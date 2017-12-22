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

#if YGL_OPENGL
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#else
#include <GL/glew.h>
#endif
#define GLFW_INCLUDE_GLCOREARB
#include <GLFW/glfw3.h>
#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw_gl3.h"
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace ygl {

namespace __impl_perlin {
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

}  // namespace __impl_perlin

// adapeted  stb_perlin.h
float perlin_noise(const vec3f& p, const vec3i& wrap) {
    return __impl_perlin::stb_perlin_noise3(
        p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_ridge_noise(const vec3f& p, float lacunarity, float gain,
    float offset, int octaves, const vec3i& wrap) {
    return __impl_perlin::stb_perlin_ridge_noise3(p.x, p.y, p.z, lacunarity,
        gain, offset, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_fbm_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
    return __impl_perlin::stb_perlin_fbm_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_turbulence_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
    return __impl_perlin::stb_perlin_turbulence_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
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
    return image4b(w, h, (vec4b*)pixels.get());
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
    return image4f(w, h, (vec4f*)pixels.get());
}

// Saves an ldr image.
bool save_image4b(const string& filename, const image4b& img) {
    if (path_extension(filename) == ".png") {
        return stbi_write_png(filename.c_str(), img.width(), img.height(), 4,
            (byte*)img.data(), img.width() * 4);
    } else if (path_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), img.width(), img.height(), 4,
            (byte*)img.data(), 75);
    } else {
        return false;
    }
}

// Saves an hdr image.
bool save_image4f(const string& filename, const image4f& img) {
    if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(
            filename.c_str(), img.width(), img.height(), 4, (float*)img.data());
    } else if (path_extension(filename) == ".exr") {
        return !SaveEXR(
            (float*)img.data(), img.width(), img.height(), 4, filename.c_str());
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

    stbir_resize_float_generic((float*)img.data(), img.width(), img.height(),
        sizeof(vec4f) * img.width(), (float*)res_img.data(), res_img.width(),
        res_img.height(), sizeof(vec4f) * res_img.width(), 4, 3,
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

    stbir_resize_uint8_generic((unsigned char*)img.data(), img.width(),
        img.height(), sizeof(vec4b) * img.width(),
        (unsigned char*)res_img.data(), res_img.width(), res_img.height(),
        sizeof(vec4b) * res_img.width(), 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

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
        if (!has_ground && j > img.height() / 2) continue;
        auto theta = pif * ((j + 0.5f) / img.height());
        theta = clamp(theta, 0.0f, pif / 2);
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
// IMPLEMENTATION FOR PATH TRACE
// -----------------------------------------------------------------------------
namespace ygl {

namespace _impl_trace {

// Phong exponent to roughness. Public API, see above.
inline float specular_exponent_to_roughness(float n) {
    return sqrtf(2 / (n + 2));
}

// Specular to fresnel eta. Public API, see above.
inline void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es = {(1 + sqrt(ks.x)) / (1 - sqrt(ks.x)),
        (1 + sqrt(ks.y)) / (1 - sqrt(ks.y)),
        (1 + sqrt(ks.z)) / (1 - sqrt(ks.z))};
    esk = {0, 0, 0};
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

// Random number smp. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
struct sampler {
    rng_pcg32& rng;        // random number state
    uint32_t pixel_hash;   // pixel hash
    int s, d;              // sample and dimension indices
    int ns, ns2;           // number of samples and its square root
    trace_rng_type rtype;  // random number type
};

// Initialize a smp ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
inline sampler make_sampler(
    rng_pcg32& rng, int i, int j, int s, int ns, trace_rng_type rtype) {
    // we use various hashes to scramble the pixel values
    return {rng, hash_uint32((uint32_t)(j + 1) << 16 | (uint32_t)(i + 1)), s, 0,
        ns, (int)round(sqrt((float)ns)), rtype};
}

// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
inline float sample_next1f(sampler& smp) {
    switch (smp.rtype) {
        case trace_rng_type::uniform: {
            return clamp(next_rand1f(smp.rng), 0.0f, 1 - flt_eps);
        } break;
        case trace_rng_type::stratified: {
            smp.d += 1;
            auto p = hash_uint64_32(
                (uint64_t)smp.pixel_hash | (uint64_t)smp.d << 32);
            auto s = hash_permute(smp.s, smp.ns, p);
            return clamp(
                (s + next_rand1f(smp.rng)) / smp.ns, 0.0f, 1 - flt_eps);
        } break;
        default: {
            assert(false);
            return 0;
        }
    }
}

// Generates a 2-dimensional sample.
//
// Implementation notes: see above. Note that using deterministic keyed
// permutaton we can use stratified sampling without preallocating samples.
inline vec2f sample_next2f(sampler& smp) {
    switch (smp.rtype) {
        case trace_rng_type::uniform: {
            return {next_rand1f(smp.rng), next_rand1f(smp.rng)};
        } break;
        case trace_rng_type::stratified: {
            smp.d += 2;
            auto p = hash_uint64_32(
                (uint64_t)smp.pixel_hash | (uint64_t)smp.d << 32);
            auto s = hash_permute(smp.s, smp.ns, p);
            return {clamp((s % smp.ns2 + next_rand1f(smp.rng)) / smp.ns2, 0.0f,
                        1 - flt_eps),
                clamp((s / smp.ns2 + next_rand1f(smp.rng)) / smp.ns2, 0.0f,
                    1 - flt_eps)};
        } break;
        default: {
            assert(false);
            return {0, 0};
        }
    }
}

// Creates a 1-dimensional sample in [0,num-1]
inline int sample_next1i(sampler& smp, int num) {
    return clamp(int(sample_next1f(smp) * num), 0, num - 1);
}

// Brdf type
enum struct brdf_type { none = 0, microfacet = 1, kajiya_kay = 2, point = 3 };

// Brdf
struct brdf {
    brdf_type type = brdf_type::none;  // type
    vec3f kd = {0, 0, 0};              // diffuse
    vec3f ks = {0, 0, 0};              // specular
    float rs = 0;                      // specular roughness
    vec3f kt = {0, 0, 0};              // transmission (thin glass)
    operator bool() const { return type != brdf_type::none; }
    vec3f rho() const { return kd + ks + kt; }
};

// Emission type
enum struct emission_type {
    none = 0,
    diffuse = 1,
    point = 2,
    line = 3,
    env = 4,
};

// Emission
struct emission {
    emission_type type = emission_type::none;
    vec3f ke = zero3f;
    operator bool() const { return type != emission_type::none; }
};

// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
struct point {
    const instance* ist = nullptr;     // instance
    const environment* env = nullptr;  // environment
    frame3f frame = identity_frame3f;  // local frame
    vec3f wo = zero3f;                 // outgoing direction
    emission em = {};                  // emission
    brdf fr = {};                      // brdf
};

// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
inline ray3f eval_camera(const camera* cam, const vec2f& uv, const vec2f& luv) {
    auto h = 2 * tan(cam->yfov / 2);
    auto w = h * cam->aspect;
    auto o = vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    auto q = vec3f{w * cam->focus * (uv.x - 0.5f),
        h * cam->focus * (uv.y - 0.5f), -cam->focus};
    return ray3f(transform_point(cam->frame, o),
        transform_direction(cam->frame, normalize(q - o)));
}

// Evaluates emission.
inline vec3f eval_emission(const point& pt) {
    auto& em = pt.em;
    auto& wo = pt.wo;
    auto& wn = pt.frame.z;

    if (!em) return zero3f;
    auto ke = zero3f;
    switch (em.type) {
        case emission_type::diffuse:
            ke += (dot(wn, wo) > 0) ? em.ke : zero3f;
            break;
        case emission_type::point: ke += em.ke; break;
        case emission_type::line: ke += em.ke; break;
        case emission_type::env: ke += em.ke; break;
        default: assert(false); break;
    }
    return ke;
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
inline vec3f eval_fresnel_dielectric(float cosw, const vec3f& eta_) {
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
inline vec3f eval_fresnel_metal(
    float cosw, const vec3f& eta, const vec3f& etak) {
    if (etak == zero3f) return eval_fresnel_dielectric(cosw, eta);

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
inline vec3f eval_fresnel_schlick(const vec3f& ks, float cosw) {
    return ks +
           (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

// Schlick approximation of Fresnel term weighted by roughness.
// This is a hack, but works better than not doing it.
inline vec3f eval_fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    auto fks = eval_fresnel_schlick(ks, cosw);
    return lerp(ks, fks, rs);
}

// Evaluates the GGX distribution and geometric term
inline float eval_ggx(float rs, float ndh, float ndi, float ndo) {
    // evaluate GGX
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pif * di * di);
#ifndef YTRACE_GGX_SMITH
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
inline float pdf_ggx(float rs, float ndh) {
    auto cos2 = ndh * ndh;
    auto tan2 = (1 - cos2) / cos2;
    auto alpha2 = rs * rs;
    auto d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
    return d;
}

// Sample the GGX distribution
inline vec3f sample_ggx(float rs, const vec2f& rn) {
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz = sqrt(1 / (tan2 + 1)), rr = sqrt(1 - rz * rz),
         rphi = 2 * pif * rn.x;
    // set to wh
    auto wh_local = vec3f{rr * cos(rphi), rr * sin(rphi), rz};
    return wh_local;
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
inline vec3f eval_brdfcos(const point& pt, const vec3f& wi) {
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
        case brdf_type::microfacet: {
            // compute wh
            auto wh = normalize(wo + wi);

            // compute dot products
            auto ndo = dot(wn, wo), ndi = dot(wn, wi),
                 ndh = clamp(dot(wh, wn), (float)-1, (float)1);

            // diffuse term
            if (fr.kd != zero3f && ndi > 0 && ndo > 0) {
                brdfcos += fr.kd * ndi / pif;
            }

            // specular term
            if (fr.ks != zero3f && ndi > 0 && ndo > 0 && ndh > 0) {
                // microfacet term
                auto dg = eval_ggx(fr.rs, ndh, ndi, ndo);

                // handle fresnel
                auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
                auto ks = eval_fresnel_schlick(fr.ks, odh, fr.rs);

                // sum up
                brdfcos += ks * ndi * dg / (4 * ndi * ndo);
            }

            // transmission hack
            if (fr.kt != zero3f && wo == -wi) brdfcos += fr.kt;
        } break;
        // hair (Kajiya-Kay)
        case brdf_type::kajiya_kay: {
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
        case brdf_type::point: {
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
inline float weight_brdfcos(const point& pt, const vec3f& wi) {
    // grab variables
    auto& fr = pt.fr;
    auto& wn = pt.frame.z;
    auto& wo = pt.wo;

    // skip if no component
    if (!fr) return 0;

    // probability of each lobe
    auto kdw = max_element(fr.kd).second, ksw = max_element(fr.ks).second,
         ktw = max_element(fr.kt).second;
    auto kaw = kdw + ksw + ktw;
    kdw /= kaw;
    ksw /= kaw;
    ktw /= kaw;

    // accumulate the probability over all lobes
    auto pdf = 0.0f;
    // sample the lobe
    switch (fr.type) {
        // reflection term
        case brdf_type::microfacet: {
            // compute wh
            auto wh = normalize(wi + wo);

            // compute dot products
            auto ndo = dot(wn, wo), ndi = dot(wn, wi), ndh = dot(wn, wh);

            // diffuse term (hemipherical cosine probability)
            if (kdw && ndo > 0 && ndi > 0) { pdf += kdw * ndi / pif; }

            // specular term (GGX)
            if (ksw && ndo > 0 && ndi > 0 && ndh > 0) {
                // probability proportional to d adjusted by wh projection
                auto d = pdf_ggx(fr.rs, ndh);
                auto hdo = dot(wo, wh);
                pdf += ksw * d / (4 * hdo);
            }

            // transmission hack
            if (ktw && wi == -wo) pdf += ktw;

            // check
            assert(isfinite(pdf));
        } break;
        // hair (Kajiya-Kay)
        case brdf_type::kajiya_kay: {
            // diffuse and specular
            pdf += (kdw + ksw) * 4 * pif;
            // transmission hack
            if (wi == -wo) pdf += ktw;
        } break;
        // point
        case brdf_type::point: {
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
inline vec3f sample_brdfcos(const point& pt, float rnl, const vec2f& rn) {
    // grab variables
    auto& fr = pt.fr;
    auto& wn = pt.frame.z;
    auto& fp = pt.frame;
    auto& wo = pt.wo;

    // skip if no component
    if (!fr) return zero3f;

    // probability of each lobe
    auto kdw = max_element(fr.kd).second, ksw = max_element(fr.ks).second,
         ktw = max_element(fr.kt).second;
    auto kaw = kdw + ksw + ktw;
    kdw /= kaw;
    ksw /= kaw;
    ktw /= kaw;

    // sample selected lobe
    switch (fr.type) {
        // reflection term
        case brdf_type::microfacet: {
            // compute cosine
            auto ndo = dot(wn, wo);

            // check to make sure we are above the surface
            if (ndo <= 0) return zero3f;

            // sample according to diffuse
            if (rnl < kdw) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                // set to wi
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(fp, wi_local);
            }
            // sample according to specular GGX
            else if (rnl < kdw + ksw) {
                // sample wh with ggx distribution
                auto wh_local = sample_ggx(fr.rs, rn);
                auto wh = transform_direction(fp, wh_local);
                // compute wi
                return normalize(wh * 2.0f * dot(wo, wh) - wo);
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return -wo;
            } else
                assert(false);
        } break;
        // hair (Kajiya-Kay)
        case brdf_type::kajiya_kay: {
            // diffuse and specular
            if (rnl < kdw + ksw) {
                // sample wi with uniform spherical distribution
                auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(fp, wi_local);
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return -wo;
            } else
                assert(false);
        } break;
        // diffuse term point
        case brdf_type::point: {
            // diffuse and specular
            if (rnl < kdw + ksw) {
                // sample wi with uniform spherical distribution
                auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(fp, wi_local);
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return -wo;
            } else
                assert(false);
        } break;
        default: assert(false); break;
    }

    // done
    return zero3f;
}

// Create a point for an environment map. Resolves material with textures.
inline point eval_envpoint(const environment* env, const vec3f& wo) {
    // set shape data
    auto pt = point();

    // env
    pt.env = env;

    // direction
    pt.wo = wo;

    // maerial
    auto ke = env->ke;
    if (env->ke_txt) {
        auto w = transform_direction_inverse(env->frame, -wo);
        auto theta = clamp(w.y, (float)-1, (float)1);
        auto phi = atan2(w.z, w.x);
        auto texcoord = vec2f{0.5f + phi / (2 * pif), 1 + 0.5f * theta / pif};
        ke *= eval_texture(env->ke_txt, texcoord).xyz();
    }

    // create emission lobe
    if (ke != zero3f) { pt.em = {emission_type::env, ke}; }

    // done
    return pt;
}

// Create a point for a shape. Resolves geometry and material with textures.
inline point eval_shapepoint(
    const instance* ist, int eid, const vec4f& euv, const vec3f& wo) {
    // set shape data
    auto pt = point();

    // instance
    pt.ist = ist;

    // direction
    pt.wo = wo;

    // shortcuts
    auto shp = ist->shp;
    auto mat = ist->shp->mat;

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
            make_frame3_fromzx({0, 0, 0}, norm, {tangsp.x, tangsp.y, tangsp.z});
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
    switch (mat->mtype) {
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
            pt.fr.type = brdf_type::point;
        if (ke != zero3f) pt.em.type = emission_type::point;
    } else if (!shp->lines.empty()) {
        if (kd.xyz() != zero3f || ks.xyz() != zero3f || kt.xyz() != zero3f)
            pt.fr.type = brdf_type::kajiya_kay;
        if (ke != zero3f) pt.em.type = emission_type::line;
    } else if (!shp->triangles.empty()) {
        if (kd.xyz() != zero3f || ks.xyz() != zero3f || kt.xyz() != zero3f)
            pt.fr.type = brdf_type::microfacet;
        if (ke != zero3f) pt.em.type = emission_type::diffuse;
    }

    // done
    return pt;
}

// Sample weight for a light point.
inline float weight_light(const point& lpt, const point& pt) {
    if (!lpt.em) return 0;
    // support only one lobe for now
    switch (lpt.em.type) {
        case emission_type::env: {
            return 4 * pif;
        } break;
        case emission_type::point: {
            auto d = length(lpt.frame.o - pt.frame.o);
            return lpt.ist->shp->elem_cdf.back() / (d * d);
        } break;
        case emission_type::line: {
            assert(false);
            return 0;
        } break;
        case emission_type::diffuse: {
            auto d = length(lpt.frame.o - pt.frame.o);
            return lpt.ist->shp->elem_cdf.back() *
                   abs(dot(lpt.frame.z, lpt.wo)) / (d * d);
        } break;
        default: {
            assert(false);
            return 0;
        } break;
    }
}

// Picks a point on a light.
inline point sample_light(
    const light* lgt, const point& pt, float rne, const vec2f& rn) {
    if (lgt->ist) {
        auto shp = lgt->ist->shp;
        auto eid = 0;
        auto euv = zero4f;
        if (!shp->triangles.empty()) {
            std::tie(eid, (vec3f&)euv) =
                sample_triangles(shp->elem_cdf, rne, rn);
        } else if (!shp->quads.empty()) {
            std::tie(eid, (vec4f&)euv) = sample_quads(shp->elem_cdf, rne, rn);
        } else if (!shp->lines.empty()) {
            std::tie(eid, (vec2f&)euv) = sample_lines(shp->elem_cdf, rne, rn.x);
        } else if (!shp->points.empty()) {
            eid = sample_points(shp->elem_cdf, rne);
            euv = {1, 0, 0, 0};
        } else {
            assert(false);
        }
        auto lpt = eval_shapepoint(lgt->ist, eid, euv, zero3f);
        lpt.wo = normalize(pt.frame.o - lpt.frame.o);
        return lpt;
    } else if (lgt->env) {
        auto z = -1 + 2 * rn.y;
        auto rr = sqrt(clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * pif * rn.x;
        auto wo = vec3f{cos(phi) * rr, z, sin(phi) * rr};
        auto lpt = eval_envpoint(lgt->env, wo);
        return lpt;
    } else {
        assert(false);
        return {};
    }
}

// Offsets a ray origin to avoid self-intersection.
inline ray3f offset_ray(
    const point& pt, const vec3f& w, const trace_params& params) {
    if (dot(w, pt.frame.z) > 0) {
        return ray3f(
            pt.frame.o + pt.frame.z * params.ray_eps, w, params.ray_eps);
    } else {
        return ray3f(
            pt.frame.o - pt.frame.z * params.ray_eps, w, params.ray_eps);
    }
}

// Offsets a ray origin to avoid self-intersection.
inline ray3f offset_ray(
    const point& pt, const point& pt2, const trace_params& params) {
    auto ray_dist = (!pt2.env) ? length(pt.frame.o - pt2.frame.o) : flt_max;
    if (dot(-pt2.wo, pt.frame.z) > 0) {
        return ray3f(pt.frame.o + pt.frame.z * params.ray_eps, -pt2.wo,
            params.ray_eps, ray_dist - 2 * params.ray_eps);
    } else {
        return ray3f(pt.frame.o - pt.frame.z * params.ray_eps, -pt2.wo,
            params.ray_eps, ray_dist - 2 * params.ray_eps);
    }
}

// Intersects a ray with the scn and return the point (or env
// point).
inline point intersect_scene(const scene* scn, const ray3f& ray) {
    auto iid = 0, eid = 0;
    auto euv = zero4f;
    auto ray_t = 0.0f;
    if (intersect_ray(scn, ray, false, ray_t, iid, eid, euv)) {
        return eval_shapepoint(scn->instances[iid], eid, euv, -ray.d);
    } else if (!scn->environments.empty()) {
        return eval_envpoint(scn->environments[0], -ray.d);
    } else {
        return {};
    }
}

// Test occlusion
inline vec3f eval_transmission(const scene* scn, const point& pt,
    const point& lpt, const trace_params& params) {
    if (params.shadow_notransmission) {
        auto shadow_ray = offset_ray(pt, lpt, params);
        // auto shadow_ray = ray3f{pt.frame.o, -lpt.wo, 0.01f, flt_max};
        return (intersect_ray(scn, shadow_ray, true)) ? zero3f : vec3f{1, 1, 1};
    } else {
        auto cpt = pt;
        auto weight = vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            cpt = intersect_scene(scn, offset_ray(cpt, lpt, params));
            if (!cpt.ist) break;
            weight *= cpt.fr.kt;
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight
inline float weight_mis(float w0, float w1) {
    if (!w0 || !w1) return 1;
    return (1 / w0) / (1 / w0 + 1 / w1);
}

// Recursive path tracing.
inline vec3f eval_li_pathtrace(const scene* scn, const ray3f& ray, sampler& smp,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    hit = pt.ist;

    // emission
    auto l = eval_emission(pt);
    if (!pt.fr || scn->lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt);

        // direct – light
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto lw = weight_light(lpt, pt) * (float)scn->lights.size();
        auto lke = eval_emission(lpt);
        auto lbc = eval_brdfcos(pt, -lpt.wo);
        auto lld = lke * lbc * lw;
        if (lld != zero3f) {
            l += weight * lld * eval_transmission(scn, pt, lpt, params) *
                 weight_mis(lw, weight_brdfcos(pt, -lpt.wo));
        }

        // direct – brdf
        auto bpt = intersect_scene(
            scn, offset_ray(pt,
                     sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp)),
                     params));
        auto bw = weight_brdfcos(pt, -bpt.wo);
        auto bke = eval_emission(bpt);
        auto bbc = eval_brdfcos(pt, -bpt.wo);
        auto bld = bke * bbc * bw;
        if (bld != zero3f) {
            l += weight * bld * weight_mis(bw, weight_light(bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;
        if (!bpt.fr) break;

        // continue path
        weight *= eval_brdfcos(pt, -bpt.wo) * weight_brdfcos(pt, -bpt.wo);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element(pt.fr.rho()).second, 0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        emission = false;
    }

    return l;
}

// Recursive path tracing.
inline vec3f eval_li_pathtrace_nomis(const scene* scn, const ray3f& ray,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    hit = pt.ist;

    // emission
    auto l = eval_emission(pt);
    if (!pt.fr) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt);

        // direct
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element(pt.fr.rho()).second, 0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto bwi = sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp));
        weight *= eval_brdfcos(pt, bwi) * weight_brdfcos(pt, bwi);
        if (weight == zero3f) break;

        auto bpt = intersect_scene(scn, offset_ray(pt, bwi, params));
        emission = false;
        if (!bpt.fr) break;

        // continue path
        pt = bpt;
    }

    return l;
}

// Recursive path tracing.
inline vec3f eval_li_pathtrace_hack(const scene* scn, const ray3f& ray,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    hit = pt.ist;

    // emission
    auto l = eval_emission(pt);
    if (!pt.fr || scn->lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // direct
        auto lgt = scn->lights[sample_next1i(smp, (int)scn->lights.size())];
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt) * (float)scn->lights.size();
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element(pt.fr.rho()).second, 0.95f);
            if (sample_next1f(smp) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto bwi = sample_brdfcos(pt, sample_next1f(smp), sample_next2f(smp));
        weight *= eval_brdfcos(pt, bwi) * weight_brdfcos(pt, bwi);
        if (weight == zero3f) break;

        auto bpt = intersect_scene(scn, offset_ray(pt, bwi, params));
        if (!bpt.fr) break;

        // continue path
        pt = bpt;
    }

    return l;
}

// Direct illumination.
inline vec3f eval_li_direct(const scene* scn, const ray3f& ray, int bounce,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    if (!bounce) hit = pt.ist;

    // emission
    auto l = eval_emission(pt);
    if (!pt.fr || scn->lights.empty()) return l;

    // ambient
    l += params.amb * pt.fr.rho();

    // direct
    for (auto& lgt : scn->lights) {
        auto lpt =
            sample_light(lgt, pt, sample_next1f(smp), sample_next2f(smp));
        auto ld = eval_emission(lpt) * eval_brdfcos(pt, -lpt.wo) *
                  weight_light(lpt, pt);
        if (ld == zero3f) continue;
        l += ld * eval_transmission(scn, pt, lpt, params);
    }

    // exit if needed
    if (bounce >= params.max_depth) return l;

    // opacity
    if (pt.fr.kt != zero3f) {
        auto ray = offset_ray(pt, -pt.wo, params);
        l += pt.fr.kt * eval_li_direct(scn, ray, bounce + 1, smp, params, hit);
    }

    // done
    return l;
}

// Direct illumination.
inline vec3f eval_li_direct(const scene* scn, const ray3f& ray, sampler& smp,
    const trace_params& params, bool& hit) {
    return eval_li_direct(scn, ray, 0, smp, params, hit);
}

// Eyelight for quick previewing.
inline vec3f eval_li_eyelight(const scene* scn, const ray3f& ray, int bounce,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    if (!bounce) hit = pt.ist;

    // emission
    auto l = eval_emission(pt);
    if (!pt.fr) return l;

    // brdf*light
    l += eval_brdfcos(pt, pt.wo) * pif;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.fr.kt != zero3f) {
        auto ray = offset_ray(pt, -pt.wo, params);
        l +=
            pt.fr.kt * eval_li_eyelight(scn, ray, bounce + 1, smp, params, hit);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
inline vec3f eval_li_eyelight(const scene* scn, const ray3f& ray, sampler& smp,
    const trace_params& params, bool& hit) {
    return eval_li_eyelight(scn, ray, 0, smp, params, hit);
}

// Debug previewing.
inline vec3f eval_li_debug_normal(const scene* scn, const ray3f& ray,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto isec = intersect_ray(scn, ray, false);
    hit = (bool)isec;
    if (!hit) return {0, 0, 0};

    // texcoord
    auto norm = eval_norm(scn->instances[isec.iid], isec.eid, isec.euv);
    return norm * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug previewing.
inline vec3f eval_li_debug_albedo(const scene* scn, const ray3f& ray,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto pt = intersect_scene(scn, ray);
    hit = pt.ist;

    return pt.fr.rho();
}

// Debug previewing.
inline vec3f eval_li_debug_texcoord(const scene* scn, const ray3f& ray,
    sampler& smp, const trace_params& params, bool& hit) {
    // intersection
    auto isec = intersect_ray(scn, ray, false);
    hit = (bool)isec;
    if (!hit) return {0, 0, 0};

    // texcoord
    auto texcoord =
        eval_texcoord(scn->instances[isec.iid]->shp, isec.eid, isec.euv);
    return {texcoord.x, texcoord.y, 0};
}

// Shader function callback.
using eval_li_fn = vec3f (*)(const scene* scn, const ray3f& ray, sampler& smp,
    const trace_params& params, bool& hit);

// Get a shader function
inline eval_li_fn get_shader(const trace_params& params) {
    switch (params.stype) {
        case trace_shader_type::eyelight: return eval_li_eyelight;
        case trace_shader_type::direct: return eval_li_direct;
        case trace_shader_type::pathtrace: return eval_li_pathtrace;
        case trace_shader_type::pathtrace_nomis: return eval_li_pathtrace_nomis;
        case trace_shader_type::debug_albedo: return eval_li_debug_albedo;
        case trace_shader_type::debug_normal: return eval_li_debug_normal;
        case trace_shader_type::debug_texcoord: return eval_li_debug_texcoord;
        default: return nullptr;
    }
}

// triangle filter (public domain from stb_image_resize)
inline float filter_triangle(float x) {
    x = (float)fabs(x);

    if (x <= 1.0f)
        return 1 - x;
    else
        return 0;
}

// cubic filter (public domain from stb_image_resize)
inline float filter_cubic(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return (4 + x * x * (3 * x - 6)) / 6;
    else if (x < 2.0f)
        return (8 + x * (-12 + x * (6 - x))) / 6;
    else
        return 0.0f;
}

// catmull-rom filter (public domain from stb_image_resize)
inline float filter_catmullrom(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return 1 - x * x * (2.5f - 1.5f * x);
    else if (x < 2.0f)
        return 2 - x * (4 + x * (0.5f * x - 2.5f));
    else
        return 0.0f;
}

// mitchell filter (public domain from stb_image_resize)
inline float filter_mitchell(float x) {
    x = (float)fabs(x);
    if (x < 1.0f)
        return (16 + x * x * (21 * x - 36)) / 18;
    else if (x < 2.0f)
        return (32 + x * (-60 + x * (36 - 7 * x))) / 18;
    else
        return 0.0f;
}

// filter function
using filter_fn = float (*)(float);

// Get a filter function
inline filter_fn get_filter(const trace_params& params) {
    switch (params.ftype) {
        case trace_filter_type::box: return nullptr;
        case trace_filter_type::triangle: return filter_triangle;
        case trace_filter_type::cubic: return filter_cubic;
        case trace_filter_type::catmull_rom: return filter_catmullrom;
        case trace_filter_type::mitchell: return filter_mitchell;
        default: return nullptr;
    }
}

// Get a filter size
inline int get_filter_size(const trace_params& params) {
    switch (params.ftype) {
        case trace_filter_type::box: return 0;
        case trace_filter_type::triangle: return 1;
        case trace_filter_type::cubic:
        case trace_filter_type::catmull_rom:
        case trace_filter_type::mitchell: return 2;
        default: return 0;
    }
}

// Renders a block of pixels. Public API, see above.
inline void trace_block(const scene* scn, image4f& img, int block_x,
    int block_y, int block_width, int block_height, int samples_min,
    int samples_max, vector<rng_pcg32>& rngs, const trace_params& params) {
    auto cam = scn->cameras[params.camera_id];
    auto shade = get_shader(params);
    for (auto j = block_y; j < block_y + block_height; j++) {
        for (auto i = block_x; i < block_x + block_width; i++) {
            auto lp = zero4f;
            for (auto s = samples_min; s < samples_max; s++) {
                auto smp = make_sampler(rngs[j * params.width + i], i, j, s,
                    params.nsamples, params.rtype);
                auto rn = sample_next2f(smp);
                auto uv = vec2f{
                    (i + rn.x) / params.width, 1 - (j + rn.y) / params.height};
                auto ray = eval_camera(cam, uv, sample_next2f(smp));
                auto hit = false;
                auto l = shade(scn, ray, smp, params, hit);
                if (!hit && params.envmap_invisible) continue;
                if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
                    log_error("NaN detected");
                    continue;
                }
                if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
                lp += {l, 1};
            }
            if (samples_min) {
                img[{i, j}] = (img[{i, j}] * (float)samples_min + lp) /
                              (float)samples_max;
            } else {
                img[{i, j}] = lp / (float)samples_max;
            }
        }
    }
}

// Trace a block of samples
inline void trace_block(const scene* scn, image4f& img, const vec2i& block_min,
    const vec2i& block_max, int samples_min, int samples_max,
    vector<rng_pcg32>& rngs, const trace_params& params) {
    auto shade = get_shader(params);
    auto cam = scn->cameras[params.camera_id];
    for (auto j = block_min.y; j < block_max.y; j++) {
        for (auto i = block_min.x; i < block_max.x; i++) {
            auto lp = zero4f;
            for (auto s = samples_min; s < samples_max; s++) {
                auto smp = make_sampler(rngs[j * params.width + i], i, j, s,
                    params.nsamples, params.rtype);
                auto rn = sample_next2f(smp);
                auto uv = vec2f{
                    (i + rn.x) / params.width, 1 - (j + rn.y) / params.height};
                auto ray = eval_camera(cam, uv, sample_next2f(smp));
                bool hit = false;
                auto l = shade(scn, ray, smp, params, hit);
                if (!hit && params.envmap_invisible) continue;
                if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
                    log_error("NaN detected");
                    continue;
                }
                if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
                lp += {l, 1};
            }
            if (samples_min) {
                img[{i, j}] = (img[{i, j}] * (float)samples_min + lp) /
                              (float)samples_max;
            } else {
                img[{i, j}] = lp / (float)samples_max;
            }
        }
    }
}

// Trace a block of samples
inline void trace_block_filtered(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, const vec2i& block_min, const vec2i& block_max,
    int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    std::mutex& image_mutex, const trace_params& params) {
    auto shade = get_shader(params);
    auto cam = scn->cameras[params.camera_id];
    auto filter = get_filter(params);
    auto filter_size = get_filter_size(params);
    static constexpr const int pad = 2;
    auto block_size =
        vec2i{block_max.x - block_min.x, block_max.y - block_min.y};
    auto acc_buffer = image4f(block_size.x + pad * 2, block_size.y + pad * 2);
    auto weight_buffer =
        image4f(block_size.x + pad * 2, block_size.y + pad * 2);
    for (auto j = block_min.y; j < block_max.y; j++) {
        for (auto i = block_min.x; i < block_max.x; i++) {
            for (auto s = samples_min; s < samples_max; s++) {
                auto smp = make_sampler(rngs[j * params.width + i], i, j, s,
                    params.nsamples, params.rtype);
                auto rn = sample_next2f(smp);
                auto uv = vec2f{
                    (i + rn.x) / params.width, 1 - (j + rn.y) / params.height};
                auto ray = eval_camera(cam, uv, sample_next2f(smp));
                auto hit = false;
                auto l = shade(scn, ray, smp, params, hit);
                if (!hit && params.envmap_invisible) continue;
                if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
                    log_error("NaN detected");
                    continue;
                }
                if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
                if (params.ftype == trace_filter_type::box) {
                    auto bi = i - block_min.x, bj = j - block_min.y;
                    acc_buffer[{bi + pad, bj + pad}] += {l, 1};
                    weight_buffer[{bi + pad, bj + pad}] += {1, 1, 1, 1};
                } else {
                    auto bi = i - block_min.x, bj = j - block_min.y;
                    for (auto fj = -filter_size; fj <= filter_size; fj++) {
                        for (auto fi = -filter_size; fi <= filter_size; fi++) {
                            auto w = filter(fi - uv.x + 0.5f) *
                                     filter(fj - uv.y + 0.5f);
                            acc_buffer[{bi + fi + pad, bj + fj + pad}] +=
                                {l * w, w};
                            weight_buffer[{bi + fi + pad, bj + fj + pad}] +=
                                {w, w, w, w};
                        }
                    }
                }
            }
        }
    }
    if (params.ftype == trace_filter_type::box) {
        for (auto j = block_min.y; j < block_max.y; j++) {
            for (auto i = block_min.x; i < block_max.x; i++) {
                auto bi = i - block_min.x, bj = j - block_min.y;
                acc[{i, j}] += acc_buffer[{bi + pad, bj + pad}];
                weight[{i, j}] += weight_buffer[{bi + pad, bj + pad}];
                img[{i, j}] = acc[{i, j}] / weight[{i, j}];
            }
        }
    } else {
        std::unique_lock<std::mutex> lock_guard(image_mutex);
        auto width = acc.width(), height = acc.height();
        for (auto j = max(block_min.y - filter_size, 0);
             j < min(block_max.y + filter_size, height); j++) {
            for (auto i = max(block_min.x - filter_size, 0);
                 i < min(block_max.x + filter_size, width); i++) {
                auto bi = i - block_min.x, bj = j - block_min.y;
                acc[{i, j}] += acc_buffer[{bi + pad, bj + pad}];
                weight[{i, j}] += weight_buffer[{bi + pad, bj + pad}];
                img[{i, j}] = acc[{i, j}] / weight[{i, j}];
            }
        }
    }
}

}  // namespace _impl_trace

// Renders a block of samples
void trace_block(const scene* scn, image4f& img, const vec2i& block_min,
    const vec2i& block_max, int samples_min, int samples_max,
    vector<rng_pcg32>& rngs, const trace_params& params) {
    _impl_trace::trace_block(
        scn, img, block_min, block_max, samples_min, samples_max, rngs, params);
}

// Renders a filtered block of samples
void trace_block_filtered(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, const vec2i& block_min, const vec2i& block_max,
    int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    std::mutex& image_mutex, const trace_params& params) {
    _impl_trace::trace_block_filtered(scn, img, acc, weight, block_min,
        block_max, samples_min, samples_max, rngs, image_mutex, params);
}

// Trace the next samples in [samples_min, samples_max) range.
// Samples have to be traced consecutively.
void trace_samples(const scene* scn, image4f& img, int samples_min,
    int samples_max, vector<rng_pcg32>& rngs, const trace_params& params) {
    auto blocks = trace_blocks(params);
    if (params.parallel) {
        parallel_for((int)blocks.size(), [&img, scn, samples_min, samples_max,
                                             &blocks, &params, &rngs](int idx) {
            trace_block(scn, img, blocks[idx].first, blocks[idx].second,
                samples_min, samples_max, rngs, params);
        });
    } else {
        for (auto idx = 0; idx < (int)blocks.size(); idx++) {
            trace_block(scn, img, blocks[idx].first, blocks[idx].second,
                samples_min, samples_max, rngs, params);
        }
    }
}

// Trace the next samples in [samples_min, samples_max) range.
// Samples have to be traced consecutively.
void trace_filtered_samples(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    const trace_params& params) {
    auto blocks = trace_blocks(params);
    std::mutex image_mutex;
    if (params.parallel) {
        parallel_for((int)blocks.size(),
            [&img, &acc, &weight, scn, samples_min, samples_max, &blocks,
                &params, &image_mutex, &rngs](int idx) {
                trace_block_filtered(scn, img, acc, weight, blocks[idx].first,
                    blocks[idx].second, samples_min, samples_max, rngs,
                    image_mutex, params);
            });
    } else {
        for (auto idx = 0; idx < (int)blocks.size(); idx++) {
            trace_block_filtered(scn, img, acc, weight, blocks[idx].first,
                blocks[idx].second, samples_min, samples_max, rngs, image_mutex,
                params);
        }
    }
}

// Starts an anyncrhounous renderer with a maximum of 256 samples.
void trace_async_start(const scene* scn, image4f& img, vector<rng_pcg32>& rngs,
    const trace_params& params, thread_pool* pool,
    const function<void(int)>& callback) {
    auto blocks = trace_blocks(params);
    for (auto sample = 0; sample < params.nsamples; sample++) {
        for (auto& block : blocks) {
            auto is_last = (block == blocks.back());
            run_async(pool, [&img, scn, sample, block, &params, callback, &rngs,
                                is_last]() {
                trace_block(scn, img, block.first, block.second, sample,
                    sample + 1, rngs, params);
                if (is_last) callback(sample);
            });
        }
    }
}

// Stop the asynchronous renderer.
void trace_async_stop(thread_pool* pool) {
    if (!pool) return;
    clear_pool(pool);
    wait_pool(pool);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

namespace _impl_obj {

// Parse a value
template <typename T>
inline void parse_val(stringstream& ss, T& v) {
    ss >> v;
}

// Parse a value
inline void parse_val(stringstream& ss, vec2f& v) {
    for (auto i = 0; i < 2; i++) parse_val(ss, v[i]);
}
// Parse a value
inline void parse_val(stringstream& ss, vec3f& v) {
    for (auto i = 0; i < 3; i++) parse_val(ss, v[i]);
}
// Parse a value
inline void parse_val(stringstream& ss, vec4f& v) {
    for (auto i = 0; i < 4; i++) parse_val(ss, v[i]);
}

// Parse a value
inline void parse_val(stringstream& ss, vec2i& v) {
    for (auto i = 0; i < 2; i++) parse_val(ss, v[i]);
}
// Parse a value
inline void parse_val(stringstream& ss, vec3i& v) {
    for (auto i = 0; i < 3; i++) parse_val(ss, v[i]);
}
// Parse a value
inline void parse_val(stringstream& ss, vec4i& v) {
    for (auto i = 0; i < 4; i++) parse_val(ss, v[i]);
}

// Parse a value
inline void parse_val(stringstream& ss, frame3f& v) {
    parse_val(ss, v.x);
    parse_val(ss, v.y);
    parse_val(ss, v.z);
    parse_val(ss, v.o);
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
            parse_val(ss, materials.back()->name);
        } else if (cmd == "illum") {
            parse_val(ss, materials.back()->illum);
        } else if (cmd == "Ke") {
            parse_val(ss, materials.back()->ke);
        } else if (cmd == "Ka") {
            parse_val(ss, materials.back()->ka);
        } else if (cmd == "Kd") {
            parse_val(ss, materials.back()->kd);
        } else if (cmd == "Ks") {
            parse_val(ss, materials.back()->ks);
        } else if (cmd == "Kr") {
            parse_val(ss, materials.back()->kr);
        } else if (cmd == "Kt" || cmd == "Tf") {
            auto vals = zero3f;
            auto ntok = 0;
            while (ss) parse_val(ss, vals[ntok++]);
            if (ntok >= 3)
                materials.back()->kt = vals;
            else
                materials.back()->kt = {vals.x, vals.x, vals.x};
        } else if (cmd == "Tr") {
            auto vals = zero3f;
            auto ntok = 0;
            while (ss) parse_val(ss, vals[ntok++]);
            if (ntok >= 3)
                materials.back()->kt = vals;
            else
                materials.back()->op = (flip_tr) ? 1 - vals.x : vals.x;
        } else if (cmd == "Ns") {
            parse_val(ss, materials.back()->ns);
        } else if (cmd == "d") {
            parse_val(ss, materials.back()->op);
        } else if (cmd == "Ni") {
            parse_val(ss, materials.back()->ior);
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
                parse_val(ss, materials.back()->unknown_props[cmd].back());
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
inline void parse_vertlist(
    stringstream& ss, vector<obj_vertex>& elems, const obj_vertex& vert_size) {
    elems.clear();
    while (true) {
        auto tok = string();
        parse_val(ss, tok);
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
            parse_val(ss, asset->pos.back());
        } else if (cmd == "vn") {
            vert_size.norm += 1;
            asset->norm.push_back({});
            parse_val(ss, asset->norm.back());
        } else if (cmd == "vt") {
            vert_size.texcoord += 1;
            asset->texcoord.push_back({});
            parse_val(ss, asset->texcoord.back());
            if (flip_texcoord)
                asset->texcoord.back().y = 1 - asset->texcoord.back().y;
        } else if (cmd == "vc") {
            vert_size.color += 1;
            asset->color.push_back({});
            parse_val(ss, asset->color.back());
        } else if (cmd == "vr") {
            vert_size.radius += 1;
            asset->radius.push_back({});
            parse_val(ss, asset->radius.back());
        } else if (cmd == "f") {
            parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::face,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "l") {
            parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::line,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "p") {
            parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::point, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "t") {
            parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::tetra, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "o") {
            auto name = string();
            parse_val(ss, name);
            asset->objects.push_back(new obj_object{name, {}});
            asset->objects.back()->groups.push_back({cur_matname, ""});
        } else if (cmd == "usemtl") {
            auto name = string();
            parse_val(ss, name);
            cur_matname = name;
            asset->objects.back()->groups.push_back({cur_matname, ""});
        } else if (cmd == "g") {
            auto name = string();
            parse_val(ss, name);
            asset->objects.back()->groups.push_back({cur_matname, name});
        } else if (cmd == "s") {
            auto name = string();
            parse_val(ss, name);
            auto smoothing = name == string("on");
            if (asset->objects.back()->groups.back().smoothing != smoothing) {
                asset->objects.back()->groups.push_back(
                    {cur_matname, name, smoothing});
            }
        } else if (cmd == "mtllib") {
            auto name = string();
            parse_val(ss, name);
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
            parse_val(ss, cam->name);
            parse_val(ss, cam->ortho);
            parse_val(ss, cam->yfov);
            parse_val(ss, cam->aspect);
            parse_val(ss, cam->aperture);
            parse_val(ss, cam->focus);
            parse_val(ss, cam->frame);
            asset->cameras.push_back(cam);
        } else if (cmd == "e") {
            auto env = new obj_environment();
            parse_val(ss, env->name);
            parse_val(ss, env->matname);
            parse_val(ss, env->frame);
            asset->environments.push_back(env);
        } else if (cmd == "i") {
            auto ist = new obj_instance();
            parse_val(ss, ist->name);
            parse_val(ss, ist->objname);
            parse_val(ss, ist->frame);
            asset->instances.push_back(ist);
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
            asset->textures.push_back(new obj_texture{txt});
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
inline void dump_val(fstream& fs, const T& v) {
    fs << v;
}

// write to stream
inline void dump_val(fstream& fs, const vec2f& v) {
    dump_val(fs, v.x);
    fs << ' ';
    dump_val(fs, v.y);
}

// write to stream
inline void dump_val(fstream& fs, const vec3f& v) {
    dump_val(fs, v.x);
    fs << ' ';
    dump_val(fs, v.y);
    fs << ' ';
    dump_val(fs, v.z);
}

// write to stream
inline void dump_val(fstream& fs, const vec4f& v) {
    dump_val(fs, v.x);
    fs << ' ';
    dump_val(fs, v.y);
    fs << ' ';
    dump_val(fs, v.z);
    fs << ' ';
    dump_val(fs, v.w);
}

// write to stream
inline void dump_val(fstream& fs, const frame3f& v) {
    dump_val(fs, v.x);
    fs << ' ';
    dump_val(fs, v.y);
    fs << ' ';
    dump_val(fs, v.z);
    fs << ' ';
    dump_val(fs, v.o);
}

// write to stream
inline void dump_val(fstream& fs, const obj_texture_info& v) {
    for (auto&& kv : v.unknown_props) {
        dump_val(fs, kv.first + " ");
        for (auto&& vv : kv.second) dump_val(fs, vv + " ");
    }
    if (v.clamp) dump_val(fs, "-clamp on ");
    dump_val(fs, v.path);
}

// write to stream
template <typename T>
inline void dump_named_val(fstream& fs, const string& name, const T& v) {
    dump_val(fs, name);
    fs << ' ';
    dump_val(fs, v);
    fs << '\n';
}

// write to stream
template <typename T>
inline void dump_opt_val(
    fstream& fs, const string& name, const T& v, const T& def = {}) {
    if (v == def) return;
    dump_named_val(fs, name, v);
}

// write an OBJ vertex triplet using only the indices that are active
inline void dump_objverts(
    fstream& fs, const char* str, int nv, const obj_vertex* verts) {
    dump_val(fs, str);
    for (auto v = 0; v < nv; v++) {
        auto& vert = verts[v];
        auto vert_ptr = &vert.pos;
        auto nto_write = 0;
        for (auto i = 0; i < 5; i++) {
            if (vert_ptr[i] >= 0) nto_write = i + 1;
        }
        for (auto i = 0; i < nto_write; i++) {
            if (vert_ptr[i] >= 0) {
                dump_val(fs, ((i == 0) ? ' ' : '/'));
                dump_val(fs, vert_ptr[i] + 1);
            } else {
                dump_val(fs, '/');
            }
        }
    }
    dump_val(fs, '\n');
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
        dump_named_val(fs, "newmtl", mat->name);
        dump_named_val(fs, "  illum", mat->illum);
        dump_opt_val(fs, "  Ke", mat->ke);
        dump_opt_val(fs, "  Ka", mat->ka);
        dump_opt_val(fs, "  Kd", mat->kd);
        dump_opt_val(fs, "  Ks", mat->ks);
        dump_opt_val(fs, "  Kr", mat->kr);
        dump_opt_val(fs, "  Tf", mat->kt);
        dump_opt_val(fs, "  Ns", mat->ns, 0.0f);
        dump_opt_val(fs, "  d", mat->op, 1.0f);
        dump_opt_val(fs, "  Ni", mat->ior, 1.0f);
        dump_opt_val(fs, "  map_Ke", mat->ke_txt);
        dump_opt_val(fs, "  map_Ka", mat->ka_txt);
        dump_opt_val(fs, "  map_Kd", mat->kd_txt);
        dump_opt_val(fs, "  map_Ks", mat->ks_txt);
        dump_opt_val(fs, "  map_Kr", mat->kr_txt);
        dump_opt_val(fs, "  map_Kt", mat->kt_txt);
        dump_opt_val(fs, "  map_Ns", mat->ns_txt);
        dump_opt_val(fs, "  map_d", mat->op_txt);
        dump_opt_val(fs, "  map_Ni", mat->ior_txt);
        dump_opt_val(fs, "  map_bump", mat->bump_txt);
        dump_opt_val(fs, "  map_disp", mat->disp_txt);
        dump_opt_val(fs, "  map_norm", mat->norm_txt);
        for (auto&& kv : mat->unknown_props) {
            dump_val(fs, kv.first);
            for (auto&& v : kv.second) {
                dump_val(fs, " ");
                dump_val(fs, v);
            }
            dump_val(fs, "\n");
        }
        dump_val(fs, "\n");
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
        dump_named_val(fs, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto cam : asset->cameras) {
        dump_val(fs, "c ");
        dump_val(fs, cam->name);
        dump_val(fs, " ");
        dump_val(fs, cam->ortho);
        dump_val(fs, " ");
        dump_val(fs, cam->yfov);
        dump_val(fs, " ");
        dump_val(fs, cam->aspect);
        dump_val(fs, " ");
        dump_val(fs, cam->aperture);
        dump_val(fs, " ");
        dump_val(fs, cam->focus);
        dump_val(fs, " ");
        dump_val(fs, cam->frame);
        dump_val(fs, '\n');
    }

    // save envs
    for (auto env : asset->environments) {
        dump_val(fs, "e ");
        dump_val(fs, env->name);
        dump_val(fs, " ");
        dump_val(fs, env->matname);
        dump_val(fs, " ");
        dump_val(fs, env->frame);
        dump_val(fs, '\n');
    }

    // save instances
    for (auto ist : asset->instances) {
        dump_val(fs, "i ");
        dump_val(fs, ist->name);
        dump_val(fs, " ");
        dump_val(fs, ist->objname);
        dump_val(fs, " ");
        dump_val(fs, ist->frame);
        dump_val(fs, '\n');
    }

    // save all vertex data
    for (auto& v : asset->pos) dump_named_val(fs, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            dump_named_val(fs, "vt", vec2f{v.x, 1 - v.y});
    } else {
        for (auto& v : asset->texcoord) dump_named_val(fs, "vt", v);
    }
    for (auto& v : asset->norm) dump_named_val(fs, "vn", v);
    for (auto& v : asset->color) dump_named_val(fs, "vc", v);
    for (auto& v : asset->radius) dump_named_val(fs, "vr", v);

    // save element data
    const char* elem_labels[] = {"", "p", "l", "f", "t"};
    for (auto object : asset->objects) {
        dump_named_val(fs, "o", object->name);
        for (auto& group : object->groups) {
            dump_opt_val(fs, "usemtl", group.matname);
            dump_opt_val(fs, "g", group.groupname);
            if (!group.smoothing) dump_named_val(fs, "s", "off");
            for (auto elem : group.elems) {
                dump_objverts(fs, elem_labels[(int)elem.type], elem.size,
                    group.verts.data() + elem.start);
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
struct vertex_hash {
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
        unordered_map<obj_vertex, int, vertex_hash> vert_map;
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

}  // namespace _impl_obj

// Loads an OBJ
obj_scene* load_obj(const string& filename, bool load_txt, bool skip_missing,
    bool flip_texcoord, bool flip_tr) {
    return _impl_obj::load_obj(
        filename, load_txt, skip_missing, flip_texcoord, flip_tr);
}

// Save an OBJ
void save_obj(const string& filename, const obj_scene* asset, bool save_txt,
    bool skip_missing, bool flip_texcoord, bool flip_tr) {
    _impl_obj::save_obj(
        filename, asset, save_txt, skip_missing, flip_texcoord, flip_tr);
}

// Flattens an scene
obj_mesh* get_mesh(
    const obj_scene* model, const obj_object& oshape, bool facet_non_smooth) {
    return _impl_obj::get_mesh(model, oshape, facet_non_smooth);
}

}  // namespace ygl

#if YGL_GLTF

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR KHRONOS GLTF
// -----------------------------------------------------------------------------
namespace ygl {

namespace _impl_gltf {

// Json alias
using json = nlohmann::json;

// #codegen begin func ---------------------------------------------------------

// Parse error
struct parse_stack {
    vector<string> path = {"glTF"};
    string pathname() {
        auto p = std::string();
        for (auto n : path) p += '/' + n;
        return p;
    }
};

// Parse support function.
template <typename T>
inline void parse(vector<T>& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        parse(v, js[i], err);
        vals[i] = v;
    }
}

// Parse int function.
inline void parse(int& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer()) throw runtime_error("integer expected");
    val = js;
}

// Parse float function.
inline void parse(float& val, const json& js, parse_stack& err) {
    if (!js.is_number()) throw runtime_error("number expected");
    val = js;
}

// Parse bool function.
inline void parse(bool& val, const json& js, parse_stack& err) {
    if (!js.is_boolean()) throw runtime_error("bool expected");
    val = js;
}

// Parse std::string function.
inline void parse(string& val, const json& js, parse_stack& err) {
    if (!js.is_string()) throw runtime_error("string expected");
    val = js;
}

// Parse json function.
inline void parse(json& val, const json& js, parse_stack& err) { val = js; }

// Parse support function.
inline void parse(vec2f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (2 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 2; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(vec3f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (3 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 3; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(vec4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (4 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 4; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(quat4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (4 != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < 4; i++) { parse(vals[i], js[i], err); }
}

// Parse support function.
inline void parse(mat4f& vals, const json& js, parse_stack& err) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (16 != js.size()) throw runtime_error("wrong array size");
    for (auto j = 0; j < 4; j++) {
        for (auto i = 0; i < 4; i++) { parse(vals[j][i], js[j * 4 + i], err); }
    }
}

// Parse support function.
template <typename T>
inline void parse(map<string, T>& vals, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        parse(vals[kv.key()], kv.value(), err);
    }
}

// Parse support function.
template <typename T>
inline void parse_attr(
    T& val, const char* name, const json& js, parse_stack& err) {
    auto iter = js.find(name);
    if (iter == js.end()) return;
    err.path.push_back(name);
    parse(val, *iter, err);
    err.path.pop_back();
}

// Parse id function.
template <typename T>
inline void parse(glTFid<T>& val, const json& js, parse_stack& err) {
    if (!js.is_number_integer()) throw runtime_error("int expected");
    val = glTFid<T>((int)js);
}

// Parses a glTFProperty object
inline void parse(glTFProperty*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFProperty();
#if YGL_GLTFJSON
    parse_attr(val->extensions, "extensions", js, err);
    parse_attr(val->extras, "extras", js, err);
#endif
}

// Parses a glTFChildOfRootProperty object
inline void parse(
    glTFChildOfRootProperty*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFChildOfRootProperty();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->name, "name", js, err);
}
// Parse a glTFAccessorSparseIndicesComponentType enum
inline void parse(glTFAccessorSparseIndicesComponentType& val, const json& js,
    parse_stack& err) {
    static map<int, glTFAccessorSparseIndicesComponentType> table = {
        {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
        {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
        {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFAccessorSparseIndices object
inline void parse(
    glTFAccessorSparseIndices*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAccessorSparseIndices();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("bufferView"))
        throw runtime_error("missing required variable");
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    if (!js.count("componentType"))
        throw runtime_error("missing required variable");
    parse_attr(val->componentType, "componentType", js, err);
}

// Parses a glTFAccessorSparseValues object
inline void parse(
    glTFAccessorSparseValues*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAccessorSparseValues();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("bufferView"))
        throw runtime_error("missing required variable");
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
}

// Parses a glTFAccessorSparse object
inline void parse(glTFAccessorSparse*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAccessorSparse();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("count")) throw runtime_error("missing required variable");
    parse_attr(val->count, "count", js, err);
    if (!js.count("indices")) throw runtime_error("missing required variable");
    parse_attr(val->indices, "indices", js, err);
    if (!js.count("values")) throw runtime_error("missing required variable");
    parse_attr(val->values, "values", js, err);
}
// Parse a glTFAccessorComponentType enum
inline void parse(
    glTFAccessorComponentType& val, const json& js, parse_stack& err) {
    static map<int, glTFAccessorComponentType> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parse a glTFAccessorType enum
inline void parse(glTFAccessorType& val, const json& js, parse_stack& err) {
    static map<string, glTFAccessorType> table = {
        {"SCALAR", glTFAccessorType::Scalar},
        {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3},
        {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2},
        {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFAccessor object
inline void parse(glTFAccessor*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAccessor();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->bufferView, "bufferView", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    if (!js.count("componentType"))
        throw runtime_error("missing required variable");
    parse_attr(val->componentType, "componentType", js, err);
    parse_attr(val->normalized, "normalized", js, err);
    if (!js.count("count")) throw runtime_error("missing required variable");
    parse_attr(val->count, "count", js, err);
    if (!js.count("type")) throw runtime_error("missing required variable");
    parse_attr(val->type, "type", js, err);
    parse_attr(val->max, "max", js, err);
    parse_attr(val->min, "min", js, err);
    parse_attr(val->sparse, "sparse", js, err);
}
// Parse a glTFAnimationChannelTargetPath enum
inline void parse(
    glTFAnimationChannelTargetPath& val, const json& js, parse_stack& err) {
    static map<string, glTFAnimationChannelTargetPath> table = {
        {"translation", glTFAnimationChannelTargetPath::Translation},
        {"rotation", glTFAnimationChannelTargetPath::Rotation},
        {"scale", glTFAnimationChannelTargetPath::Scale},
        {"weights", glTFAnimationChannelTargetPath::Weights},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFAnimationChannelTarget object
inline void parse(
    glTFAnimationChannelTarget*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAnimationChannelTarget();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("node")) throw runtime_error("missing required variable");
    parse_attr(val->node, "node", js, err);
    if (!js.count("path")) throw runtime_error("missing required variable");
    parse_attr(val->path, "path", js, err);
}

// Parses a glTFAnimationChannel object
inline void parse(
    glTFAnimationChannel*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAnimationChannel();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("sampler")) throw runtime_error("missing required variable");
    parse_attr(val->sampler, "sampler", js, err);
    if (!js.count("target")) throw runtime_error("missing required variable");
    parse_attr(val->target, "target", js, err);
}
// Parse a glTFAnimationSamplerInterpolation enum
inline void parse(
    glTFAnimationSamplerInterpolation& val, const json& js, parse_stack& err) {
    static map<string, glTFAnimationSamplerInterpolation> table = {
        {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
        {"STEP", glTFAnimationSamplerInterpolation::Step},
        {"CATMULLROMSPLINE",
            glTFAnimationSamplerInterpolation::CatmullRomSpline},
        {"CUBICSPLINE", glTFAnimationSamplerInterpolation::CubicSpline},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFAnimationSampler object
inline void parse(
    glTFAnimationSampler*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAnimationSampler();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("input")) throw runtime_error("missing required variable");
    parse_attr(val->input, "input", js, err);
    parse_attr(val->interpolation, "interpolation", js, err);
    if (!js.count("output")) throw runtime_error("missing required variable");
    parse_attr(val->output, "output", js, err);
}

// Parses a glTFAnimation object
inline void parse(glTFAnimation*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAnimation();
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("channels")) throw runtime_error("missing required variable");
    parse_attr(val->channels, "channels", js, err);
    if (!js.count("samplers")) throw runtime_error("missing required variable");
    parse_attr(val->samplers, "samplers", js, err);
}

// Parses a glTFAsset object
inline void parse(glTFAsset*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFAsset();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->copyright, "copyright", js, err);
    parse_attr(val->generator, "generator", js, err);
    if (!js.count("version")) throw runtime_error("missing required variable");
    parse_attr(val->version, "version", js, err);
    parse_attr(val->minVersion, "minVersion", js, err);
}

// Parses a glTFBuffer object
inline void parse(glTFBuffer*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFBuffer();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->uri, "uri", js, err);
    if (!js.count("byteLength"))
        throw runtime_error("missing required variable");
    parse_attr(val->byteLength, "byteLength", js, err);
}
// Parse a glTFBufferViewTarget enum
inline void parse(glTFBufferViewTarget& val, const json& js, parse_stack& err) {
    static map<int, glTFBufferViewTarget> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFBufferView object
inline void parse(glTFBufferView*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFBufferView();
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("buffer")) throw runtime_error("missing required variable");
    parse_attr(val->buffer, "buffer", js, err);
    parse_attr(val->byteOffset, "byteOffset", js, err);
    if (!js.count("byteLength"))
        throw runtime_error("missing required variable");
    parse_attr(val->byteLength, "byteLength", js, err);
    parse_attr(val->byteStride, "byteStride", js, err);
    parse_attr(val->target, "target", js, err);
}

// Parses a glTFCameraOrthographic object
inline void parse(
    glTFCameraOrthographic*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFCameraOrthographic();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("xmag")) throw runtime_error("missing required variable");
    parse_attr(val->xmag, "xmag", js, err);
    if (!js.count("ymag")) throw runtime_error("missing required variable");
    parse_attr(val->ymag, "ymag", js, err);
    if (!js.count("zfar")) throw runtime_error("missing required variable");
    parse_attr(val->zfar, "zfar", js, err);
    if (!js.count("znear")) throw runtime_error("missing required variable");
    parse_attr(val->znear, "znear", js, err);
}

// Parses a glTFCameraPerspective object
inline void parse(
    glTFCameraPerspective*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFCameraPerspective();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->aspectRatio, "aspectRatio", js, err);
    if (!js.count("yfov")) throw runtime_error("missing required variable");
    parse_attr(val->yfov, "yfov", js, err);
    parse_attr(val->zfar, "zfar", js, err);
    if (!js.count("znear")) throw runtime_error("missing required variable");
    parse_attr(val->znear, "znear", js, err);
}
// Parse a glTFCameraType enum
inline void parse(glTFCameraType& val, const json& js, parse_stack& err) {
    static map<string, glTFCameraType> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFCamera object
inline void parse(glTFCamera*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFCamera();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->orthographic, "orthographic", js, err);
    parse_attr(val->perspective, "perspective", js, err);
    if (!js.count("type")) throw runtime_error("missing required variable");
    parse_attr(val->type, "type", js, err);
}
// Parse a glTFImageMimeType enum
inline void parse(glTFImageMimeType& val, const json& js, parse_stack& err) {
    static map<string, glTFImageMimeType> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFImage object
inline void parse(glTFImage*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFImage();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->uri, "uri", js, err);
    parse_attr(val->mimeType, "mimeType", js, err);
    parse_attr(val->bufferView, "bufferView", js, err);
}

// Parses a glTFTextureInfo object
inline void parse(glTFTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFTextureInfo();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("index")) throw runtime_error("missing required variable");
    parse_attr(val->index, "index", js, err);
    parse_attr(val->texCoord, "texCoord", js, err);
}

// Parses a glTFTexture object
inline void parse(glTFTexture*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFTexture();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->sampler, "sampler", js, err);
    parse_attr(val->source, "source", js, err);
}

// Parses a glTFMaterialNormalTextureInfo object
inline void parse(
    glTFMaterialNormalTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMaterialNormalTextureInfo();
    parse((glTFTextureInfo*&)val, js, err);
    parse_attr(val->scale, "scale", js, err);
}

// Parses a glTFMaterialOcclusionTextureInfo object
inline void parse(
    glTFMaterialOcclusionTextureInfo*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMaterialOcclusionTextureInfo();
    parse((glTFTextureInfo*&)val, js, err);
    parse_attr(val->strength, "strength", js, err);
}

// Parses a glTFMaterialPbrMetallicRoughness object
inline void parse(
    glTFMaterialPbrMetallicRoughness*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMaterialPbrMetallicRoughness();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->baseColorFactor, "baseColorFactor", js, err);
    parse_attr(val->baseColorTexture, "baseColorTexture", js, err);
    parse_attr(val->metallicFactor, "metallicFactor", js, err);
    parse_attr(val->roughnessFactor, "roughnessFactor", js, err);
    parse_attr(
        val->metallicRoughnessTexture, "metallicRoughnessTexture", js, err);
}

// Parses a glTFMaterialPbrSpecularGlossiness object
inline void parse(
    glTFMaterialPbrSpecularGlossiness*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMaterialPbrSpecularGlossiness();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->diffuseFactor, "diffuseFactor", js, err);
    parse_attr(val->diffuseTexture, "diffuseTexture", js, err);
    parse_attr(val->specularFactor, "specularFactor", js, err);
    parse_attr(val->glossinessFactor, "glossinessFactor", js, err);
    parse_attr(
        val->specularGlossinessTexture, "specularGlossinessTexture", js, err);
}
// Parse a glTFMaterialAlphaMode enum
inline void parse(
    glTFMaterialAlphaMode& val, const json& js, parse_stack& err) {
    static map<string, glTFMaterialAlphaMode> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    auto v = string();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFMaterial object
inline void parse(glTFMaterial*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMaterial();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->pbrMetallicRoughness, "pbrMetallicRoughness", js, err);
    parse_attr(val->normalTexture, "normalTexture", js, err);
    parse_attr(val->occlusionTexture, "occlusionTexture", js, err);
    parse_attr(val->emissiveTexture, "emissiveTexture", js, err);
    parse_attr(val->emissiveFactor, "emissiveFactor", js, err);
    parse_attr(val->alphaMode, "alphaMode", js, err);
    parse_attr(val->alphaCutoff, "alphaCutoff", js, err);
    parse_attr(val->doubleSided, "doubleSided", js, err);
    if (js.count("extensions")) {
        auto& js_ext = js["extensions"];
        parse_attr(val->pbrSpecularGlossiness,
            "KHR_materials_pbrSpecularGlossiness", js_ext, err);
    }
}
// Parse a glTFMeshPrimitiveMode enum
inline void parse(
    glTFMeshPrimitiveMode& val, const json& js, parse_stack& err) {
    static map<int, glTFMeshPrimitiveMode> table = {
        {0, glTFMeshPrimitiveMode::Points},
        {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFMeshPrimitive object
inline void parse(glTFMeshPrimitive*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMeshPrimitive();
    parse((glTFProperty*&)val, js, err);
    if (!js.count("attributes"))
        throw runtime_error("missing required variable");
    parse_attr(val->attributes, "attributes", js, err);
    parse_attr(val->indices, "indices", js, err);
    parse_attr(val->material, "material", js, err);
    parse_attr(val->mode, "mode", js, err);
    parse_attr(val->targets, "targets", js, err);
}

// Parses a glTFMesh object
inline void parse(glTFMesh*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFMesh();
    parse((glTFChildOfRootProperty*&)val, js, err);
    if (!js.count("primitives"))
        throw runtime_error("missing required variable");
    parse_attr(val->primitives, "primitives", js, err);
    parse_attr(val->weights, "weights", js, err);
}

// Parses a glTFNode object
inline void parse(glTFNode*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFNode();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->camera, "camera", js, err);
    parse_attr(val->children, "children", js, err);
    parse_attr(val->skin, "skin", js, err);
    parse_attr(val->matrix, "matrix", js, err);
    parse_attr(val->mesh, "mesh", js, err);
    parse_attr(val->rotation, "rotation", js, err);
    parse_attr(val->scale, "scale", js, err);
    parse_attr(val->translation, "translation", js, err);
    parse_attr(val->weights, "weights", js, err);
}
// Parse a glTFSamplerMagFilter enum
inline void parse(glTFSamplerMagFilter& val, const json& js, parse_stack& err) {
    static map<int, glTFSamplerMagFilter> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parse a glTFSamplerMinFilter enum
inline void parse(glTFSamplerMinFilter& val, const json& js, parse_stack& err) {
    static map<int, glTFSamplerMinFilter> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parse a glTFSamplerWrapS enum
inline void parse(glTFSamplerWrapS& val, const json& js, parse_stack& err) {
    static map<int, glTFSamplerWrapS> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parse a glTFSamplerWrapT enum
inline void parse(glTFSamplerWrapT& val, const json& js, parse_stack& err) {
    static map<int, glTFSamplerWrapT> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    auto v = int();
    parse(v, js, err);
    if (table.find(v) == table.end()) throw runtime_error("bad enum value");
    val = table[v];
}

// Parses a glTFSampler object
inline void parse(glTFSampler*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFSampler();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->magFilter, "magFilter", js, err);
    parse_attr(val->minFilter, "minFilter", js, err);
    parse_attr(val->wrapS, "wrapS", js, err);
    parse_attr(val->wrapT, "wrapT", js, err);
}

// Parses a glTFScene object
inline void parse(glTFScene*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFScene();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->nodes, "nodes", js, err);
}

// Parses a glTFSkin object
inline void parse(glTFSkin*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTFSkin();
    parse((glTFChildOfRootProperty*&)val, js, err);
    parse_attr(val->inverseBindMatrices, "inverseBindMatrices", js, err);
    parse_attr(val->skeleton, "skeleton", js, err);
    if (!js.count("joints")) throw runtime_error("missing required variable");
    parse_attr(val->joints, "joints", js, err);
}

// Parses a glTF object
inline void parse(glTF*& val, const json& js, parse_stack& err) {
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new glTF();
    parse((glTFProperty*&)val, js, err);
    parse_attr(val->extensionsUsed, "extensionsUsed", js, err);
    parse_attr(val->extensionsRequired, "extensionsRequired", js, err);
    parse_attr(val->accessors, "accessors", js, err);
    parse_attr(val->animations, "animations", js, err);
    if (!js.count("asset")) throw runtime_error("missing required variable");
    parse_attr(val->asset, "asset", js, err);
    parse_attr(val->buffers, "buffers", js, err);
    parse_attr(val->bufferViews, "bufferViews", js, err);
    parse_attr(val->cameras, "cameras", js, err);
    parse_attr(val->images, "images", js, err);
    parse_attr(val->materials, "materials", js, err);
    parse_attr(val->meshes, "meshes", js, err);
    parse_attr(val->nodes, "nodes", js, err);
    parse_attr(val->samplers, "samplers", js, err);
    parse_attr(val->scene, "scene", js, err);
    parse_attr(val->scenes, "scenes", js, err);
    parse_attr(val->skins, "skins", js, err);
    parse_attr(val->textures, "textures", js, err);
}

// Dump support function.
template <typename T>
inline void dump(const vector<T>& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) { dump(vals[i], js[i], err); }
}

// Converts int to json.
inline void dump(const int& val, json& js, parse_stack& err) { js = val; }

// Converts float to json.
inline void dump(const float& val, json& js, parse_stack& err) { js = val; }

// Converts bool to json.
inline void dump(const bool& val, json& js, parse_stack& err) { js = val; }

// Converts string to json.
inline void dump(const string& val, json& js, parse_stack& err) { js = val; }

// Converts json to json.
inline void dump(const json& val, json& js, parse_stack& err) { js = val; }

// Dump support function.
inline void dump(const vec2f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 2; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const vec3f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 3; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const vec4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 4; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const quat4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto i = 0; i < 4; i++) { dump(vals[i], js[i], err); }
}

// Dump support function.
inline void dump(const mat4f& vals, json& js, parse_stack& err) {
    js = json::array();
    for (auto j = 0; j < 4; j++) {
        for (auto i = 0; i < 4; i++) { dump(vals[j][i], js[j * 4 + i], err); }
    }
}

// Dump support function.
template <typename T>
inline void dump(const map<string, T>& vals, json& js, parse_stack& err) {
    js = json::object();
    for (auto&& kv : vals) { dump(kv.second, js[kv.first], err); }
}

// Dump support function.
template <typename T>
inline void dump_attr(
    const T& val, const char* name, json& js, parse_stack& err) {
    err.path.push_back(name);
    dump(val, js[name], err);
    err.path.pop_back();
}

// Converts glTFid to json.
template <typename T>
inline void dump(const glTFid<T>& val, json& js, parse_stack& err) {
    js = (int)val;
}

// Converts a glTFProperty object to JSON
inline void dump(const glTFProperty* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
    if (!val->extensions.empty())
        dump_attr(val->extensions, "extensions", js, err);
    if (!val->extras.is_null()) dump_attr(val->extras, "extras", js, err);
#endif
}

// Converts a glTFChildOfRootProperty object to JSON
inline void dump(
    const glTFChildOfRootProperty* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->name != "") dump_attr(val->name, "name", js, err);
}
// Converts a glTFAccessorSparseIndicesComponentType enum to JSON
inline void dump(const glTFAccessorSparseIndicesComponentType& val, json& js,
    parse_stack& err) {
    static map<glTFAccessorSparseIndicesComponentType, int> table = {
        {glTFAccessorSparseIndicesComponentType::UnsignedByte, 5121},
        {glTFAccessorSparseIndicesComponentType::UnsignedShort, 5123},
        {glTFAccessorSparseIndicesComponentType::UnsignedInt, 5125},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessorSparseIndices object to JSON
inline void dump(
    const glTFAccessorSparseIndices* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    dump_attr(val->componentType, "componentType", js, err);
}

// Converts a glTFAccessorSparseValues object to JSON
inline void dump(
    const glTFAccessorSparseValues* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
}

// Converts a glTFAccessorSparse object to JSON
inline void dump(const glTFAccessorSparse* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->count, "count", js, err);
    dump_attr(val->indices, "indices", js, err);
    dump_attr(val->values, "values", js, err);
}
// Converts a glTFAccessorComponentType enum to JSON
inline void dump(
    const glTFAccessorComponentType& val, json& js, parse_stack& err) {
    static map<glTFAccessorComponentType, int> table = {
        {glTFAccessorComponentType::Byte, 5120},
        {glTFAccessorComponentType::UnsignedByte, 5121},
        {glTFAccessorComponentType::Short, 5122},
        {glTFAccessorComponentType::UnsignedShort, 5123},
        {glTFAccessorComponentType::UnsignedInt, 5125},
        {glTFAccessorComponentType::Float, 5126},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessorType enum to JSON
inline void dump(const glTFAccessorType& val, json& js, parse_stack& err) {
    static map<glTFAccessorType, string> table = {
        {glTFAccessorType::Scalar, "SCALAR"},
        {glTFAccessorType::Vec2, "VEC2"},
        {glTFAccessorType::Vec3, "VEC3"},
        {glTFAccessorType::Vec4, "VEC4"},
        {glTFAccessorType::Mat2, "MAT2"},
        {glTFAccessorType::Mat3, "MAT3"},
        {glTFAccessorType::Mat4, "MAT4"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAccessor object to JSON
inline void dump(const glTFAccessor* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->bufferView.is_valid())
        dump_attr(val->bufferView, "bufferView", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    dump_attr(val->componentType, "componentType", js, err);
    if (val->normalized != false)
        dump_attr(val->normalized, "normalized", js, err);
    dump_attr(val->count, "count", js, err);
    dump_attr(val->type, "type", js, err);
    if (!val->max.empty()) dump_attr(val->max, "max", js, err);
    if (!val->min.empty()) dump_attr(val->min, "min", js, err);
    if (val->sparse != nullptr) dump_attr(val->sparse, "sparse", js, err);
}
// Converts a glTFAnimationChannelTargetPath enum to JSON
inline void dump(
    const glTFAnimationChannelTargetPath& val, json& js, parse_stack& err) {
    static map<glTFAnimationChannelTargetPath, string> table = {
        {glTFAnimationChannelTargetPath::Translation, "translation"},
        {glTFAnimationChannelTargetPath::Rotation, "rotation"},
        {glTFAnimationChannelTargetPath::Scale, "scale"},
        {glTFAnimationChannelTargetPath::Weights, "weights"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAnimationChannelTarget object to JSON
inline void dump(
    const glTFAnimationChannelTarget* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->node, "node", js, err);
    dump_attr(val->path, "path", js, err);
}

// Converts a glTFAnimationChannel object to JSON
inline void dump(const glTFAnimationChannel* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->sampler, "sampler", js, err);
    dump_attr(val->target, "target", js, err);
}
// Converts a glTFAnimationSamplerInterpolation enum to JSON
inline void dump(
    const glTFAnimationSamplerInterpolation& val, json& js, parse_stack& err) {
    static map<glTFAnimationSamplerInterpolation, string> table = {
        {glTFAnimationSamplerInterpolation::Linear, "LINEAR"},
        {glTFAnimationSamplerInterpolation::Step, "STEP"},
        {glTFAnimationSamplerInterpolation::CatmullRomSpline,
            "CATMULLROMSPLINE"},
        {glTFAnimationSamplerInterpolation::CubicSpline, "CUBICSPLINE"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFAnimationSampler object to JSON
inline void dump(const glTFAnimationSampler* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->input, "input", js, err);
    if (val->interpolation != glTFAnimationSamplerInterpolation::Linear)
        dump_attr(val->interpolation, "interpolation", js, err);
    dump_attr(val->output, "output", js, err);
}

// Converts a glTFAnimation object to JSON
inline void dump(const glTFAnimation* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->channels, "channels", js, err);
    dump_attr(val->samplers, "samplers", js, err);
}

// Converts a glTFAsset object to JSON
inline void dump(const glTFAsset* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->copyright != "") dump_attr(val->copyright, "copyright", js, err);
    if (val->generator != "") dump_attr(val->generator, "generator", js, err);
    dump_attr(val->version, "version", js, err);
    if (val->minVersion != "")
        dump_attr(val->minVersion, "minVersion", js, err);
}

// Converts a glTFBuffer object to JSON
inline void dump(const glTFBuffer* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->uri != "") dump_attr(val->uri, "uri", js, err);
    dump_attr(val->byteLength, "byteLength", js, err);
}
// Converts a glTFBufferViewTarget enum to JSON
inline void dump(const glTFBufferViewTarget& val, json& js, parse_stack& err) {
    static map<glTFBufferViewTarget, int> table = {
        {glTFBufferViewTarget::ArrayBuffer, 34962},
        {glTFBufferViewTarget::ElementArrayBuffer, 34963},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFBufferView object to JSON
inline void dump(const glTFBufferView* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->buffer, "buffer", js, err);
    if (val->byteOffset != 0) dump_attr(val->byteOffset, "byteOffset", js, err);
    dump_attr(val->byteLength, "byteLength", js, err);
    if (val->byteStride != 0) dump_attr(val->byteStride, "byteStride", js, err);
    if (val->target != glTFBufferViewTarget::NotSet)
        dump_attr(val->target, "target", js, err);
}

// Converts a glTFCameraOrthographic object to JSON
inline void dump(
    const glTFCameraOrthographic* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->xmag, "xmag", js, err);
    dump_attr(val->ymag, "ymag", js, err);
    dump_attr(val->zfar, "zfar", js, err);
    dump_attr(val->znear, "znear", js, err);
}

// Converts a glTFCameraPerspective object to JSON
inline void dump(const glTFCameraPerspective* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->aspectRatio != 0)
        dump_attr(val->aspectRatio, "aspectRatio", js, err);
    dump_attr(val->yfov, "yfov", js, err);
    if (val->zfar != 0) dump_attr(val->zfar, "zfar", js, err);
    dump_attr(val->znear, "znear", js, err);
}
// Converts a glTFCameraType enum to JSON
inline void dump(const glTFCameraType& val, json& js, parse_stack& err) {
    static map<glTFCameraType, string> table = {
        {glTFCameraType::Perspective, "perspective"},
        {glTFCameraType::Orthographic, "orthographic"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFCamera object to JSON
inline void dump(const glTFCamera* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->orthographic != nullptr)
        dump_attr(val->orthographic, "orthographic", js, err);
    if (val->perspective != nullptr)
        dump_attr(val->perspective, "perspective", js, err);
    dump_attr(val->type, "type", js, err);
}
// Converts a glTFImageMimeType enum to JSON
inline void dump(const glTFImageMimeType& val, json& js, parse_stack& err) {
    static map<glTFImageMimeType, string> table = {
        {glTFImageMimeType::ImageJpeg, "image/jpeg"},
        {glTFImageMimeType::ImagePng, "image/png"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFImage object to JSON
inline void dump(const glTFImage* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->uri != "") dump_attr(val->uri, "uri", js, err);
    if (val->mimeType != glTFImageMimeType::NotSet)
        dump_attr(val->mimeType, "mimeType", js, err);
    if (val->bufferView.is_valid())
        dump_attr(val->bufferView, "bufferView", js, err);
}

// Converts a glTFTextureInfo object to JSON
inline void dump(const glTFTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->index, "index", js, err);
    if (val->texCoord != 0) dump_attr(val->texCoord, "texCoord", js, err);
}

// Converts a glTFTexture object to JSON
inline void dump(const glTFTexture* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->sampler.is_valid()) dump_attr(val->sampler, "sampler", js, err);
    if (val->source.is_valid()) dump_attr(val->source, "source", js, err);
}

// Converts a glTFMaterialNormalTextureInfo object to JSON
inline void dump(
    const glTFMaterialNormalTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFTextureInfo*)val, js, err);
    if (val->scale != 1) dump_attr(val->scale, "scale", js, err);
}

// Converts a glTFMaterialOcclusionTextureInfo object to JSON
inline void dump(
    const glTFMaterialOcclusionTextureInfo* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFTextureInfo*)val, js, err);
    if (val->strength != 1) dump_attr(val->strength, "strength", js, err);
}

// Converts a glTFMaterialPbrMetallicRoughness object to JSON
inline void dump(
    const glTFMaterialPbrMetallicRoughness* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->baseColorFactor != vec4f{1, 1, 1, 1})
        dump_attr(val->baseColorFactor, "baseColorFactor", js, err);
    if (val->baseColorTexture != nullptr)
        dump_attr(val->baseColorTexture, "baseColorTexture", js, err);
    if (val->metallicFactor != 1)
        dump_attr(val->metallicFactor, "metallicFactor", js, err);
    if (val->roughnessFactor != 1)
        dump_attr(val->roughnessFactor, "roughnessFactor", js, err);
    if (val->metallicRoughnessTexture != nullptr)
        dump_attr(
            val->metallicRoughnessTexture, "metallicRoughnessTexture", js, err);
}

// Converts a glTFMaterialPbrSpecularGlossiness object to JSON
inline void dump(
    const glTFMaterialPbrSpecularGlossiness* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (val->diffuseFactor != vec4f{1, 1, 1, 1})
        dump_attr(val->diffuseFactor, "diffuseFactor", js, err);
    if (val->diffuseTexture != nullptr)
        dump_attr(val->diffuseTexture, "diffuseTexture", js, err);
    if (val->specularFactor != vec3f{1, 1, 1})
        dump_attr(val->specularFactor, "specularFactor", js, err);
    if (val->glossinessFactor != 1)
        dump_attr(val->glossinessFactor, "glossinessFactor", js, err);
    if (val->specularGlossinessTexture != nullptr)
        dump_attr(val->specularGlossinessTexture, "specularGlossinessTexture",
            js, err);
}
// Converts a glTFMaterialAlphaMode enum to JSON
inline void dump(const glTFMaterialAlphaMode& val, json& js, parse_stack& err) {
    static map<glTFMaterialAlphaMode, string> table = {
        {glTFMaterialAlphaMode::Opaque, "OPAQUE"},
        {glTFMaterialAlphaMode::Mask, "MASK"},
        {glTFMaterialAlphaMode::Blend, "BLEND"},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFMaterial object to JSON
inline void dump(const glTFMaterial* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->pbrMetallicRoughness != nullptr)
        dump_attr(val->pbrMetallicRoughness, "pbrMetallicRoughness", js, err);
    if (val->normalTexture != nullptr)
        dump_attr(val->normalTexture, "normalTexture", js, err);
    if (val->occlusionTexture != nullptr)
        dump_attr(val->occlusionTexture, "occlusionTexture", js, err);
    if (val->emissiveTexture != nullptr)
        dump_attr(val->emissiveTexture, "emissiveTexture", js, err);
    if (val->emissiveFactor != vec3f{0, 0, 0})
        dump_attr(val->emissiveFactor, "emissiveFactor", js, err);
    if (val->alphaMode != glTFMaterialAlphaMode::Opaque)
        dump_attr(val->alphaMode, "alphaMode", js, err);
    if (val->alphaCutoff != 0.5)
        dump_attr(val->alphaCutoff, "alphaCutoff", js, err);
    if (val->doubleSided != false)
        dump_attr(val->doubleSided, "doubleSided", js, err);

    if (val->pbrSpecularGlossiness != nullptr) {
        auto& js_ext = js["extensions"];
        dump_attr(val->pbrSpecularGlossiness,
            "KHR_materials_pbrSpecularGlossiness", js_ext, err);
    }
}
// Converts a glTFMeshPrimitiveMode enum to JSON
inline void dump(const glTFMeshPrimitiveMode& val, json& js, parse_stack& err) {
    static map<glTFMeshPrimitiveMode, int> table = {
        {glTFMeshPrimitiveMode::Points, 0},
        {glTFMeshPrimitiveMode::Lines, 1},
        {glTFMeshPrimitiveMode::LineLoop, 2},
        {glTFMeshPrimitiveMode::LineStrip, 3},
        {glTFMeshPrimitiveMode::Triangles, 4},
        {glTFMeshPrimitiveMode::TriangleStrip, 5},
        {glTFMeshPrimitiveMode::TriangleFan, 6},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFMeshPrimitive object to JSON
inline void dump(const glTFMeshPrimitive* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    dump_attr(val->attributes, "attributes", js, err);
    if (val->indices.is_valid()) dump_attr(val->indices, "indices", js, err);
    if (val->material.is_valid()) dump_attr(val->material, "material", js, err);
    if (val->mode != glTFMeshPrimitiveMode::Triangles)
        dump_attr(val->mode, "mode", js, err);
    if (!val->targets.empty()) dump_attr(val->targets, "targets", js, err);
}

// Converts a glTFMesh object to JSON
inline void dump(const glTFMesh* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    dump_attr(val->primitives, "primitives", js, err);
    if (!val->weights.empty()) dump_attr(val->weights, "weights", js, err);
}

// Converts a glTFNode object to JSON
inline void dump(const glTFNode* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->camera.is_valid()) dump_attr(val->camera, "camera", js, err);
    if (!val->children.empty()) dump_attr(val->children, "children", js, err);
    if (val->skin.is_valid()) dump_attr(val->skin, "skin", js, err);
    if (val->matrix !=
        mat4f{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}})
        dump_attr(val->matrix, "matrix", js, err);
    if (val->mesh.is_valid()) dump_attr(val->mesh, "mesh", js, err);
    if (val->rotation != quat4f{0, 0, 0, 1})
        dump_attr(val->rotation, "rotation", js, err);
    if (val->scale != vec3f{1, 1, 1}) dump_attr(val->scale, "scale", js, err);
    if (val->translation != vec3f{0, 0, 0})
        dump_attr(val->translation, "translation", js, err);
    if (!val->weights.empty()) dump_attr(val->weights, "weights", js, err);
}
// Converts a glTFSamplerMagFilter enum to JSON
inline void dump(const glTFSamplerMagFilter& val, json& js, parse_stack& err) {
    static map<glTFSamplerMagFilter, int> table = {
        {glTFSamplerMagFilter::Nearest, 9728},
        {glTFSamplerMagFilter::Linear, 9729},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerMinFilter enum to JSON
inline void dump(const glTFSamplerMinFilter& val, json& js, parse_stack& err) {
    static map<glTFSamplerMinFilter, int> table = {
        {glTFSamplerMinFilter::Nearest, 9728},
        {glTFSamplerMinFilter::Linear, 9729},
        {glTFSamplerMinFilter::NearestMipmapNearest, 9984},
        {glTFSamplerMinFilter::LinearMipmapNearest, 9985},
        {glTFSamplerMinFilter::NearestMipmapLinear, 9986},
        {glTFSamplerMinFilter::LinearMipmapLinear, 9987},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerWrapS enum to JSON
inline void dump(const glTFSamplerWrapS& val, json& js, parse_stack& err) {
    static map<glTFSamplerWrapS, int> table = {
        {glTFSamplerWrapS::ClampToEdge, 33071},
        {glTFSamplerWrapS::MirroredRepeat, 33648},
        {glTFSamplerWrapS::Repeat, 10497},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSamplerWrapT enum to JSON
inline void dump(const glTFSamplerWrapT& val, json& js, parse_stack& err) {
    static map<glTFSamplerWrapT, int> table = {
        {glTFSamplerWrapT::ClampToEdge, 33071},
        {glTFSamplerWrapT::MirroredRepeat, 33648},
        {glTFSamplerWrapT::Repeat, 10497},
    };
    auto v = table.at(val);
    dump(v, js, err);
}

// Converts a glTFSampler object to JSON
inline void dump(const glTFSampler* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->magFilter != glTFSamplerMagFilter::NotSet)
        dump_attr(val->magFilter, "magFilter", js, err);
    if (val->minFilter != glTFSamplerMinFilter::NotSet)
        dump_attr(val->minFilter, "minFilter", js, err);
    if (val->wrapS != glTFSamplerWrapS::Repeat)
        dump_attr(val->wrapS, "wrapS", js, err);
    if (val->wrapT != glTFSamplerWrapT::Repeat)
        dump_attr(val->wrapT, "wrapT", js, err);
}

// Converts a glTFScene object to JSON
inline void dump(const glTFScene* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (!val->nodes.empty()) dump_attr(val->nodes, "nodes", js, err);
}

// Converts a glTFSkin object to JSON
inline void dump(const glTFSkin* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFChildOfRootProperty*)val, js, err);
    if (val->inverseBindMatrices.is_valid())
        dump_attr(val->inverseBindMatrices, "inverseBindMatrices", js, err);
    if (val->skeleton.is_valid()) dump_attr(val->skeleton, "skeleton", js, err);
    dump_attr(val->joints, "joints", js, err);
}

// Converts a glTF object to JSON
inline void dump(const glTF* val, json& js, parse_stack& err) {
    if (!js.is_object()) js = json::object();
    dump((const glTFProperty*)val, js, err);
    if (!val->extensionsUsed.empty())
        dump_attr(val->extensionsUsed, "extensionsUsed", js, err);
    if (!val->extensionsRequired.empty())
        dump_attr(val->extensionsRequired, "extensionsRequired", js, err);
    if (!val->accessors.empty())
        dump_attr(val->accessors, "accessors", js, err);
    if (!val->animations.empty())
        dump_attr(val->animations, "animations", js, err);
    dump_attr(val->asset, "asset", js, err);
    if (!val->buffers.empty()) dump_attr(val->buffers, "buffers", js, err);
    if (!val->bufferViews.empty())
        dump_attr(val->bufferViews, "bufferViews", js, err);
    if (!val->cameras.empty()) dump_attr(val->cameras, "cameras", js, err);
    if (!val->images.empty()) dump_attr(val->images, "images", js, err);
    if (!val->materials.empty())
        dump_attr(val->materials, "materials", js, err);
    if (!val->meshes.empty()) dump_attr(val->meshes, "meshes", js, err);
    if (!val->nodes.empty()) dump_attr(val->nodes, "nodes", js, err);
    if (!val->samplers.empty()) dump_attr(val->samplers, "samplers", js, err);
    if (val->scene.is_valid()) dump_attr(val->scene, "scene", js, err);
    if (!val->scenes.empty()) dump_attr(val->scenes, "scenes", js, err);
    if (!val->skins.empty()) dump_attr(val->skins, "skins", js, err);
    if (!val->textures.empty()) dump_attr(val->textures, "textures", js, err);
}
// #codegen end func

// Get directory name (including '/').
inline string _get_dirname(const string& filename) {
    auto pos = filename.rfind('/');
    if (pos == string::npos) pos = filename.rfind('\\');
    if (pos == string::npos) return "";
    return filename.substr(0, pos + 1);
}

// Get extension name
static inline string _get_extension(const string& filename) {
    auto pos = filename.rfind(".");
    if (pos == string::npos) return "";
    return filename.substr(pos);
}

// Get base name.
inline string _get_basename(const string& filename) {
    auto dirname = _get_dirname(filename);
    auto extension = _get_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

/// Encode in base64
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

/// Decode from base64
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

// Fix path
inline string _fix_path(const string& path_) {
    auto path = path_;
    for (auto& c : path)
        if (c == '\\') c = '/';
    return path;
}

// Load a binary file in memory
// http://stackoverflow.com/questions/116038/what-is-the-best-way-to-read-an-entire-file-into-a-stdstring-in-c
vector<unsigned char> load_binfile(const string& filename, bool skip_missing) {
    std::ifstream ifs(
        filename.c_str(), std::ios::in | std::ios::binary | std::ios::ate);
    if (!ifs) {
        if (skip_missing) return {};
        throw runtime_error("could not open file " + filename);
    }
    std::ifstream::pos_type fileSize = ifs.tellg();
    ifs.seekg(0, std::ios::beg);
    vector<unsigned char> bytes(fileSize);
    ifs.read((char*)&bytes[0], fileSize);
    return bytes;
}

// Saves text.
void save_textfile(const string& filename, const string& txt) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw runtime_error("cannot write file " + filename);
    fwrite(txt.c_str(), 1, (int)txt.size(), f);
    fclose(f);
}

// Saves binary.
void save_binfile(const string& filename, const vector<unsigned char>& bin,
    bool skip_missing) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f && !skip_missing)
        throw runtime_error("cannot write file " + filename);
    fwrite(bin.data(), 1, (int)bin.size(), f);
    fclose(f);
}

// Check if a string starts with a prefix
static inline bool startsiwith(const string& str, const string& prefix) {
    if (str.length() < prefix.length()) return false;
    return str.substr(0, prefix.length()) == prefix;
}

// Load buffer data.
void load_buffers(glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
        if (buffer->uri == "") continue;
        if (startsiwith(buffer->uri, "data:")) {
            // assume it is base64 and find ','
            auto pos = buffer->uri.find(',');
            if (pos == buffer->uri.npos) {
                if (skip_missing) continue;
                throw runtime_error("could not decode base64 data");
            }
            // decode
            auto data = base64_decode(buffer->uri.substr(pos + 1));
            buffer->data = vector<unsigned char>((unsigned char*)data.c_str(),
                (unsigned char*)data.c_str() + data.length());
        } else {
            buffer->data =
                load_binfile(_fix_path(dirname + buffer->uri), skip_missing);
            if (buffer->data.empty()) {
                if (skip_missing) continue;
                throw runtime_error("could not load binary file " +
                                    _fix_path(dirname + buffer->uri));
            }
        }
        if (buffer->byteLength != buffer->data.size()) {
            if (skip_missing) continue;
            throw runtime_error("mismatched buffer size");
        }
    }
}

// Loads images.
void load_images(glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto image : gltf->images) {
        image->data = image_data();
        auto filename = string();
#if YGL_IMAGEIO
        if (image->bufferView || startsiwith(image->uri, "data:")) {
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
            filename = _fix_path(dirname + image->uri);
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
    auto stack = parse_stack();
    auto gltf_ = gltf.get();
    try {
        parse(gltf_, js, stack);
    } catch (const exception& e) {
        throw runtime_error("error parsing gltf at " + stack.pathname() +
                            " with error " + string(e.what()));
    }

    // load external resources
    auto dirname = _get_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

    // done
    return gltf.release();
}

// Save buffer data.
void save_buffers(const glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
        if (startsiwith(buffer->uri, "data:")) {
            if (skip_missing) continue;
            throw runtime_error("saving of embedded data not supported");
        }
        save_binfile(dirname + buffer->uri, buffer->data, skip_missing);
    }
}

// Save images.
void save_images(const glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto image : gltf->images) {
        if (startsiwith(image->uri, "data:")) {
            if (skip_missing) continue;
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
            ok = save_imagef(filename, image->data.width, image->data.height,
                image->data.ncomp, image->data.dataf.data());
        }
#endif
        if (!ok) {
            if (skip_missing) continue;
            throw runtime_error("cannot save image " + filename);
        }
    }
}

// Saves a gltf.
void save_gltf(
    const string& filename, const glTF* gltf, bool save_bin, bool save_image) {
    // dumps json
    auto js = json();
    auto stack = parse_stack();
    dump(gltf, js, stack);

    // save json
    save_textfile(filename, js.dump(2));

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname, false);
    if (save_image) save_images(gltf, dirname, false);
}

// reading shortcut
template <typename T>
inline void read(FILE* f, T* v, int count) {
    if (fread(v, sizeof(T), count, f) != count)
        throw runtime_error("could not read binary file");
}

// writing shortcut
template <typename T>
inline void fwrite(FILE* f, const T* v, int count) {
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
    read(f, &magic, 1);
    if (magic != 0x46546C67) throw runtime_error("corrupted glb format");

    // read version
    uint32_t version;
    read(f, &version, 1);
    if (version != 1 && version != 2)
        throw runtime_error("unsupported glb version");

    // read length
    uint32_t length;
    read(f, &length, 1);

    // data
    auto json_bytes = vector<char>();
    auto buffer_bytes = vector<unsigned char>();
    uint32_t buffer_length = 0;

    if (version == 1) {
        // read content length and format
        uint32_t json_length, json_format;
        read(f, &json_length, 1);
        read(f, &json_format, 1);

        // read json bytes
        json_bytes.resize(json_length);
        read(f, json_bytes.data(), json_length);

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(length - json_length - 20);
            read(f, buffer_bytes.data(), (int)buffer_bytes.size());
            buffer_length = (int)buffer_bytes.size();
        }
    }

    if (version == 2) {
        // read content length and format
        uint32_t json_length, json_format;
        read(f, &json_length, 1);
        read(f, &json_format, 1);
        if (json_format != 0x4E4F534A) {
            throw runtime_error("corrupt binary format");
            return nullptr;
        }

        // read json bytes
        json_bytes.resize(json_length);
        read(f, json_bytes.data(), (int)json_bytes.size());

        // read content length and format
        uint32_t buffer_format;
        read(f, &buffer_length, 1);
        read(f, &buffer_format, 1);
        if (buffer_format != 0x004E4942)
            throw runtime_error("corrupt binary format");

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(buffer_length);
            read(f, buffer_bytes.data(), (int)buffer_bytes.size());
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
    auto stack = parse_stack();
    auto gltf_ = gltf.get();
    try {
        parse(gltf_, js, stack);
    } catch (const exception& e) {
        throw runtime_error("cannot parse gltf json with error at " +
                            stack.pathname() + string(" with error ") +
                            e.what());
        return nullptr;
    }

    // fix internal buffer
    auto buffer = gltf->buffers.at(0);
    buffer->byteLength = buffer_length;
    if (version == 2) buffer->uri = "";
    if (load_bin) { buffer->data = buffer_bytes; }

    // load external resources
    auto dirname = _get_dirname(filename);
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
    auto stack = parse_stack();
    dump(gltf, js, stack);

    // fix string
    auto js_str = js.dump(2);
    if (js_str.length() % 4) {
        auto count = js_str.length() % 4;
        for (auto c = 0; c < count; c++) js_str += " ";
    }
    uint32_t json_length = (uint32_t)js_str.size();

    // internal buffer
    auto buffer = gltf->buffers.at(0);
    uint32_t buffer_length = buffer->byteLength;
    if (buffer_length % 4) buffer_length += 4 - buffer_length % 4;

    // write header
    uint32_t magic = 0x46546C67;
    fwrite(f, &magic, 1);
    uint32_t version = 2;
    fwrite(f, &version, 1);
    uint32_t length = 12 + 8 + json_length + 8 + buffer_length;
    fwrite(f, &length, 1);

    // write json
    uint32_t json_type = 0x4E4F534A;
    fwrite(f, &json_length, 1);
    fwrite(f, &json_type, 1);
    fwrite(f, js_str.data(), (int)json_length);

    if (save_bin) {
        uint32_t buffer_type = 0x004E4942;
        fwrite(f, &buffer_length, 1);
        fwrite(f, &buffer_type, 1);
        fwrite(f, buffer->data.data(), (int)buffer->data.size());
        char pad = 0;
        for (auto i = 0; i < buffer_length - buffer->data.size(); i++)
            fwrite(f, &pad, 1);
    }

    // close
    fclose(f);

    // save external resources
    auto dirname = _get_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname, false);
    if (save_image) save_images(gltf, dirname, false);
}

}  // namespace _impl_gltf

// Loads a gltf file from disk
glTF* load_gltf(
    const string& filename, bool load_bin, bool load_img, bool skip_missing) {
    return _impl_gltf::load_gltf(filename, load_bin, load_img, skip_missing);
}

// Loads a binary gltf file from disk
glTF* load_binary_gltf(
    const string& filename, bool load_bin, bool load_img, bool skip_missing) {
    return _impl_gltf::load_binary_gltf(
        filename, load_bin, load_img, skip_missing);
}

// Saves a scene to disk
void save_gltf(
    const string& filename, const glTF* gltf, bool save_bin, bool save_images) {
    _impl_gltf::save_gltf(filename, gltf, save_bin, save_images);
}

// Saves a scene to disk
void save_binary_gltf(
    const string& filename, const glTF* gltf, bool save_bin, bool save_images) {
    _impl_gltf::save_binary_gltf(filename, gltf, save_bin, save_images);
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SIMPLE SCENE
// -----------------------------------------------------------------------------
namespace ygl {

namespace _impl_scn {

// Flattens an scene
inline scene* obj_to_scene(const obj_scene* obj, const load_options& opts) {
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
        mat->mtype = material_type::specular_roughness;
        mat->ke = {omat->ke.x, omat->ke.y, omat->ke.z};
        mat->kd = {omat->kd.x, omat->kd.y, omat->kd.z};
        mat->ks = {omat->ks.x, omat->ks.y, omat->ks.z};
        mat->kr = {omat->kr.x, omat->kr.y, omat->kr.z};
        mat->kt = {omat->kt.x, omat->kt.y, omat->kt.z};
        mat->rs = pow(2 / (omat->ns + 2), 1 / 4.0f);
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
                unordered_map<obj_vertex, int, obj_vertex_hash> vert_map;
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
    }

    // convert envs
    unordered_set<material*> env_mat;
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
    }

    // remove env materials
    for (auto shp : scn->shapes) env_mat.erase(shp->mat);
    for (auto mat : env_mat) {
        auto end =
            std::remove(scn->materials.begin(), scn->materials.end(), mat);
        scn->materials.erase(end, scn->materials.end());
        delete mat;
    }

    // convert instances
    for (auto oist : obj->instances) {
        for (auto shp : omap[oist->objname]) {
            auto ist = new instance();
            ist->name = oist->name;
            ist->shp = shp;
            ist->frame = oist->frame;
            scn->instances.push_back(ist);
        }
    }

    // done
    return scn;
}

// Load an obj scene
inline scene* load_obj_scene(const string& filename, const load_options& opts) {
    auto oscn = unique_ptr<obj_scene>(load_obj(filename, opts.load_textures,
        opts.skip_missing, opts.obj_flip_texcoord, opts.obj_flip_tr));
    auto scn = unique_ptr<scene>(obj_to_scene(oscn.get(), opts));
    return scn.release();
}

// Save an scene
inline obj_scene* scene_to_obj(const scene* scn) {
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
            otxt->dataf.assign((float*)txt->hdr.data(),
                (float*)txt->hdr.data() +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (txt->ldr) {
            otxt->width = txt->ldr.width();
            otxt->height = txt->ldr.height();
            otxt->ncomp = 4;
            otxt->datab.assign((uint8_t*)txt->ldr.data(),
                (uint8_t*)txt->ldr.data() +
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
        switch (mat->mtype) {
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
    for (auto shp : scn->shapes) {
        auto offset = obj_vertex{(int)obj->pos.size(),
            (int)obj->texcoord.size(), (int)obj->norm.size(),
            (int)obj->color.size(), (int)obj->radius.size()};
        for (auto& v : shp->pos) obj->pos.push_back({v.x, v.y, v.z});
        for (auto& v : shp->norm) obj->norm.push_back({v.x, v.y, v.z});
        for (auto& v : shp->texcoord) obj->texcoord.push_back({v.x, v.y});
        for (auto& v : shp->color) obj->color.push_back({v.x, v.y, v.z, v.w});
        for (auto& v : shp->radius) obj->radius.push_back(v);
        auto object = new obj_object();
        object->name = shp->name;
        object->groups.emplace_back();
        auto group = &object->groups.back();
        group->matname = (shp->mat) ? shp->mat->name : "";
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
            group->elems.push_back({(uint32_t)group->verts.size(),
                obj_element_type::face, (uint16_t)((quad.z == quad.w) ? 3 : 4)});
            if(group->elems.back().size == 3) {
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
        obj->objects.emplace_back(object);
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

    // convert instances
    for (auto ist : scn->instances) {
        auto oist = new obj_instance();
        oist->name = ist->name;
        oist->objname = (ist->shp) ? ist->shp->name : "<undefined>";
        oist->frame = ist->frame;
        obj->instances.emplace_back(oist);
    }

    return obj;
}

// Save an obj scene
inline void save_obj_scene(
    const string& filename, const scene* scn, const save_options& opts) {
    auto oscn = unique_ptr<obj_scene>(scene_to_obj(scn));
    save_obj(filename, oscn.get(), opts.save_textures, opts.skip_missing,
        opts.obj_flip_texcoord, opts.obj_flip_tr);
}

#if YGL_GLTF

// Instance gltf cameras and meshes
inline void gltf_node_to_instances(scene* scn, const vector<camera>& cameras,
    const vector<vector<shape*>>& meshes, const glTF* gltf,
    glTFid<glTFNode> nid, const mat4f& xf) {
    auto nde = gltf->get(nid);
    auto xform = xf * node_transform(nde);
    if (nde->camera) {
        auto cam = new camera(cameras[(int)nde->camera]);
        cam->frame = to_frame3f(xform);
        scn->cameras.push_back(cam);
    }
    if (nde->mesh) {
        for (auto shp : meshes[(int)nde->mesh]) {
            auto ist = new instance();
            ist->name = nde->name;
            ist->frame = to_frame3f(xform);
            ist->shp = shp;
            scn->instances.push_back(ist);
        }
    }
    for (auto cid : nde->children)
        gltf_node_to_instances(scn, cameras, meshes, gltf, cid, xform);
}

// Flattens a gltf file into a flattened asset.
inline scene* gltf_to_scene(const glTF* gltf) {
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
            mat->mtype = material_type::metallic_roughness;
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
            mat->mtype = material_type::specular_glossiness;
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

    // instance meshes and cameras
    if (gltf->scene) {
        for (auto nid : gltf->get(gltf->scene)->nodes) {
            gltf_node_to_instances(
                scn, cameras, meshes, gltf, nid, identity_mat4f);
        }
    } else if (!gltf->nodes.empty()) {
        // set up node children and root nodes
        auto is_root = vector<bool>(gltf->nodes.size(), true);
        for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
            for (auto cid : gltf->get(glTFid<glTFNode>(nid))->children)
                is_root[(int)cid] = false;
        }
        for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
            if (!is_root[nid]) continue;
            gltf_node_to_instances(scn, cameras, meshes, gltf,
                glTFid<glTFNode>(nid), identity_mat4f);
        }
    }

    return scn;
}

// Load an gltf scene
inline scene* load_gltf_scene(
    const string& filename, const load_options& opts) {
    auto gscn = unique_ptr<glTF>(
        load_gltf(filename, true, opts.load_textures, opts.skip_missing));
    auto scn = unique_ptr<scene>(gltf_to_scene(gscn.get()));
    if (!scn) {
        throw runtime_error("could not convert gltf scene");
        return nullptr;
    }
    return scn.release();
}

// Unflattnes gltf
inline glTF* scene_to_gltf(
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
            gimg->data.dataf.assign((float*)txt->hdr.data(),
                (float*)txt->hdr.data() +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (txt->ldr) {
            gimg->data.width = txt->ldr.width();
            gimg->data.height = txt->ldr.height();
            gimg->data.ncomp = 4;
            gimg->data.datab.assign((uint8_t*)txt->ldr.data(),
                (uint8_t*)txt->ldr.data() +
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
        switch (mat->mtype) {
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
                    auto bbox = make_bbox(count, (float*)data);
                    accessor->min = {bbox.min};
                    accessor->max = {bbox.max};
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
        auto gbuffer = add_opt_buffer(shp->path);
        auto gmesh = new glTFMesh();
        gmesh->name = shp->name;
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
        gmesh->primitives.push_back(gprim);
        gltf->meshes.push_back(gmesh);
    }

    // instances
    for (auto ist : scn->instances) {
        auto gnode = new glTFNode();
        gnode->name = ist->name;
        gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, ist->shp));
        gnode->matrix = to_mat4f(ist->frame);
        gltf->nodes.push_back(gnode);
    }

    // cameras
    for (auto cam : scn->cameras) {
        auto gnode = new glTFNode();
        gnode->name = cam->name;
        gnode->camera = glTFid<glTFCamera>(index(scn->cameras, cam));
        gnode->matrix = to_mat4f(cam->frame);
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

    // done
    return gltf.release();
}

// Save a gltf scene
inline void save_gltf_scene(
    const string& filename, const scene* scn, const save_options& opts) {
    auto buffer_uri = path_basename(filename) + ".bin";
    auto gscn = unique_ptr<glTF>(
        scene_to_gltf(scn, buffer_uri, opts.gltf_separate_buffers));
    save_gltf(filename, gscn.get(), true, opts.save_textures);
}

#endif

// Load a scene
inline scene* load_scene(const string& filename, const load_options& opts) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ") return load_obj_scene(filename, opts);
#if YGL_GLTF
    if (ext == ".gltf" || ext == ".GLTF")
        return load_gltf_scene(filename, opts);
#endif
    throw runtime_error("unsupported extension " + ext);
    return nullptr;
}

// Save a scene
inline void save_scene(
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

// Add missing values and elements
inline void add_elements(scene* scn, const add_elements_options& opts) {
    if (opts.smooth_normals) {
        for (auto shp : scn->shapes) {
            if (!shp->norm.empty()) continue;
            shp->norm.resize(shp->pos.size(), {0, 0, 1});
            if (!shp->lines.empty() || !shp->triangles.empty() ||
                !shp->quads.empty()) {
                shp->norm = compute_normals(
                    shp->lines, shp->triangles, shp->quads, shp->pos);
            }
        }
    }

    if (opts.tangent_space) {
        for (auto shp : scn->shapes) {
            if (!shp->tangsp.empty() || shp->triangles.empty() ||
                shp->texcoord.empty() || (shp->mat))
                continue;
            shp->tangsp = compute_tangent_frames(
                shp->triangles, shp->pos, shp->norm, shp->texcoord);
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

    if (opts.default_camera && scn->cameras.empty()) {
        update_bounds(scn);
        auto bbox = scn->bbox;
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
        cam->frame = lookat_frame3f(from, to, up);
        cam->ortho = false;
        cam->aspect = 16.0f / 9.0f;
        cam->yfov = 2 * atanf(0.5f);
        cam->aperture = 0;
        cam->focus = length(to - from);
        scn->cameras.push_back(cam);
    }

    // default environment
    if (opts.default_environment && scn->environments.empty()) {
        auto env = new environment();
        env->name = "default_environment";
        scn->environments.push_back(env);
    }
}

// Merge scene into one another
inline void merge_into(scene* merge_into, scene* merge_from) {
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

// Initialize the lights
inline void update_lights(scene* scn, bool point_only) {
    for (auto lgt : scn->lights) delete lgt;
    scn->lights.clear();

    for (auto ist : scn->instances) {
        if(!ist->shp->mat) continue;
        if (ist->shp->mat->ke == zero3f) continue;
        if (point_only && ist->shp->points.empty()) continue;
        auto lgt = new light();
        lgt->ist = ist;
        if (point_only) continue;
        auto shp = ist->shp;
        if (shp->elem_cdf.empty()) {
            if (!shp->points.empty()) {
                shp->elem_cdf = sample_points_cdf(shp->points.size());
            } else if (!ist->shp->lines.empty()) {
                shp->elem_cdf = sample_lines_cdf(shp->lines, shp->pos);
            } else if (!shp->triangles.empty()) {
                shp->elem_cdf = sample_triangles_cdf(shp->triangles, shp->pos);
            }
        }
        scn->lights.push_back(lgt);
    }

    for (auto env : scn->environments) {
        if (point_only) continue;
        if (env->ke == zero3f) continue;
        auto lgt = new light();
        lgt->env = env;
        scn->lights.push_back(lgt);
    }
}

// Print scene info (call update bounds bes before)
inline void print_info(const scene* scn) {
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

    auto bbox = scn->bbox;
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

}  // namespace _impl_scn

// Load a scene
scene* load_scene(const string& filename, const load_options& opts) {
    return _impl_scn::load_scene(filename, opts);
}

// Save a scene
void save_scene(
    const string& filename, const scene* scn, const save_options& opts) {
    _impl_scn::save_scene(filename, scn, opts);
}

// Add missing values and elements
void add_elements(scene* scn, const add_elements_options& opts) {
    _impl_scn::add_elements(scn, opts);
}

// Merge scene into one another
void merge_into(scene* merge_into, scene* merge_from) {
    _impl_scn::merge_into(merge_into, merge_from);
}

// Initialize the lights
void update_lights(scene* scn, bool point_only) {
    _impl_scn::update_lights(scn, point_only);
}

// Print scene info (call update bounds bes before)
void print_info(const scene* scn) { _impl_scn::print_info(scn); }

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvsphere(
    int usteps, int vsteps) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
    }
    return {quads, pos, norm, texcoord};
}

// Make a geodesic sphere.
tuple<vector<vec3i>, vector<vec3f>, vector<vec3f>> make_geodesicsphere(
    int level) {
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
    for (auto l = 0; l < level - 2; l++) {
        vector<vec2i> _lines;
        vector<vec4i> _quads;
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(_lines, triangles, _quads, edges, faces) =
            subdivide_elems({}, triangles, {}, (int)pos.size());
        pos = subdivide_vert(pos, edges, faces);
    }
    for (auto& p : pos) p = normalize(p);
    return {triangles, pos, pos};
}

// Make a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhemisphere(int usteps, int vsteps) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pif * uv.x, pif * 0.5f * (1 - uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
    }
    return {quads, pos, norm, texcoord};
}

// Make an inside-out sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflippedsphere(int usteps, int vsteps) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pif * uv.x, pif * uv.y};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
        texcoord[i] = {uv.x, 1 - uv.y};
    }
    return {quads, pos, norm, texcoord};
}

// Make an inside-out hemisphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflippedhemisphere(int usteps, int vsteps) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pif * uv.x, pif * (0.5f + 0.5f * uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
        texcoord[i] = {uv.x, 1 - uv.y};
    }
    return {quads, pos, norm, texcoord};
}

// Make a quad.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvquad(
    int usteps, int vsteps) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
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
tuple<vector<vec4i>, vector<vec3f>> make_cube() {
    static auto cube_pos =
        vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1}, {+1, -1, -1},
            {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    static auto cube_quads = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_quad_uv = vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    return {cube_quads, cube_pos};
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec2f>>
make_fvcube() {
    static auto cube_pos =
        vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1}, {+1, -1, -1},
            {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    static auto cube_qpos = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_texcoord = vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1},
        {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
        {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0},
        {1, 1}, {0, 1}};
    static auto cube_qtexcoord = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
    return {cube_qpos, cube_pos, cube_qtexcoord, cube_texcoord};
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
tuple<vector<vec4i>, vector<vec3f>> make_suzanne() {
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
    auto quads = suzanne_quads;
    quads.reserve(suzanne_quads.size() + suzanne_triangles.size());
    for (auto& t : suzanne_triangles) { quads.push_back({t.x, t.y, t.z, t.z}); }
    return {quads, suzanne_pos};
}

// Make a cube with uv. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvcube(
    int usteps, int vsteps) {
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
        make_uvquad(usteps, vsteps);
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
make_uvspherecube(int usteps, int vsteps) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvcube(usteps, vsteps);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = normalize(pos[i]);
        norm[i] = normalize(pos[i]);
    }
    return {quads, pos, norm, texcoord};
}

// Make a cube than stretch it towards a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherizedcube(int usteps, int vsteps, float radius) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvcube(usteps, vsteps);
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
make_uvflipcapsphere(int usteps, int vsteps, float radius) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvsphere(usteps, vsteps);
    for (auto i = 0; i < pos.size(); i++) {
        if (pos[i].z > radius) {
            pos[i].z = 2 * radius - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        } else if (pos[i].z < -radius) {
            pos[i].z = -2 * radius - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        }
    }
    return {quads, pos, norm, texcoord};
}

// Make a cutout sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvcutsphere(int usteps, int vsteps, float radius) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto p = 1 - acos(radius) / pif;
        auto a = vec2f{2 * pif * uv.x, pif * (1 - p * uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
    }
    return {quads, pos, norm, texcoord};
}

// Make a flipped and cut sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflippedcutsphere(int usteps, int vsteps, float radius) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto p = 1 - acos(radius) / pif;
        auto a = vec2f{2 * pif * uv.x, pif * ((1 - p) + p * uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
        texcoord[i] = {uv.x, (1 - uv.y)};
    }
    return {quads, pos, norm, texcoord};
}

/// Make a hair ball around a shape
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>>
make_hair(int num, int usteps, const vec2f& len, const vec2f& rad,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, uint32_t seed) {
    vector<vec3f> bpos;
    vector<vec3f> bnorm;
    vector<vec2f> btexcoord;
    tie(bpos, bnorm, btexcoord) =
        sample_triangles_points(striangles + convert_quads_to_triangles(squads),
            spos, snorm, stexcoord, num, seed);

    auto rng = init_rng(seed, 3);
    auto blen = vector<float>(bpos.size());
    for (auto& l : blen) l = lerp(len.x, len.y, next_rand1f(rng));

    auto cidx = vector<int>();
    if (clump.x > 0) {
        for (auto bidx : range(bpos.size())) {
            cidx += 0;
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

    auto lines = vector<vec2i>();
    auto texcoord = vector<vec2f>();
    tie(lines, texcoord) = make_lines(num, usteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    auto radius = vector<float>(texcoord.size());
    for (auto i : range(texcoord.size())) {
        auto u = texcoord[i].x;
        auto bidx = i / (usteps + 1);
        pos[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
        norm[i] = bnorm[bidx];
        radius[i] = lerp(rad.x, rad.y, u);
        if (clump.x > 0) {
            pos[i] = lerp(pos[i], pos[i + (cidx[bidx] - bidx) * (usteps + 1)],
                u * clump.x);
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
        cam->frame = lookat_frame3f(from, to, {0, 1, 0});
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
        ist->frame = {rotation_mat3f(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                          rotation_mat3f(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                          rotation_mat3f(vec3f{1, 0, 0}, rot[0] * pif / 180),
            pos};
        return ist;
    };

    auto make_quad = [](string name, material* mat, float scale = 1) {
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        tie(shp->quads, shp->pos, shp->norm, shp->texcoord) = make_uvquad(1, 1);
        for (auto& p : shp->pos) p *= scale;
        return shp;
    };

    auto make_box = [](string name, material* mat, vec3f scale) {
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        tie(shp->quads, shp->pos, shp->norm, shp->texcoord) = make_uvcube(1, 1);
        for (auto& p : shp->pos) p *= scale;
        return shp;
    };

    auto make_material = [](string name, vec3f kd, vec3f ke = {0, 0, 0}) {
        auto mat = new material();
        mat->mtype = material_type::specular_roughness;
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
make_uvhollowcutsphere(int usteps, int vsteps, float radius) {
    auto quads = vector<vec4i>();
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();

    vector<vec3f> mpos, mnorm;
    vector<vec2f> mtexcoord;
    vector<vec4i> mquads;
    vector<vec2i> _aux1;
    vector<vec3i> _aux2;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvcutsphere(usteps, vsteps, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvflippedcutsphere(usteps, vsteps, radius);
    for (auto& p : mpos) p *= radius;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    // dpdu = [- s r s0 s1, s r c0 s1, 0] === [- s0, c0, 0]
    // dpdv = [s c0 s1, s s0 s1, s c1] === [c0 s1, s0 s1, c1]
    // n = [c0 c1, - s0 c1, s1]
    tie(mquads, mtexcoord) = make_uvquads(usteps, vsteps);
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
make_uvhollowcutsphere1(int usteps, int vsteps, float radius) {
    auto quads = vector<vec4i>();
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();

    vector<vec3f> mpos, mnorm;
    vector<vec2f> mtexcoord;
    vector<vec4i> mquads;
    vector<vec2i> _aux1;
    vector<vec3i> _aux2;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvcutsphere(usteps, vsteps, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    for (auto i = (usteps + 1) * vsteps; i < mnorm.size(); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvflippedcutsphere(usteps, vsteps, radius * 1.05f);
    for (auto& p : mpos) p *= 0.8f;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mtexcoord) = make_uvquads(usteps, vsteps / 4);
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
    for (auto i = 0; i < (usteps + 1); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;
    return {quads, pos, norm, texcoord};
}

enum struct test_texture_type {
    none,
    grid,
    checker,
    colored,
    rcolored,
    tgrid,
    bump,
    uv,
    gamma,
    gridn,
    bumpn,
    tgridn,
    noise,
    fbm,
    ridge,
    turbulence,
    gammaf,
    sky1,
    sky2
};

inline const vector<pair<string, test_texture_type>>& test_texture_names() {
    static auto names = vector<pair<string, test_texture_type>>{
        {"none.png", test_texture_type::none},
        {"grid.png", test_texture_type::grid},
        {"checker.png", test_texture_type::checker},
        {"bump.png", test_texture_type::bumpn},
        {"tgrid.png", test_texture_type::tgrid},
        {"colored.png", test_texture_type::colored},
        {"rcolored.png", test_texture_type::rcolored},
        {"uv.png", test_texture_type::uv},
        {"gamma.png", test_texture_type::gamma},
        {"gridn.png", test_texture_type::gridn},
        {"bumpn.png", test_texture_type::bumpn},
        {"tgridn.png", test_texture_type::tgridn},
        {"noise.png", test_texture_type::noise},
        {"fbm.png", test_texture_type::fbm},
        {"ridge.png", test_texture_type::ridge},
        {"turbulence.png", test_texture_type::turbulence},
        {"gammaf.hdr", test_texture_type::gammaf},
        {"sky1.hdr", test_texture_type::sky1},
        {"sky2.hdr", test_texture_type::sky2}};
    return names;
}

texture* add_test_texture(scene* scn, test_texture_type type) {
    if (type == test_texture_type::none) return nullptr;
    auto name = ""s;
    for (auto kv : test_texture_names())
        if (kv.second == type) name = kv.first;
    for (auto txt : scn->textures)
        if (txt->path == name) return txt;
    auto txt = new texture();
    txt->path = name;
    scn->textures += txt;
    switch (type) {
        case test_texture_type::grid: {
            txt->ldr = make_grid_image(512, 512);
        } break;
        case test_texture_type::checker: {
            txt->ldr = make_checker_image(512, 512);
        } break;
        case test_texture_type::colored: {
            txt->ldr = make_uvgrid_image(512, 512);
        } break;
        case test_texture_type::rcolored: {
            txt->ldr = make_recuvgrid_image(512, 512);
        } break;
        case test_texture_type::bump: {
            txt->ldr = make_bumpdimple_image(512, 512, 32);
        } break;
        case test_texture_type::tgrid: {
            txt->ldr = make_grid_image(512, 512, 32);
        } break;
        case test_texture_type::uv: {
            txt->ldr = make_uv_image(512, 512);
        } break;
        case test_texture_type::gamma: {
            txt->ldr = make_gammaramp_image(512, 512);
        } break;
        case test_texture_type::gridn: {
            txt->ldr = bump_to_normal_map(make_grid_image(512, 512), 4);
        } break;
        case test_texture_type::bumpn: {
            txt->ldr =
                bump_to_normal_map(make_bumpdimple_image(512, 512, 32), 4);
        } break;
        case test_texture_type::tgridn: {
            txt->ldr = bump_to_normal_map(make_grid_image(512, 512, 32), 4);
        } break;
        case test_texture_type::noise: {
            txt->ldr = make_noise_image(512, 512, 8);
        } break;
        case test_texture_type::ridge: {
            txt->ldr = make_ridge_image(512, 512, 8);
        } break;
        case test_texture_type::fbm: {
            txt->ldr = make_fbm_image(512, 512, 8);
        } break;
        case test_texture_type::turbulence: {
            txt->ldr = make_turbulence_image(512, 512, 8);
        } break;
        case test_texture_type::gammaf: {
            txt->hdr = make_gammaramp_imagef(512, 512);
        } break;
        case test_texture_type::sky1: {
            txt->hdr = make_sunsky_image(512, pif / 4);
        } break;
        case test_texture_type::sky2: {
            txt->hdr = make_sunsky_image(512, pif / 2);
        } break;
        default: throw runtime_error("bad value");
    }
    return txt;
}

enum struct test_material_type {
    none,
    matte_grid,
    matte_gray,
    matte_green,
    matte_colored,
    matte_uv,
    plastic_red,
    plastic_blue,
    plastic_green,
    plastic_colored,
    plastic_bumped,
    silver_mirror,
    silver_rough,
    gold_mirror,
    gold_rough,
    transparent_red,
    transparent_green,
    transparent_blue,
    light_point,
    light_area,
    light_areal,
    light_arear,
};

inline const vector<pair<string, test_material_type>>& test_material_names() {
    static auto names = vector<pair<string, test_material_type>>{
        {"none", test_material_type::none},
        {"matte_grid", test_material_type::matte_grid},
        {"matte_gray", test_material_type::matte_gray},
        {"matte_green", test_material_type::matte_green},
        {"matte_colored", test_material_type::matte_colored},
        {"matte_uv", test_material_type::matte_uv},
        {"plastic_red", test_material_type::plastic_red},
        {"plastic_blue", test_material_type::plastic_blue},
        {"plastic_green", test_material_type::plastic_green},
        {"plastic_colored", test_material_type::plastic_colored},
        {"plastic_bumped", test_material_type::plastic_bumped},
        {"silver_mirror", test_material_type::silver_mirror},
        {"silver_rough", test_material_type::silver_rough},
        {"gold_mirror", test_material_type::gold_mirror},
        {"gold_rough", test_material_type::gold_rough},
        {"transparent_red", test_material_type::transparent_red},
        {"transparent_green", test_material_type::transparent_green},
        {"transparent_blue", test_material_type::transparent_blue},
        {"light_point", test_material_type::light_point},
        {"light_area", test_material_type::light_area},
        {"light_areal", test_material_type::light_areal},
        {"light_arear", test_material_type::light_arear},
    };

    return names;
}

inline material* add_test_material(scene* scn, test_material_type type) {
    if (type == test_material_type::none) return nullptr;
    auto name = ""s;
    for (auto kv : test_material_names())
        if (kv.second == type) name = kv.first;
    for (auto mat : scn->materials)
        if (mat->name == name) return mat;
    auto mat = new material();
    mat->name = name;
    scn->materials += mat;
    switch (type) {
        case test_material_type::matte_gray: {
            mat->kd = {0.2f, 0.2f, 0.2f};
        } break;
        case test_material_type::matte_green: {
            mat->kd = {0.2f, 0.5f, 0.2f};
        } break;
        case test_material_type::matte_grid: {
            mat->kd = {1, 1, 1};
            mat->kd_txt.txt = add_test_texture(scn, test_texture_type::grid);
        } break;
        case test_material_type::matte_colored: {
            mat->kd = {1, 1, 1};
            mat->kd_txt.txt = add_test_texture(scn, test_texture_type::colored);
        } break;
        case test_material_type::matte_uv: {
            mat->kd = {1, 1, 1};
            mat->kd_txt.txt = add_test_texture(scn, test_texture_type::uv);
        } break;
        case test_material_type::plastic_red: {
            mat->kd = {0.5f, 0.2f, 0.2f};
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = 0.25f;
        } break;
        case test_material_type::plastic_green: {
            mat->kd = {0.2f, 0.5f, 0.2f};
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = 0.1f;
        } break;
        case test_material_type::plastic_blue: {
            mat->kd = {0.2f, 0.2f, 0.5f};
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = 0.05f;
        } break;
        case test_material_type::plastic_colored: {
            mat->kd = {1, 1, 1};
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = 0.25f;
            mat->kd_txt.txt = add_test_texture(scn, test_texture_type::colored);
        } break;
        case test_material_type::plastic_bumped: {
            mat->kd = {1, 1, 1};
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = 0.25f;
            mat->kd_txt.txt = add_test_texture(scn, test_texture_type::colored);
            mat->norm_txt.txt = add_test_texture(scn, test_texture_type::bumpn);
        } break;
        case test_material_type::silver_mirror: {
            mat->ks = {0.5, 0.5, 0.5};
            mat->rs = 0.05f;
        } break;
        case test_material_type::silver_rough: {
            mat->ks = {0.5, 0.5, 0.5};
            mat->rs = 0.25f;
        } break;
        case test_material_type::gold_mirror: {
            mat->ks = {0.66f, 0.45f, 0.34f};
            mat->rs = 0.05f;
        } break;
        case test_material_type::gold_rough: {
            mat->ks = {0.66f, 0.45f, 0.34f};
            mat->rs = 0.25f;
        } break;
        case test_material_type::transparent_red: {
            mat->kd = {0.5f, 0.2f, 0.2f};
            mat->op = 0.9f;
        } break;
        case test_material_type::transparent_green: {
            mat->kd = {0.5f, 0.2f, 0.2f};
            mat->op = 0.5f;
        } break;
        case test_material_type::transparent_blue: {
            mat->kd = {0.5f, 0.2f, 0.2f};
            mat->op = 0.92f;
        } break;
        case test_material_type::light_point: {
            mat->ke = {400, 400, 400};
        } break;
        case test_material_type::light_area: {
            mat->ke = {80, 80, 80};
        } break;
        case test_material_type::light_arear: {
            mat->ke = {80, 80, 80};
        } break;
        case test_material_type::light_areal: {
            mat->ke = {20, 20, 20};
        } break;
        default: throw runtime_error("bad value");
    }
    return mat;
}

enum struct test_shape_type {
    none,
    floor,
    quad,
    cube,
    sphere,
    spherecube,
    spherizedcube,
    flipcapsphere,
    geosphere,
    geospherel,
    geospheref,
    cubep,
    cubes,
    suzanne,
    suzannes,
    cubefv,
    cubefvs,
    quads,
    matball1,
    matball2,
    matballi,
    points,
    lines1,
    lines2,
    lines3,
    linesi,
    plight,
    alight,
    alightl,
    alightr,
};

inline const vector<pair<string, test_shape_type>>& test_shape_names() {
    static auto names = vector<pair<string, test_shape_type>>{
        {"none", test_shape_type::none},
        {"floor", test_shape_type::floor},
        {"quad", test_shape_type::quad},
        {"cube", test_shape_type::cube},
        {"sphere", test_shape_type::sphere},
        {"spherecube", test_shape_type::spherecube},
        {"spherizedcube", test_shape_type::spherizedcube},
        {"flipcapsphere", test_shape_type::flipcapsphere},
        {"geosphere", test_shape_type::geosphere},
        {"geospherel", test_shape_type::geospherel},
        {"geospheref", test_shape_type::geospheref},
        {"cubep", test_shape_type::cubep},
        {"cubes", test_shape_type::cubes},
        {"suzanne", test_shape_type::suzanne},
        {"suzannes", test_shape_type::suzannes},
        {"cubefv", test_shape_type::cubefv},
        {"cubefvs", test_shape_type::cubefvs},
        {"quads", test_shape_type::quads},
        {"matball1", test_shape_type::matball1},
        {"matball2", test_shape_type::matball2},
        {"matballi", test_shape_type::matballi},
        {"points", test_shape_type::points},
        {"lines1", test_shape_type::lines1},
        {"lines2", test_shape_type::lines2},
        {"lines3", test_shape_type::lines3},
        {"linesi", test_shape_type::linesi},
        {"plight", test_shape_type::plight},
        {"alight", test_shape_type::alight},
        {"alightl", test_shape_type::alightl},
        {"alightr", test_shape_type::alightr},
    };
    return names;
}

inline shape* add_test_shape(
    scene* scn, test_shape_type stype, test_material_type mtype) {
    if (stype == test_shape_type::none) return nullptr;
    auto name = ""s;
    for (auto kv : test_shape_names())
        if (kv.second == stype) name += kv.first;
    name += "_";
    for (auto kv : test_material_names())
        if (kv.second == mtype) name += kv.first;
    for (auto shp : scn->shapes)
        if (shp->name == name) return shp;
    auto shp = new shape();
    shp->name = name;
    shp->mat = add_test_material(scn, mtype);
    scn->shapes += shp;
    switch (stype) {
        case test_shape_type::floor: {
            auto level = 6;
            auto usteps = pow2(level), vsteps = pow2(level);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(usteps, vsteps);
            for (auto& p : shp->pos) p = {-p.x, p.z, p.y};
            for (auto& n : shp->norm) n = {n.x, n.z, n.y};
            for (auto& p : shp->pos) p *= 20;
            for (auto& uv : shp->texcoord) uv *= 20;
        } break;
        case test_shape_type::quad: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(1, 1);
        } break;
        case test_shape_type::cube: {
            auto level = 0;
            auto usteps = pow2(level), vsteps = pow2(level);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvcube(usteps, vsteps);
        } break;
        case test_shape_type::sphere: {
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere(usteps, vsteps);
        } break;
        case test_shape_type::spherecube: {
            auto level = 4;
            auto usteps = pow2(level), vsteps = pow2(level);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(usteps, vsteps);
        } break;
        case test_shape_type::spherizedcube: {
            auto level = 4;
            auto usteps = pow2(level), vsteps = pow2(level);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherizedcube(usteps, vsteps, 0.75f);
        } break;
        case test_shape_type::flipcapsphere: {
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvflipcapsphere(usteps, vsteps, 0.75f);
            // for (auto& uv : shp->texcoord) uv *= uvscale;
        } break;
        case test_shape_type::geosphere: {
            auto level = 5;
            tie(shp->triangles, shp->pos, shp->norm) =
                make_geodesicsphere(level);
        } break;
        case test_shape_type::geospheref: {
            auto level = 5;
            tie(shp->triangles, shp->pos, shp->norm) =
                make_geodesicsphere(level);
            facet_shape(shp);
        } break;
        case test_shape_type::geospherel: {
            auto level = 4;
            tie(shp->triangles, shp->pos, shp->norm) =
                make_geodesicsphere(level);
            facet_shape(shp);
        } break;
        case test_shape_type::cubep: {
            tie(shp->quads, shp->pos) = make_cube();
        } break;
        case test_shape_type::cubes: {
            tie(shp->quads, shp->pos) = make_cube();
            for (auto i = 0; i < 4; i++) subdivide_shape(shp, true);
        } break;
        case test_shape_type::suzanne: {
            tie(shp->quads, shp->pos) = make_suzanne();
        } break;
        case test_shape_type::suzannes: {
            tie(shp->quads, shp->pos) = make_suzanne();
            for (auto i = 0; i < 2; i++) subdivide_shape(shp, true);
        } break;
        case test_shape_type::cubefv: {
            tie(shp->quads_pos, shp->pos, shp->quads_texcoord, shp->texcoord) =
                make_fvcube();
        } break;
        case test_shape_type::cubefvs: {
            tie(shp->quads_pos, shp->pos, shp->quads_texcoord, shp->texcoord) =
                make_fvcube();
            for (auto l = 0; l < 4; l++) subdivide_shape(shp, true);
        } break;
        case test_shape_type::quads: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(1, 1);
            for (auto i = 0; i < 4; i++) subdivide_shape(shp, true);
        } break;
        case test_shape_type::matball1: {
            auto level = 6;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvhollowcutsphere(usteps, vsteps, 0.75f);
        } break;
        case test_shape_type::matball2: {
            auto level = 6;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvhollowcutsphere1(usteps, vsteps, 0.75f);
        } break;
        case test_shape_type::matballi: {
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere(usteps, vsteps);
            for (auto& p : shp->pos) p *= 0.8f;
        } break;
        case test_shape_type::points: {
            tie(shp->points, shp->texcoord) = make_points(64 * 64 * 16);
            shp->pos.reserve(shp->texcoord.size());
            shp->norm.resize(shp->texcoord.size(), {0, 0, 1});
            shp->radius.resize(shp->texcoord.size(), 0.0025f);
            auto rn = init_rng(0);
            for (auto i = 0; i < shp->texcoord.size(); i++) {
                shp->pos += vec3f{-1 + 2 * next_rand1f(rn),
                    -1 + 2 * next_rand1f(rn), -1 + 2 * next_rand1f(rn)};
            }
        } break;
        case test_shape_type::lines1: {
            auto nhairs = 65536;
            // auto nhairs = 32768;
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(usteps, vsteps);
            tie(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius) =
                make_hair(nhairs, 4, {0.1f, 0.1f}, {0.001f, 0.0001f}, {},
                    shp->quads, shp->pos, shp->norm, shp->texcoord, {0.5f, 8});
            shp->quads.clear();
        } break;
        case test_shape_type::lines2: {
            auto nhairs = 65536;
            // auto nhairs = 32768;
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(usteps, vsteps);
            tie(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius) =
                make_hair(nhairs, 4, {0.1f, 0.1f}, {0.001f, 0.0001f}, {},
                    shp->quads, shp->pos, shp->norm, shp->texcoord, {},
                    {0.5f, 128});
            shp->quads.clear();
        } break;
        case test_shape_type::lines3: {
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(usteps, vsteps);
            tie(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius) =
                make_hair(16384, 4, {0.1f, 0.1f}, {0.01f, 0.0001f}, {},
                    shp->quads, shp->pos, shp->norm, shp->texcoord);
            shp->quads.clear();
        } break;
        case test_shape_type::linesi: {
            auto level = 5;
            auto usteps = pow2(level + 2), vsteps = pow2(level + 1);
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere(usteps, vsteps);
        } break;
        case test_shape_type::plight: {
            shp->points.push_back(0);
            shp->pos.push_back({0, 0, 0});
            shp->norm.push_back({0, 0, 1});
            shp->radius.push_back(0.001f);
        } break;
        case test_shape_type::alight: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(1, 1);
            for (auto& p : shp->points) p *= 2;
        } break;
        case test_shape_type::alightr: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(1, 1);
            for (auto& p : shp->points) p *= 2;
        } break;
        case test_shape_type::alightl: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad(1, 1);
            for (auto& p : shp->points) p *= 2;
        } break;
        default: throw runtime_error("bad value");
    }
    // for (auto& p : shp->pos) p *= scale;
    return shp;
}

enum struct test_environment_type {
    none,
    sky1,
    sky2,
};

inline const vector<pair<string, test_environment_type>>&
test_environment_names() {
    static auto names = vector<pair<string, test_environment_type>>{
        {"none", test_environment_type::none},
        {"sky1", test_environment_type::sky1},
        {"sky2", test_environment_type::sky2},
    };
    return names;
}

inline environment* add_test_environment(
    scene* scn, test_environment_type type, const frame3f& frame) {
    if (type == test_environment_type::none) return nullptr;
    auto name = ""s;
    for (auto kv : test_environment_names())
        if (kv.second == type) name += kv.first;
    for (auto env : scn->environments)
        if (env->name == name) return env;
    auto env = new environment();
    env->name = name;
    scn->environments += env;

    switch (type) {
        case test_environment_type::sky1: {
            env->ke = {1, 1, 1};
            env->ke_txt.txt = add_test_texture(scn, test_texture_type::sky1);
        } break;
        case test_environment_type::sky2: {
            env->ke = {1, 1, 1};
            env->ke_txt.txt = add_test_texture(scn, test_texture_type::sky2);
        } break;
        default: throw runtime_error("bad type");
    }

    return env;
}

instance* add_test_instance(scene* scn, test_shape_type stype,
    test_material_type mtype, const frame3f& frame) {
    auto ist = new instance();
    ist->shp = add_test_shape(scn, stype, mtype);
    auto count = 0;
    for (auto ist2 : scn->instances)
        if (ist->shp == ist2->shp) count++;
    ist->name = ist->shp->name + ((count) ? "!" + to_string(count + 1) : ""s);
    ist->frame = frame;
    scn->instances += ist;
    return ist;
}

instance* add_test_instance(scene* scn, test_shape_type stype,
    test_material_type mtype, const vec3f& pos, const vec3f& rot = {0, 0, 0}) {
    return add_test_instance(scn, stype, mtype,
        {rotation_mat3f(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                rotation_mat3f(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                rotation_mat3f(vec3f{1, 0, 0}, rot[0] * pif / 180),
            pos});
}

enum struct test_light_type {
    none,
    pointlight,
    arealight,
    arealight1,
    envlight
};

inline const vector<pair<string, test_light_type>>& test_light_names() {
    static auto names =
        vector<pair<string, test_light_type>>{{"none", test_light_type::none},
            {"pointlight", test_light_type::pointlight},
            {"arealight", test_light_type::arealight},
            {"arealight1", test_light_type::arealight1},
            {"envlight", test_light_type::envlight}};
    return names;
}

tuple<vector<instance*>, environment*> add_test_lights(
    scene* scn, test_light_type type) {
    if (type == test_light_type::none) return {{}, nullptr};
    switch (type) {
        case test_light_type::pointlight: {
            return {{add_test_instance(scn, test_shape_type::plight,
                         test_material_type::light_point, {-2, 10, 8}),
                        add_test_instance(scn, test_shape_type::plight,
                            test_material_type::light_point, {+2, 10, 8})},
                nullptr};
        } break;
        case test_light_type::arealight: {
            return {
                {add_test_instance(scn, test_shape_type::alight,
                     test_material_type::light_area,
                     lookat_frame3f({-4, 5, 8}, {0, 3, 0}, {0, 1, 0}, true)),
                    add_test_instance(scn, test_shape_type::alight,
                        test_material_type::light_area,
                        lookat_frame3f(
                            {+4, 5, 8}, {0, 3, 0}, {0, 1, 0}, true))},
                nullptr};
        } break;
        case test_light_type::arealight1: {
            return {
                {add_test_instance(scn, test_shape_type::alightr,
                     test_material_type::light_arear,
                     lookat_frame3f({+4, 10, 8}, {0, 3, 0}, {0, 1, 0}, true)),
                    add_test_instance(scn, test_shape_type::alightl,
                        test_material_type::light_areal,
                        lookat_frame3f(
                            {-8, 5, 0}, {0, 3, 0}, {0, 1, 0}, true))},
                nullptr};
        } break;
        case test_light_type::envlight: {
            return {
                {}, add_test_environment(scn, test_environment_type::sky1,
                        lookat_frame3f({0, 1, 0}, {0, 1, 1}, {0, 1, 0}, true))};
        } break;
        default: throw runtime_error("bad value");
    }
}

enum struct test_camera_type { none, cam1, cam2, cam3 };

inline const vector<pair<string, test_camera_type>>& test_camera_names() {
    static auto names = vector<pair<string, test_camera_type>>{
        {"none", test_camera_type::none}, {"cam1", test_camera_type::cam1},
        {"cam2", test_camera_type::cam2}, {"cam3", test_camera_type::cam3}};
    return names;
}

camera* add_test_camera(scene* scn, test_camera_type type) {
    if (type == test_camera_type::none) return nullptr;
    auto name = ""s;
    for (auto kv : test_camera_names())
        if (kv.second == type) name = kv.first;
    for (auto cam : scn->cameras)
        if (cam->name == name) return cam;
    auto cam = new camera();
    cam->name = name;
    scn->cameras += cam;
    switch (type) {
        case test_camera_type::cam1: {
            cam->frame = lookat_frame3f({0, 4, 10}, {0, 1, 0}, {0, 1, 0});
            cam->focus = length(vec3f{0, 4, 10} - vec3f{0, 1, 0});
            cam->yfov = 15 * pif / 180;
            cam->aperture = 0;
            cam->aspect = 1;
        } break;
        case test_camera_type::cam2: {
            cam->frame = lookat_frame3f({0, 4, 10}, {0, 1, 0}, {0, 1, 0});
            cam->focus = length(vec3f{0, 4, 10} - vec3f{0, 1, 0});
            cam->yfov = 15 * pif / 180;
            cam->aperture = 0;
            cam->aspect = 16.0f / 9.0f;
        } break;
        case test_camera_type::cam3: {
            cam->frame = lookat_frame3f({0, 4.8f, 12}, {0, 1, 0}, {0, 1, 0});
            cam->focus = length(vec3f{0, 4.8f, 12} - vec3f{0, 1, 0});
            cam->yfov = 15 * pif / 180;
            cam->aperture = 0;
            cam->aspect = 2.35f / 1.0f;  // widescreen cinema
        } break;
        default: throw runtime_error("bad value");
    }
    return cam;
}

scene* make_simple_test_scene(test_camera_type ctype,
    const vector<pair<test_shape_type, test_material_type>>& otypes,
    test_light_type ltype,
    const vector<pair<test_shape_type, test_material_type>>& itypes = {},
    test_material_type fmat = test_material_type::matte_grid) {
    auto scn = new scene();
    add_test_camera(scn, ctype);
    add_test_lights(scn, ltype);

    if (fmat != test_material_type::none) {
        add_test_instance(scn, test_shape_type::floor, fmat, identity_frame3f);
    }

    auto pos1 = vector<float>{0};
    auto pos2 = vector<float>{-1.25f, +1.25f};
    auto pos3 = vector<float>{-2.50f, 0, +2.50f};
    auto pos = pos1;
    if (otypes.size() == 2) pos = pos2;
    if (otypes.size() == 3) pos = pos3;

    for (auto i : range(otypes.size())) {
        add_test_instance(
            scn, otypes[i].first, otypes[i].second, {pos[i], 1, 0});
        if (itypes.empty()) continue;
        add_test_instance(
            scn, itypes[i].first, itypes[i].second, {pos[i], 1, 0});
    }

    return scn;
}

scene* make_instance_scene(
    const vec2i& num, const bbox2f& bbox, uint32_t seed = 13) {
    auto scn = new scene();
    add_test_camera(scn, test_camera_type::cam3);
    add_test_instance(
        scn, test_shape_type::floor, test_material_type::matte_gray, {0, 0, 0});
    add_test_lights(scn, test_light_type::pointlight);

    auto mtypes = vector<test_material_type>{test_material_type::plastic_red,
        test_material_type::plastic_green, test_material_type::plastic_blue};
    auto stypes = vector<test_shape_type>{test_shape_type::sphere,
        test_shape_type::flipcapsphere, test_shape_type::cube};

    auto rscale = 0.9f * 0.25f *
                  min((bbox.max.x - bbox.min.x) / num.x,
                      (bbox.max.x - bbox.min.x) / num.y);
    auto rng = init_rng(seed, 7);
    auto shps = unordered_set<shape*>();
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
            auto ist =
                add_test_instance(scn, stypes[next_rand1i(rng, stypes.size())],
                    mtypes[next_rand1i(rng, mtypes.size())],
                    {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, pos});
            shps.insert(ist->shp);
        }
    }

    for (auto shp : shps) {
        for (auto& p : shp->pos) p *= rscale;
    }

    return scn;
}

scene* make_test_scene(test_scene_type otype) {
    switch (otype) {
        case test_scene_type::cornell_box: {
            return make_cornell_box_scene();
        } break;
        case test_scene_type::textures: {
            auto scn = new scene();
            for (auto ttype : test_texture_names())
                add_test_texture(scn, ttype.second);
            return scn;
        } break;
        case test_scene_type::shapes: {
            auto scn = new scene();
            for (auto stype : {test_shape_type::quad, test_shape_type::cube,
                     test_shape_type::sphere, test_shape_type::spherecube,
                     test_shape_type::spherizedcube,
                     test_shape_type::flipcapsphere, test_shape_type::geosphere,
                     test_shape_type::geospheref, test_shape_type::cubes,
                     test_shape_type::suzanne, test_shape_type::suzannes,
                     test_shape_type::cubefv, test_shape_type::cubefvs,
                     test_shape_type::quads, test_shape_type::matball1,
                     test_shape_type::matball2, test_shape_type::matballi})
                add_test_shape(scn, stype, test_material_type::none);
            return scn;
        } break;
        case test_scene_type::instances_pl: {
            return make_instance_scene({10, 10}, {{-3, -3}, {3, 3}});
        } break;
        case test_scene_type::instancel_pl: {
            return make_instance_scene({100, 100}, {{-3, -3}, {3, 3}});
        } break;
        case test_scene_type::basic_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_red},
                    {test_shape_type::spherecube,
                        test_material_type::plastic_green},
                    {test_shape_type::spherizedcube,
                        test_material_type::plastic_blue},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::simple_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherecube,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherizedcube,
                        test_material_type::plastic_colored},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::simple_al: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherecube,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherizedcube,
                        test_material_type::plastic_colored},
                },
                test_light_type::arealight);
        } break;
        case test_scene_type::simple_el: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherecube,
                        test_material_type::plastic_colored},
                    {test_shape_type::spherizedcube,
                        test_material_type::plastic_colored},
                },
                test_light_type::envlight);
        } break;
        case test_scene_type::transparent_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::quad,
                        test_material_type::transparent_red},
                    {test_shape_type::quad,
                        test_material_type::transparent_green},
                    {test_shape_type::quad,
                        test_material_type::transparent_blue},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::points_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::points, test_material_type::matte_gray},
                    {test_shape_type::points, test_material_type::matte_gray},
                    {test_shape_type::points, test_material_type::matte_gray},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::lines_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::lines1, test_material_type::matte_gray},
                    {test_shape_type::lines2, test_material_type::matte_gray},
                    {test_shape_type::lines2, test_material_type::matte_gray},
                },
                test_light_type::pointlight,
                {
                    {test_shape_type::linesi, test_material_type::matte_gray},
                    {test_shape_type::linesi, test_material_type::matte_gray},
                    {test_shape_type::linesi, test_material_type::matte_gray},
                });
        } break;
        case test_scene_type::subdiv_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::cubes, test_material_type::plastic_red},
                    {test_shape_type::suzannes,
                        test_material_type::plastic_green},
                    {test_shape_type::suzannes,
                        test_material_type::plastic_blue},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::matball1_al: {
            return make_simple_test_scene(test_camera_type::cam1,
                {
                    {test_shape_type::matball1,
                        test_material_type::plastic_red},
                },
                test_light_type::arealight1,
                {{test_shape_type::matballi, test_material_type::matte_gray}});
        } break;
        case test_scene_type::matball1_el: {
            return make_simple_test_scene(test_camera_type::cam1,
                {
                    {test_shape_type::matball1,
                        test_material_type::plastic_red},
                },
                test_light_type::envlight,
                {{test_shape_type::matballi, test_material_type::matte_gray}});
        } break;
        case test_scene_type::tesselation_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::geospherel,
                        test_material_type::matte_gray},
                    {test_shape_type::geospheref,
                        test_material_type::matte_gray},
                    {test_shape_type::geosphere,
                        test_material_type::matte_gray},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::textureuv_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::matte_green},
                    {test_shape_type::flipcapsphere,
                        test_material_type::matte_colored},
                    {test_shape_type::flipcapsphere,
                        test_material_type::matte_uv},
                },
                test_light_type::pointlight);
        } break;
        case test_scene_type::normalmap_pl: {
            return make_simple_test_scene(test_camera_type::cam3,
                {
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_blue},
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_bumped},
                    {test_shape_type::flipcapsphere,
                        test_material_type::plastic_bumped},
                },
                test_light_type::pointlight);
        } break;
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
        default: throw runtime_error("bad value");
    }
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
    if (enabled)
        glEnable(GL_DEPTH_TEST);
    else
        glDisable(GL_DEPTH_TEST);
    assert(gl_check_error());
}

// Enable/disable culling
void gl_enable_culling(bool enabled) {
    assert(gl_check_error());
    if (enabled)
        glEnable(GL_CULL_FACE);
    else
        glDisable(GL_CULL_FACE);
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

// Implementation of make_texture.
void _init_texture(gl_texture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb) {
    txt._width = w;
    txt._height = h;
    txt._ncomp = nc;
    txt._float = as_float;
    txt._srgb = as_srgb;
    txt._mipmap = mipmap;
    assert(!as_srgb || !as_float);
    assert(gl_check_error());
    int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
    int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
    int* formats =
        (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
    assert(gl_check_error());
    glGenTextures(1, &txt._tid);
    glBindTexture(GL_TEXTURE_2D, txt._tid);
    glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0, formats_ub[nc - 1],
        (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE, pixels);
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
    assert(gl_check_error());
}

// Implementation of update_texture.
void _update_texture(
    gl_texture& txt, int w, int h, int nc, const void* pixels, bool floats) {
    txt._width = w;
    txt._height = h;
    assert(gl_check_error());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glBindTexture(GL_TEXTURE_2D, txt._tid);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, formats[nc - 1],
        (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE, pixels);
    if (txt._mipmap) glGenerateMipmap(GL_TEXTURE_2D);
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

// Creates a buffer with num elements of size size stored in values, where
// content is dyanamic if dynamic.
void _init_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic) {
    buf._num = n;
    buf._ncomp = nc;
    buf._float = as_float;
    assert(gl_check_error());
    buf._bid = (GLuint)0;
    glGenBuffers(1, &buf._bid);
    glBindBuffer(GL_ARRAY_BUFFER, buf._bid);
    glBufferData(GL_ARRAY_BUFFER,
        buf._num * buf._ncomp * ((as_float) ? sizeof(float) : sizeof(int)),
        values, (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    assert(gl_check_error());
}

// Updates the buffer bid with new data.
void _update_vertex_buffer(
    gl_vertex_buffer& buf, int n, int nc, const void* values, bool as_float) {
    buf._num = n;
    buf._ncomp = nc;
    buf._float = as_float;
    assert(gl_check_error());
    glBindBuffer(GL_ARRAY_BUFFER, buf._bid);
    glBufferSubData(GL_ARRAY_BUFFER, 0,
        buf._num * buf._ncomp * ((as_float) ? sizeof(float) : sizeof(int)),
        values);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
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

// Creates a buffer with num elements of size size stored in values, where
// content is dyanamic if dynamic.
// Returns the buffer id.
void _init_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values, bool dynamic) {
    buf._num = n;
    buf._ncomp = nc;
    assert(gl_check_error());
    buf._bid = (GLuint)0;
    glGenBuffers(1, &buf._bid);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf._bid);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, buf._num * buf._ncomp * sizeof(int),
        values, (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(gl_check_error());
}

// Updates the buffer bid with new data.
void _update_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values) {
    buf._num = n;
    buf._ncomp = nc;
    assert(gl_check_error());
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf._bid);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
        buf._num * buf._ncomp * sizeof(int), values);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
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
            int mtype;         // material type
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
            if(material.mtype == 0) {
                type = 0;
            } else if(material.mtype == 1) {
                type = material.etype;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt;
                rs = rs*rs;
            } else if(material.mtype == 2) {
                type = material.etype;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(material.mtype == 3) {
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
    if (!glewInit()) return nullptr;
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

// Initialize widgets
void init_widgets(gl_window* win) {
    ImGui_ImplGlfwGL3_Init(win->_gwin, false);
    ImGui::GetStyle().WindowRounding = 0;
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::SetNextWindowPos({0, 0});
    auto size = get_window_size(win);
    ImGui::SetNextWindowSize({(float)win->_widget_width, (float)size[1]});
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
bool draw_color_widget(gl_window* win, const string& lbl, vec4b& val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)&val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)&val = valb;
        return true;
    }
    return false;
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

// List widget
bool draw_list_widget(gl_window* win, const string& lbl, int& val,
    const vector<pair<string, int>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::ListBox(lbl.c_str(), &cur, _enum_widget_labels_int,
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

namespace __impl_scn_widgets {

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
        if (ist->shp) draw_tree_widgets(win, "shape: ", ist->shp, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, scene* scn, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + "cameras")) {
        for (auto cam : scn->cameras)
            draw_tree_widgets(win, "", cam, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "shapes")) {
        for (auto msh : scn->shapes) draw_tree_widgets(win, "", msh, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "instances")) {
        for (auto ist : scn->instances)
            draw_tree_widgets(win, "", ist, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "materials")) {
        for (auto mat : scn->materials)
            draw_tree_widgets(win, "", mat, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "textures")) {
        for (auto txt : scn->textures)
            draw_tree_widgets(win, "", txt, selection);
        draw_tree_widget_end(win);
    }
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
    static auto mtype_names = vector<pair<string, material_type>>{
        {"generic", material_type::specular_roughness},
        {"metallic_roughness", material_type::metallic_roughness},
        {"specular_glossiness", material_type::specular_glossiness},
    };

    auto txt_names = vector<pair<string, texture*>>{{"<none>", nullptr}};
    for (auto txt : scn->textures) txt_names.push_back({txt->path, txt});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", mat->name);
    edited += draw_value_widget(win, "mtype", mat->mtype, mtype_names);
    edited += draw_value_widget(win, "ke", mat->ke, 0, 1000);
    edited += draw_value_widget(win, "kd", mat->kd, 0, 1);
    edited += draw_value_widget(win, "ks", mat->ks, 0, 1);
    edited += draw_value_widget(win, "kt", mat->kt, 0, 1);
    edited += draw_value_widget(win, "rs", mat->rs, 0, 1);

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
    draw_label_widget(win, "verts", (int)shp->pos.size());
    if (!shp->triangles.empty())
        draw_label_widget(win, "triangles", (int)shp->triangles.size());
    if (!shp->lines.empty())
        draw_label_widget(win, "lines", (int)shp->lines.size());
    if (!shp->points.empty())
        draw_label_widget(win, "points", (int)shp->points.size());
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
    edited += draw_value_widget(win, "aperture", cam->aperture, 0, 1);
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

inline bool draw_elem_widgets(gl_window* win, scene* scn, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    for (auto cam : scn->cameras) {
        if (cam == selection)
            return draw_elem_widgets(win, scn, cam, selection, gl_txt);
    }

    for (auto shp : scn->shapes) {
        if (shp == selection)
            return draw_elem_widgets(win, scn, shp, selection, gl_txt);
    }

    for (auto ist : scn->instances) {
        if (ist == selection)
            return draw_elem_widgets(win, scn, ist, selection, gl_txt);
    }

    for (auto mat : scn->materials) {
        if (mat == selection)
            return draw_elem_widgets(win, scn, mat, selection, gl_txt);
    }

    for (auto txt : scn->textures) {
        if (txt == selection)
            return draw_elem_widgets(win, scn, txt, selection, gl_txt);
    }

    return false;
}

inline bool draw_scene_widgets(gl_window* win, const string& lbl, scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    if (draw_header_widget(win, lbl)) {
        // draw_scroll_widget_begin(win, "model", 240, false);
        draw_tree_widgets(win, "", scn, selection);
        // draw_scroll_widget_end(win);
        return draw_elem_widgets(win, scn, selection, gl_txt);
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

}  // namespace __impl_scn_widgets

bool draw_scene_widgets(gl_window* win, const string& lbl, scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    return __impl_scn_widgets::draw_scene_widgets(
        win, lbl, scn, selection, gl_txt);
}

}  // namespace ygl

#endif

// HACK to avoid compilation with MSVC2015 without dirtying code
#ifdef constexpr
#undef constexpr
#endif
