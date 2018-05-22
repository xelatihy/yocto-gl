//
// # Yocto/Image: Tiny C++ Library of Image Utilities
//
// Yocto/Image is a collection of utilities for image loading and saving,
// example image creation, tonempping, image resizing, and sunsky procedural
// images. The library contains only color conversion utilities.
// Yocto/Image is written to support the need of a minimal, but
// fully-featured, global illumination renderer, rather than the need of
// generic image editing. Image IO and resizing are wrappers of `stb_image`
// and `tinyexr` libraries.
//
// Yocto/Img support only two image types, `image4b` for low-dynamic range
// images and `image4f` for high-dynamic range images. Pixels are represented
// as RGBA values of either 8-bytes ints (LDRs) or float (HDRs) types.
//
// ## Usage
//
// 1. create images with `make_image4b()` and `make_image4f()`
// 1. load images with `load_image4b()` or `load_image4f()`
// 2. save images with `save_image4b()` or `save_image4f()`
// 3. resize images with `resize_image()`
// 4. tonemap images with `tonemap_image()`
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky_image()`
// 7. we also support low level image IO with `load_XXX()` and `save_XXX()`
//    functions.
// 8. see list below for the color conversion utilities supported
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

#ifndef _YGL_IMAGE_H_
#define _YGL_IMAGE_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

// enable image io
#ifndef YGL_IMAGEIO
#define YGL_IMAGEIO 1
#endif

#include "yocto_math.h"

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)max(0, min(int(a.x * 256), 255)),
        (byte)max(0, min(int(a.y * 256), 255)),
        (byte)max(0, min(int(a.z * 256), 255)),
        (byte)max(0, min(int(a.w * 256), 255))};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Approximate conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {pow(srgb.x / 255.0f, 2.2f), pow(srgb.y / 255.0f, 2.2f),
        pow(srgb.z / 255.0f, 2.2f), srgb.w / 255.0f};
}
// Approximate conversion to srgb.
inline vec4b linear_to_srgb(const vec4f& lin) {
    return float_to_byte({pow(lin.x, 1 / 2.2f), pow(lin.y, 1 / 2.2f),
        pow(lin.z, 1 / 2.2f), lin.w});
}

// Approximate luminance estimate
inline float luminance(const vec3f& a) { return (a.x + a.y + a.z) / 3; }

// Converts HSV to RGB.
vec3f hsv_to_rgb(const vec3f& hsv);
vec3f rgb_to_hsv(const vec3f& rgb);
// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz);
vec3f xyY_to_xyz(const vec3f& xyY);
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz);
vec3f rgb_to_xyz(const vec3f& rgb);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Approximate conversion from/to srgb.
inline std::vector<vec4f> srgb_to_linear(const std::vector<vec4b>& srgb) {
    auto lin = std::vector<vec4f>(srgb.size());
    for (auto i = 0; i < srgb.size(); i++) lin[i] = srgb_to_linear(srgb[i]);
    return lin;
}
inline std::vector<vec4b> linear_to_srgb(const std::vector<vec4f>& lin) {
    auto srgb = std::vector<vec4b>(lin.size());
    for (auto i = 0; i < lin.size(); i++) srgb[i] = linear_to_srgb(lin[i]);
    return srgb;
}

// Conversion from/to floats.
inline std::vector<vec4f> byte_to_float(const std::vector<vec4b>& bt) {
    auto fl = std::vector<vec4f>(bt.size());
    for (auto i = 0; i < bt.size(); i++) fl[i] = byte_to_float(bt[i]);
    return fl;
}
inline std::vector<vec4b> float_to_byte(const std::vector<vec4f>& fl) {
    auto bt = std::vector<vec4b>(fl.size());
    for (auto i = 0; i < fl.size(); i++) bt[i] = float_to_byte(fl[i]);
    return bt;
}

// Apply exposure to an image
std::vector<vec4f> expose_image(const std::vector<vec4f>& hdr, float exposure);

// Tone mapping linear HDR to linear LDR images in reference space.
std::vector<vec4f> filmic_tonemap_image(const std::vector<vec4f>& hdr);
std::vector<vec4f> aces_tonemap_image(const std::vector<vec4f>& hdr);

// Make example images.
std::vector<vec4f> make_grid_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1.0f},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1.0f});
std::vector<vec4f> make_checker_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1.0f},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1.0f});
std::vector<vec4f> make_bumpdimple_image(int width, int height, int tile = 8);
std::vector<vec4f> make_ramp_image(
    int width, int height, const vec4b& c0, const vec4b& c1, bool srgb = false);
std::vector<vec4f> make_gammaramp_image(int width, int height);
std::vector<vec4f> make_uv_image(int width, int height);
std::vector<vec4f> make_uvgrid_image(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
std::vector<vec4f> bump_to_normal_map(
    int width, int height, const std::vector<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
// in [1.7,10] with or without sun.
std::vector<vec4f> make_sunsky_image(int width, int height, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
std::vector<vec4f> make_lights_image(int width, int height,
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pi / 4,
    float lwidth = pi / 16, float lheight = pi / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
std::vector<vec4f> make_noise_image(
    int width, int height, float scale = 1, bool wrap = true);
std::vector<vec4f> make_fbm_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
std::vector<vec4f> make_ridge_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
std::vector<vec4f> make_turbulence_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

#if YGL_IMAGEIO

// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

// Loads/saves a 4 channel ldr/hdr image.
std::vector<vec4b> load_image4b(
    const std::string& filename, int& width, int& height);
std::vector<vec4f> load_image4f(const std::string& filename, int& width,
    int& height, bool srgb_8bit = true);
bool save_image4b(const std::string& filename, int width, int height,
    const std::vector<vec4b>& img);
bool save_image4f(const std::string& filename, int width, int height,
    const std::vector<vec4f>& img, bool srgb_8bit = true);

// Filter type and edge mode for resizing.
enum struct resize_filter {
    def,
    box,
    triangle,
    cubic_spline,
    catmull_rom,
    mitchell
};
enum struct resize_edge { def, clamp, reflect, wrap, zero };

// Resize an image.
std::vector<vec4f> resize_image(int width, int height,
    const std::vector<vec4f>& img, int res_width, int res_height,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);
std::vector<vec4b> resize_image(int width, int height,
    const std::vector<vec4b>& img, int res_width, int res_height,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);

#endif

}  // namespace ygl

#endif
