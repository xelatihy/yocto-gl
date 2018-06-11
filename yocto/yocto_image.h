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
// 1. load images with `load_image4b()` or `load_image()`
// 2. save images with `save_image4b()` or `save_image()`
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

#include "yocto_math.h"

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// IMAGE TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Float image.
struct image4f {
    int width = 0;
    int height = 0;
    std::vector<vec4f> pxl;
};

// Byte image.
struct image4b {
    int width = 0;
    int height = 0;
    std::vector<vec4b> pxl;
};

// Create image.
image4f make_image4f(int width, int height, vec4f c = {0, 0, 0, 0});
image4b make_image4b(int width, int height, vec4b c = {0, 0, 0, 0});

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Conversion from/to floats.
image4f byte_to_float(const image4b& bt);
image4b float_to_byte(const image4f& fl);

// Conversion between linear and gamma-encoded images.
image4f gamma_to_linear(const image4f& srgb, float gamma = 2.2f);
image4f linear_to_gamma(const image4f& lin, float gamma = 2.2f);

// Apply exposure and filmic tone mapping
image4f expose_image(const image4f& hdr, float exposure);
image4f filmic_tonemap_image(const image4f& hdr);

// Resize an image.
image4f resize_image(const image4f& img, int res_width, int res_height);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

// Loads/saves a 4 channel image.
image4f load_image(const std::string& filename, float ldr_gamma = 2.2f);
void save_image(
    const std::string& filename, const image4f& img, float ldr_gamma = 2.2f);
image4f load_image_from_memory(
    const byte* data, int data_size, float ldr_gamma = 2.2f);

}  // namespace ygl

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(vec4f a) {
    return {(byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255),
        (byte)clamp(int(a.z * 256), 0, 255),
        (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Conversion between linear and gamma-encoded images.
inline vec3f gamma_to_linear(vec3f srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma)};
}
inline vec3f linear_to_gamma(vec3f lin, float gamma = 2.2f) {
    return {
        pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma)};
}
inline vec4f gamma_to_linear(vec4f srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma), srgb.w};
}
inline vec4f linear_to_gamma(vec4f lin, float gamma = 2.2f) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma),
        lin.w};
}

// Approximate luminance estimate
inline float luminance(vec3f a) { return (a.x + a.y + a.z) / 3; }

// Converts HSV to RGB.
vec3f hsv_to_rgb(vec3f hsv);
vec3f rgb_to_hsv(vec3f rgb);
// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(vec3f xyz);
vec3f xyY_to_xyz(vec3f xyY);
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(vec3f xyz);
vec3f rgb_to_xyz(vec3f rgb);

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

// Make example images.
image4f make_grid_image(int width, int height, int tile = 8,
    vec4f c0 = {0.5f, 0.5f, 0.5f, 1}, vec4f c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_checker_image(int width, int height, int tile = 8,
    vec4f c0 = {0.5f, 0.5f, 0.5f, 1}, vec4f c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_bumpdimple_image(int width, int height, int tile = 8);
image4f make_ramp_image(
    int width, int height, vec4f c0, vec4f c1, float srgb = false);
image4f make_gammaramp_image(int width, int height);
image4f make_uv_image(int width, int height);
image4f make_uvgrid_image(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
// in [1.7,10] with or without sun.
image4f make_sunsky_image(int width, int height, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    vec3f ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
image4f make_lights_image(int width, int height, vec3f le = {1, 1, 1},
    int nlights = 4, float langle = pi / 4, float lwidth = pi / 16,
    float lheight = pi / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4f make_noise_image(
    int width, int height, float scale = 1, bool wrap = true);
image4f make_fbm_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image4f make_ridge_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image4f make_turbulence_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace ygl

#endif
