//
// # Yocto/Image: Tiny imaging Library mostly for rendering and color support
//
//
// Yocto/Image is a collection of image utilities useful when writing rendering
// algorithms. These include a simple image data structure, color conversion
// utilities and tone mapping.
//
//
// ## Image Utilities
//
// Yocto/Image supports a very small set is color and image utilities including
// color utilities, example image creation, tone mapping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// global illumination renderer, rather than the need of generic image editing.
// We support 4-channels float images (assumed to be in linear color) and 
// 4-channels byte images (assumed to be in sRGB).
// 
//
// 1. create images with `make_image4f()`/`make_image4b()`
// 2. resize images with `resize_image()`
// 3. tonemap images with `tonemap_image()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky_image()`
// 0. load and save image with Yocto/ImageIO
//
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

#ifndef _YOCTO_IMAGE_H_
#define _YOCTO_IMAGE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container.
struct image4f {
    int           width  = 0;
    int           height = 0;
    vector<vec4f> pixels = {};
};

// Image container.
struct image4b {
    int           width  = 0;
    int           height = 0;
    vector<vec4b> pixels = {};
};

// Image creation
inline image4f make_image(int width, int height, const vec4f& value) {
    return {width, height, {(size_t)(width * height), value}};
}
inline image4b make_image(int width, int height, const vec4b& value) {
    return {width, height, {(size_t)(width * height), value}};
}
inline image4f make_image(int width, int height, const vec4f* values) {
    return {width, height, {values, values + (size_t)(width * height)}};
}
inline image4b make_image(int width, int height, const vec4b* values) {
    return {width, height, {values, values + (size_t)(width * height)}};
}

// Pixel access
inline vec4f& at(image4f& img, const vec2i& ij) {
    return img.pixels[ij.y * img.width + ij.x];
}
inline const vec4f& at(const image4f& img, const vec2i& ij) {
    return img.pixels[ij.y * img.width + ij.x];
}
inline vec4b& at(image4b& img, const vec2i& ij) {
    return img.pixels[ij.y * img.width + ij.x];
}
inline const vec4b& at(const image4b& img, const vec2i& ij) {
    return img.pixels[ij.y * img.width + ij.x];
}
inline vec4f& at(image4f& img, int i, int j) {
    return img.pixels[j * img.width + i];
}
inline const vec4f& at(const image4f& img, int i, int j) {
    return img.pixels[j * img.width + i];
}
inline vec4b& at(image4b& img, int i, int j) {
    return img.pixels[j * img.width + i];
}
inline const vec4b& at(const image4b& img, int i, int j) {
    return img.pixels[j * img.width + i];
}

// Functions to query image data
inline vec2i  imsize(const image4f& img) { return {img.width, img.height}; }
inline vec2i  imsize(const image4b& img) { return {img.width, img.height}; }
inline bool   empty(const image4f& img) { return empty(img.pixels); }
inline bool   empty(const image4b& img) { return empty(img.pixels); }
inline size_t size(const image4f& img) { return img.pixels.size(); }
inline size_t size(const image4b& img) { return img.pixels.size(); }
inline vec4f* begin(image4f& img) { return data(img.pixels); }
inline vec4f* end(image4f& img) { return data(img.pixels) + img.pixels.size(); }
inline const vec4f* begin(const image4f& img) { return data(img.pixels); }
inline const vec4f* end(const image4f& img) {
    return data(img.pixels) + img.pixels.size();
}
inline vec4b* begin(image4b& img) { return data(img.pixels); }
inline vec4b* end(image4b& img) { return data(img.pixels) + img.pixels.size(); }
inline const vec4b* begin(const image4b& img) { return data(img.pixels); }
inline const vec4b* end(const image4b& img) {
    return data(img.pixels) + img.pixels.size();
}
inline vec4f*       data(image4f& img) { return data(img.pixels); }
inline vec4b*       data(image4b& img) { return data(img.pixels); }
inline const vec4f* data(const image4f& img) { return data(img.pixels); }
inline const vec4b* data(const image4b& img) { return data(img.pixels); }

// Image region
struct image_region {
    int offsetx = 0;
    int offsety = 0;
    int width   = 0;
    int height  = 0;
};

// Splits an image into an array of regions
vector<image_region> make_image_regions(
    int width, int height, int region_size = 32);

// Gets pixels in an image region
image4f get_image_region(const image4f& img, const image_region& region);
image4b get_image_region(const image4b& img, const image_region& region);

// Conversion from/to floats.
image4f byte_to_float(const image4b& bt);
image4b float_to_byte(const image4f& fl);

// Conversion between linear and gamma-encoded images.
image4f srgb_to_linear(const image4f& srgb);
image4f linear_to_srgb(const image4f& lin);

// Conversion between linear and gamma-encoded images.
image4f gamma_to_linear(const image4f& srgb, float gamma);
image4f linear_to_gamma(const image4f& lin, float gamma);

// Apply exposure and filmic tone mapping
image4f tonemap_image(const image4f& hdr, float exposure, bool filmic, bool srgb);
void    tonemap_image_region(image4f& ldr, const image_region& region,
       const image4f& hdr, float exposure, bool filmic, bool srgb);

// Resize an image.
image4f resize_image(const image4f& img, int width, int height);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make example images.
image4f make_grid_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_checker_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_bumpdimple_image(int width, int height, int tile = 8);
image4f make_ramp_image(int width, int height, const vec4f& c0, const vec4f& c1,
    float srgb = false);
image4f make_gammaramp_image(int width, int height);
image4f make_uvramp_image(int width, int height);
image4f make_uvgrid_image(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale = 1);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun is not simple to control
// procedurally, so we include two variables to unrealistically customize it:
// sun_size_scale and sun_emission_scale with good values in [20-] for both.
// To make editing easier, we also renormalize the (sun+sky) integral to
// the sky only integral if renormalize_sun is true. For the same reason,
// we rescale the sun dimension such that it covers at least 5 pixels in
// diameter if has_sun is enabled.
image4f make_sunsky_image(int width, int height, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_angle_scale = 1.0f,
    float        sun_emission_scale = 1.0f,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f}, bool renormalize_sun = true);
// Make an image of multiple lights.
image4f make_lights_image(int width, int height, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container.
struct volume1f {
    int           width  = 0;
    int           height = 0;
    int           depth  = 0;
    vector<float> voxels = {};
};

// Image creation
inline volume1f make_volume(int width, int height, int depth, float value) {
    return {width, height, depth,
        vector<float>((size_t)(width * height * depth), value)};
}
inline volume1f make_volume(
    int width, int height, int depth, const float* values) {
    return {width, height, depth,
        vector<float>(values, values + (size_t)(width * height * depth))};
}

// Pixel access
inline float& at(volume1f& vol, const vec3i& ijk) {
    return vol.voxels[ijk.z * vol.width * vol.height + ijk.y * vol.width + ijk.x];
}
inline const float& at(const volume1f& vol, const vec3i& ijk) {
    return vol.voxels[ijk.z * vol.width * vol.height + ijk.y * vol.width + ijk.x];
}
inline float& at(volume1f& vol, int i, int j, int k) {
    return vol.voxels[k * vol.width * vol.height + j * vol.width + i];
}
inline const float& at(const volume1f& vol, int i, int j, int k) {
    return vol.voxels[k * vol.width * vol.height + j * vol.width + i];
}

// Functions to query volume data
inline vec3i volsize(const volume1f& vol) {
    return {vol.width, vol.height, vol.depth};
}
inline bool   empty(const volume1f& vol) { return empty(vol.voxels); }
inline size_t size(const volume1f& vol) { return vol.voxels.size(); }
inline float* begin(volume1f& vol) { return data(vol.voxels); }
inline float* end(volume1f& vol) {
    return data(vol.voxels) + vol.voxels.size();
}
inline const float* begin(const volume1f& vol) { return data(vol.voxels); }
inline const float* end(const volume1f& vol) {
    return data(vol.voxels) + vol.voxels.size();
}
inline float*       data(volume1f& vol) { return data(vol.voxels); }
inline const float* data(const volume1f& vol) { return data(vol.voxels); }

// make a simple example volume
volume1f make_test_volume(
    int width, int height, int depth, float scale = 10, float exponent = 6);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {
        (byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255),
        (byte)clamp(int(a.z * 256), 0, 255),
        (byte)clamp(int(a.w * 256), 0, 255),
    };
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma), srgb.w};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma),
        lin.w};
}

// sRGB non-linear curve
inline float srgb_to_linear(float srgb) {
    if (srgb <= 0.04045) {
        return srgb / 12.92f;
    } else {
        return pow((srgb + 0.055f) / (1.0f + 0.055f), 2.4f);
    }
}
inline float linear_to_srgb(float lin) {
    if (lin <= 0.0031308f) {
        return 12.92f * lin;
    } else {
        return (1 + 0.055f) * pow(lin, 1 / 2.4f) - 0.055f;
    }
}

// Conversion between linear and srgb colors.
inline vec3f srgb_to_linear(const vec3f& srgb) {
    return {
        srgb_to_linear(srgb.x), srgb_to_linear(srgb.y), srgb_to_linear(srgb.z)};
}
inline vec3f linear_to_srgb(const vec3f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z)};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {srgb_to_linear(srgb.x), srgb_to_linear(srgb.y),
        srgb_to_linear(srgb.z), srgb.w};
}
inline vec4f linear_to_srgb(const vec4f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z),
        lin.w};
}

// Approximate luminance estimate for sRGB primaries (better relative luminance)
inline float luminance(const vec3f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a.z);
}
inline float luminance(const vec4f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a.z);
}

// Fitted ACES tonemapping curve.
inline float tonemap_filmic(float hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    // hdr *= 0.6; // brings it back to ACES range
    return (hdr * hdr * 2.51f + hdr * 0.03f) /
           (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
}
// Apply ACES fitted curve.
inline vec4f tonemap_filmic(const vec4f& hdr) {
    return {tonemap_filmic(hdr.x), tonemap_filmic(hdr.y), tonemap_filmic(hdr.z),
        hdr.w};
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_filmic(
    const vec4f& hdr, float exposure, bool filmic, bool srgb) {
    auto scale = pow(2.0f, exposure);
    auto ldr   = vec4f{hdr.x * scale, hdr.y * scale, hdr.z * scale, hdr.w};
    if (filmic) ldr = tonemap_filmic(ldr);
    if (srgb) ldr = linear_to_srgb(ldr);
    return ldr;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    return {
        +3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z,
    };
}
// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    return {
        0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z,
    };
}

// Converts HSV to RGB.
inline vec3f hsv_to_rgb(const vec3f& hsv);
inline vec3f rgb_to_hsv(const vec3f& rgb);

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Convert HSV to RGB
inline vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0.0f) return {v, v, v};

    h       = fmodf(h, 1.0f) / (60.0f / 360.0f);
    int   i = (int)h;
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
inline vec3f rgb_to_hsv(const vec3f& rgb) {
    // from Imgui.cpp
    auto  r = rgb.x, g = rgb.y, b = rgb.z;
    float K = 0.f;
    if (g < b) {
        swap(g, b);
        K = -1.f;
    }
    if (r < g) {
        swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - (g < b ? g : b);
    return {
        fabsf(K + (g - b) / (6.f * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

}  // namespace yocto

#endif
