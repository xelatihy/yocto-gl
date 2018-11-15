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
// Yocto/GL supports a very small set is color and image utilities including
// color utilities, example image creation, tone mapping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// global illumination renderer, rather than the need of generic image editing.
//
// 0. load and save image with Yocto/GLIO
// 1. create images with `image<T>` data structure
// 2. resize images with `resize_image()`
// 3. tonemap images with `tonemap_filmic()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky_image()`
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
    image4f() : size{0, 0}, pixels{} {}
    image4f(const vec2i& size, const vec4f& value = {})
        : size{size}, pixels((size_t)(size.x * size.y), value) {}
    image4f(const vec2i& size, const vec4f* values)
        : size{size}, pixels(values, values + size.x * size.y) {}

    bool  empty() const { return pixels.empty(); }

    void clear() {
        size = {0, 0};
        pixels.clear();
    }
    void resize(const vec2i& size, const vec4f& value = {});

    vec4f& operator[](const vec2i& ij) { return pixels[ij.y * size.x + ij.x]; }
    const vec4f& operator[](const vec2i& ij) const {
        return pixels[ij.y * size.x + ij.x];
    }
    vec4f&       at(const vec2i& ij) { return operator[](ij); }
    const vec4f& at(const vec2i& ij) const { return operator[](ij); }

    vec4f*       data() { return pixels.data(); }
    const vec4f* data() const { return pixels.data(); }

    vec4f*       begin() { return pixels.data(); }
    vec4f*       end() { return pixels.data() + pixels.size(); }
    const vec4f* begin() const { return pixels.data(); }

    const vec4f* end() const { return pixels.data() + pixels.size(); }

    vec2i     size   = {0, 0};
    vector<vec4f> pixels = {};
};

// Image container.
struct image4b {
    image4b() : size{0, 0}, pixels{} {}
    image4b(const vec2i& size, const vec4b& value = {})
        : size{size}, pixels((size_t)(size.x * size.y), value) {}
    image4b(const vec2i& size, const vec4b* values)
        : size{size}, pixels(values, values + size.x * size.y) {}

    bool  empty() const { return pixels.empty(); }

    void clear() {
        size = {0, 0};
        pixels.clear();
    }
    void resize(const vec2i& size, const vec4b& value = {});

    vec4b& operator[](const vec2i& ij) { return pixels[ij.y * size.x + ij.x]; }
    const vec4b& operator[](const vec2i& ij) const {
        return pixels[ij.y * size.x + ij.x];
    }
    vec4b&       at(const vec2i& ij) { return operator[](ij); }
    const vec4b& at(const vec2i& ij) const { return operator[](ij); }

    vec4b*       data() { return pixels.data(); }
    const vec4b* data() const { return pixels.data(); }

    vec4b*       begin() { return pixels.data(); }
    vec4b*       end() { return pixels.data() + pixels.size(); }
    const vec4b* begin() const { return pixels.data(); }

    const vec4b* end() const { return pixels.data() + pixels.size(); }

    vec2i     size   = {0, 0};
    vector<vec4b> pixels = {};
};

// Size
inline float get_image_aspect(const image4f& img);
inline float get_image_aspect(const image4b& img);

// Splits an image into an array of regions
vector<bbox2i> make_image_regions(const vec2i& image_size, int region_size = 32);

// Gets pixels in an image region
inline image4f get_image_region(const image4f& img, const bbox2i& region);
inline image4b get_image_region(const image4b& img, const bbox2i& region);

// Gets an image size from a suggested size and an aspect ratio. The suggested
// size may have zeros in either components. In which case, we use the aspect
// ration to compute the other.
vec2i get_image_size(const vec2i& size, float aspect);

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
image4f tonemap_image(
    const image4f& hdr, float exposure, bool filmic, bool srgb);
void tonemap_image_region(image4f& ldr, const bbox2i& region,
    const image4f& hdr, float exposure, bool filmic, bool srgb);

// Resize an image.
image4f resize_image(const image4f& img, const vec2i& size);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make example images.
image4f make_grid_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_checker_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_bumpdimple_image(const vec2i& size, int tile = 8);
image4f make_ramp_image(
    const vec2i& size, const vec4f& c0, const vec4f& c1, float srgb = false);
image4f make_gammaramp_image(const vec2i& size);
image4f make_uvramp_image(const vec2i& size);
image4f make_uvgrid_image(
    const vec2i& size, int tile = 8, bool colored = true);

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
image4f make_sunsky_image(const vec2i& size, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_angle_scale = 1.0f,
    float        sun_emission_scale = 1.0f,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f}, bool renormalize_sun = true);
// Make an image of multiple lights.
image4f make_lights_image(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4f make_noise_image(
    const vec2i& size, float scale = 1, bool wrap = true);
image4f make_fbm_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image4f make_ridge_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image4f make_turbulence_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container.
struct volume1f {
    volume1f() : size{0, 0, 0}, voxels{} {}
    volume1f(const vec3i& size, float value = 0)
        : size{size}, voxels((size_t)(size.x * size.y * size[2]), value) {}
    volume1f(const vec3i& size, const float* values)
        : size{size}, voxels(values, values + size.x * size.y + size[2]) {}

    bool  empty() const { return voxels.empty(); }

    void clear() {
        size = {0, 0, 0};
        voxels.clear();
    }
    void resize(const vec3i& size, float value = 0);

    float& operator[](const vec3i& ijk) {
        return voxels[ijk[2] * size.x * size.y + ijk.y * size.x + ijk.x];
    }
    const float& operator[](const vec3i& ijk) const {
        return voxels[ijk[2] * size.x * size.y + ijk.y * size.x + ijk.x];
    }

    float*       data() { return voxels.data(); }
    const float* data() const { return voxels.data(); }

    vec3i     size   = {0, 0, 0};
    vector<float> voxels = {};
};

// make a simple example volume
volume1f make_test_volume(
    const vec3i& size, float scale = 10, float exponent = 6);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a);
inline vec4f byte_to_float(const vec4b& a);

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma = 2.2f);
inline vec3f linear_to_gamma(const vec3f& lin, float gamma = 2.2f);
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma = 2.2f);
inline vec4f linear_to_gamma(const vec4f& lin, float gamma = 2.2f);

// sRGB non-linear curve
inline float srgb_to_linear(float srgb);
inline float linear_to_srgb(float lin);

// Conversion between linear and srgb colors.
inline vec3f srgb_to_linear(const vec3f& srgb);
inline vec3f linear_to_srgb(const vec3f& lin);
inline vec4f srgb_to_linear(const vec4f& srgb);
inline vec4f linear_to_srgb(const vec4f& lin);

// Approximate luminance estimate for sRGB primaries (better relative luminance)
inline float luminance(const vec3f& a);
inline float luminance(const vec4f& a);

// Fitted ACES tonemapping curve.
inline float tonemap_filmic(float hdr);
// Apply ACES fitted curve.
inline vec4f tonemap_filmic(const vec4f& hdr);

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_filmic(
    const vec4f& hdr, float exposure, bool filmic, bool srgb);

// Converts HSV to RGB.
inline vec3f hsv_to_rgb(const vec3f& hsv);
inline vec3f rgb_to_hsv(const vec3f& rgb);
// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz);
inline vec3f xyY_to_xyz(const vec3f& xyY);
// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz);
inline vec3f rgb_to_xyz(const vec3f& rgb);

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

inline void image4f::resize(const vec2i& size_, const vec4f& value) {
    if (size == size_) return;
    if (size == zero2i) {
        *this = image4f{size_, value};
    } else if (size_ == zero2i) {
        clear();
    } else {
        auto img = image4f{size_, value};
        for (auto j = 0; j < min(size_.y, size_.y); j++) {
            for (auto i = 0; i < min(size_.x, size_.x); i++) {
                img[{i, j}] = (*this)[{i, j}];
            }
        }
        size = size_;
        swap(pixels, img.pixels);
    }
}

inline void image4b::resize(const vec2i& size_, const vec4b& value) {
    if (size == size_) return;
    if (size == zero2i) {
        *this = image4b{size_, value};
    } else if (size_ == zero2i) {
        clear();
    } else {
        auto img = image4b{size, value};
        for (auto j = 0; j < min(size_.y, size_.y); j++) {
            for (auto i = 0; i < min(size_.x, size_.x); i++) {
                img[{i, j}] = (*this)[{i, j}];
            }
        }
        size = size_;
        swap(pixels, img.pixels);
    }
}

// Size
inline float get_image_aspect(const image4f& img) {
    return (float)img.size.x / (float)img.size.y;
}
inline float get_image_aspect(const image4b& img) {
    return (float)img.size.x / (float)img.size.y;
}

// Gets pixels in an image region
inline image4f get_image_region(const image4f& img, const bbox2i& region) {
    auto clipped = image4f{bbox_size(region)};
    for (auto j = 0; j < bbox_size(region).y; j++) {
        for (auto i = 0; i < bbox_size(region).x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
    return clipped;
}
inline image4b get_image_region(const image4b& img, const bbox2i& region) {
    auto clipped = image4b{bbox_size(region)};
    for (auto j = 0; j < bbox_size(region).y; j++) {
        for (auto i = 0; i < bbox_size(region).x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
    return clipped;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

inline void volume1f::resize(const vec3i& size_, float value) {
    if (size == size_) return;
    if (size == zero3i) {
        *this = volume1f{size_, value};
    } else if (size_ == zero3i) {
        clear();
    } else {
        auto vol = volume1f{size_, value};
        for (auto k = 0; k < min(size_[2], size_[2]); k++) {
            for (auto j = 0; j < min(size_.y, size_.y); j++) {
                for (auto i = 0; i < min(size_.x, size_.x); i++) {
                    vol[{i, j, k}] = (*this)[{i, j, k}];
                }
            }
        }
        size = size_;
        swap(voxels, vol.voxels);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255),
        (byte)clamp(int(a[2] * 256), 0, 255),
        (byte)clamp(int(a[3] * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a[2] / 255.0f, a[3] / 255.0f};
}

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb[2], gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma) {
    return {
        pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin[2], 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma) {
    return {
        pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb[2], gamma), srgb[3]};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma),
        pow(lin[2], 1 / gamma), lin[3]};
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
    return {srgb_to_linear(srgb.x), srgb_to_linear(srgb.y),
        srgb_to_linear(srgb[2])};
}
inline vec3f linear_to_srgb(const vec3f& lin) {
    return {
        linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin[2])};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {srgb_to_linear(srgb.x), srgb_to_linear(srgb.y),
        srgb_to_linear(srgb[2]), srgb[3]};
}
inline vec4f linear_to_srgb(const vec4f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y),
        linear_to_srgb(lin[2]), lin[3]};
}

// Approximate luminance estimate for sRGB primaries (better relative luminance)
inline float luminance(const vec3f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a[2]);
}
inline float luminance(const vec4f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a[2]);
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
    return {tonemap_filmic(hdr.x), tonemap_filmic(hdr.y),
        tonemap_filmic(hdr[2]), hdr[3]};
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_filmic(
    const vec4f& hdr, float exposure, bool filmic, bool srgb) {
    auto scale = pow(2.0f, exposure);
    auto ldr   = vec4f{hdr.x * scale, hdr.y * scale, hdr[2] * scale, hdr[3]};
    if (filmic) ldr = tonemap_filmic(ldr);
    if (srgb) ldr = linear_to_srgb(ldr);
    return ldr;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz[2]),
        xyz.y / (xyz.x + xyz.y + xyz[2]), xyz.y};
}
// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY[2] / xyY.y, xyY[2],
        (1 - xyY.x - xyY.y) * xyY[2] / xyY.y};
}
// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz[2],
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz[2],
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz[2]};
}
// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb[2],
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb[2],
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb[2]};
}

// Convert HSV to RGB
inline vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv[2];
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
    auto  r = rgb.x, g = rgb.y, b = rgb[2];
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
