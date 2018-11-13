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

#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::vector;

}

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container.
template <typename T>
struct image {
    image() : _size{0, 0}, _pixels{} {}
    image(const vec2i& size, const T& value = {})
        : _size{size}, _pixels((size_t)(size[0] * size[1]), value) {}
    image(const vec2i& size, const T* values)
        : _size{size}, _pixels(values, values + size[0] * size[1]) {}

    bool  empty() const { return _pixels.empty(); }
    vec2i size() const { return _size; }
    int   width() const { return _size[0]; }
    int   height() const { return _size[1]; }

    void clear() {
        _size = {0, 0};
        _pixels.clear();
    }
    void resize(const vec2i& size, const T& value = {});

    T& operator[](const vec2i& ij) { return _pixels[ij[1] * _size[0] + ij[0]]; }
    const T& operator[](const vec2i& ij) const {
        return _pixels[ij[1] * _size[0] + ij[0]];
    }
    T&       at(const vec2i& ij) { return operator[](ij); }
    const T& at(const vec2i& ij) const { return operator[](ij); }

    T*       data() { return _pixels.data(); }
    const T* data() const { return _pixels.data(); }

    T*       begin() { return _pixels.data(); }
    T*       end() { return _pixels.data() + _pixels.size(); }
    const T* begin() const { return _pixels.data(); }

    const T* end() const { return _pixels.data() + _pixels.size(); }

   private:
    vec2i     _size   = {0, 0};
    vector<T> _pixels = {};
};

// Size
template <typename T>
inline float get_image_aspect(const image<T>& img);

// Splits an image into an array of regions
vector<bbox2i> make_image_regions(const vec2i& image_size, int region_size = 32);

// Gets pixels in an image region
template <typename T>
inline void get_image_region(
    image<T>& clipped, const image<T>& img, const bbox2i& region);
template <typename T>
inline image<T> get_image_region(const image<T>& img, const bbox2i& region);

// Gets an image size from a suggested size and an aspect ratio. The suggested
// size may have zeros in either components. In which case, we use the aspect
// ration to compute the other.
vec2i get_image_size(const vec2i& size, float aspect);

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt);
image<vec4b> float_to_byte(const image<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_linear(const image<vec4f>& srgb);
image<vec4f> linear_to_srgb(const image<vec4f>& lin);

// Conversion between linear and gamma-encoded images.
image<vec4f> gamma_to_linear(const image<vec4f>& srgb, float gamma);
image<vec4f> linear_to_gamma(const image<vec4f>& lin, float gamma);

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb);
void tonemap_image_region(image<vec4f>& ldr, const bbox2i& region,
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make example images.
image<vec4f> make_grid_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_checker_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_bumpdimple_image(const vec2i& size, int tile = 8);
image<vec4f> make_ramp_image(
    const vec2i& size, const vec4f& c0, const vec4f& c1, float srgb = false);
image<vec4f> make_gammaramp_image(const vec2i& size);
image<vec4f> make_uvramp_image(const vec2i& size);
image<vec4f> make_uvgrid_image(
    const vec2i& size, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image<vec4f> bump_to_normal_map(const image<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun is not simple to control
// procedurally, so we include two variables to unrealistically customize it:
// sun_size_scale and sun_emission_scale with good values in [20-] for both.
// To make editing easier, we also renormalize the (sun+sky) integral to
// the sky only integral if renormalize_sun is true. For the same reason,
// we rescale the sun dimension such that it covers at least 5 pixels in
// diameter if has_sun is enabled.
image<vec4f> make_sunsky_image(const vec2i& size, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_angle_scale = 1.0f,
    float        sun_emission_scale = 1.0f,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f}, bool renormalize_sun = true);
// Make an image of multiple lights.
image<vec4f> make_lights_image(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image<vec4f> make_noise_image(
    const vec2i& size, float scale = 1, bool wrap = true);
image<vec4f> make_fbm_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image<vec4f> make_ridge_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image<vec4f> make_turbulence_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container.
template <typename T>
struct volume {
    volume() : _size{0, 0, 0}, _voxels{} {}
    volume(const vec3i& size, const T& value = {})
        : _size{size}, _voxels((size_t)(size[0] * size[1] * size[2]), value) {}
    volume(const vec3i& size, const T* values)
        : _size{size}, _voxels(values, values + size[0] * size[1] + size[2]) {}

    bool  empty() const { return _voxels.empty(); }
    vec3i size() const { return _size; }
    int   width() const { return _size[0]; }
    int   height() const { return _size[1]; }
    int   depth() const { return _size[2]; }

    void clear() {
        _size = {0, 0, 0};
        _voxels.clear();
    }
    void resize(const vec3i& size, const T& value = {});

    T& operator[](const vec3i& ijk) {
        return _voxels[ijk[2] * _size[0] * _size[1] + ijk[1] * _size[0] + ijk[0]];
    }
    const T& operator[](const vec3i& ijk) const {
        return _voxels[ijk[2] * _size[0] * _size[1] + ijk[1] * _size[0] + ijk[0]];
    }
    T&       at(const vec3i& ijk) { return operator[](ijk); }
    const T& at(const vec3i& ijk) const { return operator[](ijk); }

    T*       data() { return _voxels.data(); }
    const T* data() const { return _voxels.data(); }

   private:
    vec3i     _size   = {0, 0, 0};
    vector<T> _voxels = {};
};

// make a simple example volume
volume<float> make_test_volume(
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

template <typename T>
inline void image<T>::resize(const vec2i& size, const T& value) {
    if (size == _size) return;
    if (_size == zero_vec2i) {
        *this = image{size, value};
    } else if (size == zero_vec2i) {
        clear();
    } else {
        auto img = image{size, value};
        for (auto j = 0; j < min(size[1], _size[1]); j++) {
            for (auto i = 0; i < min(size[0], _size[0]); i++) {
                img[{i, j}] = (*this)[{i, j}];
            }
        }
        _size = size;
        swap(_pixels, img._pixels);
    }
}

// Size
template <typename T>
inline float get_image_aspect(const image<T>& img) {
    return (float)img.width() / (float)img.height();
}

// Gets pixels in an image region
template <typename T>
inline void get_image_region(
    image<T>& clipped, const image<T>& img, const bbox2i& region) {
    clipped.resize(bbox_size(region));
    for (auto j = 0; j < bbox_size(region)[1]; j++) {
        for (auto i = 0; i < bbox_size(region)[0]; i++) {
            clipped[{i, j}] = img[{i + region.min[0], j + region.min[1]}];
        }
    }
}
template <typename T>
inline image<T> get_image_region(const image<T>& img, const bbox2i& region) {
    auto clipped = image<T>{bbox_size(region)};
    for (auto j = 0; j < bbox_size(region)[1]; j++) {
        for (auto i = 0; i < bbox_size(region)[0]; i++) {
            clipped[{i, j}] = img[{i + region.min[0], j + region.min[1]}];
        }
    }
    return clipped;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline void volume<T>::resize(const vec3i& size, const T& value) {
    if (size == _size) return;
    if (_size == zero_vec3i) {
        *this = volume{size, value};
    } else if (size == zero_vec3i) {
        clear();
    } else {
        auto vol = volume{size, value};
        for (auto k = 0; k < min(size[2], _size[2]); k++) {
            for (auto j = 0; j < min(size[1], _size[1]); j++) {
                for (auto i = 0; i < min(size[0], _size[0]); i++) {
                    vol[{i, j, k}] = (*this)[{i, j, k}];
                }
            }
        }
        _size = size;
        swap(_voxels, vol._voxels);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)clamp(int(a[0] * 256), 0, 255),
        (byte)clamp(int(a[1] * 256), 0, 255),
        (byte)clamp(int(a[2] * 256), 0, 255),
        (byte)clamp(int(a[3] * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a[0] / 255.0f, a[1] / 255.0f, a[2] / 255.0f, a[3] / 255.0f};
}

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma) {
    return {pow(srgb[0], gamma), pow(srgb[1], gamma), pow(srgb[2], gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma) {
    return {
        pow(lin[0], 1 / gamma), pow(lin[1], 1 / gamma), pow(lin[2], 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma) {
    return {
        pow(srgb[0], gamma), pow(srgb[1], gamma), pow(srgb[2], gamma), srgb[3]};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma) {
    return {pow(lin[0], 1 / gamma), pow(lin[1], 1 / gamma),
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
    return {srgb_to_linear(srgb[0]), srgb_to_linear(srgb[1]),
        srgb_to_linear(srgb[2])};
}
inline vec3f linear_to_srgb(const vec3f& lin) {
    return {
        linear_to_srgb(lin[0]), linear_to_srgb(lin[1]), linear_to_srgb(lin[2])};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {srgb_to_linear(srgb[0]), srgb_to_linear(srgb[1]),
        srgb_to_linear(srgb[2]), srgb[3]};
}
inline vec4f linear_to_srgb(const vec4f& lin) {
    return {linear_to_srgb(lin[0]), linear_to_srgb(lin[1]),
        linear_to_srgb(lin[2]), lin[3]};
}

// Approximate luminance estimate for sRGB primaries (better relative luminance)
inline float luminance(const vec3f& a) {
    return (0.2126f * a[0] + 0.7152f * a[1] + 0.0722 * a[2]);
}
inline float luminance(const vec4f& a) {
    return (0.2126f * a[0] + 0.7152f * a[1] + 0.0722 * a[2]);
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
    return {tonemap_filmic(hdr[0]), tonemap_filmic(hdr[1]),
        tonemap_filmic(hdr[2]), hdr[3]};
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_filmic(
    const vec4f& hdr, float exposure, bool filmic, bool srgb) {
    auto scale = pow(2.0f, exposure);
    auto ldr   = vec4f{hdr[0] * scale, hdr[1] * scale, hdr[2] * scale, hdr[3]};
    if (filmic) ldr = tonemap_filmic(ldr);
    if (srgb) ldr = linear_to_srgb(ldr);
    return ldr;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero_vec3f) return zero_vec3f;
    return {xyz[0] / (xyz[0] + xyz[1] + xyz[2]),
        xyz[1] / (xyz[0] + xyz[1] + xyz[2]), xyz[1]};
}
// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY[1] == 0) return zero_vec3f;
    return {xyY[0] * xyY[2] / xyY[1], xyY[2],
        (1 - xyY[0] - xyY[1]) * xyY[2] / xyY[1]};
}
// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero_vec3f) return zero_vec3f;
    return {+3.2404542f * xyz[0] - 1.5371385f * xyz[1] - 0.4985314f * xyz[2],
        -0.9692660f * xyz[0] + 1.8760108f * xyz[1] + 0.0415560f * xyz[2],
        +0.0556434f * xyz[0] - 0.2040259f * xyz[1] + 1.0572252f * xyz[2]};
}
// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero_vec3f) return zero_vec3f;
    return {0.4124564f * rgb[0] + 0.3575761f * rgb[1] + 0.1804375f * rgb[2],
        0.2126729f * rgb[0] + 0.7151522f * rgb[1] + 0.0721750f * rgb[2],
        0.0193339f * rgb[0] + 0.1191920f * rgb[1] + 0.9503041f * rgb[2]};
}

// Convert HSV to RGB
inline vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv[0], s = hsv[1], v = hsv[2];
    if (hsv[1] == 0.0f) return {v, v, v};

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
    auto  r = rgb[0], g = rgb[1], b = rgb[2];
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
