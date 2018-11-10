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
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container.
template <typename T>
struct image {
    image() : _size{0, 0}, _pixels{} {}
    image(const vec2i& size, const T& value = {})
        : _size{size}, _pixels((size_t)(size.x * size.y), value) {}
    image(const vec2i& size, const T* values)
        : _size{size}, _pixels(values, values + size.x * size.y) {}

    bool  empty() const { return _size.x == 0 && _size.y == 0; }
    vec2i size() const { return _size; }

    void clear() {
        _size = {0, 0};
        _pixels.clear();
    }
    void resize(const vec2i& size, const T& value = {}) {
        if (size == _size) return;
        if (_size == zero2i) {
            *this = image{size, value};
        } else if (size == zero2i) {
            clear();
        } else {
            auto img = image{size, value};
            for (auto j = 0; j < min(size.y, _size.y); j++) {
                for (auto i = 0; i < min(size.x, _size.x); i++) {
                    img[{i, j}] = (*this)[{i, j}];
                }
            }
            _size = size;
            swap(_pixels, img._pixels);
        }
    }

    T& operator[](const vec2i& ij) { return _pixels[ij.y * _size.x + ij.x]; }
    const T& operator[](const vec2i& ij) const {
        return _pixels[ij.y * _size.x + ij.x];
    }
    T&       at(const vec2i& ij) { return _pixels[ij.y * _size.x + ij.x]; }
    const T& at(const vec2i& ij) const {
        return _pixels[ij.y * _size.x + ij.x];
    }

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
inline float get_image_aspect(const image<T>& img) {
    return (float)img.size().x / (float)img.size().y;
}

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

// Make a sunsky HDR model with sun at theta elevation in [0,pif/2], turbidity
// in [1.7,10] with or without sun.
image<vec4f> make_sunsky_image(const vec2i& size, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
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
        : _size{size}, _voxels((size_t)(size.x * size.y * size.z), value) {}
    volume(const vec3i& size, const T* values)
        : _size{size}, _voxels(values, values + size.x * size.y + size.z) {}

    bool  empty() const { return _size.x == 0 && _size.y == 0 && _size.z == 0; }
    vec3i size() const { return _size; }

    void clear() {
        _size = {0, 0, 0};
        _voxels.clear();
    }
    void resize(const vec3i& size, const T& value = {}) {
        if (size == _size) return;
        if (_size == zero3i) {
            *this = volume{size, value};
        } else if (size == zero3i) {
            clear();
        } else {
            auto vol = volume{size, value};
            for (auto k = 0; k < min(size.z, _size.z); k++) {
                for (auto j = 0; j < min(size.y, _size.y); j++) {
                    for (auto i = 0; i < min(size.x, _size.x); i++) {
                        vol[{i, j, k}] = (*this)[{i, j, k}];
                    }
                }
            }
            _size = size;
            swap(_voxels, vol._voxels);
        }
    }

    T& operator[](const vec3i& ijk) {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    const T& operator[](const vec3i& ijk) const {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    T& at(const vec3i& ijk) {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    const T& at(const vec3i& ijk) const {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }

    T*       data() { return _voxels.data(); }
    const T* data() const { return _voxels.data(); }

   private:
    vec3i     _size   = {0, 0};
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
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255),
        (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma = 2.2f) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma), srgb.w};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma = 2.2f) {
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

// Converts HSV to RGB.
vec3f hsv_to_rgb(const vec3f& hsv);
vec3f rgb_to_hsv(const vec3f& rgb);
// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz);
vec3f xyY_to_xyz(const vec3f& xyY);
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz);
vec3f rgb_to_xyz(const vec3f& rgb);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Gets pixels in an image region
template <typename T>
inline void get_image_region(
    image<T>& clipped, const image<T>& img, const bbox2i& region) {
    clipped.resize(bbox_size(region));
    for (auto j = 0; j < bbox_size(region).y; j++) {
        for (auto i = 0; i < bbox_size(region).x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
}
template <typename T>
inline image<T> get_image_region(const image<T>& img, const bbox2i& region) {
    auto clipped = image<T>{bbox_size(region)};
    for (auto j = 0; j < bbox_size(region).y; j++) {
        for (auto i = 0; i < bbox_size(region).x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
    return clipped;
}

}  // namespace yocto

#endif
