//
// # Yocto/Image: Tiny imaging Library mostly for rendering and color support
//
//
// Yocto/Image is a collection of image utilities useful when writing rendering
// algorithms. These include a simple image data structure, color conversion
// utilities and tone mapping. We provinde loading and saving functionality for
// images and support PNG, JPG, TGA, BMP, HDR, EXR formats.
//
// This library depends on stb_image.h, stb_image_write.h, stb_image_resize.h,
// tinyexr.h for the IO features. If thoese are not needed, it can be safely
// used without dependencies.
//
//
// ## Image Utilities
//
// Yocto/Image supports a very small set of color and image utilities including
// color utilities, example image creation, tone mapping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// global illumination renderer, rather than the need of generic image editing.
// We support 4-channels float images (assumed to be in linear color) and
// 4-channels byte images (assumed to be in sRGB).
//
//
// 1. store images using the image<T> structure
// 2. load and save images with `load_image()` and `save_image()`
// 3. resize images with `resize_image()`
// 4. tonemap images with `tonemap_image()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky_image()`
// 7. load and save images with Yocto/ImageIO
// 8. many color conversion functions are available in the code below
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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
//  LICENSE for blackbody code
//
// Copyright (c) 2015 Neil Bartlett
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
//
//

#ifndef _YOCTO_IMAGE_H_
#define _YOCTO_IMAGE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_random.h"
#include "yocto_utils.h"

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container.
template <typename T>
struct image {
    // constructors
    image() : _size{0, 0}, _pixels{} {}
    image(const vec2i& size, const T& value = {})
        : _size{size}, _pixels((size_t)size.x * (size_t)size.y, value) {}
    image(const vec2i& size, const T* value)
        : _size{size}
        , _pixels(value, value + (size_t)size.x * (size_t)size.y) {}

    // size
    bool  empty() const { return _pixels.empty(); }
    vec2i size() const { return _size; }
    bool  contains(const vec2i& ij) const {
        return ij.x > 0 && ij.x < _size.x && ij.y > 0 && ij.y < _size.y;
    }
    void resize(const vec2i& size) {
        if (size == _size) return;
        _size = size;
        _pixels.resize((size_t)size.x * (size_t)size.y);
    }

    // element access
    T& operator[](const vec2i& ij) { return _pixels[ij.y * _size.x + ij.x]; }
    const T& operator[](const vec2i& ij) const {
        return _pixels[ij.y * _size.x + ij.x];
    }

    // data access
    T*       data() { return _pixels.data(); }
    const T* data() const { return _pixels.data(); }

    // iteration
    T*       begin() { return _pixels.data(); }
    T*       end() { return _pixels.data() + _pixels.size(); }
    const T* begin() const { return _pixels.data(); }
    const T* end() const { return _pixels.data() + _pixels.size(); }

    // data
    vec2i     _size   = zero2i;
    vector<T> _pixels = {};
};

// equality
template <typename T>
inline bool operator==(const image<T>& a, const image<T>& b) {
    return a.size() == b.size() && a._pixels == b._pixels;
}
template <typename T>
inline bool operator!=(const image<T>& a, const image<T>& b) {
    return a.size() != b.size() || a._pixels != b._pixels;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
inline bool is_hdr_filename(const string& filename) {
    return get_extension(filename) == "hdr" ||
           get_extension(filename) == "exr" || get_extension(filename) == "pfm";
}
// Return the image preset type and the remaining filename
inline bool is_image_preset_filename(const string& filename) {
    return get_extension(filename) == "ypreset";
}
inline pair<string, string> get_image_preset_type(const string& filename) {
    if (get_extension(filename) == "ypreset") {
        return {get_noextension(get_noextension(get_filename(filename))),
            get_noextension(filename)};
    } else {
        return {"", filename};
    }
}

// Loads/saves a 1-4 channels float image in linear color space.
void load_image(const string& filename, image<float>& img);
void load_image(const string& filename, image<vec2f>& img);
void load_image(const string& filename, image<vec3f>& img);
void load_image(const string& filename, image<vec4f>& img);
void save_image(const string& filename, const image<float>& img);
void save_image(const string& filename, const image<vec2f>& img);
void save_image(const string& filename, const image<vec3f>& img);
void save_image(const string& filename, const image<vec4f>& img);

// Loads/saves a 1-4 byte image in sRGB color space.
void load_image(const string& filename, image<byte>& img);
void load_image(const string& filename, image<vec2b>& img);
void load_image(const string& filename, image<vec3b>& img);
void load_image(const string& filename, image<vec4b>& img);
void save_image(const string& filename, const image<byte>& img);
void save_image(const string& filename, const image<vec2b>& img);
void save_image(const string& filename, const image<vec3b>& img);
void save_image(const string& filename, const image<vec4b>& img);

// Convenience helper for loading HDR or LDR based on filename
template <int N>
inline void load_image(const string& filename, image<vec<float, N>>& hdr,
    image<vec<byte, N>>& ldr) {
    if (is_hdr_filename(filename)) {
        load_image(filename, hdr);
    } else {
        load_image(filename, ldr);
    }
}
template <int N>
inline void save_image(const string& filename, const image<vec<float, N>>& hdr,
    const image<vec<byte, N>>& ldr) {
    if (!hdr.empty()) {
        save_image(filename, hdr);
    } else {
        save_image(filename, ldr);
    }
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
template <int N>
inline void save_tonemapped_image(const string& filename,
    const image<vec<float, N>>& hdr, float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        save_image(filename, hdr);
    } else {
        auto ldr = image<vec<byte, N>>{hdr.size()};
        tonemap_image8(ldr, hdr, exposure, filmic, srgb);
        save_image(filename, ldr);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image region
struct image_region {
    vec2i min = zero2i;
    vec2i max = zero2i;
    vec2i size() const { return max - min; }
};

// Splits an image into an array of regions
inline void make_image_regions(vector<image_region>& regions, const vec2i& size,
    int region_size = 32, bool shuffled = false);

// Gets pixels in an image region
template <typename T>
inline void get_image_region(
    image<T>& clipped, const image<T>& img, const image_region& region);
template <typename T>
inline void set_image_region(
    image<T>& img, const image<T>& region, const vec2i& offset);

// Apply a function to each image pixel
template <typename T1, typename T2, typename Func>
inline void apply(
    const Func& func, image<T1>& result, const image<T2>& source) {
    result.resize(source.size());
    for (auto j = 0; j < result.size().y; j++) {
        for (auto i = 0; i < result.size().x; i++) {
            result[{i, j}] = func(source[{i, j}]);
        }
    }
}
template <typename T1, typename T2, typename Func>
inline void apply(const Func& func, image<T1>& result, const image<T2>& source,
    const image_region& region) {
    result.resize(source.size());
    for (auto j = region.min.y; j < region.max.y; j++) {
        for (auto i = region.min.x; i < region.max.x; i++) {
            result[{i, j}] = func(source[{i, j}]);
        }
    }
}

// Conversion from/to floats.
template <typename T, typename TB>
inline void byte_to_float(image<T>& fl, const image<TB>& bt) {
    return apply([](auto& a) { return byte_to_float(a); }, fl, bt);
}
template <typename T, typename TB>
inline void float_to_byte(image<TB>& bt, const image<T>& fl) {
    return apply([](auto& a) { return float_to_byte(a); }, bt, fl);
}

// Conversion between linear and gamma-encoded images.
template <typename T>
inline void srgb_to_linear(image<T>& lin, const image<T>& srgb) {
    return apply([](auto& a) { return srgb_to_linear(a); }, lin, srgb);
}
template <typename T>
inline void linear_to_srgb(image<T>& srgb, const image<T>& lin) {
    return apply([](auto& a) { return linear_to_srgb(a); }, srgb, lin);
}
template <typename T, typename TB>
inline void srgb8_to_linear(image<T>& lin, const image<TB>& srgb) {
    return apply(
        [](auto& a) { return srgb_to_linear(byte_to_float(a)); }, lin, srgb);
}
template <typename T, typename TB>
inline void linear_to_srgb8(image<TB>& srgb, const image<T>& lin) {
    return apply(
        [](auto& a) { return float_to_byte(linear_to_srgb(a)); }, srgb, lin);
}

// Conversion between linear and gamma-encoded images.
template <typename T1, typename TC>
inline void gamma_to_linear(image<T1>& lin, const image<T1>& srgb, TC gamma) {
    return apply(
        [gamma](auto& a) { return gamma_to_linear(a, gamma); }, lin, srgb);
}
template <typename T1, typename TC>
inline void linear_to_gamma(image<T1>& srgb, const image<T1>& lin, TC gamma) {
    return apply(
        [gamma](auto& a) { return linear_to_gamma(a, gamma); }, srgb, lin);
}

// Apply exposure and filmic tone mapping
template <typename T, typename TC>
inline void tonemap_image(
    image<T>& ldr, const image<T>& hdr, TC exposure, bool filmic, bool srgb) {
    return apply(
        [exposure, filmic, srgb](
            auto& a) { return tonemap_filmic(a, exposure, filmic, srgb); },
        ldr, hdr);
}
template <typename T, typename TB, typename TC>
inline void tonemap_image8(
    image<TB>& ldr, const image<T>& hdr, TC exposure, bool filmic, bool srgb) {
    return apply(
        [exposure, filmic, srgb](auto& a) {
            return float_to_byte(tonemap_filmic(a, exposure, filmic, srgb));
        },
        ldr, hdr);
}
template <typename T, typename TC>
inline void tonemap_image_region(image<T>& ldr, const image_region& region,
    const image<T>& hdr, TC exposure, bool filmic, bool srgb) {
    return apply(
        [exposure, filmic, srgb](
            auto& a) { return tonemap_filmic(a, exposure, filmic, srgb); },
        ldr, hdr, region);
}

// Resize an image.
template <typename T, int N>
inline void resize_image(image<vec<T, N>>& res, const image<vec<T, N>>& img);
inline void resize_image(image<float>& res, const image<float>& img) {
    return resize_image((image<vec1f>&)res, (const image<vec1f>&)img);
}
inline void resize_image(image<byte>& res, const image<byte>& img) {
    return resize_image((image<vec1b>&)res, (const image<vec1b>&)img);
}
template <typename T>
inline void resize_image(image<T>& res, const image<T>& img, const vec2i& size);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make example images in linear color space. Takes as input images allocated
// to the desired size and fill the pixel with expected values.
template <typename T>
inline void make_grid_image(
    image<T>& img, const vec2i& size, int tile, const T& c0, const T& c1);
template <typename T>
inline void make_checker_image(
    image<T>& img, const vec2i& size, int tile, const T& c0, const T& c1);
template <typename T>
inline void make_bumpdimple_image(
    image<T>& img, const vec2i& size, int tile, const T& c0, const T& c1);
template <typename T>
inline void make_ramp_image(
    image<T>& img, const vec2i& size, const T& c0, const T& c1);
template <typename T>
inline void make_ramp_image(image<T>& img, const vec2i& size, const T& c00,
    const T& c10, const T& c11, const T& c01);
template <typename T>
inline void make_gammaramp_image(
    image<T>& img, const vec2i& size, const T& c0, const T& c1);
template <typename T>
inline void make_uvramp_image(image<T>& img, const vec2i& size);
template <typename T, int N>
inline void make_uvgrid_image(image<vec<T, N>>& img, const vec2i& size,
    int tile = 8, bool colored = true);
template <typename T, int N>
inline void make_blackbodyramp_image(image<vec<T, N>>& img, const vec2i& size,
    float start_temperature = 1000, float end_temperature = 12000);

// Comvert a bump map to a normal map. All linear color spaces.
template <typename T, int N>
inline void bump_to_normal_map(
    image<vec<T, N>>& norm, const image<vec<T, N>>& img, T scale = 1);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
template <typename T, int N>
inline void make_sunsky_image(image<vec<T, N>>& img, const vec2i& size,
    T sun_angle, T turbidity = 3, bool has_sun = false, T sun_intensity = 1,
    T                sun_temperature = 0,
    const vec<T, 3>& ground_albedo   = {(T)0.2, (T)0.2, (T)0.2});
// Make an image of multiple lights.
template <typename T, int N>
inline void make_lights_image(image<vec<T, N>>& img, const vec2i& size,
    const vec<T, 3>& le = {1, 1, 1}, int nlights = 4, T langle = (T)pi / 4,
    T lwidth = (T)pi / 16, T lheight = (T)pi / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
template <typename T, typename T1>
inline void make_noise_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale = 1, bool wrap = true);
template <typename T, typename T1>
inline void make_fbm_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale = 1, T1 lacunarity = 2, T1 gain = 0.5f,
    int octaves = 6, bool wrap = true);
template <typename T, typename T1>
inline void make_ridge_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale = 1, T1 lacunarity = 2, T1 gain = (T1)0.5,
    T1 offset = 1, int octaves = 6, bool wrap = true);
template <typename T, typename T1>
inline void make_turbulence_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale = 1, T1 lacunarity = 2, T1 gain = (T1)0.5,
    int octaves = 6, bool wrap = true);

// Add a border to an image
template <typename T>
inline void add_image_border(
    image<T>& img, const vec2i& size, int border_width, const T& border_color);

// Make logo images. Image is resized to proper size.
template <typename T>
inline void make_logo_image(image<T>& img, const string& name);

// Make an image preset, useful for testing. See implementation for types.
inline void make_image_preset(image<vec<float, 4>>& img, const string& type);
inline void make_image_preset(image<vec<byte, 4>>& img, const string& type);
inline void make_image_preset(
    image<vec4f>& hdr, image<vec4b>& ldr, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container.
template <typename T>
struct volume {
    // constructors
    volume() : _size{0, 0, 0}, _voxels{} {}
    volume(const vec3i& size, const T& value = {})
        : _size{size}
        , _voxels((size_t)size.x * (size_t)size.y * (size_t)size.z, value) {}
    volume(const vec3i& size, const T* value)
        : _size{size}
        , _voxels(value,
              value + (size_t)size.x * (size_t)size.y * (size_t)size.z) {}

    // size
    bool  empty() const { return _voxels.empty(); }
    vec3i size() const { return _size; }
    void  resize(const vec3i& size) {
        if (size == _size) return;
        _size = size;
        _voxels.resize((size_t)size.x * (size_t)size.y * (size_t)size.z);
    }

    // element access
    T& operator[](const vec3i& ijk) {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    const T& operator[](const vec3i& ijk) const {
        return _voxels[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }

    // data access
    T*       data() { return _voxels.data(); }
    const T* data() const { return _voxels.data(); }

    // iteration
    T*       begin() { return _voxels.data(); }
    T*       end() { return _voxels.data() + _voxels.size(); }
    const T* begin() const { return _voxels.data(); }
    const T* end() const { return _voxels.data() + _voxels.size(); }

    // data
    vec3i         _size   = zero3i;
    vector<float> _voxels = {};
};

// equality
template <typename T>
inline bool operator==(const volume<T>& a, const volume<T>& b) {
    return a.size() == b.size() && a._voxels == b._voxels;
}
template <typename T>
inline bool operator!=(const volume<T>& a, const volume<T>& b) {
    return a.size() != b.size() && a._voxels != b._voxels;
}

// make a simple example volume
inline void make_test_volume(volume<float>& vol, const vec3i& size,
    float scale = 10, float exponent = 6);
inline void make_volume_preset(volume<float>& vol, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is a preset based on filename.
inline bool is_volume_preset_filename(const string& filename) {
    return get_extension(filename) == "ypreset";
}

// Loads/saves a 1 channel volume.
void load_volume(const string& filename, volume<float>& vol);
void save_volume(const string& filename, const volume<float>& vol);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

inline byte float_to_byte(float a) { return (byte)clamp(int(a * 256), 0, 255); }
inline float byte_to_float(byte a) { return a / 255.0f; }

// Element-wise float to byte conversion.
template <int N>
inline vec<byte, N> float_to_byte(const vec<float, N>& a) {
    if constexpr (N == 1) {
        return {float_to_byte(a.x)};
    } else if constexpr (N == 2) {
        return {float_to_byte(a.x), float_to_byte(a.y)};
    } else if constexpr (N == 3) {
        return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z)};
    } else if constexpr (N == 4) {
        return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z),
            float_to_byte(a.w)};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}
template <int N>
inline vec<float, N> byte_to_float(const vec<byte, N>& a) {
    if constexpr (N == 1) {
        return {byte_to_float(a.x)};
    } else if constexpr (N == 2) {
        return {byte_to_float(a.x), byte_to_float(a.y)};
    } else if constexpr (N == 3) {
        return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z)};
    } else if constexpr (N == 4) {
        return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z),
            byte_to_float(a.w)};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}

// Default alpha
template <typename T>
constexpr T _default_alpha() {
    if constexpr (std::is_same_v<T, byte>) {
        return (byte)255;
    } else {
        return (T)1;
    }
}
template <typename T>
constexpr T default_alpha = _default_alpha<T>();

// Apply an operator to a color
template <typename T>
inline vec<T, 3> color_to_rgb(T a) {
    return {a, a, a};
}
template <typename T, int N>
inline vec<T, 3> color_to_rgb(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return {a.x, a.x, a.x};
    } else if constexpr (N == 2) {
        return {a.x, a.x, a.x};
    } else if constexpr (N == 3) {
        return {a.x, a.y, a.z};
    } else if constexpr (N == 4) {
        return {a.x, a.y, a.z};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}
template <typename T>
inline vec<T, 3> color_to_rgba(T a) {
    return {a, a, a, default_alpha<T>};
}
template <typename T, int N>
inline vec<T, 4> color_to_rgba(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return {a.x, a.x, a.x, default_alpha<T>};
    } else if constexpr (N == 2) {
        return {a.x, a.x, a.x, a.y};
    } else if constexpr (N == 3) {
        return {a.x, a.y, a.z, default_alpha<T>};
    } else if constexpr (N == 4) {
        return {a.x, a.y, a.z, a.w};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}
template <typename T, int N>
inline T luminance(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return a.x;
    } else if constexpr (N == 2) {
        return a.x;
    } else if constexpr (N == 3) {
        return ((T)0.2126 * a.x + (T)0.7152 * a.y + (T)0.0722 * a.z);
    } else if constexpr (N == 4) {
        return ((T)0.2126 * a.x + (T)0.7152 * a.y + (T)0.0722 * a.z);
    } else {
        throw runtime_error("Bad number of arguments");
    }
}
template <int N>
inline byte luminance(const vec<byte, N>& a) {
    return float_to_byte(luminance(byte_to_float(a)));
}

template <typename T, int N>
inline vec<T, N> gray_to_color(T a) {
    if constexpr (N == 1) {
        return a;
    } else if constexpr (N == 2) {
        return {a, default_alpha<T>};
    } else if constexpr (N == 3) {
        return {a, a, a};
    } else if constexpr (N == 4) {
        return {a, a, a, default_alpha<T>};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}
template <typename T, int N>
inline vec<T, N> rgba_to_color(const vec<T, 4>& a) {
    if constexpr (N == 1) {
        return {luminance(a)};
    } else if constexpr (N == 2) {
        return {luminance(a), a.y};
    } else if constexpr (N == 3) {
        return {a.x, a.y, a.z};
    } else if constexpr (N == 4) {
        return {a.x, a.y, a.z, a.w};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}

template <typename T, int N>
inline void gray_to_color(image<vec<T, N>>& col, const image<T>& gray) {
    return apply([](auto& a) { return gray_to_color<T, N>(a); }, col, gray);
}
template <typename T, int N>
inline void rgba_to_color(image<vec<T, N>>& col, const image<vec<T, 4>>& rgba) {
    return apply([](auto& a) { return rgba_to_color<T, N>(a); }, col, rgba);
}

// Apply an operator to a color
template <typename T, int N, typename Func>
inline vec<T, N> apply_color(const Func& func, const vec<T, N>& a) {
    if constexpr (N == 1) {
        return {func(a.x)};
    } else if constexpr (N == 2) {
        return {func(a.x), a.y};
    } else if constexpr (N == 3) {
        return {func(a.x), func(a.y), func(a.z)};
    } else if constexpr (N == 4) {
        return {func(a.x), func(a.y), func(a.z), a.w};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}

// Lerp colors between two values
template <typename T, typename T1>
inline T lerp_color(const T& a, const T& b, T1 u) {
    return lerp(a, b, u);
}
template <typename T1>
inline byte lerp_color(byte a, byte b, T1 u) {
    return float_to_byte(lerp(byte_to_float(a), byte_to_float(b), u));
}
template <int N, typename T1>
inline vec<byte, N> lerp_color(
    const vec<byte, N>& a, const vec<byte, N> b, T1 u) {
    return float_to_byte(lerp(byte_to_float(a), byte_to_float(b), u));
}

template <typename T, typename T1>
inline T bilerp_color(
    const T& c00, const T& c10, const T& c11, const T& c01, T1 u, T1 v) {
    return bilerp(c00, c01, c11, c01, u, v);
}
template <typename T1>
inline byte bilerp_color(byte c00, byte c10, byte c11, byte c01, T1 u, T1 v) {
    return float_to_byte(bilerp(byte_to_float(c00), byte_to_float(c01),
        byte_to_float(c11), byte_to_float(c01), u, v));
}
template <int N, typename T1>
inline vec<byte, N> bilerp_color(const vec<byte, N>& c00,
    const vec<byte, N> c10, const vec<byte, N>& c11, const vec<byte, N> c01,
    T1 u, T1 v) {
    return float_to_byte(bilerp(byte_to_float(c00), byte_to_float(c01),
        byte_to_float(c11), byte_to_float(c01), u, v));
}

// Conversion between linear and gamma-encoded colors.
inline float gamma_to_linear(float srgb, float gamma) {
    return pow(srgb, gamma);
}
inline float linear_to_gamma(float srgb, float gamma) {
    return pow(srgb, 1 / gamma);
}
template <typename T, int N>
inline vec<T, N> gamma_to_linear(const vec<T, N>& srgb, float gamma) {
    return apply_color(
        [&gamma](auto& a) { return gamma_to_linear(a, gamma); }, srgb);
}
template <typename T, int N>
inline vec<T, N> linear_to_gamma(const vec<T, N>& lin, float gamma) {
    return apply_color(
        [&gamma](auto& a) { return linear_to_gamma(a, gamma); }, lin);
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
template <typename T, int N>
inline vec<T, N> srgb_to_linear(const vec<T, N>& srgb) {
    return apply_color([](auto& a) { return srgb_to_linear(a); }, srgb);
}
template <typename T, int N>
inline vec<T, N> linear_to_srgb(const vec<T, N>& lin) {
    return apply_color([](auto& a) { return linear_to_srgb(a); }, lin);
}

// Fitted ACES tonemapping curve.
template <typename T>
inline T tonemap_filmic(T hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    // hdr *= 0.6; // brings it back to ACES range
    return (hdr * hdr * (T)2.51 + hdr * (T)0.03) /
           (hdr * hdr * (T)2.43 + hdr * (T)0.59 + (T)0.14);
}
template <typename T, int N>
inline vec<T, N> tonemap_filmic(const vec<T, N>& hdr) {
    return apply_color([](auto& a) { return tonemap_filmic(a); }, hdr);
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
template <typename T, int N>
inline vec<T, N> tonemap_filmic(
    const vec<T, N>& hdr, T exposure, bool filmic, bool srgb) {
    if constexpr (N == 1) {
        auto scale = pow(2.0f, exposure);
        auto ldr   = hdr * scale;
        if (filmic) ldr = tonemap_filmic(ldr);
        if (srgb) ldr = linear_to_srgb(ldr);
        return ldr;
    } else if constexpr (N == 2) {
        auto scale = pow(2.0f, exposure);
        auto ldr   = hdr.x * scale;
        if (filmic) ldr = tonemap_filmic(ldr);
        if (srgb) ldr = linear_to_srgb(ldr);
        return {ldr, hdr.y};
    } else if constexpr (N == 3) {
        auto scale = pow(2.0f, exposure);
        auto ldr   = hdr * scale;
        if (filmic) ldr = tonemap_filmic(ldr);
        if (srgb) ldr = linear_to_srgb(ldr);
        return ldr;
    } else if constexpr (N == 4) {
        auto scale = pow(2.0f, exposure);
        auto ldr   = hdr.xyz * scale;
        if (filmic) ldr = tonemap_filmic(ldr);
        if (srgb) ldr = linear_to_srgb(ldr);
        return {ldr, hdr.w};
    } else {
        throw runtime_error("Bad number of arguments");
    }
}

// Convert between CIE XYZ and xyY
template <typename T>
inline vec<T, 3> xyz_to_xyY(const vec<T, 3>& xyz) {
    if (xyz == zero<T, 3>) return zero<T, 3>;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
template <typename T>
inline vec<T, 3> xyY_to_xyz(const vec<T, 3>& xyY) {
    if (xyY.y == 0) return zero<T, 3>;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
template <typename T>
inline vec<T, 3> xyz_to_rgb(const vec<T, 3>& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    return {
        (T) + 3.2404542 * xyz.x - (T)1.5371385 * xyz.y - (T)0.4985314 * xyz.z,
        (T)-0.9692660 * xyz.x + (T)1.8760108 * xyz.y + (T)0.0415560 * xyz.z,
        (T) + 0.0556434 * xyz.x - (T)0.2040259 * xyz.y + (T)1.0572252 * xyz.z,
    };
}
// Convert between CIE XYZ and RGB
template <typename T>
inline vec<T, 3> rgb_to_xyz(const vec<T, 3>& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    return {
        (T)0.4124564 * rgb.x + (T)0.3575761 * rgb.y + (T)0.1804375 * rgb.z,
        (T)0.2126729 * rgb.x + (T)0.7151522 * rgb.y + (T)0.0721750 * rgb.z,
        (T)0.0193339 * rgb.x + (T)0.1191920 * rgb.y + (T)0.9503041 * rgb.z,
    };
}

// Approximate color of blackbody radiation from wavelength in nm.
template <typename T>
inline vec<T, 3> blackbody_to_rgb(T temperature);

// Converts HSV to RGB.
template <typename T>
inline vec<T, 3> hsv_to_rgb(const vec<T, 3>& hsv);
template <typename T>
inline vec<T, 3> rgb_to_hsv(const vec<T, 3>& rgb);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BUILTIN IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Save with a logo embedded
template <typename T, int N>
inline void save_image_with_logo(
    const string& filename, const image<vec<T, N>>& img) {
    auto logo = image<vec<T, N>>{};
    make_logo_image(logo, "logo-render");
    auto img_copy = img;
    auto offset   = img.size() - logo.size() - 8;
    set_image_region(img_copy, logo, offset);
    save_image(filename, img_copy);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
template <int N>
inline void save_tonemapped_image_with_logo(const string& filename,
    const image<vec<float, N>>& hdr, float exposure, bool filmic, bool srgb) {
    if (is_hdr_filename(filename)) {
        save_image_with_logo(filename, hdr);
    } else {
        auto ldr = image<vec<byte, N>>{hdr.size()};
        tonemap_image8(ldr, hdr, exposure, filmic, srgb);
        save_image_with_logo(filename, ldr);
    }
}

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

#include "ext/ArHosekSkyModel.h"
#include "ext/stb_image_resize.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Approximate color of blackbody radiation from wavelength in nm.
template <typename T>
inline vec<T, 3> blackbody_to_rgb(T temperature) {
    // https://github.com/neilbartlett/color-temperature
    auto rgb = zero<T, 3>;
    if ((temperature / 100) < 66) {
        rgb.x = 255;
    } else {
        // a + b x + c Log[x] /.
        // {a -> 351.97690566805693`,
        // b -> 0.114206453784165`,
        // c -> -40.25366309332127
        // x -> (kelvin/100) - 55}
        rgb.x = (temperature / 100) - 55;
        rgb.x = (T)351.97690566805693 + (T)0.114206453784165 * rgb.x -
                (T)40.25366309332127 * log(rgb.x);
        if (rgb.x < 0) rgb.x = 0;
        if (rgb.x > 255) rgb.x = 255;
    }

    if ((temperature / 100) < 66) {
        // a + b x + c Log[x] /.
        // {a -> -155.25485562709179`,
        // b -> -0.44596950469579133`,
        // c -> 104.49216199393888`,
        // x -> (kelvin/100) - 2}
        rgb.y = (temperature / 100) - 2;
        rgb.y = (T)-155.25485562709179 - (T)0.44596950469579133 * rgb.y +
                (T)104.49216199393888 * log(rgb.y);
        if (rgb.y < 0) rgb.y = 0;
        if (rgb.y > 255) rgb.y = 255;
    } else {
        // a + b x + c Log[x] /.
        // {a -> 325.4494125711974`,
        // b -> 0.07943456536662342`,
        // c -> -28.0852963507957`,
        // x -> (kelvin/100) - 50}
        rgb.y = (temperature / 100) - 50;
        rgb.y = (T)325.4494125711974 + (T)0.07943456536662342 * rgb.y -
                (T)28.0852963507957 * log(rgb.y);
        if (rgb.y < 0) rgb.y = 0;
        if (rgb.y > 255) rgb.y = 255;
    }

    if ((temperature / 100) >= 66) {
        rgb.z = 255;
    } else {
        if ((temperature / 100) <= 20) {
            rgb.z = 0;
        } else {
            // a + b x + c Log[x] /.
            // {a -> -254.76935184120902`,
            // b -> 0.8274096064007395`,
            // c -> 115.67994401066147`,
            // x -> kelvin/100 - 10}
            rgb.z = (temperature / 100) - 10;
            rgb.z = (T)-254.76935184120902 + (T)0.8274096064007395 * rgb.z +
                    (T)115.67994401066147 * log(rgb.z);
            if (rgb.z < 0) rgb.z = 0;
            if (rgb.z > 255) rgb.z = 255;
        }
    }

    return srgb_to_linear(rgb / 255);
}

// Convert HSV to RGB
template <typename T>
inline vec<T, 3> hsv_to_rgb(const vec<T, 3>& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0) return {v, v, v};

    h       = fmod(h, (T)1) / ((T)60 / (T)360);
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
template <typename T>
inline vec<T, 3> rgb_to_hsv(const vec<T, 3>& rgb) {
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
    return {fabsf(K + (g - b) / (6.f * chroma + (T)1e-20)),
        chroma / (r + (T)1e-20), r};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Gets pixels in an image region
template <typename T>
inline void get_image_region(
    image<T>& clipped, const image<T>& img, const image_region& region) {
    clipped.resize(region.size());
    for (auto j = 0; j < region.size().y; j++) {
        for (auto i = 0; i < region.size().x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
}
template <typename T>
inline void set_image_region(
    image<T>& img, const image<T>& region, const vec2i& offset) {
    for (auto j = 0; j < region.size().y; j++) {
        for (auto i = 0; i < region.size().x; i++) {
            if (!img.contains({i, j})) continue;
            img[vec2i{i, j} + offset] = region[{i, j}];
        }
    }
}

// Splits an image into an array of regions
inline void make_image_regions(vector<image_region>& regions, const vec2i& size,
    int region_size, bool shuffled) {
    regions.clear();
    for (auto y = 0; y < size.y; y += region_size) {
        for (auto x = 0; x < size.x; x += region_size) {
            regions.push_back({{x, y},
                {min(x + region_size, size.x), min(y + region_size, size.y)}});
        }
    }
    if (shuffled) {
        auto rng = rng_state{};
        random_shuffle(regions, rng);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Resize image.
template <typename T, int N>
inline void resize_image(
    image<vec<T, N>>& res_img, const image<vec<T, N>>& img) {
    auto alpha = (N == 2 || N == 4) ? N - 1 : -1;
    if constexpr (std::is_same_v<T, float>) {
        stbir_resize_float_generic((T*)img.data(), img.size().x, img.size().y,
            sizeof(vec<T, N>) * img.size().x, (T*)res_img.data(),
            res_img.size().x, res_img.size().y,
            sizeof(vec<T, N>) * res_img.size().x, N, alpha, 0, STBIR_EDGE_CLAMP,
            STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    } else if constexpr (std::is_same_v<T, byte>) {
        stbir_resize_uint8_generic((T*)img.data(), img.size().x, img.size().y,
            sizeof(vec<T, N>) * img.size().x, (T*)res_img.data(),
            res_img.size().x, res_img.size().y,
            sizeof(vec<T, N>) * res_img.size().x, N, alpha, 0, STBIR_EDGE_CLAMP,
            STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    } else {
        throw runtime_error("type not supported");
    }
}
template <typename T>
inline void resize_image(
    image<T>& res_img, const image<T>& img, const vec2i& size_) {
    auto size = size_;
    if (size == zero2i) {
        throw std::invalid_argument("bad image size in resize_image");
    }
    if (size.y == 0) {
        size.y = (int)round(size.x * (float)img.size().y / (float)img.size().x);
    } else if (size.x == 0) {
        size.x = (int)round(size.y * (float)img.size().x / (float)img.size().y);
    }
    res_img = {size};
    resize_image(res_img, img);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image by assign values to each pixel
template <typename T, typename Func>
inline void make_image_fromij(
    image<T>& img, const vec2i& size, const Func& func) {
    img.resize(size);
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            img[{i, j}] = func(i, j);
        }
    }
}
template <typename T, typename Func>
inline void make_image_fromuv(
    image<T>& img, const vec2i& size, const Func& func) {
    img.resize(size);
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto u      = (float)i / (float)img.size().x;
            auto v      = (float)j / (float)img.size().y;
            img[{i, j}] = func(u, v);
        }
    }
}

// Make a grid image
template <typename T>
inline void make_grid_image(
    image<T>& img, const vec2i& size, int tiles, const T& c0, const T& c1) {
    make_image_fromij(
        img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            return (c) ? c0 : c1;
        });
}

// Make a checkerboard image
template <typename T>
inline void make_checker_image(
    image<T>& img, const vec2i& size, int tiles, const T& c0, const T& c1) {
    make_image_fromij(
        img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
            auto c = (i / tile + j / tile) % 2 == 0;
            return (c) ? c0 : c1;
        });
}

// Make an image with bumps and dimples.
template <typename T>
inline void make_bumpdimple_image(
    image<T>& img, const vec2i& size, int tiles, const T& c0, const T& c1) {
    make_image_fromij(
        img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
            auto c  = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r = sqrt(float(ii * ii + jj * jj)) /
                     sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) {
                h += (c) ? (0.5f - r) : -(0.5f - r);
            }
            return lerp_color(c0, c1, h);
        });
}

// Make a uv colored grid
template <typename T>
inline void make_ramp_image(
    image<T>& img, const vec2i& size, const T& c0, const T& c1) {
    make_image_fromij(img, size, [size = img.size(), &c0, &c1](int i, int j) {
        auto u = (float)i / (float)size.x;
        return lerp_color(c0, c1, u);
    });
}
template <typename T>
inline void make_ramp_image(image<T>& img, const vec2i& size, const T& c00,
    const T& c10, const T& c11, const T& c01) {
    make_image_fromij(img, size, [size, &c00, &c10, &c01, &c11](int i, int j) {
        auto u = (float)i / (float)size.x;
        auto v = (float)j / (float)size.y;
        return bilerp_color(c00, c10, c11, c01, u, v);
    });
}

// Make a gamma ramp image
template <typename T>
inline void make_gammaramp_image(
    image<T>& img, const vec2i& size, const T& c0, const T& c1) {
    make_image_fromij(img, size, [size, &c0, &c1](int i, int j) {
        auto u = j / float(size.y - 1);
        if (i < size.x / 3) u = pow(u, 2.2f);
        if (i > (size.x * 2) / 3) u = pow(u, 1 / 2.2f);
        return lerp_color(c0, c1, u);
    });
}

// Make an image color with red/green in the [0,1] range. Helpful to
// visualize uv texture coordinate application.
template <typename T, int N>
inline void make_uvramp_image(image<vec<T, N>>& img, const vec2i& size) {
    if constexpr (N == 3) {
        // FIXME: not generic
        return make_ramp_image(img, size, vec<T, N>{0, 0, 0},
            vec<T, N>{1, 0, 0}, vec<T, N>{1, 1, 0}, vec<T, N>{0, 1, 0});
    } else if constexpr (N == 4) {
        // FIXME: not generic
        return make_ramp_image(img, size, vec<T, N>{0, 0, 0, 0},
            vec<T, N>{1, 0, 0, 0}, vec<T, N>{1, 1, 0, 0},
            vec<T, N>{0, 1, 0, 0});
    } else {
        throw runtime_error("bad channels");
    }
}

// Make a uv colored grid
template <typename T, int N>
inline void make_uvgrid_image(
    image<vec<T, N>>& img, const vec2i& size, int tiles, bool colored) {
    make_image_fromij(
        img, size, [size, tile = size.x / tiles, colored](int i, int j) {
            j       = size.y - j - 1;
            auto ii = i / tile, jj = j / tile;
            auto ww = size.x / tile, hh = size.y / tile;
            auto ph = (((256 / (ww * hh)) * (ii + jj * ww) - 64 + 256) % 256) /
                      360.f;
            auto pv = 0.5f;
            auto ps = 0.8f;
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 0.05f;
                else
                    pv -= 0.05f;
            } else {
                pv = 0.8f;
                ps = 0.2f;
            }
            auto rgb = (colored) ? hsv_to_rgb(vec<T, 3>{ph, ps, pv})
                                 : vec<T, 3>{pv, pv, pv};
            if constexpr (N == 3) {
                return vec<T, 3>{rgb.x, rgb.y, rgb.z};
            } else if constexpr (N == 4) {
                return vec<T, 4>{rgb.x, rgb.y, rgb.z, 1};
            } else {
                throw runtime_error("bad number of channels");
            }
        });
}

// Makes a blackbody ramp
template <typename T, int N>
inline void make_blackbodyramp_image(image<vec<T, N>>& img, const vec2i& size,
    float start_temperature, float end_temperature) {
    make_image_fromij(
        img, size, [size, start_temperature, end_temperature](int i, int j) {
            auto temperature = start_temperature +
                               (end_temperature - start_temperature) *
                                   (float)i / (float)(size.x - 1);
            auto rgb = blackbody_to_rgb(temperature);
            if constexpr (N == 3) {
                return vec<T, 3>{rgb.x, rgb.y, rgb.z};
            } else if constexpr (N == 4) {
                return vec<T, 4>{rgb.x, rgb.y, rgb.z, 1};
            } else {
                throw runtime_error("bad number of channels");
            }
        });
}

// Make a noise image. Wrap works only if size is a power of two.
template <typename T, typename T1>
inline void make_noise_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale, bool wrap) {
    make_image_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale, &c0,
            &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g      = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            return lerp_color(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
template <typename T, typename T1>
inline void make_fbm_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale, T1 lacunarity, T1 gain, int octaves, bool wrap) {
    make_image_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g      = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            return lerp_color(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
template <typename T, typename T1>
inline void make_ridge_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale, T1 lacunarity, T1 gain, T1 offset, int octaves,
    bool wrap) {
    make_image_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, offset, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            return lerp_color(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
template <typename T, typename T1>
inline void make_turbulence_image(image<T>& img, const vec2i& size, const T& c0,
    const T& c1, T1 scale, T1 lacunarity, T1 gain, int octaves, bool wrap) {
    make_image_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_turbulence_noise(
                p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            return lerp_color(c0, c1, g);
        });
}

// Comvert a bump map to a normal map.
template <typename T, int N>
inline void bump_to_normal_map(
    image<vec<T, N>>& norm, const image<vec<T, N>>& img, T scale) {
    norm.resize(img.size());
    auto dx = 1.0f / img.size().x, dy = 1.0f / img.size().y;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto i1 = (i + 1) % img.size().x, j1 = (j + 1) % img.size().y;
            auto p00 = img[{i, j}], p10 = img[{i1, j}], p01 = img[{i, j1}];
            auto g00    = (p00.x + p00.y + p00.z) / 3;
            auto g01    = (p01.x + p01.y + p01.z) / 3;
            auto g10    = (p10.x + p10.y + p10.z) / 3;
            auto normal = vec3f{
                scale * (g00 - g10) / dx, scale * (g00 - g01) / dy, 1.0f};
            normal.y = -normal.y;  // make green pointing up, even if y axis
                                   // points down
            normal       = normalize(normal) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            norm[{i, j}] = {normal.x, normal.y, normal.z, 1};
        }
    }
}

// Add a border to an image
template <typename T>
inline void add_image_border(
    image<T>& img, int border_width, const T& border_color) {
    for (auto j = 0; j < img.size().y; j++) {
        for (auto b = 0; b < border_width; b++) {
            img[{b, j}]                    = border_color;
            img[{img.size().x - 1 - b, j}] = border_color;
        }
    }
    for (auto i = 0; i < img.size().x; i++) {
        for (auto b = 0; b < border_width; b++) {
            img[{i, b}]                    = border_color;
            img[{i, img.size().y - 1 - b}] = border_color;
        }
    }
}

#if 1

// Implementation of sunsky modified heavily from pbrt
template <typename T, int N>
inline void make_sunsky_image(image<vec<T, N>>& img, const vec2i& size,
    T theta_sun, T turbidity, bool has_sun, T sun_intensity, T sun_temperature,
    const vec<T, 3>& ground_albedo) {
    // idea adapted from pbrt

    // initialize model
    double wavelengths[9] = {630, 680, 710, 500, 530, 560, 460, 480, 490};
    ArHosekSkyModelState* skymodel_state[9];
    if (sun_temperature) {
        sun_temperature = clamp(sun_temperature, (T)2000, (T)14000);
        for (int i = 0; i < 9; ++i) {
            skymodel_state[i] = arhosekskymodelstate_alienworld_alloc_init(
                theta_sun, sun_intensity, sun_temperature, turbidity,
                ground_albedo[i / 3]);
        }
    } else {
        for (int i = 0; i < 9; ++i) {
            skymodel_state[i] = arhosekskymodelstate_alloc_init(
                theta_sun, turbidity, ground_albedo[i / 3]);
        }
    }

    // clear image
    img.resize(size);
    if constexpr (N == 3) {
        for (auto& p : img) p = {0, 0, 0};
    } else if constexpr (N == 4) {
        for (auto& p : img) p = {0, 0, 0, default_alpha<T>};
    } else {
        throw runtime_error("bad channels");
    }

    // sun-sky
    auto sun_direction = vec<T, 3>{0, sin(theta_sun), cos(theta_sun)};
    auto integral      = zero3f;
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = (j + (T)0.5) * (T)pi / img.size().y;
        if (theta > pif / 2) continue;
        for (auto i = 0; i < img.size().x; i++) {
            auto phi       = (i + (T)0.5) * 2 * (T)pi / img.size().x;
            auto direction = vec3f{
                cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(
                clamp(dot(direction, sun_direction), (T)-1, (T)1));
            for (int c = 0; c < 9; ++c) {
                auto val =
                    (has_sun)
                        ? arhosekskymodel_solar_radiance(
                              skymodel_state[c], theta, gamma, wavelengths[c])
                        : arhosekskymodel_radiance(
                              skymodel_state[c], theta, gamma, wavelengths[c]);
                // average channel over wavelengths
                img[{i, j}][c / 3] += (float)val / 3;
            }
            integral += img[{i, j}].xyz * sin(theta) /
                        (img.size().x * img.size().y / 2);
        }
    }

    // ground
    auto ground = ground_albedo * integral;
    for (auto j = img.size().y / 2; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            img[{i, j}] = {ground.x, ground.y, ground.z, 1};
        }
    }

    // cleanup
    for (auto i = 0; i < 9; i++) arhosekskymodelstate_free(skymodel_state[i]);
}

#else

// Implementation of sunsky modified heavily from pbrt
image<vec4f> make_sunsky_image(int width, int height, float theta_sun,
    float turbidity, bool has_sun, float sun_angle_scale,
    float sun_emission_scale, const vec3f& ground_albedo,
    bool renormalize_sun) {
    auto zenith_xyY = vec3f{
        (+0.00165f * pow(theta_sun, 3.f) - 0.00374f * pow(theta_sun, 2.f) +
            0.00208f * theta_sun + 0) *
                pow(turbidity, 2.f) +
            (-0.02902f * pow(theta_sun, 3.f) + 0.06377f * pow(theta_sun, 2.f) -
                0.03202f * theta_sun + 0.00394f) *
                turbidity +
            (+0.11693f * pow(theta_sun, 3.f) - 0.21196f * pow(theta_sun, 2.f) +
                0.06052f * theta_sun + 0.25885f),
        (+0.00275f * pow(theta_sun, 3.f) - 0.00610f * pow(theta_sun, 2.f) +
            0.00316f * theta_sun + 0) *
                pow(turbidity, 2.f) +
            (-0.04214f * pow(theta_sun, 3.f) + 0.08970f * pow(theta_sun, 2.f) -
                0.04153f * theta_sun + 0.00515f) *
                turbidity +
            (+0.15346f * pow(theta_sun, 3.f) - 0.26756f * pow(theta_sun, 2.f) +
                0.06669f * theta_sun + 0.26688f),
        1000 * (4.0453f * turbidity - 4.9710f) *
                tan((4.0f / 9.0f - turbidity / 120.0f) *
                    (pif - 2 * theta_sun)) -
            .2155f * turbidity + 2.4192f};

    auto perez_A_xyY = vec3f{-0.01925f * turbidity - 0.25922f,
        -0.01669f * turbidity - 0.26078f, +0.17872f * turbidity - 1.46303f};
    auto perez_B_xyY = vec3f{-0.06651f * turbidity + 0.00081f,
        -0.09495f * turbidity + 0.00921f, -0.35540f * turbidity + 0.42749f};
    auto perez_C_xyY = vec3f{-0.00041f * turbidity + 0.21247f,
        -0.00792f * turbidity + 0.21023f, -0.02266f * turbidity + 5.32505f};
    auto perez_D_xyY = vec3f{-0.06409f * turbidity - 0.89887f,
        -0.04405f * turbidity - 1.65369f, +0.12064f * turbidity - 2.57705f};
    auto perez_E_xyY = vec3f{-0.00325f * turbidity + 0.04517f,
        -0.01092f * turbidity + 0.05291f, -0.06696f * turbidity + 0.37027f};

    auto perez_f = [](vec3f A, vec3f B, vec3f C, vec3f D, vec3f E, float theta,
                       float gamma, float theta_sun, vec3f zenith) -> vec3f {
        auto den =
            ((1 + A * exp(B)) * (1 + C * exp(D * theta_sun) +
                                    E * cos(theta_sun) * cos(theta_sun)));
        auto num = ((1 + A * exp(B / cos(theta))) *
                    (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
        return zenith * num / den;
    };

    auto sky = [&perez_f, perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                   perez_E_xyY, zenith_xyY](
                   float theta, float gamma, float theta_sun) -> vec3f {
        return xyz_to_rgb(xyY_to_xyz(
                   perez_f(perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                       perez_E_xyY, theta, gamma, theta_sun, zenith_xyY))) /
               10000;
    };

    // compute sun luminance
    // TODO: how this relates to zenith intensity?
    auto sun_ko     = vec3f{0.48f, 0.75f, 0.14f};
    auto sun_kg     = vec3f{0.1f, 0.0f, 0.0f};
    auto sun_kwa    = vec3f{0.02f, 0.0f, 0.0f};
    auto sun_sol    = vec3f{20000.0f, 27000.0f, 30000.0f};
    auto sun_lambda = vec3f{680, 530, 480};
    auto sun_beta   = 0.04608365822050f * turbidity - 0.04586025928522f;
    auto sun_m =
        1.0f / (cos(theta_sun) + 0.000940f * pow(1.6386f - theta_sun, -1.253f));

    auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda / 1000, -4.08f));
    auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda / 1000, -1.3f));
    auto tauO = exp(-sun_m * sun_ko * .35f);
    auto tauG = exp(
        -1.41f * sun_kg * sun_m / pow(1 + 118.93f * sun_kg * sun_m, 0.45f));
    auto tauWA  = exp(-0.2385f * sun_kwa * 2.0f * sun_m /
                     pow(1 + 20.07f * sun_kwa * 2.0f * sun_m, 0.45f));
    auto sun_le = sun_sol * tauR * tauA * tauO * tauG * tauWA;

    // rescale by user
    sun_le *= sun_emission_scale;

    // sun scale from Wikipedia scaled by user quantity and rescaled to at
    // the minimum 5 pixel diamater
    auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
    sun_angular_radius *= sun_angle_scale;
    sun_angular_radius = max(sun_angular_radius, 5 * pif / height);

    // sun direction
    auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

    auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
        // return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f :
        //                                                zero3f;
        return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000
                                                       : zero3f;
    };

    // Make the sun sky image
    auto img          = make_image(width, height, vec4f{0, 0, 0, 1});
    auto sky_integral = 0.0f, sun_integral = 0.0f;
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = pif * ((j + 0.5f) / img.size().y);
        theta      = clamp(theta, 0.0f, pif / 2 - float_epsilon);
        for (int i = 0; i < img.size().x; i++) {
            auto phi = 2 * pif * (float(i + 0.5f) / img.size().x);
            auto w   = vec3f{
                cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
            auto sky_col = sky(theta, gamma, theta_sun);
            auto sun_col = sun(theta, gamma);
            sky_integral += mean(sky_col) * sin(theta);
            sun_integral += mean(sun_col) * sin(theta);
            auto col    = sky_col + sun_col;
            img[{i, j}] = {col.x, col.y, col.z, 1};
        }
    }

    if (renormalize_sun) {
        for (auto j = 0; j < img.size().y / 2; j++) {
            for (int i = 0; i < img.size().x; i++) {
                img[{i, j}] *= sky_integral / (sun_integral + sky_integral);
            }
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < img.size().y / 2; j++) {
            auto theta = pif * ((j + 0.5f) / img.size().y);
            for (int i = 0; i < img.size().x; i++) {
                auto pxl   = img[{i, j}];
                auto le    = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pif /
                             (img.size().x * img.size().y);
                ground += le * (ground_albedo / pif) * cos(theta) * angle;
            }
        }
        for (auto j = img.size().y / 2; j < img.size().y; j++) {
            for (int i = 0; i < img.size().x; i++) {
                img[{i, j}] = {ground.x, ground.y, ground.z, 1};
            }
        }
    }
    return img;
}

#endif

// Make an image of multiple lights.
template <typename T>
inline void make_lights_image(image<vec<T, 4>>& img, const vec2i& size,
    const vec<T, 3>& le, int nlights, T langle, T lwidth, T lheight) {
    img.resize(size);
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = (T)pi * ((j + (T)0.5) / img.size().y);
        theta      = clamp(theta, (T)0, (T)pi / 2 - float_epsilon);
        if (fabs(theta - langle) > lheight / 2) continue;
        for (int i = 0; i < img.size().x; i++) {
            auto phi     = 2 * (T)pi * (float(i + (T)0.5) / img.size().x);
            auto inlight = false;
            for (auto l = 0; l < nlights; l++) {
                auto lphi = 2 * (T)pi * (l + (T)0.5) / nlights;
                inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
            }
            img[{i, j}] = rgba(le);
        }
    }
}

template <int N>
inline void make_logo_image(image<vec<byte, N>>& img, const string& type) {
    static const auto size = vec2i{144, 28};
    // clang-format off
    static const auto logo_render = vector<byte>{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 212, 87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 255, 62, 0, 0, 0, 0, 0, 14, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 187, 245, 13, 0, 0, 0, 80, 255, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29, 88, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 88, 251, 10, 0, 0, 0, 40, 200, 253, 255, 234, 106, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90, 69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 144, 125, 0, 0, 0, 0, 0, 117, 37, 0, 0, 0, 0, 0, 0, 0, 79, 255, 101, 0, 0, 0, 178, 232, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 145, 205, 0, 0, 0, 47, 239, 210, 74, 57, 144, 232, 24, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 251, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 146, 123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 43, 35, 0, 61, 87, 0, 0, 208, 61, 0, 0, 0, 0, 0, 0, 0, 3, 224, 199, 0, 0, 24, 251, 129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 201, 149, 0, 0, 0, 180, 238, 20, 0, 0, 0, 16, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 251, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 146, 123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 119, 150, 0, 0, 208, 61, 0, 0, 0, 0, 0, 0, 0, 0, 118, 255, 41, 0, 117, 250, 26, 0, 7, 147, 223, 232, 167, 19, 0, 0, 0, 5, 143, 224, 234, 162, 20, 166, 227, 255, 210, 201, 3, 0, 7, 147, 223, 232, 167, 19, 0, 0, 0, 0, 8, 249, 93, 0, 0, 45, 255, 146, 0, 0, 0, 0, 0, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 14, 192, 114, 207, 0, 84, 224, 220, 61, 0, 68, 143, 166, 220, 78, 0, 0, 142, 234, 160, 251, 0, 0, 134, 234, 195, 25, 0, 123, 105, 211, 89, 12, 177, 236, 158, 4, 0, 44, 217, 214, 189, 123, 0, 0, 0, 149, 76, 0, 189, 108, 0, 160, 55, 123, 107, 83, 236, 240, 179, 0, 208, 128, 220, 174, 6, 0, 0, 0, 0, 0, 19, 247, 140, 0, 214, 168, 0, 0, 176, 248, 122, 111, 237, 212, 2, 0, 0, 167, 251, 128, 127, 223, 59, 103, 185, 255, 143, 117, 0, 0, 176, 248, 122, 111, 237, 212, 2, 0, 0, 0, 58, 255, 36, 0, 0, 97, 255, 77, 0, 0, 0, 0, 0, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 248, 166, 36, 13, 234, 44, 73, 215, 0, 80, 248, 65, 80, 193, 0, 66, 222, 32, 97, 251, 0, 65, 207, 21, 135, 152, 0, 144, 230, 81, 12, 129, 157, 17, 189, 88, 0, 194, 126, 17, 208, 123, 0, 0, 0, 131, 125, 7, 228, 164, 0, 224, 23, 144, 125, 5, 127, 156, 10, 0, 208, 179, 13, 205, 65, 0, 0, 0, 0, 0, 0, 158, 233, 61, 255, 59, 0, 46, 255, 121, 0, 0, 91, 255, 83, 0, 40, 254, 134, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 46, 255, 121, 0, 0, 91, 255, 83, 0, 0, 0, 114, 235, 0, 0, 0, 130, 255, 51, 0, 2, 17, 17, 17, 7, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 255, 43, 0, 76, 191, 0, 1, 242, 20, 80, 192, 0, 36, 235, 0, 138, 142, 0, 18, 251, 0, 140, 127, 0, 52, 212, 0, 144, 171, 0, 0, 204, 63, 0, 116, 148, 11, 255, 14, 0, 146, 123, 0, 0, 0, 84, 164, 46, 155, 201, 9, 231, 0, 144, 125, 0, 119, 150, 0, 0, 208, 64, 0, 164, 107, 0, 0, 0, 0, 0, 0, 49, 255, 220, 206, 0, 0, 114, 255, 46, 0, 0, 17, 255, 148, 0, 111, 255, 54, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 114, 255, 46, 0, 0, 17, 255, 148, 0, 0, 0, 171, 179, 0, 0, 0, 159, 255, 29, 0, 20, 255, 255, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 101, 242, 214, 214, 249, 40, 80, 189, 0, 34, 236, 0, 162, 121, 0, 18, 251, 0, 165, 232, 214, 219, 232, 0, 144, 126, 0, 0, 229, 222, 214, 229, 168, 34, 249, 0, 0, 146, 123, 0, 0, 0, 38, 203, 88, 107, 206, 48, 187, 0, 144, 125, 0, 119, 150, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 197, 255, 98, 0, 0, 144, 255, 22, 0, 0, 0, 249, 177, 0, 141, 255, 28, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 144, 255, 22, 0, 0, 0, 249, 177, 0, 0, 0, 227, 123, 0, 0, 0, 143, 255, 40, 0, 0, 79, 108, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 85, 188, 4, 4, 4, 0, 80, 189, 0, 34, 236, 0, 149, 132, 0, 18, 251, 0, 149, 125, 4, 4, 3, 0, 144, 125, 0, 0, 213, 62, 4, 4, 2, 21, 255, 4, 0, 146, 123, 0, 0, 0, 2, 230, 130, 67, 178, 116, 141, 0, 144, 125, 0, 119, 150, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 151, 255, 17, 0, 0, 0, 245, 183, 0, 152, 255, 21, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 151, 255, 17, 0, 0, 0, 245, 183, 0, 0, 28, 255, 67, 0, 0, 0, 116, 255, 59, 0, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 32, 236, 19, 2, 33, 0, 80, 189, 0, 34, 236, 0, 98, 195, 4, 64, 251, 0, 93, 191, 4, 11, 24, 0, 144, 125, 0, 0, 157, 131, 0, 27, 8, 1, 225, 72, 1, 190, 123, 0, 0, 0, 0, 201, 196, 26, 137, 197, 96, 0, 144, 125, 0, 109, 159, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 122, 255, 36, 0, 0, 8, 255, 152, 0, 124, 255, 42, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 122, 255, 36, 0, 0, 8, 255, 152, 0, 0, 84, 253, 13, 0, 0, 0, 79, 255, 108, 0, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 0, 123, 238, 224, 155, 0, 80, 189, 0, 34, 236, 0, 11, 201, 218, 156, 248, 0, 6, 178, 225, 234, 97, 0, 144, 125, 0, 0, 31, 216, 214, 228, 49, 0, 89, 244, 198, 179, 123, 0, 0, 0, 0, 154, 239, 0, 96, 253, 50, 0, 144, 125, 0, 34, 224, 201, 25, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 70, 255, 96, 0, 0, 67, 255, 97, 0, 76, 255, 104, 0, 0, 0, 0, 0, 113, 255, 35, 0, 0, 70, 255, 96, 0, 0, 67, 255, 97, 0, 0, 141, 210, 0, 0, 0, 0, 5, 230, 200, 2, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 19, 6, 0, 0, 0, 0, 0, 0, 0, 0, 23, 1, 0, 0, 0, 12, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 1, 211, 223, 54, 44, 205, 229, 7, 0, 2, 216, 230, 69, 74, 169, 31, 0, 55, 255, 120, 59, 26, 1, 211, 223, 54, 44, 205, 229, 7, 0, 0, 197, 154, 0, 0, 0, 0, 0, 109, 255, 166, 59, 78, 165, 255, 109, 0, 6, 255, 201, 113, 113, 113, 105, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 0, 30, 204, 255, 255, 214, 40, 0, 0, 0, 35, 208, 255, 255, 212, 47, 0, 1, 174, 252, 243, 102, 0, 30, 204, 255, 255, 214, 40, 0, 0, 6, 247, 97, 0, 0, 0, 0, 0, 0, 103, 235, 255, 255, 242, 146, 24, 0, 6, 255, 255, 255, 255, 255, 213, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 29, 0, 0, 0, 0, 0, 0, 0, 25, 31, 0, 0, 0, 0, 0, 23, 9, 0, 0, 0, 0, 24, 29, 0, 0, 0, 0, 54, 255, 41, 0, 0, 0, 0, 0, 0, 0, 0, 33, 30, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 63, 188, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    };
    // clang-format on
    if (type == "logo-render") {
        auto img1 = image<byte>{size, logo_render.data()};
        img.resize(size);
        gray_to_color(img, img1);
    } else {
        throw io_error("unknown builtin image " + type);
    }
}

template <int N>
inline void make_logo_image(image<vec<float, N>>& img, const string& type) {
    auto img8 = image<vec<byte, N>>();
    make_logo_image(img8, type);
    img.resize(img8.size());
    srgb8_to_linear(img, img8);
}

inline void make_image_preset(image<vec<float, 4>>& img, const string& type) {
    auto size = vec2i{1024, 1024};
    if (type.find("sky") != type.npos) size = {2048, 1024};
    if (type == "grid") {
        make_grid_image(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "checker") {
        make_checker_image(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "bump") {
        make_bumpdimple_image(img, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "uvramp") {
        make_uvramp_image(img, size);
    } else if (type == "gammaramp") {
        make_gammaramp_image(img, size, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "blackbodyramp") {
        make_blackbodyramp_image(img, size);
    } else if (type == "uvgrid") {
        make_uvgrid_image(img, size);
    } else if (type == "sky") {
        make_sunsky_image(img, size, pif / 4, 3.0f, false, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "sunsky") {
        make_sunsky_image(img, size, pif / 4, 3.0f, true, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "noise") {
        make_noise_image(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, true);
    } else if (type == "fbm") {
        make_fbm_image(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "ridge") {
        make_ridge_image(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f,
            0.5f, 1.0f, 6, true);
    } else if (type == "turbulence") {
        make_turbulence_image(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "bump-normal") {
        auto bump = image<vec<float, 4>>{};
        make_bumpdimple_image(bump, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
        bump_to_normal_map(img, bump, 0.05f);
    } else if (type == "logo-render") {
        make_logo_image(img, "logo-render");
    } else if (type == "images1") {
        auto sub_types = vector<string>{"grid", "uvgrid", "checker",
            "gammaramp", "bump", "bump-normal", "noise", "fbm",
            "blackbodyramp"};
        auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
        for (auto i = 0; i < sub_imgs.size(); i++) {
            sub_imgs.at(i).resize(img.size());
            make_image_preset(sub_imgs.at(i), sub_types.at(i));
        }
        auto montage_size = zero2i;
        for (auto& sub_img : sub_imgs) {
            montage_size.x += sub_img.size().x;
            montage_size.y = max(montage_size.y, sub_img.size().y);
        }
        img.resize(montage_size);
        auto pos = 0;
        for (auto& sub_img : sub_imgs) {
            set_image_region(img, sub_img, {pos, 0});
            pos += sub_img.size().x;
        }
    } else if (type == "images2") {
        auto sub_types = vector<string>{"sky", "sunsky"};
        auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
        for (auto i = 0; i < sub_imgs.size(); i++) {
            make_image_preset(sub_imgs.at(i), sub_types.at(i));
        }
        auto montage_size = zero2i;
        for (auto& sub_img : sub_imgs) {
            montage_size.x += sub_img.size().x;
            montage_size.y = max(montage_size.y, sub_img.size().y);
        }
        img.resize(montage_size);
        auto pos = 0;
        for (auto& sub_img : sub_imgs) {
            set_image_region(img, sub_img, {pos, 0});
            pos += sub_img.size().x;
        }
    } else if (type == "test-floor") {
        make_grid_image(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
        add_image_border(img, 2, {0, 0, 0, 1});
    } else if (type == "test-grid") {
        make_grid_image(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "test-checker") {
        make_checker_image(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "test-bump") {
        make_bumpdimple_image(img, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "test-uvramp") {
        make_uvramp_image(img, size);
    } else if (type == "test-gammaramp") {
        make_gammaramp_image(img, size, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "test-blackbodyramp") {
        make_blackbodyramp_image(img, size);
    } else if (type == "test-uvgrid") {
        make_uvgrid_image(img, size);
    } else if (type == "test-sky") {
        make_sunsky_image(img, size, pif / 4, 3.0f, false, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "test-sunsky") {
        make_sunsky_image(img, size, pif / 4, 3.0f, true, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "test-noise") {
        make_noise_image(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, true);
    } else if (type == "test-fbm") {
        make_fbm_image(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "test-bump-normal") {
        auto bump = image<vec<float, 4>>{};
        make_bumpdimple_image(bump, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
        bump_to_normal_map(img, bump, 0.05f);
    } else if (type == "test-fbm-displacement") {
        make_fbm_image(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 10.0f, 2.0f, 0.5f, 6, true);
    } else {
        throw std::invalid_argument("unknown image preset" + type);
    }
}

inline void make_image_preset(image<vec<byte, 4>>& img, const string& type) {
    auto imgf = image<vec4f>{};
    make_image_preset(imgf, type);
    if (type.find("-normal") == type.npos) {
        linear_to_srgb8(img, imgf);
    } else {
        float_to_byte(img, imgf);
    }
}

inline void make_image_preset(
    image<vec4f>& hdr, image<vec4b>& ldr, const string& type) {
    if (type.find("sky") == type.npos) {
        make_image_preset(ldr, type);
    } else {
        make_image_preset(hdr, type);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// make a simple example volume
inline void make_test_volume(
    volume<float>& vol, const vec3i& size, float scale, float exponent) {
    vol.resize(size);
    for (auto k = 0; k < vol.size().z; k++) {
        for (auto j = 0; j < vol.size().y; j++) {
            for (auto i = 0; i < vol.size().x; i++) {
                auto p = vec3f{i / (float)vol.size().x, j / (float)vol.size().y,
                    k / (float)vol.size().z};
                auto value = pow(
                    max(max(cos(scale * p.x), cos(scale * p.y)), 0.0f),
                    exponent);
                vol[{i, j, k}] = clamp(value, 0.0f, 1.0f);
            }
        }
    }
}

inline void make_volume_preset(volume<float>& vol, const string& type) {
    auto size = vec3i{256, 256, 256};
    if (type == "test-volume") {
        make_test_volume(vol, size, 6, 10);
    } else {
        throw runtime_error("unknown volume preset " + type);
    }
}

}  // namespace yocto

#endif
