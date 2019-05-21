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
// 3. resize images with `resize()`
// 4. tonemap images with `tonemap()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky()`
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
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Float to byte conversion.
inline byte  float_to_byte(float a);
inline float byte_to_float(byte a);
inline vec1b float_to_byte(const vec1f& a);
inline vec1f byte_to_float(const vec1b& a);
inline vec2b float_to_byte(const vec2f& a);
inline vec2f byte_to_float(const vec2b& a);
inline vec3b float_to_byte(const vec3f& a);
inline vec3f byte_to_float(const vec3b& a);
inline vec4b float_to_byte(const vec4f& a);
inline vec4f byte_to_float(const vec4b& a);

// Color luminance
inline float luminance(const vec3f& a);

// sRGB non-linear curve
inline float srgb_to_linear(float srgb);
inline float linear_to_srgb(float lin);
inline vec1f srgb_to_linear(const vec1f& srgb);
inline vec1f linear_to_srgb(const vec1f& lin);
inline vec2f srgb_to_linear(const vec2f& srgb);
inline vec2f linear_to_srgb(const vec2f& lin);
inline vec3f srgb_to_linear(const vec3f& srgb);
inline vec3f linear_to_srgb(const vec3f& lin);
inline vec4f srgb_to_linear(const vec4f& srgb);
inline vec4f linear_to_srgb(const vec4f& lin);

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f apply_contrast(const vec3f& rgb, float contrast, float grey);
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f apply_logcontrast(const vec3f& rgb, float logcontrast, float grey);
// Apply saturation.
inline vec3f apply_saturation(const vec3f& rgb, float saturation,
    const vec3f& luminance_weights = vec3f{0.25, 0.5, 0.25});

// Fitted ACES tonemapping curve.
inline vec3f tonemap_filmic(const vec3f& hdr, bool accurate_fit = true);

// Convert number of channels
template <typename T, int N1, int N2>
inline vec<T, N1> convert_channels(const vec<T, N2>& a);

// Forward declaration
struct rgb_space;

// Conversion between rgb color spaces
inline vec3f rgb_to_rgb(
    const vec3f& rgb, const rgb_space& from_space, const rgb_space& to_space);

// Conversion to/from xyz
inline vec3f rgb_to_xyz(const vec3f& rgb, const rgb_space& rgb_space);
inline vec3f xyz_to_rgb(const vec3f& xyz, const rgb_space& rgb_space);

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz);
// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY);

// Approximate color of blackbody radiation from wavelength in nm.
inline vec3f blackbody_to_rgb(float temperature);

// Converts between HSV and RGB color spaces.
inline vec3f hsv_to_rgb(const vec3f& hsv);
inline vec3f rgb_to_hsv(const vec3f& rgb);

// Tone curves
// Pure gamma tone curve: y = x^gamma
constexpr float apply_gamma(float x, float gamma);
constexpr float apply_gamma_inv(float x, float gamma);
// Pure gamma tone curve: y = (x < d) ? x * c : pow(x * a + b, gamma)
constexpr float apply_gamma(float x, float gamma, const vec4f& abcd);
constexpr float apply_gamma_inv(float x, float gamma, const vec4f& abcd);

// RGB color spaces
enum struct rgb_space_type {
    linear_srgb,   // default linear space (srgb linear)
    srgb,          // srgb color space (non-linear)
    adobe_rgb,     // Adobe rgb color space (non-linear)
    prophoto_rgb,  // ProPhoto Kodak rgb color space (non-linear)
    rec_709,       // hdtv color space (non-linear)
    rec_2020,      // uhtv color space (non-linear)
    rec_2100_pq,   // hdr color space with perceptual quantizer (non-linear)
    rec_2100_hlg,  // hdr color space with hybrid log gamma (non-linear)
    aces_2065,     // ACES storage format (linear)
    aces_cg,       // ACES CG computation (linear)
    aces_cc,       // ACES color correction (non-linear)
    aces_cct,      // ACES color correction 2 (non-linear)
    p3_dci,        // P3 DCI (non-linear)
    p3_d60,        // P3 variation for D60 (non-linear)
    p3_d65,        // P3 variation for D65 (non-linear)
    p3_display,    // Apple display P3
};

// get the color space definition for builtin color spaces
constexpr rgb_space get_rgb_space(rgb_space_type type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container.
template <typename T>
struct image {
    // constructors
    image();
    image(const vec2i& size, const T& value = {});
    image(const vec2i& size, const T* value);

    // size
    bool  empty() const;
    vec2i size() const;
    bool  contains(const vec2i& ij) const;
    void  resize(const vec2i& size);
    void  assign(const vec2i& size, const T& value = {});
    void  shrink_to_fit();

    // element access
    T&       operator[](const vec2i& ij);
    const T& operator[](const vec2i& ij) const;

    // data access
    T*       data();
    const T* data() const;

    // iteration
    T*       begin();
    T*       end();
    const T* begin() const;
    const T* end() const;

   private:
    // data
    vec2i     extent = zero2i;
    vector<T> pixels = {};
};

// equality
template <typename T>
inline bool operator==(const image<T>& a, const image<T>& b);
template <typename T>
inline bool operator!=(const image<T>& a, const image<T>& b);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear color space.
void load_image(const string& filename, image<vec4f>& img);
void save_image(const string& filename, const image<vec4f>& img);
void load_image(const string& filename, image<vec4b>& img);
void save_image(const string& filename, const image<vec4b>& img);

// Convenience helper for loading HDR or LDR based on filename
void load_image(
    const string& filename, image<vec4f>& hdr, image<vec4b>& ldr);
void save_image(
    const string& filename, const image<vec4f>& hdr, const image<vec4b>& ldr);

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
inline void make_regions(vector<image_region>& regions, const vec2i& size,
    int region_size = 32, bool shuffled = false);

// Gets pixels in an image region
template <typename T>
inline void get_region(
    image<T>& clipped, const image<T>& img, const image_region& region);
template <typename T>
inline void set_region(
    image<T>& img, const image<T>& region, const vec2i& offset);

// Conversion from/to floats.
template <typename T, typename TB>
inline void byte_to_float(image<T>& fl, const image<TB>& bt);
template <typename T, typename TB>
inline void float_to_byte(image<TB>& bt, const image<T>& fl);

// Conversion between linear and gamma-encoded images.
template <typename T>
inline void srgb_to_linear(image<T>& lin, const image<T>& srgb);
template <typename T>
inline void linear_to_srgb(image<T>& srgb, const image<T>& lin);

// Tone mapping params
struct tonemap_params {
    float exposure    = 0;
    vec3f tint        = {1, 1, 1};
    float contrast    = 0.5f;
    float logcontrast = 0.5f;
    float saturation  = 0.5f;
    bool  filmic      = false;
    bool  srgb        = true;
};

// Equality operators
inline bool operator==(const tonemap_params& a, const tonemap_params& b) {
    return memcmp(&a, &b, sizeof(a)) == 0;
}
inline bool operator!=(const tonemap_params& a, const tonemap_params& b) {
    return memcmp(&a, &b, sizeof(a)) != 0;
}

// Apply exposure and filmic tone mapping
inline void tonemap(
    image<vec4f>& ldr, const image<vec4f>& hdr, const tonemap_params& params);
inline void tonemap(
    image<vec4b>& ldr, const image<vec4f>& hdr, const tonemap_params& params);
inline void tonemap(image<vec4f>& ldr, const image<vec4f>& hdr,
    const image_region& region, const tonemap_params& params);

// minimal color grading
struct colorgrade_params {
    float contrast         = 0.5;
    float shadows          = 0.5;
    float midtones         = 0.5;
    float highlights       = 0.5;
    vec3f shadows_color    = {1, 1, 1};
    vec3f midtones_color   = {1, 1, 1};
    vec3f highlights_color = {1, 1, 1};
};

// Equality operators
inline bool operator==(const colorgrade_params& a, const colorgrade_params& b) {
    return memcmp(&a, &b, sizeof(a)) == 0;
}
inline bool operator!=(const colorgrade_params& a, const colorgrade_params& b) {
    return memcmp(&a, &b, sizeof(a)) != 0;
}

// color grade an image region
inline void colorgrade(image<vec4f>& corrected, const image<vec4f>& img,
    const image_region& region, const colorgrade_params& params);

// determine white balance colors
inline vec3f compute_white_balance(const image<vec4f>& img);

// Convert number of channels
template <int N1, int N2>
inline void convert_channels(
    image<vec<float, N1>>& result, const image<vec<float, N2>>& source);
template <int N1, int N2>
inline void convert_channels(
    image<vec<byte, N1>>& result, const image<vec<byte, N2>>& source);

// Resize an image.
inline void resize(image<float>& res, const image<float>& img);
inline void resize(image<vec1f>& res, const image<vec1f>& img);
inline void resize(image<vec2f>& res, const image<vec2f>& img);
inline void resize(image<vec3f>& res, const image<vec3f>& img);
inline void resize(image<vec4f>& res, const image<vec4f>& img);
inline void resize(image<byte>& res, const image<byte>& img);
inline void resize(image<vec1b>& res, const image<vec1b>& img);
inline void resize(image<vec2b>& res, const image<vec2b>& img);
inline void resize(image<vec3b>& res, const image<vec3b>& img);
inline void resize(image<vec4b>& res, const image<vec4b>& img);
template <typename T>
inline void resize(image<T>& res, const image<T>& img, const vec2i& size);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make example images in linear color space. Takes as input images allocated
// to the desired size and fill the pixel with expected values.
void make_grid(image<vec4f>& img, const vec2i& size, int tile, const vec4f& c0,
    const vec4f& c1);
void make_checker(image<vec4f>& img, const vec2i& size, int tile,
    const vec4f& c0, const vec4f& c1);
void make_bumpdimple(image<vec4f>& img, const vec2i& size, int tile,
    const vec4f& c0, const vec4f& c1);
void make_ramp(
    image<vec4f>& img, const vec2i& size, const vec4f& c0, const vec4f& c1);
void make_ramp(image<vec4f>& img, const vec2i& size, const vec4f& c00,
    const vec4f& c10, const vec4f& c11, const vec4f& c01);
void make_gammaramp(
    image<vec4f>& img, const vec2i& size, const vec4f& c0, const vec4f& c1);
void make_uvramp(image<vec4f>& img, const vec2i& size);
void make_uvgrid(
    image<vec4f>& img, const vec2i& size, int tile = 8, bool colored = true);
void make_blackbodyramp(image<vec4f>& img, const vec2i& size,
    float start_temperature = 1000, float end_temperature = 12000);

// Comvert a bump map to a normal map. All linear color spaces.
void bump_to_normal(
    image<vec4f>& norm, const image<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
void make_sunsky(image<vec4f>& img, const vec2i& size, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
    float sun_temperature = 0, const vec3f& ground_albedo = {0.2, 0.2, 0.2});
// Make an image of multiple lights.
void make_lights(image<vec4f>& img, const vec2i& size,
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
    float lwidth = pif / 16, float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
void make_noise(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale = 1, bool wrap = true);
void make_fbm(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale = 1, float lacunarity = 2, float gain = 0.5,
    int octaves = 6, bool wrap = true);
void make_ridge(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale = 1, float lacunarity = 2, float gain = 0.5,
    float offset = 1, int octaves = 6, bool wrap = true);
void make_turbulence(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale = 1, float lacunarity = 2, float gain = 0.5,
    int octaves = 6, bool wrap = true);

// Add a border to an image
void add_border(image<vec4f>& img, const vec2i& size, int border_width,
    const vec4f& border_color);

// Make logo images. Image is resized to proper size.
void make_logo(image<vec4f>& img, const string& name);
void make_logo(image<vec4b>& img, const string& name);

// Make an image preset, useful for testing. See implementation for types.
void make_preset(image<vec4f>& img, const string& type);
void make_preset(image<vec4b>& img, const string& type);
void make_preset(image<vec4f>& hdr, image<vec4b>& ldr, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SAVING HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
void save_tonemapped(const string& filename, const image<vec4f>& hdr,
    const tonemap_params& params);

// Save with a logo embedded
void save_image_with_logo(
    const string& filename, const image<vec4f>& img);
void save_image_with_logo(
    const string& filename, const image<vec4b>& img);

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
void save_tonemapped_with_logo(const string& filename,
    const image<vec4f>& hdr, const tonemap_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container.
template <typename T>
struct volume {
    // constructors
    volume();
    volume(const vec3i& size, const T& value = {});
    volume(const vec3i& size, const T* value);

    // size
    bool  empty() const;
    vec3i size() const;
    void  resize(const vec3i& size);
    void  assign(const vec3i& size, const T& value = {});
    void  shrink_to_fit();

    // element access
    T&       operator[](const vec3i& ijk);
    const T& operator[](const vec3i& ijk) const;

    // data access
    T*       data();
    const T* data() const;

    // iteration
    T*       begin();
    T*       end();
    const T* begin() const;
    const T* end() const;

   private:
    // data
    vec3i         extent = zero3i;
    vector<float> voxels = {};
};

// equality
template <typename T>
inline bool operator==(const volume<T>& a, const volume<T>& b);
template <typename T>
inline bool operator!=(const volume<T>& a, const volume<T>& b);

// make a simple example volume
void make_test(volume<float>& vol, const vec3i& size, float scale = 10,
    float exponent = 6);
void make_preset(volume<float>& vol, const string& type);

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// constructors
template <typename T>
inline image<T>::image() : extent{0, 0}, pixels{} {}
template <typename T>
inline image<T>::image(const vec2i& size, const T& value)
    : extent{size}, pixels((size_t)size.x * (size_t)size.y, value) {}
template <typename T>
inline image<T>::image(const vec2i& size, const T* value)
    : extent{size}, pixels(value, value + (size_t)size.x * (size_t)size.y) {}

// size
template <typename T>
inline bool image<T>::empty() const {
    return pixels.empty();
}
template <typename T>
inline vec2i image<T>::size() const {
    return extent;
}
template <typename T>
inline bool image<T>::contains(const vec2i& ij) const {
    return ij.x > 0 && ij.x < extent.x && ij.y > 0 && ij.y < extent.y;
}
template <typename T>
inline void image<T>::resize(const vec2i& size) {
    if (size == extent) return;
    extent = size;
    pixels.resize((size_t)size.x * (size_t)size.y);
}
template <typename T>
inline void image<T>::assign(const vec2i& size, const T& value) {
    extent = size;
    pixels.assign((size_t)size.x * (size_t)size.y, value);
}
template <typename T>
inline void image<T>::shrink_to_fit() {
    pixels.shrink_to_fit();
}

// element access
template <typename T>
inline T& image<T>::operator[](const vec2i& ij) {
    return pixels[ij.y * extent.x + ij.x];
}
template <typename T>
inline const T& image<T>::operator[](const vec2i& ij) const {
    return pixels[ij.y * extent.x + ij.x];
}

// data access
template <typename T>
inline T* image<T>::data() {
    return pixels.data();
}
template <typename T>
inline const T* image<T>::data() const {
    return pixels.data();
}

// iteration
template <typename T>
inline T* image<T>::begin() {
    return pixels.data();
}
template <typename T>
inline T* image<T>::end() {
    return pixels.data() + pixels.size();
}
template <typename T>
inline const T* image<T>::begin() const {
    return pixels.data();
}
template <typename T>
inline const T* image<T>::end() const {
    return pixels.data() + pixels.size();
}

// equality
template <typename T>
inline bool operator==(const image<T>& a, const image<T>& b) {
    return a.size() == b.size() && a.pixels == b.pixels;
}
template <typename T>
inline bool operator!=(const image<T>& a, const image<T>& b) {
    return a.size() != b.size() || a.pixels != b.pixels;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1 channel volume.
void load_volume(const string& filename, volume<float>& vol);
void save_volume(const string& filename, const volume<float>& vol);

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

template <typename T>
constexpr bool type_error = std::false_type::value;
template <int N>
constexpr bool channel_error = std::false_type::value;
template <int N>
constexpr bool space_error = std::false_type::value;

// Default alpha
template <typename T>
constexpr T default_alpha() {
    if constexpr (std::is_same_v<T, byte>) {
        return (byte)255;
    } else {
        return (T)1;
    }
}

inline byte float_to_byte(float a) { return (byte)clamp(int(a * 256), 0, 255); }
inline float byte_to_float(byte a) { return a / 255.0f; }
inline vec1b float_to_byte(const vec1f& a) { return {float_to_byte(a.x)}; }
inline vec2b float_to_byte(const vec2f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y)};
}
inline vec3b float_to_byte(const vec3f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z)};
}
inline vec4b float_to_byte(const vec4f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z),
        float_to_byte(a.w)};
}
inline vec1f byte_to_float(const vec1b& a) { return {byte_to_float(a.x)}; }
inline vec2f byte_to_float(const vec2b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y)};
}
inline vec3f byte_to_float(const vec3b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z),
        byte_to_float(a.w)};
}

inline float luminance(const vec3f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z);
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
inline vec1f srgb_to_linear(const vec1f& srgb) {
    return {srgb_to_linear(srgb.x)};
}
inline vec2f srgb_to_linear(const vec2f& srgb) {
    return {srgb_to_linear(srgb.x), srgb.y};
}
inline vec3f srgb_to_linear(const vec3f& srgb) {
    return {
        srgb_to_linear(srgb.x), srgb_to_linear(srgb.y), srgb_to_linear(srgb.z)};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {srgb_to_linear(srgb.x), srgb_to_linear(srgb.y),
        srgb_to_linear(srgb.z), srgb.w};
}
inline vec1f linear_to_srgb(const vec1f& lin) {
    return {linear_to_srgb(lin.x)};
}
inline vec2f linear_to_srgb(const vec2f& lin) {
    return {linear_to_srgb(lin.x), lin.y};
}
inline vec3f linear_to_srgb(const vec3f& lin) {
    return {
        linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z)};
}
inline vec4f linear_to_srgb(const vec4f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z),
        lin.w};
}

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f apply_contrast(const vec3f& rgb, float contrast, float grey) {
    return max(zero3f, grey + (rgb - grey) * (contrast * 2));
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f apply_logcontrast(const vec3f& rgb, float logcontrast, float grey) {
    auto epsilon  = (float)0.0001;
    auto log_grey = log2(grey);
    auto log_ldr  = log2(rgb + epsilon);
    auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
    return max(zero3f, exp2(adjusted) - epsilon);
}
// Apply saturation.
inline vec3f apply_saturation(
    const vec3f& rgb, float saturation, const vec3f& weights) {
    auto grey = dot(weights, rgb);
    return max(zero3f, grey + (rgb - grey) * (saturation * 2));
}

// Fitted ACES tonemapping curve.
inline float tonemap_filmic(float hdr_) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * 0.6f;  // brings it back to ACES range
    auto ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
               (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
    return max(0.0f, ldr);
}
inline vec3f tonemap_filmic(const vec3f& hdr_, bool accurate_fit) {
    if (!accurate_fit) {
        // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
        auto hdr = hdr_ * 0.6f;  // brings it back to ACES range
        auto ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
                   (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
        return max(zero3f, ldr);
    } else {
        // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
        // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
        static const auto ACESInputMat = transpose(mat3f{
            {0.59719, 0.35458, 0.04823},
            {0.07600, 0.90834, 0.01566},
            {0.02840, 0.13383, 0.83777},
        });
        // ODT_SAT => XYZ => D60_2_D65 => sRGB
        static const auto ACESOutputMat = transpose(mat3f{
            {1.60475, -0.53108, -0.07367},
            {-0.10208, 1.10813, -0.00605},
            {-0.00327, -0.07276, 1.07602},
        });
        // RRT => ODT
        auto RRTAndODTFit = [](const vec3f& v) -> vec3f {
            return (v * v + v * 0.0245786f - 0.000090537f) /
                   (v * v * 0.983729f + v * 0.4329510f + 0.238081f);
        };

        auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
        return max(zero3f, ldr);
    }
}

// Convert number of channels
template <typename T, int N1, int N2>
inline vec<T, N1> convert_channels(const vec<T, N2>& a) {
    if constexpr (N1 == 1) {
        if constexpr (N2 == 1) {
            return {a.x};
        } else if constexpr (N2 == 2) {
            return {a.x};
        } else if constexpr (N2 == 3) {
            return {(a.x + a.y + a.z) / 3};
        } else if constexpr (N2 == 4) {
            return {(a.x + a.y + a.z) / 3};
        } else {
            static_assert(channel_error<N2>, "unsupported number of channels");
        }
    } else if constexpr (N1 == 2) {
        if constexpr (N2 == 1) {
            return {a.x, default_alpha<T>()};
        } else if constexpr (N2 == 2) {
            return {a.x, a.y};
        } else if constexpr (N2 == 3) {
            return {(a.x + a.y + a.z) / 3, default_alpha<T>()};
        } else if constexpr (N2 == 4) {
            return {(a.x + a.y + a.z) / 3, a.w};
        } else {
            static_assert(channel_error<N2>, "unsupported number of channels");
        }
    } else if constexpr (N1 == 3) {
        if constexpr (N2 == 1) {
            return {a.x, a.x, a.x};
        } else if constexpr (N2 == 2) {
            return {a.x, a.x, a.x};
        } else if constexpr (N2 == 3) {
            return {a.x, a.y, a.z};
        } else if constexpr (N2 == 4) {
            return {a.x, a.y, a.z};
        } else {
            static_assert(channel_error<N2>, "unsupported number of channels");
        }
    } else if constexpr (N1 == 4) {
        if constexpr (N2 == 1) {
            return {a.x, a.x, a.x, default_alpha<T>()};
        } else if constexpr (N2 == 2) {
            return {a.x, a.x, a.x};
        } else if constexpr (N2 == 3) {
            return {a.x, a.y, a.z, default_alpha<T>()};
        } else if constexpr (N2 == 4) {
            return {a.x, a.y, a.z};
        } else {
            static_assert(channel_error<N2>, "unsupported number of channels");
        }
    } else {
        static_assert(channel_error<N1>, "unsupported number of channels");
    }
}

// https://en.wikipedia.org/wiki/SRGB
constexpr mat3f srgb_to_xyz_mat = {
    {0.4124, 0.2126, 0.0193},
    {0.3576, 0.7152, 0.1192},
    {0.1805, 0.0722, 0.9504},
};
constexpr mat3f xyz_to_srgb_mat = {
    {+3.2406, -1.5372, -0.4986},
    {-0.9689, +1.8758, +0.0415},
    {+0.0557, -0.2040, +1.0570},
};

// Curve type
enum struct rgb_curve_type {
    linear,
    gamma,
    linear_gamma,
    aces_cc,
    aces_cct,
    pq,
    hlg
};

// RGB color space definition. Various predefined color spaces are listed below.
struct rgb_space {
    // primaries
    vec2f red_chromaticity;    // xy chromaticity of the red primary
    vec2f green_chromaticity;  // xy chromaticity of the green primary
    vec2f blue_chromaticity;   // xy chromaticity of the blue primary
    vec2f white_chromaticity;  // xy chromaticity of the white point
    mat3f rgb_to_xyz_mat;      // matrix from rgb to xyz
    mat3f xyz_to_rgb_mat;      // matrix from xyz to rgb
    // tone curve
    rgb_curve_type curve_type;
    float          curve_gamma;  // gamma for power curves
    vec4f          curve_abcd;   // tone curve values for linear_gamma curves
};

// Compute the rgb -> xyz matrix from the color space definition
// Input: red, green, blue, white (x,y) chromoticities
// Algorithm from: SMPTE Recommended Practice RP 177-1993
// http://car.france3.mars.free.fr/HD/INA-%2026%20jan%2006/SMPTE%20normes%20et%20confs/rp177.pdf
constexpr mat3f rgb_to_xyz_mat(
    const vec2f& rc, const vec2f& gc, const vec2f& bc, const vec2f& wc) {
    auto rgb = mat3f{
        {rc.x, rc.y, 1 - rc.x - rc.y},
        {gc.x, gc.y, 1 - gc.x - gc.y},
        {bc.x, bc.y, 1 - bc.x - bc.y},
    };
    auto w = vec3f{wc.x, wc.y, 1 - wc.x - wc.y};
    auto c = inverse(rgb) * vec3f{w.x / w.y, 1, w.z / w.y};
    return mat3f{c.x * rgb.x, c.y * rgb.y, c.z * rgb.z};
}

// Construct an RGB color space. Predefined color spaces below
constexpr rgb_space make_linear_rgb_space(const vec2f& red, const vec2f& green,
    const vec2f& blue, const vec2f& white) {
    return rgb_space{red, green, blue, white,
        rgb_to_xyz_mat(red, green, blue, white),
        inverse(rgb_to_xyz_mat(red, green, blue, white)),
        rgb_curve_type::linear};
}
constexpr rgb_space make_gamma_rgb_space(const vec2f& red, const vec2f& green,
    const vec2f& blue, const vec2f& white, float gamma,
    const vec4f& curve_abcd = zero4f) {
    return rgb_space{red, green, blue, white,
        rgb_to_xyz_mat(red, green, blue, white),
        inverse(rgb_to_xyz_mat(red, green, blue, white)),
        curve_abcd == zero4f ? rgb_curve_type::gamma
                             : rgb_curve_type::linear_gamma};
}
constexpr rgb_space make_other_rgb_space(const vec2f& red, const vec2f& green,
    const vec2f& blue, const vec2f& white, rgb_curve_type curve_type) {
    return rgb_space{red, green, blue, white,
        rgb_to_xyz_mat(red, green, blue, white),
        inverse(rgb_to_xyz_mat(red, green, blue, white)), curve_type};
}

constexpr rgb_space get_rgb_space(rgb_space_type space) {
    switch (space) {
        // https://en.wikipedia.org/wiki/Rec._709
        case rgb_space_type::linear_srgb:
            return make_linear_rgb_space({0.6400, 0.3300}, {0.3000, 0.6000},
                {0.1500, 0.0600}, {0.3127, 0.3290});
        // https://en.wikipedia.org/wiki/Rec._709
        case rgb_space_type::srgb:
            return make_gamma_rgb_space({0.6400, 0.3300}, {0.3000, 0.6000},
                {0.1500, 0.0600}, {0.3127, 0.3290}, 2.4,
                {1.055, 0.055, 12.92, 0.0031308});
        // https://en.wikipedia.org/wiki/Academy_Color_Encoding_System
        case rgb_space_type::aces_2065:
            return make_linear_rgb_space({0.7347, 0.2653}, {0.0000, 1.0000},
                {0.0001, -0.0770}, {0.32168, 0.33767});
        // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
        case rgb_space_type::aces_cg:
            return make_linear_rgb_space({0.7130, 0.2930}, {0.1650, 0.8300},
                {0.1280, +0.0440}, {0.32168, 0.33767});
        // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
        case rgb_space_type::aces_cc:
            return make_other_rgb_space({0.7130, 0.2930}, {0.1650, 0.8300},
                {0.1280, +0.0440}, {0.32168, 0.33767}, rgb_curve_type::aces_cc);
        // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
        case rgb_space_type::aces_cct:
            return make_other_rgb_space({0.7130, 0.2930}, {0.1650, 0.8300},
                {0.1280, +0.0440}, {0.32168, 0.33767},
                rgb_curve_type::aces_cct);
        // https://en.wikipedia.org/wiki/Adobe_RGB_color_space
        case rgb_space_type::adobe_rgb:
            return make_gamma_rgb_space({0.6400, 0.3300}, {0.2100, 0.7100},
                {0.1500, 0.0600}, {0.3127, 0.3290}, 2.19921875);
        // https://en.wikipedia.org/wiki/Rec._709
        case rgb_space_type::rec_709:
            return make_gamma_rgb_space({0.6400, 0.3300}, {0.3000, 0.6000},
                {0.1500, 0.0600}, {0.3127, 0.3290}, 1 / 0.45,
                {1.099, 0.099, 4.500, 0.018});
        // https://en.wikipedia.org/wiki/Rec._2020
        case rgb_space_type::rec_2020:
            return make_gamma_rgb_space({0.7080, 0.2920}, {0.1700, 0.7970},
                {0.1310, 0.0460}, {0.3127, 0.3290}, 1 / 0.45,
                {1.09929682680944, 0.09929682680944, 4.5, 0.018053968510807});
        // https://en.wikipedia.org/wiki/Rec._2020
        case rgb_space_type::rec_2100_pq:
            return make_other_rgb_space({0.7080, 0.2920}, {0.1700, 0.7970},
                {0.1310, 0.0460}, {0.3127, 0.3290}, rgb_curve_type::pq);
        // https://en.wikipedia.org/wiki/Rec._2020
        case rgb_space_type::rec_2100_hlg:
            return make_other_rgb_space({0.7080, 0.2920}, {0.1700, 0.7970},
                {0.1310, 0.0460}, {0.3127, 0.3290}, rgb_curve_type::hlg);
        // https://en.wikipedia.org/wiki/DCI-P3
        case rgb_space_type::p3_dci:
            return make_gamma_rgb_space({0.6800, 0.3200}, {0.2650, 0.6900},
                {0.1500, 0.0600}, {0.3140, 0.3510}, 1.6);
        // https://en.wikipedia.org/wiki/DCI-P3
        case rgb_space_type::p3_d60:
            return make_gamma_rgb_space({0.6800, 0.3200}, {0.2650, 0.6900},
                {0.1500, 0.0600}, {0.32168, 0.33767}, 1.6);
        // https://en.wikipedia.org/wiki/DCI-P3
        case rgb_space_type::p3_d65:
            return make_gamma_rgb_space({0.6800, 0.3200}, {0.2650, 0.6900},
                {0.1500, 0.0600}, {0.3127, 0.3290}, 1.6);
        // https://en.wikipedia.org/wiki/DCI-P3
        case rgb_space_type::p3_display:
            return make_gamma_rgb_space({0.6800, 0.3200}, {0.2650, 0.6900},
                {0.1500, 0.0600}, {0.3127, 0.3290}, 2.4,
                {1.055, 0.055, 12.92, 0.0031308});
        // https://en.wikipedia.org/wiki/ProPhoto_RGB_color_space
        case rgb_space_type::prophoto_rgb:
            return make_gamma_rgb_space({0.7347, 0.2653}, {0.1596, 0.8404},
                {0.0366, 0.0001}, {0.3457, 0.3585}, 1.8,
                {1, 0, 16, 0.001953125});
        default: throw "unknown color space";
    }
}

// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
constexpr float acescc_display_to_linear(float x) {
    if (x < -0.3013698630f) {  // (9.72-15)/17.52
        return (exp2(x * 17.52f - 9.72f) - exp2(-16.0f)) * 2;
    } else if (x < (log2(65504.0f) + 9.72f) / 17.52f) {
        return exp2(x * 17.52f - 9.72f);
    } else {  // (in >= (log2(65504)+9.72)/17.52)
        return 65504.0f;
    }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
constexpr float acescct_display_to_linear(float x) {
    if (x < 0.155251141552511f) {
        return (x - 0.0729055341958355f) / 10.5402377416545f;
    } else {
        return exp2(x * 17.52f - 9.72f);
    }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
constexpr float acescc_linear_to_display(float x) {
    if (x <= 0) {
        return -0.3584474886f;  // =(log2( pow(2.,-16.))+9.72)/17.52
    } else if (x < exp2(-15.0f)) {
        return (log2(exp2(-16.0f) + x * 0.5f) + 9.72f) / 17.52f;
    } else {  // (in >= pow(2.,-15))
        return (log2(x) + 9.72f) / 17.52f;
    }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
constexpr float acescct_linear_to_display(float x) {
    if (x <= 0.0078125f) {
        return 10.5402377416545f * x + 0.0729055341958355f;
    } else {
        return (log2(x) + 9.72f) / 17.52f;
    }
}

// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// https://github.com/ampas/aces-dev/blob/master/transforms/ctl/lib/ACESlib.Utilities_Color.ctl
// In PQ, we assume that the linear luminance in [0,1] corresponds to
// [0,10000] cd m^2
inline float pq_display_to_linear(float x) {
    auto Np = pow(x, 1 / 78.84375f);
    auto L  = max(Np - 0.8359375f, 0.0f);
    L       = L / (18.8515625f - 18.6875f * Np);
    L       = pow(L, 1 / 0.1593017578125f);
    return L;
}
inline float pq_linear_to_display(float x) {
    return pow((0.8359375 + 18.8515625 * pow(x, 0.1593017578125f)) /
                   (1 + 18.6875f * pow(x, 0.1593017578125f)),
        78.84375f);
}
// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// In HLG, we assume that the linear luminance in [0,1] corresponds to
// [0,1000] cd m^2. Note that the version we report here is scaled in [0,1]
// range for nominal luminance. But HLG was initially defined in the [0,12]
// range where it maps 1 to 0.5 and 12 to 1. For use in HDR tonemapping that is
// likely a better range to use.
constexpr float hlg_display_to_linear(float x) {
    if (x < 0.5f) {
        return 3 * 3 * x * x;
    } else {
        return (exp((x - 0.55991073f) / 0.17883277f) + 0.28466892f) / 12;
    }
}
constexpr float hlg_linear_to_display(float x) {
    if (x < 1 / 12.0f) {
        return sqrt(3 * x);
    } else {
        return 0.17883277f * log(12 * x - 0.28466892f) + 0.55991073f;
    }
}

// Applies linear to display transformations and vice-verse
constexpr vec3f linear_to_display(const vec3f& rgb, const rgb_space& space) {
    if (space.curve_type == rgb_curve_type::linear) {
        return rgb;
    } else if (space.curve_type == rgb_curve_type::gamma) {
        return pow(rgb, 1 / space.curve_gamma);
    } else if (space.curve_type == rgb_curve_type::linear_gamma) {
        auto& [a, b, c, d] = space.curve_abcd;
        auto lim           = d;
        auto lin           = rgb * c;
        auto gamma         = a * pow(rgb, 1 / space.curve_gamma) - b;
        return {
            rgb.x < lim ? lin.x : gamma.x,
            rgb.y < lim ? lin.y : gamma.y,
            rgb.z < lim ? lin.z : gamma.z,
        };
    } else if (space.curve_type == rgb_curve_type::aces_cc) {
        return {acescc_linear_to_display(rgb.x),
            acescc_linear_to_display(rgb.y), acescc_linear_to_display(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::aces_cct) {
        return {acescct_linear_to_display(rgb.x),
            acescct_linear_to_display(rgb.y), acescct_linear_to_display(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::pq) {
        return {pq_linear_to_display(rgb.x), pq_linear_to_display(rgb.y),
            pq_linear_to_display(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::hlg) {
        return {hlg_linear_to_display(rgb.x), hlg_linear_to_display(rgb.y),
            hlg_linear_to_display(rgb.z)};
    } else {
        throw runtime_error("should not have gotten here");
    }
}
constexpr vec3f display_to_linear(const vec3f& rgb, const rgb_space& space) {
    if (space.curve_type == rgb_curve_type::linear) {
        return rgb;
    } else if (space.curve_type == rgb_curve_type::gamma) {
        return pow(rgb, space.curve_gamma);
    } else if (space.curve_type == rgb_curve_type::linear_gamma) {
        auto& [a, b, c, d] = space.curve_abcd;
        auto lim           = 1 / d;
        auto lin           = rgb / c;
        auto gamma         = pow((rgb + b) / a, space.curve_gamma);
        return {
            rgb.x < lim ? lin.x : gamma.x,
            rgb.y < lim ? lin.y : gamma.y,
            rgb.z < lim ? lin.z : gamma.z,
        };
    } else if (space.curve_type == rgb_curve_type::aces_cc) {
        return {acescc_display_to_linear(rgb.x),
            acescc_display_to_linear(rgb.y), acescc_display_to_linear(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::aces_cct) {
        return {acescct_display_to_linear(rgb.x),
            acescct_display_to_linear(rgb.y), acescct_display_to_linear(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::pq) {
        return {pq_display_to_linear(rgb.x), pq_display_to_linear(rgb.y),
            pq_display_to_linear(rgb.z)};
    } else if (space.curve_type == rgb_curve_type::hlg) {
        return {hlg_display_to_linear(rgb.x), hlg_display_to_linear(rgb.y),
            hlg_display_to_linear(rgb.z)};
    } else {
        throw runtime_error("should not have gotten here");
    }
}

// Conversion between rgb color spaces
inline vec3f rgb_to_rgb(
    const vec3f& rgb, const rgb_space& from, const rgb_space& to) {
    if (memcmp(&from, &to, sizeof(from)) == 0) {
        return rgb;
    } else if (from.rgb_to_xyz_mat == to.rgb_to_xyz_mat) {
        return linear_to_display(display_to_linear(rgb, from), to);
    } else {
        return xyz_to_rgb(rgb_to_xyz(rgb, from), to);
    }
}

// Conversion to/from xyz
inline vec3f rgb_to_xyz(const vec3f& rgb, const rgb_space& from) {
    return from.rgb_to_xyz_mat * display_to_linear(rgb, from);
}
inline vec3f xyz_to_rgb(const vec3f& xyz, const rgb_space& to) {
    return linear_to_display(to.xyz_to_rgb_mat * xyz, to);
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

// Approximate color of blackbody radiation from wavelength in nm.
inline vec3f blackbody_to_rgb(float temperature) {
    // https://github.com/neilbartlett/color-temperature
    auto rgb = zero3f;
    if ((temperature / 100) < 66) {
        rgb.x = 255;
    } else {
        // a + b x + c Log[x] /.
        // {a -> 351.97690566805693`,
        // b -> 0.114206453784165`,
        // c -> -40.25366309332127
        // x -> (kelvin/100) - 55}
        rgb.x = (temperature / 100) - 55;
        rgb.x = 351.97690566805693f + 0.114206453784165f * rgb.x -
                40.25366309332127f * log(rgb.x);
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
        rgb.y = -155.25485562709179f - 0.44596950469579133f * rgb.y +
                104.49216199393888f * log(rgb.y);
        if (rgb.y < 0) rgb.y = 0;
        if (rgb.y > 255) rgb.y = 255;
    } else {
        // a + b x + c Log[x] /.
        // {a -> 325.4494125711974`,
        // b -> 0.07943456536662342`,
        // c -> -28.0852963507957`,
        // x -> (kelvin/100) - 50}
        rgb.y = (temperature / 100) - 50;
        rgb.y = 325.4494125711974f + 0.07943456536662342f * rgb.y -
                28.0852963507957f * log(rgb.y);
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
            rgb.z = -254.76935184120902f + 0.8274096064007395f * rgb.z +
                    115.67994401066147f * log(rgb.z);
            if (rgb.z < 0) rgb.z = 0;
            if (rgb.z > 255) rgb.z = 255;
        }
    }

    return srgb_to_linear(rgb / 255);
}

// Convert HSV to RGB
inline vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0) return {v, v, v};

    h       = fmod(h, 1.0f) / (60.0f / 360.0f);
    int   i = (int)h;
    float f = h - (float)i;
    float p = v * (1 - s);
    float q = v * (1 - s * f);
    float t = v * (1 - s * (1 - f));

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
        K = -1;
    }
    if (r < g) {
        swap(r, g);
        K = -2 / 6.0f - K;
    }

    float chroma = r - (g < b ? g : b);
    return {
        fabsf(K + (g - b) / (6 * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Gets pixels in an image region
template <typename T>
inline void get_region(
    image<T>& clipped, const image<T>& img, const image_region& region) {
    clipped.resize(region.size());
    for (auto j = 0; j < region.size().y; j++) {
        for (auto i = 0; i < region.size().x; i++) {
            clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
        }
    }
}
template <typename T>
inline void set_region(
    image<T>& img, const image<T>& region, const vec2i& offset) {
    for (auto j = 0; j < region.size().y; j++) {
        for (auto i = 0; i < region.size().x; i++) {
            if (!img.contains({i, j})) continue;
            img[vec2i{i, j} + offset] = region[{i, j}];
        }
    }
}

// Splits an image into an array of regions
inline void make_regions(vector<image_region>& regions, const vec2i& size,
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
        shuffle(regions, rng);
    }
}

// Apply a function to each image pixel
template <typename T1, typename T2, typename Func>
inline void apply_image(
    image<T1>& result, const image<T2>& source, const Func& func) {
    result.resize(source.size());
    for (auto j = 0; j < result.size().y; j++) {
        for (auto i = 0; i < result.size().x; i++) {
            result[{i, j}] = func(source[{i, j}]);
        }
    }
}
template <typename T1, typename T2, typename Func>
inline void apply_image(image<T1>& result, const image<T2>& source,
    const image_region& region, const Func& func) {
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
    return apply_image(fl, bt, [](const auto& a) { return byte_to_float(a); });
}
template <typename T, typename TB>
inline void float_to_byte(image<TB>& bt, const image<T>& fl) {
    return apply_image(bt, fl, [](const auto& a) { return float_to_byte(a); });
}

// Conversion between linear and gamma-encoded images.
template <typename T>
inline void srgb_to_linear(image<T>& lin, const image<T>& srgb) {
    return apply_image(
        lin, srgb, [](const auto& a) { return srgb_to_linear(a); });
}
template <typename T>
inline void linear_to_srgb(image<T>& srgb, const image<T>& lin) {
    return apply_image(
        srgb, lin, [](const auto& a) { return linear_to_srgb(a); });
}
template <typename T, typename TB>
inline void srgb_to_linear(image<T>& lin, const image<TB>& srgb) {
    return apply_image(lin, srgb,
        [](const auto& a) { return srgb_to_linear(byte_to_float(a)); });
}
template <typename T, typename TB>
inline void linear_to_srgb(image<TB>& srgb, const image<T>& lin) {
    return apply_image(srgb, lin,
        [](const auto& a) { return float_to_byte(linear_to_srgb(a)); });
}

inline vec3f tonemap(const vec3f& hdr, const tonemap_params& params) {
    auto rgb = hdr;
    if (params.exposure != 0) rgb *= exp2(params.exposure);
    if (params.tint != vec3f{1, 1, 1}) rgb *= params.tint;
    if (params.contrast != 0.5f)
        rgb = apply_contrast(rgb, params.contrast, 0.18f);
    if (params.logcontrast != 0.5f)
        rgb = apply_logcontrast(rgb, params.logcontrast, 0.18f);
    if (params.saturation != 0.5f)
        rgb = apply_saturation(rgb, params.saturation);
    if (params.filmic) rgb = tonemap_filmic(rgb);
    if (params.srgb) rgb = linear_to_srgb(rgb);
    return rgb;
}

// Apply exposure and filmic tone mapping
inline void tonemap(
    image<vec3f>& ldr, const image<vec3f>& hdr, const tonemap_params& params) {
    return apply_image(
        ldr, hdr, [params](const vec3f& hdr) { return tonemap(hdr, params); });
}
inline void tonemap(
    image<vec4f>& ldr, const image<vec4f>& hdr, const tonemap_params& params) {
    return apply_image(ldr, hdr,
        [scale = exp2(params.exposure) * params.tint, params](
            const vec4f& hdr) {
            return vec4f{tonemap(xyz(hdr), params), hdr.w};
        });
}
inline void tonemap(
    image<vec3b>& ldr, const image<vec3f>& hdr, const tonemap_params& params) {
    return apply_image(ldr, hdr, [params](const vec3f& hdr) {
        return float_to_byte(tonemap(hdr, params));
    });
}
inline void tonemap(
    image<vec4b>& ldr, const image<vec4f>& hdr, const tonemap_params& params) {
    return apply_image(ldr, hdr, [params](const vec4f& hdr) {
        return float_to_byte(vec4f{tonemap(xyz(hdr), params), hdr.w});
    });
}
inline void tonemap(image<vec4f>& ldr, const image<vec4f>& hdr,
    const image_region& region, const tonemap_params& params) {
    return apply_image(ldr, hdr, region, [params](const vec4f& hdr) {
        return vec4f{tonemap(xyz(hdr), params), hdr.w};
    });
}

inline vec3f colorgrade(const vec3f& ldr, const colorgrade_params& params) {
    auto rgb = ldr;
    if (params.contrast != 0.5f) {
        rgb = gain(ldr, 1 - params.contrast);
    }
    if (params.shadows != 0.5f || params.midtones != 0.5f ||
        params.highlights != 0.5f || params.shadows_color != vec3f{1, 1, 1} ||
        params.midtones_color != vec3f{1, 1, 1} ||
        params.highlights_color != vec3f{1, 1, 1}) {
        auto lift  = params.shadows_color;
        auto gamma = params.midtones_color;
        auto gain  = params.highlights_color;

        lift      = lift - mean(lift) + params.shadows - (float)0.5;
        gain      = gain - mean(gain) + params.highlights + (float)0.5;
        auto grey = gamma - mean(gamma) + params.midtones;
        gamma     = log(((float)0.5 - lift) / (gain - lift)) / log(grey);

        // apply_image
        auto lerp_value = clamp01(pow(rgb, 1 / gamma));
        rgb             = gain * lerp_value + lift * (1 - lerp_value);
    }
    return rgb;
}

// Apply exposure and filmic tone mapping
inline void colorgrade(image<vec4f>& corrected, const image<vec4f>& ldr,
    const image_region& region, const colorgrade_params& params) {
    return apply_image(corrected, ldr, region, [&params](const vec4f& hdr) {
        return vec4f{colorgrade(xyz(hdr), params), hdr.w};
    });
}

// compute white balance
inline vec3f compute_white_balance(const image<vec4f>& img) {
    auto rgb = zero3f;
    for (auto& p : img) rgb += xyz(p);
    if (rgb == zero3f) return zero3f;
    return rgb / max(rgb);
}

template <int N1, int N2>
inline void convert_channels(
    image<vec<float, N1>>& result, const image<vec<float, N2>>& source) {
    return apply_image(result, source, [](const vec<float, N2>& a) {
        return convert_channels<float, N1, N2>(a);
    });
}
template <int N1, int N2>
inline void convert_channels(
    image<vec<byte, N1>>& result, const image<vec<byte, N2>>& source) {
    return apply_image(result, source, [](const vec<byte, N2>& a) {
        return convert_channels<byte, N1, N2>(a);
    });
}

// Resize image.
template <int N>
static inline void resize_impl(
    image<vec<float, N>>& res_img, const image<vec<float, N>>& img) {
    auto alpha = (N == 2 || N == 4) ? N - 1 : -1;
    stbir_resize_float_generic((float*)img.data(), img.size().x, img.size().y,
        sizeof(vec<float, N>) * img.size().x, (float*)res_img.data(),
        res_img.size().x, res_img.size().y,
        sizeof(vec<float, N>) * res_img.size().x, N, alpha, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
}
template <int N>
static inline void resize_impl(
    image<vec<byte, N>>& res_img, const image<vec<byte, N>>& img) {
    auto alpha = (N == 2 || N == 4) ? N - 1 : -1;
    stbir_resize_uint8_generic((byte*)img.data(), img.size().x, img.size().y,
        sizeof(vec<byte, N>) * img.size().x, (byte*)res_img.data(),
        res_img.size().x, res_img.size().y,
        sizeof(vec<byte, N>) * res_img.size().x, N, alpha, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
}

inline void resize(image<float>& res, const image<float>& img) {
    return resize_impl((image<vec1f>&)res, (const image<vec1f>&)img);
}
inline void resize(image<vec1f>& res, const image<vec1f>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec2f>& res, const image<vec2f>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec3f>& res, const image<vec3f>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec4f>& res, const image<vec4f>& img) {
    return resize_impl(res, img);
}

inline void resize(image<byte>& res, const image<byte>& img) {
    return resize_impl((image<vec1b>&)res, (const image<vec1b>&)img);
}
inline void resize(image<vec1b>& res, const image<vec1b>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec2b>& res, const image<vec2b>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec3b>& res, const image<vec3b>& img) {
    return resize_impl(res, img);
}
inline void resize(image<vec4b>& res, const image<vec4b>& img) {
    return resize_impl(res, img);
}

template <typename T>
inline void resize(image<T>& res_img, const image<T>& img, const vec2i& size_) {
    auto size = size_;
    if (size == zero2i) {
        throw std::invalid_argument("bad image size in resize");
    }
    if (size.y == 0) {
        size.y = (int)round(size.x * (float)img.size().y / (float)img.size().x);
    } else if (size.x == 0) {
        size.x = (int)round(size.y * (float)img.size().x / (float)img.size().y);
    }
    res_img = {size};
    resize(res_img, img);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME
// -----------------------------------------------------------------------------
namespace yocto {

// constructors
template <typename T>
inline volume<T>::volume() : extent{0, 0, 0}, voxels{} {}
template <typename T>
inline volume<T>::volume(const vec3i& size, const T& value)
    : extent{size}
    , voxels((size_t)size.x * (size_t)size.y * (size_t)size.z, value) {}
template <typename T>
inline volume<T>::volume(const vec3i& size, const T* value)
    : extent{size}
    , voxels(value, value + (size_t)size.x * (size_t)size.y * (size_t)size.z) {}

// size
template <typename T>
inline bool volume<T>::empty() const {
    return voxels.empty();
}
template <typename T>
inline vec3i volume<T>::size() const {
    return extent;
}
template <typename T>
inline void volume<T>::resize(const vec3i& size) {
    if (size == extent) return;
    extent = size;
    voxels.resize((size_t)size.x * (size_t)size.y * (size_t)size.z);
}
template <typename T>
inline void volume<T>::assign(const vec3i& size, const T& value) {
    extent = size;
    voxels.assign((size_t)size.x * (size_t)size.y * (size_t)size.z, value);
}
template <typename T>
inline void volume<T>::shrink_to_fit() {
    voxels.shrink_to_fit();
}

// element access
template <typename T>
inline T& volume<T>::operator[](const vec3i& ijk) {
    return voxels[ijk.z * extent.x * extent.y + ijk.y * extent.x + ijk.x];
}
template <typename T>
inline const T& volume<T>::operator[](const vec3i& ijk) const {
    return voxels[ijk.z * extent.x * extent.y + ijk.y * extent.x + ijk.x];
}

// data access
template <typename T>
inline T* volume<T>::data() {
    return voxels.data();
}
template <typename T>
inline const T* volume<T>::data() const {
    return voxels.data();
}

// iteration
template <typename T>
inline T* volume<T>::begin() {
    return voxels.data();
}
template <typename T>
inline T* volume<T>::end() {
    return voxels.data() + voxels.size();
}
template <typename T>
inline const T* volume<T>::begin() const {
    return voxels.data();
}
template <typename T>
inline const T* volume<T>::end() const {
    return voxels.data() + voxels.size();
}

}  // namespace yocto

#endif
