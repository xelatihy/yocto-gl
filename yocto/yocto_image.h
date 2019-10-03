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
// 5. make various image examples with the `make_proc_image()` functions
// 6. create procedural sun-sky images with `make_sunsky()`
// 7. many color conversion functions are available in the code below
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

#include "yocto_common.h"
#include "yocto_math.h"

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
  bool   empty() const;
  vec2i  size() const;
  size_t count() const;
  bool   contains(const vec2i& ij) const;
  void   resize(const vec2i& size);
  void   assign(const vec2i& size, const T& value = {});
  void   shrink_to_fit();

  // element access
  T&       operator[](int i);
  const T& operator[](int i) const;
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
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates a color image at a point `uv`.
inline vec4f eval_image(const image<vec4f>& img, const vec2f& uv,
    bool no_interpolation = false, bool clamp_to_edge = false);
inline vec4f eval_image(const image<vec4b>& img, const vec2f& uv,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image region
struct image_region {
  vec2i min = zero2i;
  vec2i max = zero2i;

  image_region();
  image_region(const vec2i& min, const vec2i& max);

  vec2i size() const;
};

// Splits an image into an array of regions
vector<image_region> make_image_regions(
    const vec2i& size, int region_size = 32, bool shuffled = false);

// Gets pixels in an image region
template <typename T>
inline image<T> get_image_region(
    const image<T>& img, const image_region& region);
template <typename T>
inline void set_region(
    image<T>& img, const image<T>& region, const vec2i& offset);
template <typename T>
inline void get_region(
    image<T>& clipped, const image<T>& img, const image_region& region);

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt);
image<vec4b> float_to_byte(const image<vec4f>& fl);
void         byte_to_float(image<vec4f>& fl, const image<vec4b>& bt);
void         float_to_byte(image<vec4b>& bt, const image<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_rgb(const image<vec4f>& srgb);
image<vec4f> rgb_to_srgb(const image<vec4f>& rgb);
image<vec4f> srgb_to_rgb(const image<vec4b>& srgb);
image<vec4b> rgb_to_srgbb(const image<vec4f>& rgb);
void         srgb_to_rgb(image<vec4f>& rgb, const image<vec4f>& srgb);
void         rgb_to_srgb(image<vec4f>& srgb, const image<vec4f>& rgb);

// Tone mapping params
struct tonemap_params {
  float exposure    = 0;
  vec3f tint        = {1, 1, 1};
  float contrast    = 0.5;
  float logcontrast = 0.5;
  float saturation  = 0.5;
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
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, const tonemap_params& params);
image<vec4b> tonemap_imageb(
    const image<vec4f>& hdr, const tonemap_params& params);
void tonemap_region(image<vec4f>& ldr, const image<vec4f>& hdr,
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
image<vec4f> colorgrade_image(
    const image<vec4f>& img, const colorgrade_params& params);
void colorgrade_region(image<vec4f>& corrected, const image<vec4f>& img,
    const image_region& region, const colorgrade_params& params);

// determine white balance colors
vec3f compute_white_balance(const image<vec4f>& img);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size);
image<vec4b> resize_image(const image<vec4b>& img, const vec2i& size);
void         resize_image(
            image<vec4f>& res, const image<vec4f>& img, const vec2i& size);
void resize_image(
    image<vec4b>& res, const image<vec4b>& img, const vec2i& size);

// Compute the difference between two images
image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool disply_diff);
void image_difference(image<vec4f>& diff, const image<vec4f>& a,
    const image<vec4f>& b, bool disply_diff);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Loads/saves a 4 channels float/byte image in linear color space.
image<vec4f> load_image(const string& filename);
void         load_image(const string& filename, image<vec4f>& img);
void         save_image(const string& filename, const image<vec4f>& img);
image<vec4b> load_imageb(const string& filename);
void         load_imageb(const string& filename, image<vec4b>& img);
void         save_imageb(const string& filename, const image<vec4b>& img);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Parameters for make_proc_image
struct proc_image_params {
  // clang-format off
  enum struct type_t {
    grid, checker, bumps, ramp, gammaramp, uvramp, uvgrid, blackbody, noise,
    turbulence, fbm, ridge };
  // clang-format on
  type_t type    = type_t::grid;
  vec2i  size    = {1024, 1024};
  float  scale   = 1;
  vec4f  color0  = {0, 0, 0, 1};
  vec4f  color1  = {1, 1, 1, 1};
  vec4f  noise   = {2, 0.5, 8, 1};  // lacunarity, gain, octaves, offset
  float  borderw = 0;
  vec4f  borderc = {0, 0, 0, 1};
};

// Make an image
image<vec4f> make_proc_image(const proc_image_params& params);
void make_proc_image(image<vec4f>& img, const proc_image_params& params);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
image<vec4f> make_sunsky(const vec2i& size, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
    float sun_radius = 1, const vec3f& ground_albedo = {0.2, 0.2, 0.2});
void         make_sunsky(image<vec4f>& img, const vec2i& size, float sun_angle,
            float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
            float sun_radius = 1, const vec3f& ground_albedo = {0.2, 0.2, 0.2});
// Make an image of multiple lights.
image<vec4f> make_lights(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);
void         make_lights(image<vec4f>& img, const vec2i& size,
            const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
            float lwidth = pif / 16, float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
image<vec4f> bump_to_normal(const image<vec4f>& img, float scale = 1);
void         bump_to_normal(
            image<vec4f>& norm, const image<vec4f>& img, float scale = 1);

// Add a border to an image
image<vec4f> add_border(const image<vec4f>& img, int width, const vec4f& color);
void         add_border(
            image<vec4f>& bordered, image<vec4f>& img, int width, const vec4f& color);

// Make logo images. Image is resized to proper size.
image<vec4b> make_logo(const string& name);
void         make_logo(image<vec4f>& img, const string& name);
void         make_logo(image<vec4b>& img, const string& name);
image<vec4f> add_logo(
    const image<vec4f>& img, const string& name = "logo-medium");
image<vec4b> add_logo(
    const image<vec4b>& img, const string& name = "logo-medium");
void add_logo(image<vec4f>& with_logo, const image<vec4f>& img,
    const string& name = "logo-medium");
void add_logo(image<vec4b>& with_logo, const image<vec4b>& img,
    const string& name = "logo-medium");

// Make an image preset, useful for testing. See implementation for types.
image<vec4f> make_image_preset(const string& type);
void         make_image_preset(image<vec4f>& img, const string& type);
void         make_image_preset(image<vec4b>& img, const string& type);
void         make_image_preset(
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
  volume();
  volume(const vec3i& size, const T& value);
  volume(const vec3i& size, const T* value);

  // size
  bool   empty() const;
  vec3i  size() const;
  size_t count() const;
  void   resize(const vec3i& size);
  void   assign(const vec3i& size, const T& value);
  void   shrink_to_fit();

  // element access
  T&       operator[](size_t i);
  const T& operator[](size_t i) const;
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
void make_voltest(volume<float>& vol, const vec3i& size, float scale = 10,
    float exponent = 6);
void make_volume_preset(volume<float>& vol, const string& type);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates a color image at a point `uv`.
inline float eval_volume(const image<float>& img, const vec3f& uvw,
    bool no_interpolation = false, bool clamp_to_edge = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME IO
// -----------------------------------------------------------------------------
namespace yocto {

// Loads/saves a 1 channel volume.
void load_volume(const string& filename, volume<float>& vol);
void save_volume(const string& filename, const volume<float>& vol);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion between flots and bytes
inline vec4b float_to_byte(const vec4f& a);
inline vec4f byte_to_float(const vec4b& a);

// Luminance
inline float luminance(const vec3f& a);

// sRGB non-linear curve
inline float srgb_to_rgb(float srgb);
inline float rgb_to_srgb(float rgb);
inline vec3f srgb_to_rgb(const vec3f& srgb);
inline vec4f srgb_to_rgb(const vec4f& srgb);
inline vec3f rgb_to_srgb(const vec3f& rgb);
inline vec4f rgb_to_srgb(const vec4f& rgb);

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f contrast(const vec3f& rgb, float contrast, float grey);
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f logcontrast(const vec3f& rgb, float logcontrast, float grey);
// Apply saturation.
inline vec3f saturate(const vec3f& rgb, float saturation,
    const vec3f& weights = vec3f{0.333333f});

// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb);
inline vec3f xyz_to_rgb(const vec3f& xyz);

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz);
inline vec3f xyY_to_xyz(const vec3f& xyY);

// Approximate color of blackbody radiation from wavelength in nm.
vec3f blackbody_to_rgb(float temperature);

// Converts between HSV and RGB color spaces.
vec3f hsv_to_rgb(const vec3f& hsv);
vec3f rgb_to_hsv(const vec3f& rgb);

// RGB color spaces
enum struct color_space {
  rgb,         // default linear space (srgb linear)
  srgb,        // srgb color space (non-linear)
  adobe,       // Adobe rgb color space (non-linear)
  prophoto,    // ProPhoto Kodak rgb color space (non-linear)
  rec709,      // hdtv color space (non-linear)
  rec2020,     // uhtv color space (non-linear)
  rec2100pq,   // hdr color space with perceptual quantizer (non-linear)
  rec2100hlg,  // hdr color space with hybrid log gamma (non-linear)
  aces2065,    // ACES storage format (linear)
  acescg,      // ACES CG computation (linear)
  acescc,      // ACES color correction (non-linear)
  acescct,     // ACES color correction 2 (non-linear)
  p3dci,       // P3 DCI (non-linear)
  p3d60,       // P3 variation for D60 (non-linear)
  p3d65,       // P3 variation for D65 (non-linear)
  p3display,   // Apple display P3
};

// Conversion between rgb color spaces
vec3f color_to_xyz(const vec3f& col, color_space from);
vec3f xyz_to_color(const vec3f& xyz, color_space to);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image container ----------------

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
inline size_t image<T>::count() const {
  return pixels.size();
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
inline T& image<T>::operator[](int i) {
  return pixels[i];
}
template <typename T>
inline const T& image<T>::operator[](int i) const {
  return pixels[i];
}
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
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup an image at coordinates `ij`
inline vec4f lookup_image(
    const image<vec4f>& img, const vec2i& ij, bool as_linear) {
  return img[ij];
}
inline vec4f lookup_image(
    const image<vec4b>& img, const vec2i& ij, bool as_linear) {
  if (as_linear) {
    return byte_to_float(img[ij]);
  } else {
    return srgb_to_rgb(byte_to_float(img[ij]));
  }
}

// Evaluate a texture
template <typename T>
inline vec4f eval_image_generic(const image<T>& img, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (img.empty()) return zero4f;

  // get image width/height
  auto size = img.size();

  // get coordinates normalized for tiling
  auto s = 0.0f, t = 0.0f;
  if (clamp_to_edge) {
    s = clamp(uv.x, 0.0f, 1.0f) * size.x;
    t = clamp(uv.y, 0.0f, 1.0f) * size.y;
  } else {
    s = fmod(uv.x, 1.0f) * size.x;
    if (s < 0) s += size.x;
    t = fmod(uv.y, 1.0f) * size.y;
    if (t < 0) t += size.y;
  }

  // get image coordinates and residuals
  auto i = clamp((int)s, 0, size.x - 1), j = clamp((int)t, 0, size.y - 1);
  auto ii = (i + 1) % size.x, jj = (j + 1) % size.y;
  auto u = s - i, v = t - j;

  if (no_interpolation) return lookup_image(img, {i, j}, as_linear);

  // handle interpolation
  return lookup_image(img, {i, j}, as_linear) * (1 - u) * (1 - v) +
         lookup_image(img, {i, jj}, as_linear) * (1 - u) * v +
         lookup_image(img, {ii, j}, as_linear) * u * (1 - v) +
         lookup_image(img, {ii, jj}, as_linear) * u * v;
}

// Evaluates a color image at a point `uv`.
inline vec4f eval_image(const image<vec4f>& img, const vec2f& uv,
    bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic(img, uv, false, no_interpolation, clamp_to_edge);
}
inline vec4f eval_image(const image<vec4b>& img, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic(
      img, uv, as_linear, no_interpolation, clamp_to_edge);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

inline image_region::image_region() {}
inline image_region::image_region(const vec2i& min, const vec2i& max)
    : min{min}, max{max} {}

inline vec2i image_region::size() const { return max - min; }

// Gets pixels in an image region
template <typename T>
inline image<T> get_image_region(
    const image<T>& img, const image_region& region) {
  auto clipped = image<T>{region.size()};
  for (auto j = 0; j < region.size().y; j++) {
    for (auto i = 0; i < region.size().x; i++) {
      clipped[{i, j}] = img[{i + region.min.x, j + region.min.y}];
    }
  }
  return clipped;
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Volume container ----------

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
inline size_t volume<T>::count() const {
  return voxels.size();
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
inline T& volume<T>::operator[](size_t i) {
  return voxels[i];
}
template <typename T>
inline const T& volume<T>::operator[](size_t i) const {
  return voxels[i];
}
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

// equality
template <typename T>
inline bool operator==(const volume<T>& a, const volume<T>& b) {
  return a.size() == b.size() && a.voxels == b.voxels;
}
template <typename T>
inline bool operator!=(const volume<T>& a, const volume<T>& b) {
  return a.size() != b.size() || a.voxels != b.voxels;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// VOLUME SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup volume
inline float lookup_volume(
    const volume<float>& vol, const vec3i& ijk, bool as_linear) {
  return vol[ijk];
}

// Evaluates a color image at a point `uv`.
inline float eval_volume(const volume<float>& vol, const vec3f& uvw,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (vol.empty()) return 0;

  // get coordinates normalized for tiling
  auto s = clamp((uvw.x + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.size().x;
  auto t = clamp((uvw.y + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.size().y;
  auto r = clamp((uvw.z + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.size().z;

  // get image coordinates and residuals
  auto i  = clamp((int)s, 0, vol.size().x - 1);
  auto j  = clamp((int)t, 0, vol.size().y - 1);
  auto k  = clamp((int)r, 0, vol.size().z - 1);
  auto ii = (i + 1) % vol.size().x, jj = (j + 1) % vol.size().y,
       kk = (k + 1) % vol.size().z;
  auto u = s - i, v = t - j, w = r - k;

  // nearest-neighbor interpolation
  if (no_interpolation) {
    i = u < 0.5 ? i : min(i + 1, vol.size().x - 1);
    j = v < 0.5 ? j : min(j + 1, vol.size().y - 1);
    k = w < 0.5 ? k : min(k + 1, vol.size().z - 1);
    return lookup_volume(vol, {i, j, k}, ldr_as_linear);
  }

  // trilinear interpolation
  return lookup_volume(vol, {i, j, k}, ldr_as_linear) * (1 - u) * (1 - v) *
             (1 - w) +
         lookup_volume(vol, {ii, j, k}, ldr_as_linear) * u * (1 - v) * (1 - w) +
         lookup_volume(vol, {i, jj, k}, ldr_as_linear) * (1 - u) * v * (1 - w) +
         lookup_volume(vol, {i, j, kk}, ldr_as_linear) * (1 - u) * (1 - v) * w +
         lookup_volume(vol, {i, jj, kk}, ldr_as_linear) * (1 - u) * v * w +
         lookup_volume(vol, {ii, j, kk}, ldr_as_linear) * u * (1 - v) * w +
         lookup_volume(vol, {ii, jj, k}, ldr_as_linear) * u * v * (1 - w) +
         lookup_volume(vol, {ii, jj, kk}, ldr_as_linear) * u * v * w;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion between flots and bytes
inline vec4b float_to_byte(const vec4f& a) {
  return {(byte)clamp(int(a.x * 256), 0, 255),
      (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255),
      (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
  return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Luminance
inline float luminance(const vec3f& a) {
  return (0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z);
}

// sRGB non-linear curve
inline float srgb_to_rgb(float srgb) {
  return (srgb <= 0.04045) ? srgb / 12.92f
                           : pow((srgb + 0.055f) / (1.0f + 0.055f), 2.4f);
}
inline float rgb_to_srgb(float rgb) {
  return (rgb <= 0.0031308f) ? 12.92f * rgb
                             : (1 + 0.055f) * pow(rgb, 1 / 2.4f) - 0.055f;
}
inline vec3f srgb_to_rgb(const vec3f& srgb) {
  return {srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z)};
}
inline vec4f srgb_to_rgb(const vec4f& srgb) {
  return {
      srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z), srgb.w};
}
inline vec3f rgb_to_srgb(const vec3f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z)};
}
inline vec4f rgb_to_srgb(const vec4f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z), rgb.w};
}

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f contrast(const vec3f& rgb, float contrast, float grey) {
  return max(zero3f, grey + (rgb - grey) * (contrast * 2));
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f logcontrast(const vec3f& rgb, float logcontrast, float grey) {
  auto epsilon  = (float)0.0001;
  auto log_grey = log2(grey);
  auto log_ldr  = log2(rgb + epsilon);
  auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
  return max(zero3f, exp2(adjusted) - epsilon);
}
// Apply saturation.
inline vec3f saturate(
    const vec3f& rgb, float saturation, const vec3f& weights) {
  auto grey = dot(weights, rgb);
  return max(zero3f, grey + (rgb - grey) * (saturation * 2));
}

// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat = mat3f{
      {0.4124, 0.2126, 0.0193},
      {0.3576, 0.7152, 0.1192},
      {0.1805, 0.0722, 0.9504},
  };
  return mat * rgb;
}
inline vec3f xyz_to_rgb(const vec3f& xyz) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat = mat3f{
      {+3.2406, -0.9689, +0.0557},
      {-1.5372, +1.8758, -0.2040},
      {-0.4986, +0.0415, +1.0570},
  };
  return mat * xyz;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
  if (xyz == zero3f) return zero3f;
  return {
      xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z), xyz.y};
}
inline vec3f xyY_to_xyz(const vec3f& xyY) {
  if (xyY.y == 0) return zero3f;
  return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}

}  // namespace yocto

#endif
