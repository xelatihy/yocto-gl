//
// # Yocto/Image : Image utilities
//
// Yocto/Image is a collection of image utilities useful when writing rendering
// algorithms. These include a simple image data structure, color conversion
// utilities and tone mapping, loading and saving functionality, and image
// resizing.
// Yocto/Image is implemented in `yocto_image.h` and `yocto_image.cpp`, and
// depends on `stb_image_resize.h` for image resizing.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef _YOCTO_IMAGE_H_
#define _YOCTO_IMAGE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <string>
#include <utility>
#include <vector>

#include "yocto_color.h"
#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image data as array of float or byte pixels. Images can be stored in linear
// or non linear color space.
template <typename T>
struct image_t {
  // iterators and references
  using iterator        = typename vector<T>::iterator;
  using const_iterator  = typename vector<T>::const_iterator;
  using reference       = typename vector<T>::reference;
  using const_reference = typename vector<T>::const_reference;

  // image api similar to vector
  image_t();
  image_t(const vec2i& size, const T& value = T{});
  image_t(const vec2i& size, const T* values);

  // size
  bool  empty() const;
  vec2i size() const;

  // iterator access
  iterator       begin();
  iterator       end();
  const_iterator begin() const;
  const_iterator end() const;

  // pixel access
  T&       operator[](vec2i ij);
  const T& operator[](vec2i ij) const;

  // data access
  T*               data();
  const T*         data() const;
  vector<T>&       pixels();
  const vector<T>& pixels() const;

 private:
  vec2i     _size = {0, 0};
  vector<T> _data = {};
};

// equality
template <typename T>
inline bool operator==(const image_t<T>& a, const image_t<T>& b);
template <typename T>
inline bool operator!=(const image_t<T>& a, const image_t<T>& b);

// Image data
template <typename T>
inline float image_aspect(const image_t<T>& image);

// Evaluates an image at a point `uv`.
template <typename T>
inline T eval_image(const image_t<T>& image, const vec2f& uv,
    bool no_interpolation = false, bool clamp_to_edge = false);

// Make an image from a function that goes from [0,1]x[0,1] to T.
template <typename Func, typename T = std::invoke_result_t<Func, vec2f>>
inline image_t<T> make_image(const vec2i& size, Func&& func);
template <typename Func, typename T = std::invoke_result_t<Func, vec2f>>
inline image_t<T> make_image(const vec2i& size, int samples, Func&& func);

// Convolve an image with a kernel
template <typename T>
inline image_t<T> convolve_image(
    const image_t<T>& img, const image_t<float>& kernel);
// Convolve an image with a separable kernel
template <typename T>
inline image_t<T> convolve_image(
    const image_t<T>& img, const vector<float>& kernel);
// Gaussian kernel
inline vector<float>  make_gaussian_kernel1d(float sigma);
inline image_t<float> make_gaussian_kernel2d(float sigma);
// Sharpening kernel
inline vector<float>  make_sharpening_kernel1d(float sigma, float scale = 1);
inline image_t<float> make_sharpening_kernel2d(float sigma, float scale = 1);

// Get/Set region
template <typename T>
inline image_t<T> get_region(
    const image_t<T>& source, const vec2i& offset, const vec2i& size);
template <typename T>
inline void get_region(image_t<T>& region, const image_t<T>& source,
    const vec2i& offset, const vec2i& size);
template <typename T>
inline void set_region(
    image_t<T>& destination, const image_t<T>& source, const vec2i& offset);

// Convenience functions
template <typename T>
inline T sum(const image_t<T>& img);
template <typename T>
inline T sum(const vector<T>& img);

// Image reconstruction
template <typename T, typename Func>
inline T reconstruct_image(const image_t<T>& img, const vec2f& uv,
    Func&& kernel, int kradius = 2, bool clamp_to_edge = false);

// Reconstruction kernels
inline float box_filter(float x);
inline float hat_filter(float x);
inline float bspline_filter(float x);
inline float catmullrom_filter(float x);
inline float mitchell_filter(float x);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
image_t<vec4f> byte_to_float(const image_t<vec4b>& bt);
image_t<vec4b> float_to_byte(const image_t<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image_t<vec4f> srgb_to_rgb(const image_t<vec4f>& srgb);
image_t<vec4f> rgb_to_srgb(const image_t<vec4f>& rgb);
image_t<vec4f> srgbb_to_rgb(const image_t<vec4b>& srgb);
image_t<vec4b> rgb_to_srgbb(const image_t<vec4f>& rgb);

// Apply exposure and filmic tone mapping
image_t<vec4f> tonemap_image(const image_t<vec4f>& hdr, float exposure = 0,
    bool filmic = false, bool srgb = true);
image_t<vec4b> tonemapb_image(const image_t<vec4f>& hdr, float exposure = 0,
    bool filmic = false, bool srgb = true);
// fast tone map for UI
void tonemap_image(image_t<vec4f>& ldr, const image_t<vec4f>& hdr,
    float exposure = 0, bool filmic = false, bool srgb = true);

// Apply exposure and filmic tone mapping
image_t<vec4f> colorgrade_image(
    const image_t<vec4f>& img, bool linear, const colorgrade_params& params);
// compute white balance
vec3f compute_white_balance(const image_t<vec4f>& img);
// fast color grade for UI
void colorgrade_image(image_t<vec4f>& graded, const image_t<vec4f>& img,
    bool linear, const colorgrade_params& params);

// image compositing
image_t<vec4f> composite_image(
    const image_t<vec4f>& foreground, const image_t<vec4f>& background);

// removes alpha
image_t<vec4f> remove_alpha(const image_t<vec4f>& img);

// turns alpha into a gray scale image
image_t<vec4f> alpha_to_gray(const image_t<vec4f>& img);

// image difference
image_t<vec4f> image_difference(
    const image_t<vec4f>& a, const image_t<vec4f>& b, bool display);

// image resizing
image_t<vec4f> resize_image(const image_t<vec4f>& img, vec2i resize);
image_t<vec4b> resize_image(const image_t<vec4b>& img, vec2i resize);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a grid image.
image_t<vec4f> make_grid(vec2i size, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a checker image.
image_t<vec4f> make_checker(vec2i size, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a bump map.
image_t<vec4f> make_bumps(vec2i size, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a ramp
image_t<vec4f> make_ramp(vec2i size, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a gamma ramp.
image_t<vec4f> make_gammaramp(vec2i size, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a uv ramp
image_t<vec4f> make_uvramp(vec2i size, float scale = 1);
// Make a uv grid
image_t<vec4f> make_uvgrid(vec2i size, float scale = 1, bool colored = true);
// Make blackbody ramp.
image_t<vec4f> make_blackbodyramp(
    vec2i size, float scale = 1, float from = 1000, float to = 12000);
// Make color map ramp.
image_t<vec4f> make_colormapramp(vec2i size, float scale = 1);
// Make a noise image. Noise parameters: lacunarity, gain, octaves, offset.
image_t<vec4f> make_noisemap(vec2i size, float scale = 1,
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});
image_t<vec4f> make_fbmmap(vec2i size, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
image_t<vec4f> make_turbulencemap(vec2i size, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
image_t<vec4f> make_ridgemap(vec2i size, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
image_t<vec4f> make_sunsky(vec2i size, float sun_angle, float turbidity = 3,
    bool has_sun = false, float sun_intensity = 1, float sun_radius = 1,
    const vec3f& ground_albedo = {0.2f, 0.2f, 0.2f});
// Make an image of multiple lights.
image_t<vec4f> make_lights(vec2i size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
image_t<vec4f> bump_to_normal(const image_t<vec4f>& image, float scale = 1);

// Add a border to an image
image_t<vec4f> add_border(
    const image_t<vec4f>& img, float width, const vec4f& color = {0, 0, 0, 1});

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
void byte_to_float(vector<vec4f>& fl, const vector<vec4b>& bt);
void float_to_byte(vector<vec4b>& bt, const vector<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4f>& srgb);
void rgb_to_srgb(vector<vec4f>& srgb, const vector<vec4f>& rgb);
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4b>& srgb);
void rgb_to_srgb(vector<vec4b>& srgb, const vector<vec4f>& rgb);

// Apply tone mapping
void tonemap_image(vector<vec4f>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic = false, bool srgb = true);
void tonemap_image(vector<vec4b>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic = false, bool srgb = true);

// Apply tone mapping using multithreading for speed
void tonemap_image_mt(vector<vec4f>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic = false, bool srgb = true);
void tonemap_image_mt(vector<vec4b>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic = false, bool srgb = true);

// Color grade a linear or srgb image to an srgb image.
// Uses multithreading for speed.
void colorgrade_image_mt(vector<vec4f>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params);
void colorgrade_image_mt(vector<vec4b>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params);

// determine white balance colors
vec3f compute_white_balance(const vector<vec4f>& img);

// Resize an image.
void resize_image(vector<vec4f>& res, const vector<vec4f>& img, int width,
    int height, int res_width, int res_height);
void resize_image(vector<vec4b>& res, const vector<vec4b>& img, int width,
    int height, int res_width, int res_height);

// Compute the difference between two images
void image_difference(vector<vec4f>& diff, const vector<vec4f>& a,
    const vector<vec4f>& b, bool disply_diff);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a grid image.
void make_grid(vector<vec4f>& pixels, int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a checker image.
void make_checker(vector<vec4f>& pixels, int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a bump map.
void make_bumps(vector<vec4f>& pixels, int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a ramp
void make_ramp(vector<vec4f>& pixels, int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a gamma ramp.
void make_gammaramp(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a uv ramp
void make_uvramp(vector<vec4f>& pixels, int width, int height, float scale = 1);
// Make a uv grid
void make_uvgrid(vector<vec4f>& pixels, int width, int height, float scale = 1,
    bool colored = true);
// Make blackbody ramp.
void make_blackbodyramp(vector<vec4f>& pixels, int width, int height,
    float scale = 1, float from = 1000, float to = 12000);
// Make color map ramp.
void make_colormapramp(
    vector<vec4f>& pixels, int width, int height, float scale = 1);
// Make a noise image. Noise parameters: lacunarity, gain, octaves, offset.
void make_noisemap(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
void make_fbmmap(vector<vec4f>& pixels, int width, int height, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
void make_turbulencemap(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& noise = {2, 0.5, 8, 1},
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});
void make_ridgemap(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& noise = {2, 0.5, 8, 1},
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});

// Make a random image.
void make_randpoints(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
void make_randlines(vector<vec4f>& pixels, int width, int height,
    float scale = 1, const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
void make_sunsky(vector<vec4f>& pixels, int width, int height, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
    float sun_radius = 1, const vec3f& ground_albedo = {0.2f, 0.2f, 0.2f});
// Make an image of multiple lights.
void make_lights(vector<vec4f>& pixels, int width, int height,
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
    float lwidth = pif / 16, float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
void bump_to_normal(vector<vec4f>& normal, const vector<vec4f>& bump, int width,
    int height, float scale = 1);

// Add a border to an image
void add_border(vector<vec4f>& pixels, const vector<vec4f>& source, int width,
    int height, float thickness, const vec4f& color = {0, 0, 0, 1});

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// image api similar to vector
template <typename T>
inline image_t<T>::image_t() : _size{0, 0}, _data{} {}
template <typename T>
inline image_t<T>::image_t(const vec2i& size, const T& value)
    : _size{size}, _data((size_t)size.x * (size_t)size.y, value) {}
template <typename T>
inline image_t<T>::image_t(const vec2i& size, const T* values)
    : _size{size}, _data{values, values + size.x * size.y} {}

// size
template <typename T>
inline bool image_t<T>::empty() const {
  return _size == zero2i;
}
template <typename T>
inline vec2i image_t<T>::size() const {
  return _size;
}

// iterator access
template <typename T>
inline typename image_t<T>::iterator image_t<T>::begin() {
  return _data.begin();
}
template <typename T>
inline typename image_t<T>::iterator image_t<T>::end() {
  return _data.end();
}
template <typename T>
inline typename image_t<T>::const_iterator image_t<T>::begin() const {
  return _data.begin();
}
template <typename T>
inline typename image_t<T>::const_iterator image_t<T>::end() const {
  return _data.end();
}

// pixel access
template <typename T>
inline T& image_t<T>::operator[](vec2i ij) {
  return _data[ij.y * _size.x + ij.x];
}
template <typename T>
inline const T& image_t<T>::operator[](vec2i ij) const {
  return _data[ij.y * _size.x + ij.x];
}

// data access
template <typename T>
inline T* image_t<T>::data() {
  return _data.data();
}
template <typename T>
inline const T* image_t<T>::data() const {
  return _data.data();
}
template <typename T>
inline vector<T>& image_t<T>::pixels() {
  return _data;
}
template <typename T>
inline const vector<T>& image_t<T>::pixels() const {
  return _data;
}

// equality
template <typename T>
inline bool operator==(const image_t<T>& a, const image_t<T>& b) {
  return a.size() == b.size() && a.pixels() == b.pixels();
}
template <typename T>
inline bool operator!=(const image_t<T>& a, const image_t<T>& b) {
  return a.size() != b.size() || a.pixels() != b.pixels();
}

// Image data
template <typename T>
inline float image_aspect(const image_t<T>& image) {
  if (image.empty()) return 0;
  return float{image.size().x} / float{image.size().y};
}

// Evaluates an image at a point `uv`.
template <typename T>
inline T eval_image(const image_t<T>& image, const vec2f& uv,
    bool no_interpolation, bool clamp_to_edge) {
  if (image.empty()) return T{};

  // get image width/height
  auto size = image.size();

  // handle interpolation
  if (no_interpolation) {
    // get coordinates normalized for tiling
    auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * (vec2f)size;

    // lookup image
    auto ij = clamp((vec2i)st, zero2i, size - 1);
    return image[ij];
  } else {
    // get coordinates normalized for tiling
    auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * (vec2f)size;
    st -= 0.5f;  // shift pixel centers to the pixel corners

    // get image coordinates
    auto ij   = clamp((vec2i)st, zero2i, size - 1);
    auto i1j  = clamp_to_edge ? min(ij + vec2i{1, 0}, size - 1)
                              : ((ij + vec2i{1, 0}) % size);
    auto ij1  = clamp_to_edge ? min(ij + vec2i{0, 1}, size - 1)
                              : ((ij + vec2i{0, 1}) % size);
    auto i1j1 = clamp_to_edge ? min(ij + vec2i{1, 1}, size - 1)
                              : ((ij + vec2i{1, 1}) % size);

    // interpolate
    auto w = st - (vec2f)ij;
    return image[ij] * (1 - w.x) * (1 - w.y) + image[ij1] * (1 - w.x) * w.y +
           image[i1j] * w.x * (1 - w.y) + image[i1j1] * w.x * w.y;
  }
}

// Make an image from a function that goes from [0,1]x[0,1] to T.
template <typename Func, typename T>
inline image_t<T> make_image(const vec2i& size, Func&& func) {
  auto img = image_t<T>(size);
  for (auto ij : range(size)) {
    img[ij] = func(((vec2f)ij + 0.5f) / (vec2f)size);
  }
  return img;
}

// Make an image from a function that goes from [0,1]x[0,1] to T.
template <typename Func, typename T>
inline image_t<T> make_image(const vec2i& size, int samples, Func&& func) {
  if (samples <= 1) return make_image(size, std::forward<Func>(func));
  auto ns  = (int)round(sqrt((float)samples));
  auto img = image_t<T>(size, T{});
  for (auto ij : range(size)) {
    for (auto sij : range({ns, ns})) {
      img[ij] += func(((vec2f)ij + ((vec2f)sij + 0.5f) / ns) / (vec2f)size);
    }
    img[ij] /= ns * ns;
  }
  return img;
}

// Evaluates a curve at a point `u`.
// This is just for the book images.
template <typename T>
inline T eval_curve(const vector<T>& curve, float u,
    bool no_interpolation = false, bool clamp_to_edge = false) {
  if (curve.empty()) return T{};

  // get curve size
  auto size = (int)curve.size();

  // handle interpolation
  if (no_interpolation) {
    auto s = (clamp_to_edge ? clamp(u, 0.0f, 1.0f) : mod(u, 1.0f)) * size;
    auto i = clamp((int)s, 0, size - 1);
    return curve[i];
  } else {
    auto s = (clamp_to_edge ? clamp(u, 0.0f, 1.0f) : mod(u, 1.0f)) * size;
    s -= 0.5f;  // shift pixel centers to the pixel corners
    auto i  = clamp((int)s, 0, size - 1);
    auto i1 = clamp_to_edge ? min(i + 1, size - 1) : ((i + 1) % size);
    auto w  = s - i;
    return curve[i] * (1 - w) + curve[i1] * w;
  }
}

// Make an image from a function that goes from [0,1] to T.
// This is just for the book images.
template <typename Func, typename T = std::invoke_result_t<Func, float>>
inline vector<T> make_curve(int size, Func&& func) {
  auto curve = vector<T>(size);
  for (auto i : range(size)) {
    curve[i] = func((i + 0.5f) / size);
  }
  return curve;
}

// Convenience functions
template <typename T>
inline T sum(const image_t<T>& img) {
  auto sum = T{};
  for (auto& value : img) sum += value;
  return sum;
}
template <typename T>
inline T sum(const vector<T>& values) {
  auto sum = T{};
  for (auto& value : values) sum += value;
  return sum;
}

// Convolve an image with a kernel
template <typename T>
inline image_t<T> convolve_image(
    const image_t<T>& img, const image_t<float>& kernel) {
  auto result  = image_t<T>(img.size());
  auto kcenter = kernel.size() / 2;
  for (auto ij : range(img.size())) {
    auto sum    = T{};
    auto weight = 0.0f;
    for (auto kij : range(kernel.size())) {
      auto iijj = ij + kij - kcenter;
      if (iijj.x < 0 || iijj.x >= img.size().x || iijj.y < 0 ||
          iijj.y >= img.size().y)
        continue;
      sum += img[iijj] * kernel[kij];
      weight += kernel[kij];
    }
    result[ij] = max(sum / weight, T{});
  }
  return result;
}
// Convolve an image with a separable kernel
template <typename T>
inline image_t<T> convolve_image(
    const image_t<T>& img, const vector<float>& kernel);
// Gaussian kernel
inline vector<float> make_gaussian_kernel1d(float sigma) {
  auto width  = (int)ceil(sigma * 3);
  auto kernel = vector<float>(width);
  for (auto i : range(kernel.size()))
    kernel[i] = exp(
        -0.5f * (i - width / 2) * (i - width / 2) / (sigma * sigma));
  for (auto& value : kernel) value /= sum(kernel);
  return kernel;
}
inline image_t<float> make_gaussian_kernel2d(float sigma) {
  auto kernel1d = make_gaussian_kernel1d(sigma);
  auto kernel   = image_t<float>({(int)kernel1d.size(), (int)kernel1d.size()});
  for (auto ij : range(kernel.size()))
    kernel[ij] = kernel1d[ij.x] * kernel1d[ij.y];
  for (auto& value : kernel) value /= sum(kernel);
  return kernel;
}
// Sharpening kernel
inline vector<float> make_sharpening_kernel1d(float sigma, float scale) {
  auto kernel = make_gaussian_kernel1d(sigma);
  for (auto& value : kernel) value = -scale * value;
  kernel[kernel.size() / 2] += 1 + scale;
  return kernel;
}
// Sharpening kernel
inline image_t<float> make_sharpening_kernel2d(float sigma, float scale) {
  auto kernel = make_gaussian_kernel2d(sigma);
  for (auto& value : kernel) value = -scale * value;
  kernel[kernel.size() / 2] += 1 + scale;
  return kernel;
}

// Get/Set region
template <typename T>
inline image_t<T> get_region(
    const image_t<T>& source, const vec2i& offset, const vec2i& size) {
  auto region = image_t<T>(size);
  for (auto ij : range(region.size())) region[ij] = source[ij + offset];
  return region;
}
template <typename T>
inline void get_region(image_t<T>& region, const image_t<T>& source,
    const vec2i& offset, const vec2i& size) {
  if (region.size() != size) region = image_t<T>(size);
  for (auto ij : range(region.size())) region[ij] = source[ij + offset];
}
template <typename T>
inline void set_region(
    image_t<T>& destination, const image_t<T>& source, const vec2i& offset) {
  for (auto ij : range(source.size())) destination[ij + offset] = source[ij];
}

inline float box_filter(float x) { return x < 0.5f && x >= -0.5f ? 1 : 0; }
inline float hat_filter(float x) {
  auto xp = abs(x);
  return xp < 1 ? 1 - xp : 0;
}
inline float bspline_filter(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return (4 + xp * xp * (3 * xp - 6)) / 6;
  else if (xp < 2.0f)
    return (8 + xp * (-12 + xp * (6 - xp))) / 6;
  else
    return 0;
}
inline float catmullrom_filter(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return 1 - xp * xp * (2.5f - 1.5f * xp);
  else if (xp < 2.0f)
    return 2 - xp * (4 + xp * (0.5f * xp - 2.5f));
  else
    return 0;
}
inline float mitchell_filter(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return (16 + xp * xp * (21 * xp - 36)) / 18;
  else if (xp < 2.0f)
    return (32 + xp * (-60 + xp * (36 - 7 * xp))) / 18;
  else
    return 0;
}

// Image reconstruction
template <typename T, typename Func>
inline T reconstruct_image(const image_t<T>& img, const vec2f& uv,
    Func&& kernel, int kradius, bool clamp_to_edge) {
  auto x    = uv * (vec2f)img.size();
  auto ij   = (vec2i)(x + 0.5f);
  auto size = img.size();
  auto sum  = T{};
  for (auto j : range(ij.y - kradius, ij.y + kradius)) {
    for (auto i : range(ij.x - kradius, ij.x + kradius)) {
      auto xi = vec2f{(float)i, (float)j} + 0.5f;  // shift centers
      auto w  = kernel(x - xi);
      sum += img[clamp(vec2i{i, j}, vec2i{0, 0}, size - 1)] * w;
    }
  }
  return sum;
}
template <typename T, typename Func>
inline T reconstruct_curve(const vector<T>& curve, float u, Func&& kernel,
    int kradius, bool clamp_to_edge) {
  auto x    = u * (float)curve.size();
  auto size = (int)curve.size();
  auto sum  = 0.0f;
  for (auto i : range(-2, size + 2)) {
    auto xi = i + 0.5f;  // shift centers
    auto w  = kernel(x - xi);
    sum += curve[clamp(i, 0, size - 1)] * w;
  }
  return sum;
}

}  // namespace yocto

#endif
