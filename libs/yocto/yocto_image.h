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
#include <tuple>
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
using std::tuple;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image data as array of float or byte pixels. Images can be stored in linear
// or non linear color space.
template <typename T>
struct image {
  // iterators and references
  using iterator        = typename vector<T>::iterator;
  using const_iterator  = typename vector<T>::const_iterator;
  using reference       = typename vector<T>::reference;
  using const_reference = typename vector<T>::const_reference;

  // image api similar to vector
  image();
  image(vec2i size, const T& value = T{});
  image(vec2i size, const T* values);

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
inline bool operator==(const image<T>& a, const image<T>& b);
template <typename T>
inline bool operator!=(const image<T>& a, const image<T>& b);

// Image data
template <typename T>
inline float image_aspect(const image<T>& image);

// Evaluates an image at a point `uv`.
template <typename T>
inline T eval_image(const image<T>& image, vec2f uv,
    bool no_interpolation = false, bool clamp_to_edge = false);

// Make an image from a function that goes from [0,1]x[0,1] to T.
template <typename Func, typename T = std::invoke_result_t<Func, vec2f>>
inline image<T> make_image(vec2i size, Func&& func);
template <typename Func, typename T = std::invoke_result_t<Func, vec2f>>
inline image<T> make_image(vec2i size, int samples, Func&& func);

// Make an image by applying a function to reach pixel of another.
template <typename S, typename Func, typename T = std::invoke_result_t<Func, S>>
inline image<T> transform_image(const image<S>& source, Func&& func);
template <typename T, typename S, typename Func>
inline void transform_image(
    image<T>& result, const image<S>& source, Func&& func);
template <typename S1, typename S2, typename Func,
    typename T = std::invoke_result_t<Func, S1, S2>>
inline image<T> transform_images(
    const image<S1>& source1, const image<S2>& source2, Func&& func);
template <typename T, typename S1, typename S2, typename Func>
inline void transform_images(image<T>& result, const image<S1>& source1,
    const image<S2>& source2, Func&& func);

// Convolve an image with a kernel
template <typename T>
inline image<T> convolve_image(const image<T>& img, const image<float>& kernel);
// Convolve an image with a separable kernel
template <typename T>
inline image<T> convolve_image(
    const image<T>& img, const vector<float>& kernel);
// Gaussian kernel
inline vector<float> make_gaussian_kernel1d(float sigma);
inline image<float>  make_gaussian_kernel2d(float sigma);
// Sharpening kernel
inline vector<float> make_sharpening_kernel1d(float sigma, float scale = 1);
inline image<float>  make_sharpening_kernel2d(float sigma, float scale = 1);

// Get/Set region
template <typename T>
inline image<T> get_region(const image<T>& source, vec2i offset, vec2i size);
template <typename T>
inline void get_region(
    image<T>& region, const image<T>& source, vec2i offset, vec2i size);
template <typename T>
inline void set_region(
    image<T>& destination, const image<T>& source, vec2i offset);

// Convenience functions
template <typename T>
inline T sum(const image<T>& img);
template <typename T>
inline T sum(const vector<T>& img);

// Reconstruction kernels
inline float box_kernel(float u);
inline float hat_kernel(float u);
inline float bspline_kernel(float u);
inline float catmullrom_kernel(float u);
inline float mitchell_kernel(float u);
inline float box_kernel(vec2f u);
inline float hat_kernel(vec2f u);
inline float bspline_kernel(vec2f u);
inline float catmullrom_kernel(vec2f u);
inline float mitchell_kernel(vec2f u);

// Generic kernel
enum struct reconstruction_type { box, hat, bspline, catmullrom, mitchell };
inline float reconstruction_kernel(reconstruction_type ktype, float u);
inline float reconstruction_kernel(reconstruction_type ktype, vec2f u);
inline int   reconstruction_radius(reconstruction_type ktype);

// Image reconstruction
template <typename T, typename Func>
inline T reconstruct_image(const image<T>& img, vec2f uv, Func&& kernel,
    int kradius = 2, bool clamp_to_edge = false);
template <typename T>
inline T reconstruct_image(const image<T>& img, vec2f uv,
    reconstruction_type ktype, bool clamp_to_edge = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt);
image<vec4b> float_to_byte(const image<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_rgb(const image<vec4f>& srgb);
image<vec4f> rgb_to_srgb(const image<vec4f>& rgb);
image<vec4f> srgbb_to_rgb(const image<vec4b>& srgb);
image<vec4b> rgb_to_srgbb(const image<vec4f>& rgb);

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_image(const image<vec4f>& hdr, float exposure = 0,
    bool filmic = false, bool srgb = true);
image<vec4b> tonemapb_image(const image<vec4f>& hdr, float exposure = 0,
    bool filmic = false, bool srgb = true);
// fast tone map for UI
void tonemap_image(image<vec4f>& ldr, const image<vec4f>& hdr,
    float exposure = 0, bool filmic = false, bool srgb = true);

// Apply exposure and filmic tone mapping
image<vec4f> colorgrade_image(
    const image<vec4f>& img, bool linear, const colorgrade_params& params);
// compute white balance
vec3f compute_white_balance(const image<vec4f>& img);
// fast color grade for UI
void colorgrade_image(image<vec4f>& graded, const image<vec4f>& img,
    bool linear, const colorgrade_params& params);

// image compositing
image<vec4f> composite_image(
    const image<vec4f>& foreground, const image<vec4f>& background);

// removes alpha
image<vec4f> remove_alpha(const image<vec4f>& img);

// turns alpha into a gray scale image
image<vec4f> alpha_to_gray(const image<vec4f>& img);

// image difference
image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool display);

// image resizing
image<vec4f> resize_image(const image<vec4f>& img, vec2i resize);
image<vec4b> resize_image(const image<vec4b>& img, vec2i resize);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a grid image.
image<vec4f> make_grid(vec2i size, float scale = 1,
    vec4f color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    vec4f color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a checker image.
image<vec4f> make_checker(vec2i size, float scale = 1,
    vec4f color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    vec4f color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a bump map.
image<vec4f> make_bumps(vec2i size, float scale = 1,
    vec4f color0 = vec4f{0, 0, 0, 1}, vec4f color1 = vec4f{1, 1, 1, 1});
// Make a ramp
image<vec4f> make_ramp(vec2i size, float scale = 1,
    vec4f color0 = vec4f{0, 0, 0, 1}, vec4f color1 = vec4f{1, 1, 1, 1});
// Make a gamma ramp.
image<vec4f> make_gammaramp(vec2i size, float scale = 1,
    vec4f color0 = vec4f{0, 0, 0, 1}, vec4f color1 = vec4f{1, 1, 1, 1});
// Make a uv ramp
image<vec4f> make_uvramp(vec2i size, float scale = 1);
// Make a uv grid
image<vec4f> make_uvgrid(vec2i size, float scale = 1, bool colored = true);
// Make color map ramp.
image<vec4f> make_colormapramp(vec2i size, float scale = 1);

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
image<vec4f> make_sunsky(vec2i size, float sun_angle, float turbidity = 3,
    bool has_sun = false, float sun_intensity = 1, float sun_radius = 1,
    vec3f ground_albedo = {0.2f, 0.2f, 0.2f});
// Make an image of multiple lights.
image<vec4f> make_lights(vec2i size, vec3f le = {1, 1, 1}, int nlights = 4,
    float langle = pif / 4, float lwidth = pif / 16, float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
image<vec4f> bump_to_normal(const image<vec4f>& image, float scale = 1);

// Add a border to an image
image<vec4f> add_border(
    const image<vec4f>& img, float width, vec4f color = {0, 0, 0, 1});

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
inline image<T>::image() : _size{0, 0}, _data{} {}
template <typename T>
inline image<T>::image(vec2i size, const T& value)
    : _size{size}, _data((size_t)size.x * (size_t)size.y, value) {}
template <typename T>
inline image<T>::image(vec2i size, const T* values)
    : _size{size}, _data{values, values + size.x * size.y} {}

// size
template <typename T>
inline bool image<T>::empty() const {
  return _size == zero2i;
}
template <typename T>
inline vec2i image<T>::size() const {
  return _size;
}

// iterator access
template <typename T>
inline typename image<T>::iterator image<T>::begin() {
  return _data.begin();
}
template <typename T>
inline typename image<T>::iterator image<T>::end() {
  return _data.end();
}
template <typename T>
inline typename image<T>::const_iterator image<T>::begin() const {
  return _data.begin();
}
template <typename T>
inline typename image<T>::const_iterator image<T>::end() const {
  return _data.end();
}

// pixel access
template <typename T>
inline T& image<T>::operator[](vec2i ij) {
  return _data[ij.y * _size.x + ij.x];
}
template <typename T>
inline const T& image<T>::operator[](vec2i ij) const {
  return _data[ij.y * _size.x + ij.x];
}

// data access
template <typename T>
inline T* image<T>::data() {
  return _data.data();
}
template <typename T>
inline const T* image<T>::data() const {
  return _data.data();
}
template <typename T>
inline vector<T>& image<T>::pixels() {
  return _data;
}
template <typename T>
inline const vector<T>& image<T>::pixels() const {
  return _data;
}

// equality
template <typename T>
inline bool operator==(const image<T>& a, const image<T>& b) {
  return a.size() == b.size() && a.pixels() == b.pixels();
}
template <typename T>
inline bool operator!=(const image<T>& a, const image<T>& b) {
  return a.size() != b.size() || a.pixels() != b.pixels();
}

// Image data
template <typename T>
inline float image_aspect(const image<T>& image) {
  if (image.empty()) return 0;
  return float{image.size().x} / float{image.size().y};
}

// Evaluates an image at a point `uv`.
template <typename T>
inline T eval_image(const image<T>& image, vec2f uv, bool no_interpolation,
    bool clamp_to_edge) {
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
inline image<T> make_image(vec2i size, Func&& func) {
  auto img = image<T>(size);
  for (auto ij : range(size)) {
    img[ij] = func(((vec2f)ij + 0.5f) / (vec2f)size);
  }
  return img;
}

// Make an image from a function that goes from [0,1]x[0,1] to T.
template <typename Func, typename T>
inline image<T> make_image(vec2i size, int samples, Func&& func) {
  if (samples <= 1) return make_image(size, std::forward<Func>(func));
  auto ns  = (int)round(sqrt((float)samples));
  auto img = image<T>(size, T{});
  for (auto ij : range(size)) {
    for (auto sij : range({ns, ns})) {
      img[ij] += func(((vec2f)ij + ((vec2f)sij + 0.5f) / ns) / (vec2f)size);
    }
    img[ij] /= ns * ns;
  }
  return img;
}

// Make an image by applying a function to reach pixel of another.
template <typename S, typename Func, typename T>
inline image<T> transform_image(const image<S>& source, Func&& func) {
  auto result = image<T>(source.size());
  for (auto&& [res, src] : zip(result, source)) res = func(src);
  return result;
}
template <typename T, typename S, typename Func>
inline void transform_image(
    image<T>& result, const image<S>& source, Func&& func) {
  if (result.size() != source.size())
    throw std::runtime_error("image size mismatch");
  for (auto&& [res, src] : zip(result, source)) res = func(src);
}
template <typename S1, typename S2, typename Func, typename T>
inline image<T> transform_images(
    const image<S1>& source1, const image<S2>& source2, Func&& func) {
  if (source1.size() != source2.size())
    throw std::runtime_error("image size mismatch");
  auto result = image<T>(source1.size());
  for (auto&& [res, src1, src2] : zip(result, source1, source2))
    res = func(src1, src2);
  return result;
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
inline T sum(const image<T>& img) {
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
inline image<T> convolve_image(
    const image<T>& img, const image<float>& kernel) {
  auto result  = image<T>(img.size());
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
inline image<T> convolve_image(
    const image<T>& img, const vector<float>& kernel);
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
inline image<float> make_gaussian_kernel2d(float sigma) {
  auto kernel1d = make_gaussian_kernel1d(sigma);
  auto kernel   = image<float>({(int)kernel1d.size(), (int)kernel1d.size()});
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
inline image<float> make_sharpening_kernel2d(float sigma, float scale) {
  auto kernel = make_gaussian_kernel2d(sigma);
  for (auto& value : kernel) value = -scale * value;
  kernel[kernel.size() / 2] += 1 + scale;
  return kernel;
}

// Get/Set region
template <typename T>
inline image<T> get_region(const image<T>& source, vec2i offset, vec2i size) {
  auto region = image<T>(size);
  for (auto ij : range(region.size())) region[ij] = source[ij + offset];
  return region;
}
template <typename T>
inline void get_region(
    image<T>& region, const image<T>& source, vec2i offset, vec2i size) {
  if (region.size() != size) region = image<T>(size);
  for (auto ij : range(region.size())) region[ij] = source[ij + offset];
}
template <typename T>
inline void set_region(
    image<T>& destination, const image<T>& source, vec2i offset) {
  for (auto ij : range(source.size())) destination[ij + offset] = source[ij];
}

// Kernels
inline float box_kernel(float x) { return x < 0.5f && x >= -0.5f ? 1 : 0; }
inline float hat_kernel(float x) {
  auto xp = abs(x);
  return xp < 1 ? 1 - xp : 0;
}
inline float bspline_kernel(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return (4 + xp * xp * (3 * xp - 6)) / 6;
  else if (xp < 2.0f)
    return (8 + xp * (-12 + xp * (6 - xp))) / 6;
  else
    return 0;
}
inline float catmullrom_kernel(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return 1 - xp * xp * (2.5f - 1.5f * xp);
  else if (xp < 2.0f)
    return 2 - xp * (4 + xp * (0.5f * xp - 2.5f));
  else
    return 0;
}
inline float mitchell_kernel(float x) {
  auto xp = abs(x);

  if (xp < 1.0f)
    return (16 + xp * xp * (21 * xp - 36)) / 18;
  else if (xp < 2.0f)
    return (32 + xp * (-60 + xp * (36 - 7 * xp))) / 18;
  else
    return 0;
}

// Kernels
inline float box_kernel(vec2f x) { return box_kernel(x.x) * box_kernel(x.y); }
inline float hat_kernel(vec2f x) { return hat_kernel(x.x) * hat_kernel(x.y); }
inline float bspline_kernel(vec2f x) {
  return bspline_kernel(x.x) * bspline_kernel(x.y);
}
inline float catmullrom_kernel(vec2f x) {
  return catmullrom_kernel(x.x) * catmullrom_kernel(x.y);
}
inline float mitchell_kernel(vec2f x) {
  return mitchell_kernel(x.x) * mitchell_kernel(x.y);
}

// Generic kernel
inline float reconstruction_kernel(reconstruction_type ktype, float u) {
  switch (ktype) {
    case reconstruction_type::box: return box_kernel(u);
    case reconstruction_type::hat: return hat_kernel(u);
    case reconstruction_type::bspline: return bspline_kernel(u);
    case reconstruction_type::catmullrom: return catmullrom_kernel(u);
    case reconstruction_type::mitchell: return mitchell_kernel(u);
    default: return 0;
  }
}
inline float reconstruction_kernel(reconstruction_type ktype, vec2f uv) {
  switch (ktype) {
    case reconstruction_type::box: return box_kernel(uv);
    case reconstruction_type::hat: return hat_kernel(uv);
    case reconstruction_type::bspline: return bspline_kernel(uv);
    case reconstruction_type::catmullrom: return catmullrom_kernel(uv);
    case reconstruction_type::mitchell: return mitchell_kernel(uv);
    default: return 0;
  }
}
inline int reconstruction_radius(reconstruction_type ktype) {
  switch (ktype) {
    case reconstruction_type::box: return 1;
    case reconstruction_type::hat: return 1;
    case reconstruction_type::bspline: return 2;
    case reconstruction_type::catmullrom: return 2;
    case reconstruction_type::mitchell: return 2;
    default: return 0;
  }
}

// Image reconstruction
template <typename T, typename Func>
inline T reconstruct_image(const image<T>& img, vec2f uv, Func&& kernel,
    int kradius, bool clamp_to_edge) {
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
template <typename T>
inline T reconstruct_image(const image<T>& img, vec2f uv,
    reconstruction_type ktype, bool clamp_to_edge) {
  switch (ktype) {
    case reconstruction_type::box:
      return reconstruct_image(
          img, uv, [](vec2f uv) { return box_kernel(uv); }, 1, clamp_to_edge);
    case reconstruction_type::hat:
      return reconstruct_image(
          img, uv, [](vec2f uv) { return hat_kernel(uv); }, 1, clamp_to_edge);
    case reconstruction_type::bspline:
      return reconstruct_image(
          img, uv, [](vec2f uv) { return bspline_kernel(uv); }, 2,
          clamp_to_edge);
    case reconstruction_type::catmullrom:
      return reconstruct_image(
          img, uv, [](vec2f uv) { return catmullrom_kernel(uv); }, 2,
          clamp_to_edge);
    case reconstruction_type::mitchell:
      return reconstruct_image(
          img, uv, [](vec2f uv) { return mitchell_kernel(uv); }, 2,
          clamp_to_edge);
    default: return T{};
  }
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
template <typename T>
inline T reconstruct_curve(const vector<T>& curve, float u,
    reconstruction_type ktype, bool clamp_to_edge) {
  return reconstruct_curve(
      curve, u, [ktype](float u) { return reconstruction_kernel(ktype, u); },
      reconstruction_radius(ktype), clamp_to_edge);
}

}  // namespace yocto

#endif
