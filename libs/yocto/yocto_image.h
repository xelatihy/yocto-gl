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

#ifndef YOCTO_IMAGE_H_
#define YOCTO_IMAGE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <stdexcept>

#include "yocto_color.h"
#include "yocto_math.h"
#include "yocto_ndarray.h"
#include "yocto_noise.h"

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image alias
template <typename T>
using image = array2d<T>;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Aspect ratio
inline float image_aspect(const vec2i& extents) {
  return (float)extents[0] / (float)extents[1];
}

// Conversion from/to floats.
inline array2d<vec3f> byte_to_float(const array2d<vec3b>& bt) {
  return fmap(bt, [](const vec3b& a) { return byte_to_float(a); });
}
inline array2d<vec4f> byte_to_float(const array2d<vec4b>& bt) {
  return fmap(bt, [](const vec4b& a) { return byte_to_float(a); });
}
inline array2d<vec3b> float_to_byte(const array2d<vec3f>& fl) {
  return fmap(fl, [](const vec3f& a) { return float_to_byte(a); });
}
inline array2d<vec4b> float_to_byte(const array2d<vec4f>& fl) {
  return fmap(fl, [](const vec4f& a) { return float_to_byte(a); });
}
inline array2d<float> byte_to_float(const array2d<byte>& bt) {
  return fmap(bt, [](byte a) { return byte_to_float(a); });
}
inline array2d<byte> float_to_byte(const array2d<float>& fl) {
  return fmap(fl, [](float a) { return float_to_byte(a); });
}

// Conversion between linear and gamma-encoded images.
inline array2d<vec3f> srgb_to_rgb(const array2d<vec3f>& srgb) {
  return fmap(srgb, [](const vec3f& a) { return srgb_to_rgb(a); });
}
inline array2d<vec4f> srgb_to_rgb(const array2d<vec4f>& srgb) {
  return fmap(srgb, [](const vec4f& a) { return srgb_to_rgb(a); });
}
inline array2d<vec3f> rgb_to_srgb(const array2d<vec3f>& rgb) {
  return fmap(rgb, [](const vec3f& a) { return rgb_to_srgb(a); });
}
inline array2d<vec4f> rgb_to_srgb(const array2d<vec4f>& rgb) {
  return fmap(rgb, [](const vec4f& a) { return rgb_to_srgb(a); });
}
inline void srgb_to_rgb(array2d<vec3f>& rgb, const array2d<vec3f>& srgb) {
  return fmap(rgb, srgb, [](const vec3f& a) { return srgb_to_rgb(a); });
}
inline void srgb_to_rgb(array2d<vec4f>& rgb, const array2d<vec4f>& srgb) {
  return fmap(rgb, srgb, [](const vec4f& a) { return srgb_to_rgb(a); });
}
inline void rgb_to_srgb(array2d<vec3f>& srgb, const array2d<vec3f>& rgb) {
  return fmap(srgb, rgb, [](const vec3f& a) { return rgb_to_srgb(a); });
}
inline void rgb_to_srgb(array2d<vec4f>& srgb, const array2d<vec4f>& rgb) {
  return fmap(srgb, rgb, [](const vec4f& a) { return rgb_to_srgb(a); });
}

// Conversion between linear and gamma-encoded images.
inline array2d<vec3f> srgbb_to_rgb(const array2d<vec3b>& srgb) {
  return fmap(srgb, [](const vec3b& a) { return srgbb_to_rgb(a); });
}
inline array2d<vec4f> srgbb_to_rgb(const array2d<vec4b>& srgb) {
  return fmap(srgb, [](const vec4b& a) { return srgbb_to_rgb(a); });
}
inline array2d<vec3b> rgb_to_srgbb(const array2d<vec3f>& rgb) {
  return fmap(rgb, [](const vec3f& a) { return rgb_to_srgbb(a); });
}
inline array2d<vec4b> rgb_to_srgbb(const array2d<vec4f>& rgb) {
  return fmap(rgb, [](const vec4f& a) { return rgb_to_srgbb(a); });
}
inline void srgbb_to_rgb(array2d<vec3f>& rgb, const array2d<vec3b>& srgb) {
  fmap(srgb, [](const vec3b& a) { return srgbb_to_rgb(a); });
}
inline void srgbb_to_rgb(array2d<vec4f>& rgb, const array2d<vec4b>& srgb) {
  fmap(srgb, [](const vec4b& a) { return srgbb_to_rgb(a); });
}
inline void rgb_to_srgbb(array2d<vec3b>& srgb, const array2d<vec3f>& rgb) {
  fmap(rgb, [](const vec3f& a) { return rgb_to_srgbb(a); });
}
inline void rgb_to_srgbb(array2d<vec4b>& srgb, const array2d<vec4f>& rgb) {
  fmap(rgb, [](const vec4f& a) { return rgb_to_srgbb(a); });
}

// Lookup pixel for evaluation
template <typename T>
inline T lookup_image(
    const array2d<T>& image, const vec2i& ij, bool as_linear = false) {
  if constexpr (!(std::is_same_v<T, byte> || std::is_same_v<T, vec3b> ||
                    std::is_same_v<T, vec4b>)) {
    return as_linear ? srgb_to_rgb(image[ij]) : image[ij];
  } else {
    return as_linear ? srgbb_to_rgb(image[ij]) : byte_to_float(image[ij]);
  }
}

// Evaluates an image at a point `uv`.
template <typename T>
inline T eval_image(const array2d<T>& image, const vec2f& uv,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  if (image.empty()) return T{0};

  // get image width/height
  auto size = image.extents();

  // get coordinates normalized for tiling
  auto st = (clamp_to_edge ? clamp(uv, 0.0f, 1.0f) : mod(uv, 1.0f)) * size;

  // handle interpolation
  if (no_interpolation) {
    auto ij = clamp((vec2i)st, 0, size - 1);
    return lookup_image(image, ij, as_linear);
  } else {
    auto ij   = clamp((vec2i)st, 0, size - 1);
    auto i1j  = (ij + vec2i{1, 0}) % size;
    auto ij1  = (ij + vec2i{0, 1}) % size;
    auto i1j1 = (ij + vec2i{1, 1}) % size;
    auto w    = st - ij;
    return lookup_image(image, ij, as_linear) * (1 - w.x) * (1 - w.y) +
           lookup_image(image, ij1, as_linear) * (1 - w.x) * w.y +
           lookup_image(image, i1j, as_linear) * w.x * (1 - w.y) +
           lookup_image(image, i1j1, as_linear) * w.x * w.y;
  }
}

// Apply tone mapping.
inline array2d<vec3f> tonemap_image(const array2d<vec3f>& image, float exposure,
    bool filmic = false, bool srgb = true) {
  return fmap(image, [=](const vec3f& pixel) {
    return tonemap(pixel, exposure, filmic, srgb);
  });
}
inline array2d<vec4f> tonemap_image(const array2d<vec4f>& image, float exposure,
    bool filmic = false, bool srgb = true) {
  return fmap(image, [=](const vec4f& pixel) {
    return tonemap(pixel, exposure, filmic, srgb);
  });
}
inline void tonemap_image(array2d<vec3f>& result, const array2d<vec3f>& image,
    float exposure, bool filmic = false, bool srgb = true) {
  return fmap(result, image, [=](const vec3f& pixel) {
    return tonemap(pixel, exposure, filmic, srgb);
  });
}
inline void tonemap_image(array2d<vec4f>& result, const array2d<vec4f>& image,
    float exposure, bool filmic = false, bool srgb = true) {
  return fmap(result, image, [=](const vec4f& pixel) {
    return tonemap(pixel, exposure, filmic, srgb);
  });
}

// Get/Set region
template <typename T>
inline array2d<T> get_region(
    const array2d<T>& image, const vec2i& offset, const vec2i& extents) {
  auto region = array2d<T>(extents);
  for (auto ij : range(region.extents())) region[ij] = image[ij + offset];
  return region;
}
template <typename T>
inline void get_region(array2d<T>& region, const array2d<T>& image,
    const vec2i& offset, const vec2i& extents) {
  if (region.extents() != extents) region = array2d<T>(extents);
  for (auto ij : range(region.extents())) region[ij] = image[ij + offset];
}
template <typename T>
inline void set_region(
    array2d<T>& image, const array2d<T>& region, const vec2i& offset) {
  for (auto ij : range(region.extents())) image[ij + offset] = region[ij];
}

// Compute the difference between two images.
inline array2d<vec4f> image_difference(
    const array2d<vec4f>& image1, const array2d<vec4f>& image2, bool display) {
  // check sizes
  check_same_size(image1, image2);

  // compute diff
  auto difference = array2d<vec4f>(image1.extents());
  for (auto idx : range(difference.extents())) {
    auto diff       = abs(image1[idx] - image2[idx]);
    difference[idx] = display ? vec4f{max(diff), max(diff), max(diff), 1}
                              : diff;
  }
  return difference;
}

// Composite two images together.
inline array2d<vec4f> composite_image(
    const array2d<vec4f>& image_a, const array2d<vec4f>& image_b) {
  if (image_a.extents() != image_b.extents())
    throw std::invalid_argument{"image should be the same size"};
  auto result = array2d<vec4f>(image_a.extents());
  for (auto ij : range(result.extents())) {
    result[ij] = composite(image_a[ij], image_b[ij]);
  }
  return result;
}

// Composite two images together.
inline void composite_image(array2d<vec4f>& result,
    const array2d<vec4f>& image_a, const array2d<vec4f>& image_b) {
  if (image_a.extents() != image_b.extents())
    throw std::invalid_argument{"image should be the same size"};
  if (image_a.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  for (auto ij : range(result.extents())) {
    result[ij] = composite(image_a[ij], image_b[ij]);
  }
}

// Color grade an hsr or ldr image to an ldr image.
inline array2d<vec4f> colorgrade_image(
    const array2d<vec4f>& image, bool linear, const colorgrade_params& params) {
  return fmap(image,
      [&](const vec4f& pixel) { return colorgrade(pixel, linear, params); });
}

// Color grade an hsr or ldr image to an ldr image.
inline void colorgrade_image(array2d<vec4f>& result,
    const array2d<vec4f>& image, bool linear, const colorgrade_params& params) {
  return fmap(result, image,
      [&](const vec4f& pixel) { return colorgrade(pixel, linear, params); });
}

// determine white balance colors
inline vec3f compute_white_balance(const array2d<vec4f>& image) {
  auto rgb = vec3f{0, 0, 0};
  for (auto ij : range(image.extents())) rgb += xyz(image[ij]);
  if (rgb == vec3f{0, 0, 0}) return vec3f{0, 0, 0};
  rgb /= max(rgb);
  return rgb;
}

// Resize an image.
array2d<vec3f> resize_image(const array2d<vec3f>& image, const vec2i& extents);
array2d<vec4f> resize_image(const array2d<vec4f>& image, const vec2i& extents);
array2d<vec3b> resize_image(const array2d<vec3b>& image, const vec2i& extents);
array2d<vec4b> resize_image(const array2d<vec4b>& image, const vec2i& extents);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image
array2d<vec4f> make_grid(const vec2i& extents = {1024, 1024}, float scale = 1,
    const vec4f& color0 = {0.5, 0.5, 0.5, 1.0},
    const vec4f& color1 = {0.7, 0.7, 0.7, 1.0});

array2d<vec4f> make_checker(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0.5, 0.5, 0.6, 1.0},
    const vec4f& color1 = {0.7, 0.7, 0.7, 1.0});

array2d<vec4f> make_bumps(const vec2i& extents = {1024, 1024}, float scale = 1,
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});

array2d<vec4f> make_ramp(const vec2i& extents = {1024, 1024}, float scale = 1,
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});

array2d<vec4f> make_gammaramp(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});

array2d<vec4f> make_uvramp(
    const vec2i& extents = {1024, 1024}, float scale = 1);

array2d<vec4f> make_orgrid(
    const vec2i& extents = {1024, 1024}, float scale = 1);

array2d<vec4f> make_uvgrid(
    const vec2i& extents = {1024, 1024}, float scale = 1, bool colored = true);

array2d<vec4f> make_colormapramp(
    const vec2i& extents = {1024, 1024}, float scale = 1);

array2d<vec4f> make_gnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
array2d<vec4f> make_vnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
array2d<vec4f> make_fnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
array2d<vec4f> make_tnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
array2d<vec4f> make_rnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});

// Add image border
array2d<vec4f> add_border(const array2d<vec4f>& image, float width,
    const vec4f& color = {0, 0, 0, 1});

// Comvert a bump map to a normal map.
void bump_to_normal(
    array2d<vec4f>& normalmap, const array2d<vec4f>& bumpmap, float scale);
array2d<vec4f> bump_to_normal(const array2d<vec4f>& bumpmap, float scale);

// Implementation of sunsky modified heavily from pbrt
array2d<vec4f> make_sunsky(const vec2i& extents = {2048, 1024},
    float theta_sun = pif / 4, float turbidity = 3, bool has_sun = false,
    float sun_intensity = 1, float sun_radius = 1,
    const vec3f& ground_albedo = vec3f{0.7, 0.7, 0.7});

// Make an image of multiple lights.
array2d<vec4f> make_lights(const vec2i& extents = {2048, 1024},
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
    float lwidth = pif / 16, float lheight = pif / 16);

}  // namespace yocto

#endif
