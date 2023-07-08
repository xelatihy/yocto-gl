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
  return extents[0] / extents[1];
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
  if constexpr (!(std::is_same_v<T, byte> || std::is_same_v<T, vec1b> ||
                    std::is_same_v<T, vec3b> || std::is_same_v<T, vec3b> ||
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
  auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * size;

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

// Apply tone mapping returning a float or byte image.
template <typename T, size_t N, typename T1>
inline array2d<vec<T, N>> tonemap_image(const array2d<vec<T, N>>& image,
    T1 exposure, bool filmic = false, bool srgb = true) {
  return fmap(image, [=](const vec<T, N>& pixel) {
    return tonemap(pixel, (T)exposure, filmic, srgb);
  });
}

// Apply tone mapping. If the input image is an ldr, does nothing.
template <typename T, size_t N, typename T1>
inline void tonemap_image(array2d<vec<T, N>>& result,
    const array2d<vec<T, N>>& image, T1 exposure, bool filmic = false,
    bool srgb = true) {
  return fmap(result, image, [=](const vec<T, N>& pixel) {
    return tonemap(pixel, (T)exposure, filmic, srgb);
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

// Convert channels
template <size_t M, typename T, size_t N>
inline array2d<vec<T, M>> convert_channels(const array2d<vec<T, N>>& image) {
  if constexpr (N == M) return image;
  return fmap(image,
      [&](const vec<T, N>& pixel) { return convert_channels<M>(pixel); });
}

// Resize an image.
template <size_t N>
array2d<vec<float, N>> resize_image(
    const array2d<vec<float, N>>& image, const vec2i& extents);
template <size_t N>
array2d<vec<byte, N>> resize_image(
    const array2d<vec<byte, N>>& image, const vec2i& extents);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

template <typename Func, typename T = result_t<Func, vec2i>>
inline array2d<T> _make_proc_image(const vec2i& extents, Func&& func) {
  auto image = array2d<T>(extents);
  for (auto ij : range(extents)) image[ij] = func(ij);
  return image;
}

// Make an image
inline array2d<vec4f> make_grid(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0.5, 0.5, 0.5, 1.0},
    const vec4f& color1 = {0.7, 0.7, 0.7, 1.0}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = fmod((4 * scale * ij) / extents, 1);
    auto thick = 0.01f / 2;
    auto c     = uv.x <= thick || uv.x >= 1 - thick || uv.y <= thick ||
             uv.y >= 1 - thick ||
             (uv.x >= 0.5f - thick && uv.x <= 0.5f + thick) ||
             (uv.y >= 0.5f - thick && uv.y <= 0.5f + thick);
    return c ? color0 : color1;
  });
}

inline array2d<vec4f> make_checker(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0.5, 0.5, 0.6, 1.0},
    const vec4f& color1 = {0.7, 0.7, 0.7, 1.0}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv = fmod((4 * scale * ij) / extents, 1);
    auto c  = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

inline array2d<vec4f> make_bumps(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv     = fmod((4 * scale * ij) / extents, 1);
    auto thick  = 0.125f;
    auto center = vec2f{
        uv.x <= 0.5f ? 0.25f : 0.75f,
        uv.y <= 0.5f ? 0.25f : 0.75f,
    };
    auto dist = clamp(length(uv - center), 0, thick) / thick;
    auto val  = uv.x <= 0.5f != uv.y <= 0.5f ? (1 + sqrt(1 - dist)) / 2
                                             : (dist * dist) / 2;
    return lerp(color0, color1, val);
  });
}

inline array2d<vec4f> make_ramp(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv = fmod((scale * ij) / extents, 1);
    return lerp(color0, color1, uv.x);
  });
}

inline array2d<vec4f> make_gammaramp(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = fmod((scale * ij) / extents, 1);
    auto gamma = 2.2f;
    if (uv.y < 1.0f / 3.0f) {
      return lerp(color0, color1, pow(uv.x, gamma));
    } else if (uv.y < 2.0f / 3.0f) {
      return lerp(color0, color1, uv.x);
    } else {
      return lerp(color0, color1, pow(uv.x, 1 / gamma));
    }
  });
}

inline array2d<vec4f> make_uvramp(
    const vec2i& extents = {1024, 1024}, float scale = 1) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv = fmod((scale * ij) / extents, 1);
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

inline array2d<vec4f> make_orgrid(
    const vec2i& extents = {1024, 1024}, float scale = 1) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv = fmod((scale * ij) / extents, 1);
    return uv.x < 0.5f ? (uv.y < 0.5f ? vec4f{0, 0, 0, 1} : vec4f{0, 1, 0, 1})
                       : (uv.y < 0.5f ? vec4f{1, 0, 0, 1} : vec4f{1, 1, 0, 1});
  });
}

inline array2d<vec4f> make_uvgrid(
    const vec2i& extents = {1024, 1024}, float scale = 1, bool colored = true) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv  = flip_v(fmod((scale * ij) / extents, 1));
    auto hue = (clamp((int)(uv.x * 8), 0, 7) +
                   (clamp((int)(uv.y * 8), 0, 7) + 5) % 8 * 8) /
               64.0f;
    auto v          = fmod(uv * 4, 1);
    auto vc         = v.x <= 0.5f != v.y <= 0.5f;
    auto value      = vc ? 0.5f - 0.05f : 0.5f + 0.05f;
    auto s          = fmod(uv * 16, 1);
    auto st         = 0.01f / 2;
    auto sc         = s.x <= st || s.x >= 1 - st || s.y <= st || s.y >= 1 - st;
    auto saturation = 0.0f;
    if (sc) {
      saturation = 0.2f;
      value      = 0.8f;
    } else {
      saturation = 0.8f;
    }
    auto hsv = vec3f{hue, saturation, value};
    auto rgb = (colored) ? hsv_to_rgb(hsv) : vec3f{value, value, value};
    return vec4f{rgb_to_srgb(rgb), 1};
  });
}

inline array2d<vec4f> make_colormapramp(
    const vec2i& extents = {1024, 1024}, float scale = 1) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv  = fmod((scale * ij) / extents, 1);
    auto rgb = vec3f{0, 0, 0};
    if (uv.y < 0.25f) {
      rgb = colormap(uv.x, colormap_type::viridis);
    } else if (uv.y < 0.50f) {
      rgb = colormap(uv.x, colormap_type::plasma);
    } else if (uv.y < 0.75f) {
      rgb = colormap(uv.x, colormap_type::magma);
    } else {
      rgb = colormap(uv.x, colormap_type::inferno);
    }
    return vec4f{rgb, 1};
  });
}

inline array2d<vec4f> make_gnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = (8 * scale * ij) / extents;
    auto value = gradient_noise(uv);
    return lerp(color0, color1, clamp(value, 0, 1));
  });
}

inline array2d<vec4f> make_vnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = (8 * scale * ij) / extents;
    auto value = value_noise(uv);
    return lerp(color0, color1, clamp(value, 0, 1));
  });
}

inline array2d<vec4f> make_fnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = (8 * scale * ij) / extents;
    auto value = fractal_noise(uv);
    return lerp(color0, color1, clamp(value, 0, 1));
  });
}

inline array2d<vec4f> make_tnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = (8 * scale * ij) / extents;
    auto value = turbulence_noise(uv);
    return lerp(color0, color1, clamp(value, 0, 1));
  });
}

inline array2d<vec4f> make_rnoisemap(const vec2i& extents = {1024, 1024},
    float scale = 1, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1}) {
  return _make_proc_image(extents, [=](vec2i ij) -> vec4f {
    auto uv    = (8 * scale * ij) / extents;
    auto value = ridge_noise(uv);
    return lerp(color0, color1, clamp(value, 0, 1));
  });
}

// Add image border
inline array2d<vec4f> add_border(const array2d<vec4f>& image, float width,
    const vec4f& color = {0, 0, 0, 1}) {
  auto result = image;
  auto size   = image.extents();
  auto scale  = 1.0f / max(image.extents());
  for (auto ij : range(image.extents())) {
    auto uv = (vec2f)ij * scale;
    if (uv.x < width || uv.y < width || uv.y > size.x * scale - width ||
        uv.y > size.y * scale - width) {
      result[ij] = color;
    }
  }
  return result;
}

// Comvert a bump map to a normal map.
inline void bump_to_normal(
    array2d<vec4f>& normalmap, const array2d<vec4f>& bumpmap, float scale) {
  if (normalmap.extents() != bumpmap.extents())
    throw std::out_of_range{"different image sizes"};
  auto dxy = 1 / (vec2f)bumpmap.extents();
  for (auto ij : range(bumpmap.extents())) {
    auto i1j = (ij + vec2i{1, 0}) % bumpmap.extents();
    auto ij1 = (ij + vec2i{0, 1}) % bumpmap.extents();
    auto p00 = bumpmap[ij], p10 = bumpmap[i1j], p01 = bumpmap[ij1];
    auto g00 = mean(p00), g10 = mean(p10), g01 = mean(p01);
    auto normal = vec3f{
        scale * (g00 - g10) / dxy[0], scale * (g00 - g01) / dxy[1], 1};
    normal[1] = -normal[1];  // make green pointing up, even if y axis
                             // points down
    normal        = normalize(normal) * 0.5f + 0.5f;
    normalmap[ij] = {normal, 1};
  }
}
inline array2d<vec4f> bump_to_normal(
    const array2d<vec4f>& bumpmap, float scale) {
  auto normalmap = array2d<vec4f>{bumpmap.extents()};
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

// Implementation of sunsky modified heavily from pbrt
inline array2d<vec4f> make_sunsky(const vec2i& extents = {2048, 1024},
    float theta_sun = pif / 4, float turbidity = 3, bool has_sun = false,
    float sun_intensity = 1, float sun_radius = 1,
    const vec3f& ground_albedo = vec3f{0.7, 0.7, 0.7}) {
  auto zenith_xyY = vec3f{
      (0.00165f * pow(theta_sun, 3) - 0.00374f * pow(theta_sun, 2) +
          0.00208f * theta_sun + 0.00000f) *
              pow(turbidity, 2) +
          (-0.02902f * pow(theta_sun, 3) + 0.06377f * pow(theta_sun, 2) -
              0.03202f * theta_sun + 0.00394f) *
              turbidity +
          (+0.11693f * pow(theta_sun, 3) - 0.21196f * pow(theta_sun, 2) +
              0.06052f * theta_sun + 0.25885f),
      (+0.00275f * pow(theta_sun, 3) - 0.00610f * pow(theta_sun, 2) +
          0.00316f * theta_sun + 0.00000f) *
              pow(turbidity, 2) +
          (-0.04214f * pow(theta_sun, 3) + 0.08970f * pow(theta_sun, 2) -
              0.04153f * theta_sun + 0.00515f) *
              turbidity +
          (+0.15346f * pow(theta_sun, 3) - 0.26756f * pow(theta_sun, 2) +
              0.06669f * theta_sun + 0.26688f),
      1000 * (4.0453f * turbidity - 4.9710f) *
              tan((4.0f / 9.0f - turbidity / 120.0f) * (pif - 2 * theta_sun)) -
          0.2155f * turbidity + 24192,
  };

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
    auto num = ((1 + A * exp(B / cos(theta))) *
                (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
    auto den = ((1 + A * exp(B)) * (1 + C * exp(D * theta_sun) +
                                       E * cos(theta_sun) * cos(theta_sun)));
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
  auto sun_ko     = vec3f{0.48f, 0.75f, 0.14f};
  auto sun_kg     = vec3f{0.1f, 0.0f, 0.0f};
  auto sun_kwa    = vec3f{0.02f, 0.0f, 0.0f};
  auto sun_sol    = vec3f{20000.0f, 27000.0f, 30000.0f};
  auto sun_lambda = vec3f{680, 530, 480};
  auto sun_beta   = 0.04608365822050f * turbidity - 0.04586025928522f;
  auto sun_m      = 1.0f /
               (cos(theta_sun) + 0.000940f * pow(1.6386f - theta_sun, -1.253f));

  auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda / 1000, -4.08f));
  auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda / 1000, -1.3f));
  auto tauO = exp(-sun_m * sun_ko * 0.35f);
  auto tauG = exp(
      -1.41f * sun_kg * sun_m / pow(1 + 118.93f * sun_kg * sun_m, 0.45f));
  auto tauWA  = exp(-0.2385f * sun_kwa * 2 * sun_m /
                    pow(1 + 20.07f * sun_kwa * 2 * sun_m, 0.45f));
  auto sun_le = sun_sol * tauR * tauA * tauO * tauG * tauWA * 10000;

  // rescale by user
  sun_le *= sun_intensity;

  // sun scale from Wikipedia scaled by user quantity and rescaled to at
  // the minimum 5 pixel diameter
  auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
  sun_angular_radius *= sun_radius;
  sun_angular_radius = max(sun_angular_radius, 2 * pif / extents.y);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000
                                                   : vec3f{0, 0, 0};
  };

  // Make the sun sky image
  auto img = array2d<vec4f>(extents);
  for (auto j : range(extents.y / 2)) {
    auto theta = pif * ((j + 0.5f) / extents.y);
    theta      = clamp(theta, 0, pif / 2 - flt_eps);
    for (auto i : range(extents.x)) {
      auto phi = 2 * pif * ((i + 0.5f) / extents.x);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1, 1));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      auto col     = sky_col + sun_col;
      img[{i, j}]  = {col, 1};
    }
  }

  if (ground_albedo != vec3f{0, 0, 0}) {
    auto ground = vec3f{0, 0, 0};
    for (auto j : range(extents.y / 2)) {
      auto theta = pif * ((j + 0.5f) / extents.y);
      for (auto i : range(extents.x)) {
        auto le    = xyz(img[{i, j}]);
        auto angle = sin(theta) * 4 * pif / (extents.x * extents.y);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j : range(extents.y / 2, extents.y)) {
      for (auto i : range(extents.x)) {
        img[{i, j}] = {ground, 1};
      }
    }
  } else {
    for (auto j : range(extents.y / 2, extents.y)) {
      for (auto i : range(extents.x)) {
        img[{i, j}] = {0, 0, 0, 1};
      }
    }
  }

  // Something seems to be wrong in the implementation. We scale down the
  // default by a factor.
  for (auto& pixel : img) pixel = {xyz(pixel) / 4, 1};

  // done
  return img;
}

// Make an image of multiple lights.
inline array2d<vec4f> make_lights(const vec2i& extents = {2048, 1024},
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
    float lwidth = pif / 16, float lheight = pif / 16) {
  auto img = array2d<vec4f>(extents);
  for (auto j : range(extents.y)) {
    auto theta = pif * ((j + 0.5f) / extents.y);
    theta      = clamp(theta, 0, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (auto i : range(extents.x)) {
      auto phi     = 2 * pif * ((i + 0.5f) / extents.x);
      auto inlight = false;
      for (auto l : range(nlights)) {
        auto lphi = 2 * pif * (l + 0.5f) / nlights;
        inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
      }
      img[{i, j}] = {le, 1};
    }
  }
  return img;
}

}  // namespace yocto

#endif
