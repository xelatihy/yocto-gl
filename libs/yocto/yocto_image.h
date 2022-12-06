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
#include "yocto_ndarray.h"
#include "yocto_noise.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
template <size_t N, typename T = float>
constexpr array2d<vec<T, N>> byte_to_float(const array2d<vec<byte, N>>& bt) {
  auto fl = array2d<vec<T, N>>{bt.extents()};
  for (auto idx : range(fl.size())) fl[idx] = byte_to_float(bt[idx]);
  return fl;
}
template <typename T, size_t N>
constexpr array2d<vec<byte, N>> float_to_byte(const array2d<vec<T, N>>& fl) {
  auto bt = array2d<vec<byte, N>>{fl.extents()};
  for (auto idx : range(fl.size())) bt[idx] = float_to_byte(fl[idx]);
  return bt;
}
template <typename T, typename Tb>
constexpr void byte_to_float(array2d<T>& fl, const array2d<Tb>& bt) {
  for (auto idx : range(fl.size())) fl[idx] = byte_to_float(bt[idx]);
}
template <typename T, typename Tb>
constexpr void float_to_byte(array2d<Tb>& bt, const array2d<T>& fl) {
  for (auto idx : range(fl.size())) bt[idx] = float_to_byte(fl[idx]);
}

// Conversion between linear and gamma-encoded images.
template <typename T>
constexpr array2d<T> srgb_to_rgb(const array2d<T>& srgb) {
  auto rgb = array2d<T>{srgb.extents()};
  for (auto idx : range(rgb.size())) rgb[idx] = srgb_to_rgb(srgb[idx]);
  return rgb;
}
template <typename T>
constexpr array2d<T> rgb_to_srgb(const array2d<T>& rgb) {
  auto srgb = array2d<T>{rgb.extents()};
  for (auto idx : range(srgb.size())) srgb[idx] = rgb_to_srgb(rgb[idx]);
  return srgb;
}
template <typename T>
constexpr void srgb_to_rgb(array2d<T>& rgb, const array2d<T>& srgb) {
  for (auto idx : range(rgb.size())) rgb[idx] = srgb_to_rgb(srgb[idx]);
}
template <typename T>
constexpr void rgb_to_srgb(array2d<T>& srgb, const array2d<T>& rgb) {
  for (auto idx : range(srgb.size())) srgb[idx] = rgb_to_srgb(rgb[idx]);
}

// Conversion between linear and gamma-encoded images.
template <size_t N, typename T = float>
constexpr array2d<vec<T, N>> srgbb_to_rgb(const array2d<vec<byte, N>>& srgb) {
  auto rgb = array2d<vec<T, N>>{srgb.extents()};
  for (auto idx : range(rgb.size())) rgb[idx] = srgbb_to_rgb<T>(srgb[idx]);
  return rgb;
}
template <typename T, size_t N>
constexpr array2d<vec<byte, N>> rgb_to_srgbb(const array2d<vec<T, N>>& rgb) {
  auto srgb = array2d<vec<byte, N>>{rgb.extents()};
  for (auto idx : range(rgb.size())) srgb[idx] = rgb_to_srgbb(rgb[idx]);
  return srgb;
}
template <typename T, typename Tb>
constexpr void srgbb_to_rgb(array2d<T>& rgb, const array2d<Tb>& srgb) {
  for (auto idx : range(rgb.size())) rgb[idx] = srgbb_to_rgb<float>(srgb[idx]);
}
template <typename T, typename Tb>
constexpr void rgb_to_srgbb(array2d<Tb>& srgb, const array2d<T>& rgb) {
  for (auto idx : range(rgb.size())) srgb[idx] = rgb_to_srgbb(rgb[idx]);
}

// Lookup pixel for evaluation
template <typename T, typename T1, size_t N>
constexpr vec<T, N> lookup_image(
    const array2d<vec<T1, N>>& image, const vec2s& ij, bool as_linear = false) {
  if constexpr (!std::is_same_v<T1, byte>) {
    return as_linear ? srgb_to_rgb(image[ij]) : image[ij];
  } else {
    return as_linear ? srgbb_to_rgb(image[ij]) : byte_to_float(image[ij]);
  }
}

// Evaluates an image at a point `uv`.
template <typename T, typename T1, size_t N>
constexpr vec<T, N> eval_image(const array2d<vec<T1, N>>& image,
    const vec<T, 2>& uv, bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  if (image.empty()) return vec<T, N>{0};

  // get image width/height
  auto size = image.extents();

  // get coordinates normalized for tiling
  auto s = (T)0, t = (T)0;
  if (clamp_to_edge) {
    s = clamp(uv.x, 0, 1) * size.x;
    t = clamp(uv.y, 0, 1) * size.y;
  } else {
    s = fmod(uv.x, 1) * size.x;
    if (s < 0) s += size.x;
    t = fmod(uv.y, 1) * size.y;
    if (t < 0) t += size.y;
  }

  // get image coordinates and residuals
  auto i  = clamp((size_t)s, (size_t)0, size.x - 1),
       j  = clamp((size_t)t, (size_t)0, size.y - 1);
  auto ii = (i + 1) % size.x, jj = (j + 1) % size.y;
  auto u = s - i, v = t - j;

  // handle interpolation
  if (no_interpolation) {
    return lookup_image(image, {i, j}, as_linear);
  } else {
    return lookup_image(image, {i, j}, as_linear) * (1 - u) * (1 - v) +
           lookup_image(image, {i, jj}, as_linear) * (1 - u) * v +
           lookup_image(image, {ii, j}, as_linear) * u * (1 - v) +
           lookup_image(image, {ii, jj}, as_linear) * u * v;
  }
}

// Apply tone mapping returning a float or byte image.
template <typename T, size_t N>
constexpr array2d<vec<T, N>> tonemap_image(const array2d<vec<T, N>>& image,
    T exposure, bool filmic = false, bool srgb = true) {
  auto result = array2d<vec<T, N>>(image.extents());
  for (auto idx : range(image.size())) {
    result[idx] = tonemap(image[idx], exposure, filmic, srgb);
  }
  return result;
}

// Apply tone mapping. If the input image is an ldr, does nothing.
template <typename T, size_t N>
constexpr void tonemap_image(array2d<vec<T, N>>& result,
    const array2d<vec<T, N>>& image, float exposure, bool filmic = false,
    bool srgb = true) {
  if (image.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  for (auto idx : range(image.size())) {
    result[idx] = tonemap(image[idx], exposure, filmic, srgb);
  }
}

// Apply tone mapping. If the input image is an ldr, does nothing.
// void tonemap_image(array2d<vec4f>& ldr, const array2d<vec4f>& image,
//     float exposure, bool filmic = false, bool srgb= true);
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(array2d<vec4f>& ldr, const array2d<vec4f>& image,
    float exposure, bool filmic = false, bool srgb = true);

// Get/Set region
template <typename T>
constexpr void get_region(array2d<T>& region, const array2d<T>& image,
    const vec2s& offset, const vec2s& extents) {
  if (region.extents() != extents) region = array2d<vec4f>(extents);
  for (auto ij : range(region.extents())) region[ij] = image[ij + offset];
}
template <typename T>
constexpr void set_region(
    array2d<T>& image, const array2d<T>& region, const vec2s& offset) {
  for (auto ij : range(region.extents())) image[ij + offset] = region[ij];
}

// Compute the difference between two images.
inline array2d<vec4f> image_difference(
    const array2d<vec4f>& image1, const array2d<vec4f>& image2, bool display) {
  // check sizes
  if (image1.extents() != image2.extents())
    throw std::invalid_argument{"image sizes are different"};

  // compute diff
  auto difference = array2d<vec4f>(image1.extents());
  for (auto idx : range(difference.size())) {
    auto diff       = abs(image1[idx] - image2[idx]);
    difference[idx] = display ? vec4f{max(diff), max(diff), max(diff), 1}
                              : diff;
  }
  return difference;
}

// Composite two images together.
template <typename T>
inline array2d<vec<T, 4>> composite_image(
    const array2d<vec<T, 4>>& image_a, const array2d<vec<T, 4>>& image_b) {
  if (image_a.extents() != image_b.extents())
    throw std::invalid_argument{"image should be the same size"};
  auto result = array2d<vec<T, 4>>(image_a.extents());
  for (auto idx : range(result.size())) {
    result[idx] = composite(image_a[idx], image_b[idx]);
  }
  return result;
}

// Composite two images together.
template <typename T>
inline void composite_image(array2d<vec<T, 4>>& result,
    const array2d<vec<T, 4>>& image_a, const array2d<vec<T, 4>>& image_b) {
  if (image_a.extents() != image_b.extents())
    throw std::invalid_argument{"image should be the same size"};
  if (image_a.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  for (auto idx : range(result.size())) {
    result[idx] = composite(image_a[idx], image_b[idx]);
  }
}

// Color grade an hsr or ldr image to an ldr image.
template <typename T>
inline array2d<vec<T, 4>> colorgrade_image(const array2d<vec<T, 4>>& image,
    bool linear, const colorgrade_gparams<T>& params) {
  auto result = array2d<vec<T, 4>>(image.extents());
  for (auto idx : range(image.size())) {
    result[idx] = colorgrade(image[idx], linear, params);
  }
  return result;
}

// Color grade an hsr or ldr image to an ldr image.
template <typename T>
inline void colorgrade_image(array2d<vec<T, 4>>& result,
    const array2d<vec<T, 4>>& image, bool linear,
    const colorgrade_gparams<T>& params) {
  if (image.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  for (auto idx : range(image.size())) {
    result[idx] = colorgrade(image[idx], linear, params);
  }
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(array2d<vec4f>& result, const array2d<vec4f>& image,
    bool linear, const colorgrade_params& params);

// determine white balance colors
template <typename T>
inline vec<T, 3> compute_white_balance(const array2d<vec<T, 4>>& image) {
  auto rgb = vec<T, 3>{0, 0, 0};
  for (auto idx : range(image.size())) rgb += xyz(image[idx]);
  if (rgb == vec<T, 3>{0, 0, 0}) return vec<T, 3>{0, 0, 0};
  rgb /= max(rgb);
  return rgb;
}

// Resize an image.
array2d<vec4f> resize_image(const array2d<vec4f>& image, const vec2s& extents);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
inline void bump_to_normal(
    array2d<vec4f>& normalmap, const array2d<vec4f>& bumpmap, float scale) {
  if (normalmap.extents() != bumpmap.extents()) 
    throw std::out_of_range{"different image sizes"};
  auto dxy = 1 / (vec2f)bumpmap.extents();
  for (auto ij : range(bumpmap.extents())) {
    auto i1j = (ij + vec2s{1, 0}) % bumpmap.extents();
    auto ij1 = (ij + vec2s{0, 1}) % bumpmap.extents();
      auto p00    = bumpmap[ij],
           p10    = bumpmap[i1j],
           p01    = bumpmap[ij1];
      auto g00    = mean(p00), g10    = mean(p10), g01    = mean(p01);
      auto normal = vec3f{
          scale * (g00 - g10) / dxy[0], scale * (g00 - g01) / dxy[1], 1.0f};
      normal.y = -normal.y;  // make green pointing up, even if y axis
                             // points down
      normal = normalize(normal) * 0.5f + 0.5f;
      normalmap[ij] = {normal.x, normal.y, normal.z, 1};
    }
}
inline array2d<vec4f> bump_to_normal(const array2d<vec4f>& bumpmap, float scale) {
  auto normalmap = array2d<vec4f>{bumpmap.extents()};
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

template <typename Shader>
inline  array2d<vec4f> make_proc_image(
    const vec2s& extents, Shader&& shader) {
  auto image = array2d<vec4f>(extents);
  for (auto ij : range(extents)) {
    auto uv = (vec2f)ij / max(extents);
    image[ij] = shader(uv);
  }
  return image;
}

// Make an image
inline array2d<vec4f> make_grid(const vec2s& extents, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto thick = 0.01f / 2;
    auto c     = uv.x <= thick || uv.x >= 1 - thick || uv.y <= thick ||
             uv.y >= 1 - thick ||
             (uv.x >= 0.5f - thick && uv.x <= 0.5f + thick) ||
             (uv.y >= 0.5f - thick && uv.y <= 0.5f + thick);
    return c ? color0 : color1;
  });
}

inline array2d<vec4f> make_checker(const vec2s& extents, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

inline array2d<vec4f> make_bumps(const vec2s& extents, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto thick  = 0.125f;
    auto center = vec2f{
        uv.x <= 0.5f ? 0.25f : 0.75f,
        uv.y <= 0.5f ? 0.25f : 0.75f,
    };
    auto dist = clamp(length(uv - center), 0.0f, thick) / thick;
    auto val  = uv.x <= 0.5f != uv.y <= 0.5f ? (1 + sqrt(1 - dist)) / 2
                                             : (dist * dist) / 2;
    return lerp(color0, color1, val);
  });
}

inline array2d<vec4f> make_ramp(const vec2s& extents, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

inline array2d<vec4f> make_gammaramp(const vec2s& extents, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    if (uv.y < 1 / 3.0f) {
      return lerp(color0, color1, pow(uv.x, 2.2f));
    } else if (uv.y < 2 / 3.0f) {
      return lerp(color0, color1, uv.x);
    } else {
      return lerp(color0, color1, pow(uv.x, 1 / 2.2f));
    }
  });
}

inline array2d<vec4f> make_uvramp(const vec2s& extents, float scale) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

inline array2d<vec4f> make_uvgrid(const vec2s& extents, float scale, bool colored) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    uv.y     = 1 - uv.y;
    auto hsv = vec3f{0, 0, 0};
    hsv.x    = (clamp((int)(uv.x * 8), 0, 7) +
                (clamp((int)(uv.y * 8), 0, 7) + 5) % 8 * 8) /
            64.0f;
    auto vuv = uv * 4;
    vuv -= vec2f{(float)(int)vuv.x, (float)(int)vuv.y};
    auto vc  = vuv.x <= 0.5f != vuv.y <= 0.5f;
    hsv.z    = vc ? 0.5f - 0.05f : 0.5f + 0.05f;
    auto suv = uv * 16;
    suv -= vec2f{(float)(int)suv.x, (float)(int)suv.y};
    auto st = 0.01f / 2;
    auto sc = suv.x <= st || suv.x >= 1 - st || suv.y <= st || suv.y >= 1 - st;
    if (sc) {
      hsv.y = 0.2f;
      hsv.z = 0.8f;
    } else {
      hsv.y = 0.8f;
    }
    auto rgb = (colored) ? hsv_to_rgb(hsv) : vec3f{hsv.z, hsv.z, hsv.z};
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

inline array2d<vec4f> make_blackbodyramp(
    const vec2s& extents, float scale, float from, float to) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

inline array2d<vec4f> make_colormapramp(const vec2s& extents, float scale) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = vec3f{0, 0, 0};
    if (uv.y < 0.25) {
      rgb = colormap(uv.x, colormap_type::viridis);
    } else if (uv.y < 0.50) {
      rgb = colormap(uv.x, colormap_type::plasma);
    } else if (uv.y < 0.75) {
      rgb = colormap(uv.x, colormap_type::magma);
    } else {
      rgb = colormap(uv.x, colormap_type::inferno);
    }
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

inline array2d<vec4f> make_noisemap(const vec2s& extents, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

inline array2d<vec4f> make_fbmmap(const vec2s& extents, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

inline array2d<vec4f> make_turbulencemap(const vec2s& extents, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

inline array2d<vec4f> make_ridgemap(const vec2s& extents, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(extents, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
inline array2d<vec4f> add_border(
    const array2d<vec4f>& image, float width, const vec4f& color) {
  auto result = image;
  auto [w, h] = image.extents();
  auto scale = 1.0f / max(image.extents());
  for (auto ij : range(image.extents())) {
      auto [u, v] = (vec2f)ij * scale;
      if (u < width || v < width || v > w * scale - width ||
          v > h * scale - width) {
        result[ij] = color;
      }
  }
  return result;
}

// Implementation of sunsky modified heavily from pbrt
inline array2d<vec4f> make_sunsky(const vec2s& extents, float theta_sun, float turbidity,
    bool has_sun, float sun_intensity, float sun_radius,
    const vec3f& ground_albedo) {
  auto zenith_xyY = vec3f{
      (+0.00165f * pow(theta_sun, 3.f) - 0.00374f * pow(theta_sun, 2.f) +
          0.00208f * theta_sun + 0.00000f) *
              pow(turbidity, 2.f) +
          (-0.02902f * pow(theta_sun, 3.f) + 0.06377f * pow(theta_sun, 2.f) -
              0.03202f * theta_sun + 0.00394f) *
              turbidity +
          (+0.11693f * pow(theta_sun, 3.f) - 0.21196f * pow(theta_sun, 2.f) +
              0.06052f * theta_sun + 0.25885f),
      (+0.00275f * pow(theta_sun, 3.f) - 0.00610f * pow(theta_sun, 2.f) +
          0.00316f * theta_sun + 0.00000f) *
              pow(turbidity, 2.f) +
          (-0.04214f * pow(theta_sun, 3.f) + 0.08970f * pow(theta_sun, 2.f) -
              0.04153f * theta_sun + 0.00515f) *
              turbidity +
          (+0.15346f * pow(theta_sun, 3.f) - 0.26756f * pow(theta_sun, 2.f) +
              0.06669f * theta_sun + 0.26688f),
      1000 * (4.0453f * turbidity - 4.9710f) *
              tan((4.0f / 9.0f - turbidity / 120.0f) * (pif - 2 * theta_sun)) -
          .2155f * turbidity + 2.4192f,
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
  auto tauO = exp(-sun_m * sun_ko * .35f);
  auto tauG = exp(
      -1.41f * sun_kg * sun_m / pow(1 + 118.93f * sun_kg * sun_m, 0.45f));
  auto tauWA  = exp(-0.2385f * sun_kwa * 2.0f * sun_m /
                    pow(1 + 20.07f * sun_kwa * 2.0f * sun_m, 0.45f));
  auto sun_le = sun_sol * tauR * tauA * tauO * tauG * tauWA * 10000;

  // rescale by user
  sun_le *= sun_intensity;

  // size
  auto [width, height] = extents;

  // sun scale from Wikipedia scaled by user quantity and rescaled to at
  // the minimum 5 pixel diameter
  auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
  sun_angular_radius *= sun_radius;
  sun_angular_radius = max(sun_angular_radius, 2 * pif / height);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000
                                                   : vec3f{0, 0, 0};
  };

  // Make the sun sky image
  auto img = array2d<vec4f>(extents);
  for (auto j : range(height / 2)) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (auto i : range(width)) {
      auto phi = 2 * pif * (float(i + 0.5f) / width);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      auto col     = sky_col + sun_col;
      img[{i, j}] = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != vec3f{0, 0, 0}) {
    auto ground = vec3f{0, 0, 0};
    for (auto j : range(height / 2)) {
      auto theta = pif * ((j + 0.5f) / height);
      for (auto i : range(width)) {
        auto pxl   = img[{i, j}];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (width * height);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j : range(height / 2, height)) {
      for (auto i : range(width)) {
        img[{i, j}] = {ground, 1};
      }
    }
  } else {
    for (auto j : range(height / 2, height)) {
      for (auto i : range(width)) {
        img[{i, j}] = {0, 0, 0, 1};
      }
    }
  }

  // done
  return img;
}

// Make an image of multiple lights.
inline array2d<vec4f> make_lights(const vec2s& extents, const vec3f& le, int nlights,
    float langle, float lwidth, float lheight) {
  auto [width, height] = extents;
  auto img = array2d<vec4f>(extents);
  for (auto j : range(height)) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (auto i : range(width)) {
      auto phi     = 2 * pif * (float(i + 0.5f) / width);
      auto inlight = false;
      for (auto l : range(nlights)) {
        auto lphi = 2 * pif * (l + 0.5f) / nlights;
        inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
      }
      img[{i, j}] = {le.x, le.y, le.z, 1};
    }
  }
  return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image data as array of float or byte pixels. Images can be stored in linear
// or non linear color space.
struct image_data {
  // image data
  int           width  = 0;
  int           height = 0;
  bool          linear = false;
  vector<vec4f> pixels = {};

  // pixel access
  vec4f&       operator[](vec2i ij);
  const vec4f& operator[](vec2i ij) const;
};

// image creation
image_data make_image(int width, int height, bool linear);

// equality
bool operator==(const image_data& a, const image_data& b);
bool operator!=(const image_data& a, const image_data& b);

// swap
void swap(image_data& a, image_data& b);

// pixel access
inline vec4f get_pixel(const image_data& image, int i, int j);
inline void  set_pixel(image_data& image, int i, int j, const vec4f& pixel);

// conversions
image_data convert_image(const image_data& image, bool linear);
void       convert_image(image_data& result, const image_data& image);

// Evaluates an image at a point `uv`.
vec4f eval_image(const image_data& image, const vec2f& uv,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// Apply tone mapping returning a float or byte image.
image_data tonemap_image(
    const image_data& image, float exposure, bool filmic = false);

// Apply tone mapping. If the input image is an ldr, does nothing.
void tonemap_image(image_data& ldr, const image_data& image, float exposure,
    bool filmic = false);
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(image_data& ldr, const image_data& image, float exposure,
    bool filmic = false);

// Resize an image.
image_data resize_image(const image_data& image, size_t width, size_t height);

// set/get region
void set_region(image_data& image, const image_data& region, int x, int y);
void get_region(image_data& region, const image_data& image, int x, int y,
    int width, int height);

// Compute the difference between two images.
image_data image_difference(
    const image_data& image_a, const image_data& image_b, bool display_diff);

// Composite two images together.
image_data composite_image(
    const image_data& image_a, const image_data& image_b);

// Composite two images together.
void composite_image(
    image_data& result, const image_data& image_a, const image_data& image_b);

// Composite two images together.
void composite_image(image_data& result, const vector<image_data>& images);

// Color grade an hsr or ldr image to an ldr image.
image_data colorgrade_image(
    const image_data& image, const colorgrade_params& params);

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image(image_data& result, const image_data& image,
    const colorgrade_params& params);

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(image_data& result, const image_data& image,
    const colorgrade_params& params);

// determine white balance colors
vec3f compute_white_balance(const image_data& image);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a grid image.
image_data make_grid(int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a checker image.
image_data make_checker(int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0.2f, 0.2f, 0.2f, 1.0f},
    const vec4f& color1 = vec4f{0.5f, 0.5f, 0.5f, 1.0f});
// Make a bump map.
image_data make_bumps(int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a ramp
image_data make_ramp(int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a gamma ramp.
image_data make_gammaramp(int width, int height, float scale = 1,
    const vec4f& color0 = vec4f{0, 0, 0, 1},
    const vec4f& color1 = vec4f{1, 1, 1, 1});
// Make a uv ramp
image_data make_uvramp(int width, int height, float scale = 1);
// Make a uv grid
image_data make_uvgrid(
    int width, int height, float scale = 1, bool colored = true);
// Make blackbody ramp.
image_data make_blackbodyramp(int width, int height, float scale = 1,
    float from = 1000, float to = 12000);
// Make color map ramp.
image_data make_colormapramp(int width, int height, float scale = 1);
// Make a noise image. Noise parameters: lacunarity, gain, octaves, offset.
image_data make_noisemap(int width, int height, float scale = 1,
    const vec4f& color0 = {0, 0, 0, 1}, const vec4f& color1 = {1, 1, 1, 1});
image_data make_fbmmap(int width, int height, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
image_data make_turbulencemap(int width, int height, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});
image_data make_ridgemap(int width, int height, float scale = 1,
    const vec4f& noise = {2, 0.5, 8, 1}, const vec4f& color0 = {0, 0, 0, 1},
    const vec4f& color1 = {1, 1, 1, 1});

// Make a sunsky HDR model with sun at sun_angle elevation in [0,pif/2],
// turbidity in [1.7,10] with or without sun. The sun can be enabled or
// disabled with has_sun. The sun parameters can be slightly modified by
// changing the sun intensity and temperature. Has a convention, a temperature
// of 0 sets the eath sun defaults (ignoring intensity too).
image_data make_sunsky(int width, int height, float sun_angle,
    float turbidity = 3, bool has_sun = false, float sun_intensity = 1,
    float sun_radius = 1, const vec3f& ground_albedo = {0.2f, 0.2f, 0.2f});
// Make an image of multiple lights.
image_data make_lights(int width, int height, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Comvert a bump map to a normal map. All linear color spaces.
image_data bump_to_normal(const image_data& image, float scale = 1);

// Add a border to an image
image_data add_border(
    const image_data& img, float width, const vec4f& color = {0, 0, 0, 1});

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
// BACKWARDS COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

using color_image = image_data;

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

// pixel access
inline vec4f& image_data::operator[](vec2i ij) {
  return pixels[ij.y * width + ij.x];
}
inline const vec4f& image_data::operator[](vec2i ij) const {
  return pixels[ij.y * width + ij.x];
}

// pixel access
inline vec4f get_pixel(const image_data& image, int i, int j) {
  return image.pixels[j * image.width + i];
}
inline void set_pixel(image_data& image, int i, int j, const vec4f& pixel) {
  image.pixels[j * image.width + i] = pixel;
}

}  // namespace yocto

#endif
