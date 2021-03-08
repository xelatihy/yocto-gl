//
// Implementation for Yocto/Image.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_image.h"

#include <memory>
#include <stdexcept>

#include "ext/stb_image_resize.h"
#include "yocto_color.h"
#include "yocto_noise.h"
#include "yocto_parallel.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// image creation
color_image make_image(int width, int height, bool linear, bool as_byte) {
  if (!as_byte) {
    return color_image{width, height, linear,
        vector<vec4f>(width * height, vec4f{0, 0, 0, 0}), {}};
  } else {
    return color_image{width, height, linear, {},
        vector<vec4b>(width * height, vec4b{0, 0, 0, 0})};
  }
}
color_image make_image(int width, int height, bool linear, const vec4f* data) {
  return color_image{
      width, height, linear, vector<vec4f>(data, data + width * height), {}};
}
color_image make_image(int width, int height, bool linear, const vec4b* data) {
  return color_image{
      width, height, linear, {}, vector<vec4b>(data, data + width * height)};
}

// equality
bool operator==(const color_image& a, const color_image& b) {
  return a.width == b.width && a.height == b.height && a.linear == b.linear &&
         a.pixelsf == b.pixelsf && a.pixelsb == b.pixelsb;
}
bool operator!=(const color_image& a, const color_image& b) {
  return a.width != b.width || a.height != b.height || a.linear != b.linear ||
         a.pixelsf != b.pixelsf || a.pixelsb != b.pixelsb;
}

// swap
void swap(color_image& a, color_image& b) {
  std::swap(a.width, b.width);
  std::swap(a.height, b.height);
  std::swap(a.linear, b.linear);
  std::swap(a.pixelsf, b.pixelsf);
  std::swap(a.pixelsb, b.pixelsb);
}

// pixel access
vec4f get_pixel(const color_image& image, int i, int j) {
  if (!image.pixelsf.empty()) {
    return image.pixelsf[j * image.width + i];
  } else {
    return byte_to_float(image.pixelsb[j * image.width + i]);
  }
}
void set_pixel(color_image& image, int i, int j, const vec4f& pixel) {
  if (!image.pixelsf.empty()) {
    image.pixelsf[j * image.width + i] = pixel;
  } else {
    image.pixelsb[j * image.width + i] = float_to_byte(pixel);
  }
}

// conversions
color_image convert_image(const color_image& image, bool linear, bool as_byte) {
  if (image.linear == linear && !image.pixelsb.empty() == as_byte) return image;
  auto result = make_image(image.width, image.height, linear, as_byte);
  convert_image(result, image);
  return result;
}
void convert_image(color_image& result, const color_image& image) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image have to be the same size"};
  if (image.linear == result.linear &&
      image.pixelsf.empty() == result.pixelsf.empty()) {
    result.pixelsf = image.pixelsf;
    result.pixelsb = image.pixelsb;
  } else if (image.linear == result.linear) {
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        set_pixel(result, i, j, get_pixel(image, i, j));
      }
    }
  } else {
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        auto color     = get_pixel(image, i, j);
        auto converted = image.linear ? rgb_to_srgb(color) : srgb_to_rgb(color);
        set_pixel(result, i, j, converted);
      }
    }
  }
}

// Evaluates an image at a point `uv`.
vec4f eval_image(const color_image& image, const vec2f& uv, bool as_linear,
    bool no_interpolation, bool clamp_to_edge) {
  if (image.width == 0 || image.height == 0) return {0, 0, 0, 0};

  // get image width/height
  auto size = vec2i{image.width, image.height};

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

  if (no_interpolation) {
    if (as_linear && !image.linear) {
      return srgb_to_rgb(get_pixel(image, i, j));
    } else {
      return get_pixel(image, i, j);
    }
  } else {
    // handle interpolation
    if (as_linear && !image.linear) {
      return srgb_to_rgb(get_pixel(image, i, j)) * (1 - u) * (1 - v) +
             srgb_to_rgb(get_pixel(image, i, jj)) * (1 - u) * v +
             srgb_to_rgb(get_pixel(image, ii, j)) * u * (1 - v) +
             srgb_to_rgb(get_pixel(image, ii, jj)) * u * v;
    } else {
      return get_pixel(image, i, j) * (1 - u) * (1 - v) +
             get_pixel(image, i, jj) * (1 - u) * v +
             get_pixel(image, ii, j) * u * (1 - v) +
             get_pixel(image, ii, jj) * u * v;
    }
  }
}

// Apply tone mapping returning a float or byte image.
color_image tonemap_image(
    const color_image& image, float exposure, bool filmic, bool as_byte) {
  if (!image.linear) return image;
  auto result = make_image(image.width, image.height, false, as_byte);
  for (auto idx = 0; idx < image.width * image.height; idx++) {
    result.pixelsb[idx] = float_to_byte(
        tonemap(image.pixelsf[idx], exposure, filmic, true));
  }
  return result;
}

// Apply tone mapping. If the input image is an ldr, does nothing.
void tonemap_image(color_image& result, const color_image& image,
    float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (image.linear) {
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        auto hdr = get_pixel(image, i, j);
        auto ldr = tonemap(hdr, exposure, filmic);
        set_pixel(result, i, j, ldr);
      }
    }
  } else {
    auto scale = vec4f{pow(2, exposure), pow(2, exposure), pow(2, exposure), 1};
    for (auto j = 0; j < image.height; j++) {
      for (auto i = 0; i < image.width; i++) {
        auto hdr = get_pixel(image, i, j);
        auto ldr = hdr * scale;
        set_pixel(result, i, j, ldr);
      }
    }
  }
}
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(color_image& result, const color_image& image,
    float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (image.linear) {
    parallel_for(image.width, image.height,
        [&result, &image, exposure, filmic](int i, int j) {
          auto hdr = get_pixel(image, i, j);
          auto ldr = tonemap(hdr, exposure, filmic);
          set_pixel(result, i, j, ldr);
        });
  } else {
    auto scale = vec4f{pow(2, exposure), pow(2, exposure), pow(2, exposure), 1};
    parallel_for(
        image.width, image.height, [&result, &image, scale](int i, int j) {
          auto hdr = get_pixel(image, i, j);
          auto ldr = hdr * scale;
          set_pixel(result, i, j, ldr);
        });
  }
}

// Resize an image.
color_image resize_image(
    const color_image& image, int res_width, int res_height) {
  if (res_width == 0 && res_height == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (res_height == 0) {
    res_height = (int)round(
        res_width * (double)image.height / (double)image.width);
  } else if (res_width == 0) {
    res_width = (int)round(
        res_height * (double)image.width / (double)image.height);
  }
  if (!image.pixelsf.empty()) {
    auto result = make_image(res_width, res_height, image.linear, false);
    stbir_resize_float_generic((float*)image.pixelsf.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4f) * image.width),
        (float*)result.pixelsf.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4f) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  } else {
    auto result = make_image(res_width, res_height, image.linear, true);
    stbir_resize_uint8_generic((byte*)image.pixelsb.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4b) * image.width),
        (byte*)result.pixelsb.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4b) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  }
}

// Compute the difference between two images.
color_image image_difference(
    const color_image& image1, const color_image& image2, bool display) {
  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    throw std::invalid_argument{"image sizes are different"};
  }

  // check types
  if (!image1.pixelsf.empty() != !image2.pixelsf.empty() ||
      !image1.pixelsf.empty() != !image2.pixelsf.empty()) {
    throw std::invalid_argument{"image types are different"};
  }

  // check types
  if (image1.linear != image2.linear || !image1.linear != !image2.linear) {
    throw std::invalid_argument{"image types are different"};
  }

  // compute diff
  auto difference = make_image(
      image1.width, image1.height, image1.linear, false);
  for (auto j = 0; j < image1.height; j++) {
    for (auto i = 0; i < image1.width; i++) {
      auto diff = abs(get_pixel(image1, i, j) - get_pixel(image2, i, j));
      if (display) {
        auto d = max(diff);
        set_pixel(difference, i, j, {d, d, d, 1});
      } else {
        set_pixel(difference, i, j, diff);
      }
    }
  }
  return difference;
}

void set_region(color_image& image, const color_image& region, int x, int y) {
  for (auto j = 0; j < region.height; j++) {
    for (auto i = 0; i < region.width; i++) {
      set_pixel(image, i + x, j + y, get_pixel(region, i, j));
    }
  }
}

void get_region(color_image& region, const color_image& image, int x, int y,
    int width, int height) {
  if (region.width != width || region.height != height) {
    region = make_image(width, height, image.linear, !image.pixelsf.empty());
  }
  for (auto j = 0; j < height; j++) {
    for (auto i = 0; i < width; i++) {
      set_pixel(region, i, j, get_pixel(region, i + x, j + y));
    }
  }
}

// Apply color grading from a linear or srgb color to an srgb color.
vec4f colorgradeb(
    const vec4f& color, bool linear, const colorgrade_params& params) {
  auto rgb   = xyz(color);
  auto alpha = color.w;
  if (linear) {
    if (params.exposure != 0) rgb *= exp2(params.exposure);
    if (params.tint != vec3f{1, 1, 1}) rgb *= params.tint;
    if (params.lincontrast != 0.5f)
      rgb = lincontrast(rgb, params.lincontrast, 0.18f);
    if (params.logcontrast != 0.5f)
      rgb = logcontrast(rgb, params.logcontrast, 0.18f);
    if (params.linsaturation != 0.5f) rgb = saturate(rgb, params.linsaturation);
    if (params.filmic) rgb = tonemap_filmic(rgb);
    if (params.srgb) rgb = rgb_to_srgb(rgb);
  }
  if (params.contrast != 0.5f) rgb = contrast(rgb, params.contrast);
  if (params.saturation != 0.5f) rgb = saturate(rgb, params.saturation);
  if (params.shadows != 0.5f || params.midtones != 0.5f ||
      params.highlights != 0.5f || params.shadows_color != vec3f{1, 1, 1} ||
      params.midtones_color != vec3f{1, 1, 1} ||
      params.highlights_color != vec3f{1, 1, 1}) {
    auto lift  = params.shadows_color;
    auto gamma = params.midtones_color;
    auto gain  = params.highlights_color;
    lift       = lift - mean(lift) + params.shadows - (float)0.5;
    gain       = gain - mean(gain) + params.highlights + (float)0.5;
    auto grey  = gamma - mean(gamma) + params.midtones;
    gamma      = log(((float)0.5 - lift) / (gain - lift)) / log(grey);
    // apply_image
    auto lerp_value = clamp(pow(rgb, 1 / gamma), 0, 1);
    rgb             = gain * lerp_value + lift * (1 - lerp_value);
  }
  return vec4f{rgb.x, rgb.y, rgb.z, alpha};
}

// Color grade an hsr or ldr image to an ldr image.
color_image colorgrade_image(
    const color_image& image, const colorgrade_params& params, bool as_byte) {
  auto result = make_image(image.width, image.height, false, as_byte);
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto color  = get_pixel(image, i, j);
      auto graded = colorgrade(color, image.linear, params);
      set_pixel(result, i, j, graded);
    }
  }
  return result;
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image(color_image& result, const color_image& image,
    const colorgrade_params& params) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"non linear expected"};
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto color  = get_pixel(image, i, j);
      auto graded = colorgrade(color, image.linear, params);
      set_pixel(result, i, j, graded);
    }
  }
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(color_image& result, const color_image& image,
    const colorgrade_params& params) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"non linear expected"};
  parallel_for(
      image.width, image.height, [&result, &image, &params](int i, int j) {
        auto color  = get_pixel(image, i, j);
        auto graded = colorgrade(color, image.linear, params);
        set_pixel(result, i, j, graded);
      });
}

// determine white balance colors
vec4f compute_white_balance(const color_image& image) {
  auto rgb = vec3f{0, 0, 0};
  for (auto j = 0; image.height; j++) {
    for (auto i = 0; image.width; i++) {
      rgb += xyz(get_pixel(image, i, j));
    }
  }
  if (rgb == vec3f{0, 0, 0}) return {0, 0, 0, 1};
  rgb /= max(rgb);
  return {rgb.x, rgb.y, rgb.z, 1};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(
    color_image& normalmap, const color_image& bumpmap, float scale) {
  auto width = bumpmap.width, height = bumpmap.height;
  if (normalmap.width != bumpmap.width || normalmap.height != bumpmap.height) {
    normalmap = make_image(
        width, height, bumpmap.linear, !bumpmap.pixelsf.empty());
  }
  auto dx = 1.0f / width, dy = 1.0f / height;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      auto i1 = (i + 1) % width, j1 = (j + 1) % height;
      auto p00 = get_pixel(bumpmap, i, j), p10 = get_pixel(bumpmap, i1, j),
           p01    = get_pixel(bumpmap, i, j1);
      auto g00    = (p00.x + p00.y + p00.z) / 3;
      auto g01    = (p01.x + p01.y + p01.z) / 3;
      auto g10    = (p10.x + p10.y + p10.z) / 3;
      auto normal = vec3f{
          scale * (g00 - g10) / dx, scale * (g00 - g01) / dy, 1.0f};
      normal.y = -normal.y;  // make green pointing up, even if y axis
                             // points down
      normal = normalize(normal) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
      set_pixel(normalmap, i, j, {normal.x, normal.y, normal.z, 1});
    }
  }
}
color_image bump_to_normal(const color_image& bumpmap, float scale) {
  auto normalmap = make_image(
      bumpmap.width, bumpmap.height, bumpmap.linear, !bumpmap.pixelsf.empty());
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

template <typename Shader>
static color_image make_proc_image(
    int width, int height, bool linear, bool as_byte, Shader&& shader) {
  auto image = make_image(width, height, linear, as_byte);
  auto scale = 1.0f / max(width, height);
  if (as_byte) {
    for (auto j = 0; j < height; j++) {
      for (auto i = 0; i < width; i++) {
        auto uv                      = vec2f{i * scale, j * scale};
        image.pixelsb[j * width + i] = float_to_byte(shader(uv));
      }
    }
  } else {
    for (auto j = 0; j < height; j++) {
      for (auto i = 0; i < width; i++) {
        auto uv                      = vec2f{i * scale, j * scale};
        image.pixelsf[j * width + i] = shader(uv);
      }
    }
  }
  return image;
}

// Make an image
color_image make_grid(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
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

color_image make_checker(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

color_image make_bumps(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
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

color_image make_ramp(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

color_image make_gammaramp(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
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

color_image make_uvramp(int width, int height, float scale) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

color_image make_uvgrid(int width, int height, float scale, bool colored) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    uv.y     = 1 - uv.y;
    auto hsv = zero3f;
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

color_image make_blackbodyramp(
    int width, int height, float scale, float from, float to) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

color_image make_colormapramp(int width, int height, float scale) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = zero3f;
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

color_image make_noisemap(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

color_image make_fbmmap(int width, int height, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

color_image make_turbulencemap(int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

color_image make_ridgemap(int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, true, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
color_image add_border(
    const color_image& image, float width, const vec4f& color) {
  auto result = image;
  auto scale  = 1.0f / max(image.width, image.height);
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto uv = vec2f{i * scale, j * scale};
      if (uv.x < width || uv.y < width || uv.x > image.width * scale - width ||
          uv.y > image.height * scale - width) {
        set_pixel(result, i, j, color);
      }
    }
  }
  return result;
}

// Implementation of sunsky modified heavily from pbrt
color_image make_sunsky(int width, int height, float theta_sun, float turbidity,
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

  // sun scale from Wikipedia scaled by user quantity and rescaled to at
  // the minimum 5 pixel diamater
  auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
  sun_angular_radius *= sun_radius;
  sun_angular_radius = max(sun_angular_radius, 2 * pif / height);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000 : zero3f;
  };

  // Make the sun sky image
  auto img          = make_image(width, height, true, false);
  auto sky_integral = 0.0f, sun_integral = 0.0f;
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < width; i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / width);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      sky_integral += mean(sky_col) * sin(theta);
      sun_integral += mean(sun_col) * sin(theta);
      auto col                   = sky_col + sun_col;
      img.pixelsf[j * width + i] = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != zero3f) {
    auto ground = zero3f;
    for (auto j = 0; j < height / 2; j++) {
      auto theta = pif * ((j + 0.5f) / height);
      for (int i = 0; i < width; i++) {
        auto pxl   = img.pixelsf[j * width + i];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (width * height);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        img.pixelsf[j * width + i] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        img.pixelsf[j * width + i] = {0, 0, 0, 1};
      }
    }
  }

  // done
  return img;
}

// Make an image of multiple lights.
color_image make_lights(int width, int height, const vec3f& le, int nlights,
    float langle, float lwidth, float lheight) {
  auto img = make_image(width, height, true, false);
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < width; i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / width);
      auto inlight = false;
      for (auto l = 0; l < nlights; l++) {
        auto lphi = 2 * pif * (l + 0.5f) / nlights;
        inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
      }
      img.pixelsf[j * width + i] = {le.x, le.y, le.z, 1};
    }
  }
  return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup an image at coordinates `ij`
static vec4f lookup_image(const vector<vec4f>& img, int width, int height,
    int i, int j, bool as_linear) {
  return img[j * width * i];
}
static vec4f lookup_image(const vector<vec4b>& img, int width, int height,
    int i, int j, bool as_linear) {
  if (as_linear) {
    return srgb_to_rgb(byte_to_float(img[j * width * i]));
  } else {
    return byte_to_float(img[j * width * i]);
  }
}

// Evaluate a texture
template <typename T>
static vec4f eval_image_generic(const vector<T>& img, int width, int height,
    const vec2f& uv, bool as_linear, bool no_interpolation,
    bool clamp_to_edge) {
  if (img.empty()) return vec4f{0, 0, 0, 0};

  // get image width/height
  auto size = vec2i{width, height};

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

  if (no_interpolation)
    return lookup_image(img, width, height, i, j, as_linear);

  // handle interpolation
  return lookup_image(img, width, height, i, j, as_linear) * (1 - u) * (1 - v) +
         lookup_image(img, width, height, i, jj, as_linear) * (1 - u) * v +
         lookup_image(img, width, height, ii, j, as_linear) * u * (1 - v) +
         lookup_image(img, width, height, ii, jj, as_linear) * u * v;
}

// Evaluates a color image at a point `uv`.
vec4f eval_image(const vector<vec4f>& img, int width, int height,
    const vec2f& uv, bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic(
      img, width, height, uv, false, no_interpolation, clamp_to_edge);
}
vec4f eval_image(const vector<vec4b>& img, int width, int height,
    const vec2f& uv, bool as_linear, bool no_interpolation,
    bool clamp_to_edge) {
  return eval_image_generic(
      img, width, height, uv, as_linear, no_interpolation, clamp_to_edge);
}

// Conversion from/to floats.
void byte_to_float(vector<vec4f>& fl, const vector<vec4b>& bt) {
  fl.resize(bt.size());
  for (auto i = 0ull; i < fl.size(); i++) fl[i] = byte_to_float(bt[i]);
}
void float_to_byte(vector<vec4b>& bt, const vector<vec4f>& fl) {
  bt.resize(fl.size());
  for (auto i = 0ull; i < bt.size(); i++) bt[i] = float_to_byte(fl[i]);
}

// Conversion between linear and gamma-encoded images.
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4f>& srgb) {
  rgb.resize(srgb.size());
  for (auto i = 0ull; i < rgb.size(); i++) rgb[i] = srgb_to_rgb(srgb[i]);
}
void rgb_to_srgb(vector<vec4f>& srgb, const vector<vec4f>& rgb) {
  srgb.resize(rgb.size());
  for (auto i = 0ull; i < srgb.size(); i++) srgb[i] = rgb_to_srgb(rgb[i]);
}
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4b>& srgb) {
  rgb.resize(srgb.size());
  for (auto i = 0ull; i < rgb.size(); i++)
    rgb[i] = srgb_to_rgb(byte_to_float(srgb[i]));
}
void rgb_to_srgbb(vector<vec4b>& srgb, const vector<vec4f>& rgb) {
  srgb.resize(rgb.size());
  for (auto i = 0ull; i < srgb.size(); i++)
    srgb[i] = float_to_byte(rgb_to_srgb(rgb[i]));
}

// Apply exposure and filmic tone mapping
void tonemap_image(vector<vec4f>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic, bool srgb) {
  ldr.resize(hdr.size());
  for (auto i = 0ull; i < hdr.size(); i++)
    ldr[i] = tonemap(hdr[i], exposure, filmic, srgb);
}
void tonemap_image(vector<vec4b>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic, bool srgb) {
  ldr.resize(hdr.size());
  for (auto i = 0ull; i < hdr.size(); i++)
    ldr[i] = float_to_byte(tonemap(hdr[i], exposure, filmic, srgb));
}

void tonemap_image_mt(vector<vec4f>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for(hdr.size(),
      [&](size_t i) { ldr[i] = tonemap(hdr[i], exposure, filmic, srgb); });
}
void tonemap_image_mt(vector<vec4b>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for(hdr.size(), [&](size_t i) {
    ldr[i] = float_to_byte(tonemap(hdr[i], exposure, filmic, srgb));
  });
}

// Apply exposure and filmic tone mapping
void colorgrade_image(vector<vec4f>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  corrected.resize(img.size());
  for (auto i = 0ull; i < img.size(); i++)
    corrected[i] = colorgrade(img[i], linear, params);
}

// Apply exposure and filmic tone mapping
void colorgrade_image_mt(vector<vec4f>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for(img.size(),
      [&](size_t i) { corrected[i] = colorgrade(img[i], linear, params); });
}
void colorgrade_image_mt(vector<vec4b>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for(img.size(), [&](size_t i) {
    corrected[i] = float_to_byte(colorgrade(img[i], linear, params));
  });
}

// compute white balance
vec3f compute_white_balance(const vector<vec4f>& img) {
  auto rgb = zero3f;
  for (auto& p : img) rgb += xyz(p);
  if (rgb == zero3f) return zero3f;
  return rgb / max(rgb);
}

void resize_image(vector<vec4f>& res, const vector<vec4f>& img, int width,
    int height, int res_width, int res_height) {
  if (res_width == 0 && res_height == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (res_height == 0) {
    res_height = (int)round(res_width * (double)height / (double)width);
  } else if (res_width == 0) {
    res_width = (int)round(res_height * (double)width / (double)height);
  }
  res.resize((size_t)res_width * (size_t)res_height);
  stbir_resize_float_generic((float*)img.data(), width, height,
      sizeof(vec4f) * width, (float*)res.data(), res_width, res_height,
      sizeof(vec4f) * res_width, 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
}
void resize_image(vector<vec4b>& res, const vector<vec4b>& img, int width,
    int height, int res_width, int res_height) {
  if (res_width == 0 && res_height == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (res_height == 0) {
    res_height = (int)round(res_width * (double)height / (double)width);
  } else if (res_width == 0) {
    res_width = (int)round(res_height * (double)width / (double)height);
  }
  res.resize((size_t)res_width * (size_t)res_height);
  stbir_resize_uint8_generic((byte*)img.data(), width, height,
      sizeof(vec4b) * width, (byte*)res.data(), res_width, res_height,
      sizeof(vec4b) * res_width, 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
}

void image_difference(vector<vec4f>& diff, const vector<vec4f>& a,
    const vector<vec4f>& b, bool display) {
  if (a.size() != b.size())
    throw std::invalid_argument{"image haev different sizes"};
  diff.resize(a.size());
  for (auto i = 0llu; i < diff.size(); i++) diff[i] = abs(a[i] - b[i]);
  if (display) {
    for (auto i = 0llu; i < diff.size(); i++) {
      auto d  = max(diff[i]);
      diff[i] = {d, d, d, 1};
    }
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(vector<vec4f>& normalmap, const vector<vec4f>& bumpmap,
    int width, int height, float scale) {
  auto dx = 1.0f / width, dy = 1.0f / height;
  for (int j = 0; j < height; j++) {
    for (int i = 0; i < width; i++) {
      auto i1 = (i + 1) % width, j1 = (j + 1) % height;
      auto p00 = bumpmap[j * width + i], p10 = bumpmap[j * width + i1],
           p01    = bumpmap[j1 * width + i];
      auto g00    = (p00.x + p00.y + p00.z) / 3;
      auto g01    = (p01.x + p01.y + p01.z) / 3;
      auto g10    = (p10.x + p10.y + p10.z) / 3;
      auto normal = vec3f{
          scale * (g00 - g10) / dx, scale * (g00 - g01) / dy, 1.0f};
      normal.y = -normal.y;  // make green pointing up, even if y axis
                             // points down
      normal = normalize(normal) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
      normalmap[j * width + i] = {normal.x, normal.y, normal.z, 1};
    }
  }
}

template <typename Shader>
static void make_proc_image(
    vector<vec4f>& pixels, int width, int height, Shader&& shader) {
  pixels.resize((size_t)width * (size_t)height);
  auto scale = 1.0f / max(width, height);
  for (auto j = 0; j < height; j++) {
    for (auto i = 0; i < width; i++) {
      auto uv               = vec2f{i * scale, j * scale};
      pixels[j * width + i] = shader(uv);
    }
  }
}

// Make an image
void make_grid(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
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

void make_checker(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

void make_bumps(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
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

void make_ramp(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

void make_gammaramp(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
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

void make_uvramp(vector<vec4f>& pixels, int width, int height, float scale) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

void make_uvgrid(
    vector<vec4f>& pixels, int width, int height, float scale, bool colored) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    uv.y     = 1 - uv.y;
    auto hsv = zero3f;
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

void make_blackbodyramp(vector<vec4f>& pixels, int width, int height,
    float scale, float from, float to) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

void make_colormapramp(
    vector<vec4f>& pixels, int width, int height, float scale) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = zero3f;
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

void make_noisemap(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

void make_fbmmap(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

void make_turbulencemap(vector<vec4f>& pixels, int width, int height,
    float scale, const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

void make_ridgemap(vector<vec4f>& pixels, int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(pixels, width, height, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
void add_border(vector<vec4f>& pixels, const vector<vec4f>& source, int width,
    int height, float thickness, const vec4f& color) {
  pixels     = source;
  auto scale = 1.0f / max(width, height);
  for (auto j = 0; j < height; j++) {
    for (auto i = 0; i < width; i++) {
      auto uv = vec2f{i * scale, j * scale};
      if (uv.x < thickness || uv.y < thickness ||
          uv.x > width * scale - thickness ||
          uv.y > height * scale - thickness) {
        pixels[j * width + i] = color;
      }
    }
  }
}

// Implementation of sunsky modified heavily from pbrt
void make_sunsky(vector<vec4f>& pixels, int width, int height, float theta_sun,
    float turbidity, bool has_sun, float sun_intensity, float sun_radius,
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

  // sun scale from Wikipedia scaled by user quantity and rescaled to at
  // the minimum 5 pixel diamater
  auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
  sun_angular_radius *= sun_radius;
  sun_angular_radius = max(sun_angular_radius, 2 * pif / height);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000 : zero3f;
  };

  // Make the sun sky image
  pixels.resize(width * height);
  auto sky_integral = 0.0f, sun_integral = 0.0f;
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < width; i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / width);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      sky_integral += mean(sky_col) * sin(theta);
      sun_integral += mean(sun_col) * sin(theta);
      auto col              = sky_col + sun_col;
      pixels[j * width + i] = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != zero3f) {
    auto ground = zero3f;
    for (auto j = 0; j < height / 2; j++) {
      auto theta = pif * ((j + 0.5f) / height);
      for (int i = 0; i < width; i++) {
        auto pxl   = pixels[j * width + i];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (width * height);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        pixels[j * width + i] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = height / 2; j < height; j++) {
      for (int i = 0; i < width; i++) {
        pixels[j * width + i] = {0, 0, 0, 1};
      }
    }
  }
}

// Make an image of multiple lights.
void make_lights(vector<vec4f>& pixels, int width, int height, const vec3f& le,
    int nlights, float langle, float lwidth, float lheight) {
  pixels.resize(width * height);
  for (auto j = 0; j < height / 2; j++) {
    auto theta = pif * ((j + 0.5f) / height);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < width; i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / width);
      auto inlight = false;
      for (auto l = 0; l < nlights; l++) {
        auto lphi = 2 * pif * (l + 0.5f) / nlights;
        inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
      }
      pixels[j * width + i] = {le.x, le.y, le.z, 1};
    }
  }
}

}  // namespace yocto
