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
image_data make_image(int width, int height, bool linear, bool as_byte) {
  if (!as_byte) {
    return image_data{width, height, linear,
        vector<vec4f>(width * height, vec4f{0, 0, 0, 0}), {}};
  } else {
    return image_data{width, height, linear, {},
        vector<vec4b>(width * height, vec4b{0, 0, 0, 0})};
  }
}
image_data make_image(int width, int height, bool linear, const vec4f* data) {
  return image_data{
      width, height, linear, vector<vec4f>(data, data + width * height), {}};
}
image_data make_image(int width, int height, bool linear, const vec4b* data) {
  return image_data{
      width, height, linear, {}, vector<vec4b>(data, data + width * height)};
}

// equality
bool operator==(const image_data& a, const image_data& b) {
  return a.width == b.width && a.height == b.height && a.linear == b.linear &&
         a.pixelsf == b.pixelsf && a.pixelsb == b.pixelsb;
}
bool operator!=(const image_data& a, const image_data& b) {
  return a.width != b.width || a.height != b.height || a.linear != b.linear ||
         a.pixelsf != b.pixelsf || a.pixelsb != b.pixelsb;
}

// swap
void swap(image_data& a, image_data& b) {
  std::swap(a.width, b.width);
  std::swap(a.height, b.height);
  std::swap(a.linear, b.linear);
  std::swap(a.pixelsf, b.pixelsf);
  std::swap(a.pixelsb, b.pixelsb);
}

// pixel access
vec4f get_pixel(const image_data& image, int i, int j) {
  if (!image.pixelsf.empty()) {
    return image.pixelsf[j * image.width + i];
  } else {
    return byte_to_float(image.pixelsb[j * image.width + i]);
  }
}
void set_pixel(image_data& image, int i, int j, const vec4f& pixel) {
  if (!image.pixelsf.empty()) {
    image.pixelsf[j * image.width + i] = pixel;
  } else {
    image.pixelsb[j * image.width + i] = float_to_byte(pixel);
  }
}

// conversions
image_data convert_image(const image_data& image, bool linear, bool as_byte) {
  if (image.linear == linear && !image.pixelsf.empty() == as_byte) return image;
  auto result = make_image(image.width, image.height, linear, as_byte);
  convert_image(result, image);
  return result;
}
void convert_image(image_data& result, const image_data& image) {
  if (image.linear == result.linear) {
    result.pixelsb = image.pixelsb;
    result.pixelsf = image.pixelsf;
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
vec4f eval_image(const image_data& image, const vec2f& uv, bool as_linear,
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
image_data tonemap_image(
    const image_data& image, float exposure, bool filmic, bool as_byte) {
  if (!image.linear) return image;
  auto result = make_image(image.width, image.height, false, as_byte);
  for (auto idx = 0; idx < image.width * image.height; idx++) {
    result.pixelsb[idx] = tonemapb(image.pixelsf[idx], exposure, filmic, true);
  }
  return result;
}

// Apply tone mapping. If the input image is an ldr, does nothing.
void tonemap_image(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"ldr expected"};
  if (!image.linear) throw std::invalid_argument{"hdr expected"};
  for (auto j = 0; j < image.height; j++) {
    for (auto i = 0; i < image.width; i++) {
      auto hdr = get_pixel(image, i, j);
      auto ldr = tonemap(hdr, exposure, filmic);
      set_pixel(result, i, j, ldr);
    }
  }
}
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"ldr expected"};
  if (!image.linear) throw std::invalid_argument{"hdr expected"};
  parallel_for(image.width, image.height,
      [&result, &image, exposure, filmic](int i, int j) {
        auto hdr = get_pixel(image, i, j);
        auto ldr = tonemap(hdr, exposure, filmic);
        set_pixel(result, i, j, ldr);
      });
}

// Resize an image.
image_data resize_image(const image_data& image, int width, int height) {
  if (width == 0 && height == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (height == 0) {
    height = (int)round(width * (double)height / (double)width);
  } else if (width == 0) {
    width = (int)round(height * (double)width / (double)height);
  }
  if (!image.pixelsf.empty()) {
    auto result = make_image(width, height, image.linear, true);
    stbir_resize_uint8_generic((byte*)image.pixelsb.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4b) * image.width),
        (byte*)result.pixelsb.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4b) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  } else {
    auto result = make_image(width, height, image.linear, false);
    stbir_resize_float_generic((float*)image.pixelsf.data(), (int)image.width,
        (int)image.height, (int)(sizeof(vec4f) * image.width),
        (float*)result.pixelsf.data(), (int)result.width, (int)result.height,
        (int)(sizeof(vec4f) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
        STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
    return result;
  }
}

// Compute the difference between two images.
image_data image_difference(
    const image_data& image1, const image_data& image2, bool display) {
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

void set_region(image_data& image, const image_data& region, int x, int y) {
  for (auto j = 0; j < region.height; j++) {
    for (auto i = 0; i < region.width; i++) {
      set_pixel(image, i + x, j + y, get_pixel(region, i, j));
    }
  }
}

void get_region(image_data& region, const image_data& image, int x, int y,
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COLOR GRADING
// -----------------------------------------------------------------------------
namespace yocto {

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
image_data colorgrade_image(
    const image_data& image, const colorgrade_params& params, bool as_byte) {
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
void colorgrade_image(image_data& result, const image_data& image,
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
void colorgrade_image_mt(image_data& result, const image_data& image,
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
vec4f compute_white_balance(const image_data& image) {
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
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(
    image_data& normalmap, const image_data& bumpmap, float scale) {
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
image_data bump_to_normal(const image_data& bumpmap, float scale) {
  auto normalmap = make_image(
      bumpmap.width, bumpmap.height, bumpmap.linear, !bumpmap.pixelsf.empty());
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

template <typename Shader>
static image_data make_proc_image(
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
image_data make_grid(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
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

image_data make_checker(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

image_data make_bumps(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
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

image_data make_ramp(int width, int height, float scale, const vec4f& color0,
    const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

image_data make_gammaramp(int width, int height, float scale,
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

image_data make_uvramp(int width, int height, float scale) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

image_data make_uvgrid(int width, int height, float scale, bool colored) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
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

image_data make_blackbodyramp(
    int width, int height, float scale, float from, float to) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

image_data make_colormapramp(int width, int height, float scale) {
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

image_data make_noisemap(int width, int height, float scale,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_fbmmap(int width, int height, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_turbulencemap(int width, int height, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image_data make_ridgemap(int width, int height, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_proc_image(width, height, false, false, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
image_data add_border(
    const image_data& image, float width, const vec4f& color) {
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
image_data make_sunsky(int width, int height, float theta_sun, float turbidity,
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
image_data make_lights(int width, int height, const vec3f& le, int nlights,
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// RGB color space definition. Various predefined color spaces are listed below.
struct color_space_params {
  // Curve type
  enum struct curve_t {
    linear,
    gamma,
    linear_gamma,
    aces_cc,
    aces_cct,
    pq,
    hlg
  };
  // primaries
  vec2f red_chromaticity;    // xy chromaticity of the red primary
  vec2f green_chromaticity;  // xy chromaticity of the green primary
  vec2f blue_chromaticity;   // xy chromaticity of the blue primary
  vec2f white_chromaticity;  // xy chromaticity of the white point
  mat3f rgb_to_xyz_mat;      // matrix from rgb to xyz
  mat3f xyz_to_rgb_mat;      // matrix from xyz to rgb
  // tone curve
  curve_t curve_type;
  float   curve_gamma;  // gamma for power curves
  vec4f   curve_abcd;   // tone curve values for linear_gamma curves
};

// Compute the rgb -> xyz matrix from the color space definition
// Input: red, green, blue, white (x,y) chromoticities
// Algorithm from: SMPTE Recommended Practice RP 177-1993
// http://car.france3.mars.free.fr/HD/INA-%2026%20jan%2006/SMPTE%20normes%20et%20confs/rp177.pdf
inline mat3f rgb_to_xyz_mat(
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
inline color_space_params get_color_scape_params(color_space space) {
  static auto make_linear_rgb_space = [](const vec2f& red, const vec2f& green,
                                          const vec2f& blue,
                                          const vec2f& white) {
    return color_space_params{red, green, blue, white,
        rgb_to_xyz_mat(red, green, blue, white),
        inverse(rgb_to_xyz_mat(red, green, blue, white)),
        color_space_params::curve_t::linear};
  };
  static auto make_gamma_rgb_space =
      [](const vec2f& red, const vec2f& green, const vec2f& blue,
          const vec2f& white, float gamma, const vec4f& curve_abcd = zero4f) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)),
            curve_abcd == zero4f ? color_space_params::curve_t::gamma
                                 : color_space_params::curve_t::linear_gamma};
      };
  static auto make_other_rgb_space =
      [](const vec2f& red, const vec2f& green, const vec2f& blue,
          const vec2f& white, color_space_params::curve_t curve_type) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)), curve_type};
      };

  // color space parameters
  // https://en.wikipedia.org/wiki/Rec._709
  static auto rgb_params = make_linear_rgb_space(
      {0.6400, 0.3300}, {0.3000, 0.6000}, {0.1500, 0.0600}, {0.3127, 0.3290});
  // https://en.wikipedia.org/wiki/Rec._709
  static auto srgb_params = make_gamma_rgb_space({0.6400, 0.3300},
      {0.3000, 0.6000}, {0.1500, 0.0600}, {0.3127, 0.3290}, 2.4,
      {1.055, 0.055, 12.92, 0.0031308});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_System
  static auto aces2065_params = make_linear_rgb_space({0.7347, 0.2653},
      {0.0000, 1.0000}, {0.0001, -0.0770}, {0.32168, 0.33767});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  static auto acescg_params = make_linear_rgb_space({0.7130, 0.2930},
      {0.1650, 0.8300}, {0.1280, +0.0440}, {0.32168, 0.33767});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  static auto acescc_params = make_other_rgb_space({0.7130, 0.2930},
      {0.1650, 0.8300}, {0.1280, +0.0440}, {0.32168, 0.33767},
      color_space_params::curve_t::aces_cc);
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  static auto acescct_params = make_other_rgb_space({0.7130, 0.2930},
      {0.1650, 0.8300}, {0.1280, +0.0440}, {0.32168, 0.33767},
      color_space_params::curve_t::aces_cct);
  // https://en.wikipedia.org/wiki/Adobe_RGB_color_space
  static auto adobe_params = make_gamma_rgb_space({0.6400, 0.3300},
      {0.2100, 0.7100}, {0.1500, 0.0600}, {0.3127, 0.3290}, 2.19921875);
  // https://en.wikipedia.org/wiki/Rec._709
  static auto rec709_params = make_gamma_rgb_space({0.6400, 0.3300},
      {0.3000, 0.6000}, {0.1500, 0.0600}, {0.3127, 0.3290}, 1 / 0.45,
      {1.099, 0.099, 4.500, 0.018});
  // https://en.wikipedia.org/wiki/Rec._2020
  static auto rec2020_params = make_gamma_rgb_space({0.7080, 0.2920},
      {0.1700, 0.7970}, {0.1310, 0.0460}, {0.3127, 0.3290}, 1 / 0.45,
      {1.09929682680944, 0.09929682680944, 4.5, 0.018053968510807});
  // https://en.wikipedia.org/wiki/Rec._2020
  static auto rec2100pq_params = make_other_rgb_space({0.7080, 0.2920},
      {0.1700, 0.7970}, {0.1310, 0.0460}, {0.3127, 0.3290},
      color_space_params::curve_t::pq);
  // https://en.wikipedia.org/wiki/Rec._2020
  static auto rec2100hlg_params = make_other_rgb_space({0.7080, 0.2920},
      {0.1700, 0.7970}, {0.1310, 0.0460}, {0.3127, 0.3290},
      color_space_params::curve_t::hlg);
  // https://en.wikipedia.org/wiki/DCI-P3
  static auto p3dci_params = make_gamma_rgb_space({0.6800, 0.3200},
      {0.2650, 0.6900}, {0.1500, 0.0600}, {0.3140, 0.3510}, 1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  static auto p3d60_params = make_gamma_rgb_space({0.6800, 0.3200},
      {0.2650, 0.6900}, {0.1500, 0.0600}, {0.32168, 0.33767}, 1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  static auto p3d65_params = make_gamma_rgb_space({0.6800, 0.3200},
      {0.2650, 0.6900}, {0.1500, 0.0600}, {0.3127, 0.3290}, 1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  static auto p3display_params = make_gamma_rgb_space({0.6800, 0.3200},
      {0.2650, 0.6900}, {0.1500, 0.0600}, {0.3127, 0.3290}, 2.4,
      {1.055, 0.055, 12.92, 0.0031308});
  // https://en.wikipedia.org/wiki/ProPhoto_RGB_color_space
  static auto prophoto_params = make_gamma_rgb_space({0.7347, 0.2653},
      {0.1596, 0.8404}, {0.0366, 0.0001}, {0.3457, 0.3585}, 1.8,
      {1, 0, 16, 0.001953125});

  // return values;
  switch (space) {
    case color_space::rgb: return rgb_params;
    case color_space::srgb: return srgb_params;
    case color_space::adobe: return adobe_params;
    case color_space::prophoto: return prophoto_params;
    case color_space::rec709: return rec709_params;
    case color_space::rec2020: return rec2020_params;
    case color_space::rec2100pq: return rec2100pq_params;
    case color_space::rec2100hlg: return rec2100hlg_params;
    case color_space::aces2065: return aces2065_params;
    case color_space::acescg: return acescg_params;
    case color_space::acescc: return acescc_params;
    case color_space::acescct: return acescct_params;
    case color_space::p3dci: return p3dci_params;
    case color_space::p3d60: return p3d60_params;
    case color_space::p3d65: return p3d65_params;
    case color_space::p3display: return p3display_params;
    default: throw std::runtime_error{"should not have gotten here"};
  }

  // return here to silence warnings
  throw std::runtime_error{"should not have gotten here"};
  return {};
}

// gamma to linear
inline float gamma_display_to_linear(float x, float gamma) {
  return pow(x, gamma);
}
inline float gamma_linear_to_display(float x, float gamma) {
  return pow(x, 1 / gamma);
}

// https://en.wikipedia.org/wiki/Rec._709
inline float gamma_display_to_linear(float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  if (x < 1 / d) {
    return x / c;
  } else {
    return pow((x + b) / a, gamma);
  }
}
inline float gamma_linear_to_display(float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  if (x < d) {
    return x * c;
  } else {
    return a * pow(x, 1 / gamma) - b;
  }
}

// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
inline float acescc_display_to_linear(float x) {
  if (x < -0.3013698630f) {  // (9.72-15)/17.52
    return (exp2(x * 17.52f - 9.72f) - exp2(-16.0f)) * 2;
  } else if (x < (log2(65504.0f) + 9.72f) / 17.52f) {
    return exp2(x * 17.52f - 9.72f);
  } else {  // (in >= (log2(65504)+9.72)/17.52)
    return 65504.0f;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
inline float acescct_display_to_linear(float x) {
  if (x < 0.155251141552511f) {
    return (x - 0.0729055341958355f) / 10.5402377416545f;
  } else {
    return exp2(x * 17.52f - 9.72f);
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
inline float acescc_linear_to_display(float x) {
  if (x <= 0) {
    return -0.3584474886f;  // =(log2( pow(2.,-16.))+9.72)/17.52
  } else if (x < exp2(-15.0f)) {
    return (log2(exp2(-16.0f) + x * 0.5f) + 9.72f) / 17.52f;
  } else {  // (in >= pow(2.,-15))
    return (log2(x) + 9.72f) / 17.52f;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
inline float acescct_linear_to_display(float x) {
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
  return pow((0.8359375f + 18.8515625f * pow(x, 0.1593017578125f)) /
                 (1 + 18.6875f * pow(x, 0.1593017578125f)),
      78.84375f);
}
// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// In HLG, we assume that the linear luminance in [0,1] corresponds to
// [0,1000] cd m^2. Note that the version we report here is scaled in [0,1]
// range for nominal luminance. But HLG was initially defined in the [0,12]
// range where it maps 1 to 0.5 and 12 to 1. For use in HDR tonemapping that is
// likely a better range to use.
inline float hlg_display_to_linear(float x) {
  if (x < 0.5f) {
    return 3 * 3 * x * x;
  } else {
    return (exp((x - 0.55991073f) / 0.17883277f) + 0.28466892f) / 12;
  }
}
inline float hlg_linear_to_display(float x) {
  if (x < 1 / 12.0f) {
    return sqrt(3 * x);
  } else {
    return 0.17883277f * log(12 * x - 0.28466892f) + 0.55991073f;
  }
}

// Conversion to/from xyz
vec3f color_to_xyz(const vec3f& col, color_space from) {
  auto space = get_color_scape_params(from);
  auto rgb   = col;
  if (space.curve_type == color_space_params::curve_t::linear) {
    // do nothing
  } else if (space.curve_type == color_space_params::curve_t::gamma) {
    rgb = {
        gamma_linear_to_display(rgb.x, space.curve_gamma),
        gamma_linear_to_display(rgb.y, space.curve_gamma),
        gamma_linear_to_display(rgb.z, space.curve_gamma),
    };
  } else if (space.curve_type == color_space_params::curve_t::linear_gamma) {
    rgb = {
        gamma_linear_to_display(rgb.x, space.curve_gamma, space.curve_abcd),
        gamma_linear_to_display(rgb.y, space.curve_gamma, space.curve_abcd),
        gamma_linear_to_display(rgb.z, space.curve_gamma, space.curve_abcd),
    };
  } else if (space.curve_type == color_space_params::curve_t::aces_cc) {
    rgb = {
        acescc_linear_to_display(rgb.x),
        acescc_linear_to_display(rgb.y),
        acescc_linear_to_display(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::aces_cct) {
    rgb = {
        acescct_linear_to_display(rgb.x),
        acescct_linear_to_display(rgb.y),
        acescct_linear_to_display(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::pq) {
    rgb = {
        pq_linear_to_display(rgb.x),
        pq_linear_to_display(rgb.y),
        pq_linear_to_display(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::hlg) {
    rgb = {
        hlg_linear_to_display(rgb.x),
        hlg_linear_to_display(rgb.y),
        hlg_linear_to_display(rgb.z),
    };
  } else {
    throw std::runtime_error{"should not have gotten here"};
  }
  return space.rgb_to_xyz_mat * rgb;
}
vec3f xyz_to_color(const vec3f& xyz, color_space to) {
  auto space = get_color_scape_params(to);
  auto rgb   = space.xyz_to_rgb_mat * xyz;
  if (space.curve_type == color_space_params::curve_t::linear) {
    // nothing
  } else if (space.curve_type == color_space_params::curve_t::gamma) {
    rgb = {
        gamma_display_to_linear(rgb.x, space.curve_gamma),
        gamma_display_to_linear(rgb.y, space.curve_gamma),
        gamma_display_to_linear(rgb.z, space.curve_gamma),
    };
  } else if (space.curve_type == color_space_params::curve_t::linear_gamma) {
    rgb = {
        gamma_display_to_linear(rgb.x, space.curve_gamma, space.curve_abcd),
        gamma_display_to_linear(rgb.y, space.curve_gamma, space.curve_abcd),
        gamma_display_to_linear(rgb.z, space.curve_gamma, space.curve_abcd),
    };
  } else if (space.curve_type == color_space_params::curve_t::aces_cc) {
    rgb = {
        acescc_display_to_linear(rgb.x),
        acescc_display_to_linear(rgb.y),
        acescc_display_to_linear(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::aces_cct) {
    rgb = {
        acescct_display_to_linear(rgb.x),
        acescct_display_to_linear(rgb.y),
        acescct_display_to_linear(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::pq) {
    rgb = {
        pq_display_to_linear(rgb.x),
        pq_display_to_linear(rgb.y),
        pq_display_to_linear(rgb.z),
    };
  } else if (space.curve_type == color_space_params::curve_t::hlg) {
    rgb = {
        hlg_display_to_linear(rgb.x),
        hlg_display_to_linear(rgb.y),
        hlg_display_to_linear(rgb.z),
    };
  } else {
    throw std::runtime_error{"should not have gotten here"};
  }
  return rgb;
}

vec3f convert_color(const vec3f& col, color_space from, color_space to) {
  if (from == to) return col;
  return xyz_to_color(color_to_xyz(col, from), to);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup an image at coordinates `ij`
vec4f lookup_image(const image<vec4f>& img, const vec2i& ij, bool as_linear) {
  return img[ij];
}
vec4f lookup_image(const image<vec4b>& img, const vec2i& ij, bool as_linear) {
  if (as_linear) {
    return byte_to_float(img[ij]);
  } else {
    return srgb_to_rgb(byte_to_float(img[ij]));
  }
}

// Lookup an image at coordinates `ij`
vec3f lookup_image(const image<vec3f>& img, const vec2i& ij, bool as_linear) {
  return img[ij];
}
vec3f lookup_image(const image<vec3b>& img, const vec2i& ij, bool as_linear) {
  if (as_linear) {
    return byte_to_float(img[ij]);
  } else {
    return srgb_to_rgb(byte_to_float(img[ij]));
  }
}

// Evaluate a texture
template <typename T, typename R>
static R eval_image_generic(const image<T>& img, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (img.empty()) return R{};

  // get image width/height
  auto size = img.imsize();

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
vec4f eval_image(const image<vec4f>& img, const vec2f& uv,
    bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic<vec4f, vec4f>(
      img, uv, false, no_interpolation, clamp_to_edge);
}
vec4f eval_image(const image<vec4b>& img, const vec2f& uv, bool as_linear,
    bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic<vec4b, vec4f>(
      img, uv, as_linear, no_interpolation, clamp_to_edge);
}

// Evaluates a color image at a point `uv`.
vec3f eval_image(const image<vec3f>& img, const vec2f& uv,
    bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic<vec3f, vec3f>(
      img, uv, false, no_interpolation, clamp_to_edge);
}
vec3f eval_image(const image<vec3b>& img, const vec2f& uv, bool as_linear,
    bool no_interpolation, bool clamp_to_edge) {
  return eval_image_generic<vec3b, vec3f>(
      img, uv, as_linear, no_interpolation, clamp_to_edge);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline void set_region(
    image<T>& img, const image<T>& region, const vec2i& offset) {
  for (auto j = 0; j < region.height(); j++) {
    for (auto i = 0; i < region.width(); i++) {
      if (!img.contains({i, j})) continue;
      img[vec2i{i, j} + offset] = region[{i, j}];
    }
  }
}

template <typename T>
inline void get_region(image<T>& clipped, const image<T>& img,
    const vec2i& offset, const vec2i& size) {
  clipped.resize(size);
  for (auto j = 0; j < size.y; j++) {
    for (auto i = 0; i < size.x; i++) {
      clipped[{i, j}] = img[{i + offset.x, j + offset.y}];
    }
  }
}

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt) {
  auto fl = image<vec4f>{bt.imsize()};
  for (auto i = 0ull; i < fl.count(); i++) fl[i] = byte_to_float(bt[i]);
  return fl;
}
image<vec4b> float_to_byte(const image<vec4f>& fl) {
  auto bt = image<vec4b>{fl.imsize()};
  for (auto i = 0ull; i < bt.count(); i++) bt[i] = float_to_byte(fl[i]);
  return bt;
}

// Conversion from/to floats.
image<vec3f> byte_to_float(const image<vec3b>& bt) {
  auto fl = image<vec3f>{bt.imsize()};
  for (auto i = 0ull; i < fl.count(); i++) fl[i] = byte_to_float(bt[i]);
  return fl;
}
image<vec3b> float_to_byte(const image<vec3f>& fl) {
  auto bt = image<vec3b>{fl.imsize()};
  for (auto i = 0ull; i < bt.count(); i++) bt[i] = float_to_byte(fl[i]);
  return bt;
}

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_rgb(const image<vec4f>& srgb) {
  auto rgb = image<vec4f>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++) rgb[i] = srgb_to_rgb(srgb[i]);
  return rgb;
}
image<vec4f> rgb_to_srgb(const image<vec4f>& rgb) {
  auto srgb = image<vec4f>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++) srgb[i] = rgb_to_srgb(rgb[i]);
  return srgb;
}
image<vec4f> srgb_to_rgb(const image<vec4b>& srgb) {
  auto rgb = image<vec4f>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++)
    rgb[i] = srgb_to_rgb(byte_to_float(srgb[i]));
  return rgb;
}
image<vec4b> rgb_to_srgbb(const image<vec4f>& rgb) {
  auto srgb = image<vec4b>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++)
    srgb[i] = float_to_byte(rgb_to_srgb(rgb[i]));
  return srgb;
}

// Conversion between linear and gamma-encoded images.
image<vec3f> srgb_to_rgb(const image<vec3f>& srgb) {
  auto rgb = image<vec3f>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++) rgb[i] = srgb_to_rgb(srgb[i]);
  return rgb;
}
image<vec3f> rgb_to_srgb(const image<vec3f>& rgb) {
  auto srgb = image<vec3f>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++) srgb[i] = rgb_to_srgb(rgb[i]);
  return srgb;
}
image<vec3f> srgb_to_rgb(const image<vec3b>& srgb) {
  auto rgb = image<vec3f>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++)
    rgb[i] = srgb_to_rgb(byte_to_float(srgb[i]));
  return rgb;
}
image<vec3b> rgb_to_srgbb(const image<vec3f>& rgb) {
  auto srgb = image<vec3b>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++)
    srgb[i] = float_to_byte(rgb_to_srgb(rgb[i]));
  return srgb;
}

// Conversion between number of channels.
image<vec4f> rgb_to_rgba(const image<vec3f>& rgb) {
  auto rgba = image<vec4f>{rgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++) rgba[i] = rgb_to_rgba(rgb[i]);
  return rgba;
}
image<vec3f> rgba_to_rgb(const image<vec4f>& rgba) {
  auto rgb = image<vec3f>{rgba.imsize()};
  for (auto i = 0ull; i < rgba.count(); i++) rgb[i] = rgba_to_rgb(rgba[i]);
  return rgb;
}
image<vec4b> rgb_to_rgba(const image<vec3b>& rgb) {
  auto rgba = image<vec4b>{rgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++) rgba[i] = rgb_to_rgba(rgb[i]);
  return rgba;
}
image<vec3b> rgba_to_rgb(const image<vec4b>& rgba) {
  auto rgb = image<vec3b>{rgba.imsize()};
  for (auto i = 0ull; i < rgba.count(); i++) rgb[i] = rgba_to_rgb(rgba[i]);
  return rgb;
}

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
  auto ldr = image<vec4f>{hdr.imsize()};
  for (auto i = 0ull; i < hdr.count(); i++)
    ldr[i] = tonemap(hdr[i], exposure, filmic, srgb);
  return ldr;
}
image<vec4b> tonemap_imageb(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
  auto ldr = image<vec4b>{hdr.imsize()};
  for (auto i = 0ull; i < hdr.count(); i++)
    ldr[i] = float_to_byte(tonemap(hdr[i], exposure, filmic, srgb));
  return ldr;
}

void tonemap_image_mt(image<vec4f>& ldr, const image<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for(hdr.width(), hdr.height(), [&](int i, int j) {
    ldr[{i, j}] = tonemap(hdr[{i, j}], exposure, filmic, srgb);
  });
}
void tonemap_image_mt(image<vec4b>& ldr, const image<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for(hdr.width(), hdr.height(), [&](int i, int j) {
    ldr[{i, j}] = float_to_byte(tonemap(hdr[{i, j}], exposure, filmic, srgb));
  });
}

vec3f colorgrade(
    const vec3f& rgb_, bool linear, const colorgrade_params& params) {
  auto rgb = rgb_;
  if (params.exposure != 0) rgb *= exp2(params.exposure);
  if (params.tint != vec3f{1, 1, 1}) rgb *= params.tint;
  if (params.lincontrast != 0.5f)
    rgb = lincontrast(rgb, params.lincontrast, linear ? 0.18f : 0.5f);
  if (params.logcontrast != 0.5f)
    rgb = logcontrast(rgb, params.logcontrast, linear ? 0.18f : 0.5f);
  if (params.linsaturation != 0.5f) rgb = saturate(rgb, params.linsaturation);
  if (params.filmic) rgb = tonemap_filmic(rgb);
  if (linear && params.srgb) rgb = rgb_to_srgb(rgb);
  if (params.contrast != 0.5f) rgb = contrast(rgb, params.contrast);
  if (params.saturation != 0.5f) rgb = saturate(rgb, params.saturation);
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
    auto lerp_value = clamp(pow(rgb, 1 / gamma), 0, 1);
    rgb             = gain * lerp_value + lift * (1 - lerp_value);
  }
  return rgb;
}
vec4f colorgrade(
    const vec4f& rgba, bool linear, const colorgrade_params& params) {
  auto graded = colorgrade(xyz(rgba), linear, params);
  return {graded.x, graded.y, graded.z, rgba.w};
}

// Apply exposure and filmic tone mapping
image<vec4f> colorgrade_image(
    const image<vec4f>& img, bool linear, const colorgrade_params& params) {
  auto corrected = image<vec4f>{img.imsize()};
  for (auto i = 0ull; i < img.count(); i++)
    corrected[i] = colorgrade(img[i], linear, params);
  return corrected;
}

// Apply exposure and filmic tone mapping
void colorgrade_image_mt(image<vec4f>& corrected, const image<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for(img.width(), img.height(), [&](int i, int j) {
    corrected[{i, j}] = colorgrade(img[{i, j}], linear, params);
  });
}
void colorgrade_image_mt(image<vec4b>& corrected, const image<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for(img.width(), img.height(), [&](int i, int j) {
    corrected[{i, j}] = float_to_byte(colorgrade(img[{i, j}], linear, params));
  });
}

// compute white balance
vec3f compute_white_balance(const image<vec4f>& img) {
  auto rgb = zero3f;
  for (auto& p : img) rgb += xyz(p);
  if (rgb == zero3f) return zero3f;
  return rgb / max(rgb);
}

static vec2i resize_size(const vec2i& img_size, const vec2i& size_) {
  auto size = size_;
  if (size == zero2i) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (size.y == 0) {
    size.y = (int)round(size.x * (float)img_size.y / (float)img_size.x);
  } else if (size.x == 0) {
    size.x = (int)round(size.y * (float)img_size.x / (float)img_size.y);
  }
  return size;
}

image<vec4f> resize_image(const image<vec4f>& img, int width, int height) {
  return resize_image(img, {width, height});
}
image<vec4b> resize_image(const image<vec4b>& img, int width, int height) {
  return resize_image(img, {width, height});
}
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size_) {
  auto size    = resize_size(img.imsize(), size_);
  auto res_img = image<vec4f>{size};
  stbir_resize_float_generic((float*)img.data(), img.width(), img.height(),
      sizeof(vec4f) * img.width(), (float*)res_img.data(), res_img.width(),
      res_img.height(), sizeof(vec4f) * res_img.width(), 4, 3, 0,
      STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return res_img;
}
image<vec4b> resize_image(const image<vec4b>& img, const vec2i& size_) {
  auto size    = resize_size(img.imsize(), size_);
  auto res_img = image<vec4b>{size};
  stbir_resize_uint8_generic((byte*)img.data(), img.width(), img.height(),
      sizeof(vec4b) * img.width(), (byte*)res_img.data(), res_img.width(),
      res_img.height(), sizeof(vec4b) * res_img.width(), 4, 3, 0,
      STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return res_img;
}

image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool display) {
  if (a.imsize() != b.imsize())
    throw std::invalid_argument{"image haev different sizes"};
  auto diff = image<vec4f>{a.imsize()};
  for (auto i = 0llu; i < diff.count(); i++) diff[i] = abs(a[i] - b[i]);
  if (display) {
    for (auto i = 0llu; i < diff.count(); i++) {
      auto d  = max(diff[i]);
      diff[i] = {d, d, d, 1};
    }
  }
  return diff;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(image<vec4f>& norm, const image<vec4f>& img, float scale) {
  norm.resize(img.imsize());
  auto dx = 1.0f / img.width(), dy = 1.0f / img.height();
  for (int j = 0; j < img.height(); j++) {
    for (int i = 0; i < img.width(); i++) {
      auto i1 = (i + 1) % img.width(), j1 = (j + 1) % img.height();
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
image<vec4f> bump_to_normal(const image<vec4f>& img, float scale) {
  auto norm = image<vec4f>{img.imsize()};
  bump_to_normal(norm, img, scale);
  return norm;
}

template <typename Shader>
image<vec4f> make_image(const vec2i& size, Shader&& shader) {
  auto img   = image<vec4f>{size};
  auto scale = 1.0f / max(size);
  for (auto j = 0; j < img.height(); j++) {
    for (auto i = 0; i < img.width(); i++) {
      auto uv     = vec2f{i * scale, j * scale};
      img[{i, j}] = shader(uv);
    }
  }
  return img;
}

// Make an image
image<vec4f> make_grid(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
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

image<vec4f> make_checker(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

image<vec4f> make_bumps(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
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

image<vec4f> make_ramp(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

image<vec4f> make_gammaramp(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
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

image<vec4f> make_uvramp(const vec2i& size, float scale) {
  return make_image(size, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

image<vec4f> make_uvgrid(const vec2i& size, float scale, bool colored) {
  return make_image(size, [=](vec2f uv) {
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

image<vec4f> make_blackbodyramp(
    const vec2i& size, float scale, float from, float to) {
  return make_image(size, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto rgb = blackbody_to_rgb(lerp(from, to, uv.x));
    return vec4f{rgb.x, rgb.y, rgb.z, 1};
  });
}

image<vec4f> make_colormapramp(const vec2i& size, float scale) {
  return make_image(size, [=](vec2f uv) {
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

image<vec4f> make_noisemap(
    const vec2i& size, float scale, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_noise(vec3f{uv.x, uv.y, 0});
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image<vec4f> make_fbmmap(const vec2i& size, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_fbm({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image<vec4f> make_turbulencemap(const vec2i& size, float scale,
    const vec4f& noise, const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_turbulence({uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z);
    v      = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

image<vec4f> make_ridgemap(const vec2i& size, float scale, const vec4f& noise,
    const vec4f& color0, const vec4f& color1) {
  return make_image(size, [=](vec2f uv) {
    uv *= 8 * scale;
    auto v = perlin_ridge(
        {uv.x, uv.y, 0}, noise.x, noise.y, (int)noise.z, noise.w);
    v = clamp(v, 0.0f, 1.0f);
    return lerp(color0, color1, v);
  });
}

// Add image border
image<vec4f> add_border(
    const image<vec4f>& source, float width, const vec4f& color) {
  auto img   = source;
  auto scale = 1.0f / max(img.imsize());
  for (auto j = 0; j < img.height(); j++) {
    for (auto i = 0; i < img.width(); i++) {
      auto uv = vec2f{i * scale, j * scale};
      if (uv.x < width || uv.y < width || uv.x > img.width() * scale - width ||
          uv.y > img.height() * scale - width) {
        img[{i, j}] = color;
      }
    }
  }
  return img;
}

// Implementation of sunsky modified heavily from pbrt
image<vec4f> make_sunsky(const vec2i& size, float theta_sun, float turbidity,
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
  sun_angular_radius = max(sun_angular_radius, 2 * pif / size.y);

  // sun direction
  auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

  auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000 : zero3f;
  };

  // Make the sun sky image
  auto img          = image<vec4f>{size};
  auto sky_integral = 0.0f, sun_integral = 0.0f;
  for (auto j = 0; j < img.height() / 2; j++) {
    auto theta = pif * ((j + 0.5f) / img.height());
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < img.width(); i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / img.width());
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      sky_integral += mean(sky_col) * sin(theta);
      sun_integral += mean(sun_col) * sin(theta);
      auto col    = sky_col + sun_col;
      img[{i, j}] = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != zero3f) {
    auto ground = zero3f;
    for (auto j = 0; j < img.height() / 2; j++) {
      auto theta = pif * ((j + 0.5f) / img.height());
      for (int i = 0; i < img.width(); i++) {
        auto pxl   = img[{i, j}];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (img.width() * img.height());
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = img.height() / 2; j < img.height(); j++) {
      for (int i = 0; i < img.width(); i++) {
        img[{i, j}] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = img.height() / 2; j < img.height(); j++) {
      for (int i = 0; i < img.width(); i++) {
        img[{i, j}] = {0, 0, 0, 1};
      }
    }
  }

  // done
  return img;
}

// Make an image of multiple lights.
image<vec4f> make_lights(const vec2i& size, const vec3f& le, int nlights,
    float langle, float lwidth, float lheight) {
  auto img = image<vec4f>{size};
  for (auto j = 0; j < img.height() / 2; j++) {
    auto theta = pif * ((j + 0.5f) / img.height());
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < img.width(); i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / img.width());
      auto inlight = false;
      for (auto l = 0; l < nlights; l++) {
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
// VOLUME SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup volume
static float lookup_volume(
    const volume<float>& vol, const vec3i& ijk, bool as_linear) {
  return vol[ijk];
}

// Evaluates a color image at a point `uv`.
float eval_volume(const volume<float>& vol, const vec3f& uvw,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (vol.empty()) return 0;

  // get coordinates normalized for tiling
  auto s = clamp((uvw.x + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.width();
  auto t = clamp((uvw.y + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.height();
  auto r = clamp((uvw.z + 1.0f) * 0.5f, 0.0f, 1.0f) * vol.depth();

  // get image coordinates and residuals
  auto i  = clamp((int)s, 0, vol.width() - 1);
  auto j  = clamp((int)t, 0, vol.height() - 1);
  auto k  = clamp((int)r, 0, vol.depth() - 1);
  auto ii = (i + 1) % vol.width(), jj = (j + 1) % vol.height(),
       kk = (k + 1) % vol.depth();
  auto u = s - i, v = t - j, w = r - k;

  // nearest-neighbor interpolation
  if (no_interpolation) {
    i = u < 0.5 ? i : min(i + 1, vol.width() - 1);
    j = v < 0.5 ? j : min(j + 1, vol.height() - 1);
    k = w < 0.5 ? k : min(k + 1, vol.depth() - 1);
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
// IMPLEMENTATION FOR VOLUME
// -----------------------------------------------------------------------------
namespace yocto {

// make a simple example volume
void make_test(
    volume<float>& vol, const vec3i& size, float scale, float exponent) {
  vol.resize(size);
  for (auto k = 0; k < vol.depth(); k++) {
    for (auto j = 0; j < vol.height(); j++) {
      for (auto i = 0; i < vol.width(); i++) {
        auto p     = vec3f{i / (float)vol.width(), j / (float)vol.height(),
            k / (float)vol.depth()};
        auto value = pow(
            max(max(cos(scale * p.x), cos(scale * p.y)), 0.0f), exponent);
        vol[{i, j, k}] = clamp(value, 0.0f, 1.0f);
      }
    }
  }
}
volume<float> make_test(const vec3i& size, float scale, float exponent) {
  auto vol = volume<float>{};
  make_test(vol, size, scale, exponent);
  return vol;
}

void make_volume_preset(volume<float>& vol, const string& type) {
  auto size = vec3i{256, 256, 256};
  if (type == "test-volume") {
    make_test(vol, size, 6, 10);
  } else {
    throw std::invalid_argument{"unknown volume preset " + type};
  }
}
volume<float> make_volume_preset(const string& type) {
  auto vol = volume<float>{};
  make_volume_preset(vol, type);
  return vol;
}

}  // namespace yocto
