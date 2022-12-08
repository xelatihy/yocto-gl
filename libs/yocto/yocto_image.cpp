//
// Implementation for Yocto/Image.
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_image.h"

#include <stb_image/stb_image_resize.h>

#include <future>
#include <memory>
#include <stdexcept>

#include "yocto_color.h"
#include "yocto_noise.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PARALLEL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for_batch(T num, T batch, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&func, &next_idx, &has_error, num, batch]() {
          try {
            while (true) {
              auto start = next_idx.fetch_add(batch);
              if (start >= num) break;
              if (has_error) break;
              auto end = std::min(num, start + batch);
              for (auto i = (T)start; i < end; i++) func(i);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for_batch(array<T, 2> num, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<T>    next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, &has_error, num]() {
          try {
            while (true) {
              auto j = next_idx.fetch_add(1);
              if (j >= num[1]) break;
              if (has_error) break;
              for (auto i = (T)0; i < num[0]; i++) func({i, j});
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(array2d<vec4f>& result, const array2d<vec4f>& image,
    float exposure, bool filmic, bool srgb) {
  if (image.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  parallel_for_batch(image.size(), image.extent(0),
      [&result, &image, exposure, filmic, srgb](size_t idx) {
        result[idx] = tonemap(image[idx], exposure, filmic, srgb);
      });
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(array2d<vec4f>& result, const array2d<vec4f>& image,
    bool linear, const colorgrade_params& params) {
  if (image.extents() != result.extents())
    throw std::invalid_argument{"image should be the same size"};
  parallel_for_batch(image.size(), image.extent(0),
      [&result, &image, &params, linear](size_t idx) {
        result[idx] = colorgrade(image[idx], linear, params);
      });
}

// Resize an image.
array2d<vec4f> resize_image(
    const array2d<vec4f>& image, const vec2s& extents_) {
  // determine new size
  auto extents = extents_;
  auto aspect  = (double)image.extent(0) / (double)image.extent(1);
  if (extents == vec2s{0, 0})
    throw std::invalid_argument{"bad image size in resize"};
  if (extents[1] == 0) {
    extents = {extents[0], (size_t)round(extents[0] / aspect)};
  } else if (extents[0] == 0) {
    extents = {(size_t)round(extents[1] * aspect), extents[1]};
  }

  // resize
  auto result = array2d<vec4f>(extents);
  stbir_resize_float_generic((float*)image.data(), (int)image.extent(0),
      (int)image.extent(1), (int)(sizeof(vec4f) * image.extent(0)),
      (float*)result.data(), (int)result.extent(0), (int)result.extent(1),
      (int)(sizeof(vec4f) * result.extent(0)), 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return result;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// image creation
image_data make_image(int width, int height, bool linear) {
  return image_data{
      width, height, linear, vector<vec4f>(width * height, vec4f{0, 0, 0, 0})};
}

// equality
bool operator==(const image_data& a, const image_data& b) {
  return a.width == b.width && a.height == b.height && a.linear == b.linear &&
         a.pixels == b.pixels;
}
bool operator!=(const image_data& a, const image_data& b) {
  return a.width != b.width || a.height != b.height || a.linear != b.linear ||
         a.pixels != b.pixels;
}

// swap
void swap(image_data& a, image_data& b) {
  std::swap(a.width, b.width);
  std::swap(a.height, b.height);
  std::swap(a.linear, b.linear);
  std::swap(a.pixels, b.pixels);
}

// conversions
image_data convert_image(const image_data& image, bool linear) {
  if (image.linear == linear) return image;
  auto result = make_image(image.width, image.height, linear);
  convert_image(result, image);
  return result;
}
void convert_image(image_data& result, const image_data& image) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image have to be the same size"};
  if (image.linear == result.linear) {
    result.pixels = image.pixels;
  } else {
    for (auto idx : range(image.pixels.size())) {
      result.pixels[idx] = image.linear ? rgb_to_srgb(image.pixels[idx])
                                        : srgb_to_rgb(image.pixels[idx]);
    }
  }
}

// Lookup pixel for evaluation
static vec4f lookup_image(
    const image_data& image, const vec2i& ij, bool as_linear) {
  auto [i, j] = ij;
  if (as_linear && !image.linear) {
    return srgb_to_rgb(image.pixels[j * image.width + i]);
  } else {
    return image.pixels[j * image.width + i];
  }
}

// Evaluates an image at a point `uv`.
vec4f eval_image(const image_data& image, const vec2f& uv, bool as_linear,
    bool no_interpolation, bool clamp_to_edge) {
  if (image.pixels.empty()) return vec4f{0, 0, 0, 0};

  // get image width/height
  auto size = vec2i{image.width, image.height};

  // get coordinates normalized for tiling
  auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * size;

  // handle interpolation
  if (no_interpolation) {
    auto ij = clamp((vec2i)st, 0, size - 1);
    return lookup_image(image, ij, as_linear);
  } else {
    auto ij     = clamp((vec2s)st, 0, size - 1);
    auto i1j    = (ij + vec2s{1, 0}) % size;
    auto ij1    = (ij + vec2s{0, 1}) % size;
    auto i1j1   = (ij + vec2s{1, 1}) % size;
    auto [u, v] = st - ij;
    return lookup_image(image, ij, as_linear) * (1 - u) * (1 - v) +
           lookup_image(image, ij1, as_linear) * (1 - u) * v +
           lookup_image(image, i1j, as_linear) * u * (1 - v) +
           lookup_image(image, i1j1, as_linear) * u * v;
  }
}

// Apply tone mapping returning a float or byte image.
image_data tonemap_image(const image_data& image, float exposure, bool filmic) {
  if (!image.linear) return image;
  auto result = make_image(image.width, image.height, false);
  for (auto idx = 0; idx < image.width * image.height; idx++) {
    result.pixels[idx] = tonemap(image.pixels[idx], exposure, filmic, true);
  }
  return result;
}

// Apply tone mapping. If the input image is an ldr, does nothing.
void tonemap_image(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (image.linear) {
    for (auto idx : range(image.pixels.size())) {
      result.pixels[idx] = tonemap(image.pixels[idx], exposure, filmic);
    }
  } else {
    auto scale = vec4f{
        pow(2.0f, exposure), pow(2.0f, exposure), pow(2.0f, exposure), 1};
    for (auto idx : range(image.pixels.size())) {
      result.pixels[idx] = image.pixels[idx] * scale;
    }
  }
}
// Apply tone mapping using multithreading for speed.
void tonemap_image_mt(
    image_data& result, const image_data& image, float exposure, bool filmic) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (result.linear) throw std::invalid_argument{"ldr expected"};
  if (image.linear) {
    parallel_for_batch((size_t)image.width * (size_t)image.height,
        (size_t)image.width, [&result, &image, exposure, filmic](size_t idx) {
          result.pixels[idx] = tonemap(image.pixels[idx], exposure, filmic);
        });
  } else {
    auto scale = vec4f{
        pow(2.0f, exposure), pow(2.0f, exposure), pow(2.0f, exposure), 1};
    parallel_for_batch((size_t)image.width * (size_t)image.height,
        (size_t)image.width, [&result, &image, scale](size_t idx) {
          result.pixels[idx] = image.pixels[idx] * scale;
        });
  }
}

// Resize an image.
image_data resize_image(
    const image_data& image, int res_width, int res_height) {
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
  auto result = make_image(res_width, res_height, image.linear);
  stbir_resize_float_generic((float*)image.pixels.data(), (int)image.width,
      (int)image.height, (int)(sizeof(vec4f) * image.width),
      (float*)result.pixels.data(), (int)result.width, (int)result.height,
      (int)(sizeof(vec4f) * result.width), 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return result;
}

// Compute the difference between two images.
image_data image_difference(
    const image_data& image1, const image_data& image2, bool display) {
  // check sizes
  if (image1.width != image2.width || image1.height != image2.height) {
    throw std::invalid_argument{"image sizes are different"};
  }

  // check types
  if (image1.linear != image2.linear) {
    throw std::invalid_argument{"image types are different"};
  }

  // compute diff
  auto difference = make_image(image1.width, image1.height, image1.linear);
  for (auto idx : range(difference.pixels.size())) {
    auto diff              = abs(image1.pixels[idx] - image2.pixels[idx]);
    difference.pixels[idx] = display ? vec4f{max(diff), max(diff), max(diff), 1}
                                     : diff;
  }
  return difference;
}

void set_region(image_data& image, const image_data& region, int x, int y) {
  for (auto j : range(region.height)) {
    for (auto i : range(region.width)) {
      image.pixels[(j + y) * image.width + (i + x)] =
          region.pixels[j * region.width + i];
    }
  }
}

void get_region(image_data& region, const image_data& image, int x, int y,
    int width, int height) {
  if (region.width != width || region.height != height) {
    region = make_image(width, height, image.linear);
  }
  for (auto j : range(height)) {
    for (auto i : range(width)) {
      region.pixels[j * region.width + i] =
          image.pixels[(j + y) * image.width + (i + x)];
    }
  }
}

// Composite two images together.
image_data composite_image(
    const image_data& image_a, const image_data& image_b) {
  if (image_a.width != image_b.width || image_a.height != image_b.height)
    throw std::invalid_argument{"image should be the same size"};
  if (image_a.linear != image_b.linear)
    throw std::invalid_argument{"image should be of the same type"};
  auto result = make_image(image_a.width, image_a.height, image_a.linear);
  for (auto idx : range(result.pixels.size())) {
    result.pixels[idx] = composite(image_a.pixels[idx], image_b.pixels[idx]);
  }
  return result;
}

// Composite two images together.
void composite_image(
    image_data& result, const image_data& image_a, const image_data& image_b) {
  if (image_a.width != image_b.width || image_a.height != image_b.height)
    throw std::invalid_argument{"image should be the same size"};
  if (image_a.linear != image_b.linear)
    throw std::invalid_argument{"image should be of the same type"};
  if (image_a.width != result.width || image_a.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (image_a.linear != result.linear)
    throw std::invalid_argument{"image should be of the same type"};
  for (auto idx : range(result.pixels.size())) {
    result.pixels[idx] = composite(image_a.pixels[idx], image_b.pixels[idx]);
  }
}

// Color grade an hsr or ldr image to an ldr image.
image_data colorgrade_image(
    const image_data& image, const colorgrade_params& params) {
  auto result = make_image(image.width, image.height, false);
  for (auto idx : range(image.pixels.size())) {
    result.pixels[idx] = colorgrade(image.pixels[idx], image.linear, params);
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
  for (auto idx : range(image.pixels.size())) {
    result.pixels[idx] = colorgrade(image.pixels[idx], image.linear, params);
  }
}

// Color grade an hsr or ldr image to an ldr image.
// Uses multithreading for speed.
void colorgrade_image_mt(image_data& result, const image_data& image,
    const colorgrade_params& params) {
  if (image.width != result.width || image.height != result.height)
    throw std::invalid_argument{"image should be the same size"};
  if (!!result.linear) throw std::invalid_argument{"non linear expected"};
  parallel_for_batch((size_t)image.width * (size_t)image.height,
      (size_t)image.width, [&result, &image, &params](size_t idx) {
        result.pixels[idx] = colorgrade(
            image.pixels[idx], image.linear, params);
      });
}

// determine white balance colors
vec3f compute_white_balance(const image_data& image) {
  auto rgb = vec3f{0, 0, 0};
  for (auto idx = (size_t)0; image.pixels.size(); idx++) {
    rgb += xyz(image.pixels[idx]);
  }
  if (rgb == vec3f{0, 0, 0}) return {0, 0, 0};
  rgb /= max(rgb);
  return rgb;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Lookup an image at coordinates `ij`
static vec4f lookup_image(
    const vector<vec4f>& img, vec2i size, vec2i ij, bool as_linear) {
  auto [width, _] = size;
  auto [i, j]     = ij;
  return img[j * width * i];
}
static vec4f lookup_image(
    const vector<vec4b>& img, vec2i size, vec2i ij, bool as_linear) {
  auto [width, _] = size;
  auto [i, j]     = ij;
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
  auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * size;

  // handle interpolation
  if (no_interpolation) {
    auto ij = clamp((vec2s)st, 0, size - 1);
    return lookup_image(img, size, ij, as_linear);
  } else {
    auto ij     = clamp((vec2s)st, 0, size - 1);
    auto i1j    = (ij + vec2s{1, 0}) % size;
    auto ij1    = (ij + vec2s{0, 1}) % size;
    auto i1j1   = (ij + vec2s{1, 1}) % size;
    auto [u, v] = st - ij;
    return lookup_image(img, size, ij, as_linear) * (1 - u) * (1 - v) +
           lookup_image(img, size, ij1, as_linear) * (1 - u) * v +
           lookup_image(img, size, i1j, as_linear) * u * (1 - v) +
           lookup_image(img, size, i1j1, as_linear) * u * v;
  }
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
  for (auto i : range(fl.size())) fl[i] = byte_to_float(bt[i]);
}
void float_to_byte(vector<vec4b>& bt, const vector<vec4f>& fl) {
  bt.resize(fl.size());
  for (auto i : range(bt.size())) bt[i] = float_to_byte(fl[i]);
}

// Conversion between linear and gamma-encoded images.
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4f>& srgb) {
  rgb.resize(srgb.size());
  for (auto i : range(rgb.size())) rgb[i] = srgb_to_rgb(srgb[i]);
}
void rgb_to_srgb(vector<vec4f>& srgb, const vector<vec4f>& rgb) {
  srgb.resize(rgb.size());
  for (auto i : range(srgb.size())) srgb[i] = rgb_to_srgb(rgb[i]);
}
void srgb_to_rgb(vector<vec4f>& rgb, const vector<vec4b>& srgb) {
  rgb.resize(srgb.size());
  for (auto i : range(rgb.size())) rgb[i] = srgb_to_rgb(byte_to_float(srgb[i]));
}
void rgb_to_srgb(vector<vec4b>& srgb, const vector<vec4f>& rgb) {
  srgb.resize(rgb.size());
  for (auto i : range(srgb.size()))
    srgb[i] = float_to_byte(rgb_to_srgb(rgb[i]));
}

// Apply exposure and filmic tone mapping
void tonemap_image(vector<vec4f>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic, bool srgb) {
  ldr.resize(hdr.size());
  for (auto i : range(hdr.size()))
    ldr[i] = tonemap(hdr[i], exposure, filmic, srgb);
}
void tonemap_image(vector<vec4b>& ldr, const vector<vec4f>& hdr, float exposure,
    bool filmic, bool srgb) {
  ldr.resize(hdr.size());
  for (auto i : range(hdr.size()))
    ldr[i] = float_to_byte(tonemap(hdr[i], exposure, filmic, srgb));
}

void tonemap_image_mt(vector<vec4f>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for_batch(hdr.size(), (size_t)1024,
      [&](size_t i) { ldr[i] = tonemap(hdr[i], exposure, filmic, srgb); });
}
void tonemap_image_mt(vector<vec4b>& ldr, const vector<vec4f>& hdr,
    float exposure, bool filmic, bool srgb) {
  parallel_for_batch(hdr.size(), (size_t)1024, [&](size_t i) {
    ldr[i] = float_to_byte(tonemap(hdr[i], exposure, filmic, srgb));
  });
}

// Apply exposure and filmic tone mapping
void colorgrade_image(vector<vec4f>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  corrected.resize(img.size());
  for (auto i : range(img.size()))
    corrected[i] = colorgrade(img[i], linear, params);
}

// Apply exposure and filmic tone mapping
void colorgrade_image_mt(vector<vec4f>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for_batch(img.size(), (size_t)1024,
      [&](size_t i) { corrected[i] = colorgrade(img[i], linear, params); });
}
void colorgrade_image_mt(vector<vec4b>& corrected, const vector<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  parallel_for_batch(img.size(), (size_t)1024, [&](size_t i) {
    corrected[i] = float_to_byte(colorgrade(img[i], linear, params));
  });
}

// compute white balance
vec3f compute_white_balance(const vector<vec4f>& img) {
  auto rgb = vec3f{0, 0, 0};
  for (auto& p : img) rgb += xyz(p);
  if (rgb == vec3f{0, 0, 0}) return {0, 0, 0};
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
  for (auto i : range(diff.size())) diff[i] = abs(a[i] - b[i]);
  if (display) {
    for (auto i : range(diff.size())) {
      auto d  = max(diff[i]);
      diff[i] = {d, d, d, 1};
    }
  }
}

}  // namespace yocto
