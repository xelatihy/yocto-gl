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
