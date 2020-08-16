//
// Implementation for Yocto/Image.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <atomic>
#include <future>
#include <memory>
#include <thread>

#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#include "yocto_color.h"
#include "yocto_noise.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::future;

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
    default: std::runtime_error("should not have gotten here");
  }

  // return here to silence warnings
  std::runtime_error("should not have gotten here");
  return {};
}

// gamma to linear
inline float gamma_display_to_linear(float x, float gamma) {
  return pow(x, gamma);
};
inline float gamma_linear_to_display(float x, float gamma) {
  return pow(x, 1 / gamma);
};

// https://en.wikipedia.org/wiki/Rec._709
inline float gamma_display_to_linear(float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  if (x < 1 / d) {
    return x / c;
  } else {
    return pow((x + b) / a, gamma);
  }
};
inline float gamma_linear_to_display(float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  if (x < d) {
    return x * c;
  } else {
    return a * pow(x, 1 / gamma) - b;
  }
};

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
    throw std::runtime_error("should not have gotten here");
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
    throw std::runtime_error("should not have gotten here");
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
inline R eval_image_generic(const image<T>& img, const vec2f& uv,
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
  for (auto j = 0; j < region.imsize().y; j++) {
    for (auto i = 0; i < region.imsize().x; i++) {
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

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto        futures  = vector<future<void>>{};
  auto        nthreads = std::thread::hardware_concurrency();
  atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, size]() {
          while (true) {
            auto j = next_idx.fetch_add(1);
            if (j >= size.y) break;
            for (auto i = 0; i < size.x; i++) func({i, j});
          }
        }));
  }
  for (auto& f : futures) f.get();
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

// Conversion from/to floats.
image<float> byte_to_float(const image<byte>& bt) {
  auto fl = image<float>{bt.imsize()};
  for (auto i = 0ull; i < fl.count(); i++) fl[i] = yocto::byte_to_float(bt[i]);
  return fl;
}
image<byte> float_to_byte(const image<float>& fl) {
  auto bt = image<byte>{fl.imsize()};
  for (auto i = 0ull; i < bt.count(); i++) bt[i] = yocto::float_to_byte(fl[i]);
  return bt;
}

// Conversion from/to floats.
image<float> ushort_to_float(const image<ushort>& bt) {
  auto fl = image<float>{bt.imsize()};
  for (auto i = 0ull; i < fl.count(); i++)
    fl[i] = yocto::ushort_to_float(bt[i]);
  return fl;
}
image<ushort> float_to_ushort(const image<float>& fl) {
  auto bt = image<ushort>{fl.imsize()};
  for (auto i = 0ull; i < bt.count(); i++)
    bt[i] = yocto::float_to_ushort(fl[i]);
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

// Conversion between linear and gamma-encoded images.
image<float> srgb_to_rgb(const image<float>& srgb) {
  auto rgb = image<float>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++)
    rgb[i] = yocto::srgb_to_rgb(srgb[i]);
  return rgb;
}
image<float> rgb_to_srgb(const image<float>& rgb) {
  auto srgb = image<float>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++)
    srgb[i] = yocto::rgb_to_srgb(rgb[i]);
  return srgb;
}
image<float> srgb_to_rgb(const image<byte>& srgb) {
  auto rgb = image<float>{srgb.imsize()};
  for (auto i = 0ull; i < rgb.count(); i++)
    rgb[i] = yocto::srgb_to_rgb(yocto::byte_to_float(srgb[i]));
  return rgb;
}
image<byte> rgb_to_srgbb(const image<float>& rgb) {
  auto srgb = image<byte>{rgb.imsize()};
  for (auto i = 0ull; i < srgb.count(); i++)
    srgb[i] = yocto::float_to_byte(yocto::rgb_to_srgb(rgb[i]));
  return srgb;
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
  parallel_for(hdr.imsize(), [&](const vec2i& ij) {
    ldr[ij] = tonemap(hdr[ij], exposure, filmic, srgb);
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
  parallel_for(img.imsize(), [&](const vec2i& ij) {
    corrected[ij] = colorgrade(img[ij], linear, params);
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
    throw std::invalid_argument("bad image size in resize");
  }
  if (size.y == 0) {
    size.y = (int)round(size.x * (float)img_size.y / (float)img_size.x);
  } else if (size.x == 0) {
    size.x = (int)round(size.y * (float)img_size.x / (float)img_size.y);
  }
  return size;
}

image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size_) {
  auto size    = resize_size(img.imsize(), size_);
  auto res_img = image<vec4f>{size};
  stbir_resize_float_generic((float*)img.data(), img.imsize().x, img.imsize().y,
      sizeof(vec4f) * img.imsize().x, (float*)res_img.data(),
      res_img.imsize().x, res_img.imsize().y,
      sizeof(vec4f) * res_img.imsize().x, 4, 3, 0, STBIR_EDGE_CLAMP,
      STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return res_img;
}
image<vec4b> resize_image(const image<vec4b>& img, const vec2i& size_) {
  auto size    = resize_size(img.imsize(), size_);
  auto res_img = image<vec4b>{size};
  stbir_resize_uint8_generic((byte*)img.data(), img.imsize().x, img.imsize().y,
      sizeof(vec4b) * img.imsize().x, (byte*)res_img.data(), res_img.imsize().x,
      res_img.imsize().y, sizeof(vec4b) * res_img.imsize().x, 4, 3, 0,
      STBIR_EDGE_CLAMP, STBIR_FILTER_DEFAULT, STBIR_COLORSPACE_LINEAR, nullptr);
  return res_img;
}

image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool display) {
  if (a.imsize() != b.imsize())
    throw std::invalid_argument("image haev different sizes");
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
  auto dx = 1.0f / img.imsize().x, dy = 1.0f / img.imsize().y;
  for (int j = 0; j < img.imsize().y; j++) {
    for (int i = 0; i < img.imsize().x; i++) {
      auto i1 = (i + 1) % img.imsize().x, j1 = (j + 1) % img.imsize().y;
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
  for (auto j = 0; j < img.imsize().y; j++) {
    for (auto i = 0; i < img.imsize().x; i++) {
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
  for (auto j = 0; j < img.imsize().y; j++) {
    for (auto i = 0; i < img.imsize().x; i++) {
      auto uv = vec2f{i * scale, j * scale};
      if (uv.x < width || uv.y < width ||
          uv.x > img.imsize().x * scale - width ||
          uv.y > img.imsize().y * scale - width) {
        img[{i, j}] = color;
      }
    }
  }
  return img;
};

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
  for (auto j = 0; j < img.imsize().y / 2; j++) {
    auto theta = pif * ((j + 0.5f) / img.imsize().y);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < img.imsize().x; i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / img.imsize().x);
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
    for (auto j = 0; j < img.imsize().y / 2; j++) {
      auto theta = pif * ((j + 0.5f) / img.imsize().y);
      for (int i = 0; i < img.imsize().x; i++) {
        auto pxl   = img[{i, j}];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (img.imsize().x * img.imsize().y);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = img.imsize().y / 2; j < img.imsize().y; j++) {
      for (int i = 0; i < img.imsize().x; i++) {
        img[{i, j}] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = img.imsize().y / 2; j < img.imsize().y; j++) {
      for (int i = 0; i < img.imsize().x; i++) {
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
  for (auto j = 0; j < img.imsize().y / 2; j++) {
    auto theta = pif * ((j + 0.5f) / img.imsize().y);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < img.imsize().x; i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / img.imsize().x);
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

image<vec4b> make_logo(const string& type) {
  static const auto logo_medium_size = vec2i{102, 36};
  static const auto logo_small_size  = vec2i{72, 28};
  // clang-format off
  static const auto logo_medium = vector<byte>{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 139, 193, 55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 209, 255, 43, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 14, 253, 239, 1, 0, 0, 0, 0, 0, 0, 0, 0, 30, 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 167, 234, 190, 0, 0, 0, 0, 0, 0, 59, 234, 233, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 69, 255, 183, 0, 0, 0, 0, 0, 0, 31, 169, 230, 255, 255, 224, 152, 15, 0, 0, 0, 0, 203, 234, 116, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 73, 255, 255, 52, 0, 0, 0, 0, 0, 164, 255, 189, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 54, 228, 253, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 127, 255, 125, 0, 0, 0, 0, 0, 66, 235, 255, 252, 206, 217, 254, 255, 225, 45, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 1, 213, 255, 157, 0, 0, 0, 0, 20, 248, 255, 73, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 185, 255, 66, 0, 0, 0, 0, 43, 251, 255, 182, 15, 0, 0, 22, 166, 188, 8, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 99, 255, 245, 16, 0, 0, 0, 117, 255, 213, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 241, 252, 12, 0, 0, 0, 0, 174, 255, 198, 4, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 7, 232, 255, 111, 0, 0, 1, 219, 255, 99, 0, 0, 0, 29, 152, 206, 238, 193, 128, 6, 0, 0, 0, 0, 0, 0, 25, 149, 207, 239, 199, 118, 4, 43, 179, 199, 255, 255, 181, 179, 136, 0, 0, 0, 0, 61, 166, 220, 231, 180, 95, 0, 0, 0, 0, 0, 0, 0, 0, 45, 255, 206, 0, 0, 0, 0, 52, 255, 255, 82, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 125, 255, 215, 0, 0, 69, 255, 232, 7, 0, 0, 56, 232, 255, 254, 235, 255, 255, 198, 21, 0, 0, 0, 0, 46, 226, 255, 254, 234, 255, 255, 170, 61, 255, 255, 255, 255, 255, 255, 166, 0, 0, 0, 101, 250, 255, 248, 241, 255, 255, 153, 3, 0, 0, 0, 0, 0, 0, 104, 255, 148, 0, 0, 0, 0, 137, 255, 232, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 20, 246, 255, 65, 0, 173, 255, 124, 0, 0, 3, 221, 255, 183, 17, 0, 49, 232, 255, 151, 0, 0, 0, 1, 213, 255, 198, 22, 0, 26, 153, 46, 5, 23, 84, 255, 255, 29, 23, 13, 0, 0, 36, 253, 255, 129, 7, 2, 90, 251, 255, 87, 0, 0, 0, 0, 0, 0, 162, 255, 90, 0, 0, 0, 0, 173, 255, 194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 151, 255, 170, 26, 251, 245, 19, 0, 0, 89, 255, 248, 15, 0, 0, 0, 91, 255, 247, 23, 0, 0, 77, 255, 253, 27, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 0, 153, 255, 199, 0, 0, 0, 0, 155, 255, 206, 0, 0, 0, 0, 0, 0, 220, 255, 32, 0, 0, 0, 0, 208, 255, 167, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 38, 254, 250, 150, 255, 149, 0, 0, 0, 174, 255, 177, 0, 0, 0, 0, 13, 250, 255, 95, 0, 0, 168, 255, 191, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 0, 238, 255, 113, 0, 0, 0, 0, 71, 255, 255, 31, 0, 0, 0, 0, 22, 255, 230, 0, 0, 0, 0, 0, 243, 255, 140, 0, 0, 0, 108, 231, 231, 231, 231, 149, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 177, 255, 255, 253, 36, 0, 0, 0, 205, 255, 140, 0, 0, 0, 0, 0, 227, 255, 125, 0, 0, 201, 255, 149, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 14, 255, 255, 76, 0, 0, 0, 0, 35, 255, 255, 61, 0, 0, 0, 0, 80, 255, 172, 0, 0, 0, 0, 1, 248, 255, 135, 0, 0, 0, 86, 255, 255, 255, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 62, 255, 255, 175, 0, 0, 0, 0, 236, 255, 119, 0, 0, 0, 0, 0, 207, 255, 156, 0, 0, 232, 255, 128, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 44, 255, 255, 55, 0, 0, 0, 0, 15, 255, 255, 92, 0, 0, 0, 0, 138, 255, 113, 0, 0, 0, 0, 0, 222, 255, 155, 0, 0, 0, 1, 4, 4, 166, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 245, 255, 105, 0, 0, 0, 0, 245, 255, 113, 0, 0, 0, 0, 0, 202, 255, 163, 0, 0, 246, 255, 120, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 53, 255, 255, 49, 0, 0, 0, 0, 10, 255, 255, 99, 0, 0, 0, 0, 196, 255, 55, 0, 0, 0, 0, 0, 194, 255, 175, 0, 0, 0, 0, 0, 0, 164, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 215, 255, 133, 0, 0, 0, 0, 0, 221, 255, 132, 0, 0, 218, 255, 141, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 23, 255, 255, 69, 0, 0, 0, 0, 29, 255, 255, 68, 0, 0, 0, 6, 248, 247, 6, 0, 0, 0, 0, 0, 166, 255, 206, 0, 0, 0, 0, 0, 0, 164, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 185, 255, 159, 0, 0, 0, 0, 3, 243, 255, 100, 0, 0, 187, 255, 167, 0, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 6, 0, 0, 0, 1, 247, 255, 95, 0, 0, 0, 0, 54, 255, 255, 36, 0, 0, 0, 57, 255, 195, 0, 0, 0, 0, 0, 0, 103, 255, 255, 32, 0, 0, 0, 0, 0, 164, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 123, 255, 232, 2, 0, 0, 0, 65, 255, 253, 36, 0, 0, 133, 255, 239, 7, 0, 0, 0, 0, 0, 0, 0, 67, 255, 255, 7, 0, 0, 0, 0, 187, 255, 170, 0, 0, 0, 0, 129, 255, 222, 3, 0, 0, 0, 115, 255, 137, 0, 0, 0, 0, 0, 0, 9, 236, 255, 127, 0, 0, 0, 0, 0, 164, 255, 165, 0, 0, 0, 222, 255, 127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 17, 243, 255, 113, 0, 0, 8, 191, 255, 170, 0, 0, 0, 23, 247, 255, 135, 0, 0, 1, 91, 18, 0, 0, 42, 255, 255, 51, 0, 1, 0, 0, 69, 255, 247, 58, 0, 0, 32, 231, 255, 106, 0, 0, 0, 0, 173, 255, 79, 0, 0, 0, 0, 0, 0, 0, 132, 255, 252, 91, 0, 0, 0, 21, 197, 255, 165, 0, 0, 0, 222, 255, 137, 20, 20, 20, 20, 20, 15, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 0, 111, 254, 255, 187, 153, 216, 255, 230, 39, 0, 0, 0, 0, 124, 255, 255, 208, 170, 219, 255, 170, 0, 0, 5, 226, 255, 233, 172, 234, 59, 0, 0, 174, 255, 244, 171, 162, 234, 255, 197, 9, 0, 0, 0, 0, 231, 255, 21, 0, 0, 0, 0, 0, 0, 0, 9, 169, 255, 255, 224, 176, 195, 246, 255, 255, 145, 0, 0, 0, 222, 255, 255, 255, 255, 255, 255, 255, 176, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 244, 255, 104, 0, 0, 0, 0, 0, 0, 87, 233, 255, 255, 253, 196, 29, 0, 0, 0, 0, 0, 0, 94, 237, 255, 255, 255, 188, 34, 0, 0, 0, 65, 233, 255, 255, 246, 103, 0, 0, 2, 138, 244, 255, 255, 247, 160, 7, 0, 0, 0, 0, 33, 255, 219, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 126, 227, 255, 255, 255, 250, 194, 73, 0, 0, 0, 0, 222, 255, 255, 255, 255, 255, 255, 255, 144, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 33, 63, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 36, 64, 18, 0, 0, 0, 0, 0, 0, 11, 56, 60, 15, 0, 0, 0, 0, 0, 3, 47, 55, 6, 0, 0, 0, 0, 0, 0, 91, 255, 160, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 67, 40, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 149, 255, 102, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 69, 140, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  };
static const auto logo_small = vector<byte> {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 165, 221, 125, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 223, 255, 139, 0, 0, 0, 0, 7, 55, 39, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 8, 236, 255, 241, 6, 0, 0, 132, 255, 255, 80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 52, 84, 88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 25, 255, 255, 81, 0, 0, 10, 160, 249, 255, 255, 242, 120, 3, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 138, 255, 255, 72, 0, 0, 214, 255, 226, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 83, 255, 255, 24, 0, 9, 207, 255, 255, 247, 249, 255, 253, 57, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 34, 254, 255, 153, 0, 40, 255, 255, 124, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 140, 255, 222, 0, 0, 115, 255, 255, 160, 9, 9, 129, 129, 0, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 182, 255, 232, 2, 122, 255, 250, 23, 0, 90, 199, 242, 214, 135, 6, 0, 0, 0, 36, 175, 231, 229, 164, 24, 173, 249, 255, 239, 176, 112, 0, 4, 130, 212, 243, 202, 95, 0, 0, 0, 0, 0, 198, 255, 165, 0, 7, 234, 255, 243, 9, 0, 0, 0, 0, 0, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 76, 255, 255, 60, 204, 255, 168, 0, 111, 255, 255, 255, 255, 255, 182, 0, 0, 41, 239, 255, 255, 255, 255, 135, 251, 255, 255, 255, 255, 130, 0, 175, 255, 255, 255, 255, 255, 118, 0, 0, 0, 7, 248, 255, 107, 0, 42, 255, 255, 172, 0, 0, 0, 0, 0, 0, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 2, 224, 255, 172, 255, 255, 61, 14, 237, 255, 208, 31, 152, 255, 255, 65, 0, 170, 255, 254, 90, 63, 143, 7, 58, 240, 255, 216, 59, 24, 60, 255, 255, 157, 31, 204, 255, 240, 16, 0, 0, 57, 255, 255, 50, 0, 77, 255, 255, 141, 0, 43, 61, 61, 61, 30, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 121, 255, 255, 255, 211, 0, 72, 255, 255, 116, 0, 47, 255, 255, 137, 5, 252, 255, 200, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 136, 255, 255, 52, 0, 111, 255, 255, 73, 0, 0, 115, 255, 244, 3, 0, 108, 255, 255, 122, 0, 160, 255, 255, 255, 129, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 22, 249, 255, 255, 105, 0, 105, 255, 255, 89, 0, 20, 255, 255, 169, 37, 255, 255, 166, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 169, 255, 255, 25, 0, 84, 255, 255, 105, 0, 0, 172, 255, 191, 0, 0, 94, 255, 255, 129, 0, 103, 211, 255, 255, 129, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 193, 255, 255, 29, 0, 115, 255, 255, 84, 0, 15, 255, 255, 179, 52, 255, 255, 159, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 179, 255, 255, 20, 0, 79, 255, 255, 115, 0, 0, 230, 255, 133, 0, 0, 67, 255, 255, 145, 0, 0, 40, 255, 255, 129, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 189, 255, 255, 26, 0, 83, 255, 255, 101, 0, 33, 255, 255, 147, 19, 255, 255, 182, 0, 0, 0, 0, 0, 235, 255, 204, 0, 0, 147, 255, 255, 37, 0, 97, 255, 255, 83, 0, 32, 255, 255, 76, 0, 0, 32, 255, 255, 195, 0, 0, 40, 255, 255, 129, 0, 213, 255, 254, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 189, 255, 255, 26, 0, 36, 254, 255, 164, 0, 98, 255, 255, 94, 0, 221, 255, 245, 28, 8, 61, 0, 0, 214, 255, 224, 2, 6, 98, 255, 255, 100, 0, 162, 255, 253, 33, 0, 89, 255, 255, 19, 0, 0, 0, 185, 255, 251, 56, 0, 58, 255, 255, 129, 0, 213, 255, 254, 63, 63, 63, 63, 7, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 189, 255, 255, 26, 0, 0, 168, 255, 255, 197, 244, 255, 222, 5, 0, 93, 255, 255, 244, 239, 255, 77, 0, 136, 255, 255, 242, 178, 6, 225, 255, 244, 197, 255, 255, 163, 0, 0, 147, 255, 217, 0, 0, 0, 0, 61, 248, 255, 251, 216, 254, 255, 255, 129, 0, 213, 255, 255, 255, 255, 255, 254, 9, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 189, 255, 255, 26, 0, 0, 17, 174, 255, 255, 255, 208, 36, 0, 0, 0, 110, 248, 255, 255, 239, 90, 0, 41, 231, 255, 255, 230, 21, 43, 212, 255, 255, 255, 168, 12, 0, 0, 204, 255, 159, 0, 0, 0, 0, 0, 55, 219, 255, 255, 255, 241, 131, 17, 0, 213, 255, 255, 255, 255, 255, 229, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 65, 34, 0, 0, 0, 0, 0, 0, 8, 54, 52, 8, 0, 0, 0, 11, 59, 56, 8, 0, 0, 1, 36, 65, 21, 0, 0, 0, 11, 251, 255, 102, 0, 0, 0, 0, 0, 0, 0, 21, 61, 35, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 161, 221, 44, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
  };
  // clang-format on
  if (type == "logo-medium") {
    auto img = image<vec4b>{logo_medium_size};
    for (auto i = 0; i < img.count(); i++)
      img[i] = vec4b{logo_medium[i], logo_medium[i], logo_medium[i], (byte)255};
    return img;
  } else if (type == "logo-small") {
    auto img = image<vec4b>{logo_small_size};
    for (auto i = 0; i < img.count(); i++)
      img[i] = vec4b{logo_small[i], logo_small[i], logo_small[i], (byte)255};
    return img;
  } else {
    throw std::runtime_error("unknown builtin image " + type);
    return {};
  }
}

image<vec4f> add_logo(const image<vec4f>& img, const string& type) {
  auto logo   = srgb_to_rgb(make_logo(type));
  auto offset = img.imsize() - logo.imsize() - 8;
  auto wlogo  = img;
  set_region(wlogo, logo, offset);
  return wlogo;
}

image<vec4b> add_logo(const image<vec4b>& img, const string& type) {
  auto logo   = make_logo(type);
  auto offset = img.imsize() - logo.imsize() - 8;
  auto wlogo  = img;
  set_region(wlogo, logo, offset);
  return wlogo;
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
// IMPLEMENTATION FOR VOLUME
// -----------------------------------------------------------------------------
namespace yocto {

// make a simple example volume
void make_test(
    volume<float>& vol, const vec3i& size, float scale, float exponent) {
  vol.resize(size);
  for (auto k = 0; k < vol.size().z; k++) {
    for (auto j = 0; j < vol.size().y; j++) {
      for (auto i = 0; i < vol.size().x; i++) {
        auto p     = vec3f{i / (float)vol.size().x, j / (float)vol.size().y,
            k / (float)vol.size().z};
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
    throw std::runtime_error("unknown volume preset " + type);
  }
}
volume<float> make_volume_preset(const string& type) {
  auto vol = volume<float>{};
  make_volume_preset(vol, type);
  return vol;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
static vector<string> split_string(const string& str) {
  auto ret = vector<string>();
  if (str.empty()) return ret;
  auto lpos = (size_t)0;
  while (lpos != str.npos) {
    auto pos = str.find_first_of(" \t\n\r", lpos);
    if (pos != str.npos) {
      if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
      lpos = pos + 1;
    } else {
      if (lpos < str.size()) ret.push_back(str.substr(lpos));
      lpos = pos;
    }
  }
  return ret;
}

// Opens a file with a utf8 file name
static FILE* fopen_utf8(const char* filename, const char* mode) {
#ifdef _Win32
  auto path8 = std::filesystem::u8path(filename);
  auto wmode = std::wstring(string{mode}.begin(), string{mode}.end());
  return _wfopen(path.c_str(), wmode.c_str());
#else
  return fopen(filename, mode);
#endif
}

// Pfm load
static float* load_pfm(const char* filename, int* w, int* h, int* nc, int req) {
  auto fs = fopen_utf8(filename, "rb");
  if (!fs) return nullptr;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  // buffer
  char buffer[4096];
  auto toks = vector<string>();

  // read magic
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  if (toks[0] == "Pf")
    *nc = 1;
  else if (toks[0] == "PF")
    *nc = 3;
  else
    return nullptr;

  // read w, h
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  *w   = atoi(toks[0].c_str());
  *h   = atoi(toks[1].c_str());

  // read scale
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks   = split_string(buffer);
  auto s = atof(toks[0].c_str());

  // read the data (flip y)
  auto npixels = (size_t)(*w) * (size_t)(*h);
  auto nvalues = npixels * (size_t)(*nc);
  auto nrow    = (size_t)(*w) * (size_t)(*nc);
  auto pixels  = std::unique_ptr<float[]>(new float[nvalues]);
  for (auto j = *h - 1; j >= 0; j--) {
    if (fread(pixels.get() + j * nrow, sizeof(float), nrow, fs) != nrow)
      return nullptr;
  }

  // endian conversion
  if (s > 0) {
    for (auto i = 0; i < nvalues; ++i) {
      auto dta = (uint8_t*)(pixels.get() + i);
      std::swap(dta[0], dta[3]);
      std::swap(dta[1], dta[2]);
    }
  }

  // scale
  auto scl = (s > 0) ? s : -s;
  if (scl != 1) {
    for (auto i = 0; i < nvalues; i++) pixels[i] *= scl;
  }

  // proper number of channels
  if (!req || *nc == req) return pixels.release();

  // pack into channels
  if (req < 0 || req > 4) {
    return nullptr;
  }
  auto cpixels = std::unique_ptr<float[]>(new float[req * npixels]);
  for (auto i = 0ull; i < npixels; i++) {
    auto vp = pixels.get() + i * (*nc);
    auto cp = cpixels.get() + i * req;
    if (*nc == 1) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[0];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[0];
          cp[2] = vp[0];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[0];
          cp[2] = vp[0];
          cp[3] = 1;
          break;
      }
    } else {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          cp[3] = 1;
          break;
      }
    }
  }
  return cpixels.release();
}

// save pfm
static bool save_pfm(
    const char* filename, int w, int h, int nc, const float* pixels) {
  auto fs = fopen_utf8(filename, "wb");
  if (!fs) return false;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  if (fprintf(fs, "%s\n", (nc == 1) ? "Pf" : "PF") < 0) return false;
  if (fprintf(fs, "%d %d\n", w, h) < 0) return false;
  if (fprintf(fs, "-1\n") < 0) return false;
  if (nc == 1 || nc == 3) {
    if (fwrite(pixels, sizeof(float), w * h * nc, fs) != w * h * nc)
      return false;
  } else {
    for (auto i = 0; i < w * h; i++) {
      auto vz = 0.0f;
      auto v  = pixels + i * nc;
      if (fwrite(v + 0, sizeof(float), 1, fs) != 1) return false;
      if (fwrite(v + 1, sizeof(float), 1, fs) != 1) return false;
      if (nc == 2) {
        if (fwrite(&vz, sizeof(float), 1, fs) != 1) return false;
      } else {
        if (fwrite(v + 2, sizeof(float), 1, fs) != 1) return false;
      }
    }
  }

  return true;
}

// Get extension (not including '.').
static string path_extension(const string& filename) {
  auto pos = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename) {
  auto ext = path_extension(filename);
  return ext == ".hdr" || ext == ".exr" || ext == ".pfm";
}

// Loads an hdr image.
[[nodiscard]] bool load_image(
    const string& filename, image<vec4f>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0)
      return read_error();
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec4f*)pixels};
    delete[] pixels;
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec4f*)pixels};
    free(pixels);
    return true;
  } else if (!is_hdr_filename(filename)) {
    auto imgb = image<vec4b>{};
    if (!load_image(filename, imgb, error)) return false;
    img = srgb_to_rgb(imgb);
    return true;
  } else {
    return format_error();
  }
}

// Saves an hdr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<vec4f>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    if (!stbi_write_hdr(filename.c_str(), img.imsize().x, img.imsize().y, 4,
            (float*)img.data()))
      return write_error();
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), img.imsize().x, img.imsize().y, 4,
            (float*)img.data()))
      return write_error();
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((float*)img.data(), img.imsize().x, img.imsize().y, 4,
            filename.c_str()) < 0)
      return write_error();
    return true;
  } else if (!is_hdr_filename(filename)) {
    return save_image(filename, rgb_to_srgbb(img), error);
  } else {
    return format_error();
  }
}

// Loads an ldr image.
[[nodiscard]] bool load_image(
    const string& filename, image<vec4b>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
      ext == ".tga" || ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 4);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec4b*)pixels};
    free(pixels);
    return true;
  } else if (is_hdr_filename(filename)) {
    auto imgf = image<vec4f>{};
    if (!load_image(filename, imgf, error)) return false;
    img = rgb_to_srgbb(imgf);
    return true;
  } else {
    return format_error();
  }
}

// Saves an ldr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<vec4b>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), img.imsize().x, img.imsize().y, 4,
            img.data(), img.imsize().x * 4))
      return write_error();
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(filename.c_str(), img.imsize().x, img.imsize().y, 4,
            img.data(), 75))
      return write_error();
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(
            filename.c_str(), img.imsize().x, img.imsize().y, 4, img.data()))
      return write_error();
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(
            filename.c_str(), img.imsize().x, img.imsize().y, 4, img.data()))
      return write_error();
    return true;
  } else if (is_hdr_filename(filename)) {
    return save_image(filename, srgb_to_rgb(img), error);
  } else {
    return format_error();
  }
}

// Loads an hdr image.
[[nodiscard]] bool load_image(
    const string& filename, image<vec3f>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0)
      return read_error();
    if (!pixels) return read_error();
    auto cpixels = std::unique_ptr<float[]>(new float[width * height * 3]);
    for (auto i = (size_t)0; i < (size_t)width * (size_t)height; i++) {
      cpixels[i * 3 + 0] = pixels[i * 4 + 0];
      cpixels[i * 3 + 1] = pixels[i * 4 + 1];
      cpixels[i * 3 + 2] = pixels[i * 4 + 2];
    }
    img = image{{width, height}, (const vec3f*)cpixels.get()};
    free(pixels);
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 3);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec3f*)pixels};
    delete[] pixels;
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 3);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec3f*)pixels};
    free(pixels);
    return true;
  } else if (!is_hdr_filename(filename)) {
    auto imgb = image<vec3b>{};
    if (!load_image(filename, imgb, error)) return false;
    img = srgb_to_rgb(imgb);
    return true;
  } else {
    return format_error();
  }
}

// Saves an hdr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<vec3f>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    if (!stbi_write_hdr(filename.c_str(), img.imsize().x, img.imsize().y, 3,
            (float*)img.data()))
      return write_error();
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), img.imsize().x, img.imsize().y, 3,
            (float*)img.data()))
      return write_error();
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((float*)img.data(), img.imsize().x, img.imsize().y, 3,
            filename.c_str()) < 0)
      return write_error();
    return true;
  } else if (!is_hdr_filename(filename)) {
    return save_image(filename, rgb_to_srgbb(img), error);
  } else {
    return format_error();
  }
}

// Loads an ldr image.
[[nodiscard]] bool load_image(
    const string& filename, image<vec3b>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
      ext == ".tga" || ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 3);
    if (!pixels) return read_error();
    img = image{{width, height}, (const vec3b*)pixels};
    free(pixels);
    return true;
  } else if (is_hdr_filename(filename)) {
    auto imgf = image<vec3f>{};
    if (!load_image(filename, imgf, error)) return false;
    img = rgb_to_srgbb(imgf);
    return true;
  } else {
    return format_error();
  }
}

// Saves an ldr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<vec3b>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), img.imsize().x, img.imsize().y, 3,
            img.data(), img.imsize().x * 3))
      return write_error();
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(filename.c_str(), img.imsize().x, img.imsize().y, 3,
            img.data(), 75))
      return write_error();
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(
            filename.c_str(), img.imsize().x, img.imsize().y, 3, img.data()))
      return write_error();
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(
            filename.c_str(), img.imsize().x, img.imsize().y, 3, img.data()))
      return write_error();
    return true;
  } else if (is_hdr_filename(filename)) {
    return save_image(filename, srgb_to_rgb(img), error);
  } else {
    return format_error();
  }
}

// Loads an hdr image.
[[nodiscard]] bool load_image(
    const string& filename, image<float>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".exr" || ext == ".EXR") {
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0)
      return read_error();
    if (!pixels) return read_error();
    auto cpixels = std::unique_ptr<float[]>(new float[width * height]);
    for (auto i = (size_t)0; i < (size_t)width * (size_t)height; i++) {
      cpixels[i + 0] =
          (pixels[i * 4 + 0] + pixels[i * 4 + 1] + pixels[i * 4 + 2]) / 3;
    }
    img = image{{width, height}, (const float*)cpixels.get()};
    free(pixels);
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 1);
    if (!pixels) return read_error();
    img = image{{width, height}, (const float*)pixels};
    delete[] pixels;
    return true;
  } else if (ext == ".hdr" || ext == ".HDR") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 1);
    if (!pixels) return read_error();
    img = image{{width, height}, (const float*)pixels};
    free(pixels);
    return true;
  } else if (!is_hdr_filename(filename)) {
    auto imgb = image<byte>{};
    if (!load_image(filename, imgb, error)) return false;
    img = srgb_to_rgb(imgb);
    return true;
  } else {
    return format_error();
  }
}

// Saves an hdr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<float>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".hdr" || ext == ".HDR") {
    auto npixels = (size_t)img.imsize().x * (size_t)img.imsize().y;
    auto cpixels = std::unique_ptr<float[]>(new float[npixels * 3]);
    for (auto idx = (size_t)0; idx < npixels; idx++) {
      cpixels[idx * 3 + 0] = *(img.begin() + idx);
      cpixels[idx * 3 + 1] = *(img.begin() + idx);
      cpixels[idx * 3 + 2] = *(img.begin() + idx);
    }
    if (!stbi_write_hdr(filename.c_str(), img.imsize().x, img.imsize().y, 3,
            (float*)cpixels.get()))
      return write_error();
    return true;
  } else if (ext == ".pfm" || ext == ".PFM") {
    if (!save_pfm(filename.c_str(), img.imsize().x, img.imsize().y, 1,
            (float*)img.data()))
      return write_error();
    return true;
  } else if (ext == ".exr" || ext == ".EXR") {
    if (SaveEXR((float*)img.data(), img.imsize().x, img.imsize().y, 1,
            filename.c_str()) < 0)
      return write_error();
    return true;
  } else if (!is_hdr_filename(filename)) {
    return save_image(filename, rgb_to_srgbb(img), error);
  } else {
    return format_error();
  }
}

// Loads an ldr image.
[[nodiscard]] bool load_image(
    const string& filename, image<byte>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
      ext == ".tga" || ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 1);
    if (!pixels) return read_error();
    img = image{{width, height}, (const byte*)pixels};
    free(pixels);
    return true;
  } else if (is_hdr_filename(filename)) {
    auto imgf = image<float>{};
    if (!load_image(filename, imgf, error)) return false;
    img = rgb_to_srgbb(imgf);
    return true;
  } else {
    return format_error();
  }
}

// Saves an ldr image.
[[nodiscard]] bool save_image(
    const string& filename, const image<byte>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG") {
    if (!stbi_write_png(filename.c_str(), img.imsize().x, img.imsize().y, 1,
            img.data(), img.imsize().x * 1))
      return write_error();
    return true;
  } else if (ext == ".jpg" || ext == ".JPG") {
    if (!stbi_write_jpg(filename.c_str(), img.imsize().x, img.imsize().y, 1,
            img.data(), 75))
      return write_error();
    return true;
  } else if (ext == ".tga" || ext == ".TGA") {
    if (!stbi_write_tga(
            filename.c_str(), img.imsize().x, img.imsize().y, 1, img.data()))
      return write_error();
    return true;
  } else if (ext == ".bmp" || ext == ".BMP") {
    if (!stbi_write_bmp(
            filename.c_str(), img.imsize().x, img.imsize().y, 1, img.data()))
      return write_error();
    return true;
  } else if (is_hdr_filename(filename)) {
    return save_image(filename, srgb_to_rgb(img), error);
  } else {
    return format_error();
  }
}

// Loads a 16 bit image.
[[nodiscard]] bool load_image(
    const string& filename, image<ushort>& img, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".png" || ext == ".PNG" || ext == ".jpg" || ext == ".JPG" ||
      ext == ".tga" || ext == ".TGA" || ext == ".bmp" || ext == ".BMP") {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_16(filename.c_str(), &width, &height, &ncomp, 1);
    if (!pixels) return read_error();
    img = image{{width, height}, (const uint16_t*)pixels};
    free(pixels);
    return true;
  } else if (is_hdr_filename(filename)) {
    auto imgf = image<float>{};
    if (!load_image(filename, imgf, error)) return false;
    img = float_to_ushort(imgf);
    return true;
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Volume load
static float* load_yvol(
    const char* filename, int* w, int* h, int* d, int* nc, int req) {
  auto fs = fopen_utf8(filename, "rb");
  if (!fs) return nullptr;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  // buffer
  char buffer[4096];
  auto toks = vector<string>();

  // read magic
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  if (toks[0] != "YVOL") return nullptr;

  // read w, h
  if (!fgets(buffer, sizeof(buffer), fs)) return nullptr;
  toks = split_string(buffer);
  *w   = atoi(toks[0].c_str());
  *h   = atoi(toks[1].c_str());
  *d   = atoi(toks[2].c_str());
  *nc  = atoi(toks[3].c_str());

  // read data
  auto nvoxels = (size_t)(*w) * (size_t)(*h) * (size_t)(*d);
  auto nvalues = nvoxels * (size_t)(*nc);
  auto voxels  = std::unique_ptr<float[]>(new float[nvalues]);
  if (fread(voxels.get(), sizeof(float), nvalues, fs) != nvalues)
    return nullptr;

  // proper number of channels
  if (!req || *nc == req) return voxels.release();

  // pack into channels
  if (req < 0 || req > 4) {
    return nullptr;
  }
  auto cvoxels = std::unique_ptr<float[]>(new float[req * nvoxels]);
  for (auto i = 0; i < nvoxels; i++) {
    auto vp = voxels.get() + i * (*nc);
    auto cp = cvoxels.get() + i * req;
    if (*nc == 1) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[0];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[0];
          cp[2] = vp[0];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[0];
          cp[2] = vp[0];
          cp[3] = 1;
          break;
      }
    } else if (*nc == 2) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
      }
    } else if (*nc == 3) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          cp[3] = 1;
          break;
      }
    } else if (*nc == 4) {
      switch (req) {
        case 1: cp[0] = vp[0]; break;
        case 2:
          cp[0] = vp[0];
          cp[1] = vp[1];
          break;
        case 3:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          break;
        case 4:
          cp[0] = vp[0];
          cp[1] = vp[1];
          cp[2] = vp[2];
          cp[3] = vp[3];
          break;
      }
    }
  }
  return cvoxels.release();
}

// save pfm
static bool save_yvol(
    const char* filename, int w, int h, int d, int nc, const float* voxels) {
  auto fs = fopen_utf8(filename, "wb");
  if (!fs) return false;
  auto fs_guard = std::unique_ptr<FILE, void (*)(FILE*)>{
      fs, [](FILE* f) { fclose(f); }};

  if (fprintf(fs, "YVOL\n") < 0) return false;
  if (fprintf(fs, "%d %d %d %d\n", w, h, d, nc) < 0) return false;
  auto nvalues = (size_t)w * (size_t)h * (size_t)d * (size_t)nc;
  if (fwrite(voxels, sizeof(float), nvalues, fs) != nvalues) return false;

  return true;
}

// Loads volume data from binary format.
bool load_volume(const string& filename, volume<float>& vol, string& error) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto width = 0, height = 0, depth = 0, ncomp = 0;
  auto voxels = load_yvol(filename.c_str(), &width, &height, &depth, &ncomp, 1);
  if (!voxels) return read_error();
  vol = volume{{width, height, depth}, (const float*)voxels};
  delete[] voxels;
  return true;
}

// Saves volume data in binary format.
bool save_volume(
    const string& filename, const volume<float>& vol, string& error) {
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };
  if (!save_yvol(filename.c_str(), vol.size().x, vol.size().y, vol.size().z, 1,
          vol.data()))
    return write_error();
  return true;
}

}  // namespace yocto
