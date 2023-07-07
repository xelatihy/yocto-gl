//
// # Yocto/Color: Color operations
//
// Yocto/Color provides basic color utilities for writing graphics applications.
// In particular, we support color conversion to/from linear rgb, srgb, hsv,
// xyz, byte to float color conversions, colormaps, and a few color
// manipulations like contrast and saturation.
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
//
// LICENSE for colormap code
//
// License CC0 (public domain)
//
// Code by Matt Zucker from https://www.shadertoy.com/view/WlfXRN
// Data fitted from https://github.com/BIDS/colormap/blob/master/colormaps.py
// (which is licensed CC0).
// Magma, Inferno, Plasma, Viridis are released under CC0 by Nathaniel J. Smith,
// Stefan van der Walt, and (in the case of Viridis) Eric Firing:
// https://github.com/BIDS/colormap/blob/master/colormaps.py.
//
//

#ifndef YOCTO_COLOR_H_
#define YOCTO_COLOR_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <stdexcept>
#include <utility>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifndef kernel
#ifdef __CUDACC__
#define kernel __device__
#else
#define kernel
#endif
#endif

// -----------------------------------------------------------------------------
// COLOR OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Extracting components
inline kernel const vec3f& rgb(const vec4f& color) {
  return (const vec3f&)color;
}
inline kernel const float& alpha(const vec4f& color) { return color[3]; }

// Conversion between floats and bytes
constexpr kernel byte float_to_byte(float a) {
  return (byte)clamp(int(a * 256), 0, 255);
}
constexpr kernel float byte_to_float(byte a) { return a / 255.0f; }
constexpr kernel vec3b float_to_byte(const vec3f& a) {
  return (vec3b)clamp(vec3i(a * 256), 0, 255);
}
constexpr kernel vec4b float_to_byte(const vec4f& a) {
  return (vec4b)clamp(vec4i(a * 256), 0, 255);
}
constexpr kernel vec3f byte_to_float(const vec3b& a) { return a / 255.0f; }
constexpr kernel vec4f byte_to_float(const vec4b& a) { return a / 255.0f; }

// Luminance
constexpr kernel float luminance(const vec3f& a) {
  return 0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z;
}

// sRGB non-linear curve
constexpr kernel float srgb_to_rgb(float srgb) {
  return (srgb <= 0.04045f) ? srgb / 12.92f
                            : pow((srgb + 0.055f) / (1 + 0.055f), 2.4f);
}

constexpr kernel float rgb_to_srgb(float rgb) {
  return (rgb <= 0.0031308f) ? 12.92f * rgb
                             : (1 + 0.055f) * pow(rgb, 1 / 2.4f) - 0.055f;
}
constexpr kernel vec3f srgb_to_rgb(const vec3f& srgb) {
  return {srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z)};
}
constexpr kernel vec4f srgb_to_rgb(const vec4f& srgb) {
  return {
      srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z), srgb.w};
}
constexpr kernel vec3f rgb_to_srgb(const vec3f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z)};
}
constexpr kernel vec4f rgb_to_srgb(const vec4f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z), rgb.w};
}

// sRGB non-linear curve
constexpr kernel float srgbb_to_rgb(byte srgb) {
  return srgb_to_rgb(byte_to_float(srgb));
}

constexpr kernel byte rgb_to_srgbb(float rgb) {
  return float_to_byte(srgb_to_rgb(rgb));
}
constexpr kernel vec3f srgbb_to_rgb(const vec3b& srgb) {
  return srgb_to_rgb(byte_to_float(srgb));
}
constexpr kernel vec4f srgbb_to_rgb(const vec4b& srgb) {
  return srgb_to_rgb(byte_to_float(srgb));
}
constexpr kernel vec3b rgb_to_srgbb(const vec3f& rgb) {
  return float_to_byte(rgb_to_srgb(rgb));
}
constexpr kernel vec4b rgb_to_srgbb(const vec4f& rgb) {
  return float_to_byte(rgb_to_srgb(rgb));
}

// Conversion between number of channels.
constexpr kernel vec4f rgb_to_rgba(const vec3f& rgb) { return {rgb, 1.0f}; }
constexpr kernel vec4b rgb_to_rgba(const vec3b& rgb) {
  return {rgb, (byte)255};
}
constexpr kernel vec3f rgba_to_rgb(const vec4f& rgba) { return xyz(rgba); }

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
constexpr kernel vec3f lincontrast(
    const vec3f& rgb, float contrast, float grey) {
  return max(grey + (rgb - grey) * (contrast * 2), 0);
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
constexpr kernel vec3f logcontrast(
    const vec3f& rgb, float logcontrast, float grey) {
  auto epsilon  = 0.0001f;
  auto log_grey = log2(grey);
  auto log_ldr  = log2(rgb + epsilon);
  auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
  return max(exp2(adjusted) - epsilon, 0);
}
// Apply an s-shaped contrast.
constexpr kernel vec3f contrast(const vec3f& rgb, float contrast) {
  return gain(rgb, 1 - contrast);
}
// Apply saturation.
constexpr kernel vec3f saturate(const vec3f& rgb, float saturation,
    const vec3f& weights = vec3f{1.0f / 3.0f, 1.0f / 3.0f, 1.0f / 3.0f}) {
  auto grey = dot(weights, rgb);
  return max(grey + (rgb - grey) * (saturation * 2), 0);
}

// Convert channels
template <size_t M, typename T, size_t N>
constexpr kernel vec<T, M> convert_channels(const vec<T, N>& color) {
  constexpr auto a1 = std::is_same_v<T, byte> ? (byte)255 : 1.0f;
  if constexpr (N == M) {
    return color;
  } else if constexpr (N == 1) {
    if constexpr (M == 1) return {color[0]};
    if constexpr (M == 2) return {color[0], a1};
    if constexpr (M == 3) return {color[0], color[0], color[0]};
    if constexpr (M == 4) return {color[0], color[0], color[0], a1};
  } else if constexpr (N == 2) {
    if constexpr (M == 1) return {color[0]};
    if constexpr (M == 2) return {color[0], color[1]};
    if constexpr (M == 3) return {color[0], color[0], color[0]};
    if constexpr (M == 4) return {color[0], color[0], color[0], color[1]};
  } else if constexpr (N == 3) {
    if constexpr (M == 1) return {mean(color[0])};
    if constexpr (M == 2) return {mean(color[0]), a1};
    if constexpr (M == 3) return {color[0], color[1], color[2]};
    if constexpr (M == 4) return {color[0], color[1], color[2], a1};
  } else if constexpr (N == 4) {
    if constexpr (M == 1) return {mean(color[0])};
    if constexpr (M == 2) return {mean(color[0]), color[3]};
    if constexpr (M == 3) return {color[0], color[1], color[2]};
    if constexpr (M == 4) return {color[0], color[1], color[2], color[3]};
  } else {
    return vec<T, M>{0};
  }
}

#ifndef __CUDACC__

// Filmic tonemapping
constexpr kernel vec3f tonemap_filmic(
    const vec3f& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * 0.6f;  // brings it back to ACES range
    auto ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
               (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
    return max(vec3f{0, 0, 0}, ldr);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    constexpr auto ACESInputMat = transpose(mat3f{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    constexpr auto ACESOutputMat = transpose(mat3f{
        {1.60475, -0.53108, -0.07367},
        {-0.10208, 1.10813, -0.00605},
        {-0.00327, -0.07276, 1.07602},
    });
    // RRT => ODT
    auto RRTAndODTFit = [](const vec3f& v) -> vec3f {
      return (v * v + v * 0.0245786f - 0.000090537f) /
             (v * v * 0.983729f + v * 0.4329510f + 0.238081f);
    };

    auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
    return max(vec3f{0, 0, 0}, ldr);
  }
}

#else

// Filmic tonemapping
inline kernel vec3f tonemap_filmic(
    const vec3f& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * 0.6f;  // brings it back to ACES range
    auto ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
               (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
    return max(vec3f{0, 0, 0}, ldr);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    auto ACESInputMat = transpose(mat3f{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    auto ACESOutputMat = transpose(mat3f{
        {1.60475, -0.53108, -0.07367},
        {-0.10208, 1.10813, -0.00605},
        {-0.00327, -0.07276, 1.07602},
    });
    // RRT => ODT
    auto RRTAndODTFit = [](const vec3f& v) -> vec3f {
      return (v * v + v * 0.0245786f - 0.000090537f) /
             (v * v * 0.983729f + v * 0.4329510f + 0.238081f);
    };

    auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
    return max(vec3f{0, 0, 0}, ldr);
  }
}

#endif

// Tonemap
constexpr kernel vec3f tonemap(
    const vec3f& hdr, float exposure, bool filmic, bool srgb = true) {
  auto rgb = hdr;
  if (exposure != 0) rgb *= exp2(exposure);
  if (filmic) rgb = tonemap_filmic(rgb);
  if (srgb) rgb = rgb_to_srgb(rgb);
  return rgb;
}

// Tonemap
constexpr kernel vec4f tonemap(
    const vec4f& hdr, float exposure, bool filmic, bool srgb = true) {
  auto ldr = tonemap(xyz(hdr), exposure, filmic, srgb);
  return {ldr, alpha(hdr)};
}

// Composite colors
constexpr kernel vec4f composite(const vec4f& a, const vec4f& b) {
  if (a.w == 0 && b.w == 0) return {0, 0, 0, 0};
  auto cc = xyz(a) * alpha(a) + xyz(b) * alpha(b) * (1 - alpha(a));
  auto ca = alpha(a) + alpha(b) * (1 - alpha(a));
  return {cc / ca, ca};
}

// Convert between CIE XYZ and RGB
constexpr kernel vec3f rgb_to_xyz(const vec3f& rgb) {
  // https://en.wikipedia.org/wiki/SRGB
  constexpr auto m = mat3f{
      {0.4124f, 0.2126f, 0.0193f},
      {0.3576f, 0.7152f, 0.1192f},
      {0.1805f, 0.0722f, 0.9504f},
  };
  return m * rgb;
}

constexpr kernel vec3f xyz_to_rgb(const vec3f& xyz) {
  // https://en.wikipedia.org/wiki/SRGB
  constexpr auto m = mat3f{
      {+3.2406f, -0.9689f, +0.0557f},
      {-1.5372f, +1.8758f, -0.2040f},
      {-0.4986f, +0.0415f, +1.0570f},
  };
  return m * xyz;
}

// Convert between CIE XYZ and xyY

constexpr kernel vec3f xyz_to_xyY(const vec3f& xyz) {
  if (xyz == 0) return {0, 0, 0};
  auto [x, y, z] = xyz;
  return {x / (x + y + z), y / (x + y + z), y};
}

constexpr kernel vec3f xyY_to_xyz(const vec3f& xyY) {
  auto [x, y, Y] = xyY;
  if (y == 0) return {0, 0, 0};
  return {x * Y / y, Y, (1 - x - y) * Y / y};
}

// Convert HSV to RGB
constexpr kernel vec3f hsv_to_rgb(const vec3f& hsv) {
  // from Imgui.cpp
  auto [h, s, v] = hsv;
  if (s == 0) return {v, v, v};

  h      = fmod(h, 1.0f) / (60.0f / 360.0f);
  auto i = (int)h;
  auto f = h - (float)i;
  auto p = v * (1 - s);
  auto q = v * (1 - s * f);
  auto t = v * (1 - s * (1 - f));

  switch (i) {
    case 0: return {v, t, p};
    case 1: return {q, v, p};
    case 2: return {p, v, t};
    case 3: return {p, q, v};
    case 4: return {t, p, v};
    case 5: return {v, p, q};
    default: return {};
  }
}

constexpr kernel vec3f rgb_to_hsv(const vec3f& rgb) {
  // from Imgui.cpp
  auto [r, g, b] = rgb;
  auto K         = 0.0f;
  if (g < b) {
    std::swap(g, b);
    K = -1;
  }
  if (r < g) {
    std::swap(r, g);
    K = -2 / 6.0f - K;
  }

  auto c = r - (g < b ? g : b);  // chroma
  return {abs(K + (g - b) / (6 * c + 1e-20f)), c / (r + 1e-20f), r};
}

// Approximate color of blackbody radiation from wavelength in nm.
constexpr kernel vec3f blackbody_to_rgb(float temperature) {
  // clamp to valid range
  auto t = clamp(temperature, 1667.0f, 25000.0f) / 1000.0f;
  // compute x
  auto x = 0.0f;
  if (temperature < 4000.0f) {
    x = -0.2661239f * 1 / (t * t * t) - 0.2343589f * 1 / (t * t) +
        0.8776956f * (1 / t) + 0.179910f;
  } else {
    x = -3.0258469f * 1 / (t * t * t) + 2.1070379f * 1 / (t * t) +
        0.2226347f * (1 / t) + 0.240390f;
  }
  // compute y
  auto y = 0.0f;
  if (temperature < 2222.0f) {
    y = -1.1063814f * (x * x * x) - 1.34811020f * (x * x) + 2.18555832f * x -
        0.20219683f;
  } else if (temperature < 4000.0f) {
    y = -0.9549476f * (x * x * x) - 1.37418593f * (x * x) + 2.09137015f * x -
        0.16748867f;
  } else {
    y = 3.0817580f * (x * x * x) - 5.87338670f * (x * x) + 3.75112997f * x -
        0.37001483f;
  }
  return xyz_to_rgb(xyY_to_xyz(vec3f{x, y, 1}));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR MAPS
// -----------------------------------------------------------------------------
namespace yocto {

constexpr kernel vec3f colormap_viridis(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec3f{
      0.2777273272234177f, 0.005407344544966578f, 0.3340998053353061f};
  constexpr auto c1 = vec3f{
      0.1050930431085774f, 1.404613529898575f, 1.384590162594685f};
  constexpr auto c2 = vec3f{
      -0.3308618287255563f, 0.214847559468213f, 0.09509516302823659f};
  constexpr auto c3 = vec3f{
      -4.634230498983486f, -5.799100973351585f, -19.33244095627987f};
  constexpr auto c4 = vec3f{
      6.228269936347081f, 14.17993336680509f, 56.69055260068105f};
  constexpr auto c5 = vec3f{
      4.776384997670288f, -13.74514537774601f, -65.35303263337234f};
  constexpr auto c6 = vec3f{
      -5.435455855934631f, 4.645852612178535f, 26.3124352495832f};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

constexpr kernel vec3f colormap_plasma(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec3f{
      0.05873234392399702f, 0.02333670892565664f, 0.5433401826748754f};
  constexpr auto c1 = vec3f{
      2.176514634195958f, 0.2383834171260182f, 0.7539604599784036f};
  constexpr auto c2 = vec3f{
      -2.689460476458034f, -7.455851135738909f, 3.110799939717086f};
  constexpr auto c3 = vec3f{
      6.130348345893603f, 42.3461881477227f, -28.51885465332158f};
  constexpr auto c4 = vec3f{
      -11.10743619062271f, -82.66631109428045f, 60.13984767418263f};
  constexpr auto c5 = vec3f{
      10.02306557647065f, 71.41361770095349f, -54.07218655560067f};
  constexpr auto c6 = vec3f{
      -3.658713842777788f, -22.93153465461149f, 18.19190778539828f};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

constexpr kernel vec3f colormap_magma(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec3f{
      -0.002136485053939582f, -0.000749655052795221f, -0.005386127855323933f};
  constexpr auto c1 = vec3f{
      0.2516605407371642f, 0.6775232436837668f, 2.494026599312351f};
  constexpr auto c2 = vec3f{
      8.353717279216625f, -3.577719514958484f, 0.3144679030132573f};
  constexpr auto c3 = vec3f{
      -27.66873308576866f, 14.26473078096533f, -13.64921318813922f};
  constexpr auto c4 = vec3f{
      52.17613981234068f, -27.94360607168351f, 12.94416944238394f};
  constexpr auto c5 = vec3f{
      -50.76852536473588f, 29.04658282127291f, 4.23415299384598f};
  constexpr auto c6 = vec3f{
      18.65570506591883f, -11.48977351997711f, -5.601961508734096f};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

constexpr kernel vec3f colormap_inferno(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec3f{
      0.0002189403691192265f, 0.001651004631001012f, -0.01948089843709184f};
  constexpr auto c1 = vec3f{
      0.1065134194856116f, 0.5639564367884091, 3.932712388889277f};
  constexpr auto c2 = vec3f{
      11.60249308247187f, -3.972853965665698, -15.9423941062914f};
  constexpr auto c3 = vec3f{
      -41.70399613139459f, 17.43639888205313, 44.35414519872813f};
  constexpr auto c4 = vec3f{
      77.162935699427f, -33.40235894210092, -81.80730925738993f};
  constexpr auto c5 = vec3f{
      -71.31942824499214f, 32.62606426397723f, 73.20951985803202f};
  constexpr auto c6 = vec3f{
      25.13112622477341f, -12.24266895238567f, -23.07032500287172f};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

// Colormap type
enum struct colormap_type { viridis, plasma, magma, inferno };

// Colormaps from {0,1] to color

constexpr kernel vec3f colormap(float t, colormap_type type) {
  t = clamp(t, 0.0f, 1.0f);
  switch (type) {
    case colormap_type::viridis: return colormap_viridis(t);
    case colormap_type::magma: return colormap_magma(t);
    case colormap_type::inferno: return colormap_inferno(t);
    case colormap_type::plasma: return colormap_plasma(t);
  }
  return {0, 0, 0};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR GRADING
// -----------------------------------------------------------------------------
namespace yocto {

// minimal color grading
struct colorgrade_params {
  float exposure         = 0;
  vec3f tint             = {1, 1, 1};
  float lincontrast      = 0.5;
  float logcontrast      = 0.5;
  float linsaturation    = 0.5;
  bool  filmic           = false;
  bool  srgb             = true;
  float contrast         = 0.5;
  float saturation       = 0.5;
  float shadows          = 0.5;
  float midtones         = 0.5;
  float highlights       = 0.5;
  vec3f shadows_color    = {1, 1, 1};
  vec3f midtones_color   = {1, 1, 1};
  vec3f highlights_color = {1, 1, 1};
};

constexpr kernel vec3f colorgrade(
    const vec3f& rgb_, bool linear, const colorgrade_params& params) {
  auto rgb = rgb_;
  if (params.exposure != 0) rgb *= exp2(params.exposure);
  if (params.tint != 1) rgb *= params.tint;
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
      params.highlights != 0.5f || params.shadows_color != 1 ||
      params.midtones_color != 1 || params.highlights_color != 1) {
    auto lift  = params.shadows_color;
    auto gamma = params.midtones_color;
    auto gain  = params.highlights_color;

    lift      = lift - mean(lift) + params.shadows - 0.5f;
    gain      = gain - mean(gain) + params.highlights + 0.5f;
    auto grey = gamma - mean(gamma) + params.midtones;
    gamma     = log((0.5f - lift) / (gain - lift)) / log(grey);

    // apply_image
    auto lerp_value = clamp(pow(rgb, 1 / gamma), 0.0f, 1.0f);
    rgb             = gain * lerp_value + lift * (1 - lerp_value);
  }
  return rgb;
}

constexpr kernel vec4f colorgrade(
    const vec4f& rgba, bool linear, const colorgrade_params& params) {
  auto graded = colorgrade(xyz(rgba), linear, params);
  return {graded, alpha(rgba)};
}

}  // namespace yocto

#ifndef __CUDACC__

// -----------------------------------------------------------------------------
// COLOR SPACE CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// RGB color spaces
enum struct color_space {
  rgb,         // default linear space (srgb linear)
  srgb,        // srgb color space (non-linear)
  adobe,       // Adobe rgb color space (non-linear)
  prophoto,    // ProPhoto Kodak rgb color space (non-linear)
  rec709,      // hdtv color space (non-linear)
  rec2020,     // uhtv color space (non-linear)
  rec2100pq,   // hdr color space with perceptual quantizer (non-linear)
  rec2100hlg,  // hdr color space with hybrid log gamma (non-linear)
  aces2065,    // ACES storage format (linear)
  acescg,      // ACES CG computation (linear)
  acescc,      // ACES color correction (non-linear)
  acescct,     // ACES color correction 2 (non-linear)
  p3dci,       // P3 DCI (non-linear)
  p3d60,       // P3 variation for D60 (non-linear)
  p3d65,       // P3 variation for D65 (non-linear)
  p3display,   // Apple display P3
};

// Conversion between rgb color spaces

constexpr vec3f color_to_xyz(const vec3f& col, color_space from);

constexpr vec3f xyz_to_color(const vec3f& xyz, color_space to);

// Conversion between rgb color spaces

constexpr vec3f convert_color(
    const vec3f& col, color_space from, color_space to);

}  // namespace yocto

#endif

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#ifndef __CUDACC__

// -----------------------------------------------------------------------------
// COLOR SPACES
// -----------------------------------------------------------------------------
namespace yocto {

// RGB color space definition. Various predefined color spaces are listed below.

struct color_space_params {
  // Curve type
  enum struct curve_t {
    // clang-format off
    linear, gamma, linear_gamma, aces_cc, aces_cct, pq, hlg
    // clang-format on
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
// Input: red, green, blue, white (x,y) chromaticities
// Algorithm from: SMPTE Recommended Practice RP 177-1993
// http://car.france3.mars.free.fr/HD/INA-%2026%20jan%2006/SMPTE%20normes%20et%20confs/rp177.pdf

constexpr mat3f rgb_to_xyz_mat(
    const vec2f& rc, const vec2f& gc, const vec2f& bc, const vec2f& wc) {
  auto r = vec3f{rc, 1 - sum(rc)}, g = vec3f{gc, 1 - sum(gc)},
       b = vec3f{bc, 1 - sum(bc)}, w = vec3f{wc, 1 - sum(wc)};
  auto [cr, cg, cb] = inverse(mat3f{r, g, b}) * w / w.y;
  return {cr * r, cg * g, cb * b};
}

// Construct an RGB color space. Predefined color spaces below
constexpr color_space_params get_color_scape_params(color_space space) {
  constexpr auto make_linear_rgb_space =
      [](const vec2f& red, const vec2f& green, const vec2f& blue,
          const vec2f& white) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)),
            color_space_params::curve_t::linear};
      };
  constexpr auto make_gamma_rgb_space =
      [](const vec2f& red, const vec2f& green, const vec2f& blue,
          const vec2f& white, float gamma,
          const vec4f& curve_abcd = vec4f{0, 0, 0, 0}) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)),
            curve_abcd == vec4f{0, 0, 0, 0}
                ? color_space_params::curve_t::gamma
                : color_space_params::curve_t::linear_gamma};
      };
  constexpr auto make_other_rgb_space =
      [](const vec2f& red, const vec2f& green, const vec2f& blue,
          const vec2f& white, color_space_params::curve_t curve_type) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)), curve_type};
      };

  // color space parameters
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto rgb_params = make_linear_rgb_space({0.6400f, 0.3300f},
      {0.3000f, 0.6000f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f});
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto srgb_params = make_gamma_rgb_space({0.6400f, 0.3300f},
      {0.3000f, 0.6000f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f}, 2.4f,
      {1.055f, 0.055f, 12.92f, 0.0031308f});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_System
  constexpr auto aces2065_params = make_linear_rgb_space({0.7347f, 0.2653f},
      {0.0000f, 1.0000f}, {0.0001f, -0.0770f}, {0.32168f, 0.33767f});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescg_params = make_linear_rgb_space({0.7130f, 0.2930f},
      {0.1650f, 0.8300f}, {0.1280f, 0.0440f}, {0.32168f, 0.33767f});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescc_params = make_other_rgb_space({0.7130f, 0.2930f},
      {0.1650f, 0.8300f}, {0.1280f, 0.0440f}, {0.32168f, 0.33767f},
      color_space_params::curve_t::aces_cc);
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescct_params = make_other_rgb_space({0.7130f, 0.2930f},
      {0.1650f, 0.8300f}, {0.1280f, 0.0440f}, {0.32168f, 0.33767f},
      color_space_params::curve_t::aces_cct);
  // https://en.wikipedia.org/wiki/Adobe_RGB_color_space
  constexpr auto adobe_params = make_gamma_rgb_space({0.6400f, 0.3300f},
      {0.2100f, 0.7100f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f}, 2.19921875f);
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto rec709_params = make_gamma_rgb_space({0.6400f, 0.3300f},
      {0.3000f, 0.6000f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f}, 1 / 0.45f,
      {1.099f, 0.099f, 4.500f, 0.018f});
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2020_params = make_gamma_rgb_space({0.7080f, 0.2920f},
      {0.1700f, 0.7970f}, {0.1310f, 0.0460f}, {0.3127f, 0.3290f}, 1 / 0.45f,
      {1.09929682680944f, 0.09929682680944f, 4.5f, 0.018053968510807f});
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2100pq_params = make_other_rgb_space({0.7080f, 0.2920f},
      {0.1700f, 0.7970f}, {0.1310f, 0.0460f}, {0.3127f, 0.3290f},
      color_space_params::curve_t::pq);
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2100hlg_params = make_other_rgb_space({0.7080f, 0.2920f},
      {0.1700f, 0.7970f}, {0.1310f, 0.0460f}, {0.3127f, 0.3290f},
      color_space_params::curve_t::hlg);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3dci_params = make_gamma_rgb_space({0.6800f, 0.3200f},
      {0.2650f, 0.6900f}, {0.1500f, 0.0600f}, {0.3140f, 0.3510f}, 1.6f);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3d60_params = make_gamma_rgb_space({0.6800f, 0.3200f},
      {0.2650f, 0.6900f}, {0.1500f, 0.0600f}, {0.32168f, 0.33767f}, 1.6f);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3d65_params = make_gamma_rgb_space({0.6800f, 0.3200f},
      {0.2650f, 0.6900f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f}, 1.6f);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3display_params = make_gamma_rgb_space({0.6800f, 0.3200f},
      {0.2650f, 0.6900f}, {0.1500f, 0.0600f}, {0.3127f, 0.3290f}, 2.4f,
      {1.055f, 0.055f, 12.92f, 0.0031308f});
  // https://en.wikipedia.org/wiki/ProPhoto_RGB_color_space
  constexpr auto prophoto_params = make_gamma_rgb_space({0.7347f, 0.2653f},
      {0.1596f, 0.8404f}, {0.0366f, 0.0001f}, {0.3457f, 0.3585f}, 1.8f,
      {1.0f, 0.0f, 16.0f, 0.001953125f});

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

constexpr float gamma_display_to_linear(float x, float gamma) {
  return pow(x, gamma);
}

constexpr float gamma_linear_to_display(float x, float gamma) {
  return pow(x, 1 / gamma);
}

constexpr vec3f gamma_display_to_linear(const vec3f& rgb, float gamma) {
  return {gamma_display_to_linear(rgb.x, gamma),
      gamma_display_to_linear(rgb.y, gamma),
      gamma_display_to_linear(rgb.z, gamma)};
}

constexpr vec3f gamma_linear_to_display(const vec3f& rgb, float gamma) {
  return {gamma_linear_to_display(rgb.x, gamma),
      gamma_linear_to_display(rgb.y, gamma),
      gamma_linear_to_display(rgb.z, gamma)};
}

// https://en.wikipedia.org/wiki/Rec._709

constexpr float gamma_display_to_linear(
    float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  return (x < 1 / d) ? (x / c) : pow((x + b) / a, gamma);
}

constexpr float gamma_linear_to_display(
    float x, float gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  return (x < d) ? (x * c) : (a * pow(x, 1 / gamma) - b);
}

constexpr vec3f gamma_display_to_linear(
    const vec3f& rgb, float gamma, const vec4f& abcd) {
  return {gamma_display_to_linear(rgb.x, gamma, abcd),
      gamma_display_to_linear(rgb.y, gamma, abcd),
      gamma_display_to_linear(rgb.z, gamma, abcd)};
}

constexpr vec3f gamma_linear_to_display(
    const vec3f& rgb, float gamma, const vec4f& abcd) {
  return {gamma_linear_to_display(rgb.x, gamma, abcd),
      gamma_linear_to_display(rgb.y, gamma, abcd),
      gamma_linear_to_display(rgb.z, gamma, abcd)};
}

// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx

constexpr float acescc_display_to_linear(float x) {
  if (x < -0.3013698630f) {  // (9.72-15)/17.52
    return (exp2(x * 17.52f - 9.72f) - exp2(-16.0f)) * 2;
  } else if (x < (log2(65504.0f) + 9.72f) / 17.52f) {
    return exp2(x * 17.52f - 9.72f);
  } else {  // (in >= (log2(65504)+9.72)/17.52)
    return 65504.0f;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx

constexpr float acescct_display_to_linear(float x) {
  if (x < 0.155251141552511f) {
    return (x - 0.0729055341958355f) / 10.5402377416545f;
  } else {
    return exp2(x * 17.52f - 9.72f);
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx

constexpr float acescc_linear_to_display(float x) {
  if (x <= 0) {
    return -0.3584474886f;  // =(log2( pow(2.,-16.))+9.72)/17.52
  } else if (x < exp2(-15.0f)) {
    return (log2(exp2(-16.0f) + x * 0.5f) + 9.72f) / 17.52f;
  } else {  // (in >= pow(2.,-15))
    return (log2(x) + 9.72f) / 17.52f;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx

constexpr float acescct_linear_to_display(float x) {
  if (x <= 0.0078125f) {
    return 10.5402377416545f * x + 0.0729055341958355f;
  } else {
    return (log2(x) + 9.72f) / 17.52f;
  }
}

constexpr vec3f acescc_display_to_linear(const vec3f& rgb) {
  return {acescc_display_to_linear(rgb.x), acescc_display_to_linear(rgb.y),
      acescc_display_to_linear(rgb.z)};
}

constexpr vec3f acescct_display_to_linear(const vec3f& rgb) {
  return {acescct_display_to_linear(rgb.x), acescct_display_to_linear(rgb.y),
      acescct_display_to_linear(rgb.z)};
}

constexpr vec3f acescc_linear_to_display(const vec3f& rgb) {
  return {acescc_linear_to_display(rgb.x), acescc_linear_to_display(rgb.y),
      acescc_linear_to_display(rgb.z)};
}

constexpr vec3f acescct_linear_to_display(const vec3f& rgb) {
  return {acescct_linear_to_display(rgb.x), acescct_linear_to_display(rgb.y),
      acescct_linear_to_display(rgb.z)};
}

// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// https://github.com/ampas/aces-dev/blob/master/transforms/ctl/lib/ACESlib.Utilities_Color.ctl
// In PQ, we assume that the linear luminance in [0,1] corresponds to
// [0,10000] cd m^2

constexpr float pq_display_to_linear(float x) {
  auto Np = pow(x, 1 / 78.84375f);
  auto L  = max(Np - 0.8359375f, 0.0f);
  L       = L / (18.8515625f - 18.6875f * Np);
  L       = pow(L, 1 / 0.1593017578125f);
  return L;
}

constexpr float pq_linear_to_display(float x) {
  return pow((0.8359375f + 18.8515625f * pow(x, 0.1593017578125f)) /
                 (1 + 18.6875f * pow(x, 0.1593017578125f)),
      78.84375f);
}

constexpr vec3f pq_display_to_linear(const vec3f& rgb) {
  return {pq_display_to_linear(rgb.x), pq_display_to_linear(rgb.y),
      pq_display_to_linear(rgb.z)};
}

constexpr vec3f pq_linear_to_display(const vec3f& rgb) {
  return {pq_linear_to_display(rgb.x), pq_linear_to_display(rgb.y),
      pq_linear_to_display(rgb.z)};
}

// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// In HLG, we assume that the linear luminance in [0,1] corresponds to
// [0,1000] cd m^2. Note that the version we report here is scaled in [0,1]
// range for nominal luminance. But HLG was initially defined in the [0,12]
// range where it maps 1 to 0.5 and 12 to 1. For use in HDR tonemapping that is
// likely a better range to use.

constexpr float hlg_display_to_linear(float x) {
  if (x < 0.5f) {
    return 3 * 3 * x * x;
  } else {
    return (exp((x - 0.55991073f) / 0.17883277f) + 0.28466892f) / 12;
  }
}

constexpr float hlg_linear_to_display(float x) {
  if (x < 1 / 12.0f) {
    return sqrt(3 * x);
  } else {
    return 0.17883277f * log(12 * x - 0.28466892f) + 0.55991073f;
  }
}

constexpr vec3f hlg_display_to_linear(const vec3f& rgb) {
  return {hlg_display_to_linear(rgb.x), hlg_display_to_linear(rgb.y),
      hlg_display_to_linear(rgb.z)};
}

constexpr vec3f hlg_linear_to_display(const vec3f& rgb) {
  return {hlg_linear_to_display(rgb.x), hlg_linear_to_display(rgb.y),
      hlg_linear_to_display(rgb.z)};
}

// Conversion to/from xyz
constexpr vec3f color_to_xyz(const vec3f& col, color_space from) {
  auto space = get_color_scape_params(from);
  auto rgb   = col;
  if (space.curve_type == color_space_params::curve_t::linear) {
    // do nothing
  } else if (space.curve_type == color_space_params::curve_t::gamma) {
    rgb = gamma_linear_to_display(rgb, space.curve_gamma);
  } else if (space.curve_type == color_space_params::curve_t::linear_gamma) {
    rgb = gamma_linear_to_display(rgb, space.curve_gamma, space.curve_abcd);
  } else if (space.curve_type == color_space_params::curve_t::aces_cc) {
    rgb = acescc_linear_to_display(rgb);
  } else if (space.curve_type == color_space_params::curve_t::aces_cct) {
    rgb = acescct_linear_to_display(rgb);
  } else if (space.curve_type == color_space_params::curve_t::pq) {
    rgb = pq_linear_to_display(rgb);
  } else if (space.curve_type == color_space_params::curve_t::hlg) {
    rgb = hlg_linear_to_display(rgb);
  } else {
    throw std::runtime_error{"should not have gotten here"};
  }
  return space.rgb_to_xyz_mat * rgb;
}

constexpr vec3f xyz_to_color(const vec3f& xyz, color_space to) {
  auto space = get_color_scape_params(to);
  auto rgb   = space.xyz_to_rgb_mat * xyz;
  if (space.curve_type == color_space_params::curve_t::linear) {
    // nothing
  } else if (space.curve_type == color_space_params::curve_t::gamma) {
    rgb = gamma_display_to_linear(rgb, space.curve_gamma);
  } else if (space.curve_type == color_space_params::curve_t::linear_gamma) {
    rgb = gamma_display_to_linear(rgb, space.curve_gamma, space.curve_abcd);
  } else if (space.curve_type == color_space_params::curve_t::aces_cc) {
    rgb = acescc_display_to_linear(rgb);
  } else if (space.curve_type == color_space_params::curve_t::aces_cct) {
    rgb = acescct_display_to_linear(rgb);
  } else if (space.curve_type == color_space_params::curve_t::pq) {
    rgb = pq_display_to_linear(rgb);
  } else if (space.curve_type == color_space_params::curve_t::hlg) {
    rgb = hlg_display_to_linear(rgb);
  } else {
    throw std::runtime_error{"should not have gotten here"};
  }
  return rgb;
}

constexpr vec3f convert_color(
    const vec3f& col, color_space from, color_space to) {
  if (from == to) return col;
  return xyz_to_color(color_to_xyz(col, from), to);
}

}  // namespace yocto

#endif

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
