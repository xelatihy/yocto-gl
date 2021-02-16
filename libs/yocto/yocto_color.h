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
//
// LICENSE for blackbody code
//
// Copyright (c) 2015 Neil Bartlett
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
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
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

#ifndef _YOCTO_COLOR_H_
#define _YOCTO_COLOR_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <stdexcept>
#include <utility>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// COLOR OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion between flots and bytes
inline vec4b float_to_byte(const vec4f& a);
inline vec4f byte_to_float(const vec4b& a);
inline byte  float_to_byte(float a);
inline float byte_to_float(byte a);

// Luminance
inline float luminance(const vec3f& a);

// sRGB non-linear curve
inline float srgb_to_rgb(float srgb);
inline float rgb_to_srgb(float rgb);
inline vec3f srgb_to_rgb(const vec3f& srgb);
inline vec4f srgb_to_rgb(const vec4f& srgb);
inline vec3f rgb_to_srgb(const vec3f& rgb);
inline vec4f rgb_to_srgb(const vec4f& rgb);

// Conversion between number of channels.
inline vec4f rgb_to_rgba(const vec3f& rgb);
inline vec3f rgba_to_rgb(const vec4f& rgba);

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f lincontrast(const vec3f& rgb, float contrast, float grey);
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f logcontrast(const vec3f& rgb, float logcontrast, float grey);
// Apply an s-shaped contrast.
inline vec3f contrast(const vec3f& rgb, float contrast);
// Apply saturation.
inline vec3f saturate(const vec3f& rgb, float saturation,
    const vec3f& weights = vec3f{0.333333, 0.333333, 0.333333});

// Apply tone mapping
inline vec3f tonemap(
    const vec3f& hdr, float exposure, bool filmic = false, bool srgb = true);
inline vec4f tonemap(
    const vec4f& hdr, float exposure, bool filmic = false, bool srgb = true);

// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb);
inline vec3f xyz_to_rgb(const vec3f& xyz);

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz);
inline vec3f xyY_to_xyz(const vec3f& xyY);

// Converts between HSV and RGB color spaces.
inline vec3f hsv_to_rgb(const vec3f& hsv);
inline vec3f rgb_to_hsv(const vec3f& rgb);

// Approximate color of blackbody radiation from wavelength in nm.
inline vec3f blackbody_to_rgb(float temperature);

// Colormap type
enum struct colormap_type { viridis, plasma, magma, inferno };

// Colormaps from [0,1] to color
inline vec3f colormap(float t, colormap_type type = colormap_type::viridis);

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

// Apply color grading from a linear or srgb color to an srgb color.
inline vec3f colorgrade(
    const vec3f& color, bool linear, const colorgrade_params& params);
inline vec4f colorgrade(
    const vec4f& color, bool linear, const colorgrade_params& params);

}  // namespace yocto

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
inline vec3f color_to_xyz(const vec3f& col, color_space from);
inline vec3f xyz_to_color(const vec3f& xyz, color_space to);

// Conversion between rgb color spaces
inline vec3f convert_color(const vec3f& col, color_space from, color_space to);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion between floats and bytes
inline vec4b float_to_byte(const vec4f& a) {
  return {(byte)clamp(int(a.x * 256), 0, 255),
      (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255),
      (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
  return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}
inline byte float_to_byte(float a) { return (byte)clamp(int(a * 256), 0, 255); }
inline float byte_to_float(byte a) { return a / 255.0f; }

// Luminance
inline float luminance(const vec3f& a) {
  return (0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z);
}

// sRGB non-linear curve
inline float srgb_to_rgb(float srgb) {
  return (srgb <= 0.04045) ? srgb / 12.92f
                           : pow((srgb + 0.055f) / (1.0f + 0.055f), 2.4f);
}
inline float rgb_to_srgb(float rgb) {
  return (rgb <= 0.0031308f) ? 12.92f * rgb
                             : (1 + 0.055f) * pow(rgb, 1 / 2.4f) - 0.055f;
}
inline vec3f srgb_to_rgb(const vec3f& srgb) {
  return {srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z)};
}
inline vec4f srgb_to_rgb(const vec4f& srgb) {
  return {
      srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z), srgb.w};
}
inline vec3f rgb_to_srgb(const vec3f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z)};
}
inline vec4f rgb_to_srgb(const vec4f& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z), rgb.w};
}

// Conversion between number of channels.
inline vec4f rgb_to_rgba(const vec3f& rgb) { return {rgb.x, rgb.y, rgb.z, 1}; }
inline vec3f rgba_to_rgb(const vec4f& rgba) { return xyz(rgba); }

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f lincontrast(const vec3f& rgb, float contrast, float grey) {
  return max(zero3f, grey + (rgb - grey) * (contrast * 2));
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
inline vec3f logcontrast(const vec3f& rgb, float logcontrast, float grey) {
  auto epsilon  = (float)0.0001;
  auto log_grey = log2(grey);
  auto log_ldr  = log2(rgb + epsilon);
  auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
  return max(zero3f, exp2(adjusted) - epsilon);
}
// Apply an s-shaped contrast.
inline vec3f contrast(const vec3f& rgb, float contrast) {
  return gain(rgb, 1 - contrast);
}
// Apply saturation.
inline vec3f saturate(
    const vec3f& rgb, float saturation, const vec3f& weights) {
  auto grey = dot(weights, rgb);
  return max(zero3f, grey + (rgb - grey) * (saturation * 2));
}

// Filmic tonemapping
inline vec3f tonemap_filmic(const vec3f& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * 0.6f;  // brings it back to ACES range
    auto ldr = (hdr * hdr * 2.51f + hdr * 0.03f) /
               (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
    return max(zero3f, ldr);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    static const auto ACESInputMat = transpose(mat3f{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    static const auto ACESOutputMat = transpose(mat3f{
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
    return max(zero3f, ldr);
  }
}

inline vec3f tonemap(const vec3f& hdr, float exposure, bool filmic, bool srgb) {
  auto rgb = hdr;
  if (exposure != 0) rgb *= exp2(exposure);
  if (filmic) rgb = tonemap_filmic(rgb);
  if (srgb) rgb = rgb_to_srgb(rgb);
  return rgb;
}
inline vec4f tonemap(const vec4f& hdr, float exposure, bool filmic, bool srgb) {
  auto ldr = tonemap(xyz(hdr), exposure, filmic, srgb);
  return {ldr.x, ldr.y, ldr.z, hdr.w};
}

// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat = mat3f{
      {0.4124, 0.2126, 0.0193},
      {0.3576, 0.7152, 0.1192},
      {0.1805, 0.0722, 0.9504},
  };
  return mat * rgb;
}
inline vec3f xyz_to_rgb(const vec3f& xyz) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat = mat3f{
      {+3.2406, -0.9689, +0.0557},
      {-1.5372, +1.8758, -0.2040},
      {-0.4986, +0.0415, +1.0570},
  };
  return mat * xyz;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
  if (xyz == zero3f) return zero3f;
  return {
      xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z), xyz.y};
}
inline vec3f xyY_to_xyz(const vec3f& xyY) {
  if (xyY.y == 0) return zero3f;
  return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}

// Convert HSV to RGB
inline vec3f hsv_to_rgb(const vec3f& hsv) {
  // from Imgui.cpp
  auto h = hsv.x, s = hsv.y, v = hsv.z;
  if (hsv.y == 0) return {v, v, v};

  h       = fmod(h, 1.0f) / (60.0f / 360.0f);
  int   i = (int)h;
  float f = h - (float)i;
  float p = v * (1 - s);
  float q = v * (1 - s * f);
  float t = v * (1 - s * (1 - f));

  switch (i) {
    case 0: return {v, t, p};
    case 1: return {q, v, p};
    case 2: return {p, v, t};
    case 3: return {p, q, v};
    case 4: return {t, p, v};
    case 5: return {v, p, q};
    default: return {v, p, q};
  }
}

inline vec3f rgb_to_hsv(const vec3f& rgb) {
  // from Imgui.cpp
  auto r = rgb.x, g = rgb.y, b = rgb.z;
  auto K = 0.f;
  if (g < b) {
    std::swap(g, b);
    K = -1;
  }
  if (r < g) {
    std::swap(r, g);
    K = -2 / 6.0f - K;
  }

  auto chroma = r - (g < b ? g : b);
  return {abs(K + (g - b) / (6 * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

// Approximate color of blackbody radiation from wavelength in nm.
inline vec3f blackbody_to_rgb(float temperature) {
  // https://github.com/neilbartlett/color-temperature
  auto rgb = zero3f;
  if ((temperature / 100) < 66) {
    rgb.x = 255;
  } else {
    // a + b x + c Log[x] /.
    // {a -> 351.97690566805693`,
    // b -> 0.114206453784165`,
    // c -> -40.25366309332127
    // x -> (kelvin/100) - 55}
    rgb.x = (temperature / 100) - 55;
    rgb.x = 351.97690566805693f + 0.114206453784165f * rgb.x -
            40.25366309332127f * log(rgb.x);
    if (rgb.x < 0) rgb.x = 0;
    if (rgb.x > 255) rgb.x = 255;
  }

  if ((temperature / 100) < 66) {
    // a + b x + c Log[x] /.
    // {a -> -155.25485562709179`,
    // b -> -0.44596950469579133`,
    // c -> 104.49216199393888`,
    // x -> (kelvin/100) - 2}
    rgb.y = (temperature / 100) - 2;
    rgb.y = -155.25485562709179f - 0.44596950469579133f * rgb.y +
            104.49216199393888f * log(rgb.y);
    if (rgb.y < 0) rgb.y = 0;
    if (rgb.y > 255) rgb.y = 255;
  } else {
    // a + b x + c Log[x] /.
    // {a -> 325.4494125711974`,
    // b -> 0.07943456536662342`,
    // c -> -28.0852963507957`,
    // x -> (kelvin/100) - 50}
    rgb.y = (temperature / 100) - 50;
    rgb.y = 325.4494125711974f + 0.07943456536662342f * rgb.y -
            28.0852963507957f * log(rgb.y);
    if (rgb.y < 0) rgb.y = 0;
    if (rgb.y > 255) rgb.y = 255;
  }

  if ((temperature / 100) >= 66) {
    rgb.z = 255;
  } else {
    if ((temperature / 100) <= 20) {
      rgb.z = 0;
    } else {
      // a + b x + c Log[x] /.
      // {a -> -254.76935184120902`,
      // b -> 0.8274096064007395`,
      // c -> 115.67994401066147`,
      // x -> kelvin/100 - 10}
      rgb.z = (temperature / 100) - 10;
      rgb.z = -254.76935184120902f + 0.8274096064007395f * rgb.z +
              115.67994401066147f * log(rgb.z);
      if (rgb.z < 0) rgb.z = 0;
      if (rgb.z > 255) rgb.z = 255;
    }
  }

  return srgb_to_rgb(rgb / 255);
}

inline vec3f colormap_viridis(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  static const auto c0 = vec3f{
      0.2777273272234177, 0.005407344544966578, 0.3340998053353061};
  static const auto c1 = vec3f{
      0.1050930431085774, 1.404613529898575, 1.384590162594685};
  static const auto c2 = vec3f{
      -0.3308618287255563, 0.214847559468213, 0.09509516302823659};
  static const auto c3 = vec3f{
      -4.634230498983486, -5.799100973351585, -19.33244095627987};
  static const auto c4 = vec3f{
      6.228269936347081, 14.17993336680509, 56.69055260068105};
  static const auto c5 = vec3f{
      4.776384997670288, -13.74514537774601, -65.35303263337234};
  static const auto c6 = vec3f{
      -5.435455855934631, 4.645852612178535, 26.3124352495832};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

inline vec3f colormap_plasma(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  static const auto c0 = vec3f{
      0.05873234392399702, 0.02333670892565664, 0.5433401826748754};
  static const auto c1 = vec3f{
      2.176514634195958, 0.2383834171260182, 0.7539604599784036};
  static const auto c2 = vec3f{
      -2.689460476458034, -7.455851135738909, 3.110799939717086};
  static const auto c3 = vec3f{
      6.130348345893603, 42.3461881477227, -28.51885465332158};
  static const auto c4 = vec3f{
      -11.10743619062271, -82.66631109428045, 60.13984767418263};
  static const auto c5 = vec3f{
      10.02306557647065, 71.41361770095349, -54.07218655560067};
  static const auto c6 = vec3f{
      -3.658713842777788, -22.93153465461149, 18.19190778539828};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

inline vec3f colormap_magma(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  static const auto c0 = vec3f{
      -0.002136485053939582, -0.000749655052795221, -0.005386127855323933};
  static const auto c1 = vec3f{
      0.2516605407371642, 0.6775232436837668, 2.494026599312351};
  static const auto c2 = vec3f{
      8.353717279216625, -3.577719514958484, 0.3144679030132573};
  static const auto c3 = vec3f{
      -27.66873308576866, 14.26473078096533, -13.64921318813922};
  static const auto c4 = vec3f{
      52.17613981234068, -27.94360607168351, 12.94416944238394};
  static const auto c5 = vec3f{
      -50.76852536473588, 29.04658282127291, 4.23415299384598};
  static const auto c6 = vec3f{
      18.65570506591883, -11.48977351997711, -5.601961508734096};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

inline vec3f colormap_inferno(float t) {
  // https://www.shadertoy.com/view/WlfXRN
  static const auto c0 = vec3f{
      0.0002189403691192265, 0.001651004631001012, -0.01948089843709184};
  static const auto c1 = vec3f{
      0.1065134194856116, 0.5639564367884091, 3.932712388889277};
  static const auto c2 = vec3f{
      11.60249308247187, -3.972853965665698, -15.9423941062914};
  static const auto c3 = vec3f{
      -41.70399613139459, 17.43639888205313, 44.35414519872813};
  static const auto c4 = vec3f{
      77.162935699427, -33.40235894210092, -81.80730925738993};
  static const auto c5 = vec3f{
      -71.31942824499214, 32.62606426397723, 73.20951985803202};
  static const auto c6 = vec3f{
      25.13112622477341, -12.24266895238567, -23.07032500287172};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

// Colormaps from {0,1] to color
inline vec3f colormap(float t, colormap_type type) {
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

inline vec3f colorgrade(
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
inline vec4f colorgrade(
    const vec4f& rgba, bool linear, const colorgrade_params& params) {
  auto graded = colorgrade(xyz(rgba), linear, params);
  return {graded.x, graded.y, graded.z, rgba.w};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR SPACES
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
inline vec3f color_to_xyz(const vec3f& col, color_space from) {
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
inline vec3f xyz_to_color(const vec3f& xyz, color_space to) {
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

inline vec3f convert_color(const vec3f& col, color_space from, color_space to) {
  if (from == to) return col;
  return xyz_to_color(color_to_xyz(col, from), to);
}

}  // namespace yocto

#endif
