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

#include <utility>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// COLOR OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion between flots and bytes
inline vec3b  float_to_byte(const vec3f& a);
inline vec3f  byte_to_float(const vec3b& a);
inline vec4b  float_to_byte(const vec4f& a);
inline vec4f  byte_to_float(const vec4b& a);
inline byte   float_to_byte(float a);
inline float  byte_to_float(byte a);
inline ushort float_to_ushort(float a);
inline float  ushort_to_float(ushort a);

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
inline vec4b rgb_to_rgba(const vec3b& rgb);
inline vec3b rgba_to_rgb(const vec4b& rgba);
inline vec4f red_to_rgba(float red);
inline float rgba_to_red(const vec4f& rgba);
inline vec4b red_to_rgba(byte red);
inline byte  rgba_to_red(const vec4b& rgba);
inline vec4f gray_to_rgba(float gray);
inline float rgba_to_gray(const vec4f& rgba);
inline vec4b gray_to_rgba(byte gray);
inline byte  rgba_to_gray(const vec4b& rgba);

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

// Conversion between flots and bytes
inline vec3b float_to_byte(const vec3f& a) {
  return {(byte)clamp(int(a.x * 256), 0, 255),
      (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255)};
}
inline vec3f byte_to_float(const vec3b& a) {
  return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f};
}
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
inline vec4s float_to_ushort(const vec4f& a) {
  return {(ushort)clamp(int(a.x * 65536), 0, 65535),
      (ushort)clamp(int(a.y * 65536), 0, 65535),
      (ushort)clamp(int(a.z * 65536), 0, 65535),
      (ushort)clamp(int(a.w * 65536), 0, 65535)};
}
inline vec4f ushort_to_float(const vec4s& a) {
  return {a.x / 65535.0f, a.y / 65535.0f, a.z / 65535.0f, a.w / 65535.0f};
}
inline ushort float_to_ushort(float a) {
  return (ushort)clamp(int(a * 65536), 0, 65535);
}
inline float ushort_to_float(ushort a) { return a / 65535.0f; }

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
inline vec4b rgb_to_rgba(const vec3b& rgb) {
  return {rgb.x, rgb.y, rgb.z, 255};
}
inline vec3b rgba_to_rgb(const vec4b& rgba) { return xyz(rgba); }
inline vec4f red_to_rgba(float red) { return {red, red, red, 1}; }
inline float rgba_to_red(const vec4f& rgba) { return rgba.x; }
inline vec4b red_to_rgba(byte red) { return {red, red, red, 255}; }
inline byte  rgba_to_red(const vec4b& rgba) {
  return (byte)(((int)rgba.x + (int)rgba.y + (int)rgba.z) / 3);
}
inline vec4f gray_to_rgba(float gray) { return {gray, gray, gray, 1}; }
inline float rgba_to_gray(const vec4f& rgba) {
  return (rgba.x + rgba.y + rgba.z) / 3;
}
inline vec4b gray_to_rgba(byte gray) { return {gray, gray, gray, 255}; }
inline byte  rgba_to_gray(const vec4b& rgba) {
  return (byte)(((int)rgba.x + (int)rgba.y + (int)rgba.z) / 3);
}

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
}

}  // namespace yocto

#endif
