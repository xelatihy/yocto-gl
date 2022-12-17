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

#ifndef _YOCTO_COLOR_H_
#define _YOCTO_COLOR_H_

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
template <typename T>
constexpr kernel vec<T, 3> rgb(const vec<T, 4>& color) {
  return (vec<T, 3>&)color;
}
template <typename T>
constexpr kernel T alpha(const vec<T, 4>& color) {
  return color[3];
}

// Conversion between floats and bytes
template <typename T>
constexpr kernel byte float_to_byte(T a) {
  return (byte)clamp(int(a * 256), 0, 255);
}
template <typename T = float>
constexpr kernel T byte_to_float(byte a) {
  return a / (T)255;
}
template <size_t N, typename T>
constexpr kernel vec<byte, N> float_to_byte(const vec<T, N>& a) {
  return (vec<byte, N>)clamp(vec<int, N>(a * 256), 0, 255);
}
template <size_t N, typename T = float>
constexpr kernel vec<T, N> byte_to_float(const vec<byte, N>& a) {
  return a / (T)255;
}

// Luminance
template <typename T>
constexpr kernel T luminance(const vec<T, 3>& a) {
  auto [r, g, b] = a;
  return (T)0.2126 * r + (T)0.7152 * g + (T)0.0722 * b;
}

// sRGB non-linear curve
template <typename T>
constexpr kernel T srgb_to_rgb(T srgb) {
  return (srgb <= (T)0.04045) ? srgb / (T)12.92
                              : pow((srgb + (T)0.055) / (1 + (T)0.055), (T)2.4);
}
template <typename T>
constexpr kernel T rgb_to_srgb(T rgb) {
  return (rgb <= (T)0.0031308)
             ? (T)12.92 * rgb
             : (1 + (T)0.055) * pow(rgb, 1 / (T)2.4) - (T)0.055;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> srgb_to_rgb(const vec<T, N>& srgb) {
  if constexpr (N == 4) {
    return {srgb_to_rgb(xyz(srgb)), alpha(srgb)};
  } else {
    return map(srgb, [](T srgb) { return srgb_to_rgb(srgb); });
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> rgb_to_srgb(const vec<T, N>& rgb) {
  if constexpr (N == 4) {
    return {rgb_to_srgb(xyz(rgb)), alpha(rgb)};
  } else {
    return map(rgb, [](T rgb) { return rgb_to_srgb(rgb); });
  }
}

// sRGB non-linear curve
template <typename T>
constexpr kernel T srgbb_to_rgb(byte srgb) {
  return srgb_to_rgb(byte_to_float<T>(srgb));
}
template <typename T>
constexpr kernel byte rgb_to_srgbb(T rgb) {
  return float_to_byte(srgb_to_rgb(rgb));
}
template <typename T, size_t N>
constexpr kernel vec<T, N> srgbb_to_rgb(const vec<byte, N>& srgb) {
  return srgb_to_rgb(byte_to_float<N, T>(srgb));
}
template <typename T, size_t N>
constexpr kernel vec<byte, N> rgb_to_srgbb(const vec<T, N>& rgb) {
  return float_to_byte(rgb_to_srgb(rgb));
}

// Conversion between number of channels.
template <typename T>
constexpr kernel vec<T, 4> rgb_to_rgba(const vec<T, 3>& rgb) {
  if constexpr (!std::is_same_v<T, byte>) {
    return {rgb, (T)1};
  } else {
    return {rgb, (byte)255};
  }
}
template <typename T>
constexpr kernel vec<T, 3> rgba_to_rgb(const vec<T, 4>& rgba) {
  return xyz(rgba);
}

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
constexpr kernel vec<T, 3> lincontrast(
    const vec<T, 3>& rgb, T contrast, T grey) {
  return max(grey + (rgb - grey) * (contrast * 2), 0);
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
constexpr kernel vec<T, 3> logcontrast(
    const vec<T, 3>& rgb, T logcontrast, T grey) {
  auto epsilon  = (T)0.0001;
  auto log_grey = log2(grey);
  auto log_ldr  = log2(rgb + epsilon);
  auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
  return max(exp2(adjusted) - epsilon, 0);
}
// Apply an s-shaped contrast.
template <typename T>
constexpr kernel vec<T, 3> contrast(const vec<T, 3>& rgb, T contrast) {
  return gain(rgb, 1 - contrast);
}
// Apply saturation.
template <typename T>
constexpr kernel vec<T, 3> saturate(const vec<T, 3>& rgb, T saturation,
    const vec<T, 3>& weights = vec<T, 3>{
        (T)(1.0 / 3.0), (T)(1.0 / 3.0), (T)(1.0 / 3.0)}) {
  auto grey = dot(weights, rgb);
  return max(grey + (rgb - grey) * (saturation * 2), 0);
}

#ifndef __CUDACC__

// Filmic tonemapping
template <typename T>
constexpr kernel vec<T, 3> tonemap_filmic(
    const vec<T, 3>& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * (T)0.6;  // brings it back to ACES range
    auto ldr = (hdr * hdr * (T)2.51 + hdr * (T)0.03) /
               (hdr * hdr * (T)2.43 + hdr * (T)0.59 + (T)0.14);
    return max(vec<T, 3>{0, 0, 0}, ldr);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    constexpr auto ACESInputMat = transpose(mat<T, 3, 3>{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    constexpr auto ACESOutputMat = transpose(mat<T, 3, 3>{
        {1.60475, -0.53108, -0.07367},
        {-0.10208, 1.10813, -0.00605},
        {-0.00327, -0.07276, 1.07602},
    });
    // RRT => ODT
    auto RRTAndODTFit = [](const vec<T, 3>& v) -> vec<T, 3> {
      return (v * v + v * (T)0.0245786 - (T)0.000090537) /
             (v * v * (T)0.983729 + v * (T)0.4329510 + (T)0.238081);
    };

    auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
    return max(vec<T, 3>{0, 0, 0}, ldr);
  }
}

#else

// Filmic tonemapping
template <typename T>
inline kernel vec<T, 3> tonemap_filmic(
    const vec<T, 3>& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * (T)0.6;  // brings it back to ACES range
    auto ldr = (hdr * hdr * (T)2.51 + hdr * (T)0.03) /
               (hdr * hdr * (T)2.43 + hdr * (T)0.59 + (T)0.14);
    return max(vec<T, 3>{0, 0, 0}, ldr);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    auto ACESInputMat = transpose(mat<T, 3, 3>{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    auto ACESOutputMat = transposemat<T, 3, 3>{
        {1.60475, -0.53108, -0.07367},
        {-0.10208, 1.10813, -0.00605},
        {-0.00327, -0.07276, 1.07602},
    });
    // RRT => ODT
    auto RRTAndODTFit = [](const vec<T, 3>& v) -> vec<T, 3> {
      return (v * v + v * (T)0.0245786 - (T)0.000090537) /
             (v * v * (T)0.983729(T) + v * (T)0.4329510 + (T)0.238081);
    };

    auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
    return max(vec<T, 3>{0, 0, 0}, ldr);
  }
}

#endif

template <typename T>
constexpr kernel vec<T, 3> tonemap(
    const vec<T, 3>& hdr, T exposure, bool filmic, bool srgb = true) {
  auto rgb = hdr;
  if (exposure != 0) rgb *= exp2(exposure);
  if (filmic) rgb = tonemap_filmic(rgb);
  if (srgb) rgb = rgb_to_srgb(rgb);
  return rgb;
}
template <typename T>
constexpr kernel vec<T, 4> tonemap(
    const vec<T, 4>& hdr, T exposure, bool filmic, bool srgb = true) {
  auto ldr = tonemap(xyz(hdr), exposure, filmic, srgb);
  return {ldr, alpha(hdr)};
}

// Composite colors
template <typename T>
constexpr kernel vec<T, 4> composite(const vec<T, 4>& a, const vec<T, 4>& b) {
  if (a.w == 0 && b.w == 0) return {0, 0, 0, 0};
  auto cc = xyz(a) * alpha(a) + xyz(b) * alpha(b) * (1 - alpha(a));
  auto ca = alpha(a) + alpha(b) * (1 - alpha(a));
  return {cc / ca, ca};
}

// Convert between CIE XYZ and RGB
template <typename T>
constexpr kernel vec<T, 3> rgb_to_xyz(const vec<T, 3>& rgb) {
  // https://en.wikipedia.org/wiki/SRGB
  constexpr auto m = mat<T, 3, 3>{
      {(T)0.4124, (T)0.2126, (T)0.0193},
      {(T)0.3576, (T)0.7152, (T)0.1192},
      {(T)0.1805, (T)0.0722, (T)0.9504},
  };
  return m * rgb;
}
template <typename T>
constexpr kernel vec<T, 3> xyz_to_rgb(const vec<T, 3>& xyz) {
  // https://en.wikipedia.org/wiki/SRGB
  constexpr auto m = mat<T, 3, 3>{
      {(T) + 3.2406, (T)-0.9689, (T) + 0.0557},
      {(T)-1.5372, (T) + 1.8758, (T)-0.2040},
      {(T)-0.4986, (T) + 0.0415, (T) + 1.0570},
  };
  return m * xyz;
}

// Convert between CIE XYZ and xyY
template <typename T>
constexpr kernel vec<T, 3> xyz_to_xyY(const vec<T, 3>& xyz) {
  if (xyz == 0) return {0, 0, 0};
  auto [x, y, z] = xyz;
  return {x / (x + y + z), y / (x + y + z), y};
}
template <typename T>
constexpr kernel vec<T, 3> xyY_to_xyz(const vec<T, 3>& xyY) {
  auto [x, y, Y] = xyY;
  if (y == 0) return {0, 0, 0};
  return {x * Y / y, Y, (1 - x - y) * Y / y};
}

// Convert HSV to RGB
template <typename T>
constexpr kernel vec<T, 3> hsv_to_rgb(const vec<T, 3>& hsv) {
  // from Imgui.cpp
  auto [h, s, v] = hsv;
  if (s == 0) return {v, v, v};

  h      = fmod(h, (T)1) / ((T)60 / (T)360);
  auto i = (int)h;
  auto f = h - (T)i;
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
    default: return {v, p, q};
  }
}

template <typename T>
constexpr kernel vec<T, 3> rgb_to_hsv(const vec<T, 3>& rgb) {
  // from Imgui.cpp
  auto [r, g, b] = rgb;
  auto K         = (T)0;
  if (g < b) {
    std::swap(g, b);
    K = -1;
  }
  if (r < g) {
    std::swap(r, g);
    K = -2 / (T)6 - K;
  }

  auto c = r - (g < b ? g : b);  // chroma
  return {abs(K + (g - b) / (6 * c + (T)1e-20)), c / (r + (T)1e-20), r};
}

// Approximate color of blackbody radiation from wavelength in nm.
template <typename T>
constexpr kernel vec<T, 3> blackbody_to_rgb(T temperature) {
  // clamp to valid range
  auto t = clamp(temperature, (T)1667.0, (T)25000.0) / (T)1000.0;
  // compute x
  auto x = (T)0;
  if (temperature < (T)4000.0) {
    x = (T)-0.2661239 * 1 / (t * t * t) - (T)0.2343589 * 1 / (t * t) +
        (T)0.8776956 * (1 / t) + (T)0.179910;
  } else {
    x = (T)-3.0258469 * 1 / (t * t * t) + (T)2.1070379 * 1 / (t * t) +
        (T)0.2226347 * (1 / t) + (T)0.240390;
  }
  // compute y
  auto y = (T)0;
  if (temperature < (T)2222) {
    y = (T)-1.1063814 * (x * x * x) - (T)1.34811020 * (x * x) +
        (T)2.18555832 * x - (T)0.20219683;
  } else if (temperature < (T)4000.0) {
    y = (T)-0.9549476 * (x * x * x) - (T)1.37418593 * (x * x) +
        (T)2.09137015 * x - (T)0.16748867;
  } else {
    y = (T)3.0817580 * (x * x * x) - (T)5.87338670 * (x * x) +
        (T)3.75112997 * x - (T)0.37001483;
  }
  return xyz_to_rgb(xyY_to_xyz(vec<T, 3>{x, y, 1}));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COLOR MAPS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
constexpr kernel vec<T, 3> colormap_viridis(T t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec<T, 3>{
      (T)0.2777273272234177, (T)0.005407344544966578, (T)0.3340998053353061};
  constexpr auto c1 = vec<T, 3>{
      (T)0.1050930431085774, (T)1.404613529898575, (T)1.384590162594685};
  constexpr auto c2 = vec<T, 3>{
      (T)-0.3308618287255563, (T)0.214847559468213, (T)0.09509516302823659};
  constexpr auto c3 = vec<T, 3>{
      (T)-4.634230498983486, (T)-5.799100973351585, (T)-19.33244095627987};
  constexpr auto c4 = vec<T, 3>{
      (T)6.228269936347081, (T)14.17993336680509, (T)56.69055260068105};
  constexpr auto c5 = vec<T, 3>{
      (T)4.776384997670288, (T)-13.74514537774601, (T)-65.35303263337234};
  constexpr auto c6 = vec<T, 3>{
      (T)-5.435455855934631, (T)4.645852612178535, (T)26.3124352495832};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

template <typename T>
constexpr kernel vec<T, 3> colormap_plasma(T t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec<T, 3>{
      (T)0.05873234392399702, (T)0.02333670892565664, (T)0.5433401826748754};
  constexpr auto c1 = vec<T, 3>{
      (T)2.176514634195958, (T)0.2383834171260182, (T)0.7539604599784036};
  constexpr auto c2 = vec<T, 3>{
      (T)-2.689460476458034, (T)-7.455851135738909, (T)3.110799939717086};
  constexpr auto c3 = vec<T, 3>{
      (T)6.130348345893603, (T)42.3461881477227, (T)-28.51885465332158};
  constexpr auto c4 = vec<T, 3>{
      (T)-11.10743619062271, (T)-82.66631109428045, (T)60.13984767418263};
  constexpr auto c5 = vec<T, 3>{
      (T)10.02306557647065, (T)71.41361770095349, (T)-54.07218655560067};
  constexpr auto c6 = vec<T, 3>{
      (T)-3.658713842777788, (T)-22.93153465461149, (T)18.19190778539828};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

template <typename T>
constexpr kernel vec<T, 3> colormap_magma(T t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec<T, 3>{(T)-0.002136485053939582,
      (T)-0.000749655052795221, (T)-0.005386127855323933};
  constexpr auto c1 = vec<T, 3>{
      (T)0.2516605407371642, (T)0.6775232436837668, (T)2.494026599312351};
  constexpr auto c2 = vec<T, 3>{
      (T)8.353717279216625, (T)-3.577719514958484, (T)0.3144679030132573};
  constexpr auto c3 = vec<T, 3>{
      (T)-27.66873308576866, (T)14.26473078096533, (T)-13.64921318813922};
  constexpr auto c4 = vec<T, 3>{
      (T)52.17613981234068, (T)-27.94360607168351, (T)12.94416944238394};
  constexpr auto c5 = vec<T, 3>{
      (T)-50.76852536473588, (T)29.04658282127291, (T)4.23415299384598};
  constexpr auto c6 = vec<T, 3>{
      (T)18.65570506591883, (T)-11.48977351997711, (T)-5.601961508734096};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

template <typename T>
constexpr kernel vec<T, 3> colormap_inferno(T t) {
  // https://www.shadertoy.com/view/WlfXRN
  constexpr auto c0 = vec<T, 3>{(T)0.0002189403691192265,
      (T)0.001651004631001012, (T)-0.01948089843709184};
  constexpr auto c1 = vec<T, 3>{
      (T)0.1065134194856116, 0.5639564367884091, (T)3.932712388889277};
  constexpr auto c2 = vec<T, 3>{
      (T)11.60249308247187, -3.972853965665698, (T)-15.9423941062914};
  constexpr auto c3 = vec<T, 3>{
      (T)-41.70399613139459, 17.43639888205313, (T)44.35414519872813};
  constexpr auto c4 = vec<T, 3>{
      (T)77.162935699427, -33.40235894210092, (T)-81.80730925738993};
  constexpr auto c5 = vec<T, 3>{
      (T)-71.31942824499214, (T)(T)32.62606426397723, (T)73.20951985803202};
  constexpr auto c6 = vec<T, 3>{
      (T)25.13112622477341, (T)-12.24266895238567, (T)-23.07032500287172};
  return c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))));
}

// Colormap type
enum struct colormap_type { viridis, plasma, magma, inferno };

// Colormaps from {0,1] to color
template <typename T>
constexpr kernel vec<T, 3> colormap(T t, colormap_type type) {
  t = clamp(t, (T)0, (T)1);
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
template <typename T = float>
struct colorgrade_gparams {
  T         exposure         = 0;
  vec<T, 3> tint             = {1, 1, 1};
  T         lincontrast      = 0.5;
  T         logcontrast      = 0.5;
  T         linsaturation    = 0.5;
  bool      filmic           = false;
  bool      srgb             = true;
  T         contrast         = 0.5;
  T         saturation       = 0.5;
  T         shadows          = 0.5;
  T         midtones         = 0.5;
  T         highlights       = 0.5;
  vec<T, 3> shadows_color    = {1, 1, 1};
  vec<T, 3> midtones_color   = {1, 1, 1};
  vec<T, 3> highlights_color = {1, 1, 1};
};

// Alias
using colorgrade_params = colorgrade_gparams<float>;

template <typename T>
constexpr kernel vec<T, 3> colorgrade(
    const vec<T, 3>& rgb_, bool linear, const colorgrade_gparams<T>& params) {
  auto rgb = rgb_;
  if (params.exposure != 0) rgb *= exp2(params.exposure);
  if (params.tint != 1) rgb *= params.tint;
  if (params.lincontrast != (T)0.5)
    rgb = lincontrast(rgb, params.lincontrast, linear ? (T)0.18 : (T)0.5);
  if (params.logcontrast != (T)0.5)
    rgb = logcontrast(rgb, params.logcontrast, linear ? (T)0.18 : (T)0.5);
  if (params.linsaturation != (T)0.5) rgb = saturate(rgb, params.linsaturation);
  if (params.filmic) rgb = tonemap_filmic(rgb);
  if (linear && params.srgb) rgb = rgb_to_srgb(rgb);
  if (params.contrast != (T)0.5) rgb = contrast(rgb, params.contrast);
  if (params.saturation != (T)0.5) rgb = saturate(rgb, params.saturation);
  if (params.shadows != (T)0.5 || params.midtones != (T)0.5 ||
      params.highlights != (T)0.5 || params.shadows_color != 1 ||
      params.midtones_color != 1 || params.highlights_color != 1) {
    auto lift  = params.shadows_color;
    auto gamma = params.midtones_color;
    auto gain  = params.highlights_color;

    lift      = lift - mean(lift) + params.shadows - (T)0.5;
    gain      = gain - mean(gain) + params.highlights + (T)0.5;
    auto grey = gamma - mean(gamma) + params.midtones;
    gamma     = log(((T)0.5 - lift) / (gain - lift)) / log(grey);

    // apply_image
    auto lerp_value = clamp(pow(rgb, 1 / gamma), (T)0, (T)1);
    rgb             = gain * lerp_value + lift * (1 - lerp_value);
  }
  return rgb;
}

template <typename T>
constexpr kernel vec<T, 4> colorgrade(
    const vec<T, 4>& rgba, bool linear, const colorgrade_params& params) {
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
template <typename T>
constexpr vec<T, 3> color_to_xyz(const vec<T, 3>& col, color_space from);
template <typename T>
constexpr vec<T, 3> xyz_to_color(const vec<T, 3>& xyz, color_space to);

// Conversion between rgb color spaces
template <typename T>
constexpr vec<T, 3> convert_color(
    const vec<T, 3>& col, color_space from, color_space to);

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
template <typename T>
struct color_space_gparams {
  // Curve type
  enum struct curve_t {
    // clang-format off
    linear, gamma, linear_gamma, aces_cc, aces_cct, pq, hlg
    // clang-format on
  };
  // primaries
  vec<T, 2>    red_chromaticity;    // xy chromaticity of the red primary
  vec<T, 2>    green_chromaticity;  // xy chromaticity of the green primary
  vec<T, 2>    blue_chromaticity;   // xy chromaticity of the blue primary
  vec<T, 2>    white_chromaticity;  // xy chromaticity of the white point
  mat<T, 3, 3> rgb_to_xyz_mat;      // matrix from rgb to xyz
  mat<T, 3, 3> xyz_to_rgb_mat;      // matrix from xyz to rgb
  // tone curve
  curve_t   curve_type;
  T         curve_gamma;  // gamma for power curves
  vec<T, 4> curve_abcd;   // tone curve values for linear_gamma curves
};

// Alias
using color_space_params = color_space_gparams<float>;

// Compute the rgb -> xyz matrix from the color space definition
// Input: red, green, blue, white (x,y) chromaticities
// Algorithm from: SMPTE Recommended Practice RP 177-1993
// http://car.france3.mars.free.fr/HD/INA-%2026%20jan%2006/SMPTE%20normes%20et%20confs/rp177.pdf
template <typename T>
constexpr mat<T, 3, 3> rgb_to_xyz_mat(const vec<T, 2>& rc, const vec<T, 2>& gc,
    const vec<T, 2>& bc, const vec<T, 2>& wc) {
  auto r = vec<T, 3>{rc, 1 - sum(rc)}, g = vec<T, 3>{gc, 1 - sum(gc)},
       b = vec<T, 3>{bc, 1 - sum(bc)}, w = vec<T, 3>{wc, 1 - sum(wc)};
  auto [_, wy, __]  = w;
  auto [cr, cg, cb] = inverse({r, g, b}) * w / wy;
  return mat<T, 3, 3>{cr * r, cg * g, cb * b};
}

// Construct an RGB color space. Predefined color spaces below
template <typename T>
constexpr color_space_params get_color_scape_params(color_space space) {
  constexpr auto make_linear_rgb_space =
      [](const vec<T, 2>& red, const vec<T, 2>& green, const vec<T, 2>& blue,
          const vec<T, 2>& white) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)),
            color_space_params::curve_t::linear};
      };
  constexpr auto make_gamma_rgb_space =
      [](const vec<T, 2>& red, const vec<T, 2>& green, const vec<T, 2>& blue,
          const vec<T, 2>& white, T gamma,
          const vec4f& curve_abcd = vec4f{0, 0, 0, 0}) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)),
            curve_abcd == vec4f{0, 0, 0, 0}
                ? color_space_params::curve_t::gamma
                : color_space_params::curve_t::linear_gamma};
      };
  constexpr auto make_other_rgb_space =
      [](const vec<T, 2>& red, const vec<T, 2>& green, const vec<T, 2>& blue,
          const vec<T, 2>& white, color_space_params::curve_t curve_type) {
        return color_space_params{red, green, blue, white,
            rgb_to_xyz_mat(red, green, blue, white),
            inverse(rgb_to_xyz_mat(red, green, blue, white)), curve_type};
      };

  // color space parameters
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto rgb_params = make_linear_rgb_space({(T)0.6400, (T)0.3300},
      {(T)0.3000, (T)0.6000}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290});
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto srgb_params = make_gamma_rgb_space({(T)0.6400, (T)0.3300},
      {(T)0.3000, (T)0.6000}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290},
      (T)2.4, {(T)1.055, (T)0.055, (T)12.92, (T)0.0031308});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_System
  constexpr auto aces2065_params = make_linear_rgb_space({(T)0.7347, (T)0.2653},
      {(T)0.0000, (T)1.0000}, {(T)0.0001, (T)-0.0770},
      {(T)0.32168, (T)0.33767});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescg_params = make_linear_rgb_space({(T)0.7130, (T)0.2930},
      {(T)0.1650, (T)0.8300}, {(T)0.1280, (T)0.0440}, {(T)0.32168, (T)0.33767});
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescc_params = make_other_rgb_space({(T)0.7130, (T)0.2930},
      {(T)0.1650, (T)0.8300}, {(T)0.1280, (T)0.0440}, {(T)0.32168, (T)0.33767},
      color_space_params::curve_t::aces_cc);
  // https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
  constexpr auto acescct_params = make_other_rgb_space({(T)0.7130, (T)0.2930},
      {(T)0.1650, (T)0.8300}, {(T)0.1280, (T)0.0440}, {(T)0.32168, (T)0.33767},
      color_space_params::curve_t::aces_cct);
  // https://en.wikipedia.org/wiki/Adobe_RGB_color_space
  constexpr auto adobe_params = make_gamma_rgb_space({(T)0.6400, (T)0.3300},
      {(T)0.2100, (T)0.7100}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290},
      (T)2.19921875);
  // https://en.wikipedia.org/wiki/Rec._709
  constexpr auto rec709_params = make_gamma_rgb_space({(T)0.6400, (T)0.3300},
      {(T)0.3000, (T)0.6000}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290},
      1 / (T)0.45, {(T)1.099, (T)0.099, (T)4.500, (T)0.018});
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2020_params = make_gamma_rgb_space({(T)0.7080, (T)0.2920},
      {(T)0.1700, (T)0.7970}, {(T)0.1310, (T)0.0460}, {(T)0.3127, (T)0.3290},
      1 / (T)0.45,
      {(T)1.09929682680944, (T)0.09929682680944, (T)4.5, (T)0.018053968510807});
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2100pq_params = make_other_rgb_space({(T)0.7080, (T)0.2920},
      {(T)0.1700, (T)0.7970}, {(T)0.1310, (T)0.0460}, {(T)0.3127, (T)0.3290},
      color_space_params::curve_t::pq);
  // https://en.wikipedia.org/wiki/Rec._2020
  constexpr auto rec2100hlg_params = make_other_rgb_space(
      {(T)0.7080, (T)0.2920}, {(T)0.1700, (T)0.7970}, {(T)0.1310, (T)0.0460},
      {(T)0.3127, (T)0.3290}, color_space_params::curve_t::hlg);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3dci_params = make_gamma_rgb_space({(T)0.6800, (T)0.3200},
      {(T)0.2650, (T)0.6900}, {(T)0.1500, (T)0.0600}, {(T)0.3140, (T)0.3510},
      (T)1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3d60_params = make_gamma_rgb_space({(T)0.6800, (T)0.3200},
      {(T)0.2650, (T)0.6900}, {(T)0.1500, (T)0.0600}, {(T)0.32168, (T)0.33767},
      (T)1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3d65_params = make_gamma_rgb_space({(T)0.6800, (T)0.3200},
      {(T)0.2650, (T)0.6900}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290},
      (T)1.6);
  // https://en.wikipedia.org/wiki/DCI-P3
  constexpr auto p3display_params = make_gamma_rgb_space({(T)0.6800, (T)0.3200},
      {(T)0.2650, (T)0.6900}, {(T)0.1500, (T)0.0600}, {(T)0.3127, (T)0.3290},
      (T)2.4, {(T)1.055, (T)0.055, (T)12.92, (T)0.0031308});
  // https://en.wikipedia.org/wiki/ProPhoto_RGB_color_space
  constexpr auto prophoto_params = make_gamma_rgb_space({(T)0.7347, (T)0.2653},
      {(T)0.1596, (T)0.8404}, {(T)0.0366, (T)0.0001}, {(T)0.3457, (T)0.3585},
      (T)1.8, {(T)1.0, (T)0.0, (T)16.0, (T)0.001953125});

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
template <typename T>
constexpr T gamma_display_to_linear(T x, T gamma) {
  return pow(x, gamma);
}
template <typename T>
constexpr T gamma_linear_to_display(T x, T gamma) {
  return pow(x, 1 / gamma);
}
template <typename T>
constexpr vec<T, 3> gamma_display_to_linear(const vec<T, 3>& rgb, T gamma) {
  return map(rgb, [gamma](T a) { return gamma_display_to_linear(a, gamma); });
}
template <typename T>
constexpr vec<T, 3> gamma_linear_to_display(const vec<T, 3>& rgb, T gamma) {
  return map(rgb, [gamma](T a) { return gamma_linear_to_display(a, gamma); });
}

// https://en.wikipedia.org/wiki/Rec._709
template <typename T>
constexpr T gamma_display_to_linear(T x, T gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  return (x < 1 / d) ? (x / c) : pow((x + b) / a, gamma);
}
template <typename T>
constexpr T gamma_linear_to_display(T x, T gamma, const vec4f& abcd) {
  auto& [a, b, c, d] = abcd;
  return (x < d) ? (x * c) : (a * pow(x, 1 / gamma) - b);
}
template <typename T>
constexpr vec<T, 3> gamma_display_to_linear(
    const vec<T, 3>& rgb, T gamma, const vec4f& abcd) {
  return map(rgb, [&](T a) { return gamma_display_to_linear(a, gamma, abcd); });
}
template <typename T>
constexpr vec<T, 3> gamma_linear_to_display(
    const vec<T, 3>& rgb, T gamma, const vec4f& abcd) {
  return map(rgb, [&](T a) { return gamma_linear_to_display(a, gamma, abcd); });
}

// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
template <typename T>
constexpr T acescc_display_to_linear(T x) {
  if (x < (T)-0.3013698630) {  // (9.72-15)/17.52
    return (exp2(x * (T)17.52 - (T)9.72) - exp2((T)-16.0)) * 2;
  } else if (x < (log2((T)65504.0) + (T)9.72) / (T)17.52) {
    return exp2(x * (T)17.52 - (T)9.72);
  } else {  // (in >= (log2(65504)+9.72)/17.52)
    return (T)65504.0;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
template <typename T>
constexpr T acescct_display_to_linear(T x) {
  if (x < (T)0.155251141552511) {
    return (x - (T)0.0729055341958355) / (T)10.5402377416545;
  } else {
    return exp2(x * (T)17.52 - (T)9.72);
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
template <typename T>
constexpr T acescc_linear_to_display(T x) {
  if (x <= 0) {
    return (T)-0.3584474886;  // =(log2( pow(2.,-16.))+9.72)/17.52
  } else if (x < exp2((T)-15.0)) {
    return (log2(exp2((T)-16.0) + x * (T)0.5) + (T)9.72) / (T)17.52;
  } else {  // (in >= pow(2.,-15))
    return (log2(x) + (T)9.72) / (T)17.52;
  }
}
// https://en.wikipedia.org/wiki/Academy_Color_Encoding_Systemx
template <typename T>
constexpr T acescct_linear_to_display(T x) {
  if (x <= (T)0.0078125) {
    return (T)10.5402377416545 * x + (T)0.0729055341958355;
  } else {
    return (log2(x) + (T)9.72) / (T)17.52;
  }
}
template <typename T>
constexpr vec<T, 3> acescc_display_to_linear(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return acescc_display_to_linear(a); });
}
template <typename T>
constexpr vec<T, 3> acescct_display_to_linear(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return acescct_display_to_linear(a); });
}
template <typename T>
constexpr vec<T, 3> acescc_linear_to_display(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return acescc_linear_to_display(a); });
}
template <typename T>
constexpr vec<T, 3> acescct_linear_to_display(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return acescct_linear_to_display(a); });
}

// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// https://github.com/ampas/aces-dev/blob/master/transforms/ctl/lib/ACESlib.Utilities_Color.ctl
// In PQ, we assume that the linear luminance in [0,1] corresponds to
// [0,10000] cd m^2
template <typename T>
constexpr T pq_display_to_linear(T x) {
  auto Np = pow(x, 1 / (T)78.84375);
  auto L  = max(Np - (T)0.8359375, (T)0.0);
  L       = L / ((T)18.8515625 - (T)18.6875 * Np);
  L       = pow(L, 1 / (T)0.1593017578125);
  return L;
}
template <typename T>
constexpr T pq_linear_to_display(T x) {
  return pow(((T)0.8359375 + (T)18.8515625 * pow(x, (T)0.1593017578125)) /
                 (1 + (T)18.6875 * pow(x, (T)0.1593017578125)),
      (T)78.84375);
}
template <typename T>
constexpr vec<T, 3> pq_display_to_linear(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return pq_display_to_linear(a); });
}
template <typename T>
constexpr vec<T, 3> pq_linear_to_display(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return pq_linear_to_display(a); });
}

// https://en.wikipedia.org/wiki/High-dynamic-range_video#Perceptual_Quantizer
// In HLG, we assume that the linear luminance in [0,1] corresponds to
// [0,1000] cd m^2. Note that the version we report here is scaled in [0,1]
// range for nominal luminance. But HLG was initially defined in the [0,12]
// range where it maps 1 to 0.5 and 12 to 1. For use in HDR tonemapping that is
// likely a better range to use.
template <typename T>
constexpr T hlg_display_to_linear(T x) {
  if (x < (T)0.5) {
    return 3 * 3 * x * x;
  } else {
    return (exp((x - (T)0.55991073) / (T)0.17883277) + (T)0.28466892) / 12;
  }
}
template <typename T>
constexpr T hlg_linear_to_display(T x) {
  if (x < 1 / (T)12.0) {
    return sqrt(3 * x);
  } else {
    return (T)0.17883277 * log(12 * x - (T)0.28466892) + (T)0.55991073;
  }
}
template <typename T>
constexpr vec<T, 3> hlg_display_to_linear(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return hlg_display_to_linear(a); });
}
template <typename T>
constexpr vec<T, 3> hlg_linear_to_display(const vec<T, 3>& rgb) {
  return map(rgb, [](T a) { return hlg_linear_to_display(a); });
}

// Conversion to/from xyz
template <typename T>
constexpr vec<T, 3> color_to_xyz(const vec<T, 3>& col, color_space from) {
  constexpr auto space = get_color_scape_params<T>(from);
  auto           rgb   = col;
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
template <typename T>
constexpr vec<T, 3> xyz_to_color(const vec<T, 3>& xyz, color_space to) {
  constexpr auto space = get_color_scape_params<T>(to);
  auto           rgb   = space.xyz_to_rgb_mat * xyz;
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

template <typename T>
constexpr vec<T, 3> convert_color(
    const vec<T, 3>& col, color_space from, color_space to) {
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
