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

#include <stb_image/stb_image_resize2.h>

#include <memory>
#include <stdexcept>

#include "yocto_color.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt) {
  return transform_image(bt, [](vec4b a) { return byte_to_float(a); });
}
image<vec4b> float_to_byte(const image<vec4f>& fl) {
  return transform_image(fl, [](vec4f a) { return float_to_byte(a); });
}

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_rgb(const image<vec4f>& srgb) {
  return transform_image(srgb, [](vec4f a) { return srgb_to_rgb(a); });
}
image<vec4f> rgb_to_srgb(const image<vec4f>& rgb) {
  return transform_image(rgb, [](vec4f a) { return rgb_to_srgb(a); });
}
image<vec4f> srgbb_to_rgb(const image<vec4b>& srgb) {
  return transform_image(srgb, [](vec4b a) { return srgbb_to_rgb(a); });
}
image<vec4b> rgb_to_srgbb(const image<vec4f>& rgb) {
  return transform_image(rgb, [](vec4f a) { return rgb_to_srgbb(a); });
}

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
  return transform_image(hdr, [exposure, filmic, srgb](vec4f a) {
    return tonemap(a, exposure, filmic, srgb);
  });
}
image<vec4b> tonemapb_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
  return transform_image(hdr, [exposure, filmic, srgb](vec4f a) {
    return float_to_byte(tonemap(a, exposure, filmic, srgb));
  });
}
void tonemap_image(image<vec4f>& ldr, const image<vec4f>& hdr, float exposure,
    bool filmic, bool srgb) {
  return transform_image(ldr, hdr, [exposure, filmic, srgb](vec4f a) {
    return tonemap(a, exposure, filmic, srgb);
  });
}

// Apply exposure and filmic tone mapping
image<vec4f> colorgrade_image(
    const image<vec4f>& img, bool linear, const colorgrade_params& params) {
  return transform_image(
      img, [linear, params](vec4f a) { return colorgrade(a, linear, params); });
}
void colorgrade_image(image<vec4f>& graded, const image<vec4f>& img,
    bool linear, const colorgrade_params& params) {
  return transform_image(graded, img,
      [linear, params](vec4f a) { return colorgrade(a, linear, params); });
}

// compute white balance
vec3f compute_white_balance(const image<vec4f>& img) {
  auto rgb = vec3f{0, 0, 0};
  for (auto& p : img) rgb += xyz(p);
  if (rgb == vec3f{0, 0, 0}) return {0, 0, 0};
  return rgb / max(rgb);
}

// image compositing
image<vec4f> composite_image(
    const image<vec4f>& foreground, const image<vec4f>& background) {
  return transform_images(
      foreground, background, [](vec4f a, vec4f b) { return composite(a, b); });
}

// removes alpha
image<vec4f> remove_alpha(const image<vec4f>& img) {
  return transform_image(img, [](vec4f a) -> vec4f { return {xyz(a), 1}; });
}

// turns alpha into a gray scale image
image<vec4f> alpha_to_gray(const image<vec4f>& img) {
  return transform_image(img, [](vec4f a) -> vec4f {
    auto g = a.w;
    return {g, g, g, 1};
  });
}

image<vec4f> image_difference(
    const image<vec4f>& a, const image<vec4f>& b, bool display) {
  return transform_images(a, b, [display](vec4f a, vec4f b) -> vec4f {
    auto diff = abs(a - b);
    if (display) {
      auto d = max(diff);
      return {d, d, d, 1};
    } else {
      return diff;
    }
  });
}

image<vec4f> resize_image(const image<vec4f>& img, vec2i resize) {
  auto size = img.size();
  if (resize.x == 0 && resize.y == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (resize.y == 0) {
    resize.y = (int)round(resize.x * (double)size.y / (double)size.x);
  } else if (resize.x == 0) {
    resize.x = (int)round(resize.y * (double)size.x / (double)size.y);
  }
  auto res = image<vec4f>{resize};
  stbir_resize_float_linear((float*)img.data(), size.x, size.y,
      sizeof(vec4f) * size.x, (float*)res.data(), resize.x, resize.y,
      sizeof(vec4f) * resize.x, STBIR_RGBA);
  return res;
}
image<vec4b> resize_image(const image<vec4b>& img, vec2i resize) {
  auto size = img.size();
  if (resize.x == 0 && resize.y == 0) {
    throw std::invalid_argument{"bad image size in resize"};
  }
  if (resize.y == 0) {
    resize.y = (int)round(resize.x * (double)size.y / (double)size.x);
  } else if (resize.x == 0) {
    resize.x = (int)round(resize.y * (double)size.x / (double)size.y);
  }
  auto res = image<vec4b>{resize};
  stbir_resize_uint8_linear((byte*)img.data(), size.x, size.y,
      sizeof(vec4b) * size.x, (byte*)res.data(), resize.x, resize.y,
      sizeof(vec4b) * resize.x, STBIR_RGBA);
  return res;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Comvert a bump map to a normal map.
void bump_to_normal(
    image<vec4f>& normalmap, const image<vec4f>& bumpmap, float scale) {
  if (normalmap.size() != bumpmap.size()) normalmap = {bumpmap.size()};
  auto dx = 1.0f / bumpmap.size().x, dy = 1.0f / bumpmap.size().y;
  auto size = bumpmap.size();
  for (auto [i, j] : range(size)) {
    auto i1 = (i + 1) % size.x, j1 = (j + 1) % size.y;
    auto p00 = bumpmap[{i, j}], p10 = bumpmap[{i1, j}], p01 = bumpmap[{i, j1}];
    auto g00    = (p00.x + p00.y + p00.z) / 3;
    auto g01    = (p01.x + p01.y + p01.z) / 3;
    auto g10    = (p10.x + p10.y + p10.z) / 3;
    auto normal = vec3f{
        scale * (g00 - g10) / dx, scale * (g00 - g01) / dy, 1.0f};
    normal.y = -normal.y;  // make green pointing up, even if y axis
                           // points down
    normal            = normalize(normal) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
    normalmap[{i, j}] = {normal.x, normal.y, normal.z, 1};
  }
}
image<vec4f> bump_to_normal(const image<vec4f>& bumpmap, float scale) {
  auto normalmap = image<vec4f>{bumpmap.size()};
  bump_to_normal(normalmap, bumpmap, scale);
  return normalmap;
}

template <typename Shader>
static image<vec4f> make_proc_image(vec2i size, Shader&& shader) {
  auto image_ = image<vec4f>{size};
  auto scale  = 1.0f / max(size);
  for (auto ij : range(size)) {
    image_[ij] = shader((vec2f)ij * scale);
  }
  return image_;
}

// Make an image
image<vec4f> make_grid(vec2i size, float scale, vec4f color0, vec4f color1) {
  return make_proc_image(size, [=](vec2f uv) {
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

image<vec4f> make_checker(vec2i size, float scale, vec4f color0, vec4f color1) {
  return make_proc_image(size, [=](vec2f uv) {
    uv *= 4 * scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    auto c = uv.x <= 0.5f != uv.y <= 0.5f;
    return c ? color0 : color1;
  });
}

image<vec4f> make_bumps(vec2i size, float scale, vec4f color0, vec4f color1) {
  return make_proc_image(size, [=](vec2f uv) {
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

image<vec4f> make_ramp(vec2i size, float scale, vec4f color0, vec4f color1) {
  return make_proc_image(size, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return lerp(color0, color1, uv.x);
  });
}

image<vec4f> make_gammaramp(
    vec2i size, float scale, vec4f color0, vec4f color1) {
  return make_proc_image(size, [=](vec2f uv) {
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

image<vec4f> make_uvramp(vec2i size, float scale) {
  return make_proc_image(size, [=](vec2f uv) {
    uv *= scale;
    uv -= vec2f{(float)(int)uv.x, (float)(int)uv.y};
    return vec4f{uv.x, uv.y, 0, 1};
  });
}

image<vec4f> make_uvgrid(vec2i size, float scale, bool colored) {
  return make_proc_image(size, [=](vec2f uv) {
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

image<vec4f> make_colormapramp(vec2i size, float scale) {
  return make_proc_image(size, [=](vec2f uv) {
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

// Add image border
image<vec4f> add_border(const image<vec4f>& image, float width, vec4f color) {
  auto result = image;
  auto scale  = 1.0f / max(image.size());
  for (auto ij : range(image.size())) {
    auto uv = (vec2f)ij * scale;
    if (uv.x < width || uv.y < width || uv.x > image.size().x * scale - width ||
        uv.y > image.size().y * scale - width) {
      result[ij] = color;
    }
  }
  return result;
}

// Implementation of sunsky modified heavily from pbrt
image<vec4f> make_sunsky(vec2i size, float theta_sun, float turbidity,
    bool has_sun, float sun_intensity, float sun_radius, vec3f ground_albedo) {
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
    return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000
                                                   : vec3f{0, 0, 0};
  };

  // Make the sun sky image
  auto img = image<vec4f>{size};
  for (auto j = 0; j < size.y / 2; j++) {
    auto theta = pif * ((j + 0.5f) / size.y);
    theta      = clamp(theta, 0.0f, pif / 2 - flt_eps);
    for (int i = 0; i < size.x; i++) {
      auto phi = 2 * pif * (float(i + 0.5f) / size.x);
      auto w = vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
      auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
      auto sky_col = sky(theta, gamma, theta_sun);
      auto sun_col = sun(theta, gamma);
      auto col     = sky_col + sun_col;
      img[{i, j}]  = {col.x, col.y, col.z, 1};
    }
  }

  if (ground_albedo != vec3f{0, 0, 0}) {
    auto ground = vec3f{0, 0, 0};
    for (auto j = 0; j < size.y / 2; j++) {
      auto theta = pif * ((j + 0.5f) / size.y);
      for (int i = 0; i < size.x; i++) {
        auto pxl   = img[{i, j}];
        auto le    = vec3f{pxl.x, pxl.y, pxl.z};
        auto angle = sin(theta) * 4 * pif / (size.x * size.y);
        ground += le * (ground_albedo / pif) * cos(theta) * angle;
      }
    }
    for (auto j = size.y / 2; j < size.y; j++) {
      for (int i = 0; i < size.x; i++) {
        img[{i, j}] = {ground.x, ground.y, ground.z, 1};
      }
    }
  } else {
    for (auto j = size.y / 2; j < size.y; j++) {
      for (int i = 0; i < size.x; i++) {
        img[{i, j}] = {0, 0, 0, 1};
      }
    }
  }

  // done
  return img;
}

// Make an image of multiple lights.
image<vec4f> make_lights(vec2i size, vec3f le, int nlights, float langle,
    float lwidth, float lheight) {
  auto img = image<vec4f>{size};
  for (auto j = 0; j < size.y / 2; j++) {
    auto theta = pif * ((j + 0.5f) / size.y);
    theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
    if (fabs(theta - langle) > lheight / 2) continue;
    for (int i = 0; i < size.x; i++) {
      auto phi     = 2 * pif * (float(i + 0.5f) / size.x);
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
