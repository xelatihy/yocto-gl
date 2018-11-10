//
// Implementation for Yocto/Image.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "yocto_image.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
// Convert between CIE XYZ and RGB
vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

// Convert HSV to RGB
vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0.0f) return {v, v, v};

    h       = fmodf(h, 1.0f) / (60.0f / 360.0f);
    int   i = (int)h;
    float f = h - (float)i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - s * f);
    float t = v * (1.0f - s * (1.0f - f));

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
vec3f rgb_to_hsv(const vec3f& rgb) {
    // from Imgui.cpp
    auto  r = rgb.x, g = rgb.y, b = rgb.z;
    float K = 0.f;
    if (g < b) {
        swap(g, b);
        K = -1.f;
    }
    if (r < g) {
        swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - (g < b ? g : b);
    return {
        fabsf(K + (g - b) / (6.f * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Gets an image size from a suggested size and an aspect ratio. The suggested
// size may have zeros in either components. In which case, we use the aspect
// ration to compute the other.
vec2i get_image_size(const vec2i& size, float aspect) {
    if (size == zero2i) {
        return {(int)round(720 * aspect), 720};
    } else if (size.y == 0) {
        return {size.x, (int)round(size.x / aspect)};
    } else if (size.x == 0) {
        return {(int)round(size.y * aspect), size.y};
    } else {
        return size;
    }
}

// Splits an image into an array of regions
vector<bbox2i> make_image_regions(const vec2i& image_size, int region_size) {
    auto regions = vector<bbox2i>{};
    for (auto y = 0; y < image_size.y; y += region_size) {
        for (auto x = 0; x < image_size.x; x += region_size) {
            regions.push_back({{x, y}, {min(x + region_size, image_size.x),
                                           min(y + region_size, image_size.y)}});
        }
    }
    return regions;
}

// Conversion between linear and gamma-encoded images.
image<vec4f> gamma_to_linear(const image<vec4f>& srgb, float gamma) {
    auto lin = image<vec4f>{srgb.size()};
    for (auto j = 0; j < srgb.size().y; j++) {
        for (auto i = 0; i < srgb.size().x; i++) {
            lin[{i, j}] = gamma_to_linear(srgb[{i, j}], gamma);
        }
    }
    return lin;
}
image<vec4f> linear_to_gamma(const image<vec4f>& lin, float gamma) {
    auto srgb = image<vec4f>{lin.size()};
    for (auto j = 0; j < srgb.size().y; j++) {
        for (auto i = 0; i < srgb.size().x; i++) {
            srgb[{i, j}] = linear_to_gamma(lin[{i, j}], gamma);
        }
    }
    return srgb;
}

// Conversion between linear and gamma-encoded images.
image<vec4f> srgb_to_linear(const image<vec4f>& srgb) {
    auto lin = image<vec4f>{srgb.size()};
    for (auto j = 0; j < srgb.size().y; j++) {
        for (auto i = 0; i < srgb.size().x; i++) {
            lin[{i, j}] = srgb_to_linear(srgb[{i, j}]);
        }
    }
    return lin;
}
image<vec4f> linear_to_srgb(const image<vec4f>& lin) {
    auto srgb = image<vec4f>{lin.size()};
    for (auto j = 0; j < srgb.size().y; j++) {
        for (auto i = 0; i < srgb.size().x; i++) {
            srgb[{i, j}] = linear_to_srgb(lin[{i, j}]);
        }
    }
    return srgb;
}

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt) {
    auto fl = image<vec4f>{bt.size()};
    for (auto j = 0; j < bt.size().y; j++) {
        for (auto i = 0; i < bt.size().x; i++) {
            fl[{i, j}] = byte_to_float(bt[{i, j}]);
        }
    }
    return fl;
}
image<vec4b> float_to_byte(const image<vec4f>& fl) {
    auto bt = image<vec4b>{fl.size()};
    for (auto j = 0; j < bt.size().y; j++) {
        for (auto i = 0; i < bt.size().x; i++) {
            bt[{i, j}] = float_to_byte(fl[{i, j}]);
        }
    }
    return bt;
}

// Tonemap image
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
    auto ldr = image<vec4f>{hdr.size()};
    for (auto j = 0; j < hdr.size().y; j++) {
        for (auto i = 0; i < hdr.size().x; i++) {
            ldr[{i, j}] = tonemap_filmic(hdr[{i, j}], exposure, filmic, srgb);
        }
    }
    return ldr;
}

// Tonemap image
void tonemap_image_region(image<vec4f>& ldr, const bbox2i& region,
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb) {
    for (auto j = region.min.y; j < region.max.y; j++) {
        for (auto i = region.min.x; i < region.max.x; i++) {
            ldr[{i, j}] = tonemap_filmic(hdr[{i, j}], exposure, filmic, srgb);
        }
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a grid image
image<vec4f> make_grid_image(
    const vec2i& size, int tiles, const vec4f& c0, const vec4f& c1) {
    auto img  = image<vec4f>{size};
    auto tile = img.size().x / tiles;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            img[{i, j}] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make a checkerboard image
image<vec4f> make_checker_image(
    const vec2i& size, int tiles, const vec4f& c0, const vec4f& c1) {
    auto img  = image<vec4f>{size};
    auto tile = img.size().x / tiles;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto c      = (i / tile + j / tile) % 2 == 0;
            img[{i, j}] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make an image with bumps and dimples.
image<vec4f> make_bumpdimple_image(const vec2i& size, int tiles) {
    auto img  = image<vec4f>{size};
    auto tile = img.size().x / tiles;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto c  = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r = sqrt(float(ii * ii + jj * jj)) /
                     sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) {
                h += (c) ? (0.5f - r) : -(0.5f - r);
            }
            img[{i, j}] = {h, h, h, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image<vec4f> make_ramp_image(const vec2i& size, const vec4f& c0, const vec4f& c1) {
    auto img = image<vec4f>{size};
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto u      = (float)i / (float)img.size().x;
            img[{i, j}] = c0 * (1 - u) + c1 * u;
        }
    }
    return img;
}

// Make a gamma ramp image
image<vec4f> make_gammaramp_imagef(const vec2i& size) {
    auto img = image<vec4f>{size};
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto u = j / float(img.size().y - 1);
            if (i < img.size().x / 3) u = pow(u, 2.2f);
            if (i > (img.size().x * 2) / 3) u = pow(u, 1 / 2.2f);
            img[{i, j}] = {u, u, u, 1};
        }
    }
    return img;
}

// Make an image color with red/green in the [0,1] range. Helpful to
// visualize uv texture coordinate application.
image<vec4f> make_uvramp_image(const vec2i& size) {
    auto img = image<vec4f>{size};
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            img[{i, j}] = {i / (float)(img.size().x - 1),
                j / (float)(img.size().y - 1), 0, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image<vec4f> make_uvgrid_image(const vec2i& size, int tiles, bool colored) {
    auto img  = image<vec4f>{size};
    auto tile = img.size().x / tiles;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto ii = i / tile, jj = j / tile;
            auto ww = img.size().x / tile, hh = img.size().y / tile;
            auto ph = (((256 / (ww * hh)) * (ii + jj * ww) - 64 + 256) % 256) /
                      360.f;
            auto pv = 0.5f;
            auto ps = 0.8f;
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 0.05f;
                else
                    pv -= 0.05f;
            } else {
                pv = 0.8f;
                ps = 0.2f;
            }
            auto rgb = (colored) ? hsv_to_rgb({ph, ps, pv}) : vec3f{pv, pv, pv};
            img[{i, img.size().y - j - 1}] = {rgb.x, rgb.y, rgb.z, 1};
        }
    }
    return img;
}

// Comvert a bump map to a normal map.
image<vec4f> bump_to_normal_map(const image<vec4f>& img, float scale) {
    auto norm = image<vec4f>{img.size()};
    auto dx = 1.0f / img.size().x, dy = 1.0f / img.size().y;
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto i1 = (i + 1) % img.size().x, j1 = (j + 1) % img.size().y;
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
    return img;
}

// Implementation of sunsky modified heavily from pbrt
image<vec4f> make_sunsky_image(const vec2i& size, float thetaSun,
    float turbidity, bool has_sun, const vec3f& ground_albedo) {
    auto wSun = vec3f{0, cos(thetaSun), sin(thetaSun)};

    // sunSpectralRad =  ComputeAttenuatedSunlight(thetaS, turbidity);
    auto sunAngularRadius = 9.35e-03f / 2;  // Wikipedia
    auto thetaS           = thetaSun;

    auto t1 = thetaSun, t2 = thetaSun * thetaSun,
         t3 = thetaSun * thetaSun * thetaSun;
    auto T  = turbidity;
    auto T2 = turbidity * turbidity;

    auto zenith_xyY = vec3f{
        (+0.00165f * t3 - 0.00374f * t2 + 0.00208f * t1 + 0) * T2 +
            (-0.02902f * t3 + 0.06377f * t2 - 0.03202f * t1 + 0.00394f) * T +
            (+0.11693f * t3 - 0.21196f * t2 + 0.06052f * t1 + 0.25885f),
        (+0.00275f * t3 - 0.00610f * t2 + 0.00316f * t1 + 0) * T2 +
            (-0.04214f * t3 + 0.08970f * t2 - 0.04153f * t1 + 0.00515f) * T +
            (+0.15346f * t3 - 0.26756f * t2 + 0.06669f * t1 + 0.26688f),
        1000 * (4.0453f * T - 4.9710f) *
                tan((4.0f / 9.0f - T / 120.0f) * (pif - 2 * t1)) -
            .2155f * T + 2.4192f};

    auto perez_A_xyY = vec3f{-0.01925f * T - 0.25922f, -0.01669f * T - 0.26078f,
        +0.17872f * T - 1.46303f};
    auto perez_B_xyY = vec3f{-0.06651f * T + 0.00081f, -0.09495f * T + 0.00921f,
        -0.35540f * T + 0.42749f};
    auto perez_C_xyY = vec3f{-0.00041f * T + 0.21247f, -0.00792f * T + 0.21023f,
        -0.02266f * T + 5.32505f};
    auto perez_D_xyY = vec3f{-0.06409f * T - 0.89887f, -0.04405f * T - 1.65369f,
        +0.12064f * T - 2.57705f};
    auto perez_E_xyY = vec3f{-0.00325f * T + 0.04517f, -0.01092f * T + 0.05291f,
        -0.06696f * T + 0.37027f};

    auto perez_f = [thetaS](float A, float B, float C, float D, float E,
                       float theta, float gamma, float zenith) -> float {
        auto den = ((1 + A * exp(B)) *
                    (1 + C * exp(D * thetaS) + E * cos(thetaS) * cos(thetaS)));
        auto num = ((1 + A * exp(B / cos(theta))) *
                    (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
        return zenith * num / den;
    };

    auto sky = [&perez_f, perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                   perez_E_xyY, zenith_xyY](auto theta, auto gamma) -> vec3f {
        auto x = perez_f(perez_A_xyY.x, perez_B_xyY.x, perez_C_xyY.x,
            perez_D_xyY.x, perez_E_xyY.x, theta, gamma, zenith_xyY.x);
        auto y = perez_f(perez_A_xyY.y, perez_B_xyY.y, perez_C_xyY.y,
            perez_D_xyY.y, perez_E_xyY.y, theta, gamma, zenith_xyY.y);
        auto Y = perez_f(perez_A_xyY.z, perez_B_xyY.z, perez_C_xyY.z,
            perez_D_xyY.z, perez_E_xyY.z, theta, gamma, zenith_xyY.z);
        return xyz_to_rgb(xyY_to_xyz({x, y, Y})) / 10000.0f;
    };

    // compute sun luminance
    // TODO: how this relates to zenith intensity?
    auto sun_ko     = vec3f{0.48f, 0.75f, 0.14f};
    auto sun_kg     = vec3f{0.1f, 0.0f, 0.0f};
    auto sun_kwa    = vec3f{0.02f, 0.0f, 0.0f};
    auto sun_sol    = vec3f{20000.0f, 27000.0f, 30000.0f};
    auto sun_lambda = vec3f{680, 530, 480};
    auto sun_beta   = 0.04608365822050f * turbidity - 0.04586025928522f;
    auto sun_m      = 1.0f /
                 (cos(thetaSun) + 0.000940f * pow(1.6386f - thetaSun, -1.253f));

    auto sun_le = zero3f;
    for (auto i = 0; i < 3; i++) {
        auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda[i] / 1000, -4.08f));
        auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda[i] / 1000, -1.3f));
        auto tauO = exp(-sun_m * sun_ko[i] * .35f);
        auto tauG = exp(-1.41f * sun_kg[i] * sun_m /
                        pow(1 + 118.93f * sun_kg[i] * sun_m, 0.45f));
        auto tauWA = exp(-0.2385f * sun_kwa[i] * 2.0f * sun_m /
                         pow(1 + 20.07f * sun_kwa[i] * 2.0f * sun_m, 0.45f));
        sun_le[i]  = sun_sol[i] * tauR * tauA * tauO * tauG * tauWA;
    }

    auto sun = [has_sun, sunAngularRadius, sun_le](auto theta, auto gamma) {
        return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f : zero3f;
    };

    auto img = image<vec4f>{size, {0, 0, 0, 1}};
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = pif * ((j + 0.5f) / img.size().y);
        theta      = clamp(theta, 0.0f, pif / 2 - epsf);
        for (int i = 0; i < img.size().x; i++) {
            auto phi = 2 * pif * (float(i + 0.5f) / img.size().x);
            auto w   = vec3f{
                cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma  = acos(clamp(dot(w, wSun), -1.0f, 1.0f));
            auto col    = sky(theta, gamma) + sun(theta, gamma);
            img[{i, j}] = {col.x, col.y, col.z, 1};
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < img.size().y / 2; j++) {
            auto theta = pif * ((j + 0.5f) / img.size().y);
            for (int i = 0; i < img.size().x; i++) {
                auto pxl   = img[{i, j}];
                auto le    = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pif / (img.size().x * img.size().y);
                ground += le * (ground_albedo / pif) * cos(theta) * angle;
            }
        }
        for (auto j = img.size().y / 2; j < img.size().y; j++) {
            for (int i = 0; i < img.size().x; i++) {
                img[{i, j}] = {ground.x, ground.y, ground.z, 1};
            }
        }
    }
    return img;
}

// Make an image of multiple lights.
image<vec4f> make_lights_image(const vec2i& size, const vec3f& le, int nlights,
    float langle, float lwidth, float lheight) {
    auto img = image<vec4f>{size, {0, 0, 0, 1}};
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = pif * ((j + 0.5f) / img.size().y);
        theta      = clamp(theta, 0.0f, pif / 2 - epsf);
        if (fabs(theta - langle) > lheight / 2) continue;
        for (int i = 0; i < img.size().x; i++) {
            auto phi     = 2 * pif * (float(i + 0.5f) / img.size().x);
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

// Make a noise image. Wrap works only if size is a power of two.
image<vec4f> make_noise_image(const vec2i& size, float scale, bool wrap) {
    auto img    = image<vec4f>{size};
    auto wrap3i = (wrap) ? vec3i{img.size().x, img.size().y, 2} : zero3i;
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto p = vec3f{i / (float)img.size().x, j / (float)img.size().y,
                         0.5f} *
                     scale;
            auto g      = perlin_noise(p, wrap3i);
            g           = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if size is a power of two.
image<vec4f> make_fbm_image(const vec2i& size, float scale, float lacunarity,
    float gain, int octaves, bool wrap) {
    auto img    = image<vec4f>{size};
    auto wrap3i = (wrap) ? vec3i{img.size().x, img.size().y, 2} : zero3i;
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto p = vec3f{i / (float)img.size().x, j / (float)img.size().y,
                         0.5f} *
                     scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g      = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if size is a power of two.
image<vec4f> make_ridge_image(const vec2i& size, float scale, float lacunarity,
    float gain, float offset, int octaves, bool wrap) {
    auto img    = image<vec4f>{size};
    auto wrap3i = (wrap) ? vec3i{img.size().x, img.size().y, 2} : zero3i;
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto p = vec3f{i / (float)img.size().x, j / (float)img.size().y,
                         0.5f} *
                     scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g           = clamp(g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if size is a power of two.
image<vec4f> make_turbulence_image(const vec2i& size, float scale,
    float lacunarity, float gain, int octaves, bool wrap) {
    auto img    = image<vec4f>{size};
    auto wrap3i = (wrap) ? vec3i{img.size().x, img.size().y, 2} : zero3i;
    for (auto j = 0; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            auto p = vec3f{i / (float)img.size().x, j / (float)img.size().y,
                         0.5f} *
                     scale;
            auto g = perlin_turbulence_noise(
                p, lacunarity, gain, octaves, wrap3i);
            g           = clamp(g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// make a simple example volume
volume<float> make_test_volume(const vec3i& size, float scale, float exponent) {
    auto vol = volume<float>{size};
    for (auto k = 0; k < vol.size().z; k++) {
        for (auto j = 0; j < vol.size().y; j++) {
            for (auto i = 0; i < vol.size().x; i++) {
                auto p = vec3f{
                    i / (float)size.x, j / (float)size.y, k / (float)size.z};
                float value = pow(
                    max(max(cos(scale * p.x), cos(scale * p.y)), 0.0f), exponent);
                vol[{i, j, k}] = clamp(value, 0.0f, 1.0f);
            }
        }
    }
    return vol;
}

}  // namespace yocto
