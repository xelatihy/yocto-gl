//
// Implementation for Yocto/Image.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// #ifndef _clang_analyzer__

#define STB_IMAGE_IMPLEMENTATION
#include "ext/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "ext/stb_image_resize.h"

#define TINYEXR_IMPLEMENTATION
#include "ext/tinyexr.h"

// #endif

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic pop
#endif

#include "ext/ArHosekSkyModel.cpp"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make an image by assign values to each pixel
template <typename T, typename Func>
void make_fromij(image<T>& img, const vec2i& size, const Func& func) {
    img.resize(size);
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            img[{i, j}] = func(i, j);
        }
    }
}
template <typename T, typename Func>
void make_fromuv(image<T>& img, const vec2i& size, const Func& func) {
    img.resize(size);
    for (int j = 0; j < img.size().y; j++) {
        for (int i = 0; i < img.size().x; i++) {
            auto u      = (float)i / (float)img.size().x;
            auto v      = (float)j / (float)img.size().y;
            img[{i, j}] = func(u, v);
        }
    }
}

// Make a grid image
void make_grid(image<vec4f>& img, const vec2i& size, int tiles, const vec4f& c0,
    const vec4f& c1) {
    make_fromij(img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
        auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                 j % tile == tile - 1;
        return (c) ? c0 : c1;
    });
}

// Make a checkerboard image
void make_checker(image<vec4f>& img, const vec2i& size, int tiles,
    const vec4f& c0, const vec4f& c1) {
    make_fromij(img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
        auto c = (i / tile + j / tile) % 2 == 0;
        return (c) ? c0 : c1;
    });
}

// Make an image with bumps and dimples.
void make_bumpdimple(image<vec4f>& img, const vec2i& size, int tiles,
    const vec4f& c0, const vec4f& c1) {
    make_fromij(img, size, [tile = size.x / tiles, &c0, &c1](int i, int j) {
        auto c  = (i / tile + j / tile) % 2 == 0;
        auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
        auto r = sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
        auto h = 0.5f;
        if (r < 0.5f) {
            h += (c) ? (0.5f - r) : -(0.5f - r);
        }
        return lerp(c0, c1, h);
    });
}

// Make a uv colored grid
void make_ramp(
    image<vec4f>& img, const vec2i& size, const vec4f& c0, const vec4f& c1) {
    make_fromij(img, size, [size = img.size(), &c0, &c1](int i, int j) {
        auto u = (float)i / (float)size.x;
        return lerp(c0, c1, u);
    });
}
void make_ramp(image<vec4f>& img, const vec2i& size, const vec4f& c00,
    const vec4f& c10, const vec4f& c11, const vec4f& c01) {
    make_fromij(img, size, [size, &c00, &c10, &c01, &c11](int i, int j) {
        auto u = (float)i / (float)size.x;
        auto v = (float)j / (float)size.y;
        return bilerp(c00, c10, c11, c01, u, v);
    });
}

// Make a gamma ramp image
void make_gammaramp(
    image<vec4f>& img, const vec2i& size, const vec4f& c0, const vec4f& c1) {
    make_fromij(img, size, [size, &c0, &c1](int i, int j) {
        auto u = j / float(size.y - 1);
        if (i < size.x / 3) u = pow(u, 2.2f);
        if (i > (size.x * 2) / 3) u = pow(u, 1 / 2.2f);
        return lerp(c0, c1, u);
    });
}

// Make an image color with red/green in the [0,1] range. Helpful to
// visualize uv texture coordinate application.
void make_uvramp(image<vec4f>& img, const vec2i& size) {
    return make_ramp(
        img, size, {0, 0, 0, 0}, {1, 0, 0, 0}, {1, 1, 0, 0}, {0, 1, 0, 0});
}

// Make a uv colored grid
void make_uvgrid(
    image<vec4f>& img, const vec2i& size, int tiles, bool colored) {
    make_fromij(
        img, size, [size, tile = size.x / tiles, colored](int i, int j) {
            j       = size.y - j - 1;
            auto ii = i / tile, jj = j / tile;
            auto ww = size.x / tile, hh = size.y / tile;
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
            auto rgb = (colored) ? hsv_to_rgb(vec3f{ph, ps, pv})
                                 : vec3f{pv, pv, pv};
            return vec4f{rgb, 1};
        });
}

// Makes a blackbody ramp
void make_blackbodyramp(image<vec4f>& img, const vec2i& size,
    float start_temperature, float end_temperature) {
    make_fromij(
        img, size, [size, start_temperature, end_temperature](int i, int j) {
            auto temperature = start_temperature +
                               (end_temperature - start_temperature) *
                                   (float)i / (float)(size.x - 1);
            auto rgb = blackbody_to_rgb(temperature);
            return vec4f{rgb, 1};
        });
}

// Make a noise image. Wrap works only if size is a power of two.
void make_noise(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale, bool wrap) {
    make_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale, &c0,
            &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g      = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            return lerp(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
void make_fbm(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale, float lacunarity, float gain, int octaves, bool wrap) {
    make_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g      = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            return lerp(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
void make_ridge(image<vec4f>& img, const vec2i& size, const vec4f& c0,
    const vec4f& c1, float scale, float lacunarity, float gain, float offset, int octaves,
    bool wrap) {
    make_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, offset, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            return lerp(c0, c1, g);
        });
}

// Make a noise image. Wrap works only if size is a power of two.
void make_turbulence(image<vec4f>& img, const vec2i& size,
    const vec4f& c0, const vec4f& c1, float scale, float lacunarity, float gain,
    int octaves, bool wrap) {
    make_fromij(img, size,
        [wrap3i = (wrap) ? vec3i{size.x, size.y, 2} : zero3i, size, scale,
            lacunarity, gain, octaves, &c0, &c1](int i, int j) {
            auto p = vec3f{i / (float)size.x, j / (float)size.y, 0.5f} * scale;
            auto g = perlin_turbulence_noise(
                p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            return lerp(c0, c1, g);
        });
}

// Comvert a bump map to a normal map.
void bump_to_normal(
    image<vec4f>& norm, const image<vec4f>& img, float scale) {
    norm.resize(img.size());
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
}

// Add a border to an image
void add_border(
    image<vec4f>& img, int border_width, const vec4f& border_color) {
    for (auto j = 0; j < img.size().y; j++) {
        for (auto b = 0; b < border_width; b++) {
            img[{b, j}]                    = border_color;
            img[{img.size().x - 1 - b, j}] = border_color;
        }
    }
    for (auto i = 0; i < img.size().x; i++) {
        for (auto b = 0; b < border_width; b++) {
            img[{i, b}]                    = border_color;
            img[{i, img.size().y - 1 - b}] = border_color;
        }
    }
}

#if 1

// Implementation of sunsky modified heavily from pbrt
void make_sunsky(image<vec4f>& img, const vec2i& size, float theta_sun,
    float turbidity, bool has_sun, float sun_intensity, float sun_temperature,
    const vec3f& ground_albedo) {
    // idea adapted from pbrt

    // initialize model
    double wavelengths[9] = {630, 680, 710, 500, 530, 560, 460, 480, 490};
    ArHosekSkyModelState* skymodel_state[9];
    if (sun_temperature) {
        sun_temperature = clamp(sun_temperature, 2000.0f, 14000.0f);
        for (int i = 0; i < 9; ++i) {
            skymodel_state[i] = arhosekskymodelstate_alienworld_alloc_init(
                theta_sun, sun_intensity, sun_temperature, turbidity,
                ground_albedo[i / 3]);
        }
    } else {
        for (int i = 0; i < 9; ++i) {
            skymodel_state[i] = arhosekskymodelstate_alloc_init(
                theta_sun, turbidity, ground_albedo[i / 3]);
        }
    }

    // clear image
    img.resize(size);
    for (auto& p : img) p = {0, 0, 0, 1};

    // sun-sky
    auto sun_direction = vec3f{0, sin(theta_sun), cos(theta_sun)};
    auto integral      = zero3f;
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = (j + 0.5f) * pif / img.size().y;
        if (theta > pif / 2) continue;
        for (auto i = 0; i < img.size().x; i++) {
            auto phi       = (i + 0.5f) * 2 * pif / img.size().x;
            auto direction = vec3f{
                cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(clamp(dot(direction, sun_direction), -1.0f, 1.0f));
            for (int c = 0; c < 9; ++c) {
                auto val =
                    (has_sun)
                        ? arhosekskymodel_solar_radiance(
                              skymodel_state[c], theta, gamma, wavelengths[c])
                        : arhosekskymodel_radiance(
                              skymodel_state[c], theta, gamma, wavelengths[c]);
                // average channel over wavelengths
                img[{i, j}][c / 3] += (float)val / 3;
            }
            integral += xyz(img[{i, j}]) * sin(theta) /
                        (img.size().x * img.size().y / 2);
        }
    }

    // ground
    auto ground = ground_albedo * integral;
    for (auto j = img.size().y / 2; j < img.size().y; j++) {
        for (auto i = 0; i < img.size().x; i++) {
            img[{i, j}] = {ground.x, ground.y, ground.z, 1};
        }
    }

    // cleanup
    for (auto i = 0; i < 9; i++) arhosekskymodelstate_free(skymodel_state[i]);
}

#else

// Implementation of sunsky modified heavily from pbrt
image<vec4f> make_sunsky(int width, int height, float theta_sun,
    float turbidity, bool has_sun, float sun_angle_scale,
    float sun_emission_scale, const vec3f& ground_albedo,
    bool renormalize_sun) {
    auto zenith_xyY = vec3f{
        (+0.00165f * pow(theta_sun, 3.f) - 0.00374f * pow(theta_sun, 2.f) +
            0.00208f * theta_sun + 0) *
                pow(turbidity, 2.f) +
            (-0.02902f * pow(theta_sun, 3.f) + 0.06377f * pow(theta_sun, 2.f) -
                0.03202f * theta_sun + 0.00394f) *
                turbidity +
            (+0.11693f * pow(theta_sun, 3.f) - 0.21196f * pow(theta_sun, 2.f) +
                0.06052f * theta_sun + 0.25885f),
        (+0.00275f * pow(theta_sun, 3.f) - 0.00610f * pow(theta_sun, 2.f) +
            0.00316f * theta_sun + 0) *
                pow(turbidity, 2.f) +
            (-0.04214f * pow(theta_sun, 3.f) + 0.08970f * pow(theta_sun, 2.f) -
                0.04153f * theta_sun + 0.00515f) *
                turbidity +
            (+0.15346f * pow(theta_sun, 3.f) - 0.26756f * pow(theta_sun, 2.f) +
                0.06669f * theta_sun + 0.26688f),
        1000 * (4.0453f * turbidity - 4.9710f) *
                tan((4.0f / 9.0f - turbidity / 120.0f) *
                    (pif - 2 * theta_sun)) -
            .2155f * turbidity + 2.4192f};

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
        auto den =
            ((1 + A * exp(B)) * (1 + C * exp(D * theta_sun) +
                                    E * cos(theta_sun) * cos(theta_sun)));
        auto num = ((1 + A * exp(B / cos(theta))) *
                    (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
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
    // TODO: how this relates to zenith intensity?
    auto sun_ko     = vec3f{0.48f, 0.75f, 0.14f};
    auto sun_kg     = vec3f{0.1f, 0.0f, 0.0f};
    auto sun_kwa    = vec3f{0.02f, 0.0f, 0.0f};
    auto sun_sol    = vec3f{20000.0f, 27000.0f, 30000.0f};
    auto sun_lambda = vec3f{680, 530, 480};
    auto sun_beta   = 0.04608365822050f * turbidity - 0.04586025928522f;
    auto sun_m =
        1.0f / (cos(theta_sun) + 0.000940f * pow(1.6386f - theta_sun, -1.253f));

    auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda / 1000, -4.08f));
    auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda / 1000, -1.3f));
    auto tauO = exp(-sun_m * sun_ko * .35f);
    auto tauG = exp(
        -1.41f * sun_kg * sun_m / pow(1 + 118.93f * sun_kg * sun_m, 0.45f));
    auto tauWA  = exp(-0.2385f * sun_kwa * 2.0f * sun_m /
                     pow(1 + 20.07f * sun_kwa * 2.0f * sun_m, 0.45f));
    auto sun_le = sun_sol * tauR * tauA * tauO * tauG * tauWA;

    // rescale by user
    sun_le *= sun_emission_scale;

    // sun scale from Wikipedia scaled by user quantity and rescaled to at
    // the minimum 5 pixel diamater
    auto sun_angular_radius = 9.35e-03f / 2;  // Wikipedia
    sun_angular_radius *= sun_angle_scale;
    sun_angular_radius = max(sun_angular_radius, 5 * pif / height);

    // sun direction
    auto sun_direction = vec3f{0, cos(theta_sun), sin(theta_sun)};

    auto sun = [has_sun, sun_angular_radius, sun_le](auto theta, auto gamma) {
        // return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f :
        //                                                zero3f;
        return (has_sun && gamma < sun_angular_radius) ? sun_le / 10000
                                                       : zero3f;
    };

    // Make the sun sky image
    auto img          = make_image(width, height, vec4f{0, 0, 0, 1});
    auto sky_integral = 0.0f, sun_integral = 0.0f;
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = pif * ((j + 0.5f) / img.size().y);
        theta      = clamp(theta, 0.0f, pif / 2 - float_epsilon);
        for (int i = 0; i < img.size().x; i++) {
            auto phi = 2 * pif * (float(i + 0.5f) / img.size().x);
            auto w   = vec3f{
                cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma   = acos(clamp(dot(w, sun_direction), -1.0f, 1.0f));
            auto sky_col = sky(theta, gamma, theta_sun);
            auto sun_col = sun(theta, gamma);
            sky_integral += mean(sky_col) * sin(theta);
            sun_integral += mean(sun_col) * sin(theta);
            auto col    = sky_col + sun_col;
            img[{i, j}] = {col.x, col.y, col.z, 1};
        }
    }

    if (renormalize_sun) {
        for (auto j = 0; j < img.size().y / 2; j++) {
            for (int i = 0; i < img.size().x; i++) {
                img[{i, j}] *= sky_integral / (sun_integral + sky_integral);
            }
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < img.size().y / 2; j++) {
            auto theta = pif * ((j + 0.5f) / img.size().y);
            for (int i = 0; i < img.size().x; i++) {
                auto pxl   = img[{i, j}];
                auto le    = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pif /
                             (img.size().x * img.size().y);
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

#endif

// Make an image of multiple lights.
void make_lights(image<vec4f>& img, const vec2i& size, const vec3f& le,
    int nlights, float langle, float lwidth, float lheight) {
    img.resize(size);
    for (auto j = 0; j < img.size().y / 2; j++) {
        auto theta = pif * ((j + 0.5f) / img.size().y);
        theta      = clamp(theta, 0.0f, pif / 2 - 0.00001f);
        if (fabs(theta - langle) > lheight / 2) continue;
        for (int i = 0; i < img.size().x; i++) {
            auto phi     = 2 * pif * (float(i + 0.5f) / img.size().x);
            auto inlight = false;
            for (auto l = 0; l < nlights; l++) {
                auto lphi = 2 * pif * (l + 0.5f) / nlights;
                inlight   = inlight || fabs(phi - lphi) < lwidth / 2;
            }
            img[{i, j}] = {le, 1};
        }
    }
}

void make_logo(image<vec4b>& img, const string& type) {
    static const auto size = vec2i{144, 28};
    // clang-format off
    static const auto logo_render = vector<byte>{
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 212, 87, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 255, 62, 0, 0, 0, 0, 0, 14, 27, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 187, 245, 13, 0, 0, 0, 80, 255, 90, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29, 88, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 88, 251, 10, 0, 0, 0, 40, 200, 253, 255, 234, 106, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 12, 147, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 90, 69, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 144, 125, 0, 0, 0, 0, 0, 117, 37, 0, 0, 0, 0, 0, 0, 0, 79, 255, 101, 0, 0, 0, 178, 232, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 145, 205, 0, 0, 0, 47, 239, 210, 74, 57, 144, 232, 24, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 251, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 146, 123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 43, 35, 0, 61, 87, 0, 0, 208, 61, 0, 0, 0, 0, 0, 0, 0, 3, 224, 199, 0, 0, 24, 251, 129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 201, 149, 0, 0, 0, 180, 238, 20, 0, 0, 0, 16, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 251, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 146, 123, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 119, 150, 0, 0, 208, 61, 0, 0, 0, 0, 0, 0, 0, 0, 118, 255, 41, 0, 117, 250, 26, 0, 7, 147, 223, 232, 167, 19, 0, 0, 0, 5, 143, 224, 234, 162, 20, 166, 227, 255, 210, 201, 3, 0, 7, 147, 223, 232, 167, 19, 0, 0, 0, 0, 8, 249, 93, 0, 0, 45, 255, 146, 0, 0, 0, 0, 0, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 14, 192, 114, 207, 0, 84, 224, 220, 61, 0, 68, 143, 166, 220, 78, 0, 0, 142, 234, 160, 251, 0, 0, 134, 234, 195, 25, 0, 123, 105, 211, 89, 12, 177, 236, 158, 4, 0, 44, 217, 214, 189, 123, 0, 0, 0, 149, 76, 0, 189, 108, 0, 160, 55, 123, 107, 83, 236, 240, 179, 0, 208, 128, 220, 174, 6, 0, 0, 0, 0, 0, 19, 247, 140, 0, 214, 168, 0, 0, 176, 248, 122, 111, 237, 212, 2, 0, 0, 167, 251, 128, 127, 223, 59, 103, 185, 255, 143, 117, 0, 0, 176, 248, 122, 111, 237, 212, 2, 0, 0, 0, 58, 255, 36, 0, 0, 97, 255, 77, 0, 0, 0, 0, 0, 0, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 248, 166, 36, 13, 234, 44, 73, 215, 0, 80, 248, 65, 80, 193, 0, 66, 222, 32, 97, 251, 0, 65, 207, 21, 135, 152, 0, 144, 230, 81, 12, 129, 157, 17, 189, 88, 0, 194, 126, 17, 208, 123, 0, 0, 0, 131, 125, 7, 228, 164, 0, 224, 23, 144, 125, 5, 127, 156, 10, 0, 208, 179, 13, 205, 65, 0, 0, 0, 0, 0, 0, 158, 233, 61, 255, 59, 0, 46, 255, 121, 0, 0, 91, 255, 83, 0, 40, 254, 134, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 46, 255, 121, 0, 0, 91, 255, 83, 0, 0, 0, 114, 235, 0, 0, 0, 130, 255, 51, 0, 2, 17, 17, 17, 7, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 255, 43, 0, 76, 191, 0, 1, 242, 20, 80, 192, 0, 36, 235, 0, 138, 142, 0, 18, 251, 0, 140, 127, 0, 52, 212, 0, 144, 171, 0, 0, 204, 63, 0, 116, 148, 11, 255, 14, 0, 146, 123, 0, 0, 0, 84, 164, 46, 155, 201, 9, 231, 0, 144, 125, 0, 119, 150, 0, 0, 208, 64, 0, 164, 107, 0, 0, 0, 0, 0, 0, 49, 255, 220, 206, 0, 0, 114, 255, 46, 0, 0, 17, 255, 148, 0, 111, 255, 54, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 114, 255, 46, 0, 0, 17, 255, 148, 0, 0, 0, 171, 179, 0, 0, 0, 159, 255, 29, 0, 20, 255, 255, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 101, 242, 214, 214, 249, 40, 80, 189, 0, 34, 236, 0, 162, 121, 0, 18, 251, 0, 165, 232, 214, 219, 232, 0, 144, 126, 0, 0, 229, 222, 214, 229, 168, 34, 249, 0, 0, 146, 123, 0, 0, 0, 38, 203, 88, 107, 206, 48, 187, 0, 144, 125, 0, 119, 150, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 197, 255, 98, 0, 0, 144, 255, 22, 0, 0, 0, 249, 177, 0, 141, 255, 28, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 144, 255, 22, 0, 0, 0, 249, 177, 0, 0, 0, 227, 123, 0, 0, 0, 143, 255, 40, 0, 0, 79, 108, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 85, 188, 4, 4, 4, 0, 80, 189, 0, 34, 236, 0, 149, 132, 0, 18, 251, 0, 149, 125, 4, 4, 3, 0, 144, 125, 0, 0, 213, 62, 4, 4, 2, 21, 255, 4, 0, 146, 123, 0, 0, 0, 2, 230, 130, 67, 178, 116, 141, 0, 144, 125, 0, 119, 150, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 151, 255, 17, 0, 0, 0, 245, 183, 0, 152, 255, 21, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 151, 255, 17, 0, 0, 0, 245, 183, 0, 0, 28, 255, 67, 0, 0, 0, 116, 255, 59, 0, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 32, 236, 19, 2, 33, 0, 80, 189, 0, 34, 236, 0, 98, 195, 4, 64, 251, 0, 93, 191, 4, 11, 24, 0, 144, 125, 0, 0, 157, 131, 0, 27, 8, 1, 225, 72, 1, 190, 123, 0, 0, 0, 0, 201, 196, 26, 137, 197, 96, 0, 144, 125, 0, 109, 159, 0, 0, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 122, 255, 36, 0, 0, 8, 255, 152, 0, 124, 255, 42, 0, 0, 0, 0, 0, 115, 255, 33, 0, 0, 122, 255, 36, 0, 0, 8, 255, 152, 0, 0, 84, 253, 13, 0, 0, 0, 79, 255, 108, 0, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 16, 253, 0, 0, 0, 123, 238, 224, 155, 0, 80, 189, 0, 34, 236, 0, 11, 201, 218, 156, 248, 0, 6, 178, 225, 234, 97, 0, 144, 125, 0, 0, 31, 216, 214, 228, 49, 0, 89, 244, 198, 179, 123, 0, 0, 0, 0, 154, 239, 0, 96, 253, 50, 0, 144, 125, 0, 34, 224, 201, 25, 208, 61, 0, 162, 108, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 70, 255, 96, 0, 0, 67, 255, 97, 0, 76, 255, 104, 0, 0, 0, 0, 0, 113, 255, 35, 0, 0, 70, 255, 96, 0, 0, 67, 255, 97, 0, 0, 141, 210, 0, 0, 0, 0, 5, 230, 200, 2, 0, 0, 39, 255, 109, 0, 6, 255, 157, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 14, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 18, 0, 0, 0, 0, 0, 19, 6, 0, 0, 0, 0, 0, 0, 0, 0, 23, 1, 0, 0, 0, 12, 6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 1, 211, 223, 54, 44, 205, 229, 7, 0, 2, 216, 230, 69, 74, 169, 31, 0, 55, 255, 120, 59, 26, 1, 211, 223, 54, 44, 205, 229, 7, 0, 0, 197, 154, 0, 0, 0, 0, 0, 109, 255, 166, 59, 78, 165, 255, 109, 0, 6, 255, 201, 113, 113, 113, 105, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 255, 33, 0, 0, 0, 30, 204, 255, 255, 214, 40, 0, 0, 0, 35, 208, 255, 255, 212, 47, 0, 1, 174, 252, 243, 102, 0, 30, 204, 255, 255, 214, 40, 0, 0, 6, 247, 97, 0, 0, 0, 0, 0, 0, 103, 235, 255, 255, 242, 146, 24, 0, 6, 255, 255, 255, 255, 255, 213, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 29, 0, 0, 0, 0, 0, 0, 0, 25, 31, 0, 0, 0, 0, 0, 23, 9, 0, 0, 0, 0, 24, 29, 0, 0, 0, 0, 54, 255, 41, 0, 0, 0, 0, 0, 0, 0, 0, 33, 30, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 63, 188, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
    };
    // clang-format on
    if (type == "logo-render") {
        auto img1 = image<vec<byte, 1>>{size, (vec1b*)logo_render.data()};
        img.resize(size);
        convert_channels(img, img1);
    } else {
        throw io_error("unknown builtin image " + type);
    }
}

void make_logo(image<vec4f>& img, const string& type) {
    auto img8 = image<vec4b>();
    make_logo(img8, type);
    img.resize(img8.size());
    srgb_to_linear(img, img8);
}

void make_preset(image<vec4f>& img, const string& type) {
    auto size = vec2i{1024, 1024};
    if (type.find("sky") != type.npos) size = {2048, 1024};
    if (type == "grid") {
        make_grid(img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "checker") {
        make_checker(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "bump") {
        make_bumpdimple(img, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "uvramp") {
        make_uvramp(img, size);
    } else if (type == "gammaramp") {
        make_gammaramp(img, size, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "blackbodyramp") {
        make_blackbodyramp(img, size);
    } else if (type == "uvgrid") {
        make_uvgrid(img, size);
    } else if (type == "sky") {
        make_sunsky(img, size, pif / 4, 3.0f, false, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "sunsky") {
        make_sunsky(img, size, pif / 4, 3.0f, true, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "noise") {
        make_noise(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, true);
    } else if (type == "fbm") {
        make_fbm(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "ridge") {
        make_ridge(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f,
            1.0f, 6, true);
    } else if (type == "turbulence") {
        make_turbulence(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "bump-normal") {
        auto bump = image<vec4f>{};
        make_bumpdimple(bump, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
        bump_to_normal(img, bump, 0.05f);
    } else if (type == "logo-render") {
        make_logo(img, "logo-render");
    } else if (type == "images1") {
        auto sub_types = vector<string>{"grid", "uvgrid", "checker",
            "gammaramp", "bump", "bump-normal", "noise", "fbm",
            "blackbodyramp"};
        auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
        for (auto i = 0; i < sub_imgs.size(); i++) {
            sub_imgs.at(i).resize(img.size());
            make_preset(sub_imgs.at(i), sub_types.at(i));
        }
        auto montage_size = zero2i;
        for (auto& sub_img : sub_imgs) {
            montage_size.x += sub_img.size().x;
            montage_size.y = max(montage_size.y, sub_img.size().y);
        }
        img.resize(montage_size);
        auto pos = 0;
        for (auto& sub_img : sub_imgs) {
            set_region(img, sub_img, {pos, 0});
            pos += sub_img.size().x;
        }
    } else if (type == "images2") {
        auto sub_types = vector<string>{"sky", "sunsky"};
        auto sub_imgs  = vector<image<vec4f>>(sub_types.size());
        for (auto i = 0; i < sub_imgs.size(); i++) {
            make_preset(sub_imgs.at(i), sub_types.at(i));
        }
        auto montage_size = zero2i;
        for (auto& sub_img : sub_imgs) {
            montage_size.x += sub_img.size().x;
            montage_size.y = max(montage_size.y, sub_img.size().y);
        }
        img.resize(montage_size);
        auto pos = 0;
        for (auto& sub_img : sub_imgs) {
            set_region(img, sub_img, {pos, 0});
            pos += sub_img.size().x;
        }
    } else if (type == "test-floor") {
        make_grid(img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
        add_border(img, 2, {0, 0, 0, 1});
    } else if (type == "test-grid") {
        make_grid(img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "test-checker") {
        make_checker(
            img, size, 8, {0.2f, 0.2f, 0.2f, 1}, {0.5f, 0.5f, 0.5f, 1});
    } else if (type == "test-bump") {
        make_bumpdimple(img, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "test-uvramp") {
        make_uvramp(img, size);
    } else if (type == "test-gammaramp") {
        make_gammaramp(img, size, {0, 0, 0, 1}, {1, 1, 1, 1});
    } else if (type == "test-blackbodyramp") {
        make_blackbodyramp(img, size);
    } else if (type == "test-uvgrid") {
        make_uvgrid(img, size);
    } else if (type == "test-sky") {
        make_sunsky(img, size, pif / 4, 3.0f, false, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "test-sunsky") {
        make_sunsky(img, size, pif / 4, 3.0f, true, 1.0f, 0.0f,
            vec<float, 3>{0.7f, 0.7f, 0.7f});
    } else if (type == "test-noise") {
        make_noise(img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, true);
    } else if (type == "test-fbm") {
        make_fbm(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 1.0f, 2.0f, 0.5f, 6, true);
    } else if (type == "test-bump-normal") {
        auto bump = image<vec4f>{};
        make_bumpdimple(bump, size, 8, {0, 0, 0, 1}, {1, 1, 1, 1});
        bump_to_normal(img, bump, 0.05f);
    } else if (type == "test-fbm-displacement") {
        make_fbm(
            img, size, {0, 0, 0, 1}, {1, 1, 1, 1}, 10.0f, 2.0f, 0.5f, 6, true);
    } else {
        throw std::invalid_argument("unknown image preset " + type);
    }
}

void make_preset(image<vec<byte, 4>>& img, const string& type) {
    auto imgf = image<vec4f>{};
    make_preset(imgf, type);
    if (type.find("-normal") == type.npos) {
        linear_to_srgb(img, imgf);
    } else {
        float_to_byte(img, imgf);
    }
}

void make_preset(
    image<vec4f>& hdr, image<vec4b>& ldr, const string& type) {
    if (type.find("sky") == type.npos) {
        make_preset(ldr, type);
    } else {
        make_preset(hdr, type);
    }
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
                auto p = vec3f{i / (float)vol.size().x, j / (float)vol.size().y,
                    k / (float)vol.size().z};
                auto value = pow(
                    max(max(cos(scale * p.x), cos(scale * p.y)), 0.0f),
                    exponent);
                vol[{i, j, k}] = clamp(value, 0.0f, 1.0f);
            }
        }
    }
}

void make_preset(volume<float>& vol, const string& type) {
    auto size = vec3i{256, 256, 256};
    if (type == "test-volume") {
        make_test(vol, size, 6, 10);
    } else {
        throw runtime_error("unknown volume preset " + type);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace yocto {

// Split a string
static inline vector<string> split_string(const string& str) {
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

// Pfm load
static inline float* load_pfm(
    const char* filename, int* w, int* h, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
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
    auto pixels  = unique_ptr<float[]>(new float[nvalues]);
    for (auto j = *h - 1; j >= 0; j--) {
        if (fread(pixels.get() + j * nrow, sizeof(float), nrow, fs) != nrow)
            return nullptr;
    }

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels.get() + i);
            swap(dta[0], dta[3]);
            swap(dta[1], dta[2]);
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
    auto cpixels = unique_ptr<float[]>(new float[req * npixels]);
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
static inline bool save_pfm(
    const char* filename, int w, int h, int nc, const float* pixels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
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

// load pfm image
template <int N>
static inline void load_pfm(const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    delete[] pixels;
}
template <int N>
static inline void save_pfm(
    const string& filename, const image<vec<float, N>>& img) {
    if (!save_pfm(filename.c_str(), img.size().x, img.size().y, N,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load exr image weith tiny exr
static inline const char* get_tinyexr_error(int error) {
    switch (error) {
        case TINYEXR_ERROR_INVALID_MAGIC_NUMBER: return "INVALID_MAGIC_NUMBER";
        case TINYEXR_ERROR_INVALID_EXR_VERSION: return "INVALID_EXR_VERSION";
        case TINYEXR_ERROR_INVALID_ARGUMENT: return "INVALID_ARGUMENT";
        case TINYEXR_ERROR_INVALID_DATA: return "INVALID_DATA";
        case TINYEXR_ERROR_INVALID_FILE: return "INVALID_FILE";
        // case TINYEXR_ERROR_INVALID_PARAMETER: return "INVALID_PARAMETER";
        case TINYEXR_ERROR_CANT_OPEN_FILE: return "CANT_OPEN_FILE";
        case TINYEXR_ERROR_UNSUPPORTED_FORMAT: return "UNSUPPORTED_FORMAT";
        case TINYEXR_ERROR_INVALID_HEADER: return "INVALID_HEADER";
        default: throw io_error("unknown tinyexr error");
    }
}

template <int N>
static inline void load_exr(const string& filename, image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    auto width = 0, height = 0;
    auto pixels = (float*)nullptr;
    if (auto error = LoadEXR(
            &pixels, &width, &height, filename.c_str(), nullptr);
        error < 0) {
        throw io_error("error loading image " + filename + "("s +
                       get_tinyexr_error(error) + ")"s);
    }
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void save_exr(
    const string& filename, const image<vec<float, N>>& img) {
    // TODO
    if (N != 4) throw runtime_error("bad number of channels");
    if (!SaveEXR((float*)img.data(), img.size().x, img.size().y, N,
            filename.c_str())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
static inline void load_stb(const string& filename, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void load_stb(const string& filename, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, N);
    if (!pixels) {
        throw io_error("error loading image " + filename);
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

// save an image with stbi
template <int N>
static inline void save_png(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_png(filename.c_str(), img.size().x, img.size().y, N,
            img.data(), img.size().x * 4)) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_jpg(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_jpg(
            filename.c_str(), img.size().x, img.size().y, 4, img.data(), 75)) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_tga(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_tga(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_bmp(
    const string& filename, const image<vec<byte, N>>& img) {
    if (!stbi_write_bmp(
            filename.c_str(), img.size().x, img.size().y, 4, img.data())) {
        throw io_error("error saving image " + filename);
    }
}
template <int N>
static inline void save_hdr(
    const string& filename, const image<vec<float, N>>& img) {
    if (!stbi_write_hdr(filename.c_str(), img.size().x, img.size().y, 4,
            (float*)img.data())) {
        throw io_error("error saving image " + filename);
    }
}

// load an image using stbi library
template <int N>
static inline void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<byte, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_load_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image");
    }
    img = image{{width, height}, (const vec<byte, N>*)pixels};
    free(pixels);
}
template <int N>
static inline void load_stb_image_from_memory(
    const byte* data, int data_size, image<vec<float, N>>& img) {
    auto width = 0, height = 0, ncomp = 0;
    auto pixels = stbi_loadf_from_memory(
        data, data_size, &width, &height, &ncomp, 4);
    if (!pixels) {
        throw io_error("error loading in-memory image {}");
    }
    img = image{{width, height}, (const vec<float, N>*)pixels};
    free(pixels);
}

template <int N>
static inline void load_image_preset(
    const string& filename, image<vec<float, N>>& img) {
    auto [type, nfilename] = get_preset_type(filename);
    if constexpr (N == 4) {
        img.resize({1024, 1024});
        if (type == "images2") img.resize({2048, 1024});
        make_preset(img, type);
    } else {
        auto img4 = image<vec<float, 4>>({1024, 1024});
        if (type == "images2") img4.resize({2048, 1024});
        make_preset(img4, type);
        img.resize(img4.size());
        convert_channels(img, img4);
    }
}
template <int N>
static inline void load_image_preset(
    const string& filename, image<vec<byte, N>>& img) {
    auto imgf = image<vec<float, N>>{};
    load_image_preset(filename, imgf);
    img.resize(imgf.size());
    linear_to_srgb(img, imgf);
}

// Forward declarations
template <int N>
static inline void load_image_impl(
    const string& filename, image<vec<byte, N>>& img);
template <int N>
static inline void save_image_impl(
    const string& filename, const image<vec<byte, N>>& img);

// Loads an hdr image.
template <int N>
static inline void load_image_impl(
    const string& filename, image<vec<float, N>>& img) {
    if (is_preset_filename(filename)) {
        return load_image_preset(filename, img);
    }
    auto ext = get_extension(filename);
    if (ext == "exr" || ext == "EXR") {
        load_exr(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        load_pfm(filename, img);
    } else if (ext == "hdr" || ext == "HDR") {
        load_stb(filename, img);
    } else if (!is_hdr_filename(filename)) {
        auto img8 = image<vec<byte, N>>{};
        load_image_impl(filename, img8);
        srgb_to_linear(img, img8);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an hdr image.
template <int N>
static inline void save_image_impl(
    const string& filename, const image<vec<float, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "hdr" || ext == "HDR") {
        save_hdr(filename, img);
    } else if (ext == "pfm" || ext == "PFM") {
        save_pfm(filename, img);
    } else if (ext == "exr" || ext == "EXR") {
        save_exr(filename, img);
    } else if (!is_hdr_filename(filename)) {
        auto img8 = image<vec<byte, N>>{img.size()};
        linear_to_srgb(img8, img);
        save_image_impl(filename, img8);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads an hdr image.
template <int N>
static inline void load_image_impl(
    const string& filename, image<vec<byte, N>>& img) {
    if (is_preset_filename(filename)) {
        return load_image_preset(filename, img);
    }
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        load_stb(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        load_stb(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        load_stb(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        load_stb(filename, img);
    } else if (is_hdr_filename(filename)) {
        auto imgf = image<vec<float, N>>{};
        load_image_impl(filename, imgf);
        linear_to_srgb(img, imgf);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Saves an ldr image.
template <int N>
static inline void save_image_impl(
    const string& filename, const image<vec<byte, N>>& img) {
    auto ext = get_extension(filename);
    if (ext == "png" || ext == "PNG") {
        save_png(filename, img);
    } else if (ext == "jpg" || ext == "JPG") {
        save_jpg(filename, img);
    } else if (ext == "tga" || ext == "TGA") {
        save_tga(filename, img);
    } else if (ext == "bmp" || ext == "BMP") {
        save_bmp(filename, img);
    } else if (is_hdr_filename(filename)) {
        auto imgf = image<vec<float, N>>{img.size()};
        srgb_to_linear(imgf, img);
        save_image_impl(filename, imgf);
    } else {
        throw io_error("unsupported image format " + ext);
    }
}

// Loads/saves a 1-4 channels float image in linear color space.
void load_image(const string& filename, image<float>& img) {
    load_image_impl(filename, (image<vec1f>&)img);
}
void load_image(const string& filename, image<vec2f>& img) {
    load_image_impl(filename, img);
}
void load_image(const string& filename, image<vec3f>& img) {
    load_image_impl(filename, img);
}
void load_image(const string& filename, image<vec4f>& img) {
    load_image_impl(filename, img);
}
void save_image(const string& filename, const image<float>& img) {
    save_image_impl(filename, (const image<vec1f>&)img);
}
void save_image(const string& filename, const image<vec2f>& img) {
    save_image_impl(filename, img);
}
void save_image(const string& filename, const image<vec3f>& img) {
    save_image_impl(filename, img);
}
void save_image(const string& filename, const image<vec4f>& img) {
    save_image_impl(filename, img);
}

// Loads/saves a 1-4 byte image in sRGB color space.
void load_image(const string& filename, image<byte>& img) {
    load_image_impl(filename, (image<vec1b>&)img);
}
void load_image(const string& filename, image<vec2b>& img) {
    load_image_impl(filename, img);
}
void load_image(const string& filename, image<vec3b>& img) {
    load_image_impl(filename, img);
}
void load_image(const string& filename, image<vec4b>& img) {
    load_image_impl(filename, img);
}
void save_image(const string& filename, const image<byte>& img) {
    save_image_impl(filename, (const image<vec1b>&)img);
}
void save_image(const string& filename, const image<vec2b>& img) {
    save_image_impl(filename, img);
}
void save_image(const string& filename, const image<vec3b>& img) {
    save_image_impl(filename, img);
}
void save_image(const string& filename, const image<vec4b>& img) {
    save_image_impl(filename, img);
}

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename) {
    return get_extension(filename) == "hdr" ||
           get_extension(filename) == "exr" || get_extension(filename) == "pfm";
}

// Convenience helper for loading HDR or LDR based on filename
void load_image(
    const string& filename, image<vec4f>& hdr, image<vec4b>& ldr) {
    if (is_hdr_filename(filename)) {
        load_image(filename, hdr);
    } else {
        load_image(filename, ldr);
    }
}
void save_image(
    const string& filename, const image<vec4f>& hdr, const image<vec4b>& ldr) {
    if (!hdr.empty()) {
        save_image(filename, hdr);
    } else {
        save_image(filename, ldr);
    }
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
void save_tonemapped(const string& filename,
    const image<vec4f>& hdr, const tonemap_params& params) {
    if (is_hdr_filename(filename)) {
        save_image(filename, hdr);
    } else {
        auto ldr = image<vec4b>{hdr.size()};
        tonemap(ldr, hdr, params);
        save_image(filename, ldr);
    }
}

// Save with a logo embedded
void save_image_with_logo(
    const string& filename, const image<vec4f>& img) {
    auto logo = image<vec4f>{};
    make_logo(logo, "logo-render");
    auto img_copy = img;
    auto offset   = img.size() - logo.size() - 8;
    set_region(img_copy, logo, offset);
    save_image(filename, img_copy);
}
void save_image_with_logo(
    const string& filename, const image<vec4b>& img) {
    auto logo = image<vec4b>{};
    make_logo(logo, "logo-render");
    auto img_copy = img;
    auto offset   = img.size() - logo.size() - 8;
    set_region(img_copy, logo, offset);
    save_image(filename, img_copy);
}

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
void save_tonemapped_with_logo(const string& filename,
    const image<vec4f>& hdr, const tonemap_params& params) {
    if (is_hdr_filename(filename)) {
        save_image_with_logo(filename, hdr);
    } else {
        auto ldr = image<vec4b>{hdr.size()};
        tonemap(ldr, hdr, params);
        save_image_with_logo(filename, ldr);
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace yocto {

namespace impl {

// Volume load
static inline float* load_yvol(
    const char* filename, int* w, int* h, int* d, int* nc, int req) {
    auto fs = fopen(filename, "rb");
    if (!fs) return nullptr;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
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
    auto voxels  = unique_ptr<float[]>(new float[nvalues]);
    if (fread(voxels.get(), sizeof(float), nvalues, fs) != nvalues)
        return nullptr;

    // proper number of channels
    if (!req || *nc == req) return voxels.release();

    // pack into channels
    if (req < 0 || req > 4) {
        return nullptr;
    }
    auto cvoxels = unique_ptr<float[]>(new float[req * nvoxels]);
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
static inline bool save_yvol(
    const char* filename, int w, int h, int d, int nc, const float* voxels) {
    auto fs = fopen(filename, "wb");
    if (!fs) return false;
    auto fs_guard = unique_ptr<FILE, void (*)(FILE*)>{
        fs, [](FILE* f) { fclose(f); }};

    if (fprintf(fs, "YVOL\n") < 0) return false;
    if (fprintf(fs, "%d %d %d %d\n", w, h, d, nc) < 0) return false;
    auto nvalues = (size_t)w * (size_t)h * (size_t)d * (size_t)nc;
    if (fwrite(voxels, sizeof(float), nvalues, fs) != nvalues) return false;

    return true;
}

// Loads volume data from binary format.
void load_volume(const string& filename, volume<float>& vol) {
    auto width = 0, height = 0, depth = 0, ncomp = 0;
    auto voxels = load_yvol(
        filename.c_str(), &width, &height, &depth, &ncomp, 1);
    if (!voxels) {
        throw io_error("error loading volume " + filename);
    }
    vol = volume{{width, height, depth}, (const float*)voxels};
    delete[] voxels;
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume<float>& vol) {
    if (!save_yvol(filename.c_str(), vol.size().x, vol.size().y, vol.size().z,
            1, vol.data())) {
        throw io_error("error saving volume " + filename);
    }
}

}  // namespace impl

// Loads volume data from binary format.
void load_volume(const string& filename, volume<float>& vol) {
    impl::load_volume(filename, vol);
}

// Saves volume data in binary format.
void save_volume(const string& filename, const volume<float>& vol) {
    impl::save_volume(filename, vol);
}

}  // namespace yocto
