//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

// general includes ------------
#include "yapp.h"

#include <map>

#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"

#include "sunsky/ArHosekSkyModel.c"
#include "sunsky/ArHosekSkyModel.h"

ym::frame3f xform(const ym::vec3f& pos, const ym::vec3f& rot) {
    ym::frame3f xf = ym::identity_frame3f;
    xf = ym::rotation_frame3(ym::vec3f(1, 0, 0), rot[0] * ym::pif / 180) * xf;
    xf = ym::rotation_frame3(ym::vec3f(0, 1, 0), rot[1] * ym::pif / 180) * xf;
    xf = ym::rotation_frame3(ym::vec3f(0, 0, 1), rot[2] * ym::pif / 180) * xf;
    xf = ym::translation_frame3(pos) * xf;
    return xf;
}

yapp::shape transform_shape(const yapp::shape& shape_, const ym::frame3f& xf,
                            const ym::vec3f& s = ym::vec3f{1, 1, 1}) {
    auto shape = shape_;
    for (auto& p : shape.pos) p = ym::transform_point(xf, s * p);
    for (auto& n : shape.norm) n = ym::transform_direction(xf, s * n);
    return shape;
}

yapp::shape transform_shape(const yapp::shape& shape_, const ym::vec3f& pos,
                            const ym::vec3f& rot = ym::zero3f,
                            const ym::vec3f& s = ym::vec3f{1, 1, 1}) {
    return transform_shape(shape_, xform(pos, rot), s);
}

yapp::shape xform_shape(const yapp::shape& shape_, const ym::vec3f& pos,
                        const ym::vec3f& rot,
                        const ym::vec3f& s = ym::vec3f{1, 1, 1}) {
    return transform_shape(shape_, xform(pos, rot), s);
}

yapp::shape lookat_shape(const yapp::shape& shape, const ym::vec3f& pos,
                         const ym::vec3f& to = {0, 0, 0},
                         const ym::vec3f& s = ym::vec3f{1, 1, 1}) {
    auto xf = ym::lookat_frame3(pos, to, {0, 1, 0});
    xf[2] = -xf[2];
    xf[0] = -xf[0];
    return transform_shape(shape, xf, s);
}

yapp::shape make_shape(const std::string& name, int matid, int l,
                       yshape::stype stype) {
    ym::vec4f params = {0.75f, 0.75f, 0, 0};
    yapp::shape shape;
    shape.name = name;
    shape.matid = matid;
    yshape::make_stdsurface(stype, l, params, shape.triangles, shape.pos,
                            shape.norm, shape.texcoord);
    return shape;
}

yapp::shape make_floor(const std::string& name, int matid, float s, float p,
                       int l) {
    auto n = (int)round(powf(2, (float)l));
    yapp::shape shape;
    shape.name = name;
    shape.matid = matid;
    yshape::make_uvsurface(
        n, n, shape.triangles, shape.pos, shape.norm, shape.texcoord,
        [p](const ym::vec2f& uv) {
            auto pos = ym::zero3f;
            auto x = 2 * uv[0] - 1;
            auto y = 2 * (1 - uv[1]) - 1;
            if (y >= 0 || !p) {
                pos = {x, 0, y};
            } else {
                pos = {x, std::pow(-y, p), y};
            }
            return pos;
        },
        [](const ym::vec2f& uv) { return ym::vec3f(0, 1, 0); },
        [s](const ym::vec2f& uv) { return uv * s; });
    if (p) {
        yshape::compute_normals(shape.points, shape.lines, shape.triangles,
                                shape.pos, shape.norm);
    }
    return shape;
}

yapp::material make_material(const std::string& name, const ym::vec3f& ke,
                             const ym::vec3f& kd, const ym::vec3f& ks, float n,
                             int ke_txt, int kd_txt, int ks_txt) {
    yapp::material mat;
    mat.name = name;
    mat.ke = ke;
    mat.kd = kd;
    mat.ks = ks;
    mat.rs = sqrtf(2 / (n + 2));
    mat.ke_txt = ke_txt;
    mat.kd_txt = kd_txt;
    mat.ks_txt = ks_txt;
    return mat;
}

yapp::material make_emission(const std::string& name, const ym::vec3f& ke,
                             int txt = -1) {
    return make_material(name, ke, ym::zero3f, ym::zero3f, 0, txt, -1, -1);
}

yapp::material make_diffuse(const std::string& name, const ym::vec3f& kd,
                            int txt = -1) {
    return make_material(name, ym::zero3f, kd, ym::zero3f, 0, -1, txt, -1);
}

yapp::material make_plastic(const std::string& name, const ym::vec3f& kd,
                            float n, int txt = -1) {
    return make_material(name, ym::zero3f, kd, {0.04f, 0.04f, 0.04f}, n, -1,
                         txt, -1);
}

yapp::material make_metal(const std::string& name, const ym::vec3f& kd, float n,
                          int txt = -1) {
    return make_material(name, ym::zero3f, ym::zero3f, kd, n, 1, 1, txt);
}

yapp::camera make_camera(const std::string& name, const ym::vec3f& from,
                         const ym::vec3f& to, float h, float a) {
    yapp::camera cam;
    cam.name = name;
    cam.frame = ym::lookat_frame3(from, to, {0, 1, 0});
    cam.aperture = a;
    cam.focus = ym::dist(from, to);
    cam.yfov = 2 * std::atan(h / 2);
    cam.aspect = 16.0f / 9.0f;
    return cam;
}

yapp::environment make_env(const std::string& name, int matid,
                           const ym::vec3f& from, const ym::vec3f& to) {
    yapp::environment env;
    env.name = name;
    env.matid = matid;
    env.frame = ym::lookat_frame3(from, to, {0, 1, 0});
    return env;
}

yapp::shape make_points(const std::string& name, int matid, int num) {
    yapp::shape shape;
    shape.name = name;
    shape.matid = matid;

    ym::rng_pcg32 rn;
    yshape::make_points(
        num, shape.points, shape.pos, shape.norm, shape.texcoord, shape.radius,
        [&rn](float u) {
            return ym::vec3f(ym::rng_nextf(rn), ym::rng_nextf(rn),
                             ym::rng_nextf(rn));
        },
        [](float u) {
            return ym::vec3f{0, 0, 1};
        },
        [](float u) {
            return ym::vec2f{u, 0};
        },
        [](float u) { return 0.0025f; });
    return shape;
}

yapp::shape make_lines(const std::string& name, int matid, int num, int n,
                       float r, float c, float s) {
    yapp::shape shape;
    shape.name = name;
    shape.matid = matid;

    ym::rng_pcg32 rn;
    std::vector<ym::vec3f> base(num + 1), dir(num + 1);
    std::vector<float> ln(num + 1);
    for (auto i = 0; i <= num; i++) {
        auto z = -1 + 2 * ym::rng_nextf(rn);
        auto r = std::sqrt(ym::clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * ym::pif * ym::rng_nextf(rn);
        base[i] = ym::vec3f{r * cosf(phi), r * sinf(phi), z};
        dir[i] = base[i];
        ln[i] = 0.15f + 0.15f * ym::rng_nextf(rn);
    }

    yshape::make_lines(
        n, num, shape.lines, shape.pos, shape.norm, shape.texcoord,
        shape.radius,
        [num, base, dir, ln, r, s, c, &rn](const ym::vec2f& uv) {
            auto i = ym::clamp((int)(uv[1] * (num + 1)), 0, num);
            auto pos = base[i] * (1 + uv[0] * ln[i]);
            if (r) {
                pos += ym::vec3f{r * (0.5f - ym::rng_nextf(rn)),
                                 r * (0.5f - ym::rng_nextf(rn)),
                                 r * (0.5f - ym::rng_nextf(rn))};
            }
            if (s && uv[0]) {
                ym::frame3f rotation =
                    ym::rotation_frame3(ym::vec3f(0, 1, 0), s * uv[0] * uv[0]);
                pos = ym::transform_point(rotation, pos);
            }
            auto nc = 128;
            if (c && i > nc) {
                int cc = 0;
                float md = HUGE_VALF;
                for (int k = 0; k < nc; k++) {
                    float d = ym::dist(base[i], base[k]);
                    if (d < md) {
                        md = d;
                        cc = k;
                    }
                }
                ym::vec3f cpos = base[cc] * (1 + uv[0] * ln[cc]);
                pos =
                    pos * (1 - c * uv[0] * uv[0]) + cpos * (c * uv[0] * uv[0]);
            }
            return pos;
        },
        [](const ym::vec2f& uv) {
            return ym::vec3f{0, 0, 1};
        },
        [](const ym::vec2f& uv) { return uv; },
        [](const ym::vec2f& uv) { return 0.001f + 0.001f * (1 - uv[0]); });

    yshape::compute_normals(shape.points, shape.lines, shape.triangles,
                            shape.pos, shape.norm);

    return shape;
}

std::vector<yapp::shape> make_random_shapes(int nshapes, int l,
                                            bool use_xform) {
    std::vector<yapp::shape> shapes(nshapes);
    shapes[0] = transform_shape(make_floor("floor", 0, 6, 4, 6), {0, 0, -4},
                                ym::zero3f, {6, 6, 6});
    ym::vec3f pos[1024];
    float radius[1024];
    int levels[1024];

    ym::rng_pcg32 rn;
    for (auto i = 1; i < nshapes; i++) {
        auto done = false;
        while (!done) {
            auto x = -2 + 4 * ym::rng_nextf(rn);
            auto z = 1 - 3 * ym::rng_nextf(rn);
            radius[i] = 0.15f + ((1 - z) / 3) * ((1 - z) / 3) * 0.5f;
            pos[i] = ym::vec3f{x, radius[i], z};
            levels[i] = (int)round(log2f(powf(2, (float)l) * radius[i] / 0.5f));
            done = true;
            for (int j = 1; j < i && done; j++) {
                if (ym::dist(pos[i], pos[j]) < radius[i] + radius[j])
                    done = false;
            }
        }
    }

    for (auto i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        yshape::stype stypes[3] = {
            yshape::stype::uvspherecube, yshape::stype::uvspherizedcube,
            yshape::stype::uvflipcapsphere,
        };
        auto stype = stypes[(int)(ym::rng_nextf(rn) * 3)];
        if (stype == yshape::stype::uvflipcapsphere) levels[i]++;
        shapes[i] = make_shape(name, i, levels[i], stype);
        if (use_xform) {
            shapes[i] = xform_shape(shapes[i], pos[i], ym::zero3f,
                                    {radius[i], radius[i], radius[i]});
        } else {
            shapes[i] = transform_shape(shapes[i], pos[i], ym::zero3f,
                                        {radius[i], radius[i], radius[i]});
        }
    }

    return shapes;
}

std::vector<yapp::texture> make_random_textures() {
    const std::string txts[5] = {"grid.png", "checker.png", "rchecker.png",
                                 "colored.png", "rcolored.png"};
    std::vector<yapp::texture> textures;
    for (auto txt : txts) {
        textures.push_back({});
        textures.back().path = txt;
    }
    return textures;
}

std::vector<yapp::material> make_random_materials(int nshapes) {
    std::vector<yapp::material> materials(nshapes);
    materials[0] = make_diffuse("floor", {1, 1, 1}, 0);

    ym::rng_pcg32 rn;
    for (auto i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        auto txt = -1;
        if (ym::rng_nextf(rn) < 0.5f) {
            txt = (int)(ym::rng_nextf(rn) * 6) - 1;
        }
        auto c = (txt >= 0) ? ym::vec3f{1, 1, 1}
                            : ym::vec3f{0.2f + 0.3f * ym::rng_nextf(rn),
                                        0.2f + 0.3f * ym::rng_nextf(rn),
                                        0.2f + 0.3f * ym::rng_nextf(rn)};
        auto rs = 0.01f + 0.25f * ym::rng_nextf(rn);
        auto ns = 2 / (rs * rs) - 2;
        auto mt = (int)(ym::rng_nextf(rn) * 4);
        if (mt == 0) {
            materials[i] = make_diffuse(name, c, txt);
        } else if (mt == 1) {
            materials[i] = make_metal(name, c, ns, txt);
        } else {
            materials[i] = make_plastic(name, c, ns, txt);
        }
    }

    return materials;
}

std::vector<yapp::shape> make_random_rigid_shapes(int nshapes, int l) {
    std::vector<yapp::shape> shapes(nshapes);
    shapes[0] = xform_shape(make_shape("floor", 0, 2, yshape::stype::uvcube),
                            {0, -0.5, 0}, ym::zero3f, {6, 0.5, 6});
    ym::vec3f pos[1024];
    float radius[1024];
    int levels[1024];

    ym::rng_pcg32 rn;
    for (int i = 1; i < nshapes; i++) {
        bool done = false;
        while (!done) {
            radius[i] = 0.1f + 0.4f * ym::rng_nextf(rn);
            pos[i] =
                ym::vec3f{-2 + 4 * ym::rng_nextf(rn), 1 + 4 * ym::rng_nextf(rn),
                          -2 + 4 * ym::rng_nextf(rn)};
            levels[i] = (int)round(log2f(powf(2, (float)l) * radius[i] / 0.5f));
            done = true;
            for (int j = 1; j < i && done; j++) {
                if (ym::dist(pos[i], pos[j]) < radius[i] + radius[j])
                    done = false;
            }
        }
    }

    for (int i = 1; i < nshapes; i++) {
        auto name = "obj" + std::to_string(i);
        yshape::stype stypes[2] = {yshape::stype::uvspherecube,
                                   yshape::stype::uvcube};
        auto stype = stypes[(int)(ym::rng_nextf(rn) * 2)];
        shapes[i] = make_shape(name, i, levels[i], stype);
        shapes[i] = xform_shape(shapes[i], pos[i], ym::zero3f,
                                {radius[i], radius[i], radius[i]});
    }

    return shapes;
}

yapp::scene make_scene(const std::vector<yapp::camera>& cameras,
                       const std::vector<yapp::shape>& shapes0,
                       const std::vector<yapp::shape>& shapes1,
                       const std::vector<yapp::material>& materials0,
                       const std::vector<yapp::material>& materials1,
                       const std::vector<yapp::texture>& textures,
                       const std::vector<yapp::environment>& envs = {}) {
    yapp::scene scene;
    scene.cameras = cameras;
    scene.shapes = shapes0;
    for (auto&& s : shapes1) scene.shapes.push_back(s);
    scene.materials = materials0;
    for (auto&& m : materials1) scene.materials.push_back(m);
    scene.textures = textures;
    scene.environments = envs;
    return scene;
}

using ubyte = unsigned char;
struct rgba {
    ubyte r, g, b, a;
};

std::vector<rgba> make_grid(int s) {
    std::vector<rgba> pixels(s * s);
    int g = 64;
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            if (i % g == 0 || i % g == g - 1 || j % g == 0 || j % g == g - 1)
                pixels[j * s + i] = rgba{90, 90, 90, 255};
            else
                pixels[j * s + i] = rgba{128, 128, 128, 255};
        }
    }
    return pixels;
}

std::vector<rgba> make_checker(int s) {
    std::vector<rgba> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            if ((i / 64 + j / 64) % 2)
                pixels[j * s + i] = rgba{90, 90, 90, 255};
            else
                pixels[j * s + i] = rgba{128, 128, 128, 255};
        }
    }
    return pixels;
}

// http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
rgba hsv_to_rgb(ubyte h, ubyte s, ubyte v) {
    rgba rgb = {0, 0, 0, 255};
    ubyte region, remainder, p, q, t;

    if (s == 0) {
        rgb.r = v;
        rgb.g = v;
        rgb.b = v;
        return rgb;
    }

    region = h / 43;
    remainder = (h - (region * 43)) * 6;

    p = (v * (255 - s)) >> 8;
    q = (v * (255 - ((s * remainder) >> 8))) >> 8;
    t = (v * (255 - ((s * (255 - remainder)) >> 8))) >> 8;

    switch (region) {
        case 0:
            rgb.r = v;
            rgb.g = t;
            rgb.b = p;
            break;
        case 1:
            rgb.r = q;
            rgb.g = v;
            rgb.b = p;
            break;
        case 2:
            rgb.r = p;
            rgb.g = v;
            rgb.b = t;
            break;
        case 3:
            rgb.r = p;
            rgb.g = q;
            rgb.b = v;
            break;
        case 4:
            rgb.r = t;
            rgb.g = p;
            rgb.b = v;
            break;
        default:
            rgb.r = v;
            rgb.g = p;
            rgb.b = q;
            break;
    }

    return rgb;
}

std::vector<rgba> make_rcolored(int s) {
    std::vector<rgba> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            ubyte ph = 32 * (i / (s / 8));
            ubyte pv = 128;
            ubyte ps = 64 + 16 * (7 - j / (s / 8));
            if (i % 32 && j % 32) {
                if ((i / 64 + j / 64) % 2)
                    pv += 16;
                else
                    pv -= 16;
                if ((i / 16 + j / 16) % 2)
                    pv += 4;
                else
                    pv -= 4;
                if ((i / 4 + j / 4) % 2)
                    pv += 1;
                else
                    pv -= 1;
            } else {
                pv = 196;
                ps = 32;
            }
            pixels[j * s + i] = hsv_to_rgb(ph, ps, pv);
        }
    }
    return pixels;
}

std::vector<rgba> make_gammaramp(int s) {
    std::vector<rgba> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            auto u = j / float(s - 1);
            if (i < s / 3) u = std::pow(u, 2.2f);
            if (i > (s * 2) / 3) u = std::pow(u, 1 / 2.2f);
            auto c = (unsigned char)(u * 255);
            pixels[j * s + i] = {c, c, c, 255};
        }
    }
    return pixels;
}

std::vector<ym::vec4f> make_gammarampf(int s) {
    std::vector<ym::vec4f> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            auto u = j / float(s - 1);
            if (i < s / 3) u = std::pow(u, 2.2f);
            if (i > (s * 2) / 3) u = std::pow(u, 1 / 2.2f);
            pixels[j * s + i] = {u, u, u, 1};
        }
    }
    return pixels;
}

std::vector<rgba> make_colored(int s) {
    std::vector<rgba> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            ubyte ph = 32 * (i / (s / 8));
            ubyte pv = 128;
            ubyte ps = 64 + 16 * (7 - j / (s / 8));
            if (i % 32 && j % 32) {
                if ((i / 64 + j / 64) % 2)
                    pv += 16;
                else
                    pv -= 16;
            } else {
                pv = 196;
                ps = 32;
            }
            pixels[j * s + i] = hsv_to_rgb(ph, ps, pv);
        }
    }
    return pixels;
}

std::vector<rgba> make_rchecker(int s) {
    std::vector<rgba> pixels(s * s);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            ubyte pv = 128;
            if (i % 32 && j % 32) {
                if ((i / 64 + j / 64) % 2)
                    pv += 16;
                else
                    pv -= 16;
                if ((i / 16 + j / 16) % 2)
                    pv += 4;
                else
                    pv -= 4;
                if ((i / 4 + j / 4) % 2)
                    pv += 1;
                else
                    pv -= 1;
            } else {
                pv = 196;
            }
            pixels[j * s + i] = rgba{pv, pv, pv, 255};
        }
    }
    return pixels;
}

#define sqr(x) ((x) * (x))

std::vector<ym::vec4f> make_sunsky_hdr(int w, int h, float sun_theta,
                                       float turbidity, ym::vec3f ground,
                                       float scale, bool include_ground) {
    std::vector<ym::vec4f> rgba(w * h);
    ArHosekSkyModelState* skymodel_state[3] = {
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground[0], sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground[0], sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground[0], sun_theta),
    };
    auto sun_phi = ym::pif;
    auto sun_w = ym::vec3f{cosf(sun_phi) * sinf(sun_theta),
                           sinf(sun_phi) * sinf(sun_theta), cosf(sun_theta)};
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            auto theta = ym::pif * (j + 0.5f) / h;
            auto phi = 2 * ym::pif * (i + 0.5f) / w;
            if (include_ground)
                theta = ym::clamp(theta, 0.0f, ym::pif / 2 - 0.001f);
            auto pw = ym::vec3f{cosf(phi) * sinf(theta),
                                sinf(phi) * sinf(theta), cosf(theta)};
            auto gamma =
                std::acos(ym::clamp(ym::dot(sun_w, pw), (float)-1, (float)1));
            auto sky = ym::vec3f{(float)(arhosek_tristim_skymodel_radiance(
                                     skymodel_state[0], theta, gamma, 0)),
                                 (float)(arhosek_tristim_skymodel_radiance(
                                     skymodel_state[1], theta, gamma, 1)),
                                 (float)(arhosek_tristim_skymodel_radiance(
                                     skymodel_state[2], theta, gamma, 2))};
            rgba[j * w + i] = {scale * sky[0], scale * sky[1], scale * sky[2],
                               1};
        }
    }
    arhosekskymodelstate_free(skymodel_state[0]);
    arhosekskymodelstate_free(skymodel_state[1]);
    arhosekskymodelstate_free(skymodel_state[2]);
    return rgba;
}

void save_image(const std::string& filename, const std::string& dirname,
                const rgba* pixels, int s) {
    std::string path = std::string(dirname) + "/" + std::string(filename);
    stbi_write_png(path.c_str(), s, s, 4, pixels, s * 4);
}

void save_image_hdr(const std::string& filename, const std::string& dirname,
                    const ym::vec4f* pixels, int w, int h) {
    std::string path = std::string(dirname) + "/" + std::string(filename);
    stbi_write_hdr(path.c_str(), w, h, 4, (float*)pixels);
}

void save_scene(const std::string& filename, const std::string& dirname,
                const yapp::scene& scene) {
    std::string errmsg;
    yapp::save_obj_scene(dirname + "/" + filename, scene, errmsg);
    yapp::save_gltf_scene(
        dirname + "/" + ycmd::get_basename(filename) + ".gltf", scene, errmsg);
}

yapp::texture make_texture(const std::string& path) {
    auto txt = yapp::texture();
    txt.path = path;
    return txt;
}

yapp::shape make_point(const std::string& name, int matid,
                       const ym::vec3f& pos = {0, 0, 0},
                       float radius = 0.001f) {
    yapp::shape shape;
    shape.name = name;
    shape.matid = matid;
    shape.points.push_back(0);
    shape.pos.push_back(pos);
    shape.norm.push_back({0, 0, 1});
    shape.radius.push_back(radius);
    return shape;
}

std::vector<yapp::camera> make_simple_cameras() {
    return {make_camera("cam", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0),
            make_camera("cam_dof", {0, 1.5f, 5}, {0, 0.5, 0}, 0.5f, 0.1f)};
}

std::vector<yapp::shape> make_simple_lightshapes(int matoffset,
                                                 bool arealights) {
    if (!arealights) {
        return {
            transform_shape(make_point("light01", matoffset + 0), {0.7f, 4, 3}),
            transform_shape(make_point("light02", matoffset + 1),
                            {-0.7f, 4, 3})};
    } else {
        return {lookat_shape(make_shape("light01", matoffset + 0, 0,
                                        yshape::stype::uvquad),
                             {2, 2, 4}, {0, 1, 0}),
                lookat_shape(make_shape("light02", matoffset + 1, 0,
                                        yshape::stype::uvquad),
                             {-2, 2, 4}, {0, 1, 0})};
    }
}

std::vector<yapp::material> make_simple_lightmaterials(bool arealights) {
    if (!arealights) {
        return {
            make_emission("light01", {100, 100, 100}),
            make_emission("light02", {100, 100, 100}),
        };
    } else {
        return {
            make_emission("light01", {40, 40, 40}),
            make_emission("light02", {40, 40, 40}),
        };
    }
}

yapp::scene make_simple_scene(bool textured, bool arealights) {
    std::vector<yapp::shape> shapes = {
        transform_shape(make_floor("floor", 0, 6, 4, 6), {0, 0, -4}, ym::zero3f,
                        {6, 6, 6}),
        transform_shape(
            make_shape("obj01", 1, 5, yshape::stype::uvflipcapsphere),
            {-1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}),
        transform_shape(
            make_shape("obj02", 2, 4, yshape::stype::uvspherizedcube),
            {0, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}),
        transform_shape(make_shape("obj03", 2, 4, yshape::stype::uvspherecube),
                        {1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f})};
    std::vector<yapp::material> materials;
    std::vector<yapp::texture> textures;
    if (!textured) {
        materials = std::vector<yapp::material>{
            make_diffuse("floor", {0.2f, 0.2f, 0.2f}, -1),
            make_plastic("obj01", {0.5f, 0.2f, 0.2f}, 50, -1),
            make_plastic("obj02", {0.2f, 0.5f, 0.2f}, 100, -1),
            make_plastic("obj03", {0.2f, 0.2f, 0.5f}, 500, -1)};
    } else {
        materials = std::vector<yapp::material>{
            make_diffuse("floor", {1, 1, 1}, 0),
            make_plastic("obj01", {1, 1, 1}, 50, 1),
            make_plastic("obj02", {1, 1, 1}, 100, 2),
            make_plastic("obj03", {1, 1, 1}, 500, 3)};
        textures = std::vector<yapp::texture>{
            make_texture("grid.png"), make_texture("rcolored.png"),
            make_texture("checker.png"), make_texture("colored.png"),
        };
    }
    return make_scene(
        make_simple_cameras(), shapes,
        make_simple_lightshapes((int)materials.size(), arealights), materials,
        make_simple_lightmaterials(arealights), textures);
}

yapp::scene make_pointslines_scene(bool lines, bool arealights) {
    std::vector<yapp::shape> shapes;
    std::vector<yapp::material> materials;
    std::vector<yapp::texture> textures;
    shapes.push_back(transform_shape(make_floor("floor", 0, 6, 4, 6),
                                     {0, 0, -4}, ym::zero3f, {6, 6, 6}));
    materials =
        std::vector<yapp::material>{make_diffuse("floor", {0.2f, 0.2f, 0.2f}),
                                    make_diffuse("obj", {0.2f, 0.2f, 0.2f}),
                                    make_diffuse("points", {0.2f, 0.2f, 0.2f}),
                                    make_diffuse("lines", {0.2f, 0.2f, 0.2f})};
    if (!lines) {
        shapes.push_back(
            transform_shape(make_points("points01", 2, 64 * 64 * 16),
                            {0, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
    } else {
        shapes.push_back(
            transform_shape(make_shape("obj01", 1, 6, yshape::stype::uvsphere),
                            {1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
        shapes.push_back(transform_shape(
            make_lines("lines01", 3, 64 * 64 * 16, 4, 0.1f, 0, 0),
            {1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
        shapes.push_back(
            transform_shape(make_shape("obj02", 1, 6, yshape::stype::uvsphere),
                            {0, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
        shapes.push_back(transform_shape(
            make_lines("lines02", 3, 64 * 64 * 16, 4, 0, 0.75f, 0),
            {0, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
        shapes.push_back(
            transform_shape(make_shape("obj03", 1, 6, yshape::stype::uvsphere),
                            {-1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
        shapes.push_back(transform_shape(
            make_lines("lines03", 3, 64 * 64 * 16, 4, 0, 0, 0.5f),
            {-1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}));
    }

    return make_scene(
        make_simple_cameras(), shapes,
        make_simple_lightshapes((int)materials.size(), arealights), materials,
        make_simple_lightmaterials(arealights), textures);
}

yapp::scene make_random_scene(int nshapes, bool arealights) {
    std::vector<yapp::camera> cameras = {
        make_camera("cam", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {0, 1.5f, 5}, {0, 0.5, 0}, 0.5f, 0.1f)};
    std::vector<yapp::shape> shapes = make_random_shapes(nshapes, 5, false);
    std::vector<yapp::material> materials = make_random_materials(nshapes);
    std::vector<yapp::texture> textures = make_random_textures();
    return make_scene(
        cameras, shapes,
        make_simple_lightshapes((int)materials.size(), arealights), materials,
        make_simple_lightmaterials(arealights), textures);
}

// http://graphics.cs.williams.edu/data
// http://www.graphics.cornell.edu/online/box/data.html
yapp::scene make_cornell_box_scene() {
    std::vector<yapp::camera> cameras = {
        make_camera("cam", {0, 1, 4}, {0, 1, 0}, 0.7f, 0)};
    std::vector<yapp::shape> shapes = {
        transform_shape(make_shape("floor", 0, 0, yshape::stype::uvquad),
                        ym::zero3f, {-90, 0, 0}),
        transform_shape(make_shape("ceiling", 0, 0, yshape::stype::uvquad),
                        {0, 2, 0}, {90, 0, 0}),
        transform_shape(make_shape("back", 0, 0, yshape::stype::uvquad),
                        {0, 1, -1}, ym::zero3f),
        transform_shape(make_shape("back", 2, 0, yshape::stype::uvquad),
                        {+1, 1, 0}, {0, -90, 0}),
        transform_shape(make_shape("back", 1, 0, yshape::stype::uvquad),
                        {-1, 1, 0}, {0, 90, 0}),
        transform_shape(make_shape("tallbox", 0, 0, yshape::stype::uvcube),
                        {-0.33f, 0.6f, -0.29f}, {0, 15, 0}, {0.3f, 0.6f, 0.3f}),
        transform_shape(make_shape("shortbox", 0, 0, yshape::stype::uvcube),
                        {0.33f, 0.3f, 0.33f}, {0, -15, 0}, {0.3f, 0.3f, 0.3f}),
        transform_shape(make_shape("light", 3, 0, yshape::stype::uvquad),
                        {0, 1.999f, 0}, {90, 0, 0}, {0.25f, 0.25f, 0.25f})};
    std::vector<yapp::material> materials = {
        make_diffuse("white", {0.725f, 0.71f, 0.68f}),
        make_diffuse("red", {0.63f, 0.065f, 0.05f}),
        make_diffuse("green", {0.14f, 0.45f, 0.091f}),
        make_emission("light", {17, 12, 4}),
    };
    return make_scene(cameras, shapes, {}, materials, {}, {});
}

yapp::scene make_envmap_scene(bool as_shape, bool use_map) {
    std::vector<yapp::camera> cameras = {
        make_camera("cam", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0.1f)};
    std::vector<yapp::shape> shapes = {
        transform_shape(make_floor("floor", 0, 6, 4, 6), {0, 0, -4}, ym::zero3f,
                        {6, 6, 6}),
        transform_shape(
            make_shape("obj01", 1, 5, yshape::stype::uvflipcapsphere),
            {-1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}),
        transform_shape(
            make_shape("obj02", 2, 4, yshape::stype::uvspherizedcube),
            {0, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f}),
        transform_shape(make_shape("obj03", 3, 4, yshape::stype::uvspherecube),
                        {1.25f, 0.5f, 0}, ym::zero3f, {0.5f, 0.5f, 0.5f})};
    std::vector<yapp::material> materials = {
        make_diffuse("floor", {0.2f, 0.2f, 0.2f}),
        make_plastic("obj01", {0.5f, 0.2f, 0.2f}, 50),
        make_plastic("obj02", {0.2f, 0.5f, 0.2f}, 100),
        make_plastic("obj03", {0.2f, 0.2f, 0.5f}, 500),
        make_emission("env", {1, 1, 1}, (use_map) ? 0 : -1)};
    std::vector<yapp::texture> textures;
    std::vector<yapp::environment> environments;
    if (as_shape) {
        shapes.push_back(transform_shape(
            make_shape("env_sphere", 4, 6, yshape::stype::uvflippedsphere),
            {0, 0.5f, 0}, {-90, 0, 0}, {10000, 10000, 10000}));
    } else {
        environments.push_back(
            make_env("env", 4, {0, 0.5f, 0}, {-1.5f, 0.5f, 0}));
    }
    if (use_map) {
        textures.push_back(make_texture("env.hdr"));
    }

    return make_scene(cameras, shapes, {}, materials, {}, textures,
                      environments);
}

yapp::scene make_rigid_scene(int config) {
    std::vector<yapp::camera> cameras = {
        make_camera("cam", {5, 5, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {5, 5, 5}, {0, 0.5f, 0}, 0.5f, 0.1f)};
    std::vector<yapp::shape> shapes;
    std::vector<yapp::material> materials = {
        make_diffuse("floor", {1, 1, 1}, 0),
        make_plastic("obj", {1, 1, 1}, 50, 1)};
    std::vector<yapp::texture> textures = {make_texture("grid.png"),
                                           make_texture("checker.png")};

    if (config == 0 || config == 1) {
        shapes = {
            (config)
                ? xform_shape(make_shape("floor", 0, 2, yshape::stype::uvcube),
                              {0, -2.5, 0}, {30, 0, 0}, {6, 0.5f, 6})
                : xform_shape(make_shape("floor", 0, 4, yshape::stype::uvcube),
                              {0, -0.5f, 0}, {0, 0, 0}, {6, 0.5f, 6}),
            xform_shape(make_shape("obj01", 1, 2, yshape::stype::uvcube),
                        {-1.25f, 0.5f, 0}, {0, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj02", 1, 3, yshape::stype::uvspherecube),
                        {0, 1, 0}, {0, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj03", 1, 2, yshape::stype::uvcube),
                        {1.25f, 1.5f, 0}, {0, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj11", 1, 2, yshape::stype::uvcube),
                        {-1.25f, 0.5f, 1.5f}, {0, 45, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj12", 1, 3, yshape::stype::uvspherecube),
                        {0, 1, 1.5f}, {45, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj13", 1, 2, yshape::stype::uvcube),
                        {1.25f, 1.5f, 1.5f}, {45, 0, 45}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj21", 1, 2, yshape::stype::uvcube),
                        {-1.25f, 0.5f, -1.5f}, {0, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj22", 1, 3, yshape::stype::uvspherecube),
                        {0, 1, -1.5f}, {22.5, 0, 0}, {0.5f, 0.5f, 0.5f}),
            xform_shape(make_shape("obj23", 1, 2, yshape::stype::uvcube),
                        {1.25f, 1.5f, -1.5f}, {22.5f, 0, 22.5f},
                        {0.5f, 0.5f, 0.5f})};
    } else if (config == 2) {
        shapes = make_random_rigid_shapes(128, 1);
    } else {
        assert(false);
    }

    shapes.push_back(transform_shape(make_point("light01", 2), {0.7f, 4, 3},
                                     {0, 0, 0}, {1, 1, 1}));
    shapes.push_back(transform_shape(make_point("light02", 3), {-0.7f, 4, 3},
                                     {0, 0, 0}, {1, 1, 1}));
    materials.push_back(make_emission("light01", {100, 100, 100}));
    materials.push_back(make_emission("light02", {100, 100, 100}));

    return make_scene(cameras, shapes, {}, materials, {}, textures);
}

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ycmd::make_parser(argc, argv, "make tests");
    auto dirname = ycmd::parse_arg<std::string>(parser, "dirname",
                                                "directory name", ".", true);
    ycmd::check_parser(parser);

// make directories
#ifndef _MSC_VER
    auto cmd = "mkdir -p " + dirname;
#else
    auto cmd = "mkdir " + dirname;
#endif
    system(cmd.c_str());

    // simple scene ------------------------------
    printf("generating simple scenes ...\n");
    save_scene("basic_pointlight.obj", dirname,
               make_simple_scene(false, false));
    save_scene("simple_pointlight.obj", dirname,
               make_simple_scene(true, false));
    save_scene("simple_arealight.obj", dirname, make_simple_scene(true, true));

    // point and lines scene ------------------------------
    printf("generating points and lines scenes ...\n");
    save_scene("points_pointlight.obj", dirname,
               make_pointslines_scene(false, false));
    save_scene("points_arealight.obj", dirname,
               make_pointslines_scene(false, true));
    save_scene("lines_pointlight.obj", dirname,
               make_pointslines_scene(true, false));
    save_scene("lines_arealight.obj", dirname,
               make_pointslines_scene(true, true));

    // random obj scene --------------------------
    printf("generating random shapes scenes ...\n");
    save_scene("random_pointlight.obj", dirname, make_random_scene(32, false));
    save_scene("random_arealight.obj", dirname, make_random_scene(32, true));

    // env scene ------------------------------
    printf("generating envmaps scenes ...\n");
    save_scene("env_shape_const.obj", dirname, make_envmap_scene(true, false));
    save_scene("env_shape_map.obj", dirname, make_envmap_scene(true, true));
    save_scene("env_inf_const.obj", dirname, make_envmap_scene(false, false));
    save_scene("env_inf_map.obj", dirname, make_envmap_scene(false, true));

    // cornell box ------------------------------
    printf("generating cornell box scenes ...\n");

    // save scene
    save_scene("cornell_box.obj", dirname, make_cornell_box_scene());

    // rigid body scenes ------------------------
    printf("generating rigid body scenes ...\n");
    save_scene("rigid_01.obj", dirname, make_rigid_scene(0));
    save_scene("rigid_02.obj", dirname, make_rigid_scene(1));
    // save_scene("rigid_03.obj", dirname, make_rigid_scene(2));

    // textures ---------------------------------
    printf("generating simple textures ...\n");
    save_image("grid.png", dirname, make_grid(512).data(), 512);
    save_image("checker.png", dirname, make_checker(512).data(), 512);
    save_image("rchecker.png", dirname, make_rchecker(512).data(), 512);
    save_image("colored.png", dirname, make_colored(512).data(), 512);
    save_image("rcolored.png", dirname, make_rcolored(512).data(), 512);
    save_image("gamma.png", dirname, make_gammaramp(512).data(), 512);
    save_image_hdr("gamma.hdr", dirname, make_gammarampf(512).data(), 512, 512);
    printf("generating envmaps textures ...\n");
    save_image_hdr("env.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8f, 8,
                                   ym::vec3f{0.2f, 0.2f, 0.2f}, 1 / powf(2, 6),
                                   true)
                       .data(),
                   1024, 512);
    save_image_hdr("env01.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8f, 8,
                                   ym::vec3f{0.2f, 0.2f, 0.2f}, 1 / powf(2, 6),
                                   true)
                       .data(),
                   1024, 512);
}
