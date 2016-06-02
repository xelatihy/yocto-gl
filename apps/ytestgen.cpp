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
#define YGL_DECLARATION

#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#include "ext/sunsky/ArHosekSkyModel.c"
#include "ext/sunsky/ArHosekSkyModel.h"

ym_frame3f xform(const ym_vec3f& pos, const ym_vec3f& rot) {
    ym_frame3f xf = ym_identity_frame3f;
    xf = ym_rotation_xform3(ym_vec3f(1, 0, 0), rot.x * ym_pif / 180) * xf;
    xf = ym_rotation_xform3(ym_vec3f(0, 1, 0), rot.y * ym_pif / 180) * xf;
    xf = ym_rotation_xform3(ym_vec3f(0, 0, 1), rot.z * ym_pif / 180) * xf;
    xf = ym_translation_xform3(pos) * xf;
    return xf;
}

yo_shape transform_shape(yo_shape shape, const ym_vec3f& pos,
                         const ym_vec3f& rot, const ym_vec3f& s) {
    ym_frame3f xf = xform(pos, rot);
    for (int i = 0; i < shape.nverts; i++) {
        shape.pos[i] = ym_transform_point(xf, s * shape.pos[i]);
        shape.norm[i] = ym_transform_direction(xf, s * shape.norm[i]);
    }
    return shape;
}

yo_shape transform_shape(yo_shape shape, const ym_vec3f& pos,
                         const ym_vec3f& rot, float s) {
    return transform_shape(shape, pos, rot, {s, s, s});
}

yo_shape xform_shape(yo_shape shape, const ym_vec3f& pos, const ym_vec3f& rot,
                     const ym_vec3f& s) {
    for (int i = 0; i < shape.nverts; i++) shape.pos[i] *= s;
    ym_affine3f xf = (ym_affine3f)xform(pos, rot);
    if (shape.xformed) xf = xf * shape.xform;
    shape.xform = xf;
    shape.xformed = true;
    return shape;
}

yo_shape xform_shape(yo_shape shape, const ym_vec3f& pos, const ym_vec3f& rot,
                     float s) {
    return xform_shape(shape, pos, rot, {s, s, s});
}

yo_shape make_shape(const char* name, const char* matname, int l, int stype) {
    int etype = yo_etype_triangle;
    ym_vec4f params = {0.75f, 0.75f, 0, 0};
    if (stype == ys_stype_points || stype == ys_stype_rnpoints) {
        etype = yo_etype_point;
        params[0] = 0.001f;
        params[0] = 0.0025f;
    }
    yo_shape shape;
    shape.name = name;
    shape.matname = matname;
    shape.etype = etype;
    ys_make_stdshape(stype, l, etype, params, &shape.nelems, &shape.elem,
                     &shape.nverts, &shape.pos, &shape.norm, &shape.texcoord,
                     0);
    return shape;
}

yo_shape make_floor(const char* name, const char* matname, float s, float p,
                    int l) {
    int n = (int)round(powf(2, (float)l));
    yo_shape shape;
    shape.name = name;
    shape.matname = matname;
    shape.etype = yo_etype_triangle;
    ys_get_uvgrid_size(n, n, ys_etype_triangle, &shape.nelems, &shape.nverts);
    shape.elem.resize(shape.nelems * 3);
    shape.texcoord.resize(shape.nverts);
    shape.pos.resize(shape.nverts);
    shape.norm.resize(shape.nverts);
    ys_make_uvgrid(n, n, ys_etype_triangle, shape.elem.data(),
                   shape.texcoord.data());
    for (int i = 0; i < shape.nverts; i++) {
        ym_vec2f uv = shape.texcoord[i];
        float x = 2 * uv[0] - 1;
        float y = 2 * (1 - uv[1]) - 1;
        if (y >= 0 || !p) {
            shape.pos[i] = {x, 0, y};
        } else {
            shape.pos[i] = {x, powf(-y, p), y};
        }

        shape.norm[i] = {0, 1, 0};
        shape.texcoord[i] = uv * s;
    }
    if (p) {
        ys_compute_normals(shape.nelems, shape.elem.data(), shape.etype,
                           shape.nverts, shape.pos.data(), shape.norm.data(),
                           true);
    }
    return shape;
}

yo_material make_material(const char* name, const ym_vec3f& ke,
                          const ym_vec3f& kd, const ym_vec3f& ks, float ns,
                          const char* ke_txt, const char* kd_txt,
                          const char* ks_txt) {
    yo_material mat;
    mat.name = name;
    mat.ke = ke;
    mat.kd = kd;
    mat.ks = ks;
    mat.ns = ns;
    mat.ke_txt = (ke_txt) ? ke_txt : "";
    mat.kd_txt = (kd_txt) ? kd_txt : "";
    mat.ks_txt = (ks_txt) ? ks_txt : "";
    return mat;
}

yo_material make_emission(const char* name, const ym_vec3f& ke,
                          const char* txt) {
    return make_material(name, ke, ym_zero3f, ym_zero3f, 0, txt, 0, 0);
}

yo_material make_diffuse(const char* name, const ym_vec3f& kd,
                         const char* txt) {
    return make_material(name, ym_zero3f, kd, ym_zero3f, 0, 0, txt, 0);
}

yo_material make_plastic(const char* name, const ym_vec3f& kd, float n,
                         const char* txt) {
    return make_material(name, ym_zero3f, kd, {0.04f, 0.04f, 0.04f}, n, 0, txt,
                         0);
}

yo_material make_metal(const char* name, const ym_vec3f& kd, float n,
                       const char* txt) {
    return make_material(name, ym_zero3f, ym_zero3f, kd, n, 0, 0, txt);
}

yo_camera make_camera(const char* name, const ym_vec3f& from,
                      const ym_vec3f& to, float h, float a) {
    yo_camera cam;
    cam.name = name;
    cam.from = from;
    cam.to = to;
    cam.up = ym_vec3f{0, 1, 0};
    cam.aperture = a;
    cam.height = h;
    cam.width = cam.height * 16.0f / 9.0f;
    return cam;
}

yo_env make_env(const char* name, const char* matname, const ym_vec3f& from,
                const ym_vec3f& to) {
    yo_env env;
    env.name = name;
    env.matname = matname;
    env.from = from;
    env.to = to;
    env.up = ym_vec3f{0, 1, 0};
    return env;
}

yo_shape make_lines(const char* name, const char* matname, int num, int n,
                    float r, float c, float s) {
    yo_shape shape;
    shape.name = name;
    shape.matname = matname;
    shape.etype = yo_etype_line;

    ym_rng_pcg32 rn;
    ym_vector<ym_vec3f> base(num + 1), dir(num + 1);
    ym_vector<float> ln(num + 1);
    for (int i = 0; i <= num; i++) {
        float z = -1 + 2 * ym_rng_nextf(&rn);
        float r = sqrtf(ym_clamp(1 - z * z, 0, 1));
        float phi = 2 * ym_pif * ym_rng_nextf(&rn);
        base[i] = ym_vec3f{r * cosf(phi), r * sinf(phi), z};
        dir[i] = base[i];
        ln[i] = 0.15f + 0.15f * ym_rng_nextf(&rn);
    }

    ys_get_uvgrid_size(n, num, shape.etype, &shape.nelems, &shape.nverts);
    shape.elem.resize(shape.nelems * shape.etype);
    shape.texcoord.resize(shape.nverts);
    shape.pos.resize(shape.nverts);
    shape.norm.resize(shape.nverts);
    shape.radius.resize(shape.nverts);
    ys_make_uvgrid(n, num, shape.etype, shape.elem.data(),
                   shape.texcoord.data());

    for (int j = 0; j < shape.nverts; j++) {
        ym_vec2f uv = shape.texcoord[j];

        int i = ym_clamp((int)(uv[1] * (num + 1)), 0, num);
        ym_vec3f pos = base[i] * (1 + uv[0] * ln[i]);
        if (r) {
            pos += ym_vec3f{r * (0.5f - ym_rng_nextf(&rn)),
                            r * (0.5f - ym_rng_nextf(&rn)),
                            r * (0.5f - ym_rng_nextf(&rn))};
        }
        if (s && uv[0]) {
            ym_frame3f rotation =
                ym_rotation_xform3(ym_vec3f(0, 1, 0), s * uv[0] * uv[0]);
            pos = ym_transform_point(rotation, pos);
        }
        int nc = 128;
        if (c && i > nc) {
            int cc = 0;
            float md = HUGE_VALF;
            for (int k = 0; k < nc; k++) {
                float d = ym_dist(base[i], base[k]);
                if (d < md) {
                    md = d;
                    cc = k;
                }
            }
            ym_vec3f cpos = base[cc] * (1 + uv[0] * ln[cc]);
            pos = pos * (1 - c * uv[0] * uv[0]) + cpos * (c * uv[0] * uv[0]);
        }

        shape.pos[j] = pos;
        shape.radius[j] = 0.001f + 0.001f * (1 - shape.texcoord[j][0]);
    }

    ys_compute_normals(shape.nelems, shape.elem.data(), shape.etype,
                       shape.nverts, shape.pos.data(), shape.norm.data(), true);

    return shape;
}

ym_vector<yo_shape> make_random_shapes(int nshapes, int l, bool use_xform) {
    ym_vector<yo_shape> shapes(nshapes);
    shapes[0] = transform_shape(make_floor("floor", "floor", 6, 4, 6),
                                {0, 0, -4}, ym_zero3f, 6);
    ym_vec3f pos[1024];
    float radius[1024];
    int levels[1024];

    ym_rng_pcg32 rn;
    for (int i = 1; i < nshapes; i++) {
        bool done = false;
        while (!done) {
            float x = -2 + 4 * ym_rng_nextf(&rn);
            float z = 1 - 3 * ym_rng_nextf(&rn);
            radius[i] = 0.15f + ((1 - z) / 3) * ((1 - z) / 3) * 0.5f;
            pos[i] = ym_vec3f{x, radius[i], z};
            levels[i] = (int)round(log2f(powf(2, (float)l) * radius[i] / 0.5f));
            done = true;
            for (int j = 1; j < i && done; j++) {
                if (ym_dist(pos[i], pos[j]) < radius[i] + radius[j])
                    done = false;
            }
        }
    }

    for (int i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        int stypes[3] = {
            ys_stype_uvspherecube, ys_stype_uvspherizedcube,
            ys_stype_uvflipcapsphere,
        };
        int stype = stypes[(int)(ym_rng_nextf(&rn) * 3)];
        if (stype == ys_stype_uvflipcapsphere) levels[i]++;
        shapes[i] = make_shape(name, name, levels[i], stype);
        if (use_xform) {
            shapes[i] = xform_shape(shapes[i], pos[i], ym_zero3f, radius[i]);
        } else {
            shapes[i] =
                transform_shape(shapes[i], pos[i], ym_zero3f, radius[i]);
        }
    }

    return shapes;
}

ym_vector<yo_material> make_random_materials(int nshapes) {
    ym_vector<yo_material> materials(nshapes);
    materials[0] = make_diffuse("floor", {1, 1, 1}, "grid.png");

    const char* txts[6] = {0,
                           "grid.png",
                           "checker.png",
                           "rchecker.png",
                           "colored.png",
                           "rcolored.png"};
    ym_rng_pcg32 rn;
    for (int i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        const char* txt = 0;
        if (ym_rng_nextf(&rn) < 0.5f) {
            txt = txts[(int)(ym_rng_nextf(&rn) * 6)];
        }
        ym_vec3f c = (txt) ? ym_vec3f{1, 1, 1}
                           : ym_vec3f{0.2f + 0.3f * ym_rng_nextf(&rn),
                                      0.2f + 0.3f * ym_rng_nextf(&rn),
                                      0.2f + 0.3f * ym_rng_nextf(&rn)};
        float rs = 0.01f + 0.25f * ym_rng_nextf(&rn);
        float ns = 2 / (rs * rs) - 2;
        int mt = (int)(ym_rng_nextf(&rn) * 4);
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

ym_vector<yo_shape> make_random_rigid_shapes(int nshapes, int l) {
    ym_vector<yo_shape> shapes(nshapes);
    shapes[0] = xform_shape(make_shape("floor", "floor", 2, ys_stype_uvcube),
                            {0, -0.5, 0}, ym_zero3f, {6, 0.5, 6});
    ym_vec3f pos[1024];
    float radius[1024];
    int levels[1024];

    ym_rng_pcg32 rn;
    for (int i = 1; i < nshapes; i++) {
        bool done = false;
        while (!done) {
            radius[i] = 0.1f + 0.4f * ym_rng_nextf(&rn);
            pos[i] =
                ym_vec3f{-2 + 4 * ym_rng_nextf(&rn), 1 + 4 * ym_rng_nextf(&rn),
                         -2 + 4 * ym_rng_nextf(&rn)};
            levels[i] = (int)round(log2f(powf(2, (float)l) * radius[i] / 0.5f));
            done = true;
            for (int j = 1; j < i && done; j++) {
                if (ym_dist(pos[i], pos[j]) < radius[i] + radius[j])
                    done = false;
            }
        }
    }

    for (int i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        int stypes[2] = {ys_stype_uvspherecube, ys_stype_uvcube};
        int stype = stypes[(int)(ym_rng_nextf(&rn) * 2)];
        shapes[i] = make_shape(name, "obj", levels[i], stype);
        shapes[i] = xform_shape(shapes[i], pos[i], ym_zero3f, radius[i]);
    }

    return shapes;
}

yo_scene make_scene(const ym_vector<yo_camera>& cameras,
                    const ym_vector<yo_shape>& shapes,
                    const ym_vector<yo_material>& materials,
                    const ym_vector<yo_env>& envs = {}) {
    yo_scene scene;
    scene.cameras = cameras;
    scene.shapes = shapes;
    scene.materials = materials;
    scene.envs = envs;
    return scene;
}

typedef unsigned char ubyte;
typedef struct { ubyte r, g, b, a; } rgba;

ym_vector<rgba> make_grid(int s) {
    ym_vector<rgba> pixels(s * s);
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

ym_vector<rgba> make_checker(int s) {
    ym_vector<rgba> pixels(s * s);
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

ym_vector<rgba> make_rcolored(int s) {
    ym_vector<rgba> pixels(s * s);
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

ym_vector<rgba> make_colored(int s) {
    ym_vector<rgba> pixels(s * s);
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

ym_vector<rgba> make_rchecker(int s) {
    ym_vector<rgba> pixels(s * s);
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

ym_vector<ym_vec4f> make_sunsky_hdr(int w, int h, float sun_theta,
                                    float turbidity, ym_vec3f ground,
                                    float scale, bool include_ground) {
    ym_vector<ym_vec4f> rgba(w * h);
    ArHosekSkyModelState* skymodel_state[3] = {
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
    };
    float sun_phi = ym_pif;
    ym_vec3f sun_w = ym_vec3f{cosf(sun_phi) * sinf(sun_theta),
                              sinf(sun_phi) * sinf(sun_theta), cosf(sun_theta)};
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            float theta = ym_pif * (j + 0.5f) / h;
            float phi = 2 * ym_pif * (i + 0.5f) / w;
            if (include_ground)
                theta = ym_clamp(theta, 0.0f, ym_pif / 2 - 0.001f);
            ym_vec3f pw = ym_vec3f{cosf(phi) * sinf(theta),
                                   sinf(phi) * sinf(theta), cosf(theta)};
            float gamma = acosf(ym_clamp(ym_dot(sun_w, pw), -1, 1));
            ym_vec3f sky = {(float)(arhosek_tristim_skymodel_radiance(
                                skymodel_state[0], theta, gamma, 0)),
                            (float)(arhosek_tristim_skymodel_radiance(
                                skymodel_state[1], theta, gamma, 1)),
                            (float)(arhosek_tristim_skymodel_radiance(
                                skymodel_state[2], theta, gamma, 2))};
            rgba[j * w + i] = {scale * sky.x, scale * sky.y, scale * sky.z, 1};
        }
    }
    arhosekskymodelstate_free(skymodel_state[0]);
    arhosekskymodelstate_free(skymodel_state[1]);
    arhosekskymodelstate_free(skymodel_state[2]);
    return rgba;
}

void save_image(const char* filename, const char* dirname, const rgba* pixels,
                int s) {
    ym_string path = ym_string(dirname) + "/" + ym_string(filename);
    stbi_write_png(path.c_str(), s, s, 4, pixels, s * 4);
}

void save_image_hdr(const char* filename, const char* dirname,
                    const ym_vec4f* pixels, int w, int h) {
    ym_string path = ym_string(dirname) + "/" + ym_string(filename);
    stbi_write_hdr(path.c_str(), w, h, 4, (float*)pixels);
}

void save_scene(const char* filename, const char* dirname,
                const yo_scene scene) {
    ym_string path = ym_string(dirname) + "/" + ym_string(filename);
    yo_save_obj(path.c_str(), &scene, true);
}

int main(int argc, const char** argv) {
    // command line params
    yc_parser* parser = yc_init_parser(argc, argv, "make tests");
    const char* dirname =
        yc_parse_args(parser, "dirname", "directory name", ".", true);
    yc_done_parser(parser);

    // make directories
    ym_string cmd;
#ifndef _MSC_VER
    cmd = ym_string("mkdir -p ") + ym_string(dirname);
#else
    cmd = ym_string("mkdir ") + ym_string(dirname);
#endif
    system(cmd.c_str());

    // floor scene ------------------------------
    printf("generating simple scenes ...\n");
    ym_vector<yo_camera> sh_cams = {
        make_camera("cam", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {0, 1.5f, 5}, {0, 0.5, 0}, 0.5f, 0.1f)};
    ym_vector<yo_shape> sh_shapes = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), {0, 0, -4},
                        ym_zero3f, 6),
        transform_shape(
            make_shape("obj01", "obj01", 5, ys_stype_uvflipcapsphere),
            {-1.25f, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(
            make_shape("obj02", "obj02", 4, ys_stype_uvspherizedcube),
            {0, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(make_shape("obj03", "obj03", 4, ys_stype_uvspherecube),
                        {1.25f, 0.5f, 0}, ym_zero3f, 0.5)};
    ym_vector<yo_material> sh_materials = {
        make_diffuse("floor", {0.2f, 0.2f, 0.2f}, 0),
        make_plastic("obj01", {0.5f, 0.2f, 0.2f}, 50, 0),
        make_plastic("obj02", {0.2f, 0.5f, 0.2f}, 100, 0),
        make_plastic("obj03", {0.2f, 0.2f, 0.5f}, 500, 0)};
    ym_vector<yo_shape> sh_light_shapes = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_points),
                        {0.7f, 4, 3}, ym_zero3f, 1),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_points),
                        {-0.7f, 4, 3}, ym_zero3f, 1)};
    ym_vector<yo_material> sh_light_materials = {
        make_emission("light01", {100, 100, 100}, 0),
        make_emission("light02", {100, 100, 100}, 0),
    };

    // textures
    ym_vector<yo_material> sh_textured_materials = {
        make_diffuse("floor", {1, 1, 1}, "grid.png"),
        make_plastic("obj01", {1, 1, 1}, 50, "rcolored.png"),
        make_plastic("obj02", {1, 1, 1}, 100, "checker.png"),
        make_plastic("obj03", {1, 1, 1}, 500, "colored.png")};

    // area lights
    ym_vector<yo_shape> sh_area_shapes = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                        {4, 2, 4}, {0, 180 + 45, 0}, 2),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                        {-4, 2, 4}, {0, 180 - 45, 0}, 2)};
    ym_vector<yo_material> sh_area_materials = {
        make_emission("light01", {10, 10, 10}, 0),
        make_emission("light02", {10, 10, 10}, 0),
    };

    // colored area lights
    ym_vector<yo_shape> sh_area1_shapes = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                        {4, 2, 4}, {0, 180 + 45, 0}, 2),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                        {-4, 2, 4}, {0, 180 - 45, 0}, 2)};
    ym_vector<yo_material> sh_area1_materials = {
        make_emission("light01", {10, 6, 3}, 0),
        make_emission("light02", {3, 6, 10}, 0),
    };

    // save plane scenes
    save_scene("sh01.obj", dirname,
               make_scene(sh_cams, sh_shapes + sh_light_shapes,
                          sh_materials + sh_light_materials));
    save_scene("sh02.obj", dirname,
               make_scene(sh_cams, sh_shapes + sh_light_shapes,
                          sh_textured_materials + sh_light_materials));
    save_scene("sh03.obj", dirname,
               make_scene(sh_cams, sh_shapes + sh_area_shapes,
                          sh_textured_materials + sh_area_materials));
    save_scene("sh04.obj", dirname,
               make_scene(sh_cams, sh_shapes + sh_area1_shapes,
                          sh_textured_materials + sh_area1_materials));

    // point and lines scene ------------------------------
    printf("generating points and lines scenes ...\n");
    ym_vector<yo_shape> ps_shapes = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), {0, 0, -4},
                        {0, 0, 0}, 6),
        transform_shape(make_shape("points01", "points", 18, ys_stype_rnpoints),
                        {0, 0.5f, 0}, ym_zero3f, 0.5)};
    ym_vector<yo_shape> ls_shapes = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), {0, 0, -4},
                        ym_zero3f, 6),
        transform_shape(make_shape("obj01", "obj", 6, ys_stype_uvsphere),
                        {1.25f, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(
            make_lines("lines01", "lines", 64 * 64 * 16, 4, 0.1f, 0, 0),
            {1.25f, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(make_shape("obj02", "obj", 6, ys_stype_uvsphere),
                        {0, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(
            make_lines("lines02", "lines", 64 * 64 * 16, 4, 0, 0.75f, 0),
            {0, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(make_shape("obj03", "obj", 6, ys_stype_uvsphere),
                        {-1.25f, 0.5f, 0}, ym_zero3f, 0.5),
        transform_shape(
            make_lines("lines03", "lines", 64 * 64 * 16, 4, 0, 0, 0.5f),
            {-1.25f, 0.5f, 0}, ym_zero3f, 0.5)};
    ym_vector<yo_material> pl_materials = {
        make_diffuse("floor", {0.2f, 0.2f, 0.2f}, 0),
        make_diffuse("obj", {0.2f, 0.2f, 0.2f}, 0),
        make_diffuse("points", {0.2f, 0.2f, 0.2f}, 0),
        make_diffuse("lines", {0.2f, 0.2f, 0.2f}, 0)};

    // save scenes
    save_scene("ps01.obj", dirname,
               make_scene(sh_cams, ps_shapes + sh_light_shapes,
                          pl_materials + sh_light_materials));
    save_scene("ps02.obj", dirname,
               make_scene(sh_cams, ps_shapes + sh_area_shapes,
                          pl_materials + sh_area_materials));
    save_scene("ps03.obj", dirname,
               make_scene(sh_cams, ps_shapes + sh_area1_shapes,
                          pl_materials + sh_area1_materials));
    save_scene("ls01.obj", dirname,
               make_scene(sh_cams, ls_shapes + sh_light_shapes,
                          pl_materials + sh_light_materials));
    save_scene("ls02.obj", dirname,
               make_scene(sh_cams, ls_shapes + sh_area_shapes,
                          pl_materials + sh_area_materials));
    save_scene("ls03.obj", dirname,
               make_scene(sh_cams, ls_shapes + sh_area1_shapes,
                          pl_materials + sh_area1_materials));

    // random obj scene --------------------------
    printf("generating random shapes scenes ...\n");
    ym_vector<yo_shape> rsh_shapes = make_random_shapes(32, 5, false);
    ym_vector<yo_material> rsh_materials = make_random_materials(32);

    // save plane scenes
    save_scene("rs01.obj", dirname,
               make_scene(sh_cams, rsh_shapes + sh_light_shapes,
                          rsh_materials + sh_light_materials));
    save_scene("rs02.obj", dirname,
               make_scene(sh_cams, rsh_shapes + sh_area_shapes,
                          rsh_materials + sh_area_materials));
    save_scene("rs03.obj", dirname,
               make_scene(sh_cams, rsh_shapes + sh_area1_shapes,
                          rsh_materials + sh_area1_materials));

    // xform obj scene ---------------------------
    printf("generating transform tests scenes ...\n");
    ym_vector<yo_shape> xsh_shapes = make_random_shapes(32, 5, true);
    ym_vector<yo_shape> xsh_light_shapes = {
        xform_shape(make_shape("light01", "light01", 0, ys_stype_points),
                    {0.7f, 4, 3}, ym_zero3f, 1),
        xform_shape(make_shape("light02", "light02", 0, ys_stype_points),
                    {-0.7f, 4, 3}, ym_zero3f, 1)};
    ym_vector<yo_shape> xsh_area_shapes = {
        xform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                    {4, 2, 4}, {0, 180 + 45, 0}, 2),
        xform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                    {-4, 2, 4}, {0, 180 - 45, 0}, 2)};

    // save plane scenes
    save_scene("xs01.obj", dirname,
               make_scene(sh_cams, xsh_shapes + xsh_light_shapes,
                          rsh_materials + sh_light_materials));
    save_scene("xs02.obj", dirname,
               make_scene(sh_cams, xsh_shapes + xsh_area_shapes,
                          rsh_materials + sh_area_materials));

    // env scene ------------------------------
    printf("generating envmaps scenes ...\n");
    ym_vector<yo_camera> env_cams = {
        make_camera("cam", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {0, 1.5f, 5}, {0, 0.5f, 0}, 0.5f, 0.1f)};
    ym_vector<yo_env> env_envs = {
        make_env("env", "env", {0, 0.5f, 0}, {0, 1.5f, 0})};
    ym_vector<yo_env> env_txt_envs = {
        make_env("env", "env_txt", {0, 0.5f, 0}, {-1.5f, 0.5f, 0})};
    ym_vector<yo_shape> env_shapes = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), {0, 0, -4},
                        ym_zero3f, 6),
        transform_shape(
            make_shape("obj01", "obj01", 5, ys_stype_uvflipcapsphere),
            {-1.25f, 0.5f, 0}, ym_zero3f, 0.5f),
        transform_shape(
            make_shape("obj02", "obj02", 4, ys_stype_uvspherizedcube),
            {0, 0.5f, 0}, ym_zero3f, 0.5f),
        transform_shape(make_shape("obj03", "obj03", 4, ys_stype_uvspherecube),
                        {1.25f, 0.5f, 0}, ym_zero3f, 0.5f)};
    ym_vector<yo_material> env_materials = {
        make_diffuse("floor", {0.2f, 0.2f, 0.2f}, 0),
        make_plastic("obj01", {0.5f, 0.2f, 0.2f}, 50, 0),
        make_plastic("obj02", {0.2f, 0.5f, 0.2f}, 100, 0),
        make_plastic("obj03", {0.2f, 0.2f, 0.5f}, 500, 0),
        make_emission("env", {1, 1, 1}, 0),
        make_emission("env_txt", {1, 1, 1}, "env.hdr")};
    ym_vector<yo_shape> env_light_shapes = {transform_shape(
        make_shape("env_sphere", "env", 6, ys_stype_uvflippedsphere),
        {0, 0.5f, 0}, {-90, 0, 0}, 10000)};
    ym_vector<yo_material> env_light_materials = {
        make_emission("env", {1, 1, 1}, 0)};
    ym_vector<yo_shape> env_txt_light_shapes = {transform_shape(
        make_shape("env_sphere", "env_txt", 6, ys_stype_uvflippedsphere),
        {0, 0.5f, 0}, {-90, 0, 0}, 10000)};
    ym_vector<yo_material> env_txt_light_materials = {
        make_emission("env_txt", {1, 1, 1}, "env.hdr")};

    // save plane scenes
    save_scene("es01.obj", dirname,
               make_scene(env_cams, env_shapes, env_materials, env_envs));
    save_scene("es02.obj", dirname,
               make_scene(env_cams, env_shapes + env_light_shapes,
                          env_materials + env_light_materials));
    save_scene("es03.obj", dirname,
               make_scene(env_cams, env_shapes, env_materials, env_txt_envs));
    save_scene("es04.obj", dirname,
               make_scene(env_cams, env_shapes + env_txt_light_shapes,
                          env_materials + env_txt_light_materials));

    // cornell box ------------------------------
    // http://graphics.cs.williams.edu/data
    // http://www.graphics.cornell.edu/online/box/data.html
    printf("generating cornell box scenes ...\n");
    ym_vector<yo_camera> cb_cams = {
        make_camera("cam", {0, 1, 4}, {0, 1, 0}, 0.7f, 0)};
    ym_vector<yo_shape> cb_box = {
        transform_shape(make_shape("floor", "white", 0, ys_stype_uvquad),
                        ym_zero3f, {-90, 0, 0}, 1),
        transform_shape(make_shape("ceiling", "white", 0, ys_stype_uvquad),
                        {0, 2, 0}, {90, 0, 0}, 1),
        transform_shape(make_shape("back", "white", 0, ys_stype_uvquad),
                        {0, 1, -1}, ym_zero3f, 1),
        transform_shape(make_shape("back", "green", 0, ys_stype_uvquad),
                        {+1, 1, 0}, {0, -90, 0}, 1),
        transform_shape(make_shape("back", "red", 0, ys_stype_uvquad),
                        {-1, 1, 0}, {0, 90, 0}, 1),
        transform_shape(make_shape("tallbox", "white", 0, ys_stype_uvcube),
                        {-0.33f, 0.6f, -0.29f}, {0, 15, 0}, {0.3f, 0.6f, 0.3f}),
        transform_shape(make_shape("shortbox", "white", 0, ys_stype_uvcube),
                        {0.33f, 0.3f, 0.33f}, {0, -15, 0}, {0.3f, 0.3f, 0.3f}),
        transform_shape(make_shape("light", "light", 0, ys_stype_uvquad),
                        {0, 1.999f, 0}, {90, 0, 0}, {0.25f, 0.25f, 0.25f})};
    ym_vector<yo_material> cb_materials = {
        make_diffuse("white", {0.725f, 0.71f, 0.68f}, 0),
        make_diffuse("red", {0.63f, 0.065f, 0.05f}, 0),
        make_diffuse("green", {0.14f, 0.45f, 0.091f}, 0),
        make_emission("light", {17, 12, 4}, 0),
    };

    // save scene
    save_scene("cb01.obj", dirname, make_scene(cb_cams, cb_box, cb_materials));

    // rigid body scenes ------------------------
    printf("generating rigid body scenes ...\n");
    ym_vector<yo_camera> rb_cams = {
        make_camera("cam", {5, 5, 5}, {0, 0.5f, 0}, 0.5f, 0),
        make_camera("cam_dof", {5, 5, 5}, {0, 0.5f, 0}, 0.5f, 0.1f)};
    ym_vector<yo_shape> rb_shapes = {
        xform_shape(make_shape("floor", "floor", 4, ys_stype_uvcube),
                    {0, -0.5f, 0}, {0, 0, 0}, {6, 0.5f, 6}),
        xform_shape(make_shape("obj01", "obj", 2, ys_stype_uvcube),
                    {-1.25f, 0.5f, 0}, {0, 0, 0}, 0.5f),
        xform_shape(make_shape("obj02", "obj", 3, ys_stype_uvspherecube),
                    {0, 1, 0}, {0, 0, 0}, 0.5f),
        xform_shape(make_shape("obj03", "obj", 2, ys_stype_uvcube),
                    {1.25f, 1.5f, 0}, {0, 0, 0}, 0.5f),
        xform_shape(make_shape("obj11", "obj", 2, ys_stype_uvcube),
                    {-1.25f, 0.5f, 1.5f}, {0, 45, 0}, 0.5f),
        xform_shape(make_shape("obj12", "obj", 3, ys_stype_uvspherecube),
                    {0, 1, 1.5f}, {45, 0, 0}, 0.5f),
        xform_shape(make_shape("obj13", "obj", 2, ys_stype_uvcube),
                    {1.25f, 1.5f, 1.5f}, {45, 0, 45}, 0.5f),
        xform_shape(make_shape("obj21", "obj", 2, ys_stype_uvcube),
                    {-1.25f, 0.5f, -1.5f}, {0, 0, 0}, 0.5f),
        xform_shape(make_shape("obj22", "obj", 3, ys_stype_uvspherecube),
                    {0, 1, -1.5f}, {22.5, 0, 0}, 0.5),
        xform_shape(make_shape("obj23", "obj", 2, ys_stype_uvcube),
                    {1.25f, 1.5f, -1.5f}, {22.5f, 0, 22.5f}, 0.5f)};
    ym_vector<yo_material> rb_materials = {
        make_diffuse("floor", {1, 1, 1}, "grid.png"),
        make_plastic("obj", {1, 1, 1}, 50, "checker.png")};
    ym_vector<yo_shape> rb_light_shapes = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_points),
                        {0.7f, 4, 3}, {0, 0, 0}, 1),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_points),
                        {-0.7f, 4, 3}, {0, 0, 0}, 1)};
    ym_vector<yo_material> rb_light_materials = {
        make_emission("light01", {100, 100, 100}, 0),
        make_emission("light02", {100, 100, 100}, 0),
    };

    ym_vector<yo_shape> rb2_shapes = {
        xform_shape(make_shape("floor", "floor", 2, ys_stype_uvcube),
                    {0, -2.5, 0}, {30, 0, 0}, {6, 0.5f, 6}),
        rb_shapes[1], rb_shapes[2], rb_shapes[3], rb_shapes[4], rb_shapes[5],
        rb_shapes[6], rb_shapes[7], rb_shapes[8], rb_shapes[9],
    };

    ym_vector<yo_shape> rb3_shapes = make_random_rigid_shapes(128, 1);

    // save plane scenes
    save_scene("rb01.obj", dirname,
               make_scene(rb_cams, rb_shapes + rb_light_shapes,
                          rb_materials + rb_light_materials));
    save_scene("rb02.obj", dirname,
               make_scene(rb_cams, rb2_shapes + rb_light_shapes,
                          rb_materials + rb_light_materials));
    save_scene("rb03.obj", dirname,
               make_scene(rb_cams, rb3_shapes + rb_light_shapes,
                          rb_materials + rb_light_materials));

    // textures ---------------------------------
    printf("generating simple textures ...\n");
    save_image("grid.png", dirname, make_grid(512).data(), 512);
    save_image("checker.png", dirname, make_checker(512).data(), 512);
    save_image("rchecker.png", dirname, make_rchecker(512).data(), 512);
    save_image("colored.png", dirname, make_colored(512).data(), 512);
    save_image("rcolored.png", dirname, make_rcolored(512).data(), 512);
    printf("generating envmaps textures ...\n");
    save_image_hdr("env.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8f, 8,
                                   ym_vec3f{0.2f, 0.2f, 0.2f}, 1 / powf(2, 6),
                                   true)
                       .data(),
                   1024, 512);
    save_image_hdr("env01.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8f, 8,
                                   ym_vec3f{0.2f, 0.2f, 0.2f}, 1 / powf(2, 6),
                                   true)
                       .data(),
                   1024, 512);
}
