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
#include "../yocto/yocto_cmdline.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_math.h"
#include "../yocto/yocto_obj.h"
#include "../yocto/yocto_shape.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "ext/stb_image_write.h"

#include "ext/sunsky/ArHosekSkyModel.c"
#include "ext/sunsky/ArHosekSkyModel.h"

ym_vec3f
v(float x, float y, float z) {
    return (ym_vec3f){ x, y, z };
}
ym_vec3f
v1(float x) {
    return (ym_vec3f){ x, x, x };
}
ym_vec3f
zero() {
    return (ym_vec3f){ 0, 0, 0 };
}

ym_mat4f
xform(ym_vec3f pos, ym_vec3f rot, ym_vec3f s) {
    ym_mat4f xf = ym_identity4f();
    xf = ym_mmul4f(ym_scaling4f(s), xf);
    xf = ym_mmul4f(ym_rotationx4f(rot.x * ym_pif / 180), xf);
    xf = ym_mmul4f(ym_rotationy4f(rot.y * ym_pif / 180), xf);
    xf = ym_mmul4f(ym_rotationz4f(rot.z * ym_pif / 180), xf);
    xf = ym_mmul4f(ym_translation4f(pos), xf);
    return xf;
}

yo_shape
transform_shape(yo_shape shape, ym_vec3f pos, ym_vec3f rot, ym_vec3f s) {
    ym_mat4f xf = xform(pos, rot, s);
    for (int i = 0; i < shape.nverts; i++) {
        ym_vec3f* pos = (ym_vec3f*)shape.pos + i;
        *pos = ym_transform_point3f(xf, *pos);
        ym_vec3f* norm = (ym_vec3f*)shape.norm + i;
        *norm = ym_transform_direction3f(xf, *norm);
    }
    return shape;
}

yo_shape
xform_shape(yo_shape shape, ym_vec3f pos, ym_vec3f rot, ym_vec3f s) {
    ym_mat4f scale = ym_scaling4f(s);
    for (int i = 0; i < shape.nverts; i++) {
        ym_vec3f* pos = (ym_vec3f*)shape.pos + i;
        *pos = ym_transform_point3f(scale, *pos);
    }
    ym_mat4f xf = xform(pos, rot, v1(1));
    if (shape.xform) {
        xf = ym_mmul4f(xf, *(ym_mat4f*)shape.xform);
    }
    shape.xform = (float*)calloc(16, sizeof(float));
    *(ym_mat4f*)shape.xform = xf;
    return shape;
}

yo_shape
make_shape(char* name, char* matname, int l, int stype) {
    int etype = yo_etype_quad;
    float params[2] = { 0.75f, 0.75f };
    if (stype == ys_stype_points || stype == ys_stype_rnpoints) {
        etype = yo_etype_point;
        params[0] = 0.001;
        params[0] = 0.0025;
    }
    yo_shape shape = { 0 };
    shape.name = (char*)strdup(name);
    shape.matname = (char*)strdup(matname);
    shape.etype = etype;
    ys_make_stdshape(stype, l, etype, params, &shape.nelems, &shape.elem,
                     &shape.nverts, &shape.pos, &shape.norm, &shape.texcoord,
                     0);
    return shape;
}

struct floor_params {
    float size;
    float cur_power;
};
typedef struct floor_params floor_params;

void
floor_pos(void* ctx, int vid, const float uv[2], float pos[3]) {
    floor_params* params = (floor_params*)ctx;
    float x = 2 * uv[0] - 1;
    float y = 2 * (1 - uv[1]) - 1;
    if (y >= 0 || !params->cur_power) {
        pos[0] = x;
        pos[1] = 0;
        pos[2] = y;
    } else {
        pos[0] = x;
        pos[1] = powf(-y, params->cur_power);
        pos[2] = y;
    }
}

void
floor_norm(void* ctx, int vid, const float uv[2], float norm[3]) {
    norm[0] = 0;
    norm[1] = 1;
    norm[2] = 0;
}

void
floor_texcoord(void* ctx, int vid, const float uv[2], float texcoord[2]) {
    floor_params* params = (floor_params*)ctx;
    texcoord[0] = uv[0] * params->size;
    texcoord[1] = uv[1] * params->size;
}

yo_shape
make_floor(char* name, char* matname, float s, float p, int l) {
    int n = round(powf(2, l));
    yo_shape shape = { 0 };
    shape.name = (char*)strdup(name);
    shape.matname = (char*)strdup(matname);
    shape.etype = yo_etype_quad;
    floor_params params = { s, p };
    ys_make_uvshape(n, n, &params, ys_etype_quad, 3, (int[3]){ 3, 3, 2 },
                    (vfunc[3]){ floor_pos, floor_norm, floor_texcoord },
                    &shape.nelems, &shape.elem, &shape.nverts,
                    (float* * [3]){ &shape.pos, &shape.norm, &shape.texcoord });
    if (p) {
        ys_compute_normals(shape.nelems, shape.elem, shape.etype, shape.nverts,
                           shape.pos, shape.norm, true);
    }
    return shape;
}

yo_material
make_material(char* name, ym_vec3f ke, ym_vec3f kd, ym_vec3f ks, float ns,
              char* ke_txt, char* kd_txt, char* ks_txt) {
    yo_material mat = { 0 };
    mat.name = (char*)strdup(name);
    *(ym_vec3f*)mat.ke = ke;
    *(ym_vec3f*)mat.kd = kd;
    *(ym_vec3f*)mat.ks = ks;
    mat.ns = ns;
    mat.ke_txt = ke_txt;
    mat.kd_txt = kd_txt;
    mat.ks_txt = ks_txt;
    return mat;
}

yo_material
make_emission(char* name, ym_vec3f ke, char* txt) {
    return make_material(name, ke, ym_zero3f(), ym_zero3f(), 0, txt, 0, 0);
}

yo_material
make_diffuse(char* name, ym_vec3f kd, char* txt) {
    return make_material(name, ym_zero3f(), kd, ym_zero3f(), 0, 0, txt, 0);
}

yo_material
make_plastic(char* name, ym_vec3f kd, float n, char* txt) {
    return make_material(name, ym_zero3f(), kd, v(0.04, 0.04, 0.04), n, 0, txt,
                         0);
}

yo_material
make_metal(char* name, ym_vec3f kd, float n, char* txt) {
    return make_material(name, ym_zero3f(), ym_zero3f(), kd, n, 0, 0, txt);
}

yo_camera
make_camera(char* name, ym_vec3f from, ym_vec3f to, float h, float a) {
    yo_camera cam = { 0 };
    cam.name = (char*)strdup(name);
    *(ym_vec3f*)cam.from = from;
    *(ym_vec3f*)cam.to = to;
    *(ym_vec3f*)cam.up = (ym_vec3f){ 0, 1, 0 };
    cam.aperture = a;
    cam.height = h;
    cam.width = cam.height * 16.0f / 9.0f;
    return cam;
}

yo_env
make_env(char* name, char* matname, ym_vec3f from, ym_vec3f to) {
    yo_env env = { 0 };
    env.name = (char*)strdup(name);
    env.matname = (char*)strdup(matname);
    *(ym_vec3f*)env.from = from;
    *(ym_vec3f*)env.to = to;
    *(ym_vec3f*)env.up = (ym_vec3f){ 0, 1, 0 };
    return env;
}

static inline float
rng_nextf(unsigned int* state) {
    *state = (1103515245U * *state + 12345U) % 2147483648U;
    return fminf((float)*state / (float)2147483648U, 0.999999f);
}

typedef struct {
    int n;
    float r, c, s;
    ym_vec3f* base;
    ym_vec3f* dir;
    float* ln;
} line_params;

void
lines_pos(void* ctx, int vid, const float uv[2], float pos_[3]) {
    unsigned int rn = 0;
    line_params* params = (line_params*)ctx;
    int i = ym_clamp(uv[1] * (params->n), 0, params->n - 1);
    ym_vec3f* pos = (ym_vec3f*)pos_;
    *pos = ym_smul3f(params->base[i], 1 + uv[0] * params->ln[i]);
    if (params->r) {
        *pos =
            ym_sum3f(*pos, (ym_vec3f){ params->r * (0.5f - rng_nextf(&rn)),
                                       params->r * (0.5f - rng_nextf(&rn)),
                                       params->r * (0.5f - rng_nextf(&rn)) });
    }
    if (params->s && uv[0]) {
        ym_mat4f rotation = ym_rotationy4f(params->s * uv[0] * uv[0]);
        *pos = ym_affine_transform_point3f(rotation, *pos);
    }
    int nc = 128;
    if (params->c && i > nc) {
        int c = 0;
        float md = HUGE_VALF;
        for (int k = 0; k < nc; k++) {
            float d = ym_dist3f(params->base[i], params->base[k]);
            if (d < md) {
                md = d;
                c = k;
            }
        }
        ym_vec3f cpos = ym_smul3f(params->base[c], 1 + uv[0] * params->ln[c]);
        *pos = ym_sum3f(ym_smul3f(*pos, 1 - params->c * uv[0] * uv[0]),
                        ym_smul3f(cpos, params->c * uv[0] * uv[0]));
    }
}

void
lines_norm(void* ctx, int vid, const float uv[2], float norm[3]) {
    *(ym_vec3f*)norm = (ym_vec3f){ 0, 1, 0 };
}

void
lines_radius(void* ctx, int vid, const float uv[2], float* radius) {
    *radius = 0.001 + 0.001 * (1 - uv[0]);
}

yo_shape
make_lines(char* name, char* matname, int num, int n, float r, float c,
           float s) {
    yo_shape shape = { 0 };
    shape.name = (char*)strdup(name);
    shape.matname = (char*)strdup(matname);
    shape.etype = yo_etype_line;

    unsigned int rn = 0;
    line_params params = { num + 1,
                           r,
                           c,
                           s,
                           (ym_vec3f*)calloc(num + 1, sizeof(ym_vec3f)),
                           (ym_vec3f*)calloc(num + 1, sizeof(ym_vec3f)),
                           (float*)calloc(num + 1, sizeof(float)) };
    for (int i = 0; i <= num; i++) {
        float z = -1 + 2 * rng_nextf(&rn);
        float r = sqrtf(ym_fclampf(1 - z * z, 0, 1));
        float phi = 2 * ym_pif * rng_nextf(&rn);
        params.base[i] = (ym_vec3f){ r * cosf(phi), r * sinf(phi), z };
        params.dir[i] = params.base[i];
        params.ln[i] = 0.15f + 0.15f * rng_nextf(&rn);
    }

    ys_make_uvshape(n, num, &params, shape.etype, 4, (int[4]){ 3, 3, 2, 1 },
                    (vfunc[4]){ lines_pos, lines_norm, 0, lines_radius },
                    &shape.nelems, &shape.elem, &shape.nverts,
                    (float* * [4]){ &shape.pos, &shape.norm, &shape.texcoord,
                                    &shape.radius });

    ys_compute_normals(shape.nelems, shape.elem, shape.etype, shape.nverts,
                       shape.pos, shape.norm, true);

    return shape;
}

yo_shape*
make_random_shapes(int nshapes, int l, bool use_xform) {
    yo_shape* shapes = (yo_shape*)calloc(nshapes, sizeof(yo_shape));
    shapes[0] = transform_shape(make_floor("floor", "floor", 6, 4, 6),
                                v(0, 0, -4), v(0, 0, 0), v1(6));
    ym_vec3f pos[1024];
    float radius[1024];
    int levels[1024];

    unsigned int rn = 0;
    for (int i = 1; i < nshapes; i++) {
        bool done = false;
        while (!done) {
            float x = -2 + 4 * rng_nextf(&rn);
            float z = 1 - 3 * rng_nextf(&rn);
            radius[i] = 0.15f + ((1 - z) / 3) * ((1 - z) / 3) * 0.5f;
            pos[i] = (ym_vec3f){ x, radius[i], z };
            levels[i] = (int)round(log2f(powf(2, l) * radius[i] / 0.5));
            done = true;
            for (int j = 1; j < i && done; j++) {
                if (ym_dist3f(pos[i], pos[j]) < radius[i] + radius[j])
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
        int stype = stypes[(int)(rng_nextf(&rn) * 3)];
        if (stype == ys_stype_uvflipcapsphere) levels[i]++;
        shapes[i] = make_shape(name, name, levels[i], stype);
        if (use_xform) {
            shapes[i] = xform_shape(shapes[i], pos[i], (ym_vec3f){ 0, 0, 0 },
                                    v1(radius[i]));
        } else {
            shapes[i] = transform_shape(shapes[i], pos[i],
                                        (ym_vec3f){ 0, 0, 0 }, v1(radius[i]));
        }
    }

    return shapes;
}

yo_material*
make_random_materials(int nshapes) {
    yo_material* materials = (yo_material*)calloc(nshapes, sizeof(yo_material));
    materials[0] = make_diffuse("floor", v(1, 1, 1), "grid.png");

    unsigned int rn = 0;
    for (int i = 1; i < nshapes; i++) {
        char name[1024];
        sprintf(name, "obj%02d", i);
        char* txt = 0;
        if (rng_nextf(&rn) < 0.5f) {
            char* txts[6] = { 0,
                              "grid.png",
                              "checker.png",
                              "rchecker.png",
                              "colored.png",
                              "rcolored.png" };
            txt = txts[(int)(rng_nextf(&rn) * 6)];
        }
        ym_vec3f c = (txt) ? (ym_vec3f){ 1, 1, 1 }
                           : (ym_vec3f){ 0.2f + 0.3f * rng_nextf(&rn),
                                         0.2f + 0.3f * rng_nextf(&rn),
                                         0.2f + 0.3f * rng_nextf(&rn) };
        float rs = 0.01f + 0.25f * rng_nextf(&rn);
        float ns = 2 / (rs * rs) - 2;
        int mt = rng_nextf(&rn) * 4;
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

#define arraydup(a, n) memcpy(malloc(sizeof(*a) * n), a, sizeof(*a) * n)

yo_scene
make_scene(int ncameras, yo_camera* cameras, int nenvs, yo_env* envs, int nobjs,
           yo_shape* obj_shapes, int nmaterials, yo_material* obj_materials,
           int nlights, yo_shape* light_shapes, yo_material* light_materials) {
    yo_scene scene = { 0 };
    scene.ncameras = ncameras;
    scene.cameras = (yo_camera*)calloc(ncameras, sizeof(yo_camera));
    for (int i = 0; i < ncameras; i++) scene.cameras[i] = cameras[i];
    scene.nenvs = nenvs;
    scene.envs = (yo_env*)calloc(nenvs, sizeof(yo_env));
    for (int i = 0; i < nenvs; i++) scene.envs[i] = envs[i];
    scene.nshapes = nobjs + nlights;
    scene.shapes = (yo_shape*)calloc(scene.nshapes, sizeof(yo_shape));
    for (int i = 0; i < nobjs; i++) scene.shapes[i] = obj_shapes[i];
    for (int i = 0; i < nlights; i++) scene.shapes[i + nobjs] = light_shapes[i];
    scene.nmaterials = nmaterials + nlights;
    scene.materials =
        (yo_material*)calloc(scene.nmaterials, sizeof(yo_material));
    for (int i = 0; i < nmaterials; i++) scene.materials[i] = obj_materials[i];
    for (int i = 0; i < nlights; i++)
        scene.materials[i + nmaterials] = light_materials[i];
    return scene;
}

typedef unsigned char ubyte;
typedef struct { ubyte r, g, b, a; } rgba;

rgba*
make_grid(int s) {
    rgba* pixels = (rgba*)calloc(s * s, 4);
    int g = 64;
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            if (i % g == 0 || i % g == g - 1 || j % g == 0 || j % g == g - 1)
                pixels[j * s + i] = (rgba){ 90, 90, 90, 255 };
            else
                pixels[j * s + i] = (rgba){ 128, 128, 128, 255 };
        }
    }
    return pixels;
}

rgba*
make_checker(int s) {
    rgba* pixels = (rgba*)calloc(s * s, 4);
    for (int j = 0; j < s; j++) {
        for (int i = 0; i < s; i++) {
            if ((i / 64 + j / 64) % 2)
                pixels[j * s + i] = (rgba){ 90, 90, 90, 255 };
            else
                pixels[j * s + i] = (rgba){ 128, 128, 128, 255 };
        }
    }
    return pixels;
}

// http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
rgba
hsv_to_rgb(ubyte h, ubyte s, ubyte v) {
    rgba rgb = { 0, 0, 0, 255 };
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

rgba*
make_rcolored(int s) {
    rgba* pixels = (rgba*)calloc(s * s, 4);
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

rgba*
make_colored(int s) {
    rgba* pixels = (rgba*)calloc(s * s, 4);
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

rgba*
make_rchecker(int s) {
    rgba* pixels = (rgba*)calloc(s * s, 4);
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
            pixels[j * s + i] = (rgba){ pv, pv, pv, 255 };
        }
    }
    return pixels;
}

#define sqr(x) ((x) * (x))

static inline float
sunsky_PerezFunctionFromDensity(const float* lambdas, float invDenLum,
                                float cosTheta, float gamma, float cosGamma) {
    float num = (1.0 + lambdas[0] * exp(lambdas[1] / cosTheta)) *
                (1.0 + lambdas[2] * exp(lambdas[3] * gamma) +
                 lambdas[4] * sqr(cosGamma));

    return num * invDenLum;
}

float*
make_sunsky2_hdr(int w, int h, float sun_theta, float turbidity, float scale) {
    float* rgba = calloc(4 * w * h, sizeof(float));

    float theta2 = sun_theta * sun_theta;
    float theta3 = sun_theta * sun_theta * sun_theta;
    float T = turbidity;
    float T2 = T * T;

    float chi = (4.0 / 9.0 - T / 120.0) * (ym_pi - 2.0 * sun_theta);

    // mZenith stored as xyY
    float mZenith[3] = {
        (0.00165 * theta3 - 0.00374 * theta2 + 0.00208 * sun_theta + 0.0) * T2 +
            (-0.02902 * theta3 + 0.06377 * theta2 - 0.03202 * sun_theta +
             0.00394) *
                T +
            (0.11693 * theta3 - 0.21196 * theta2 + 0.06052 * sun_theta +
             0.25885),
        (0.00275 * theta3 - 0.00610 * theta2 + 0.00316 * sun_theta + 0.0) * T2 +
            (-0.04214 * theta3 + 0.08970 * theta2 - 0.04153 * sun_theta +
             0.00515) *
                T +
            (0.15346 * theta3 - 0.26756 * theta2 + 0.06669 * sun_theta +
             0.26688),
        (4.0453 * T - 4.9710) * tan(chi) - .2155 * T + 2.4192
    };

    // Perez constants
    float mPerez_Y[5] = { 0.17872 * T - 1.46303 - 0.35540 * T + 0.42749,
                          -0.02266 * T + 5.32505, 0.12064 * T - 2.57705,
                          -0.06696 * T + 0.37027 };

    float mPerez_x[5] = { -0.01925 * T - 0.25922, -0.06651 * T + 0.00081,
                          -0.00041 * T + 0.21247, -0.06409 * T - 0.89887,
                          -0.00325 * T + 0.04517 };

    float mPerez_y[5] = { -0.01669 * T - 0.26078, -0.09495 * T + 0.00921,
                          -0.00792 * T + 0.21023, -0.04405 * T - 1.65369,
                          -0.01092 * T + 0.05291 };

    // initialize sun-constant parts of the Perez functions
    float mPerezInvDen[3] = {
        mZenith[0] / ((1.0 + mPerez_x[0] * exp(mPerez_x[1])) *
                      (1.0 + mPerez_x[2] * exp(mPerez_x[3] * sun_theta) +
                       mPerez_x[4] * sqr(cos(sun_theta)))),
        mZenith[1] / ((1.0 + mPerez_y[0] * exp(mPerez_y[1])) *
                      (1.0 + mPerez_y[2] * exp(mPerez_y[3] * sun_theta) +
                       mPerez_y[4] * sqr(cos(sun_theta)))),
        mZenith[2] / ((1.0 + mPerez_Y[0] * exp(mPerez_Y[1])) *
                      (1.0 + mPerez_Y[2] * exp(mPerez_Y[3] * sun_theta) +
                       mPerez_Y[4] * sqr(cos(sun_theta))))
    };

    float sun_phi = ym_pif;
    ym_vec3f sun_w =
        (ym_vec3f){ cosf(sun_phi) * sinf(sun_theta),
                    sinf(sun_phi) * sinf(sun_theta), cosf(sun_theta) };

    const ym_vec3f kXYZToR = { 2.80298, -1.18735, -0.437286 };
    const ym_vec3f kXYZToG = { -1.07909, 1.97927, 0.0423073 };
    const ym_vec3f kXYZToB = { 0.0746015, -0.248513, 1.08196 };

    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            float theta = ym_pif * (j + 0.5f) / h;
            float phi = 2 * ym_pif * (i + 0.5f) / w;
            ym_vec3f pw = (ym_vec3f){ cosf(phi) * sinf(theta),
                                      sinf(phi) * sinf(theta), cosf(theta) };
            float cosTheta = pw.z;
            float cosGamma = ym_dot3f(sun_w, pw);
            float gamma = acos(cosGamma);

            if (cosTheta < 0.0f) cosTheta = 0.0f;

            float xyY[3] = {
                sunsky_PerezFunctionFromDensity(mPerez_x, mPerezInvDen[0],
                                                cosTheta, gamma, cosGamma),
                sunsky_PerezFunctionFromDensity(mPerez_y, mPerezInvDen[1],
                                                cosTheta, gamma, cosGamma),
                sunsky_PerezFunctionFromDensity(mPerez_Y, mPerezInvDen[2],
                                                cosTheta, gamma, cosGamma)
            };

            ym_vec3f cieXYZ = { xyY[0] * xyY[2] / xyY[1],
                                xyY[1] * xyY[2] / xyY[1],
                                (1.0 - xyY[0] - xyY[1]) * xyY[2] / xyY[1] };
            ym_vec3f rgb = { ym_dot3f(kXYZToR, cieXYZ),
                             ym_dot3f(kXYZToG, cieXYZ),
                             ym_dot3f(kXYZToB, cieXYZ) };

            rgba[(j * w + i) * 4 + 0] = scale * rgb.x;
            rgba[(j * w + i) * 4 + 1] = scale * rgb.y;
            rgba[(j * w + i) * 4 + 2] = scale * rgb.z;
            rgba[(j * w + i) * 4 + 3] = 1;
        }
    }

    return rgba;
}

float*
make_sunsky_hdr(int w, int h, float sun_theta, float turbidity, ym_vec3f ground,
                float scale, bool include_ground) {
    float* rgba = calloc(4 * w * h, sizeof(float));
    ArHosekSkyModelState* skymodel_state[3] = {
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
        arhosek_rgb_skymodelstate_alloc_init(turbidity, ground.x, sun_theta),
    };
    float sun_phi = ym_pif;
    ym_vec3f sun_w =
        (ym_vec3f){ cosf(sun_phi) * sinf(sun_theta),
                    sinf(sun_phi) * sinf(sun_theta), cosf(sun_theta) };
    for (int j = 0; j < h; j++) {
        for (int i = 0; i < w; i++) {
            float theta = ym_pif * (j + 0.5f) / h;
            float phi = 2 * ym_pif * (i + 0.5f) / w;
            if (include_ground)
                theta = ym_fclampf(theta, 0, ym_pif / 2 - 0.001);
            ym_vec3f pw = (ym_vec3f){ cosf(phi) * sinf(theta),
                                      sinf(phi) * sinf(theta), cosf(theta) };
            float gamma = acosf(ym_fclampf(ym_dot3f(sun_w, pw), -1, 1));
            ym_vec3f sky = { arhosek_tristim_skymodel_radiance(
                                 skymodel_state[0], theta, gamma, 0),
                             arhosek_tristim_skymodel_radiance(
                                 skymodel_state[1], theta, gamma, 1),
                             arhosek_tristim_skymodel_radiance(
                                 skymodel_state[2], theta, gamma, 2) };
            rgba[(j * w + i) * 4 + 0] = scale * sky.x;
            rgba[(j * w + i) * 4 + 1] = scale * sky.y;
            rgba[(j * w + i) * 4 + 2] = scale * sky.z;
            rgba[(j * w + i) * 4 + 3] = 1;
        }
    }
    arhosekskymodelstate_free(skymodel_state[0]);
    arhosekskymodelstate_free(skymodel_state[1]);
    arhosekskymodelstate_free(skymodel_state[2]);
    return rgba;
}

void
save_image(const char* filename, const char* dirname, const rgba* pixels,
           int s) {
    char path[4096];
    sprintf(path, "%s/%s", dirname, filename);
    stbi_write_png(path, s, s, 4, pixels, s * 4);
}

void
save_image_hdr(const char* filename, const char* dirname, const float* pixels,
               int w, int h) {
    char path[4096];
    sprintf(path, "%s/%s", dirname, filename);
    stbi_write_hdr(path, w, h, 4, pixels);
}

void
save_scene(const char* filename, const char* dirname, const yo_scene scene) {
    char path[4096];
    sprintf(path, "%s/%s", dirname, filename);
    yo_save_obj(path, &scene, true);
}

int
main(int argc, const char** argv) {
    // command line params
    yc_parser* parser = yc_init_parser(argc, argv, "make tests");
    const char* dirname =
        yc_parse_args(parser, "dirname", "directory name", ".", true);
    yc_done_parser(parser);

    // make directories
    char cmd[4096];
#ifndef _MSC_VER
    sprintf(cmd, "mkdir -p %s", dirname);
#else
    sprintf(cmd, "mkdir %s", dirname);
#endif
    system(cmd);

    // floor scene ------------------------------
    printf("generating simple scenes ...\n");
    yo_camera sh_cams[2] = {
        make_camera("cam", v(0, 1.5, 5), v(0, 0.5, 0), 0.5, 0),
        make_camera("cam_dof", v(0, 1.5, 5), v(0, 0.5, 0), 0.5, 0.1)
    };
    yo_shape sh_shapes[4] = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), v(0, 0, -4),
                        v(0, 0, 0), v1(6)),
        transform_shape(
            make_shape("obj01", "obj01", 5, ys_stype_uvflipcapsphere),
            v(-1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(
            make_shape("obj02", "obj02", 4, ys_stype_uvspherizedcube),
            v(0, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(make_shape("obj03", "obj03", 4, ys_stype_uvspherecube),
                        v(1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5))
    };
    yo_material sh_materials[4] = {
        make_diffuse("floor", v(0.2, 0.2, 0.2), 0),
        make_plastic("obj01", v(0.5, 0.2, 0.2), 50, 0),
        make_plastic("obj02", v(0.2, 0.5, 0.2), 100, 0),
        make_plastic("obj03", v(0.2, 0.2, 0.5), 500, 0)
    };
    yo_shape sh_light_shapes[2] = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_points),
                        v(0.7, 4, 3), v(0, 0, 0), v1(1)),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_points),
                        v(-0.7, 4, 3), v(0, 0, 0), v1(1))
    };
    yo_material sh_light_materials[2] = {
        make_emission("light01", v(100, 100, 100), 0),
        make_emission("light02", v(100, 100, 100), 0),
    };

    // textures
    yo_material sh_textured_materials[4] = {
        make_diffuse("floor", v(1, 1, 1), "grid.png"),
        make_plastic("obj01", v(1, 1, 1), 50, "rcolored.png"),
        make_plastic("obj02", v(1, 1, 1), 100, "checker.png"),
        make_plastic("obj03", v(1, 1, 1), 500, "colored.png")
    };

    // area lights
    yo_shape sh_area_shapes[2] = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                        v(4, 2, 4), v(0, 180 + 45, 0), v1(2)),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                        v(-4, 2, 4), v(0, 180 - 45, 0), v1(2))
    };
    yo_material sh_area_materials[2] = {
        make_emission("light01", v(10, 10, 10), 0),
        make_emission("light02", v(10, 10, 10), 0),
    };

    // colored area lights
    yo_shape sh_area1_shapes[2] = {
        transform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                        v(4, 2, 4), v(0, 180 + 45, 0), v1(2)),
        transform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                        v(-4, 2, 4), v(0, 180 - 45, 0), v1(2))
    };
    yo_material sh_area1_materials[2] = {
        make_emission("light01", v(10, 6, 3), 0),
        make_emission("light02", v(3, 6, 10), 0),
    };

    // save plane scenes
    save_scene("sh01.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 4, sh_shapes, 4, sh_materials, 2,
                          sh_light_shapes, sh_light_materials));
    save_scene("sh02.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 4, sh_shapes, 4,
                          sh_textured_materials, 2, sh_light_shapes,
                          sh_light_materials));
    save_scene("sh03.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 4, sh_shapes, 4,
                          sh_textured_materials, 2, sh_area_shapes,
                          sh_area_materials));
    save_scene("sh04.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 4, sh_shapes, 4,
                          sh_textured_materials, 2, sh_area1_shapes,
                          sh_area1_materials));

    // point and lines scene ------------------------------
    printf("generating points and lines scenes ...\n");
    yo_shape ps_shapes[2] = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), v(0, 0, -4),
                        v(0, 0, 0), v1(6)),
        transform_shape(make_shape("points01", "points", 18, ys_stype_rnpoints),
                        v(0, 0.5f, 0), v(0, 0, 0), v1(0.5))
    };
    yo_shape ls_shapes[7] = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), v(0, 0, -4),
                        v(0, 0, 0), v1(6)),
        transform_shape(make_shape("obj01", "obj", 6, ys_stype_uvsphere),
                        v(1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(
            make_lines("lines01", "lines", 64 * 64 * 16, 4, 0.1, 0, 0),
            v(1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(make_shape("obj02", "obj", 6, ys_stype_uvsphere),
                        v(0, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(
            make_lines("lines02", "lines", 64 * 64 * 16, 4, 0, 0.75, 0),
            v(0, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(make_shape("obj03", "obj", 6, ys_stype_uvsphere),
                        v(-1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(
            make_lines("lines03", "lines", 64 * 64 * 16, 4, 0, 0, 0.5),
            v(-1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5))
    };
    yo_material pl_materials[4] = { make_diffuse("floor", v(0.2, 0.2, 0.2), 0),
                                    make_diffuse("obj", v(0.2, 0.2, 0.2), 0),
                                    make_diffuse("points", v(0.2, 0.2, 0.2), 0),
                                    make_diffuse("lines", v(0.2, 0.2, 0.2),
                                                 0) };

    // save scenes
    save_scene("ps01.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 2, ps_shapes, 4, pl_materials, 2,
                          sh_light_shapes, sh_light_materials));
    save_scene("ps02.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 2, ps_shapes, 4, pl_materials, 2,
                          sh_area_shapes, sh_area_materials));
    save_scene("ps03.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 2, ps_shapes, 4, pl_materials, 2,
                          sh_area1_shapes, sh_area1_materials));
    save_scene("ls01.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 7, ls_shapes, 4, pl_materials, 2,
                          sh_light_shapes, sh_light_materials));
    save_scene("ls02.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 7, ls_shapes, 4, pl_materials, 2,
                          sh_area_shapes, sh_area_materials));
    save_scene("ls03.obj", dirname,
               make_scene(2, sh_cams, 0, 0, 7, ls_shapes, 4, pl_materials, 2,
                          sh_area1_shapes, sh_area1_materials));

    // random obj scene --------------------------
    printf("generating random shapes scenes ...\n");
    int rsh_nshapes = 32;
    yo_shape* rsh_shapes = make_random_shapes(rsh_nshapes, 5, false);
    yo_material* rsh_materials = make_random_materials(rsh_nshapes);

    // save plane scenes
    save_scene("rs01.obj", dirname,
               make_scene(2, sh_cams, 0, 0, rsh_nshapes, rsh_shapes,
                          rsh_nshapes, rsh_materials, 2, sh_light_shapes,
                          sh_light_materials));
    save_scene("rs02.obj", dirname,
               make_scene(2, sh_cams, 0, 0, rsh_nshapes, rsh_shapes,
                          rsh_nshapes, rsh_materials, 2, sh_area_shapes,
                          sh_area_materials));
    save_scene("rs03.obj", dirname,
               make_scene(2, sh_cams, 0, 0, rsh_nshapes, rsh_shapes,
                          rsh_nshapes, rsh_materials, 2, sh_area1_shapes,
                          sh_area1_materials));

    // xform obj scene ---------------------------
    printf("generating transform tests scenes ...\n");
    int xsh_nshapes = 32;
    yo_shape* xsh_shapes = make_random_shapes(xsh_nshapes, 5, true);
    yo_shape xsh_light_shapes[2] = {
        xform_shape(make_shape("light01", "light01", 0, ys_stype_points),
                    v(0.7, 4, 3), v1(0), v1(1)),
        xform_shape(make_shape("light02", "light02", 0, ys_stype_points),
                    v(-0.7, 4, 3), v1(0), v1(1))
    };
    yo_shape xsh_area_shapes[2] = {
        xform_shape(make_shape("light01", "light01", 0, ys_stype_uvquad),
                    v(4, 2, 4), v(0, 180 + 45, 0), v1(2)),
        xform_shape(make_shape("light02", "light02", 0, ys_stype_uvquad),
                    v(-4, 2, 4), v(0, 180 - 45, 0), v1(2))
    };

    // save plane scenes
    save_scene("xs01.obj", dirname,
               make_scene(2, sh_cams, 0, 0, xsh_nshapes, xsh_shapes,
                          rsh_nshapes, rsh_materials, 2, xsh_light_shapes,
                          sh_light_materials));
    save_scene("xs02.obj", dirname,
               make_scene(2, sh_cams, 0, 0, xsh_nshapes, xsh_shapes,
                          rsh_nshapes, rsh_materials, 2, xsh_area_shapes,
                          sh_area_materials));

    // env scene ------------------------------
    printf("generating envmaps scenes ...\n");
    yo_camera env_cams[2] = {
        make_camera("cam", v(0, 1.5, 5), v(0, 0.5, 0), 0.5, 0),
        make_camera("cam_dof", v(0, 1.5, 5), v(0, 0.5, 0), 0.5, 0.1)
    };
    yo_env env_envs[1] = { make_env("env", "env", v(0, 0.5, 0), v(0, 1.5, 0)) };
    yo_env env_txt_envs[1] = { make_env("env", "env_txt", v(0, 0.5, 0),
                                        v(-1.5, 0.5, 0)) };
    yo_shape env_shapes[4] = {
        transform_shape(make_floor("floor", "floor", 6, 4, 6), v(0, 0, -4),
                        v(0, 0, 0), v1(6)),
        transform_shape(
            make_shape("obj01", "obj01", 5, ys_stype_uvflipcapsphere),
            v(-1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(
            make_shape("obj02", "obj02", 4, ys_stype_uvspherizedcube),
            v(0, 0.5f, 0), v(0, 0, 0), v1(0.5)),
        transform_shape(make_shape("obj03", "obj03", 4, ys_stype_uvspherecube),
                        v(1.25f, 0.5f, 0), v(0, 0, 0), v1(0.5))
    };
    yo_material env_materials[6] = {
        make_diffuse("floor", v(0.2, 0.2, 0.2), 0),
        make_plastic("obj01", v(0.5, 0.2, 0.2), 50, 0),
        make_plastic("obj02", v(0.2, 0.5, 0.2), 100, 0),
        make_plastic("obj03", v(0.2, 0.2, 0.5), 500, 0),
        make_emission("env", v(1, 1, 1), 0),
        make_emission("env_txt", v(1, 1, 1), "env.hdr")
    };
    yo_shape env_light_shapes[1] = { transform_shape(
        make_shape("env_sphere", "env", 6, ys_stype_uvflippedsphere),
        v(0, 0.5, 0), v(-90, 0, 0), v1(10000)) };
    yo_material env_light_materials[1] = { make_emission("env", v(1, 1, 1),
                                                         0) };
    yo_shape env_txt_light_shapes[1] = { transform_shape(
        make_shape("env_sphere", "env_txt", 6, ys_stype_uvflippedsphere),
        v(0, 0.5, 0), v(-90, 0, 0), v1(10000)) };
    yo_material env_txt_light_materials[1] = { make_emission(
        "env_txt", v(1, 1, 1), "env.hdr") };

    // save plane scenes
    save_scene("es01.obj", dirname,
               make_scene(2, env_cams, 1, env_envs, 4, env_shapes, 6,
                          env_materials, 0, 0, 0));
    save_scene("es02.obj", dirname,
               make_scene(2, env_cams, 0, 0, 4, env_shapes, 6, env_materials, 1,
                          env_light_shapes, env_light_materials));
    save_scene("es03.obj", dirname,
               make_scene(2, env_cams, 1, env_txt_envs, 4, env_shapes, 6,
                          env_materials, 0, 0, 0));
    save_scene("es04.obj", dirname,
               make_scene(2, env_cams, 0, 0, 4, env_shapes, 6, env_materials, 1,
                          env_txt_light_shapes, env_txt_light_materials));

    // cornell box ------------------------------
    // http://graphics.cs.williams.edu/data
    // http://www.graphics.cornell.edu/online/box/data.html
    printf("generating cornell box scenes ...\n");
    yo_camera cb_cams[1] = { make_camera("cam", v(0, 1, 4), v(0, 1, 0), 0.7,
                                         0) };
    yo_shape cb_box[8] = {
        transform_shape(make_shape("floor", "white", 0, ys_stype_uvquad),
                        v(0, 0, 0), v(-90, 0, 0), v1(1)),
        transform_shape(make_shape("ceiling", "white", 0, ys_stype_uvquad),
                        v(0, 2, 0), v(90, 0, 0), v1(1)),
        transform_shape(make_shape("back", "white", 0, ys_stype_uvquad),
                        v(0, 1, -1), v(0, 0, 0), v1(1)),
        transform_shape(make_shape("back", "green", 0, ys_stype_uvquad),
                        v(+1, 1, 0), v(0, -90, 0), v1(1)),
        transform_shape(make_shape("back", "red", 0, ys_stype_uvquad),
                        v(-1, 1, 0), v(0, 90, 0), v1(1)),
        transform_shape(make_shape("tallbox", "white", 0, ys_stype_uvcube),
                        v(-0.33, 0.6, -0.29), v(0, 15, 0), v(0.3, 0.6, 0.3)),
        transform_shape(make_shape("shortbox", "white", 0, ys_stype_uvcube),
                        v(0.33, 0.3, 0.33), v(0, -15, 0), v(0.3, 0.3, 0.3)),
        transform_shape(make_shape("light", "light", 0, ys_stype_uvquad),
                        v(0, 1.999, 0), v(90, 0, 0), v(0.25, 0.25, 0.25))
    };
    yo_material cb_materials[4] = {
        make_diffuse("white", v(0.725, 0.71, 0.68), 0),
        make_diffuse("red", v(0.63, 0.065, 0.05), 0),
        make_diffuse("green", v(0.14, 0.45, 0.091), 0),
        make_emission("light", v(17, 12, 4), 0),
    };

    // save scene
    save_scene("cb01.obj", dirname, make_scene(1, cb_cams, 0, 0, 8, cb_box, 4,
                                               cb_materials, 0, 0, 0));

    // textures ---------------------------------
    printf("generating simple textures ...\n");
    save_image("grid.png", dirname, make_grid(512), 512);
    save_image("checker.png", dirname, make_checker(512), 512);
    save_image("rchecker.png", dirname, make_rchecker(512), 512);
    save_image("colored.png", dirname, make_colored(512), 512);
    save_image("rcolored.png", dirname, make_rcolored(512), 512);
    printf("generating envmaps textures ...\n");
    save_image_hdr("env.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8, 8,
                                   (ym_vec3f){ 0.2, 0.2, 0.2 }, 1 / powf(2, 6),
                                   true),
                   1024, 512);
    save_image_hdr("env01.hdr", dirname,
                   make_sunsky_hdr(1024, 512, 0.8, 8,
                                   (ym_vec3f){ 0.2, 0.2, 0.2 }, 1 / powf(2, 6),
                                   true),
                   1024, 512);
    save_image_hdr("env02.hdr", dirname,
                   make_sunsky2_hdr(1024, 512, 0.4, 8, 1 / powf(2, 12)), 1024,
                   512);
}
