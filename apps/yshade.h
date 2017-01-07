//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YSHADE_APP_H_
#define _YSHADE_APP_H_

// clang-format off
#ifndef __APPLE__
#include <GL/glew.h>
#else
#include <OpenGL/gl.h>
#include <OpenGL/gl3.h>
#endif
#include <GLFW/glfw3.h>
// clang-format on

#include "yapp.h"
#include "../yocto/yocto_glu.h"

namespace yshade_app {

//
// Init shading
//
inline void init_shade(const yapp::scene* sc, yglu::uint& shade_prog,
                       yglu::uint& shade_vao,
                       std::vector<yglu::uint>& shade_txt,
                       std::vector<std::array<yglu::uint, 7>>& shade_vbo) {
    yglu::stdshader::make_program(&shade_prog, &shade_vao);
    for (auto txt : sc->textures) {
        if (!txt->hdr.empty()) {
            shade_txt.push_back(
                yglu::modern::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->hdr.data(), true, true, true));
        } else if (!txt->ldr.empty()) {
            shade_txt.push_back(
                yglu::modern::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->ldr.data(), true, true, true));
        } else
            assert(false);
    }
    for (auto shape : sc->shapes) {
        shade_vbo.push_back({0, 0, 0, 0});
        auto& vbos = shade_vbo.back();
        if (!shape->pos.empty())
            vbos[0] = yglu::modern::make_buffer(
                (int)shape->pos.size(), 3 * sizeof(float), shape->pos.data(),
                false, false);
        if (!shape->norm.empty())
            vbos[1] = yglu::modern::make_buffer(
                (int)shape->norm.size(), 3 * sizeof(float), shape->norm.data(),
                false, false);
        if (!shape->texcoord.empty())
            vbos[2] = yglu::modern::make_buffer(
                (int)shape->texcoord.size(), 2 * sizeof(float),
                shape->texcoord.data(), false, false);
        if (!shape->color.empty())
            vbos[3] = yglu::modern::make_buffer(
                (int)shape->color.size(), 3 * sizeof(float),
                shape->color.data(), false, false);
        if (!shape->points.empty())
            vbos[4] = yglu::modern::make_buffer(
                (int)shape->points.size(), sizeof(int), shape->points.data(),
                true, false);
        if (!shape->lines.empty())
            vbos[5] = yglu::modern::make_buffer(
                (int)shape->lines.size(), 2 * sizeof(int), shape->lines.data(),
                true, false);
        if (!shape->triangles.empty())
            vbos[6] = yglu::modern::make_buffer(
                (int)shape->triangles.size(), 3 * sizeof(int),
                shape->triangles.data(), true, false);
    }
}

//
// Display a scene
//
inline void shade(const yapp::scene* sc, int cur_camera, yglu::uint prog,
                  yglu::uint vao, const std::vector<yglu::uint>& txt,
                  std::vector<std::array<yglu::uint, 7>>& vbo,
                  const ym::vec4f& background, float exposure, float gamma_,
                  bool srgb, bool wireframe, bool edges, bool camera_lights,
                  const ym::vec3f& amb) {
    // begin frame
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    auto cam = sc->cameras[cur_camera];
    auto camera_xform = ym::to_mat((ym::frame3f)cam->frame);
    auto camera_view = ym::to_mat(ym::inverse((ym::frame3f)cam->frame));
    auto camera_proj =
        ym::perspective_mat4(cam->yfov, cam->aspect, 0.1f, 100000.0f);

    yglu::stdshader::begin_frame(prog, vao, camera_lights, exposure, gamma_,
                                 srgb, camera_xform, camera_view, camera_proj);

    if (!camera_lights) {
        auto nlights = 0;
        std::array<ym::vec3f, 16> light_pos, light_ke;
        std::array<yglu::ltype, 16> light_type;
        for (auto shape : sc->shapes) {
            if (shape->matid < 0) continue;
            auto mat = sc->materials[shape->matid];
            if (mat->ke == ym::zero3f) continue;
            for (auto p : shape->points) {
                if (nlights >= 16) continue;
                light_pos[nlights] = shape->pos[p];
                light_pos[nlights] = ym::transform_point(
                    (ym::frame3f)shape->frame, light_pos[nlights]);
                light_ke[nlights] = mat->ke;
                light_type[nlights] = yglu::ltype::point;
                nlights++;
            }
        }
        yglu::stdshader::set_lights(
            prog, amb, nlights, (yglu::float3*)light_pos.data(),
            (yglu::float3*)light_ke.data(), light_type.data());
    }

    for (auto sid = 0; sid < sc->shapes.size(); sid++) {
        auto shape = sc->shapes[sid];
        yglu::stdshader::begin_shape(prog,
                                     ym::to_mat((ym::frame3f)shape->frame));

        if (shape->matid >= 0) {
            auto mat = sc->materials[shape->matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yglu::stdshader::set_material(
                prog, mat->ke, mat->kd, mat->ks, mat->rs, __txt(mat->ke_txt),
                __txt(mat->kd_txt), __txt(mat->ks_txt), __txt(mat->rs_txt),
                false);
        } else {
            auto kd = ym::vec3f{0.8f, 0.8f, 0.8f};
            yglu::stdshader::set_material(prog, {0, 0, 0}, kd, {0, 0, 0}, 0, 0,
                                          0, 0, 0, false);
        }

        auto bids = vbo[sid];
        yglu::stdshader::set_vert(prog, bids[0], bids[1], bids[2], bids[3]);

        yglu::stdshader::draw_points(prog, (int)shape->points.size(), bids[4]);
        yglu::stdshader::draw_lines(prog, (int)shape->lines.size(), bids[5]);
        yglu::stdshader::draw_triangles(prog, (int)shape->triangles.size(),
                                        bids[6]);

        if (edges && !wireframe) {
            yglu::stdshader::set_material(prog, {0, 0, 0}, {0, 0, 0}, {0, 0, 0},
                                          0, 0, 0, 0, 0, false);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yglu::stdshader::draw_triangles(prog, (int)shape->triangles.size(),
                                            bids[6]);
            glDepthRange(0, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glLineWidth(1);
        }

        yglu::stdshader::end_shape();
    }

    yglu::stdshader::end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

//
// Init shading
//
inline void init_draw(const yapp::scene* sc,
                      std::vector<yglu::uint>& shade_txt) {
    for (auto txt : sc->textures) {
        if (!txt->hdr.empty()) {
            shade_txt.push_back(
                yglu::legacy::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->hdr.data(), true, true));
        } else if (!txt->ldr.empty()) {
            shade_txt.push_back(
                yglu::legacy::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->ldr.data(), true, true));
        } else
            assert(false);
    }
}

//
// Draw a sc
//
inline void draw(const yapp::scene* sc, int cur_camera,
                 const std::vector<yglu::uint>& txt,
                 const ym::vec4f& background, float exposure, float gamma_,
                 bool wireframe, bool edges, bool camera_lights,
                 const ym::vec3f& amb) {
    // begin frame
    glEnable(GL_DEPTH_TEST);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glDisable(GL_CULL_FACE);

    if (wireframe) glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    auto cam = sc->cameras[cur_camera];
    auto camera_xform = ym::to_mat((ym::frame3f)cam->frame);
    auto camera_view = ym::to_mat(ym::inverse((ym::frame3f)cam->frame));
    auto camera_proj =
        ym::perspective_mat4(cam->yfov, cam->aspect, 0.1f, 100000.0f);

    auto nlights = 0;
    std::array<ym::vec3f, 16> light_pos, light_ke;
    std::array<yglu::ltype, 16> light_type;

    yglu::legacy::begin_frame(camera_xform, camera_view, camera_proj,
                              camera_lights, true);

    if (!camera_lights) {
        for (auto shape : sc->shapes) {
            if (shape->matid < 0) continue;
            auto mat = sc->materials[shape->matid];
            if (mat->ke == ym::zero3f) continue;
            for (auto p : shape->points) {
                if (nlights >= 16) continue;
                light_pos[nlights] = shape->pos[p];
                light_pos[nlights] = ym::transform_point(
                    (ym::frame3f)shape->frame, light_pos[nlights]);
                light_ke[nlights] = mat->ke;
                light_type[nlights] = yglu::ltype::point;
                nlights++;
            }
        }
        yglu::legacy::set_lights(amb, nlights, (yglu::float3*)light_pos.data(),
                                 (yglu::float3*)light_ke.data(),
                                 light_type.data());
    }

    for (auto shape : sc->shapes) {
        yglu::legacy::begin_shape(ym::to_mat((ym::frame3f)shape->frame));

        if (shape->matid >= 0) {
            auto mat = sc->materials[shape->matid];

#define __txt(i) ((i >= 0) ? txt[i] : 0)
            yglu::legacy::set_material(
                mat->ke, mat->kd, mat->ks,
                yglu::legacy::specular_roughness_to_exponent(mat->rs),
                __txt(mat->kd_txt), true);
        } else {
            auto kd = ym::vec3f{0.8f, 0.8f, 0.8f};
            yglu::legacy::set_material({0, 0, 0}, kd, {0, 0, 0}, 0, 0, true);
        }

        yglu::legacy::draw_points(
            (int)shape->points.size(), shape->points.data(),
            (yglu::float3*)shape->pos.data(), (yglu::float3*)shape->norm.data(),
            (yglu::float2*)shape->texcoord.data(),
            (yglu::float3*)shape->color.data());
        yglu::legacy::draw_lines(
            (int)shape->lines.size(), (yglu::int2*)shape->lines.data(),
            (yglu::float3*)shape->pos.data(), (yglu::float3*)shape->norm.data(),
            (yglu::float2*)shape->texcoord.data(),
            (yglu::float3*)shape->color.data());
        yglu::legacy::draw_triangles(
            (int)shape->triangles.size(), (yglu::int3*)shape->triangles.data(),
            (yglu::float3*)shape->pos.data(), (yglu::float3*)shape->norm.data(),
            (yglu::float2*)shape->texcoord.data(),
            (yglu::float3*)shape->color.data());

        if (edges && !wireframe) {
            yglu::legacy::set_material({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 0, 0,
                                       true);

            glLineWidth(2);
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
            glDepthRange(0, 0.999999);
            yglu::legacy::draw_triangles((int)shape->triangles.size(),
                                         (yglu::int3*)shape->triangles.data(),
                                         (yglu::float3*)shape->pos.data(),
                                         (yglu::float3*)shape->norm.data(),
                                         (yglu::float2*)shape->texcoord.data(),
                                         (yglu::float3*)shape->color.data());
            glDepthRange(0, 1);
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glLineWidth(1);
        }

        yglu::legacy::end_shape();
    }

    yglu::legacy::end_frame();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

//
// Shade state
//
struct params : virtual yapp::params {
    bool legacy_gl, no_ui;

    int width, height;

    int camera_id;
    float background, amb;
    float exposure, gamma;
    bool srgb;

    bool wireframe, edges;
    bool camera_lights;

    yglu::uint shade_prog = 0;
    yglu::uint shade_vao = 0;
    std::vector<yglu::uint> shade_txt;
    std::vector<std::array<yglu::uint, 7>> shade_vbo;

    void* widget_ctx = nullptr;
};

inline void init_params(params* pars, ycmd::parser* parser) {
    // parse cmdline
    pars->exposure =
        ycmd::parse_opt<float>(parser, "--exposure", "-e", "image exposure", 0);
    pars->gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "image gamma", 1);
    pars->srgb =
        ycmd::parse_opt<bool>(parser, "--srgb", "", "image srgb", true);
    pars->amb =
        ycmd::parse_opt<float>(parser, "--ambient", "", "ambient factor", 0);
    pars->camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                           "enable camera lights", false);
    pars->camera_id =
        ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    pars->no_ui =
        ycmd::parse_flag(parser, "--no-ui", "", "runs offline", false);
    pars->legacy_gl = ycmd::parse_flag(parser, "--legacy_opengl", "-L",
                                       "uses legacy OpenGL", false);
    auto aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a",
                                         "image aspect", 16.0f / 9.0f);
    auto res = ycmd::parse_opt<int>(parser, "--resolution", "-r",
                                    "image resolution", 720);

    // check up
    if (!pars->scene) yapp::init_params(pars, parser);
    if (!pars->scene) return;

    // fixing scene
    for (auto cam : pars->scene->cameras) cam->aspect = aspect;

    // rendering params
    pars->width = (int)std::round(aspect * res);
    pars->height = res;
}

inline void init(params* pars) {
    if (pars->legacy_gl) {
        init_draw(pars->scene, pars->shade_txt);
    } else {
        init_shade(pars->scene, pars->shade_prog, pars->shade_vao,
                   pars->shade_txt, pars->shade_vbo);
    }
}

inline void render(params* pars) {
    if (pars->legacy_gl) {
        draw(pars->scene, pars->camera_id, pars->shade_txt,
             {pars->background, pars->background, pars->background, 0},
             pars->exposure, pars->gamma, pars->wireframe, pars->edges,
             pars->camera_lights, {pars->amb, pars->amb, pars->amb});
    } else {
        shade(pars->scene, pars->camera_id, pars->shade_prog, pars->shade_vao,
              pars->shade_txt, pars->shade_vbo,
              {pars->background, pars->background, pars->background, 0},
              pars->exposure, pars->gamma, pars->srgb, pars->wireframe,
              pars->edges, pars->camera_lights,
              {pars->amb, pars->amb, pars->amb});
    }
}

}  // namespace

#endif
