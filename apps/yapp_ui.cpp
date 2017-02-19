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

//
// YAPP OPENGL IMPLEMENTATION
//

#include "yapp.h"

#include "../yocto/yocto_glu.h"
#include "../yocto/yocto_math.h"

namespace yapp {

//
// Init shading
//
void init_shade_state(const yapp::scene* sc, yglu::uint& shade_prog,
                      yglu::uint& shade_vao,
                      std::map<texture*, yglu::uint>& shade_txt,
                      std::vector<std::array<yglu::uint, 7>>& shade_vbo) {
    yglu::stdshader::make_program(&shade_prog, &shade_vao);
    for (auto txt : sc->textures) {
        if (!txt->hdr.empty()) {
            shade_txt[txt] =
                yglu::modern::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->hdr.data(), true, true, true);
        } else if (!txt->ldr.empty()) {
            shade_txt[txt] =
                yglu::modern::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->ldr.data(), true, true, true);
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
void shade_scene(const yapp::scene* sc, int cur_camera, uint prog, uint vao,
                 const std::map<texture*, uint>& txts,
                 const std::vector<std::array<uint, 7>>& vbo,
                 const float4& background, float exposure, float gamma_,
                 bool srgb, bool wireframe, bool edges, bool camera_lights,
                 const float3& amb) {
    // begin frame
    yglu::enable_depth_test(true);
    yglu::enable_culling(false);
    yglu::clear_buffers();

    yglu::enable_wireframe(wireframe);

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
            if (!shape->mat) continue;
            auto mat = shape->mat;
            if (mat->ke == std::array<float, 3>{0, 0, 0}) continue;
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

        if (shape->mat) {
            auto mat = shape->mat;

#define __txt(txt) ((txt) ? txts.at(txt) : 0)
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

            yglu::line_width(2);
            yglu::enable_edges(true);
            yglu::stdshader::draw_points(prog, (int)shape->points.size(),
                                         bids[4]);
            yglu::stdshader::draw_lines(prog, (int)shape->lines.size(),
                                        bids[5]);
            yglu::stdshader::draw_triangles(prog, (int)shape->triangles.size(),
                                            bids[6]);
            yglu::enable_edges(false);
            yglu::line_width(1);
        }

        yglu::stdshader::end_shape();
    }

    yglu::stdshader::end_frame();

    yglu::enable_wireframe(false);
}

//
// Init shading
//
void init_draw_state(const yapp::scene* sc,
                     std::map<texture*, yglu::uint>& shade_txt) {
    for (auto txt : sc->textures) {
        if (!txt->hdr.empty()) {
            shade_txt[txt] =
                yglu::legacy::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->hdr.data(), true, true);
        } else if (!txt->ldr.empty()) {
            shade_txt[txt] =
                yglu::legacy::make_texture(txt->width, txt->height, txt->ncomp,
                                           txt->ldr.data(), true, true);
        } else
            assert(false);
    }
}

//
// Draw a sc
//
void draw_scene(const yapp::scene* sc, int cur_camera,
                const std::map<texture*, yglu::uint>& txts,
                const ym::vec4f& background, float exposure, float gamma_,
                bool wireframe, bool edges, bool camera_lights,
                const ym::vec3f& amb) {
    // begin frame
    yglu::enable_depth_test(true);
    yglu::enable_culling(false);
    yglu::clear_buffers();
    yglu::enable_wireframe(wireframe);

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
            if (!shape->mat) continue;
            auto mat = shape->mat;
            if (mat->ke == std::array<float, 3>{0, 0, 0}) continue;
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

        if (shape->mat) {
            auto mat = shape->mat;

#define __txt(txt) ((txt) ? txts.at(txt) : 0)
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

            yglu::line_width(2);
            yglu::enable_edges(true);
            yglu::legacy::draw_triangles((int)shape->triangles.size(),
                                         (yglu::int3*)shape->triangles.data(),
                                         (yglu::float3*)shape->pos.data(),
                                         (yglu::float3*)shape->norm.data(),
                                         (yglu::float2*)shape->texcoord.data(),
                                         (yglu::float3*)shape->color.data());
            yglu::line_width(1);
            yglu::enable_edges(false);
        }

        yglu::legacy::end_shape();
    }

    yglu::legacy::end_frame();

    yglu::enable_wireframe(false);
}

//
// OpenGL view
//
struct shade_state {
    // shade state
    uint shade_prog = 0;
    uint shade_vao = 0;
    std::map<texture*, uint> shade_txt;
    std::vector<std::array<uint, 7>> shade_vbo;
};

shade_state* init_shade_state(const yapp::scene* scn, const params* pars) {
    auto st = new shade_state();
    if (pars->legacy_gl) {
        init_draw_state(scn, st->shade_txt);
    } else {
        init_shade_state(scn, st->shade_prog, st->shade_vao, st->shade_txt,
                         st->shade_vbo);
    }
    return st;
}

//
// Cleanup state
//
void free_shade_state(shade_state* st) {
    if (st) delete st;
}

void shade_scene(const scene* scn, const params* pars, const shade_state* st) {
    if (pars->legacy_gl) {
        draw_scene(scn, pars->render_params.camera_id, st->shade_txt,
                   pars->background, pars->exposure, pars->gamma,
                   pars->wireframe, pars->edges,
                   pars->render_params.stype == ytrace::shader_type::eyelight,
                   pars->render_params.amb);
    } else {
        shade_scene(scn, pars->render_params.camera_id, st->shade_prog,
                    st->shade_vao, st->shade_txt, st->shade_vbo,
                    pars->background, pars->exposure, pars->gamma, pars->srgb,
                    pars->wireframe, pars->edges,
                    pars->render_params.stype == ytrace::shader_type::eyelight,
                    pars->render_params.amb);
    }
}

}  // namespace
