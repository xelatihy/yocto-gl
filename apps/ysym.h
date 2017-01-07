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

#ifndef _YSYM_APP_H_
#define _YSYM_APP_H_

#include "yapp.h"

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_sym.h"

namespace ysym_app {

inline ybvh::scene* make_bvh(const yapp::scene* scene) {
    auto scene_bvh = ybvh::make_scene((int)scene->shapes.size());
    auto sid = 0;
    for (auto shape : scene->shapes) {
        if (!shape->points.empty()) {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->points.size(),
                                  shape->points.data(), (int)shape->pos.size(),
                                  shape->pos.data(), shape->radius.data());
        } else if (!shape->lines.empty()) {
            ybvh::set_line_shape(scene_bvh, sid++, shape->frame,
                                 (int)shape->lines.size(), shape->lines.data(),
                                 (int)shape->pos.size(), shape->pos.data(),
                                 shape->radius.data());

        } else if (!shape->triangles.empty()) {
            ybvh::set_triangle_shape(
                scene_bvh, sid++, shape->frame, (int)shape->triangles.size(),
                shape->triangles.data(), (int)shape->pos.size(),
                shape->pos.data(), shape->radius.data());

        } else {
            ybvh::set_point_shape(scene_bvh, sid++, shape->frame,
                                  (int)shape->pos.size(), shape->pos.data(),
                                  shape->radius.data());
        }
    }
    ybvh::build_bvh(scene_bvh);
    return scene_bvh;
}

inline ysym::scene* make_rigid_scene(const yapp::scene* scene,
                                     ybvh::scene* scene_bvh) {
    // allocate scene
    auto rigid_scene = ysym::make_scene((int)scene->shapes.size());

    // add each shape
    auto sid = 0;
    for (auto shape : scene->shapes) {
        auto mat = scene->materials[shape->matid];
        auto density =
            (shape->name != "floor" && ym::length((ym::vec3f)mat->ke) == 0 &&
             !shape->triangles.empty())
                ? 1.0f
                : 0.0f;
        ysym::set_body(rigid_scene, sid++, shape->frame, {0, 0, 0}, {0, 0, 0},
                       density, (int)shape->triangles.size(),
                       shape->triangles.data(), (int)shape->pos.size(),
                       shape->pos.data());
    }

    // set up final bvh
    scene_bvh = make_bvh(scene);

    // setup collisions
    ysym::set_overlap_callbacks(
        rigid_scene, scene_bvh,
        [](auto ctx, std::vector<ysym::int2>* overlaps) {
            auto scene_bvh = (ybvh::scene*)ctx;
            ybvh::overlap_shape_bounds(scene_bvh, scene_bvh, false, true, true,
                                       overlaps);
        },
        [](auto ctx, int sid, const ysym::float3& pt, float max_dist) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto overlap =
                ybvh::overlap_point(scene_bvh, sid, pt, max_dist, false);
            return *(ysym::overlap_point*)&overlap;
        },
        [](auto ctx, int sid1, int sid2, float max_dist,
           std::vector<std::pair<ysym::overlap_point, ysym::int2>>* overlaps_) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto overlaps =
                (std::vector<std::pair<ybvh::point, ybvh::int2>>*)overlaps_;
            ybvh::overlap_verts(scene_bvh, scene_bvh, sid1, sid2, true,
                                max_dist, true, overlaps);
        },
        [](auto ctx, auto rigid_scene) {
            auto scene_bvh = (ybvh::scene*)ctx;
            for (auto sid = 0; sid < rigid_scene->shapes.size(); sid++) {
                ybvh::set_shape_frame(scene_bvh, sid,
                                      ysym::get_body_frame(rigid_scene, sid));
            }
            ybvh::refit_bvh(scene_bvh);
        });

    // initialize
    ysym::init_simulation(rigid_scene);

    return rigid_scene;
}

inline void simulate_step(yapp::scene* scene, ysym::scene* rigid_scene,
                          float dt) {
    ysym::advance_simulation(rigid_scene, dt);
    for (auto sid = 0; sid < scene->shapes.size(); sid++) {
        scene->shapes[sid]->frame = ysym::get_body_frame(rigid_scene, sid);
    }
}

struct params : virtual yapp::params {
    std::string outfilename;

    ybvh::scene* scene_bvh;
    ysym::scene* rigid_scene;

    float dt;
    int nframes;

    ~params() {
        if (scene_bvh) delete scene_bvh;
        if (rigid_scene) delete rigid_scene;
    }
};

inline void init_params(params* pars, ycmd::parser* parser) {
    pars->dt = ycmd::parse_opt<float>(parser, "--delta_time", "-dt",
                                      "delta time", 1 / 60.0f);
    pars->nframes = ycmd::parse_opt<int>(parser, "--nframes", "-n",
                                         "number of frames", 1000);
    pars->outfilename = ycmd::parse_opt<std::string>(
        parser, "--output", "-o", "output filename", "out.%04d.obj");

    // check up
    if (!pars->scene) yapp::init_params(pars, parser);
    if (!pars->scene) return;

    // init rigid simulation
    pars->scene_bvh = make_bvh(pars->scene);
    pars->rigid_scene = make_rigid_scene(pars->scene, pars->scene_bvh);
}

}  // namespace

#endif
