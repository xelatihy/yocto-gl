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

#ifndef _YTRACE_APP_H_
#define _YTRACE_APP_H_

#include "yapp.h"

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_cmd.h"
#include "../yocto/yocto_trace.h"

#include "ThreadPool.h"

namespace ytrace_app {

inline std::vector<ym::vec4i> make_image_blocks(int w, int h, int bs) {
    std::vector<ym::vec4i> blocks;
    for (int j = 0; j < h; j += bs) {
        for (int i = 0; i < w; i += bs) {
            blocks.push_back({i, j, ym::min(bs, w - i), ym::min(bs, h - j)});
        }
    }
    return blocks;
}

inline void save_image(const std::string& filename, int width, int height,
                       const ym::vec4f* hdr, float exposure, float gamma,
                       bool srgb_output) {
    auto ext = ycmd::get_extension(filename);
    auto tone_mapped = (const ym::vec4f*)nullptr;
    auto tone_mapped_buffer = std::vector<ym::vec4f>();
    if (exposure != 0 || gamma != 1) {
        tone_mapped_buffer =
            std::vector<ym::vec4f>(width * height, {0, 0, 0, 0});
        ym::exposure_gamma(width, height, 4, (const float*)hdr,
                           (float*)tone_mapped_buffer.data(), exposure, gamma,
                           false);
        tone_mapped = tone_mapped_buffer.data();
    } else {
        tone_mapped = hdr;
    }
    if (ext == ".hdr") {
        stbi_write_hdr(filename.c_str(), width, height, 4,
                       (float*)tone_mapped->data());
    } else if (ext == ".png") {
        auto ldr = std::vector<ym::vec4b>(width * height, {0, 0, 0, 0});
        if (srgb_output) {
            ym::linear_to_srgb(width, height, 4, (const float*)tone_mapped,
                               (unsigned char*)ldr.data());
        } else {
            ym::linear_to_byte(width, height, 4, (const float*)tone_mapped,
                               (unsigned char*)ldr.data());
        }
        stbi_write_png(filename.c_str(), width, height, 4, ldr.data(),
                       width * 4);
    } else {
        printf("supports only hdr and png for image writing\n");
        return;
    }
}

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

inline ytrace::scene* make_trace_scene(const yapp::scene* scene,
                                       const ybvh::scene* scene_bvh,
                                       int camera) {
    auto trace_scene = ytrace::make_scene(
        (int)scene->cameras.size(), (int)scene->shapes.size(),
        (int)scene->materials.size(), (int)scene->textures.size(),
        (int)scene->environments.size());

    auto cid = 0;
    for (auto cam : scene->cameras) {
        ytrace::set_camera(trace_scene, cid++, cam->frame, cam->yfov,
                           cam->aspect, cam->aperture, cam->focus);
    }

    auto tid = 0;
    for (auto txt : scene->textures) {
        if (!txt->hdr.empty()) {
            ytrace::set_texture(trace_scene, tid++, txt->width, txt->height,
                                txt->ncomp, txt->hdr.data());
        } else if (!txt->ldr.empty()) {
            ytrace::set_texture(trace_scene, tid++, txt->width, txt->height,
                                txt->ncomp, txt->ldr.data());
        } else
            assert(false);
    }

    auto eid = 0;
    for (auto env : scene->environments) {
        auto mat = scene->materials[env->matid];
        ytrace::set_environment(trace_scene, eid++, env->frame, mat->ke,
                                mat->ke_txt);
    }

    auto mid = 0;
    for (auto mat : scene->materials) {
        ytrace::set_material(trace_scene, mid++, mat->ke, mat->kd, mat->ks,
                             mat->rs, mat->ke_txt, mat->kd_txt, mat->ks_txt,
                             mat->ks_txt);
    }

    auto sid = 0;
    for (auto shape : scene->shapes) {
        if (!shape->points.empty()) {
            ytrace::set_point_shape(trace_scene, sid++, shape->frame,
                                    shape->matid, (int)shape->points.size(),
                                    shape->points.data(),
                                    (int)shape->pos.size(), shape->pos.data(),
                                    shape->norm.data(), shape->texcoord.data(),
                                    shape->color.data(), shape->radius.data());
        } else if (!shape->lines.empty()) {
            ytrace::set_line_shape(trace_scene, sid++, shape->frame,
                                   shape->matid, (int)shape->lines.size(),
                                   shape->lines.data(), (int)shape->pos.size(),
                                   shape->pos.data(), shape->norm.data(),
                                   shape->texcoord.data(), shape->color.data(),
                                   shape->radius.data());

        } else if (!shape->triangles.empty()) {
            ytrace::set_triangle_shape(
                trace_scene, sid++, shape->frame, shape->matid,
                (int)shape->triangles.size(), shape->triangles.data(),
                (int)shape->pos.size(), shape->pos.data(), shape->norm.data(),
                shape->texcoord.data(), shape->color.data());

        } else {
        }
    }

    ytrace::set_intersection_callbacks(
        trace_scene, (void*)scene_bvh,
        [](auto ctx, auto o, auto d, auto tmin, auto tmax) {
            auto scene_bvh = (ybvh::scene*)ctx;
            auto isec = ybvh::intersect_ray(scene_bvh, o, d, tmin, tmax, false);
            auto ipt = ytrace::intersect_point();
            ipt.dist = isec.dist;
            ipt.sid = isec.sid;
            ipt.eid = isec.eid;
            ipt.euv = {isec.euv[0], isec.euv[1], isec.euv[2]};
            return ipt;
        },
        [](auto ctx, auto o, auto d, auto tmin, auto tmax) {
            auto scene_bvh = (ybvh::scene*)ctx;
            return (bool)ybvh::intersect_ray(scene_bvh, o, d, tmin, tmax, true);
        });

    ytrace::init_lights(trace_scene);

    return trace_scene;
}

struct params : virtual yapp::params {
    ybvh::scene* scene_bvh;
    ytrace::scene* trace_scene;

    int width, height;
    std::vector<ym::vec4f> hdr;

    int samples;
    ytrace::render_params render_params;

    int block_size;
    std::vector<ym::vec4i> blocks;
    ThreadPool* pool;

    float exposure = 0;
    float gamma = 1;
    bool srgb = true;

    ~params() {
        if (scene_bvh) delete scene_bvh;
        if (trace_scene) delete trace_scene;
    }
};

inline void init_params(params* pars, ycmd::parser* parser) {
    auto rtype_names = std::unordered_map<std::string, ytrace::rng_type>{
        {"default", ytrace::rng_type::def},
        {"uniform", ytrace::rng_type::uniform},
        {"stratified", ytrace::rng_type::stratified},
        {"cmjs", ytrace::rng_type::cmjs}};
    auto stype_names = std::unordered_map<std::string, ytrace::shader_type>{
        {"default", ytrace::shader_type::def},
        {"eye", ytrace::shader_type::eyelight},
        {"direct", ytrace::shader_type::direct},
        {"path", ytrace::shader_type::pathtrace}};

    // params
    auto exposure = ycmd::parse_opt<float>(parser, "--exposure", "-e",
                                           "hdr image exposure", 0);
    auto gamma =
        ycmd::parse_opt<float>(parser, "--gamma", "-g", "hdr image gamma", 1);
    auto srgb =
        ycmd::parse_opt<bool>(parser, "--srgb", "", "hdr srgb output", true);
    auto rtype = ycmd::parse_opte<ytrace::rng_type>(
        parser, "--random", "", "random type", ytrace::rng_type::def,
        rtype_names);
    auto stype = ycmd::parse_opte<ytrace::shader_type>(
        parser, "--integrator", "-i", "integrator type",
        ytrace::shader_type::def, stype_names);
    auto camera_lights = ycmd::parse_flag(parser, "--camera_lights", "-c",
                                          "enable camera lights", false);
    auto camera = ycmd::parse_opt<int>(parser, "--camera", "-C", "camera", 0);
    auto nthreads = ycmd::parse_opt<int>(
        parser, "--threads", "-t", "number of threads [0 for default]", 0);
    auto block_size =
        ycmd::parse_opt<int>(parser, "--block_size", "", "block size", 32);
    auto samples =
        ycmd::parse_opt<int>(parser, "--samples", "-s", "image samples", 256);
    auto aspect = ycmd::parse_opt<float>(parser, "--aspect", "-a",
                                         "image aspect", 16.0f / 9.0f);
    auto res = ycmd::parse_opt<int>(parser, "--resolution", "-r",
                                    "image resolution", 720);

    // check up
    if (!pars->scene) yapp::init_params(pars, parser);
    if (!pars->scene) return;

    // setting up multithreading
    if (!nthreads) nthreads = std::thread::hardware_concurrency();
    pars->pool = new ThreadPool(nthreads);

    // fixing scene
    for (auto cam : pars->scene->cameras) cam->aspect = aspect;

    // preparing raytracer
    pars->scene_bvh = make_bvh(pars->scene);
    stype = (camera_lights) ? ytrace::shader_type::eyelight : stype;
    pars->trace_scene = make_trace_scene(pars->scene, pars->scene_bvh, camera);

    // rendering params
    pars->render_params = ytrace::render_params();
    pars->render_params.camera_id = camera;
    pars->render_params.nsamples = samples;
    pars->render_params.stype = stype;
    pars->render_params.rtype = rtype;
    pars->width = (int)std::round(aspect * res);
    pars->height = res;
    pars->samples = samples;
    pars->hdr =
        std::vector<ym::vec4f>(pars->width * pars->height, {0, 0, 0, 0});
    pars->block_size = block_size;
    pars->blocks = make_image_blocks(pars->width, pars->height, block_size);
    pars->exposure = exposure;
    pars->gamma = gamma;
    pars->srgb = srgb;
}

inline void render(params* pars) {
    std::vector<std::future<void>> futures;
    printf("tracing %s to %s\n", pars->filename.c_str(),
           pars->imfilename.c_str());
    printf("rendering ...");
    fflush(stdout);
    for (auto cur_sample = 0; cur_sample < pars->render_params.nsamples;
         cur_sample++) {
        printf("\rrendering sample %d/%d", cur_sample + 1,
               pars->render_params.nsamples);
        fflush(stdout);
        futures.clear();
        for (auto cur_block = 0; cur_block < pars->blocks.size(); cur_block++) {
            auto block = pars->blocks[cur_block];
            futures.push_back(pars->pool->enqueue([=]() {
                ytrace::trace_block(pars->trace_scene, pars->width,
                                    pars->height,
                                    (ytrace::float4*)pars->hdr.data(), block[0],
                                    block[1], block[2], block[3], cur_sample,
                                    cur_sample + 1, pars->render_params, true);
            }));
        }
        for (auto& future : futures) future.wait();
    }
    printf("\rrendering done\n");
    fflush(stdout);
    save_image(pars->imfilename, pars->width, pars->height, pars->hdr.data(),
               pars->exposure, pars->gamma, pars->srgb);
}

}  // namespace

#endif
