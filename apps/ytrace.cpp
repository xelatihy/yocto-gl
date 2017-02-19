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

#include "yapp.h"

int main(int argc, char* argv[]) {
    // logging
    yapp::set_default_loggers();

    // params
    auto pars = yapp::init_params("render scene with path tracing", argc, argv,
                                  true, false, false, false);

    // setting up rendering
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "loading scene %s",
                   pars->filenames[0].c_str());
    auto scene = yapp::load_scenes(pars->filenames, pars->scene_scale);
    auto scene_bvh = yapp::make_bvh(scene);
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "setting up tracer");
    auto trace_scene =
        yapp::make_trace_scene(scene, scene_bvh, pars->render_params.camera_id);

    // fixing camera
    for (auto cam : scene->cameras)
        cam->aspect = (float)pars->width / (float)pars->height;

    // init renderer
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "initializing tracer");
    ytrace::init_lights(trace_scene);

    // allocate image
    auto hdr = new std::array<float, 4>[pars->width * pars->height];
    for (auto i = 0; i < pars->width * pars->height; i++) hdr[i] = {0, 0, 0, 0};

    // render
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "starting renderer");
    auto blocks =
        yapp::make_trace_blocks(pars->width, pars->height, pars->block_size);
    auto pool = ycmd::make_thread_pool(pars->nthreads);
    std::vector<std::future<void>> futures;
    for (auto cur_sample = 0; cur_sample < pars->render_params.nsamples;
         cur_sample++) {
        ycmd::log_msgf(ycmd::log_level_info, "ytrace",
                       "rendering sample %4d/%d", cur_sample + 1,
                       pars->render_params.nsamples);
        futures.clear();
        for (auto cur_block = 0; cur_block < blocks.size(); cur_block++) {
            auto block = blocks[cur_block];
            futures.push_back(ycmd::pool_enqueue(pool, [=]() {
                ytrace::trace_block(trace_scene, pars->width, pars->height,
                                    (ytrace::float4*)hdr, block[0], block[1],
                                    block[2], block[3], cur_sample,
                                    cur_sample + 1, pars->render_params);
            }));
        }
        for (auto& future : futures) future.wait();
        if (pars->save_progressive &&
            (cur_sample + 1) % pars->save_progressive == 0) {
            auto imfilename = ycmd::get_dirname(pars->imfilename) +
                              ycmd::get_basename(pars->imfilename) +
                              ycmd::format_str(".%04d", cur_sample + 1) +
                              ycmd::get_extension(pars->imfilename);
            ycmd::log_msgf(ycmd::log_level_info, "ytrace", "saving image %s",
                           imfilename.c_str());
            yapp::save_image(imfilename, pars->width, pars->height, hdr,
                             pars->exposure, pars->gamma, pars->srgb);
        }
    }
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "rendering done");

    // save image
    ycmd::log_msgf(ycmd::log_level_info, "ytrace", "saving image %s",
                   pars->imfilename.c_str());
    yapp::save_image(pars->imfilename, pars->width, pars->height, hdr,
                     pars->exposure, pars->gamma, pars->srgb);

    // done
    // cleanup
    delete scene;
    ybvh::free_scene(scene_bvh);
    ytrace::free_scene(trace_scene);
    delete[] hdr;
    return EXIT_SUCCESS;
}
