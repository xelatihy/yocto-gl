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

#include "../yocto/yocto_gl.h"
using namespace ygl;

// Application state
struct app_state {
    // scene data
    scene* scn = nullptr;

    // filenames
    string filename;
    string imfilename;

    // render
    int resolution = 0;
    float exposure = 0, gamma = 2.2f;
    bool filmic = false;
    vec4f background = {0, 0, 0, 0};

    // trace
    trace_params trace_params_;
    bool trace_save_progressive = false;
    int trace_block_size = 32;
    int trace_batch_size = 16;

    // rendered images and buffers
    image4f trace_img;
    image4f trace_acc;
    image4f trace_weight;
    vector<rng_pcg32> trace_rngs;

    ~app_state() {
        if (scn) delete scn;
    }
};

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = make_parser(argc, argv, "ytrace", "offline oath tracing");
    app->trace_params_.camera_id = 0;
    app->trace_save_progressive =
        parse_flag(parser, "--save-progressive", "", "save progressive images");
    app->trace_params_.rtype = parse_opt(parser, "--random", "", "random type",
        trace_rng_names(), trace_rng_type::stratified);
    app->trace_params_.ftype = parse_opt(parser, "--filter", "", "filter type",
        trace_filter_names(), trace_filter_type::box);
    app->trace_params_.stype =
        parse_opt(parser, "--shader", "-S", "path estimator type",
            trace_shader_names(), trace_shader_type::pathtrace);
    app->trace_params_.envmap_invisible =
        parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
    app->trace_params_.shadow_notransmission = parse_flag(
        parser, "--shadow-notransmission", "", "shadow without transmission");
    app->trace_block_size =
        parse_opt(parser, "--block-size", "", "block size", 32);
    app->trace_batch_size =
        parse_opt(parser, "--batch-size", "", "batch size", 16);
    app->trace_params_.nsamples =
        parse_opt(parser, "--samples", "-s", "image samples", 256);
    app->trace_params_.parallel =
        !parse_flag(parser, "--no-parallel", "", "so not run in parallel");
    app->exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->gamma = parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->filmic = parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    app->resolution =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    auto camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
    app->trace_params_.amb = {amb, amb, amb};
    if (camera_lights) {
        app->trace_params_.stype = trace_shader_type::eyelight;
    }
    auto log_filename = parse_opt(parser, "--log", "", "log to disk", ""s);
    if (log_filename != "") add_file_stream(log_filename, true);
    app->imfilename =
        parse_opt(parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = parse_arg(parser, "scene", "scene filename", ""s);
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setting up rendering
    log_info("loading scene {}", app->filename);
    try {
        app->scn = load_scene(app->filename);
    } catch (exception e) {
        log_fatal("cannot load scene {}", app->filename);
        return 1;
    }
    auto opts = add_elements_options();
    opts.pointline_radius = 0.001f;
    add_elements(app->scn, opts);

    // build bvh
    log_info("building bvh");
    build_bvh(app->scn);

    // init renderer
    log_info("initializing tracer");
    update_lights(app->scn, false);

    // initialize rendering objects
    auto width =
        (int)round(app->scn->cameras[app->trace_params_.camera_id]->aspect *
                   app->resolution);
    auto height = app->resolution;
    app->trace_params_.width = width;
    app->trace_params_.height = height;
    app->trace_img = image4f(width, height);
    app->trace_acc = image4f(width, height);
    app->trace_weight = image4f(width, height);
    app->trace_rngs = trace_rngs(app->trace_params_);

    // render
    log_info("starting renderer");
    for (auto cur_sample = 0; cur_sample < app->trace_params_.nsamples;
         cur_sample += app->trace_batch_size) {
        if (app->trace_save_progressive && cur_sample) {
            auto imfilename = format("{}{}.{}{}", path_dirname(app->imfilename),
                path_basename(app->imfilename), cur_sample,
                path_extension(app->imfilename));
            log_info("saving image {}", imfilename);
            save_image(imfilename, app->trace_img, app->exposure, app->gamma,
                app->filmic);
        }
        log_info(
            "rendering sample {}/{}", cur_sample, app->trace_params_.nsamples);
        if (app->trace_params_.ftype == trace_filter_type::box) {
            trace_samples(app->scn, app->trace_img, cur_sample,
                min(cur_sample + app->trace_batch_size,
                    app->trace_params_.nsamples),
                app->trace_rngs, app->trace_params_);
        } else {
            trace_filtered_samples(app->scn, app->trace_img, app->trace_acc,
                app->trace_weight, cur_sample,
                min(cur_sample + app->trace_batch_size,
                    app->trace_params_.nsamples),
                app->trace_rngs, app->trace_params_);
        }
    }
    log_info("rendering done");

    // save image
    log_info("saving image {}", app->imfilename);
    save_image(app->imfilename, app->trace_img, app->exposure, app->gamma,
        app->filmic);

    // cleanup
    delete app;

    // done
    return 0;
}
