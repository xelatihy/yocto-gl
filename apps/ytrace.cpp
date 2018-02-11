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
using namespace std::literals;

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* view = nullptr;
    ygl::camera* cam = nullptr;
    ygl::bvh_tree* bvh = nullptr;
    std::string filename;
    std::string imfilename;
    ygl::image4f img;
    ygl::image<ygl::trace_pixel> pixels;
    ygl::trace_params params;
    ygl::trace_lights lights;
    float exposure = 0, gamma = 2.2f;
    bool filmic = false;
    ygl::vec4f background = {0, 0, 0, 0};
    bool save_progressive = false;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
        if (bvh) delete bvh;
    }
};

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser =
        ygl::make_parser(argc, argv, "ytrace", "offline oath tracing");
    app->save_progressive = ygl::parse_flag(
        parser, "--save-progressive", "", "save progressive images");
    app->params.rtype = ygl::parse_opt(parser, "--random", "", "random type",
        ygl::enum_names<ygl::trace_rng_type>(),
        ygl::trace_rng_type::stratified);
    app->params.ftype = ygl::parse_opt(parser, "--filter", "", "filter type",
        ygl::enum_names<ygl::trace_filter_type>(), ygl::trace_filter_type::box);
    app->params.stype = ygl::parse_opt(parser, "--shader", "-S",
        "path estimator type", ygl::enum_names<ygl::trace_shader_type>(),
        ygl::trace_shader_type::pathtrace);
    app->params.envmap_invisible =
        ygl::parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
    app->params.shadow_notransmission = ygl::parse_flag(
        parser, "--shadow-notransmission", "", "shadow without transmission");
    app->params.block_size =
        ygl::parse_opt(parser, "--block-size", "", "block size", 32);
    app->params.batch_size =
        ygl::parse_opt(parser, "--batch-size", "", "batch size", 16);
    app->params.nsamples =
        ygl::parse_opt(parser, "--samples", "-s", "image samples", 256);
    app->params.parallel =
        !ygl::parse_flag(parser, "--no-parallel", "", "so not run in parallel");
    app->exposure =
        ygl::parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->gamma =
        ygl::parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->filmic =
        ygl::parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    app->params.height =
        ygl::parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = ygl::parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    auto camera_lights = ygl::parse_flag(
        parser, "--camera-lights", "-c", "enable camera lights");
    app->params.ambient = {amb, amb, amb};
    if (camera_lights) { app->params.stype = ygl::trace_shader_type::eyelight; }
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setting up rendering
    ygl::log_info("loading scene {}", app->filename);
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
        return 1;
    }

    // add elements
    auto opts = ygl::add_elements_options();
    add_elements(app->scn, opts);

    // view camera
    app->view = make_view_camera(app->scn, 0);
    app->cam = app->view;

    // build bvh
    ygl::log_info("building bvh");
    app->bvh = make_bvh(app->scn);

    // init renderer
    ygl::log_info("initializing tracer");
    app->lights = make_trace_lights(app->scn);

    // initialize rendering objects
    app->params.width = (int)round(app->cam->aspect * app->params.height);
    app->img = ygl::image4f(app->params.width, app->params.height);
    app->pixels = make_trace_pixels(app->params);

    // render
    ygl::log_info("starting renderer");
    for (auto cur_sample = 0; cur_sample < app->params.nsamples;
         cur_sample += app->params.batch_size) {
        if (app->save_progressive && cur_sample) {
            auto imfilename =
                ygl::format("{}{}.{}{}", ygl::path_dirname(app->imfilename),
                    ygl::path_basename(app->imfilename), cur_sample,
                    ygl::path_extension(app->imfilename));
            ygl::log_info("saving image {}", imfilename);
            save_image(
                imfilename, app->img, app->exposure, app->gamma, app->filmic);
        }
        ygl::log_info(
            "rendering sample {}/{}", cur_sample, app->params.nsamples);
        trace_samples(app->scn, app->cam, app->bvh, app->lights, app->img,
            app->pixels, app->params.batch_size, app->params);
    }
    ygl::log_info("rendering done");

    // save image
    ygl::log_info("saving image {}", app->imfilename);
    ygl::save_image(
        app->imfilename, app->img, app->exposure, app->gamma, app->filmic);

    // cleanup
    delete app;

    // done
    return 0;
}
