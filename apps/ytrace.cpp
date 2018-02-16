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
    bool save_batch = false;
    int batch_size = 16;
    bool quiet = false;

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
        ygl::make_parser(argc, argv, "ytrace", "Offline oath tracing");
    app->params = ygl::parse_params(parser, "", app->params);
    app->batch_size = ygl::parse_opt(parser, "--batch-size", "",
        "Compute images in <val> samples batches", 16);
    app->save_batch = ygl::parse_flag(
        parser, "--save-batch", "", "Save images progressively");
    app->quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (app->quiet) ygl::get_default_logger()->verbose = false;

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
    app->img =
        ygl::image4f((int)round(app->cam->aspect * app->params.resolution),
            app->params.resolution);
    app->pixels = ygl::make_trace_pixels(app->img, app->params);

    // render
    ygl::log_info("starting renderer");
    for (auto cur_sample = 0; cur_sample < app->params.nsamples;
         cur_sample += app->batch_size) {
        if (app->save_batch && cur_sample) {
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
            app->pixels, app->batch_size, app->params);
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
