//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
using namespace std::literals;

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* view = nullptr;
    ygl::camera* cam = nullptr;
    ygl::bvh_tree* bvh = nullptr;
    std::string filename;
    std::string imfilename;
    ygl::load_options loadopts;
    ygl::image4f img;
    std::vector<ygl::trace_pixel> pixels;
    ygl::trace_params params;
    ygl::trace_lights lights;
    ygl::tonemap_type tonemapper = ygl::tonemap_type::gamma;
    float exposure = 0;
    bool save_batch = false;
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
    app->params.resolution = ygl::parse_opt(parser, "--resolution", "-r",
        "Image vertical resolution.", app->params.resolution);
    app->params.nsamples = ygl::parse_opt(
        parser, "--nsamples", "-s", "Number of samples.", app->params.nsamples);
    app->params.tracer = ygl::parse_opt(parser, "--tracer", "-T", "Trace type.",
        ygl::trace_type_names(), app->params.tracer);
    app->params.notransmission = ygl::parse_opt(parser, "--notransmission", "",
        "Whether to test transmission in shadows.", app->params.notransmission);
    app->params.double_sided = ygl::parse_opt(parser, "--double-sided", "-D",
        "Force double sided rendering.", app->params.double_sided);
    app->params.ambient = ygl::parse_opt(
        parser, "--ambient", "", "Ambient lighting.", app->params.ambient);
    app->params.envmap_invisible = ygl::parse_opt(parser, "--envmap-invisible",
        "", "View environment map.", app->params.envmap_invisible);
    app->params.min_depth = ygl::parse_opt(
        parser, "--min-depth", "", "Minimum ray depth.", app->params.min_depth);
    app->params.max_depth = ygl::parse_opt(
        parser, "--max-depth", "", "Maximum ray depth.", app->params.max_depth);
    app->params.pixel_clamp = ygl::parse_opt(parser, "--pixel-clamp", "",
        "Final pixel clamping.", app->params.pixel_clamp);
    app->params.ray_eps = ygl::parse_opt(parser, "--ray-eps", "",
        "Ray intersection epsilon.", app->params.ray_eps);
    app->params.parallel = ygl::parse_opt(
        parser, "--parallel", "", "Parallel execution.", app->params.parallel);
    app->params.seed = ygl::parse_opt(parser, "--seed", "",
        "Seed for the random number generators.", app->params.seed);
    app->params.preview_resolution = ygl::parse_opt(parser,
        "--preview-resolution", "", "Preview resolution for async rendering.",
        app->params.preview_resolution);
    app->params.batch_size = ygl::parse_opt(parser, "--batch-size", "",
        "Sample batch size.", app->params.batch_size);
    app->save_batch = ygl::parse_flag(
        parser, "--save-batch", "", "Save images progressively");
    app->tonemapper = ygl::parse_opt(parser, "--tonemapper", "t",
        "Tonemapper type.", ygl::tonemap_type_names(), app->tonemapper);
    app->exposure = ygl::parse_opt(
        parser, "--exposure", "t", "Hdr exposure", app->exposure);
    app->quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    auto filenames = ygl::parse_args(
        parser, "scenes", "Scene filenames", std::vector<std::string>());
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (app->quiet) ygl::get_default_logger()->verbose = false;

    // scene loading
    app->scn = new ygl::scene();
    for (auto filename : filenames) {
        try {
            ygl::log_info("loading scene {}", filename);
            auto scn = load_scene(filename, app->loadopts);
            ygl::merge_into(app->scn, scn);
            delete scn;
        } catch (std::exception e) {
            ygl::log_fatal("cannot load scene {}", filename);
        }
    }
    app->filename = filenames.front();

    // add elements
    ygl::add_names(app->scn);
    ygl::add_tangent_space(app->scn);

    // validate
    ygl::validate(app->scn, false, true);

    // view camera
    app->view = make_view_camera(app->scn, 0);
    app->cam = app->view;

    // build bvh
    ygl::log_info("building bvh");
    app->bvh = ygl::build_scene_bvh(app->scn);

    // init renderer
    ygl::log_info("initializing tracer");
    app->lights = make_trace_lights(app->scn);

    // initialize rendering objects
    app->img =
        ygl::make_image4f((int)round(app->cam->aspect * app->params.resolution),
            app->params.resolution);
    app->pixels = ygl::make_trace_pixels(app->img, app->params);

    // render
    ygl::log_info("starting renderer");
    ygl::trace_image(app->scn, app->cam, app->bvh, app->lights, app->img,
        app->pixels, app->params, [app](int cur_sample) {
            if (app->save_batch && cur_sample) {
                auto imfilename =
                    ygl::format("{}{}.{}{}", ygl::path_dirname(app->imfilename),
                        ygl::path_basename(app->imfilename), cur_sample,
                        ygl::path_extension(app->imfilename));
                ygl::log_info("saving image {}", imfilename);
                save_image(
                    imfilename, app->img, app->tonemapper, app->exposure);
            }
            ygl::log_info(
                "rendering sample {}/{}", cur_sample, app->params.nsamples);
        });
    ygl::log_info("rendering done");

    // save image
    ygl::log_info("saving image {}", app->imfilename);
    ygl::save_image(app->imfilename, app->img, app->tonemapper, app->exposure);

    // cleanup
    delete app;

    // done
    return 0;
}
