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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
#include "CLI11.hpp"
using namespace std::literals;

struct app_state {
    std::string filename = "scene.json"s;
    std::string imfilename = "out.hdr"s;
    int resolution = 512;
    int nsamples = 256;
    std::string tracer = "pathtrace"s;
    ygl::trace_func tracef = &ygl::trace_path;
    int nbounces = 8;
    float pixel_clamp = 100.0f;
    bool noparallel = false;
    int seed = 7;
    int batch_size = 16;
    bool save_batch = false;
    float exposure = 0.0f;
    bool double_sided = false;
    bool add_skyenv = false;
    bool quiet = false;
};

auto tracer_names = std::unordered_map<std::string, ygl::trace_func>{
    {"pathtrace", ygl::trace_path},
    {"direct", ygl::trace_direct},
    {"environment", ygl::trace_environment},
    {"eyelight", ygl::trace_eyelight},
    {"pathtrace-nomis", ygl::trace_path_nomis},
    {"pathtrace-naive", ygl::trace_path_naive},
    {"direct-nomis", ygl::trace_direct_nomis},
    {"debug-normal", ygl::trace_debug_normal},
    {"debug-albedo", ygl::trace_debug_albedo},
    {"debug-texcoord", ygl::trace_debug_texcoord},
    {"debug-frontfacing", ygl::trace_debug_frontfacing},
};

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = std::make_shared<app_state>();

    // parse command line
    CLI::App parser("Offline path tracing", "ytrace");
    parser.add_option(
        "--resolution,-r", app->resolution, "Image vertical resolution.");
    parser.add_option("--nsamples,-s", app->nsamples, "Number of samples.");
    parser.add_option("--tracer,-t", app->tracer, "Trace type.")
        ->check([](const std::string& s) -> std::string {
            if (tracer_names.find(s) == tracer_names.end())
                throw CLI::ValidationError("unknown tracer name");
            return s;
        });
    parser.add_option(
        "--nbounces", app->nbounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", app->pixel_clamp, "Final pixel clamping.");
    parser.add_flag(
        "--noparallel", app->noparallel, "Disable parallel execution.");
    parser.add_option(
        "--seed", app->seed, "Seed for the random number generators.");
    parser.add_option("--batch-size", app->batch_size, "Sample batch size.");
    parser.add_flag(
        "--save-batch", app->save_batch, "Save images progressively");
    parser.add_option("--exposure,-e", app->exposure, "Hdr exposure");
    parser.add_flag(
        "--double-sided,-D", app->double_sided, "Double-sided rendering.");
    parser.add_flag("--add-skyenv,-E", app->add_skyenv, "add missing env map");
    parser.add_flag("--quiet,-q", app->quiet, "Print only errors messages");
    parser.add_option("--output-image,-o", app->imfilename, "Image filename");
    parser.add_option("scene", app->filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }
    app->tracef = tracer_names.at(app->tracer);

    // setup logger
    if (app->quiet) ygl::log_verbose() = false;

    // scene loading
    ygl::log_info_begin("loading scene {}", app->filename);
    auto scn = std::shared_ptr<ygl::scene>();
    try {
        scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        ygl::log_error("error during scene loading: "s + e.what());
        ygl::log_fatal("cannot load scene {}", app->filename);
    }
    ygl::log_info_end();

    // tesselate
    ygl::log_info("tesselating scene elements");
    ygl::update_tesselation(scn);

    // update bbox and transforms
    ygl::update_transforms(scn);
    ygl::update_bbox(scn);

    // add components
    ygl::log_info("adding scene elements");
    if (app->add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_environment("sky", {1, 1, 1},
            ygl::make_texture("sky", "sky.exr",
                ygl::make_sunsky_image(1024, 512, ygl::pi / 4))));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (app->double_sided) {
        for (auto mat : scn->materials) mat->double_sided = true;
    }
    if (scn->cameras.empty()) {
        scn->cameras.push_back(ygl::make_bbox_camera("<view>", scn->bbox));
    }
    auto cam = scn->cameras[0];
    ygl::add_missing_names(scn);
    for (auto err : ygl::validate(scn)) ygl::log_warning(err);

    // build bvh
    ygl::log_info_begin("building bvh");
    ygl::update_bvh(scn);
    ygl::log_info_end();

    // init renderer
    ygl::log_info("initializing lights");
    ygl::update_lights(scn);

    // initialize rendering objects
    ygl::log_info("initializing tracer data");
    auto img =
        ygl::image4f{(int)round(app->resolution * cam->width / cam->height),
            app->resolution};
    auto rngs = ygl::make_trace_rngs(img.width(), img.height(), app->seed);

    // render
    ygl::log_info_begin("rendering image");
    for (auto sample = 0; sample < app->nsamples; sample += app->batch_size) {
        if (app->save_batch && sample) {
            auto filename = ygl::replace_path_extension(app->imfilename, 
                std::to_string(sample) + "." + ygl::get_extension(app->imfilename));
            ygl::log_info("saving image {}", app->filename);
            ygl::save_image(filename, ygl::expose_image(img, app->exposure));
        }
        ygl::log_info_begin("rendering sample {}/{}", sample, app->nsamples);
        if (app->noparallel) {
            ygl::trace_samples(scn, cam, img, rngs, sample,
                std::min(app->batch_size, app->nsamples - sample), app->tracef,
                app->nbounces, app->pixel_clamp);
        } else {
            ygl::trace_samples_mt(scn, cam, img, rngs, sample,
                std::min(app->batch_size, app->nsamples - sample), app->tracef,
                app->nbounces, app->pixel_clamp);
        }
        ygl::log_info_end();
    }
    ygl::log_info_end();

    // save image
    ygl::log_info("saving image {}", app->imfilename);
    ygl::save_image(app->imfilename, ygl::expose_image(img, app->exposure));

    // done
    return 0;
}
