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
using namespace std::literals;

#include <map>

int main(int argc, char* argv[]) {
    static auto trace_names = std::map<ygl::trace_func, std::string>{
        {ygl::trace_path, "pathtrace"},
        {ygl::trace_direct, "direct"},
        {ygl::trace_environment, "environment"},
        {ygl::trace_eyelight, "eyelight"},
        {ygl::trace_path_nomis, "pathtrace-nomis"},
        {ygl::trace_path_naive, "pathtrace-naive"},
        {ygl::trace_direct_nomis, "direct-nomis"},
        {ygl::trace_debug_normal, "debug-normal"},
        {ygl::trace_debug_albedo, "debug-albedo"},
        {ygl::trace_debug_texcoord, "debug-texcoord"},
        {ygl::trace_debug_frontfacing, "debug-frontfacing"},
    };

    // parse command line
    auto parser =
        ygl::make_parser(argc, argv, "ytrace", "Offline oath tracing");
    auto resolution = ygl::parse_opt(
        parser, "--resolution", "-r", "Image vertical resolution.", 512);
    auto nsamples =
        ygl::parse_opt(parser, "--nsamples", "-s", "Number of samples.", 256);
    auto tracer = ygl::parse_opte(
        parser, "--tracer", "-t", "Trace type.", trace_names, &ygl::trace_path);
    auto max_depth = ygl::parse_opt(
        parser, "--nbounces", "", "Maximum number of bounces.", 8);
    auto pixel_clamp = ygl::parse_opt(
        parser, "--pixel-clamp", "", "Final pixel clamping.", 100.0f);
    auto noparallel = ygl::parse_flag(
        parser, "--noparallel", "", "Disable parallel execution.", false);
    auto seed = ygl::parse_opt(
        parser, "--seed", "", "Seed for the random number generators.", 7);
    auto batch_size =
        ygl::parse_opt(parser, "--batch-size", "", "Sample batch size.", 16);
    auto save_batch = ygl::parse_flag(
        parser, "--save-batch", "", "Save images progressively");
    auto exposure =
        ygl::parse_opt(parser, "--exposure", "-e", "Hdr exposure", 0.0f);
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided", "-D", "Double-sided rendering.", false);
    auto add_skyenv = ygl::parse_flag(
        parser, "--add-skyenv", "-E", "add missing env map", false);
    auto quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    auto imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    auto filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (quiet) ygl::log_verbose() = false;

    // scene loading
    ygl::log_info_begin("loading scene {}", filename);
    auto scn = std::shared_ptr<ygl::scene>();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        ygl::log_error("error during scene loading: "s + e.what());
        ygl::log_fatal("cannot load scene {}", filename);
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
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_environment("sky", {1, 1, 1},
            ygl::make_texture("sky", "sky.exr",
                ygl::make_sunsky_image(1024, 512, ygl::pi / 4))));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided) {
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
    auto img = ygl::make_image4f(
        (int)round(resolution * cam->width / cam->height), resolution);
    auto rngs = ygl::make_rng_seq(img.width * img.height, seed);

    // render
    ygl::log_info_begin("rendering image");
    for (auto sample = 0; sample < nsamples; sample += batch_size) {
        if (save_batch && sample) {
            auto filename = ygl::format("{}{}.{}{}",
                ygl::path_dirname(imfilename), ygl::path_basename(imfilename),
                sample, ygl::path_extension(imfilename));
            ygl::log_info("saving image {}", filename);
            ygl::save_image(filename, ygl::expose_image(img, exposure));
        }
        ygl::log_info_begin("rendering sample {}/{}", sample, nsamples);
        if (noparallel) {
            ygl::trace_samples(scn, cam, img, rngs, sample,
                std::min(batch_size, nsamples - sample), tracer, max_depth,
                pixel_clamp);
        } else {
            ygl::trace_samples_mt(scn, cam, img, rngs, sample,
                std::min(batch_size, nsamples - sample), tracer, max_depth,
                pixel_clamp);
        }
        ygl::log_info_end();
    }
    ygl::log_info_end();

    // save image
    ygl::log_info("saving image {}", imfilename);
    ygl::save_image(imfilename, ygl::expose_image(img, exposure));

    // done
    return 0;
}
