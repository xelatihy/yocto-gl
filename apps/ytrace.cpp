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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"
#include "CLI11.hpp"
using namespace std::literals;

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
    // command line parameters
    auto filename = "scene.json"s;  // scene filename
    auto imfilename = "out.hdr"s;   // image output filename
    auto camid = 0;                 // camera index
    auto resolution = 512;          // image vertical resolution
    auto nsamples = 256;            // image samples
    auto tracer = "pathtrace"s;     // tracer algorithm
    auto nbounces = 8;              // number of bounces
    auto pixel_clamp = 100.0f;      // pixel clamping
    auto noparallel = false;        // disable parallel
    auto seed = 7;                  // random seed
    auto nbatch = 16;               // batch size
    auto save_batch = false;        // whether to save bacthes
    auto exposure = 0.0f;           // exposure
    auto gamma = 2.2;               // gamma
    auto filmic = false;            // filmic
    auto double_sided = false;      // double sided
    auto add_skyenv = false;        // add environment
    auto quiet = false;             // quiet mode

    // parse command line
    CLI::App parser("Offline path tracing", "ytrace");
    parser.add_option("--camera", camid, "Camera index.");
    parser.add_option(
        "--resolution,-r", resolution, "Image vertical resolution.");
    parser.add_option("--nsamples,-s", nsamples, "Number of samples.");
    parser.add_option("--tracer,-t", tracer, "Trace type.")
        ->check([](const std::string& s) -> std::string {
            if (tracer_names.find(s) == tracer_names.end())
                throw CLI::ValidationError("unknown tracer name");
            return s;
        });
    parser.add_option("--nbounces", nbounces, "Maximum number of bounces.");
    parser.add_option("--pixel-clamp", pixel_clamp, "Final pixel clamping.");
    parser.add_flag("--noparallel", noparallel, "Disable parallel execution.");
    parser.add_option("--seed", seed, "Seed for the random number generators.");
    parser.add_option("--nbatch", nbatch, "Sample batch size.");
    parser.add_flag("--save-batch", save_batch, "Save images progressively");
    parser.add_option("--exposure,-e", exposure, "Hdr exposure");
    parser.add_flag(
        "--double-sided,-D", double_sided, "Double-sided rendering.");
    parser.add_flag("--add-skyenv,-E", add_skyenv, "add missing env map");
    parser.add_flag("--quiet,-q", quiet, "Print only errors messages");
    parser.add_option("--output-image,-o", imfilename, "Image filename");
    parser.add_option("scene", filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // scene loading
    auto scn = std::shared_ptr<ygl::scene>();
    if (!quiet) std::cout << "loading scene" << filename << "\n";
    auto load_start = ygl::get_time();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    if (!quiet)
        std::cout << "loading in "
                  << ygl::format_duration(ygl::get_time() - load_start) << "\n";

    // tesselate
    if (!quiet) std::cout << "tesselating scene elements\n";
    ygl::update_tesselation(scn);

    // update bbox and transforms
    ygl::update_transforms(scn);
    ygl::update_bbox(scn);

    // add components
    if (!quiet) std::cout << "adding scene elements\n";
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_sky_environment("sky"));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : scn->materials) mat->double_sided = true;
    if (scn->cameras.empty())
        scn->cameras.push_back(ygl::make_bbox_camera("<view>", scn->bbox));
    ygl::add_missing_names(scn);
    for (auto err : ygl::validate(scn)) std::cout << "warning: " << err << "\n";

    // build bvh
    if (!quiet) std::cout << "building bvh\n";
    auto bvh_start = ygl::get_time();
    ygl::update_bvh(scn);
    if (!quiet)
        std::cout << "building bvh in "
                  << ygl::format_duration(ygl::get_time() - bvh_start) << "\n";

    // init renderer
    if (!quiet) std::cout << "initializing lights\n";
    ygl::update_lights(scn);

    // initialize rendering objects
    if (!quiet) std::cout << "initializing tracer data\n";
    auto tracef = tracer_names.at(tracer);
    auto st = ygl::trace_state();

    // render
    if (!quiet) std::cout << "rendering image\n";
    auto render_start = ygl::get_time();
    auto done = false;
    while (!done) {
        if (!quiet)
            std::cout << "rendering sample " << st.sample << "/" << nsamples
                      << "\n";
        auto block_start = ygl::get_time();
        done = ygl::trace_samples(st, scn, camid, resolution, nsamples, tracef,
            nbatch, nbounces, pixel_clamp, noparallel, seed);
        if (!quiet)
            std::cout << "rendering block in "
                      << ygl::format_duration(ygl::get_time() - block_start)
                      << "\n";
        if (save_batch) {
            auto filename = ygl::replace_extension(
                imfilename, std::to_string(st.sample) + "." +
                                ygl::get_extension(imfilename));
            if (!quiet) std::cout << "saving image " << filename << "\n";
            if (ygl::is_hdr_filename(filename)) {
                ygl::save_image(filename, st.img);
            } else {
                ygl::save_image(filename,
                    ygl::tonemap_image(st.img, exposure, gamma, filmic));
            }
        }
    }
    if (!quiet)
        std::cout << "rendering image in "
                  << ygl::format_duration(ygl::get_time() - render_start)
                  << "\n";
    // save image
    if (!quiet) std::cout << "saving image " << imfilename << "\n";
    if (ygl::is_hdr_filename(imfilename)) {
        ygl::save_image(imfilename, st.img);
    } else {
        ygl::save_image(
            imfilename, ygl::tonemap_image(st.img, exposure, gamma, filmic));
    }

    // done
    return 0;
}
