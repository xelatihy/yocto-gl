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
    {"pathtrace", ygl::trace_path}, {"direct", ygl::trace_direct},
    {"environment", ygl::trace_environment}, {"eyelight", ygl::trace_eyelight},
    {"pathtrace-nomis", ygl::trace_path_nomis},
    {"pathtrace-naive", ygl::trace_path_naive},
    {"direct-nomis", ygl::trace_direct_nomis},
    {"debug_normal", ygl::trace_debug_normal},
    {"debug_albedo", ygl::trace_debug_albedo},
    {"debug_texcoord", ygl::trace_debug_texcoord},
    {"debug_frontfacing", ygl::trace_debug_frontfacing},
    {"debug_diffuse", ygl::trace_debug_diffuse},
    {"debug_specular", ygl::trace_debug_specular},
    {"debug_roughness", ygl::trace_debug_roughness}};

int main(int argc, char* argv[]) {
    // command line parameters
    auto filename = "scene.json"s;        // scene filename
    auto imfilename = "out.hdr"s;         // image output filename
    auto camid = 0;                       // camera index
    auto resolution = 512;                // image vertical resolution
    auto nsamples = 256;                  // image samples
    auto tracer = "pathtrace"s;           // tracer algorithm
    auto nbounces = 4;                    // number of bounces
    auto pixel_clamp = 100.0f;            // pixel clamping
    auto noparallel = false;              // disable parallel
    auto seed = ygl::trace_default_seed;  // random seed
    auto nbatch = 16;                     // batch size
    auto save_batch = false;              // whether to save bacthes
    auto exposure = 0.0f;                 // exposure
    auto gamma = 2.2;                     // gamma
    auto filmic = false;                  // filmic
    auto double_sided = false;            // double sided
    auto add_skyenv = false;              // add environment
    auto embree = false;                  // use embree
    auto quiet = false;                   // quiet mode

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
    parser.add_flag("--embree", embree, "Use Embree ratracer");
    parser.add_flag("--quiet,-q", quiet, "Print only errors messages");
    parser.add_option("--output-image,-o", imfilename, "Image filename");
    parser.add_option("scene", filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // scene loading
    auto scn = (ygl::scene*)nullptr;
    if (!quiet) printf("loading scene %s\n", filename.c_str());
    auto load_start = ygl::get_time();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\n", filename.c_str());
        printf("error: %s\n", e.what());
        exit(1);
    }
    if (!quiet)
        printf("loading in %s\n",
            ygl::format_duration(ygl::get_time() - load_start).c_str());

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_sky_environment("sky"));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : scn->materials) mat->double_sided = true;
    if (scn->cameras.empty())
        scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", ygl::compute_bbox(scn)));
    ygl::add_missing_names(scn);
    for (auto err : ygl::validate(scn)) printf("warning: %s\n", err.c_str());

    // build bvh
    if (!quiet) printf("building bvh\n");
    auto bvh_start = ygl::get_time();
    ygl::build_bvh(scn);
#if YGL_EMBREE
    if (embree) ygl::build_bvh_embree(scn);
#endif
    if (!quiet)
        printf("building bvh in %s\n",
            ygl::format_duration(ygl::get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    ygl::init_lights(scn);

    // initialize rendering objects
    if (!quiet) printf("initializing tracer data\n");
    auto tracef = tracer_names.at(tracer);
    auto cam = scn->cameras.at(camid);
    auto img = ygl::make_image4f(
        ygl::image_width(cam, resolution), ygl::image_height(cam, resolution));
    auto rng = ygl::make_trace_rngs(ygl::image_width(cam, resolution),
        ygl::image_height(cam, resolution), seed);

    // render
    if (!quiet) printf("rendering image\n");
    auto render_start = ygl::get_time();
    for (auto sample = 0; sample < nsamples; sample += nbatch) {
        if (!quiet) printf("rendering sample %04d/%04d", sample, nsamples);
        auto block_start = ygl::get_time();
        ygl::trace_samples(scn, cam, nbatch, tracef, &img, &rng, sample,
            nbounces, pixel_clamp, noparallel);
        if (!quiet)
            printf("rendering block in %s\n",
                ygl::format_duration(ygl::get_time() - block_start).c_str());
        if (save_batch) {
            auto filename = ygl::replace_extension(imfilename,
                std::to_string(sample) + "." + ygl::get_extension(imfilename));
            if (!quiet) printf("saving image %s\n", filename.c_str());
            if (ygl::is_hdr_filename(filename)) {
                ygl::save_image(filename, img);
            } else {
                ygl::save_image(
                    filename, ygl::tonemap_image(img, exposure, gamma, filmic));
            }
        }
    }
    if (!quiet)
        printf("rendering image in %s\n",
            ygl::format_duration(ygl::get_time() - render_start).c_str());

    // stata
    if (!quiet) {
        printf("using %s rays in %s paths\n",
            ygl::format_num(ygl::get_trace_stats().first).c_str(),
            ygl::format_num(ygl::get_trace_stats().second).c_str());
    }

    // save image
    if (!quiet) printf("saving image %s\n", imfilename.c_str());
    if (ygl::is_hdr_filename(imfilename)) {
        ygl::save_image(imfilename, img);
    } else {
        ygl::save_image(
            imfilename, ygl::tonemap_image(img, exposure, gamma, filmic));
    }

    // cleanup
    delete scn;

    // done
    return 0;
}
