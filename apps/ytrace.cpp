//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

int main(int argc, char* argv[]) {
    // options
    auto load_options  = load_scene_options{};
    auto bvh_options   = build_bvh_options{};
    auto trace_options = trace_image_options{};
    auto no_parallel   = false;
    auto save_batch    = false;
    auto exposure      = 0.0f;
    auto filmic        = false;
    auto srgb          = true;
    auto add_skyenv    = false;
    auto imfilename    = "out.hdr"s;
    auto filename      = "scene.json"s;

    // names for enums
    auto trace_sampler_type_namemap = map<string, trace_sampler_type>{};
    for (auto type = 0; type < trace_sampler_type_names.size(); type++) {
        trace_sampler_type_namemap[trace_sampler_type_names[type]] =
            (trace_sampler_type)type;
    }

    // parse command line
    auto parser = CLI::App{"Offline path tracing"};
    parser.add_option("--camera", trace_options.camera_id, "Camera index.");
    parser.add_option(
        "--hres,-R", trace_options.image_width, "Image horizontal resolution.");
    parser.add_option(
        "--vres,-r", trace_options.image_height, "Image vertical resolution.");
    parser.add_option(
        "--nsamples,-s", trace_options.num_samples, "Number of samples.");
    parser.add_option("--tracer,-t", trace_options.sampler_type, "Trace type.")
        ->transform(CLI::IsMember(trace_sampler_type_namemap));
    parser.add_option(
        "--nbounces", trace_options.max_bounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", trace_options.pixel_clamp, "Final pixel clamping.");
    parser.add_flag("--parallel,!--no-parallel", no_parallel,
        "Disable parallel execution.");
    parser.add_option("--seed", trace_options.random_seed,
        "Seed for the random number generators.");
    parser.add_option(
        "--nbatch,-b", trace_options.samples_per_batch, "Samples per batch.");
    parser.add_flag("--env-hidden,!--no-env-hidden",
        trace_options.environments_hidden,
        "Environments are hidden in renderer");
    parser.add_flag("--double-sided,!--no-double-sided,-D",
        trace_options.double_sided, "Double-sided rendering.");
    parser.add_option("--save-batch", save_batch, "Save images progressively");
    parser.add_option("--exposure,-e", exposure, "Hdr exposure");
    parser.add_option("--filmic", filmic, "Hdr filmic");
    parser.add_option("--no-srgb", srgb, "No srgb");
    parser.add_flag(
        "--embree,!--no-embree", bvh_options.use_embree, "Use Embree ratracer");
    parser.add_flag("--flatten-embree,!--no-flatten-embree",
        bvh_options.flatten_embree, "Flatten embree scene");
    parser.add_flag(
        "--add-skyenv,!--no-add-skyenv", add_skyenv, "Add sky envmap");
    parser.add_option("--output-image,-o", imfilename, "Image filename");
    parser.add_option("scene", filename, "Scene filename", true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix parallel code
    if (no_parallel) {
        bvh_options.run_serially  = true;
        load_options.run_serially = true;
    }

    // scene loading
    printf("loading scene");
    auto start_load = get_time();
    auto scene      = yocto_scene{};
    try {
        load_scene(filename, scene, load_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }
    printf(" [%s]\n", format_duration(get_time() - start_load).c_str());

    // tesselate
    tesselate_shapes(scene);

    // add components
    print_validation_errors(scene);

    // add sky
    if (add_skyenv) add_sky_environment(scene);

    // build bvh
    printf("building bvh");
    auto start_bvh = get_time();
    auto bvh       = bvh_scene{};
    build_scene_bvh(scene, bvh, bvh_options);
    printf(" [%s]\n", format_duration(get_time() - start_bvh).c_str());

    // init renderer
    auto lights = trace_lights{};
    init_trace_lights(lights, scene);

    // fix renderer type if no lights
    if ((lights.instances.empty() && lights.environments.empty()) &&
        is_trace_sampler_lit(trace_options)) {
        printf("no lights presents, switching to eyelight shader\n");
        trace_options.sampler_type = trace_sampler_type::eyelight;
    }

    // allocate buffers
    auto [width, height] = get_camera_image_size(
        scene.cameras[trace_options.camera_id], trace_options.image_width,
        trace_options.image_height);
    auto image = yocto::image{{width, height}, zero4f};
    auto state = trace_state{};
    init_trace_state(state, width, height, trace_options.random_seed);

    // render
    for (auto sample = 0; sample < trace_options.num_samples;
         sample += trace_options.samples_per_batch) {
        auto nsamples = min(trace_options.samples_per_batch,
            trace_options.num_samples - sample);
        printf("rendering image [%d/%d]", sample, trace_options.num_samples);
        auto start_batch = get_time();
        trace_image_samples(
            image, state, scene, bvh, lights, sample, trace_options);
        printf(" [%s]\n", format_duration(get_time() - start_batch).c_str());
        if (save_batch) {
            auto filename = replace_extension(
                imfilename, std::to_string(sample + nsamples) + "." +
                                get_extension(imfilename));
            try {
                save_tonemapped_image(filename, image, exposure, filmic, srgb);
            } catch (const std::exception& e) {
                exit_error(e.what());
            }
        }
    }

    // save image
    printf("saving image");
    auto start_save = get_time();
    try {
        save_tonemapped_image(imfilename, image, exposure, filmic, srgb);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }
    printf(" [%s]\n", format_duration(get_time() - start_save).c_str());

    // done
    return 0;
}
