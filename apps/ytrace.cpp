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

int main(int argc, char* argv[]) {
    // options
    load_scene_options  load_options  = {};
    build_bvh_options   bvh_options   = {};
    trace_image_options trace_options = {};

    // parse command line
    auto parser = cmdline_parser{};
    init_cmdline_parser(parser, argc, argv, "Offline path tracing", "ytrace");
    trace_options.camera_id = parse_cmdline_argument(
        parser, "--camera", 0, "Camera index.");
    trace_options.image_width = parse_cmdline_argument(
        parser, "--hres,-R", 1280, "Image horizontal resolution.");
    trace_options.image_height = parse_cmdline_argument(
        parser, "--vres,-r", 720, "Image vertical resolution.");
    trace_options.num_samples = parse_cmdline_argument(
        parser, "--nsamples,-s", 256, "Number of samples.");
    trace_options.sampler_type = parse_cmdline_argument(parser, "--tracer,-t",
        trace_sampler_type::path, "Trace type.", trace_sampler_type_names);
    trace_options.max_bounces  = parse_cmdline_argument(
        parser, "--nbounces", 8, "Maximum number of bounces.");
    trace_options.pixel_clamp = parse_cmdline_argument(
        parser, "--pixel-clamp", 10.0f, "Final pixel clamping.");
    auto no_parallel          = parse_cmdline_argument(parser,
        "--parallel/--no-parallel", false, "Disable parallel execution.");
    trace_options.random_seed = parse_cmdline_argument(
        parser, "--seed", 13, "Seed for the random number generators.");
    trace_options.samples_per_batch = parse_cmdline_argument(
        parser, "--nbatch,-b", 16, "Samples per batch.");
    trace_options.environments_hidden = parse_cmdline_argument(parser,
        "--env-hidden/--no-env-hidden", false,
        "Environments are hidden in renderer");
    trace_options.double_sided        = parse_cmdline_argument(parser,
        "--double-sided/--no-double-sided,-D", false,
        "Double-sided rendering.");
    auto save_batch                   = parse_cmdline_argument(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure = parse_cmdline_argument(
        parser, "--exposure,-e", 0.0f, "Hdr exposure");
    auto filmic = parse_cmdline_argument(
        parser, "--filmic", false, "Hdr filmic");
    auto srgb = parse_cmdline_argument(parser, "--no-srgb", true, "No srgb");
    bvh_options.use_embree = parse_cmdline_argument(
        parser, "--embree/--no-embree", false, "Use Embree ratracer");
    bvh_options.flatten_embree = parse_cmdline_argument(parser,
        "--flatten-embree/--no-flatten-embree", true, "Flatten embree scene");
    auto add_skyenv            = parse_cmdline_argument(
        parser, "--add-skyenv/--no-add-skyenv", false, "Add sky envmap");
    auto imfilename = parse_cmdline_argument(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    auto filename = parse_cmdline_argument(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline_parser(parser);

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
    tesselate_shapes_and_surfaces(scene);

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
    if ((empty(lights.instances) && empty(lights.environments)) &&
        is_trace_sampler_lit(trace_options)) {
        printf("no lights presents, switching to eyelight shader\n");
        trace_options.sampler_type = trace_sampler_type::eyelight;
    }

    // allocate buffers
    auto [width, height] = get_camera_image_size(
        scene.cameras[trace_options.camera_id], trace_options.image_width,
        trace_options.image_height);
    auto image = yocto::image{width, height, zero4f};
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
