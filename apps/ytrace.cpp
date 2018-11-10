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
using namespace ygl;

int main(int argc, char* argv[]) {
    // options
    load_scene_options  load_options  = {};
    build_bvh_options   bvh_options   = {};
    trace_image_options trace_options = {};

    // parse command line
    auto parser = cmdline_parser{};
    init_cmdline_parser(parser, argc, argv, "Offline path tracing", "ytrace");
    trace_options.camera_id = parse_argument(
        parser, "--camera", 0, "Camera index.");
    trace_options.image_size = vec2i{0, parse_argument(parser, "--resolution,-r",
                                            512, "Image vertical resolution.")};
    trace_options.num_samples = parse_argument(
        parser, "--nsamples,-s", 256, "Number of samples.");
    trace_options.sampler_type = parse_argument(parser, "--tracer,-t",
        trace_sampler_type::path, "Trace type.", trace_sampler_type_names);
    trace_options.max_bounces  = parse_argument(
        parser, "--nbounces", 8, "Maximum number of bounces.");
    trace_options.pixel_clamp = parse_argument(
        parser, "--pixel-clamp", 100.0f, "Final pixel clamping.");
    auto no_parallel = parse_argument(
        parser, "--noparallel", false, "Disable parallel execution.");
    trace_options.random_seed = parse_argument(
        parser, "--seed", 13, "Seed for the random number generators.");
    trace_options.samples_per_batch = parse_argument(
        parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch = parse_argument(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure = parse_argument(parser, "--exposure,-e", 0.0f, "Hdr exposure");
    auto filmic   = parse_argument(parser, "--filmic", false, "Hdr filmic");
    auto srgb     = parse_argument(parser, "--no-srgb", true, "No srgb");
    bvh_options.use_embree = parse_argument(
        parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = parse_argument(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto add_skyenv = parse_argument(
        parser, "--add-skyenv,-E", false, "add missing environment map");
    auto imfilename = parse_argument(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    auto filename = parse_argument(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // fix parallel code
    if (no_parallel) {
        bvh_options.run_serially  = true;
        load_options.run_serially = true;
    }

    // scene loading
    auto scene = yocto_scene{};
    if (!load_scene(filename, scene, load_options))
        log_fatal("cannot load scene {}", filename);

    // tesselate
    tesselate_shapes_and_surfaces(scene);

    // add components
    if (add_skyenv && scene.environments.empty()) add_sky_environment(scene);
    if (double_sided)
        for (auto& material : scene.materials) material.double_sided = true;
    log_validation_errors(scene);

    // build bvh
    auto bvh = bvh_scene{};
    build_scene_bvh(scene, bvh, bvh_options);

    // init renderer
    auto lights = make_trace_lights(scene);

    // fix renderer type if no lights
    if (empty(lights) &&
        trace_options.sampler_type != trace_sampler_type::eyelight) {
        log_info("no lights presents, switching to eyelight shader");
        trace_options.sampler_type = trace_sampler_type::eyelight;
    }

    // allocate buffers
    auto image_size = get_camera_image_size(
        scene.cameras[trace_options.camera_id], trace_options.image_size);
    auto rendered_image = image<vec4f>{image_size};
    auto trace_pixels = make_trace_pixels(image_size, trace_options.random_seed);

    // render
    auto scope = log_trace_begin("rendering image");
    for (auto sample = 0; sample < trace_options.num_samples;
         sample += trace_options.samples_per_batch) {
        auto nsamples = min(trace_options.samples_per_batch,
            trace_options.num_samples - sample);
        trace_image_samples(rendered_image, trace_pixels, scene, bvh, lights,
            sample, trace_options);
        if (save_batch) {
            auto filename = replace_extension(imfilename,
                to_string(sample + nsamples) + "." + get_extension(imfilename));
            if (!save_tonemapped_image(
                    filename, rendered_image, exposure, filmic, srgb))
                log_fatal("cannot save image " + filename);
        }
    }
    log_trace_end(scope);

    // save image
    if (!save_tonemapped_image(imfilename, rendered_image, exposure, filmic, srgb))
        log_fatal("cannot save image " + imfilename);

    // done
    return 0;
}
