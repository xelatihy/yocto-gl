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
    // trace options
    auto params = trace_params();

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "Offline path tracing", "ytrace");
    params.camera_id   = parse_arg(parser, "--camera", 0, "Camera index.");
    params.image_size  = {0,
        parse_arg(parser, "--resolution,-r", 512, "Image vertical resolution.")};
    params.num_samples = parse_arg(
        parser, "--nsamples,-s", 256, "Number of samples.");
    params.sample_tracer = parse_arge(parser, "--tracer,-t", trace_type::path,
        "Trace type.", trace_type_names);
    params.max_bounces   = parse_arg(
        parser, "--nbounces", 8, "Maximum number of bounces.");
    params.pixel_clamp = parse_arg(
        parser, "--pixel-clamp", 100.0f, "Final pixel clamping.");
    params.no_parallel = parse_arg(
        parser, "--noparallel", false, "Disable parallel execution.");
    params.random_seed = parse_arg(
        parser, "--seed", 13, "Seed for the random number generators.");
    params.samples_per_batch = parse_arg(
        parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch = parse_arg(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure = parse_arg(parser, "--exposure,-e", 0.0f, "Hdr exposure");
    auto filmic   = parse_arg(parser, "--filmic", false, "Hdr filmic");
    auto srgb     = parse_arg(parser, "--no-srgb", true, "No srgb");
    auto embree   = parse_arg(parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = parse_arg(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto add_skyenv = parse_arg(
        parser, "--add-skyenv,-E", false, "add missing environment map");
    auto imfilename = parse_arg(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    auto filename = parse_arg(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    auto scene = yocto_scene{};
    if (!load_scene(filename, scene))
        log_fatal("cannot load scene {}", filename);

    // tesselate
    tesselate_shapes_and_surfaces(scene);

    // add components
    if (add_skyenv && scene.environments.empty()) add_sky_environment(scene);
    if (double_sided)
        for (auto& material : scene.materials) material.double_sided = true;
    log_validation_errors(scene);

    // build bvh
    auto bvh = make_scene_bvh(scene, true, embree);

    // init renderer
    auto lights = make_trace_lights(scene);

    // fix renderer type if no lights
    if (empty(lights) && params.sample_tracer != trace_type::eyelight) {
        log_info("no lights presents, switching to eyelight shader");
        params.sample_tracer = trace_type::eyelight;
    }

    // initialize rendering objects
    auto image_size = get_camera_image_size(
        scene.cameras[params.camera_id], params.image_size);
    auto rendered_image = image<vec4f>{image_size};
    auto trace_rngs     = make_trace_rngs(image_size, params.random_seed);

    // render
    auto scope = log_trace_begin("rendering image");
    for (auto sample = 0; sample < params.num_samples;
         sample += params.samples_per_batch) {
        trace_samples(rendered_image, scene, bvh, lights, sample,
            params.samples_per_batch, trace_rngs, params);
        if (save_batch) {
            auto filename = replace_extension(
                imfilename, to_string(sample + params.samples_per_batch) + "." +
                                get_extension(imfilename));
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
