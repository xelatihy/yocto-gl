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
    params.camera_id = parse_arg(parser, "--camera", 0, "Camera index.");
    params.vertical_resolution = parse_arg(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    params.num_samples = parse_arg(
        parser, "--nsamples,-s", 256, "Number of samples.");
    params.sample_tracer = parse_arge(parser, "--tracer,-t", trace_type::path,
        "Trace type.", trace_type_names);
    params.max_bounces   = parse_arg(
        parser, "--nbounces", 8, "Maximum number of bounces.");
    params.no_parallel = parse_arg(
        parser, "--noparallel", false, "Disable parallel execution.");
    params.samples_per_batch = parse_arg(
        parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch = parse_arg(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure = parse_arg(parser, "--exposure,-e", 0.0f, "Hdr exposure");
    auto filmic   = parse_arg(parser, "--filmic", false, "Hdr filmic");
    auto srgb     = parse_arg(parser, "--no-srgb", true, "No srgb");
    auto embree   = parse_arg(parser, "--embree", false, "Use Embree ratracer");
    auto imfilename = parse_arg(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    auto filename = parse_arg(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    printf("loading scene %s\n", filename.c_str());
    auto scene = load_scene(filename);
    if (!scene) {
        printf("could not load %s\n", filename.c_str());
        exit(1);
    }

    // tesselate
    printf("tesselating scene elements\n");
    tesselate_subdivs(scene);

    // build bvh
    printf("building bvh\n");
    auto bvh = make_scene_bvh(scene, true, embree);

    // init renderer
    printf("initializing lights\n");
    auto lights = make_trace_lights(scene, params);

    // fix renderer type if no lights
    if (!lights && params.sample_tracer != trace_type::eyelight) {
        printf("no lights presents, switching to eyelight shader\n");
        params.sample_tracer = trace_type::eyelight;
    }

    // initialize rendering objects
    printf("initializing tracer data\n");
    auto state = make_trace_state(scene, params);

    // render
    printf("rendering image\n");
    auto done = false;
    while (!done) {
        printf("rendering sample %d/%d\n", state->current_sample,
            params.num_samples);
        done = trace_samples(state, scene, bvh, lights, params);
        if (save_batch) {
            auto filename = replace_extension(
                imfilename, to_string(state->current_sample) + "." +
                                get_extension(imfilename));
            printf("saving image %s\n", filename.c_str());
            save_tonemapped_image(
                filename, state->rendered_image, exposure, filmic, srgb);
        }
    }

    // save image
    printf("saving image %s\n", imfilename.c_str());
    save_tonemapped_image(
        imfilename, state->rendered_image, exposure, filmic, srgb);

    // cleanup
    delete scene;
    delete state;
    delete bvh;
    delete lights;

    // done
    return 0;
}
