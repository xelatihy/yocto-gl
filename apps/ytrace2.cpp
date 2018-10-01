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
    auto parser =
        make_cmdline_parser(argc, argv, "Offline path tracing", "ytrace");
    params.camid = parse_int(parser, "--camera", 0, "Camera index.");
    params.yresolution =
        parse_int(parser, "--resolution,-r", 512, "Image vertical resolution.");
    params.nsamples =
        parse_int(parser, "--nsamples,-s", 256, "Number of samples.");
    params.tracer = (trace_type)parse_enum(
        parser, "--tracer,-t", 0, "Trace type.", trace_type_names);
    params.nbounces =
        parse_int(parser, "--nbounces", 8, "Maximum number of bounces.");
    params.noparallel = parse_flag(
        parser, "--noparallel", false, "Disable parallel execution.");
    params.nbatch = parse_int(parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch =
        parse_flag(parser, "--save-batch", false, "Save images progressively");
    auto exposure = parse_float(parser, "--exposure,-e", 0, "Hdr exposure");
    auto gamma = parse_float(parser, "--gamma,-g", 2.2f, "Hdr gamma");
    auto filmic = parse_flag(parser, "--filmic", false, "Hdr filmic");
    auto embree = parse_flag(parser, "--embree", false, "Use Embree ratracer");
    auto imfilename =
        parse_string(parser, "--output-image,-o", "out.hdr", "Image filename");
    auto filename =
        parse_string(parser, "scene", "scene.json", "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    printf("loading scene %s\n", filename.c_str());
    auto scn = load_scene(filename);

    // tesselate
    printf("tesselating scene elements\n");
    tesselate_subdivs(scn);

    // build bvh
    printf("building bvh\n");
    auto bvh = build_bvh(scn, true, embree);

    // init renderer
    printf("initializing lights\n");
    auto lights = make_trace_lights(scn, params);

    // initialize rendering objects
    printf("initializing tracer data\n");
    auto state = make_trace_state(scn, params);

    // render
    printf("rendering image\n");
    auto done = false;
    while (!done) {
        printf("rendering sample %d/%d\n", state->sample, params.nsamples);
        done = trace_samples(state, scn, bvh, lights, params);
        if (save_batch) {
            auto filename = replace_extension(
                imfilename, std::to_string(state->sample) + "." +
                                get_extension(imfilename));
            printf("saving image %s\n", filename.c_str());
            save_tonemapped_image(
                filename, state->img, exposure, gamma, filmic);
        }
    }

    // save image
    printf("saving image %s\n", imfilename.c_str());
    save_tonemapped_image(imfilename, state->img, exposure, gamma, filmic);

    // cleanup
    delete scn;
    delete state;
    delete bvh;
    delete lights;

    // done
    return 0;
}
