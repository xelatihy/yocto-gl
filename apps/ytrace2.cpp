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

auto tracer_names = std::vector<std::string>{"pathtrace", "direct",
    "environment", "eyelight", "pathtrace-nomis", "pathtrace-naive",
    "direct-nomis", "debug_normal", "debug_albedo", "debug_texcoord",
    "debug_frontfacing", "debug_diffuse", "debug_specular", "debug_roughness"};

auto tracer_funcs = std::vector<ygl::trace_func>{ygl::trace_path,
    ygl::trace_direct, ygl::trace_environment, ygl::trace_eyelight,
    ygl::trace_path_nomis, ygl::trace_path_naive, ygl::trace_direct_nomis,
    ygl::trace_debug_normal, ygl::trace_debug_albedo, ygl::trace_debug_texcoord,
    ygl::trace_debug_frontfacing, ygl::trace_debug_diffuse,
    ygl::trace_debug_specular, ygl::trace_debug_roughness};

int main(int argc, char* argv[]) {
    // parse command line
    auto parser =
        ygl::make_cmdline_parser(argc, argv, "Offline path tracing", "ytrace");
    auto camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    auto resolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    auto nsamples =
        ygl::parse_int(parser, "--nsamples,-s", 256, "Number of samples.");
    auto tracer_id =
        ygl::parse_enum(parser, "--tracer,-t", 0, "Trace type.", tracer_names);
    auto nbounces =
        ygl::parse_int(parser, "--nbounces", 8, "Maximum number of bounces.");
    auto noparallel = ygl::parse_flag(
        parser, "--noparallel", false, "Disable parallel execution.");
    auto nbatch =
        ygl::parse_int(parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch = ygl::parse_flag(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure =
        ygl::parse_float(parser, "--exposure,-e", 0, "Hdr exposure");
    auto gamma = ygl::parse_float(parser, "--gamma,-g", 2.2f, "Hdr gamma");
    auto filmic = ygl::parse_flag(parser, "--filmic", false, "Hdr filmic");
    auto embree =
        ygl::parse_flag(parser, "--embree", false, "Use Embree ratracer");
    auto imfilename = ygl::parse_string(
        parser, "--output-image,-o", "out.hdr", "Image filename");
    auto filename = ygl::parse_string(
        parser, "scene", "scene.json", "Scene filename", true);
    ygl::check_cmdline(parser);

    // scene loading
    printf("loading scene %s\n", filename.c_str());
    auto scn = ygl::load_scene(filename);

    // tesselate
    printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(scn);

    // build bvh
    printf("building bvh\n");
    ygl::build_bvh(scn);
#if YGL_EMBREE
    if (embree) ygl::build_bvh_embree(scn);
#endif

    // init renderer
    printf("initializing lights\n");
    ygl::init_lights(scn);

    // initialize rendering objects
    printf("initializing tracer data\n");
    auto tracer_func = tracer_funcs.at(tracer_id);
    auto cam = scn->cameras.at(camid);
    auto img = ygl::make_image4f(
        ygl::image_width(cam, resolution), ygl::image_height(cam, resolution));
    auto rng = ygl::make_trace_rngs(
        ygl::image_width(cam, resolution), ygl::image_height(cam, resolution));

    // render
    printf("rendering image\n");
    for (auto sample = 0; sample < nsamples; sample += nbatch) {
        printf("rendering sample %d/%d\n", sample, nsamples);
        ygl::trace_samples(scn, cam, nbatch, tracer_func, img, rng, sample,
            nbounces, 100, noparallel);
        if (save_batch) {
            auto filename = ygl::replace_extension(imfilename,
                std::to_string(sample) + "." + ygl::get_extension(imfilename));
            printf("saving image %s\n", filename.c_str());
            ygl::save_tonemapped_image4f(
                filename, img, exposure, gamma, filmic);
        }
    }

    // save image
    printf("saving image %s\n", imfilename.c_str());
    ygl::save_tonemapped_image4f(imfilename, img, exposure, gamma, filmic);

    // cleanup
    delete scn;

    // done
    return 0;
}
