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

int main(int argc, char* argv[]) {
    // trace options
    auto prm = ygl::trace_params();

    // parse command line
    auto parser =
        ygl::make_cmdline_parser(argc, argv, "Offline path tracing", "ytrace");
    prm.camid = ygl::parse_int(parser, "--camera", 0, "Camera index.");
    prm.yresolution = ygl::parse_int(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    prm.nsamples =
        ygl::parse_int(parser, "--nsamples,-s", 256, "Number of samples.");
    prm.tracer = (ygl::trace_type)ygl::parse_enum(
        parser, "--tracer,-t", 0, "Trace type.", ygl::trace_type_names);
    prm.nbounces =
        ygl::parse_int(parser, "--nbounces", 8, "Maximum number of bounces.");
    prm.pixel_clamp =
        ygl::parse_float(parser, "--pixel-clamp", 100, "Final pixel clamping.");
    prm.noparallel = ygl::parse_flag(
        parser, "--noparallel", false, "Disable parallel execution.");
    prm.seed = ygl::parse_int(
        parser, "--seed", 13, "Seed for the random number generators.");
    prm.nbatch =
        ygl::parse_int(parser, "--nbatch,-b", 16, "Samples per batch.");
    auto save_batch = ygl::parse_flag(
        parser, "--save-batch", false, "Save images progressively");
    auto exposure =
        ygl::parse_float(parser, "--exposure,-e", 0, "Hdr exposure");
    auto gamma = ygl::parse_float(parser, "--gamma,-g", 2.2f, "Hdr gamma");
    auto filmic = ygl::parse_flag(parser, "--filmic", false, "Hdr filmic");
    auto embree =
        ygl::parse_flag(parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided,-D", false, "Double-sided rendering.");
    auto add_skyenv = ygl::parse_flag(
        parser, "--add-skyenv,-E", false, "add missing env map");
    auto imfilename = ygl::parse_string(
        parser, "--output-image,-o", "out.hdr", "Image filename");
    auto filename = ygl::parse_string(
        parser, "scene", "scene.json", "Scene filename", true);
    ygl::check_cmdline(parser);

    // scene loading
    auto scn = (ygl::scene*)nullptr;
    printf("loading scene %s\n", filename.c_str());
    auto load_start = ygl::get_time();
    try {
        scn = ygl::load_scene(filename);
    } catch (const std::exception& e) {
        printf("cannot load scene %s\nerror: %s\n", filename.c_str(), e.what());
        exit(1);
    }
    printf("loading in %s\n",
        ygl::format_duration(ygl::get_time() - load_start).c_str());

    // tesselate
    printf("tesselating scene elements\n");
    ygl::tesselate_subdivs(scn);

    // add components
    printf("adding scene elements\n");
    if (add_skyenv && scn->environments.empty()) {
        scn->environments.push_back(ygl::make_sky_environment("sky"));
        scn->textures.push_back(scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : scn->materials) mat->double_sided = true;
    for (auto& err : ygl::validate(scn)) printf("warning: %s\n", err.c_str());

    // build bvh
    printf("building bvh\n");
    auto bvh_start = ygl::get_time();
    ygl::build_bvh(scn);
#if YGL_EMBREE
    if (embree) ygl::build_bvh_embree(scn);
#endif
    printf("building bvh in %s\n",
        ygl::format_duration(ygl::get_time() - bvh_start).c_str());

    // init renderer
    printf("initializing lights\n");
    ygl::init_lights(scn);

    // initialize rendering objects
    printf("initializing tracer data\n");
    auto stt = ygl::make_trace_state(scn, prm);

    // render
    printf("rendering image\n");
    auto render_start = ygl::get_time();
    auto done = false;
    while (!done) {
        printf("rendering sample %04d/%04d\n", stt->sample, prm.nsamples);
        auto block_start = ygl::get_time();
        done = ygl::trace_samples(stt, scn, prm);
        printf("rendering block in %s\n",
            ygl::format_duration(ygl::get_time() - block_start).c_str());
        if (save_batch) {
            auto filename = ygl::replace_extension(
                imfilename, std::to_string(stt->sample) + "." +
                                ygl::get_extension(imfilename));
            printf("saving image %s\n", filename.c_str());
            ygl::save_tonemapped_image4f(
                filename, stt->img, exposure, gamma, filmic);
        }
    }
    printf("rendering image in %s\n",
        ygl::format_duration(ygl::get_time() - render_start).c_str());

    // save image
    printf("saving image %s\n", imfilename.c_str());
    ygl::save_tonemapped_image4f(imfilename, stt->img, exposure, gamma, filmic);

    // cleanup
    delete stt;
    delete scn;

    // done
    return 0;
}
