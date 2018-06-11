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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_utils.h"
using namespace std::literals;

void mkdir(const std::string& dir) {
    if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
        return;
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
#else
    system(("mkdir " + dir).c_str());
#endif
}

int main(int argc, char** argv) {
    // command line params
    auto parser =
        ygl::make_parser(argc, argv, "yscnproc", "converts obj to gltf");
    auto notextures =
        ygl::parse_flag(parser, "--notextures", "", "Disable textures.", false);
    auto output = ygl::parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filename = ygl::parse_arg(parser, "scene", "input scene", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load obj
    auto scn = (ygl::scene*)nullptr;
    try {
        scn = ygl::load_scene(filename, !notextures);

    } catch (const std::exception& e) {
        ygl::log_error("error during scene loading: "s + e.what());
        ygl::log_fatal("unable to load " + filename);
    }

    // make a directory if needed
    try {
        mkdir(ygl::path_dirname(output));
    } catch (const std::exception& e) {
        ygl::log_fatal("unable to make directory {} with error {}",
            ygl::path_dirname(output), e.what());
    }
    // save scene
    try {
        ygl::save_scene(output, scn);
    } catch (const std::exception& e) {
        ygl::log_fatal(
            "unable to save scene {} with error {}", output, e.what());
    }

    // done
    return 0;
}
