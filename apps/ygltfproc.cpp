//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#include "../yocto/yocto_gltf.h"
#include "../yocto/yocto_utils.h"

#include "../yocto/ext/json.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <string>

using namespace nlohmann;

int main(int argc, char** argv) {
    // command line params
    auto parser = yu::cmdline::make_parser(
        argc, argv, "ygltfproc", "prints info about gltf");
    auto use_scene = parse_flag(
        parser, "--use-flag", "-s", "convert to/from scene internally");
    auto save_out =
        parse_opts(parser, "--out-filename", "-o", "output filename", "");
    auto filename = parse_args(parser, "filename", "input filename", "", true);
    check_parser(parser);

    // load gltf
    auto gltf = ygltf::load_gltf(filename, true);

    // scene
    if (use_scene) {
        auto scns = ygltf::gltf_to_scenes(gltf);
        delete gltf;
        auto buffer_uri = yu::path::get_basename(filename) + ".bin";
        gltf = ygltf::scenes_to_gltf(scns, buffer_uri);
        delete scns;
    }

    // saving output
    if (!save_out.empty()) ygltf::save_gltf(save_out, gltf, true, false);

    // cleanup
    delete gltf;
    return 0;
}
