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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

#include "ext/CLI11.hpp"

bool mkdir(const string& dir) {
    if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
        return true;
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
    return true;
#else
    system(("mkdir " + dir).c_str());
    return true;
#endif
}

int main(int argc, char** argv) {
    // command line parameters
    auto scene_postfix     = false;
    auto skip_textures     = false;
    auto mesh_filenames    = true;
    auto mesh_directory    = "models/"s;
    auto texture_filenames = true;
    auto texture_directory = "textures/"s;
    auto info              = false;
    auto output            = "out.json"s;
    auto filenames         = vector<string>{};

    // parse command line
    auto parser = CLI::App{"Merge scenes"};
    parser.add_flag("--scene-postfix,!--no-scene-postfix", scene_postfix,
        "Append unique scene postfix to each name");
    parser.add_flag("--skip-textures,!--no-skip-textures", skip_textures,
        "Disable textures.");
    parser.add_flag("--mesh-filenames,!--no-mesh-filenames", mesh_filenames,
        "Add mesh filenames.");
    parser.add_option("--mesh-directory", mesh_directory,
        "Mesh directory when adding names.");
    parser.add_flag("--texture-filenames,!--no-texture-filenames",
        texture_filenames, "Add texture filenames.");
    parser.add_option("--texture-directory", texture_directory,
        "Texture directory when adding names.");
    parser.add_option("--print-info,-i", info, "print scene info");
    parser.add_option("--output,-o", output, "output scene")->required(true);
    parser.add_option("scenes", filenames, "scene filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix options
    auto load_prms          = sceneio_params();
    auto save_prms          = sceneio_params();
    load_prms.skip_textures = skip_textures;
    save_prms.skip_textures = skip_textures;

    // load scene
    auto scene = yocto_scene{};
    for (auto& filename : filenames) {
        auto to_merge = yocto_scene{};
        try {
            load_scene(filename, to_merge, load_prms);
        } catch (const std::exception& e) {
            print_fatal(e.what());
        }
        print_validation(to_merge, true);
        if (scene_postfix) {
            auto postfix = "{" + get_filename(filename) + "}";
            for (auto& val : to_merge.cameras) val.uri += postfix;
            for (auto& val : to_merge.textures) val.uri += postfix;
            for (auto& val : to_merge.voltextures) val.uri += postfix;
            for (auto& val : to_merge.materials) val.uri += postfix;
            for (auto& val : to_merge.shapes) val.uri += postfix;
            for (auto& val : to_merge.instances) val.uri += postfix;
            for (auto& val : to_merge.environments) val.uri += postfix;
            for (auto& val : to_merge.nodes) val.uri += postfix;
            for (auto& val : to_merge.animations) val.uri += postfix;
        }
        merge_scene(scene, to_merge);
    }

    // validate scene
    print_validation(scene, true);

    // print info
    if (info) print_info("{}", format_stats(scene));

    // add missing mesh names if necessary
    if (!mesh_directory.empty() && mesh_directory.back() != '/')
        mesh_directory += '/';
    if (get_extension(output) == "yaml") {
        for (auto& shape : scene.shapes) {
            shape.uri = mesh_directory + get_filename(shape.uri);
        }
    }
    // gltf does not support embedded data
    if (get_extension(output) == "gltf") {
        for (auto& shape : scene.shapes) {
            shape.uri = mesh_directory + shape.uri + ".bin";
        }
    }

    // add missing textures names if necessary
    if (!texture_directory.empty() && texture_directory.back() != '/')
        texture_directory += '/';
    if (get_extension(output) == "yaml") {
        auto tid = 0;
        for (auto& texture : scene.textures) {
            if (texture.uri.empty()) {
                texture.uri = texture_directory + "ytexture#" +
                              std::to_string(tid) + ".png";
            } else if (texture_filenames) {
                texture.uri = texture_directory + get_filename(texture.uri);
            }
            tid++;
        }
    }

    // make a directory if needed
    if (!mkdir(get_dirname(output))) {
        print_fatal("cannot create directory " + get_dirname(output));
    }

    // save scene
    try {
        save_scene(output, scene, save_prms);
    } catch (const std::exception& e) {
        print_fatal(e.what());
    }

    // done
    return 0;
}
