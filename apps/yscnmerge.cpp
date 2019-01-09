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

#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_utils.h"
using namespace yocto;

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
    // parse command line
    auto parser = make_cmdline_parser(argc, argv, "Process scene", "yscnproc");
    auto scene_postfix  = parse_argument(parser,
        "--scene-postfix/no-scene-postfix", false,
        "Append unique scene postfix to each name");
    auto skip_textures  = parse_argument(parser,
        "--skip-textures/--no-skip-textures", false, "Disable textures.");
    auto mesh_filenames = parse_argument(parser,
        "--mesh-filenames/--no-mesh-filenames", true, "Add mesh filenames.");
    auto mesh_directory = parse_argument(parser, "--mesh-directory", "models/"s,
        "Mesh directory when adding names.");
    auto texture_filenames = parse_argument(parser,
        "--texture-filenames/--no-texture-filenames", true,
        "Add texture filenames.");
    auto texture_directory = parse_argument(parser, "--texture-directory",
        "textures/"s, "Texture directory when adding names.");
    auto print_info        = parse_argument(
        parser, "--print-info,-i", false, "print scene info");
    auto output = parse_argument(
        parser, "--output,-o", "out.json"s, "output scene", true);
    auto filenames = parse_arguments(
        parser, "scenes", vector<string>{}, "scene filenames", true);
    check_cmdline(parser);

    // fix options
    auto load_options          = load_scene_options();
    auto save_options          = save_scene_options();
    load_options.skip_textures = skip_textures;
    save_options.skip_textures = skip_textures;

    // load scene
    auto scene = yocto_scene{};
    for (auto& filename : filenames) {
        auto to_merge = yocto_scene{};
        if (!load_scene(filename, to_merge, load_options))
            log_fatal("cannot load scene {}", filename);
        log_validation_errors(to_merge, true);
        if (scene_postfix) {
            auto postfix = "{" + get_filename(filename) + "}";
            for (auto& val : to_merge.cameras) val.name += postfix;
            for (auto& val : to_merge.textures) val.name += postfix;
            for (auto& val : to_merge.voltextures) val.name += postfix;
            for (auto& val : to_merge.materials) val.name += postfix;
            for (auto& val : to_merge.shapes) val.name += postfix;
            for (auto& val : to_merge.surfaces) val.name += postfix;
            for (auto& val : to_merge.instances) val.name += postfix;
            for (auto& val : to_merge.environments) val.name += postfix;
            for (auto& val : to_merge.nodes) val.name += postfix;
            for (auto& val : to_merge.animations) val.name += postfix;
        }
        merge_scene_into(scene, to_merge);
    }

    // validate scene
    log_validation_errors(scene, true);

    // print info
    if (print_info) cout << print_scene_stats(scene) << "\n";

    // add missing mesh names if necessary
    if (!mesh_directory.empty() && mesh_directory.back() != '/')
        mesh_directory += '/';
    if (get_extension(output) == "json") {
        for (auto& shape : scene.shapes) {
            if (shape.filename.empty()) {
                shape.filename = mesh_directory + shape.name + ".ply";
            } else if (mesh_filenames) {
                shape.filename = mesh_directory + get_filename(shape.filename);
            }
        }
        for (auto& surface : scene.surfaces) {
            surface.filename = mesh_directory + surface.name + ".obj";
        }
    }
    // gltf does not support embedded data
    if (get_extension(output) == "gltf") {
        for (auto& shape : scene.shapes) {
            shape.filename = mesh_directory + shape.name + ".bin";
        }
    }

    // add missing textures names if necessary
    if (!texture_directory.empty() && texture_directory.back() != '/')
        texture_directory += '/';
    if (get_extension(output) == "json") {
        for (auto& texture : scene.textures) {
            if (texture.filename.empty()) {
                texture.filename = texture_directory + texture.name + ".ply";
            } else if (texture_filenames) {
                texture.filename = texture_directory +
                                   get_filename(texture.filename);
            }
        }
    }

    // make a directory if needed
    if (!mkdir(get_dirname(output))) {
        log_fatal("cannot create directory {}", get_dirname(output));
    }

    // save scene
    if (!save_scene(output, scene, save_options))
        log_fatal("cannot save scene {}", output);

    // done
    return 0;
}
