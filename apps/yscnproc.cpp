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
    auto skip_textures    = false;
    auto skip_meshes      = false;
    auto mesh_filenames   = false;
    auto shape_directory  = "shapes/"s;
    auto subdiv_directory = "subdivs/"s;
    auto uniform_txt      = false;
    auto validate         = false;
    auto print_info       = false;
    auto output           = "out.json"s;
    auto filename         = "scene.json"s;

    // parse command line
    auto parser = CLI::App{"Process scene"};
    parser.add_flag("--skip-textures,!--no-skip-textures", skip_textures,
        "Disable textures.");
    parser.add_flag(
        "--skip-meshes,!--no-skip-meshes", skip_meshes, "Disable meshes.");
    parser.add_flag("--mesh-filenames,!--no-mesh-filenames", mesh_filenames,
        "Add mesh filenames.");
    parser.add_option("--shape-directory", shape_directory,
        "Shape directory when adding names.");
    parser.add_option("--subdiv-directory", subdiv_directory,
        "Subdiv directory when adding names.");
    parser.add_flag("--uniform-texture,!--no-uniform-textures", uniform_txt,
        "uniform texture formats");
    parser.add_flag("--print-info,-i", print_info, "print scene info");
    parser.add_flag("--validate,!--no-validate", validate, "Validate scene");
    parser.add_option("--output,-o", output, "output scene")->required(true);
    parser.add_option("scene", filename, "input scene")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix options
    auto load_options          = load_scene_options();
    auto save_options          = save_scene_options();
    load_options.skip_textures = skip_textures;
    save_options.skip_textures = skip_textures;
    load_options.skip_meshes   = skip_meshes;
    save_options.skip_meshes   = skip_meshes;

    // load scene
    printf("loading scene ...\n");
    auto start_load = get_time();
    auto scene      = yocto_scene{};
    try {
        load_scene(filename, scene, load_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }
    printf("loading scene [%s]\n",
        format_duration(get_time() - start_load).c_str());

    // validate scene
    if (validate) {
        printf("validating scene ...\n");
        auto start_val = get_time();
        print_validation_errors(scene);
        printf("validating scene [%s]\n",
            format_duration(get_time() - start_val).c_str());
    }

    // print info
    if (print_info) printf("%s\n", print_scene_stats(scene).c_str());

    // change texture names
    if (uniform_txt) {
        for (auto& texture : scene.textures) {
            auto ext = get_extension(texture.uri);
            if (is_hdr_filename(texture.uri)) {
                if (ext == "hdr") continue;
                get_noextension(filename) + ".hdr";
            } else {
                if (ext == "png") continue;
                get_noextension(filename) + ".png";
            }
        }
    }

    printf("tesselating scene ...\n");
    auto start_tess = get_time();
    tesselate_subdivs(scene);
    printf("tesselating scene [%s]\n",
        format_duration(get_time() - start_tess).c_str());

    // add missing mesh names if necessary
    if (!shape_directory.empty() && shape_directory.back() != '/')
        shape_directory += '/';
    if (mesh_filenames) {
        auto sid = 0;
        for (auto& shape : scene.shapes) {
            if (!shape.quads_positions.empty()) {
                shape.uri = shape_directory + "shape_" + std::to_string(sid) +
                            ".obj";
            } else {
                shape.uri = shape_directory + "shape_" + std::to_string(sid) +
                            ".ply";
            }
            sid++;
        }
    }
    // add missing mesh names if necessary
    if (!subdiv_directory.empty() && subdiv_directory.back() != '/')
        subdiv_directory += '/';
    if (mesh_filenames) {
        auto sid = 0;
        for (auto& subdiv : scene.subdivs) {
            if (!subdiv.quads_positions.empty()) {
                subdiv.uri = subdiv_directory + "subdiv_" +
                             std::to_string(sid) + ".obj";
            } else {
                subdiv.uri = subdiv_directory + "subdiv_" +
                             std::to_string(sid) + ".ply";
            }
            sid++;
        }
    }

    // make a directory if needed
    if (!mkdir(get_dirname(output))) {
        exit_error("cannot create directory " + get_dirname(output));
    }

    // save scene
    printf("saving scene ...\n");
    auto start_save = get_time();
    try {
        save_scene(output, scene, save_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }
    printf("saving scene [%s]\n",
        format_duration(get_time() - start_save).c_str());

    // done
    return 0;
}
