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

#include "../yocto/yocto_imageio.h"
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
    auto skip_textures  = false;
    auto skip_meshes    = false;
    auto mesh_filenames = true;
    auto mesh_directory = "models/"s;
    auto uniform_txt    = false;
    auto print_info     = false;
    auto output         = "out.json"s;
    auto filename       = "scene.json"s;

    // parse command line
    auto parser = CLI::App{"Process scene"};
    parser.add_flag("--skip-textures,!--no-skip-textures", skip_textures,
        "Disable textures.");
    parser.add_flag(
        "--skip-meshes,!--no-skip-meshes", skip_meshes, "Disable meshes.");
    parser.add_flag("--mesh-filenames,!--no-mesh-filenames", mesh_filenames,
        "Add mesh filenames.");
    parser.add_option("--mesh-directory", mesh_directory,
        "Mesh directory when adding names.");
    parser.add_flag("--uniform-texture,!--no-uniform-textures", uniform_txt,
        "uniform texture formats");
    parser.add_flag("--print-info,-i", print_info, "print scene info");
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
    auto scene = yocto_scene{};
    try {
        load_scene(filename, scene, load_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // validate scene
    print_validation_errors(scene, true);

    // print info
    if (print_info) printf("%s\n", print_scene_stats(scene).c_str());

    // change texture names
    if (uniform_txt) {
        for (auto& texture : scene.textures) {
            auto ext = get_extension(texture.filename);
            if (is_hdr_filename(texture.filename)) {
                if (ext == "hdr" || ext == "exr") continue;
                if (ext == "pfm")
                    replace_extension(filename, "hdr");
                else
                    throw io_error("unknown texture format " + ext);
            } else {
                if (ext == "png" || ext == "jpg") continue;
                if (ext == "tga" || ext == "bmp")
                    replace_extension(filename, "png");
                else
                    throw io_error("unknown texture format " + ext);
            }
        }
    }

    // add missing mesh names if necessary
    if (!mesh_directory.empty() && mesh_directory.back() != '/')
        mesh_directory += '/';
    if (mesh_filenames && get_extension(output) == "json") {
        for (auto& shape : scene.shapes) {
            shape.filename = "";
            if (shape.positions.size() <= 16) continue;
            if (shape.preserve_facevarying) {
                shape.filename = mesh_directory + shape.name + ".obj";
            } else {
                shape.filename = mesh_directory + shape.name + ".ply";
            }
        }
    }
    // gltf does not support embedded data
    if (get_extension(output) == "gltf") {
        for (auto& shape : scene.shapes) {
            shape.filename = mesh_directory + shape.name + ".bin";
        }
    }

    // make a directory if needed
    if (!mkdir(get_dirname(output))) {
        exit_error("cannot create directory " + get_dirname(output));
    }

    // save scene
    try {
        save_scene(output, scene, save_options);
    } catch (const std::exception& e) {
        exit_error(e.what());
    }

    // done
    return 0;
}
