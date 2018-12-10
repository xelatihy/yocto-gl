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
    auto skip_textures = parse_argument(parser,
        "--skip-textures/--no-skip-textures", false, "Disable textures.");
    auto skip_meshes   = parse_argument(
        parser, "--skip-meshes/--no-skip-meshes", false, "Disable meshes.");
    auto mesh_filenames = parse_argument(parser,
        "--mesh-filenames/--no-mesh-filenames", true, "Add mesh filenames.");
    auto mesh_directory = parse_argument(parser, "--mesh-directory", "models/"s,
        "Mesh directory when adding names.");
    auto uniform_txt    = parse_argument(parser,
        "--uniform-texture/--no-uniform-textures", false,
        "uniform texture formats");
    auto print_info     = parse_argument(
        parser, "--print-info,-i", false, "print scene info");
    auto output = parse_argument(
        parser, "--output,-o", "out.json"s, "output scene", true);
    auto filename = parse_argument(
        parser, "scene", "scene.json"s, "input scene", true);
    check_cmdline(parser);

    // fix options
    auto load_options          = load_scene_options();
    auto save_options          = save_scene_options();
    load_options.skip_textures = skip_textures;
    save_options.skip_textures = skip_textures;
    load_options.skip_meshes   = skip_meshes;
    save_options.skip_meshes   = skip_meshes;

    // load scene
    auto scene = yocto_scene{};
    if (!load_scene(filename, scene, load_options))
        log_fatal("cannot load scene {}", filename);

    // validate scene
    log_validation_errors(scene, true);

    // print info
    if (print_info) cout << print_scene_stats(scene) << "\n";

    // change texture names
    if (uniform_txt) {
        for (auto& texture : scene.textures) {
            auto ext = get_extension(texture.filename);
            if (is_hdr_filename(texture.filename)) {
                if (ext == "hdr" || ext == "exr") continue;
                if (ext == "pfm")
                    replace_extension(filename, "hdr");
                else
                    log_error("unknown texture format {}", ext);
            } else {
                if (ext == "png" || ext == "jpg") continue;
                if (ext == "tga" || ext == "bmp")
                    replace_extension(filename, "png");
                else
                    log_error("unknown texture format {}", ext);
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
            shape.filename = mesh_directory + shape.name + ".ply";
        }
        for (auto& surface : scene.surfaces) {
            surface.filename = "";
            if (surface.positions.size() <= 16) continue;
            surface.filename = mesh_directory + surface.name + ".obj";
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
        log_fatal("cannot create directory {}", get_dirname(output));
    }

    // save scene
    if (!save_scene(output, scene, save_options))
        log_fatal("cannot save scene {}", output);

    // done
    return 0;
}
