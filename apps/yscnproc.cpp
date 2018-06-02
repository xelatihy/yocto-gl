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

void validate_texture_paths(ygl::scene* obj, const std::string& filename) {
    auto dirname = ygl::path_dirname(filename);
    for (auto txt : obj->textures) {
        auto f = fopen((dirname + txt->path).c_str(), "rb");
        if (!f)
            printf("Missing texture: %s\n", txt->path.c_str());
        else
            fclose(f);
    }
}

int main(int argc, char** argv) {
    // command line params
    auto parser =
        ygl::make_parser(argc, argv, "yscnproc", "converts obj to gltf");
    auto textures =
        ygl::parse_flag(parser, "--textures", "-t", "process textures");
    // auto no_flip_opacity =
    //     ygl::parse_flag(parser, "--no-flip-opacity", "", "flip opacity");
    // auto facet_non_smooth = ygl::parse_flag(
    //     parser, "--facet-non-smooth", "", "facet non smooth surfaces");
    auto scale = ygl::parse_opt(parser, "--scale", "", "scale the model", 1.0f);
    auto flipyz =
        ygl::parse_flag(parser, "--flipyz", "", "flip y and z coords");
    auto add_normals = ygl::parse_flag(parser, "--normals", "", "add normals");
    //    auto add_specgloss =
    //        ygl::parse_flag( "--specgloss", "", "add spec gloss");
    auto save_separate_buffers = ygl::parse_flag(
        parser, "--separate-buffers", "", "save separate buffers");
    auto info =
        ygl::parse_flag(parser, "--print-info", "", "print information");
    auto validate_textures = ygl::parse_flag(
        parser, "--validate-textures", "", "validate texture paths");
    auto validate =
        ygl::parse_flag(parser, "--validate", "", "validate after saving");
    auto output = ygl::parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filenames = ygl::parse_args(
        parser, "scenes", "input scene filenames", std::vector<std::string>{});
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load obj
    auto scn = new ygl::scene();
    for (auto filename : filenames) {
        auto to_merge = (ygl::scene*)nullptr;
        try {
            to_merge = ygl::load_scene(filename, textures);

        } catch (const std::exception& e) {
            ygl::log_fatal("unable to load file %s with error {}\n",
                filename.c_str(), e.what());
        }

        // check missing texture
        if (validate_textures) validate_texture_paths(to_merge, filename);

        // print information
        if (info) {
            printf("information ------------------------\n");
            printf("filename: %s\n", filename.c_str());
            ygl::print_stats(to_merge);
        }

        // merge the scene into the other
        ygl::merge_into(scn, to_merge);
        delete to_merge;
    }

    // print information
    if (info && filenames.size() > 1) {
        printf("information ------------------------\n");
        printf("merged\n");
        ygl::print_stats(scn);
    }

    // add missing elements
    if (add_normals) ygl::add_missing_normals(scn);

    // process geometry
    if (scale != 1.0f) {
        for (auto shp : scn->shapes) {
            for (auto& p : shp->pos) p *= scale;
        }
    }

    if (flipyz) {
        for (auto shp : scn->shapes) {
            for (auto& p : shp->pos) std::swap(p.y, p.z);
            for (auto& n : shp->norm) std::swap(n.y, n.z);
        }
    }

    // print infomation again if needed
    if (info && scale != 1.0f) {
        printf("post-correction information -------\n");
        printf("output: %s\n", output.c_str());
        ygl::print_stats(scn);
    }

    // make a directory if needed
    try {
        mkdir(ygl::path_dirname(output));
    } catch (const std::exception& e) {
        ygl::log_fatal("unable to make directory %s with error {}\n",
            ygl::path_dirname(output), e.what());
    }
    // save scene
    try {
        ygl::save_scene(output, scn);
    } catch (const std::exception& e) {
        ygl::log_fatal("unable to save scene %s with error {}\n",
            output.c_str(), e.what());
    }

    // validate
    if (validate) {
        auto vscn = ygl::load_scene(output);
        if (info) {
            printf("validate information -----------------\n");
            ygl::print_stats(vscn);
        }
    }

    // done
    return 0;
}
