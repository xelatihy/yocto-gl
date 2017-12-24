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

#include "../yocto/yocto_gl.h"
using namespace ygl;

#include <algorithm>
#include <array>
#include <memory>

void mkdir(const string& dir) {
    if (dir == "" || dir == "." || dir == ".." || dir == "./" || dir == "../")
        return;
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
#else
    system(("mkdir " + dir).c_str());
#endif
}

void validate_texture_paths(scene* obj, const string& filename) {
    auto dirname = path_dirname(filename);
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
    using namespace string_literals;
    auto parser = make_parser(argc, argv, "yobj2gltf", "converts obj to gltf");
    auto textures = parse_flag(parser, "--textures", "-t", "process textures");
    auto no_flipy_texcoord = parse_flag(
        parser, "--no-flipy-texcoord", "", "texcoord vertical flipping");
    auto no_flip_opacity =
        parse_flag(parser, "--no-flip-opacity", "", "flip opacity");
    auto facet_non_smooth = parse_flag(
        parser, "--facet-non-smooth", "", "facet non smooth surfaces");
    auto scale = parse_opt(parser, "--scale", "", "scale the model", 1.0f);
    auto flipyz = parse_flag(parser, "--flipyz", "", "flip y and z coords");
    auto add_scene = parse_flag(parser, "--scene", "", "add scene");
    auto add_normals = parse_flag(parser, "--normals", "", "add normals");
    auto tesselation =
        parse_opt(parser, "--tesselation", "-T", "tesselation level", 0);
    auto subdiv =
        parse_opt(parser, "--subdiv", "", "Catmull-Clark subdivision", 0);
    //    auto add_specgloss =
    //        parse_flag( "--specgloss", "", "add spec gloss");
    auto save_separate_buffers =
        parse_flag(parser, "--separate-buffers", "", "save separate buffers");
    auto info = parse_flag(parser, "--print-info", "", "print information");
    auto validate_textures =
        parse_flag(parser, "--validate-textures", "", "validate texture paths");
    auto validate =
        parse_flag(parser, "--validate", "", "validate after saving");
    auto output = parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filenames =
        parse_args(parser, "scenes", "input scene filenames", vector<string>{});
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load obj
    auto scn = unique_ptr<scene>(new scene());
    for (auto filename : filenames) {
        auto to_merge = unique_ptr<scene>(nullptr);
        try {
            auto opts = load_options();
            opts.load_textures = textures;
            opts.obj_flip_texcoord = !no_flipy_texcoord;
            opts.obj_flip_tr = !no_flip_opacity;
            opts.obj_facet_non_smooth = facet_non_smooth;
            opts.preserve_quads = true;
            to_merge = unique_ptr<scene>(load_scene(filename, opts));

        } catch (const exception& e) {
            log_fatal("unable to load file %s with error {}\n",
                filename.c_str(), e.what());
        }

        // check missing texture
        if (validate_textures) validate_texture_paths(to_merge.get(), filename);

        // print information
        if (info) {
            printf("information ------------------------\n");
            printf("filename: %s\n", filename.c_str());
            print_info(to_merge.get());
        }

        // merge the scene into the other
        if (filenames.size() > 1) {
            merge_into(scn.get(), to_merge.get());
        } else {
            swap(scn, to_merge);
        }
    }

    // print information
    if (info && filenames.size() > 1) {
        printf("information ------------------------\n");
        printf("merged\n");
        print_info(scn.get());
    }

    // add missing elements
    {
        auto opts = add_elements_options::none();
        opts.default_names = true;
        opts.smooth_normals = add_normals;
        opts.shape_instances = add_scene;
        opts.default_paths = true;
        add_elements(scn.get(), opts);
    }

    // process geometry
    if (scale != 1.0f) {
        for (auto shp : scn->shapes) {
            for (auto& p : shp->pos) p *= scale;
        }
    }

    if (flipyz) {
        for (auto shp : scn->shapes) {
            for (auto& p : shp->pos) swap(p.y, p.z);
            for (auto& n : shp->norm) swap(n.y, n.z);
        }
    }

    if (tesselation) {
        for (auto shp : scn->shapes) {
            for (auto l = 0; l < tesselation; l++) { subdivide_shape(shp); }
        }
    }

    if (subdiv) {
        for (auto shp : scn->shapes) {
            for (auto l = 0; l < subdiv; l++) { subdivide_shape(shp, true); }
        }
    }

    // print infomation again if needed
    if (info && scale != 1.0f) {
        printf("post-correction information -------\n");
        printf("output: %s\n", output.c_str());
        print_info(scn.get());
    }

    // make a directory if needed
    try {
        mkdir(path_dirname(output));
    } catch (const exception& e) {
        log_fatal("unable to make directory %s with error {}\n",
            path_dirname(output), e.what());
    }
    // save scene
    try {
        auto opts = save_options();
        opts.save_textures = textures;
        opts.gltf_separate_buffers = save_separate_buffers;
        save_scene(output, scn.get(), opts);
    } catch (const exception& e) {
        log_fatal("unable to save scene %s with error {}\n", output.c_str(),
            e.what());
    }

    // validate
    if (validate) {
        auto vscn = unique_ptr<scene>(load_scene(output));
        if (info) {
            printf("validate information -----------------\n");
            print_info(vscn.get());
        }
    }

    // done
    return 0;
}
