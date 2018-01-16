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

void mkdir(const string& dir) {
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
#else
    system(("mkdir " + dir).c_str());
#endif
}

void rmdir(const string& dir) {
#ifndef _MSC_VER
    auto rcmd = "rm -rf " + dir;
    system(rcmd.c_str());
#else
    auto rcmd = "del " + dir + "\\*.*; rmdir " + dir;
    system(rcmd.c_str());
#endif
}

void save_test_scene(const string& sname, const string& basedir) {
    auto dirname = basedir + "/" + sname + "/";
    printf("generating %s scenes ...\n", sname.c_str());
    try {
        auto scn = (sname == "cornell_box") ? make_cornell_box_scene() : make_prim_scene(test_scene_presets().at(sname));
        mkdir(dirname);
        if (sname == "textures") {
            for (auto txt : scn->textures) {
                if (txt->hdr) save_image4f(dirname + txt->path, txt->hdr);
                if (txt->ldr) save_image4b(dirname + txt->path, txt->ldr);
            }
        } else if (sname == "shapes") {
            for (auto shp : scn->shapes) {
                auto sscn = new scene();
                sscn->shapes += shp;
                auto mat = shp->mat;
                shp->mat = nullptr;
                auto shp_name = partition(shp->name, "_")[0];
                auto opts = save_options();
                save_scene(dirname + shp_name + ".obj", sscn, opts);
                shp->mat = mat;
                sscn->shapes.clear();
                delete sscn;
            }
        } else {
            auto facevarying = false;
            for (auto shp : scn->shapes)
                facevarying = facevarying || !shp->quads_pos.empty();

            auto opts = save_options();
            opts.save_textures = true;
            if (!facevarying) save_scene(dirname + sname + ".gltf", scn, opts);
            if (!startswith(sname, "instance")) flatten_instances(scn);
            save_scene(dirname + sname + ".obj", scn, opts);
        }
        delete scn;
    } catch (exception& e) { log_fatal("error {}", e.what()); }
}

int main(int argc, char* argv[]) {
    // put together scene names
    auto scene_names = vector<string>{"cornell_box"};
    for (auto& kv : test_scene_presets()) scene_names += kv.first;

    // command line params
    auto parser = make_parser(argc, argv, "ytestgen", "make tests");
    auto scene = parse_opt(
        parser, "--scene", "-s", "scene name", ""s, false, scene_names);
    auto clean = parse_flag(parser, "--clean", "-c", "clean directory");
    auto dirname =
        parse_opt(parser, "--dirname", "-d", "directory name", "tests"s);
    auto no_parallel =
        parse_flag(parser, "--no-parallel", "", "do not run in parallel");
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // make directories
    if (clean) rmdir(dirname);
    mkdir(dirname);

    if(scene != "") {
        save_test_scene(scene, dirname);
    } else if (no_parallel) {
        for (auto scn : scene_names) save_test_scene(scn, dirname);
    } else {
        parallel_for(scene_names.size(), [&scene_names,dirname](int idx) {
                save_test_scene(scene_names[idx], dirname);
        });
    }
}
