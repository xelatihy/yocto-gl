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
using namespace std::literals;

void mkdir(const std::string& dir) {
#ifndef _MSC_VER
    system(("mkdir -p " + dir).c_str());
#else
    auto fdir = dir;
    for (auto& c : fdir)
        if (c == '/') c = '\\';
    system(("mkdir " + fdir).c_str());
#endif
}

void rmdir(const std::string& dir) {
#ifndef _MSC_VER
    auto rcmd = "rm -rf " + dir;
    system(rcmd.c_str());
#else
    auto rcmd = "del " + dir + "\\*.*; rmdir " + dir;
    system(rcmd.c_str());
#endif
}

void save_scene(ygl::scene* scn, const std::string& sname,
    const std::string& dirname, bool flatten_obj) {
    auto facevarying = false;
    for (auto sgr : scn->shapes)
        for (auto shp : sgr->shapes)
            facevarying = facevarying || !shp->quads_pos.empty();

    auto opts = ygl::save_options();
    opts.save_textures = true;
    if (!facevarying) save_scene(dirname + sname + ".gltf", scn, opts);
    if (flatten_obj) ygl::flatten_instances(scn);
    ygl::save_scene(dirname + sname + ".obj", scn, opts);
}

void save_test_scene(const std::string& sname, const std::string& basedir) {
    auto dirname = basedir + "/" + sname + "/";
    try {
        mkdir(dirname);
        ygl::log_info("generating scene {}", sname);
        if (sname == "cornell_box") {
            auto scn = ygl::make_cornell_box_scene();
            save_scene(scn, sname, dirname, true);
            return;
        }
        auto test_scns =
            std::vector<std::pair<ygl::proc_scene*, ygl::scene*>>();
        for (auto& preset : ygl::proc_split_scene_presets()) {
            if (preset->scn->name != sname) continue;
            test_scns.push_back(
                {preset->scn, ygl::make_proc_elems(preset->scn)});
            for (auto& view : preset->views) {
                test_scns.push_back({view, ygl::make_proc_elems(view)});
            }
        }
        ygl::log_info("saving scene {}", sname);
        if (sname == "textures") {
            auto scn = test_scns.front().second;
            for (auto txt : scn->textures) {
                if (!txt->hdr.empty())
                    ygl::save_image4f(dirname + txt->path, txt->hdr);
                if (!txt->ldr.empty())
                    ygl::save_image4b(dirname + txt->path, txt->ldr);
            }
        } else if (sname == "shapes") {
            auto scn = test_scns.front().second;
            for (auto sgr : scn->shapes) {
                auto sscn = new ygl::scene();
                auto ssgr = new ygl::shape_group();
                ssgr->name = sgr->name;
                for (auto shp : sgr->shapes) {
                    auto sshp = new ygl::shape(*shp);
                    sshp->materials.clear();
                    ssgr->shapes.push_back(sshp);
                }
                sscn->shapes.push_back(ssgr);
                auto shp_name = ygl::partition(sgr->name, "_")[0];
                auto opts = ygl::save_options();
                ygl::save_scene(dirname + shp_name + ".obj", sscn, opts);
                delete sscn;
            }
        } else {
            auto count = 0;
            for (auto& test_scn : test_scns) {
                auto name = "scene"s;
                if (count) name = "view"s + std::to_string(count);
                count += 1;
                save_scene(test_scn.second, name, dirname,
                    !ygl::startswith(sname, "instance"));
                ygl::save_proc_scene(dirname + name + ".json", test_scn.first);
            }
        }
    } catch (std::exception& e) { ygl::log_fatal("error {}", e.what()); }
}

int main(int argc, char* argv[]) {
    // put together scene names
    auto scene_names = std::vector<std::string>{"cornell_box"};
    for (auto preset : ygl::proc_split_scene_presets())
        scene_names.push_back(preset->scn->name);

    // command line params
    auto parser = ygl::make_parser(argc, argv, "ytestgen", "make tests");
    auto scene = ygl::parse_opt(
        parser, "--scene", "-s", "scene name", ""s, false, scene_names);
    auto clean = ygl::parse_flag(parser, "--clean", "-c", "clean directory");
    auto dirname =
        ygl::parse_opt(parser, "--dirname", "-d", "directory name", "tests"s);
    auto no_parallel =
        ygl::parse_flag(parser, "--no-parallel", "", "do not run in parallel");
    auto quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (quiet) ygl::get_default_logger()->verbose = false;

    // make directories
    if (clean) {
        rmdir(dirname);
        ygl::log_info("cleaning directory {}", dirname);
    }
    ygl::log_info("creating directory {}", dirname);
    mkdir(dirname);

    if (scene != "") {
        save_test_scene(scene, dirname);
    } else if (no_parallel) {
        for (auto scn : scene_names) save_test_scene(scn, dirname);
    } else {
        auto threads = std::vector<std::thread>();
        for (auto scene_name : scene_names) {
            threads.push_back(std::thread([scene_name, dirname]() {
                save_test_scene(scene_name, dirname);
            }));
        }
        for (auto& t : threads) t.join();
    }
}
