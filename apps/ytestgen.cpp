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

#include <thread>
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_utils.h"
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

void save_scene(const ygl::scene* scn, const std::string& sname,
    const std::string& dirname, bool flatten_obj) {
    try {
        auto gltf = true;
        for (auto shp : scn->shapes) {
            if(!shp->quads_pos.empty()) gltf = false;
            if(!shp->beziers.empty()) gltf = false;
        }

        if (gltf) save_scene(dirname + sname + ".gltf", scn, true);
        ygl::save_scene(dirname + sname + ".obj", scn, true, !flatten_obj);
    } catch (std::exception& e) { ygl::log_fatal("error {}", e.what()); }
}

inline std::vector<std::string>& proc_scenes_types() {
    static auto names =
        std::vector<std::string>{"cornellbox", "textures", "shapes", "basic", "simple", "lines", "subdivs"};
    return names;
}

std::vector<ygl::scene*> make_proc_scenes(const std::string& name) {
    auto scenes = std::vector<ygl::scene*>();

    // cornell box
    if (name == "cornellbox") {
        scenes.push_back(
            ygl::make_cornellbox_scene("cornellbox_al", "arealight"));
        scenes.push_back(
            ygl::make_cornellbox_scene("cornellbox_pl", "pointlight"));
        scenes.push_back(
            ygl::make_cornellbox_scene("cornellbox_el", "envlight"));
    }

    // textures
    if (name == "textures") {
        for (auto txt : ygl::proc_texture_types()) {
            auto scn = ygl::make_simple_scene(
                "texture_" + txt, "quad", "matte", "arealight", "");
            ygl::log_error("fix texture scene");
            scenes.push_back(scn);
        }
    }

    // shapes
    if (name == "shapes") {
        for (auto shp : ygl::proc_shape_types()) {
            scenes.push_back(ygl::make_simple_scene(
                "shape_" + shp, shp, "plastic_colored", "arealight", ""));
        }
    }

    // light types
    auto all_lights = {"pointlight"s, "arealight"s, "envlight"s};
    auto area_lights = {"arealight"s, "envlight"s};
    auto lights_names = std::map<std::string, std::string>{
        {"pointlight", "pl"}, {"arealight", "al"}, {"envlight", "el"}
    };

    // basic shapes
    if (name == "basic") {
        for (auto lgt : all_lights) {
            scenes.push_back(ygl::make_simple_scene("basic_" + lights_names.at(lgt),
                {"sphere_flipcap", "sphere_cube", "cube_rounded"},
                {"plastic_red", "plastic_green", "plastic_blue"}, lgt, "matte"));
        }
    }

    // simple shapes
    if (name == "simple") {
        for (auto lgt : all_lights) {
            scenes.push_back(ygl::make_simple_scene("simple_" + lights_names.at(lgt),
                {"sphere_flipcap", "sphere_cube", "cube_rounded"},
                {"plastic_colored", "plastic_colored", "plastic_colored"},
                lgt));
        }
    }

    // lines
    if (name == "lines") {
        for (auto lgt : area_lights) {
            scenes.push_back(ygl::make_simple_scene("lines_" + lights_names.at(lgt),
                {"hairball_noise", "hairball_clump", "hairball"},
                {"matte_gray", "matte_gray", "matte_gray"}, lgt));
        }
    }

    // subdivs
    if (name == "subdivs") {
        for (auto lgt : area_lights) {
            scenes.push_back(ygl::make_simple_scene("subdivs_" + lights_names.at(lgt),
                {"cube_subdiv", "suzanne_subdiv", "fvcube_subdiv"},
                {"plastic_red", "plastic_green", "plastic_colored"}, lgt));
        }
    }

#if 0
    // transparent shapes
    presets["transparent_al"] = ygl::make_simple_scene({"quad", "quad", "quad"},
        {"transparent_red", "transparent_green", "transparent_blue"},
        "arealight");

    // plastics shapes
    presets["plastics_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "arealight");
    presets["plastics_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "envlight");

    // metals shapes
    presets["metals_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"metal_gold_rough", "metal_gold_sharp", "metal_silver_mirror"},
            "arealight");
    presets["metals_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"metal_gold_rough", "metal_gold_sharp", "metal_silver_mirror"},
            "envlight");

    // glass shapes
    presets["glass_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"glass_mirror", "glass_colored", "glass_rough"}, "arealight");
    presets["glass_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"glass_mirror", "glass_colored", "glass_rough"}, "envlight");

    // car paints shapes
    presets["paints_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"carpaint_black", "carpaint_blue", "metal_blue"}, "arealight");
    presets["paints_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"carpaint_black", "carpaint_blue", "metal_blue"}, "envlight");

    // lihght testing
    presets["arealights_nl"] =
        ygl::make_simple_scene({"quad", "sphere", "quad_stack"},
            {"emission", "emission", "emission"}, "nolights");
    presets["arealight_nl"] =
        ygl::make_simple_scene({"quad"}, {"emission"}, "nolights");

    // matball
    // presets["mattball_al"] = ygl::make_simple_scene( {"matball"},
    // {"carpaint_black"}, "arealight", ""); last = add_simple_scene(presets,
    // "mattball_el", {"matball"}, {"carpaint_black"}, "envlight", "");
    //    presets["matball_el"] =
    //        ygl::make_simple_scene({"sphere"}, {"transparent_none"},
    //        "envlight",
    //        "");
    //    presets["matball_al"] =
    //        ygl::make_simple_scene({"sphere"}, {"transparent_none"},
    //        "arealight",
    //        "");

    // tesselation shapes
    //    presets["tesselation_pl"] = ygl::make_simple_scene(
    //        {"geodesic_spherel", "geodesic_spheref", "geodesic_sphere"},
    //        {"matte_gray", "matte_gray", "matte_gray"}, "pointlight");

    // textureuv shapes
    presets["textureuv_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"matte_green", "matte_colored", "matte_uv"}, "pointlight");

    // normalmap shapes
    presets["normalmap_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"plastic_blue", "plastic_blue_grid_norm", "plastic_colored_bump_norm"},
        "pointlight");

    // animated shapes
    // presets["animated_pl"] = ygl::make_simple_scene(
    //     {"sphere_flipcap",
    //         "sphere_cube",
    //         "cube_rounded"},
    //     {"plastic_colored",
    //         "plastic_colored",
    //         "plastic_colored"},
    //     "pointlight", true, {"scale", "bounce", "rotation"});

    // instances
    presets["instances_pl"] =
        ygl::make_random_instances_scene({10, 10}, {{-3, -3, 0}, {3, 3, 0}});
    presets["instancel_pl"] =
        ygl::make_random_instances_scene({100, 100}, {{-3, -3, 0}, {3, 3, 0}});
#endif

    return scenes;
}

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ygl::make_parser(argc, argv, "ytestgen", "make tests");
    auto sname = ygl::parse_opt(parser, "--scene", "-s", "scenes name", ""s,
        false, proc_scenes_types());
    auto noclean = ygl::parse_flag(parser, "--noclean", "", "do not clean directory");
    auto dirname =
        ygl::parse_opt(parser, "--dirname", "-d", "directory name", "tests"s);
    auto noparallel =
        ygl::parse_flag(parser, "--noparallel", "", "do not run in parallel");
    auto quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // scenes
    auto snames =
        (sname != "") ? std::vector<std::string>{sname} : proc_scenes_types();

    // setup logger
    if (quiet) ygl::log_verbose() = false;

    // make directories
    if (!noclean) {
        rmdir(dirname);
        ygl::log_info("cleaning directory {}", dirname);
    }
    ygl::log_info("creating directory {}", dirname);
    mkdir(dirname);

    auto sgroups = std::map<std::string, std::vector<ygl::scene*>>();
    for (auto sname : snames) {
        ygl::log_info("generating scenes {}", sname);
        sgroups[sname] = make_proc_scenes(sname);
    }
    auto threads = std::vector<std::thread>();
    for(auto sname : snames) {
        auto scns = sgroups[sname];
        ygl::log_info("creating directory {}", dirname + "/" + sname);
        mkdir(dirname + "/" + sname);
        for (auto scn : scns) {
            if (noparallel) {
                ygl::log_info("saving scene {}/{}", sname, scn->name);
                save_scene(scn, scn->name, dirname + "/" + sname + "/",
                    sname == "instances");
            } else {
                threads.push_back(std::thread([=]() {
                    ygl::log_info("saving scene {}/{}", sname, scn->name);
                    save_scene(scn, scn->name, dirname + "/" + sname + "/",
                        sname == "instances");
                }));
            }
        }
    }
    for (auto& t : threads) t.join();
    for(auto sname : snames) {
        auto scns = sgroups[sname];
        for(auto scn : scns) delete scn;
    }
}
