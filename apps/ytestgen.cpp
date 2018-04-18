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
#include "../yocto/yocto_utils.h"
#include <thread>
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
    auto facevarying = false;
    for (auto shp : scn->shapes)
        facevarying = facevarying || !shp->quads_pos.empty();

    if (!facevarying) save_scene(dirname + sname + ".gltf", scn, true);
    ygl::save_scene(dirname + sname + ".obj", scn, true, !flatten_obj);
}

void save_test_scene(const ygl::scene* scn, const std::string& sname,
    const std::string& basedir) {
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };

    auto dirname = basedir + "/" + sname + "/";
    try {
        mkdir(dirname);
        ygl::log_info("saving scene {}", sname);
        if (sname == "textures") {
            for (auto txt : scn->textures) {
                if (!txt->hdr.pixels.empty())
                    ygl::save_image4f(dirname + txt->path, txt->hdr);
                if (!txt->ldr.pixels.empty())
                    ygl::save_image4b(dirname + txt->path, txt->ldr);
            }
        } else if (sname == "shapes") {
            for (auto shp : scn->shapes) {
                auto sscn = new ygl::scene();
                auto sshp = new ygl::shape(*shp);
                sscn->shapes.push_back(sshp);
                auto shp_name = sshp->name;
                ygl::save_scene(dirname + shp_name + ".obj", sscn);
                delete sscn;
            }
        } else {
            save_scene(scn, "scene", dirname, !startswith(sname, "instance"));
        }
    } catch (std::exception& e) { ygl::log_fatal("error {}", e.what()); }
}

std::map<std::string, ygl::scene*> make_scene_presets() {
    auto presets = std::map<std::string, ygl::scene*>();

    // cornell box
    presets["cornell_box"] = ygl::make_cornell_box_scene();

    // textures
    presets["textures"] = new ygl::scene();
    for (auto txt : ygl::proc_texture_types()) {
        presets["textures"]->textures.push_back(
            ygl::make_proc_texture(txt, txt, 512));
    }

    // shapes
    presets["shapes"] = new ygl::scene();
    for (auto shp : ygl::proc_shape_types()) {
        presets["shapes"]->shapes.push_back(ygl::make_proc_shape(shp, shp));
    }

    // basic shapes
    presets["basic_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_red", "plastic_green", "plastic_blue"}, "pointlights");
    presets["basic_al"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_red", "plastic_green", "plastic_blue"}, "arealights");
    presets["basic_el"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_red", "plastic_green", "plastic_blue"}, "envlights");

    // simple shapes
    presets["simple_al"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_colored", "plastic_colored", "plastic_colored"},
        "arealights");
    presets["simple_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_colored", "plastic_colored", "plastic_colored"},
        "pointlights");
    presets["simple_el"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_colored", "plastic_colored", "plastic_colored"}, "envlights");

    // simple shapes 1
    presets["spheres_al"] =
        ygl::make_simple_scene({"sphere_flipcap", "sphere_cube", "sphere"},
            {"plastic_colored", "plastic_colored", "plastic_colored"},
            "arealights");

    // simple shapes 2
    presets["cubes_al"] =
        ygl::make_simple_scene({"quad", "cube", "cube_rounded"},
            {"plastic_colored", "plastic_colored", "plastic_colored"},
            "arealights");

    // simple shapes 3
    presets["cylinders_al"] =
        ygl::make_simple_scene({"disk", "cylinder", "cylinder_rounded"},
            {"plastic_colored", "plastic_colored", "plastic_colored"},
            "arealights");

    // simple shapes 3
    presets["disks_al"] =
        ygl::make_simple_scene({"disk", "disk_quad", "disk_bulged"},
            {"plastic_colored", "plastic_colored", "plastic_colored"},
            "arealights");

    // transparent shapes
    presets["transparent_al"] = ygl::make_simple_scene({"quad", "quad", "quad"},
        {"transparent_red", "transparent_green", "transparent_blue"},
        "arealights");

    // lines shapes
    presets["lines_al"] =
        ygl::make_simple_scene({"hairball_noise", "hairball_clump", "hairball"},
            {"matte_gray", "matte_gray", "matte_gray"}, "arealights");

    // subdiv shapes
    presets["subdiv_al"] = ygl::make_simple_scene(
        {"cube_subdiv", "suzanne_subdiv", "fvcube_subdiv"},
        {"plastic_red", "plastic_green", "plastic_colored"}, "arealights");

    // plastics shapes
    presets["plastics_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "arealights");
    presets["plastics_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "envlights");

    // metals shapes
    presets["metals_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"metal_gold_rough", "metal_gold_sharp", "metal_silver_mirror"},
            "arealights");
    presets["metals_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"metal_gold_rough", "metal_gold_sharp", "metal_silver_mirror"},
            "envlights");

    // glass shapes
    presets["glass_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"glass_mirror", "glass_colored", "glass_rough"}, "arealights");
    presets["glass_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"glass_mirror", "glass_colored", "glass_rough"}, "envlights");

    // car paints shapes
    presets["paints_al"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"carpaint_black", "carpaint_blue", "metal_blue"}, "arealights");
    presets["paints_el"] =
        ygl::make_simple_scene({"matball", "matball", "matball"},
            {"carpaint_black", "carpaint_blue", "metal_blue"}, "envlights");

    // matball
    // presets["mattball_al"] = ygl::make_simple_scene( {"matball"},
    // {"carpaint_black"}, "arealights", ""); last = add_simple_scene(presets,
    // "mattball_el", {"matball"}, {"carpaint_black"}, "envlights", "");
    //    presets["matball_el"] =
    //        ygl::make_simple_scene({"sphere"}, {"transparent_none"},
    //        "envlights",
    //        "");
    //    presets["matball_al"] =
    //        ygl::make_simple_scene({"sphere"}, {"transparent_none"},
    //        "arealights",
    //        "");

    // tesselation shapes
    //    presets["tesselation_pl"] = ygl::make_simple_scene(
    //        {"geodesic_spherel", "geodesic_spheref", "geodesic_sphere"},
    //        {"matte_gray", "matte_gray", "matte_gray"}, "pointlights");

    // textureuv shapes
    presets["textureuv_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"matte_green", "matte_colored", "matte_uv"}, "pointlights");

    // normalmap shapes
    presets["normalmap_pl"] = ygl::make_simple_scene(
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"plastic_blue", "plastic_blue_grid_norm", "plastic_colored_bump_norm"},
        "pointlights");

    // animated shapes
    // presets["animated_pl"] = ygl::make_simple_scene(
    //     {"sphere_flipcap",
    //         "sphere_cube",
    //         "cube_rounded"},
    //     {"plastic_colored",
    //         "plastic_colored",
    //         "plastic_colored"},
    //     "pointlights", true, {"scale", "bounce", "rotation"});

    // instances
    presets["instances_pl"] =
        ygl::make_random_scene({10, 10}, {{-3, -3, 0}, {3, 3, 0}});
    presets["instancel_pl"] =
        ygl::make_random_scene({100, 100}, {{-3, -3, 0}, {3, 3, 0}});

    return presets;
}

int main(int argc, char* argv[]) {
    // put together scene names
    auto presets = make_scene_presets();
    auto scene_names = std::vector<std::string>{};
    for (auto preset : presets) scene_names.push_back(preset.first);

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
        save_test_scene(presets.at(scene), scene, dirname);
    } else if (no_parallel) {
        for (auto& scn_kv : presets)
            save_test_scene(scn_kv.second, scn_kv.first, dirname);
    } else {
        auto threads = std::vector<std::thread>();
        for (auto& scn_kv : presets) {
            threads.push_back(std::thread([scn_kv, dirname]() {
                save_test_scene(scn_kv.second, scn_kv.first, dirname);
            }));
        }
        for (auto& t : threads) t.join();
    }
}
