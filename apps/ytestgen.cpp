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
    const std::string& dirname, bool preserve_instances) {
    try {
        save_scene(dirname + sname + ".gltf", scn, true);
        ygl::save_scene(
            dirname + sname + ".obj", scn, true, preserve_instances);
    } catch (std::exception& e) { ygl::log_fatal("error {}", e.what()); }
}

inline std::vector<std::string>& proc_scenes_types() {
    static auto names = std::vector<std::string>{"cornellbox", "textures",
        "shapes", "environments", "basic", "simple", "simplemat", "lines",
        "subdivs", "materials", "lighttests", "transparent", "instances",
        "animated", "examples"};
    return names;
}

std::vector<ygl::scene*> make_proc_scenes(const std::string& name) {
    auto scenes = std::vector<ygl::scene*>();

    auto nists = std::vector<ygl::instance*>{};

    // cornell box
    if (name == "cornellbox") {
        scenes.push_back(ygl::make_cornellbox_scene("cornellbox_al", false));
        scenes.push_back(ygl::make_cornellbox_scene("cornellbox_el", true));
    }

    // textures
    if (name == "textures") {
        auto txts = std::vector<ygl::texture*>{
            ygl::make_grid_texture("grid"),
            ygl::make_grid_texture("uvgrid"),
        };
        for (auto txt : txts) {
            scenes.push_back(ygl::make_simple_scene(txt->name,
                ygl::make_quad_shape("obj"),
                ygl::make_matte_material(txt->name, {1, 1, 1}, txt), false));
        }
    }

    // shapes
    if (name == "shapes") {
        auto shps = std::vector<ygl::shape*>{
            ygl::make_cube_shape("cube"),
            ygl::make_sphere_shape("sphere"),
            ygl::make_sphereflipcap_shape("sphereflipcap"),
            ygl::make_spherecube_shape("spherecube"),
            ygl::make_quad_shape("quad"),
        };
        for (auto shp : shps) {
            scenes.push_back(ygl::make_simple_scene(shp->name, shp,
                ygl::make_plastic_material(
                    "mat", {1, 1, 1}, 0.1f, ygl::make_uvgrid_texture("uvgrid")),
                false));
        }
    }

    // environments
    if (name == "environments") {
        auto envs = std::vector<ygl::environment*>{
            ygl::make_environment("const-env", {1, 1, 1}),
            ygl::make_environment(
                "sky-env", {1, 1, 1}, ygl::make_sky_texture("sky"))};
        for (auto env : envs) {
            scenes.push_back(ygl::make_environment_scene(env->name, env));
        }
    }

    // basic shapes
    if (name == "basic") {
        for (auto env : {0, 1}) {
            scenes.push_back(ygl::make_simple_scene(
                "basic"s + ((env) ? "-el" : "-al"),
                {
                    ygl::make_instance("obj1",
                        ygl::make_sphereflipcap_shape("obj1"),
                        ygl::make_plastic_material("obj1", {0.7f, 0.5f, 0.5f})),
                    ygl::make_instance("obj2",
                        ygl::make_spherecube_shape("obj2"),
                        ygl::make_plastic_material("obj2", {0.5f, 0.7f, 0.5f})),
                    ygl::make_instance("obj3",
                        ygl::make_cuberounded_shape("obj3"),
                        ygl::make_plastic_material("obj3", {0.5f, 0.5f, 0.7f})),
                },
                env));
        }
    }

    // simple shapes
    if (name == "simple") {
        for (auto env : {0, 1}) {
            scenes.push_back(
                ygl::make_simple_scene("simple"s + ((env) ? "-el" : "-al"),
                    {
                        ygl::make_instance("obj1",
                            ygl::make_sphereflipcap_shape("obj1"),
                            ygl::make_plastic_material("obj1", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                        ygl::make_instance("obj2",
                            ygl::make_spherecube_shape("obj2"),
                            ygl::make_plastic_material("obj2", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                        ygl::make_instance("obj3",
                            ygl::make_cuberounded_shape("obj3"),
                            ygl::make_plastic_material("obj3", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                    },
                    env));
        }
    }

    // simple materials
    if (name == "simplemat") {
        for (auto env : {0, 1}) {
            scenes.push_back(ygl::make_simple_scene(
                "simplemat"s + ((env) ? "-el" : "-al"),
                {
                    ygl::make_instance("obj1",
                        ygl::make_spherecube_shape("obj1"),
                        ygl::make_matte_material("obj1", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("obj2",
                        ygl::make_spherecube_shape("obj2"),
                        ygl::make_metal_material(
                            "obj2", {0.7f, 0.7f, 0.7f}, 0)),
                    ygl::make_instance("obj3",
                        ygl::make_spherecube_shape("obj3"),
                        ygl::make_plastic_material(
                            "obj3", {0.5f, 0.5f, 0.7f}, 0.01f)),
                },
                env));
        }
    }

    // lines
    if (name == "lines") {
        for (auto env : {0, 1}) {
            scenes.push_back(ygl::make_simple_scene(
                "lines"s + ((env) ? "-el" : "-al"),
                {
                    ygl::make_instance("obj1",
                        ygl::make_hairball_shape("obj1", 16, 2, 2,
                            {0.15f, 0.2f}, {0.05f, 100}, {0, 0}),
                        ygl::make_matte_material("obj2", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("obj1",
                        ygl::make_hairball_shape(
                            "obj2", 16, 2, 2, {0.15f, 0.2f}, {0.5f, 8}, {0, 0}),
                        ygl::make_matte_material("obj3", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("obj1",
                        ygl::make_hairball_shape("obj3", 16, 2, 2,
                            {0.15f, 0.2f}, {0, 0}, {0.5f, 128}),
                        ygl::make_matte_material("obj3", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("int1",
                        ygl::make_sphere_shape("int1", 4, 1.6f),
                        ygl::make_matte_material("int1", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("int2",
                        ygl::make_sphere_shape("int2", 4, 1.6f),
                        ygl::make_matte_material("int2", {0.7f, 0.7f, 0.7f})),
                    ygl::make_instance("int3",
                        ygl::make_sphere_shape("int3", 4, 1.6f),
                        ygl::make_matte_material("int3", {0.7f, 0.7f, 0.7f})),
                },
                env));
        }
    }

    // subdivs
    if (name == "subdivs") {
        for (auto env : {0, 1}) {
            scenes.push_back(
                ygl::make_simple_scene("subdivs"s + ((env) ? "-el" : "-al"),
                    {
                        ygl::make_instance("obj1",
                            ygl::make_cube_subdiv_shape("obj1"),
                            ygl::make_plastic_material("obj1", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid")),
                            ygl::make_cube_subdiv("obj1")),
                        ygl::make_instance("obj1",
                            ygl::make_suzanne_subdiv_shape("obj2"),
                            ygl::make_plastic_material("obj2", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid")),
                            ygl::make_suzanne_subdiv("obj2")),
                        ygl::make_instance("obj1",
                            ygl::make_fvcube_subdiv_shape("obj3"),
                            ygl::make_plastic_material("obj3", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid")),
                            ygl::make_fvcube_subdiv("obj3")),
                    },
                    env));
        }
        for (auto env : {0, 1}) {
            scenes.push_back(
                ygl::make_simple_scene("tsubdivs"s + ((env) ? "-el" : "-al"),
                    {
                        ygl::make_instance("obj1",
                            ygl::make_cube_subdiv_shape("obj1"),
                            ygl::make_plastic_material("obj1", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                        ygl::make_instance("obj2",
                            ygl::make_suzanne_subdiv_shape("obj2"),
                            ygl::make_plastic_material("obj2", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                        ygl::make_instance("obj3",
                            ygl::make_fvcube_subdiv_shape("obj3"),
                            ygl::make_plastic_material("obj3", {1, 1, 1}, 0.1f,
                                ygl::make_uvgrid_texture("uvgrid"))),
                    },
                    env));
        }
    }

    // materials
    if (name == "materials") {
        auto make_matgroups = []() {
            return std::vector<std::vector<ygl::material*>>{
                {
                    ygl::make_plastic_material(
                        "plastic1", {0.7f, 0.5f, 0.5f}, 0.1f),
                    ygl::make_plastic_material(
                        "plastic2", {0.5f, 0.7f, 0.5f}, 0.05f),
                    ygl::make_plastic_material(
                        "plastic3", {0.5f, 0.5f, 0.7f}, 0),
                },
                {
                    ygl::make_metal_material(
                        "metal1", {0.66f, 0.45f, 0.34f}, 0.1f),
                    ygl::make_metal_material(
                        "metal2", {0.66f, 0.45f, 0.34f}, 0.05f),
                    ygl::make_metal_material("metal3", {0.7f, 0.7f, 0.7f}, 0),
                },
                {
                    ygl::make_glass_material("glass1", {1, 1, 1}, 0),
                    ygl::make_glass_material("glass2", {1, 0.7f, 0.7f}, 0),
                    ygl::make_glass_material("glass3", {1, 1, 1}, 0.1f),
                },
                {
                    ygl::make_solidglass_material("refractive1", {1, 1, 1}, 0),
                    ygl::make_solidglass_material(
                        "refractive2", {1, 0.7f, 0.7f}, 0),
                    ygl::make_solidglass_material(
                        "refractive3", {1, 1, 1}, 0.1f),
                },
                {
                    ygl::make_transparent_material(
                        "transparent1", {0.5f, 0.2f, 0.2f}, 0.9f),
                    ygl::make_transparent_material(
                        "transparent2", {0.2f, 0.5f, 0.2f}, 0.5f),
                    ygl::make_transparent_material(
                        "transparent3", {0.2f, 0.2f, 0.5f}, 0.2f),
                },
                {
                    ygl::make_plastic_material("bump1", {0.2f, 0.5f, 0.2f},
                        0.05f, nullptr, ygl::make_bumpnorm_texture("norm1")),
                    ygl::make_plastic_material("bump2", {0.2f, 0.5f, 0.2f},
                        0.05f, nullptr, ygl::make_bumpnorm_texture("norm2")),
                    ygl::make_plastic_material("bump3", {0.2f, 0.5f, 0.2f},
                        0.05f, nullptr, ygl::make_bumpnorm_texture("norm3")),
                },
            };
        };
        for (auto env : {0, 1}) {
            for (auto& mats : make_matgroups()) {
                for (auto mat : mats) {
                    scenes.push_back(ygl::make_simple_scene(
                        mat->name + ((env) ? "-el" : "-al"),
                        ygl::make_matball_shape("obj"), mat, env));
                }
            }
        }
        for (auto env : {0, 1}) {
            for (auto& mats : make_matgroups()) {
                auto name = mats[0]->name;
                name.back() = 's';
                scenes.push_back(
                    ygl::make_simple_scene(name + ((env) ? "-el" : "-al"),
                        {
                            ygl::make_instance("obj1",
                                ygl::make_sphere_shape("obj1"), mats[0]),
                            ygl::make_instance("obj2",
                                ygl::make_sphere_shape("obj2"), mats[1]),
                            ygl::make_instance("obj3",
                                ygl::make_sphere_shape("obj3"), mats[2]),
                        },
                        env));
            }
        }
    }

    // light tests
    if (name == "lighttests") {
        auto shps = std::vector<ygl::shape*>{
            ygl::make_quad_shape("quad"),
            ygl::make_sphere_shape("sphere"),
            ygl::make_quadstack_shape("quadstack"),
        };
        for (auto shp : shps) {
            scenes.push_back(ygl::make_simple_scene(shp->name + "light", shp,
                ygl::make_emission_material("light", {1, 1, 1}), true));
        }
    }

    // instances
    if (name == "instances") {
        scenes.push_back(ygl::make_random_instances_scene(
            "instances-small", {10, 10}, {{-3, -3, 0}, {3, 3, 0}}));
        scenes.push_back(ygl::make_random_instances_scene(
            "instances-large", {100, 100}, {{-3, -3, 0}, {3, 3, 0}}));
    }

    // instances
    if (name == "animated") {
        scenes.push_back(ygl::make_simple_scene("animated-al",
            {
                ygl::make_instance("obj1",
                    ygl::make_sphereflipcap_shape("obj1"),
                    ygl::make_plastic_material("obj1", {0.5f, 0.2f, 0.2f})),
                ygl::make_instance("obj2", ygl::make_spherecube_shape("obj2"),
                    ygl::make_plastic_material("obj2", {0.2f, 0.5f, 0.2f})),
                ygl::make_instance("obj3", ygl::make_cuberounded_shape("obj3"),
                    ygl::make_plastic_material("obj3", {0.2f, 0.2f, 0.5f})),
            },
            false,
            {
                ygl::make_animation("obj1", {0, 1, 2}, {}, {},
                    {{1, 1, 1}, {0.7f, 0.7f, 0.7f}, {1, 1, 1}}),
                ygl::make_animation(
                    "obj2", {0, 1, 2}, {{0, 1, 0}, {0, 2, 0}, {0, 1, 0}}),
                ygl::make_animation("obj3", {0, 1, 2}, {},
                    {ygl::rotation_quat({0, 1, 0}, 0),
                        ygl::rotation_quat({0, 1, 0}, ygl::pi),
                        ygl::rotation_quat({0, 1, 0}, 0)}),
            }));
    }

    return scenes;
}

int main(int argc, char* argv[]) {
    // command line params
    auto parser = ygl::make_parser(argc, argv, "ytestgen", "make tests");
    auto sname = ygl::parse_opt(parser, "--scene", "-s", "scenes name", ""s,
        false, proc_scenes_types());
    auto noclean =
        ygl::parse_flag(parser, "--noclean", "", "do not clean directory");
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

    auto threads = std::vector<std::thread>();
    auto sgroups = std::map<std::string, std::vector<ygl::scene*>>();
    for (auto sname : snames) {
        sgroups[sname] = {};
        if (noparallel) {
            ygl::log_info("generating scenes {}", sname);
            sgroups[sname] = make_proc_scenes(sname);
        } else {
            threads.push_back(std::thread([=, &sgroups]() {
                ygl::log_info("start generating scenes {}", sname);
                sgroups[sname] = make_proc_scenes(sname);
                ygl::log_info("end generating scenes {}", sname);
            }));
        }
    }
    for (auto& t : threads) t.join();
    threads.clear();
    for (auto sname : snames) {
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
                    ygl::log_info("start saving scene {}/{}", sname, scn->name);
                    save_scene(scn, scn->name, dirname + "/" + sname + "/",
                        sname == "instances");
                    ygl::log_info("end saving scene {}/{}", sname, scn->name);
                }));
            }
        }
    }
    for (auto& t : threads) t.join();
    for (auto sname : snames) {
        auto scns = sgroups[sname];
        for (auto scn : scns) delete scn;
    }
}
