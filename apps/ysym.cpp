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

#include "yscene.h"

// ---------------------------------------------------------------------------
// INTERACTIVE FUNCTIONS
// ---------------------------------------------------------------------------

void simulate_step(yscene* scene);
void simulate_reset(yscene* scene);

#ifndef YOCTO_NO_OPENGL

bool update(yscene* scn) {
    // advance if simulating
    if (scn->animate) {
        simulate_step(scn);
        scn->time += scn->simulation_params.dt;
        scn->time_range.y += scn->simulation_params.dt;
        return true;
    }
    return false;
}

void draw_custom_widgets(yscene* scn) {}

#endif

// ---------------------------------------------------------------------------
// SIMULATION SCENE
// ---------------------------------------------------------------------------

ysym::scene* make_simulation_scene(const yobj::scene* scene) {
    // allocate scene
    auto simulation_scene = ysym::make_scene();

    // add each material
    auto material_map = std::map<yobj::material*, int>{{nullptr, -1}};
    for (auto mat : scene->materials) {
        auto simulated = (mat->name != "floor" && mat->name != "floor_txt" &&
                          mat->ke == ym::zero3f);
        auto density = (simulated) ? 1.0f : 0.0f;
        material_map[mat] = ysym::add_rigid_material(simulation_scene, density);
    }

    // add each shape
    auto shape_map = std::map<yobj::shape*, int>{{nullptr, -1}};
    for (auto msh : scene->meshes) {
        for (auto shp : msh->shapes) {
            shape_map[shp] = ysym::add_rigid_shape(simulation_scene,
                (int)shp->triangles.size(), (ym::vec3i*)shp->triangles.data(),
                (int)shp->pos.size(), (ym::vec3f*)shp->pos.data());
        }
    }

    // add each instance (support only one shape per instance)
    for (auto ist : scene->instances) {
        auto shp = ist->msh->shapes[0];
        ysym::add_rigid_body(simulation_scene, to_frame(ym::mat4f(ist->xform)),
            shape_map.at(shp), material_map.at(shp->mat), {0, 0, 0}, {0, 0, 0});
    }

    // initialize
    ysym::init_simulation(simulation_scene);

    return simulation_scene;
}

ysym::scene* make_simulation_scene(const ygltf::scene_group* scene) {
    // allocate scene
    auto simulation_scene = ysym::make_scene();

    // add each material
    auto material_map = std::map<ygltf::material*, int>{{nullptr, -1}};
    for (auto mat : scene->materials) {
        auto simulated = (mat->name != "floor" && mat->name != "floor_txt" &&
                          mat->emission == ym::zero3f);
        auto density = (simulated) ? 1.0f : 0.0f;
        material_map[mat] = ysym::add_rigid_material(simulation_scene, density);
    }

    // add each shape
    auto shape_map = std::map<ygltf::shape*, int>{{nullptr, -1}};
    for (auto msh : scene->meshes) {
        for (auto shp : msh->shapes) {
            shape_map[shp] = ysym::add_rigid_shape(simulation_scene,
                (int)shp->triangles.size(), (ym::vec3i*)shp->triangles.data(),
                (int)shp->pos.size(), (ym::vec3f*)shp->pos.data());
        }
    }

    // add each instance (support only one shape per instance)
    auto instances = ygltf::get_mesh_nodes(scene->scenes[0]);
    for (auto ist : instances) {
        auto shp = ist->msh->shapes[0];
        ysym::add_rigid_body(simulation_scene, to_frame(ym::mat4f(ist->xform)),
            shape_map.at(shp), material_map.at(shp->mat), {0, 0, 0},
            {0, 0, 0});
    }

    // initialize
    ysym::init_simulation(simulation_scene);

    return simulation_scene;
}

// ---------------------------------------------------------------------------
// ADVANCE SIMULATION
// ---------------------------------------------------------------------------

void simulate_step(yscene* scn) {
    ysym::advance_simulation(scn->simulation_scene, scn->simulation_params);
    if (scn->oscn) {
        for (auto iid = 0; iid < scn->oscn->instances.size(); iid++) {
            scn->oscn->instances[iid]->xform = to_mat(ym::frame3f(
                ysym::get_rigid_body_frame(scn->simulation_scene, iid)));
        }
    } else if (scn->gscn) {
        auto instances = ygltf::get_mesh_nodes(scn->gscn->default_scene);
        for (auto iid = 0; iid < instances.size(); iid++) {
            instances[iid]->xform = to_mat(ym::frame3f(
                ysym::get_rigid_body_frame(scn->simulation_scene, iid)));
        }
    }
}

void simulate_reset(yscene* scn) {
    if (scn->oscn) {
        for (int iid = 0; iid < scn->oscn->instances.size(); iid++) {
            scn->oscn->instances[iid]->xform =
                to_mat(ym::frame3f(scn->simulation_initial_state[iid]));
            ysym::set_rigid_body_frame(
                scn->simulation_scene, iid, scn->simulation_initial_state[iid]);
        }
    } else if (scn->gscn) {
        auto instances = ygltf::get_mesh_nodes(scn->gscn->default_scene);
        for (int iid = 0; iid < instances.size(); iid++) {
            instances[iid]->xform =
                to_mat(ym::frame3f(scn->simulation_initial_state[iid]));
            ysym::set_rigid_body_frame(
                scn->simulation_scene, iid, scn->simulation_initial_state[iid]);
        }
    }
}

// ---------------------------------------------------------------------------
// OFFLINE SIMULATION
// ---------------------------------------------------------------------------

void simulate_offline(yscene* scene) {
    // simulate each frame and save the results to a new scene
    printf("rigid body simulation for %s to %s\n", scene->filename.c_str(),
        scene->outfilename.c_str());
    printf("simulating ...");
    for (auto i = 0; i < scene->simulation_nframes; i++) {
        printf("\rsimulating frame %d/%d", i, scene->simulation_nframes);
        simulate_step(scene);
        std::string errmsg;
        char frame_filename[4096];
        sprintf(frame_filename, scene->outfilename.c_str(), i);
        if (scene->oscn) {
            yobj::save_scene(frame_filename, scene->oscn, false);
        } else if (scene->gscn) {
            ygltf::save_scenes(frame_filename, scene->gscn, false);
        }
    }
    printf("\rsimulating done\n");
}

// ---------------------------------------------------------------------------
// MAIN
// ---------------------------------------------------------------------------

int main(int argc, char* argv[]) {
    // create empty scene
    auto scn = new yscene();

    // command line
    parse_cmdline(
        scn, argc, argv, "rigid body simulation of a scene", true, false);

    // setting up rendering
    load_scene(scn, scn->filename, true, false);
    scn->simulation_scene = (scn->oscn) ? make_simulation_scene(scn->oscn) :
                                          make_simulation_scene(scn->gscn);

    // initialize simulation
    ysym::init_simulation(scn->simulation_scene);

    // initialize overlap
    ysym::init_overlap(scn->simulation_scene);

#ifndef YOCTO_NO_OPENGL
    // render offline or online
    if (!scn->interactive) {
        simulate_offline(scn);
    } else {
        // save init values
        if (scn->oscn) {
            scn->simulation_initial_state.resize(scn->oscn->instances.size());
            for (auto i = 0; i < scn->simulation_initial_state.size(); i++)
                scn->simulation_initial_state[i] =
                    ym::to_frame(scn->oscn->instances[i]->xform);

        } else if (scn->gscn) {
            auto instances = ygltf::get_mesh_nodes(scn->gscn->default_scene);
            scn->simulation_initial_state.resize(instances.size());
            for (auto i = 0; i < scn->simulation_initial_state.size(); i++)
                scn->simulation_initial_state[i] =
                    ym::to_frame(instances[i]->xform);
        }

        // run ui
        auto width = (int)std::round(scn->view_cam->aspect * scn->resolution);
        auto height = scn->resolution;
        run_ui(scn, width, height, "ysym", shade_init, shade_draw, update);
    }
#else
    simulate_offline(scn);
#endif

    // cleanup
    delete scn;

    // done
    return 0;
}
