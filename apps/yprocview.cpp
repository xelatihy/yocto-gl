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

#define YGL_OPENGL 1
#include "../yocto/yocto_gl.h"
using namespace ygl;

// Procedural camera
struct proc_camera {
    string name = "cam";
    vec3f from = {0, 0, -4};
    vec3f to = {0, 0, 0};
    vec3f up = {0, 1, 0};
    float yfov = 1;
    float aspect = 1;

    vector<camera*> cams;

    ~proc_camera() {
        for (auto e : cams) delete e;
    }
};

// Procedural shape type
enum struct proc_shape_type { none, floor, prim };

// Procedural floor shape params
struct proc_floor_shape_params {
    float size = 10;
    int tesselation = 5;
    bool textured = false;
};

// Procesural shape
struct proc_shape {
    string name = "shp";
    proc_shape_type type = proc_shape_type::none;
    
    proc_floor_shape_params floor_params = {};

    vector<shape*> shps;
    vector<material*> mats;

    ~proc_shape() {
        for (auto e : shps) delete e;
        for (auto e : mats) delete e;
    }
};

// Procedural instance type
enum struct proc_instance_type { single };

// Procedural instances
struct proc_instance {
    string name = "ist";
    proc_instance_type type = proc_instance_type::single;
    frame3f frame = identity_frame3f;
    proc_shape* shp = nullptr;

    vector<instance*> ists;

    ~proc_instance() {
        for (auto e : ists) delete e;
    }
};

// Procedural scene
struct proc_scene {
    vector<proc_camera*> cameras;
    vector<proc_shape*> shapes;
    vector<proc_instance*> instances;

    scene* scn = nullptr;
    
    ~proc_scene() {
        if(scn) {
            scn->cameras.clear();
            scn->shapes.clear();
            scn->materials.clear();
            scn->instances.clear();
            delete scn;
        }
        for(auto e : cameras) delete e;
        for(auto e : shapes) delete e;
        for(auto e : instances) delete e;
    }
};

// init a procedural scene with default objects
proc_scene* init_proc_scene() {
    auto pscn = new proc_scene();

    auto pcam = new proc_camera();
    pscn->cameras += pcam;
    
    auto pshp = new proc_shape();
    pscn->shapes += pshp;
    pshp->type = proc_shape_type::floor;

    auto pist = new proc_instance();
    pscn->instances += pist;
    pist->name = "floor";
    pist->shp = pshp;
    
    return pscn;
}

// add a certain number of objects
template <typename T>
bool add_proc_objects(vector<T*>& objs, int num) {
    if (objs.size() == num) return false;
    if (objs.size() < num) {
        for (auto i = 0; i < num - objs.size(); i++) objs += new T();
        return true;
    } else {
        throw runtime_error("removal not yet supported");
        return false;
    }
}

// Updates the proc camera
bool update_proc_camera(proc_camera* pcam) {
    auto num_changed = add_proc_objects(pcam->cams, 1);
    auto cam = pcam->cams.front();
    cam->name = pcam->name;
    cam->frame = lookat_frame3f(pcam->from, pcam->to, pcam->up);
    cam->yfov = pcam->yfov;
    cam->aspect = pcam->aspect;
    cam->near = 0.01f;
    cam->far = 10.0f;
    cam->aperture = 0;
    cam->focus = length(pcam->from - pcam->to);
    return num_changed;
}

// Updates for floor proc shape
bool update_proc_floor_shape(proc_shape* pshp) {
    auto num_changed_shp = add_proc_objects(pshp->shps, 1);
    auto num_changed_mat = add_proc_objects(pshp->mats, 1);
    
    auto& params = pshp->floor_params;
    auto shp = pshp->shps[0];
    auto mat = pshp->mats[0];

    shp->name = pshp->name;
    tie(shp->quads, shp->pos, shp->norm, shp->texcoord) = make_uvquad(params.tesselation);
    
    mat->name = pshp->name;
    mat->mtype = material_type::specular_roughness;
    mat->kd = {0.2f,0.2f,0.2f};
    mat->ks = zero3f;
    mat->rs = 1;
    
    return num_changed_shp || num_changed_mat;
}

// Updates for prim proc shape
bool update_proc_prim_shape(proc_shape* pshp) {
    auto num_changed_shp = add_proc_objects(pshp->shps, 1);
    auto num_changed_mat = add_proc_objects(pshp->mats, 1);
    return num_changed_shp || num_changed_mat;
}

// Updates the proc shape
bool update_proc_shape(proc_shape* pshp) {
    switch (pshp->type) {
        case proc_shape_type::none: {
            auto num_changed_shp = add_proc_objects(pshp->shps, 0);
            auto num_changed_mat = add_proc_objects(pshp->mats, 0);
            return num_changed_shp || num_changed_mat;
        } break;
        case proc_shape_type::floor: {
            return update_proc_floor_shape(pshp);
        } break;
        case proc_shape_type::prim: {
            return update_proc_prim_shape(pshp);
        } break;
        default: throw runtime_error("should not have entered here");
    }
}

// Updates the proc single instance
bool update_proc_single_instance(proc_instance* pist) {
    auto num_changed = add_proc_objects(pist->ists, pist->shp->shps.size());
    for (auto sid : range(pist->shp->shps.size())) {
        pist->ists[sid]->name = pist->name + "_" + to_string(sid);
        pist->ists[sid]->frame = pist->frame;
        pist->ists[sid]->shp = pist->shp->shps[sid];
    }
    return num_changed;
}

// Updates the proc instance
bool update_proc_instance(proc_instance* pist) {
    switch (pist->type) {
        case proc_instance_type::single: {
            return update_proc_single_instance(pist);
        } break;
        default: throw runtime_error("should not have entered here");
    }
}

// Updates the view scene from the procedural one
void update_proc_scene(proc_scene* pscn) {
    if (!pscn->scn) pscn->scn = new scene();
    auto changed_cams = false;
    for (auto pcam : pscn->cameras)
        if (update_proc_camera(pcam)) changed_cams = true;
    if (changed_cams) {
        pscn->scn->cameras.clear();
        for (auto pcam : pscn->cameras) pscn->scn->cameras += pcam->cams;
    }
    auto changed_shps = false;
    for (auto pshp : pscn->shapes)
        if (update_proc_shape(pshp)) changed_shps = true;
    if (changed_shps) {
        pscn->scn->shapes.clear();
        pscn->scn->materials.clear();
        for (auto pshp : pscn->shapes) pscn->scn->shapes += pshp->shps;
        for (auto pshp : pscn->shapes) pscn->scn->materials += pshp->mats;
    }
    auto changed_ists = false;
    for (auto pist : pscn->instances)
        if (update_proc_instance(pist)) changed_ists = true;
    if (changed_shps) {
        pscn->scn->instances.clear();
        for (auto pist : pscn->instances) pscn->scn->instances += pist->ists;
    }
}

// Application state
struct app_state {
    proc_scene* pscn = nullptr;
    scene* scn = nullptr;
    string filename;
    string imfilename;
    string outfilename;
    gl_stdsurface_params shparams = {};
    gl_stdsurface_state* shstate = nullptr;
    bool navigation_fps = false;
    void* selection = nullptr;

    ~app_state() {
        if (shstate) delete shstate;
        if (pscn) delete pscn;
        // scn is owned by pscn
    }
};

// draw with shading
inline void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    auto framebuffer_size = get_framebuffer_size(win);
    app->shparams.width = framebuffer_size.x;
    app->shparams.height = framebuffer_size.y;
    auto cam = app->scn->cameras[app->shparams.camera_id];
    cam->aspect = (float)framebuffer_size.x / (float)framebuffer_size.y;

    update_lights(app->scn, false, false);
    update_stdsurface_state(app->shstate, app->scn, app->shparams);
    if (app->shstate->lights_pos.empty()) app->shparams.camera_lights = true;

    gl_clear_buffers();
    draw_stdsurface_scene(app->shstate, app->scn, app->shparams);

    if (begin_widgets(win, "yprocview")) {
        draw_label_widget(win, "scene", app->filename);
        draw_camera_widget(win, "camera", app->scn, app->shparams.camera_id);
        draw_value_widget(win, "wire", app->shparams.wireframe);
        draw_continue_widget(win);
        draw_value_widget(win, "edges", app->shparams.edges);
        draw_continue_widget(win);
        draw_value_widget(win, "cutout", app->shparams.cutout);
        draw_continue_widget(win);
        draw_value_widget(win, "fps", app->navigation_fps);
        draw_tonemap_widgets(win, "", app->shparams.exposure,
            app->shparams.gamma, app->shparams.filmic);
        draw_scene_widgets(
            win, "scene", app->scn, app->selection, app->shstate->txt);
    }
    end_widgets(win);

    swap_buffers(win);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = make_window(app->shparams.width, app->shparams.height,
        "yprocview | " + app->filename, app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // load textures and vbos
    app->shstate = make_stdsurface_state();
    update_stdsurface_state(app->shstate, app->scn, app->shparams);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        // handle mouse and keyboard for navigation
        auto cam = app->scn->cameras[app->shparams.camera_id];
        handle_camera_navigation(win, cam, app->navigation_fps);

        // draw
        draw(win);

        // event hadling
        poll_events(win);
    }

    clear_window(win);
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = make_parser(argc, argv, "yview", "views scenes inteactively");
    app->shparams.exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->shparams.gamma =
        parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->shparams.filmic =
        parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    app->shparams.height =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    app->shparams.ambient = {amb, amb, amb};
    app->shparams.camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
    auto log_filename = parse_opt(parser, "--log", "", "log to disk", ""s);
    if (log_filename != "") add_file_stream(log_filename, true);
    auto preserve_quads =
        parse_flag(parser, "--preserve-quads", "-q", "preserve quads on load");
    auto preserve_facevarying = parse_flag(
        parser, "--preserve-facevarying", "-f", "preserve facevarying on load");
    app->imfilename =
        parse_opt(parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = parse_arg(parser, "scene", "scene filename", ""s);
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // scene loading
    // log_info("loading scene {}", app->filename);
    // DISABLED FOR TESTING
    
    // make a test scene
    app->pscn = init_proc_scene();
    
    // convert to scene
    update_proc_scene(app->pscn);
    
    // fix scene pointer
    app->scn = app->pscn->scn;

    // run ui
    auto cam = app->scn->cameras[app->shparams.camera_id];
    app->shparams.width = (int)round(cam->aspect * app->shparams.height);
    run_ui(app);

    // clear
    delete app;

    // done
    return 0;
}
