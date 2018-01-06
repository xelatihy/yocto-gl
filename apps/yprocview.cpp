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
    vec3f from = {0, 4, 10};
    vec3f to = {0, 1, 0};
    vec3f up = {0, 1, 0};
    float yfov = 15 * pif / 180;
    float aspect = 16.0f / 9.0f;

    vector<camera*> cams;

    ~proc_camera() {
        for (auto e : cams) delete e;
    }
};

// Procedural shape type
enum struct proc_shape_type { none, floor, prim };

// Names for enumeration
inline const vector<pair<string, proc_shape_type>>& proc_shape_names() {
    static auto names = vector<pair<string, proc_shape_type>>{
        {"none", proc_shape_type::none},
        {"floor", proc_shape_type::floor},
        {"prim", proc_shape_type::prim},
    };
    return names;
}

// Procedural floor shape params
struct proc_floor_shape_params {
    float size = 10;
    int tesselation = 5;
    vec3f color = {0.2f, 0.2f, 0.2f};
    bool textured = false;
};

// Procedural prim shape type
enum struct proc_prim_shape_type { sphere };

// Names for enumeration
inline const vector<pair<string, proc_prim_shape_type>>&
proc_prim_shape_names() {
    static auto names = vector<pair<string, proc_prim_shape_type>>{
        {"sphere", proc_prim_shape_type::sphere},
    };
    return names;
}

// Procedural prim shape params
struct proc_prim_shape_params {
    proc_prim_shape_type ptype = proc_prim_shape_type::sphere;
    float size = 1;
    int tesselation = 5;
    int subdivision = 0;
    vec3f color = {0.2f, 0.2f, 0.2f};
    bool textured = false;
};

// Procesural shape
struct proc_shape {
    string name = "shp";
    proc_shape_type type = proc_shape_type::none;

    proc_floor_shape_params floor_params = {};
    proc_prim_shape_params prim_params = {};

    vector<shape*> shps;
    vector<material*> mats;

    ~proc_shape() {
        for (auto e : shps) delete e;
        for (auto e : mats) delete e;
    }
};

// Procedural instance type
enum struct proc_instance_type { single };

// Names for enumeration
inline const vector<pair<string, proc_instance_type>>& proc_instance_names() {
    static auto names = vector<pair<string, proc_instance_type>>{
        {"single", proc_instance_type::single},
    };
    return names;
}

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
        if (scn) {
            scn->cameras.clear();
            scn->shapes.clear();
            scn->materials.clear();
            scn->instances.clear();
            delete scn;
        }
        for (auto e : cameras) delete e;
        for (auto e : shapes) delete e;
        for (auto e : instances) delete e;
    }
};

// init a procedural scene with default objects
proc_scene* init_proc_scene() {
    auto pscn = new proc_scene();

    auto pcam = new proc_camera();
    pscn->cameras += pcam;

    auto pshp = new proc_shape();
    pscn->shapes += pshp;
    pshp->name = "floor";
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
    cam->far = 100.0f;
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
    tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
        make_uvquad(params.tesselation);
    for (auto& p : shp->pos) {
        p *= params.size / 2;
        swap(p.y, p.z);
    }
    for (auto& n : shp->norm) n = {0, 1, 0};

    mat->name = pshp->name;
    mat->mtype = material_type::specular_roughness;
    mat->kd = {0.2f, 0.2f, 0.2f};
    mat->ks = zero3f;
    mat->rs = 1;

    return num_changed_shp || num_changed_mat;
}

// Updates for prim proc shape
bool update_proc_prim_shape(proc_shape* pshp) {
    auto num_changed_shp = add_proc_objects(pshp->shps, 1);
    auto num_changed_mat = add_proc_objects(pshp->mats, 1);

    auto& params = pshp->prim_params;
    auto shp = pshp->shps[0];
    auto mat = pshp->mats[0];

    shp->name = pshp->name;
    switch (params.ptype) {
        case proc_prim_shape_type::sphere: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere(params.tesselation);
        } break;
        default: throw runtime_error("should not have gotten here");
    }
    for (auto& p : shp->pos) p *= params.size / 2;

    mat->name = pshp->name;
    mat->mtype = material_type::specular_roughness;
    mat->kd = {0.2f, 0.2f, 0.2f};
    mat->ks = zero3f;
    mat->rs = 1;

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

inline void draw_proctree_widgets(
    gl_window* win, const string& lbl, proc_camera* cam, void*& selection) {
    draw_tree_widget_leaf(win, lbl + cam->name, selection, cam);
}

inline void draw_proctree_widgets(
    gl_window* win, const string& lbl, proc_shape* shp, void*& selection) {
    draw_tree_widget_leaf(win, lbl + shp->name, selection, shp);
}

inline void draw_proctree_widgets(
    gl_window* win, const string& lbl, proc_instance* ist, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + ist->name, selection, ist)) {
        if (ist->shp)
            draw_proctree_widgets(win, "shape: ", ist->shp, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_proctree_widgets(
    gl_window* win, const string& lbl, proc_scene* scn, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + "cameras")) {
        for (auto cam : scn->cameras)
            draw_proctree_widgets(win, "", cam, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "shapes")) {
        for (auto msh : scn->shapes)
            draw_proctree_widgets(win, "", msh, selection);
        draw_tree_widget_end(win);
    }

    if (draw_tree_widget_begin(win, lbl + "instances")) {
        for (auto ist : scn->instances)
            draw_proctree_widgets(win, "", ist, selection);
        draw_tree_widget_end(win);
    }
}

inline bool draw_procelem_widgets(gl_window* win, proc_scene* scn,
    proc_shape* shp, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", shp->name);
    edited += draw_value_widget(win, "type", shp->type, proc_shape_names());
    switch (shp->type) {
        case proc_shape_type::none: break;
        case proc_shape_type::floor: {
            auto& params = shp->floor_params;
            edited += draw_value_widget(win, "size", params.size, 0, 20);
            edited +=
                draw_value_widget(win, "tesselation", params.tesselation, 0, 8);
            edited += draw_color_widget(win, "color", params.color);
            edited += draw_value_widget(win, "textured", params.textured);
        } break;
        case proc_shape_type::prim: {
            auto& params = shp->prim_params;
            edited += draw_value_widget(win, "size", params.size, 0, 20);
            edited +=
                draw_value_widget(win, "tesselation", params.tesselation, 0, 8);
            edited +=
                draw_value_widget(win, "subdivision", params.subdivision, 0, 8);
            edited += draw_color_widget(win, "color", params.color);
            edited += draw_value_widget(win, "textured", params.textured);
        } break;
        default: throw runtime_error("should not have gotten here");
    }
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_procelem_widgets(gl_window* win, proc_scene* scn,
    proc_camera* cam, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", cam->name);
    edited += draw_value_widget(win, "from", cam->yfov, -10, 10);
    edited += draw_value_widget(win, "to", cam->yfov, -10, 10);
    edited += draw_value_widget(win, "up", cam->yfov, -10, 10);
    edited += draw_value_widget(win, "yfov", cam->yfov, 0.1, 4);
    edited += draw_value_widget(win, "aspect", cam->aspect, 0.1, 4);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_procelem_widgets(gl_window* win, proc_scene* scn,
    proc_instance* ist, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto shp_names = vector<pair<string, proc_shape*>>{{"<none>", nullptr}};
    for (auto shp : scn->shapes) shp_names.push_back({shp->name, shp});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", ist->name);
    edited += draw_value_widget(win, "frame", ist->frame, -10, 10);
    edited += draw_value_widget(win, "shape", ist->shp, shp_names);
    edited += draw_value_widget(win, "type", ist->type, proc_instance_names());
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_procelem_widgets(gl_window* win, proc_scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    for (auto cam : scn->cameras) {
        if (cam == selection)
            return draw_procelem_widgets(win, scn, cam, selection, gl_txt);
    }

    for (auto shp : scn->shapes) {
        if (shp == selection)
            return draw_procelem_widgets(win, scn, shp, selection, gl_txt);
    }

    for (auto ist : scn->instances) {
        if (ist == selection)
            return draw_procelem_widgets(win, scn, ist, selection, gl_txt);
    }

    return false;
}

inline bool draw_edit_widgets(gl_window* win, proc_scene* scn) {
    static auto stype = proc_shape_type::prim;
    static auto cur_shape = 0;
    draw_value_widget(win, "add type", stype, proc_shape_names());
    if (draw_button_widget(win, "add")) {
        auto shp = new proc_shape();
        shp->type = stype;
        shp->name = "shape" + to_string(cur_shape++);
        scn->shapes += shp;
        auto ist = new proc_instance();
        ist->name = shp->name;
        ist->shp = shp;
        scn->instances += ist;
        update_proc_scene(scn);
    }
    return false;
}

inline bool draw_procscene_widgets(gl_window* win, const string& lbl,
    proc_scene* scn, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    if (draw_header_widget(win, lbl)) {
        // draw_scroll_widget_begin(win, "model", 240, false);
        draw_proctree_widgets(win, "", scn, selection);
        draw_edit_widgets(win, scn);
        // draw_scroll_widget_end(win);
        return draw_procelem_widgets(win, scn, selection, gl_txt);
    } else
        return false;
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
    bool scene_updated = false;

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
        app->scene_updated = draw_procscene_widgets(
            win, "proc scene", app->pscn, app->selection, app->shstate->txt);
        draw_scene_widgets(
            win, "scene", app->scn, app->selection, app->shstate->txt);
    }
    end_widgets(win);

    swap_buffers(win);
}

// update the scene from a changed to the proc scene
void update(app_state* app) {
    if (!app->scene_updated) return;
    for (auto pcam : app->pscn->cameras) {
        if (pcam != app->selection) continue;
        update_proc_camera(pcam);
    }
    for (auto pshp : app->pscn->shapes) {
        if (pshp != app->selection) continue;
        update_proc_shape(pshp);
        update_stdsurface_state(app->shstate, app->scn, app->shparams,
            unordered_set<shape*>{pshp->shps.begin(), pshp->shps.end()}, {});
    }
    for (auto pist : app->pscn->instances) {
        if (pist != app->selection) continue;
        update_proc_instance(pist);
    }
    app->scene_updated = false;
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

        // update
        update(app);

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
