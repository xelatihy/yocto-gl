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

#include "../yocto/yocto_gl.h"
#include "../yocto/yocto_glio.h"
#include "../yocto/yocto_glutils.h"
#include "CLI11.hpp"
#include "yapp_ui.h"
using namespace std::literals;

// Application state
struct app_state {
    std::shared_ptr<ygl::scene> scn = nullptr;
    std::shared_ptr<ygl::camera> cam = nullptr;
    std::string filename = "scene.json";
    std::string imfilename = "out.png";
    std::string outfilename = "scene.json";

    int resolution = 512;                           // image resolution
    bool wireframe = false;                         // wireframe drawing
    bool edges = false;                             // draw edges
    float edge_offset = 0.01f;                      // offset for edges
    bool eyelight = false;                          // camera light mode
    float exposure = 0;                             // exposure
    float gamma = 2.2f;                             // gamma
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};  // background
    ygl::vec3f ambient = {0, 0, 0};                 // ambient lighting

    ygl::glsurface_program prog;

    bool navigation_fps = false;
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    float time = 0;
    std::string anim_group = "";
    ygl::vec2f time_range = ygl::zero2f;
    bool animate = false;
    bool quiet = false;
    bool screenshot_and_exit = false;
    bool no_glwidgets = false;
    std::unordered_map<std::string, std::string> inspector_highlights;
    bool double_sided = false;
    std::string highlight_filename = ""s;
};

// draw with shading
inline void draw(const std::shared_ptr<ygl::glwindow>& win,
    const std::shared_ptr<app_state>& app) {
    auto framebuffer_size = get_glwindow_framebuffer_size(win);
    app->resolution = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app->update_list) {
        if (sel.as<ygl::texture>()) {
            ygl::update_gldata(sel.as<ygl::texture>());
        }
        if (sel.as<ygl::subdiv>()) {
            auto sbd = sel.as<ygl::subdiv>();
            for (auto ist : app->scn->instances) {
                if (ist->sbd != sbd) continue;
                ygl::update_tesselation(ist->sbd, ist->shp);
                ygl::update_gldata(ist->shp);
                break;
            }
        }
        if (sel.as<ygl::shape>()) { ygl::update_gldata(sel.as<ygl::shape>()); }
        if (sel.as<ygl::node>() || sel.as<ygl::animation>() ||
            app->time != last_time) {
            ygl::update_transforms(app->scn, app->time, app->anim_group);
            last_time = app->time;
        }
        if (sel.as<ygl::shape>() || sel.as<ygl::material>() ||
            sel.as<ygl::node>()) {
            ygl::update_lights(app->scn);
            if (app->scn->lights.empty()) app->eyelight = true;
        }
    }
    app->update_list.clear();

    ygl::clear_glbuffers(app->background);
    ygl::enable_gldepth_test(true);
    ygl::draw_glscene(app->scn, app->cam, app->prog, framebuffer_size,
        app->selection.ptr, app->eyelight, app->wireframe, app->edges,
        app->exposure, app->gamma);

    if (app->no_glwidgets) {
        ygl::swap_glwindow_buffers(win);
        return;
    }

    if (ygl::begin_glwidgets_frame(win, "yview")) {
        ygl::draw_glwidgets_text(win, "scene", app->filename);
        if (ygl::draw_glwidgets_button(win, "new")) {
            ygl::clear_gldata(app->scn);
            app->scn = std::make_shared<ygl::scene>();
            ygl::update_gldata(app->scn);
        }
        ygl::continue_glwidgets_line(win);
        if (ygl::draw_glwidgets_button(win, "load")) {
            ygl::clear_gldata(app->scn);
            app->scn = ygl::load_scene(app->filename);
            ygl::update_gldata(app->scn);
        }
        ygl::continue_glwidgets_line(win);
        if (ygl::draw_glwidgets_button(win, "save")) {
            ygl::save_scene(app->filename, app->scn);
        }
        ygl::continue_glwidgets_line(win);
        if (ygl::draw_glwidgets_button(win, "print stats"))
            ygl::print_stats(app->scn);
        if (app->time_range != ygl::zero2f) {
            ygl::draw_glwidgets_dragbox(
                win, "time", app->time, app->time_range.x, app->time_range.y);
            ygl::draw_glwidgets_text(win, "anim group", app->anim_group);
            ygl::draw_glwidgets_checkbox(win, "animate", app->animate);
        }
        if (ygl::begin_glwidgets_tree(win, "render settings")) {
            ygl::draw_glwidgets_combobox(
                win, "camera", app->cam, app->scn->cameras);
            ygl::draw_glwidgets_dragbox(
                win, "resolution", app->resolution, 256, 4096);
            ygl::draw_glwidgets_checkbox(win, "eyelight", app->eyelight);
            ygl::continue_glwidgets_line(win);
            ygl::draw_glwidgets_checkbox(win, "wireframe", app->wireframe);
            ygl::continue_glwidgets_line(win);
            ygl::draw_glwidgets_checkbox(win, "edges", app->edges);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "view settings")) {
            ygl::draw_glwidgets_dragbox(
                win, "exposure", app->exposure, -10, 10);
            ygl::draw_glwidgets_dragbox(win, "gamma", app->gamma, 0.1f, 4);
            ygl::draw_glwidgets_colorbox(win, "background", app->background);
            ygl::draw_glwidgets_checkbox(win, "fps", app->navigation_fps);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "scene tree")) {
            ygl::draw_glwidgets_scene_tree(win, "", app->scn, app->selection,
                app->update_list, 200, app->inspector_highlights);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "scene object")) {
            ygl::draw_glwidgets_scene_inspector(win, "", app->scn,
                app->selection, app->update_list, 200,
                app->inspector_highlights);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "log")) {
            ygl::draw_glwidgets_log(win, 200);
            ygl::end_glwidgets_tree(win);
        }
    }
    ygl::end_glwidgets_frame(win);

    ygl::swap_glwindow_buffers(win);
}

// run ui loop
void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win = ygl::make_glwindow(
        (int)std::round(app->resolution * app->cam->width / app->cam->height) +
            ygl::default_glwidgets_width,
        app->resolution, "yview | " + app->filename);
    ygl::set_glwindow_callbacks(
        win, nullptr, nullptr, [win, app]() { draw(win, app); });

    // load textures and vbos
    app->prog = ygl::make_glsurface_program();
    ygl::update_gldata(app->scn);
    ygl::update_transforms(app->scn, app->time);
    ygl::update_lights(app->scn);
    if (app->scn->lights.empty()) app->eyelight = true;

    // init widget
    ygl::init_glwidgets(win);

    // loop
    while (!should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->navigation_fps) {
            ygl::handle_glcamera_fps(win, app->cam);
        } else {
            ygl::handle_glcamera_turntable(win, app->cam,
                ygl::handle_glcamera_turntable_dist(
                    win, app->cam, app->scn, app->selection));
        }

        // animation
        if (app->animate) {
            app->time += 1 / 60.0f;
            if (app->time < app->time_range.x || app->time > app->time_range.y)
                app->time = app->time_range.x;
            ygl::update_transforms(app->scn, app->time);
        }

        // draw
        draw(win, app);

        // check if exiting is needed
        if (app->screenshot_and_exit) {
            if (!app->quiet) std::cout << "taking screenshot and exiting\n";
            auto img = ygl::take_glwindow_screenshot(win);
            ygl::save_image(app->imfilename, img);
            break;
        }

        // event hadling
        if (ygl::get_glwindow_mouse_button(win) ||
            ygl::get_glwidgets_active(win) || app->animate) {
            ygl::poll_glwindow_events(win);
        } else {
            ygl::wait_glwindow_events(win);
        }
    }
}

// Load INI file. The implementation does not handle escaping.
std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
load_ini(const std::string& filename) {
    auto f = fopen(filename.c_str(), "rt");
    if (!f) throw std::runtime_error("cannot open " + filename);
    auto ret = std::unordered_map<std::string,
        std::unordered_map<std::string, std::string>>();
    auto cur_group = ""s;
    ret[""] = {};

    char buf[4096];
    while (fgets(buf, 4096, f)) {
        auto line = std::string(buf);
        if (line.empty()) continue;
        if (line.front() == ';') continue;
        if (line.front() == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') throw std::runtime_error("bad INI format");
            cur_group = line.substr(1, line.length() - 2);
            ret[cur_group] = {};
        } else if (line.find('=') != line.npos) {
            auto var = line.substr(0, line.find('='));
            auto val = line.substr(line.find('=') + 1);
            ret[cur_group][var] = val;
        } else {
            throw std::runtime_error("bad INI format");
        }
    }

    fclose(f);

    return ret;
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = std::make_shared<app_state>();

    // parse command line
    CLI::App parser("views scenes inteactively", "yview");
    parser.add_flag("--eyelight,-c", app->eyelight, "Eyelight rendering.");
    parser.add_flag(
        "--double-sided,-D", app->double_sided, "Double-sided rendering.");
    parser.add_flag("--quiet,-q", app->quiet, "Print only errors messages");
    parser.add_flag("--screenshot-and-exit", app->screenshot_and_exit,
        "Take a screenshot and exit");
    parser.add_flag("--no-widgets", app->no_glwidgets, "Disable widgets");
    parser.add_option(
        "--highlights", app->highlight_filename, "Highlight filename");
    parser.add_option("--output-image,-o", app->imfilename, "Image filename");
    parser.add_option("scene", app->filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }

    // fix hilights
    if (!app->highlight_filename.empty()) {
        try {
            app->inspector_highlights =
                load_ini(app->highlight_filename).at("");
        } catch (const std::exception& e) {
            std::cout << "cannot load highlight " << app->highlight_filename
                      << "\n";
            std::cout << "error: " << e.what() << "\n";
            exit(1);
        }
    }

    // scene loading
    if (!app->quiet) std::cout << "loading scene" << app->filename << "\n";
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << app->filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }

    // tesselate
    if (!app->quiet) std::cout << "tesselating scene elements\n";
    ygl::update_tesselation(app->scn);

    // update bbox and transforms
    ygl::update_transforms(app->scn);
    ygl::update_bbox(app->scn);

    // add components
    if (!app->quiet) std::cout << "adding scene elements\n";
    if (app->double_sided) {
        for (auto mat : app->scn->materials) mat->double_sided = true;
    }
    if (app->scn->cameras.empty()) {
        app->scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", app->scn->bbox));
    }
    app->cam = app->scn->cameras[0];
    ygl::add_missing_names(app->scn);
    for (auto err : ygl::validate(app->scn))
        std::cout << "warning: " << err << "\n";

    // animation
    app->time_range = ygl::compute_animation_range(app->scn);
    app->time = app->time_range.x;

    // lights
    ygl::update_bbox(app->scn);
    ygl::update_lights(app->scn);

    // run ui
    run_ui(app);

    // done
    return 0;
}
