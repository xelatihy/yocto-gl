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

#include "../yocto/yocto_glutils.h"
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_utils.h"
#include "yapp_ui.h"
using namespace std::literals;

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* cam = nullptr;
    std::string filename;
    std::string imfilename;
    std::string outfilename;

    int resolution = 512;       // image resolution
    bool wireframe = false;     // wireframe drawing
    bool edges = false;         // draw edges
    float edge_offset = 0.01f;  // offset for edges
    bool cutout = false;        // draw with binary transparency
    bool eyelight = false;      // camera light mode
    float exposure = 0;         // exposure
    float gamma = 2.2f;         // gamma
    ygl::vec4b background = {222, 222, 222, 0};  // background
    ygl::vec3f ambient = {0, 0, 0};              // ambient lighting
    bool cull_backface = false;                  // culling back face

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

    ~app_state() {
        if (scn) delete scn;
    }
};

// draw with shading
inline void draw(ygl::glwindow* win, app_state* app) {
    auto framebuffer_size = get_glwindow_framebuffer_size(win);
    app->resolution = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app->update_list) {
        if (sel.as<ygl::texture>()) {
            ygl::update_gldata(sel.as<ygl::texture>());
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
    ygl::enable_glculling(app->cull_backface);
    ygl::draw_glscene(app->scn, app->cam, app->prog, framebuffer_size,
        app->selection.ptr, app->eyelight, app->wireframe, app->edges,
        app->cutout, app->exposure, app->gamma, app->cull_backface);

    if (app->no_glwidgets) {
        ygl::swap_glwindow_buffers(win);
        return;
    }

    if (ygl::begin_glwidgets_frame(win, "yview")) {
        ygl::draw_glwidgets_text(win, "scene", app->filename);
        if (ygl::draw_glwidgets_button(win, "new")) {
            delete app->scn;
            ygl::clear_gldata(app->scn);
            app->scn = new ygl::scene();
            ygl::update_gldata(app->scn);
        }
        ygl::continue_glwidgets_line(win);
        if (ygl::draw_glwidgets_button(win, "load")) {
            ygl::clear_gldata(app->scn);
            delete app->scn;
            app->scn = ygl::load_scene(app->filename, {});
            ygl::update_gldata(app->scn);
        }
        ygl::continue_glwidgets_line(win);
        if (ygl::draw_glwidgets_button(win, "save")) {
            ygl::save_scene(app->filename, app->scn, {});
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
            draw_glwidgets_dragbox(
                win, "resolution", app->resolution, 256, 4096);
            draw_glwidgets_checkbox(win, "eyelight", app->eyelight);
            draw_glwidgets_checkbox(win, "wireframe", app->wireframe);
            ygl::continue_glwidgets_line(win);
            draw_glwidgets_checkbox(win, "edges", app->edges);
            ygl::continue_glwidgets_line(win);
            draw_glwidgets_checkbox(win, "cutout", app->cutout);
            ygl::continue_glwidgets_line(win);
            draw_glwidgets_checkbox(win, "cull", app->cull_backface);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "view settings")) {
            draw_glwidgets_dragbox(win, "exposure", app->exposure, -10, 10);
            draw_glwidgets_dragbox(win, "gamma", app->gamma, 0.1f, 4);
            draw_glwidgets_colorbox(win, "background", app->background);
            ygl::draw_glwidgets_checkbox(win, "fps", app->navigation_fps);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "scene")) {
            ygl::draw_glwidgets_scene_tree(win, "", app->scn, app->selection,
                app->update_list, app->inspector_highlights);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "object")) {
            ygl::draw_glwidgets_scene_inspector(win, "", app->scn,
                app->selection, app->update_list, app->inspector_highlights);
            ygl::end_glwidgets_tree(win);
        }
    }
    ygl::end_glwidgets_frame(win);

    ygl::swap_glwindow_buffers(win);
}

inline void refresh(ygl::glwindow* win) {
    return draw(win, (app_state*)ygl::get_glwindow_user_pointer(win));
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win =
        ygl::make_glwindow((int)std::round(app->cam->aspect * app->resolution) +
                               ygl::default_glwidgets_width,
            app->resolution, "yview | " + app->filename, app);
    ygl::set_glwindow_callbacks(win, nullptr, nullptr, refresh);

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
        ygl::handle_glcamera_navigation(win, app->cam, app->navigation_fps);

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
            ygl::log_info("taking screenshot and exiting...");
            auto width = 0, height = 0;
            auto img = std::vector<ygl::vec4b>();
            ygl::take_glwindow_screenshot4b(win, width, height, img);
            ygl::save_image4b(app->imfilename, width, height, img);
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

    // cleanup
    delete win;
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
    auto app = new app_state();

    // parse command line
    auto parser =
        ygl::make_parser(argc, argv, "yview", "views scenes inteactively");
    app->eyelight = ygl::parse_flag(
        parser, "--eyelight", "-c", "Eyelight rendering.", false);
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided", "-D", "Force double sided rendering.", false);
    app->quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    app->screenshot_and_exit = ygl::parse_flag(
        parser, "--screenshot-and-exit", "", "Take a screenshot and exit");
    app->no_glwidgets =
        ygl::parse_flag(parser, "--no-widgets", "", "Disable widgets");
    auto highlight_filename =
        ygl::parse_opt(parser, "--highlights", "", "Highlight filename", ""s);
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.png"s);
    app->filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (app->quiet) ygl::log_verbose() = false;

    // fix hilights
    if (!highlight_filename.empty()) {
        try {
            app->inspector_highlights = load_ini(highlight_filename).at("");
        } catch (std::exception e) {
            ygl::log_fatal("cannot load highlihgt file {}", highlight_filename);
        }
    }

    // scene loading
    try {
        ygl::log_info("loading scene {}", app->filename);
        app->scn = ygl::load_scene(app->filename, true, true);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }

    // tesselate input shapes
    ygl::tesselate_shapes(app->scn, true, false, false, false);

    // fix scene
    ygl::update_bbox(app->scn);
    ygl::add_missing_camera(app->scn);
    ygl::add_missing_names(app->scn);
    ygl::add_missing_tangent_space(app->scn);
    app->cam = app->scn->cameras[0];
    if (double_sided) {
        for (auto mat : app->scn->materials) mat->double_sided = true;
    }

    // validate
    for (auto err : ygl::validate(app->scn)) ygl::log_warning(err);

    // animation
    app->time_range = ygl::compute_animation_range(app->scn);
    app->time = app->time_range.x;

    // lights
    ygl::update_bbox(app->scn);
    ygl::update_lights(app->scn);

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
