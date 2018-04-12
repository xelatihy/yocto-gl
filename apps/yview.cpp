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

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* view = nullptr;
    ygl::camera* cam = nullptr;
    std::string filename;
    std::string imfilename;
    std::string outfilename;
    ygl::load_options loadopts;
    ygl::glsurface_params params = {};
    ygl::glsurface_program prog;
    std::unordered_map<ygl::texture*, ygl::gltexture> textures;
    std::unordered_map<ygl::shape*, ygl::glshape> shapes;
    ygl::gllights lights;
    bool navigation_fps = false;
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    float time = 0;
    std::string anim_group = "";
    ygl::vec2f time_range = ygl::zero2f;
    bool animate = false;
    bool quiet = false;
    bool screenshot_and_exit = false;
    bool no_widgets = false;
    std::unordered_map<std::string, std::string> inspector_highlights;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
    }
};

// draw with shading
inline void draw(ygl::glwindow* win, app_state* app) {
    auto framebuffer_size = get_glframebuffer_size(win);
    app->params.resolution = framebuffer_size.y;

    static auto last_time = 0.0f;
    for (auto& sel : app->update_list) {
        if (sel.is<ygl::texture>()) {
            ygl::update_gltexture(sel.get<ygl::texture>(),
                app->textures[sel.get<ygl::texture>()]);
        }
        if (sel.is<ygl::shape>()) {
            ygl::update_glshape(
                sel.get<ygl::shape>(), app->shapes[sel.get<ygl::shape>()]);
        }
        if (sel.is<ygl::node>() || sel.is<ygl::animation>() ||
            app->time != last_time) {
            ygl::update_transforms(app->scn, app->time, app->anim_group);
            last_time = app->time;
        }
        if (sel.is<ygl::shape>() || sel.is<ygl::material>() ||
            sel.is<ygl::node>()) {
            app->lights = ygl::make_gllights(app->scn);
            if (app->lights.pos.empty()) app->params.eyelight = true;
        }
    }
    app->update_list.clear();

    ygl::clear_glbuffers(app->params.background);
    ygl::enable_gldepth_test(true);
    ygl::enable_glculling(app->params.cull_backface);
    ygl::draw_glsurface_scene(app->scn, app->cam, app->prog, app->shapes,
        app->textures, app->lights, framebuffer_size,
        app->selection.get_untyped(), app->params);

    if (app->no_widgets) {
        ygl::swap_glwindow_buffers(win);
        return;
    }

    if (ygl::begin_imgui_frame(win, "yview")) {
        if (ygl::draw_imgui_header(win, "view")) {
            ygl::push_imgui_groupid(win, app);
            ygl::draw_imgui_text(win, "scene", app->filename);
            if (ygl::draw_imgui_button(win, "new")) {
                delete app->scn;
                ygl::clear_gltextures(app->textures);
                ygl::clear_glshapes(app->shapes);
                app->scn = new ygl::scene();
                app->textures = ygl::make_gltextures(app->scn);
                app->shapes = ygl::make_glshapes(app->scn);
            }
            ygl::continue_imgui_line(win);
            if (ygl::draw_imgui_button(win, "load")) {
                delete app->scn;
                app->scn = ygl::load_scene(app->filename, {});
                ygl::clear_gltextures(app->textures);
                ygl::clear_glshapes(app->shapes);
                app->textures = ygl::make_gltextures(app->scn);
                app->shapes = ygl::make_glshapes(app->scn);
            }
            ygl::continue_imgui_line(win);
            if (ygl::draw_imgui_button(win, "save")) {
                ygl::save_scene(app->filename, app->scn, {});
            }
            ygl::draw_imgui_camera_selector(
                win, "camera", app->cam, app->scn, app->view);
            ygl::draw_imgui_camera_inspector(win, "camera", app->cam);
            ygl::draw_imgui_checkbox(win, "fps", app->navigation_fps);
            if (app->time_range != ygl::zero2f) {
                ygl::draw_imgui_dragbox(win, "time", app->time,
                    app->time_range.x, app->time_range.y);
                ygl::draw_imgui_text(win, "anim group", app->anim_group);
                ygl::draw_imgui_checkbox(win, "animate", app->animate);
            }
            if (ygl::draw_imgui_button(win, "print stats"))
                ygl::print_stats(app->scn);
            ygl::pop_imgui_groupid(win);
        }
        if (ygl::draw_imgui_header(win, "params")) {
            ygl::draw_imgui_stdsurface_inspector(win, "", app->params);
        }
        if (ygl::draw_imgui_header(win, "scene")) {
            ygl::draw_imgui_scene_tree(win, "", app->scn, app->selection,
                app->update_list, app->inspector_highlights);
            ygl::draw_imgui_scene_inspector(win, "", app->scn, app->selection,
                app->update_list, app->inspector_highlights);
        }
    }
    ygl::end_imgui_frame(win);

    ygl::swap_glwindow_buffers(win);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = ygl::make_glwindow(
        (int)std::round(app->cam->aspect * app->params.resolution),
        app->params.resolution, "yview | " + app->filename);
    ygl::set_glwindow_callbacks(
        win, nullptr, nullptr, [win, app]() { draw(win, app); });

    // load textures and vbos
    app->prog = ygl::make_glsurface_program();
    app->textures = ygl::make_gltextures(app->scn);
    app->shapes = ygl::make_glshapes(app->scn);
    ygl::update_transforms(app->scn, app->time);
    app->lights = ygl::make_gllights(app->scn);
    if (app->lights.pos.empty()) app->params.eyelight = true;

    // init widget
    ygl::init_imgui(win);

    // loop
    while (!should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->cam == app->view) {
            ygl::handle_glcamera_navigation(
                win, app->view, app->navigation_fps);
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
            ygl::log_info("taking screenshot and exiting...");
            ygl::save_image4b(app->imfilename, ygl::take_glscreenshot4b(win));
            break;
        }

        // event hadling
        if (ygl::get_glmouse_button(win) || ygl::get_imgui_active(win) ||
            app->animate) {
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
    app->quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    app->screenshot_and_exit = ygl::parse_flag(
        parser, "--screenshot-and-exit", "", "Take a screenshot and exit");
    app->no_widgets =
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
    if (app->quiet) ygl::get_default_logger()->verbose = false;

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
        app->scn = load_scene(app->filename, app->loadopts);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }

    // tesselate input shapes
    ygl::tesselate_shapes(app->scn, true,
        !app->loadopts.obj_preserve_facevarying,
        !app->loadopts.obj_preserve_quads &&
            !app->loadopts.obj_preserve_facevarying,
        false);

    // add missing data
    ygl::add_names(app->scn);
    ygl::add_tangent_space(app->scn);

    // validate
    ygl::validate(app->scn, false, true);

    // view camera
    app->view = ygl::make_view_camera(app->scn, 0);
    app->cam = app->view;

    // animation
    app->time_range = ygl::compute_animation_range(app->scn);
    app->time = app->time_range.x;

    // lights
    app->lights = ygl::make_gllights(app->scn);

    // run ui
    run_ui(app);

    // cleanup
    delete app;

    // done
    return 0;
}
