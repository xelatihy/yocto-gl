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
    ygl::gl_stdsurface_params params = {};
    ygl::gl_stdsurface_program prog;
    std::unordered_map<ygl::texture*, ygl::gl_texture> textures;
    std::unordered_map<ygl::shape*, ygl::gl_shape> shapes;
    ygl::gl_lights lights;
    bool navigation_fps = false;
    void* selection = nullptr;
    ygl::test_scene_params edit_params;
    float time = 0;
    ygl::vec2f time_range = ygl::zero2f;
    bool animate = false;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
    }
};

// draw with shading
inline void draw(ygl::gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    auto framebuffer_size = get_framebuffer_size(win);
    app->params.resolution = framebuffer_size.y;

    ygl::update_transforms(app->scn, app->time);
    app->lights = ygl::make_gl_lights(app->scn);
    ygl::update_textures(app->scn, app->textures);
    ygl::update_shapes(app->scn, app->shapes);
    if (app->lights.pos.empty()) app->params.eyelight = true;

    app->params.highlighted = app->selection;

    ygl::gl_clear_buffers(app->params.background);
    ygl::gl_enable_depth_test(true);
    ygl::gl_enable_culling(app->params.cull_backface);
    ygl::draw_stdsurface_scene(app->scn, app->cam, app->prog, app->shapes,
        app->textures, app->lights, framebuffer_size, app->params);

    if (ygl::begin_widgets(win, "yview")) {
        if (ygl::draw_header_widget(win, "view")) {
            ygl::draw_value_widget(win, "scene", app->filename);
            if (ygl::draw_button_widget(win, "new")) {
                app->edit_params = ygl::test_scene_presets().at("plane_al");
                delete app->scn;
                app->scn = new ygl::scene();
                ygl::update_test_scene(app->scn, app->edit_params);
                ygl::update_textures(app->scn, app->textures, {}, true);
                ygl::update_shapes(app->scn, app->shapes, {}, {}, true);
            }
            ygl::draw_continue_widget(win);
            if (ygl::draw_button_widget(win, "load")) {
                app->scn = ygl::load_scene(app->filename, {});
                ygl::update_textures(app->scn, app->textures, {}, true);
                ygl::update_shapes(app->scn, app->shapes, {}, {}, true);
            }
            ygl::draw_continue_widget(win);
            if (ygl::draw_button_widget(win, "save")) {
                ygl::save_scene(app->filename, app->scn, {});
            }
            ygl::draw_continue_widget(win);
            if (draw_button_widget(win, "save proc")) {
                ygl::save_test_scene(
                    ygl::replace_path_extension(app->filename, ".json"),
                    app->edit_params);
            }
            ygl::draw_camera_selection_widget(
                win, "camera", app->cam, app->scn, app->view);
            ygl::draw_value_widget(win, "fps", app->navigation_fps);
            if (app->time_range != ygl::zero2f) {
                ygl::draw_value_widget(win, "time", app->time,
                    app->time_range.x, app->time_range.y);
                ygl::draw_value_widget(win, "animate", app->animate);
            }
        }
        if (ygl::draw_header_widget(win, "params")) {
            ygl::draw_params_widgets(win, "", app->params);
        }
        if (ygl::draw_header_widget(win, "scene")) {
            if (ygl::draw_scene_widgets(win, "", app->scn, app->selection,
                    app->textures, &app->edit_params)) {
                ygl::update_textures(
                    app->scn, app->textures, {(ygl::texture*)app->selection});
                ygl::update_shapes(app->scn, app->shapes,
                    {(ygl::shape*)app->selection},
                    {(ygl::shape_group*)app->selection});
            }
        }
    }
    ygl::end_widgets(win);

    ygl::swap_buffers(win);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = ygl::make_window(
        (int)std::round(app->cam->aspect * app->params.resolution),
        app->params.resolution, "yview | " + app->filename, app);
    ygl::set_window_callbacks(win, nullptr, nullptr, draw);

    // load textures and vbos
    app->prog = ygl::make_stdsurface_program();
    ygl::update_textures(app->scn, app->textures);
    ygl::update_shapes(app->scn, app->shapes);

    // init widget
    ygl::init_widgets(win);

    // loop
    while (!should_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->cam == app->view) {
            ygl::handle_camera_navigation(win, app->view, app->navigation_fps);
        }

        // animation
        if (app->animate) {
            app->time += 1 / 60.0f;
            if (app->time < app->time_range.x || app->time > app->time_range.y)
                app->time = app->time_range.x;
        }

        // draw
        draw(win);

        // event hadling
        ygl::poll_events(win);
    }

    ygl::clear_window(win);
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser =
        ygl::make_parser(argc, argv, "yview", "views scenes inteactively");
    app->params = ygl::parse_params(parser, "", app->params);
    auto preserve_quads = ygl::parse_flag(
        parser, "--preserve-quads", "-q", "Preserve quads on load");
    auto preserve_facevarying = ygl::parse_flag(
        parser, "--preserve-facevarying", "-f", "Preserve facevarying on load");
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // scene loading
    ygl::log_info("loading scene {}", app->filename);
    try {
        auto opts = ygl::load_options();
        opts.preserve_quads = preserve_quads;
        opts.preserve_facevarying = preserve_facevarying;
        opts.preserve_hierarchy = true;
        app->scn = load_scene(app->filename, opts);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }

    // tesselate input shapes
    ygl::tesselate_shapes(app->scn, true, !preserve_facevarying,
        !preserve_quads && !preserve_facevarying, false);

    // add missing data
    ygl::add_elements(app->scn);

    // view camera
    app->view = ygl::make_view_camera(app->scn, 0);
    app->cam = app->view;

    // animation
    app->time_range = ygl::compute_animation_range(app->scn);
    app->time = app->time_range.x;

    // lights
    app->lights = ygl::make_gl_lights(app->scn);

    // run ui
    run_ui(app);

    // clear
    delete app;

    // done
    return 0;
}
