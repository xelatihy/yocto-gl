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
using namespace ygl;

// Application state
struct app_state {
    scene* scn = nullptr;
    camera* view = nullptr;
    camera* cam = nullptr;
    string filename;
    string imfilename;
    string outfilename;
    gl_stdsurface_params params = {};
    gl_stdsurface_program prog;
    unordered_map<texture*, gl_texture> textures;
    unordered_map<shape*, gl_shape> shapes;
    gl_lights lights;
    bool navigation_fps = false;
    void* selection = nullptr;
    test_scene_params edit_params;
    float time = 0;
    vec2f time_range = zero2f;
    bool animate = false;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
    }
};

// draw with shading
inline void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    auto framebuffer_size = get_framebuffer_size(win);
    app->params.width = framebuffer_size.x;
    app->params.height = framebuffer_size.y;

    update_transforms(app->scn, app->time);
    app->lights = make_gl_lights(app->scn);
    update_textures(app->scn, app->textures);
    update_shapes(app->scn, app->shapes);
    if (app->lights.pos.empty()) app->params.camera_lights = true;

    app->params.highlighted = app->selection;

    gl_clear_buffers(app->params.background);
    gl_enable_depth_test(true);
    gl_enable_culling(app->params.cull_backface);
    draw_stdsurface_scene(app->scn, app->cam, app->prog, app->shapes,
        app->textures, app->lights, app->params);

    if (begin_widgets(win, "yview")) {
        if (draw_header_widget(win, "file")) {
            draw_value_widget(win, "scene", app->filename);
            if (draw_button_widget(win, "new")) {
                app->edit_params = test_scene_presets().at("plane_al");
                delete app->scn;
                app->scn = new scene();
                update_test_scene(app->scn, app->edit_params);
                update_textures(app->scn, app->textures, {}, true);
                update_shapes(app->scn, app->shapes, {}, {}, true);
            }
            draw_continue_widget(win);
            if (draw_button_widget(win, "load")) {
                app->scn = load_scene(app->filename, {});
                update_textures(app->scn, app->textures, {}, true);
                update_shapes(app->scn, app->shapes, {}, {}, true);
            }
            draw_continue_widget(win);
            if (draw_button_widget(win, "save")) {
                save_scene(app->filename, app->scn, {});
            }
            draw_continue_widget(win);
            if (draw_button_widget(win, "save proc")) {
                save_test_scene(replace_path_extension(app->filename, ".json"),
                    app->edit_params);
            }
        }
        if (draw_header_widget(win, "view")) {
            draw_camera_widget(win, "camera", app->cam, app->scn, app->view);
            draw_value_widget(win, "wire", app->params.wireframe);
            draw_continue_widget(win);
            draw_value_widget(win, "edges", app->params.edges);
            draw_continue_widget(win);
            draw_value_widget(win, "cutout", app->params.cutout);
            draw_continue_widget(win);
            draw_value_widget(win, "fps", app->navigation_fps);
            draw_tonemap_widgets(win, "", app->params.exposure,
                app->params.gamma, app->params.filmic);
            if (app->time_range != zero2f) {
                draw_value_widget(win, "time", app->time, app->time_range.x,
                    app->time_range.y);
                draw_value_widget(win, "animate", app->animate);
            }
        }
        if (draw_scene_widgets(win, "scene", app->scn, app->selection,
                app->textures, &app->edit_params)) {
            update_textures(
                app->scn, app->textures, {(texture*)app->selection});
            update_shapes(app->scn, app->shapes, {(shape*)app->selection},
                {(shape_group*)app->selection});
        }
    }
    end_widgets(win);

    swap_buffers(win);
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = make_window(
        app->params.width, app->params.height, "yview | " + app->filename, app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // load textures and vbos
    app->prog = make_stdsurface_program();
    update_textures(app->scn, app->textures);
    update_shapes(app->scn, app->shapes);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->cam == app->view) {
            handle_camera_navigation(win, app->view, app->navigation_fps);
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
        poll_events(win);
    }

    clear_window(win);
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = make_parser(argc, argv, "yview", "views scenes inteactively");
    app->params.exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->params.gamma =
        parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->params.filmic =
        parse_flag(parser, "--filmic", "-F", "hdr filmic output");
    app->params.height =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    app->params.ambient = {amb, amb, amb};
    app->params.camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
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
    log_info("loading scene {}", app->filename);
    try {
        auto opts = load_options();
        opts.preserve_quads = preserve_quads;
        opts.preserve_facevarying = preserve_facevarying;
        opts.preserve_hierarchy = true;
        app->scn = load_scene(app->filename, opts);
    } catch (exception e) { log_fatal("cannot load scene {}", app->filename); }

    // tesselate input shapes
    tesselate_shapes(app->scn, true, !preserve_facevarying,
        !preserve_quads && !preserve_facevarying, false);

    // add missing data
    add_elements(app->scn);

    // HACK
    print_info_visit(app->scn);

    // view camera
    app->view = make_view_camera(app->scn, 0);
    app->cam = app->view;

    // animation
    app->time_range = compute_animation_range(app->scn);
    app->time = app->time_range.x;

    // lights
    app->lights = make_gl_lights(app->scn);

    // run ui
    app->params.width = (int)round(app->cam->aspect * app->params.height);
    run_ui(app);

    // clear
    delete app;

    // done
    return 0;
}
