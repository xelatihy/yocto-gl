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

// Application state
struct app_state {
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
        if (scn) delete scn;
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

    app->shparams.highlighted = app->selection;

    gl_clear_buffers(app->shparams.background);
    gl_enable_depth_test(true);
    gl_enable_culling(app->shparams.cull_backface);
    draw_stdsurface_scene(app->shstate, app->scn, app->shparams);

    if (begin_widgets(win, "yview")) {
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
        "yview | " + app->filename, app);
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
    log_info("loading scene {}", app->filename);
    try {
        auto opts = load_options();
        opts.preserve_quads = preserve_quads;
        opts.preserve_facevarying = preserve_facevarying;
        app->scn = load_scene(app->filename, opts);
    } catch (exception e) { log_fatal("cannot load scene {}", app->filename); }

    // tesselate input shapes
    tesselate_shapes(app->scn, true, !preserve_facevarying,
        !preserve_quads && !preserve_facevarying, false);

    // add missing data
    add_elements(app->scn);

    // light
    update_lights(app->scn, false, false);

    // run ui
    auto cam = app->scn->cameras[app->shparams.camera_id];
    app->shparams.width = (int)round(cam->aspect * app->shparams.height);
    run_ui(app);

    // clear
    delete app;

    // done
    return 0;
}
