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
    // scene data
    scene* scn = nullptr;

    // filenames
    string filename;
    string imfilename;

    // ui
    bool scene_updated = false;
    bool update_bvh = false;

    // navigation
    bool navigation_fps = false;

    // trace
    trace_params params;
    bool save_progressive = false;
    trace_state* state = nullptr;
    bool rendering = false;

    // image view
    gl_stdimage_params imparams = {};
    gl_texture trace_texture = {};
    gl_stdimage_program gl_prog = {};

    // editing support
    void* selection = nullptr;

    ~app_state() {
        if (state) delete state;
        if (scn) delete scn;
    }
};

void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    // update texture
    update_texture(app->trace_texture, get_trace_image(app->state));

    // draw image
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    gl_set_viewport(framebuffer_size);
    app->imparams.win_size = window_size;
    draw_image(app->gl_prog, app->trace_texture, app->imparams);

    auto edited = vector<bool>();
    if (begin_widgets(win, "yitrace")) {
        draw_label_widget(win, "scene", app->filename);
        draw_label_widget(
            win, "size", "{} x {}", app->params.width, app->params.height);
        draw_label_widget(win, "sample", get_trace_sample(app->state));
        if (draw_header_widget(win, "trace")) {
            draw_value_widget(win, "samples", app->params.nsamples, 1, 4096, 1);
            edited += draw_value_widget(
                win, "shader type", app->params.stype, trace_shader_names());
            edited += draw_value_widget(
                win, "random type", app->params.rtype, trace_rng_names());
            edited += draw_value_widget(
                win, "filter type", app->params.ftype, trace_filter_names());
            edited += draw_camera_widget(
                win, "camera", app->scn, app->params.camera_id);
            edited += draw_value_widget(win, "update bvh", app->update_bvh);
            draw_value_widget(win, "fps", app->navigation_fps);
        }
        if (draw_header_widget(win, "image")) {
            draw_imageview_widgets(win, "", app->imparams);
            draw_imageinspect_widgets(win, "", get_trace_image(app->state), {},
                get_mouse_posf(win), app->imparams);
        }
        edited +=
            draw_scene_widgets(win, "scene", app->scn, app->selection, {});
    }
    end_widgets(win);
    app->scene_updated =
        app->scene_updated ||
        std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });

    swap_buffers(win);
}

bool update(app_state* app) {
    if (app->scene_updated) {
        trace_async_stop(app->state);
        app->rendering = false;

        // update BVH
        if (app->update_bvh) refit_bvh(app->scn);

        // render preview
        auto pparams = app->params;
        pparams.width = app->params.width / app->params.block_size;
        pparams.height = app->params.height / app->params.block_size;
        pparams.nsamples = 1;
        pparams.ftype = trace_filter_type::box;
        auto preview_state = make_trace_state(pparams);
        trace_samples(preview_state, app->scn, 1, pparams);
        resize_image(get_trace_image(preview_state),
            (image4f&)get_trace_image(app->state), resize_filter::box);
        update_texture(app->trace_texture, get_trace_image(app->state));
        delete preview_state;

        app->scene_updated = false;
    } else if (!app->rendering) {
        trace_async_start(app->state, app->scn, app->params);
        app->rendering = true;
    }
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = make_window(app->params.width, app->params.height,
        "yitrace | " + app->filename, app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // load textures
    app->gl_prog = make_stdimage_program();
    app->trace_texture =
        make_texture(get_trace_image(app->state), false, false, true);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        // handle mouse and keyboard for navigation
        auto cam = app->scn->cameras[app->params.camera_id];
        if (handle_camera_navigation(win, cam, app->navigation_fps))
            app->scene_updated = true;

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
    auto parser =
        make_parser(argc, argv, "yitrace", "path trace images interactively");
    app->params.camera_id = 0;
    app->save_progressive =
        parse_flag(parser, "--save-progressive", "", "save progressive images");
    app->params.rtype = parse_opt(parser, "--random", "", "random type",
        trace_rng_names(), trace_rng_type::stratified);
    app->params.ftype = parse_opt(parser, "--filter", "", "filter type",
        trace_filter_names(), trace_filter_type::box);
    app->params.stype =
        parse_opt(parser, "--shader", "-S", "path estimator type",
            trace_shader_names(), trace_shader_type::pathtrace);
    app->params.envmap_invisible =
        parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
    app->params.shadow_notransmission = parse_flag(
        parser, "--shadow-notransmission", "", "shadow without transmission");
    app->params.block_size =
        parse_opt(parser, "--block-size", "", "block size", 32);
    app->params.batch_size =
        parse_opt(parser, "--batch-size", "", "batch size", 16);
    app->params.nsamples =
        parse_opt(parser, "--samples", "-s", "image samples", 256);
    app->imparams.exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->imparams.gamma =
        parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->imparams.filmic =
        parse_flag(parser, "--filmic", "-F", "hdr image filmic");
    app->params.height =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    auto camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
    app->params.ambient = {amb, amb, amb};
    if (camera_lights) { app->params.stype = trace_shader_type::eyelight; }
    auto log_filename = parse_opt(parser, "--log", "", "log to disk", ""s);
    if (log_filename != "") add_file_stream(log_filename, true);
    app->imfilename =
        parse_opt(parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = parse_arg(parser, "scene", "scene filename", ""s);
    if (should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setting up rendering
    log_info("loading scene {}", app->filename);
    try {
        app->scn = load_scene(app->filename);
    } catch (exception e) { log_fatal("cannot load scene {}", app->filename); }
    add_elements(app->scn);

    // build bvh
    log_info("building bvh");
    build_bvh(app->scn);

    // init renderer
    log_info("initializing tracer");
    update_lights(app->scn, true, true);

    // fix renderer type if no lights
    if (app->scn->lights.empty() &&
        app->params.stype != trace_shader_type::eyelight) {
        log_info("no lights presents, switching to eyelight shader");
        app->params.stype = trace_shader_type::eyelight;
    }

    // initialize rendering objects
    auto cam = app->scn->cameras[app->params.camera_id];
    app->params.width = (int)round(cam->aspect * app->params.height);
    app->state = make_trace_state(app->params);
    app->scene_updated = true;

    // run interactive
    run_ui(app);

    // cleanup
    trace_async_stop(app->state);
    delete app;

    // done
    return 0;
}
