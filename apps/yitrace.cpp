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
    // scene data
    scene* scn = nullptr;

    // camera selection
    camera* scam = nullptr;

    // filenames
    string filename;
    string imfilename;

    // render
    int resolution = 0;
    float exposure = 0, gamma = 2.2f;
    bool filmic = false;
    vec4f background = {0, 0, 0, 0};

    // ui
    bool scene_updated = false;
    bool update_bvh = false;

    // navigation
    bool navigation_fps = false;

    // trace
    trace_params trace_params_;
    bool trace_save_progressive = false;
    int trace_block_size = 32;
    int trace_batch_size = 16;
    int trace_nthreads = 0;

    // interactive trace
    gl_texture trace_texture = {};
    gl_stdimage_program gl_prog = {};
    int trace_blocks_per_update = 8;
    image4f trace_img;
    image4f preview_img;
    bool trace_async_rendering = false;
    vector<rng_pcg32> trace_rngs;
    thread_pool* trace_pool = nullptr;
    int trace_cur_sample = 0;

    // editing support
    void* selection = nullptr;

    ~app_state() {
        if (trace_pool) delete trace_pool;
        if (scn) delete scn;
    }
};

void draw(gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    // begin frame
    gl_clear_buffers(app->background);

    // update texture
    update_texture(app->trace_texture, app->trace_img);

    // draw image
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    gl_set_viewport(framebuffer_size);
    draw_image(app->gl_prog, app->trace_texture, window_size, zero2f, 1,
        app->exposure, app->gamma, app->filmic);

    auto edited = vector<bool>();
    if (begin_widgets(win, "yitrace")) {
        draw_label_widget(win, "scene", app->filename);
        draw_label_widget(win, "sample", app->trace_cur_sample);
        draw_value_widget(
            win, "samples", app->trace_params_.nsamples, 1, 4096, 1);
        edited += draw_value_widget(
            win, "shader type", app->trace_params_.stype, trace_shader_names());
        edited += draw_value_widget(
            win, "random type", app->trace_params_.rtype, trace_rng_names());
        edited += draw_value_widget(
            win, "filter type", app->trace_params_.ftype, trace_filter_names());
        edited += draw_camera_widget(win, "camera", app->scn, app->scam);
        edited += draw_value_widget(win, "update bvh", app->update_bvh);
        draw_value_widget(win, "fps", app->navigation_fps);
        draw_tonemap_widgets(win, "", app->exposure, app->gamma, app->filmic);
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
        trace_async_stop(app->trace_pool);
        app->trace_cur_sample = 0;
        app->trace_async_rendering = false;

        // update BVH
        if (app->update_bvh) refit_bvh(app->scn);

        // render preview
        auto pparams = app->trace_params_;
        pparams.width = app->trace_params_.width / app->trace_block_size;
        pparams.height = app->trace_params_.height / app->trace_block_size;
        pparams.nsamples = 1;
        pparams.ftype = trace_filter_type::box;
        app->preview_img = image4f(pparams.width, pparams.height);
        auto prngs = trace_rngs(pparams);
        trace_samples(app->scn, app->preview_img, 0, 1, prngs, pparams);
        resize_image(app->preview_img, app->trace_img, resize_filter::box);
        update_texture(app->trace_texture, app->trace_img);

        app->scene_updated = false;
    } else if (!app->trace_async_rendering) {
        trace_async_start(app->scn, app->trace_img, app->trace_rngs,
            app->trace_params_, app->trace_pool,
            [app](int s) { app->trace_cur_sample = s; });
        app->trace_async_rendering = true;
    }
    return true;
}

// run ui loop
void run_ui(app_state* app, int w, int h) {
    // window
    auto win = make_window(w, h, "yitrace", app);
    set_window_callbacks(win, nullptr, nullptr, draw);

    // window values
    int mouse_button = 0;
    vec2f mouse_pos, mouse_last;

    // load textures
    app->gl_prog = make_stdimage_program();
    app->trace_texture = make_texture(app->trace_img, false, false, true);

    // init widget
    init_widgets(win);

    // loop
    while (!should_close(win)) {
        mouse_last = mouse_pos;
        mouse_pos = get_mouse_posf(win);
        mouse_button = get_mouse_button(win);

        set_window_title(win, ("yshade | " + app->filename));

        // handle mouse and keyboard for navigation
        if (mouse_button && !get_widget_active(win)) {
            if (app->navigation_fps) {
                auto dolly = 0.0f;
                auto pan = zero2f;
                auto rotate = zero2f;
                switch (mouse_button) {
                    case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                    case 2:
                        dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f;
                        break;
                    case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                    default: break;
                }
                camera_fps(app->scam->frame, {0, 0, 0}, rotate);
            } else {
                auto dolly = 0.0f;
                auto pan = zero2f;
                auto rotate = zero2f;
                switch (mouse_button) {
                    case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                    case 2:
                        dolly = (mouse_pos[0] - mouse_last[0]) / 100.0f;
                        break;
                    case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                    default: break;
                }

                camera_turntable(
                    app->scam->frame, app->scam->focus, rotate, dolly, pan);
            }
            app->scene_updated = true;
        }

        // handle keytboard for navigation
        if (!get_widget_active(win) && app->navigation_fps) {
            auto transl = zero3f;
            if (get_key(win, 'a')) transl.x -= 1;
            if (get_key(win, 'd')) transl.x += 1;
            if (get_key(win, 's')) transl.z += 1;
            if (get_key(win, 'w')) transl.z -= 1;
            if (get_key(win, 'e')) transl.y += 1;
            if (get_key(win, 'q')) transl.y -= 1;
            if (transl != zero3f) {
                camera_fps(app->scam->frame, transl, {0, 0});
                app->scene_updated = true;
            }
        }

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
    app->trace_params_.camera_id = 0;
    app->trace_save_progressive =
        parse_flag(parser, "--save-progressive", "", "save progressive images");
    app->trace_params_.rtype = parse_opt(parser, "--random", "", "random type",
        trace_rng_names(), trace_rng_type::stratified);
    app->trace_params_.ftype = parse_opt(parser, "--filter", "", "filter type",
        trace_filter_names(), trace_filter_type::box);
    app->trace_params_.stype =
        parse_opt(parser, "--shader", "-S", "path estimator type",
            trace_shader_names(), trace_shader_type::pathtrace);
    app->trace_params_.envmap_invisible =
        parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
    app->trace_params_.shadow_notransmission = parse_flag(
        parser, "--shadow-notransmission", "", "shadow without transmission");
    app->trace_block_size =
        parse_opt(parser, "--block-size", "", "block size", 32);
    app->trace_batch_size =
        parse_opt(parser, "--batch-size", "", "batch size", 16);
    app->trace_params_.nsamples =
        parse_opt(parser, "--samples", "-s", "image samples", 256);
    app->exposure =
        parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->gamma = parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->filmic = parse_flag(parser, "--filmic", "-F", "hdr image filmic");
    app->resolution =
        parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    auto camera_lights =
        parse_flag(parser, "--camera-lights", "-c", "enable camera lights");
    app->trace_params_.amb = {amb, amb, amb};
    if (camera_lights) {
        app->trace_params_.stype = trace_shader_type::eyelight;
    }
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
    app->scam = app->scn->cameras[0];

    // build bvh
    log_info("building bvh");
    build_bvh(app->scn);

    // init renderer
    log_info("initializing tracer");
    update_lights(app->scn, false);

    // initialize rendering objects
    auto width = (int)round(app->scam->aspect * app->resolution);
    auto height = app->resolution;
    app->trace_params_.width = width;
    app->trace_params_.height = height;
    app->trace_img = image4f(width, height);
    app->trace_rngs = trace_rngs(app->trace_params_);
    app->trace_pool = new thread_pool();
    app->preview_img = image4f();
    app->scene_updated = true;

    // run interactive
    run_ui(app, width, height);

    // cleanup
    trace_async_stop(app->trace_pool);
    delete app;

    // done
    return 0;
}
