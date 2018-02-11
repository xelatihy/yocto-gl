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
    ygl::bvh_tree* bvh = nullptr;
    std::string filename;
    std::string imfilename;
    ygl::image4f img;
    ygl::image<ygl::trace_pixel> pixels;
    ygl::trace_lights lights;
    ygl::trace_params params;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    bool scene_updated = false;
    bool update_bvh = false;
    bool navigation_fps = false;
    bool save_progressive = false;
    bool rendering = false;
    ygl::gl_stdimage_params imparams = {};
    ygl::gl_texture trace_texture = {};
    ygl::gl_stdimage_program gl_prog = {};
    void* selection = nullptr;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
        if (bvh) delete bvh;
    }
};

void draw(ygl::gl_window* win) {
    auto app = (app_state*)get_user_pointer(win);

    // update texture
    update_texture(app->trace_texture, app->img);

    // draw image
    auto window_size = get_window_size(win);
    auto framebuffer_size = get_framebuffer_size(win);
    ygl::gl_set_viewport(framebuffer_size);
    app->imparams.win_size = window_size;
    ygl::draw_image(app->gl_prog, app->trace_texture, app->imparams);

    auto edited = 0;
    if (ygl::begin_widgets(win, "yitrace")) {
        ygl::draw_label_widget(win, "scene", app->filename);
        ygl::draw_label_widget(
            win, "size", "{} x {}", app->params.width, app->params.height);
        ygl::draw_label_widget(win, "sample", app->pixels.at(0, 0).sample);
        if (ygl::draw_header_widget(win, "trace")) {
            ygl::draw_value_widget(
                win, "samples", app->params.nsamples, 1, 4096, 1);
            edited += ygl::draw_value_widget(win, "shader type",
                app->params.stype, ygl::enum_names<ygl::trace_shader_type>());
            edited += ygl::draw_value_widget(win, "random type",
                app->params.rtype, ygl::enum_names<ygl::trace_rng_type>());
            edited += ygl::draw_value_widget(win, "filter type",
                app->params.ftype, ygl::enum_names<ygl::trace_filter_type>());
            edited += ygl::draw_camera_widget(
                win, "camera", app->cam, app->scn, app->view);
            edited +=
                ygl::draw_value_widget(win, "update bvh", app->update_bvh);
            ygl::draw_value_widget(win, "fps", app->navigation_fps);
        }
        if (ygl::draw_header_widget(win, "image")) {
            ygl::draw_imageview_widgets(win, "", app->imparams);
            ygl::draw_imageinspect_widgets(
                win, "", app->img, {}, get_mouse_posf(win), app->imparams);
        }
        edited +=
            ygl::draw_scene_widgets(win, "scene", app->scn, app->selection, {});
    }
    ygl::end_widgets(win);
    app->scene_updated = app->scene_updated || (bool)edited;

    ygl::swap_buffers(win);
}

bool update(app_state* app) {
    if (app->scene_updated) {
        ygl::trace_async_stop(app->async_threads, app->async_stop);
        app->rendering = false;

        // update BVH
        if (app->update_bvh) ygl::refit_bvh(app->bvh, app->scn, false);

        // render preview
        auto pparams = app->params;
        pparams.width = app->params.width / app->params.block_size;
        pparams.height = app->params.height / app->params.block_size;
        pparams.nsamples = 1;
        pparams.ftype = ygl::trace_filter_type::box;
        auto preview_pixels = make_trace_pixels(pparams);
        auto preview_img = ygl::image4f(pparams.width, pparams.height);
        ygl::trace_samples(app->scn, app->cam, app->bvh, app->lights,
            preview_img, preview_pixels, 1, pparams);
        ygl::resize_image(preview_img, app->img, ygl::resize_filter::box);
        ygl::update_texture(app->trace_texture, app->img);

        app->scene_updated = false;
    } else if (!app->rendering) {
        ygl::trace_async_start(app->scn, app->cam, app->bvh, app->lights,
            app->img, app->pixels, app->async_threads, app->async_stop,
            app->params);
        app->rendering = true;
    }
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = ygl::make_window(app->params.width, app->params.height,
        "yitrace | " + app->filename, app);
    ygl::set_window_callbacks(win, nullptr, nullptr, draw);

    // load textures
    app->gl_prog = ygl::make_stdimage_program();
    app->trace_texture = ygl::make_texture(app->img, false, false, true);

    // init widget
    ygl::init_widgets(win);

    // loop
    while (!ygl::should_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->cam == app->view) {
            if (ygl::handle_camera_navigation(
                    win, app->view, app->navigation_fps))
                app->scene_updated = true;
        }

        // draw
        draw(win);

        // update
        update(app);

        // event hadling
        ygl::poll_events(win);
    }

    ygl::clear_window(win);
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = ygl::make_parser(
        argc, argv, "yitrace", "path trace images interactively");
    app->save_progressive = ygl::parse_flag(
        parser, "--save-progressive", "", "save progressive images");
    app->params.rtype = ygl::parse_opt(parser, "--random", "", "random type",
        ygl::enum_names<ygl::trace_rng_type>(),
        ygl::trace_rng_type::stratified);
    app->params.ftype = ygl::parse_opt(parser, "--filter", "", "filter type",
        ygl::enum_names<ygl::trace_filter_type>(), ygl::trace_filter_type::box);
    app->params.stype = ygl::parse_opt(parser, "--shader", "-S",
        "path estimator type", ygl::enum_names<ygl::trace_shader_type>(),
        ygl::trace_shader_type::pathtrace);
    app->params.envmap_invisible =
        ygl::parse_flag(parser, "--envmap-invisible", "", "envmap invisible");
    app->params.shadow_notransmission = ygl::parse_flag(
        parser, "--shadow-notransmission", "", "shadow without transmission");
    app->params.block_size =
        ygl::parse_opt(parser, "--block-size", "", "block size", 32);
    app->params.batch_size =
        ygl::parse_opt(parser, "--batch-size", "", "batch size", 16);
    app->params.nsamples =
        ygl::parse_opt(parser, "--samples", "-s", "image samples", 256);
    app->imparams.exposure =
        ygl::parse_opt(parser, "--exposure", "-e", "hdr image exposure", 0.0f);
    app->imparams.gamma =
        ygl::parse_opt(parser, "--gamma", "-g", "hdr image gamma", 2.2f);
    app->imparams.filmic =
        ygl::parse_flag(parser, "--filmic", "-F", "hdr image filmic");
    app->params.height =
        ygl::parse_opt(parser, "--resolution", "-r", "image resolution", 540);
    auto amb = ygl::parse_opt(parser, "--ambient", "", "ambient factor", 0.0f);
    auto camera_lights = ygl::parse_flag(
        parser, "--camera-lights", "-c", "enable camera lights");
    app->params.ambient = {amb, amb, amb};
    if (camera_lights) { app->params.stype = ygl::trace_shader_type::eyelight; }
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setting up rendering
    ygl::log_info("loading scene {}", app->filename);
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }

    // add elements
    auto opts = ygl::add_elements_options();
    ygl::add_elements(app->scn, opts);

    // view camera
    app->view = ygl::make_view_camera(app->scn, 0);
    app->cam = app->view;

    // build bvh
    ygl::log_info("building bvh");
    app->bvh = ygl::make_bvh(app->scn);

    // init renderer
    ygl::log_info("initializing tracer");
    app->lights = ygl::make_trace_lights(app->scn);

    // fix renderer type if no lights
    if (app->lights.empty() &&
        app->params.stype != ygl::trace_shader_type::eyelight) {
        ygl::log_info("no lights presents, switching to eyelight shader");
        app->params.stype = ygl::trace_shader_type::eyelight;
    }

    // initialize rendering objects
    app->params.width = (int)round(app->cam->aspect * app->params.height);
    app->img = ygl::image4f(app->params.width, app->params.height);
    app->pixels = ygl::make_trace_pixels(app->params);
    app->scene_updated = true;

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->async_threads, app->async_stop);
    delete app;

    // done
    return 0;
}
