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
    int resolution = 512;
    ygl::image4f img;
    ygl::image<ygl::trace_pixel> pixels;
    ygl::trace_lights lights;
    ygl::trace_params params;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    bool scene_updated = false;
    bool navigation_fps = false;
    int preview_res = 64;
    bool rendering = false;
    ygl::gl_stdimage_params imparams = {};
    ygl::gl_texture trace_texture = {};
    ygl::gl_stdimage_program gl_prog = {};
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;

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
    ygl::draw_image(
        app->gl_prog, app->trace_texture, window_size, app->imparams);

    auto edited = 0;
    if (ygl::begin_widgets(win, "yitrace")) {
        if (ygl::draw_header_widget(win, "trace")) {
            ygl::draw_label_widget(win, "scene", app->filename);
            ygl::draw_label_widget(
                win, "size", "{} x {}", app->img.width(), app->img.height());
            ygl::draw_label_widget(win, "sample", app->pixels.at(0, 0).sample);
            edited += ygl::draw_camera_selection_widget(
                win, "camera", app->cam, app->scn, app->view);
            ygl::draw_value_widget(win, "fps", app->navigation_fps);
        }
        if (ygl::draw_header_widget(win, "params")) {
            if (ygl::draw_params_widgets(win, "", app->params)) {
                app->scene_updated = true;
            }
        }
        if (ygl::draw_header_widget(win, "image")) {
            ygl::draw_params_widgets(win, "", app->imparams);
            ygl::draw_imageinspect_widgets(
                win, "", app->img, {}, get_mouse_posf(win), app->imparams);
        }
        if (ygl::draw_header_widget(win, "scene")) {
            if (ygl::draw_scene_widgets(
                    win, "", app->scn, app->selection, app->update_list, {}))
                app->scene_updated = true;
        }
    }
    ygl::end_widgets(win);

    ygl::swap_buffers(win);
}

bool update(app_state* app) {
    if (app->scene_updated || !app->update_list.empty()) {
        ygl::trace_async_stop(app->async_threads, app->async_stop);
        app->rendering = false;

        // update BVH
        for (auto sel : app->update_list) {
            if (sel.ist || sel.sgr) {
                ygl::refit_bvh(app->bvh, app->scn, false);
            }
            if (sel.nde) {
                ygl::update_transforms(app->scn, 0);
                ygl::refit_bvh(app->bvh, app->scn, false);
            }
        }
        app->update_list.clear();

        // render preview
        auto pparams = app->params;
        pparams.nsamples = 1;
        pparams.filter = ygl::trace_filter_type::box;
        auto pimg =
            ygl::image4f((int)std::round(app->cam->aspect * app->preview_res),
                app->preview_res);
        auto ppixels = make_trace_pixels(pimg, pparams);
        ygl::trace_samples(app->scn, app->cam, app->bvh, app->lights, pimg,
            ppixels, 1, pparams);
        ygl::resize_image(pimg, app->img, ygl::resize_filter::box);
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
    auto win = ygl::make_window(
        app->img.width(), app->img.height(), "yitrace | " + app->filename, app);
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
        argc, argv, "yitrace", "Path trace images interactively");
    app->params = ygl::parse_params(parser, "", app->params);
    app->imparams = ygl::parse_params(parser, "", app->imparams);
    app->preview_res =
        ygl::parse_opt(parser, "--preview-res", "", "preview resolution", 32);
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
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
        app->params.shader != ygl::trace_shader_type::eyelight) {
        ygl::log_info("no lights presents, switching to eyelight shader");
        app->params.shader = ygl::trace_shader_type::eyelight;
    }

    // initialize rendering objects
    app->img =
        ygl::image4f((int)round(app->cam->aspect * app->params.resolution),
            app->params.resolution);
    app->pixels = ygl::make_trace_pixels(app->img, app->params);
    app->scene_updated = true;

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->async_threads, app->async_stop);
    delete app;

    // done
    return 0;
}
