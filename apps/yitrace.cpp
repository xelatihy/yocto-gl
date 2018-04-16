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

#include "../yocto/yocto_bvh.h"
#include "../yocto/yocto_glutils.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
#include "yapp_ui.h"
using namespace std::literals;

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* view = nullptr;
    ygl::camera* cam = nullptr;
    ygl::bvh_tree* bvh = nullptr;
    std::string filename;
    std::string imfilename;
    ygl::load_options loadopts;
    int resolution = 512;
    ygl::image4f img;
    std::vector<ygl::trace_pixel> pixels;
    ygl::trace_lights lights;
    ygl::trace_params params;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    bool update_texture = true;
    bool navigation_fps = false;
    bool rendering = false;
    ygl::glimage_params imparams = {};
    ygl::gltexture trace_texture = {};
    ygl::glimage_program gl_prog = {};
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    int cur_sample = 0;
    bool quiet = false;

    ~app_state() {
        if (scn) delete scn;
        if (view) delete view;
        if (bvh) delete bvh;
    }
};

namespace ygl {

bool draw_imgui_trace_inspector(
    glwindow* win, const std::string& lbl, trace_params& params) {
    auto edited = 0;
    edited +=
        draw_imgui_dragbox(win, "resolution", params.resolution, 256, 4096);
    edited += draw_imgui_dragbox(win, "nsamples", params.nsamples, 16, 4096);
    edited +=
        draw_imgui_combobox(win, "tracer", params.tracer, trace_type_names());
    edited += draw_imgui_checkbox(win, "notransmission", params.notransmission);
    edited += draw_imgui_checkbox(win, "double_sided", params.double_sided);
    edited += draw_imgui_colorbox(win, "ambient", params.ambient);
    edited +=
        draw_imgui_checkbox(win, "envmap_invisible", params.envmap_invisible);
    edited += draw_imgui_dragbox(win, "min_depth", params.min_depth, 1, 10);
    edited += draw_imgui_dragbox(win, "max_depth", params.max_depth, 1, 10);
    edited +=
        draw_imgui_dragbox(win, "pixel_clamp", params.pixel_clamp, 10, 1000);
    edited +=
        draw_imgui_dragbox(win, "ray_eps", params.ray_eps, 0.0001f, 0.001f);
    edited += draw_imgui_checkbox(win, "parallel", params.parallel);
    edited += draw_imgui_dragbox(win, "seed", (int&)params.seed, 0, 1000);
    edited +=
        draw_imgui_dragbox(win, "preview", params.preview_resolution, 64, 1080);
    return edited;
}

}

void draw(ygl::glwindow* win, app_state* app) {
    // draw image
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glframebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::draw_glimage(
        app->gl_prog, app->trace_texture, window_size, app->imparams);

    if (ygl::begin_imgui_frame(win, "yitrace")) {
        if (ygl::draw_imgui_header(win, "trace")) {
            ygl::push_imgui_groupid(win, app);
            ygl::draw_imgui_label(win, "scene", app->filename);
            ygl::draw_imgui_label(
                win, "size", ygl::format("{} x {}", app->img.width, app->img.height));
            ygl::draw_imgui_label(win, "sample", std::to_string(app->cur_sample));
            if (ygl::draw_imgui_camera_selector(
                    win, "camera", app->cam, app->scn, app->view)) {
                app->update_list.push_back(app->cam);
            }
            if (ygl::draw_imgui_camera_inspector(win, "camera", app->cam)) {
                app->update_list.push_back(app->cam);
            }
            ygl::draw_imgui_checkbox(win, "fps", app->navigation_fps);
            if (ygl::draw_imgui_button(win, "print stats"))
                ygl::print_stats(app->scn);
            ygl::pop_imgui_groupid(win);
        }
        if (ygl::draw_imgui_header(win, "params")) {
            if (ygl::draw_imgui_trace_inspector(win, "", app->params)) {
                app->update_list.push_back(&app->params);
            }
        }
        if (ygl::draw_imgui_header(win, "image")) {
            ygl::draw_imgui_stdimage_inspector(win, "", app->imparams);
            ygl::draw_imgui_image_inspector(
                win, "", app->img, {}, get_glmouse_posf(win), app->imparams);
        }
        if (ygl::draw_imgui_header(win, "scene")) {
            ygl::draw_imgui_scene_tree(
                win, "", app->scn, app->selection, app->update_list);
        }
        if (ygl::draw_imgui_header(win, "inspect")) {
            ygl::draw_imgui_scene_inspector(
                win, "", app->scn, app->selection, app->update_list);
        }
    }
    ygl::end_imgui_frame(win);

    ygl::swap_glwindow_buffers(win);
}

bool update(ygl::glwindow* win, app_state* app) {
    if (!app->update_list.empty()) {
        ygl::trace_async_stop(app->async_threads, app->async_stop);

        // update BVH
        for (auto sel : app->update_list) {
            if (sel.is<ygl::shape>()) {
                ygl::refit_shape_bvh(app->bvh, sel.get<ygl::shape>(), false);
            }
            if (sel.is<ygl::node>()) {
                ygl::update_transforms(app->scn, 0);
                ygl::refit_scene_bvh(app->bvh, app->scn, false);
            }
        }
        app->update_list.clear();

        ygl::trace_async_start(app->scn, app->cam, app->bvh, app->lights,
            app->img, app->pixels, app->async_threads, app->async_stop,
            app->params, [win, app](int s, int j) {
                if (j % app->params.preview_resolution) return;
                app->update_texture = true;
#ifdef __APPLE__
                ygl::post_glwindow_event(win);
#endif
            });
    }
    if (app->update_texture) {
        ygl::update_gltexture(app->trace_texture, app->img, false, false, true);
        app->update_texture = false;
        return true;
    }
    return false;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = ygl::make_glwindow(
        app->img.width, app->img.height, "yitrace | " + app->filename);
    ygl::set_glwindow_callbacks(
        win, nullptr, nullptr, [win, app]() { draw(win, app); });

    // load textures
    app->gl_prog = ygl::make_glimage_program();
    ygl::update_gltexture(app->trace_texture, app->img, false, false, true);

    // init widget
    ygl::init_imgui(win);

    // loop
    while (!ygl::should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->cam == app->view) {
            if (ygl::handle_glcamera_navigation(
                    win, app->view, app->navigation_fps)) {
                app->update_list.push_back(app->view);
            }
        }
        ygl::handle_glscene_selection(win, app->scn, app->cam, app->bvh,
            app->params.resolution, app->imparams, app->selection);

        // draw
        draw(win, app);

        // update
        update(win, app);

#ifdef __APPLE__
        // event hadling
        if (ygl::get_glmouse_button(win) || ygl::get_imgui_active(win)) {
            ygl::poll_glwindow_events(win);
        } else {
            ygl::wait_glwindow_events(win);
        }
#else
        // event hadling
        ygl::poll_glwindow_events(win);
#endif
    }

    // cleanup
    delete win;
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = new app_state();

    // parse command line
    auto parser = ygl::make_parser(
        argc, argv, "yitrace", "Path trace images interactively");
    app->params.resolution = ygl::parse_opt(parser, "--resolution", "-r",
        "Image vertical resolution.", app->params.resolution);
    app->params.nsamples = ygl::parse_opt(
        parser, "--nsamples", "-s", "Number of samples.", app->params.nsamples);
    app->params.tracer = ygl::parse_opt(parser, "--tracer", "-T", "Trace type.",
        ygl::trace_type_names(), app->params.tracer);
    app->params.notransmission = ygl::parse_opt(parser, "--notransmission", "",
        "Whether to test transmission in shadows.", app->params.notransmission);
    app->params.double_sided = ygl::parse_opt(parser, "--double-sided", "-D",
        "Force double sided rendering.", app->params.double_sided);
    app->params.ambient = ygl::parse_opt(
        parser, "--ambient", "", "Ambient lighting.", app->params.ambient);
    app->params.envmap_invisible = ygl::parse_opt(parser, "--envmap-invisible",
        "", "View environment map.", app->params.envmap_invisible);
    app->params.min_depth = ygl::parse_opt(
        parser, "--min-depth", "", "Minimum ray depth.", app->params.min_depth);
    app->params.max_depth = ygl::parse_opt(
        parser, "--max-depth", "", "Maximum ray depth.", app->params.max_depth);
    app->params.pixel_clamp = ygl::parse_opt(parser, "--pixel-clamp", "",
        "Final pixel clamping.", app->params.pixel_clamp);
    app->params.ray_eps = ygl::parse_opt(parser, "--ray-eps", "",
        "Ray intersection epsilon.", app->params.ray_eps);
    app->params.parallel = ygl::parse_opt(
        parser, "--parallel", "", "Parallel execution.", app->params.parallel);
    app->params.seed = ygl::parse_opt(parser, "--seed", "",
        "Seed for the random number generators.", app->params.seed);
    app->params.preview_resolution = ygl::parse_opt(parser,
        "--preview-resolution", "", "Preview resolution for async rendering.",
        app->params.preview_resolution);
    app->params.batch_size = ygl::parse_opt(parser, "--batch-size", "",
        "Sample batch size.", app->params.batch_size);
    app->quiet =
        ygl::parse_flag(parser, "--quiet", "-q", "Print only errors messages");
    app->imfilename = ygl::parse_opt(
        parser, "--output-image", "-o", "Image filename", "out.hdr"s);
    app->filename = ygl::parse_arg(parser, "scene", "Scene filename", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // setup logger
    if (app->quiet) ygl::get_default_logger()->verbose = false;

    // scene loading
    try {
        ygl::log_info("loading scene {}", app->filename);
        app->scn = load_scene(app->filename, app->loadopts);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }

    // add elements
    ygl::add_names(app->scn);
    ygl::add_tangent_space(app->scn);

    // validate
    ygl::validate(app->scn, false, true);

    // view camera
    app->view = ygl::make_view_camera(app->scn, 0);
    app->cam = app->view;

    // build bvh
    ygl::log_info("building bvh");
    app->bvh = ygl::build_scene_bvh(app->scn);

    // init renderer
    ygl::log_info("initializing tracer");
    app->lights = ygl::make_trace_lights(app->scn);

    // fix renderer type if no lights
    if (app->lights.lights.empty() &&
        app->params.tracer != ygl::trace_type::eyelight) {
        ygl::log_info("no lights presents, switching to eyelight shader");
        app->params.tracer = ygl::trace_type::eyelight;
    }

    // initialize rendering objects
    app->img =
        ygl::make_image4f((int)round(app->cam->aspect * app->params.resolution),
            app->params.resolution);
    app->pixels = ygl::make_trace_pixels(app->img, app->params);
    app->update_list.push_back(app->scn);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->async_threads, app->async_stop);
    delete app;

    // done
    return 0;
}
