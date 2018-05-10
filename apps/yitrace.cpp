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

#include <map>

// Application state
struct app_state {
    ygl::scene* scn = nullptr;
    ygl::camera* cam = nullptr;
    std::string filename;
    std::string imfilename;

    // rendering params
    int resolution = 512;                      // image vertical resolution
    int nsamples = 256;                        // number of samples
    ygl::trace_func tracer = ygl::trace_path;  // tracer
    int nbounces = 8;                          // max depth
    int seed = 7;                              // seed
    float pixel_clamp = 100.0f;                // pixel clamping

    // rendered image
    int width = 0, height = 512;
    std::vector<ygl::vec4f> img;

    // rendering state
    std::vector<ygl::rng_state> rngs;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    int cur_sample = 0;

    // view image
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = true;
    float exposure = 0;
    float gamma = 2.2f;
    ygl::vec4b background = {222, 222, 222, 0};
    ygl::gltexture gl_txt = {};
    ygl::glimage_program gl_prog = {};
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    bool quiet = false;
    int preview_resolution = 64;
    bool navigation_fps = false;
    bool rendering = false;

    ~app_state() {
        if (scn) delete scn;
    }
};

auto trace_names = std::map<ygl::trace_func, std::string>{
    {ygl::trace_path, "pathtrace"},
    {ygl::trace_direct, "direct"},
    {ygl::trace_eyelight, "eyelight"},
    {ygl::trace_path_nomis, "pathtrace_nomis"},
    {ygl::trace_path_naive, "pathtrace_naive"},
    {ygl::trace_direct_nomis, "direct_nomis"},
    {ygl::trace_debug_normal, "debug_normal"},
    {ygl::trace_debug_albedo, "debug_albedo"},
    {ygl::trace_debug_texcoord, "debug_texcoord"},
    {ygl::trace_debug_frontfacing, "debug_frontfacing"},
};

void draw(ygl::glwindow* win, app_state* app) {
    // update image
    ygl::update_gltexture(
        app->gl_txt, app->width, app->height, app->img, false, false, true);
    // draw image
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glwindow_framebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::clear_glbuffers(app->background);
    ygl::draw_glimage(app->gl_prog, app->gl_txt, window_size, app->imframe,
        app->exposure, app->gamma);

    if (ygl::begin_glwidgets_frame(win, "yitrace")) {
        ygl::draw_glwidgets_label(win, "scene", app->filename);
        ygl::draw_glwidgets_label(win, "image",
            ygl::format("{} x {} @ {} samples", app->width, app->height,
                app->cur_sample));
        if (ygl::begin_glwidgets_tree(win, "render settings")) {
            auto edited = 0;
            edited += ygl::draw_glwidgets_combobox(
                win, "camera", app->cam, app->scn->cameras);
            edited += draw_glwidgets_dragbox(
                win, "resolution", app->resolution, 256, 4096);
            edited += draw_glwidgets_dragbox(
                win, "nsamples", app->nsamples, 16, 4096);
            edited += draw_glwidgets_combobox(
                win, "tracer", app->tracer, trace_names);
            edited +=
                draw_glwidgets_dragbox(win, "nbounces", app->nbounces, 1, 10);
            edited +=
                draw_glwidgets_dragbox(win, "seed", (int&)app->seed, 0, 1000);
            edited += draw_glwidgets_dragbox(
                win, "preview", app->preview_resolution, 64, 1080);
            if (edited) app->update_list.push_back(nullptr);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "view settings")) {
            ygl::draw_glwidgets_dragbox(win, "exposure", app->exposure, -5, 5);
            ygl::draw_glwidgets_dragbox(win, "gamma", app->gamma, 0.2f, 4);
            ygl::draw_glwidgets_colorbox(win, "background", app->background);
            auto zoom = app->imframe.x.x;
            if (ygl::draw_glwidgets_dragbox(win, "zoom", zoom, 0.1, 10))
                app->imframe.x.x = app->imframe.y.y = zoom;
            ygl::draw_glwidgets_checkbox(win, "zoom to fit", app->zoom_to_fit);
            ygl::continue_glwidgets_line(win);
            ygl::draw_glwidgets_checkbox(win, "fps", app->navigation_fps);
            auto ij = ygl::get_glimage_coords(get_glwidnow_mouse_posf(win),
                app->imframe, {app->width, app->height});
            ygl::draw_glwidgets_dragbox(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app->width && ij.y >= 0 &&
                ij.y < app->height) {
                ygl::draw_glwidgets_colorbox(
                    win, "pixel", app->img[ij.x + ij.y * app->width]);
            } else {
                auto zero4f_ = ygl::zero4f;
                ygl::draw_glwidgets_colorbox(win, "pixel", zero4f_);
            }
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "scene tree")) {
            if (ygl::draw_glwidgets_button(win, "print stats"))
                ygl::print_stats(app->scn);
            ygl::draw_glwidgets_scene_tree(
                win, "", app->scn, app->selection, app->update_list);
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "scene object")) {
            ygl::draw_glwidgets_scene_inspector(
                win, "", app->scn, app->selection, app->update_list);
            ygl::end_glwidgets_tree(win);
        }
    }
    ygl::end_glwidgets_frame(win);

    ygl::swap_glwindow_buffers(win);
}

bool update(ygl::glwindow* win, app_state* app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->async_threads, app->async_stop);

    // update BVH
    for (auto sel : app->update_list) {
        if (sel.as<ygl::shape>()) { ygl::refit_bvh(sel.as<ygl::shape>()); }
        if (sel.as<ygl::instance>()) { ygl::refit_bvh(app->scn, false); }
        if (sel.as<ygl::node>()) {
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn, false);
        }
    }
    app->update_list.clear();

    // render preview image
    if (app->preview_resolution) {
        auto pwidth =
            (int)std::round(app->cam->aspect * app->preview_resolution);
        auto pheight = app->preview_resolution;
        auto pimg = std::vector<ygl::vec4f>(pwidth * pheight);
        auto prngs = ygl::make_rng_seq(pwidth * pheight, 7);
        trace_samples(app->scn, app->cam, pwidth, pheight, pimg, prngs, 0, 1,
            app->tracer, app->nbounces);
        auto pratio = app->resolution / app->preview_resolution;
        for (auto j = 0; j < app->height; j++) {
            for (auto i = 0; i < app->width; i++) {
                auto pi = ygl::clamp(i / pratio, 0, pwidth - 1),
                     pj = ygl::clamp(j / pratio, 0, pheight - 1);
                app->img[i + j * app->width] = pimg[pi + pwidth * pj];
            }
        }
    } else {
        for (auto& p : app->img) p = ygl::zero4f;
    }

    // restart renderer
    app->rngs = ygl::make_rng_seq(app->img.size(), app->seed);
    ygl::trace_async_start(app->scn, app->cam, app->width, app->height,
        app->img, app->rngs, app->nsamples, app->tracer, app->nbounces,
        app->async_threads, app->async_stop, app->cur_sample, app->pixel_clamp);

    // updated
    return true;
}

void refresh(ygl::glwindow* win) {
    auto app = (app_state*)ygl::get_glwindow_user_pointer(win);
    ygl::center_glimage(app->imframe, {app->width, app->height},
        ygl::get_glwindow_size(win), app->zoom_to_fit);
    draw(win, (app_state*)ygl::get_glwindow_user_pointer(win));
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win_width = app->width + ygl::default_glwidgets_width;
    auto win_height = ygl::clamp(app->height, 512, 1024);
    auto win = ygl::make_glwindow(win_width, win_height, "yitrace", app);
    ygl::set_glwindow_callbacks(win, nullptr, nullptr, refresh);
    ygl::center_glimage(app->imframe, {app->width, app->height},
        ygl::get_glwindow_size(win), app->zoom_to_fit);

    // load textures
    app->gl_prog = ygl::make_glimage_program();
    ygl::update_gltexture(
        app->gl_txt, app->width, app->height, app->img, false, false, true);

    // init widget
    ygl::init_glwidgets(win);

    // loop
    while (!ygl::should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (ygl::handle_glcamera_navigation(
                win, app->cam, app->navigation_fps)) {
            app->update_list.push_back(app->cam);
        }
        ygl::handle_glscene_selection(win, app->scn, app->cam, app->resolution,
            app->imframe, app->selection);

        // draw
        ygl::center_glimage(app->imframe, {app->width, app->height},
            ygl::get_glwindow_size(win), app->zoom_to_fit);
        draw(win, app);

        // update
        update(win, app);

        // event hadling
        if (!ygl::get_glwindow_mouse_button(win) &&
            !ygl::get_glwidgets_active(win)) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
        // event hadling
        ygl::poll_glwindow_events(win);
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
    app->resolution = ygl::parse_opt(parser, "--resolution", "-r",
        "Image vertical resolution.", app->resolution);
    app->nsamples = ygl::parse_opt(
        parser, "--nsamples", "-s", "Number of samples.", app->nsamples);
    app->tracer = ygl::parse_opte(
        parser, "--tracer", "-t", "Trace type.", trace_names, app->tracer);
    app->nbounces = ygl::parse_opt(
        parser, "--nbounces", "", "Maximum number of bounces.", app->nbounces);
    app->pixel_clamp = ygl::parse_opt(
        parser, "--pixel-clamp", "", "Final pixel clamping.", 100.0f);
    app->seed = ygl::parse_opt(
        parser, "--seed", "", "Seed for the random number generators.", 7);
    app->preview_resolution = ygl::parse_opt(parser, "--preview-resolution", "",
        "Preview resolution for async rendering.", app->preview_resolution);
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
    if (app->quiet) ygl::log_verbose() = false;

    // scene loading
    ygl::log_info_begin("loading scene {}", app->filename);
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (std::exception e) {
        ygl::log_fatal("cannot load scene {}", app->filename);
    }
    ygl::log_info_end();

    // tesselate
    ygl::log_info("tesselating scene elements");
    ygl::update_tesselation(app->scn);

    // fix scene
    ygl::log_info("adding missing scene elements");
    ygl::update_bbox(app->scn);
    ygl::add_missing_camera(app->scn);
    ygl::add_missing_names(app->scn);
    ygl::add_missing_tangent_space(app->scn);
    app->cam = app->scn->cameras[0];
    for (auto err : ygl::validate(app->scn)) ygl::log_warning(err);

    // build bvh
    ygl::log_info_begin("building bvh");
    ygl::update_bvh(app->scn);
    ygl::log_info_end();

    // init renderer
    ygl::log_info("initializing lights");
    ygl::update_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->scn->environments.empty() &&
        app->tracer != ygl::trace_eyelight) {
        ygl::log_info("no lights presents, switching to eyelight shader");
        app->tracer = ygl::trace_eyelight;
    }

    // initialize rendering objects
    ygl::log_info("initializing tracer data");
    app->width = (int)round(app->cam->aspect * app->resolution);
    app->height = app->resolution;
    app->img = std::vector<ygl::vec4f>(app->width * app->height);
    app->rngs = ygl::make_rng_seq(app->width * app->height, app->seed);
    app->update_list.push_back(app->scn);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->async_threads, app->async_stop);
    delete app;

    // done
    return 0;
}
