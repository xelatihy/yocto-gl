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
    int resolution = 512;  // image vertical resolution
    int nsamples = 256;    // number of samples
    ygl::trace_type tracer = ygl::trace_type::pathtrace;  // tracer
    int nbounces = 8;                                     // max depth
    int seed = 7;                                         // seed
    float pixel_clamp = 100.0f;                           // pixel clamping

    // rendered image
    int width = 0, height = 512;
    std::vector<ygl::vec4f> img;

    // rendering state
    std::vector<ygl::rng_state> rngs;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    bool update_texture = true;

    // view image
    ygl::vec2f offset = {0, 0};
    float zoom = 1;
    float exposure = 0;
    float gamma = 2.2f;
    ygl::vec4f background = {0, 0, 0, 0};
    ygl::gltexture gl_txt = {};
    ygl::glimage_program gl_prog = {};
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    int cur_sample = 0;
    bool quiet = false;
    int preview_resolution = 64;
    bool navigation_fps = false;
    bool rendering = false;

    ~app_state() {
        if (scn) delete scn;
    }
};

auto trace_names = std::map<ygl::trace_type, std::string>{
    {ygl::trace_type::pathtrace, "pathtrace"},
    {ygl::trace_type::eyelight, "eyelight"},
    {ygl::trace_type::direct, "direct"},
    {ygl::trace_type::pathtrace_nomis, "pathtrace_nomis"},
    {ygl::trace_type::debug_normal, "debug_normal"},
    {ygl::trace_type::debug_albedo, "debug_albedo"},
    {ygl::trace_type::debug_texcoord, "debug_texcoord"},
    {ygl::trace_type::debug_frontfacing, "debug_frontfacing"},
};

namespace ygl {

bool draw_imgui_trace_inspector(glwindow* win, app_state* app) {
    auto edited = 0;
    edited += draw_imgui_dragbox(win, "resolution", app->resolution, 256, 4096);
    edited += draw_imgui_dragbox(win, "nsamples", app->nsamples, 16, 4096);
    edited += draw_imgui_combobox(win, "tracer", app->tracer, trace_names);
    edited += draw_imgui_dragbox(win, "nbounces", app->nbounces, 1, 10);
    edited += draw_imgui_dragbox(win, "seed", (int&)app->seed, 0, 1000);
    edited +=
        draw_imgui_dragbox(win, "preview", app->preview_resolution, 64, 1080);
    return edited;
}

}  // namespace ygl

void draw(ygl::glwindow* win, app_state* app) {
    // draw image
    auto window_size = get_glwindow_size(win);
    auto framebuffer_size = get_glframebuffer_size(win);
    ygl::set_glviewport(framebuffer_size);
    ygl::clear_glbuffers(app->background);
    ygl::draw_glimage(app->gl_prog, app->gl_txt, window_size, app->offset,
        app->zoom, app->exposure, app->gamma);

    if (ygl::begin_imgui_frame(win, "yitrace")) {
        if (ygl::draw_imgui_header(win, "trace")) {
            ygl::push_imgui_groupid(win, app);
            ygl::draw_imgui_label(win, "scene", app->filename);
            ygl::draw_imgui_label(
                win, "size", ygl::format("{} x {}", app->width, app->height));
            ygl::draw_imgui_label(
                win, "sample", std::to_string(app->cur_sample));
            if (ygl::draw_imgui_combobox(
                    win, "camera", app->cam, app->scn->cameras)) {
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
            if (ygl::draw_imgui_trace_inspector(win, app)) {
                app->update_list.push_back(nullptr);
            }
        }
        if (ygl::draw_imgui_header(win, "image")) {
            ygl::draw_imgui_dragbox(win, "offset", app->offset, -4096, 4096);
            ygl::draw_imgui_dragbox(win, "zoom", app->zoom, 0.01, 10);
            ygl::draw_imgui_colorbox(win, "background", app->background);
            ygl::draw_imgui_dragbox(win, "exposure", app->exposure, -5, 5);
            ygl::draw_imgui_dragbox(win, "gamma", app->gamma, 0.2f, 4);
            auto ij = ygl::get_glimage_coords(
                get_glmouse_posf(win), app->offset, app->zoom);
            ygl::draw_imgui_dragbox(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app->width && ij.y >= 0 &&
                ij.y < app->height) {
                ygl::draw_imgui_colorbox(
                    win, "pixel", app->img[ij.x + ij.y * app->width]);
            } else {
                auto zero4f_ = ygl::zero4f;
                ygl::draw_imgui_colorbox(win, "pixel", zero4f_);
            }
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
            trace_samples(app->scn, app->cam, pwidth, pheight, pimg, prngs, 0,
                1, app->tracer, app->nbounces);
            auto pratio = app->resolution / app->preview_resolution;
            for (auto j = 0; j < app->height; j++) {
                for (auto i = 0; i < app->width; i++) {
                    auto pi = i / pratio, pj = j / pratio;
                    app->img[i + j * app->width] = pimg[pi + pwidth * pj];
                }
            }
            // app->img = resize_image(
            //     pwidth, pheight, pimg, app->width, app->height,
            //     ygl::resize_filter::box);
        } else {
            for (auto& p : app->img) p = ygl::zero4f;
        }
        app->update_texture = true;

        // restart renderer
        app->rngs = ygl::make_rng_seq(app->img.size(), app->seed);
        ygl::trace_async_start(app->scn, app->cam, app->width, app->height,
            app->img, app->rngs, app->nsamples, app->tracer, app->nbounces,
            app->async_threads, app->async_stop, app->pixel_clamp, false,
            [win, app](int s, int j) {
                if (j % app->preview_resolution) return;
                app->update_texture = true;
#ifdef __APPLE__
                ygl::post_glwindow_event(win);
#endif
            });
    }
    if (app->update_texture) {
        ygl::update_gltexture(
            app->gl_txt, app->width, app->height, app->img, false, false, true);
        app->update_texture = false;
        return true;
    }
    return false;
}

void refresh(ygl::glwindow* win) {
    return draw(win, (app_state*)ygl::get_glwindow_user_pointer(win));
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win = ygl::make_glwindow(
        app->width, app->height, "yitrace | " + app->filename, app);
    ygl::set_glwindow_callbacks(win, nullptr, nullptr, refresh);

    // load textures
    app->gl_prog = ygl::make_glimage_program();
    ygl::update_gltexture(
        app->gl_txt, app->width, app->height, app->img, false, false, true);

    // init widget
    ygl::init_imgui(win);

    // loop
    while (!ygl::should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (ygl::handle_glcamera_navigation(
                win, app->cam, app->navigation_fps)) {
            app->update_list.push_back(app->cam);
        }
        ygl::handle_glscene_selection(win, app->scn, app->cam, app->resolution,
            app->offset, app->zoom, app->selection);

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
    app->resolution = ygl::parse_opt(parser, "--resolution", "-r",
        "Image vertical resolution.", app->resolution);
    app->nsamples = ygl::parse_opt(
        parser, "--nsamples", "-s", "Number of samples.", app->nsamples);
    app->tracer = ygl::parse_opt(
        parser, "--tracer", "-T", "Trace type.", trace_names, app->tracer);
    auto double_sided = ygl::parse_flag(
        parser, "--double-sided", "-D", "Force double sided rendering.", false);
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

    // fix scene
    ygl::log_info("adding missing scene elements");
    ygl::update_bbox(app->scn);
    ygl::add_missing_camera(app->scn);
    ygl::add_missing_names(app->scn);
    ygl::add_missing_tangent_space(app->scn);
    app->cam = app->scn->cameras[0];
    if (double_sided) {
        for (auto mat : app->scn->materials) mat->double_sided = true;
    }
    for (auto err : ygl::validate(app->scn)) ygl::log_warning(err);

    // build bvh
    ygl::log_info_begin("building bvh");
    ygl::update_bvh(app->scn);
    ygl::log_info_end();

    // init renderer
    ygl::log_info("initializing lights");
    ygl::update_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->tracer != ygl::trace_type::eyelight) {
        ygl::log_info("no lights presents, switching to eyelight shader");
        app->tracer = ygl::trace_type::eyelight;
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
