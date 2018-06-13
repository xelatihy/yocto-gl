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

#include "../yocto/yocto_gl.h"
#include "../yocto/yocto_glio.h"
#include "../yocto/yocto_glutils.h"
#include "CLI11.hpp"
#include "yapp_ui.h"
using namespace std::literals;

#include <map>

// Application state
struct app_state {
    std::shared_ptr<ygl::scene> scn = nullptr;
    std::shared_ptr<ygl::camera> cam = nullptr;
    std::string filename = "scene.json"s;
    std::string imfilename = "out.obj"s;

    // rendering params
    int resolution = 512;                      // image vertical resolution
    int nsamples = 256;                        // number of samples
    std::string tracer = "pathtrace"s;         // tracer name
    ygl::trace_func tracef = ygl::trace_path;  // tracer
    int nbounces = 8;                          // max depth
    int seed = 7;                              // seed
    float pixel_clamp = 100.0f;                // pixel clamping
    bool double_sided = false;                 // double sided
    bool add_skyenv = false;                   // add sky environment

    // rendered image
    ygl::image4f img = {};

    // rendering state
    ygl::image<ygl::rng_state> rngs;
    std::vector<std::thread> async_threads;
    bool async_stop = false;
    int cur_sample = 0;

    // view image
    ygl::frame2f imframe = ygl::identity_frame2f;
    bool zoom_to_fit = true;
    float exposure = 0;
    float gamma = 2.2f;
    ygl::vec4f background = {0.8f, 0.8f, 0.8f, 0};
    ygl::gltexture gl_txt = {};
    ygl::glimage_program gl_prog = {};
    ygl::scene_selection selection = {};
    std::vector<ygl::scene_selection> update_list;
    bool quiet = false;
    int preview_resolution = 64;
    bool navigation_fps = false;
    bool rendering = false;
};

auto trace_names = std::vector<std::string>{
    "pathtrace",
    "direct",
    "environment",
    "eyelight",
    "pathtrace_nomis",
    "pathtrace_naive",
    "direct_nomis",
    "debug_normal",
    "debug_albedo",
    "debug_diffuse",
    "debug_specular",
    "debug_roughness",
    "debug_texcoord",
    "debug_frontfacing",
};

auto tracer_names = std::unordered_map<std::string, ygl::trace_func>{
    {"pathtrace", ygl::trace_path},
    {"direct", ygl::trace_direct},
    {"environment", ygl::trace_environment},
    {"eyelight", ygl::trace_eyelight},
    {"pathtrace-nomis", ygl::trace_path_nomis},
    {"pathtrace-naive", ygl::trace_path_naive},
    {"direct-nomis", ygl::trace_direct_nomis},
    {"debug-normal", ygl::trace_debug_normal},
    {"debug-albedo", ygl::trace_debug_albedo},
    {"debug-texcoord", ygl::trace_debug_texcoord},
    {"debug-frontfacing", ygl::trace_debug_frontfacing},
};

void draw(const std::shared_ptr<ygl::glwindow>& win,
    const std::shared_ptr<app_state>& app) {
    // update image
    ygl::center_glimage(app->imframe, {app->img.width(), app->img.height()},
        ygl::get_glwindow_size(win), app->zoom_to_fit);
    ygl::update_gltexture(app->gl_txt, app->img.width(), app->img.height(),
        app->img.pixels(), false, false, true, false);
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
            std::to_string(app->img.width()) + " x " +
                std::to_string(app->img.height()) + " @ " +
                std::to_string(app->cur_sample));
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
            app->tracef = tracer_names.at(app->tracer);
            edited +=
                draw_glwidgets_dragbox(win, "nbounces", app->nbounces, 1, 10);
            edited +=
                draw_glwidgets_dragbox(win, "seed", (int&)app->seed, 0, 1000);
            edited += draw_glwidgets_dragbox(
                win, "preview", app->preview_resolution, 64, 1080);
            if (edited) app->update_list.push_back(ygl::scene_selection());
            ygl::end_glwidgets_tree(win);
        }
        if (ygl::begin_glwidgets_tree(win, "view settings")) {
            ygl::draw_glwidgets_dragbox(win, "exposure", app->exposure, -5, 5);
            ygl::draw_glwidgets_dragbox(win, "gamma", app->gamma, 1, 3);
            ygl::draw_glwidgets_colorbox(win, "background", app->background);
            auto zoom = app->imframe.x.x;
            if (ygl::draw_glwidgets_dragbox(win, "zoom", zoom, 0.1, 10))
                app->imframe.x.x = app->imframe.y.y = zoom;
            ygl::draw_glwidgets_checkbox(win, "zoom to fit", app->zoom_to_fit);
            ygl::continue_glwidgets_line(win);
            ygl::draw_glwidgets_checkbox(win, "fps", app->navigation_fps);
            auto ij = ygl::get_glimage_coords(get_glwidnow_mouse_posf(win),
                app->imframe, {app->img.width(), app->img.height()});
            ygl::draw_glwidgets_dragbox(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app->img.width() && ij.y >= 0 &&
                ij.y < app->img.height()) {
                ygl::draw_glwidgets_colorbox(win, "pixel", app->img[ij]);
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
        if (ygl::begin_glwidgets_tree(win, "log")) {
            ygl::draw_glwidgets_log(win, 200);
            ygl::end_glwidgets_tree(win);
        }
    }
    ygl::end_glwidgets_frame(win);

    ygl::swap_glwindow_buffers(win);
}

bool update(const std::shared_ptr<ygl::glwindow>& win,
    const std::shared_ptr<app_state>& app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    ygl::trace_async_stop(app->async_threads, app->async_stop);

    // update BVH
    for (auto sel : app->update_list) {
        if (sel.as<ygl::shape>()) {
            if (!app->quiet) std::cout << "refit shape bvh\n";
            ygl::refit_bvh(sel.as<ygl::shape>());
        }
        if (sel.as<ygl::instance>()) {
            if (!app->quiet) std::cout << "refit scene bvh\n";
            ygl::refit_bvh(app->scn, false);
        }
        if (sel.as<ygl::node>()) {
            if (!app->quiet) std::cout << "refit scene bvh\n";
            ygl::update_transforms(app->scn, 0);
            ygl::refit_bvh(app->scn, false);
        }
    }
    app->update_list.clear();

    // render preview image
    if (app->preview_resolution) {
        auto pwidth =
            (int)std::round(app->preview_resolution * app->img.width() /
                            (float)app->img.height());
        auto pheight = app->preview_resolution;
        auto pimg = ygl::image4f{pwidth, pheight};
        auto prngs = ygl::make_trace_rngs(pwidth, pheight, 7);
        trace_samples(
            app->scn, app->cam, pimg, prngs, 0, 1, app->tracef, app->nbounces);
        auto pratio = app->resolution / app->preview_resolution;
        for (auto j = 0; j < app->img.height(); j++) {
            for (auto i = 0; i < app->img.width(); i++) {
                auto pi = ygl::clamp(i / pratio, 0, pwidth - 1),
                     pj = ygl::clamp(j / pratio, 0, pheight - 1);
                app->img[{i, j}] = pimg[{pi, pj}];
            }
        }
    } else {
        for (auto& p : app->img) p = ygl::zero4f;
    }

    // restart renderer
    if (!app->quiet) std::cout << "restart renderer\n";
    app->rngs =
        ygl::make_trace_rngs(app->img.width(), app->img.height(), app->seed);
    ygl::trace_async_start(app->scn, app->cam, app->img, app->rngs,
        app->nsamples, app->tracef, app->nbounces, app->async_threads,
        app->async_stop, app->cur_sample, app->pixel_clamp);

    // updated
    return true;
}

// run ui loop
void run_ui(const std::shared_ptr<app_state>& app) {
    // window
    auto win_width = app->img.width() + ygl::default_glwidgets_width;
    auto win_height = ygl::clamp(app->img.height(), 512, 1024);
    auto win = ygl::make_glwindow(win_width, win_height, "yitrace");
    ygl::set_glwindow_callbacks(
        win, nullptr, nullptr, [app, win]() { draw(win, app); });
    ygl::center_glimage(app->imframe, {app->img.width(), app->img.height()},
        ygl::get_glwindow_size(win), app->zoom_to_fit);

    // load textures
    app->gl_prog = ygl::make_glimage_program();
    ygl::update_gltexture(app->gl_txt, app->img.width(), app->img.height(),
        app->img.pixels(), false, false, true, false);

    // init widget
    ygl::init_glwidgets(win);

    // loop
    while (!ygl::should_glwindow_close(win)) {
        // handle mouse and keyboard for navigation
        if (app->navigation_fps) {
            if (ygl::handle_glcamera_fps(win, app->cam)) {
                app->update_list.push_back(app->cam);
            }
        } else {
            if (ygl::handle_glcamera_turntable(win, app->cam,
                    ygl::handle_glcamera_turntable_dist(
                        win, app->cam, app->scn, app->selection))) {
                app->update_list.push_back(app->cam);
            }
        }
        ygl::handle_glscene_selection(win, app->scn, app->cam,
            {app->img.width(), app->img.height()}, app->imframe,
            app->selection);

        // draw
        ygl::center_glimage(app->imframe, {app->img.width(), app->img.height()},
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
}

int main(int argc, char* argv[]) {
    // create empty scene
    auto app = std::make_shared<app_state>();

    // parse command line
    CLI::App parser("progressive path tracing", "yitrace");
    parser.add_option(
        "--resolution,-r", app->resolution, "Image vertical resolution.");
    parser.add_option("--nsamples,-s", app->nsamples, "Number of samples.");
    parser.add_option("--tracer,-t", app->tracer, "Trace type.")
        ->check([](const std::string& s) -> std::string {
            if (tracer_names.find(s) == tracer_names.end())
                throw CLI::ValidationError("unknown tracer name");
            return s;
        });
    parser.add_option(
        "--nbounces", app->nbounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", app->pixel_clamp, "Final pixel clamping.");
    parser.add_option(
        "--seed", app->seed, "Seed for the random number generators.");
    parser.add_option("--exposure,-e", app->exposure, "Hdr exposure");
    parser.add_flag(
        "--double-sided,-D", app->double_sided, "Double-sided rendering.");
    parser.add_flag("--add-skyenv,-E", app->add_skyenv, "add missing env map");
    parser.add_flag("--quiet,-q", app->quiet, "Print only errors messages");
    parser.add_option("--preview-resolution", app->preview_resolution,
        "Preview resolution for async rendering.");
    parser.add_option("--output-image,-o", app->imfilename, "Image filename");
    parser.add_option("scene", app->filename, "Scene filename")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) { return parser.exit(e); }
    app->tracef = tracer_names.at(app->tracer);

    // scene loading
    if (!app->quiet) std::cout << "loading scene" << app->filename << "\n";
    auto load_start = ygl::get_time();
    try {
        app->scn = ygl::load_scene(app->filename);
    } catch (const std::exception& e) {
        std::cout << "cannot load scene " << app->filename << "\n";
        std::cout << "error: " << e.what() << "\n";
        exit(1);
    }
    if (!app->quiet)
        std::cout << "loading in "
                  << ygl::format_duration(ygl::get_time() - load_start) << "\n";

    // tesselate
    if (!app->quiet) std::cout << "tesselating scene elements\n";
    ygl::update_tesselation(app->scn);

    // update bbox and transforms
    ygl::update_transforms(app->scn);
    ygl::update_bbox(app->scn);

    // add components
    if (!app->quiet) std::cout << "adding scene elements\n";
    if (app->add_skyenv && app->scn->environments.empty()) {
        app->scn->environments.push_back(ygl::make_sky_environment("sky"));
        app->scn->textures.push_back(app->scn->environments.back()->ke_txt);
    }
    if (app->double_sided)
        for (auto mat : app->scn->materials) mat->double_sided = true;
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            ygl::make_bbox_camera("<view>", app->scn->bbox));
    app->cam = app->scn->cameras[0];
    ygl::add_missing_names(app->scn);
    for (auto err : ygl::validate(app->scn))
        std::cout << "warning: " << err << "\n";

    // build bvh
    if (!app->quiet) std::cout << "building bvh\n";
    auto bvh_start = ygl::get_time();
    ygl::update_bvh(app->scn);
    if (!app->quiet)
        std::cout << "building bvh in "
                  << ygl::format_duration(ygl::get_time() - bvh_start) << "\n";

    // init renderer
    if (!app->quiet) std::cout << "initializing lights\n";
    ygl::update_lights(app->scn);

    // fix renderer type if no lights
    if (app->scn->lights.empty() && app->scn->environments.empty() &&
        app->tracer != "eyelight") {
        if (!app->quiet)
            std::cout << "no lights presents, switching to eyelight shader\n";
        app->tracer = "eyelight";
        app->tracef = ygl::trace_eyelight;
    }

    // initialize rendering objects
    if (!app->quiet) std::cout << "initializing tracer data\n";
    app->img = ygl::image4f{
        (int)round(app->resolution * app->cam->width / app->cam->height),
        app->resolution};
    app->rngs =
        ygl::make_trace_rngs(app->img.width(), app->img.height(), app->seed);
    app->update_list.push_back(app->scn);

    // run interactive
    run_ui(app);

    // cleanup
    ygl::trace_async_stop(app->async_threads, app->async_stop);

    // done
    return 0;
}
