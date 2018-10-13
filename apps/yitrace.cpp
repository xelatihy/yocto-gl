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

#include "../yocto/ygl.h"
#include "../yocto/yglio.h"
#include "yglutils.h"
#include "ysceneui.h"

// Application state
struct app_state {
    // scene
    scene*    scn = nullptr;
    bvh_tree* bvh = nullptr;

    // rendering params
    std::string  filename   = "scene.json";
    std::string  imfilename = "out.obj";
    trace_params params     = {};

    // rendering state
    trace_state*  state  = nullptr;
    trace_lights* lights = nullptr;

    // view image
    vec2f                                      imcenter     = zero2f;
    float                                      imscale      = 1;
    bool                                       zoom_to_fit  = true;
    bool                                       widgets_open = false;
    void*                                      selection    = nullptr;
    std::vector<std::pair<std::string, void*>> update_list;
    bool                                       navigation_fps = false;
    bool                                       quiet          = false;
    int64_t                                    trace_start    = 0;
    uint                                       gl_txt         = 0;

    ~app_state() {
        if (scn) delete scn;
        if (bvh) delete bvh;
        if (state) delete state;
        if (lights) delete lights;
    }
};

void draw_glwidgets(glwindow* win) {
    auto app = (app_state*)get_user_pointer(win);
    begin_glwidgets_frame(win);
    if (begin_glwidgets_window(win, "yitrace")) {
        if (begin_header_glwidget(win, "scene")) {
            draw_label_glwidgets(win, "scene", app->filename);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "trace")) {
            draw_label_glwidgets(win, "image", "%d x %d @ %d",
                app->state->img.size().x, app->state->img.size().y,
                app->state->sample);
            auto cam_names = std::vector<std::string>();
            for (auto cam : app->scn->cameras) cam_names.push_back(cam->name);
            auto edited = 0;
            edited += draw_combobox_glwidget(
                win, "camera", app->params.camid, cam_names);
            edited += draw_slider_glwidget(
                win, "resolution", app->params.yresolution, 256, 4096);
            edited += draw_slider_glwidget(
                win, "nsamples", app->params.nsamples, 16, 4096);
            edited += draw_combobox_glwidget(
                win, "tracer", (int&)app->params.tracer, trace_type_names);
            edited += draw_slider_glwidget(
                win, "nbounces", app->params.nbounces, 1, 10);
            edited += draw_slider_glwidget(
                win, "seed", (int&)app->params.seed, 0, 1000);
            edited += draw_slider_glwidget(
                win, "pratio", app->params.preview_ratio, 1, 64);
            if (edited) app->update_list.push_back({"app", app});
            draw_label_glwidgets(win, "time/sample", "%0.3lf",
                (app->state->sample) ? (get_time() - app->trace_start) /
                                           (1000000000.0 * app->state->sample) :
                                       0.0);
            draw_slider_glwidget(win, "exposure", app->params.exposure, -5, 5);
            draw_checkbox_glwidget(win, "filmic", app->params.filmic);
            draw_checkbox_glwidget(win, "srgb", app->params.srgb);
            draw_slider_glwidget(win, "zoom", app->imscale, 0.1, 10);
            draw_checkbox_glwidget(win, "zoom to fit", app->zoom_to_fit);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "fps", app->navigation_fps);
            auto mouse_pos = get_glmouse_pos(win);
            auto ij        = get_image_coords(
                mouse_pos, app->imcenter, app->imscale, app->state->img.size());
            draw_dragger_glwidget(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app->state->img.size().x && ij.y >= 0 &&
                ij.y < app->state->img.size().y) {
                draw_coloredit_glwidget(
                    win, "pixel", app->state->img[{ij.x, ij.y}]);
            } else {
                auto zero4f_ = zero4f;
                draw_coloredit_glwidget(win, "pixel", zero4f_);
            }
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "navigate")) {
            draw_glwidgets_scene_tree(
                win, "", app->scn, app->selection, app->update_list, 200);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "inspec")) {
            draw_glwidgets_scene_inspector(
                win, "", app->scn, app->selection, app->update_list, 200);
            end_header_glwidget(win);
        }
    }
    end_glwidgets_frame(win);
}

void draw(glwindow* win) {
    auto app      = (app_state*)get_user_pointer(win);
    auto win_size = get_glwindow_size(win);
    auto fb_size  = get_glframebuffer_size(win);
    set_glviewport(fb_size);
    clear_glframebuffer(vec4f{0.8f, 0.8f, 0.8f, 1.0f});
    center_image4f(app->imcenter, app->imscale, app->state->display.size(),
        win_size, app->zoom_to_fit);
    if (!app->gl_txt) {
        app->gl_txt = make_gltexture(app->state->display, false, false, false);
    } else {
        update_gltexture(app->gl_txt, app->state->display, false, false, false);
    }
    draw_glimage(app->gl_txt, app->state->display.size(), win_size,
        app->imcenter, app->imscale);
    draw_glwidgets(win);
    swap_glbuffers(win);
}

bool update(app_state* app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    trace_async_stop(app->state);

    // update BVH
    for (auto& sel : app->update_list) {
        if (sel.first == "shape") {
            for (auto sid = 0; sid < app->scn->shapes.size(); sid++) {
                if (app->scn->shapes[sid] == sel.second) {
                    refit_bvh((shape*)sel.second, app->bvh->shape_bvhs[sid]);
                    break;
                }
            }
            refit_bvh(app->scn, app->bvh);
        }
        if (sel.first == "instance") { refit_bvh(app->scn, app->bvh); }
        if (sel.first == "node") {
            update_transforms(app->scn, 0);
            refit_bvh(app->scn, app->bvh);
        }
    }
    app->update_list.clear();

    delete app->state;
    app->trace_start = get_time();
    app->state       = make_trace_state(app->scn, app->params);
    trace_async_start(app->state, app->scn, app->bvh, app->lights, app->params);

    // updated
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto win_size = clamp(app->state->img.size(), 256, 1440);
    auto win      = make_glwindow(win_size, "yitrace", app, draw);

    // init widgets
    init_glwidgets(win);

    // loop
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_glwindow_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_glmouse_pos(win);
        auto mouse_left     = get_glmouse_left(win);
        auto mouse_right    = get_glmouse_right(win);
        auto alt_down       = get_glalt_key(win);
        auto shift_down     = get_glshift_key(win);
        auto widgets_active = get_glwidgets_active(win);

        // handle mouse and keyboard for navigation
        if ((mouse_left || mouse_right) && !alt_down && !widgets_active) {
            auto dolly  = 0.0f;
            auto pan    = zero2f;
            auto rotate = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto cam = app->scn->cameras.at(app->params.camid);
            camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            app->update_list.push_back({"camera", cam});
        }

        // selection
        if ((mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto ij = get_image_coords(mouse_pos, app->imcenter, app->imscale,
                {app->state->img.size().x, app->state->img.size().y});
            if (ij.x < 0 || ij.x >= app->state->img.size().x || ij.y < 0 ||
                ij.y >= app->state->img.size().y) {
                auto cam = app->scn->cameras.at(app->params.camid);
                auto ray = eval_camera_ray(
                    cam, ij, app->state->img.size(), {0.5f, 0.5f}, zero2f);
                auto isec = intersect_ray(app->scn, app->bvh, ray);
                if (isec.ist) app->selection = isec.ist;
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        if (!(mouse_left || mouse_right) && !widgets_active)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        process_glevents(win);
    }

    // clear
    delete_glwindow(win);
}

int main(int argc, char* argv[]) {
    // application
    auto app = new app_state();

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "progressive path tracing", "yitrace");
    app->params.camid       = parse_arg(parser, "--camera", 0, "Camera index.");
    app->params.yresolution = parse_arg(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->params.nsamples = parse_arg(
        parser, "--nsamples,-s", 4096, "Number of samples.");
    app->params.tracer = (trace_type)parse_arge(
        parser, "--tracer,-t", 0, "Tracer type.", trace_type_names);
    app->params.nbounces = parse_arg(
        parser, "--nbounces", 4, "Maximum number of bounces.");
    app->params.pixel_clamp = parse_arg(
        parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app->params.seed = parse_arg(
        parser, "--seed", 7, "Seed for the random number generators.");
    auto embree = parse_arg(parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = parse_arg(
        parser, "--double-sided", false, "Double-sided rendering.");
    auto add_skyenv = parse_arg(
        parser, "--add-skyenv", false, "Add missing env map");
    auto quiet = parse_arg(
        parser, "--quiet", false, "Print only errors messages");
    app->imfilename = parse_arg(
        parser, "--output-image,-o", "out.hdr", "Image filename");
    app->filename = parse_arg(
        parser, "scene", "scene.json", "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    if (!quiet) printf("loading scene %s\n", app->filename.c_str());
    auto load_start = get_time();
    app->scn = load_scene(app->filename);
    if(!app->scn)
        log_fatal("cannot load scene " + app->filename);
    if (!quiet)
        printf("loading in %s\n",
            format_duration(get_time() - load_start).c_str());

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    tesselate_subdivs(app->scn);

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (add_skyenv && app->scn->environments.empty()) {
        app->scn->environments.push_back(make_sky_environment("sky"));
        app->scn->textures.push_back(app->scn->environments.back()->ke_txt);
    }
    if (double_sided)
        for (auto mat : app->scn->materials) mat->double_sided = true;
    if (app->scn->cameras.empty())
        app->scn->cameras.push_back(
            make_bbox_camera("<view>", compute_bbox(app->scn)));
    add_missing_names(app->scn);
    for (auto& err : validate(app->scn)) printf("warning: %s\n", err.c_str());

    // build bvh
    if (!quiet) printf("building bvh\n");
    auto bvh_start = get_time();
    app->bvh       = build_bvh(app->scn, true, embree);
    if (!quiet)
        printf("building bvh in %s\n",
            format_duration(get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    app->lights = make_trace_lights(app->scn, app->params);

    // fix renderer type if no lights
    if (app->lights->lights.empty() && app->lights->environments.empty() &&
        app->params.tracer != trace_type::eyelight) {
        if (!quiet)
            printf("no lights presents, switching to eyelight shader\n");
        app->params.tracer = trace_type::eyelight;
    }

    // prepare renderer
    app->state = make_trace_state(app->scn, app->params);

    // initialize rendering objects
    if (!quiet) printf("starting async renderer\n");
    app->trace_start = get_time();
    trace_async_start(app->state, app->scn, app->bvh, app->lights, app->params);

    // run interactive
    run_ui(app);

    // cleanup
    trace_async_stop(app->state);
    delete app;

    // done
    return 0;
}
