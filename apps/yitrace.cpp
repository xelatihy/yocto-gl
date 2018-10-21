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
    unique_ptr<yocto_scene> scene = nullptr;
    unique_ptr<bvh_tree>    bvh   = nullptr;

    // rendering params
    string       filename   = "scene.json";
    string       imfilename = "out.obj";
    trace_params params     = {};

    // rendering state
    unique_ptr<trace_state>  state  = nullptr;
    unique_ptr<trace_lights> lights = nullptr;

    // view image
    vec2f                        imcenter     = zero2f;
    float                        imscale      = 1;
    bool                         zoom_to_fit  = true;
    bool                         widgets_open = false;
    void*                        selection    = nullptr;
    vector<tuple<string, void*>> update_list;
    bool                         navigation_fps = false;
    bool                         quiet          = false;
    int64_t                      trace_start    = 0;
    gltexture                    gl_txt         = {};
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
                app->state->rendered_image.width,
                app->state->rendered_image.height, app->state->current_sample);
            auto cam_names = vector<string>();
            for (auto camera : app->scene->cameras)
                cam_names.push_back(camera->name);
            auto edited = 0;
            edited += draw_combobox_glwidget(
                win, "camera", app->params.camera_id, cam_names);
            edited += draw_slider_glwidget(
                win, "resolution", app->params.vertical_resolution, 256, 4096);
            edited += draw_slider_glwidget(
                win, "nsamples", app->params.num_samples, 16, 4096);
            edited += draw_combobox_glwidget(win, "tracer",
                (int&)app->params.sample_tracer, trace_type_names);
            edited += draw_slider_glwidget(
                win, "nbounces", app->params.max_bounces, 1, 10);
            edited += draw_slider_glwidget(
                win, "seed", (int&)app->params.random_seed, 0, 1000);
            edited += draw_slider_glwidget(
                win, "pratio", app->params.preview_ratio, 1, 64);
            if (edited) app->update_list.push_back({"app", app});
            draw_label_glwidgets(win, "time/sample", "%0.3lf",
                (app->state->current_sample) ?
                    (get_time() - app->trace_start) /
                        (1000000000.0 * app->state->current_sample) :
                    0.0);
            draw_slider_glwidget(
                win, "exposure", app->params.display_exposure, -5, 5);
            draw_checkbox_glwidget(win, "filmic", app->params.display_filmic);
            draw_checkbox_glwidget(win, "srgb", app->params.display_srgb);
            draw_slider_glwidget(win, "zoom", app->imscale, 0.1, 10);
            draw_checkbox_glwidget(win, "zoom to fit", app->zoom_to_fit);
            continue_glwidgets_line(win);
            draw_checkbox_glwidget(win, "fps", app->navigation_fps);
            auto mouse_pos = get_glmouse_pos(win);
            auto ij = get_image_coords(mouse_pos, app->imcenter, app->imscale,
                {app->state->rendered_image.width,
                    app->state->rendered_image.height});
            draw_dragger_glwidget(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app->state->rendered_image.width &&
                ij.y >= 0 && ij.y < app->state->rendered_image.height) {
                draw_coloredit_glwidget(
                    win, "pixel", at(app->state->rendered_image, ij.x, ij.y));
            } else {
                auto zero4f_ = zero4f;
                draw_coloredit_glwidget(win, "pixel", zero4f_);
            }
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "navigate")) {
            draw_glwidgets_scene_tree(win, "", app->scene.get(), app->selection,
                app->update_list, 200);
            end_header_glwidget(win);
        }
        if (begin_header_glwidget(win, "inspec")) {
            draw_glwidgets_scene_inspector(win, "", app->scene.get(),
                app->selection, app->update_list, 200);
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
    center_image(app->imcenter, app->imscale,
        {app->state->display_image.width, app->state->display_image.height},
        win_size, app->zoom_to_fit);
    if (!app->gl_txt) {
        app->gl_txt = make_gltexture(
            app->state->display_image, false, false, false);
    } else {
        update_gltexture(
            app->gl_txt, app->state->display_image, false, false, false);
    }
    draw_glimage(app->gl_txt,
        {app->state->display_image.width, app->state->display_image.height},
        win_size, app->imcenter, app->imscale);
    draw_glwidgets(win);
    swap_glbuffers(win);
}

bool update(app_state* app) {
    // exit if no updated
    if (app->update_list.empty()) return false;

    // stop renderer
    trace_async_stop(app->state.get());

    // update BVH
    for (auto& sel : app->update_list) {
        if (get<0>(sel) == "shape") {
            for (auto sid = 0; sid < app->scene->shapes.size(); sid++) {
                if (app->scene->shapes[sid] == get<1>(sel)) {
                    refit_shape_bvh(
                        (yocto_shape*)get<1>(sel), app->bvh->shape_bvhs[sid]);
                    break;
                }
            }
            refit_scene_bvh(app->scene.get(), app->bvh.get());
        }
        if (get<0>(sel) == "instance") {
            refit_scene_bvh(app->scene.get(), app->bvh.get());
        }
        if (get<0>(sel) == "node") {
            update_transforms(app->scene.get(), 0);
            refit_scene_bvh(app->scene.get(), app->bvh.get());
        }
    }
    app->update_list.clear();

    app->state       = {};
    app->trace_start = get_time();
    app->state       = unique_ptr<trace_state>(
        make_trace_state(app->scene.get(), app->params));
    trace_async_start(app->state.get(), app->scene.get(), app->bvh.get(),
        app->lights.get(), app->params);

    // updated
    return true;
}

// run ui loop
void run_ui(app_state* app) {
    // window
    auto width  = clamp(app->state->rendered_image.width, 256, 1440);
    auto height = clamp(app->state->rendered_image.height, 256, 1440);
    auto win    = make_glwindow(width, height, "yitrace", app, draw);

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
            auto camera = app->scene->cameras.at(app->params.camera_id);
            camera_turntable(
                camera->frame, camera->focus_distance, rotate, dolly, pan);
            app->update_list.push_back({"camera", camera});
        }

        // selection
        if ((mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto ij = get_image_coords(mouse_pos, app->imcenter, app->imscale,
                {app->state->rendered_image.width,
                    app->state->rendered_image.height});
            if (ij.x < 0 || ij.x >= app->state->rendered_image.width ||
                ij.y < 0 || ij.y >= app->state->rendered_image.height) {
                auto camera = app->scene->cameras.at(app->params.camera_id);
                auto ray    = evaluate_camera_ray(camera, ij,
                    {app->state->rendered_image.width,
                        app->state->rendered_image.height},
                    {0.5f, 0.5f}, zero2f);
                auto isec = intersect_sene(app->scene.get(), app->bvh.get(), ray);
                if (isec.instance) app->selection = isec.instance;
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
    app->params.camera_id = parse_arg(parser, "--camera", 0, "Camera index.");
    app->params.vertical_resolution = parse_arg(
        parser, "--resolution,-r", 512, "Image vertical resolution.");
    app->params.num_samples = parse_arg(
        parser, "--nsamples,-s", 4096, "Number of samples.");
    app->params.sample_tracer = parse_arge(parser, "--tracer,-t",
        trace_type::path, "Tracer type.", trace_type_names);
    app->params.max_bounces   = parse_arg(
        parser, "--nbounces", 4, "Maximum number of bounces.");
    app->params.pixel_clamp = parse_arg(
        parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app->params.random_seed = parse_arg(
        parser, "--seed", 7, "Seed for the random number generators.");
    auto embree = parse_arg(parser, "--embree", false, "Use Embree ratracer");
    auto double_sided = parse_arg(
        parser, "--double-sided", false, "Double-sided rendering.");
    auto add_skyenv = parse_arg(
        parser, "--add-skyenv", false, "Add missing environment map");
    auto quiet = parse_arg(
        parser, "--quiet", false, "Print only errors messages");
    app->imfilename = parse_arg(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    app->filename = parse_arg(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // scene loading
    if (!quiet) printf("loading scene %s\n", app->filename.c_str());
    auto load_start = get_time();
    app->scene      = unique_ptr<yocto_scene>(load_scene(app->filename));
    if (!app->scene) log_fatal("cannot load scene " + app->filename);
    if (!quiet)
        printf("loading in %s\n",
            format_duration(get_time() - load_start).c_str());

    // tesselate
    if (!quiet) printf("tesselating scene elements\n");
    tesselate_subdivs(app->scene.get());

    // add components
    if (!quiet) printf("adding scene elements\n");
    if (add_skyenv && app->scene->environments.empty()) {
        app->scene->environments.push_back(make_sky_environment("sky"));
        app->scene->textures.push_back(
            app->scene->environments.back()->emission_texture);
    }
    if (double_sided)
        for (auto mat : app->scene->materials) mat->double_sided = true;
    if (app->scene->cameras.empty())
        app->scene->cameras.push_back(
            make_bbox_camera("<view>", compute_scene_bounds(app->scene.get())));
    add_missing_names(app->scene.get());
    for (auto& err : validate_scene(app->scene.get()))
        printf("warning: %s\n", err.c_str());

    // build bvh
    if (!quiet) printf("building bvh\n");
    auto bvh_start = get_time();
    app->bvh       = unique_ptr<bvh_tree>(
        make_scene_bvh(app->scene.get(), true, embree));
    if (!quiet)
        printf("building bvh in %s\n",
            format_duration(get_time() - bvh_start).c_str());

    // init renderer
    if (!quiet) printf("initializing lights\n");
    app->lights = unique_ptr<trace_lights>(
        make_trace_lights(app->scene.get(), app->params));

    // fix renderer type if no lights
    if (!app->lights && app->params.sample_tracer != trace_type::eyelight) {
        if (!quiet)
            printf("no lights presents, switching to eyelight shader\n");
        app->params.sample_tracer = trace_type::eyelight;
    }

    // prepare renderer
    app->state = unique_ptr<trace_state>(
        make_trace_state(app->scene.get(), app->params));

    // initialize rendering objects
    if (!quiet) printf("starting async renderer\n");
    app->trace_start = get_time();
    trace_async_start(app->state.get(), app->scene.get(), app->bvh.get(),
        app->lights.get(), app->params);

    // run interactive
    run_ui(app);

    // cleanup
    trace_async_stop(app->state.get());

    // done
    return 0;
}
