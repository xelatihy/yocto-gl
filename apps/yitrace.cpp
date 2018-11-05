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
    // loading params
    string filename       = "scene.json";
    string imfilename     = "out.obj";
    bool   use_embree_bvh = false;
    bool   double_sided   = false;
    bool   add_skyenv     = false;

    // scene
    yocto_scene scene = {};
    bvh_scene   bvh   = {};

    // rendering state
    trace_params                   params         = {};
    trace_lights                   lights         = {};
    image<rng_state>               trace_rngs     = {};
    image<vec4f>                   rendered_image = {};
    image<vec4f>                   display_image  = {};
    bool                           trace_stop     = false;
    int                            trace_sample   = 0;
    vector<thread>                 trace_threads  = {};
    concurrent_queue<image_region> trace_queue    = {};

    // view image
    vec2f                      image_center = zero2f;
    float                      image_scale  = 1;
    bool                       zoom_to_fit  = true;
    bool                       widgets_open = false;
    tuple<string, int>         selection    = {"", -1};
    vector<tuple<string, int>> update_list;
    bool                       navigation_fps = false;
    bool                       quiet          = false;
    int64_t                    trace_start    = 0;
    opengl_texture             gl_txt         = {};

    // app status
    bool   load_done = false, load_running = false;
    string status = "";
};

void start_rendering_async(app_state& app) {
    trace_async_stop(app.trace_threads, app.trace_stop, app.trace_queue);
    app.status      = "rendering image";
    app.trace_start = get_time();
    auto image_size = get_camera_image_size(
        app.scene.cameras[app.params.camera_id], app.params.image_size);
    app.rendered_image = make_image<vec4f>(image_size);
    app.display_image  = make_image<vec4f>(image_size);
    app.trace_rngs     = make_trace_rngs(image_size, app.params.random_seed);
    trace_async_start(app.rendered_image, app.display_image, app.scene, app.bvh,
        app.lights, app.trace_rngs, app.trace_threads, app.trace_stop,
        app.trace_sample, app.trace_queue, app.params);
}

void stop_rendering_async(app_state& app) {
    trace_async_stop(app.trace_threads, app.trace_stop, app.trace_queue);
}

bool load_scene_sync(app_state& app) {
    // scene loading
    app.status = "loading scene";
    if (!load_scene(app.filename, app.scene)) {
        log_fatal("cannot load scene " + app.filename);
        return false;
    }

    // tesselate
    app.status = "tesselating surfaces";
    tesselate_shapes_and_surfaces(app.scene);

    // add components
    if (app.add_skyenv && app.scene.environments.empty())
        add_sky_environment(app.scene);
    if (app.double_sided)
        for (auto& material : app.scene.materials) material.double_sided = true;
    add_missing_cameras(app.scene);
    add_missing_names(app.scene);
    log_validation_errors(app.scene);

    // build bvh
    app.status = "computing bvh";
    app.bvh    = make_scene_bvh(app.scene, true, app.use_embree_bvh);

    // init renderer
    app.status = "initializing lights";
    app.lights = make_trace_lights(app.scene);

    // fix renderer type if no lights
    if (empty(app.lights) && app.params.sample_tracer != trace_type::eyelight) {
        log_info("no lights presents, switching to eyelight shader\n");
        app.params.sample_tracer = trace_type::eyelight;
    }

    // set flags
    app.load_done    = true;
    app.load_running = false;
    app.status       = "loading done";

    // start rendering
    start_rendering_async(app);

    // done
    return false;
}

void load_scene_async(app_state& app) {
    if (app.load_running) {
        log_error("already loading");
        return;
    }
    app.load_done    = false;
    app.load_running = true;
    app.status       = "uninitialized";
    app.scene        = {};
    app.bvh          = {};
    app.lights       = {};
    auto load_thread = thread([&app]() { load_scene_sync(app); });
    load_thread.detach();
}

void draw_opengl_widgets(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    begin_opengl_widgets_frame(win);
    if (begin_opengl_widgets_window(win, "yitrace")) {
        if (begin_header_opengl_widget(win, "scene")) {
            draw_label_opengl_widget(win, "scene", get_filename(app.filename));
            if (draw_button_opengl_widget(win, "load")) {
                stop_rendering_async(app);
                load_scene_async(app);
            }
            draw_label_opengl_widget(win, "filename", app.filename);
            draw_label_opengl_widget(win, "status", app.status);
            end_header_opengl_widget(win);
        }
        if (begin_header_opengl_widget(win, "trace")) {
            draw_label_opengl_widget(win, "image", "%d x %d @ %d",
                app.rendered_image.size.x, app.rendered_image.size.y,
                app.trace_sample);
            auto cam_names = vector<string>();
            for (auto& camera : app.scene.cameras)
                cam_names.push_back(camera.name);
            auto edited = 0;
            if (app.load_done) {
                edited += draw_combobox_opengl_widget(
                    win, "camera", app.params.camera_id, cam_names);
            }
            edited += draw_slider_opengl_widget(
                win, "size", app.params.image_size, 256, 4096);
            edited += draw_slider_opengl_widget(
                win, "nsamples", app.params.num_samples, 16, 4096);
            edited += draw_combobox_opengl_widget(win, "tracer",
                (int&)app.params.sample_tracer, trace_type_names);
            edited += draw_slider_opengl_widget(
                win, "nbounces", app.params.max_bounces, 1, 10);
            edited += draw_slider_opengl_widget(
                win, "seed", (int&)app.params.random_seed, 0, 1000);
            edited += draw_slider_opengl_widget(
                win, "pratio", app.params.preview_ratio, 1, 64);
            if (edited) app.update_list.push_back({"app", -1});
            draw_label_opengl_widget(win, "time/sample", "%0.3lf",
                (app.trace_sample) ? (get_time() - app.trace_start) /
                                         (1000000000.0 * app.trace_sample) :
                                     0.0);
            draw_slider_opengl_widget(
                win, "exposure", app.params.display_exposure, -5, 5);
            draw_checkbox_opengl_widget(win, "filmic", app.params.display_filmic);
            draw_checkbox_opengl_widget(win, "srgb", app.params.display_srgb);
            draw_slider_opengl_widget(win, "zoom", app.image_scale, 0.1, 10);
            draw_checkbox_opengl_widget(win, "zoom to fit", app.zoom_to_fit);
            continue_opengl_widget_line(win);
            draw_checkbox_opengl_widget(win, "fps", app.navigation_fps);
            auto mouse_pos = get_opengl_mouse_pos(win);
            auto ij        = get_image_coords(mouse_pos, app.image_center,
                app.image_scale, app.rendered_image.size);
            draw_dragger_opengl_widget(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app.rendered_image.size.x && ij.y >= 0 &&
                ij.y < app.rendered_image.size.y) {
                draw_coloredit_opengl_widget(
                    win, "pixel", at(app.rendered_image, ij));
            } else {
                auto zero4f_ = zero4f;
                draw_coloredit_opengl_widget(win, "pixel", zero4f_);
            }
            end_header_opengl_widget(win);
        }
        if (app.load_done && begin_header_opengl_widget(win, "navigate")) {
            draw_opengl_widgets_scene_tree(
                win, "", app.scene, app.selection, app.update_list, 200);
            end_header_opengl_widget(win);
        }
        if (app.load_done && begin_header_opengl_widget(win, "inspec")) {
            draw_opengl_widgets_scene_inspector(
                win, "", app.scene, app.selection, app.update_list, 200);
            end_header_opengl_widget(win);
        }
    }
    end_opengl_widgets_frame(win);
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    set_glviewport(get_opengl_framebuffer_size(win));
    clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    if (app.load_done) {
        center_image(app.image_center, app.image_scale, app.display_image.size,
            win_size, app.zoom_to_fit);
        if (!app.gl_txt) {
            init_opengl_texture(
                app.gl_txt, app.display_image.size, false, false, false, false);
        } else {
            auto region = image_region{};
            while (app.trace_queue.try_pop(region))
                update_opengl_texture_region(
                    app.gl_txt, app.display_image, region, false);
        }
        set_glblending(true);
        draw_glimage_background(app.display_image.size, win_size,
            app.image_center, app.image_scale);
        draw_glimage(app.gl_txt, app.display_image.size, win_size,
            app.image_center, app.image_scale);
        set_glblending(false);
    }
    draw_opengl_widgets(win);
    swap_opengl_buffers(win);
}

bool update(app_state& app) {
    // exit if no updated
    if (!app.load_done || app.update_list.empty()) return false;

    // stop renderer
    stop_rendering_async(app);

    // update BVH
    for (auto& sel : app.update_list) {
        if (get<0>(sel) == "shape") {
            refit_shape_bvh(
                app.scene.shapes[get<1>(sel)], app.bvh.shape_bvhs[get<1>(sel)]);
            refit_scene_bvh(app.scene, app.bvh);
        }
        if (get<0>(sel) == "instance") {
            refit_scene_bvh(app.scene, app.bvh);
        }
        if (get<0>(sel) == "node") {
            update_transforms(app.scene, 0);
            refit_scene_bvh(app.scene, app.bvh);
        }
    }
    app.update_list.clear();

    // start rendering
    start_rendering_async(app);

    // updated
    return true;
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    trace_async_stop(app.trace_threads, app.trace_stop, app.trace_queue);
    app.filename = paths.front();
    load_scene_async(app);
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280, 720},
        "yitrace | " + get_filename(app.filename), &app, draw);
    set_drop_opengl_callback(win, drop_callback);

    // init widgets
    init_opengl_widgets(win);

    // loop
    auto mouse_pos = zero2f, last_pos = zero2f;
    while (!should_opengl_window_close(win)) {
        last_pos            = mouse_pos;
        mouse_pos           = get_opengl_mouse_pos(win);
        auto mouse_left     = get_opengl_mouse_left(win);
        auto mouse_right    = get_opengl_mouse_right(win);
        auto alt_down       = get_opengl_alt_key(win);
        auto shift_down     = get_opengl_shift_key(win);
        auto widgets_active = get_opengl_widgets_active(win);

        // handle mouse and keyboard for navigation
        if (app.load_done && (mouse_left || mouse_right) && !alt_down &&
            !widgets_active) {
            auto dolly  = 0.0f;
            auto pan    = zero2f;
            auto rotate = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
            auto& camera = app.scene.cameras.at(app.params.camera_id);
            camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            app.update_list.push_back({"camera", app.params.camera_id});
        }

        // selection
        if (app.load_done && (mouse_left || mouse_right) && alt_down &&
            !widgets_active) {
            auto ij = get_image_coords(mouse_pos, app.image_center,
                app.image_scale, app.rendered_image.size);
            if (ij.x < 0 || ij.x >= app.rendered_image.size.x || ij.y < 0 ||
                ij.y >= app.rendered_image.size.y) {
                auto& camera = app.scene.cameras.at(app.params.camera_id);
                auto  ray    = evaluate_camera_ray(
                    camera, ij, app.rendered_image.size, {0.5f, 0.5f}, zero2f);
                auto isec = intersect_scene(app.scene, app.bvh, ray);
                if (isec.instance_id >= 0)
                    app.selection = {"instance", isec.instance_id};
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        if (!(mouse_left || mouse_right) && !widgets_active)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        process_opengl_events(win);
    }

    // clear
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // application
    auto app = app_state();

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "progressive path tracing", "yitrace");
    app.params.camera_id   = parse_arg(parser, "--camera", 0, "Camera index.");
    app.params.image_size  = {0,
        parse_arg(parser, "--resolution,-r", 512, "Image vertical resolution.")};
    app.params.num_samples = parse_arg(
        parser, "--nsamples,-s", 4096, "Number of samples.");
    app.params.sample_tracer = parse_arge(parser, "--tracer,-t",
        trace_type::path, "Tracer type.", trace_type_names);
    app.params.max_bounces   = parse_arg(
        parser, "--nbounces", 4, "Maximum number of bounces.");
    app.params.pixel_clamp = parse_arg(
        parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app.params.random_seed = parse_arg(
        parser, "--seed", 7, "Seed for the random number generators.");
    app.use_embree_bvh = parse_arg(
        parser, "--embree", false, "Use Embree ratracer");
    app.double_sided = parse_arg(
        parser, "--double-sided", false, "Double-sided rendering.");
    app.add_skyenv = parse_arg(
        parser, "--add-skyenv", false, "Add missing environment map");
    app.imfilename = parse_arg(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    app.filename = parse_arg(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // load scene
    load_scene_async(app);

    // run interactive
    run_ui(app);

    // cleanup
    trace_async_stop(app.trace_threads, app.trace_stop, app.trace_queue);

    // done
    return 0;
}
