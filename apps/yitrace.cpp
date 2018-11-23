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

#include "../yocto/yocto_imageio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
#include "ysceneui.h"

// Application state
struct app_state {
    // loading options
    string filename   = "scene.json";
    string imfilename = "out.obj";

    // options
    load_scene_options  load_options  = {};
    build_bvh_options   bvh_options   = {};
    trace_image_options trace_options = {};
    float               exposure      = 0;
    bool                filmic        = false;
    bool                srgb          = true;
    int                 preview_ratio = 8;
    int                 image_width   = 0;
    int                 image_height  = 0;

    // scene
    yocto_scene scene = {};
    bvh_scene   bvh   = {};

    // rendering state
    trace_lights                   lights  = {};
    trace_state                    state   = {};
    image4f                        image   = {};
    image4f                        display = {};
    image4f                        preview = {};
    atomic<bool>                   trace_stop;
    atomic<int>                    trace_sample;
    vector<thread>                 trace_threads = {};
    concurrent_queue<image_region> trace_queue   = {};

    // view image
    vec2f                         image_center = zero2f;
    float                         image_scale  = 1;
    bool                          zoom_to_fit  = true;
    bool                          widgets_open = false;
    pair<type_index, int>         selection    = {typeid(void), -1};
    vector<pair<type_index, int>> update_list;
    bool                          navigation_fps  = false;
    bool                          quiet           = false;
    int64_t                       trace_start     = 0;
    opengl_texture                display_texture = {};

    // app status
    atomic<bool> load_done, load_running;
    string       status = "";
};

void stop_rendering_async(app_state& app) {
    trace_image_async_stop(app.trace_threads, app.trace_queue, app.trace_options);
}

void start_rendering_async(app_state& app) {
    stop_rendering_async(app);
    app.status       = "rendering image";
    app.trace_start  = get_time();
    app.trace_stop   = false;
    app.trace_sample = 0;

    tie(app.image_width, app.image_height) = get_camera_image_size(
        app.scene.cameras[app.trace_options.camera_id],
        app.trace_options.image_width, app.trace_options.image_height);
    app.image   = make_image(app.image_width, app.image_height, zero4f);
    app.display = make_image(app.image_width, app.image_height, zero4f);
    app.state   = make_trace_state(
        app.image_width, app.image_height, app.trace_options.random_seed);

    auto preview_options = app.trace_options;
    preview_options.image_width /= app.preview_ratio;
    preview_options.image_height /= app.preview_ratio;
    preview_options.num_samples = 1;
    app.preview = trace_image(app.scene, app.bvh, app.lights, preview_options);
    auto display_preview = tonemap_image(
        app.preview, app.exposure, app.filmic, app.srgb);
    auto large_preview = make_image(app.image_width, app.image_height, zero4f);
    for (auto j = 0; j < app.image_height; j++) {
        for (auto i = 0; i < app.image_width; i++) {
            auto pi = clamp(i / app.preview_ratio, 0, display_preview.width - 1),
                 pj = clamp(j / app.preview_ratio, 0, display_preview.height - 1);
            at(large_preview, i, j) = at(display_preview, pi, pj);
        }
    }
    app.preview = large_preview;
    app.trace_queue.push({0, 0, 0, 0});

    app.trace_options.cancel_flag = &app.trace_stop;
    trace_image_async_start(app.image, app.state, app.scene, app.bvh, app.lights,
        app.trace_threads, app.trace_sample, app.trace_queue, app.trace_options);
}

bool load_scene_sync(app_state& app) {
    // scene loading
    app.status = "loading scene";
    if (!load_scene(app.filename, app.scene, app.load_options)) {
        log_fatal("cannot load scene " + app.filename);
        return false;
    }

    // tesselate
    app.status = "tesselating surfaces";
    tesselate_shapes_and_surfaces(app.scene);

    // add components
    add_missing_cameras(app.scene);
    add_missing_names(app.scene);
    log_validation_errors(app.scene);

    // build bvh
    app.status = "computing bvh";
    app.bvh    = make_scene_bvh(app.scene, app.bvh_options);

    // init renderer
    app.status = "initializing lights";
    app.lights = make_trace_lights(app.scene);

    // fix renderer type if no lights
    if (empty(app.lights.instances) && empty(app.lights.environments) &&
        app.trace_options.sampler_type != trace_sampler_type::eyelight) {
        log_info("no lights presents, switching to eyelight shader\n");
        app.trace_options.sampler_type = trace_sampler_type::eyelight;
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
                app.image.width, app.image.height, (int)app.trace_sample);
            auto cam_names = vector<string>();
            for (auto& camera : app.scene.cameras)
                cam_names.push_back(camera.name);
            auto edited = 0;
            if (app.load_done) {
                edited += draw_combobox_opengl_widget(
                    win, "camera", app.trace_options.camera_id, cam_names);
            }
            edited += draw_slider_opengl_widget(
                win, "width", app.trace_options.image_width, 0, 4096);
            edited += draw_slider_opengl_widget(
                win, "height", app.trace_options.image_height, 0, 4096);
            edited += draw_slider_opengl_widget(
                win, "nsamples", app.trace_options.num_samples, 16, 4096);
            edited += draw_combobox_opengl_widget(win, "tracer",
                (int&)app.trace_options.sampler_type, trace_sampler_type_names);
            edited += draw_slider_opengl_widget(
                win, "nbounces", app.trace_options.max_bounces, 1, 10);
            edited += draw_checkbox_opengl_widget(
                win, "double sided", app.trace_options.double_sided);
            edited += draw_slider_opengl_widget(
                win, "seed", (int&)app.trace_options.random_seed, 0, 1000000);
            edited += draw_slider_opengl_widget(
                win, "pratio", app.preview_ratio, 1, 64);
            if (edited) app.update_list.push_back({typeid(app_state), -1});
            draw_label_opengl_widget(win, "time/sample", "%0.3lf",
                (app.trace_sample) ? (get_time() - app.trace_start) /
                                         (1000000000.0 * app.trace_sample) :
                                     0.0);
            draw_slider_opengl_widget(win, "exposure", app.exposure, -5, 5);
            draw_checkbox_opengl_widget(win, "filmic", app.filmic);
            draw_checkbox_opengl_widget(win, "srgb", app.srgb);
            draw_slider_opengl_widget(win, "zoom", app.image_scale, 0.1, 10);
            draw_checkbox_opengl_widget(win, "zoom to fit", app.zoom_to_fit);
            continue_opengl_widget_line(win);
            draw_checkbox_opengl_widget(win, "fps", app.navigation_fps);
            continue_opengl_widget_line(win);
            if (draw_button_opengl_widget(win, "print cams")) {
                for (auto& camera : app.scene.cameras) {
                    print("c {} {} {} {} {} {} {} {}\n", camera.name,
                        (int)camera.orthographic, camera.film_width,
                        camera.film_height, camera.focal_length,
                        camera.focus_distance, camera.lens_aperture,
                        camera.frame);
                }
            }
            auto mouse_pos = get_opengl_mouse_pos(win);
            auto ij        = get_image_coords(mouse_pos, app.image_center,
                app.image_scale, {app.image.width, app.image.height});
            draw_dragger_opengl_widget(win, "mouse", ij);
            if (ij.x >= 0 && ij.x < app.image.width && ij.y >= 0 &&
                ij.y < app.image.height) {
                draw_coloredit_opengl_widget(
                    win, "pixel", at(app.image, ij.x, ij.y));
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
    set_opengl_viewport(get_opengl_framebuffer_size(win));
    clear_opengl_lframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    if (app.load_done) {
        update_image_view(app.image_center, app.image_scale,
            {app.image.width, app.image.height}, win_size, app.zoom_to_fit);
        if (!app.display_texture) {
            if (app.image_width != 0 || app.image_height != 0) {
                init_opengl_texture(app.display_texture, app.image_width,
                    app.image_height, false, false, false, false);
            }
        } else {
            auto region = image_region{};
            auto size   = 0;
            while (app.trace_queue.try_pop(region)) {
                if (region.width == 0) {
                    update_opengl_texture(
                        app.display_texture, app.preview, false);
                    break;
                } else {
                    tonemap_image_region(app.display, region, app.image,
                        app.exposure, app.filmic, app.srgb);
                    update_opengl_texture_region(
                        app.display_texture, app.display, region, false);
                    size += region.width * region.height;
                    if (size >= app.image.width * app.image.height) break;
                }
            }
        }
        set_opengl_blending(true);
        draw_glimage_background(app.display_texture, win_size.x, win_size.y,
            app.image_center, app.image_scale);
        draw_glimage(app.display_texture, win_size.x, win_size.y,
            app.image_center, app.image_scale);
        set_opengl_blending(false);
    }
    draw_opengl_widgets(win);
    swap_opengl_buffers(win);
}

bool update(app_state& app) {
    // exit if no updated
    if (!app.load_done || empty(app.update_list)) return false;

    // stop renderer
    stop_rendering_async(app);

    // update BVH
    auto updated_instances = vector<int>{}, updated_shapes = vector<int>{},
         updated_surfaces = vector<int>{};
    for (auto [type, index] : app.update_list) {
        if (type == typeid(yocto_shape)) updated_shapes.push_back(index);
        if (type == typeid(yocto_instance)) updated_surfaces.push_back(index);
        if (type == typeid(yocto_scene_node))
            updated_instances.push_back(index);
    }
    if (!empty(updated_instances) || !empty(updated_shapes) ||
        !empty(updated_surfaces))
        refit_scene_bvh(app.scene, app.bvh, updated_instances, updated_shapes,
            updated_surfaces);
    app.update_list.clear();

    // start rendering
    start_rendering_async(app);

    // updated
    return true;
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    trace_image_async_stop(app.trace_threads, app.trace_queue, app.trace_options);
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
            auto& camera = app.scene.cameras.at(app.trace_options.camera_id);
            auto  dolly  = 0.0f;
            auto  pan    = zero2f;
            auto  rotate = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down)
                pan = (mouse_pos - last_pos) * camera.focus_distance / 200.0f;
            pan.x = -pan.x;
            update_camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            app.update_list.push_back(
                {typeid(yocto_camera), app.trace_options.camera_id});
        }

        // selection
        if (app.load_done && (mouse_left || mouse_right) && alt_down &&
            !widgets_active) {
            auto ij = get_image_coords(mouse_pos, app.image_center,
                app.image_scale, {app.image.width, app.image.height});
            if (ij.x < 0 || ij.x >= app.image.width || ij.y < 0 ||
                ij.y >= app.image.height) {
                auto& camera = app.scene.cameras.at(app.trace_options.camera_id);
                auto  ray    = evaluate_camera_ray(camera, ij,
                    {app.image.width, app.image.height}, {0.5f, 0.5f}, zero2f);
                auto  isec   = intersect_scene_bvh(app.bvh, ray);
                if (isec.instance_id >= 0)
                    app.selection = {typeid(yocto_instance), isec.instance_id};
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        process_opengl_events(win);
    }

    // clear
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // application
    app_state app{};
    app.trace_options.samples_per_batch = 1;

    // parse command line
    auto parser = make_cmdline_parser(
        argc, argv, "progressive path tracing", "yitrace");
    app.trace_options.camera_id = parse_argument(
        parser, "--camera", 0, "Camera index.");
    app.trace_options.image_width = parse_argument(
        parser, "--hres,-R", 1280, "Image horizontal resolution.");
    app.trace_options.image_height = parse_argument(
        parser, "--vres,-r", 720, "Image vertical resolution.");
    app.trace_options.num_samples = parse_argument(
        parser, "--nsamples,-s", 4096, "Number of samples.");
    app.trace_options.sampler_type = parse_argument(parser, "--tracer,-t",
        trace_sampler_type::path, "Tracer type.", trace_sampler_type_names);
    app.trace_options.max_bounces  = parse_argument(
        parser, "--nbounces", 8, "Maximum number of bounces.");
    app.trace_options.pixel_clamp = parse_argument(
        parser, "--pixel-clamp", 100, "Final pixel clamping.");
    app.trace_options.random_seed = parse_argument(
        parser, "--seed", 7, "Seed for the random number generators.");
    app.trace_options.environments_hidden = parse_argument(parser,
        "--env-hidden/--no-env-hidden", false,
        "Environments are hidden in renderer");
    auto no_parallel = parse_argument(parser, "--parallel/--no-parallel", false,
        "Disable parallel execution.");
    app.bvh_options.use_embree = parse_argument(
        parser, "--embree/--no-embree", false, "Use Embree ratracer");
    app.trace_options.double_sided = parse_argument(parser,
        "--double-sided/--no-double-sided", false, "Double-sided rendering.");
    app.imfilename                 = parse_argument(
        parser, "--output-image,-o", "out.hdr"s, "Image filename");
    app.filename = parse_argument(
        parser, "scene", "scene.json"s, "Scene filename", true);
    check_cmdline(parser);

    // fix parallel code
    if (no_parallel) {
        app.bvh_options.run_serially   = true;
        app.load_options.run_serially  = true;
        app.trace_options.run_serially = true;
    }

    // init app
    app.load_done    = false;
    app.load_running = false;

    // load scene
    load_scene_async(app);

    // run interactive
    run_ui(app);

    // cleanup
    trace_image_async_stop(app.trace_threads, app.trace_queue, app.trace_options);

    // done
    return 0;
}
