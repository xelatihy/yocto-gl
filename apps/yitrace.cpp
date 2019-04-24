//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "../yocto/yocto_utils.h"
#include "yocto_opengl.h"
#include "ysceneui.h"
using namespace yocto;

#include "ext/CLI11.hpp"

#include <map>

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};  // namespace yocto

// Application task
enum struct app_task_type { none, load, save, bvh, lights, render, save_scene };

struct app_task {
    app_task_type                  type;
    future<void>                   result;
    atomic<bool>                   stop;
    concurrent_queue<image_region> queue;

    app_task(app_task_type type) : type{type}, result{}, stop{false}, queue{} {}
    ~app_task() {
        stop = true;
        if (result.valid()) {
            try {
                result.get();
            } catch (...) {
            }
        }
    }
};

// Application scene
struct app_scene {
    // loading options
    string filename   = "scene.json";
    string imfilename = "out.obj";

    // options
    load_scene_options    load_options          = {};
    bvh_build_options     bvh_options           = {};
    trace_image_options   trace_options         = {};
    tonemap_image_options tonemap_options = {};
    int                   preview_ratio         = 8;
    vec2i                 image_size            = {0, 0};

    // scene
    yocto_scene scene      = {};
    bvh_scene   bvh        = {};
    bool        add_skyenv = false;

    // rendering state
    trace_lights                   lights  = {};
    trace_state                    state   = {};
    image<vec4f>                   render  = {};
    image<vec4f>                   display = {};
    image<vec4f>                   preview = {};
    atomic<bool>*                  trace_stop = new atomic<bool>{false};
    atomic<int>*                   trace_sample = new atomic<int>{0};
    vector<future<void>>*           trace_futures = new vector<future<void>>{};
    concurrent_queue<image_region> trace_queue   = {};

    // view image
    vec2f            image_center = zero2f;
    float            image_scale  = 1;
    bool             zoom_to_fit  = true;
    bool             widgets_open = false;
    app_selection    selection    = {typeid(void), -1};
    vector<app_edit> update_list;
    bool             navigation_fps  = false;
    bool             quiet           = false;
    opengl_texture   display_texture = {};

    // app status
    bool load_done, load_running;
    string       status = "";
};

// Application state
struct app_state {
    vector<app_scene> scenes;
    int               selected = -1;
    deque<string>     errors;
};

void stop_rendering_async(app_scene& scn) {
    trace_image_async_stop(
        *scn.trace_futures, scn.trace_queue, scn.trace_options);
}

void start_rendering_async(app_scene& scn) {
    stop_rendering_async(scn);
    scn.status       = "rendering image";
    *scn.trace_stop   = false;
    *scn.trace_sample = 0;

    scn.image_size = get_camera_image_size(
        scn.scene.cameras[scn.trace_options.camera_id],
        scn.trace_options.image_size);
    scn.render  = {scn.image_size, zero4f};
    scn.display = {scn.image_size, zero4f};
    init_trace_state(scn.state, scn.image_size, scn.trace_options.random_seed);

    auto preview_options = scn.trace_options;
    preview_options.image_size /= scn.preview_ratio;
    preview_options.num_samples = 1;
    scn.preview = trace_image(scn.scene, scn.bvh, scn.lights, preview_options);
    auto display_preview = scn.preview;
    tonemap_image(display_preview, scn.preview, scn.tonemap_options);
    auto large_preview = image{scn.image_size, zero4f};
    for (auto j = 0; j < scn.image_size.y; j++) {
        for (auto i = 0; i < scn.image_size.x; i++) {
            auto pi = clamp(
                     i / scn.preview_ratio, 0, display_preview.size().x - 1),
                 pj = clamp(
                     j / scn.preview_ratio, 0, display_preview.size().y - 1);
            large_preview[{i, j}] = display_preview[{pi, pj}];
        }
    }
    scn.preview = large_preview;
    scn.trace_queue.push({{0, 0}, {0, 0}});

    scn.trace_options.cancel_flag = scn.trace_stop;
    trace_image_async_start(scn.render, scn.state, scn.scene, scn.bvh,
        scn.lights, *scn.trace_futures, *scn.trace_sample, scn.trace_queue,
        scn.trace_options);
}

bool load_scene_sync(app_scene& scn) {
    // scene loading
    scn.status = "loading scene";
    try {
        load_scene(scn.filename, scn.scene, scn.load_options);
    } catch (const std::exception& e) {
        print_fatal(e.what());
        return false;
    }

    // tesselate
    scn.status = "tesselating shapes";
    tesselate_subdivs(scn.scene);

    // add sky
    if (scn.add_skyenv) add_sky_environment(scn.scene);

    // build bvh
    scn.status = "computing bvh";
    build_scene_bvh(scn.scene, scn.bvh, scn.bvh_options);

    // init renderer
    scn.status = "initializing lights";
    init_trace_lights(scn.lights, scn.scene);

    // fix renderer type if no lights
    if (scn.lights.instances.empty() && scn.lights.environments.empty() &&
        is_trace_sampler_lit(scn.trace_options)) {
        print_info("no lights presents, switching to eyelight shader");
        scn.trace_options.sampler_type = trace_sampler_type::eyelight;
    }

    // set flags
    scn.load_done    = true;
    scn.load_running = false;
    scn.status       = "loading done";

    // start rendering
    start_rendering_async(scn);

    // done
    return false;
}

void load_scene_async(app_scene& scn) {
    if (scn.load_running) {
        print_info("error: already loading");
        return;
    }
    scn.load_done    = false;
    scn.load_running = true;
    scn.status       = "uninitialized";
    scn.scene        = {};
    scn.bvh          = {};
    scn.lights       = {};
    auto load_thread = thread([&scn]() { load_scene_sync(scn); });
    load_thread.detach();
}

void add_new_scene(app_state& app, const string& filename,
    const load_scene_options& load_options,
    const trace_image_options&   trace_options,
    const bvh_build_options&     bvh_options,
    const tonemap_image_options& tonemap_options,
    bool validate = false, bool add_skyenv = false) {
    auto& scn                           = app.scenes.emplace_back();
    scn.filename                        = filename;
    scn.imfilename                      = get_noextension(filename) + ".png";
    scn.load_options = load_options;
    scn.trace_options                   = trace_options;
    scn.trace_options.samples_per_batch = 1;
    scn.bvh_options                     = bvh_options;
    scn.tonemap_options                 = tonemap_options;
    scn.add_skyenv = add_skyenv;
    load_scene_async(scn);
}

void draw_opengl_widgets(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    begin_opengl_widgets_frame(win);
    if (begin_opengl_widgets_window(win, "yitrace")) {
        auto& scn = app.scenes[app.selected];
        if (begin_tabbar_opengl_widget(win, "tabs")) {
            if (begin_tabitem_opengl_widget(win, "trace")) {
                if (draw_button_opengl_widget(win, "load")) {
                    stop_rendering_async(scn);
                    load_scene_async(scn);
                }
                draw_label_opengl_widget(
                    win, "scene", get_filename(scn.filename));
                draw_label_opengl_widget(win, "filename", scn.filename);
                draw_label_opengl_widget(win, "status", scn.status);
                draw_label_opengl_widget(win, "image", "%d x %d @ %d",
                    scn.render.size().x, scn.render.size().y,
                    (int)*scn.trace_sample);
                auto cam_names = vector<string>();
                for (auto& camera : scn.scene.cameras)
                    cam_names.push_back(camera.uri);
                {
                    auto edited = false;
                    if (scn.load_done) {
                        if (draw_combobox_opengl_widget(win, "camera",
                                scn.trace_options.camera_id, cam_names)) {
                            edited = true;
                        }
                    }
                    if (draw_slider_opengl_widget(win, "width",
                            scn.trace_options.image_size.x, 0, 4096)) {
                        edited = true;
                    }
                    if (draw_slider_opengl_widget(win, "height",
                            scn.trace_options.image_size.y, 0, 4096)) {
                        edited = true;
                    }
                    if (draw_slider_opengl_widget(win, "nsamples",
                            scn.trace_options.num_samples, 16, 4096)) {
                        edited = true;
                    }
                    if (draw_combobox_opengl_widget(win, "tracer",
                            (int&)scn.trace_options.sampler_type,
                            trace_sampler_type_names)) {
                        edited = true;
                    }
                    if (draw_combobox_opengl_widget(win, "false color",
                            (int&)scn.trace_options.falsecolor_type,
                            trace_falsecolor_type_names)) {
                        edited = true;
                    }
                    if (draw_slider_opengl_widget(win, "nbounces",
                            scn.trace_options.max_bounces, 1, 10)) {
                        edited = true;
                    }
                    if (draw_checkbox_opengl_widget(win, "double sided",
                            scn.trace_options.double_sided)) {
                        edited = true;
                    }
                    if (draw_slider_opengl_widget(win, "seed",
                            (int&)scn.trace_options.random_seed, 0, 1000000)) {
                        edited = true;
                    }
                    if (draw_slider_opengl_widget(
                            win, "pratio", scn.preview_ratio, 1, 64)) {
                        edited = true;
                    }
                    if (edited) {
                        scn.update_list.push_back({typeid(trace_image_options),
                            -1, scn.trace_options, false});
                    }
                }
                draw_slider_opengl_widget(
                    win, "exposure", scn.tonemap_options.exposure, -5, 5);
                draw_checkbox_opengl_widget(
                    win, "filmic", scn.tonemap_options.filmic);
                continue_opengl_widget_line(win);
                draw_checkbox_opengl_widget(
                    win, "srgb", scn.tonemap_options.srgb);
                draw_slider_opengl_widget(
                    win, "zoom", scn.image_scale, 0.1, 10);
                draw_checkbox_opengl_widget(
                    win, "zoom to fit", scn.zoom_to_fit);
                continue_opengl_widget_line(win);
                draw_checkbox_opengl_widget(win, "fps", scn.navigation_fps);
                if (draw_button_opengl_widget(win, "print cams")) {
                    for (auto& camera : scn.scene.cameras) {
                        print_obj_camera(camera);
                    }
                }
                continue_opengl_widget_line(win);
                if (draw_button_opengl_widget(win, "print stats")) {
                    print_info("{}", format_scene_stats(scn.scene).c_str());
                    print_info("{}", print_scene_bvh_stats(scn.bvh).c_str());
                }
                auto mouse_pos = get_opengl_mouse_pos(win);
                auto ij        = get_image_coords(mouse_pos, scn.image_center,
                    scn.image_scale, scn.render.size());
                draw_dragger_opengl_widget(win, "mouse", ij);
                if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
                    ij.y < scn.render.size().y) {
                    draw_coloredit_opengl_widget(
                        win, "pixel", scn.render[{ij.x, ij.y}]);
                } else {
                    auto zero4f_ = zero4f;
                    draw_coloredit_opengl_widget(win, "pixel", zero4f_);
                }
                end_tabitem_opengl_widget(win);
            }
            if (scn.load_done && begin_tabitem_opengl_widget(win, "navigate")) {
                draw_opengl_widgets_scene_tree(
                    win, "", scn.scene, scn.selection, scn.update_list, 200);
                end_tabitem_opengl_widget(win);
            }
            if (scn.load_done && begin_tabitem_opengl_widget(win, "inspec")) {
                draw_opengl_widgets_scene_inspector(
                    win, "", scn.scene, scn.selection, scn.update_list, 200);
                end_tabitem_opengl_widget(win);
            }
            end_tabbar_opengl_widget(win);
        }
    }
    end_opengl_widgets_frame(win);
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    set_opengl_viewport(get_opengl_framebuffer_size(win));
    clear_opengl_lframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    for (auto& scn : app.scenes) {
        if (!scn.load_done) break;
        update_image_view(scn.image_center, scn.image_scale, scn.render.size(),
            win_size, scn.zoom_to_fit);
        if (!scn.display_texture ||
            scn.display_texture.size != scn.image_size) {
            if (scn.display_texture) delete_opengl_texture(scn.display_texture);
            if (scn.image_size != zero2i) {
                init_opengl_texture(scn.display_texture, scn.image_size, false,
                    false, false, false);
            }
        } else {
            auto region = image_region{};
            auto size   = 0;
            while (scn.trace_queue.try_pop(region)) {
                if (region.size() == zero2i) {
                    update_opengl_texture(
                        scn.display_texture, scn.preview, false);
                    break;
                } else {
                    tonemap_image_region(scn.display, region, scn.render,
                        scn.tonemap_options);
                    update_opengl_texture_region(
                        scn.display_texture, scn.display, region, false);
                    size += region.size().x * region.size().y;
                    if (size >= scn.render.size().x * scn.render.size().y)
                        break;
                }
            }
        }
    }
    if (!app.scenes.empty() && app.selected >= 0) {
        auto& scn = app.scenes[app.selected];
        set_opengl_blending(true);
        draw_opengl_image_background(scn.display_texture, win_size.x,
            win_size.y, scn.image_center, scn.image_scale);
        draw_opengl_image(scn.display_texture, win_size.x, win_size.y,
            scn.image_center, scn.image_scale);
        set_opengl_blending(false);
    }
    draw_opengl_widgets(win);
    swap_opengl_buffers(win);
}

bool update(app_state& app) {
    auto updated = false;
    for (auto& scn : app.scenes) {
        // exit if no updated
        if (!scn.load_done || scn.update_list.empty()) return false;

        // stop renderer
        stop_rendering_async(scn);

        // update data
        auto updated_instances = vector<int>{}, updated_shapes = vector<int>{};
        auto updated_lights = false;
        for (auto& [type, index, data, reload] : scn.update_list) {
            if (type == typeid(yocto_camera)) {
                scn.scene.cameras[index] = any_cast<yocto_camera>(data);
            } else if (type == typeid(yocto_texture)) {
                scn.scene.textures[index] = any_cast<yocto_texture>(data);
                if (reload) {
                    auto& texture = scn.scene.textures[index];
                    load_image(get_dirname(scn.filename) + texture.uri,
                        texture.hdr_image, texture.ldr_image);
                }
            } else if (type == typeid(yocto_voltexture)) {
                scn.scene.voltextures[index] = any_cast<yocto_voltexture>(data);
                if (reload) {
                    auto& texture = scn.scene.voltextures[index];
                    load_volume(get_dirname(scn.filename) + texture.uri,
                        texture.volume_data);
                }
            } else if (type == typeid(yocto_shape)) {
                scn.scene.shapes[index] = any_cast<yocto_shape>(data);
                if (reload) {
                    auto& shape = scn.scene.shapes[index];
                    load_shape(get_dirname(scn.filename) + shape.uri,
                        shape.points, shape.lines, shape.triangles, shape.quads,
                        shape.quads_positions, shape.quads_normals,
                        shape.quads_texturecoords, shape.positions,
                        shape.normals, shape.texturecoords, shape.colors,
                        shape.radius, false);
                }
                updated_shapes.push_back(index);
            } else if (type == typeid(yocto_subdiv)) {
                // TODO: this needs more fixing?
                scn.scene.subdivs[index] = any_cast<yocto_subdiv>(data);
                if (reload) {
                    auto& subdiv = scn.scene.subdivs[index];
                    load_shape(get_dirname(scn.filename) + subdiv.uri,
                        subdiv.points, subdiv.lines, subdiv.triangles,
                        subdiv.quads, subdiv.quads_positions,
                        subdiv.quads_normals, subdiv.quads_texturecoords,
                        subdiv.positions, subdiv.normals, subdiv.texturecoords,
                        subdiv.colors, subdiv.radius,
                        subdiv.preserve_facevarying);
                }
                tesselate_subdiv(scn.scene, scn.scene.subdivs[index]);
                updated_shapes.push_back(
                    scn.scene.subdivs[index].tesselated_shape);
            } else if (type == typeid(yocto_material)) {
                auto old_emission = scn.scene.materials[index].emission;
                scn.scene.materials[index] = any_cast<yocto_material>(data);
                if (old_emission != scn.scene.materials[index].emission) {
                    updated_lights = true;
                }
            } else if (type == typeid(yocto_instance)) {
                scn.scene.instances[index] = any_cast<yocto_instance>(data);
                updated_instances.push_back(index);
            } else if (type == typeid(yocto_environment)) {
                auto old_emission = scn.scene.materials[index].emission;
                scn.scene.environments[index] = any_cast<yocto_environment>(
                    data);
                if (old_emission != scn.scene.materials[index].emission) {
                    updated_lights = true;
                }
            } else if (type == typeid(trace_image_options)) {
                scn.trace_options = any_cast<trace_image_options>(data);
            } else {
                throw runtime_error("unsupported type "s + type.name());
            }
        }
        // update bvh
        if (!updated_instances.empty() || !updated_shapes.empty()) {
            refit_scene_bvh(scn.scene, scn.bvh, updated_instances,
                updated_shapes, scn.bvh_options);
        }
        // update lights
        if (updated_lights) {
            init_trace_lights(scn.lights, scn.scene);
        }

        // clear
        scn.update_list.clear();

        // start rendering
        start_rendering_async(scn);
    }

    // updated
    return true;
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    for (auto& path : paths) add_new_scene(app, path, {}, {}, {}, {});
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280, 720}, "yitrace", &app, draw);
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
        if (app.selected >= 0 && app.scenes[app.selected].load_done &&
            (mouse_left || mouse_right) && !alt_down && !widgets_active) {
            auto& scn    = app.scenes[app.selected];
            auto  camera = scn.scene.cameras.at(scn.trace_options.camera_id);
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
            scn.update_list.push_back({typeid(yocto_camera),
                scn.trace_options.camera_id, camera, false});
        }

        // selection
        if (app.selected >= 0 && app.scenes[app.selected].load_done &&
            (mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto& scn = app.scenes[app.selected];
            auto  ij  = get_image_coords(mouse_pos, scn.image_center,
                scn.image_scale, scn.render.size());
            if (ij.x < 0 || ij.x >= scn.render.size().x || ij.y < 0 ||
                ij.y >= scn.render.size().y) {
                auto& camera = scn.scene.cameras.at(
                    scn.trace_options.camera_id);
                auto ray = evaluate_camera_ray(
                    camera, ij, scn.render.size(), {0.5f, 0.5f}, zero2f);
                if (auto isec = bvh_intersection{};
                    intersect_scene_bvh(scn.scene, scn.bvh, ray, isec)) {
                    scn.selection = {typeid(yocto_instance), isec.instance_id};
                }
            }
        }

        // update
        update(app);

        // draw
        draw(win);

        // event hadling
        process_opengl_events(win);
    }

    // stop
    for(auto& scn : app.scenes) stop_rendering_async(scn);

    // clear
    delete_opengl_window(win);
}

int main(int argc, char* argv[]) {
    // application
    app_state app{};
    auto      load_options          = load_scene_options{};
    auto      bvh_options           = bvh_build_options{};
    auto      trace_options         = trace_image_options{};
    auto      tonemap_options       = tonemap_image_options{};
    auto      no_parallel           = false;
    auto      add_skyenv            = false;
    auto      validate              = false;
    auto      filenames             = vector<string>{};
    trace_options.samples_per_batch = 1;

    // names for enums
    auto trace_sampler_type_namemap = std::map<string, trace_sampler_type>{};
    for (auto type = 0; type < trace_sampler_type_names.size(); type++) {
        trace_sampler_type_namemap[trace_sampler_type_names[type]] =
            (trace_sampler_type)type;
    }
    auto trace_falsecolor_type_namemap =
        std::map<string, trace_falsecolor_type>{};
    for (auto type = 0; type < trace_falsecolor_type_names.size(); type++) {
        trace_falsecolor_type_namemap[trace_falsecolor_type_names[type]] =
            (trace_falsecolor_type)type;
    }

    // parse command line
    auto parser = CLI::App{"progressive path tracing"};
    parser.add_option("--camera", trace_options.camera_id, "Camera index.");
    parser.add_option("--hres,-R", trace_options.image_size.x,
        "Image horizontal resolution.");
    parser.add_option(
        "--vres,-r", trace_options.image_size.y, "Image vertical resolution.");
    parser.add_option(
        "--nsamples,-s", trace_options.num_samples, "Number of samples.");
    parser
        .add_option("--tracer,-t", trace_options.sampler_type, "Tracer type.")
        ->transform(CLI::IsMember(trace_sampler_type_namemap));
    parser
        .add_option("--falsecolor,-F", trace_options.falsecolor_type,
            "Tracer false color type.")
        ->transform(CLI::IsMember(trace_falsecolor_type_namemap));
    parser.add_option(
        "--nbounces", trace_options.max_bounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", trace_options.pixel_clamp, "Final pixel clamping.");
    parser.add_option("--seed", trace_options.random_seed,
        "Seed for the random number generators.");
    parser.add_flag("--env-hidden,!--no-env-hidden",
        trace_options.environments_hidden,
        "Environments are hidden in renderer");
    parser.add_flag("--parallel,!--no-parallel", no_parallel,
        "Disable parallel execution.");
    parser.add_flag("--bvh-high-quality,!--no-bvh-high-quality",
        bvh_options.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
    parser.add_flag("--bvh-embree,!--no-bvh-embree", bvh_options.use_embree,
        "Use Embree ratracer");
    parser.add_flag("--bvh-embree-flatten,!--no-bvh-embree-flatten",
        bvh_options.embree_flatten, "Flatten embree scene");
    parser.add_flag("--bvh-embree-compact,!--no-bvh-embree-compact",
        bvh_options.embree_compact, "Embree runs in compact memory");
#endif
    parser.add_flag("--double-sided,!--no-double-sided",
        trace_options.double_sided, "Double-sided rendering.");
    parser.add_flag(
        "--add-skyenv,!--no-add-skyenv", add_skyenv, "Add sky envmap");
    parser.add_option("scenes", filenames, "Scene filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix parallel code
    if (no_parallel) {
        bvh_options.run_serially   = true;
        load_options.run_serially  = true;
        trace_options.run_serially = true;
    }

    // loading images
    for (auto filename : filenames)
        add_new_scene(app, filename, load_options, trace_options, bvh_options,
            tonemap_options, validate, add_skyenv);
    app.selected = 0;

    // run interactive
    run_ui(app);

    // done
    return 0;
}
