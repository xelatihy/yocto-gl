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
enum struct app_task_type {
    none,
    load,
    bvh,
    lights,
    render,
    param_edit,
    scene_edit,
    save_image,
    save_scene
};

struct app_task {
    app_task_type                  type;
    future<void>                   result;
    atomic<bool>                   stop;
    atomic<int>                    current;
    concurrent_queue<image_region> queue;
    app_edit                       edit;

    app_task(app_task_type type, const app_edit& edit = {})
        : type{type}, result{}, stop{false}, current{-1}, queue{}, edit{edit} {}
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
    string filename  = "scene.yaml";
    string imagename = "out.png";
    string outname   = "out.yaml";
    string name      = "";

    // options
    load_scene_options    load_options    = {};
    save_scene_options    save_options    = {};
    bvh_build_options     bvh_options     = {};
    trace_image_options   trace_options   = {};
    tonemap_image_options tonemap_options = {};
    int                   preview_ratio   = 8;
    vec2i                 image_size      = {0, 0};

    // scene
    yocto_scene scene      = {};
    bvh_scene   bvh        = {};
    bool        add_skyenv = false;

    // rendering state
    trace_lights lights        = {};
    trace_state  state         = {};
    image<vec4f> render        = {};
    image<vec4f> display       = {};
    image<vec4f> preview       = {};
    int          render_sample = 0;

    // view scene
    vec2f          image_center   = zero2f;
    float          image_scale    = 1;
    bool           zoom_to_fit    = true;
    bool           navigation_fps = false;
    opengl_texture gl_txt         = {};

    // tasks
    bool load_done = false, bvh_done = false, lights_done = false,
         render_done = false;
    deque<app_task> task_queue;
    app_selection   selection = {typeid(void), -1};
};

// Application state
struct app_state {
    deque<app_scene> scenes;
    int              selected = -1;
    deque<string>    errors;
};

void update_app_render(const string& filename, image<vec4f>& render,
    image<vec4f>& display, image<vec4f>& preview, trace_state& state,
    const yocto_scene& scene, const trace_lights& lights, const bvh_scene& bvh,
    const trace_image_options&   trace_options,
    const tonemap_image_options& tonemap_options, int preview_ratio,
    atomic<bool>& stop, atomic<int>& current_sample,
    concurrent_queue<image_region>& queue) {
    auto preview_options = trace_options;
    preview_options.image_size /= preview_ratio;
    preview_options.num_samples = 1;
    auto small_preview   = trace_image(scene, bvh, lights, preview_options);
    auto display_preview = small_preview;
    tonemap_image(display_preview, small_preview, tonemap_options);
    for (auto j = 0; j < preview.size().y; j++) {
        for (auto i = 0; i < preview.size().x; i++) {
            auto pi = clamp(i / preview_ratio, 0, display_preview.size().x - 1),
                 pj = clamp(j / preview_ratio, 0, display_preview.size().y - 1);
            preview[{i, j}] = display_preview[{pi, pj}];
        }
    }
    queue.push({{0, 0}, {0, 0}});
    current_sample = 0;

    auto& camera     = scene.cameras.at(trace_options.camera_id);
    auto  image_size = get_camera_image_size(camera, trace_options.image_size);
    state            = trace_state{};
    init_trace_state(state, image_size, trace_options.random_seed);
    auto regions = vector<image_region>{};
    make_image_regions(regions, render.size(), trace_options.region_size, true);

    for (auto sample = 0; sample < trace_options.num_samples;
         sample += trace_options.samples_per_batch) {
        if (stop) return;
        current_sample   = sample;
        auto num_samples = min(trace_options.samples_per_batch,
            trace_options.num_samples - current_sample);
        parallel_foreach(
            regions,
            [num_samples, &trace_options, &tonemap_options, &render, &display,
                &scene, &lights, &bvh, &state,
                &queue](const image_region& region) {
                trace_image_region(render, state, scene, bvh, lights, region,
                    num_samples, trace_options);
                tonemap_image_region(display, region, render, tonemap_options);
                queue.push(region);
            },
            &stop);
    }
    current_sample = trace_options.num_samples;
}

void add_new_scene(app_state& app, const string& filename,
    const load_scene_options&    load_options,
    const trace_image_options&   trace_options,
    const bvh_build_options&     bvh_options,
    const tonemap_image_options& tonemap_options, bool validate = false,
    bool add_skyenv = false) {
    auto& scn         = app.scenes.emplace_back();
    scn.filename      = filename;
    scn.imagename     = get_noextension(filename) + ".png";
    scn.outname       = get_noextension(filename) + ".edited.yaml";
    scn.name          = get_filename(scn.filename);
    scn.load_options  = load_options;
    scn.trace_options = trace_options;
    scn.trace_options.samples_per_batch = 1;
    scn.bvh_options                     = bvh_options;
    scn.tonemap_options                 = tonemap_options;
    scn.add_skyenv                      = add_skyenv;
    scn.task_queue.emplace_back(app_task_type::load);
    app.selected = (int)app.scenes.size() - 1;
}

void apply_scene_edit(const string& filename, yocto_scene& scene,
    bvh_scene& bvh, trace_lights& lights, const app_edit& edit,
    const bvh_build_options& bvh_options) {
    // update data
    auto updated_instances = vector<int>{}, updated_shapes = vector<int>{};
    auto updated_lights = false;

    auto& [type, index, data, reload] = edit;

    if (type == typeid(yocto_camera)) {
        scene.cameras[index] = any_cast<yocto_camera>(data);
    } else if (type == typeid(yocto_texture)) {
        scene.textures[index] = any_cast<yocto_texture>(data);
        if (reload) {
            auto& texture = scene.textures[index];
            load_image(get_dirname(filename) + texture.uri, texture.hdr_image,
                texture.ldr_image);
        }
    } else if (type == typeid(yocto_voltexture)) {
        scene.voltextures[index] = any_cast<yocto_voltexture>(data);
        if (reload) {
            auto& texture = scene.voltextures[index];
            load_volume(
                get_dirname(filename) + texture.uri, texture.volume_data);
        }
    } else if (type == typeid(yocto_shape)) {
        scene.shapes[index] = any_cast<yocto_shape>(data);
        if (reload) {
            auto& shape = scene.shapes[index];
            load_shape(get_dirname(filename) + shape.uri, shape.points,
                shape.lines, shape.triangles, shape.quads,
                shape.quads_positions, shape.quads_normals,
                shape.quads_texturecoords, shape.positions, shape.normals,
                shape.texturecoords, shape.colors, shape.radius, false);
        }
        updated_shapes.push_back(index);
    } else if (type == typeid(yocto_subdiv)) {
        // TODO: this needs more fixing?
        scene.subdivs[index] = any_cast<yocto_subdiv>(data);
        if (reload) {
            auto& subdiv = scene.subdivs[index];
            load_shape(get_dirname(filename) + subdiv.uri, subdiv.points,
                subdiv.lines, subdiv.triangles, subdiv.quads,
                subdiv.quads_positions, subdiv.quads_normals,
                subdiv.quads_texturecoords, subdiv.positions, subdiv.normals,
                subdiv.texturecoords, subdiv.colors, subdiv.radius,
                subdiv.preserve_facevarying);
        }
        tesselate_subdiv(scene, scene.subdivs[index]);
        updated_shapes.push_back(scene.subdivs[index].tesselated_shape);
    } else if (type == typeid(yocto_material)) {
        auto old_emission      = scene.materials[index].emission;
        scene.materials[index] = any_cast<yocto_material>(data);
        if (old_emission != scene.materials[index].emission) {
            updated_lights = true;
        }
    } else if (type == typeid(yocto_instance)) {
        scene.instances[index] = any_cast<yocto_instance>(data);
        updated_instances.push_back(index);
    } else if (type == typeid(yocto_environment)) {
        auto old_emission         = scene.materials[index].emission;
        scene.environments[index] = any_cast<yocto_environment>(data);
        if (old_emission != scene.materials[index].emission) {
            updated_lights = true;
        }
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }

    // update bvh
    if (!updated_instances.empty() || !updated_shapes.empty()) {
        refit_scene_bvh(
            scene, bvh, updated_instances, updated_shapes, bvh_options);
    }
    // update lights
    if (updated_lights) init_trace_lights(lights, scene);
}

void apply_param_edit(const string& filename,
    trace_image_options& trace_options, tonemap_image_options& tonemap_options,
    const app_edit& edit) {
    // update data
    auto& [type, index, data, reload] = edit;
    if (type == typeid(trace_image_options)) {
        trace_options = any_cast<trace_image_options>(data);
    } else if (type == typeid(tonemap_image_options)) {
        tonemap_options = any_cast<tonemap_image_options>(data);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }
}

void draw_opengl_widgets(const opengl_window& win) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    if (begin_opengl_widgets_window(win, "yitrace")) {
        auto& scn = app.scenes[app.selected];
        if (begin_tabbar_opengl_widget(win, "tabs")) {
            if (begin_tabitem_opengl_widget(win, "trace")) {
                if (draw_button_opengl_widget(win, "load")) {
                    // TODO
                }
                draw_label_opengl_widget(
                    win, "scene", get_filename(scn.filename));
                draw_label_opengl_widget(win, "filename", scn.filename);
                draw_label_opengl_widget(win, "image", "%d x %d @ %d",
                    scn.render.size().x, scn.render.size().y,
                    scn.render_sample);
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
                        scn.task_queue.emplace_back(app_task_type::param_edit,
                            app_edit{typeid(trace_image_options), -1,
                                scn.trace_options, false});
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
                    win, "", scn.scene, scn.selection, 200);
                end_tabitem_opengl_widget(win);
            }
            if (scn.load_done && begin_tabitem_opengl_widget(win, "inspec")) {
                auto edit = app_edit{};
                if (draw_opengl_widgets_scene_inspector(
                        win, "", scn.scene, scn.selection, edit, 200)) {
                    scn.task_queue.emplace_back(
                        app_task_type::scene_edit, edit);
                }
                end_tabitem_opengl_widget(win);
            }
            end_tabbar_opengl_widget(win);
        }
    }
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    set_opengl_viewport(get_opengl_framebuffer_size(win));
    clear_opengl_lframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
    if (!app.scenes.empty() && app.selected >= 0) {
        auto& scn = app.scenes.at(app.selected);
        if (scn.load_done && scn.gl_txt) {
            update_image_view(scn.image_center, scn.image_scale,
                scn.display.size(), win_size, scn.zoom_to_fit);
            draw_opengl_image_background(scn.gl_txt, win_size.x, win_size.y,
                scn.image_center, scn.image_scale);
            set_opengl_blending(true);
            draw_opengl_image(scn.gl_txt, win_size.x, win_size.y,
                scn.image_center, scn.image_scale);
            set_opengl_blending(false);
        }
    }
    begin_opengl_widgets_frame(win);
    draw_opengl_widgets(win);
    end_opengl_widgets_frame(win);
    swap_opengl_buffers(win);
}

void update(app_state& app) {
    // consume partial results
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (task.type != app_task_type::render || task.queue.empty()) continue;
        auto region  = image_region{};
        auto updated = false;
        while (scn.task_queue.front().queue.try_pop(region)) {
            if (region.size() == zero2i) {
                update_opengl_texture_region(scn.gl_txt, scn.preview,
                    {zero2i, scn.preview.size()}, false);
            } else {
                update_opengl_texture_region(
                    scn.gl_txt, scn.display, region, false);
            }
            updated = true;
        }
        if (updated) {
            scn.render_sample = max(scn.render_sample, (int)task.current);
            scn.name          = format(
                "{} [{}x{}@{}]", get_filename(scn.filename), scn.render.size().x,
                                     scn.render.size().y, scn.render_sample);
        }
    }
    // remove unneeded tasks
    for (auto& scn : app.scenes) {
        while (scn.task_queue.size() > 1 &&
               scn.task_queue.at(0).type == app_task_type::render &&
               (scn.task_queue.at(1).type == app_task_type::render ||
                   scn.task_queue.at(1).type == app_task_type::scene_edit ||
                   scn.task_queue.at(1).type == app_task_type::param_edit)) {
            log_info("cancel rendering {}", scn.filename);
            auto& task = scn.task_queue.front();
            task.stop  = true;
            if (task.result.valid()) {
                try {
                    task.result.get();
                } catch (...) {
                }
            }
            scn.task_queue.pop_front();
        }
    }
    // schedule tasks not running
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (task.result.valid()) continue;
        task.stop = false;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::load: {
                log_info("start loading {}", scn.filename);
                scn.load_done   = false;
                scn.bvh_done    = false;
                scn.lights_done = false;
                task.result     = async([&scn]() {
                    load_scene(scn.filename, scn.scene, scn.load_options);
                    tesselate_subdivs(scn.scene);
                    if (scn.add_skyenv) add_sky_environment(scn.scene);
                });
            } break;
            case app_task_type::bvh: {
                log_info("start building bvh {}", scn.filename);
                scn.bvh_done = false;
                task.result  = async([&scn]() {
                    build_scene_bvh(scn.scene, scn.bvh, scn.bvh_options);
                });
            } break;
            case app_task_type::lights: {
                log_info("start building lights {}", scn.filename);
                scn.lights_done = false;
                task.result     = async(
                    [&scn]() { init_trace_lights(scn.lights, scn.scene); });
            } break;
            case app_task_type::save_image: {
                log_info("start saving {}", scn.imagename);
                task.result = async([&scn]() {
                    save_tonemapped_image(
                        scn.imagename, scn.render, scn.tonemap_options);
                });
            } break;
            case app_task_type::save_scene: {
                log_info("start saving {}", scn.outname);
                task.result = async([&scn]() {
                    save_scene(scn.outname, scn.scene, scn.save_options);
                });
            } break;
            case app_task_type::render: {
                log_info("start rendering {}", scn.filename);
                scn.render_done = false;
                scn.image_size  = get_camera_image_size(
                    scn.scene.cameras[scn.trace_options.camera_id],
                    scn.trace_options.image_size);
                scn.render.assign(scn.image_size, zero4f);
                scn.display.assign(scn.image_size, zero4f);
                scn.preview.assign(scn.image_size, zero4f);
                if (scn.lights.instances.empty() &&
                    scn.lights.environments.empty() &&
                    is_trace_sampler_lit(scn.trace_options)) {
                    log_info(
                        "no lights presents, switching to eyelight shader");
                    scn.trace_options.sampler_type =
                        trace_sampler_type::eyelight;
                }
                scn.render_sample = scn.trace_options.num_samples;
                scn.name          = format("{} [{}x{}@{}]",
                    get_filename(scn.filename), scn.render.size().x,
                        scn.render.size().y, scn.render_sample);
                task.result       = async([&scn, &task]() {
                    update_app_render(scn.filename, scn.render, scn.display,
                        scn.preview, scn.state, scn.scene, scn.lights, scn.bvh,
                        scn.trace_options, scn.tonemap_options,
                        scn.preview_ratio, task.stop, task.current, task.queue);
                });
                init_opengl_texture(
                    scn.gl_txt, scn.display, false, false, false);
            } break;
            case app_task_type::scene_edit: {
                log_info("start editing {}", scn.filename);
                scn.render_done = false;
                task.result     = async([&scn, &task]() {
                    apply_scene_edit(scn.filename, scn.scene, scn.bvh,
                        scn.lights, task.edit, scn.bvh_options);
                });
            } break;
            case app_task_type::param_edit: {
                log_info("start editing params for {}", scn.filename);
                scn.render_done = false;
                task.result     = async([&scn, &task]() {
                    apply_param_edit(scn.filename, scn.trace_options,
                        scn.tonemap_options, task.edit);
                });
            } break;
        }
    }
    // grab result of finished tasks
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (!task.result.valid() || !task.queue.empty() ||
            task.result.wait_for(std::chrono::nanoseconds(10)) !=
                std::future_status::ready)
            continue;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::load: {
                try {
                    task.result.get();
                    scn.load_done = true;
                    scn.name      = format("{}", get_filename(scn.filename));
                    log_info("done loading {}", scn.filename);
                    scn.task_queue.emplace_back(app_task_type::bvh);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot load " + scn.filename);
                }
            } break;
            case app_task_type::bvh: {
                try {
                    task.result.get();
                    scn.bvh_done = true;
                    scn.name     = format("{}", get_filename(scn.filename));
                    log_info("done building bvh {}", scn.filename);
                    scn.task_queue.emplace_back(app_task_type::lights);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot build bvh " + scn.filename);
                }
            } break;
            case app_task_type::lights: {
                try {
                    task.result.get();
                    scn.lights_done = true;
                    scn.name        = format("{}", get_filename(scn.filename));
                    log_info("done building lights {}", scn.filename);
                    scn.task_queue.emplace_back(app_task_type::render);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot build lights " + scn.filename);
                }
            } break;
            case app_task_type::save_image: {
                try {
                    task.result.get();
                    log_info("done saving {}", scn.imagename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot save " + scn.imagename);
                }
            } break;
            case app_task_type::save_scene: {
                try {
                    task.result.get();
                    log_info("done saving {}", scn.outname);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot save " + scn.outname);
                }
            } break;
            case app_task_type::render: {
                try {
                    task.result.get();
                    scn.render_done = true;
                    log_info("done rendering {}", scn.filename);
                    scn.render_sample = scn.trace_options.num_samples;
                    scn.name          = format("{} [{}x{}@{}]",
                        get_filename(scn.filename), scn.render.size().x,
                            scn.render.size().y, scn.render_sample);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot render " + scn.filename);
                }
            } break;
            case app_task_type::scene_edit: {
                try {
                    task.result.get();
                    log_info("done editing {}", scn.filename);
                    scn.task_queue.emplace_back(app_task_type::render);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot edit " + scn.filename);
                }
            } break;
            case app_task_type::param_edit: {
                try {
                    task.result.get();
                    log_info("done editing params for {}", scn.filename);
                    scn.task_queue.emplace_back(app_task_type::render);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot edit " + scn.filename);
                }
            } break;
        }
        scn.task_queue.pop_front();
    }

    // for (auto& scn : app.scenes) {
    //     // exit if no updated
    //     if (!scn.load_done || scn.update_list.empty()) return false;

    //     // stop renderer
    //     stop_rendering_async(scn);

    //     // update data
    //     auto updated_instances = vector<int>{}, updated_shapes =
    //     vector<int>{}; auto updated_lights = false; for (auto& [type, index,
    //     data, reload] : scn.update_list) {
    //         if (type == typeid(yocto_camera)) {
    //             scn.scene.cameras[index] = any_cast<yocto_camera>(data);
    //         } else if (type == typeid(yocto_texture)) {
    //             scn.scene.textures[index] = any_cast<yocto_texture>(data);
    //             if (reload) {
    //                 auto& texture = scn.scene.textures[index];
    //                 load_image(get_dirname(scn.filename) + texture.uri,
    //                     texture.hdr_image, texture.ldr_image);
    //             }
    //         } else if (type == typeid(yocto_voltexture)) {
    //             scn.scene.voltextures[index] =
    //             any_cast<yocto_voltexture>(data); if (reload) {
    //                 auto& texture = scn.scene.voltextures[index];
    //                 load_volume(get_dirname(scn.filename) + texture.uri,
    //                     texture.volume_data);
    //             }
    //         } else if (type == typeid(yocto_shape)) {
    //             scn.scene.shapes[index] = any_cast<yocto_shape>(data);
    //             if (reload) {
    //                 auto& shape = scn.scene.shapes[index];
    //                 load_shape(get_dirname(scn.filename) + shape.uri,
    //                     shape.points, shape.lines, shape.triangles,
    //                     shape.quads, shape.quads_positions,
    //                     shape.quads_normals, shape.quads_texturecoords,
    //                     shape.positions, shape.normals, shape.texturecoords,
    //                     shape.colors, shape.radius, false);
    //             }
    //             updated_shapes.push_back(index);
    //         } else if (type == typeid(yocto_subdiv)) {
    //             // TODO: this needs more fixing?
    //             scn.scene.subdivs[index] = any_cast<yocto_subdiv>(data);
    //             if (reload) {
    //                 auto& subdiv = scn.scene.subdivs[index];
    //                 load_shape(get_dirname(scn.filename) + subdiv.uri,
    //                     subdiv.points, subdiv.lines, subdiv.triangles,
    //                     subdiv.quads, subdiv.quads_positions,
    //                     subdiv.quads_normals, subdiv.quads_texturecoords,
    //                     subdiv.positions, subdiv.normals,
    //                     subdiv.texturecoords, subdiv.colors, subdiv.radius,
    //                     subdiv.preserve_facevarying);
    //             }
    //             tesselate_subdiv(scn.scene, scn.scene.subdivs[index]);
    //             updated_shapes.push_back(
    //                 scn.scene.subdivs[index].tesselated_shape);
    //         } else if (type == typeid(yocto_material)) {
    //             auto old_emission = scn.scene.materials[index].emission;
    //             scn.scene.materials[index] = any_cast<yocto_material>(data);
    //             if (old_emission != scn.scene.materials[index].emission) {
    //                 updated_lights = true;
    //             }
    //         } else if (type == typeid(yocto_instance)) {
    //             scn.scene.instances[index] = any_cast<yocto_instance>(data);
    //             updated_instances.push_back(index);
    //         } else if (type == typeid(yocto_environment)) {
    //             auto old_emission = scn.scene.materials[index].emission;
    //             scn.scene.environments[index] = any_cast<yocto_environment>(
    //                 data);
    //             if (old_emission != scn.scene.materials[index].emission) {
    //                 updated_lights = true;
    //             }
    //         } else if (type == typeid(trace_image_options)) {
    //             scn.trace_options = any_cast<trace_image_options>(data);
    //         } else {
    //             throw runtime_error("unsupported type "s + type.name());
    //         }
    //     }
    //     // update bvh
    //     if (!updated_instances.empty() || !updated_shapes.empty()) {
    //         refit_scene_bvh(scn.scene, scn.bvh, updated_instances,
    //             updated_shapes, scn.bvh_options);
    //     }
    //     // update lights
    //     if (updated_lights) {
    //         init_trace_lights(scn.lights, scn.scene);
    //     }

    //     // clear
    //     scn.update_list.clear();

    //     // start rendering
    //     start_rendering_async(scn);
    // }
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

    // setup logging
    set_log_callback(
        [&win](const string& msg) { add_log_opengl_widget(win, msg.c_str()); });

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
            scn.task_queue.emplace_back(app_task_type::scene_edit,
                app_edit{typeid(yocto_camera), scn.trace_options.camera_id,
                    camera, false});
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
