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
    load_scene,
    load_element,
    build_bvh,
    refit_bvh,
    init_lights,
    render_image,
    apply_edit,
    save_image,
    save_scene,
    close_scene
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
    sceneio_params sceneio_prms  = {};
    bvh_params     bvh_prms      = {};
    trace_params   trace_prms    = {};
    tonemap_params tonemap_prms  = {};
    int            preview_ratio = 8;
    vec2i          image_size    = {0, 0};

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
    // data
    deque<app_scene> scenes;
    int              selected = -1;
    deque<string>    errors;

    // default options
    sceneio_params sceneio_prms = {};
    bvh_params     bvh_prms     = {};
    trace_params   trace_prms   = {};
    tonemap_params tonemap_prms = {};
    bool           add_skyenv   = false;
};

void update_app_render(const string& filename, image<vec4f>& render,
    image<vec4f>& display, image<vec4f>& preview, trace_state& state,
    const yocto_scene& scene, const trace_lights& lights, const bvh_scene& bvh,
    const trace_params& trace_prms, const tonemap_params& tonemap_prms,
    int preview_ratio, atomic<bool>& stop, atomic<int>& current_sample,
    concurrent_queue<image_region>& queue) {
    auto preview_options = trace_prms;
    preview_options.image_size /= preview_ratio;
    preview_options.num_samples = 1;
    auto small_preview   = trace_image(scene, bvh, lights, preview_options);
    auto display_preview = small_preview;
    tonemap(display_preview, small_preview, tonemap_prms);
    for (auto j = 0; j < preview.size().y; j++) {
        for (auto i = 0; i < preview.size().x; i++) {
            auto pi = clamp(i / preview_ratio, 0, display_preview.size().x - 1),
                 pj = clamp(j / preview_ratio, 0, display_preview.size().y - 1);
            preview[{i, j}] = display_preview[{pi, pj}];
        }
    }
    queue.push({{0, 0}, {0, 0}});
    current_sample = 0;

    auto& camera     = scene.cameras.at(trace_prms.camera_id);
    auto  image_size = camera_image_size(camera, trace_prms.image_size);
    state            = trace_state{};
    init_trace_state(state, image_size, trace_prms.random_seed);
    auto regions = vector<image_region>{};
    make_regions(regions, render.size(), trace_prms.region_size, true);

    for (auto sample = 0; sample < trace_prms.num_samples;
         sample += trace_prms.samples_per_batch) {
        if (stop) return;
        current_sample   = sample;
        auto num_samples = min(trace_prms.samples_per_batch,
            trace_prms.num_samples - current_sample);
        parallel_foreach(
            regions,
            [num_samples, &trace_prms, &tonemap_prms, &render, &display, &scene,
                &lights, &bvh, &state, &queue](const image_region& region) {
                trace_region(render, state, scene, bvh, lights, region,
                    num_samples, trace_prms);
                tonemap(display, render, region, tonemap_prms);
                queue.push(region);
            },
            &stop);
    }
    current_sample = trace_prms.num_samples;
}

void add_new_scene(app_state& app, const string& filename) {
    auto& scn        = app.scenes.emplace_back();
    scn.filename     = filename;
    scn.imagename    = get_noextension(filename) + ".png";
    scn.outname      = get_noextension(filename) + ".edited.yaml";
    scn.name         = get_filename(scn.filename);
    scn.sceneio_prms = app.sceneio_prms;
    scn.trace_prms   = app.trace_prms;
    scn.bvh_prms     = app.bvh_prms;
    scn.tonemap_prms = app.tonemap_prms;
    scn.add_skyenv   = app.add_skyenv;
    scn.task_queue.emplace_back(app_task_type::load_scene);
    app.selected = (int)app.scenes.size() - 1;
}

void draw_opengl_widgets(const opengl_window& win) {
    static string load_path = "", save_path = "", error_message = "";
    auto&         app = *(app_state*)get_opengl_user_pointer(win);
    if (!begin_opengl_widgets_window(win, "yitrace")) return;
    if (!app.errors.empty() && error_message.empty()) {
        error_message = app.errors.front();
        app.errors.pop_front();
        open_modal_opengl_widget(win, "error");
    }
    if (!draw_modal_message_opengl_window(win, "error", error_message)) {
        error_message = "";
    }
    if (draw_modal_fileialog_opengl_widgets(
            win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
        add_new_scene(app, load_path);
    }
    if (draw_modal_fileialog_opengl_widgets(win, "save", save_path, true,
            get_dirname(save_path), get_filename(save_path),
            "*.yaml;*.obj;*.pbrt")) {
        app.scenes[app.selected].outname = save_path;
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::save_scene);
        save_path = "";
    }
    if (draw_modal_fileialog_opengl_widgets(win, "save image", save_path, true,
            get_dirname(save_path), get_filename(save_path),
            "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
        app.scenes[app.selected].imagename = save_path;
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::save_image);
        save_path = "";
    }
    if (draw_button_opengl_widget(win, "load")) {
        open_modal_opengl_widget(win, "load");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save",
            app.selected >= 0 && app.scenes[app.selected].task_queue.empty())) {
        save_path = app.scenes[app.selected].outname;
        open_modal_opengl_widget(win, "save");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "save image",
            app.selected >= 0 && app.scenes[app.selected].render_done)) {
        save_path = app.scenes[app.selected].imagename;
        open_modal_opengl_widget(win, "save image");
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "close", app.selected >= 0)) {
        app.scenes[app.selected].task_queue.emplace_back(
            app_task_type::close_scene);
    }
    continue_opengl_widget_line(win);
    if (draw_button_opengl_widget(win, "quit")) {
        set_close_opengl_window(win, true);
    }
    if (app.scenes.empty()) return;
    draw_combobox_opengl_widget(
        win, "scene", app.selected, (int)app.scenes.size(),
        [&app](int idx) { return app.scenes[idx].name.c_str(); }, false);
    auto& scn = app.scenes[app.selected];
    if (begin_header_opengl_widget(win, "trace")) {
        auto cam_names = vector<string>();
        for (auto& camera : scn.scene.cameras) cam_names.push_back(camera.uri);
        auto trace_prms = scn.trace_prms;
        if (scn.load_done) {
            if (draw_combobox_opengl_widget(
                    win, "camera", trace_prms.camera_id, cam_names)) {
            }
        }
        draw_slider_opengl_widget(
            win, "width", trace_prms.image_size.x, 0, 4096);
        draw_slider_opengl_widget(
            win, "height", trace_prms.image_size.y, 0, 4096);
        draw_slider_opengl_widget(
            win, "nsamples", trace_prms.num_samples, 16, 4096);
        draw_combobox_opengl_widget(
            win, "tracer", (int&)trace_prms.sampler_type, trace_sampler_names);
        draw_combobox_opengl_widget(win, "false color",
            (int&)trace_prms.falsecolor_type, trace_falsecolor_names);
        draw_slider_opengl_widget(
            win, "nbounces", trace_prms.max_bounces, 1, 128);
        draw_slider_opengl_widget(
            win, "seed", (int&)trace_prms.random_seed, 0, 1000000);
        draw_slider_opengl_widget(win, "pratio", scn.preview_ratio, 1, 64);
        auto tonemap_prms = scn.tonemap_prms;
        draw_slider_opengl_widget(
            win, "exposure", tonemap_prms.exposure, -5, 5);
        draw_checkbox_opengl_widget(win, "filmic", tonemap_prms.filmic);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "srgb", tonemap_prms.srgb);
        if (trace_prms != scn.trace_prms) {
            scn.task_queue.emplace_back(app_task_type::apply_edit,
                app_edit{typeid(trace_params), -1, trace_prms, false});
        }
        if (tonemap_prms != scn.tonemap_prms) {
            scn.task_queue.emplace_back(app_task_type::apply_edit,
                app_edit{typeid(tonemap_params), -1, tonemap_prms, false});
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "inspect")) {
        draw_label_opengl_widget(win, "scene", get_filename(scn.filename));
        draw_label_opengl_widget(win, "filename", scn.filename);
        draw_label_opengl_widget(win, "outname", scn.outname);
        draw_label_opengl_widget(win, "imagename", scn.imagename);
        draw_label_opengl_widget(win, "image", "%d x %d @ %d",
            scn.render.size().x, scn.render.size().y, scn.render_sample);
        draw_slider_opengl_widget(win, "zoom", scn.image_scale, 0.1, 10);
        draw_checkbox_opengl_widget(win, "zoom to fit", scn.zoom_to_fit);
        continue_opengl_widget_line(win);
        draw_checkbox_opengl_widget(win, "fps", scn.navigation_fps);
        if (draw_button_opengl_widget(win, "print cams")) {
            for (auto& camera : scn.scene.cameras) {
                print_obj_camera(camera);
            }
        }
        continue_opengl_widget_line(win);
        if (draw_button_opengl_widget(win, "print stats")) {
            print_info("{}", format_stats(scn.scene).c_str());
            print_info("{}", print_stats(scn.bvh).c_str());
        }
        auto mouse_pos = get_opengl_mouse_pos(win);
        auto ij        = get_image_coords(
            mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
        draw_dragger_opengl_widget(win, "mouse", ij);
        if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
            ij.y < scn.render.size().y) {
            draw_coloredit_opengl_widget(
                win, "pixel", scn.render[{ij.x, ij.y}]);
        } else {
            auto zero4f_ = zero4f;
            draw_coloredit_opengl_widget(win, "pixel", zero4f_);
        }
        end_header_opengl_widget(win);
    }
    if (scn.load_done && begin_header_opengl_widget(win, "scene tree")) {
        draw_opengl_widgets_scene_tree(win, "", scn.scene, scn.selection, 200);
        end_header_opengl_widget(win);
    }
    if (scn.load_done && begin_header_opengl_widget(win, "scene object")) {
        auto edit = app_edit{};
        if (draw_opengl_widgets_scene_inspector(
                win, "", scn.scene, scn.selection, edit, 200)) {
            scn.task_queue.emplace_back(app_task_type::apply_edit, edit);
        }
        end_header_opengl_widget(win);
    }
    if (begin_header_opengl_widget(win, "log")) {
        draw_log_opengl_widget(win);
        end_header_opengl_widget(win);
    }
}

void draw(const opengl_window& win) {
    auto& app      = *(app_state*)get_opengl_user_pointer(win);
    auto  win_size = get_opengl_window_size(win);
    set_opengl_viewport(get_opengl_framebuffer_viewport(win));
    clear_opengl_framebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
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

void apply_edit(const string& filename, yocto_scene& scene,
    trace_params& trace_prms, tonemap_params& tonemap_prms,
    bool& reload_element, bool& updated_lights, bool& updated_bvh,
    const app_edit& edit) {
    auto& [type, index, data, reload] = edit;

    if (type == typeid(yocto_camera)) {
        scene.cameras[index] = any_cast<yocto_camera>(data);
    } else if (type == typeid(yocto_texture)) {
        scene.textures[index] = any_cast<yocto_texture>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_voltexture)) {
        scene.voltextures[index] = any_cast<yocto_voltexture>(data);
        if (reload) reload_element = true;
    } else if (type == typeid(yocto_shape)) {
        scene.shapes[index] = any_cast<yocto_shape>(data);
        if (reload) {
            reload_element = true;
            updated_bvh    = true;
        }
    } else if (type == typeid(yocto_subdiv)) {
        scene.subdivs[index] = any_cast<yocto_subdiv>(data);
        if (reload) {
            reload_element = true;
            updated_bvh    = true;
        }
    } else if (type == typeid(yocto_material)) {
        auto old_emission      = scene.materials[index].emission_factor;
        scene.materials[index] = any_cast<yocto_material>(data);
        if (old_emission != scene.materials[index].emission_factor) {
            updated_lights = true;
        }
    } else if (type == typeid(yocto_instance)) {
        auto old_instance      = scene.instances[index];
        scene.instances[index] = any_cast<yocto_instance>(data);
        if (old_instance.shape != scene.instances[index].shape ||
            old_instance.frame != scene.instances[index].frame) {
            updated_bvh = true;
        }
    } else if (type == typeid(yocto_environment)) {
        auto old_emission         = scene.materials[index].emission_factor;
        scene.environments[index] = any_cast<yocto_environment>(data);
        if (old_emission != scene.materials[index].emission_factor) {
            updated_lights = true;
        }
    } else if (type == typeid(trace_params)) {
        trace_prms = any_cast<trace_params>(data);
    } else if (type == typeid(tonemap_params)) {
        tonemap_prms = any_cast<tonemap_params>(data);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }
}

// reload an element
void load_element(
    const string& filename, yocto_scene& scene, const app_edit& edit) {
    auto& [type, index, data, reload] = edit;

    if (type == typeid(yocto_texture)) {
        auto& texture = scene.textures[index];
        load_image(get_dirname(filename) + texture.uri, texture.hdr_image,
            texture.ldr_image);
    } else if (type == typeid(yocto_voltexture)) {
        auto& texture = scene.voltextures[index];
        load_volume(get_dirname(filename) + texture.uri, texture.volume_data);
    } else if (type == typeid(yocto_shape)) {
        auto& shape = scene.shapes[index];
        load_shape(get_dirname(filename) + shape.uri, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.quads_positions,
            shape.quads_normals, shape.quads_texcoords, shape.positions,
            shape.normals, shape.texcoords, shape.colors, shape.radius, false);
    } else if (type == typeid(yocto_subdiv)) {
        // TODO: this needs more fixing?
        auto& subdiv = scene.subdivs[index];
        load_shape(get_dirname(filename) + subdiv.uri, subdiv.points,
            subdiv.lines, subdiv.triangles, subdiv.quads,
            subdiv.quads_positions, subdiv.quads_normals,
            subdiv.quads_texcoords, subdiv.positions, subdiv.normals,
            subdiv.texcoords, subdiv.colors, subdiv.radius,
            subdiv.preserve_facevarying);
        tesselate_subdiv(scene, scene.subdivs[index]);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }
}

void refit_bvh(const string& filename, yocto_scene& scene, bvh_scene& bvh,
    const bvh_params& bvh_prms, const app_edit& edit) {
    auto& [type, index, data, reload] = edit;

    auto updated_shapes    = vector<int>{};
    auto updated_instances = vector<int>{};
    if (type == typeid(yocto_shape)) {
        updated_shapes.push_back(index);
    } else if (type == typeid(yocto_subdiv)) {
        auto& subdiv = scene.subdivs[index];
        updated_shapes.push_back(subdiv.tesselated_shape);
    } else if (type == typeid(yocto_instance)) {
        updated_instances.push_back(index);
    } else {
        throw runtime_error("unsupported type "s + type.name());
    }

    refit_bvh(scene, bvh, updated_instances, updated_shapes);
}

void update(app_state& app) {
    // close if needed
    while (!app.scenes.empty()) {
        auto pos = -1;
        for (auto idx = 0; idx < app.scenes.size(); idx++) {
            for (auto& task : app.scenes[idx].task_queue) {
                if (task.type == app_task_type::close_scene) pos = idx;
            }
        }
        if (pos < 0) break;
        app.scenes.erase(app.scenes.begin() + pos);
        app.selected = app.scenes.empty() ? -1 : 0;
    }

    // consume partial results
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (task.type != app_task_type::render_image || task.queue.empty())
            continue;
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
            scn.name = format("{} [{}x{}@{}]", get_filename(scn.filename),
                scn.render.size().x, scn.render.size().y, scn.render_sample);
        }
    }

    // remove unneeded tasks
    for (auto& scn : app.scenes) {
        while (scn.task_queue.size() > 1) {
            auto& task = scn.task_queue.at(0);
            auto& next = scn.task_queue.at(1);
            if (task.type == app_task_type::render_image) {
                if (next.type != app_task_type::render_image &&
                    next.type != app_task_type::apply_edit)
                    break;
                log_info("cancel rendering {}", scn.filename);
            } else if (task.type == app_task_type::apply_edit) {
                if (next.type != app_task_type::apply_edit ||
                    task.edit.type != next.edit.type ||
                    task.edit.index != next.edit.index)
                    break;
                log_info("cancel editing {}", scn.filename);
            } else {
                break;
            }
            task.stop = true;
            if (task.result.valid()) {
                try {
                    task.result.get();
                } catch (...) {
                }
            }
            scn.task_queue.pop_front();
        }
    }

    // apply synchronous edit
    for (auto& scn : app.scenes) {
        while (!scn.task_queue.empty()) {
            auto& task = scn.task_queue.front();
            if (task.type != app_task_type::apply_edit) break;
            log_info("start editing {}", scn.filename);
            try {
                scn.render_done     = false;
                auto reload_element = false, update_bvh = false,
                     update_lights = false;
                apply_edit(scn.filename, scn.scene, scn.trace_prms,
                    scn.tonemap_prms, reload_element, update_lights, update_bvh,
                    task.edit);
                log_info("done editing {}", scn.filename);
                if (reload_element) {
                    scn.load_done = false;
                    scn.task_queue.emplace_back(
                        app_task_type::load_element, task.edit);
                }
                if (update_bvh) {
                    scn.bvh_done = false;
                    scn.task_queue.emplace_back(
                        app_task_type::refit_bvh, task.edit);
                }
                if (update_lights) {
                    scn.lights_done = false;
                    scn.task_queue.emplace_back(app_task_type::init_lights);
                }
                scn.task_queue.emplace_back(app_task_type::render_image);
            } catch (std::exception& e) {
                log_error(e.what());
                app.errors.push_back("cannot edit " + scn.filename);
            }
            scn.task_queue.pop_front();
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
            case app_task_type::close_scene: break;
            case app_task_type::load_scene: {
                try {
                    task.result.get();
                    scn.load_done  = true;
                    scn.image_size = camera_image_size(
                        scn.scene.cameras[scn.trace_prms.camera_id],
                        scn.trace_prms.image_size);
                    scn.render.resize(scn.image_size);
                    scn.display.resize(scn.image_size);
                    scn.preview.resize(scn.image_size);
                    scn.name = format("{} [{}x{}@0]",
                        get_filename(scn.filename), scn.render.size().x,
                        scn.render.size().y);
                    log_info("done loading {}", scn.filename);
                    init_opengl_texture(
                        scn.gl_txt, scn.display, false, false, false);
                    scn.task_queue.emplace_back(app_task_type::build_bvh);
                    scn.task_queue.emplace_back(app_task_type::init_lights);
                    scn.task_queue.emplace_back(app_task_type::render_image);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot load " + scn.filename);
                }
            } break;
            case app_task_type::load_element: {
                try {
                    task.result.get();
                    scn.load_done = true;
                    log_info("done loading element from {}", scn.filename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back(
                        "cannot load element from " + scn.filename);
                }
            } break;
            case app_task_type::build_bvh: {
                try {
                    task.result.get();
                    scn.bvh_done = true;
                    scn.name     = format("{}", get_filename(scn.filename));
                    log_info("done building bvh {}", scn.filename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot build bvh " + scn.filename);
                }
            } break;
            case app_task_type::refit_bvh: {
                try {
                    task.result.get();
                    scn.bvh_done = true;
                    scn.name     = format("{}", get_filename(scn.filename));
                    log_info("done refitting bvh {}", scn.filename);
                } catch (std::exception& e) {
                    log_error(e.what());
                    scn.name = format("{} [error]", get_filename(scn.filename));
                    app.errors.push_back("cannot refit bvh " + scn.filename);
                }
            } break;
            case app_task_type::init_lights: {
                try {
                    task.result.get();
                    scn.lights_done = true;
                    scn.name        = format("{}", get_filename(scn.filename));
                    log_info("done building lights {}", scn.filename);
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
            case app_task_type::render_image: {
                try {
                    task.result.get();
                    scn.render_done = true;
                    log_info("done rendering {}", scn.filename);
                    scn.render_sample = scn.trace_prms.num_samples;
                    scn.name          = format("{} [{}x{}@{}]",
                        get_filename(scn.filename), scn.render.size().x,
                        scn.render.size().y, scn.render_sample);
                } catch (std::exception& e) {
                    log_error(e.what());
                    app.errors.push_back("cannot render " + scn.filename);
                }
            } break;
            case app_task_type::apply_edit: break;
        }
        scn.task_queue.pop_front();
    }
    // schedule tasks not running
    for (auto& scn : app.scenes) {
        if (scn.task_queue.empty()) continue;
        auto& task = scn.task_queue.front();
        if (task.result.valid()) continue;
        task.stop = false;
        switch (task.type) {
            case app_task_type::none: break;
            case app_task_type::close_scene: break;
            case app_task_type::load_scene: {
                log_info("start loading {}", scn.filename);
                scn.load_done   = false;
                scn.bvh_done    = false;
                scn.lights_done = false;
                task.result     = async([&scn]() {
                    load_scene(scn.filename, scn.scene, scn.sceneio_prms);
                    tesselate_subdivs(scn.scene);
                    if (scn.add_skyenv) add_sky(scn.scene);
                });
            } break;
            case app_task_type::load_element: {
                log_info("start loading element for {}", scn.filename);
                scn.load_done = false;
                task.result   = async([&scn, &task]() {
                    load_element(scn.filename, scn.scene, task.edit);
                });
            } break;
            case app_task_type::build_bvh: {
                log_info("start building bvh {}", scn.filename);
                scn.bvh_done = false;
                task.result  = async(
                    [&scn]() { build_bvh(scn.scene, scn.bvh, scn.bvh_prms); });
            } break;
            case app_task_type::refit_bvh: {
                log_info("start refitting bvh {}", scn.filename);
                scn.bvh_done = false;
                task.result  = async([&scn, &task]() {
                    refit_bvh(scn.filename, scn.scene, scn.bvh, scn.bvh_prms,
                        task.edit);
                });
            } break;
            case app_task_type::init_lights: {
                log_info("start building lights {}", scn.filename);
                scn.lights_done = false;
                task.result     = async(
                    [&scn]() { init_trace_lights(scn.lights, scn.scene); });
            } break;
            case app_task_type::save_image: {
                log_info("start saving {}", scn.imagename);
                task.result = async([&scn]() {
                    save_tonemapped(
                        scn.imagename, scn.render, scn.tonemap_prms);
                });
            } break;
            case app_task_type::save_scene: {
                log_info("start saving {}", scn.outname);
                task.result = async([&scn]() {
                    save_scene(scn.outname, scn.scene, scn.sceneio_prms);
                });
            } break;
            case app_task_type::render_image: {
                log_info("start rendering {}", scn.filename);
                scn.render_done = false;
                scn.image_size  = camera_image_size(
                    scn.scene.cameras[scn.trace_prms.camera_id],
                    scn.trace_prms.image_size);
                if (scn.lights.instances.empty() &&
                    scn.lights.environments.empty() &&
                    is_sampler_lit(scn.trace_prms)) {
                    log_info(
                        "no lights presents, switching to eyelight shader");
                    scn.trace_prms.sampler_type = trace_sampler_type::eyelight;
                }
                scn.render_sample = 0;
                scn.name = format("{} [{}x{}@{}]", get_filename(scn.filename),
                    scn.render.size().x, scn.render.size().y,
                    scn.render_sample);
                task.result = async([&scn, &task]() {
                    update_app_render(scn.filename, scn.render, scn.display,
                        scn.preview, scn.state, scn.scene, scn.lights, scn.bvh,
                        scn.trace_prms, scn.tonemap_prms, scn.preview_ratio,
                        task.stop, task.current, task.queue);
                });
                if (scn.render.size() != scn.image_size) {
                    scn.render.resize(scn.image_size);
                    scn.display.assign(scn.image_size);
                    scn.preview.assign(scn.image_size);
                    if (scn.gl_txt) {
                        delete_opengl_texture(scn.gl_txt);
                        scn.gl_txt = {};
                    }
                    init_opengl_texture(
                        scn.gl_txt, scn.display, false, false, false);
                }
            } break;
            case app_task_type::apply_edit: break;
        }
    }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
    auto& app = *(app_state*)get_opengl_user_pointer(win);
    for (auto& path : paths) add_new_scene(app, path);
}

// run ui loop
void run_ui(app_state& app) {
    // window
    auto win = opengl_window();
    init_opengl_window(win, {1280 + 320, 720}, "yitrace", &app, draw);
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
            auto& scn        = app.scenes[app.selected];
            auto& old_camera = scn.scene.cameras.at(scn.trace_prms.camera_id);
            auto  camera     = scn.scene.cameras.at(scn.trace_prms.camera_id);
            auto  dolly      = 0.0f;
            auto  pan        = zero2f;
            auto  rotate     = zero2f;
            if (mouse_left && !shift_down)
                rotate = (mouse_pos - last_pos) / 100.0f;
            if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
            if (mouse_left && shift_down)
                pan = (mouse_pos - last_pos) * camera.focus_distance / 200.0f;
            pan.x = -pan.x;
            update_camera_turntable(
                camera.frame, camera.focus_distance, rotate, dolly, pan);
            if (camera.frame != old_camera.frame ||
                camera.focus_distance != old_camera.focus_distance) {
                scn.task_queue.emplace_back(app_task_type::apply_edit,
                    app_edit{typeid(yocto_camera), scn.trace_prms.camera_id,
                        camera, false});
            }
        }

        // selection
        if (app.selected >= 0 && app.scenes[app.selected].load_done &&
            (mouse_left || mouse_right) && alt_down && !widgets_active) {
            auto& scn = app.scenes[app.selected];
            auto  ij  = get_image_coords(mouse_pos, scn.image_center,
                scn.image_scale, scn.render.size());
            if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
                ij.y < scn.render.size().y) {
                auto& camera = scn.scene.cameras.at(scn.trace_prms.camera_id);
                auto  ray    = eval_camera(
                    camera, ij, scn.render.size(), {0.5f, 0.5f}, zero2f);
                if (auto isec = bvh_intersection{};
                    intersect_bvh(scn.scene, scn.bvh, ray, isec)) {
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
    app.trace_prms.samples_per_batch = 1;
    auto no_parallel                 = false;
    auto filenames                   = vector<string>{};

    // names for enums
    auto trace_sampler_type_namemap = std::map<string, trace_sampler_type>{};
    for (auto type = 0; type < trace_sampler_names.size(); type++) {
        trace_sampler_type_namemap[trace_sampler_names[type]] =
            (trace_sampler_type)type;
    }
    auto trace_falsecolor_type_namemap =
        std::map<string, trace_falsecolor_type>{};
    for (auto type = 0; type < trace_falsecolor_names.size(); type++) {
        trace_falsecolor_type_namemap[trace_falsecolor_names[type]] =
            (trace_falsecolor_type)type;
    }

    // parse command line
    auto parser = CLI::App{"progressive path tracing"};
    parser.add_option("--camera", app.trace_prms.camera_id, "Camera index.");
    parser.add_option("--hres,-R", app.trace_prms.image_size.x,
        "Image horizontal resolution.");
    parser.add_option(
        "--vres,-r", app.trace_prms.image_size.y, "Image vertical resolution.");
    parser.add_option(
        "--nsamples,-s", app.trace_prms.num_samples, "Number of samples.");
    parser
        .add_option("--tracer,-t", app.trace_prms.sampler_type, "Tracer type.")
        ->transform(CLI::IsMember(trace_sampler_type_namemap));
    parser
        .add_option("--falsecolor,-F", app.trace_prms.falsecolor_type,
            "Tracer false color type.")
        ->transform(CLI::IsMember(trace_falsecolor_type_namemap));
    parser.add_option(
        "--nbounces", app.trace_prms.max_bounces, "Maximum number of bounces.");
    parser.add_option(
        "--pixel-clamp", app.trace_prms.pixel_clamp, "Final pixel clamping.");
    parser.add_option("--seed", app.trace_prms.random_seed,
        "Seed for the random number generators.");
    parser.add_flag("--env-hidden,!--no-env-hidden",
        app.trace_prms.environments_hidden,
        "Environments are hidden in renderer");
    parser.add_flag("--parallel,!--no-parallel", no_parallel,
        "Disable parallel execution.");
    parser.add_option(
        "--exposure,-e", app.tonemap_prms.exposure, "Hdr exposure");
    parser.add_flag(
        "--filmic,!--no-filmic", app.tonemap_prms.filmic, "Hdr filmic");
    parser.add_flag("--srgb,!--no-srgb", app.tonemap_prms.srgb, "Hdr srgb");
    parser.add_flag("--bvh-high-quality,!--no-bvh-high-quality",
        app.bvh_prms.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
    parser.add_flag("--bvh-embree,!--no-bvh-embree", app.bvh_prms.use_embree,
        "Use Embree ratracer");
    parser.add_flag("--bvh-embree-flatten,!--no-bvh-embree-flatten",
        app.bvh_prms.embree_flatten, "Flatten embree scene");
    parser.add_flag("--bvh-embree-compact,!--no-bvh-embree-compact",
        app.bvh_prms.embree_compact, "Embree runs in compact memory");
#endif
    parser.add_flag(
        "--add-skyenv,!--no-add-skyenv", app.add_skyenv, "Add sky envmap");
    parser.add_option("scenes", filenames, "Scene filenames")->required(true);
    try {
        parser.parse(argc, argv);
    } catch (const CLI::ParseError& e) {
        return parser.exit(e);
    }

    // fix parallel code
    if (no_parallel) {
        app.bvh_prms.run_serially     = true;
        app.sceneio_prms.run_serially = true;
        app.trace_prms.run_serially   = true;
    }

    // loading images
    for (auto filename : filenames) add_new_scene(app, filename);
    app.selected = app.scenes.empty() ? -1 : 0;

    // run interactive
    run_ui(app);

    // done
    return 0;
}
