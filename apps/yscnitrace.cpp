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

#include <atomic>
#include <future>
#include <map>
#include <thread>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

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
  app_task_type            type;
  std::future<void>        result;
  std::atomic<bool>        stop;
  std::atomic<int>         current;
  std::deque<image_region> queue;
  std::mutex               queuem;
  app_edit                 edit;

  app_task(app_task_type type, const app_edit& edit = {})
      : type{type}, result{}, stop{false}, current{-1}, edit{edit} {}
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
  load_params    load_prms     = {};
  save_params    save_prms     = {};
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
  std::deque<app_task> task_queue;
  app_selection        selection = {typeid(void), -1};
};

// Application state
struct app_state {
  // data
  std::deque<app_scene> scenes;
  int                   selected = -1;
  std::deque<string>    errors;

  // default options
  load_params    load_prms    = {};
  save_params    save_prms    = {};
  bvh_params     bvh_prms     = {};
  trace_params   trace_prms   = {};
  tonemap_params tonemap_prms = {};
  bool           add_skyenv   = false;
};

void update_app_render(const string& filename, image<vec4f>& render,
    image<vec4f>& display, image<vec4f>& preview, trace_state& state,
    const yocto_scene& scene, const trace_lights& lights, const bvh_scene& bvh,
    const trace_params& trace_prms, const tonemap_params& tonemap_prms,
    int preview_ratio, std::atomic<bool>* cancel,
    std::atomic<int>& current_sample, std::deque<image_region>& queue,
    std::mutex& queuem) {
  auto preview_prms = trace_prms;
  preview_prms.resolution /= preview_ratio;
  preview_prms.samples = 1;
  auto small_preview   = trace_image(scene, bvh, lights, preview_prms);
  auto display_preview = tonemap(small_preview, tonemap_prms);
  for (auto j = 0; j < preview.size().y; j++) {
    for (auto i = 0; i < preview.size().x; i++) {
      auto pi = clamp(i / preview_ratio, 0, display_preview.size().x - 1),
           pj = clamp(j / preview_ratio, 0, display_preview.size().y - 1);
      preview[{i, j}] = display_preview[{pi, pj}];
    }
  }
  {
    std::lock_guard guard{queuem};
    queue.push_back({{0, 0}, {0, 0}});
  }
  current_sample = 0;

  auto& camera     = scene.cameras.at(trace_prms.camera);
  auto  image_size = camera_resolution(camera, trace_prms.resolution);
  state            = make_trace_state(image_size, trace_prms.seed);
  auto regions     = make_regions(render.size(), trace_prms.region, true);

  for (auto sample = 0; sample < trace_prms.samples;
       sample += trace_prms.batch) {
    if (cancel && *cancel) return;
    current_sample   = sample;
    auto num_samples = min(
        trace_prms.batch, trace_prms.samples - current_sample);
    auto                futures  = vector<std::future<void>>{};
    auto                nthreads = std::thread::hardware_concurrency();
    std::atomic<size_t> next_idx(0);
    for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
      futures.emplace_back(std::async(std::launch::async,
          [num_samples, &trace_prms, &tonemap_prms, &render, &display, &scene,
              &lights, &bvh, &state, &queue, &queuem, &next_idx, cancel,
              &regions]() {
            while (true) {
              if (cancel && *cancel) break;
              auto idx = next_idx.fetch_add(1);
              if (idx >= regions.size()) break;
              auto region = regions[idx];
              trace_region(render, state, scene, bvh, lights, region,
                  num_samples, trace_prms);
              tonemap(display, render, region, tonemap_prms);
              {
                std::lock_guard guard{queuem};
                queue.push_back(region);
              }
            }
          }));
    }
    current_sample = trace_prms.samples;
  }
}

void add_new_scene(app_state& app, const string& filename) {
  auto& scn        = app.scenes.emplace_back();
  scn.filename     = filename;
  scn.imagename    = fs::path(filename).replace_extension(".png");
  scn.outname      = fs::path(filename).replace_extension(".edited.yaml");
  scn.name         = fs::path(scn.filename).filename();
  scn.load_prms    = app.load_prms;
  scn.save_prms    = app.save_prms;
  scn.trace_prms   = app.trace_prms;
  scn.bvh_prms     = app.bvh_prms;
  scn.tonemap_prms = app.tonemap_prms;
  scn.add_skyenv   = app.add_skyenv;
  scn.task_queue.emplace_back(app_task_type::load_scene);
  app.selected = (int)app.scenes.size() - 1;
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app = *(app_state*)get_gluser_pointer(win);
  if (!begin_glwidgets_window(win, "yscnitrace")) return;
  if (!app.errors.empty() && error_message.empty()) {
    error_message = app.errors.front();
    app.errors.pop_front();
    open_glmodal(win, "error");
  }
  if (!draw_glmessage(win, "error", error_message)) {
    error_message = "";
  }
  if (draw_glfiledialog(
          win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
    add_new_scene(app, load_path);
  }
  if (draw_glfiledialog(win, "save", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.yaml;*.obj;*.pbrt")) {
    app.scenes[app.selected].outname = save_path;
    app.scenes[app.selected].task_queue.emplace_back(app_task_type::save_scene);
    save_path = "";
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    app.scenes[app.selected].imagename = save_path;
    app.scenes[app.selected].task_queue.emplace_back(app_task_type::save_image);
    save_path = "";
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save",
          app.selected >= 0 && app.scenes[app.selected].task_queue.empty())) {
    save_path = app.scenes[app.selected].outname;
    open_glmodal(win, "save");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save image",
          app.selected >= 0 && app.scenes[app.selected].render_done)) {
    save_path = app.scenes[app.selected].imagename;
    open_glmodal(win, "save image");
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", app.selected >= 0)) {
    app.scenes[app.selected].task_queue.emplace_back(
        app_task_type::close_scene);
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.scenes.empty()) return;
  draw_glcombobox(
      win, "scene", app.selected, (int)app.scenes.size(),
      [&app](int idx) { return app.scenes[idx].name.c_str(); }, false);
  auto& scn = app.scenes[app.selected];
  if (begin_glheader(win, "trace")) {
    auto cam_names = vector<string>();
    for (auto& camera : scn.scene.cameras) cam_names.push_back(camera.uri);
    auto trace_prms = scn.trace_prms;
    if (scn.load_done) {
      if (draw_glcombobox(win, "camera", trace_prms.camera, cam_names)) {
      }
    }
    draw_glslider(win, "resolution", trace_prms.resolution, 180, 4096);
    draw_glslider(win, "nsamples", trace_prms.samples, 16, 4096);
    draw_glcombobox(
        win, "tracer", (int&)trace_prms.sampler, trace_sampler_names);
    draw_glcombobox(win, "false color", (int&)trace_prms.falsecolor,
        trace_falsecolor_names);
    draw_glslider(win, "nbounces", trace_prms.bounces, 1, 128);
    draw_glcheckbox(win, "env hidden", trace_prms.envhidden);
    continue_glline(win);
    draw_glcheckbox(win, "filter", trace_prms.tentfilter);
    draw_glslider(win, "seed", (int&)trace_prms.seed, 0, 1000000);
    draw_glslider(win, "pratio", scn.preview_ratio, 1, 64);
    auto tonemap_prms = scn.tonemap_prms;
    draw_glslider(win, "exposure", tonemap_prms.exposure, -5, 5);
    draw_glcheckbox(win, "filmic", tonemap_prms.filmic);
    continue_glline(win);
    draw_glcheckbox(win, "srgb", tonemap_prms.srgb);
    if (trace_prms != scn.trace_prms) {
      scn.task_queue.emplace_back(app_task_type::apply_edit,
          app_edit{typeid(trace_params), -1, trace_prms, false});
    }
    if (tonemap_prms != scn.tonemap_prms) {
      scn.task_queue.emplace_back(app_task_type::apply_edit,
          app_edit{typeid(tonemap_params), -1, tonemap_prms, false});
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    draw_gllabel(win, "scene", fs::path(scn.filename).filename());
    draw_gllabel(win, "filename", scn.filename);
    draw_gllabel(win, "outname", scn.outname);
    draw_gllabel(win, "imagename", scn.imagename);
    draw_gllabel(win, "image", "%d x %d @ %d", scn.render.size().x,
        scn.render.size().y, scn.render_sample);
    draw_glslider(win, "zoom", scn.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", scn.zoom_to_fit);
    continue_glline(win);
    draw_glcheckbox(win, "fps", scn.navigation_fps);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : scn.scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      printf("%s\n", format_stats(scn.scene).c_str());
      printf("%s\n", format_stats(scn.bvh).c_str());
    }
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
    draw_gldragger(win, "mouse", ij);
    if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
        ij.y < scn.render.size().y) {
      draw_glcoloredit(win, "pixel", scn.render[{ij.x, ij.y}]);
    } else {
      auto zero4f_ = zero4f;
      draw_glcoloredit(win, "pixel", zero4f_);
    }
    end_glheader(win);
  }
  if (scn.load_done && begin_glheader(win, "scene tree")) {
    draw_glscenetree(win, "", scn.scene, scn.selection, 200);
    end_glheader(win);
  }
  if (scn.load_done && begin_glheader(win, "scene object")) {
    auto edit = app_edit{};
    if (draw_glsceneinspector(win, "", scn.scene, scn.selection, edit, 200)) {
      scn.task_queue.emplace_back(app_task_type::apply_edit, edit);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win) {
  auto& app      = *(app_state*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  set_glviewport(get_glframebuffer_viewport(win));
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.scenes.empty() && app.selected >= 0) {
    auto& scn = app.scenes.at(app.selected);
    if (scn.load_done && scn.gl_txt) {
      update_imview(scn.image_center, scn.image_scale, scn.display.size(),
          win_size, scn.zoom_to_fit);
      draw_glimage_background(scn.gl_txt, win_size.x, win_size.y,
          scn.image_center, scn.image_scale);
      set_glblending(true);
      draw_glimage(scn.gl_txt, win_size.x, win_size.y, scn.image_center,
          scn.image_scale);
      set_glblending(false);
    }
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
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
    auto old_emission      = scene.materials[index].emission;
    scene.materials[index] = any_cast<yocto_material>(data);
    if (old_emission != scene.materials[index].emission) {
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
    auto old_emission         = scene.materials[index].emission;
    scene.environments[index] = any_cast<yocto_environment>(data);
    if (old_emission != scene.materials[index].emission) {
      updated_lights = true;
    }
  } else if (type == typeid(trace_params)) {
    trace_prms = any_cast<trace_params>(data);
  } else if (type == typeid(tonemap_params)) {
    tonemap_prms = any_cast<tonemap_params>(data);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }
}

// reload an element
void load_element(
    const string& filename, yocto_scene& scene, const app_edit& edit) {
  auto& [type, index, data, reload] = edit;

  if (type == typeid(yocto_texture)) {
    auto& texture = scene.textures[index];
    if (is_hdr_filename(texture.uri)) {
      load_image(fs::path(filename).parent_path() / texture.uri, texture.hdr);
    } else {
      load_imageb(fs::path(filename).parent_path() / texture.uri, texture.ldr);
    }
  } else if (type == typeid(yocto_voltexture)) {
    auto& texture = scene.voltextures[index];
    load_volume(fs::path(filename).parent_path() / texture.uri, texture.vol);
  } else if (type == typeid(yocto_shape)) {
    auto& shape = scene.shapes[index];
    load_shape(fs::path(filename).parent_path() / shape.uri, shape.points,
        shape.lines, shape.triangles, shape.quads, shape.quadspos,
        shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
        shape.texcoords, shape.colors, shape.radius, false);
  } else if (type == typeid(yocto_subdiv)) {
    // TODO: this needs more fixing?
    auto& subdiv = scene.subdivs[index];
    load_shape(fs::path(filename).parent_path() / subdiv.uri, subdiv.points,
        subdiv.lines, subdiv.triangles, subdiv.quads, subdiv.quadspos,
        subdiv.quadsnorm, subdiv.quadstexcoord, subdiv.positions,
        subdiv.normals, subdiv.texcoords, subdiv.colors, subdiv.radius,
        subdiv.facevarying);
    tesselate_subdiv(scene, scene.subdivs[index]);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
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
    updated_shapes.push_back(subdiv.shape);
  } else if (type == typeid(yocto_instance)) {
    updated_instances.push_back(index);
  } else {
    throw std::runtime_error("unsupported type "s + type.name());
  }

  refit_bvh(bvh, scene, updated_instances, updated_shapes, bvh_prms);
}

void update(const opengl_window& win, app_state& app) {
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
    if (task.type != app_task_type::render_image) continue;
    auto updated = false;
    while (true) {
      auto region = image_region{};
      {
        std::lock_guard guard{task.queuem};
        if (task.queue.empty()) break;
        region = task.queue.front();
        task.queue.pop_front();
      }
      if (region.size() == zero2i) {
        update_gltexture_region(
            scn.gl_txt, scn.preview, {zero2i, scn.preview.size()}, false);
      } else {
        update_gltexture_region(scn.gl_txt, scn.display, region, false);
      }
      updated = true;
    }
    if (updated) {
      scn.render_sample = max(scn.render_sample, (int)task.current);
      scn.name          = fs::path(scn.filename).filename().string() + "[" +
                 std::to_string(scn.render.size().x) + "x" +
                 std::to_string(scn.render.size().y) + " @ " +
                 std::to_string(scn.render_sample) + "]";
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
        log_glinfo(win, "cancel rendering " + scn.filename);
      } else if (task.type == app_task_type::apply_edit) {
        if (next.type != app_task_type::apply_edit ||
            task.edit.type != next.edit.type ||
            task.edit.index != next.edit.index)
          break;
        log_glinfo(win, "cancel editing " + scn.filename);
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
      log_glinfo(win, "start editing " + scn.filename);
      try {
        scn.render_done     = false;
        auto reload_element = false, update_bvh = false, update_lights = false;
        apply_edit(scn.filename, scn.scene, scn.trace_prms, scn.tonemap_prms,
            reload_element, update_lights, update_bvh, task.edit);
        log_glinfo(win, "done editing " + scn.filename);
        if (reload_element) {
          scn.load_done = false;
          scn.task_queue.emplace_back(app_task_type::load_element, task.edit);
        }
        if (update_bvh) {
          scn.bvh_done = false;
          scn.task_queue.emplace_back(app_task_type::refit_bvh, task.edit);
        }
        if (update_lights) {
          scn.lights_done = false;
          scn.task_queue.emplace_back(app_task_type::init_lights);
        }
        scn.task_queue.emplace_back(app_task_type::render_image);
      } catch (std::exception& e) {
        log_glerror(win, e.what());
        app.errors.push_back("cannot edit " + scn.filename);
      }
      scn.task_queue.pop_front();
    }
  }

  // grab result of finished tasks
  for (auto& scn : app.scenes) {
    if (scn.task_queue.empty()) continue;
    auto& task = scn.task_queue.front();
    if (!task.result.valid()) continue;
    if (task.type == app_task_type::render_image) {
      std::lock_guard guard{task.queuem};
      if (!task.queue.empty()) continue;
    }
    if (task.result.wait_for(std::chrono::nanoseconds(10)) !=
        std::future_status::ready)
      continue;
    switch (task.type) {
      case app_task_type::none: break;
      case app_task_type::close_scene: break;
      case app_task_type::load_scene: {
        try {
          task.result.get();
          scn.load_done  = true;
          scn.image_size = camera_resolution(
              scn.scene.cameras[scn.trace_prms.camera],
              scn.trace_prms.resolution);
          scn.render.resize(scn.image_size);
          scn.display.resize(scn.image_size);
          scn.preview.resize(scn.image_size);
          scn.name = fs::path(scn.filename).filename().string() + " [" +
                     std::to_string(scn.render.size().x) + "x" +
                     std::to_string(scn.render.size().y) + " @ 0]";
          log_glinfo(win, "done loading " + scn.filename);
          init_gltexture(scn.gl_txt, scn.display, false, false, false);
          scn.task_queue.emplace_back(app_task_type::build_bvh);
          scn.task_queue.emplace_back(app_task_type::init_lights);
          scn.task_queue.emplace_back(app_task_type::render_image);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [errpr]";
          app.errors.push_back("cannot load " + scn.filename);
        }
      } break;
      case app_task_type::load_element: {
        try {
          task.result.get();
          scn.load_done = true;
          log_glinfo(win, "done loading element from " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [errpr]";
          app.errors.push_back("cannot load element from " + scn.filename);
        }
      } break;
      case app_task_type::build_bvh: {
        try {
          task.result.get();
          scn.bvh_done = true;
          scn.name     = fs::path(scn.filename).filename().string();
          log_glinfo(win, "done building bvh " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [errpr]";
          app.errors.push_back("cannot build bvh " + scn.filename);
        }
      } break;
      case app_task_type::refit_bvh: {
        try {
          task.result.get();
          scn.bvh_done = true;
          scn.name     = fs::path(scn.filename).filename().string();
          log_glinfo(win, "done refitting bvh " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [errpr]";
          app.errors.push_back("cannot refit bvh " + scn.filename);
        }
      } break;
      case app_task_type::init_lights: {
        try {
          task.result.get();
          scn.lights_done = true;
          scn.name        = fs::path(scn.filename).filename().string();
          log_glinfo(win, "done building lights " + scn.filename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          scn.name = fs::path(scn.filename).filename().string() + " [errpr]";
          app.errors.push_back("cannot build lights " + scn.filename);
        }
      } break;
      case app_task_type::save_image: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.imagename);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.imagename);
        }
      } break;
      case app_task_type::save_scene: {
        try {
          task.result.get();
          log_glinfo(win, "done saving " + scn.outname);
        } catch (std::exception& e) {
          log_glerror(win, e.what());
          app.errors.push_back("cannot save " + scn.outname);
        }
      } break;
      case app_task_type::render_image: {
        try {
          task.result.get();
          scn.render_done = true;
          log_glinfo(win, "done rendering " + scn.filename);
          scn.render_sample = scn.trace_prms.samples;
          scn.name = fs::path(scn.filename).filename().string() + " [" +
                     std::to_string(scn.render.size().x) + "x" +
                     std::to_string(scn.render.size().y) + " @ " +
                     std::to_string(scn.render_sample) + "]";
        } catch (std::exception& e) {
          log_glerror(win, e.what());
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
        log_glinfo(win, "start loading " + scn.filename);
        scn.load_done   = false;
        scn.bvh_done    = false;
        scn.lights_done = false;
        task.result     = std::async(std::launch::async, [&scn]() {
          load_scene(scn.filename, scn.scene, scn.load_prms);
          tesselate_subdivs(scn.scene);
          if (scn.add_skyenv) add_sky(scn.scene);
        });
      } break;
      case app_task_type::load_element: {
        log_glinfo(win, "start loading element for " + scn.filename);
        scn.load_done = false;
        task.result   = std::async(std::launch::async, [&scn, &task]() {
          load_element(scn.filename, scn.scene, task.edit);
        });
      } break;
      case app_task_type::build_bvh: {
        log_glinfo(win, "start building bvh " + scn.filename);
        scn.bvh_done = false;
        task.result  = std::async(std::launch::async,
            [&scn]() { make_bvh(scn.bvh, scn.scene, scn.bvh_prms); });
      } break;
      case app_task_type::refit_bvh: {
        log_glinfo(win, "start refitting bvh " + scn.filename);
        scn.bvh_done = false;
        task.result  = std::async(std::launch::async, [&scn, &task]() {
          refit_bvh(scn.filename, scn.scene, scn.bvh, scn.bvh_prms, task.edit);
        });
      } break;
      case app_task_type::init_lights: {
        log_glinfo(win, "start building lights " + scn.filename);
        scn.lights_done = false;
        task.result     = std::async(std::launch::async,
            [&scn]() { scn.lights = make_trace_lights(scn.scene); });
      } break;
      case app_task_type::save_image: {
        log_glinfo(win, "start saving " + scn.imagename);
        task.result = std::async(std::launch::async, [&scn]() {
          if (is_hdr_filename(scn.imagename)) {
            save_image(scn.imagename, scn.render);
          } else {
            save_imageb(scn.imagename, tonemapb(scn.render, scn.tonemap_prms));
          }
        });
      } break;
      case app_task_type::save_scene: {
        log_glinfo(win, "start saving " + scn.outname);
        task.result = std::async(std::launch::async,
            [&scn]() { save_scene(scn.outname, scn.scene, scn.save_prms); });
      } break;
      case app_task_type::render_image: {
        log_glinfo(win, "start rendering " + scn.filename);
        scn.render_done = false;
        scn.image_size  = camera_resolution(
            scn.scene.cameras[scn.trace_prms.camera],
            scn.trace_prms.resolution);
        if (scn.lights.instances.empty() && scn.lights.environments.empty() &&
            is_sampler_lit(scn.trace_prms)) {
          log_glinfo(win, "no lights presents, switching to eyelight shader");
          scn.trace_prms.sampler = trace_params::sampler_type::eyelight;
        }
        scn.render_sample = 0;
        scn.name          = fs::path(scn.filename).filename().string() + " [" +
                   std::to_string(scn.render.size().x) + "x" +
                   std::to_string(scn.render.size().y) + " @ " +
                   std::to_string(scn.render_sample) + "]";
        task.result = std::async(std::launch::async, [&scn, &task]() {
          update_app_render(scn.filename, scn.render, scn.display, scn.preview,
              scn.state, scn.scene, scn.lights, scn.bvh, scn.trace_prms,
              scn.tonemap_prms, scn.preview_ratio, &task.stop, task.current,
              task.queue, task.queuem);
        });
        if (scn.render.size() != scn.image_size) {
          scn.render.resize(scn.image_size);
          scn.display.assign(scn.image_size);
          scn.preview.assign(scn.image_size);
          if (scn.gl_txt) {
            delete_gltexture(scn.gl_txt);
            scn.gl_txt = {};
          }
          init_gltexture(scn.gl_txt, scn.display, false, false, false);
        }
      } break;
      case app_task_type::apply_edit: break;
    }
  }
}

void drop_callback(const opengl_window& win, const vector<string>& paths) {
  auto& app = *(app_state*)get_gluser_pointer(win);
  for (auto& path : paths) add_new_scene(app, path);
}

// run ui loop
void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", &app, draw);
  set_drop_glcallback(win, drop_callback);

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
    if (app.selected >= 0 && app.scenes[app.selected].load_done &&
        (mouse_left || mouse_right) && !alt_down && !widgets_active) {
      auto& scn        = app.scenes[app.selected];
      auto& old_camera = scn.scene.cameras.at(scn.trace_prms.camera);
      auto  camera     = scn.scene.cameras.at(scn.trace_prms.camera);
      auto  dolly      = 0.0f;
      auto  pan        = zero2f;
      auto  rotate     = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      if (camera.frame != old_camera.frame ||
          camera.focus != old_camera.focus) {
        scn.task_queue.emplace_back(app_task_type::apply_edit,
            app_edit{
                typeid(yocto_camera), scn.trace_prms.camera, camera, false});
      }
    }

    // selection
    if (app.selected >= 0 && app.scenes[app.selected].load_done &&
        (mouse_left || mouse_right) && alt_down && !widgets_active) {
      auto& scn = app.scenes[app.selected];
      auto  ij  = get_image_coords(
          mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
      if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
          ij.y < scn.render.size().y) {
        auto& camera = scn.scene.cameras.at(scn.trace_prms.camera);
        auto  ray    = eval_camera(
            camera, ij, scn.render.size(), {0.5f, 0.5f}, zero2f);
        if (auto isec = intersect_bvh(scn.bvh, ray); isec.hit) {
          scn.selection = {typeid(yocto_instance), isec.instance};
        }
      }
    }

    // update
    update(win, app);

    // draw
    draw(win);

    // event hadling
    process_glevents(win);
  }

  // clear
  delete_glwindow(win);
}

int main(int argc, const char* argv[]) {
  // application
  app_state app{};
  app.trace_prms.batch = 1;
  auto no_parallel     = false;
  auto filenames       = vector<string>{};

  // names for enums
  auto sampler_namemap = std::map<string, trace_params::sampler_type>{};
  for (auto type = 0; type < trace_sampler_names.size(); type++) {
    sampler_namemap[trace_sampler_names[type]] =
        (trace_params::sampler_type)type;
  }
  auto falsecolor_namemap = std::map<string, trace_params::falsecolor_type>{};
  for (auto type = 0; type < trace_falsecolor_names.size(); type++) {
    falsecolor_namemap[trace_falsecolor_names[type]] =
        (trace_params::falsecolor_type)type;
  }

  // parse command line
  auto cli = make_cli("yscnitrace", "progressive path tracing");
  add_cli_option(cli, "--camera", app.trace_prms.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app.trace_prms.resolution, "Image resolution.");
  add_cli_option(
      cli, "--samples,-s", app.trace_prms.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)app.trace_prms.sampler,
      "Tracer type.", trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)app.trace_prms.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", app.trace_prms.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", app.trace_prms.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", app.trace_prms.tentfilter, "Filter image.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", app.trace_prms.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(cli, "--parallel,/--no-parallel", no_parallel,
      "Disable parallel execution.");
  add_cli_option(
      cli, "--exposure,-e", app.tonemap_prms.exposure, "Hdr exposure");
  add_cli_option(
      cli, "--filmic/--no-filmic", app.tonemap_prms.filmic, "Hdr filmic");
  add_cli_option(cli, "--srgb/--no-srgb", app.tonemap_prms.srgb, "Hdr srgb");
  add_cli_option(cli, "--bvh-high-quality/--no-bvh-high-quality",
      app.bvh_prms.high_quality, "Use high quality bvh mode");
#if YOCTO_EMBREE
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", app.bvh_prms.use_embree,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-flatten/--no-bvh-embree-flatten",
      app.bvh_prms.embree_flatten, "Flatten embree scene");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      app.bvh_prms.embree_compact, "Embree runs in compact memory");
#endif
  add_cli_option(cli, "--add-skyenv", app.add_skyenv, "Add sky envmap");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (no_parallel) {
    app.bvh_prms.noparallel   = true;
    app.load_prms.noparallel  = true;
    app.save_prms.noparallel  = true;
    app.trace_prms.noparallel = true;
  }

  // loading images
  for (auto filename : filenames) add_new_scene(app, filename);
  app.selected = app.scenes.empty() ? -1 : 0;

  // run interactive
  run_ui(app);

  // done
  return 0;
}
