//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_gui.h>
using namespace yocto;

#include <atomic>
#include <deque>
#include <future>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(scene_camera* camera);
};

// Application state
struct app_state {
  // loading parameters
  string filename  = "shape.obj";
  string imagename = "out.png";
  string outname   = "out.obj";
  string name      = "";

  // options
  ogl_scene_params drawgl_prms = {};

  // scene
  generic_shape* ioshape = new generic_shape{};

  // rendering state
  ogl_scene*  glscene  = new ogl_scene{};
  ogl_camera* glcamera = nullptr;

  // loading status
  std::atomic<bool> ok           = false;
  std::future<void> loader       = {};
  string            status       = "";
  string            error        = "";
  std::atomic<int>  current      = 0;
  std::atomic<int>  total        = 0;
  string            loader_error = "";

  ~app_state() {
    if (glscene) delete glscene;
    if (ioshape) delete ioshape;
  }
};

// Application state
struct app_states {
  // data
  vector<app_state*>     states   = {};
  app_state*             selected = nullptr;
  std::deque<app_state*> loading  = {};

  // default options
  ogl_scene_params drawgl_prms = {};

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

void load_shape_async(
    app_states* apps, const string& filename, const string& camera_name = "") {
  auto app         = apps->states.emplace_back(new app_state{});
  app->filename    = filename;
  app->imagename   = sfs::path(filename).replace_extension(".png");
  app->outname     = sfs::path(filename).replace_extension(".edited.obj");
  app->name        = sfs::path(app->filename).filename();
  app->drawgl_prms = apps->drawgl_prms;
  app->status      = "load";
  app->loader      = std::async(std::launch::async, [app, camera_name]() {
    auto progress_cb = [app](const string& message, int current, int total) {
      app->current = current;
      app->total   = total;
    };
    if (!load_shape(app->filename, *app->ioshape, app->loader_error)) return;
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = app;
}

// TODO(fabio): move this function to math
frame3f camera_frame(
    const bbox3f& bbox, float lens, float aspect, float film = 0.036) {
  auto center      = (bbox.max + bbox.min) / 2;
  auto bbox_radius = length(bbox.max - bbox.min) / 2;
  auto camera_dir  = vec3f{0, 0, 1};
  auto camera_dist = bbox_radius * lens / (film / aspect);
  camera_dist *= 2.0f;  // correction
  return lookat_frame(camera_dir * camera_dist + center, center, {0, 1, 0});
}

void init_glscene(ogl_scene* glscene, const generic_shape* ioshape,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 4};

  // create scene
  init_scene(glscene);

  // compute bounding box
  auto bbox = invalidb3f;
  for (auto& pos : ioshape->positions) bbox = merge(bbox, pos);
  // TODO(fabio): this should be a math function

  // camera
  if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
  auto glcamera = add_camera(glscene);
  set_frame(glcamera, camera_frame(bbox, 0.050, 16.0f / 9.0f, 0.036));
  set_lens(glcamera, 0.050, 16.0f / 9.0f, 0.036);
  set_nearfar(glcamera, 0.001, 10000);
  glcamera->focus = length(glcamera->frame.o - center(bbox));

  // material
  if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
  auto glmaterial = add_material(glscene);
  set_emission(glmaterial, {0, 0, 0}, nullptr);
  set_color(glmaterial, {0.5, 1, 0.5}, nullptr);
  set_specular(glmaterial, 1, nullptr);
  set_metallic(glmaterial, 0, nullptr);
  set_roughness(glmaterial, 0.2, nullptr);
  set_opacity(glmaterial, 1, nullptr);
  set_normalmap(glmaterial, nullptr);

  // shapes
  if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
  auto glshape = add_shape(glscene);
  set_positions(glshape, ioshape->positions);
  set_normals(glshape, ioshape->normals);
  set_texcoords(glshape, ioshape->texcoords);
  set_colors(glshape, ioshape->colors);
  set_points(glshape, ioshape->points);
  set_lines(glshape, ioshape->lines);
  set_triangles(glshape, ioshape->triangles);
  set_quads(glshape, ioshape->quads);
  set_edges(glshape, ioshape->triangles, ioshape->quads);

  // shapes
  if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
  auto globject = add_object(glscene);
  set_frame(globject, identity3x4f);
  set_shape(globject, glshape);
  set_material(globject, glmaterial);

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);
}

bool draw_widgets(gui_window* win, scene_model* ioscene, scene_shape* ioshape) {
  if (!ioshape) return false;
  auto edited = 0;
  draw_label(win, "name", ioshape->name);
  draw_label(win, "points", std::to_string(ioshape->points.size()));
  draw_label(win, "lines", std::to_string(ioshape->lines.size()));
  draw_label(win, "triangles", std::to_string(ioshape->triangles.size()));
  draw_label(win, "quads", std::to_string(ioshape->quads.size()));
  draw_label(win, "positions", std::to_string(ioshape->positions.size()));
  draw_label(win, "normals", std::to_string(ioshape->normals.size()));
  draw_label(win, "texcoords", std::to_string(ioshape->texcoords.size()));
  draw_label(win, "colors", std::to_string(ioshape->colors.size()));
  draw_label(win, "radius", std::to_string(ioshape->radius.size()));
  draw_label(win, "tangents", std::to_string(ioshape->tangents.size()));
  draw_label(win, "quads pos", std::to_string(ioshape->quadspos.size()));
  draw_label(win, "quads norm", std::to_string(ioshape->quadsnorm.size()));
  draw_label(
      win, "quads texcoord", std::to_string(ioshape->quadstexcoord.size()));
  edited += draw_slider(win, "subdivisions", ioshape->subdivisions, 0, 5);
  edited += draw_checkbox(win, "catmull-clark", ioshape->catmullclark);
  edited += draw_slider(win, "displacement", ioshape->displacement, 0, 1);
  edited += draw_combobox(win, "displacement_tex", ioshape->displacement_tex,
      ioscene->textures, true);
  return edited;
}

// draw with shading
void draw_widgets(gui_window* win, app_states* apps, const gui_input& input) {
  static auto load_path = ""s, save_path = ""s, error_message = ""s;
  if (draw_filedialog_button(win, "load", true, "load", load_path, false, "./",
          "", "*.ply;*.obj")) {
    load_shape_async(apps, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", apps->selected && apps->selected->ok,
          "save", save_path, true, sfs::path(save_path).parent_path(),
          sfs::path(save_path).filename(), "*.ply;*.obj")) {
    auto app     = apps->selected;
    app->outname = save_path;
    auto ok      = save_shape(app->outname, *app->ioshape, app->error);
    save_path    = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) return;
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    set_close(win, true);
  }
  if (apps->states.empty()) return;
  draw_combobox(win, "shape", apps->selected, apps->states, false);
  if (!apps->selected) return;
  draw_progressbar(win, apps->selected->status.c_str(), apps->selected->current,
      apps->selected->total);
  if (apps->selected->error != "") {
    draw_label(win, "error", apps->selected->error);
    return;
  }
  if (!apps->selected->ok) return;
  auto app = apps->selected;
  if (begin_header(win, "view")) {
    auto& params = app->drawgl_prms;
    draw_slider(win, "resolution", params.resolution, 0, 4096);
    draw_combobox(win, "shading", (int&)params.shading, ogl_shading_names);
    draw_checkbox(win, "wireframe", params.wireframe);
    continue_line(win);
    draw_checkbox(win, "edges", params.edges);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
    end_header(win);
  }
  if (begin_header(win, "inspect")) {
    draw_label(win, "shape", sfs::path(app->filename).filename());
    draw_label(win, "filename", app->filename);
    draw_label(win, "outname", app->outname);
    draw_label(win, "imagename", app->imagename);
    end_header(win);
  }
}

// draw with shading
void draw(gui_window* win, app_states* apps, const gui_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app = apps->selected;
  draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
      app->drawgl_prms);
}

// update
void update(gui_window* win, app_states* apps) {
  auto is_ready = [](const std::future<void>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    auto progress_cb = [app](const string& message, int current, int total) {
      app->current = current;
      app->total   = total;
    };
    app->loader.get();
    if (app->loader_error.empty()) {
      init_glscene(app->glscene, app->ioshape, progress_cb);
      app->glcamera = app->glscene->cameras.front();
      app->ok       = true;
      app->status   = "ok";
    } else {
      app->status = "error";
      app->error  = app->loader_error;
    }
  }
}

int main(int argc, const char* argv[]) {
  // initialize app
  auto apps_guard  = std::make_unique<app_states>();
  auto apps        = apps_guard.get();
  auto filenames   = vector<string>{};
  auto camera_name = ""s;

  // parse command line
  auto cli = make_cli("yshapeview", "views shapes inteactively");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(cli, "--resolution,-r", apps->drawgl_prms.resolution,
      "Image resolution.");
  add_option(cli, "--shading", apps->drawgl_prms.shading, "Shading type.",
      ogl_shading_names);
  add_option(cli, "shapes", filenames, "Shape filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_shape_async(apps, filename, camera_name);

  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [apps](gui_window* win, const gui_input& input) {
    for (auto app : apps->states) clear_scene(app->glscene);
  };
  callbacks.draw_cb = [apps](gui_window* win, const gui_input& input) {
    draw(win, apps, input);
  };
  callbacks.widgets_cb = [apps](gui_window* win, const gui_input& input) {
    draw_widgets(win, apps, input);
  };
  callbacks.drop_cb = [apps](gui_window* win, const vector<string>& paths,
                          const gui_input& input) {
    for (auto& path : paths) load_shape_async(apps, path);
  };
  callbacks.update_cb = [apps](gui_window* win, const gui_input& input) {
    update(win, apps);
  };
  callbacks.uiupdate_cb = [apps](gui_window* win, const gui_input& input) {
    if (!apps->selected || !apps->selected->ok) return;
    auto app = apps->selected;

    // handle mouse and keyboard for navigation
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      update_turntable(
          app->glcamera->frame, app->glcamera->focus, rotate, dolly, pan);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "ysceneview", callbacks);

  // done
  return 0;
}
