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
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_gui.h>
using namespace yocto::math;
namespace sio = yocto::sceneio;
namespace cli = yocto::commonio;
namespace gui = yocto::gui;

#include <atomic>
#include <deque>
#include <future>
using namespace std::string_literals;

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto::sceneio {
void print_obj_camera(sio::camera* camera);
};

// Application state
struct app_state {
  // loading parameters
  std::string filename  = "scene.json";
  std::string imagename = "out.png";
  std::string outname   = "scene.json";
  std::string name      = "";

  // options
  gui::scene_params drawgl_prms = {};

  // scene
  sio::model*  ioscene  = new sio::model{};
  sio::camera* iocamera = nullptr;

  // rendering state
  gui::scene*  glscene  = new gui::scene{};
  gui::camera* glcamera = nullptr;

  // editing
  sio::camera*      selected_camera      = nullptr;
  sio::object*      selected_object      = nullptr;
  sio::instance*    selected_instance    = nullptr;
  sio::shape*       selected_shape       = nullptr;
  sio::subdiv*      selected_subdiv      = nullptr;
  sio::material*    selected_material    = nullptr;
  sio::environment* selected_environment = nullptr;
  sio::texture*     selected_texture     = nullptr;

  // loading status
  std::atomic<bool> ok           = false;
  std::future<void> loader       = {};
  std::string       status       = "";
  std::string       error        = "";
  std::atomic<int>  current      = 0;
  std::atomic<int>  total        = 0;
  std::string       loader_error = "";

  ~app_state() {
    if (ioscene) delete ioscene;
    if (glscene) delete glscene;
  }
};

void update_lights(gui::scene* glscene, sio::model* ioscene) {
  clear_lights(glscene);
  for (auto ioobject : ioscene->objects) {
    if (has_max_lights(glscene)) break;
    if (ioobject->material->emission == zero3f) continue;
    auto ioshape = ioobject->shape;
    auto bbox    = invalidb3f;
    for (auto p : ioshape->positions) bbox = merge(bbox, p);
    auto pos  = (bbox.max + bbox.min) / 2;
    auto area = 0.0f;
    if (!ioshape->triangles.empty()) {
      for (auto t : ioshape->triangles)
        area += triangle_area(ioshape->positions[t.x], ioshape->positions[t.y],
            ioshape->positions[t.z]);
    } else if (!ioshape->quads.empty()) {
      for (auto q : ioshape->quads)
        area += quad_area(ioshape->positions[q.x], ioshape->positions[q.y],
            ioshape->positions[q.z], ioshape->positions[q.w]);
    } else if (!ioshape->lines.empty()) {
      for (auto l : ioshape->lines)
        area += line_length(ioshape->positions[l.x], ioshape->positions[l.y]);
    } else {
      area += ioshape->positions.size();
    }
    auto ke = ioobject->material->emission * area;
    set_light(add_light(glscene), transform_point(ioobject->frame, pos), ke,
        gui::light_type::point, false);
  }
}

void init_glscene(gui::scene* glscene, sio::model* ioscene,
    gui::camera*& glcamera, sio::camera* iocamera,
    sio::progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->materials.size() +
             (int)ioscene->textures.size() + (int)ioscene->shapes.size() +
             (int)ioscene->subdivs.size() + (int)ioscene->instances.size() +
             (int)ioscene->objects.size()};

  // create scene
  init_scene(glscene);

  // camera
  auto camera_map     = std::unordered_map<sio::camera*, gui::camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto camera = add_camera(glscene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_nearfar(camera, 0.001, 10000);
    camera_map[iocamera] = camera;
  }

  // textures
  auto texture_map     = std::unordered_map<sio::texture*, gui::texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto gltexture = add_texture(glscene);
    if (!iotexture->colorf.empty()) {
      set_texture(gltexture, iotexture->colorf);
    } else if (!iotexture->colorb.empty()) {
      set_texture(gltexture, iotexture->colorb);
    } else if (!iotexture->scalarf.empty()) {
      set_texture(gltexture, iotexture->scalarf);
    } else if (!iotexture->scalarb.empty()) {
      set_texture(gltexture, iotexture->scalarb);
    }
    texture_map[iotexture] = gltexture;
  }

  // material
  auto material_map     = std::unordered_map<sio::material*, gui::material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto glmaterial = add_material(glscene);
    set_emission(glmaterial, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
        texture_map.at(iomaterial->color_tex));
    set_specular(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_metallic(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_roughness(glmaterial, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(glmaterial, iomaterial->opacity,
        texture_map.at(iomaterial->opacity_tex));
    set_normalmap(glmaterial, texture_map.at(iomaterial->normal_tex));
    material_map[iomaterial] = glmaterial;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    if (progress_cb) progress_cb("convert subdiv", progress.x++, progress.y);
    tesselate_subdiv(ioscene, iosubdiv);
  }

  // shapes
  auto shape_map     = std::unordered_map<sio::shape*, gui::shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
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
    shape_map[ioshape] = glshape;
  }

  // instances
  auto instance_map     = std::unordered_map<sio::instance*, gui::instance*>{};
  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto glinstance = add_instance(glscene);
    set_frames(glinstance, ioinstance->frames);
    instance_map[ioinstance] = glinstance;
  }

  // shapes
  for (auto ioobject : ioscene->objects) {
    if (progress_cb) progress_cb("convert object", progress.x++, progress.y);
    auto globject = add_object(glscene);
    set_frame(globject, ioobject->frame);
    set_shape(globject, shape_map.at(ioobject->shape));
    set_material(globject, material_map.at(ioobject->material));
    set_instance(globject, instance_map.at(ioobject->instance));
  }

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);

  // get cmmera
  glcamera = camera_map.at(iocamera);
}

int main(int argc, const char* argv[]) {
  // initialize app
  auto app_guard   = std::make_unique<app_state>();
  auto app         = app_guard.get();
  auto camera_name = ""s;

  // parse command line
  auto cli = cli::make_cli("ysceneviews", "views scene inteactively");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(
      cli, "--resolution,-r", app->drawgl_prms.resolution, "Image resolution.");
  add_option(cli, "--shading", app->drawgl_prms.shading, "Eyelight rendering.",
      gui::shading_names);
  add_option(cli, "scene", app->filename, "Scene filename", true);
  parse_cli(cli, argc, argv);

  // loading scene
  auto ioerror = ""s;
  if (!load_scene(app->filename, app->ioscene, ioerror, cli::print_progress))
    cli::print_fatal(ioerror);

  // get camera
  app->iocamera = get_camera(app->ioscene, camera_name);

  // callbacks
  auto callbacks    = gui::ui_callbacks{};
  callbacks.init_cb = [app](gui::window* win, const gui::input& input) {
    init_glscene(app->glscene, app->ioscene, app->glcamera, app->iocamera,
        [app](const std::string& message, int current, int total) {
          app->status  = "init scene";
          app->current = current;
          app->total   = total;
        });
  };
  callbacks.clear_cb = [app](gui::window* win, const gui::input& input) {
    clear_scene(app->glscene);
  };
  callbacks.draw_cb = [app](gui::window* win, const gui::input& input) {
    if (app->drawgl_prms.shading == gui::shading_type::lights)
      update_lights(app->glscene, app->ioscene);
    draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
        app->drawgl_prms);
  };
  callbacks.widgets_cb = [app](gui::window* win, const gui::input& input) {
    draw_progressbar(win, app->status.c_str(), app->current, app->total);
    if (draw_combobox(win, "camera", app->iocamera, app->ioscene->cameras)) {
      for (auto idx = 0; idx < app->ioscene->cameras.size(); idx++) {
        if (app->ioscene->cameras[idx] == app->iocamera)
          app->glcamera = app->glscene->cameras[idx];
      }
    }
    auto& params = app->drawgl_prms;
    draw_slider(win, "resolution", params.resolution, 0, 4096);
    draw_checkbox(win, "wireframe", params.wireframe);
    draw_combobox(win, "shading", (int&)params.shading, gui::shading_names);
    continue_line(win);
    draw_checkbox(win, "edges", params.edges);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
  };
  callbacks.update_cb = [app](gui::window* win, const gui::input& input) {
    // update(win, apps);
  };
  callbacks.uiupdate_cb = [app](gui::window* win, const gui::input& input) {
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
          app->iocamera->frame, app->iocamera->focus, rotate, dolly, pan);
      set_frame(app->glcamera, app->iocamera->frame);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "ysceneviews", callbacks);

  // done
  return 0;
}
