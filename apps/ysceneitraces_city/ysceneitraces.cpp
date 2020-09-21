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
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_gui.h>
#include <yocto_pathtrace/yocto_pathtrace.h>
#include "ext/filesystem.hpp"

using namespace yocto::math;
namespace ptr = yocto::pathtrace;
namespace cli = yocto::commonio;
namespace img = yocto::image;
namespace gui = yocto::gui;
namespace sio = yocto::sceneio;
namespace sfs = ghc::filesystem;
namespace shp = yocto::shape;

#include <future>
#include <memory>
using namespace std::string_literals;

// Application state
struct app_state {
  // loading options
  std::string filename  = "scene.yaml";
  std::string imagename = "out.png";
  std::string outname   = "out.yaml";
  std::string name      = "";

  // options
  ptr::trace_params params = {};
  bool add_lens  = false;
  std::string error = "";

  // scene
  ptr::scene*              scene        = new ptr::scene{};
  ptr::camera*             camera       = nullptr;
  std::vector<std::string> camera_names = {};

  // rendering state
  img::image<vec4f> render   = {};
  img::image<vec4f> display  = {};
  float             exposure = 0;

  // view scene
  gui::image*       glimage  = new gui::image{};
  gui::image_params glparams = {};

  // computation
  int               render_sample  = 0;
  int               render_counter = 0;
  ptr::state*       render_state   = new ptr::state{};
  std::future<void> render_worker  = {};
  std::atomic<bool> render_stop    = {};

  // status
  std::atomic<float> progress = 0;

  ~app_state() {
    if (render_worker.valid()) {
      render_stop = true;
      render_worker.get();
    }
    if (render_state) delete render_state;
    if (scene) delete scene;
    if (glimage) delete glimage;
  }
};

// construct a scene from io
void init_scene(ptr::scene* scene, sio::model* ioscene, ptr::camera*& camera,
    sio::camera* iocamera, sio::progress_callback print_progress = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->subdivs.size() +
             (int)ioscene->objects.size()};

  auto camera_map     = std::unordered_map<sio::camera*, ptr::camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (print_progress)
      print_progress("convert camera", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
    camera_map[iocamera] = camera;
  }

  auto texture_map     = std::unordered_map<sio::texture*, ptr::texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (print_progress)
      print_progress("convert texture", progress.x++, progress.y);
    auto texture = add_texture(scene);
    if (!iotexture->colorf.empty()) {
      set_texture(texture, iotexture->colorf);
    } else if (!iotexture->colorb.empty()) {
      set_texture(texture, iotexture->colorb);
    } else if (!iotexture->scalarf.empty()) {
      set_texture(texture, iotexture->scalarf);
    } else if (!iotexture->scalarb.empty()) {
      set_texture(texture, iotexture->scalarb);
    }
    texture_map[iotexture] = texture;
  }

  auto material_map     = std::unordered_map<sio::material*, ptr::material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (print_progress)
      print_progress("convert material", progress.x++, progress.y);
    auto material = add_material(scene);
    set_emission(material, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(
        material, iomaterial->color, texture_map.at(iomaterial->color_tex));
    set_specular(material, iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_ior(material, iomaterial->ior);
    set_metallic(material, iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_transmission(material, iomaterial->transmission, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->transmission_tex));
    set_roughness(material, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(
        material, iomaterial->opacity, texture_map.at(iomaterial->opacity_tex));
    set_thin(material, iomaterial->thin);
    set_scattering(material, iomaterial->scattering, iomaterial->scanisotropy,
        texture_map.at(iomaterial->scattering_tex));
    set_normalmap(material, texture_map.at(iomaterial->normal_tex));
    material_map[iomaterial] = material;
  }

  auto shape_map     = std::unordered_map<sio::shape*, ptr::shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (print_progress)
      print_progress("convert shape", progress.x++, progress.y);
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    if (!ioshape->quads.empty())
      set_triangles(shape, shp::quads_to_triangles(ioshape->quads));
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_radius(shape, ioshape->radius);
    shape_map[ioshape] = shape;
  }

  auto subdiv_map = std::unordered_map<sio::subdiv*, ptr::shape*>{};
  for (auto iosubdiv : ioscene->subdivs) {
    if (print_progress)
      print_progress("convert subdiv", progress.x++, progress.y);
    auto subdiv = add_shape(scene);
    set_subdiv_quadspos(subdiv, iosubdiv->quadspos);
    set_subdiv_quadstexcoord(subdiv, iosubdiv->quadstexcoord);
    set_subdiv_positions(subdiv, iosubdiv->positions);
    set_subdiv_texcoords(subdiv, iosubdiv->texcoords);
    subdiv_map[iosubdiv] = subdiv;
  }

  for (auto ioobject : ioscene->objects) {
    if (print_progress)
      print_progress("convert object", progress.x++, progress.y);
    if (ioobject->instance) {
      for (auto frame : ioobject->instance->frames) {
        auto object = add_object(scene);
        set_frame(object, frame * ioobject->frame);
        if (ioobject->shape) set_shape(object, shape_map.at(ioobject->shape));
        if (ioobject->subdiv) {
          set_shape(object, subdiv_map.at(ioobject->subdiv));
          set_subdiv_subdivision(subdiv_map.at(ioobject->subdiv),
              ioobject->material->subdivisions, true);
          set_subdiv_displacement(subdiv_map.at(ioobject->subdiv),
              ioobject->material->displacement,
              texture_map.at(ioobject->material->displacement_tex));
        }
        set_material(object, material_map.at(ioobject->material));
      }
    } else {
      auto object = add_object(scene);
      set_frame(object, ioobject->frame);
      if (ioobject->shape) set_shape(object, shape_map.at(ioobject->shape));
      if (ioobject->subdiv) {
        set_shape(object, subdiv_map.at(ioobject->subdiv));
        set_subdiv_subdivision(subdiv_map.at(ioobject->subdiv),
            ioobject->material->subdivisions, true);
        set_subdiv_displacement(subdiv_map.at(ioobject->subdiv),
            ioobject->material->displacement,
            texture_map.at(ioobject->material->displacement_tex));
      }
      set_material(object, material_map.at(ioobject->material));
    }
  }

  for (auto ioenvironment : ioscene->environments) {
    if (print_progress)
      print_progress("convert environment", progress.x++, progress.y);
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // done
  if (print_progress) print_progress("convert done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

// init camera names
void init_camera_names(std::vector<std::string>& names,
    const std::vector<sio::camera*>&             iocameras) {
  for (auto iocamera : iocameras) {
    names.push_back(iocamera->name);
  }
}

void reset_display(app_state* app) {
  // stop render
  app->render_stop = true;
  if (app->render_worker.valid()) app->render_worker.get();

  // init state
  init_state(app->render_state, app->scene, app->camera, app->params);
  if (app->render.size() != app->render_state->render.size()) {
    app->render.resize(app->render_state->render.size());
    app->display.resize(app->render_state->render.size());
  }

  // render preview
  auto pstate_guard = std::make_unique<ptr::state>();
  auto pstate       = pstate_guard.get();
  auto pprms        = app->params;
  pprms.resolution /= app->params.pratio;
  pprms.samples = 1;
  init_state(pstate, app->scene, app->camera, pprms);
  trace_samples(pstate, app->scene, app->camera, pprms);
  for (auto j = 0; j < app->render.size().y; j++) {
    for (auto i = 0; i < app->render.size().x; i++) {
      auto pi = clamp(i / app->params.pratio, 0, pstate->render.size().x - 1),
           pj = clamp(j / app->params.pratio, 0, pstate->render.size().y - 1);
      app->render[{i, j}]  = pstate->render[{pi, pj}];
      app->display[{i, j}] = tonemap(app->render[{i, j}], app->exposure);
    }
  }

  // start render
  app->render_counter = 0;
  app->render_stop    = false;
  app->render_worker  = std::async(std::launch::async, [app]() {
    for (auto sample = 0; sample < app->params.samples; sample++) {
      app->progress = (float)sample / (float)app->params.samples;
      trace_samples(app->render_state, app->scene, app->camera, app->params,
          &app->render_stop);
      for (auto j = 0; j < app->render_state->render.size().y; j++) {
        for (auto i = 0; i < app->render_state->render.size().x; i++) {
          if (app->render_stop) return;
          app->render[{i, j}]  = app->render_state->render[{i, j}];
          app->display[{i, j}] = tonemap(app->render[{i, j}], app->exposure);
        }
      }
    }
  });
  app->progress       = 1;
}

int main(int argc, const char* argv[]) {
  // application
  auto app_guard = std::make_unique<app_state>();
  auto app       = app_guard.get();

  // command line options
  auto camera_name = ""s;
  bool add_skyenv = false;

  // parse command line
  auto cli = cli::make_cli("yscnitraces", "progressive path tracing");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(
      cli, "--resolution,-r", app->params.resolution, "Image resolution.");
  add_option(cli, "--samples,-s", app->params.samples, "Number of samples.");
  add_option(cli, "--shader,-t", app->params.shader, "Tracer type.",
      ptr::shader_names);
  add_option(
      cli, "--bounces,-b", app->params.bounces, "Maximum number of bounces.");
  add_option(cli, "--clamp", app->params.clamp, "Final pixel clamping.");
  add_option(cli, "--skyenv/--no-skyenv", add_skyenv, "Add sky envmap");
  add_option(cli, "--output,-o", app->imagename, "Image output");
  add_option(cli, "--lens/--no-lens", app->add_lens, "Add lens");
  add_option(cli, "scene", app->filename, "Scene filename", true);

  parse_cli(cli, argc, argv);

  // scene loading
  auto ioscene_guard = std::make_unique<sio::model>();
  auto ioscene       = ioscene_guard.get();
  auto ioerror       = ""s;
  if (!load_scene(app->filename, ioscene, ioerror, cli::print_progress))
    cli::print_fatal(ioerror);

  // get camera
  auto iocamera = get_camera(ioscene, camera_name);

  // conversion
  init_scene(app->scene, ioscene, app->camera, iocamera, cli::print_progress);

  // camera names
  init_camera_names(app->camera_names, ioscene->cameras);

  // cleanup
  ioscene_guard.reset();

  // init subdivs
  init_subdivs(app->scene, app->params, cli::print_progress);

  // build bvh
  init_bvh(app->scene, app->params, cli::print_progress);

  // build lights
  init_lights(app->scene, app->params, cli::print_progress);

  // allocate buffers
  reset_display(app);

  // callbacks
  auto callbacks     = gui::ui_callbacks{};
  callbacks.clear_cb = [app](gui::window* win, const gui::input& input) {
    if (app->glimage) clear_image(app->glimage);
  };
  callbacks.draw_cb = [app](gui::window* win, const gui::input& input) {
    if (!is_initialized(app->glimage)) init_image(app->glimage);
    if (!app->render_counter)
      set_image(app->glimage, app->display, false, false);
    app->glparams.window      = input.window_size;
    app->glparams.framebuffer = input.framebuffer_viewport;
    update_imview(app->glparams.center, app->glparams.scale,
        app->display.size(), app->glparams.window, app->glparams.fit);
    draw_image(app->glimage, app->glparams);
    app->render_counter++;
    if (app->render_counter > 10) app->render_counter = 0;
  };
  callbacks.widgets_cb = [app](gui::window* win, const gui::input& input) {
    std::string save_path = "";
    if (draw_filedialog_button(win, "save image",
          true, "save image", save_path, true,
          sfs::path(save_path).parent_path(), sfs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
      app->outname = save_path;
      save_image(app->imagename, app->display, app->error);
      save_path = "";
    }
    auto  edited  = 0;
    auto& tparams = app->params;
    draw_progressbar(win, "render", app->progress);
    edited += draw_combobox(
        win, "camera", app->camera, app->scene->cameras, app->camera_names);
    edited += draw_slider(win, "resolution", tparams.resolution, 180, 4096);
    edited += draw_slider(win, "nsamples", tparams.samples, 16, 4096);
    edited += draw_combobox(
        win, "shader", (int&)tparams.shader, ptr::shader_names);
    edited += draw_slider(win, "nbounces", tparams.bounces, 1, 128);
    edited += draw_slider(win, "pratio", tparams.pratio, 1, 64);
    edited += draw_slider(win, "exposure", app->exposure, -5, 5);    
    if (app->add_lens)
      edited += draw_slider(win, "lens", app->camera->lens, 0.01f, 1);
    if (edited) reset_display(app);
  };
  callbacks.uiupdate_cb = [app](gui::window* win, const gui::input& input) {
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
        pan = (input.mouse_pos - input.mouse_last) * app->camera->focus /
              200.0f;
      pan.x = -pan.x;
      update_turntable(
          app->camera->frame, app->camera->focus, rotate, dolly, pan);
      reset_display(app);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yscnitraces", callbacks);

  // done
  return 0;
}
