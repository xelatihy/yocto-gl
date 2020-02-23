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

#include "../yocto/yocto_cli.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "../yocto_gui/yocto_gui.h"
using namespace ym;

#include <atomic>
#include <deque>
#include <future>
using namespace std::string_literals;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

namespace ysio {
void print_obj_camera(ysio::camera* camera);
};  // namespace ysio

// Application scene
struct app_state {
  // loading options
  std::string filename  = "app->yaml";
  std::string imagename = "out.png";
  std::string outname   = "out.yaml";
  std::string name      = "";

  // scene
  ysio::model*  ioscene  = new ysio::model{};
  ytrc::scene*  scene    = new ytrc::scene{};
  ysio::camera* iocamera = nullptr;
  ytrc::camera* camera   = nullptr;

  // options
  ytrc::trace_params params = {};

  // rendering state
  yimg::image<vec4f> render   = {};
  yimg::image<vec4f> display  = {};
  float              exposure = 0;

  // view scene
  ygui::image*       glimage  = new ygui::image{};
  ygui::image_params glparams = {};

  // editing
  ysio::camera*      selected_camera      = nullptr;
  ysio::object*      selected_object      = nullptr;
  ysio::instance*    selected_instance    = nullptr;
  ysio::shape*       selected_shape       = nullptr;
  ysio::subdiv*      selected_subdiv      = nullptr;
  ysio::material*    selected_material    = nullptr;
  ysio::environment* selected_environment = nullptr;
  ysio::texture*     selected_texture     = nullptr;

  // computation
  int          render_sample  = 0;
  int          render_counter = 0;
  ytrc::state* render_state   = new ytrc::state{};

  // loading status
  std::atomic<bool>  ok           = false;
  std::future<void>  loader       = {};
  std::string        status       = "";
  std::string        error        = "";
  std::atomic<float> progress     = 0.5;
  std::string        loader_error = "";

  ~app_state() {
    if (render_state) {
      ytrc::trace_stop(render_state);
      delete render_state;
    }
    if (scene) delete scene;
    if (ioscene) delete ioscene;
    if (glimage) delete glimage;
  }
};

// Application state
struct app_states {
  // data
  std::vector<app_state*> states   = {};
  app_state*              selected = nullptr;
  std::deque<app_state*>  loading  = {};

  // default options
  ytrc::trace_params params     = {};
  bool               add_skyenv = false;

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

// Construct a scene from io
void init_scene(ytrc::scene* scene, ysio::model* ioscene, ytrc::camera*& camera,
    ysio::camera* iocamera, ysio::progress_callback progress_cb = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->subdivs.size() +
             (int)ioscene->instances.size() + (int)ioscene->objects.size()};

  auto camera_map     = std::unordered_map<ysio::camera*, ytrc::camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb)
      progress_cb("converting cameras", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
    camera_map[iocamera] = camera;
  }

  auto texture_map     = std::unordered_map<ysio::texture*, ytrc::texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb)
      progress_cb("converting textures", progress.x++, progress.y);
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

  auto material_map = std::unordered_map<ysio::material*, ytrc::material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb)
      progress_cb("converting materials", progress.x++, progress.y);
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
    set_normalmap(material, texture_map.at(iomaterial->normal_tex));
    set_scattering(material, iomaterial->scattering, iomaterial->scanisotropy,
        texture_map.at(iomaterial->scattering_tex));
    material_map[iomaterial] = material;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    if (progress_cb)
      progress_cb("converting subdivs", progress.x++, progress.y);
    tesselate_subdiv(ioscene, iosubdiv);
  }

  auto shape_map     = std::unordered_map<ysio::shape*, ytrc::shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("converting shapes", progress.x++, progress.y);
    auto shape = add_shape(scene);
    set_points(shape, ioshape->points);
    set_lines(shape, ioshape->lines);
    set_triangles(shape, ioshape->triangles);
    set_quads(shape, ioshape->quads);
    set_positions(shape, ioshape->positions);
    set_normals(shape, ioshape->normals);
    set_texcoords(shape, ioshape->texcoords);
    set_colors(shape, ioshape->colors);
    set_radius(shape, ioshape->radius);
    set_tangents(shape, ioshape->tangents);
    shape_map[ioshape] = shape;
  }

  auto instance_map = std::unordered_map<ysio::instance*, ytrc::instance*>{};
  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    if (progress_cb)
      progress_cb("converting instances", progress.x++, progress.y);
    auto instance = add_instance(scene);
    set_frames(instance, ioinstance->frames);
    instance_map[ioinstance] = instance;
  }

  for (auto ioobject : ioscene->objects) {
    if (progress_cb)
      progress_cb("converting objects", progress.x++, progress.y);
    auto object = add_object(scene);
    set_frame(object, ioobject->frame);
    set_shape(object, shape_map.at(ioobject->shape));
    set_material(object, material_map.at(ioobject->material));
    set_instance(object, instance_map.at(ioobject->instance));
  }

  for (auto ioenvironment : ioscene->environments) {
    if (progress_cb)
      progress_cb("converting environments", progress.x++, progress.y);
    auto environment = add_environment(scene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // done
  if (progress_cb) progress_cb("converting done", progress.x++, progress.y);

  // get camera
  camera = camera_map.at(iocamera);
}

void stop_display(app_state* app) {
  // stop render
  ytrc::trace_stop(app->render_state);
}

void reset_display(app_state* app) {
  // stop render
  ytrc::trace_stop(app->render_state);

  // start render
  app->status         = "render";
  app->render_counter = 0;
  ytrc::trace_start(
      app->render_state, app->scene, app->camera, app->params,
      [app](const std::string& message, int sample, int nsamples) {
        app->progress = (float)sample / (float)nsamples;
      },
      [app](const yimg::image<vec4f>& render, int current, int total) {
        if (current > 0) return;
        app->render  = render;
        app->display = tonemap_image(app->render, app->exposure);
      },
      [app](const yimg::image<vec4f>& render, int current, int total,
          const vec2i& ij) {
        app->render[ij]  = render[ij];
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
}

void load_scene_async(app_states* apps, const std::string& filename,
    const std::string& camera_name = "") {
  auto app       = apps->states.emplace_back(new app_state{});
  app->name      = fs::path(filename).filename().string() + " [loading]";
  app->filename  = filename;
  app->imagename = fs::path(filename).replace_extension(".png");
  app->outname   = fs::path(filename).replace_extension(".edited.yaml");
  app->params    = apps->params;
  app->status    = "load";
  app->loader    = std::async(std::launch::async, [app, camera_name]() {
    auto progress_cb = [app](
                           const std::string& message, int current, int total) {
      app->progress = (float)current / (float)total;
    };
    if (!load_scene(
            app->filename, app->ioscene, app->loader_error, progress_cb))
      return;
    app->progress = 1;
    app->iocamera = get_camera(app->ioscene, camera_name);
    init_scene(
        app->scene, app->ioscene, app->camera, app->iocamera, progress_cb);
    init_bvh(app->scene, app->params);
    init_lights(app->scene);
    if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
      app->params.sampler = ytrc::sampler_type::eyelight;
    }
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = app;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::camera* iocamera) {
  if (!iocamera) return false;
  auto edited = 0;
  draw_label(win, "name", iocamera->name);
  edited += draw_slider(win, "frame.x", iocamera->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", iocamera->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", iocamera->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", iocamera->frame.o, -10, 10);
  edited += draw_checkbox(win, "ortho", iocamera->orthographic);
  edited += draw_slider(win, "lens", iocamera->lens, 0.01f, 1);
  edited += draw_slider(win, "film", iocamera->film, 0.01f, 0.1f);
  edited += draw_slider(win, "focus", iocamera->focus, 0.01f, 1000);
  edited += draw_slider(win, "aperture", iocamera->aperture, 0, 5);
  auto from         = iocamera->frame.o,
       to           = iocamera->frame.o - iocamera->focus * iocamera->frame.z;
  auto from_changed = draw_slider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_slider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    iocamera->frame = lookat_frame(from, to, {0, 1, 0});
    iocamera->focus = length(from - to);
    edited += 1;
  }
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::texture* iotexture) {
  if (!iotexture) return false;
  auto edited = 0;
  draw_label(win, "name", iotexture->name);
  draw_label(win, "colorf",
      std::to_string(iotexture->colorf.size().x) + " x " +
          std::to_string(iotexture->colorf.size().y));
  draw_label(win, "colorb",
      std::to_string(iotexture->colorb.size().x) + " x " +
          std::to_string(iotexture->colorb.size().y));
  draw_label(win, "scalarf",
      std::to_string(iotexture->scalarf.size().x) + " x " +
          std::to_string(iotexture->scalarf.size().y));
  draw_label(win, "scalarb",
      std::to_string(iotexture->scalarb.size().x) + " x " +
          std::to_string(iotexture->scalarb.size().y));
  // TODO: load texture
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::material* iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
  draw_label(win, "name", iomaterial->name);
  edited += draw_hdrcoloredit(win, "emission", iomaterial->emission);
  edited += draw_coloredit(win, "color", iomaterial->color);
  edited += draw_slider(win, "opacity", iomaterial->opacity, 0, 1);
  edited += draw_slider(win, "metallic", iomaterial->metallic, 0, 1);
  edited += draw_slider(win, "roughness", iomaterial->roughness, 0, 1);
  edited += draw_slider(win, "specular", iomaterial->specular, 0, 1);
  edited += draw_slider(win, "coat", iomaterial->coat, 0, 1);
  edited += draw_slider(win, "transmission", iomaterial->transmission, 0, 1);
  edited += draw_coloredit(win, "spectint", iomaterial->spectint);
  edited += draw_checkbox(win, "thin", iomaterial->thin);
  edited += draw_coloredit(win, "scattering", iomaterial->scattering);
  edited += draw_slider(win, "trdepth", iomaterial->trdepth, 0, 1);
  edited += draw_slider(win, "scanisotropy", iomaterial->scanisotropy, -1, 1);
  edited += draw_combobox(
      win, "emission_tex", iomaterial->emission_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "color_tex", iomaterial->color_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "opacity_tex", iomaterial->opacity_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "metallic_tex", iomaterial->metallic_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "roughness_tex", iomaterial->roughness_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "specular_tex", iomaterial->specular_tex, ioscene->textures, true);
  edited += draw_combobox(win, "transmission_tex", iomaterial->transmission_tex,
      ioscene->textures, true);
  edited += draw_combobox(win, "scattering_tex", iomaterial->scattering_tex,
      ioscene->textures, true);
  edited += draw_combobox(
      win, "spectint_tex", iomaterial->spectint_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "normal_tex", iomaterial->normal_tex, ioscene->textures, true);
  edited += draw_combobox(win, "displacement_tex", iomaterial->displacement_tex,
      ioscene->textures, true);
  edited += draw_slider(win, "subdivisions", iomaterial->subdivisions, 0, 5);
  edited += draw_checkbox(win, "smooth", iomaterial->smooth);
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::shape* ioshape) {
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
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::instance* ioinstance) {
  if (!ioinstance) return false;
  auto edited = 0;
  draw_label(win, "name", ioinstance->name);
  draw_label(win, "frames", std::to_string(ioinstance->frames.size()));
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::object* ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  draw_label(win, "name", ioobject->name);
  edited += draw_slider(win, "frame.x", ioobject->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", ioobject->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", ioobject->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_combobox(win, "shape", ioobject->shape, ioscene->shapes);
  edited += draw_combobox(
      win, "material", ioobject->material, ioscene->materials);
  edited += draw_combobox(
      win, "instance", ioobject->instance, ioscene->instances);
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::subdiv* iosubdiv) {
  if (!iosubdiv) return false;
  auto edited = 0;
  draw_label(win, "name", iosubdiv->name);
  draw_label(win, "quads pos", std::to_string(iosubdiv->quadspos.size()));
  draw_label(win, "quads norm", std::to_string(iosubdiv->quadsnorm.size()));
  draw_label(
      win, "quads texcoord", std::to_string(iosubdiv->quadstexcoord.size()));
  draw_label(win, "pos", std::to_string(iosubdiv->positions.size()));
  draw_label(win, "norm", std::to_string(iosubdiv->normals.size()));
  draw_label(win, "texcoord", std::to_string(iosubdiv->texcoords.size()));
  return edited;
}

bool draw_widgets(
    ygui::window* win, ysio::model* ioscene, ysio::environment* ioenvironment) {
  if (!ioenvironment) return false;
  auto edited = 0;
  draw_label(win, "name", ioenvironment->name);
  edited += draw_slider(win, "frame.x", ioenvironment->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", ioenvironment->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", ioenvironment->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_hdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_combobox(win, "emission texture", ioenvironment->emission_tex,
      ioscene->textures, true);
  return edited;
}

template <typename T, typename T1>
T1* get_element(T* ioelement, const std::vector<T*>& ioelements,
    const std::vector<T1*>& elements) {
  if (!ioelement) return nullptr;
  for (auto pos = 0; pos < ioelements.size(); pos++) {
    if (ioelements[pos] == ioelement) return elements[pos];
  }
  throw std::runtime_error("element not found");
}

void draw_widgets(
    ygui::window* win, app_states* apps, const ygui::input& input) {
  static std::string load_path = "", save_path = "", error_message = "";
  if (draw_filedialog_button(win, "load", true, "load", load_path, false, "./",
          "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", apps->selected && apps->selected->ok,
          "save", save_path, true, fs::path(save_path).parent_path(),
          fs::path(save_path).filename(), "*.yaml;*.obj;*.pbrt")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_scene(app->outname, app->ioscene, app->error);
    save_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save image",
          apps->selected && apps->selected->ok, "save image", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_image(app->imagename, app->display, app->error);
    save_path = "";
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
  draw_combobox(win, "scene", apps->selected, apps->states, false);
  if (!apps->selected) return;
  draw_progressbar(
      win, apps->selected->status.c_str(), apps->selected->progress);
  if (apps->selected->error != "") {
    draw_label(win, "error", apps->selected->error);
    return;
  }
  if (!apps->selected->ok) return;
  auto app = apps->selected;
  if (begin_header(win, "trace")) {
    auto edited = 0;
    if (draw_combobox(win, "camera", app->iocamera, app->ioscene->cameras)) {
      app->camera = get_element(
          app->iocamera, app->ioscene->cameras, app->scene->cameras);
      edited += 1;
    }
    auto& tparams = app->params;
    edited += draw_slider(win, "resolution", tparams.resolution, 180, 4096);
    edited += draw_slider(win, "nsamples", tparams.samples, 16, 4096);
    edited += draw_combobox(
        win, "tracer", (int&)tparams.sampler, ytrc::sampler_names);
    edited += draw_combobox(
        win, "false color", (int&)tparams.falsecolor, ytrc::falsecolor_names);
    edited += draw_slider(win, "nbounces", tparams.bounces, 1, 128);
    edited += draw_checkbox(win, "envhidden", tparams.envhidden);
    continue_line(win);
    edited += draw_checkbox(win, "filter", tparams.tentfilter);
    edited += draw_slider(win, "seed", (int&)tparams.seed, 0, 1000000);
    edited += draw_slider(win, "pratio", tparams.pratio, 1, 64);
    edited += draw_slider(win, "exposure", app->exposure, -5, 5);
    if (edited) reset_display(app);
    end_header(win);
  }
  if (begin_header(win, "inspect")) {
    draw_label(win, "scene", fs::path(app->filename).filename());
    draw_label(win, "filename", app->filename);
    draw_label(win, "outname", app->outname);
    draw_label(win, "imagename", app->imagename);
    if (app->ok) {
      draw_label(win, "image",
          std::to_string(app->render.size().x) + " x " +
              std::to_string(app->render.size().y) + " @ " +
              std::to_string(app->render_sample));
      draw_slider(win, "zoom", app->glparams.scale, 0.1, 10);
      draw_checkbox(win, "zoom to fit", app->glparams.fit);
      continue_line(win);
      if (draw_button(win, "print cams")) {
        for (auto iocamera : app->ioscene->cameras) {
          print_obj_camera(iocamera);
        }
      }
      continue_line(win);
      if (draw_button(win, "print stats")) {
        for (auto stat : scene_stats(app->ioscene)) ycli::print_info(stat);
      }
      auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->render.size());
      draw_dragger(win, "mouse", ij);
      if (ij.x >= 0 && ij.x < app->render.size().x && ij.y >= 0 &&
          ij.y < app->render.size().y) {
        draw_coloredit(win, "pixel", app->render[{ij.x, ij.y}]);
      } else {
        auto zero4f_ = zero4f;
        draw_coloredit(win, "pixel", zero4f_);
      }
    }
    end_header(win);
  }
  auto get_texture = [app](ysio::texture* iotexture) {
    return get_element(iotexture, app->ioscene->textures, app->scene->textures);
  };
  if (!app->ioscene->cameras.empty() && begin_header(win, "cameras")) {
    draw_combobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras, true);
    if (draw_widgets(win, app->ioscene, app->selected_camera)) {
      stop_display(app);
      auto iocamera = app->selected_camera;
      auto camera   = get_element(
          iocamera, app->ioscene->cameras, app->scene->cameras);
      set_frame(camera, iocamera->frame);
      set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_focus(camera, iocamera->aperture, iocamera->focus);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->environments.empty() &&
      begin_header(win, "environments")) {
    draw_combobox(win, "environment##2", app->selected_environment,
        app->ioscene->environments, true);
    if (draw_widgets(win, app->ioscene, app->selected_environment)) {
      stop_display(app);
      auto ioenvironment = app->selected_environment;
      auto environment   = get_element(
          ioenvironment, app->ioscene->environments, app->scene->environments);
      set_frame(environment, ioenvironment->frame);
      set_emission(environment, ioenvironment->emission,
          get_texture(ioenvironment->emission_tex));
      init_lights(app->scene);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->objects.empty() && begin_header(win, "objects")) {
    draw_combobox(
        win, "object##2", app->selected_object, app->ioscene->objects, true);
    if (draw_widgets(win, app->ioscene, app->selected_object)) {
      stop_display(app);
      auto ioobject = app->selected_object;
      auto object   = get_element(
          ioobject, app->ioscene->objects, app->scene->objects);
      set_frame(object, ioobject->frame);
      set_shape(object, get_element(ioobject->shape, app->ioscene->shapes,
                            app->scene->shapes));
      set_material(object, get_element(ioobject->material,
                               app->ioscene->materials, app->scene->materials));
      set_instance(object, get_element(ioobject->instance,
                               app->ioscene->instances, app->scene->instances));
      update_bvh(app->scene, {object}, {}, {}, app->params);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->shapes.empty() && begin_header(win, "shapes")) {
    draw_combobox(
        win, "shape##2", app->selected_shape, app->ioscene->shapes, true);
    if (draw_widgets(win, app->ioscene, app->selected_shape)) {
      stop_display(app);
      auto ioshape = app->selected_shape;
      auto shape   = get_element(
          ioshape, app->ioscene->shapes, app->scene->shapes);
      set_points(shape, ioshape->points);
      set_lines(shape, ioshape->lines);
      set_triangles(shape, ioshape->triangles);
      set_quads(shape, ioshape->quads);
      set_positions(shape, ioshape->positions);
      set_normals(shape, ioshape->normals);
      set_texcoords(shape, ioshape->texcoords);
      set_colors(shape, ioshape->colors);
      set_radius(shape, ioshape->radius);
      set_tangents(shape, ioshape->tangents);
      update_bvh(app->scene, {}, {shape}, {}, app->params);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->materials.empty() && begin_header(win, "materials")) {
    draw_combobox(win, "material##2", app->selected_material,
        app->ioscene->materials, true);
    if (draw_widgets(win, app->ioscene, app->selected_material)) {
      stop_display(app);
      auto iomaterial = app->selected_material;
      auto material   = get_element(
          iomaterial, app->ioscene->materials, app->scene->materials);
      set_emission(material, iomaterial->emission,
          get_texture(iomaterial->emission_tex));
      set_color(
          material, iomaterial->color, get_texture(iomaterial->color_tex));
      set_specular(material, iomaterial->specular,
          get_texture(iomaterial->specular_tex));
      set_ior(material, iomaterial->ior);
      set_metallic(material, iomaterial->metallic,
          get_texture(iomaterial->metallic_tex));
      set_transmission(material, iomaterial->transmission, iomaterial->thin,
          iomaterial->trdepth, get_texture(iomaterial->transmission_tex));
      set_roughness(material, iomaterial->roughness,
          get_texture(iomaterial->roughness_tex));
      set_opacity(
          material, iomaterial->opacity, get_texture(iomaterial->opacity_tex));
      set_thin(material, iomaterial->thin);
      set_normalmap(material, get_texture(iomaterial->normal_tex));
      set_scattering(material, iomaterial->scattering, iomaterial->scanisotropy,
          get_texture(iomaterial->scattering_tex));
      init_lights(app->scene);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->textures.empty() && begin_header(win, "textures")) {
    draw_combobox(win, "textures##2", app->selected_texture,
        app->ioscene->textures, true);
    if (draw_widgets(win, app->ioscene, app->selected_texture)) {
      stop_display(app);
      auto iotexture = app->selected_texture;
      auto texture   = get_element(
          iotexture, app->ioscene->textures, app->scene->textures);
      if (!iotexture->colorf.empty()) {
        set_texture(texture, iotexture->colorf);
      } else if (!iotexture->colorb.empty()) {
        set_texture(texture, iotexture->colorb);
      } else if (!iotexture->scalarf.empty()) {
        set_texture(texture, iotexture->scalarf);
      } else if (!iotexture->scalarb.empty()) {
        set_texture(texture, iotexture->scalarb);
      }
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->instances.empty() && begin_header(win, "instances")) {
    draw_combobox(win, "instance##2", app->selected_instance,
        app->ioscene->instances, true);
    if (draw_widgets(win, app->ioscene, app->selected_instance)) {
      stop_display(app);
      auto ioinstance = app->selected_instance;
      auto instance   = get_element(
          ioinstance, app->ioscene->instances, app->scene->instances);
      set_frames(instance, ioinstance->frames);
      update_bvh(app->scene, {}, {}, {instance}, app->params);
      reset_display(app);
    }
    end_header(win);
  }
  if (!app->ioscene->subdivs.empty() && begin_header(win, "subdivs")) {
    draw_combobox(
        win, "selection##2", app->selected_subdiv, app->ioscene->subdivs, true);
    if (draw_widgets(win, app->ioscene, app->selected_subdiv)) {
      stop_display(app);
      // TODO: subdiv not implemented yet
      reset_display(app);
    }
    end_header(win);
  }
}

void draw(ygui::window* win, app_states* apps, const ygui::input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app                  = apps->selected;
  app->glparams.window      = input.window_size;
  app->glparams.framebuffer = input.framebuffer_viewport;
  if (!is_initialized(app->glimage)) init_glimage(app->glimage);
  if (!app->render_counter)
    set_glimage(app->glimage, app->display, false, false);
  update_imview(app->glparams.center, app->glparams.scale, app->display.size(),
      app->glparams.window, app->glparams.fit);
  draw_glimage(app->glimage, app->glparams);
  app->render_counter++;
  if (app->render_counter > 10) app->render_counter = 0;
}

void update(ygui::window* win, app_states* apps) {
  auto is_ready = [](const std::future<void>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    app->loader.get();
    if (app->loader_error.empty()) {
      app->status = "done";
      app->ok     = true;
      reset_display(app);
    } else {
      app->error  = app->loader_error;
      app->status = "error";
    }
  }
}

int main(int argc, const char* argv[]) {
  // application
  auto apps_guard  = std::make_unique<app_states>();
  auto apps        = apps_guard.get();
  auto filenames   = std::vector<std::string>{};
  auto camera_name = ""s;

  // parse command line
  auto cli = ycli::make_cli("yscnitrace", "progressive path tracing");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(
      cli, "--resolution,-r", apps->params.resolution, "Image resolution.");
  add_option(cli, "--samples,-s", apps->params.samples, "Number of samples.");
  add_option(cli, "--tracer,-t", apps->params.sampler, "Tracer type.",
      ytrc::sampler_names);
  add_option(cli, "--falsecolor,-F", apps->params.falsecolor,
      "Tracer false color type.", ytrc::falsecolor_names);
  add_option(
      cli, "--bounces", apps->params.bounces, "Maximum number of bounces.");
  add_option(cli, "--clamp", apps->params.clamp, "Final pixel clamping.");
  add_option(
      cli, "--filter/--no-filter", apps->params.tentfilter, "Filter image.");
  add_option(cli, "--env-hidden/--no-env-hidden", apps->params.envhidden,
      "Environments are hidden in renderer");
  add_option(cli, "--bvh", apps->params.bvh, "Bvh type", ytrc::bvh_names);
  add_option(cli, "scenes", filenames, "Scene filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename, camera_name);

  // callbacks
  auto callbacks    = ygui::ui_callbacks{};
  callbacks.draw_cb = [apps](ygui::window* win, const ygui::input& input) {
    draw(win, apps, input);
  };
  callbacks.widgets_cb = [apps](ygui::window* win, const ygui::input& input) {
    draw_widgets(win, apps, input);
  };
  callbacks.drop_cb = [apps](ygui::window*                win,
                          const std::vector<std::string>& paths,
                          const ygui::input&              input) {
    for (auto& path : paths) load_scene_async(apps, path);
  };
  callbacks.update_cb = [apps](ygui::window* win, const ygui::input& input) {
    update(win, apps);
  };
  callbacks.uiupdate_cb = [apps](ygui::window* win, const ygui::input& input) {
    if (!apps->selected) return;
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
        pan = (input.mouse_pos - input.mouse_last) * app->iocamera->focus /
              200.0f;
      pan.x = -pan.x;
      stop_display(app);
      update_turntable(
          app->iocamera->frame, app->iocamera->focus, rotate, dolly, pan);
      set_frame(app->camera, app->iocamera->frame);
      set_lens(app->camera, app->iocamera->lens, app->iocamera->aspect,
          app->iocamera->film);
      set_focus(app->camera, app->iocamera->aperture, app->iocamera->focus);
      reset_display(app);
    }

    // selection
    if ((input.mouse_left || input.mouse_right) && input.modifier_alt &&
        !input.widgets_active) {
      auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->render.size());
      if (ij.x >= 0 && ij.x < app->render.size().x && ij.y >= 0 &&
          ij.y < app->render.size().y) {
        auto ray = camera_ray(app->camera->frame, app->camera->lens,
            app->camera->film,
            vec2f{ij.x + 0.5f, ij.y + 0.5f} / vec2f{(float)app->render.size().x,
                                                  (float)app->render.size().y});
        if (auto isec = intersect_scene_bvh(app->scene, ray); isec.hit) {
          app->selected_object = app->ioscene->objects[isec.object];
        }
      }
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yscnitrace", callbacks);

  // done
  return 0;
}
