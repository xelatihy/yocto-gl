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

#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <deque>
#include <future>
using namespace std;

#include "ext/CLI11.hpp"
#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

namespace yocto {
void print_obj_camera(sceneio_camera* camera);
};  // namespace yocto

// Application scene
struct app_state {
  // loading options
  string filename  = "app->yaml";
  string imagename = "out.png";
  string outname   = "out.yaml";
  string name      = "";

  // options
  trace_params params = {};
  int          pratio = 8;

  // scene
  sceneio_model* ioscene = new sceneio_model{};
  trace_scene*   scene   = new trace_scene{};

  // rendering state
  image<vec4f> render   = {};
  image<vec4f> display  = {};
  float        exposure = 0;

  // view scene
  opengl_image*       glimage  = new opengl_image{};
  draw_glimage_params glparams = {};

  // editing
  sceneio_camera*      selected_camera      = nullptr;
  sceneio_object*      selected_object      = nullptr;
  sceneio_instance*    selected_instance    = nullptr;
  sceneio_shape*       selected_shape       = nullptr;
  sceneio_subdiv*      selected_subdiv      = nullptr;
  sceneio_material*    selected_material    = nullptr;
  sceneio_environment* selected_environment = nullptr;
  sceneio_texture*     selected_texture     = nullptr;

  // computation
  int          render_sample  = 0;
  int          render_counter = 0;
  trace_state* render_state   = new trace_state{};

  // loading status
  atomic<bool>       ok       = false;
  future<void>       loader   = {};
  string             error    = "";
  std::atomic<float> progress = 0.5;

  ~app_state() {
    if (render_state) {
      trace_async_stop(render_state);
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
  vector<app_state*> states   = {};
  app_state*         selected = nullptr;
  deque<app_state*>  loading  = {};

  // default options
  trace_params params     = {};
  bool         add_skyenv = false;

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

// Construct a scene from io
void init_scene(trace_scene* scene, sceneio_model* ioscene,
    sceneio_progress progress_cb = {}) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->environments.size() +
             (int)ioscene->materials.size() + (int)ioscene->textures.size() +
             (int)ioscene->shapes.size() + (int)ioscene->subdivs.size() +
             (int)ioscene->instances.size() + (int)ioscene->objects.size()};

  for (auto iocamera : ioscene->cameras) {
    if (progress_cb)
      progress_cb("converting cameras", progress.x++, progress.y);
    auto camera = add_camera(scene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_focus(camera, iocamera->aperture, iocamera->focus);
  }

  auto texture_map     = unordered_map<sceneio_texture*, trace_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb)
      progress_cb("converting textures", progress.x++, progress.y);
    auto texture = add_texture(scene);
    if (!iotexture->hdr.empty()) {
      set_texture(texture, std::move(iotexture->hdr));
    } else if (!iotexture->ldr.empty()) {
      set_texture(texture, std::move(iotexture->ldr));
    }
    texture_map[iotexture] = texture;
  }

  auto material_map     = unordered_map<sceneio_material*, trace_material*>{};
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

  auto shape_map     = unordered_map<sceneio_shape*, trace_shape*>{};
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

  auto instance_map     = unordered_map<sceneio_instance*, trace_instance*>{};
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
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto        futures  = vector<future<void>>{};
  auto        nthreads = thread::hardware_concurrency();
  atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(async(launch::async, [&func, &next_idx, size]() {
      while (true) {
        auto j = next_idx.fetch_add(1);
        if (j >= size.y) break;
        for (auto i = 0; i < size.x; i++) func({i, j});
      }
    }));
  }
  for (auto& f : futures) f.get();
}

void stop_display(app_state* app) {
  // stop render
  trace_async_stop(app->render_state);
  app->render_state = nullptr;
}

void reset_display(app_state* app) {
  // stop render
  trace_async_stop(app->render_state);

  // start render
  app->render_counter = 0;
  trace_async_start(
      app->render_state, app->scene, app->params,
      [app](const string& message, int sample, int nsamples) {
        app->progress = (float)sample / (float)nsamples;
      },
      [app](const image<vec4f>& render, int current, int total) {
        if (current > 0) return;
        app->render  = render;
        app->display = tonemap_image(app->render, app->exposure);
      },
      [app](
          const image<vec4f>& render, int current, int total, const vec2i& ij) {
        app->render[ij]  = render[ij];
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
}

void load_scene_async(app_states* apps, const string& filename) {
  auto app       = apps->states.emplace_back(new app_state{});
  app->name      = fs::path(filename).filename().string() + " [loading]";
  app->filename  = filename;
  app->imagename = fs::path(filename).replace_extension(".png");
  app->outname   = fs::path(filename).replace_extension(".edited.yaml");
  app->params    = app->params;
  app->loader    = async(launch::async, [app]() {
    auto progress_cb = [app](const string& message, int current, int total) {
      app->progress = (float)current / (float)total;
    };
    load_scene(app->filename, app->ioscene, progress_cb);
    app->progress = 1;
    init_scene(app->scene, app->ioscene, progress_cb);
    init_bvh(app->scene, app->params);
    init_lights(app->scene);
    if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
      app->params.sampler = trace_sampler_type::eyelight;
    }
    // app->render.resize(app->state->size());
    // app->display.resize(app->state->size());
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = app;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_camera* iocamera) {
  if (!iocamera) return false;
  auto edited = 0;
  draw_gllabel(win, "name", iocamera->name);
  edited += draw_glslider(win, "frame.x", iocamera->frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", iocamera->frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", iocamera->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", iocamera->frame.o, -10, 10);
  edited += draw_glcheckbox(win, "ortho", iocamera->orthographic);
  edited += draw_glslider(win, "lens", iocamera->lens, 0.01f, 1);
  edited += draw_glslider(win, "film", iocamera->film, 0.01f, 0.1f);
  edited += draw_glslider(win, "focus", iocamera->focus, 0.01f, 1000);
  edited += draw_glslider(win, "aperture", iocamera->aperture, 0, 5);
  auto from         = iocamera->frame.o,
       to           = iocamera->frame.o - iocamera->focus * iocamera->frame.z;
  auto from_changed = draw_glslider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_glslider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    iocamera->frame = lookat_frame(from, to, {0, 1, 0});
    iocamera->focus = length(from - to);
    edited += 1;
  }
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_texture* iotexture) {
  if (!iotexture) return false;
  auto edited = 0;
  draw_gllabel(win, "name", iotexture->name);
  draw_gllabel(win, "hdr",
      to_string(iotexture->hdr.size().x) + " x " +
          to_string(iotexture->hdr.size().y));
  draw_gllabel(win, "ldr",
      to_string(iotexture->ldr.size().x) + " x " +
          to_string(iotexture->ldr.size().y));
  // TODO: load texture
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_material* iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
  draw_gllabel(win, "name", iomaterial->name);
  edited += draw_glhdrcoloredit(win, "emission", iomaterial->emission);
  edited += draw_glcoloredit(win, "color", iomaterial->color);
  edited += draw_glslider(win, "specular", iomaterial->specular, 0, 1);
  edited += draw_glslider(win, "metallic", iomaterial->metallic, 0, 1);
  edited += draw_glslider(win, "roughness", iomaterial->roughness, 0, 1);
  edited += draw_glslider(win, "coat", iomaterial->coat, 0, 1);
  edited += draw_glslider(win, "transmission", iomaterial->transmission, 0, 1);
  edited += draw_glcoloredit(win, "spectint", iomaterial->spectint);
  edited += draw_glcheckbox(win, "thin", iomaterial->thin);
  edited += draw_glcoloredit(win, "scattering", iomaterial->scattering);
  edited += draw_glslider(win, "trdepth", iomaterial->trdepth, 0, 1);
  edited += draw_glslider(win, "scanisotropy", iomaterial->scanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", iomaterial->opacity, 0, 1);
  edited += draw_glcombobox(
      win, "emission_tex", iomaterial->emission_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "color_tex", iomaterial->color_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", iomaterial->metallic_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", iomaterial->specular_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "transmission_tex",
      iomaterial->transmission_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "scattering_tex", iomaterial->scattering_tex,
      ioscene->textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", iomaterial->roughness_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "spectint_tex", iomaterial->spectint_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", iomaterial->normal_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "displacement_tex",
      iomaterial->displacement_tex, ioscene->textures, true);
  edited += draw_glcheckbox(win, "glTF textures", iomaterial->gltf_textures);
  edited += draw_glslider(win, "subdivisions", iomaterial->subdivisions, 0, 5);
  edited += draw_glcheckbox(win, "smooth", iomaterial->smooth);
  edited += draw_glcheckbox(win, "glTF textures", iomaterial->gltf_textures);
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_shape* ioshape) {
  if (!ioshape) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioshape->name);
  draw_gllabel(win, "points", to_string(ioshape->points.size()));
  draw_gllabel(win, "lines", to_string(ioshape->lines.size()));
  draw_gllabel(win, "triangles", to_string(ioshape->triangles.size()));
  draw_gllabel(win, "quads", to_string(ioshape->quads.size()));
  draw_gllabel(win, "positions", to_string(ioshape->positions.size()));
  draw_gllabel(win, "normals", to_string(ioshape->normals.size()));
  draw_gllabel(win, "texcoords", to_string(ioshape->texcoords.size()));
  draw_gllabel(win, "colors", to_string(ioshape->colors.size()));
  draw_gllabel(win, "radius", to_string(ioshape->radius.size()));
  draw_gllabel(win, "tangents", to_string(ioshape->tangents.size()));
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_instance* ioinstance) {
  if (!ioinstance) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioinstance->name);
  draw_gllabel(win, "frames", to_string(ioinstance->frames.size()));
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_object* ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioobject->name);
  edited += draw_glslider(win, "frame.x", ioobject->frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", ioobject->frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", ioobject->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_glcombobox(win, "shape", ioobject->shape, ioscene->shapes);
  edited += draw_glcombobox(
      win, "material", ioobject->material, ioscene->materials);
  edited += draw_glcombobox(
      win, "instance", ioobject->instance, ioscene->instances);
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_subdiv* iosubdiv) {
  if (!iosubdiv) return false;
  auto edited = 0;
  draw_gllabel(win, "name", iosubdiv->name);
  draw_gllabel(win, "quads pos", to_string(iosubdiv->quadspos.size()));
  draw_gllabel(win, "quads norm", to_string(iosubdiv->quadsnorm.size()));
  draw_gllabel(
      win, "quads texcoord", to_string(iosubdiv->quadstexcoord.size()));
  draw_gllabel(win, "pos", to_string(iosubdiv->positions.size()));
  draw_gllabel(win, "norm", to_string(iosubdiv->normals.size()));
  draw_gllabel(win, "texcoord", to_string(iosubdiv->texcoords.size()));
  return edited;
}

bool draw_glwidgets(opengl_window* win, sceneio_model* ioscene,
    sceneio_environment* ioenvironment) {
  if (!ioenvironment) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioenvironment->name);
  edited += draw_glslider(win, "frame.x", ioenvironment->frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", ioenvironment->frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", ioenvironment->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_glcombobox(win, "emission texture",
      ioenvironment->emission_tex, ioscene->textures, true);
  return edited;
}

template <typename T, typename T1>
T1* get_element(
    T* ioelement, const vector<T*>& ioelements, const vector<T1*>& elements) {
  if (!ioelement) return nullptr;
  for (auto pos = 0; pos < ioelements.size(); pos++) {
    if (ioelements[pos] == ioelement) return elements[pos];
  }
  throw std::runtime_error("element not found");
}

void draw_glwidgets(
    opengl_window* win, app_states* apps, const opengl_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  if (draw_glfiledialog_button(win, "load", true, "load", load_path, false,
          "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save",
          apps->selected && apps->selected->ok, "save", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.yaml;*.obj;*.pbrt")) {
    auto app     = apps->selected;
    app->outname = save_path;
    try {
      save_scene(app->outname, app->ioscene);
    } catch (std::exception& e) {
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save image",
          apps->selected && apps->selected->ok, "save image", save_path, true,
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->selected;
    app->outname = save_path;
    try {
      save_image(app->imagename, app->display);
    } catch (std::exception& e) {
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) return;
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_close(win, true);
  }
  draw_glcombobox(win, "scene", apps->selected, apps->states, false);
  if (!apps->selected) return;
  if (apps->selected->error != "") {
    draw_gllabel(win, "error", apps->selected->error);
    return;
  }
  if (apps->selected->progress) {
    draw_glprogressbar(win, apps->selected->progress);
  }
  if (!apps->selected->ok) return;
  auto app = apps->selected;
  if (begin_glheader(win, "trace")) {
    auto  edited  = 0;
    auto& tparams = app->params;
    edited += draw_glcombobox(
        win, "camera", tparams.camera, app->ioscene->cameras);
    edited += draw_glslider(win, "resolution", tparams.resolution, 180, 4096);
    edited += draw_glslider(win, "nsamples", tparams.samples, 16, 4096);
    edited += draw_glcombobox(
        win, "tracer", (int&)tparams.sampler, trace_sampler_names);
    edited += draw_glcombobox(
        win, "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
    edited += draw_glslider(win, "nbounces", tparams.bounces, 1, 128);
    edited += draw_glcheckbox(win, "envhidden", tparams.envhidden);
    continue_glline(win);
    edited += draw_glcheckbox(win, "filter", tparams.tentfilter);
    edited += draw_glslider(win, "seed", (int&)tparams.seed, 0, 1000000);
    edited += draw_glslider(win, "pratio", app->pratio, 1, 64);
    edited += draw_glslider(win, "exposure", app->exposure, -5, 5);
    if (edited) reset_display(app);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    draw_gllabel(win, "scene", fs::path(app->filename).filename());
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "imagename", app->imagename);
    if (app->ok) {
      draw_gllabel(win, "image",
          to_string(app->render.size().x) + " x " +
              to_string(app->render.size().y) + " @ " +
              to_string(app->render_sample));
      draw_glslider(win, "zoom", app->glparams.scale, 0.1, 10);
      draw_glcheckbox(win, "zoom to fit", app->glparams.fit);
      continue_glline(win);
      if (draw_glbutton(win, "print cams")) {
        for (auto iocamera : app->ioscene->cameras) {
          print_obj_camera(iocamera);
        }
      }
      continue_glline(win);
      if (draw_glbutton(win, "print stats")) {
        for (auto stat : scene_stats(app->ioscene))
          printf("%s\n", stat.c_str());
      }
      auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->render.size());
      draw_gldragger(win, "mouse", ij);
      if (ij.x >= 0 && ij.x < app->render.size().x && ij.y >= 0 &&
          ij.y < app->render.size().y) {
        draw_glcoloredit(win, "pixel", app->render[{ij.x, ij.y}]);
      } else {
        auto zero4f_ = zero4f;
        draw_glcoloredit(win, "pixel", zero4f_);
      }
    }
    end_glheader(win);
  }
  auto get_texture = [app](sceneio_texture* iotexture) {
    return get_element(iotexture, app->ioscene->textures, app->scene->textures);
  };
  if (!app->ioscene->cameras.empty() && begin_glheader(win, "cameras")) {
    draw_glcombobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_camera)) {
      stop_display(app);
      auto iocamera = app->selected_camera;
      auto camera   = get_element(
          iocamera, app->ioscene->cameras, app->scene->cameras);
      set_frame(camera, iocamera->frame);
      set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_focus(camera, iocamera->aperture, iocamera->focus);
      reset_display(app);
    }
    end_glheader(win);
  }
  if (!app->ioscene->environments.empty() &&
      begin_glheader(win, "environments")) {
    draw_glcombobox(win, "environment##2", app->selected_environment,
        app->ioscene->environments, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_environment)) {
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
    end_glheader(win);
  }
  if (!app->ioscene->objects.empty() && begin_glheader(win, "objects")) {
    draw_glcombobox(
        win, "object##2", app->selected_object, app->ioscene->objects, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_object)) {
      stop_display(app);
      // auto ioobject = app->ioscene->shapes[app->selected_shape];
      // TODO: update editing
      // update_bvh(app->scene, {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (!app->ioscene->shapes.empty() && begin_glheader(win, "shapes")) {
    draw_glcombobox(
        win, "shape##2", app->selected_shape, app->ioscene->shapes, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_shape)) {
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
      // TODO: bvh update
      // update_bvh(app->scene, {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (!app->ioscene->materials.empty() && begin_glheader(win, "materials")) {
    draw_glcombobox(win, "material##2", app->selected_material,
        app->ioscene->materials, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_material)) {
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
    end_glheader(win);
  }
  if (!app->ioscene->textures.empty() && begin_glheader(win, "textures")) {
    draw_glcombobox(win, "textures##2", app->selected_texture,
        app->ioscene->textures, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_texture)) {
      stop_display(app);
      auto iotexture = app->selected_texture;
      auto texture   = get_element(
          iotexture, app->ioscene->textures, app->scene->textures);
      if (!iotexture->hdr.empty()) {
        set_texture(texture, iotexture->hdr);
      } else if (!iotexture->ldr.empty()) {
        set_texture(texture, iotexture->ldr);
      }
      reset_display(app);
    }
    end_glheader(win);
  }
  if (!app->ioscene->instances.empty() && begin_glheader(win, "instances")) {
    draw_glcombobox(win, "instance##2", app->selected_instance,
        app->ioscene->instances, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_instance)) {
      stop_display(app);
      // auto ioinstance = app->ioscene->instances[app->selected_instance];
      // TODO: update editing
      // update_bvh(app->scene, {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (!app->ioscene->subdivs.empty() && begin_glheader(win, "subdivs")) {
    draw_glcombobox(
        win, "selection##2", app->selected_subdiv, app->ioscene->subdivs, true);
    if (draw_glwidgets(win, app->ioscene, app->selected_subdiv)) {
      stop_display(app);
      // TODO: this is bogus
      // auto iosubdiv = app->ioscene->subdivs[app->selected_subdiv];
      // tesselate_subdiv(app->ioscene, iosubdiv);
      // TODO: this is bogus
      // auto shape   = app->scene->shapes[app->selected_subdiv];
      // auto ioshape = iosubdiv->shape;
      // set_points(shape, ioshape->points);
      // set_lines(shape, ioshape->lines);
      // set_triangles(shape, ioshape->triangles);
      // set_quads(shape, ioshape->quads);
      // set_positions(shape, ioshape->positions);
      // set_normals(shape, ioshape->normals);
      // set_texcoords(shape, ioshape->texcoords);
      // set_colors(shape, ioshape->colors);
      // set_radius(shape, ioshape->radius);
      // set_tangents(shape, ioshape->tangents);
      // TODO: bvh update
      // update_bvh(app->scene, {}, {app->selected_subdiv}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
}

void draw(opengl_window* win, app_states* apps, const opengl_input& input) {
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

void update(opengl_window* win, app_states* apps) {
  auto is_ready = [](const future<void>& result) -> bool {
    return result.valid() &&
           result.wait_for(chrono::microseconds(0)) == future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    try {
      app->loader.get();
      reset_display(app);
      app->name = fs::path(app->filename).filename().string() + " [ok]";
      app->ok   = true;
    } catch (std::exception& e) {
      app->name  = fs::path(app->filename).filename().string() + " [error]";
      app->error = e.what();
    }
  }
}

int run_app(int argc, const char* argv[]) {
  // application
  auto apps_guard = make_unique<app_states>();
  auto apps       = apps_guard.get();
  auto filenames  = vector<string>{};

  // maps for getting param
  auto trace_sampler_map = map<string, trace_sampler_type>{};
  for (auto idx = 0; idx < trace_sampler_names.size(); idx++) {
    trace_sampler_map[trace_sampler_names[idx]] = (trace_sampler_type)idx;
  }
  auto trace_falsecolor_map = map<string, trace_falsecolor_type>{};
  for (auto idx = 0; idx < trace_falsecolor_names.size(); idx++) {
    trace_falsecolor_map[trace_falsecolor_names[idx]] =
        (trace_falsecolor_type)idx;
  }
  auto trace_bvh_map = map<string, trace_bvh_type>{};
  for (auto idx = 0; idx < trace_bvh_names.size(); idx++) {
    trace_bvh_map[trace_bvh_names[idx]] = (trace_bvh_type)idx;
  }

  // parse command line
  auto cli = CLI::App{"progressive path tracing"};
  cli.add_option("--camera", apps->params.camera, "Camera index.");
  cli.add_option(
      "--resolution,-r", apps->params.resolution, "Image resolution.");
  cli.add_option("--samples,-s", apps->params.samples, "Number of samples.");
  cli.add_option("--tracer,-t", apps->params.sampler, "Tracer type.")
      ->transform(CLI::CheckedTransformer(trace_sampler_map));
  cli.add_option(
         "--falsecolor,-F", apps->params.falsecolor, "Tracer false color type.")
      ->transform(CLI::CheckedTransformer(trace_falsecolor_map));
  cli.add_option(
      "--bounces", apps->params.bounces, "Maximum number of bounces.");
  cli.add_option("--clamp", apps->params.clamp, "Final pixel clamping.");
  cli.add_flag("--filter", apps->params.tentfilter, "Filter image.");
  cli.add_flag("--env-hidden,!--no-env-hidden", apps->params.envhidden,
      "Environments are hidden in renderer");
  cli.add_option("--bvh", apps->params.bvh, "Bvh type")
      ->transform(CLI::CheckedTransformer(trace_bvh_map));
  cli.add_option("scenes", filenames, "Scene filenames")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  // window
  auto win_guard = make_glwindow({1280 + 320, 720}, "yscnitrace", true);
  auto win       = win_guard.get();

  // callbacks
  set_draw_glcallback(
      win, [apps](opengl_window* win, const opengl_input& input) {
        draw(win, apps, input);
      });
  set_widgets_glcallback(
      win, [apps](opengl_window* win, const opengl_input& input) {
        draw_glwidgets(win, apps, input);
      });
  set_drop_glcallback(
      win, [apps](opengl_window* win, const vector<string>& paths,
               const opengl_input& input) {
        for (auto& path : paths) load_scene_async(apps, path);
      });
  set_update_glcallback(
      win, [apps](opengl_window* win, const opengl_input& input) {
        update(win, apps);
      });
  set_uiupdate_glcallback(win, [apps](opengl_window*   win,
                                   const opengl_input& input) {
    if (!apps->selected) return;
    auto app = apps->selected;

    // handle mouse and keyboard for navigation
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
      auto iocamera = app->ioscene->cameras.at(app->params.camera);
      auto dolly    = 0.0f;
      auto pan      = zero2f;
      auto rotate   = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * iocamera->focus / 200.0f;
      pan.x = -pan.x;
      stop_display(app);
      update_turntable(iocamera->frame, iocamera->focus, rotate, dolly, pan);
      set_frame(app->scene->cameras[app->params.camera], iocamera->frame);
      set_lens(app->scene->cameras[app->params.camera], iocamera->lens,
          iocamera->aspect, iocamera->film);
      set_focus(app->scene->cameras[app->params.camera], iocamera->aperture,
          iocamera->focus);
      reset_display(app);
    }

    // selection
    if ((input.mouse_left || input.mouse_right) && input.modifier_alt &&
        !input.widgets_active) {
      auto ij = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->render.size());
      if (ij.x >= 0 && ij.x < app->render.size().x && ij.y >= 0 &&
          ij.y < app->render.size().y) {
        auto camera = app->scene->cameras.at(app->params.camera);
        auto ray    = camera_ray(camera->frame, camera->lens, camera->film,
            vec2f{ij.x + 0.5f, ij.y + 0.5f} / vec2f{(float)app->render.size().x,
                                                  (float)app->render.size().y});
        if (auto isec = intersect_scene_bvh(app->scene, ray); isec.hit) {
          app->selected_object = app->ioscene->objects[isec.object];
        }
      }
    }
  });

  // run ui
  run_ui(win);

  // clear
  clear_glwindow(win);

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    fprintf(stderr, "%s\n", e.what());
    return 1;
  }
}
