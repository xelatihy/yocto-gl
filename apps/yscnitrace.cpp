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

#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <deque>
#include <future>
using namespace std;

namespace yocto {
void print_obj_camera(const sceneio_camera* camera);
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
  shared_ptr<sceneio_model> ioscene = nullptr;
  shared_ptr<trace_scene>   scene   = nullptr;

  // rendering state
  shared_ptr<trace_state> state    = nullptr;
  image<vec4f>            render   = {};
  image<vec4f>            display  = {};
  float                   exposure = 0;

  // view scene
  opengl_image        glimage  = {};
  draw_glimage_params glparams = {};

  // editing
  int selected_camera      = -1;
  int selected_object      = -1;
  int selected_instance    = -1;
  int selected_shape       = -1;
  int selected_subdiv      = -1;
  int selected_material    = -1;
  int selected_environment = -1;
  int selected_texture     = -1;

  // editing maps
  unordered_map<sceneio_texture*, trace_texture*>   texture_map  = {};
  unordered_map<sceneio_material*, trace_material*> material_map = {};
  unordered_map<sceneio_shape*, trace_shape*>       shape_map    = {};
  unordered_map<sceneio_instance*, trace_instance*> instance_map = {};

  // computation
  int          render_sample  = 0;
  atomic<bool> render_stop    = {};
  future<void> render_future  = {};
  int          render_counter = 0;

  ~app_state() {
    render_stop = true;
    if (render_future.valid()) render_future.get();
  }
};

// Application state
struct app_states {
  // data
  vector<shared_ptr<app_state>> states   = {};
  int                           selected = -1;

  // loading
  deque<future<shared_ptr<app_state>>> loaders = {};

  // default options
  trace_params params     = {};
  bool         add_skyenv = false;
};

// Construct a scene from io
void init_scene(shared_ptr<app_state> app) {
  auto  scene        = app->scene.get();
  auto  ioscene      = app->ioscene.get();
  auto& texture_map  = app->texture_map;
  auto& material_map = app->material_map;
  auto& shape_map    = app->shape_map;
  auto& instance_map = app->instance_map;

  for (auto iocamera : ioscene->cameras) {
    auto camera = add_camera(scene);
    set_camera_frame(camera, iocamera->frame);
    set_camera_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_camera_focus(camera, iocamera->aperture, iocamera->focus);
  }

  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    auto texture = add_texture(scene);
    if (!iotexture->hdr.empty()) {
      set_texture(texture, std::move(iotexture->hdr));
    } else if (!iotexture->ldr.empty()) {
      set_texture(texture, std::move(iotexture->ldr));
    }
    texture_map[iotexture] = texture;
  }

  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    auto material = add_material(scene);
    set_shape_emission(material, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_shape_color(
        material, iomaterial->color, texture_map.at(iomaterial->color_tex));
    set_shape_specular(material, iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_shape_ior(material, iomaterial->ior);
    set_shape_metallic(material, iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_shape_transmission(material, iomaterial->transmission, iomaterial->thin,
        iomaterial->trdepth, texture_map.at(iomaterial->transmission_tex));
    set_shape_roughness(material, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_shape_opacity(
        material, iomaterial->opacity, texture_map.at(iomaterial->opacity_tex));
    set_shape_thin(material, iomaterial->thin);
    set_shape_normalmap(material, texture_map.at(iomaterial->normal_tex));
    set_shape_scattering(material, iomaterial->scattering,
        iomaterial->scanisotropy, texture_map.at(iomaterial->scattering_tex));
    material_map[iomaterial] = material;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
  }

  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    auto shape = add_shape(scene);
    set_shape_points(shape, ioshape->points);
    set_shape_lines(shape, ioshape->lines);
    set_shape_triangles(shape, ioshape->triangles);
    set_shape_quads(shape, ioshape->quads);
    set_shape_positions(shape, ioshape->positions);
    set_shape_normals(shape, ioshape->normals);
    set_shape_texcoords(shape, ioshape->texcoords);
    set_shape_colors(shape, ioshape->colors);
    set_shape_radius(shape, ioshape->radius);
    set_shape_tangents(shape, ioshape->tangents);
    shape_map[ioshape] = shape;
  }

  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    auto instance = add_instance(scene);
    set_frames(instance, ioinstance->frames);
    instance_map[ioinstance] = instance;
  }

  for (auto ioobject : ioscene->objects) {
    auto object = add_object(scene);
    set_frame(object, ioobject->frame);
    set_shape(object, shape_map.at(ioobject->shape));
    set_material(object, material_map.at(ioobject->material));
    set_instance(object, instance_map.at(ioobject->instance));
  }

  for (auto ioenvironment : ioscene->environments) {
    auto environment = add_environment(scene);
    set_environment_frame(environment, ioenvironment->frame);
    set_environment_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }
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

void stop_display(shared_ptr<app_state> app) {
  // stop render
  app->render_stop = true;
  if (app->render_future.valid()) app->render_future.get();
}

void reset_display(shared_ptr<app_state> app) {
  // stop render
  app->render_stop = true;
  if (app->render_future.valid()) app->render_future.get();

  // reset state
  init_state(app->state.get(), app->scene.get(), app->params);
  app->render.resize(app->state->size());
  app->display.resize(app->state->size());

  // render preview
  auto preview_prms = app->params;
  preview_prms.resolution /= app->pratio;
  preview_prms.samples = 1;
  auto preview         = trace_image(app->scene.get(), preview_prms);
  preview              = tonemap_image(preview, app->exposure);
  for (auto j = 0; j < app->display.size().y; j++) {
    for (auto i = 0; i < app->display.size().x; i++) {
      auto pi              = clamp(i / app->pratio, 0, preview.size().x - 1),
           pj              = clamp(j / app->pratio, 0, preview.size().y - 1);
      app->display[{i, j}] = preview[{pi, pj}];
    }
  }

  // start renderer
  app->render_counter = 0;
  app->render_stop    = false;
  app->render_future  = async(launch::async, [app]() {
    for (auto sample = 0; sample < app->params.samples; sample++) {
      if (app->render_stop) return;
      parallel_for(app->render.size(), [app](const vec2i& ij) {
        if (app->render_stop) return;
        app->render[ij] = trace_sample(
            app->state.get(), app->scene.get(), ij, app->params);
        app->display[ij] = tonemap(app->render[ij], app->exposure);
      });
    }
  });
}

void load_scene_async(shared_ptr<app_states> apps, const string& filename) {
  apps->loaders.push_back(
      async(launch::async, [apps, filename]() -> shared_ptr<app_state> {
        auto app       = make_shared<app_state>();
        app->filename  = filename;
        app->imagename = replace_extension(filename, ".png");
        app->outname   = replace_extension(filename, ".edited.yaml");
        app->name      = get_filename(app->filename);
        app->params    = app->params;
        app->ioscene   = make_shared<sceneio_model>();
        load_scene(app->filename, app->ioscene.get());
        app->scene = make_shared<trace_scene>();
        init_scene(app);
        init_bvh(app->scene.get(), app->params);
        init_lights(app->scene.get());
        if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
          app->params.sampler = trace_sampler_type::eyelight;
        }
        app->state = make_shared<trace_state>();
        init_state(app->state.get(), app->scene.get(), app->params);
        app->render.resize(app->state->size());
        app->display.resize(app->state->size());
        app->name = get_filename(app->filename) + " [" +
                    to_string(app->render.size().x) + "x" +
                    to_string(app->render.size().y) + " @ 0]";
        app->selected_camera      = app->ioscene->cameras.empty() ? -1 : 0;
        app->selected_shape       = app->ioscene->shapes.empty() ? -1 : 0;
        app->selected_subdiv      = app->ioscene->subdivs.empty() ? -1 : 0;
        app->selected_material    = app->ioscene->shapes.empty() ? -1 : 0;
        app->selected_texture     = app->ioscene->textures.empty() ? -1 : 0;
        app->selected_environment = app->ioscene->environments.empty() ? -1 : 0;
        return app;
      }));
}

bool draw_glwidgets_camera(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iocamera = app->ioscene->cameras[id];
  auto edited   = 0;
  edited += (int)draw_gltextinput(win, "name", iocamera->name);
  edited += (int)draw_glslider(win, "frame.x", iocamera->frame.x, -1, 1);
  edited += (int)draw_glslider(win, "frame.y", iocamera->frame.y, -1, 1);
  edited += (int)draw_glslider(win, "frame.z", iocamera->frame.z, -1, 1);
  edited += (int)draw_glslider(win, "frame.o", iocamera->frame.o, -10, 10);
  edited += (int)draw_glcheckbox(win, "ortho", iocamera->orthographic);
  edited += (int)draw_glslider(win, "lens", iocamera->lens, 0.01f, 1);
  edited += (int)draw_glslider(win, "film", iocamera->film, 0.01f, 0.1f);
  edited += (int)draw_glslider(win, "focus", iocamera->focus, 0.01f, 1000);
  edited += (int)draw_glslider(win, "aperture", iocamera->aperture, 0, 5);
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

bool draw_glwidgets_texture(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iotexture = app->ioscene->textures[id];
  auto edited    = 0;
  edited += draw_gltextinput(win, "name", iotexture->name);
  draw_gllabel(win, "hdr",
      to_string(iotexture->hdr.size().x) + " x " +
          to_string(iotexture->hdr.size().y));
  draw_gllabel(win, "ldr",
      to_string(iotexture->ldr.size().x) + " x " +
          to_string(iotexture->ldr.size().y));
  // TODO: load texture
  return edited;
}

bool draw_glwidgets_material(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iomaterial = app->ioscene->materials[id];
  auto edited     = 0;
  edited += draw_gltextinput(win, "name", iomaterial->name);
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
  edited += draw_glcombobox(win, "emission_tex", iomaterial->emission_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(
      win, "color_tex", iomaterial->color_tex, app->ioscene->textures, true);
  edited += draw_glcombobox(win, "metallic_tex", iomaterial->metallic_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(win, "specular_tex", iomaterial->specular_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(win, "transmission_tex",
      iomaterial->transmission_tex, app->ioscene->textures, true);
  edited += draw_glcombobox(win, "scattering_tex", iomaterial->scattering_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(win, "roughness_tex", iomaterial->roughness_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(win, "spectint_tex", iomaterial->spectint_tex,
      app->ioscene->textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", iomaterial->normal_tex, app->ioscene->textures, true);
  edited += draw_glcheckbox(win, "glTF textures", iomaterial->gltf_textures);
  return edited;
}

bool draw_glwidgets_shape(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioshape = app->ioscene->shapes[id];
  auto edited  = 0;
  edited += draw_gltextinput(win, "name", ioshape->name);
  draw_gllabel(win, "points", to_string(ioshape->points.size()));
  draw_gllabel(win, "lines", to_string(ioshape->lines.size()));
  draw_gllabel(win, "triangles", to_string(ioshape->triangles.size()));
  draw_gllabel(win, "quads", to_string(ioshape->quads.size()));
  draw_gllabel(win, "pos", to_string(ioshape->positions.size()));
  draw_gllabel(win, "norm", to_string(ioshape->normals.size()));
  draw_gllabel(win, "texcoord", to_string(ioshape->texcoords.size()));
  draw_gllabel(win, "color", to_string(ioshape->colors.size()));
  draw_gllabel(win, "radius", to_string(ioshape->radius.size()));
  draw_gllabel(win, "tangsp", to_string(ioshape->tangents.size()));
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioinstance = app->ioscene->instances[id];
  auto edited     = 0;
  edited += draw_gltextinput(win, "name", ioinstance->name);
  draw_gllabel(win, "frames", to_string(ioinstance->frames.size()));
  return edited;
}

bool draw_glwidgets_object(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioobject = app->ioscene->objects[id];
  auto edited   = 0;
  edited += draw_gltextinput(win, "name", ioobject->name);
  edited += draw_glslider(win, "frame.x", ioobject->frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", ioobject->frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", ioobject->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", ioobject->shape, app->ioscene->shapes);
  edited += draw_glcombobox(
      win, "material", ioobject->material, app->ioscene->materials);
  edited += draw_glcombobox(
      win, "instance", ioobject->instance, app->ioscene->instances);
  return edited;
}

bool draw_glwidgets_subdiv(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iosubdiv = app->ioscene->subdivs[id];
  auto edited   = 0;
  edited += draw_gltextinput(win, "name", iosubdiv->name);
  draw_gllabel(win, "points", to_string(iosubdiv->points.size()));
  draw_gllabel(win, "lines", to_string(iosubdiv->lines.size()));
  draw_gllabel(win, "triangles", to_string(iosubdiv->triangles.size()));
  draw_gllabel(win, "quads", to_string(iosubdiv->quads.size()));
  draw_gllabel(win, "quads pos", to_string(iosubdiv->quadspos.size()));
  draw_gllabel(win, "quads norm", to_string(iosubdiv->quadsnorm.size()));
  draw_gllabel(
      win, "quads texcoord", to_string(iosubdiv->quadstexcoord.size()));
  draw_gllabel(win, "pos", to_string(iosubdiv->positions.size()));
  draw_gllabel(win, "norm", to_string(iosubdiv->normals.size()));
  draw_gllabel(win, "texcoord", to_string(iosubdiv->texcoords.size()));
  draw_gllabel(win, "color", to_string(iosubdiv->colors.size()));
  draw_gllabel(win, "radius", to_string(iosubdiv->radius.size()));
  draw_gllabel(win, "tangsp", to_string(iosubdiv->tangents.size()));
  // TODO: load
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioenvironment = app->ioscene->environments[id];
  auto edited        = 0;
  edited += draw_gltextinput(win, "name", ioenvironment->name);
  edited += draw_glslider(win, "frame.x", ioenvironment->frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", ioenvironment->frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", ioenvironment->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_glcombobox(win, "emission texture",
      ioenvironment->emission_tex, app->ioscene->textures, true);
  return edited;
}

void draw_glwidgets(const opengl_window& win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  auto          app = (!apps->states.empty() && apps->selected >= 0)
                 ? apps->states[apps->selected]
                 : nullptr;
  if (draw_glfiledialog_button(win, "load", true, "load", load_path, false,
          "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save", (bool)app, "save", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.yaml;*.obj;*.pbrt")) {
    app->outname = save_path;
    try {
      save_scene(app->outname, app->ioscene.get());
    } catch (std::exception& e) {
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save image", (bool)app, "save image",
          save_path, true, get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
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
  if (draw_glbutton(win, "close", (bool)app)) {
    apps->states.erase(apps->states.begin() + apps->selected);
    apps->selected = apps->states.empty() ? -1 : 0;
    app            = (!apps->states.empty() && apps->selected >= 0)
              ? apps->states[apps->selected]
              : nullptr;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_close(win, true);
  }
  draw_glcombobox(
      win, "scene", apps->selected, (int)apps->states.size(),
      [apps](int idx) { return apps->states[apps->selected]->name.c_str(); },
      false);
  if (app && begin_glheader(win, "trace")) {
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
  if (app && begin_glheader(win, "inspect")) {
    draw_gllabel(win, "scene", get_filename(app->filename));
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "imagename", app->imagename);
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
      for (auto stat : scene_stats(app->ioscene.get())) print_info(stat);
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
    end_glheader(win);
  }
  if (app && !app->ioscene->cameras.empty() && begin_glheader(win, "cameras")) {
    draw_glcombobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras);
    if (draw_glwidgets_camera(win, app, app->selected_camera)) {
      stop_display(app);
      auto iocamera = app->ioscene->cameras[app->selected_camera];
      set_camera_frame(
          app->scene->cameras[app->selected_camera], iocamera->frame);
      set_camera_lens(app->scene->cameras[app->selected_camera], iocamera->lens,
          iocamera->aspect, iocamera->film);
      set_camera_focus(app->scene->cameras[app->selected_camera],
          iocamera->aperture, iocamera->focus);
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->environments.empty() &&
      begin_glheader(win, "environments")) {
    draw_glcombobox(win, "environment##2", app->selected_environment,
        app->ioscene->environments);
    if (draw_glwidgets_environment(win, app, app->selected_environment)) {
      stop_display(app);
      auto ioenvironment =
          app->ioscene->environments[app->selected_environment];
      set_environment_frame(app->scene->environments[app->selected_environment],
          ioenvironment->frame);
      set_environment_emission(
          app->scene->environments[app->selected_environment],
          ioenvironment->emission,
          app->texture_map.at(ioenvironment->emission_tex));
      init_lights(app->scene.get());
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->objects.empty() && begin_glheader(win, "objects")) {
    draw_glcombobox(
        win, "object##2", app->selected_object, app->ioscene->objects);
    if (draw_glwidgets_shape(win, app, app->selected_object)) {
      stop_display(app);
      // auto ioobject = app->ioscene->shapes[app->selected_shape];
      // TODO: update editing
      // update_bvh(app->scene, {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->instances.empty() &&
      begin_glheader(win, "instances")) {
    draw_glcombobox(
        win, "instance##2", app->selected_instance, app->ioscene->instances);
    if (draw_glwidgets_shape(win, app, app->selected_instance)) {
      stop_display(app);
      // auto ioinstance = app->ioscene->instances[app->selected_instance];
      // TODO: update editing
      // update_bvh(app->scene, {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->shapes.empty() && begin_glheader(win, "shapes")) {
    draw_glcombobox(win, "shape##2", app->selected_shape, app->ioscene->shapes);
    if (draw_glwidgets_shape(win, app, app->selected_shape)) {
      stop_display(app);
      auto shape   = app->scene->shapes[app->selected_shape];
      auto ioshape = app->ioscene->shapes[app->selected_shape];
      set_shape_points(shape, ioshape->points);
      set_shape_lines(shape, ioshape->lines);
      set_shape_triangles(shape, ioshape->triangles);
      set_shape_quads(shape, ioshape->quads);
      set_shape_positions(shape, ioshape->positions);
      set_shape_normals(shape, ioshape->normals);
      set_shape_texcoords(shape, ioshape->texcoords);
      set_shape_colors(shape, ioshape->colors);
      set_shape_radius(shape, ioshape->radius);
      set_shape_tangents(shape, ioshape->tangents);
      update_bvh(app->scene.get(), {}, {app->selected_shape}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->shapes.empty() &&
      begin_glheader(win, "materials")) {
    draw_glcombobox(
        win, "material##2", app->selected_material, app->ioscene->materials);
    if (draw_glwidgets_material(win, app, app->selected_material)) {
      stop_display(app);
      auto material   = app->scene->materials[app->selected_material];
      auto iomaterial = app->ioscene->materials[app->selected_material];
      set_shape_emission(material, iomaterial->emission,
          app->texture_map.at(iomaterial->emission_tex));
      set_shape_color(material, iomaterial->color,
          app->texture_map.at(iomaterial->color_tex));
      set_shape_specular(material, iomaterial->specular,
          app->texture_map.at(iomaterial->specular_tex));
      set_shape_ior(material, iomaterial->ior);
      set_shape_metallic(material, iomaterial->metallic,
          app->texture_map.at(iomaterial->metallic_tex));
      set_shape_transmission(material, iomaterial->transmission,
          iomaterial->thin, iomaterial->trdepth,
          app->texture_map.at(iomaterial->transmission_tex));
      set_shape_roughness(material, iomaterial->roughness,
          app->texture_map.at(iomaterial->roughness_tex));
      set_shape_opacity(material, iomaterial->opacity,
          app->texture_map.at(iomaterial->opacity_tex));
      set_shape_thin(material, iomaterial->thin);
      set_shape_normalmap(
          material, app->texture_map.at(iomaterial->normal_tex));
      set_shape_scattering(material, iomaterial->scattering,
          iomaterial->scanisotropy,
          app->texture_map.at(iomaterial->scattering_tex));
      init_lights(app->scene.get());
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->textures.empty() &&
      begin_glheader(win, "textures")) {
    draw_glcombobox(
        win, "textures##2", app->selected_texture, app->ioscene->textures);
    if (draw_glwidgets_texture(win, app, app->selected_texture)) {
      stop_display(app);
      auto texture   = app->scene->textures[app->selected_texture];
      auto iotexture = app->ioscene->textures[app->selected_texture];
      if (!iotexture->hdr.empty()) {
        set_texture(texture, iotexture->hdr);
      } else if (!iotexture->ldr.empty()) {
        set_texture(texture, iotexture->ldr);
      }
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->subdivs.empty() && begin_glheader(win, "subdivs")) {
    draw_glcombobox(
        win, "selection##2", app->selected_subdiv, app->ioscene->subdivs);
    if (draw_glwidgets_subdiv(win, app, app->selected_subdiv)) {
      stop_display(app);
      auto iosubdiv = app->ioscene->subdivs[app->selected_subdiv];
      tesselate_subdiv(app->ioscene.get(), iosubdiv);
      // TODO: this is bogus
      auto shape   = app->scene->shapes[app->selected_subdiv];
      auto ioshape = iosubdiv->shape;
      set_shape_points(shape, ioshape->points);
      set_shape_lines(shape, ioshape->lines);
      set_shape_triangles(shape, ioshape->triangles);
      set_shape_quads(shape, ioshape->quads);
      set_shape_positions(shape, ioshape->positions);
      set_shape_normals(shape, ioshape->normals);
      set_shape_texcoords(shape, ioshape->texcoords);
      set_shape_colors(shape, ioshape->colors);
      set_shape_radius(shape, ioshape->radius);
      set_shape_tangents(shape, ioshape->tangents);
      update_bvh(app->scene.get(), {}, {app->selected_subdiv}, app->params);
      // TODO: maybe we should update lights for this
      reset_display(app);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->states.empty() && apps->selected >= 0) {
    auto app                  = apps->states[apps->selected];
    app->glparams.window      = input.window_size;
    app->glparams.framebuffer = input.framebuffer_viewport;
    if (!is_initialized(app->glimage)) init_glimage(app->glimage);
    if (!app->render_counter)
      set_glimage(app->glimage, app->display, false, false);
    update_imview(app->glparams.center, app->glparams.scale,
        app->display.size(), app->glparams.window, app->glparams.fit);
    draw_glimage(app->glimage, app->glparams);
    app->render_counter++;
    if (app->render_counter > 10) app->render_counter = 0;
  }
}

void update(const opengl_window& win, shared_ptr<app_states> apps) {
  auto is_ready = [](const future<shared_ptr<app_state>>& result) -> bool {
    return result.valid() &&
           result.wait_for(chrono::microseconds(0)) == future_status::ready;
  };

  while (!apps->loaders.empty() && is_ready(apps->loaders.front())) {
    try {
      auto app = apps->loaders.front().get();
      apps->loaders.pop_front();
      apps->states.push_back(app);
      reset_display(app);
      if (apps->selected < 0) apps->selected = (int)apps->states.size() - 1;
    } catch (std::exception& e) {
      apps->loaders.pop_front();
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
  }
}

void run_app(int argc, const char* argv[]) {
  // application
  auto apps      = make_shared<app_states>();
  auto filenames = vector<string>{};

  // parse command line
  auto cli = make_cli("yscnitrace", "progressive path tracing");
  add_cli_option(cli, "--camera", apps->params.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", apps->params.resolution, "Image resolution.");
  add_cli_option(
      cli, "--samples,-s", apps->params.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)apps->params.sampler, "Tracer type.",
      trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)apps->params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", apps->params.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", apps->params.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", apps->params.tentfilter, "Filter image.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", apps->params.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(
      cli, "--bvh", (int&)apps->params.bvh, "Bvh type", trace_bvh_names);
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  // window
  auto win = opengl_window{};
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", true);

  // callbacks
  set_draw_glcallback(
      win, [apps](const opengl_window& win, const opengl_input& input) {
        draw(win, apps, input);
      });
  set_widgets_glcallback(
      win, [apps](const opengl_window& win, const opengl_input& input) {
        draw_glwidgets(win, apps, input);
      });
  set_drop_glcallback(
      win, [apps](const opengl_window& win, const vector<string>& paths,
               const opengl_input& input) {
        for (auto& path : paths) load_scene_async(apps, path);
      });
  set_update_glcallback(
      win, [apps](const opengl_window& win, const opengl_input& input) {
        update(win, apps);
      });
  set_uiupdate_glcallback(win, [apps](const opengl_window& win,
                                   const opengl_input&     input) {
    auto scene_ok = !apps->states.empty() && apps->selected >= 0;
    if (!scene_ok) return;

    // handle mouse and keyboard for navigation
    if (scene_ok && (input.mouse_left || input.mouse_right) &&
        !input.modifier_alt && !input.widgets_active) {
      auto app      = apps->states[apps->selected];
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
      set_camera_frame(
          app->scene->cameras[app->params.camera], iocamera->frame);
      set_camera_lens(app->scene->cameras[app->params.camera], iocamera->lens,
          iocamera->aspect, iocamera->film);
      set_camera_focus(app->scene->cameras[app->params.camera],
          iocamera->aperture, iocamera->focus);
      reset_display(app);
    }

    // selection
    if (scene_ok && (input.mouse_left || input.mouse_right) &&
        input.modifier_alt && !input.widgets_active) {
      auto app = apps->states[apps->selected];
      auto ij  = get_image_coords(input.mouse_pos, app->glparams.center,
          app->glparams.scale, app->render.size());
      if (ij.x >= 0 && ij.x < app->render.size().x && ij.y >= 0 &&
          ij.y < app->render.size().y) {
        auto camera = app->scene->cameras.at(app->params.camera);
        auto ray    = camera_ray(camera->frame, camera->lens, camera->film,
            vec2f{ij.x + 0.5f, ij.y + 0.5f} / vec2f{(float)app->render.size().x,
                                                  (float)app->render.size().y});
        if (auto isec = intersect_scene_bvh(app->scene.get(), ray); isec.hit) {
          app->selected_object = isec.object;
        }
      }
    }
  });

  // run ui
  run_ui(win);

  // clear
  clear_glwindow(win);
}

int main(int argc, const char* argv[]) {
  try {
    run_app(argc, argv);
    return 0;
  } catch (std::exception& e) {
    print_fatal(e.what());
    return 1;
  }
}
