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
void print_obj_camera(const sceneio_camera& camera);
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
  sceneio_model           ioscene    = {};
  unique_ptr<trace_scene> scene      = {};
  bool                    add_skyenv = false;

  // rendering state
  unique_ptr<trace_state> state    = {};
  image<vec4f>            render   = {};
  image<vec4f>            display  = {};
  float                   exposure = 0;

  // view scene
  unique_ptr<opengl_image> glimage  = {};
  draw_glimage_params      glparams = {};

  // editing
  pair<string, int> selection = {"camera", 0};

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

// Load state
struct load_state {
  string                filename = "";
  shared_ptr<app_state> app      = nullptr;
  sceneio_status        status   = {};
};

// Application state
struct app_states {
  // data
  vector<shared_ptr<app_state>> states   = {};
  int                           selected = -1;

  // loading
  deque<future<load_state>> loaders = {};

  // default options
  trace_params params     = {};
  bool         add_skyenv = false;
};

// set material values
void set_material(
    trace_scene* scene, int idx, const sceneio_material& iomaterial) {
  set_material_emission(
      scene, idx, iomaterial.emission, iomaterial.emission_tex);
  set_material_opacity(scene, idx, iomaterial.opacity, iomaterial.opacity_tex);
  set_material_normalmap(scene, idx, iomaterial.normal_tex);
  switch (iomaterial.type) {
    case sceneio_material_type::standard:
    case sceneio_material_type::substrate: {
      set_material_diffuse(
          scene, idx, iomaterial.base, iomaterial.base_tex);
      set_material_specular(
          scene, idx, iomaterial.specular, iomaterial.specular_tex);
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
    }; break;
    case sceneio_material_type::matte: {
      set_material_diffuse(
          scene, idx, iomaterial.base, iomaterial.base_tex);
      set_material_specular(scene, idx, zero3f);
      set_material_roughness(scene, idx, 1);
    } break;
    case sceneio_material_type::reflective: {
      set_material_diffuse(scene, idx, zero3f);
      set_material_specular(
          scene, idx, iomaterial.base, iomaterial.base_tex);
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
    } break;
    case sceneio_material_type::metallic: {
      set_material_diffuse(
          scene, idx, iomaterial.base, iomaterial.base_tex);
      set_material_metallic(
          scene, idx, mean(iomaterial.specular), iomaterial.specular_tex);
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
    }; break;
    case sceneio_material_type::transparent: {
      set_material_specular(
          scene, idx, iomaterial.specular, iomaterial.specular_tex);
      set_material_transmission(
          scene, idx, iomaterial.base, iomaterial.base_tex);
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
    }; break;
    case sceneio_material_type::refractive: {
      set_material_specular(
          scene, idx, iomaterial.specular, iomaterial.specular_tex);
      set_material_transmission(scene, idx, {1, 1, 1});
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
      set_material_volume(scene, idx, iomaterial.volemission,
          iomaterial.base, zero3f,
          iomaterial.volume, iomaterial.volscale, iomaterial.volanisotropy,
          iomaterial.volume_tex);
      set_material_refract(scene, idx, true);
    }; break;
    case sceneio_material_type::subsurface: {
      set_material_specular(
          scene, idx, iomaterial.specular, iomaterial.specular_tex);
      set_material_transmission(scene, idx, {1, 1, 1});
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
      set_material_volume(scene, idx, iomaterial.volemission,
          iomaterial.base, zero3f,
          iomaterial.volume, iomaterial.volscale, iomaterial.volanisotropy,
          iomaterial.volume_tex);
      set_material_refract(scene, idx, true);
    }; break;
    case sceneio_material_type::volume: {
      set_material_specular(
          scene, idx, iomaterial.specular, iomaterial.specular_tex);
      set_material_transmission(scene, idx, {1, 1, 1});
      set_material_roughness(
          scene, idx, iomaterial.roughness, iomaterial.roughness_tex);
      set_material_volume(scene, idx, iomaterial.volemission,
          iomaterial.base, zero3f,
          iomaterial.volume, iomaterial.volscale, iomaterial.volanisotropy,
          iomaterial.volume_tex);
      set_material_refract(scene, idx, false);
    }; break;
  }
}

// Construct a scene from io
trace_scene* make_scene(const sceneio_model& ioscene) {
  auto scene = make_unique<trace_scene>();

  for (auto& iocamera : ioscene.cameras) {
    add_camera(scene.get(), iocamera.frame, iocamera.lens, iocamera.aspect,
        iocamera.film, iocamera.aperture, iocamera.focus);
  }
  for (auto& iotexture : ioscene.textures) {
    if (!iotexture.hdr.empty()) {
      add_texture(scene.get(), std::move(iotexture.hdr));
    } else if (!iotexture.ldr.empty()) {
      add_texture(scene.get(), std::move(iotexture.ldr));
    }
  }
  for (auto& iomaterial : ioscene.materials) {
    auto id = add_material(scene.get());
    set_material(scene.get(), id, iomaterial);
  }
  for (auto& ioshape_ : ioscene.shapes) {
    auto tshape = (needs_tesselation(ioscene, ioshape_))
                      ? tesselate_shape(ioscene, ioshape_)
                      : sceneio_shape{};
    auto& ioshape = (needs_tesselation(ioscene, ioshape_)) ? tshape : ioshape_;
    if (!ioshape.points.empty()) {
      add_shape(scene.get(), ioshape.points, ioshape.positions, ioshape.normals,
          ioshape.texcoords, ioshape.colors, ioshape.radius);
    } else if (!ioshape.lines.empty()) {
      add_shape(scene.get(), ioshape.lines, ioshape.positions, ioshape.normals,
          ioshape.texcoords, ioshape.colors, ioshape.radius);
    } else if (!ioshape.triangles.empty()) {
      add_shape(scene.get(), ioshape.triangles, ioshape.positions,
          ioshape.normals, ioshape.texcoords, ioshape.colors, ioshape.tangents);
    } else if (!ioshape.quads.empty()) {
      add_shape(scene.get(), ioshape.quads, ioshape.positions, ioshape.normals,
          ioshape.texcoords, ioshape.colors, ioshape.tangents);
    } else if (!ioshape.quadspos.empty()) {
      add_shape(scene.get(), ioshape.quadspos, ioshape.quadsnorm,
          ioshape.quadstexcoord, ioshape.positions, ioshape.normals,
          ioshape.texcoords);
    }
    tshape = {};
  }
  for (auto& ioinstance : ioscene.instances) {
    add_instance(
        scene.get(), ioinstance.frame, ioinstance.shape, ioinstance.material);
  }
  for (auto& ioenvironment : ioscene.environments) {
    add_environment(scene.get(), ioenvironment.frame, ioenvironment.emission,
        ioenvironment.emission_tex);
  }

  return scene.release();
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
  app->state = unique_ptr<trace_state>{
      make_state(app->scene.get(), app->params)};
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
      async(launch::async, [apps, filename]() -> load_state {
        auto app        = make_shared<app_state>();
        app->filename   = filename;
        app->imagename  = replace_extension(filename, ".png");
        app->outname    = replace_extension(filename, ".edited.yaml");
        app->name       = get_filename(app->filename);
        app->params     = app->params;
        app->add_skyenv = app->add_skyenv;
        if (auto ret = load_scene(app->filename, app->ioscene); !ret)
          return {filename, nullptr, ret};
        app->scene = unique_ptr<trace_scene>{make_scene(app->ioscene)};
        init_bvh(app->scene.get(), app->params);
        init_lights(app->scene.get());
        if (app->scene->lights.empty() && is_sampler_lit(app->params)) {
          app->params.sampler = trace_sampler_type::eyelight;
        }
        app->state = unique_ptr<trace_state>{
            make_state(app->scene.get(), app->params)};
        app->render.resize(app->state->size());
        app->display.resize(app->state->size());
        app->name = get_filename(app->filename) + " [" +
                    to_string(app->render.size().x) + "x" +
                    to_string(app->render.size().y) + " @ 0]";
        return {filename, app, {}};
      }));
}

bool draw_glwidgets_camera(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& camera = app->ioscene.cameras[id];
  auto  edited = 0;
  edited += (int)draw_gltextinput(win, "name", camera.name);
  edited += (int)draw_glslider(win, "frame.x", camera.frame.x, -1, 1);
  edited += (int)draw_glslider(win, "frame.y", camera.frame.y, -1, 1);
  edited += (int)draw_glslider(win, "frame.z", camera.frame.z, -1, 1);
  edited += (int)draw_glslider(win, "frame.o", camera.frame.o, -10, 10);
  edited += (int)draw_glcheckbox(win, "ortho", camera.orthographic);
  edited += (int)draw_glslider(win, "lens", camera.lens, 0.01f, 1);
  edited += (int)draw_glslider(win, "film", camera.film, 0.01f, 0.1f);
  edited += (int)draw_glslider(win, "focus", camera.focus, 0.01f, 1000);
  edited += (int)draw_glslider(win, "aperture", camera.aperture, 0, 5);
  auto from         = camera.frame.o,
       to           = camera.frame.o - camera.focus * camera.frame.z;
  auto from_changed = draw_glslider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_glslider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    camera.frame = lookat_frame(from, to, {0, 1, 0});
    camera.focus = length(from - to);
    edited += 1;
  }
  return edited;
}

bool draw_glwidgets_texture(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& texture      = app->ioscene.textures[id];
  auto  old_filename = texture.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", texture.name);
  edited += draw_gltextinput(win, "filename", texture.filename);
  draw_gllabel(win, "hdr",
      to_string(texture.hdr.size().x) + " x " +
          to_string(texture.hdr.size().y));
  draw_gllabel(win, "ldr",
      to_string(texture.ldr.size().x) + " x " +
          to_string(texture.ldr.size().y));
  if (edited && old_filename != texture.filename) {
    if (auto ret = load_texture(app->filename, texture); !ret) {
      push_glmessage(win, ret.error);
      log_glinfo(win, ret.error);
    }
  }
  return edited;
}

bool draw_glwidgets_material(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& material = app->ioscene.materials[id];
  auto  edited   = 0;
  edited += draw_gltextinput(win, "name", material.name);
  edited += draw_glhdrcoloredit(win, "emission", material.emission);
  edited += draw_glcoloredit(win, "diffuse", material.diffuse);
  edited += draw_glcoloredit(win, "specular", material.specular);
  edited += draw_glslider(win, "metallic", material.metallic, 0, 1);
  edited += draw_glslider(win, "roughness", material.roughness, 0, 1);
  edited += draw_glcoloredit(win, "transmission", material.transmission);
  edited += draw_glcoloredit(win, "vol scatter", material.volscatter);
  edited += draw_glcoloredit(win, "vol emission", material.volemission);
  edited += draw_glslider(win, "vol scale", material.volscale, 0, 1);
  edited += draw_glslider(win, "vol anisotropy", material.volanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", material.opacity, 0, 1);
  edited += draw_glcombobox(
      win, "emission_tex", material.emission_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      app->ioscene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      app->ioscene.textures, true);
  edited += draw_glcombobox(win, "roughness_tex", material.roughness_tex,
      app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, app->ioscene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  return edited;
}

bool draw_glwidgets_shape(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& shape        = app->ioscene.shapes[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  draw_gllabel(win, "points", to_string(shape.points.size()));
  draw_gllabel(win, "lines", to_string(shape.lines.size()));
  draw_gllabel(win, "triangles", to_string(shape.triangles.size()));
  draw_gllabel(win, "quads", to_string(shape.quads.size()));
  draw_gllabel(win, "quads pos", to_string(shape.quadspos.size()));
  draw_gllabel(win, "quads norm", to_string(shape.quadsnorm.size()));
  draw_gllabel(win, "quads texcoord", to_string(shape.quadstexcoord.size()));
  draw_gllabel(win, "pos", to_string(shape.positions.size()));
  draw_gllabel(win, "norm", to_string(shape.normals.size()));
  draw_gllabel(win, "texcoord", to_string(shape.texcoords.size()));
  draw_gllabel(win, "color", to_string(shape.colors.size()));
  draw_gllabel(win, "radius", to_string(shape.radius.size()));
  draw_gllabel(win, "tangsp", to_string(shape.tangents.size()));
  if (edited && old_filename != shape.filename) {
    if (auto ret = load_shape(app->filename, shape); !ret) {
      push_glmessage(win, "cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
    }
  }
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& instance = app->ioscene.instances[id];
  auto  edited   = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame.x", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", instance.shape, app->ioscene.shapes, true);
  edited += draw_glcombobox(
      win, "material", instance.material, app->ioscene.materials, true);
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& environment = app->ioscene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame.x", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      app->ioscene.textures, true);
  return edited;
}

void draw_glwidgets(const opengl_window* win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  static string load_path = "", save_path = "", error_message = "";
  auto          scene_ok = !apps->states.empty() && apps->selected >= 0;
  if (draw_glfiledialog_button(win, "load", true, "load", load_path, false,
          "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save", scene_ok, "save", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.yaml;*.obj;*.pbrt")) {
    auto app     = apps->states[apps->selected];
    app->outname = save_path;
    if (auto ret = save_scene(app->outname, app->ioscene); !ret) {
      push_glmessage(win, "cannot save " + app->outname);
      log_glinfo(win, "cannot save " + app->outname);
      log_glinfo(win, ret.error);
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save image", scene_ok, "save image",
          save_path, true, get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto app     = apps->states[apps->selected];
    app->outname = save_path;
    if (auto ret = save_image(app->imagename, app->display); !ret) {
      push_glmessage(win, "cannot save " + app->outname);
      log_glinfo(win, "cannot save " + app->outname);
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", scene_ok)) {
    auto it = apps->states.begin();
    advance(it, apps->selected);
    apps->states.erase(it);
    apps->selected = apps->states.empty() ? -1 : 0;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_close(win, true);
  }
  draw_glcombobox(win, "scene", apps->selected, (int)apps->states.size(),
      [apps](int idx) { return apps->states[apps->selected]->name.c_str(); },
      false);
  if (scene_ok && begin_glheader(win, "trace")) {
    auto  edited  = 0;
    auto  app     = apps->states[apps->selected];
    auto& tparams = app->params;
    edited += draw_glcombobox(
        win, "camera", tparams.camera, app->ioscene.cameras);
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
  if (scene_ok && begin_glheader(win, "inspect")) {
    auto app = apps->states[apps->selected];
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
      for (auto& camera : app->ioscene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : scene_stats(app->ioscene)) print_info(stat);
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
  if (scene_ok && begin_glheader(win, "edit")) {
    static auto labels = vector<string>{
        "camera", "shape", "environment", "instance", "material", "texture"};
    auto app = apps->states[apps->selected];
    if (draw_glcombobox(win, "selection##1", app->selection.first, labels))
      app->selection.second = 0;
    if (app->selection.first == "camera") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->ioscene.cameras);
      if (draw_glwidgets_camera(win, app, app->selection.second)) {
        stop_display(app);
        auto& iocamera = app->ioscene.cameras[app->selection.second];
        set_camera(app->scene.get(), app->selection.second, iocamera.frame,
            iocamera.lens, iocamera.aspect, iocamera.film, iocamera.aperture,
            iocamera.focus);
        reset_display(app);
      }
    } else if (app->selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->ioscene.textures);
      if (draw_glwidgets_texture(win, app, app->selection.second)) {
        stop_display(app);
        auto& iotexture = app->ioscene.textures[app->selection.second];
        if (!iotexture.hdr.empty()) {
          set_texture(app->scene.get(), app->selection.second, iotexture.hdr);
        } else if (!iotexture.ldr.empty()) {
          set_texture(app->scene.get(), app->selection.second, iotexture.ldr);
        }
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app->selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->ioscene.materials);
      if (draw_glwidgets_material(win, app, app->selection.second)) {
        stop_display(app);
        auto& iomaterial = app->ioscene.materials[app->selection.second];
        set_material(app->scene.get(), app->selection.second, iomaterial);
        init_lights(app->scene.get());
        reset_display(app);
      }
    } else if (app->selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->ioscene.shapes);
      if (draw_glwidgets_shape(win, app, app->selection.second)) {
        stop_display(app);
        auto& ioshape_ = app->ioscene.shapes[app->selection.second];
        auto  tshape   = (needs_tesselation(app->ioscene, ioshape_))
                          ? tesselate_shape(app->ioscene, ioshape_)
                          : sceneio_shape{};
        auto& ioshape = (needs_tesselation(app->ioscene, ioshape_)) ? tshape
                                                                    : ioshape_;
        if (!ioshape.points.empty()) {
          set_shape(app->scene.get(), app->selection.second, ioshape.points,
              ioshape.positions, ioshape.normals, ioshape.texcoords,
              ioshape.colors, ioshape.radius);
        } else if (!ioshape.lines.empty()) {
          set_shape(app->scene.get(), app->selection.second, ioshape.lines,
              ioshape.positions, ioshape.normals, ioshape.texcoords,
              ioshape.colors, ioshape.radius);
        } else if (!ioshape.triangles.empty()) {
          set_shape(app->scene.get(), app->selection.second, ioshape.triangles,
              ioshape.positions, ioshape.normals, ioshape.texcoords,
              ioshape.colors, ioshape.tangents);
        } else if (!ioshape.quads.empty()) {
          set_shape(app->scene.get(), app->selection.second, ioshape.quads,
              ioshape.positions, ioshape.normals, ioshape.texcoords,
              ioshape.colors, ioshape.tangents);
        } else if (!ioshape.quadspos.empty()) {
          set_shape(app->scene.get(), app->selection.second, ioshape.quadspos,
              ioshape.quadsnorm, ioshape.quadstexcoord, ioshape.positions,
              ioshape.normals, ioshape.texcoords);
        }
        update_bvh(app->scene.get(), {}, {app->selection.second}, app->params);
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app->selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->ioscene.instances);
      if (draw_glwidgets_instance(win, app, app->selection.second)) {
        stop_display(app);
        auto& ioinstance = app->ioscene.instances[app->selection.second];
        set_instance(app->scene.get(), app->selection.second, ioinstance.frame,
            ioinstance.shape, ioinstance.material);
        update_bvh(app->scene.get(), {app->selection.second}, {}, app->params);
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app->selection.first == "environment") {
      draw_glcombobox(win, "selection##2", app->selection.second,
          app->ioscene.environments);
      if (draw_glwidgets_environment(win, app, app->selection.second)) {
        stop_display(app);
        auto& ioenvironment = app->ioscene.environments[app->selection.second];
        set_environment(app->scene.get(), app->selection.second,
            ioenvironment.frame, ioenvironment.emission,
            ioenvironment.emission_tex);
        init_lights(app->scene.get());
        reset_display(app);
      }
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window* win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->states.empty() && apps->selected >= 0) {
    auto app                  = apps->states[apps->selected];
    app->glparams.window      = input.window_size;
    app->glparams.framebuffer = input.framebuffer_viewport;
    if (!app->glimage) app->glimage = unique_ptr<opengl_image>(make_glimage());
    if (!app->render_counter)
      set_glimage(app->glimage.get(), app->display, false, false);
    update_imview(app->glparams.center, app->glparams.scale,
        app->display.size(), app->glparams.window, app->glparams.fit);
    draw_glimage(app->glimage.get(), app->glparams);
    app->render_counter++;
    if (app->render_counter > 10) app->render_counter = 0;
  }
}

void update(const opengl_window* win, shared_ptr<app_states> apps) {
  auto is_ready = [](const future<load_state>& result) -> bool {
    return result.valid() &&
           result.wait_for(chrono::microseconds(0)) == future_status::ready;
  };

  while (!apps->loaders.empty() && is_ready(apps->loaders.front())) {
    auto [filename, app, status] = apps->loaders.front().get();
    apps->loaders.pop_front();
    if (!status) {
      push_glmessage(win, "cannot load scene " + filename);
      log_glinfo(win, "cannot load scene " + filename);
      log_glinfo(win, status.error);
    } else {
      apps->states.push_back(app);
      reset_display(app);
      if (apps->selected < 0) apps->selected = (int)apps->states.size() - 1;
    }
  }
}

int main(int argc, const char* argv[]) {
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
  add_cli_option(cli, "--add-skyenv", apps->add_skyenv, "Add sky envmap");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  // window
  auto win = make_glwindow({1280 + 320, 720}, "yscnitrace", true);

  // callbacks
  set_draw_glcallback(
      win, [apps](const opengl_window* win, const opengl_input& input) {
        draw(win, apps, input);
      });
  set_widgets_glcallback(
      win, [apps](const opengl_window* win, const opengl_input& input) {
        draw_glwidgets(win, apps, input);
      });
  set_drop_glcallback(
      win, [apps](const opengl_window* win, const vector<string>& paths,
               const opengl_input& input) {
        for (auto& path : paths) load_scene_async(apps, path);
      });
  set_update_glcallback(
      win, [apps](const opengl_window* win, const opengl_input& input) {
        update(win, apps);
      });
  set_uiupdate_glcallback(win, [apps](const opengl_window* win,
                                   const opengl_input&     input) {
    auto scene_ok = !apps->states.empty() && apps->selected >= 0;
    if (!scene_ok) return;

    // handle mouse and keyboard for navigation
    if (scene_ok && (input.mouse_left || input.mouse_right) &&
        !input.modifier_alt && !input.widgets_active) {
      auto  app    = apps->states[apps->selected];
      auto& camera = app->ioscene.cameras.at(app->params.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) * camera.focus / 200.0f;
      pan.x = -pan.x;
      stop_display(app);
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      set_camera(app->scene.get(), app->params.camera, camera.frame,
          camera.lens, camera.aspect, camera.film, camera.aperture,
          camera.focus);
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
        auto& camera = app->scene->cameras.at(app->params.camera);
        auto  ray    = camera_ray(camera.frame, camera.lens, camera.film,
            vec2f{ij.x + 0.5f, ij.y + 0.5f} / vec2f{(float)app->render.size().x,
                                                  (float)app->render.size().y});
        if (auto isec = intersect_scene_bvh(app->scene.get(), ray); isec.hit) {
          app->selection = {"instance", isec.instance};
        }
      }
    }
  });

  // run ui
  run_ui(win);

  // clear
  delete_glwindow(win);

  // done
  return 0;
}
