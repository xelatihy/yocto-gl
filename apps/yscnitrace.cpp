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

#include <future>
#include <list>

namespace yocto {
void print_obj_camera(const sceneio_camera& camera);
};  // namespace yocto

// Application scene
struct app_state {
  // loading options
  string filename  = "app.yaml";
  string imagename = "out.png";
  string outname   = "out.yaml";
  string name      = "";

  // options
  trace_params params = {};
  int          pratio = 8;

  // scene
  sceneio_model ioscene    = {};
  trace_scene   scene      = {};
  bool          add_skyenv = false;

  // rendering state
  trace_state  state    = {};
  image<vec4f> render   = {};
  image<vec4f> display  = {};
  float        exposure = 0;

  // view scene
  opengl_image        glimage  = {};
  draw_glimage_params glparams = {};

  // editing
  pair<string, int> selection = {"camera", 0};

  // computation
  int               render_sample  = 0;
  std::atomic<bool> render_stop    = {};
  std::future<void> render_future  = {};
  int               render_counter = 0;

  ~app_state() {
    render_stop = true;
    if (render_future.valid()) render_future.get();
  }
};

// Application state
struct app_states {
  // data
  std::list<app_state>                   states;
  int                                    selected = -1;
  std::list<app_state>                   loading;
  std::list<std::future<sceneio_status>> loaders;

  // get image
  app_state& get_selected() {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }
  const app_state& get_selected() const {
    auto it = states.begin();
    std::advance(it, selected);
    return *it;
  }

  // default options
  trace_params params     = {};
  bool         add_skyenv = false;
};

// convert scene objects
void update_trace_camera(trace_camera& camera, const sceneio_camera& iocamera) {
  camera.frame = iocamera.frame;
  camera.film  = iocamera.aspect >= 1
                    ? vec2f{iocamera.film, iocamera.film / iocamera.aspect}
                    : vec2f{iocamera.film * iocamera.aspect, iocamera.film};
  camera.lens     = iocamera.lens;
  camera.focus    = iocamera.focus;
  camera.aperture = iocamera.aperture;
}
void update_trace_texture(
    trace_texture& texture, const sceneio_texture& iotexture) {
  texture.hdr = iotexture.hdr;
  texture.ldr = iotexture.ldr;
}
void update_trace_material(
    trace_material& material, const sceneio_material& iomaterial) {
  material.emission         = iomaterial.emission;
  material.diffuse          = iomaterial.diffuse;
  material.specular         = iomaterial.specular;
  material.transmission     = iomaterial.transmission;
  material.roughness        = iomaterial.roughness;
  material.opacity          = iomaterial.opacity;
  material.refract          = iomaterial.refract;
  material.volemission      = iomaterial.volemission;
  material.voltransmission  = iomaterial.voltransmission;
  material.volmeanfreepath  = iomaterial.volmeanfreepath;
  material.volscatter       = iomaterial.volscatter;
  material.volscale         = iomaterial.volscale;
  material.volanisotropy    = iomaterial.volanisotropy;
  material.emission_tex     = iomaterial.emission_tex;
  material.diffuse_tex      = iomaterial.diffuse_tex;
  material.specular_tex     = iomaterial.specular_tex;
  material.transmission_tex = iomaterial.transmission_tex;
  material.roughness_tex    = iomaterial.roughness_tex;
  material.opacity_tex      = iomaterial.opacity_tex;
  material.subsurface_tex   = iomaterial.subsurface_tex;
  material.normal_tex       = iomaterial.normal_tex;
}
void update_trace_shape(trace_shape& shape, const sceneio_shape& ioshape,
    const sceneio_model& ioscene) {
  if (needs_tesselation(ioscene, ioshape)) {
    return update_trace_shape(
        shape, tesselate_shape(ioscene, ioshape), ioscene);
  }
  shape.points        = ioshape.points;
  shape.lines         = ioshape.lines;
  shape.triangles     = ioshape.triangles;
  shape.quads         = ioshape.quads;
  shape.quadspos      = ioshape.quadspos;
  shape.quadsnorm     = ioshape.quadsnorm;
  shape.quadstexcoord = ioshape.quadstexcoord;
  shape.positions     = ioshape.positions;
  shape.normals       = ioshape.normals;
  shape.texcoords     = ioshape.texcoords;
  shape.colors        = ioshape.colors;
  shape.radius        = ioshape.radius;
  shape.tangents      = ioshape.tangents;
}
void update_trace_instance(
    trace_instance& instance, const sceneio_instance& ioinstance) {
  instance.frame    = ioinstance.frame;
  instance.shape    = ioinstance.shape;
  instance.material = ioinstance.material;
}
void update_trace_environment(
    trace_environment& environment, const sceneio_environment& ioenvironment) {
  environment.frame        = ioenvironment.frame;
  environment.emission     = ioenvironment.emission;
  environment.emission_tex = ioenvironment.emission_tex;
}

// Construct a scene from io
trace_scene make_scene(const sceneio_model& ioscene) {
  auto scene = trace_scene{};

  for (auto& iocamera : ioscene.cameras) {
    update_trace_camera(scene.cameras.emplace_back(), iocamera);
  }
  for (auto& iotexture : ioscene.textures) {
    update_trace_texture(scene.textures.emplace_back(), iotexture);
  }
  for (auto& iomaterial : ioscene.materials) {
    update_trace_material(scene.materials.emplace_back(), iomaterial);
  }
  for (auto& ioshape : ioscene.shapes) {
    update_trace_shape(scene.shapes.emplace_back(), ioshape, ioscene);
  }
  for (auto& ioinstance : ioscene.instances) {
    update_trace_instance(scene.instances.emplace_back(), ioinstance);
  }
  for (auto& ioenvironment : ioscene.environments) {
    update_trace_environment(scene.environments.emplace_back(), ioenvironment);
  }

  return scene;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto             futures  = vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, size]() {
          while (true) {
            auto j = next_idx.fetch_add(1);
            if (j >= size.y) break;
            for (auto i = 0; i < size.x; i++) func({i, j});
          }
        }));
  }
  for (auto& f : futures) f.get();
}

void stop_display(app_state& app) {
  // stop render
  app.render_stop = true;
  if (app.render_future.valid()) app.render_future.get();
}

void reset_display(app_state& app) {
  // stop render
  app.render_stop = true;
  if (app.render_future.valid()) app.render_future.get();

  // reset state
  app.state = make_state(app.scene, app.params);
  app.render.resize(app.state.size());
  app.display.resize(app.state.size());

  // render preview
  auto preview_prms = app.params;
  preview_prms.resolution /= app.pratio;
  preview_prms.samples = 1;
  auto preview         = trace_image(app.scene, preview_prms);
  preview              = tonemap_image(preview, app.exposure);
  for (auto j = 0; j < app.display.size().y; j++) {
    for (auto i = 0; i < app.display.size().x; i++) {
      auto pi             = clamp(i / app.pratio, 0, preview.size().x - 1),
           pj             = clamp(j / app.pratio, 0, preview.size().y - 1);
      app.display[{i, j}] = preview[{pi, pj}];
    }
  }

  // start renderer
  app.render_counter = 0;
  app.render_stop    = false;
  app.render_future  = std::async(std::launch::async, [&app]() {
    for (auto sample = 0; sample < app.params.samples; sample++) {
      if (app.render_stop) return;
      parallel_for(app.render.size(), [&app](const vec2i& ij) {
        if (app.render_stop) return;
        app.render[ij]  = trace_sample(app.state, app.scene, ij, app.params);
        app.display[ij] = tonemap(app.render[ij], app.exposure);
      });
    }
  });
}

void load_scene_async(app_states& apps, const string& filename) {
  auto& app      = apps.loading.emplace_back();
  app.filename   = filename;
  app.imagename  = replace_extension(filename, ".png");
  app.outname    = replace_extension(filename, ".edited.yaml");
  app.name       = get_filename(app.filename);
  app.params     = app.params;
  app.add_skyenv = app.add_skyenv;
  apps.loaders.push_back(
      std::async(std::launch::async, [&app]() -> sceneio_status {
        if (auto ret = load_scene(app.filename, app.ioscene); !ret) return ret;
        app.scene = make_scene(app.ioscene);
        init_bvh(app.scene, app.params);
        init_lights(app.scene);
        if (app.scene.lights.empty() && is_sampler_lit(app.params)) {
          app.params.sampler = trace_sampler_type::eyelight;
        }
        app.state = make_state(app.scene, app.params);
        app.render.resize(app.state.size());
        app.display.resize(app.state.size());
        app.name = get_filename(app.filename) + " [" +
                   std::to_string(app.render.size().x) + "x" +
                   std::to_string(app.render.size().y) + " @ 0]";
        return {};
      }));
}

bool draw_glwidgets_camera(const opengl_window& win, app_state& app, int id) {
  auto& camera = app.ioscene.cameras[id];
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

bool draw_glwidgets_texture(const opengl_window& win, app_state& app, int id) {
  auto& texture      = app.ioscene.textures[id];
  auto  old_filename = texture.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", texture.name);
  edited += draw_gltextinput(win, "filename", texture.filename);
  draw_gllabel(win, "hdr",
      std::to_string(texture.hdr.size().x) + " x " +
          std::to_string(texture.hdr.size().y));
  draw_gllabel(win, "ldr",
      std::to_string(texture.ldr.size().x) + " x " +
          std::to_string(texture.ldr.size().y));
  if (edited && old_filename != texture.filename) {
    if(auto ret = load_texture(app.filename, texture); !ret) {
      push_glmessage(ret.error);
      log_glinfo(win, ret.error);
    }
  }
  return edited;
}

bool draw_glwidgets_material(const opengl_window& win, app_state& app, int id) {
  auto& material = app.ioscene.materials[id];
  auto  edited   = 0;
  edited += draw_gltextinput(win, "name", material.name);
  edited += draw_glhdrcoloredit(win, "emission", material.emission);
  edited += draw_glcoloredit(win, "diffuse", material.diffuse);
  edited += draw_glcoloredit(win, "specular", material.specular);
  edited += draw_glslider(win, "metallic", material.metallic, 0, 1);
  edited += draw_glslider(win, "roughness", material.roughness, 0, 1);
  edited += draw_glcoloredit(win, "coat", material.coat);
  edited += draw_glcoloredit(win, "transmission", material.transmission);
  edited += draw_glcheckbox(win, "refract", material.refract);
  edited += draw_glcoloredit(win, "vol transmission", material.voltransmission);
  edited += draw_glcoloredit(win, "vol meanfreepath", material.volmeanfreepath);
  edited += draw_glcoloredit(win, "vol scatter", material.volscatter);
  edited += draw_glcoloredit(win, "vol emission", material.volemission);
  edited += draw_glslider(win, "vol scale", material.volscale, 0, 1);
  edited += draw_glslider(win, "vol anisotropy", material.volanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", material.opacity, 0, 1);
  edited += draw_glcombobox(
      win, "emission_tex", material.emission_tex, app.ioscene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, app.ioscene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, app.ioscene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, app.ioscene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      app.ioscene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      app.ioscene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, app.ioscene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, app.ioscene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_state& app, int id) {
  auto& shape        = app.ioscene.shapes[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  draw_gllabel(win, "points", std::to_string(shape.points.size()));
  draw_gllabel(win, "lines", std::to_string(shape.lines.size()));
  draw_gllabel(win, "triangles", std::to_string(shape.triangles.size()));
  draw_gllabel(win, "quads", std::to_string(shape.quads.size()));
  draw_gllabel(win, "quads pos", std::to_string(shape.quadspos.size()));
  draw_gllabel(win, "quads norm", std::to_string(shape.quadsnorm.size()));
  draw_gllabel(win, "quads texcoord", std::to_string(shape.quadstexcoord.size()));
  draw_gllabel(win, "pos", std::to_string(shape.positions.size()));
  draw_gllabel(win, "norm", std::to_string(shape.normals.size()));
  draw_gllabel(win, "texcoord", std::to_string(shape.texcoords.size()));
  draw_gllabel(win, "color", std::to_string(shape.colors.size()));
  draw_gllabel(win, "radius", std::to_string(shape.radius.size()));
  draw_gllabel(win, "tangsp", std::to_string(shape.tangents.size()));
  if (edited && old_filename != shape.filename) {
    if(auto ret = load_shape(app.filename, shape); !ret) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
    }
  }
  return edited;
}

bool draw_glwidgets_instance(const opengl_window& win, app_state& app, int id) {
  auto& instance     = app.ioscene.instances[id];
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame.x", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(win, "shape", instance.shape, app.ioscene.shapes, true);
  edited += draw_glcombobox(win, "material", instance.material, app.ioscene.materials, true);
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, app_state& app, int id) {
  auto& environment = app.ioscene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame.x", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame.y", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame.z", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      app.ioscene.textures, true);
  return edited;
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         apps     = *(app_states*)get_gluser_pointer(win);
  auto          scene_ok = !apps.states.empty() && apps.selected >= 0;
  if (!begin_glwidgets_window(win, "yscnitrace")) return;
  draw_glmessages(win);
  if (draw_glfiledialog_button(win, "load", true, "load", load_path, false,
          "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save", scene_ok, "save", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.yaml;*.obj;*.pbrt")) {
    auto& app   = apps.get_selected();
    app.outname = save_path;
    if (auto ret = save_scene(app.outname, app.ioscene); !ret) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
      log_glinfo(win, ret.error);
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save image", scene_ok, "save image",
          save_path, true, get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto& app   = apps.get_selected();
    app.outname = save_path;
    if (auto ret = save_image(app.imagename, app.display); !ret) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", scene_ok)) {
    auto it = apps.states.begin();
    std::advance(it, apps.selected);
    apps.states.erase(it);
    apps.selected = apps.states.empty() ? -1 : 0;
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  draw_glcombobox(win, "scene", apps.selected, (int)apps.states.size(),
      [&apps](int idx) {
        auto it = apps.states.begin();
        std::advance(it, idx);
        return it->name.c_str();
      },
      false);
  if (scene_ok && begin_glheader(win, "trace")) {
    auto  edited  = 0;
    auto& app     = apps.get_selected();
    auto& tparams = app.params;
    edited += draw_glcombobox(
        win, "camera", tparams.camera, app.ioscene.cameras);
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
    edited += draw_glslider(win, "pratio", app.pratio, 1, 64);
    edited += draw_glslider(win, "exposure", app.exposure, -5, 5);
    if (edited) reset_display(app);
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "inspect")) {
    auto& app = apps.get_selected();
    draw_gllabel(win, "scene", get_filename(app.filename));
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    draw_gllabel(win, "imagename", app.imagename);
    draw_gllabel(win, "image",
        std::to_string(app.render.size().x) + " x " +
            std::to_string(app.render.size().y) + " @ " +
            std::to_string(app.render_sample));
    draw_glslider(win, "zoom", app.glparams.scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", app.glparams.fit);
    continue_glline(win);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : app.ioscene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : scene_stats(app.ioscene)) print_info(stat);
    }
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, app.glparams.center, app.glparams.scale, app.render.size());
    draw_gldragger(win, "mouse", ij);
    if (ij.x >= 0 && ij.x < app.render.size().x && ij.y >= 0 &&
        ij.y < app.render.size().y) {
      draw_glcoloredit(win, "pixel", app.render[{ij.x, ij.y}]);
    } else {
      auto zero4f_ = zero4f;
      draw_glcoloredit(win, "pixel", zero4f_);
    }
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "edit")) {
    static auto labels = vector<string>{
        "camera", "shape", "environment", "instance", "material", "texture"};
    auto& app = apps.get_selected();
    if (draw_glcombobox(win, "selection##1", app.selection.first, labels))
      app.selection.second = 0;
    if (app.selection.first == "camera") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.cameras);
      if(draw_glwidgets_camera(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_camera(app.scene.cameras[app.selection.second],
            app.ioscene.cameras[app.selection.second]);
        reset_display(app);
      }
    } else if (app.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.textures);
      if(draw_glwidgets_texture(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_texture(app.scene.textures[app.selection.second],
            app.ioscene.textures[app.selection.second]);
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app.selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.materials);
      if(draw_glwidgets_material(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_material(app.scene.materials[app.selection.second],
            app.ioscene.materials[app.selection.second]);
        init_lights(app.scene);
        reset_display(app);
      }
    } else if (app.selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.shapes);
      if(draw_glwidgets_shape(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_shape(app.scene.shapes[app.selection.second],
            app.ioscene.shapes[app.selection.second], app.ioscene);
        update_bvh(app.scene, {}, {app.selection.second}, app.params);
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app.selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.instances);
      if(draw_glwidgets_instance(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_instance(app.scene.instances[app.selection.second],
            app.ioscene.instances[app.selection.second]);
        update_bvh(app.scene, {app.selection.second}, {}, app.params);
        // TODO: maybe we should update lights for this
        reset_display(app);
      }
    } else if (app.selection.first == "environment") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.ioscene.environments);
      if(draw_glwidgets_environment(win, app, app.selection.second)) {
        stop_display(app);
        update_trace_environment(app.scene.environments[app.selection.second],
            app.ioscene.environments[app.selection.second]);
        init_lights(app.scene);
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

void draw(const opengl_window& win) {
  auto& apps = *(app_states*)get_gluser_pointer(win);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!apps.states.empty() && apps.selected >= 0) {
    auto& app                = apps.get_selected();
    app.glparams.window      = get_glwindow_size(win);
    app.glparams.framebuffer = get_glframebuffer_viewport(win);
    if (!app.glimage || app.glimage.size() != app.display.size() ||
        !app.render_counter) {
      update_glimage(app.glimage, app.display, false, false);
    }
    update_imview(app.glparams.center, app.glparams.scale, app.display.size(),
        app.glparams.window, app.glparams.fit);
    draw_glimage(app.glimage, app.glparams);
    app.render_counter++;
    if (app.render_counter > 10) app.render_counter = 0;
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_states& app) {
  auto is_ready = [](const std::future<sceneio_status>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!app.loaders.empty() && is_ready(app.loaders.front())) {
    if (auto ret = app.loaders.front().get(); !ret) {
      push_glmessage(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, ret.error);
      break;
    }
    app.states.splice(app.states.end(), app.loading, app.loading.begin());
    app.loaders.pop_front();
    reset_display(app.states.back());
    if (app.selected < 0) app.selected = (int)app.states.size() - 1;
  }
}

// run ui loop
void run_ui(app_states& apps) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", &apps, draw);
  set_drop_glcallback(
      win, [](const opengl_window& win, const vector<string>& paths) {
        auto& app = *(app_states*)get_gluser_pointer(win);
        for (auto& path : paths) load_scene_async(app, path);
      });

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
    auto scene_ok       = !apps.states.empty() && apps.selected >= 0;

    // handle mouse and keyboard for navigation
    if (scene_ok && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& app    = apps.get_selected();
      auto& camera = app.ioscene.cameras.at(app.params.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      stop_display(app);
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      update_trace_camera(app.scene.cameras.at(app.params.camera), camera);
      reset_display(app);
    }

    // selection
    if (scene_ok && (mouse_left || mouse_right) && alt_down &&
        !widgets_active) {
      auto& app = apps.get_selected();
      auto  ij  = get_image_coords(mouse_pos, app.glparams.center,
          app.glparams.scale, app.render.size());
      if (ij.x >= 0 && ij.x < app.render.size().x && ij.y >= 0 &&
          ij.y < app.render.size().y) {
        auto& camera = app.scene.cameras.at(app.params.camera);
        auto  ray    = camera_ray(camera.frame, camera.lens, camera.film,
            vec2f{ij.x + 0.5f, ij.y + 0.5f} /
                vec2f{(float)app.render.size().x, (float)app.render.size().y});
        if (auto isec = intersect_scene_bvh(app.scene, ray); isec.hit) {
          app.selection = {"instance", isec.instance};
        }
      }
    }

    // update
    update(win, apps);

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
  app_states app{};
  auto       filenames = vector<string>{};

  // parse command line
  auto cli = make_cli("yscnitrace", "progressive path tracing");
  add_cli_option(cli, "--camera", app.params.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app.params.resolution, "Image resolution.");
  add_cli_option(cli, "--samples,-s", app.params.samples, "Number of samples.");
  add_cli_option(cli, "--tracer,-t", (int&)app.params.sampler, "Tracer type.",
      trace_sampler_names);
  add_cli_option(cli, "--falsecolor,-F", (int&)app.params.falsecolor,
      "Tracer false color type.", trace_falsecolor_names);
  add_cli_option(
      cli, "--bounces", app.params.bounces, "Maximum number of bounces.");
  add_cli_option(cli, "--clamp", app.params.clamp, "Final pixel clamping.");
  add_cli_option(cli, "--filter", app.params.tentfilter, "Filter image.");
  add_cli_option(cli, "--env-hidden/--no-env-hidden", app.params.envhidden,
      "Environments are hidden in renderer");
  add_cli_option(
      cli, "--bvh", (int&)app.params.bvh, "Bvh type", trace_bvh_names);
  add_cli_option(cli, "--add-skyenv", app.add_skyenv, "Add sky envmap");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) load_scene_async(app, filename);

  // run interactive
  run_ui(app);

  // done
  return 0;
}
