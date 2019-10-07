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

#include "../yocto/yocto_common.h"
#include "../yocto/yocto_commonio.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "../yocto/yocto_trace.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <list>
#include <map>

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};  // namespace yocto

// Application scene
struct app_state {
  // loading options
  string filename  = "app.yaml";
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

  // scene
  yocto_scene scene      = {};
  trace_bvh   bvh        = {};
  bool        add_skyenv = false;

  // rendering state
  trace_lights lights  = {};
  trace_state  state   = {};
  image<vec4f> render  = {};
  image<vec4f> display = {};

  // view scene
  vec2f        image_center   = zero2f;
  float        image_scale    = 1;
  bool         zoom_to_fit    = true;
  bool         navigation_fps = false;
  opengl_image gl_image       = {};

  // editing
  pair<string, int> selection = {"camera", 0};

  // computation
  bool                 render_preview = true;
  int                  render_sample  = 0;
  int                  render_region  = 0;
  vector<image_region> render_regions = {};
};

// Application state
struct app_states {
  // data
  std::list<app_state>                   states;
  int                                    selected = -1;
  std::list<string>                      errors;
  std::list<app_state>                   loading;
  std::list<std::future<void>>           load_workers;
  std::deque<std::unique_ptr<app_state>> minchia;

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
  load_params    load_prms    = {};
  save_params    save_prms    = {};
  bvh_params     bvh_prms     = {};
  trace_params   trace_prms   = {};
  tonemap_params tonemap_prms = {};
  bool           add_skyenv   = false;
};

void reset_display(app_state& app) {
  auto image_size = camera_resolution(
      app.scene.cameras[app.trace_prms.camera], app.trace_prms.resolution);
  app.render.resize(image_size);
  app.display.resize(image_size);
  app.render_preview = true;
  app.render_sample  = 0;
  app.render_region  = 0;
  app.state          = make_trace_state(app.render.size(), app.trace_prms.seed);
  app.render_regions = make_image_regions(
      app.render.size(), app.trace_prms.region, true);
}

void load_scene_async(app_states& apps, const string& filename) {
  auto& app        = apps.loading.emplace_back();
  app.filename     = filename;
  app.imagename    = replace_extension(filename, ".png");
  app.outname      = replace_extension(filename, ".edited.yaml");
  app.name         = get_filename(app.filename);
  app.load_prms    = app.load_prms;
  app.save_prms    = app.save_prms;
  app.trace_prms   = app.trace_prms;
  app.bvh_prms     = app.bvh_prms;
  app.tonemap_prms = app.tonemap_prms;
  app.add_skyenv   = app.add_skyenv;
  apps.load_workers.push_back(run_async([&app]() {
    load_scene(app.filename, app.scene);
    make_bvh(app.bvh, app.scene, app.bvh_prms);
    make_trace_lights(app.lights, app.scene);
    if (app.lights.instances.empty() && app.lights.environments.empty() &&
        is_sampler_lit(app.trace_prms)) {
      app.trace_prms.sampler = trace_params::sampler_type::eyelight;
    }
    auto image_size = camera_resolution(
        app.scene.cameras[app.trace_prms.camera], app.trace_prms.resolution);
    app.render.resize(image_size);
    app.display.resize(image_size);
    app.name = get_filename(app.filename) + " [" +
               std::to_string(app.render.size().x) + "x" +
               std::to_string(app.render.size().y) + " @ 0]";
  }));
}

bool draw_glwidgets_camera(const opengl_window& win, app_state& app, int id) {
  auto& camera = app.scene.cameras[id];
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

/// Visit struct elements.
bool draw_glwidgets_texture(const opengl_window& win, app_state& app, int id) {
  auto& texture      = app.scene.textures[id];
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
    try {
      if (is_hdr_filename(texture.filename)) {
        load_image(texture.filename, texture.hdr);
      } else {
        load_imageb(texture.filename, texture.ldr);
      }
    } catch (std::exception& e) {
      push_glmessage("cannot load " + texture.filename);
      log_glinfo(win, "cannot load " + texture.filename);
      log_glinfo(win, e.what());
    }
    // TODO: update lights
  }
  return edited;
}

bool draw_glwidgets_material(const opengl_window& win, app_state& app, int id) {
  auto& material = app.scene.materials[id];
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
      win, "emission_tex", material.emission_tex, app.scene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, app.scene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, app.scene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, app.scene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      app.scene.textures, true);
  edited += draw_glcombobox(
      win, "subsurface_tex", material.subsurface_tex, app.scene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, app.scene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, app.scene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_state& app, int id) {
  auto& shape        = app.scene.shapes[id];
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
  draw_gllabel(
      win, "quads texcoord", std::to_string(shape.quadstexcoord.size()));
  draw_gllabel(win, "pos", std::to_string(shape.positions.size()));
  draw_gllabel(win, "norm", std::to_string(shape.normals.size()));
  draw_gllabel(win, "texcoord", std::to_string(shape.texcoords.size()));
  draw_gllabel(win, "color", std::to_string(shape.colors.size()));
  draw_gllabel(win, "radius", std::to_string(shape.radius.size()));
  draw_gllabel(win, "tangsp", std::to_string(shape.tangents.size()));
  if (edited && old_filename != shape.filename) {
    try {
      load_shape(shape.filename, shape.points, shape.lines, shape.triangles,
          shape.quads, shape.positions, shape.normals, shape.texcoords,
          shape.colors, shape.radius);
    } catch (std::exception& e) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
      log_glinfo(win, e.what());
    }
    update_bvh(app.bvh, app.scene, {}, {id}, app.bvh_prms);
    // TODO: update lights
  }
  return edited;
}

bool draw_glwidgets_instance(const opengl_window& win, app_state& app, int id) {
  auto& instance     = app.scene.instances[id];
  auto  old_instance = instance;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame[0]", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", instance.shape, app.scene.shapes, true);
  edited += draw_glcombobox(
      win, "material", instance.material, app.scene.materials, true);
  if (edited && instance.shape != old_instance.shape)
    update_bvh(app.bvh, app.scene, {}, {id}, app.bvh_prms);
  if (edited && instance.frame != old_instance.frame)
    update_bvh(app.bvh, app.scene, {}, {id}, app.bvh_prms);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, app_state& app, int id) {
  auto& environment = app.scene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame[0]", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      app.scene.textures, true);
  if (edited) {
    // TODO: update lights
  }
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
    try {
      save_scene(app.outname, app.scene);
    } catch (std::exception& e) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
      log_glinfo(win, e.what());
    }
    save_path = "";
  }
  continue_glline(win);
  if (draw_glfiledialog_button(win, "save image", scene_ok, "save image",
          save_path, true, get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    auto& app   = apps.get_selected();
    app.outname = save_path;
    try {
      save_image(app.imagename, app.display);
    } catch (std::exception& e) {
      push_glmessage("cannot save " + app.outname);
      log_glinfo(win, "cannot save " + app.outname);
      log_glinfo(win, e.what());
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
  draw_glcombobox(
      win, "scene", apps.selected, (int)apps.states.size(),
      [&apps](int idx) {
        auto it = apps.states.begin();
        std::advance(it, apps.selected);
        return it->name.c_str();
      },
      false);
  if (scene_ok && begin_glheader(win, "trace")) {
    auto  edited  = 0;
    auto& app     = apps.get_selected();
    auto& tparams = app.trace_prms;
    edited += draw_glcombobox(win, "camera", tparams.camera, app.scene.cameras);
    edited += draw_glslider(win, "resolution", tparams.resolution, 180, 4096);
    edited += draw_glslider(win, "nsamples", tparams.samples, 16, 4096);
    edited += draw_glcombobox(
        win, "tracer", (int&)tparams.sampler, trace_sampler_names);
    edited += draw_glcombobox(
        win, "false color", (int&)tparams.falsecolor, trace_falsecolor_names);
    edited += draw_glslider(win, "nbounces", tparams.bounces, 1, 128);
    edited += draw_glcheckbox(win, "env hidden", tparams.envhidden);
    continue_glline(win);
    edited += draw_glcheckbox(win, "filter", tparams.tentfilter);
    edited += draw_glslider(win, "seed", (int&)tparams.seed, 0, 1000000);
    edited += draw_glslider(win, "pratio", app.preview_ratio, 1, 64);
    auto& dparams = app.tonemap_prms;
    edited += draw_glslider(win, "exposure", dparams.exposure, -5, 5);
    draw_glcheckbox(win, "filmic", dparams.filmic);
    continue_glline(win);
    edited += draw_glcheckbox(win, "srgb", dparams.srgb);
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
    draw_glslider(win, "zoom", app.image_scale, 0.1, 10);
    draw_glcheckbox(win, "zoom to fit", app.zoom_to_fit);
    continue_glline(win);
    draw_glcheckbox(win, "fps", app.navigation_fps);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : app.scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : format_stats(app.scene)) print_info(stat);
      // for (auto stat : format_stats(app.bvh)) print_info(stat);
    }
    auto mouse_pos = get_glmouse_pos(win);
    auto ij        = get_image_coords(
        mouse_pos, app.image_center, app.image_scale, app.render.size());
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
        "camera", "shape", "environment", "instance", "materials", "textures"};
    auto& app = apps.get_selected();
    if (draw_glcombobox(win, "selection##1", app.selection.first, labels))
      app.selection.second = 0;
    auto edited = 0;
    if (app.selection.first == "camera") {
      edited += draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.cameras);
      edited += draw_glwidgets_camera(win, app, app.selection.second);
    } else if (app.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.textures);
      edited += draw_glwidgets_texture(win, app, app.selection.second);
    } else if (app.selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.materials);
      edited += draw_glwidgets_material(win, app, app.selection.second);
    } else if (app.selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.shapes);
      edited += draw_glwidgets_shape(win, app, app.selection.second);
    } else if (app.selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.instances);
      edited += draw_glwidgets_instance(win, app, app.selection.second);
    } else if (app.selection.first == "environment") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.environments);
      edited += draw_glwidgets_environment(win, app, app.selection.second);
    }
    if (edited) reset_display(app);
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

void draw(const opengl_window& win) {
  auto& apps     = *(app_states*)get_gluser_pointer(win);
  auto  win_size = get_glwindow_size(win);
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!apps.states.empty() && apps.selected >= 0) {
    auto& app = apps.get_selected();
    if (!app.gl_image || app.gl_image.size() != app.display.size())
      update_glimage(app.gl_image, app.display, false, false);
    update_imview(app.image_center, app.image_scale, app.display.size(),
        win_size, app.zoom_to_fit);
    draw_glimage(
        app.gl_image, win_size, app.image_center, app.image_scale, true);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_states& app) {
  while (!app.load_workers.empty() && is_ready(app.load_workers.front())) {
    try {
      app.load_workers.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    app.states.splice(app.states.end(), app.loading, app.loading.begin());
    app.load_workers.pop_front();
    reset_display(app.states.back());
    if (app.selected < 0) app.selected = (int)app.states.size() - 1;
  }
  for (auto& app : app.states) {
    if (app.render_preview) {
      // rendering preview
      auto preview_prms = app.trace_prms;
      preview_prms.resolution /= app.preview_ratio;
      preview_prms.samples = 1;
      auto preview = trace_image(app.scene, app.bvh, app.lights, preview_prms);
      preview      = tonemap_image(preview, app.tonemap_prms);
      for (auto j = 0; j < app.display.size().y; j++) {
        for (auto i = 0; i < app.display.size().x; i++) {
          auto pi = clamp(i / app.preview_ratio, 0, preview.size().x - 1),
               pj = clamp(j / app.preview_ratio, 0, preview.size().y - 1);
          app.display[{i, j}] = preview[{pi, pj}];
        }
      }
      if (!app.gl_image || app.gl_image.size() != app.display.size()) {
        update_glimage(app.gl_image, app.display, false, false);
      } else {
        update_glimage(app.gl_image, app.display, false, false);
      }
      app.render_preview = false;
    } else if (app.render_sample < app.trace_prms.samples) {
      // rendering blocks
      auto num_regions = min(
          128, app.render_regions.size() - app.render_region);
      parallel_for(app.render_region, app.render_region + num_regions,
          [&app](int region_id) {
            trace_region(app.render, app.state, app.scene, app.bvh, app.lights,
                app.render_regions[region_id], 1, app.trace_prms);
            tonemap_region(app.display, app.render,
                app.render_regions[region_id], app.tonemap_prms);
          });
      if (!app.gl_image || app.gl_image.size() != app.display.size()) {
        update_glimage(app.gl_image, app.display, false, false);
      } else {
        for (auto idx = 0; idx < num_regions; idx++)
          update_glimage_region(app.gl_image, app.display,
              app.render_regions[app.render_region + idx]);
      }
      app.render_region += num_regions;
      if (app.render_region >= app.render_regions.size()) {
        app.render_region = 0;
        app.render_sample += 1;
      }
    }
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
      auto& camera = app.scene.cameras.at(app.trace_prms.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      reset_display(app);
    }

    // selection
    if (scene_ok && (mouse_left || mouse_right) && alt_down &&
        !widgets_active) {
      auto& app = apps.get_selected();
      auto  ij  = get_image_coords(
          mouse_pos, app.image_center, app.image_scale, app.render.size());
      if (ij.x >= 0 && ij.x < app.render.size().x && ij.y >= 0 &&
          ij.y < app.render.size().y) {
        auto& camera = app.scene.cameras.at(app.trace_prms.camera);
        auto  ray    = eval_camera(
            camera, ij, app.render.size(), {0.5f, 0.5f}, zero2f);
        if (auto isec = intersect_scene_bvh(app.bvh, ray); isec.hit) {
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
  add_cli_option(cli, "--bvh-embree/--no-bvh-embree", app.bvh_prms.embree,
      "Use Embree ratracer");
  add_cli_option(cli, "--bvh-embree-compact/--no-bvh-embree-compact",
      app.bvh_prms.compact, "Embree runs in compact memory");
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
  for (auto filename : filenames) load_scene_async(app, filename);

  // run interactive
  run_ui(app);

  // done
  return 0;
}
