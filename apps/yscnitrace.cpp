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
using namespace yocto;

#include <map>
#include <list>

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};  // namespace yocto

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

  // scene
  yocto_scene scene      = {};
  bvh_scene   bvh        = {};
  bool        add_skyenv = false;

  // rendering state
  trace_lights lights  = {};
  trace_state  state   = {};
  image<vec4f> render  = {};
  image<vec4f> display = {};

  // view scene
  vec2f          image_center   = zero2f;
  float          image_scale    = 1;
  bool           zoom_to_fit    = true;
  bool           navigation_fps = false;
  opengl_texture gl_txt         = {};

  // editing
  pair<string, int> selection = {"camera", 0};

  // computation
  bool                 render_preview = true;
  int                  render_sample  = 0;
  int                  render_region  = 0;
  vector<image_region> render_regions = {};
};

// Application state
struct app_state {
  // data
  std::list<app_scene>    scenes;
  int                     selected = -1;
  std::list<string>       errors;
  std::list<app_scene>    loading;
  std::list<future<void>> load_workers;

  // get image
  app_scene& get_selected() {
    auto it = scenes.begin();
    std::advance(it, selected);
    return *it;
  }
  const app_scene& get_selected() const {
    auto it = scenes.begin();
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

void reset_display(app_scene& scene) {
  auto image_size = camera_resolution(
      scene.scene.cameras[scene.trace_prms.camera],
      scene.trace_prms.resolution);
  scene.render.resize(image_size);
  scene.display.resize(image_size);
  scene.render_preview = true;
  scene.render_sample  = 0;
  scene.render_region  = 0;
  scene.state = make_trace_state(scene.render.size(), scene.trace_prms.seed);
  scene.render_regions = make_regions(
      scene.render.size(), scene.trace_prms.region, true);
}

void load_scene_async(app_state& app, const string& filename) {
  auto& scene        = app.loading.emplace_back();
  scene.filename     = filename;
  scene.imagename    = replace_extension(filename, ".png");
  scene.outname      = replace_extension(filename, ".edited.yaml");
  scene.name         = get_filename(scene.filename);
  scene.load_prms    = app.load_prms;
  scene.save_prms    = app.save_prms;
  scene.trace_prms   = app.trace_prms;
  scene.bvh_prms     = app.bvh_prms;
  scene.tonemap_prms = app.tonemap_prms;
  scene.add_skyenv   = app.add_skyenv;
  app.load_workers.push_back(run_async([&scene]() {
    load_scene(scene.filename, scene.scene);
    make_bvh(scene.bvh, scene.scene, scene.bvh_prms);
    make_trace_lights(scene.lights, scene.scene);
    auto image_size = camera_resolution(
        scene.scene.cameras[scene.trace_prms.camera],
        scene.trace_prms.resolution);
    scene.render.resize(image_size);
    scene.display.resize(image_size);
    scene.name = get_filename(scene.filename) + " [" +
                 std::to_string(scene.render.size().x) + "x" +
                 std::to_string(scene.render.size().y) + " @ 0]";
  }));
}

bool draw_glwidgets_camera(const opengl_window& win, app_scene& scene, int id) {
  auto& camera = scene.scene.cameras[id];
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
bool draw_glwidgets_texture(
    const opengl_window& win, app_scene& scene, int id) {
  auto& texture      = scene.scene.textures[id];
  auto  old_filename = texture.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", texture.name);
  edited += draw_gltextinput(win, "filename", texture.filename);
  draw_gllabel(
      win, "hdr", "%d x %d", texture.hdr.size().x, texture.hdr.size().y);
  draw_gllabel(
      win, "ldr", "%d x %d", texture.ldr.size().x, texture.ldr.size().y);
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

bool draw_glwidgets_material(
    const opengl_window& win, app_scene& scene, int id) {
  auto& material = scene.scene.materials[id];
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
      win, "emission_tex", material.emission_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, scene.scene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, scene.scene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, scene.scene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_scene& scene, int id) {
  auto& shape        = scene.scene.shapes[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  draw_gllabel(win, "points", "%ld", shape.points.size());
  draw_gllabel(win, "lines", "%ld", shape.lines.size());
  draw_gllabel(win, "triangles", "%ld", shape.triangles.size());
  draw_gllabel(win, "quads", "%ld", shape.quads.size());
  draw_gllabel(win, "quads pos", "%ld", shape.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", shape.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", shape.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", shape.positions.size());
  draw_gllabel(win, "norm", "%ld", shape.normals.size());
  draw_gllabel(win, "texcoord", "%ld", shape.texcoords.size());
  draw_gllabel(win, "color", "%ld", shape.colors.size());
  draw_gllabel(win, "radius", "%ld", shape.radius.size());
  draw_gllabel(win, "tangsp", "%ld", shape.tangents.size());
  if (edited && old_filename != shape.filename) {
    try {
      load_shape(shape.filename, shape.points, shape.lines, shape.triangles,
          shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords, shape.colors,
          shape.radius, false);
    } catch (std::exception& e) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
      log_glinfo(win, e.what());
    }
    refit_bvh(scene.bvh, scene.scene, {}, {id}, scene.bvh_prms);
    // TODO: update lights
  }
  return edited;
}

inline bool draw_glwidgets_subdiv(
    const opengl_window& win, app_scene& scene, int id) {
  auto& shape        = scene.scene.subdivs[id];
  auto  old_filename = shape.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", shape.name);
  edited += draw_gltextinput(win, "filename", shape.filename);
  edited += draw_glslider(win, "subdivisions", shape.subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", shape.catmullclark);
  edited += draw_glcheckbox(win, "smooth", shape.smooth);
  edited += draw_glcheckbox(win, "facevarying", shape.facevarying);
  edited += draw_glcombobox(
      win, "shape", shape.shape, scene.scene.shapes, true);
  edited += draw_glcombobox(win, "displacement_tex", shape.displacement_tex,
      scene.scene.textures, true);
  edited += draw_glslider(win, "displacement", shape.displacement, 0, 1);
  draw_gllabel(win, "points", "%ld", shape.points.size());
  draw_gllabel(win, "lines", "%ld", shape.lines.size());
  draw_gllabel(win, "triangles", "%ld", shape.triangles.size());
  draw_gllabel(win, "quads", "%ld", shape.quads.size());
  draw_gllabel(win, "quads pos", "%ld", shape.quadspos.size());
  draw_gllabel(win, "quads norm", "%ld", shape.quadsnorm.size());
  draw_gllabel(win, "quads texcoord", "%ld", shape.quadstexcoord.size());
  draw_gllabel(win, "pos", "%ld", shape.positions.size());
  draw_gllabel(win, "norm", "%ld", shape.normals.size());
  draw_gllabel(win, "texcoord", "%ld", shape.texcoords.size());
  draw_gllabel(win, "color", "%ld", shape.colors.size());
  draw_gllabel(win, "radius", "%ld", shape.radius.size());
  if (edited && old_filename != shape.filename) {
    try {
      load_shape(shape.filename, shape.points, shape.lines, shape.triangles,
          shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords, shape.colors,
          shape.radius, false);
    } catch (std::exception& e) {
      push_glmessage("cannot load " + shape.filename);
      log_glinfo(win, "cannot load " + shape.filename);
      log_glinfo(win, e.what());
    }
    tesselate_subdiv(scene.scene, shape);
    refit_bvh(scene.bvh, scene.scene, {}, {id}, scene.bvh_prms);
    // TODO: update lights
  }
  if (edited && old_filename == shape.filename) {
    tesselate_subdiv(scene.scene, shape);
    refit_bvh(scene.bvh, scene.scene, {}, {shape.shape}, scene.bvh_prms);
    // TODO: update lights
  }
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window& win, app_scene& scene, int id) {
  auto& instance     = scene.scene.instances[id];
  auto  old_instance = instance;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame[0]", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", instance.shape, scene.scene.shapes, true);
  edited += draw_glcombobox(
      win, "material", instance.material, scene.scene.materials, true);
  if (edited && instance.shape != old_instance.shape)
    refit_bvh(scene.bvh, scene.scene, {}, {id}, scene.bvh_prms);
  if (edited && instance.frame != old_instance.frame)
    refit_bvh(scene.bvh, scene.scene, {}, {id}, scene.bvh_prms);
  // TODO: update lights
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, app_scene& scene, int id) {
  auto& environment = scene.scene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame[0]", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      scene.scene.textures, true);
  if (edited) {
    // TODO: update lights
  }
  return edited;
}

void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         app      = *(app_state*)get_gluser_pointer(win);
  auto          scene_ok = !app.scenes.empty() && app.selected >= 0;
  if (!begin_glwidgets_window(win, "yscnitrace")) return;
  draw_glmessages(win);
  if (draw_glfiledialog(
          win, "load", load_path, false, "./", "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(app, load_path);
  }
  if (draw_glfiledialog(win, "save", save_path, true, get_dirname(save_path),
          get_filename(save_path), "*.yaml;*.obj;*.pbrt")) {
    // app.scenes[app.selected].outname = save_path;
    // TODO: support save
    save_path = "";
  }
  if (draw_glfiledialog(win, "save image", save_path, true,
          get_dirname(save_path), get_filename(save_path),
          "*.png;*.jpg;*.tga;*.bmp;*.hdr;*.exr")) {
    // app.scenes[app.selected].imagename = save_path;
    // TODO: support save
    save_path = "";
  }
  if (draw_glbutton(win, "load")) {
    open_glmodal(win, "load");
  }
  continue_glline(win);
  if (draw_glbutton(win, "save", scene_ok)) {
    // save_path = app.scenes[app.selected].outname;
    // open_glmodal(win, "save");
    // TODO: support save
  }
  continue_glline(win);
  if (draw_glbutton(win, "save image", scene_ok)) {
    // save_path = app.scenes[app.selected].imagename;
    // open_glmodal(win, "save image");
    // TODO: support save
  }
  continue_glline(win);
  if (draw_glbutton(win, "close", scene_ok)) {
    // TODO: support close
  }
  continue_glline(win);
  if (draw_glbutton(win, "quit")) {
    set_glwindow_close(win, true);
  }
  if (app.scenes.empty()) return;
  draw_glcombobox(
      win, "scene", app.selected, (int)app.scenes.size(),
      [&app](int idx) { 
        auto it = app.scenes.begin();
        std::advance(it, app.selected);
        return it->name.c_str(); 
  }, false);
  if (scene_ok && begin_glheader(win, "trace")) {
    auto  edited  = 0;
    auto& scene   = app.get_selected();
    auto& tparams = scene.trace_prms;
    edited += draw_glcombobox(
        win, "camera", tparams.camera, scene.scene.cameras);
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
    edited += draw_glslider(win, "pratio", scene.preview_ratio, 1, 64);
    auto& dparams = scene.tonemap_prms;
    edited += draw_glslider(win, "exposure", dparams.exposure, -5, 5);
    draw_glcheckbox(win, "filmic", dparams.filmic);
    continue_glline(win);
    edited += draw_glcheckbox(win, "srgb", dparams.srgb);
    if (edited) reset_display(scene);
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "inspect")) {
    auto& scn = app.get_selected();
    draw_gllabel(win, "scene", get_filename(scn.filename));
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
      for (auto stat : format_stats(scn.scene)) print_info(stat);
      for (auto stat : format_stats(scn.bvh)) print_info(stat);
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
  if (scene_ok && begin_glheader(win, "edit")) {
    static auto labels = vector<string>{"camera", "shape", "environment",
        "instance", "materials", "textures", "subdivs"};
    auto&       scene  = app.get_selected();
    if (draw_glcombobox(win, "selection##1", scene.selection.first, labels))
      scene.selection.second = 0;
    auto edited = 0;
    if (scene.selection.first == "camera") {
      edited += draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.cameras);
      edited += draw_glwidgets_camera(win, scene, scene.selection.second);
    } else if (scene.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.textures);
      edited += draw_glwidgets_texture(win, scene, scene.selection.second);
    } else if (scene.selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.materials);
      edited += draw_glwidgets_material(win, scene, scene.selection.second);
    } else if (scene.selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.subdivs);
      edited += draw_glwidgets_shape(win, scene, scene.selection.second);
    } else if (scene.selection.first == "subdiv") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.subdivs);
      edited += draw_glwidgets_subdiv(win, scene, scene.selection.second);
    } else if (scene.selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", scene.selection.second, scene.scene.instances);
      edited += draw_glwidgets_instance(win, scene, scene.selection.second);
    } else if (scene.selection.first == "environment") {
      draw_glcombobox(win, "selection##2", scene.selection.second,
          scene.scene.environments);
      edited += draw_glwidgets_environment(win, scene, scene.selection.second);
    }
    if (edited) reset_display(scene);
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
  auto  fb_view  = get_glframebuffer_viewport(win);
  set_glviewport(fb_view);
  clear_glframebuffer(vec4f{0.15f, 0.15f, 0.15f, 1.0f});
  if (!app.scenes.empty() && app.selected >= 0) {
    auto& scene = app.get_selected();
    if (!scene.gl_txt || scene.gl_txt.size != scene.display.size())
      init_gltexture(scene.gl_txt, scene.display, false, false, false);
    update_imview(scene.image_center, scene.image_scale, scene.display.size(),
        win_size, scene.zoom_to_fit);
    draw_glimage_background(scene.gl_txt, win_size.x, win_size.y,
        scene.image_center, scene.image_scale);
    set_glblending(true);
    draw_glimage(scene.gl_txt, win_size.x, win_size.y, scene.image_center,
        scene.image_scale);
    set_glblending(false);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

void update(const opengl_window& win, app_state& app) {
  while(!app.load_workers.empty() && is_ready(app.load_workers.front())) {
    try {
      app.load_workers.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, "cannot load scene " + app.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    app.scenes.splice(app.scenes.end(), app.loading, app.loading.begin());
    reset_display(app.scenes.back());
    if(app.selected < 0) app.selected = (int)app.scenes.size() - 1;
  }
  for (auto& scene : app.scenes) {
    if (scene.render_preview) {
      // rendering preview
      auto preview_prms = scene.trace_prms;
      preview_prms.resolution /= scene.preview_ratio;
      preview_prms.samples = 1;
      auto preview         = trace_image(
          scene.scene, scene.bvh, scene.lights, preview_prms);
      preview = tonemap(preview, scene.tonemap_prms);
      for (auto j = 0; j < scene.display.size().y; j++) {
        for (auto i = 0; i < scene.display.size().x; i++) {
          auto pi = clamp(i / scene.preview_ratio, 0, preview.size().x - 1),
               pj = clamp(j / scene.preview_ratio, 0, preview.size().y - 1);
          scene.display[{i, j}] = preview[{pi, pj}];
        }
      }
      if (!scene.gl_txt || scene.gl_txt.size != scene.display.size()) {
        init_gltexture(scene.gl_txt, scene.display, false, false, false);
      } else {
        update_gltexture(scene.gl_txt, scene.display, false);
      }
      scene.render_preview = false;
    } else if (scene.render_sample < scene.trace_prms.samples) {
      // rendering blocks
      auto num_regions = min(
          128, scene.render_regions.size() - scene.render_region);
      parallel_for(scene.render_region, scene.render_region + num_regions,
          [&scene](int region_id) {
            trace_region(scene.render, scene.state, scene.scene, scene.bvh,
                scene.lights, scene.render_regions[region_id], 1,
                scene.trace_prms);
            tonemap(scene.display, scene.render,
                scene.render_regions[region_id], scene.tonemap_prms);
          });
      if (!scene.gl_txt || scene.gl_txt.size != scene.display.size()) {
        init_gltexture(scene.gl_txt, scene.display, false, false, false);
      } else {
        for(auto idx = 0; idx < num_regions; idx++)
        update_gltexture_region(scene.gl_txt, scene.display, 
          scene.render_regions[scene.render_region + idx], false);
      }
      scene.render_region += num_regions;
      if (scene.render_region >= scene.render_regions.size()) {
        scene.render_region = 0;
        scene.render_sample += 1;
      }
    }
  }
}

// run ui loop
void run_ui(app_state& app) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnitrace", &app, draw);
  set_drop_glcallback(win, [](const opengl_window& win, const vector<string>& paths){
      auto& app = *(app_state*)get_gluser_pointer(win);
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
    auto scene_ok = !app.scenes.empty() && app.selected >= 0;

    // handle mouse and keyboard for navigation
    if (scene_ok && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& scene  = app.get_selected();
      auto& camera = scene.scene.cameras.at(scene.trace_prms.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down)
        pan = (mouse_pos - last_pos) * camera.focus / 200.0f;
      pan.x = -pan.x;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      reset_display(scene);
    }

    // selection
    if (scene_ok && (mouse_left || mouse_right) && alt_down &&
        !widgets_active) {
      auto& scn = app.get_selected();
      auto  ij  = get_image_coords(
          mouse_pos, scn.image_center, scn.image_scale, scn.render.size());
      if (ij.x >= 0 && ij.x < scn.render.size().x && ij.y >= 0 &&
          ij.y < scn.render.size().y) {
        auto& camera = scn.scene.cameras.at(scn.trace_prms.camera);
        auto  ray    = eval_camera(
            camera, ij, scn.render.size(), {0.5f, 0.5f}, zero2f);
        if (auto isec = intersect_bvh(scn.bvh, ray); isec.hit) {
          scn.selection = {"instance", isec.instance};
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
  for (auto filename : filenames) load_scene_async(app, filename);

  // run interactive
  run_ui(app);

  // done
  return 0;
}
