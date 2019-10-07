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
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_scene.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <list>

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(const yocto_camera& camera);
};

// Application state
struct app_state {
  // loading parameters
  string filename  = "scene.json";
  string imagename = "out.png";
  string outname   = "scene.json";
  string name      = "";

  // options
  load_params         load_prms   = {};
  save_params         save_prms   = {};
  draw_glscene_params drawgl_prms = {};

  // scene
  yocto_scene scene = {};

  // rendering state
  opengl_scene state = {};

  // view image
  bool   navigation_fps = false;
  float  time           = 0;
  string anim_group     = "";
  vec2f  time_range     = zero2f;
  bool   animate        = false;

  // editing
  pair<string, int> selection = {"camera", 0};
};

// Application state
struct app_states {
  // data
  std::list<app_state>         states;
  int                          selected = -1;
  std::list<string>            errors;
  std::list<app_state>         loading;
  std::list<std::future<void>> load_workers;

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
  load_params         load_prms   = {};
  save_params         save_prms   = {};
  draw_glscene_params drawgl_prms = {};
};

void load_scene_async(app_states& apps, const string& filename) {
  auto& app       = apps.loading.emplace_back();
  app.filename    = filename;
  app.imagename   = replace_extension(filename, ".png");
  app.outname     = replace_extension(filename, ".edited.yaml");
  app.name        = get_filename(app.filename);
  app.load_prms   = app.load_prms;
  app.save_prms   = app.save_prms;
  app.drawgl_prms = app.drawgl_prms;
  apps.load_workers.push_back(run_async([&app]() {
    load_scene(app.filename, app.scene);
    update_tesselation(app.scene);
    app.time_range = compute_animation_range(app.scene);
    app.time       = app.time_range.x;
  }));
}

void update_glcamera(opengl_camera& glcamera, const yocto_camera& camera) {
  glcamera.frame  = camera.frame;
  glcamera.yfov   = camera_yfov(camera);
  glcamera.asepct = camera_aspect(camera);
  glcamera.near   = 0.001f;
  glcamera.far    = 10000;
}

void update_gltexture(opengl_texture& gltexture, const yocto_texture& texture) {
  if (!texture.hdr.empty()) {
    init_gltexture(gltexture, texture.hdr, true, true, true);
  } else if (!texture.ldr.empty()) {
    init_gltexture(gltexture, texture.ldr, true, true, true);
  } else {
    throw std::runtime_error("bad texture");
  }
}

void update_glmaterial(
    opengl_material& glmaterial, const yocto_material& material) {
  glmaterial.emission     = material.emission;
  glmaterial.diffuse      = material.diffuse;
  glmaterial.specular     = material.specular;
  glmaterial.metallic     = material.metallic;
  glmaterial.roughness    = material.roughness;
  glmaterial.opacity      = material.opacity;
  glmaterial.emission_map = material.emission_tex;
  glmaterial.diffuse_map  = material.diffuse_tex;
  glmaterial.specular_map = material.specular_tex;
  glmaterial.metallic_map = material.metallic_tex;
  glmaterial.normal_map   = material.normal_tex;
}

void update_glshape(opengl_shape& glshape, const yocto_shape& shape) {
  if (shape.quadspos.empty()) {
    if (!shape.positions.empty())
      init_glarraybuffer(glshape.positions, shape.positions, false);
    if (!shape.normals.empty())
      init_glarraybuffer(glshape.normals, shape.normals, false);
    if (!shape.texcoords.empty())
      init_glarraybuffer(glshape.texcoords, shape.texcoords, false);
    if (!shape.colors.empty())
      init_glarraybuffer(glshape.colors, shape.colors, false);
    if (!shape.tangents.empty())
      init_glarraybuffer(glshape.tangentsps, shape.tangents, false);
    if (!shape.points.empty())
      init_glelementbuffer(glshape.points, shape.points, false);
    if (!shape.lines.empty())
      init_glelementbuffer(glshape.lines, shape.lines, false);
    if (!shape.triangles.empty())
      init_glelementbuffer(glshape.triangles, shape.triangles, false);
    if (!shape.quads.empty()) {
      auto triangles = quads_to_triangles(shape.quads);
      init_glelementbuffer(glshape.quads, triangles, false);
    }
  } else {
    auto quads     = vector<vec4i>{};
    auto positions = vector<vec3f>{};
    auto normals   = vector<vec3f>{};
    auto texcoords = vector<vec2f>{};
    split_facevarying(quads, positions, normals, texcoords, shape.quadspos,
        shape.quadsnorm, shape.quadstexcoord, shape.positions, shape.normals,
        shape.texcoords);
    if (!positions.empty())
      init_glarraybuffer(glshape.positions, positions, false);
    if (!normals.empty()) init_glarraybuffer(glshape.normals, normals, false);
    if (!texcoords.empty())
      init_glarraybuffer(glshape.texcoords, texcoords, false);
    if (!quads.empty()) {
      auto triangles = quads_to_triangles(quads);
      init_glelementbuffer(glshape.quads, triangles, false);
    }
  }
}

void update_glinstance(
    opengl_instance& glinstance, const yocto_instance& instance) {
  glinstance.frame    = instance.frame;
  glinstance.shape    = instance.shape;
  glinstance.material = instance.material;
}

void update_gllights(opengl_scene& state, const yocto_scene& scene) {
  state.lights = {};
  for (auto& instance : scene.instances) {
    if (state.lights.size() >= 16) break;
    if (instance.shape < 0) continue;
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    if (material.emission == zero3f) continue;
    auto bbox = compute_bounds(shape);
    auto pos  = (bbox.max + bbox.min) / 2;
    auto area = 0.0f;
    if (!shape.triangles.empty()) {
      for (auto t : shape.triangles)
        area += triangle_area(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
    } else if (!shape.quads.empty()) {
      for (auto q : shape.quads)
        area += quad_area(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w]);
    } else if (!shape.lines.empty()) {
      for (auto l : shape.lines)
        area += line_length(shape.positions[l.x], shape.positions[l.y]);
    } else {
      area += shape.positions.size();
    }
    auto  ke       = material.emission * area;
    auto& light    = state.lights.emplace_back();
    light.position = transform_point(instance.frame, pos);
    light.emission = ke;
    light.type     = 0;
  }
}

void make_glscene(opengl_scene& glscene, const yocto_scene& scene) {
  // load program
  make_glscene(glscene);

  // camera
  for (auto& camera : scene.cameras) {
    update_glcamera(glscene.cameras.emplace_back(), camera);
  }

  // textures
  for (auto& texture : scene.textures) {
    update_gltexture(glscene.textures.emplace_back(), texture);
  }

  // materials
  for (auto& material : scene.materials) {
    update_glmaterial(glscene.materials.emplace_back(), material);
  }

  // shapes
  for (auto& shape : scene.shapes) {
    update_glshape(glscene.shapes.emplace_back(), shape);
  }

  // instances
  for (auto& instance : scene.instances) {
    update_glinstance(glscene.instances.emplace_back(), instance);
  }
}

bool draw_glwidgets_camera(const opengl_window& win, app_state& scene, int id) {
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
  if (edited) update_glcamera(scene.state.cameras[id], camera);
  return edited;
}

/// Visit struct elements.
bool draw_glwidgets_texture(
    const opengl_window& win, app_state& scene, int id) {
  auto& texture      = scene.scene.textures[id];
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
  }
  if (edited) update_gltexture(scene.state.textures[id], texture);
  return edited;
}

bool draw_glwidgets_material(
    const opengl_window& win, app_state& scene, int id) {
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
  if (edited) update_glmaterial(scene.state.materials[id], material);
  return edited;
}

bool draw_glwidgets_shape(const opengl_window& win, app_state& scene, int id) {
  auto& shape        = scene.scene.shapes[id];
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
  edited += draw_glslider(win, "subdivisions", shape.subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", shape.catmullclark);
  edited += draw_glcheckbox(win, "smooth", shape.smooth);
  edited += draw_glcheckbox(win, "facevarying", shape.facevarying);
  edited += draw_glcombobox(win, "displacement_tex", shape.displacement_tex,
      scene.scene.textures, true);
  edited += draw_glslider(win, "displacement", shape.displacement, 0, 1);
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
    // TODO: update lights
  }
  if (edited) update_glshape(scene.state.shapes[id], shape);
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window& win, app_state& scene, int id) {
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
  // TODO: update lights
  if (edited) update_glinstance(scene.state.instances[id], instance);
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, app_state& scene, int id) {
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

// draw with shading
void draw_glwidgets(const opengl_window& win) {
  static string load_path = "", save_path = "", error_message = "";
  auto&         apps     = *(app_states*)get_gluser_pointer(win);
  auto          scene_ok = !apps.states.empty() && apps.selected >= 0;
  if (!begin_glwidgets_window(win, "yscnview")) return;
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
  if (apps.states.empty()) return;
  draw_glcombobox(
      win, "scene", apps.selected, (int)apps.states.size(),
      [&apps](int idx) {
        auto it = apps.states.begin();
        std::advance(it, apps.selected);
        return it->name.c_str();
      },
      false);
  if (scene_ok && begin_glheader(win, "view")) {
    auto& app    = apps.get_selected();
    auto& params = app.drawgl_prms;
    draw_glcombobox(win, "camera", params.camera, app.scene.cameras);
    draw_glslider(win, "resolution", params.resolution, 0, 4096);
    draw_glcheckbox(win, "eyelight", params.eyelight);
    continue_glline(win);
    draw_glcheckbox(win, "wireframe", params.wireframe);
    continue_glline(win);
    draw_glcheckbox(win, "edges", params.edges);
    if (app.time_range != zero2f) {
      draw_glslider(win, "time", app.time, app.time_range.x, app.time_range.y);
      draw_gltextinput(win, "anim group", app.anim_group);
      draw_glcheckbox(win, "animate", app.animate);
    }
    draw_glslider(win, "exposure", params.exposure, -10, 10);
    draw_glslider(win, "gamma", params.gamma, 0.1f, 4);
    draw_glcheckbox(win, "double sided", params.double_sided);
    draw_glslider(win, "near", params.near, 0.01f, 1.0f);
    draw_glslider(win, "far", params.far, 1000.0f, 10000.0f);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    auto& app = apps.get_selected();
    draw_gllabel(win, "scene", get_filename(app.filename));
    draw_gllabel(win, "filename", app.filename);
    draw_gllabel(win, "outname", app.outname);
    draw_gllabel(win, "imagename", app.imagename);
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
          win, "selection##2", app.selection.second, app.scene.cameras);
      draw_glwidgets_camera(win, app, app.selection.second);
    } else if (app.selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.textures);
      draw_glwidgets_texture(win, app, app.selection.second);
    } else if (app.selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.materials);
      draw_glwidgets_material(win, app, app.selection.second);
    } else if (app.selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.shapes);
      draw_glwidgets_shape(win, app, app.selection.second);
    } else if (app.selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.instances);
      draw_glwidgets_instance(win, app, app.selection.second);
    } else if (app.selection.first == "environment") {
      draw_glcombobox(
          win, "selection##2", app.selection.second, app.scene.environments);
      draw_glwidgets_environment(win, app, app.selection.second);
    }
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(const opengl_window& win) {
  auto& apps = *(app_states*)get_gluser_pointer(win);

  if (!apps.states.empty() && apps.selected >= 0) {
    auto& app = apps.get_selected();
    draw_glscene(app.state, get_glframebuffer_viewport(win), app.drawgl_prms);
  }
  begin_glwidgets(win);
  draw_glwidgets(win);
  end_glwidgets(win);
  swap_glbuffers(win);
}

// update
void update(const opengl_window& win, app_states& apps) {
  while (!apps.load_workers.empty() && is_ready(apps.load_workers.front())) {
    try {
      apps.load_workers.front().get();
    } catch (const std::exception& e) {
      push_glmessage(win, "cannot load scene " + apps.loading.front().filename);
      log_glinfo(win, "cannot load scene " + apps.loading.front().filename);
      log_glinfo(win, e.what());
      break;
    }
    apps.states.splice(apps.states.end(), apps.loading, apps.loading.begin());
    apps.load_workers.pop_front();
    make_glscene(apps.states.back().state, apps.states.back().scene);
    update_gllights(apps.states.back().state, apps.states.back().scene);
    if (apps.selected < 0) apps.selected = (int)apps.states.size() - 1;
  }
}

// run ui loop
void run_ui(app_states& apps) {
  // window
  auto win = opengl_window();
  init_glwindow(win, {1280 + 320, 720}, "yscnview", &apps, draw);
  set_drop_glcallback(
      win, [](const opengl_window& win, const vector<string>& paths) {
        auto& app = *(app_states*)get_gluser_pointer(win);
        for (auto& path : paths) load_scene_async(app, path);
      });

  // init widget
  init_glwidgets(win);

  // loop
  auto mouse_pos = zero2f, last_pos = zero2f;
  auto last_time = std::chrono::high_resolution_clock::now();
  while (!should_glwindow_close(win)) {
    last_pos            = mouse_pos;
    mouse_pos           = get_glmouse_pos(win);
    auto mouse_left     = get_glmouse_left(win);
    auto mouse_right    = get_glmouse_right(win);
    auto alt_down       = get_glalt_key(win);
    auto shift_down     = get_glshift_key(win);
    auto widgets_active = get_glwidgets_active(win);
    auto scene_ok       = !apps.states.empty() && apps.selected >= 0;

    // update trasforms
    if (scene_ok) {
      auto& app = apps.get_selected();
      update_transforms(app.scene, app.time);
    }

    // handle mouse and keyboard for navigation
    if (scene_ok && (mouse_left || mouse_right) && !alt_down &&
        !widgets_active) {
      auto& app    = apps.get_selected();
      auto& camera = app.scene.cameras.at(app.drawgl_prms.camera);
      auto  dolly  = 0.0f;
      auto  pan    = zero2f;
      auto  rotate = zero2f;
      if (mouse_left && !shift_down) rotate = (mouse_pos - last_pos) / 100.0f;
      if (mouse_right) dolly = (mouse_pos.x - last_pos.x) / 100.0f;
      if (mouse_left && shift_down) pan = (mouse_pos - last_pos) / 100.0f;
      update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
      update_glcamera(app.state.cameras[app.drawgl_prms.camera], camera);
    }

    // animation
    if (scene_ok && apps.get_selected().animate) {
      auto& app     = apps.get_selected();
      auto  now     = std::chrono::high_resolution_clock::now();
      auto  elapsed = now - last_time;
      auto  time    = (double)(elapsed.count()) / 1000000000.0;
      app.time += min(1 / 60.0f, (float)time);
      if (app.time < app.time_range.x || app.time > app.time_range.y)
        app.time = app.time_range.x;
      update_transforms(app.scene, app.time);
      last_time = now;
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
  // initialize app
  app_states app{};
  auto       filenames  = vector<string>{};
  auto       noparallel = false;

  // parse command line
  auto cli = make_cli("yscnview", "views scenes inteactively");
  add_cli_option(cli, "--camera", app.drawgl_prms.camera, "Camera index.");
  add_cli_option(
      cli, "--resolution,-r", app.drawgl_prms.resolution, "Image resolution.");
  add_cli_option(cli, "--eyelight/--no-eyelight,-c", app.drawgl_prms.eyelight,
      "Eyelight rendering.");
  add_cli_option(
      cli, "--noparallel", noparallel, "Disable parallel execution.");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // fix parallel code
  if (noparallel) {
    app.load_prms.noparallel = true;
    app.save_prms.noparallel = true;
  }

  // loading images
  for (auto filename : filenames) load_scene_async(app, filename);

  // run ui
  run_ui(app);

  // done
  return 0;
}
