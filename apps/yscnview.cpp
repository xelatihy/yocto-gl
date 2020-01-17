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
#include "../yocto/yocto_image.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <deque>
#include <future>
#include <memory>
using namespace std;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(const sceneio_camera& camera);
};

// Application state
struct app_state {
  // loading parameters
  string filename  = "scene.json";
  string imagename = "out.png";
  string outname   = "scene.json";
  string name      = "";

  // options
  draw_glscene_params drawgl_prms = {};

  // scene
  sceneio_model scene = {};

  // rendering state
  unique_ptr<opengl_scene> glscene = {};

  // view image
  float  time       = 0;
  string anim_group = "";
  vec2f  time_range = zero2f;
  bool   animate    = false;

  // editing
  pair<string, int> selection = {"camera", 0};

  // error
  string error = "";
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
  draw_glscene_params drawgl_prms = {};
};

// Compute animation range
vec2f compute_animation_range(
    const sceneio_model& scene, const string& anim_group = "") {
  if (scene.animations.empty()) return zero2f;
  auto range = vec2f{+flt_max, -flt_max};
  for (auto& animation : scene.animations) {
    if (anim_group != "" && animation.group != anim_group) continue;
    range.x = min(range.x, animation.times.front());
    range.y = max(range.y, animation.times.back());
  }
  if (range.y < range.x) return zero2f;
  return range;
}

void load_scene_async(shared_ptr<app_states> apps, const string& filename) {
  apps->loaders.push_back(
      async(launch::async, [apps, filename]() -> load_state {
        auto app         = make_shared<app_state>();
        app->filename    = filename;
        app->imagename   = replace_extension(filename, ".png");
        app->outname     = replace_extension(filename, ".edited.yaml");
        app->name        = get_filename(app->filename);
        app->drawgl_prms = apps->drawgl_prms;
        if (auto ret = load_scene(app->filename, app->scene); !ret)
          return {filename, nullptr, ret};
        app->time_range = compute_animation_range(app->scene);
        app->time       = app->time_range.x;
        return {filename, app, {}};
      }));
}

void update_lights(opengl_scene* glscene, const sceneio_model& scene) {
  clear_lights(glscene);
  for (auto& instance : scene.instances) {
    if (has_max_lights(glscene)) break;
    if (instance.shape < 0) continue;
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    if (material.emission == zero3f) continue;
    auto bbox = invalidb3f;
    for (auto& p : shape.positions) bbox = merge(bbox, p);
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
    auto ke = material.emission * area;
    add_light(glscene, transform_point(instance.frame, pos), ke, false);
  }
}

opengl_scene* make_scene(sceneio_model& scene) {
  // load program
  auto glscene = make_glscene();

  // camera
  for (auto& camera : scene.cameras) {
    add_camera(glscene, camera.frame, camera.lens, camera.aspect, camera.film,
        0.001, 10000);
  }

  // textures
  for (auto& texture : scene.textures) {
    if (!texture.hdr.empty()) {
      add_texture(glscene, texture.hdr);
    } else if (!texture.ldr.empty()) {
      add_texture(glscene, texture.ldr);
    }
  }

  // materials
  for (auto& material : scene.materials) {
    auto id = add_material(glscene);
    set_material_emission(
        glscene, id, material.emission, material.emission_tex);
    set_material_diffuse(glscene, id, material.diffuse, material.diffuse_tex);
    set_material_specular(
        glscene, id, material.specular, material.specular_tex);
    set_material_metallic(
        glscene, id, material.metallic, material.metallic_tex);
    set_material_roughness(
        glscene, id, material.roughness, material.roughness_tex);
    set_material_opacity(glscene, id, material.opacity, material.opacity_tex);
    set_material_normalmap(glscene, id, material.normal_tex);
  }

  for (auto& subdiv : scene.subdivs) {
    tesselate_subdiv(scene, subdiv);
  }

  // shapes
  for (auto& shape : scene.shapes) {
    auto id = add_shape(glscene);
    if (!shape.points.empty()) {
      set_shape(glscene, id, shape.points, shape.positions, shape.normals,
          shape.texcoords, shape.colors);
    } else if (!shape.lines.empty()) {
      set_shape(glscene, id, shape.lines, shape.positions, shape.normals,
          shape.texcoords, shape.colors);
    } else if (!shape.triangles.empty()) {
      set_shape(glscene, id, shape.triangles, shape.positions, shape.normals,
          shape.texcoords, shape.colors, shape.tangents);
    } else if (!shape.quads.empty()) {
      set_shape(glscene, id, shape.quads, shape.positions, shape.normals,
          shape.texcoords, shape.colors, shape.tangents);
    } else if (!shape.quadspos.empty()) {
      auto [quads, positions, normals, texcoords] = split_facevarying(
          shape.quadspos, shape.quadsnorm, shape.quadstexcoord, shape.positions,
          shape.normals, shape.texcoords);
      set_shape(glscene, id, quads, positions, normals, texcoords);
    }
  }

  // instances
  for (auto& instance : scene.instances) {
    add_instance(glscene, instance.frame, instance.shape, instance.material);
  }

  return glscene;
}

bool draw_glwidgets_camera(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& camera = app->scene.cameras[id];
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
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& texture      = app->scene.textures[id];
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
  auto& material = app->scene.materials[id];
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
      win, "emission_tex", material.emission_tex, app->scene.textures, true);
  edited += draw_glcombobox(
      win, "diffuse_tex", material.diffuse_tex, app->scene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", material.metallic_tex, app->scene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", material.specular_tex, app->scene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", material.transmission_tex,
      app->scene.textures, true);
  edited += draw_glcombobox(win, "subsurface_tex", material.subsurface_tex,
      app->scene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", material.roughness_tex, app->scene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", material.normal_tex, app->scene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", material.gltf_textures);
  return edited;
}

bool draw_glwidgets_shape(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& shape        = app->scene.shapes[id];
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
  edited += draw_glslider(win, "subdivisions", shape.subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", shape.catmullclark);
  edited += draw_glcheckbox(win, "smooth", shape.smooth);
  edited += draw_glcheckbox(win, "facevarying", shape.facevarying);
  edited += draw_glcombobox(win, "displacement_tex", shape.displacement_tex,
      app->scene.textures, true);
  edited += draw_glslider(win, "displacement", shape.displacement, 0, 1);
  if (edited && old_filename != shape.filename) {
    if (auto ret = load_shape(app->filename, shape); !ret) {
      push_glmessage(win, ret.error);
      log_glinfo(win, ret.error);
    }
  }
  return edited;
}

bool draw_glwidgets_subdiv(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& subdiv       = app->scene.subdivs[id];
  auto  old_filename = subdiv.filename;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", subdiv.name);
  edited += draw_gltextinput(win, "filename", subdiv.filename);
  draw_gllabel(win, "points", to_string(subdiv.points.size()));
  draw_gllabel(win, "lines", to_string(subdiv.lines.size()));
  draw_gllabel(win, "triangles", to_string(subdiv.triangles.size()));
  draw_gllabel(win, "quads", to_string(subdiv.quads.size()));
  draw_gllabel(win, "quads pos", to_string(subdiv.quadspos.size()));
  draw_gllabel(win, "quads norm", to_string(subdiv.quadsnorm.size()));
  draw_gllabel(win, "quads texcoord", to_string(subdiv.quadstexcoord.size()));
  draw_gllabel(win, "pos", to_string(subdiv.positions.size()));
  draw_gllabel(win, "norm", to_string(subdiv.normals.size()));
  draw_gllabel(win, "texcoord", to_string(subdiv.texcoords.size()));
  draw_gllabel(win, "color", to_string(subdiv.colors.size()));
  draw_gllabel(win, "radius", to_string(subdiv.radius.size()));
  draw_gllabel(win, "tangsp", to_string(subdiv.tangents.size()));
  edited += draw_glslider(win, "subdivisions", subdiv.subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", subdiv.catmullclark);
  edited += draw_glcheckbox(win, "smooth", subdiv.smooth);
  edited += draw_glcheckbox(win, "facevarying", subdiv.facevarying);
  edited += draw_glcombobox(win, "displacement_tex", subdiv.displacement_tex,
      app->scene.textures, true);
  edited += draw_glslider(win, "displacement", subdiv.displacement, 0, 1);
  if (edited && old_filename != subdiv.filename) {
    if (auto ret = load_subdiv(app->filename, subdiv); !ret) {
      push_glmessage(win, ret.error);
      log_glinfo(win, ret.error);
    }
  }
  return edited;
}

bool draw_glwidgets_instance(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& instance     = app->scene.instances[id];
  auto  old_instance = instance;
  auto  edited       = 0;
  edited += draw_gltextinput(win, "name", instance.name);
  edited += draw_glslider(win, "frame[0]", instance.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", instance.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", instance.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", instance.frame.o, -10, 10);
  edited += draw_glcombobox(
      win, "shape", instance.shape, app->scene.shapes, true);
  edited += draw_glcombobox(
      win, "material", instance.material, app->scene.materials, true);
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window* win, shared_ptr<app_state> app, int id) {
  auto& environment = app->scene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", environment.name);
  edited += draw_glslider(win, "frame[0]", environment.frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", environment.frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", environment.frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", environment.frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", environment.emission);
  edited += draw_glcombobox(win, "emission texture", environment.emission_tex,
      app->scene.textures, true);
  return edited;
}

// draw with shading
void draw_glwidgets(const opengl_window* win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  static auto load_path = ""s, save_path = ""s, error_message = ""s;
  auto        scene_ok = !apps->states.empty() && apps->selected >= 0;
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
    if (auto ret = save_scene(app->outname, app->scene); !ret) {
      push_glmessage(win, "cannot save " + app->outname);
      log_glinfo(win, "cannot save " + app->outname);
      log_glinfo(win, ret.error);
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
  if (apps->states.empty()) return;
  draw_glcombobox(
      win, "scene", apps->selected, (int)apps->states.size(),
      [apps](int idx) { return apps->states[idx]->name.c_str(); }, false);
  if (scene_ok && begin_glheader(win, "view")) {
    auto  app    = apps->states[apps->selected];
    auto& params = app->drawgl_prms;
    draw_glcombobox(win, "camera", params.camera, app->scene.cameras);
    draw_glslider(win, "resolution", params.resolution, 0, 4096);
    draw_glcheckbox(win, "eyelight", params.eyelight);
    continue_glline(win);
    draw_glcheckbox(win, "wireframe", params.wireframe);
    continue_glline(win);
    draw_glcheckbox(win, "edges", params.edges);
    if (app->time_range != zero2f) {
      draw_glslider(
          win, "time", app->time, app->time_range.x, app->time_range.y);
      draw_gltextinput(win, "anim group", app->anim_group);
      draw_glcheckbox(win, "animate", app->animate);
    }
    draw_glslider(win, "exposure", params.exposure, -10, 10);
    draw_glslider(win, "gamma", params.gamma, 0.1f, 4);
    draw_glcheckbox(win, "double sided", params.double_sided);
    draw_glslider(win, "near", params.near, 0.01f, 1.0f);
    draw_glslider(win, "far", params.far, 1000.0f, 10000.0f);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
    auto app = apps->states[apps->selected];
    draw_gllabel(win, "scene", get_filename(app->filename));
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "imagename", app->imagename);
    continue_glline(win);
    if (draw_glbutton(win, "print cams")) {
      for (auto& camera : app->scene.cameras) {
        print_obj_camera(camera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : scene_stats(app->scene)) print_info(stat);
    }
    end_glheader(win);
  }
  if (scene_ok && begin_glheader(win, "edit")) {
    static auto labels  = vector<string>{"camera", "shape", "subdiv",
        "environment", "instance", "material", "texture"};
    auto        app     = apps->states[apps->selected];
    auto        glscene = app->glscene.get();
    if (draw_glcombobox(win, "selection##1", app->selection.first, labels))
      app->selection.second = 0;
    if (app->selection.first == "camera") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.cameras);
      if (draw_glwidgets_camera(win, app, app->selection.second)) {
        auto& camera = app->scene.cameras[app->selection.second];
        set_camera(glscene, app->selection.second, camera.frame, camera.lens,
            camera.aspect, camera.film, 0.001, 10000);
      }
    } else if (app->selection.first == "texture") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.textures);
      if (draw_glwidgets_texture(win, app, app->selection.second)) {
        auto& texture = app->scene.textures[app->selection.second];
        if (!texture.hdr.empty()) {
          set_texture(glscene, app->selection.second, texture.hdr);
        } else if (!texture.hdr.empty()) {
          set_texture(glscene, app->selection.second, texture.ldr);
        }
      }
    } else if (app->selection.first == "material") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.materials);
      if (draw_glwidgets_material(win, app, app->selection.second)) {
        auto& material = app->scene.materials[app->selection.second];
        set_material_emission(glscene, app->selection.second, material.emission,
            material.emission_tex);
        set_material_diffuse(glscene, app->selection.second, material.diffuse,
            material.diffuse_tex);
        set_material_specular(glscene, app->selection.second, material.specular,
            material.specular_tex);
        set_material_metallic(glscene, app->selection.second, material.metallic,
            material.metallic_tex);
        set_material_roughness(glscene, app->selection.second,
            material.roughness, material.roughness_tex);
        set_material_opacity(glscene, app->selection.second, material.opacity,
            material.opacity_tex);
        set_material_normalmap(
            glscene, app->selection.second, material.normal_tex);
      }
    } else if (app->selection.first == "shape") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.shapes);
      if (!draw_glwidgets_shape(win, app, app->selection.second)) {
        auto& shape = app->scene.shapes[app->selection.second];
        auto  idx   = app->selection.second;
        if (!shape.points.empty()) {
          set_shape(glscene, idx, shape.points, shape.positions, shape.normals,
              shape.texcoords, shape.colors);
        } else if (!shape.lines.empty()) {
          set_shape(glscene, idx, shape.lines, shape.positions, shape.normals,
              shape.texcoords, shape.colors);
        } else if (!shape.triangles.empty()) {
          set_shape(glscene, idx, shape.triangles, shape.positions,
              shape.normals, shape.texcoords, shape.colors, shape.tangents);
        } else if (!shape.quads.empty()) {
          set_shape(glscene, idx, shape.quads, shape.positions, shape.normals,
              shape.texcoords, shape.colors, shape.tangents);
        } else if (!shape.quadspos.empty()) {
          auto [quads, positions, normals, texcoords] = split_facevarying(
              shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
              shape.positions, shape.normals, shape.texcoords);
          set_shape(glscene, idx, quads, positions, normals, texcoords);
        }
      }
    } else if (app->selection.first == "subdiv") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.subdivs);
      if (!draw_glwidgets_subdiv(win, app, app->selection.second)) {
        auto& subdiv = app->scene.subdivs[app->selection.second];
        tesselate_subdiv(app->scene, subdiv);
        auto& shape = app->scene.shapes[subdiv.shape];
        auto  idx   = app->selection.second;
        if (!shape.points.empty()) {
          set_shape(glscene, idx, shape.points, shape.positions, shape.normals,
              shape.texcoords, shape.colors);
        } else if (!shape.lines.empty()) {
          set_shape(glscene, idx, shape.lines, shape.positions, shape.normals,
              shape.texcoords, shape.colors);
        } else if (!shape.triangles.empty()) {
          set_shape(glscene, idx, shape.triangles, shape.positions,
              shape.normals, shape.texcoords, shape.colors, shape.tangents);
        } else if (!shape.quads.empty()) {
          set_shape(glscene, idx, shape.quads, shape.positions, shape.normals,
              shape.texcoords, shape.colors, shape.tangents);
        } else if (!shape.quadspos.empty()) {
          auto [quads, positions, normals, texcoords] = split_facevarying(
              shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
              shape.positions, shape.normals, shape.texcoords);
          set_shape(glscene, idx, quads, positions, normals, texcoords);
        }
      }
    } else if (app->selection.first == "instance") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.instances);
      if (draw_glwidgets_instance(win, app, app->selection.second)) {
        auto& instance = app->scene.instances[app->selection.second];
        set_instance(glscene, app->selection.second, instance.frame,
            instance.shape, instance.material);
      }
    } else if (app->selection.first == "environment") {
      draw_glcombobox(
          win, "selection##2", app->selection.second, app->scene.environments);
      draw_glwidgets_environment(win, app, app->selection.second);
    }
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(const opengl_window* win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->states.empty() && apps->selected >= 0) {
    auto app = apps->states[apps->selected];
    draw_glscene(
        app->glscene.get(), input.framebuffer_viewport, app->drawgl_prms);
  }
}

// update
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
      break;
    } else {
      apps->states.push_back(app);
      app->glscene = unique_ptr<opengl_scene>{make_scene(app->scene)};
      update_lights(app->glscene.get(), app->scene);
      if (apps->selected < 0) apps->selected = (int)apps->states.size() - 1;
    }
  }
}

int main(int argc, const char* argv[]) {
  // initialize app
  auto apps       = make_shared<app_states>();
  auto filenames  = vector<string>{};
  auto noparallel = false;

  // parse command line
  auto cli = make_cli("yscnview", "views scenes inteactively");
  add_cli_option(cli, "--camera", apps->drawgl_prms.camera, "Camera index.");
  add_cli_option(cli, "--resolution,-r", apps->drawgl_prms.resolution,
      "Image resolution.");
  add_cli_option(cli, "--eyelight/--no-eyelight,-c", apps->drawgl_prms.eyelight,
      "Eyelight rendering.");
  add_cli_option(
      cli, "--noparallel", noparallel, "Disable parallel execution.");
  add_cli_option(cli, "scenes", filenames, "Scene filenames", true);
  if (!parse_cli(cli, argc, argv)) exit(1);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  auto win = make_glwindow({1280 + 320, 720}, "yscnview", true);

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
  set_uiupdate_glcallback(
      win, [apps](const opengl_window* win, const opengl_input& input) {
        auto scene_ok = !apps->states.empty() && apps->selected >= 0;
        if (!scene_ok) return;

        // update trasforms
        if (scene_ok) {
          auto app = apps->states[apps->selected];
          update_transforms(app->scene, app->time);
        }

        // handle mouse and keyboard for navigation
        if (scene_ok && (input.mouse_left || input.mouse_right) &&
            !input.modifier_alt && !input.widgets_active) {
          auto  app    = apps->states[apps->selected];
          auto& camera = app->scene.cameras.at(app->drawgl_prms.camera);
          auto  dolly  = 0.0f;
          auto  pan    = zero2f;
          auto  rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_right)
            dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) / 100.0f;
          update_turntable(camera.frame, camera.focus, rotate, dolly, pan);
          set_camera(app->glscene.get(), app->drawgl_prms.camera, camera.frame,
              camera.lens, camera.aspect, camera.film, 0.001, 10000);
        }

        // animation
        if (scene_ok && apps->states[apps->selected]->animate) {
          auto app = apps->states[apps->selected];
          app->time += min(1 / 60.0f, (float)input.time_delta);
          if (app->time < app->time_range.x || app->time > app->time_range.y)
            app->time = app->time_range.x;
          update_transforms(app->scene, app->time);
        }
      });

  // run ui
  run_ui(win);

  // clear
  delete_glwindow(win);

  // done
  return 0;
}
