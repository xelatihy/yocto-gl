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
void print_obj_camera(const sceneio_camera* camera);
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
  sceneio_model ioscene = {};

  // rendering state
  opengl_scene glscene = {};

  // view image
  float  time       = 0;
  string anim_group = "";
  vec2f  time_range = zero2f;
  bool   animate    = false;

  // editing
  int selected_camera      = -1;
  int selected_shape       = -1;
  int selected_subdiv      = -1;
  int selected_material    = -1;
  int selected_environment = -1;
  int selected_texture     = -1;

  // error
  string error = "";
};

// Application state
struct app_states {
  // data
  vector<shared_ptr<app_state>> states   = {};
  int                           selected = -1;

  // loading
  deque<future<shared_ptr<app_state>>> loaders = {};

  // default options
  draw_glscene_params drawgl_prms = {};
};

// Compute animation range
vec2f compute_animation_range(
    const sceneio_model& scene, const string& anim_group = "") {
  if (scene.animations.empty()) return zero2f;
  auto range = vec2f{+flt_max, -flt_max};
  for (auto animation : scene.animations) {
    if (anim_group != "" && animation->group != anim_group) continue;
    range.x = min(range.x, animation->times.front());
    range.y = max(range.y, animation->times.back());
  }
  if (range.y < range.x) return zero2f;
  return range;
}

void load_scene_async(shared_ptr<app_states> apps, const string& filename) {
  apps->loaders.push_back(
      async(launch::async, [apps, filename]() -> shared_ptr<app_state> {
        auto app         = make_shared<app_state>();
        app->filename    = filename;
        app->imagename   = replace_extension(filename, ".png");
        app->outname     = replace_extension(filename, ".edited.yaml");
        app->name        = get_filename(app->filename);
        app->drawgl_prms = apps->drawgl_prms;
        load_scene(app->filename, app->ioscene);
        app->time_range           = compute_animation_range(app->ioscene);
        app->time                 = app->time_range.x;
        app->selected_camera      = app->ioscene.cameras.empty() ? -1 : 0;
        app->selected_shape       = app->ioscene.shapes.empty() ? -1 : 0;
        app->selected_subdiv      = app->ioscene.subdivs.empty() ? -1 : 0;
        app->selected_material    = app->ioscene.shapes.empty() ? -1 : 0;
        app->selected_texture     = app->ioscene.textures.empty() ? -1 : 0;
        app->selected_environment = app->ioscene.environments.empty() ? -1 : 0;
        return app;
      }));
}

void update_lights(opengl_scene& glscene, const sceneio_model& ioscene) {
  clear_lights(glscene);
  for (auto ioshape : ioscene.shapes) {
    if (has_max_lights(glscene)) break;
    if (ioshape->emission == zero3f) continue;
    auto bbox = invalidb3f;
    for (auto p : ioshape->positions) bbox = merge(bbox, p);
    auto pos  = (bbox.max + bbox.min) / 2;
    auto area = 0.0f;
    if (!ioshape->triangles.empty()) {
      for (auto t : ioshape->triangles)
        area += triangle_area(
            ioshape->positions[t.x], ioshape->positions[t.y], ioshape->positions[t.z]);
    } else if (!ioshape->quads.empty()) {
      for (auto q : ioshape->quads)
        area += quad_area(ioshape->positions[q.x], ioshape->positions[q.y],
            ioshape->positions[q.z], ioshape->positions[q.w]);
    } else if (!ioshape->lines.empty()) {
      for (auto l : ioshape->lines)
        area += line_length(ioshape->positions[l.x], ioshape->positions[l.y]);
    } else {
      area += ioshape->positions.size();
    }
    auto ke  = ioshape->emission * area;
    auto lid = add_light(glscene);
    set_light(glscene, lid, transform_point(ioshape->frame, pos), ke, false);
  }
}

void init_scene(opengl_scene& glscene, sceneio_model& ioscene) {
  // load program
  init_glscene(glscene);

  // camera
  for (auto iocamera : ioscene.cameras) {
    auto id = add_camera(glscene);
    set_camera_frame(glscene, id, iocamera->frame);
    set_camera_lens(glscene, id, iocamera->lens, iocamera->aspect, iocamera->film);
    set_camera_nearfar(glscene, id, 0.001, 10000);
  }

  // textures
  for (auto iotexture : ioscene.textures) {
    auto id = add_texture(glscene);
    if (!iotexture->hdr.empty()) {
      set_texture(glscene, id, iotexture->hdr);
    } else if (!iotexture->ldr.empty()) {
      set_texture(glscene, id, iotexture->ldr);
    }
  }

  for (auto iosubdiv : ioscene.subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
  }

  // shapes
  for (auto ioshape : ioscene.shapes) {
    auto id = add_shape(glscene);
    set_shape_positions(glscene, id, ioshape->positions);
    set_shape_normals(glscene, id, ioshape->normals);
    set_shape_texcoords(glscene, id, ioshape->texcoords);
    set_shape_colors(glscene, id, ioshape->colors);
    set_shape_points(glscene, id, ioshape->points);
    set_shape_lines(glscene, id, ioshape->lines);
    set_shape_triangles(glscene, id, ioshape->triangles);
    set_shape_quads(glscene, id, ioshape->quads);
    set_shape_frame(glscene, id, ioshape->frame);
    set_shape_instances(glscene, id, ioshape->instances);
    set_shape_emission(glscene, id, ioshape->emission, ioshape->emission_tex);
    set_shape_color(
        glscene, id, (1 - ioshape->transmission) * ioshape->color, ioshape->color_tex);
    set_shape_specular(glscene, id, (1 - ioshape->transmission) * ioshape->specular,
        ioshape->specular_tex);
    set_shape_metallic(glscene, id, (1 - ioshape->transmission) * ioshape->metallic,
        ioshape->metallic_tex);
    set_shape_roughness(glscene, id, ioshape->roughness, ioshape->roughness_tex);
    set_shape_opacity(glscene, id, ioshape->opacity, ioshape->opacity_tex);
    set_shape_normalmap(glscene, id, ioshape->normal_tex);
  }
}

bool draw_glwidgets_camera(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iocamera = app->ioscene.cameras[id];
  auto  edited = 0;
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

/// Visit struct elements.
bool draw_glwidgets_texture(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iotexture = app->ioscene.textures[id];
  draw_gllabel(win, "name", iotexture->name);
  draw_gllabel(win, "hdr",
      to_string(iotexture->hdr.size().x) + " x " +
          to_string(iotexture->hdr.size().y));
  draw_gllabel(win, "ldr",
      to_string(iotexture->ldr.size().x) + " x " +
          to_string(iotexture->ldr.size().y));
  return false;
}

bool draw_glwidgets_material(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioshape  = app->ioscene.shapes[id];
  auto  edited = 0;
  edited += draw_gltextinput(win, "name", app->ioscene.shapes[id]->name);
  edited += draw_glhdrcoloredit(win, "emission", ioshape->emission);
  edited += draw_glcoloredit(win, "color", ioshape->color);
  edited += draw_glslider(win, "specular", ioshape->specular, 0, 1);
  edited += draw_glslider(win, "metallic", ioshape->metallic, 0, 1);
  edited += draw_glslider(win, "roughness", ioshape->roughness, 0, 1);
  edited += draw_glslider(win, "coat", ioshape->coat, 0, 1);
  edited += draw_glslider(win, "transmission", ioshape->transmission, 0, 1);
  edited += draw_glcoloredit(win, "spectint", ioshape->spectint);
  edited += draw_glcheckbox(win, "thin", ioshape->thin);
  edited += draw_glcoloredit(win, "scattering", ioshape->scattering);
  edited += draw_glslider(win, "trdepth", ioshape->trdepth, 0, 1);
  edited += draw_glslider(win, "scanisotropy", ioshape->scanisotropy, -1, 1);
  edited += draw_glslider(win, "opacity", ioshape->opacity, 0, 1);
  edited += draw_glcombobox(
      win, "emission_tex", ioshape->emission_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "color_tex", ioshape->color_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", ioshape->metallic_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", ioshape->specular_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(win, "transmission_tex", ioshape->transmission_tex,
      app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "scattering_tex", ioshape->scattering_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", ioshape->roughness_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "spectint_tex", ioshape->spectint_tex, app->ioscene.textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", ioshape->normal_tex, app->ioscene.textures, true);
  edited += draw_glcheckbox(win, "glTF textures", ioshape->gltf_textures);
  return edited;
}

bool draw_glwidgets_shape(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioshape  = app->ioscene.shapes[id];
  auto  edited = 0;
  edited += draw_gltextinput(win, "name", ioshape->name);
  edited += draw_glslider(win, "frame[0]", ioshape->frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", ioshape->frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", ioshape->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioshape->frame.o, -10, 10);
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
  // TODO: load
  return edited;
}

bool draw_glwidgets_subdiv(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto iosubdiv = app->ioscene.subdivs[id];
  auto  edited = 0;
  edited += draw_gltextinput(win, "name", iosubdiv->name);
  draw_gllabel(win, "points", to_string(iosubdiv->points.size()));
  draw_gllabel(win, "lines", to_string(iosubdiv->lines.size()));
  draw_gllabel(win, "triangles", to_string(iosubdiv->triangles.size()));
  draw_gllabel(win, "quads", to_string(iosubdiv->quads.size()));
  draw_gllabel(win, "quads pos", to_string(iosubdiv->quadspos.size()));
  draw_gllabel(win, "quads norm", to_string(iosubdiv->quadsnorm.size()));
  draw_gllabel(win, "quads texcoord", to_string(iosubdiv->quadstexcoord.size()));
  draw_gllabel(win, "pos", to_string(iosubdiv->positions.size()));
  draw_gllabel(win, "norm", to_string(iosubdiv->normals.size()));
  draw_gllabel(win, "texcoord", to_string(iosubdiv->texcoords.size()));
  draw_gllabel(win, "color", to_string(iosubdiv->colors.size()));
  draw_gllabel(win, "radius", to_string(iosubdiv->radius.size()));
  draw_gllabel(win, "tangsp", to_string(iosubdiv->tangents.size()));
  edited += draw_glslider(win, "subdivisions", iosubdiv->subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", iosubdiv->catmullclark);
  edited += draw_glcheckbox(win, "smooth", iosubdiv->smooth);
  edited += draw_glcombobox(win, "displacement_tex", iosubdiv->displacement_tex,
      app->ioscene.textures, true);
  edited += draw_glslider(win, "displacement", iosubdiv->displacement, 0, 1);
  // TODO: load
  return edited;
}

bool draw_glwidgets_environment(
    const opengl_window& win, shared_ptr<app_state> app, int id) {
  auto ioenvironment = app->ioscene.environments[id];
  auto  edited      = 0;
  edited += draw_gltextinput(win, "name", ioenvironment->name);
  edited += draw_glslider(win, "frame[0]", ioenvironment->frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", ioenvironment->frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", ioenvironment->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_glcombobox(win, "emission texture", ioenvironment->emission_tex,
      app->ioscene.textures, true);
  return edited;
}

// draw with shading
void draw_glwidgets(const opengl_window& win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  static auto load_path = ""s, save_path = ""s, error_message = ""s;
  auto        app = (!apps->states.empty() && apps->selected >= 0)
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
      save_scene(app->outname, app->ioscene);
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
  if (apps->states.empty()) return;
  draw_glcombobox(
      win, "scene", apps->selected, (int)apps->states.size(),
      [apps](int idx) { return apps->states[idx]->name.c_str(); }, false);
  if (app && begin_glheader(win, "view")) {
    auto  app    = apps->states[apps->selected];
    auto& params = app->drawgl_prms;
    draw_glcombobox(win, "camera", params.camera, app->ioscene.cameras);
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
  if (app && begin_glheader(win, "inspect")) {
    draw_gllabel(win, "scene", get_filename(app->filename));
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "imagename", app->imagename);
    continue_glline(win);
    if (draw_glbutton(win, "print cams")) {
      for (auto iocamera : app->ioscene.cameras) {
        print_obj_camera(iocamera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : scene_stats(app->ioscene)) print_info(stat);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.cameras.empty() && begin_glheader(win, "cameras")) {
    draw_glcombobox(
        win, "camera##2", app->selected_camera, app->ioscene.cameras);
    if (draw_glwidgets_camera(win, app, app->selected_camera)) {
      auto iocamera = app->ioscene.cameras[app->selected_camera];
      set_camera_frame(app->glscene, app->selected_camera, iocamera->frame);
      set_camera_lens(app->glscene, app->selected_camera, iocamera->lens,
          iocamera->aspect, iocamera->film);
      set_camera_nearfar(app->glscene, app->selected_camera, 0.001, 10000);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.environments.empty() &&
      begin_glheader(win, "environments")) {
    draw_glcombobox(win, "environments##2", app->selected_environment,
        app->ioscene.environments);
    if (draw_glwidgets_environment(win, app, app->selected_environment)) {
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.shapes.empty() && begin_glheader(win, "shapes")) {
    draw_glcombobox(win, "shape##2", app->selected_shape, app->ioscene.shapes);
    if (!draw_glwidgets_shape(win, app, app->selected_shape)) {
      auto ioshape = app->ioscene.shapes[app->selected_shape];
      auto  idx   = app->selected_shape;
      set_shape_positions(app->glscene, idx, ioshape->positions);
      set_shape_normals(app->glscene, idx, ioshape->normals);
      set_shape_texcoords(app->glscene, idx, ioshape->texcoords);
      set_shape_colors(app->glscene, idx, ioshape->colors);
      set_shape_points(app->glscene, idx, ioshape->points);
      set_shape_lines(app->glscene, idx, ioshape->lines);
      set_shape_triangles(app->glscene, idx, ioshape->triangles);
      set_shape_quads(app->glscene, idx, ioshape->quads);
      set_shape_frame(app->glscene, idx, ioshape->frame);
      set_shape_instances(app->glscene, idx, ioshape->instances);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.shapes.empty() && begin_glheader(win, "materials")) {
    draw_glcombobox(
        win, "material##2", app->selected_material, app->ioscene.shapes);
    if (draw_glwidgets_material(win, app, app->selected_material)) {
      auto ioshape = app->ioscene.shapes[app->selected_material];
      set_shape_emission(app->glscene, app->selected_material, ioshape->emission,
          ioshape->emission_tex);
      set_shape_color(app->glscene, app->selected_material,
          (1 - ioshape->transmission) * ioshape->color, ioshape->color_tex);
      set_shape_specular(app->glscene, app->selected_material,
          (1 - ioshape->transmission) * ioshape->specular, ioshape->specular_tex);
      set_shape_metallic(app->glscene, app->selected_material,
          (1 - ioshape->transmission) * ioshape->metallic, ioshape->metallic_tex);
      set_shape_roughness(app->glscene, app->selected_material, ioshape->roughness,
          ioshape->roughness_tex);
      set_shape_opacity(app->glscene, app->selected_material, ioshape->opacity,
          ioshape->opacity_tex);
      set_shape_normalmap(
          app->glscene, app->selected_material, ioshape->normal_tex);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.textures.empty() &&
      begin_glheader(win, "textures")) {
    draw_glcombobox(
        win, "texture##2", app->selected_texture, app->ioscene.textures);
    if (draw_glwidgets_texture(win, app, app->selected_texture)) {
      auto iotexture = app->ioscene.textures[app->selected_texture];
      if (!iotexture->hdr.empty()) {
        set_texture(app->glscene, app->selected_texture, iotexture->hdr);
      } else if (!iotexture->hdr.empty()) {
        set_texture(app->glscene, app->selected_texture, iotexture->ldr);
      }
    }
    end_glheader(win);
  }
  if (app && !app->ioscene.subdivs.empty() && begin_glheader(win, "subdivs")) {
    draw_glcombobox(
        win, "subdiv##2", app->selected_subdiv, app->ioscene.subdivs);
    if (!draw_glwidgets_subdiv(win, app, app->selected_subdiv)) {
      auto iosubdiv = app->ioscene.subdivs[app->selected_subdiv];
      tesselate_subdiv(app->ioscene, iosubdiv);
      auto ioshape = app->ioscene.shapes[iosubdiv->shape];
      auto  idx   = app->selected_subdiv;
      set_shape_positions(app->glscene, idx, ioshape->positions);
      set_shape_normals(app->glscene, idx, ioshape->normals);
      set_shape_texcoords(app->glscene, idx, ioshape->texcoords);
      set_shape_colors(app->glscene, idx, ioshape->colors);
      set_shape_points(app->glscene, idx, ioshape->points);
      set_shape_lines(app->glscene, idx, ioshape->lines);
      set_shape_triangles(app->glscene, idx, ioshape->triangles);
      set_shape_quads(app->glscene, idx, ioshape->quads);
      set_shape_frame(app->glscene, idx, ioshape->frame);
      set_shape_instances(app->glscene, idx, ioshape->instances);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(const opengl_window& win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->states.empty() && apps->selected >= 0) {
    auto app = apps->states[apps->selected];
    draw_glscene(app->glscene, input.framebuffer_viewport, app->drawgl_prms);
  }
}

// update
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
      init_scene(app->glscene, app->ioscene);
      update_lights(app->glscene, app->ioscene);
      if (apps->selected < 0) apps->selected = (int)apps->states.size() - 1;
    } catch (std::exception& e) {
      apps->loaders.pop_front();
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
  }
}

void run_app(int argc, const char* argv[]) {
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
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  auto win = opengl_window{};
  init_glwindow(win, {1280 + 320, 720}, "yscnview", true);

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
  set_uiupdate_glcallback(
      win, [apps](const opengl_window& win, const opengl_input& input) {
        auto scene_ok = !apps->states.empty() && apps->selected >= 0;
        if (!scene_ok) return;

        // update trasforms
        if (scene_ok) {
          auto app = apps->states[apps->selected];
          update_transforms(app->ioscene, app->time);
        }

        // handle mouse and keyboard for navigation
        if (scene_ok && (input.mouse_left || input.mouse_right) &&
            !input.modifier_alt && !input.widgets_active) {
          auto  app    = apps->states[apps->selected];
          auto  iocamera = app->ioscene.cameras.at(app->drawgl_prms.camera);
          auto  dolly  = 0.0f;
          auto  pan    = zero2f;
          auto  rotate = zero2f;
          if (input.mouse_left && !input.modifier_shift)
            rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
          if (input.mouse_right)
            dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
          if (input.mouse_left && input.modifier_shift)
            pan = (input.mouse_pos - input.mouse_last) / 100.0f;
          update_turntable(iocamera->frame, iocamera->focus, rotate, dolly, pan);
          set_camera_frame(app->glscene, app->drawgl_prms.camera, iocamera->frame);
        }

        // animation
        if (scene_ok && apps->states[apps->selected]->animate) {
          auto app = apps->states[apps->selected];
          app->time += min(1 / 60.0f, (float)input.time_delta);
          if (app->time < app->time_range.x || app->time > app->time_range.y)
            app->time = app->time_range.x;
          update_transforms(app->ioscene, app->time);
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
