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

#include "../yocto/yocto_image.h"
#include "../yocto/yocto_sceneio.h"
#include "../yocto/yocto_shape.h"
#include "yocto_opengl.h"
using namespace yocto;

#include <deque>
#include <future>
#include <memory>
using namespace std;

#include "ext/CLI11.hpp"
#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(shared_ptr<sceneio_camera> camera);
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
  shared_ptr<sceneio_model> ioscene = make_shared<sceneio_model>();

  // rendering state
  shared_ptr<opengl_scene> glscene = nullptr;

  // view image
  float  time       = 0;
  string anim_group = "";
  vec2f  time_range = zero2f;
  bool   animate    = false;

  // editing
  shared_ptr<sceneio_camera>      selected_camera      = nullptr;
  shared_ptr<sceneio_object>      selected_object      = nullptr;
  shared_ptr<sceneio_instance>    selected_instance    = nullptr;
  shared_ptr<sceneio_shape>       selected_shape       = nullptr;
  shared_ptr<sceneio_subdiv>      selected_subdiv      = nullptr;
  shared_ptr<sceneio_material>    selected_material    = nullptr;
  shared_ptr<sceneio_environment> selected_environment = nullptr;
  shared_ptr<sceneio_texture>     selected_texture     = nullptr;

  // editing maps
  unordered_map<shared_ptr<sceneio_camera>, shared_ptr<opengl_camera>>
      camera_map = {};
  unordered_map<shared_ptr<sceneio_texture>, shared_ptr<opengl_texture>>
      texture_map = {};
  unordered_map<shared_ptr<sceneio_material>, shared_ptr<opengl_material>>
                                                                     material_map = {};
  unordered_map<shared_ptr<sceneio_shape>, shared_ptr<opengl_shape>> shape_map =
      {};
  unordered_map<shared_ptr<sceneio_instance>, shared_ptr<opengl_instance>>
      instance_map = {};
  unordered_map<shared_ptr<sceneio_object>, shared_ptr<opengl_object>>
      object_map = {};
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
    shared_ptr<sceneio_model> ioscene, const string& anim_group = "") {
  if (ioscene->animations.empty()) return zero2f;
  auto range = vec2f{+flt_max, -flt_max};
  for (auto animation : ioscene->animations) {
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
        app->imagename   = fs::path(filename).replace_extension(".png");
        app->outname     = fs::path(filename).replace_extension(".edited.yaml");
        app->name        = fs::path(app->filename).filename();
        app->drawgl_prms = apps->drawgl_prms;
        app->ioscene     = load_scene(app->filename);
        app->time_range  = compute_animation_range(app->ioscene);
        app->time        = app->time_range.x;
        return app;
      }));
}

void update_lights(
    shared_ptr<opengl_scene> glscene, shared_ptr<sceneio_model> ioscene) {
  clear_lights(glscene);
  for (auto ioobject : ioscene->objects) {
    if (has_max_lights(glscene)) break;
    if (ioobject->material->emission == zero3f) continue;
    auto ioshape = ioobject->shape;
    auto bbox    = invalidb3f;
    for (auto p : ioshape->positions) bbox = merge(bbox, p);
    auto pos  = (bbox.max + bbox.min) / 2;
    auto area = 0.0f;
    if (!ioshape->triangles.empty()) {
      for (auto t : ioshape->triangles)
        area += triangle_area(ioshape->positions[t.x], ioshape->positions[t.y],
            ioshape->positions[t.z]);
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
    auto ke = ioobject->material->emission * area;
    set_light(
        add_light(glscene), transform_point(ioobject->frame, pos), ke, false);
  }
}

void init_scene(shared_ptr<app_state> app) {
  // create scene
  app->glscene = make_glscene();

  auto glscene = app->glscene;
  auto ioscene = app->ioscene;

  // camera
  app->camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    auto camera = add_camera(glscene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  app->texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    auto gltexture = add_texture(glscene);
    if (!iotexture->hdr.empty()) {
      set_texture(gltexture, iotexture->hdr);
    } else if (!iotexture->ldr.empty()) {
      set_texture(gltexture, iotexture->ldr);
    }
    app->texture_map[iotexture] = gltexture;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    tesselate_subdiv(ioscene, iosubdiv);
  }

  // material
  app->material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    auto glmaterial = add_material(glscene);
    set_emission(glmaterial, iomaterial->emission,
        app->texture_map.at(iomaterial->emission_tex));
    set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
        app->texture_map.at(iomaterial->color_tex));
    set_specular(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->specular,
        app->texture_map.at(iomaterial->specular_tex));
    set_metallic(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->metallic,
        app->texture_map.at(iomaterial->metallic_tex));
    set_roughness(glmaterial, iomaterial->roughness,
        app->texture_map.at(iomaterial->roughness_tex));
    set_opacity(glmaterial, iomaterial->opacity,
        app->texture_map.at(iomaterial->opacity_tex));
    set_normalmap(glmaterial, app->texture_map.at(iomaterial->normal_tex));
    app->material_map[iomaterial] = glmaterial;
  }

  // shapes
  app->shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    auto glshape = add_shape(glscene);
    set_positions(glshape, ioshape->positions);
    set_normals(glshape, ioshape->normals);
    set_texcoords(glshape, ioshape->texcoords);
    set_colors(glshape, ioshape->colors);
    set_points(glshape, ioshape->points);
    set_lines(glshape, ioshape->lines);
    set_triangles(glshape, ioshape->triangles);
    set_quads(glshape, ioshape->quads);
    app->shape_map[ioshape] = glshape;
  }

  // instances
  app->instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    auto glinstance = add_instance(glscene);
    set_frames(glinstance, ioinstance->frames);
    app->instance_map[ioinstance] = glinstance;
  }

  // shapes
  app->object_map[nullptr] = nullptr;
  for (auto ioobject : ioscene->objects) {
    auto globject = add_object(glscene);
    set_frame(globject, ioobject->frame);
    set_shape(globject, app->shape_map.at(ioobject->shape));
    set_material(globject, app->material_map.at(ioobject->material));
    set_instance(globject, app->instance_map.at(ioobject->instance));
    app->object_map[ioobject] = globject;
  }
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model> ioscene, shared_ptr<sceneio_camera> iocamera) {
  if (!iocamera) return false;
  auto edited = 0;
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
bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model> ioscene, shared_ptr<sceneio_texture> iotexture) {
  if (!iotexture) return false;
  draw_gllabel(win, "name", iotexture->name);
  draw_gllabel(win, "hdr",
      to_string(iotexture->hdr.size().x) + " x " +
          to_string(iotexture->hdr.size().y));
  draw_gllabel(win, "ldr",
      to_string(iotexture->ldr.size().x) + " x " +
          to_string(iotexture->ldr.size().y));
  return false;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model>                 ioscene,
    shared_ptr<sceneio_material>              iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
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
  edited += draw_glcheckbox(win, "glTF textures", iomaterial->gltf_textures);
  return edited;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model> ioscene, shared_ptr<sceneio_shape> ioshape) {
  if (!ioshape) return false;
  auto edited = 0;
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
  // TODO: load
  return edited;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model>                 ioscene,
    shared_ptr<sceneio_instance>              ioinstance) {
  if (!ioinstance) return false;
  auto edited = 0;
  edited += draw_gltextinput(win, "name", ioinstance->name);
  draw_gllabel(win, "frames", to_string(ioinstance->frames.size()));
  // TODO: load
  return edited;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model> ioscene, shared_ptr<sceneio_object> ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  edited += draw_gltextinput(win, "name", ioobject->name);
  edited += draw_glslider(win, "frame[0]", ioobject->frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", ioobject->frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", ioobject->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_glcombobox(win, "shape", ioobject->shape, ioscene->shapes);
  edited += draw_glcombobox(
      win, "material", ioobject->material, ioscene->materials);
  edited += draw_glcombobox(
      win, "instance", ioobject->instance, ioscene->instances, true);
  // TODO: load
  return edited;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model> ioscene, shared_ptr<sceneio_subdiv> iosubdiv) {
  if (!iosubdiv) return false;
  auto edited = 0;
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
  edited += draw_glslider(win, "subdivisions", iosubdiv->subdivisions, 0, 10);
  edited += draw_glcheckbox(win, "catmullclark", iosubdiv->catmullclark);
  edited += draw_glcheckbox(win, "smooth", iosubdiv->smooth);
  edited += draw_glcombobox(win, "displacement_tex", iosubdiv->displacement_tex,
      ioscene->textures, true);
  edited += draw_glslider(win, "displacement", iosubdiv->displacement, 0, 1);
  // TODO: load
  return edited;
}

bool draw_glwidgets(shared_ptr<opengl_window> win,
    shared_ptr<sceneio_model>                 ioscene,
    shared_ptr<sceneio_environment>           ioenvironment) {
  if (!ioenvironment) return false;
  auto edited = 0;
  edited += draw_gltextinput(win, "name", ioenvironment->name);
  edited += draw_glslider(win, "frame[0]", ioenvironment->frame.x, -1, 1);
  edited += draw_glslider(win, "frame[1]", ioenvironment->frame.y, -1, 1);
  edited += draw_glslider(win, "frame[2]", ioenvironment->frame.z, -1, 1);
  edited += draw_glslider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_glhdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_glcombobox(win, "emission texture",
      ioenvironment->emission_tex, ioscene->textures, true);
  return edited;
}

// draw with shading
void draw_glwidgets(shared_ptr<opengl_window> win, shared_ptr<app_states> apps,
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
          fs::path(save_path).parent_path(), fs::path(save_path).filename(),
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
    draw_glcombobox(win, "camera", params.camera, app->ioscene->cameras);
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
    draw_gllabel(win, "scene", fs::path(app->filename).filename());
    draw_gllabel(win, "filename", app->filename);
    draw_gllabel(win, "outname", app->outname);
    draw_gllabel(win, "imagename", app->imagename);
    continue_glline(win);
    if (draw_glbutton(win, "print cams")) {
      for (auto iocamera : app->ioscene->cameras) {
        print_obj_camera(iocamera);
      }
    }
    continue_glline(win);
    if (draw_glbutton(win, "print stats")) {
      for (auto stat : scene_stats(app->ioscene)) std::cout << stat << "\n";
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->cameras.empty() && begin_glheader(win, "cameras")) {
    draw_glcombobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras);
    if (draw_glwidgets(win, app->ioscene, app->selected_camera)) {
      auto iocamera = app->selected_camera;
      auto glcamera = app->camera_map.at(iocamera);
      set_frame(glcamera, iocamera->frame);
      set_lens(glcamera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_nearfar(glcamera, 0.001, 10000);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->environments.empty() &&
      begin_glheader(win, "environments")) {
    draw_glcombobox(win, "environments##2", app->selected_environment,
        app->ioscene->environments);
    if (draw_glwidgets(win, app->ioscene, app->selected_environment)) {
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->objects.empty() && begin_glheader(win, "objects")) {
    draw_glcombobox(
        win, "object##2", app->selected_object, app->ioscene->objects);
    if (!draw_glwidgets(win, app->ioscene, app->selected_object)) {
      auto ioobject = app->selected_object;
      auto globject = app->object_map.at(ioobject);
      set_frame(globject, ioobject->frame);
      set_shape(globject, app->shape_map.at(ioobject->shape));
      set_material(globject, app->material_map.at(ioobject->material));
      set_instance(globject, app->instance_map.at(ioobject->instance));
      // TODO: add the rest
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->shapes.empty() && begin_glheader(win, "shapes")) {
    draw_glcombobox(win, "shape##2", app->selected_shape, app->ioscene->shapes);
    if (!draw_glwidgets(win, app->ioscene, app->selected_shape)) {
      auto ioshape = app->selected_shape;
      auto glshape = app->shape_map.at(ioshape);
      set_positions(glshape, ioshape->positions);
      set_normals(glshape, ioshape->normals);
      set_texcoords(glshape, ioshape->texcoords);
      set_colors(glshape, ioshape->colors);
      set_points(glshape, ioshape->points);
      set_lines(glshape, ioshape->lines);
      set_triangles(glshape, ioshape->triangles);
      set_quads(glshape, ioshape->quads);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->materials.empty() &&
      begin_glheader(win, "materials")) {
    draw_glcombobox(
        win, "material##2", app->selected_material, app->ioscene->materials);
    if (draw_glwidgets(win, app->ioscene, app->selected_material)) {
      auto iomaterial = app->selected_material;
      auto glmaterial = app->material_map.at(iomaterial);
      set_emission(glmaterial, iomaterial->emission,
          app->texture_map.at(iomaterial->emission_tex));
      set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
          app->texture_map.at(iomaterial->color_tex));
      set_specular(glmaterial,
          (1 - iomaterial->transmission) * iomaterial->specular,
          app->texture_map.at(iomaterial->specular_tex));
      set_metallic(glmaterial,
          (1 - iomaterial->transmission) * iomaterial->metallic,
          app->texture_map.at(iomaterial->metallic_tex));
      set_roughness(glmaterial, iomaterial->roughness,
          app->texture_map.at(iomaterial->roughness_tex));
      set_opacity(glmaterial, iomaterial->opacity,
          app->texture_map.at(iomaterial->opacity_tex));
      set_normalmap(glmaterial, app->texture_map.at(iomaterial->normal_tex));
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->instances.empty() &&
      begin_glheader(win, "instances")) {
    draw_glcombobox(
        win, "instance##2", app->selected_instance, app->ioscene->instances);
    if (!draw_glwidgets(win, app->ioscene, app->selected_instance)) {
      auto ioinstance = app->selected_instance;
      auto glinstance = app->instance_map.at(ioinstance);
      set_frames(glinstance, ioinstance->frames);
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->textures.empty() &&
      begin_glheader(win, "textures")) {
    draw_glcombobox(
        win, "texture##2", app->selected_texture, app->ioscene->textures);
    if (draw_glwidgets(win, app->ioscene, app->selected_texture)) {
      auto iotexture = app->selected_texture;
      auto gltexture = app->texture_map.at(iotexture);
      if (!iotexture->hdr.empty()) {
        set_texture(gltexture, iotexture->hdr);
      } else if (!iotexture->hdr.empty()) {
        set_texture(gltexture, iotexture->ldr);
      }
    }
    end_glheader(win);
  }
  if (app && !app->ioscene->subdivs.empty() && begin_glheader(win, "subdivs")) {
    draw_glcombobox(
        win, "subdiv##2", app->selected_subdiv, app->ioscene->subdivs);
    if (!draw_glwidgets(win, app->ioscene, app->selected_subdiv)) {
      // auto iosubdiv = app->ioscene->subdivs[app->selected_subdiv];
      // tesselate_subdiv(app->ioscene, iosubdiv);
      // TODO: FIX SUBDIVS
      // auto glshape = app->glscene->shapes[app->selected_shape];
      // auto ioshape = iosubdiv->shape;
      // set_positions(glshape, ioshape->positions);
      // set_normals(glshape, ioshape->normals);
      // set_texcoords(glshape, ioshape->texcoords);
      // set_colors(glshape, ioshape->colors);
      // set_points(glshape, ioshape->points);
      // set_lines(glshape, ioshape->lines);
      // set_triangles(glshape, ioshape->triangles);
      // set_quads(glshape, ioshape->quads);
    }
    end_glheader(win);
  }
  if (begin_glheader(win, "log")) {
    draw_gllog(win);
    end_glheader(win);
  }
}

// draw with shading
void draw(shared_ptr<opengl_window> win, shared_ptr<app_states> apps,
    const opengl_input& input) {
  if (!apps->states.empty() && apps->selected >= 0) {
    auto app = apps->states[apps->selected];
    draw_glscene(app->glscene, input.framebuffer_viewport, app->drawgl_prms);
  }
}

// update
void update(shared_ptr<opengl_window> win, shared_ptr<app_states> apps) {
  auto is_ready = [](const future<shared_ptr<app_state>>& result) -> bool {
    return result.valid() &&
           result.wait_for(chrono::microseconds(0)) == future_status::ready;
  };

  while (!apps->loaders.empty() && is_ready(apps->loaders.front())) {
    try {
      auto app = apps->loaders.front().get();
      apps->loaders.pop_front();
      apps->states.push_back(app);
      init_scene(app);
      update_lights(app->glscene, app->ioscene);
      if (apps->selected < 0) apps->selected = (int)apps->states.size() - 1;
    } catch (std::exception& e) {
      apps->loaders.pop_front();
      push_glmessage(win, e.what());
      log_glinfo(win, e.what());
    }
  }
}

int run_app(int argc, const char* argv[]) {
  // initialize app
  auto apps       = make_shared<app_states>();
  auto filenames  = vector<string>{};
  auto noparallel = false;

  // parse command line
  auto cli = CLI::App{"views scenes inteactively"};
  cli.add_option("--camera", apps->drawgl_prms.camera, "Camera index.");
  cli.add_option(
      "--resolution,-r", apps->drawgl_prms.resolution, "Image resolution.");
  cli.add_flag("--eyelight!,--no-eyelight,-c", apps->drawgl_prms.eyelight,
      "Eyelight rendering.");
  cli.add_flag("--noparallel", noparallel, "Disable parallel execution.");
  cli.add_option("scenes", filenames, "Scene filenames")->required();
  try {
    cli.parse(argc, argv);
  } catch (CLI::ParseError& e) {
    return cli.exit(e);
  }

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  auto win = make_glwindow({1280 + 320, 720}, "yscnview", true);

  // callbacks
  set_draw_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        draw(win, apps, input);
      });
  set_widgets_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        draw_glwidgets(win, apps, input);
      });
  set_drop_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const vector<string>& paths,
               const opengl_input& input) {
        for (auto& path : paths) load_scene_async(apps, path);
      });
  set_update_glcallback(
      win, [apps](shared_ptr<opengl_window> win, const opengl_input& input) {
        update(win, apps);
      });
  set_uiupdate_glcallback(win, [apps](shared_ptr<opengl_window> win,
                                   const opengl_input&          input) {
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
      auto app      = apps->states[apps->selected];
      auto iocamera = app->ioscene->cameras.at(app->drawgl_prms.camera);
      auto dolly    = 0.0f;
      auto pan      = zero2f;
      auto rotate   = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      update_turntable(iocamera->frame, iocamera->focus, rotate, dolly, pan);
      set_frame(
          app->glscene->cameras[app->drawgl_prms.camera], iocamera->frame);
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

  // done
  return 0;
}

int main(int argc, const char* argv[]) {
  try {
    return run_app(argc, argv);
  } catch (std::exception& e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
}
