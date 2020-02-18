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

#include <atomic>
#include <deque>
#include <future>
using std::atomic;
using std::deque;
using std::future;
using std::to_string;

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(sceneio_camera* camera);
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
  sceneio_model* ioscene = new sceneio_model{};

  // rendering state
  opengl_scene* glscene = new opengl_scene{};

  // editing
  sceneio_camera*      selected_camera      = nullptr;
  sceneio_object*      selected_object      = nullptr;
  sceneio_instance*    selected_instance    = nullptr;
  sceneio_shape*       selected_shape       = nullptr;
  sceneio_subdiv*      selected_subdiv      = nullptr;
  sceneio_material*    selected_material    = nullptr;
  sceneio_environment* selected_environment = nullptr;
  sceneio_texture*     selected_texture     = nullptr;

  // loading status
  atomic<bool>       ok           = false;
  future<void>       loader       = {};
  string             status       = "";
  string             error        = "";
  std::atomic<float> progress     = 0.5;
  string             loader_error = "";

  ~app_state() {
    if (ioscene) delete ioscene;
    if (glscene) delete glscene;
  }
};

// Application state
struct app_states {
  // data
  vector<app_state*> states   = {};
  app_state*         selected = nullptr;
  deque<app_state*>  loading  = {};

  // default options
  draw_glscene_params drawgl_prms = {};

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

void load_scene_async(app_states* apps, const string& filename) {
  auto app         = apps->states.emplace_back(new app_state{});
  app->filename    = filename;
  app->imagename   = fs::path(filename).replace_extension(".png");
  app->outname     = fs::path(filename).replace_extension(".edited.yaml");
  app->name        = fs::path(app->filename).filename();
  app->drawgl_prms = apps->drawgl_prms;
  app->status      = "load";
  app->loader      = std::async(std::launch::async, [app]() {
    auto progress_cb = [app](const string& message, int current, int total) {
      app->progress = (float)current / (float)total;
    };
    if (!load_scene(
            app->filename, app->ioscene, app->loader_error, progress_cb))
      return;
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = app;
}

void update_lights(opengl_scene* glscene, sceneio_model* ioscene) {
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

void init_glscene(opengl_scene* glscene, sceneio_model* ioscene,
    sceneio_progress progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->materials.size() +
             (int)ioscene->textures.size() + (int)ioscene->shapes.size() +
             (int)ioscene->subdivs.size() + (int)ioscene->instances.size() +
             (int)ioscene->objects.size()};

  // create scene
  init_glscene(glscene);

  // camera
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto camera = add_camera(glscene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  auto texture_map     = unordered_map<sceneio_texture*, opengl_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto gltexture = add_texture(glscene);
    if (!iotexture->colorf.empty()) {
      set_texture(gltexture, iotexture->colorf);
    } else if (!iotexture->colorb.empty()) {
      set_texture(gltexture, iotexture->colorb);
    } else if (!iotexture->scalarf.empty()) {
      set_texture(gltexture, iotexture->scalarf);
    } else if (!iotexture->scalarb.empty()) {
      set_texture(gltexture, iotexture->scalarb);
    }
    texture_map[iotexture] = gltexture;
  }

  // material
  auto material_map     = unordered_map<sceneio_material*, opengl_material*>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto glmaterial = add_material(glscene);
    set_emission(glmaterial, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
        texture_map.at(iomaterial->color_tex));
    set_specular(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_metallic(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_roughness(glmaterial, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(glmaterial, iomaterial->opacity,
        texture_map.at(iomaterial->opacity_tex));
    set_normalmap(glmaterial, texture_map.at(iomaterial->normal_tex));
    material_map[iomaterial] = glmaterial;
  }

  for (auto iosubdiv : ioscene->subdivs) {
    if (progress_cb) progress_cb("convert subdiv", progress.x++, progress.y);
    tesselate_subdiv(ioscene, iosubdiv);
  }

  // shapes
  auto shape_map     = unordered_map<sceneio_shape*, opengl_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    auto glshape = add_shape(glscene);
    set_positions(glshape, ioshape->positions);
    set_normals(glshape, ioshape->normals);
    set_texcoords(glshape, ioshape->texcoords);
    set_colors(glshape, ioshape->colors);
    set_points(glshape, ioshape->points);
    set_lines(glshape, ioshape->lines);
    set_triangles(glshape, ioshape->triangles);
    set_quads(glshape, ioshape->quads);
    shape_map[ioshape] = glshape;
  }

  // instances
  auto instance_map     = unordered_map<sceneio_instance*, opengl_instance*>{};
  instance_map[nullptr] = nullptr;
  for (auto ioinstance : ioscene->instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto glinstance = add_instance(glscene);
    set_frames(glinstance, ioinstance->frames);
    instance_map[ioinstance] = glinstance;
  }

  // shapes
  for (auto ioobject : ioscene->objects) {
    if (progress_cb) progress_cb("convert object", progress.x++, progress.y);
    auto globject = add_object(glscene);
    set_frame(globject, ioobject->frame);
    set_shape(globject, shape_map.at(ioobject->shape));
    set_material(globject, material_map.at(ioobject->material));
    set_instance(globject, instance_map.at(ioobject->instance));
  }

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);
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

/// Visit struct elements.
bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_texture* iotexture) {
  if (!iotexture) return false;
  draw_gllabel(win, "name", iotexture->name);
  draw_gllabel(win, "colorf",
      std::to_string(iotexture->colorf.size().x) + " x " +
          std::to_string(iotexture->colorf.size().y));
  draw_gllabel(win, "colorb",
      std::to_string(iotexture->colorb.size().x) + " x " +
          std::to_string(iotexture->colorb.size().y));
  draw_gllabel(win, "scalarf",
      std::to_string(iotexture->scalarf.size().x) + " x " +
          std::to_string(iotexture->scalarf.size().y));
  draw_gllabel(win, "scalarb",
      std::to_string(iotexture->scalarb.size().x) + " x " +
          std::to_string(iotexture->scalarb.size().y));
  return false;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_material* iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
  draw_gllabel(win, "name", iomaterial->name);
  edited += draw_glhdrcoloredit(win, "emission", iomaterial->emission);
  edited += draw_glcoloredit(win, "color", iomaterial->color);
  edited += draw_glslider(win, "opacity", iomaterial->opacity, 0, 1);
  edited += draw_glslider(win, "metallic", iomaterial->metallic, 0, 1);
  edited += draw_glslider(win, "roughness", iomaterial->roughness, 0, 1);
  edited += draw_glslider(win, "specular", iomaterial->specular, 0, 1);
  edited += draw_glslider(win, "coat", iomaterial->coat, 0, 1);
  edited += draw_glslider(win, "transmission", iomaterial->transmission, 0, 1);
  edited += draw_glcoloredit(win, "spectint", iomaterial->spectint);
  edited += draw_glcheckbox(win, "thin", iomaterial->thin);
  edited += draw_glcoloredit(win, "scattering", iomaterial->scattering);
  edited += draw_glslider(win, "trdepth", iomaterial->trdepth, 0, 1);
  edited += draw_glslider(win, "scanisotropy", iomaterial->scanisotropy, -1, 1);
  edited += draw_glslider(win, "displacement", iomaterial->displacement, 0, 1);
  edited += draw_glcombobox(
      win, "emission_tex", iomaterial->emission_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "color_tex", iomaterial->color_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "opacity_tex", iomaterial->opacity_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "metallic_tex", iomaterial->metallic_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "roughness_tex", iomaterial->roughness_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "specular_tex", iomaterial->specular_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "transmission_tex",
      iomaterial->transmission_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "scattering_tex", iomaterial->scattering_tex,
      ioscene->textures, true);
  edited += draw_glcombobox(
      win, "spectint_tex", iomaterial->spectint_tex, ioscene->textures, true);
  edited += draw_glcombobox(
      win, "normal_tex", iomaterial->normal_tex, ioscene->textures, true);
  edited += draw_glcombobox(win, "displacement_tex",
      iomaterial->displacement_tex, ioscene->textures, true);
  edited += draw_glslider(win, "subdivisions", iomaterial->subdivisions, 0, 5);
  edited += draw_glcheckbox(win, "smooth", iomaterial->smooth);
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
  // TODO: load
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_instance* ioinstance) {
  if (!ioinstance) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioinstance->name);
  draw_gllabel(win, "frames", to_string(ioinstance->frames.size()));
  // TODO: load
  return edited;
}

bool draw_glwidgets(
    opengl_window* win, sceneio_model* ioscene, sceneio_object* ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  draw_gllabel(win, "name", ioobject->name);
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
  // TODO: load
  return edited;
}

bool draw_glwidgets(opengl_window* win, sceneio_model* ioscene,
    sceneio_environment* ioenvironment) {
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

template <typename T, typename T1>
T1* get_element(
    T* ioelement, const vector<T*>& ioelements, const vector<T1*>& elements) {
  if (!ioelement) return nullptr;
  for (auto pos = 0; pos < ioelements.size(); pos++) {
    if (ioelements[pos] == ioelement) return elements[pos];
  }
  throw std::runtime_error("element not found");
}

// draw with shading
void draw_glwidgets(
    opengl_window* win, app_states* apps, const opengl_input& input) {
  static auto load_path = ""s, save_path = ""s, error_message = ""s;
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
    save_scene(app->outname, app->ioscene, app->error);
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
  if (apps->states.empty()) return;
  draw_glcombobox(win, "scene", apps->selected, apps->states, false);
  if (!apps->selected) return;
  draw_glprogressbar(
      win, apps->selected->status.c_str(), apps->selected->progress);
  if (apps->selected->error != "") {
    draw_gllabel(win, "error", apps->selected->error);
    return;
  }
  if (!apps->selected->ok) return;
  auto app = apps->selected;
  if (begin_glheader(win, "view")) {
    auto& params = app->drawgl_prms;
    draw_glcombobox(win, "camera", params.camera, app->ioscene->cameras);
    draw_glslider(win, "resolution", params.resolution, 0, 4096);
    draw_glcheckbox(win, "eyelight", params.eyelight);
    continue_glline(win);
    draw_glcheckbox(win, "wireframe", params.wireframe);
    continue_glline(win);
    draw_glcheckbox(win, "edges", params.edges);
    draw_glslider(win, "exposure", params.exposure, -10, 10);
    draw_glslider(win, "gamma", params.gamma, 0.1f, 4);
    draw_glcheckbox(win, "double sided", params.double_sided);
    draw_glslider(win, "near", params.near, 0.01f, 1.0f);
    draw_glslider(win, "far", params.far, 1000.0f, 10000.0f);
    end_glheader(win);
  }
  if (begin_glheader(win, "inspect")) {
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
      for (auto stat : scene_stats(app->ioscene)) print_info(stat);
    }
    end_glheader(win);
  }
  auto get_texture = [app](sceneio_texture* iotexture) {
    return get_element(
        iotexture, app->ioscene->textures, app->glscene->textures);
  };
  if (!app->ioscene->cameras.empty() && begin_glheader(win, "cameras")) {
    draw_glcombobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras);
    if (draw_glwidgets(win, app->ioscene, app->selected_camera)) {
      auto iocamera = app->selected_camera;
      auto glcamera = get_element(
          iocamera, app->ioscene->cameras, app->glscene->cameras);
      set_frame(glcamera, iocamera->frame);
      set_lens(glcamera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_nearfar(glcamera, 0.001, 10000);
    }
    end_glheader(win);
  }
  if (!app->ioscene->environments.empty() &&
      begin_glheader(win, "environments")) {
    draw_glcombobox(win, "environments##2", app->selected_environment,
        app->ioscene->environments);
    if (draw_glwidgets(win, app->ioscene, app->selected_environment)) {
    }
    end_glheader(win);
  }
  if (!app->ioscene->objects.empty() && begin_glheader(win, "objects")) {
    draw_glcombobox(
        win, "object##2", app->selected_object, app->ioscene->objects);
    if (!draw_glwidgets(win, app->ioscene, app->selected_object)) {
      auto ioobject = app->selected_object;
      auto globject = get_element(
          ioobject, app->ioscene->objects, app->glscene->objects);
      set_frame(globject, ioobject->frame);
      set_shape(globject, get_element(ioobject->shape, app->ioscene->shapes,
                              app->glscene->shapes));
      set_material(
          globject, get_element(ioobject->material, app->ioscene->materials,
                        app->glscene->materials));
      set_instance(
          globject, get_element(ioobject->instance, app->ioscene->instances,
                        app->glscene->instances));
    }
    end_glheader(win);
  }
  if (!app->ioscene->shapes.empty() && begin_glheader(win, "shapes")) {
    draw_glcombobox(win, "shape##2", app->selected_shape, app->ioscene->shapes);
    if (!draw_glwidgets(win, app->ioscene, app->selected_shape)) {
      auto ioshape = app->selected_shape;
      auto glshape = get_element(
          ioshape, app->ioscene->shapes, app->glscene->shapes);
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
  if (!app->ioscene->materials.empty() && begin_glheader(win, "materials")) {
    draw_glcombobox(
        win, "material##2", app->selected_material, app->ioscene->materials);
    if (draw_glwidgets(win, app->ioscene, app->selected_material)) {
      auto iomaterial = app->selected_material;
      auto glmaterial = get_element(
          iomaterial, app->ioscene->materials, app->glscene->materials);
      set_emission(glmaterial, iomaterial->emission,
          get_texture(iomaterial->emission_tex));
      set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
          get_texture(iomaterial->color_tex));
      set_specular(glmaterial,
          (1 - iomaterial->transmission) * iomaterial->specular,
          get_texture(iomaterial->specular_tex));
      set_metallic(glmaterial,
          (1 - iomaterial->transmission) * iomaterial->metallic,
          get_texture(iomaterial->metallic_tex));
      set_roughness(glmaterial, iomaterial->roughness,
          get_texture(iomaterial->roughness_tex));
      set_opacity(glmaterial, iomaterial->opacity,
          get_texture(iomaterial->opacity_tex));
      set_normalmap(glmaterial, get_texture(iomaterial->normal_tex));
    }
    end_glheader(win);
  }
  if (!app->ioscene->instances.empty() && begin_glheader(win, "instances")) {
    draw_glcombobox(
        win, "instance##2", app->selected_instance, app->ioscene->instances);
    if (!draw_glwidgets(win, app->ioscene, app->selected_instance)) {
      auto ioinstance = app->selected_instance;
      auto glinstance = get_element(
          ioinstance, app->ioscene->instances, app->glscene->instances);
      set_frames(glinstance, ioinstance->frames);
    }
    end_glheader(win);
  }
  if (!app->ioscene->textures.empty() && begin_glheader(win, "textures")) {
    draw_glcombobox(
        win, "texture##2", app->selected_texture, app->ioscene->textures);
    if (draw_glwidgets(win, app->ioscene, app->selected_texture)) {
      auto iotexture = app->selected_texture;
      auto gltexture = get_element(
          iotexture, app->ioscene->textures, app->glscene->textures);
      if (!iotexture->colorf.empty()) {
        set_texture(gltexture, iotexture->colorf);
      } else if (!iotexture->colorb.empty()) {
        set_texture(gltexture, iotexture->colorb);
      } else if (!iotexture->scalarf.empty()) {
        set_texture(gltexture, iotexture->scalarf);
      } else if (!iotexture->scalarb.empty()) {
        set_texture(gltexture, iotexture->scalarb);
      }
    }
    end_glheader(win);
  }
  if (!app->ioscene->subdivs.empty() && begin_glheader(win, "subdivs")) {
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
}

// draw with shading
void draw(opengl_window* win, app_states* apps, const opengl_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app = apps->selected;
  draw_glscene(app->glscene, input.framebuffer_viewport, app->drawgl_prms);
}

// update
void update(opengl_window* win, app_states* apps) {
  auto is_ready = [](const future<void>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    auto progress_cb = [app](const string& message, int current, int total) {
      app->progress = (float)current / (float)total;
    };
    app->loader.get();
    if (app->loader_error.empty()) {
      init_glscene(app->glscene, app->ioscene, progress_cb);
      update_lights(app->glscene, app->ioscene);
      app->ok     = true;
      app->status = "ok";
    } else {
      app->status = "error";
      app->error  = app->loader_error;
    }
  }
}

int main(int argc, const char* argv[]) {
  // initialize app
  auto apps_guard = make_unique<app_states>();
  auto apps       = apps_guard.get();
  auto filenames  = vector<string>{};

  // parse command line
  auto cli = make_cli("yscnview", "views scenes inteactively");
  add_option(cli, "--camera", apps->drawgl_prms.camera, "Camera index.");
  add_option(cli, "--resolution,-r", apps->drawgl_prms.resolution,
      "Image resolution.");
  add_option(cli, "--eyelight/--no-eyelight", apps->drawgl_prms.eyelight,
      "Eyelight rendering.");
  add_option(cli, "scenes", filenames, "Scene filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename);

  auto win_guard = make_glwindow({1280 + 320, 720}, "yscnview", true);
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
    if (!apps->selected || !apps->selected->ok) return;
    auto app = apps->selected;

    // handle mouse and keyboard for navigation
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
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
  });

  // run ui
  run_ui(win);

  // clear
  clear_glwindow(win);

  // done
  return 0;
}
