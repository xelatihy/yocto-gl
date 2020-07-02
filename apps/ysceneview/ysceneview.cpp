//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_gui.h>
using namespace yocto;

#include <atomic>
#include <deque>
#include <future>

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto {
void print_obj_camera(scene_camera* camera);
};

// Application state
struct app_state {
  // loading parameters
  string filename  = "scene.json";
  string imagename = "out.png";
  string outname   = "scene.json";
  string name      = "";

  // options
  ogl_scene_params drawgl_prms = {};

  // scene
  scene_model*  ioscene  = new scene_model{};
  scene_camera* iocamera = nullptr;

  // rendering state
  ogl_scene*  glscene  = new ogl_scene{};
  ogl_camera* glcamera = nullptr;

  // editing
  scene_camera*      selected_camera      = nullptr;
  scene_instance*    selected_instance    = nullptr;
  scene_shape*       selected_shape       = nullptr;
  scene_material*    selected_material    = nullptr;
  scene_environment* selected_environment = nullptr;
  scene_texture*     selected_texture     = nullptr;

  // loading status
  std::atomic<bool> ok           = false;
  std::future<void> loader       = {};
  string            status       = "";
  string            error        = "";
  std::atomic<int>  current      = 0;
  std::atomic<int>  total        = 0;
  string            loader_error = "";

  ~app_state() {
    if (ioscene) delete ioscene;
    if (glscene) delete glscene;
  }
};

// Application state
struct app_states {
  // data
  vector<app_state*>     states   = {};
  app_state*             selected = nullptr;
  std::deque<app_state*> loading  = {};

  // default options
  ogl_scene_params drawgl_prms = {};

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

void load_scene_async(
    app_states* apps, const string& filename, const string& camera_name = "") {
  auto app         = apps->states.emplace_back(new app_state{});
  app->filename    = filename;
  app->imagename   = sfs::path(filename).replace_extension(".png");
  app->outname     = sfs::path(filename).replace_extension(".edited.yaml");
  app->name        = sfs::path(app->filename).filename();
  app->drawgl_prms = apps->drawgl_prms;
  app->status      = "load";
  app->loader      = std::async(std::launch::async, [app, camera_name]() {
    auto progress_cb = [app](const string& message, int current, int total) {
      app->current = current;
      app->total   = total;
    };
    if (!load_scene(
            app->filename, app->ioscene, app->loader_error, progress_cb))
      return;
    app->iocamera = get_camera(app->ioscene, camera_name);
    tesselate_shapes(app->ioscene, progress_cb);
  });
  apps->loading.push_back(app);
  if (!apps->selected) apps->selected = app;
}

void update_lights(ogl_scene* glscene, scene_model* ioscene) {
  clear_lights(glscene);
  for (auto ioobject : ioscene->instances) {
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
    set_light(add_light(glscene), transform_point(ioobject->frame, pos), ke,
        ogl_light_type::point, false);
  }
}

void init_glscene(ogl_scene* glscene, scene_model* ioscene,
    ogl_camera*& glcamera, scene_camera* iocamera,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->materials.size() +
             (int)ioscene->textures.size() + (int)ioscene->shapes.size() +
             (int)ioscene->instances.size()};

  // create scene
  init_scene(glscene);

  // camera
  auto camera_map     = unordered_map<scene_camera*, ogl_camera*>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto camera = add_camera(glscene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_nearfar(camera, 0.001, 10000);
    camera_map[iocamera] = camera;
  }

  // textures
  auto texture_map     = unordered_map<scene_texture*, ogl_texture*>{};
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
  auto material_map     = unordered_map<scene_material*, ogl_material*>{};
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

  // shapes
  auto shape_map     = unordered_map<scene_shape*, ogl_shape*>{};
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
    set_edges(glshape, ioshape->triangles, ioshape->quads);
    shape_map[ioshape] = glshape;
  }

  // shapes
  for (auto ioinstance : ioscene->instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto globject = add_object(glscene);
    set_frame(globject, ioinstance->frame);
    set_shape(globject, shape_map.at(ioinstance->shape));
    set_material(globject, material_map.at(ioinstance->material));
  }

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);

  // get cmmera
  glcamera = camera_map.at(iocamera);
}

bool draw_widgets(
    gui_window* win, scene_model* ioscene, scene_camera* iocamera) {
  if (!iocamera) return false;
  auto edited = 0;
  draw_label(win, "name", iocamera->name);
  edited += draw_slider(win, "frame.x", iocamera->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", iocamera->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", iocamera->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", iocamera->frame.o, -10, 10);
  edited += draw_checkbox(win, "ortho", iocamera->orthographic);
  edited += draw_slider(win, "lens", iocamera->lens, 0.01f, 1);
  edited += draw_slider(win, "film", iocamera->film, 0.01f, 0.1f);
  edited += draw_slider(win, "focus", iocamera->focus, 0.01f, 1000);
  edited += draw_slider(win, "aperture", iocamera->aperture, 0, 5);
  auto from         = iocamera->frame.o,
       to           = iocamera->frame.o - iocamera->focus * iocamera->frame.z;
  auto from_changed = draw_slider(win, "!!from", from, -10, 10);
  auto to_changed   = draw_slider(win, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    iocamera->frame = lookat_frame(from, to, {0, 1, 0});
    iocamera->focus = length(from - to);
    edited += 1;
  }
  return edited;
}

/// Visit struct elements.
bool draw_widgets(
    gui_window* win, scene_model* ioscene, scene_texture* iotexture) {
  if (!iotexture) return false;
  draw_label(win, "name", iotexture->name);
  draw_label(win, "colorf",
      std::to_string(iotexture->colorf.imsize().x) + " x " +
          std::to_string(iotexture->colorf.imsize().y));
  draw_label(win, "colorb",
      std::to_string(iotexture->colorb.imsize().x) + " x " +
          std::to_string(iotexture->colorb.imsize().y));
  draw_label(win, "scalarf",
      std::to_string(iotexture->scalarf.imsize().x) + " x " +
          std::to_string(iotexture->scalarf.imsize().y));
  draw_label(win, "scalarb",
      std::to_string(iotexture->scalarb.imsize().x) + " x " +
          std::to_string(iotexture->scalarb.imsize().y));
  return false;
}

bool draw_widgets(
    gui_window* win, scene_model* ioscene, scene_material* iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
  draw_label(win, "name", iomaterial->name);
  edited += draw_hdrcoloredit(win, "emission", iomaterial->emission);
  edited += draw_coloredit(win, "color", iomaterial->color);
  edited += draw_slider(win, "opacity", iomaterial->opacity, 0, 1);
  edited += draw_slider(win, "metallic", iomaterial->metallic, 0, 1);
  edited += draw_slider(win, "roughness", iomaterial->roughness, 0, 1);
  edited += draw_slider(win, "specular", iomaterial->specular, 0, 1);
  edited += draw_slider(win, "coat", iomaterial->coat, 0, 1);
  edited += draw_slider(win, "transmission", iomaterial->transmission, 0, 1);
  edited += draw_coloredit(win, "spectint", iomaterial->spectint);
  edited += draw_checkbox(win, "thin", iomaterial->thin);
  edited += draw_coloredit(win, "scattering", iomaterial->scattering);
  edited += draw_slider(win, "trdepth", iomaterial->trdepth, 0, 1);
  edited += draw_slider(win, "scanisotropy", iomaterial->scanisotropy, -1, 1);
  edited += draw_combobox(
      win, "emission_tex", iomaterial->emission_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "color_tex", iomaterial->color_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "opacity_tex", iomaterial->opacity_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "metallic_tex", iomaterial->metallic_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "roughness_tex", iomaterial->roughness_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "specular_tex", iomaterial->specular_tex, ioscene->textures, true);
  edited += draw_combobox(win, "transmission_tex", iomaterial->transmission_tex,
      ioscene->textures, true);
  edited += draw_combobox(win, "scattering_tex", iomaterial->scattering_tex,
      ioscene->textures, true);
  edited += draw_combobox(
      win, "spectint_tex", iomaterial->spectint_tex, ioscene->textures, true);
  edited += draw_combobox(
      win, "normal_tex", iomaterial->normal_tex, ioscene->textures, true);
  return edited;
}

bool draw_widgets(gui_window* win, scene_model* ioscene, scene_shape* ioshape) {
  if (!ioshape) return false;
  auto edited = 0;
  draw_label(win, "name", ioshape->name);
  draw_label(win, "points", std::to_string(ioshape->points.size()));
  draw_label(win, "lines", std::to_string(ioshape->lines.size()));
  draw_label(win, "triangles", std::to_string(ioshape->triangles.size()));
  draw_label(win, "quads", std::to_string(ioshape->quads.size()));
  draw_label(win, "positions", std::to_string(ioshape->positions.size()));
  draw_label(win, "normals", std::to_string(ioshape->normals.size()));
  draw_label(win, "texcoords", std::to_string(ioshape->texcoords.size()));
  draw_label(win, "colors", std::to_string(ioshape->colors.size()));
  draw_label(win, "radius", std::to_string(ioshape->radius.size()));
  draw_label(win, "tangents", std::to_string(ioshape->tangents.size()));
  draw_label(win, "quads pos", std::to_string(ioshape->quadspos.size()));
  draw_label(win, "quads norm", std::to_string(ioshape->quadsnorm.size()));
  draw_label(
      win, "quads texcoord", std::to_string(ioshape->quadstexcoord.size()));
  edited += draw_slider(win, "subdivisions", ioshape->subdivisions, 0, 5);
  edited += draw_checkbox(win, "catmull-clark", ioshape->catmullclark);
  edited += draw_slider(win, "displacement", ioshape->displacement, 0, 1);
  edited += draw_combobox(win, "displacement_tex", ioshape->displacement_tex,
      ioscene->textures, true);
  return edited;
}

bool draw_widgets(
    gui_window* win, scene_model* ioscene, scene_instance* ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  draw_label(win, "name", ioobject->name);
  edited += draw_slider(win, "frame.x", ioobject->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", ioobject->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", ioobject->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_combobox(win, "shape", ioobject->shape, ioscene->shapes);
  edited += draw_combobox(
      win, "material", ioobject->material, ioscene->materials);
  return edited;
}

bool draw_widgets(
    gui_window* win, scene_model* ioscene, scene_environment* ioenvironment) {
  if (!ioenvironment) return false;
  auto edited = 0;
  edited += draw_textinput(win, "name", ioenvironment->name);
  edited += draw_slider(win, "frame.x", ioenvironment->frame.x, -1, 1);
  edited += draw_slider(win, "frame.y", ioenvironment->frame.y, -1, 1);
  edited += draw_slider(win, "frame.z", ioenvironment->frame.z, -1, 1);
  edited += draw_slider(win, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_hdrcoloredit(win, "emission", ioenvironment->emission);
  edited += draw_combobox(win, "emission texture", ioenvironment->emission_tex,
      ioscene->textures, true);
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
void draw_widgets(gui_window* win, app_states* apps, const gui_input& input) {
  static auto load_path = ""s, save_path = ""s, error_message = ""s;
  if (draw_filedialog_button(win, "load", true, "load", load_path, false, "./",
          "", "*.yaml;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_line(win);
  if (draw_filedialog_button(win, "save", apps->selected && apps->selected->ok,
          "save", save_path, true, sfs::path(save_path).parent_path(),
          sfs::path(save_path).filename(), "*.yaml;*.obj;*.pbrt")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_scene(app->outname, app->ioscene, app->error);
    save_path = "";
  }
  continue_line(win);
  if (draw_button(win, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) return;
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_line(win);
  if (draw_button(win, "quit")) {
    set_close(win, true);
  }
  if (apps->states.empty()) return;
  draw_combobox(win, "scene", apps->selected, apps->states, false);
  if (!apps->selected) return;
  draw_progressbar(win, apps->selected->status.c_str(), apps->selected->current,
      apps->selected->total);
  if (apps->selected->error != "") {
    draw_label(win, "error", apps->selected->error);
    return;
  }
  if (!apps->selected->ok) return;
  auto app = apps->selected;
  if (begin_header(win, "view")) {
    if (draw_combobox(win, "camera", app->iocamera, app->ioscene->cameras)) {
      app->glcamera = get_element(
          app->iocamera, app->ioscene->cameras, app->glscene->cameras);
    }
    auto& params = app->drawgl_prms;
    draw_slider(win, "resolution", params.resolution, 0, 4096);
    draw_combobox(win, "shading", (int&)params.shading, ogl_shading_names);
    draw_checkbox(win, "wireframe", params.wireframe);
    continue_line(win);
    draw_checkbox(win, "edges", params.edges);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
    end_header(win);
  }
  if (begin_header(win, "inspect")) {
    draw_label(win, "scene", sfs::path(app->filename).filename());
    draw_label(win, "filename", app->filename);
    draw_label(win, "outname", app->outname);
    draw_label(win, "imagename", app->imagename);
    continue_line(win);
    if (draw_button(win, "print cams")) {
      for (auto iocamera : app->ioscene->cameras) {
        print_obj_camera(iocamera);
      }
    }
    continue_line(win);
    if (draw_button(win, "print stats")) {
      for (auto stat : scene_stats(app->ioscene)) print_info(stat);
    }
    end_header(win);
  }
  auto get_texture = [app](scene_texture* iotexture) {
    return get_element(
        iotexture, app->ioscene->textures, app->glscene->textures);
  };
  if (!app->ioscene->cameras.empty() && begin_header(win, "cameras")) {
    draw_combobox(
        win, "camera##2", app->selected_camera, app->ioscene->cameras);
    if (draw_widgets(win, app->ioscene, app->selected_camera)) {
      auto iocamera = app->selected_camera;
      auto glcamera = get_element(
          iocamera, app->ioscene->cameras, app->glscene->cameras);
      set_frame(glcamera, iocamera->frame);
      set_lens(glcamera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_nearfar(glcamera, 0.001, 10000);
    }
    end_header(win);
  }
  if (!app->ioscene->environments.empty() &&
      begin_header(win, "environments")) {
    draw_combobox(win, "environments##2", app->selected_environment,
        app->ioscene->environments);
    if (draw_widgets(win, app->ioscene, app->selected_environment)) {
    }
    end_header(win);
  }
  if (!app->ioscene->instances.empty() && begin_header(win, "objects")) {
    draw_combobox(
        win, "instance##2", app->selected_instance, app->ioscene->instances);
    if (!draw_widgets(win, app->ioscene, app->selected_instance)) {
      auto ioobject = app->selected_instance;
      auto globject = get_element(
          ioobject, app->ioscene->instances, app->glscene->instances);
      set_frame(globject, ioobject->frame);
      set_shape(globject, get_element(ioobject->shape, app->ioscene->shapes,
                              app->glscene->shapes));
      set_material(
          globject, get_element(ioobject->material, app->ioscene->materials,
                        app->glscene->materials));
    }
    end_header(win);
  }
  if (!app->ioscene->shapes.empty() && begin_header(win, "shapes")) {
    draw_combobox(win, "shape##2", app->selected_shape, app->ioscene->shapes);
    if (!draw_widgets(win, app->ioscene, app->selected_shape)) {
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
    end_header(win);
  }
  if (!app->ioscene->materials.empty() && begin_header(win, "materials")) {
    draw_combobox(
        win, "material##2", app->selected_material, app->ioscene->materials);
    if (draw_widgets(win, app->ioscene, app->selected_material)) {
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
    end_header(win);
  }
  if (!app->ioscene->textures.empty() && begin_header(win, "textures")) {
    draw_combobox(
        win, "texture##2", app->selected_texture, app->ioscene->textures);
    if (draw_widgets(win, app->ioscene, app->selected_texture)) {
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
    end_header(win);
  }
}

// draw with shading
void draw(gui_window* win, app_states* apps, const gui_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app = apps->selected;
  if (app->drawgl_prms.shading == ogl_shading_type::lights)
    update_lights(app->glscene, app->ioscene);
  draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
      app->drawgl_prms);
}

// update
void update(gui_window* win, app_states* apps) {
  auto is_ready = [](const std::future<void>& result) -> bool {
    return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                                 std::future_status::ready;
  };

  while (!apps->loading.empty()) {
    auto app = apps->loading.front();
    if (!is_ready(app->loader)) break;
    apps->loading.pop_front();
    auto progress_cb = [app](const string& message, int current, int total) {
      app->current = current;
      app->total   = total;
    };
    app->loader.get();
    if (app->loader_error.empty()) {
      init_glscene(app->glscene, app->ioscene, app->glcamera, app->iocamera,
          progress_cb);
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
  auto apps_guard  = std::make_unique<app_states>();
  auto apps        = apps_guard.get();
  auto filenames   = vector<string>{};
  auto camera_name = ""s;

  // parse command line
  auto cli = make_cli("ysceneview", "views scenes inteactively");
  add_option(cli, "--camera", camera_name, "Camera name.");
  add_option(cli, "--resolution,-r", apps->drawgl_prms.resolution,
      "Image resolution.");
  add_option(cli, "--shading", apps->drawgl_prms.shading, "Shading type.",
      ogl_shading_names);
  add_option(cli, "scenes", filenames, "Scene filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename, camera_name);

  // callbacks
  auto callbacks     = gui_callbacks{};
  callbacks.clear_cb = [apps](gui_window* win, const gui_input& input) {
    for (auto app : apps->states) clear_scene(app->glscene);
  };
  callbacks.draw_cb = [apps](gui_window* win, const gui_input& input) {
    draw(win, apps, input);
  };
  callbacks.widgets_cb = [apps](gui_window* win, const gui_input& input) {
    draw_widgets(win, apps, input);
  };
  callbacks.drop_cb = [apps](gui_window* win, const vector<string>& paths,
                          const gui_input& input) {
    for (auto& path : paths) load_scene_async(apps, path);
  };
  callbacks.update_cb = [apps](gui_window* win, const gui_input& input) {
    update(win, apps);
  };
  callbacks.uiupdate_cb = [apps](gui_window* win, const gui_input& input) {
    if (!apps->selected || !apps->selected->ok) return;
    auto app = apps->selected;

    // handle mouse and keyboard for navigation
    if ((input.mouse_left || input.mouse_right) && !input.modifier_alt &&
        !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      update_turntable(
          app->iocamera->frame, app->iocamera->focus, rotate, dolly, pan);
      set_frame(app->glcamera, app->iocamera->frame);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "ysceneview", callbacks);

  // done
  return 0;
}
