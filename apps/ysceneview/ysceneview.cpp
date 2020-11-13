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
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>
#include <yocto_gui/yocto_window.h>
using namespace yocto;

#include <deque>

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
  shade_params drawgl_prms = {};

  // scene
  sceneio_scene*  ioscene  = new sceneio_scene{};
  sceneio_camera* iocamera = nullptr;

  // rendering state
  shade_scene*  glscene  = new shade_scene{};
  shade_camera* glcamera = nullptr;

  // editing
  sceneio_camera*      selected_camera      = nullptr;
  sceneio_instance*    selected_instance    = nullptr;
  sceneio_shape*       selected_shape       = nullptr;
  sceneio_material*    selected_material    = nullptr;
  sceneio_environment* selected_environment = nullptr;
  sceneio_texture*     selected_texture     = nullptr;

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
  shade_params drawgl_prms = {};

  // imgui
  gui_widgets widgets = {};

  // cleanup
  ~app_states() {
    for (auto state : states) delete state;
  }
};

void load_scene_async(
    app_states* apps, const string& filename, const string& camera_name = "") {
  auto app         = apps->states.emplace_back(new app_state{});
  app->filename    = filename;
  app->imagename   = replace_extension(filename, ".png");
  app->outname     = replace_extension(filename, ".edited.json");
  app->name        = path_filename(app->filename);
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

void init_glscene(shade_scene* glscene, sceneio_scene* ioscene,
    shade_camera*& glcamera, sceneio_camera* iocamera,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->materials.size() +
             (int)ioscene->textures.size() + (int)ioscene->shapes.size() +
             (int)ioscene->instances.size()};

  // init scene
  init_scene(glscene);

  // camera
  auto camera_map     = unordered_map<sceneio_camera*, shade_camera*>{};
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
  auto texture_map     = unordered_map<sceneio_texture*, shade_texture*>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto gltexture = add_texture(glscene);
    if (!iotexture->hdr.empty()) {
      set_texture(gltexture, iotexture->hdr);
    } else if (!iotexture->ldr.empty()) {
      set_texture(gltexture, iotexture->ldr);
    }
    texture_map[iotexture] = gltexture;
  }

  // material
  auto material_map     = unordered_map<sceneio_material*, shade_material*>{};
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
  auto shape_map     = unordered_map<sceneio_shape*, shade_shape*>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    auto glshape       = add_shape(glscene, ioshape->points, ioshape->lines,
        ioshape->triangles, ioshape->quads, ioshape->positions,
        ioshape->normals, ioshape->texcoords, ioshape->colors);
    shape_map[ioshape] = glshape;
  }

  // shapes
  for (auto ioinstance : ioscene->instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto instance = add_instance(glscene);
    set_frame(instance, ioinstance->frame);
    set_shape(instance, shape_map.at(ioinstance->shape));
    set_material(instance, material_map.at(ioinstance->material));
  }

  // environments
  for (auto ioenvironment : ioscene->environments) {
    auto environment = add_environment(glscene);
    set_frame(environment, ioenvironment->frame);
    set_emission(environment, ioenvironment->emission,
        texture_map.at(ioenvironment->emission_tex));
  }

  // init environments
  init_environments(glscene);

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);

  // get cmmera
  glcamera = camera_map.at(iocamera);
}

bool draw_widgets(
    gui_widgets* widgets, sceneio_scene* ioscene, sceneio_camera* iocamera) {
  if (!iocamera) return false;
  auto edited = 0;
  draw_label(widgets, "name", iocamera->name);
  edited += draw_slider(widgets, "frame.x", iocamera->frame.x, -1, 1);
  edited += draw_slider(widgets, "frame.y", iocamera->frame.y, -1, 1);
  edited += draw_slider(widgets, "frame.z", iocamera->frame.z, -1, 1);
  edited += draw_slider(widgets, "frame.o", iocamera->frame.o, -10, 10);
  edited += draw_checkbox(widgets, "ortho", iocamera->orthographic);
  edited += draw_slider(widgets, "lens", iocamera->lens, 0.01f, 1);
  edited += draw_slider(widgets, "film", iocamera->film, 0.01f, 0.1f);
  edited += draw_slider(widgets, "focus", iocamera->focus, 0.01f, 1000);
  edited += draw_slider(widgets, "aperture", iocamera->aperture, 0, 5);
  auto from         = iocamera->frame.o,
       to           = iocamera->frame.o - iocamera->focus * iocamera->frame.z;
  auto from_changed = draw_slider(widgets, "!!from", from, -10, 10);
  auto to_changed   = draw_slider(widgets, "!!to", to, -10, 10);
  if (from_changed || to_changed) {
    iocamera->frame = lookat_frame(from, to, {0, 1, 0});
    iocamera->focus = length(from - to);
    edited += 1;
  }
  return edited;
}

/// Visit struct elements.
bool draw_widgets(
    gui_widgets* widgets, sceneio_scene* ioscene, sceneio_texture* iotexture) {
  if (!iotexture) return false;
  draw_label(widgets, "name", iotexture->name);
  draw_label(widgets, "hdr",
      std::to_string(iotexture->hdr.width()) + " x " +
          std::to_string(iotexture->hdr.height()));
  draw_label(widgets, "ldr",
      std::to_string(iotexture->ldr.width()) + " x " +
          std::to_string(iotexture->ldr.height()));
  return false;
}

bool draw_widgets(gui_widgets* widgets, sceneio_scene* ioscene,
    sceneio_material* iomaterial) {
  if (!iomaterial) return false;
  auto edited = 0;
  draw_label(widgets, "name", iomaterial->name);
  edited += draw_hdrcoloredit(widgets, "emission", iomaterial->emission);
  edited += draw_coloredit(widgets, "color", iomaterial->color);
  edited += draw_slider(widgets, "opacity", iomaterial->opacity, 0, 1);
  edited += draw_slider(widgets, "metallic", iomaterial->metallic, 0, 1);
  edited += draw_slider(widgets, "roughness", iomaterial->roughness, 0, 1);
  edited += draw_slider(widgets, "specular", iomaterial->specular, 0, 1);
  edited += draw_slider(widgets, "coat", iomaterial->coat, 0, 1);
  edited += draw_slider(
      widgets, "transmission", iomaterial->transmission, 0, 1);
  edited += draw_coloredit(widgets, "spectint", iomaterial->spectint);
  edited += draw_checkbox(widgets, "thin", iomaterial->thin);
  edited += draw_coloredit(widgets, "scattering", iomaterial->scattering);
  edited += draw_slider(widgets, "trdepth", iomaterial->trdepth, 0, 1);
  edited += draw_slider(
      widgets, "scanisotropy", iomaterial->scanisotropy, -1, 1);
  edited += draw_combobox(widgets, "emission_tex", iomaterial->emission_tex,
      ioscene->textures, true);
  edited += draw_combobox(
      widgets, "color_tex", iomaterial->color_tex, ioscene->textures, true);
  edited += draw_combobox(
      widgets, "opacity_tex", iomaterial->opacity_tex, ioscene->textures, true);
  edited += draw_combobox(widgets, "metallic_tex", iomaterial->metallic_tex,
      ioscene->textures, true);
  edited += draw_combobox(widgets, "roughness_tex", iomaterial->roughness_tex,
      ioscene->textures, true);
  edited += draw_combobox(widgets, "specular_tex", iomaterial->specular_tex,
      ioscene->textures, true);
  edited += draw_combobox(widgets, "transmission_tex",
      iomaterial->transmission_tex, ioscene->textures, true);
  edited += draw_combobox(widgets, "scattering_tex", iomaterial->scattering_tex,
      ioscene->textures, true);
  edited += draw_combobox(widgets, "spectint_tex", iomaterial->spectint_tex,
      ioscene->textures, true);
  edited += draw_combobox(
      widgets, "normal_tex", iomaterial->normal_tex, ioscene->textures, true);
  return edited;
}

bool draw_widgets(
    gui_widgets* widgets, sceneio_scene* ioscene, sceneio_shape* ioshape) {
  if (!ioshape) return false;
  auto edited = 0;
  draw_label(widgets, "name", ioshape->name);
  draw_label(widgets, "points", std::to_string(ioshape->points.size()));
  draw_label(widgets, "lines", std::to_string(ioshape->lines.size()));
  draw_label(widgets, "triangles", std::to_string(ioshape->triangles.size()));
  draw_label(widgets, "quads", std::to_string(ioshape->quads.size()));
  draw_label(widgets, "positions", std::to_string(ioshape->positions.size()));
  draw_label(widgets, "normals", std::to_string(ioshape->normals.size()));
  draw_label(widgets, "texcoords", std::to_string(ioshape->texcoords.size()));
  draw_label(widgets, "colors", std::to_string(ioshape->colors.size()));
  draw_label(widgets, "radius", std::to_string(ioshape->radius.size()));
  draw_label(widgets, "tangents", std::to_string(ioshape->tangents.size()));
  draw_label(widgets, "quads pos", std::to_string(ioshape->quadspos.size()));
  draw_label(widgets, "quads norm", std::to_string(ioshape->quadsnorm.size()));
  draw_label(
      widgets, "quads texcoord", std::to_string(ioshape->quadstexcoord.size()));
  edited += draw_slider(widgets, "subdivisions", ioshape->subdivisions, 0, 5);
  edited += draw_checkbox(widgets, "catmull-clark", ioshape->catmullclark);
  edited += draw_slider(widgets, "displacement", ioshape->displacement, 0, 1);
  edited += draw_combobox(widgets, "displacement_tex",
      ioshape->displacement_tex, ioscene->textures, true);
  return edited;
}

bool draw_widgets(
    gui_widgets* widgets, sceneio_scene* ioscene, sceneio_instance* ioobject) {
  if (!ioobject) return false;
  auto edited = 0;
  draw_label(widgets, "name", ioobject->name);
  edited += draw_slider(widgets, "frame.x", ioobject->frame.x, -1, 1);
  edited += draw_slider(widgets, "frame.y", ioobject->frame.y, -1, 1);
  edited += draw_slider(widgets, "frame.z", ioobject->frame.z, -1, 1);
  edited += draw_slider(widgets, "frame.o", ioobject->frame.o, -10, 10);
  edited += draw_combobox(widgets, "shape", ioobject->shape, ioscene->shapes);
  edited += draw_combobox(
      widgets, "material", ioobject->material, ioscene->materials);
  return edited;
}

bool draw_widgets(gui_widgets* widgets, sceneio_scene* ioscene,
    sceneio_environment* ioenvironment) {
  if (!ioenvironment) return false;
  auto edited = 0;
  edited += draw_textinput(widgets, "name", ioenvironment->name);
  edited += draw_slider(widgets, "frame.x", ioenvironment->frame.x, -1, 1);
  edited += draw_slider(widgets, "frame.y", ioenvironment->frame.y, -1, 1);
  edited += draw_slider(widgets, "frame.z", ioenvironment->frame.z, -1, 1);
  edited += draw_slider(widgets, "frame.o", ioenvironment->frame.o, -10, 10);
  edited += draw_hdrcoloredit(widgets, "emission", ioenvironment->emission);
  edited += draw_combobox(widgets, "emission texture",
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
  print_fatal("element not found");
  return nullptr;
}

// draw with shading
void draw_widgets(app_states* apps, const gui_input& input) {
  auto widgets = &apps->widgets;
  begin_imgui(widgets, "ysceneview");

  static auto load_path = ""s, save_path = ""s, error_message = ""s;
  if (draw_filedialog_button(widgets, "load", true, "load", load_path, false,
          "./", "", "*.json;*.obj;*.pbrt")) {
    load_scene_async(apps, load_path);
    load_path = "";
  }
  continue_line(widgets);
  if (draw_filedialog_button(widgets, "save",
          apps->selected && apps->selected->ok, "save", save_path, true,
          path_dirname(save_path), path_filename(save_path),
          "*.json;*.obj;*.pbrt")) {
    auto app     = apps->selected;
    app->outname = save_path;
    save_scene(app->outname, app->ioscene, app->error);
    save_path = "";
  }
  continue_line(widgets);
  if (draw_button(widgets, "close", (bool)apps->selected)) {
    if (apps->selected->loader.valid()) {
      end_imgui(widgets);
      return;
    }
    delete apps->selected;
    apps->states.erase(
        std::find(apps->states.begin(), apps->states.end(), apps->selected));
    apps->selected = apps->states.empty() ? nullptr : apps->states.front();
  }
  continue_line(widgets);
  if (draw_button(widgets, "quit")) {
    set_close(widgets->window, true);
  }
  if (apps->states.empty()) {
    end_imgui(widgets);
    return;
  }
  draw_combobox(widgets, "scene", apps->selected, apps->states, false);
  if (!apps->selected) {
    end_imgui(widgets);
    return;
  }
  draw_progressbar(widgets, apps->selected->status.c_str(),
      apps->selected->current, apps->selected->total);
  if (apps->selected->error != "") {
    draw_label(widgets, "error", apps->selected->error);
    end_imgui(widgets);
    return;
  }
  if (!apps->selected->ok) {
    end_imgui(widgets);
    return;
  }
  auto app = apps->selected;
  if (begin_header(widgets, "view")) {
    if (draw_combobox(
            widgets, "camera", app->iocamera, app->ioscene->cameras)) {
      app->glcamera = get_element(
          app->iocamera, app->ioscene->cameras, app->glscene->cameras);
    }
    auto& params = app->drawgl_prms;
    draw_slider(widgets, "resolution", params.resolution, 0, 4096);
    draw_combobox(
        widgets, "lighting", (int&)params.lighting, shade_lighting_names);
    draw_checkbox(widgets, "wireframe", params.wireframe);
    continue_line(widgets);
    draw_checkbox(widgets, "faceted", params.faceted);
    continue_line(widgets);
    draw_checkbox(widgets, "double sided", params.double_sided);
    draw_slider(widgets, "exposure", params.exposure, -10, 10);
    draw_slider(widgets, "gamma", params.gamma, 0.1f, 4);
    draw_slider(widgets, "near", params.near, 0.01f, 1.0f);
    draw_slider(widgets, "far", params.far, 1000.0f, 10000.0f);
    end_header(widgets);
  }
  if (begin_header(widgets, "inspect")) {
    draw_label(widgets, "scene", path_filename(app->filename));
    draw_label(widgets, "filename", app->filename);
    draw_label(widgets, "outname", app->outname);
    draw_label(widgets, "imagename", app->imagename);
    continue_line(widgets);
    if (draw_button(widgets, "print cams")) {
      for (auto iocamera : app->ioscene->cameras) {
        print_obj_camera(iocamera);
      }
    }
    continue_line(widgets);
    if (draw_button(widgets, "print stats")) {
      for (auto stat : scene_stats(app->ioscene)) print_info(stat);
    }
    end_header(widgets);
  }
  auto get_texture = [app](sceneio_texture* iotexture) {
    return get_element(
        iotexture, app->ioscene->textures, app->glscene->textures);
  };
  if (!app->ioscene->cameras.empty() && begin_header(widgets, "cameras")) {
    draw_combobox(
        widgets, "camera##2", app->selected_camera, app->ioscene->cameras);
    if (draw_widgets(widgets, app->ioscene, app->selected_camera)) {
      auto iocamera = app->selected_camera;
      auto glcamera = get_element(
          iocamera, app->ioscene->cameras, app->glscene->cameras);
      set_frame(glcamera, iocamera->frame);
      set_lens(glcamera, iocamera->lens, iocamera->aspect, iocamera->film);
      set_nearfar(glcamera, 0.001, 10000);
    }
    end_header(widgets);
  }
  if (!app->ioscene->environments.empty() &&
      begin_header(widgets, "environments")) {
    draw_combobox(widgets, "environments##2", app->selected_environment,
        app->ioscene->environments);
    if (draw_widgets(widgets, app->ioscene, app->selected_environment)) {
      auto ioenvironment = app->selected_environment;
      auto glenvironment = get_element(ioenvironment,
          app->ioscene->environments, app->glscene->environments);
      set_emission(glenvironment, ioenvironment->emission);
    }
    end_header(widgets);
  }
  if (!app->ioscene->instances.empty() && begin_header(widgets, "instances")) {
    draw_combobox(widgets, "instance##2", app->selected_instance,
        app->ioscene->instances);
    if (draw_widgets(widgets, app->ioscene, app->selected_instance)) {
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
    end_header(widgets);
  }
  if (!app->ioscene->shapes.empty() && begin_header(widgets, "shapes")) {
    draw_combobox(
        widgets, "shape##2", app->selected_shape, app->ioscene->shapes);
    if (!draw_widgets(widgets, app->ioscene, app->selected_shape)) {
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
    end_header(widgets);
  }
  if (!app->ioscene->materials.empty() && begin_header(widgets, "materials")) {
    draw_combobox(widgets, "material##2", app->selected_material,
        app->ioscene->materials);
    if (draw_widgets(widgets, app->ioscene, app->selected_material)) {
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
    end_header(widgets);
  }
  if (!app->ioscene->textures.empty() && begin_header(widgets, "textures")) {
    draw_combobox(
        widgets, "texture##2", app->selected_texture, app->ioscene->textures);
    if (draw_widgets(widgets, app->ioscene, app->selected_texture)) {
      auto iotexture = app->selected_texture;
      auto gltexture = get_element(
          iotexture, app->ioscene->textures, app->glscene->textures);
      if (!iotexture->hdr.empty()) {
        set_texture(gltexture, iotexture->hdr);
      } else if (!iotexture->ldr.empty()) {
        set_texture(gltexture, iotexture->ldr);
      }
    }
    end_header(widgets);
  }
  end_imgui(widgets);
}

void drop(app_states* apps, const gui_input& input) {
  if (input.dropped.size()) {
    for (auto& path : input.dropped) load_scene_async(apps, path);
  }
}

// draw with shading
void draw(app_states* apps, const gui_input& input) {
  if (!apps->selected || !apps->selected->ok) return;
  auto app = apps->selected;
  draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
      app->drawgl_prms);
}

// update
void update(app_states* apps) {
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
      app->ok     = true;
      app->status = "ok";
    } else {
      app->status = "error";
      app->error  = app->loader_error;
    }
  }
}

void update_camera(app_states* apps, const gui_input& input) {
  if (is_active(&apps->widgets)) return;

  if (!apps->selected || !apps->selected->ok) return;
  auto app = apps->selected;

  // handle mouse and keyboard for navigation
  if ((input.mouse_left || input.mouse_right) && !input.modifier_alt) {
    auto dolly  = 0.0f;
    auto pan    = zero2f;
    auto rotate = zero2f;
    if (input.mouse_left && !input.modifier_shift)
      rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
    if (input.mouse_right)
      dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
    if (input.mouse_left && input.modifier_shift)
      pan = (input.mouse_pos - input.mouse_last) / 100.0f;
    std::tie(app->iocamera->frame, app->iocamera->focus) = camera_turntable(
        app->iocamera->frame, app->iocamera->focus, rotate, dolly, pan);
    set_frame(app->glcamera, app->iocamera->frame);
  }
};

void update_app(const gui_input& input, void* data) {
  auto apps = (app_states*)data;

  update_camera(apps, input);
  drop(apps, input);
  update(apps);

  draw(apps, input);
  draw_widgets(apps, input);
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
  add_option(cli, "--lighting", apps->drawgl_prms.lighting, "Lighting type.",
      shade_lighting_names);
  add_option(cli, "scenes", filenames, "Scene filenames", true);
  parse_cli(cli, argc, argv);

  // loading images
  for (auto filename : filenames) load_scene_async(apps, filename, camera_name);

  auto window = new gui_window{};
  init_window(window, {1280 + 320, 720}, "ysceneviews", true);
  window->user_data = apps;
  apps->widgets     = create_imgui(window);

  // run ui
  run_ui(window, update_app);
  for (auto app : apps->states) clear_scene(app->glscene);

  // done
  return 0;
}
