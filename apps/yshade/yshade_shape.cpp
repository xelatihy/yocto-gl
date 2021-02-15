//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <yocto/yocto_cli.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_parallel.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>
using namespace yocto;

#include "yshade_shape.h"

#ifdef _WIN32
#undef near
#undef far
#endif

// Application state
struct shade_shape_state {
  // loading parameters
  string filename  = "shape.ply";
  string imagename = "out.png";
  string outname   = "out.ply";
  string name      = "";

  // options
  shade_params drawgl_prms = {};

  // scene
  sceneio_scene ioscene  = sceneio_scene{};
  camera_handle iocamera = invalid_handle;

  // rendering state
  shade_scene glscene = {};
};

// Create a shape with small spheres for each point
static shape_data make_spheres(
    const vector<vec3f>& positions, float radius, int steps) {
  auto shape = shape_data{};
  for (auto position : positions) {
    auto sphere = make_sphere(steps, radius);
    for (auto& p : sphere.positions) p += position;
    merge_quads(shape.quads, shape.positions, shape.normals, shape.texcoords,
        sphere.quads, sphere.positions, sphere.normals, sphere.texcoords);
  }
  return shape;
}
static shape_data make_cylinders(const vector<vec2i>& lines,
    const vector<vec3f>& positions, float radius, const vec3i& steps) {
  auto shape = shape_data{};
  for (auto line : lines) {
    auto len      = length(positions[line.x] - positions[line.y]);
    auto dir      = normalize(positions[line.x] - positions[line.y]);
    auto center   = (positions[line.x] + positions[line.y]) / 2;
    auto cylinder = make_uvcylinder({4, 1, 1}, {radius, len / 2});
    auto frame    = frame_fromz(center, dir);
    for (auto& p : cylinder.positions) p = transform_point(frame, p);
    for (auto& n : cylinder.normals) n = transform_direction(frame, n);
    merge_quads(shape.quads, shape.positions, shape.normals, shape.texcoords,
        cylinder.quads, cylinder.positions, cylinder.normals,
        cylinder.texcoords);
  }
  return shape;
}

static frame3f camera_frame(float lens, float aspect, float film = 0.036) {
  auto camera_dir  = normalize(vec3f{0, 0.5, 1});
  auto bbox_radius = 2.0f;
  auto camera_dist = bbox_radius * lens / (film / aspect);
  return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
}

static void convert_scene(scene_scene& scene, const scene_shape& ioshape_,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 5};
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);

  // init scene
  scene = {};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));
  // TODO(fabio): this should be a math function

  // camera
  if (progress_cb) progress_cb("create camera", progress.x++, progress.y);
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  if (progress_cb) progress_cb("create material", progress.x++, progress.y);
  auto& shape_material     = scene.materials.emplace_back();
  shape_material.type      = material_type::plastic;
  shape_material.color     = {0.5, 1, 0.5};
  shape_material.roughness = 0.2;
  auto& edges_material     = scene.materials.emplace_back();
  edges_material.type      = material_type::matte;
  edges_material.color     = {0, 0, 0};
  auto& vertices_material  = scene.materials.emplace_back();
  vertices_material.type   = material_type::matte;
  vertices_material.color  = {0, 0, 0};

  // shapes
  if (progress_cb) progress_cb("create shape", progress.x++, progress.y);
  scene.shapes.emplace_back(ioshape);
  auto& edges_shape        = scene.shapes.emplace_back();
  edges_shape.positions    = ioshape.positions;
  edges_shape.lines        = get_edges(ioshape.triangles, ioshape.quads);
  auto& vertices_shape     = scene.shapes.emplace_back(ioshape);
  vertices_shape.positions = ioshape.positions;
  vertices_shape.points    = vector<int>(ioshape.positions.size());
  for (auto idx = 0; idx < (int)vertices_shape.points.size(); idx++)
    vertices_shape.points[idx] = idx;

#if 0
  auto froms = vector<vec3f>();
  auto tos   = vector<vec3f>();
  froms.reserve(edges.size());
  tos.reserve(edges.size());
  float avg_edge_length = 0;
  for (auto& edge : edges) {
    auto from = ioshape.positions[edge.x];
    auto to   = ioshape.positions[edge.y];
    froms.push_back(from);
    tos.push_back(to);
    avg_edge_length += length(from - to);
  }
  avg_edge_length /= edges.size();
  auto cylinder_radius = 0.05f * avg_edge_length;
  auto cylinder        = make_uvcylinder({4, 1, 1}, {cylinder_radius, 1});
  for (auto& p : cylinder.positions) {
    p.z = p.z * 0.5 + 0.5;
  }
  auto edges_shape = add_shape(glscene, {}, {}, {}, cylinder.quads,
      cylinder.positions, cylinder.normals, cylinder.texcoords, {});
  set_instances(glscene.shapes[edges_shape], froms, tos);

  auto vertices_radius = 3.0f * cylinder_radius;
  auto vertices        = make_spheres(ioshape.positions, vertices_radius, 2);
  auto vertices_shape  = add_shape(glscene, {}, {}, {}, vertices.quads,
      vertices.positions, vertices.normals, vertices.texcoords, {});
  set_instances(glscene.shapes[vertices_shape], {}, {});
#endif

  // instances
  if (progress_cb) progress_cb("create instance", progress.x++, progress.y);
  auto& shape_instance       = scene.instances.emplace_back();
  shape_instance.shape       = 0;
  shape_instance.material    = 0;
  auto& edges_instance       = scene.instances.emplace_back();
  edges_instance.shape       = 1;
  edges_instance.material    = 1;
  auto& vertices_instance    = scene.instances.emplace_back();
  vertices_instance.shape    = 2;
  vertices_instance.material = 2;

  // done
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);
}

static void init_glscene(shade_scene& glscene, const sceneio_scene& ioscene,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene.cameras.size() + (int)ioscene.materials.size() +
             (int)ioscene.textures.size() + (int)ioscene.shapes.size() +
             (int)ioscene.instances.size()};

  // init scene
  init_scene(glscene);

  // camera
  for (auto& iocamera : ioscene.cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto& camera = glscene.cameras.at(add_camera(glscene));
    set_frame(camera, iocamera.frame);
    set_lens(camera, iocamera.lens, iocamera.aspect, iocamera.film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  for (auto& iotexture : ioscene.textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto  handle    = add_texture(glscene);
    auto& gltexture = glscene.textures[handle];
    if (!iotexture.pixelsf.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsf);
    } else if (!iotexture.pixelsb.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsb);
    }
  }

  // material
  for (auto& iomaterial : ioscene.materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto  handle     = add_material(glscene);
    auto& glmaterial = glscene.materials[handle];
    set_emission(glmaterial, iomaterial.emission, iomaterial.emission_tex);
    set_opacity(glmaterial, iomaterial.opacity, invalid_handle);
    set_normalmap(glmaterial, iomaterial.normal_tex);
    switch (iomaterial.type) {
      case material_type::matte: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, 0, invalid_handle);
      } break;
      case material_type::plastic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      case material_type::metal: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 1, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
      case material_type::metallic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, iomaterial.metallic, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      default: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
    }
  }

  // shapes
  for (auto& ioshape : ioscene.shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    add_shape(glscene, ioshape.points, ioshape.lines, ioshape.triangles,
        ioshape.quads, ioshape.positions, ioshape.normals, ioshape.texcoords,
        ioshape.colors);
  }

  // shapes
  for (auto& ioinstance : ioscene.instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto  handle     = add_instance(glscene);
    auto& glinstance = glscene.instances[handle];
    set_frame(glinstance, ioinstance.frame);
    set_shape(glinstance, ioinstance.shape);
    set_material(glinstance, ioinstance.material);
  }

  // environments
  for (auto& ioenvironment : ioscene.environments) {
    auto  handle        = add_environment(glscene);
    auto& glenvironment = glscene.environments[handle];
    set_frame(glenvironment, ioenvironment.frame);
    set_emission(
        glenvironment, ioenvironment.emission, ioenvironment.emission_tex);
  }

  // init environments
  init_environments(glscene);

  // done
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);
}

int run_shade_shape(const shade_shape_params& params) {
  // initialize app
  auto app = shade_shape_state();

  // copy command line
  app.filename = params.shape;

  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  print_progress("load shape", 0, 1);
  if (!load_shape(app.filename, ioshape, ioerror)) print_fatal(ioerror);
  print_progress("load shape", 1, 1);

  // create scene
  convert_scene(app.ioscene, ioshape, print_progress);

  // get camera
  app.iocamera = find_camera(app.ioscene, "");

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&app](gui_window* win, const gui_input& input) {
    init_glscene(app.glscene, app.ioscene, print_progress);
    app.glscene.instances[1].hidden = true;
    app.glscene.instances[2].hidden = true;
  };
  callbacks.clear_cb = [&app](gui_window* win, const gui_input& input) {
    clear_scene(app.glscene);
  };
  callbacks.draw_cb = [&app](gui_window* win, const gui_input& input) {
    draw_scene(app.glscene, app.glscene.cameras.at(0),
        input.framebuffer_viewport, app.drawgl_prms);
  };
  callbacks.widgets_cb = [&app](gui_window* win, const gui_input& input) {
    auto& params = app.drawgl_prms;
    draw_checkbox(win, "wireframe", params.wireframe);
    continue_line(win);
    draw_checkbox(win, "faceted", params.faceted);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_combobox(win, "lighting", (int&)params.lighting, shade_lighting_names);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
    // draw_label(win, "shape", app.name);
    // draw_label(win, "filename", app.filename);
    // draw_label(win, "outname", app.outname);
    // draw_label(win, "imagename", app.imagename);
    // auto& ioshape = app.ioshape;
    // draw_label(win, "points", std::to_string(ioshape.points.size()));
    // draw_label(win, "lines", std::to_string(ioshape.lines.size()));
    // draw_label(win, "triangles", std::to_string(ioshape.triangles.size()));
    // draw_label(win, "quads", std::to_string(ioshape.quads.size()));
    // draw_label(win, "positions", std::to_string(ioshape.positions.size()));
    // draw_label(win, "normals", std::to_string(ioshape.normals.size()));
    // draw_label(win, "texcoords", std::to_string(ioshape.texcoords.size()));
    // draw_label(win, "colors", std::to_string(ioshape.colors.size()));
    // draw_label(win, "radius", std::to_string(ioshape.radius.size()));
  };
  callbacks.update_cb = [](gui_window* win, const gui_input& input) {
    // update(win, apps);
  };
  callbacks.uiupdate_cb = [&app](gui_window* win, const gui_input& input) {
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
      auto& camera = app.ioscene.cameras.at(app.iocamera);
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, pan);
      set_frame(app.glscene.cameras.at(app.iocamera), camera.frame);
    }
  };

  // run ui
  run_ui({1280 + 320, 720}, "yshade", callbacks);

  // done
  return 0;
}
