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
#include <yocto/yocto_math.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_scene.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#if YOCTO_OPENGL == 1
#include <yocto_gui/yocto_glview.h>
#endif
using namespace yocto;

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

void add_command(cli_command& cli, const string& name, view_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
  add_option(cmd, "output", value.output, "Output shape.", {}, "o");
  add_option(cmd, "addsky", value.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view shapes
int run_view(const view_params& params) {
  // shape data
  auto shape = shape_data{};

  // load mesh
  auto ioerror = ""s;
  print_progress("load shape", 0, 1);
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
  }
  print_progress("load shape", 1, 1);

  // run view
  view_shape("yshape", params.shape, shape, params.addsky, print_progress);

  // done
  return 0;
}

#endif

struct glview_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glview_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_scene make_shapescene(
    const scene_shape& ioshape_, progress_callback progress_cb) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // handle progress
  auto progress = vec2i{0, 5};
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);

  // init scene
  auto scene = scene_scene{};

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

  // shapes
  if (progress_cb) progress_cb("create shape", progress.x++, progress.y);
  scene.shapes.emplace_back(ioshape);

  // instances
  if (progress_cb) progress_cb("create instance", progress.x++, progress.y);
  auto& shape_instance    = scene.instances.emplace_back();
  shape_instance.shape    = 0;
  shape_instance.material = 0;

  // done
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);
  return scene;
}

int run_glview(const glview_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto shape   = scene_shape{};
  print_progress("load shape", 0, 1);
  if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
  print_progress("load shape", 1, 1);

  // create scene
  auto scene = make_shapescene(shape, print_progress);

  // run viewer
  glview_scene(scene, params.shape, "", print_progress);

  // done
  return 0;
}

#endif

struct glpath_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glpath_params& value,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", value.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_scene make_pathscene(
    const scene_shape& ioshape_, progress_callback progress_cb) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // handle progress
  auto progress = vec2i{0, 5};
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);

  // init scene
  auto scene = scene_scene{};

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
  auto& shape_material      = scene.materials.emplace_back();
  shape_material.type       = material_type::plastic;
  shape_material.color      = {0.5, 1, 0.5};
  shape_material.roughness  = 0.2;
  auto& points_material     = scene.materials.emplace_back();
  points_material.type      = material_type::plastic;
  points_material.color     = {1, 0.5, 0.5};
  points_material.roughness = 0.2;
  auto& lines_material      = scene.materials.emplace_back();
  lines_material.type       = material_type::plastic;
  lines_material.color      = {0.5, 0.5, 1};
  lines_material.roughness  = 0.2;

  // shapes
  if (progress_cb) progress_cb("create shape", progress.x++, progress.y);
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(points_to_spheres({{0, 0, 0}}));
  scene.shapes.emplace_back(lines_to_cylinders({{0, 0, 0}, {0, 0, 0}}));

  // instances
  if (progress_cb) progress_cb("create instance", progress.x++, progress.y);
  auto& shape_instance     = scene.instances.emplace_back();
  shape_instance.shape     = 0;
  shape_instance.material  = 0;
  auto& points_instance    = scene.instances.emplace_back();
  points_instance.shape    = 1;
  points_instance.material = 1;
  auto& lines_instance     = scene.instances.emplace_back();
  lines_instance.shape     = 2;
  lines_instance.material  = 2;

  // done
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);
  return scene;
}

int run_glpath(const glpath_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  print_progress("load shape", 0, 1);
  if (!load_shape(params.shape, ioshape, ioerror, true)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads     = {};
  }
  print_progress("load shape", 1, 1);

  // create scene
  auto scene = make_pathscene(ioshape, print_progress);

  // bvh
  auto& shape = scene.shapes.at(0);
  auto  bvh   = make_triangles_bvh(shape.triangles, shape.positions, {});

  // stroke
  auto stroke = vector<shape_point>{};

  // geodesic solver
  auto adjacencies = face_adjacencies(shape.triangles);
  auto solver      = make_dual_geodesic_solver(
      shape.triangles, shape.positions, adjacencies);

  // run viewer
  glview_scene(
      scene, params.shape, "", print_progress,
      [&](gui_window* win, const gui_input& input, scene_scene& scene,
          shade_scene& glscene) {},
      [&](gui_window* win, const gui_input& input, scene_scene& scene,
          shade_scene& glscene) {
        auto& shape   = scene.shapes.at(0);
        auto& camera  = scene.cameras.at(0);
        auto  updated = false;
        if (input.mouse_left && input.modifier_ctrl) {
          if (input.modifier_shift) {
            stroke.clear();
            updated = true;
          } else {
            auto mouse_uv = vec2f{
                input.mouse_pos.x / float(input.window_size.x),
                input.mouse_pos.y / float(input.window_size.y)};
            auto ray  = camera_ray(camera.frame, camera.lens, camera.aspect,
                camera.film, mouse_uv);
            auto isec = intersect_triangles_bvh(
                bvh, shape.triangles, shape.positions, ray, false);
            if (isec.hit) {
              if (stroke.empty() || stroke.back().element != isec.element ||
                  stroke.back().uv != isec.uv) {
                stroke.push_back({isec.element, isec.uv});
                updated = true;
              }
            }
          }
        }
        if (updated) {
          auto positions = vector<vec3f>{};
          for (auto [element, uv] : stroke) {
            positions.push_back(eval_position(shape, element, uv));
          }
          auto& points = scene.shapes.at(1);
          points       = points_to_spheres(positions);
          set_positions(glscene.shapes.at(1), points.positions);
          set_normals(glscene.shapes.at(1), points.normals);
          set_texcoords(glscene.shapes.at(1), points.texcoords);
          set_quads(glscene.shapes.at(1), points.quads);
          glscene.shapes.at(1).point_size = 10;
          auto path       = compute_shortest_path(solver, shape.triangles,
              shape.positions, adjacencies, (vector<mesh_point>&)stroke);
          auto ppositions = vector<vec3f>{};
          for (auto [element, uv] : path) {
            ppositions.push_back(eval_position(shape, element, uv));
          }
          auto& lines = scene.shapes.at(2);
          lines       = lines_to_cylinders(ppositions);
          set_positions(glscene.shapes.at(2), lines.positions);
          set_normals(glscene.shapes.at(2), lines.normals);
          set_texcoords(glscene.shapes.at(2), lines.texcoords);
          set_quads(glscene.shapes.at(2), lines.quads);
        }
      });

  // done
  return 0;
}

#endif

struct app_params {
  string        command = "view";
  view_params   view    = {};
  glview_params glview  = {};
  glpath_params glpath  = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& value,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", value.command, "Command.");
  add_command(cli, "view", value.view, "View shapes.");
  add_command(cli, "glview", value.glview, "View shapes with OpenGL.");
  add_command(cli, "glpath", value.glpath, "Trace paths with OpenGL.");
}

// Parse cli
void parse_cli(app_params& params, int argc, const char** argv) {
  auto cli = cli_command{};
  add_commands(cli, "ymesh", params, "Process and view meshes.");
  parse_cli(cli, argc, argv);
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto params = app_params{};
  parse_cli(params, argc, argv);

  // dispatch commands
  if (params.command == "view") {
    return run_view(params.view);
  } else if (params.command == "glview") {
    return run_glview(params.glview);
  } else if (params.command == "glpath") {
    return run_glpath(params.glpath);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
