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

#include <queue>

// view params
struct view_params {
  string shape  = "shape.ply";
  string output = "out.ply";
  bool   addsky = false;
};

void add_command(cli_command& cli, const string& name, view_params& params,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "output", params.output, "Output shape.", {});
  add_option(cmd, "addsky", params.addsky, "Add sky.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_view(const view_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

// view shapes
int run_view(const view_params& params) {
  // load mesh
  auto shape   = scene_shape{};
  auto ioerror = ""s;
  if (path_filename(params.shape) == ".ypreset") {
    if (!make_shape_preset(shape, path_basename(params.shape), ioerror))
      print_fatal(ioerror);
  } else {
    if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);
  }

  // make scene
  auto scene = make_shape_scene(shape, params.addsky);

  // run view
  view_scene("ymesh", params.shape, scene);

  // done
  return 0;
}

#endif

struct glview_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glview_params& params,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glview(const glview_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_model make_shapescene(const scene_shape& ioshape_) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // init scene
  auto scene = scene_model{};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));
  // TODO(fabio): this should be a math function

  // camera
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  auto& shape_material     = scene.materials.emplace_back();
  shape_material.type      = scene_material_type::glossy;
  shape_material.color     = {0.5, 1, 0.5};
  shape_material.roughness = 0.2;

  // shapes
  scene.shapes.emplace_back(ioshape);

  // instances
  auto& shape_instance    = scene.instances.emplace_back();
  shape_instance.shape    = 0;
  shape_instance.material = 0;

  // done
  return scene;
}

int run_glview(const glview_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto shape   = scene_shape{};
  if (!load_shape(params.shape, shape, ioerror, true)) print_fatal(ioerror);

  // create scene
  auto scene = make_shapescene(shape);

  // run viewer
  glview_scene("ymesh", params.shape, scene, {});

  // done
  return 0;
}

#endif

struct glpath_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glpath_params& params,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glpath(const glpath_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_model make_pathscene(const scene_shape& ioshape_) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // init scene
  auto scene = scene_model{};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));
  // TODO(fabio): this should be a math function

  // camera
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  auto& shape_material      = scene.materials.emplace_back();
  shape_material.type       = scene_material_type::glossy;
  shape_material.color      = {0.5, 1, 0.5};
  shape_material.roughness  = 0.2;
  auto& points_material     = scene.materials.emplace_back();
  points_material.type      = scene_material_type::glossy;
  points_material.color     = {1, 0.5, 0.5};
  points_material.roughness = 0.2;
  auto& lines_material      = scene.materials.emplace_back();
  lines_material.type       = scene_material_type::glossy;
  lines_material.color      = {0.5, 0.5, 1};
  lines_material.roughness  = 0.2;

  // shapes
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(points_to_spheres({{0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));

  // instances
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
  return scene;
}

int run_glpath(const glpath_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  if (!load_shape(params.shape, ioshape, ioerror, true)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads     = {};
  }

  // create scene
  auto scene = make_pathscene(ioshape);

  // bvh
  auto& shape = scene.shapes.at(0);
  auto  bvh   = make_triangles_bvh(shape.triangles, shape.positions, {});

  // stroke
  auto stroke = vector<shape_point>{};

  // geodesic solver
  auto adjacencies = face_adjacencies(shape.triangles);
  auto solver      = make_dual_geodesic_solver(
      shape.triangles, shape.positions, adjacencies);
  auto bezier = true;

  // bezier algos
  auto params1      = spline_params{};
  params1.algorithm = spline_algorithm::de_casteljau_uniform;
  auto params2      = spline_params{};
  params2.algorithm = spline_algorithm::de_casteljau_adaptive;
  auto params3      = spline_params{};
  params3.algorithm = spline_algorithm::lane_riesenfeld_uniform;
  auto params4      = spline_params{};
  params4.algorithm = spline_algorithm::lane_riesenfeld_adaptive;

  // run viewer
  glview_scene(
      "ymesh", params.shape, scene, {},
      [&](const glinput_state&, vector<int>&, vector<int>&) {},
      [&](const glinput_state& input, vector<int>& updated_shapes,
          vector<int>&) {
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
          scene.shapes.at(1) = points_to_spheres(positions);
          updated_shapes.push_back(1);
          auto path = bezier ? compute_bezier_path(solver, shape.triangles,
                                   shape.positions, adjacencies,
                                   (vector<mesh_point>&)stroke, params1)
                             : compute_shortest_path(solver, shape.triangles,
                                   shape.positions, adjacencies,
                                   (vector<mesh_point>&)stroke);
          auto ppositions = vector<vec3f>{};
          for (auto [element, uv] : path) {
            ppositions.push_back(eval_position(shape, element, uv));
          }
          scene.shapes.at(2) = polyline_to_cylinders(ppositions);
          updated_shapes.push_back(2);
        }
      });

  // done
  return 0;
}

#endif

struct glpathd_params {
  string shape = "shape.ply";
};

// Cli
void add_command(cli_command& cli, const string& name, glpathd_params& params,
    const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
}

#ifndef YOCTO_OPENGL

// view shapes
int run_glpathd(const glpathd_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

static scene_model make_pathdscene(const scene_shape& ioshape) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // init scene
  auto scene = scene_model{};

  // camera
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o);

  // material
  auto& shape_material     = scene.materials.emplace_back();
  shape_material.type      = scene_material_type::glossy;
  shape_material.color     = {0.5, 1, 0.5};
  shape_material.roughness = 0.2;
  auto& points_material    = scene.materials.emplace_back();
  points_material.type     = scene_material_type::matte;
  points_material.color    = {1, 0.5, 0.5};
  auto& lines1_material    = scene.materials.emplace_back();
  lines1_material.type     = scene_material_type::matte;
  lines1_material.color    = {0.5, 0.5, 1};
  auto& lines2_material    = scene.materials.emplace_back();
  lines2_material.type     = scene_material_type::matte;
  lines2_material.color    = {1, 1, 0.5};
  auto& lines3_material    = scene.materials.emplace_back();
  lines3_material.type     = scene_material_type::matte;
  lines3_material.color    = {1, 0.5, 1};
  auto& lines4_material    = scene.materials.emplace_back();
  lines4_material.type     = scene_material_type::matte;
  lines4_material.color    = {0.5, 0.5, 0.5};
  auto& edges_material     = scene.materials.emplace_back();
  edges_material.type      = scene_material_type::matte;
  edges_material.color     = {0, 0, 0};

  // shapes
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(points_to_spheres({{0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));
  scene.shapes.emplace_back(polyline_to_cylinders({{0, 0, 0}, {0, 0, 0}}));

// make edges
#if 0
  auto edges         = get_edges(ioshape.triangles);
  auto edge_vertices = vector<vec3f>{};
  for (auto [i, j] : edges) {
    edge_vertices.push_back(ioshape.positions[i]);
    edge_vertices.push_back(ioshape.positions[j]);
  }
  scene.shapes.back() = lines_to_cylinders(edge_vertices, 4, 0.001);
#endif

  // instances
  auto& shape_instance     = scene.instances.emplace_back();
  shape_instance.shape     = 0;
  shape_instance.material  = 0;
  auto& points_instance    = scene.instances.emplace_back();
  points_instance.shape    = 1;
  points_instance.material = 1;
  auto& lines1_instance    = scene.instances.emplace_back();
  lines1_instance.shape    = 2;
  lines1_instance.material = 2;
  auto& lines2_instance    = scene.instances.emplace_back();
  lines2_instance.shape    = 3;
  lines2_instance.material = 3;
  auto& lines3_instance    = scene.instances.emplace_back();
  lines3_instance.shape    = 4;
  lines3_instance.material = 4;
  auto& lines4_instance    = scene.instances.emplace_back();
  lines4_instance.shape    = 5;
  lines4_instance.material = 5;
  auto& edges_instance     = scene.instances.emplace_back();
  edges_instance.shape     = 6;
  edges_instance.material  = 6;

  // done
  return scene;
}

int run_glpathd(const glpathd_params& params) {
  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  if (!load_shape(params.shape, ioshape, ioerror, true)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads     = {};
  }

  // rescale shape to unit
  auto bbox = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));

  // reordina indici nei triangoli
  auto rng = rng_state{};
  for (auto& t : ioshape.triangles) {
    auto o = rand1i(rng, 3);
    t      = {t[(0 + o) % 3], t[(1 + o) % 3], t[(2 + o) % 3]};
  }

  // create scene
  auto scene = make_pathdscene(ioshape);

  // bvh
  auto& shape = scene.shapes.at(0);
  auto  bvh   = make_triangles_bvh(shape.triangles, shape.positions, {});

  // geodesic solver
  auto adjacencies = face_adjacencies(shape.triangles);
  auto solver      = make_dual_geodesic_solver(
      shape.triangles, shape.positions, adjacencies);

  // other solver
  //  auto v2t = vertex_to_triangles(shape.triangles, shape.positions,
  //  adjacencies); auto solver2 = make_geodesic_solver(
  //      shape.triangles, shape.positions, adjacencies, v2t);
  //  auto total_angles = vector<float>{};
  //  auto angles       = compute_angles(
  //      shape.triangles, shape.positions, adjacencies, v2t, total_angles,
  //      true);

  // points at random
  auto point1 = mesh_point{0, {0.5, 0.5}};
  auto point2 = mesh_point{1, {0.5, 0.5}};

  // run viewer
  glview_scene(
      "ymesh", params.shape, scene, {},
      [&](const glinput_state&, vector<int>&, vector<int>&) {},
      [&](const glinput_state& input, vector<int>& updated_shapes,
          vector<int>&) {
        auto& shape   = scene.shapes.at(0);
        auto& camera  = scene.cameras.at(0);
        auto  updated = false;
        if (input.mouse_left && input.modifier_ctrl) {
          auto mouse_uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
              input.mouse_pos.y / float(input.window_size.y)};
          auto ray      = camera_ray(
              camera.frame, camera.lens, camera.aspect, camera.film, mouse_uv);
          auto isec = intersect_triangles_bvh(
              bvh, shape.triangles, shape.positions, ray, false);
          if (isec.hit) {
            if (input.modifier_shift) {
              point2 = {isec.element, isec.uv};
            } else {
              point1 = {isec.element, isec.uv};
            }
            updated = true;
          }
        }
        if (updated) {
          auto positions = vector<vec3f>{};
          positions.push_back(eval_position(shape, point1.face, point1.uv));
          positions.push_back(eval_position(shape, point2.face, point2.uv));
          scene.shapes.at(1) = points_to_spheres(positions, 2, 0.002);
          updated_shapes.push_back(1);
          auto path1      = compute_shortest_path(solver, shape.triangles,
              shape.positions, adjacencies, point1, point2);
          auto positions1 = vector<vec3f>{};
          for (auto [element, uv] : path1) {
            positions1.push_back(eval_position(shape, element, uv));
          }
          scene.shapes.at(2) = polyline_to_cylinders(positions1, 4, 0.002);
          updated_shapes.push_back(2);
          auto positions2    = visualize_shortest_path(solver, shape.triangles,
              shape.positions, adjacencies, point1, point2, true);
          scene.shapes.at(3) = polyline_to_cylinders(positions2, 4, 0.002);
          updated_shapes.push_back(3);
          // auto path3 = visualize_shortest_path(solver2, shape.triangles,
          //     shape.positions, adjacencies, v2t, angles, point1, point2,
          //     false);
          // auto positions3 = vector<vec3f>{};
          // for (auto [element, uv] : path3) {
          //   positions3.push_back(eval_position(shape, element, uv));
          // }
          // auto &lines3 = scene.shapes.at(4);
          // lines3       = polyline_to_cylinders(positions3, 4, 0.002);
          // set_positions(glscene.shapes.at(4), lines3.positions);
          // set_normals(glscene.shapes.at(4), lines3.normals);
          // set_texcoords(glscene.shapes.at(4), lines3.texcoords);
          // set_quads(glscene.shapes.at(4), lines3.quads);
          // auto path4      = visualize_shortest_path(solver2, shape.triangles,
          //     shape.positions, adjacencies, v2t, angles, point1, point2,
          //     true);
          // auto positions4 = vector<vec3f>{};
          // for (auto [element, uv] : path2) {
          //   positions4.push_back(eval_position(shape, element, uv));
          // }
          // auto &lines4 = scene.shapes.at(5);
          // lines4       = polyline_to_cylinders(positions4, 4, 0.002);
          // set_positions(glscene.shapes.at(5), lines4.positions);
          // set_normals(glscene.shapes.at(5), lines4.normals);
          // set_texcoords(glscene.shapes.at(5), lines4.texcoords);
          // set_quads(glscene.shapes.at(5), lines4.quads);
        }
      });

  // done
  return 0;
}

#endif

struct glsculpt_params {
  string shape   = "shape.ply"s;
  string texture = "";
};

// Cli
inline void add_command(cli_command& cli, const string& name,
    glsculpt_params& params, const string& usage) {
  auto& cmd = add_command(cli, name, usage);
  add_argument(cmd, "shape", params.shape, "Input shape.");
  add_option(cmd, "texture", params.texture, "Brush texture.");
}

#ifndef YOCTO_OPENGL

// view scene
int run_glsculpt(const glsculpt_params& params) {
  return print_fatal("Opengl not compiled");
}

#else

enum struct sculpt_brush_type { gaussian, texture, smooth };
auto const sculpt_brush_names = vector<std::string>{
    "gaussian brush", "texture brush", "smooth brush"};

struct sculpt_params {
  sculpt_brush_type type     = sculpt_brush_type::gaussian;
  float             radius   = 0.35f;
  float             strength = 1.0f;
  bool              negative = false;
};

struct sculpt_state {
  // data structures
  shape_bvh           bvh         = {};
  hash_grid           grid        = {};
  vector<vector<int>> adjacencies = {};
  geodesic_solver     solver      = {};
  // brush
  scene_texture tex_image = {};
  // stroke
  vector<shape_point> stroke   = {};
  vec2f               last_uv  = {};
  bool                instroke = false;
  // shape at the beginning of the stroke
  scene_shape base_shape = {};
};

// Initialize all sculpting parameters.
sculpt_state make_sculpt_state(
    const scene_shape& shape, const scene_texture& texture) {
  auto state = sculpt_state{};
  state.bvh  = make_triangles_bvh(
      shape.triangles, shape.positions, shape.radius);
  state.grid       = make_hash_grid(shape.positions, 0.05f);
  auto adjacencies = face_adjacencies(shape.triangles);
  state.solver     = make_geodesic_solver(
      shape.triangles, adjacencies, shape.positions);
  state.adjacencies = vertex_adjacencies(shape.triangles, adjacencies);
  state.tex_image   = texture;
  state.base_shape  = shape;
  state.base_shape.texcoords.assign(shape.positions.size(), {0, 0});
  return state;
}

scene_shape make_circle(
    const vec3f& center, const mat3f& basis, float radius, int steps) {
  // 4 initial vertices
  auto  lines    = make_lines({1, 4});
  vec3f next_dir = basis.x;
  for (auto line : lines.lines) {
    auto v1             = line.x;
    auto v2             = line.y;
    lines.positions[v1] = center + next_dir * radius;
    lines.normals[v1]   = basis.z;
    next_dir            = cross(basis.z, next_dir);
    lines.positions[v2] = center + next_dir * radius;
    lines.normals[v2]   = basis.z;
  }

  // create polylines and fix lenght to radius
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto lins      = vector<vec2i>{};
  steps          = int(steps / 4);
  for (auto line : lines.lines) {
    auto v1  = line.x;
    auto v2  = line.y;
    auto pos = lines.positions[v1];
    positions.push_back(pos);
    normals.push_back(lines.normals[v1]);
    auto dir    = lines.positions[v2] - lines.positions[v1];
    auto lenght = length(dir) / steps;
    dir         = normalize(dir);
    for (int i = 0; i < steps; i++) {
      auto new_pos = pos + dir * lenght * (i + 1);
      auto new_dir = normalize(new_pos - center);
      new_pos      = center + new_dir * radius;
      positions.push_back(new_pos);
      normals.push_back(basis.z);
      lins.push_back({int(positions.size() - 2), int(positions.size() - 1)});
    }
  }

  // apply
  lines.positions = positions;
  lines.normals   = normals;
  lines.lines     = lins;
  return lines;
}

// To visualize mouse intersection on mesh
scene_shape make_cursor(const vec3f& position, const vec3f& normal,
    float radius, float height = 0.05f) {
  auto basis  = basis_fromz(normal);
  auto cursor = make_circle(position, basis, radius, 32);
  cursor.normals.clear();
  cursor.texcoords.clear();
  cursor.positions.push_back(position);
  cursor.positions.push_back(position + normal * 0.05f);
  cursor.lines.push_back(
      {int(cursor.positions.size() - 2), int(cursor.positions.size() - 1)});
  return cursor;
}

// To take shape positions indices associate with planar coordinates
vector<int> stroke_parameterization(vector<vec2f>& coords,
    const geodesic_solver& solver, const vector<int>& sampling,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    float radius) {
  if (normals.empty()) return vector<int>{};

  // Planar coordinates by local computation between neighbors
  auto compute_coordinates = [](vector<vec2f>&        coords,
                                 const vector<mat3f>& frames,
                                 const vector<vec3f>& positions, int node,
                                 int neighbor, float weight) -> void {
    // Project a vector on a plane, maintaining vector length
    auto project_onto_plane = [](const mat3f& basis, const vec3f& p) -> vec2f {
      auto v  = p - dot(p, basis.z) * basis.z;
      auto v1 = vec2f{dot(v, basis.x), dot(v, basis.y)};
      return v1 * (length(p) / length(v1));
    };

    auto current_coord = coords[node];
    auto edge          = positions[node] - positions[neighbor];
    auto projection    = project_onto_plane(frames[neighbor], edge);
    auto new_coord     = coords[neighbor] + projection;
    auto avg_lenght    = (length(current_coord) + length(new_coord)) / 2;
    auto new_dir       = normalize(current_coord + new_coord);
    coords[node] = current_coord == zero2f ? new_coord : new_dir * avg_lenght;

    // following doesn't work
    // coords[node] = current_coord + (weight * (coords[neighbor] +
    // projection));
  };

  // Frame by local computation between neighbors
  auto compute_frame = [](vector<mat3f>& frames, const vector<vec3f>& normals,
                           int node, int neighbor, float weight) -> void {
    auto current_dir = frames[node].x;
    auto rotation    = basis_fromz(normals[neighbor]) *
                    transpose(basis_fromz(normals[node]));
    auto neighbor_dir = frames[neighbor].x;
    current_dir       = current_dir + (rotation * weight * neighbor_dir);
    frames[node].z    = normals[node];
    frames[node].y    = cross(frames[node].z, normalize(current_dir));
    frames[node].x    = cross(frames[node].y, frames[node].z);
  };

  // Classic Dijkstra
  auto dijkstra = [](const geodesic_solver& solver, const vector<int>& sources,
                      vector<float>& distances, float max_distance,
                      auto&& update) -> void {
    auto compare = [&](int i, int j) { return distances[i] > distances[j]; };
    std::priority_queue<int, vector<int>, decltype(compare)> queue(compare);

    // setup queue
    for (auto source : sources) queue.push(source);

    while (!queue.empty()) {
      int node = queue.top();
      queue.pop();

      auto distance = distances[node];
      if (distance > max_distance) continue;  // early exit

      for (auto arc : solver.graph[node]) {
        auto new_distance = distance + arc.length;

        update(node, arc.node, new_distance);

        if (new_distance < distances[arc.node]) {
          distances[arc.node] = new_distance;
          queue.push(arc.node);
        }
      }
    }
  };

  // Compute initial frames of stroke sampling vertices
  auto compute_stroke_frames = [](vector<mat3f>&        frames,
                                   const vector<vec3f>& positions,
                                   const vector<vec3f>& normals,
                                   const vector<int>& stroke_sampling) -> void {
    // frames follow stroke direction
    for (int i = 0; i < (int)stroke_sampling.size() - 1; i++) {
      int  curr    = stroke_sampling[i];
      int  next    = stroke_sampling[i + 1];
      auto dir     = positions[next] - positions[curr];
      auto z       = normals[curr];
      auto y       = cross(z, normalize(dir));
      auto x       = cross(y, z);
      frames[curr] = {x, y, z};
    }

    if (stroke_sampling.size() == 1) {
      frames[stroke_sampling[0]] = basis_fromz(normals[stroke_sampling[0]]);
    } else {
      int  final    = stroke_sampling[stroke_sampling.size() - 1];
      int  prev     = stroke_sampling[stroke_sampling.size() - 2];
      auto dir      = positions[final] - positions[prev];
      auto z        = normals[final];
      auto y        = cross(z, normalize(dir));
      auto x        = cross(y, z);
      frames[final] = {x, y, z};
    }

    // average frames direction of middle stroke vertices
    for (int i = 1; i < stroke_sampling.size() - 1; i++) {
      int  curr    = stroke_sampling[i];
      int  next    = stroke_sampling[i + 1];
      int  prev    = stroke_sampling[i - 1];
      auto dir     = frames[prev].x + frames[next].x;
      auto z       = normals[curr];
      auto y       = cross(z, normalize(dir));
      auto x       = cross(y, z);
      frames[curr] = {x, y, z};
    }
  };

  // init params
  auto vertices = std::unordered_set<int>{};  // to avoid duplicates

  auto visited = vector<bool>(positions.size(), false);
  for (auto sample : sampling) visited[sample] = true;

  coords              = vector<vec2f>(solver.graph.size(), zero2f);
  coords[sampling[0]] = {radius, radius};
  vertices.insert(sampling[0]);
  for (int i = 1; i < sampling.size(); i++) {
    auto edge           = positions[sampling[i]] - positions[sampling[i - 1]];
    coords[sampling[i]] = {coords[sampling[i - 1]].x + length(edge), radius};
    vertices.insert(sampling[i]);
  }

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : sampling) distances[sample] = 0.0f;

  auto frames = vector<mat3f>(positions.size(), identity3x3f);
  compute_stroke_frames(frames, positions, normals, sampling);

  auto update = [&](int node, int neighbor, float new_distance) {
    vertices.insert(node);
    if (!visited[neighbor]) return;
    float weight = length(positions[neighbor] - positions[node]) + flt_eps;
    weight       = 1.0f / weight;
    compute_coordinates(coords, frames, positions, node, neighbor, weight);
    compute_frame(frames, normals, node, neighbor, weight);
    visited[node] = true;
  };
  dijkstra(solver, sampling, distances, radius, update);

  auto vec_vertices = vector<int>(vertices.begin(), vertices.end());

  // conversion in [0, 1]
  for (int i = 0; i < vertices.size(); i++)
    coords[vec_vertices[i]] /= radius * 2.0f;

  return vec_vertices;
}

// Compute gaussian function
float gaussian_distribution(const vec3f& origin, const vec3f& position,
    float standard_dev, float scale_factor, float strength, float radius) {
  auto scaled_strength = strength / ((((radius - 0.1f) * 0.5f) / 0.7f) + 0.2f);
  auto N               = 1.0f / (((standard_dev * scaled_strength) *
                       (standard_dev * scaled_strength) *
                       (standard_dev * scaled_strength)) *
                      yocto::sqrt((2.0f * pif) * (2.0f * pif) * (2.0f * pif)));
  auto d               = (origin - position) * scale_factor;
  auto E               = dot(d, d) / (2.0f * standard_dev * standard_dev);
  return N * yocto::exp(-E);
}

// To apply brush on intersected points' neighbors
bool gaussian_brush(vector<vec3f>& positions, const hash_grid& grid,
    const vector<vec3i>& triangles, const vector<vec3f>& base_positions,
    const vector<vec3f>& base_normals, const vector<shape_point>& stroke,
    const sculpt_params& params) {
  if (stroke.empty()) return false;

  // helpers
  auto eval_position = [](const vector<vec3i>&  triangles,
                           const vector<vec3f>& positions, int element,
                           const vec2f& uv) {
    auto& triangle = triangles[element];
    return interpolate_triangle(positions[triangle.x], positions[triangle.y],
        positions[triangle.z], uv);
  };
  auto eval_normal = [](const vector<vec3i>&  triangles,
                         const vector<vec3f>& normals, int element,
                         const vec2f& uv) {
    auto& triangle = triangles[element];
    return normalize(interpolate_triangle(
        normals[triangle.x], normals[triangle.y], normals[triangle.z], uv));
  };

  // for a correct gaussian distribution
  float scale_factor = 3.5f / params.radius;

  auto neighbors = vector<int>{};
  for (auto [element, uv] : stroke) {
    auto position = eval_position(triangles, base_positions, element, uv);
    auto normal   = eval_normal(triangles, base_normals, element, uv);
    find_neighbors(grid, neighbors, position, params.radius);
    if (params.negative) normal = -normal;
    for (auto neighbor : neighbors) {
      auto gauss_height = gaussian_distribution(position, positions[neighbor],
          0.7f, scale_factor, params.strength, params.radius);
      positions[neighbor] += normal * gauss_height;
    }
    neighbors.clear();
  }

  return true;
}

// Compute texture values through the parameterization
bool texture_brush(vector<vec3f>& positions, vector<vec2f>& texcoords,
    const geodesic_solver& solver, const scene_texture& texture,
    const vector<vec3i>& triangles, const vector<vec3f>& base_positions,
    const vector<vec3f>& base_normals, const vector<shape_point>& stroke,
    const sculpt_params& params) {
  if (texture.pixelsf.empty() && texture.pixelsb.empty()) return false;

  // Taking closest vertex of an intersection
  auto closest_vertex = [](const vector<vec3i>& triangles, int element,
                            const vec2f& uv) -> int {
    auto& triangle = triangles[element];
    if (uv.x < 0.5f && uv.y < 0.5f) return triangle.x;
    if (uv.x > uv.y) return triangle.y;
    return triangle.z;
  };

  auto sampling = vector<int>{};
  for (auto [element, uv] : stroke) {
    sampling.push_back(closest_vertex(triangles, element, uv));
  }

  auto vertices = stroke_parameterization(
      texcoords, solver, sampling, base_positions, base_normals, params.radius);
  if (vertices.empty()) return false;

  positions = base_positions;

  auto scale_factor = 3.5f / params.radius;
  auto max_height   = gaussian_distribution(
      zero3f, zero3f, 0.7f, scale_factor, params.strength, params.radius);

  for (auto idx : vertices) {
    auto uv     = texcoords[idx];
    auto height = max(xyz(eval_texture(texture, uv)));
    auto normal = base_normals[idx];
    if (params.negative) normal = -normal;
    height *= max_height;
    positions[idx] += normal * height;
  }

  return true;
}

// Cotangent operator

// Compute edge cotangents weights
float laplacian_weight(const vector<vec3f>& positions,
    const vector<vector<int>>& adjacencies, int node, int neighbor) {
  auto cotan = [](vec3f& a, vec3f& b) {
    return dot(a, b) / length(cross(a, b));
  };

  auto num_neighbors = int(adjacencies[node].size());

  int ind = -1;
  for (int i = 0; i < adjacencies[node].size(); i++) {
    if (adjacencies[node][i] == neighbor) {
      ind = i;
      break;
    }
  }

  auto prev = adjacencies[node][(ind - 1) >= 0 ? ind - 1 : num_neighbors - 1];
  auto next = adjacencies[node][(ind + 1) % num_neighbors];

  auto v1 = positions[node] - positions[prev];
  auto v2 = positions[neighbor] - positions[prev];
  auto v3 = positions[node] - positions[next];
  auto v4 = positions[neighbor] - positions[next];

  float cot_alpha = cotan(v1, v2);
  float cot_beta  = cotan(v3, v4);
  float weight    = (cot_alpha + cot_beta) / 2;

  float cotan_max = yocto::cos(flt_min) / yocto::sin(flt_min);
  weight          = clamp(weight, cotan_max, -cotan_max);
  return weight;
}

// Smooth brush with Laplace Operator Discretization and Cotangents Weights
bool smooth_brush(vector<vec3f>& positions, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vector<int>>& adjacencies,
    vector<shape_point>& stroke, const sculpt_params& params) {
  if (stroke.empty()) return false;

  // Taking closest vertex of an intersection
  auto closest_vertex = [](const vector<vec3i>& triangles, int element,
                            const vec2f& uv) -> int {
    auto& triangle = triangles[element];
    if (uv.x < 0.5f && uv.y < 0.5f) return triangle.x;
    if (uv.x > uv.y) return triangle.y;
    return triangle.z;
  };

  auto stroke_vertices = vector<int>{};
  for (auto [element, uv] : stroke)
    stroke_vertices.push_back(closest_vertex(triangles, element, uv));

  // Classic Dijkstra
  auto dijkstra = [](const geodesic_solver& solver, const vector<int>& sources,
                      vector<float>& distances, float max_distance,
                      auto&& update) -> void {
    auto compare = [&](int i, int j) { return distances[i] > distances[j]; };
    std::priority_queue<int, vector<int>, decltype(compare)> queue(compare);

    // setup queue
    for (auto source : sources) queue.push(source);

    while (!queue.empty()) {
      int node = queue.top();
      queue.pop();

      auto distance = distances[node];
      if (distance > max_distance) continue;  // early exit

      for (auto arc : solver.graph[node]) {
        auto new_distance = distance + arc.length;

        update(node, arc.node, new_distance);

        if (new_distance < distances[arc.node]) {
          distances[arc.node] = new_distance;
          queue.push(arc.node);
        }
      }
    }
  };

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : stroke_vertices) distances[sample] = 0.0f;

  auto current_node = -1;
  auto neighbors    = vector<int>{};
  auto weights      = vector<float>{};
  auto update       = [&](int node, int neighbor, float new_distance) {
    if (current_node == -1) current_node = node;
    if (node != current_node) {
      vec3f sum1 = zero3f;
      float sum2 = 0.0f;
      for (int i = 0; i < neighbors.size(); i++) {
        sum1 += positions[neighbors[i]] * weights[i];
        sum2 += weights[i];
      }
      positions[current_node] += 0.5f *
                                 ((sum1 / sum2) - positions[current_node]);
      current_node = node;
      neighbors.clear();
      weights.clear();
    }
    neighbors.push_back(neighbor);
    weights.push_back(laplacian_weight(positions, adjacencies, node, neighbor));
  };
  dijkstra(solver, stroke_vertices, distances, params.radius, update);

  return true;
}

static scene_model make_sculptscene(const scene_shape& ioshape_) {
  // Frame camera
  auto camera_frame = [](float lens, float aspect,
                          float film = 0.036) -> frame3f {
    auto camera_dir  = normalize(vec3f{0, 0.5, 1});
    auto bbox_radius = 2.0f;
    auto camera_dist = bbox_radius * lens / (film / aspect);
    return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
  };

  // init scene
  auto scene = scene_model{};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto& pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto& pos : ioshape.positions) pos -= center(bbox);
  for (auto& pos : ioshape.positions) pos /= max(size(bbox));

  // camera
  auto& camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  auto& shape_material  = scene.materials.emplace_back();
  shape_material.type   = scene_material_type::matte;
  shape_material.color  = {0.78f, 0.31f, 0.23f};
  auto& cursor_material = scene.materials.emplace_back();
  cursor_material.type  = scene_material_type::matte;

  // shapes
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(make_cursor({0, 0, 0}, {0, 0, 1}, 1));

  // instances
  auto& shape_instance     = scene.instances.emplace_back();
  shape_instance.shape     = 0;
  shape_instance.material  = 0;
  auto& cursor_instance    = scene.instances.emplace_back();
  cursor_instance.shape    = 1;
  cursor_instance.material = 1;

  // done
  return scene;
}

// To make the stroke sampling (position, normal) following the mouse
static pair<vector<shape_point>, vec2f> sample_stroke(const shape_bvh& bvh,
    const scene_shape& shape, const vec2f& last_uv, const vec2f& mouse_uv,
    const scene_camera& camera, const sculpt_params& params) {
  // helper
  auto intersect_shape = [&](const vec2f& uv) {
    auto ray = camera_ray(
        camera.frame, camera.lens, camera.aspect, camera.film, uv);
    return intersect_triangles_bvh(
        bvh, shape.triangles, shape.positions, ray, false);
  };

  // eval current intersection
  auto last  = intersect_shape(last_uv);
  auto mouse = intersect_shape(mouse_uv);
  if (!mouse.hit || !last.hit) return {{}, last_uv};

  // sample
  auto delta_pos   = distance(eval_position(shape, last.element, last.uv),
      eval_position(shape, mouse.element, mouse.uv));
  auto stroke_dist = params.radius * 0.2f;
  auto steps       = int(delta_pos / stroke_dist);
  if (steps == 0) return {};
  auto update_uv = (mouse_uv - last_uv) * stroke_dist / delta_pos;
  auto cur_uv    = last_uv;
  auto samples   = vector<shape_point>{};
  for (auto step = 0; step < steps; step++) {
    cur_uv += update_uv;
    auto isec = intersect_shape(cur_uv);
    if (!isec.hit) continue;
    samples.push_back({isec.element, isec.uv});
  }

  return {samples, cur_uv};
}

static pair<bool, bool> sculpt_update(sculpt_state& state, scene_shape& shape,
    scene_shape& cursor, const scene_camera& camera, const vec2f& mouse_uv,
    bool mouse_pressed, const sculpt_params& params) {
  auto updated_shape = false, updated_cursor = false;

  auto ray = camera_ray(
      camera.frame, camera.lens, camera.aspect, camera.film, mouse_uv);
  auto isec = intersect_triangles_bvh(
      state.bvh, shape.triangles, shape.positions, ray, false);
  if (isec.hit) {
    cursor         = make_cursor(eval_position(shape, isec.element, isec.uv),
        eval_normal(shape, isec.element, isec.uv),
        params.radius *
            (params.type == sculpt_brush_type::gaussian ? 0.5f : 1.0f));
    updated_cursor = true;
  }

  // sculpting
  if (isec.hit && mouse_pressed) {
    if (!state.instroke) {
      state.instroke = true;
      state.last_uv  = mouse_uv;
    } else {
      auto last_uv           = state.last_uv;
      auto [samples, cur_uv] = sample_stroke(
          state.bvh, shape, last_uv, mouse_uv, camera, params);
      if (!samples.empty()) {
        state.last_uv = cur_uv;
        if (params.type == sculpt_brush_type::gaussian) {
          state.stroke  = samples;
          updated_shape = gaussian_brush(shape.positions, state.grid,
              shape.triangles, state.base_shape.positions,
              state.base_shape.normals, state.stroke, params);
        } else if (params.type == sculpt_brush_type::smooth) {
          state.stroke  = samples;
          updated_shape = smooth_brush(shape.positions, state.solver,
              shape.triangles, state.adjacencies, state.stroke, params);
        } else if (params.type == sculpt_brush_type::texture) {
          state.stroke.insert(
              state.stroke.end(), samples.begin(), samples.end());
          updated_shape = texture_brush(shape.positions,
              state.base_shape.texcoords, state.solver, state.tex_image,
              state.base_shape.triangles, state.base_shape.positions,
              state.base_shape.normals, state.stroke, params);
        }
      }
    }
    if (updated_shape) {
      triangles_normals(shape.normals, shape.triangles, shape.positions);
      update_triangles_bvh(state.bvh, shape.triangles, shape.positions);
      state.grid = make_hash_grid(shape.positions, state.grid.cell_size);
    }
  } else {
    if (state.instroke) {
      state.instroke = false;
      state.stroke.clear();
      state.base_shape.positions = shape.positions;
      state.base_shape.normals   = shape.normals;
      state.base_shape.texcoords.assign(
          state.base_shape.texcoords.size(), {0, 0});
    }
  }

  return {updated_shape, updated_cursor};
}

int run_glsculpt(const glsculpt_params& params_) {
  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  if (!load_shape(params_.shape, ioshape, ioerror, true)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads.clear();
  }

  // loading texture
  auto texture = scene_texture{};
  if (!params_.texture.empty()) {
    if (!load_texture(params_.texture, texture, ioerror)) print_fatal(ioerror);
  }

  // setup app
  auto scene = make_sculptscene(ioshape);

  // sculpt params
  auto params = sculpt_params{};
  auto state  = make_sculpt_state(scene.shapes.front(), texture);

  // callbacks
  glview_scene(
      "ymesh", params_.shape, scene, {},
      [&](const glinput_state&, vector<int>&, vector<int>&) {
        draw_glcombobox("brush type", (int&)params.type, sculpt_brush_names);
        if (params.type == sculpt_brush_type::gaussian) {
          if (params.strength < 0.8f || params.strength > 1.5f)
            params.strength = 1.0f;
          draw_glslider("radius", params.radius, 0.1f, 0.8f);
          draw_glslider("strength", params.strength, 1.5f, 0.9f);
          draw_glcheckbox("negative", params.negative);
        } else if (params.type == sculpt_brush_type::texture) {
          if (params.strength < 0.8f || params.strength > 1.5f)
            params.strength = 1.0f;
          draw_glslider("radius", params.radius, 0.1f, 0.8f);
          draw_glslider("strength", params.strength, 1.5f, 0.9f);
          draw_glcheckbox("negative", params.negative);
        } else if (params.type == sculpt_brush_type::smooth) {
          draw_glslider("radius", params.radius, 0.1f, 0.8f);
          draw_glslider("strength", params.strength, 0.1f, 1.0f);
        }
      },
      [&](const glinput_state& input, vector<int>& updated_shapes,
          vector<int>&) {
        auto  mouse_uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
            input.mouse_pos.y / float(input.window_size.y)};
        auto& shape    = scene.shapes.at(0);
        auto& cursor   = scene.shapes.at(1);
        auto& camera   = scene.cameras.at(0);
        auto [updated_shape, updated_cursor] = sculpt_update(state, shape,
            cursor, camera, mouse_uv, input.mouse_left && input.modifier_ctrl,
            params);
        if (updated_cursor) {
          updated_shapes.push_back(1);
          // glscene.instances.at(1).hidden = false;
        } else {
          // glscene.instances.at(1).hidden = true;
        }
        if (updated_shape) {
          updated_shapes.push_back(0);
        }
      });

  // done
  return 0;
}

#endif

struct app_params {
  string          command  = "view";
  view_params     view     = {};
  glview_params   glview   = {};
  glpath_params   glpath   = {};
  glpathd_params  glpathd  = {};
  glsculpt_params glsculpt = {};
};

// Cli
void add_commands(cli_command& cli, const string& name, app_params& params,
    const string& usage) {
  cli = make_cli(name, usage);
  add_command_name(cli, "command", params.command, "Command.");
  add_command(cli, "view", params.view, "View shapes.");
  add_command(cli, "glview", params.glview, "View shapes with OpenGL.");
  add_command(cli, "glpath", params.glpath, "Trace paths with OpenGL.");
  add_command(cli, "glpathd", params.glpathd, "Trace debug paths with OpenGL.");
  add_command(cli, "glsculpt", params.glsculpt, "Sculpt meshes with OpenGL.");
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
  } else if (params.command == "glpathd") {
    return run_glpathd(params.glpathd);
  } else if (params.command == "glsculpt") {
    return run_glsculpt(params.glsculpt);
  } else {
    return print_fatal("unknown command " + params.command);
  }
}
