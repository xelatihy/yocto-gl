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
#include <yocto/yocto_math.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_sampling.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
using namespace yocto;

#include "ext/json.hpp"

using json = nlohmann::json;

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

using json = nlohmann::json;
using std::array;

// support for json conversions
inline void to_json(json& j, const vec2f& value) {
  nlohmann::to_json(j, (const array<float, 2>&)value);
}
inline void to_json(json& j, const vec3f& value) {
  nlohmann::to_json(j, (const array<float, 3>&)value);
}
inline void to_json(json& j, const vec4f& value) {
  nlohmann::to_json(j, (const array<float, 4>&)value);
}
inline void to_json(json& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}
inline void to_json(json& j, const mat4f& value) {
  nlohmann::to_json(j, (const array<float, 16>&)value);
}

inline void from_json(const json& j, vec3f& value) {
  nlohmann::from_json(j, (array<float, 3>&)value);
}
inline void from_json(const json& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const json& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}

inline void to_json(json& j, const mesh_point& value) {
  nlohmann::to_json(j, pair{value.face, value.uv});
}

}  // namespace yocto

// load/save json
static bool load_json(const string& filename, json& js, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    js = json::parse(text);
    return true;
  } catch (std::exception&) {
    return parse_error();
  }
}

static bool save_json(const string& filename, const json& js, string& error) {
  return save_text(filename, js.dump(2), error);
}

vector<mesh_point> sample_points(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const shape_bvh& bvh, int num_points = 4) {
  // pick based on area
  auto points = vector<mesh_point>{};
  auto rng    = make_rng(9867198237913);
  auto cdf    = sample_triangles_cdf(triangles, positions);
  for (auto idx = 0; idx < num_points; idx++) {
    auto [triangle, uv] = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    points.push_back({mesh_point{triangle, uv}});
  }
  return points;
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto smooth    = false;
  auto faceted   = false;
  auto validate  = false;
  auto info      = false;
  auto pathname  = "path.ply"s;
  auto statsname = "stats.json"s;
  auto scenename = "scene.json"s;
  auto meshname  = "mesh.ply"s;

  // parse command line
  auto cli = make_cli("ymshproc", "Applies operations on a triangle mesh");
  add_option(cli, "--smooth", smooth, "Compute smooth normals");
  add_option(cli, "--faceted", faceted, "Remove normals");
  add_option(cli, "--validate,-v", validate, "validate mesh");
  add_option(cli, "--info,-i", info, "print mesh info");
  add_option(cli, "--path,-p", pathname, "output path");
  add_option(cli, "--stats,-s", statsname, "output stats");
  add_option(cli, "--scene,-S", scenename, "output scene");
  add_option(cli, "mesh", meshname, "input mesh", true);
  parse_cli(cli, argc, argv);

  // mesh data
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto texcoords = vector<vec2f>{};
  auto colors    = vector<vec3f>{};
  auto triangles = vector<vec3i>{};

  // stats
  auto stats                = json{};
  stats["mesh"]["filename"] = meshname;
  stats["mesh"]["valid"]    = false;

  // load mesh
  auto ioerror    = ""s;
  auto load_timer = print_timed("load mesh");
  if (!load_mesh(
          meshname, triangles, positions, normals, texcoords, colors, ioerror))
    print_fatal(ioerror);
  stats["time"]["load"] = print_elapsed(load_timer);

  // check if valid
  if (validate) {
  } else {
    stats["mesh"]["valid"] = true;
  }

  // print info
  stats["mesh"]["triangles"] = triangles.size();
  stats["mesh"]["vertices"]  = positions.size();

  // transform
  auto rescale_timer = print_timed("rescale bbox");
  auto bbox          = invalidb3f;
  for (auto& position : positions) bbox = merge(bbox, position);
  stats["time"]["rescale"] = print_elapsed(rescale_timer);

  // build bvh
  auto bvh_timer       = print_timed("build bvh");
  auto bvh             = make_triangles_bvh(triangles, positions, {});
  stats["time"]["bvh"] = print_elapsed(bvh_timer);

  // pick points
  auto points_timer       = print_timed("sample points");
  auto points             = sample_points(triangles, positions, bvh);
  stats["time"]["points"] = print_elapsed(points_timer);

  // stats
  stats["points"]["positions"] = points;

  // build graph
  auto graph_timer       = print_timed("graph bvh");
  auto graph             = make_triangles_bvh(triangles, positions, {});
  stats["time"]["graph"] = print_elapsed(graph_timer);

  // trace path
  auto path_timer = print_timed("trace path");
  // auto path             = trace_path(graph, triangles, positions, points);
  stats["time"]["path"] = print_elapsed(path_timer);

  // save path
  stats["path"]["name"] = pathname;
  // if (!save_path(pathname, path, ioerror)) return log_error(ioerror);

  // save scene
  auto scene_guard = std::make_unique<sceneio_scene>();
  auto scene       = scene_guard.get();
  // if (!path_to_scene(scene, triangles, positions, path, ioerror))
  //   return log_error(ioerror);
  stats["scene"]["name"] = scenename;
  if (!save_scene(scenename, scene, ioerror)) print_fatal(ioerror);

  // save stats
  if (!save_json(statsname, stats, ioerror)) print_fatal(ioerror);

  if (info) {
    print_info("mesh stats ------------");
    auto stats = mesh_stats(triangles, positions, normals, texcoords, colors);
    for (auto& stat : stats) print_info(stat);
  }

  // done
  return 0;
}
