//
// Implementation for Yocto/Shape Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

#include "yocto_shapeio.h"

#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>

#include "yocto_commonio.h"
#include "yocto_geometry.h"
#include "yocto_modelio.h"

// -----------------------------------------------------------------------------
// SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply mesh
shape_data load_shape(const string& filename, bool flip_texcoord) {
  auto shape = shape_data{};
  load_shape(filename, shape, flip_texcoord);
  return shape;
}
void load_shape(const string& filename, shape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    load_ply(filename, ply);
    get_positions(ply, shape.positions);
    get_normals(ply, shape.normals);
    get_texcoords(ply, shape.texcoords, flip_texcoord);
    get_colors(ply, shape.colors);
    get_radius(ply, shape.radius);
    get_faces(ply, shape.triangles, shape.quads);
    get_lines(ply, shape.lines);
    get_points(ply, shape.points);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      throw io_error::shape_error(filename);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    load_obj(filename, obj, false);
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_faces(obj, shape.triangles, shape.quads, materials);
    get_lines(obj, shape.lines, materials);
    get_points(obj, shape.points, materials);
    if (shape.points.empty() && shape.lines.empty() &&
        shape.triangles.empty() && shape.quads.empty())
      throw io_error::shape_error(filename);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    load_stl(filename, stl, true);
    if (stl.shapes.size() != 1) throw io_error::shape_error(filename);
    auto fnormals = vector<vec3f>{};
    if (!get_triangles(stl, 0, shape.triangles, shape.positions, fnormals))
      throw io_error::shape_error(filename);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_shape_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
  }
}

// Save ply mesh
void save_shape(const string& filename, const shape_data& shape,
    bool flip_texcoord, bool ascii) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_positions(ply, shape.positions);
    add_normals(ply, shape.normals);
    add_texcoords(ply, shape.texcoords, flip_texcoord);
    add_colors(ply, shape.colors);
    add_radius(ply, shape.radius);
    add_faces(ply, shape.triangles, shape.quads);
    add_lines(ply, shape.lines);
    add_points(ply, shape.points);
    save_ply(filename, ply);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.normals);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_triangles(obj, shape.triangles, 0, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_quads(
        obj, shape.quads, 0, !shape.normals.empty(), !shape.texcoords.empty());
    add_lines(
        obj, shape.lines, 0, !shape.normals.empty(), !shape.texcoords.empty());
    add_points(
        obj, shape.points, 0, !shape.normals.empty(), !shape.texcoords.empty());
    save_obj(filename, obj);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.lines.empty()) throw io_error{filename, "lines not supported"};
    if (!shape.points.empty()) throw io_error{filename, "points not supported"};
    if (!shape.triangles.empty()) {
      add_triangles(stl, shape.triangles, shape.positions, {});
    } else if (!shape.quads.empty()) {
      add_triangles(stl, quads_to_triangles(shape.quads), shape.positions, {});
    } else {
      throw io_error::shape_error(filename);
    }
    save_stl(filename, stl);
  } else if (ext == ".cpp" || ext == ".CPP") {
    auto to_cpp = [](const string& name, const string& vname,
                      const auto& values) -> string {
      using T = typename std::remove_const_t<
          std::remove_reference_t<decltype(values)>>::value_type;
      if (values.empty()) return ""s;
      auto str = "auto " + name + "_" + vname + " = ";
      if constexpr (std::is_same_v<int, T>) str += "vector<int>{\n";
      if constexpr (std::is_same_v<float, T>) str += "vector<float>{\n";
      if constexpr (std::is_same_v<vec2i, T>) str += "vector<vec2i>{\n";
      if constexpr (std::is_same_v<vec2f, T>) str += "vector<vec2f>{\n";
      if constexpr (std::is_same_v<vec3i, T>) str += "vector<vec3i>{\n";
      if constexpr (std::is_same_v<vec3f, T>) str += "vector<vec3f>{\n";
      if constexpr (std::is_same_v<vec4i, T>) str += "vector<vec4i>{\n";
      if constexpr (std::is_same_v<vec4f, T>) str += "vector<vec4f>{\n";
      for (auto& value : values) {
        if constexpr (std::is_same_v<int, T> || std::is_same_v<float, T>) {
          str += std::to_string(value) + ",\n";
        } else if constexpr (std::is_same_v<vec2i, T> ||
                             std::is_same_v<vec2f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "},\n";
        } else if constexpr (std::is_same_v<vec3i, T> ||
                             std::is_same_v<vec3f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "},\n";
        } else if constexpr (std::is_same_v<vec4i, T> ||
                             std::is_same_v<vec4f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "," + std::to_string(value.w) +
                 "},\n";
        } else {
          throw std::invalid_argument{"cannot print this"};
        }
      }
      str += "};\n\n";
      return str;
    };

    auto name = string{"shape"};
    auto str  = ""s;
    str += to_cpp(name, "positions", shape.positions);
    str += to_cpp(name, "normals", shape.normals);
    str += to_cpp(name, "texcoords", shape.texcoords);
    str += to_cpp(name, "colors", shape.colors);
    str += to_cpp(name, "radius", shape.radius);
    str += to_cpp(name, "points", shape.points);
    str += to_cpp(name, "lines", shape.lines);
    str += to_cpp(name, "triangles", shape.triangles);
    str += to_cpp(name, "quads", shape.quads);
    save_text(filename, str);
  } else {
    throw io_error::format_error(filename);
  }
}

// Load face-varying mesh
fvshape_data load_fvshape(const string& filename, bool flip_texcoord) {
  auto shape = fvshape_data{};
  load_fvshape(filename, shape, flip_texcoord);
  return shape;
}
void load_fvshape(
    const string& filename, fvshape_data& shape, bool flip_texcoord) {
  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    load_ply(filename, ply);
    get_positions(ply, shape.positions);
    get_normals(ply, shape.normals);
    get_texcoords(ply, shape.texcoords, flip_texcoord);
    get_quads(ply, shape.quadspos);
    if (!shape.normals.empty()) shape.quadsnorm = shape.quadspos;
    if (!shape.texcoords.empty()) shape.quadstexcoord = shape.quadspos;
    if (shape.quadspos.empty()) throw io_error::shape_error(filename);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    load_obj(filename, obj, true);
    auto materials = vector<int>{};
    get_positions(obj, shape.positions);
    get_normals(obj, shape.normals);
    get_texcoords(obj, shape.texcoords, flip_texcoord);
    get_fvquads(
        obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, materials);
    if (shape.quadspos.empty()) throw io_error::shape_error(filename);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    load_stl(filename, stl, true);
    if (stl.shapes.empty()) throw io_error::shape_error(filename);
    if (stl.shapes.size() > 1) throw io_error::shape_error(filename);
    auto fnormals  = vector<vec3f>{};
    auto triangles = vector<vec3i>{};
    if (!get_triangles(stl, 0, triangles, shape.positions, fnormals))
      throw io_error::shape_error(filename);
    shape.quadspos = triangles_to_quads(triangles);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    shape = make_fvshape_preset(path_basename(filename));
  } else {
    throw io_error::format_error(filename);
  }
}

// Save ply mesh
void save_fvshape(const string& filename, const fvshape_data& shape,
    bool flip_texcoord, bool ascii) {
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply             = ply_model{};
    auto split_quads     = vector<vec4i>{};
    auto split_positions = vector<vec3f>{};
    auto split_normals   = vector<vec3f>{};
    auto split_texcoords = vector<vec2f>{};
    split_facevarying(split_quads, split_positions, split_normals,
        split_texcoords, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
        shape.positions, shape.normals, shape.texcoords);
    add_positions(ply, split_positions);
    add_normals(ply, split_normals);
    add_texcoords(ply, split_texcoords, flip_texcoord);
    add_faces(ply, {}, split_quads);
    save_ply(filename, ply);
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, shape.positions);
    add_normals(obj, shape.positions);
    add_texcoords(obj, shape.texcoords, flip_texcoord);
    add_fvquads(obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, 0);
    save_obj(filename, obj);
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (!shape.quadspos.empty()) {
      auto split_quads     = vector<vec4i>{};
      auto split_positions = vector<vec3f>{};
      auto split_normals   = vector<vec3f>{};
      auto split_texcoords = vector<vec2f>{};
      split_facevarying(split_quads, split_positions, split_normals,
          split_texcoords, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords);
      add_triangles(stl, quads_to_triangles(split_quads), split_positions, {});
    } else {
      throw io_error::shape_error(filename);
    }
    save_stl(filename, stl);
  } else if (ext == ".cpp" || ext == ".CPP") {
    auto to_cpp = [](const string& name, const string& vname,
                      const auto& values) -> string {
      using T = typename std::remove_const_t<
          std::remove_reference_t<decltype(values)>>::value_type;
      if (values.empty()) return ""s;
      auto str = "auto " + name + "_" + vname + " = ";
      if constexpr (std::is_same_v<int, T>) str += "vector<int>{\n";
      if constexpr (std::is_same_v<float, T>) str += "vector<float>{\n";
      if constexpr (std::is_same_v<vec2i, T>) str += "vector<vec2i>{\n";
      if constexpr (std::is_same_v<vec2f, T>) str += "vector<vec2f>{\n";
      if constexpr (std::is_same_v<vec3i, T>) str += "vector<vec3i>{\n";
      if constexpr (std::is_same_v<vec3f, T>) str += "vector<vec3f>{\n";
      if constexpr (std::is_same_v<vec4i, T>) str += "vector<vec4i>{\n";
      if constexpr (std::is_same_v<vec4f, T>) str += "vector<vec4f>{\n";
      for (auto& value : values) {
        if constexpr (std::is_same_v<int, T> || std::is_same_v<float, T>) {
          str += std::to_string(value) + ",\n";
        } else if constexpr (std::is_same_v<vec2i, T> ||
                             std::is_same_v<vec2f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "},\n";
        } else if constexpr (std::is_same_v<vec3i, T> ||
                             std::is_same_v<vec3f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "},\n";
        } else if constexpr (std::is_same_v<vec4i, T> ||
                             std::is_same_v<vec4f, T>) {
          str += "{" + std::to_string(value.x) + "," + std::to_string(value.y) +
                 "," + std::to_string(value.z) + "," + std::to_string(value.w) +
                 "},\n";
        } else {
          throw std::invalid_argument{"cannot print this"};
        }
      }
      str += "};\n\n";
      return str;
    };
    auto name = string{"shape"};
    auto str  = ""s;
    str += to_cpp(name, "positions", shape.positions);
    str += to_cpp(name, "normals", shape.normals);
    str += to_cpp(name, "texcoords", shape.texcoords);
    str += to_cpp(name, "quadspos", shape.quadspos);
    str += to_cpp(name, "quadsnorm", shape.quadsnorm);
    str += to_cpp(name, "quadstexcoord", shape.quadstexcoord);
    save_text(filename, str);
  } else {
    throw io_error::format_error(filename);
  }
}

// Shape presets used for testing.
shape_data make_shape_preset(const string& type) {
  if (type == "default-quad") {
    return make_rect();
  } else if (type == "default-quady") {
    return make_recty();
  } else if (type == "default-cube") {
    return make_box();
  } else if (type == "default-cube-rounded") {
    return make_rounded_box();
  } else if (type == "default-sphere") {
    return make_sphere();
  } else if (type == "default-matcube") {
    return make_rounded_box();
  } else if (type == "default-matsphere") {
    return make_uvspherey();
  } else if (type == "default-disk") {
    return make_disk();
  } else if (type == "default-disk-bulged") {
    return make_bulged_disk();
  } else if (type == "default-quad-bulged") {
    return make_bulged_rect();
  } else if (type == "default-uvsphere") {
    return make_uvsphere();
  } else if (type == "default-uvsphere-flipcap") {
    return make_capped_uvsphere();
  } else if (type == "default-uvspherey") {
    return make_uvspherey();
  } else if (type == "default-uvspherey-flipcap") {
    return make_capped_uvspherey();
  } else if (type == "default-uvdisk") {
    return make_uvdisk();
  } else if (type == "default-uvcylinder") {
    return make_uvcylinder();
  } else if (type == "default-uvcylinder-rounded") {
    return make_rounded_uvcylinder({32, 32, 32});
  } else if (type == "default-geosphere") {
    return make_geosphere();
  } else if (type == "default-floor") {
    return make_floor();
  } else if (type == "default-floor-bent") {
    return make_bent_floor();
  } else if (type == "default-matball") {
    return make_sphere();
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8f);
    return make_hair(base, {4, 65536}, {0.2f, 0.2f}, {0.002f, 0.001f});
  } else if (type == "default-hairball-interior") {
    return make_sphere(pow2(5), 0.8f);
  } else if (type == "default-suzanne") {
    return make_monkey();
  } else if (type == "default-cube-facevarying") {
    return fvshape_to_shape(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    return fvshape_to_shape(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    return make_recty({256, 256});
  } else if (type == "default-sphere-displaced") {
    return make_sphere(128);
  } else if (type == "test-cube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvsphere") {
    auto shape = make_uvsphere({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvsphere-flipcap") {
    auto shape = make_capped_uvsphere({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvspherey") {
    auto shape = make_uvspherey({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvspherey-flipcap") {
    auto shape = make_capped_uvspherey({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-sphere") {
    auto shape = make_sphere(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-matcube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-matsphere") {
    auto shape = make_uvspherey({32, 32}, 0.075f, {2, 1});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-sphere-displaced") {
    auto shape = make_sphere(128, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-smallsphere") {
    auto shape = make_sphere(32, 0.015f, 1);
    for (auto& p : shape.positions) p += {0, 0.015f, 0};
    return shape;
  } else if (type == "test-disk") {
    auto shape = make_disk(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-uvcylinder") {
    auto shape = make_rounded_uvcylinder(
        {32, 32, 32}, {0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-floor") {
    return make_floor({1, 1}, {2, 2}, {20, 20});
  } else if (type == "test-smallfloor") {
    return make_floor({1, 1}, {0.5f, 0.5f}, {1, 1});
  } else if (type == "test-quad") {
    return make_rect({1, 1}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quady") {
    return make_recty({1, 1}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quad-displaced") {
    return make_rect({256, 256}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-quady-displaced") {
    return make_recty({256, 256}, {0.075f, 0.075f}, {1, 1});
  } else if (type == "test-matball") {
    auto shape = make_sphere(32, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-geosphere") {
    auto shape = make_geosphere(0.075f, 3);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-geosphere-flat") {
    auto shape = make_geosphere(0.075f, 3);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    shape.normals = {};
    return shape;
  } else if (type == "test-geosphere-subdivided") {
    auto shape = make_geosphere(0.075f, 6);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03f, 100});
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f});
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075f, 0};
    return make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128});
  } else if (type == "test-hairball-interior") {
    auto shape = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-suzanne-subdiv") {
    auto shape = make_monkey(0.075f * 0.8f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-cube-subdiv") {
    auto fvshape    = make_fvcube(0.075f);
    auto shape      = shape_data{};
    shape.quads     = fvshape.quadspos;
    shape.positions = fvshape.positions;
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-arealight1") {
    return make_rect({1, 1}, {0.2f, 0.2f});
  } else if (type == "test-arealight2") {
    return make_rect({1, 1}, {0.2f, 0.2f});
  } else if (type == "test-largearealight1") {
    return make_rect({1, 1}, {0.4f, 0.4f});
  } else if (type == "test-largearealight2") {
    return make_rect({1, 1}, {0.4f, 0.4f});
  } else if (type == "test-pointlight1") {
    return make_point(0);
  } else if (type == "test-pointlight2") {
    return make_point(0);
  } else if (type == "test-point") {
    return make_points(1);
  } else if (type == "test-points") {
    return make_points(4096);
  } else if (type == "test-points-random") {
    auto shape = make_random_points(4096, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape;
  } else if (type == "test-points-grid") {
    auto shape = make_points({256, 256}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f;
    return shape;
  } else if (type == "test-lines-grid") {
    auto shape = make_lines({256, 256}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f;
    return shape;
  } else if (type == "test-thickpoints-grid") {
    auto shape = make_points({16, 16}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f * 10;
    return shape;
  } else if (type == "test-thicklines-grid") {
    auto shape = make_lines({16, 16}, {0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    for (auto& r : shape.radius) r *= 0.075f * 10;
    return shape;
  } else if (type == "test-particles") {
    return make_points(4096);
  } else if (type == "test-cloth") {
    return make_rect({64, 64}, {0.2f, 0.2f});
  } else if (type == "test-clothy") {
    return make_recty({64, 64}, {0.2f, 0.2f});
  } else {
    throw io_error::preset_error(type);
  }
}

// Shape presets used for testing.
fvshape_data make_fvshape_preset(const string& type) {
  if (type == "default-quad") {
    return shape_to_fvshape(make_rect());
  } else if (type == "default-quady") {
    return shape_to_fvshape(make_recty());
  } else if (type == "default-cube") {
    return shape_to_fvshape(make_box());
  } else if (type == "default-cube-rounded") {
    return shape_to_fvshape(make_rounded_box());
  } else if (type == "default-sphere") {
    return shape_to_fvshape(make_sphere());
  } else if (type == "default-matcube") {
    return shape_to_fvshape(make_rounded_box());
  } else if (type == "default-matsphere") {
    return shape_to_fvshape(make_uvspherey());
  } else if (type == "default-disk") {
    return shape_to_fvshape(make_disk());
  } else if (type == "default-disk-bulged") {
    return shape_to_fvshape(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    return shape_to_fvshape(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    return shape_to_fvshape(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    return shape_to_fvshape(make_capped_uvsphere());
  } else if (type == "default-uvspherey") {
    return shape_to_fvshape(make_uvspherey());
  } else if (type == "default-uvspherey-flipcap") {
    return shape_to_fvshape(make_capped_uvspherey());
  } else if (type == "default-uvdisk") {
    return shape_to_fvshape(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    return shape_to_fvshape(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    return shape_to_fvshape(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    return shape_to_fvshape(make_geosphere());
  } else if (type == "default-floor") {
    return shape_to_fvshape(make_floor());
  } else if (type == "default-floor-bent") {
    return shape_to_fvshape(make_bent_floor());
  } else if (type == "default-matball") {
    return shape_to_fvshape(make_sphere());
  } else if (type == "default-hairball-interior") {
    return shape_to_fvshape(make_sphere(pow2(5), 0.8f));
  } else if (type == "default-suzanne") {
    return shape_to_fvshape(make_monkey());
  } else if (type == "default-cube-facevarying") {
    return make_fvbox();
  } else if (type == "default-sphere-facevarying") {
    return make_fvsphere();
  } else if (type == "default-quady-displaced") {
    return shape_to_fvshape(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    return shape_to_fvshape(make_sphere(128));
  } else if (type == "test-cube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-matsphere") {
    auto shape = make_uvspherey({32, 32}, 0.075f, {2, 1});
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvsphere") {
    auto shape = make_uvsphere({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvsphere-flipcap") {
    auto shape = make_capped_uvsphere({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvspherey") {
    auto shape = make_uvspherey({32, 32}, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvspherey-flipcap") {
    auto shape = make_capped_uvspherey({32, 32}, 0.075f, {1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-sphere") {
    auto shape = make_sphere(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-sphere-displaced") {
    auto shape = make_sphere(128, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-matcube") {
    auto shape = make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-disk") {
    auto shape = make_disk(32, 0.075f, 1);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-uvcylinder") {
    auto shape = make_rounded_uvcylinder(
        {32, 32, 32}, {0.075f, 0.075f}, {1, 1, 1}, 0.3f * 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-floor") {
    return shape_to_fvshape(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-smallfloor") {
    return shape_to_fvshape(make_floor({1, 1}, {0.5f, 0.5f}, {1, 1}));
  } else if (type == "test-quad") {
    return shape_to_fvshape(make_rect({1, 1}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quady") {
    return shape_to_fvshape(make_recty({1, 1}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    return shape_to_fvshape(make_rect({256, 256}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    return shape_to_fvshape(make_recty({256, 256}, {0.075f, 0.075f}, {1, 1}));
  } else if (type == "test-matball") {
    auto shape = make_sphere(32, 0.075f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-suzanne-subdiv") {
    auto shape = make_monkey(0.075f * 0.8f);
    for (auto& p : shape.positions) p += {0, 0.075f, 0};
    return shape_to_fvshape(shape);
  } else if (type == "test-cube-subdiv") {
    auto fvshape = make_fvcube(0.075f);
    for (auto& p : fvshape.positions) p += {0, 0.075f, 0};
    return fvshape;
  } else if (type == "test-arealight1") {
    return shape_to_fvshape(make_rect({1, 1}, {0.2f, 0.2f}));
  } else if (type == "test-arealight2") {
    return shape_to_fvshape(make_rect({1, 1}, {0.2f, 0.2f}));
  } else if (type == "test-largearealight1") {
    return shape_to_fvshape(make_rect({1, 1}, {0.4f, 0.4f}));
  } else if (type == "test-largearealight2") {
    return shape_to_fvshape(make_rect({1, 1}, {0.4f, 0.4f}));
  } else if (type == "test-cloth") {
    return shape_to_fvshape(make_rect({64, 64}, {0.2f, 0.2f}));
  } else if (type == "test-clothy") {
    return shape_to_fvshape(make_recty({64, 64}, {0.2f, 0.2f}));
  } else {
    throw io_error::preset_error(type);
  }
}

// Load ply mesh
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoord) {
  try {
    load_shape(filename, shape, flip_texcoord);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save ply mesh
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoord, bool ascii) {
  try {
    save_shape(filename, shape, flip_texcoord, ascii);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Load ply mesh
bool load_fvshape(const string& filename, fvshape_data& fvshape, string& error,
    bool flip_texcoord) {
  try {
    load_fvshape(filename, fvshape, flip_texcoord);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Save ply mesh
bool save_fvshape(const string& filename, const fvshape_data& fvshape,
    string& error, bool flip_texcoord, bool ascii) {
  try {
    save_fvshape(filename, fvshape, flip_texcoord, ascii);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Shape presets used ofr testing.
bool make_shape_preset(shape_data& shape, const string& type, string& error) {
  try {
    shape = make_shape_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

// Shape presets used for testing.
bool make_fvshape_preset(
    fvshape_data& fvshape, const string& type, string& error) {
  try {
    fvshape = make_fvshape_preset(type);
    return true;
  } catch (const io_error& exception) {
    error = exception.what();
    return false;
  }
}

}  // namespace yocto
