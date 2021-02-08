//
// Implementation for Yocto/Shape
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_shape.h"

#include <algorithm>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>

#include "yocto_commonio.h"
#include "yocto_geometry.h"
#include "yocto_modelio.h"
#include "yocto_noise.h"
#include "yocto_sampling.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::deque;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const shape_data& shape, int element, const vec2f& uv) {
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.positions[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.positions[line.x], shape.positions[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return normalize(shape.normals[point]);
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return normalize(
        interpolate_line(shape.normals[line.x], shape.normals[line.y], uv.x));
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return normalize(interpolate_triangle(shape.normals[triangle.x],
        shape.normals[triangle.y], shape.normals[triangle.z], uv));
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec3f eval_tangent(const shape_data& shape, int element, const vec2f& uv) {
  return eval_normal(shape, element, uv);
}

vec2f eval_texcoord(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.texcoords[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.texcoords[line.x], shape.texcoords[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.texcoords[triangle.x],
        shape.texcoords[triangle.y], shape.texcoords[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

vec4f eval_color(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.colors[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.colors[line.x], shape.colors[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.colors[triangle.x],
        shape.colors[triangle.y], shape.colors[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.colors[quad.x], shape.colors[quad.y],
        shape.colors[quad.z], shape.colors[quad.w], uv);
  } else {
    return {0, 0};
  }
}

float eval_radius(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.radius.empty()) return 0;
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.radius[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.radius[line.x], shape.radius[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.radius[triangle.x],
        shape.radius[triangle.y], shape.radius[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.radius[quad.x], shape.radius[quad.y],
        shape.radius[quad.z], shape.radius[quad.w], uv);
  } else {
    return 0;
  }
}

// Evaluate element normals
vec3f eval_element_normal(const shape_data& shape, int element) {
  if (!shape.points.empty()) {
    return {0, 0, 1};
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return line_tangent(shape.positions[line.x], shape.positions[line.y]);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return triangle_normal(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z]);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const shape_data& shape) {
  if (!shape.points.empty()) {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    return lines_tangents(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return triangles_normals(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return quads_normals(shape.quads, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const shape_data& shape) {
  if (!shape.points.empty()) {
    normals.assign(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    lines_tangents(normals, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    triangles_normals(normals, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    quads_normals(normals, shape.quads, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Shape sampling
vector<float> sample_shape_cdf(const shape_data& shape) {
  if (!shape.points.empty()) {
    return sample_points_cdf((int)shape.points.size());
  } else if (!shape.lines.empty()) {
    return sample_lines_cdf(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return sample_triangles_cdf(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return sample_quads_cdf(shape.quads, shape.positions);
  } else {
    return sample_points_cdf((int)shape.positions.size());
  }
}

void sample_shape_cdf(vector<float>& cdf, const shape_data& shape) {
  if (!shape.points.empty()) {
    sample_points_cdf(cdf, (int)shape.points.size());
  } else if (!shape.lines.empty()) {
    sample_lines_cdf(cdf, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    sample_triangles_cdf(cdf, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    sample_quads_cdf(cdf, shape.quads, shape.positions);
  } else {
    sample_points_cdf(cdf, (int)shape.positions.size());
  }
}

shape_point sample_shape(const shape_data& shape, const vector<float>& cdf,
    float rn, const vec2f& ruv) {
  if (!shape.points.empty()) {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  } else if (!shape.lines.empty()) {
    auto [element, u] = sample_lines(cdf, rn, ruv.x);
    return {element, {u, 0}};
  } else if (!shape.triangles.empty()) {
    auto [element, uv] = sample_triangles(cdf, rn, ruv);
    return {element, uv};
  } else if (!shape.quads.empty()) {
    auto [element, uv] = sample_quads(cdf, rn, ruv);
    return {element, uv};
  } else {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  }
}

vector<shape_point> sample_shape(
    const shape_data& shape, int num_samples, uint64_t seed) {
  auto cdf    = sample_shape_cdf(shape);
  auto points = vector<shape_point>(num_samples);
  auto rng    = make_rng(seed);
  for (auto& point : points) {
    point = sample_shape(shape, cdf, rand1f(rng), rand2f(rng));
  }
  return points;
}

// Conversions
shape_data quads_to_triangles(const shape_data& shape) {
  auto result = shape;
  quads_to_triangles(result, result);
  return result;
}
void quads_to_triangles(shape_data& result, const shape_data& shape) {
  result.triangles = quads_to_triangles(shape.quads);
  result.quads     = {};
}

// Subdivision
shape_data subdivide_shape(
    const shape_data& shape, int subdivisions, bool catmullclark) {
  // This should probably be reimplemented in a faster fashion,
  // but how it is not obvious
  auto subdivided = shape_data{};
  if (!shape.points.empty()) {
    // nothing to do
  } else if (!shape.lines.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_lines(
        shape.lines, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_lines(
        shape.lines, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_lines(
        shape.lines, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_lines(
        shape.lines, shape.radius, subdivisions);
    std::tie(subdivided.lines, subdivided.positions) = subdivide_lines(
        shape.lines, shape.positions, subdivisions);
  } else if (!shape.triangles.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_triangles(
        shape.triangles, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_triangles(
        shape.triangles, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_triangles(
        shape.triangles, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_triangles(
        shape.triangles, shape.radius, subdivisions);
    std::tie(subdivided.triangles, subdivided.positions) = subdivide_triangles(
        shape.triangles, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && !catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_quads(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_quads(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_quads(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_quads(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_quads(
        shape.quads, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_catmullclark(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_catmullclark(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_catmullclark(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_catmullclark(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_catmullclark(
        shape.quads, shape.positions, subdivisions);
  } else {
    // empty shape
  }
  return subdivided;
}

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, const vec2f& uv) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadsnorm[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec2f eval_texcoord(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadstexcoord[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

// Evaluate element normals
vec3f eval_element_normal(const fvshape_data& shape, int element) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const fvshape_data& shape) {
  if (!shape.quadspos.empty()) {
    return quads_normals(shape.quadspos, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const fvshape_data& shape) {
  if (!shape.quadspos.empty()) {
    quads_normals(normals, shape.quadspos, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Conversions
shape_data fvshape_to_shape(const fvshape_data& fvshape, bool as_triangles) {
  auto shape = shape_data{};
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, fvshape.quadspos, fvshape.quadsnorm,
      fvshape.quadstexcoord, fvshape.positions, fvshape.normals,
      fvshape.texcoords);
  return shape;
}
fvshape_data shape_to_fvshape(const shape_data& shape) {
  if (!shape.points.empty() || !shape.lines.empty())
    throw std::invalid_argument{"cannor convert shape"};
  auto fvshape          = fvshape_data{};
  fvshape.positions     = shape.positions;
  fvshape.normals       = shape.normals;
  fvshape.texcoords     = shape.texcoords;
  fvshape.quadspos      = !shape.quads.empty() ? shape.quads
                                               : triangles_to_quads(shape.triangles);
  fvshape.quadsnorm     = !shape.normals.empty() ? fvshape.quadspos
                                                 : vector<vec4i>{};
  fvshape.quadstexcoord = !shape.texcoords.empty() ? fvshape.quadspos
                                                   : vector<vec4i>{};
  return fvshape;
}

// Subdivision
fvshape_data subdivide_fvshape(
    const fvshape_data& shape, int subdivisions, bool catmullclark) {
  auto subdivided = fvshape_data{};
  if (!catmullclark) {
    std::tie(subdivided.quadspos, subdivided.positions) = subdivide_quads(
        shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_quads(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) = subdivide_quads(
        shape.quadstexcoord, shape.texcoords, subdivisions);
  } else {
    std::tie(subdivided.quadspos, subdivided.positions) =
        subdivide_catmullclark(shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_catmullclark(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) =
        subdivide_catmullclark(
            shape.quadstexcoord, shape.texcoords, subdivisions, true);
  }
  return subdivided;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape from/to OBJ.
static bool load_obj_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord);
static bool save_obj_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load/save a shape from/to PLY.
static bool load_ply_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord);
static bool save_ply_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load/save a shape from/to STL.
static bool load_stl_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord);
static bool save_stl_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Save a shape to C++.
static bool save_cpp_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load ply mesh
bool load_shape(const string& filename, shape_data& shape, string& error,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    return load_ply_shape(filename, shape, error, flip_texcoord);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_shape(filename, shape, error, flip_texcoord);
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_shape(filename, shape, error, flip_texcoord);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    if (!make_shape_preset(shape, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_shape(const string& filename, const shape_data& shape, string& error,
    bool flip_texcoord, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    return save_ply_shape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_shape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_shape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".cpp" || ext == ".CPP") {
    return save_cpp_shape(filename, shape, error, flip_texcoord, ascii);
  } else {
    return format_error();
  }
}

// Load/save a shape from/to OBJ.
static bool load_obj_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  shape = {};

  // load obj
  auto obj = obj_shape{};
  if (!load_obj(filename, obj, error, true)) return false;

  // decide what to do and get properties
  auto materials = vector<int>{};
  get_positions(obj, shape.positions);
  get_normals(obj, shape.normals);
  get_texcoords(obj, shape.texcoords, flip_texcoord);
  get_faces(obj, shape.triangles, shape.quads, materials);
  get_lines(obj, shape.lines, materials);
  get_points(obj, shape.points, materials);

  if (shape.points.empty() && shape.lines.empty() && shape.triangles.empty() &&
      shape.quads.empty())
    return shape_error();
  return true;
}
static bool save_obj_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
    return false;
  };

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
  return save_obj(filename, obj, error);
}

// Load/save a shape from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  shape = {};

  // open ply
  auto ply = ply_model{};
  if (!load_ply(filename, ply, error)) return false;

  // get vertex
  get_positions(ply, shape.positions);
  get_normals(ply, shape.normals);
  get_texcoords(ply, shape.texcoords, flip_texcoord);
  get_colors(ply, shape.colors);
  get_radius(ply, shape.radius);

  // get faces
  get_faces(ply, shape.triangles, shape.quads);
  get_lines(ply, shape.lines);
  get_points(ply, shape.points);

  if (shape.points.empty() && shape.lines.empty() && shape.triangles.empty() &&
      shape.quads.empty())
    return shape_error();
  return true;
}
static bool save_ply_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  // create ply
  auto ply = ply_model{};
  add_positions(ply, shape.positions);
  add_normals(ply, shape.normals);
  add_texcoords(ply, shape.texcoords, flip_texcoord);
  add_colors(ply, shape.colors);
  add_radius(ply, shape.radius);
  add_faces(ply, shape.triangles, shape.quads);
  add_lines(ply, shape.lines);
  add_points(ply, shape.points);
  return save_ply(filename, ply, error);
}

// Load/save a shape from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_shape(const string& filename, shape_data& shape,
    string& error, bool flip_texcoord) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  // load stl
  auto stl = stl_model{};
  if (!load_stl(filename, stl, error, true)) return false;

  // get shape
  if (stl.shapes.size() != 1) return shape_error();
  auto fnormals = vector<vec3f>{};
  if (!get_triangles(stl, 0, shape.triangles, shape.positions, fnormals))
    return shape_error();
  return true;
}
static bool save_stl_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  // create stl
  auto stl = stl_model{};
  if (!shape.lines.empty()) return line_error();
  if (!shape.points.empty()) return point_error();
  if (!shape.triangles.empty()) {
    add_triangles(stl, shape.triangles, shape.positions, {});
  } else if (!shape.quads.empty()) {
    add_triangles(stl, quads_to_triangles(shape.quads), shape.positions, {});
  } else {
    return shape_error();
  }
  return save_stl(filename, stl, error);
}

// Save a shape to C++.
static bool save_cpp_shape(const string& filename, const shape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
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
  return save_text(filename, str, error);
}

// Load/save a shape from/to OBJ.
static bool load_obj_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord);
static bool save_obj_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load/save a shape from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord);
static bool save_ply_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load/save a shape from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord);
static bool save_stl_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Save a shape to C++.
static bool save_cpp_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii);

// Load ply mesh
bool load_fvshape(const string& filename, fvshape_data& shape, string& error,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  shape = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    return load_ply_fvshape(filename, shape, error, flip_texcoord);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_fvshape(filename, shape, error, flip_texcoord);
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_fvshape(filename, shape, error, flip_texcoord);
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    // create preset
    if (!make_fvshape_preset(shape, path_basename(filename), error))
      return preset_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    return save_ply_fvshape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_fvshape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_fvshape(filename, shape, error, flip_texcoord, ascii);
  } else if (ext == ".cpp" || ext == ".CPP") {
    return save_cpp_fvshape(filename, shape, error, flip_texcoord, ascii);
  } else {
    return format_error();
  }
}

// Load/save a shape from/to OBJ.
static bool load_obj_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  shape = {};

  // load obj
  auto obj = obj_shape{};
  if (!load_obj(filename, obj, error, true)) return false;
  auto materials = vector<int>{};
  get_positions(obj, shape.positions);
  get_normals(obj, shape.normals);
  get_texcoords(obj, shape.texcoords, flip_texcoord);
  get_fvquads(
      obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, materials);

  if (shape.quadspos.empty()) return shape_error();
  return true;
}
static bool save_obj_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
    return false;
  };

  auto obj = obj_shape{};
  add_positions(obj, shape.positions);
  add_normals(obj, shape.positions);
  add_texcoords(obj, shape.texcoords, flip_texcoord);
  add_fvquads(obj, shape.quadspos, shape.quadsnorm, shape.quadstexcoord, 0);
  return save_obj(filename, obj, error);
}

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other
// data.
static bool load_ply_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  shape = {};

  // open ply
  auto ply = ply_model{};
  if (!load_ply(filename, ply, error)) return false;

  get_positions(ply, shape.positions);
  get_normals(ply, shape.normals);
  get_texcoords(ply, shape.texcoords, flip_texcoord);
  get_quads(ply, shape.quadspos);
  if (!shape.normals.empty()) shape.quadsnorm = shape.quadspos;
  if (!shape.texcoords.empty()) shape.quadstexcoord = shape.quadspos;

  if (shape.positions.empty()) return shape_error();
  return true;
}
static bool save_ply_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  // create ply
  auto ply = ply_model{};
  // split data
  auto [split_quads, split_positions, split_normals, split_texcoords] =
      split_facevarying(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
          shape.positions, shape.normals, shape.texcoords);
  add_positions(ply, split_positions);
  add_normals(ply, split_normals);
  add_texcoords(ply, split_texcoords, flip_texcoord);
  add_faces(ply, {}, split_quads);
  return save_ply(filename, ply, error);
}

// Load/save a scene from/to STL. Loads/saves only one mesh with no other
// data.
static bool load_stl_fvshape(const string& filename, fvshape_data& shape,
    string& error, bool flip_texcoord) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  shape = {};

  // load obj
  auto stl = stl_model{};
  if (!load_stl(filename, stl, error, true)) return false;

  // get shape
  if (stl.shapes.empty()) return shape_error();
  if (stl.shapes.size() > 1) return shape_error();
  auto fnormals  = vector<vec3f>{};
  auto triangles = vector<vec3i>{};
  if (!get_triangles(stl, 0, triangles, shape.positions, fnormals))
    return shape_error();
  shape.quadspos = triangles_to_quads(triangles);
  return true;
}
static bool save_stl_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto line_error = [filename, &error]() {
    error = filename + ": unsupported lines";
    return false;
  };
  auto point_error = [filename, &error]() {
    error = filename + ": unsupported points";
    return false;
  };

  // create ply
  auto stl = stl_model{};
  if (!shape.quadspos.empty()) {
    // split data
    auto [split_quads, split_positions, split_normals, split_texcoords] =
        split_facevarying(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords);
    add_triangles(stl, quads_to_triangles(split_quads), split_positions, {});
  } else {
    return shape_error();
  }
  return save_stl(filename, stl, error);
}

// Save a shape to C++.
static bool save_cpp_fvshape(const string& filename, const fvshape_data& shape,
    string& error, bool flip_texcoord, bool ascii) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

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
  return save_text(filename, str, error);
}

// Shape presets used ofr testing.
bool make_shape_preset(shape_data& shape, const string& type, string& error) {
  auto set_quads = [&](quads_shape&& shape_) {
    shape.quads     = shape_.quads;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };
  auto set_triangles = [&](triangles_shape&& shape_) {
    shape.triangles = shape_.triangles;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };
  auto set_lines = [&](lines_shape&& shape_) {
    shape.lines     = shape_.lines;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
    shape.radius    = shape_.radius;
  };
  auto set_points = [&](points_shape&& shape_) {
    shape.points    = shape_.points;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
    shape.radius    = shape_.radius;
  };
  auto set_fvquads = [&](quads_fvshape&& shape_) {
    shape.quads     = shape_.quadspos;
    shape.positions = shape_.positions;
    shape.normals   = shape_.normals;
    shape.texcoords = shape_.texcoords;
  };

  if (type == "default-quad") {
    set_quads(make_rect());
  } else if (type == "default-quady") {
    set_quads(make_recty());
  } else if (type == "default-cube") {
    set_quads(make_box());
  } else if (type == "default-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "default-sphere") {
    set_quads(make_sphere());
  } else if (type == "default-disk") {
    set_quads(make_disk());
  } else if (type == "default-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "default-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "default-floor") {
    set_quads(make_floor());
  } else if (type == "default-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "default-matball") {
    set_quads(make_sphere());
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8);
    set_lines(make_hair(base, {4, 65536}, {0.2, 0.2}, {0.002, 0.001}));
  } else if (type == "default-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "default-suzanne") {
    set_quads(make_monkey());
  } else if (type == "default-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    set_quads(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    set_quads(make_sphere(128));
  } else if (type == "test-cube") {
    set_quads(make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    set_quads(make_uvsphere({32, 32}, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    set_quads(make_sphere(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    set_quads(make_sphere(128, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    set_quads(make_disk(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    set_quads(make_rounded_uvcylinder(
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    set_quads(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-quad") {
    set_quads(make_rect({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady") {
    set_quads(make_recty({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    set_quads(make_rect({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    set_quads(make_recty({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-matball") {
    set_quads(make_sphere(32, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100}));
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}));
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128}));
  } else if (type == "test-hairball-interior") {
    set_quads(make_sphere(32, 0.075f * 0.8f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    set_quads(make_monkey(0.075f * 0.8f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // set_quads(make_cube( 0.075f);
    set_fvquads(make_fvcube(0.075f));
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-arealight2") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-largearealight1") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-largearealight2") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-pointlight1") {
    set_points(make_point(0));
  } else if (type == "test-pointlight2") {
    set_points(make_point(0));
  } else if (type == "test-point") {
    set_points(make_points(1));
  } else if (type == "test-points") {
    set_points(make_points(4096));
  } else if (type == "test-points-random") {
    set_points(make_random_points(4096, {0.2, 0.2, 0.2}));
  } else if (type == "test-particles") {
    set_points(make_points(4096));
  } else if (type == "test-cloth") {
    set_quads(make_rect({64, 64}, {0.2, 0.2}));
  } else if (type == "test-clothy") {
    set_quads(make_recty({64, 64}, {0.2, 0.2}));
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

// Shape presets used for testing.
bool make_fvshape_preset(
    fvshape_data& shape, const string& type, string& error) {
  auto set_quads = [&](quads_shape&& shape_) {
    shape.quadspos  = shape_.quads;
    shape.positions = shape_.positions;
    if (!shape_.normals.empty()) shape.quadsnorm = shape_.quads;
    shape.normals = shape_.normals;
    if (!shape_.texcoords.empty()) shape.quadstexcoord = shape_.quads;
    shape.texcoords = shape_.texcoords;
  };
  auto set_triangles = [&](triangles_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_lines = [&](lines_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_points = [&](points_shape&& shape) {
    throw std::invalid_argument{"bad shape type"};
  };
  auto set_fvquads = [&](quads_fvshape&& shape_) {
    shape.quadspos      = shape_.quadspos;
    shape.quadsnorm     = shape_.quadsnorm;
    shape.quadstexcoord = shape_.quadstexcoord;
    shape.positions     = shape_.positions;
    shape.normals       = shape_.normals;
    shape.texcoords     = shape_.texcoords;
  };

  if (type == "default-quad") {
    set_quads(make_rect());
  } else if (type == "default-quady") {
    set_quads(make_recty());
  } else if (type == "default-cube") {
    set_quads(make_box());
  } else if (type == "default-cube-rounded") {
    set_quads(make_rounded_box());
  } else if (type == "default-sphere") {
    set_quads(make_sphere());
  } else if (type == "default-disk") {
    set_quads(make_disk());
  } else if (type == "default-disk-bulged") {
    set_quads(make_bulged_disk());
  } else if (type == "default-quad-bulged") {
    set_quads(make_bulged_rect());
  } else if (type == "default-uvsphere") {
    set_quads(make_uvsphere());
  } else if (type == "default-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere());
  } else if (type == "default-uvdisk") {
    set_quads(make_uvdisk());
  } else if (type == "default-uvcylinder") {
    set_quads(make_uvcylinder());
  } else if (type == "default-uvcylinder-rounded") {
    set_quads(make_rounded_uvcylinder({32, 32, 32}));
  } else if (type == "default-geosphere") {
    set_triangles(make_geosphere());
  } else if (type == "default-floor") {
    set_quads(make_floor());
  } else if (type == "default-floor-bent") {
    set_quads(make_bent_floor());
  } else if (type == "default-matball") {
    set_quads(make_sphere());
  } else if (type == "default-hairball") {
    auto base = make_sphere(pow2(5), 0.8);
    set_lines(make_hair(base, {4, 65536}, {0.2, 0.2}, {0.002, 0.001}));
  } else if (type == "default-hairball-interior") {
    set_quads(make_sphere(pow2(5), 0.8));
  } else if (type == "default-suzanne") {
    set_quads(make_monkey());
  } else if (type == "default-cube-facevarying") {
    set_fvquads(make_fvbox());
  } else if (type == "default-sphere-facevarying") {
    set_fvquads(make_fvsphere());
  } else if (type == "default-quady-displaced") {
    set_quads(make_recty({256, 256}));
  } else if (type == "default-sphere-displaced") {
    set_quads(make_sphere(128));
  } else if (type == "test-cube") {
    set_quads(make_rounded_box(
        {32, 32, 32}, {0.075f, 0.075f, 0.075f}, {1, 1, 1}, 0.3 * 0.075f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere") {
    set_quads(make_uvsphere({32, 32}, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvsphere-flipcap") {
    set_quads(make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere") {
    set_quads(make_sphere(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-sphere-displaced") {
    set_quads(make_sphere(128, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-disk") {
    set_quads(make_disk(32, 0.075f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-uvcylinder") {
    set_quads(make_rounded_uvcylinder(
        {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-floor") {
    set_quads(make_floor({1, 1}, {2, 2}, {20, 20}));
  } else if (type == "test-quad") {
    set_quads(make_rect({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady") {
    set_quads(make_recty({1, 1}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quad-displaced") {
    set_quads(make_rect({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-quady-displaced") {
    set_quads(make_recty({256, 256}, {0.075, 0.075}, {1, 1}));
  } else if (type == "test-matball") {
    set_quads(make_sphere(32, 0.075));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-hairball1") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100}));
  } else if (type == "test-hairball2") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}));
  } else if (type == "test-hairball3") {
    auto base = make_sphere(32, 0.075f * 0.8f, 1);
    for (auto& p : base.positions) p += {0, 0.075, 0};
    set_lines(make_hair(base, {4, 65536}, {0.1f * 0.15f, 0.1f * 0.15f},
        {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128}));
  } else if (type == "test-hairball-interior") {
    set_quads(make_sphere(32, 0.075f * 0.8f, 1));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-suzanne-subdiv") {
    set_quads(make_monkey(0.075f * 0.8f));
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-cube-subdiv") {
    // set_quads(make_cube( 0.075f);
    set_fvquads(make_fvcube(0.075f));
    // make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals,
    //      texcoords, {1, 1, 1}, {0.075f, 0.075f, 0.075f});
    for (auto& p : shape.positions) p += {0, 0.075, 0};
  } else if (type == "test-arealight1") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-arealight2") {
    set_quads(make_rect({1, 1}, {0.2, 0.2}));
  } else if (type == "test-largearealight1") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-largearealight2") {
    set_quads(make_rect({1, 1}, {0.4, 0.4}));
  } else if (type == "test-pointlight1") {
    set_points(make_point(0));
  } else if (type == "test-pointlight2") {
    set_points(make_point(0));
  } else if (type == "test-point") {
    set_points(make_points(1));
  } else if (type == "test-points") {
    set_points(make_points(4096));
  } else if (type == "test-points-random") {
    set_points(make_random_points(4096, {0.2, 0.2, 0.2}));
  } else if (type == "test-particles") {
    set_points(make_points(4096));
  } else if (type == "test-cloth") {
    set_quads(make_rect({64, 64}, {0.2, 0.2}));
  } else if (type == "test-clothy") {
    set_quads(make_recty({64, 64}, {0.2, 0.2}));
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> shape_stats(const shape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("points:       " + format(shape.points.size()));
  stats.push_back("lines:        " + format(shape.lines.size()));
  stats.push_back("triangles:    " + format(shape.triangles.size()));
  stats.push_back("quads:        " + format(shape.quads.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("colors:       " + format(shape.colors.size()));
  stats.push_back("radius:       " + format(shape.radius.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

vector<string> fvshape_stats(const fvshape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("fvquads:      " + format(shape.quadspos.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMPUTATION OF PER-VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex tangents for lines.
vector<vec3f> lines_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto tangents = vector<vec3f>{positions.size()};
  for (auto& tangent : tangents) tangent = zero3f;
  for (auto& l : lines) {
    auto tangent = line_tangent(positions[l.x], positions[l.y]);
    auto length  = line_length(positions[l.x], positions[l.y]);
    tangents[l.x] += tangent * length;
    tangents[l.y] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
  return tangents;
}

// Compute per-vertex normals for triangles.
vector<vec3f> triangles_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = zero3f;
  for (auto& t : triangles) {
    auto normal = triangle_normal(
        positions[t.x], positions[t.y], positions[t.z]);
    auto area = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
    normals[t.x] += normal * area;
    normals[t.y] += normal * area;
    normals[t.z] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Compute per-vertex normals for quads.
vector<vec3f> quads_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = zero3f;
  for (auto& q : quads) {
    auto normal = quad_normal(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    auto area = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    normals[q.x] += normal * area;
    normals[q.y] += normal * area;
    normals[q.z] += normal * area;
    if (q.z != q.w) normals[q.w] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Compute per-vertex tangents for lines.
void lines_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions) {
  if (tangents.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& tangent : tangents) tangent = zero3f;
  for (auto& l : lines) {
    auto tangent = line_tangent(positions[l.x], positions[l.y]);
    auto length  = line_length(positions[l.x], positions[l.y]);
    tangents[l.x] += tangent * length;
    tangents[l.y] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
}

// Compute per-vertex normals for triangles.
void triangles_normals(vector<vec3f>& normals, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = zero3f;
  for (auto& t : triangles) {
    auto normal = triangle_normal(
        positions[t.x], positions[t.y], positions[t.z]);
    auto area = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
    normals[t.x] += normal * area;
    normals[t.y] += normal * area;
    normals[t.z] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex normals for quads.
void quads_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = zero3f;
  for (auto& q : quads) {
    auto normal = quad_normal(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    auto area = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    normals[q.x] += normal * area;
    normals[q.y] += normal * area;
    normals[q.z] += normal * area;
    if (q.z != q.w) normals[q.w] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> triangles_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  auto tangu = vector<vec3f>(positions.size(), zero3f);
  auto tangv = vector<vec3f>(positions.size(), zero3f);
  for (auto t : triangles) {
    auto tutv = triangle_tangents_fromuv(positions[t.x], positions[t.y],
        positions[t.z], texcoords[t.x], texcoords[t.y], texcoords[t.z]);
    for (auto vid : {t.x, t.y, t.z}) tangu[vid] += normalize(tutv.first);
    for (auto vid : {t.x, t.y, t.z}) tangv[vid] += normalize(tutv.second);
  }
  for (auto& t : tangu) t = normalize(t);
  for (auto& t : tangv) t = normalize(t);

  auto tangent_spaces = vector<vec4f>(positions.size());
  for (auto& tangent : tangent_spaces) tangent = zero4f;
  for (auto i = 0; i < positions.size(); i++) {
    tangu[i] = orthonormalize(tangu[i], normals[i]);
    auto s   = (dot(cross(normals[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
    tangent_spaces[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
  }
  return tangent_spaces;
}

// Apply skinning
pair<vector<vec3f>, vector<vec3f>> skin_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i = 0; i < positions.size(); i++) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
  }
  for (auto i = 0; i < normals.size(); i++) {
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning as specified in Khronos glTF
pair<vector<vec3f>, vector<vec3f>> skin_matrices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i = 0; i < positions.size(); i++) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning
void skin_vertices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i = 0; i < positions.size(); i++) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
  }
  for (auto i = 0; i < normals.size(); i++) {
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
}

// Apply skinning as specified in Khronos glTF
void skin_matrices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i = 0; i < positions.size(); i++) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  return lines_tangents(lines, positions);
}
vector<vec3f> compute_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  return triangles_normals(triangles, positions);
}
vector<vec3f> compute_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  return quads_normals(quads, positions);
}
// Update normals and tangents
void update_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions) {
  return lines_tangents(tangents, lines, positions);
}
void update_normals(vector<vec3f>& normals, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  return triangles_normals(normals, triangles, positions);
}
void update_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  return quads_normals(normals, quads, positions);
}

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  return triangles_tangent_spaces(triangles, positions, normals, texcoords);
}

// Apply skinning to vertex position and normals.
pair<vector<vec3f>, vector<vec3f>> compute_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms) {
  return skin_vertices(positions, normals, weights, joints, xforms);
}
// Apply skinning as specified in Khronos glTF.
pair<vector<vec3f>, vector<vec3f>> compute_matrix_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms) {
  return skin_matrices(positions, normals, weights, joints, xforms);
}
// Update skinning
void udpate_skinning(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms) {
  return skin_vertices(skinned_positions, skinned_normals, positions, normals,
      weights, joints, xforms);
}
void update_matrix_skinning(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms) {
  return skin_matrices(skinned_positions, skinned_normals, positions, normals,
      weights, joints, xforms);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
vector<vec3f> flip_normals(const vector<vec3f>& normals) {
  auto flipped = normals;
  for (auto& n : flipped) n = -n;
  return flipped;
}
// Flip face orientation
vector<vec3i> flip_triangles(const vector<vec3i>& triangles) {
  auto flipped = triangles;
  for (auto& t : flipped) swap(t.y, t.z);
  return flipped;
}
vector<vec4i> flip_quads(const vector<vec4i>& quads) {
  auto flipped = quads;
  for (auto& q : flipped) {
    if (q.z != q.w) {
      swap(q.y, q.w);
    } else {
      swap(q.y, q.z);
      q.w = q.z;
    }
  }
  return flipped;
}

// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
vector<vec3f> align_vertices(
    const vector<vec3f>& positions, const vec3i& alignment) {
  auto bounds = invalidb3f;
  for (auto& p : positions) bounds = merge(bounds, p);
  auto offset = vec3f{0, 0, 0};
  switch (alignment.x) {
    case 1: offset.x = bounds.min.x; break;
    case 2: offset.x = (bounds.min.x + bounds.max.x) / 2; break;
    case 3: offset.x = bounds.max.x; break;
  }
  switch (alignment.y) {
    case 1: offset.y = bounds.min.y; break;
    case 2: offset.y = (bounds.min.y + bounds.max.y) / 2; break;
    case 3: offset.y = bounds.max.y; break;
  }
  switch (alignment.z) {
    case 1: offset.z = bounds.min.z; break;
    case 2: offset.z = (bounds.min.z + bounds.max.z) / 2; break;
    case 3: offset.z = bounds.max.z; break;
  }
  auto aligned = positions;
  for (auto& p : aligned) p -= offset;
  return aligned;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EDGEA AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize an edge map with elements.
edge_map make_edge_map(const vector<vec3i>& triangles) {
  auto emap = edge_map{};
  for (auto& t : triangles) {
    insert_edge(emap, {t.x, t.y});
    insert_edge(emap, {t.y, t.z});
    insert_edge(emap, {t.z, t.x});
  }
  return emap;
}
edge_map make_edge_map(const vector<vec4i>& quads) {
  auto emap = edge_map{};
  for (auto& q : quads) {
    insert_edge(emap, {q.x, q.y});
    insert_edge(emap, {q.y, q.z});
    if (q.z != q.w) insert_edge(emap, {q.z, q.w});
    insert_edge(emap, {q.w, q.x});
  }
  return emap;
}
void insert_edges(edge_map& emap, const vector<vec3i>& triangles) {
  for (auto& t : triangles) {
    insert_edge(emap, {t.x, t.y});
    insert_edge(emap, {t.y, t.z});
    insert_edge(emap, {t.z, t.x});
  }
}
void insert_edges(edge_map& emap, const vector<vec4i>& quads) {
  for (auto& q : quads) {
    insert_edge(emap, {q.x, q.y});
    insert_edge(emap, {q.y, q.z});
    if (q.z != q.w) insert_edge(emap, {q.z, q.w});
    insert_edge(emap, {q.w, q.x});
  }
}
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge) {
  auto es = edge.x < edge.y ? edge : vec2i{edge.y, edge.x};
  auto it = emap.index.find(es);
  if (it == emap.index.end()) {
    auto idx = (int)emap.edges.size();
    emap.index.insert(it, {es, idx});
    emap.edges.push_back(es);
    emap.nfaces.push_back(1);
    return idx;
  } else {
    auto idx = it->second;
    emap.nfaces[idx] += 1;
    return idx;
  }
}
// Get number of edges
int num_edges(const edge_map& emap) { return (int)emap.edges.size(); }
// Get the edge index
int edge_index(const edge_map& emap, const vec2i& edge) {
  auto es       = edge.x < edge.y ? edge : vec2i{edge.y, edge.x};
  auto iterator = emap.index.find(es);
  if (iterator == emap.index.end()) return -1;
  return iterator->second;
}
// Get a list of edges, boundary edges, boundary vertices
vector<vec2i> get_edges(const edge_map& emap) { return emap.edges; }
vector<vec2i> get_boundary(const edge_map& emap) {
  auto boundary = vector<vec2i>{};
  for (auto idx = 0; idx < emap.edges.size(); idx++) {
    if (emap.nfaces[idx] < 2) boundary.push_back(emap.edges[idx]);
  }
  return boundary;
}
vector<vec2i> get_edges(const vector<vec3i>& triangles) {
  return get_edges(make_edge_map(triangles));
}
vector<vec2i> get_edges(const vector<vec4i>& quads) {
  return get_edges(make_edge_map(quads));
}
vector<vec2i> get_edges(
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
  auto edges      = get_edges(triangles);
  auto more_edges = get_edges(quads);
  edges.insert(edges.end(), more_edges.begin(), more_edges.end());
  return edges;
}

// Build adjacencies between faces (sorted counter-clockwise)
vector<vec3i> face_adjacencies(const vector<vec3i>& triangles) {
  auto get_edge = [](const vec3i& triangle, int i) -> vec2i {
    auto x = triangle[i], y = triangle[i < 2 ? i + 1 : 0];
    return x < y ? vec2i{x, y} : vec2i{y, x};
  };
  auto adjacencies = vector<vec3i>{triangles.size(), vec3i{-1, -1, -1}};
  auto edge_map    = unordered_map<vec2i, int>();
  edge_map.reserve((size_t)(triangles.size() * 1.5));
  for (int i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      auto edge = get_edge(triangles[i], k);
      auto it   = edge_map.find(edge);
      if (it == edge_map.end()) {
        edge_map.insert(it, {edge, i});
      } else {
        auto neighbor     = it->second;
        adjacencies[i][k] = neighbor;
        for (int kk = 0; kk < 3; ++kk) {
          auto edge2 = get_edge(triangles[neighbor], kk);
          if (edge2 == edge) {
            adjacencies[neighbor][kk] = i;
            break;
          }
        }
      }
    }
  }
  return adjacencies;
}

// Build adjacencies between vertices (sorted counter-clockwise)
vector<vector<int>> vertex_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies) {
  auto find_index = [](const vec3i& v, int x) {
    if (v.x == x) return 0;
    if (v.y == x) return 1;
    if (v.z == x) return 2;
    return -1;
  };

  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<int>(triangles.size() * 3, -1);

  for (int i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; k++) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<int>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (int i = 0; i < num_vertices; ++i) {
    result[i].reserve(6);
    auto first_face = face_from_vertex[i];
    if (first_face == -1) continue;

    auto face = first_face;
    while (true) {
      auto k = find_index(triangles[face], i);
      k      = k != 0 ? k - 1 : 2;
      result[i].push_back(triangles[face][k]);
      face = adjacencies[face][k];
      if (face == -1) break;
      if (face == first_face) break;
    }
  }

  return result;
}

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
vector<vector<int>> vertex_to_faces_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies) {
  auto find_index = [](const vec3i& v, int x) {
    if (v.x == x) return 0;
    if (v.y == x) return 1;
    if (v.z == x) return 2;
    return -1;
  };

  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<int>(triangles.size() * 3, -1);

  for (int i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; k++) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<int>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (int i = 0; i < num_vertices; ++i) {
    result[i].reserve(6);
    auto first_face = face_from_vertex[i];
    if (first_face == -1) continue;

    auto face = first_face;
    while (true) {
      auto k = find_index(triangles[face], i);
      k      = k != 0 ? k - 1 : 2;
      face   = adjacencies[face][k];
      result[i].push_back(face);
      if (face == -1) break;
      if (face == first_face) break;
    }
  }

  return result;
}

// Compute boundaries as a list of loops (sorted counter-clockwise)
vector<vector<int>> ordered_boundaries(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int num_vertices) {
  // map every boundary vertex to its next one
  auto next_vert = vector<int>(num_vertices, -1);
  for (int i = 0; i < triangles.size(); ++i) {
    for (int k = 0; k < 3; ++k) {
      if (adjacency[i][k] == -1)
        next_vert[triangles[i][k]] = triangles[i][(k + 1) % 3];
    }
  }

  // result
  auto boundaries = vector<vector<int>>();

  // arrange boundary vertices in loops
  for (int i = 0; i < next_vert.size(); i++) {
    if (next_vert[i] == -1) continue;

    // add new empty boundary
    boundaries.emplace_back();
    auto current = i;

    while (true) {
      auto next = next_vert[current];
      if (next == -1) {
        return {};
      }
      next_vert[current] = -1;
      boundaries.back().push_back(current);

      // close loop if necessary
      if (next == i)
        break;
      else
        current = next;
    }
  }

  return boundaries;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH
// -----------------------------------------------------------------------------
namespace yocto {

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-variable"
#ifndef __clang__
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#endif
#endif

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
static pair<int, int> split_sah(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto split_axis = 0;
  auto mid        = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, split_axis};

  // consider N bins, compute their cost and keep the minimum
  const int nbins    = 16;
  auto      middle   = 0.0f;
  auto      min_cost = flt_max;
  auto      area     = [](auto& b) {
    auto size = b.max - b.min;
    return 1e-12f + 2 * size.x * size.y + 2 * size.x * size.z +
           2 * size.y * size.z;
  };
  for (auto saxis = 0; saxis < 3; saxis++) {
    for (auto b = 1; b < nbins; b++) {
      auto split     = cbbox.min[saxis] + b * csize[saxis] / nbins;
      auto left_bbox = invalidb3f, right_bbox = invalidb3f;
      auto left_nprims = 0, right_nprims = 0;
      for (auto i = start; i < end; i++) {
        if (centers[primitives[i]][saxis] < split) {
          left_bbox = merge(left_bbox, bboxes[primitives[i]]);
          left_nprims += 1;
        } else {
          right_bbox = merge(right_bbox, bboxes[primitives[i]]);
          right_nprims += 1;
        }
      }
      auto cost = 1 + left_nprims * area(left_bbox) / area(cbbox) +
                  right_nprims * area(right_bbox) / area(cbbox);
      if (cost < min_cost) {
        min_cost   = cost;
        middle     = split;
        split_axis = saxis;
      }
    }
  }
  // split
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [split_axis, middle, &centers](
                      auto a) { return centers[a][split_axis] < middle; }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    split_axis = 0;
    mid        = (start + end) / 2;
    // throw std::runtime_error("bad bvh split");
  }

  return {mid, split_axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and
// axis.
static pair<int, int> split_balanced(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // balanced tree split: find the largest axis of the
  // bounding box and split along this one right in the middle
  mid = (start + end) / 2;
  std::nth_element(primitives.data() + start, primitives.data() + mid,
      primitives.data() + end, [axis, &centers](auto a, auto b) {
        return centers[a][axis] < centers[b][axis];
      });

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    axis = 0;
    mid  = (start + end) / 2;
    // throw std::runtime_error("bad bvh split");
  }

  return {mid, axis};
}

// Splits a BVH node using the middle heutirtic. Returns split position and
// axis.
static pair<int, int> split_middle(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // split the space in the middle along the largest axis
  auto cmiddle = (cbbox.max + cbbox.min) / 2;
  auto middle  = cmiddle[axis];
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [axis, middle, &centers](
                      auto a) { return centers[a][axis] < middle; }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    axis = 0;
    mid  = (start + end) / 2;
    // throw std::runtime_error("bad bvh split");
  }

  return {mid, axis};
}

#if !defined(_WIN32) && !defined(_WIN64)
#pragma GCC diagnostic pop
#endif

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh(shape_bvh& bvh, vector<bbox3f>& bboxes) {
  // get values
  auto& nodes      = bvh.nodes;
  auto& primitives = bvh.primitives;

  // prepare to build nodes
  nodes.clear();
  nodes.reserve(bboxes.size() * 2);

  // prepare primitives
  bvh.primitives.resize(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) bvh.primitives[idx] = idx;

  // prepare centers
  auto centers = vector<vec3f>(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++)
    centers[idx] = center(bboxes[idx]);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)bboxes.size()}};
  nodes.emplace_back();

  // create nodes until the queue is empty
  while (!queue.empty()) {
    // grab node to work on
    auto next = queue.front();
    queue.pop_front();
    auto nodeid = next.x, start = next.y, end = next.z;

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox = merge(node.bbox, bboxes[primitives[i]]);

    // split into two children
    if (end - start > bvh_max_prims) {
      // get split
      auto [mid, axis] = split_middle(primitives, bboxes, centers, start, end);

      // make an internal node
      node.internal = true;
      node.axis     = (int8_t)axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = (int16_t)(end - start);
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

// Update bvh
static void update_bvh(shape_bvh& bvh, const vector<bbox3f>& bboxes) {
  for (auto nodeid = (int)bvh.nodes.size() - 1; nodeid >= 0; nodeid--) {
    auto& node = bvh.nodes[nodeid];
    node.bbox  = invalidb3f;
    if (node.internal) {
      for (auto idx = 0; idx < 2; idx++) {
        node.bbox = merge(node.bbox, bvh.nodes[node.start + idx].bbox);
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        node.bbox = merge(node.bbox, bboxes[bvh.primitives[node.start + idx]]);
      }
    }
  }
}

// Build shape bvh
shape_bvh make_points_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // build nodes
  auto bvh = shape_bvh{};
  build_bvh(bvh, bboxes);
  return bvh;
}
shape_bvh make_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& l     = lines[idx];
    bboxes[idx] = line_bounds(
        positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
  }

  // build nodes
  auto bvh = shape_bvh{};
  build_bvh(bvh, bboxes);
  return bvh;
}
shape_bvh make_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& t     = triangles[idx];
    bboxes[idx] = triangle_bounds(
        positions[t.x], positions[t.y], positions[t.z]);
  }

  // build nodes
  auto bvh = shape_bvh{};
  build_bvh(bvh, bboxes);
  return bvh;
}
shape_bvh make_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& q     = quads[idx];
    bboxes[idx] = quad_bounds(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
  }

  // build nodes
  auto bvh = shape_bvh{};
  build_bvh(bvh, bboxes);
  return bvh;
}

void update_points_bvh(shape_bvh& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_lines_bvh(shape_bvh& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& l     = lines[idx];
    bboxes[idx] = line_bounds(
        positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_triangles_bvh(shape_bvh& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& t     = triangles[idx];
    bboxes[idx] = triangle_bounds(
        positions[t.x], positions[t.y], positions[t.z]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_quads_bvh(shape_bvh& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& q     = quads[idx];
    bboxes[idx] = quad_bounds(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}

// Intersect ray with a bvh.
template <typename Intersect>
static bool intersect_elements_bvh(const shape_bvh& bvh,
    Intersect&& intersect_element, const ray3f& ray_, int& element, vec2f& uv,
    float& distance, bool find_any) {
  // check empty
  if (bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 128>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto hit = false;

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        if (intersect_element(primitive, ray, uv, distance)) {
          hit      = true;
          element  = primitive;
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh.
shape_intersection intersect_points_bvh(const shape_bvh& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = intersect_elements_bvh(
      bvh,
      [&points, &positions, &radius](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& p = points[idx];
        return intersect_point(ray, positions[p], radius[p], uv, distance);
      },
      ray, intersection.element, intersection.uv, intersection.distance,
      find_any);
  return intersection;
}
shape_intersection intersect_lines_bvh(const shape_bvh& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = intersect_elements_bvh(
      bvh,
      [&lines, &positions, &radius](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& l = lines[idx];
        return intersect_line(ray, positions[l.x], positions[l.y], radius[l.x],
            radius[l.y], uv, distance);
      },
      ray, intersection.element, intersection.uv, intersection.distance,
      find_any);
  return intersection;
}
shape_intersection intersect_triangles_bvh(const shape_bvh& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = intersect_elements_bvh(
      bvh,
      [&triangles, &positions](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& t = triangles[idx];
        return intersect_triangle(
            ray, positions[t.x], positions[t.y], positions[t.z], uv, distance);
      },
      ray, intersection.element, intersection.uv, intersection.distance,
      find_any);
  return intersection;
}
shape_intersection intersect_quads_bvh(const shape_bvh& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = intersect_elements_bvh(
      bvh,
      [&quads, &positions](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& t = quads[idx];
        return intersect_quad(ray, positions[t.x], positions[t.y],
            positions[t.z], positions[t.w], uv, distance);
      },
      ray, intersection.element, intersection.uv, intersection.distance,
      find_any);
  return intersection;
}

// Intersect ray with a bvh.
template <typename Overlap>
static bool overlap_elements_bvh(const shape_bvh& bvh,
    Overlap&& overlap_element, const vec3f& pos, float max_distance,
    int& element, vec2f& uv, float& distance, bool find_any) {
  // check if empty
  if (bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 128>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        if (overlap_element(primitive, pos, max_distance, uv, distance)) {
          hit          = true;
          element      = primitive;
          max_distance = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
shape_intersection overlap_points_bvh(const shape_bvh& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = overlap_elements_bvh(
      bvh,
      [&points, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& p = points[idx];
        return overlap_point(
            pos, max_distance, positions[p], radius[p], uv, distance);
      },
      pos, max_distance, intersection.element, intersection.uv,
      intersection.distance, find_any);
  return intersection;
}
shape_intersection overlap_lines_bvh(const shape_bvh& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = overlap_elements_bvh(
      bvh,
      [&lines, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& l = lines[idx];
        return overlap_line(pos, max_distance, positions[l.x], positions[l.y],
            radius[l.x], radius[l.y], uv, distance);
      },
      pos, max_distance, intersection.element, intersection.uv,
      intersection.distance, find_any);
  return intersection;
}
shape_intersection overlap_triangles_bvh(const shape_bvh& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = overlap_elements_bvh(
      bvh,
      [&triangles, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& t = triangles[idx];
        return overlap_triangle(pos, max_distance, positions[t.x],
            positions[t.y], positions[t.z], radius[t.x], radius[t.y],
            radius[t.z], uv, distance);
      },
      pos, max_distance, intersection.element, intersection.uv,
      intersection.distance, find_any);
  return intersection;
}
shape_intersection overlap_quads_bvh(const shape_bvh& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  auto intersection = shape_intersection{};
  intersection.hit  = overlap_elements_bvh(
      bvh,
      [&quads, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& q = quads[idx];
        return overlap_quad(pos, max_distance, positions[q.x], positions[q.y],
            positions[q.z], positions[q.w], radius[q.x], radius[q.y],
            radius[q.z], radius[q.w], uv, distance);
      },
      pos, max_distance, intersection.element, intersection.uv,
      intersection.distance, find_any);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHTBORS
// -----------------------------------------------------------------------------

namespace yocto {

// Gets the cell index
vec3i get_cell_index(const hash_grid& grid, const vec3f& position) {
  auto scaledpos = position * grid.cell_inv_size;
  return vec3i{(int)scaledpos.x, (int)scaledpos.y, (int)scaledpos.z};
}

// Create a hash_grid
hash_grid make_hash_grid(float cell_size) {
  auto grid          = hash_grid{};
  grid.cell_size     = cell_size;
  grid.cell_inv_size = 1 / cell_size;
  return grid;
}
hash_grid make_hash_grid(const vector<vec3f>& positions, float cell_size) {
  auto grid          = hash_grid{};
  grid.cell_size     = cell_size;
  grid.cell_inv_size = 1 / cell_size;
  for (auto& position : positions) insert_vertex(grid, position);
  return grid;
}
// Inserts a point into the grid
int insert_vertex(hash_grid& grid, const vec3f& position) {
  auto vertex_id = (int)grid.positions.size();
  auto cell      = get_cell_index(grid, position);
  grid.cells[cell].push_back(vertex_id);
  grid.positions.push_back(position);
  return vertex_id;
}
// Finds the nearest neighbors within a given radius
void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    const vec3f& position, float max_radius, int skip_id) {
  auto cell        = get_cell_index(grid, position);
  auto cell_radius = (int)(max_radius * grid.cell_inv_size) + 1;
  neighbors.clear();
  auto max_radius_squared = max_radius * max_radius;
  for (auto k = -cell_radius; k <= cell_radius; k++) {
    for (auto j = -cell_radius; j <= cell_radius; j++) {
      for (auto i = -cell_radius; i <= cell_radius; i++) {
        auto ncell         = cell + vec3i{i, j, k};
        auto cell_iterator = grid.cells.find(ncell);
        if (cell_iterator == grid.cells.end()) continue;
        auto& ncell_vertices = cell_iterator->second;
        for (auto vertex_id : ncell_vertices) {
          if (distance_squared(grid.positions[vertex_id], position) >
              max_radius_squared)
            continue;
          if (vertex_id == skip_id) continue;
          neighbors.push_back(vertex_id);
        }
      }
    }
  }
}
void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    const vec3f& position, float max_radius) {
  find_neighbors(grid, neighbors, position, max_radius, -1);
}
void find_neighbors(const hash_grid& grid, vector<int>& neighbors, int vertex,
    float max_radius) {
  find_neighbors(grid, neighbors, grid.positions[vertex], max_radius, vertex);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
vector<vec3i> quads_to_triangles(const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  return triangles;
}

// Convert triangles to quads by creating degenerate quads
vector<vec4i> triangles_to_quads(const vector<vec3i>& triangles) {
  auto quads = vector<vec4i>{};
  quads.reserve(triangles.size());
  for (auto& t : triangles) quads.push_back({t.x, t.y, t.z, t.z});
  return quads;
}

// Convert beziers to lines using 3 lines for each bezier.
vector<vec2i> bezier_to_lines(const vector<vec4i>& beziers) {
  auto lines = vector<vec2i>{};
  lines.reserve(beziers.size() * 3);
  for (auto b : beziers) {
    lines.push_back({b.x, b.y});
    lines.push_back({b.y, b.z});
    lines.push_back({b.z, b.w});
  }
  return lines;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
split_facevarying(const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords) {
  auto split =
      std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>{};
  auto& [split_quads, split_positions, split_normals, split_texcoords] = split;
  // make faces unique
  unordered_map<vec3i, int> vert_map;
  split_quads.resize(quadspos.size());
  for (auto fid = 0; fid < quadspos.size(); fid++) {
    for (auto c = 0; c < 4; c++) {
      auto v = vec3i{
          (&quadspos[fid].x)[c],
          (!quadsnorm.empty()) ? (&quadsnorm[fid].x)[c] : -1,
          (!quadstexcoord.empty()) ? (&quadstexcoord[fid].x)[c] : -1,
      };
      auto it = vert_map.find(v);
      if (it == vert_map.end()) {
        auto s = (int)vert_map.size();
        vert_map.insert(it, {v, s});
        (&split_quads[fid].x)[c] = s;
      } else {
        (&split_quads[fid].x)[c] = it->second;
      }
    }
  }

  // fill vert data
  split_positions.clear();
  if (!positions.empty()) {
    split_positions.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_positions[index] = positions[vert.x];
    }
  }
  split_normals.clear();
  if (!normals.empty()) {
    split_normals.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_normals[index] = normals[vert.y];
    }
  }
  split_texcoords.clear();
  if (!texcoords.empty()) {
    split_texcoords.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_texcoords[index] = texcoords[vert.z];
    }
  }

  return split;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
void split_facevarying(vector<vec4i>& split_quads,
    vector<vec3f>& split_positions, vector<vec3f>& split_normals,
    vector<vec2f>& split_texcoords, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  // make faces unique
  unordered_map<vec3i, int> vert_map;
  split_quads.resize(quadspos.size());
  for (auto fid = 0; fid < quadspos.size(); fid++) {
    for (auto c = 0; c < 4; c++) {
      auto v = vec3i{
          (&quadspos[fid].x)[c],
          (!quadsnorm.empty()) ? (&quadsnorm[fid].x)[c] : -1,
          (!quadstexcoord.empty()) ? (&quadstexcoord[fid].x)[c] : -1,
      };
      auto it = vert_map.find(v);
      if (it == vert_map.end()) {
        auto s = (int)vert_map.size();
        vert_map.insert(it, {v, s});
        (&split_quads[fid].x)[c] = s;
      } else {
        (&split_quads[fid].x)[c] = it->second;
      }
    }
  }

  // fill vert data
  split_positions.clear();
  if (!positions.empty()) {
    split_positions.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_positions[index] = positions[vert.x];
    }
  }
  split_normals.clear();
  if (!normals.empty()) {
    split_normals.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_normals[index] = normals[vert.y];
    }
  }
  split_texcoords.clear();
  if (!texcoords.empty()) {
    split_texcoords.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_texcoords[index] = texcoords[vert.z];
    }
  }
}

// Split primitives per id
template <typename T>
vector<vector<T>> ungroup_elems_impl(
    const vector<T>& elems, const vector<int>& ids) {
  auto max_id      = *max_element(ids.begin(), ids.end());
  auto split_elems = vector<vector<T>>(max_id + 1);
  for (auto elem_id = 0; elem_id < elems.size(); elem_id++) {
    split_elems[ids[elem_id]].push_back(elems[elem_id]);
  }
  return split_elems;
}
vector<vector<vec2i>> ungroup_lines(
    const vector<vec2i>& lines, const vector<int>& ids) {
  return ungroup_elems_impl(lines, ids);
}
vector<vector<vec3i>> ungroup_triangles(
    const vector<vec3i>& triangles, const vector<int>& ids) {
  return ungroup_elems_impl(triangles, ids);
}
vector<vector<vec4i>> ungroup_quads(
    const vector<vec4i>& quads, const vector<int>& ids) {
  return ungroup_elems_impl(quads, ids);
}

// Weld vertices within a threshold.
pair<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold) {
  auto indices   = vector<int>(positions.size());
  auto welded    = vector<vec3f>{};
  auto grid      = make_hash_grid(threshold);
  auto neighbors = vector<int>{};
  for (auto vertex = 0; vertex < positions.size(); vertex++) {
    auto& position = positions[vertex];
    find_neighbors(grid, neighbors, position, threshold);
    if (neighbors.empty()) {
      welded.push_back(position);
      indices[vertex] = (int)welded.size() - 1;
      insert_vertex(grid, position);
    } else {
      indices[vertex] = neighbors.front();
    }
  }
  return {welded, indices};
}
pair<vector<vec3i>, vector<vec3f>> weld_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    float threshold) {
  auto [wpositions, indices] = weld_vertices(positions, threshold);
  auto wtriangles            = triangles;
  for (auto& t : wtriangles) t = {indices[t.x], indices[t.y], indices[t.z]};
  return {wtriangles, wpositions};
}
pair<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold) {
  auto [wpositions, indices] = weld_vertices(positions, threshold);
  auto wquads                = quads;
  for (auto& q : wquads)
    q = {
        indices[q.x],
        indices[q.y],
        indices[q.z],
        indices[q.w],
    };
  return {wquads, wpositions};
}

// Merge shape elements
void merge_lines(
    vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts) {
  for (auto& l : merge_lines)
    lines.push_back({l.x + num_verts, l.y + num_verts});
}
void merge_triangles(vector<vec3i>& triangles,
    const vector<vec3i>& merge_triangles, int num_verts) {
  for (auto& t : merge_triangles)
    triangles.push_back({t.x + num_verts, t.y + num_verts, t.z + num_verts});
}
void merge_quads(
    vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts) {
  for (auto& q : merge_quads)
    quads.push_back(
        {q.x + num_verts, q.y + num_verts, q.z + num_verts, q.w + num_verts});
}
void merge_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& tangents, vector<vec2f>& texcoords, vector<float>& radius,
    const vector<vec2i>& merge_lines, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_tangents,
    const vector<vec2f>& merge_texturecoords,
    const vector<float>& merge_radius) {
  auto merge_verts = (int)positions.size();
  for (auto& l : merge_lines)
    lines.push_back({l.x + merge_verts, l.y + merge_verts});
  positions.insert(
      positions.end(), merge_positions.begin(), merge_positions.end());
  tangents.insert(tangents.end(), merge_tangents.begin(), merge_tangents.end());
  texcoords.insert(
      texcoords.end(), merge_texturecoords.begin(), merge_texturecoords.end());
  radius.insert(radius.end(), merge_radius.begin(), merge_radius.end());
}
void merge_triangles(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vector<vec3i>& merge_triangles, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals,
    const vector<vec2f>& merge_texturecoords) {
  auto merge_verts = (int)positions.size();
  for (auto& t : merge_triangles)
    triangles.push_back(
        {t.x + merge_verts, t.y + merge_verts, t.z + merge_verts});
  positions.insert(
      positions.end(), merge_positions.begin(), merge_positions.end());
  normals.insert(normals.end(), merge_normals.begin(), merge_normals.end());
  texcoords.insert(
      texcoords.end(), merge_texturecoords.begin(), merge_texturecoords.end());
}
void merge_quads(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vector<vec4i>& merge_quads, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals,
    const vector<vec2f>& merge_texturecoords) {
  auto merge_verts = (int)positions.size();
  for (auto& q : merge_quads)
    quads.push_back({q.x + merge_verts, q.y + merge_verts, q.z + merge_verts,
        q.w + merge_verts});
  positions.insert(
      positions.end(), merge_positions.begin(), merge_positions.end());
  normals.insert(normals.end(), merge_normals.begin(), merge_normals.end());
  texcoords.insert(
      texcoords.end(), merge_texturecoords.begin(), merge_texturecoords.end());
}

void merge_triangles_and_quads(
    vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles) {
  if (quads.empty()) return;
  if (force_triangles) {
    auto qtriangles = quads_to_triangles(quads);
    triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
    quads = {};
  } else {
    auto tquads = triangles_to_quads(triangles);
    quads.insert(quads.end(), tquads.begin(), tquads.end());
    triangles = {};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines.
template <typename T>
void subdivide_lines_impl(vector<vec2i>& lines, vector<T>& vert,
    const vector<vec2i>& lines_, const vector<T>& vert_, int level) {
  // initialization
  lines = lines_;
  vert  = vert_;
  // early exit
  if (lines.empty() || vert.empty()) return;
  // loop over levels
  for (auto l = 0; l < level; l++) {
    // sizes
    auto nverts = (int)vert.size();
    auto nlines = (int)lines.size();
    // create vertices
    auto tvert = vector<T>(nverts + nlines);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nlines; i++) {
      auto l            = lines[i];
      tvert[nverts + i] = (vert[l.x] + vert[l.y]) / 2;
    }
    // create lines
    auto tlines = vector<vec2i>(nlines * 2);
    for (auto i = 0; i < nlines; i++) {
      auto l            = lines[i];
      tlines[i * 2 + 0] = {l.x, nverts + i};
      tlines[i * 2 + 0] = {nverts + i, l.y};
    }
    swap(tlines, lines);
    swap(tvert, vert);
  }
}
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_lines_impl(
    const vector<vec2i>& lines, const vector<T>& vert, int level) {
  auto tess = pair<vector<vec2i>, vector<T>>{};
  subdivide_lines_impl(tess.first, tess.second, lines, vert, level);
  return tess;
}

// Subdivide triangle.
template <typename T>
void subdivide_triangles_impl(vector<vec3i>& triangles, vector<T>& vert,
    const vector<vec3i>& triangles_, const vector<T>& vert_, int level) {
  // initialization
  triangles = triangles_;
  vert      = vert_;
  // early exit
  if (triangles.empty() || vert.empty()) return;
  // loop over levels
  for (auto l = 0; l < level; l++) {
    // get edges
    auto emap  = make_edge_map(triangles);
    auto edges = get_edges(emap);
    // number of elements
    auto nverts = (int)vert.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)triangles.size();
    // create vertices
    auto tvert = vector<T>(nverts + nedges);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
      auto e            = edges[i];
      tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    // create triangles
    auto ttriangles = vector<vec3i>(nfaces * 4);
    for (auto i = 0; i < nfaces; i++) {
      auto t                = triangles[i];
      ttriangles[i * 4 + 0] = {t.x, nverts + edge_index(emap, {t.x, t.y}),
          nverts + edge_index(emap, {t.z, t.x})};
      ttriangles[i * 4 + 1] = {t.y, nverts + edge_index(emap, {t.y, t.z}),
          nverts + edge_index(emap, {t.x, t.y})};
      ttriangles[i * 4 + 2] = {t.z, nverts + edge_index(emap, {t.z, t.x}),
          nverts + edge_index(emap, {t.y, t.z})};
      ttriangles[i * 4 + 3] = {nverts + edge_index(emap, {t.x, t.y}),
          nverts + edge_index(emap, {t.y, t.z}),
          nverts + edge_index(emap, {t.z, t.x})};
    }
    swap(ttriangles, triangles);
    swap(tvert, vert);
  }
}
template <typename T>
pair<vector<vec3i>, vector<T>> subdivide_triangles_impl(
    const vector<vec3i>& triangles, const vector<T>& vert, int level) {
  auto tess = pair<vector<vec3i>, vector<T>>{};
  subdivide_triangles_impl(tess.first, tess.second, triangles, vert, level);
  return tess;
}

// Subdivide quads.
template <typename T>
void subdivide_quads_impl(vector<vec4i>& quads, vector<T>& vert,
    const vector<vec4i>& quads_, const vector<T>& vert_, int level) {
  // initialization
  quads = quads_;
  vert  = vert_;
  // early exit
  if (quads.empty() || vert.empty()) return;
  // loop over levels
  for (auto l = 0; l < level; l++) {
    // get edges
    auto emap  = make_edge_map(quads);
    auto edges = get_edges(emap);
    // number of elements
    auto nverts = (int)vert.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();
    // create vertices
    auto tvert = vector<T>(nverts + nedges + nfaces);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
      auto e            = edges[i];
      tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    for (auto i = 0; i < nfaces; i++) {
      auto q = quads[i];
      if (q.z != q.w) {
        tvert[nverts + nedges + i] =
            (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4;
      } else {
        tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.z]) / 3;
      }
    }
    // create quads
    auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
    auto qi     = 0;
    for (auto i = 0; i < nfaces; i++) {
      auto q = quads[i];
      if (q.z != q.w) {
        tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
            nverts + nedges + i, nverts + edge_index(emap, {q.w, q.x})};
        tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
            nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.w}),
            nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
        tquads[qi++] = {q.w, nverts + edge_index(emap, {q.w, q.x}),
            nverts + nedges + i, nverts + edge_index(emap, {q.z, q.w})};
      } else {
        tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
            nverts + nedges + i, nverts + edge_index(emap, {q.z, q.x})};
        tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
            nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.x}),
            nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
      }
    }
    tquads.resize(qi);
    // done
    swap(tquads, quads);
    swap(tvert, vert);
  }
}
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_quads_impl(
    const vector<vec4i>& quads, const vector<T>& vert, int level) {
  auto tess = pair<vector<vec4i>, vector<T>>{};
  subdivide_quads_impl(tess.first, tess.second, quads, vert, level);
  return tess;
}

// Subdivide beziers.
template <typename T>
void subdivide_beziers_impl(vector<vec4i>& beziers, vector<T>& vert,
    const vector<vec4i>& beziers_, const vector<T>& vert_, int level) {
  // initialization
  beziers = beziers_;
  vert    = vert_;
  // early exit
  if (beziers.empty() || vert.empty()) return;
  // loop over levels
  for (auto l = 0; l < level; l++) {
    // get edges
    auto vmap     = unordered_map<int, int>();
    auto tvert    = vector<T>();
    auto tbeziers = vector<vec4i>();
    for (auto b : beziers) {
      if (vmap.find(b.x) == vmap.end()) {
        vmap[b.x] = (int)tvert.size();
        tvert.push_back(vert[b.x]);
      }
      if (vmap.find(b.w) == vmap.end()) {
        vmap[b.w] = (int)tvert.size();
        tvert.push_back(vert[b.w]);
      }
      auto bo = (int)tvert.size();
      tbeziers.push_back({vmap.at(b.x), bo + 0, bo + 1, bo + 2});
      tbeziers.push_back({bo + 2, bo + 3, bo + 4, vmap.at(b.w)});
      tvert.push_back(vert[b.x] / 2 + vert[b.y] / 2);
      tvert.push_back(vert[b.x] / 4 + vert[b.y] / 2 + vert[b.z] / 4);
      tvert.push_back(vert[b.x] / 8 + vert[b.y] * ((float)3 / (float)8) +
                      vert[b.z] * ((float)3 / (float)8) + vert[b.w] / 8);
      tvert.push_back(vert[b.y] / 4 + vert[b.z] / 2 + vert[b.w] / 4);
      tvert.push_back(vert[b.z] / 2 + vert[b.w] / 2);
    }

    // done
    swap(tbeziers, beziers);
    swap(tvert, vert);
  }
}
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_beziers_impl(
    const vector<vec4i>& beziers, const vector<T>& vert, int level) {
  auto tess = pair<vector<vec4i>, vector<T>>{};
  subdivide_beziers_impl(tess.first, tess.second, beziers, vert, level);
  return tess;
}

// Subdivide catmullclark.
template <typename T>
void subdivide_catmullclark_impl(vector<vec4i>& quads, vector<T>& vert,
    const vector<vec4i>& quads_, const vector<T>& vert_, int level,
    bool lock_boundary) {
  // initialization
  quads = quads_;
  vert  = vert_;
  // early exit
  if (quads.empty() || vert.empty()) return;
  // loop over levels
  for (auto l = 0; l < level; l++) {
    // get edges
    auto emap     = make_edge_map(quads);
    auto edges    = get_edges(emap);
    auto boundary = get_boundary(emap);
    // number of elements
    auto nverts    = (int)vert.size();
    auto nedges    = (int)edges.size();
    auto nboundary = (int)boundary.size();
    auto nfaces    = (int)quads.size();

    // split elements ------------------------------------
    // create vertices
    auto tvert = vector<T>(nverts + nedges + nfaces);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
      auto e            = edges[i];
      tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    for (auto i = 0; i < nfaces; i++) {
      auto q = quads[i];
      if (q.z != q.w) {
        tvert[nverts + nedges + i] =
            (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4;
      } else {
        tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.z]) / 3;
      }
    }
    // create quads
    auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
    auto qi     = 0;
    for (auto i = 0; i < nfaces; i++) {
      auto q = quads[i];
      if (q.z != q.w) {
        tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
            nverts + nedges + i, nverts + edge_index(emap, {q.w, q.x})};
        tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
            nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.w}),
            nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
        tquads[qi++] = {q.w, nverts + edge_index(emap, {q.w, q.x}),
            nverts + nedges + i, nverts + edge_index(emap, {q.z, q.w})};
      } else {
        tquads[qi++] = {q.x, nverts + edge_index(emap, {q.x, q.y}),
            nverts + nedges + i, nverts + edge_index(emap, {q.z, q.x})};
        tquads[qi++] = {q.y, nverts + edge_index(emap, {q.y, q.z}),
            nverts + nedges + i, nverts + edge_index(emap, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + edge_index(emap, {q.z, q.x}),
            nverts + nedges + i, nverts + edge_index(emap, {q.y, q.z})};
      }
    }
    tquads.resize(qi);

    // split boundary
    auto tboundary = vector<vec2i>(nboundary * 2);
    for (auto i = 0; i < nboundary; i++) {
      auto e               = boundary[i];
      tboundary[i * 2 + 0] = {e.x, nverts + edge_index(emap, e)};
      tboundary[i * 2 + 1] = {nverts + edge_index(emap, e), e.y};
    }

    // setup creases -----------------------------------
    auto tcrease_edges = vector<vec2i>();
    auto tcrease_verts = vector<int>();
    if (lock_boundary) {
      for (auto& b : tboundary) {
        tcrease_verts.push_back(b.x);
        tcrease_verts.push_back(b.y);
      }
    } else {
      for (auto& b : tboundary) tcrease_edges.push_back(b);
    }

    // define vertex valence ---------------------------
    auto tvert_val = vector<int>(tvert.size(), 2);
    for (auto& e : tboundary) {
      tvert_val[e.x] = (lock_boundary) ? 0 : 1;
      tvert_val[e.y] = (lock_boundary) ? 0 : 1;
    }

    // averaging pass ----------------------------------
    auto avert  = vector<T>(tvert.size(), T());
    auto acount = vector<int>(tvert.size(), 0);
    for (auto p : tcrease_verts) {
      if (tvert_val[p] != 0) continue;
      avert[p] += tvert[p];
      acount[p] += 1;
    }
    for (auto& e : tcrease_edges) {
      auto c = (tvert[e.x] + tvert[e.y]) / 2;
      for (auto vid : {e.x, e.y}) {
        if (tvert_val[vid] != 1) continue;
        avert[vid] += c;
        acount[vid] += 1;
      }
    }
    for (auto& q : tquads) {
      auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4;
      for (auto vid : {q.x, q.y, q.z, q.w}) {
        if (tvert_val[vid] != 2) continue;
        avert[vid] += c;
        acount[vid] += 1;
      }
    }
    for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i = 0; i < tvert.size(); i++) {
      if (tvert_val[i] != 2) continue;
      avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4 / (float)acount[i]);
    }
    tvert = avert;

    // done
    swap(tquads, quads);
    swap(tvert, vert);
  }
}
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_catmullclark_impl(
    const vector<vec4i>& quads, const vector<T>& vert, int level,
    bool lock_boundary) {
  auto tess = pair<vector<vec4i>, vector<T>>{};
  subdivide_catmullclark_impl(
      tess.first, tess.second, quads, vert, level, lock_boundary);
  return tess;
}

pair<vector<vec2i>, vector<float>> subdivide_lines(
    const vector<vec2i>& lines, const vector<float>& vert, int level) {
  return subdivide_lines_impl(lines, vert, level);
}
pair<vector<vec2i>, vector<vec2f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec2f>& vert, int level) {
  return subdivide_lines_impl(lines, vert, level);
}
pair<vector<vec2i>, vector<vec3f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec3f>& vert, int level) {
  return subdivide_lines_impl(lines, vert, level);
}
pair<vector<vec2i>, vector<vec4f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec4f>& vert, int level) {
  return subdivide_lines_impl(lines, vert, level);
}

pair<vector<vec3i>, vector<float>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<float>& vert, int level) {
  return subdivide_triangles_impl(triangles, vert, level);
}
pair<vector<vec3i>, vector<vec2f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec2f>& vert, int level) {
  return subdivide_triangles_impl(triangles, vert, level);
}
pair<vector<vec3i>, vector<vec3f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& vert, int level) {
  return subdivide_triangles_impl(triangles, vert, level);
}
pair<vector<vec3i>, vector<vec4f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec4f>& vert, int level) {
  return subdivide_triangles_impl(triangles, vert, level);
}

pair<vector<vec4i>, vector<float>> subdivide_quads(
    const vector<vec4i>& quads, const vector<float>& vert, int level) {
  return subdivide_quads_impl(quads, vert, level);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec2f>& vert, int level) {
  return subdivide_quads_impl(quads, vert, level);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec3f>& vert, int level) {
  return subdivide_quads_impl(quads, vert, level);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec4f>& vert, int level) {
  return subdivide_quads_impl(quads, vert, level);
}

pair<vector<vec4i>, vector<float>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<float>& vert, int level) {
  return subdivide_beziers_impl(beziers, vert, level);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec2f>& vert, int level) {
  return subdivide_beziers_impl(beziers, vert, level);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec3f>& vert, int level) {
  return subdivide_beziers_impl(beziers, vert, level);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec4f>& vert, int level) {
  return subdivide_beziers_impl(beziers, vert, level);
}

pair<vector<vec4i>, vector<float>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<float>& vert, int level,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vert, level, lock_boundary);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec2f>& vert, int level,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vert, level, lock_boundary);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec3f>& vert, int level,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vert, level, lock_boundary);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec4f>& vert, int level,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vert, level, lock_boundary);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
int sample_points(int npoints, float re) { return sample_uniform(npoints, re); }
int sample_points(const vector<float>& cdf, float re) {
  return sample_discrete_cdf(cdf, re);
}
vector<float> sample_points_cdf(int npoints) {
  auto cdf = vector<float>(npoints);
  for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i != 0 ? cdf[i - 1] : 0);
  return cdf;
}
void sample_points_cdf(vector<float>& cdf, int npoints) {
  for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i != 0 ? cdf[i - 1] : 0);
}

// Pick a point on lines uniformly.
pair<int, float> sample_lines(const vector<float>& cdf, float re, float ru) {
  return {sample_discrete_cdf(cdf, re), ru};
}
vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto cdf = vector<float>(lines.size());
  for (auto i = 0; i < cdf.size(); i++) {
    auto& l = lines[i];
    auto  w = line_length(positions[l.x], positions[l.y]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_lines_cdf(vector<float>& cdf, const vector<vec2i>& lines,
    const vector<vec3f>& positions) {
  for (auto i = 0; i < cdf.size(); i++) {
    auto& l = lines[i];
    auto  w = line_length(positions[l.x], positions[l.y]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
}

// Pick a point on a triangle mesh uniformly.
pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv) {
  return {sample_discrete_cdf(cdf, re), sample_triangle(ruv)};
}
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto cdf = vector<float>(triangles.size());
  for (auto i = 0; i < cdf.size(); i++) {
    auto& t = triangles[i];
    auto  w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_triangles_cdf(vector<float>& cdf, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  for (auto i = 0; i < cdf.size(); i++) {
    auto& t = triangles[i];
    auto  w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
}

// Pick a point on a quad mesh uniformly.
pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv) {
  return {sample_discrete_cdf(cdf, re), ruv};
}
pair<int, vec2f> sample_quads(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv) {
  auto element = sample_discrete_cdf(cdf, re);
  if (quads[element].z == quads[element].w) {
    return {element, sample_triangle(ruv)};
  } else {
    return {element, ruv};
  }
}
vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto cdf = vector<float>(quads.size());
  for (auto i = 0; i < cdf.size(); i++) {
    auto& q = quads[i];
    auto  w = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    cdf[i] = w + (i ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_quads_cdf(vector<float>& cdf, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  for (auto i = 0; i < cdf.size(); i++) {
    auto& q = quads[i];
    auto  w = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    cdf[i] = w + (i ? cdf[i - 1] : 0);
  }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
void sample_triangles(vector<vec3f>& sampled_positions,
    vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
    int seed) {
  sampled_positions.resize(npoints);
  sampled_normals.resize(npoints);
  sampled_texcoords.resize(npoints);
  auto cdf = sample_triangles_cdf(triangles, positions);
  auto rng = make_rng(seed);
  for (auto i = 0; i < npoints; i++) {
    auto  sample         = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto& t              = triangles[sample.first];
    auto  uv             = sample.second;
    sampled_positions[i] = interpolate_triangle(
        positions[t.x], positions[t.y], positions[t.z], uv);
    if (!sampled_normals.empty()) {
      sampled_normals[i] = normalize(
          interpolate_triangle(normals[t.x], normals[t.y], normals[t.z], uv));
    } else {
      sampled_normals[i] = triangle_normal(
          positions[t.x], positions[t.y], positions[t.z]);
    }
    if (!sampled_texcoords.empty()) {
      sampled_texcoords[i] = interpolate_triangle(
          texcoords[t.x], texcoords[t.y], texcoords[t.z], uv);
    } else {
      sampled_texcoords[i] = zero2f;
    }
  }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
void sample_quads(vector<vec3f>& sampled_positions,
    vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
    int seed) {
  sampled_positions.resize(npoints);
  sampled_normals.resize(npoints);
  sampled_texcoords.resize(npoints);
  auto cdf = sample_quads_cdf(quads, positions);
  auto rng = make_rng(seed);
  for (auto i = 0; i < npoints; i++) {
    auto  sample         = sample_quads(cdf, rand1f(rng), rand2f(rng));
    auto& q              = quads[sample.first];
    auto  uv             = sample.second;
    sampled_positions[i] = interpolate_quad(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w], uv);
    if (!sampled_normals.empty()) {
      sampled_normals[i] = normalize(interpolate_quad(
          normals[q.x], normals[q.y], normals[q.z], normals[q.w], uv));
    } else {
      sampled_normals[i] = quad_normal(
          positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    }
    if (!sampled_texcoords.empty()) {
      sampled_texcoords[i] = interpolate_quad(
          texcoords[q.x], texcoords[q.y], texcoords[q.z], texcoords[q.w], uv);
    } else {
      sampled_texcoords[i] = zero2f;
    }
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a quad.
void make_rect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale) {
  positions.resize((steps.x + 1) * (steps.y + 1));
  normals.resize((steps.x + 1) * (steps.y + 1));
  texcoords.resize((steps.x + 1) * (steps.y + 1));
  for (auto j = 0; j <= steps.y; j++) {
    for (auto i = 0; i <= steps.x; i++) {
      auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
      positions[j * (steps.x + 1) + i] = {
          (2 * uv.x - 1) * scale.x, (2 * uv.y - 1) * scale.y, 0};
      normals[j * (steps.x + 1) + i]   = {0, 0, 1};
      texcoords[j * (steps.x + 1) + i] = vec2f{uv.x, 1 - uv.y} * uvscale;
    }
  }

  quads.resize(steps.x * steps.y);
  for (auto j = 0; j < steps.y; j++) {
    for (auto i = 0; i < steps.x; i++) {
      quads[j * steps.x + i] = {j * (steps.x + 1) + i,
          j * (steps.x + 1) + i + 1, (j + 1) * (steps.x + 1) + i + 1,
          (j + 1) * (steps.x + 1) + i};
    }
  }
}

void make_bulged_rect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale, float height) {
  make_rect(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (height != 0) {
    height      = min(height, min(scale));
    auto radius = (1 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i = 0; i < positions.size(); i++) {
      auto pn      = normalize(positions[i] - center);
      positions[i] = center + pn * radius;
      normals[i]   = pn;
    }
  }
}

// Make a quad.
void make_recty(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale) {
  make_rect(quads, positions, normals, texcoords, steps, scale, uvscale);
  for (auto& p : positions) {
    std::swap(p.y, p.z);
    p.z = -p.z;
  }
  for (auto& n : normals) std::swap(n.y, n.z);
}

void make_bulged_recty(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale, float height) {
  make_bulged_rect(
      quads, positions, normals, texcoords, steps, scale, uvscale, height);
  for (auto& p : positions) {
    std::swap(p.y, p.z);
    p.z = -p.z;
  }
  for (auto& n : normals) std::swap(n.y, n.z);
}

// Make a cube.
void make_box(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec3f& scale, const vec3f& uvscale) {
  quads.clear();
  positions.clear();
  normals.clear();
  texcoords.clear();
  auto qquads         = vector<vec4i>{};
  auto qpositions     = vector<vec3f>{};
  auto qnormals       = vector<vec3f>{};
  auto qtexturecoords = vector<vec2f>{};
  // + z
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.x, steps.y},
      {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qpositions) p = {p.x, p.y, scale.z};
  for (auto& n : qnormals) n = {0, 0, 1};
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
  // - z
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.x, steps.y},
      {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qpositions) p = {-p.x, p.y, -scale.z};
  for (auto& n : qnormals) n = {0, 0, -1};
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
  // + x
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.z, steps.y},
      {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qpositions) p = {scale.x, p.y, -p.x};
  for (auto& n : qnormals) n = {1, 0, 0};
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
  // - x
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.z, steps.y},
      {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qpositions) p = {-scale.x, p.y, p.x};
  for (auto& n : qnormals) n = {-1, 0, 0};
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
  // + y
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.x, steps.z},
      {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i = 0; i < qpositions.size(); i++) {
    qpositions[i] = {qpositions[i].x, scale.y, -qpositions[i].y};
    qnormals[i]   = {0, 1, 0};
  }
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
  // - y
  make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.x, steps.z},
      {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i = 0; i < qpositions.size(); i++) {
    qpositions[i] = {qpositions[i].x, -scale.y, qpositions[i].y};
    qnormals[i]   = {0, -1, 0};
  }
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexturecoords);
}

void make_rounded_box(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec3f& scale, const vec3f& uvscale, float radius) {
  make_box(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (radius != 0) {
    radius = min(radius, min(scale));
    auto c = scale - radius;
    for (auto i = 0; i < positions.size(); i++) {
      auto pc = vec3f{
          abs(positions[i].x), abs(positions[i].y), abs(positions[i].z)};
      auto ps = vec3f{positions[i].x < 0 ? -1.0f : 1.0f,
          positions[i].y < 0 ? -1.0f : 1.0f, positions[i].z < 0 ? -1.0f : 1.0f};
      if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
        auto pn      = normalize(pc - c);
        positions[i] = c + radius * pn;
        normals[i]   = pn;
      } else if (pc.x >= c.x && pc.y >= c.y) {
        auto pn      = normalize((pc - c) * vec3f{1, 1, 0});
        positions[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
        normals[i]   = pn;
      } else if (pc.x >= c.x && pc.z >= c.z) {
        auto pn      = normalize((pc - c) * vec3f{1, 0, 1});
        positions[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
        normals[i]   = pn;
      } else if (pc.y >= c.y && pc.z >= c.z) {
        auto pn      = normalize((pc - c) * vec3f{0, 1, 1});
        positions[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
        normals[i]   = pn;
      } else {
        continue;
      }
      positions[i] *= ps;
      normals[i] *= ps;
    }
  }
}

void make_rect_stack(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec3f& scale, const vec2f& uvscale) {
  auto qquads         = vector<vec4i>{};
  auto qpositions     = vector<vec3f>{};
  auto qnormals       = vector<vec3f>{};
  auto qtexturecoords = vector<vec2f>{};
  for (auto i = 0; i <= steps.z; i++) {
    make_rect(qquads, qpositions, qnormals, qtexturecoords, {steps.x, steps.y},
        {scale.x, scale.y}, uvscale);
    for (auto& p : qpositions) p.z = (-1 + 2 * (float)i / steps.z) * scale.z;
    merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
        qnormals, qtexturecoords);
  }
}

void make_floor(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale) {
  make_rect(quads, positions, normals, texcoords, steps, scale, uvscale);
  for (auto& p : positions) {
    std::swap(p.y, p.z);
    p.z = -p.z;
  }
  for (auto& n : normals) std::swap(n.y, n.z);
}

void make_bent_floor(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& scale, const vec2f& uvscale, float radius) {
  make_floor(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (radius != 0) {
    radius     = min(radius, scale.y);
    auto start = (scale.y - radius) / 2;
    auto end   = start + radius;
    for (auto i = 0; i < positions.size(); i++) {
      if (positions[i].z < -end) {
        positions[i] = {positions[i].x, -positions[i].z - end + radius, -end};
        normals[i]   = {0, 0, 1};
      } else if (positions[i].z < -start && positions[i].z >= -end) {
        auto phi     = (pif / 2) * (-positions[i].z - start) / radius;
        positions[i] = {positions[i].x, -cos(phi) * radius + radius,
            -sin(phi) * radius - start};
        normals[i]   = {0, cos(phi), sin(phi)};
      } else {
      }
    }
  }
}

// Generate a sphere
void make_sphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps, float scale,
    float uvscale) {
  make_box(quads, positions, normals, texcoords, {steps, steps, steps},
      {scale, scale, scale}, {uvscale, uvscale, uvscale});
  for (auto& p : positions) p = normalize(p) * scale;
  normals = positions;
  for (auto& n : normals) n = normalize(n);
}

// Generate a uvsphere
void make_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float scale, const vec2f& uvscale) {
  make_rect(quads, positions, normals, texcoords, steps, {1, 1}, {1, 1});
  for (auto i = 0; i < positions.size(); i++) {
    auto uv      = texcoords[i];
    auto a       = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
    positions[i] = vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)} *
                   scale;
    normals[i]   = normalize(positions[i]);
    texcoords[i] = uv * uvscale;
  }
}

void make_capped_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float scale, const vec2f& uvscale, float cap) {
  make_uvsphere(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (cap != 0) {
    cap        = min(cap, scale / 2);
    auto zflip = (scale - cap);
    for (auto i = 0; i < positions.size(); i++) {
      if (positions[i].z > zflip) {
        positions[i].z = 2 * zflip - positions[i].z;
        normals[i].x   = -normals[i].x;
        normals[i].y   = -normals[i].y;
      } else if (positions[i].z < -zflip) {
        positions[i].z = 2 * (-zflip) - positions[i].z;
        normals[i].x   = -normals[i].x;
        normals[i].y   = -normals[i].y;
      }
    }
  }
}

// Generate a disk
void make_disk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps, float scale,
    float uvscale) {
  make_rect(quads, positions, normals, texcoords, {steps, steps}, {1, 1},
      {uvscale, uvscale});
  for (auto& position : positions) {
    // Analytical Methods for Squaring the Disc, by C. Fong
    // https://arxiv.org/abs/1509.06344
    auto xy = vec2f{position.x, position.y};
    auto uv = vec2f{
        xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
    position = vec3f{uv.x, uv.y, 0} * scale;
  }
}

void make_bulged_disk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps, float scale,
    float uvscale, float height) {
  make_disk(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (height != 0) {
    height      = min(height, scale);
    auto radius = (1 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i = 0; i < positions.size(); i++) {
      auto pn      = normalize(positions[i] - center);
      positions[i] = center + pn * radius;
      normals[i]   = pn;
    }
  }
}

// Generate a uvdisk
void make_uvdisk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float scale, const vec2f& uvscale) {
  make_rect(quads, positions, normals, texcoords, steps, {1, 1}, {1, 1});
  for (auto i = 0; i < positions.size(); i++) {
    auto uv      = texcoords[i];
    auto phi     = 2 * pif * uv.x;
    positions[i] = vec3f{cos(phi) * uv.y, sin(phi) * uv.y, 0} * scale;
    normals[i]   = {0, 0, 1};
    texcoords[i] = uv * uvscale;
  }
}

// Generate a uvcylinder
void make_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec2f& scale, const vec3f& uvscale) {
  auto qquads     = vector<vec4i>{};
  auto qpositions = vector<vec3f>{};
  auto qnormals   = vector<vec3f>{};
  auto qtexcoords = vector<vec2f>{};
  // side
  make_rect(qquads, qpositions, qnormals, qtexcoords, {steps.x, steps.y},
      {1, 1}, {1, 1});
  for (auto i = 0; i < qpositions.size(); i++) {
    auto uv       = qtexcoords[i];
    auto phi      = 2 * pif * uv.x;
    qpositions[i] = {
        cos(phi) * scale.x, sin(phi) * scale.x, (2 * uv.y - 1) * scale.y};
    qnormals[i]   = {cos(phi), sin(phi), 0};
    qtexcoords[i] = uv * vec2f{uvscale.x, uvscale.y};
  }
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexcoords);
  // top
  make_rect(qquads, qpositions, qnormals, qtexcoords, {steps.x, steps.z},
      {1, 1}, {1, 1});
  for (auto i = 0; i < qpositions.size(); i++) {
    auto uv         = qtexcoords[i];
    auto phi        = 2 * pif * uv.x;
    qpositions[i]   = {cos(phi) * uv.y * scale.x, sin(phi) * uv.y * scale.x, 0};
    qnormals[i]     = {0, 0, 1};
    qtexcoords[i]   = uv * vec2f{uvscale.x, uvscale.z};
    qpositions[i].z = scale.y;
  }
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexcoords);
  // bottom
  make_rect(qquads, qpositions, qnormals, qtexcoords, {steps.x, steps.z},
      {1, 1}, {1, 1});
  for (auto i = 0; i < qpositions.size(); i++) {
    auto uv         = qtexcoords[i];
    auto phi        = 2 * pif * uv.x;
    qpositions[i]   = {cos(phi) * uv.y * scale.x, sin(phi) * uv.y * scale.x, 0};
    qnormals[i]     = {0, 0, 1};
    qtexcoords[i]   = uv * vec2f{uvscale.x, uvscale.z};
    qpositions[i].z = -scale.y;
    qnormals[i]     = -qnormals[i];
  }
  for (auto& qquad : qquads) swap(qquad.x, qquad.z);
  merge_quads(quads, positions, normals, texcoords, qquads, qpositions,
      qnormals, qtexcoords);
}

// Generate a uvcylinder
void make_rounded_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec2f& scale, const vec3f& uvscale, float radius) {
  make_uvcylinder(quads, positions, normals, texcoords, steps, scale, uvscale);
  if (radius != 0) {
    radius = min(radius, min(scale));
    auto c = scale - radius;
    for (auto i = 0; i < positions.size(); i++) {
      auto phi = atan2(positions[i].y, positions[i].x);
      auto r   = length(vec2f{positions[i].x, positions[i].y});
      auto z   = positions[i].z;
      auto pc  = vec2f{r, abs(z)};
      auto ps  = (z < 0) ? -1.0f : 1.0f;
      if (pc.x >= c.x && pc.y >= c.y) {
        auto pn      = normalize(pc - c);
        positions[i] = {cos(phi) * (c.x + radius * pn.x),
            sin(phi) * (c.x + radius * pn.x), ps * (c.y + radius * pn.y)};
        normals[i]   = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
      } else {
        continue;
      }
    }
  }
}

// Make a plane.
shape_data make_rect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_rect(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
shape_data make_bulged_rect(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = quads_shape{};
  make_bulged_rect(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a plane in the xz plane.
shape_data make_recty(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_recty(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
shape_data make_bulged_recty(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = quads_shape{};
  make_bulged_recty(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Make a box.
shape_data make_box(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = quads_shape{};
  make_box(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
shape_data make_rounded_box(const vec3i& steps, const vec3f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = quads_shape{};
  make_rounded_box(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a quad stack
shape_data make_rect_stack(
    const vec3i& steps, const vec3f& scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_rect_stack(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a floor.
shape_data make_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
shape_data make_bent_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale, float bent) {
  auto shape = quads_shape{};
  make_bent_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, bent);
  return shape;
}

// Make a sphere.
shape_data make_sphere(int steps, float scale, float uvscale) {
  auto shape = quads_shape{};
  make_sphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere.
shape_data make_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_uvsphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere with slipped caps.
shape_data make_capped_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale, float height) {
  auto shape = quads_shape{};
  make_capped_uvsphere(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, height);
  return shape;
}
// Make a disk
shape_data make_disk(int steps, float scale, float uvscale) {
  auto shape = quads_shape{};
  make_disk(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}

// Make a bulged disk
shape_data make_bulged_disk(
    int steps, float scale, float uvscale, float height) {
  auto shape = quads_shape{};
  make_bulged_disk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, height);
  return shape;
}

// Make a uv disk
shape_data make_uvdisk(const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = quads_shape{};
  make_uvdisk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a uv cylinder
shape_data make_uvcylinder(
    const vec3i& steps, const vec2f& scale, const vec3f& uvscale) {
  auto shape = quads_shape{};
  make_uvcylinder(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(const vec3i& steps, const vec2f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = quads_shape{};
  make_rounded_uvcylinder(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Generate lines set along a quad.
void make_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    const vec2i& steps, const vec2f& size, const vec2f& uvscale,
    const vec2f& rad) {
  auto nverts = (steps.x + 1) * steps.y;
  auto nlines = steps.x * steps.y;
  auto vid    = [steps](int i, int j) { return j * (steps.x + 1) + i; };
  auto fid    = [steps](int i, int j) { return j * steps.x + i; };

  positions.resize(nverts);
  normals.resize(nverts);
  texcoords.resize(nverts);
  radius.resize(nverts);
  if (steps.y > 1) {
    for (auto j = 0; j < steps.y; j++) {
      for (auto i = 0; i <= steps.x; i++) {
        auto uv = vec2f{
            i / (float)steps.x, j / (float)(steps.y > 1 ? steps.y - 1 : 1)};
        positions[vid(i, j)] = {
            (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
        normals[vid(i, j)]   = {1, 0, 0};
        texcoords[vid(i, j)] = uv * uvscale;
      }
    }
  } else {
    for (auto i = 0; i <= steps.x; i++) {
      auto uv              = vec2f{i / (float)steps.x, 0};
      positions[vid(i, 0)] = {(uv.x - 0.5f) * size.x, 0, 0};
      normals[vid(i, 0)]   = {1, 0, 0};
      texcoords[vid(i, 0)] = uv * uvscale;
    }
  }

  lines.resize(nlines);
  for (int j = 0; j < steps.y; j++) {
    for (int i = 0; i < steps.x; i++) {
      lines[fid(i, j)] = {vid(i, j), vid(i + 1, j)};
    }
  }
}

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
shape_data make_lines(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, const vec2f& rad) {
  auto shape = lines_shape{};
  make_lines(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, steps, scale, uvscale, rad);
  return shape;
}

// Generate a point at the origin.
void make_point(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    float point_radius) {
  points    = {0};
  positions = {{0, 0, 0}};
  normals   = {{0, 0, 1}};
  texcoords = {{0, 0}};
  radius    = {point_radius};
}

// Generate a point set with points placed at the origin with texcoords
// varying along u.
void make_points(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    int num, float uvscale, float point_radius) {
  points.resize(num);
  for (auto i = 0; i < num; i++) points[i] = i;
  positions.assign(num, {0, 0, 0});
  normals.assign(num, {0, 0, 1});
  texcoords.assign(num, {0, 0});
  radius.assign(num, point_radius);
  for (auto i = 0; i < texcoords.size(); i++)
    texcoords[i] = {(float)i / (float)num, 0};
}

// Generate a point set.
void make_random_points(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    int num, const vec3f& size, float uvscale, float point_radius,
    uint64_t seed) {
  make_points(points, positions, normals, texcoords, radius, num, uvscale,
      point_radius);
  auto rng = make_rng(seed);
  for (auto& position : positions) {
    position = (rand3f(rng) - vec3f{0.5f, 0.5f, 0.5f}) * size;
  }
}

// Make point primitives. Returns points, pos, norm, texcoord, radius.
shape_data make_point(float radius) {
  auto shape = points_shape{};
  make_point(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, radius);
  return shape;
}

shape_data make_points(int num, float uvscale, float radius) {
  auto shape = points_shape{};
  make_points(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, num, uvscale, radius);
  return shape;
}

shape_data make_random_points(
    int num, const vec3f& size, float uvscale, float radius, uint64_t seed) {
  auto shape = points_shape{};
  make_random_points(shape.points, shape.positions, shape.normals,
      shape.texcoords, shape.radius, num, size, uvscale, radius, seed);
  return shape;
}

// Make a bezier circle. Returns bezier, pos.
void make_bezier_circle(
    vector<vec4i>& beziers, vector<vec3f>& positions, float size) {
  // constant from http://spencermortensen.com/articles/bezier-circle/
  const auto  c              = 0.551915024494f;
  static auto circle_pos     = vector<vec3f>{{1, 0, 0}, {1, c, 0}, {c, 1, 0},
      {0, 1, 0}, {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
      {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
  static auto circle_beziers = vector<vec4i>{
      {0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
  positions = circle_pos;
  beziers   = circle_beziers;
  for (auto& p : positions) p *= size;
}

// Make fvquad
void make_fvrect(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const vec2f& size, const vec2f& uvscale) {
  make_rect(quadspos, positions, normals, texcoords, steps, size, uvscale);
  quadsnorm     = quadspos;
  quadstexcoord = quadspos;
}
void make_fvbox(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec3i& steps,
    const vec3f& size, const vec3f& uvscale) {
  make_box(quadspos, positions, normals, texcoords, steps, size, uvscale);
  quadsnorm                     = quadspos;
  quadstexcoord                 = quadspos;
  std::tie(quadspos, positions) = weld_quads(
      quadspos, positions, 0.1f * min(size) / max(steps));
}
void make_fvsphere(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps, float size,
    float uvscale) {
  make_fvbox(quadspos, quadsnorm, quadstexcoord, positions, normals, texcoords,
      {steps, steps, steps}, {size, size, size}, {uvscale, uvscale, uvscale});
  quadsnorm = quadspos;
  normals   = positions;
  for (auto& n : normals) n = normalize(n);
}

// Make a facevarying rect
fvshape_data make_fvrect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = quads_fvshape{};
  make_fvrect(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying box
fvshape_data make_fvbox(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = quads_fvshape{};
  make_fvbox(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying sphere
fvshape_data make_fvsphere(int steps, float scale, float uvscale) {
  auto shape = quads_fvshape{};
  make_fvsphere(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

extern vector<vec3f> suzanne_positions;
extern vector<vec4i> suzanne_quads;

// Predefined meshes
void make_monkey(vector<vec4i>& quads, vector<vec3f>& positions, float scale) {
  quads     = suzanne_quads;
  positions = suzanne_positions;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

extern vector<vec3f> bunny_positions;
extern vector<vec3f> bunny_normals;
extern vector<vec2f> bunny_texcoords;
extern vector<vec3i> bunny_triangles;

triangles_shape make_bunny(float scale, bool align_middle) {
  auto shape      = triangles_shape{};
  shape.triangles = bunny_triangles;
  shape.positions = bunny_positions;
  shape.normals   = bunny_normals;
  shape.texcoords = bunny_texcoords;
  // scale to height 1
  auto bbox = invalidb3f;
  for (auto& t : shape.triangles) {
    bbox = merge(bbox, triangle_bounds(shape.positions[t.x],
                           shape.positions[t.y], shape.positions[t.z]));
  }
  auto yscale = 2 / size(bbox).y;
  for (auto& p : shape.positions) p *= yscale;
  if (align_middle) {
    for (auto& p : shape.positions) p.y -= 1;
  }
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}

void make_quad(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale) {
  static const auto quad_positions = vector<vec3f>{
      {-1, -1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, +1, 0}};
  static const auto quad_normals = vector<vec3f>{
      {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
  static const auto quad_texcoords = vector<vec2f>{
      {0, 1}, {1, 1}, {1, 0}, {0, 0}};
  static const auto quad_quads = vector<vec4i>{{0, 1, 2, 3}};
  quads                        = quad_quads;
  positions                    = quad_positions;
  normals                      = quad_normals;
  texcoords                    = quad_texcoords;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

void make_quady(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale) {
  static const auto quady_positions = vector<vec3f>{
      {-1, 0, -1}, {-1, 0, +1}, {+1, 0, +1}, {+1, 0, -1}};
  static const auto quady_normals = vector<vec3f>{
      {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
  static const auto quady_texcoords = vector<vec2f>{
      {0, 0}, {1, 0}, {1, 1}, {0, 1}};
  static const auto quady_quads = vector<vec4i>{{0, 1, 2, 3}};
  quads                         = quady_quads;
  positions                     = quady_positions;
  normals                       = quady_normals;
  texcoords                     = quady_texcoords;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

void make_cube(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale) {
  static const auto cube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
      {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
      {+1, +1, -1}, {+1, -1, +1}, {+1, -1, -1}, {+1, +1, -1}, {+1, +1, +1},
      {-1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {-1, +1, -1}, {-1, +1, +1},
      {+1, +1, +1}, {+1, +1, -1}, {-1, +1, -1}, {+1, -1, +1}, {-1, -1, +1},
      {-1, -1, -1}, {+1, -1, -1}};
  static const auto cube_normals   = vector<vec3f>{{0, 0, +1}, {0, 0, +1},
      {0, 0, +1}, {0, 0, +1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1},
      {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {-1, 0, 0}, {-1, 0, 0},
      {-1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0},
      {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
  static const auto cube_texcoords = vector<vec2f>{{0, 1}, {1, 1}, {1, 0},
      {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0},
      {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1},
      {1, 1}, {1, 0}, {0, 0}};
  static const auto cube_quads     = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
      {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
  quads                            = cube_quads;
  positions                        = cube_positions;
  normals                          = cube_normals;
  texcoords                        = cube_texcoords;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

void make_fvcube(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale) {
  static const auto fvcube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
      {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
      {+1, +1, -1}};
  static const auto fvcube_normals   = vector<vec3f>{{0, 0, +1}, {0, 0, +1},
      {0, 0, +1}, {0, 0, +1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1},
      {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {-1, 0, 0}, {-1, 0, 0},
      {-1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0},
      {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
  static const auto fvcube_texcoords = vector<vec2f>{{0, 1}, {1, 1}, {1, 0},
      {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0},
      {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1},
      {1, 1}, {1, 0}, {0, 0}};
  static const auto fvcube_quadspos  = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
      {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};
  static const auto fvcube_quadsnorm = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
      {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
  static const auto fvcube_quadstexcoord = vector<vec4i>{{0, 1, 2, 3},
      {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19},
      {20, 21, 22, 23}};
  quadspos                               = fvcube_quadspos;
  quadsnorm                              = fvcube_quadsnorm;
  quadstexcoord                          = fvcube_quadstexcoord;
  positions                              = fvcube_positions;
  normals                                = fvcube_normals;
  texcoords                              = fvcube_texcoords;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

void make_geosphere(
    vector<vec3i>& triangles, vector<vec3f>& positions, float scale) {
  // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
  const float X                   = 0.525731112119133606f;
  const float Z                   = 0.850650808352039932f;
  static auto geosphere_positions = vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z},
      {-X, 0.0, -Z}, {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X},
      {0.0, -Z, -X}, {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
  static auto geosphere_triangles = vector<vec3i>{{0, 1, 4}, {0, 4, 9},
      {9, 4, 5}, {4, 8, 5}, {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3},
      {5, 3, 2}, {2, 3, 7}, {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0},
      {0, 6, 1}, {6, 10, 1}, {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
  triangles                       = geosphere_triangles;
  positions                       = geosphere_positions;
  if (scale != 1) {
    for (auto& p : positions) p *= scale;
  }
}

// Predefined meshes
shape_data make_monkey(float scale) {
  auto shape = quads_shape{};
  make_monkey(shape.quads, shape.positions, scale);
  return shape;
}
shape_data make_quad(float scale) {
  auto shape = quads_shape{};
  make_quad(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_quady(float scale) {
  auto shape = quads_shape{};
  make_quady(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_cube(float scale) {
  auto shape = quads_shape{};
  make_cube(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
fvshape_data make_fvcube(float scale) {
  auto shape = quads_fvshape{};
  make_fvcube(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
shape_data make_geosphere(float scale) {
  auto shape = triangles_shape{};
  make_geosphere(shape.triangles, shape.positions, scale);
  return shape;
}

// Make a hair ball around a shape
void make_hair(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2i& steps, const vec2f& len,
    const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
  auto alltriangles    = striangles;
  auto quads_triangles = quads_to_triangles(squads);
  alltriangles.insert(
      alltriangles.end(), quads_triangles.begin(), quads_triangles.end());
  auto bpos      = vector<vec3f>{};
  auto bnorm     = vector<vec3f>{};
  auto btexcoord = vector<vec2f>{};
  sample_triangles(bpos, bnorm, btexcoord, alltriangles, spos, snorm, stexcoord,
      steps.y, seed);

  auto rng  = make_rng(seed, 3);
  auto blen = vector<float>(bpos.size());
  for (auto& l : blen) {
    l = lerp(len.x, len.y, rand1f(rng));
  }

  auto cidx = vector<int>();
  if (clump.x > 0) {
    for (auto bidx = 0; bidx < bpos.size(); bidx++) {
      cidx.push_back(0);
      auto cdist = flt_max;
      for (auto c = 0; c < clump.y; c++) {
        auto d = length(bpos[bidx] - bpos[c]);
        if (d < cdist) {
          cdist       = d;
          cidx.back() = c;
        }
      }
    }
  }

  make_lines(lines, positions, normals, texcoords, radius, steps, {1, 1},
      {1, 1}, {1, 1});
  for (auto i = 0; i < positions.size(); i++) {
    auto u       = texcoords[i].x;
    auto bidx    = i / (steps.x + 1);
    positions[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
    normals[i]   = bnorm[bidx];
    radius[i]    = lerp(rad.x, rad.y, u);
    if (clump.x > 0) {
      positions[i] =
          positions[i] +
          (positions[i + (cidx[bidx] - bidx) * (steps.x + 1)] - positions[i]) *
              u * clump.x;
    }
    if (noise.x > 0) {
      auto nx =
          (perlin_noise(positions[i] * noise.y + vec3f{0, 0, 0}) * 2 - 1) *
          noise.x;
      auto ny =
          (perlin_noise(positions[i] * noise.y + vec3f{3, 7, 11}) * 2 - 1) *
          noise.x;
      auto nz =
          (perlin_noise(positions[i] * noise.y + vec3f{13, 17, 19}) * 2 - 1) *
          noise.x;
      positions[i] += {nx, ny, nz};
    }
  }

  if (clump.x > 0 || noise.x > 0 || rotation.x > 0) {
    normals = lines_tangents(lines, positions);
  }
}

// Make a hair ball around a shape
lines_shape make_hair(const triangles_shape& base, const vec2i& steps,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
  auto shape = lines_shape{};
  make_hair(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, base.triangles, base.quads, base.positions, base.normals,
      base.texcoords, steps, len, rad, noise, clump, rotation, seed);
  return shape;
}

// Thickens a shape by copy9ing the shape content, rescaling it and flipping
// its normals. Note that this is very much not robust and only useful for
// trivial cases.
void make_shell(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float thickness) {
  auto bbox = invalidb3f;
  for (auto p : positions) bbox = merge(bbox, p);
  auto center              = yocto::center(bbox);
  auto inner_quads         = quads;
  auto inner_positions     = positions;
  auto inner_normals       = normals;
  auto inner_texturecoords = texcoords;
  for (auto& p : inner_positions) p = (1 - thickness) * (p - center) + center;
  for (auto& n : inner_normals) n = -n;
  merge_quads(quads, positions, normals, texcoords, inner_quads,
      inner_positions, inner_normals, inner_texturecoords);
}

// Make a heightfield mesh.
shape_data make_heightfield(const vec2i& size, const vector<float>& height) {
  auto shape = make_recty(
      size - 1, vec2f{(float)size.x, (float)size.y} / max(size), {1, 1});
  for (auto j = 0; j < size.y; j++)
    for (auto i = 0; i < size.x; i++)
      shape.positions[j * size.x + i].y = height[j * size.x + i];
  shape.normals = quads_normals(shape.quads, shape.positions);
  return shape;
}
shape_data make_heightfield(const vec2i& size, const vector<vec4f>& color) {
  auto shape = make_recty(
      size - 1, vec2f{(float)size.x, (float)size.y} / max(size), {1, 1});
  for (auto j = 0; j < size.y; j++)
    for (auto i = 0; i < size.x; i++)
      shape.positions[j * size.x + i].y = mean(xyz(color[j * size.x + i]));
  shape.normals = quads_normals(shape.quads, shape.positions);
  return shape;
}

}  // namespace yocto
