//
// Implementation for Yocto/Shape
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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
// IMPLEMENTATION FO SHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const shape_data& shape, int element, const vec2f& uv) {
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.positions[point];
  } else if (!shape.lines.empty()) {
    auto [u, _] = uv;
    auto& line  = shape.lines[element];
    return interpolate_line(shape.positions, line, u);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.positions, triangle, uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.positions, quad, uv);
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
    auto [u, _] = uv;
    auto& line  = shape.lines[element];
    return normalize(interpolate_line(shape.normals, line, u));
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return normalize(interpolate_triangle(shape.normals, triangle, uv));
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return normalize(interpolate_quad(shape.normals, quad, uv));
  } else {
    return {0, 0, 1};
  }
}

vec3f eval_tangent(const shape_data& shape, int element, const vec2f& uv) {
  return eval_normal(shape, element, uv);
}

vec2f eval_texcoord(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return uv;
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.texcoords[point];
  } else if (!shape.lines.empty()) {
    auto [u, _] = uv;
    auto& line  = shape.lines[element];
    return interpolate_line(shape.texcoords, line, u);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.texcoords, triangle, uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.texcoords, quad, uv);
  } else {
    return uv;
  }
}

vec4f eval_color(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.colors[point];
  } else if (!shape.lines.empty()) {
    auto [u, _] = uv;
    auto& line  = shape.lines[element];
    return interpolate_line(shape.colors, line, u);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.colors, triangle, uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.colors, quad, uv);
  } else {
    return {0, 0, 0, 0};
  }
}

float eval_radius(const shape_data& shape, int element, const vec2f& uv) {
  if (shape.radius.empty()) return 0;
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.radius[point];
  } else if (!shape.lines.empty()) {
    auto [u, _] = uv;
    auto& line  = shape.lines[element];
    return interpolate_line(shape.radius, line, u);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.radius, triangle, uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.radius, quad, uv);
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
    return line_tangent(shape.positions, line);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return triangle_normal(shape.positions, triangle);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return quad_normal(shape.positions, quad);
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
    auto [ru, _]      = ruv;
    auto [element, u] = sample_lines(cdf, rn, ru);
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
  if (!shape.quads.empty()) {
    result.triangles = quads_to_triangles(shape.quads);
    result.quads     = {};
  }
  return result;
}
void quads_to_triangles_inplace(shape_data& shape) {
  if (shape.quads.empty()) return;
  shape.triangles = quads_to_triangles(shape.quads);
  shape.quads     = {};
}

// Subdivision
shape_data subdivide_shape(
    const shape_data& shape, int subdivisions, bool catmullclark) {
  // This should probably be re-implemented in a faster fashion,
  // but how it is not obvious
  if (subdivisions == 0) return shape;
  auto subdivided = shape_data{};
  if (!subdivided.points.empty()) {
    subdivided = shape;
  } else if (!subdivided.lines.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_lines(
        shape.lines, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_lines(
        shape.lines, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_lines(
        shape.lines, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_lines(
        subdivided.lines, shape.radius, subdivisions);
    std::tie(subdivided.lines, subdivided.positions) = subdivide_lines(
        shape.lines, shape.positions, subdivisions);
  } else if (!subdivided.triangles.empty()) {
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
  } else if (!subdivided.quads.empty() && !catmullclark) {
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
  } else if (!subdivided.quads.empty() && catmullclark) {
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

// Transform shape
shape_data transform_shape(
    const frame3f& frame, const shape_data& shape, bool non_rigid) {
  return transform_shape(frame, shape, 1, non_rigid);
}
shape_data transform_shape(const frame3f& frame, const shape_data& shape,
    float radius_scale, bool non_rigid) {
  auto transformed = shape;
  for (auto& position : transformed.positions)
    position = transform_point(frame, position);
  for (auto& normal : transformed.normals)
    normal = transform_normal(frame, normal, non_rigid);
  for (auto& radius : transformed.radius) radius *= radius_scale;
  return transformed;
}
shape_data remove_normals(const shape_data& shape) {
  auto transformed    = shape;
  transformed.normals = {};
  return transformed;
}

vector<string> shape_stats(const shape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto nums) {
    auto str = string{};
    for (auto num : nums) {
      if (!str.empty()) str += " ";
      str += std::to_string(num);
    }
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
  stats.push_back("center:       " + format3(bbox_center(bbox)));
  stats.push_back("diagonal:     " + format3(bbox_diagonal(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR FVSHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, const vec2f& uv) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return interpolate_quad(shape.positions, quad, uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadsnorm[element];
    return normalize(interpolate_quad(shape.normals, quad, uv));
  } else {
    return {0, 0, 1};
  }
}

vec2f eval_texcoord(const fvshape_data& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return uv;
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadstexcoord[element];
    return interpolate_quad(shape.texcoords, quad, uv);
  } else {
    return uv;
  }
}

// Evaluate element normals
vec3f eval_element_normal(const fvshape_data& shape, int element) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return quad_normal(shape.positions, quad);
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
  // This should be probably re-implemeneted in a faster fashion.
  if (subdivisions == 0) return shape;
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

// Transform shape
fvshape_data transform_fvshape(
    const frame3f& frame, const fvshape_data& shape, bool non_rigid) {
  auto transformed = shape;
  for (auto& position : transformed.positions)
    position = transform_point(frame, position);
  for (auto& normal : transformed.normals)
    normal = transform_normal(frame, normal, non_rigid);
  return transformed;
}

vector<string> fvshape_stats(const fvshape_data& shape, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto nums) {
    auto str = string{};
    for (auto num : nums) {
      if (!str.empty()) str += " ";
      str += std::to_string(num);
    }
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
  stats.push_back("center:       " + format3(bbox_center(bbox)));
  stats.push_back("diagonal:     " + format3(bbox_diagonal(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Predefined meshes
shape_data make_quad(int subdivisions) {
  static const auto quad_positions = vector<vec3f>{
      {-1, -1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, +1, 0}};
  static const auto quad_normals = vector<vec3f>{
      {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
  static const auto quad_texcoords = vector<vec2f>{
      {0, 1}, {1, 1}, {1, 0}, {0, 0}};
  static const auto quad_quads = vector<vec4i>{{0, 1, 2, 3}};
  auto              shape      = shape_data{};
  if (subdivisions == 0) {
    shape.quads     = quad_quads;
    shape.positions = quad_positions;
    shape.normals   = quad_normals;
    shape.texcoords = quad_texcoords;
  } else {
    std::tie(shape.quads, shape.positions) = subdivide_quads(
        quad_quads, quad_positions, subdivisions);
    std::tie(shape.quads, shape.normals) = subdivide_quads(
        quad_quads, quad_normals, subdivisions);
    std::tie(shape.quads, shape.texcoords) = subdivide_quads(
        quad_quads, quad_texcoords, subdivisions);
  }
  return shape;
}
shape_data make_quady(int subdivisions) {
  static const auto quady_positions = vector<vec3f>{
      {-1, 0, -1}, {-1, 0, +1}, {+1, 0, +1}, {+1, 0, -1}};
  static const auto quady_normals = vector<vec3f>{
      {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
  static const auto quady_texcoords = vector<vec2f>{
      {0, 0}, {1, 0}, {1, 1}, {0, 1}};
  static const auto quady_quads = vector<vec4i>{{0, 1, 2, 3}};
  auto              shape       = shape_data{};
  if (subdivisions == 0) {
    shape.quads     = quady_quads;
    shape.positions = quady_positions;
    shape.normals   = quady_normals;
    shape.texcoords = quady_texcoords;
  } else {
    std::tie(shape.quads, shape.positions) = subdivide_quads(
        quady_quads, quady_positions, subdivisions);
    std::tie(shape.quads, shape.normals) = subdivide_quads(
        quady_quads, quady_normals, subdivisions);
    std::tie(shape.quads, shape.texcoords) = subdivide_quads(
        quady_quads, quady_texcoords, subdivisions);
  }
  return shape;
}
shape_data make_cube(int subdivisions) {
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

  auto shape = shape_data{};
  if (subdivisions == 0) {
    shape.quads     = cube_quads;
    shape.positions = cube_positions;
    shape.normals   = cube_normals;
    shape.texcoords = cube_texcoords;
  } else {
    std::tie(shape.quads, shape.positions) = subdivide_quads(
        cube_quads, cube_positions, subdivisions);
    std::tie(shape.quads, shape.normals) = subdivide_quads(
        cube_quads, cube_normals, subdivisions);
    std::tie(shape.quads, shape.texcoords) = subdivide_quads(
        cube_quads, cube_texcoords, subdivisions);
  }
  return shape;
}
shape_data make_sphere(int subdivisions) {
  auto shape = make_cube(subdivisions);
  for (auto& p : shape.positions) p = normalize(p);
  shape.normals = shape.positions;
  return shape;
}
shape_data make_geosphere(int subdivisions) {
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

  auto shape = shape_data{};
  if (subdivisions == 0) {
    shape.triangles = geosphere_triangles;
    shape.positions = geosphere_positions;
    shape.normals   = geosphere_positions;
  } else {
    std::tie(shape.triangles, shape.positions) = subdivide_triangles(
        geosphere_triangles, geosphere_positions, subdivisions);
    for (auto& position : shape.positions) position = normalize(position);
    shape.normals = shape.positions;
  }
  return shape;
}
shape_data make_disk(int subdivisions) {
  auto shape = make_quad(subdivisions);
  for (auto& position : shape.positions) {
    // Analytical Methods for Squaring the Disc, by C. Fong
    // https://arxiv.org/abs/1509.06344
    auto [x, y, _] = position;
    auto uv        = vec2f{x * sqrt(1 - y * y / 2), y * sqrt(1 - x * x / 2)};
    position       = vec3f{uv, 0};
  }
  return shape;
}
shape_data make_wtcube(int subdivisions) {
  static const auto wtcube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
      {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
      {+1, +1, -1}};
  static const auto wtcube_quads     = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
          {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};

  auto shape = shape_data{};
  if (subdivisions == 0) {
    shape.quads     = wtcube_quads;
    shape.positions = wtcube_positions;
  } else {
    std::tie(shape.quads, shape.positions) = subdivide_quads(
        wtcube_quads, wtcube_positions, subdivisions);
  }
  return shape;
}
shape_data make_wtsphere(int subdivisions) {
  auto shape = make_wtcube(subdivisions);
  for (auto& p : shape.positions) p = normalize(p);
  shape.normals = shape.positions;
  return shape;
}
shape_data make_floor(int subdivisions, float scale) {
  auto shape = make_quady(subdivisions);
  for (auto& position : shape.positions) position *= scale;
  for (auto& texcoord : shape.texcoords) texcoord *= scale;
  return shape;
}
shape_data make_monkey(int subdivisions) {
  extern vector<vec3f> suzanne_positions;
  extern vector<vec4i> suzanne_quads;

  auto shape = shape_data{};
  if (subdivisions == 0) {
    shape.quads     = suzanne_quads;
    shape.positions = suzanne_positions;
  } else {
    std::tie(shape.quads, shape.positions) = subdivide_quads(
        suzanne_quads, suzanne_positions, subdivisions);
  }
  return shape;
}

// Deformed shapes
shape_data make_bulged_quad(int subdivisions, float height) {
  if (height == 0) return make_quad(subdivisions);
  auto shape  = make_quad(subdivisions);
  height      = min(height, 1);
  auto radius = (1 + height * height) / (2 * height);
  auto center = vec3f{0, 0, -radius + height};
  for (auto i : range(shape.positions.size())) {
    auto pn            = normalize(shape.positions[i] - center);
    shape.positions[i] = center + pn * radius;
    shape.normals[i]   = pn;
  }
  return shape;
}
shape_data make_bulged_disk(int subdivisions, float height) {
  if (height == 0) return make_disk(subdivisions);
  auto shape  = make_disk(subdivisions);
  height      = min(height, 1);
  auto radius = (1 + height * height) / (2 * height);
  auto center = vec3f{0, 0, -radius + height};
  for (auto i : range(shape.positions.size())) {
    auto pn            = normalize(shape.positions[i] - center);
    shape.positions[i] = center + pn * radius;
    shape.normals[i]   = pn;
  }
  return shape;
}
shape_data make_rounded_cube(int subdivisions, float radius) {
  if (radius == 0) return make_cube(subdivisions);
  auto shape = make_cube(subdivisions);
  radius     = min(radius, 1);
  auto c     = 1 - radius;
  for (auto i : range(shape.positions.size())) {
    auto pc = abs(shape.positions[i]);
    auto ps = select(component_less(shape.positions[i], 0), -1.0f, 1.0f);
    auto [pc_x, pc_y, pc_z] = pc;
    if (pc_x >= c && pc_y >= c && pc_z >= c) {
      auto pn            = normalize(pc - c);
      shape.positions[i] = c + radius * pn;
      shape.normals[i]   = pn;
    } else if (pc_x >= c && pc_y >= c) {
      auto pn              = normalize((pc - c) * vec3f{1, 1, 0});
      auto [pn_x, pn_y, _] = pn;
      shape.positions[i]   = {c + radius * pn_x, c + radius * pn_y, pc_z};
      shape.normals[i]     = pn;
    } else if (pc_x >= c && pc_z >= c) {
      auto pn              = normalize((pc - c) * vec3f{1, 0, 1});
      auto [pn_x, _, pn_z] = pn;
      shape.positions[i]   = {c + radius * pn_x, pc_y, c + radius * pn_z};
      shape.normals[i]     = pn;
    } else if (pc_y >= c && pc_z >= c) {
      auto pn              = normalize((pc - c) * vec3f{0, 1, 1});
      auto [_, pn_y, pn_z] = pn;
      shape.positions[i]   = {pc_x, c + radius * pn_y, c + radius * pn_z};
      shape.normals[i]     = pn;
    } else {
      continue;
    }
    shape.positions[i] *= ps;
    shape.normals[i] *= ps;
  }
  return shape;
}
shape_data make_bent_floor(int subdivisions, float scale, float radius) {
  if (radius == 0) return make_floor(subdivisions, scale);
  auto shape = make_floor(subdivisions, scale);
  radius     = min(radius * scale, scale);
  auto start = (scale - radius) / 2;
  auto end   = start + radius;
  for (auto i : range(shape.positions.size())) {
    auto [px, py, pz] = shape.positions[i];
    if (pz < -end) {
      shape.positions[i] = {px, -pz - end + radius, -end};
      shape.normals[i]   = {0, 0, 1};
    } else if (pz < -start && pz >= -end) {
      auto phi           = (pif / 2) * (-pz - start) / radius;
      shape.positions[i] = {
          px, -cos(phi) * radius + radius, -sin(phi) * radius - start};
      shape.normals[i] = {0, cos(phi), sin(phi)};
    } else {
    }
  }
  return shape;
}

// Predefined meshes
fvshape_data make_fvcube(int subdivisions) {
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

  auto shape = fvshape_data{};
  if (subdivisions == 0) {
    shape.quadspos      = fvcube_quadspos;
    shape.quadsnorm     = fvcube_quadsnorm;
    shape.quadstexcoord = fvcube_quadstexcoord;
    shape.positions     = fvcube_positions;
    shape.normals       = fvcube_normals;
    shape.texcoords     = fvcube_texcoords;
  } else {
    std::tie(shape.quadspos, shape.positions) = subdivide_quads(
        fvcube_quadspos, fvcube_positions, subdivisions);
    std::tie(shape.quadsnorm, shape.normals) = subdivide_quads(
        fvcube_quadsnorm, fvcube_normals, subdivisions);
    std::tie(shape.quadstexcoord, shape.texcoords) = subdivide_quads(
        fvcube_quadstexcoord, fvcube_texcoords, subdivisions);
  }
  return shape;
}

// Make a tesselated rectangle. Useful in other subdivisions.
struct make_quads_vertex {
  vec3f position;
  vec3f normal;
  vec2f texcoord;
};
template <typename Vertex>
static shape_data make_quads(const vec2i& steps, Vertex&& vertex) {
  auto shape = shape_data{};

  auto [stride, _] = steps;
  auto vid = [stride = stride](int i, int j) { return j * (stride + 1) + i; };
  auto fid = [stride = stride](int i, int j) { return j * stride + i; };

  shape.positions.resize(prod(steps + 1));
  shape.normals.resize(prod(steps + 1));
  shape.texcoords.resize(prod(steps + 1));
  for (auto [i, j] : range(steps + 1)) {
    auto [position, normal, texcoord] = vertex(vec2i(i, j) / (vec2f)steps);
    shape.positions[vid(i, j)]        = position;
    shape.normals[vid(i, j)]          = normal;
    shape.texcoords[vid(i, j)]        = texcoord;
  }

  shape.quads.resize(prod(steps));
  for (auto [i, j] : range(steps)) {
    shape.quads[fid(i, j)] = {
        vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
  }

  return shape;
}

// Merge shape elements
void merge_shape_inplace(shape_data& shape, const shape_data& merge) {
  auto offset = (int)shape.positions.size();
  for (auto& p : merge.points) shape.points.push_back(p + offset);
  for (auto& l : merge.lines) shape.lines.push_back(l + offset);
  for (auto& t : merge.triangles) shape.triangles.push_back(t + offset);
  for (auto& q : merge.quads) shape.quads.push_back(q + offset);
  shape.positions.insert(
      shape.positions.end(), merge.positions.begin(), merge.positions.end());
  shape.normals.insert(
      shape.normals.end(), merge.normals.begin(), merge.normals.end());
  shape.tangents.insert(
      shape.tangents.end(), merge.tangents.begin(), merge.tangents.end());
  shape.texcoords.insert(
      shape.texcoords.end(), merge.texcoords.begin(), merge.texcoords.end());
  shape.colors.insert(
      shape.colors.end(), merge.colors.begin(), merge.colors.end());
  shape.radius.insert(
      shape.radius.end(), merge.radius.begin(), merge.radius.end());
}

// Flip the y and z axis for example shapes
shape_data _flip_yz(const shape_data& shape) {
  auto transformed = shape;
  for (auto& position : transformed.positions) {
    auto [px, py, pz] = position;
    position          = {px, pz, -py};
  }
  for (auto& normal : transformed.normals) {
    auto [nx, ny, nz] = normal;
    normal            = {nx, nz, -ny};
  }
  return transformed;
}

// Make a plane.
shape_data make_rect(const vec2i& steps, const vec2f& scale) {
  return make_quads(steps, [=](vec2f uv) {
    return make_quads_vertex{
        vec3f(scale * (uv * 2 - 1), 0), vec3f(0, 0, 1), flip_v(uv)};
  });
}
shape_data make_bulged_rect(
    const vec2i& steps, const vec2f& scale, float height) {
  if (height == 0) return make_rect(steps, scale);
  auto shape  = make_rect(steps, scale);
  height      = min(height, min(scale));
  auto radius = (1 + height * height) / (2 * height);
  auto center = vec3f{0, 0, -radius + height};
  for (auto i : range(shape.positions.size())) {
    auto pn            = normalize(shape.positions[i] - center);
    shape.positions[i] = center + pn * radius;
    shape.normals[i]   = pn;
  }
  return shape;
}

// Make a plane in the xz plane.
shape_data make_recty(const vec2i& steps, const vec2f& scale) {
  return _flip_yz(make_rect(steps, scale));
}
shape_data make_bulged_recty(
    const vec2i& steps, const vec2f& scale, float height) {
  return _flip_yz(make_bulged_rect(steps, scale, height));
}

// Make a box.
shape_data make_box(const vec3i& steps, const vec3f& scale) {
  auto shape = shape_data{};
  // steps
  auto [steps_x, steps_y, steps_z] = steps;
  // + z
  merge_shape_inplace(shape, make_quads({steps_x, steps_y}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(+sx * (u * 2 - 1), +sy * (v * 2 - 1), +sz),
        vec3f(0, 0, 1), flip_v(uv)};
  }));
  // - z
  merge_shape_inplace(shape, make_quads({steps_x, steps_y}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(-sx * (u * 2 - 1), +sy * (v * 2 - 1), -sz),
        vec3f(0, 0, -1), flip_v(uv)};
  }));
  // + x
  merge_shape_inplace(shape, make_quads({steps_z, steps_y}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(+sx, +sy * (v * 2 - 1), -sz * (u * 2 - 1)),
        vec3f(1, 0, 0), flip_v(uv)};
  }));
  // - x
  merge_shape_inplace(shape, make_quads({steps_z, steps_y}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(-sx, +sy * (v * 2 - 1), +sz * (u * 2 - 1)),
        vec3f(-1, 0, 0), flip_v(uv)};
  }));
  // + y
  merge_shape_inplace(shape, make_quads({steps_x, steps_z}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(+sx * (u * 2 - 1), +sy, -sz * (v * 2 - 1)),
        vec3f(0, +1, 0), flip_v(uv)};
  }));
  // - y
  merge_shape_inplace(shape, make_quads({steps_x, steps_z}, [=](vec2f uv) {
    auto [sx, sy, sz] = scale;
    auto [u, v]       = uv;
    return make_quads_vertex{vec3f(+sx * (u * 2 - 1), -sy, +sz * (v * 2 - 1)),
        vec3f(0, -1, 0), flip_v(uv)};
  }));
  return shape;
}
shape_data make_rounded_box(
    const vec3i& steps, const vec3f& scale, float radius) {
  if (radius == 0) make_box(steps, scale);
  auto shape = make_box(steps, scale);
  radius     = min(radius, min(scale));
  auto c     = scale - radius;
  for (auto i : range(shape.positions.size())) {
    auto pc = abs(shape.positions[i]);
    auto ps = select(component_less(shape.positions[i], 0), -1.0f, 1.0f);
    auto [pc_x, pc_y, pc_z] = pc;
    auto [c_x, c_y, c_z]    = pc;
    if (pc_x >= c_x && pc_y >= c_y && pc_z >= c_z) {
      auto pn            = normalize(pc - c);
      shape.positions[i] = c + radius * pn;
      shape.normals[i]   = pn;
    } else if (pc_x >= c_x && pc_y >= c_y) {
      auto pn              = normalize((pc - c) * vec3f{1, 1, 0});
      auto [pn_x, pn_y, _] = pn;
      shape.positions[i]   = {c_x + radius * pn_x, c_y + radius * pn_y, pc_z};
      shape.normals[i]     = pn;
    } else if (pc_x >= c_x && pc_z >= c_z) {
      auto pn              = normalize((pc - c) * vec3f{1, 0, 1});
      auto [pn_x, _, pn_z] = pn;
      shape.positions[i]   = {c_x + radius * pn_x, pc_y, c_z + radius * pn_z};
      shape.normals[i]     = pn;
    } else if (pc_y >= c_y && pc_z >= c_z) {
      auto pn              = normalize((pc - c) * vec3f{0, 1, 1});
      auto [_, pn_y, pn_z] = pn;
      shape.positions[i]   = {pc_x, c_y + radius * pn_y, c_z + radius * pn_z};
      shape.normals[i]     = pn;
    } else {
      continue;
    }
    shape.positions[i] *= ps;
    shape.normals[i] *= ps;
  }
  return shape;
}

// Make a quad stack
shape_data make_rect_stack(int vsteps, const vec2i& steps, const vec3f& scale) {
  auto shape = shape_data{};
  for (auto i : range(vsteps + 1)) {
    auto w = i / (float)vsteps;
    merge_shape_inplace(shape, make_quads(steps, [=](vec2f uv) {
      auto uvw = vec3f{uv, w};
      return make_quads_vertex{
          scale * (uvw * 2 - 1), vec3f(0, 0, 1), flip_v(uv)};
    }));
  }
  return shape;
}

// Make a sphere.
shape_data make_tsphere(const vec3i& steps, float scale) {
  auto shape = make_box(steps, {scale, scale, scale});
  for (auto& p : shape.positions) p = normalize(p) * scale;
  shape.normals = shape.positions;
  for (auto& n : shape.normals) n = normalize(n);
  return shape;
}

// Make a sphere.
shape_data make_uvsphere(const vec2i& steps, float scale) {
  return make_quads(steps, [=](vec2f uv) {
    auto [phi, theta] = flip_v(uv) * vec2f{2 * pif, pif};
    return make_quads_vertex{
        vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)} * scale,
        vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)}, uv};
  });
}

// Make a sphere.
shape_data make_uvspherey(const vec2i& steps, float scale) {
  return _flip_yz(make_uvsphere(steps, scale));
}

// Make a sphere with slipped caps.
shape_data make_capped_uvsphere(const vec2i& steps, float scale, float cap) {
  if (cap != 0) return make_uvsphere(steps, scale);
  auto shape = make_uvsphere(steps, scale);
  cap        = min(cap, scale / 2);
  auto zflip = (scale - cap);
  for (auto i : range(shape.positions.size())) {
    auto [px, py, pz] = shape.positions[i];
    auto [nx, ny, nz] = shape.normals[i];
    if (pz > zflip) {
      shape.positions[i] = {px, py, 2 * zflip - pz};
      shape.normals[i]   = {-nx, -ny, nz};
    } else if (pz < -zflip) {
      shape.positions[i] = {px, py, 2 * (-zflip) - pz};
      shape.normals[i]   = {-nx, -ny, nz};
    }
  }
  return shape;
}

// Make a sphere with slipped caps.
shape_data make_capped_uvspherey(const vec2i& steps, float scale, float cap) {
  return _flip_yz(make_capped_uvsphere(steps, scale, cap));
}

// Make a uv disk
shape_data make_uvdisk(const vec2i& steps, float scale) {
  return make_quads(steps, [=](vec2f uv) {
    auto [phi, r] = flip_v(uv) * vec2f{2 * pif, 1};
    return make_quads_vertex{vec3f{r * cos(phi), r * sin(phi), 0} * scale,
        vec3f(0, 0, 1), flip_v(uv)};
  });
}

// Make a uv cylinder
shape_data make_uvcylinder(const vec3i& steps, const vec2f& scale) {
  auto shape = shape_data{};
  // steps
  auto [steps_circle, steps_sides, steps_ends] = steps;
  // side
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_sides}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, h]         = uv * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{cos(phi) * radius, sin(phi) * radius, height * (2 * h - 1)},
            vec3f{cos(phi), sin(phi), 0}, uv};
      }));
  // top
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_ends}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, r]         = flip_v(uv) * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{r * radius * cos(phi), r * radius * sin(phi), +height},
            vec3f(0, 0, 1), flip_v(uv)};
      }));
  // bottom
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_ends}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, r]         = uv * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{r * radius * cos(phi), r * radius * sin(phi), -height},
            vec3f(0, 0, -1), uv};
      }));
  return shape;
}

// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(
    const vec3i& steps, const vec2f& scale, float radius) {
  auto shape = make_uvcylinder(steps, scale);
  if (radius != 0) {
    radius = min(radius, min(scale));
    auto c = scale - radius;
    for (auto i : range(shape.positions.size())) {
      auto [px, py, pz] = shape.positions[i];
      auto phi          = atan2(py, px);
      auto r            = length(vec2f{px, py});
      auto z            = pz;
      auto pc           = vec2f{r, abs(z)};
      auto [pc_x, pc_y] = pc;
      auto [c_x, c_y]   = c;
      auto ps           = (z < 0) ? -1.0f : 1.0f;
      if (pc_x >= c_x && pc_y >= c_y) {
        auto pn            = normalize(pc - c);
        auto [pn_x, pn_y]  = pn;
        shape.positions[i] = {cos(phi) * (c_x + radius * pn_x),
            sin(phi) * (c_x + radius * pn_x), ps * (c_y + radius * pn_y)};
        shape.normals[i]   = {cos(phi) * pn_x, sin(phi) * pn_x, ps * pn_y};
      } else {
        continue;
      }
    }
  }
  return shape;
}

// Make a uv capsule
shape_data make_uvcapsule(const vec3i& steps, const vec2f& scale) {
  auto shape = shape_data{};
  // steps
  auto [steps_circle, steps_sides, steps_ends] = steps;
  // side
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_sides}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, h]         = uv * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{cos(phi) * radius, sin(phi) * radius, height * (2 * h - 1)},
            vec3f{cos(phi), sin(phi), 0}, uv};
      }));
  // top
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_ends}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, theta]     = flip_v(uv) * vec2f{2 * pif, pif / 2};
        return make_quads_vertex{
            vec3f{cos(phi) * sin(theta) * radius,
                sin(phi) * sin(theta) * radius, cos(theta) * radius + height},
            vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)},
            flip_v(uv)};
      }));
  // bottom
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_ends}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, theta]     = uv * vec2f{2 * pif, pif / 2};
        return make_quads_vertex{
            vec3f{cos(phi) * sin(theta) * radius,
                sin(phi) * sin(theta) * radius, -cos(theta) * radius - height},
            vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), -cos(theta)},
            uv};
      }));
  // done
  return shape;
}

// Make a uv cone
shape_data make_uvcone(const vec3i& steps, const vec2f& scale) {
  auto shape = shape_data{};
  // steps
  auto [steps_circle, steps_sides, steps_ends] = steps;
  // side
  merge_shape_inplace(
      shape, make_quads({steps_circle, steps_sides}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [nr, nh]         = normalize(scale * vec2f{1, 2});
        auto [phi, h]         = uv * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{cos(phi) * (1 - h) * radius, sin(phi) * (1 - h) * radius,
                height * (2 * h - 1)},
            vec3f{cos(phi) * nh, sin(phi) * nh, nr}, uv};
      }));
  // bottom
  merge_shape_inplace(
      shape, make_quads({steps_sides, steps_ends}, [=](vec2f uv) {
        auto [radius, height] = scale;
        auto [phi, r]         = uv * vec2f{2 * pif, 1};
        return make_quads_vertex{
            vec3f{r * radius * cos(phi), r * radius * sin(phi), -height},
            vec3f(0, 0, -1), uv};
      }));
  return shape;
}

// Make a face-varying rect
fvshape_data make_fvrect(const vec2i& steps, const vec2f& scale) {
  auto rect           = make_rect(steps, scale);
  auto shape          = fvshape_data{};
  shape.positions     = rect.positions;
  shape.normals       = rect.normals;
  shape.texcoords     = rect.texcoords;
  shape.quadspos      = rect.quads;
  shape.quadsnorm     = rect.quads;
  shape.quadstexcoord = rect.quads;
  return shape;
}

// Make a face-varying box
fvshape_data make_fvbox(const vec3i& steps, const vec3f& scale) {
  auto box                                  = make_box(steps, scale);
  auto shape                                = fvshape_data{};
  shape.quadsnorm                           = box.quads;
  shape.quadstexcoord                       = box.quads;
  std::tie(shape.quadspos, shape.positions) = weld_quads(
      box.quads, box.positions, 0.1f * min(scale) / max(steps));
  return shape;
}

// Make a face-varying sphere
fvshape_data make_fvsphere(int steps, float scale) {
  auto shape      = make_fvbox({steps, steps, steps}, {scale, scale, scale});
  shape.quadsnorm = shape.quadspos;
  shape.normals   = shape.positions;
  for (auto& n : shape.normals) n = normalize(n);
  return shape;
}

// Generate lines.
struct make_lines_vertex {
  vec3f position;
  vec3f normal;
  vec2f texcoord;
  float radius;
};
// Generate lines set along a quad.
template <typename Func>
static shape_data make_lines(int num, int steps, Func&& func) {
  auto shape = shape_data{};

  auto vid = [stride = steps](int i, int j) { return j * (stride + 1) + i; };
  auto lid = [stride = steps](int i, int j) { return j * stride + i; };

  shape.positions.resize((steps + 1) * num);
  shape.normals.resize((steps + 1) * num);
  shape.texcoords.resize((steps + 1) * num);
  shape.radius.resize((steps + 1) * num);
  for (auto [i, j] : range(vec2i(steps + 1, num))) {
    auto uv = vec2f{i / (float)steps, j / (float)(num - 1)};
    auto [position, normal, texcoord, radius] = func(uv);
    shape.positions[vid(i, j)]                = position;
    shape.normals[vid(i, j)]                  = normal;
    shape.texcoords[vid(i, j)]                = texcoord;
    shape.radius[vid(i, j)]                   = radius;
  }

  shape.lines.resize(steps * num);
  for (auto [i, j] : range(vec2i(steps, num))) {
    shape.lines[lid(i, j)] = {vid(i, j), vid(i + 1, j)};
  }

  return shape;
}
// Generate lines set along a quad.
shape_data make_lines(
    int num, int steps, const vec2f& scale, const vec2f& rad) {
  return make_lines(num, steps, [&](const vec2f& uv) {
    auto [len, width]              = scale;
    auto [base_radius, tip_radius] = rad;
    auto [u, v]                    = uv;
    return make_lines_vertex{vec3f{(2 * u - 1) * len, (2 * v - 1) * width, 0},
        vec3f{1, 0, 0}, vec2f{u, v}, lerp(base_radius, tip_radius, u)};
  });
}

// Make point primitives. Returns points, pos, norm, texcoord, radius.
shape_data make_point(float radius, bool generate_uv) {
  auto shape      = shape_data{};
  shape.points    = {0};
  shape.positions = {{0, 0, 0}};
  shape.normals   = {{0, 0, 1}};
  if (generate_uv) shape.texcoords = {{0, 0}};
  shape.radius = {radius};
  return shape;
}

// Generate a point set with points placed at the origin with texcoords
// varying along u.
shape_data make_points(int num, float radius, bool generate_uv) {
  auto shape = shape_data{};
  shape.points.resize(num);
  for (auto i : range(num)) shape.points[i] = i;
  shape.positions.assign(num, {0, 0, 0});
  shape.normals.assign(num, {0, 0, 1});
  shape.texcoords.assign(num, {0, 0});
  shape.radius.assign(num, radius);
  if (generate_uv) {
    for (auto i : range(shape.texcoords.size()))
      shape.texcoords[i] = {(float)i / (float)num, 0};
  }
  return shape;
}

shape_data make_quad_grid(const vec2i& steps, float scale_) {
  auto scale = scale_ * vec2f{aspect_ratio((vec2f)steps), 1};
  return make_quads(steps, [=](vec2f uv) {
    return make_quads_vertex{
        vec3f(scale * (uv * 2 - 1), 0), vec3f(0, 0, 1), flip_v(uv)};
  });
}

shape_data make_point_grid(const vec2i& steps, float scale, float radius) {
  auto shape  = make_quad_grid(steps, scale);
  shape.quads = {};
  shape.points.resize(shape.positions.size());
  shape.radius.resize(shape.positions.size());
  for (auto i : range(shape.positions.size())) shape.points[i] = (int)i;
  for (auto i : range(shape.texcoords.size())) shape.radius[i] = radius;
  return shape;
}

shape_data make_line_grid(const vec2i& steps, float scale, float radius) {
  auto shape  = make_quad_grid(steps, scale);
  shape.lines = get_edges(shape.quads);
  shape.quads = {};
  shape.radius.resize(shape.positions.size());
  for (auto i : range(shape.texcoords.size())) shape.radius[i] = radius;
  return shape;
}

// Make points inside a cube
shape_data make_random_points(
    int num, float radius, bool generate_uv, uint64_t seed) {
  auto shape = make_points(num, radius);
  auto rng   = make_rng(seed);
  for (auto& position : shape.positions) position = (2 * rand3f(rng) - 1);
  for (auto& texcoord : shape.texcoords) texcoord = rand2f(rng);
  return shape;
}

// Grow lines around a shape
shape_data make_random_points(
    const shape_data& base, int num, float radius, bool generate_uv, int seed) {
  auto points     = sample_shape(base, num, seed);
  auto bpositions = vector<vec3f>{};
  auto btexcoord  = vector<vec2f>{};
  for (auto& point : points) {
    bpositions.push_back(eval_position(base, point.element, point.uv));
    btexcoord.push_back(eval_texcoord(base, point.element, point.uv));
  }
  auto shape      = make_points(num, radius, generate_uv);
  shape.positions = bpositions;
  return shape;
}

// Grow lines around a shape
shape_data make_random_lines(const shape_data& base, int num, int steps,
    const vec2f& len, const vec2f& radius, bool generate_uv, int seed) {
  auto points     = sample_shape(base, num, seed);
  auto bpositions = vector<vec3f>{};
  auto bnormals   = vector<vec3f>{};
  auto btexcoord  = vector<vec2f>{};
  for (auto& point : points) {
    bpositions.push_back(eval_position(base, point.element, point.uv));
    bnormals.push_back(eval_normal(base, point.element, point.uv));
    btexcoord.push_back(eval_texcoord(base, point.element, point.uv));
  }

  auto [min_len, max_len] = len;

  auto shape = make_lines(num, steps, {1, 1}, radius);
  auto rng   = make_rng(seed);
  for (auto idx : range(num)) {
    auto offset = idx * (steps + 1);
    auto length = lerp(min_len, max_len, rand1f(rng));
    for (auto iidx : range(steps + 1)) {
      auto u                         = iidx / (float)steps;
      shape.positions[offset + iidx] = bpositions[idx] +
                                       u * length * bnormals[idx];
    }
  }

  shape.normals = lines_tangents(shape.lines, shape.positions);

  return shape;
}

// Grow hairs around a shape
shape_data make_random_hairs(const shape_data& base, int num, int steps,
    const vec2f& len, const vec2f& radius, float noise, float gravity,
    bool generate_uv, int seed) {
  auto points     = sample_shape(base, num, seed);
  auto bpositions = vector<vec3f>{};
  auto bnormals   = vector<vec3f>{};
  auto btexcoord  = vector<vec2f>{};
  for (auto& point : points) {
    bpositions.push_back(eval_position(base, point.element, point.uv));
    bnormals.push_back(eval_normal(base, point.element, point.uv));
    btexcoord.push_back(eval_texcoord(base, point.element, point.uv));
  }

  auto [min_len, max_len] = len;

  auto shape = make_lines(num, steps, {1, 1}, radius);
  auto rng   = make_rng(seed);
  for (auto idx : range(num)) {
    auto offset             = idx * (steps + 1);
    auto position           = bpositions[idx];
    auto direction          = bnormals[idx];
    auto length             = lerp(min_len, max_len, rand1f(rng));
    shape.positions[offset] = position;
    for (auto iidx = 1; iidx <= steps; iidx++) {
      shape.positions[offset + iidx] = position;
      shape.positions[offset + iidx] += direction * length / (float)steps;
      shape.positions[offset + iidx] += (2 * rand3f(rng) - 1) * noise;
      shape.positions[offset + iidx] += vec3f{0, -gravity, 0};
      direction = normalize(shape.positions[offset + iidx] - position);
      position  = shape.positions[offset + iidx];
    }
  }

  shape.normals = lines_tangents(shape.lines, shape.positions);

  return shape;
}

// Make a heightfield mesh.
shape_data make_heightfield(const array2d<float>& height) {
  auto size       = (vec2i)height.extents();
  auto shape      = make_recty(size + 1, (vec2f)size / (float)max(size));
  auto [width, _] = size;
  for (auto [i, j] : range(size)) {
    auto& position = shape.positions[j * width + i];
    auto [x, y, z] = position;
    position       = {x, height[vec2i{i, j}], z};
  }
  shape.normals = quads_normals(shape.quads, shape.positions);
  return shape;
}
shape_data make_heightfield(const array2d<vec4f>& height) {
  auto size       = (vec2i)height.extents();
  auto shape      = make_recty(size + 1, (vec2f)size / (float)max(size));
  auto [width, _] = size;
  for (auto [i, j] : range(size)) {
    auto& position = shape.positions[j * width + i];
    auto [x, y, z] = position;
    position       = {x, mean(xyz(height[vec2i{i, j}])), z};
  }
  shape.normals = quads_normals(shape.quads, shape.positions);
  return shape;
}

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res.
shape_data points_to_spheres(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  for (auto& vertex : vertices) {
    auto sphere = make_tsphere({steps, steps, steps}, scale);
    for (auto& position : sphere.positions) position += vertex;
    merge_shape_inplace(shape, sphere);
  }
  return shape;
}
shape_data polyline_to_cylinders(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  for (auto idx = 0; idx < (int)vertices.size() - 1; idx++) {
    auto cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1});
    auto frame    = frame_fromz((vertices[idx] + vertices[idx + 1]) / 2,
           vertices[idx] - vertices[idx + 1]);
    auto length   = distance(vertices[idx], vertices[idx + 1]);
    for (auto& position : cylinder.positions)
      position = transform_point(frame, position * vec3f{1, 1, length / 2});
    for (auto& normal : cylinder.normals)
      normal = transform_direction(frame, normal);
    merge_shape_inplace(shape, cylinder);
  }
  return shape;
}
shape_data lines_to_cylinders(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  for (auto idx = 0; idx < (int)vertices.size(); idx += 2) {
    auto cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1});
    auto frame    = frame_fromz((vertices[idx + 0] + vertices[idx + 1]) / 2,
           vertices[idx + 0] - vertices[idx + 1]);
    auto length   = distance(vertices[idx + 0], vertices[idx + 1]);
    for (auto& position : cylinder.positions)
      position = transform_point(frame, position * vec3f{1, 1, length / 2});
    for (auto& normal : cylinder.normals)
      normal = transform_direction(frame, normal);
    merge_shape_inplace(shape, cylinder);
  }
  return shape;
}
shape_data lines_to_cylinders(const vector<vec2i>& lines,
    const vector<vec3f>& positions, int steps, float scale) {
  auto shape = shape_data{};
  for (auto& [v1, v2] : lines) {
    auto &p1 = positions[v1], &p2 = positions[v2];
    auto  cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1});
    auto  frame    = frame_fromz((p1 + p2) / 2, p1 - p2);
    auto  length   = distance(p1, p2);
    for (auto& position : cylinder.positions)
      position = transform_point(frame, position * vec3f{1, 1, length / 2});
    for (auto& normal : cylinder.normals)
      normal = transform_direction(frame, normal);
    merge_shape_inplace(shape, cylinder);
  }
  return shape;
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
  for (auto& tangent : tangents) tangent = {0, 0, 0};
  for (auto& [v1, v2] : lines) {
    auto &p1 = positions[v1], &p2 = positions[v2];
    auto  tangent = line_tangent(p1, p2);
    auto  length  = line_length(p1, p2);
    tangents[v1] += tangent * length;
    tangents[v2] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
  return tangents;
}

// Compute per-vertex normals for triangles.
vector<vec3f> triangles_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& triangle : triangles) {
    auto normal = triangle_normal(positions, triangle);
    auto area   = triangle_area(positions, triangle);
    for (auto vid : triangle) normals[vid] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Compute per-vertex normals for quads.
vector<vec3f> quads_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& quad : quads) {
    auto normal = quad_normal(positions, quad);
    auto area   = quad_area(positions, quad);
    if (!is_triangle(quad)) {
      for (auto vid : quad) normals[vid] += normal * area;
    } else {
      for (auto vid : as_triangle(quad)) normals[vid] += normal * area;
    }
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
  for (auto& tangent : tangents) tangent = {0, 0, 0};
  for (auto& line : lines) {
    auto tangent = line_tangent(positions, line);
    auto length  = line_length(positions, line);
    for (auto vid : line) tangents[vid] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
}

// Compute per-vertex normals for triangles.
void triangles_normals(vector<vec3f>& normals, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& triangle : triangles) {
    auto normal = triangle_normal(positions, triangle);
    auto area   = triangle_area(positions, triangle);
    for (auto vid : triangle) normals[vid] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex normals for quads.
void quads_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& quad : quads) {
    auto normal = quad_normal(positions, quad);
    auto area   = quad_area(positions, quad);
    if (!is_triangle(quad)) {
      for (auto vid : quad) normals[vid] += normal * area;
    } else {
      for (auto vid : as_triangle(quad)) normals[vid] += normal * area;
    }
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
  auto tangu = vector<vec3f>(positions.size(), vec3f{0, 0, 0});
  auto tangv = vector<vec3f>(positions.size(), vec3f{0, 0, 0});
  for (auto& triangle : triangles) {
    auto [v1, v2, v3] = triangle;
    auto tutv         = triangle_tangents_fromuv(positions[v1], positions[v2],
                positions[v3], texcoords[v1], texcoords[v2], texcoords[v3]);
    for (auto vid : triangle) tangu[vid] += normalize(tutv.first);
    for (auto vid : triangle) tangv[vid] += normalize(tutv.second);
  }
  for (auto& tangent : tangu) tangent = normalize(tangent);
  for (auto& tangent : tangv) tangent = normalize(tangent);

  auto tangent_spaces = vector<vec4f>(positions.size());
  for (auto& tangent : tangent_spaces) tangent = vec4f{0, 0, 0, 0};
  for (auto vid : range(positions.size())) {
    tangu[vid] = orthonormalize(tangu[vid], normals[vid]);
    auto s     = dot(cross(normals[vid], tangu[vid]), tangv[vid]) < 0 ? -1.0f
                                                                      : 1.0f;
    tangent_spaces[vid] = {tangu[vid], s};
  }
  return tangent_spaces;
}

// Apply skinning
pair<vector<vec3f>, vector<vec3f>> skin_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i : range(positions.size())) {
    auto [j1, j2, j3, j4] = joints[i];
    auto [w1, w2, w3, w4] = weights[i];
    skinned_positions[i]  = transform_point(xforms[j1], positions[i]) * w1 +
                           transform_point(xforms[j2], positions[i]) * w2 +
                           transform_point(xforms[j3], positions[i]) * w3 +
                           transform_point(xforms[j4], positions[i]) * w4;
    skinned_normals[i] = normalize(
        transform_direction(xforms[j1], normals[i]) * w1 +
        transform_direction(xforms[j2], normals[i]) * w2 +
        transform_direction(xforms[j3], normals[i]) * w3 +
        transform_direction(xforms[j4], normals[i]) * w4);
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning as specified in Khronos glTF
pair<vector<vec3f>, vector<vec3f>> skin_matrices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms) {
  auto skinned_positions = vector<vec3f>{positions.size()};
  auto skinned_normals   = vector<vec3f>{positions.size()};
  for (auto i : range(positions.size())) {
    auto [j1, j2, j3, j4] = joints[i];
    auto [w1, w2, w3, w4] = weights[i];
    auto xform = xforms[j1] * w1 + xforms[j2] * w2 + xforms[j3] * w3 +
                 xforms[j4] * w4;
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
  for (auto i : range(positions.size())) {
    auto [j1, j2, j3, j4] = joints[i];
    auto [w1, w2, w3, w4] = weights[i];
    skinned_positions[i]  = transform_point(xforms[j1], positions[i]) * w1 +
                           transform_point(xforms[j2], positions[i]) * w2 +
                           transform_point(xforms[j3], positions[i]) * w3 +
                           transform_point(xforms[j4], positions[i]) * w4;
    skinned_normals[i] = normalize(
        transform_direction(xforms[j1], normals[i]) * w1 +
        transform_direction(xforms[j2], normals[i]) * w2 +
        transform_direction(xforms[j3], normals[i]) * w3 +
        transform_direction(xforms[j4], normals[i]) * w4);
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
  for (auto i : range(positions.size())) {
    auto [j1, j2, j3, j4] = joints[i];
    auto [w1, w2, w3, w4] = weights[i];
    auto xform = xforms[j1] * w1 + xforms[j2] * w2 + xforms[j3] * w3 +
                 xforms[j4] * w4;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
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
  for (auto& [v1, v2, v3] : flipped) swap(v2, v3);
  return flipped;
}
vector<vec4i> flip_quads(const vector<vec4i>& quads) {
  auto flipped = quads;
  for (auto& [v1, v2, v3, v4] : flipped) {
    if (v3 != v4) {
      swap(v2, v4);
    } else {
      swap(v2, v3);
      v4 = v3;
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
  for (auto a : range(3)) {
    switch (alignment[a]) {
      case 0: break;
      case 1: offset[a] = bounds.min[a]; break;
      case 2: offset[a] = (bounds.min[a] + bounds.max[a]) / 2; break;
      case 3: offset[a] = bounds.max[a]; break;
      default: throw std::invalid_argument{"invalid alignment"};
    }
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
  for (auto& [v1, v2, v3] : triangles) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    insert_edge(emap, {v3, v1});
  }
  return emap;
}
edge_map make_edge_map(const vector<vec4i>& quads) {
  auto emap = edge_map{};
  for (auto& [v1, v2, v3, v4] : quads) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    if (v3 != v4) insert_edge(emap, {v3, v4});
    insert_edge(emap, {v4, v1});
  }
  return emap;
}
void insert_edges(edge_map& emap, const vector<vec3i>& triangles) {
  for (auto& [v1, v2, v3] : triangles) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    insert_edge(emap, {v3, v1});
  }
}
void insert_edges(edge_map& emap, const vector<vec4i>& quads) {
  for (auto& [v1, v2, v3, v4] : quads) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    if (v3 != v4) insert_edge(emap, {v3, v4});
    insert_edge(emap, {v4, v1});
  }
}
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge) {
  auto es = vec2i{min(edge), max(edge)};
  auto it = emap.edges.find(es);
  if (it == emap.edges.end()) {
    auto data = edge_map::edge_data{(int)emap.edges.size(), 1};
    emap.edges.insert(it, {es, data});
    return data.index;
  } else {
    auto& data = it->second;
    data.nfaces += 1;
    return data.index;
  }
}
// Get number of edges
int num_edges(const edge_map& emap) { return (int)emap.edges.size(); }
// Get the edge index
int edge_index(const edge_map& emap, const vec2i& edge) {
  auto es       = vec2i{min(edge), max(edge)};
  auto iterator = emap.edges.find(es);
  if (iterator == emap.edges.end()) return -1;
  return iterator->second.index;
}
// Get a list of edges, boundary edges, boundary vertices
vector<vec2i> get_edges(const edge_map& emap) {
  auto edges = vector<vec2i>(emap.edges.size());
  for (auto& [edge, data] : emap.edges) edges[data.index] = edge;
  return edges;
}
vector<vec2i> get_boundary(const edge_map& emap) {
  auto boundary = vector<vec2i>{};
  for (auto& [edge, data] : emap.edges) {
    if (data.nfaces < 2) boundary.push_back(edge);
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
  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto k = 0; k < 3; ++k) {
      auto edge = get_edge(triangles[i], k);
      auto it   = edge_map.find(edge);
      if (it == edge_map.end()) {
        edge_map.insert(it, {edge, i});
      } else {
        auto neighbor     = it->second;
        adjacencies[i][k] = neighbor;
        for (auto kk = 0; kk < 3; ++kk) {
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
    for (auto vid : range(3))
      if (v[vid] == x) return vid;
    return -1;
  };

  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<int>(triangles.size() * 3, -1);

  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto k : range(3)) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<int>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (auto i = 0; i < num_vertices; ++i) {
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
    for (auto vid : range(3))
      if (v[vid] == x) return vid;
    return -1;
  };

  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<int>(triangles.size() * 3, -1);

  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto k : range(3)) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<int>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (auto i = 0; i < num_vertices; ++i) {
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
  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto k = 0; k < 3; ++k) {
      if (adjacency[i][k] == -1)
        next_vert[triangles[i][k]] = triangles[i][(k + 1) % 3];
    }
  }

  // result
  auto boundaries = vector<vector<int>>();

  // arrange boundary vertices in loops
  for (auto i : range(next_vert.size())) {
    if (next_vert[i] == -1) continue;

    // add new empty boundary
    boundaries.emplace_back();
    auto current = (int)i;

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

// Splits a BVH node using the middle heuristic_ Returns split position and
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
  if (csize == vec3f{0, 0, 0}) return {mid, axis};

  // split along largest
  axis = (int)argmax(csize);

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

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static bvh_tree make_bvh(vector<bbox3f>& bboxes) {
  // bvh
  auto bvh = bvh_tree{};

  // prepare to build nodes
  bvh.nodes.clear();
  bvh.nodes.reserve(bboxes.size() * 2);

  // prepare primitives
  bvh.primitives.resize(bboxes.size());
  for (auto idx : range(bboxes.size())) bvh.primitives[idx] = (int)idx;

  // prepare centers
  auto centers = vector<vec3f>(bboxes.size());
  for (auto idx : range(bboxes.size())) centers[idx] = bbox_center(bboxes[idx]);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)bboxes.size()}};
  bvh.nodes.emplace_back();

  // create nodes until the queue is empty
  while (!queue.empty()) {
    // grab node to work on
    auto [nodeid, start, end] = queue.front();
    queue.pop_front();

    // grab node
    auto& node = bvh.nodes[nodeid];

    // compute bounds
    node.bbox = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox = merge(node.bbox, bboxes[bvh.primitives[i]]);

    // split into two children
    if (end - start > bvh_max_prims) {
      // get split
      auto [mid, axis] = split_middle(
          bvh.primitives, bboxes, centers, start, end);

      // make an internal node
      node.internal = true;
      node.axis     = (int8_t)axis;
      node.num      = 2;
      node.start    = (int)bvh.nodes.size();
      bvh.nodes.emplace_back();
      bvh.nodes.emplace_back();
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
  bvh.nodes.shrink_to_fit();

  // done
  return bvh;
}

// Update bvh
static void update_bvh(bvh_tree& bvh, const vector<bbox3f>& bboxes) {
  for (auto nodeid = (int)bvh.nodes.size() - 1; nodeid >= 0; nodeid--) {
    auto& node = bvh.nodes[nodeid];
    node.bbox  = invalidb3f;
    if (node.internal) {
      for (auto idx : range(2)) {
        node.bbox = merge(node.bbox, bvh.nodes[node.start + idx].bbox);
      }
    } else {
      for (auto idx : range(node.num)) {
        node.bbox = merge(node.bbox, bboxes[bvh.primitives[node.start + idx]]);
      }
    }
  }
}

// Build shape bvh
bvh_tree make_points_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx : range(bboxes.size())) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // build nodes
  return make_bvh(bboxes);
}
bvh_tree make_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2] = lines[idx];
    bboxes[idx]    = line_bounds(
        positions[v1], positions[v2], radius[v1], radius[v2]);
  }

  // build nodes
  return make_bvh(bboxes);
}
bvh_tree make_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2, v3] = triangles[idx];
    bboxes[idx] = triangle_bounds(positions[v1], positions[v2], positions[v3]);
  }

  // build nodes
  return make_bvh(bboxes);
}
bvh_tree make_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2, v3, v4] = quads[idx];
    bboxes[idx]            = quad_bounds(
        positions[v1], positions[v2], positions[v3], positions[v4]);
  }

  // build nodes
  return make_bvh(bboxes);
}

void update_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx : range(bboxes.size())) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2] = lines[idx];
    bboxes[idx]    = line_bounds(
        positions[v1], positions[v2], radius[v1], radius[v2]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2, v3] = triangles[idx];
    bboxes[idx] = triangle_bounds(positions[v1], positions[v2], positions[v3]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_quads_bvh(
    bvh_tree& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx : range(bboxes.size())) {
    auto& [v1, v2, v3, v4] = quads[idx];
    bboxes[idx]            = quad_bounds(
        positions[v1], positions[v2], positions[v3], positions[v4]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}

// Intersect ray with a bvh.
template <typename Intersect>
static shape_intersection intersect_elements_bvh(const bvh_tree& bvh,
    Intersect&& intersect_element, const ray3f& ray_, bool find_any) {
  // check empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 128>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto intersection = shape_intersection{};

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = 1 / ray.d;
  auto ray_dsign = component_less(ray_dinv, 0);

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
      for (auto idx : range(node.num)) {
        auto primitive     = bvh.primitives[node.start + idx];
        auto eintersection = intersect_element(primitive, ray);
        if (!eintersection.hit) continue;
        intersection = {
            primitive, eintersection.uv, eintersection.distance, true};
        ray.tmax = eintersection.distance;
      }
    }

    // check for early exit
    if (find_any && intersection.hit) return intersection;
  }

  return intersection;
}

// Intersect ray with a bvh.
shape_intersection intersect_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&points, &positions, &radius](int idx, const ray3f& ray) {
        auto& p = points[idx];
        return intersect_point(ray, positions[p], radius[p]);
      },
      ray, find_any);
}
shape_intersection intersect_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&lines, &positions, &radius](int idx, const ray3f& ray) {
        auto& [v1, v2] = lines[idx];
        return intersect_line(
            ray, positions[v1], positions[v2], radius[v1], radius[v2]);
      },
      ray, find_any);
}
shape_intersection intersect_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&triangles, &positions](int idx, const ray3f& ray) {
        auto& [v1, v2, v3] = triangles[idx];
        return intersect_triangle(
            ray, positions[v1], positions[v2], positions[v3]);
      },
      ray, find_any);
}
shape_intersection intersect_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&quads, &positions](int idx, const ray3f& ray) {
        auto& [v1, v2, v3, v4] = quads[idx];
        return intersect_quad(
            ray, positions[v1], positions[v2], positions[v3], positions[v4]);
      },
      ray, find_any);
}

// Intersect ray with a bvh.
template <typename Overlap>
static shape_intersection overlap_elements_bvh(const bvh_tree& bvh,
    Overlap&& overlap_element, const vec3f& pos, float max_distance,
    bool find_any) {
  // check if empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 128>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto intersection = shape_intersection{};

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
      for (auto idx : range(node.num)) {
        auto primitive     = bvh.primitives[node.start + idx];
        auto eintersection = overlap_element(primitive, pos, max_distance);
        if (!eintersection.hit) continue;
        intersection = {
            primitive, eintersection.uv, eintersection.distance, true};
        max_distance = eintersection.distance;
      }
    }

    // check for early exit
    if (find_any && intersection.hit) return intersection;
  }

  return intersection;
}

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
shape_intersection overlap_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&points, &positions, &radius](
          int idx, const vec3f& pos, float max_distance) {
        auto& p = points[idx];
        return overlap_point(pos, max_distance, positions[p], radius[p]);
      },
      pos, max_distance, find_any);
}
shape_intersection overlap_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&lines, &positions, &radius](
          int idx, const vec3f& pos, float max_distance) {
        auto& [v1, v2] = lines[idx];
        return overlap_line(pos, max_distance, positions[v1], positions[v2],
            radius[v1], radius[v2]);
      },
      pos, max_distance, find_any);
}
shape_intersection overlap_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&triangles, &positions, &radius](
          int idx, const vec3f& pos, float max_distance) {
        auto& [v1, v2, v3] = triangles[idx];
        return overlap_triangle(pos, max_distance, positions[v1], positions[v2],
            positions[v3], radius[v1], radius[v2], radius[v3]);
      },
      pos, max_distance, find_any);
}
shape_intersection overlap_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&quads, &positions, &radius](
          int idx, const vec3f& pos, float max_distance) {
        auto& [v1, v2, v3, v4] = quads[idx];
        return overlap_quad(pos, max_distance, positions[v1], positions[v2],
            positions[v3], positions[v4], radius[v1], radius[v2], radius[v3],
            radius[v4]);
      },
      pos, max_distance, find_any);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------

namespace yocto {

// Gets the cell index
vec3i get_cell_index(const hash_grid& grid, const vec3f& position) {
  auto scaledpos = position * grid.cell_inv_size;
  return (vec3i)scaledpos;
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
  auto vid  = (int)grid.positions.size();
  auto cell = get_cell_index(grid, position);
  grid.cells[cell].push_back(vid);
  grid.positions.push_back(position);
  return vid;
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
        for (auto vid : ncell_vertices) {
          if (distance2(grid.positions[vid], position) > max_radius_squared)
            continue;
          if (vid == skip_id) continue;
          neighbors.push_back(vid);
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
  for (auto& [v1, v2, v3, v4] : quads) {
    triangles.push_back({v1, v2, v4});
    if (v3 != v4) triangles.push_back({v3, v4, v2});
  }
  return triangles;
}

// Convert triangles to quads by creating degenerate quads
vector<vec4i> triangles_to_quads(const vector<vec3i>& triangles) {
  auto quads = vector<vec4i>{};
  quads.reserve(triangles.size());
  for (auto& [v1, v2, v3] : triangles) quads.push_back({v1, v2, v3, v3});
  return quads;
}

// Convert beziers to lines using 3 lines for each bezier.
vector<vec2i> bezier_to_lines(const vector<vec4i>& beziers) {
  auto lines = vector<vec2i>{};
  lines.reserve(beziers.size() * 3);
  for (auto [v1, v2, v3, v4] : beziers) {
    lines.push_back({v1, v2});
    lines.push_back({v2, v3});
    lines.push_back({v3, v4});
  }
  return lines;
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
  for (auto fid : range(quadspos.size())) {
    for (auto c : range(4)) {
      auto v = vec3i{
          quadspos[fid][c],
          (!quadsnorm.empty()) ? quadsnorm[fid][c] : -1,
          (!quadstexcoord.empty()) ? quadstexcoord[fid][c] : -1,
      };
      auto it = vert_map.find(v);
      if (it == vert_map.end()) {
        auto s = (int)vert_map.size();
        vert_map.insert(it, {v, s});
        split_quads[fid][c] = s;
      } else {
        split_quads[fid][c] = it->second;
      }
    }
  }

  // fill vert data
  split_positions.clear();
  if (!positions.empty()) {
    split_positions.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      auto [vid, _, __]      = vert;
      split_positions[index] = positions[vid];
    }
  }
  split_normals.clear();
  if (!normals.empty()) {
    split_normals.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      auto [_, vid, __]    = vert;
      split_normals[index] = normals[vid];
    }
  }
  split_texcoords.clear();
  if (!texcoords.empty()) {
    split_texcoords.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      auto [_, __, vid]      = vert;
      split_texcoords[index] = texcoords[vid];
    }
  }
}

// Weld vertices within a threshold.
pair<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold) {
  auto indices   = vector<int>(positions.size());
  auto welded    = vector<vec3f>{};
  auto grid      = make_hash_grid(threshold);
  auto neighbors = vector<int>{};
  for (auto vertex : range(positions.size())) {
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
  for (auto& triangle : wtriangles) {
    auto [v1, v2, v3] = triangle;
    triangle          = {indices[v1], indices[v2], indices[v3]};
  }
  return {wtriangles, wpositions};
}
pair<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold) {
  auto [wpositions, indices] = weld_vertices(positions, threshold);
  auto wquads                = quads;
  for (auto& quad : wquads) {
    auto [v1, v2, v3, v4] = quad;
    quad = {indices[v1], indices[v2], indices[v3], indices[v4]};
  }
  return {wquads, wpositions};
}

// Merge shape elements
void merge_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& tangents, vector<vec2f>& texcoords, vector<float>& radius,
    const vector<vec2i>& merge_lines, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_tangents,
    const vector<vec2f>& merge_texturecoords,
    const vector<float>& merge_radius) {
  auto merge_verts = (int)positions.size();
  for (auto& [v1, v2] : merge_lines)
    lines.push_back({v1 + merge_verts, v2 + merge_verts});
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
  for (auto& triangle : merge_triangles)
    triangles.push_back(triangle + merge_verts);
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
  for (auto& quad : merge_quads) quads.push_back(quad + merge_verts);
  positions.insert(
      positions.end(), merge_positions.begin(), merge_positions.end());
  normals.insert(normals.end(), merge_normals.begin(), merge_normals.end());
  texcoords.insert(
      texcoords.end(), merge_texturecoords.begin(), merge_texturecoords.end());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines.
template <typename T>
static pair<vector<vec2i>, vector<T>> subdivide_lines_impl(
    const vector<vec2i>& lines, const vector<T>& vertices) {
  // early exit
  if (lines.empty() || vertices.empty()) return {lines, vertices};
  // create vertices
  auto tvertices = vector<T>{};
  tvertices.reserve(vertices.size() + lines.size());
  for (auto& vertex : vertices) tvertices.push_back(vertex);
  for (auto& [v1, v2] : lines) {
    tvertices.push_back((vertices[v1] + vertices[v2]) / 2);
  }
  // create lines
  auto tlines = vector<vec2i>{};
  tlines.reserve(lines.size() * 2);
  auto line_vertex = [nverts = (int)vertices.size()](
                         size_t line_id) { return nverts + (int)line_id; };
  for (auto&& [line_id, line] : enumerate(lines)) {
    auto& [v1, v2] = line;
    tlines.push_back({v1, line_vertex(line_id)});
    tlines.push_back({line_vertex(line_id), v2});
  }
  // done
  return {tlines, tvertices};
}

// Subdivide triangle.
template <typename T>
static pair<vector<vec3i>, vector<T>> subdivide_triangles_impl(
    const vector<vec3i>& triangles, const vector<T>& vertices) {
  // early exit
  if (triangles.empty() || vertices.empty()) return {triangles, vertices};
  // get edges
  auto emap  = make_edge_map(triangles);
  auto edges = get_edges(emap);
  // create vertices
  auto tvertices = vector<T>{};
  tvertices.reserve(vertices.size() + edges.size());
  for (auto& vertex : vertices) tvertices.push_back(vertex);
  for (auto& [v1, v2] : edges)
    tvertices.push_back((vertices[v1] + vertices[v2]) / 2);
  // create triangles
  auto ttriangles = vector<vec3i>{};
  ttriangles.reserve(triangles.size() * 4);
  auto edge_vertex = [&emap, nverts = (int)vertices.size()](const vec2i& edge) {
    return nverts + edge_index(emap, edge);
  };
  for (auto& triangle : triangles) {
    auto [v1, v2, v3] = triangle;
    ttriangles.push_back({v1, edge_vertex({v1, v2}), edge_vertex({v3, v1})});
    ttriangles.push_back({v2, edge_vertex({v2, v3}), edge_vertex({v1, v2})});
    ttriangles.push_back({v3, edge_vertex({v3, v1}), edge_vertex({v2, v3})});
    ttriangles.push_back(
        {edge_vertex({v1, v2}), edge_vertex({v2, v3}), edge_vertex({v3, v1})});
  }
  // done
  return {ttriangles, tvertices};
}

// Subdivide quads.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_quads_impl(
    const vector<vec4i>& quads, const vector<T>& vertices) {
  // early exit
  if (quads.empty() || vertices.empty()) return {quads, vertices};
  // get edges
  auto emap  = make_edge_map(quads);
  auto edges = get_edges(emap);
  // create vertices
  auto tvertices = vector<T>{};
  tvertices.reserve(vertices.size() + edges.size() + quads.size());
  for (auto& vertex : vertices) tvertices.push_back(vertex);
  for (auto& [v1, v2] : edges)
    tvertices.push_back((vertices[v1] + vertices[v2]) / 2);
  for (auto& [v1, v2, v3, v4] : quads) {
    if (v3 != v4) {
      tvertices.push_back(
          (vertices[v1] + vertices[v2] + vertices[v3] + vertices[v4]) / 4);
    } else {
      tvertices.push_back((vertices[v1] + vertices[v2] + vertices[v3]) / 3);
    }
  }
  // create quads
  auto tquads = vector<vec4i>{};
  tquads.reserve(quads.size() * 4);
  auto edge_vertex = [&emap, nverts = (int)vertices.size()](const vec2i& edge) {
    return nverts + edge_index(emap, edge);
  };
  auto quad_vertex = [nverts    = (int)vertices.size(),
                         nedges = (int)edges.size()](size_t quad_id) {
    return nverts + nedges + (int)quad_id;
  };
  for (auto&& [quad_id, quad] : enumerate(quads)) {
    auto [v1, v2, v3, v4] = quad;
    if (v3 != v4) {
      tquads.push_back({v1, edge_vertex({v1, v2}), quad_vertex(quad_id),
          edge_vertex({v4, v1})});
      tquads.push_back({v2, edge_vertex({v2, v3}), quad_vertex(quad_id),
          edge_vertex({v1, v2})});
      tquads.push_back({v3, edge_vertex({v3, v4}), quad_vertex(quad_id),
          edge_vertex({v2, v3})});
      tquads.push_back({v4, edge_vertex({v4, v1}), quad_vertex(quad_id),
          edge_vertex({v3, v4})});
    } else {
      tquads.push_back({v1, edge_vertex({v1, v2}), quad_vertex(quad_id),
          edge_vertex({v3, v1})});
      tquads.push_back({v2, edge_vertex({v2, v3}), quad_vertex(quad_id),
          edge_vertex({v1, v2})});
      tquads.push_back({v3, edge_vertex({v3, v1}), quad_vertex(quad_id),
          edge_vertex({v2, v3})});
    }
  }
  // done
  return {tquads, tvertices};
}

// Subdivide beziers.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_beziers_impl(
    const vector<vec4i>& beziers, const vector<T>& vertices) {
  // early exit
  if (beziers.empty() || vertices.empty()) return {beziers, vertices};
  // get edges
  auto vmap      = unordered_map<int, int>();
  auto tvertices = vector<T>();
  auto tbeziers  = vector<vec4i>();
  for (auto& bezier : beziers) {
    auto [v1, v2, v3, v4] = bezier;
    for (auto vid : {v1, v4}) {
      if (vmap.find(vid) == vmap.end()) {
        vmap[vid] = (int)tvertices.size();
        tvertices.push_back(vertices[vid]);
      }
    }
    auto bo = (int)tvertices.size();
    tbeziers.push_back({vmap.at(v1), bo + 0, bo + 1, bo + 2});
    tbeziers.push_back({bo + 2, bo + 3, bo + 4, vmap.at(v4)});
    tvertices.push_back(vertices[v1] / 2 + vertices[v2] / 2);
    tvertices.push_back(vertices[v1] / 4 + vertices[v2] / 2 + vertices[v3] / 4);
    tvertices.push_back(
        vertices[v1] / 8 + vertices[v2] * ((float)3 / (float)8) +
        vertices[v3] * ((float)3 / (float)8) + vertices[v4] / 8);
    tvertices.push_back(vertices[v2] / 4 + vertices[v3] / 2 + vertices[v4] / 4);
    tvertices.push_back(vertices[v3] / 2 + vertices[v4] / 2);
  }

  // done
  return {tbeziers, tvertices};
}

// Subdivide catmullclark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_impl(
    const vector<vec4i>& quads, const vector<T>& vertices, bool lock_boundary) {
  // early exit
  if (quads.empty() || vertices.empty()) return {quads, vertices};
  // get edges
  auto emap     = make_edge_map(quads);
  auto edges    = get_edges(emap);
  auto boundary = get_boundary(emap);

  // split elements ------------------------------------
  // create vertices
  auto tvertices = vector<T>{};
  tvertices.reserve(vertices.size() + edges.size() + quads.size());
  for (auto& vertex : vertices) tvertices.push_back(vertex);
  for (auto& [v1, v2] : edges)
    tvertices.push_back((vertices[v1] + vertices[v2]) / 2);
  for (auto& [v1, v2, v3, v4] : quads) {
    if (v3 != v4) {
      tvertices.push_back(
          (vertices[v1] + vertices[v2] + vertices[v3] + vertices[v4]) / 4);
    } else {
      tvertices.push_back((vertices[v1] + vertices[v2] + vertices[v3]) / 3);
    }
  }
  // create quads
  auto tquads = vector<vec4i>{};
  tquads.reserve(quads.size() * 4);
  auto edge_vertex = [&emap, nverts = (int)vertices.size()](const vec2i& edge) {
    return nverts + edge_index(emap, edge);
  };
  auto quad_vertex = [nverts    = (int)vertices.size(),
                         nedges = (int)edges.size()](size_t quad_id) {
    return nverts + nedges + (int)quad_id;
  };
  for (auto&& [quad_id, quad] : enumerate(quads)) {
    auto& [v1, v2, v3, v4] = quad;
    if (v3 != v4) {
      tquads.push_back({v1, edge_vertex({v1, v2}), quad_vertex(quad_id),
          edge_vertex({v4, v1})});
      tquads.push_back({v2, edge_vertex({v2, v3}), quad_vertex(quad_id),
          edge_vertex({v1, v2})});
      tquads.push_back({v3, edge_vertex({v3, v4}), quad_vertex(quad_id),
          edge_vertex({v2, v3})});
      tquads.push_back({v4, edge_vertex({v4, v1}), quad_vertex(quad_id),
          edge_vertex({v3, v4})});
    } else {
      tquads.push_back({v1, edge_vertex({v1, v2}), quad_vertex(quad_id),
          edge_vertex({v3, v1})});
      tquads.push_back({v2, edge_vertex({v2, v3}), quad_vertex(quad_id),
          edge_vertex({v1, v2})});
      tquads.push_back({v3, edge_vertex({v3, v1}), quad_vertex(quad_id),
          edge_vertex({v2, v3})});
    }
  }

  // split boundary
  auto tboundary = vector<vec2i>{};
  tboundary.reserve(boundary.size());
  for (auto& [v1, v2] : boundary) {
    tboundary.push_back({v1, edge_vertex({v1, v2})});
    tboundary.push_back({edge_vertex({v1, v2}), v2});
  }

  // setup creases -----------------------------------
  auto tcrease_edges = vector<vec2i>{};
  auto tcrease_verts = vector<int>{};
  if (lock_boundary) {
    for (auto& [v1, v2] : tboundary) {
      tcrease_verts.push_back(v1);
      tcrease_verts.push_back(v2);
    }
  } else {
    for (auto& b : tboundary) tcrease_edges.push_back(b);
  }

  // define vertices valence ---------------------------
  auto tvert_val = vector<int>(tvertices.size(), 2);
  for (auto& [v1, v2] : tboundary) {
    tvert_val[v1] = (lock_boundary) ? 0 : 1;
    tvert_val[v2] = (lock_boundary) ? 0 : 1;
  }

  // averaging pass ----------------------------------
  auto avert  = vector<T>(tvertices.size(), T());
  auto acount = vector<int>(tvertices.size(), 0);
  for (auto& point : tcrease_verts) {
    if (tvert_val[point] != 0) continue;
    avert[point] += tvertices[point];
    acount[point] += 1;
  }
  for (auto& edge : tcrease_edges) {
    auto [v1, v2] = edge;
    auto centroid = (tvertices[v1] + tvertices[v2]) / 2;
    for (auto vid : edge) {
      if (tvert_val[vid] != 1) continue;
      avert[vid] += centroid;
      acount[vid] += 1;
    }
  }
  for (auto& quad : tquads) {
    auto [v1, v2, v3, v4] = quad;
    auto centroid =
        (tvertices[v1] + tvertices[v2] + tvertices[v3] + tvertices[v4]) / 4;
    for (auto vid : quad) {
      if (tvert_val[vid] != 2) continue;
      avert[vid] += centroid;
      acount[vid] += 1;
    }
  }
  for (auto i : range(tvertices.size())) avert[i] /= (float)acount[i];

  // correction pass ----------------------------------
  // p = p + (avg_p - p) * (4/avg_count)
  for (auto i : range(tvertices.size())) {
    if (tvert_val[i] != 2) continue;
    avert[i] = tvertices[i] +
               (avert[i] - tvertices[i]) * (4 / (float)acount[i]);
  }
  tvertices = avert;

  // done
  return {tquads, tvertices};
}

pair<vector<vec2i>, vector<float>> subdivide_lines(
    const vector<vec2i>& lines, const vector<float>& vertices) {
  return subdivide_lines_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec2f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec2f>& vertices) {
  return subdivide_lines_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec3f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec3f>& vertices) {
  return subdivide_lines_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec4f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec4f>& vertices) {
  return subdivide_lines_impl(lines, vertices);
}

pair<vector<vec3i>, vector<float>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<float>& vertices) {
  return subdivide_triangles_impl(triangles, vertices);
}
pair<vector<vec3i>, vector<vec2f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec2f>& vertices) {
  return subdivide_triangles_impl(triangles, vertices);
}
pair<vector<vec3i>, vector<vec3f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& vertices) {
  return subdivide_triangles_impl(triangles, vertices);
}
pair<vector<vec3i>, vector<vec4f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec4f>& vertices) {
  return subdivide_triangles_impl(triangles, vertices);
}

pair<vector<vec4i>, vector<float>> subdivide_quads(
    const vector<vec4i>& quads, const vector<float>& vertices) {
  return subdivide_quads_impl(quads, vertices);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec2f>& vertices) {
  return subdivide_quads_impl(quads, vertices);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec3f>& vertices) {
  return subdivide_quads_impl(quads, vertices);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec4f>& vertices) {
  return subdivide_quads_impl(quads, vertices);
}

pair<vector<vec4i>, vector<float>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<float>& vertices) {
  return subdivide_beziers_impl(beziers, vertices);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec2f>& vertices) {
  return subdivide_beziers_impl(beziers, vertices);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec3f>& vertices) {
  return subdivide_beziers_impl(beziers, vertices);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec4f>& vertices) {
  return subdivide_beziers_impl(beziers, vertices);
}

pair<vector<vec4i>, vector<float>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<float>& vertices,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vertices, lock_boundary);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec2f>& vertices,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vertices, lock_boundary);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec3f>& vertices,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vertices, lock_boundary);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec4f>& vertices,
    bool lock_boundary) {
  return subdivide_catmullclark_impl(quads, vertices, lock_boundary);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
int sample_points(int npoints, float re) { return sample_uniform(npoints, re); }
int sample_points(const vector<float>& cdf, float re) {
  return sample_discrete(cdf, re);
}
vector<float> sample_points_cdf(int npoints) {
  auto cdf = vector<float>(npoints);
  for (auto i : range(cdf.size())) cdf[i] = 1 + (i != 0 ? cdf[i - 1] : 0);
  return cdf;
}
void sample_points_cdf(vector<float>& cdf, int npoints) {
  for (auto i : range(cdf.size())) cdf[i] = 1 + (i != 0 ? cdf[i - 1] : 0);
}

// Pick a point on lines uniformly.
pair<int, float> sample_lines(const vector<float>& cdf, float re, float ru) {
  return {sample_discrete(cdf, re), ru};
}
vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto cdf = vector<float>(lines.size());
  for (auto i : range(cdf.size())) {
    auto& [v1, v2] = lines[i];
    auto w         = line_length(positions[v1], positions[v2]);
    cdf[i]         = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_lines_cdf(vector<float>& cdf, const vector<vec2i>& lines,
    const vector<vec3f>& positions) {
  for (auto i : range(cdf.size())) {
    auto& [v1, v2] = lines[i];
    auto w         = line_length(positions[v1], positions[v2]);
    cdf[i]         = w + (i != 0 ? cdf[i - 1] : 0);
  }
}

// Pick a point on a triangle mesh uniformly.
pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv) {
  return {sample_discrete(cdf, re), sample_triangle(ruv)};
}
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto cdf = vector<float>(triangles.size());
  for (auto i : range(cdf.size())) {
    auto& [v1, v2, v3] = triangles[i];
    auto w = triangle_area(positions[v1], positions[v2], positions[v3]);
    cdf[i] = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_triangles_cdf(vector<float>& cdf, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  for (auto i : range(cdf.size())) {
    auto& [v1, v2, v3] = triangles[i];
    auto w = triangle_area(positions[v1], positions[v2], positions[v3]);
    cdf[i] = w + (i != 0 ? cdf[i - 1] : 0);
  }
}

// Pick a point on a quad mesh uniformly.
pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv) {
  return {sample_discrete(cdf, re), ruv};
}
pair<int, vec2f> sample_quads(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv) {
  auto element = sample_discrete(cdf, re);
  if (!is_triangle(quads[element])) {
    return {element, sample_triangle(ruv)};
  } else {
    return {element, ruv};
  }
}
vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto cdf = vector<float>(quads.size());
  for (auto i : range(cdf.size())) {
    auto& [v1, v2, v3, v4] = quads[i];
    auto w                 = quad_area(
        positions[v1], positions[v2], positions[v3], positions[v4]);
    cdf[i] = w + (i ? cdf[i - 1] : 0);
  }
  return cdf;
}
void sample_quads_cdf(vector<float>& cdf, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  for (auto i : range(cdf.size())) {
    auto& [v1, v2, v3, v4] = quads[i];
    auto w                 = quad_area(
        positions[v1], positions[v2], positions[v3], positions[v4]);
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
  for (auto i : range(npoints)) {
    auto sample          = sample_triangles(cdf, rand1f(rng), rand2f(rng));
    auto& [v1, v2, v3]   = triangles[sample.first];
    auto uv              = sample.second;
    sampled_positions[i] = interpolate_triangle(
        positions[v1], positions[v2], positions[v3], uv);
    if (!sampled_normals.empty()) {
      sampled_normals[i] = normalize(
          interpolate_triangle(normals[v1], normals[v2], normals[v3], uv));
    } else {
      sampled_normals[i] = triangle_normal(
          positions[v1], positions[v2], positions[v3]);
    }
    if (!sampled_texcoords.empty()) {
      sampled_texcoords[i] = interpolate_triangle(
          texcoords[v1], texcoords[v2], texcoords[v3], uv);
    } else {
      sampled_texcoords[i] = vec2f{0, 0};
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
  for (auto i : range(npoints)) {
    auto sample            = sample_quads(cdf, rand1f(rng), rand2f(rng));
    auto& [v1, v2, v3, v4] = quads[sample.first];
    auto uv                = sample.second;
    sampled_positions[i]   = interpolate_quad(
        positions[v1], positions[v2], positions[v3], positions[v4], uv);
    if (!sampled_normals.empty()) {
      sampled_normals[i] = normalize(interpolate_quad(
          normals[v1], normals[v2], normals[v3], normals[v4], uv));
    } else {
      sampled_normals[i] = quad_normal(
          positions[v1], positions[v2], positions[v3], positions[v4]);
    }
    if (!sampled_texcoords.empty()) {
      sampled_texcoords[i] = interpolate_quad(
          texcoords[v1], texcoords[v2], texcoords[v3], texcoords[v4], uv);
    } else {
      sampled_texcoords[i] = vec2f{0, 0};
    }
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DATA
// -----------------------------------------------------------------------------
namespace yocto {

vector<vec3f> suzanne_positions = vector<vec3f>{{0.4375, 0.1640625, 0.765625},
    {-0.4375, 0.1640625, 0.765625}, {0.5, 0.09375, 0.6875},
    {-0.5, 0.09375, 0.6875}, {0.546875, 0.0546875, 0.578125},
    {-0.546875, 0.0546875, 0.578125}, {0.3515625, -0.0234375, 0.6171875},
    {-0.3515625, -0.0234375, 0.6171875}, {0.3515625, 0.03125, 0.71875},
    {-0.3515625, 0.03125, 0.71875}, {0.3515625, 0.1328125, 0.78125},
    {-0.3515625, 0.1328125, 0.78125}, {0.2734375, 0.1640625, 0.796875},
    {-0.2734375, 0.1640625, 0.796875}, {0.203125, 0.09375, 0.7421875},
    {-0.203125, 0.09375, 0.7421875}, {0.15625, 0.0546875, 0.6484375},
    {-0.15625, 0.0546875, 0.6484375}, {0.078125, 0.2421875, 0.65625},
    {-0.078125, 0.2421875, 0.65625}, {0.140625, 0.2421875, 0.7421875},
    {-0.140625, 0.2421875, 0.7421875}, {0.2421875, 0.2421875, 0.796875},
    {-0.2421875, 0.2421875, 0.796875}, {0.2734375, 0.328125, 0.796875},
    {-0.2734375, 0.328125, 0.796875}, {0.203125, 0.390625, 0.7421875},
    {-0.203125, 0.390625, 0.7421875}, {0.15625, 0.4375, 0.6484375},
    {-0.15625, 0.4375, 0.6484375}, {0.3515625, 0.515625, 0.6171875},
    {-0.3515625, 0.515625, 0.6171875}, {0.3515625, 0.453125, 0.71875},
    {-0.3515625, 0.453125, 0.71875}, {0.3515625, 0.359375, 0.78125},
    {-0.3515625, 0.359375, 0.78125}, {0.4375, 0.328125, 0.765625},
    {-0.4375, 0.328125, 0.765625}, {0.5, 0.390625, 0.6875},
    {-0.5, 0.390625, 0.6875}, {0.546875, 0.4375, 0.578125},
    {-0.546875, 0.4375, 0.578125}, {0.625, 0.2421875, 0.5625},
    {-0.625, 0.2421875, 0.5625}, {0.5625, 0.2421875, 0.671875},
    {-0.5625, 0.2421875, 0.671875}, {0.46875, 0.2421875, 0.7578125},
    {-0.46875, 0.2421875, 0.7578125}, {0.4765625, 0.2421875, 0.7734375},
    {-0.4765625, 0.2421875, 0.7734375}, {0.4453125, 0.3359375, 0.78125},
    {-0.4453125, 0.3359375, 0.78125}, {0.3515625, 0.375, 0.8046875},
    {-0.3515625, 0.375, 0.8046875}, {0.265625, 0.3359375, 0.8203125},
    {-0.265625, 0.3359375, 0.8203125}, {0.2265625, 0.2421875, 0.8203125},
    {-0.2265625, 0.2421875, 0.8203125}, {0.265625, 0.15625, 0.8203125},
    {-0.265625, 0.15625, 0.8203125}, {0.3515625, 0.2421875, 0.828125},
    {-0.3515625, 0.2421875, 0.828125}, {0.3515625, 0.1171875, 0.8046875},
    {-0.3515625, 0.1171875, 0.8046875}, {0.4453125, 0.15625, 0.78125},
    {-0.4453125, 0.15625, 0.78125}, {0.0, 0.4296875, 0.7421875},
    {0.0, 0.3515625, 0.8203125}, {0.0, -0.6796875, 0.734375},
    {0.0, -0.3203125, 0.78125}, {0.0, -0.1875, 0.796875},
    {0.0, -0.7734375, 0.71875}, {0.0, 0.40625, 0.6015625},
    {0.0, 0.5703125, 0.5703125}, {0.0, 0.8984375, -0.546875},
    {0.0, 0.5625, -0.8515625}, {0.0, 0.0703125, -0.828125},
    {0.0, -0.3828125, -0.3515625}, {0.203125, -0.1875, 0.5625},
    {-0.203125, -0.1875, 0.5625}, {0.3125, -0.4375, 0.5703125},
    {-0.3125, -0.4375, 0.5703125}, {0.3515625, -0.6953125, 0.5703125},
    {-0.3515625, -0.6953125, 0.5703125}, {0.3671875, -0.890625, 0.53125},
    {-0.3671875, -0.890625, 0.53125}, {0.328125, -0.9453125, 0.5234375},
    {-0.328125, -0.9453125, 0.5234375}, {0.1796875, -0.96875, 0.5546875},
    {-0.1796875, -0.96875, 0.5546875}, {0.0, -0.984375, 0.578125},
    {0.4375, -0.140625, 0.53125}, {-0.4375, -0.140625, 0.53125},
    {0.6328125, -0.0390625, 0.5390625}, {-0.6328125, -0.0390625, 0.5390625},
    {0.828125, 0.1484375, 0.4453125}, {-0.828125, 0.1484375, 0.4453125},
    {0.859375, 0.4296875, 0.59375}, {-0.859375, 0.4296875, 0.59375},
    {0.7109375, 0.484375, 0.625}, {-0.7109375, 0.484375, 0.625},
    {0.4921875, 0.6015625, 0.6875}, {-0.4921875, 0.6015625, 0.6875},
    {0.3203125, 0.7578125, 0.734375}, {-0.3203125, 0.7578125, 0.734375},
    {0.15625, 0.71875, 0.7578125}, {-0.15625, 0.71875, 0.7578125},
    {0.0625, 0.4921875, 0.75}, {-0.0625, 0.4921875, 0.75},
    {0.1640625, 0.4140625, 0.7734375}, {-0.1640625, 0.4140625, 0.7734375},
    {0.125, 0.3046875, 0.765625}, {-0.125, 0.3046875, 0.765625},
    {0.203125, 0.09375, 0.7421875}, {-0.203125, 0.09375, 0.7421875},
    {0.375, 0.015625, 0.703125}, {-0.375, 0.015625, 0.703125},
    {0.4921875, 0.0625, 0.671875}, {-0.4921875, 0.0625, 0.671875},
    {0.625, 0.1875, 0.6484375}, {-0.625, 0.1875, 0.6484375},
    {0.640625, 0.296875, 0.6484375}, {-0.640625, 0.296875, 0.6484375},
    {0.6015625, 0.375, 0.6640625}, {-0.6015625, 0.375, 0.6640625},
    {0.4296875, 0.4375, 0.71875}, {-0.4296875, 0.4375, 0.71875},
    {0.25, 0.46875, 0.7578125}, {-0.25, 0.46875, 0.7578125},
    {0.0, -0.765625, 0.734375}, {0.109375, -0.71875, 0.734375},
    {-0.109375, -0.71875, 0.734375}, {0.1171875, -0.8359375, 0.7109375},
    {-0.1171875, -0.8359375, 0.7109375}, {0.0625, -0.8828125, 0.6953125},
    {-0.0625, -0.8828125, 0.6953125}, {0.0, -0.890625, 0.6875},
    {0.0, -0.1953125, 0.75}, {0.0, -0.140625, 0.7421875},
    {0.1015625, -0.1484375, 0.7421875}, {-0.1015625, -0.1484375, 0.7421875},
    {0.125, -0.2265625, 0.75}, {-0.125, -0.2265625, 0.75},
    {0.0859375, -0.2890625, 0.7421875}, {-0.0859375, -0.2890625, 0.7421875},
    {0.3984375, -0.046875, 0.671875}, {-0.3984375, -0.046875, 0.671875},
    {0.6171875, 0.0546875, 0.625}, {-0.6171875, 0.0546875, 0.625},
    {0.7265625, 0.203125, 0.6015625}, {-0.7265625, 0.203125, 0.6015625},
    {0.7421875, 0.375, 0.65625}, {-0.7421875, 0.375, 0.65625},
    {0.6875, 0.4140625, 0.7265625}, {-0.6875, 0.4140625, 0.7265625},
    {0.4375, 0.546875, 0.796875}, {-0.4375, 0.546875, 0.796875},
    {0.3125, 0.640625, 0.8359375}, {-0.3125, 0.640625, 0.8359375},
    {0.203125, 0.6171875, 0.8515625}, {-0.203125, 0.6171875, 0.8515625},
    {0.1015625, 0.4296875, 0.84375}, {-0.1015625, 0.4296875, 0.84375},
    {0.125, -0.1015625, 0.8125}, {-0.125, -0.1015625, 0.8125},
    {0.2109375, -0.4453125, 0.7109375}, {-0.2109375, -0.4453125, 0.7109375},
    {0.25, -0.703125, 0.6875}, {-0.25, -0.703125, 0.6875},
    {0.265625, -0.8203125, 0.6640625}, {-0.265625, -0.8203125, 0.6640625},
    {0.234375, -0.9140625, 0.6328125}, {-0.234375, -0.9140625, 0.6328125},
    {0.1640625, -0.9296875, 0.6328125}, {-0.1640625, -0.9296875, 0.6328125},
    {0.0, -0.9453125, 0.640625}, {0.0, 0.046875, 0.7265625},
    {0.0, 0.2109375, 0.765625}, {0.328125, 0.4765625, 0.7421875},
    {-0.328125, 0.4765625, 0.7421875}, {0.1640625, 0.140625, 0.75},
    {-0.1640625, 0.140625, 0.75}, {0.1328125, 0.2109375, 0.7578125},
    {-0.1328125, 0.2109375, 0.7578125}, {0.1171875, -0.6875, 0.734375},
    {-0.1171875, -0.6875, 0.734375}, {0.078125, -0.4453125, 0.75},
    {-0.078125, -0.4453125, 0.75}, {0.0, -0.4453125, 0.75},
    {0.0, -0.328125, 0.7421875}, {0.09375, -0.2734375, 0.78125},
    {-0.09375, -0.2734375, 0.78125}, {0.1328125, -0.2265625, 0.796875},
    {-0.1328125, -0.2265625, 0.796875}, {0.109375, -0.1328125, 0.78125},
    {-0.109375, -0.1328125, 0.78125}, {0.0390625, -0.125, 0.78125},
    {-0.0390625, -0.125, 0.78125}, {0.0, -0.203125, 0.828125},
    {0.046875, -0.1484375, 0.8125}, {-0.046875, -0.1484375, 0.8125},
    {0.09375, -0.15625, 0.8125}, {-0.09375, -0.15625, 0.8125},
    {0.109375, -0.2265625, 0.828125}, {-0.109375, -0.2265625, 0.828125},
    {0.078125, -0.25, 0.8046875}, {-0.078125, -0.25, 0.8046875},
    {0.0, -0.2890625, 0.8046875}, {0.2578125, -0.3125, 0.5546875},
    {-0.2578125, -0.3125, 0.5546875}, {0.1640625, -0.2421875, 0.7109375},
    {-0.1640625, -0.2421875, 0.7109375}, {0.1796875, -0.3125, 0.7109375},
    {-0.1796875, -0.3125, 0.7109375}, {0.234375, -0.25, 0.5546875},
    {-0.234375, -0.25, 0.5546875}, {0.0, -0.875, 0.6875},
    {0.046875, -0.8671875, 0.6875}, {-0.046875, -0.8671875, 0.6875},
    {0.09375, -0.8203125, 0.7109375}, {-0.09375, -0.8203125, 0.7109375},
    {0.09375, -0.7421875, 0.7265625}, {-0.09375, -0.7421875, 0.7265625},
    {0.0, -0.78125, 0.65625}, {0.09375, -0.75, 0.6640625},
    {-0.09375, -0.75, 0.6640625}, {0.09375, -0.8125, 0.640625},
    {-0.09375, -0.8125, 0.640625}, {0.046875, -0.8515625, 0.6328125},
    {-0.046875, -0.8515625, 0.6328125}, {0.0, -0.859375, 0.6328125},
    {0.171875, 0.21875, 0.78125}, {-0.171875, 0.21875, 0.78125},
    {0.1875, 0.15625, 0.7734375}, {-0.1875, 0.15625, 0.7734375},
    {0.3359375, 0.4296875, 0.7578125}, {-0.3359375, 0.4296875, 0.7578125},
    {0.2734375, 0.421875, 0.7734375}, {-0.2734375, 0.421875, 0.7734375},
    {0.421875, 0.3984375, 0.7734375}, {-0.421875, 0.3984375, 0.7734375},
    {0.5625, 0.3515625, 0.6953125}, {-0.5625, 0.3515625, 0.6953125},
    {0.5859375, 0.2890625, 0.6875}, {-0.5859375, 0.2890625, 0.6875},
    {0.578125, 0.1953125, 0.6796875}, {-0.578125, 0.1953125, 0.6796875},
    {0.4765625, 0.1015625, 0.71875}, {-0.4765625, 0.1015625, 0.71875},
    {0.375, 0.0625, 0.7421875}, {-0.375, 0.0625, 0.7421875},
    {0.2265625, 0.109375, 0.78125}, {-0.2265625, 0.109375, 0.78125},
    {0.1796875, 0.296875, 0.78125}, {-0.1796875, 0.296875, 0.78125},
    {0.2109375, 0.375, 0.78125}, {-0.2109375, 0.375, 0.78125},
    {0.234375, 0.359375, 0.7578125}, {-0.234375, 0.359375, 0.7578125},
    {0.1953125, 0.296875, 0.7578125}, {-0.1953125, 0.296875, 0.7578125},
    {0.2421875, 0.125, 0.7578125}, {-0.2421875, 0.125, 0.7578125},
    {0.375, 0.0859375, 0.7265625}, {-0.375, 0.0859375, 0.7265625},
    {0.4609375, 0.1171875, 0.703125}, {-0.4609375, 0.1171875, 0.703125},
    {0.546875, 0.2109375, 0.671875}, {-0.546875, 0.2109375, 0.671875},
    {0.5546875, 0.28125, 0.671875}, {-0.5546875, 0.28125, 0.671875},
    {0.53125, 0.3359375, 0.6796875}, {-0.53125, 0.3359375, 0.6796875},
    {0.4140625, 0.390625, 0.75}, {-0.4140625, 0.390625, 0.75},
    {0.28125, 0.3984375, 0.765625}, {-0.28125, 0.3984375, 0.765625},
    {0.3359375, 0.40625, 0.75}, {-0.3359375, 0.40625, 0.75},
    {0.203125, 0.171875, 0.75}, {-0.203125, 0.171875, 0.75},
    {0.1953125, 0.2265625, 0.75}, {-0.1953125, 0.2265625, 0.75},
    {0.109375, 0.4609375, 0.609375}, {-0.109375, 0.4609375, 0.609375},
    {0.1953125, 0.6640625, 0.6171875}, {-0.1953125, 0.6640625, 0.6171875},
    {0.3359375, 0.6875, 0.59375}, {-0.3359375, 0.6875, 0.59375},
    {0.484375, 0.5546875, 0.5546875}, {-0.484375, 0.5546875, 0.5546875},
    {0.6796875, 0.453125, 0.4921875}, {-0.6796875, 0.453125, 0.4921875},
    {0.796875, 0.40625, 0.4609375}, {-0.796875, 0.40625, 0.4609375},
    {0.7734375, 0.1640625, 0.375}, {-0.7734375, 0.1640625, 0.375},
    {0.6015625, 0.0, 0.4140625}, {-0.6015625, 0.0, 0.4140625},
    {0.4375, -0.09375, 0.46875}, {-0.4375, -0.09375, 0.46875},
    {0.0, 0.8984375, 0.2890625}, {0.0, 0.984375, -0.078125},
    {0.0, -0.1953125, -0.671875}, {0.0, -0.4609375, 0.1875},
    {0.0, -0.9765625, 0.4609375}, {0.0, -0.8046875, 0.34375},
    {0.0, -0.5703125, 0.3203125}, {0.0, -0.484375, 0.28125},
    {0.8515625, 0.234375, 0.0546875}, {-0.8515625, 0.234375, 0.0546875},
    {0.859375, 0.3203125, -0.046875}, {-0.859375, 0.3203125, -0.046875},
    {0.7734375, 0.265625, -0.4375}, {-0.7734375, 0.265625, -0.4375},
    {0.4609375, 0.4375, -0.703125}, {-0.4609375, 0.4375, -0.703125},
    {0.734375, -0.046875, 0.0703125}, {-0.734375, -0.046875, 0.0703125},
    {0.59375, -0.125, -0.1640625}, {-0.59375, -0.125, -0.1640625},
    {0.640625, -0.0078125, -0.4296875}, {-0.640625, -0.0078125, -0.4296875},
    {0.3359375, 0.0546875, -0.6640625}, {-0.3359375, 0.0546875, -0.6640625},
    {0.234375, -0.3515625, 0.40625}, {-0.234375, -0.3515625, 0.40625},
    {0.1796875, -0.4140625, 0.2578125}, {-0.1796875, -0.4140625, 0.2578125},
    {0.2890625, -0.7109375, 0.3828125}, {-0.2890625, -0.7109375, 0.3828125},
    {0.25, -0.5, 0.390625}, {-0.25, -0.5, 0.390625},
    {0.328125, -0.9140625, 0.3984375}, {-0.328125, -0.9140625, 0.3984375},
    {0.140625, -0.7578125, 0.3671875}, {-0.140625, -0.7578125, 0.3671875},
    {0.125, -0.5390625, 0.359375}, {-0.125, -0.5390625, 0.359375},
    {0.1640625, -0.9453125, 0.4375}, {-0.1640625, -0.9453125, 0.4375},
    {0.21875, -0.28125, 0.4296875}, {-0.21875, -0.28125, 0.4296875},
    {0.2109375, -0.2265625, 0.46875}, {-0.2109375, -0.2265625, 0.46875},
    {0.203125, -0.171875, 0.5}, {-0.203125, -0.171875, 0.5},
    {0.2109375, -0.390625, 0.1640625}, {-0.2109375, -0.390625, 0.1640625},
    {0.296875, -0.3125, -0.265625}, {-0.296875, -0.3125, -0.265625},
    {0.34375, -0.1484375, -0.5390625}, {-0.34375, -0.1484375, -0.5390625},
    {0.453125, 0.8671875, -0.3828125}, {-0.453125, 0.8671875, -0.3828125},
    {0.453125, 0.9296875, -0.0703125}, {-0.453125, 0.9296875, -0.0703125},
    {0.453125, 0.8515625, 0.234375}, {-0.453125, 0.8515625, 0.234375},
    {0.4609375, 0.5234375, 0.4296875}, {-0.4609375, 0.5234375, 0.4296875},
    {0.7265625, 0.40625, 0.3359375}, {-0.7265625, 0.40625, 0.3359375},
    {0.6328125, 0.453125, 0.28125}, {-0.6328125, 0.453125, 0.28125},
    {0.640625, 0.703125, 0.0546875}, {-0.640625, 0.703125, 0.0546875},
    {0.796875, 0.5625, 0.125}, {-0.796875, 0.5625, 0.125},
    {0.796875, 0.6171875, -0.1171875}, {-0.796875, 0.6171875, -0.1171875},
    {0.640625, 0.75, -0.1953125}, {-0.640625, 0.75, -0.1953125},
    {0.640625, 0.6796875, -0.4453125}, {-0.640625, 0.6796875, -0.4453125},
    {0.796875, 0.5390625, -0.359375}, {-0.796875, 0.5390625, -0.359375},
    {0.6171875, 0.328125, -0.5859375}, {-0.6171875, 0.328125, -0.5859375},
    {0.484375, 0.0234375, -0.546875}, {-0.484375, 0.0234375, -0.546875},
    {0.8203125, 0.328125, -0.203125}, {-0.8203125, 0.328125, -0.203125},
    {0.40625, -0.171875, 0.1484375}, {-0.40625, -0.171875, 0.1484375},
    {0.4296875, -0.1953125, -0.2109375}, {-0.4296875, -0.1953125, -0.2109375},
    {0.890625, 0.40625, -0.234375}, {-0.890625, 0.40625, -0.234375},
    {0.7734375, -0.140625, -0.125}, {-0.7734375, -0.140625, -0.125},
    {1.0390625, -0.1015625, -0.328125}, {-1.0390625, -0.1015625, -0.328125},
    {1.28125, 0.0546875, -0.4296875}, {-1.28125, 0.0546875, -0.4296875},
    {1.3515625, 0.3203125, -0.421875}, {-1.3515625, 0.3203125, -0.421875},
    {1.234375, 0.5078125, -0.421875}, {-1.234375, 0.5078125, -0.421875},
    {1.0234375, 0.4765625, -0.3125}, {-1.0234375, 0.4765625, -0.3125},
    {1.015625, 0.4140625, -0.2890625}, {-1.015625, 0.4140625, -0.2890625},
    {1.1875, 0.4375, -0.390625}, {-1.1875, 0.4375, -0.390625},
    {1.265625, 0.2890625, -0.40625}, {-1.265625, 0.2890625, -0.40625},
    {1.2109375, 0.078125, -0.40625}, {-1.2109375, 0.078125, -0.40625},
    {1.03125, -0.0390625, -0.3046875}, {-1.03125, -0.0390625, -0.3046875},
    {0.828125, -0.0703125, -0.1328125}, {-0.828125, -0.0703125, -0.1328125},
    {0.921875, 0.359375, -0.21875}, {-0.921875, 0.359375, -0.21875},
    {0.9453125, 0.3046875, -0.2890625}, {-0.9453125, 0.3046875, -0.2890625},
    {0.8828125, -0.0234375, -0.2109375}, {-0.8828125, -0.0234375, -0.2109375},
    {1.0390625, 0.0, -0.3671875}, {-1.0390625, 0.0, -0.3671875},
    {1.1875, 0.09375, -0.4453125}, {-1.1875, 0.09375, -0.4453125},
    {1.234375, 0.25, -0.4453125}, {-1.234375, 0.25, -0.4453125},
    {1.171875, 0.359375, -0.4375}, {-1.171875, 0.359375, -0.4375},
    {1.0234375, 0.34375, -0.359375}, {-1.0234375, 0.34375, -0.359375},
    {0.84375, 0.2890625, -0.2109375}, {-0.84375, 0.2890625, -0.2109375},
    {0.8359375, 0.171875, -0.2734375}, {-0.8359375, 0.171875, -0.2734375},
    {0.7578125, 0.09375, -0.2734375}, {-0.7578125, 0.09375, -0.2734375},
    {0.8203125, 0.0859375, -0.2734375}, {-0.8203125, 0.0859375, -0.2734375},
    {0.84375, 0.015625, -0.2734375}, {-0.84375, 0.015625, -0.2734375},
    {0.8125, -0.015625, -0.2734375}, {-0.8125, -0.015625, -0.2734375},
    {0.7265625, 0.0, -0.0703125}, {-0.7265625, 0.0, -0.0703125},
    {0.71875, -0.0234375, -0.171875}, {-0.71875, -0.0234375, -0.171875},
    {0.71875, 0.0390625, -0.1875}, {-0.71875, 0.0390625, -0.1875},
    {0.796875, 0.203125, -0.2109375}, {-0.796875, 0.203125, -0.2109375},
    {0.890625, 0.2421875, -0.265625}, {-0.890625, 0.2421875, -0.265625},
    {0.890625, 0.234375, -0.3203125}, {-0.890625, 0.234375, -0.3203125},
    {0.8125, -0.015625, -0.3203125}, {-0.8125, -0.015625, -0.3203125},
    {0.8515625, 0.015625, -0.3203125}, {-0.8515625, 0.015625, -0.3203125},
    {0.828125, 0.078125, -0.3203125}, {-0.828125, 0.078125, -0.3203125},
    {0.765625, 0.09375, -0.3203125}, {-0.765625, 0.09375, -0.3203125},
    {0.84375, 0.171875, -0.3203125}, {-0.84375, 0.171875, -0.3203125},
    {1.0390625, 0.328125, -0.4140625}, {-1.0390625, 0.328125, -0.4140625},
    {1.1875, 0.34375, -0.484375}, {-1.1875, 0.34375, -0.484375},
    {1.2578125, 0.2421875, -0.4921875}, {-1.2578125, 0.2421875, -0.4921875},
    {1.2109375, 0.0859375, -0.484375}, {-1.2109375, 0.0859375, -0.484375},
    {1.046875, 0.0, -0.421875}, {-1.046875, 0.0, -0.421875},
    {0.8828125, -0.015625, -0.265625}, {-0.8828125, -0.015625, -0.265625},
    {0.953125, 0.2890625, -0.34375}, {-0.953125, 0.2890625, -0.34375},
    {0.890625, 0.109375, -0.328125}, {-0.890625, 0.109375, -0.328125},
    {0.9375, 0.0625, -0.3359375}, {-0.9375, 0.0625, -0.3359375},
    {1.0, 0.125, -0.3671875}, {-1.0, 0.125, -0.3671875},
    {0.9609375, 0.171875, -0.3515625}, {-0.9609375, 0.171875, -0.3515625},
    {1.015625, 0.234375, -0.375}, {-1.015625, 0.234375, -0.375},
    {1.0546875, 0.1875, -0.3828125}, {-1.0546875, 0.1875, -0.3828125},
    {1.109375, 0.2109375, -0.390625}, {-1.109375, 0.2109375, -0.390625},
    {1.0859375, 0.2734375, -0.390625}, {-1.0859375, 0.2734375, -0.390625},
    {1.0234375, 0.4375, -0.484375}, {-1.0234375, 0.4375, -0.484375},
    {1.25, 0.46875, -0.546875}, {-1.25, 0.46875, -0.546875},
    {1.3671875, 0.296875, -0.5}, {-1.3671875, 0.296875, -0.5},
    {1.3125, 0.0546875, -0.53125}, {-1.3125, 0.0546875, -0.53125},
    {1.0390625, -0.0859375, -0.4921875}, {-1.0390625, -0.0859375, -0.4921875},
    {0.7890625, -0.125, -0.328125}, {-0.7890625, -0.125, -0.328125},
    {0.859375, 0.3828125, -0.3828125}, {-0.859375, 0.3828125, -0.3828125}};

vector<vec4i> suzanne_quads = vector<vec4i>{{46, 0, 2, 44}, {3, 1, 47, 45},
    {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4}, {7, 9, 3, 5}, {0, 10, 8, 2},
    {9, 11, 1, 3}, {10, 12, 14, 8}, {15, 13, 11, 9}, {8, 14, 16, 6},
    {17, 15, 9, 7}, {14, 20, 18, 16}, {19, 21, 15, 17}, {12, 22, 20, 14},
    {21, 23, 13, 15}, {22, 24, 26, 20}, {27, 25, 23, 21}, {20, 26, 28, 18},
    {29, 27, 21, 19}, {26, 32, 30, 28}, {31, 33, 27, 29}, {24, 34, 32, 26},
    {33, 35, 25, 27}, {34, 36, 38, 32}, {39, 37, 35, 33}, {32, 38, 40, 30},
    {41, 39, 33, 31}, {38, 44, 42, 40}, {43, 45, 39, 41}, {36, 46, 44, 38},
    {45, 47, 37, 39}, {46, 36, 50, 48}, {51, 37, 47, 49}, {36, 34, 52, 50},
    {53, 35, 37, 51}, {34, 24, 54, 52}, {55, 25, 35, 53}, {24, 22, 56, 54},
    {57, 23, 25, 55}, {22, 12, 58, 56}, {59, 13, 23, 57}, {12, 10, 62, 58},
    {63, 11, 13, 59}, {10, 0, 64, 62}, {65, 1, 11, 63}, {0, 46, 48, 64},
    {49, 47, 1, 65}, {88, 173, 175, 90}, {175, 174, 89, 90}, {86, 171, 173, 88},
    {174, 172, 87, 89}, {84, 169, 171, 86}, {172, 170, 85, 87},
    {82, 167, 169, 84}, {170, 168, 83, 85}, {80, 165, 167, 82},
    {168, 166, 81, 83}, {78, 91, 145, 163}, {146, 92, 79, 164},
    {91, 93, 147, 145}, {148, 94, 92, 146}, {93, 95, 149, 147},
    {150, 96, 94, 148}, {95, 97, 151, 149}, {152, 98, 96, 150},
    {97, 99, 153, 151}, {154, 100, 98, 152}, {99, 101, 155, 153},
    {156, 102, 100, 154}, {101, 103, 157, 155}, {158, 104, 102, 156},
    {103, 105, 159, 157}, {160, 106, 104, 158}, {105, 107, 161, 159},
    {162, 108, 106, 160}, {107, 66, 67, 161}, {67, 66, 108, 162},
    {109, 127, 159, 161}, {160, 128, 110, 162}, {127, 178, 157, 159},
    {158, 179, 128, 160}, {125, 155, 157, 178}, {158, 156, 126, 179},
    {123, 153, 155, 125}, {156, 154, 124, 126}, {121, 151, 153, 123},
    {154, 152, 122, 124}, {119, 149, 151, 121}, {152, 150, 120, 122},
    {117, 147, 149, 119}, {150, 148, 118, 120}, {115, 145, 147, 117},
    {148, 146, 116, 118}, {113, 163, 145, 115}, {146, 164, 114, 116},
    {113, 180, 176, 163}, {176, 181, 114, 164}, {109, 161, 67, 111},
    {67, 162, 110, 112}, {111, 67, 177, 182}, {177, 67, 112, 183},
    {176, 180, 182, 177}, {183, 181, 176, 177}, {134, 136, 175, 173},
    {175, 136, 135, 174}, {132, 134, 173, 171}, {174, 135, 133, 172},
    {130, 132, 171, 169}, {172, 133, 131, 170}, {165, 186, 184, 167},
    {185, 187, 166, 168}, {130, 169, 167, 184}, {168, 170, 131, 185},
    {143, 189, 188, 186}, {188, 189, 144, 187}, {184, 186, 188, 68},
    {188, 187, 185, 68}, {129, 130, 184, 68}, {185, 131, 129, 68},
    {141, 192, 190, 143}, {191, 193, 142, 144}, {139, 194, 192, 141},
    {193, 195, 140, 142}, {138, 196, 194, 139}, {195, 197, 138, 140},
    {137, 70, 196, 138}, {197, 70, 137, 138}, {189, 143, 190, 69},
    {191, 144, 189, 69}, {69, 190, 205, 207}, {206, 191, 69, 207},
    {70, 198, 199, 196}, {200, 198, 70, 197}, {196, 199, 201, 194},
    {202, 200, 197, 195}, {194, 201, 203, 192}, {204, 202, 195, 193},
    {192, 203, 205, 190}, {206, 204, 193, 191}, {198, 203, 201, 199},
    {202, 204, 198, 200}, {198, 207, 205, 203}, {206, 207, 198, 204},
    {138, 139, 163, 176}, {164, 140, 138, 176}, {139, 141, 210, 163},
    {211, 142, 140, 164}, {141, 143, 212, 210}, {213, 144, 142, 211},
    {143, 186, 165, 212}, {166, 187, 144, 213}, {80, 208, 212, 165},
    {213, 209, 81, 166}, {208, 214, 210, 212}, {211, 215, 209, 213},
    {78, 163, 210, 214}, {211, 164, 79, 215}, {130, 129, 71, 221},
    {71, 129, 131, 222}, {132, 130, 221, 219}, {222, 131, 133, 220},
    {134, 132, 219, 217}, {220, 133, 135, 218}, {136, 134, 217, 216},
    {218, 135, 136, 216}, {216, 217, 228, 230}, {229, 218, 216, 230},
    {217, 219, 226, 228}, {227, 220, 218, 229}, {219, 221, 224, 226},
    {225, 222, 220, 227}, {221, 71, 223, 224}, {223, 71, 222, 225},
    {223, 230, 228, 224}, {229, 230, 223, 225}, {182, 180, 233, 231},
    {234, 181, 183, 232}, {111, 182, 231, 253}, {232, 183, 112, 254},
    {109, 111, 253, 255}, {254, 112, 110, 256}, {180, 113, 251, 233},
    {252, 114, 181, 234}, {113, 115, 249, 251}, {250, 116, 114, 252},
    {115, 117, 247, 249}, {248, 118, 116, 250}, {117, 119, 245, 247},
    {246, 120, 118, 248}, {119, 121, 243, 245}, {244, 122, 120, 246},
    {121, 123, 241, 243}, {242, 124, 122, 244}, {123, 125, 239, 241},
    {240, 126, 124, 242}, {125, 178, 235, 239}, {236, 179, 126, 240},
    {178, 127, 237, 235}, {238, 128, 179, 236}, {127, 109, 255, 237},
    {256, 110, 128, 238}, {237, 255, 257, 275}, {258, 256, 238, 276},
    {235, 237, 275, 277}, {276, 238, 236, 278}, {239, 235, 277, 273},
    {278, 236, 240, 274}, {241, 239, 273, 271}, {274, 240, 242, 272},
    {243, 241, 271, 269}, {272, 242, 244, 270}, {245, 243, 269, 267},
    {270, 244, 246, 268}, {247, 245, 267, 265}, {268, 246, 248, 266},
    {249, 247, 265, 263}, {266, 248, 250, 264}, {251, 249, 263, 261},
    {264, 250, 252, 262}, {233, 251, 261, 279}, {262, 252, 234, 280},
    {255, 253, 259, 257}, {260, 254, 256, 258}, {253, 231, 281, 259},
    {282, 232, 254, 260}, {231, 233, 279, 281}, {280, 234, 232, 282},
    {66, 107, 283, 72}, {284, 108, 66, 72}, {107, 105, 285, 283},
    {286, 106, 108, 284}, {105, 103, 287, 285}, {288, 104, 106, 286},
    {103, 101, 289, 287}, {290, 102, 104, 288}, {101, 99, 291, 289},
    {292, 100, 102, 290}, {99, 97, 293, 291}, {294, 98, 100, 292},
    {97, 95, 295, 293}, {296, 96, 98, 294}, {95, 93, 297, 295},
    {298, 94, 96, 296}, {93, 91, 299, 297}, {300, 92, 94, 298},
    {307, 308, 327, 337}, {328, 308, 307, 338}, {306, 307, 337, 335},
    {338, 307, 306, 336}, {305, 306, 335, 339}, {336, 306, 305, 340},
    {88, 90, 305, 339}, {305, 90, 89, 340}, {86, 88, 339, 333},
    {340, 89, 87, 334}, {84, 86, 333, 329}, {334, 87, 85, 330},
    {82, 84, 329, 331}, {330, 85, 83, 332}, {329, 335, 337, 331},
    {338, 336, 330, 332}, {329, 333, 339, 335}, {340, 334, 330, 336},
    {325, 331, 337, 327}, {338, 332, 326, 328}, {80, 82, 331, 325},
    {332, 83, 81, 326}, {208, 341, 343, 214}, {344, 342, 209, 215},
    {80, 325, 341, 208}, {342, 326, 81, 209}, {78, 214, 343, 345},
    {344, 215, 79, 346}, {78, 345, 299, 91}, {300, 346, 79, 92},
    {76, 323, 351, 303}, {352, 324, 76, 303}, {303, 351, 349, 77},
    {350, 352, 303, 77}, {77, 349, 347, 304}, {348, 350, 77, 304},
    {304, 347, 327, 308}, {328, 348, 304, 308}, {325, 327, 347, 341},
    {348, 328, 326, 342}, {295, 297, 317, 309}, {318, 298, 296, 310},
    {75, 315, 323, 76}, {324, 316, 75, 76}, {301, 357, 355, 302},
    {356, 358, 301, 302}, {302, 355, 353, 74}, {354, 356, 302, 74},
    {74, 353, 315, 75}, {316, 354, 74, 75}, {291, 293, 361, 363},
    {362, 294, 292, 364}, {363, 361, 367, 365}, {368, 362, 364, 366},
    {365, 367, 369, 371}, {370, 368, 366, 372}, {371, 369, 375, 373},
    {376, 370, 372, 374}, {313, 377, 373, 375}, {374, 378, 314, 376},
    {315, 353, 373, 377}, {374, 354, 316, 378}, {353, 355, 371, 373},
    {372, 356, 354, 374}, {355, 357, 365, 371}, {366, 358, 356, 372},
    {357, 359, 363, 365}, {364, 360, 358, 366}, {289, 291, 363, 359},
    {364, 292, 290, 360}, {73, 359, 357, 301}, {358, 360, 73, 301},
    {283, 285, 287, 289}, {288, 286, 284, 290}, {283, 289, 359, 73},
    {360, 290, 284, 73}, {293, 295, 309, 361}, {310, 296, 294, 362},
    {309, 311, 367, 361}, {368, 312, 310, 362}, {311, 381, 369, 367},
    {370, 382, 312, 368}, {313, 375, 369, 381}, {370, 376, 314, 382},
    {347, 349, 385, 383}, {386, 350, 348, 384}, {317, 383, 385, 319},
    {386, 384, 318, 320}, {297, 299, 383, 317}, {384, 300, 298, 318},
    {299, 343, 341, 383}, {342, 344, 300, 384}, {313, 321, 379, 377},
    {380, 322, 314, 378}, {315, 377, 379, 323}, {380, 378, 316, 324},
    {319, 385, 379, 321}, {380, 386, 320, 322}, {349, 351, 379, 385},
    {380, 352, 350, 386}, {399, 387, 413, 401}, {414, 388, 400, 402},
    {399, 401, 403, 397}, {404, 402, 400, 398}, {397, 403, 405, 395},
    {406, 404, 398, 396}, {395, 405, 407, 393}, {408, 406, 396, 394},
    {393, 407, 409, 391}, {410, 408, 394, 392}, {391, 409, 411, 389},
    {412, 410, 392, 390}, {409, 419, 417, 411}, {418, 420, 410, 412},
    {407, 421, 419, 409}, {420, 422, 408, 410}, {405, 423, 421, 407},
    {422, 424, 406, 408}, {403, 425, 423, 405}, {424, 426, 404, 406},
    {401, 427, 425, 403}, {426, 428, 402, 404}, {401, 413, 415, 427},
    {416, 414, 402, 428}, {317, 319, 443, 441}, {444, 320, 318, 442},
    {319, 389, 411, 443}, {412, 390, 320, 444}, {309, 317, 441, 311},
    {442, 318, 310, 312}, {381, 429, 413, 387}, {414, 430, 382, 388},
    {411, 417, 439, 443}, {440, 418, 412, 444}, {437, 445, 443, 439},
    {444, 446, 438, 440}, {433, 445, 437, 435}, {438, 446, 434, 436},
    {431, 447, 445, 433}, {446, 448, 432, 434}, {429, 447, 431, 449},
    {432, 448, 430, 450}, {413, 429, 449, 415}, {450, 430, 414, 416},
    {311, 447, 429, 381}, {430, 448, 312, 382}, {311, 441, 445, 447},
    {446, 442, 312, 448}, {415, 449, 451, 475}, {452, 450, 416, 476},
    {449, 431, 461, 451}, {462, 432, 450, 452}, {431, 433, 459, 461},
    {460, 434, 432, 462}, {433, 435, 457, 459}, {458, 436, 434, 460},
    {435, 437, 455, 457}, {456, 438, 436, 458}, {437, 439, 453, 455},
    {454, 440, 438, 456}, {439, 417, 473, 453}, {474, 418, 440, 454},
    {427, 415, 475, 463}, {476, 416, 428, 464}, {425, 427, 463, 465},
    {464, 428, 426, 466}, {423, 425, 465, 467}, {466, 426, 424, 468},
    {421, 423, 467, 469}, {468, 424, 422, 470}, {419, 421, 469, 471},
    {470, 422, 420, 472}, {417, 419, 471, 473}, {472, 420, 418, 474},
    {457, 455, 479, 477}, {480, 456, 458, 478}, {477, 479, 481, 483},
    {482, 480, 478, 484}, {483, 481, 487, 485}, {488, 482, 484, 486},
    {485, 487, 489, 491}, {490, 488, 486, 492}, {463, 475, 485, 491},
    {486, 476, 464, 492}, {451, 483, 485, 475}, {486, 484, 452, 476},
    {451, 461, 477, 483}, {478, 462, 452, 484}, {457, 477, 461, 459},
    {462, 478, 458, 460}, {453, 473, 479, 455}, {480, 474, 454, 456},
    {471, 481, 479, 473}, {480, 482, 472, 474}, {469, 487, 481, 471},
    {482, 488, 470, 472}, {467, 489, 487, 469}, {488, 490, 468, 470},
    {465, 491, 489, 467}, {490, 492, 466, 468}, {391, 389, 503, 501},
    {504, 390, 392, 502}, {393, 391, 501, 499}, {502, 392, 394, 500},
    {395, 393, 499, 497}, {500, 394, 396, 498}, {397, 395, 497, 495},
    {498, 396, 398, 496}, {399, 397, 495, 493}, {496, 398, 400, 494},
    {387, 399, 493, 505}, {494, 400, 388, 506}, {493, 501, 503, 505},
    {504, 502, 494, 506}, {493, 495, 499, 501}, {500, 496, 494, 502},
    {313, 381, 387, 505}, {388, 382, 314, 506}, {313, 505, 503, 321},
    {504, 506, 314, 322}, {319, 321, 503, 389}, {504, 322, 320, 390},
    // ttriangles
    {60, 64, 48, 48}, {49, 65, 61, 61}, {62, 64, 60, 60}, {61, 65, 63, 63},
    {60, 58, 62, 62}, {63, 59, 61, 61}, {60, 56, 58, 58}, {59, 57, 61, 61},
    {60, 54, 56, 56}, {57, 55, 61, 61}, {60, 52, 54, 54}, {55, 53, 61, 61},
    {60, 50, 52, 52}, {53, 51, 61, 61}, {60, 48, 50, 50}, {51, 49, 61, 61},
    {224, 228, 226, 226}, {227, 229, 225, 255}, {72, 283, 73, 73},
    {73, 284, 72, 72}, {341, 347, 383, 383}, {384, 348, 342, 342},
    {299, 345, 343, 343}, {344, 346, 300, 300}, {323, 379, 351, 351},
    {352, 380, 324, 324}, {441, 443, 445, 445}, {446, 444, 442, 442},
    {463, 491, 465, 465}, {466, 492, 464, 464}, {495, 497, 499, 499},
    {500, 498, 496, 496}};

}  // namespace yocto
