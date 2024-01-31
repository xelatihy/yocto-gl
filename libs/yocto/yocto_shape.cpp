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
#include "yocto_modeling.h"
#include "yocto_modelio.h"
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
vec3f eval_position(const shape_data& shape, int element, vec2f uv) {
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

vec3f eval_normal(const shape_data& shape, int element, vec2f uv) {
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

vec3f eval_tangent(const shape_data& shape, int element, vec2f uv) {
  return eval_normal(shape, element, uv);
}

vec2f eval_texcoord(const shape_data& shape, int element, vec2f uv) {
  if (shape.texcoords.empty()) return uv;
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
    return uv;
  }
}

vec4f eval_color(const shape_data& shape, int element, vec2f uv) {
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
    return {1, 1, 1, 1};
  }
}

float eval_radius(const shape_data& shape, int element, vec2f uv) {
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

shape_point sample_shape(
    const shape_data& shape, const vector<float>& cdf, float rn, vec2f ruv) {
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
    const shape_data& shape, const frame3f& frame, bool non_rigid) {
  return transform_shape(shape, frame, 1, non_rigid);
}
shape_data transform_shape(const shape_data& shape, const frame3f& frame,
    float radius_scale, bool non_rigid) {
  auto transformed = shape;
  for (auto& position : transformed.positions)
    position = transform_point(frame, position);
  for (auto& normal : transformed.normals)
    normal = transform_normal(frame, normal, non_rigid);
  for (auto& radius : transformed.radius) radius *= radius_scale;
  return transformed;
}
shape_data scale_shape(const shape_data& shape, float scale, float uvscale) {
  if (scale == 1 && uvscale == 1) return shape;
  auto transformed = shape;
  for (auto& position : transformed.positions) position *= scale;
  for (auto& texcoord : transformed.texcoords) texcoord *= uvscale;
  return transformed;
}
shape_data scale_shape(shape_data&& shape, float scale, float uvscale) {
  if (scale == 1 && uvscale == 1) return std::move(shape);
  auto transformed = std::move(shape);
  for (auto& position : transformed.positions) position *= scale;
  for (auto& texcoord : transformed.texcoords) texcoord *= uvscale;
  return transformed;
}

// Manipulate vertex data
shape_data remove_normals(const shape_data& shape) {
  auto transformed    = shape;
  transformed.normals = {};
  return transformed;
}
shape_data add_normals(const shape_data& shape) {
  auto transformed    = shape;
  transformed.normals = compute_normals(shape);
  return transformed;
}

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR FVSHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, vec2f uv) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const fvshape_data& shape, int element, vec2f uv) {
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

vec2f eval_texcoord(const fvshape_data& shape, int element, vec2f uv) {
  if (shape.texcoords.empty()) return uv;
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadstexcoord[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return uv;
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

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
void split_facevarying(vector<vec4i>& split_quads,
    vector<vec3f>& split_positions, vector<vec3f>& split_normals,
    vector<vec2f>& split_texcoords, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  // make vertices
  auto vertices = vector<vec3i>{};
  vertices.reserve(quadspos.size() * 4);
  for (auto fid : range(quadspos.size())) {
    for (auto c : range(4)) {
      auto vertex = vec3i{
          (&quadspos[fid].x)[c],
          (!quadsnorm.empty()) ? (&quadsnorm[fid].x)[c] : -1,
          (!quadstexcoord.empty()) ? (&quadstexcoord[fid].x)[c] : -1,
      };
      vertices.push_back(vertex);
    }
  }

  // sort vertices and remove duplicates
  auto compare_vertices = [](vec3i a, vec3i b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y) ||
           (a.x == b.x && a.y == b.y && a.z < b.z);
  };
  std::sort(vertices.begin(), vertices.end(), compare_vertices);
  vertices.erase(std::unique(vertices.begin(), vertices.end()), vertices.end());

  // fill vert data
  if (!positions.empty()) {
    split_positions.resize(vertices.size());
    for (auto&& [index, vertex] : enumerate(vertices)) {
      split_positions[index] = positions[vertex.x];
    }
  } else {
    split_positions.clear();
  }
  if (!normals.empty()) {
    split_normals.resize(vertices.size());
    for (auto&& [index, vertex] : enumerate(vertices)) {
      split_normals[index] = normals[vertex.y];
    }
  } else {
    split_normals.clear();
  }
  if (!texcoords.empty()) {
    split_texcoords.resize(vertices.size());
    for (auto&& [index, vertex] : enumerate(vertices)) {
      split_texcoords[index] = texcoords[vertex.z];
    }
  } else {
    split_texcoords.clear();
  }

  // fill quad data
  split_quads.resize(quadspos.size());
  for (auto fid : range(quadspos.size())) {
    for (auto c : range(4)) {
      auto vertex = vec3i{
          (&quadspos[fid].x)[c],
          (!quadsnorm.empty()) ? (&quadsnorm[fid].x)[c] : -1,
          (!quadstexcoord.empty()) ? (&quadstexcoord[fid].x)[c] : -1,
      };
      split_quads[fid][c] =
          (int)(std::lower_bound(vertices.begin(), vertices.end(), vertex,
                    compare_vertices) -
                vertices.begin());
    }
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
    throw std::invalid_argument{"cannot convert shape"};
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
    const fvshape_data& shape, const frame3f& frame, bool non_rigid) {
  auto transformed = shape;
  for (auto& position : transformed.positions)
    position = transform_point(frame, position);
  for (auto& normal : transformed.normals)
    normal = transform_normal(frame, normal, non_rigid);
  return transformed;
}
fvshape_data scale_fvshape(
    const fvshape_data& shape, float scale, float uvscale) {
  if (scale == 1 && uvscale == 1) return shape;
  auto transformed = shape;
  for (auto& position : transformed.positions) position *= scale;
  for (auto& texcoord : transformed.texcoords) texcoord *= uvscale;
  return transformed;
}
fvshape_data scale_fvshape(fvshape_data&& shape, float scale, float uvscale) {
  if (scale == 1 && uvscale == 1) return std::move(shape);
  auto transformed = std::move(shape);
  for (auto& position : transformed.positions) position *= scale;
  for (auto& texcoord : transformed.texcoords) texcoord *= uvscale;
  return transformed;
}

// Vertex properties
fvshape_data remove_normals(const fvshape_data& shape) {
  auto transformed      = shape;
  transformed.quadsnorm = {};
  transformed.normals   = {};
  return transformed;
}
fvshape_data add_normals(const fvshape_data& shape) {
  auto transformed      = shape;
  transformed.quadsnorm = transformed.quadspos;
  transformed.normals   = compute_normals(shape);
  return transformed;
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
// IMPLEMENTATION FOR SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a tesselated rectangle. Useful in other subdivisions.
static shape_data make_quads(vec2i steps, vec2f scale, vec2f uvscale) {
  auto shape = shape_data{};

  shape.positions.resize((steps.x + 1) * (steps.y + 1));
  shape.normals.resize((steps.x + 1) * (steps.y + 1));
  shape.texcoords.resize((steps.x + 1) * (steps.y + 1));
  for (auto j : range(steps.y + 1)) {
    for (auto i : range(steps.x + 1)) {
      auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
      shape.positions[j * (steps.x + 1) + i] = {
          (2 * uv.x - 1) * scale.x, (2 * uv.y - 1) * scale.y, 0};
      shape.normals[j * (steps.x + 1) + i]   = {0, 0, 1};
      shape.texcoords[j * (steps.x + 1) + i] = vec2f{uv.x, 1 - uv.y} * uvscale;
    }
  }

  shape.quads.resize(steps.x * steps.y);
  for (auto j : range(steps.y)) {
    for (auto i : range(steps.x)) {
      shape.quads[j * steps.x + i] = {j * (steps.x + 1) + i,
          j * (steps.x + 1) + i + 1, (j + 1) * (steps.x + 1) + i + 1,
          (j + 1) * (steps.x + 1) + i};
    }
  }

  return shape;
}

// Merge shape elements
void merge_shape_inplace(shape_data& shape, const shape_data& merge) {
  auto offset = (int)shape.positions.size();
  for (auto& p : merge.points) shape.points.push_back(p + offset);
  for (auto& l : merge.lines)
    shape.lines.push_back({l.x + offset, l.y + offset});
  for (auto& t : merge.triangles)
    shape.triangles.push_back({t.x + offset, t.y + offset, t.z + offset});
  for (auto& q : merge.quads)
    shape.quads.push_back(
        {q.x + offset, q.y + offset, q.z + offset, q.w + offset});
  shape.positions.insert(
      shape.positions.end(), merge.positions.begin(), merge.positions.end());
  shape.tangents.insert(
      shape.tangents.end(), merge.tangents.begin(), merge.tangents.end());
  shape.texcoords.insert(
      shape.texcoords.end(), merge.texcoords.begin(), merge.texcoords.end());
  shape.colors.insert(
      shape.colors.end(), merge.colors.begin(), merge.colors.end());
  shape.radius.insert(
      shape.radius.end(), merge.radius.begin(), merge.radius.end());
}

// Make a plane.
shape_data make_rect(vec2i steps, vec2f scale, vec2f uvscale) {
  return make_quads(steps, scale, uvscale);
}
shape_data make_bulged_rect(
    vec2i steps, vec2f scale, vec2f uvscale, float height) {
  auto shape = make_rect(steps, scale, uvscale);
  if (height != 0) {
    height      = min(height, min(scale));
    auto radius = (1 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i : range(shape.positions.size())) {
      auto pn            = normalize(shape.positions[i] - center);
      shape.positions[i] = center + pn * radius;
      shape.normals[i]   = pn;
    }
  }
  return shape;
}

// Make a plane in the xz plane.
shape_data make_recty(vec2i steps, vec2f scale, vec2f uvscale) {
  auto shape = make_rect(steps, scale, uvscale);
  for (auto& position : shape.positions)
    position = {position.x, position.z, -position.y};
  for (auto& normal : shape.normals) normal = {normal.x, normal.z, normal.y};
  return shape;
}
shape_data make_bulged_recty(
    vec2i steps, vec2f scale, vec2f uvscale, float height) {
  auto shape = make_bulged_rect(steps, scale, uvscale, height);
  for (auto& position : shape.positions)
    position = {position.x, position.z, -position.y};
  for (auto& normal : shape.normals) normal = {normal.x, normal.z, normal.y};
  return shape;
}

// Make a box.
shape_data make_box(vec3i steps, vec3f scale, vec3f uvscale) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // + z
  qshape = make_rect(
      {steps.x, steps.y}, {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qshape.positions) p = {p.x, p.y, scale.z};
  for (auto& n : qshape.normals) n = {0, 0, 1};
  merge_shape_inplace(shape, qshape);
  // - z
  qshape = make_rect(
      {steps.x, steps.y}, {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qshape.positions) p = {-p.x, p.y, -scale.z};
  for (auto& n : qshape.normals) n = {0, 0, -1};
  merge_shape_inplace(shape, qshape);
  // + x
  qshape = make_rect(
      {steps.z, steps.y}, {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qshape.positions) p = {scale.x, p.y, -p.x};
  for (auto& n : qshape.normals) n = {1, 0, 0};
  merge_shape_inplace(shape, qshape);
  // - x
  qshape = make_rect(
      {steps.z, steps.y}, {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qshape.positions) p = {-scale.x, p.y, p.x};
  for (auto& n : qshape.normals) n = {-1, 0, 0};
  merge_shape_inplace(shape, qshape);
  // + y
  qshape = make_rect(
      {steps.x, steps.z}, {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i : range(qshape.positions.size())) {
    qshape.positions[i] = {
        qshape.positions[i].x, scale.y, -qshape.positions[i].y};
    qshape.normals[i] = {0, 1, 0};
  }
  merge_shape_inplace(shape, qshape);
  // - y
  qshape = make_rect(
      {steps.x, steps.z}, {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i : range(qshape.positions.size())) {
    qshape.positions[i] = {
        qshape.positions[i].x, -scale.y, qshape.positions[i].y};
    qshape.normals[i] = {0, -1, 0};
  }
  merge_shape_inplace(shape, qshape);
  return shape;
}
shape_data make_rounded_box(
    vec3i steps, vec3f scale, vec3f uvscale, float radius) {
  auto shape = make_box(steps, scale, uvscale);
  if (radius != 0) {
    radius = min(radius, min(scale));
    auto c = scale - radius;
    for (auto i : range(shape.positions.size())) {
      auto pc = vec3f{abs(shape.positions[i].x), abs(shape.positions[i].y),
          abs(shape.positions[i].z)};
      auto ps = vec3f{shape.positions[i].x < 0 ? -1.0f : 1.0f,
          shape.positions[i].y < 0 ? -1.0f : 1.0f,
          shape.positions[i].z < 0 ? -1.0f : 1.0f};
      if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
        auto pn            = normalize(pc - c);
        shape.positions[i] = c + radius * pn;
        shape.normals[i]   = pn;
      } else if (pc.x >= c.x && pc.y >= c.y) {
        auto pn            = normalize((pc - c) * vec3f{1, 1, 0});
        shape.positions[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
        shape.normals[i]   = pn;
      } else if (pc.x >= c.x && pc.z >= c.z) {
        auto pn            = normalize((pc - c) * vec3f{1, 0, 1});
        shape.positions[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
        shape.normals[i]   = pn;
      } else if (pc.y >= c.y && pc.z >= c.z) {
        auto pn            = normalize((pc - c) * vec3f{0, 1, 1});
        shape.positions[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
        shape.normals[i]   = pn;
      } else {
        continue;
      }
      shape.positions[i] *= ps;
      shape.normals[i] *= ps;
    }
  }
  return shape;
}

// Make a quad stack
shape_data make_rect_stack(vec3i steps, vec3f scale, vec2f uvscale) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  for (auto i : range(steps.z + 1)) {
    qshape = make_rect({steps.x, steps.y}, {scale.x, scale.y}, uvscale);
    for (auto& p : qshape.positions)
      p.z = (-1 + 2 * (float)i / steps.z) * scale.z;
    merge_shape_inplace(shape, qshape);
  }
  return shape;
}

// Make a floor.
shape_data make_floor(vec2i steps, vec2f scale, vec2f uvscale) {
  auto shape = make_rect(steps, scale, uvscale);
  for (auto& position : shape.positions)
    position = {position.x, position.z, -position.y};
  for (auto& normal : shape.normals) normal = {normal.x, normal.z, normal.y};
  return shape;
}
shape_data make_bent_floor(
    vec2i steps, vec2f scale, vec2f uvscale, float radius) {
  auto shape = make_floor(steps, scale, uvscale);
  if (radius != 0) {
    radius     = min(radius, scale.y);
    auto start = (scale.y - radius) / 2;
    auto end   = start + radius;
    for (auto i : range(shape.positions.size())) {
      if (shape.positions[i].z < -end) {
        shape.positions[i] = {
            shape.positions[i].x, -shape.positions[i].z - end + radius, -end};
        shape.normals[i] = {0, 0, 1};
      } else if (shape.positions[i].z < -start &&
                 shape.positions[i].z >= -end) {
        auto phi = (pif / 2) * (-shape.positions[i].z - start) / radius;
        shape.positions[i] = {shape.positions[i].x, -cos(phi) * radius + radius,
            -sin(phi) * radius - start};
        shape.normals[i]   = {0, cos(phi), sin(phi)};
      } else {
      }
    }
  }
  return shape;
}

// Make a sphere.
shape_data make_sphere(int steps, float scale, float uvscale) {
  auto shape = make_box({steps, steps, steps}, {scale, scale, scale},
      {uvscale, uvscale, uvscale});
  for (auto& p : shape.positions) p = normalize(p) * scale;
  shape.normals = shape.positions;
  for (auto& n : shape.normals) n = normalize(n);
  return shape;
}

// Make a sphere.
shape_data make_uvsphere(vec2i steps, float scale, vec2f uvscale) {
  auto shape = make_rect(steps, {1, 1});
  for (auto i : range(shape.positions.size())) {
    auto uv = shape.texcoords[i];
    auto a  = vec2f{2 * pif * uv.x, pif * uv.y};
    shape.positions[i] =
        vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)} * scale;
    shape.normals[i]   = normalize(shape.positions[i]);
    shape.texcoords[i] = uv * uvscale;
  }
  return shape;
}

// Make a sphere.
shape_data make_uvspherey(vec2i steps, float scale, vec2f uvscale) {
  auto shape = make_uvsphere(steps, scale, uvscale);
  for (auto& position : shape.positions)
    position = {position.x, position.z, position.y};
  for (auto& normal : shape.normals) normal = {normal.x, normal.z, normal.y};
  for (auto& texcoord : shape.texcoords) texcoord = {texcoord.x, texcoord.y};
  for (auto& quad : shape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  return shape;
}

// Make a sphere with slipped caps.
shape_data make_capped_uvsphere(
    vec2i steps, float scale, vec2f uvscale, float cap) {
  auto shape = make_uvsphere(steps, scale, uvscale);
  if (cap != 0) {
    cap        = min(cap, scale / 2);
    auto zflip = (scale - cap);
    for (auto i : range(shape.positions.size())) {
      if (shape.positions[i].z > zflip) {
        shape.positions[i].z = 2 * zflip - shape.positions[i].z;
        shape.normals[i].x   = -shape.normals[i].x;
        shape.normals[i].y   = -shape.normals[i].y;
      } else if (shape.positions[i].z < -zflip) {
        shape.positions[i].z = 2 * (-zflip) - shape.positions[i].z;
        shape.normals[i].x   = -shape.normals[i].x;
        shape.normals[i].y   = -shape.normals[i].y;
      }
    }
  }
  return shape;
}

// Make a sphere with slipped caps.
shape_data make_capped_uvspherey(
    vec2i steps, float scale, vec2f uvscale, float cap) {
  auto shape = make_capped_uvsphere(steps, scale, uvscale, cap);
  for (auto& position : shape.positions)
    position = {position.x, position.z, position.y};
  for (auto& normal : shape.normals) normal = {normal.x, normal.z, normal.y};
  for (auto& texcoord : shape.texcoords)
    texcoord = {texcoord.x, 1 - texcoord.y};
  for (auto& quad : shape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  return shape;
}

// Make a disk
shape_data make_disk(int steps, float scale, float uvscale) {
  auto shape = make_rect({steps, steps}, {1, 1}, {uvscale, uvscale});
  for (auto& position : shape.positions) {
    // Analytical Methods for Squaring the Disc, by C. Fong
    // https://arxiv.org/abs/1509.06344
    auto xy = vec2f{position.x, position.y};
    auto uv = vec2f{
        xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
    position = vec3f{uv.x, uv.y, 0} * scale;
  }
  return shape;
}

// Make a bulged disk
shape_data make_bulged_disk(
    int steps, float scale, float uvscale, float height) {
  auto shape = make_disk(steps, scale, uvscale);
  if (height != 0) {
    height      = min(height, scale);
    auto radius = (1 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i : range(shape.positions.size())) {
      auto pn            = normalize(shape.positions[i] - center);
      shape.positions[i] = center + pn * radius;
      shape.normals[i]   = pn;
    }
  }
  return shape;
}

// Make a uv disk
shape_data make_uvdisk(vec2i steps, float scale, vec2f uvscale) {
  auto shape = make_rect(steps, {1, 1}, {1, 1});
  for (auto i : range(shape.positions.size())) {
    auto uv            = shape.texcoords[i];
    auto phi           = 2 * pif * uv.x;
    shape.positions[i] = vec3f{cos(phi) * uv.y, sin(phi) * uv.y, 0} * scale;
    shape.normals[i]   = {0, 0, 1};
    shape.texcoords[i] = uv * uvscale;
  }
  return shape;
}

// Make a uv cylinder
shape_data make_uvcylinder(vec3i steps, vec2f scale, vec3f uvscale) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // side
  qshape = make_rect({steps.x, steps.y}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    qshape.positions[i] = {
        cos(phi) * scale.x, sin(phi) * scale.x, -(2 * uv.y - 1) * scale.y};
    qshape.normals[i]   = {cos(phi), sin(phi), 0};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.y};
  }
  for (auto& quad : qshape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  merge_shape_inplace(shape, qshape);
  // top
  qshape = make_rect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    qshape.positions[i] = {
        cos(phi) * uv.y * scale.x, sin(phi) * uv.y * scale.x, 0};
    qshape.normals[i]     = {0, 0, 1};
    qshape.texcoords[i]   = uv * vec2f{uvscale.x, uvscale.z};
    qshape.positions[i].z = scale.y;
  }
  merge_shape_inplace(shape, qshape);
  // bottom
  qshape = make_rect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    qshape.positions[i] = {
        cos(phi) * uv.y * scale.x, sin(phi) * uv.y * scale.x, 0};
    qshape.normals[i]     = {0, 0, 1};
    qshape.texcoords[i]   = uv * vec2f{uvscale.x, uvscale.z};
    qshape.positions[i].z = -scale.y;
    qshape.normals[i]     = -qshape.normals[i];
  }
  for (auto& qquad : qshape.quads) swap(qquad.x, qquad.z);
  merge_shape_inplace(shape, qshape);
  return shape;
}

// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(
    vec3i steps, vec2f scale, vec3f uvscale, float radius) {
  auto shape = make_uvcylinder(steps, scale, uvscale);
  if (radius != 0) {
    radius = min(radius, min(scale));
    auto c = scale - radius;
    for (auto i : range(shape.positions.size())) {
      auto phi = atan2(shape.positions[i].y, shape.positions[i].x);
      auto r   = length(vec2f{shape.positions[i].x, shape.positions[i].y});
      auto z   = shape.positions[i].z;
      auto pc  = vec2f{r, abs(z)};
      auto ps  = (z < 0) ? -1.0f : 1.0f;
      if (pc.x >= c.x && pc.y >= c.y) {
        auto pn            = normalize(pc - c);
        shape.positions[i] = {cos(phi) * (c.x + radius * pn.x),
            sin(phi) * (c.x + radius * pn.x), ps * (c.y + radius * pn.y)};
        shape.normals[i]   = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
      } else {
        continue;
      }
    }
  }
  return shape;
}

// Make a uv capsule
shape_data make_uvcapsule(vec3i steps, vec2f scale, vec3f uvscale) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // side
  qshape = make_rect({steps.x, steps.y}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    qshape.positions[i] = {
        cos(phi) * scale.x, sin(phi) * scale.x, (2 * uv.y - 1) * scale.y};
    qshape.normals[i]   = {cos(phi), sin(phi), 0};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.y};
  }
  for (auto& quad : qshape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  merge_shape_inplace(shape, qshape);
  // top
  qshape = make_rect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv  = qshape.texcoords[i];
    auto phi = uv.x * 2 * pif, theta = (1 - uv.y) * pif / 2;
    qshape.positions[i] = {cos(phi) * sin(theta) * scale.x,
        sin(phi) * sin(theta) * scale.x, cos(theta) * scale.x + scale.y};
    qshape.normals[i]   = {
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.z};
  }
  merge_shape_inplace(shape, qshape);
  // bottom
  qshape = make_rect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv  = qshape.texcoords[i];
    auto phi = uv.x * 2 * pif, theta = (1 - uv.y) * pif / 2;
    qshape.positions[i] = {cos(phi) * sin(theta) * scale.x,
        sin(phi) * sin(theta) * scale.x, -cos(theta) * scale.x - scale.y};
    qshape.normals[i]   = {
        cos(phi) * sin(theta), sin(phi) * sin(theta), -cos(theta)};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.z};
    qshape.normals[i]   = -qshape.normals[i];
  }
  for (auto& qquad : qshape.quads) swap(qquad.x, qquad.z);
  merge_shape_inplace(shape, qshape);
  return shape;
}

// Make a uv cone
shape_data make_uvcone(vec3i steps, vec2f scale, vec3f uvscale) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // side
  qshape = make_rect({steps.x, steps.y}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    auto normal2d       = normalize(scale * vec2f{1, 2});
    qshape.positions[i] = {cos(phi) * (1 - uv.y) * scale.x,
        sin(phi) * (1 - uv.y) * scale.x, scale.y * (2 * uv.y - 1)};
    qshape.normals[i]   = {
        cos(phi) * normal2d.y, sin(phi) * normal2d.y, normal2d.x};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.y};
  }
  for (auto& quad : qshape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  merge_shape_inplace(shape, qshape);
  // bottom
  qshape = make_rect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = uv.x * 2 * pif;
    qshape.positions[i] = {
        uv.y * scale.x * cos(phi), uv.y * scale.x * sin(phi), -scale.y};
    qshape.normals[i]   = {0, 0, -1};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.z};
  }
  for (auto& qquad : qshape.quads) swap(qquad.x, qquad.z);
  merge_shape_inplace(shape, qshape);
  return shape;
}

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
shape_data make_lines(
    int num, int steps, vec2f scale, vec2f uvscale, vec2f rad) {
  auto shape = shape_data{};
  shape.positions.resize((steps + 1) * num);
  shape.normals.resize((steps + 1) * num);
  shape.texcoords.resize((steps + 1) * num);
  shape.radius.resize((steps + 1) * num);
  if (num > 1) {
    for (auto j : range(num)) {
      for (auto i : range(steps + 1)) {
        auto uv = vec2f{i / (float)steps, j / (float)(num - 1)};
        shape.positions[j * (steps + 1) + i] = {
            (uv.x - 0.5f) * scale.x, (uv.y - 0.5f) * scale.y, 0};
        shape.normals[j * (steps + 1) + i]   = {1, 0, 0};
        shape.texcoords[j * (steps + 1) + i] = uv * uvscale;
        shape.radius[j * (steps + 1) + i]    = lerp(rad.x, rad.y, uv.x);
      }
    }
  } else {
    for (auto i : range(steps + 1)) {
      auto uv            = vec2f{i / (float)steps, 0};
      shape.positions[i] = {(uv.x - 0.5f) * scale.x, 0, 0};
      shape.normals[i]   = {1, 0, 0};
      shape.texcoords[i] = uv * uvscale;
      shape.radius[i]    = lerp(rad.x, rad.y, uv.x);
    }
  }

  shape.lines.resize(steps * num);
  for (auto j = 0; j < num; j++) {
    for (auto i = 0; i < steps; i++) {
      shape.lines[j * steps + i] = {
          j * (steps + 1) + i, j * (steps + 1) + i + 1};
    }
  }

  return shape;
}

// Make point primitives. Returns points, pos, norm, texcoord, radius.
shape_data make_point(float radius) {
  auto shape      = shape_data{};
  shape.points    = {0};
  shape.positions = {{0, 0, 0}};
  shape.normals   = {{0, 0, 1}};
  shape.texcoords = {{0, 0}};
  shape.radius    = {radius};
  return shape;
}

// Generate a point set with points placed at the origin with texcoords
// varying along u.
shape_data make_points(int num, float uvscale, float radius) {
  auto shape = shape_data{};
  shape.points.resize(num);
  for (auto i : range(num)) shape.points[i] = i;
  shape.positions.assign(num, {0, 0, 0});
  shape.normals.assign(num, {0, 0, 1});
  shape.texcoords.assign(num, {0, 0});
  shape.radius.assign(num, radius);
  for (auto i : range(shape.texcoords.size()))
    shape.texcoords[i] = {(float)i / (float)num, 0};
  return shape;
}

shape_data make_points(vec2i steps, vec2f size, vec2f uvscale, vec2f radius) {
  auto shape  = make_rect(steps, size, uvscale);
  shape.quads = {};
  shape.points.resize(shape.positions.size());
  for (auto i : range(shape.positions.size())) shape.points[i] = (int)i;
  shape.radius.resize(shape.positions.size());
  for (auto i : range(shape.texcoords.size())) {
    shape.radius[i] = lerp(
        radius.x, radius.y, shape.texcoords[i].y / uvscale.y);
  }
  return shape;
}

shape_data make_random_points(
    int num, vec3f size, float uvscale, float radius, uint64_t seed) {
  auto shape = make_points(num, uvscale, radius);
  auto rng   = make_rng(seed);
  for (auto& position : shape.positions)
    position = (2 * rand3f(rng) - 1) * size;
  for (auto& texcoord : shape.texcoords) texcoord = rand2f(rng);
  return shape;
}

// Predefined meshes
shape_data make_monkey(int subdivisions, float scale) {
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}
shape_data make_quad(int subdivisions, float scale) {
  static const auto quad_positions = vector<vec3f>{
      {-1, -1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, +1, 0}};
  static const auto quad_normals = vector<vec3f>{
      {0, 0, 1}, {0, 0, 1}, {0, 0, 1}, {0, 0, 1}};
  static const auto quad_texcoords = vector<vec2f>{
      {0, 0}, {1, 0}, {1, 1}, {0, 1}};
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}
shape_data make_quady(int subdivisions, float scale) {
  static const auto quady_positions = vector<vec3f>{
      {-1, 0, -1}, {-1, 0, +1}, {+1, 0, +1}, {+1, 0, -1}};
  static const auto quady_normals = vector<vec3f>{
      {0, 1, 0}, {0, 1, 0}, {0, 1, 0}, {0, 1, 0}};
  static const auto quady_texcoords = vector<vec2f>{
      {0, 1}, {1, 1}, {1, 0}, {0, 0}};
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}
shape_data make_cube(int subdivisions, float scale) {
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
  static const auto cube_texcoords = vector<vec2f>{{0, 0}, {1, 0}, {1, 1},
      {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1},
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
      {1, 0}, {1, 1}, {0, 1}};
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}
fvshape_data make_fvcube(int subdivisions, float scale) {
  static const auto fvcube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
      {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
      {+1, +1, -1}};
  static const auto fvcube_normals   = vector<vec3f>{{0, 0, +1}, {0, 0, +1},
        {0, 0, +1}, {0, 0, +1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1},
        {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {-1, 0, 0}, {-1, 0, 0},
        {-1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0},
        {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
  static const auto fvcube_texcoords = vector<vec2f>{{0, 0}, {1, 0}, {1, 1},
      {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1},
      {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
      {1, 0}, {1, 1}, {0, 1}};
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}
shape_data make_geosphere(int subdivisions, float scale) {
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
  if (scale != 1) {
    for (auto& p : shape.positions) p *= scale;
  }
  return shape;
}

// Make a hair ball around a shape
shape_data make_hair(const shape_data& base, vec2i steps, vec2f len, vec2f rad,
    vec2f noise, vec2f clump, vec2f rotation, int seed) {
  auto points    = sample_shape(base, steps.y, seed);
  auto bpos      = vector<vec3f>{};
  auto bnorm     = vector<vec3f>{};
  auto btexcoord = vector<vec2f>{};
  for (auto& point : points) {
    bpos.push_back(eval_position(base, point.element, point.uv));
    bnorm.push_back(eval_normal(base, point.element, point.uv));
    btexcoord.push_back(eval_texcoord(base, point.element, point.uv));
  }

  auto rng  = make_rng(seed, 3);
  auto blen = vector<float>(bpos.size());
  for (auto& l : blen) {
    l = lerp(len.x, len.y, rand1f(rng));
  }

  auto cidx = vector<int>();
  if (clump.x > 0) {
    for (auto bidx : range((int)bpos.size())) {
      cidx.push_back(0);
      auto cdist = flt_max;
      for (auto c : range((int)clump.y)) {
        auto d = length(bpos[bidx] - bpos[c]);
        if (d < cdist) {
          cdist       = d;
          cidx.back() = c;
        }
      }
    }
  }

  auto shape = make_lines(steps.y, steps.x, {1, 1}, {1, 1}, {1, 1});
  for (auto i : range((int)shape.positions.size())) {
    auto u             = shape.texcoords[i].x;
    auto bidx          = i / (steps.x + 1);
    shape.positions[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
    shape.normals[i]   = bnorm[bidx];
    shape.radius[i]    = lerp(rad.x, rad.y, u);
    if (clump.x > 0) {
      shape.positions[i] =
          shape.positions[i] +
          (shape.positions[i + (cidx[bidx] - bidx) * (steps.x + 1)] -
              shape.positions[i]) *
              u * clump.x;
    }
    if (noise.x > 0) {
      auto nx =
          (gradient_noise(shape.positions[i] * noise.y + vec3f{0, 0, 0}) * 2 -
              1) *
          noise.x;
      auto ny =
          (gradient_noise(shape.positions[i] * noise.y + vec3f{3, 7, 11}) * 2 -
              1) *
          noise.x;
      auto nz =
          (gradient_noise(shape.positions[i] * noise.y + vec3f{13, 17, 19}) *
                  2 -
              1) *
          noise.x;
      shape.positions[i] += {nx, ny, nz};
    }
  }

  if (clump.x > 0 || noise.x > 0 || rotation.x > 0) {
    shape.normals = lines_tangents(shape.lines, shape.positions);
  }

  return shape;
}

// Grow points around a shape
shape_data make_random_points(
    const shape_data& shape, int num, float radius, uint64_t seed) {
  auto samples = sample_shape(shape, num, seed);
  auto points  = make_points(num, 1, radius);
  for (auto idx : range(num)) {
    points.positions[idx] = eval_position(
        shape, samples[idx].element, samples[idx].uv);
    points.normals[idx] = eval_normal(
        shape, samples[idx].element, samples[idx].uv);
  }
  return points;
}

// Grow hairs around a shape
shape_data make_random_hairs(const shape_data& shape, int num, int steps,
    vec2f len, vec2f radius, float noise, float gravity, uint64_t seed) {
  auto samples = sample_shape(shape, num, seed);
  auto hairs   = make_lines(num, steps, {1, 1}, {1, 1}, radius);
  auto rng     = make_rng(seed);
  for (auto idx : range(num)) {
    auto offset   = idx * (steps + 1);
    auto position = eval_position(shape, samples[idx].element, samples[idx].uv);
    auto direction = eval_normal(shape, samples[idx].element, samples[idx].uv);
    auto length    = lerp(len.x, len.y, rand1f(rng));
    hairs.positions[offset] = position;
    for (auto iidx = 1; iidx <= steps; iidx++) {
      hairs.positions[offset + iidx] = position;
      hairs.positions[offset + iidx] += direction * length / (float)steps;
      hairs.positions[offset + iidx] += (2 * rand3f(rng) - 1) * noise;
      hairs.positions[offset + iidx] += vec3f{0, -gravity, 0};
      direction = normalize(hairs.positions[offset + iidx] - position);
      position  = hairs.positions[offset + iidx];
    }
  }

  hairs.normals = lines_tangents(hairs.lines, hairs.positions);

  return hairs;
}

// Grow hairs around a shape
shape_data make_hair2(const shape_data& base, vec2i steps, vec2f len,
    vec2f radius, float noise, float gravity, int seed) {
  auto points     = sample_shape(base, steps.y, seed);
  auto bpositions = vector<vec3f>{};
  auto bnormals   = vector<vec3f>{};
  auto btexcoord  = vector<vec2f>{};
  for (auto& point : points) {
    bpositions.push_back(eval_position(base, point.element, point.uv));
    bnormals.push_back(eval_normal(base, point.element, point.uv));
    btexcoord.push_back(eval_texcoord(base, point.element, point.uv));
  }

  auto shape = make_lines(steps.y, steps.x, {1, 1}, {1, 1}, radius);
  auto rng   = make_rng(seed);
  for (auto idx : range(steps.y)) {
    auto offset             = idx * (steps.x + 1);
    auto position           = bpositions[idx];
    auto direction          = bnormals[idx];
    auto length             = rand1f(rng) * (len.y - len.x) + len.x;
    shape.positions[offset] = position;
    for (auto iidx = 1; iidx <= steps.x; iidx++) {
      shape.positions[offset + iidx] = position;
      shape.positions[offset + iidx] += direction * length / (float)steps.x;
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
shape_data make_heightfield(vec2i size, const vector<float>& height) {
  auto shape = make_recty({size.x - 1, size.y - 1},
      vec2f{(float)size.x, (float)size.y} / (float)max(size), {1, 1});
  for (auto j : range(size.y))
    for (auto i : range(size.x))
      shape.positions[j * size.x + i].y = height[j * size.x + i];
  shape.normals = quads_normals(shape.quads, shape.positions);
  return shape;
}
shape_data make_heightfield(vec2i size, const vector<vec4f>& color) {
  auto shape = make_recty({size.x - 1, size.y - 1},
      vec2f{(float)size.x, (float)size.y} / (float)max(size), {1, 1});
  for (auto j : range(size.y))
    for (auto i : range(size.x))
      shape.positions[j * size.x + i].y = mean(xyz(color[j * size.x + i]));
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
    auto sphere = make_sphere(steps, scale, 1);
    for (auto& position : sphere.positions) position += vertex;
    merge_shape_inplace(shape, sphere);
  }
  return shape;
}
shape_data polyline_to_cylinders(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = shape_data{};
  for (auto idx = 0; idx < (int)vertices.size() - 1; idx++) {
    auto cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1}, {1, 1, 1});
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
    auto cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1}, {1, 1, 1});
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
  for (auto& line : lines) {
    auto cylinder = make_uvcylinder({steps, 1, 1}, {scale, 1}, {1, 1, 1});
    auto frame    = frame_fromz((positions[line.x] + positions[line.y]) / 2,
           positions[line.x] - positions[line.y]);
    auto length   = distance(positions[line.x], positions[line.y]);
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