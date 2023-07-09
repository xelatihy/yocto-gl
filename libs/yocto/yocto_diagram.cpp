//
// Implementation for Yocto/Diagram
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2023 Fabio Pellacini
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

#include "yocto_diagram.h"

#include <yocto/yocto_sceneio.h>

#include <unordered_set>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Init diagram
diagram_data make_diagram() {
  auto diagram = diagram_data{};
  diagram.scenes.reserve(16);
  return diagram;
}

// Shape transformations
shape_data scale_shape(const shape_data& shape, float size) {
  auto tshape = shape;
  for (auto& position : tshape.positions) position *= size;
  return tshape;
}
shape_data transform_shape(
    const shape_data& shape, const frame3f& frame, float size) {
  return transform_shape(shape, frame * scaling_frame(vec3f{size, size, size}));
}

// Shapes
shape_data load_dshape(const string& filename) {
  auto loaded = load_shape(filename);
  auto shape  = shape_data{
       .points    = loaded.points,
       .lines     = loaded.lines,
       .triangles = loaded.triangles,
       .quads     = loaded.quads,
       .positions = loaded.positions,
       .texcoords = loaded.texcoords,
  };
  return shape;
}

// Helper
static vector<vec3f> transform_points(const frame3f& frame,
    const vector<vec3f>& positions, float scale = 1,
    const vec3f& center = {0, 0, 0}) {
  auto tpositions = vector<vec3f>{};
  for (auto& position : positions)
    tpositions.push_back(transform_point(frame, position * scale + center));
  return tpositions;
}

// Helper
template <typename T>
struct zip_labels {
  struct iterator {
    iterator(
        size_t idx, const vector<string>& labels, const vector<T>& others) :
        idx{idx}, labels{labels}, others{others} {}

    bool operator!=(const iterator& other) const { return idx != other.idx; }
    void operator++() { idx++; }
    pair<const string&, const T&> operator*() const {
      return {labels[idx], others[idx]};
    }

   private:
    size_t                idx;
    const vector<string>& labels;
    const vector<T>&      others;
  };

  zip_labels(const vector<string>& labels, const vector<T>& others) :
      labels{labels}, others{others} {}
  iterator begin() { return {(size_t)0, labels, others}; }
  iterator end() {
    return {std::min(labels.size(), others.size()), labels, others};
  }

 private:
  const vector<string>& labels;
  const vector<T>&      others;
};

// Automatically set offsets
static void update_offsets(diagram_data& diagram) {
  auto total_size = vec2f{0, 0};
  for (auto& scene : diagram.scenes) total_size += scene.size + scene.margin;
  auto offset = vec2f{-total_size.x / 2, 0};
  for (auto& scene : diagram.scenes) {
    scene.offset.x = offset.x + scene.size.x / 2 + scene.margin.x / 2;
    scene.offset.y = 0;
    offset.x += scene.size.x + scene.margin.x;
  }
}

// Scene
diagram_scene& add_scene(
    diagram_data& diagram, const diagram_scene_params& params) {
  auto& scene  = diagram.scenes.emplace_back();
  scene.size   = params.size;
  scene.margin = params.margin;
  scene.frame  = dtranslation({-params.center.x, -params.center.y, 0}) *
                params.frame;

  if (!params.title.empty()) {
    scene.objects.push_back({
        .labels     = {params.title + "!!t"},
        .lpositions = {{0, scene.size.y / 2, 0}},
    });
  }

  if (!params.subtitle.empty()) {
    scene.objects.push_back({
        .labels     = {params.subtitle + "!!b"},
        .lpositions = {{0, -scene.size.y / 2, 0}},
    });
  }

  update_offsets(diagram);

  return scene;
}

// Add labels
void add_labels(diagram_scene& diagram, const diagram_labels_params& params) {
  auto positions = transform_points(
      diagram.frame * params.frame, params.positions);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }

  diagram.objects.push_back({
      .labels     = labels,
      .lpositions = lpositions,
  });
}

// Add points
void add_points(diagram_scene& diagram, const diagram_points_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);

  auto points = vector<int>{};
  for (auto idx : range((int)params.positions.size())) points.push_back(idx);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }

  diagram.objects.push_back({
      .points     = points,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add lines
void add_lines(diagram_scene& diagram, const diagram_lines_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);

  auto lines = vector<vec2i>{};
  if (params.lines.empty()) {
    for (auto idx : range((int)params.positions.size() / 2))
      lines.push_back({idx * 2 + 0, idx * 2 + 1});
  } else {
    lines = params.lines;
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, lines)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
      .connect    = params.connector ? 0.3f : 0.0f,
  });
}

// Add arrows
void add_arrows(diagram_scene& diagram, const diagram_arrows_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);

  auto arrows = vector<vec2i>{};
  if (params.arrows.empty()) {
    for (auto idx : range((int)params.positions.size() / 2))
      arrows.push_back({idx * 2 + 0, idx * 2 + 1});
  } else {
    arrows = params.arrows;
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, arrows)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .arrows     = arrows,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
      .connect    = params.connector ? 0.3f : 0.0f,
  });
}

// Add vectors
void add_vectors(diagram_scene& diagram, const diagram_vectors_params& params) {
  auto positions_ = params.vectors;
  positions_.push_back({0, 0, 0});
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto arrows = vector<vec2i>{};
  for (auto idx : range((int)params.vectors.size()))
    arrows.push_back({(int)params.vectors.size(), idx});

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, arrows)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .arrows     = arrows,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add axes
void add_axes(diagram_scene& diagram, const diagram_axes_params& params) {
  auto positions_ = vector<vec3f>{params.axes.o,
      params.axes.o + params.aspect.x * params.axes.x,
      params.axes.o + params.aspect.y * params.axes.y,
      params.axes.o + params.aspect.z * params.axes.z};
  auto positions  = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto arrows = vector<vec2i>{{0, 1}, {0, 2}, {0, 3}};

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, arrows)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .arrows     = arrows,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add rays
void add_rays(diagram_scene& diagram, const diagram_rays_params& params) {
  auto positions = vector<vec3f>{};
  for (auto& ray : params.rays) {
    positions.push_back(ray.o);
    positions.push_back(ray.o + ray.d);
  }

  add_arrows(diagram, {.frame        = params.frame,
                          .positions = positions,
                          .position  = params.position,
                          .scale     = params.scale,
                          .connector = params.connector,
                          .labels    = params.labels,
                          .clabels   = params.clabels,
                          .stroke    = params.stroke,
                          .thickness = params.thickness});

  if (params.llength <= 0.001f) return;

  auto lpositions = vector<vec3f>{};
  for (auto& ray : params.rays) {
    lpositions.push_back(ray.o + ray.d);
    lpositions.push_back(ray.o + ray.d * params.llength);
  }

  add_lines(diagram, {.frame        = params.frame,
                         .positions = lpositions,
                         .position  = params.position,
                         .scale     = params.scale,
                         .connector = true,
                         .labels    = params.llabels,
                         .clabels   = params.lclabels,
                         .stroke    = params.lstroke,
                         .thickness = params.lthickness});
}

// Add quads
void add_quads(diagram_scene& diagram, const diagram_quads_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);
  auto texcoords = params.texcoords;

  auto quads = vector<vec4i>{};
  if (params.quads.empty()) {
    for (auto idx : range((int)params.positions.size() / 4))
      quads.push_back({idx * 4 + 0, idx * 4 + 1, idx * 4 + 2, idx * 4 + 3});
  } else {
    quads = params.quads;
  }

  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, quad] : zip_labels(params.clabels, quads)) {
    labels.push_back(label);
    lpositions.push_back((positions[quad.x] + positions[quad.y] +
                             positions[quad.z] + positions[quad.w]) /
                         4);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .texcoords  = texcoords,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add triangles
void add_triangles(
    diagram_scene& diagram, const diagram_triangles_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);
  auto texcoords = params.texcoords;

  auto triangles = vector<vec3i>{};
  if (params.triangles.empty()) {
    for (auto idx : range((int)params.positions.size() / 3))
      triangles.push_back({idx * 3 + 0, idx * 3 + 1, idx * 3 + 2});
  } else {
    triangles = params.triangles;
  }

  auto lines = get_edges(triangles);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, triangle] : zip_labels(params.clabels, triangles)) {
    labels.push_back(label);
    lpositions.push_back((positions[triangle.x] + positions[triangle.y] +
                             positions[triangle.z]) /
                         3);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .triangles  = triangles,
      .positions  = positions,
      .texcoords  = texcoords,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add polyline
void add_polyline(
    diagram_scene& diagram, const diagram_polyline_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      params.positions, params.scale, params.position);

  auto lines = vector<vec2i>{};
  for (auto idx : range((int)params.positions.size() - 1))
    lines.push_back({idx + 0, idx + 1});

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, lines)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add polygon
void add_polygon(diagram_scene& diagram, const diagram_polygon_params& params) {
  auto bbox = invalidb3f;
  for (auto& position : params.positions) bbox = merge(bbox, position);
  auto center = (bbox.min + bbox.max) / 2;

  auto steps      = (int)params.positions.size();
  auto positions_ = params.positions;
  positions_.push_back(center);
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto triangles = vector<vec3i>{};
  for (auto idx : range(steps)) {
    triangles.push_back({idx, (idx + 1) % steps, steps});
  }

  auto lines = vector<vec2i>{};
  for (auto idx : range(steps)) {
    lines.push_back({idx, (idx + 1) % steps});
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  if (!params.clabels.empty()) {
    labels.push_back(params.clabels.front());
    lpositions.push_back(positions.back());
  }

  diagram.objects.push_back({
      .lines      = lines,
      .triangles  = triangles,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add rect
void add_rect(diagram_scene& diagram, const diagram_rect_params& params) {
  auto positions_ = dconstants::quad_positions;
  auto scale      = vec3f{params.aspect.x / params.aspect.y, 1, 1};
  for (auto& position : positions_) position *= scale;
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);
  auto texcoords = dconstants::quad_texcoords;

  auto quads = dconstants::quad_quads;
  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, quad] : zip_labels(params.clabels, quads)) {
    labels.push_back(label);
    lpositions.push_back((positions[quad.x] + positions[quad.y] +
                             positions[quad.z] + positions[quad.w]) /
                         4);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .texcoords  = texcoords,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add disk
void add_disk(diagram_scene& diagram, const diagram_disk_params& params) {
  auto positions_ = vector<vec3f>{};
  for (auto idx : range(params.steps)) {
    auto theta = 2 * pif * idx / (float)params.steps;
    positions_.push_back({cos(theta), sin(theta), 0});
  }
  positions_.push_back(params.position);
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto triangles = vector<vec3i>{};
  for (auto idx : range(params.steps)) {
    triangles.push_back({idx, (idx + 1) % params.steps, params.steps});
  }

  auto lines = vector<vec2i>{};
  for (auto idx : range(params.steps)) {
    lines.push_back({idx, (idx + 1) % params.steps});
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  if (!params.clabels.empty()) {
    labels.push_back(params.clabels.front());
    lpositions.push_back(positions.back());
  }

  diagram.objects.push_back({
      .lines      = lines,
      .triangles  = triangles,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add arc
void add_arc(diagram_scene& diagram, const diagram_arc_params& params) {
  // refence frame
  auto x     = normalize(params.from - params.position);
  auto y     = normalize(params.to - params.position);
  auto angle = acos(dot(x, y));
  auto z     = normalize(cross(x, y));
  y          = normalize(cross(z, x));
  auto frame = frame3f{x, y, z, {0, 0, 0}};

  auto positions_ = vector<vec3f>{};
  for (auto idx : range(params.steps + 1)) {
    auto theta = angle * idx / (float)params.steps;
    positions_.push_back(transform_point(frame, {cos(theta), sin(theta), 0}));
  }
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto lines = vector<vec2i>{};
  for (auto idx : range(params.steps)) {
    lines.push_back({idx, idx + 1});
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  if (!params.clabels.empty()) {
    labels.push_back(params.clabels.front());
    lpositions.push_back(positions[positions.size() / 2]);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add qbezier
void add_qbeziers(
    diagram_scene& diagram, const diagram_qbeziers_params& params) {
  auto qbezier_point = [](vec3f p1, vec3f p2, vec3f p3, float u) {
    return (1 - u) * (1 - u) * p1 + 2 * (1 - u) * u * p2 + u * u * p3;
  };
  auto positions = vector<vec3f>();
  for (auto segment : range(params.positions.size() / 3)) {
    for (auto idx : range(params.steps + 1)) {
      positions.push_back(qbezier_point(params.positions[0 + segment * 3],
          params.positions[1 + segment * 3], params.positions[2 + segment * 3],
          idx / (float)params.steps));
    }
  }

  add_polyline(diagram, {.frame        = params.frame,
                            .positions = positions,
                            .position  = params.position,
                            .scale     = params.scale,
                            .labels    = params.labels,
                            .clabels   = params.clabels,
                            .stroke    = params.stroke,
                            .thickness = params.thickness});
}

// Add bezier
void add_beziers(diagram_scene& diagram, const diagram_beziers_params& params) {
  auto bezier_point = [](vec3f p1, vec3f p2, vec3f p3, vec3f p4, float u) {
    return (1 - u) * (1 - u) * (1 - u) * p1 + 3 * (1 - u) * (1 - u) * u * p2 +
           3 * (1 - u) * u * u * p3 + u * u * u * p4;
  };
  auto positions = vector<vec3f>();
  for (auto segment : range(params.positions.size() / 4)) {
    for (auto idx : range(params.steps + 1)) {
      positions.push_back(bezier_point(params.positions[0 + segment * 4],
          params.positions[1 + segment * 4], params.positions[2 + segment * 4],
          params.positions[3 + segment * 4], idx / (float)params.steps));
    }
  }

  add_polyline(diagram, {.frame        = params.frame,
                            .positions = positions,
                            .position  = params.position,
                            .scale     = params.scale,
                            .labels    = params.labels,
                            .clabels   = params.clabels,
                            .stroke    = params.stroke,
                            .thickness = params.thickness});
}

// Add measure
void add_measure(diagram_scene& diagram, const diagram_measure_params& params) {
  // refence frame
  auto a = params.from + params.offset, b = params.to + params.offset;

  auto positions_ = vector<vec3f>{a, b};
  auto positions  = transform_points(diagram.frame * params.frame, positions_);

  auto lines = vector<vec2i>{{0, 1}};

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  if (!params.clabels.empty()) {
    labels.push_back(params.clabels.front());
    lpositions.push_back((a + b) / 2);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add cube
void add_cube(diagram_scene& diagram, const diagram_cube_params& params) {
  auto positions = transform_points(diagram.frame * params.frame,
      dconstants::cube_positions, params.scale, params.position);

  auto quads = dconstants::cube_quads;

  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add cube
void add_sphere(diagram_scene& diagram, const diagram_sphere_params& params) {
  auto shape = params.uvsphere ? make_uvsphere({params.steps * 2, params.steps})
                               : make_sphere(params.steps);
  auto positions = transform_points(diagram.frame * params.frame,
      shape.positions, params.scale, params.position);
  // auto normals = transform_normals(diagram.frame * params.frame,
  // shape.normals);
  auto texcoords = shape.texcoords;

  auto quads = shape.quads;

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }

  diagram.objects.push_back({
      .quads     = quads,
      .positions = positions,
      // .normals    = normals,
      .texcoords  = texcoords,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
      .texture    = params.texture,
  });
}

// Add grid
void add_grid(diagram_scene& diagram, const diagram_grid_params& params) {
  auto ratio     = (float)params.steps.x / (float)params.steps.y;
  auto shape     = make_rect({params.steps.x, params.steps.y}, {ratio, 1});
  auto positions = transform_points(diagram.frame * params.frame,
      shape.positions, params.scale, params.position);

  auto quads = shape.quads;
  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, quad] : zip_labels(params.clabels, quads)) {
    labels.push_back(label);
    lpositions.push_back((positions[quad.x] + positions[quad.y] +
                             positions[quad.z] + positions[quad.w]) /
                         4);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add affine grid
void add_affinegrid(
    diagram_scene& diagram, const diagram_affinegrid_params& params) {
  auto ratio      = (float)params.steps.x / (float)params.steps.y;
  auto shape      = make_rect({params.steps.x, params.steps.y}, {ratio, 1});
  auto positions_ = shape.positions;
  for (auto& position : positions_)
    position = position.x * params.axes_a + position.y * params.axes_b;
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto quads = shape.quads;
  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, quad] : zip_labels(params.clabels, quads)) {
    labels.push_back(label);
    lpositions.push_back((positions[quad.x] + positions[quad.y] +
                             positions[quad.z] + positions[quad.w]) /
                         4);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
  });
}

// Add image
void add_image(diagram_scene& diagram, const diagram_image_params& params) {
  auto ratio = (float)params.image.extent(0) / (float)params.image.extent(1);
  auto scale = !params.scaley ? vec3f{ratio, 1, 1} : vec3f{1, 1 / ratio, 1};
  auto positions_ = dconstants::quad_positions;
  for (auto& position : positions_) position *= scale;
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto texcoords = dconstants::quad_texcoords;

  auto quads = dconstants::quad_quads;

  auto lines = get_edges(quads);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  if (!params.clabels.empty()) {
    labels.push_back(params.clabels.front());
    lpositions.push_back(
        (positions[0] + positions[1] + positions[2] + positions[3]) / 4);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .quads      = quads,
      .positions  = positions,
      .texcoords  = texcoords,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .fill       = params.fill,
      .thickness  = params.thickness,
      .texture    = params.image,
      .nearest    = !params.interpolate,
  });

  if (!params.olabels.empty()) {
    auto oscale     = vec3f{3, 1, 1} * 0.24f;
    auto ocenter    = positions[1] + vec3f{-0.9, +0.4, +0.1};
    auto opositions = dconstants::quad_positions;
    for (auto& oposition : opositions) oposition = oposition * oscale + ocenter;
    diagram.objects.push_back({
        .quads      = dconstants::quad_quads,
        .positions  = opositions,
        .labels     = {params.olabels.front()},
        .lpositions = {ocenter},
        .stroke     = dcolors::transparent,
        .fill       = dcolors::black,
        .text       = dcolors::white,
    });
  }
}

// Add image grid
void add_imagegrid(
    diagram_scene& diagram, const diagram_imagegrid_params& params) {
  auto ratio = (float)params.image.extent(0) / (float)params.image.extent(1);
  add_image(diagram, {.frame          = params.frame,
                         .image       = params.image,
                         .position    = params.position,
                         .scale       = params.scale,
                         .interpolate = params.interpolate,
                         .fill        = params.fill});
  add_grid(diagram, {.frame        = params.frame,
                        .position  = params.position,
                        .steps     = (vec2i)params.image.extents(),
                        .scale     = params.scale,
                        .labels    = params.labels,
                        .clabels   = params.clabels,
                        .stroke    = params.stroke,
                        .fill      = dcolors::transparent,
                        .thickness = params.thickness});
}

// Add image mosaic
void add_imagemosaic(
    diagram_scene& diagram, const diagram_imagemosaic_params& params) {
  if (params.images.size() != 4)
    throw std::out_of_range{"supports only 4 images for now"};
  auto ratio = (float)params.images.at(0).extent(0) /
               (float)params.images.at(0).extent(1);
  auto positions = vector<vec3f>{{-ratio / 2, 0.5, 0}, {+ratio / 2, 0.5, 0},
      {-ratio / 2, -0.5, 0}, {+ratio / 2, -0.5, 0}};
  for (auto idx : range((int)size(params.images))) {
    add_image(
        diagram, {.frame          = params.frame,
                     .image       = params.images[idx],
                     .position    = positions[idx] * params.scale,
                     .scale       = params.scale / 2,
                     .interpolate = params.interpolate,
                     .labels      = params.labels,
                     .clabels     = idx < params.clabels.size()
                                        ? vector<string>{params.clabels.at(idx)}
                                        : vector<string>{},
                     .olabels     = idx < params.olabels.size()
                                        ? vector<string>{params.olabels.at(idx)}
                                        : vector<string>{},
                     .stroke      = params.stroke,
                     .fill        = params.fill,
                     .thickness   = params.thickness});
  }
}

// Add random points
void add_randompoints(
    diagram_scene& diagram, const diagram_randompoints_params& params) {
  auto rng        = make_rng(187291);
  auto positions_ = vector<vec3f>{};
  auto ssteps     = int(round(sqrt((float)params.steps)));
  for (auto idx : range(params.steps)) {
    auto uv = !params.stratified ? rand2f(rng)
                                 : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                       (idx / ssteps + rand1f(rng)) / ssteps};
    switch (params.type) {
      case diagram_randompoints_type::quad: {
        positions_.push_back({2 * uv.x - 1, 2 * uv.y - 1, 0});
      } break;
      case diagram_randompoints_type::disk: {
        positions_.push_back({sample_disk(uv).x, sample_disk(uv).y, 0});
      } break;
      case diagram_randompoints_type::disknu: {
        auto uv_ = vec2f{uv.x, uv.y * uv.y};
        positions_.push_back({sample_disk(uv_).x, sample_disk(uv_).y, 0});
      } break;
      case diagram_randompoints_type::triangle: {
        positions_.push_back(interpolate_triangle(params.triangle[0],
            params.triangle[1], params.triangle[2], sample_triangle(uv)));
      } break;
    }
  }
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto points = vector<int>{};
  for (auto idx : range(params.steps)) points.push_back(idx);

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }

  diagram.objects.push_back({
      .points     = points,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add random lines
void add_randomlines(
    diagram_scene& diagram, const diagram_randomlines_params& params) {
  auto rng        = make_rng(187291);
  auto ssteps     = int(round(sqrt((float)params.steps)));
  auto positions_ = vector<vec3f>{};
  for (auto idx : range(params.steps)) {
    auto uv = !params.stratified ? rand2f(rng)
                                 : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                       (idx / ssteps + rand1f(rng)) / ssteps};
    switch (params.type) {
      case diagram_randomlines_type::hemi:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere({0, 0, 1}, uv));
        break;
      case diagram_randomlines_type::hemicos:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere_cos({0, 0, 1}, uv));
        break;
      case diagram_randomlines_type::hemicospower:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere_cospower(64, {0, 0, 1}, uv));
        break;
      case diagram_randomlines_type::beam:
        auto rdisk = sample_disk(uv);
        positions_.push_back({rdisk.x, rdisk.y, 0});
        positions_.push_back({rdisk.x, rdisk.y, 1});
        break;
    }
  }
  auto positions = transform_points(
      diagram.frame * params.frame, positions_, params.scale, params.position);

  auto lines = vector<vec2i>{};
  for (auto idx : range(params.steps)) {
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  }

  auto labels     = vector<string>{};
  auto lpositions = vector<vec3f>{};
  for (auto&& [label, position] : zip_labels(params.labels, positions)) {
    labels.push_back(label);
    lpositions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(params.clabels, lines)) {
    labels.push_back(label);
    lpositions.push_back((positions[line.x] + positions[line.y]) / 2);
  }

  diagram.objects.push_back({
      .lines      = lines,
      .positions  = positions,
      .labels     = labels,
      .lpositions = lpositions,
      .stroke     = params.stroke,
      .thickness  = params.thickness,
  });
}

// Add plot
diagram_scene& add_plot(
    diagram_data& diagram, const diagram_plot_params& params) {
  return add_scene(diagram, {.size           = params.size,
                                .center      = params.center,
                                .title       = params.title,
                                .subtitle    = params.subtitle,
                                .frame       = params.frame,
                                .margin      = params.margin,
                                .spacing     = params.spacing,
                                .auto_offset = params.auto_offset});
}

// Add plot axes
diagram_scene& add_plotaxes(
    diagram_scene& diagram, const diagram_plotaxes_params& params) {
  // add frame
  auto ratio = diagram.size.x / diagram.size.y;
  add_quads(diagram, {.positions = transform_points(
                          dscaling({diagram.size.x / 2, diagram.size.y / 2, 1}),
                          dconstants::quad_positions),
                         .quads     = dconstants::quad_quads,
                         .stroke    = params.stroke,
                         .fill      = dcolors::transparent,
                         .thickness = params.thickness});
  if (!params.title.empty()) {
    add_labels(diagram, {.positions = {{0, diagram.size.y, 0}},
                            .labels = {params.title + "!!t"}});
  }

  auto min    = vec2f{params.xbounds.x, params.ybounds.x};
  auto max    = vec2f{params.xbounds.y, params.ybounds.y};
  auto center = (max + min) / 2, size = (max - min);
  auto pframe = dscaling(
                    {diagram.size.x / size.x, diagram.size.y / size.y, 1.0f}) *
                dtranslation({-center.x, -center.y, 0.0f});

  for (auto&& [label, tick] : zip_labels(params.xlabels, params.xticks)) {
    add_labels(diagram, {.positions = {{transform_point(pframe, {tick, 0, 0}).x,
                             -diagram.size.y / 2, 0}},
                            .labels = {label + "!!b"}});
  }
  for (auto&& [label, tick] : zip_labels(params.ylabels, params.yticks)) {
    add_labels(diagram, {.positions = {{-diagram.size.x / 2,
                             transform_point(pframe, {0, tick, 0}).y, 0}},
                            .labels = {label + "!!l"}});
  }

  diagram.frame = diagram.frame * pframe;

  return diagram;
}

// Plot functions
void add_plotdata(
    diagram_scene& diagram, const diagram_plotdata_params& params) {
  auto& points = !params.points.empty() ? params.points
                                        : sample_function(params.function,
                                              params.range, params.steps);

  auto clip = invalidb2f;  // TODO: save on context
  if (clip == invalidb2f) {
    auto positions = vector<vec3f>{};
    for (auto point : points) positions.push_back({point.x, point.y, 0});
    return add_polyline(diagram, {.positions    = positions,
                                     .labels    = params.labels,
                                     .stroke    = params.stroke,
                                     .thickness = params.thickness});
  } else {
    auto clamp = [](vec2f x, vec2f min, vec2f max) -> vec2f {
      return {yocto::clamp(x.x, min.x, max.x), yocto::clamp(x.y, min.y, max.y)};
    };
    auto positions = vector<vec3f>{};
    for (auto idx : range((int)points.size() - 1)) {
      auto point1 = points[idx + 0], point2 = points[idx + 1];
      auto clipped1 = clamp(point1, clip.min, clip.max),
           clipped2 = clamp(point2, clip.min, clip.max);
      if (point1 != clipped1 && point2 != clipped2) continue;
      positions.push_back({clipped1.x, clipped1.y, 0});
      positions.push_back({clipped2.x, clipped2.y, 0});
    }
    return add_lines(diagram, {.positions    = positions,
                                  .labels    = params.labels,
                                  .stroke    = params.stroke,
                                  .thickness = params.thickness});
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM RENDERING
// -----------------------------------------------------------------------------
namespace yocto {

// Make default camera
// 100 : 36 = 25 :  x => x = 36 * 25 / 100 = 9
// 100 : 36 =  d : 10 => d = 100 / 36 * 10 = 27.77777
//   l : 36 = 25 : 10 => l = 25 * 36 / 10 = 90
static camera_data make_default_camera() {
  return {
      .frame  = lookat_frame(vec3f{0, 0, 25}, vec3f{0, 0, 0}, vec3f{0, 1, 0}),
      .lens   = 0.09f,
      .aspect = 16.0f / 9.0f,
      .focus  = 25,
  };
}

// Render text
static texture_data make_text_texture(const string& text) {
  // Rendered text location
  static string text_location = "./textcache";

  // Latex template
  static string latex_template_header =
      "\\documentclass[10pt]{article}\n"
      "\\usepackage[cochineal]{newtx}\n"
      "\\usepackage[T1]{fontenc}\n"
      "\\usepackage{parskip}\n"
      "\n"
      "\\begin{document}\n"
      "\\pagenumbering{gobble}\n\n";
  static string latex_template_footer =
      "\n\n"
      "\\end{document}\n\n";

  // djb2 algorithm from http://www.cse.yorku.ca/~oz/hash.html
  auto hash = [](const string& text) -> string {
    auto str  = text.c_str();
    auto hash = (unsigned long)5381;
    auto c    = (int)*str++;
    while (c) {
      hash = ((hash << 5) + hash) + c;  // hash * 33 + c
      c    = (int)*str++;
    }
    return std::to_string(hash);
  };

  auto text_name    = hash(text);
  auto texture_name = text_location + "/" + text_name + ".png";
  auto latex_name   = text_location + "/" + text_name + ".tex";
  if (!path_exists(texture_name)) {
    // save latex
    printf("missing text: %s\n", text.c_str());
    auto latex = latex_template_header + text + latex_template_footer;
    save_text(latex_name, latex);
    // render latex
    auto command = "cd " + text_location + " && pdflatex " + text_name + ".tex";
    printf("%s\n", command.c_str());
    system(command.c_str());
    command = "cd " + text_location + " && pdfcrop " + text_name + ".pdf " +
              text_name + "-cropped.pdf";
    printf("%s\n", command.c_str());
    system(command.c_str());
    command = "cd " + text_location +
              " && gs -sDEVICE=pnggray -sBATCH -sOutputFile=" + text_name +
              ".png" + " -dNOPAUSE -r1200 " + text_name + "-cropped.pdf";
    printf("%s\n", command.c_str());
    system(command.c_str());
    system(("rm " + text_location + "/*.aux").c_str());
    system(("rm " + text_location + "/*.fls").c_str());
    system(("rm " + text_location + "/*.log").c_str());
    system(("rm " + text_location + "/*.gz").c_str());
    system(("rm " + text_location + "/*.fdb_latexmk").c_str());
  }

  auto texture = load_texture(texture_name);
  for (auto& c : texture.pixelsb) {
    c = {255, 255, 255, (byte)(255 - c.x)};
  }
  return texture;
}

// Text shape
static shape_data make_text_shape(const string& text, const vec3f& offset,
    const texture_data& texture, float scale) {
  auto extents = max(texture.pixelsf.extents(), texture.pixelsb.extents());
  auto [width, height] = scale * extents / 72.0f;
  auto offset_         = vec3f{scale * 144 / 72.0f, scale * 144 / 72.0f, 1};
  auto z               = 0.1f;  // small offset to avoid z fighting
  auto shape           = shape_data{};
  shape.quads          = {{0, 1, 2, 3}};
  shape.positions      = {{-width, -height, z}, {+width, -height, z},
           {+width, +height, z}, {-width, +height, z}};
  for (auto& p : shape.positions) {
    if (offset.x > 0) p.x += +width + offset_.x * offset.x;
    if (offset.x < 0) p.x += -width + offset_.x * offset.x;
    if (offset.y > 0) p.y += +height + offset_.y * offset.y;
    if (offset.y < 0) p.y += -height + offset_.y * offset.y;
    if (offset.z > 0) p.z += offset.z;
  }
  shape.texcoords = {{0, 1}, {1, 1}, {1, 0}, {0, 0}};
  return shape;
}

// Text offset
static frame3f make_text_frame(const string& text, const texture_data& texture,
    const vec3f& camera, const vec3f& position, float scale) {
  return lookat_frame(position, camera, {0, 1, 0}, true);
}

// Make text material
static material_data make_text_material(const vec4f& color, int textureid) {
  return {
      .type      = material_type::matte,
      .emission  = {color.x, color.y, color.z},
      .opacity   = color.w,
      .color_tex = textureid,
  };
}

// split label
[[maybe_unused]] static pair<string, vec3f> split_label(const string& label) {
  auto pos = label.find("!!");
  if (pos == string::npos) return {label, {0, 0, 0}};
  auto text   = label.substr(0, pos);
  auto offset = vec3f{0, 0, 0};
  auto scale  = 1.0f;
  for (auto c : label.substr(pos + 2)) {
    auto reset_scale = true;
    switch (c) {
      case 'h':
        scale       = 0.5f;
        reset_scale = false;
        break;
      case 'a':
        scale       = 0.001f;
        reset_scale = false;
        break;
      case 'l': offset.x -= scale; break;
      case 'r': offset.x += scale; break;
      case 't': offset.y += scale; break;
      case 'b': offset.y -= scale; break;
      case 'n': offset.z += scale; break;
      case 'f': offset.z -= scale; break;
    }
    if (reset_scale) scale = 1.0f;
  }
  return {text, offset};
}

// Make text labels
static vector<string> make_text_labels(const vector<string>& labels) {
  if (labels.empty()) return {};
  auto texts = vector<string>{};
  for (auto& label : labels) texts.push_back(split_label(label).first);
  return texts;
}

// Make text offsets
static vector<vec3f> make_text_offsets(const vector<string>& labels) {
  if (labels.empty()) return {};
  auto offsets = vector<vec3f>{};
  for (auto& label : labels) offsets.push_back(split_label(label).second);
  return offsets;
}

// Make surface
static shape_data make_triangles_shape(const vector<vec3i>& triangles,
    const vector<vec3f>& positions_, const vector<vec2f>& texcoords,
    const frame3f& frame) {
  auto positions = positions_;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape      = shape_data{};
  nshape.triangles = triangles;
  nshape.positions = positions;
  nshape.texcoords = texcoords;
  return nshape;
}

// Make surface
static shape_data make_quads_shape(const vector<vec4i>& quads,
    const vector<vec3f>& positions_, const vector<vec2f>& texcoords,
    const frame3f& frame) {
  auto positions = positions_;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape      = shape_data{};
  nshape.quads     = quads;
  nshape.positions = positions;
  nshape.texcoords = texcoords;
  return nshape;
}

// Make fill material
static material_data make_fill_material(const vec4f& fill, int textureid) {
  return {
      .type         = material_type::matte,
      .emission     = {fill.x, fill.y, fill.z},
      .opacity      = fill.w,
      .emission_tex = textureid,
  };
}

// Make a shape of lines and points
static shape_data make_points_shape(const vector<int>& points,
    const vector<vec3f>& positions_, const frame3f& frame, float radius) {
  auto make_sphere = [](vec3f position, float radius) {
    auto frame = translation_frame(position) *
                 scaling_frame(vec3f{radius * 3, radius * 3, radius * 3});
    auto sphere = yocto::make_sphere(8);
    for (auto& p : sphere.positions) p = transform_point(frame, p);
    return sphere;
  };

  auto positions = positions_;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& point : points) {
    merge_shape_inplace(nshape, make_sphere(positions[point], radius));
  }
  return nshape;
}

// Make a shape of lines
static shape_data make_lines_shape(const vector<vec2i>& lines,
    const vector<vec3f>& positions_, const frame3f& frame, float radius,
    float connect) {
  auto make_capsule = [](vec3f position1, vec3f position2, float radius) {
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    return capsule;
  };

  auto positions = positions_;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& line : lines) {
    if (connect <= 0) {
      merge_shape_inplace(
          nshape, make_capsule(positions[line.x], positions[line.y], radius));
    } else {
      auto position1 = positions[line.x], position2 = positions[line.y];
      auto middle = (position1 + position2) / 2;
      auto offset = max(0.0f, length(position2 - position1) - connect) *
                    normalize(position2 - position1) / 2;
      merge_shape_inplace(
          nshape, make_capsule(middle - offset, middle + offset, radius));
    }
  }
  return nshape;
}

// Make a shape of arrowa
static shape_data make_arrows_shape(const vector<vec2i>& lines,
    const vector<vec3f>& positions_, const frame3f& frame, float radius,
    float connect) {
  auto make_capsule = [](vec3f position1, vec3f position2, float radius) {
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    return capsule;
  };
  auto make_arrow = [](vec3f position1, vec3f position2, float radius) {
    auto nshape         = shape_data{};
    auto axis           = position2 - position1;
    auto rotation_axis  = normalize(normalize(axis) + vec3f{0, 0, 1});
    auto rotation_angle = pif;
    auto frame          = translation_frame((position1 + position2) / 2) *
                 rotation_frame(rotation_axis, rotation_angle);
    auto capsule = make_uvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    merge_shape_inplace(nshape, capsule);
    auto frame2 = translation_frame(position2) *
                  rotation_frame(rotation_axis, rotation_angle);
    auto cone = make_uvcone({8, 8, 8}, {radius * 3, radius * 5});
    for (auto& p : cone.positions) p = transform_point(frame2, p);
    merge_shape_inplace(nshape, cone);
    return nshape;
  };

  auto positions = positions_;
  for (auto& position : positions) position = transform_point(frame, position);

  auto nshape = shape_data{};
  for (auto& line : lines) {
    if (connect <= 0) {
      merge_shape_inplace(
          nshape, make_arrow(positions[line.x], positions[line.y], radius));
    } else {
      auto position1 = positions[line.x], position2 = positions[line.y];
      auto middle = (position1 + position2) / 2;
      auto offset = max(0.0f, length(position2 - position1) - connect) *
                    normalize(position2 - position1) / 2;
      merge_shape_inplace(
          nshape, make_arrow(middle - offset, middle + offset, radius));
    }
  }
  return nshape;
}

// Make stroke material
static material_data make_stroke_material(
    const vec4f& stroke, int textureid = invalidid) {
  return {
      .type         = material_type::matte,
      .emission     = {stroke.x, stroke.y, stroke.z},
      .opacity      = stroke.w,
      .emission_tex = textureid,
  };
}

// Make an image texture
static texture_data make_image_texture(
    const array2d<vec4f>& image, bool nearest) {
  return texture_data{.pixelsb = float_to_byte(image), .nearest = nearest};
}

// Add an object to a scene
static void add_object(scene_data& scene, const diagram_object& object) {
  // triangles
  if (!object.triangles.empty() && object.fill.w > 0) {
    scene.shapes.push_back(make_triangles_shape(
        object.triangles, object.positions, object.texcoords, object.frame));
    if (object.texture.empty()) {
      scene.materials.push_back(make_fill_material(object.fill, invalidid));
    } else {
      scene.textures.emplace_back(
          make_image_texture(object.texture, object.nearest));
      scene.materials.push_back(
          make_fill_material(object.fill, (int)scene.textures.size() - 1));
    }
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // quads
  if (!object.quads.empty() && object.fill.w > 0) {
    scene.shapes.push_back(make_quads_shape(
        object.quads, object.positions, object.texcoords, object.frame));
    if (object.texture.empty()) {
      scene.materials.push_back(make_fill_material(object.fill, invalidid));
    } else {
      scene.textures.emplace_back(
          make_image_texture(object.texture, object.nearest));
      scene.materials.push_back(
          make_fill_material(object.fill, (int)scene.textures.size() - 1));
    }
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // points
  if (!object.points.empty() && object.stroke.w > 0) {
    scene.shapes.push_back(make_points_shape(
        object.points, object.positions, object.frame, object.thickness));
    scene.materials.push_back(make_stroke_material(object.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // lines
  if (!object.lines.empty() && object.stroke.w > 0) {
    scene.shapes.push_back(make_lines_shape(object.lines, object.positions,
        object.frame, object.thickness, object.connect));
    scene.materials.push_back(make_stroke_material(object.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // arrows
  if (!object.arrows.empty() && object.stroke.w > 0) {
    scene.shapes.push_back(make_arrows_shape(object.arrows, object.positions,
        object.frame, object.thickness, object.connect));
    scene.materials.push_back(make_stroke_material(object.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // text
  if ((!object.labels.empty()) && object.text.w > 0) {
    auto positions = object.lpositions;
    auto offsets   = make_text_offsets(object.labels);
    auto labels    = make_text_labels(object.labels);
    for (auto idx = 0; idx < (int)labels.size(); idx++) {
      if (labels[idx].empty()) continue;
      auto& texture = scene.textures.emplace_back(
          make_text_texture(labels[idx]));
      scene.materials.push_back(
          make_text_material(object.text, (int)scene.textures.size() - 1));
      scene.shapes.push_back(make_text_shape(
          labels[idx], offsets[idx], texture, object.textscale));
      auto& ninstance    = scene.instances.emplace_back();
      ninstance.frame    = make_text_frame(labels[idx], texture,
             scene.cameras.front().frame.o,
             transform_point(object.frame, positions[idx]), object.textscale);
      ninstance.material = (int)scene.materials.size() - 1;
      ninstance.shape    = (int)scene.shapes.size() - 1;
    }
  }
}

// Convert diagram to scene
static scene_data convert_scene(const diagram_scene& diagram) {
  auto scene = scene_data{};
  scene.cameras.push_back(make_default_camera());
  for (auto& object : diagram.objects) add_object(scene, object);
  return scene;
}

// crop image vertically
static array2d<vec4f> crop_image(const array2d<vec4f>& image) {
  // find min and max
  auto [width, height] = image.extents();
  auto min = 0, max = height;
  for (auto j : range(height)) {
    auto empty = true;
    for (auto i : range(width)) {
      if (image[{i, j}] != vec4f{1, 1, 1, 1}) empty = false;
    }
    if (empty) {
      min = j;
    } else {
      break;
    }
  }
  for (auto j : range(height)) {
    auto empty = true;
    for (auto i : range(width)) {
      if (image[{i, height - j - 1}] != vec4f{1, 1, 1, 1}) empty = false;
    }
    if (empty) {
      max = height - j;
    } else {
      break;
    }
  }

  // crop
  auto cropped = array2d<vec4f>({width, max - min});
  for (auto j : range(min, max)) {
    for (auto i : range(width)) {
      cropped[{i, j - min}] = image[{i, j}];
    }
  }

  // done
  return cropped;
}

// Image rendering
static array2d<vec4f> render_image(const scene_data& scene, int resolution,
    int samples, bool noparallel = false);

// Render a diagram to an image
array2d<vec4f> render_diagram(
    const diagram_data& diagram, bool draft, bool crop) {
  // final image
  auto composite = array2d<vec4f>{{1440, 960}, {1, 1, 1, 1}};

  // render scenes
  for (auto& scene : diagram.scenes) {
    // convert diagram to scene
    auto yscene = convert_scene(scene);

    // image
    auto image = render_image(yscene, 1440, draft ? 16 : 256, false);

    // helpers
    auto draw_quad = [](array2d<vec4f>& image, const vec2i& center,
                         const vec2i& size, const vec4f& color) {
      auto extents = (vec2i)image.extents();
      for (auto i : range(-size.x / 2, +size.x / 2)) {
        image[clamp(center + vec2i{i, -size.y / 2}, vec2i{0, 0}, extents - 1)] =
            color;
        image[clamp(center + vec2i{i, +size.y / 2}, vec2i{0, 0}, extents - 1)] =
            color;
      }
      for (auto j : range(-size.y / 2, +size.y / 2)) {
        image[clamp(center + vec2i{-size.x / 2, j}, vec2i{0, 0}, extents - 1)] =
            color;
        image[clamp(center + vec2i{+size.x / 2, j}, vec2i{0, 0}, extents - 1)] =
            color;
      }
    };
    auto copy_quad = [](array2d<vec4f>& image, const array2d<vec4f>& source,
                         const vec2i& icenter, const vec2i& scenter,
                         const vec2i& size, const vec4f& color) {
      auto iextents = (vec2i)image.extents();
      auto sextents = (vec2i)source.extents();
      for (auto j : range(-size.y / 2, +size.y / 2)) {
        for (auto i : range(-size.x / 2, +size.x / 2)) {
          auto ij = vec2i{i, j};
          image[clamp(ij + icenter, vec2i{0, 0}, iextents - 1)] =
              source[clamp(ij + scenter, vec2i{0, 0}, sextents - 1)];
        }
      }
    };

    // composite
    auto size    = (vec2i)(scene.size * 1440 / 10);
    auto margin  = (vec2i)(scene.margin * 1440 / 10);
    auto icenter = (vec2i)image.extents() / 2;
    auto ccenter = (vec2i)composite.extents() / 2 +
                   (vec2i)(scene.offset * 1440 / 10);

    copy_quad(composite, image, ccenter, icenter, size + margin, {1, 0, 0, 1});
    draw_quad(composite, ccenter, size, {0, 1, 0, 1});
    draw_quad(composite, ccenter, size + margin, {0, 0, 1, 1});

    /*
        // composite
        auto crop_size        = (vec2i)((scene.size + scene.margin) * 1440 /
      10); auto crop_offset      = (vec2i)image.extents() / 2 - crop_size / 2;
        auto composite_offset = (vec2i)composite.extents() / 2 +
                                (vec2i)(scene.offset * 1440 / 10) - crop_size
      / 2; for (auto ij : range(crop_size)) { composite[clamp(ij +
      composite_offset, vec2s{0, 0}, composite.extents() - 1)] = image[ij +
      crop_offset]; if (ij.x == 0 || ij.y == 0 || ij.x == crop_size.x - 1 ||
              ij.y == crop_size.y - 1)
            composite[clamp(ij + composite_offset, vec2s{0, 0},
                composite.extents() - 1)] = {1, 0, 0, 1};
        }
      }
      */
  }

  // crop
  // if (crop) composite = crop_image_(composite);

  // done
  return composite;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide catmullclark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_impl_(
    const vector<vec4i>& quads, const vector<T>& vert, bool uncorrected) {
  // get edges
  auto emap     = make_edge_map(quads);
  auto edges    = get_edges(emap);
  auto boundary = get_boundary(emap);
  // number of elements
  auto nverts = (int)vert.size();
  auto nedges = (int)edges.size();
  auto nfaces = (int)quads.size();

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

  // averaging pass ----------------------------------
  auto avert  = vector<T>(tvert.size(), T{});
  auto acount = vector<int>(tvert.size(), 0);
  for (auto& q : tquads) {
    auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4;
    for (auto vid : {q.x, q.y, q.z, q.w}) {
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];

  // correction pass ----------------------------------
  // p = p + (avg_p - p) * (4/avg_count)
  if (!uncorrected) {
    for (auto i = 0; i < tvert.size(); i++) {
      avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4 / (float)acount[i]);
    }
  }
  tvert = avert;

  return {tquads, tvert};
}

// Bspline subdivision
pair<vector<vec4i>, vector<float>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<float>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec2f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec3f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec4f>& vertices,
    bool uncorrected) {
  return subdivide_catmullclark_impl_(quads, vertices, uncorrected);
}

// Subdivide bspline.
template <typename T>
static pair<vector<vec2i>, vector<T>> subdivide_bspline_impl(
    const vector<vec2i>& lines, const vector<T>& vert) {
  // split elements ------------------------------------
  auto [tlines, tvert] = subdivide_lines(lines, vert);

  // boundary ------------------------------------------
  auto valence = vector<bool>(tvert.size(), 1);
  // TODO: fixme
  auto tcreases = vector<int>{};

  // averaging pass ----------------------------------
  auto avert  = vector<T>(tvert.size(), T{});
  auto acount = vector<int>(tvert.size(), 0);
  for (auto& p : tcreases) {
    auto c = tvert[p];
    for (auto vid : {p}) {
      if (valence[vid] != 0) continue;
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto& l : tlines) {
    auto c = (tvert[l.x] + tvert[l.y]) / 2;
    for (auto vid : {l.x, l.y}) {
      if (valence[vid] != 1) continue;
      avert[vid] += c;
      acount[vid] += 1;
    }
  }
  for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];
  tvert = avert;

  return {tlines, tvert};
}

pair<vector<vec2i>, vector<float>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<float>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec2f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec2f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec3f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec3f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}
pair<vector<vec2i>, vector<vec4f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec4f>& vertices) {
  return subdivide_bspline_impl(lines, vertices);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM RENDERER
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the two integer indices.
template <typename Func>
inline void parallel_for_batch(vec2i num, Func&& func) {
  auto              futures  = vector<std::future<void>>{};
  auto              nthreads = std::thread::hardware_concurrency();
  std::atomic<int>  next_idx(0);
  std::atomic<bool> has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, &has_error, num]() {
          try {
            while (true) {
              auto j = next_idx.fetch_add(1);
              if (j >= num[1]) break;
              if (has_error) break;
              for (auto i = 0; i < num[0]; i++) func(vec2i{i, j});
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Convenience functions
[[maybe_unused]] static vec3f eval_position(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec2f eval_texcoord(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_texcoord(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static material_point eval_material(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_material(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}

// Generates a ray from a camera.
static ray3f eval_camera(const camera_data& camera, const vec2f& uv_) {
  auto film = vec2f{camera.film, camera.film / camera.aspect};
  if (!camera.orthographic) {
    // point on film
    auto uv = flip_u(uv_);
    auto q  = vec3f{film * (uv - 0.5f), camera.lens};
    // point on the lens
    auto e = vec3f{0, 0, 0};
    // ray direction through the lens center
    auto d = normalize(e - q);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  } else {
    // point on film
    auto uv    = flip_u(uv_);
    auto scale = 1 / camera.lens;
    auto q     = vec3f{film * (uv - 0.5f) * scale, camera.lens};
    // point on the lens
    auto e = vec3f{-xy(q), 0};
    // correct ray direction to account for camera focusing
    auto d = normalize(e - q);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  }
}

// Diagram previewing.
static vec4f render_ray(
    const scene_data& scene, const scene_bvh& bvh, const ray3f& ray_) {
  // initialize
  auto radiance = vec3f{0, 0, 0};
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;

  // trace  path
  for (auto bounce = 0; bounce < 16; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      radiance += weight * vec3f{1, 1, 1};
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_position(scene, intersection);
    auto normal   = eval_normal(scene, intersection);
    auto material = eval_material(scene, intersection);

    // accumulate color
    radiance += weight * material.emission * material.opacity;

    // handle opacity
    if (material.opacity >= 1) break;
    weight *= 1 - material.opacity;
    ray = {position + ray.d * 1e-2f, ray.d};
  }

  return {radiance, 1.0f};
}

// Progressively computes an image.
static array2d<vec4f> render_image(
    const scene_data& scene, int resolution, int samples, bool noparallel) {
  // Bvh
  auto bvh = make_scene_bvh(scene, false, noparallel);
  // Camera
  auto& camera = scene.cameras[0];
  // Image
  auto size  = vec2i{resolution, (int)round(resolution / camera.aspect)};
  auto image = array2d<vec4f>{size, {0, 0, 0, 0}};
  // Start rendering
  auto nsamples = (int)round(sqrt((float)samples));
  if (noparallel) {
    for (auto ij : range(size)) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto ray = eval_camera(
            camera, (ij + (sij + 0.5f) / (float)nsamples) / image.extents());
        auto color = render_ray(scene, bvh, ray);
        image[ij] += isfinite(color) ? color : vec4f{1, 1, 1, 1};
      }
      image[ij] /= nsamples * nsamples;
    }
  } else {
    parallel_for_batch(size, [&](vec2i ij) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto ray = eval_camera(
            camera, (ij + (sij + 0.5f) / (float)nsamples) / image.extents());
        auto color = render_ray(scene, bvh, ray);
        image[ij] += isfinite(color) ? color : vec4f{1, 1, 1, 1};
      }
      image[ij] /= nsamples * nsamples;
    });
  }

  return image;
}

}  // namespace yocto
