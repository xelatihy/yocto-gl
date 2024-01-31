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

#include <future>
#include <unordered_set>

#include "yocto_bvh.h"
#include "yocto_sampling.h"
#include "yocto_sceneio.h"

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
// SHAPE HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Merge shape elements
static void merge_dshape_inplace(shape_data& shape, const shape_data& merge) {
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

// Make a tesselated rectangle. Useful in other subdivisions.
static shape_data make_dquads(vec2i steps, vec2f scale, vec2f uvscale) {
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

// Make a rect
static shape_data make_drect(vec2i steps, vec2f scale, vec2f uvscale = {1, 1}) {
  return make_dquads(steps, scale, uvscale);
}

// Make a box.
static shape_data make_dbox(
    vec3i steps, vec3f scale, vec3f uvscale = {1, 1, 1}) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // + z
  qshape = make_drect(
      {steps.x, steps.y}, {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qshape.positions) p = {p.x, p.y, scale.z};
  for (auto& n : qshape.normals) n = {0, 0, 1};
  merge_dshape_inplace(shape, qshape);
  // - z
  qshape = make_drect(
      {steps.x, steps.y}, {scale.x, scale.y}, {uvscale.x, uvscale.y});
  for (auto& p : qshape.positions) p = {-p.x, p.y, -scale.z};
  for (auto& n : qshape.normals) n = {0, 0, -1};
  merge_dshape_inplace(shape, qshape);
  // + x
  qshape = make_drect(
      {steps.z, steps.y}, {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qshape.positions) p = {scale.x, p.y, -p.x};
  for (auto& n : qshape.normals) n = {1, 0, 0};
  merge_dshape_inplace(shape, qshape);
  // - x
  qshape = make_drect(
      {steps.z, steps.y}, {scale.z, scale.y}, {uvscale.z, uvscale.y});
  for (auto& p : qshape.positions) p = {-scale.x, p.y, p.x};
  for (auto& n : qshape.normals) n = {-1, 0, 0};
  merge_dshape_inplace(shape, qshape);
  // + y
  qshape = make_drect(
      {steps.x, steps.z}, {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i : range(qshape.positions.size())) {
    qshape.positions[i] = {
        qshape.positions[i].x, scale.y, -qshape.positions[i].y};
    qshape.normals[i] = {0, 1, 0};
  }
  merge_dshape_inplace(shape, qshape);
  // - y
  qshape = make_rect(
      {steps.x, steps.z}, {scale.x, scale.z}, {uvscale.x, uvscale.z});
  for (auto i : range(qshape.positions.size())) {
    qshape.positions[i] = {
        qshape.positions[i].x, -scale.y, qshape.positions[i].y};
    qshape.normals[i] = {0, -1, 0};
  }
  merge_dshape_inplace(shape, qshape);
  return shape;
}

// Make a sphere.
static shape_data make_dsphere(int steps, float scale, float uvscale = 1) {
  auto shape = make_dbox({steps, steps, steps}, {scale, scale, scale},
      {uvscale, uvscale, uvscale});
  for (auto& p : shape.positions) p = normalize(p) * scale;
  shape.normals = shape.positions;
  for (auto& n : shape.normals) n = normalize(n);
  return shape;
}

// Make a sphere.
static shape_data make_duvsphere(
    vec2i steps, float scale, vec2f uvscale = {1, 1}) {
  auto shape = make_drect({steps.x, steps.y}, {1, 1});
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

// Make a hemisphere.
static shape_data make_duvhemisphere(
    vec2i steps, float scale, vec2f uvscale = {1, 1}) {
  auto shape = make_drect({steps.x, steps.y}, {1, 1});
  for (auto i : range(shape.positions.size())) {
    auto uv = shape.texcoords[i];
    auto a  = vec2f{pif * uv.x, pif * uv.y};
    shape.positions[i] =
        vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)} * scale;
    shape.normals[i]   = normalize(shape.positions[i]);
    shape.texcoords[i] = uv * uvscale;
  }
  return shape;
}

// Make a uv cylinder
shape_data make_duvcylinder(
    vec3i steps, vec2f scale, vec3f uvscale = {1, 1, 1}) {
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

// Make a uv capsule
static shape_data make_duvcapsule(
    vec3i steps, vec2f scale, vec3f uvscale = {1, 1, 1}) {
  auto shape  = shape_data{};
  auto qshape = shape_data{};
  // side
  qshape = make_drect({steps.x, steps.y}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv             = qshape.texcoords[i];
    auto phi            = 2 * pif * uv.x;
    qshape.positions[i] = {
        cos(phi) * scale.x, sin(phi) * scale.x, (2 * uv.y - 1) * scale.y};
    qshape.normals[i]   = {cos(phi), sin(phi), 0};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.y};
  }
  for (auto& quad : qshape.quads) quad = {quad.x, quad.w, quad.z, quad.y};
  merge_dshape_inplace(shape, qshape);
  // top
  qshape = make_drect({steps.x, steps.z}, {1, 1}, {1, 1});
  for (auto i : range(qshape.positions.size())) {
    auto uv  = qshape.texcoords[i];
    auto phi = uv.x * 2 * pif, theta = (1 - uv.y) * pif / 2;
    qshape.positions[i] = {cos(phi) * sin(theta) * scale.x,
        sin(phi) * sin(theta) * scale.x, cos(theta) * scale.x + scale.y};
    qshape.normals[i]   = {
        cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
    qshape.texcoords[i] = uv * vec2f{uvscale.x, uvscale.z};
  }
  merge_dshape_inplace(shape, qshape);
  // bottom
  qshape = make_drect({steps.x, steps.z}, {1, 1}, {1, 1});
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
  merge_dshape_inplace(shape, qshape);
  return shape;
}

// Make a uv cone
static shape_data make_duvcone(
    vec3i steps, vec2f scale, vec3f uvscale = {1, 1, 1}) {
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
// TODO: remove
template <typename Params>
static frame3f params_frame(const Params& params) {
  return params.frame;
}

// Helper
template <typename T>
struct zip_labels {
  struct iterator {
    iterator(size_t idx, const vector<string>& labels, const vector<T>& others)
        : idx{idx}, labels{labels}, others{others} {}

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

  zip_labels(const vector<string>& labels, const vector<T>& others)
      : labels{labels}, others{others} {}
  iterator begin() { return {(size_t)0, labels, others}; }
  iterator end() { return {min(labels.size(), others.size()), labels, others}; }

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
diagram_scene& add_scene(diagram_data& diagram, const string& title, vec2f size,
    const frame3f& frame, vec2f margin) {
  return add_scenews(diagram, title, "", size, frame, margin);
}
diagram_scene& add_scene(diagram_data& diagram, const string& title,
    const frame3f& frame, vec2f margin) {
  return add_scenews(diagram, title, "", {2, 2}, frame, margin);
}
diagram_scene& add_scenews(diagram_data& diagram, const string& title,
    const string& subtitle, vec2f size, const frame3f& frame, vec2f margin) {
  auto& scene  = diagram.scenes.emplace_back();
  scene.size   = size;
  scene.margin = margin;
  // scene.frame  = dscaling({params.scale, params.scale, params.scale}) *
  //               dtranslation({-params.center.x, -params.center.y, 0}) *
  //               params.frame;
  scene.frame = frame;

  if (!title.empty()) {
    scene.objects.push_back({.labels = {
                                 .labels    = {title + "!!t"},
                                 .positions = {{0, scene.size.y / 2, 0}},
                             }});
  }

  if (!subtitle.empty()) {
    scene.objects.push_back({.labels = {
                                 .labels    = {subtitle + "!!b"},
                                 .positions = {{0, -scene.size.y / 2, 0}},
                             }});
  }

  update_offsets(diagram);

  return scene;
}

// Labels
diagram_labels dlabels(const vector<string>& labels) {
  return dlabels({}, labels);
}
diagram_labels dlabels(
    const vector<vec3f>& positions, const vector<string>& labels) {
  return {.labels = labels, .positions = positions};
}
diagram_labels dllabels(const vector<string>& labels) {
  return {.llabels = labels};
}
diagram_labels dflabels(const vector<string>& labels) {
  return {.flabels = labels};
}
diagram_labels dclabels(const vector<string>& labels) {
  return {.clabels = labels};
}
diagram_labels dglabels(const vector<string>& labels,
    const vector<string>& llabels, const vector<string>& flabels,
    const vector<string>& clabels) {
  return {.labels = labels,
      .llabels    = llabels,
      .flabels    = flabels,
      .clabels    = clabels,
      .positions  = {}};
}
diagram_labels dglabels(const vector<vec3f>& positions,
    const vector<string>& labels, const vector<string>& llabels,
    const vector<string>& flabels, const vector<string>& clabels) {
  return {.labels = labels,
      .llabels    = llabels,
      .flabels    = flabels,
      .clabels    = clabels,
      .positions  = positions};
}

// Add labels
void add_labels(diagram_scene& diagram, const diagram_labels& labels,
    const diagram_style& style) {
  add_labels(diagram, identity3x4f, labels, style);
}

void add_labels(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_labels& labels, const diagram_style& style) {
  diagram.objects.push_back({diagram.frame * frame.frame, {}, labels, style});
}

// Add shape
void add_shape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_style& style) {
  return add_shape(diagram, identity3x4f, shape, diagram_labels{}, style);
}
void add_shape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_labels& labels, const diagram_style& style) {
  return add_shape(diagram, identity3x4f, shape, labels, style);
}
void add_shape(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_shape& shape, const diagram_style& style) {
  return add_shape(diagram, frame.frame, shape, diagram_labels{}, style);
}

void add_shape(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_shape& shape_, const diagram_labels& labels_,
    const diagram_style& style) {
  // Copy back
  auto shape = shape_;

  // TODO: wireframe
  auto wireframe = true;
  if (wireframe && shape.lines.empty() && !shape.triangles.empty()) {
    shape.lines = triangles_edges(shape.triangles);
  }
  if (wireframe && shape.lines.empty() && !shape.quads.empty()) {
    shape.lines = quads_edges(shape.quads);
  }

  auto labels = diagram_labels{};
  for (auto&& [label, position] : zip_labels(labels_.labels, shape.positions)) {
    labels.labels.push_back(label);
    labels.positions.push_back(position);
  }
  for (auto&& [label, line] : zip_labels(labels_.llabels, shape.lines)) {
    labels.labels.push_back(label);
    labels.positions.push_back(
        (shape.positions[line.x] + shape.positions[line.y]) / 2);
  }
  for (auto&& [label, arrow] : zip_labels(labels_.llabels, shape.arrows)) {
    labels.labels.push_back(label);
    labels.positions.push_back(
        (shape.positions[arrow.x] + shape.positions[arrow.y]) / 2);
  }
  for (auto&& [label, triangle] :
      zip_labels(labels_.flabels, shape.triangles)) {
    labels.labels.push_back(label);
    labels.positions.push_back(
        (shape.positions[triangle.x] + shape.positions[triangle.y] +
            shape.positions[triangle.z]) /
        3);
  }
  for (auto&& [label, quad] : zip_labels(labels_.flabels, shape.quads)) {
    labels.labels.push_back(label);
    labels.positions.push_back(
        (shape.positions[quad.x] + shape.positions[quad.y] +
            shape.positions[quad.z] + shape.positions[quad.w]) /
        4);
  }
  if (!labels_.clabels.empty()) {
    auto bbox = invalidb3f;
    for (auto& position : shape.positions) bbox = merge(bbox, position);
    auto center = (bbox.min + bbox.max) / 2;
    for (auto&& label : labels_.clabels) {
      labels.labels.push_back(label);
      labels.positions.push_back(center);
    }
  }

  diagram.objects.push_back(
      {diagram.frame * frame.frame, shape, labels, style});
}

// Add points
diagram_shape dpoints(const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  for (auto idx : range((int)positions.size())) shape.points.push_back(idx);
  return shape;
}
diagram_shape dpoints(
    const vector<int>& points, const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.points    = points;
  return shape;
}

// Add lines
diagram_shape dlines(const vector<vec3f>& positions) {
  auto lines = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  return dlines(lines, positions);
}
diagram_shape dlines(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.lines     = lines;
  return shape;
}
diagram_shape dclines(const vector<vec3f>& positions) {
  auto lines = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    lines.push_back({idx * 2 + 0, idx * 2 + 1});
  return dclines(lines, positions);
}
diagram_shape dclines(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.lines     = lines;
  shape.connect   = 0.3f;
  return shape;
}

// Add arrows
diagram_shape darrows(const vector<vec3f>& positions) {
  auto arrows = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    arrows.push_back({idx * 2 + 0, idx * 2 + 1});
  return darrows(arrows, positions);
}
diagram_shape darrows(
    const vector<vec2i>& arrows, const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.arrows    = arrows;
  return shape;
}
diagram_shape dcarrows(const vector<vec3f>& positions) {
  auto arrows = vector<vec2i>{};
  for (auto idx : range((int)positions.size() / 2))
    arrows.push_back({idx * 2 + 0, idx * 2 + 1});
  return dcarrows(arrows, positions);
}
diagram_shape dcarrows(
    const vector<vec2i>& arrows, const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.arrows    = arrows;
  shape.connect   = 0.3f;
  return shape;
}

// Add vectors
diagram_shape dvectors(const vector<vec3f>& vectors) {
  auto shape      = diagram_shape{};
  shape.positions = vectors;
  shape.positions.push_back({0, 0, 0});
  for (auto idx : range((int)vectors.size()))
    shape.arrows.push_back({(int)vectors.size(), idx});
  return shape;
}
diagram_shape dcvectors(const vector<vec3f>& vectors) {
  auto shape      = diagram_shape{};
  shape.positions = vectors;
  shape.positions.push_back({0, 0, 0});
  for (auto idx : range((int)vectors.size()))
    shape.arrows.push_back({(int)vectors.size(), idx});
  shape.connect = 0.3;
  return shape;
}

// Add axes
diagram_shape daxes(const frame3f& axes, vec3f aspect) {
  auto shape      = diagram_shape{};
  shape.positions = {axes.o, axes.o + aspect.x * axes.x,
      axes.o + aspect.y * axes.y, axes.o + aspect.z * axes.z};
  shape.arrows    = {{0, 1}, {0, 2}, {0, 3}};
  return shape;
}
diagram_shape daxes2(const frame3f& axes, vec3f aspect) {
  auto shape      = diagram_shape{};
  shape.positions = {
      axes.o, axes.o + aspect.x * axes.x, axes.o + aspect.y * axes.y};
  shape.arrows = {{0, 1}, {0, 2}};
  return shape;
}

// Add rays
diagram_shape drays(const vector<ray3f>& rays) {
  auto positions = vector<vec3f>{};
  for (auto& ray : rays) {
    positions.push_back(ray.o);
    positions.push_back(ray.o + ray.d);
  }
  return darrows(positions);
}
diagram_shape drayscont(const vector<ray3f>& rays, float length) {
  auto positions = vector<vec3f>{};
  for (auto& ray : rays) {
    positions.push_back(ray.o + ray.d);
    positions.push_back(ray.o + ray.d * length);
  }
  return dclines(positions);
}

// Add quads
diagram_shape dquads(
    const vector<vec3f>& positions, const vector<vec2f>& texcoords) {
  auto quads = vector<vec4i>{};
  for (auto idx : range((int)positions.size() / 4))
    quads.push_back({idx * 4 + 0, idx * 4 + 1, idx * 4 + 2, idx * 4 + 3});
  return dquads(quads, positions, texcoords);
}
diagram_shape dquads(const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec2f>& texcoords) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.texcoords = texcoords;
  if (texcoords.empty() && positions.size() == 4) {
    shape.texcoords = dconstants::quad_texcoords;
  }

  if (quads.empty()) {
    for (auto idx : range((int)shape.positions.size() / 4))
      shape.quads.push_back(
          {idx * 4 + 0, idx * 4 + 1, idx * 4 + 2, idx * 4 + 3});
  } else {
    shape.quads = quads;
  }
  shape.lines = quads_edges(shape.quads);

  return shape;
}

// Add triangles
diagram_shape dtriangles(
    const vector<vec3f>& positions, const vector<vec2f>& texcoords) {
  auto triangles = vector<vec3i>{};
  for (auto idx : range((int)positions.size() / 3))
    triangles.push_back({idx * 3 + 0, idx * 3 + 1, idx * 3 + 2});
  return dtriangles(triangles, positions, {});
}
diagram_shape dtriangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec2f>& texcoords) {
  auto shape      = diagram_shape{};
  shape.positions = positions;
  shape.texcoords = texcoords;

  if (triangles.empty()) {
    for (auto idx : range((int)positions.size() / 3))
      shape.triangles.push_back({idx * 3 + 0, idx * 3 + 1, idx * 3 + 2});
  } else {
    shape.triangles = triangles;
  }
  shape.lines = triangles_edges(shape.triangles);

  return shape;
}

// Add shape
diagram_shape dshape(const shape_data& shape_, bool wireframe) {
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.points    = shape_.points;
  shape.lines     = shape_.lines;
  shape.triangles = shape_.triangles;
  shape.quads     = shape_.quads;

  if (wireframe && shape.lines.empty() && !shape.triangles.empty()) {
    shape.lines = triangles_edges(shape.triangles);
  }
  if (wireframe && shape.lines.empty() && !shape.quads.empty()) {
    shape.lines = quads_edges(shape.quads);
  }

  return shape;
}

// Add polyline
diagram_shape dpolyline(const vector<vec3f>& positions) {
  auto shape      = diagram_shape{};
  shape.positions = positions;

  for (auto idx : range((int)positions.size() - 1))
    shape.lines.push_back({idx + 0, idx + 1});

  return shape;
}

// Add polygon
diagram_shape dpolygon(const vector<vec3f>& positions) {
  auto bbox = invalidb3f;
  for (auto& position : positions) bbox = merge(bbox, position);
  auto center = (bbox.min + bbox.max) / 2;

  auto steps      = (int)positions.size();
  auto positions_ = positions;
  positions_.push_back(center);

  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(steps)) {
    shape.triangles.push_back({idx, (idx + 1) % steps, steps});
  }

  for (auto idx : range(steps)) {
    shape.lines.push_back({idx, (idx + 1) % steps});
  }

  return shape;
}

// Add rect
diagram_shape drect(vec2f aspect) {
  auto positions_ = dconstants::quad_positions;
  auto scale      = vec3f{aspect.x / aspect.y, 1, 1};
  for (auto& position : positions_) position *= scale;
  auto shape      = diagram_shape{};
  shape.positions = positions_;
  shape.texcoords = dconstants::quad_texcoords;

  shape.quads = dconstants::quad_quads;
  shape.lines = quads_edges(shape.quads);

  return shape;
}

// Add disk
diagram_shape ddisk(int steps) {
  auto positions_ = vector<vec3f>{};
  for (auto idx : range(steps)) {
    auto theta = 2 * pif * idx / (float)steps;
    positions_.push_back({cos(theta), sin(theta), 0});
  }
  positions_.push_back({0, 0, 0});
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(steps)) {
    shape.triangles.push_back({idx, (idx + 1) % steps, steps});
  }

  for (auto idx : range(steps)) {
    shape.lines.push_back({idx, (idx + 1) % steps});
  }

  return shape;
}

// Add half disk
diagram_shape dhalfdisk(int steps) {
  auto positions_ = vector<vec3f>{};
  for (auto idx : range(steps + 1)) {
    auto theta = pif * idx / (float)steps;
    positions_.push_back({cos(theta), sin(theta), 0});
  }
  positions_.push_back({0, 0, 0});
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(steps)) {
    shape.triangles.push_back({idx, idx, steps + 1});
  }

  for (auto idx : range(steps)) {
    shape.lines.push_back({idx, idx + 1});
  }
  shape.lines.push_back({0, steps + 1});
  shape.lines.push_back({steps, steps + 1});

  return shape;
}

// Add arc
diagram_shape darc(vec3f from, vec3f to, vec3f center, int steps) {
  // refence frame
  auto x     = normalize(from - center);
  auto y     = normalize(to - center);
  auto angle = acos(dot(x, y));
  auto z     = normalize(cross(x, y));
  y          = normalize(cross(z, x));
  auto frame = frame3f{x, y, z, center};

  auto positions_ = vector<vec3f>{};
  for (auto idx : range(steps + 1)) {
    auto theta = angle * idx / (float)steps;
    positions_.push_back(transform_point(frame, {cos(theta), sin(theta), 0}));
  }
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(steps)) {
    shape.lines.push_back({idx, idx + 1});
  }

  return shape;
}

// Add qbezier
diagram_shape dqbeziers(const vector<vec3f>& positions, int steps) {
  auto qbezier_point = [](vec3f p1, vec3f p2, vec3f p3, float u) {
    return (1 - u) * (1 - u) * p1 + 2 * (1 - u) * u * p2 + u * u * p3;
  };
  auto positions_ = vector<vec3f>();
  for (auto segment : range(positions.size() / 3)) {
    for (auto idx : range(steps + 1)) {
      positions_.push_back(
          qbezier_point(positions[0 + segment * 3], positions[1 + segment * 3],
              positions[2 + segment * 3], idx / (float)steps));
    }
  }

  return dpolyline(positions_);
}

// Add bezier
diagram_shape dbeziers(const vector<vec3f>& positions, int steps) {
  auto bezier_point = [](vec3f p1, vec3f p2, vec3f p3, vec3f p4, float u) {
    return (1 - u) * (1 - u) * (1 - u) * p1 + 3 * (1 - u) * (1 - u) * u * p2 +
           3 * (1 - u) * u * u * p3 + u * u * u * p4;
  };
  auto positions_ = vector<vec3f>();
  for (auto segment : range(positions.size() / 4)) {
    for (auto idx : range(steps + 1)) {
      positions_.push_back(bezier_point(positions[0 + segment * 4],
          positions[1 + segment * 4], positions[2 + segment * 4],
          positions[3 + segment * 4], idx / (float)steps));
    }
  }

  return dpolyline(positions_);
}

// Add cube
diagram_shape dcube() {
  auto shape      = diagram_shape{};
  shape.positions = dconstants::cube_positions;

  shape.quads = dconstants::cube_quads;

  shape.lines = quads_edges(shape.quads);

  return shape;
}
diagram_shape dcube(int steps) {
  auto shape_     = make_dbox({steps, steps, steps}, {1, 1, 1});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;
  shape.quads     = shape_.quads;
  shape.lines     = quads_edges(shape.quads);
  return shape;
}

// Sphere
diagram_shape dsphere(int steps, bool isolines) {
  auto shape_     = make_dsphere(steps, 1);
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;
  if (isolines) {
    using dconstants::xup3x4f;
    using dconstants::yup3x4f;
    for (auto frame : {identity3x4f, xup3x4f, yup3x4f}) {
      auto offset = (int)shape.positions.size();
      for (auto idx = 0; idx < steps * 2; idx++) {
        auto u = 2 * pif * idx / float(steps * 2);
        shape.positions.push_back(
            transform_point(frame, vec3f{cos(u), sin(u), 0}));
      }
      for (auto idx = 0; idx < steps * 2; idx++) {
        shape.lines.push_back({offset + idx, offset + (idx + 1) % (steps * 2)});
      }
    }
  } else {
    shape.lines = quads_edges(shape.quads);
  }

  return shape;
}
diagram_shape duvsphere(int steps, bool isolines) {
  auto shape_     = make_duvsphere({steps * 2, steps}, 1);
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;
  if (isolines) {
    using dconstants::xup3x4f;
    using dconstants::yup3x4f;
    for (auto frame : {identity3x4f, xup3x4f, yup3x4f}) {
      auto offset = (int)shape.positions.size();
      for (auto idx = 0; idx < steps * 2; idx++) {
        auto u = 2 * pif * idx / float(steps * 2);
        shape.positions.push_back(
            transform_point(frame, vec3f{cos(u), sin(u), 0}));
      }
      for (auto idx = 0; idx < steps * 2; idx++) {
        shape.lines.push_back({offset + idx, offset + (idx + 1) % (steps * 2)});
      }
    }
  } else {
    shape.lines = quads_edges(shape.quads);
  }

  return shape;
}

// Add hemisphere
diagram_shape dhemisphere(int steps, bool isolines) {
  auto shape_     = make_duvhemisphere({steps * 2, steps / 2}, 1);
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;
  if (isolines) {
    using dconstants::xup3x4f;
    using dconstants::yup3x4f;
    auto offset = (int)shape.positions.size();
    for (auto idx = 0; idx < steps * 2; idx++) {
      auto u = 2 * pif * idx / float(steps * 2);
      shape.positions.push_back({cos(u), sin(u), 0});
      shape.lines.push_back({offset + idx, offset + (idx + 1) % (steps * 2)});
    }
    for (auto frame : {yup3x4f, drotation({1, 0, 0}, 90) * xup3x4f}) {
      auto offset = (int)shape.positions.size();
      for (auto idx = 0; idx < steps + 1; idx++) {
        auto u = pif * idx / float(steps);
        shape.positions.push_back(transform_point(frame, {cos(u), sin(u), 0}));
        if (idx != steps)
          shape.lines.push_back({offset + idx, offset + (idx + 1)});
      }
    }
  } else {
    shape.lines = quads_edges(shape.quads);
  }

  return shape;
}

// Add cone
diagram_shape dcone(int steps, bool isolines) {
  auto shape_     = make_duvcone({steps, 1, 1}, {1, 1});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;
  if (isolines) {
    auto offset = (int)shape.positions.size();
    for (auto idx = 0; idx < steps * 2; idx++) {
      auto u = 2 * pif * idx / float(steps * 2);
      shape.positions.push_back({cos(u), sin(u), 0});
    }
    for (auto idx = 0; idx < steps * 2; idx++) {
      shape.lines.push_back({offset + idx, offset + (idx + 1) % steps});
    }
  } else {
    shape.lines = quads_edges(shape.quads);
  }
  return shape;
}

// Cylinder
diagram_shape duvcylinder(
    float hscale, int rsteps, int hsteps, int csteps, bool isolines) {
  auto shape_     = make_duvcylinder({rsteps, hsteps, csteps}, {1, hscale});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;

  // TODO: isolines

  return shape;
}

// Capsule
diagram_shape duvcapsule(
    float hscale, int rsteps, int hsteps, int csteps, bool isolines) {
  auto shape_     = make_duvcapsule({rsteps, hsteps, csteps}, {1, hscale});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.normals   = shape_.normals;
  shape.texcoords = shape_.texcoords;

  shape.quads = shape_.quads;

  // TODO: isolines

  return shape;
}

// Add bbox
diagram_shape dbbox(const bbox3f& bbox_, float epsilon) {
  auto bbox       = bbox3f{bbox_.min - epsilon, bbox_.max + epsilon};
  auto positions_ = dconstants::cube_positions;
  for (auto& position : positions_) {
    position.x = position.x < 0 ? bbox.min.x : bbox.max.x;
    position.y = position.y < 0 ? bbox.min.y : bbox.max.y;
    position.z = position.z < 0 ? bbox.min.z : bbox.max.z;
  }

  auto shape      = diagram_shape{};
  shape.positions = positions_;

  shape.quads = dconstants::cube_quads;
  shape.lines = quads_edges(shape.quads);

  return shape;
}

// Add grid
diagram_shape ddotgrid(vec2i steps) {
  auto ratio      = (float)steps.x / (float)steps.y;
  auto shape_     = make_drect({steps.x, steps.y}, {ratio, 1});
  auto positions_ = vector<vec3f>{};
  for (auto& quad : shape_.quads) {
    positions_.push_back(
        (shape_.positions[quad.x] + shape_.positions[quad.y] +
            shape_.positions[quad.z] + shape_.positions[quad.w]) /
        4);
  }
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range((int)steps.x * (int)steps.y))
    shape.points.push_back(idx);

  return shape;
}

// Add grid
diagram_shape dgrid(vec2i steps) {
  auto ratio      = (float)steps.x / (float)steps.y;
  auto shape_     = make_drect({steps.x, steps.y}, {ratio, 1});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;

  shape.quads = shape_.quads;
  shape.lines = quads_edges(shape.quads);

  return shape;
}

// Add grid
diagram_shape ddiskgrid(vec2i steps, int dsteps) {
  auto shape = diagram_shape{};
  for (auto idxr : range(steps.x)) {
    auto offset = (int)shape.positions.size();
    auto radius = (idxr + 1) / (float)steps.x;
    for (auto idx : range(dsteps * 2)) {
      auto u = 2 * pif * idx / float(dsteps * 2);
      shape.positions.push_back({radius * cos(u), radius * sin(u), 0});
      shape.lines.push_back({offset + idx, offset + (idx + 1) % (dsteps * 2)});
    }
  }

  for (auto idx : range(steps.y)) {
    auto offset = (int)shape.positions.size();
    auto theta  = 2 * pif * idx / (float)steps.y;
    shape.positions.push_back({cos(theta), sin(theta), 0});
    shape.positions.push_back({cos(theta + pif), sin(theta + pif), 0});
    shape.lines.push_back({offset + 0, offset + 1});
  }

  return shape;
}
diagram_shape dudiskgrid(vec2i steps, int dsteps) {
  auto shape = diagram_shape{};
  for (auto idxr : range(steps.x)) {
    auto offset = (int)shape.positions.size();
    auto radius = sqrt((idxr + 1) / (float)steps.x);
    for (auto idx : range(dsteps * 2)) {
      auto u = 2 * pif * idx / float(dsteps * 2);
      shape.positions.push_back({radius * cos(u), radius * sin(u), 0});
      shape.lines.push_back({offset + idx, offset + (idx + 1) % (dsteps * 2)});
    }
  }

  for (auto idx : range(steps.y)) {
    auto offset = (int)shape.positions.size();
    auto theta  = 2 * pif * idx / (float)steps.y;
    shape.positions.push_back({cos(theta), sin(theta), 0});
    shape.positions.push_back({cos(theta + pif), sin(theta + pif), 0});
    shape.lines.push_back({offset + 0, offset + 1});
  }

  return shape;
}

// Add affine grid
diagram_shape daffinegrid(vec3f axes_a, vec3f axes_b, vec2i steps) {
  auto ratio      = (float)steps.x / (float)steps.y;
  auto shape_     = make_drect({steps.x, steps.y}, {ratio, 1});
  auto positions_ = shape_.positions;
  for (auto& position : positions_)
    position = position.x * axes_a + position.y * axes_b;
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  shape.quads = shape_.quads;
  shape.lines = quads_edges(shape.quads);

  return shape;
}

// Add image
diagram_shape dimagerect(vec2i size) {
  auto ratio      = (float)size.x / (float)size.y;
  auto scale      = vec3f{ratio, 1, 1};
  auto positions_ = dconstants::quad_positions;
  for (auto& position : positions_) position *= scale;
  auto shape      = diagram_shape{};
  shape.positions = positions_;
  shape.texcoords = dconstants::quad_texcoords;

  shape.quads = dconstants::quad_quads;
  shape.lines = quads_edges(shape.quads);

  return shape;
}
diagram_shape dimagerect(const image<vec4f>& img) {
  return dimagerect(img.size());
}

// Add image grid
diagram_shape dimagegrid(vec2i size) {
  auto ratio      = (float)size.x / (float)size.y;
  auto shape_     = make_dquads({size.x, size.y}, {ratio, 1}, {1, 1});
  auto shape      = diagram_shape{};
  shape.positions = shape_.positions;
  shape.texcoords = shape_.texcoords;
  shape.quads     = shape_.quads;
  return shape;
}
diagram_shape dimagegrid(const image<vec4f>& img) {
  return dimagegrid(img.size());
}

// Add image label
diagram_shape dimagelabel(vec2i size, float scale) {
  auto ratio      = (float)size.x / (float)size.y;
  auto oscale     = vec3f{3 * scale, 1, 1} * 0.1f;
  auto ocenter    = vec3f{ratio, -1, 0} + vec3f{-0.40, +0.175, +0.01};
  auto opositions = dconstants::quad_positions;
  for (auto& oposition : opositions) oposition = oposition * oscale + ocenter;
  return {
      .positions = opositions,
      .texcoords = dconstants::quad_texcoords,
      .quads     = dconstants::quad_quads,
  };
}
diagram_shape dimagelabel(const image<vec4f>& image, float scale) {
  return dimagelabel(image.size(), scale);
}

// Add random points
diagram_shape drandompoints(int num, bool stratified) {
  return drandompoints(drandompoints_type::quad, num, stratified);
}
diagram_shape drandompoints(drandompoints_type type, int num, bool stratified) {
  auto rng        = make_rng(187291);
  auto positions_ = vector<vec3f>{};
  auto ssteps     = int(round(sqrt((float)num)));
  for (auto idx : range(num)) {
    auto uv = !stratified ? rand2f(rng)
                          : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                (idx / ssteps + rand1f(rng)) / ssteps};
    switch (type) {
      case drandompoints_type::quad: {
        positions_.push_back({2 * uv.x - 1, 2 * uv.y - 1, 0});
      } break;
      case drandompoints_type::disk: {
        positions_.push_back({sample_disk(uv).x, sample_disk(uv).y, 0});
      } break;
      case drandompoints_type::disknu: {
        auto uv_ = vec2f{uv.x, uv.y * uv.y};
        positions_.push_back({sample_disk(uv_).x, sample_disk(uv_).y, 0});
      } break;
      case drandompoints_type::triangle: {
        auto& triangle = dconstants::triangle_positions;
        positions_.push_back(interpolate_triangle(
            triangle[0], triangle[1], triangle[2], sample_triangle(uv)));
      } break;
      case drandompoints_type::hemi: {
        positions_.push_back(sample_hemisphere(uv));
      } break;
      case drandompoints_type::hemicos: {
        positions_.push_back(sample_hemisphere_cos(uv));
      } break;
      case drandompoints_type::hemicospower: {
        positions_.push_back(sample_hemisphere_cospower(64, uv));
      } break;
      case drandompoints_type::sphere: {
        positions_.push_back(sample_sphere(uv));
      } break;
      case drandompoints_type::linex: {
        positions_.push_back({2 * uv.x - 1, 0, 0});
      } break;
      case drandompoints_type::liney: {
        positions_.push_back({0, 2 * uv.y - 1, 0});
      } break;
    }
  }
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(num)) shape.points.push_back(idx);

  return shape;
}

// Add random points
diagram_shape drandompoints3(int num, bool stratified) {
  return drandompoints3(drandompoints3_type::cube, num, stratified);
}
diagram_shape drandompoints3(
    drandompoints3_type type, int num, bool stratified) {
  auto rng        = make_rng(187291);
  auto positions_ = vector<vec3f>{};
  if (!stratified) {
    for (auto idx : range(num)) {
      auto uv = rand3f(rng);
      switch (type) {
        case drandompoints3_type::cube: {
          positions_.push_back(2 * uv - 1);
        } break;
        case drandompoints3_type::linex: {
          positions_.push_back({2 * uv.x - 1, 0, 0});
        } break;
        case drandompoints3_type::liney: {
          positions_.push_back({0, 2 * uv.y - 1, 0});
        } break;
        case drandompoints3_type::linez: {
          positions_.push_back({0, 0, 2 * uv.z - 1});
        } break;
      }
    }
  } else {
    auto ssteps = int(round(pow((float)num, 1 / 3.0f)));
    for (auto k : range(ssteps)) {
      for (auto j : range(ssteps)) {
        for (auto i : range(ssteps)) {
          auto idx = vec3i{i, j, k};
          auto uv  = ((vec3f)idx + rand3f(rng)) / ssteps;
          switch (type) {
            case drandompoints3_type::cube: {
              positions_.push_back(2 * uv - 1);
            } break;
            case drandompoints3_type::linex: {
              positions_.push_back({2 * uv.x - 1, 0, 0});
            } break;
            case drandompoints3_type::liney: {
              positions_.push_back({0, 2 * uv.y - 1, 0});
            } break;
            case drandompoints3_type::linez: {
              positions_.push_back({0, 0, 2 * uv.z - 1});
            } break;
          }
        }
      }
    }
  }
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(num)) shape.points.push_back(idx);

  return shape;
}

// Add random lines
diagram_shape drandomlines(int num, bool stratified) {
  return drandomlines(drandomlines_type::hemi, num, stratified);
}
diagram_shape drandomlines(drandomlines_type type, int num, bool stratified) {
  auto rng        = make_rng(187291);
  auto ssteps     = int(round(sqrt((float)num)));
  auto positions_ = vector<vec3f>{};
  for (auto idx : range(num)) {
    auto uv = !stratified ? rand2f(rng)
                          : vec2f{(idx % ssteps + rand1f(rng)) / ssteps,
                                (idx / ssteps + rand1f(rng)) / ssteps};
    switch (type) {
      case drandomlines_type::hemi:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere({0, 0, 1}, uv));
        break;
      case drandomlines_type::hemicos:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere_cos({0, 0, 1}, uv));
        break;
      case drandomlines_type::hemicospower:
        positions_.push_back({0, 0, 0});
        positions_.push_back(sample_hemisphere_cospower(64, {0, 0, 1}, uv));
        break;
      case drandomlines_type::beam:
        auto rdisk = sample_disk(uv);
        positions_.push_back({rdisk.x, rdisk.y, 0});
        positions_.push_back({rdisk.x, rdisk.y, 1});
        break;
    }
  }
  auto shape      = diagram_shape{};
  shape.positions = positions_;

  for (auto idx : range(num)) {
    shape.lines.push_back({idx * 2 + 0, idx * 2 + 1});
  }

  return shape;
}

// Add plot
diagram_scene& add_plot(diagram_data& diagram, const string& title, vec2f size,
    const frame3f& frame, vec2f margin) {
  return add_scene(diagram, title, size, frame, margin);
}

// Add plot axes
diagram_scene& add_plotaxes(diagram_scene& diagram, const bbox2f& bounds,
    const vector<pair<float, string>>& xticks,
    const vector<pair<float, string>>& yticks, const diagram_style& style,
    const diagram_style& lstyle) {
  // set frame
  auto min = bounds.min, max = bounds.max;
  auto center = (max + min) / 2, diagonal = (max - min);
  auto pframe = dscaling({diagram.size.x / diagonal.x,
                    diagram.size.y / diagonal.y, 1.0f}) *
                dtranslation({-center.x, -center.y, 0.0f});
  diagram.frame = diagram.frame * pframe;

  // add frame
  add_shape(diagram,
      dquads({{min.x, min.y, 0}, {max.x, min.y, 0}, {max.x, max.y, 0},
          {min.x, max.y, 0}}),
      style);

  auto lpositions = vector<vec3f>{};
  for (auto& [x, _] : xticks) {
    lpositions.push_back({x, min.y, 0});
    lpositions.push_back({x, max.y, 0});
  }
  for (auto& [y, _] : yticks) {
    lpositions.push_back({min.x, y, 0});
    lpositions.push_back({max.x, y, 0});
  }
  add_shape(diagram, dlines(lpositions), lstyle);

  auto tpositions = vector<vec3f>{};
  auto tlabels    = vector<string>{};
  for (auto&& [x, label] : xticks) {
    tpositions.push_back({x, min.y, 0});
    tlabels.push_back(label + "!!b");
  }
  for (auto&& [y, label] : yticks) {
    tpositions.push_back({min.x, y, 0});
    tlabels.push_back(label + "!!l");
  }
  add_labels(diagram, {.positions = tpositions, .labels = tlabels});

  return diagram;
}

// Add plot shape
void add_plotshape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_style& style) {
  return add_plotshape(diagram, shape, {}, style);
}
void add_plotshape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_labels& labels, const diagram_style& style) {
  return add_shape(diagram, shape, labels, style);
}

// Plot functions
diagram_shape dplotline(const vector<vec2f>& points) {
  auto shape = diagram_shape{};
  for (auto& point : points) shape.positions.push_back({point.x, point.y, 0});
  for (auto idx : range((int)points.size() - 1))
    shape.lines.push_back({idx, idx + 1});
  return shape;
}
diagram_shape dplotpoints(const vector<vec2f>& points) {
  auto shape = diagram_shape{};
  for (auto& point : points) shape.positions.push_back({point.x, point.y, 0});
  for (auto idx : range((int)points.size())) shape.points.push_back(idx);
  return shape;
}
diagram_shape dplotarrows(const vector<vec2f>& points) {
  auto shape = diagram_shape{};
  for (auto& point : points) shape.positions.push_back({point.x, point.y, 0});
  for (auto idx : range((int)points.size() - 1))
    shape.arrows.push_back({idx, idx + 1});
  return shape;
}
vector<vec2f> dplotfunc(
    const function<float(float)>& func, vec2f range_, int samples) {
  auto points = vector<vec2f>(samples);
  for (auto idx : range(samples)) {
    auto x      = lerp(range_.x, range_.y, (float)idx / (float)(samples - 1));
    points[idx] = {x, func(x)};
  }
  return points;
}
vector<vec2f> dplotpfunc(const function<float(float)>& func, int samples) {
  auto points = vector<vec2f>(samples + 1);
  for (auto idx : range(samples)) {
    auto theta  = lerp(0, 2 * pif, (float)idx / (float)(samples));
    auto radius = func(theta);
    points[idx] = {radius * cos(theta), radius * sin(theta)};
  }
  points[samples] = points[0];
  return points;
}
vector<vec2f> dplotcurve(const vector<float>& curve, bool center) {
  auto size   = (int)curve.size();
  auto points = vector<vec2f>(size);
  for (auto i : range(size))
    points[i] = {center ? ((i + 0.5f) / size) : (float(i) / size), curve[i]};
  return points;
}
vector<vec2f> dplotcurve(const vector<float>& curve, vec2f range) {
  auto size   = (int)curve.size();
  auto points = vector<vec2f>(curve.size());
  for (auto i : yocto::range(size))
    points[i] = {lerp(range.x, range.y, i / (float)(size - 1)), curve[i]};
  return points;
}

// Add plot
diagram_scene& add_plot3(diagram_data& diagram, const string& title, vec2f size,
    const frame3f& frame, vec2f margin) {
  return add_scene(diagram, title, size, frame, margin);
}

// Add plot axes
diagram_scene& add_plotaxes3(diagram_scene& diagram, const bbox3f& bounds,
    vec3f size, const vector<pair<float, string>>& xticks,
    const vector<pair<float, string>>& yticks,
    const vector<pair<float, string>>& zticks, const diagram_style& style,
    const diagram_style& lstyle) {
  // set frame
  auto min = bounds.min, max = bounds.max;
  auto center = (max + min) / 2, diagonal = (max - min);
  auto pframe = dscaling(size / diagonal) * dtranslation(-center);
  // diagram.frame = diagram.frame * pframe * dconstants::yup3x4f;
  diagram.frame = diagram.frame * pframe;

  // add frame
  add_shape(diagram,
      dquads({
          // z quad
          {min.x, min.y, min.z}, {max.x, min.y, min.z}, {max.x, max.y, min.z},
          {min.x, max.y, min.z},  // z quad
          // x quad
          {min.x, min.y, min.z}, {min.x, max.y, min.z}, {min.x, max.y, max.z},
          {min.x, min.y, max.z},  // x quad
          // y quad
          {min.x, min.y, min.z}, {max.x, min.y, min.z}, {max.x, min.y, max.z},
          {min.x, min.y, max.z}  // y quad
      }),
      style);

#if 0
  auto lpositions = vector<vec3f>{};
  for (auto& [x, _] : xticks) {
    lpositions.push_back({x, min.y, 0});
    lpositions.push_back({x, max.y, 0});
  }
  for (auto& [y, _] : yticks) {
    lpositions.push_back({min.x, y, 0});
    lpositions.push_back({max.x, y, 0});
  }
  add_shape(diagram, dlines(lpositions), lstyle);
#endif

  auto tpositions = vector<vec3f>{};
  auto tlabels    = vector<string>{};
  for (auto&& [x, label] : xticks) {
    tpositions.push_back({x, min.y, min.z});
    tlabels.push_back(label + "!!l");
  }
  for (auto&& [y, label] : yticks) {
    tpositions.push_back({min.x, y, min.z});
    tlabels.push_back(label + "!!r");
  }
  for (auto&& [z, label] : zticks) {
    tpositions.push_back({min.x, min.y, z});
    tlabels.push_back(label + "!!r");
  }
  add_labels(diagram, {.positions = tpositions, .labels = tlabels});

  return diagram;
}

// Plot surface
diagram_shape dplotsurface(const function<float(vec2f)>& func, vec2f xrange,
    vec2f yrange, vec2i steps, bool isolines) {
  auto shape = diagram_shape{};

  // steps
  auto stepsp1 = vec2i{steps.x + 1, steps.y + 1};

  // vertices
  shape.positions = vector<vec3f>(stepsp1.x * stepsp1.y);
  shape.texcoords = vector<vec2f>(stepsp1.x * stepsp1.y);
  for (auto j : range(stepsp1.y)) {
    for (auto i : range(stepsp1.x)) {
      auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
      auto x  = lerp(xrange.x, xrange.y, uv.x);
      auto y  = lerp(yrange.x, yrange.y, uv.y);
      auto z  = func({x, y});
      shape.positions[j * stepsp1.x + i] = {x, y, z};
      shape.texcoords[j * stepsp1.x + i] = uv;
    }
  }

  // faces
  shape.quads = vector<vec4i>(steps.x * steps.y);
  for (auto j : range(steps.y)) {
    for (auto i : range(steps.x)) {
      shape.quads[j * steps.x + i] = {j * stepsp1.x + i, j * stepsp1.x + i + 1,
          (j + 1) * stepsp1.x + i + 1, (j + 1) * stepsp1.x + i};
    }
  }

  // no isolines
  if (!isolines) return shape;

  // isolines
  for (auto j : range(0, steps.y + 1, steps.y / 4)) {
    for (auto i : range(steps.x)) {
      shape.lines.push_back({j * stepsp1.x + i, j * stepsp1.x + i + 1});
    }
  }
  for (auto i : range(0, steps.x + 1, steps.x / 4)) {
    for (auto j : range(steps.y)) {
      shape.lines.push_back({j * stepsp1.x + i, (j + 1) * stepsp1.x + i});
    }
  }

  // done
  return shape;
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
//   l : 36 = 25 :  9 => l = 25 * 36 /  9 = 100
static camera_data make_default_camera() {
  return {
      .frame = lookat_frame(vec3f{0, 0, 25}, vec3f{0, 0, 0}, vec3f{0, 1, 0}),
      .lens  = 0.1f,
      // .aspect = 16.0f / 9.0f,
      .aspect = 1,
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
static shape_data make_text_shape(const string& text, vec3f offset,
    const texture_data& texture, float scale) {
  auto size            = max(texture.pixelsf.size(), texture.pixelsb.size());
  auto [width, height] = scale * (vec2f)size / 72.0f;
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
    vec3f camera, vec3f position, float scale) {
  return lookat_frame(position, camera, {0, 1, 0}, true);
}

// Make text material
static material_data make_text_material(vec4f color, int textureid) {
  return {
      .type      = material_type::matte,
      .color     = {color.x, color.y, color.z},
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
static material_data make_fill_material(
    vec4f fill, bool highlight, int textureid) {
  return {
      .type      = highlight ? material_type::glossy : material_type::matte,
      .color     = {fill.x, fill.y, fill.z},
      .opacity   = fill.w,
      .color_tex = textureid,
  };
}

// Make a shape of lines and points
static shape_data make_points_shape(const vector<int>& points,
    const vector<vec3f>& positions_, const frame3f& frame, float radius) {
  auto make_sphere = [](vec3f position, float radius) {
    auto frame = translation_frame(position) *
                 scaling_frame(vec3f{radius * 3, radius * 3, radius * 3});
    auto sphere = make_dsphere(8, 1);
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
    auto capsule = make_duvcapsule({16, 8, 8}, {radius, length(axis) / 2});
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
    auto capsule = make_duvcapsule({16, 8, 8}, {radius, length(axis) / 2});
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
    auto capsule = make_duvcapsule({16, 8, 8}, {radius, length(axis) / 2});
    for (auto& p : capsule.positions) p = transform_point(frame, p);
    merge_shape_inplace(nshape, capsule);
    auto frame2 = translation_frame(position2) *
                  rotation_frame(rotation_axis, rotation_angle);
    auto cone = make_duvcone({8, 8, 8}, {radius * 3, radius * 5});
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
    vec4f stroke, int textureid = invalidid) {
  return {
      .type      = material_type::matte,
      .color     = {stroke.x, stroke.y, stroke.z},
      .opacity   = stroke.w,
      .color_tex = textureid,
  };
}

// Make an image texture
static texture_data make_image_texture(
    const image<vec4f>& image, bool linear, bool nearest) {
  if (linear) {
    return texture_data{.pixelsf = image, .nearest = nearest};
  } else {
    return texture_data{.pixelsb = float_to_byte(image), .nearest = nearest};
  }
}

// Add an object to a scene
static void add_object(scene_data& scene, const diagram_object& object) {
  // unpack
  auto& [frame, shape, labels, style] = object;

  // triangles
  if (!shape.triangles.empty() && style.fill.w > 0) {
    scene.shapes.push_back(make_triangles_shape(
        shape.triangles, shape.positions, shape.texcoords, frame.frame));
    if (style.texture.empty()) {
      scene.materials.push_back(
          make_fill_material(style.fill, style.highlight, invalidid));
    } else {
      scene.textures.emplace_back(
          make_image_texture(style.texture, style.linear, style.nearest));
      scene.materials.push_back(make_fill_material(
          style.fill, style.highlight, (int)scene.textures.size() - 1));
    }
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // quads
  if (!shape.quads.empty() && style.fill.w > 0) {
    scene.shapes.push_back(make_quads_shape(
        shape.quads, shape.positions, shape.texcoords, frame.frame));
    if (style.texture.empty()) {
      scene.materials.push_back(
          make_fill_material(style.fill, style.highlight, invalidid));
    } else {
      scene.textures.emplace_back(
          make_image_texture(style.texture, style.linear, style.nearest));
      scene.materials.push_back(make_fill_material(
          style.fill, style.highlight, (int)scene.textures.size() - 1));
    }
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // points
  if (!shape.points.empty() && style.stroke.w > 0) {
    scene.shapes.push_back(make_points_shape(
        shape.points, shape.positions, frame.frame, style.thickness));
    scene.materials.push_back(make_stroke_material(style.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // lines
  if (!shape.lines.empty() && style.stroke.w > 0) {
    scene.shapes.push_back(make_lines_shape(shape.lines, shape.positions,
        frame.frame, style.thickness, max(shape.connect, style.connect)));
    scene.materials.push_back(make_stroke_material(style.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // arrows
  if (!shape.arrows.empty() && style.stroke.w > 0) {
    scene.shapes.push_back(make_arrows_shape(shape.arrows, shape.positions,
        frame.frame, style.thickness, max(shape.connect, style.connect)));
    scene.materials.push_back(make_stroke_material(style.stroke));
    auto& ninstance    = scene.instances.emplace_back();
    ninstance.shape    = (int)scene.shapes.size() - 1;
    ninstance.material = (int)scene.materials.size() - 1;
    ninstance.frame    = identity3x4f;
  }

  // text
  if ((!labels.labels.empty()) && style.text.w > 0) {
    auto positions = labels.positions;
    auto offsets   = make_text_offsets(labels.labels);
    auto texts     = make_text_labels(labels.labels);
    for (auto idx = 0; idx < (int)texts.size(); idx++) {
      if (texts[idx].empty()) continue;
      auto& texture = scene.textures.emplace_back(
          make_text_texture(texts[idx]));
      scene.materials.push_back(
          make_text_material(style.text, (int)scene.textures.size() - 1));
      scene.shapes.push_back(
          make_text_shape(texts[idx], offsets[idx], texture, style.textscale));
      auto& ninstance    = scene.instances.emplace_back();
      ninstance.frame    = make_text_frame(texts[idx], texture,
             scene.cameras.front().frame.o,
             transform_point(frame.frame, positions[idx]), style.textscale);
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
static image<vec4f> crop_image(const image<vec4f>& source) {
  // find min and max
  auto [width, height] = source.size();
  auto min = 0, max = height;
  for (auto j : range(height)) {
    auto empty = true;
    for (auto i : range(width)) {
      if (source[{i, j}] != vec4f{1, 1, 1, 1}) empty = false;
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
      if (source[{i, height - j - 1}] != vec4f{1, 1, 1, 1}) empty = false;
    }
    if (empty) {
      max = height - j;
    } else {
      break;
    }
  }

  // no content
  if (max < min) return source;

  // crop
  auto cropped = image<vec4f>({width, max - min});
  for (auto j : range(min, max)) {
    for (auto i : range(width)) {
      cropped[{i, j - min}] = source[{i, j}];
    }
  }

  // done
  return cropped;
}

// Image rendering
static image<vec4f> render_image(const scene_data& scene, int resolution,
    int samples, bool noparallel = false);

// Render a diagram to an image
image<vec4f> render_diagram(const diagram_data& diagram, int resolution,
    int samples, bool boxes, bool crop) {
  // final image
  auto composite = image<vec4f>{{resolution, resolution}, {1, 1, 1, 1}};

  // render scenes
  for (auto& scene : diagram.scenes) {
    // convert diagram to scene
    auto yscene = convert_scene(scene);

    // image
    auto render = render_image(yscene, resolution, samples, false);

    // helpers
    auto draw_quad = [](image<vec4f>& composite, vec2i center, vec2i size,
                         vec4f color) {
      auto extents = composite.size();
      for (auto i : range(-size.x / 2, +size.x / 2)) {
        composite[clamp(
            center + vec2i{i, -size.y / 2}, vec2i{0, 0}, extents - 1)] = color;
        composite[clamp(
            center + vec2i{i, +size.y / 2}, vec2i{0, 0}, extents - 1)] = color;
      }
      for (auto j : range(-size.y / 2, +size.y / 2)) {
        composite[clamp(
            center + vec2i{-size.x / 2, j}, vec2i{0, 0}, extents - 1)] = color;
        composite[clamp(
            center + vec2i{+size.x / 2, j}, vec2i{0, 0}, extents - 1)] = color;
      }
    };
    auto copy_quad = [](image<vec4f>& composite, const image<vec4f>& source,
                         vec2i icenter, vec2i scenter, vec2i size) {
      auto csize = composite.size();
      auto ssize = source.size();
      for (auto j : range(-size.y / 2, +size.y / 2)) {
        for (auto i : range(-size.x / 2, +size.x / 2)) {
          auto cij = vec2i{i, j} + icenter, sij = vec2i{i, j} + scenter;
          if (cij.x < 0 || cij.y < 0 || cij.x >= csize.x || cij.y >= csize.y)
            continue;
          if (sij.x < 0 || sij.y < 0 || sij.x >= ssize.x || sij.y >= ssize.y)
            continue;
          composite[cij] = source[sij];
        }
      }
    };

    // composite
    auto size    = (vec2i)(scene.size * resolution / 9.0f);
    auto margin  = (vec2i)(scene.margin * resolution / 9.0f);
    auto icenter = render.size() / 2;
    auto ccenter = composite.size() / 2 +
                   (vec2i)(scene.offset * resolution / 9.0f);

    copy_quad(composite, render, ccenter, icenter, size + margin);
    if (boxes) {
      draw_quad(composite, ccenter, size, {0, 1, 0, 1});
      draw_quad(composite, ccenter, size + margin, {0, 0, 1, 1});
    }
  }

  // crop
  if (crop) composite = crop_image(composite);

  // done
  return composite;
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

// Generates a ray from a camera.
static ray3f eval_camera(const camera_data& camera, vec2f uv_) {
  auto film = vec2f{camera.film, camera.film / camera.aspect};
  if (!camera.orthographic) {
    // point on film
    auto uv = vec2f{1 - uv_.x, uv_.y};
    auto q = vec3f{film.x * (uv.x - 0.5f), film.y * (uv.y - 0.5f), camera.lens};
    // point on the lens
    auto e = vec3f{0, 0, 0};
    // ray direction through the lens center
    auto d = normalize(e - q);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  } else {
    // point on film
    auto uv    = vec2f{1 - uv_.x, uv_.y};
    auto scale = 1 / camera.lens;
    auto q     = vec3f{film.x * (uv.x - 0.5f) * scale,
        film.y * (uv.y - 0.5f) * scale, camera.lens};
    // point on the lens
    auto e = vec3f{-q.x, -q.y, 0};
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
  // lights
  auto light = normalize(vec3f{1, 1, 2});

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
    radiance += weight * material.color * material.opacity;
    if (material.type == material_type::glossy) {
      radiance += weight * pow(max(dot(normal, light), 0.0f), 10) *
                  material.opacity;
    }

    // handle opacity
    if (material.opacity >= 1) break;
    weight *= 1 - material.opacity;
    ray = {position + ray.d * 1e-2f, ray.d};
  }

  return {radiance.x, radiance.y, radiance.z, 1.0f};
}

// Progressively computes an image.
static image<vec4f> render_image(
    const scene_data& scene, int resolution, int samples, bool noparallel) {
  // Bvh
  auto bvh = make_scene_bvh(scene, false, noparallel);
  // Camera
  auto& camera = scene.cameras[0];
  // Image
  auto size   = vec2i{resolution, (int)round(resolution / camera.aspect)};
  auto render = image<vec4f>{size, {0, 0, 0, 0}};
  // Start rendering
  auto nsamples = (int)round(sqrt((float)samples));
  if (noparallel) {
    for (auto ij : range(size)) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto ray = eval_camera(
            camera, ((vec2f)ij + ((vec2f)sij + 0.5f) / (float)nsamples) /
                        (vec2f)render.size());
        auto color = render_ray(scene, bvh, ray);
        render[ij] += isfinite(color) ? color : vec4f{1, 1, 1, 1};
      }
      render[ij] /= nsamples * nsamples;
    }
  } else {
    parallel_for_batch(size, [&](vec2i ij) {
      for (auto sij : range(vec2i{nsamples, nsamples})) {
        auto ray = eval_camera(
            camera, ((vec2f)ij + ((vec2f)sij + 0.5f) / (float)nsamples) /
                        (vec2f)render.size());
        auto color = render_ray(scene, bvh, ray);
        render[ij] += isfinite(color) ? color : vec4f{1, 1, 1, 1};
      }
      render[ij] /= nsamples * nsamples;
    });
  }

  return render;
}

// Helpers to clip lines
diagram_shape clip_lines(const diagram_shape& shape, const bbox3f& bbox) {
  return clip_lines(identity3x4f, shape, bbox);
}

diagram_shape clip_lines(
    const frame3f& frame, const diagram_shape& shape, const bbox3f& bbox) {
  auto inside = vector<bool>(shape.positions.size(), false);
  for (auto index : range(shape.positions.size())) {
    inside[index] = contains(
        bbox, transform_point(frame, shape.positions[index]));
  }
  auto nshape = shape;
  nshape.lines.clear();
  for (auto line : shape.lines) {
    if (inside.at(line.x) && inside.at(line.y)) nshape.lines.push_back(line);
  }
  return nshape;
}

}  // namespace yocto
