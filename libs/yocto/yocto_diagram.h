//
// # Yocto/Diagram: Diagram representation
//
// Yocto/Diagram defines diagram representations.
// Yocto/Diagram is implemented in `yocto_diagram.h` and `yocto_diagram.cpp`.
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
//

#ifndef _YOCTO_DIAGRAM_H_
#define _YOCTO_DIAGRAM_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <yocto/yocto_scene.h>

#include <array>
#include <functional>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::function;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM REPRESENTATION
// -----------------------------------------------------------------------------
namespace yocto {

// Diagram frame
struct diagram_frame {
  frame3f frame = identity3x4f;

  diagram_frame() {}
  diagram_frame(const diagram_frame&) = default;
  diagram_frame(const frame3f& frame_) : frame{frame_} {}
};

// Diagram shape
struct diagram_shape {
  // Element data
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec2i> arrows    = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // Vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};

  // Options
  float connect = 0;
};

// Diagram labels
struct diagram_labels {
  // Labels data
  vector<string> labels    = {};
  vector<string> llabels   = {};
  vector<string> flabels   = {};
  vector<string> clabels   = {};
  vector<vec3f>  positions = {};
};

// Diagram style
struct diagram_style {
  // Style data
  vec4f        stroke    = {0, 0, 0, 1};
  vec4f        fill      = {0.4, 1.0, 1.0, 1};
  vec4f        text      = {0, 0, 0, 1};
  image<vec4f> texture   = {};
  bool         highlight = false;
  bool         arrow     = false;
  bool         nearest   = false;
  bool         linear    = false;
  bool         wireframe = true;
  float        thickness = 0.015f;
  float        textscale = 0.05f;
  float        connect   = 0;
};

// Diagram object
struct diagram_object {
  diagram_frame  frame  = {};
  diagram_shape  shape  = {};
  diagram_labels labels = {};
  diagram_style  style  = {};
};

// Diagram scene
struct diagram_scene {
  // scene rendering
  vec2f size   = {0, 0};
  vec2f offset = {0, 0};
  vec2f margin = {0.1, 0.3};

  // scene transforms
  frame3f frame = identity3x4f;

  // scene objects
  vector<diagram_object> objects = {};
};

// Diagram data
struct diagram_data {
  // Scene data
  vector<diagram_scene> scenes = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM INITIALIZATION AND RENDERING
// -----------------------------------------------------------------------------
namespace yocto {

// Init diagram
diagram_data make_diagram();

// Rendering a diagram
image<vec4f> render_diagram(const diagram_data& diagram, int resolution = 1440,
    int samples = 64, bool boxes = false, bool crop = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Useful constants
namespace dcolors {
inline const auto black       = vec4f{0.0, 0.0, 0.0, 1.0};
inline const auto gray        = vec4f{0.5, 0.5, 0.5, 1.0};
inline const auto white       = vec4f{1.0, 1.0, 1.0, 1.0};
inline const auto transparent = vec4f{1.0, 1.0, 1.0, 0.0};
inline const auto stroke1     = vec4f{0.0, 0.7, 0.7, 1.0};
inline const auto stroke2     = vec4f{0.0, 0.7, 0.0, 1.0};
inline const auto stroke3     = vec4f{0.7, 0.3, 0.0, 1.0};
inline const auto stroke4     = vec4f{0.7, 0.3, 0.7, 1.0};
inline const auto stroke5     = vec4f{0.7, 0.7, 0.0, 1.0};
inline const auto stroke6     = vec4f{0.5, 0.5, 0.5, 1.0};
inline const auto fill1       = vec4f{0.4, 1.0, 1.0, 1.0};
inline const auto fill2       = vec4f{0.4, 1.0, 0.4, 1.0};
inline const auto fill3       = vec4f{1.0, 0.6, 0.2, 1.0};
inline const auto fill4       = vec4f{1.0, 0.6, 1.0, 1.0};
inline const auto fill5       = vec4f{1.0, 1.0, 0.4, 1.0};
inline const auto fill6       = vec4f{0.7, 0.7, 0.7, 1.0};
inline const auto tfill1      = vec4f{0.0, 1.0, 1.0, 0.4};
inline const auto tfill2      = vec4f{0.0, 1.0, 0.0, 0.4};
inline const auto tfill3      = vec4f{1.0, 0.3, 0.0, 0.4};
inline const auto tfill4      = vec4f{1.0, 0.3, 1.0, 0.4};
inline const auto tfill5      = vec4f{1.0, 1.0, 0.0, 0.4};
inline const auto tfill6      = vec4f{0.8, 0.8, 0.8, 0.4};
inline const auto etfill1     = vec4f{0.0, 1.0, 1.0, 0.1};
inline const auto etfill2     = vec4f{0.0, 1.0, 0.0, 0.1};
inline const auto etfill3     = vec4f{1.0, 0.3, 0.0, 0.1};
inline const auto etfill4     = vec4f{1.0, 0.3, 1.0, 0.1};
inline const auto etfill5     = vec4f{1.0, 1.0, 0.0, 0.1};
};  // namespace dcolors

// Thickness
namespace dthickness {
inline const auto default_    = 0.015f;
inline const auto regular     = 0.0144f;
inline const auto medium      = 0.0156f;
inline const auto thin        = 0.010f;
inline const auto thick       = 0.020f;
inline const auto ethin       = 0.005f;
inline const auto extra_thin  = 0.005f;
inline const auto extra_thick = 0.035f;
};  // namespace dthickness

// Points and frames
namespace dconstants {
inline const auto xup3x4f = frame3f{{0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, 0}};
inline const auto yup3x4f = frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
inline const auto zup3x4f = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
inline const auto quad_quads     = vector<vec4i>{{0, 1, 2, 3}};
inline const auto quad_positions = vector<vec3f>{
    {-1, -1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, +1, 0}};
inline const auto quad_texcoords = vector<vec2f>{
    {0, 1}, {1, 1}, {1, 0}, {0, 0}};
inline const auto triangle_triangles = vector<vec3i>{{0, 1, 2}};
inline const auto triangle_positions = vector<vec3f>{
    {-1, -1, 0}, {+1, -1, 0}, {0, +1, 0}};
inline const auto triangle_texcoords = vector<vec2f>{{0, 1}, {1, 1}, {0.5, 0}};
inline const auto cube_quads         = vector<vec4i>{{0, 3, 2, 1}, {4, 5, 6, 7},
            {1, 2, 6, 5}, {0, 4, 7, 3}, {0, 1, 5, 4}, {2, 3, 7, 6}};
inline const auto cube_positions     = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
        {+1, -1, -1}, {-1, -1, -1}, {-1, +1, +1}, {+1, +1, +1}, {+1, +1, -1},
        {-1, +1, -1}};
inline const auto tetra_triangles    = vector<vec3i>{
    {0, 2, 1}, {0, 1, 3}, {0, 3, 2}, {1, 3, 2}};
inline const auto tetra_positions         = vector<vec3f>{{sqrt(2.0f), -1, 0},
            {-sqrt(1.0f / 2.0f), -1, -sqrt(3.0f / 2.0f)},
            {-sqrt(1.0f / 2.0f), -1, +sqrt(3.0f / 2.0f)}, {0, +1, 0}};
inline const auto fourtriangles_triangles = vector<vec3i>{
    {0, 1, 3}, {3, 1, 4}, {4, 1, 2}, {4, 2, 5}};
inline const auto fourtriangles_positions  = vector<vec3f>{{-2.5, -1, 0},
     {-0.5, -1, 0}, {+1.5, -1, 0}, {-1.5, +1, 0}, {+0.5, +1, 0}, {+2.5, +1, 0}};
inline const auto fourtriangles2_triangles = vector<vec3i>{
    {0, 1, 2}, {2, 1, 3}, {2, 3, 4}, {4, 3, 5}};
inline const auto fourtriangles2_positions = vector<vec3f>{
    {-2, -2, 0}, {0, -2, 0}, {-1, 0, 0}, {+1, 0, 0}, {0, +2, 0}, {+2, +2, 0}};
inline const auto slabs_quads     = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
inline const auto slabs_positions = vector<vec3f>{{-2, -2, +1}, {+2, -2, +1},
    {+2, +2, +1}, {-2, +2, +1}, {+2, -2, -1}, {-2, -2, -1}, {-2, +2, -1},
    {+2, +2, -1}, {+1, -2, +2}, {+1, -2, -2}, {+1, +2, -2}, {+1, +2, +2},
    {-1, -2, -2}, {-1, -2, +2}, {-1, +2, +2}, {-1, +2, -2}, {-2, +1, +2},
    {+2, +1, +2}, {+2, +1, -2}, {-2, +1, -2}, {+2, -1, +2}, {-2, -1, +2},
    {-2, -1, -2}, {+2, -1, -2}};
inline const auto slabs_texcoords = vector<vec2f>{{0, 1}, {1, 1}, {1, 0},
    {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0},
    {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1},
    {1, 1}, {1, 0}, {0, 0}};
inline const auto slab_quads      = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7}};
inline const auto slab_positions  = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
     {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
     {+1, +1, -1}};
inline const auto slab_texcoords  = vector<vec2f>{
    {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}};
inline const auto fvcube_quads     = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};
inline const auto fvcube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
    {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
    {+1, +1, -1}};
inline const auto fvcube_normals   = vector<vec3f>{{0, 0, +1}, {0, 0, +1},
      {0, 0, +1}, {0, 0, +1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1},
      {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {-1, 0, 0}, {-1, 0, 0},
      {-1, 0, 0}, {-1, 0, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0},
      {0, -1, 0}, {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
inline const auto fvcube_texcoords = vector<vec2f>{{0, 0}, {1, 0}, {1, 1},
    {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1},
    {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
    {1, 0}, {1, 1}, {0, 1}};
inline const auto fvcube_quadspos  = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
     {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};
inline const auto fvcube_quadsnorm = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
    {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
inline const auto fvcube_quadstexcoord = vector<vec4i>{{0, 1, 2, 3},
    {4, 5, 6, 7}, {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19},
    {20, 21, 22, 23}};
};  // namespace dconstants

// Easy checking in lists
inline bool contains(const string& tag, const vector<string>& tags) {
  for (auto tag_ : tags) {
    if (tag_ == tag) return true;
  }
  return false;
}

// Frames
inline frame3f dtranslation(vec3f translation) {
  return translation_frame(translation);
}
inline frame3f drotation(vec3f rotation) {
  return rotation_frame(vec3f{1, 0, 0}, radians(rotation.x)) *
         rotation_frame(vec3f{0, 1, 0}, radians(rotation.y)) *
         rotation_frame(vec3f{0, 0, 1}, radians(rotation.z));
}
inline frame3f drotation(vec3f axis, float angle) {
  return rotation_frame(axis, radians(angle));
}
inline frame3f dscaling(vec3f scaling) { return scaling_frame(scaling); }
inline frame3f dscaling(float scale) {
  return scaling_frame(vec3f{scale, scale, scale});
}
inline frame3f dscale(float scale) {
  return scaling_frame(vec3f{scale, scale, scale});
}
inline frame3f dframez(vec3f z, vec3f pos = {0, 0, 0}) {
  return frame_fromz(pos, z);
}
inline frame3f dlookat(vec3f from, vec3f to = {0, 0, 0}, vec3f up = {0, 1, 0}) {
  return lookat_frame(from, to, up, true);
}
inline frame3f dtransform(vec3f translation, vec3f rotation) {
  return dtranslation(translation) * drotation(rotation);
}
inline frame3f dgtransform(vec3f translation, vec3f rotation, vec3f scaling) {
  return dtranslation(translation) * drotation(rotation) * dscaling(scaling);
}
inline frame3f dreflectx() {
  return frame3f{{-1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
}
inline frame3f dcentering(vec2f center) {
  return dtranslation({-center.x, -center.y, 0});
}

inline vector<vec3f> transform_points(
    const frame3f& frame, const vector<vec3f>& points) {
  auto tpoints = vector<vec3f>{};
  tpoints.reserve(points.size());
  for (auto& point : points) tpoints.push_back(transform_point(frame, point));
  return tpoints;
}
inline vector<vec3f> transform_normals(
    const frame3f& frame, const vector<vec3f>& normals) {
  auto tnormals = vector<vec3f>{};
  tnormals.reserve(normals.size());
  for (auto& normal : normals)
    tnormals.push_back(transform_normal(frame, normal, true));
  return tnormals;
}
inline diagram_shape transform_shape(
    const frame3f& frame, const diagram_shape& shape) {
  auto tshape = shape;
  for (auto& position : tshape.positions)
    position = transform_point(frame, position);
  return tshape;
}

// Styles
inline diagram_style dstroke(
    vec4f stroke = dcolors::black, float thickness = dthickness::default_) {
  return {
      .stroke = stroke, .fill = dcolors::transparent, .thickness = thickness};
}
inline diagram_style dfill(vec4f fill = dcolors::fill1) {
  return {.stroke = dcolors::transparent, .fill = fill};
}
inline diagram_style dfilled(vec4f fill = dcolors::fill1,
    vec4f stroke = dcolors::black, float thickness = dthickness::default_) {
  return {.stroke = stroke, .fill = fill, .thickness = thickness};
}
inline diagram_style dtextured(const image<vec4f>& texture,
    vec4f stroke = dcolors::black, float thickness = dthickness::default_) {
  return {.stroke = stroke,
      .fill       = {1, 1, 1, 1},
      .thickness  = thickness,
      .texture    = texture};
}
inline diagram_style dtextured(const image<vec4f>& texture, bool interpolate,
    vec4f stroke = dcolors::black, float thickness = dthickness::default_) {
  return {.stroke = stroke,
      .fill       = {1, 1, 1, 1},
      .thickness  = thickness,
      .texture    = texture,
      .nearest    = !interpolate};
}
inline diagram_style dimtextured(const image<vec4f>& texture,
    vec4f stroke    = dcolors::transparent,
    float thickness = dthickness::default_) {
  return {
      .stroke    = stroke,
      .fill      = {1, 1, 1, 1},
      .thickness = thickness,
      .texture   = texture,
      .nearest   = true,
  };
}
inline diagram_style dimtextured(const image<vec4f>& texture, bool interpolate,
    vec4f stroke    = dcolors::transparent,
    float thickness = dthickness::default_) {
  return {.stroke = stroke,
      .fill       = {1, 1, 1, 1},
      .thickness  = thickness,
      .texture    = texture,
      .nearest    = !interpolate};
}
inline diagram_style dtextcolor(const vec4f textcolor = dcolors::black,
    vec4f fill = dcolors::fill1, vec4f stroke = dcolors::black,
    float thickness = dthickness::default_) {
  return {.stroke = stroke,
      .fill       = fill,
      .thickness  = thickness,
      .text       = textcolor};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Add scene
diagram_scene& add_scene(diagram_data& diagram, const string& title,
    const frame3f& frame = identity3x4f, vec2f margin = {0.2, 1.0});
diagram_scene& add_scene(diagram_data& diagram, const string& title, vec2f size,
    const frame3f& frame = identity3x4f, vec2f margin = {0.2, 1.0});
diagram_scene& add_scenews(diagram_data& diagram, const string& title,
    const string& subtitle, vec2f size = {2, 2},
    const frame3f& frame = identity3x4f, vec2f margin = {0.2, 1.0});

// Labels
diagram_labels dlabels(const vector<string>& labels);
diagram_labels dlabels(
    const vector<vec3f>& positions, const vector<string>& labels);
diagram_labels dllabels(const vector<string>& labels);
diagram_labels dflabels(const vector<string>& labels);
diagram_labels dclabels(const vector<string>& labels);
diagram_labels dglabels(const vector<vec3f>& positions,
    const vector<string>& labels, const vector<string>& llabels,
    const vector<string>& flabels, const vector<string>& clabels);
diagram_labels dglabels(const vector<string>& labels,
    const vector<string>& llabels, const vector<string>& flabels,
    const vector<string>& clabels);

// Add labels
void add_labels(diagram_scene& diagram, const diagram_labels& labels,
    const diagram_style& style = dstroke());
void add_labels(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_labels& labels, const diagram_style& style = dstroke());

// Add shapes
void add_shape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_style& style = dfilled());
void add_shape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_labels& labels, const diagram_style& style = dfilled());
void add_shape(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_shape& shape, const diagram_style& style = dfilled());
void add_shape(diagram_scene& diagram, const diagram_frame& frame,
    const diagram_shape& shape, const diagram_labels& labels,
    const diagram_style& style = dfilled());

// Points
diagram_shape dpoints(const vector<vec3f>& positions = {{0, 0, 0}});
diagram_shape dpoints(
    const vector<int>& points, const vector<vec3f>& positions);

// Lines
diagram_shape dlines(const vector<vec3f>& positions);
diagram_shape dlines(
    const vector<vec2i>& lines, const vector<vec3f>& positions);
diagram_shape dclines(const vector<vec3f>& positions);
diagram_shape dclines(
    const vector<vec2i>& lines, const vector<vec3f>& positions);

// Arrows
diagram_shape darrows(const vector<vec3f>& positions);
diagram_shape darrows(
    const vector<vec2i>& arrows, const vector<vec3f>& positions);
diagram_shape dcarrows(const vector<vec3f>& positions);
diagram_shape dcarrows(
    const vector<vec2i>& arrows, const vector<vec3f>& positions);

// Vectors
diagram_shape dvectors(const vector<vec3f>& positions);
diagram_shape dcvectors(const vector<vec3f>& positions);

// Axes
diagram_shape daxes(
    const frame3f& axes = identity3x4f, vec3f aspect = {1, 1, 1});
diagram_shape daxes2(
    const frame3f& axes = identity3x4f, vec3f aspect = {1, 1, 1});

// Rays
diagram_shape drays(const vector<ray3f>& rays = {{{0, 0, 0}, {0, 0, 1}}});
diagram_shape drayscont(
    const vector<ray3f>& rays = {{{0, 0, 0}, {0, 0, 1}}}, float length = 5);

// Quads
diagram_shape dquads(const vector<vec3f>& positions = {{-1, -1, 0}, {+1, -1, 0},
                         {+1, +1, 0}, {-1, +1, 0}},
    const vector<vec2f>&                  texcoords = {});
diagram_shape dquads(const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec2f>& texcoords = {});

// Triangles
diagram_shape dtriangles(
    const vector<vec3f>& positions = {{-1, -1, 0}, {+1, -1, 0}, {0, +1, 0}},
    const vector<vec2f>& texcoords = {});
diagram_shape dtriangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec2f>& texcoords = {});

// Shape
diagram_shape dshape(const shape_data& shape, bool wireframe = true);

// Polyline
diagram_shape dpolyline(const vector<vec3f>& positions);

// Polygon
diagram_shape dpolygon(const vector<vec3f>& positions);  // TODO: remove

// Rect
diagram_shape drect(vec2f aspect = {1, 1});  // TODO: remove

// Disk
diagram_shape ddisk(int steps = 64);

// Half disk
diagram_shape dhalfdisk(int steps = 32);

// Arc
diagram_shape darc(
    vec3f from, vec3f to, vec3f center = {0, 0, 0}, int steps = 16);

// Qbezier
diagram_shape dqbeziers(const vector<vec3f>& positions, int steps = 32);

// Bezier
diagram_shape dbeziers(const vector<vec3f>& positions, int steps = 32);

// Cube
diagram_shape dcube();
diagram_shape dcube(int steps);

// Sphere
diagram_shape dsphere(int steps = 32, bool isolines = true);
diagram_shape duvsphere(int steps = 32, bool isolines = true);

// Bbox
diagram_shape dbbox(
    const bbox3f& bbox = {{-1, -1, -1}, {+1, +1, +1}}, float epsilon = 0.01);

// Hemisphere
diagram_shape dhemisphere(int steps = 32, bool isolines = true);

// Cone
diagram_shape dcone(int steps = 32, bool isolines = true);

// Cylinder
diagram_shape duvcylinder(float hscale = 2, int rsteps = 32, int hsteps = 2,
    int csteps = 4, bool isolines = true);

// Capsule
diagram_shape duvcapsule(float hscale = 2, int rsteps = 32, int hsteps = 2,
    int csteps = 4, bool isolines = true);

// Dot grid
diagram_shape ddotgrid(vec2i steps = {4, 4});

// Grid
diagram_shape dgrid(vec2i steps = {4, 4});

// Disk grid
diagram_shape ddiskgrid(vec2i steps = {4, 4}, int dstep = 32);
diagram_shape dudiskgrid(vec2i steps = {4, 4}, int dstep = 32);

// Affine grid
diagram_shape daffinegrid(
    vec3f axes_a = {1, 0, 0}, vec3f axes_b = {0, 1, 0}, vec2i steps = {4, 4});

// Image
diagram_shape dimagerect(vec2i size);
diagram_shape dimagerect(const image<vec4f>& image);
diagram_shape dimagegrid(vec2i size);
diagram_shape dimagegrid(const image<vec4f>& image);
diagram_shape dimagelabel(vec2i size, float scale = 1);
diagram_shape dimagelabel(const image<vec4f>& image, float scale = 1);

// Random points
enum struct drandompoints_type {
  // clang-format off
  quad, disk, disknu, triangle, hemi, hemicos, hemicospower, sphere,
  linex, liney
  // clang-format on
};
diagram_shape drandompoints(int num = 16, bool stratified = true);
diagram_shape drandompoints(
    drandompoints_type type, int num = 16, bool stratified = true);

// Random points
enum struct drandompoints3_type {
  // clang-format off
  cube, linex, liney, linez
  // clang-format on
};
diagram_shape drandompoints3(int num = 16, bool stratified = true);
diagram_shape drandompoints3(
    drandompoints3_type type, int num = 16, bool stratified = true);

// Random lines
enum struct drandomlines_type { hemi, hemicos, hemicospower, beam };
diagram_shape drandomlines(int num = 16, bool stratified = true);
diagram_shape drandomlines(
    drandomlines_type type, int num = 16, bool stratified = true);

// Add plot
diagram_scene& add_plot(diagram_data& diagram, const string& title,
    vec2f size = {4, 2}, const frame3f& frame = identity3x4f,
    vec2f margin = {0.8, 1.0});

// Add plot axes
diagram_scene& add_plotaxes(diagram_scene& diagram,
    const bbox2f&                          bounds = {{0, 0}, {1, 1}},
    const vector<pair<float, string>>&     xticks = {},
    const vector<pair<float, string>>&     yticks = {},
    const diagram_style&                   style  = dstroke(),
    const diagram_style& lstyle = dstroke(dcolors::transparent));

// Add plot shape
void add_plotshape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_style& style = dstroke(dcolors::black));
void add_plotshape(diagram_scene& diagram, const diagram_shape& shape,
    const diagram_labels& labels,
    const diagram_style&  style = dstroke(dcolors::black));

// Plot function
diagram_shape dplotline(const vector<vec2f>& points);
diagram_shape dplotpoints(const vector<vec2f>& points);
diagram_shape dplotarrows(const vector<vec2f>& points);
vector<vec2f> dplotfunc(const function<float(float)>& func,
    vec2f range = {0, 1}, int samples = 100);
vector<vec2f> dplotpfunc(const function<float(float)>& func, int samples = 100);
vector<vec2f> dplotcurve(const vector<float>& curve, bool center = true);
vector<vec2f> dplotcurve(const vector<float>& curve, vec2f range);
inline diagram_shape dplotline(const function<float(float)>& func,
    vec2f range = {0, 1}, int samples = 100) {
  return dplotline(dplotfunc(func, range, samples));
}
inline diagram_shape dplotpoints(const function<float(float)>& func,
    vec2f range = {0, 1}, int samples = 100) {
  return dplotpoints(dplotfunc(func, range, samples));
}
inline diagram_shape dplotline(const vector<float>& curve, bool center = true) {
  return dplotline(dplotcurve(curve, center));
}
inline diagram_shape dplotpoints(
    const vector<float>& curve, bool center = true) {
  return dplotpoints(dplotcurve(curve, center));
}

// Add plot3
diagram_scene& add_plot3(diagram_data& diagram, vec2f size = {4, 3},
    const frame3f& frame = drotation({25, 0, 0}) * drotation({0, 45, 0}) *
                           drotation({270, 0, 0}) * drotation({0, 0, 180}),
    vec2f margin = {0.8, 1.0});
diagram_scene& add_plot3(diagram_data& diagram, const string& title,
    vec2f          size  = {4, 3},
    const frame3f& frame = drotation({25, 0, 0}) * drotation({0, 45, 0}) *
                           drotation({270, 0, 0}) * drotation({0, 0, 180}),
    vec2f margin = {0.8, 1.0});

// Add plot axes
diagram_scene& add_plotaxes3(diagram_scene& diagram,
    const bbox3f& bounds = {{0, 0, 0}, {1, 1, 1}}, vec3f size = {1, 1, 1},
    const vector<pair<float, string>>& xticks = {},
    const vector<pair<float, string>>& yticks = {},
    const vector<pair<float, string>>& zticks = {},
    const diagram_style&               style  = dstroke(),
    const diagram_style&               lstyle = dstroke(dcolors::transparent));

// Plot surface
diagram_shape dplotsurface(const function<float(vec2f)>& func,
    vec2f xrange = {0, 1}, vec2f yrange = {0, 1}, vec2i steps = {64, 64},
    bool wireframe = true);

// Helpers to clip lines
diagram_shape clip_lines(const diagram_shape& shape, const bbox3f& bbox);
diagram_shape clip_lines(
    const frame3f& frame, const diagram_shape& shape, const bbox3f& bbox);

}  // namespace yocto

#endif
