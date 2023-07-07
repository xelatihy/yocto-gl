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
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM REPRESENTATION
// -----------------------------------------------------------------------------
namespace yocto {

// Diagram object
struct diagram_object {
  // Object data
  frame3f frame = identity3x4f;

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

  // Labels data
  vector<string> labels     = {};
  vector<vec3f>  lpositions = {};

  // Style data
  vec4f          stroke    = {0, 0, 0, 1};
  vec4f          fill      = {0.4, 1.0, 1.0, 1};
  vec4f          text      = {0, 0, 0, 1};
  array2d<vec4f> texture   = {};
  bool           arrow     = false;
  bool           nearest   = false;
  bool           wireframe = true;
  float          thickness = 0.015f;
  float          textscale = 0.05f;
  float          connect   = 0;
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

// Diagram function ref
template <typename Ret, typename Arg>
struct diagram_func {
  diagram_func(const diagram_func&) = default;
  template <typename Func>
  diagram_func(Func&& func) {
    _ptr      = (void*)&func;
    _callback = [](void* _ptr, Arg x) -> Ret { return (*(Func*)_ptr)(x); };
  }
  Ret operator()(Arg x) const { return _callback(_ptr, x); }

 private:
  void* _ptr                           = nullptr;
  Ret (*_callback)(void* ptr, Arg arg) = nullptr;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM INITIALIZATION AND RENDERING
// -----------------------------------------------------------------------------
namespace yocto {

// Init diagram
diagram_data make_diagram();

// Rendering a diagram
array2d<vec4f> render_diagram(
    const diagram_data& diagram, bool draft = false, bool crop = true);

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
inline const auto fill1       = vec4f{0.4, 1.0, 1.0, 1.0};
inline const auto fill2       = vec4f{0.4, 1.0, 0.4, 1.0};
inline const auto fill3       = vec4f{1.0, 0.6, 0.2, 1.0};
inline const auto fill4       = vec4f{1.0, 0.6, 1.0, 1.0};
inline const auto fill5       = vec4f{1.0, 1.0, 0.4, 1.0};
inline const auto tfill1      = vec4f{0.0, 1.0, 1.0, 0.4};
inline const auto tfill2      = vec4f{0.0, 1.0, 0.0, 0.4};
inline const auto tfill3      = vec4f{1.0, 0.3, 0.0, 0.4};
inline const auto tfill4      = vec4f{1.0, 0.3, 1.0, 0.4};
inline const auto tfill5      = vec4f{1.0, 1.0, 0.0, 0.4};
};  // namespace dcolors

// Thickness
namespace dthickness {
inline const auto default_    = 0.015f;
inline const auto medium      = 0.02f;
inline const auto thin        = 0.01f;
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
inline const auto triangle_texcoords = vector<vec2f>{{0, 0}, {1, 0}, {0.5, 1}};
inline const auto cube_quads         = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
            {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};
inline const auto cube_positions     = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
        {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
        {+1, +1, -1}};
inline const auto tetra_triangles    = vector<vec3i>{
    {0, 1, 2}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};
inline const auto tetra_positions = vector<vec3f>{{sqrt(2.0f), -1, 0},
    {-sqrt(1.0f / 2.0f), -1, -sqrt(3.0f / 2.0f)},
    {-sqrt(1.0f / 2.0f), -1, +sqrt(3.0f / 2.0f)}, {0, +1, 0}};
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
inline const auto slab_positions  = vector<vec3f>{{-2, -2, +1}, {+2, -2, +1},
     {+2, +2, +1}, {-2, +2, +1}, {+2, -2, -1}, {-2, -2, -1}, {-2, +2, -1},
     {+2, +2, -1}};
inline const auto slab_texcoords  = vector<vec2f>{
    {0, 1}, {1, 1}, {1, 0}, {0, 0}, {0, 1}, {1, 1}, {1, 0}, {0, 0}};
inline const auto fvcube_quads     = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {1, 4, 7, 2}, {5, 0, 3, 6}, {3, 2, 7, 6}, {1, 0, 5, 4}};
inline const auto fvcube_positions = vector<vec3f>{{-1, -1, +1}, {+1, -1, +1},
    {+1, +1, +1}, {-1, +1, +1}, {+1, -1, -1}, {-1, -1, -1}, {-1, +1, -1},
    {+1, +1, -1}};
};  // namespace dconstants

// Easy checking in lists
inline bool contains(const string& tag, const vector<string>& tags) {
  for (auto tag_ : tags) {
    if (tag_ == tag) return true;
  }
  return false;
}

// Frames
inline frame3f dtranslation(const vec3f& translation) {
  return translation_frame(translation);
}
inline frame3f drotation(const vec3f& rotation) {
  return rotation_frame(vec3f{1, 0, 0}, radians(rotation.x)) *
         rotation_frame(vec3f{0, 1, 0}, radians(rotation.y)) *
         rotation_frame(vec3f{0, 0, 1}, radians(rotation.z));
}
inline frame3f dscaling(const vec3f& scaling) { return scaling_frame(scaling); }
inline frame3f dscaling(float scale) {
  return scaling_frame(vec3f{scale, scale, scale});
}
inline frame3f dlookat(const vec3f& from, const vec3f& to = {0, 0, 0},
    const vec3f& up = {0, 1, 0}) {
  return lookat_frame(from, to, up, true);
}
inline frame3f dtransform(const vec3f& translation, const vec3f& rotation) {
  return dtranslation(translation) * drotation(rotation);
}
inline frame3f dgtransform(
    const vec3f& translation, const vec3f& rotation, const vec3f& scaling) {
  return dtranslation(translation) * drotation(rotation) * dscaling(scaling);
}

// Sampling functions
inline vector<vec2f> sample_function(const diagram_func<float, float>& function,
    const vec2f& range_, int samples = 100) {
  auto data = vector<vec2f>(samples);
  for (auto idx : range(samples)) {
    auto x    = lerp(range_.x, range_.y, (float)idx / (float)(samples - 1));
    data[idx] = {x, function(x)};
  }
  return data;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// DIAGRAM CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Add scene
struct diagram_scene_params {
  vec2f   size        = {2, 2};
  vec2f   center      = {0, 0};
  string  title       = "";
  string  subtitle    = "";
  frame3f frame       = identity3x4f;
  vec2f   margin      = {0.2, 1.0};
  vec2f   spacing     = {0.15, 0.15};
  bool    auto_offset = true;
};
diagram_scene& add_scene(
    diagram_data& diagram, const diagram_scene_params& params);

// Add labels
struct diagram_labels_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {};
  vector<string> labels    = {};
};
void add_labels(diagram_scene& diagram, const diagram_labels_params& params);

// Add points
struct diagram_points_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, 0, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_points(diagram_scene& diagram, const diagram_points_params& params);

// Add lines
struct diagram_lines_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, -1, 0}, {0, 1, 0}};
  vector<vec2i>  lines     = {};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  bool           connector = false;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_lines(diagram_scene& diagram, const diagram_lines_params& params);

// Add arrows
struct diagram_arrows_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, -1, 0}, {0, 1, 0}};
  vector<vec2i>  arrows    = {};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  bool           connector = false;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_arrows(diagram_scene& diagram, const diagram_arrows_params& params);

// Add vectors
struct diagram_vectors_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  vectors   = {{0, 1, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_vectors(diagram_scene& diagram, const diagram_vectors_params& params);

// Add axes
struct diagram_axes_params {
  frame3f        frame     = identity3x4f;
  frame3f        axes      = identity3x4f;
  vec3f          aspect    = {1, 1, 1};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_axes(diagram_scene& diagram, const diagram_axes_params& params);

// Add rays
struct diagram_rays_params {
  frame3f        frame      = identity3x4f;
  vector<ray3f>  rays       = {{{0, 0, 0}, {0, 0, 1}}};
  vec3f          position   = {0, 0, 0};
  float          scale      = 1;
  bool           connector  = false;
  vector<string> labels     = {};
  vector<string> clabels    = {};
  vec4f          stroke     = dcolors::black;
  float          thickness  = dthickness::default_;
  float          llength    = 0;
  vec4f          lstroke    = dcolors::gray;
  float          lthickness = dthickness::default_;
  vector<string> llabels    = {};
  vector<string> lclabels   = {};
};
void add_rays(diagram_scene& diagram, const diagram_rays_params& params);

// Add quads
struct diagram_quads_params {
  frame3f       frame     = identity3x4f;
  vector<vec3f> positions = {
      {-1, -1, 0}, {+1, -1, 0}, {+1, +1, 0}, {-1, +1, 0}};
  vector<vec2f>  texcoords = {};
  vector<vec4i>  quads     = {};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
  array2d<vec4f> texture   = {};
};
void add_quads(diagram_scene& diagram, const diagram_quads_params& params);

// Add triangles
struct diagram_triangles_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{-1, -1, 0}, {+1, -1, 0}, {0, +1, 0}};
  vector<vec2f>  texcoords = {};
  vector<vec3i>  triangles = {};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
  array2d<vec4f> texture   = {};
};
void add_triangles(
    diagram_scene& diagram, const diagram_triangles_params& params);

// Add polyline
struct diagram_polyline_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, -1, 0}, {0, +1, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_polyline(
    diagram_scene& diagram, const diagram_polyline_params& params);

// Add polygon
struct diagram_polygon_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{-1, -1, 0}, {+1, -1, 0}, {0, +1, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
};
void add_polygon(diagram_scene& diagram, const diagram_polygon_params& params);

// Add rect
struct diagram_rect_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  vec2f          aspect    = {1, 1};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
  array2d<vec4f> texture   = {};
};
void add_rect(diagram_scene& diagram, const diagram_rect_params& params);

// Add disk
struct diagram_disk_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  vec2f          aspect    = {1, 1};
  float          scale     = 1;
  int            steps     = 32;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
};
void add_disk(diagram_scene& diagram, const diagram_disk_params& params);

// Add arc
struct diagram_arc_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  vec3f          from      = {1, 0, 0};
  vec3f          to        = {0, 1, 0};
  float          scale     = 1;
  int            steps     = 16;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
};
void add_arc(diagram_scene& diagram, const diagram_arc_params& params);

// Add qbezier
struct diagram_qbeziers_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, -1, 0}, {0, 1, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  int            steps     = 32;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_qbeziers(
    diagram_scene& diagram, const diagram_qbeziers_params& params);

// Add bezier
struct diagram_beziers_params {
  frame3f        frame     = identity3x4f;
  vector<vec3f>  positions = {{0, -1, 0}, {0, 1, 0}};
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  int            steps     = 32;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
void add_beziers(diagram_scene& diagram, const diagram_beziers_params& params);

// Add measure
struct diagram_measure_params {
  frame3f        frame     = identity3x4f;
  vec3f          from      = {0, 0, 0};
  vec3f          to        = {1, 0, 0};
  vec3f          offset    = {0, 1, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
};
void add_measure(diagram_scene& diagram, const diagram_measure_params& params);

// Add cube
struct diagram_cube_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
};
void add_cube(diagram_scene& diagram, const diagram_cube_params& params);

// Add sphere
struct diagram_sphere_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  bool           uvsphere  = false;
  int            steps     = 32;
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::fill1;
  float          thickness = dthickness::default_;
  array2d<vec4f> texture   = {};
};
void add_sphere(diagram_scene& diagram, const diagram_sphere_params& params);

// Add grid
struct diagram_grid_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  float          scale     = 1;
  vec2i          steps     = {4, 4};
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::transparent;
  float          thickness = dthickness::default_;
};
void add_grid(diagram_scene& diagram, const diagram_grid_params& params);

// Add affine grid
struct diagram_affinegrid_params {
  frame3f        frame     = identity3x4f;
  vec3f          position  = {0, 0, 0};
  vec3f          axes_a    = {1, 0, 0};
  vec3f          axes_b    = {0, 1, 0};
  float          scale     = 1;
  vec2i          steps     = {4, 4};
  vector<string> labels    = {};
  vector<string> clabels   = {};
  vec4f          stroke    = dcolors::black;
  vec4f          fill      = dcolors::transparent;
  float          thickness = dthickness::default_;
};
void add_affinegrid(
    diagram_scene& diagram, const diagram_affinegrid_params& params);

// Add image
struct diagram_image_params {
  frame3f        frame       = identity3x4f;
  array2d<vec4f> image       = {};
  vec3f          position    = {0, 0, 0};
  float          scale       = 1;
  bool           scaley      = false;
  bool           interpolate = false;
  vector<string> labels      = {};
  vector<string> clabels     = {};
  vector<string> olabels     = {};
  vec4f          stroke      = dcolors::transparent;
  vec4f          fill        = dcolors::white;
  float          thickness   = dthickness::default_;
};
void add_image(diagram_scene& diagram, const diagram_image_params& params);

// Add image grid
struct diagram_imagegrid_params {
  frame3f        frame       = identity3x4f;
  array2d<vec4f> image       = {};
  vec3f          position    = {0, 0, 0};
  float          scale       = 1;
  bool           interpolate = false;
  vector<string> labels      = {};
  vector<string> clabels     = {};
  vec4f          stroke      = dcolors::black;
  vec4f          fill        = dcolors::white;
  float          thickness   = dthickness::default_;
};
void add_imagegrid(
    diagram_scene& diagram, const diagram_imagegrid_params& params);

// Add image mosaic
struct diagram_imagemosaic_params {
  frame3f                frame       = identity3x4f;
  vector<array2d<vec4f>> images      = {};
  vec3f                  position    = {0, 0, 0};
  float                  scale       = 1;
  bool                   scaley      = false;
  bool                   interpolate = false;
  vector<string>         labels      = {};
  vector<string>         clabels     = {};
  vector<string>         olabels     = {};
  vec4f                  stroke      = dcolors::transparent;
  vec4f                  fill        = dcolors::white;
  float                  thickness   = dthickness::default_;
};
void add_imagemosaic(
    diagram_scene& diagram, const diagram_imagemosaic_params& params);

// Add random points
enum struct diagram_randompoints_type { quad, disk, disknu, triangle };
struct diagram_randompoints_params {
  frame3f                   frame      = identity3x4f;
  vec3f                     position   = {0, 0, 0};
  float                     scale      = 1;
  diagram_randompoints_type type       = diagram_randompoints_type::quad;
  int                       steps      = 16;
  bool                      stratified = false;
  vector<vec3f>             triangle   = dconstants::triangle_positions;
  vector<string>            labels     = {};
  vec4f                     stroke     = dcolors::black;
  float                     thickness  = dthickness::default_;
};
void add_randompoints(
    diagram_scene& diagram, const diagram_randompoints_params& params);

// Add random lines
enum struct diagram_randomlines_type { hemi, hemicos, hemicospower, beam };
struct diagram_randomlines_params {
  frame3f                  frame      = identity3x4f;
  vec3f                    position   = {0, 0, 0};
  float                    scale      = 1;
  diagram_randomlines_type type       = diagram_randomlines_type::hemi;
  int                      steps      = 16;
  bool                     stratified = false;
  vector<string>           labels     = {};
  vector<string>           clabels    = {};
  vec4f                    stroke     = dcolors::black;
  float                    thickness  = dthickness::default_;
};
void add_randomlines(
    diagram_scene& diagram, const diagram_randomlines_params& params);

// Add plot
struct diagram_plot_params {
  vec2f   size        = {4, 2};
  vec2f   center      = {0, 0};
  string  title       = "";
  string  subtitle    = "";
  frame3f frame       = identity3x4f;
  vec2f   margin      = {0.8, 1.0};
  vec2f   spacing     = {0.15, 0.15};
  bool    auto_offset = true;
};
diagram_scene& add_plot(
    diagram_data& diagram, const diagram_plot_params& params);

// Add plot axes
struct diagram_plotaxes_params {
  string         title     = "";
  vec2f          xbounds   = {0, 1};
  vec2f          ybounds   = {0, 1};
  vector<string> xlabels   = {};
  vector<float>  xticks    = {};
  vector<string> ylabels   = {};
  vector<float>  yticks    = {};
  vec4f          stroke    = dcolors::black;
  float          thickness = dthickness::default_;
};
diagram_scene& add_plotaxes(
    diagram_scene& diagram, const diagram_plotaxes_params& params);

// Plot functions
struct diagram_plotdata_params {
  frame3f                    frame     = identity3x4f;
  vector<vec2f>              points    = {};
  diagram_func<float, float> function  = [](float x) { return x; };
  vec2f                      range     = {0, 1};
  int                        steps     = 100;
  vector<string>             labels    = {};
  vec4f                      stroke    = dcolors::black;
  float                      thickness = dthickness::default_;
};
void add_plotdata(diagram_scene& diagram, const diagram_plotdata_params& frame);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Bspline subdivision
pair<vector<vec2i>, vector<float>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<float>& vertices);
pair<vector<vec2i>, vector<vec2f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec2f>& vertices);
pair<vector<vec2i>, vector<vec3f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec3f>& vertices);
pair<vector<vec2i>, vector<vec4f>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<vec4f>& vertices);

// Subdivide lines by splitting each line in half.
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline_(
    const vector<vec2i>& lines, const vector<T>& vertices, int level) {
  if (level < 1) return {lines, vertices};
  auto tess = pair{lines, vertices};
  for (auto idx : range(level))
    tess = subdivide_bspline(tess.first, tess.second);
  return tess;
}

// Bspline subdivision
pair<vector<vec4i>, vector<float>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<float>& vertices,
    bool uncorrected = false);
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec2f>& vertices,
    bool uncorrected = false);
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec3f>& vertices,
    bool uncorrected = false);
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<vec4f>& vertices,
    bool uncorrected = false);

template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<T>& vertices, int level,
    bool uncorrected = false);

// Subdivide lines by splitting each line in half.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_catmullclark_(
    const vector<vec4i>& quads, const vector<T>& vertices, int level,
    bool uncorrected) {
  if (level < 1) return {quads, vertices};
  auto tess = pair{quads, vertices};
  for (auto idx : range(level))
    tess = subdivide_catmullclark_(tess.first, tess.second, uncorrected);
  return tess;
}

}  // namespace yocto

#endif
