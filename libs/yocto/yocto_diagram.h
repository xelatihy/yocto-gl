//
// # Yocto/Diagram: Diagram representation
//
// Yocto/Diagram defines diagram representations.
// Yocto/Diagram is implemented in `yocto_diagram.h` and `yocto_diagram.cpp`.
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
// PATTERN REPRESENTATION
// -----------------------------------------------------------------------------
namespace yocto {

// Diagram style
struct style_data {
  vec4f          stroke    = {0, 0, 0, 1};
  vec4f          fill      = {0.4, 1.0, 1.0, 1};
  vec4f          text      = {0, 0, 0, 1};
  vector<vec4f>  fcolors   = {};
  array2d<vec4f> texture   = {};
  bool           arrow     = false;
  bool           nearest   = false;
  bool           wireframe = true;
  float          thickness = 0.015f;
  float          textscale = 0.05f;
  float          connect   = 0;
};

// Diagram text
struct label_data {
  vector<vec3f>  positions = {};
  vector<string> labels    = {};
};

// Diagram data
struct diagram_data {
  scene_data scene = {};
};

// Diagram context
struct diagram_context {
  diagram_data& diagram;
  scene_data&   scene;
  frame3f       frame = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
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

// Diagram callback
using diagram_callback = diagram_func<void, const diagram_context&>;

// Render a diagram to an image
array2d<vec4f> render_diagram(
    const diagram_data& diagram, bool draft = false, bool crop = true);
array2d<vec4f> render_diagram(
    const diagram_callback& callback, bool draft = false, bool crop = true);

#ifdef YOCTO_OPENGL

// Render a diagram interactively
void view_diagram(const string& name, const diagram_data& diagram);
void view_diagram(const string& name, const diagram_callback& callback);

#endif

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

// Group nodes
frame3f add_group(diagram_data& diagram, const vec3f& offset = {0, 0, 0},
    const frame3f& frame = identity3x4f);
frame3f add_tgroup(diagram_data& diagram, const frame3f& frame);

// Point transformations
inline vector<vec3f> translate(
    const vec3f& translation, const vector<vec3f>& positions) {
  auto transformed = positions;
  for (auto& position : transformed) position += translation;
  return transformed;
}

// Shape nodes
void add_label(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, const string& label,
    const vec4f& textcolor = dcolors::black);
void add_labels(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels,
    const vec4f& textcolor = dcolors::black);
void add_point(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, const string& label = "",
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_points(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& color     = dcolors::black,
    float        thickness = dthickness::default_);
void add_points(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& positions, const vector<string>& labels = {},
    const vec4f& color     = dcolors::black,
    float        thickness = dthickness::default_);
void add_line(diagram_data& diagram, const frame3f& frame,
    const vec3f& positions1, const vec3f& position2,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_lines(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_lines(diagram_data& diagram, const frame3f& frame,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_arrow(diagram_data& diagram, const frame3f& frame,
    const vec3f& positions1, const vec3f& position2,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_arrows(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_cline(diagram_data& diagram, const frame3f& frame,
    const vec3f& positions1, const vec3f& position2,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_carrow(diagram_data& diagram, const frame3f& frame,
    const vec3f& positions1, const vec3f& position2,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_vector(diagram_data& diagram, const frame3f& frame,
    const vec3f& direction, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_axes(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_axes(diagram_data& diagram, const frame3f& frame,
    const frame3f& frame1, float scale, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_ray(diagram_data& diagram, const frame3f& frame, const vec3f& position,
    const vec3f& direction, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_polyline(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_quad(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_quadv(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_tquad(diagram_data& diagram, const frame3f& frame,
    const array2d<vec4f>& texture, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rect(diagram_data& diagram, const frame3f& frame, const vec2f& aspect,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_triangle(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_trianglev(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_polygon(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_triangles(diagram_data& diagram, const frame3f& frame,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<string>& labels = {}, const vec4f& fill = dcolors::fill1,
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_quads(diagram_data& diagram, const frame3f& frame,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<string>& labels = {}, const vec4f& fill = dcolors::fill1,
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_disk(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_circle(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_arc(diagram_data& diagram, const frame3f& frame, const vec3f& position,
    float scale, float angle, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_cube(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_sphere(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill      = dcolors::fill1,
    const vec4f& stroke    = dcolors::transparent,
    float        thickness = dthickness::default_);
void add_uvsphere(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill      = dcolors::fill1,
    const vec4f& stroke    = dcolors::transparent,
    float        thickness = dthickness::default_);
void add_uvspheret(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& texture,
    const vector<string>& labels    = {},
    const vec4f&          stroke    = dcolors::transparent,
    float                 thickness = dthickness::default_);
void add_grid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_cgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_cgridnu(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_tgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_tgridv(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_bezier(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_qbezier(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& positions, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_slabs(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_slab(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_image(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& image,
    const vector<string>& labels    = {},
    const vec4f&          stroke    = dcolors::transparent,
    float                 thickness = dthickness::default_);
void add_imagegrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const array2d<vec4f>& image,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_imagefull(diagram_data& diagram, const frame3f& frame,
    const array2d<vec4f>& image, const string& label = "");
void add_imagemosaic(diagram_data& diagram, const frame3f& frame,
    const vector<array2d<vec4f>>& images, const vector<string>& labels = {});
void add_rpoints(diagram_data& diagram, const frame3f& frame,
    const diagram_func<vec3f, vec2f>& func, int steps, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rdpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rdpointsnu(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rtpoints(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, int steps, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_rtpoints(diagram_data& diagram, const frame3f& frame,
    const vector<vec3f>& triangle, int steps = 32, bool stratified = false,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
enum struct rlines_type { hemi, hemicos, hemicospower, beam };
void add_rlines(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, rlines_type type, int steps,
    bool stratified = false, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_atgrid(diagram_data& diagram, const frame3f& frame,
    const vec3f& position, float scale, const vec2i& steps,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);

// Shape nodes
void add_lines(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_shape(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels = {},
    const vec4f& fill = dcolors::fill1, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_shape(diagram_data& diagram, const frame3f& frame,
    const shape_data& shape, const vector<string>& labels,
    const array2d<vec4f>& texture, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);

// Easy checking in lists
bool contains(const string& tag, const vector<string>& tags);

// Frames
frame3f dtranslation(const vec3f& translation);
frame3f drotation(const vec3f& rotation);
frame3f dscaling(const vec3f& scale);
frame3f dscaling(float scale);
frame3f dlookat(const vec3f& from, const vec3f& to = {0, 0, 0},
    const vec3f& up = {0, 1, 0});
frame3f dtransform(const vec3f& translation, const vec3f& rotation);
frame3f dgtransform(
    const vec3f& translation, const vec3f& rotation, const vec3f& scaling);

// Plot axes
struct dplotaxes_params {
  vec2f          aspect  = {2, 1};
  vec2f          xbounds = {0, 1};
  vec2f          ybounds = {0, 1};
  bool           polar   = false;
  vector<string> titles  = {};
  vector<string> xlabels = {};
  vector<float>  xticks  = {};
  vector<string> ylabels = {};
  vector<float>  yticks  = {};
};

// Plot functions
frame3f add_plot(
    diagram_data& diagram, const frame3f& frame, const dplotaxes_params& axes);
void add_lineplot(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& points, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_lineplot(diagram_data& diagram, const frame3f& frame,
    const diagram_func<float, float>& function, const vec2f& range,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);
void add_scatterplot(diagram_data& diagram, const frame3f& frame,
    const vector<vec2f>& points, const vector<string>& labels = {},
    const vec4f& stroke    = dcolors::black,
    float        thickness = dthickness::default_);
void add_scatterplot(diagram_data& diagram, const frame3f& frame,
    const diagram_func<float, float>& function, const vec2f& range,
    const vector<string>& labels = {}, const vec4f& stroke = dcolors::black,
    float thickness = dthickness::default_);

// Sampling functions
vector<vec2f> sample_function(const diagram_func<float, float>& function,
    const vec2f& range, int samples = 100);

// Color shapes
array2d<vec4f> random_image(const vec2i& steps);
array2d<vec4f> gamma_image(const vec2i& steps);

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

template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vertices, int level);

// Subdivide lines by splitting each line in half.
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
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

// -----------------------------------------------------------------------------
// EXAMPLE DIAGRAMS
// -----------------------------------------------------------------------------
namespace yocto {

inline auto  placeholder_tags = vector<string>{"missing"};
diagram_data placeholder_diagrams(const string& tag);

inline auto frame_tags = vector<string>{
    "composition", "coordinates", "definition"};
diagram_data frame_diagrams(const string& tag);

inline auto  objtransform_tags = vector<string>{"rotationx", "rotationy",
     "rotationz", "scalingnu", "scalingu", "translation"};
diagram_data objtransform_diagrams(const string& tag);

inline auto  primitive_tags = vector<string>{"types"};
diagram_data primitive_diagrams(const string& tag);

inline auto transform_tags = vector<string>{
    "rotation", "scaling", "translation"};
diagram_data transform_diagrams(const string& tag);

inline auto  vector_tags = vector<string>{"operations", "products"};
diagram_data vector_diagrams(const string& tag);

inline auto compositing_tags = vector<string>{
    "alpha1", "alpha2", "alpha3", "alpha4"};
diagram_data compositing_diagrams(const string& tag);

inline auto  image_tags = vector<string>{"approximation", "gamma", "grid"};
diagram_data image_diagrams(const string& tag);

}  // namespace yocto

#endif
