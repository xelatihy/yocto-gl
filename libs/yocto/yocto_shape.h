//
// # Yocto/Shape: Shape utilities
//
// Yocto/Shape is a collection of utilities for manipulating shapes in 3D
// graphics, with a focus on triangle and quad meshes.
// Yocto/Shape is implemented in `yocto_shape.h` and `yocto_shape.cpp`.
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

#ifndef _YOCTO_SHAPE_H_
#define _YOCTO_SHAPE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct shape_data {
  // element data
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};
};

// Shape creation
template <typename PFunc, typename TFunc>
inline shape_data make_lines(int steps, PFunc&& position, TFunc&& tangent);
template <typename PFunc, typename NFunc>
inline shape_data make_quads(int steps, PFunc&& position, NFunc&& normal);

// Interpolate vertex data
vec3f eval_position(const shape_data& shape, int element, vec2f uv);
vec3f eval_normal(const shape_data& shape, int element, vec2f uv);
vec3f eval_tangent(const shape_data& shape, int element, vec2f uv);
vec2f eval_texcoord(const shape_data& shape, int element, vec2f uv);
vec4f eval_color(const shape_data& shape, int element, vec2f uv);
float eval_radius(const shape_data& shape, int element, vec2f uv);

// Evaluate element normals
vec3f eval_element_normal(const shape_data& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const shape_data& shape);
void          compute_normals(vector<vec3f>& normals, const shape_data& shape);

// An unevaluated location on a shape
struct shape_point {
  int   element = 0;
  vec2f uv      = {0, 0};
};

// Shape sampling
vector<float> sample_shape_cdf(const shape_data& shape);
shape_point   sample_shape(
      const shape_data& shape, const vector<float>& cdf, float rn, vec2f ruv);
vector<shape_point> sample_shape(
    const shape_data& shape, int num_samples, uint64_t seed = 98729387);

// Conversions
shape_data quads_to_triangles(const shape_data& shape);
void       quads_to_triangles_inplace(shape_data& shape);

// Subdivision
shape_data subdivide_shape(
    const shape_data& shape, int subdivisions, bool catmullclark);

// Transform shape
shape_data transform_shape(
    const shape_data& shape, const frame3f& frame, bool non_rigid = false);
shape_data transform_shape(const shape_data& shape, const frame3f& frame,
    float radius_scale, bool non_rigid = false);
shape_data scale_shape(const shape_data& shape, float scale, float uvscale = 1);
shape_data scale_shape(shape_data&& shape, float scale, float uvscale = 1);
shape_data flipyz_shape(const shape_data& shape);

// Manipulate vertex data
shape_data remove_normals(const shape_data& shape);
shape_data add_normals(const shape_data& shape);

// Merge a shape into another
void merge_shape_inplace(shape_data& shape, const shape_data& merge);

// Shape statistics
vector<string> shape_stats(const shape_data& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FACE-VARYING SHAPE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Shape data stored as a face-varying mesh
struct fvshape_data {
  // element data
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
};

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, vec2f uv);
vec3f eval_normal(const fvshape_data& shape, int element, vec2f uv);
vec2f eval_texcoord(const shape_data& shape, int element, vec2f uv);

// Evaluate element normals
vec3f eval_element_normal(const fvshape_data& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const fvshape_data& shape);
void compute_normals(vector<vec3f>& normals, const fvshape_data& shape);

// Conversions
shape_data fvshape_to_shape(
    const fvshape_data& shape, bool as_triangles = false);
fvshape_data shape_to_fvshape(const shape_data& shape);

// Subdivision
fvshape_data subdivide_fvshape(
    const fvshape_data& shape, int subdivisions, bool catmullclark);

// Transform shape
fvshape_data transform_fvshape(
    const fvshape_data& shape, const frame3f& frame, bool non_rigid = false);
fvshape_data scale_fvshape(
    const fvshape_data& shape, float scale, float uvscale = 1);
fvshape_data scale_fvshape(
    fvshape_data&& shape, float scale, float uvscale = 1);

// Vertex properties
fvshape_data remove_normals(const fvshape_data& shape);
fvshape_data add_normals(const fvshape_data& shape);

// Shape statistics
vector<string> fvshape_stats(const fvshape_data& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a plane.
shape_data make_rect(
    vec2i steps = {1, 1}, vec2f scale = {1, 1}, vec2f uvscale = {1, 1});
shape_data make_bulged_rect(vec2i steps = {1, 1}, vec2f scale = {1, 1},
    vec2f uvscale = {1, 1}, float radius = 0.3f);
// Make a plane in the xz plane.
shape_data make_recty(
    vec2i steps = {1, 1}, vec2f scale = {1, 1}, vec2f uvscale = {1, 1});
shape_data make_bulged_recty(vec2i steps = {1, 1}, vec2f scale = {1, 1},
    vec2f uvscale = {1, 1}, float radius = 0.3f);
// Make a box.
shape_data make_box(vec3i steps = {1, 1, 1}, vec3f scale = {1, 1, 1},
    vec3f uvscale = {1, 1, 1});
shape_data make_rounded_box(vec3i steps = {1, 1, 1}, vec3f scale = {1, 1, 1},
    vec3f uvscale = {1, 1, 1}, float radius = 0.3f);
// Make a quad stack
shape_data make_rect_stack(
    vec3i steps = {1, 1, 1}, vec3f scale = {1, 1, 1}, vec2f uvscale = {1, 1});
// Make a floor.
shape_data make_floor(
    vec2i steps = {1, 1}, vec2f scale = {10, 10}, vec2f uvscale = {10, 10});
shape_data make_bent_floor(vec2i steps = {1, 1}, vec2f scale = {10, 10},
    vec2f uvscale = {10, 10}, float bent = 0.5f);
// Make a sphere.
shape_data make_sphere(int steps = 32, float scale = 1, float uvscale = 1);
// Make a sphere.
shape_data make_uvsphere(
    vec2i steps = {32, 32}, float scale = 1, vec2f uvscale = {1, 1});
shape_data make_uvspherey(
    vec2i steps = {32, 32}, float scale = 1, vec2f uvscale = {1, 1});
// Make a sphere with slipped caps.
shape_data make_capped_uvsphere(vec2i steps = {32, 32}, float scale = 1,
    vec2f uvscale = {1, 1}, float height = 0.3f);
shape_data make_capped_uvspherey(vec2i steps = {32, 32}, float scale = 1,
    vec2f uvscale = {1, 1}, float height = 0.3f);
// Make a disk
shape_data make_disk(int steps = 32, float scale = 1, float uvscale = 1);
// Make a bulged disk
shape_data make_bulged_disk(
    int steps = 32, float scale = 1, float uvscale = 1, float height = 0.3f);
// Make a uv disk
shape_data make_uvdisk(
    vec2i steps = {32, 32}, float scale = 1, vec2f uvscale = {1, 1});
// Make a uv cylinder
shape_data make_uvcylinder(vec3i steps = {32, 32, 32}, vec2f scale = {1, 1},
    vec3f uvscale = {1, 1, 1});
// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(vec3i steps = {32, 32, 32},
    vec2f scale = {1, 1}, vec3f uvscale = {1, 1, 1}, float radius = 0.3f);
// Make a uv capsule
shape_data make_uvcapsule(vec3i steps = {32, 32, 32}, vec2f scale = {1, 1},
    vec3f uvscale = {1, 1, 1});
// Make a uv cone
shape_data make_uvcone(vec3i steps = {32, 32, 32}, vec2f scale = {1, 1},
    vec3f uvscale = {1, 1, 1});

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
shape_data make_lines(int num = 65536, int steps = 4, vec2f scale = {1, 1},
    vec2f uvscale = {1, 1}, vec2f radius = {0.001f, 0.001f});

// Make a point primitive. Returns points, pos, norm, texcoord, radius.
shape_data make_point(float radius = 0.001f);
// Make a point set on a grid. Returns points, pos, norm, texcoord, radius.
shape_data make_points(
    int num = 65536, float uvscale = 1, float radius = 0.001f);
shape_data make_points(vec2i steps = {256, 256}, vec2f size = {1, 1},
    vec2f uvscale = {1, 1}, vec2f radius = {0.001f, 0.001f});
// Make random points in a cube. Returns points, pos, norm, texcoord, radius.
shape_data make_random_points(int num = 65536, vec3f size = {1, 1, 1},
    float uvscale = 1, float radius = 0.001f, uint64_t seed = 17);

// Predefined meshes
shape_data   make_monkey(int subdivisions = 0, float scale = 1);
shape_data   make_quad(int subdivisions = 0, float scale = 1);
shape_data   make_quady(int subdivisions = 0, float scale = 1);
shape_data   make_cube(int subdivisions = 0, float scale = 1);
fvshape_data make_fvcube(int subdivisions = 0, float scale = 1);
shape_data   make_geosphere(int subdivisions = 0, float scale = 1);

// Make a hair ball around a shape.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (strength/number)
// rotation: rotation added to hair (angle/strength)
shape_data make_hair(const shape_data& shape, vec2i steps = {8, 65536},
    vec2f length = {0.1f, 0.1f}, vec2f radius = {0.001f, 0.001f},
    vec2f noise = {0, 10}, vec2f clump = {0, 128}, vec2f rotation = {0, 0},
    int seed = 7);

// Grow hairs around a shape
shape_data make_hair2(const shape_data& shape, vec2i steps = {8, 65536},
    vec2f length = {0.1f, 0.1f}, vec2f radius = {0.001f, 0.001f},
    float noise = 0, float gravity = 0.001f, int seed = 7);

// Grow hairs around a shape
shape_data make_random_hairs(const shape_data& shape, int num = 65536,
    int steps = 8, vec2f length = {1, 1}, vec2f radius = {0.01, 0.01},
    float noise = 0, float gravity = 0.05, uint64_t seed = 7);

// Grow points around a shape
shape_data make_random_points(const shape_data& shape, int num = 65536,
    float radius = 0.01f, uint64_t seed = 7);

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res.
shape_data points_to_spheres(
    const vector<vec3f>& vertices, int steps = 2, float scale = 0.01f);
shape_data polyline_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);
shape_data lines_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);
shape_data lines_to_cylinders(const vector<vec2i>& lines,
    const vector<vec3f>& positions, int steps = 4, float scale = 0.01f);

// Make a heightfield mesh.
shape_data make_heightfield(vec2i size, const vector<float>& height);
shape_data make_heightfield(vec2i size, const vector<vec4f>& color);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// SHAPE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Shape creation
template <typename PFunc, typename TFunc>
inline shape_data make_lines(int steps, PFunc&& position, TFunc&& tangent) {
  auto shape      = shape_data{};
  shape.positions = vector<vec3f>(steps + 1);
  shape.normals   = vector<vec3f>(steps + 1);
  shape.texcoords = vector<vec2f>(steps + 1);
  for (auto idx : range(steps + 1)) {
    auto u               = (float)idx / (float)steps;
    shape.positions[idx] = position(u);
    shape.normals[idx]   = tangent(u);
    shape.texcoords[idx] = {u, 0};
  }
  shape.lines = vector<vec2i>(steps);
  for (auto idx : range(steps)) shape.lines[idx] = {idx, idx + 1};
  return shape;
}
template <typename PFunc, typename NFunc>
inline shape_data make_quads(vec2i steps, PFunc&& position, NFunc&& normal) {
  auto shape      = shape_data{};
  shape.positions = vector<vec3f>((steps.x + 1) * (steps.y + 1));
  shape.normals   = vector<vec3f>((steps.x + 1) * (steps.y + 1));
  shape.texcoords = vector<vec2f>((steps.x + 1) * (steps.y + 1));
  for (auto j : range(steps.y + 1)) {
    for (auto i : range(steps.x + 1)) {
      auto uv              = vec2f{i / (float)steps.x, j / (float)steps.y};
      auto idx             = j * (steps.x + 1) + i;
      shape.positions[idx] = position(uv);
      shape.normals[idx]   = normal(uv);
      shape.texcoords[idx] = {uv.x, 1 - uv.y};
    }
  }
  shape.quads = vector<vec4i>(steps.x * steps.y);
  for (auto j : range(steps.y)) {
    for (auto i : range(steps.x)) {
      auto idx         = j * steps.x + i;
      shape.quads[idx] = {j * (steps.x + 1) + i, j * (steps.x + 1) + i + 1,
          (j + 1) * (steps.x + 1) + i + 1, (j + 1) * (steps.x + 1) + i};
    }
  }
  return shape;
}

}  // namespace yocto

#endif
