//
// # Yocto/Modeling: Curve and surface operations
//
// Yocto/Modeling provides basic utilities for working with smooth curves and
// surfaces.
//

//
// LICENSE:
//
// Copyright (c) 2024 -- 2024 Fabio Pellacini
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

#ifndef _YOCTO_MODELING_H_
#define _YOCTO_MODELING_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::tuple;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// BÉZIER SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate Bernstein polynomials (compile-time version)
template <int Index, int Degree>
inline float eval_bezier_basis(float u);

// Evaluate a Bézier curve at a given u (compile-time version)
template <typename T, size_t N>
inline T eval_bezier(const array<T, N>& control_points, float u);

// Evaluate a Bézier tangent at a given u (compile-time version)
template <typename T, size_t N>
inline T eval_bezier_tangent(const array<T, N>& control_points, float u);

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T, size_t N>
inline vector<T> tesselate_bezier(
    const array<T, N>& control_points, int steps = 100);

// Evaluate Bernstein polynomials
inline float eval_bezier_basis(float u, int index, int degree);

// Evaluate Bézier curve made of multiple segments
template <typename T>
inline T eval_bezier(
    const vector<T>& control_points, int degree, bool closed, float u);

// Evaluate Bézier curve tangents for a curve
template <typename T>
inline T eval_bezier_tangent(
    const vector<T>& control_points, int degree, bool closed, float u);

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T>
inline vector<T> tesselate_bezier(
    const vector<T>& control_points, int degree, bool closed, int steps = 100);

// Subdivide a Bézier curve made of multiple segments using De Casteljau
template <typename T>
inline vector<vector<T>> subdivide_bezier(
    const vector<T>& control_points, int levels = 6);

// Evaluate a Bézier curve
template <typename T>
inline T eval_bezier(const vector<T>& control_points, float u);

// Evaluate Bézier curve tangents for a curve
template <typename T>
inline T eval_bezier_tangent(const vector<T>& control_points, float u);

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T>
inline vector<T> tesselate_bezier(
    const vector<T>& control_points, int steps = 100);

// Flatten a set of Bezier control points to a polyline
template <typename T>
inline vector<T> flatten_beziers(const vector<vector<T>>& beziers);

// Convert a long Bézier spline to Bézier segments
template <typename T>
inline vector<vector<T>> bezier_to_beziers(
    const vector<T>& bezier, int degree, bool closed);

// Subdivide a Bézier curve made of multiple segments using De Casteljau
template <typename T>
inline vector<vector<T>> subdivide_beziers(
    const vector<vector<T>>& control_points, int levels = 6);

}  // namespace yocto

// -----------------------------------------------------------------------------
// CATMULL-ROM SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate a Catmull-Rom spline
template <typename T>
inline T eval_catmullrom(
    const vector<T>& control_points, int degree, bool closed, float u);

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_catmullrom(
    const vector<T>& control_points, int degree, bool closed, int steps = 100);

// Convert a Catmull-Rom spline to Bézier segments
template <typename T>
inline vector<array<T, 4>> catmullrom_to_beziers(
    const vector<T>& catmullrom, int degree, bool closed);

}  // namespace yocto

// -----------------------------------------------------------------------------
// B-SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate the bspline basis up to degree 3 (compile-time version)
template <int Degree>
inline float eval_bspline_basis(float u);

// Evaluate B-spline by convolution
template <typename T, size_t N>
inline T eval_bspline(const array<T, N>& control_points, bool closed, float u);

// Evaluate the bspline basis up to degree 3
inline float eval_bspline_basis(float u, int degree);

// Evaluate B-spline by convolution
template <typename T>
inline T eval_bspline(
    const vector<T>& control_points, int degree, bool closed, float u);

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_bspline(
    const vector<T>& control_points, int degree, bool closed, int steps = 100);

// Evaluate B-spline by convolution
template <typename T>
inline T eval_bspline(const vector<T>& control_points, bool closed, float u);

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_bspline(
    const vector<T>& control_points, bool closed, int steps = 100);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BÉZIER PATCHES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate Bernstein polynomials
inline float eval_patch_basis(vec2f uv, vec2i index, vec2i degree);

// Evaluate Bézier patch made of multiple segments.
// For now, we only support one patch.
template <typename T>
inline T eval_patch(const vector<vector<T>>& control_points, vec2i degree,
    bool closed, vec2f uv);

// Extract an isoline from a Bézier patch made of multiple segments.
// For now, we only support one patch.
template <typename T>
inline vector<T> patch_isoline(const vector<vector<T>>& control_points,
    vec2i degree, bool closed, float param, bool u_isoline);

// Tesselate a Bézier patch made of multiple segments using uniform sampling.
// For now, we only support one patch.
template <typename T>
inline pair<vector<vec4i>, vector<T>> tesselate_patch(
    const vector<vector<T>>& control_points, vec2i degree, bool closed,
    vec2i steps = {100, 100});

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION CURVES
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines.
template <typename T>
static pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vertices,
    int subdivisions = 1);

// Check if it is a boundary vertex
inline vector<int> get_boundary_vertices(
    int npoints, const vector<vec2i>& lines);

// B-spline subdivision
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vert, int subdivisions = 1);
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vert,
    const vector<int>& corners, int subdivisions = 1);

// Subdivide beziers.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vertices,
    int subdivisions = 1);

}  // namespace yocto

// -----------------------------------------------------------------------------
// MESH EDGES AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Binary search for edge index
inline int search_edge(const vector<vec2i>& edges, vec2i edge);

// Get edges of meshes
inline vector<vec2i> triangles_edges(const vector<vec3i>& triangles);
inline vector<vec2i> quads_edges(const vector<vec4i>& quads);
inline vector<vec2i> faces_edges(
    const vector<vec3i>& triangles, const vector<vec4i>& quads);

// Get boundaries of meshes
inline vector<vec2i> triangles_boundaries(const vector<vec3i>& triangles);
inline vector<vec2i> quads_boundaries(const vector<vec4i>& quads);
inline vector<vec2i> faces_boundaries(
    const vector<vec3i>& triangles, const vector<vec4i>& quads);

// Get edges and boundaries of meshes
inline pair<vector<vec2i>, vector<vec2i>> triangles_edges_boundaries(
    const vector<vec3i>& triangles);
inline pair<vector<vec2i>, vector<vec2i>> quads_edges_boundaries(
    const vector<vec4i>& quads);
inline pair<vector<vec2i>, vector<vec2i>> faces_edges_boundaries(
    const vector<vec3i>& triangles, const vector<vec4i>& quads);

// Build adjacencies between faces (sorted counter-clockwise)
inline vector<vec3i> face_adjacencies(const vector<vec3i>& triangles);

// Build adjacencies between vertices (sorted counter-clockwise)
inline vector<vector<int>> vertex_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies);

// Compute boundaries as a list of loops (sorted counter-clockwise)
inline vector<vector<int>> ordered_boundaries(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int num_vertices);

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
inline vector<vector<int>> vertex_to_faces_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION SURFACES
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide triangles.
template <typename T>
static pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vertices,
    int subdivisions = 1);

// Subdivide quads.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vertices,
    int subdivisions = 1);

// Subdivide Catmull-Clark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vertices, int subdivisions = 1,
    bool lock_boundary = false);

// Subdivide Catmull-Clark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_creased(
    const vector<vec4i>& quads, const vector<T>& vertices,
    const vector<vec2i>& creases, const vector<int>& corners,
    int subdivisions = 1, bool lock_boundary = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DISPLACEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// Displace vertices
inline vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const image<float>& displacement, float scale = 1.0f, float offset = 0.5f);
inline vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const image<vec4f>& displacement, float scale = 1.0f, float offset = 0.5f);

// Shape displacement
inline shape_data displace_shape(const shape_data& shape,
    const image<float>& displacement, float height = 1, float offset = 0.5f);
inline shape_data displace_shape(const shape_data& shape,
    const image<vec4f>& displacement, float height = 1, float offset = 0.5f);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIGNED-DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Here are some examples of analytic SDFs. You can find more at
// https://iquilezles.org/articles/distfunctions/

// Sphere SDF
inline float sdf_sphere(vec3f position, vec3f center, float radius);
// Box SDF
inline float sdf_box(vec3f position, vec3f center, vec3f size);
// Capsule SDF
inline float sdf_capsule(vec3f position, vec3f a, vec3f b, float radius);

// Union SDF
inline float sdf_union(float sdf1, float sdf2);
// Intersection SDF
inline float sdf_intersection(float sdf1, float sdf2);
// Difference SDF
inline float sdf_difference(float sdf1, float sdf2);

// Soft union SDF
inline float sdf_smooth_union(float sdf1, float sdf2, float k = 0.1f);
// Soft intersection SDF
inline float sdf_smooth_intersection(float sdf1, float sdf2, float k = 0.1f);
// Soft difference SDF
inline float sdf_smooth_difference(float sdf1, float sdf2, float k = 0.1f);

}  // namespace yocto

// -----------------------------------------------------------------------------
// NOISE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// the implementation follows ideas from "Building up Perlin Noise"
// http://eastfarthing.com/blog/2015-04-21-noise/

// Value noise
inline float value_noise(vec2f position);
inline float value_noise(vec3f position);

// Gradient noise
inline float gradient_noise(vec2f position);
inline float gradient_noise(vec3f position);

// Cell noise
inline float cell_noise(vec2f position, bool inverted = false);
inline float cell_noise(vec3f position, bool inverted = false);

// Voronoi-like noise
inline float voronoi_noise(vec2f position);
inline float voronoi_noise(vec3f position);

// Fractal noise
template <typename Noise>
inline float fractal_noise(Noise&& base_noise, vec2f position, int num = 8);
// Turbulence noise
template <typename Noise>
inline float turbulence_noise(Noise&& base_noise, vec2f position, int num = 8);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTOR HASHING
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::vec2i> {
  size_t operator()(const yocto::vec2i& v) const {
    static const auto hasher = std::hash<int>();
    auto              h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};
template <>
struct hash<yocto::vec3i> {
  size_t operator()(const yocto::vec3i& v) const {
    static const auto hasher = std::hash<int>();
    auto              h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};
template <>
struct hash<yocto::vec4i> {
  size_t operator()(const yocto::vec4i& v) const {
    static const auto hasher = std::hash<int>();
    auto              h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.w) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------
namespace yocto {

// A sparse grid of cells, containing list of points. Cells are stored in
// a dictionary to get sparsity. Helpful for nearest neighbor lookups.
struct hash_grid {
  float                             cell_size     = 0;
  float                             cell_inv_size = 0;
  vector<vec3f>                     positions     = {};
  unordered_map<vec3i, vector<int>> cells         = {};
};

// Create a hash_grid
inline hash_grid make_hash_grid(float cell_size);
inline hash_grid make_hash_grid(
    const vector<vec3f>& positions, float cell_size);
// Inserts a point into the grid
inline int insert_vertex(hash_grid& grid, vec3f position);
// Finds the nearest neighbors within a given radius
inline void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    vec3f position, float max_radius);
inline void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    int vertex, float max_radius);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// ARRAY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Concatenation
template <typename T>
inline vector<T> join(const vector<T>& a, const vector<T>& b) {
  auto c = vector<T>{};
  c.reserve(a.size() + b.size());
  c.insert(c.end(), a.begin(), a.end());
  c.insert(c.end(), b.begin(), b.end());
  return c;
}

// Concatenation
template <typename T>
inline vector<T>& append(vector<T>& a, const vector<T>& b) {
  a.reserve(a.size() + b.size());
  a.insert(a.end(), b.begin(), b.end());
  return a;
}

// Remove duplicates
template <typename T>
inline vector<T> remove_duplicates(const vector<T>& v) {
  auto result = v;
  std::sort(result.begin(), result.end());
  result.erase(std::unique(result.begin(), result.end()), result.end());
  return result;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BÉZIER SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate Bernstein polynomials (compile-time version)
template <int Index, int Degree>
inline float eval_bezier_basis(float u) {
  // Check ranges
  static_assert(
      Index < 0 || Index > Degree || Degree < 0, "index out of range");

  // Implementation using recursion
  if constexpr (Degree == 0) {
    return 1;
  } else if constexpr (Degree == 1) {
    if constexpr (Index == 0) return 1 - u;
    if constexpr (Index == 1) return u;
  } else if constexpr (Degree == 2) {
    if constexpr (Index == 0) return (1 - u) * (1 - u);
    if constexpr (Index == 1) return 2 * u * (1 - u);
    if constexpr (Index == 2) return u * u;
  } else if constexpr (Degree == 3) {
    if constexpr (Index == 0) return (1 - u) * (1 - u) * (1 - u);
    if constexpr (Index == 1) return 3 * u * (1 - u) * (1 - u);
    if constexpr (Index == 2) return 3 * u * u * (1 - u);
    if constexpr (Index == 3) return u * u * u;
  } else {
    // Recursive definition of Bernstein polynomials
    return eval_bezier_basis<Index - 1, Degree - 1>(u) * u +
           eval_bezier_basis<Index, Degree - 1>(u) * (1 - u);
  }
}

// Evaluate a Bézier curve at a given u (compile-time version)
template <typename T, size_t N>
inline T eval_bezier(const array<T, N>& control_points, float u) {
  // get degree
  constexpr int Degree = (int)N - 1;

  // check supported degree
  static_assert(Degree >= 0 && Degree <= 3, "unsupported degree");

  // evaluate bernstein polynomials directly
  if constexpr (Degree == 0) {
    return control_points[0];
  } else if constexpr (Degree == 1) {
    return eval_bezier_basis<0, Degree>(u) * control_points[0] +
           eval_bezier_basis<1, Degree>(u) * control_points[1];
  } else if constexpr (Degree == 2) {
    return eval_bezier_basis<0, Degree>(u) * control_points[0] +
           eval_bezier_basis<1, Degree>(u) * control_points[1] +
           eval_bezier_basis<2, Degree>(u) * control_points[2];
  } else if constexpr (Degree == 3) {
    return eval_bezier_basis<0, Degree>(u) * control_points[0] +
           eval_bezier_basis<1, Degree>(u) * control_points[1] +
           eval_bezier_basis<2, Degree>(u) * control_points[2] +
           eval_bezier_basis<3, Degree>(u) * control_points[3];
  } else {
    // use recursive definition of Bézier curves
    // TODO: recurve definition since we do not have a constexpr for loop
    return {};
  }
}

// Evaluate a Bézier tangent at a given u (compile-time version)
template <typename T, size_t N>
inline T eval_bezier_tangent(const array<T, N>& control_points, float u) {
  // get degree
  constexpr int Degree = (int)N - 1;

  // check supported degree
  static_assert(Degree >= 0 && Degree <= 3, "unsupported degree");

  // evaluate bernstein polynomials directly
  if constexpr (Degree == 0) {
    return T{};
  } else if constexpr (Degree == 1) {
    return Degree * eval_bezier_basis<0, Degree - 1>(u) *
           (control_points[1] - control_points[0]);
  } else if constexpr (Degree == 2) {
    return Degree * eval_bezier_basis<0, Degree - 1>(u) *
               (control_points[1] - control_points[0]) +
           Degree * eval_bezier_basis<1, Degree - 1>(u) *
               (control_points[2] - control_points[1]);
  } else if constexpr (Degree == 3) {
    return Degree * eval_bezier_basis<0, Degree - 1>(u) *
               (control_points[1] - control_points[0]) +
           Degree * eval_bezier_basis<1, Degree - 1>(u) *
               (control_points[2] - control_points[1]) +
           Degree * eval_bezier_basis<1, Degree - 1>(u) *
               (control_points[3] - control_points[2]);
  } else {
    // use recursive definition of Bézier curves
    // TODO: recurve definition since we do not have a constexpr for loop
    return {};
  }
}

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T, size_t N>
inline vector<T> tesselate_bezier(
    const array<T, N>& control_points, int steps) {
  // tesselate each segment
  auto tessellated_points = vector<T>();
  for (auto step : range(steps + 1)) {
    auto u = (float)step / (float)steps;
    tessellated_points.push_back(eval_bezier(control_points, u));
  }

  // done
  return tessellated_points;
}

// Evaluate Bernstein polynomials
inline float eval_bezier_basis(float u, int index, int degree) {
  // Implementation using recursion
  if (index < 0 || index > degree || degree < 0) return 0;
  if (degree == 0 && index == 0) return 1;
  if (degree == 1) {
    if (index == 0) return 1 - u;
    if (index == 1) return u;
    return 0;
  }
  if (degree == 2) {
    if (index == 0) return (1 - u) * (1 - u);
    if (index == 1) return 2 * u * (1 - u);
    if (index == 2) return u * u;
    return 0;
  }
  if (degree == 3) {
    if (index == 0) return (1 - u) * (1 - u) * (1 - u);
    if (index == 1) return 3 * u * (1 - u) * (1 - u);
    if (index == 2) return 3 * u * u * (1 - u);
    if (index == 3) return u * u * u;
    return 0;
  }
  return eval_bezier_basis(u, index - 1, degree - 1) * u +
         eval_bezier_basis(u, index, degree - 1) * (1 - u);
}

// Evaluate Bézier curve made of multiple segments
template <typename T>
inline T eval_bezier(
    const vector<T>& control_points, int degree, bool closed, float u) {
  // get number of segments
  auto segments = ((int)control_points.size() - (closed ? 0 : 1)) / degree;

  // get segment and number of points
  auto segment = clamp((int)u, 0, segments - 1);
  auto npoints = (int)control_points.size();

  // compute local u
  auto usegment = u - segment;

  // evaluate bernstein polynomials
  auto point = T();
  for (auto index : range(degree + 1)) {
    point += eval_bezier_basis(usegment, index, degree) *
             control_points[(segment * degree + index) % npoints];
  }
  return point;
}

// Evaluate Bézier curve tangents for a curve
template <typename T>
inline T eval_bezier_tangent(
    const vector<T>& control_points, int degree, bool closed, float u) {
  // get number of segments
  auto segments = ((int)control_points.size() - (closed ? 0 : 1)) / degree;

  // get segment and number of points
  auto segment = clamp((int)u, 0, segments - 1);
  auto npoints = (int)control_points.size();

  // compute local u
  auto usegment = u - segment;

  // evaluate bernstein polynomials
  auto tangent = T();
  for (auto index : range(degree)) {
    tangent += eval_bezier_basis(usegment, index, degree - 1) *
               (control_points[(segment * degree + index + 1) % npoints] -
                   control_points[(segment * degree + index) % npoints]);
  }
  tangent = degree * tangent;
  return tangent;
}

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T>
inline vector<T> tesselate_bezier(
    const vector<T>& control_points, int degree, bool closed, int steps) {
  // get number of segments
  auto segments = ((int)control_points.size() - (closed ? 0 : 1)) / degree;

  // tesselate each segment with n steps
  auto tessellated_points = vector<T>();
  for (auto segment : range(segments)) {
    auto npoints = (int)control_points.size();
    auto points  = vector<T>(degree + 1);
    for (auto index : range(degree + 1)) {
      points[index] = control_points[(segment * degree + index) % npoints];
    }
    for (auto step : range(steps + 1)) {
      auto u = segment + (float)step / (float)steps;
      tessellated_points.push_back(
          eval_bezier(control_points, degree, closed, u));
    }
  }
  return tessellated_points;
}

// Subdivide a Bézier curve made of multiple segments using De Casteljau
template <typename T>
inline vector<vector<T>> subdivide_bezier(
    const vector<T>& control_points, int levels) {
  // get degree
  auto degree = (int)control_points.size() - 1;

  // allocate buffers
  auto points = vector<vector<vec3f>>(degree + 1);
  for (auto sdegree : range(degree + 1)) points[sdegree].resize(sdegree + 1);

  // tesselate each segment recursively
  auto tessellated_points = vector<vector<T>>{control_points};
  for (auto level : range(levels)) {
    // get number of segments
    auto segments = (int)tessellated_points.size();

    // subdivide each segment
    auto subdivided_points = vector<vector<T>>{};
    for (auto segment : range(segments)) {
      // copy control points
      points[degree] = tessellated_points[segment];

      // evaluate midpoints points
      for (auto sdegree = degree - 1; sdegree >= 0; sdegree--) {
        for (auto index : range(sdegree + 1)) {
          points[sdegree][index] = 0.5f * points[sdegree + 1][index + 0] +
                                   0.5f * points[sdegree + 1][index + 1];
        }
      }

      // copy vertices
      subdivided_points.push_back({});
      for (auto sdegree : range(degree + 1)) {
        subdivided_points.back().push_back(points[degree - sdegree].front());
      }
      subdivided_points.push_back({});
      for (auto sdegree : range(degree + 1)) {
        subdivided_points.back().push_back(points[sdegree].back());
      }
    }

    // swap for next level
    tessellated_points = subdivided_points;
  }

  // done
  return tessellated_points;
}

// Evaluate a Bézier curve
template <typename T>
inline T eval_bezier(const vector<T>& control_points, float u) {
  return eval_bezier(control_points, (int)control_points.size() - 1, false, u);
}

// Evaluate Bézier curve tangents for a curve
template <typename T>
inline T eval_bezier_tangent(const vector<T>& control_points, float u) {
  return eval_bezier_tangent(
      control_points, (int)control_points.size() - 1, false, u);
}

// Tesselate a Bézier curve made of multiple segments using uniform sampling.
template <typename T>
inline vector<T> tesselate_bezier(const vector<T>& control_points, int steps) {
  return tesselate_bezier(
      control_points, (int)control_points.size() - 1, false, steps);
}

// Flatten a set of Bezier control points to a polyline
template <typename T>
inline vector<T> flatten_beziers(const vector<vector<T>>& beziers) {
  auto points = vector<T>{};
  for (auto& bezier : beziers) append(points, bezier);
  return points;
}

// Convert a long Bézier spline to Bézier segments
template <typename T>
inline vector<vector<T>> bezier_to_beziers(
    const vector<T>& bezier, int degree, bool closed) {
  // check for minimal length
  if (bezier.size() < degree + 1)
    throw std::runtime_error("spline is too short");

  // number of segments
  auto segments = ((int)bezier.size() - (closed ? 0 : 1)) / degree;
  auto npoints  = (int)bezier.size();

  // compute bezier control points
  auto beziers = vector<vector<T>>(segments);
  for (auto segment : range(segments)) {
    beziers[segment].resize(degree + 1);
    for (auto index : range(degree + 1)) {
      beziers[segment][index] = bezier[(segment * degree + index) % npoints];
    }
  }

  // done
  return beziers;
}

// Subdivide a Bézier curve made of multiple segments using De Casteljau
template <typename T>
inline vector<vector<T>> subdivide_beziers(
    const vector<vector<T>>& control_points, int levels) {
  // tesselate each segment
  auto subdivided_beziers = vector<vector<T>>();
  for (auto& control_points_ : control_points) {
    append(subdivided_beziers, subdivide_bezier(control_points_, levels));
  }
  return subdivided_beziers;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CATMULL-ROM SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate a Catmull-Rom spline
template <typename T>
inline T eval_catmullrom(
    const vector<T>& control_points, int degree, bool closed, float u) {
  // check for degree
  if (degree != 3) throw std::runtime_error("only cubic Catmull-Rom supported");

  // check for minimal length
  if (control_points.size() < 4)
    throw std::runtime_error("spline is too short");

  // get number of segments
  auto segments = (int)control_points.size() - (closed ? 0 : 1);

  // get segment and number of points
  auto segment = clamp((int)u, 0, segments - 1);
  auto npoints = (int)control_points.size();

  // save control points for convenience
  auto &p1 = control_points[segment],
       &p2 = control_points[(segment + 1) % npoints];
  auto &p0 = control_points[closed ? (segment - 1 + npoints) % npoints
                                   : max(segment - 1, 0)],
       &p3 = control_points[closed ? (segment + 2) % npoints
                                   : min(segment + 2, npoints - 1)];

  // construct Bezier segment
  auto bezier = array<T, 4>{
      p1,
      p1 + 1.0f / 6.0f * (p2 - p0),
      p2 - 1.0f / 6.0f * (p3 - p1),
      p2,
  };

  // evaluate bezier
  auto usegment = u - segment;
  auto point    = T();
  for (auto index : range(4)) {
    point += eval_bezier_basis(usegment, index, 3) * bezier[index];
  }
  return point;
}

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_catmullrom(
    const vector<T>& control_points, int degree, bool closed, int steps) {
  // get number of segments
  auto segments = (int)control_points.size() - (closed ? 0 : 1);

  // tesselate each segment with n steps
  auto tessellated_points = vector<T>();
  for (auto segment : range(segments)) {
    for (auto step : range(steps + 1)) {
      auto u = segment + (float)step / (float)steps;
      tessellated_points.push_back(
          eval_catmullrom(control_points, degree, closed, u));
    }
  }
  return tessellated_points;
}

// Convert a Catmull-Rom spline to Bézier segments
template <typename T>
inline vector<array<T, 4>> catmullrom_to_beziers(
    const vector<T>& catmullrom, int degree, bool closed) {
  // check for degree
  if (degree != 3) throw std::runtime_error("only cubic Catmull-Rom supported");

  // check for minimal length
  if (catmullrom.size() < 4) throw std::runtime_error("spline is too short");

  // number of segments
  auto segments = (int)catmullrom.size() - (closed ? 0 : 1);
  auto npoints  = (int)catmullrom.size();

  // compute bezier control points
  auto beziers = vector<array<T, 4>>(segments);
  for (auto index : range(segments)) {
    // save points for convenience
    auto &p1 = catmullrom[index], &p2 = catmullrom[(index + 1) % npoints];
    auto &p0 = catmullrom[closed ? (index - 1 + npoints) % npoints
                                 : max(index - 1, 0)],
         &p3 = catmullrom[closed ? (index + 2) % npoints
                                 : min(index + 2, npoints - 1)];
    // end points
    beziers[index][0] = p1;
    beziers[index][3] = p2;
    // tangents
    beziers[index][1] = p1 + 1.0f / 6.0f * (p2 - p0);
    beziers[index][2] = p2 - 1.0f / 6.0f * (p3 - p1);
  }

  // done
  return beziers;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// B-SPLINES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate the bspline basis up to degree 3 (compile-time version)
template <int Degree>
inline float eval_bspline_basis(float u) {
  // check degree
  static_assert(Degree >= 0 || Degree <= 3, "degree not supported");

  // get parameter
  auto abs_u = abs(u);

  if constexpr (Degree == 0) {
    if (abs_u < 0.5f) return 1;
    return 0;
  } else if constexpr (Degree == 1) {
    if (abs_u < 1.0f) return 1 - abs_u;
    return 0;
  } else if constexpr (Degree == 2) {
    if (abs_u < 1.0f) return 0.5f * abs_u * abs_u;
    return 0;
  } else if constexpr (Degree == 3) {
    if (abs_u < 1.0f) return 2.0f / 3.0f - 0.5f * abs_u * abs_u * (2 - abs_u);
    if (abs_u < 2.0f)
      return 1.0f / 6.0f * (2 - abs_u) * (2 - abs_u) * (2 - abs_u);
    return 0;
  } else {
    return 0;
  }
}

// Evaluate B-spline by convolution
template <typename T, size_t N>
inline T eval_bspline(const array<T, N>& control_points, bool closed, float u) {
  // get degree
  constexpr int Degree = (int)N - 1;

  // get segment and number of points
  auto npoints = (int)control_points.size();
  auto segment = clamp((int)u, 0, npoints - (closed ? 0 : 1));

  // accumulate point by discrete-continuous convolution
  auto point = T();
  for (auto index : range(segment - Degree, segment + Degree + 1)) {
    auto pindex = closed ? ((index + npoints) % npoints)
                         : clamp(index, 0, npoints - 1);
    point += eval_bspline_basis<Degree>(u - index) * control_points[pindex];
  }
  return point;
}

// Evaluate the bspline basis up to degree 3
inline float eval_bspline_basis(float u, int degree) {
  // check degree
  if (degree < 0 || degree > 3)
    throw std::runtime_error("degree not supported");

  // get parameter
  auto abs_u = abs(u);

  if (degree == 0) {
    if (abs_u < 0.5f) return 1;
    return 0;
  } else if (degree == 1) {
    if (abs_u < 1.0f) return 1 - abs_u;
    return 0;
  } else if (degree == 2) {
    if (abs_u < 0.5f) return -abs_u * abs_u + 0.75f;
    if (abs_u < 1.5f) return 0.5f * (1.5f - abs_u) * (1.5f - abs_u);
    return 0;
  } else if (degree == 3) {
    if (abs_u < 1.0f) return 2.0f / 3.0f - 0.5f * abs_u * abs_u * (2 - abs_u);
    if (abs_u < 2.0f)
      return 1.0f / 6.0f * (2 - abs_u) * (2 - abs_u) * (2 - abs_u);
    return 0;
  } else {
    throw std::runtime_error("degree not supported");
  }
}

// Evaluate B-spline by convolution
template <typename T>
inline T eval_bspline(
    const vector<T>& control_points, int degree, bool closed, float u) {
  // get segment and number of points
  auto segment = clamp(
      (int)u, 0, (int)control_points.size() - (closed ? 0 : 1));
  auto npoints = (int)control_points.size();

  // accumulate point by discrete-continuous convolution
  auto point = T();
  // for (auto index : range(segment - degree, segment + degree + 1)) {
  //   auto pindex = closed ? ((index + npoints) % npoints)
  //                        : clamp(index, 0, npoints - 1);
  //   point += eval_bspline_basis(u - index, degree) * control_points[pindex];
  // }
  for (auto index : range(segment - degree - 1, segment + degree + 2)) {
    auto pindex = closed ? ((index + npoints) % npoints)
                         : clamp(index, 0, npoints - 1);
    point += eval_bspline_basis(u - index, degree) * control_points[pindex];
  }
  return point;
}

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_bspline(
    const vector<T>& control_points, int degree, bool closed, int steps) {
  // get number of segments
  auto segments = (int)control_points.size() - (closed ? 0 : 1);

  // tesselate each segment with n steps
  auto tessellated_points = vector<T>();
  for (auto segment : range(segments)) {
    for (auto step : range(steps + 1)) {
      auto u = segment + (float)step / (float)steps;
      tessellated_points.push_back(
          eval_bspline(control_points, degree, closed, u));
    }
  }
  return tessellated_points;
}

// Evaluate B-spline by convolution
template <typename T>
inline T eval_bspline(const vector<T>& control_points, bool closed, float u) {
  return eval_bspline(control_points, 3, closed, u);
}

// Tessellate a B-spline curve
template <typename T>
inline vector<T> tesselate_bspline(
    const vector<T>& control_points, bool closed, int steps) {
  return tesselate_bspline(control_points, 3, closed, steps);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BÉZIER PATCHES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate Bernstein polynomials
inline float eval_patch_basis(vec2f uv, vec2i index, vec2i degree) {
  return eval_bezier_basis(uv.x, index.x, degree.x) *
         eval_bezier_basis(uv.y, index.y, degree.y);
}

// Evaluate Bézier patch made of multiple segments.
// For now, we only support one patch.
template <typename T>
inline T eval_patch(const vector<vector<T>>& control_points, vec2i degree,
    bool closed, vec2f uv) {
  // evaluate bernstein polynomials
  auto point = T();
  for (auto index : range(degree + 1)) {
    point += eval_patch_basis(uv, index, degree) *
             control_points[index.x][index.y];
  }
  return point;
}

// Extract an isoline from a Bézier patch made of multiple segments.
// For now, we only support one patch.
template <typename T>
inline vector<T> patch_isoline(const vector<vector<T>>& control_points,
    vec2i degree, bool closed, float param, bool u_isoline) {
  auto bezier_points = vector<T>{};
  if (u_isoline) {
    bezier_points = vector<T>(degree.y + 1, T{});
    for (auto index : range(degree + 1)) {
      bezier_points[index.y] += eval_bezier_basis(param, index.x, degree.x) *
                                control_points[index.x][index.y];
    }
  } else {
    bezier_points = vector<T>(degree.x + 1, T{});
    for (auto index : range(degree + 1)) {
      bezier_points[index.x] += eval_bezier_basis(param, index.y, degree.y) *
                                control_points[index.x][index.y];
    }
  }
  return bezier_points;
}

// Tesselate a Bézier patch made of multiple segments using uniform sampling.
// For now, we only support one patch.
template <typename T>
inline pair<vector<vec4i>, vector<T>> tesselate_patch(
    const vector<vector<T>>& control_points, vec2i degree, bool closed,
    vec2i steps) {
  // tesselate each patch with n x n steps
  auto tessellated_points = vector<T>();
  for (auto step : range(steps + 1)) {
    auto uv = (vec2f)step / (vec2f)steps;
    tessellated_points.push_back(
        eval_patch(control_points, degree, closed, uv));
  }
  auto tessellated_quads = vector<vec4i>();
  for (auto [i, j] : range(steps)) {
    tessellated_quads.push_back({
        (j + 0) * (steps.x + 1) + (i + 0),
        (j + 0) * (steps.x + 1) + (i + 1),
        (j + 1) * (steps.x + 1) + (i + 1),
        (j + 1) * (steps.x + 1) + (i + 0),
    });
  }
  return {tessellated_quads, tessellated_points};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SUBDIVISION CURVES
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines.
template <typename T>
static pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines_, const vector<T>& vertices_, int subdivisions) {
  // early exit
  if (lines_.empty() || vertices_.empty() || subdivisions < 1)
    return {lines_, vertices_};

  // copy data
  auto lines    = lines_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // create vertices
    auto tvertices = vector<T>{};
    tvertices.reserve(vertices.size() + lines.size());
    for (auto& vertex : vertices) tvertices.push_back(vertex);
    for (auto& line : lines) {
      tvertices.push_back((vertices[line.x] + vertices[line.y]) / 2);
    }
    // create lines
    auto tlines = vector<vec2i>{};
    tlines.reserve(lines.size() * 2);
    auto line_vertex = [nverts = (int)vertices.size()](
                           size_t line_id) { return nverts + (int)line_id; };
    for (auto&& [line_id, line] : enumerate(lines)) {
      tlines.push_back({line.x, line_vertex(line_id)});
      tlines.push_back({line_vertex(line_id), line.y});
    }

    // swap
    lines.swap(tlines);
    vertices.swap(tvertices);
  }

  // done
  return {lines, vertices};
}

// Check if it is a boundary vertex
inline vector<int> get_boundary_vertices(
    int npoints, const vector<vec2i>& lines) {
  // count edge degrees
  auto degree = vector<int>(npoints, 0);
  for (auto& [i, j] : lines) {
    degree[i] += 1;
    degree[j] += 1;
  }

  // find boundary vertices
  auto boundary = vector<int>{};
  for (auto index : range(npoints)) {
    if (degree[index] <= 1) boundary.push_back(index);
  }

  // done
  return boundary;
}

// B-spline subdivision
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vertices, int subdivisions) {
  return subdivide_bspline(lines, vertices, {}, subdivisions);
}

// B-spline subdivision
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines_, const vector<T>& vertices_,
    const vector<int>& corners, int subdivisions) {
  // early exit
  if (lines_.empty() || vertices_.empty() || subdivisions < 1)
    return {lines_, vertices_};

  // copy data
  auto lines    = lines_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // split elements ------------------------------------
    auto [tlines, tvertices] = subdivide_lines(lines, vertices);

    // boundary ------------------------------------------
    auto boundary = get_boundary_vertices((int)tvertices.size(), tlines);
    auto valence  = vector<int>(tvertices.size(), 1);
    auto tcorners = vector<int>{};
    for (auto index : remove_duplicates(join(boundary, corners))) {
      tcorners.push_back(index);
      valence[index] = 0;
    }

    // averaging pass ----------------------------------
    auto avertices = vector<T>(tvertices.size(), T{});
    auto acounts   = vector<int>(tvertices.size(), 0);
    for (auto& vid : tcorners) {
      auto& vertex = tvertices[vid];
      if (valence[vid] != 0) continue;
      avertices[vid] += vertex;
      acounts[vid] += 1;
    }
    for (auto& line : tlines) {
      auto midpoint = (tvertices[line.x] + tvertices[line.y]) / 2;
      for (auto vid : line) {
        if (valence[vid] != 1) continue;
        avertices[vid] += midpoint;
        acounts[vid] += 1;
      }
    }
    for (auto index : range(tvertices.size()))
      avertices[index] /= (float)acounts[index];
    tvertices = avertices;

    // swap
    lines.swap(tlines);
    vertices.swap(tvertices);
  }

  // done
  return {lines, vertices};
}

// Subdivide beziers.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers_, const vector<T>& vertices_,
    int subdivisions) {
  // early exit
  if (beziers_.empty() || vertices_.empty() || subdivisions < 1)
    return {beziers_, vertices_};

  // setup subdivision
  auto beziers  = beziers_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges
    auto vmap      = unordered_map<int, int>();
    auto tvertices = vector<T>();
    auto tbeziers  = vector<vec4i>();
    for (auto& bezier : beziers) {
      if (vmap.find(bezier.x) == vmap.end()) {
        vmap[bezier.x] = (int)tvertices.size();
        tvertices.push_back(vertices[bezier.x]);
      }
      if (vmap.find(bezier.w) == vmap.end()) {
        vmap[bezier.w] = (int)tvertices.size();
        tvertices.push_back(vertices[bezier.w]);
      }
      auto bo = (int)tvertices.size();
      tbeziers.push_back({vmap.at(bezier.x), bo + 0, bo + 1, bo + 2});
      tbeziers.push_back({bo + 2, bo + 3, bo + 4, vmap.at(bezier.w)});
      tvertices.push_back(vertices[bezier.x] / 2 + vertices[bezier.y] / 2);
      tvertices.push_back(vertices[bezier.x] / 4 + vertices[bezier.y] / 2 +
                          vertices[bezier.z] / 4);
      tvertices.push_back(
          vertices[bezier.x] / 8 + vertices[bezier.y] * ((float)3 / (float)8) +
          vertices[bezier.z] * ((float)3 / (float)8) + vertices[bezier.w] / 8);
      tvertices.push_back(vertices[bezier.y] / 4 + vertices[bezier.z] / 2 +
                          vertices[bezier.w] / 4);
      tvertices.push_back(vertices[bezier.z] / 2 + vertices[bezier.w] / 2);
    }

    // swap
    beziers.swap(tbeziers);
    vertices.swap(tvertices);
  }

  // done
  return {beziers, vertices};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EDGE AND ADJACENCIES FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Edge comparator
inline bool compare_edges(vec2i a, vec2i b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

// Search for edge index
inline vec2i make_edge(vec2i e) { return e.x < e.y ? e : vec2i{e.y, e.x}; }
inline vec2i make_edge(int a, int b) { return make_edge(vec2i{a, b}); }
inline int   search_edge(const vector<vec2i>& edges, vec2i edge_) {
  auto edge = make_edge(edge_);
  auto pos = std::lower_bound(edges.begin(), edges.end(), edge, compare_edges);
  if (pos == edges.end() || *pos != edge) return -1;
  return (int)(pos - edges.begin());
}

// Get edges of meshes
inline vector<vec2i> triangles_edges(const vector<vec3i>& triangles) {
  return faces_edges(triangles, {});
}
inline vector<vec2i> quads_edges(const vector<vec4i>& quads) {
  return faces_edges({}, quads);
}
inline vector<vec2i> faces_edges(
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
  // allocate edges
  auto edges = vector<vec2i>();
  edges.reserve(triangles.size() * 3 + quads.size() * 4);

  // add edges
  for (auto& triangle : triangles) {
    edges.push_back(make_edge(triangle.x, triangle.y));
    edges.push_back(make_edge(triangle.y, triangle.z));
    edges.push_back(make_edge(triangle.z, triangle.x));
  }
  for (auto& quad : quads) {
    edges.push_back(make_edge(quad.x, quad.y));
    edges.push_back(make_edge(quad.y, quad.z));
    edges.push_back(make_edge(quad.z, quad.w));
    edges.push_back(make_edge(quad.w, quad.x));
  }

  // remove degenerates
  std::erase_if(edges, [](vec2i edge) { return edge.x == edge.y; });

  // sort edges
  std::sort(edges.begin(), edges.end(), compare_edges);

  // remove duplicates
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  // done
  return edges;
}

// Get edges and boundaries of meshes
inline vector<vec2i> triangles_boundaries(const vector<vec3i>& triangles) {
  return faces_boundaries(triangles, {});
}
inline vector<vec2i> quads_boundaries(const vector<vec4i>& quads) {
  return faces_boundaries({}, quads);
}
inline vector<vec2i> faces_boundaries(
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
  // allocate edges
  auto edges = vector<vec2i>();
  edges.reserve(triangles.size() * 3 + quads.size() * 4);

  // add edges
  for (auto& triangle : triangles) {
    edges.push_back(make_edge(triangle.x, triangle.y));
    edges.push_back(make_edge(triangle.y, triangle.z));
    edges.push_back(make_edge(triangle.z, triangle.x));
  }
  for (auto& quad : quads) {
    edges.push_back(make_edge(quad.x, quad.y));
    edges.push_back(make_edge(quad.y, quad.z));
    edges.push_back(make_edge(quad.z, quad.w));
    edges.push_back(make_edge(quad.w, quad.x));
  }

  // remove degenerates
  std::erase_if(edges, [](vec2i edge) { return edge.x == edge.y; });

  // sort edges
  std::sort(edges.begin(), edges.end(), compare_edges);

  // construct boundaries
  auto boundaries = vector<vec2i>{};
  for (auto index : range(edges.size())) {
    auto cur  = edges[index];
    auto prev = index == 0 ? vec2i{-1, -1} : edges[index - 1];
    auto next = index == edges.size() - 1 ? vec2i{-1, -1} : edges[index + 1];
    if (cur != next && cur != prev) boundaries.push_back(cur);
  }

  // done
  return boundaries;
}

// Get edges and boundaries of meshes
inline pair<vector<vec2i>, vector<vec2i>> triangles_edges_boundaries(
    const vector<vec3i>& triangles) {
  return faces_edges_boundaries(triangles, {});
}
inline pair<vector<vec2i>, vector<vec2i>> quads_edges_boundaries(
    const vector<vec4i>& quads) {
  return faces_edges_boundaries({}, quads);
}
inline pair<vector<vec2i>, vector<vec2i>> faces_edges_boundaries(
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
  // allocate edges
  auto edges = vector<vec2i>();
  edges.reserve(triangles.size() * 3 + quads.size() * 4);

  // add edges
  for (auto& triangle : triangles) {
    edges.push_back(make_edge(triangle.x, triangle.y));
    edges.push_back(make_edge(triangle.y, triangle.z));
    edges.push_back(make_edge(triangle.z, triangle.x));
  }
  for (auto& quad : quads) {
    edges.push_back(make_edge(quad.x, quad.y));
    edges.push_back(make_edge(quad.y, quad.z));
    edges.push_back(make_edge(quad.z, quad.w));
    edges.push_back(make_edge(quad.w, quad.x));
  }

  // remove degenerates
  std::erase_if(edges, [](vec2i edge) { return edge.x == edge.y; });

  // sort edges
  std::sort(edges.begin(), edges.end(), compare_edges);

  // construct boundaries
  auto boundaries = vector<vec2i>{};
  for (auto index : range(edges.size())) {
    auto cur  = edges[index];
    auto prev = index == 0 ? vec2i{-1, -1} : edges[index - 1];
    auto next = index == edges.size() - 1 ? vec2i{-1, -1} : edges[index + 1];
    if (cur != next && cur != prev) boundaries.push_back(cur);
  }

  // remove duplicates
  edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

  // done
  return {edges, boundaries};
}

// Build adjacencies between faces (sorted counter-clockwise)
inline vector<vec3i> face_adjacencies(const vector<vec3i>& triangles) {
  struct face_edge {
    vec2i edge;
    int   face;
    int   edge_index;
  };
  auto get_edge = [](vec3i triangle, int i) -> vec2i {
    auto x = triangle[i], y = triangle[i < 2 ? i + 1 : 0];
    return x < y ? vec2i{x, y} : vec2i{y, x};
  };
  auto edges = vector<face_edge>{};
  edges.reserve(triangles.size() * 3);
  for (auto tid : range((int)triangles.size())) {
    for (auto k = 0; k < 3; ++k) {
      auto edge = get_edge(triangles[tid], k);
      edges.push_back({edge, tid, k});
    }
  }
  sort(edges.begin(), edges.end(), [](auto a_, auto b_) {
    auto a = a_.edge, b = b_.edge;
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  });

  auto adjacencies = vector<vec3i>{triangles.size(), vec3i{-1, -1, -1}};
  for (auto eid : range((int)edges.size() - 1)) {
    if (edges[eid].edge == edges[eid + 1].edge) continue;
    adjacencies[edges[eid].face][edges[eid].edge_index] = edges[eid + 1].face;
    adjacencies[edges[eid + 1].face][edges[eid + 1].edge_index] =
        edges[eid].face;
  }

  return adjacencies;
}

// Build adjacencies between vertices (sorted counter-clockwise)
inline vector<vector<int>> vertex_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies) {
  auto find_index = [](vec3i v, int x) {
    if (v.x == x) return 0;
    if (v.y == x) return 1;
    if (v.z == x) return 2;
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
inline vector<vector<int>> vertex_to_faces_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies) {
  auto find_index = [](vec3i v, int x) {
    if (v.x == x) return 0;
    if (v.y == x) return 1;
    if (v.z == x) return 2;
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
inline vector<vector<int>> ordered_boundaries(const vector<vec3i>& triangles,
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
// SUBDIVISION SURFACES
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide triangle.
template <typename T>
static pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles_, const vector<T>& vertices_,
    int subdivisions) {
  // early exit
  if (triangles_.empty() || vertices_.empty() || subdivisions < 1)
    return {triangles_, vertices_};

  // copy data
  auto triangles = triangles_;
  auto vertices  = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges
    auto edges = triangles_edges(triangles);

    // number of elements
    auto nverts = (int)vertices.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)triangles.size();

    // indices
    auto edge_vertex = [&edges, nverts](vec2i edge) {
      return nverts + search_edge(edges, edge);
    };

    // create vertices
    auto tvertices = vector<T>{};
    tvertices.reserve(nverts + nedges);
    for (auto& vertex : vertices) tvertices.push_back(vertex);
    for (auto& edge : edges)
      tvertices.push_back((vertices[edge.x] + vertices[edge.y]) / 2);
    // create triangles
    auto ttriangles = vector<vec3i>{};
    ttriangles.reserve(triangles.size() * 4);
    for (auto& triangle : triangles) {
      ttriangles.push_back({triangle.x, edge_vertex({triangle.x, triangle.y}),
          edge_vertex({triangle.z, triangle.x})});
      ttriangles.push_back({triangle.y, edge_vertex({triangle.y, triangle.z}),
          edge_vertex({triangle.x, triangle.y})});
      ttriangles.push_back({triangle.z, edge_vertex({triangle.z, triangle.x}),
          edge_vertex({triangle.y, triangle.z})});
      ttriangles.push_back({edge_vertex({triangle.x, triangle.y}),
          edge_vertex({triangle.y, triangle.z}),
          edge_vertex({triangle.z, triangle.x})});
    }

    // swap for next iteration
    triangles.swap(ttriangles);
    vertices.swap(tvertices);
  }

  // done
  return {triangles, vertices};
}

// Subdivide quads.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads_, const vector<T>& vertices_, int subdivisions) {
  // early exit
  if (quads_.empty() || vertices_.empty() || subdivisions < 1)
    return {quads_, vertices_};

  // setup subdivision
  auto quads    = quads_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges
    auto edges = quads_edges(quads);

    // number of elements
    auto nverts = (int)vertices.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();

    // indices
    auto edge_vertex = [&edges, nverts](vec2i edge) {
      return nverts + search_edge(edges, edge);
    };
    auto quad_vertex = [nverts, nedges](size_t quad_id) {
      return nverts + nedges + (int)quad_id;
    };

    // create vertices
    auto tvertices = vector<T>{};
    tvertices.reserve(nverts + nedges + nfaces);
    for (auto& vertex : vertices) tvertices.push_back(vertex);
    for (auto& edge : edges)
      tvertices.push_back((vertices[edge.x] + vertices[edge.y]) / 2);
    for (auto& quad : quads) {
      if (quad.z != quad.w) {
        tvertices.push_back((vertices[quad.x] + vertices[quad.y] +
                                vertices[quad.z] + vertices[quad.w]) /
                            4);
      } else {
        tvertices.push_back(
            (vertices[quad.x] + vertices[quad.y] + vertices[quad.z]) / 3);
      }
    }
    // create quads
    auto tquads = vector<vec4i>{};
    tquads.reserve(quads.size() * 4);
    for (auto&& [quad_id, quad] : enumerate(quads)) {
      if (quad.z != quad.w) {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.w, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.w}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
        tquads.push_back({quad.w, edge_vertex({quad.w, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.w})});
      } else {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
      }
    }

    // swap for next iteration
    quads.swap(tquads);
    vertices.swap(tvertices);
  }

  // done
  return {quads, vertices};
}

// Subdivide Catmull-Clark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads_, const vector<T>& vertices_, int subdivisions,
    bool lock_boundary) {
  // early exit
  if (quads_.empty() || vertices_.empty() || subdivisions < 1)
    return {quads_, vertices_};

  // setup subdivision
  auto quads    = quads_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges and boundaries
    auto [edges, boundary] = quads_edges_boundaries(quads);

    // number of elements
    auto nverts = (int)vertices.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();

    // indices
    auto edge_vertex = [&edges, nverts](vec2i edge) {
      return nverts + search_edge(edges, edge);
    };
    auto quad_vertex = [nverts, nedges](size_t quad_id) {
      return nverts + nedges + (int)quad_id;
    };

    // split elements ------------------------------------
    // create vertices
    auto tvertices = vector<T>{};
    tvertices.reserve(nverts + nedges + nfaces);
    for (auto& vertex : vertices) tvertices.push_back(vertex);
    for (auto& edge : edges)
      tvertices.push_back((vertices[edge.x] + vertices[edge.y]) / 2);
    for (auto& quad : quads) {
      if (quad.z != quad.w) {
        tvertices.push_back((vertices[quad.x] + vertices[quad.y] +
                                vertices[quad.z] + vertices[quad.w]) /
                            4);
      } else {
        tvertices.push_back(
            (vertices[quad.x] + vertices[quad.y] + vertices[quad.z]) / 3);
      }
    }
    // create quads
    auto tquads = vector<vec4i>{};
    tquads.reserve(quads.size() * 4);
    for (auto&& [quad_id, quad] : enumerate(quads)) {
      if (quad.z != quad.w) {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.w, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.w}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
        tquads.push_back({quad.w, edge_vertex({quad.w, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.w})});
      } else {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
      }
    }

    // split boundary
    auto tboundary = vector<vec2i>{};
    tboundary.reserve(boundary.size());
    for (auto& edge : boundary) {
      tboundary.push_back({edge.x, edge_vertex(edge)});
      tboundary.push_back({edge_vertex(edge), edge.y});
    }

    // setup creases -----------------------------------
    auto tcrease_edges = vector<vec2i>{};
    auto tcrease_verts = vector<int>{};
    if (lock_boundary) {
      for (auto& b : tboundary) {
        tcrease_verts.push_back(b.x);
        tcrease_verts.push_back(b.y);
      }
    } else {
      for (auto& b : tboundary) tcrease_edges.push_back(b);
    }

    // define vertices valence ---------------------------
    auto tvalences = vector<int>(tvertices.size(), 2);
    for (auto& edge : tboundary) {
      tvalences[edge.x] = (lock_boundary) ? 0 : 1;
      tvalences[edge.y] = (lock_boundary) ? 0 : 1;
    }

    // averaging pass ----------------------------------
    auto avertices = vector<T>(tvertices.size(), T());
    auto acounts   = vector<int>(tvertices.size(), 0);
    for (auto& point : tcrease_verts) {
      if (tvalences[point] != 0) continue;
      avertices[point] += tvertices[point];
      acounts[point] += 1;
    }
    for (auto& edge : tcrease_edges) {
      auto centroid = (tvertices[edge.x] + tvertices[edge.y]) / 2;
      for (auto vid : {edge.x, edge.y}) {
        if (tvalences[vid] != 1) continue;
        avertices[vid] += centroid;
        acounts[vid] += 1;
      }
    }
    for (auto& quad : tquads) {
      auto centroid = (tvertices[quad.x] + tvertices[quad.y] +
                          tvertices[quad.z] + tvertices[quad.w]) /
                      4;
      for (auto vid : {quad.x, quad.y, quad.z, quad.w}) {
        if (tvalences[vid] != 2) continue;
        avertices[vid] += centroid;
        acounts[vid] += 1;
      }
    }
    for (auto i : range(tvertices.size())) avertices[i] /= (float)acounts[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i : range(tvertices.size())) {
      if (tvalences[i] != 2) continue;
      avertices[i] = tvertices[i] +
                     (avertices[i] - tvertices[i]) * (4 / (float)acounts[i]);
    }
    tvertices = avertices;

    // swap for next iteration
    quads.swap(tquads);
    vertices.swap(tvertices);
  }

  // done
  return {quads, vertices};
}

// Subdivide Catmull-Clark.
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_creased(
    const vector<vec4i>& quads_, const vector<T>& vertices_,
    const vector<vec2i>& creases_, const vector<int>& corners, int subdivisions,
    bool lock_boundary) {
  // early exit
  if (quads_.empty() || vertices_.empty() || subdivisions < 1)
    return {quads_, vertices_};

  // setup subdivision
  auto quads    = quads_;
  auto vertices = vertices_;
  auto creases  = creases_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges and boundaries
    auto [edges, boundaries] = quads_edges_boundaries(quads);

    // number of elements
    auto nverts = (int)vertices.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();

    // indices
    auto edge_vertex = [&edges, nverts](vec2i edge) {
      return nverts + search_edge(edges, edge);
    };
    auto quad_vertex = [nverts, nedges](size_t quad_id) {
      return nverts + nedges + (int)quad_id;
    };

    // split elements ------------------------------------
    // create vertices
    auto tvertices = vector<T>{};
    tvertices.reserve(vertices.size() + edges.size() + quads.size());
    for (auto& vertex : vertices) tvertices.push_back(vertex);
    for (auto& edge : edges)
      tvertices.push_back((vertices[edge.x] + vertices[edge.y]) / 2);
    for (auto& quad : quads) {
      if (quad.z != quad.w) {
        tvertices.push_back((vertices[quad.x] + vertices[quad.y] +
                                vertices[quad.z] + vertices[quad.w]) /
                            4);
      } else {
        tvertices.push_back(
            (vertices[quad.x] + vertices[quad.y] + vertices[quad.z]) / 3);
      }
    }
    // create quads
    auto tquads = vector<vec4i>{};
    tquads.reserve(quads.size() * 4);
    for (auto&& [quad_id, quad] : enumerate(quads)) {
      if (quad.z != quad.w) {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.w, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.w}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
        tquads.push_back({quad.w, edge_vertex({quad.w, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.w})});
      } else {
        tquads.push_back({quad.x, edge_vertex({quad.x, quad.y}),
            quad_vertex(quad_id), edge_vertex({quad.z, quad.x})});
        tquads.push_back({quad.y, edge_vertex({quad.y, quad.z}),
            quad_vertex(quad_id), edge_vertex({quad.x, quad.y})});
        tquads.push_back({quad.z, edge_vertex({quad.z, quad.x}),
            quad_vertex(quad_id), edge_vertex({quad.y, quad.z})});
      }
    }

    // join boundary and creases
    boundaries.insert(boundaries.end(), creases.begin(), creases.end());
    boundaries.erase(
        std::unique(boundaries.begin(), boundaries.end()), boundaries.end());

    // split boundary and creases
    auto tboundaries = vector<vec2i>{};
    tboundaries.reserve(boundaries.size());
    for (auto& edge : boundaries) {
      tboundaries.push_back({edge.x, edge_vertex(edge)});
      tboundaries.push_back({edge_vertex(edge), edge.y});
    }

    // split creases to return them
    auto tcreases = vector<vec2i>{};
    for (auto& edge : creases) {
      tcreases.push_back({edge.x, edge_vertex(edge)});
      tcreases.push_back({edge_vertex(edge), edge.y});
    }
    tcreases.erase(
        std::unique(tcreases.begin(), tcreases.end()), tcreases.end());

    // split corners
    auto tcorners = corners;
    if (lock_boundary) {
      for (auto& edge : tboundaries) {
        tcorners.push_back(edge.x);
        tcorners.push_back(edge.y);
      }
    }

    // define vertices valence ---------------------------
    auto tvalences = vector<int>(tvertices.size(), 2);
    for (auto& edge : tboundaries) {
      tvalences[edge.x] = (lock_boundary) ? 0 : 1;
      tvalences[edge.y] = (lock_boundary) ? 0 : 1;
    }
    for (auto& corner : tcorners) {  // apply later to override boundaries
      tvalences[corner] = 0;
    }

    // join boundary and creases

    // averaging pass ----------------------------------
    auto avertices = vector<T>(tvertices.size(), T());
    auto acounts   = vector<int>(tvertices.size(), 0);
    for (auto& point : tcorners) {
      if (tvalences[point] != 0) continue;
      avertices[point] += tvertices[point];
      acounts[point] += 1;
    }
    for (auto& edge : tboundaries) {
      auto centroid = (tvertices[edge.x] + tvertices[edge.y]) / 2;
      for (auto vid : {edge.x, edge.y}) {
        if (tvalences[vid] != 1) continue;
        avertices[vid] += centroid;
        acounts[vid] += 1;
      }
    }
    for (auto& quad : tquads) {
      auto centroid = (tvertices[quad.x] + tvertices[quad.y] +
                          tvertices[quad.z] + tvertices[quad.w]) /
                      4;
      for (auto vid : {quad.x, quad.y, quad.z, quad.w}) {
        if (tvalences[vid] != 2) continue;
        avertices[vid] += centroid;
        acounts[vid] += 1;
      }
    }
    for (auto i : range(tvertices.size())) avertices[i] /= (float)acounts[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i : range(tvertices.size())) {
      if (tvalences[i] != 2) continue;
      avertices[i] = tvertices[i] +
                     (avertices[i] - tvertices[i]) * (4 / (float)acounts[i]);
    }
    tvertices = avertices;

    // swap for next iteration
    quads.swap(tquads);
    vertices.swap(tvertices);
    creases.swap(tcreases);
  }

  // done
  return {quads, vertices};
}

// Subdivide quads using Carmull-Clark subdivision rules.
// VERSION USED FOR ILLUSTRATION PURPOSES
template <typename T>
static pair<vector<vec4i>, vector<T>> subdivide_catmullclark_illustration(
    const vector<vec4i>& quads_, const vector<T>& vertices_, bool uncorrected,
    int subdivisions) {
  // early exit
  if (quads_.empty() || vertices_.empty() || subdivisions < 1)
    return {quads_, vertices_};

  // setup subdivision
  auto quads    = quads_;
  auto vertices = vertices_;

  // subdivide
  for (auto subdivision : range(subdivisions)) {
    // get edges and boundaries
    auto [edges, boundary] = quads_edges_boundaries(quads);

    // number of elements
    auto nverts = (int)vertices.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();

    // indices
    auto edge_vertex = [&edges, nverts](vec2i edge) {
      return nverts + search_edge(edges, edge);
    };
    auto quad_vertex = [nverts, nedges](size_t quad_id) {
      return nverts + nedges + (int)quad_id;
    };

    // split elements ------------------------------------
    // create vertices
    auto tvertices = vector<T>(nverts + nedges + nfaces);
    for (auto i = 0; i < nverts; i++) tvertices[i] = vertices[i];
    for (auto i = 0; i < nedges; i++) {
      auto e                = edges[i];
      tvertices[nverts + i] = (vertices[e.x] + vertices[e.y]) / 2;
    }
    for (auto i = 0; i < nfaces; i++) {
      auto q = quads[i];
      if (q.z != q.w) {
        tvertices[nverts + nedges + i] =
            (vertices[q.x] + vertices[q.y] + vertices[q.z] + vertices[q.w]) / 4;
      } else {
        tvertices[nverts + nedges + i] =
            (vertices[q.x] + vertices[q.y] + vertices[q.z]) / 3;
      }
    }
    // create quads
    auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
    auto qi     = 0;
    for (auto i : range(nfaces)) {
      auto q = quads[i];
      if (q.z != q.w) {
        tquads[qi++] = {q.x, nverts + search_edge(edges, {q.x, q.y}),
            nverts + nedges + i, nverts + search_edge(edges, {q.w, q.x})};
        tquads[qi++] = {q.y, nverts + search_edge(edges, {q.y, q.z}),
            nverts + nedges + i, nverts + search_edge(edges, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + search_edge(edges, {q.z, q.w}),
            nverts + nedges + i, nverts + search_edge(edges, {q.y, q.z})};
        tquads[qi++] = {q.w, nverts + search_edge(edges, {q.w, q.x}),
            nverts + nedges + i, nverts + search_edge(edges, {q.z, q.w})};
      } else {
        tquads[qi++] = {q.x, nverts + search_edge(edges, {q.x, q.y}),
            nverts + nedges + i, nverts + search_edge(edges, {q.z, q.x})};
        tquads[qi++] = {q.y, nverts + search_edge(edges, {q.y, q.z}),
            nverts + nedges + i, nverts + search_edge(edges, {q.x, q.y})};
        tquads[qi++] = {q.z, nverts + search_edge(edges, {q.z, q.x}),
            nverts + nedges + i, nverts + search_edge(edges, {q.y, q.z})};
      }
    }
    tquads.resize(qi);

    // averaging pass ----------------------------------
    auto avertices = vector<T>(tvertices.size(), T{});
    auto acounts   = vector<int>(tvertices.size(), 0);
    for (auto& q : tquads) {
      auto c =
          (tvertices[q.x] + tvertices[q.y] + tvertices[q.z] + tvertices[q.w]) /
          4;
      for (auto vid : {q.x, q.y, q.z, q.w}) {
        avertices[vid] += c;
        acounts[vid] += 1;
      }
    }
    for (auto i = 0; i < tvertices.size(); i++)
      avertices[i] /= (float)acounts[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    if (!uncorrected) {
      for (auto i = 0; i < tvertices.size(); i++) {
        avertices[i] = tvertices[i] +
                       (avertices[i] - tvertices[i]) * (4 / (float)acounts[i]);
      }
    }
    tvertices = avertices;

    // swap for next iteration
    quads.swap(tquads);
    vertices.swap(tvertices);
  }

  // done
  return {quads, vertices};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DISPLACEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// Displace vertices
inline vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const image<float>& displacement, float scale, float offset) {
  if (texcoords.empty()) return positions;
  auto displaced = positions;
  for (auto idx : range(displaced.size())) {
    displaced[idx] += normals[idx] *
                      (eval_image(displacement, texcoords[idx]) - offset) *
                      scale;
  }
  return displaced;
}
inline vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const image<vec4f>& displacement, float scale, float offset) {
  if (texcoords.empty()) return positions;
  auto displaced = positions;
  for (auto idx : range(displaced.size())) {
    displaced[idx] +=
        normals[idx] *
        (mean(xyz(eval_image(displacement, texcoords[idx]))) - offset) * scale;
  }
  return displaced;
}

// Displacement
inline shape_data displace_shape(const shape_data& shape,
    const image<float>& displacement, float height, float offset) {
  if (displacement.empty() || shape.texcoords.empty() ||
      shape.normals.empty() || (shape.triangles.empty() && shape.quads.empty()))
    return shape;
  auto displaced      = shape;
  displaced.positions = displace_vertices(displaced.positions,
      displaced.normals, displaced.texcoords, displacement, height, offset);
  displaced.normals   = compute_normals(displaced);
  return displaced;
}
inline shape_data displace_shape(const shape_data& shape,
    const image<vec4f>& displacement, float height, float offset) {
  if (displacement.empty() || shape.texcoords.empty() ||
      shape.normals.empty() || (shape.triangles.empty() && shape.quads.empty()))
    return shape;
  auto displaced      = shape;
  displaced.positions = displace_vertices(displaced.positions,
      displaced.normals, displaced.texcoords, displacement, height, offset);
  displaced.normals   = compute_normals(displaced);
  return displaced;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIGNED-DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Sphere SDF
inline float sdf_sphere(vec3f position, vec3f center, float radius) {
  return length(position - center) - radius;
}
// Box SDF
inline float sdf_box(vec3f position, vec3f center, vec3f size) {
  auto q = abs(position - center) - size;
  return length(max(q, 0.0f)) + min(max(q), 0.0f);
}
// Capsule SDF
inline float sdf_capsule(vec3f position, vec3f a, vec3f b, float radius) {
  auto  pos_a = position - a, axis = b - a;
  float height = clamp(dot(pos_a, axis) / dot(axis, axis), 0.0f, 1.0f);
  return length(pos_a - axis * height) - radius;
}

// Union SDF
inline float sdf_union(float sdf1, float sdf2) { return min(sdf1, sdf2); }
// Intersection SDF
inline float sdf_intersection(float sdf1, float sdf2) {
  return max(sdf1, sdf2);
}
// Difference SDF
inline float sdf_difference(float sdf1, float sdf2) { return max(sdf1, -sdf2); }

// Smooth min/max
inline float smin(float a, float b, float k = 32) {
  return -log2(exp2(-k * a) + exp2(-k * b)) / k;
}
inline float smax(float a, float b, float k = 32) {
  return +log2(exp2(+k * a) + exp2(+k * b)) / k;
}

// Soft union SDF
inline float sdf_smooth_union(float sdf1, float sdf2, float k) {
  return smin(sdf1, sdf2);
}
// Soft intersection SDF
inline float sdf_smooth_intersection(float sdf1, float sdf2, float k) {
  return smax(sdf1, sdf2);
}
// Soft difference SDF
inline float sdf_smooth_difference(float sdf1, float sdf2, float k) {
  return smax(sdf1, -sdf2);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// NOISE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// works for 2 and 3 dimensions
// http://www.jcgt.org/published/0009/03/02/
template <size_t N>
inline array<uint, N> hash_pcg(array<uint, N> v_) {
  auto v = v_;
  for (auto i : range(N)) v[i] = v[i] * 1664525u + 1013904223u;
  for (auto i : range(N)) v[i] += v[(i + (size_t)1) % N] * 1664525u;
  for (auto i : range(N)) v[i] = v[i] ^ (v[i] >> 16u);
  for (auto i : range(N)) v[i] += v[(i + (size_t)1) % N] * 1664525u;
  for (auto i : range(N)) v[i] = v[i] ^ (v[i] >> 16u);
  return v;
}
template <size_t N>
inline array<uint, N> hash_pcg(array<int, N> v_) {
  auto vu = array<uint, N>{};
  for (auto i : range(N)) vu[i] = *(uint*)&v_[i];
  return hash_pcg(vu);
}
inline vec2f hash_pcg(vec2i v) {
  auto hashed = hash_pcg<2>(array<int, 2>{(int)v.x, (int)v.y});
  return {
      (float)hashed[0] / (float)UINT_MAX, (float)hashed[1] / (float)UINT_MAX};
}
inline vec3f hash_pcg(vec3i v) {
  auto hashed = hash_pcg<3>(array<int, 3>{(int)v.x, (int)v.y, (int)v.z});
  return {(float)hashed[0] / (float)UINT_MAX,
      (float)hashed[1] / (float)UINT_MAX, (float)hashed[2] / (float)UINT_MAX};
}

// the implementation follows ideas from "Building up Perlin Noise"
// http://eastfarthing.com/blog/2015-04-21-noise/

// Interpolants and surflets
inline float noise_kernel_cubic(float x_) {
  auto x = abs(x_);
  return x >= 1 ? 0 : (1 - (3 * x * x - 2 * x * x * x));
}
inline float noise_kernel_cubic(vec2f x) {
  return noise_kernel_cubic(x.x) * noise_kernel_cubic(x.y);
}
inline float noise_kernel_cubic(vec3f x) {
  return noise_kernel_cubic(x.x) * noise_kernel_cubic(x.y) *
         noise_kernel_cubic(x.z);
}
inline float noise_kernel_quintic(float x_) {
  auto x = abs(x_);
  return x >= 1 ? 0 : (1 - ((6 * x - 15) * x + 10) * x * x * x);
}
inline float noise_kernel_quintic(vec2f x) {
  return noise_kernel_quintic(x.x) * noise_kernel_quintic(x.y);
}
inline float noise_kernel_quintic(vec3f x) {
  return noise_kernel_quintic(x.x) * noise_kernel_quintic(x.y) *
         noise_kernel_quintic(x.z);
}

// Value noise
inline float value_noise_rvalue(vec2i cell) { return hash_pcg(cell).x; }
inline float value_noise_rvalue(vec3i cell) { return hash_pcg(cell).x; }
inline float value_noise(vec2f position) {
  auto cell  = (vec2i)floor(position);
  auto noise = 0.0f;
  for (auto offset : range(vec2i{2, 2})) {
    auto index = cell + offset;
    noise += value_noise_rvalue(index) * noise_kernel_cubic(position - index);
  }
  return noise;
}
inline float value_noise(vec3f position) {
  auto cell  = (vec3i)floor(position);
  auto noise = 0.0f;
  for (auto offset : range(vec3i{2, 2, 2})) {
    auto index = cell + offset;
    noise += value_noise_rvalue(index) * noise_kernel_cubic(position - index);
  }
  return noise;
}

// Gradient noise
inline vec2f gradient_noise_rvector(vec2i cell) {
  return normalize(2 * hash_pcg(cell) - 1);
}
inline vec3f gradient_noise_rvector(vec3i cell) {
  return normalize(2 * hash_pcg(cell) - 1);
}
inline float gradient_noise_rgradient(vec2f position, vec2i cell) {
  return dot(position - cell, gradient_noise_rvector(cell));
}
inline float gradient_noise_rgradient(vec3f position, vec3i cell) {
  return dot(position - cell, gradient_noise_rvector(cell));
}
inline float gradient_noise(vec2f position) {
  auto cell  = (vec2i)floor(position);
  auto noise = 0.0f;
  for (auto offset : range(vec2i{2, 2})) {
    auto index = cell + offset;
    noise += gradient_noise_rgradient(position, index) *
             noise_kernel_cubic(position - index);
  }
  noise = (noise + 1) / 2;
  return noise;
}
inline float gradient_noise(vec3f position) {
  auto cell  = (vec3i)floor(position);
  auto noise = 0.0f;
  for (auto offset : range(vec3i{2, 2, 2})) {
    auto index = cell + offset;
    noise += gradient_noise_rgradient(position, index) *
             noise_kernel_cubic(position - index);
  }
  noise = (noise + 1) / 2;
  return noise;
}

// Cell noise
inline vec2f cell_noise_rpoint(vec2i cell) { return cell + hash_pcg(cell); }
inline vec3f cell_noise_rpoint(vec3i cell) { return cell + hash_pcg(cell); }
inline float cell_noise(vec2f position, bool inverted) {
  auto cell  = (vec2i)floor(position);
  auto noise = 2.0f;
  for (auto offset_ : range(vec2i{3, 3})) {
    auto offset = offset_ - 1;
    auto index  = cell + offset;
    auto dist   = length(position - cell_noise_rpoint(index));
    noise       = min(noise, dist);
  }
  return !inverted ? noise : (1 - noise);
}
inline float cell_noise(vec3f position, bool inverted) {
  auto cell  = (vec3i)floor(position);
  auto noise = 2.0f;
  for (auto offset_ : range(vec3i{3, 3, 3})) {
    auto offset = offset_ - 1;
    auto index  = cell + offset;
    auto dist   = length(position - cell_noise_rpoint(index));
    noise       = min(noise, dist);
  }
  return !inverted ? noise : (1 - noise);
}

// Cell noise
inline float voronoi_noise(vec2f position) {
  auto cell  = (vec2i)floor(position);
  auto noise = value_noise_rvalue(cell);
  auto mdist = 2.0f;
  for (auto offset_ : range(vec2i{3, 3})) {
    auto offset = offset_ - 1;
    auto index  = cell + offset;
    auto dist   = length(position - cell_noise_rpoint(index));
    if (dist < mdist) {
      noise = value_noise_rvalue(index);
      mdist = dist;
    }
  }
  return noise;
}
inline float voronoi_noise(vec3f position) {
  auto cell  = (vec3i)floor(position);
  auto noise = value_noise_rvalue(cell);
  auto mdist = 2.0f;
  for (auto offset_ : range(vec3i{3, 3, 3})) {
    auto offset = offset_ - 1;
    auto index  = cell + offset;
    auto dist   = length(position - cell_noise_rpoint(index));
    if (dist < mdist) {
      noise = value_noise_rvalue(index);
      mdist = dist;
    }
  }
  return noise;
}

// Fractal noise
template <typename Noise>
inline float fractal_noise(Noise&& base_noise, vec2f position, int num) {
  auto noise = 0.0f;
  for (auto index : range(num)) {
    noise += base_noise(position * pow(2, index)) / pow(2, index);
  }
  return noise / 2;
}

// Turbulence noise
template <typename Noise>
inline float turbulence_noise(Noise&& base_noise, vec2f position, int num) {
  auto noise = 0.0f;
  for (auto index : range(num)) {
    noise += abs(2 * base_noise(position * pow(2, index)) - 1) / pow(2, index);
  }
  return noise;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------

namespace yocto {

// Gets the cell index
inline vec3i get_cell_index(const hash_grid& grid, vec3f position) {
  auto scaledpos = position * grid.cell_inv_size;
  return vec3i{(int)scaledpos.x, (int)scaledpos.y, (int)scaledpos.z};
}

// Create a hash_grid
inline hash_grid make_hash_grid(float cell_size) {
  auto grid          = hash_grid{};
  grid.cell_size     = cell_size;
  grid.cell_inv_size = 1 / cell_size;
  return grid;
}
inline hash_grid make_hash_grid(
    const vector<vec3f>& positions, float cell_size) {
  auto grid          = hash_grid{};
  grid.cell_size     = cell_size;
  grid.cell_inv_size = 1 / cell_size;
  for (auto& position : positions) insert_vertex(grid, position);
  return grid;
}
// Inserts a point into the grid
inline int insert_vertex(hash_grid& grid, vec3f position) {
  auto vertex_id = (int)grid.positions.size();
  auto cell      = get_cell_index(grid, position);
  grid.cells[cell].push_back(vertex_id);
  grid.positions.push_back(position);
  return vertex_id;
}
// Finds the nearest neighbors within a given radius
inline void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    vec3f position, float max_radius, int skip_id) {
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
inline void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    vec3f position, float max_radius) {
  find_neighbors(grid, neighbors, position, max_radius, -1);
}
inline void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    int vertex, float max_radius) {
  find_neighbors(grid, neighbors, grid.positions[vertex], max_radius, vertex);
}

}  // namespace yocto

#endif
