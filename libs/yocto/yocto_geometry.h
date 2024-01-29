//
// # Yocto/Geometry: Geometry operations
//
// Yocto/Geometry defines basic geometry operations, including computation of
// basic geometry quantities, ray-primitive intersection, point-primitive
// distance, primitive bounds, and several interpolation functions.
// Yocto/Geometry is implemented in `yocto_geometry.h`.
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

#ifndef _YOCTO_GEOMETRY_H_
#define _YOCTO_GEOMETRY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <utility>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef __CUDACC__
#define inline inline __device__ __forceinline__
#endif

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2f {
  vec2f min = {flt_max, flt_max};
  vec2f max = {flt_min, flt_min};

  constexpr bbox2f() : min{flt_max, flt_max}, max{flt_min, flt_min} {}
  constexpr bbox2f(vec2f min, vec2f max) : min{min}, max{max} {}

  inline vec2f&       operator[](int i);
  inline const vec2f& operator[](int i) const;
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
  vec3f min = {flt_max, flt_max, flt_max};
  vec3f max = {flt_min, flt_min, flt_min};

  constexpr bbox3f()
      : min{flt_max, flt_max, flt_max}, max{flt_min, flt_min, flt_min} {}
  constexpr bbox3f(vec3f min, vec3f max) : min{min}, max{max} {}

  inline vec3f&       operator[](int i);
  inline const vec3f& operator[](int i) const;
};

// Empty bbox constant.
constexpr auto invalidb2f = bbox2f{};
constexpr auto invalidb3f = bbox3f{};

// Bounding box properties
inline vec2f center(const bbox2f& a);
inline vec2f size(const bbox2f& a);

// Bounding box tests
inline bool contains(const bbox2f& a, vec2f b);
inline bool contains(const bbox2f& a, const bbox2f& b);

// Bounding box comparisons.
inline bool operator==(const bbox2f& a, const bbox2f& b);
inline bool operator!=(const bbox2f& a, const bbox2f& b);

// Bounding box expansions with points and other boxes.
inline bbox2f merge(const bbox2f& a, vec2f b);
inline bbox2f merge(const bbox2f& a, const bbox2f& b);
inline void   expand(bbox2f& a, vec2f b);
inline void   expand(bbox2f& a, const bbox2f& b);

// Bounding box properties
inline vec3f center(const bbox3f& a);
inline vec3f size(const bbox3f& a);

// Bounding box tests
inline bool contains(const bbox3f& a, vec3f b);
inline bool contains(const bbox3f& a, const bbox3f& b);

// Bounding box comparisons.
inline bool operator==(const bbox3f& a, const bbox3f& b);
inline bool operator!=(const bbox3f& a, const bbox3f& b);

// Bounding box expansions with points and other boxes.
inline bbox3f merge(const bbox3f& a, vec3f b);
inline bbox3f merge(const bbox3f& a, const bbox3f& b);
inline void   expand(bbox3f& a, vec3f b);
inline void   expand(bbox3f& a, const bbox3f& b);

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Ray epsilon
constexpr auto ray_eps = 1e-4f;

struct ray2f {
  vec2f o    = {0, 0};
  vec2f d    = {0, 1};
  float tmin = ray_eps;
  float tmax = flt_max;

  constexpr ray2f() : o{0, 0}, d{0, 1}, tmin{ray_eps}, tmax{flt_max} {}
  constexpr ray2f(
      vec2f o_, vec2f d_, float tmin_ = ray_eps, float tmax_ = flt_max)
      : o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Rays with origin, direction and min/max t value.
struct ray3f {
  vec3f o    = {0, 0, 0};
  vec3f d    = {0, 0, 1};
  float tmin = ray_eps;
  float tmax = flt_max;

  constexpr ray3f() : o{0, 0, 0}, d{0, 0, 1}, tmin{ray_eps}, tmax{flt_max} {}
  constexpr ray3f(
      vec3f o_, vec3f d_, float tmin_ = ray_eps, float tmax_ = flt_max)
      : o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Computes a point on a ray
inline vec2f ray_point(const ray2f& ray, float t);
inline vec3f ray_point(const ray3f& ray, float t);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms rays.
inline ray3f transform_ray(const mat4f& a, const ray3f& b);
inline ray3f transform_ray(const frame3f& a, const ray3f& b);

// Transforms bounding boxes by matrices.
inline bbox3f transform_bbox(const mat4f& a, const bbox3f& b);
inline bbox3f transform_bbox(const frame3f& a, const bbox3f& b);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRIMITIVE BOUNDS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive bounds.
inline bbox3f point_bounds(vec3f p);
inline bbox3f point_bounds(vec3f p, float r);
inline bbox3f line_bounds(vec3f p0, vec3f p1);
inline bbox3f line_bounds(vec3f p0, vec3f p1, float r0, float r1);
inline bbox3f triangle_bounds(vec3f p0, vec3f p1, vec3f p2);
inline bbox3f quad_bounds(vec3f p0, vec3f p1, vec3f p2, vec3f p3);
inline bbox3f sphere_bounds(vec3f p, float r);
inline bbox3f capsule_bounds(vec3f p0, vec3f p1, float r0, float r1);

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
inline vec3f line_tangent(vec3f p0, vec3f p1);
inline float line_length(vec3f p0, vec3f p1);

// Triangle properties.
inline vec3f triangle_normal(vec3f p0, vec3f p1, vec3f p2);
inline float triangle_area(vec3f p0, vec3f p1, vec3f p2);

// Quad properties.
inline vec3f quad_normal(vec3f p0, vec3f p1, vec3f p2, vec3f p3);
inline float quad_area(vec3f p0, vec3f p1, vec3f p2, vec3f p3);

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
inline T interpolate_line(const T& p0, const T& p1, float u);

// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(const T& p0, const T& p1, const T& p2, vec2f uv);

// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, vec2f uv);

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, float u);

// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, float u);

// Interpolated line properties.
inline vec3f line_point(vec3f p0, vec3f p1, float u);
inline vec3f line_tangent(vec3f t0, vec3f t1, float u);

// Interpolated triangle properties.
inline vec3f triangle_point(vec3f p0, vec3f p1, vec3f p2, vec2f uv);
inline vec3f triangle_normal(vec3f n0, vec3f n1, vec3f n2, vec2f uv);

// Interpolated quad properties.
inline vec3f quad_point(vec3f p0, vec3f p1, vec3f p2, vec2f uv);
inline vec3f quad_normal(vec3f n0, vec3f n1, vec3f n2, vec3f n3, vec2f uv);

// Interpolated sphere properties.
inline vec3f sphere_point(const vec3f p, float r, vec2f uv);
inline vec3f sphere_normal(const vec3f p, float r, vec2f uv);

// Triangle tangent and bitangent from uv
inline pair<vec3f, vec3f> triangle_tangents_fromuv(
    vec3f p0, vec3f p1, vec3f p2, vec2f uv0, vec2f uv1, vec2f uv2);

// Quad tangent and bitangent from uv. Note that we pass a current_uv since
// internally we may want to split the quad in two and we need to known where
// to do it. If not interested in the split, just pass vec2f{0,0} here.
inline pair<vec3f, vec3f> quad_tangents_fromuv(vec3f p0, vec3f p1, vec3f p2,
    vec3f p3, vec2f uv0, vec2f uv1, vec2f uv2, vec2f uv3, vec2f current_uv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generate a ray from a camera
inline ray3f camera_ray(
    const frame3f& frame, float lens, vec2f film, vec2f image_uv);

// Generate a ray from a camera
inline ray3f camera_ray(
    const frame3f& frame, float lens, float aspect, float film, vec2f image_uv);

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive intersection
struct prim_intersection {
  vec2f uv       = {0, 0};
  float distance = flt_max;
  bool  hit      = false;
};

// Intersect a ray with a point (approximate)
inline prim_intersection intersect_point(const ray3f& ray, vec3f p, float r);

// Intersect a ray with a line
inline prim_intersection intersect_line(
    const ray3f& ray, vec3f p0, vec3f p1, float r0, float r1);

// Intersect a ray with a triangle
inline prim_intersection intersect_triangle(
    const ray3f& ray, vec3f p0, vec3f p1, vec3f p2);

// Intersect a ray with a quad.
inline prim_intersection intersect_quad(
    const ray3f& ray, vec3f p0, vec3f p1, vec3f p2, vec3f p3);

// Intersect a ray with a sphere
inline prim_intersection intersect_sphere(const ray3f& ray, vec3f p, float r);

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(
    const ray3f& ray, vec3f ray_dinv, const bbox3f& bbox);

}  // namespace yocto

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_point(
    vec3f pos, float dist_max, vec3f p, float r);

// Compute the closest line uv to a give position pos.
inline float closestuv_line(vec3f pos, vec3f p0, vec3f p1);

// Check if a line overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_line(
    vec3f pos, float dist_max, vec3f p0, vec3f p1, float r0, float r1);

// Compute the closest triangle uv to a give position pos.
inline vec2f closestuv_triangle(vec3f pos, vec3f p0, vec3f p1, vec3f p2);

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
inline prim_intersection overlap_triangle(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, float r0, float r1, float r2);

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_quad(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, vec3f p3, float r0, float r1, float r2, float r3);

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
inline bool overlap_bbox(vec3f pos, float dist_max, const bbox3f& bbox);

// Check if two bboxes overlap.
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex normals/tangents for lines/triangles/quads.
inline vector<vec3f> lines_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions);
inline vector<vec3f> triangles_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
inline vector<vec3f> quads_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions);
// Update normals and tangents
inline void lines_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions);
inline void triangles_normals(vector<vec3f>& normals,
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
inline void quads_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
inline vector<vec4f> triangle_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
inline vector<vec3f> flip_normals(const vector<vec3f>& normals);
// Flip face orientation
inline vector<vec3i> flip_triangles(const vector<vec3i>& triangles);
inline vector<vec4i> flip_quads(const vector<vec4i>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
inline vector<vec3f> align_vertices(
    const vector<vec3f>& positions, vec3i alignment);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE ELEMENT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
inline vector<vec3i> quads_to_triangles(const vector<vec4i>& quads);
// Convert triangles to quads by creating degenerate quads
inline vector<vec4i> triangles_to_quads(const vector<vec3i>& triangles);
// Convert beziers to lines using 3 lines for each bezier.
inline vector<vec4i> bezier_to_lines(vector<vec2i>& lines);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Axis aligned bounding box represented as a min/max vector pairs.
inline vec2f&       bbox2f::operator[](int i) { return (&min)[i]; }
inline const vec2f& bbox2f::operator[](int i) const { return (&min)[i]; }

// Axis aligned bounding box represented as a min/max vector pairs.
inline vec3f&       bbox3f::operator[](int i) { return (&min)[i]; }
inline const vec3f& bbox3f::operator[](int i) const { return (&min)[i]; }

// Bounding box properties
inline vec2f center(const bbox2f& a) { return (a.min + a.max) / 2; }
inline vec2f size(const bbox2f& a) { return a.max - a.min; }

// Bounding box tests
inline bool contains(const bbox2f& a, vec2f b) {
  return b.x >= a.min.x && b.x <= a.max.x && b.y >= a.min.y && b.y <= a.max.y;
}
inline bool contains(const bbox2f& a, const bbox2f& b) {
  return contains(a, b.min) && contains(a, b.max);
}

// Bounding box comparisons.
inline bool operator==(const bbox2f& a, const bbox2f& b) {
  return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox2f& a, const bbox2f& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
inline bbox2f merge(const bbox2f& a, vec2f b) {
  return {min(a.min, b), max(a.max, b)};
}
inline bbox2f merge(const bbox2f& a, const bbox2f& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
inline void expand(bbox2f& a, vec2f b) { a = merge(a, b); }
inline void expand(bbox2f& a, const bbox2f& b) { a = merge(a, b); }

// Bounding box properties
inline vec3f center(const bbox3f& a) { return (a.min + a.max) / 2; }
inline vec3f size(const bbox3f& a) { return a.max - a.min; }

// Bounding box tests
inline bool contains(const bbox3f& a, vec3f b) {
  return b.x >= a.min.x && b.x <= a.max.x && b.y >= a.min.y && b.y <= a.max.y &&
         b.z >= a.min.z && b.z <= a.max.z;
}
inline bool contains(const bbox3f& a, const bbox3f& b) {
  return contains(a, b.min) && contains(a, b.max);
}

// Bounding box comparisons.
inline bool operator==(const bbox3f& a, const bbox3f& b) {
  return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox3f& a, const bbox3f& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
inline bbox3f merge(const bbox3f& a, vec3f b) {
  return {min(a.min, b), max(a.max, b)};
}
inline bbox3f merge(const bbox3f& a, const bbox3f& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
inline void expand(bbox3f& a, vec3f b) { a = merge(a, b); }
inline void expand(bbox3f& a, const bbox3f& b) { a = merge(a, b); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Computes a point on a ray
inline vec2f ray_point(const ray2f& ray, float t) { return ray.o + ray.d * t; }
inline vec3f ray_point(const ray3f& ray, float t) { return ray.o + ray.d * t; }

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms rays and bounding boxes by matrices.
inline ray3f transform_ray(const mat4f& a, const ray3f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
inline ray3f transform_ray(const frame3f& a, const ray3f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
inline bbox3f transform_bbox(const mat4f& a, const bbox3f& b) {
  auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
      vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
      vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
      vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
      vec3f{b.max.x, b.max.y, b.max.z}};
  auto xformed = bbox3f();
  for (auto& corner : corners)
    xformed = merge(xformed, transform_point(a, corner));
  return xformed;
}
inline bbox3f transform_bbox(const frame3f& a, const bbox3f& b) {
  auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
      vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
      vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
      vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
      vec3f{b.max.x, b.max.y, b.max.z}};
  auto xformed = bbox3f();
  for (auto& corner : corners)
    xformed = merge(xformed, transform_point(a, corner));
  return xformed;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRIMITIVE BOUNDS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive bounds.
inline bbox3f point_bounds(vec3f p) { return {p, p}; }
inline bbox3f point_bounds(vec3f p, float r) {
  return {min(p - r, p + r), max(p - r, p + r)};
}
inline bbox3f line_bounds(vec3f p0, vec3f p1) {
  return {min(p0, p1), max(p0, p1)};
}
inline bbox3f line_bounds(vec3f p0, vec3f p1, float r0, float r1) {
  return {min(p0 - r0, p1 - r1), max(p0 + r0, p1 + r1)};
}
inline bbox3f triangle_bounds(vec3f p0, vec3f p1, vec3f p2) {
  return {min(p0, min(p1, p2)), max(p0, max(p1, p2))};
}
inline bbox3f quad_bounds(vec3f p0, vec3f p1, vec3f p2, vec3f p3) {
  return {min(p0, min(p1, min(p2, p3))), max(p0, max(p1, max(p2, p3)))};
}
inline bbox3f sphere_bounds(vec3f p, float r) { return {p - r, p + r}; }
inline bbox3f capsule_bounds(vec3f p0, vec3f p1, float r0, float r1) {
  return {min(p0 - r0, p1 - r1), max(p0 + r0, p1 + r1)};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
inline vec3f line_tangent(vec3f p0, vec3f p1) { return normalize(p1 - p0); }
inline float line_length(vec3f p0, vec3f p1) { return length(p1 - p0); }

// Triangle properties.
inline vec3f triangle_normal(vec3f p0, vec3f p1, vec3f p2) {
  return normalize(cross(p1 - p0, p2 - p0));
}
inline float triangle_area(vec3f p0, vec3f p1, vec3f p2) {
  return length(cross(p1 - p0, p2 - p0)) / 2;
}

// Quad propeties.
inline vec3f quad_normal(vec3f p0, vec3f p1, vec3f p2, vec3f p3) {
  return normalize(triangle_normal(p0, p1, p3) + triangle_normal(p2, p3, p1));
}
inline float quad_area(vec3f p0, vec3f p1, vec3f p2, vec3f p3) {
  return triangle_area(p0, p1, p3) + triangle_area(p2, p3, p1);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
inline T interpolate_line(const T& p0, const T& p1, float u) {
  return p0 * (1 - u) + p1 * u;
}

// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(const T& p0, const T& p1, const T& p2, vec2f uv) {
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}

// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, vec2f uv) {
  if (uv.x + uv.y <= 1) {
    return interpolate_triangle(p0, p1, p3, uv);
  } else {
    return interpolate_triangle(p2, p3, p1, 1 - uv);
  }
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
  return p0 * (1 - u) * (1 - u) * (1 - u) + p1 * 3 * u * (1 - u) * (1 - u) +
         p2 * 3 * u * u * (1 - u) + p3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
  return (p1 - p0) * 3 * (1 - u) * (1 - u) + (p2 - p1) * 6 * u * (1 - u) +
         (p3 - p2) * 3 * u * u;
}

// Interpolated line properties.
inline vec3f line_point(vec3f p0, vec3f p1, float u) {
  return p0 * (1 - u) + p1 * u;
}
inline vec3f line_tangent(vec3f t0, vec3f t1, float u) {
  return normalize(t0 * (1 - u) + t1 * u);
}

// Interpolated triangle properties.
inline vec3f triangle_point(vec3f p0, vec3f p1, vec3f p2, vec2f uv) {
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
inline vec3f triangle_normal(vec3f n0, vec3f n1, vec3f n2, vec2f uv) {
  return normalize(n0 * (1 - uv.x - uv.y) + n1 * uv.x + n2 * uv.y);
}

// Interpolated quad properties.
inline vec3f quad_point(vec3f p0, vec3f p1, vec3f p2, vec3f p3, vec2f uv) {
  if (uv.x + uv.y <= 1) {
    return triangle_point(p0, p1, p3, uv);
  } else {
    return triangle_point(p2, p3, p1, 1 - uv);
  }
}
inline vec3f quad_normal(vec3f n0, vec3f n1, vec3f n2, vec3f n3, vec2f uv) {
  if (uv.x + uv.y <= 1) {
    return triangle_normal(n0, n1, n3, uv);
  } else {
    return triangle_normal(n2, n3, n1, 1 - uv);
  }
}

// Interpolated sphere properties.
inline vec3f sphere_point(const vec3f p, float r, vec2f uv) {
  return p + r * vec3f{cos(uv.x * 2 * pif) * sin(uv.y * pif),
                     sin(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif)};
}
inline vec3f sphere_normal(const vec3f p, float r, vec2f uv) {
  return normalize(vec3f{cos(uv.x * 2 * pif) * sin(uv.y * pif),
      sin(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif)});
}

// Triangle tangent and bitangent from uv
inline pair<vec3f, vec3f> triangle_tangents_fromuv(
    vec3f p0, vec3f p1, vec3f p2, vec2f uv0, vec2f uv1, vec2f uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  auto p   = p1 - p0;
  auto q   = p2 - p0;
  auto s   = vec2f{uv1.x - uv0.x, uv2.x - uv0.x};
  auto t   = vec2f{uv1.y - uv0.y, uv2.y - uv0.y};
  auto div = s.x * t.y - s.y * t.x;

  if (div != 0) {
    auto tu = vec3f{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                  t.y * p.z - t.x * q.z} /
              div;
    auto tv = vec3f{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                  s.x * q.z - s.y * p.z} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}

// Quad tangent and bitangent from uv.
inline pair<vec3f, vec3f> quad_tangents_fromuv(vec3f p0, vec3f p1, vec3f p2,
    vec3f p3, vec2f uv0, vec2f uv1, vec2f uv2, vec2f uv3, vec2f current_uv) {
  if (current_uv.x + current_uv.y <= 1) {
    return triangle_tangents_fromuv(p0, p1, p3, uv0, uv1, uv3);
  } else {
    return triangle_tangents_fromuv(p2, p3, p1, uv2, uv3, uv1);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generate a ray from a camera
inline ray3f camera_ray(
    const frame3f& frame, float lens, vec2f film, vec2f image_uv) {
  auto e = vec3f{0, 0, 0};
  auto q = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), lens};
  auto q1  = -q;
  auto d   = normalize(q1 - e);
  auto ray = ray3f{transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

// Generate a ray from a camera
inline ray3f camera_ray(const frame3f& frame, float lens, float aspect,
    float film_, vec2f image_uv) {
  auto film = aspect >= 1 ? vec2f{film_, film_ / aspect}
                          : vec2f{film_ * aspect, film_};
  auto e    = vec3f{0, 0, 0};
  auto q    = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), lens};
  auto q1  = -q;
  auto d   = normalize(q1 - e);
  auto ray = ray3f{transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate)
inline prim_intersection intersect_point(const ray3f& ray, vec3f p, float r) {
  // find parameter for line-point minimum distance
  auto w = p - ray.o;
  auto t = dot(w, ray.d) / dot(ray.d, ray.d);

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return {};

  // test for line-point distance vs point radius
  auto rp  = ray.o + ray.d * t;
  auto prp = p - rp;
  if (dot(prp, prp) > r * r) return {};

  // intersection occurred: set params and exit
  return {{0, 0}, t, true};
}

// Intersect a ray with a line
inline prim_intersection intersect_line(
    const ray3f& ray, vec3f p0, vec3f p1, float r0, float r1) {
  // setup intersection params
  auto u = ray.d;
  auto v = p1 - p0;
  auto w = ray.o - p0;

  // compute values to solve a linear system
  auto a   = dot(u, u);
  auto b   = dot(u, v);
  auto c   = dot(v, v);
  auto d   = dot(u, w);
  auto e   = dot(v, w);
  auto det = a * c - b * b;

  // check determinant and exit if lines are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return {};

  // compute Parameters on both ray and segment
  auto t = (b * e - c * d) / det;
  auto s = (a * e - b * d) / det;

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return {};

  // clamp segment param to segment corners
  s = clamp(s, (float)0, (float)1);

  // compute segment-segment distance on the closest points
  auto pr  = ray.o + ray.d * t;
  auto pl  = p0 + (p1 - p0) * s;
  auto prl = pr - pl;

  // check with the line radius at the same point
  auto d2 = dot(prl, prl);
  auto r  = r0 * (1 - s) + r1 * s;
  if (d2 > r * r) return {};

  // intersection occurred: set params and exit
  return {{s, sqrt(d2) / r}, t, true};
}

// Intersect a ray with a triangle
inline prim_intersection intersect_triangle(
    const ray3f& ray, vec3f p0, vec3f p1, vec3f p2) {
  // compute triangle edges
  auto edge1 = p1 - p0;
  auto edge2 = p2 - p0;

  // compute determinant to solve a linear system
  auto pvec = cross(ray.d, edge2);
  auto det  = dot(edge1, pvec);

  // check determinant and exit if triangle and ray are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return {};
  auto inv_det = 1.0f / det;

  // compute and check first bricentric coordinated
  auto tvec = ray.o - p0;
  auto u    = dot(tvec, pvec) * inv_det;
  if (u < 0 || u > 1) return {};

  // compute and check second bricentric coordinated
  auto qvec = cross(tvec, edge1);
  auto v    = dot(ray.d, qvec) * inv_det;
  if (v < 0 || u + v > 1) return {};

  // compute and check ray parameter
  auto t = dot(edge2, qvec) * inv_det;
  if (t < ray.tmin || t > ray.tmax) return {};

  // intersection occurred: set params and exit
  return {{u, v}, t, true};
}

// Intersect a ray with a quad.
inline prim_intersection intersect_quad(
    const ray3f& ray, vec3f p0, vec3f p1, vec3f p2, vec3f p3) {
  if (p2 == p3) return intersect_triangle(ray, p0, p1, p3);
  auto isec1 = intersect_triangle(ray, p0, p1, p3);
  auto isec2 = intersect_triangle(ray, p2, p3, p1);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}

// Intersect a ray with a sphere
inline prim_intersection intersect_sphere(const ray3f& ray, vec3f p, float r) {
  // compute parameters
  auto a = dot(ray.d, ray.d);
  auto b = 2 * dot(ray.o - p, ray.d);
  auto c = dot(ray.o - p, ray.o - p) - r * r;

  // check discriminant
  auto dis = b * b - 4 * a * c;
  if (dis < 0) return {};

  // compute ray parameter
  auto t = (-b - sqrt(dis)) / (2 * a);

  // try other ray parameter
  if (t < ray.tmin || t > ray.tmax) t = (-b + sqrt(dis)) / (2 * a);

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return {};

  // compute local point for uvs
  auto plocal = ((ray.o + ray.d * t) - p) / r;
  auto u      = atan2(plocal.y, plocal.x) / (2 * pif);
  if (u < 0) u += 1;
  auto v = acos(clamp(plocal.z, -1.0f, 1.0f)) / pif;

  // intersection occurred: set params and exit
  return {{u, v}, t, true};
}

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(const ray3f& ray, const bbox3f& bbox) {
  // determine intersection ranges
  auto invd = 1.0f / ray.d;
  auto t0   = (bbox.min - ray.o) * invd;
  auto t1   = (bbox.max - ray.o) * invd;
  // flip based on range directions
  if (invd.x < 0.0f) swap(t0.x, t1.x);
  if (invd.y < 0.0f) swap(t0.y, t1.y);
  if (invd.z < 0.0f) swap(t0.z, t1.z);
  auto tmin = max(t0.z, max(t0.y, max(t0.x, ray.tmin)));
  auto tmax = min(t1.z, min(t1.y, min(t1.x, ray.tmax)));
  tmax *= 1.00000024f;  // for double: 1.0000000000000004
  return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(
    const ray3f& ray, vec3f ray_dinv, const bbox3f& bbox) {
  auto it_min = (bbox.min - ray.o) * ray_dinv;
  auto it_max = (bbox.max - ray.o) * ray_dinv;
  auto tmin   = min(it_min, it_max);
  auto tmax   = max(it_min, it_max);
  auto t0     = max(max(tmin), ray.tmin);
  auto t1     = min(min(tmax), ray.tmax);
  t1 *= 1.00000024f;  // for double: 1.0000000000000004
  return t0 <= t1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_point(
    vec3f pos, float dist_max, vec3f p, float r) {
  auto d2 = dot(pos - p, pos - p);
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  return {{0, 0}, sqrt(d2), true};
}

// Compute the closest line uv to a give position pos.
inline float closestuv_line(vec3f pos, vec3f p0, vec3f p1) {
  auto ab = p1 - p0;
  auto d  = dot(ab, ab);
  // Project c onto ab, computing parameterized position d(t) = a + t*(b –
  // a)
  auto u = dot(pos - p0, ab) / d;
  u      = clamp(u, (float)0, (float)1);
  return u;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_line(
    vec3f pos, float dist_max, vec3f p0, vec3f p1, float r0, float r1) {
  auto u = closestuv_line(pos, p0, p1);
  // Compute projected position from the clamped t d = a + t * ab;
  auto p  = p0 + (p1 - p0) * u;
  auto r  = r0 + (r1 - r0) * u;
  auto d2 = dot(pos - p, pos - p);
  // check distance
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  // done
  return {{u, 0}, sqrt(d2), true};
}

// Compute the closest triangle uv to a give position pos.
inline vec2f closestuv_triangle(vec3f pos, vec3f p0, vec3f p1, vec3f p2) {
  // this is a complicated test -> I probably "--"+prefix to use a sequence of
  // test (triangle body, and 3 edges)
  auto ab = p1 - p0;
  auto ac = p2 - p0;
  auto ap = pos - p0;

  auto d1 = dot(ab, ap);
  auto d2 = dot(ac, ap);

  // corner and edge cases
  if (d1 <= 0 && d2 <= 0) return {0, 0};

  auto bp = pos - p1;
  auto d3 = dot(ab, bp);
  auto d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return {1, 0};

  auto vc = d1 * d4 - d3 * d2;
  if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return {d1 / (d1 - d3), 0};

  auto cp = pos - p2;
  auto d5 = dot(ab, cp);
  auto d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) return {0, 1};

  auto vb = d5 * d2 - d1 * d6;
  if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) return {0, d2 / (d2 - d6)};

  auto va = d3 * d6 - d5 * d4;
  if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
    auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return {1 - w, w};
  }

  // face case
  auto denom = 1 / (va + vb + vc);
  auto u     = vb * denom;
  auto v     = vc * denom;
  return {u, v};
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
inline prim_intersection overlap_triangle(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, float r0, float r1, float r2) {
  auto cuv = closestuv_triangle(pos, p0, p1, p2);
  auto p   = p0 * (1 - cuv.x - cuv.y) + p1 * cuv.x + p2 * cuv.y;
  auto r   = r0 * (1 - cuv.x - cuv.y) + r1 * cuv.x + r2 * cuv.y;
  auto dd  = dot(p - pos, p - pos);
  if (dd > (dist_max + r) * (dist_max + r)) return {};
  return {cuv, sqrt(dd), true};
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
inline prim_intersection overlap_quad(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, vec3f p3, float r0, float r1, float r2, float r3) {
  if (p2 == p3) return overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r2);
  auto isec1 = overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r2);
  auto isec2 = overlap_triangle(pos, dist_max, p2, p3, p1, r2, r3, r1);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
inline bool overlap_bbox(vec3f pos, float dist_max, const bbox3f& bbox) {
  // computing distance
  auto dd = 0.0f;

  // For each axis count any excess distance outside box extents
  if (pos.x < bbox.min.x) dd += (bbox.min.x - pos.x) * (bbox.min.x - pos.x);
  if (pos.x > bbox.max.x) dd += (pos.x - bbox.max.x) * (pos.x - bbox.max.x);
  if (pos.y < bbox.min.y) dd += (bbox.min.y - pos.y) * (bbox.min.y - pos.y);
  if (pos.y > bbox.max.y) dd += (pos.y - bbox.max.y) * (pos.y - bbox.max.y);
  if (pos.z < bbox.min.z) dd += (bbox.min.z - pos.z) * (bbox.min.z - pos.z);
  if (pos.z > bbox.max.z) dd += (pos.z - bbox.max.z) * (pos.z - bbox.max.z);

  // check distance
  return dd < dist_max * dist_max;
}

// Check if two bboxe overlap.
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
  if (bbox1.max.x < bbox2.min.x || bbox1.min.x > bbox2.max.x) return false;
  if (bbox1.max.y < bbox2.min.y || bbox1.min.y > bbox2.max.y) return false;
  if (bbox1.max.z < bbox2.min.z || bbox1.min.z > bbox2.max.z) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARD COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate)
[[deprecated]] inline bool intersect_point(
    const ray3f& ray, vec3f p, float r, vec2f& uv, float& dist) {
  auto intersection = intersect_point(ray, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a line
[[deprecated]] inline bool intersect_line(const ray3f& ray, vec3f p0, vec3f p1,
    float r0, float r1, vec2f& uv, float& dist) {
  auto intersection = intersect_line(ray, p0, p1, r0, r1);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a sphere
[[deprecated]] inline bool intersect_sphere(
    const ray3f& ray, vec3f p, float r, vec2f& uv, float& dist) {
  auto intersection = intersect_sphere(ray, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a triangle
[[deprecated]] inline bool intersect_triangle(
    const ray3f& ray, vec3f p0, vec3f p1, vec3f p2, vec2f& uv, float& dist) {
  auto intersection = intersect_triangle(ray, p0, p1, p2);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a quad.
[[deprecated]] inline bool intersect_quad(const ray3f& ray, vec3f p0, vec3f p1,
    vec3f p2, vec3f p3, vec2f& uv, float& dist) {
  auto intersection = intersect_quad(ray, p0, p1, p2, p3);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a point overlaps a position pos withint a maximum distance dist_max.
[[deprecated]] inline bool overlap_point(
    vec3f pos, float dist_max, vec3f p, float r, vec2f& uv, float& dist) {
  auto intersection = overlap_point(pos, dist_max, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
[[deprecated]] inline bool overlap_line(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, float r0, float r1, vec2f& uv, float& dist) {
  auto intersection = overlap_line(pos, dist_max, p0, p1, r0, r1);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
[[deprecated]] inline bool overlap_triangle(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, float r0, float r1, float r2, vec2f& uv, float& dist) {
  auto intersection = overlap_triangle(pos, dist_max, p0, p1, p2, r0, r1, r2);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
[[deprecated]] inline bool overlap_quad(vec3f pos, float dist_max, vec3f p0,
    vec3f p1, vec3f p2, vec3f p3, float r0, float r1, float r2, float r3,
    vec2f& uv, float& dist) {
  auto intersection = overlap_quad(
      pos, dist_max, p0, p1, p2, p3, r0, r1, r2, r3);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMPUTATION OF PER-VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex tangents for lines.
inline vector<vec3f> lines_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto tangents = vector<vec3f>{positions.size()};
  for (auto& tangent : tangents) tangent = {0, 0, 0};
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
inline vector<vec3f> triangles_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
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
inline vector<vec3f> quads_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
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
inline void lines_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions) {
  if (tangents.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& tangent : tangents) tangent = {0, 0, 0};
  for (auto& l : lines) {
    auto tangent = line_tangent(positions[l.x], positions[l.y]);
    auto length  = line_length(positions[l.x], positions[l.y]);
    tangents[l.x] += tangent * length;
    tangents[l.y] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
}

// Compute per-vertex normals for triangles.
inline void triangles_normals(vector<vec3f>& normals,
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = {0, 0, 0};
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
inline void quads_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions) {
  if (normals.size() != positions.size()) {
    throw std::out_of_range("array should be the same length");
  }
  for (auto& normal : normals) normal = {0, 0, 0};
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
inline vector<vec4f> triangles_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  auto tangu = vector<vec3f>(positions.size(), vec3f{0, 0, 0});
  auto tangv = vector<vec3f>(positions.size(), vec3f{0, 0, 0});
  for (auto t : triangles) {
    auto tutv = triangle_tangents_fromuv(positions[t.x], positions[t.y],
        positions[t.z], texcoords[t.x], texcoords[t.y], texcoords[t.z]);
    for (auto vid : {t.x, t.y, t.z}) tangu[vid] += normalize(tutv.first);
    for (auto vid : {t.x, t.y, t.z}) tangv[vid] += normalize(tutv.second);
  }
  for (auto& t : tangu) t = normalize(t);
  for (auto& t : tangv) t = normalize(t);

  auto tangent_spaces = vector<vec4f>(positions.size());
  for (auto& tangent : tangent_spaces) tangent = vec4f{0, 0, 0, 0};
  for (auto i : range(positions.size())) {
    tangu[i] = orthonormalize(tangu[i], normals[i]);
    auto s   = (dot(cross(normals[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
    tangent_spaces[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
  }
  return tangent_spaces;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
inline vector<vec3f> flip_normals(const vector<vec3f>& normals) {
  auto flipped = normals;
  for (auto& n : flipped) n = -n;
  return flipped;
}
// Flip face orientation
inline vector<vec3i> flip_triangles(const vector<vec3i>& triangles) {
  auto flipped = triangles;
  for (auto& t : flipped) swap(t.y, t.z);
  return flipped;
}
inline vector<vec4i> flip_quads(const vector<vec4i>& quads) {
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
inline vector<vec3f> align_vertices(
    const vector<vec3f>& positions, vec3i alignment) {
  auto bounds = invalidb3f;
  for (auto& p : positions) bounds = merge(bounds, p);
  auto offset = vec3f{0, 0, 0};
  switch (alignment.x) {
    case 0: break;
    case 1: offset.x = bounds.min.x; break;
    case 2: offset.x = (bounds.min.x + bounds.max.x) / 2; break;
    case 3: offset.x = bounds.max.x; break;
    default: throw std::invalid_argument{"invalid alignment"};
  }
  switch (alignment.y) {
    case 0: break;
    case 1: offset.y = bounds.min.y; break;
    case 2: offset.y = (bounds.min.y + bounds.max.y) / 2; break;
    case 3: offset.y = bounds.max.y; break;
    default: throw std::invalid_argument{"invalid alignment"};
  }
  switch (alignment.z) {
    case 0: break;
    case 1: offset.z = bounds.min.z; break;
    case 2: offset.z = (bounds.min.z + bounds.max.z) / 2; break;
    case 3: offset.z = bounds.max.z; break;
    default: throw std::invalid_argument{"invalid alignment"};
  }
  auto aligned = positions;
  for (auto& p : aligned) p -= offset;
  return aligned;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE ELEMENT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
inline vector<vec3i> quads_to_triangles(const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  return triangles;
}

// Convert triangles to quads by creating degenerate quads
inline vector<vec4i> triangles_to_quads(const vector<vec3i>& triangles) {
  auto quads = vector<vec4i>{};
  quads.reserve(triangles.size());
  for (auto& t : triangles) quads.push_back({t.x, t.y, t.z, t.z});
  return quads;
}

// Convert beziers to lines using 3 lines for each bezier.
inline vector<vec2i> bezier_to_lines(const vector<vec4i>& beziers) {
  auto lines = vector<vec2i>{};
  lines.reserve(beziers.size() * 3);
  for (auto b : beziers) {
    lines.push_back({b.x, b.y});
    lines.push_back({b.y, b.z});
    lines.push_back({b.z, b.w});
  }
  return lines;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef __CUDACC__
#undef inline
#endif

#endif
