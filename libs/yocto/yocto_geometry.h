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

#ifndef YOCTO_GEOMETRY_H_
#define YOCTO_GEOMETRY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <utility>
#include <vector>

#include "yocto_math.h"
#include "yocto_views.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifndef kernel
#ifdef __CUDACC__
#define kernel __device__
#else
#define kernel
#endif
#endif

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2f {
  vec2f min = {flt_max, flt_max};
  vec2f max = {flt_min, flt_min};

  constexpr kernel bbox2f() : min{flt_max, flt_max}, max{flt_min, flt_min} {}
  constexpr kernel bbox2f(const vec2f& min_, const vec2f& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec2f& operator[](size_t i) { return i == 0 ? min : max; }
  constexpr kernel const vec2f& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
  vec3f min = {flt_max, flt_max, flt_max};
  vec3f max = {flt_min, flt_min, flt_min};

  constexpr kernel bbox3f() :
      min{flt_max, flt_max, flt_max}, max{flt_min, flt_min, flt_min} {}
  constexpr kernel bbox3f(const vec3f& min_, const vec3f& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec3f& operator[](size_t i) { return i == 0 ? min : max; }
  constexpr kernel const vec3f& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Empty bbox constant.
constexpr auto invalidb2f = bbox2f{};
constexpr auto invalidb3f = bbox3f{};

// Bounding box properties
constexpr kernel bool valid(const bbox2f& a) {
  return a.min.x <= a.max.x && a.min.y <= a.max.y;
}
constexpr kernel bool  empty(const bbox2f& a) { return a.min == a.max; }
constexpr kernel vec2f center(const bbox2f& a) { return (a.min + a.max) / 2; }
constexpr kernel vec2f diagonal(const bbox2f& a) { return a.max - a.min; }
constexpr kernel float area(const bbox2f& a) {
  auto d = a.max - a.min;
  return 2 * d.x * d.y;
}

// Bounding box comparisons.
constexpr kernel bool operator==(const bbox2f& a, const bbox2f& b) {
  return a.min == b.min && a.max == b.max;
}
constexpr kernel bool operator!=(const bbox2f& a, const bbox2f& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
constexpr kernel bbox2f merge(const bbox2f& a, const vec2f& b) {
  return {min(a.min, b), max(a.max, b)};
}
constexpr kernel bbox2f merge(const bbox2f& a, const bbox2f& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
constexpr kernel void expand(bbox2f& a, const vec2f& b) { a = merge(a, b); }
constexpr kernel void expand(bbox2f& a, const bbox2f& b) { a = merge(a, b); }

// Bounding box properties
constexpr kernel bool valid(const bbox3f& a) {
  return a.min.x <= a.max.x && a.min.y <= a.max.y && a.min.z <= a.max.z;
}
constexpr kernel bool  empty(const bbox3f& a) { return a.min == a.max; }
constexpr kernel vec3f center(const bbox3f& a) { return (a.min + a.max) / 2; }
constexpr kernel vec3f diagonal(const bbox3f& a) { return a.max - a.min; }
constexpr kernel float area(const bbox3f& a) {
  auto d = a.max - a.min;
  return 2 * d.x * d.y + 2 * d.x * d.z + 2 * d.y * d.z;
}

// Bounding box comparisons.
constexpr kernel bool operator==(const bbox3f& a, const bbox3f& b) {
  return a.min == b.min && a.max == b.max;
}
constexpr kernel bool operator!=(const bbox3f& a, const bbox3f& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
constexpr kernel bbox3f merge(const bbox3f& a, const vec3f& b) {
  return {min(a.min, b), max(a.max, b)};
}
constexpr kernel bbox3f merge(const bbox3f& a, const bbox3f& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
constexpr kernel void expand(bbox3f& a, const vec3f& b) { a = merge(a, b); }
constexpr kernel void expand(bbox3f& a, const bbox3f& b) { a = merge(a, b); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Ray epsilon
constexpr float ray_eps = 1e-4f;

// Rays with origin, direction and min/max t value.
struct ray2f {
  vec2f o    = {0, 0};
  vec2f d    = {0, 1};
  float tmin = ray_eps;
  float tmax = flt_max;

  constexpr kernel ray2f() : o{0, 0}, d{0, 1}, tmin{ray_eps}, tmax{flt_max} {}
  constexpr kernel ray2f(const vec2f& o_, const vec2f& d_,
      float tmin_ = ray_eps, float tmax_ = flt_max) :
      o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Rays with origin, direction and min/max t value.
struct ray3f {
  vec3f o    = {0, 0, 0};
  vec3f d    = {0, 0, 1};
  float tmin = ray_eps;
  float tmax = flt_max;

  constexpr kernel ray3f() : o{0}, d{0}, tmin{ray_eps}, tmax{flt_max} {}
  constexpr kernel ray3f(const vec3f& o_, const vec3f& d_,
      float tmin_ = ray_eps, float tmax_ = flt_max) :
      o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Computes a point on a ray
constexpr kernel vec2f ray_point(const ray2f& ray, float t) {
  return ray.o + ray.d * t;
}
constexpr kernel vec2f eval_ray(const ray2f& ray, float t) {
  return ray.o + ray.d * t;
}

// Computes a point on a ray
constexpr kernel vec3f ray_point(const ray3f& ray, float t) {
  return ray.o + ray.d * t;
}
constexpr kernel vec3f eval_ray(const ray3f& ray, float t) {
  return ray.o + ray.d * t;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms rays.
constexpr kernel ray2f transform_ray(const mat3f& a, const ray2f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
constexpr kernel ray2f transform_ray(const frame2f& a, const ray2f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}

// Transforms rays.
constexpr kernel ray3f transform_ray(const mat4f& a, const ray3f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
constexpr kernel ray3f transform_ray(const frame3f& a, const ray3f& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}

// Transforms bboxes.
constexpr kernel bbox2f transform_bbox(const mat3f& a, const bbox2f& b) {
  auto xformed = bbox2f{};
  for (auto ij : range(vec2i(2, 2))) {
    auto corner = vec2f{b[ij.x][0], b[ij.y][1]};
    xformed     = merge(xformed, transform_point(a, corner));
  }
  return xformed;
}
constexpr kernel bbox2f transform_bbox(const frame2f& a, const bbox2f& b) {
  auto xformed = bbox2f{};
  for (auto ij : range(vec2i(2, 2))) {
    auto corner = vec2f{b[ij.x][0], b[ij.y][1]};
    xformed     = merge(xformed, transform_point(a, corner));
  }
  return xformed;
}

constexpr kernel bbox3f transform_bbox(const mat4f& a, const bbox3f& b) {
  auto xformed = bbox3f{};
  for (auto ijk : range(vec3i(2, 2, 2))) {
    auto corner = vec3f{b[ijk.x][0], b[ijk.y][1], b[ijk.z][2]};
    xformed     = merge(xformed, transform_point(a, corner));
  }
  return xformed;
}
constexpr kernel bbox3f transform_bbox(const frame3f& a, const bbox3f& b) {
  auto xformed = bbox3f{};
  for (auto ijk : range(vec3i(2, 2, 2))) {
    auto corner = vec3f{b[ijk.x][0], b[ijk.y][1], b[ijk.z][2]};
    xformed     = merge(xformed, transform_point(a, corner));
  }
  return xformed;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUAD CONVENTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a quad is a triangle
constexpr kernel bool is_triangle(const vec4i& quad) {
  return quad.z == quad.w;
}

// Get the triangle from the quad
constexpr kernel vec3i as_triangle(const vec4i& quad) { return xyz(quad); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRIMITIVE BOUNDS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive bounds.
constexpr kernel bbox3f point_bounds(const vec3f& p) { return {p, p}; }
constexpr kernel bbox3f point_bounds(const vec3f& p, float r) {
  return {min(p - r, p + r), max(p - r, p + r)};
}
constexpr kernel bbox3f line_bounds(const vec3f& p1, const vec3f& p2) {
  return {min(p1, p2), max(p1, p2)};
}
constexpr kernel bbox3f line_bounds(
    const vec3f& p1, const vec3f& p2, float r1, float r2) {
  return {min(p1 - r1, p2 - r2), max(p1 + r1, p2 + r2)};
}
constexpr kernel bbox3f triangle_bounds(
    const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  return {min(p1, min(p2, p3)), max(p1, max(p2, p3))};
}
constexpr kernel bbox3f quad_bounds(
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec3f& p4) {
  return {min(p1, min(p2, min(p3, p4))), max(p1, max(p2, max(p3, p4)))};
}
constexpr kernel bbox3f sphere_bounds(const vec3f& p, float r) {
  return {p - r, p + r};
}
constexpr kernel bbox3f capsule_bounds(
    const vec3f& p1, const vec3f& p2, float r1, float r2) {
  return {min(p1 - r1, p2 - r2), max(p1 + r1, p2 + r2)};
}

// Primitive bounds.
constexpr kernel bbox3f point_bounds(
    const vector<vec3f>& positions, int point) {
  auto v1 = point;
  return point_bounds(positions[v1]);
}
constexpr kernel bbox3f point_bounds(
    const vector<vec3f>& positions, const vector<float>& radius, int point) {
  auto v1 = point;
  return point_bounds(positions[v1], radius[v1]);
}
constexpr kernel bbox3f line_bounds(
    const vector<vec3f>& positions, const vec2i& line) {
  return line_bounds(positions[line.x], positions[line.y]);
}
constexpr kernel bbox3f line_bounds(const vector<vec3f>& positions,
    const vector<float>& radius, const vec2i& line) {
  return line_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}
constexpr kernel bbox3f triangle_bounds(
    const vector<vec3f>& positions, const vec3i& triangle) {
  return triangle_bounds(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
constexpr kernel bbox3f quad_bounds(
    const vector<vec3f>& positions, const vec4i& quad) {
  return quad_bounds(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}

// Primitive bounds in indexed arrays.
constexpr kernel bbox3f point_bounds(cspan<vec3f> positions, int point) {
  return point_bounds(positions[point]);
}
constexpr kernel bbox3f point_bounds(
    cspan<vec3f> positions, cspan<float> radius, int point) {
  return point_bounds(positions[point], radius[point]);
}
constexpr kernel bbox3f line_bounds(cspan<vec3f> positions, vec2i line) {
  return line_bounds(positions[line.x], positions[line.y]);
}
constexpr kernel bbox3f line_bounds(
    cspan<vec3f> positions, cspan<float> radius, vec2i line) {
  return line_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}
constexpr kernel bbox3f triangle_bounds(
    cspan<vec3f> positions, vec3i triangle) {
  return triangle_bounds(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
constexpr kernel bbox3f quad_bounds(cspan<vec3f> positions, vec4i quad) {
  return quad_bounds(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}
constexpr kernel bbox3f sphere_bounds(
    cspan<vec3f> positions, cspan<float> radius, int point) {
  return sphere_bounds(positions[point], radius[point]);
}
constexpr kernel bbox3f capsule_bounds(
    cspan<vec3f> positions, cspan<float> radius, vec2i line) {
  return capsule_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
constexpr kernel vec3f line_tangent(const vec3f& p1, const vec3f& p2) {
  return normalize(p2 - p1);
}
constexpr kernel vec3f line_tangent(
    const vector<vec3f>& positions, const vec2i& line) {
  return line_tangent(positions[line.x], positions[line.y]);
}
constexpr kernel float line_length(const vec3f& p1, const vec3f& p2) {
  return length(p2 - p1);
}
constexpr kernel float line_length(
    const vector<vec3f>& positions, const vec2i& line) {
  return line_length(positions[line.x], positions[line.y]);
}

// Triangle properties.
constexpr kernel vec3f triangle_normal(
    const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  return normalize(cross(p2 - p1, p3 - p1));
}
constexpr kernel vec3f triangle_normal(
    const vector<vec3f>& positions, const vec3i& triangle) {
  return triangle_normal(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
constexpr kernel float triangle_area(
    const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  return length(cross(p2 - p1, p3 - p1)) / 2;
}
constexpr kernel float triangle_area(
    const vector<vec3f>& positions, const vec3i& triangle) {
  return triangle_area(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}

// Quad properties.
constexpr kernel vec3f quad_normal(
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec3f& p4) {
  return normalize(triangle_normal(p1, p2, p4) + triangle_normal(p3, p4, p2));
}
constexpr kernel vec3f quad_normal(
    const vector<vec3f>& positions, const vec4i& quad) {
  return quad_normal(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}
constexpr kernel float quad_area(
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec3f& p4) {
  return triangle_area(p1, p2, p4) + triangle_area(p3, p4, p2);
}
constexpr kernel float quad_area(
    const vector<vec3f>& positions, const vec4i& quad) {
  return quad_area(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
constexpr kernel T interpolate_line(const T& p1, const T& p2, float u) {
  return p1 * (1 - u) + p2 * u;
}
template <typename T>
constexpr kernel T interpolate_line(
    const vector<T>& vertices, const vec2i& line, float u) {
  return interpolate_line(vertices[line.x], vertices[line.y], u);
}
// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
constexpr kernel T interpolate_line(const T& p1, const T& p2, const vec2f& uv) {
  return interpolate_line(p1, p2, uv.x);
}
template <typename T>
constexpr kernel T interpolate_line(
    const vector<T>& vertices, const vec2i& line, const vec2f& uv) {
  return interpolate_line(vertices, line, uv.x);
}
// Interpolates values over a triangle parameterized by u and v along the
// (p2-p1) and (p3-p1) directions. Same as barycentric interpolation.
template <typename T>
constexpr kernel T interpolate_triangle(
    const T& p1, const T& p2, const T& p3, const vec2f& uv) {
  return p1 * (1 - uv.x - uv.y) + p2 * uv.x + p3 * uv.y;
}
template <typename T>
constexpr kernel T interpolate_triangle(
    const vector<T>& vertices, const vec3i& triangle, const vec2f& uv) {
  return interpolate_triangle(
      vertices[triangle.x], vertices[triangle.y], vertices[triangle.z], uv);
}
// Interpolates values over a quad parameterized by u and v along the
// (p2-p1) and (p3-p2) directions. Same as bilinear interpolation.
template <typename T>
constexpr kernel T interpolate_quad(
    const T& p1, const T& p2, const T& p3, const T& p4, const vec2f& uv) {
  if (sum(uv) <= 1) {
    return interpolate_triangle(p1, p2, p4, uv);
  } else {
    return interpolate_triangle(p3, p4, p2, 1 - uv);
  }
}
template <typename T>
constexpr kernel T interpolate_quad(
    const vector<T>& vertices, const vec4i& quad, const vec2f& uv) {
  return interpolate_quad(vertices[quad.x], vertices[quad.y], vertices[quad.z],
      vertices[quad.w], uv);
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
constexpr kernel T interpolate_bezier(
    const T& p1, const T& p2, const T& p3, const T& p4, float u) {
  return p1 * (1 - u) * (1 - u) * (1 - u) + p2 * 3 * u * (1 - u) * (1 - u) +
         p3 * 3 * u * u * (1 - u) + p4 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
constexpr kernel T interpolate_bezier_derivative(
    const T& p1, const T& p2, const T& p3, const T& p4, float u) {
  return (p2 - p1) * 3 * (1 - u) * (1 - u) + (p3 - p2) * 6 * u * (1 - u) +
         (p4 - p3) * 3 * u * u;
}

// Interpolated line properties.
constexpr kernel vec3f line_point(const vec3f& p1, const vec3f& p2, float u) {
  return p1 * (1 - u) + p2 * u;
}
constexpr kernel vec3f line_point(
    const vector<vec3f>& positions, const vec2i& line, float u) {
  return line_point(positions[line.x], positions[line.y], u);
}
constexpr kernel vec3f line_tangent(const vec3f& t0, const vec3f& t1, float u) {
  return normalize(t0 * (1 - u) + t1 * u);
}
constexpr kernel vec3f line_tangent(
    const vector<vec3f>& tangents, const vec2i& line, float u) {
  return line_tangent(tangents[line.x], tangents[line.y], u);
}

// Interpolated triangle properties.
constexpr kernel vec3f triangle_point(
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec2f& uv) {
  return p1 * (1 - uv.x - uv.y) + p2 * uv.x + p3 * uv.y;
}
constexpr kernel vec3f triangle_point(
    const vector<vec3f>& positions, const vec3i& triangle, const vec2f& uv) {
  return triangle_point(
      positions[triangle.x], positions[triangle.y], positions[triangle.z], uv);
}
constexpr kernel vec3f triangle_normal(
    const vec3f& n0, const vec3f& n1, const vec3f& n2, const vec2f& uv) {
  return normalize(n0 * (1 - uv.x - uv.y) + n1 * uv.x + n2 * uv.y);
}
constexpr kernel vec3f triangle_normal(
    const vector<vec3f>& normals, const vec3i& triangle, const vec2f& uv) {
  return triangle_normal(
      normals[triangle.x], normals[triangle.y], normals[triangle.z], uv);
}

// Interpolated quad properties.
constexpr kernel vec3f quad_point(const vec3f& p1, const vec3f& p2,
    const vec3f& p3, const vec3f& p4, const vec2f& uv) {
  if (sum(uv) <= 1) {
    return triangle_point(p1, p2, p4, uv);
  } else {
    return triangle_point(p3, p4, p2, 1 - uv);
  }
}
constexpr kernel vec3f quad_point(
    const vector<vec3f>& positions, const vec4i& quad, const vec2f& uv) {
  return quad_point(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w], uv);
}
constexpr kernel vec3f quad_normal(const vec3f& n0, const vec3f& n1,
    const vec3f& n2, const vec3f& n3, const vec2f& uv) {
  if (sum(uv) <= 1) {
    return triangle_normal(n0, n1, n3, uv);
  } else {
    return triangle_normal(n2, n3, n1, 1 - uv);
  }
}
constexpr kernel vec3f quad_normal(
    const vector<vec3f>& normals, const vec4i& quad, const vec2f& uv) {
  return quad_normal(
      normals[quad.x], normals[quad.y], normals[quad.z], normals[quad.w], uv);
}

// Interpolated sphere properties.
constexpr kernel vec3f sphere_point(const vec3f p, float r, const vec2f& uv) {
  auto phi = uv.x * 2 * pif, theta = uv.y * pif;
  return p +
         r * vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
constexpr kernel vec3f sphere_normal(const vec3f p, float r, const vec2f& uv) {
  auto phi = uv.x * 2 * pif, theta = uv.y * pif;
  return normalize(
      vec3f{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)});
}

// Triangle tangent and bi-tangent from uv
constexpr kernel pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& p1,
    const vec3f& p2, const vec3f& p3, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  // TODO: do this without indices
  auto p = p2 - p1, q = p3 - p1;
  auto s   = vec2f{uv1[0] - uv0[0], uv2[0] - uv0[0]};
  auto t   = vec2f{uv1[1] - uv0[1], uv2[1] - uv0[1]};
  auto div = cross(s, t);

  if (div != 0) {
    auto tu = vec3f{t[1] * p[0] - t[0] * q[0], t[1] * p[1] - t[0] * q[1],
                  t[1] * p[2] - t[0] * q[2]} /
              div;
    auto tv = vec3f{s[0] * q[0] - s[1] * p[0], s[0] * q[1] - s[1] * p[1],
                  s[0] * q[2] - s[1] * p[2]} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}
constexpr kernel pair<vec3f, vec3f> triangle_tangents_fromuv(
    const vector<vec3f>& positions, const vector<vec2f>& texcoords,
    const vec3i& triangle) {
  return triangle_tangents_fromuv(positions[triangle.x], positions[triangle.y],
      positions[triangle.z], texcoords[triangle.x], texcoords[triangle.y],
      texcoords[triangle.z]);
}

// Quad tangent and bi-tangent from uv.
constexpr kernel pair<vec3f, vec3f> quad_tangents_fromuv(const vec3f& p1,
    const vec3f& p2, const vec3f& p3, const vec3f& p4, const vec2f& uv0,
    const vec2f& uv1, const vec2f& uv2, const vec2f& uv3,
    const vec2f& current_uv) {
  if (sum(current_uv) <= 1) {
    return triangle_tangents_fromuv(p1, p2, p4, uv0, uv1, uv3);
  } else {
    return triangle_tangents_fromuv(p3, p4, p2, uv2, uv3, uv1);
  }
}
constexpr kernel pair<vec3f, vec3f> quad_tangents_fromuv(
    const vector<vec3f>& positions, const vector<vec2f>& texcoords,
    const vec4i& quad, const vec2f& current_uv) {
  return quad_tangents_fromuv(positions[quad.x], positions[quad.y],
      positions[quad.z], positions[quad.w], texcoords[quad.x],
      texcoords[quad.y], texcoords[quad.z], texcoords[quad.w], current_uv);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes the aspect ratio.
constexpr kernel float aspect_ratio(const vec2i& size) {
  return size.x / size.y;
}

// Flip u from [0,1] to [1,0]
constexpr kernel vec2f flip_u(const vec2f& uv) { return {1 - uv.x, uv.y}; }
// Flip v from [0,1] to [1,0]
constexpr kernel vec2f flip_v(const vec2f& uv) { return {uv.x, 1 - uv.y}; }

// Generate a ray from a camera
constexpr kernel ray3f camera_ray(const frame3f& frame, float lens,
    const vec2f& film, const vec2f& image_uv) {
  auto uv  = flip_u(image_uv);
  auto e   = vec3f{0, 0, 0};
  auto q   = vec3f{film * (uv - 0.5f), lens};
  auto d   = normalize(e - q);
  auto ray = yocto::ray3f{
      transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

// Generate a ray from a camera
constexpr kernel ray3f camera_ray(const frame3f& frame, float lens,
    float aspect, float film_, const vec2f& image_uv) {
  auto uv   = flip_u(image_uv);
  auto film = aspect >= 1 ? vec2f{film_, film_ / aspect}
                          : vec2f{film_ * aspect, film_};
  auto e    = vec3f{0, 0, 0};
  auto q    = vec3f{film * (uv - 0.5f), lens};
  auto d    = normalize(e - q);
  auto ray  = yocto::ray3f{
      transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive intersection
struct pintersection3f {
  vec2f uv       = {0, 0};
  float distance = flt_max;  // TODO: flt_max
  bool  hit      = false;

  constexpr kernel pintersection3f() : hit{false} {}
  constexpr kernel pintersection3f(const vec2f& uv_, float distance_) :
      uv{uv_}, distance{distance_}, hit{true} {}
};

// Intersect a ray with a point (approximate)
constexpr kernel pintersection3f intersect_point(
    const ray3f& ray, const vec3f& p, float r) {
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
  return {{0, 0}, t};
}
constexpr kernel pintersection3f intersect_point(const ray3f& ray,
    const vector<vec3f>& positions, const vector<float>& radius, int point) {
  auto v1 = point;
  return intersect_point(ray, positions[v1], radius[v1]);
}

// Intersect a ray with a line
constexpr kernel pintersection3f intersect_line(
    const ray3f& ray, const vec3f& p1, const vec3f& p2, float r1, float r2) {
  // setup intersection params
  auto u = ray.d;
  auto v = p2 - p1;
  auto w = ray.o - p1;

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
  auto pl  = p1 + (p2 - p1) * s;
  auto prl = pr - pl;

  // check with the line radius at the same point
  auto d2 = dot(prl, prl);
  auto r  = r1 * (1 - s) + r2 * s;
  if (d2 > r * r) return {};

  // intersection occurred: set params and exit
  return {{s, sqrt(d2) / r}, t};
}
constexpr kernel pintersection3f intersect_line(const ray3f& ray,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec2i& line) {
  return intersect_line(ray, positions[line.x], positions[line.y],
      radius[line.x], radius[line.y]);
}

// Intersect a ray with a sphere
constexpr kernel pintersection3f intersect_sphere(
    const ray3f& ray, const vec3f& p, float r) {
  // compute parameters
  auto a = dot(ray.d, ray.d);
  auto b = 2 * dot(ray.o - p, ray.d);
  auto c = dot(ray.o - p, ray.o - p) - r * r;

  // check discriminant
  auto dis = b * b - 4 * a * c;
  if (dis < 0) return {};

  // compute ray parameter
  auto t = (-b - sqrt(dis)) / (2 * a);

  // exit if not within boundpif
  if (t < ray.tmin || t > ray.tmax) return {};

  // try other ray parameter
  t = (-b + sqrt(dis)) / (2 * a);

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return {};

  // compute local point for uvs
  auto uv = cartesian_to_sphericaluv(((ray.o + ray.d * t) - p) / r);

  // intersection occurred: set params and exit
  return {uv, t};
}

// Intersect a ray with a triangle
constexpr kernel pintersection3f intersect_triangle(
    const ray3f& ray, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  // compute triangle edges
  auto edge1 = p2 - p1;
  auto edge2 = p3 - p1;

  // compute determinant to solve a linear system
  auto pvec = cross(ray.d, edge2);
  auto det  = dot(edge1, pvec);

  // check determinant and exit if triangle and ray are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return {};
  auto inv_det = 1.0f / det;

  // compute and check first barycentric coordinate
  auto tvec = ray.o - p1;
  auto u    = dot(tvec, pvec) * inv_det;
  if (u < 0 || u > 1) return {};

  // compute and check second barycentric coordinate
  auto qvec = cross(tvec, edge1);
  auto v    = dot(ray.d, qvec) * inv_det;
  if (v < 0 || u + v > 1) return {};

  // compute and check ray parameter
  auto t = dot(edge2, qvec) * inv_det;
  if (t < ray.tmin || t > ray.tmax) return {};

  // intersection occurred: set params and exit
  return {{u, v}, t};
}
constexpr kernel pintersection3f intersect_triangle(
    const ray3f& ray, const vector<vec3f>& positions, const vec3i& triangle) {
  return intersect_triangle(
      ray, positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}

// Intersect a ray with a quad.
constexpr kernel pintersection3f intersect_quad(const ray3f& ray,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec3f& p4) {
  if (p3 == p4) return intersect_triangle(ray, p1, p2, p4);
  auto isec1 = intersect_triangle(ray, p1, p2, p4);
  auto isec2 = intersect_triangle(ray, p3, p4, p2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
constexpr kernel pintersection3f intersect_quad(
    const ray3f& ray, const vector<vec3f>& positions, const vec4i& quad) {
  return intersect_quad(ray, positions[quad.x], positions[quad.y],
      positions[quad.z], positions[quad.w]);
}

// Intersect a ray with a axis-aligned bounding box
constexpr kernel bool intersect_bbox(const ray3f& ray, const bbox3f& bbox) {
  auto ray_dinv = 1 / ray.d;
  auto it_min   = (bbox.min - ray.o) * ray_dinv;
  auto it_max   = (bbox.max - ray.o) * ray_dinv;
  auto tmin     = min(it_min, it_max);
  auto tmax     = max(it_min, it_max);
  auto t0       = max(max(tmin), ray.tmin);
  auto t1       = min(min(tmax), ray.tmax);
  t1 *= 1.00000024f;
  return t0 <= t1;
}

// Intersect a ray with a axis-aligned bounding box
constexpr kernel bool intersect_bbox(
    const ray3f& ray, const vec3f& ray_dinv, const bbox3f& bbox) {
  auto it_min = (bbox.min - ray.o) * ray_dinv;
  auto it_max = (bbox.max - ray.o) * ray_dinv;
  auto tmin   = min(it_min, it_max);
  auto tmax   = max(it_min, it_max);
  auto t0     = max(max(tmin), ray.tmin);
  auto t1     = min(min(tmax), ray.tmax);
  t1 *= 1.00000024f;
  return t0 <= t1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
constexpr kernel pintersection3f overlap_point(
    const vec3f& pos, float dist_max, const vec3f& p, float r) {
  auto d2 = dot(pos - p, pos - p);
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  return {{0, 0}, sqrt(d2)};
}
constexpr kernel pintersection3f overlap_point(const vec3f& pos, float dist_max,
    const vector<vec3f>& positions, const vector<float>& radius, int point) {
  auto v1 = point;
  return overlap_point(pos, dist_max, positions[v1], radius[v1]);
}

// Compute the closest line uv to a give position pos.
constexpr kernel float closestuv_line(
    const vec3f& pos, const vec3f& p1, const vec3f& p2) {
  auto ab = p2 - p1;
  auto d  = dot(ab, ab);
  // Project c onto ab, computing parameterized position d(t) = a + t*(b –
  // a)
  auto u = dot(pos - p1, ab) / d;
  u      = clamp(u, 0.0f, 1.0f);
  return u;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
constexpr kernel pintersection3f overlap_line(const vec3f& pos, float dist_max,
    const vec3f& p1, const vec3f& p2, float r1, float r2) {
  auto u = closestuv_line(pos, p1, p2);
  // Compute projected position from the clamped t d = a + t * ab;
  auto p  = p1 + (p2 - p1) * u;
  auto r  = r1 + (r2 - r1) * u;
  auto d2 = dot(pos - p, pos - p);
  // check distance
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  // done
  return {{u, 0}, sqrt(d2)};
}
constexpr kernel pintersection3f overlap_line(const vec3f& pos, float dist_max,
    const vector<vec3f>& positions, const vector<float> radius,
    const vec2i& line) {
  return overlap_line(pos, dist_max, positions[line.x], positions[line.y],
      radius[line.x], radius[line.y]);
}

// Compute the closest triangle uv to a give position pos.
constexpr kernel vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  // this is a complicated test -> I probably "--"+prefix to use a sequence of
  // test (triangle body, and 3 edges)
  auto ab = p2 - p1;
  auto ac = p3 - p1;
  auto ap = pos - p1;

  auto d1 = dot(ab, ap);
  auto d2 = dot(ac, ap);

  // corner and edge cases
  if (d1 <= 0 && d2 <= 0) return {0, 0};

  auto bp = pos - p2;
  auto d3 = dot(ab, bp);
  auto d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return {1, 0};

  auto vc = d1 * d4 - d3 * d2;
  if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return {d1 / (d1 - d3), 0};

  auto cp = pos - p3;
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
constexpr kernel pintersection3f overlap_triangle(const vec3f& pos,
    float dist_max, const vec3f& p1, const vec3f& p2, const vec3f& p3, float r1,
    float r2, float r3) {
  auto uv = closestuv_triangle(pos, p1, p2, p3);
  auto p  = interpolate_triangle(p1, p2, p3, uv);
  auto r  = interpolate_triangle(r1, r2, r3, uv);
  auto dd = dot(p - pos, p - pos);
  if (dd > (dist_max + r) * (dist_max + r)) return {};
  return {uv, sqrt(dd)};
}
constexpr kernel pintersection3f overlap_triangle(const vec3f& pos,
    float dist_max, const vector<vec3f>& positions, const vector<float> radius,
    const vec3i& triangle) {
  return overlap_triangle(pos, dist_max, positions[triangle.x],
      positions[triangle.y], positions[triangle.z], radius[triangle.x],
      radius[triangle.y], radius[triangle.z]);
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
constexpr kernel pintersection3f overlap_quad(const vec3f& pos, float dist_max,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, const vec3f& p4,
    float r1, float r2, float r3, float r4) {
  if (p3 == p4) return overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec1 = overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec2 = overlap_triangle(pos, dist_max, p3, p4, p2, r3, r4, r2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
constexpr kernel pintersection3f overlap_quad(const vec3f& pos, float dist_max,
    const vector<vec3f>& positions, const vector<float> radius,
    const vec4i& quad) {
  return overlap_quad(pos, dist_max, positions[quad.x], positions[quad.y],
      positions[quad.z], positions[quad.w], radius[quad.x], radius[quad.y],
      radius[quad.z], radius[quad.w]);
}

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
constexpr kernel bool overlap_bbox(
    const vec3f& pos, float dist_max, const bbox3f& bbox) {
  // computing distance
  auto dd = 0.0f;

  // For each axis count any excess distance outside box extents
  for (auto a : range(3)) {
    if (pos[a] < bbox.min[a])
      dd += (bbox.min[a] - pos[a]) * (bbox.min[a] - pos[a]);
    if (pos[a] > bbox.max[a])
      dd += (pos[a] - bbox.max[a]) * (pos[a] - bbox.max[a]);
  }

  // check distance
  return dd < dist_max * dist_max;
}

// Check if two bboxes overlap.
constexpr kernel bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
  for (auto a : range(3)) {
    if (bbox1.max[a] < bbox2.min[a] || bbox1.min[a] > bbox2.max[a])
      return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
