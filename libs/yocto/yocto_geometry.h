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
template <typename T, size_t N>
struct bbox;

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 1> {
  vec<T, 1> min = {num_max<T>};
  vec<T, 1> max = {num_min<T>};

  constexpr kernel bbox() : min{num_max<T>}, max{num_min<T>} {}
  constexpr kernel bbox(const vec<T, 1>& min_, const vec<T, 1>& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec<T, 1>& operator[](size_t i) {
    return i == 0 ? min : max;
  }
  constexpr kernel const vec<T, 1>& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 2> {
  vec<T, 2> min = {num_max<T>, num_max<T>};
  vec<T, 2> max = {num_min<T>, num_min<T>};

  constexpr kernel bbox() :
      min{num_max<T>, num_max<T>}, max{num_min<T>, num_min<T>} {}
  constexpr kernel bbox(const vec<T, 2>& min_, const vec<T, 2>& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec<T, 2>& operator[](size_t i) {
    return i == 0 ? min : max;
  }
  constexpr kernel const vec<T, 2>& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 3> {
  vec<T, 3> min = {num_max<T>, num_max<T>, num_max<T>};
  vec<T, 3> max = {num_min<T>, num_min<T>, num_min<T>};

  constexpr kernel bbox() :
      min{num_max<T>, num_max<T>, num_max<T>},
      max{num_min<T>, num_min<T>, num_min<T>} {}
  constexpr kernel bbox(const vec<T, 3>& min_, const vec<T, 3>& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec<T, 3>& operator[](size_t i) {
    return i == 0 ? min : max;
  }
  constexpr kernel const vec<T, 3>& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 4> {
  vec<T, 4> min = {num_max<T>, num_max<T>, num_max<T>, num_max<T>};
  vec<T, 4> max = {num_min<T>, num_min<T>, num_min<T>, num_min<T>};

  constexpr kernel bbox() :
      min{num_max<T>, num_max<T>, num_max<T>, num_max<T>},
      max{num_min<T>, num_min<T>, num_min<T>, num_min<T>} {}
  constexpr kernel bbox(const vec<T, 4>& min_, const vec<T, 4>& max_) :
      min{min_}, max{max_} {}

  constexpr kernel vec<T, 4>& operator[](size_t i) {
    return i == 0 ? min : max;
  }
  constexpr kernel const vec<T, 4>& operator[](size_t i) const {
    return i == 0 ? min : max;
  }
};

// Bbox aliases
using bbox1f = bbox<float, 1>;
using bbox2f = bbox<float, 2>;
using bbox3f = bbox<float, 3>;
using bbox4f = bbox<float, 4>;

// Empty bbox constant.
constexpr auto invalidb1f = bbox1f{};
constexpr auto invalidb2f = bbox2f{};
constexpr auto invalidb3f = bbox3f{};
constexpr auto invalidb4f = bbox4f{};

// Bounding box properties
template <typename T, size_t N>
constexpr kernel vec<T, N> bbox_center(const bbox<T, N>& a) {
  return (a.min + a.max) / 2;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> bbox_diagonal(const bbox<T, N>& a) {
  return a.max - a.min;
}
template <typename T>
constexpr kernel T bbox_area(const bbox<T, 3>& a) {
  auto d = a.max - a.min;
  return 2 * d.x * d.y + 2 * d.x * d.z + 2 * d.y * d.z;
}

// Bounding box comparisons.
template <typename T, size_t N>
constexpr kernel bool operator==(const bbox<T, N>& a, const bbox<T, N>& b) {
  return a.min == b.min && a.max == b.max;
}
template <typename T, size_t N>
constexpr kernel bool operator!=(const bbox<T, N>& a, const bbox<T, N>& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T, size_t N>
constexpr kernel bbox<T, N> merge(const bbox<T, N>& a, const vec<T, N>& b) {
  return {min(a.min, b), max(a.max, b)};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> merge(const bbox<T, N>& a, const bbox<T, N>& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
template <typename T, size_t N>
constexpr kernel void expand(bbox<T, N>& a, const vec<T, N>& b) {
  a = merge(a, b);
}
template <typename T, size_t N>
constexpr kernel void expand(bbox<T, N>& a, const bbox<T, N>& b) {
  a = merge(a, b);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Ray epsilon
template <typename T>
constexpr auto ray_eps = (T)1e-4;

// Rays with origin, direction and min/max t value.
template <typename T, size_t N>
struct ray;

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray<T, 2> {
  vec<T, 3> o    = {0, 0};
  vec<T, 3> d    = {0, 1};
  T         tmin = ray_eps<T>;
  T         tmax = num_max<T>;

  constexpr kernel ray() :
      o{0, 0}, d{0, 1}, tmin{ray_eps<T>}, tmax{num_max<T>} {}
  constexpr kernel ray(const vec<T, 2>& o_, const vec<T, 2>& d_,
      T tmin_ = ray_eps<T>, T tmax_ = num_max<T>) :
      o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray<T, 3> {
  vec<T, 3> o    = {0, 0, 0};
  vec<T, 3> d    = {0, 0, 1};
  T         tmin = ray_eps<T>;
  T         tmax = num_max<T>;

  constexpr kernel ray() : o{0}, d{0}, tmin{ray_eps<T>}, tmax{num_max<T>} {}
  constexpr kernel ray(const vec<T, 3>& o_, const vec<T, 3>& d_,
      T tmin_ = ray_eps<T>, T tmax_ = num_max<T>) :
      o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
};

// Ray aliases
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// Computes a point on a ray
template <typename T, size_t N>
constexpr kernel vec<T, N> ray_point(const ray<T, N>& ray, T t) {
  return ray.o + ray.d * t;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> eval_ray(const ray<T, N>& ray, T t) {
  return ray.o + ray.d * t;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms rays.
template <typename T, size_t N>
constexpr kernel ray<T, N> transform_ray(
    const mat<T, N + 1, N + 1>& a, const ray<T, N>& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}

template <typename T, size_t N>
constexpr kernel ray<T, N> transform_ray(
    const frame<T, N>& a, const ray<T, N>& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}

// Transforms bboxes.
template <typename T, size_t N>
constexpr kernel bbox<T, N> transform_bbox(
    const mat<T, N + 1, N + 1>& a, const bbox<T, N>& b) {
  if constexpr (N == 2) {
    auto xformed = bbox<T, N>();
    for (auto ij : range(vec<int, 2>(2, 2))) {
      auto corner = vec<T, 2>{b[ij.x][0], b[ij.y][1]};
      xformed     = merge(xformed, transform_point(a, corner));
    }
    return xformed;
  } else if constexpr (N == 3) {
    auto xformed = bbox<T, N>();
    for (auto ijk : range(vec<int, 3>(2, 2, 2))) {
      auto corner = vec<T, 3>{b[ijk.x][0], b[ijk.y][1], b[ijk.z][2]};
      xformed     = merge(xformed, transform_point(a, corner));
    }
    return xformed;
  }
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> transform_bbox(
    const frame<T, N>& a, const bbox<T, N>& b) {
  if constexpr (N == 2) {
    auto xformed = bbox<T, N>();
    for (auto ij : range(vec<int, 2>(2, 2))) {
      auto corner = vec<T, 2>{b[ij.x][0], b[ij.y][1]};
      xformed     = merge(xformed, transform_point(a, corner));
    }
    return xformed;
  } else if constexpr (N == 3) {
    auto xformed = bbox<T, N>();
    for (auto ijk : range(vec<int, 3>(2, 2, 2))) {
      auto corner = vec<T, 3>{b[ijk.x][0], b[ijk.y][1], b[ijk.z][2]};
      xformed     = merge(xformed, transform_point(a, corner));
    }
    return xformed;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUAD CONVENTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a quad is a triangle
template <typename I>
constexpr kernel bool is_triangle(const vec<I, 4>& quad) {
  return quad.z == quad.w;
}

// Get the triangle from the quad
template <typename I>
constexpr kernel vec<I, 3> as_triangle(const vec<I, 4>& quad) {
  return xyz(quad);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRIMITIVE BOUNDS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive bounds.
template <typename T, size_t N>
constexpr kernel bbox<T, N> point_bounds(const vec<T, N>& p) {
  return {p, p};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> point_bounds(const vec<T, N>& p, T r) {
  return {min(p - r, p + r), max(p - r, p + r)};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> line_bounds(
    const vec<T, N>& p1, const vec<T, N>& p2) {
  return {min(p1, p2), max(p1, p2)};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> line_bounds(
    const vec<T, N>& p1, const vec<T, N>& p2, T r1, T r2) {
  return {min(p1 - r1, p2 - r2), max(p1 + r1, p2 + r2)};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> triangle_bounds(
    const vec<T, N>& p1, const vec<T, N>& p2, const vec<T, N>& p3) {
  return {min(p1, min(p2, p3)), max(p1, max(p2, p3))};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> quad_bounds(const vec<T, N>& p1,
    const vec<T, N>& p2, const vec<T, N>& p3, const vec<T, N>& p4) {
  return {min(p1, min(p2, min(p3, p4))), max(p1, max(p2, max(p3, p4)))};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> sphere_bounds(const vec<T, N>& p, T r) {
  return {p - r, p + r};
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> capsule_bounds(
    const vec<T, N>& p1, const vec<T, N>& p2, T r1, T r2) {
  return {min(p1 - r1, p2 - r2), max(p1 + r1, p2 + r2)};
}

// Primitive bounds.
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> point_bounds(
    const vector<vec<T, N>>& positions, I point) {
  auto v1 = point;
  return point_bounds(positions[v1]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> point_bounds(
    const vector<vec<T, N>>& positions, const vector<T>& radius, I point) {
  auto v1 = point;
  return point_bounds(positions[v1], radius[v1]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(
    const vector<vec<T, N>>& positions, const vec<I, 2>& line) {
  return line_bounds(positions[line.x], positions[line.y]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(const vector<vec<T, N>>& positions,
    const vector<T>& radius, const vec<I, 2>& line) {
  return line_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> triangle_bounds(
    const vector<vec<T, N>>& positions, const vec<I, 3>& triangle) {
  return triangle_bounds(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> quad_bounds(
    const vector<vec<T, N>>& positions, const vec<I, 4>& quad) {
  return quad_bounds(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}

// Primitive bounds in indexed arrays.
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> point_bounds(cspan<vec<T, N>> positions, I point) {
  return point_bounds(positions[point]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> point_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, I point) {
  return point_bounds(positions[point], radius[point]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(
    cspan<vec<T, N>> positions, vec<I, 2> line) {
  return line_bounds(positions[line.x], positions[line.y]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, vec<I, 2> line) {
  return line_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> triangle_bounds(
    cspan<vec<T, N>> positions, vec<I, 3> triangle) {
  return triangle_bounds(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> quad_bounds(
    cspan<vec<T, N>> positions, vec<I, 4> quad) {
  return quad_bounds(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.q]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> sphere_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, I point) {
  return sphere_bounds(positions[point], radius[point]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> capsule_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, vec<I, 2> line) {
  return capsule_bounds(
      positions[line.x], positions[line.y], radius[line.x], radius[line.y]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
template <typename T>
constexpr kernel vec<T, 3> line_tangent(
    const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return normalize(p2 - p1);
}
template <typename T, typename I>
constexpr kernel vec<T, 3> line_tangent(
    const vector<vec<T, 3>>& positions, const vec<I, 2>& line) {
  return line_tangent(positions[line.x], positions[line.y]);
}
template <typename T>
constexpr kernel T line_length(const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return length(p2 - p1);
}
template <typename T, typename I>
constexpr kernel T line_length(
    const vector<vec<T, 3>>& positions, const vec<I, 2>& line) {
  return line_length(positions[line.x], positions[line.y]);
}

// Triangle properties.
template <typename T>
constexpr kernel vec<T, 3> triangle_normal(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
  return normalize(cross(p2 - p1, p3 - p1));
}
template <typename T, typename I>
constexpr kernel vec<T, 3> triangle_normal(
    const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  return triangle_normal(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}
template <typename T>
constexpr kernel T triangle_area(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
  return length(cross(p2 - p1, p3 - p1)) / 2;
}
template <typename T, typename I>
constexpr kernel T triangle_area(
    const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  return triangle_area(
      positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}

// Quad properties.
template <typename T>
constexpr kernel vec<T, 3> quad_normal(const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 3>& p3, const vec<T, 3>& p4) {
  return normalize(triangle_normal(p1, p2, p4) + triangle_normal(p3, p4, p2));
}
template <typename T, typename I>
constexpr kernel vec<T, 3> quad_normal(
    const vector<vec<T, 3>>& positions, const vec<I, 4>& quad) {
  return quad_normal(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}
template <typename T>
constexpr kernel T quad_area(const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 3>& p3, const vec<T, 3>& p4) {
  return triangle_area(p1, p2, p4) + triangle_area(p3, p4, p2);
}
template <typename T, typename I>
constexpr kernel T quad_area(
    const vector<vec<T, 3>>& positions, const vec<I, 4>& quad) {
  return quad_area(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w]);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
constexpr kernel T interpolate_line(const T& p1, const T& p2, T1 u) {
  return p1 * (1 - u) + p2 * u;
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_line(
    const vector<T>& vertices, const vec<I, 2>& line, T1 u) {
  return interpolate_line(vertices[line.x], vertices[line.y], u);
}
// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
constexpr kernel T interpolate_line(
    const T& p1, const T& p2, const vec<T1, 2>& uv) {
  return interpolate_line(p1, p2, uv.x);
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_line(
    const vector<T>& vertices, const vec<I, 2>& line, const vec<T1, 2>& uv) {
  return interpolate_line(vertices, line, uv.x);
}
// Interpolates values over a triangle parameterized by u and v along the
// (p2-p1) and (p3-p1) directions. Same as barycentric interpolation.
template <typename T, typename T1>
constexpr kernel T interpolate_triangle(
    const T& p1, const T& p2, const T& p3, const vec<T1, 2>& uv) {
  return p1 * (1 - uv.x - uv.y) + p2 * uv.x + p3 * uv.y;
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_triangle(
    const vector<T>& vertices, const vec<I, 3>& triangle, T1 u) {
  return interpolate_triangle(
      vertices[triangle.x], vertices[triangle.y], vertices[triangle.z], u);
}
// Interpolates values over a quad parameterized by u and v along the
// (p2-p1) and (p3-p2) directions. Same as bilinear interpolation.
template <typename T, typename T1>
constexpr kernel T interpolate_quad(
    const T& p1, const T& p2, const T& p3, const T& p4, const vec<T1, 2>& uv) {
  if (sum(uv) <= 1) {
    return interpolate_triangle(p1, p2, p4, uv);
  } else {
    return interpolate_triangle(p3, p4, p2, 1 - uv);
  }
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_quad(
    const vector<T>& vertices, const vec<I, 4>& quad, T1 u) {
  return interpolate_quad(vertices[quad.x], vertices[quad.y], vertices[quad.z],
      vertices[quad.w], u);
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
constexpr kernel T interpolate_bezier(
    const T& p1, const T& p2, const T& p3, const T& p4, T1 u) {
  return p1 * (1 - u) * (1 - u) * (1 - u) + p2 * 3 * u * (1 - u) * (1 - u) +
         p3 * 3 * u * u * (1 - u) + p4 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
constexpr kernel T interpolate_bezier_derivative(
    const T& p1, const T& p2, const T& p3, const T& p4, T1 u) {
  return (p2 - p1) * 3 * (1 - u) * (1 - u) + (p3 - p2) * 6 * u * (1 - u) +
         (p4 - p3) * 3 * u * u;
}

// Interpolated line properties.
template <typename T>
constexpr kernel vec<T, 3> line_point(
    const vec<T, 3>& p1, const vec<T, 3>& p2, T u) {
  return p1 * (1 - u) + p2 * u;
}
template <typename T, typename I>
constexpr kernel vec<T, 3> line_point(
    const vector<vec<T, 3>>& positions, const vec<I, 2>& line, T u) {
  return line_point(positions[line.x], positions[line.y], u);
}
template <typename T>
constexpr kernel vec<T, 3> line_tangent(
    const vec<T, 3>& t0, const vec<T, 3>& t1, T u) {
  return normalize(t0 * (1 - u) + t1 * u);
}
template <typename T, typename I>
constexpr kernel vec<T, 3> line_tangent(
    const vector<vec<T, 3>>& tangents, const vec<I, 2>& line, T u) {
  return line_tangent(tangents[line.x], tangents[line.y], u);
}

// Interpolated triangle properties.
template <typename T>
constexpr kernel vec<T, 3> triangle_point(const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3, const vec<T, 2>& uv) {
  return p1 * (1 - uv.x - uv.y) + p2 * uv.x + p3 * uv.y;
}
template <typename T, typename I>
constexpr kernel vec<T, 3> triangle_point(const vector<vec<T, 3>>& positions,
    const vec<I, 3>& triangle, const vec<T, 2>& uv) {
  return triangle_point(
      positions[triangle.x], positions[triangle.y], positions[triangle.z], uv);
}
template <typename T>
constexpr kernel vec<T, 3> triangle_normal(const vec<T, 3>& n0,
    const vec<T, 3>& n1, const vec<T, 3>& n2, const vec<T, 2>& uv) {
  return normalize(n0 * (1 - uv.x - uv.y) + n1 * uv.x + n2 * uv.y);
}
template <typename T, typename I>
constexpr kernel vec<T, 3> triangle_normal(const vector<vec<T, 3>>& normals,
    const vec<I, 3>& triangle, const vec<T, 2>& uv) {
  return triangle_normal(
      normals[triangle.x], normals[triangle.y], normals[triangle.z], uv);
}

// Interpolated quad properties.
template <typename T>
constexpr kernel vec<T, 3> quad_point(const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 3>& p3, const vec<T, 3>& p4, const vec<T, 2>& uv) {
  if (sum(uv) <= 1) {
    return triangle_point(p1, p2, p4, uv);
  } else {
    return triangle_point(p3, p4, p2, 1 - uv);
  }
}
template <typename T, typename I>
constexpr kernel vec<T, 3> quad_point(const vector<vec<T, 3>>& positions,
    const vec<I, 4>& quad, const vec<T, 2>& uv) {
  return quad_point(positions[quad.x], positions[quad.y], positions[quad.z],
      positions[quad.w], uv);
}
template <typename T>
constexpr kernel vec<T, 3> quad_normal(const vec<T, 3>& n0, const vec<T, 3>& n1,
    const vec<T, 3>& n2, const vec<T, 3>& n3, const vec<T, 2>& uv) {
  if (sum(uv) <= 1) {
    return triangle_normal(n0, n1, n3, uv);
  } else {
    return triangle_normal(n2, n3, n1, 1 - uv);
  }
}
template <typename T, typename I>
constexpr kernel vec<T, 3> quad_normal(const vector<vec<T, 3>>& normals,
    const vec<I, 4>& quad, const vec<T, 2>& uv) {
  return quad_normal(
      normals[quad.x], normals[quad.y], normals[quad.z], normals[quad.w], uv);
}

// Interpolated sphere properties.
template <typename T>
constexpr kernel vec<T, 3> sphere_point(
    const vec<T, 3> p, T r, const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return p + r * vec<T, 3>{
                     cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
template <typename T>
constexpr kernel vec<T, 3> sphere_normal(
    const vec<T, 3> p, T r, const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return normalize(
      vec<T, 3>{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)});
}

// Triangle tangent and bi-tangent from uv
template <typename T>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  // TODO: do this without indices
  auto p = p2 - p1, q = p3 - p1;
  auto s   = vec<T, 2>{uv1[0] - uv0[0], uv2[0] - uv0[0]};
  auto t   = vec<T, 2>{uv1[1] - uv0[1], uv2[1] - uv0[1]};
  auto div = cross(s, t);

  if (div != 0) {
    auto tu = vec<T, 3>{t[1] * p[0] - t[0] * q[0], t[1] * p[1] - t[0] * q[1],
                  t[1] * p[2] - t[0] * q[2]} /
              div;
    auto tv = vec<T, 3>{s[0] * q[0] - s[1] * p[0], s[0] * q[1] - s[1] * p[1],
                  s[0] * q[2] - s[1] * p[2]} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}
template <typename T, typename I>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 2>>& texcoords,
    const vec<I, 3>& triangle) {
  return triangle_tangents_fromuv(positions[triangle.x], positions[triangle.y],
      positions[triangle.z], texcoords[triangle.x], texcoords[triangle.y],
      texcoords[triangle.z]);
}

// Quad tangent and bi-tangent from uv.
template <typename T>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4, const vec<T, 2>& uv0, const vec<T, 2>& uv1,
    const vec<T, 2>& uv2, const vec<T, 2>& uv3, const vec<T, 2>& current_uv) {
  if (sum(current_uv) <= 1) {
    return triangle_tangents_fromuv(p1, p2, p4, uv0, uv1, uv3);
  } else {
    return triangle_tangents_fromuv(p3, p4, p2, uv2, uv3, uv1);
  }
}
template <typename T, typename I>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 2>>& texcoords,
    const vec<I, 4>& quad, const vec<T, 2>& current_uv) {
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
template <typename T>
constexpr kernel T aspect_ratio(const vec<T, 2>& size) {
  return size.x / size.y;
}

// Flip u from [0,1] to [1,0]
template <typename T>
constexpr kernel vec<T, 2> flip_u(const vec<T, 2>& uv) {
  return {1 - uv.x, uv.y};
}
// Flip v from [0,1] to [1,0]
template <typename T>
constexpr kernel vec<T, 2> flip_v(const vec<T, 2>& uv) {
  return {uv.x, 1 - uv.y};
}

// Generate a ray from a camera
template <typename T>
constexpr kernel ray<T, 3> camera_ray(const frame<T, 3>& frame, T lens,
    const vec<T, 2>& film, const vec<T, 2>& image_uv) {
  auto uv  = flip_u(image_uv);
  auto e   = vec<T, 3>{0, 0, 0};
  auto q   = vec<T, 3>{film * (uv - (T)0.5), lens};
  auto d   = normalize(e - q);
  auto ray = yocto::ray<T, 3>{
      transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

// Generate a ray from a camera
template <typename T>
constexpr kernel ray<T, 3> camera_ray(const frame<T, 3>& frame, T lens,
    T aspect, T film_, const vec<T, 2>& image_uv) {
  auto uv   = flip_u(image_uv);
  auto film = aspect >= 1 ? vec<T, 2>{film_, film_ / aspect}
                          : vec<T, 2>{film_ * aspect, film_};
  auto e    = vec<T, 3>{0, 0, 0};
  auto q    = vec<T, 3>{film * (uv - (T)0.5), lens};
  auto d    = normalize(e - q);
  auto ray  = yocto::ray<T, 3>{
      transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Primitive intersection
template <typename T, size_t N>
struct pintersection {
  vec<T, 2> uv       = {0, 0};
  T         distance = num_max<T>;  // TODO: num_max<T>
  bool      hit      = false;

  constexpr kernel pintersection() : hit{false} {}
  constexpr kernel pintersection(const vec<T, 2>& uv_, T distance_) :
      uv{uv_}, distance{distance_}, hit{true} {}
};

// Typedefs
using pintersection3f = pintersection<float, 3>;

// Intersect a ray with a point (approximate)
template <typename T>
constexpr kernel pintersection<T, 3> intersect_point(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r) {
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
template <typename T, typename I>
constexpr kernel pintersection<T, 3> intersect_point(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vector<T>& radius, I point) {
  auto v1 = point;
  return intersect_point(ray, positions[v1], radius[v1]);
}

// Intersect a ray with a line
template <typename T>
constexpr kernel pintersection<T, 3> intersect_line(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, T r1, T r2) {
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
  s = clamp(s, (T)0, (T)1);

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
template <typename T, typename I>
constexpr kernel pintersection<T, 3> intersect_line(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vector<T>& radius,
    const vec<I, 2>& line) {
  return intersect_line(ray, positions[line.x], positions[line.y],
      radius[line.x], radius[line.y]);
}

// Intersect a ray with a sphere
template <typename T>
constexpr kernel pintersection<T, 3> intersect_sphere(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r) {
  // compute parameters
  auto a = dot(ray.d, ray.d);
  auto b = 2 * dot(ray.o - p, ray.d);
  auto c = dot(ray.o - p, ray.o - p) - r * r;

  // check discriminant
  auto dis = b * b - 4 * a * c;
  if (dis < 0) return {};

  // compute ray parameter
  auto t = (-b - sqrt(dis)) / (2 * a);

  // exit if not within bound(T)pi
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
template <typename T>
constexpr kernel pintersection<T, 3> intersect_triangle(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
  // compute triangle edges
  auto edge1 = p2 - p1;
  auto edge2 = p3 - p1;

  // compute determinant to solve a linear system
  auto pvec = cross(ray.d, edge2);
  auto det  = dot(edge1, pvec);

  // check determinant and exit if triangle and ray are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return {};
  auto inv_det = (T)1.0 / det;

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
template <typename T, typename I>
constexpr kernel pintersection<T, 3> intersect_triangle(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  return intersect_triangle(
      ray, positions[triangle.x], positions[triangle.y], positions[triangle.z]);
}

// Intersect a ray with a quad.
template <typename T>
constexpr kernel pintersection<T, 3> intersect_quad(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4) {
  if (p3 == p4) return intersect_triangle(ray, p1, p2, p4);
  auto isec1 = intersect_triangle(ray, p1, p2, p4);
  auto isec2 = intersect_triangle(ray, p3, p4, p2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
template <typename T, typename I>
constexpr kernel pintersection<T, 3> intersect_quad(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vec<I, 4>& quad) {
  return intersect_quad(ray, positions[quad.x], positions[quad.y],
      positions[quad.z], positions[quad.w]);
}

// Intersect a ray with a axis-aligned bounding box
template <typename T>
constexpr kernel bool intersect_bbox(
    const ray<T, 3>& ray, const bbox<T, 3>& bbox) {
  auto ray_dinv = 1 / ray.d;
  auto it_min   = (bbox.min - ray.o) * ray_dinv;
  auto it_max   = (bbox.max - ray.o) * ray_dinv;
  auto tmin     = min(it_min, it_max);
  auto tmax     = max(it_min, it_max);
  auto t0       = max(max(tmin), ray.tmin);
  auto t1       = min(min(tmax), ray.tmax);
  if constexpr (std::is_same_v<T, float>) t1 *= 1.00000024f;
  if constexpr (std::is_same_v<T, double>) t1 *= 1.0000000000000004;
  return t0 <= t1;
}

// Intersect a ray with a axis-aligned bounding box
template <typename T>
constexpr kernel bool intersect_bbox(
    const ray<T, 3>& ray, const vec<T, 3>& ray_dinv, const bbox<T, 3>& bbox) {
  auto it_min = (bbox.min - ray.o) * ray_dinv;
  auto it_max = (bbox.max - ray.o) * ray_dinv;
  auto tmin   = min(it_min, it_max);
  auto tmax   = max(it_min, it_max);
  auto t0     = max(max(tmin), ray.tmin);
  auto t1     = min(min(tmax), ray.tmax);
  if constexpr (std::is_same_v<T, float>) t1 *= 1.00000024f;
  if constexpr (std::is_same_v<T, double>) t1 *= 1.0000000000000004;
  return t0 <= t1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel pintersection<T, 3> overlap_point(
    const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p, T r) {
  auto d2 = dot(pos - p, pos - p);
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  return {{0, 0}, sqrt(d2)};
}
template <typename T, typename I>
constexpr kernel pintersection<T, 3> overlap_point(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    I point) {
  auto v1 = point;
  return overlap_point(pos, dist_max, positions[v1], radius[v1]);
}

// Compute the closest line uv to a give position pos.
template <typename T>
constexpr kernel T closestuv_line(
    const vec<T, 3>& pos, const vec<T, 3>& p1, const vec<T, 3>& p2) {
  auto ab = p2 - p1;
  auto d  = dot(ab, ab);
  // Project c onto ab, computing parameterized position d(t) = a + t*(b –
  // a)
  auto u = dot(pos - p1, ab) / d;
  u      = clamp(u, (T)0, (T)1);
  return u;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel pintersection<T, 3> overlap_line(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, T r1, T r2) {
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
template <typename T, typename I>
constexpr kernel pintersection<T, 3> overlap_line(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 2>& line) {
  return overlap_line(pos, dist_max, positions[line.x], positions[line.y],
      radius[line.x], radius[line.y]);
}

// Compute the closest triangle uv to a give position pos.
template <typename T>
constexpr kernel vec<T, 2> closestuv_triangle(const vec<T, 3>& pos,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
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
template <typename T, typename I>
constexpr kernel pintersection<T, 3> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  return overlap_triangle(pos, dist_max, positions[triangle.x],
      positions[triangle.y], positions[triangle.z]);
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
template <typename T>
constexpr kernel pintersection<T, 3> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    T r1, T r2, T r3) {
  auto uv = closestuv_triangle(pos, p1, p2, p3);
  auto p  = interpolate_triangle(p1, p2, p3, uv);
  auto r  = interpolate_triangle(r1, r2, r3, uv);
  auto dd = dot(p - pos, p - pos);
  if (dd > (dist_max + r) * (dist_max + r)) return {};
  return {uv, sqrt(dd)};
}
template <typename T, typename I>
constexpr kernel pintersection<T, 3> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 3>& triangle) {
  return overlap_triangle(pos, dist_max, positions[triangle.x],
      positions[triangle.y], positions[triangle.z], radius[triangle.x],
      radius[triangle.y], radius[triangle.z]);
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel pintersection<T, 3> overlap_quad(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4, T r1, T r2, T r3, T r4) {
  if (p3 == p4) return overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec1 = overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec2 = overlap_triangle(pos, dist_max, p3, p4, p2, r3, r4, r2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
template <typename T, typename I>
constexpr kernel pintersection<T, 3> overlap_quad(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 4>& quad) {
  return overlap_quad(pos, dist_max, positions[quad.x], positions[quad.y],
      positions[quad.z], positions[quad.w], radius[quad.x], radius[quad.y],
      radius[quad.z], radius[quad.w]);
}

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel bool overlap_bbox(
    const vec<T, 3>& pos, T dist_max, const bbox<T, 3>& bbox) {
  // computing distance
  auto dd = (T)0.0;

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
template <typename T>
constexpr kernel bool overlap_bbox(
    const bbox<T, 3>& bbox1, const bbox<T, 3>& bbox2) {
  for (auto a : range(3)) {
    if (bbox1.max[a] < bbox2.min[a] || bbox1.min[a] > bbox2.max[a])
      return false;
  }
  return true;
}

}  // namespace yocto

#ifndef __CUDACC__

// -----------------------------------------------------------------------------
// GENERIC BVH FOR INTERSECTION AND OVERLAP
// -----------------------------------------------------------------------------
namespace yocto {

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
template <typename T, size_t N>
struct bvh_gnode {
  bbox<T, N> bbox_    = {};
  int32_t    start    = 0;
  int16_t    num      = 0;
  int8_t     axis     = 0;
  bool       internal = false;
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// Application data is not stored explicitly.
template <typename T, size_t N>
struct bvh_gdata {
  vector<bvh_gnode<T, N>> nodes      = {};
  vector<int>             primitives = {};
};

// Typedefs
using bvh_node = bvh_gnode<float, 3>;
using bvh_data = bvh_gdata<float, 3>;

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
template <typename T, size_t N>
struct intersection {
  int       instance = -1;
  int       element  = -1;
  vec<T, 2> uv       = {0, 0};
  T         distance = 0;
  bool      hit      = false;

  intersection() : hit{false} {}
  intersection(int element_, const vec<T, 2>& uv_, T distance_) :
      element{element_}, uv{uv_}, distance{distance_}, hit{true} {}
  intersection(int instance_, int element_, const vec<T, 2>& uv_, T distance_) :
      instance{instance_},
      element{element_},
      uv{uv_},
      distance{distance_},
      hit{true} {}
  intersection(int element_, const pintersection<T, 3>& intersection_) :
      element{element_},
      uv{intersection_.uv},
      distance{intersection_.distance},
      hit{intersection_.hit} {}
  intersection(int instance_, const intersection& intersection_) :
      instance{instance_},
      element{intersection_.element},
      uv{intersection_.uv},
      distance{intersection_.distance},
      hit{intersection_.hit} {}
};

// Typedefs
using intersection3f = intersection<float, 3>;

// Build BVH nodes
template <typename E, typename Func,
    typename T = decltype(std::declval<result_t<Func, E>>().min.x)>
inline bvh_gdata<T, 3> make_bvh(
    const vector<E>& elements, bool highquality, Func&& bbox_func);

// Update bvh
template <typename T, typename E, typename Func>
inline void refit_bvh(
    bvh_gdata<T, 3>& bvh, const vector<E>& elements, Func&& bbox_func);

template <typename T, typename E, typename Func>
inline intersection<T, 3> intersect_bvh(const bvh_gdata<T, 3>& bvh,
    const vector<E>& elements, const ray<T, 3>& ray_, bool find_any,
    Func&& intersect_element);

template <typename T, typename E, typename Func>
inline intersection<T, 3> overlap_bvh(const bvh_gdata<T, 3>& bvh,
    const vector<E>& elements, const ray<T, 3>& pos, float max_distance,
    bool find_any, Func&& overlap_element);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH FOR COLLECTION OF PRIMITIVES
// -----------------------------------------------------------------------------
namespace yocto {

// Build elements bvh
inline bvh_data make_points_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
inline bvh_data make_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
inline bvh_data make_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius);
inline bvh_data make_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius);

// Refit elements bvh
inline void update_points_bvh(bvh_data& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
inline void update_lines_bvh(bvh_data& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
inline void update_triangles_bvh(bvh_data& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
inline void update_quads_bvh(
    bvh_data& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions);

// Bvh intersection
inline intersection3f intersect_points_bvh(const bvh_data& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
inline intersection3f intersect_lines_bvh(const bvh_data& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
inline intersection3f intersect_triangles_bvh(const bvh_data& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);
inline intersection3f intersect_quads_bvh(const bvh_data& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR GENERIC BVH
// -----------------------------------------------------------------------------
namespace yocto {

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
inline vec2i split_sah(
    vector<int>& primitives, const vector<bbox3f>& bboxes, int start, int end) {
  // compute primitive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, bbox_center(bboxes[primitives[i]]));
  auto csize = bbox_diagonal(cbbox);
  if (csize == vec3f{0, 0, 0}) return {(start + end) / 2, 0};

  // consider N bins, compute their cost and keep the minimum
  auto      axis     = 0;
  const int nbins    = 16;
  auto      split    = 0.0f;
  auto      min_cost = flt_max;
  for (auto saxis : range(3)) {
    for (auto b = 1; b < nbins; b++) {
      auto bsplit    = cbbox.min[saxis] + b * csize[saxis] / nbins;
      auto left_bbox = invalidb3f, right_bbox = invalidb3f;
      auto left_nprims = 0, right_nprims = 0;
      for (auto i = start; i < end; i++) {
        if (bbox_center(bboxes[primitives[i]])[saxis] < bsplit) {
          left_bbox = merge(left_bbox, bboxes[primitives[i]]);
          left_nprims += 1;
        } else {
          right_bbox = merge(right_bbox, bboxes[primitives[i]]);
          right_nprims += 1;
        }
      }
      auto cost =
          1 + left_nprims * bbox_area(left_bbox) / (bbox_area(cbbox) + 1e-12f) +
          right_nprims * bbox_area(right_bbox) / (bbox_area(cbbox) + 1e-12f);
      if (cost < min_cost) {
        min_cost = cost;
        split    = bsplit;
        axis     = saxis;
      }
    }
  }
  // split
  auto middle =
      (int)(std::partition(primitives.data() + start, primitives.data() + end,
                [axis, split, &bboxes](auto primitive) {
                  return bbox_center(bboxes[primitive])[axis] < split;
                }) -
            primitives.data());

  // if we were not able to split, just break the primitives in half
  if (middle == start || middle == end) return {(start + end) / 2, axis};

  // done
  return {middle, axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and
// axis.
inline vec2i split_balanced(
    vector<int>& primitives, const vector<bbox3f>& bboxes, int start, int end) {
  // compute primitives bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, bbox_center(bboxes[primitives[i]]));
  auto csize = bbox_diagonal(cbbox);
  if (csize == vec3f{0, 0, 0}) return {(start + end) / 2, 0};

  // split along largest
  auto axis = (int)argmax(csize);

  // balanced tree split: find the largest axis of the
  // bounding box and split along this one right in the middle
  auto middle = (start + end) / 2;
  std::nth_element(primitives.data() + start, primitives.data() + middle,
      primitives.data() + end,
      [axis, &bboxes](auto primitive_a, auto primitive_b) {
        return bbox_center(bboxes[primitive_a])[axis] <
               bbox_center(bboxes[primitive_b])[axis];
      });

  // if we were not able to split, just break the primitives in half
  if (middle == start || middle == end) return {(start + end) / 2, axis};

  // done
  return {middle, axis};
}

// Splits a BVH node using the middle heuristic. Returns split position and
// axis.
inline vec2i split_middle(
    vector<int>& primitives, const vector<bbox3f>& bboxes, int start, int end) {
  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, bbox_center(bboxes[primitives[i]]));
  auto csize = bbox_diagonal(cbbox);
  if (csize == vec3f{0, 0, 0}) return {(start + end) / 2, 0};

  // split along largest
  auto axis = (int)argmax(csize);

  // split the space in the middle along the largest axis
  auto split = bbox_center(cbbox)[axis];
  auto middle =
      (int)(std::partition(primitives.data() + start, primitives.data() + end,
                [axis, split, &bboxes](auto primitive) {
                  return bbox_center(bboxes[primitive])[axis] < split;
                }) -
            primitives.data());

  // if we were not able to split, just break the primitives in half
  if (middle == start || middle == end) return {(start + end) / 2, axis};

  // done
  return {middle, axis};
}

// Maximum number of primitives per BVH node.
constexpr auto bvh_max_prims = 4;

// Build BVH nodes
template <typename E, typename Func, typename T>
inline bvh_gdata<T, 3> make_bvh(
    const vector<E>& elements, bool highquality, Func&& bbox_func) {
  // create primitives

  // bvh
  auto bvh = bvh_gdata<T, 3>{};

  // prepare to build nodes
  bvh.nodes.clear();
  bvh.nodes.reserve(elements.size() * 2);

  // prepare bboxes
  auto bboxes = vector<bbox<T, 3>>(elements.size());
  for (auto&& [bbox, element] : zip(bboxes, elements))
    bbox = bbox_func(element);

  // prepare primitives
  bvh.primitives = vector<int>(bboxes.size());
  for (auto&& [idx, primitive] : enumerate(bvh.primitives))
    primitive = (int)idx;

  // push first node onto the stack
  auto stack = vector<vec3i>{{0, 0, (int)bboxes.size()}};
  bvh.nodes.emplace_back();

  // create nodes until the stack is empty
  while (!stack.empty()) {
    // grab node to work on
    auto [nodeid, start, end] = stack.back();
    stack.pop_back();

    // grab node
    auto& node = bvh.nodes[nodeid];

    // compute bounds
    node.bbox_ = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox_ = merge(node.bbox_, bboxes[bvh.primitives[i]]);

    // split into two children
    if (end - start > bvh_max_prims) {
      // get split
      auto [mid, axis] = highquality
                             ? split_sah(bvh.primitives, bboxes, start, end)
                             : split_middle(bvh.primitives, bboxes, start, end);

      // make an internal node
      node.internal = true;
      node.axis     = (uint8_t)axis;
      node.num      = 2;
      node.start    = (int)bvh.nodes.size();
      bvh.nodes.emplace_back();
      bvh.nodes.emplace_back();
      stack.push_back({node.start + 0, start, mid});
      stack.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = (int16_t)(end - start);
      node.start    = start;
    }
  }

  // cleanup
  bvh.nodes.shrink_to_fit();

  // done
  return bvh;
}

// Update bvh
template <typename T, typename E, typename Func>
inline void refit_bvh(
    bvh_gdata<T, 3>& bvh, const vector<E>& elements, Func&& bbox_func) {
  for (auto nodeid = (int)bvh.nodes.size() - 1; nodeid >= 0; nodeid--) {
    auto& node = bvh.nodes[nodeid];
    node.bbox_ = invalidb3f;
    if (node.internal) {
      for (auto idx : range(2)) {
        node.bbox_ = merge(node.bbox_, bvh.nodes[node.start + idx].bbox_);
      }
    } else {
      for (auto idx : range(node.start, node.start + node.num)) {
        node.bbox_ = merge(
            node.bbox_, bbox_func(elements[bvh.primitives[idx]]));
      }
    }
  }
}

template <typename T, typename E, typename Func>
inline intersection<T, 3> intersect_bvh(const bvh_gdata<T, 3>& bvh,
    const vector<E>& elements, const ray<T, 3>& ray_, bool find_any,
    Func&& intersect_element) {
  // check empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 128>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto sintersection = intersection<T, 3>{};

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = 1 / ray.d;
  auto ray_dsign = component_less(ray_dinv, 0);

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox_)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto eintersection = intersect_element(
            ray, elements[bvh.primitives[idx]]);
        if (!eintersection.hit) continue;
        sintersection = {bvh.primitives[idx], eintersection};
        ray.tmax      = eintersection.distance;
      }
    }

    // check for early exit
    if (find_any && sintersection.hit) return sintersection;
  }

  return sintersection;
}

// Check if a point overlaps some elements.
template <typename T, typename E, typename Func>
inline intersection<T, 3> overlap_bvh(const bvh_gdata<T, 3>& bvh,
    const vector<E>& elements, const vec<T, 3>& pos, float max_distance,
    bool find_any, Func&& overlap_element) {
  // check if empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 64>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // intersection
  auto sintersection = intersection<T, 3>{};

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox_)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx : range(node.start, node.start + node.num)) {
        auto eintersection = overlap_element(
            pos, max_distance, elements[bvh.primitives[idx]]);
        if (!eintersection.hit) continue;
        sintersection = {bvh.primitives[idx], eintersection};
        max_distance  = eintersection.distance;
      }
    }

    // check for early exit
    if (find_any && sintersection.hit) return sintersection;
  }

  return sintersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR ELEMENTS BVH
// -----------------------------------------------------------------------------
namespace yocto {

// Build shape bvh
inline bvh_data make_points_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool highquality = false) {
  return make_bvh(points, highquality,
      [&](int point) { return point_bounds(positions, radius, point); });
}
inline bvh_data make_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool highquality = false) {
  return make_bvh(lines, highquality,
      [&](const vec2i& line) { return line_bounds(positions, radius, line); });
}
inline bvh_data make_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, bool highquality = false) {
  return make_bvh(triangles, highquality, [&](const vec3i& triangle) {
    return triangle_bounds(positions, triangle);
  });
}
inline bvh_data make_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, bool highquality = false) {
  return make_bvh(quads, highquality,
      [&](const vec4i& quad) { return quad_bounds(positions, quad); });
}

inline void update_points_bvh(bvh_data& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  refit_bvh(bvh, points,
      [&](int point) { return point_bounds(positions, radius, point); });
}
inline void update_lines_bvh(bvh_data& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  refit_bvh(bvh, lines,
      [&](const vec2i& line) { return line_bounds(positions, radius, line); });
}
inline void update_triangles_bvh(bvh_data& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  refit_bvh(bvh, triangles, [&](const vec3i& triangle) {
    return triangle_bounds(positions, triangle);
  });
}
inline void update_quads_bvh(
    bvh_data& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions) {
  refit_bvh(bvh, quads,
      [&](const vec4i& quad) { return quad_bounds(positions, quad); });
}

inline intersection3f intersect_points_bvh(const bvh_data& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  return intersect_bvh(
      bvh, points, ray, find_any, [&](const ray3f& ray, const int& point) {
        return intersect_point(ray, positions, radius, point);
      });
}
inline intersection3f intersect_lines_bvh(const bvh_data& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any) {
  return intersect_bvh(
      bvh, lines, ray, find_any, [&](const ray3f& ray, const vec2i& line) {
        return intersect_line(ray, positions, radius, line);
      });
}
inline intersection3f intersect_triangles_bvh(const bvh_data& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  return intersect_bvh(bvh, triangles, ray, find_any,
      [&](const ray3f& ray, const vec3i& triangle) {
        return intersect_triangle(ray, positions, triangle);
      });
}
inline intersection3f intersect_quads_bvh(const bvh_data& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any) {
  return intersect_bvh(
      bvh, quads, ray, find_any, [&](const ray3f& ray, const vec4i& quad) {
        return intersect_quad(ray, positions, quad);
      });
}

}  // namespace yocto

#endif

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
