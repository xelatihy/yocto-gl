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
struct bbox {
  vec<T, N> min = {num_max<T>};
  vec<T, N> max = {num_min<T>};

  constexpr kernel bbox() : min{num_max<T>}, max{num_min<T>} {}
  constexpr kernel bbox(const vec<T, N>& min_, const vec<T, N>& max_)
      : min{min_}, max{max_} {}

  constexpr kernel vec<T, N>& operator[](size_t i) {
    return i == 0 ? min : max;
  }
  constexpr kernel const vec<T, N>& operator[](size_t i) const {
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
[[deprecated]] constexpr kernel vec<T, N> center(const bbox<T, N>& a) {
  return (a.min + a.max) / 2;
}
template <typename T, size_t N>
[[deprecated]] constexpr kernel vec<T, N> size(const bbox<T, N>& a) {
  return a.max - a.min;
}
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
  auto diagonal = a.max - a.min;
  return 2 * diagonal.x * diagonal.y + 2 * diagonal.x * diagonal.z +
         2 * diagonal.y * diagonal.z;
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
struct ray {
  vec<T, N> o    = {0};
  vec<T, N> d    = {0};  // TODO: initialize this properly
  T         tmin = ray_eps<T>;
  T         tmax = num_max<T>;

  constexpr kernel ray() : o{0}, d{0}, tmin{ray_eps<T>}, tmax{num_max<T>} {}
  constexpr kernel ray(const vec<T, N>& o_, const vec<T, N>& d_,
      T tmin_ = ray_eps<T>, T tmax_ = num_max<T>)
      : o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}
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
  if constexpr (N == 3) {
    auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
        vec<T, 3>{b.min.x, b.min.y, b.max.z},
        vec<T, 3>{b.min.x, b.max.y, b.min.z},
        vec<T, 3>{b.min.x, b.max.y, b.max.z},
        vec<T, 3>{b.max.x, b.min.y, b.min.z},
        vec<T, 3>{b.max.x, b.min.y, b.max.z},
        vec<T, 3>{b.max.x, b.max.y, b.min.z},
        vec<T, 3>{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox<T, N>();
    for (auto& corner : corners)
      xformed = merge(xformed, transform_point(a, corner));
    return xformed;
  }
}
template <typename T, size_t N>
constexpr kernel bbox<T, N> transform_bbox(
    const frame<T, N>& a, const bbox<T, N>& b) {
  if constexpr (N == 3) {
    auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
        vec<T, 3>{b.min.x, b.min.y, b.max.z},
        vec<T, 3>{b.min.x, b.max.y, b.min.z},
        vec<T, 3>{b.min.x, b.max.y, b.max.z},
        vec<T, 3>{b.max.x, b.min.y, b.min.z},
        vec<T, 3>{b.max.x, b.min.y, b.max.z},
        vec<T, 3>{b.max.x, b.max.y, b.min.z},
        vec<T, 3>{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox<T, N>();
    for (auto& corner : corners)
      xformed = merge(xformed, transform_point(a, corner));
    return xformed;
  }
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
  auto [v1, v2] = line;
  return line_bounds(positions[v1], positions[v2]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(const vector<vec<T, N>>& positions,
    const vector<T>& radius, const vec<I, 2>& line) {
  auto [v1, v2] = line;
  return line_bounds(positions[v1], positions[v2], radius[v1], radius[v2]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> triangle_bounds(
    const vector<vec<T, N>>& positions, const vec<I, 3>& triangle) {
  auto [v1, v2, v3] = triangle;
  return triangle_bounds(positions[v1], positions[v2], positions[v3]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> quad_bounds(
    const vector<vec<T, N>>& positions, const vec<I, 4>& quad) {
  auto [v1, v2, v3, v4] = quad;
  return quad_bounds(
      positions[v1], positions[v2], positions[v3], positions[v4]);
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
  auto [v1, v2] = line;
  return line_bounds(positions[v1], positions[v2]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> line_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, vec<I, 2> line) {
  auto [v1, v2] = line;
  return line_bounds(positions[v1], positions[v2], radius[v1], radius[v2]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> triangle_bounds(
    cspan<vec<T, N>> positions, vec<I, 3> triangle) {
  auto [v1, v2, v3] = triangle;
  return triangle_bounds(positions[v1], positions[v2], positions[v3]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> quad_bounds(
    cspan<vec<T, N>> positions, vec<I, 4> quad) {
  auto [v1, v2, v3, v4] = quad;
  return quad_bounds(
      positions[v1], positions[v2], positions[v3], positions[v4]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> sphere_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, I point) {
  return sphere_bounds(positions[point], radius[point]);
}
template <typename T, size_t N, typename I>
constexpr kernel bbox<T, N> capsule_bounds(
    cspan<vec<T, N>> positions, cspan<T> radius, vec<I, 2> line) {
  auto [v1, v2] = line;
  return capsule_bounds(positions[v1], positions[v2], radius[v1], radius[v2]);
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
  auto [v1, v2] = line;
  return line_tangent(positions[v1], positions[v2]);
}
template <typename T>
constexpr kernel T line_length(const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return length(p2 - p1);
}
template <typename T, typename I>
constexpr kernel T line_length(
    const vector<vec<T, 3>>& positions, const vec<I, 2>& line) {
  auto [v1, v2] = line;
  return line_length(positions[v1], positions[v2]);
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
  auto [v1, v2, v3] = triangle;
  return triangle_normal(positions[v1], positions[v2], positions[v3]);
}
template <typename T>
constexpr kernel T triangle_area(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
  return length(cross(p2 - p1, p3 - p1)) / 2;
}
template <typename T, typename I>
constexpr kernel T triangle_area(
    const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  auto [v1, v2, v3] = triangle;
  return triangle_area(positions[v1], positions[v2], positions[v3]);
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
  auto [v1, v2, v3, v4] = quad;
  return quad_normal(
      positions[v1], positions[v2], positions[v3], positions[v4]);
}
template <typename T>
constexpr kernel T quad_area(const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 3>& p3, const vec<T, 3>& p4) {
  return triangle_area(p1, p2, p4) + triangle_area(p3, p4, p2);
}
template <typename T, typename I>
constexpr kernel T quad_area(
    const vector<vec<T, 3>>& positions, const vec<I, 4>& quad) {
  auto [v1, v2, v3, v4] = quad;
  return quad_area(positions[v1], positions[v2], positions[v3], positions[v4]);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
constexpr kernel T interpolate_line(const T& p1, const T& p2, T1 u) {
  return p1 * (1 - u) + p2 * u;
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_line(
    const vector<T>& vertices, const vec<I, 2>& line, T1 u) {
  auto [v1, v2] = line;
  return interpolate_line(vertices[v1], vertices[v2], u);
}
// Interpolates values over a triangle parameterized by u and v along the
// (p2-p1) and (p3-p1) directions. Same as barycentric interpolation.
template <typename T, typename T1>
constexpr kernel T interpolate_triangle(
    const T& p1, const T& p2, const T& p3, const vec<T1, 2>& uv) {
  auto [u, v] = uv;
  return p1 * (1 - u - v) + p2 * u + p3 * v;
}
template <typename T, typename T1, typename I>
constexpr kernel T interpolate_triangle(
    const vector<T>& vertices, const vec<I, 3>& triangle, T1 u) {
  auto [v1, v2, v3] = triangle;
  return interpolate_triangle(vertices[v1], vertices[v2], vertices[v3], u);
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
  auto [v1, v2, v3, v4] = quad;
  return interpolate_quad(
      vertices[v1], vertices[v2], vertices[v3], vertices[v4], u);
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
  auto [v1, v2] = line;
  return line_point(positions[v1], positions[v2], u);
}
template <typename T>
constexpr kernel vec<T, 3> line_tangent(
    const vec<T, 3>& t0, const vec<T, 3>& t1, T u) {
  return normalize(t0 * (1 - u) + t1 * u);
}
template <typename T, typename I>
constexpr kernel vec<T, 3> line_tangent(
    const vector<vec<T, 3>>& tangents, const vec<I, 2>& line, T u) {
  auto [v1, v2] = line;
  return line_tangent(tangents[v1], tangents[v2], u);
}

// Interpolated triangle properties.
template <typename T>
constexpr kernel vec<T, 3> triangle_point(const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3, const vec<T, 2>& uv) {
  auto [u, v] = uv;
  return p1 * (1 - u - v) + p2 * u + p3 * v;
}
template <typename T, typename I>
constexpr kernel vec<T, 3> triangle_point(const vector<vec<T, 3>>& positions,
    const vec<I, 3>& triangle, const vec<T, 2>& uv) {
  auto [v1, v2, v3] = triangle;
  return triangle_point(positions[v1], positions[v2], positions[v3], uv);
}
template <typename T>
constexpr kernel vec<T, 3> triangle_normal(const vec<T, 3>& n0,
    const vec<T, 3>& n1, const vec<T, 3>& n2, const vec<T, 2>& uv) {
  auto [u, v] = uv;
  return normalize(n0 * (1 - u - v) + n1 * u + n2 * v);
}
template <typename T, typename I>
constexpr kernel vec<T, 3> triangle_normal(const vector<vec<T, 3>>& normals,
    const vec<I, 3>& triangle, const vec<T, 2>& uv) {
  auto [v1, v2, v3] = triangle;
  return triangle_normal(normals[v1], normals[v2], normals[v3], uv);
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
  auto [v1, v2, v3, v4] = quad;
  return quad_point(
      positions[v1], positions[v2], positions[v3], positions[v4], uv);
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
  auto [v1, v2, v3, v4] = quad;
  return quad_normal(normals[v1], normals[v2], normals[v3], normals[v4], uv);
}

// Interpolated sphere properties.
template <typename T>
constexpr kernel vec<T, 3> sphere_point(
    const vec<T, 3> p, T r, const vec<T, 2>& uv) {
  return p + r * vec<T, 3>{cos(uv.x * 2 * (T)pi) * sin(uv.y * (T)pi),
                     sin(uv.x * 2 * (T)pi) * sin(uv.y * (T)pi),
                     cos(uv.y * (T)pi)};
}
template <typename T>
constexpr kernel vec<T, 3> sphere_normal(
    const vec<T, 3> p, T r, const vec<T, 2>& uv) {
  return normalize(vec<T, 3>{cos(uv.x * 2 * (T)pi) * sin(uv.y * (T)pi),
      sin(uv.x * 2 * (T)pi) * sin(uv.y * (T)pi), cos(uv.y * (T)pi)});
}

// Triangle tangent and bi-tangent from uv
template <typename T>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  auto p = p2 - p1, q = p3 - p1;
  auto s   = vec<T, 2>{uv1.x - uv0.x, uv2.x - uv0.x};
  auto t   = vec<T, 2>{uv1.y - uv0.y, uv2.y - uv0.y};
  auto div = cross(s, t);

  if (div != 0) {
    auto tu = vec<T, 3>{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                  t.y * p.z - t.x * q.z} /
              div;
    auto tv = vec<T, 3>{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                  s.x * q.z - s.y * p.z} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}
template <typename T, typename I>
constexpr kernel pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 2>>& texcoords,
    const vec<I, 4>& triangle, const vec<T, 2>& uv) {
  auto [v1, v2, v3] = triangle;
  return triangle_tangents_fromuv(positions[v1], positions[v2], positions[v3],
      texcoords[v1], texcoords[v2], texcoords[v3], uv);
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
    const vec<I, 4>& quad, const vec<T, 2>& uv) {
  auto [v1, v2, v3, v4] = quad;
  return quad_tangents_fromuv(positions[v1], positions[v2], positions[v3],
      positions[v4], texcoords[v1], texcoords[v2], texcoords[v3], texcoords[v4],
      uv);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generate a ray from a camera
template <typename T>
constexpr kernel ray<T, 3> camera_ray(const frame<T, 3>& frame, T lens,
    const vec<T, 2>& film, const vec<T, 2>& image_uv) {
  auto uv  = image_uv * vec{-1, 1};  // flip x
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
  auto uv   = image_uv * vec{-1, 1};  // flip x
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
template <typename T = float>
struct prim_gintersection {
  vec<T, 2> uv       = {0, 0};
  T         distance = num_max<T>;  // TODO: num_max<T>
  bool      hit      = false;
};

// Typedefs
using prim_intersection = prim_gintersection<float>;

// Intersect a ray with a point (approximate)
template <typename T>
constexpr kernel prim_gintersection<T> intersect_point(
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
  return {{0, 0}, t, true};
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> intersect_point(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vector<T>& radius, I point) {
  auto v1 = point;
  return intersect_point(ray, positions[v1], radius[v1]);
}

// Intersect a ray with a line
template <typename T>
constexpr kernel prim_gintersection<T> intersect_line(const ray<T, 3>& ray,
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
  return {{s, sqrt(d2) / r}, t, true};
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> intersect_line(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vector<T>& radius,
    const vec<I, 2>& line) {
  auto [v1, v2] = line;
  return intersect_line(
      ray, positions[v1], positions[v2], radius[v1], radius[v2]);
}

// Intersect a ray with a sphere
template <typename T>
constexpr kernel prim_gintersection<T> intersect_sphere(
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
  auto plocal = ((ray.o + ray.d * t) - p) / r;
  auto u      = atan2(plocal.y, plocal.x) / (2 * (T)pi);
  if (u < 0) u += 1;
  auto v = acos(clamp(plocal.z, (T)-1, (T)1)) / (T)pi;

  // intersection occurred: set params and exit
  return {{u, v}, t, true};
}

// Intersect a ray with a triangle
template <typename T>
constexpr kernel prim_gintersection<T> intersect_triangle(const ray<T, 3>& ray,
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

  // compute and check first bricentric coordinated
  auto tvec = ray.o - p1;
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
template <typename T, typename I>
constexpr kernel prim_gintersection<T> intersect_triangle(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  auto [v1, v2, v3] = triangle;
  return intersect_triangle(ray, positions[v1], positions[v2], positions[v3]);
}

// Intersect a ray with a quad.
template <typename T>
constexpr kernel prim_gintersection<T> intersect_quad(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4) {
  if (p3 == p4) return intersect_triangle(ray, p1, p2, p4);
  auto isec1 = intersect_triangle(ray, p1, p2, p4);
  auto isec2 = intersect_triangle(ray, p3, p4, p2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> intersect_quad(const ray<T, 3>& ray,
    const vector<vec<T, 3>>& positions, const vec<I, 4>& quad) {
  auto [v1, v2, v3, v4] = quad;
  return intersect_quad(
      ray, positions[v1], positions[v2], positions[v3], positions[v4]);
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
constexpr kernel prim_gintersection<T> overlap_point(
    const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p, T r) {
  auto d2 = dot(pos - p, pos - p);
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  return {{0, 0}, sqrt(d2), true};
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> overlap_point(const vec<T, 3>& pos,
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
constexpr kernel prim_gintersection<T> overlap_line(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, T r1, T r2) {
  auto u = closestuv_line(pos, p1, p2);
  // Compute projected position from the clamped t d = a + t * ab;
  auto p  = p1 + (p2 - p1) * u;
  auto r  = r1 + (r2 - r1) * u;
  auto d2 = dot(pos - p, pos - p);
  // check distance
  if (d2 > (dist_max + r) * (dist_max + r)) return {};
  // done
  return {{u, 0}, sqrt(d2), true};
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> overlap_line(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 2>& line) {
  auto [v1, v2] = line;
  return overlap_line(
      pos, dist_max, positions[v1], positions[v2], radius[v1], radius[v2]);
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
constexpr kernel prim_gintersection<T> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vec<I, 3>& triangle) {
  auto [v1, v2, v3] = triangle;
  return overlap_triangle(
      pos, dist_max, positions[v1], positions[v2], positions[v3]);
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
template <typename T>
constexpr kernel prim_gintersection<T> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    T r1, T r2, T r3) {
  auto uv = closestuv_triangle(pos, p1, p2, p3);
  auto p  = interpolate_triangle(p1, p2, p3, uv);
  auto r  = interpolate_triangle(r1, r2, r3, uv);
  auto dd = dot(p - pos, p - pos);
  if (dd > (dist_max + r) * (dist_max + r)) return {};
  return {uv, sqrt(dd), true};
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 3>& triangle) {
  auto [v1, v2, v3] = triangle;
  return overlap_triangle(pos, dist_max, positions[v1], positions[v2],
      positions[v3], radius[v1], radius[v2], radius[v3]);
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel prim_gintersection<T> overlap_quad(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4, T r1, T r2, T r3, T r4) {
  if (p3 == p4) return overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec1 = overlap_triangle(pos, dist_max, p1, p2, p4, r1, r2, r3);
  auto isec2 = overlap_triangle(pos, dist_max, p3, p4, p2, r3, r4, r2);
  if (isec2.hit) isec2.uv = 1 - isec2.uv;
  return isec1.distance < isec2.distance ? isec1 : isec2;
}
template <typename T, typename I>
constexpr kernel prim_gintersection<T> overlap_quad(const vec<T, 3>& pos,
    T dist_max, const vector<vec<T, 3>>& positions, const vector<T> radius,
    const vec<I, 4>& quad) {
  auto [v1, v2, v3, v4] = quad;
  return overlap_quad(pos, dist_max, positions[v1], positions[v2],
      positions[v3], positions[v4], radius[v1], radius[v2], radius[v3],
      radius[v4]);
}

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
template <typename T>
constexpr kernel bool overlap_bbox(
    const vec<T, 3>& pos, T dist_max, const bbox<T, 3>& bbox) {
  // computing distance
  auto dd = (T)0.0;

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
template <typename T>
constexpr kernel bool overlap_bbox(
    const bbox<T, 3>& bbox1, const bbox<T, 3>& bbox2) {
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
template <typename T>
[[deprecated]] constexpr kernel bool intersect_point(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r, vec<T, 2>& uv, T& dist) {
  auto intersection = intersect_point(ray, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a line
template <typename T>
[[deprecated]] constexpr kernel bool intersect_line(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, T r1, T r2, vec<T, 2>& uv,
    T& dist) {
  auto intersection = intersect_line(ray, p1, p2, r1, r2);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a sphere
template <typename T>
[[deprecated]] constexpr kernel bool intersect_sphere(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r, vec<T, 2>& uv, T& dist) {
  auto intersection = intersect_sphere(ray, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a triangle
template <typename T>
[[deprecated]] constexpr kernel bool intersect_triangle(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    vec<T, 2>& uv, T& dist) {
  auto intersection = intersect_triangle(ray, p1, p2, p3);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Intersect a ray with a quad.
template <typename T>
[[deprecated]] constexpr kernel bool intersect_quad(const ray<T, 3>& ray,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4, vec<T, 2>& uv, T& dist) {
  auto intersection = intersect_quad(ray, p1, p2, p3, p4);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a point overlaps a position pos withint a maximum distance dist_max.
template <typename T>
[[deprecated]] constexpr kernel bool overlap_point(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p, T r, vec<T, 2>& uv, T& dist) {
  auto intersection = overlap_point(pos, dist_max, p, r);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
template <typename T>
[[deprecated]] constexpr kernel bool overlap_line(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, T r1, T r2,
    vec<T, 2>& uv, T& dist) {
  auto intersection = overlap_line(pos, dist_max, p1, p2, r1, r2);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
template <typename T>
[[deprecated]] constexpr kernel bool overlap_triangle(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    T r1, T r2, T r3, vec<T, 2>& uv, T& dist) {
  auto intersection = overlap_triangle(pos, dist_max, p1, p2, p3, r1, r2, r3);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
  return true;
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
template <typename T>
[[deprecated]] constexpr kernel bool overlap_quad(const vec<T, 3>& pos,
    T dist_max, const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 3>& p4, T r1, T r2, T r3, T r4, vec<T, 2>& uv, T& dist) {
  auto intersection = overlap_quad(
      pos, dist_max, p1, p2, p3, p4, r1, r2, r3, r4);
  if (!intersection.hit) return false;
  uv   = intersection.uv;
  dist = intersection.distance;
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
