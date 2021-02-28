//
// Implementation for Yocto/Mesh
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

// TODO(fabio): remove asserts
// TODO(fabio): better name for v2t
// TODO(fabio): better name for adjacencies

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_mesh.h"

#include <cassert>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>

#include "yocto_geometry.h"
#include "yocto_modelio.h"

// #define TESTS_MAY_FAIL

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::deque;
using std::pair;
using std::unordered_set;
using namespace std::string_literals;

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

}  // namespace std

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// find a value in a vector or vecs
static int find_in_vec(const vector<int>& vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x) return i;
  return -1;
}

static int find_in_vec(const vec3i& vec, int x) {
  for (auto i = 0; i < 3; i++)
    if (vec[i] == x) return i;
  return -1;
}

static int mod3(int i) { return (i > 2) ? i - 3 : i; }

static vec3f get_bary(const vec2f& uv) {
  return vec3f{1 - uv.x - uv.y, uv.x, uv.y};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

pair<bool, int> bary_is_vert(const vec3f& bary, float tol) {
  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

pair<bool, int> bary_is_edge(const vec3f& bary, float tol) {
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

pair<bool, int> point_is_vert(const mesh_point& p, float tol) {
  auto bary = vec3f{1 - p.uv.x - p.uv.y, p.uv.x, p.uv.y};
  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

pair<bool, int> point_is_edge(const mesh_point& p, float tol) {
  auto bary = vec3f{1 - p.uv.x - p.uv.y, p.uv.x, p.uv.y};
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

pair<bool, vec2f> point_in_triangle(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int tid, const vec3f& point, float tol) {
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  auto b  = vec3f{0, 0, 0};
  auto v0 = positions[triangles[tid].x];
  auto v1 = positions[triangles[tid].y];
  auto v2 = positions[triangles[tid].z];

  auto u = v1 - v0, v = v2 - v0, w = point - v0;
  auto d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
       d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0) return {false, zero2f};

  b[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(b[2]));
  b[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(b[1]));
  b[0] = 1 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (auto i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol) return {false, zero2f};
  }
  auto uv = vec2f{b.y, b.z};
  uv      = clamp(uv, 0.f, 1.f);
  return {true, uv};
}

vec3f eval_position(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& sample) {
  auto [x, y, z] = triangles[sample.face];
  return interpolate_triangle(
      positions[x], positions[y], positions[z], sample.uv);
}

vec3f eval_normal(const vector<vec3i>& triangles, const vector<vec3f>& normals,
    const mesh_point& point) {
  auto [x, y, z] = triangles[point.face];
  return normalize(
      interpolate_triangle(normals[x], normals[y], normals[z], point.uv));
}

vec2f intersect_circles(const vec2f& c2, float R2, const vec2f& c1, float R1) {
  auto R = length_squared(c2 - c1);
  assert(R > 0);
  auto invR   = 1 / R;
  auto result = (c1 + c2);
  result += (c2 - c1) * ((R1 - R2) * invR);
  auto A = 2 * (R1 + R2) * invR;
  auto B = (R1 - R2) * invR;
  auto s = A - B * B - 1;
  assert(s >= 0);
  result += vec2f{c2.y - c1.y, c1.x - c2.x} * sqrt(s);
  return result / 2;
}

[[maybe_unused]] static vec2f intersect_circles_double(
    float c2x, float c2y, float R2, double c1x, double c1y, float R1) {
  auto R = (c2x - c1x) * (c2x - c1x) + (c2y - c1y) * (c2y - c1y);
  assert(R > 0);
  auto invR    = 1 / R;
  auto resultx = (c1x + c2x);
  auto resulty = (c1y + c2y);
  resultx += (c2x - c1x) * ((R1 - R2) * invR);
  resulty += (c2y - c1y) * ((R1 - R2) * invR);

  auto A = 2 * (R1 + R2) * invR;
  auto B = (R1 - R2) * invR;
  auto s = A - B * B - 1;
  assert(s >= 0);
  resultx += c2y - c1y * yocto::sqrt(s);
  resulty += c1x - c2x * yocto::sqrt(s);
  return vec2f{(float)(resultx / 2), (float)(resulty / 2)};
}

unfold_triangle init_flat_triangle(
    const vector<vec3f>& positions, const vec3i& tr) {
  auto tr2d = unfold_triangle{};
  tr2d[0]   = {0, 0};
  tr2d[1]   = {0, length(positions[tr.x] - positions[tr.y])};
  auto rx   = length_squared(positions[tr.x] - positions[tr.z]);
  auto ry   = length_squared(positions[tr.y] - positions[tr.z]);
  tr2d[2]   = intersect_circles(tr2d[0], rx, tr2d[1], ry);
  return tr2d;
}

static int find_adjacent_triangle(
    const vec3i& triangle, const vec3i& adjacent) {
  for (int i = 0; i < 3; i++) {
    auto k = find_in_vec(adjacent, triangle[i]);
    if (k != -1) {
      if (find_in_vec(adjacent, triangle[mod3(i + 1)]) != -1) {
        return i;
      } else {
        return mod3(i + 2);
      }
    }
  }
  assert(0 && "input triangles are not adjacent");
  return -1;
}

static int find_adjacent_triangle(
    const vector<vec3i>& triangles, int face, int neighbor) {
  return find_adjacent_triangle(triangles[face], triangles[neighbor]);
}

[[maybe_unused]] static vec2f unfold_point(const vec3f& pa, const vec3f& pb,
    const vec3f& pv, const vec2f& ca, const vec2f& cb) {
  // Unfold position of vertex v
  auto ex = normalize(ca - cb);
  auto ey = vec2f{-ex.y, ex.x};

  auto result = vec2f{};

  // Ortogonal projection
  auto pv_pb = pv - pb;
  auto x     = dot(pa - pb, pv_pb) / length(pa - pb);
  result     = x * ex;

  // Pythagorean theorem
  auto y = dot(pv_pb, pv_pb) - x * x;
  assert(y > 0);
  y = yocto::sqrt(y);
  result += y * ey;
  return result;
}

// given the 2D coordinates in tanget space of a triangle, find the coordinates
// of the k-th neighbor triangle
unfold_triangle unfold_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const unfold_triangle& tr, int face,
    int neighbor) {
  auto k = find_adjacent_triangle(triangles, face, neighbor);
  auto j = find_adjacent_triangle(triangles, neighbor, face);
  assert(j != -1);
  assert(k != -1);
  auto v = triangles[neighbor][mod3(j + 2)];
  auto a = triangles[face][k];
  auto b = triangles[face][mod3(k + 1)];

  auto result         = unfold_triangle{};
  result[j]           = tr[mod3(k + 1)];
  result[mod3(j + 1)] = tr[k];

  // TODO(splinesurf): check which unfolding method is better.
#if 1
  // old method
  auto r0             = length_squared(positions[v] - positions[a]);
  auto r1             = length_squared(positions[v] - positions[b]);
  result[mod3(j + 2)] = intersect_circles(
      result[j], r1, result[mod3(j + 1)], r0);
#else
  // new method
  auto& pa    = positions[a];
  auto& pb    = positions[b];
  auto& pv    = positions[v];
  auto  point = unfold_point(pa, pb, pv, result[mod3(j + 1)], result[j]);
  result[mod3(j + 2)] = result[j] + point;
#endif

  assert(result[0] != result[1]);
  assert(result[1] != result[2]);
  assert(result[2] != result[0]);
  return result;
}

static unfold_triangle unfold_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const unfold_triangle& tr, int face, int k) {
  return unfold_face(triangles, positions, tr, face, adjacencies[face][k]);
}

// assign 2D coordinates to vertices of the triangle containing the mesh
// point, putting the point at (0, 0)
unfold_triangle triangle_coordinates(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& point) {
  auto result = unfold_triangle{};
  auto tr     = triangles[point.face];
  result[0]   = {0, 0};
  result[1]   = {0, length(positions[tr.x] - positions[tr.y])};
  auto rx     = length_squared(positions[tr.x] - positions[tr.z]);
  auto ry     = length_squared(positions[tr.y] - positions[tr.z]);
  result[2]   = intersect_circles(result[0], rx, result[1], ry);

  // Transform coordinates such that point = (0, 0)
  auto point_coords = interpolate_triangle(
      result[0], result[1], result[2], point.uv);
  result[0] -= point_coords;
  result[1] -= point_coords;
  result[2] -= point_coords;

  assert(result[0] != result[1]);
  assert(result[1] != result[2]);
  assert(result[2] != result[0]);
  return result;
}

// assign 2D coordinates to a strip of triangles. point start is at (0, 0)
vector<unfold_triangle> unfold_strip(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<int>& strip,
    const mesh_point& start) {
  auto coords = vector<unfold_triangle>(strip.size());
  assert(start.face == strip[0]);
  coords[0] = triangle_coordinates(triangles, positions, start);

  for (auto i = 1; i < strip.size(); i++) {
    assert(coords[i - 1][0] != coords[i - 1][1]);
    assert(coords[i - 1][1] != coords[i - 1][2]);
    assert(coords[i - 1][2] != coords[i - 1][0]);
    coords[i] = unfold_face(
        triangles, positions, coords[i - 1], strip[i - 1], strip[i]);
  }

  return coords;
}

// Create sequence of 2D segments (portals) needed for funneling.
static vector<pair<vec2f, vec2f>> make_funnel_portals(
    const vector<vec3i>& triangles, const vector<unfold_triangle>& coords,
    const vector<int>& strip, const mesh_point& to) {
  auto portals = vector<pair<vec2f, vec2f>>(strip.size());
  for (auto i = 0; i < strip.size() - 1; i++) {
    auto curr = strip[i], next = strip[i + 1];
    auto k     = find_adjacent_triangle(triangles, curr, next);
    auto tr    = coords[i];
    portals[i] = {tr[k], tr[mod3(k + 1)]};
  }
  auto end = interpolate_triangle(
      coords.back()[0], coords.back()[1], coords.back()[2], to.uv);
  portals.back() = {end, end};
  return portals;
}

static vector<pair<vec2f, vec2f>> unfold_funnel_portals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<int>& strip, const mesh_point& start, const mesh_point& end) {
  auto coords = unfold_strip(triangles, positions, strip, start);
  return make_funnel_portals(triangles, coords, strip, end);
}

vec2i get_edge(const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int f0, int f1) {
  auto k = find_in_vec(adjacencies[f0], f1);
  if (k == -1) return {-1, -1};
  auto tr = triangles[f0];
  return vec2i{tr[k], tr[mod3(k + 1)]};
}

// "strip" must be such that strip.back()=p.face and strip[0] must share at
// least one vertex with p.face.
// Under this hypoteses, this function gives back the distance between the
// opposite vertex to the edge between strip[0] and strip[1] and p handling
// concave paths.
// first_sample_pos is the position in the 2D-reference system defined in
// "init_flat_triangle" where the path intersect the edge between strip[0] and
// strip[1].
float length_by_flattening(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& p, const vector<int>& strip) {
  auto opp_pid = strip[0];
  auto h       = find_in_vec(adjacencies[opp_pid], strip[1]);
  auto coords  = vector<unfold_triangle>(strip.size());
  coords[0]    = init_flat_triangle(positions, triangles[opp_pid]);
  for (auto i = 1; i < strip.size(); i++) {
    auto k = find_in_vec(adjacencies[strip[i - 1]], strip[i]);
    assert(k != -1);
    auto tr = unfold_face(
        triangles, positions, adjacencies, coords[i - 1], strip[i - 1], k);
    coords[i] = tr;
  }

  auto last  = coords.back();
  auto pos2d = interpolate_triangle(last[0], last[1], last[2], p.uv);
  auto v     = pos2d - coords[0][(h + 2) % 3];
  auto w0    = coords[0][h] - coords[0][(h + 2) % 3];
  auto w1    = coords[0][(h + 1) % 3] - coords[0][(h + 2) % 3];
  auto phi = angle(w0, w1), theta0 = angle(v, w0), theta1 = angle(v, w1);

  if (theta0 < phi && theta1 < phi) {
    // first_sample_direction = pos2d;
    return length(v);
  } else if (theta1 > phi) {
    // first_sample_direction = coords[0][h];
    return length(w0) + length(coords[0][h] - pos2d);
  } else {
    // first_sample_direction = coords[0][(h + 1) % 3];
    return length(w1) + length(coords[0][(h + 1) % 3] - pos2d);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Triangle fan starting from a face and going towards the k-th neighbor face.
vector<int> triangle_fan(
    const vector<vec3i>& adjacencies, int face, int k, bool clockwise) {
  auto result = vector<int>{};
  result.push_back(face);
  auto prev   = face;
  auto node   = adjacencies[face][k];
  auto offset = 2 - (int)clockwise;
  // for (int i = 0; i < 256; i++) {
  while (true) {
    if (node == -1) break;
    if (node == face) break;
    result.push_back(node);
    auto kk = find_in_vec(adjacencies[node], prev);
    assert(kk != -1);
    kk   = mod3(kk + offset);
    prev = node;
    node = adjacencies[node][kk];
  }
  return result;
}

// returns the list of triangles incident at each vertex in ccw order
vector<vector<int>> vertex_to_triangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies) {
  auto v2t    = vector<vector<int>>{positions.size(), vector<int>{}};
  auto offset = 0;
  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      auto curr = triangles[i][j];
      if (v2t[curr].size() == 0) {
        offset    = find_in_vec(triangles[i], curr);
        v2t[curr] = triangle_fan(adjacencies, i, mod3(offset + 2));
      }
    }
  }
  return v2t;
}

// Face adjacent to t and opposite to vertex vid
int opposite_face(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, int t, int vid) {
  auto triangle = adjacencies[t];
  for (auto i = 0; i < 3; ++i) {
    if (find_in_vec(triangles[triangle[i]], vid) < 0) return triangle[i];
  }
  return -1;
}

// Finds the opposite vertex of an edge
int opposite_vertex(const vec3i& triangle, const vec2i& edge) {
  for (auto i = 0; i < 3; ++i) {
    if (triangle[i] != edge.x && triangle[i] != edge.y) return triangle[i];
  }
  return -1;
}

int opposite_vertex(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, int face, int k) {
  int neighbor = adjacencies[face][k];
  int j        = find_in_vec(adjacencies[neighbor], face);
  assert(j != -1);
  auto tt = triangles[neighbor];
  return tt[mod3(j + 2)];
}

// Finds common edge between triangles
vec2i common_edge(const vec3i& triangle0, const vec3i& triangle1) {
  for (auto i = 0; i < 3; i++) {
    for (auto k = 0; k < 3; k++) {
      if (triangle0[i] == triangle1[k] &&
          triangle0[mod3(i + 1)] == triangle1[mod3(k + 2)])
        return {triangle0[i], triangle0[mod3(i + 1)]};
    }
  }
  return {-1, -1};
}

vec2i opposite_edge(const vec3i& t, int vid) {
  auto offset = find_in_vec(t, vid);
  auto v0     = t[mod3(offset + 1)];
  auto v1     = t[mod3(offset + 2)];
  return vec2i{v0, v1};
}

// TODO: cleanup
vec2i common_edge(const vector<vec3i>& triangles, int pid0, int pid1) {
  auto& poly0 = triangles[pid0];
  auto& poly1 = triangles[pid1];
  for (auto i = 0; i < 3; ++i) {
    auto& vid    = poly0[i];
    auto  offset = find_in_vec(poly1, vid);
    if (offset < 0) {
      offset  = find_in_vec(poly0, vid);
      auto e0 = poly0[(offset + 1) % 3];
      auto e1 = poly0[(offset + 2) % 3];
      if (find_in_vec(poly1, e0) != -1 && find_in_vec(poly1, e1) != -1)
        return {e0, e1};
      else
        return {-1, -1};
    }
  }
  return {-1, -1};
}

// TODO: cleanup
int common_vertex(const vector<vec3i>& triangles, int pid0, int pid1) {
  auto& poly0 = triangles[pid0];
  auto& poly1 = triangles[pid1];
  for (auto i = 0; i < 3; ++i) {
    auto& vid    = poly0[i];
    auto  offset = find_in_vec(poly1, vid);
    if (offset != -1) return vid;
  }
  return -1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

// Extract isoline from surface scalar field.
void meandering_triangles(const vector<float>& field, float isoline,
    int selected_tag, int t0, int t1, vector<vec3i>& triangles,
    vector<int>& tags, vector<vec3f>& positions, vector<vec3f>& normals) {
  auto num_triangles = triangles.size();

  // Edgemap to keep track of the added vertex on each splitted edge.
  // key: edge (ordered std::pair), value: vertex index
  auto emap = unordered_map<vec2i, int>();

  // Helper procedures.
  auto make_edge = [](int a, int b) -> vec2i {
    return a < b ? vec2i{a, b} : vec2i{b, a};
  };
  auto add_vertex = [&](int a, int b, float coeff) -> int {
    auto position = coeff * positions[a] + (1 - coeff) * positions[b];
    auto normal   = normalize(coeff * normals[a] + (1 - coeff) * normals[b]);
    auto index    = (int)positions.size();
    positions.push_back(position);
    normals.push_back(normal);
    return index;
  };
  auto get_tag = [&](int v) { return field[v] > isoline ? t1 : t0; };

  for (auto i = 0; i < num_triangles; ++i) {
    if (tags[i] != selected_tag) continue;

    auto tr = triangles[i];

    // Find which vertex has different tag, if any.
    int j = -1;
    for (auto k = 0; k < 3; k++) {
      if (get_tag(tr[k]) != get_tag(tr[mod3(k + 1)]) &&
          get_tag(tr[k]) != get_tag(tr[mod3(k + 2)])) {
        j = k;
      }
    }

    // If all vertices have same tag, then tag the whole triangle and continue.
    if (j == -1) {
      tags[i] = get_tag(tr.x);
      continue;
    }

    // Reorder the triangle such that vertex with different tag is always z
    tr = {tr[mod3(j + 2)], tr[mod3(j + 1)], tr[j]};

    auto values = vec3f{field[tr.x], field[tr.y], field[tr.z]};
    values -= isoline;

    // Create or retrieve new vertices that split the edges
    auto new_verts = std::array<int, 2>{-1, -1};
    for (auto k : {0, 1}) {
      auto vert = tr[k];
      auto a    = values[k];
      auto b    = values.z;
      auto edge = make_edge(tr.z, vert);
      auto it   = emap.find(edge);

      if (it != emap.end()) {
        // Edge already processed.
        new_verts[k] = it->second;
      } else {
        // Compute new vertex via interpolation.
        auto alpha   = abs(a / (b - a));
        new_verts[k] = add_vertex(tr.z, vert, alpha);
        emap.insert(it, {edge, new_verts[k]});
      }
    }

    /*
                     tr.z
                       /\
                      /  \
                     /    \
                    /  i   \
                   /        \
    new_verts[1]  /..........\  new_verts[0]
                 /   . n_f[0] \
                /       .      \
               /           .    \
              / new_faces[1]  .  \
             /___________________.\
       tr.y                      tr.x

    */

    // Add two new faces.
    triangles.push_back({new_verts[0], new_verts[1], tr.x});
    tags.push_back(get_tag(tr.x));

    triangles.push_back({tr.y, tr.x, new_verts[1]});
    tags.push_back(get_tag(tr.y));

    // Edit old face.
    triangles[i] = vec3i{new_verts[0], tr.z, new_verts[1]};
    tags[i]      = get_tag(tr.z);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto {

static void connect_nodes(geodesic_solver& solver, int a, int b, float length) {
  solver.graph[a].push_back({b, length});
  solver.graph[b].push_back({a, length});
}

static float opposite_nodes_arc_length(
    const vector<vec3f>& positions, int a, int c, const vec2i& edge) {
  // Triangles (a, b, d) and (b, d, c) are connected by (b, d) edge
  // Nodes a and c must be connected.

  auto b = edge.x, d = edge.y;
  auto ba = positions[a] - positions[b];
  auto bc = positions[c] - positions[b];
  auto bd = positions[d] - positions[b];

  auto cos_alpha = dot(normalize(ba), normalize(bd));
  auto cos_beta  = dot(normalize(bc), normalize(bd));
  auto sin_alpha = sqrt(max(0.0f, 1 - cos_alpha * cos_alpha));
  auto sin_beta  = sqrt(max(0.0f, 1 - cos_beta * cos_beta));

  // cos(alpha + beta)
  auto cos_alpha_beta = cos_alpha * cos_beta - sin_alpha * sin_beta;
  if (cos_alpha_beta <= -1) return flt_max;

  // law of cosines (generalized Pythagorean theorem)
  auto len = dot(ba, ba) + dot(bc, bc) -
             length(ba) * length(bc) * 2 * cos_alpha_beta;

  if (len <= 0)
    return flt_max;
  else
    return sqrt(len);
}

static void connect_opposite_nodes(geodesic_solver& solver,
    const vector<vec3f>& positions, const vec3i& tr0, const vec3i& tr1,
    const vec2i& edge) {
  auto opposite_vertex = [](const vec3i& tr, const vec2i& edge) -> int {
    for (auto i = 0; i < 3; ++i) {
      if (tr[i] != edge.x && tr[i] != edge.y) return tr[i];
    }
    return -1;
  };

  auto v0 = opposite_vertex(tr0, edge);
  auto v1 = opposite_vertex(tr1, edge);
  if (v0 == -1 || v1 == -1) return;
  auto length = opposite_nodes_arc_length(positions, v0, v1, edge);
  connect_nodes(solver, v0, v1, length);
}

geodesic_solver make_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vec3f>& positions) {
  auto solver = geodesic_solver{};
  solver.graph.resize(positions.size());
  for (auto face = 0; face < triangles.size(); face++) {
    for (auto k = 0; k < 3; k++) {
      auto a = triangles[face][k];
      auto b = triangles[face][mod3(k + 1)];

      // connect mesh edges
      auto len = length(positions[a] - positions[b]);
      if (a < b) connect_nodes(solver, a, b, len);

      // connect opposite nodes
      auto neighbor = adjacencies[face][k];
      if (face < neighbor) {
        connect_opposite_nodes(
            solver, positions, triangles[face], triangles[neighbor], {a, b});
      }
    }
  }
  return solver;
}

// Construct a graph to compute geodesic distances
dual_geodesic_solver make_dual_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies) {
  auto get_triangle_center = [](const vector<vec3i>&  triangles,
                                 const vector<vec3f>& positions, int face) {
    return (positions[triangles[face].x] + positions[triangles[face].y] +
               positions[triangles[face].z]) /
           3;
  };

  auto solver = dual_geodesic_solver{};
  solver.graph.resize(triangles.size());
  for (auto i = 0; i < solver.graph.size(); ++i) {
    for (auto k = 0; k < 3; ++k) {
      solver.graph[i][k].node = adjacencies[i][k];

      if (adjacencies[i][k] == -1) {
        solver.graph[i][k].length = flt_max;
      } else {
        // TODO(splinesurf): check which solver is better/faster.
#if 1
        solver.graph[i][k].length = length(
            get_triangle_center(triangles, positions, i) -
            get_triangle_center(triangles, positions, adjacencies[i][k]));
#else
        auto p0    = mesh_point{i, {1 / 3.f, 1 / 3.f}};
        auto p1    = mesh_point{adjacencies[i][k], {1 / 3.f, 1 / 3.f}};
        auto path  = geodesic_path{};
        path.strip = {i, adjacencies[i][k]};
        path.start = p0;
        path.end   = p1;
        straighten_path(path, triangles, positions, adjacencies);
        auto len = path_length(path, triangles, positions, adjacencies);
        solver.graph[i][k].length = len;
#endif
      }
    }
  }
  return solver;
}

// Builds a graph-based geodesic solver with arcs arranged in counterclockwise
// order by using the vertex-to-face adjacencies
geodesic_solver make_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t) {
  auto solver = geodesic_solver{};
  solver.graph.resize(positions.size());
  for (auto i = 0; i < positions.size(); ++i) {
    auto& star = v2t[i];
    auto& vert = positions[i];
    if (star.size() == 0) continue;
    for (auto j = 0; j < star.size(); ++j) {
      auto tid    = star[j];
      auto offset = find_in_vec(triangles[tid], i);
      auto p      = triangles[tid][(offset + 1) % 3];
      auto e      = positions[p] - vert;
      solver.graph[i].push_back({p, length(e)});
      auto opp   = opposite_face(triangles, adjacencies, tid, i);
      auto strip = vector<int>{tid, opp};
      auto k     = find_in_vec(
          adjacencies[tid], opp);  // TODO(fabio): this is not needeed
      assert(k != -1);
      auto a       = opposite_vertex(triangles, adjacencies, tid, k);
      offset       = find_in_vec(triangles[opp], a);
      auto bary    = zero3f;
      bary[offset] = 1;
      auto l       = length_by_flattening(
          triangles, positions, adjacencies, {opp, {bary.y, bary.z}}, strip);
      solver.graph[i].push_back({a, l});
    }
  }
  return solver;
}

// TODO: cleanup
vector<vector<float>> compute_angles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, vector<float>& total_angles,
    bool with_opposite) {
  auto angles = vector<vector<float>>(positions.size());
  total_angles.resize(positions.size());
  for (auto i = 0; i < positions.size(); ++i) {
    auto& star = v2t[i];
    if (star.size() == 0) continue;
    auto theta       = 0.0f;
    auto curr_angles = vector<float>{};
    auto tid         = 0;
    auto tr2d        = unfold_triangle{};
    auto tr3d        = zero3i;
    for (auto j = 0; j < star.size(); ++j) {
      tid  = star[j];
      tr3d = triangles[tid];
      if (with_opposite) {
        auto opp    = opposite_face(triangles, adjacencies, tid, i);
        auto k      = find_in_vec(adjacencies[tid], opp);
        tr2d        = init_flat_triangle(positions, tr3d);
        auto tr_opp = unfold_face(
            triangles, positions, adjacencies, tr2d, tid, k);
        auto h        = find_in_vec(adjacencies[opp], tid);
        auto flat_opp = tr_opp[(h + 2) % 3];
        auto w0       = tr2d[k] - tr2d[(k + 2) % 3];
        auto w1       = tr2d[(k + 1) % 3] - tr2d[(k + 2) % 3];
        auto v        = flat_opp - tr2d[(k + 2) % 3];
        auto phi = angle(w0, w1), theta0 = angle(v, w0), theta1 = angle(v, w1);

        if (theta0 < phi && theta1 < phi) {
          // no concave
          curr_angles.push_back(theta0);
          curr_angles.push_back(theta1);
          theta += phi;
        } else if (theta0 > theta1) {
          curr_angles.push_back(phi);
          curr_angles.push_back(phi);
          theta += phi;
        } else {
          curr_angles.push_back(0);
          curr_angles.push_back(phi);
          theta += phi;
        }
      } else {
        auto k = find_in_vec(triangles[tid], i);
        auto v = positions[tr3d[(k + 1) % 3]] - positions[tr3d[k]];
        auto w = positions[tr3d[(k + 2) % 3]] - positions[tr3d[k]];
        theta += angle(v, w);
      }
    }

    auto scale_factor = 2 * pif / theta;
    total_angles[i]   = theta;
    angles[i].push_back(0);
    theta = 0;
    for (auto j = 0; j < curr_angles.size() - 1; ++j) {
      auto curr_angle = curr_angles[j] * scale_factor;
      theta += curr_angle;
      angles[i].push_back(theta);
    }
  }
  return angles;
}

// `update` is a function that is executed during expansion, every time a node
// is put into queue. `exit` is a function that tells whether to expand the
// current node or perform early exit.
template <typename Update, typename Stop, typename Exit>
void visit_geodesic_graph(vector<float>& field, const geodesic_solver& solver,
    const vector<int>& sources, Update&& update, Stop&& stop, Exit&& exit) {
  /*
     This algortithm uses the heuristic Small Label Fisrt and Large Label Last
     https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

     Large Label Last (LLL): When extracting nodes from the queue, pick the
     front one. If it weights more than the average weight of the queue, put
     on the back and check the next node. Continue this way.
     Sometimes average_weight is less than every value due to floating point
     errors (doesn't happen with double precision).

     Small Label First (SLF): When adding a new node to queue, instead of
     always pushing it to the end of the queue, if it weights less than the
     front node of the queue, it is put on front. Otherwise the node is put at
     the end of the queue.
  */

  auto in_queue = vector<bool>(solver.graph.size(), false);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  auto cumulative_weight = 0.0;

  // setup queue
  auto queue = deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node)) break;
    if (stop(node)) continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      // Distance of neighbor through this node
      auto new_distance = field[node] + solver.graph[node][i].length;
      auto neighbor     = solver.graph[node][i].node;

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      update(node, neighbor, new_distance);
    }
  }
}

// Compute geodesic distances
void update_geodesic_distances(vector<float>& distances,
    const geodesic_solver& solver, const vector<int>& sources,
    float max_distance) {
  auto update = [](int node, int neighbor, float new_distance) {};
  auto stop   = [&](int node) { return distances[node] > max_distance; };
  auto exit   = [](int node) { return false; };
  visit_geodesic_graph(distances, solver, sources, update, stop, exit);
}

vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<int>& sources, float max_distance) {
  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto source : sources) distances[source] = 0.0f;
  update_geodesic_distances(distances, solver, sources, max_distance);
  return distances;
}

// Compute all shortest paths from source vertices to any other vertex.
// Paths are implicitly represented: each node is assigned its previous node
// in the path. Graph search early exits when reching end_vertex.
vector<int> compute_geodesic_parents(
    const geodesic_solver& solver, const vector<int>& sources, int end_vertex) {
  auto parents   = vector<int>(solver.graph.size(), -1);
  auto distances = vector<float>(solver.graph.size(), flt_max);
  auto update    = [&parents](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
  };
  auto stop = [end_vertex](int node) { return node == end_vertex; };
  auto exit = [](int node) { return false; };
  for (auto source : sources) distances[source] = 0.0f;
  visit_geodesic_graph(distances, solver, sources, update, stop, exit);
  return parents;
}

// Sample vertices with a Poisson distribution using geodesic distances
// Sampling strategy is farthest point sampling (FPS): at every step
// take the farthers point from current sampled set until done.
vector<int> sample_vertices_poisson(
    const geodesic_solver& solver, int num_samples) {
  auto verts = vector<int>{};
  verts.reserve(num_samples);
  auto distances = vector<float>(solver.graph.size(), flt_max);
  while (true) {
    auto max_index =
        (int)(std::max_element(distances.begin(), distances.end()) -
              distances.begin());
    verts.push_back(max_index);
    if (verts.size() >= num_samples) break;
    distances[max_index] = 0;
    update_geodesic_distances(distances, solver, {max_index}, flt_max);
  }
  return verts;
}

// Compute the distance field needed to compute a voronoi diagram
vector<vector<float>> compute_voronoi_fields(
    const geodesic_solver& solver, const vector<int>& generators) {
  auto fields = vector<vector<float>>(generators.size());

  // Find max distance from a generator to set an early exit condition for the
  // following distance field computations. This optimization makes
  // computation time weakly dependant on the number of generators.
  auto total = compute_geodesic_distances(solver, generators);
  auto max   = *std::max_element(total.begin(), total.end());
  for (auto i = 0; i < generators.size(); ++i) {
    fields[i]                = vector<float>(solver.graph.size(), flt_max);
    fields[i][generators[i]] = 0;
    fields[i] = compute_geodesic_distances(solver, {generators[i]}, max);
  }
  return fields;
}

vector<vec3f> colors_from_field(
    const vector<float>& field, float scale, const vec3f& c0, const vec3f& c1) {
  auto colors = vector<vec3f>{field.size()};
  for (auto i = 0; i < colors.size(); i++) {
    colors[i] = ((int64_t)(field[i] * scale)) % 2 ? c0 : c1;
  }
  return colors;
}

// `update` is a function that is executed during expansion, every time a node
// is put into queue. `exit` is a function that tells whether to expand the
// current node or perform early exit.
// DO NOT TOUCH THIS FUNCTION!!!
template <typename Update, typename Stop, typename Exit>
void visit_geodesic_graph(vector<float>& field,
    const dual_geodesic_solver& solver, const vector<int>& sources,
    Update&& update, Stop&& stop, Exit&& exit) {
  // DO NOT TOUCH THIS FUNCTION!!!
  /*
     This algortithm uses the heuristic Small Label Fisrt and Large Label Last
     https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm

     Large Label Last (LLL): When extracting nodes from the queue, pick the
     front one. If it weights more than the average weight of the queue, put
     on the back and check the next node. Continue this way.
     Sometimes average_weight is less than every value due to floating point
     errors (doesn't happen with double precision).

     Small Label First (SLF): When adding a new node to queue, instead of
     always pushing it to the end of the queue, if it weights less than the
     front node of the queue, it is put on front. Otherwise the node is put at
     the end of the queue.
  */

  auto in_queue = vector<bool>(solver.graph.size(), false);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  auto cumulative_weight = 0.0;

  // setup queue
  auto queue = std::deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node)) break;
    if (stop(node)) continue;

    for (auto i = 0; i < solver.graph[node].size(); i++) {
      // Distance of neighbor through this node
      auto new_distance = field[node] + solver.graph[node][i].length;
      auto neighbor     = solver.graph[node][i].node;

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      if (update(node, neighbor, new_distance)) return;
    }
  }
}

vector<mesh_point> compute_shortest_path(const dual_geodesic_solver& graph,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& start,
    const mesh_point& end) {
  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end   = end;
    path.strip = {start.face};
  } else {
    auto strip = strip_on_dual_graph(
        graph, triangles, positions, end.face, start.face);
    path = shortest_path(triangles, positions, adjacencies, start, end, strip);
  }
  // get mesh points
  return convert_mesh_path(
      triangles, adjacencies, path.strip, path.lerps, path.start, path.end)
      .points;
};

vector<mesh_point> compute_shortest_path(const dual_geodesic_solver& graph,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& points) {
  // geodesic path
  auto path = vector<mesh_point>{};
  for (auto idx = 0; idx < (int)points.size() - 1; idx++) {
    auto segment = compute_shortest_path(
        graph, triangles, positions, adjacencies, points[idx], points[idx + 1]);
    path.insert(path.end(), segment.begin(), segment.end());
  }
  return path;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF INTEGRAL PATHS
// -----------------------------------------------------------------------------
namespace yocto {
vec3f compute_gradient(const vec3i& triangle, const vector<vec3f>& positions,
    const vector<float>& field) {
  auto xy     = positions[triangle.y] - positions[triangle.x];
  auto yz     = positions[triangle.z] - positions[triangle.y];
  auto zx     = positions[triangle.x] - positions[triangle.z];
  auto normal = normalize(cross(zx, xy));
  auto result = zero3f;
  result += field[triangle.x] * cross(normal, yz);
  result += field[triangle.y] * cross(normal, zx);
  result += field[triangle.z] * cross(normal, xy);
  return result;
}

// Description of a discrete path along the surface of a triangle
struct surface_path {
  struct vertex {
    vec2i edge  = {0, 0};
    int   face  = 0;
    float alpha = 0;
  };
  int            start, end;
  vector<vertex> vertices;
};

// TODO: cleanup
// The functionalities under namespace `integral_paths` are now obsolete and can
// be deleted.
namespace integral_paths {

static int adjacent_face(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int face, const vec2i& edge) {
  // Given a face and an edge, return the adjacent face
  for (auto k = 0; k < 3; ++k) {
    auto x = triangles[face][k];
    auto y = triangles[face][k < 2 ? k + 1 : 0];
    auto e = vec2i{x, y};
    if (e == edge) return adjacency[face][k];
    if (vec2i{e.y, e.x} == edge) return adjacency[face][k];
  }
  return -1;
}

static vector<int> get_face_ring(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int face, int vertex) {
  auto contains = [](const vec3i& v, int a) -> bool {
    return v.x == a || v.y == a || v.z == a;
  };
  auto containsv = [](const vector<int>& v, int a) -> bool {
    return find(v.begin(), v.end(), a) != v.end();
  };

  auto result = vector<int>();
  result.reserve(12);

  if (face == -1) {
    for (int i = 0; i < triangles.size(); ++i) {
      if (contains(triangles[i], vertex)) {
        face = i;
        break;
      }
    }
  }

  auto queue = vector<int>();
  queue.push_back(face);

  while (!queue.empty()) {
    auto f = queue.back();
    queue.pop_back();
    result.push_back(f);

    for (auto k = 0; k < 3; ++k) {
      auto edge = vec2i{triangles[f][k], triangles[f][mod3(k + 1)]};
      if (edge.x == vertex || edge.y == vertex) {
        int neighbor_face = adjacency[f][k];
        if (neighbor_face == -1) continue;
        if (!containsv(result, neighbor_face)) {
          queue.push_back(neighbor_face);
          result.push_back(neighbor_face);
        }
      }
    }
  }

  return result;
}

static float step_from_point_to_edge(const vec3f& right, const vec3f& left,
    const vec3f& direction, float epsilon = 0.0001f) {
  auto normal            = cross(right, left);
  auto inverse_transform = mat3f{right - left, left, normal};
  auto transform         = inverse(inverse_transform);
  auto dir               = transform * direction;
  return clamp(1 - dir.x / dir.y, 0.0f + epsilon, 1.0f - epsilon);
}

static std::pair<float, bool> step_from_edge_to_edge(const vec3f& point,
    const vec3f& a, const vec3f& b, const vec3f& c, const vec3f& direction,
    float epsilon = 0.0001f) {
  //      b
  //     /\
  //    /  \
  //   /    \
  //  /______\
  // c        a

  auto right = a - point, left = c - point, front = b - point;
  auto normal = triangle_normal(a, b, c);

  auto right_side = dot(cross(direction, front), normal) > 0;
  if (right_side) {
    auto below_edge = dot(cross(direction, right), normal) > 0;
    if (below_edge) return {epsilon, true};
    auto x = step_from_point_to_edge(right, front, direction);
    return {x, true};
  } else {
    auto below_edge = dot(cross(left, direction), normal) > 0;
    if (below_edge) return {1.0f - epsilon, false};
    auto x = step_from_point_to_edge(front, left, direction);
    return {x, false};
  }
}

static surface_path::vertex step_from_point(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, const vector<float>& field, int vertex,
    int start_face, int tag = -1, float epsilon = 0.0001f) {
  auto opposite_edge = [](int vertex, const vec3i& tr) -> vec2i {
    for (auto i = 0; i < 3; ++i) {
      if (tr[i] == vertex) return {tr[mod3(i + 1)], tr[mod3(i + 2)]};
    }
    return {-1, -1};
  };
  auto is_direction_inbetween = [](const vec3f& right, const vec3f& left,
                                    const vec3f& direction,
                                    float        threshold = 0) -> bool {
    auto normal      = cross(right, left);
    auto cross_right = cross(right, direction);
    auto cross_left  = cross(direction, left);
    return dot(cross_right, normal) > threshold &&
           dot(cross_left, normal) > threshold;
  };

  auto triangle_fan = get_face_ring(triangles, adjacency, start_face, vertex);

  auto best_alignment = 0.0;
  auto fallback_lerp  = surface_path::vertex{vec2i{-1, -1}, -1, 0};

  for (auto i = 0; i < triangle_fan.size(); ++i) {
    auto face = triangle_fan[i];
    if (tag != -1 && tags[face] != tag) continue;
    auto edge = opposite_edge(vertex, triangles[face]);
    if (edge == vec2i{-1, -1}) throw std::runtime_error("edge is {-1, -1}\n");

    auto a = positions[vertex], b = positions[edge.x], c = positions[edge.y];
    auto left = c - a, right = b - a;
    auto direction = normalize(
        compute_gradient(triangles[face], positions, field));

    auto right_dot = dot(normalize(right), direction);
    auto left_dot  = dot(normalize(left), direction);

    // Look for the most aligned edge with the gradient
    // direction. If no other face is suitable for selection, we return the
    // most aligned edge as result (fallback_lerp)
    if (max(right_dot, left_dot) > best_alignment) {
      best_alignment = max(right_dot, left_dot);
      auto alpha     = right_dot > left_dot ? 0.0f + epsilon : 1.0f - epsilon;
      fallback_lerp  = surface_path::vertex{edge, face, alpha};
    }

    // Check if gradient direction is in between ac and ab.
    if (is_direction_inbetween(right, left, direction)) {
      auto alpha = step_from_point_to_edge(right, left, direction);
      return surface_path::vertex{edge, face, alpha};
    }
  }

  return fallback_lerp;
}

surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from) {
  auto opposite_vertex = [](const vec3i& tr, const vec2i& edge) -> int {
    for (auto i = 0; i < 3; ++i) {
      if (tr[i] != edge.x && tr[i] != edge.y) return tr[i];
    }
    return -1;
  };
  auto find_index = [](const vec3i& v, int x) -> int {
    if (v.x == x) return 0;
    if (v.y == x) return 1;
    if (v.z == x) return 2;
    return -1;
  };

  // trace function
  auto lerps = vector<surface_path::vertex>();
  auto lerp  = step_from_point(
      triangles, positions, adjacency, tags, field, from, -1, tag);
  if (lerp.face == -1) return {};
  lerps.push_back(lerp);

  const int num_steps = 10000;

  for (auto i = 0; i < num_steps; i++) {
    auto [old_edge, old_face, old_alpha] = lerps.back();
    if (old_face == -1) throw std::runtime_error("programmer error");
    auto point = (1 - old_alpha) * positions[old_edge.x] +
                 old_alpha * positions[old_edge.y];

    auto face = adjacent_face(triangles, adjacency, old_face, old_edge);
    if (face == -1) {
      lerps.push_back({vec2i{-1, -1}, face, 0});
      return surface_path{from, -1, lerps};
    }

    if (tags[face] != tag) {
      auto k  = find_index(triangles[face], old_edge.x);
      auto to = triangles[face][mod3(k + 1)];

      // @Hack!: We store the tag of the reached region in edge.x
      auto edge = vec2i{to, tags[face]};
      lerps.push_back({edge, face, 0});
      return surface_path{from, to, lerps};
    }

    auto direction = normalize(
        compute_gradient(triangles[face], positions, field));

    auto front_idx = opposite_vertex(triangles[face], old_edge);
    if (front_idx == -1) {
      throw std::runtime_error("programmer error");
      return {};
    }

    if (old_alpha < 0 || old_alpha > 1)
      throw std::runtime_error("programmer error");

    auto& a = positions[old_edge.x];
    auto& b = positions[front_idx];
    auto& c = positions[old_edge.y];

    auto [x, step_right] = step_from_edge_to_edge(point, a, b, c, direction);

    auto edge = old_edge;
    if (step_right) {
      point  = (1 - x) * a + x * b;
      edge.y = front_idx;
    } else {
      point  = (1 - x) * b + x * c;
      edge.x = front_idx;
    }
    lerps.push_back({edge, face, x});
  }

  throw std::runtime_error("integral path ended nowhere");
  return surface_path{from, 0, lerps};
}

surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from,
    int to) {
  auto opposite_vertex = [](const vec3i& tr, const vec2i& edge) -> int {
    for (int i = 0; i < 3; ++i) {
      if (tr[i] != edge.x && tr[i] != edge.y) return tr[i];
    }
    return -1;
  };
  auto contains = [](const vec3i& v, int a) -> bool {
    return v.x == a || v.y == a || v.z == a;
  };

  auto lerps = vector<surface_path::vertex>();

  lerps.push_back(
      step_from_point(triangles, positions, adjacency, tags, field, from, -1));

  const int num_steps = 10000;

  for (auto i = 0; i < num_steps; i++) {
    auto [old_edge, old_face, old_alpha] = lerps.back();
    auto point = (1 - old_alpha) * positions[old_edge.x] +
                 old_alpha * positions[old_edge.y];

    auto face = adjacent_face(triangles, adjacency, old_face, old_edge);
    if (face == -1) break;

    if (contains(triangles[face], to)) {
      for (int k = 0; k < 3; ++k) {
        auto edge = vec2i{triangles[face][k], triangles[face][mod3(k + 1)]};
        if (edge.x == to) {
          lerps.push_back({edge, face, 0});
          return {from, to, lerps};
        }
      }
    }
    auto direction = normalize(
        compute_gradient(triangles[face], positions, field));

    int front_idx = opposite_vertex(triangles[face], old_edge);
    if (front_idx == -1) {
      throw std::runtime_error("programmer error: front_idx is -1");
      break;
    }

    if (old_alpha < 0 || old_alpha > 1)
      throw std::runtime_error("programmer error");
    if (old_alpha == 0 || old_alpha == 1) {
      int  vertex = old_alpha == 0 ? old_edge.x : old_edge.y;
      auto lerp   = step_from_point(
          triangles, positions, adjacency, tags, field, vertex, old_face);
      lerps.push_back(lerp);
      if (lerp.alpha == 0 && lerp.edge.x == to) break;
      if (lerp.alpha == 1 && lerp.edge.y == to) break;
      continue;
    }

    auto& a = positions[old_edge.x];
    auto& b = positions[front_idx];
    auto& c = positions[old_edge.y];

    auto [x, step_right] = step_from_edge_to_edge(point, a, b, c, direction);

    auto edge = old_edge;
    if (step_right) {
      point  = (1 - x) * a + x * b;
      edge.y = front_idx;
    } else {
      point  = (1 - x) * b + x * c;
      edge.x = front_idx;
    }
    if (opposite_vertex(triangles[face], edge) == -1) {
      throw std::runtime_error("opposite vertex == -1");
    }

    lerps.push_back({edge, face, x});
    if (x == 0 && edge.x == to) break;
    if (x == 1 && edge.y == to) break;
  }

  return {from, to, lerps};
}

vector<vec3f> make_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions) {
  if (path.vertices.empty()) return {};

  auto positions = vector<vec3f>();
  positions.reserve(path.vertices.size() + 1);
  positions.push_back(mesh_positions[path.start]);

  for (auto i = 0; i < path.vertices.size() - 1; ++i) {
    auto [edge, face, x] = path.vertices[i];
    auto p0              = mesh_positions[edge.x];
    auto p1              = mesh_positions[edge.y];
    auto position        = (1 - x) * p0 + x * p1;
    positions.push_back(position);
  }

  if (path.end != -1) {
    positions.push_back(mesh_positions[path.end]);
  }
  return positions;
}

}  // namespace integral_paths

// TODO: cleanup
// Trace integral path following the gradient of a scalar field
surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from) {
  return integral_paths::integrate_field(
      triangles, positions, adjacency, tags, tag, field, from);
}

surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from,
    int to) {
  return integral_paths::integrate_field(
      triangles, positions, adjacency, tags, tag, field, from, to);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto {

vector<vec3f> make_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions) {
  return integral_paths::make_positions_from_path(path, mesh_positions);
}

static vector<pair<int, float>> nodes_around_point(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& p) {
  auto nodes               = vector<pair<int, float>>{};
  auto [is_vertex, offset] = bary_is_vert(get_bary(p.uv));
  if (is_vertex) {
    auto vid = triangles[p.face][offset];
    nodes.push_back({vid, 0});
  } else {
    auto pid = p.face;
    auto pos = eval_position(triangles, positions, p);
    for (int i = 0; i < 3; ++i) {
      int   p0 = triangles[pid][i], p1 = triangles[pid][(i + 1) % 3];
      float d = length(positions[p0] - pos);
      nodes.push_back(std::make_pair(p0, d));
      int         CW_pid = adjacencies[pid][i];
      int         opp    = opposite_vertex(triangles, adjacencies, pid, i);
      vector<int> strip  = {CW_pid, pid};
      float       l      = length_by_flattening(
          triangles, positions, adjacencies, p, strip);
      nodes.push_back(std::make_pair(opp, l));
      int opp_pid = opposite_face(triangles, adjacencies, CW_pid, p0);
      strip       = {opp_pid, CW_pid, pid};
      int k       = find_in_vec(adjacencies[CW_pid], opp_pid);
      int q       = opposite_vertex(triangles, adjacencies, CW_pid, k);
      d = length_by_flattening(triangles, positions, adjacencies, p, strip);
      nodes.push_back(std::make_pair(q, d));
      opp_pid = opposite_face(triangles, adjacencies, CW_pid, p1);
      strip   = {opp_pid, CW_pid, pid};
      k       = find_in_vec(adjacencies[CW_pid], opp_pid);
      q       = opposite_vertex(triangles, adjacencies, CW_pid, k);
      d = length_by_flattening(triangles, positions, adjacencies, p, strip);
      nodes.push_back(std::make_pair(q, d));
    }
  }
  return nodes;
}

#if 0
// compute the distance between a point p and some vertices around him
// handling concave path
static vector<pair<int, float>> nodes_around_point(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& p) {
  auto nodes = vector<pair<int, float>>{};
  if (auto [is_vert, offset] = point_is_vert(p); is_vert) {
    auto vid = triangles[p.face][offset];
    nodes.push_back({vid, 0});
  } else {
    auto pid = p.face;
    auto pos = eval_position(triangles, positions, p);
    for (auto i = 0; i < 3; ++i) {
      auto p0 = triangles[pid][i];
      // auto p1 = triangles[pid][(i + 1) % 3];
      auto d = length(positions[p0] - pos);
      nodes.push_back({p0, d});
      auto cw_pid = adjacencies[pid][i];
      auto opp    = opposite_vertex(triangles, adjacencies, pid, i);
      auto strip  = vector<int>{cw_pid, pid};
      auto l      = length_by_flattening(
          triangles, positions, adjacencies, p, strip);
      nodes.push_back({opp, l});
    }
  }

  return nodes;
}
#endif

vector<float> solve_with_parents(const geodesic_solver& solver,
    const vector<pair<int, float>>&                     sources_and_dist,
    const vector<pair<int, float>>& targets, vector<int>& parents,
    bool with_parents = false) {
  parents.assign(solver.graph.size(), -1);
  auto update = [&parents](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
  };
  auto stop = [](int node) { return false; };

  auto exit_verts = vector<int>(targets.size());
  for (auto i = 0; i < targets.size(); ++i) {
    exit_verts[i] = targets[i].first;
  }
  auto exit = [&exit_verts](int node) {
    auto it = find(exit_verts.begin(), exit_verts.end(), node);
    if (it != exit_verts.end()) {
      exit_verts.erase(it);
    }
    return exit_verts.empty();
  };

  auto distances  = vector<float>(solver.graph.size(), flt_max);
  auto sources_id = vector<int>(sources_and_dist.size());
  for (auto i = 0; i < sources_and_dist.size(); ++i) {
    sources_id[i]                        = sources_and_dist[i].first;
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }

  visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);
  return distances;
}

// given a set of vertices and distances (nbr) computed with
// "nodes_around_point" and a scalar field (f), returns the parent of the
// point having nbr as neighborhood
int set_target_parent(const vector<pair<int, float>>& nbr,
    const vector<int>& parents, const vector<float>& f) {
  if (nbr.size() == 1) {
    return parents[nbr[0].first];
  } else {
    auto vid    = -1;
    auto lambda = flt_max;
    for (auto i = 0; i < nbr.size(); ++i) {
      auto val = f[nbr[i].first] + nbr[i].second;
      if (val < lambda) {
        lambda = val;
        vid    = nbr[i].first;
      }
    }
    return vid;
  }
}

// given a vector of parents(parents) and a starting vertex (target_parent),
// return the vertex v such that parents[v]=-1; all the vertices visited
// during the navigation are stored in "path"
int set_source_child(
    const vector<int>& parents, int target_parent, vector<int>& path) {
  auto stop          = false;
  auto prev          = target_parent;
  auto source_parent = parents[prev];
  path               = {target_parent};
  if (source_parent == -1) return target_parent;
  while (!stop) {
    prev = source_parent;
    path.push_back(prev);
    source_parent = parents[source_parent];
    if (source_parent == -1) stop = true;
  }
  return prev;
}

// utilities:nodes_around_point-->length_by_flattening
// returns the shortest path starting from target to source as a list of
// indices of nodes of the graph(solver). note: the list does not contains
// source if it is a vertex
vector<int> compute_geodesic_parents(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& source,
    const mesh_point& target) {
  auto source_nodes = nodes_around_point(
      triangles, positions, adjacencies, source);
  auto target_nodes = nodes_around_point(
      triangles, positions, adjacencies, target);
  auto parents   = vector<int>{};
  auto distances = solve_with_parents(
      solver, source_nodes, target_nodes, parents);
  auto target_parent = -1, source_child = -1;
  if (target_nodes.size() == 1)
    target_parent = parents[target_nodes[0].first];
  else
    target_parent = set_target_parent(target_nodes, parents, distances);
  auto path = vector<int>{};
  if (target_parent == -1) return path;
  source_child = set_source_child(parents, target_parent, path);
  if (source_nodes.size() == 1)
    path.pop_back();  // we remove "source" from the list of parents if it is
                      // a vertex
  return path;
}

// same function of "compute_geodesic_parents" of yocto_mesh.cpp that makes
// early exit when reaching end_vertex, the name is changed because the input
// parameters are the same
vector<int> compute_pruned_geodesic_parents(
    const geodesic_solver& solver, const vector<int>& sources, int end_vertex) {
  auto parents   = vector<int>(solver.graph.size(), -1);
  auto distances = vector<float>(solver.graph.size(), flt_max);
  auto update    = [&parents](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
  };
  auto stop = [](int node) { return false; };
  auto exit = [end_vertex](int node) { return node == end_vertex; };
  for (auto source : sources) distances[source] = 0.0f;
  visit_geodesic_graph(distances, solver, sources, update, stop, exit);
  return parents;
}

// returns the shortest path starting from target to source as a list of
// indices of nodes of the graph(solver). note: the list does not contains
// source and target.
vector<int> compute_geodesic_parents(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int source, int target) {
  vector<int> sources = {source};
  vector<int> parents = compute_pruned_geodesic_parents(
      solver, sources, target);
  vector<int> path = {};
  set_source_child(parents, target, path);
  path.pop_back();  // we remove source from the list
  return path;
}

static vector<pair<int, float>> check_nodes(vector<pair<int, float>>& nodes) {
  sort(nodes.begin(), nodes.end());
  auto new_nodes = vector<pair<int, float>>{};
  for (auto i = 1; i < nodes.size(); ++i) {
    auto prev = nodes[i - 1];
    auto curr = nodes[i];
    if (prev.first == curr.first) {
      auto d0 = prev.second, d1 = curr.second;
      if (d0 <= d1)
        new_nodes.push_back(prev);
      else
        new_nodes.push_back(curr);
      ++i;
    } else {
      new_nodes.push_back(prev);
    }
  }
  nodes = new_nodes;
  return nodes;
}

static vector<float> solve(const geodesic_solver& solver,
    const vector<pair<int, float>>&               sources_and_dist) {
  auto update = [](int node, int neighbor, float new_distance) {};
  auto stop   = [](int node) { return false; };
  auto exit   = [](int node) { return false; };

  auto distances = vector<float>{};
  distances.assign(solver.graph.size(), flt_max);
  auto sources_id = vector<int>(sources_and_dist.size());
  for (auto i = 0; i < sources_and_dist.size(); ++i) {
    sources_id[i]                        = sources_and_dist[i].first;
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }
  visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);

  return distances;
}

// TODO: cleanup
#if 0
static vector<float> solve_with_parents(const geodesic_solver& solver,
    const vector<pair<int, float>>&                     sources_and_dist,
    const vector<pair<int, float>>& targets, vector<int>& parents,
    bool with_parents = false) {
  parents.assign(solver.graph.size(), -1);
  auto update = [&parents](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
  };
  auto stop = [](int node) { return false; };

  vector<int> exit_verts(targets.size());
  for (int i = 0; i < targets.size(); ++i) {
    exit_verts[i] = targets[i].first;
  }
  auto exit = [&exit_verts](int node) {
    auto it = find(exit_verts.begin(), exit_verts.end(), node);
    if (it != exit_verts.end()) {
      exit_verts.erase(it);
    }
    return exit_verts.empty();
  };
  vector<float> distances;
  distances.assign(solver.graph.size(), flt_max);
  vector<int> sources_id(sources_and_dist.size());
  for (int i = 0; i < sources_and_dist.size(); ++i) {
    sources_id[i]                        = sources_and_dist[i].first;
    distances[sources_and_dist[i].first] = sources_and_dist[i].second;
  }

  visit_geodesic_graph(distances, solver, sources_id, update, stop, exit);

  return distances;
}
#endif

vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& sources) {
  auto source_nodes = vector<pair<int, float>>{};
  for (auto i = 0; i < sources.size(); ++i) {
    auto curr_nodes = nodes_around_point(
        triangles, positions, adjacencies, sources[i]);
    source_nodes.insert(
        source_nodes.end(), curr_nodes.begin(), curr_nodes.end());
  }
  if (source_nodes.size() > 1) check_nodes(source_nodes);
  return solve(solver, source_nodes);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STRIPS
// -----------------------------------------------------------------------------
namespace yocto {

// `update` is a function that is executed during expansion, every time a node
// is put into queue. `exit` is a function that tells whether to expand the
// current node or perform early exit.
// TODO(fabio): this needs a lof of cleaning
template <typename Update, typename Stop, typename Exit>
void heuristic_visit_geodesic_graph(vector<float>& field,
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int start, int end, Update&& update,
    Stop&& stop, Exit&& exit) {
  auto destination_pos = eval_position(
      triangles, positions, {end, {1.0f / 3, 1.0f / 3}});

  auto estimate_dist = [&](int face) {
    auto p = eval_position(triangles, positions, {face, {1.0f / 3, 1.0f / 3}});
    return length(p - destination_pos);
  };
  field[start] = estimate_dist(start);

  auto in_queue = vector<bool>(solver.graph.size(), false);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue      = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node)) break;
    if (stop(node)) continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1) continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      if (update(node, neighbor, new_distance)) return;
    }
  }
}

// Get strip by ascending geodesic distance field computed by the graph
vector<int> strip_ascending_distance_field(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int start, int end) {
  auto argmin = [](const vec3f& vec) -> int {
    auto min_index = 0;
    for (auto i = 1; i < 3; ++i)
      if (vec[i] < vec[min_index]) min_index = i;
    return min_index;
  };

  if (start == end) return {start};
  assert(start != -1);
  assert(end != -1);

  if (find_in_vec(adjacencies[start], end) != -1) {
    return {start, end};
  }

  auto [x, y, z] = triangles[end];
  auto tr        = triangles[start];
  auto distance  = vector<float>(positions.size(), flt_max);
  auto update    = [](int node, int neighbor, float new_distance) {};
  auto stop      = [](int node) { return false; };
  auto exit      = [&](int node) {
    return distance[tr.x] != flt_max && distance[tr.y] != flt_max &&
           distance[tr.z] != flt_max;
  };
  visit_geodesic_graph(distance, solver, {x, y, z}, update, stop, exit);

  auto result    = vector<int>{start};
  auto node      = start;
  auto visited   = vector<bool>(triangles.size(), false);
  visited[start] = true;

  // Descend geodesic field.
  while (node != end) {
    auto dist = vec3f{flt_max, flt_max, flt_max};
    for (auto i = 0; i < 3; ++i) {
      if (adjacencies[node][i] == -1) continue;
      if (visited[adjacencies[node][i]]) continue;
      auto v  = opposite_vertex(triangles, adjacencies, node, i);
      dist[i] = distance[v];
    }

    node = adjacencies[node][argmin(dist)];
    result.push_back(node);
    visited[node] = true;
  }

  // assert(check_strip(adjacencies, result));
  return result;
}

// TODO(fabio): this is copied somewhere else. Make sure it is th proper
// version.
#if 0
bool check_strip(const vector<vec3i>& adjacencies, const vector<int>& strip) {
  auto faces = std::unordered_set<int>{};
  faces.insert(strip[0]);
  for (auto i = 1; i < strip.size(); ++i) {
    assert(faces.count(strip[i]) == 0);  // face appears twice in the strip
    faces.insert(strip[i]);
    assert(find_in_vec(adjacencies[strip[i - 1]], strip[i]) != -1);
    assert(find_in_vec(adjacencies[strip[i]], strip[i - 1]) != -1);
  }
  return true;
}
#endif

vector<int> strip_on_dual_graph(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions, int start,
    int end) {
  if (start == end) return {start};

  // initialize once for all and sparsely cleanup at the end of every solve
  auto parents = vector<int>(solver.graph.size(), -1);
  auto field   = vector<float>(solver.graph.size(), flt_max);

  auto sources = vector<int>{start};
  auto update  = [&parents, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  heuristic_visit_geodesic_graph(
      field, solver, triangles, positions, start, end, update, stop, exit);
  // auto visited = sources;
  // visit_dual_graph(field, solver, triangles, positions, in_queue, start, end,
  //     update, visited);

  // extract_strip
  auto strip = vector<int>{};
  auto node  = end;
  assert(parents[end] != -1);
  strip.reserve((int)sqrt(parents.size()));
  // for (int i = 0; i < parents.size() && node != -1; i++) {
  while (node != -1) {
    assert(find_in_vec(strip, node) != 1);
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  // for (auto& v : visited) {
  //   parents[v]  = -1;
  //   field[v]    = flt_max;
  //   in_queue[v] = false;
  // }
  // assert(check_strip(adjacencies, strip));
  return strip;
}

// TODO: cleanup
static int node_is_neighboor(const geodesic_solver& solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}

// TODO: cleanup
static bool set_ord(int s, int prev_entry, int next_entry, bool nei_is_dual) {
  auto ccw_count = -1, cw_count = -1;
  if (prev_entry < next_entry) {
    ccw_count = next_entry - prev_entry;
    cw_count  = prev_entry + s - next_entry;
  } else {
    ccw_count = s - prev_entry + next_entry;
    cw_count  = prev_entry - next_entry;
  }
  if (!nei_is_dual) --ccw_count;
  if (ccw_count < cw_count)
    return true;
  else
    return false;
}

static void fill_strip(vector<int>& strip, const vector<vector<int>>& v2t,
    int vid, int first, int last, bool nei_is_dual, bool ccw) {
  auto  start = first, end = last;
  auto& star = v2t[vid];
  auto  s    = (int)star.size();
  if (ccw && !nei_is_dual)
    end = (s - 1 + end) % s;  // I can stop one face earlier;
  if (start == end) {
    if (strip.back() != star[start]) strip.push_back(star[start]);
  } else if (ccw) {
    if (strip.back() == star[start]) start = (start + 1) % s;
    if (start > end) end += s;
    for (auto i = start; i <= end; ++i) {
      strip.push_back(star[i % s]);
    }
  } else {
    if (strip.back() == star[start % s]) start = (s - 1 + start) % s;
    if (start < end) start += s;
    for (auto i = start; i >= end; --i) {
      strip.push_back(star[i % s]);
    }
  }
}

static int get_entry(vector<int>& strip, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, int parent, const mesh_point& p) {
  strip = {p.face};
  if (auto [is_vert, offset] = point_is_vert(p); is_vert) {
    auto vid   = triangles[p.face][offset];
    auto entry = node_is_neighboor(solver, vid, parent);
    if (entry < 0) {
      assert(vid == parent);
      return -1;
    }
    auto star        = v2t[vid];
    auto s           = (int)star.size();
    auto it          = find(star.begin(), star.end(), p.face);
    auto first       = (int)distance(star.begin(), it);
    auto last        = (entry % 2) ? (entry - 1) / 2 : (entry / 2) % s;
    auto ccw         = set_ord(s, first, last, entry % 2);
    auto nei_is_dual = (bool)((entry + 2) % 2);
    fill_strip(strip, v2t, vid, first, last, nei_is_dual, ccw);
    if (entry % 2) {
      first    = (entry - 1) / 2;
      auto tid = opposite_face(triangles, adjacencies, v2t[vid][first], vid);
      strip.push_back(tid);
    }
    entry = node_is_neighboor(solver, parent, vid);
    return entry;
  } else {
    auto        h = find_in_vec(triangles[p.face], parent);
    vector<int> adj_tri(3);
    if (h == -1) {
      for (auto i = 0; i < 3; ++i) {
        auto adj   = adjacencies[p.face][i];
        adj_tri[i] = adj;
        h          = find_in_vec(triangles[adj], parent);
        if (h != -1) {
          strip.push_back(adj);
          auto entry = node_is_neighboor(
              solver, parent, triangles[adj][(h + 1) % 3]);
          assert(entry >= 0);
          return entry + 1;
        }
      }
    } else {
      auto entry = node_is_neighboor(
          solver, parent, triangles[p.face][(h + 1) % 3]);
      assert(entry >= 0);
      return entry + 1;
    }
    for (auto i = 0; i < 3; ++i) {
      auto p0 = triangles[p.face][i], p1 = triangles[p.face][(i + 1) % 3];
      auto adj = adj_tri[i];
      auto opp = opposite_face(triangles, adjacencies, adj, p0);
      h        = find_in_vec(triangles[opp], parent);
      if (h != -1) {
        strip.push_back(adj);
        strip.push_back(opp);
        auto entry = node_is_neighboor(
            solver, parent, triangles[opp][(h + 1) % 3]);
        assert(entry >= 0);
        return entry + 1;
      }
      opp = opposite_face(triangles, adjacencies, adj, p1);
      h   = find_in_vec(triangles[opp], parent);
      if (h != -1) {
        strip.push_back(adj);
        strip.push_back(opp);
        auto entry = node_is_neighboor(
            solver, parent, triangles[opp][(h + 1) % 3]);
        assert(entry >= 0);
        return entry + 1;
      }
    }
  }
  assert(false);
  return -1;
}

#if 0
static int get_entry(vector<int>& strip, const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, int parent, const mesh_point& p) {
  strip = {p.face};
  if (auto [is_vert, offset] = point_is_vert(p); is_vert) {
    auto vid   = triangles[p.face][offset];
    auto entry = node_is_neighboor(solver, vid, parent);
    assert(entry >= 0);
    auto star        = v2t[vid];
    auto s           = (int)star.size();
    auto it          = find(star.begin(), star.end(), p.face);
    auto first       = (int)distance(star.begin(), it);
    auto last        = (entry % 2) ? (entry - 1) / 2 : (entry / 2) % s;
    auto ccw         = set_ord(s, first, last, entry % 2);
    auto nei_is_dual = (bool)(entry % 2);  // TODO(fabio): ma qusto e' giusto?
    fill_strip(strip, v2t, vid, first, last, nei_is_dual, ccw);
    if (entry % 2) {
      first    = (entry - 1) / 2;
      auto tid = opposite_face(triangles, adjacencies, v2t[vid][first], vid);
      strip.push_back(tid);
    }
    entry = node_is_neighboor(solver, parent, vid);
    return entry;
  } else {
    auto h = find_in_vec(triangles[p.face], parent);
    if (h == -1) {
      for (auto i = 0; i < 3; ++i) {
        auto adj = adjacencies[p.face][i];
        h        = find_in_vec(triangles[adj], parent);
        if (h != -1) {
          strip.push_back(adj);
          auto entry = node_is_neighboor(
              solver, parent, triangles[adj][(h + 1) % 3]);
          assert(entry >= 0);
          return entry + 1;
        }
      }
    } else {
      auto entry = node_is_neighboor(
          solver, parent, triangles[p.face][(h + 1) % 3]);
      assert(entry >= 0);
      return entry + 1;
    }
  }
  return 0;  // TODO(fabio): cosa deve fare qui?
}
#endif

void close_strip(vector<int>& strip, const vector<vector<int>>& v2t, int vid,
    int prev_tri, int last_tri) {
  auto star    = v2t[vid];
  auto s       = (int)star.size();
  auto prev_it = find(star.begin(), star.end(), prev_tri);
  assert(prev_it != star.end());
  auto first   = (int)distance(star.begin(), prev_it);
  auto next_it = find(star.begin(), star.end(), last_tri);
  assert(next_it != star.end());
  auto last = (int)distance(star.begin(), next_it);
  auto ccw  = set_ord(s, first, last, true);
  fill_strip(strip, v2t, vid, first, last, true, ccw);
}

// particular case of "get strip" when one of the two point is the parent of
// the other so the size of the strip is one or two
static vector<int> short_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  auto entry = 0;
  auto strip = vector<int>{};
  if (auto [is_vert, offset] = point_is_vert(target); is_vert) {
    auto vid = triangles[target.face][offset];
    entry    = get_entry(strip, solver, triangles, positions, adjacencies, v2t,
        angles, vid, source);
    if (strip.back() != target.face) {
      close_strip(strip, v2t, vid, strip.back(), target.face);
    }
    reverse(strip.begin(), strip.end());
  } else if (auto [is_vert, offset] = point_is_vert(source); is_vert) {
    auto vid = triangles[source.face][offset];
    entry    = get_entry(strip, solver, triangles, positions, adjacencies, v2t,
        angles, vid, target);
    if (strip.back() != source.face) {
      close_strip(strip, v2t, vid, strip.back(), source.face);
    }
  } else {
    assert(false);
  }
  return strip;
}

// returns a strip of triangles such target belongs to the first one and
// source to the last one
// TODO(fabio): better name
vector<int> get_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  if (target.face == source.face) return {target.face};
  auto parents = compute_geodesic_parents(
      solver, triangles, positions, adjacencies, source, target);
  auto N     = (int)parents.size();
  auto first = 0, last = 0, prev_entry = 0, next_entry = 0;
  auto strip_to_point = vector<int>{}, strip = vector<int>{};
  auto ccw = false, nei_is_dual = false;
  if (N == 0) {
    return short_strip(
        solver, triangles, positions, adjacencies, v2t, angles, source, target);
  } else if (N == 1) {
    auto v      = parents[0];
    prev_entry  = get_entry(strip, solver, triangles, positions, adjacencies,
        v2t, angles, v, target);
    next_entry  = get_entry(strip_to_point, solver, triangles, positions,
        adjacencies, v2t, angles, v, source);
    first       = find_in_vec(v2t[v], strip.back());
    nei_is_dual = next_entry % 2;
    last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
    ccw         = set_ord((int)v2t[v].size(), first, last, nei_is_dual);
    fill_strip(strip, v2t, v, first, last, nei_is_dual, ccw);
    close_strip(strip, v2t, v, strip.back(), strip_to_point.back());
    if (strip.back() == strip_to_point.back()) strip_to_point.pop_back();
    reverse(strip_to_point.begin(), strip_to_point.end());
    strip.insert(strip.end(), strip_to_point.begin(), strip_to_point.end());
    return strip;
  } else {
    prev_entry = get_entry(strip, solver, triangles, positions, adjacencies,
        v2t, angles, parents[0], target);
  }

  for (auto i = 0; i < N; ++i) {
    auto v = parents[i];
    if (i == N - 1) {
      first      = find_in_vec(v2t[v], strip.back());
      next_entry = get_entry(strip_to_point, solver, triangles, positions,
          adjacencies, v2t, angles, v, source);
      last       = find_in_vec(v2t[v], strip_to_point.back());
      ccw        = set_ord((int)v2t[v].size(), first, last, next_entry % 2);
    } else {
      first = find_in_vec(v2t[v], strip.back());
      assert(first != -1);
      next_entry  = node_is_neighboor(solver, v, parents[i + 1]);
      nei_is_dual = next_entry % 2;
      last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
      ccw         = set_ord((int)v2t[v].size(), first, last, nei_is_dual);
    }

    fill_strip(strip, v2t, v, first, last, nei_is_dual, ccw);

    if (nei_is_dual && i != N - 1) {
      auto tid = opposite_face(triangles, adjacencies, strip.back(), v);
      strip.push_back(tid);
    }
  }

  close_strip(strip, v2t, parents.back(), strip.back(), strip_to_point.back());
  if (strip.back() == strip_to_point.back()) strip_to_point.pop_back();
  reverse(strip_to_point.begin(), strip_to_point.end());
  strip.insert(strip.end(), strip_to_point.begin(), strip_to_point.end());
  return strip;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEODESIC PATH
// -----------------------------------------------------------------------------
namespace yocto {

vector<vec3f> path_positions(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies) {
  auto result = vector<vec3f>(path.lerps.size() + 2);
  result[0]   = eval_position(triangles, positions, path.start);
  for (auto i = 0; i < path.lerps.size(); i++) {
    auto e = get_edge(
        triangles, positions, adjacencies, path.strip[i], path.strip[i + 1]);
    if (e == vec2i{-1, -1}) continue;
    auto x        = path.lerps[i];
    auto p        = lerp(positions[e.x], positions[e.y], x);
    result[i + 1] = p;
  }
  result.back() = eval_position(triangles, positions, path.end);
  return result;
}

vector<float> path_parameters(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies) {
  return path_parameters(
      path_positions(path, triangles, positions, adjacencies));
}

float path_length(const geodesic_path& path, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies) {
  return path_length(path_positions(path, triangles, positions, adjacencies));
}

vector<float> path_parameters(const vector<vec3f>& positions) {
  auto len         = 0.0f;
  auto parameter_t = vector<float>(positions.size());
  for (auto i = 0; i < positions.size(); i++) {
    if (i) len += length(positions[i] - positions[i - 1]);
    parameter_t[i] = len;
  }
  for (auto& t : parameter_t) t /= len;
  return parameter_t;
}

float path_length(const vector<vec3f>& positions) {
  auto len = 0.0f;
  for (auto i = 0; i < positions.size(); i++) {
    if (i) len += length(positions[i] - positions[i - 1]);
  }
  return len;
}

// Find barycentric coordinates of a point inside a triangle (a, b, c).
vec2f barycentric_coordinates(
    const vec2f& point, const vec2f& a, const vec2f& b, const vec2f& c) {
  auto  v0 = b - a, v1 = c - a, v2 = point - a;
  float d00   = dot(v0, v0);
  float d01   = dot(v0, v1);
  float d11   = dot(v1, v1);
  float d20   = dot(v2, v0);
  float d21   = dot(v2, v1);
  float denom = d00 * d11 - d01 * d01;
  return vec2f{d11 * d20 - d01 * d21, d00 * d21 - d01 * d20} / denom;
};

// given a direction expressed in tangent space of the face start,
// continue the path as straight as possible.
geodesic_path straightest_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& start, const vec2f& direction, float path_length) {
  auto path  = geodesic_path{};
  path.start = start;
  path.strip.push_back(start.face);

  auto coords    = triangle_coordinates(triangles, positions, start);
  auto prev_face = -2, face = start.face;
  auto len = 0.0f;

  // https://rootllama.wordpress.com/2014/06/20/ray-line-segment-intersection-test-in-2d/
  auto intersect = [](const vec2f& direction, const vec2f& left,
                       const vec2f& right) {
    auto v1 = -left;
    auto v2 = right - left;
    auto v3 = vec2f{-direction.y, direction.x};
    auto t0 = cross(v2, v1) / dot(v2, v3);
    auto t1 = -dot(left, v3) / dot(v2, v3);
    return pair<float, float>{t0, t1};
  };

  while (len < path_length) {
    // Given the triangle, find which edge is intersected by the line.
    for (auto k = 0; k < 3; ++k) {
      auto neighbor = adjacencies[face][k];
      if (neighbor == prev_face) continue;
      auto left     = coords[k];
      auto right    = coords[mod3(k + 1)];
      auto [t0, t1] = intersect(direction, left, right);
      if (t0 > 0 && t1 >= 0 && t1 <= 1) {
        len = t0;
        if (t0 < path_length) {
          path.lerps.push_back(t1);
          // Step to next face.
          prev_face = face;
          if (neighbor == -1) {
            path_length = len;
            break;
          }
          coords = unfold_face(triangles, positions, coords, face, neighbor);
          face   = adjacencies[face][k];
          path.strip.push_back(face);
        }
        break;
      }
    }
  }

  auto p   = direction * path_length;
  auto uv  = barycentric_coordinates(p, coords[0], coords[1], coords[2]);
  path.end = {face, uv};
  return path;
}

mat2f parallel_transport_rotation(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_path& path) {
  if (path.start.face == path.end.face) return identity2x2f;

  auto coords   = unfold_strip(triangles, positions, path.strip, path.start);
  auto a        = coords.back()[0];
  auto b        = coords.back()[1];
  auto y_axis   = normalize(b - a);
  auto rotation = mat2f{};
  rotation.y    = y_axis;
  rotation.x    = {y_axis.y, -y_axis.x};
  rotation      = transpose(rotation);
  return rotation;
}

static float intersect_segments(const vec2f& start1, const vec2f& end1,
    const vec2f& start2, const vec2f& end2) {
  if (end1 == start2) return 0;
  if (end2 == start1) return 1;
  if (start2 == start1) return 0;
  if (end2 == end1) return 1;
  auto a   = end1 - start1;    // direction of line a
  auto b   = start2 - end2;    // direction of line b, reversed
  auto d   = start2 - start1;  // right-hand side
  auto det = a.x * b.y - a.y * b.x;
  assert(det);
  return (a.x * d.y - a.y * d.x) / det;
}

bool path_check_strip(
    const vector<vec3i>& adjacencies, const vector<int>& strip) {
  auto faces = unordered_set<int>{};
  faces.insert(strip[0]);
  for (auto i = 1; i < strip.size(); ++i) {
    assert(faces.count(strip[i]) == 0);  // face appears twice in the strip
    faces.insert(strip[i]);
    assert(find_in_vec(adjacencies[strip[i - 1]], strip[i]) != -1);
    assert(find_in_vec(adjacencies[strip[i]], strip[i - 1]) != -1);
  }
  return true;
}

struct funnel_point {
  int   face = 0;
  vec2f pos  = {0, 0};
};

static int max_curvature_point(const vector<funnel_point>& path) {
  // Among vertices around which the path curves, find the vertex
  // with maximum angle. We are going to fix that vertex. Actually, max_index is
  // the index of the first face containing that vertex.
  auto max_index = -1;
  auto max_angle = 0.0f;
  for (auto i = 1; i < path.size() - 1; ++i) {
    auto pos   = path[i].pos;
    auto prev  = path[i - 1].pos;
    auto next  = path[i + 1].pos;
    auto v0    = normalize(pos - prev);
    auto v1    = normalize(next - pos);
    auto angle = 1 - dot(v0, v1);
    if (angle > max_angle) {
      max_index = path[i].face;
      max_angle = angle;
    }
  }
  return max_index;
}

static vector<float> funnel(
    const vector<pair<vec2f, vec2f>>& portals, int& max_index) {
  // Find straight path.
  auto start       = vec2f{0, 0};
  auto apex_index  = 0;
  auto left_index  = 0;
  auto right_index = 0;
  auto apex        = start;
  auto left_bound  = portals[0].first;
  auto right_bound = portals[0].second;

  // Add start point.
  auto points = vector<funnel_point>{{apex_index, apex}};
  points.reserve(portals.size());

  // @Speed: is this slower than an inlined function?
  auto area = [](const vec2f a, const vec2f b, const vec2f c) {
    return cross(b - a, c - a);
  };

  for (auto i = 1; i < portals.size(); ++i) {
    auto left = portals[i].first, right = portals[i].second;
    // Update right vertex.
    if (area(apex, right_bound, right) <= 0) {
      if (apex == right_bound || area(apex, left_bound, right) > 0) {
        // Tighten the funnel.
        right_bound = right;
        right_index = i;
      } else {
        // Right over left, insert left to path and restart scan from
        // portal left point.
        if (left_bound != apex) {
          points.push_back({left_index, left_bound});
          // Make current left the new apex.
          apex       = left_bound;
          apex_index = left_index;
          // Reset portal
          left_bound  = apex;
          right_bound = apex;
          left_index  = apex_index;
          right_index = apex_index;
          // Restart scan
          i = apex_index;
          continue;
        }
      }
    }

    // Update left vertex.
    if (area(apex, left_bound, left) >= 0) {
      if (apex == left_bound || area(apex, right_bound, left) < 0) {
        // Tighten the funnel.
        left_bound = left;
        left_index = i;
      } else {
        if (right_bound != apex) {
          points.push_back({right_index, right_bound});
          // Make current right the new apex.
          apex       = right_bound;
          apex_index = right_index;
          // Reset portal
          left_bound  = apex;
          right_bound = apex;
          left_index  = apex_index;
          right_index = apex_index;
          // Restart scan
          i = apex_index;
          continue;
        }
      }
    }
  }

  // This happens when we got an apex on the last edge of the strip
  if (points.back().pos != portals.back().first) {
    points.push_back({(int)portals.size() - 1, portals.back().first});
  }
  assert(points.back().pos == portals.back().first);
  assert(points.back().pos == portals.back().second);

  auto lerps = vector<float>();
  lerps.reserve(portals.size());
  for (auto i = 0; i < points.size() - 1; i++) {
    auto a = points[i].pos;
    auto b = points[i + 1].pos;
    for (auto k = points[i].face; k < points[i + 1].face; k++) {
      auto portal = portals[k];
      auto s      = intersect_segments(a, b, portal.first, portal.second);
      assert(s >= -0.01f && s <= 1.01f);
      auto p = clamp(s, 0.0f, 1.0f);
      lerps.push_back(p);
    }
  }

  auto index = 1;
  for (auto i = 1; i < portals.size(); ++i) {
    if ((portals[i].first == points[index].pos) ||
        (portals[i].second == points[index].pos)) {
      points[index].face = i;
      index += 1;
    }
  }
  max_index = max_curvature_point(points);
  assert(lerps.size() == portals.size() - 1);
  return lerps;
}

bool check_point(const mesh_point& point) {
  assert(point.face != -1);
  assert(point.uv.x >= 0);
  assert(point.uv.y >= 0);
  assert(point.uv.x <= 1);
  assert(point.uv.y <= 1);
  return true;
}

static vector<int> fix_strip(const vector<vec3i>& adjacencies,
    const vector<int>& strip, int index, int k, bool left) {
  assert(path_check_strip(adjacencies, strip));
  assert(index < strip.size() - 1);
  auto face = strip[index];
  if (!left) k = mod3(k + 2);

  // Create triangle fan that starts at face, walks backward along the strip for
  // a while, exits and then re-enters back.
  auto fan = triangle_fan(adjacencies, face, k, left);
  if (fan.empty()) return strip;

  // strip in the array of faces and fan is a loop of faces which has partial
  // intersection with strip. We wan to remove the intersection from strip and
  // insert there the remaining part of fan, so that we have a new valid strip.
  auto first_strip_intersection = index;
  auto first_fan_intersection   = 0;
  for (auto i = 1; i < fan.size(); i++) {
    auto fan_index   = i;
    auto strip_index = max(index - i, 0);
    if (strip_index < 0) break;
    if (fan[fan_index] == strip[strip_index]) {
      first_strip_intersection = strip_index;
      first_fan_intersection   = fan_index;
    } else {
      break;
    }
  }
  auto second_strip_intersection = index;
  auto second_fan_intersection   = 0;
  for (auto i = 0; i < fan.size(); i++) {
    auto fan_index   = (int)fan.size() - 1 - i;
    auto strip_index = index + i + 1;
    if (strip_index >= (int)strip.size()) break;
    if (fan[fan_index] == strip[strip_index]) {
      second_strip_intersection = strip_index;
      second_fan_intersection   = fan_index;
    } else {
      break;
    }
  }

  if (first_strip_intersection >= second_strip_intersection) return strip;
  if (first_fan_intersection >= second_fan_intersection) return strip;

  auto result = vector<int>{};
  result.reserve(strip.size() + 12);

  // Initial part of original strip, up until intersection with fan.
  result.insert(
      result.end(), strip.begin(), strip.begin() + first_strip_intersection);

  // Append out-flanking part of fan.
  result.insert(result.end(), fan.begin() + first_fan_intersection,
      fan.begin() + second_fan_intersection);

  // Append remaining part of strip after intersection with fan.
  for (auto i = second_strip_intersection; i < strip.size(); ++i)
    result.push_back(strip[i]);

  assert(path_check_strip(adjacencies, result));
  return result;
}

static void straighten_path(geodesic_path& path, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies) {
  auto index = -1, vertex = -1;
  auto init_portals = unfold_funnel_portals(
      triangles, positions, path.strip, path.start, path.end);
  path.lerps = funnel(init_portals, index);

  // while(true) { this may never break...
  auto max_iterations = path.strip.size() * 2;
  for (auto i = 0; i < max_iterations && index != -1; i++) {
    auto new_vertex = -1;
    auto face       = path.strip[index];
    auto next       = path.strip[index + 1];
    auto edge       = common_edge(triangles[face], triangles[next]);
    auto flank_left = false;
    if (path.lerps[index] == 0) {
      new_vertex = edge.x;
      flank_left = false;
    } else if (path.lerps[index] == 1) {
      new_vertex = edge.y;
      flank_left = true;
    }
    if (new_vertex == vertex) break;
    vertex = new_vertex;

    path.strip = fix_strip(adjacencies, path.strip, index,
        find_in_vec(triangles[face], vertex), flank_left);

    auto portals = unfold_funnel_portals(
        triangles, positions, path.strip, path.start, path.end);
    path.lerps = funnel(portals, index);
  }
}

geodesic_path shortest_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& start, const mesh_point& end, const vector<int>& strip) {
  auto path  = geodesic_path{};
  path.start = start;
  path.end   = end;
  path.strip = strip;
  straighten_path(path, triangles, positions, adjacencies);
  return path;
}

mesh_path convert_mesh_path(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<int>& strip,
    const vector<float>& lerps, const mesh_point& start,
    const mesh_point& end) {
  auto result = mesh_path{};
  result.points.resize(lerps.size() + 2);
  result.points[0] = start;
  for (int i = 0; i < lerps.size(); ++i) {
    auto  k              = find_in_vec(adjacencies[strip[i]], strip[i + 1]);
    vec2f uvw[3]         = {{0, 0}, {1, 0}, {0, 1}};
    auto  a              = uvw[k];
    auto  b              = uvw[mod3(k + 1)];
    auto  uv             = lerp(a, b, lerps[i]);
    result.points[i + 1] = {strip[i], uv};
  }
  result.points.back() = end;
  return result;
}

mesh_point eval_path_point(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<float>& parameter_t,
    float t) {
  if (t <= 0) return path.start;
  if (t >= 1) return path.end;

  // strip with 1 triangle are trivial, just average the uvs
  if (path.start.face == path.end.face) {
    return mesh_point{path.start.face, lerp(path.start.uv, path.end.uv, t)};
  }
  // util function
  auto rotate = [](const vec3f& v, int k) {
    if (mod3(k) == 0)
      return v;
    else if (mod3(k) == 1)
      return vec3f{v.z, v.x, v.y};
    else
      return vec3f{v.y, v.z, v.x};
  };

  // TODO(splinesurf): check which is better betwee linear saerch, binary search
  // and "linear fit" search.
#if 1
  // linear search
  auto i = 0;
  for (; i < parameter_t.size() - 1; i++) {
    if (parameter_t[i + 1] >= t) break;
  }
#else
  // binary search
  // (but "linear fit search" should be faster here)
  // https://blog.demofox.org/2019/03/22/linear-fit-search/
  auto i      = -1;
  auto i_low  = 0;
  auto i_high = (int)parameter_t.size() - 1;
  while (true) {
    i           = (i_high + i_low) / 2;
    auto t_low  = parameter_t[i];
    auto t_high = parameter_t[i + 1];
    if (t_low <= t && t_high >= t) break;

    if (t_low > t) {
      i_high = i;
    } else {
      i_low = i;
    }
  }
#endif

  auto t_low  = parameter_t[i];
  auto t_high = parameter_t[i + 1];
  auto alpha  = (t - t_low) / (t_high - t_low);
  auto face   = path.strip[i];
  auto uv_low = vec2f{0, 0};
  if (i == 0) {
    uv_low = path.start.uv;
  } else {
    auto uvw = lerp(vec3f{1, 0, 0}, vec3f{0, 1, 0}, 1 - path.lerps[i - 1]);
    auto prev_face = path.strip[i - 1];
    auto k         = find_in_vec(adjacencies[face], prev_face);
    uvw            = rotate(uvw, k + 2);
    uv_low         = {uvw.x, uvw.y};
  }
  auto uv_high = vec2f{0, 0};
  if (i == parameter_t.size() - 2) {
    uv_high = path.end.uv;
  } else {
    auto uvw       = lerp(vec3f{1, 0, 0}, vec3f{0, 1, 0}, path.lerps[i]);
    auto next_face = path.strip[i + 1];
    auto k         = find_in_vec(adjacencies[face], next_face);
    uvw            = rotate(uvw, k + 2);
    uv_high        = {uvw.x, uvw.y};
  }
  auto uv = lerp(uv_low, uv_high, alpha);
  return mesh_point{face, uv};
}

mesh_point eval_path_point(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, float t) {
  auto parameter_t = path_parameters(
      path_positions(path, triangles, positions, adjacencies));
  return eval_path_point(
      path, triangles, positions, adjacencies, parameter_t, t);
}

[[maybe_unused]] static mesh_point eval_path_midpoint(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies) {
  return eval_path_point(path, triangles, positions, adjacencies, 0.5);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MESH BEZIER
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
static int find_in_vector(const T& vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x) return i;
  return -1;
}

[[maybe_unused]] static bool check_strip(
    const vector<vec3i>& adjacencies, const vector<int>& strip) {
  auto faces = std::unordered_set<int>{};
  faces.insert(strip[0]);
  for (auto i = 1; i < strip.size(); ++i) {
    assert(faces.count(strip[i]) == 0);  // face appears twice in the strip
    faces.insert(strip[i]);
    assert(find_in_vector(adjacencies[strip[i - 1]], strip[i]) != -1);
    assert(find_in_vector(adjacencies[strip[i]], strip[i - 1]) != -1);
  }
  return true;
}

template <typename Update, typename Stop, typename Exit>
static void search_strip(vector<float>& field, vector<bool>& in_queue,
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int start, int end, Update&& update,
    Stop&& stop, Exit&& exit) {
  auto destination_pos = eval_position(
      triangles, positions, {end, {1.0f / 3, 1.0f / 3}});

  auto estimate_dist = [&](int face) {
    auto p = eval_position(triangles, positions, {face, {1.0f / 3, 1.0f / 3}});
    return length(p - destination_pos);
  };
  field[start] = estimate_dist(start);

  // Cumulative weights of elements in queue. Used to keep track of the
  // average weight of the queue.
  double cumulative_weight = 0.0;

  // setup queue
  auto queue      = std::deque<int>{};
  in_queue[start] = true;
  cumulative_weight += field[start];
  queue.push_back(start);

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)(cumulative_weight / queue.size());

    // Large Label Last (see comment at the beginning)
    for (auto tries = 0; tries < queue.size() + 1; tries++) {
      if (field[node] <= average_weight) break;
      queue.pop_front();
      queue.push_back(node);
      node = queue.front();
    }

    // Remove node from queue.
    queue.pop_front();
    in_queue[node] = false;
    cumulative_weight -= field[node];

    // Check early exit condition.
    if (exit(node)) break;
    if (stop(node)) continue;

    for (auto i = 0; i < (int)solver.graph[node].size(); i++) {
      auto neighbor = solver.graph[node][i].node;
      if (neighbor == -1) continue;

      // Distance of neighbor through this node
      auto new_distance = field[node];
      new_distance += solver.graph[node][i].length;
      new_distance += estimate_dist(neighbor);
      new_distance -= estimate_dist(node);

      auto old_distance = field[neighbor];
      if (new_distance >= old_distance) continue;

      if (in_queue[neighbor]) {
        // If neighbor already in queue, don't add it.
        // Just update cumulative weight.
        cumulative_weight += new_distance - old_distance;
      } else {
        // If neighbor not in queue, add node to queue using Small Label
        // First (see comment at the beginning).
        if (queue.empty() || (new_distance < field[queue.front()]))
          queue.push_front(neighbor);
        else
          queue.push_back(neighbor);

        // Update queue information.
        in_queue[neighbor] = true;
        cumulative_weight += new_distance;
      }

      // Update distance of neighbor.
      field[neighbor] = new_distance;
      if (update(node, neighbor, new_distance)) return;
    }
  }
}

static vector<int> compute_strip(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions, int start,
    int end) {
  if (start == end) return {start};

  thread_local static auto parents  = vector<int>{};
  thread_local static auto field    = vector<float>{};
  thread_local static auto in_queue = vector<bool>{};

  if (parents.size() != solver.graph.size()) {
    parents.assign(solver.graph.size(), -1);
    field.assign(solver.graph.size(), flt_max);
    in_queue.assign(solver.graph.size(), false);
  }

  // initialize once for all and sparsely cleanup at the end of every solve
  auto visited = vector<int>{start};
  auto sources = vector<int>{start};
  auto update  = [&visited, end](int node, int neighbor, float new_distance) {
    parents[neighbor] = node;
    visited.push_back(neighbor);
    return neighbor == end;
  };
  auto stop = [](int node) { return false; };
  auto exit = [](int node) { return false; };

  search_strip(field, in_queue, solver, triangles, positions, start, end,
      update, stop, exit);

  // extract_strip
  auto strip = vector<int>{};
  auto node  = end;
  strip.reserve((int)yocto::sqrt((float)parents.size()));

  // for (int i = 0; i < parents.size() && node != -1; i++) {
  while (node != -1) {
    assert(find_in_vector(strip, node) != 1);
    assert(strip.size() < parents.size());
    strip.push_back(node);
    node = parents[node];
  }

  // cleanup buffers
  for (auto& v : visited) {
    parents[v]  = -1;
    field[v]    = flt_max;
    in_queue[v] = false;
  }
  return strip;
}

static geodesic_path compute_geodesic_path(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& start,
    const mesh_point& end) {
  check_point(start);
  check_point(end);
  auto path = geodesic_path{};
  if (start.face == end.face) {
    path.start = start;
    path.end   = end;
    path.strip = {start.face};
    return path;
  }
  auto strip = compute_strip(
      solver, triangles, positions, end.face, start.face);
  path = shortest_path(triangles, positions, adjacencies, start, end, strip);
  return path;
}
static mesh_point geodesic_lerp(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& start,
    const mesh_point& end, float t) {
  // profile_function();
  if (start.face == end.face) {
    return mesh_point{start.face, lerp(start.uv, end.uv, t)};
  }
  auto path = compute_geodesic_path(
      solver, triangles, positions, adjacencies, start, end);
  auto point = eval_path_point(path, triangles, positions, adjacencies, t);
  assert(check_point(point));
  return point;
}

using spline_polygon = array<mesh_point, 4>;

static pair<spline_polygon, spline_polygon> subdivide_bezier_polygon(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& input, float t) {
  auto Q0 = geodesic_lerp(
      solver, triangles, positions, adjacencies, input[0], input[1], t);
  auto Q1 = geodesic_lerp(
      solver, triangles, positions, adjacencies, input[1], input[2], t);
  auto Q2 = geodesic_lerp(
      solver, triangles, positions, adjacencies, input[2], input[3], t);
  auto R0 = geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, t);
  auto R1 = geodesic_lerp(solver, triangles, positions, adjacencies, Q1, Q2, t);
  auto S  = geodesic_lerp(solver, triangles, positions, adjacencies, R0, R1, t);
  return {{input[0], Q0, R0, S}, {S, R1, Q2, input[3]}};
}

static vector<mesh_point> de_casteljau_uniform(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, int subdivisions) {
  auto segments = vector<spline_polygon>{control_points};
  auto result   = vector<spline_polygon>();
  for (auto subdivision = 0; subdivision < subdivisions; subdivision++) {
    result.resize(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      auto [split0, split1] = subdivide_bezier_polygon(
          solver, triangles, positions, adjacencies, segments[i], 0.5);
      result[2 * i]     = split0;
      result[2 * i + 1] = split1;
    }
    swap(segments, result);
  }
  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 4};
}

vector<mesh_point> compute_bezier_path(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>&        adjacencies,
    const array<mesh_point, 4>& control_points, int subdivisions) {
  return de_casteljau_uniform(
      solver, triangles, positions, adjacencies, control_points, subdivisions);
}

vector<mesh_point> compute_bezier_path(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& control_points,
    int subdivisions) {
  auto path = vector<mesh_point>{};
  for (auto idx = 0; idx < (int)control_points.size() - 3; idx += 3) {
    auto polygon = spline_polygon{control_points[idx + 0],
        control_points[idx + 1], control_points[idx + 2],
        control_points[idx + 3]};
    auto segment = de_casteljau_uniform(
        solver, triangles, positions, adjacencies, polygon, subdivisions);
    path.insert(path.end(), segment.begin(), segment.end());
  }
  return path;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MESH BEZIER
// -----------------------------------------------------------------------------
namespace yocto {

static mesh_point geodesic_midpoint(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& start,
    const mesh_point& end) {
  // profile_function();

  if (start.face == end.face) {
    return mesh_point{start.face, (start.uv + end.uv) * 0.5};
  }

  auto path = compute_geodesic_path(
      solver, triangles, positions, adjacencies, start, end);
  auto midpoint = eval_path_midpoint(path, triangles, positions, adjacencies);
  //  assert(check_point(midpoint));
  return midpoint;
}

static mesh_point geodesic_lerp(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& a, const mesh_point& b,
    const mesh_point& c, float t0, float t1) {
  // den := (1-t0-t1) + t0 = 1 - t1;
  auto t  = t0 / (1 - t1);
  auto ab = geodesic_lerp(solver, triangles, positions, adjacencies, a, b, t);
  return geodesic_lerp(solver, triangles, positions, adjacencies, ab, c, t1);
}

using spline_polygon = array<mesh_point, 4>;

struct bezier_polygon {
  spline_polygon               segment;
  std::array<geodesic_path, 3> lines;
};

static mesh_point eval_path_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path, float t) {
  return eval_path_point(path, triangles, positions, adjacencies, t);
}

static mesh_point eval_path_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path,
    const vector<float>& path_parameter_t, const float& t) {
  return eval_path_point(
      path, triangles, positions, adjacencies, path_parameter_t, t);
}

static mesh_point eval_path_midpoint(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path) {
  return eval_path_point(path, triangles, positions, adjacencies, 0.5);
}

static vector<vec3f> path_positions(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path) {
  return path_positions(path, triangles, positions, adjacencies);
}

static float path_length(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path) {
  auto positions_ = path_positions(
      solver, triangles, positions, adjacencies, path);

  auto result = 0.0f;
  for (int i = 1; i < positions_.size(); ++i)
    result += length(positions_[i] - positions_[i - 1]);

  return result;
}

// TODO(giacomo): put in yocto_mesh (explode params)
static vec2f tangent_path_direction(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& path,
    bool start = true) {
  auto find = [](const vec3i& vec, int x) {
    for (int i = 0; i < size(vec); i++)
      if (vec[i] == x) return i;
    return -1;
  };

  auto direction = vec2f{};

  if (start) {
    auto start_tr = triangle_coordinates(triangles, positions, path.start);

    if (path.lerps.empty()) {
      direction = interpolate_triangle(
          start_tr[0], start_tr[1], start_tr[2], path.end.uv);
    } else {
      auto x    = path.lerps[0];
      auto k    = find(adjacencies[path.strip[0]], path.strip[1]);
      direction = lerp(start_tr[k], start_tr[(k + 1) % 3], x);
    }
  } else {
    auto end_tr = triangle_coordinates(triangles, positions, path.end);
    if (path.lerps.empty()) {
      direction = interpolate_triangle(
          end_tr[0], end_tr[1], end_tr[2], path.start.uv);
    } else {
      auto x = path.lerps.back();
      auto k = find(
          adjacencies[path.strip.rbegin()[0]], path.strip.rbegin()[1]);
      direction = lerp(end_tr[k], end_tr[(k + 1) % 3], 1 - x);
    }
  }
  return normalize(direction);
}

// TODO(fabio): compare with subdivided_bezier_polygon
static pair<bezier_polygon, bezier_polygon> subdivide_bezier(
    const bezier_polygon& input, const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies) {
  auto& segment = input.segment;
  auto  Q0      = eval_path_midpoint(
      solver, triangles, positions, adjacencies, input.lines[0]);
  auto Q1 = eval_path_midpoint(
      solver, triangles, positions, adjacencies, input.lines[1]);
  auto Q2 = eval_path_midpoint(
      solver, triangles, positions, adjacencies, input.lines[2]);

  auto R0 = geodesic_midpoint(
      solver, triangles, positions, adjacencies, Q0, Q1);
  auto R1 = geodesic_midpoint(
      solver, triangles, positions, adjacencies, Q1, Q2);
  auto S = geodesic_midpoint(solver, triangles, positions, adjacencies, R0, R1);
  auto result           = pair<bezier_polygon, bezier_polygon>{};
  result.first.lines[0] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, segment[0], Q0);
  result.first.lines[1] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, Q0, R0);
  result.first.lines[2] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, R0, S);
  result.first.segment   = {segment[0], Q0, R0, S};
  result.second.lines[0] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, S, R1);
  result.second.lines[1] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, R1, Q2);
  result.second.lines[2] = compute_geodesic_path(
      solver, triangles, positions, adjacencies, Q2, segment[3]);
  result.second.segment = {S, R1, Q2, segment[1]};
  return result;
  // return {{P0, Q0, R0, S}, {S, R1, Q2, P3}};
}

static bool is_control_polygon_unfoldable(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& segment) {
  if (segment[0].face != segment[1].face) return false;
  if (segment[1].face != segment[2].face) return false;
  if (segment[2].face != segment[3].face) return false;
  return true;
};

static mesh_point de_casteljau_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& segment, float t0,
    float t_start, float t_end) {
  auto is_control_polygon_unfoldable =
      [](const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
          const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
          const spline_polygon& segment) -> bool {
    if (segment[0].face != segment[1].face) return false;
    if (segment[1].face != segment[2].face) return false;
    if (segment[2].face != segment[3].face) return false;
    return true;
  };

  auto eval_bezier_point_cheap =
      [](const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
          const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
          const spline_polygon& polygon, float t) -> mesh_point {
    auto Q0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[0], polygon[1], t);
    auto Q1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[1], polygon[2], t);
    auto Q2 = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[2], polygon[3], t);
    auto R0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q0, Q1, t);
    auto R1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q1, Q2, t);
    return geodesic_lerp(solver, triangles, positions, adjacencies, R0, R1, t);
  };

  auto points = segment;

  // Recursively subdivide, shrinking the active control polygon around the
  // point of interest, until it is unfoldable.
  while (true) {
    auto Q0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[0], points[1], 0.5);
    auto Q1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[1], points[2], 0.5);
    auto Q2 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[2], points[3], 0.5);
    auto R0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
    auto R1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q1, Q2, 0.5);
    auto S = geodesic_lerp(
        solver, triangles, positions, adjacencies, R0, R1, 0.5);
    auto mid_t = (t_start + t_end) / 2;
    if (t0 < mid_t) {
      points[1] = Q0;
      points[2] = R0;
      points[3] = S;
      t_end     = mid_t;
    } else {
      points[0] = S;
      points[1] = R1;
      points[2] = Q2;
      t_start   = mid_t;
    }

    if (is_control_polygon_unfoldable(
            solver, triangles, positions, adjacencies, points))
      break;
  }

  // Exact evaluation.
  float t = (t0 - t_start) / (t_end - t_start);
  return eval_bezier_point_cheap(
      solver, triangles, positions, adjacencies, points, t);
}

static mesh_point lane_riesenfeld_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& polygon, float t0,
    float precision);

mesh_point eval_bezier_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& segment, float t,
    bool lane_riesenfeld, float precision) {
  if (lane_riesenfeld) {
    return lane_riesenfeld_point(
        solver, triangles, positions, adjacencies, segment, t, precision);
  } else {
    return de_casteljau_point(
        solver, triangles, positions, adjacencies, segment, t, 0, 1);
  }
}

static spline_polygon from_spline_to_bezier(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<float>& interval,
    const spline_polygon& polygon) {
  vector<float>  deltas(5);
  spline_polygon spline_polygon;
  for (auto i = 1; i < 6; ++i) {
    deltas[i] = interval[i] - interval[i - 1];
  }
  auto a = (yocto::pow(deltas[2], 2)) /
           ((deltas[0] + deltas[1] + deltas[2]) * (deltas[1] + deltas[2]));
  auto b = (yocto::pow(deltas[2], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  auto c = (yocto::pow(deltas[2], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[2] + deltas[3]));
  auto d = (yocto::pow(deltas[2], 2)) /
           ((deltas[2] + deltas[3] + deltas[4]) * (deltas[3] + deltas[4]));
  auto e = (deltas[1] * deltas[2]) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  auto f = (yocto::pow(deltas[1], 2)) /
           ((deltas[1] + deltas[2] + deltas[3]) * (deltas[1] + deltas[2]));
  if (a != 0 && 1 - a - f != 0 && f != 0)
    spline_polygon[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[0], polygon[1], polygon[2], 1 - a - f, f);
  else if (a != 0 && 1 - a != 0)
    spline_polygon[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[0], polygon[1], 1 - a);
  else if (a != 0 && f != 0)
    spline_polygon[0] = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[0], polygon[2], f);
  else if (a != 0)
    spline_polygon[0] = polygon[0];
  else if (f != 0 && 1 - f != 0)
    spline_polygon[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[1], polygon[2], 1 - f);
  else if (f != 0)
    spline_polygon[0] = polygon[2];
  else
    spline_polygon[0] = polygon[1];

  if (1 - e - f != 0 && e + f != 0)
    spline_polygon[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[1], polygon[2], e + f);
  else if (1 - e - f != 0)
    spline_polygon[1] = polygon[1];
  else
    spline_polygon[1] = polygon[2];

  if (1 - 2 * e - b - f != 0 && b + 2 * e + f != 0)
    spline_polygon[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[1], polygon[2], b + 2 * e + f);
  else if (1 - 2 * e - b - f != 0)
    spline_polygon[2] = polygon[1];
  else
    spline_polygon[2] = polygon[2];

  if (1 - 3 * e - 2 * b + c - f != 0 && 2 * b - c - d + 3 * e + f != 0 &&
      d != 0)
    spline_polygon[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[1], polygon[2], polygon[3], 2 * c - c - d + 3 * e + f, d);
  else if (1 - 3 * e - 2 * b + c - f != 0 && 2 * b - c + 3 * e + f != 0)
    spline_polygon[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
        polygon[1], polygon[2], 2 * b - c + 3 * e + f);
  else if (1 - 3 * e - 2 * b + c - f != 0 && d != 0)
    spline_polygon[3] = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[1], polygon[2], d);
  else if (1 - 3 * e - 2 * b + c - f != 0)
    spline_polygon[3] = polygon[1];
  else if (2 * b - c - d + 3 * e + f != 0 && d != 0)
    spline_polygon[3] = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[2], polygon[3], d);
  else if (d != 0)
    spline_polygon[3] = polygon[3];
  else
    spline_polygon[3] = polygon[2];

  return spline_polygon;
}

static mesh_point de_boor_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& polygon,
    const vector<float>& knot_vector, const float& t0) {
  spline_polygon old_points  = polygon;
  spline_polygon curr_points = {};
  for (auto j = 1; j <= 3; ++j) {
    for (auto i = 3; i >= j; --i) {
      auto alpha = (t0 - knot_vector[i]) /
                   (knot_vector[i - j + 4] - knot_vector[i]);
      curr_points[i] = geodesic_lerp(solver, triangles, positions, adjacencies,
          old_points[i - 1], old_points[i], alpha);
    }
    old_points = curr_points;
  }
  return curr_points.back();
}

static bool spline_stop_criterion(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& polygon,
    const float& treshold) {
  auto L01 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, polygon[0], polygon[1]);
  auto L12 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, polygon[1], polygon[2]);
  auto L23 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, polygon[2], polygon[3]);
  auto L03 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, polygon[0], polygon[3]);
  auto perimeter = path_length(solver, triangles, positions, adjacencies, L01) +
                   path_length(solver, triangles, positions, adjacencies, L12) +
                   path_length(solver, triangles, positions, adjacencies, L23);
  if (perimeter - path_length(solver, triangles, positions, adjacencies, L03) <=
      treshold)
    return true;

  return false;
}

static mesh_point lane_riesenfeld_regular(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& a, const mesh_point& b,
    const mesh_point& c) {
  auto Q0 = geodesic_lerp(
      solver, triangles, positions, adjacencies, a, b, 0.75);
  auto Q1 = geodesic_lerp(
      solver, triangles, positions, adjacencies, b, c, 0.25);
  return geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
}

static mesh_point lane_riesenfeld_regular(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& l0,
    const geodesic_path& l1) {
  auto Q0 = eval_path_point(
      solver, triangles, positions, adjacencies, l0, 0.75);
  auto Q1 = eval_path_point(
      solver, triangles, positions, adjacencies, l1, 0.25);
  return geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
}

static mesh_point lane_riesenfeld_init(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& a, const mesh_point& b,
    const mesh_point& c) {
  auto Q0 = geodesic_lerp(
      solver, triangles, positions, adjacencies, a, b, 5 / 8.f);
  auto Q1 = geodesic_lerp(
      solver, triangles, positions, adjacencies, b, c, 3 / 8.f);
  return geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
}

static mesh_point lane_riesenfeld_boundary(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& a, const mesh_point& b,
    const mesh_point& c, const bool& left) {
  auto Q0 = mesh_point{};
  auto Q1 = mesh_point{};
  if (left) {
    Q0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, a, b, 5 / 8.f);
    Q1 = geodesic_lerp(solver, triangles, positions, adjacencies, b, c, 0.25);

  } else {
    Q0 = geodesic_lerp(solver, triangles, positions, adjacencies, a, b, 0.75);
    Q1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, b, c, 3 / 8.f);
  }
  return geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
}

static mesh_point lane_riesenfeld_boundary(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const geodesic_path& l0,
    const geodesic_path& l1, const bool& left) {
  auto Q0 = mesh_point{};
  auto Q1 = mesh_point{};
  if (left) {
    Q0 = eval_path_point(
        solver, triangles, positions, adjacencies, l0, 5 / 8.f);
    Q1 = eval_path_point(solver, triangles, positions, adjacencies, l1, 0.25);

  } else {
    Q0 = eval_path_point(solver, triangles, positions, adjacencies, l0, 0.75);
    Q1 = eval_path_point(
        solver, triangles, positions, adjacencies, l1, 3 / 8.f);
  }
  return geodesic_lerp(solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
}

static pair<bool, spline_polygon> handle_boundary(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, const vector<int>& new_ones_entries,
    const int k) {
  spline_polygon new_ones = {};
  if (new_ones_entries[0] > 3 && new_ones_entries.back() < yocto::pow(2, k) - 1)
    return {false, {}};
  else if (new_ones_entries[0] <= 3) {
    switch (new_ones_entries[0]) {
      case 0: {
        new_ones[0] = control_points[0];
        new_ones[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[0], control_points[1], 0.5);
        new_ones[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[1], control_points[2], 0.25);
        new_ones[3] = lane_riesenfeld_boundary(solver, triangles, positions,
            adjacencies, control_points[1], control_points[2],
            control_points[3], true);
      } break;
      case 1: {
        new_ones[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[0], control_points[1], 0.5);
        new_ones[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[1], control_points[2], 0.25);
        new_ones[2] = lane_riesenfeld_boundary(solver, triangles, positions,
            adjacencies, control_points[1], control_points[2],
            control_points[3], true);
        new_ones[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[2], control_points[3], 0.5);

      } break;
      case 2: {
        new_ones[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[0], control_points[1], 0.25);
        new_ones[1] = lane_riesenfeld_boundary(solver, triangles, positions,
            adjacencies, control_points[0], control_points[1],
            control_points[2], true);
        new_ones[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[1], control_points[2], 0.5);
        new_ones[3] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, control_points[1], control_points[2],
            control_points[3]);

      } break;

      case 3: {
        new_ones[0] = lane_riesenfeld_boundary(solver, triangles, positions,
            adjacencies, control_points[0], control_points[1],
            control_points[2], true);
        new_ones[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[1], control_points[2], 0.5);
        new_ones[2] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, control_points[1], control_points[2],
            control_points[3]);
        new_ones[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
            control_points[2], control_points[3], 0.5);

      } break;
    }
  } else {
    if (new_ones_entries.back() == yocto::pow(2, k) + 2) {
      new_ones[0] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, control_points[0], control_points[1], control_points[2],
          false);
      new_ones[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[1], control_points[2], 0.75);
      new_ones[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[2], control_points[3], 0.5);
      new_ones[3] = control_points[3];
    } else if (new_ones_entries.back() == yocto::pow(2, k) + 1) {
      new_ones[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[0], control_points[1], 0.5);
      new_ones[1] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, control_points[0], control_points[1], control_points[2],
          false);
      new_ones[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[1], control_points[2], 0.75);
      new_ones[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[2], control_points[3], 0.5);
    } else if (new_ones_entries.back() == yocto::pow(2, k)) {
      new_ones[0] = lane_riesenfeld_regular(solver, triangles, positions,
          adjacencies, control_points[0], control_points[1], control_points[2]);
      new_ones[1] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[1], control_points[2], 0.5);
      new_ones[2] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, control_points[1], control_points[2], control_points[3],
          false);
      new_ones[3] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[2], control_points[3], 0.75);
    } else if (new_ones_entries.back() == yocto::pow(2, k) - 1) {
      new_ones[0] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[0], control_points[1], 0.5);
      new_ones[1] = lane_riesenfeld_regular(solver, triangles, positions,
          adjacencies, control_points[0], control_points[1], control_points[2]);
      new_ones[2] = geodesic_lerp(solver, triangles, positions, adjacencies,
          control_points[1], control_points[2], 0.5);
      new_ones[3] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, control_points[1], control_points[2], control_points[3],
          false);
    }
  }
  return {true, new_ones};
}

static std::tuple<int, vec2f, spline_polygon> find_leaf(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, const float& t0,
    const float& treshold) {
  auto q = vector<mesh_point>(7);

  auto& p = control_points;
  q[0]    = p[0];
  q[1]    = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[0], p[1], 1 / 4.f);
  auto p0p1 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[0], p[1], 1 / 2.f);
  auto p1p2 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[1], p[2], 1 / 2.f);
  q[2] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p0p1, p1p2, 1 / 4.f);
  auto p2p3 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[2], p[3], 1 / 2.f);
  q[3] = lane_riesenfeld_init(
      solver, triangles, positions, adjacencies, p0p1, p1p2, p2p3);
  q[4] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p1p2, p2p3, 3 / 4.f);
  q[5] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p2p3, p[3], 1 / 2.f);
  q[6] = p[3];

  spline_polygon leaf = {};
  auto           t    = zero2f;
  if (t0 <= 0.25) {
    leaf = {q[0], q[1], q[2], q[3]};
    t    = {0, 0.25};
  } else if (t0 <= 0.5) {
    leaf = {q[1], q[2], q[3], q[4]};
    t    = {0.25, 0.5};
  } else if (t0 <= 0.75) {
    leaf = {q[2], q[3], q[4], q[5]};
    t    = {0.5, 0.75};
  } else if (t0 < 1) {
    leaf = {q[3], q[4], q[5], q[6]};
    t    = {0.75, 1};
  } else
    assert(false);
  auto k = 3;
  while (true) {
    auto curr_t     = (t.x + t.y) / 2;
    int  curr_entry = (int)pow(2.0f, (float)k) * curr_t;
    curr_entry += 3;
    assert(curr_entry % 2 == 0);
    if (t0 < curr_t) {
      auto curr_entries = vector<int>{
          curr_entry - 4, curr_entry - 3, curr_entry - 2, curr_entry - 1};
      auto [are_boundaries, new_ones] = handle_boundary(
          solver, triangles, positions, adjacencies, leaf, curr_entries, k);
      if (!are_boundaries) {
        new_ones[0] = geodesic_lerp(
            solver, triangles, positions, adjacencies, leaf[0], leaf[1], 0.5);
        new_ones[1] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, leaf[0], leaf[1], leaf[2]);
        new_ones[2] = geodesic_lerp(
            solver, triangles, positions, adjacencies, leaf[1], leaf[2], 0.5);
        new_ones[3] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, leaf[1], leaf[2], leaf[3]);
      }
      leaf = new_ones;
      t    = {t.x, curr_t};
    } else {
      auto curr_entries = vector<int>{
          curr_entry - 3, curr_entry - 2, curr_entry - 1, curr_entry};
      auto [are_boundaries, new_ones] = handle_boundary(
          solver, triangles, positions, adjacencies, leaf, curr_entries, k);
      if (!are_boundaries) {
        new_ones[0] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, leaf[0], leaf[1], leaf[2]);
        new_ones[1] = geodesic_lerp(
            solver, triangles, positions, adjacencies, leaf[1], leaf[2], 0.5);
        new_ones[2] = lane_riesenfeld_regular(solver, triangles, positions,
            adjacencies, leaf[1], leaf[2], leaf[3]);
        new_ones[3] = geodesic_lerp(
            solver, triangles, positions, adjacencies, leaf[2], leaf[3], 0.5);
      }
      leaf = new_ones;
      t    = {curr_t, t.y};
    }

    if (spline_stop_criterion(
            solver, triangles, positions, adjacencies, leaf, treshold))
      return {k, t, leaf};
    else
      ++k;
  }
}

static mesh_point eval_spline_point_cheap(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& control_polygon,
    const vector<float>& knot_vector, const float& t0, const int entry) {
  assert(entry >= 3);
  auto control_points = {control_polygon[entry], control_polygon[entry - 1],
      control_polygon[entry - 2], control_polygon[entry - 3]};
  vector<mesh_point> old_points = control_points;
  vector<mesh_point> curr_points(4);
  auto               offset = entry - 3;
  for (auto j = 1; j <= 3; ++j) {
    for (auto i = entry - 3 + j; i <= entry; ++i) {
      auto alpha = (t0 - knot_vector[i]) /
                   (knot_vector[i - j + 4] - knot_vector[i]);
      curr_points[i - offset] = geodesic_lerp(solver, triangles, positions,
          adjacencies, old_points[i - 1 - offset], old_points[i - offset],
          alpha);
    }
    old_points = curr_points;
  }
  return curr_points.back();
}

static mesh_point lane_riesenfeld_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& polygon, float t0,
    float precision) {
  if (t0 <= 0) return polygon[0];
  if (t0 >= 1) return polygon[3];
  auto treshold     = precision;
  auto [k, t, leaf] = find_leaf(
      solver, triangles, positions, adjacencies, polygon, t0, treshold);
  float         step  = 1 / yocto::pow(2, k);
  vector<float> knots = {
      yocto::max(t.x - 3 * step, 0.f),
      yocto::max(t.x - 2 * step, 0.f),
      yocto::max(t.x - step, 0.f),
      t.x,
      t.y,
      yocto::min(t.y + step, 1.f),
      yocto::min(t.y + 2 * step, 1.f),
      yocto::min(t.y + 3 * step, 1.f),
  };

  return de_boor_point(
      solver, triangles, positions, adjacencies, leaf, knots, t0);
}

static std::array<spline_polygon, 2> de_casteljau_insert(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& polygon, float t0) {
  auto t_start = 0.f;
  auto t_end   = 1.f;
  auto points  = polygon;
  while (true) {
    auto Q0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[0], points[1], 0.5);
    auto Q1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[1], points[2], 0.5);
    auto Q2 = geodesic_lerp(
        solver, triangles, positions, adjacencies, points[2], points[3], 0.5);
    auto R0 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q0, Q1, 0.5);
    auto R1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, Q1, Q2, 0.5);
    auto S = geodesic_lerp(
        solver, triangles, positions, adjacencies, R0, R1, 0.5);
    auto mid_t = (t_start + t_end) / 2.f;
    if (t0 < mid_t) {
      points[1] = Q0;
      points[2] = R0;
      points[3] = S;
      t_end     = mid_t;
    } else {
      points[0] = S;
      points[1] = R1;
      points[2] = Q2;
      t_start   = mid_t;
    }
    if (is_control_polygon_unfoldable(
            solver, triangles, positions, adjacencies, points))
      break;
  }
  // Compute the parameter t local to the leaf control polygon.
  auto tP_local = (t0 - t_start) / (t_end - t_start);
  // Subdivide the leaf control with De Castljeau creating two new control
  // polygons. They are segment_left and segment_right.
  auto [segment_left, segment_right] = subdivide_bezier_polygon(
      solver, triangles, positions, adjacencies, points, tP_local);
  auto left_side  = compute_geodesic_path(solver, triangles, positions,
      adjacencies, segment_left.back(), segment_left[2]);
  auto right_side = compute_geodesic_path(solver, triangles, positions,
      adjacencies, segment_right[0], segment_right[1]);
  // P is the inserted mesh point that sepraters segment_left and
  // segment_right.
  // assert(segment_left[3] == segment_right[0]);
  auto P = segment_right[0];
  // left part
  {
    auto Pp2_len = path_length(
        solver, triangles, positions, adjacencies, left_side);
    auto Pp2_dir = tangent_path_direction(
        solver, triangles, positions, adjacencies, left_side);
    //    assert(left_leaf.start == P);
    auto delta_len = t_start * Pp2_len / (t0 - t_start);
    auto path      = straightest_path(
        triangles, positions, adjacencies, P, Pp2_dir, delta_len + Pp2_len);
    auto Pp2 = path.end;
    auto Pp1 = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[0], polygon[1], t0);
    segment_left = {polygon[0], Pp1, Pp2, P};
  }
  // right part
  {
    auto Pp1_len = path_length(
        solver, triangles, positions, adjacencies, right_side);
    auto Pp1_dir = tangent_path_direction(
        solver, triangles, positions, adjacencies, right_side);
    auto delta_len = (1 - t_end) / (t_end - t0) * Pp1_len;
    auto path      = straightest_path(
        triangles, positions, adjacencies, P, Pp1_dir, delta_len + Pp1_len);
    auto Pp1 = path.end;
    auto Pp2 = geodesic_lerp(
        solver, triangles, positions, adjacencies, polygon[2], polygon[3], t0);
    segment_right = {P, Pp1, Pp2, polygon[3]};
  }
  return {segment_left, segment_right};
}

static std::array<spline_polygon, 2> lane_riesenfeld_insert(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& polygon, float t0, float precision) {
  // Go down the tree and find the leaf node containing the point.
  std::array<spline_polygon, 2> result;

  if (t0 <= 0) {
    result[0] = {polygon[0], polygon[0], polygon[0], polygon[0]};
    result[1] = polygon;

  } else if (t0 >= 1) {
    result[0] = polygon;
    result[1] = {polygon[3], polygon[3], polygon[3], polygon[3]};

  } else {
    auto treshold         = precision;
    auto [depth, t, leaf] = find_leaf(
        solver, triangles, positions, adjacencies, polygon, t0, treshold);
    float         step           = 1 / yocto::pow(2, depth);
    vector<float> knots          = {yocto::max(t.x - 3 * step, 0.f),
        yocto::max(t.x - 2 * step, 0.f), yocto::max(t.x - step, 0.f), t.x, t.y,
        yocto::min(t.y + step, 1.f), yocto::min(t.y + 2 * step, 1.f),
        yocto::min(t.y + 3 * step, 1.f)};
    auto          spline_polygon = from_spline_to_bezier(
        solver, triangles, positions, adjacencies, knots, leaf);
    auto t_rel         = (t0 - t.x) / (t.y - t.x);
    auto [left, right] = subdivide_bezier_polygon(
        solver, triangles, positions, adjacencies, spline_polygon, t_rel);

    {
      auto L32 = compute_geodesic_path(
          solver, triangles, positions, adjacencies, left[3], left[2]);
      auto Pp2_len = path_length(
          solver, triangles, positions, adjacencies, L32);
      auto Pp2_dir = tangent_path_direction(
          solver, triangles, positions, adjacencies, L32);
      auto delta_len = t0 * Pp2_len / (t0 - t.x);
      auto path      = straightest_path(
          triangles, positions, adjacencies, left[3], Pp2_dir, delta_len);
      auto Pp2 = path.end;

      auto Pp1 = geodesic_lerp(solver, triangles, positions, adjacencies,
          polygon[0], polygon[1], t0);

      result[0] = {polygon[0], Pp1, Pp2, left[3]};
    }

    // right part
    {
      auto L01 = compute_geodesic_path(
          solver, triangles, positions, adjacencies, right[0], right[1]);
      auto Pp1_len = path_length(
          solver, triangles, positions, adjacencies, L01);
      auto Pp1_dir = tangent_path_direction(
          solver, triangles, positions, adjacencies, L01);
      auto t0_local = (t.y - t0) / (1 - t0);
      auto path = straightest_path(triangles, positions, adjacencies, right[0],
          Pp1_dir, Pp1_len / t0_local);
      auto Pp1  = path.end;
      auto Pp2  = geodesic_lerp(solver, triangles, positions, adjacencies,
          polygon[2], polygon[3], t.y);
      result[1] = {right[0], Pp1, Pp2, polygon[3]};
    }
  }
  return result;
}

array<spline_polygon, 2> insert_bezier_point(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& polygon, float t,
    bool lane_riesenfeld, float precision) {
  if (lane_riesenfeld) {
    return lane_riesenfeld_insert(
        solver, triangles, positions, adjacencies, polygon, t, precision);
  } else {
    throw de_casteljau_insert(
        solver, triangles, positions, adjacencies, polygon, t);
  }
}

static bool is_bezier_straight_enough(const geodesic_path& a,
    const geodesic_path& b, const geodesic_path& c,
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_params& params) {
  // TODO(giacomo): we don't need all positions!
  // auto a_positions = path_positions(solver, triangles, positions,
  // adjacencies,  a); auto b_positions = path_positions(solver, triangles,
  // positions, adjacencies,  b); auto c_positions = path_positions(solver,
  // triangles, positions, adjacencies,  c);

  {
    // On curve apex we may never reach straightess, so we check curve
    // length.
    auto pos  = array<vec3f, 4>{};
    pos[0]    = eval_position(triangles, positions, a.start);
    pos[1]    = eval_position(triangles, positions, b.start);
    pos[2]    = eval_position(triangles, positions, c.start);
    pos[3]    = eval_position(triangles, positions, c.end);
    float len = 0;
    for (int i = 0; i < 3; i++) {
      len += length(pos[i] - pos[i + 1]);
    }
    if (len < params.min_curve_size) return true;
  }

  {
    auto dir0 = tangent_path_direction(
        solver, triangles, positions, adjacencies, a, false);  // end
    auto dir1 = tangent_path_direction(
        solver, triangles, positions, adjacencies, b, true);  // start
    auto angle1 = cross(dir0, dir1);
    if (fabs(angle1) > params.precision) {
      // printf("a1: %f > %f\n", angle1, params.precision);
      return false;
    }
  }

  {
    auto dir0 = tangent_path_direction(
        solver, triangles, positions, adjacencies, b, false);  // end
    auto dir1 = tangent_path_direction(
        solver, triangles, positions, adjacencies, c, true);  // start
    auto angle1 = cross(dir0, dir1);
    if (fabs(angle1) > params.precision) {
      // printf("a2: %f > %f\n", angle1, params.precision);
      return false;
    }
  }

  return true;
}

static void subdivide_compute_bezier_de_casteljau_adaptive(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& input, spline_params params,
    vector<mesh_point>& result, int depth = 0) {
  // resulting beziers: (P0, Q0, R0, S) (S, R1, Q2, P3)

  if (depth > params.max_depth) {
    return;
  }
  auto [P0, P1, P2, P3] = input;

  auto P0_P1 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, P0, P1);
  auto P1_P2 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, P1, P2);
  auto P2_P3 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, P2, P3);

  if (is_bezier_straight_enough(P0_P1, P1_P2, P2_P3, solver, triangles,
          positions, adjacencies, params)) {
    result.push_back(P0);
    result.push_back(P1);
    result.push_back(P2);
    result.push_back(P3);
    return;
  }

  auto Q0    = eval_path_midpoint(P0_P1, triangles, positions, adjacencies);
  auto Q1    = eval_path_midpoint(P1_P2, triangles, positions, adjacencies);
  auto Q2    = eval_path_midpoint(P2_P3, triangles, positions, adjacencies);
  auto Q0_Q1 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, Q0, Q1);
  auto Q1_Q2 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, Q1, Q2);

  auto R0    = eval_path_midpoint(Q0_Q1, triangles, positions, adjacencies);
  auto R1    = eval_path_midpoint(Q1_Q2, triangles, positions, adjacencies);
  auto R0_R1 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, R0, R1);

  auto S = eval_path_midpoint(R0_R1, triangles, positions, adjacencies);

  subdivide_compute_bezier_de_casteljau_adaptive(solver, triangles, positions,
      adjacencies, {P0, Q0, R0, S}, params, result, depth + 1);
  subdivide_compute_bezier_de_casteljau_adaptive(solver, triangles, positions,
      adjacencies, {S, R1, Q2, P3}, params, result, depth + 1);
}

static vector<mesh_point> de_casteljau_adaptive(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, const spline_params& params) {
  auto result = vector<mesh_point>{};
  subdivide_compute_bezier_de_casteljau_adaptive(solver, triangles, positions,
      adjacencies, control_points, params, result);
  return result;
}

[[maybe_unused]] static vector<mesh_point> de_casteljau_uniform_mt(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, const spline_params& params) {
  auto segments = vector<spline_polygon>{control_points};
  auto result   = vector<spline_polygon>();
  auto threads  = vector<std::thread>((int)yocto::pow(2, params.subdivisions));

  auto f = [&](int k) {
    auto [split0, split1] = subdivide_bezier_polygon(
        solver, triangles, positions, adjacencies, segments[k], 0.5);
    result[k * 2]     = split0;
    result[k * 2 + 1] = split1;
  };

  for (auto subdivision = 0; subdivision < params.subdivisions; subdivision++) {
    result.resize(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      threads[i] = std::thread(f, i);
    }
    for (auto i = 0; i < segments.size(); i++) threads[i].join();

    swap(segments, result);
  }

  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 4};
}

static vector<mesh_point> de_casteljau_uniform(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_polygon& control_points, const spline_params& params) {
  auto segments = vector<spline_polygon>{control_points};
  auto result   = vector<spline_polygon>();

  for (auto subdivision = 0; subdivision < params.subdivisions; subdivision++) {
    result.resize(segments.size() * 2);
    for (auto i = 0; i < segments.size(); i++) {
      auto [split0, split1] = subdivide_bezier_polygon(
          solver, triangles, positions, adjacencies, segments[i], 0.5);
      result[2 * i]     = split0;
      result[2 * i + 1] = split1;
    }
    swap(segments, result);
  }
  return {(mesh_point*)segments.data(),
      (mesh_point*)segments.data() + segments.size() * 4};
}

static vector<mesh_point> lane_riesenfeld_uniform(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 4>& control_points, int num_subdivisions);
static vector<mesh_point> lane_riesenfeld_adaptive(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 4>& polygon, const spline_params& params);

[[maybe_unused]] static vector<vec3f> polyline_positions(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<mesh_point>& points) {
  auto result = vector<vec3f>();
  result.reserve(points.size() * 4);
  for (int i = 0; i < points.size() - 1; i++) {
    auto a    = points[i];
    auto b    = points[i + 1];
    auto path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, a, b);
    auto path_pos = path_positions(
        solver, triangles, positions, adjacencies, path);
    result.insert(result.end(), path_pos.begin(), path_pos.end());
    // append(result, path_pos);
  }
  return result;
}

vector<mesh_point> compute_bezier_path(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const spline_polygon& control_points,
    const spline_params& params) {
  switch (params.algorithm) {
    case spline_algorithm::de_casteljau_uniform: {
      return de_casteljau_uniform(
          solver, triangles, positions, adjacencies, control_points, params);
    }
    case spline_algorithm::de_casteljau_adaptive: {
      return de_casteljau_adaptive(
          solver, triangles, positions, adjacencies, control_points, params);
    }
    case spline_algorithm::lane_riesenfeld_uniform: {
      return lane_riesenfeld_uniform(solver, triangles, positions, adjacencies,
          control_points, params.subdivisions);
    }
    case spline_algorithm::lane_riesenfeld_adaptive: {
      return lane_riesenfeld_adaptive(
          solver, triangles, positions, adjacencies, control_points, params);
    }
  }
}

vector<mesh_point> compute_bezier_path(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& control_points,
    const spline_params& params) {
  auto path = vector<mesh_point>{};
  for (auto idx = 0; idx < (int)control_points.size() - 3; idx += 3) {
    auto polygon = spline_polygon{control_points[idx + 0],
        control_points[idx + 1], control_points[idx + 2],
        control_points[idx + 3]};
    auto segment = compute_bezier_path(
        solver, triangles, positions, adjacencies, polygon, params);
    path.insert(path.end(), segment.begin(), segment.end());
  }
  return path;
}

// weighted averages (Note:gradients needs to be a vector such that at the
// i-th entry constains the gradient field of the squared distance field from
// the i-th control points)
[[maybe_unused]] static vector<mesh_point> weighted_average(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<mesh_point>& rectangle, const vector<vector<float>>& weights) {
#if HEAVY
  vector<vector<vec3f>> gradients(rectangle.size());
  vector<vector<float>> field(rectangle.size());
  vector<mesh_point>    points_on_surface(weights.size(), {-1, zero2f});
  for (auto i = 0; i < rectangle.size(); ++i) {
    auto f = compute_distance_field(
        triangles, positions, adjacencies, v2t, solver, i, rectangle, graph);
    field[i] = f;
    std::transform(f.begin(), f.end(), f.begin(),
        [](float lambda) { return lambda * lambda; });
    gradients[i] = compute_grad(solver, triangles, positions, normals, Grad, f);
  }
  // set from as the control points associate to greatest weight
  vector<pair<vector<float>, int>> badones;
  for (auto i = 0; i < weights.size(); ++i) {
    auto from = optimal_seed(rectangle, weights[i]);
    auto p = gradient_descent(triangles, positions, adjacencies, v2t, normals,
        solver, angles, total_angles, gradients, weights[i], from);
    if (p.face != -1)
      points_on_surface[i] = p;
    else
      badones.push_back({weights[i], i});
  }
  if (points_on_surface.size() == 0) {
    std::cout << "rectangle is too big" << std::endl;
    return {};
  }
  bool stop              = badones.size() == 0;
  bool needs_subdivision = false;
  auto prev_proccessed   = badones.size();
  auto curr_processed    = 0;
  while (!stop) {
    stop = true;
    for (auto i = 0; i < badones.size(); ++i) {
      if (points_on_surface[badones[i].second].face == -1) {
        stop      = false;
        auto from = optimal_seed(points_on_surface, badones[i].first, weights);
        auto p    = gradient_descent(triangles, positions, adjacencies, v2t,
            normals, solver, angles, total_angles, gradients, badones[i].first,
            from);
        if (p.face != -1)
          points_on_surface[badones[i].second] = p;
        else
          ++curr_processed;
      }
      if (curr_processed == prev_proccessed) {
        stop              = true;
        needs_subdivision = true;

      } else {
        prev_proccessed = curr_processed;
        curr_processed  = 0;
      }
    }
  }
  // this part doesn't work yet
  if (needs_subdivision) {
    for (auto i = 0; i < badones.size(); ++i) {
      if (points_on_surface[badones[i].second].face == -1) {
        auto from = optimal_seed_with_weights(
            points_on_surface, badones[i].first, weights);
        auto curr_gradients = gradients_inside_cell(triangles, positions,
            adjacencies, v2t, solver, angles, normals, Grad, field, rectangle,
            from.first, from.second, badones[i].first);
        auto p = gradient_descent(triangles, positions, adjacencies, v2t,
            normals, solver, angles, total_angles, curr_gradients,
            badones[i].first, from.first);
        if (p.face != -1) {
          points_on_surface[badones[i].second] = p;
          std::cout << "we subdivide" << std::endl;
        } else {
          std::cout << "rectangle is too big" << std::endl;
        }
      }
    }
  }

  return points_on_surface;
#else
  return {};
#endif
}

struct spline_node {
  std::array<mesh_point, 4>    points  = {};
  std::array<geodesic_path, 3> lines   = {};
  vec2f                        t       = {};
  int                          depth   = 0;
  bool                         is_good = false;
};

static pair<bool, vector<mesh_point>> handle_boundary_node(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_node& leaf, const vector<int>& new_ones_entries) {
  vector<mesh_point> new_ones(5);
  auto               k = leaf.depth + 1;
  if (new_ones_entries[0] > 3 && new_ones_entries.back() < yocto::pow(2, k) - 1)
    return {false, new_ones};
  else if (new_ones_entries[0] <= 3) {
    if (new_ones_entries[0] == 0) {
      new_ones[0] = leaf.points[0];
      new_ones[1] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[0], 0.5);
      new_ones[2] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[1], 0.25);
      new_ones[3] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, leaf.lines[1], leaf.lines[2], true);
      new_ones[4] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[2], 0.5);
    } else if (new_ones_entries[0] == 2) {
      new_ones[0] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[0], 0.25);
      new_ones[1] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, leaf.lines[0], leaf.lines[1], true);
      new_ones[2] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[1], 0.5);
      new_ones[3] = lane_riesenfeld_regular(solver, triangles, positions,
          adjacencies, leaf.lines[1], leaf.lines[2]);
      new_ones[4] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[2], 0.5);
    } else
      assert(false);
  } else {
    if (new_ones_entries.back() == yocto::pow(2, k) + 2) {
      new_ones[0] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[0], 0.5);
      new_ones[1] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, leaf.lines[0], leaf.lines[1], false);
      new_ones[2] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[1], 0.75);
      new_ones[3] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[2], 0.5);
      new_ones[4] = leaf.points[3];
    } else if (new_ones_entries.back() == yocto::pow(2, k)) {
      new_ones[0] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[0], 0.5);
      new_ones[1] = lane_riesenfeld_regular(solver, triangles, positions,
          adjacencies, leaf.lines[0], leaf.lines[1]);
      new_ones[2] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[1], 0.5);
      new_ones[3] = lane_riesenfeld_boundary(solver, triangles, positions,
          adjacencies, leaf.lines[1], leaf.lines[2], false);
      new_ones[4] = eval_path_point(
          solver, triangles, positions, adjacencies, leaf.lines[2], 0.75);
    } else
      assert(false);
  }
  return {true, new_ones};
}

static pair<spline_node, spline_node> lane_riesenfeld_split_node(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const spline_node& leaf, bool& max_depth_reached) {
  auto curr_t     = (leaf.t.x + leaf.t.y) / 2;
  int  curr_entry = (int)(yocto::pow(2, leaf.depth + 1) * curr_t);
  curr_entry += 3;
  if (curr_entry % 2) {
    max_depth_reached = true;
    return {leaf, {}};
  }
  auto curr_entries               = vector<int>{curr_entry - 4, curr_entry - 3,
      curr_entry - 2, curr_entry - 1, curr_entry};
  auto [are_boundaries, new_ones] = handle_boundary_node(
      solver, triangles, positions, adjacencies, leaf, curr_entries);
  if (!are_boundaries) {
    new_ones[0] = eval_path_point(
        solver, triangles, positions, adjacencies, leaf.lines[0], 0.5);
    new_ones[1] = lane_riesenfeld_regular(solver, triangles, positions,
        adjacencies, leaf.lines[0], leaf.lines[1]);
    new_ones[2] = eval_path_point(
        solver, triangles, positions, adjacencies, leaf.lines[1], 0.5);
    new_ones[3] = lane_riesenfeld_regular(solver, triangles, positions,
        adjacencies, leaf.lines[1], leaf.lines[2]);
    new_ones[4] = eval_path_point(
        solver, triangles, positions, adjacencies, leaf.lines[2], 0.5);
  }
  auto L01 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, new_ones[0], new_ones[1]);
  auto L12 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, new_ones[1], new_ones[2]);
  auto L23 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, new_ones[2], new_ones[3]);
  spline_node P0 = {{new_ones[0], new_ones[1], new_ones[2], new_ones[3]},
      {L01, L12, L23}, {leaf.t.x, curr_t}, leaf.depth + 1};
  L01            = compute_geodesic_path(
      solver, triangles, positions, adjacencies, new_ones[3], new_ones[4]);
  spline_node P1 = {{new_ones[1], new_ones[2], new_ones[3], new_ones[4]},
      {L12, L23, L01}, {curr_t, leaf.t.y}, leaf.depth + 1};
  return {P0, P1};
}

static vector<mesh_point> lane_riesenfeld_uniform(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 4>& control_points, int num_subdivisions) {
  auto size = 7;
  struct parametric_path {
    geodesic_path path = {};
    vector<float> t    = {};
  };
  parametric_path curr_path = {};
  parametric_path gamma01   = {};
  parametric_path gamma32   = {};
  auto            prev      = mesh_point{};
  auto            curr      = mesh_point{};
  auto            q         = vector<mesh_point>(size);

  {
    auto& p      = control_points;
    gamma01.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[0], p[1]);
    gamma01.t = path_parameters(
        gamma01.path, triangles, positions, adjacencies);
    gamma32.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[3], p[2]);
    gamma32.t = path_parameters(
        gamma32.path, triangles, positions, adjacencies);
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[1], p[2]);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);
    q[0]           = p[0];
    q[1]           = eval_path_point(solver, triangles, positions, adjacencies,
        gamma01.path, gamma01.t, 0.25);
    auto p0p1      = eval_path_point(solver, triangles, positions, adjacencies,
        gamma01.path, gamma01.t, 0.5);
    auto p1p2      = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.5);
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p0p1, p1p2);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);
    q[2]           = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.25);
    auto p2p3      = eval_path_point(solver, triangles, positions, adjacencies,
        gamma32.path, gamma32.t, 0.5);
    prev           = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p1p2, p2p3);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);
    curr = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 3 / 8.f);
    q[3] = geodesic_lerp(
        solver, triangles, positions, adjacencies, prev, curr, 0.5);
    q[4] = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.75);
    q[5] = eval_path_point(solver, triangles, positions, adjacencies,
        gamma32.path, gamma32.t, 0.25);
    q[6] = p[3];
  }

  auto p = vector<mesh_point>{};

  for (int subdiv = 0; subdiv < num_subdivisions; subdiv++) {
    std::swap(p, q);

    auto new_size = 2 * size - 3;
    q.resize(new_size);
    q[0]           = p[0];
    q[1]           = eval_path_point(solver, triangles, positions, adjacencies,
        gamma01.path, gamma01.t, 1.f / yocto::pow(2, 3 + subdiv));
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[1], p[2]);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);
    q[2]           = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.25);
    prev           = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 5 / 8.f);
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[2], p[3]);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);
    curr = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.25);
    q[3] = geodesic_lerp(
        solver, triangles, positions, adjacencies, prev, curr, 0.5);
    prev = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.5);
    for (int j = 4; j < 2 * size - 8; j += 2) {
      q[j] = prev;
      prev = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(solver, triangles, positions,
          adjacencies, p[j / 2 + 1], p[j / 2 + 2]);
      curr_path.t    = path_parameters(
          curr_path.path, triangles, positions, adjacencies);
      curr     = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.25);
      q[j + 1] = geodesic_lerp(
          solver, triangles, positions, adjacencies, prev, curr, 1 / 2.f);
      prev = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.5);
    }
    q[2 * size - 8] = prev;
    {
      auto qq = &q[new_size - 4];
      auto pp = &p[size - 4];
      prev    = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(
          solver, triangles, positions, adjacencies, pp[1], pp[2]);
      curr_path.t = path_parameters(
          curr_path.path, triangles, positions, adjacencies);
      curr  = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 3 / 8.f);
      qq[0] = geodesic_lerp(
          solver, triangles, positions, adjacencies, prev, curr, 0.5);
      qq[1] = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.75);
      qq[2] = eval_path_point(solver, triangles, positions, adjacencies,
          gamma32.path, gamma32.t, 1.f / yocto::pow(2, 3 + subdiv));
      qq[3] = pp[3];
    }
    size = new_size;
  }
  return q;
}

[[maybe_unused]] static vector<mesh_point> spline_subdivision_uniform_quadric(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 3>& control_points, int num_subdivisions) {
  auto size = 6;
  struct parametric_path {
    geodesic_path path = {};
    vector<float> t    = {};
  };
  parametric_path curr_path = {};
  parametric_path gamma01   = {};
  parametric_path gamma21   = {};

  auto q = vector<mesh_point>(size);

  {
    auto& p      = control_points;
    gamma01.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[0], p[1]);
    gamma01.t = path_parameters(
        gamma01.path, triangles, positions, adjacencies);
    gamma21.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[2], p[1]);
    gamma21.t = path_parameters(
        gamma21.path, triangles, positions, adjacencies);
    q[0] = p[0];
    q[1] = eval_path_point(solver, triangles, positions, adjacencies,
        gamma01.path, gamma01.t, 0.25);
    q[4] = eval_path_point(solver, triangles, positions, adjacencies,
        gamma21.path, gamma21.t, 0.25);
    q[5] = p[2];

    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, q[2], q[5]);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);

    q[2] = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.25);
    q[3] = eval_path_point(solver, triangles, positions, adjacencies,
        curr_path.path, curr_path.t, 0.75);
  }

  auto p         = vector<mesh_point>{};
  curr_path.path = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[1], q[2]);
  curr_path.t = path_parameters(
      curr_path.path, triangles, positions, adjacencies);
  for (int subdiv = 0; subdiv < num_subdivisions; subdiv++) {
    std::swap(p, q);

    auto new_size = 2 * size - 2;
    q.resize(new_size);
    q[0]           = p[0];
    q[1]           = eval_path_point(solver, triangles, positions, adjacencies,
        gamma01.path, gamma01.t, 1.f / yocto::pow(2, 3 + subdiv));
    curr_path.path = compute_geodesic_path(
        solver, triangles, positions, adjacencies, p[1], p[2]);
    curr_path.t = path_parameters(
        curr_path.path, triangles, positions, adjacencies);

    for (int j = 2; j < 2 * size - 4; j += 2) {
      q[j]     = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.25);
      q[j + 1] = eval_path_point(solver, triangles, positions, adjacencies,
          curr_path.path, curr_path.t, 0.75);
      curr_path.path = compute_geodesic_path(solver, triangles, positions,
          adjacencies, p[j / 2 + 1], p[j / 2 + 2]);
      curr_path.t    = path_parameters(
          curr_path.path, triangles, positions, adjacencies);
    }

    q[2 * size - 4] = eval_path_point(solver, triangles, positions, adjacencies,
        gamma21.path, gamma21.t, 1.f / yocto::pow(2, 3 + subdiv));
    q[2 * size - 3] = p.back();
    size            = new_size;
  }
  return q;
}

static vector<mesh_point> lane_riesenfeld_adaptive(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 4>& polygon, const spline_params& params) {
  auto q = vector<mesh_point>(7);

  auto& p = polygon;
  q[0]    = p[0];
  q[1]    = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[0], p[1], 1 / 4.f);
  auto p0p1 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[0], p[1], 1 / 2.f);
  auto p1p2 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[1], p[2], 1 / 2.f);
  q[2] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p0p1, p1p2, 1 / 4.f);
  auto p2p3 = geodesic_lerp(
      solver, triangles, positions, adjacencies, p[2], p[3], 1 / 2.f);
  q[3] = lane_riesenfeld_init(
      solver, triangles, positions, adjacencies, p0p1, p1p2, p2p3);
  q[4] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p1p2, p2p3, 3 / 4.f);
  q[5] = geodesic_lerp(
      solver, triangles, positions, adjacencies, p2p3, p[3], 1 / 2.f);
  q[6]     = p[3];
  auto L01 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[0], q[1]);
  auto L12 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[1], q[2]);
  auto L23 = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[2], q[3]);
  spline_node P0 = {{q[0], q[1], q[2], q[3]}, {L01, L12, L23}, {0, 0.25}, 2};
  L01            = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[3], q[4]);
  spline_node P1 = {{q[1], q[2], q[3], q[4]}, {L12, L23, L01}, {0.25, 0.5}, 2};
  L12            = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[4], q[5]);
  spline_node P2 = {{q[2], q[3], q[4], q[5]}, {L23, L01, L12}, {0.5, 0.75}, 2};
  L23            = compute_geodesic_path(
      solver, triangles, positions, adjacencies, q[5], q[6]);
  spline_node P3 = {{q[3], q[4], q[5], q[6]}, {L01, L12, L23}, {0.75, 1}, 2};

  P0.is_good = is_bezier_straight_enough(P0.lines[0], P0.lines[1], P0.lines[2],
      solver, triangles, positions, adjacencies, params);
  P1.is_good = is_bezier_straight_enough(P1.lines[0], P1.lines[1], P1.lines[2],
      solver, triangles, positions, adjacencies, params);
  P2.is_good = is_bezier_straight_enough(P2.lines[0], P2.lines[1], P2.lines[2],
      solver, triangles, positions, adjacencies, params);
  P3.is_good = is_bezier_straight_enough(P3.lines[0], P3.lines[1], P3.lines[2],
      solver, triangles, positions, adjacencies, params);
  // auto                    count_path = 16;
  // auto                    count_eval = 8;
  std::deque<spline_node> Q;
  Q.push_back(P3);
  Q.push_back(P2);
  Q.push_back(P1);
  Q.push_back(P0);
  auto P                 = vector<spline_node>{};
  bool max_depth_reached = false;
  while (!Q.empty()) {
    auto curr = Q.back();
    Q.pop_back();
    if (P.size() > 0 && P.back().depth == curr.depth && curr.is_good) {
      P.push_back(curr);

    } else {
      auto [left, right] = lane_riesenfeld_split_node(
          solver, triangles, positions, adjacencies, curr, max_depth_reached);
      if (max_depth_reached) {
        curr.is_good = true;
        if (P.size() > 0) {
          auto last = P.back();
          auto L    = compute_geodesic_path(solver, triangles, positions,
              adjacencies, last.points[0], curr.points[0]);
          if (is_bezier_straight_enough(L, curr.lines[0], curr.lines[1], solver,
                  triangles, positions, adjacencies, params))
            P.push_back(curr);
          else {
            P.pop_back();
            last = {{last.points[0], curr.points[0], curr.points[1],
                        curr.points[2]},
                {L, curr.lines[0], curr.lines[1]}, last.t, last.depth, false};
            Q.push_back(curr);
            Q.push_back(last);
          }
        } else
          P.push_back(curr);

        max_depth_reached = false;
      } else {
        left.is_good  = is_bezier_straight_enough(left.lines[0], left.lines[1],
            left.lines[2], solver, triangles, positions, adjacencies, params);
        right.is_good = is_bezier_straight_enough(right.lines[0],
            right.lines[1], right.lines[2], solver, triangles, positions,
            adjacencies, params);
        if (left.is_good && right.is_good) {
          if (P.size() == 0) {
            P.push_back(left);
            P.push_back(right);
          } else {
            auto last = P.back();
            auto L    = compute_geodesic_path(solver, triangles, positions,
                adjacencies, last.points[0], left.points[0]);
            if (is_bezier_straight_enough(L, left.lines[0], left.lines[1],
                    solver, triangles, positions, adjacencies, params)) {
              last     = {{last.points[0], left.points[0], left.points[1],
                          left.points[2]},
                  {L, left.lines[0], left.lines[1]}, last.t, last.depth, true};
              P.back() = last;
              P.push_back(left);
              P.push_back(right);
            } else if (left.depth < last.depth) {
              left.is_good  = false;
              right.is_good = false;
              Q.push_back(right);
              Q.push_back(left);
            } else {
              P.pop_back();
              last = {{last.points[0], left.points[0], left.points[1],
                          left.points[2]},
                  {L, left.lines[0], left.lines[1]}, last.t, last.depth, false};
              Q.push_back(right);
              Q.push_back(left);
              Q.push_back(last);
            }
          }
        } else {
          Q.push_back(right);
          Q.push_back(left);
        }
      }
    }
  }
  // std::cout << "Number of paths";
  // std::cout << count_path << std::endl;
  // std::cout << "Number of eval";
  // std::cout << count_eval << std::endl;
  // std::cout << "Precision";
  // std::cout << yocto::pow(2, -params.precision) / pif * 180 << std::endl;
  auto polyline = vector<mesh_point>{};
  for (auto i = 0; i < P.size(); ++i) {
    if (i == 0) {
      for (auto j = 0; j < 4; ++j) {
        polyline.push_back(P[i].points[j]);
      }
    } else
      polyline.push_back(P[i].points.back());
  }

  return polyline;
}

[[maybe_unused]] static vector<mesh_point> degree_elevation(
    const dual_geodesic_solver& solver, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const array<mesh_point, 4>& control_points, int num_subdivisions) {
  auto p = vector<mesh_point>{control_points[0], control_points[1],
      control_points[2], control_points[3]};
  auto q = p;

  for (auto i = 0; i < num_subdivisions; ++i) {
    std::swap(p, q);
    auto new_size = p.size() + 1;
    q.resize(new_size);
    q[0] = p[0];
    for (auto j = 1; j < new_size - 1; ++j) {
      float alpha = (new_size - 1 - j) / (float)(new_size - 1);
      q[j]        = geodesic_lerp(
          solver, triangles, positions, adjacencies, p[j - 1], p[j], alpha);
    }
    q.back() = p.back();
  }

  return q;
}

}  // namespace yocto
