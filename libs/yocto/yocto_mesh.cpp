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
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>

#include "yocto_commonio.h"
#include "yocto_geometry.h"
#include "yocto_modelio.h"

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// TODO: cleanup
pair<bool, int> bary_is_vert(const vec3f& bary, float tol) {
  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

// TODO: cleanup
pair<bool, int> bary_is_edge(const vec3f& bary, float tol) {
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

// TODO: cleanup
pair<bool, int> point_is_vert(const mesh_point& p, float tol) {
  auto bary = vec3f{1 - p.uv.x - p.uv.y, p.uv.x, p.uv.y};
  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

// TODO: cleanup
pair<bool, int> point_is_edge(const mesh_point& p, float tol) {
  auto bary = vec3f{1 - p.uv.x - p.uv.y, p.uv.x, p.uv.y};
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol) return {true, 0};
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol) return {true, 1};
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol) return {true, 2};
  return {false, -1};
}

// TODO: cleanup
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

// TODO: cleanup
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

inline int find_adjacent_triangle(
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

inline int find_adjacent_triangle(
    const vector<vec3i>& triangles, int face, int neighbor) {
  return find_adjacent_triangle(triangles[face], triangles[neighbor]);
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
  auto v  = triangles[neighbor][mod3(j + 2)];
  auto a  = triangles[face][k];
  auto b  = triangles[face][mod3(k + 1)];
  auto r0 = length_squared(positions[v] - positions[a]);
  auto r1 = length_squared(positions[v] - positions[b]);

  auto res         = unfold_triangle{};
  res[j]           = tr[mod3(k + 1)];
  res[mod3(j + 1)] = tr[k];
  res[mod3(j + 2)] = intersect_circles(res[j], r1, res[mod3(j + 1)], r0);
  return res;
}

// TODO: cleanup
static unfold_triangle unfold_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const unfold_triangle& tr, int face, int k) {
  return unfold_face(triangles, positions, tr, face, adjacencies[face][k]);
}

// assign 2D coordinates to vertices of the triangle containing the mesh
// point, putting the point at (0, 0)
unfold_triangle triangle_coordinates(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& point) {
  auto first = unfold_triangle{};
  auto tr    = triangles[point.face];
  first[0]   = {0, 0};
  first[1]   = {0, length(positions[tr.x] - positions[tr.y])};
  auto rx    = length_squared(positions[tr.x] - positions[tr.z]);
  auto ry    = length_squared(positions[tr.y] - positions[tr.z]);
  first[2]   = intersect_circles(first[0], rx, first[1], ry);

  // Transform coordinates such that point = (0, 0)
  auto point_coords = interpolate_triangle(
      first[0], first[1], first[2], point.uv);
  first[0] -= point_coords;
  first[1] -= point_coords;
  first[2] -= point_coords;
  return first;
}

// assign 2D coordinates to a strip of triangles. point start is at (0, 0)
vector<unfold_triangle> unfold_strip(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<int>& strip,
    const mesh_point& start) {
  auto coords = vector<unfold_triangle>(strip.size());
  assert(start.face == strip[0]);
  coords[0] = triangle_coordinates(triangles, positions, start);

  for (auto i = 1; i < strip.size(); i++) {
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

inline vector<pair<vec2f, vec2f>> unfold_funnel_portals(
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

// TODO: cleanup
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
        solver.graph[i][k].length = length(
            get_triangle_center(triangles, positions, i) -
            get_triangle_center(triangles, positions, adjacencies[i][k]));
      }
    }
  }
  return solver;
}

// TODO: cleanup
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
vector<int> compute_geodesic_paths(
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
  // @Speed: use parallel_for
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

// TODO: cleanup
surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from,
    int to) {
  return integral_paths::integrate_field(
      triangles, positions, adjacency, tags, tag, field, from, to);
}

// TODO: cleanup
vector<vec3f> make_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions) {
  return integral_paths::make_positions_from_path(path, mesh_positions);
}

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
// utilities:nodes_around_point-->length_by_flattening
// returns the shortest path starting from target to source as a list of
// indices of nodes of the graph(solver). note: the list does not contains
// source if it is a vertex
vector<int> point_to_point_geodesic_path(const geodesic_solver& solver,
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

// TODO: cleanup
// same function of "compute_geodesic_paths" of yocto_mesh.cpp that makes
// early exit when reaching end_vertex, the name is changed because the input
// parameters are the same
vector<int> compute_pruned_geodesic_paths(
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

// TODO: cleanup
// returns the shortest path starting from target to source as a list of
// indices of nodes of the graph(solver). note: the list does not contains
// source and target.
vector<int> point_to_point_geodesic_path(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int source, int target) {
  vector<int> sources = {source};
  vector<int> parents = compute_pruned_geodesic_paths(solver, sources, target);
  vector<int> path    = {};
  set_source_child(parents, target, path);
  path.pop_back();  // we remove source from the list
  return path;
}

// TODO: cleanup
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

// TODO: cleanup
// TODO(fabio): better name
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

// TODO: cleanup
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
  // assert(check_strip(mesh.adjacencies, strip));
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

// TODO: cleanup
// static bool set_ord(float theta_next, float theta_prev) {
//   if (theta_next > theta_prev) {
//     return theta_next - theta_prev < pif;
//   } else {
//     return theta_prev - theta_next > pif;
//   }
// }

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
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

// TODO: cleanup
// returns a strip of triangles such target belongs to the first one and
// source to the last one
// TODO(fabio_): may be the names could change in order to get the call
// more consistent with the output)
vector<int> get_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  if (target.face == source.face) return {target.face};
  auto parents = point_to_point_geodesic_path(
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

inline float intersect_segments(const vec2f& start1, const vec2f& end1,
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

  auto result = vector<int>{};
  result.reserve(strip.size() + 12);
  // Initial part of original strip, until intersection with fan.
  for (auto i = 0; i < first_strip_intersection; ++i)
    result.push_back(strip[i]);

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
  for (auto i = 0; i < path.strip.size() * 2 && index != -1; i++) {
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
    const vector<vec3i>& adjacencies, float t) {
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

  auto parameter_t = path_parameters(
      path_positions(path, triangles, positions, adjacencies));

  // find the point in the middle
  auto i = 0;
  for (; i < parameter_t.size() - 1; i++) {
    if (parameter_t[i + 1] >= 0.5) break;
  }
  auto t_low = parameter_t[i], t_high = parameter_t[i + 1];
  auto alpha = (t - t_low) / (t_high - t_low);
  // alpha == 0 -> t_low, alpha == 1 -> t_high
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MESH IO
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
static vector<vec3i> quads_to_triangles(const vector<vec4i>& quads) {
  auto triangles = vector<vec3i>{};
  triangles.reserve(quads.size() * 2);
  for (auto& q : quads) {
    triangles.push_back({q.x, q.y, q.w});
    if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
  }
  return triangles;
}

// Load ply mesh
bool load_mesh(const string& filename, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec4f>& colors, string& error, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  triangles = {};
  positions = {};
  normals   = {};
  texcoords = {};
  colors    = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // open ply
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    // gets vertex
    get_positions(ply, positions);
    get_normals(ply, normals);
    get_texcoords(ply, texcoords, flip_texcoord);
    get_colors(ply, colors);
    // get faces
    get_triangles(ply, triangles);
    if (triangles.empty()) return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    // load obj
    auto obj = obj_shape{};
    if (!load_obj(filename, obj, error, true)) return false;
    // decide what to do and get properties
    auto materials = vector<int>{};
    get_positions(obj, positions);
    get_normals(obj, normals);
    get_texcoords(obj, texcoords, flip_texcoord);
    get_triangles(obj, triangles, materials);
    if (triangles.empty()) return shape_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_mesh(const string& filename, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<vec4f>& colors, string& error,
    bool ascii, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // create ply
    auto ply = ply_model{};
    add_positions(ply, positions);
    add_normals(ply, normals);
    add_texcoords(ply, texcoords, flip_texcoord);
    add_colors(ply, colors);
    add_triangles(ply, triangles);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, positions);
    add_normals(obj, normals);
    add_texcoords(obj, texcoords, flip_texcoord);
    add_triangles(obj, triangles, 0, !normals.empty(), !texcoords.empty());
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (triangles.empty()) return shape_error();
    add_triangles(stl, triangles, positions, {});
    if (!save_stl(filename, stl, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

// Load ply mesh
bool load_mesh(const string& filename, vector<vec3i>& triangles,
    vector<vec3f>& positions, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  triangles = {};
  positions = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // open ply
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    get_positions(ply, positions);
    get_triangles(ply, triangles);
    if (positions.empty()) return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    // load obj
    auto obj = obj_shape{};
    if (!load_obj(filename, obj, error, true)) return false;
    auto materials = vector<int>{};
    get_positions(obj, positions);
    get_triangles(obj, triangles, materials);
    if (triangles.empty()) return shape_error();
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    // open ply
    auto stl = stl_model{};
    if (!load_stl(filename, stl, error)) return false;
    if (stl.shapes.empty()) return shape_error();
    if (stl.shapes.size() > 1) return shape_error();
    auto fnormals = vector<vec3f>{};
    if (!get_triangles(stl, 0, triangles, positions, fnormals))
      return shape_error();
    if (positions.empty()) return shape_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_mesh(const string& filename, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, string& error, bool ascii) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // create ply
    auto ply = ply_model{};
    if (triangles.empty()) return shape_error();
    add_positions(ply, positions);
    add_triangles(ply, triangles);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, positions);
    add_triangles(obj, triangles, 0, false, false);
    auto err = ""s;
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else if (ext == ".stl" || ext == ".STL") {
    auto stl = stl_model{};
    if (triangles.empty()) return shape_error();
    add_triangles(stl, triangles, positions, {});
    if (!save_stl(filename, stl, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

// Load ply mesh
bool load_lines(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec4f>& colors, string& error, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  lines     = {};
  positions = {};
  normals   = {};
  texcoords = {};
  colors    = {};

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // open ply
    auto ply = ply_model{};
    if (!load_ply(filename, ply, error)) return false;
    // gets vertex
    get_positions(ply, positions);
    get_normals(ply, normals);
    get_texcoords(ply, texcoords, flip_texcoord);
    get_colors(ply, colors);
    // get faces
    get_lines(ply, lines);
    if (positions.empty()) return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    // load obj
    auto obj = obj_shape{};
    if (!load_obj(filename, obj, error, true)) return false;
    auto materials = vector<int>{};
    get_positions(obj, positions);
    get_normals(obj, normals);
    get_texcoords(obj, texcoords, flip_texcoord);
    get_lines(obj, lines, materials);
    if (lines.empty()) return shape_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
bool save_lines(const string& filename, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<vec4f>& colors, string& error,
    bool ascii, bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // create ply
    auto ply = ply_model{};
    add_positions(ply, positions);
    add_normals(ply, normals);
    add_texcoords(ply, texcoords, flip_texcoord);
    add_colors(ply, colors);
    add_lines(ply, lines);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj = obj_shape{};
    add_positions(obj, positions);
    add_normals(obj, normals);
    add_texcoords(obj, texcoords, flip_texcoord);
    add_lines(obj, lines, 0, !normals.empty(), !texcoords.empty());
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> mesh_stats(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<vec4f>& colors, bool verbose) {
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = invalidb3f;
  for (auto& pos : positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("triangles:    " + format(triangles.size()));
  stats.push_back("positions:    " + format(positions.size()));
  stats.push_back("normals:      " + format(normals.size()));
  stats.push_back("texcoords:    " + format(texcoords.size()));
  stats.push_back("colors:       " + format(colors.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto
