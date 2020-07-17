//
// Implementation for Yocto/Mesh
//

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_mesh.h"

#include <atomic>
#include <cassert>
#include <deque>
#include <memory>
#include <string>
#include <utility>

#include "yocto_geometry.h"
#include "yocto_modelio.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::deque;
using std::pair;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// find a value in a vector or vecs
template <typename T>
inline int find_in_vec(const T& vec, int x) {
  for (int i = 0; i < size(vec); i++)
    if (vec[i] == x) return i;
  return -1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STRIPS
// -----------------------------------------------------------------------------
namespace yocto {

using unfold_triangle = std::array<vec2f, 3>;

inline vec2f intersect_circles(
    const vec2f& c2, float R2, const vec2f& c1, float R1) {
  auto R = length_squared(c2 - c1);
  assert(R > 0);
  auto result = (c1 + c2) / 2;
  result += (c2 - c1) * ((R1 - R2) / (2 * R));
  auto A = 2 * (R1 + R2) / R;
  auto B = (R1 - R2) / R;
  auto s = A - B * B - 1;
  assert(s >= 0);
  result += vec2f{c2.y - c1.y, c1.x - c2.x} * (0.5f * sqrt(s));
  return result;
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

inline unfold_triangle unfold_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const unfold_triangle& tr, int face, int k) {
  auto v = opposite_vertex(triangles, adjacencies, face, k);
  assert(v != -1);

  auto a  = triangles[face][k];
  auto b  = triangles[face][(k + 1) % 3];
  auto r0 = length_squared(positions[v] - positions[a]);
  auto r1 = length_squared(positions[v] - positions[b]);

  auto neighbor = adjacencies[face][k];
  auto j        = find_in_vec(adjacencies[neighbor], face);
  assert(j != -1);

  auto res         = unfold_triangle{};
  res[j]           = tr[(k + 1) % 3];
  res[(j + 1) % 3] = tr[k];
  res[(j + 2) % 3] = intersect_circles(res[j], r1, res[(j + 1) % 3], r0);
  return res;
}

// Utilities
inline vec3f make_bary(const vec2f& bary) {
  return vec3f{1 - bary.x - bary.y, bary.x, bary.y};
}

inline mesh_point make_point(int tid, const vec3f& bary) {
  return {tid, vec2f{bary.y, bary.z}};
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
    const mesh_point& p, vector<int>& strip, vec2f& first_sample_direction) {
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

  auto last   = coords.back();
  auto bary   = make_bary(p.uv);
  auto pos2d  = last[0] * bary.x + last[1] * bary.y + last[2] * bary.z;
  auto v      = pos2d - coords[0][(h + 2) % 3];
  auto w0     = coords[0][h] - coords[0][(h + 2) % 3];
  auto w1     = coords[0][(h + 1) % 3] - coords[0][(h + 2) % 3];
  auto phi    = angle(w0, w1);
  auto theta0 = angle(v, w0);
  auto theta1 = angle(v, w1);

  if (theta0 < phi && theta1 < phi) {
    first_sample_direction = pos2d;
    return length(v);
  } else if (theta1 > phi) {
    first_sample_direction = coords[0][h];
    float len              = length(w0);
    len += length(coords[0][h] - pos2d);
    return len;
  } else {
    first_sample_direction = coords[0][(h + 1) % 3];
    float len              = length(w1);
    len += length(coords[0][(h + 1) % 3] - pos2d);
    return len;
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
  auto offset = 2 - int(clockwise);
  while (true) {
    if (node == -1) break;
    if (node == face) break;
    result.push_back(node);
    auto kk = find_in_vec(adjacencies[node], prev);
    assert(kk != -1);
    kk   = (kk + offset) % 3;
    prev = node;
    node = adjacencies[node][kk];
  }
  return result;
}

// returns the list of triangles incident at each vertex in CCW order
vector<vector<int>> vertex_to_triangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies) {
  auto v2t    = vector<vector<int>>{positions.size(), vector<int>{}};
  auto offset = 0;
  for (auto i = 0; i < (int)triangles.size(); ++i) {
    for (auto j = 0; j < 3; ++j) {
      auto curr = triangles[i][j];
      if (v2t[curr].size() == 0) {
        offset    = find_in_vec(triangles[i], curr);
        v2t[curr] = triangle_fan(adjacencies, i, (offset + 2) % 3);
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
  return tt[(j + 2) % 3];
}

// Finds common edge between triangles
vec2i common_edge(const vec3i& triangle0, const vec3i& triangle1) {
  for (auto i = 0; i < 3; i++) {
    for (auto k = 0; k < 3; k++) {
      if (triangle0[i] == triangle1[k] &&
          triangle0[(i + 1) % 3] == triangle1[(k + 2) % 3])
        return {triangle0[i], triangle0[(i + 1) % 3]};
    }
  }
  return {-1, -1};
}

vec2i opposite_edge(const vec3i& t, int vid) {
  auto offset = find_in_vec(t, vid);
  auto v0     = t[(offset + 1) % 3];
  auto v1     = t[(offset + 2) % 3];
  return vec2i{v0, v1};
}

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

bool point_in_triangle(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int tid, const vec3f& point, vec2f& uv,
    float tol) {
  // http://www.r-5.org/files/books/computers/algo-list/realtime-3d/Christer_Ericson-Real-Time_Collision_Detection-EN.pdf
  // pag.48
  auto b  = vec3f{0.0, 0.0, 0.0};
  auto v0 = positions[triangles[tid].x];
  auto v1 = positions[triangles[tid].y];
  auto v2 = positions[triangles[tid].z];

  auto u = v1 - v0, v = v2 - v0, w = point - v0;
  auto d00 = dot(u, u), d01 = dot(u, v), d11 = dot(v, v), d20 = dot(w, u),
       d21 = dot(w, v), d = d00 * d11 - d01 * d01;

  if (d == 0) return false;

  b[2] = (d00 * d21 - d01 * d20) / d;
  assert(!isnan(b[2]));
  b[1] = (d11 * d20 - d01 * d21) / d;
  assert(!isnan(b[1]));
  b[0] = 1.0 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (int i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol) return false;
  }

  uv = vec2f{b.y, b.z};
  uv = clamp(uv, 0.f, 1.f);

  return true;
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
  b[0] = 1.0 - b[1] - b[2];
  assert(!isnan(b[0]));

  for (auto i = 0; i < 3; ++i) {
    if (b[i] < -tol || b[i] > 1.0 + tol) return {false, zero2f};
  }
  return {true, vec2f{b.y, b.z}};
}

// compute angles in tangent space, if opposite is false we do not consider
// opposite vertices
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

        auto phi    = angle(w0, w1);
        auto theta0 = angle(v, w0);
        auto theta1 = angle(v, w1);

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
        auto k   = find_in_vec(triangles[tid], i);
        auto v   = positions[tr3d[(k + 1) % 3]] - positions[tr3d[k]];
        auto w   = positions[tr3d[(k + 2) % 3]] - positions[tr3d[k]];
        auto phi = angle(v, w);
        theta += phi;
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
      if (get_tag(tr[k]) != get_tag(tr[(k + 1) % 3]) &&
          get_tag(tr[k]) != get_tag(tr[(k + 2) % 3])) {
        j = k;
      }
    }

    // If all vertices have same tag, then tag the whole triangle and continue.
    if (j == -1) {
      tags[i] = get_tag(tr.x);
      continue;
    }

    // Reorder the triangle such that vertex with different tag is always z
    tr = {tr[(j + 2) % 3], tr[(j + 1) % 3], tr[j]};

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

static inline void connect_nodes(
    geodesic_solver& solver, int a, int b, float length) {
  solver.graph[a].push_back({b, length});
  solver.graph[b].push_back({a, length});
}

static inline float opposite_nodes_arc_length(
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

static inline void connect_opposite_nodes(geodesic_solver& solver,
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
      auto b = triangles[face][(k + 1) % 3];

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
      auto k = find_in_vec(adjacencies[tid], opp);  // TODO: this is not needeed
      assert(k != -1);
      auto a       = opposite_vertex(triangles, adjacencies, tid, k);
      offset       = find_in_vec(triangles[opp], a);
      auto bary    = zero3f;
      bary[offset] = 1;
      auto pos     = zero2f;
      auto l       = length_by_flattening(
          triangles, positions, adjacencies, make_point(opp, bary), strip, pos);
      solver.graph[i].push_back({a, l});
    }
  }
  return solver;
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
  double cumulative_weight = 0.0;

  // setup queue
  auto queue = deque<int>();
  for (auto source : sources) {
    in_queue[source] = true;
    cumulative_weight += field[source];
    queue.push_back(source);
  }

  while (!queue.empty()) {
    auto node           = queue.front();
    auto average_weight = (float)cumulative_weight / queue.size();

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
  };
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF INTEGRAL PATHS
// -----------------------------------------------------------------------------
namespace yocto {

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
      auto edge = vec2i{triangles[f][k], triangles[f][(k + 1) % 3]};
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
      if (tr[i] == vertex) return {tr[(i + 1) % 3], tr[(i + 2) % 3]};
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
      auto to = triangles[face][(k + 1) % 3];

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
        auto edge = vec2i{triangles[face][k], triangles[face][(k + 1) % 3]};
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

vector<vec3f> make_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions) {
  return integral_paths::make_positions_from_path(path, mesh_positions);
}

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

inline bool is_vert(const mesh_point& p, int& offset, float tol = 1e-2) {
  auto bary = make_bary(p.uv);
  if (bary[0] > tol && bary[1] <= tol && bary[2] <= tol) {
    offset = 0;
    return true;
  }
  if (bary[1] > tol && bary[0] <= tol && bary[2] <= tol) {
    offset = 1;
    return true;
  }
  if (bary[2] > tol && bary[0] <= tol && bary[1] <= tol) {
    offset = 2;
    return true;
  }
  return false;
}

inline bool is_edge(const mesh_point& p, int& offset, float tol = 5 * 1e-3) {
  auto bary = make_bary(p.uv);
  if (bary[0] > tol && bary[1] > tol && bary[2] <= tol) {
    offset = 0;
    return true;
  }
  if (bary[1] > tol && bary[2] > tol && bary[0] <= tol) {
    offset = 1;
    return true;
  }
  if (bary[2] > tol && bary[0] > tol && bary[1] <= tol) {
    offset = 2;
    return true;
  }
  return false;
}

inline vec3f eval_position(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& sample) {
  auto [x, y, z] = triangles[sample.face];
  return interpolate_triangle(
      positions[x], positions[y], positions[z], sample.uv);
}

inline vec3f eval_normal(const vector<vec3i>& triangles,
    const vector<vec3f>& normals, const mesh_point& point) {
  auto [x, y, z] = triangles[point.face];
  return normalize(
      interpolate_triangle(normals[x], normals[y], normals[z], point.uv));
}

inline int node_is_neighboor(const geodesic_solver& solver, int vid, int node) {
  auto nbr = solver.graph[vid];
  for (auto i = 0; i < nbr.size(); ++i) {
    if (nbr[i].node == node) {
      return i;
    }
  }
  return -1;
}
// vid is the parent, p is the point and k is adjacency coeffient of the face
// containing p with the one in the star of vid, if k==-1, then vid is a vertex
// of p.face
static float get_rescaled_angle(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<float>>& angles,
    int tid, int vid, const mesh_point& p, int k, int& entry) {
  if (k == -1) {
    auto h    = find_in_vec(triangles[p.face], vid);
    auto pos  = eval_position(triangles, positions, p);
    auto v    = normalize(pos - positions[vid]);
    auto vid0 = triangles[p.face][(h + 1) % 3];
    auto vid1 = triangles[p.face][(h + 2) % 3];
    entry     = node_is_neighboor(solver, vid, vid0);

    auto s = angles[vid].size();

    auto v0     = positions[vid0] - positions[vid];
    auto v1     = positions[vid1] - positions[vid];
    auto theta0 = angles[vid][entry];
    auto phi3D  = angle(v0, v1);

    auto phi2D = (entry + 2 == s) ? 2 * pif - theta0
                                  : angles[vid][entry + 2] - theta0;
    auto scale_factor = phi3D / phi2D;
    auto theta        = angle(v0, v);
    ++entry;
    return theta0 + theta * scale_factor;
  } else {
    auto tr2d      = init_flat_triangle(positions, triangles[tid]);
    auto flat_face = unfold_face(
        triangles, positions, adjacencies, tr2d, tid, k);
    auto bary = make_bary(p.uv);
    auto flat = bary.x * flat_face[0] + bary.y * flat_face[1] +
                bary.z * flat_face[2];

    auto v      = flat - tr2d[(k + 2) % 3];
    auto w      = tr2d[k] - tr2d[(k + 2) % 3];
    auto u      = tr2d[(k + 1) % 3] - tr2d[(k + 2) % 3];
    auto theta  = angle(v, w);
    auto phi3D  = angle(w, u);
    auto vid0   = triangles[tid][k];
    auto s      = angles[vid].size();
    entry       = node_is_neighboor(solver, vid, vid0);
    auto theta0 = angles[vid][entry];
    auto phi2D  = (entry + 2 == s) ? 2 * pif - theta0
                                  : angles[vid][entry + 2] - theta0;
    auto scale_factor = phi3D / phi2D;
    ++entry;

    return theta0 + theta * scale_factor;
  }
}

inline float get_angle(const geodesic_solver& solver,
    const vector<vector<float>>& angles, int vid0, int vid1, int& entry) {
  entry = node_is_neighboor(solver, vid0, vid1);
  assert(entry >= 0);
  return angles[vid0][entry];
}

inline bool set_ord(int s, int prev_entry, int next_entry, bool nei_is_dual) {
  int CCW_count, CW_count;
  if (prev_entry < next_entry) {
    CCW_count = next_entry - prev_entry;
    CW_count  = prev_entry + s - next_entry;
  } else {
    CCW_count = s - prev_entry + next_entry;
    CW_count  = prev_entry - next_entry;
  }
  if (!nei_is_dual) --CCW_count;
  if (CCW_count < CW_count)
    return true;
  else
    return false;
}
inline bool set_ord(float teta_next, float teta_prev) {
  if (teta_next > teta_prev) {
    if (teta_next - teta_prev < pif)
      return true;
    else
      return false;
  } else {
    if (teta_prev - teta_next > pif)
      return true;
    else
      return false;
  }
}
void fill_strip(const vector<vector<int>>& v2t, int vid, int first, int last,
    bool nei_is_dual, vector<int>& strip, bool CCW) {
  auto start = first;
  auto end   = last;

  auto star = v2t[vid];
  auto s    = star.size();

  if (CCW && !nei_is_dual)
    end = (s - 1 + end) % s;  // I can stop one face earlier;

  if (start == end) {
    if (strip.back() != star[start]) strip.push_back(star[start]);

  } else if (CCW) {
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
static float get_angle(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, int parent, const mesh_point& p,
    vector<int>& strip, int& entry) {
  auto k = -1;
  strip  = {p.face};
  if (is_vert(p, k)) {
    auto vid   = triangles[p.face][k];
    auto entry = node_is_neighboor(solver, vid, parent);
    assert(entry >= 0);
    auto star        = v2t[vid];
    auto s           = star.size();
    auto it          = find(star.begin(), star.end(), p.face);
    auto first       = distance(star.begin(), it);
    auto last        = (entry % 2) ? (entry - 1) / 2 : (entry / 2) % s;
    bool CCW         = set_ord(s, first, last, entry % 2);
    bool nei_is_dual = entry % 2;
    fill_strip(v2t, vid, first, last, nei_is_dual, strip, CCW);
    if (entry % 2) {
      first    = (entry - 1) / 2;
      auto tid = opposite_face(triangles, adjacencies, v2t[vid][first], vid);
      strip.push_back(tid);
    }
    entry = node_is_neighboor(solver, parent, vid);
    return angles[parent][entry];

  } else {
    auto h = find_in_vec(triangles[p.face], parent);
    if (h == -1) {
      for (auto i = 0; i < 3; ++i) {
        auto adj = adjacencies[p.face][i];
        h        = find_in_vec(triangles[adj], parent);
        if (h != -1) {
          strip.push_back(adj);
          h = find_in_vec(adjacencies[adj], p.face);
          return get_rescaled_angle(solver, triangles, positions, adjacencies,
              angles, adj, parent, p, h, entry);
        }
      }
    } else {
      return get_rescaled_angle(solver, triangles, positions, adjacencies,
          angles, p.face, parent, p, -1, entry);
    }
  }
  return {};
}

void close_strip(const vector<vector<int>>& v2t, int vid, int prev_tri,
    int last_tri, vector<int>& strip) {
  auto star    = v2t[vid];
  auto s       = star.size();
  auto prev_it = find(star.begin(), star.end(), prev_tri);
  assert(prev_it != star.end());
  auto first   = distance(star.begin(), prev_it);
  auto next_it = find(star.begin(), star.end(), last_tri);
  assert(next_it != star.end());
  auto last = distance(star.begin(), next_it);
  auto CCW  = set_ord(s, first, last, true);
  fill_strip(v2t, vid, first, last, true, strip, CCW);
}

// particular case of "get strip" when one of the two point is the parent of
// the other so the size of the strip is one or two
static vector<int> short_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  auto offset = 0, entry = 0;
  auto strip = vector<int>{};

  if (is_vert(target, offset)) {
    auto vid = triangles[target.face][offset];
    get_angle(solver, triangles, positions, adjacencies, v2t, angles, vid,
        source, strip, entry);
    if (strip.back() != target.face) {
      close_strip(v2t, vid, strip.back(), target.face, strip);
    }
    reverse(strip.begin(), strip.end());
  } else if (is_vert(source, offset)) {
    auto vid = triangles[source.face][offset];
    get_angle(solver, triangles, positions, adjacencies, v2t, angles, vid,
        target, strip, entry);
    if (strip.back() != source.face) {
      close_strip(v2t, vid, strip.back(), source.face, strip);
    }
  } else {
    assert(false);
  }

  return strip;
}

// returns a strip of triangles such target belongs to the first one and
// source to the last one
//(TO DO:may be the names could change in order to get the call
// more consistent with the output)
vector<int> get_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target) {
  // TODO: ensure initialization of all variables
  // profile_function();
  if (target.face == source.face) return {target.face};
  auto parents = point_to_point_geodesic_path(
      solver, triangles, positions, adjacencies, source, target);
  auto N     = parents.size();
  auto first = 0, last = 0, prev_entry = 0, next_entry = 0;
  auto strip_to_point = vector<int>{}, strip = vector<int>{};
  auto teta_prev = 0.0f, teta_next = 0.0f;
  auto CCW = false, nei_is_dual = false;
  if (N == 0) {
    return short_strip(
        solver, triangles, positions, adjacencies, v2t, angles, source, target);
  } else if (N == 1) {
    auto v      = parents[0];
    teta_prev   = get_angle(solver, triangles, positions, adjacencies, v2t,
        angles, v, target, strip, prev_entry);
    teta_next   = get_angle(solver, triangles, positions, adjacencies, v2t,
        angles, v, source, strip_to_point, next_entry);
    first       = find_in_vec(v2t[v], strip.back());
    nei_is_dual = next_entry % 2;
    last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
    CCW         = set_ord(v2t[v].size(), first, last, nei_is_dual);
    fill_strip(v2t, v, first, last, nei_is_dual, strip, CCW);
    close_strip(v2t, v, strip.back(), strip_to_point.back(), strip);
    if (strip.back() == strip_to_point.back()) strip_to_point.pop_back();
    reverse(strip_to_point.begin(), strip_to_point.end());
    strip.insert(strip.end(), strip_to_point.begin(), strip_to_point.end());
    return strip;
  } else {
    teta_prev = get_angle(solver, triangles, positions, adjacencies, v2t,
        angles, parents[0], target, strip, prev_entry);
  }

  for (auto i = 0; i < N; ++i) {
    auto v = parents[i];

    if (i == N - 1) {
      first     = find_in_vec(v2t[v], strip.back());
      teta_next = get_angle(solver, triangles, positions, adjacencies, v2t,
          angles, v, source, strip_to_point, next_entry);
      last      = find_in_vec(v2t[v], strip_to_point.back());
      CCW       = set_ord(v2t[v].size(), first, last, next_entry % 2);
    } else {
      first = find_in_vec(v2t[v], strip.back());
      assert(first != -1);
      next_entry  = node_is_neighboor(solver, v, parents[i + 1]);
      nei_is_dual = next_entry % 2;
      last        = (nei_is_dual) ? (next_entry - 1) / 2 : next_entry / 2;
      CCW         = set_ord(v2t[v].size(), first, last, nei_is_dual);
    }

    fill_strip(v2t, v, first, last, nei_is_dual, strip, CCW);

    if (nei_is_dual && i != N - 1) {
      auto tid = opposite_face(triangles, adjacencies, strip.back(), v);
      strip.push_back(tid);
    }
  }

  close_strip(v2t, parents.back(), strip.back(), strip_to_point.back(), strip);
  if (strip.back() == strip_to_point.back()) strip_to_point.pop_back();
  reverse(strip_to_point.begin(), strip_to_point.end());
  strip.insert(strip.end(), strip_to_point.begin(), strip_to_point.end());
  return strip;
}

// compute the distance between a point p and some vertices around him
// handling concave path
static vector<pair<int, float>> nodes_around_point(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& p) {
  auto offset = 0;
  auto nodes  = vector<pair<int, float>>{};
  if (is_vert(p, offset)) {
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

      auto CW_pid = adjacencies[pid][i];
      auto opp    = opposite_vertex(triangles, adjacencies, pid, i);
      auto strip  = vector<int>{CW_pid, pid};
      auto pos2d  = zero2f;
      auto l      = length_by_flattening(
          triangles, positions, adjacencies, p, strip, pos2d);

      nodes.push_back({opp, l});
    }
  }

  return nodes;
}

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
  auto vid    = -1;
  auto lambda = flt_max;
  if (nbr.size() == 1)
    vid = parents[nbr[0].first];
  else {
    for (auto i = 0; i < nbr.size(); ++i) {
      auto val = f[nbr[i].first] + nbr[i].second;
      if (val < lambda) {
        lambda = val;
        vid    = nbr[i].first;
      }
    }
  }
  return vid;
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

vector<float> compute_geodesic_distances(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const geodesic_solver& solver, const vector<mesh_point>& sources) {
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
// IMPLEMENTATION OF MESH IO
// -----------------------------------------------------------------------------
namespace yocto {

// Get extension (not including '.').
static string get_extension(const string& filename) {
  auto pos = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Load ply mesh
[[nodiscard]] bool load_mesh(const string& filename, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, string& error, bool flip_texcoord) {
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

  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // open ply
    auto ply_guard = std::make_unique<ply_model>();
    auto ply       = ply_guard.get();
    if (!load_ply(filename, ply, error)) return false;

    // gets vertex
    get_positions(ply, positions);
    get_normals(ply, normals);
    get_texcoords(ply, texcoords, flip_texcoord);
    get_colors(ply, colors);

    // get faces
    get_triangles(ply, triangles);

    if (positions.empty()) return shape_error();
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    // load obj
    auto obj_guard = std::make_unique<obj_model>();
    auto obj       = obj_guard.get();
    if (!load_obj(filename, obj, error, true)) return false;

    // get shape
    if (obj->shapes.empty()) return shape_error();
    if (obj->shapes.size() > 1) return shape_error();
    auto shape = obj->shapes.front();
    if (shape->points.empty() && shape->lines.empty() && shape->faces.empty())
      return shape_error();

    // decide what to do and get properties
    auto materials  = vector<obj_material*>{};
    auto ematerials = vector<int>{};
    if (!shape->faces.empty()) {
      get_triangles(shape, triangles, positions, normals, texcoords, materials,
          ematerials, flip_texcoord);
    } else {
      return shape_error();
    }

    if (positions.empty()) return shape_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
[[nodiscard]] bool save_mesh(const string& filename,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, string& error, bool ascii,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // create ply
    auto ply_guard = std::make_unique<ply_model>();
    auto ply       = ply_guard.get();
    add_positions(ply, positions);
    add_normals(ply, normals);
    add_texcoords(ply, texcoords, flip_texcoord);
    add_colors(ply, colors);
    add_triangles(ply, triangles);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj_guard = std::make_unique<obj_model>();
    auto obj       = obj_guard.get();
    auto oshape    = add_shape(obj);
    if (!triangles.empty()) {
      set_triangles(
          oshape, triangles, positions, normals, texcoords, {}, flip_texcoord);
    } else {
      return shape_error();
    }
    auto err = ""s;
    if (!save_obj(filename, obj, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

// Load ply mesh
[[nodiscard]] bool load_lines(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, string& error, bool flip_texcoord) {
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

  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // open ply
    auto ply_guard = std::make_unique<ply_model>();
    auto ply       = ply_guard.get();
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
    auto obj_guard = std::make_unique<obj_model>();
    auto obj       = obj_guard.get();
    if (!load_obj(filename, obj, error, true)) return false;

    // get shape
    if (obj->shapes.empty()) return shape_error();
    if (obj->shapes.size() > 1) return shape_error();
    auto shape = obj->shapes.front();
    if (shape->points.empty() && shape->lines.empty() && shape->faces.empty())
      return shape_error();

    // decide what to do and get properties
    auto materials  = vector<obj_material*>{};
    auto ematerials = vector<int>{};
    if (!shape->faces.empty()) {
      get_lines(shape, lines, positions, normals, texcoords, materials,
          ematerials, flip_texcoord);
    } else {
      return shape_error();
    }

    if (positions.empty()) return shape_error();
    return true;
  } else {
    return format_error();
  }
}

// Save ply mesh
[[nodiscard]] bool save_lines(const string& filename,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, string& error, bool ascii,
    bool flip_texcoord) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    // create ply
    auto ply_guard = std::make_unique<ply_model>();
    auto ply       = ply_guard.get();
    add_positions(ply, positions);
    add_normals(ply, normals);
    add_texcoords(ply, texcoords, flip_texcoord);
    add_colors(ply, colors);
    add_lines(ply, lines);
    if (!save_ply(filename, ply, error)) return false;
    return true;
  } else if (ext == ".obj" || ext == ".OBJ") {
    auto obj_guard = std::make_unique<obj_model>();
    auto obj       = obj_guard.get();
    auto oshape    = add_shape(obj);
    if (!lines.empty()) {
      set_lines(
          oshape, lines, positions, normals, texcoords, {}, flip_texcoord);
    } else {
      return shape_error();
    }
    auto err = ""s;
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
    const vector<vec2f>& texcoords, const vector<vec3f>& colors, bool verbose) {
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
