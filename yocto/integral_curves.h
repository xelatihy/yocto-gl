#include "yocto_shape.h"

using std::vector;
using namespace yocto;

// @Duplicate: function already defined in yocto_shape.cpp
inline int opposite_vertex(const vec3i& tr, const vec2i& edge) {
  for (int i = 0; i < 3; ++i) {
    if (tr[i] != edge.x && tr[i] != edge.y) return tr[i];
  }
  return -1;
}

inline vec2i opposite_edge(int vertex, const vec3i& tr) {
  for (int i = 0; i < 3; ++i) {
    if (tr[i] == vertex) return {tr[(i + 1) % 3], tr[(i + 2) % 3]};
  }
  return {-1, -1};
}

inline vec2i get_edge(const vec3i& triangle, int i) {
  auto x = triangle[i], y = triangle[i < 2 ? i + 1 : 0];
  return vec2i{x, y};
};

inline vec2i make_edge(int x, int y) {
  return x < y ? vec2i{x, y} : vec2i{y, x};
}

inline vec2i make_edge(const vec2i& e) { return make_edge(e.x, e.y); }

inline bool same_edge(const vec2i& a, const vec2i& b) {
  return make_edge(a) == make_edge(b);
}

inline int adjacent_face(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int face, const vec2i& edge) {
  // Given a face and an edge, return the adjacent face
  for (int k = 0; k < 3; ++k) {
    auto e = get_edge(triangles[face], k);
    if (e == edge) return adjacency[face][k];
    if (vec2i{e.y, e.x} == edge) return adjacency[face][k];
  }
  return -1;
}

inline int find(const vec3i& v, int x) {
  if (v.x == x) return 0;
  if (v.y == x) return 1;
  if (v.z == x) return 2;
  return -1;
}

struct lerp_point {
  vec2i edge;
  int   face;
  float alpha;

  lerp_point() {}
};

// Description of a discrete path along the surface of the mesh.
struct Path {
  int                start, end;
  vector<lerp_point> lerps;
};

template <typename Type>
inline bool contains(const vector<Type>& v, const Type& value) {
  return find(v.begin(), v.end(), value) != v.end();
}

inline bool contains(const vec3i& v, int a) {
  return v.x == a or v.y == a or v.z == a;
}

inline vector<int> get_face_ring(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int face, int vertex) {
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

  while (not queue.empty()) {
    int f = queue.back();
    queue.pop_back();
    result.push_back(f);

    // for (auto[edge, k] : edges(triangles[f])) {
    for (int k = 0; k < 3; ++k) {
      auto edge = vec2i{triangles[f][k], triangles[f][(k + 1) % 3]};
      if (edge.x == vertex or edge.y == vertex) {
        int neighbor_face = adjacency[f][k];
        if (neighbor_face == -1) continue;
        if (not contains(result, neighbor_face)) {
          queue.push_back(neighbor_face);
          result.push_back(neighbor_face);
        }
      }
    }
  }

  return result;
}

inline vec3f gradient_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& field, int face) {
  auto  result = zero3f;
  auto& tr     = triangles[face];
  auto& a      = positions[tr.x];
  auto& b      = positions[tr.y];
  auto& c      = positions[tr.z];
  auto  ab     = b - a;
  auto  bc     = c - b;
  auto  ca     = a - c;

  auto normal = normalize(cross(ca, ab));
  // auto normal = cross(ca, ab);
  result += field[tr.x] * cross(normal, bc);
  result += field[tr.y] * cross(normal, ca);
  result += field[tr.z] * cross(normal, ab);
  return normalize(result);
}

const float eps = 0.0001f;

float step_from_point_to_edge(
    const vec3f& right, const vec3f& left, const vec3f& direction) {
  vec3f normal = cross(right, left);
  mat3f inverse_transform;
  inverse_transform.x = right - left;
  inverse_transform.y = left;
  inverse_transform.z = normal;
  auto transform      = inverse(inverse_transform);
  auto dir            = transform * direction;
  return clamp(1 - dir.x / dir.y, 0.0f + eps, 1.0f - eps);
}

inline bool is_direction_inbetween(const vec3f& right, const vec3f& left,
    const vec3f& direction, float threshold = 0) {
  vec3f normal      = cross(right, left);
  auto  cross_right = cross(right, direction);
  auto  cross_left  = cross(direction, left);
  return dot(cross_right, normal) > threshold and
         dot(cross_left, normal) > threshold;
}

inline pair<float, bool> step_from_edge_to_edge(const vec3f& point,
    const vec3f& a, const vec3f& b, const vec3f& c, const vec3f& direction) {
  //      b
  //     /\
    //    /  \
    //   /    \
    //  /______\
    // c        a

  vec3f right  = a - point;
  vec3f left   = c - point;
  vec3f front  = b - point;
  vec3f normal = triangle_normal(a, b, c);

  bool right_side = dot(cross(direction, front), normal) > 0;
  if (right_side) {
    bool below_edge = dot(cross(direction, right), normal) > 0;
    if (below_edge) return {eps, true};

    auto x = step_from_point_to_edge(right, front, direction);
    return {x, true};
  } else {
    bool below_edge = dot(cross(left, direction), normal) > 0;
    if (below_edge) return {1.0f - eps, false};

    auto x = step_from_point_to_edge(front, left, direction);
    return {x, false};
  }
}

inline lerp_point make_lerp(vec2i e, int f, float a) {
  lerp_point lerp;
  lerp.edge  = e;
  lerp.face  = f;
  lerp.alpha = a;
  //    lerp.loc   = l;

  return lerp;
}

inline lerp_point step_from_point(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, const vector<float>& field, int vertex,
    int start_face, int tag = -1) {
  auto triangle_fan = get_face_ring(triangles, adjacency, start_face, vertex);

  auto best_alignment = 0.0;
  auto fallback_lerp  = make_lerp((vec2i{-1, -1}), -1, 0);

  for (int i = 0; i < triangle_fan.size(); ++i) {
    int face = triangle_fan[i];
    if (tag != -1 and tags[face] != tag) continue;

    auto edge = opposite_edge(vertex, triangles[face]);
    //        assert(edge != (vec2i{-1, -1}));
    if (edge == vec2i{-1, -1}) {
      printf("edge is {-1, -1}\n");
      continue;
    }

    auto a         = positions[vertex];
    auto b         = positions[edge.x];
    auto c         = positions[edge.y];
    auto left      = c - a;
    auto right     = b - a;
    auto direction = gradient_face(triangles, positions, field, face);

    auto right_dot = dot(normalize(right), direction);
    auto left_dot  = dot(normalize(left), direction);

    // Look for the most aligned edge with the gradient
    // direction. If no other face is suitable for selection, we return the
    // most aligned edge as result (fallback_lerp)
    if (max(right_dot, left_dot) > best_alignment) {
      best_alignment = max(right_dot, left_dot);
      float alpha    = right_dot > left_dot ? 0.0f + eps : 1.0f - eps;
      fallback_lerp  = make_lerp(edge, face, alpha);
    }

    // Check if gradient direction is in between ac and ab.
    if (is_direction_inbetween(right, left, direction)) {
      float alpha = step_from_point_to_edge(right, left, direction);
      return make_lerp(edge, face, alpha);
    }
  }

  // assert(fallback_lerp.face != -1);
  return fallback_lerp;
}

inline Path follow_gradient_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from) {
  // TRACE_FUNCTION
  auto lerps = vector<lerp_point>();
  {
    auto lerp = step_from_point(
        triangles, positions, adjacency, tags, field, from, -1, tag);
    if (lerp.face == -1) return {};
    lerps.push_back(lerp);
  }

  const int num_steps = 10000;

  for (int i = 0; i < num_steps; i++) {
    auto[old_edge, old_face, old_alpha] = lerps.back();
    assert(old_face != -1);
    vec3f point = (1.0f - old_alpha) * positions[old_edge.x] +
                  old_alpha * positions[old_edge.y];

    int face = adjacent_face(triangles, adjacency, old_face, old_edge);
    if (face == -1) {
      lerps.push_back(make_lerp(vec2i{-1, -1}, face, 0));
      return Path{from, -1, lerps};
    }

    if (tags[face] != tag) {
      int k  = find(triangles[face], old_edge.x);
      int to = triangles[face][(k + 1) % 3];
      // int   to   = opposite_vertex(old_edge, triangles[face]);

      // @Hack!: We store the tag of the reached region in edge.x
      vec2i edge = {to, tags[face]};
      lerps.push_back(make_lerp(edge, face, 0));
      return Path{from, to, lerps};
    }

    vec3f direction = gradient_face(triangles, positions, field, face);

    int front_idx = opposite_vertex(triangles[face], old_edge);
    if (front_idx == -1) {
      printf("front_idx is -1!\n");
      return {};
    }

    assert(old_alpha > 0);
    assert(old_alpha < 1);

    auto& a = positions[old_edge.x];
    auto& b = positions[front_idx];
    auto& c = positions[old_edge.y];

    auto[x, step_right] = step_from_edge_to_edge(point, a, b, c, direction);

    vec2i edge = old_edge;
    if (step_right) {
      point  = (1 - x) * a + x * b;
      edge.y = front_idx;
    } else {
      point  = (1 - x) * b + x * c;
      edge.x = front_idx;
    }
    lerps.push_back(make_lerp(edge, face, x));
  }

  assert(0 && "path ended in nowhere");
  return Path{from, 0, lerps};
}

inline Path follow_gradient_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from,
    int to) {
  // TRACE_FUNCTION
  auto lerps = vector<lerp_point>();

  lerps.push_back(
      step_from_point(triangles, positions, adjacency, tags, field, from, -1));
  //  assert(
  //      opposite_vertex(lerps.back().edge, triangles[lerps.back().face]) ==
  //      from);

  const int num_steps = 10000;

  for (int i = 0; i < num_steps; i++) {
    auto[old_edge, old_face, old_alpha] = lerps.back();
    assert(old_face != -1);
    vec3f point = (1.0f - old_alpha) * positions[old_edge.x] +
                  old_alpha * positions[old_edge.y];

    int face = adjacent_face(triangles, adjacency, old_face, old_edge);
    if (face == -1) {
      //          assert(face != -1);
      break;
    }

    if (contains(triangles[face], to)) {
      // for (auto[edge, k] : edges(triangles[face])) {
      for (int k = 0; k < 3; ++k) {
        auto edge = vec2i{triangles[face][k], triangles[face][(k + 1) % 3]};
        if (edge.x == to) {
          lerps.push_back(make_lerp(edge, face, 0));
          return {from, to, lerps};
        }
      }
    }
    vec3f direction = gradient_face(triangles, positions, field, face);

    int front_idx = opposite_vertex(triangles[face], old_edge);
    if (front_idx == -1) {
      printf("front_idx is -1!\n");
      break;
    }

    assert(old_alpha >= 0);
    assert(old_alpha <= 1);
    if (old_alpha == 0 or old_alpha == 1) {
      int  vertex = old_alpha == 0 ? old_edge.x : old_edge.y;
      auto lerp   = step_from_point(
          triangles, positions, adjacency, tags, field, vertex, old_face);
      lerps.push_back(lerp);
      if (lerp.alpha == 0 and lerp.edge.x == to) break;
      if (lerp.alpha == 1 and lerp.edge.y == to) break;
      continue;
    }

    auto& a = positions[old_edge.x];
    auto& b = positions[front_idx];
    auto& c = positions[old_edge.y];

    auto[x, step_right] = step_from_edge_to_edge(point, a, b, c, direction);

    vec2i edge = old_edge;
    if (step_right) {
      point  = (1 - x) * a + x * b;
      edge.y = front_idx;
    } else {
      point  = (1 - x) * b + x * c;
      edge.x = front_idx;
    }
    if (opposite_vertex(triangles[face], edge) == -1) {
      assert(0 && "opposite vertex == -1");
      break;
    }

    lerps.push_back(make_lerp(edge, face, x));
    if (x == 0 and edge.x == to) break;
    if (x == 1 and edge.y == to) break;
  }

  return {from, to, lerps};
}

struct State {
  // solver used for scalar field computation
  geodesic_solver solver;

  // mesh data
  vector<vec3f> positions;
  vector<vec3f> normals;
  vector<vec3i> triangles;

  // triangle adjacency used for fast boundary computation
  vector<vec3i> triangle_graph;

  // editing data
  vector<int> tags;  // per-triangle region tag
};

inline int add_triangle(State& state, const vec3i& triangle = {}, int tag = -1,
    const vec3i& neighbors = {-1, -1, -1}) {
  assert(state.triangles.size() == state.tags.size());
  assert(state.triangles.size() == state.triangle_graph.size());
  int index = state.triangles.size();
  state.triangles.push_back(triangle);
  state.tags.push_back(tag);
  state.triangle_graph.push_back(neighbors);
  return index;
}

vector<vec3f> points_from_lerps(
    const vector<vec3f>& positions, const Path& path) {
  auto result = vector<vec3f>();
  result.reserve(path.lerps.size() + 1);
  result.push_back(positions[path.start]);
  // + 0.0001 * state.normals[path.start]);
  if (path.lerps.empty()) return {};
  for (int i = 0; i < path.lerps.size() - 1; ++i) {
    auto[edge, face, x] = path.lerps[i];
    auto a        = positions[edge.x];  // + 0.0001 * state.normals[edge.x];
    auto b        = positions[edge.y];  // + 0.0001 * state.normals[edge.y];
    auto position = (1 - x) * a + x * b;
    result.push_back(position);
  }
  if (path.end != -1) {
    result.push_back(positions[path.end]);
    //+ 0.0001 * state.normals[path.end]);
  }
  return result;
}

inline int add_node(
    geodesic_solver& solver, const vec3f& position, const vec3f& normal) {
  // assert(solver.positions.size() == solver.graph.size());
  // assert(solver.normals.size() == solver.graph.size());
  // solver.positions.push_back(position);
  // solver.normals.push_back(normal);
  solver.graph.emplace_back();
  solver.graph.back().reserve(solver.min_arcs);
  return (int)solver.graph.size() - 1;
}

inline void delete_arc(geodesic_solver& solver, int a, int b) {
  // search and delete b in a's adjacency
  for (int k = 0; k < solver.graph[a].size(); ++k) {
    if (solver.graph[a][k].node == b) {
      solver.graph[a].erase(solver.graph[a].begin() + k);
      break;
    }
    // assert(k != solver.graph[a].size() - 1);
  }
  // search and delete a in b's adjacency
  for (int k = 0; k < solver.graph[b].size(); ++k) {
    if (solver.graph[b][k].node == a) {
      solver.graph[b].erase(solver.graph[b].begin() + k);
      break;
    }
    // assert(k != solver.graph[b].size() - 1);
  }
}

inline int add_vertex(
    State& state, const vec3f& position, const vec3f& normal) {
  assert(state.positions.size() == state.normals.size());
  int index = state.positions.size();
  state.positions.push_back(position);
  state.normals.push_back(normal);
  add_node(state.solver, position, normal);
  return index;
}

inline int add_vertex(State& state, int a, int b, float coeff) {
  auto position = coeff * state.positions[a] + (1 - coeff) * state.positions[b];
  auto normal   = normalize(
      coeff * state.normals[a] + (1 - coeff) * state.normals[b]);
  return add_vertex(state, position, normal);
}

inline void connect_opposite_nodes(
    State& state, int f0, int f1, const vec2i& edge) {
  if (f0 != -1 and f1 != -1)
    connect_opposite_nodes(state.solver, state.positions, state.triangles[f0],
        state.triangles[f1], edge);
}

inline void disconnect_opposite_nodes(geodesic_solver& solver, const vec3i& f0,
    const vec3i& f1, const vec2i& edge) {
  /*
       a
      /\
     /  \
  b /____\
    \    /
     \  /
      \/
      c
  */

  int a = opposite_vertex(f0, edge);
  int c = opposite_vertex(f1, edge);
  if (a == -1 or c == -1) return;
  delete_arc(solver, a, c);
}

inline void disconnect_triangle(State& state, int face) {
  for (int k = 0; k < 3; ++k) {
    auto edge = vec2i{
        state.triangles[face][k], state.triangles[face][(k + 1) % 3]};
    int neighbor = state.triangle_graph[face][k];
    disconnect_opposite_nodes(
        state.solver, state.triangles[face], state.triangles[neighbor], edge);
  }
}

inline int adjacent_face(const State& state, int face, const vec2i& edge) {
  adjacent_face(state.triangles, state.triangle_graph, face, edge);
  return -1;
}

inline void update_adjacency(State& state, const vector<int>& updated_faces) {
  for (auto face : updated_faces) {
    //        for (auto [edge, k] : edges(state.triangles[face])) {
    for (int k : {0, 1, 2}) {
      vec2i edge = {
          state.triangles[face][k], state.triangles[face][(k + 1) % 3]};

      int neighbor = state.triangle_graph[face][k];
      if (neighbor == -1) continue;

      for (int j = 0; j < 3; ++j) {
        auto neighbor_edge = vec2i{state.triangles[neighbor][j],
            state.triangles[neighbor][(j + 1) % 3]};
        //      for (auto[neighbor_edge, j] : edges(state.triangles[neighbor]))
        //      {
        if (same_edge(edge, neighbor_edge)) {
          state.triangle_graph[neighbor][j] = face;
        }
      }
    }
  }
}

inline void connect_adjacent_nodes(State& state, int v0, int v1) {
  auto len = length(state.positions[v0] - state.positions[v1]);
  connect_nodes(state.solver, v0, v1, len);
}

inline void disconnect_adjacent_nodes(State& state, int v0, int v1) {
  if (v0 != -1 and v1 != -1) delete_arc(state.solver, v0, v1);
}

// Refer to slice_triangle.png to understand the names of the following
// functions.
inline void split_triangle(State& state, const vec3i& triangle, int fu, int fl,
    int fr, int vr, int vl, int erd, int eru, int elu, int eld, int cd) {
  auto[x, y, z] = triangle;

  disconnect_adjacent_nodes(state, y, z);
  disconnect_adjacent_nodes(state, z, x);

  state.triangles[fu]      = {vl, vr, z};
  state.triangle_graph[fu] = {fl, eru, elu};

  state.triangles[fl]      = {x, vr, vl};
  state.triangle_graph[fl] = {fr, fu, eld};

  state.triangles[fr]      = {x, y, vr};
  state.triangle_graph[fr] = {cd, erd, fl};

  connect_adjacent_nodes(state, y, vr);
  connect_adjacent_nodes(state, vr, z);
  connect_adjacent_nodes(state, z, vl);
  connect_adjacent_nodes(state, vl, x);
  connect_adjacent_nodes(state, x, vr);
  connect_adjacent_nodes(state, vl, vr);

  connect_opposite_nodes(state, fu, fl, {vl, vr});
  connect_opposite_nodes(state, fu, eru, {z, vr});
  connect_opposite_nodes(state, fu, elu, {z, vl});

  connect_opposite_nodes(state, fl, eld, {vl, x});
  connect_opposite_nodes(state, fl, fr, {x, vr});

  connect_opposite_nodes(state, fr, cd, {x, y});
  connect_opposite_nodes(state, fr, erd, {y, vr});
}

inline void split_triangle_(State& state, const vec3i& triangle, int face,
    int fl, int fr, int vr, int vl, int erd, int eru, int elu, int eld,
    int cd) {
  auto[x, y, z] = triangle;

  disconnect_adjacent_nodes(state, y, z);
  disconnect_adjacent_nodes(state, z, x);

  state.triangles[face]      = {vl, vr, z};
  state.triangle_graph[face] = {fr, eru, elu};

  state.triangles[fl]      = {vl, x, y};
  state.triangle_graph[fl] = {eld, cd, fr};

  state.triangles[fr]      = {vl, y, vr};
  state.triangle_graph[fr] = {fl, erd, face};

  connect_adjacent_nodes(state, y, vr);
  connect_adjacent_nodes(state, vr, z);
  connect_adjacent_nodes(state, z, vl);
  connect_adjacent_nodes(state, vl, x);
  connect_adjacent_nodes(state, y, vl);
  connect_adjacent_nodes(state, vl, vr);

  connect_opposite_nodes(state, face, fr, {vl, vr});
  connect_opposite_nodes(state, face, eru, {z, vr});
  connect_opposite_nodes(state, face, elu, {z, vl});

  connect_opposite_nodes(state, fl, eld, {vl, x});
  connect_opposite_nodes(state, fl, fr, {y, vl});
  connect_opposite_nodes(state, fl, cd, {x, y});

  connect_opposite_nodes(state, fr, erd, {y, vr});
}

inline void split_triangle(State& state, const vec3i& triangle, int fu, int fd,
    int v, int efd, int erd, int eru, int efu) {
  auto[x, y, z] = triangle;

  disconnect_adjacent_nodes(state, y, z);

  state.triangles[fu]      = {x, v, z};
  state.triangle_graph[fu] = {fd, eru, efu};

  state.triangles[fd]      = {x, y, v};
  state.triangle_graph[fd] = {efd, erd, fu};

  connect_adjacent_nodes(state, v, z);
  connect_adjacent_nodes(state, z, x);
  connect_adjacent_nodes(state, x, v);
  connect_adjacent_nodes(state, x, y);
  connect_adjacent_nodes(state, y, v);

  connect_opposite_nodes(state, fu, fd, {x, v});
  connect_opposite_nodes(state, fu, eru, {v, z});
  connect_opposite_nodes(state, fu, efu, {z, x});

  connect_opposite_nodes(state, fd, efd, {x, y});
  connect_opposite_nodes(state, fd, erd, {y, v});
}

inline bool slice_path(State& state, int tag, const Path& path, int tag_left,
    int tag_right, vector<int>& left_faces, vector<int>& right_faces) {
  auto& lerps = path.lerps;
  auto  start = path.start;
  left_faces.reserve(lerps.size());
  right_faces.reserve(lerps.size());

  {
    auto[edge, face, alpha] = lerps[0];
    int x                   = start;
    int y                   = edge.x;
    int z                   = edge.y;
    //    assert(opposite_vertex(edge, state.triangles[face]) == start);

    disconnect_triangle(state, face);

    auto v   = add_vertex(state, edge.x, edge.y, 1 - alpha);
    int  efd = adjacent_face(state, face, {x, y});
    int  efu = adjacent_face(state, face, {x, z});
    int  eru = -1;
    int  erd = -1;

    int fd = add_triangle(state);
    split_triangle(state, {x, y, z}, face, fd, v, efd, erd, eru, efu);
    // state.tags[fd]   = tag_right;
    // state.tags[face] = tag_left;
    right_faces.push_back(fd);
    left_faces.push_back(face);
  }

  int last_face_left  = lerps[0].face;
  int last_face_right = state.triangles.size() - 1;

  for (int i = 1; i < lerps.size() - 1; ++i) {
    auto[edge, face, alpha] = lerps[i];
    bool step_right         = lerps[i - 1].edge.x == edge.x;

    // int old_face_left  = lerps[i - 1].face;
    // int old_face_right = state.triangles.size() - 1;
    //        if (step_right) std::swap(old_face_left, old_face_right);
    auto triangle = state.triangles[face];

    if (!step_right) {
      int x = opposite_vertex(triangle, edge);
      if (x == -1) return false;
      assert(x != -1);
      assert(x == lerps[i - 1].edge.x);
      int  y  = edge.x;
      int  z  = edge.y;
      int  vl = state.positions.size() - 1;
      auto vr = add_vertex(state, edge.x, edge.y, 1 - alpha);

      disconnect_triangle(state, face);
      int erd = -1;
      int eru = -1;
      int elu = last_face_left;
      int eld = last_face_right;
      int cd  = state.triangle_graph[face][find(triangle, x)];
      //            assert(cd != -1);

      state.triangles[face] = {x, y, z};
      int fl                = add_triangle(state);
      int fr                = add_triangle(state);

      split_triangle(
          state, {x, y, z}, face, fl, fr, vr, vl, erd, eru, elu, eld, cd);
      // state.tags[face] = tag_left;
      // state.tags[fl]   = tag_right;
      // state.tags[fr]   = tag_right;
      left_faces.push_back(face);
      right_faces.push_back(fl);
      right_faces.push_back(fr);

      last_face_left  = face;
      last_face_right = fr;
    } else {
      int x = edge.y;
      int y = opposite_vertex(triangle, edge);
      if (y == -1) return false;
      assert(y != -1);
      int z = edge.x;

      int  vr = state.positions.size() - 1;
      auto vl = add_vertex(state, edge.x, edge.y, 1 - alpha);

      disconnect_triangle(state, face);
      int eld = -1;
      int elu = -1;
      int erd = last_face_left;
      int eru = last_face_right;
      int cd  = state.triangle_graph[face][find(triangle, x)];
      //            assert(cd != -1);

      state.triangles[face] = {x, y, z};
      int fr                = add_triangle(state);
      int fl                = add_triangle(state);

      // split_triangle_(editing_state& state, const vec3i& triangle, int
      // fu,
      //                            int fl, int fr, int vr, int vl, int
      //                            erd, int eru, int elu, int eld, int
      //                            cd) {
      split_triangle_(
          state, {x, y, z}, face, fl, fr, vr, vl, erd, eru, elu, eld, cd);
      // state.tags[face] = tag_right;
      // state.tags[fl]   = tag_left;
      // state.tags[fr]   = tag_left;
      right_faces.push_back(face);
      left_faces.push_back(fl);
      left_faces.push_back(fr);

      last_face_left  = fl;
      last_face_right = face;
    }
  }

  int end = lerps.back().edge.x;
  assert(end == path.end);

  if (lerps.back().face != -1) {
    auto[edge, face, alpha] = lerps.back();
    //        assert((state.tags[face] != tag) == (arrival == -1));

    //        auto prev_lerp = lerps[lerps.size() - 2];
    // int  old_face_left  = prev_lerp.face;
    // int  old_face_right = state.triangles.size() - 1;

    // bool step_right = prev_lerp.edge.x == edge.x;
    // if (step_right) std::swap(old_face_left, old_face_right);
    auto tr = state.triangles[face];

    int x = end;
    int y = tr[(find(tr, x) + 1) % 3];  // prev_lerp.edge.y;
    int z = tr[(find(tr, x) + 2) % 3];  // prev_lerp.edge.x;
    // assert(opposite_vertex(edge, state.triangles[face]) == start);
    disconnect_triangle(state, face);

    {
      auto v   = state.positions.size() - 1;
      auto adj = state.triangle_graph[face];
      // int efd = adjacent_face(state, face, {x, y});
      // int efu = adjacent_face(state, face, {x, z});
      int efd               = adj[find(tr, x)];
      int efu               = adj[find(tr, z)];
      int eru               = last_face_right;
      int erd               = last_face_left;
      state.triangles[face] = {x, y, z};

      connect_adjacent_nodes(state, v, y);
      //          length(state.positions[v] - state.positions[y]));
      connect_adjacent_nodes(state, v, z);
      //          length(state.positions[v] - state.positions[z]));

      int fd = add_triangle(state);

      split_triangle(state, {x, y, z}, face, fd, v, efd, erd, eru, efu);
      if (state.tags[face] == tag) {
        // state.tags[fd]   = tag_left;
        // state.tags[face] = tag_right;
        left_faces.push_back(fd);
        right_faces.push_back(face);

      } else {
        state.tags[fd] = state.tags[face];
        update_adjacency(state, {face, fd});
      }
    }
  }

  // for (int i = 0; i < state.triangles.size(); ++i) {
  //     if (state.tags[i] == tag_right) right_faces.push_back(i);
  //     if (state.tags[i] == tag_left) left_faces.push_back(i);
  // }

  // for (auto& f : left_faces) state.tags[f] = tag_left;
  // for (auto& f : right_faces) state.tags[f] = tag_right;

  update_adjacency(state, right_faces);
  update_adjacency(state, left_faces);

  // assert(state_is_sane(state));
  return true;
}
