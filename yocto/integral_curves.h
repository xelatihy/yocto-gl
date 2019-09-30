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
