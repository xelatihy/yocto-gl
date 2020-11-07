#pragma once
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_shape.h>

#include <deque>
#include <queue>

using namespace yocto;
using namespace std;

// Taking closest vertex of an intersection
int closest_vertex(const vector<vec3i> &triangles, vec2f uv, int element) {
  auto tr = triangles[element];
  if (uv.x < 0.5f && uv.y < 0.5f) return tr.x;
  if (uv.x > uv.y) return tr.y;
  return tr.z;
}

vec3f eval_position(const generic_shape *shape, int element, const vec2f &uv) {
  return eval_position(shape->triangles, shape->positions, {element, uv});
}
vec3f eval_normal(const generic_shape *shape, int element, const vec2f &uv) {
  return eval_normal(shape->triangles, shape->normals, {element, uv});
}
enum struct axes { x, y, z };
auto const symmetric_axes = vector<std::string>{"asse x", "asse y", "asse z"};

enum struct brush_type { gaussian, texture, smooth };
auto const brushes_names = vector<std::string>{
    "gaussian brush", "texture brush", "smooth brush"};

// To obtain symmetric from stroke result
inline vector<pair<vec3f, vec3f>> symmetric_stroke(
    vector<pair<vec3f, vec3f>> &pairs, generic_shape *shape, shape_bvh &tree,
    vector<int> &symmetric_stroke_sampling, axes &axis) {
  vector<pair<vec3f, vec3f>> symmetric_pairs;
  if (pairs.empty()) return symmetric_pairs;
  for (int i = 0; i < pairs.size(); i++) {
    auto ray = ray3f{};
    ray.d    = pairs[i].first;

    if (axis == axes::x) {
      ray.o   = vec3f{0.0f, ray.d.y, 0.0f};
      ray.d   = normalize(ray.d - ray.o);
      ray.d.x = -ray.d.x;
    }

    if (axis == axes::y) {
      ray.o   = vec3f{ray.d.x, 0.0f, 0.0f};
      ray.d   = normalize(ray.d - ray.o);
      ray.d.y = -ray.d.y;
    }

    if (axis == axes::z) {
      ray.o   = vec3f{0.0f, ray.d.y, 0.0f};
      ray.d   = normalize(ray.d - ray.o);
      ray.d.z = -ray.d.z;
    }

    auto inter = intersect_triangles_bvh(
        tree, shape->triangles, shape->positions, ray);
    if (!inter.hit) continue;
    auto pos  = eval_position(shape, inter.element, inter.uv);
    auto nor  = eval_normal(shape, inter.element, inter.uv);
    auto pair = std::pair<vec3f, vec3f>{pos, nor};
    symmetric_stroke_sampling.push_back(
        closest_vertex(shape->triangles, inter.uv, inter.element));
    symmetric_pairs.push_back(pair);
  }
  return symmetric_pairs;
}

// Project a vector on a plane, maintaining vector lenght
inline vec2f project_onto_plane(const mat3f &basis, const vec3f &p) {
  auto v  = p - dot(p, basis.z) * basis.z;
  auto v1 = vec2f{dot(v, basis.x), dot(v, basis.y)};
  return v1 * (length(p) / length(v1));
}

// Planar coordinates by local computation between neighbors
inline void compute_coordinates(vector<vec2f> &coords, vector<mat3f> &frames,
    vector<vec3f> &positions, int node, int neighbor, float weight) {
  auto current_coord = coords[node];
  auto edge          = positions[node] - positions[neighbor];
  auto projection    = project_onto_plane(frames[neighbor], edge);
  auto new_coord     = coords[neighbor] + projection;
  auto avg_lenght    = (length(current_coord) + length(new_coord)) / 2;
  auto new_dir       = normalize(current_coord + new_coord);
  coords[node] = current_coord == zero2f ? new_coord : new_dir * avg_lenght;

  // following doesn't work
  // coords[node] = current_coord + (weight * (coords[neighbor] + projection));
}

// Frame by local computation between neighbors
inline void compute_frame(vector<mat3f> &frames, vector<vec3f> &normals,
    int node, int neighbor, float weight) {
  auto current_dir = frames[node].x;
  auto rotation    = basis_fromz(normals[neighbor]) *
                  transpose(basis_fromz(normals[node]));
  auto neighbor_dir = frames[neighbor].x;
  current_dir       = current_dir + (rotation * weight * neighbor_dir);
  frames[node].z    = normals[node];
  frames[node].y    = cross(frames[node].z, normalize(current_dir));
  frames[node].x    = cross(frames[node].y, frames[node].z);
}

// Classic Dijkstra
template <typename Update>
void dijkstra(const geodesic_solver &solver, const vector<int> &sources,
    vector<float> &distances, float max_distance, Update &&update) {
  auto compare = [&](int i, int j) { return distances[i] > distances[j]; };
  std::priority_queue<int, vector<int>, decltype(compare)> queue(compare);

  // setup queue
  for (auto source : sources) queue.push(source);

  while (!queue.empty()) {
    int node = queue.top();
    queue.pop();

    auto distance = distances[node];
    if (distance > max_distance) continue;  // early exit

    for (auto arc : solver.graph[node]) {
      auto new_distance = distance + arc.length;

      update(node, arc.node, new_distance);

      if (new_distance < distances[arc.node]) {
        distances[arc.node] = new_distance;
        queue.push(arc.node);
      }
    }
  }
}

// Compute initial frames of stroke sampling vertices
inline void compute_stroke_frames(vector<mat3f> &frames,
    vector<vec3f> &positions, vector<vec3f> &normals,
    vector<int> &stroke_sampling) {
  // frames follow stroke direction
  for (int i = 0; i < stroke_sampling.size() - 1; i++) {
    int  curr    = stroke_sampling[i];
    int  next    = stroke_sampling[i + 1];
    auto dir     = positions[next] - positions[curr];
    auto z       = normals[curr];
    auto y       = cross(z, normalize(dir));
    auto x       = cross(y, z);
    frames[curr] = {x, y, z};
  }

  if (stroke_sampling.size() == 1) {
    frames[stroke_sampling[0]] = basis_fromz(normals[stroke_sampling[0]]);
  } else {
    int  final    = stroke_sampling[stroke_sampling.size() - 1];
    int  prev     = stroke_sampling[stroke_sampling.size() - 2];
    auto dir      = positions[final] - positions[prev];
    auto z        = normals[final];
    auto y        = cross(z, normalize(dir));
    auto x        = cross(y, z);
    frames[final] = {x, y, z};
  }

  // average frames direction of middle stroke vertices
  for (int i = 1; i < stroke_sampling.size() - 1; i++) {
    int  curr    = stroke_sampling[i];
    int  next    = stroke_sampling[i + 1];
    int  prev    = stroke_sampling[i - 1];
    auto dir     = frames[prev].x + frames[next].x;
    auto z       = normals[curr];
    auto y       = cross(z, normalize(dir));
    auto x       = cross(y, z);
    frames[curr] = {x, y, z};
  }
}

// To take shape positions indices associate with planar coordinates
inline vector<int> stroke_parameterization(geodesic_solver &solver,
    vector<vec2f> &coords, vector<int> &stroke_sampling,
    vector<vec3f> &positions, vector<vec3f> &normals, float radius) {
  if (stroke_sampling.empty()) return vector<int>{};

  // init params
  std::set<int> vertices;  // to avoid duplicates

  auto visited = vector<bool>(positions.size(), false);
  for (auto sample : stroke_sampling) visited[sample] = true;

  coords                     = vector<vec2f>(solver.graph.size(), zero2f);
  coords[stroke_sampling[0]] = {radius, radius};
  vertices.insert(stroke_sampling[0]);
  for (int i = 1; i < stroke_sampling.size(); i++) {
    auto edge = positions[stroke_sampling[i]] -
                positions[stroke_sampling[i - 1]];
    coords[stroke_sampling[i]] = {
        coords[stroke_sampling[i - 1]].x + length(edge), radius};
    vertices.insert(stroke_sampling[i]);
  }

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : stroke_sampling) distances[sample] = 0.0f;

  auto frames = vector<mat3f>(positions.size(), identity3x3f);
  compute_stroke_frames(frames, positions, normals, stroke_sampling);

  auto update = [&](int node, int neighbor, float new_distance) {
    vertices.insert(node);
    if (!visited[neighbor]) return;
    float weight = length(positions[neighbor] - positions[node]) + flt_eps;
    weight       = 1.0f / weight;
    compute_coordinates(coords, frames, positions, node, neighbor, weight);
    compute_frame(frames, normals, node, neighbor, weight);
    visited[node] = true;
  };
  dijkstra(solver, stroke_sampling, distances, radius, update);

  auto vec_vertices = vector<int>(vertices.begin(), vertices.end());

  // conversion in [0, 1]
  for (int i = 0; i < vertices.size(); i++)
    coords[vec_vertices[i]] /= radius * 2.0f;

  return vec_vertices;
}
