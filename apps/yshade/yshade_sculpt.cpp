#include "yshade_sculpt.h"

#include <yocto/yocto_cli.h>
#include <yocto/yocto_color.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>

#include <deque>
#include <queue>
#include <set>

using namespace yocto;
using namespace std;

#ifdef _WIN32
#undef near
#undef far
#undef max
#endif

// Application state
struct shade_sculpt_state {
  // loading parameters
  string shapename   = "";
  string filename    = "shape.ply";
  string imagename   = "";
  string outname     = "out.ply";
  string texturename = "";
  string name        = "";

  // scene
  sceneio_scene ioscene  = sceneio_scene{};
  camera_handle iocamera = invalid_handle;

  // options
  shade_params drawgl_prms = {};

  // rendering state
  shade_scene glscene = {};
};

enum struct brush_type { gaussian, texture, smooth };
auto const brushes_names = vector<std::string>{
    "gaussian brush", "texture brush", "smooth brush"};

struct sculpt_params {
  brush_type type     = brush_type::gaussian;
  float      radius   = 0.35f;
  float      strength = 1.0f;
  bool       negative = false;
};

struct sculpt_state {
  shape_bvh           bvh         = {};
  hash_grid           grid        = {};
  vector<vector<int>> adjacencies = {};
  geodesic_solver     solver      = {};
  image_data          tex_image   = {};
};

// sculpt stroke
struct sculpt_stroke {
  vector<pair<vec3f, vec3f>> pairs     = {};
  vector<int>                sampling  = {};
  vec2f                      locked_uv = {};
  bool                       lock      = false;
};

// sculpt buffers
struct sculpt_buffers {
  vector<vec2f> coords        = {};
  vector<vec3f> old_positions = {};
  vector<vec3f> old_normals   = {};
};

// Initialize all sculpting parameters.
sculpt_state make_sculpt_state(
    const shape_data &shape, const scene_texture &texture) {
  auto state = sculpt_state{};
  state.bvh  = make_triangles_bvh(
      shape.triangles, shape.positions, shape.radius);
  state.grid       = make_hash_grid(shape.positions, 0.05f);
  auto adjacencies = face_adjacencies(shape.triangles);
  state.solver     = make_geodesic_solver(
      shape.triangles, adjacencies, shape.positions);
  state.adjacencies = vertex_adjacencies(shape.triangles, adjacencies);
  state.tex_image   = texture;
  return state;
}

sculpt_buffers make_sculpt_buffers(const shape_data &shape) {
  auto buffers          = sculpt_buffers{};
  buffers.old_positions = shape.positions;
  buffers.old_normals   = shape.normals;
  buffers.coords        = vector<vec2f>(shape.positions.size(), vec2f{0, 0});
  return buffers;
}

shape_data make_circle(
    const vec3f &center, const mat3f &basis, float radius, int steps) {
  // 4 initial vertices
  auto  lines    = make_lines({1, 4});
  vec3f next_dir = basis.x;
  for (auto line : lines.lines) {
    auto v1             = line.x;
    auto v2             = line.y;
    lines.positions[v1] = center + next_dir * radius;
    lines.normals[v1]   = basis.z;
    next_dir            = cross(basis.z, next_dir);
    lines.positions[v2] = center + next_dir * radius;
    lines.normals[v2]   = basis.z;
  }

  // create polylines and fix lenght to radius
  auto positions = vector<vec3f>{};
  auto normals   = vector<vec3f>{};
  auto lins      = vector<vec2i>{};
  steps          = int(steps / 4);
  for (auto line : lines.lines) {
    auto v1  = line.x;
    auto v2  = line.y;
    auto pos = lines.positions[v1];
    positions.push_back(pos);
    normals.push_back(lines.normals[v1]);
    auto dir    = lines.positions[v2] - lines.positions[v1];
    auto lenght = length(dir) / steps;
    dir         = normalize(dir);
    for (int i = 0; i < steps; i++) {
      auto new_pos = pos + dir * lenght * (i + 1);
      auto new_dir = normalize(new_pos - center);
      new_pos      = center + new_dir * radius;
      positions.push_back(new_pos);
      normals.push_back(basis.z);
      lins.push_back({int(positions.size() - 2), int(positions.size() - 1)});
    }
  }

  // apply
  lines.positions = positions;
  lines.normals   = normals;
  lines.lines     = lins;
  return lines;
}

// To visualize mouse intersection on mesh
shape_data make_cursor(const vec3f &position, const vec3f &normal, float radius,
    float height = 0.05f) {
  auto basis  = basis_fromz(normal);
  auto cursor = make_circle(position, basis, radius, 32);
  cursor.normals.clear();
  cursor.texcoords.clear();
  cursor.positions.push_back(position);
  cursor.positions.push_back(position + normal * 0.05f);
  cursor.lines.push_back(
      {int(cursor.positions.size() - 2), int(cursor.positions.size() - 1)});
  return cursor;
}

// Taking closest vertex of an intersection
int closest_vertex(
    const vector<vec3i> &triangles, int element, const vec2f &uv) {
  auto &triangle = triangles[element];
  if (uv.x < 0.5f && uv.y < 0.5f) return triangle.x;
  if (uv.x > uv.y) return triangle.y;
  return triangle.z;
}

// Project a vector on a plane, maintaining vector length
inline vec2f project_onto_plane(const mat3f &basis, const vec3f &p) {
  auto v  = p - dot(p, basis.z) * basis.z;
  auto v1 = vec2f{dot(v, basis.x), dot(v, basis.y)};
  return v1 * (length(p) / length(v1));
}

// Planar coordinates by local computation between neighbors
inline void compute_coordinates(vector<vec2f> &coords,
    const vector<mat3f> &frames, const vector<vec3f> &positions, int node,
    int neighbor, float weight) {
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
inline void compute_frame(vector<mat3f> &frames, const vector<vec3f> &normals,
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
    const vector<vec3f> &positions, const vector<vec3f> &normals,
    const vector<int> &stroke_sampling) {
  // frames follow stroke direction
  for (int i = 0; i < (int)stroke_sampling.size() - 1; i++) {
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
vector<int> stroke_parameterization(vector<vec2f> &coords,
    const geodesic_solver &solver, const vector<int> &sampling,
    const vector<vec3f> &positions, const vector<vec3f> &normals,
    float radius) {
  if (normals.empty()) return vector<int>{};

  // init params
  std::set<int> vertices;  // to avoid duplicates

  auto visited = vector<bool>(positions.size(), false);
  for (auto sample : sampling) visited[sample] = true;

  coords              = vector<vec2f>(solver.graph.size(), zero2f);
  coords[sampling[0]] = {radius, radius};
  vertices.insert(sampling[0]);
  for (int i = 1; i < sampling.size(); i++) {
    auto edge           = positions[sampling[i]] - positions[sampling[i - 1]];
    coords[sampling[i]] = {coords[sampling[i - 1]].x + length(edge), radius};
    vertices.insert(sampling[i]);
  }

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : sampling) distances[sample] = 0.0f;

  auto frames = vector<mat3f>(positions.size(), identity3x3f);
  compute_stroke_frames(frames, positions, normals, sampling);

  auto update = [&](int node, int neighbor, float new_distance) {
    vertices.insert(node);
    if (!visited[neighbor]) return;
    float weight = length(positions[neighbor] - positions[node]) + flt_eps;
    weight       = 1.0f / weight;
    compute_coordinates(coords, frames, positions, node, neighbor, weight);
    compute_frame(frames, normals, node, neighbor, weight);
    visited[node] = true;
  };
  dijkstra(solver, sampling, distances, radius, update);

  auto vec_vertices = vector<int>(vertices.begin(), vertices.end());

  // conversion in [0, 1]
  for (int i = 0; i < vertices.size(); i++)
    coords[vec_vertices[i]] /= radius * 2.0f;

  return vec_vertices;
}

// End stroke settings
void end_stroke(sculpt_params &params, sculpt_stroke &stroke,
    sculpt_buffers &buffers, scene_shape &shape) {
  stroke.lock = false;
  stroke.sampling.clear();
  buffers.coords.assign(buffers.coords.size(), {0, 0});
  buffers.old_positions = shape.positions;
  buffers.old_normals   = shape.normals;
}

// Compute gaussian function
float gaussian_distribution(vec3f origin, vec3f position, float standard_dev,
    float scale_factor, float strength, float radius) {
  float scaled_strength = strength / ((((radius - 0.1f) * 0.5f) / 0.7f) + 0.2f);
  float N               = 1.0f / (((standard_dev * scaled_strength) *
                        (standard_dev * scaled_strength) *
                        (standard_dev * scaled_strength)) *
                       sqrt((2.0f * pi) * (2.0f * pi) * (2.0f * pi)));
  float dx              = (origin.x - position.x) * scale_factor;
  float dy              = (origin.y - position.y) * scale_factor;
  float dz              = (origin.z - position.z) * scale_factor;
  float E               = ((dx * dx) + (dy * dy) + (dz * dz)) /
            (2.0f * standard_dev * standard_dev);
  return N * yocto::exp(-E);
}

// Change positions, normals, bounding volume hierarchy and hash grid
void apply_brush(shape_data &shape, const vector<vec3f> &positions,
    shape_bvh &tree, hash_grid &grid) {
  shape.positions = positions;
  triangles_normals(shape.normals, shape.triangles, shape.positions);
  update_triangles_bvh(tree, shape.triangles, shape.positions);
  grid = make_hash_grid(positions, grid.cell_size);
}

// To apply brush on intersected points' neighbors
bool gaussian_brush(vector<vec3f> &positions, const hash_grid &grid,
    const vector<pair<vec3f, vec3f>> &stroke_samples,
    const sculpt_params &             params) {
  if (stroke_samples.empty()) return false;

  // for a correct gaussian distribution
  float scale_factor = 3.5f / params.radius;

  auto neighbors = vector<int>{};
  for (auto [position, normal] : stroke_samples) {
    find_neighbors(grid, neighbors, position, params.radius);
    if (params.negative) normal = -normal;
    for (auto neighbor : neighbors) {
      auto gauss_height = gaussian_distribution(position, positions[neighbor],
          0.7f, scale_factor, params.strength, params.radius);
      positions[neighbor] += normal * gauss_height;
    }
    neighbors.clear();
  }

  return true;
}

// Compute texture values through the parameterization
bool texture_brush(vector<vec3f> &positions, const vector<int> &vertices,
    const image_data &texture, const vector<vec2f> &coords,
    const vector<vec3f> &base_positions, const vector<vec3f> &base_normals,
    const sculpt_params &params) {
  if (vertices.empty()) return false;
  if (texture.pixelsf.empty() && texture.pixelsb.empty()) return false;

  positions = base_positions;

  auto scale_factor = 3.5f / params.radius;
  auto max_height   = gaussian_distribution(
      zero3f, zero3f, 0.7f, scale_factor, params.strength, params.radius);

  for (auto idx : vertices) {
    auto uv     = coords[idx];
    auto height = max(xyz(eval_image(texture, uv)));
    auto normal = base_normals[idx];
    if (params.negative) normal = -normal;
    height *= max_height;
    positions[idx] += normal * height;
  }

  return true;
}

// Cotangent operator
float cotan(vec3f &a, vec3f &b) { return dot(a, b) / length(cross(a, b)); }

// Compute edge cotangents weights
float laplacian_weight(const vector<vec3f> &positions,
    const vector<vector<int>> &adjacencies, int node, int neighbor) {
  auto num_neighbors = int(adjacencies[node].size());

  int ind = -1;
  for (int i = 0; i < adjacencies[node].size(); i++) {
    if (adjacencies[node][i] == neighbor) {
      ind = i;
      break;
    }
  }

  auto prev = adjacencies[node][(ind - 1) >= 0 ? ind - 1 : num_neighbors - 1];
  auto next = adjacencies[node][(ind + 1) % num_neighbors];

  auto v1 = positions[node] - positions[prev];
  auto v2 = positions[neighbor] - positions[prev];
  auto v3 = positions[node] - positions[next];
  auto v4 = positions[neighbor] - positions[next];

  float cot_alpha = cotan(v1, v2);
  float cot_beta  = cotan(v3, v4);
  float weight    = (cot_alpha + cot_beta) / 2;

  float cotan_max = yocto::cos(flt_min) / yocto::sin(flt_min);
  weight          = clamp(weight, cotan_max, -cotan_max);
  return weight;
}

// Smooth brush with Laplace Operator Discretization and Cotangents Weights
bool smooth_brush(vector<vec3f> &positions, const geodesic_solver &solver,
    const vector<vector<int>> &adjacencies, vector<int> &stroke_sampling,
    const sculpt_params &params) {
  if (stroke_sampling.empty()) return false;

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : stroke_sampling) distances[sample] = 0.0f;

  auto current_node = -1;
  auto neighbors    = vector<int>{};
  auto weights      = vector<float>{};
  auto update       = [&](int node, int neighbor, float new_distance) {
    if (current_node == -1) current_node = node;
    if (node != current_node) {
      vec3f sum1 = zero3f;
      float sum2 = 0.0f;
      for (int i = 0; i < neighbors.size(); i++) {
        sum1 += positions[neighbors[i]] * weights[i];
        sum2 += weights[i];
      }
      positions[current_node] += 0.5f *
                                 (((sum1 / sum2) - positions[current_node]));
      current_node = node;
      neighbors.clear();
      weights.clear();
    }
    neighbors.push_back(neighbor);
    weights.push_back(laplacian_weight(positions, adjacencies, node, neighbor));
  };
  dijkstra(solver, stroke_sampling, distances, params.radius, update);

  stroke_sampling.clear();

  return true;
}

// TODO(fabio): move this function to math
static frame3f camera_frame(float lens, float aspect, float film = 0.036) {
  auto camera_dir  = normalize(vec3f{0, 0.5, 1});
  auto bbox_radius = 2.0f;
  auto camera_dist = bbox_radius * lens / (film / aspect);
  return lookat_frame(camera_dir * camera_dist, {0, 0, 0}, {0, 1, 0});
}

static void convert_scene(scene_scene &scene, const scene_shape &ioshape_,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 5};
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);

  // init scene
  scene = {};

  // rescale shape to unit
  auto ioshape = ioshape_;
  auto bbox    = invalidb3f;
  for (auto &pos : ioshape.positions) bbox = merge(bbox, pos);
  for (auto &pos : ioshape.positions) pos -= center(bbox);
  for (auto &pos : ioshape.positions) pos /= max(size(bbox));

  // camera
  if (progress_cb) progress_cb("create camera", progress.x++, progress.y);
  auto &camera  = scene.cameras.emplace_back();
  camera.frame  = camera_frame(0.050, 16.0f / 9.0f, 0.036);
  camera.lens   = 0.050;
  camera.aspect = 16.0f / 9.0f;
  camera.film   = 0.036;
  camera.focus  = length(camera.frame.o - center(bbox));

  // material
  if (progress_cb) progress_cb("create material", progress.x++, progress.y);
  auto &shape_material  = scene.materials.emplace_back();
  shape_material.type   = material_type::matte;
  shape_material.color  = {0.78f, 0.31f, 0.23f};
  auto &cursor_material = scene.materials.emplace_back();
  cursor_material.type  = material_type::matte;

  // shapes
  if (progress_cb) progress_cb("create shape", progress.x++, progress.y);
  scene.shapes.emplace_back(ioshape);
  scene.shapes.emplace_back(make_cursor({0, 0, 0}, {0, 0, 1}, 1));

  // instances
  if (progress_cb) progress_cb("create instance", progress.x++, progress.y);
  auto &shape_instance     = scene.instances.emplace_back();
  shape_instance.shape     = 0;
  shape_instance.material  = 0;
  auto &cursor_instance    = scene.instances.emplace_back();
  cursor_instance.shape    = 1;
  cursor_instance.material = 1;

  // done
  if (progress_cb) progress_cb("create scene", progress.x++, progress.y);
}

static void init_glscene(shade_scene &glscene, const sceneio_scene &ioscene,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene.cameras.size() + (int)ioscene.materials.size() +
             (int)ioscene.textures.size() + (int)ioscene.shapes.size() +
             (int)ioscene.instances.size()};

  // init scene
  init_scene(glscene);

  // camera
  for (auto &iocamera : ioscene.cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto &camera = glscene.cameras.at(add_camera(glscene));
    set_frame(camera, iocamera.frame);
    set_lens(camera, iocamera.lens, iocamera.aspect, iocamera.film);
    set_nearfar(camera, 0.001, 10000);
  }

  // textures
  for (auto &iotexture : ioscene.textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto  handle    = add_texture(glscene);
    auto &gltexture = glscene.textures[handle];
    if (!iotexture.pixelsf.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsf);
    } else if (!iotexture.pixelsb.empty()) {
      set_texture(
          gltexture, iotexture.width, iotexture.height, iotexture.pixelsb);
    }
  }

  // material
  for (auto &iomaterial : ioscene.materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto  handle     = add_material(glscene);
    auto &glmaterial = glscene.materials[handle];
    set_emission(glmaterial, iomaterial.emission, iomaterial.emission_tex);
    set_opacity(glmaterial, iomaterial.opacity, invalid_handle);
    set_normalmap(glmaterial, iomaterial.normal_tex);
    switch (iomaterial.type) {
      case material_type::matte: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, 0, invalid_handle);
      } break;
      case material_type::plastic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      case material_type::metal: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 1, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
      case material_type::metallic: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 1, invalid_handle);
        set_metallic(glmaterial, iomaterial.metallic, invalid_handle);
        set_roughness(glmaterial, iomaterial.roughness, invalid_handle);
      } break;
      default: {
        set_color(glmaterial, iomaterial.color, iomaterial.color_tex);
        set_specular(glmaterial, 0, invalid_handle);
        set_metallic(glmaterial, 0, invalid_handle);
        set_roughness(
            glmaterial, iomaterial.roughness, iomaterial.roughness_tex);
      } break;
    }
  }

  // shapes
  for (auto &ioshape : ioscene.shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    add_shape(glscene, ioshape.points, ioshape.lines, ioshape.triangles,
        ioshape.quads, ioshape.positions, ioshape.normals, ioshape.texcoords,
        ioshape.colors);
  }

  // shapes
  for (auto &ioinstance : ioscene.instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto  handle     = add_instance(glscene);
    auto &glinstance = glscene.instances[handle];
    set_frame(glinstance, ioinstance.frame);
    set_shape(glinstance, ioinstance.shape);
    set_material(glinstance, ioinstance.material);
  }

  // environments
  for (auto &ioenvironment : ioscene.environments) {
    auto  handle        = add_environment(glscene);
    auto &glenvironment = glscene.environments[handle];
    set_frame(glenvironment, ioenvironment.frame);
    set_emission(
        glenvironment, ioenvironment.emission, ioenvironment.emission_tex);
  }

  // init environments
  init_environments(glscene);

  // done
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);
}

// To make the stroke sampling (position, normal) following the mouse
bool sample_stroke(sculpt_stroke &stroke, const shape_bvh &bvh,
    const scene_shape &shape, const vec2f &mouse_uv, const scene_camera &camera,
    bool clear, const sculpt_params &params) {
  // clear
  if (clear) {
    stroke.pairs.clear();
    stroke.sampling.clear();
  }

  // helper
  auto intersect_shape = [&](const vec2f &uv) {
    auto ray = camera_ray(
        camera.frame, camera.lens, camera.aspect, camera.film, uv);
    return intersect_triangles_bvh(
        bvh, shape.triangles, shape.positions, ray, false);
  };

  // setup uv
  auto last_uv = stroke.locked_uv;

  // eval current intersection
  auto last  = intersect_shape(last_uv);
  auto mouse = intersect_shape(mouse_uv);
  if (!mouse.hit || !last.hit) return false;

  // sample
  auto delta_pos   = distance(eval_position(shape, last.element, last.uv),
      eval_position(shape, mouse.element, mouse.uv));
  auto stroke_dist = params.radius * 0.2f;
  auto delta_uv    = distance(mouse_uv, last_uv);
  auto steps       = int(delta_pos / stroke_dist);
  if (steps == 0) return true;
  auto stroke_uv = delta_uv * stroke_dist / delta_pos;
  auto mouse_dir = normalize(mouse_uv - stroke.locked_uv);
  for (auto step = 0; step < steps; step++) {
    last_uv += stroke_uv * mouse_dir;
    auto isec = intersect_shape(stroke.locked_uv);
    if (!isec.hit) continue;
    stroke.pairs.push_back({eval_position(shape, isec.element, isec.uv),
        eval_normal(shape, isec.element, isec.uv)});
    stroke.sampling.push_back(
        closest_vertex(shape.triangles, isec.element, isec.uv));
  }

  // update
  stroke.locked_uv = mouse_uv;

  return true;
}

pair<bool, bool> update_stroke(sculpt_stroke &stroke, sculpt_state &state,
    sculpt_params &params, sculpt_buffers &buffers, scene_shape &shape,
    scene_shape &cursor, const scene_camera &camera, const vec2f &mouse_uv,
    bool mouse_pressed) {
  auto updated_shape = false, updated_cursor = false;

  auto ray = camera_ray(
      camera.frame, camera.lens, camera.aspect, camera.film, mouse_uv);
  auto isec = intersect_triangles_bvh(
      state.bvh, shape.triangles, shape.positions, ray, false);
  if (isec.hit) {
    cursor         = make_cursor(eval_position(shape, isec.element, isec.uv),
        eval_normal(shape, isec.element, isec.uv),
        params.radius * (params.type == brush_type::gaussian ? 0.5f : 1.0f));
    updated_cursor = true;
  }

  // sculpting
  if (isec.hit && mouse_pressed) {
    if (!stroke.lock) {
      stroke.lock      = true;
      stroke.locked_uv = mouse_uv;
    } else {
      sample_stroke(stroke, state.bvh, shape, mouse_uv, camera,
          params.type != brush_type::texture, params);
      if (params.type == brush_type::gaussian) {
        updated_shape = gaussian_brush(
            shape.positions, state.grid, stroke.pairs, params);
      } else if (params.type == brush_type::smooth) {
        updated_shape = smooth_brush(shape.positions, state.solver,
            state.adjacencies, stroke.sampling, params);
      } else if (params.type == brush_type::texture && !stroke.pairs.empty() &&
                 !stroke.sampling.empty()) {
        auto vertices = stroke_parameterization(buffers.coords, state.solver,
            stroke.sampling, buffers.old_positions, buffers.old_normals,
            params.radius);
        updated_shape = texture_brush(shape.positions, vertices,
            state.tex_image, buffers.coords, buffers.old_positions,
            buffers.old_normals, params);
      }
    }
    if (updated_shape) {
      triangles_normals(shape.normals, shape.triangles, shape.positions);
      update_triangles_bvh(state.bvh, shape.triangles, shape.positions);
      state.grid = make_hash_grid(shape.positions, state.grid.cell_size);
    }
  } else {
    if (stroke.lock) end_stroke(params, stroke, buffers, shape);
  }

  return {updated_shape, updated_cursor};
}

int run_shade_sculpt(const shade_sculpt_params &params_) {
  // initialize app
  auto app = shade_sculpt_state();

  // copy command line
  app.filename    = params_.shape;
  app.texturename = params_.texture;
  app.imagename   = params_.texture;
  app.shapename   = params_.shape;

  // loading shape
  auto ioerror = ""s;
  auto ioshape = scene_shape{};
  print_progress("load shape", 0, 1);
  if (!load_shape(app.filename, ioshape, ioerror)) print_fatal(ioerror);
  if (!ioshape.quads.empty()) {
    ioshape.triangles = quads_to_triangles(ioshape.quads);
    ioshape.quads.clear();
  }
  print_progress("load shape", 1, 1);

  // loading texture
  auto iotexture = scene_texture{};
  if (!app.texturename.empty()) {
    print_progress("load texture", 0, 1);
    if (!load_texture(app.texturename, iotexture, ioerror))
      print_fatal(ioerror);
    print_progress("load texture", 1, 1);
  }

  // setup app
  convert_scene(app.ioscene, ioshape, print_progress);

  // sculpt params
  auto params  = sculpt_params{};
  auto stroke  = sculpt_stroke{};
  auto state   = make_sculpt_state(app.ioscene.shapes.front(), iotexture);
  auto buffers = make_sculpt_buffers(app.ioscene.shapes.front());

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [&app](gui_window *win, const gui_input &input) {
    init_glscene(app.glscene, app.ioscene, print_progress);
  };
  callbacks.clear_cb = [&app](gui_window *win, const gui_input &input) {
    clear_scene(app.glscene);
  };
  callbacks.draw_cb = [&app](gui_window *win, const gui_input &input) {
    draw_scene(app.glscene, app.glscene.cameras.at(0),
        input.framebuffer_viewport, app.drawgl_prms);
  };
  callbacks.widgets_cb = [&app, &params](
                             gui_window *win, const gui_input &input) {
    auto &glparams = app.drawgl_prms;
    draw_slider(win, "resolution", glparams.resolution, 0, 4096);
    draw_checkbox(win, "wireframe", glparams.wireframe);
    draw_combobox(
        win, "lighting", (int &)glparams.lighting, shade_lighting_names);
    continue_line(win);
    draw_checkbox(win, "double sided", glparams.double_sided);
    draw_slider(win, "exposure", glparams.exposure, -10, 10);
    draw_slider(win, "gamma", glparams.gamma, 0.1f, 4);
    draw_slider(win, "near", glparams.near, 0.01f, 1.0f);
    draw_slider(win, "far", glparams.far, 1000.0f, 10000.0f);
    draw_label(win, "", "");
    draw_label(win, "", "sculpt params");
    draw_combobox(win, "brush type", (int &)params.type, brushes_names);
    if (params.type == brush_type::gaussian) {
      if (params.strength < 0.8f || params.strength > 1.5f)
        params.strength = 1.0f;
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", params.negative);
    } else if (params.type == brush_type::texture) {
      if (params.strength < 0.8f || params.strength > 1.5f)
        params.strength = 1.0f;
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", params.negative);
    } else if (params.type == brush_type::smooth) {
      draw_slider(win, "radius", params.radius, 0.1f, 0.8f);
      draw_slider(win, "strength", params.strength, 0.1f, 1.0f);
    }
  };
  callbacks.update_cb = [](gui_window *win, const gui_input &input) {
    // pass
  };
  callbacks.uiupdate_cb = [&app, &params, &stroke, &buffers, &state](
                              gui_window *win, const gui_input &input) {
    // intersect mouse position and shape
    auto  mouse_uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
        input.mouse_pos.y / float(input.window_size.y)};
    auto &shape    = app.ioscene.shapes.at(0);
    auto &cursor   = app.ioscene.shapes.at(1);
    auto &camera   = app.ioscene.cameras.at(0);
    auto &glcamera = app.glscene.cameras.at(0);
    auto &glshape  = app.glscene.shapes.at(0);
    auto &glscene  = app.glscene;
    auto [updated_shape, updated_cursor] = update_stroke(stroke, state, params,
        buffers, shape, cursor, camera, mouse_uv,
        input.mouse_left && !input.modifier_ctrl);
    if (updated_cursor) {
      set_positions(glscene.shapes.at(1), cursor.positions);
      set_lines(glscene.shapes.at(1), cursor.lines);
      glscene.instances.at(1).hidden = false;
    } else {
      glscene.instances.at(1).hidden = true;
    }
    if (updated_shape) {
      set_positions(glshape, shape.positions);
      set_normals(glshape, shape.normals);
    }
    if (input.modifier_ctrl && !input.widgets_active) {
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      std::tie(camera.frame, camera.focus) = camera_turntable(
          camera.frame, camera.focus, rotate, dolly, -pan);
      glcamera.frame = camera.frame;
      glcamera.focus = camera.focus;
    }
  };

  // run ui
  run_ui({1900, 876}, "yshade", callbacks);

  // done
  return 0;
}
