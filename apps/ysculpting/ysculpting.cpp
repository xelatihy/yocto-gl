
#include <yocto/yocto_commonio.h>
#include <yocto/yocto_geometry.h>
#include <yocto/yocto_image.h>
#include <yocto/yocto_mesh.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto_gui/yocto_imgui.h>
#include <yocto_gui/yocto_shade.h>

#include <set>

using namespace yocto;

#include <atomic>
#include <deque>
#include <future>
#include <queue>
using namespace std::string_literals;

#include "ext/filesystem.hpp"
namespace sfs = ghc::filesystem;

#ifdef _WIN32
#undef near
#undef far
#endif

namespace yocto::sceneio {
void print_obj_camera(sceneio_camera *camera);
};

// Application state
struct app_state {
  // loading parameters
  string shapename = "";
  string filename  = "scene.json";
  string imagename = "";
  string outname   = "scene.json";
  string name      = "";

  // options
  shade_params drawgl_prms = {};

  // scene
  sceneio_scene * ioscene  = new sceneio_scene{};
  sceneio_camera *iocamera = nullptr;

  // rendering state
  shade_scene * glscene  = new shade_scene{};
  shade_camera *glcamera = nullptr;

  // editing
  sceneio_camera *     selected_camera      = nullptr;
  sceneio_instance *   selected_instance    = nullptr;
  sceneio_shape *      selected_shape       = nullptr;
  sceneio_material *   selected_material    = nullptr;
  sceneio_environment *selected_environment = nullptr;
  sceneio_texture *    selected_texture     = nullptr;

  // loading status
  std::atomic<bool> ok           = false;
  std::future<void> loader       = {};
  string            status       = "";
  string            error        = "";
  std::atomic<int>  current      = 0;
  std::atomic<int>  total        = 0;
  string            loader_error = "";

  ~app_state() {
    if (ioscene) delete ioscene;
    if (glscene) delete glscene;
    // pool->~ThreadPool();
  }
};

enum struct axes { x, y, z };
auto const symmetric_axes = vector<std::string>{"asse x", "asse y", "asse z"};

enum struct brush_type { gaussian, texture, smooth };
auto const brushes_names = vector<std::string>{
    "gaussian brush", "texture brush", "smooth brush"};

struct sculpt_params {
  // brush type
  brush_type type = brush_type::gaussian;

  // intersection
  ray3f              camera_ray       = {};
  sceneio_instance * shape_instance   = nullptr;
  shape_bvh          bvh_shape_tree   = {};
  shape_intersection bvh_intersection = {};

  // hash grid
  hash_grid hash_grid = {};

  // stroke
  vec3f locked_position = {};
  vec2f locked_uv       = {};
  bool  lock            = false;
  bool  continuous      = false;
  bool  symmetric       = false;
  axes  symmetric_axis  = axes::x;

  // brush
  float               radius      = 0.35f;
  float               strength    = 1.0f;
  bool                negative    = false;
  vector<float>       opacity     = {};
  bool                saturation  = false;
  vector<vector<int>> adjacencies = {};

  // stroke parameterization
  vector<int>            stroke_sampling           = {};
  vector<int>            symmetric_stroke_sampling = {};
  yocto::geodesic_solver solver                    = {};
  vector<vec2f>          coords                    = {};
  yocto::image<vec3f>    tex_image                 = {};
  vector<vec3f>          old_positions             = {};
  vector<vec3f>          old_normals               = {};
};

void init_glscene(shade_scene *glscene, sceneio_scene *ioscene,
    shade_camera *&glcamera, sceneio_camera *iocamera,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{
      0, (int)ioscene->cameras.size() + (int)ioscene->materials.size() +
             (int)ioscene->textures.size() + (int)ioscene->shapes.size() +
             (int)ioscene->instances.size()};

  // create scene
  init_scene(glscene);

  // camera
  auto camera_map     = unordered_map<sceneio_camera *, shade_camera *>{};
  camera_map[nullptr] = nullptr;
  for (auto iocamera : ioscene->cameras) {
    if (progress_cb) progress_cb("convert camera", progress.x++, progress.y);
    auto camera = add_camera(glscene);
    set_frame(camera, iocamera->frame);
    set_lens(camera, iocamera->lens, iocamera->aspect, iocamera->film);
    set_nearfar(camera, 0.001, 10000);
    camera_map[iocamera] = camera;
  }

  // textures
  auto texture_map     = unordered_map<sceneio_texture *, shade_texture *>{};
  texture_map[nullptr] = nullptr;
  for (auto iotexture : ioscene->textures) {
    if (progress_cb) progress_cb("convert texture", progress.x++, progress.y);
    auto gltexture = add_texture(glscene);
    if (!iotexture->hdr.empty()) {
      set_texture(gltexture, iotexture->hdr);
    } else if (!iotexture->ldr.empty()) {
      set_texture(gltexture, iotexture->ldr);
    }
    texture_map[iotexture] = gltexture;
  }

  // material
  auto material_map     = unordered_map<sceneio_material *, shade_material *>{};
  material_map[nullptr] = nullptr;
  for (auto iomaterial : ioscene->materials) {
    if (progress_cb) progress_cb("convert material", progress.x++, progress.y);
    auto glmaterial = add_material(glscene);
    set_emission(glmaterial, iomaterial->emission,
        texture_map.at(iomaterial->emission_tex));
    set_color(glmaterial, (1 - iomaterial->transmission) * iomaterial->color,
        texture_map.at(iomaterial->color_tex));
    set_specular(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->specular,
        texture_map.at(iomaterial->specular_tex));
    set_metallic(glmaterial,
        (1 - iomaterial->transmission) * iomaterial->metallic,
        texture_map.at(iomaterial->metallic_tex));
    set_roughness(glmaterial, iomaterial->roughness,
        texture_map.at(iomaterial->roughness_tex));
    set_opacity(glmaterial, iomaterial->opacity,
        texture_map.at(iomaterial->opacity_tex));
    set_normalmap(glmaterial, texture_map.at(iomaterial->normal_tex));
    material_map[iomaterial] = glmaterial;
  }

  // shapes
  auto shape_map     = unordered_map<sceneio_shape *, shade_shape *>{};
  shape_map[nullptr] = nullptr;
  for (auto ioshape : ioscene->shapes) {
    if (progress_cb) progress_cb("convert shape", progress.x++, progress.y);
    auto glshape = add_shape(glscene);
    set_positions(glshape, ioshape->positions);
    set_normals(glshape, ioshape->normals);
    set_texcoords(glshape, ioshape->texcoords);
    set_colors(glshape, ioshape->colors);
    set_points(glshape, ioshape->points);
    set_lines(glshape, ioshape->lines);
    set_triangles(glshape, ioshape->triangles);
    set_quads(glshape, ioshape->quads);
    shape_map[ioshape] = glshape;
  }

  // shapes
  for (auto ioobject : ioscene->instances) {
    if (progress_cb) progress_cb("convert instance", progress.x++, progress.y);
    auto globject = add_instance(glscene);
    set_frame(globject, ioobject->frame);
    set_shape(globject, shape_map.at(ioobject->shape));
    set_material(globject, material_map.at(ioobject->material));
  }

  // done
  if (progress_cb) progress_cb("convert done", progress.x++, progress.y);

  // get cmmera
  glcamera = camera_map.at(iocamera);
}

// Initialize all sculpting parameters.
sculpt_params *init_sculpt_tool(yocto::sceneio_scene *ioscene,
    yocto::sceneio_camera *iocamera, std::string shapename,
    std::string texture_name, sceneio_material *material = nullptr) {
  sculpt_params *params = new sculpt_params();

  // create sceneio_shape and sceneio_instance
  sceneio_instance *shape_instance = add_instance(ioscene, "object");
  sceneio_shape *   shape          = add_shape(ioscene, "object");

  // loading shape
  auto ioerror = ""s;
  if (!load_shape(shapename, shape->points, shape->lines, shape->triangles,
          shape->quads, shape->quadspos, shape->quadsnorm, shape->quadstexcoord,
          shape->positions, shape->normals, shape->texcoords, shape->colors,
          shape->radius, ioerror))
    print_fatal(ioerror);

  // convert quads meshes to triangles meshes
  // (geodesic solver works only on triangles)
  if (!shape->quads.empty()) {
    shape->triangles = quads_to_triangles(shape->quads);
    shape->quads.clear();
  }

  // set shape instance
  shape_instance->shape  = shape;
  shape_instance->frame  = identity3x4f;
  params->shape_instance = shape_instance;

  // save positions
  params->old_positions = shape_instance->shape->positions;

  // set/create sceneio_material
  if (material == nullptr) {
    material                 = add_material(ioscene, "clay");
    material->color          = {0.78f, 0.31f, 0.23f};
    material->specular       = 0.0f;
    shape_instance->material = material;
  } else {
    shape_instance->material = material;
  }
  shape_instance->shape->colors = vector<vec4f>(
      shape->positions.size(), vec4f{1.0f, 1.0f, 1.0f, 1.0f});
  // sceneio_texture* tex = add_texture(app->ioscene);

  // create camera
  add_cameras(ioscene);

  // create bvh structure
  params->bvh_shape_tree = make_triangles_bvh(shape_instance->shape->triangles,
      shape_instance->shape->positions, shape_instance->shape->radius);

  // create an hash grid
  params->hash_grid = make_hash_grid(shape_instance->shape->positions, 0.05f);

  // init saturation buffer
  params->opacity = vector<float>(
      shape_instance->shape->positions.size(), 0.0f);

  // create geodesic distance graph (ONLY TRIANGLES MESHES!)
  auto adjacencies    = face_adjacencies(shape_instance->shape->triangles);
  params->solver      = make_geodesic_solver(shape_instance->shape->triangles,
      adjacencies, shape_instance->shape->positions);
  params->adjacencies = vertex_adjacencies(
      shape_instance->shape->triangles, adjacencies);

  // init texture
  if (texture_name != "") {
    auto img = image<vec4f>{};
    if (!load_image(texture_name, img, ioerror)) print_fatal(ioerror);
    params->tex_image.resize(img.imsize());
    for (auto idx = 0; idx < img.count(); idx++)
      params->tex_image[idx] = xyz(img[idx]);
  }

  return params;
}

lines_shape make_circle(vec3f center, mat3f basis, float radius, int steps) {
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
void view_pointer(sceneio_instance *instance, shade_shape *glshape,
    shape_intersection intersection, float radius, int definition,
    brush_type &type) {
  if (intersection.hit) {
    if (type == brush_type::gaussian) radius *= 0.5f;
    auto pos   = eval_position(instance, intersection.element, intersection.uv);
    auto nor   = eval_normal(instance, intersection.element, intersection.uv);
    auto basis = basis_fromz(nor);
    auto circle    = make_circle(pos, basis, radius, definition);
    auto positions = instance->shape->positions;
    positions.insert(
        positions.end(), circle.positions.begin(), circle.positions.end());
    auto lines = instance->shape->lines;
    for (auto &line : circle.lines)
      line += int(instance->shape->positions.size());
    lines.insert(lines.end(), circle.lines.begin(), circle.lines.end());
    auto nor_line         = make_lines({1, 1});
    nor_line.positions[0] = pos;
    nor_line.positions[1] = pos + nor * 0.05f;
    positions.insert(
        positions.end(), nor_line.positions.begin(), nor_line.positions.end());
    lines.push_back({int(positions.size() - 2), int(positions.size() - 1)});
    set_positions(glshape, positions);
    set_lines(glshape, lines);
  } else {
    set_positions(glshape, instance->shape->positions);
    set_lines(glshape, instance->shape->lines);
  }
}

// Taking closest vertex of an intersection
int closest_vertex(sceneio_shape *shape, vec2f uv, int element) {
  auto tr = shape->triangles[element];
  if (uv.x < 0.5f && uv.y < 0.5f) return tr.x;
  if (uv.x > uv.y) return tr.y;
  return tr.z;
}

// To make the stroke sampling (position, normal) following the mouse
vector<pair<vec3f, vec3f>> stroke(
    sculpt_params *params, vec2f mouse_uv, shade_camera *glcamera) {
  // eval current intersection
  vector<pair<vec3f, vec3f>> pairs;
  auto &                     inter = params->bvh_intersection;
  auto  pos = eval_position(params->shape_instance, inter.element, inter.uv);
  auto  nor = eval_normal(params->shape_instance, inter.element, inter.uv);
  float delta_pos   = distance(pos, params->locked_position);
  float stroke_dist = params->radius * 0.2f;

  // handle continuous stroke
  if (params->continuous && (delta_pos < stroke_dist)) {
    pair<vec3f, vec3f> pair;
    pair.first  = pos;
    pair.second = nor;
    pairs.push_back(pair);
    params->stroke_sampling.push_back(
        closest_vertex(params->shape_instance->shape, inter.uv, inter.element));
    return pairs;
  }

  // handle first stroke intersection
  if (!params->lock) {
    params->locked_position = pos;
    params->locked_uv       = mouse_uv;
    params->lock            = true;
    return pairs;
  }

  float delta_uv = distance(mouse_uv, params->locked_uv);
  int   steps    = int(delta_pos / stroke_dist);
  if (steps == 0) return pairs;
  pairs           = vector<pair<vec3f, vec3f>>(steps);
  float stroke_uv = delta_uv * stroke_dist / delta_pos;
  auto  mouse_dir = normalize(mouse_uv - params->locked_uv);
  for (int step = 0; step < steps; step++) {
    params->locked_uv += stroke_uv * mouse_dir;
    auto ray = camera_ray(glcamera->frame, glcamera->lens, glcamera->aspect,
        glcamera->film, params->locked_uv);
    inter    = intersect_triangles_bvh(params->bvh_shape_tree,
        params->shape_instance->shape->triangles,
        params->shape_instance->shape->positions, ray, false);
    pos      = eval_position(params->shape_instance, inter.element, inter.uv);
    nor      = eval_normal(params->shape_instance, inter.element, inter.uv);
    params->stroke_sampling.push_back(
        closest_vertex(params->shape_instance->shape, inter.uv, inter.element));
    auto pair               = std::pair<vec3f, vec3f>{pos, nor};
    params->locked_position = pos;
    pairs.push_back(pair);
  }
  return pairs;
}

// To obtain symmetric from stroke result
vector<pair<vec3f, vec3f>> symmetric_stroke(vector<pair<vec3f, vec3f>> &pairs,
    sceneio_instance *instance, shape_bvh &tree,
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
        tree, instance->shape->triangles, instance->shape->positions, ray);
    if (!inter.hit) continue;
    auto pos  = eval_position(instance, inter.element, inter.uv);
    auto nor  = eval_normal(instance, inter.element, inter.uv);
    auto pair = std::pair<vec3f, vec3f>{pos, nor};
    symmetric_stroke_sampling.push_back(
        closest_vertex(instance->shape, inter.uv, inter.element));
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
void compute_coordinates(vector<vec2f> &coords, vector<mat3f> &frames,
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
void compute_frame(vector<mat3f> &frames, vector<vec3f> &normals, int node,
    int neighbor, float weight) {
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
void compute_stroke_frames(vector<mat3f> &frames, vector<vec3f> &positions,
    vector<vec3f> &normals, vector<int> &stroke_sampling) {
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
vector<int> stroke_parameterization(yocto::geodesic_solver &solver,
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
        coords[stroke_sampling[i - 1]].x + yocto::length(edge), radius};
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

// End stroke settings
void end_stroke(sculpt_params *params) {
  params->lock = false;
  std::fill(params->opacity.begin(), params->opacity.end(), 0.0f);
  params->stroke_sampling.clear();
  params->coords.clear();
  params->symmetric_stroke_sampling.clear();
  params->old_positions = params->shape_instance->shape->positions;
  params->old_normals   = params->shape_instance->shape->normals;
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

// Change positions, normals, boundig volume hierarchy and hash grid
void apply_brush(sceneio_shape *shape, vector<vec3f> &positions,
    shade_shape *glshape, shape_bvh &tree, hash_grid &grid) {
  shape->positions = positions;
  set_positions(glshape, positions);
  update_normals(shape->normals, shape->triangles, shape->positions);
  set_normals(glshape, shape->normals);
  update_triangles_bvh(tree, shape->triangles, shape->positions);
  grid = make_hash_grid(positions, grid.cell_size);
}

// To apply brush on intersected points' neighbors
void brush(sculpt_params *params, shade_shape *glshape,
    vector<pair<vec3f, vec3f>> &pairs) {
  if (pairs.empty()) return;
  auto &positions = params->shape_instance->shape->positions;
  auto  neighbors = vector<int>{};

  // for a correct gaussian distribution
  float scale_factor = 3.5f / params->radius;

  for (int p = 0; p < pairs.size(); p++) {
    find_neighbors(
        params->hash_grid, neighbors, pairs[p].first, params->radius);
    auto center = pairs[p].first;
    auto normal = pairs[p].second;
    if (params->negative) normal = -normal;
    if (params->saturation) {
      float max_height = gaussian_distribution(
          center, center, 0.7f, scale_factor, params->strength, params->radius);
      for (int v = 0; v < neighbors.size(); v++) {
        float opacity_val = params->opacity[neighbors[v]];
        if (opacity_val == 1.0f) continue;
        float gauss_height = gaussian_distribution(center,
            positions[neighbors[v]], 0.7f, scale_factor, params->strength,
            params->radius);
        float gauss_ratio  = gauss_height / max_height;
        if ((opacity_val + gauss_ratio) > 1.0f) {
          gauss_height = max_height * (1.0f - opacity_val);
          gauss_ratio  = 1.0f - opacity_val;
        }
        params->opacity[neighbors[v]] += gauss_ratio;
        positions[neighbors[v]] += normal * gauss_height;
      }
    } else {
      for (int v = 0; v < neighbors.size(); v++) {
        float gauss_height = gaussian_distribution(center,
            positions[neighbors[v]], 0.7f, scale_factor, params->strength,
            params->radius);
        positions[neighbors[v]] += normal * gauss_height;
      }
    }
    neighbors.clear();
  }

  apply_brush(params->shape_instance->shape, positions, glshape,
      params->bvh_shape_tree, params->hash_grid);
}

// Compute texture values through the parameterization
void texture_brush(vector<int> &vertices, yocto::image<vec3f> &texture,
    vector<vec2f> &coords, sculpt_params *params, shade_shape *glshape,
    vector<vec3f> positions, vector<vec3f> normals) {
  if (vertices.empty()) return;
  if (texture.empty()) return;

  auto scale_factor = 3.5f / params->radius;
  auto max_height   = gaussian_distribution(
      zero3f, zero3f, 0.7f, scale_factor, params->strength, params->radius);

  for (auto i : vertices) {
    auto uv     = coords[i];
    auto color  = eval_image(texture, uv);
    auto height = yocto::length(color);
    auto normal = normals[i];
    if (params->negative) normal = -normal;
    height *= max_height;
    positions[i] += normal * height;
  }

  apply_brush(params->shape_instance->shape, positions, glshape,
      params->bvh_shape_tree, params->hash_grid);
}

// Cotangent operator
float cotan(vec3f &a, vec3f &b) {
  return dot(a, b) / yocto::length(cross(a, b));
}

// Compute edge cotangents weights
float laplacian_weight(vector<vec3f> &positions,
    vector<vector<int>> &adjacencies, int node, int neighbor) {
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
void smooth(geodesic_solver &solver, vector<int> &stroke_sampling,
    vector<vec3f> &positions, sculpt_params *params, shade_shape *glshape) {
  if (stroke_sampling.empty()) return;

  auto distances = vector<float>(solver.graph.size(), flt_max);
  for (auto sample : stroke_sampling) distances[sample] = 0.0f;

  int  current_node = -1;
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
    weights.push_back(
        laplacian_weight(positions, params->adjacencies, node, neighbor));
  };
  dijkstra(solver, stroke_sampling, distances, params->radius, update);

  apply_brush(params->shape_instance->shape, positions, glshape,
      params->bvh_shape_tree, params->hash_grid);

  stroke_sampling.clear();
}

int main(int argc, const char *argv[]) {
  // initialize app
  auto app_guard = std::make_unique<app_state>();
  auto app       = app_guard.get();
  // auto camera_name = "default";

  // parse command line
  auto cli = make_cli("ysculpting", " sculpt a mesh interactively");
  add_option(
      cli, "--resolution,-r", app->drawgl_prms.resolution, "Image resolution.");
  add_option(cli, "--lighting", app->drawgl_prms.lighting, "Lighting mode.",
      shade_lighting_names);
  add_option(cli, "shape", app->shapename, "Shape filename.", true);
  add_option(cli, "texture", app->imagename, "Texture filename.");
  parse_cli(cli, argc, argv);

  auto ioerror              = ""s;
  app->drawgl_prms.lighting = shade_lighting_type::eyelight;

  sculpt_params *sculpt_params = init_sculpt_tool(
      app->ioscene, app->iocamera, app->shapename, app->imagename);

  app->iocamera = get_camera(app->ioscene);

  // tesselation
  tesselate_shapes(app->ioscene, print_progress);

  // callbacks
  auto callbacks    = gui_callbacks{};
  callbacks.init_cb = [app](gui_window *win, const gui_input &input) {
    init_glscene(app->glscene, app->ioscene, app->glcamera, app->iocamera,
        [app](const string &message, int current, int total) {
          app->status  = "init scene";
          app->current = current;
          app->total   = total;
        });
  };
  callbacks.clear_cb = [app](gui_window *win, const gui_input &input) {
    clear_scene(app->glscene);
  };
  callbacks.draw_cb = [app](gui_window *win, const gui_input &input) {
    draw_scene(app->glscene, app->glcamera, input.framebuffer_viewport,
        app->drawgl_prms);
  };
  callbacks.widgets_cb = [app, &sculpt_params](
                             gui_window *win, const gui_input &input) {
    draw_progressbar(win, app->status.c_str(), app->current, app->total);
    if (draw_combobox(win, "camera", app->iocamera, app->ioscene->cameras)) {
      for (auto idx = 0; idx < app->ioscene->cameras.size(); idx++) {
        if (app->ioscene->cameras[idx] == app->iocamera)
          app->glcamera = app->glscene->cameras[idx];
      }
    }
    auto &params = app->drawgl_prms;
    // float mouse_x = input.mouse_pos.x;
    // float mouse_y = input.mouse_pos.y;
    draw_slider(win, "resolution", params.resolution, 0, 4096);
    draw_checkbox(win, "wireframe", params.wireframe);
    draw_combobox(
        win, "lighting", (int &)params.lighting, shade_lighting_names);
    continue_line(win);
    draw_checkbox(win, "double sided", params.double_sided);
    draw_slider(win, "exposure", params.exposure, -10, 10);
    draw_slider(win, "gamma", params.gamma, 0.1f, 4);
    draw_slider(win, "near", params.near, 0.01f, 1.0f);
    draw_slider(win, "far", params.far, 1000.0f, 10000.0f);
    draw_label(win, "", "");
    draw_label(win, "", "sculpt params");
    draw_combobox(win, "brush type", (int &)sculpt_params->type, brushes_names);
    if (sculpt_params->type == brush_type::gaussian) {
      if (sculpt_params->strength < 0.8f || sculpt_params->strength > 1.5f)
        sculpt_params->strength = 1.0f;
      draw_slider(win, "radius", sculpt_params->radius, 0.1f, 0.8f);
      draw_slider(win, "strength", sculpt_params->strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", sculpt_params->negative);
      draw_checkbox(win, "continuous", sculpt_params->continuous);
      draw_checkbox(win, "saturation", sculpt_params->saturation);
      draw_checkbox(win, "symmetric", sculpt_params->symmetric);
      if (sculpt_params->symmetric)
        draw_combobox(win, "symmetric axis",
            (int &)sculpt_params->symmetric_axis, symmetric_axes);
    } else if (sculpt_params->type == brush_type::texture) {
      if (sculpt_params->strength < 0.8f || sculpt_params->strength > 1.5f)
        sculpt_params->strength = 1.0f;
      draw_slider(win, "radius", sculpt_params->radius, 0.1f, 0.8f);
      draw_slider(win, "strength", sculpt_params->strength, 1.5f, 0.9f);
      draw_checkbox(win, "negative", sculpt_params->negative);
      draw_checkbox(win, "symmetric", sculpt_params->symmetric);
      if (sculpt_params->symmetric)
        draw_combobox(win, "symmetric axis",
            (int &)sculpt_params->symmetric_axis, symmetric_axes);
      /*
      capire come funziona
      draw_filedialog_button(win, "Load", true, "Load Texture", app->imagename,
          false, "./", "", "");
          */
    } else if (sculpt_params->type == brush_type::smooth) {
      draw_slider(win, "radius", sculpt_params->radius, 0.1f, 0.8f);
      draw_slider(win, "strength", sculpt_params->strength, 0.1f, 1.0f);
      draw_checkbox(win, "symmetric", sculpt_params->symmetric);
      if (sculpt_params->symmetric)
        draw_combobox(win, "symmetric axis",
            (int &)sculpt_params->symmetric_axis, symmetric_axes);
    }
  };
  callbacks.update_cb = [](gui_window *win, const gui_input &input) {
    // update(win, apps);
  };
  callbacks.uiupdate_cb = [app, &sculpt_params](
                              gui_window *win, const gui_input &input) {
    // handle mouse and keyboard for navigation
    // if (input.modifier_alt) printf("%s\n", app->imagename);

    // intersect mouse position and shape
    vec2f mouse_uv = vec2f{input.mouse_pos.x / float(input.window_size.x),
        input.mouse_pos.y / float(input.window_size.y)};

    sculpt_params->camera_ray       = camera_ray(app->glcamera->frame,
        app->glcamera->lens, app->glcamera->aspect, app->glcamera->film,
        mouse_uv);
    sculpt_params->bvh_intersection = intersect_triangles_bvh(
        sculpt_params->bvh_shape_tree,
        sculpt_params->shape_instance->shape->triangles,
        sculpt_params->shape_instance->shape->positions,
        sculpt_params->camera_ray, false);
    view_pointer(sculpt_params->shape_instance,
        app->glscene->instances[0]->shape, sculpt_params->bvh_intersection,
        sculpt_params->radius, 20, sculpt_params->type);

    // sculpting
    if (input.mouse_left && sculpt_params->bvh_intersection.hit &&
        !input.modifier_ctrl) {
      auto        pairs = stroke(sculpt_params, mouse_uv, app->glcamera);
      vector<int> vertices;
      if (sculpt_params->type == brush_type::gaussian) {
        brush(sculpt_params, app->glscene->instances[0]->shape, pairs);
      } else if (sculpt_params->type == brush_type::smooth) {
        smooth(sculpt_params->solver, sculpt_params->stroke_sampling,
            sculpt_params->shape_instance->shape->positions, sculpt_params,
            app->glscene->instances[0]->shape);
      } else if (sculpt_params->type == brush_type::texture && !pairs.empty()) {
        vertices = stroke_parameterization(sculpt_params->solver,
            sculpt_params->coords, sculpt_params->stroke_sampling,
            sculpt_params->old_positions, sculpt_params->old_normals,
            sculpt_params->radius);
        texture_brush(vertices, sculpt_params->tex_image, sculpt_params->coords,
            sculpt_params, app->glscene->instances[0]->shape,
            sculpt_params->old_positions, sculpt_params->old_normals);
      }
      if (sculpt_params->symmetric) {
        pairs = symmetric_stroke(pairs, sculpt_params->shape_instance,
            sculpt_params->bvh_shape_tree,
            sculpt_params->symmetric_stroke_sampling,
            sculpt_params->symmetric_axis);
        if (sculpt_params->type == brush_type::gaussian) {
          brush(sculpt_params, app->glscene->instances[0]->shape, pairs);
        } else if (sculpt_params->type == brush_type::smooth) {
          smooth(sculpt_params->solver,
              sculpt_params->symmetric_stroke_sampling,
              sculpt_params->shape_instance->shape->positions, sculpt_params,
              app->glscene->instances[0]->shape);
        } else if (sculpt_params->type == brush_type::texture &&
                   !pairs.empty()) {
          vertices = stroke_parameterization(sculpt_params->solver,
              sculpt_params->coords, sculpt_params->symmetric_stroke_sampling,
              sculpt_params->old_positions, sculpt_params->old_normals,
              sculpt_params->radius);
          texture_brush(vertices, sculpt_params->tex_image,
              sculpt_params->coords, sculpt_params,
              app->glscene->instances[0]->shape,
              sculpt_params->shape_instance->shape->positions,
              sculpt_params->shape_instance->shape->normals);
        }
      }
    } else {
      end_stroke(sculpt_params);
    }
    if (input.modifier_ctrl && !input.widgets_active) {
      // printf("%d\n", input.window_size.x);
      auto dolly  = 0.0f;
      auto pan    = zero2f;
      auto rotate = zero2f;
      if (input.mouse_left && !input.modifier_shift)
        rotate = (input.mouse_pos - input.mouse_last) / 100.0f;
      if (input.mouse_right)
        dolly = (input.mouse_pos.x - input.mouse_last.x) / 100.0f;
      if (input.mouse_left && input.modifier_shift)
        pan = (input.mouse_pos - input.mouse_last) / 100.0f;
      update_turntable(
          app->iocamera->frame, app->iocamera->focus, rotate, dolly, -pan);
      set_frame(app->glcamera, app->iocamera->frame);
    }
  };

  // run ui
  run_ui({1900, 876}, "ysculpting", callbacks);

  // done
  return 0;
}
