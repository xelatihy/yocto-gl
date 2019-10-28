//
// Implementation for Yocto/GL Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#include "yocto_sceneio.h"
#include "yocto_commonio.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_obj.h"
#include "yocto_ply.h"
#include "yocto_pbrt.h"
#include "yocto_shape.h"

#include <atomic>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <deque>
#include <future>

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, Func&& func) {
  auto             futures  = vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(begin);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, end]() {
          while (true) {
            auto idx = next_idx.fetch_add(1);
            if (idx >= end) break;
            func(idx);
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Find the first keyframe value that is greater than the argument.
inline int keyframe_index(const vector<float>& times, const float& time) {
  for (auto i = 0; i < times.size(); i++)
    if (times[i] > time) return i;
  return (int)times.size();
}

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T keyframe_step(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f keyframe_slerp(
    const vector<float>& times, const vector<vec4f>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline T keyframe_linear(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evaluates a keyframed value using Bezier interpolation.
template <typename T>
inline T keyframe_bezier(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return interpolate_bezier(
      vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_model& scene) {
  auto shape_bbox = vector<bbox3f>{};
  for (auto& shape : scene.shapes) {
    auto& sbvh = shape_bbox.emplace_back();
    sbvh       = invalidb3f;
    for (auto p : shape.positions) sbvh = merge(sbvh, p);
  }
  auto bbox = invalidb3f;
  for (auto& instance : scene.instances) {
    bbox = merge(
        bbox, transform_bbox(instance.frame, shape_bbox[instance.shape]));
  }
  return bbox;
}

// Add missing cameras.
void add_cameras(scene_model& scene) {
  if (!scene.cameras.empty()) return;
  auto& camera = scene.cameras.emplace_back();
  camera.name  = "default";
  // TODO: fix me
  // FIXME: error in camera.lens and camera.film
  camera.orthographic = false;
  camera.film         = 0.036;
  camera.aperture     = 0;
  camera.lens         = 0.050;
  auto bbox           = compute_bounds(scene);
  auto center         = (bbox.max + bbox.min) / 2;
  auto bbox_radius    = length(bbox.max - bbox.min) / 2;
  auto camera_dir     = camera.frame.o - center;
  if (camera_dir == zero3f) camera_dir = {0, 0, 1};
  auto camera_dist = bbox_radius / camera.film;
  auto from        = camera_dir * (camera_dist * 1) + center;
  auto to          = center;
  auto up          = vec3f{0, 1, 0};
  camera.frame     = lookat_frame(from, to, up);
  camera.focus     = length(from - to);
}

// Add missing materials.
void add_materials(scene_model& scene) {
  auto material_id = -1;
  for (auto& instance : scene.instances) {
    if (instance.material >= 0) continue;
    if (material_id < 0) {
      auto material    = scene_material{};
      material.name    = "default";
      material.diffuse = {0.2f, 0.2f, 0.2f};
      scene.materials.push_back(material);
      material_id = (int)scene.materials.size() - 1;
    }
    instance.material = material_id;
  }
}

// Add missing radius.
void add_radius(scene_model& scene, float radius = 0.001f) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add a sky environment
void add_sky(scene_model& scene, float sun_angle) {
  auto texture     = scene_texture{};
  texture.name     = "sky";
  texture.filename = "textures/sky.hdr";
  make_sunsky(texture.hdr, {1024, 512}, sun_angle);
  scene.textures.push_back(texture);
  auto environment         = scene_environment{};
  environment.name         = "sky";
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
  scene.environments.push_back(environment);
}

vector<string> format_stats(const scene_model& scene, bool verbose) {
  auto accumulate = [](const auto& values, const auto& func) -> size_t {
    auto sum = (size_t)0;
    for (auto& value : values) sum += func(value);
    return sum;
  };
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

  auto bbox = compute_bounds(scene);

  auto stats = vector<string>{};
  stats.push_back("cameras:      " + format(scene.cameras.size()));
  stats.push_back("shapes:       " + format(scene.shapes.size()));
  stats.push_back("instances:    " + format(scene.instances.size()));
  stats.push_back("environments: " + format(scene.environments.size()));
  stats.push_back("textures:     " + format(scene.textures.size()));
  stats.push_back("materials:    " + format(scene.materials.size()));
  stats.push_back("nodes:        " + format(scene.nodes.size()));
  stats.push_back("animations:   " + format(scene.animations.size()));
  stats.push_back(
      "points:       " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.points.size(); })));
  stats.push_back(
      "lines:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.lines.size(); })));
  stats.push_back("triangles:    " +
                  format(accumulate(scene.shapes,
                      [](auto& shape) { return shape.triangles.size(); })));
  stats.push_back(
      "quads:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.quads.size(); })));
  stats.push_back("fvquads:      " +
                  format(accumulate(scene.shapes,
                      [](auto& shape) { return shape.quadspos.size(); })));
  stats.push_back(
      "texels4b:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.ldr.size().x * (size_t)texture.ldr.size().x;
      })));
  stats.push_back(
      "texels4f:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.hdr.size().x * (size_t)texture.hdr.size().y;
      })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Reduce memory usage
void trim_memory(scene_model& scene) {
  for (auto& shape : scene.shapes) {
    shape.points.shrink_to_fit();
    shape.lines.shrink_to_fit();
    shape.triangles.shrink_to_fit();
    shape.quads.shrink_to_fit();
    shape.quadspos.shrink_to_fit();
    shape.quadsnorm.shrink_to_fit();
    shape.quadstexcoord.shrink_to_fit();
    shape.positions.shrink_to_fit();
    shape.normals.shrink_to_fit();
    shape.texcoords.shrink_to_fit();
    shape.colors.shrink_to_fit();
    shape.radius.shrink_to_fit();
    shape.tangents.shrink_to_fit();
  }
  for (auto& texture : scene.textures) {
    texture.ldr.shrink_to_fit();
    texture.hdr.shrink_to_fit();
  }
  scene.cameras.shrink_to_fit();
  scene.shapes.shrink_to_fit();
  scene.instances.shrink_to_fit();
  scene.materials.shrink_to_fit();
  scene.textures.shrink_to_fit();
  scene.environments.shrink_to_fit();
  scene.nodes.shrink_to_fit();
  scene.animations.shrink_to_fit();
}

// Checks for validity of the scene.
vector<string> format_validation(const scene_model& scene, bool notextures) {
  auto errs        = vector<string>();
  auto check_names = [&errs](const auto& vals, const string& base) {
    auto used = unordered_map<string, int>();
    used.reserve(vals.size());
    for (auto& value : vals) used[value.name] += 1;
    for (auto& [name, used] : used) {
      if (name == "") {
        errs.push_back("empty " + base + " name");
      } else if (used > 1) {
        errs.push_back("duplicated " + base + " name " + name);
      }
    }
  };
  auto check_empty_textures = [&errs](const vector<scene_texture>& vals) {
    for (auto& value : vals) {
      if (value.hdr.empty() && value.ldr.empty()) {
        errs.push_back("empty texture " + value.name);
      }
    }
  };

  check_names(scene.cameras, "camera");
  check_names(scene.shapes, "shape");
  check_names(scene.textures, "texture");
  check_names(scene.materials, "material");
  check_names(scene.instances, "instance");
  check_names(scene.environments, "environment");
  check_names(scene.nodes, "node");
  check_names(scene.animations, "animation");
  if (!notextures) check_empty_textures(scene.textures);

  return errs;
}

// Apply subdivision and displacement rules.
scene_shape subdivide_shape(const scene_shape& shape) {
  if (!shape.subdivisions) return shape;
  auto subdiv         = shape;
  subdiv.subdivisions = 0;
  if (!shape.points.empty()) {
    throw std::runtime_error("point subdivision not supported");
  } else if (!shape.lines.empty()) {
    subdivide_lines(subdiv.lines, subdiv.positions, subdiv.normals,
        subdiv.texcoords, subdiv.colors, subdiv.radius, shape.lines,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius, shape.subdivisions);
    if (shape.smooth)
      subdiv.normals = compute_tangents(subdiv.lines, subdiv.positions);
  } else if (!shape.triangles.empty()) {
    subdivide_triangles(subdiv.triangles, subdiv.positions, subdiv.normals,
        subdiv.texcoords, subdiv.colors, subdiv.radius, shape.triangles,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius, shape.subdivisions);
    if (shape.smooth)
      subdiv.normals = compute_normals(subdiv.triangles, subdiv.positions);
  } else if (!shape.quads.empty() && !shape.catmullclark) {
    subdivide_quads(subdiv.quads, subdiv.positions, subdiv.normals,
        subdiv.texcoords, subdiv.colors, subdiv.radius, shape.quads,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius, shape.subdivisions);
    if (subdiv.smooth)
      subdiv.normals = compute_normals(subdiv.quads, subdiv.positions);
  } else if (!shape.quads.empty() && shape.catmullclark) {
    subdivide_catmullclark(subdiv.quads, subdiv.positions, subdiv.normals,
        subdiv.texcoords, subdiv.colors, subdiv.radius, shape.quads,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius, shape.subdivisions);
    if (subdiv.smooth)
      subdiv.normals = compute_normals(subdiv.quads, subdiv.positions);
  } else if (!shape.quadspos.empty() && !shape.catmullclark) {
    subdivide_quads(subdiv.quadspos, subdiv.positions, subdiv.quadspos,
        shape.positions, shape.subdivisions);
    subdivide_quads(subdiv.quadsnorm, subdiv.normals, subdiv.quadsnorm,
        shape.normals, shape.subdivisions);
    subdivide_quads(subdiv.quadstexcoord, subdiv.texcoords,
        subdiv.quadstexcoord, shape.texcoords, shape.subdivisions);
    if (subdiv.smooth) {
      subdiv.normals   = compute_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    }
  } else if (!shape.quadspos.empty() && shape.catmullclark) {
    subdivide_catmullclark(subdiv.quadspos, subdiv.positions, shape.quadspos,
        shape.positions, subdiv.subdivisions);
    subdivide_catmullclark(subdiv.quadstexcoord, subdiv.texcoords,
        shape.quadstexcoord, shape.texcoords, shape.subdivisions, true);
    if (shape.smooth) {
      subdiv.normals   = compute_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    }
  } else {
    throw std::runtime_error("empty shape");
  }
  return subdiv;
}
// Apply displacement to a shape
scene_shape displace_shape(const scene_model& scene, const scene_shape& shape) {
  // Evaluate a texture
  auto eval_texture = [](const scene_texture& texture, const vec2f& texcoord,
                          bool ldr_as_linear, bool no_interpolation = false,
                          bool clamp_to_edge = false) -> vec4f {
    if (!texture.hdr.empty()) {
      return eval_image(texture.hdr, texcoord, no_interpolation, clamp_to_edge);
    } else if (!texture.ldr.empty()) {
      return eval_image(texture.ldr, texcoord, ldr_as_linear, no_interpolation,
          clamp_to_edge);
    } else {
      return {1, 1, 1, 1};
    }
  };

  if (!shape.displacement || shape.displacement_tex < 0) return shape;
  auto& displacement = scene.textures[shape.displacement_tex];
  if (shape.texcoords.empty()) {
    throw std::runtime_error("missing texture coordinates");
    return {};
  }

  auto subdiv             = shape;
  subdiv.displacement     = 0;
  subdiv.displacement_tex = -1;

  // simple case
  if (!shape.triangles.empty()) {
    auto normals = shape.normals;
    if (normals.empty())
      normals = compute_normals(shape.triangles, shape.positions);
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
      auto disp = mean(
          xyz(eval_texture(displacement, shape.texcoords[vid], true)));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      subdiv.positions[vid] += normals[vid] * shape.displacement * disp;
    }
    if (shape.smooth || !shape.normals.empty())
      subdiv.normals = compute_normals(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    auto normals = shape.normals;
    if (normals.empty())
      normals = compute_normals(shape.triangles, shape.positions);
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
      auto disp = mean(
          xyz(eval_texture(displacement, shape.texcoords[vid], true)));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      subdiv.positions[vid] += shape.normals[vid] * shape.displacement * disp;
    }
    if (shape.smooth || !shape.normals.empty())
      subdiv.normals = compute_normals(shape.quads, shape.positions);
  } else if (!shape.quadspos.empty()) {
    // facevarying case
    auto offset = vector<float>(shape.positions.size(), 0);
    auto count  = vector<int>(shape.positions.size(), 0);
    for (auto fid = 0; fid < shape.quadspos.size(); fid++) {
      auto qpos = shape.quadspos[fid];
      auto qtxt = shape.quadstexcoord[fid];
      for (auto i = 0; i < 4; i++) {
        auto disp = mean(
            xyz(eval_texture(displacement, shape.texcoords[qtxt[i]], true)));
        if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
        offset[qpos[i]] += shape.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = vector<vec3f>{shape.positions.size()};
    compute_normals(normals, shape.quadspos, shape.positions);
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
      subdiv.positions[vid] += normals[vid] * offset[vid] / count[vid];
    }
    if (shape.smooth || !shape.normals.empty()) {
      subdiv.quadsnorm = shape.quadspos;
      subdiv.normals   = compute_normals(shape.quadspos, shape.positions);
    }
  }
  return subdiv;
}

// Update animation transforms
void update_transforms(scene_model& scene, scene_animation& animation,
    float time, const string& anim_group) {
  if (anim_group != "" && anim_group != animation.group) return;

  if (!animation.translations.empty()) {
    auto value = vec3f{0, 0, 0};
    switch (animation.interpolation) {
      case scene_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.translations, time);
        break;
      case scene_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.translations, time);
        break;
      case scene_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.translations, time);
        break;
      default: throw std::runtime_error("should not have been here");
    }
    for (auto target : animation.targets)
      scene.nodes[target].translation = value;
  }
  if (!animation.rotations.empty()) {
    auto value = vec4f{0, 0, 0, 1};
    switch (animation.interpolation) {
      case scene_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.rotations, time);
        break;
      case scene_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.rotations, time);
        break;
      case scene_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.rotations, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].rotation = value;
  }
  if (!animation.scales.empty()) {
    auto value = vec3f{1, 1, 1};
    switch (animation.interpolation) {
      case scene_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.scales, time);
        break;
      case scene_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.scales, time);
        break;
      case scene_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.scales, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].scale = value;
  }
}

// Update node transforms
void update_transforms(scene_model& scene, scene_node& node,
    const frame3f& parent = identity3x4f) {
  auto frame = parent * node.local * translation_frame(node.translation) *
               rotation_frame(node.rotation) * scaling_frame(node.scale);
  if (node.instance >= 0) scene.instances[node.instance].frame = frame;
  if (node.camera >= 0) scene.cameras[node.camera].frame = frame;
  if (node.environment >= 0) scene.environments[node.environment].frame = frame;
  for (auto child : node.children)
    update_transforms(scene, scene.nodes[child], frame);
}

// Update node transforms
void update_transforms(
    scene_model& scene, float time, const string& anim_group) {
  for (auto& agr : scene.animations)
    update_transforms(scene, agr, time, anim_group);
  for (auto& node : scene.nodes) node.children.clear();
  for (auto node_id = 0; node_id < scene.nodes.size(); node_id++) {
    auto& node = scene.nodes[node_id];
    if (node.parent >= 0) scene.nodes[node.parent].children.push_back(node_id);
  }
  for (auto& node : scene.nodes)
    if (node.parent < 0) update_transforms(scene, node);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin YAML format.
static void load_yaml_scene(
    const string& filename, scene_model& scene, const load_params& params);
static void save_yaml_scene(const string& filename, const scene_model& scene,
    const save_params& params);

// Load/save a scene from/to OBJ.
static void load_obj_scene(
    const string& filename, scene_model& scene, const load_params& params);
static void save_obj_scene(const string& filename, const scene_model& scene,
    const save_params& params);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static void load_ply_scene(
    const string& filename, scene_model& scene, const load_params& params);
static void save_ply_scene(const string& filename, const scene_model& scene,
    const save_params& params);

// Load/save a scene from/to glTF.
static void load_gltf_scene(
    const string& filename, scene_model& scene, const load_params& params);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static void load_pbrt_scene(
    const string& filename, scene_model& scene, const load_params& params);
static void save_pbrt_scene(const string& filename, const scene_model& scene,
    const save_params& params);

// Load a scene
void load_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  try {
    auto ext = get_extension(filename);
    if (ext == ".yaml" || ext == ".YAML") {
      load_yaml_scene(filename, scene, params);
    } else if (ext == ".obj" || ext == ".OBJ") {
      load_obj_scene(filename, scene, params);
    } else if (ext == ".gltf" || ext == ".GLTF") {
      load_gltf_scene(filename, scene, params);
    } else if (ext == ".pbrt" || ext == ".PBRT") {
      load_pbrt_scene(filename, scene, params);
    } else if (ext == ".ply" || ext == ".PLY") {
      load_ply_scene(filename, scene, params);
    } else {
      scene = {};
      throw std::runtime_error("unsupported scene format " + ext);
    }
  } catch (std::exception& e) {
    throw std::runtime_error("cannot load scene " + filename + "\n" + e.what());
  }
}

// Save a scene
void save_scene(const string& filename, const scene_model& scene,
    const save_params& params) {
  try {
    auto ext = get_extension(filename);
    if (ext == ".yaml" || ext == ".YAML") {
      save_yaml_scene(filename, scene, params);
    } else if (ext == ".obj" || ext == ".OBJ") {
      save_obj_scene(filename, scene, params);
    } else if (ext == ".pbrt" || ext == ".PBRT") {
      save_pbrt_scene(filename, scene, params);
    } else if (ext == ".ply" || ext == ".PLY") {
      save_ply_scene(filename, scene, params);
    } else {
      throw std::runtime_error("unsupported scene format " + ext);
    }
  } catch (std::exception& e) {
    throw std::runtime_error("cannot load scene " + filename + "\n" + e.what());
  }
}

void load_texture(scene_texture& texture, const string& dirname) {
  if (is_hdr_filename(texture.filename)) {
    load_image(dirname + texture.filename, texture.hdr);
  } else {
    load_imageb(dirname + texture.filename, texture.ldr);
  }
}

void load_textures(
    scene_model& scene, const string& dirname, const load_params& params) {
  if (params.notextures) return;

  // load images
  if (params.noparallel) {
    for (auto& texture : scene.textures) {
      if (!texture.hdr.empty() || !texture.ldr.empty()) return;
      load_texture(texture, dirname);
    }
  } else {
    parallel_foreach(scene.textures, [&dirname](scene_texture& texture) {
      if (!texture.hdr.empty() || !texture.ldr.empty()) return;
      load_texture(texture, dirname);
    });
  }
}

void save_texture(const scene_texture& texture, const string& dirname) {
  if (!texture.hdr.empty()) {
    save_image(dirname + texture.filename, texture.hdr);
  } else {
    save_imageb(dirname + texture.filename, texture.ldr);
  }
}

// helper to save textures
void save_textures(const scene_model& scene, const string& dirname,
    const save_params& params) {
  if (params.notextures) return;

  // save images
  if (params.noparallel) {
    for (auto& texture : scene.textures) {
      save_texture(texture, dirname);
    }
  } else {
    parallel_foreach(scene.textures, [&dirname](const scene_texture& texture) {
      save_texture(texture, dirname);
    });
  }
}

// Load json meshes
void load_shapes(
    scene_model& scene, const string& dirname, const load_params& params) {
  // load shapes
  if (params.noparallel) {
    for (auto& shape : scene.shapes) {
      if (!shape.positions.empty()) continue;
      if (!shape.facevarying) {
        load_shape(dirname + shape.filename, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius);
      } else {
        load_fvshape(dirname + shape.filename, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords);
      }
    }
  } else {
    parallel_foreach(scene.shapes, [&dirname](scene_shape& shape) {
      if (!shape.positions.empty()) return;
      if (!shape.facevarying) {
        load_shape(dirname + shape.filename, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius);
      } else {
        load_fvshape(dirname + shape.filename, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords);
      }
    });
  }
}

// Save json meshes
void save_shapes(const scene_model& scene, const string& dirname,
    const save_params& params) {
  // save shapes
  if (params.noparallel) {
    for (auto& shape : scene.shapes) {
      if (shape.quadspos.empty()) {
        save_shape(dirname + shape.filename, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius);
      } else {
        save_fvshape(dirname + shape.filename, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords);
      }
    }
  } else {
    parallel_foreach(scene.shapes, [&dirname](const scene_shape& shape) {
      if (shape.quadspos.empty()) {
        save_shape(dirname + shape.filename, shape.points, shape.lines,
            shape.triangles, shape.quads, shape.positions, shape.normals,
            shape.texcoords, shape.colors, shape.radius);
      } else {
        save_fvshape(dirname + shape.filename, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords);
      }
    });
  }
}

// create and cleanup names and filenames
static string make_safe_name(
    const string& name_, const string& base, int count) {
  auto name = name_;
  if (name.empty()) name = base + std::to_string(count);
  if (name.front() == '-') name = "_" + name;
  if (name.front() >= '0' && name.front() <= '9') name = "_" + name;
  for (auto& c : name) {
    if (c == '-' || c == '_') continue;
    if (c >= '0' && c <= '9') continue;
    if (c >= 'a' && c <= 'z') continue;
    if (c >= 'A' && c <= 'Z') continue;
    c = '_';
  }
  std::transform(name.begin(), name.end(), name.begin(),
      [](unsigned char c) { return std::tolower(c); });
  return name;
}
static inline string make_safe_filename(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename) {
    if (c == ' ') c = '_';
  }
  return filename;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL YAML DECLARATIONS
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;
using namespace std::literals::string_view_literals;

// Yaml value type
enum struct yaml_value_type { number, boolean, string, array };

// Yaml value
struct yaml_value {
  yaml_value_type   type    = yaml_value_type::number;
  double            number  = 0;
  bool              boolean = false;
  string            string_ = "";
  array<double, 16> array_  = {};
};

// Yaml element
struct yaml_element {
  string                           name       = "";
  vector<pair<string, yaml_value>> key_values = {};
};

// Yaml model
struct yaml_model {
  vector<string>       comments = {};
  vector<yaml_element> elements = {};
};

// Load/save yaml
void load_yaml(const string& filename, yaml_model& yaml);
void save_yaml(const string& filename, const yaml_model& yaml);

// Load Yaml properties
bool read_yaml_property(file_wrapper& fs, string& group, string& key,
    bool& newobj, yaml_value& value);

// Write Yaml properties
void write_yaml_comment(file_wrapper& fs, const string& comment);
void write_yaml_property(file_wrapper& fs, const string& object,
    const string& key, bool newobj, const yaml_value& value);
void write_yaml_object(file_wrapper& fs, const string& object);

// type-cheked yaml value access
void get_yaml_value(const yaml_value& yaml, string& value);
void get_yaml_value(const yaml_value& yaml, bool& value);
void get_yaml_value(const yaml_value& yaml, int& value);
void get_yaml_value(const yaml_value& yaml, float& value);
void get_yaml_value(const yaml_value& yaml, vec2f& value);
void get_yaml_value(const yaml_value& yaml, vec3f& value);
void get_yaml_value(const yaml_value& yaml, mat3f& value);
void get_yaml_value(const yaml_value& yaml, frame3f& value);

// yaml value construction
yaml_value make_yaml_value(const string& value);
yaml_value make_yaml_value(bool value);
yaml_value make_yaml_value(int value);
yaml_value make_yaml_value(float value);
yaml_value make_yaml_value(const vec2f& value);
yaml_value make_yaml_value(const vec3f& value);
yaml_value make_yaml_value(const mat3f& value);
yaml_value make_yaml_value(const frame3f& value);

}

// -----------------------------------------------------------------------------
// LOW-LEVEL YAML IMPLEMENTATION
// -----------------------------------------------------------------------------
namespace yocto {

// Read a line
static inline bool read_yaml_line(file_wrapper& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

static inline bool is_yaml_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_yaml_newline(char c) { return c == '\r' || c == '\n'; }
static inline bool is_yaml_digit(char c) { return c >= '0' && c <= '9'; }
static inline bool is_yaml_alpha(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

static inline void skip_yaml_whitespace(string_view& str) {
  while (!str.empty() && is_yaml_space(str.front())) str.remove_prefix(1);
}
static inline void trim_yaml_whitespace(string_view& str) {
  while (!str.empty() && is_yaml_space(str.front())) str.remove_prefix(1);
  while (!str.empty() && is_yaml_space(str.back())) str.remove_suffix(1);
}

static inline bool is_yaml_whitespace(string_view str) {
  while (!str.empty()) {
    if (!is_yaml_space(str.front())) return false;
    str.remove_prefix(1);
  }
  return true;
}

static inline void remove_yaml_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_yaml_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

static inline void parse_yaml_varname(string_view& str, string_view& value) {
  skip_yaml_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (!is_yaml_alpha(str.front())) throw std::runtime_error("cannot parse value");
  auto pos = 0;
  while (is_yaml_alpha(str[pos]) || str[pos] == '_' || is_yaml_digit(str[pos])) {
    pos += 1;
    if (pos >= str.size()) break;
  }
  value = str.substr(0, pos);
  str.remove_prefix(pos);
}
static inline void parse_yaml_varname(string_view& str, string& value) {
  auto view = ""sv;
  parse_yaml_varname(str, view);
  value = string{view};
}

inline void parse_yaml_value(string_view& str, string_view& value) {
  skip_yaml_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_yaml_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
    if (str.front() != '"') throw std::runtime_error("cannot parse value");
    str.remove_prefix(1);
    if (str.empty()) throw std::runtime_error("cannot parse value");
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
    if (cpy.empty()) throw std::runtime_error("cannot parse value");
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    str.remove_prefix(1);
  }
}
inline void parse_yaml_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_yaml_value(str, valuev);
  value = string{valuev};
}
inline void parse_yaml_value(string_view& str, int& value) {
  skip_yaml_whitespace(str);
  char* end = nullptr;
  value     = (int)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_yaml_value(string_view& str, float& value) {
  skip_yaml_whitespace(str);
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_yaml_value(string_view& str, double& value) {
  skip_yaml_whitespace(str);
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}

// parse yaml value
void get_yaml_value(const yaml_value& yaml, string& value) {
  if (yaml.type != yaml_value_type::string)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.string_;
}
void get_yaml_value(const yaml_value& yaml, bool& value) {
  if (yaml.type != yaml_value_type::boolean)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.boolean;
}
void get_yaml_value(const yaml_value& yaml, int& value) {
  if (yaml.type != yaml_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (int)yaml.number;
}
void get_yaml_value(const yaml_value& yaml, float& value) {
  if (yaml.type != yaml_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (float)yaml.number;
}
void get_yaml_value(const yaml_value& yaml, vec2f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 2)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1]};
}
void get_yaml_value(const yaml_value& yaml, vec3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 3)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1], (float)yaml.array_[2]};
}
void get_yaml_value(const yaml_value& yaml, mat3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 9)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 9; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}
void get_yaml_value(const yaml_value& yaml, frame3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 12)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 12; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}

// construction
yaml_value make_yaml_value(const string& value) {
  return {yaml_value_type::string, 0, false, value};
}
yaml_value make_yaml_value(bool value) {
  return {yaml_value_type::boolean, 0, value};
}
yaml_value make_yaml_value(int value) {
  return {yaml_value_type::number, (double)value};
}
yaml_value make_yaml_value(float value) {
  return {yaml_value_type::number, (double)value};
}
yaml_value make_yaml_value(const vec2f& value) {
  return {
      yaml_value_type::array, 2, false, "", {(double)value.x, (double)value.y}};
}
yaml_value make_yaml_value(const vec3f& value) {
  return {yaml_value_type::array, 3, false, "",
      {(double)value.x, (double)value.y, (double)value.z}};
}
yaml_value make_yaml_value(const mat3f& value) {
  auto yaml = yaml_value{yaml_value_type::array, 9};
  for (auto i = 0; i < 9; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
}
yaml_value make_yaml_value(const frame3f& value) {
  auto yaml = yaml_value{yaml_value_type::array, 12};
  for (auto i = 0; i < 12; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
}

void parse_yaml_value(string_view& str, yaml_value& value) {
  trim_yaml_whitespace(str);
  if (str.empty()) throw std::runtime_error("bad yaml");
  if (str.front() == '[') {
    str.remove_prefix(1);
    value.type   = yaml_value_type::array;
    value.number = 0;
    while (!str.empty()) {
      skip_yaml_whitespace(str);
      if (str.empty()) throw std::runtime_error("bad yaml");
      if (str.front() == ']') {
        str.remove_prefix(1);
        break;
      }
      if (value.number >= 16) throw std::runtime_error("array too large");
      parse_yaml_value(str, value.array_[(int)value.number]);
      value.number += 1;
      skip_yaml_whitespace(str);
      if (str.front() == ',') {
        str.remove_prefix(1);
        continue;
      } else if (str.front() == ']') {
        str.remove_prefix(1);
        break;
      } else {
        throw std::runtime_error("bad yaml");
      }
    }
  } else if (is_yaml_digit(str.front()) || str.front() == '-' ||
             str.front() == '+') {
    value.type = yaml_value_type::number;
    parse_yaml_value(str, value.number);
  } else {
    value.type = yaml_value_type::string;
    parse_yaml_value(str, value.string_);
    if (value.string_ == "true" || value.string_ == "false") {
      value.type    = yaml_value_type::boolean;
      value.boolean = value.string_ == "true";
    }
  }
  skip_yaml_whitespace(str);
  if (!str.empty() && !is_yaml_whitespace(str)) throw std::runtime_error("bad yaml");
}

// Load/save yaml
void load_yaml(const string& filename, yaml_model& yaml) {
  auto fs = open_file(filename, "rt");

  // read the file line by line
  auto group = ""s;
  auto key   = ""s;
  auto value = yaml_value{};
  char buffer[4096];
  while (read_yaml_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_yaml_comment(line);
    if (line.empty()) continue;
    if (is_yaml_whitespace(line)) continue;

    // peek commands
    if (is_yaml_space(line.front())) {
      // indented property
      if (group == "") throw std::runtime_error("bad yaml");
      skip_yaml_whitespace(line);
      if (line.empty()) throw std::runtime_error("bad yaml");
      if (line.front() == '-') {
        auto& element = yaml.elements.emplace_back();
        element.name  = group;
        line.remove_prefix(1);
        skip_yaml_whitespace(line);
      } else if (yaml.elements.empty() || yaml.elements.back().name != group) {
        auto& element = yaml.elements.emplace_back();
        element.name  = group;
      }
      parse_yaml_varname(line, key);
      skip_yaml_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      parse_yaml_value(line, value);
      yaml.elements.back().key_values.push_back({key, value});
    } else if (is_yaml_alpha(line.front())) {
      // new group
      parse_yaml_varname(line, key);
      skip_yaml_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      if (!line.empty() && !is_yaml_whitespace(line)) {
        group = "";
        if (yaml.elements.empty() || yaml.elements.back().name != group) {
          auto& element = yaml.elements.emplace_back();
          element.name  = group;
        }
        parse_yaml_value(line, value);
        yaml.elements.back().key_values.push_back({key, value});
      } else {
        group = key;
        key   = "";
      }
    } else {
      throw std::runtime_error("bad yaml");
    }
  }
}

static inline void format_yaml_value(string& str, const string& value) {
  str += value;
}
static inline void format_yaml_value(string& str, const char* value) {
  str += value;
}
static inline void format_yaml_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}

static inline void format_yaml_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static inline void format_yaml_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_yaml_value(str, arg);
  format_yaml_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
static inline void format_yaml_values(
    file_wrapper& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_yaml_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
static inline void format_yaml_value(file_wrapper& fs, const T& value) {
  auto str = ""s;
  format_yaml_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

static inline void format_yaml_value(string& str, const yaml_value& value) {
  switch (value.type) {
    case yaml_value_type::number: format_yaml_value(str, value.number); break;
    case yaml_value_type::boolean:
      format_yaml_value(str, value.boolean ? "true" : "false");
      break;
    case yaml_value_type::string:
      if (value.string_.empty() || is_yaml_digit(value.string_.front())) {
        format_yaml_values(str, "\"{}\"", value.string_);
      } else {
        format_yaml_values(str, "{}", value.string_);
      }
      break;
    case yaml_value_type::array:
      format_yaml_value(str, "[ ");
      for (auto i = 0; i < value.number; i++) {
        if (i) format_yaml_value(str, ", ");
        format_yaml_value(str, value.array_[i]);
      }
      format_yaml_value(str, " ]");
      break;
  }
}

void save_yaml(const string& filename, const yaml_model& yaml) {
  auto fs = open_file(filename, "wt");

  // save comments
  format_yaml_values(fs, "#\n");
  format_yaml_values(fs, "# Written by Yocto/GL\n");
  format_yaml_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_yaml_values(fs, "#\n\n");
  for (auto& comment : yaml.comments) {
    format_yaml_values(fs, "# {}\n", comment);
  }
  format_yaml_values(fs, "\n");

  auto group = ""s;
  for (auto& element : yaml.elements) {
    if (group != element.name) {
      group = element.name;
      if (group != "") {
        format_yaml_values(fs, "\n{}:\n", group);
      } else {
        format_yaml_values(fs, "\n");
      }
      auto first = true;
      for (auto& [key, value] : element.key_values) {
        if (group != "") {
          format_yaml_values(fs, "  {} {}: {}\n", first ? "-" : " ", key, value);
          first = false;
        } else {
          format_yaml_values(fs, "{}: {}\n", key, value);
        }
      }
    }
  }
}

bool read_yaml_property(file_wrapper& fs, string& group, string& key,
    bool& newobj, yaml_value& value) {
  // read the file line by line
  char buffer[4096];
  while (read_yaml_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_yaml_comment(line);
    if (line.empty()) continue;
    if (is_yaml_whitespace(line)) continue;

    // peek commands
    if (is_yaml_space(line.front())) {
      // indented property
      if (group == "") throw std::runtime_error("bad yaml");
      skip_yaml_whitespace(line);
      if (line.empty()) throw std::runtime_error("bad yaml");
      if (line.front() == '-') {
        newobj = true;
        line.remove_prefix(1);
        skip_yaml_whitespace(line);
      } else {
        newobj = false;
      }
      parse_yaml_varname(line, key);
      skip_yaml_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      parse_yaml_value(line, value);
      return true;
    } else if (is_yaml_alpha(line.front())) {
      // new group
      parse_yaml_varname(line, key);
      skip_yaml_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      if (!line.empty() && !is_yaml_whitespace(line)) {
        group = "";
        parse_yaml_value(line, value);
        return true;
      } else {
        group = key;
        key   = "";
        return true;
      }
    } else {
      throw std::runtime_error("bad yaml");
    }
  }
  return false;
}

static inline vector<string> split_yaml_string(
    const string& str, const string& delim) {
  auto tokens = vector<string>{};
  auto last = (size_t)0, next = (size_t)0;
  while ((next = str.find(delim, last)) != string::npos) {
    tokens.push_back(str.substr(last, next - last));
    last = next + delim.size();
  }
  if (last < str.size()) tokens.push_back(str.substr(last));
  return tokens;
}

void write_yaml_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_yaml_string(comment, "\n");
  for (auto& line : lines) {
    format_yaml_values(fs, "# {}\n", line);
  }
  format_yaml_values(fs, "\n");
}

// Save yaml property
void write_yaml_property(file_wrapper& fs, const string& object,
    const string& key, bool newobj, const yaml_value& value) {
  if (key.empty()) {
    format_yaml_values(fs, "\n{}:\n", object);
  } else {
    if (!object.empty()) {
      format_yaml_values(fs, "  {} {}: {}\n", newobj ? "-" : " ", key, value);
    } else {
      format_yaml_values(fs, "{}: {}\n", key, value);
    }
  }
}

void write_yaml_object(file_wrapper& fs, const string& object) {
  format_yaml_values(fs, "\n{}:\n", object);
}

}

// -----------------------------------------------------------------------------
// YAML SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

void load_yaml(
    const string& filename, scene_model& scene, const load_params& params) {
  // open file
  auto fs = open_file(filename);

  // parse state
  enum struct parsing_type {
    // clang-format off
    none, camera, texture, material, shape, instance, environment
    // clang-format on
  };
  auto type = parsing_type::none;

  auto tmap = unordered_map<string, int>{{"", -1}};
  auto vmap = unordered_map<string, int>{{"", -1}};
  auto mmap = unordered_map<string, int>{{"", -1}};
  auto smap = unordered_map<string, int>{{"", -1}};

  // parse yaml reference
  auto get_yaml_ref = [](const yaml_value& yaml, int& value,
                          const unordered_map<string, int>& refs) {
    if (yaml.type != yaml_value_type::string)
      throw std::runtime_error("error parsing yaml value");
    if (yaml.string_ == "") return;
    try {
      value = refs.at(yaml.string_);
    } catch (...) {
      throw std::runtime_error("reference not found " + yaml.string_);
    }
  };

  // load yaml
  auto group  = ""s;
  auto key    = ""s;
  auto newobj = false;
  auto value  = yaml_value{};
  while (read_yaml_property(fs, group, key, newobj, value)) {
    if (group.empty()) {
      throw std::runtime_error("bad yaml");
    }
    if (key.empty()) {
      type = parsing_type::none;
      continue;
    }
    if (newobj) {
      if (group == "cameras") {
        type = parsing_type::camera;
        scene.cameras.push_back({});
      } else if (group == "textures") {
        type = parsing_type::texture;
        scene.textures.push_back({});
      } else if (group == "materials") {
        type = parsing_type::material;
        scene.materials.push_back({});
      } else if (group == "shapes") {
        type = parsing_type::shape;
        scene.shapes.push_back({});
      } else if (group == "instances") {
        type = parsing_type::instance;
        scene.instances.push_back({});
      } else if (group == "environments") {
        type = parsing_type::environment;
        scene.environments.push_back({});
      } else {
        type = parsing_type::none;
        throw std::runtime_error("unknown object type " + string(group));
      }
    }
    if (type == parsing_type::none) {
      throw std::runtime_error("bad yaml");
    } else if (type == parsing_type::camera) {
      auto& camera = scene.cameras.back();
      if (key == "name") {
        get_yaml_value(value, camera.name);
      } else if (key == "uri") {
        get_yaml_value(value, camera.name);
        camera.name = get_basename(camera.name);
      } else if (key == "frame") {
        get_yaml_value(value, camera.frame);
      } else if (key == "orthographic") {
        get_yaml_value(value, camera.orthographic);
      } else if (key == "lens") {
        get_yaml_value(value, camera.lens);
      } else if (key == "aspect") {
        get_yaml_value(value, camera.aspect);
      } else if (key == "film") {
        get_yaml_value(value, camera.film);
      } else if (key == "focus") {
        get_yaml_value(value, camera.focus);
      } else if (key == "aperture") {
        get_yaml_value(value, camera.aperture);
      } else if (key == "lookat") {
        auto lookat = identity3x3f;
        get_yaml_value(value, lookat);
        camera.frame = lookat_frame(lookat.x, lookat.y, lookat.z);
        camera.focus = length(lookat.x - lookat.y);
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else if (type == parsing_type::texture) {
      auto& texture = scene.textures.back();
      if (key == "name") {
        get_yaml_value(value, texture.name);
        tmap[texture.name] = (int)scene.textures.size() - 1;
      } else if (key == "filename") {
        get_yaml_value(value, texture.filename);
      } else if (key == "preset") {
        auto preset = ""s;
        get_yaml_value(value, preset);
        make_image_preset(texture.hdr, texture.ldr, preset);
        if (texture.filename.empty()) {
          texture.filename = "textures/ypreset-" + preset +
                             (texture.hdr.empty() ? ".png" : ".hdr");
        }
      } else if (key == "uri") {
        get_yaml_value(value, texture.filename);
        texture.name           = get_basename(texture.filename);
        tmap[texture.filename] = (int)scene.textures.size() - 1;
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else if (type == parsing_type::material) {
      auto& material = scene.materials.back();
      if (key == "name") {
        get_yaml_value(value, material.name);
        mmap[material.name] = (int)scene.materials.size() - 1;
      } else if (key == "uri") {
        get_yaml_value(value, material.name);
        mmap[material.name] = (int)scene.materials.size() - 1;
        material.name       = get_basename(material.name);
      } else if (key == "emission") {
        get_yaml_value(value, material.emission);
      } else if (key == "diffuse") {
        get_yaml_value(value, material.diffuse);
      } else if (key == "metallic") {
        get_yaml_value(value, material.metallic);
      } else if (key == "specular") {
        get_yaml_value(value, material.specular);
      } else if (key == "roughness") {
        get_yaml_value(value, material.roughness);
      } else if (key == "coat") {
        get_yaml_value(value, material.coat);
      } else if (key == "transmission") {
        get_yaml_value(value, material.transmission);
      } else if (key == "refract") {
        get_yaml_value(value, material.refract);
      } else if (key == "voltransmission") {
        get_yaml_value(value, material.voltransmission);
      } else if (key == "volmeanfreepath") {
        get_yaml_value(value, material.volmeanfreepath);
      } else if (key == "volscatter") {
        get_yaml_value(value, material.volscatter);
      } else if (key == "volemission") {
        get_yaml_value(value, material.volemission);
      } else if (key == "volanisotropy") {
        get_yaml_value(value, material.volanisotropy);
      } else if (key == "volscale") {
        get_yaml_value(value, material.volscale);
      } else if (key == "opacity") {
        get_yaml_value(value, material.opacity);
      } else if (key == "coat") {
        get_yaml_value(value, material.coat);
      } else if (key == "emission_tex") {
        get_yaml_ref(value, material.emission_tex, tmap);
      } else if (key == "diffuse_tex") {
        get_yaml_ref(value, material.diffuse_tex, tmap);
      } else if (key == "metallic_tex") {
        get_yaml_ref(value, material.metallic_tex, tmap);
      } else if (key == "specular_tex") {
        get_yaml_ref(value, material.specular_tex, tmap);
      } else if (key == "transmission_tex") {
        get_yaml_ref(value, material.transmission_tex, tmap);
      } else if (key == "roughness_tex") {
        get_yaml_ref(value, material.roughness_tex, tmap);
      } else if (key == "subsurface_tex") {
        get_yaml_ref(value, material.subsurface_tex, tmap);
      } else if (key == "opacity_tex") {
        get_yaml_ref(value, material.normal_tex, tmap);
      } else if (key == "normal_tex") {
        get_yaml_ref(value, material.normal_tex, tmap);
      } else if (key == "gltf_textures") {
        get_yaml_value(value, material.gltf_textures);
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else if (type == parsing_type::shape) {
      auto& shape = scene.shapes.back();
      if (key == "name") {
        get_yaml_value(value, shape.name);
        smap[shape.name] = (int)scene.shapes.size() - 1;
      } else if (key == "filename") {
        get_yaml_value(value, shape.filename);
      } else if (key == "preset") {
        auto preset = ""s;
        get_yaml_value(value, preset);
        make_shape_preset(shape.points, shape.lines, shape.triangles,
            shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords, shape.colors,
            shape.radius, preset);
        if (shape.filename.empty()) {
          shape.filename = "shapes/ypreset-" + preset + ".yvol";
        }
      } else if (key == "subdivisions") {
        get_yaml_value(value, shape.subdivisions);
      } else if (key == "catmullclark") {
        get_yaml_value(value, shape.catmullclark);
      } else if (key == "smooth") {
        get_yaml_value(value, shape.smooth);
      } else if (key == "facevarying") {
        get_yaml_value(value, shape.facevarying);
      } else if (key == "displacement_tex") {
        get_yaml_ref(value, shape.displacement_tex, tmap);
      } else if (key == "displacement") {
        get_yaml_value(value, shape.displacement);
      } else if (key == "uri") {
        get_yaml_value(value, shape.filename);
        shape.name           = get_basename(shape.filename);
        smap[shape.filename] = (int)scene.shapes.size() - 1;
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else if (type == parsing_type::instance) {
      auto& instance = scene.instances.back();
      if (key == "name") {
        get_yaml_value(value, instance.name);
      } else if (key == "uri") {
        get_yaml_value(value, instance.name);
      } else if (key == "frame") {
        get_yaml_value(value, instance.frame);
      } else if (key == "shape") {
        get_yaml_ref(value, instance.shape, smap);
      } else if (key == "material") {
        get_yaml_ref(value, instance.material, mmap);
      } else if (key == "lookat") {
        auto lookat = identity3x3f;
        get_yaml_value(value, lookat);
        instance.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else if (type == parsing_type::environment) {
      auto& environment = scene.environments.back();
      if (key == "name") {
        get_yaml_value(value, environment.name);
      } else if (key == "uri") {
        get_yaml_value(value, environment.name);
      } else if (key == "frame") {
        get_yaml_value(value, environment.frame);
      } else if (key == "emission") {
        get_yaml_value(value, environment.emission);
      } else if (key == "emission_tex") {
        get_yaml_ref(value, environment.emission_tex, tmap);
      } else if (key == "lookat") {
        auto lookat = identity3x3f;
        get_yaml_value(value, lookat);
        environment.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      } else {
        throw std::runtime_error("unknown property " + string(key));
      }
    } else {
      assert(false);  // should not get here
    }
  }
}

// Save a scene in the builtin YAML format.
static void load_yaml_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  scene = {};

  // Parse yaml
  load_yaml(filename, scene, params);

  // load shape and textures
  auto dirname = get_dirname(filename);
  load_shapes(scene, dirname, params);
  load_textures(scene, dirname, params);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);
  trim_memory(scene);
}

// Save yaml
static void save_yaml(const string& filename, const scene_model& scene,
    bool ply_instances = false, const string& instances_name = "") {
  // open file
  auto fs = open_file(filename, "w");

  // write_yaml_comment(fs, get_save_scene_message(scene, ""));

  static const auto def_camera      = scene_camera{};
  static const auto def_texture     = scene_texture{};
  static const auto def_material    = scene_material{};
  static const auto def_shape       = scene_shape{};
  static const auto def_instance    = scene_instance{};
  static const auto def_environment = scene_environment{};

  auto yvalue = yaml_value{};

  if (!scene.cameras.empty()) write_yaml_object(fs, "cameras");
  for (auto& camera : scene.cameras) {
    write_yaml_property(
        fs, "cameras", "name", true, make_yaml_value(camera.name));
    if (camera.frame != identity3x4f)
      write_yaml_property(
          fs, "cameras", "frame", false, make_yaml_value(camera.frame));
    if (camera.orthographic)
      write_yaml_property(fs, "cameras", "orthographic", false,
          make_yaml_value(camera.orthographic));
    write_yaml_property(
        fs, "cameras", "lens", false, make_yaml_value(camera.lens));
    write_yaml_property(
        fs, "cameras", "aspect", false, make_yaml_value(camera.aspect));
    write_yaml_property(
        fs, "cameras", "film", false, make_yaml_value(camera.film));
    write_yaml_property(
        fs, "cameras", "focus", false, make_yaml_value(camera.focus));
    if (camera.aperture)
      write_yaml_property(
          fs, "cameras", "aperture", false, make_yaml_value(camera.aperture));
  }

  if (!scene.textures.empty()) write_yaml_object(fs, "textures");
  for (auto& texture : scene.textures) {
    write_yaml_property(
        fs, "textures", "name", true, make_yaml_value(texture.name));
    if (!texture.filename.empty())
      write_yaml_property(
          fs, "textures", "filename", false, make_yaml_value(texture.filename));
  }

  if (!scene.materials.empty()) write_yaml_object(fs, "materials");
  for (auto& material : scene.materials) {
    write_yaml_property(
        fs, "materials", "name", true, make_yaml_value(material.name));
    if (material.emission != zero3f)
      write_yaml_property(fs, "materials", "emission", false,
          make_yaml_value(material.emission));
    if (material.diffuse != zero3f)
      write_yaml_property(
          fs, "materials", "diffuse", false, make_yaml_value(material.diffuse));
    if (material.specular != zero3f)
      write_yaml_property(fs, "materials", "specular", false,
          make_yaml_value(material.specular));
    if (material.metallic)
      write_yaml_property(fs, "materials", "metallic", false,
          make_yaml_value(material.metallic));
    if (material.transmission != zero3f)
      write_yaml_property(fs, "materials", "transmission", false,
          make_yaml_value(material.transmission));
    write_yaml_property(fs, "materials", "roughness", false,
        make_yaml_value(material.roughness));
    if (material.refract)
      write_yaml_property(
          fs, "materials", "refract", false, make_yaml_value(material.refract));
    if (material.voltransmission != zero3f)
      write_yaml_property(fs, "materials", "voltransmission", false,
          make_yaml_value(material.voltransmission));
    if (material.volmeanfreepath != zero3f)
      write_yaml_property(fs, "materials", "volmeanfreepath", false,
          make_yaml_value(material.volmeanfreepath));
    if (material.volscatter != zero3f)
      write_yaml_property(fs, "materials", "volscatter", false,
          make_yaml_value(material.volscatter));
    if (material.volemission != zero3f)
      write_yaml_property(fs, "materials", "volemission", false,
          make_yaml_value(material.volemission));
    if (material.volanisotropy)
      write_yaml_property(fs, "materials", "volanisotropy", false,
          make_yaml_value(material.volanisotropy));
    if (material.voltransmission != zero3f ||
        material.volmeanfreepath != zero3f)
      write_yaml_property(fs, "materials", "volscale", false,
          make_yaml_value(material.volscale));
    if (material.coat != zero3f)
      write_yaml_property(
          fs, "materials", "coat", false, make_yaml_value(material.coat));
    if (material.opacity != 1)
      write_yaml_property(
          fs, "materials", "opacity", false, make_yaml_value(material.opacity));
    if (material.emission_tex >= 0)
      write_yaml_property(fs, "materials", "emission_tex", false,
          make_yaml_value(scene.textures[material.emission_tex].name));
    if (material.diffuse_tex >= 0)
      write_yaml_property(fs, "materials", "diffuse_tex", false,
          make_yaml_value(scene.textures[material.diffuse_tex].name));
    if (material.metallic_tex >= 0)
      write_yaml_property(fs, "materials", "metallic_tex", false,
          make_yaml_value(scene.textures[material.metallic_tex].name));
    if (material.specular_tex >= 0)
      write_yaml_property(fs, "materials", "specular_tex", false,
          make_yaml_value(scene.textures[material.specular_tex].name));
    if (material.roughness_tex >= 0)
      write_yaml_property(fs, "materials", "roughness_tex", false,
          make_yaml_value(scene.textures[material.roughness_tex].name));
    if (material.transmission_tex >= 0)
      write_yaml_property(fs, "materials", "transmission_tex", false,
          make_yaml_value(scene.textures[material.transmission_tex].name));
    if (material.subsurface_tex >= 0)
      write_yaml_property(fs, "materials", "subsurface_tex", false,
          make_yaml_value(scene.textures[material.subsurface_tex].name));
    if (material.coat_tex >= 0)
      write_yaml_property(fs, "materials", "coat_tex", false,
          make_yaml_value(scene.textures[material.coat_tex].name));
    if (material.opacity_tex >= 0)
      write_yaml_property(fs, "materials", "opacity_tex", false,
          make_yaml_value(scene.textures[material.opacity_tex].name));
    if (material.normal_tex >= 0)
      write_yaml_property(fs, "materials", "normal_tex", false,
          make_yaml_value(scene.textures[material.normal_tex].name));
    if (material.gltf_textures)
      write_yaml_property(fs, "materials", "gltf_textures", false,
          make_yaml_value(material.gltf_textures));
  }

  if (!scene.shapes.empty()) write_yaml_object(fs, "shapes");
  for (auto& shape : scene.shapes) {
    write_yaml_property(
        fs, "shapes", "name", true, make_yaml_value(shape.name));
    if (!shape.filename.empty())
      write_yaml_property(
          fs, "shapes", "filename", false, make_yaml_value(shape.filename));
    write_yaml_property(fs, "subdivs", "subdivisions", false,
        make_yaml_value(shape.subdivisions));
    write_yaml_property(fs, "subdivs", "catmullclark", false,
        make_yaml_value(shape.catmullclark));
    write_yaml_property(
        fs, "subdivs", "smooth", false, make_yaml_value(shape.smooth));
    if (shape.facevarying)
      write_yaml_property(fs, "subdivs", "facevarying", false,
          make_yaml_value(shape.facevarying));
    if (shape.displacement_tex >= 0)
      write_yaml_property(fs, "subdivs", "displacement_tex", false,
          make_yaml_value(scene.textures[shape.displacement_tex].name));
    if (shape.displacement_tex >= 0)
      write_yaml_property(fs, "subdivs", "displacement", false,
          make_yaml_value(shape.displacement));
  }

  if (!ply_instances) {
    if (!scene.instances.empty()) write_yaml_object(fs, "instances");
    for (auto& instance : scene.instances) {
      write_yaml_property(
          fs, "instances", "name", true, make_yaml_value(instance.name));
      if (instance.frame != identity3x4f)
        write_yaml_property(
            fs, "instances", "frame", false, make_yaml_value(instance.frame));
      if (instance.shape >= 0)
        write_yaml_property(fs, "instances", "shape", false,
            make_yaml_value(scene.shapes[instance.shape].name));
      if (instance.material >= 0)
        write_yaml_property(fs, "instances", "material", false,
            make_yaml_value(scene.materials[instance.material].name));
    }
  } else {
    if (!scene.instances.empty()) write_yaml_object(fs, "ply_instances");
    write_yaml_property(
        fs, "ply_instances", "filename", true, make_yaml_value(instances_name));
  }

  if (!scene.environments.empty()) write_yaml_object(fs, "environments");
  for (auto& environment : scene.environments) {
    write_yaml_property(
        fs, "environments", "name", true, make_yaml_value(environment.name));
    if (environment.frame != identity3x4f)
      write_yaml_property(fs, "environments", "frame", false,
          make_yaml_value(environment.frame));
    write_yaml_property(fs, "environments", "emission", false,
        make_yaml_value(environment.emission));
    if (environment.emission_tex >= 0)
      write_yaml_property(fs, "environments", "emission_tex", false,
          make_yaml_value(scene.textures[environment.emission_tex].name));
  }
}

// Save a scene in the builtin YAML format.
static void save_yaml_scene(const string& filename, const scene_model& scene,
    const save_params& params) {
  try {
    // save yaml file
    save_yaml(filename, scene);

    // save meshes and textures
    auto dirname = get_dirname(filename);
    save_shapes(scene, dirname, params);
    save_textures(scene, dirname, params);
  } catch (const std::exception& e) {
    throw std::runtime_error("cannot save scene " + filename + "\n" + e.what());
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

void load_obj(
    const string& filename, scene_model& scene, const load_params& params) {
  // load obj
  auto obj = obj_model{};
  load_obj(filename, obj, false, true, true);

  // convert cameras
  for (auto& ocam : obj.cameras) {
    auto& camera = scene.cameras.emplace_back();
    camera.name  = make_safe_name(ocam.name, "cam", (int)scene.cameras.size());
    camera.frame = ocam.frame;
    camera.orthographic = ocam.ortho;
    camera.film         = max(ocam.width, ocam.height);
    camera.aspect       = ocam.width / ocam.height;
    camera.focus        = ocam.focus;
    camera.lens         = ocam.lens;
    camera.aperture     = ocam.aperture;
  }

  // helper to create texture maps
  auto texture_map = unordered_map<string, int>{{"", -1}};
  auto get_texture = [&texture_map, &scene](const obj_texture_info& info) {
    if (info.path == "") return -1;
    auto it = texture_map.find(info.path);
    if (it != texture_map.end()) return it->second;
    auto& texture = scene.textures.emplace_back();
    texture.name  = make_safe_name(
        get_basename(info.path), "texture", (int)scene.textures.size());
    texture.filename       = info.path;
    texture_map[info.path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert materials and textures
  auto material_map = unordered_map<string, int>{{"", -1}};
  for (auto& omat : obj.materials) {
    auto& material = scene.materials.emplace_back();
    material.name  = make_safe_name(
        omat.name, "material", (int)scene.materials.size());
    material.emission         = omat.emission;
    material.diffuse          = omat.diffuse;
    material.specular         = omat.specular;
    material.roughness        = obj_exponent_to_roughness(omat.exponent);
    material.metallic         = omat.pbr_metallic;
    material.coat             = omat.reflection;
    material.transmission     = omat.transmission;
    material.voltransmission  = omat.vol_transmission;
    material.volmeanfreepath  = omat.vol_meanfreepath;
    material.volemission      = omat.vol_emission;
    material.volscatter       = omat.vol_scattering;
    material.volanisotropy    = omat.vol_anisotropy;
    material.volscale         = omat.vol_scale;
    material.opacity          = omat.opacity;
    material.emission_tex     = get_texture(omat.emission_map);
    material.diffuse_tex      = get_texture(omat.diffuse_map);
    material.specular_tex     = get_texture(omat.specular_map);
    material.metallic_tex     = get_texture(omat.pbr_metallic_map);
    material.roughness_tex    = get_texture(omat.pbr_roughness_map);
    material.transmission_tex = get_texture(omat.transmission_map);
    material.coat_tex         = get_texture(omat.reflection_map);
    material.opacity_tex      = get_texture(omat.opacity_map);
    material.normal_tex       = get_texture(omat.normal_map);
    material_map[omat.name]   = (int)scene.materials.size() - 1;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto& oshape : obj.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shape.name  = oshape.name;
    if (shape.name == "") shape.name = "shape";
    shape_name_counts[shape.name] += 1;
    if (shape_name_counts[shape.name] > 1)
      shape.name += std::to_string(shape_name_counts[shape.name]);
    shape.name = make_safe_name(shape.name, "shape", (int)scene.shapes.size());
    shape.filename  = make_safe_filename("shapes/" + shape.name + ".ply");
    auto materials  = vector<string>{};
    auto ematerials = vector<int>{};
    auto has_quads  = has_obj_quads(oshape);
    if (!oshape.faces.empty() && !params.facevarying && !has_quads) {
      get_obj_triangles(obj, oshape, shape.triangles, shape.positions,
          shape.normals, shape.texcoords, materials, ematerials, true);
    } else if (!oshape.faces.empty() && !params.facevarying && has_quads) {
      get_obj_quads(obj, oshape, shape.quads, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.lines.empty()) {
      get_obj_lines(obj, oshape, shape.lines, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.points.empty()) {
      get_obj_points(obj, oshape, shape.points, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.faces.empty() && params.facevarying) {
      get_obj_fvquads(obj, oshape, shape.quadspos, shape.quadsnorm,
          shape.quadstexcoord, shape.positions, shape.normals, shape.texcoords,
          materials, ematerials, true);
    } else {
      throw std::runtime_error("should not have gotten here");
    }
    // get material
    if (oshape.materials.size() != 1) {
      throw std::runtime_error("missing material for " + oshape.name);
    }
    auto material = -1;
    try {
      material = material_map.at(oshape.materials.at(0));
    } catch (...) {
      throw std::runtime_error(
          "cannot find material " + oshape.materials.at(0));
    }
    // make instances
    if (oshape.instances.empty()) {
      auto& instance    = scene.instances.emplace_back();
      instance.name     = shape.name;
      instance.material = material;
      instance.shape    = (int)scene.shapes.size() - 1;
    } else {
      for (auto& frame : oshape.instances) {
        auto& instance    = scene.instances.emplace_back();
        instance.name     = shape.name;
        instance.frame    = frame;
        instance.material = material;
        instance.shape    = (int)scene.shapes.size() - 1;
      }
    }
  }

  // convert environments
  for (auto& oenvironment : obj.environments) {
    auto& environment = scene.environments.emplace_back();
    environment.name  = make_safe_name(
        oenvironment.name, "environment", scene.environments.size());
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = get_texture(oenvironment.emission_map);
  }
}

// Loads an OBJ
static void load_obj_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  scene = {};

  // Parse obj
  load_obj(filename, scene, params);

  // load textures
  auto dirname = get_dirname(filename);
  load_textures(scene, dirname, params);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);
}

static void save_obj(const string& filename, const scene_model& scene,
    const save_params& params) {
  auto obj = obj_model{};

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = camera.name;
    ocamera.frame    = camera.frame;
    ocamera.ortho    = camera.orthographic;
    ocamera.width    = camera.film;
    ocamera.height   = camera.film / camera.aspect;
    ocamera.focus    = camera.focus;
    ocamera.lens     = camera.lens;
    ocamera.aperture = camera.aperture;
  }

  // textures
  auto get_texture = [&scene](int tex) {
    if (tex < 0) return obj_texture_info{};
    auto info = obj_texture_info{};
    info.path = scene.textures[tex].filename;
    return info;
  };

  // convert materials and textures
  for (auto& material : scene.materials) {
    auto& omaterial             = obj.materials.emplace_back();
    omaterial.name              = material.name;
    omaterial.illum             = 2;
    omaterial.emission          = material.emission;
    omaterial.diffuse           = material.diffuse;
    omaterial.specular          = material.specular;
    omaterial.exponent          = obj_roughness_to_exponent(material.roughness);
    omaterial.pbr_metallic      = material.metallic;
    omaterial.reflection        = material.coat;
    omaterial.transmission      = material.transmission;
    omaterial.opacity           = material.opacity;
    omaterial.emission_map      = get_texture(material.emission_tex);
    omaterial.diffuse_map       = get_texture(material.diffuse_tex);
    omaterial.specular_map      = get_texture(material.specular_tex);
    omaterial.pbr_metallic_map  = get_texture(material.metallic_tex);
    omaterial.pbr_roughness_map = get_texture(material.roughness_tex);
    omaterial.transmission_map  = get_texture(material.transmission_tex);
    omaterial.reflection_map    = get_texture(material.coat_tex);
    omaterial.opacity_map       = get_texture(material.opacity_tex);
    omaterial.normal_map        = get_texture(material.normal_tex);
    if (material.voltransmission != zero3f ||
        material.volmeanfreepath != zero3f) {
      omaterial.vol_transmission = material.voltransmission;
      omaterial.vol_meanfreepath = material.volmeanfreepath;
      omaterial.vol_emission     = material.volemission;
      omaterial.vol_scattering   = material.volscatter;
      omaterial.vol_anisotropy   = material.volanisotropy;
      omaterial.vol_scale        = material.volscale;
    }
  }

  // convert shapes
  if (params.objinstances) {
    for (auto& shape : scene.shapes) {
      if (!shape.triangles.empty()) {
        add_obj_triangles(obj, shape.name, shape.triangles, shape.positions,
            shape.normals, shape.texcoords, {}, {}, true);
      } else if (!shape.quads.empty()) {
        add_obj_quads(obj, shape.name, shape.quads, shape.positions,
            shape.normals, shape.texcoords, {}, {}, true);
      } else if (!shape.lines.empty()) {
        add_obj_lines(obj, shape.name, shape.lines, shape.positions,
            shape.normals, shape.texcoords, {}, {}, true);
      } else if (!shape.points.empty()) {
        add_obj_points(obj, shape.name, shape.points, shape.positions,
            shape.normals, shape.texcoords, {}, {}, true);
      } else if (!shape.quadspos.empty()) {
        add_obj_fvquads(obj, shape.name, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, shape.positions, shape.normals,
            shape.texcoords, {}, {}, true);
      } else {
        throw std::runtime_error("do not support empty shapes");
      }
    }
    for (auto& instance : scene.instances) {
      obj.shapes[instance.shape].instances.push_back(instance.frame);
    }
  } else {
    for (auto& instance : scene.instances) {
      auto& shape     = scene.shapes[instance.shape];
      auto  materials = vector{scene.materials[instance.material].name};
      auto  positions = shape.positions, normals = shape.normals;
      for (auto& p : positions) p = transform_point(instance.frame, p);
      for (auto& n : normals) n = transform_normal(instance.frame, n);
      if (!shape.triangles.empty()) {
        add_obj_triangles(obj, instance.name, shape.triangles, positions,
            normals, shape.texcoords, materials, {}, true);
      } else if (!shape.quads.empty()) {
        add_obj_quads(obj, instance.name, shape.quads, positions, normals,
            shape.texcoords, materials, {}, true);
      } else if (!shape.lines.empty()) {
        add_obj_lines(obj, instance.name, shape.lines, positions, normals,
            shape.texcoords, materials, {}, true);
      } else if (!shape.points.empty()) {
        add_obj_points(obj, instance.name, shape.points, positions, normals,
            shape.texcoords, materials, {}, true);
      } else if (!shape.quadspos.empty()) {
        add_obj_fvquads(obj, instance.name, shape.quadspos, shape.quadsnorm,
            shape.quadstexcoord, positions, normals, shape.texcoords, materials,
            {}, true);
      } else {
        throw std::runtime_error("do not support empty shapes");
      }
    }
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = environment.name;
    oenvironment.frame        = environment.frame;
    oenvironment.emission     = environment.emission;
    oenvironment.emission_map = get_texture(environment.emission_tex);
  }

  // save obj
  save_obj(filename, obj);
}

static void save_obj_scene(const string& filename, const scene_model& scene,
    const save_params& params) {
  save_obj(filename, scene, params);
  auto dirname = get_dirname(filename);
  save_textures(scene, dirname, params);
}

void print_obj_camera(const scene_camera& camera) {
  printf("c %s %d %g %g %g %g %g %g %g %g %g %g%g %g %g %g %g %g %g\n",
      camera.name.c_str(), (int)camera.orthographic, camera.film,
      camera.film / camera.aspect, camera.lens, camera.focus, camera.aperture,
      camera.frame.x.x, camera.frame.x.y, camera.frame.x.z, camera.frame.y.x,
      camera.frame.y.y, camera.frame.y.z, camera.frame.z.x, camera.frame.z.y,
      camera.frame.z.z, camera.frame.o.x, camera.frame.o.y, camera.frame.o.z);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_ply_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  scene = {};

  // load ply mesh
  scene.shapes.push_back({});
  auto& shape    = scene.shapes.back();
  shape.name     = "shape";
  shape.filename = get_filename(filename);
  load_shape(filename, shape.points, shape.lines, shape.triangles, shape.quads,
      shape.positions, shape.normals, shape.texcoords, shape.colors,
      shape.radius);

  // add instance
  auto instance  = scene_instance{};
  instance.name  = shape.name;
  instance.shape = 0;
  scene.instances.push_back(instance);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);
}

static void save_ply_scene(const string& filename, const scene_model& scene,
    const save_params& params) {
  if (scene.shapes.empty()) {
    throw std::runtime_error("cannot save empty scene " + filename);
  }
  auto& shape = scene.shapes.front();
  if (shape.quadspos.empty()) {
    save_shape(filename, shape.points, shape.lines, shape.triangles,
        shape.quads, shape.positions, shape.normals, shape.texcoords,
        shape.colors, shape.radius);
  } else {
    save_fvshape(filename, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
        shape.positions, shape.normals, shape.texcoords);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// convert gltf to scene
static void load_gltf(const string& filename, scene_model& scene) {
  auto gltf = gltf_model{};
  load_gltf(filename, gltf);

  // convert textures
  for (auto& gtexture : gltf.textures) {
    auto& texture = scene.textures.emplace_back();
    if (!gtexture.name.empty()) {
      texture.name = make_safe_name(
          gtexture.name, "texture", (int)scene.textures.size());
    } else {
      texture.name = make_safe_name(get_basename(gtexture.filename), "texture",
          (int)scene.textures.size());
    }
    texture.filename = gtexture.filename;
  }

  // convert materials
  for (auto& gmaterial : gltf.materials) {
    auto& material = scene.materials.emplace_back();
    material.name  = make_safe_name(
        gmaterial.name, "material", (int)scene.materials.size());
    material.emission     = gmaterial.emission;
    material.emission_tex = gmaterial.emission_tex;
    if (gmaterial.has_specgloss) {
      material.diffuse      = xyz(gmaterial.sg_diffuse);
      material.opacity      = gmaterial.sg_diffuse.w;
      material.specular     = gmaterial.sg_specular;
      material.diffuse_tex  = gmaterial.sg_diffuse_tex;
      material.specular_tex = gmaterial.sg_specular_tex;
    } else if (gmaterial.has_metalrough) {
      material.diffuse      = xyz(gmaterial.mr_base);
      material.opacity      = gmaterial.mr_base.w;
      material.specular     = vec3f{0.04f};
      material.diffuse_tex  = gmaterial.mr_base_tex;
      material.metallic_tex = gmaterial.mr_metallic_tex;
    }
    material.normal_tex = gmaterial.normal_tex;
  }

  // convert shapes
  auto shape_indices = vector<vector<vec2i>>{};
  for (auto& gmesh : gltf.meshes) {
    shape_indices.push_back({});
    for (auto& gprim : gmesh.primitives) {
      auto& shape = scene.shapes.emplace_back();
      shape_indices.back().push_back(
          {(int)scene.shapes.size() - 1, gprim.material});
      shape.name =
          gmesh.name.empty()
              ? ""s
              : (gmesh.name + std::to_string(shape_indices.back().size()));
      make_safe_name(shape.name, "shape", (int)scene.shapes.size());
      shape.filename = make_safe_filename(
          "shapes/shape" + std::to_string(scene.shapes.size()));
      shape.positions = gprim.positions;
      shape.normals   = gprim.normals;
      shape.texcoords = gprim.texcoords;
      shape.colors    = gprim.colors;
      shape.radius    = gprim.radius;
      shape.tangents  = gprim.tangents;
      shape.triangles = gprim.triangles;
      shape.lines     = gprim.lines;
      shape.points    = gprim.points;
    }
  }

  // convert cameras
  auto cameras = vector<scene_camera>{};
  for (auto& gcamera : gltf.cameras) {
    auto& camera  = cameras.emplace_back();
    camera.name   = gcamera.name;
    camera.aspect = gcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = gcamera.aspect >= 1
                      ? (2 * camera.aspect * tan(gcamera.yfov / 2))
                      : (2 * tan(gcamera.yfov / 2));
    camera.focus = 10;
  }

  // convert scene nodes
  for (auto& gnode : gltf.nodes) {
    if (gnode.camera >= 0) {
      auto& camera = scene.cameras.emplace_back(cameras[gnode.camera]);
      camera.name  = make_safe_name(
          camera.name, "caemra", (int)scene.cameras.size());
      camera.frame = gnode.frame;
    }
    if (gnode.mesh >= 0) {
      for (auto [shape, material] : shape_indices[gnode.mesh]) {
        auto& instance = scene.instances.emplace_back();
        instance.name  = make_safe_name(
            scene.shapes[shape].name, "instance", (int)scene.instances.size());
        instance.frame    = gnode.frame;
        instance.shape    = shape;
        instance.material = material;
      }
    }
  }
}

// Load a scene
static void load_gltf_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  // initialization
  scene = {};

  // load gltf
  load_gltf(filename, scene);

  // load textures
  auto dirname = get_dirname(filename);
  load_textures(scene, dirname, params);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);

  // fix cameras
  auto bbox = compute_bounds(scene);
  for (auto& camera : scene.cameras) {
    auto center   = (bbox.min + bbox.max) / 2;
    auto distance = dot(-camera.frame.z, center - camera.frame.o);
    if (distance > 0) camera.focus = distance;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

static void load_pbrt(
    const string& filename, scene_model& scene, const load_params& params) {
  // load pbrt
  auto pbrt = pbrt_model{};
  load_pbrt(filename, pbrt);

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.name   = make_safe_name("", "camera", (int)scene.cameras.size());
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = pcamera.lens;
    camera.focus  = pcamera.focus;
  }

  // convert textures
  auto texture_map = unordered_map<string, int>{{"", -1}};
  for (auto& ptexture : pbrt.textures) {
    if (ptexture.filename.empty()) continue;
    auto& texture = scene.textures.emplace_back();
    texture.name  = make_safe_name(
        ptexture.name, "texture", (int)scene.textures.size());
    texture.filename           = ptexture.filename;
    texture_map[ptexture.name] = (int)scene.textures.size() - 1;
  }

  // convert materials
  auto get_texture = [&texture_map](const string& name) {
    if (name == "") return -1;
    if (texture_map.find(name) == texture_map.end())
      throw std::runtime_error("cannot find texture " + name);
    return texture_map.at(name);
  };
  auto material_map = unordered_map<string, int>{{"", -1}};
  for (auto& pmaterial : pbrt.materials) {
    auto& material = scene.materials.emplace_back();
    material.name  = make_safe_name(
        pmaterial.name, "material", (int)scene.materials.size());
    material.diffuse      = pmaterial.diffuse;
    material.specular     = pmaterial.sspecular;
    material.transmission = pmaterial.transmission;
    material.roughness    = mean(pmaterial.roughness);
    material.opacity      = pmaterial.opacity == vec3f{1} ? 1
                                                     : mean(pmaterial.opacity);
    material.diffuse_tex         = get_texture(pmaterial.diffuse_map);
    material_map[pmaterial.name] = (int)scene.materials.size() - 1;
  }

  // convert arealights
  auto arealight_map = unordered_map<string, int>{{"", -1}};
  for (auto& parealight : pbrt.arealights) {
    auto& material = scene.materials.emplace_back();
    material.name  = make_safe_name(
        parealight.name, "arealight", (int)arealight_map.size());
    material.emission              = parealight.emission;
    arealight_map[parealight.name] = (int)scene.materials.size() - 1;
  }

  // convert shapes
  for (auto& pshape : pbrt.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shape.name  = make_safe_name(
        get_basename(shape.filename), "shape", (int)scene.shapes.size());
    if (pshape.filename.empty()) {
      shape.name     = make_safe_name("", "shape", (int)scene.shapes.size());
      shape.filename = make_safe_filename(
          "shapes/shape" + std::to_string(scene.shapes.size()) + ".ply");
    } else {
      shape.filename = pshape.filename;
      shape.name     = make_safe_name(
          get_basename(pshape.filename), "shape", (int)scene.shapes.size());
    }
    shape.positions = pshape.positions;
    shape.normals   = pshape.normals;
    shape.texcoords = pshape.texcoords;
    shape.triangles = pshape.triangles;
    for (auto& uv : shape.texcoords) uv.y = 1 - uv.y;
    auto material_id  = material_map.at(pshape.material);
    auto arealight_id = arealight_map.at(pshape.arealight);
    auto instance_id  = 0;
    for (auto& frame : pshape.instance_frames) {
      auto& instance    = scene.instances.emplace_back();
      instance.name     = shape.name + (pshape.instance_frames.empty()
                                           ? ""s
                                           : std::to_string(instance_id++));
      instance.frame    = frame * pshape.frame;
      instance.material = arealight_id >= 0 ? arealight_id : material_id;
      instance.shape    = (int)scene.shapes.size() - 1;
    }
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment = scene.environments.emplace_back();
    environment.name  = make_safe_name(
        "", "environment", (int)scene.environments.size());
    environment.frame    = penvironment.frame;
    environment.emission = penvironment.emission;
    if (!penvironment.filename.empty()) {
      auto& texture    = scene.textures.emplace_back();
      texture.name     = make_safe_name(get_basename(penvironment.filename),
          "environment", (int)scene.environments.size());
      texture.filename = penvironment.filename;
      environment.emission_tex = (int)scene.textures.size() - 1;
    } else {
      environment.emission_tex = -1;
    }
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape       = scene.shapes.emplace_back();
    shape.name        = make_safe_name("", "light", (int)scene.shapes.size());
    shape.filename    = make_safe_filename("shapes/" + shape.name + ".ply");
    shape.triangles   = plight.area_triangles;
    shape.positions   = plight.area_positions;
    shape.normals     = plight.area_normals;
    auto& material    = scene.materials.emplace_back();
    material.name     = shape.name;
    material.emission = plight.area_emission;
    auto& instance    = scene.instances.emplace_back();
    instance.name     = shape.name;
    instance.frame    = plight.area_frame;
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = (int)scene.materials.size() - 1;
  }
}

// load pbrt scenes
static void load_pbrt_scene(
    const string& filename, scene_model& scene, const load_params& params) {
  scene = scene_model{};

  // Parse pbrt
  load_pbrt(filename, scene, params);

  // load textures
  auto dirname = get_dirname(filename);
  load_shapes(scene, dirname, params);
  load_textures(scene, dirname, params);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);
}

// Convert a scene to pbrt format
static void save_pbrt(const string& filename, const scene_model& scene) {
  auto pbrt = pbrt_model{};

  // embed data
  for (auto stat : format_stats(scene)) pbrt.comments.push_back(stat);

  // convert camera
  auto& camera     = scene.cameras.front();
  auto& pcamera    = pbrt.cameras.emplace_back();
  pcamera.frame    = camera.frame;
  pcamera.lens     = camera.lens;
  pcamera.aspect   = camera.aspect;
  auto& pfilm      = pbrt.films.emplace_back();
  pfilm.filename   = "out.png";
  pfilm.resolution = {1280, (int)(1280 / pcamera.aspect)};

  // convert textures
  for (auto& texture : scene.textures) {
    auto& ptexture    = pbrt.textures.emplace_back();
    ptexture.name     = texture.name;
    ptexture.filename = texture.filename;
  }

  // convert materials
  for (auto& material : scene.materials) {
    auto& pmaterial        = pbrt.materials.emplace_back();
    pmaterial.name         = material.name;
    pmaterial.diffuse      = material.diffuse;
    pmaterial.specular     = material.specular;
    pmaterial.transmission = material.transmission;
    pmaterial.roughness    = {material.roughness, material.roughness};
    pmaterial.diffuse_map  = material.diffuse_tex >= 0
                                ? scene.textures[material.diffuse_tex].name
                                : ""s;
    auto& parealight    = pbrt.arealights.emplace_back();
    parealight.name     = material.name;
    parealight.emission = material.emission;
  }

  // convert instances
  for (auto& instance : scene.instances) {
    auto& shape      = scene.shapes[instance.shape];
    auto& material   = scene.materials[instance.material];
    auto& pshape     = pbrt.shapes.emplace_back();
    pshape.filename  = replace_extension(shape.filename, ".ply");
    pshape.frame     = instance.frame;
    pshape.material  = material.name;
    pshape.arealight = material.emission == zero3f ? ""s : material.name;
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment    = pbrt.environments.emplace_back();
    penvironment.emission = environment.emission;
    if (environment.emission_tex >= 0) {
      penvironment.filename = scene.textures[environment.emission_tex].filename;
    }
  }

  save_pbrt(filename, pbrt);
}

// Save a pbrt scene
void save_pbrt_scene(const string& filename, const scene_model& scene,
    const save_params& params) {
  // save pbrt
  save_pbrt(filename, scene);

  // save meshes
  auto dirname = get_dirname(filename);
  for (auto& shape : scene.shapes) {
    if (shape.quadspos.empty()) {
      save_shape(replace_extension(dirname + shape.filename, ".ply"),
          shape.points, shape.lines, shape.triangles, shape.quads,
          shape.positions, shape.normals, shape.texcoords, shape.colors,
          shape.radius);
    } else {
      save_fvshape(replace_extension(dirname + shape.filename, ".ply"),
          shape.quadspos, shape.quadsnorm, shape.quadstexcoord, shape.positions,
          shape.normals, shape.texcoords);
    }
  }

  // save textures
  save_textures(scene, dirname, params);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox_scene(scene_model& scene) {
  scene.name              = "cornellbox";
  auto& camera            = scene.cameras.emplace_back();
  camera.name             = "camera";
  camera.frame            = frame3f{{0, 1, 3.9}};
  camera.lens             = 0.035;
  camera.aperture         = 0.0;
  camera.film             = 0.024;
  camera.aspect           = 1;
  auto& floor_mat         = scene.materials.emplace_back();
  floor_mat.name          = "floor";
  floor_mat.diffuse       = {0.725, 0.71, 0.68};
  auto& ceiling_mat       = scene.materials.emplace_back();
  ceiling_mat.name        = "ceiling";
  ceiling_mat.diffuse     = {0.725, 0.71, 0.68};
  auto& backwall_mat      = scene.materials.emplace_back();
  backwall_mat.name       = "backwall";
  backwall_mat.diffuse    = {0.725, 0.71, 0.68};
  auto& rightwall_mat     = scene.materials.emplace_back();
  rightwall_mat.name      = "rightwall";
  rightwall_mat.diffuse   = {0.14, 0.45, 0.091};
  auto& leftwall_mat      = scene.materials.emplace_back();
  leftwall_mat.name       = "leftwall";
  leftwall_mat.diffuse    = {0.63, 0.065, 0.05};
  auto& shortbox_mat      = scene.materials.emplace_back();
  shortbox_mat.name       = "shortbox";
  shortbox_mat.diffuse    = {0.725, 0.71, 0.68};
  auto& tallbox_mat       = scene.materials.emplace_back();
  tallbox_mat.name        = "tallbox";
  tallbox_mat.diffuse     = {0.725, 0.71, 0.68};
  auto& light_mat         = scene.materials.emplace_back();
  light_mat.name          = "light";
  light_mat.emission      = {17, 12, 4};
  auto& floor_shp         = scene.shapes.emplace_back();
  floor_shp.name          = "floor";
  floor_shp.filename      = "shapes/floor.obj";
  floor_shp.positions     = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor_shp.triangles     = {{0, 1, 2}, {2, 3, 0}};
  auto& ceiling_shp       = scene.shapes.emplace_back();
  ceiling_shp.name        = "ceiling";
  ceiling_shp.name        = "shapes/ceiling.obj";
  ceiling_shp.positions   = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling_shp.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& backwall_shp      = scene.shapes.emplace_back();
  backwall_shp.name       = "backwall";
  backwall_shp.filename   = "shapes/backwall.obj";
  backwall_shp.positions  = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall_shp.triangles  = {{0, 1, 2}, {2, 3, 0}};
  auto& rightwall_shp     = scene.shapes.emplace_back();
  rightwall_shp.name      = "rightwall";
  rightwall_shp.filename  = "shapes/rightwall.obj";
  rightwall_shp.positions = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& leftwall_shp      = scene.shapes.emplace_back();
  leftwall_shp.name       = "leftwall";
  leftwall_shp.filename   = "shapes/leftwall.obj";
  leftwall_shp.positions  = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall_shp.triangles  = {{0, 1, 2}, {2, 3, 0}};
  auto& shortbox_shp      = scene.shapes.emplace_back();
  shortbox_shp.name       = "shortbox";
  shortbox_shp.filename   = "shapes/shortbox.obj";
  shortbox_shp.positions  = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  shortbox_shp.triangles  = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& tallbox_shp       = scene.shapes.emplace_back();
  tallbox_shp.name        = "tallbox";
  tallbox_shp.filename    = "shapes/tallbox.obj";
  tallbox_shp.positions   = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
      {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
      {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
      {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
      {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
      {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
      {-0.71, 0.0, -0.49}};
  tallbox_shp.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& light_shp         = scene.shapes.emplace_back();
  light_shp.name          = "light";
  light_shp.filename      = "shapes/light.obj";
  light_shp.positions     = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
      {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  light_shp.triangles     = {{0, 1, 2}, {2, 3, 0}};
  scene.instances.push_back({"floor", identity3x4f, 0, 0});
  scene.instances.push_back({"ceiling", identity3x4f, 1, 1});
  scene.instances.push_back({"backwall", identity3x4f, 2, 2});
  scene.instances.push_back({"rightwall", identity3x4f, 3, 3});
  scene.instances.push_back({"leftwall", identity3x4f, 4, 4});
  scene.instances.push_back({"shortbox", identity3x4f, 5, 5});
  scene.instances.push_back({"tallbox", identity3x4f, 6, 6});
  scene.instances.push_back({"light", identity3x4f, 7, 7});
}

}  // namespace yocto
