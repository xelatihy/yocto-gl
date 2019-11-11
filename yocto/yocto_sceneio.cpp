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

//
// TODO: update transforms -> should be compute transforms?
// TODO: update tesselation -> should be compute tesselation
// TODO: move out animation utilities
//

#include "yocto_sceneio.h"

#include <atomic>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <deque>
#include <future>

#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
static inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)

    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get directory name (including '/').
static inline string get_dirname(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
static inline string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
static inline string get_filename(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
static inline string get_noextension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
static inline string get_basename(const string& filename) {
  return get_noextension(get_filename(filename));
}

// Replaces extensions
static inline string replace_extension(
    const string& filename, const string& ext) {
  return get_noextension(filename) + ext;
}

}  // namespace yocto

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
bbox3f compute_bounds(const sceneio_model& scene) {
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
void add_cameras(sceneio_model& scene) {
  if (!scene.cameras.empty()) return;
  auto& camera = scene.cameras.emplace_back();
  camera.name  = "default";
  // TODO: error in camera.lens and camera.film
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
void add_materials(sceneio_model& scene) {
  auto material_id = -1;
  for (auto& instance : scene.instances) {
    if (instance.material >= 0) continue;
    if (material_id < 0) {
      auto material    = sceneio_material{};
      material.name    = "default";
      material.diffuse = {0.2f, 0.2f, 0.2f};
      scene.materials.push_back(material);
      material_id = (int)scene.materials.size() - 1;
    }
    instance.material = material_id;
  }
}

// Add missing radius.
void add_radius(sceneio_model& scene, float radius = 0.001f) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add a sky environment
void add_sky(sceneio_model& scene, float sun_angle) {
  auto texture     = sceneio_texture{};
  texture.name     = "sky";
  texture.filename = "textures/sky.hdr";
  texture.hdr      = make_sunsky({1024, 512}, sun_angle);
  scene.textures.push_back(texture);
  auto environment         = sceneio_environment{};
  environment.name         = "sky";
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
  scene.environments.push_back(environment);
}

vector<string> scene_stats(const sceneio_model& scene, bool verbose) {
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
void trim_memory(sceneio_model& scene) {
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
vector<string> scene_validation(const sceneio_model& scene, bool notextures) {
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
  auto check_empty_textures = [&errs](const vector<sceneio_texture>& vals) {
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
sceneio_shape subdivide_shape(const sceneio_shape& shape) {
  using std::ignore;
  if (!shape.subdivisions) return shape;
  auto subdiv         = shape;
  subdiv.subdivisions = 0;
  if (!shape.points.empty()) {
    throw std::runtime_error("point subdivision not supported");
  } else if (!shape.lines.empty()) {
    tie(ignore, subdiv.normals) = subdivide_lines(
        subdiv.lines, subdiv.normals, shape.subdivisions);
    tie(ignore, subdiv.texcoords) = subdivide_lines(
        subdiv.lines, subdiv.texcoords, shape.subdivisions);
    tie(ignore, subdiv.colors) = subdivide_lines(
        subdiv.lines, subdiv.colors, shape.subdivisions);
    tie(ignore, subdiv.radius) = subdivide_lines(
        subdiv.lines, subdiv.radius, shape.subdivisions);
    tie(subdiv.lines, subdiv.positions) = subdivide_lines(
        subdiv.lines, subdiv.positions, shape.subdivisions);
    if (shape.smooth)
      subdiv.normals = compute_tangents(subdiv.lines, subdiv.positions);
  } else if (!shape.triangles.empty()) {
    tie(ignore, subdiv.normals) = subdivide_triangles(
        subdiv.triangles, subdiv.normals, shape.subdivisions);
    tie(ignore, subdiv.texcoords) = subdivide_triangles(
        subdiv.triangles, subdiv.texcoords, shape.subdivisions);
    tie(ignore, subdiv.colors) = subdivide_triangles(
        subdiv.triangles, subdiv.colors, shape.subdivisions);
    tie(ignore, subdiv.radius) = subdivide_triangles(
        subdiv.triangles, subdiv.radius, shape.subdivisions);
    tie(subdiv.triangles, subdiv.positions) = subdivide_triangles(
        subdiv.triangles, subdiv.positions, shape.subdivisions);
    if (shape.smooth)
      subdiv.normals = compute_normals(subdiv.triangles, subdiv.positions);
  } else if (!shape.quads.empty() && !shape.catmullclark) {
    tie(ignore, subdiv.normals) = subdivide_quads(
        subdiv.quads, subdiv.normals, shape.subdivisions);
    tie(ignore, subdiv.texcoords) = subdivide_quads(
        subdiv.quads, subdiv.texcoords, shape.subdivisions);
    tie(ignore, subdiv.colors) = subdivide_quads(
        subdiv.quads, subdiv.colors, shape.subdivisions);
    tie(ignore, subdiv.radius) = subdivide_quads(
        subdiv.quads, subdiv.radius, shape.subdivisions);
    tie(subdiv.quads, subdiv.positions) = subdivide_quads(
        subdiv.quads, subdiv.positions, shape.subdivisions);
    if (subdiv.smooth)
      subdiv.normals = compute_normals(subdiv.quads, subdiv.positions);
  } else if (!shape.quads.empty() && shape.catmullclark) {
    tie(ignore, subdiv.normals) = subdivide_catmullclark(
        subdiv.quads, subdiv.normals, shape.subdivisions);
    tie(ignore, subdiv.texcoords) = subdivide_catmullclark(
        subdiv.quads, subdiv.texcoords, shape.subdivisions);
    tie(ignore, subdiv.colors) = subdivide_catmullclark(
        subdiv.quads, subdiv.colors, shape.subdivisions);
    tie(ignore, subdiv.radius) = subdivide_catmullclark(
        subdiv.quads, subdiv.radius, shape.subdivisions);
    tie(subdiv.quads, subdiv.positions) = subdivide_catmullclark(
        subdiv.quads, subdiv.positions, shape.subdivisions);
    if (subdiv.smooth)
      subdiv.normals = compute_normals(subdiv.quads, subdiv.positions);
  } else if (!shape.quadspos.empty() && !shape.catmullclark) {
    std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_quads(
        subdiv.quadsnorm, subdiv.normals, shape.subdivisions);
    std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_quads(
        subdiv.quadstexcoord, subdiv.texcoords, shape.subdivisions);
    std::tie(subdiv.quadspos, subdiv.positions) = subdivide_quads(
        subdiv.quadspos, subdiv.positions, shape.subdivisions);
    if (subdiv.smooth) {
      subdiv.normals   = compute_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    }
  } else if (!shape.quadspos.empty() && shape.catmullclark) {
    std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_catmullclark(
        subdiv.quadstexcoord, subdiv.texcoords, shape.subdivisions, true);
    std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_catmullclark(
        subdiv.quadsnorm, subdiv.normals, shape.subdivisions, true);
    std::tie(subdiv.quadspos, subdiv.positions) = subdivide_catmullclark(
        subdiv.quadspos, subdiv.positions, shape.subdivisions);
    if (shape.smooth) {
      subdiv.normals   = compute_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    } else {
      subdiv.normals   = {};
      subdiv.quadsnorm = {};
    }
  } else {
    throw std::runtime_error("empty shape");
  }
  return subdiv;
}
// Apply displacement to a shape
sceneio_shape displace_shape(
    const sceneio_model& scene, const sceneio_shape& shape) {
  // Evaluate a texture
  auto eval_texture = [](const sceneio_texture& texture,
                          const vec2f&          texcoord) -> vec4f {
    if (!texture.hdr.empty()) {
      return eval_image(texture.hdr, texcoord, false, false);
    } else if (!texture.ldr.empty()) {
      return eval_image(texture.ldr, texcoord, true, false, false);
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
      auto disp = mean(xyz(eval_texture(displacement, shape.texcoords[vid])));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      subdiv.positions[vid] += normals[vid] * shape.displacement * disp;
    }
    if (shape.smooth || !shape.normals.empty()) {
      subdiv.normals = compute_normals(subdiv.triangles, subdiv.positions);
    }
  } else if (!shape.quads.empty()) {
    auto normals = shape.normals;
    if (normals.empty())
      normals = compute_normals(shape.triangles, shape.positions);
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
      auto disp = mean(xyz(eval_texture(displacement, shape.texcoords[vid])));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      subdiv.positions[vid] += normals[vid] * shape.displacement * disp;
    }
    if (shape.smooth || !shape.normals.empty()) {
      subdiv.normals = compute_normals(subdiv.quads, subdiv.positions);
    }
  } else if (!shape.quadspos.empty()) {
    // facevarying case
    auto offset = vector<float>(shape.positions.size(), 0);
    auto count  = vector<int>(shape.positions.size(), 0);
    for (auto fid = 0; fid < shape.quadspos.size(); fid++) {
      auto qpos = shape.quadspos[fid];
      auto qtxt = shape.quadstexcoord[fid];
      for (auto i = 0; i < 4; i++) {
        auto disp = mean(
            xyz(eval_texture(displacement, shape.texcoords[qtxt[i]])));
        if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
        offset[qpos[i]] += shape.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = compute_normals(shape.quadspos, shape.positions);
    for (auto vid = 0; vid < shape.positions.size(); vid++) {
      subdiv.positions[vid] += normals[vid] * offset[vid] / count[vid];
    }
    if (shape.smooth || !shape.normals.empty()) {
      subdiv.quadsnorm = shape.quadspos;
      subdiv.normals   = compute_normals(subdiv.quadspos, subdiv.positions);
    }
  }
  return subdiv;
}
sceneio_shape tesselate_shape(const sceneio_model& scene,
    const sceneio_shape& shape, bool no_quads, bool no_facevarying) {
  if (!needs_tesselation(scene, shape, no_quads, no_facevarying)) return shape;
  auto subdiv = shape;
  if (subdiv.subdivisions) subdiv = subdivide_shape(subdiv);
  if (subdiv.displacement) subdiv = displace_shape(scene, subdiv);
  if (!subdiv.quadspos.empty() && no_facevarying) {
    std::tie(subdiv.quads, subdiv.positions, subdiv.normals, subdiv.texcoords) =
        split_facevarying(subdiv.quadspos, subdiv.quadsnorm,
            subdiv.quadstexcoord, subdiv.positions, subdiv.normals,
            subdiv.texcoords);
  }
  if (!subdiv.quads.empty() && no_quads) {
    subdiv.triangles = quads_to_triangles(subdiv.quads);
    subdiv.quads     = {};
  }
  return subdiv;
}
bool needs_tesselation(const sceneio_model& scene, const sceneio_shape& shape,
    bool no_quads, bool no_facevarying) {
  if (!shape.quadspos.empty() && (no_facevarying || no_quads)) return true;
  if (!shape.quads.empty() && no_quads) return true;
  if (shape.subdivisions || shape.displacement) return true;
  return false;
}

// Update animation transforms
void update_transforms(sceneio_model& scene, sceneio_animation& animation,
    float time, const string& anim_group) {
  if (anim_group != "" && anim_group != animation.group) return;

  if (!animation.translations.empty()) {
    auto value = vec3f{0, 0, 0};
    switch (animation.interpolation) {
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.translations, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.translations, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
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
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.rotations, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.rotations, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.rotations, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].rotation = value;
  }
  if (!animation.scales.empty()) {
    auto value = vec3f{1, 1, 1};
    switch (animation.interpolation) {
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.scales, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.scales, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.scales, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].scale = value;
  }
}

// Update node transforms
void update_transforms(sceneio_model& scene, sceneio_node& node,
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
    sceneio_model& scene, float time, const string& anim_group) {
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
static sceneio_status load_yaml_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static sceneio_status save_yaml_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load/save a scene from/to OBJ.
static sceneio_status load_obj_scene(const string& filename,
    sceneio_model& scene, bool facevarying, bool noparallel);
static sceneio_status save_obj_scene(const string& filename,
    const sceneio_model& scene, bool instances, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static sceneio_status load_ply_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static sceneio_status save_ply_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load/save a scene from/to glTF.
static sceneio_status load_gltf_scene(
    const string& filename, sceneio_model& scene, bool noparallel);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static sceneio_status load_pbrt_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static sceneio_status save_pbrt_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load a scene
sceneio_status load_scene(const string& filename, sceneio_model& scene,
    bool obj_facevarying, bool noparallel) {
  auto ext = get_extension(filename);
  if (ext == ".yaml" || ext == ".YAML") {
    return load_yaml_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, obj_facevarying, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, noparallel);
  } else {
    scene = {};
    return {filename + ": unsupported format"};
  }
}

// Save a scene
sceneio_status save_scene(const string& filename, const sceneio_model& scene,
    bool obj_instances, bool noparallel) {
  auto ext = get_extension(filename);
  if (ext == ".yaml" || ext == ".YAML") {
    return save_yaml_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, obj_instances, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, noparallel);
  } else {
    return {filename + ": unsupported format"};
  }
}

sceneio_status load_texture(const string& filename, sceneio_texture& texture) {
  if (is_hdr_filename(texture.filename)) {
    if (auto ret = load_image(
            get_dirname(filename) + texture.filename, texture.hdr);
        !ret) {
      return {filename + ": missing texture (" + ret.error + ")"};
    } else {
      return {};
    }
  } else {
    if (auto ret = load_imageb(
            get_dirname(filename) + texture.filename, texture.ldr);
        !ret) {
      return {filename + ": missing texture (" + ret.error + ")"};
    } else {
      return {};
    }
  }
}

sceneio_status save_texture(
    const string& filename, const sceneio_texture& texture) {
  if (!texture.hdr.empty()) {
    if (auto ret = save_image(
            get_dirname(filename) + texture.filename, texture.hdr);
        !ret) {
      return {filename + ": missing texture (" + ret.error + ")"};
    } else {
      return {};
    }
  } else {
    if (auto ret = save_imageb(
            get_dirname(filename) + texture.filename, texture.ldr);
        !ret) {
      return {filename + ": missing texture (" + ret.error + ")"};
    } else {
      return {};
    }
  }
}

sceneio_status load_textures(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // load images
  if (noparallel) {
    for (auto& texture : scene.textures) {
      if (!texture.hdr.empty() || !texture.ldr.empty()) continue;
      if (auto ret = load_texture(filename, texture); !ret) return ret;
    }
    return {};
  } else {
    auto mutex  = std::mutex{};
    auto status = sceneio_status{};
    parallel_foreach(
        scene.textures, [&filename, &mutex, &status](sceneio_texture& texture) {
          if (!status) return;
          if (!texture.hdr.empty() || !texture.ldr.empty()) return;
          if (auto ret = load_texture(filename, texture); !ret) {
            auto lock = std::lock_guard{mutex};
            status    = ret;
          }
        });
    return status;
  }
}

// helper to save textures
sceneio_status save_textures(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // save images
  if (noparallel) {
    for (auto& texture : scene.textures) {
      if (auto ret = save_texture(filename, texture); !ret) return ret;
    }
    return {};
  } else {
    auto mutex  = std::mutex{};
    auto status = sceneio_status{};
    parallel_foreach(scene.textures,
        [&filename, &mutex, &status](const sceneio_texture& texture) {
          if (!status) return;
          if (auto ret = save_texture(filename, texture); !ret) {
            auto lock = std::lock_guard{mutex};
            status    = ret;
          }
        });
    return status;
  }
}

sceneio_status load_shape(const string& filename, sceneio_shape& shape) {
  if (!shape.facevarying) {
    if (auto ret = load_shape(get_dirname(filename) + shape.filename,
            shape.points, shape.lines, shape.triangles, shape.quads,
            shape.positions, shape.normals, shape.texcoords, shape.colors,
            shape.radius);
        !ret) {
      return {filename + ": missing shape (" + ret.error + ")"};
    } else {
      return {};
    }
  } else {
    if (auto ret = load_fvshape(get_dirname(filename) + shape.filename,
            shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords);
        !ret) {
      return {filename + ": missing shape (" + ret.error + ")"};
    } else {
      return {};
    }
  }
}

sceneio_status save_shape(const string& filename, const sceneio_shape& shape) {
  if (shape.quadspos.empty()) {
    if (auto ret = save_shape(get_dirname(filename) + shape.filename,
            shape.points, shape.lines, shape.triangles, shape.quads,
            shape.positions, shape.normals, shape.texcoords, shape.colors,
            shape.radius);
        !ret) {
      return {filename + ": missing shape (" + ret.error + ")"};
    } else {
      return {};
    }
  } else {
    if (auto ret = save_fvshape(get_dirname(filename) + shape.filename,
            shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords);
        !ret) {
      return {filename + ": missing shape (" + ret.error + ")"};
    } else {
      return {};
    }
  }
}

// Load json meshes
sceneio_status load_shapes(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // load shapes
  if (noparallel) {
    for (auto& shape : scene.shapes) {
      if (!shape.positions.empty()) continue;
      if (auto ret = load_shape(filename, shape); !ret) return ret;
    }
    return {};
  } else {
    auto mutex  = std::mutex{};
    auto status = sceneio_status{};
    parallel_foreach(
        scene.shapes, [&filename, &status, &mutex](sceneio_shape& shape) {
          if (!status) return;
          if (!shape.positions.empty()) return;
          if (auto ret = load_shape(filename, shape); !ret) {
            auto lock = std::lock_guard{mutex};
            status    = ret;
          }
        });
    return status;
  }
}

// Save json meshes
sceneio_status save_shapes(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // save shapes
  if (noparallel) {
    for (auto& shape : scene.shapes) {
      if (auto ret = save_shape(filename, shape); !ret)
        return {filename + ": missing shape (" + ret.error + ")"};
    }
    return {};
  } else {
    auto mutex  = std::mutex{};
    auto status = sceneio_status{};
    parallel_foreach(
        scene.shapes, [&filename, &status, &mutex](const sceneio_shape& shape) {
          if (!status) return;
          if (auto ret = save_shape(filename, shape); !ret) {
            auto lock = std::lock_guard{mutex};
            status    = ret;
          }
        });
    return {};
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
// YAML SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

static bool make_image_preset(
    image<vec4f>& hdr, image<vec4b>& ldr, const string& type) {
  if (type.find("sky") == type.npos) {
    auto imgf = make_image_preset(type);
    if (imgf.empty()) return false;
    if (type.find("-normal") == type.npos &&
        type.find("-displacement") == type.npos) {
      ldr = rgb_to_srgbb(imgf);
    } else {
      ldr = float_to_byte(imgf);
    }
    return true;
  } else {
    hdr = make_image_preset(type);
    return true;
  }
}

sceneio_status load_yaml(const string& filename, sceneio_model& scene) {
  // open file
  auto yaml = yaml_model{};
  if (auto ret = load_yaml(filename, yaml); !ret) return {ret.error};

  auto tmap = unordered_map<string, int>{{"", -1}};
  auto vmap = unordered_map<string, int>{{"", -1}};
  auto mmap = unordered_map<string, int>{{"", -1}};
  auto smap = unordered_map<string, int>{{"", -1}};

  // parse yaml reference
  auto get_yaml_ref = [](const yaml_element& yelment, const string& name,
                          int&                              value,
                          const unordered_map<string, int>& refs) -> bool {
    auto ref = ""s;
    if (!get_yaml_value(yelment, name, ref)) return false;
    if (ref == "") {
      value = -1;
      return true;
    } else {
      if (refs.find(ref) == refs.end()) return false;
      value = refs.at(ref);
      return true;
    }
  };

  // cameras
  for (auto& yelement : yaml.elements) {
    if (yelement.name == "cameras") {
      auto& camera = scene.cameras.emplace_back();
      if (!get_yaml_value(yelement, "name", camera.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "uri", camera.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "frame", camera.frame))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "orthographic", camera.orthographic))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "lens", camera.lens))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "aspect", camera.aspect))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "film", camera.film))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "focus", camera.focus))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "aperture", camera.aperture))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "uri")) {
        auto uri = ""s;
        if (!get_yaml_value(yelement, "uri", uri))
          return {filename + ": parse error"};
        camera.name = get_basename(uri);
      }
      if (has_yaml_value(yelement, "lookat")) {
        auto lookat = identity3x3f;
        if (!get_yaml_value(yelement, "lookat", lookat))
          return {filename + ": parse error"};
        camera.frame = lookat_frame(lookat.x, lookat.y, lookat.z);
        camera.focus = length(lookat.x - lookat.y);
      }
    } else if (yelement.name == "textures") {
      auto& texture = scene.textures.emplace_back();
      if (!get_yaml_value(yelement, "name", texture.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "filename", texture.filename))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "preset")) {
        auto preset = ""s;
        if (!get_yaml_value(yelement, "preset", preset))
          return {filename + ": parse error"};
        make_image_preset(texture.hdr, texture.ldr, preset);
        if (texture.filename.empty()) {
          texture.filename = "textures/ypreset-" + preset +
                             (texture.hdr.empty() ? ".png" : ".hdr");
        }
      }
      if (has_yaml_value(yelement, "uri")) {
        if (!get_yaml_value(yelement, "uri", texture.filename))
          return {filename + ": parse error"};
        texture.name           = get_basename(texture.filename);
        tmap[texture.filename] = (int)scene.textures.size() - 1;
      }
      tmap[texture.name] = (int)scene.textures.size() - 1;
    } else if (yelement.name == "materials") {
      auto& material = scene.materials.emplace_back();
      if (!get_yaml_value(yelement, "name", material.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "emission", material.emission))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "diffuse", material.diffuse))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "metallic", material.metallic))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "specular", material.specular))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "roughness", material.roughness))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "coat", material.coat))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "transmission", material.transmission))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "refract", material.refract))
        return {filename + ": parse error"};
      if (!get_yaml_value(
              yelement, "voltransmission", material.voltransmission))
        return {filename + ": parse error"};
      if (!get_yaml_value(
              yelement, "volmeanfreepath", material.volmeanfreepath))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "volscatter", material.volscatter))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "volemission", material.volemission))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "volanisotropy", material.volanisotropy))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "volscale", material.volscale))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "opacity", material.opacity))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "coat", material.coat))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "emission_tex", material.emission_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "diffuse_tex", material.diffuse_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "metallic_tex", material.metallic_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "specular_tex", material.specular_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(
              yelement, "transmission_tex", material.transmission_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(
              yelement, "roughness_tex", material.roughness_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(
              yelement, "subsurface_tex", material.subsurface_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "normal_tex", material.normal_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "normal_tex", material.normal_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "gltf_textures", material.gltf_textures))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "uri")) {
        if (!get_yaml_value(yelement, "uri", material.name))
          return {filename + ": parse error"};
        mmap[material.name] = (int)scene.materials.size() - 1;
        material.name       = get_basename(material.name);
      }
      mmap[material.name] = (int)scene.materials.size() - 1;
    } else if (yelement.name == "shapes") {
      auto& shape = scene.shapes.emplace_back();
      if (!get_yaml_value(yelement, "name", shape.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "filename", shape.filename))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "subdivisions", shape.subdivisions))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "catmullclark", shape.catmullclark))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "smooth", shape.smooth))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "facevarying", shape.facevarying))
        return {filename + ": parse error"};
      if (!get_yaml_ref(
              yelement, "displacement_tex", shape.displacement_tex, tmap))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "displacement", shape.displacement))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "uri")) {
        if (!get_yaml_value(yelement, "uri", shape.filename))
          return {filename + ": parse error"};
        shape.name           = get_basename(shape.filename);
        smap[shape.filename] = (int)scene.shapes.size() - 1;
      }
      if (has_yaml_value(yelement, "preset")) {
        auto preset = ""s;
        if (!get_yaml_value(yelement, "preset", preset))
          return {filename + ": parse error"};
        make_shape_preset(shape.points, shape.lines, shape.triangles,
            shape.quads, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords, shape.colors,
            shape.radius, preset);
        if (shape.filename.empty()) {
          shape.filename = "shapes/ypreset-" + preset + ".yvol";
        }
      }
      smap[shape.name] = (int)scene.shapes.size() - 1;
    } else if (yelement.name == "instances") {
      auto& instance = scene.instances.emplace_back();
      if (!get_yaml_value(yelement, "name", instance.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "frame", instance.frame))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "shape", instance.shape, smap))
        return {filename + ": parse error"};
      if (!get_yaml_ref(yelement, "material", instance.material, mmap))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "uri")) {
        auto uri = ""s;
        if (!get_yaml_value(yelement, "uri", uri))
          return {filename + ": parse error"};
        instance.name = get_basename(uri);
      }
      if (has_yaml_value(yelement, "lookat")) {
        auto lookat = identity3x3f;
        if (!get_yaml_value(yelement, "lookat", lookat))
          return {filename + ": parse error"};
        instance.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
    } else if (yelement.name == "environments") {
      auto& environment = scene.environments.emplace_back();
      if (!get_yaml_value(yelement, "name", environment.name))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "frame", environment.frame))
        return {filename + ": parse error"};
      if (!get_yaml_value(yelement, "emission", environment.emission))
        return {filename + ": parse error"};
      if (!get_yaml_ref(
              yelement, "emission_tex", environment.emission_tex, tmap))
        return {filename + ": parse error"};
      if (has_yaml_value(yelement, "uri")) {
        auto uri = ""s;
        if (!get_yaml_value(yelement, "uri", uri))
          return {filename + ": parse error"};
        environment.name = get_basename(uri);
      }
      if (has_yaml_value(yelement, "lookat")) {
        auto lookat = identity3x3f;
        if (!get_yaml_value(yelement, "lookat", lookat))
          return {filename + ": parse error"};
        environment.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
    }
  }

  return {};
}

// Save a scene in the builtin YAML format.
static sceneio_status load_yaml_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  scene = {};

  // Parse yaml
  if (auto ret = load_yaml(filename, scene); !ret) return ret;

  // load shape and textures
  if (auto ret = load_shapes(filename, scene, noparallel); !ret) return ret;
  if (auto ret = load_textures(filename, scene, noparallel); !ret) return ret;

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);
  trim_memory(scene);

  return {};
}

// Save yaml
static sceneio_status save_yaml(const string& filename,
    const sceneio_model& scene, bool ply_instances = false,
    const string& instances_name = "") {
  auto yaml = yaml_model{};

  for (auto stat : scene_stats(scene)) yaml.comments.push_back(stat);

  for (auto& camera : scene.cameras) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "cameras";
    add_yaml_value(yelement, "name", camera.name);
    add_yaml_value(yelement, "frame", camera.frame);
    if (camera.orthographic)
      add_yaml_value(yelement, "orthographic", camera.orthographic);
    add_yaml_value(yelement, "lens", camera.lens);
    add_yaml_value(yelement, "aspect", camera.aspect);
    add_yaml_value(yelement, "film", camera.film);
    add_yaml_value(yelement, "focus", camera.focus);
    add_yaml_value(yelement, "aperture", camera.aperture);
  }

  for (auto& texture : scene.textures) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "textures";
    add_yaml_value(yelement, "name", texture.name);
    add_yaml_value(yelement, "filename", texture.filename);
  }

  for (auto& material : scene.materials) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "materials";
    add_yaml_value(yelement, "name", material.name);
    add_yaml_value(yelement, "emission", material.emission);
    add_yaml_value(yelement, "diffuse", material.diffuse);
    add_yaml_value(yelement, "specular", material.specular);
    if (material.metallic)
      add_yaml_value(yelement, "metallic", material.metallic);
    if (material.transmission != zero3f)
      add_yaml_value(yelement, "transmission", material.transmission);
    add_yaml_value(yelement, "roughness", material.roughness);
    if (material.refract) add_yaml_value(yelement, "refract", material.refract);
    if (material.voltransmission != zero3f)
      add_yaml_value(yelement, "voltransmission", material.voltransmission);
    if (material.volmeanfreepath != zero3f)
      add_yaml_value(yelement, "volmeanfreepath", material.volmeanfreepath);
    if (material.volscatter != zero3f)
      add_yaml_value(yelement, "volscatter", material.volscatter);
    if (material.volemission != zero3f)
      add_yaml_value(yelement, "volemission", material.volemission);
    if (material.volanisotropy)
      add_yaml_value(yelement, "volanisotropy", material.volanisotropy);
    if (material.voltransmission != zero3f ||
        material.volmeanfreepath != zero3f)
      add_yaml_value(yelement, "volscale", material.volscale);
    if (material.coat != zero3f)
      add_yaml_value(yelement, "coat", material.coat);
    if (material.opacity != 1)
      add_yaml_value(yelement, "opacity", material.opacity);
    if (material.emission_tex >= 0)
      add_yaml_value(
          yelement, "emission_tex", scene.textures[material.emission_tex].name);
    if (material.diffuse_tex >= 0)
      add_yaml_value(
          yelement, "diffuse_tex", scene.textures[material.diffuse_tex].name);
    if (material.metallic_tex >= 0)
      add_yaml_value(
          yelement, "metallic_tex", scene.textures[material.metallic_tex].name);
    if (material.specular_tex >= 0)
      add_yaml_value(
          yelement, "specular_tex", scene.textures[material.specular_tex].name);
    if (material.roughness_tex >= 0)
      add_yaml_value(yelement, "roughness_tex",
          scene.textures[material.roughness_tex].name);
    if (material.transmission_tex >= 0)
      add_yaml_value(yelement, "transmission_tex",
          scene.textures[material.transmission_tex].name);
    if (material.subsurface_tex >= 0)
      add_yaml_value(yelement, "subsurface_tex",
          scene.textures[material.subsurface_tex].name);
    if (material.coat_tex >= 0)
      add_yaml_value(
          yelement, "coat_tex", scene.textures[material.coat_tex].name);
    if (material.opacity_tex >= 0)
      add_yaml_value(
          yelement, "opacity_tex", scene.textures[material.opacity_tex].name);
    if (material.normal_tex >= 0)
      add_yaml_value(
          yelement, "normal_tex", scene.textures[material.normal_tex].name);
    if (material.gltf_textures)
      add_yaml_value(yelement, "gltf_textures", material.gltf_textures);
  }

  for (auto& shape : scene.shapes) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "shapes";
    add_yaml_value(yelement, "name", shape.name);
    add_yaml_value(yelement, "filename", shape.filename);
    add_yaml_value(yelement, "subdivisions", shape.subdivisions);
    add_yaml_value(yelement, "catmullclark", shape.catmullclark);
    add_yaml_value(yelement, "smooth", shape.smooth);
    if (shape.facevarying)
      add_yaml_value(yelement, "facevarying", shape.facevarying);
    if (shape.displacement_tex >= 0)
      add_yaml_value(yelement, "displacement_tex",
          scene.textures[shape.displacement_tex].name);
    if (shape.displacement_tex >= 0)
      add_yaml_value(yelement, "displacement", shape.displacement);
  }

  for (auto& instance : scene.instances) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "instances";
    add_yaml_value(yelement, "name", instance.name);
    add_yaml_value(yelement, "frame", instance.frame);
    if (instance.shape >= 0)
      add_yaml_value(yelement, "shape", scene.shapes[instance.shape].name);
    if (instance.material >= 0)
      add_yaml_value(
          yelement, "material", scene.materials[instance.material].name);
  }

  for (auto& environment : scene.environments) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "environments";
    add_yaml_value(yelement, "name", environment.name);
    add_yaml_value(yelement, "frame", environment.frame);
    add_yaml_value(yelement, "emission", environment.emission);
    if (environment.emission_tex >= 0)
      add_yaml_value(yelement, "emission_tex",
          scene.textures[environment.emission_tex].name);
  }

  if (auto ret = save_yaml(filename, yaml); !ret) return {ret.error};

  return {};
}

// Save a scene in the builtin YAML format.
static sceneio_status save_yaml_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // save yaml file
  if (auto ret = save_yaml(filename, scene); !ret) return ret;
  if (auto ret = save_shapes(filename, scene, noparallel); !ret) return ret;
  if (auto ret = save_textures(filename, scene, noparallel); !ret) return ret;

  return {};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static sceneio_status load_obj(
    const string& filename, sceneio_model& scene, bool facevarying) {
  // load obj
  auto obj = obj_model{};
  if (auto ret = load_obj(filename, obj, false, true, true); !ret)
    return {ret.error};

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
    if (!oshape.faces.empty() && !facevarying && !has_quads) {
      get_obj_triangles(obj, oshape, shape.triangles, shape.positions,
          shape.normals, shape.texcoords, materials, ematerials, true);
    } else if (!oshape.faces.empty() && !facevarying && has_quads) {
      get_obj_quads(obj, oshape, shape.quads, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.lines.empty()) {
      get_obj_lines(obj, oshape, shape.lines, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.points.empty()) {
      get_obj_points(obj, oshape, shape.points, shape.positions, shape.normals,
          shape.texcoords, materials, ematerials, true);
    } else if (!oshape.faces.empty() && facevarying) {
      get_obj_fvquads(obj, oshape, shape.quadspos, shape.quadsnorm,
          shape.quadstexcoord, shape.positions, shape.normals, shape.texcoords,
          materials, ematerials, true);
    } else {
      return {filename + ": empty shape"};
    }
    // get material
    if (oshape.materials.size() != 1) {
      return {filename + ": missing material for " + oshape.name};
    }
    if (material_map.find(oshape.materials.at(0)) == material_map.end()) {
      return {filename + ": missing material " + oshape.materials.at(0)};
    }
    auto material = material_map.at(oshape.materials.at(0));
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

  return {};
}

// Loads an OBJ
static sceneio_status load_obj_scene(const string& filename,
    sceneio_model& scene, bool facevarying, bool noparallel) {
  scene = {};

  // Parse obj
  if (auto ret = load_obj(filename, scene, facevarying); !ret) return ret;
  if (auto ret = load_textures(filename, scene, noparallel); !ret) return ret;

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);

  return {};
}

static sceneio_status save_obj(
    const string& filename, const sceneio_model& scene, bool instances) {
  auto obj = obj_model{};

  for (auto stat : scene_stats(scene)) obj.comments.push_back(stat);

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
  if (instances) {
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
        return {filename + ": do not support empty shapes"};
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
  if (auto ret = save_obj(filename, obj); !ret) return {ret.error};

  return {};
}

static sceneio_status save_obj_scene(const string& filename,
    const sceneio_model& scene, bool instances, bool noparallel) {
  if (auto ret = save_obj(filename, scene, instances); !ret) return ret;
  if (auto ret = save_textures(filename, scene, noparallel)) return ret;
  return {};
}

void print_obj_camera(const sceneio_camera& camera) {
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

static sceneio_status load_ply_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  scene = {};

  // load ply mesh
  scene.shapes.push_back({});
  auto& shape    = scene.shapes.back();
  shape.name     = "shape";
  shape.filename = get_filename(filename);
  if (auto ret = load_shape(filename, shape.points, shape.lines,
          shape.triangles, shape.quads, shape.positions, shape.normals,
          shape.texcoords, shape.colors, shape.radius);
      !ret) {
    return {ret.error};
  }

  // add instance
  auto instance  = sceneio_instance{};
  instance.name  = shape.name;
  instance.shape = 0;
  scene.instances.push_back(instance);

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);

  return {};
}

static sceneio_status save_ply_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  if (scene.shapes.empty()) return {filename + ": empty scene"};
  auto& shape = scene.shapes.front();
  if (shape.quadspos.empty()) {
    save_shape(filename, shape.points, shape.lines, shape.triangles,
        shape.quads, shape.positions, shape.normals, shape.texcoords,
        shape.colors, shape.radius);
  } else {
    save_fvshape(filename, shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
        shape.positions, shape.normals, shape.texcoords);
  }
  return {};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// convert gltf to scene
static sceneio_status load_gltf(const string& filename, sceneio_model& scene) {
  auto gltf = gltf_model{};
  if (auto ret = load_gltf(filename, gltf); !ret) return {ret.error};

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
  auto cameras = vector<sceneio_camera>{};
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

  return {};
}

// Load a scene
static sceneio_status load_gltf_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // initialization
  scene = {};

  // load gltf
  if (auto ret = load_gltf(filename, scene); !ret) return ret;
  if (auto ret = load_textures(filename, scene, noparallel)) return ret;

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

  return {};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

static sceneio_status load_pbrt(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // load pbrt
  auto pbrt = pbrt_model{};
  if (auto ret = load_pbrt(filename, pbrt); !ret) return {ret.error};

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
    if (pshape.instance_frames.empty()) {
      auto& instance    = scene.instances.emplace_back();
      instance.name     = shape.name;
      instance.frame    = pshape.frame;
      instance.material = arealight_id >= 0 ? arealight_id : material_id;
      instance.shape    = (int)scene.shapes.size() - 1;
    } else {
      auto instance_id = 0;
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

  return {};
}

// load pbrt scenes
static sceneio_status load_pbrt_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  scene = sceneio_model{};

  // Parse pbrt
  if (auto ret = load_pbrt(filename, scene, noparallel); !ret) return ret;
  if (auto ret = load_shapes(filename, scene, noparallel); !ret) return ret;
  if (auto ret = load_textures(filename, scene, noparallel); !ret) return ret;

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_materials(scene);
  add_radius(scene);

  return {};
}

// Convert a scene to pbrt format
static sceneio_status save_pbrt(
    const string& filename, const sceneio_model& scene) {
  auto pbrt = pbrt_model{};

  // embed data
  for (auto stat : scene_stats(scene)) pbrt.comments.push_back(stat);

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

  if (auto ret = save_pbrt(filename, pbrt); !ret) return {ret.error};

  return {};
}

// Save a pbrt scene
sceneio_status save_pbrt_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // save pbrt
  if (auto ret = save_pbrt(filename, scene); !ret) return ret;

  // save meshes
  auto dirname = get_dirname(filename);
  for (auto& shape : scene.shapes) {
    if (shape.quadspos.empty()) {
      if (auto ret = save_shape(
              replace_extension(dirname + shape.filename, ".ply"), shape.points,
              shape.lines, shape.triangles, shape.quads, shape.positions,
              shape.normals, shape.texcoords, shape.colors, shape.radius);
          !ret) {
        return {filename + ": missing shape (" + ret.error + ")"};
      }
    } else {
      if (auto ret = save_fvshape(
              replace_extension(dirname + shape.filename, ".ply"),
              shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
              shape.positions, shape.normals, shape.texcoords);
          !ret) {
        return {filename + ": missing shape (" + ret.error + ")"};
      }
    }
  }

  // save textures
  if (auto ret = save_textures(filename, scene, noparallel); !ret) return ret;

  return {};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox_scene(sceneio_model& scene) {
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
