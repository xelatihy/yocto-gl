//
// Implementation for Yocto/Scene Input and Output functions.
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

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "ext/cgltf.h"
#include "yocto_color.h"
#include "yocto_commonio.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_parallel.h"
#include "yocto_sceneio.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::deque;
using std::unique_ptr;
using namespace std::string_literals;

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
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> scene_stats(const scene_scene* scene, bool verbose) {
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
  stats.push_back("cameras:      " + format(scene->cameras.size()));
  stats.push_back("shapes:       " + format(scene->shapes.size()));
  stats.push_back("environments: " + format(scene->environments.size()));
  stats.push_back("textures:     " + format(scene->textures.size()));
  stats.push_back(
      "points:       " + format(accumulate(scene->shapes,
                             [](auto& shape) { return shape.points.size(); })));
  stats.push_back(
      "lines:        " + format(accumulate(scene->shapes,
                             [](auto& shape) { return shape.lines.size(); })));
  stats.push_back("triangles:    " +
                  format(accumulate(scene->shapes,
                      [](auto& shape) { return shape.triangles.size(); })));
  stats.push_back(
      "quads:        " + format(accumulate(scene->shapes,
                             [](auto& shape) { return shape.quads.size(); })));
  stats.push_back("fvquads:     " +
                  format(accumulate(scene->shapes,
                      [](auto& shape) { return shape.quadspos.size(); })));
  stats.push_back(
      "texels4b:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->ldr.width() * (size_t)texture->ldr.width();
      })));
  stats.push_back(
      "texels4f:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->hdr.width() * (size_t)texture->hdr.height();
      })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene->
vector<string> scene_validation(const scene_scene* scene, bool notextures) {
  auto errs        = vector<string>();
  auto check_names = [&errs](const vector<string>& names, const string& base) {
    auto used = unordered_map<string, int>();
    used.reserve(names.size());
    for (auto& name : names) used[name] += 1;
    for (auto& [name, used] : used) {
      if (name.empty()) {
        errs.push_back("empty " + base + " name");
      } else if (used > 1) {
        errs.push_back("duplicated " + base + " name " + name);
      }
    }
  };
  auto check_empty_textures = [&errs](const scene_scene* scene) {
    for (auto idx = 0; idx < (int)scene->textures.size(); idx++) {
      auto& texture = scene->textures[idx];
      if (texture->hdr.empty() && texture->ldr.empty()) {
        errs.push_back("empty texture " + scene->texture_names[idx]);
      }
    }
  };

  check_names(scene->camera_names, "camera");
  check_names(scene->shape_names, "shape");
  check_names(scene->material_names, "material");
  check_names(scene->instance_names, "instance");
  check_names(scene->texture_names, "texture");
  check_names(scene->environment_names, "environment");
  if (!notextures) check_empty_textures(scene);

  return errs;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

scene_scene::~scene_scene() {
  for (auto texture : textures) delete texture;
}

// add an element
template <typename T>
static element_handle add_element(vector<T*>& elements, vector<string>& names,
    unordered_map<const T*, string>& name_map, const string& name,
    const string& base) {
  names.push_back(
      !name.empty() ? name : (base + std::to_string(elements.size())));
  auto element      = elements.emplace_back(new T{});
  name_map[element] = name;
  return (int)elements.size() - 1;
}
template <typename T>
static element_handle add_element(vector<T>& elements, vector<string>& names,
    const string& name, const string& base) {
  names.push_back(
      !name.empty() ? name : (base + std::to_string(elements.size())));
  elements.emplace_back();
  return (int)elements.size() - 1;
}

// add element
camera_handle add_camera(scene_scene* scene, const string& name) {
  return add_element(scene->cameras, scene->camera_names, name, "camera");
}
environment_handle add_environment(scene_scene* scene, const string& name) {
  return add_element(
      scene->environments, scene->environment_names, name, "environment");
}
shape_handle add_shape(scene_scene* scene, const string& name) {
  return add_element(scene->shapes, scene->shape_names, name, "shape");
}
texture_handle add_texture(scene_scene* scene, const string& name) {
  return add_element(scene->textures, scene->texture_names, scene->texture_map,
      name, "texture");
}
instance_handle add_instance(scene_scene* scene, const string& name) {
  return add_element(scene->instances, scene->instance_names, name, "instance");
}
material_handle add_material(scene_scene* scene, const string& name) {
  return add_element(scene->materials, scene->material_names, name, "material");
}
instance_handle add_complete_instance(scene_scene* scene, const string& name) {
  auto  handle      = add_instance(scene, name);
  auto& instance    = get_instance(scene, handle);
  instance.shape    = add_shape(scene, name);
  instance.material = add_material(scene, name);
  return handle;
}

// get element from a scene
scene_camera& get_camera(scene_scene* scene, camera_handle handle) {
  return scene->cameras.at(handle);
}
scene_environment& get_environment(
    scene_scene* scene, environment_handle handle) {
  return scene->environments.at(handle);
}
scene_instance& get_instance(scene_scene* scene, instance_handle handle) {
  return scene->instances.at(handle);
}
scene_material& get_material(scene_scene* scene, material_handle handle) {
  return scene->materials.at(handle);
}
scene_shape& get_shape(scene_scene* scene, shape_handle handle) {
  return scene->shapes.at(handle);
}
scene_texture* get_texture(scene_scene* scene, texture_handle handle) {
  return scene->textures.at(handle);
}
scene_instance& get_complete_instance(
    scene_scene* scene, instance_handle handle) {
  return scene->instances.at(handle);
}

// get element from a scene
const scene_camera& get_camera(const scene_scene* scene, camera_handle handle) {
  return scene->cameras.at(handle);
}
const scene_environment& get_environment(
    const scene_scene* scene, environment_handle handle) {
  return scene->environments.at(handle);
}
const scene_instance& get_instance(
    const scene_scene* scene, instance_handle handle) {
  return scene->instances.at(handle);
}
const scene_material& get_material(
    const scene_scene* scene, material_handle handle) {
  return scene->materials.at(handle);
}
const scene_shape& get_shape(const scene_scene* scene, shape_handle handle) {
  return scene->shapes.at(handle);
}
const scene_texture* get_texture(
    const scene_scene* scene, texture_handle handle) {
  return scene->textures.at(handle);
}
const scene_instance& get_complete_instance(
    const scene_scene* scene, instance_handle handle) {
  return scene->instances.at(handle);
}

// get name
string get_camera_name(const scene_scene* scene, int idx) {
  return scene->camera_names[idx];
}
string get_environment_name(const scene_scene* scene, int idx) {
  return scene->environment_names[idx];
}
string get_shape_name(const scene_scene* scene, int idx) {
  return scene->shape_names[idx];
}
string get_texture_name(const scene_scene* scene, int idx) {
  return scene->texture_names[idx];
}
string get_instance_name(const scene_scene* scene, int idx) {
  return scene->instance_names[idx];
}
string get_material_name(const scene_scene* scene, int idx) {
  return scene->material_names[idx];
}

string get_camera_name(const scene_scene* scene, const scene_camera& camera) {
  return scene->camera_names.at(&camera - scene->cameras.data());
}
string get_environment_name(
    const scene_scene* scene, const scene_environment& environment) {
  return scene->environment_names.at(&environment - scene->environments.data());
}
string get_shape_name(const scene_scene* scene, const scene_shape& shape) {
  return scene->shape_names.at(&shape - scene->shapes.data());
}
string get_texture_name(
    const scene_scene* scene, const scene_texture* texture) {
  return scene->texture_map.at(texture);
}
string get_instance_name(
    const scene_scene* scene, const scene_instance& instance) {
  return scene->instance_names.at(&instance - scene->instances.data());
}
string get_material_name(
    const scene_scene* scene, const scene_material& material) {
  return scene->material_names.at(&material - scene->materials.data());
}

// Add missing cameras.
void add_cameras(scene_scene* scene) {
  if (!scene->cameras.empty()) return;
  auto& camera        = get_camera(scene, add_camera(scene, "camera"));
  camera.orthographic = false;
  camera.film         = 0.036;
  camera.aspect       = (float)16 / (float)9;
  camera.aperture     = 0;
  camera.lens         = 0.050;
  auto bbox           = compute_bounds(scene);
  auto center         = (bbox.max + bbox.min) / 2;
  auto bbox_radius    = length(bbox.max - bbox.min) / 2;
  auto camera_dir     = vec3f{0, 0, 1};
  auto camera_dist = bbox_radius * camera.lens / (camera.film / camera.aspect);
  camera_dist *= 2.0f;  // correction for tracer camera implementation
  auto from    = camera_dir * camera_dist + center;
  auto to      = center;
  auto up      = vec3f{0, 1, 0};
  camera.frame = lookat_frame(from, to, up);
  camera.focus = length(from - to);
}

// Add missing radius.
void add_radius(scene_scene* scene, float radius) {
  for (auto& shape : scene->shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add missing materials.
void add_materials(scene_scene* scene) {
  auto default_material = invalid_handle;
  for (auto& instance : scene->instances) {
    if (instance.material != invalid_handle) continue;
    if (default_material == invalid_handle) {
      default_material = add_material(scene);
      auto& material   = get_material(scene, default_material);
      material.color   = {0.8, 0.8, 0.8};
    }
    instance.material = default_material;
  }
}

// Add a sky environment
void add_sky(scene_scene* scene, float sun_angle) {
  auto thandle         = add_texture(scene, "sky");
  auto texture         = get_texture(scene, thandle);
  texture->hdr         = make_sunsky({1024, 512}, sun_angle);
  auto& environment    = get_environment(scene, add_environment(scene, "sky"));
  environment.emission = {1, 1, 1};
  environment.emission_tex = thandle;
}

// get named camera or default if camera is empty
scene_camera& get_camera(scene_scene* scene, const string& name) {
  return get_camera(scene, get_camera_handle(scene, name));
}
const scene_camera& get_camera(const scene_scene* scene, const string& name) {
  return get_camera(scene, get_camera_handle(scene, name));
}
camera_handle get_camera_handle(const scene_scene* scene, const string& name) {
  if (scene->cameras.empty()) return invalid_handle;
  for (auto idx = 0; idx < (int)scene->camera_names.size(); idx++) {
    if (get_camera_name(scene, idx) == name) return idx;
  }
  for (auto idx = 0; idx < (int)scene->camera_names.size(); idx++) {
    if (get_camera_name(scene, idx) == "default") return idx;
  }
  for (auto idx = 0; idx < (int)scene->camera_names.size(); idx++) {
    if (get_camera_name(scene, idx) == "camera") return idx;
  }
  for (auto idx = 0; idx < (int)scene->camera_names.size(); idx++) {
    if (get_camera_name(scene, idx) == "camera1") return idx;
  }
  return 0;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_scene* scene) {
  auto shape_bbox = vector<bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto& shape : scene->shapes) {
    auto sbvh = shape_bbox.emplace_back();
    for (auto p : shape.positions) sbvh = merge(sbvh, p);
  }
  for (auto& instance : scene->instances) {
    auto& sbvh = shape_bbox[instance.shape];
    bbox       = merge(bbox, transform_bbox(instance.frame, sbvh));
  }
  return bbox;
}

// Clone a scene
void clone_scene(scene_scene* dest, const scene_scene* scene) {
  dest->name      = scene->name;
  dest->copyright = scene->copyright;
  throw std::runtime_error("not implemented yet");
}

// Reduce memory usage
void trim_memory(scene_scene* scene) {
  for (auto& shape : scene->shapes) {
    shape.points.shrink_to_fit();
    shape.lines.shrink_to_fit();
    shape.triangles.shrink_to_fit();
    shape.quads.shrink_to_fit();
    shape.positions.shrink_to_fit();
    shape.normals.shrink_to_fit();
    shape.texcoords.shrink_to_fit();
    shape.colors.shrink_to_fit();
    shape.radius.shrink_to_fit();
    shape.tangents.shrink_to_fit();
    shape.quadspos.shrink_to_fit();
    shape.quadsnorm.shrink_to_fit();
    shape.quadstexcoord.shrink_to_fit();
  }
  for (auto texture : scene->textures) {
    texture->hdr.shrink_to_fit();
    texture->ldr.shrink_to_fit();
  }
  scene->cameras.shrink_to_fit();
  scene->shapes.shrink_to_fit();
  scene->textures.shrink_to_fit();
  scene->environments.shrink_to_fit();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

void tesselate_shape(scene_scene* scene, scene_shape& shape) {
  if (shape.subdivisions > 0) {
    if (!shape.points.empty()) {
      throw std::runtime_error("cannot subdivide points");
    } else if (!shape.lines.empty()) {
      std::tie(std::ignore, shape.texcoords) = subdivide_lines(
          shape.lines, shape.texcoords, shape.subdivisions);
      std::tie(std::ignore, shape.normals) = subdivide_lines(
          shape.lines, shape.normals, shape.subdivisions);
      std::tie(std::ignore, shape.colors) = subdivide_lines(
          shape.lines, shape.colors, shape.subdivisions);
      std::tie(std::ignore, shape.radius) = subdivide_lines(
          shape.lines, shape.radius, shape.subdivisions);
      std::tie(shape.lines, shape.positions) = subdivide_lines(
          shape.lines, shape.positions, shape.subdivisions);
    } else if (!shape.triangles.empty()) {
      std::tie(std::ignore, shape.texcoords) = subdivide_triangles(
          shape.triangles, shape.texcoords, shape.subdivisions);
      std::tie(std::ignore, shape.normals) = subdivide_triangles(
          shape.triangles, shape.normals, shape.subdivisions);
      std::tie(std::ignore, shape.colors) = subdivide_triangles(
          shape.triangles, shape.colors, shape.subdivisions);
      std::tie(std::ignore, shape.radius) = subdivide_triangles(
          shape.triangles, shape.radius, shape.subdivisions);
      std::tie(shape.triangles, shape.positions) = subdivide_triangles(
          shape.triangles, shape.positions, shape.subdivisions);
    } else if (!shape.quads.empty()) {
      if (shape.catmullclark) {
        std::tie(std::ignore, shape.texcoords) = subdivide_catmullclark(
            shape.quads, shape.texcoords, shape.subdivisions, true);
        std::tie(std::ignore, shape.normals) = subdivide_catmullclark(
            shape.quads, shape.normals, shape.subdivisions, true);
        std::tie(std::ignore, shape.colors) = subdivide_catmullclark(
            shape.quads, shape.colors, shape.subdivisions);
        std::tie(std::ignore, shape.radius) = subdivide_catmullclark(
            shape.quads, shape.radius, shape.subdivisions);
        std::tie(std::ignore, shape.positions) = subdivide_catmullclark(
            shape.quads, shape.positions, shape.subdivisions);
      } else {
        std::tie(std::ignore, shape.texcoords) = subdivide_quads(
            shape.quads, shape.texcoords, shape.subdivisions);
        std::tie(std::ignore, shape.normals) = subdivide_quads(
            shape.quads, shape.normals, shape.subdivisions);
        std::tie(std::ignore, shape.colors) = subdivide_quads(
            shape.quads, shape.colors, shape.subdivisions);
        std::tie(std::ignore, shape.radius) = subdivide_quads(
            shape.quads, shape.radius, shape.subdivisions);
        std::tie(shape.quads, shape.positions) = subdivide_quads(
            shape.quads, shape.positions, shape.subdivisions);
      }
    } else if (!shape.quadspos.empty()) {
      if (shape.catmullclark) {
        std::tie(shape.quadstexcoord, shape.texcoords) = subdivide_catmullclark(
            shape.quadstexcoord, shape.texcoords, shape.subdivisions, true);
        std::tie(shape.quadsnorm, shape.normals) = subdivide_catmullclark(
            shape.quadsnorm, shape.normals, shape.subdivisions, true);
        std::tie(shape.quadspos, shape.positions) = subdivide_catmullclark(
            shape.quadspos, shape.positions, shape.subdivisions);
      } else {
        std::tie(shape.quadstexcoord, shape.texcoords) = subdivide_quads(
            shape.quadstexcoord, shape.texcoords, shape.subdivisions);
        std::tie(shape.quadsnorm, shape.normals) = subdivide_quads(
            shape.quadsnorm, shape.normals, shape.subdivisions);
        std::tie(shape.quadspos, shape.positions) = subdivide_quads(
            shape.quadspos, shape.positions, shape.subdivisions);
      }
      if (shape.smooth) {
        shape.normals   = quads_normals(shape.quadspos, shape.positions);
        shape.quadsnorm = shape.quadspos;
      } else {
        shape.normals   = {};
        shape.quadsnorm = {};
      }
    } else {
      throw std::runtime_error("not supported yet");
    }
    shape.subdivisions = 0;
  }

  if (shape.displacement != 0 && shape.displacement_tex != invalid_handle) {
    if (shape.texcoords.empty())
      throw std::runtime_error("missing texture coordinates");

    if (!shape.triangles.empty() || !shape.quads.empty()) {
      auto no_normals = shape.normals.empty();
      if (shape.normals.empty())
        shape.normals = !shape.triangles.empty()
                            ? triangles_normals(
                                  shape.triangles, shape.positions)
                            : quads_normals(shape.quads, shape.positions);
      for (auto idx = 0; idx < shape.positions.size(); idx++) {
        auto diplacement_tex = get_texture(scene, shape.displacement_tex);
        auto disp            = mean(
            eval_texture(diplacement_tex, shape.texcoords[idx], true));
        if (!diplacement_tex->ldr.empty()) disp -= 0.5f;
        shape.positions[idx] += shape.normals[idx] * shape.displacement * disp;
      }
      if (shape.smooth) {
        shape.normals = !shape.triangles.empty()
                            ? triangles_normals(
                                  shape.triangles, shape.positions)
                            : quads_normals(shape.quads, shape.positions);
      } else if (no_normals) {
        shape.normals = {};
      }
    } else if (!shape.quadspos.empty()) {
      // facevarying case
      auto offset = vector<float>(shape.positions.size(), 0);
      auto count  = vector<int>(shape.positions.size(), 0);
      for (auto fid = 0; fid < shape.quadspos.size(); fid++) {
        auto qpos = shape.quadspos[fid];
        auto qtxt = shape.quadstexcoord[fid];
        for (auto i = 0; i < 4; i++) {
          auto displacement_tex = get_texture(scene, shape.displacement_tex);
          auto disp             = mean(
              eval_texture(displacement_tex, shape.texcoords[qtxt[i]], true));
          if (!displacement_tex->ldr.empty()) disp -= 0.5f;
          offset[qpos[i]] += shape.displacement * disp;
          count[qpos[i]] += 1;
        }
      }
      auto normals = quads_normals(shape.quadspos, shape.positions);
      for (auto vid = 0; vid < shape.positions.size(); vid++) {
        shape.positions[vid] += normals[vid] * offset[vid] / count[vid];
      }
      if (shape.smooth || !shape.normals.empty()) {
        shape.quadsnorm = shape.quadspos;
        shape.normals   = quads_normals(shape.quadspos, shape.positions);
      }
    }

    shape.displacement     = 0;
    shape.displacement_tex = invalid_handle;
  }
  if (!shape.quadspos.empty()) {
    std::tie(shape.quads, shape.positions, shape.normals, shape.texcoords) =
        split_facevarying(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
            shape.positions, shape.normals, shape.texcoords);
    shape.points    = {};
    shape.lines     = {};
    shape.triangles = {};
    shape.colors    = {};
    shape.radius    = {};
  }
}

void tesselate_shapes(
    scene_scene* scene, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, (int)scene->shapes.size() + 1};
  if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);

  // tesselate shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
    tesselate_shape(scene, shape);
  }

  // done
  if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Check texture size
vec2i texture_size(const scene_texture* texture) {
  if (!texture->hdr.empty()) {
    return texture->hdr.imsize();
  } else if (!texture->ldr.empty()) {
    return texture->ldr.imsize();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
vec4f lookup_texture(
    const scene_texture* texture, const vec2i& ij, bool ldr_as_linear) {
  if (!texture->hdr.empty()) {
    return texture->hdr[ij];
  } else if (!texture->ldr.empty()) {
    return ldr_as_linear ? byte_to_float(texture->ldr[ij])
                         : srgb_to_rgb(byte_to_float(texture->ldr[ij]));
  } else {
    return {1, 1, 1, 1};
  }
}

// Evaluate a texture
vec4f eval_texture(const scene_texture* texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  // get texture
  if (texture == nullptr) return {1, 1, 1, 1};

  // get image width/height
  auto size = texture_size(texture);

  // get coordinates normalized for tiling
  auto s = 0.0f, t = 0.0f;
  if (clamp_to_edge) {
    s = clamp(uv.x, 0.0f, 1.0f) * size.x;
    t = clamp(uv.y, 0.0f, 1.0f) * size.y;
  } else {
    s = fmod(uv.x, 1.0f) * size.x;
    if (s < 0) s += size.x;
    t = fmod(uv.y, 1.0f) * size.y;
    if (t < 0) t += size.y;
  }

  // get image coordinates and residuals
  auto i = clamp((int)s, 0, size.x - 1), j = clamp((int)t, 0, size.y - 1);
  auto ii = (i + 1) % size.x, jj = (j + 1) % size.y;
  auto u = s - i, v = t - j;

  if (no_interpolation) return lookup_texture(texture, {i, j}, ldr_as_linear);

  // handle interpolation
  return lookup_texture(texture, {i, j}, ldr_as_linear) * (1 - u) * (1 - v) +
         lookup_texture(texture, {i, jj}, ldr_as_linear) * (1 - u) * v +
         lookup_texture(texture, {ii, j}, ldr_as_linear) * u * (1 - v) +
         lookup_texture(texture, {ii, jj}, ldr_as_linear) * u * v;
}

// Helpers
vec4f eval_texture(const scene_scene* scene, texture_handle texture,
    const vec2f& uv, bool ldr_as_linear, bool no_interpolation,
    bool clamp_to_edge) {
  if (texture == invalid_handle) return {1, 1, 1, 1};
  return eval_texture(
      get_texture(scene, texture), uv, ldr_as_linear, no_interpolation);
}

// Generates a ray from a camera for yimg::image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera(
    const scene_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  if (!camera.orthographic) {
    auto q = vec3f{film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f),
        camera.lens};
    // ray direction through the lens center
    auto dc = -normalize(q);
    // point on the lens
    auto e = vec3f{
        lens_uv.x * camera.aperture / 2, lens_uv.y * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = dc * camera.focus / abs(dc.z);
    // correct ray direction to account for camera focusing
    auto d = normalize(p - e);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  } else {
    auto scale = 1 / camera.lens;
    auto q     = vec3f{film.x * (0.5f - image_uv.x) * scale,
        film.y * (image_uv.y - 0.5f) * scale, camera.lens};
    // point on the lens
    auto e = vec3f{-q.x, -q.y, 0} + vec3f{lens_uv.x * camera.aperture / 2,
                                        lens_uv.y * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = vec3f{-q.x, -q.y, -camera.focus};
    // correct ray direction to account for camera focusing
    auto d = normalize(p - e);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  }
}

// Eval position
vec3f eval_position(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = get_shape(scene, instance.shape);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_point(
        instance.frame, interpolate_triangle(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z], uv));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_point(instance.frame,
        interpolate_quad(shape.positions[q.x], shape.positions[q.y],
            shape.positions[q.z], shape.positions[q.w], uv));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_point(instance.frame,
        interpolate_line(shape.positions[l.x], shape.positions[l.y], uv.x));
  } else if (!shape.points.empty()) {
    return transform_point(
        instance.frame, shape.positions[shape.points[element]]);
  } else {
    return zero3f;
  }
}

// Shape element normal.
vec3f eval_element_normal(
    const scene_scene* scene, const scene_instance& instance, int element) {
  auto& shape = get_shape(scene, instance.shape);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, triangle_normal(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z]));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_normal(
        instance.frame, quad_normal(shape.positions[q.x], shape.positions[q.y],
                            shape.positions[q.z], shape.positions[q.w]));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_normal(instance.frame,
        line_tangent(shape.positions[l.x], shape.positions[l.y]));
  } else if (!shape.points.empty()) {
    return {0, 0, 1};
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
vec3f eval_normal(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = get_shape(scene, instance.shape);
  if (shape.normals.empty())
    return eval_element_normal(scene, instance, element);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, normalize(interpolate_triangle(shape.normals[t.x],
                            shape.normals[t.y], shape.normals[t.z], uv)));
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return transform_normal(instance.frame,
        normalize(interpolate_quad(shape.normals[q.x], shape.normals[q.y],
            shape.normals[q.z], shape.normals[q.w], uv)));
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return transform_normal(instance.frame,
        normalize(
            interpolate_line(shape.normals[l.x], shape.normals[l.y], uv.x)));
  } else if (!shape.points.empty()) {
    return transform_normal(
        instance.frame, normalize(shape.normals[shape.points[element]]));
  } else {
    return zero3f;
  }
}

// Eval texcoord
vec2f eval_texcoord(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = get_shape(scene, instance.shape);
  if (shape.texcoords.empty()) return uv;
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.texcoords[t.x], shape.texcoords[t.y], shape.texcoords[t.z], uv);
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return interpolate_quad(shape.texcoords[q.x], shape.texcoords[q.y],
        shape.texcoords[q.z], shape.texcoords[q.w], uv);
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return interpolate_line(shape.texcoords[l.x], shape.texcoords[l.y], uv.x);
  } else if (!shape.points.empty()) {
    return shape.texcoords[shape.points[element]];
  } else {
    return zero2f;
  }
}

#if 0
// Shape element normal.
static pair<vec3f, vec3f> eval_tangents(
    const trace_shape& shape, int element, const vec2f& uv) {
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    if (shape.texcoords.empty()) {
      return triangle_tangents_fromuv(shape.positions[t.x],
          shape.positions[t.y], shape.positions[t.z], {0, 0}, {1, 0}, {0, 1});
    } else {
      return triangle_tangents_fromuv(shape.positions[t.x],
          shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
          shape.texcoords[t.y], shape.texcoords[t.z]);
    }
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    if (shape.texcoords.empty()) {
      return quad_tangents_fromuv(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w], {0, 0}, {1, 0}, {0, 1},
          {1, 1}, uv);
    } else {
      return quad_tangents_fromuv(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w], shape.texcoords[q.x],
          shape.texcoords[q.y], shape.texcoords[q.z], shape.texcoords[q.w],
          uv);
    }
  } else {
    return {zero3f, zero3f};
  }
}
#endif

// Shape element normal.
pair<vec3f, vec3f> eval_element_tangents(
    const scene_scene* scene, const scene_instance& instance, int element) {
  auto& shape = get_shape(scene, instance.shape);
  if (!shape.triangles.empty() && !shape.texcoords.empty()) {
    auto t        = shape.triangles[element];
    auto [tu, tv] = triangle_tangents_fromuv(shape.positions[t.x],
        shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
        shape.texcoords[t.y], shape.texcoords[t.z]);
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else if (!shape.quads.empty() && !shape.texcoords.empty()) {
    auto q        = shape.quads[element];
    auto [tu, tv] = quad_tangents_fromuv(shape.positions[q.x],
        shape.positions[q.y], shape.positions[q.z], shape.positions[q.w],
        shape.texcoords[q.x], shape.texcoords[q.y], shape.texcoords[q.z],
        shape.texcoords[q.w], {0, 0});
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else {
    return {};
  }
}

vec3f eval_normalmap(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape    = get_shape(scene, instance.shape);
  auto& material = get_material(scene, instance.material);
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalid_handle &&
      (!shape.triangles.empty() || !shape.quads.empty())) {
    auto normal_tex = get_texture(scene, material.normal_tex);
    auto normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, true));
    auto [tu, tv]   = eval_element_tangents(scene, instance, element);
    auto frame      = frame3f{tu, tv, normal, zero3f};
    frame.x         = orthonormalize(frame.x, frame.z);
    frame.y         = normalize(cross(frame.z, frame.x));
    auto flip_v     = dot(frame.y, tv) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading normal
vec3f eval_shading_normal(const scene_scene* scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = get_shape(scene, instance.shape);
  auto& material = get_material(scene, instance.shape);
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalid_handle) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (!material.thin) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else if (!shape.lines.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    return orthonormalize(outgoing, normal);
  } else if (!shape.points.empty()) {
    return outgoing;
  } else {
    return zero3f;
  }
}

// Eval color
vec4f eval_color(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = get_shape(scene, instance.shape);
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.colors[t.x], shape.colors[t.y], shape.colors[t.z], uv);
  } else if (!shape.quads.empty()) {
    auto q = shape.quads[element];
    return interpolate_quad(shape.colors[q.x], shape.colors[q.y],
        shape.colors[q.z], shape.colors[q.w], uv);
  } else if (!shape.lines.empty()) {
    auto l = shape.lines[element];
    return interpolate_line(shape.colors[l.x], shape.colors[l.y], uv.x);
  } else if (!shape.points.empty()) {
    return shape.colors[shape.points[element]];
  } else {
    return {0, 0, 0, 0};
  }
}

// Evaluate environment color.
vec3f eval_environment(const scene_scene* scene,
    const scene_environment& environment, const vec3f& direction) {
  auto wl       = transform_direction(inverse(environment.frame), direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_scene* scene, const vec3f& direction) {
  auto emission = zero3f;
  for (auto environment : scene->environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

// Evaluate point
scene_material_sample eval_material(const scene_scene* scene,
    const scene_material& material, const vec2f& texcoord) {
  auto mat = scene_material_sample{};
  mat.emission =
      material.emission *
      xyz(eval_texture(scene, material.emission_tex, texcoord, false));
  mat.color = material.color *
              xyz(eval_texture(scene, material.color_tex, texcoord, false));
  mat.specular = material.specular *
                 eval_texture(scene, material.specular_tex, texcoord, true).x;
  mat.metallic = material.metallic *
                 eval_texture(scene, material.metallic_tex, texcoord, true).x;
  mat.roughness = material.roughness *
                  eval_texture(scene, material.roughness_tex, texcoord, true).x;
  mat.ior  = material.ior;
  mat.coat = material.coat *
             eval_texture(scene, material.coat_tex, texcoord, true).x;
  mat.transmission =
      material.transmission *
      eval_texture(scene, material.emission_tex, texcoord, true).x;
  mat.translucency =
      material.translucency *
      eval_texture(scene, material.translucency_tex, texcoord, true).x;
  mat.opacity = material.opacity *
                eval_texture(scene, material.opacity_tex, texcoord, true).x;
  mat.thin = material.thin || material.transmission == 0;
  mat.scattering =
      material.scattering *
      xyz(eval_texture(scene, material.scattering_tex, texcoord, false));
  mat.scanisotropy = material.scanisotropy;
  mat.trdepth      = material.trdepth;
  mat.normalmap    = material.normal_tex != invalid_handle
                         ? -1 + 2 * xyz(eval_texture(scene, material.normal_tex,
                                     texcoord, true))
                         : vec3f{0, 0, 1};
  return mat;
}

// constant values
static const auto coat_ior       = 1.5f;
static const auto coat_roughness = 0.03f * 0.03f;

// Eval material to obtain emission, brdf and opacity.
vec3f eval_emission(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv, const vec3f& normal, const vec3f& outgoing) {
  auto& material = get_material(scene, instance.material);
  auto  texcoord = eval_texcoord(scene, instance, element, uv);
  return material.emission *
         xyz(eval_texture(scene, material.emission_tex, texcoord));
}

// Eval material to obtain emission, brdf and opacity.
float eval_opacity(const scene_scene* scene, const scene_instance& instance,
    int element, const vec2f& uv, const vec3f& normal, const vec3f& outgoing) {
  auto& material = get_material(scene, instance.material);
  auto  texcoord = eval_texcoord(scene, instance, element, uv);
  auto  opacity  = material.opacity *
                 eval_texture(scene, material.opacity_tex, texcoord, true).x;
  if (opacity > 0.999f) opacity = 1;
  return opacity;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox(scene_scene* scene) {
  scene->name     = "cornellbox";
  auto& camera    = get_camera(scene, add_camera(scene));
  camera.frame    = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 3.9}};
  camera.lens     = 0.035;
  camera.aperture = 0.0;
  camera.focus    = 3.9;
  camera.film     = 0.024;
  camera.aspect   = 1;
  auto& floor     = get_instance(scene, add_complete_instance(scene, "floor"));
  get_shape(scene, floor.shape).positions = {
      {-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  get_shape(scene, floor.shape).triangles   = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, floor.material).color = {0.725, 0.71, 0.68};
  auto& ceiling = get_instance(scene, add_complete_instance(scene, "ceiling"));
  get_shape(scene, ceiling.shape).positions = {
      {-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  get_shape(scene, ceiling.shape).triangles   = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, ceiling.material).color = {0.725, 0.71, 0.68};
  auto& backwall                              = get_instance(
      scene, add_complete_instance(scene, "backwall"));
  get_shape(scene, backwall.shape).positions = {
      {-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  get_shape(scene, backwall.shape).triangles   = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, backwall.material).color = {0.725, 0.71, 0.68};
  auto& rightwall                              = get_instance(
      scene, add_complete_instance(scene, "rightwall"));
  get_shape(scene, rightwall.shape).positions = {
      {1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  get_shape(scene, rightwall.shape).triangles   = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, rightwall.material).color = {0.14, 0.45, 0.091};
  auto& leftwall                                = get_instance(
      scene, add_complete_instance(scene, "leftwall"));
  get_shape(scene, leftwall.shape).positions = {
      {-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  get_shape(scene, leftwall.shape).triangles   = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, leftwall.material).color = {0.63, 0.065, 0.05};
  auto& shortbox                               = get_instance(
      scene, add_complete_instance(scene, "shortbox"));
  get_shape(scene, shortbox.shape).positions = {{0.53, 0.6, 0.75},
      {0.7, 0.6, 0.17}, {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57},
      {-0.05, 0.0, 0.57}, {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0},
      {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75}, {0.53, 0.6, 0.75},
      {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57}, {0.7, 0.0, 0.17},
      {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75}, {0.13, 0.0, 0.0},
      {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17}, {0.53, 0.0, 0.75},
      {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0}, {-0.05, 0.0, 0.57}};
  get_shape(scene, shortbox.shape).triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6},
      {6, 7, 4}, {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12},
      {16, 17, 18}, {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  get_material(scene, shortbox.material).color = {0.725, 0.71, 0.68};
  auto& tallbox = get_instance(scene, add_complete_instance(scene, "tallbox"));
  get_shape(scene, tallbox.shape).positions = {{-0.53, 1.2, 0.09},
      {0.04, 1.2, -0.09}, {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49},
      {-0.53, 0.0, 0.09}, {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49},
      {-0.14, 1.2, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67},
      {-0.14, 1.2, -0.67}, {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 0.0, -0.09}, {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09},
      {-0.53, 0.0, 0.09}, {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09},
      {-0.14, 0.0, -0.67}, {-0.71, 0.0, -0.49}};
  get_shape(scene, tallbox.shape).triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6},
      {6, 7, 4}, {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12},
      {16, 17, 18}, {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  get_material(scene, tallbox.material).color = {0.725, 0.71, 0.68};
  auto& light = get_instance(scene, add_complete_instance(scene, "light"));
  get_shape(scene, light.shape).positions      = {{-0.25, 1.99, 0.25},
      {-0.25, 1.99, -0.25}, {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  get_shape(scene, light.shape).triangles      = {{0, 1, 2}, {2, 3, 0}};
  get_material(scene, light.material).emission = {17, 12, 4};
}

}  // namespace yocto
