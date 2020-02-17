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
// TODO: make_safe_name
// TODO: pbrt load/save ply
//

//
// TODO: update transforms -> should be compute transforms?
// TODO: update tesselation -> should be compute tesselation
// TODO: move out animation utilities
//

#include "yocto_sceneio.h"

#include <atomic>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <deque>
#include <future>
#include <memory>
using std::make_unique;

#include "ext/filesystem.hpp"
#include "ext/json.hpp"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_shape.h"
namespace fs = ghc::filesystem;

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
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> scene_stats(const sceneio_model* scene, bool verbose) {
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
  stats.push_back("subdivs:      " + format(scene->subdivs.size()));
  stats.push_back("environments: " + format(scene->environments.size()));
  stats.push_back("textures:     " + format(scene->textures.size()));
  stats.push_back(
      "points:       " + format(accumulate(scene->shapes,
                             [](auto shape) { return shape->points.size(); })));
  stats.push_back(
      "lines:        " + format(accumulate(scene->shapes,
                             [](auto shape) { return shape->lines.size(); })));
  stats.push_back("triangles:    " +
                  format(accumulate(scene->shapes,
                      [](auto shape) { return shape->triangles.size(); })));
  stats.push_back(
      "quads:        " + format(accumulate(scene->shapes,
                             [](auto shape) { return shape->quads.size(); })));
  stats.push_back("sfvquads:     " +
                  format(accumulate(scene->subdivs,
                      [](auto shape) { return shape->quadspos.size(); })));
  stats.push_back(
      "texels4b:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->ldr.size().x * (size_t)texture->ldr.size().x;
      })));
  stats.push_back(
      "texels4f:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->hdr.size().x * (size_t)texture->hdr.size().y;
      })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene->
vector<string> scene_validation(const sceneio_model* scene, bool notextures) {
  auto errs        = vector<string>();
  auto check_names = [&errs](const auto& vals, const string& base) {
    auto used = unordered_map<string, int>();
    used.reserve(vals.size());
    for (auto& value : vals) used[value->name] += 1;
    for (auto& [name, used] : used) {
      if (name == "") {
        errs.push_back("empty " + base + " name");
      } else if (used > 1) {
        errs.push_back("duplicated " + base + " name " + name);
      }
    }
  };
  auto check_empty_textures = [&errs](const vector<sceneio_texture*>& vals) {
    for (auto value : vals) {
      if (value->hdr.empty() && value->ldr.empty()) {
        errs.push_back("empty texture " + value->name);
      }
    }
  };

  check_names(scene->cameras, "camera");
  check_names(scene->shapes, "shape");
  check_names(scene->subdivs, "subdiv");
  // check_names(scene->shapes, "instance");
  // check_names(scene->shapes, "object");
  check_names(scene->textures, "texture");
  check_names(scene->environments, "environment");
  if (!notextures) check_empty_textures(scene->textures);

  return errs;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

sceneio_model::~sceneio_model() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto subdiv : subdivs) delete subdiv;
  for (auto material : materials) delete material;
  for (auto instance : instances) delete instance;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
}

unique_ptr<sceneio_model> make_sceneio_model() {
  return make_unique<sceneio_model>();
}

// add element
sceneio_camera* add_camera(sceneio_model* scene) {
  auto camera  = scene->cameras.emplace_back(new sceneio_camera{});
  camera->name = "cameras/camera" + std::to_string(scene->cameras.size()) +
                 ".json";
  return camera;
}
sceneio_environment* add_environment(sceneio_model* scene) {
  auto environment = scene->environments.emplace_back(
      new sceneio_environment{});
  environment->name = "environments/environment" +
                      std::to_string(scene->environments.size()) + ".json";
  return environment;
}
sceneio_shape* add_shape(sceneio_model* scene) {
  auto shape  = scene->shapes.emplace_back(new sceneio_shape{});
  shape->name = "shapes/shape" + std::to_string(scene->shapes.size()) + ".ply";
  return shape;
}
sceneio_subdiv* add_subdiv(sceneio_model* scene) {
  auto subdiv  = scene->subdivs.emplace_back(new sceneio_subdiv{});
  subdiv->name = "subdivs/subdiv" + std::to_string(scene->subdivs.size()) +
                 ".obj";
  return subdiv;
}
sceneio_texture* add_texture(sceneio_model* scene) {
  auto texture  = scene->textures.emplace_back(new sceneio_texture{});
  texture->name = "textures/texture" + std::to_string(scene->textures.size()) +
                  ".png";
  return texture;
}
sceneio_object* add_object(sceneio_model* scene) {
  auto object  = scene->objects.emplace_back(new sceneio_object{});
  object->name = "objects/object" + std::to_string(scene->objects.size()) +
                 ".json";
  return object;
}
sceneio_instance* add_instance(sceneio_model* scene) {
  auto instance  = scene->instances.emplace_back(new sceneio_instance{});
  instance->name = "instances/instance" +
                   std::to_string(scene->instances.size()) + ".ply";
  return instance;
}
sceneio_material* add_material(sceneio_model* scene) {
  auto material  = scene->materials.emplace_back(new sceneio_material{});
  material->name = "materials/material" +
                   std::to_string(scene->materials.size()) + ".json";
  return material;
}
sceneio_object* add_complete_object(
    sceneio_model* scene, const string& basename) {
  auto object      = add_object(scene);
  object->shape    = add_shape(scene);
  object->material = add_material(scene);
  if (!basename.empty()) {
    object->name           = "objects/" + basename + ".json";
    object->shape->name    = "materials/" + basename + ".json";
    object->material->name = "shapes/" + basename + ".ply";
  }
  return object;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const sceneio_model* scene) {
  auto shape_bbox = unordered_map<sceneio_shape*, bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto shape : scene->shapes) {
    auto sbvh = invalidb3f;
    for (auto p : shape->positions) sbvh = merge(sbvh, p);
    shape_bbox[shape] = sbvh;
  }
  for (auto object : scene->objects) {
    if (object->instance) {
      for (auto& frame : object->instance->frames) {
        auto sbvh = shape_bbox[object->shape];
        bbox      = merge(bbox, transform_bbox(frame * object->frame, sbvh));
      }
    } else {
      auto sbvh = shape_bbox[object->shape];
      bbox      = merge(bbox, transform_bbox(object->frame, sbvh));
    }
  }
  return bbox;
}

// Add missing cameras.
void add_cameras(sceneio_model* scene) {
  if (!scene->cameras.empty()) return;
  auto camera          = add_camera(scene);
  camera->name         = "cameras/default.json";
  camera->orthographic = false;
  camera->film         = 0.036;
  camera->aspect       = (float)16 / (float)9;
  camera->aperture     = 0;
  camera->lens         = 0.050;
  auto bbox            = compute_bounds(scene);
  auto center          = (bbox.max + bbox.min) / 2;
  auto bbox_radius     = length(bbox.max - bbox.min) / 2;
  auto camera_dir      = vec3f{0, 0, 1};
  auto camera_dist     = bbox_radius * camera->lens /
                     (camera->film / camera->aspect);
  camera_dist *= 2.0f;  // correction for tracer camera implementation
  auto from     = camera_dir * camera_dist + center;
  auto to       = center;
  auto up       = vec3f{0, 1, 0};
  camera->frame = lookat_frame(from, to, up);
  camera->focus = length(from - to);
}

// Add missing radius.
void add_radius(sceneio_model* scene, float radius = 0.001f) {
  for (auto shape : scene->shapes) {
    if (shape->points.empty() && shape->lines.empty()) continue;
    if (!shape->radius.empty()) continue;
    shape->radius.assign(shape->positions.size(), radius);
  }
}

// Add missing materials.
void add_materials(sceneio_model* scene) {
  auto default_material = (sceneio_material*)nullptr;
  for (auto& object : scene->objects) {
    if (object->material) continue;
    if (!default_material) {
      default_material        = add_material(scene);
      default_material->color = {0.8, 0.8, 0.8};
    }
    object->material = default_material;
  }
}

// Add a sky environment
void add_sky(sceneio_model* scene, float sun_angle) {
  auto texture              = add_texture(scene);
  texture->name             = "environments/sky.hdr";
  texture->hdr              = make_sunsky({1024, 512}, sun_angle);
  auto environment          = add_environment(scene);
  environment->name         = "environments/sky.yaml";
  environment->emission     = {1, 1, 1};
  environment->emission_tex = texture;
}

// Reduce memory usage
void trim_memory(sceneio_model* scene) {
  for (auto shape : scene->shapes) {
    shape->points.shrink_to_fit();
    shape->lines.shrink_to_fit();
    shape->triangles.shrink_to_fit();
    shape->quads.shrink_to_fit();
    shape->positions.shrink_to_fit();
    shape->normals.shrink_to_fit();
    shape->texcoords.shrink_to_fit();
    shape->colors.shrink_to_fit();
    shape->radius.shrink_to_fit();
    shape->tangents.shrink_to_fit();
  }
  for (auto subdiv : scene->subdivs) {
    subdiv->quadspos.shrink_to_fit();
    subdiv->quadsnorm.shrink_to_fit();
    subdiv->quadstexcoord.shrink_to_fit();
    subdiv->positions.shrink_to_fit();
    subdiv->normals.shrink_to_fit();
    subdiv->texcoords.shrink_to_fit();
  }
  for (auto texture : scene->textures) {
    texture->ldr.shrink_to_fit();
    texture->hdr.shrink_to_fit();
  }
  scene->cameras.shrink_to_fit();
  scene->shapes.shrink_to_fit();
  scene->textures.shrink_to_fit();
  scene->environments.shrink_to_fit();
}

// Apply subdivision and displacement rules.
unique_ptr<sceneio_subdiv> subdivide_subdiv(
    sceneio_subdiv* shape, int subdivisions, bool smooth) {
  using std::ignore;
  auto tesselated = make_unique<sceneio_subdiv>(*shape);
  if (!subdivisions) return tesselated;
  std::tie(tesselated->quadstexcoord, tesselated->texcoords) =
      subdivide_catmullclark(
          tesselated->quadstexcoord, tesselated->texcoords, subdivisions, true);
  std::tie(tesselated->quadsnorm, tesselated->normals) = subdivide_catmullclark(
      tesselated->quadsnorm, tesselated->normals, subdivisions, true);
  std::tie(tesselated->quadspos, tesselated->positions) =
      subdivide_catmullclark(
          tesselated->quadspos, tesselated->positions, subdivisions);
  if (smooth) {
    tesselated->normals = compute_normals(
        tesselated->quadspos, tesselated->positions);
    tesselated->quadsnorm = tesselated->quadspos;
  } else {
    tesselated->normals   = {};
    tesselated->quadsnorm = {};
  }
  return tesselated;
}
// Apply displacement to a shape
unique_ptr<sceneio_subdiv> displace_subdiv(sceneio_subdiv* subdiv,
    float displacement, sceneio_texture* displacement_tex, bool smooth) {
  // Evaluate a texture
  auto eval_texture = [](sceneio_texture* texture,
                          const vec2f&    texcoord) -> vec4f {
    if (!texture->hdr.empty()) {
      return eval_image(texture->hdr, texcoord, false, false);
    } else if (!texture->ldr.empty()) {
      return eval_image(texture->ldr, texcoord, true, false, false);
    } else {
      return {1, 1, 1, 1};
    }
  };

  auto displaced = make_unique<sceneio_subdiv>(*subdiv);

  if (!displacement || !displacement_tex) return displaced;
  if (subdiv->texcoords.empty())
    throw std::runtime_error("missing texture coordinates");

  // facevarying case
  auto offset = vector<float>(subdiv->positions.size(), 0);
  auto count  = vector<int>(subdiv->positions.size(), 0);
  for (auto fid = 0; fid < subdiv->quadspos.size(); fid++) {
    auto qpos = subdiv->quadspos[fid];
    auto qtxt = subdiv->quadstexcoord[fid];
    for (auto i = 0; i < 4; i++) {
      auto disp = mean(
          xyz(eval_texture(displacement_tex, subdiv->texcoords[qtxt[i]])));
      if (!displacement_tex->ldr.empty()) disp -= 0.5f;
      offset[qpos[i]] += displacement * disp;
      count[qpos[i]] += 1;
    }
  }
  auto normals = compute_normals(subdiv->quadspos, subdiv->positions);
  for (auto vid = 0; vid < subdiv->positions.size(); vid++) {
    displaced->positions[vid] += normals[vid] * offset[vid] / count[vid];
  }
  if (smooth || !subdiv->normals.empty()) {
    displaced->quadsnorm = subdiv->quadspos;
    displaced->normals   = compute_normals(
        displaced->quadspos, displaced->positions);
  }

  return displaced;
}

void tesselate_subdiv(sceneio_model* scene, sceneio_subdiv* subdiv) {
  auto material = (sceneio_material*)nullptr;
  auto shape    = (sceneio_shape*)nullptr;
  for (auto object : scene->objects) {
    if (object->subdiv == subdiv) {
      material = object->material;
      shape    = object->shape;
      break;
    }
  }
  auto tesselated = subdivide_subdiv(
      subdiv, material->subdivisions, material->smooth);
  auto displaced = displace_subdiv(tesselated.get(), material->displacement,
      material->displacement_tex, material->smooth);
  std::tie(shape->quads, shape->positions, shape->normals, shape->texcoords) =
      split_facevarying(displaced->quadspos, displaced->quadsnorm,
          displaced->quadstexcoord, displaced->positions, displaced->normals,
          displaced->texcoords);
  shape->points    = {};
  shape->lines     = {};
  shape->triangles = {};
  shape->colors    = {};
  shape->radius    = {};
}

void tesselate_subdivs(sceneio_model* scene, sceneio_progress progress_cb) {
  if (scene->subdivs.empty()) return;

  // handle progress
  auto progress = vec2i{0, (int)scene->subdivs.size()};

  // tesselate subdivs
  for (auto subdiv : scene->subdivs) {
    if (progress_cb) progress_cb("tesseleate subdiv", progress.x++, progress.y);
    tesselate_subdiv(scene, subdiv);
  }

  // done
  if (progress_cb) progress_cb("tesseleate subdiv", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
[[nodiscard]] static bool load_json_scene(const string& filename,
    sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);
[[nodiscard]] static bool save_json_scene(const string& filename,
    const sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);

// Load/save a scene from/to OBJ.
[[nodiscard]] static bool load_obj_scene(const string& filename,
    sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);
[[nodiscard]] static bool save_obj_scene(const string& filename,
    const sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
[[nodiscard]] static bool load_ply_scene(const string& filename,
    sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);
[[nodiscard]] static bool save_ply_scene(const string& filename,
    const sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);

// Load/save a scene from/to glTF.
[[nodiscard]] static bool load_gltf_scene(const string& filename,
    sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
[[nodiscard]] static bool load_pbrt_scene(const string& filename,
    sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);
[[nodiscard]] static bool save_pbrt_scene(const string& filename,
    const sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel);

// Load a scene
unique_ptr<sceneio_model> load_scene(const string& filename, string& error,
    sceneio_progress progress_cb, bool noparallel) {
  auto scene = make_unique<sceneio_model>();
  if (!load_scene(filename, scene.get(), error, progress_cb, noparallel))
    return nullptr;
  return scene;
}

// Load a scene
[[nodiscard]] bool load_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto ext = fs::path(filename).extension();
  if (ext == ".json" || ext == ".JSON") {
    return load_json_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    throw std::runtime_error{filename + ": unknown format"};
  }
}

// Save a scene
[[nodiscard]] bool save_scene(const string& filename,
    const sceneio_model* scene, string& error, sceneio_progress progress_cb,
    bool noparallel) {
  auto ext = fs::path(filename).extension();
  if (ext == ".json" || ext == ".JSON") {
    return save_json_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    throw std::runtime_error{filename + ": unknown format"};
  }
}

// create a name from another
static string make_safe_name(const string& prefix, const string& name,
    const string& ext = ".json", size_t count = 0) {
  if (name.empty()) return prefix + "s/" + prefix + std::to_string(count) + ext;
  auto lname = name;
  for (auto& c : lname) {
    c = tolower((int)c);
    if (c == ' ') c = '-';
  }
  if (count) lname += "[" + std::to_string(count) + "]";
  return prefix + "s/" + lname + ext;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto {

// Get extension (not including '.').
static string get_extension(const string& filename) {
  auto pos = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// load instances
[[nodiscard]] static bool load_instances(
    const string& filename, vector<frame3f>& frames, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, &ply, error)) return false;
    frames = get_values(&ply, "frame",
        array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz",
            "ox", "oy", "oz"});
    return true;
  } else {
    return format_error();
  }
}

// save instances
[[nodiscard]] static bool save_instances(const string& filename,
    const vector<frame3f>& frames, string& error, bool ascii = false) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(&ply, frames, "frame",
        array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz",
            "ox", "oy", "oz"});
    if (!save_ply(filename, &ply, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Load a text file
[[nodiscard]] inline bool load_text(
    const string& filename, string& str, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  str.resize(length);
  if (fread(str.data(), 1, length, fs) != length) return read_error();
  return true;
}

// Save a text file
[[nodiscard]] inline bool save_text(
    const string& filename, const string& str, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  if (fprintf(fs, "%s", str.c_str()) < 0) return write_error();
  return true;
}

// Load a binary file
inline string load_text(const string& filename, string& error) {
  auto text = string{};
  if (!load_text(filename, text, error)) return {};
  return text;
}

// Load a binary file
[[nodiscard]] inline bool load_binary(
    const string& filename, vector<byte>& data, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
  auto fs = fopen(filename.c_str(), "rb");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  fseek(fs, 0, SEEK_END);
  auto length = ftell(fs);
  fseek(fs, 0, SEEK_SET);
  data.resize(length);
  if (fread(data.data(), 1, length, fs) != length) return read_error();
  return true;
}

// Save a binary file
[[nodiscard]] inline bool save_binary(
    const string& filename, const vector<byte>& data, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  auto fs = fopen(filename.c_str(), "wb");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};
  if (fwrite(data.data(), 1, data.size(), fs) != data.size())
    return write_error();
  return true;
}

// Load a binary file
inline vector<byte> load_binary(const string& filename, string& error) {
  auto data = vector<byte>{};
  if (!load_binary(filename, data, error)) return {};
  return data;
}

using json = nlohmann::json;

// support for json conversions
inline void to_json(json& j, const vec3f& value) {
  nlohmann::to_json(j, (const array<float, 3>&)value);
}
inline void to_json(json& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}

inline void from_json(const json& j, vec3f& value) {
  nlohmann::from_json(j, (array<float, 3>&)value);
}
inline void from_json(const json& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const json& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}

// load/save json
[[nodiscard]] inline bool load_json(
    const string& filename, json& js, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    js = json::parse(text);
    return true;
  } catch (std::exception& e) {
    return parse_error();
  }
}

[[nodiscard]] inline bool save_json(
    const string& filename, const json& js, string& error) {
  return save_text(filename, js.dump(2), error);
}

inline json load_json(const string& filename, string& error) {
  auto js = json{};
  if (!load_json(filename, js, error)) return {};
  return js;
}

// Save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto material_error = [filename, &error](const string& name) {
    error = filename + ": missing material " + name;
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // open file
  auto js = json{};
  if (!load_json(filename, js, error)) return false;

  // gets a json value
  auto get_value = [](const json& ejs, const string& name,
                       auto& value) -> bool {
    if (!ejs.contains(name)) return true;
    try {
      ejs.at(name).get_to(value);
      return true;
    } catch (...) {
      return false;
    }
  };

  // parse yaml reference
  auto get_ref = [&material_error, &get_value](const json& ejs,
                     const string& name, auto& value,
                     const auto& refs) -> bool {
    if (!ejs.contains(name)) return true;
    auto ref = ""s;
    if (!get_value(ejs, name, ref)) return false;
    if (ref == "") {
      value = nullptr;
    } else {
      if (refs.find(ref) == refs.end()) return material_error(ref);
      value = refs.at(ref);
    }
    return true;
  };

  // parse json reference
  auto texture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_texture = [scene, &texture_map, &get_value](const json& ejs,
                         const string& name, sceneio_texture*& value,
                         const string& dirname = "textures/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) {
      value = it->second;
      return true;
    }
    auto texture      = add_texture(scene);
    texture->name     = path;
    texture_map[path] = texture;
    value             = texture;
    return true;
  };

  // parse json reference
  auto shape_map = unordered_map<string, sceneio_shape*>{{"", nullptr}};
  auto get_shape = [scene, &shape_map, &get_value](const json& ejs,
                       const string& name, sceneio_shape*& value,
                       const string& dirname = "shapes/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = shape_map.find(path);
    if (it != shape_map.end()) {
      value = it->second;
      return it->second;
    }
    auto shape      = add_shape(scene);
    shape->name     = path;
    shape_map[path] = shape;
    value           = shape;
    return true;
  };

  // parse json reference
  auto subdiv_map = unordered_map<string, sceneio_subdiv*>{{"", nullptr}};
  auto get_subdiv = [scene, &subdiv_map, &get_value](const json& ejs,
                        const string& name, sceneio_subdiv*& value,
                        const string& dirname = "subdivs/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = subdiv_map.find(path);
    if (it != subdiv_map.end()) {
      value = it->second;
      return true;
    }
    auto subdiv      = add_subdiv(scene);
    subdiv->name     = path;
    subdiv_map[path] = subdiv;
    value            = subdiv;
    return true;
  };

  // load json instance
  auto instance_map = unordered_map<string, sceneio_instance*>{{"", nullptr}};
  auto get_instance = [scene, &instance_map, &get_value](const json& ejs,
                          const string& name, sceneio_instance*& value,
                          const string& dirname = "instances/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = instance_map.find(path);
    if (it != instance_map.end()) {
      value = it->second;
      return true;
    }
    auto instance      = add_instance(scene);
    instance->name     = path;
    instance_map[path] = instance;
    value              = instance;
    return true;
  };

  // material map
  auto material_map = unordered_map<string, sceneio_material*>{{"", nullptr}};

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // check for conversion errors
  // cameras
  if (js.contains("cameras")) {
    for (auto& ejs : js.at("cameras")) {
      auto camera = add_camera(scene);
      if (!get_value(ejs, "name", camera->name)) return false;
      if (!get_value(ejs, "frame", camera->frame)) return false;
      if (!get_value(ejs, "orthographic", camera->orthographic)) return false;
      if (!get_value(ejs, "lens", camera->lens)) return false;
      if (!get_value(ejs, "aspect", camera->aspect)) return false;
      if (!get_value(ejs, "film", camera->film)) return false;
      if (!get_value(ejs, "focus", camera->focus)) return false;
      if (!get_value(ejs, "aperture", camera->aperture)) return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        camera->frame = lookat_frame(lookat.x, lookat.y, lookat.z);
        camera->focus = length(lookat.x - lookat.y);
      }
    }
  }
  if (js.contains("environments")) {
    for (auto& ejs : js.at("environments")) {
      auto environment = add_environment(scene);
      if (!get_value(ejs, "name", environment->name)) return false;
      if (!get_value(ejs, "frame", environment->frame)) return false;
      if (!get_value(ejs, "emission", environment->emission)) return false;
      if (!get_texture(
              ejs, "emission_tex", environment->emission_tex, "environments/"))
        return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        environment->frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
    }
  }
  if (js.contains("materials")) {
    for (auto& ejs : js.at("materials")) {
      auto material = add_material(scene);
      if (!get_value(ejs, "name", material->name)) return false;
      if (!get_value(ejs, "emission", material->emission)) return false;
      if (!get_value(ejs, "color", material->color)) return false;
      if (!get_value(ejs, "metallic", material->metallic)) return false;
      if (!get_value(ejs, "specular", material->specular)) return false;
      if (!get_value(ejs, "roughness", material->roughness)) return false;
      if (!get_value(ejs, "coat", material->coat)) return false;
      if (!get_value(ejs, "transmission", material->transmission)) return false;
      if (!get_value(ejs, "thin", material->thin)) return false;
      if (!get_value(ejs, "ior", material->ior)) return false;
      if (!get_value(ejs, "trdepth", material->trdepth)) return false;
      if (!get_value(ejs, "scattering", material->scattering)) return false;
      if (!get_value(ejs, "scanisotropy", material->scanisotropy)) return false;
      if (!get_value(ejs, "opacity", material->opacity)) return false;
      if (!get_value(ejs, "coat", material->coat)) return false;
      if (!get_value(ejs, "displacement", material->displacement)) return false;
      if (!get_texture(ejs, "emission_tex", material->emission_tex))
        return false;
      if (!get_texture(ejs, "color_tex", material->color_tex)) return false;
      if (!get_texture(ejs, "metallic_tex", material->metallic_tex))
        return false;
      if (!get_texture(ejs, "specular_tex", material->specular_tex))
        return false;
      if (!get_texture(ejs, "transmission_tex", material->transmission_tex))
        return false;
      if (!get_texture(ejs, "roughness_tex", material->roughness_tex))
        return false;
      if (!get_texture(ejs, "scattering_tex", material->scattering_tex))
        return false;
      if (!get_texture(ejs, "opacity_tex", material->opacity_tex)) return false;
      if (!get_texture(ejs, "normal_tex", material->normal_tex)) return false;
      if (!get_texture(ejs, "displacement_tex", material->displacement_tex))
        return false;
      if (!get_value(ejs, "subdivisions", material->subdivisions))
        return false;  // hack fir subd
      if (!get_value(ejs, "smooth", material->smooth))
        return false;  // hack for subd
      if (!get_value(ejs, "gltf_textures", material->gltf_textures))
        return false;
      material_map[material->name] = material;
    }
  }
  if (js.contains("objects")) {
    for (auto& ejs : js.at("objects")) {
      auto object = add_object(scene);
      if (!get_value(ejs, "name", object->name)) return false;
      if (!get_value(ejs, "frame", object->frame)) return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        object->frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
      if (!get_ref(ejs, "material", object->material, material_map))
        return false;
      if (!get_shape(ejs, "shape", object->shape)) return false;
      if (!get_subdiv(ejs, "subdiv", object->subdiv)) return false;
      if (!get_instance(ejs, "instance", object->instance)) return false;
    }
  }

  // handle progress
  progress.y += scene->shapes.size();
  progress.y += scene->subdivs.size();
  progress.y += scene->textures.size();
  progress.y += scene->instances.size();

  // load shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
    if (!load_shape(fs::path(filename).parent_path() / shape->name,
            shape->points, shape->lines, shape->triangles, shape->quads,
            shape->positions, shape->normals, shape->texcoords, shape->colors,
            shape->radius, error))
      return dependent_error();
  }
  // load subdivs
  for (auto subdiv : scene->subdivs) {
    if (progress_cb) progress_cb("load subdiv", progress.x++, progress.y);
    if (!load_fvshape(fs::path(filename).parent_path() / subdiv->name,
            subdiv->quadspos, subdiv->quadsnorm, subdiv->quadstexcoord,
            subdiv->positions, subdiv->normals, subdiv->texcoords, error))
      return dependent_error();
  }
  // load textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (is_hdr_filename(texture->name)) {
      if (!load_image(fs::path(filename).parent_path() / texture->name,
              texture->hdr, error))
        return dependent_error();
    } else {
      if (!load_imageb(fs::path(filename).parent_path() / texture->name,
              texture->ldr, error))
        return dependent_error();
    }
  }
  // load instances
  for (auto instance : scene->instances) {
    if (progress_cb) progress_cb("load instance", progress.x++, progress.y);
    if (!load_instances(fs::path(filename).parent_path() / instance->name,
            instance->frames, error))
      return dependent_error();
  }

  // fix scene
  scene->name = filename;
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);
  trim_memory(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // helper
  auto add_val = [](json& ejs, const string& name, const auto& value) {
    ejs[name] = value;
  };
  auto add_opt = [](json& ejs, const string& name, const auto& value,
                     const auto& def) {
    if (value == def) return;
    ejs[name] = value;
  };
  auto add_tex = [](json& ejs, const string& name, sceneio_texture* texture) {
    if (!texture) return;
    ejs[name] = texture->name;
  };
  auto add_ref = [](json& ejs, const string& name, auto ref) {
    if (!ref) return;
    ejs[name] = ref->name;
  };

  // handle progress
  auto progress = vec2i{
      0, 2 + (int)scene->shapes.size() + (int)scene->subdivs.size() +
             (int)scene->textures.size() + (int)scene->instances.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save yaml file
  auto js     = json::object();
  js["asset"] = json::object();

  auto def_cam = sceneio_camera{};
  if (!scene->cameras.empty()) js["cameras"] = json::array();
  for (auto& camera : scene->cameras) {
    auto& ejs = js["cameras"].emplace_back();
    add_val(ejs, "name", camera->name);
    add_opt(ejs, "frame", camera->frame, def_cam.frame);
    add_opt(ejs, "ortho", camera->orthographic, def_cam.orthographic);
    add_opt(ejs, "lens", camera->lens, def_cam.lens);
    add_opt(ejs, "aspect", camera->aspect, def_cam.aspect);
    add_opt(ejs, "film", camera->film, def_cam.film);
    add_opt(ejs, "focus", camera->focus, def_cam.focus);
    add_opt(ejs, "aperture", camera->aperture, def_cam.aperture);
  }

  auto def_env = sceneio_environment{};
  if (!scene->environments.empty()) js["environments"] = json::array();
  for (auto environment : scene->environments) {
    auto& ejs = js["environments"].emplace_back();
    add_val(ejs, "name", environment->name);
    add_opt(ejs, "frame", environment->frame, def_env.frame);
    add_opt(ejs, "emission", environment->emission, def_env.emission);
    add_tex(ejs, "emission_tex", environment->emission_tex);
  }

  auto def_material = sceneio_material{};
  if (!scene->materials.empty()) js["materials"] = json::array();
  for (auto material : scene->materials) {
    auto& ejs = js["materials"].emplace_back();
    add_val(ejs, "name", material->name);
    add_opt(ejs, "emission", material->emission, def_material.emission);
    add_opt(ejs, "color", material->color, def_material.color);
    add_opt(ejs, "specular", material->specular, def_material.specular);
    add_opt(ejs, "metallic", material->metallic, def_material.metallic);
    add_opt(ejs, "coat", material->coat, def_material.coat);
    add_opt(ejs, "roughness", material->roughness, def_material.roughness);
    add_opt(ejs, "ior", material->ior, def_material.ior);
    add_opt(
        ejs, "transmission", material->transmission, def_material.transmission);
    add_opt(ejs, "trdepth", material->trdepth, def_material.trdepth);
    add_opt(ejs, "scattering", material->scattering, def_material.scattering);
    add_opt(
        ejs, "scanisotropy", material->scanisotropy, def_material.scanisotropy);
    add_opt(ejs, "opacity", material->opacity, def_material.opacity);
    add_opt(ejs, "displacement", material->opacity, def_material.displacement);
    add_opt(ejs, "thin", material->thin, def_material.thin);
    add_tex(ejs, "emission_tex", material->emission_tex);
    add_tex(ejs, "color_tex", material->color_tex);
    add_tex(ejs, "metallic_tex", material->metallic_tex);
    add_tex(ejs, "specular_tex", material->specular_tex);
    add_tex(ejs, "roughness_tex", material->roughness_tex);
    add_tex(ejs, "transmission_tex", material->transmission_tex);
    add_tex(ejs, "scattering_tex", material->scattering_tex);
    add_tex(ejs, "coat_tex", material->coat_tex);
    add_tex(ejs, "opacity_tex", material->opacity_tex);
    add_tex(ejs, "normal_tex", material->normal_tex);
    add_tex(ejs, "displacement_tex", material->displacement_tex);
    add_opt(ejs, "subdivisions", material->subdivisions,
        def_material.subdivisions);  // hack fir subd
    add_opt(
        ejs, "smooth", material->smooth, def_material.smooth);  // hack for subd
    add_opt(ejs, "gltf_textures", material->gltf_textures,
        def_material.gltf_textures);
  }

  auto def_object = sceneio_object{};
  auto def_subdiv = sceneio_subdiv{};
  if (!scene->objects.empty()) js["objects"] = json::array();
  for (auto object : scene->objects) {
    auto& ejs = js["objects"].emplace_back();
    add_val(ejs, "name", object->name);
    add_opt(ejs, "frame", object->frame, def_object.frame);
    add_ref(ejs, "shape", object->shape);
    add_ref(ejs, "subdiv", object->subdiv);
    add_ref(ejs, "material", object->material);
    add_ref(ejs, "instance", object->instance);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json
  if (!save_json(filename, js, error)) return false;

  // save shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    if (!shape->positions.empty()) {
      if (!save_shape(fs::path(filename).parent_path() / shape->name,
              shape->points, shape->lines, shape->triangles, shape->quads,
              shape->positions, shape->normals, shape->texcoords, shape->colors,
              shape->radius, error))
        return dependent_error();
    }
  }

  // save subdivs
  for (auto subdiv : scene->subdivs) {
    if (progress_cb) progress_cb("save subdiv", progress.x++, progress.y);
    if (!subdiv->positions.empty()) {
      if (!save_fvshape(fs::path(filename).parent_path() / subdiv->name,
              subdiv->quadspos, subdiv->quadsnorm, subdiv->quadstexcoord,
              subdiv->positions, subdiv->normals, subdiv->texcoords, error))
        return dependent_error();
    }
  }

  // save instances
  for (auto instance : scene->instances) {
    if (progress_cb) progress_cb("save instance", progress.x++, progress.y);
    if (!instance->frames.empty()) {
      if (!save_instances(fs::path(filename).parent_path() / instance->name,
              instance->frames, error))
        return dependent_error();
    }
  }

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    if (!texture->ldr.empty() || !texture->hdr.empty()) {
      if (!texture->hdr.empty()) {
        if (!save_image(fs::path(filename).parent_path() / texture->name,
                texture->hdr, error))
          return dependent_error();
      } else {
        if (!save_imageb(fs::path(filename).parent_path() / texture->name,
                texture->ldr, error))
          return dependent_error();
      }
    }
  }

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Loads an OBJ
static bool load_obj_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load obj
  auto obj_guard = make_unique<obj_model>();
  auto obj       = obj_guard.get();
  if (!load_obj(filename, obj, error, false, true, false)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto ocam : obj->cameras) {
    auto camera = add_camera(scene);
    // camera->name         = make_safe_name("camera", ocam->name);
    camera->frame        = ocam->frame;
    camera->orthographic = ocam->ortho;
    camera->film         = max(ocam->width, ocam->height);
    camera->aspect       = ocam->width / ocam->height;
    camera->focus        = ocam->focus;
    camera->lens         = ocam->lens;
    camera->aperture     = ocam->aperture;
  }

  // helper to create texture maps
  auto texture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_texture = [&texture_map, scene](
                         const obj_texture_info& info) -> sceneio_texture* {
    auto path = info.path;
    if (path == "") return nullptr;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    auto texture = add_texture(scene);
    if (is_hdr_filename(path))
      texture->name = fs::path(texture->name).replace_extension(".hdr");
    // texture->name = make_safe_name(
    //     "texture", get_basename(path), is_hdr_filename(path) ? ".hdr" :
    //     ".png");
    texture_map[path] = texture;
    return texture;
  };

  // handler for materials
  auto material_map = unordered_map<obj_material*, sceneio_material*>{};
  for (auto omat : obj->materials) {
    auto material = add_material(scene);
    // material->name             = make_safe_name("material", omat->name);
    material->emission         = omat->pbr_emission;
    material->color            = omat->pbr_base;
    material->specular         = omat->pbr_specular;
    material->roughness        = omat->pbr_roughness;
    material->ior              = omat->pbr_ior;
    material->metallic         = omat->pbr_metallic;
    material->coat             = omat->pbr_coat;
    material->transmission     = omat->pbr_transmission;
    material->scattering       = omat->pbr_volscattering;
    material->scanisotropy     = omat->pbr_volanisotropy;
    material->trdepth          = omat->pbr_volscale;
    material->opacity          = omat->pbr_opacity;
    material->thin             = true;
    material->emission_tex     = get_texture(omat->pbr_emission_tex);
    material->color_tex        = get_texture(omat->pbr_base_tex);
    material->specular_tex     = get_texture(omat->pbr_specular_tex);
    material->metallic_tex     = get_texture(omat->pbr_metallic_tex);
    material->roughness_tex    = get_texture(omat->pbr_roughness_tex);
    material->transmission_tex = get_texture(omat->pbr_transmission_tex);
    material->coat_tex         = get_texture(omat->pbr_coat_tex);
    material->opacity_tex      = get_texture(omat->pbr_opacity_tex);
    material->normal_tex       = get_texture(omat->normal_tex);
    material_map[omat]         = material;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto oshape : obj->shapes) {
    auto materials = get_materials(obj, oshape);
    if (materials.empty()) materials.push_back(nullptr);
    for (auto material_idx = 0; material_idx < materials.size();
         material_idx++) {
      auto object      = add_object(scene);
      object->shape    = add_shape(scene);
      object->material = material_map.at(materials[material_idx]);
      // if (!oshape->name.empty()) {
      //   object->name        = make_safe_name("object", oshape->name, ".json",
      //       materials.size() > 1 ? material_idx + 1 : 0);
      //   object->shape->name = make_safe_name("shape", oshape->name, ".ply",
      //       materials.size() > 1 ? material_idx + 1 : 0);
      // }
      auto has_quads_ = has_quads(oshape);
      if (!oshape->faces.empty() && !has_quads_) {
        get_triangles(obj, oshape, material_idx, object->shape->triangles,
            object->shape->positions, object->shape->normals,
            object->shape->texcoords, true);
      } else if (!oshape->faces.empty() && has_quads_) {
        get_quads(obj, oshape, material_idx, object->shape->quads,
            object->shape->positions, object->shape->normals,
            object->shape->texcoords, true);
      } else if (!oshape->lines.empty()) {
        get_lines(obj, oshape, material_idx, object->shape->lines,
            object->shape->positions, object->shape->normals,
            object->shape->texcoords, true);
      } else if (!oshape->points.empty()) {
        get_points(obj, oshape, material_idx, object->shape->points,
            object->shape->positions, object->shape->normals,
            object->shape->texcoords, true);
      } else {
        return shape_error();
      }
      if (!oshape->instances.empty()) {
        object->instance = add_instance(scene);
        // if (!oshape->name.empty()) {
        //   object->instance->name = make_safe_name("object", oshape->name,
        //       ".json", materials.size() > 1 ? material_idx + 1 : 0);
        // }
        object->instance->frames = oshape->instances;
      }
    }
  }

  // convert environments
  for (auto oenvironment : obj->environments) {
    auto environment = add_environment(scene);
    // environment->name     = make_safe_name("environment", oenvironment.name);
    environment->frame        = oenvironment->frame;
    environment->emission     = oenvironment->emission;
    environment->emission_tex = get_texture(oenvironment->emission_tex);
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // load textures
  texture_map.erase("");
  for (auto [path, texture] : texture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (is_hdr_filename(path)) {
      if (!load_image(
              fs::path(filename).parent_path() / path, texture->hdr, error))
        return dependent_error();
    } else {
      if (!load_imageb(
              fs::path(filename).parent_path() / path, texture->ldr, error))
        return dependent_error();
    }
  }

  // fix scene
  scene->name = filename;
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_obj_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2 + (int)scene->textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  auto obj_guard = make_obj();
  auto obj       = obj_guard.get();

  // convert cameras
  for (auto camera : scene->cameras) {
    auto ocamera      = add_camera(obj);
    ocamera->name     = fs::path(camera->name).stem();
    ocamera->frame    = camera->frame;
    ocamera->ortho    = camera->orthographic;
    ocamera->width    = camera->film;
    ocamera->height   = camera->film / camera->aspect;
    ocamera->focus    = camera->focus;
    ocamera->lens     = camera->lens;
    ocamera->aperture = camera->aperture;
  }

  // textures
  auto get_texture = [](sceneio_texture* texture) {
    if (!texture) return obj_texture_info{};
    auto info = obj_texture_info{};
    info.path = texture->name;
    return info;
  };

  // convert materials and textures
  for (auto material : scene->materials) {
    auto omaterial                  = add_material(obj);
    omaterial->name                 = fs::path(material->name).stem();
    omaterial->illum                = 2;
    omaterial->as_pbr               = true;
    omaterial->pbr_emission         = material->emission;
    omaterial->pbr_base             = material->color;
    omaterial->pbr_specular         = material->specular;
    omaterial->pbr_roughness        = material->roughness;
    omaterial->pbr_metallic         = material->metallic;
    omaterial->pbr_coat             = material->coat;
    omaterial->pbr_transmission     = material->transmission;
    omaterial->pbr_opacity          = material->opacity;
    omaterial->pbr_emission_tex     = get_texture(material->emission_tex);
    omaterial->pbr_base_tex         = get_texture(material->color_tex);
    omaterial->pbr_specular_tex     = get_texture(material->specular_tex);
    omaterial->pbr_metallic_tex     = get_texture(material->metallic_tex);
    omaterial->pbr_roughness_tex    = get_texture(material->roughness_tex);
    omaterial->pbr_transmission_tex = get_texture(material->transmission_tex);
    omaterial->pbr_coat_tex         = get_texture(material->coat_tex);
    omaterial->pbr_opacity_tex      = get_texture(material->opacity_tex);
    omaterial->normal_tex           = get_texture(material->normal_tex);
  }

  // convert objects
  for (auto object : scene->objects) {
    auto shape     = object->shape;
    auto positions = shape->positions, normals = shape->normals;
    for (auto& p : positions) p = transform_point(object->frame, p);
    for (auto& n : normals) n = transform_normal(object->frame, n);
    if (!shape->triangles.empty()) {
      add_triangles(obj, shape->name, shape->triangles, positions, normals,
          shape->texcoords, {}, {},
          object->instance ? object->instance->frames : vector<frame3f>{},
          true);
    } else if (!shape->quads.empty()) {
      add_quads(obj, shape->name, shape->quads, positions, normals,
          shape->texcoords, {}, {},
          object->instance ? object->instance->frames : vector<frame3f>{},
          true);
    } else if (!shape->lines.empty()) {
      add_lines(obj, shape->name, shape->lines, positions, normals,
          shape->texcoords, {}, {},
          object->instance ? object->instance->frames : vector<frame3f>{},
          true);
    } else if (!shape->points.empty()) {
      add_points(obj, shape->name, shape->points, positions, normals,
          shape->texcoords, {}, {},
          object->instance ? object->instance->frames : vector<frame3f>{},
          true);
    } else {
      return shape_error();
    }
  }

  // convert environments
  for (auto environment : scene->environments) {
    auto oenvironment          = add_environment(obj);
    oenvironment->name         = fs::path(environment->name).stem();
    oenvironment->frame        = environment->frame;
    oenvironment->emission     = environment->emission;
    oenvironment->emission_tex = get_texture(environment->emission_tex);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save obj
  if (!save_obj(filename, obj, error)) return false;

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    if (texture->ldr.empty() && texture->hdr.empty()) continue;
    if (!texture->hdr.empty()) {
      if (!save_image(fs::path(filename).parent_path() / texture->name,
              texture->hdr, error))
        return dependent_error();
    } else {
      if (!save_imageb(fs::path(filename).parent_path() / texture->name,
              texture->ldr, error))
        return dependent_error();
    }
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

void print_obj_camera(sceneio_camera* camera) {
  printf("c %s %d %g %g %g %g %g %g %g %g %g %g%g %g %g %g %g %g %g\n",
      camera->name.c_str(), (int)camera->orthographic, camera->film,
      camera->film / camera->aspect, camera->lens, camera->focus,
      camera->aperture, camera->frame.x.x, camera->frame.x.y, camera->frame.x.z,
      camera->frame.y.x, camera->frame.y.y, camera->frame.y.z,
      camera->frame.z.x, camera->frame.z.y, camera->frame.z.z,
      camera->frame.o.x, camera->frame.o.y, camera->frame.o.z);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_ply_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load ply mesh
  auto shape = add_shape(scene);
  if (!load_shape(filename, shape->points, shape->lines, shape->triangles,
          shape->quads, shape->positions, shape->normals, shape->texcoords,
          shape->colors, shape->radius, error))
    return false;

  // fix scene
  scene->name = filename;
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_ply_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  if (scene->shapes.empty())
    throw std::runtime_error{filename + ": empty shape"};

  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save shape
  auto shape = scene->shapes.front();
  if (!save_shape(filename, shape->points, shape->lines, shape->triangles,
          shape->quads, shape->positions, shape->normals, shape->texcoords,
          shape->colors, shape->radius, error))
    return false;

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static bool load_gltf_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load gltf
  auto gltf_guard = make_unique<gltf_model>();
  auto gltf       = gltf_guard.get();
  if (!load_gltf(filename, gltf, error)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto gcamera : gltf->cameras) {
    for (auto frame : gcamera->frames) {
      auto camera    = add_camera(scene);
      camera->frame  = frame;
      camera->aspect = gcamera->aspect;
      camera->film   = 0.036;
      camera->lens   = gcamera->aspect >= 1
                         ? (2 * camera->aspect * tan(gcamera->yfov / 2))
                         : (2 * tan(gcamera->yfov / 2));
      camera->focus = 10;
    }
  }

  // convert textures
  auto texture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_texture = [&scene, &texture_map](
                         gltf_texture* gtexture) -> sceneio_texture* {
    if (!gtexture) return nullptr;
    auto path = gtexture->filename;
    if (path == "") return nullptr;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    auto texture      = add_texture(scene);
    texture->name     = make_safe_name("texture", fs::path(path).stem(),
        (!texture->ldr.empty() ? ".png" : ".hdr"));
    texture_map[path] = texture;
    return texture;
  };

  // convert materials
  auto material_map = unordered_map<gltf_material*, sceneio_material*>{
      {nullptr, nullptr}};
  for (auto gmaterial : gltf->materials) {
    auto material           = add_material(scene);
    material->emission      = gmaterial->emission;
    material->emission_tex  = get_texture(gmaterial->emission_tex);
    material->color         = gmaterial->color;
    material->opacity       = gmaterial->opacity;
    material->specular      = 1;
    material->color_tex     = get_texture(gmaterial->color_tex);
    material->metallic_tex  = get_texture(gmaterial->metallic_tex);
    material->normal_tex    = get_texture(gmaterial->normal_tex);
    material_map[gmaterial] = material;
  }

  // convert shapes
  auto shape_map = unordered_map<gltf_primitive*, sceneio_shape*>{
      {nullptr, nullptr}};
  for (auto gprim : gltf->primitives) {
    auto shape       = add_shape(scene);
    shape_map[gprim] = shape;
    shape->positions = gprim->positions;
    shape->normals   = gprim->normals;
    shape->texcoords = gprim->texcoords;
    shape->colors    = gprim->colors;
    shape->radius    = gprim->radius;
    shape->tangents  = gprim->tangents;
    shape->triangles = gprim->triangles;
    shape->lines     = gprim->lines;
    shape->points    = gprim->points;
  }

  // convert object
  for (auto gmesh : gltf->meshes) {
    for (auto gprim : gmesh->primitives) {
      auto object   = add_object(scene);
      object->frame = gmesh->frames.empty() ? identity3x4f
                                            : gmesh->frames.front();
      object->shape    = shape_map.at(gprim);
      object->material = material_map.at(gprim->material);
      object->frame    = identity3x4f;
      if (gmesh->frames.size() == 1) {
        object->frame = gmesh->frames.front();
      } else {
        object->instance         = add_instance(scene);
        object->instance->frames = gmesh->frames;
      }
    }
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // load texture
  texture_map.erase("");
  for (auto [path, texture] : texture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (is_hdr_filename(path)) {
      if (!load_image(
              fs::path(filename).parent_path() / path, texture->hdr, error))
        return dependent_error();
    } else {
      if (!load_imageb(
              fs::path(filename).parent_path() / path, texture->ldr, error))
        return dependent_error();
    }
  }

  // fix scene
  scene->name = filename;
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // fix cameras
  auto bbox = compute_bounds(scene);
  for (auto camera : scene->cameras) {
    auto center   = (bbox.min + bbox.max) / 2;
    auto distance = dot(-camera->frame.z, center - camera->frame.o);
    if (distance > 0) camera->focus = distance;
  }

  // load done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static bool load_pbrt_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load pbrt
  auto pbrt_guard = make_unique<pbrt_model>();
  auto pbrt       = pbrt_guard.get();
  if (!load_pbrt(filename, pbrt, error)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto pcamera : pbrt->cameras) {
    auto camera    = add_camera(scene);
    camera->frame  = pcamera->frame;
    camera->aspect = pcamera->aspect;
    camera->film   = 0.036;
    camera->lens   = pcamera->lens;
    camera->focus  = pcamera->focus;
  }

  // convert materials
  auto texture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_texture = [&scene, &texture_map](
                         const string& path) -> sceneio_texture* {
    if (path == "") return nullptr;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    auto texture      = add_texture(scene);
    texture->name     = make_safe_name("texture", fs::path(path).stem(),
        (!is_hdr_filename(path) ? ".png" : ".hdr"));
    texture_map[path] = texture;
    return texture;
  };
  auto alpha_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_alpha = [&scene, &alpha_map](
                       const string& path) -> sceneio_texture* {
    if (path == "") return nullptr;
    auto it = alpha_map.find(path);
    if (it != alpha_map.end()) return it->second;
    auto texture  = add_texture(scene);
    texture->name = make_safe_name(
        "texture", fs::path(path).stem(), "-alpha.png");
    alpha_map[path] = texture;
    return texture;
  };

  // convert material
  auto material_map = unordered_map<pbrt_material*, sceneio_material*>{};
  for (auto pmaterial : pbrt->materials) {
    auto material          = add_material(scene);
    material->emission     = pmaterial->emission;
    material->color        = pmaterial->color;
    material->metallic     = pmaterial->metallic;
    material->specular     = pmaterial->specular;
    material->transmission = pmaterial->transmission;
    material->ior          = pmaterial->ior;
    material->roughness    = pmaterial->roughness;
    material->opacity      = pmaterial->opacity;
    material->thin         = pmaterial->thin;
    material->color_tex    = get_texture(pmaterial->color_tex);
    material->opacity_tex  = get_texture(pmaterial->opacity_tex);
    if (!material->opacity_tex)
      material->opacity_tex = get_alpha(pmaterial->alpha_tex);
    material_map[pmaterial] = material;
  }

  // hack for pbrt empty material
  material_map[nullptr] = add_material(scene);

  // convert shapes
  for (auto pshape : pbrt->shapes) {
    auto object   = add_object(scene);
    object->shape = add_shape(scene);
    object->frame = pshape->frame;
    if (!pshape->instances.empty()) {
      object->instance         = add_instance(scene);
      object->instance->frames = pshape->instances;
    }
    object->shape->positions = pshape->positions;
    object->shape->normals   = pshape->normals;
    object->shape->texcoords = pshape->texcoords;
    object->shape->triangles = pshape->triangles;
    for (auto& uv : object->shape->texcoords) uv.y = 1 - uv.y;
    object->material = material_map.at(pshape->material);
  }

  // convert environments
  for (auto penvironment : pbrt->environments) {
    auto environment          = add_environment(scene);
    environment->frame        = penvironment->frame;
    environment->emission     = penvironment->emission;
    environment->emission_tex = get_texture(penvironment->emission_tex);
  }

  // lights
  for (auto plight : pbrt->lights) {
    auto object                = add_object(scene);
    object->shape              = add_shape(scene);
    object->frame              = plight->area_frame;
    object->shape->triangles   = plight->area_triangles;
    object->shape->positions   = plight->area_positions;
    object->shape->normals     = plight->area_normals;
    object->material           = add_material(scene);
    object->material->emission = plight->area_emission;
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // load texture
  texture_map.erase("");
  for (auto [path, texture] : texture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (is_hdr_filename(path)) {
      if (!load_image(
              fs::path(filename).parent_path() / path, texture->hdr, error))
        return dependent_error();
    } else {
      if (!load_imageb(
              fs::path(filename).parent_path() / path, texture->ldr, error))
        return dependent_error();
    }
  }

  // load alpha
  alpha_map.erase("");
  for (auto [path, texture] : alpha_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (is_hdr_filename(path)) {
      if (!load_image(
              fs::path(filename).parent_path() / path, texture->hdr, error))
        return dependent_error();
      for (auto& c : texture->hdr)
        c = (max(xyz(c)) < 0.01) ? vec4f{0, 0, 0, 0} : vec4f{1, 1, 1, 1};
    } else {
      if (!load_imageb(
              fs::path(filename).parent_path() / path, texture->ldr, error))
        return dependent_error();
      for (auto& c : texture->ldr) {
        c = (max(max(c.x, c.y), c.z) < 2) ? vec4b{0, 0, 0, 0}
                                          : vec4b{255, 255, 255, 255};
      }
    }
  }

  // fix scene
  scene->name = filename;
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt_guard = make_pbrt();
  auto pbrt       = pbrt_guard.get();

  // convert camera
  auto camera         = scene->cameras.front();
  auto pcamera        = add_camera(pbrt);
  pcamera->frame      = camera->frame;
  pcamera->lens       = camera->lens;
  pcamera->aspect     = camera->aspect;
  pcamera->resolution = {1280, (int)(1280 / pcamera->aspect)};

  // get texture name
  auto get_texture = [](const sceneio_texture* texture) {
    return texture ? texture->name : "";
  };

  // convert materials
  auto material_map = unordered_map<sceneio_material*, pbrt_material*>{};
  for (auto material : scene->materials) {
    auto pmaterial          = add_material(pbrt);
    pmaterial->name         = fs::path(material->name).stem();
    pmaterial->emission     = material->emission;
    pmaterial->color        = material->color;
    pmaterial->metallic     = material->metallic;
    pmaterial->specular     = material->specular;
    pmaterial->transmission = material->transmission;
    pmaterial->roughness    = material->roughness;
    pmaterial->ior          = material->ior;
    pmaterial->opacity      = material->opacity;
    pmaterial->color_tex    = get_texture(material->color_tex);
    pmaterial->opacity_tex  = get_texture(material->opacity_tex);
    material_map[material]  = pmaterial;
  }

  // convert instances
  for (auto object : scene->objects) {
    auto pshape       = add_shape(pbrt);
    pshape->filename_ = fs::path(object->shape->name).replace_extension(".ply");
    pshape->frame     = object->frame;
    pshape->frend     = object->frame;
    pshape->material  = material_map.at(object->material);
    if (object->instance) {
      pshape->instances = object->instance->frames;
      pshape->instances = object->instance->frames;
    }
  }

  // convert environments
  for (auto environment : scene->environments) {
    auto penvironment          = add_environment(pbrt);
    penvironment->emission     = environment->emission;
    penvironment->emission_tex = get_texture(environment->emission_tex);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  if (!save_pbrt(filename, pbrt, error)) return false;

  // handle progress
  progress.y += (int)scene->shapes.size() + (int)scene->textures.size();

  // save meshes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    if (shape->positions.empty()) continue;
    if (!save_shape((fs::path(filename).parent_path() / shape->name)
                        .replace_extension(".ply"),
            shape->points, shape->lines, shape->triangles, shape->quads,
            shape->positions, shape->normals, shape->texcoords, shape->colors,
            shape->radius, error))
      return dependent_error();
  }

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    if (texture->ldr.empty() && texture->hdr.empty()) continue;
    if (!texture->hdr.empty()) {
      if (!save_image(fs::path(filename).parent_path() / texture->name,
              texture->hdr, error))
        return dependent_error();
    } else {
      if (!save_imageb(fs::path(filename).parent_path() / texture->name,
              texture->ldr, error))
        return dependent_error();
    }
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox_scene(sceneio_model* scene) {
  scene->name                = "cornellbox";
  auto camera                = add_camera(scene);
  camera->frame              = frame3f{{0, 1, 3.9}};
  camera->lens               = 0.035;
  camera->aperture           = 0.0;
  camera->film               = 0.024;
  camera->aspect             = 1;
  auto floor                 = add_complete_object(scene, "floor");
  floor->shape->positions    = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor->shape->triangles    = {{0, 1, 2}, {2, 3, 0}};
  floor->material->color     = {0.725, 0.71, 0.68};
  auto ceiling               = add_complete_object(scene, "ceiling");
  ceiling->shape->positions  = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling->shape->triangles  = {{0, 1, 2}, {2, 3, 0}};
  ceiling->material->color   = {0.725, 0.71, 0.68};
  auto backwall              = add_complete_object(scene, "backwall");
  backwall->shape->positions = {
      {-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall->shape->triangles  = {{0, 1, 2}, {2, 3, 0}};
  backwall->material->color   = {0.725, 0.71, 0.68};
  auto rightwall              = add_complete_object(scene, "rightwall");
  rightwall->shape->positions = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall->shape->triangles = {{0, 1, 2}, {2, 3, 0}};
  rightwall->material->color  = {0.14, 0.45, 0.091};
  auto leftwall               = add_complete_object(scene, "leftwall");
  leftwall->shape->positions  = {
      {-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall->shape->triangles = {{0, 1, 2}, {2, 3, 0}};
  leftwall->material->color  = {0.63, 0.065, 0.05};
  auto shortbox              = add_complete_object(scene, "shortbox");
  shortbox->shape->positions = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  shortbox->shape->triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  shortbox->material->color  = {0.725, 0.71, 0.68};
  auto tallbox               = add_complete_object(scene, "tallbox");
  tallbox->shape->positions  = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
      {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
      {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
      {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
      {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
      {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
      {-0.71, 0.0, -0.49}};
  tallbox->shape->triangles  = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  tallbox->material->color   = {0.725, 0.71, 0.68};
  auto light                 = add_complete_object(scene, "light");
  light->shape->positions    = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
      {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  light->shape->triangles    = {{0, 1, 2}, {2, 3, 0}};
  light->material->emission  = {17, 12, 4};
}

}  // namespace yocto
