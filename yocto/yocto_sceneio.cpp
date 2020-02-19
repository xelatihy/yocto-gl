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

#ifdef __clang__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#endif

#define CGLTF_IMPLEMENTATION
#include "ext/cgltf.h"

#ifdef __clang__
#pragma clang diagnostic pop
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

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
      "texels3b:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->colorb.size().x *
               (size_t)texture->colorb.size().x;
      })));
  stats.push_back(
      "texels3f:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->colorf.size().x *
               (size_t)texture->colorf.size().y;
      })));
  stats.push_back(
      "texels1b:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->scalarb.size().x *
               (size_t)texture->scalarb.size().x;
      })));
  stats.push_back(
      "texels1f:     " + format(accumulate(scene->textures, [](auto texture) {
        return (size_t)texture->scalarf.size().x *
               (size_t)texture->scalarf.size().y;
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
      if (value->colorf.empty() && value->colorb.empty() &&
          value->scalarf.empty() && value->scalarb.empty()) {
        errs.push_back("empty texture " + value->name);
      }
    }
  };

  check_names(scene->cameras, "camera");
  check_names(scene->shapes, "shape");
  check_names(scene->subdivs, "subdiv");
  check_names(scene->shapes, "instance");
  check_names(scene->shapes, "object");
  check_names(scene->textures, "texture");
  check_names(scene->environments, "environment");
  if (!notextures) check_empty_textures(scene->textures);

  return errs;
}

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

sceneio_model::~sceneio_model() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto subdiv : subdivs) delete subdiv;
  for (auto material : materials) delete material;
  for (auto instance : instances) delete instance;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
}

// add element
sceneio_camera* add_camera(sceneio_model* scene, const string& name) {
  auto camera = scene->cameras.emplace_back(new sceneio_camera{});
  camera->name =
      name != "" ? name : ("camera" + std::to_string(scene->cameras.size()));
  return camera;
}
sceneio_environment* add_environment(sceneio_model* scene, const string& name) {
  auto environment = scene->environments.emplace_back(
      new sceneio_environment{});
  environment->name =
      name != "" ? name
                 : ("environment" + std::to_string(scene->environments.size()));
  return environment;
}
sceneio_shape* add_shape(sceneio_model* scene, const string& name) {
  auto shape  = scene->shapes.emplace_back(new sceneio_shape{});
  shape->name = name != "" ? name
                           : ("shape" + std::to_string(scene->shapes.size()));
  return shape;
}
sceneio_subdiv* add_subdiv(sceneio_model* scene, const string& name) {
  auto subdiv = scene->subdivs.emplace_back(new sceneio_subdiv{});
  subdiv->name =
      name != "" ? name : ("subdiv" + std::to_string(scene->subdivs.size()));
  return subdiv;
}
sceneio_texture* add_texture(sceneio_model* scene, const string& name) {
  auto texture = scene->textures.emplace_back(new sceneio_texture{});
  texture->name =
      name != "" ? name : ("texture" + std::to_string(scene->textures.size()));
  return texture;
}
sceneio_object* add_object(sceneio_model* scene, const string& name) {
  auto object = scene->objects.emplace_back(new sceneio_object{});
  object->name =
      name != "" ? name : ("object" + std::to_string(scene->objects.size()));
  return object;
}
sceneio_instance* add_instance(sceneio_model* scene, const string& name) {
  auto instance  = scene->instances.emplace_back(new sceneio_instance{});
  instance->name = name != ""
                       ? name
                       : ("instance" + std::to_string(scene->instances.size()));
  return instance;
}
sceneio_material* add_material(sceneio_model* scene, const string& name) {
  auto material  = scene->materials.emplace_back(new sceneio_material{});
  material->name = name != ""
                       ? name
                       : ("material" + std::to_string(scene->materials.size()));
  return material;
}
sceneio_object* add_complete_object(sceneio_model* scene, const string& name) {
  auto object      = add_object(scene, name);
  object->shape    = add_shape(scene, name);
  object->material = add_material(scene, name);
  return object;
}

// get named camera or default if camera is empty
sceneio_camera* get_camera(const sceneio_model* scene, const string& name) {
  if (scene->cameras.empty()) return nullptr;
  for (auto camera : scene->cameras) {
    if (camera->name == name) return camera;
  }
  for (auto camera : scene->cameras) {
    if (camera->name == "default") return camera;
  }
  for (auto camera : scene->cameras) {
    if (camera->name == "camera") return camera;
  }
  return scene->cameras.front();
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
  auto camera          = add_camera(scene, "camera");
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
  using yocto::image::make_sunsky;
  auto texture = add_texture(scene, "sky");
  auto sunsky  = make_sunsky({1024, 512}, sun_angle);
  texture->colorf.resize(sunsky.size());
  for (auto j = 0; j < sunsky.size().y; j++)
    for (auto i = 0; j < sunsky.size().x; i++)
      texture->colorf[{i, j}] = xyz(sunsky[{i, j}]);
  auto environment          = add_environment(scene, "sky");
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
    texture->colorf.shrink_to_fit();
    texture->colorb.shrink_to_fit();
    texture->scalarf.shrink_to_fit();
    texture->scalarb.shrink_to_fit();
  }
  scene->cameras.shrink_to_fit();
  scene->shapes.shrink_to_fit();
  scene->textures.shrink_to_fit();
  scene->environments.shrink_to_fit();
}

// Check texture size
static vec2i texture_size(const sceneio_texture* texture) {
  if (!texture->colorf.empty()) {
    return texture->colorf.size();
  } else if (!texture->colorb.empty()) {
    return texture->colorb.size();
  } else if (!texture->scalarf.empty()) {
    return texture->scalarf.size();
  } else if (!texture->scalarb.empty()) {
    return texture->scalarb.size();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
static vec3f lookup_texture(const sceneio_texture* texture, const vec2i& ij,
    bool ldr_as_linear = false) {
  if (!texture->colorf.empty()) {
    return texture->colorf[ij];
  } else if (!texture->colorb.empty()) {
    return ldr_as_linear ? byte_to_float(texture->colorb[ij])
                         : srgb_to_rgb(byte_to_float(texture->colorb[ij]));
  } else if (!texture->scalarf.empty()) {
    return vec3f{texture->scalarf[ij]};
  } else if (!texture->scalarb.empty()) {
    return ldr_as_linear
               ? byte_to_float(vec3b{texture->scalarb[ij]})
               : srgb_to_rgb(byte_to_float(vec3b{texture->scalarb[ij]}));
  } else {
    return {1, 1, 1};
  }
}

// Evaluate a texture
static vec3f eval_texture(const sceneio_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  // get texture
  if (!texture) return {1, 1, 1};

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

// Compute per-vertex normals for quads.
static vector<vec3f> compute_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto normals = vector<vec3f>{positions.size()};
  for (auto& normal : normals) normal = zero3f;
  for (auto& q : quads) {
    auto normal = quad_normal(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    auto area = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    normals[q.x] += normal * area;
    normals[q.y] += normal * area;
    normals[q.z] += normal * area;
    if (q.z != q.w) normals[q.w] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>,
    vector<vec2f>> static split_facevarying(const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords) {
  auto split =
      std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>{};
  auto& [split_quads, split_positions, split_normals, split_texcoords] = split;
  // make faces unique
  unordered_map<vec3i, int> vert_map;
  split_quads.resize(quadspos.size());
  for (auto fid = 0; fid < quadspos.size(); fid++) {
    for (auto c = 0; c < 4; c++) {
      auto v = vec3i{
          (&quadspos[fid].x)[c],
          (!quadsnorm.empty()) ? (&quadsnorm[fid].x)[c] : -1,
          (!quadstexcoord.empty()) ? (&quadstexcoord[fid].x)[c] : -1,
      };
      auto it = vert_map.find(v);
      if (it == vert_map.end()) {
        auto s = (int)vert_map.size();
        vert_map.insert(it, {v, s});
        (&split_quads[fid].x)[c] = s;
      } else {
        (&split_quads[fid].x)[c] = it->second;
      }
    }
  }

  // fill vert data
  split_positions.clear();
  if (!positions.empty()) {
    split_positions.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_positions[index] = positions[vert.x];
    }
  }
  split_normals.clear();
  if (!normals.empty()) {
    split_normals.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_normals[index] = normals[vert.y];
    }
  }
  split_texcoords.clear();
  if (!texcoords.empty()) {
    split_texcoords.resize(vert_map.size());
    for (auto& [vert, index] : vert_map) {
      split_texcoords[index] = texcoords[vert.z];
    }
  }

  return split;
}

// using
using std::unique_ptr;
using yocto::shape::subdivide_catmullclark;

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
          eval_texture(displacement_tex, subdiv->texcoords[qtxt[i]], true));
      if (!displacement_tex->scalarb.empty() ||
          !displacement_tex->colorb.empty())
        disp -= 0.5f;
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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Load/save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);
static bool save_json_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);
static bool save_obj_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);
static bool save_ply_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);
static bool save_pbrt_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel);

// Load a scene
bool load_scene(const string& filename, sceneio_model* scene, string& error,
    sceneio_progress progress_cb, bool noparallel) {
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
bool save_scene(const string& filename, const sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Get extension (not including '.').
static string get_extension(const string& filename) {
  auto pos = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Loads/saves a  channel float/byte image in linear/srgb color space.
static bool load_image(const string& filename, image<vec4f>& colorf,
    image<vec4b>& colorb, string& error) {
  using yocto::image::is_hdr_filename;
  using yocto::image::load_image;
  if (is_hdr_filename(filename)) {
    return load_image(filename, colorf, error);
  } else {
    return load_image(filename, colorb, error);
  }
}

// Loads/saves a 3 channel float/byte image in linear/srgb color space.
static bool load_image(const string& filename, image<vec3f>& colorf,
    image<vec3b>& colorb, string& error) {
  using yocto::image::is_hdr_filename;
  using yocto::image::load_image;
  if (is_hdr_filename(filename)) {
    return load_image(filename, colorf, error);
  } else {
    return load_image(filename, colorb, error);
  }
}
static bool save_image(const string& filename, const image<vec3f>& colorf,
    const image<vec3b>& colorb, string& error) {
  using yocto::image::is_hdr_filename;
  using yocto::image::save_image;
  if (is_hdr_filename(filename)) {
    return save_image(filename, colorf, error);
  } else {
    return save_image(filename, colorb, error);
  }
}

// Loads/saves a 1 channel float/byte image in linear/srgb color space.
static bool load_image(const string& filename, image<float>& scalarf,
    image<byte>& scalarb, string& error) {
  using yocto::image::is_hdr_filename;
  using yocto::image::load_image;
  if (is_hdr_filename(filename)) {
    return load_image(filename, scalarf, error);
  } else {
    return load_image(filename, scalarb, error);
  }
}
static bool save_image(const string& filename, const image<float>& scalarf,
    const image<byte>& scalarb, string& error) {
  using yocto::image::is_hdr_filename;
  using yocto::image::save_image;
  if (is_hdr_filename(filename)) {
    return save_image(filename, scalarf, error);
  } else {
    return save_image(filename, scalarb, error);
  }
}

// load instances
static bool load_instance(
    const string& filename, vector<frame3f>& frames, string& error) {
  using namespace yocto::modelio;
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, &ply, error)) return false;
    frames = get_values(&ply, "instance",
        array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz",
            "ox", "oy", "oz"});
    return true;
  } else {
    return format_error();
  }
}

// save instances
static bool save_instance(const string& filename, const vector<frame3f>& frames,
    string& error, bool ascii = false) {
  using namespace yocto::modelio;
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = get_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(&ply, frames, "instance",
        array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz",
            "ox", "oy", "oz"});
    if (!save_ply(filename, &ply, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto::math {

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

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

using json = nlohmann::json;

// Load a text file
inline bool load_text(const string& filename, string& str, string& error) {
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
inline bool save_text(
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
inline bool load_binary(
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
inline bool save_binary(
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

// load/save json
inline bool load_json(const string& filename, json& js, string& error) {
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

inline bool save_json(const string& filename, const json& js, string& error) {
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [scene, &ctexture_map, &get_value](const json& ejs,
                          const string& name, sceneio_texture*& value,
                          const string& dirname = "textures/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) {
      value = it->second;
      return true;
    }
    auto texture       = add_texture(scene, path);
    ctexture_map[path] = texture;
    value              = texture;
    return true;
  };

  // parse json reference
  auto stexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_stexture = [scene, &stexture_map, &get_value](const json& ejs,
                          const string& name, sceneio_texture*& value,
                          const string& dirname = "textures/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) {
      value = it->second;
      return true;
    }
    auto texture       = add_texture(scene, path);
    stexture_map[path] = texture;
    value              = texture;
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
    auto shape      = add_shape(scene, path);
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
    auto subdiv      = add_subdiv(scene, path);
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
    auto instance      = add_instance(scene, path);
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
    for (auto& [name, ejs] : js.at("cameras").items()) {
      auto camera  = add_camera(scene);
      camera->name = name;
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
    for (auto& [name, ejs] : js.at("environments").items()) {
      auto environment  = add_environment(scene);
      environment->name = name;
      if (!get_value(ejs, "frame", environment->frame)) return false;
      if (!get_value(ejs, "emission", environment->emission)) return false;
      if (!get_ctexture(
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
    for (auto& [name, ejs] : js.at("materials").items()) {
      auto material  = add_material(scene);
      material->name = name;
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
      if (!get_ctexture(ejs, "emission_tex", material->emission_tex))
        return false;
      if (!get_ctexture(ejs, "color_tex", material->color_tex)) return false;
      if (!get_stexture(ejs, "metallic_tex", material->metallic_tex))
        return false;
      if (!get_stexture(ejs, "specular_tex", material->specular_tex))
        return false;
      if (!get_stexture(ejs, "transmission_tex", material->transmission_tex))
        return false;
      if (!get_stexture(ejs, "roughness_tex", material->roughness_tex))
        return false;
      if (!get_stexture(ejs, "scattering_tex", material->scattering_tex))
        return false;
      if (!get_stexture(ejs, "opacity_tex", material->opacity_tex))
        return false;
      if (!get_ctexture(ejs, "normal_tex", material->normal_tex)) return false;
      if (!get_stexture(ejs, "displacement_tex", material->displacement_tex))
        return false;
      if (!get_value(ejs, "subdivisions", material->subdivisions))
        return false;  // hack fir subd
      if (!get_value(ejs, "smooth", material->smooth))
        return false;  // hack for subd
      material_map[material->name] = material;
    }
  }
  if (js.contains("objects")) {
    for (auto& [name, ejs] : js.at("objects").items()) {
      auto object  = add_object(scene);
      object->name = name;
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

  // get filename from name
  auto get_filename = [filename](const string& name, const string& group,
                          const vector<string>& extensions) {
    for (auto& extension : extensions) {
      auto filepath = fs::path(filename).parent_path() / group /
                      (name + extension);
      if (fs::exists(filepath)) return filepath;
    }
    return fs::path(filename).parent_path() / group /
           (name + extensions.front());
  };

  // using
  using yocto::shape::load_fvshape;
  using yocto::shape::load_shape;

  // load shapes
  shape_map.erase("");
  for (auto [name, shape] : shape_map) {
    if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
    auto path = get_filename(name, "shapes", {".ply", ".obj"});
    if (!load_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->positions, shape->normals, shape->texcoords,
            shape->colors, shape->radius, error))
      return dependent_error();
  }
  // load subdivs
  subdiv_map.erase("");
  for (auto [name, subdiv] : subdiv_map) {
    if (progress_cb) progress_cb("load subdiv", progress.x++, progress.y);
    auto path = get_filename(name, "subdivs", {".obj"});
    if (!load_fvshape(path, subdiv->quadspos, subdiv->quadsnorm,
            subdiv->quadstexcoord, subdiv->positions, subdiv->normals,
            subdiv->texcoords, error))
      return dependent_error();
  }
  // load textures
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto path = get_filename(
        name, "textures", {".hdr", ".exr", ".png", ".jpg"});
    if (!load_image(path, texture->colorf, texture->colorb, error))
      return dependent_error();
  }
  // load textures
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto path = get_filename(
        name, "textures", {".hdr", ".exr", ".png", ".jpg"});
    if (!load_image(path, texture->scalarf, texture->scalarb, error))
      return dependent_error();
  }
  // load instances
  instance_map.erase("");
  for (auto [name, instance] : instance_map) {
    if (progress_cb) progress_cb("load instance", progress.x++, progress.y);
    auto path = get_filename(name, "instances", {".ply"});
    if (!load_instance(path, instance->frames, error)) return dependent_error();
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
  if (!scene->cameras.empty()) js["cameras"] = json::object();
  for (auto& camera : scene->cameras) {
    auto& ejs = js["cameras"][camera->name];
    add_opt(ejs, "frame", camera->frame, def_cam.frame);
    add_opt(ejs, "ortho", camera->orthographic, def_cam.orthographic);
    add_opt(ejs, "lens", camera->lens, def_cam.lens);
    add_opt(ejs, "aspect", camera->aspect, def_cam.aspect);
    add_opt(ejs, "film", camera->film, def_cam.film);
    add_opt(ejs, "focus", camera->focus, def_cam.focus);
    add_opt(ejs, "aperture", camera->aperture, def_cam.aperture);
  }

  auto def_env = sceneio_environment{};
  if (!scene->environments.empty()) js["environments"] = json::object();
  for (auto environment : scene->environments) {
    auto& ejs = js["environments"][environment->name];
    add_opt(ejs, "frame", environment->frame, def_env.frame);
    add_opt(ejs, "emission", environment->emission, def_env.emission);
    add_tex(ejs, "emission_tex", environment->emission_tex);
  }

  auto def_material = sceneio_material{};
  if (!scene->materials.empty()) js["materials"] = json::object();
  for (auto material : scene->materials) {
    auto& ejs = js["materials"][material->name];
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
        def_material.subdivisions);  // hack for subd
    add_opt(
        ejs, "smooth", material->smooth, def_material.smooth);  // hack for subd
  }

  auto def_object = sceneio_object{};
  auto def_subdiv = sceneio_subdiv{};
  if (!scene->objects.empty()) js["objects"] = json::object();
  for (auto object : scene->objects) {
    auto& ejs = js["objects"][object->name];
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

  // get filename from name
  auto get_filename = [filename](const string& name, const string& group,
                          const string& extension) {
    return fs::path(filename).parent_path() / group / (name + extension);
  };

  // using
  using yocto::shape::save_fvshape;
  using yocto::shape::save_shape;

  // save shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    auto path = get_filename(shape->name, "shapes", ".ply");
    if (!save_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->positions, shape->normals, shape->texcoords,
            shape->colors, shape->radius, error))
      return dependent_error();
  }

  // save subdivs
  for (auto subdiv : scene->subdivs) {
    if (progress_cb) progress_cb("save subdiv", progress.x++, progress.y);
    auto path = get_filename(subdiv->name, "subdivs", ".obj");
    if (!save_fvshape(path, subdiv->quadspos, subdiv->quadsnorm,
            subdiv->quadstexcoord, subdiv->positions, subdiv->normals,
            subdiv->texcoords, error))
      return dependent_error();
  }

  // save instances
  for (auto instance : scene->instances) {
    if (progress_cb) progress_cb("save instance", progress.x++, progress.y);
    auto path = get_filename(instance->name, "instances", ".ply");
    if (!save_instance(path, instance->frames, error)) return dependent_error();
  }

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = get_filename(texture->name, "textures",
        (!texture->colorf.empty() || !texture->scalarf.empty()) ? ".hdr"
                                                                : ".png");
    if (!texture->colorf.empty() || !texture->colorb.empty()) {
      if (!save_image(path, texture->colorf, texture->colorb, error))
        return dependent_error();
    } else {
      if (!save_image(path, texture->scalarf, texture->scalarb, error))
        return dependent_error();
    }
  }

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Loads an OBJ
static bool load_obj_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  using namespace yocto::modelio;
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&ctexture_map, scene](
                          const obj_texture_info& info) -> sceneio_texture* {
    auto path = info.path;
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };

  // helper to create texture maps
  auto stexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_stexture = [&stexture_map, scene](
                          const obj_texture_info& info) -> sceneio_texture* {
    auto path = info.path;
    if (path == "") return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
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
    material->emission_tex     = get_ctexture(omat->pbr_emission_tex);
    material->color_tex        = get_ctexture(omat->pbr_base_tex);
    material->specular_tex     = get_stexture(omat->pbr_specular_tex);
    material->metallic_tex     = get_stexture(omat->pbr_metallic_tex);
    material->roughness_tex    = get_stexture(omat->pbr_roughness_tex);
    material->transmission_tex = get_stexture(omat->pbr_transmission_tex);
    material->coat_tex         = get_stexture(omat->pbr_coat_tex);
    material->opacity_tex      = get_stexture(omat->pbr_opacity_tex);
    material->normal_tex       = get_ctexture(omat->normal_tex);
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
      auto has_quads_  = has_quads(oshape);
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
        object->instance         = add_instance(scene);
        object->instance->frames = oshape->instances;
      }
    }
  }

  // convert environments
  for (auto oenvironment : obj->environments) {
    auto environment          = add_environment(scene);
    environment->frame        = oenvironment->frame;
    environment->emission     = oenvironment->emission;
    environment->emission_tex = get_ctexture(oenvironment->emission_tex);
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // get filename from name
  auto get_filename = [filename](const string& name) {
    return fs::path(filename).parent_path() / name;
  };

  // load textures
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(
            get_filename(name), texture->colorf, texture->colorb, error))
      return dependent_error();
  }

  // load textures
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(
            get_filename(name), texture->scalarf, texture->scalarb, error))
      return dependent_error();
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
  using namespace yocto::modelio;
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

  auto obj_guard = make_unique<obj_model>();
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

  // get filename from name
  auto get_filename = [filename](const string& name, const string& group,
                          const string& extension) {
    return fs::path(filename).parent_path() / group / (name + extension);
  };

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = get_filename(texture->name, "textures",
        (!texture->colorf.empty() || !texture->scalarf.empty()) ? ".hdr"
                                                                : ".png");
    if (!texture->colorf.empty() || !texture->colorb.empty()) {
      if (!save_image(path, texture->colorf, texture->colorb, error))
        return dependent_error();
    } else {
      if (!save_image(path, texture->scalarf, texture->scalarb, error))
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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

static bool load_ply_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // using
  using yocto::shape::load_shape;

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

  // using
  using yocto::shape::save_shape;

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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// Load a scene
static bool load_gltf_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto primitive_error = [filename, &error]() {
    error = filename + ": primitive error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 3};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load gltf
  auto params = cgltf_options{};
  memset(&params, 0, sizeof(params));
  auto data   = (cgltf_data*)nullptr;
  auto result = cgltf_parse_file(&params, filename.c_str(), &data);
  if (result != cgltf_result_success) return read_error();
  auto gltf = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      data, cgltf_free};

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load buffers
  auto dirname = fs::path(filename).parent_path().string();
  if (dirname != "") dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success)
    return read_error();

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde = &gltf->nodes[nid];
    if (!gnde->camera) continue;
    auto mat = mat4f{};
    cgltf_node_transform_world(gnde, &mat.x.x);
    auto gcam            = gnde->camera;
    auto camera          = add_camera(scene);
    camera->frame        = (frame3f)mat;
    camera->orthographic = gcam->type == cgltf_camera_type_orthographic;
    if (camera->orthographic) {
      auto ortho     = &gcam->data.orthographic;
      camera->aspect = ortho->xmag / ortho->ymag;
      camera->lens   = ortho->ymag;  // this is probably bogus
      camera->film   = 0.036;
    } else {
      auto persp     = &gcam->data.perspective;
      camera->aspect = persp->aspect_ratio;
      camera->film   = 0.036;
      camera->lens   = camera->aspect >= 1
                         ? (2 * camera->aspect * tan(persp->yfov / 2))
                         : (2 * tan(persp->yfov / 2));
    }
  }

  // convert color textures
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const cgltf_texture_view& ginfo) -> sceneio_texture* {
    if (!ginfo.texture || !ginfo.texture->image) return nullptr;
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  // convert color opacity textures
  auto cotexture_map =
      unordered_map<string, pair<sceneio_texture*, sceneio_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_cotexture = [&scene, &cotexture_map](const cgltf_texture_view& ginfo)
      -> pair<sceneio_texture*, sceneio_texture*> {
    if (!ginfo.texture || !ginfo.texture->image) return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return {nullptr, nullptr};
    auto it = cotexture_map.find(path);
    if (it != cotexture_map.end()) return it->second;
    auto color_texture   = add_texture(scene);
    auto opacity_texture = add_texture(scene);
    cotexture_map[path]  = {color_texture, opacity_texture};
    return {color_texture, opacity_texture};
  };
  // convert textures
  auto mrtexture_map =
      unordered_map<string, pair<sceneio_texture*, sceneio_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_mrtexture = [&scene, &mrtexture_map](const cgltf_texture_view& ginfo)
      -> pair<sceneio_texture*, sceneio_texture*> {
    if (!ginfo.texture || !ginfo.texture->image) return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return {nullptr, nullptr};
    auto it = mrtexture_map.find(path);
    if (it != mrtexture_map.end()) return it->second;
    auto metallic_texture  = add_texture(scene);
    auto roughness_texture = add_texture(scene);
    mrtexture_map[path]    = {metallic_texture, roughness_texture};
    return {metallic_texture, roughness_texture};
  };

  // convert materials
  auto material_map = unordered_map<cgltf_material*, sceneio_material*>{
      {nullptr, nullptr}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto gmaterial         = &gltf->materials[mid];
    auto material          = add_material(scene);
    material->emission     = {gmaterial->emissive_factor[0],
        gmaterial->emissive_factor[1], gmaterial->emissive_factor[2]};
    material->emission_tex = get_ctexture(gmaterial->emissive_texture);
    if (gmaterial->has_pbr_metallic_roughness) {
      auto gmr          = &gmaterial->pbr_metallic_roughness;
      material->color   = {gmr->base_color_factor[0], gmr->base_color_factor[1],
          gmr->base_color_factor[2]};
      material->opacity = gmr->base_color_factor[3];
      material->metallic  = gmr->metallic_factor;
      material->roughness = gmr->roughness_factor;
      material->specular  = 1;
      std::tie(material->color_tex, material->opacity_tex) = get_cotexture(
          gmr->base_color_texture);
      std::tie(material->metallic_tex, material->roughness_tex) = get_mrtexture(
          gmr->metallic_roughness_texture);
    }
    material->normal_tex    = get_ctexture(gmaterial->normal_texture);
    material_map[gmaterial] = material;
  }

  // convert meshes
  auto mesh_map = unordered_map<cgltf_mesh*, vector<sceneio_object*>>{
      {nullptr, {}}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh = &gltf->meshes[mid];
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (!gprim->attributes_count) continue;
      auto object = add_object(scene);
      mesh_map[gmesh].push_back(object);
      auto shape       = add_shape(scene);
      object->shape    = shape;
      object->material = material_map.at(gprim->material);
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name ? gattr->name : "");
        auto gacc     = gattr->data;
        if (semantic == "POSITION") {
          shape->positions.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->positions[i].x, 3);
        } else if (semantic == "NORMAL") {
          shape->normals.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->normals[i].x, 3);
        } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
          shape->texcoords.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->texcoords[i].x, 2);
        } else if (semantic == "COLOR" || semantic == "COLOR_0") {
          shape->colors.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->colors[i].x, 3);
        } else if (semantic == "TANGENT") {
          shape->tangents.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->tangents[i].x, 4);
          for (auto& t : shape->tangents) t.w = -t.w;
        } else if (semantic == "RADIUS") {
          shape->radius.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->radius[i], 1);
        } else {
          // ignore
        }
      }
      // indices
      if (!gprim->indices) {
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape->triangles.resize(shape->positions.size() / 3);
          for (auto i = 0; i < shape->positions.size() / 3; i++)
            shape->triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape->triangles.resize(shape->positions.size() - 2);
          for (auto i = 2; i < shape->positions.size(); i++)
            shape->triangles[i - 2] = {0, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape->triangles.resize(shape->positions.size() - 2);
          for (auto i = 2; i < shape->positions.size(); i++)
            shape->triangles[i - 2] = {i - 2, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape->lines.resize(shape->positions.size() / 2);
          for (auto i = 0; i < shape->positions.size() / 2; i++)
            shape->lines[i] = {i * 2 + 0, i * 2 + 1};
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape->lines.resize(shape->positions.size());
          for (auto i = 1; i < shape->positions.size(); i++)
            shape->lines[i - 1] = {i - 1, i};
          shape->lines.back() = {(int)shape->positions.size() - 1, 0};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape->lines.resize(shape->positions.size() - 1);
          for (auto i = 1; i < shape->positions.size(); i++)
            shape->lines[i - 1] = {i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          return primitive_error();
        } else {
          return primitive_error();
        }
      } else {
        auto giacc = gprim->indices;
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape->triangles.resize(giacc->count / 3);
          for (auto i = 0; i < giacc->count / 3; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 3 + 0, (uint*)&shape->triangles[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 1, (uint*)&shape->triangles[i].y, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 2, (uint*)&shape->triangles[i].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape->triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, 0 + 0, (uint*)&shape->triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape->triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape->triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape->triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, i - 2, (uint*)&shape->triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape->triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape->triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape->lines.resize(giacc->count / 2);
          for (auto i = 0; i < giacc->count / 2; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 2 + 0, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 2 + 1, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape->lines.resize(giacc->count);
          for (auto i = 0; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape->lines.resize(giacc->count - 1);
          for (auto i = 0; i < giacc->count - 1; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          return primitive_error();
        } else {
          return primitive_error();
        }
      }
    }
  }

  // convert nodes
  auto instance_map = unordered_map<cgltf_mesh*, vector<frame3f>>{};
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde = &gltf->nodes[nid];
    if (!gnde->mesh) continue;
    auto mat = mat4f{};
    cgltf_node_transform_world(gnde, &mat.x.x);
    auto frame = (frame3f)mat;
    instance_map[gnde->mesh].push_back(frame);
  }
  for (auto& [gmsh, frames] : instance_map) {
    if (frames.size() == 1) {
      for (auto object : mesh_map.at(gmsh)) object->frame = frames.front();
    } else {
      auto instance    = add_instance(scene);
      instance->frames = frames;
      for (auto object : mesh_map.at(gmsh)) object->instance = instance;
    }
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // load texture
  ctexture_map.erase("");
  for (auto [path, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(fs::path(filename).parent_path() / path, texture->colorf,
            texture->colorb, error))
      return dependent_error();
  }

  // load texture
  cotexture_map.erase("");
  for (auto [path, textures] : cotexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto color_opacityf = image<vec4f>{};
    auto color_opacityb = image<vec4b>{};
    if (!load_image(fs::path(filename).parent_path() / path, color_opacityf,
            color_opacityb, error))
      return dependent_error();
    if (!color_opacityf.empty()) {
      auto [ctexture, otexture] = textures;
      ctexture->colorf.resize(color_opacityf.size());
      otexture->scalarf.resize(color_opacityf.size());
      for (auto j = 0; j < color_opacityf.size().y; j++) {
        for (auto i = 0; i < color_opacityf.size().x; i++) {
          ctexture->colorf[{i, j}]  = xyz(color_opacityf[{i, j}]);
          otexture->scalarf[{i, j}] = color_opacityf[{i, j}].w;
        }
      }
    }
    if (!color_opacityb.empty()) {
      auto [ctexture, otexture] = textures;
      ctexture->colorb.resize(color_opacityb.size());
      otexture->scalarb.resize(color_opacityb.size());
      for (auto j = 0; j < color_opacityb.size().y; j++) {
        for (auto i = 0; i < color_opacityb.size().x; i++) {
          ctexture->colorb[{i, j}]  = xyz(color_opacityb[{i, j}]);
          otexture->scalarb[{i, j}] = color_opacityb[{i, j}].w;
        }
      }
    }
  }

  // load texture
  mrtexture_map.erase("");
  for (auto [path, textures] : mrtexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto metallic_roughnessf = image<vec3f>{};
    auto metallic_roughnessb = image<vec3b>{};
    if (!load_image(fs::path(filename).parent_path() / path,
            metallic_roughnessf, metallic_roughnessb, error))
      return dependent_error();
    if (!metallic_roughnessf.empty()) {
      auto [mtexture, rtexture] = textures;
      mtexture->scalarf.resize(metallic_roughnessf.size());
      rtexture->scalarf.resize(metallic_roughnessf.size());
      for (auto j = 0; j < metallic_roughnessf.size().y; j++) {
        for (auto i = 0; i < metallic_roughnessf.size().x; i++) {
          mtexture->scalarf[{i, j}] = metallic_roughnessf[{i, j}].z;
          rtexture->scalarf[{i, j}] = metallic_roughnessf[{i, j}].y;
        }
      }
    }
    if (!metallic_roughnessb.empty()) {
      auto [mtexture, rtexture] = textures;
      mtexture->scalarb.resize(metallic_roughnessb.size());
      rtexture->scalarb.resize(metallic_roughnessb.size());
      for (auto j = 0; j < metallic_roughnessb.size().y; j++) {
        for (auto i = 0; i < metallic_roughnessb.size().x; i++) {
          mtexture->scalarb[{i, j}] = metallic_roughnessb[{i, j}].z;
          rtexture->scalarb[{i, j}] = metallic_roughnessb[{i, j}].y;
        }
      }
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

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

// load pbrt scenes
static bool load_pbrt_scene(const string& filename, sceneio_model* scene,
    string& error, sceneio_progress progress_cb, bool noparallel) {
  using namespace yocto::modelio;
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const string& path) -> sceneio_texture* {
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  auto stexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_stexture = [&scene, &stexture_map](
                          const string& path) -> sceneio_texture* {
    if (path == "") return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
    return texture;
  };
  auto atexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_atexture = [&scene, &atexture_map](
                          const string& path) -> sceneio_texture* {
    if (path == "") return nullptr;
    auto it = atexture_map.find(path);
    if (it != atexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    atexture_map[path] = texture;
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
    material->color_tex    = get_ctexture(pmaterial->color_tex);
    material->opacity_tex  = get_stexture(pmaterial->opacity_tex);
    if (!material->opacity_tex)
      material->opacity_tex = get_atexture(pmaterial->alpha_tex);
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
    environment->emission_tex = get_ctexture(penvironment->emission_tex);
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

  // get filename from name
  auto get_filename = [filename](const string& name) {
    return fs::path(filename).parent_path() / name;
  };

  // load texture
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(
            get_filename(name), texture->colorf, texture->colorb, error))
      return dependent_error();
  }

  // load texture
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(
            get_filename(name), texture->scalarf, texture->scalarb, error))
      return dependent_error();
  }

  // load alpha
  atexture_map.erase("");
  for (auto [name, texture] : atexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(
            get_filename(name), texture->scalarf, texture->scalarb, error))
      return dependent_error();
    for (auto& c : texture->scalarf) c = (c < 0.01) ? 1 : 1;
    for (auto& c : texture->scalarb) c = (c < 2) ? 0 : 255;
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
  using namespace yocto::modelio;
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt_guard = make_unique<pbrt_model>();
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

  // get filename from name
  auto get_filename = [filename](const string& name, const string& group,
                          const string& extension) {
    return fs::path(filename).parent_path() / group / (name + extension);
  };

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = get_filename(texture->name, "textures",
        (!texture->colorf.empty() || !texture->scalarf.empty()) ? ".hdr"
                                                                : ".png");
    if (!texture->colorf.empty() || !texture->colorb.empty()) {
      if (!save_image(path, texture->colorf, texture->colorb, error))
        return dependent_error();
    } else {
      if (!save_image(path, texture->scalarf, texture->scalarb, error))
        return dependent_error();
    }
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto::sceneio

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto::sceneio {

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

}  // namespace yocto::sceneio
