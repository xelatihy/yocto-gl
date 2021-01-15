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

vector<string> scene_stats(const scene_scene& scene, bool verbose) {
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
  stats.push_back("environments: " + format(scene.environments.size()));
  stats.push_back("textures:     " + format(scene.textures.size()));
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
  stats.push_back("fvquads:     " +
                  format(accumulate(scene.subdivs,
                      [](auto& subdiv) { return subdiv.quadspos.size(); })));
  stats.push_back(
      "texels4b:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.ldr.width() * (size_t)texture.ldr.width();
      })));
  stats.push_back(
      "texels4f:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.hdr.width() * (size_t)texture.hdr.height();
      })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene.
vector<string> scene_validation(const scene_scene& scene, bool notextures) {
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
  auto check_empty_textures = [&errs](const scene_scene& scene) {
    for (auto idx = 0; idx < (int)scene.textures.size(); idx++) {
      auto& texture = scene.textures[idx];
      if (texture.hdr.empty() && texture.ldr.empty()) {
        errs.push_back("empty texture " + scene.texture_names[idx]);
      }
    }
  };

  check_names(scene.camera_names, "camera");
  check_names(scene.shape_names, "shape");
  check_names(scene.material_names, "material");
  check_names(scene.instance_names, "instance");
  check_names(scene.texture_names, "texture");
  check_names(scene.environment_names, "environment");
  if (!notextures) check_empty_textures(scene);

  return errs;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

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
camera_handle add_camera(scene_scene& scene, const string& name) {
  return add_element(scene.cameras, scene.camera_names, name, "camera");
}
environment_handle add_environment(scene_scene& scene, const string& name) {
  return add_element(
      scene.environments, scene.environment_names, name, "environment");
}
shape_handle add_shape(scene_scene& scene, const string& name) {
  return add_element(scene.shapes, scene.shape_names, name, "shape");
}
texture_handle add_texture(scene_scene& scene, const string& name) {
  return add_element(scene.textures, scene.texture_names, name, "texture");
}
instance_handle add_instance(scene_scene& scene, const string& name) {
  return add_element(scene.instances, scene.instance_names, name, "instance");
}
material_handle add_material(scene_scene& scene, const string& name) {
  return add_element(scene.materials, scene.material_names, name, "material");
}
subdiv_handle add_subdiv(scene_scene& scene, const string& name) {
  return add_element(scene.subdivs, scene.subdiv_names, name, "subdiv");
}
instance_handle add_complete_instance(scene_scene& scene, const string& name) {
  auto  handle      = add_instance(scene, name);
  auto& instance    = scene.instances[handle];
  instance.shape    = add_shape(scene, name);
  instance.material = add_material(scene, name);
  return handle;
}

camera_handle add_camera(scene_scene& scene, const string& name,
    const vec3f& from, const vec3f& to, const vec3f& up, float lens,
    float aspect, float aperture, bool orthographic, float film) {
  scene.camera_names.emplace_back(name);
  auto& camera        = scene.cameras.emplace_back();
  camera.frame        = lookat_frame(from, to, up);
  camera.lens         = lens;
  camera.aspect       = aspect;
  camera.film         = film;
  camera.orthographic = orthographic;
  camera.aperture     = aperture;
  camera.focus        = length(from - to);
  return (int)scene.cameras.size() - 1;
}
camera_handle add_camera(scene_scene& scene, const string& name,
    const frame3f& frame, float lens, float aspect, float aperture, float focus,
    bool orthographic, float film) {
  scene.camera_names.emplace_back(name);
  auto& camera        = scene.cameras.emplace_back();
  camera.frame        = frame;
  camera.lens         = lens;
  camera.aspect       = aspect;
  camera.film         = film;
  camera.orthographic = orthographic;
  camera.aperture     = aperture;
  camera.focus        = focus;
  return (int)scene.cameras.size() - 1;
}
instance_handle add_instance(scene_scene& scene, const string& name,
    const frame3f& frame, shape_handle shape, material_handle material) {
  scene.instance_names.emplace_back(name);
  auto& instance    = scene.instances.emplace_back();
  instance.frame    = frame;
  instance.shape    = shape;
  instance.material = material;
  return (int)scene.instances.size() - 1;
}
environment_handle add_environment(scene_scene& scene, const string& name,
    const frame3f& frame, const vec3f& emission, texture_handle emission_tex) {
  scene.environment_names.emplace_back(name);
  auto& environment        = scene.environments.emplace_back();
  environment.frame        = frame;
  environment.emission     = emission;
  environment.emission_tex = emission_tex;
  return (int)scene.environments.size() - 1;
  ;
}
texture_handle add_texture(scene_scene& scene, const string& name,
    const image<vec4f>& img, bool hdr, bool ldr_linear) {
  scene.texture_names.emplace_back(name);
  auto& texture = scene.textures.emplace_back();
  if (hdr) {
    texture.hdr = img;
  } else {
    texture.ldr = ldr_linear ? float_to_byte(img) : rgb_to_srgbb(img);
  }
  return (int)scene.textures.size() - 1;
}
shape_handle add_shape(
    scene_scene& scene, const string& name, const quads_shape& shape_data) {
  scene.shape_names.emplace_back(name);
  auto& shape     = scene.shapes.emplace_back();
  shape.points    = shape_data.points;
  shape.lines     = shape_data.lines;
  shape.triangles = shape_data.triangles;
  shape.quads     = shape_data.quads;
  shape.positions = shape_data.positions;
  shape.normals   = shape_data.normals;
  shape.texcoords = shape_data.texcoords;
  shape.colors    = shape_data.colors;
  shape.radius    = shape_data.radius;
  return (int)scene.shapes.size() - 1;
}
shape_handle add_shape(
    scene_scene& scene, const string& name, const fvshape_data& shape_data) {
  scene.shape_names.emplace_back(name);
  auto& shape = scene.shapes.emplace_back();
  std::tie(shape.quads, shape.positions, shape.normals, shape.texcoords) =
      split_facevarying(shape_data.quadspos, shape_data.quadsnorm,
          shape_data.quadstexcoord, shape_data.positions, shape_data.normals,
          shape_data.texcoords);
  return (int)scene.shapes.size() - 1;
}
subdiv_handle add_subdiv(scene_scene& scene, const string& name,
    const quads_shape& shape_data, shape_handle shape, int subdivisions,
    float displacement, texture_handle displacement_tex) {
  scene.subdiv_names.emplace_back(name);
  auto& subdiv    = scene.subdivs.emplace_back();
  auto& quads     = (!shape_data.quads.empty())
                        ? shape_data.quads
                        : triangles_to_quads(shape_data.triangles);
  subdiv.quadspos = quads;
  if (!shape_data.normals.empty()) subdiv.quadsnorm = quads;
  if (!shape_data.texcoords.empty()) subdiv.quadstexcoord = quads;
  subdiv.positions        = shape_data.positions;
  subdiv.normals          = shape_data.normals;
  subdiv.texcoords        = shape_data.texcoords;
  subdiv.shape            = shape;
  subdiv.subdivisions     = subdivisions;
  subdiv.smooth           = subdivisions > 0 || displacement_tex;
  subdiv.displacement     = displacement;
  subdiv.displacement_tex = displacement_tex;
  return (int)scene.subdivs.size() - 1;
}
subdiv_handle add_subdiv(scene_scene& scene, const string& name,
    const quads_fvshape& subdiv_data, shape_handle shape, int subdivisions,
    float displacement, texture_handle displacement_tex) {
  scene.subdiv_names.emplace_back(name);
  auto& subdiv            = scene.subdivs.emplace_back();
  subdiv.quadspos         = subdiv_data.quadspos;
  subdiv.quadsnorm        = subdiv_data.quadsnorm;
  subdiv.quadstexcoord    = subdiv_data.quadstexcoord;
  subdiv.positions        = subdiv_data.positions;
  subdiv.normals          = subdiv_data.normals;
  subdiv.texcoords        = subdiv_data.texcoords;
  subdiv.shape            = shape;
  subdiv.subdivisions     = subdivisions;
  subdiv.smooth           = subdivisions > 0 || displacement_tex;
  subdiv.displacement     = displacement;
  subdiv.displacement_tex = displacement_tex;
  return (int)scene.subdivs.size() - 1;
}
material_handle add_emission_material(scene_scene& scene, const string& name,
    const vec3f& emission, texture_handle emission_tex) {
  scene.material_names.emplace_back(name);
  auto& material        = scene.materials.emplace_back();
  material.emission     = emission;
  material.emission_tex = emission_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_matte_material(scene_scene& scene, const string& name,
    const vec3f& color, texture_handle color_tex, texture_handle normal_tex) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::matte;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = 1;
  material.metallic      = 0;
  material.ior           = 1.5;
  material.color_tex     = color_tex;
  material.roughness_tex = invalid_handle;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_plastic_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, texture_handle color_tex,
    texture_handle roughness_tex, texture_handle normal_tex, float ior) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::plastic;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = roughness;
  material.ior           = ior;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_metal_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, texture_handle color_tex,
    texture_handle roughness_tex, texture_handle normal_tex) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::metal;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = roughness;
  material.metallic      = 1;
  material.ior           = 1.5;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_metallic_material(scene_scene& scene, const string& name,
    const vec3f& color, texture_handle color_tex, float roughness,
    float metallic, texture_handle roughness_tex, texture_handle normal_tex,
    texture_handle metallic_tex, float ior) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::metallic;
  material.color         = color;
  material.opacity       = 1;
  material.metallic      = metallic;
  material.roughness     = roughness;
  material.ior           = ior;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_thinglass_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, texture_handle color_tex,
    texture_handle roughness_tex, texture_handle normal_tex, float ior) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::thinglass;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = roughness;
  material.metallic      = 0;
  material.ior           = ior;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_glass_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, texture_handle color_tex,
    texture_handle roughness_tex, texture_handle normal_tex, float ior) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::glass;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = roughness;
  material.metallic      = 0;
  material.ior           = ior;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_glass_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, const vec3f& scattering,
    texture_handle color_tex, texture_handle roughness_tex,
    texture_handle normal_tex, float ior, float scanisotropy, float trdepth) {
  scene.material_names.emplace_back(name);
  auto& material         = scene.materials.emplace_back();
  material.type          = material_type::glass;
  material.color         = color;
  material.opacity       = 1;
  material.roughness     = roughness;
  material.metallic      = 0;
  material.ior           = ior;
  material.scattering    = scattering;
  material.scanisotropy  = scanisotropy;
  material.trdepth       = trdepth;
  material.color_tex     = color_tex;
  material.roughness_tex = roughness_tex;
  material.normal_tex    = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_subsurface_material(scene_scene& scene, const string& name,
    const vec3f& color, float roughness, const vec3f& scattering,
    texture_handle color_tex, texture_handle roughness_tex,
    const texture_handle scattering_tex, texture_handle normal_tex, float ior,
    float scanisotropy, float trdepth) {
  scene.material_names.emplace_back(name);
  auto& material          = scene.materials.emplace_back();
  material.type           = material_type::subsurface;
  material.color          = color;
  material.color_tex      = color_tex;
  material.roughness      = roughness;
  material.roughness_tex  = roughness_tex;
  material.scattering     = scattering;
  material.scattering_tex = scattering_tex;
  material.ior            = ior;
  material.scanisotropy   = scanisotropy;
  material.trdepth        = trdepth;
  material.normal_tex     = normal_tex;
  return (int)scene.materials.size() - 1;
}
material_handle add_volume_material(scene_scene& scene, const string& name,
    const vec3f& color, const vec3f& scattering, float scanisotropy,
    float trdepth) {
  scene.material_names.emplace_back(name);
  auto& material        = scene.materials.emplace_back();
  material.type         = material_type::volume;
  material.color        = color;
  material.scattering   = scattering;
  material.scanisotropy = scanisotropy;
  material.trdepth      = trdepth;
  material.roughness    = 0;
  material.ior          = 1;
  material.opacity      = 1;
  return (int)scene.materials.size() - 1;
}
material_handle add_transparent_material(scene_scene& scene, const string& name,
    const vec3f& color, float opacity, texture_handle color_tex,
    texture_handle normal_tex) {
  scene.material_names.emplace_back(name);
  auto& material      = scene.materials.emplace_back();
  material.type       = material_type::matte;
  material.color      = color;
  material.color_tex  = color_tex;
  material.roughness  = 1;
  material.opacity    = opacity;
  material.normal_tex = normal_tex;
  return (int)scene.materials.size() - 1;
}

// Add missing cameras.
void add_cameras(scene_scene& scene) {
  if (!scene.cameras.empty()) return;
  scene.camera_names.emplace_back("camera");
  auto& camera        = scene.cameras.emplace_back();
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
void add_radius(scene_scene& scene, float radius) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add missing materials.
void add_materials(scene_scene& scene) {
  auto default_material = invalid_handle;
  for (auto& instance : scene.instances) {
    if (instance.material != invalid_handle) continue;
    if (default_material == invalid_handle) {
      if (!scene.instance_names.empty())
        scene.material_names.emplace_back("default");
      auto& material   = scene.materials.emplace_back();
      material.color   = {0.8, 0.8, 0.8};
      default_material = (int)scene.materials.size() - 1;
    }
    instance.material = default_material;
  }
}

// Add a sky environment
void add_sky(scene_scene& scene, float sun_angle) {
  scene.texture_names.emplace_back("sky");
  auto& texture = scene.textures.emplace_back();
  texture.hdr   = make_sunsky({1024, 512}, sun_angle);
  scene.environment_names.emplace_back("sky");
  auto& environment        = scene.environments.emplace_back();
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
}

// get named camera or default if camera is empty
camera_handle find_camera(const scene_scene& scene, const string& name) {
  if (scene.cameras.empty()) return invalid_handle;
  if (scene.camera_names.empty()) return 0;
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == name) return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "default") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera1") return idx;
  }
  return 0;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_scene& scene) {
  auto shape_bbox = vector<bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto& shape : scene.shapes) {
    auto sbvh = shape_bbox.emplace_back();
    for (auto p : shape.positions) sbvh = merge(sbvh, p);
  }
  for (auto& instance : scene.instances) {
    auto& sbvh = shape_bbox[instance.shape];
    bbox       = merge(bbox, transform_bbox(instance.frame, sbvh));
  }
  return bbox;
}

// Reduce memory usage
void trim_memory(scene_scene& scene) {
  for (auto& shape : scene.shapes) {
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
  }
  for (auto& subdiv : scene.subdivs) {
    subdiv.positions.shrink_to_fit();
    subdiv.normals.shrink_to_fit();
    subdiv.texcoords.shrink_to_fit();
    subdiv.quadspos.shrink_to_fit();
    subdiv.quadsnorm.shrink_to_fit();
    subdiv.quadstexcoord.shrink_to_fit();
  }
  for (auto& texture : scene.textures) {
    texture.hdr.shrink_to_fit();
    texture.ldr.shrink_to_fit();
  }
  scene.cameras.shrink_to_fit();
  scene.shapes.shrink_to_fit();
  scene.subdivs.shrink_to_fit();
  scene.instances.shrink_to_fit();
  scene.materials.shrink_to_fit();
  scene.textures.shrink_to_fit();
  scene.environments.shrink_to_fit();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

void tesselate_subdiv(
    scene_shape& shape, scene_subdiv& subdiv_, const scene_scene& scene) {
  auto subdiv = subdiv_;

  if (subdiv.subdivisions > 0) {
    if (subdiv.catmullclark) {
      std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_catmullclark(
          subdiv.quadstexcoord, subdiv.texcoords, subdiv.subdivisions, true);
      std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_catmullclark(
          subdiv.quadsnorm, subdiv.normals, subdiv.subdivisions, true);
      std::tie(subdiv.quadspos, subdiv.positions) = subdivide_catmullclark(
          subdiv.quadspos, subdiv.positions, subdiv.subdivisions);
    } else {
      std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_quads(
          subdiv.quadstexcoord, subdiv.texcoords, subdiv.subdivisions);
      std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_quads(
          subdiv.quadsnorm, subdiv.normals, subdiv.subdivisions);
      std::tie(subdiv.quadspos, subdiv.positions) = subdivide_quads(
          subdiv.quadspos, subdiv.positions, subdiv.subdivisions);
    }
    if (subdiv.smooth) {
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    } else {
      subdiv.normals   = {};
      subdiv.quadsnorm = {};
    }
  }

  if (subdiv.displacement != 0 && subdiv.displacement_tex != invalid_handle) {
    if (subdiv.texcoords.empty())
      throw std::runtime_error("missing texture coordinates");

    // facevarying case
    auto offset = vector<float>(subdiv.positions.size(), 0);
    auto count  = vector<int>(subdiv.positions.size(), 0);
    for (auto fid = 0; fid < subdiv.quadspos.size(); fid++) {
      auto qpos = subdiv.quadspos[fid];
      auto qtxt = subdiv.quadstexcoord[fid];
      for (auto i = 0; i < 4; i++) {
        auto& displacement_tex = scene.textures[subdiv.displacement_tex];
        auto  disp             = mean(
            eval_texture(displacement_tex, subdiv.texcoords[qtxt[i]], true));
        if (!displacement_tex.ldr.empty()) disp -= 0.5f;
        offset[qpos[i]] += subdiv.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = quads_normals(subdiv.quadspos, subdiv.positions);
    for (auto vid = 0; vid < subdiv.positions.size(); vid++) {
      subdiv.positions[vid] += normals[vid] * offset[vid] / count[vid];
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      subdiv.quadsnorm = subdiv.quadspos;
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
    }
  }

  shape = {};
  std::tie(shape.quads, shape.positions, shape.normals, shape.texcoords) =
      split_facevarying(subdiv.quadspos, subdiv.quadsnorm, subdiv.quadstexcoord,
          subdiv.positions, subdiv.normals, subdiv.texcoords);
}

void tesselate_shapes(
    scene_scene& scene, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, (int)scene.subdivs.size() + 1};
  if (progress_cb) progress_cb("tesselate subdivs", progress.x++, progress.y);

  // tesselate shapes
  for (auto& subdiv : scene.subdivs) {
    if (progress_cb) progress_cb("tesselate subdiv", progress.x++, progress.y);
    tesselate_subdiv(scene.shapes[subdiv.shape], subdiv, scene);
  }

  // done
  if (progress_cb) progress_cb("tesselate subdivs", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Check texture size
vec2i texture_size(const scene_texture& texture) {
  if (!texture.hdr.empty()) {
    return texture.hdr.imsize();
  } else if (!texture.ldr.empty()) {
    return texture.ldr.imsize();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
vec4f lookup_texture(
    const scene_texture& texture, const vec2i& ij, bool ldr_as_linear) {
  if (!texture.hdr.empty()) {
    return texture.hdr[ij];
  } else if (!texture.ldr.empty()) {
    return ldr_as_linear ? byte_to_float(texture.ldr[ij])
                         : srgb_to_rgb(byte_to_float(texture.ldr[ij]));
  } else {
    return {1, 1, 1, 1};
  }
}

// Evaluate a texture
vec4f eval_texture(const scene_texture& texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
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
vec4f eval_texture(const scene_scene& scene, texture_handle texture,
    const vec2f& uv, bool ldr_as_linear, bool no_interpolation,
    bool clamp_to_edge) {
  if (texture == invalid_handle) return {1, 1, 1, 1};
  return eval_texture(
      scene.textures[texture], uv, ldr_as_linear, no_interpolation);
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
vec3f eval_position(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
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
    const scene_scene& scene, const scene_instance& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
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
vec3f eval_normal(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
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
vec2f eval_texcoord(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
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
    const scene_scene& scene, const scene_instance& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
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

vec3f eval_normalmap(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalid_handle &&
      (!shape.triangles.empty() || !shape.quads.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, true));
    auto [tu, tv]    = eval_element_tangents(scene, instance, element);
    auto frame       = frame3f{tu, tv, normal, zero3f};
    frame.x          = orthonormalize(frame.x, frame.z);
    frame.y          = normalize(cross(frame.z, frame.x));
    auto flip_v      = dot(frame.y, tv) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading normal
vec3f eval_shading_normal(const scene_scene& scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalid_handle) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::glass) return normal;
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
vec4f eval_color(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
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
vec3f eval_environment(const scene_scene& scene,
    const scene_environment& environment, const vec3f& direction) {
  auto wl       = transform_direction(inverse(environment.frame), direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_scene& scene, const vec3f& direction) {
  auto emission = zero3f;
  for (auto environment : scene.environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

// constant values
static const auto coat_ior       = 1.5f;
static const auto coat_roughness = 0.03f * 0.03f;

// Evaluate material
material_point eval_material(const scene_scene& scene,
    const scene_instance& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, false);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, false);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, true);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, false);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::glass ||
      material.type == material_type::volume ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::metallic ||
      point.type == material_type::plastic) {
    point.roughness = clamp(point.roughness, coat_roughness, 1.0f);
  }

  return point;
}

// Evaluate material
material_point eval_material(const scene_scene& scene,
    const scene_material& material, const vec2f& texcoord,
    const vec4f& color_shp) {
  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, false);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, false);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, true);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, false);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::glass ||
      material.type == material_type::volume ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::metallic ||
      point.type == material_type::plastic) {
    point.roughness = clamp(point.roughness, coat_roughness, 1.0f);
  }

  return point;
}

// check if a material is a delta
bool is_delta(const scene_material& material) {
  return (material.type == material_type::metal && material.roughness == 0) ||
         (material.type == material_type::glass && material.roughness == 0) ||
         (material.type == material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == material_type::volume);
}
bool is_volumetric(const scene_material& material) {
  return material.type == material_type::glass ||
         material.type == material_type::volume ||
         material.type == material_type::subsurface;
}
bool is_volumetric(const scene_scene& scene, const scene_instance& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

// check if a brdf is a delta
bool is_delta(const material_point& material) {
  return (material.type == material_type::metal && material.roughness == 0) ||
         (material.type == material_type::glass && material.roughness == 0) ||
         (material.type == material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == material_type::volume);
}
bool has_volume(const material_point& material) {
  return material.type == material_type::glass ||
         material.type == material_type::volume ||
         material.type == material_type::subsurface;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox(scene_scene& scene) {
  scene.asset.name = "cornellbox";
  auto& camera     = scene.cameras[add_camera(scene, "camera")];
  camera.frame     = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 3.9}};
  camera.lens      = 0.035;
  camera.aperture  = 0.0;
  camera.focus     = 3.9;
  camera.film      = 0.024;
  camera.aspect    = 1;
  auto& floor      = scene.instances[add_complete_instance(scene, "floor")];
  scene.shapes[floor.shape].positions = {
      {-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  scene.shapes[floor.shape].triangles   = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[floor.material].color = {0.725, 0.71, 0.68};
  auto& ceiling = scene.instances[add_complete_instance(scene, "ceiling")];
  scene.shapes[ceiling.shape].positions = {
      {-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  scene.shapes[ceiling.shape].triangles   = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[ceiling.material].color = {0.725, 0.71, 0.68};
  auto& backwall = scene.instances[add_complete_instance(scene, "backwall")];
  scene.shapes[backwall.shape].positions = {
      {-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  scene.shapes[backwall.shape].triangles   = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[backwall.material].color = {0.725, 0.71, 0.68};
  auto& rightwall = scene.instances[add_complete_instance(scene, "rightwall")];
  scene.shapes[rightwall.shape].positions = {
      {1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  scene.shapes[rightwall.shape].triangles   = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[rightwall.material].color = {0.14, 0.45, 0.091};
  auto& leftwall = scene.instances[add_complete_instance(scene, "leftwall")];
  scene.shapes[leftwall.shape].positions = {
      {-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  scene.shapes[leftwall.shape].triangles   = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[leftwall.material].color = {0.63, 0.065, 0.05};
  auto& shortbox = scene.instances[add_complete_instance(scene, "shortbox")];
  scene.shapes[shortbox.shape].positions = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  scene.shapes[shortbox.shape].triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6},
      {6, 7, 4}, {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12},
      {16, 17, 18}, {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  scene.materials[shortbox.material].color = {0.725, 0.71, 0.68};
  auto& tallbox = scene.instances[add_complete_instance(scene, "tallbox")];
  scene.shapes[tallbox.shape].positions   = {{-0.53, 1.2, 0.09},
      {0.04, 1.2, -0.09}, {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49},
      {-0.53, 0.0, 0.09}, {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49},
      {-0.14, 1.2, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67},
      {-0.14, 1.2, -0.67}, {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 0.0, -0.09}, {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09},
      {-0.53, 0.0, 0.09}, {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09},
      {-0.14, 0.0, -0.67}, {-0.71, 0.0, -0.49}};
  scene.shapes[tallbox.shape].triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6},
      {6, 7, 4}, {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12},
      {16, 17, 18}, {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  scene.materials[tallbox.material].color = {0.725, 0.71, 0.68};
  auto& light = scene.instances[add_complete_instance(scene, "light")];
  scene.shapes[light.shape].positions      = {{-0.25, 1.99, 0.25},
      {-0.25, 1.99, -0.25}, {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  scene.shapes[light.shape].triangles      = {{0, 1, 2}, {2, 3, 0}};
  scene.materials[light.material].emission = {17, 12, 4};
}

}  // namespace yocto
