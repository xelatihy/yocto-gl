//
// Implementation for Yocto/Scene Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include "yocto_scene.h"

#include <atomic>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <deque>
#include <future>
#include <memory>

#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::deque;
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

vector<string> scene_stats(const scene_model* scene, bool verbose) {
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
  stats.push_back("fvquads:     " +
                  format(accumulate(scene->shapes,
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
vector<string> scene_validation(const scene_model* scene, bool notextures) {
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
  auto check_empty_textures = [&errs](const vector<scene_texture*>& vals) {
    for (auto value : vals) {
      if (value->colorf.empty() && value->colorb.empty() &&
          value->scalarf.empty() && value->scalarb.empty()) {
        errs.push_back("empty texture " + value->name);
      }
    }
  };

  check_names(scene->cameras, "camera");
  check_names(scene->shapes, "shape");
  check_names(scene->shapes, "instance");
  check_names(scene->shapes, "object");
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

scene_shape::~scene_shape() {
  if (bvh) delete bvh;
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene(embree_bvh);
#endif
}

scene_model::~scene_model() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto object : objects) delete object;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
  for (auto light : lights) delete light;
  if (bvh) delete bvh;
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene(embree_bvh);
#endif
}

// add an element
template <typename T>
static T* add_element(
    vector<T*>& elements, const string& name, const string& base) {
  auto element  = elements.emplace_back(new T{});
  element->name = name != "" ? name : (base + std::to_string(elements.size()));
  return element;
}

// add element
scene_camera* add_camera(scene_model* scene, const string& name) {
  return add_element(scene->cameras, name, "camera");
}
scene_environment* add_environment(scene_model* scene, const string& name) {
  return add_element(scene->environments, name, "environment");
}
scene_shape* add_shape(scene_model* scene, const string& name) {
  return add_element(scene->shapes, name, "shape");
}
scene_texture* add_texture(scene_model* scene, const string& name) {
  return add_element(scene->textures, name, "texture");
}
scene_object* add_object(scene_model* scene, const string& name) {
  return add_element(scene->objects, name, "object");
}
scene_material* add_material(scene_model* scene, const string& name) {
  return add_element(scene->materials, name, "material");
}
scene_object* add_complete_object(scene_model* scene, const string& name) {
  auto object      = add_object(scene, name);
  object->shape    = add_shape(scene, name);
  object->material = add_material(scene, name);
  return object;
}

// Set cameras
void set_frame(scene_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(
    scene_camera* camera, float lens, float aspect, float film, bool ortho) {
  camera->lens         = lens;
  camera->aspect       = aspect;
  camera->film         = film;
  camera->orthographic = ortho;
}
void set_focus(scene_camera* camera, float aperture, float focus) {
  camera->aperture = aperture;
  camera->focus    = focus;
}

// Add texture
void set_texture(scene_texture* texture, const image<vec3b>& img) {
  texture->colorb  = img;
  texture->colorf  = {};
  texture->scalarb = {};
  texture->scalarf = {};
}
void set_texture(scene_texture* texture, const image<vec3f>& img) {
  texture->colorb  = {};
  texture->colorf  = img;
  texture->scalarb = {};
  texture->scalarf = {};
}
void set_texture(scene_texture* texture, const image<byte>& img) {
  texture->colorb  = {};
  texture->colorf  = {};
  texture->scalarb = img;
  texture->scalarf = {};
}
void set_texture(scene_texture* texture, const image<float>& img) {
  texture->colorb  = {};
  texture->colorf  = {};
  texture->scalarb = {};
  texture->scalarf = img;
}

// Add shape
void set_points(scene_shape* shape, const vector<int>& points) {
  shape->points = points;
}
void set_lines(scene_shape* shape, const vector<vec2i>& lines) {
  shape->lines = lines;
}
void set_triangles(scene_shape* shape, const vector<vec3i>& triangles) {
  shape->triangles = triangles;
}
void set_quads(scene_shape* shape, const vector<vec4i>& quads) {
  shape->quads = quads;
}
void set_positions(scene_shape* shape, const vector<vec3f>& positions) {
  shape->positions = positions;
}
void set_normals(scene_shape* shape, const vector<vec3f>& normals) {
  shape->normals = normals;
}
void set_texcoords(scene_shape* shape, const vector<vec2f>& texcoords) {
  shape->texcoords = texcoords;
}
void set_colors(scene_shape* shape, const vector<vec3f>& colors) {
  shape->colors = colors;
}
void set_radius(scene_shape* shape, const vector<float>& radius) {
  shape->radius = radius;
}
void set_tangents(scene_shape* shape, const vector<vec4f>& tangents) {
  shape->tangents = tangents;
}

// Add object
void set_frame(scene_object* object, const frame3f& frame) {
  object->frame = frame;
}
void set_shape(scene_object* object, scene_shape* shape) {
  object->shape = shape;
}
void set_material(scene_object* object, scene_material* material) {
  object->material = material;
}

// Add material
void set_emission(scene_material* material, const vec3f& emission,
    scene_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    scene_material* material, const vec3f& color, scene_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    scene_material* material, float specular, scene_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_metallic(
    scene_material* material, float metallic, scene_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_ior(scene_material* material, float ior) { material->ior = ior; }
void set_transmission(scene_material* material, float transmission, bool thin,
    float trdepth, scene_texture* transmission_tex) {
  material->transmission     = transmission;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->transmission_tex = transmission_tex;
}
void set_translucency(scene_material* material, float translucency, bool thin,
    float trdepth, scene_texture* translucency_tex) {
  material->translucency     = translucency;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->translucency_tex = translucency_tex;
}
void set_thin(scene_material* material, bool thin) { material->thin = thin; }
void set_roughness(
    scene_material* material, float roughness, scene_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    scene_material* material, float opacity, scene_texture* opacity_tex) {
  material->opacity     = opacity;
  material->opacity_tex = opacity_tex;
}
void set_scattering(scene_material* material, const vec3f& scattering,
    float scanisotropy, scene_texture* scattering_tex) {
  material->scattering     = scattering;
  material->scanisotropy   = scanisotropy;
  material->scattering_tex = scattering_tex;
}
void set_normalmap(scene_material* material, scene_texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// Add environment
void set_frame(scene_environment* environment, const frame3f& frame) {
  environment->frame = frame;
}
void set_emission(scene_environment* environment, const vec3f& emission,
    scene_texture* emission_tex) {
  environment->emission     = emission;
  environment->emission_tex = emission_tex;
}

// Add missing cameras.
void add_cameras(scene_model* scene) {
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
void add_radius(scene_model* scene, float radius) {
  for (auto shape : scene->shapes) {
    if (shape->points.empty() && shape->lines.empty()) continue;
    if (!shape->radius.empty()) continue;
    shape->radius.assign(shape->positions.size(), radius);
  }
}

// Add missing materials.
void add_materials(scene_model* scene) {
  auto default_material = (scene_material*)nullptr;
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
void add_sky(scene_model* scene, float sun_angle) {
  auto texture = add_texture(scene, "sky");
  auto sunsky  = image<vec4f>{{1024, 512}};
  make_sunsky(sunsky, sunsky.size(), sun_angle);
  texture->colorf.resize(sunsky.size());
  for (auto j = 0; j < sunsky.size().y; j++)
    for (auto i = 0; i < sunsky.size().x; i++)
      texture->colorf[{i, j}] = xyz(sunsky[{i, j}]);
  auto environment          = add_environment(scene, "sky");
  environment->emission     = {1, 1, 1};
  environment->emission_tex = texture;
}

// get named camera or default if camera is empty
scene_camera* get_camera(const scene_model* scene, const string& name) {
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
  for (auto camera : scene->cameras) {
    if (camera->name == "camera1") return camera;
  }
  return scene->cameras.front();
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_model* scene) {
  auto shape_bbox = unordered_map<scene_shape*, bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto shape : scene->shapes) {
    auto sbvh = invalidb3f;
    for (auto p : shape->positions) sbvh = merge(sbvh, p);
    shape_bbox[shape] = sbvh;
  }
  for (auto object : scene->objects) {
    auto sbvh = shape_bbox[object->shape];
    bbox      = merge(bbox, transform_bbox(object->frame, sbvh));
  }
  return bbox;
}

// Clone a scene
void clone_scene(scene_model* dest, const scene_model* scene) {
  dest->name      = scene->name;
  dest->copyright = scene->copyright;
  throw std::runtime_error("not implemented yet");
}

// Reduce memory usage
void trim_memory(scene_model* scene) {
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
    shape->quadspos.shrink_to_fit();
    shape->quadsnorm.shrink_to_fit();
    shape->quadstexcoord.shrink_to_fit();
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
static vec2i texture_size(const scene_texture* texture) {
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
static vec3f lookup_texture(
    const scene_texture* texture, const vec2i& ij, bool ldr_as_linear = false) {
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
static vec3f eval_texture(const scene_texture* texture, const vec2f& uv,
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

void tesselate_shape(scene_shape* shape) {
  if (shape->subdivisions) {
    if (!shape->points.empty()) {
      throw std::runtime_error("cannot subdivide points");
    } else if (!shape->lines.empty()) {
      std::tie(std::ignore, shape->texcoords) = subdivide_lines(
          shape->lines, shape->texcoords, shape->subdivisions);
      std::tie(std::ignore, shape->normals) = subdivide_lines(
          shape->lines, shape->normals, shape->subdivisions);
      std::tie(std::ignore, shape->colors) = subdivide_lines(
          shape->lines, shape->colors, shape->subdivisions);
      std::tie(std::ignore, shape->radius) = subdivide_lines(
          shape->lines, shape->radius, shape->subdivisions);
      std::tie(shape->lines, shape->positions) = subdivide_lines(
          shape->lines, shape->positions, shape->subdivisions);
    } else if (!shape->triangles.empty()) {
      std::tie(std::ignore, shape->texcoords) = subdivide_triangles(
          shape->triangles, shape->texcoords, shape->subdivisions);
      std::tie(std::ignore, shape->normals) = subdivide_triangles(
          shape->triangles, shape->normals, shape->subdivisions);
      std::tie(std::ignore, shape->colors) = subdivide_triangles(
          shape->triangles, shape->colors, shape->subdivisions);
      std::tie(std::ignore, shape->radius) = subdivide_triangles(
          shape->triangles, shape->radius, shape->subdivisions);
      std::tie(shape->triangles, shape->positions) = subdivide_triangles(
          shape->triangles, shape->positions, shape->subdivisions);
    } else if (!shape->quads.empty()) {
      if (shape->catmullclark) {
        std::tie(std::ignore, shape->texcoords) = subdivide_catmullclark(
            shape->quads, shape->texcoords, shape->subdivisions, true);
        std::tie(std::ignore, shape->normals) = subdivide_catmullclark(
            shape->quads, shape->normals, shape->subdivisions, true);
        std::tie(std::ignore, shape->colors) = subdivide_catmullclark(
            shape->quads, shape->colors, shape->subdivisions);
        std::tie(std::ignore, shape->radius) = subdivide_catmullclark(
            shape->quads, shape->radius, shape->subdivisions);
        std::tie(std::ignore, shape->positions) = subdivide_catmullclark(
            shape->quads, shape->positions, shape->subdivisions);
      } else {
        std::tie(std::ignore, shape->texcoords) = subdivide_quads(
            shape->quads, shape->texcoords, shape->subdivisions);
        std::tie(std::ignore, shape->normals) = subdivide_quads(
            shape->quads, shape->normals, shape->subdivisions);
        std::tie(std::ignore, shape->colors) = subdivide_quads(
            shape->quads, shape->colors, shape->subdivisions);
        std::tie(std::ignore, shape->radius) = subdivide_quads(
            shape->quads, shape->radius, shape->subdivisions);
        std::tie(shape->quads, shape->positions) = subdivide_quads(
            shape->quads, shape->positions, shape->subdivisions);
      }
    } else if (!shape->quadspos.empty()) {
      if (shape->catmullclark) {
        std::tie(shape->quadstexcoord, shape->texcoords) =
            subdivide_catmullclark(shape->quadstexcoord, shape->texcoords,
                shape->subdivisions, true);
        std::tie(shape->quadsnorm, shape->normals) = subdivide_catmullclark(
            shape->quadsnorm, shape->normals, shape->subdivisions, true);
        std::tie(shape->quadspos, shape->positions) = subdivide_catmullclark(
            shape->quadspos, shape->positions, shape->subdivisions);
      } else {
        std::tie(shape->quadstexcoord, shape->texcoords) = subdivide_quads(
            shape->quadstexcoord, shape->texcoords, shape->subdivisions);
        std::tie(shape->quadsnorm, shape->normals) = subdivide_quads(
            shape->quadsnorm, shape->normals, shape->subdivisions);
        std::tie(shape->quadspos, shape->positions) = subdivide_quads(
            shape->quadspos, shape->positions, shape->subdivisions);
      }
      if (shape->smooth) {
        shape->normals   = compute_normals(shape->quadspos, shape->positions);
        shape->quadsnorm = shape->quadspos;
      } else {
        shape->normals   = {};
        shape->quadsnorm = {};
      }
    } else {
      throw std::runtime_error("not supported yet");
    }
    shape->subdivisions = 0;
  }

  if (shape->displacement && shape->displacement_tex) {
    if (shape->texcoords.empty())
      throw std::runtime_error("missing texture coordinates");

    if (!shape->triangles.empty() || !shape->quads.empty()) {
      auto no_normals = shape->normals.empty();
      if (shape->normals.empty())
        shape->normals = !shape->triangles.empty()
                             ? compute_normals(
                                   shape->triangles, shape->positions)
                             : compute_normals(shape->quads, shape->positions);
      for (auto idx = 0; idx < shape->positions.size(); idx++) {
        auto disp = mean(
            eval_texture(shape->displacement_tex, shape->texcoords[idx], true));
        if (!shape->displacement_tex->scalarb.empty() ||
            !shape->displacement_tex->colorb.empty())
          disp -= 0.5f;
        shape->positions[idx] += shape->normals[idx] * shape->displacement *
                                 disp;
      }
      if (shape->smooth) {
        shape->normals = !shape->triangles.empty()
                             ? compute_normals(
                                   shape->triangles, shape->positions)
                             : compute_normals(shape->quads, shape->positions);
      } else if (no_normals) {
        shape->normals = {};
      }
    } else if (!shape->quadspos.empty()) {
      // facevarying case
      auto offset = vector<float>(shape->positions.size(), 0);
      auto count  = vector<int>(shape->positions.size(), 0);
      for (auto fid = 0; fid < shape->quadspos.size(); fid++) {
        auto qpos = shape->quadspos[fid];
        auto qtxt = shape->quadstexcoord[fid];
        for (auto i = 0; i < 4; i++) {
          auto disp = mean(eval_texture(
              shape->displacement_tex, shape->texcoords[qtxt[i]], true));
          if (!shape->displacement_tex->scalarb.empty() ||
              !shape->displacement_tex->colorb.empty())
            disp -= 0.5f;
          offset[qpos[i]] += shape->displacement * disp;
          count[qpos[i]] += 1;
        }
      }
      auto normals = compute_normals(shape->quadspos, shape->positions);
      for (auto vid = 0; vid < shape->positions.size(); vid++) {
        shape->positions[vid] += normals[vid] * offset[vid] / count[vid];
      }
      if (shape->smooth || !shape->normals.empty()) {
        shape->quadsnorm = shape->quadspos;
        shape->normals   = compute_normals(shape->quadspos, shape->positions);
      }
    }

    shape->displacement     = 0;
    shape->displacement_tex = nullptr;
  }
  if (!shape->quadspos.empty()) {
    std::tie(shape->quads, shape->positions, shape->normals, shape->texcoords) =
        split_facevarying(shape->quadspos, shape->quadsnorm,
            shape->quadstexcoord, shape->positions, shape->normals,
            shape->texcoords);
    shape->points    = {};
    shape->lines     = {};
    shape->triangles = {};
    shape->colors    = {};
    shape->radius    = {};
  }
}  // namespace yocto

void tesselate_shapes(scene_model* scene, progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, (int)scene->shapes.size()};

  // tesselate shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
    tesselate_shape(shape);
  }

  // done
  if (progress_cb) progress_cb("tesselate shape", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

#ifdef YOCTO_EMBREE
// Get Embree device
static atomic<ssize_t> embree_memory = 0;
static RTCDevice       embree_device() {
  static RTCDevice device = nullptr;
  if (!device) {
    device = rtcNewDevice("");
    rtcSetDeviceErrorFunction(
        device,
        [](void* ctx, RTCError code, const char* str) {
          switch (code) {
            case RTC_ERROR_UNKNOWN:
              throw std::runtime_error("RTC_ERROR_UNKNOWN: "s + str);
            case RTC_ERROR_INVALID_ARGUMENT:
              throw std::runtime_error("RTC_ERROR_INVALID_ARGUMENT: "s + str);
            case RTC_ERROR_INVALID_OPERATION:
              throw std::runtime_error("RTC_ERROR_INVALID_OPERATION: "s + str);
            case RTC_ERROR_OUT_OF_MEMORY:
              throw std::runtime_error("RTC_ERROR_OUT_OF_MEMORY: "s + str);
            case RTC_ERROR_UNSUPPORTED_CPU:
              throw std::runtime_error("RTC_ERROR_UNSUPPORTED_CPU: "s + str);
            case RTC_ERROR_CANCELLED:
              throw std::runtime_error("RTC_ERROR_CANCELLED: "s + str);
            default: throw std::runtime_error("invalid error code");
          }
        },
        nullptr);
    rtcSetDeviceMemoryMonitorFunction(
        device,
        [](void* userPtr, ssize_t bytes, bool post) {
          embree_memory += bytes;
          return true;
        },
        nullptr);
  }
  return device;
}

// Initialize Embree BVH
static void init_embree_bvh(
    scene_shape* shape, const scene_bvh_params& params) {
  auto edevice = embree_device();
  if (shape->embree_bvh) rtcReleaseScene(shape->embree_bvh);
  shape->embree_bvh = rtcNewScene(edevice);
  auto escene       = shape->embree_bvh;
  if (params.bvh == scene_bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == scene_bvh_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  if (!shape->points.empty()) {
    throw std::runtime_error("embree does not support points");
  } else if (!shape->lines.empty()) {
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto& l : shape->lines) {
      if (last_index == l.x) {
        elines.push_back((int)epositions.size() - 1);
        epositions.push_back({shape->positions[l.y], shape->radius[l.y]});
      } else {
        elines.push_back((int)epositions.size());
        epositions.push_back({shape->positions[l.x], shape->radius[l.x]});
        epositions.push_back({shape->positions[l.y], shape->radius[l.y]});
      }
      last_index = l.y;
    }
    auto egeometry = rtcNewGeometry(
        edevice, RTC_GEOMETRY_TYPE_FLAT_LINEAR_CURVE);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, 4 * 4, epositions.size());
    auto embree_lines     = rtcSetNewGeometryBuffer(
        egeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, 4, elines.size());
    memcpy(embree_positions, epositions.data(), epositions.size() * 16);
    memcpy(embree_lines, elines.data(), elines.size() * 4);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else if (!shape->triangles.empty()) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    if (params.bvh == scene_bvh_type::embree_compact) {
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape->positions.data(), 0, 3 * 4,
          shape->positions.size());
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT3, shape->triangles.data(), 0, 3 * 4,
          shape->triangles.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape->positions.size());
      auto embree_triangles = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
          shape->triangles.size());
      memcpy(embree_positions, shape->positions.data(),
          shape->positions.size() * 12);
      memcpy(embree_triangles, shape->triangles.data(),
          shape->triangles.size() * 12);
    }
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else if (!shape->quads.empty()) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    if (params.bvh == scene_bvh_type::embree_compact) {
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape->positions.data(), 0, 3 * 4,
          shape->positions.size());
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT4, shape->quads.data(), 0, 4 * 4, shape->quads.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape->positions.size());
      auto embree_quads     = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4,
          shape->quads.size());
      memcpy(embree_positions, shape->positions.data(),
          shape->positions.size() * 12);
      memcpy(embree_quads, shape->quads.data(), shape->quads.size() * 16);
    }
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else {
    throw std::runtime_error("empty shapes not supported");
  }
  rtcCommitScene(escene);
}

static void init_embree_bvh(
    scene_model* scene, const scene_bvh_params& params) {
  // scene bvh
  auto edevice = embree_device();
  if (scene->embree_bvh) rtcReleaseScene(scene->embree_bvh);
  scene->embree_bvh = rtcNewScene(edevice);
  auto escene       = scene->embree_bvh;
  if (params.bvh == scene_bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == scene_bvh_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  auto object_id = 0;
  for (auto object : scene->objects) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_INSTANCE);
    rtcSetGeometryInstancedScene(egeometry, object->shape->embree_bvh);
    rtcSetGeometryTransform(
        egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &object->frame);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, object_id++);
  }
  rtcCommitScene(escene);
}

static void update_embree_bvh(scene_model* scene,
    const vector<scene_object*>&           updated_objects,
    const vector<scene_shape*>&            updated_shapes,
    const scene_bvh_params&                params) {
  throw std::runtime_error("not implemented yet");
  // // scene bvh
  // auto escene = scene->embree_bvh;
  // for (auto& [object_id, instance_id] : scene->embree_instances) {
  //   auto object    = scene->objects[object_id];
  //   auto frame     = scene->objects[instance_id]->frame * object->frame;
  //   auto egeometry = rtcGetGeometry(escene, instance_id);
  //   rtcSetGeometryInstancedScene(egeometry, object->shape->embree_bvh);
  //   rtcSetGeometryTransform(
  //       egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &frame);
  //   rtcCommitGeometry(egeometry);
  // }
  // rtcCommitScene(escene);
}

static bool intersect_shape_embree_bvh(scene_shape* shape, const ray3f& ray,
    int& element, vec2f& uv, float& distance, bool find_any) {
  RTCRayHit embree_ray;
  embree_ray.ray.org_x     = ray.o.x;
  embree_ray.ray.org_y     = ray.o.y;
  embree_ray.ray.org_z     = ray.o.z;
  embree_ray.ray.dir_x     = ray.d.x;
  embree_ray.ray.dir_y     = ray.d.y;
  embree_ray.ray.dir_z     = ray.d.z;
  embree_ray.ray.tnear     = ray.tmin;
  embree_ray.ray.tfar      = ray.tmax;
  embree_ray.ray.flags     = 0;
  embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
  embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  RTCIntersectContext embree_ctx;
  rtcInitIntersectContext(&embree_ctx);
  rtcIntersect1(shape->embree_bvh, &embree_ctx, &embree_ray);
  if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}

static bool intersect_scene_embree_bvh(const scene_model* scene,
    const ray3f& ray, int& object, int& element, vec2f& uv, float& distance,
    bool find_any) {
  RTCRayHit embree_ray;
  embree_ray.ray.org_x     = ray.o.x;
  embree_ray.ray.org_y     = ray.o.y;
  embree_ray.ray.org_z     = ray.o.z;
  embree_ray.ray.dir_x     = ray.d.x;
  embree_ray.ray.dir_y     = ray.d.y;
  embree_ray.ray.dir_z     = ray.d.z;
  embree_ray.ray.tnear     = ray.tmin;
  embree_ray.ray.tfar      = ray.tmax;
  embree_ray.ray.flags     = 0;
  embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
  embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  RTCIntersectContext embree_ctx;
  rtcInitIntersectContext(&embree_ctx);
  rtcIntersect1(scene->embree_bvh, &embree_ctx, &embree_ray);
  if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
  object   = (int)embree_ray.hit.instID[0];
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}
#endif

// primitive used to sort bvh entries
struct scene_bvh_primitive {
  bbox3f bbox      = invalidb3f;
  vec3f  center    = zero3f;
  int    primitive = 0;
};

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
static std::pair<int, int> split_sah(
    vector<scene_bvh_primitive>& primitives, int start, int end) {
  // initialize split axis and position
  auto split_axis = 0;
  auto mid        = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++) cbbox = merge(cbbox, primitives[i].center);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, split_axis};

  // consider N bins, compute their cost and keep the minimum
  const int nbins    = 16;
  auto      middle   = 0.0f;
  auto      min_cost = flt_max;
  auto      area     = [](auto& b) {
    auto size = b.max - b.min;
    return 1e-12f + 2 * size.x * size.y + 2 * size.x * size.z +
           2 * size.y * size.z;
  };
  for (auto saxis = 0; saxis < 3; saxis++) {
    for (auto b = 1; b < nbins; b++) {
      auto split     = cbbox.min[saxis] + b * csize[saxis] / nbins;
      auto left_bbox = invalidb3f, right_bbox = invalidb3f;
      auto left_nprims = 0, right_nprims = 0;
      for (auto i = start; i < end; i++) {
        if (primitives[i].center[saxis] < split) {
          left_bbox = merge(left_bbox, primitives[i].bbox);
          left_nprims += 1;
        } else {
          right_bbox = merge(right_bbox, primitives[i].bbox);
          right_nprims += 1;
        }
      }
      auto cost = 1 + left_nprims * area(left_bbox) / area(cbbox) +
                  right_nprims * area(right_bbox) / area(cbbox);
      if (cost < min_cost) {
        min_cost   = cost;
        middle     = split;
        split_axis = saxis;
      }
    }
  }
  // split
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [split_axis, middle](auto& primitive) {
                    return primitive.center[split_axis] < middle;
                  }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    throw std::runtime_error("bad bvh split");
    split_axis = 0;
    mid        = (start + end) / 2;
  }

  return {mid, split_axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and
// axis.
static std::pair<int, int> split_balanced(
    vector<scene_bvh_primitive>& primitives, int start, int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++) cbbox = merge(cbbox, primitives[i].center);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // balanced tree split: find the largest axis of the
  // bounding box and split along this one right in the middle
  mid = (start + end) / 2;
  std::nth_element(primitives.data() + start, primitives.data() + mid,
      primitives.data() + end, [axis](auto& primitive_a, auto& primitive_b) {
        return primitive_a.center[axis] < primitive_b.center[axis];
      });

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    // throw std::runtime_error("bad bvh split");
    mid = (start + end) / 2;
  }

  return {mid, axis};
}

// Splits a BVH node using the middle heuristic. Returns split position and
// axis.
static std::pair<int, int> split_middle(
    vector<scene_bvh_primitive>& primitives, int start, int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++) cbbox = merge(cbbox, primitives[i].center);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // split the space in the middle along the largest axis
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [axis, middle = center(cbbox)[axis]](auto& primitive) {
                    return primitive.center[axis] < middle;
                  }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    // throw std::runtime_error("bad bvh split");
    mid = (start + end) / 2;
  }

  return {mid, axis};
}

// Split bvh nodes according to a type
static std::pair<int, int> split_nodes(vector<scene_bvh_primitive>& primitives,
    int start, int end, scene_bvh_type type) {
  switch (type) {
    case scene_bvh_type::default_: return split_middle(primitives, start, end);
    case scene_bvh_type::highquality: return split_sah(primitives, start, end);
    case scene_bvh_type::middle: return split_middle(primitives, start, end);
    case scene_bvh_type::balanced:
      return split_balanced(primitives, start, end);
    default: throw std::runtime_error("should not have gotten here");
  }
}

// Maximum number of primitives per BVH node.
const int scene_bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh_serial(vector<scene_bvh_node>& nodes,
    vector<scene_bvh_primitive>& primitives, scene_bvh_type type) {
  // prepare to build nodes
  nodes.clear();
  nodes.reserve(primitives.size() * 2);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // create nodes until the queue is empty
  while (!queue.empty()) {
    // grab node to work on
    auto next = queue.front();
    queue.pop_front();
    auto nodeid = next.x, start = next.y, end = next.z;

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox = merge(node.bbox, primitives[i].bbox);

    // split into two children
    if (end - start > scene_bvh_max_prims) {
      // get split
      auto [mid, axis] = split_nodes(primitives, start, end, type);

      // make an internal node
      node.internal = true;
      node.axis     = axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = end - start;
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

#if 0

// Build BVH nodes
static void build_bvh_parallel(
    const shared_ptr<bvh_tree>& bvh, vector<bbox3f>& bboxes, bvh_type type) {
  // get values
  auto& nodes      = bvh->nodes;
  auto& primitives = bvh->primitives;

  // prepare to build nodes
  nodes.clear();
  nodes.reserve(bboxes.size() * 2);

  // prepare primitives
  bvh->primitives.resize(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) bvh->primitives[idx] = idx;

  // prepare centers
  auto centers = vector<vec3f>(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++)
    centers[idx] = center(bboxes[idx]);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // synchronization
  atomic<int>          num_processed_prims(0);
  std::mutex                queue_mutex;
  vector<future<void>> futures;
  auto                      nthreads = std::thread::hardware_concurrency();

  // create nodes until the queue is empty
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&nodes, &primitives, &bboxes, &centers, &type,
                                &num_processed_prims, &queue_mutex, &queue] {
          while (true) {
            // exit if needed
            if (num_processed_prims >= primitives.size()) return;

            // grab node to work on
            auto next = zero3i;
            {
              std::lock_guard<std::mutex> lock{queue_mutex};
              if (!queue.empty()) {
                next = queue.front();
                queue.pop_front();
              }
            }

            // wait a bit if needed
            if (next == zero3i) {
              std::this_thread::sleep_for(std::chrono::microseconds(10));
              continue;
            }

            // grab node
            auto  nodeid = next.x, start = next.y, end = next.z;
            auto& node = nodes[nodeid];

            // compute bounds
            node.bbox = invalidb3f;
            for (auto i = start; i < end; i++)
              node.bbox = merge(node.bbox, bboxes[primitives[i]]);

            // split into two children
            if (end - start > bvh_max_prims) {
              // get split
              auto [mid, axis] = split_nodes(
                  primitives, bboxes, centers, start, end, type);

              // make an internal node
              {
                std::lock_guard<std::mutex> lock{queue_mutex};
                node.internal = true;
                node.axis     = axis;
                node.num      = 2;
                node.start    = (int)nodes.size();
                nodes.emplace_back();
                nodes.emplace_back();
                queue.push_back({node.start + 0, start, mid});
                queue.push_back({node.start + 1, mid, end});
              }
            } else {
              // Make a leaf node
              node.internal = false;
              node.num      = end - start;
              node.start    = start;
              num_processed_prims += node.num;
            }
          }
        }));
  }
  for (auto& f : futures) f.get();

  // cleanup
  nodes.shrink_to_fit();
}

#endif

// Update bvh
static void update_bvh(scene_bvh* bvh, const vector<bbox3f>& bboxes) {
  for (auto nodeid = (int)bvh->nodes.size() - 1; nodeid >= 0; nodeid--) {
    auto& node = bvh->nodes[nodeid];
    node.bbox  = invalidb3f;
    if (node.internal) {
      for (auto idx = 0; idx < 2; idx++) {
        node.bbox = merge(node.bbox, bvh->nodes[node.start + idx].bbox);
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        node.bbox = merge(node.bbox, bboxes[node.start + idx]);
      }
    }
  }
}

static void init_bvh(scene_shape* shape, const scene_bvh_params& params) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (params.bvh == scene_bvh_type::embree_default ||
      params.bvh == scene_bvh_type::embree_highquality ||
      params.bvh == scene_bvh_type::embree_compact) {
    return init_embree_bvh(shape, params);
  }
#endif

  // build primitives
  auto primitives = vector<scene_bvh_primitive>{};
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < shape->points.size(); idx++) {
      auto& p             = shape->points[idx];
      auto& primitive     = primitives.emplace_back();
      primitive.bbox      = point_bounds(shape->positions[p], shape->radius[p]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  } else if (!shape->lines.empty()) {
    for (auto idx = 0; idx < shape->lines.size(); idx++) {
      auto& l         = shape->lines[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  } else if (!shape->triangles.empty()) {
    for (auto idx = 0; idx < shape->triangles.size(); idx++) {
      auto& primitive = primitives.emplace_back();
      auto& t         = shape->triangles[idx];
      primitive.bbox  = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  } else if (!shape->quads.empty()) {
    for (auto idx = 0; idx < shape->quads.size(); idx++) {
      auto& q         = shape->quads[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = idx;
    }
  }

  // build nodes
  if (shape->bvh) delete shape->bvh;
  shape->bvh = new scene_bvh{};
  build_bvh_serial(shape->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  shape->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    shape->bvh->primitives.push_back(primitive.primitive);
  }
}

void init_bvh(scene_model* scene, const scene_bvh_params& params,
    progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1 + (int)scene->shapes.size()};

  // shapes
  for (auto idx = 0; idx < scene->shapes.size(); idx++) {
    if (progress_cb) progress_cb("build shape bvh", progress.x++, progress.y);
    init_bvh(scene->shapes[idx], params);
  }

  // embree
#ifdef YOCTO_EMBREE
  if (params.bvh == scene_bvh_type::embree_default ||
      params.bvh == scene_bvh_type::embree_highquality ||
      params.bvh == scene_bvh_type::embree_compact) {
    return init_embree_bvh(scene, params);
  }
#endif

  // handle progress
  if (progress_cb) progress_cb("build scene bvh", progress.x++, progress.y);

  // instance bboxes
  auto primitives            = vector<scene_bvh_primitive>{};
  auto object_id             = 0;
  auto empty_instance_frames = vector<frame3f>{identity3x4f};
  for (auto object : scene->objects) {
    auto& primitive = primitives.emplace_back();
    primitive.bbox =
        object->shape->bvh->nodes.empty()
            ? invalidb3f
            : transform_bbox(object->frame, object->shape->bvh->nodes[0].bbox);
    primitive.center    = center(primitive.bbox);
    primitive.primitive = object_id++;
  }

  // build nodes
  if (scene->bvh) delete scene->bvh;
  scene->bvh = new scene_bvh{};
  build_bvh_serial(scene->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  scene->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    scene->bvh->primitives.push_back(primitive.primitive);
  }

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

static void update_bvh(scene_shape* shape, const scene_bvh_params& params) {
#ifdef YOCTO_EMBREE
  if (shape->embree_bvh) {
    throw std::runtime_error("embree shape update not implemented");
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>(shape->bvh->primitives.size());
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& p     = shape->points[shape->bvh->primitives[idx]];
      bboxes[idx] = point_bounds(shape->positions[p], shape->radius[p]);
    }
  } else if (!shape->lines.empty()) {
    bboxes = vector<bbox3f>(shape->lines.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& l     = shape->lines[shape->bvh->primitives[idx]];
      bboxes[idx] = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
    }
  } else if (!shape->triangles.empty()) {
    bboxes = vector<bbox3f>(shape->triangles.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& t     = shape->triangles[shape->bvh->primitives[idx]];
      bboxes[idx] = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
    }
  } else if (!shape->quads.empty()) {
    bboxes = vector<bbox3f>(shape->quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape->quads[shape->bvh->primitives[idx]];
      bboxes[idx] = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
    }
  }

  // update nodes
  update_bvh(shape->bvh, bboxes);
}

void update_bvh(scene_model*     scene,
    const vector<scene_object*>& updated_objects,
    const vector<scene_shape*>&  updated_shapes,
    const scene_bvh_params&      params) {
  for (auto shape : updated_shapes) update_bvh(shape, params);

#ifdef YOCTO_EMBREE
  if (scene->embree_bvh) {
    update_embree_bvh(scene, updated_objects, updated_shapes, params);
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>(scene->bvh->primitives.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto instance = scene->bvh->primitives[idx];
    auto object   = scene->objects[instance];
    auto sbvh     = object->shape->bvh;
    bboxes[idx]   = transform_bbox(object->frame, sbvh->nodes[0].bbox);
  }

  // update nodes
  update_bvh(scene->bvh, bboxes);
}

// Intersect ray with a bvh->
static bool intersect_shape_bvh(scene_shape* shape, const ray3f& ray_,
    int& element, vec2f& uv, float& distance, bool find_any) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (shape->embree_bvh) {
    return intersect_shape_embree_bvh(
        shape, ray_, element, uv, distance, find_any);
  }
#endif

  // get bvh and shape pointers for fast access
  auto bvh = shape->bvh;

  // check empty
  if (bvh->nodes.empty()) return false;

  // node stack
  int  node_stack[128];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto hit = false;

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh->nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else if (!shape->points.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& p = shape->points[shape->bvh->primitives[idx]];
        if (intersect_point(
                ray, shape->positions[p], shape->radius[p], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape->lines.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& l = shape->lines[shape->bvh->primitives[idx]];
        if (intersect_line(ray, shape->positions[l.x], shape->positions[l.y],
                shape->radius[l.x], shape->radius[l.y], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape->triangles.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& t = shape->triangles[shape->bvh->primitives[idx]];
        if (intersect_triangle(ray, shape->positions[t.x],
                shape->positions[t.y], shape->positions[t.z], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape->quads.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape->quads[shape->bvh->primitives[idx]];
        if (intersect_quad(ray, shape->positions[q.x], shape->positions[q.y],
                shape->positions[q.z], shape->positions[q.w], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh->
static bool intersect_scene_bvh(const scene_model* scene, const ray3f& ray_,
    int& object, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (scene->embree_bvh) {
    return intersect_scene_embree_bvh(
        scene, ray_, object, element, uv, distance, find_any);
  }
#endif

  // get bvh and scene pointers for fast access
  auto bvh = scene->bvh;

  // check empty
  if (bvh->nodes.empty()) return false;

  // node stack
  int  node_stack[128];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // shared variables
  auto hit = false;

  // copy ray to modify it
  auto ray = ray_;

  // prepare ray for fast queries
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh->nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto object_ = scene->objects[scene->bvh->primitives[idx]];
        auto inv_ray = transform_ray(
            inverse(object_->frame, non_rigid_frames), ray);
        if (intersect_shape_bvh(
                object_->shape, inv_ray, element, uv, distance, find_any)) {
          hit      = true;
          object   = scene->bvh->primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh->
static bool intersect_instance_bvh(const scene_object* object, const ray3f& ray,
    int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto inv_ray = transform_ray(inverse(object->frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      object->shape, inv_ray, element, uv, distance, find_any);
}

scene_intersection intersect_scene_bvh(const scene_model* scene,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = scene_intersection{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.object,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}
scene_intersection intersect_instance_bvh(const scene_object* object,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = scene_intersection{};
  intersection.hit  = intersect_instance_bvh(object, ray, intersection.element,
      intersection.uv, intersection.distance, find_any, non_rigid_frames);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox(scene_model* scene) {
  scene->name                = "cornellbox";
  auto camera                = add_camera(scene);
  camera->frame              = frame3f{{0, 1, 3.9}};
  camera->lens               = 0.035;
  camera->aperture           = 0.0;
  camera->focus              = 3.9;
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
