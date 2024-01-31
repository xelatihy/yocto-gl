//
// Implementation for Yocto/Scene.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modeling.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::unique_ptr;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CAMERA PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Image resolution from camera
vec2i camera_resolution(const camera_data& camera, int resolution) {
  return (camera.aspect >= 1)
             ? vec2i{resolution, (int)round(resolution / camera.aspect)}
             : vec2i{(int)round(resolution * camera.aspect), resolution};
}

// Generates a ray from a camera for yimg::image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera(const camera_data& camera, vec2f image_uv, vec2f lens_uv) {
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// pixel access
vec4f lookup_texture(
    const texture_data& texture, vec2i ij, bool ldr_as_linear) {
  if (!texture.pixelsf.empty()) return texture.pixelsf[ij];
  if (!texture.pixelsb.empty())
    return ldr_as_linear ? byte_to_float(texture.pixelsb[ij])
                         : srgb_to_rgb(byte_to_float(texture.pixelsb[ij]));
  return vec4f{0, 0, 0, 0};
}

// Evaluates an image at a point `uv`.
vec4f eval_texture(const texture_data& texture, vec2f uv, bool ldr_as_linear) {
  if (texture.pixelsf.empty() && texture.pixelsb.empty()) return {0, 0, 0, 0};

  // get texture width/height
  auto size = max(texture.pixelsf.size(), texture.pixelsb.size());

  // get coordinates normalized for tiling
  auto st = (texture.clamp ? clamp(uv, 0.0f, 1.0f) : mod(uv, 1.0f)) *
            (vec2f)size;

  // handle interpolation
  if (texture.nearest) {
    auto ij = clamp((vec2i)st, {0, 0}, size - 1);
    return lookup_texture(texture, ij, ldr_as_linear);
  } else {
    auto ij   = clamp((vec2i)st, zero2i, size - 1);
    auto i1j  = (ij + vec2i{1, 0}) % size;
    auto ij1  = (ij + vec2i{0, 1}) % size;
    auto i1j1 = (ij + vec2i{1, 1}) % size;
    auto w    = st - (vec2f)ij;
    return lookup_texture(texture, ij, ldr_as_linear) * (1 - w.x) * (1 - w.y) +
           lookup_texture(texture, ij1, ldr_as_linear) * (1 - w.x) * w.y +
           lookup_texture(texture, i1j, ldr_as_linear) * w.x * (1 - w.y) +
           lookup_texture(texture, i1j1, ldr_as_linear) * w.x * w.y;
  }
}

// Helpers
vec4f eval_texture(
    const scene_data& scene, int texture, vec2f uv, bool ldr_as_linear) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(scene.textures[texture], uv, ldr_as_linear);
}

// conversion from image
texture_data image_to_texture(const image<vec4f>& image, bool linear) {
  auto texture = texture_data{};
  if (linear) {
    texture.pixelsf = image;
  } else {
    texture.pixelsb = float_to_byte(image);
  }
  return texture;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATERIAL PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// constant values
static const auto min_roughness = 0.03f * 0.03f;

// Evaluate material
material_point eval_material(const scene_data& scene,
    const material_data& material, vec2f texcoord, vec4f color_shp) {
  // evaluate textures
  auto emission_tex  = eval_texture(scene, material.emission_tex, texcoord);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, true);
  auto scattering_tex = eval_texture(scene, material.scattering_tex, texcoord);

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
  if (material.type == material_type::refractive ||
      material.type == material_type::volumetric ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::gltfpbr ||
      point.type == material_type::glossy) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// check if a material is a delta or volumetric
bool is_delta(const material_data& material) {
  return (material.type == material_type::reflective &&
             material.roughness == 0) ||
         (material.type == material_type::refractive &&
             material.roughness == 0) ||
         (material.type == material_type::transparent &&
             material.roughness == 0) ||
         (material.type == material_type::volumetric);
}
bool is_volumetric(const material_data& material) {
  return material.type == material_type::refractive ||
         material.type == material_type::volumetric ||
         material.type == material_type::subsurface;
}

// check if a brdf is a delta
bool is_delta(const material_point& material) {
  return (material.type == material_type::reflective &&
             material.roughness == 0) ||
         (material.type == material_type::refractive &&
             material.roughness == 0) ||
         (material.type == material_type::transparent &&
             material.roughness == 0) ||
         (material.type == material_type::volumetric);
}
bool has_volume(const material_point& material) {
  return material.type == material_type::refractive ||
         material.type == material_type::volumetric ||
         material.type == material_type::subsurface;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INSTANCE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Eval position
vec3f eval_position(const scene_data& scene, const instance_data& instance,
    int element, vec2f uv) {
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
    return {0, 0, 0};
  }
}

// Shape element normal.
vec3f eval_element_normal(
    const scene_data& scene, const instance_data& instance, int element) {
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
vec3f eval_normal(const scene_data& scene, const instance_data& instance,
    int element, vec2f uv) {
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
    return {0, 0, 0};
  }
}

// Eval texcoord
vec2f eval_texcoord(const scene_data& scene, const instance_data& instance,
    int element, vec2f uv) {
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
    return vec2f{0, 0};
  }
}

#if 0
// Shape element normal.
static pair<vec3f, vec3f> eval_tangents(
    const trace_shape& shape, int element, vec2f uv) {
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
    return {{0,0,0}, {0,0,0}};
  }
}
#endif

// Shape element normal.
pair<vec3f, vec3f> eval_element_tangents(
    const scene_data& scene, const instance_data& instance, int element) {
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

vec3f eval_normalmap(const scene_data& scene, const instance_data& instance,
    int element, vec2f uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalidid &&
      (!shape.triangles.empty() || !shape.quads.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, true));
    auto [tu, tv]    = eval_element_tangents(scene, instance, element);
    auto frame       = frame3f{tu, tv, normal, {0, 0, 0}};
    frame.x          = orthonormalize(frame.x, frame.z);
    frame.y          = normalize(cross(frame.z, frame.x));
    auto flip_v      = dot(frame.y, tv) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading position
vec3f eval_shading_position(const scene_data& scene,
    const instance_data& instance, int element, vec2f uv, vec3f outgoing) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    return eval_position(scene, instance, element, uv);
  } else if (!shape.lines.empty()) {
    return eval_position(scene, instance, element, uv);
  } else if (!shape.points.empty()) {
    return eval_position(shape, element, uv);
  } else {
    return {0, 0, 0};
  }
}

// Eval shading normal
vec3f eval_shading_normal(const scene_data& scene,
    const instance_data& instance, int element, vec2f uv, vec3f outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalidid) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::refractive) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else if (!shape.lines.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    return orthonormalize(outgoing, normal);
  } else if (!shape.points.empty()) {
    return outgoing;
  } else {
    return {0, 0, 0};
  }
}

// Eval color
vec4f eval_color(const scene_data& scene, const instance_data& instance,
    int element, vec2f uv) {
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

// Evaluate material
material_point eval_material(const scene_data& scene,
    const instance_data& instance, int element, vec2f uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex  = eval_texture(scene, material.emission_tex, texcoord);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, true);
  auto scattering_tex = eval_texture(scene, material.scattering_tex, texcoord);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex) * xyz(color_shp);
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
  if (material.type == material_type::refractive ||
      material.type == material_type::volumetric ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::gltfpbr ||
      point.type == material_type::glossy) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  } else if (material.type == material_type::volumetric) {
    point.roughness = 0;
  } else {
    if (point.roughness < min_roughness) point.roughness = 0;
  }

  return point;
}

// check if an instance is volumetric
bool is_volumetric(const scene_data& scene, const instance_data& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENVIRONMENT PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate environment color.
vec3f eval_environment(const scene_data& scene,
    const environment_data& environment, vec3f direction) {
  auto wl       = transform_direction(inverse(environment.frame), direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_data& scene, vec3f direction) {
  auto emission = vec3f{0, 0, 0};
  for (auto& environment : scene.environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Add missing cameras.
void add_camera(scene_data& scene) {
  scene.camera_names.emplace_back("camera");
  auto& camera        = scene.cameras.emplace_back();
  camera.orthographic = false;
  camera.film         = 0.036f;
  camera.aspect       = (float)16 / (float)9;
  camera.aperture     = 0;
  camera.lens         = 0.050f;
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

// Add a sky environment
void add_sky(scene_data& scene, float sun_angle) {
  scene.texture_names.emplace_back("sky");
  auto& texture = scene.textures.emplace_back();
  texture       = image_to_texture(make_sunsky({1024, 512}, sun_angle), true);
  scene.environment_names.emplace_back("sky");
  auto& environment        = scene.environments.emplace_back();
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
}

// get named camera or default if camera is empty
int find_camera(const scene_data& scene, const string& name) {
  if (scene.cameras.empty()) return invalidid;
  if (scene.camera_names.empty()) return 0;
  for (auto idx : range(scene.camera_names.size())) {
    if (scene.camera_names[idx] == name) return (int)idx;
  }
  for (auto idx : range(scene.camera_names.size())) {
    if (scene.camera_names[idx] == "default") return (int)idx;
  }
  for (auto idx : range(scene.camera_names.size())) {
    if (scene.camera_names[idx] == "camera") return (int)idx;
  }
  for (auto idx : range(scene.camera_names.size())) {
    if (scene.camera_names[idx] == "camera0") return (int)idx;
  }
  for (auto idx : range(scene.camera_names.size())) {
    if (scene.camera_names[idx] == "camera1") return (int)idx;
  }
  return 0;
}

// check if it has lights
bool has_lights(const scene_data& scene) {
  for (auto& environment : scene.environments) {
    if (environment.emission != vec3f{0, 0, 0}) return true;
  }
  for (auto& instance : scene.instances) {
    auto& shape = scene.shapes[instance.shape];
    if (shape.triangles.empty() || shape.quads.empty()) continue;
    auto& material = scene.materials[instance.material];
    if (material.emission != vec3f{0, 0, 0}) return true;
  }
  return false;
}

// create a scene from a shape
scene_data make_shape_scene(const shape_data& shape, bool addsky) {
  // scene
  auto scene = scene_data{};
  // shape
  scene.shape_names.emplace_back("shape");
  scene.shapes.push_back(shape);
  // material
  scene.material_names.emplace_back("material");
  auto& shape_material     = scene.materials.emplace_back();
  shape_material.type      = material_type::glossy;
  shape_material.color     = {0.5f, 1.0f, 0.5f};
  shape_material.roughness = 0.2f;
  // instance
  scene.instance_names.emplace_back("instance");
  auto& shape_instance    = scene.instances.emplace_back();
  shape_instance.shape    = 0;
  shape_instance.material = 0;
  // camera
  add_camera(scene);
  // environment
  if (addsky) add_sky(scene);
  // done
  return scene;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_data& scene) {
  auto shape_bbox = vector<bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto& shape : scene.shapes) {
    auto& sbvh = shape_bbox.emplace_back();
    for (auto p : shape.positions) sbvh = merge(sbvh, p);
  }
  for (auto& instance : scene.instances) {
    auto& sbvh = shape_bbox[instance.shape];
    bbox       = merge(bbox, transform_bbox(instance.frame, sbvh));
  }
  return bbox;
}

// Scene creation helpers
scene_data make_scene() { return scene_data{}; }

// Scene creation helpers
int add_camera(
    scene_data& scene, const string& name, const camera_data& camera) {
  scene.camera_names.push_back(name);
  scene.cameras.push_back(camera);
  return (int)scene.cameras.size() - 1;
}
int add_texture(
    scene_data& scene, const string& name, const texture_data& texture) {
  scene.texture_names.push_back(name);
  scene.textures.push_back(texture);
  return (int)scene.textures.size() - 1;
}
int add_material(
    scene_data& scene, const string& name, const material_data& material) {
  scene.material_names.push_back(name);
  scene.materials.push_back(material);
  return (int)scene.materials.size() - 1;
}
int add_shape(scene_data& scene, const string& name, const shape_data& shape) {
  scene.shape_names.push_back(name);
  scene.shapes.push_back(shape);
  return (int)scene.shapes.size() - 1;
}
int add_instance(
    scene_data& scene, const string& name, const instance_data& instance) {
  scene.instance_names.push_back(name);
  scene.instances.push_back(instance);
  return (int)scene.instances.size() - 1;
}
int add_environment(scene_data& scene, const string& name,
    const environment_data& environment) {
  scene.environment_names.push_back(name);
  scene.environments.push_back(environment);
  return (int)scene.environments.size() - 1;
}

int add_camera(scene_data& scene, const string& name, vec3f from, vec3f to,
    float lens, float aspect, float aperture, float focus_offset) {
  return add_camera(scene, name,
      {.frame       = lookat_frame(from, to, {0, 1, 0}),
          .lens     = lens,
          .aspect   = aspect,
          .aperture = aperture,
          .focus    = length(from - to) + focus_offset});
}

int add_camera(scene_data& scene, const string& name, const frame3f& frame,
    float lens, float aspect, float aperture, float focus) {
  return add_camera(scene, name,
      {.frame       = frame,
          .lens     = lens,
          .aspect   = aspect,
          .aperture = aperture,
          .focus    = focus});
}

// Scene creation helpers
int add_material(scene_data& scene, const string& name, vec3f emission,
    material_type type, vec3f color, float roughness, int color_tex,
    int roughness_tex, int normal_tex) {
  return add_material(scene, name,
      {.type             = type,
          .emission      = emission,
          .color         = color,
          .roughness     = roughness,
          .color_tex     = color_tex,
          .roughness_tex = roughness_tex,
          .normal_tex    = normal_tex});
}

// Scene creation helpers
int add_material(scene_data& scene, const string& name, material_type type,
    vec3f color, float roughness, int color_tex, int roughness_tex,
    int normal_tex) {
  return add_material(scene, name,
      {.type             = type,
          .color         = color,
          .roughness     = roughness,
          .color_tex     = color_tex,
          .roughness_tex = roughness_tex,
          .normal_tex    = normal_tex});
}

// Scene creation helpers
int add_material(scene_data& scene, const string& name, material_type type,
    vec3f color, float roughness, vec3f scattering, float scanisotropy,
    float trdepth, int color_tex, int roughness_tex, int scattering_tex,
    int normal_tex) {
  return add_material(scene, name,
      {.type              = type,
          .color          = color,
          .roughness      = roughness,
          .scattering     = scattering,
          .scanisotropy   = scanisotropy,
          .trdepth        = trdepth,
          .color_tex      = color_tex,
          .roughness_tex  = roughness_tex,
          .scattering_tex = scattering_tex,
          .normal_tex     = normal_tex});
}

// Scene creation helpers
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    int shape, int material) {
  return add_instance(
      scene, name, {.frame = frame, .shape = shape, .material = material});
}

// Scene creation helpers
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    const shape_data& shape, const material_data& material) {
  return add_instance(scene, name,
      {.frame       = frame,
          .shape    = add_shape(scene, name, shape),
          .material = add_material(scene, name, material)});
}

// Scene creation helpers
int add_environment(scene_data& scene, const string& name, const frame3f& frame,
    vec3f emission, int emission_tex) {
  return add_environment(scene, name,
      {.frame = frame, .emission = emission, .emission_tex = emission_tex});
  return (int)scene.environments.size() - 1;
}

// Scene creation helpers
int add_texture(scene_data& scene, const string& name,
    const image<vec4f>& texture, bool linear) {
  if (texture.empty()) return invalidid;
  if (linear) {
    scene.textures.push_back({.pixelsf = texture});
  } else {
    scene.textures.push_back({.pixelsb = float_to_byte(texture)});
  }
  return (int)scene.textures.size() - 1;
}

// Scene creation helpers
int add_texture(scene_data& scene, const string& name,
    const image<vec4b>& texture, bool linear) {
  if (texture.empty()) return invalidid;
  if (linear) {
    scene.textures.push_back({.pixelsf = srgbb_to_rgb(texture)});
  } else {
    scene.textures.push_back({.pixelsb = texture});
  }
  return (int)scene.textures.size() - 1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE TESSELATION
// -----------------------------------------------------------------------------
namespace yocto {

void tesselate_subdiv(
    shape_data& shape, subdiv_data& subdiv_, const scene_data& scene) {
  auto subdiv = subdiv_;

  if (subdiv.subdivisions > 0) {
    if (subdiv.catmullclark) {
      for ([[maybe_unused]] auto subdivision : range(subdiv.subdivisions)) {
        std::tie(subdiv.quadstexcoord, subdiv.texcoords) =
            subdivide_catmullclark(
                subdiv.quadstexcoord, subdiv.texcoords, true);
        std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_catmullclark(
            subdiv.quadsnorm, subdiv.normals, true);
        std::tie(subdiv.quadspos, subdiv.positions) = subdivide_catmullclark(
            subdiv.quadspos, subdiv.positions);
      }
    } else {
      for ([[maybe_unused]] auto subdivision : range(subdiv.subdivisions)) {
        std::tie(subdiv.quadstexcoord, subdiv.texcoords) = subdivide_quads(
            subdiv.quadstexcoord, subdiv.texcoords);
        std::tie(subdiv.quadsnorm, subdiv.normals) = subdivide_quads(
            subdiv.quadsnorm, subdiv.normals);
        std::tie(subdiv.quadspos, subdiv.positions) = subdivide_quads(
            subdiv.quadspos, subdiv.positions);
      }
    }
    if (subdiv.smooth) {
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
      subdiv.quadsnorm = subdiv.quadspos;
    } else {
      subdiv.normals   = {};
      subdiv.quadsnorm = {};
    }
  }

  if (subdiv.displacement != 0 && subdiv.displacement_tex != invalidid) {
    if (subdiv.texcoords.empty())
      throw std::runtime_error("missing texture coordinates");

    // facevarying case
    auto offset = vector<float>(subdiv.positions.size(), 0);
    auto count  = vector<int>(subdiv.positions.size(), 0);
    for (auto fid = (size_t)0; fid < subdiv.quadspos.size(); fid++) {
      auto qpos = subdiv.quadspos[fid];
      auto qtxt = subdiv.quadstexcoord[fid];
      for (auto i : range(4)) {
        auto& displacement_tex = scene.textures[subdiv.displacement_tex];
        auto  disp             = mean(
            eval_texture(displacement_tex, subdiv.texcoords[qtxt[i]], false));
        if (!displacement_tex.pixelsb.empty()) disp -= 0.5f;
        offset[qpos[i]] += subdiv.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = quads_normals(subdiv.quadspos, subdiv.positions);
    for (auto vid = (size_t)0; vid < subdiv.positions.size(); vid++) {
      subdiv.positions[vid] += normals[vid] * offset[vid] / (float)count[vid];
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      subdiv.quadsnorm = subdiv.quadspos;
      subdiv.normals   = quads_normals(subdiv.quadspos, subdiv.positions);
    }
  }

  shape = fvshape_to_shape(
      {subdiv.quadspos, subdiv.quadsnorm, subdiv.quadstexcoord,
          subdiv.positions, subdiv.normals, subdiv.texcoords});
}

void tesselate_subdivs(scene_data& scene) {
  // tesselate shapes
  for (auto& subdiv : scene.subdivs) {
    tesselate_subdiv(scene.shapes[subdiv.shape], subdiv, scene);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

scene_data make_cornellbox() {
  auto scene = scene_data{};

  auto& camera    = scene.cameras.emplace_back();
  camera.frame    = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 3.9f}};
  camera.lens     = 0.035f;
  camera.aperture = 0.0;
  camera.focus    = 3.9f;
  camera.film     = 0.024f;
  camera.aspect   = 1;

  auto& floor_shape       = scene.shapes.emplace_back();
  floor_shape.positions   = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& floor_material    = scene.materials.emplace_back();
  floor_material.color    = {0.725f, 0.71f, 0.68f};
  auto& floor_instance    = scene.instances.emplace_back();
  floor_instance.shape    = (int)scene.shapes.size() - 1;
  floor_instance.material = (int)scene.materials.size() - 1;

  auto& ceiling_shape       = scene.shapes.emplace_back();
  ceiling_shape.positions   = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& ceiling_material    = scene.materials.emplace_back();
  ceiling_material.color    = {0.725f, 0.71f, 0.68f};
  auto& ceiling_instance    = scene.instances.emplace_back();
  ceiling_instance.shape    = (int)scene.shapes.size() - 1;
  ceiling_instance.material = (int)scene.materials.size() - 1;

  auto& backwall_shape     = scene.shapes.emplace_back();
  backwall_shape.positions = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& backwall_material  = scene.materials.emplace_back();
  backwall_material.color  = {0.725f, 0.71f, 0.68f};
  auto& backwall_instance  = scene.instances.emplace_back();
  backwall_instance.shape  = (int)scene.shapes.size() - 1;
  backwall_instance.material = (int)scene.materials.size() - 1;

  auto& rightwall_shape       = scene.shapes.emplace_back();
  rightwall_shape.positions   = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& rightwall_material    = scene.materials.emplace_back();
  rightwall_material.color    = {0.14f, 0.45f, 0.091f};
  auto& rightwall_instance    = scene.instances.emplace_back();
  rightwall_instance.shape    = (int)scene.shapes.size() - 1;
  rightwall_instance.material = (int)scene.materials.size() - 1;

  auto& leftwall_shape     = scene.shapes.emplace_back();
  leftwall_shape.positions = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& leftwall_material  = scene.materials.emplace_back();
  leftwall_material.color  = {0.63f, 0.065f, 0.05f};
  auto& leftwall_instance  = scene.instances.emplace_back();
  leftwall_instance.shape  = (int)scene.shapes.size() - 1;
  leftwall_instance.material = (int)scene.materials.size() - 1;

  auto& shortbox_shape       = scene.shapes.emplace_back();
  shortbox_shape.positions   = {{0.53f, 0.6f, 0.75f}, {0.7f, 0.6f, 0.17f},
        {0.13f, 0.6f, 0.0f}, {-0.05f, 0.6f, 0.57f}, {-0.05f, 0.0f, 0.57f},
        {-0.05f, 0.6f, 0.57f}, {0.13f, 0.6f, 0.0f}, {0.13f, 0.0f, 0.0f},
        {0.53f, 0.0f, 0.75f}, {0.53f, 0.6f, 0.75f}, {-0.05f, 0.6f, 0.57f},
        {-0.05f, 0.0f, 0.57f}, {0.7f, 0.0f, 0.17f}, {0.7f, 0.6f, 0.17f},
        {0.53f, 0.6f, 0.75f}, {0.53f, 0.0f, 0.75f}, {0.13f, 0.0f, 0.0f},
        {0.13f, 0.6f, 0.0f}, {0.7f, 0.6f, 0.17f}, {0.7f, 0.0f, 0.17f},
        {0.53f, 0.0f, 0.75f}, {0.7f, 0.0f, 0.17f}, {0.13f, 0.0f, 0.0f},
        {-0.05f, 0.0f, 0.57f}};
  shortbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
        {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
        {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& shortbox_material    = scene.materials.emplace_back();
  shortbox_material.color    = {0.725f, 0.71f, 0.68f};
  auto& shortbox_instance    = scene.instances.emplace_back();
  shortbox_instance.shape    = (int)scene.shapes.size() - 1;
  shortbox_instance.material = (int)scene.materials.size() - 1;

  auto& tallbox_shape       = scene.shapes.emplace_back();
  tallbox_shape.positions   = {{-0.53f, 1.2f, 0.09f}, {0.04f, 1.2f, -0.09f},
        {-0.14f, 1.2f, -0.67f}, {-0.71f, 1.2f, -0.49f}, {-0.53f, 0.0f, 0.09f},
        {-0.53f, 1.2f, 0.09f}, {-0.71f, 1.2f, -0.49f}, {-0.71f, 0.0f, -0.49f},
        {-0.71f, 0.0f, -0.49f}, {-0.71f, 1.2f, -0.49f}, {-0.14f, 1.2f, -0.67f},
        {-0.14f, 0.0f, -0.67f}, {-0.14f, 0.0f, -0.67f}, {-0.14f, 1.2f, -0.67f},
        {0.04f, 1.2f, -0.09f}, {0.04f, 0.0f, -0.09f}, {0.04f, 0.0f, -0.09f},
        {0.04f, 1.2f, -0.09f}, {-0.53f, 1.2f, 0.09f}, {-0.53f, 0.0f, 0.09f},
        {-0.53f, 0.0f, 0.09f}, {0.04f, 0.0f, -0.09f}, {-0.14f, 0.0f, -0.67f},
        {-0.71f, 0.0f, -0.49f}};
  tallbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
        {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
        {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& tallbox_material    = scene.materials.emplace_back();
  tallbox_material.color    = {0.725f, 0.71f, 0.68f};
  auto& tallbox_instance    = scene.instances.emplace_back();
  tallbox_instance.shape    = (int)scene.shapes.size() - 1;
  tallbox_instance.material = (int)scene.materials.size() - 1;

  auto& light_shape       = scene.shapes.emplace_back();
  light_shape.positions   = {{-0.25f, 1.99f, 0.25f}, {-0.25f, 1.99f, -0.25f},
        {0.25f, 1.99f, -0.25f}, {0.25f, 1.99f, 0.25f}};
  light_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& light_material    = scene.materials.emplace_back();
  light_material.emission = {17, 12, 4};
  auto& light_instance    = scene.instances.emplace_back();
  light_instance.shape    = (int)scene.shapes.size() - 1;
  light_instance.material = (int)scene.materials.size() - 1;

  return scene;
}

}  // namespace yocto
