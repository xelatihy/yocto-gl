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
#include "yocto_noise.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CAMERA PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera for image plane coordinate image_uv and
// the lens coordinates lens_uv.
ray3f eval_camera(
    const camera_data& camera, const vec2f& image_uv, const vec2f& lens_uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  if (!camera.orthographic) {
    auto uv = flip_u(image_uv);
    auto q  = vec3f{film * (uv - 0.5f), camera.lens};
    // ray direction through the lens center
    auto dc = -normalize(q);
    // point on the lens
    auto e = vec3f{lens_uv * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = dc * camera.focus / abs(dc.z);
    // correct ray direction to account for camera focusing
    auto d = normalize(p - e);
    // done
    return ray3f{
        transform_point(camera.frame, e), transform_direction(camera.frame, d)};
  } else {
    auto uv    = flip_u(image_uv);
    auto scale = 1 / camera.lens;
    auto q     = vec3f{film * (uv - 0.5f) * scale, camera.lens};
    // point on the lens
    auto e = vec3f{-xy(q), 0} + vec3f{lens_uv * camera.aperture / 2, 0};
    // point on the focus plane
    auto p = vec3f{-xy(q), -camera.focus};
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
    const texture_data& texture, const vec2s& ij, bool ldr_as_linear) {
  if (!texture.pixelsf.empty()) return texture.pixelsf[ij];
  if (!texture.pixelsb.empty())
    return ldr_as_linear ? byte_to_float(texture.pixelsb[ij])
                         : srgb_to_rgb(byte_to_float(texture.pixelsb[ij]));
  return vec4f{0, 0, 0, 0};
}

// Evaluates an image at a point `uv`.
vec4f eval_texture(const texture_data& texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (texture.pixelsf.empty() && texture.pixelsb.empty()) return {0, 0, 0, 0};

  // get texture width/height
  auto size = max(texture.pixelsf.extents(), texture.pixelsb.extents());

  // get coordinates normalized for tiling
  auto st = (clamp_to_edge ? clamp(uv, 0, 1) : mod(uv, 1)) * size;

  // handle interpolation
  if (no_interpolation) {
    auto ij = clamp((vec2s)st, 0, size - 1);
    return lookup_texture(texture, ij, ldr_as_linear);
  } else {
    auto ij   = clamp((vec2s)st, 0, size - 1);
    auto i1j  = (ij + vec2s{1, 0}) % size;
    auto ij1  = (ij + vec2s{0, 1}) % size;
    auto i1j1 = (ij + vec2s{1, 1}) % size;
    auto w    = st - ij;
    return lookup_texture(texture, ij, ldr_as_linear) * (1 - w.x) * (1 - w.y) +
           lookup_texture(texture, ij1, ldr_as_linear) * (1 - w.x) * w.y +
           lookup_texture(texture, i1j, ldr_as_linear) * w.x * (1 - w.y) +
           lookup_texture(texture, i1j1, ldr_as_linear) * w.x * w.y;
  }
}
vec4f eval_texture(
    const texture_data& texture, const vec2f& uv, bool ldr_as_linear) {
  return eval_texture(
      texture, uv, ldr_as_linear, texture.nearest, texture.clamp);
}

// Helpers
vec4f eval_texture(
    const scene_data& scene, int texture, const vec2f& uv, bool ldr_as_linear) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(scene.textures[texture], uv, ldr_as_linear,
      scene.textures[texture].nearest, scene.textures[texture].clamp);
}
vec4f eval_texture(const scene_data& scene, int texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(scene.textures[texture], uv, ldr_as_linear,
      no_interpolation, clamp_to_edge);
}

// conversion from image
texture_data image_to_texture(const array2d<vec4f>& image, bool linear) {
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
    const material_data& material, const vec2f& texcoord,
    const vec4f& shape_color) {
  // evaluate textures
  auto emission_tex   = eval_texture(scene, material.emission_tex, texcoord);
  auto color_tex      = eval_texture(scene, material.color_tex, texcoord);
  auto roughness_tex  = eval_texture(scene, material.roughness_tex, texcoord);
  auto scattering_tex = eval_texture(scene, material.scattering_tex, texcoord);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex) * xyz(shape_color);
  point.color        = material.color * xyz(color_tex) * xyz(shape_color);
  point.opacity      = material.opacity * alpha(color_tex) * alpha(shape_color);
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
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (is_triangles(shape)) {
    auto triangle = shape.triangles[element];
    return transform_point(
        instance.frame, interpolate_triangle(shape.positions, triangle, uv));
  } else if (is_quads(shape)) {
    auto quad = shape.quads[element];
    return transform_point(
        instance.frame, interpolate_quad(shape.positions, quad, uv));
  } else if (is_lines(shape)) {
    auto line = shape.lines[element];
    return transform_point(
        instance.frame, interpolate_line(shape.positions, line, uv));
  } else if (is_points(shape)) {
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
  if (is_triangles(shape)) {
    auto triangle = shape.triangles[element];
    return transform_normal(
        instance.frame, triangle_normal(shape.positions, triangle));
  } else if (is_quads(shape)) {
    auto quad = shape.quads[element];
    return transform_normal(instance.frame, quad_normal(shape.positions, quad));
  } else if (is_lines(shape)) {
    auto line = shape.lines[element];
    return transform_normal(
        instance.frame, line_tangent(shape.positions, line));
  } else if (is_points(shape)) {
    return {0, 0, 1};
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
vec3f eval_normal(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.normals.empty())
    return eval_element_normal(scene, instance, element);
  if (is_triangles(shape)) {
    auto triangle = shape.triangles[element];
    return transform_normal(instance.frame,
        normalize(interpolate_triangle(shape.normals, triangle, uv)));
  } else if (is_quads(shape)) {
    auto quad = shape.quads[element];
    return transform_normal(
        instance.frame, normalize(interpolate_quad(shape.normals, quad, uv)));
  } else if (is_lines(shape)) {
    auto line = shape.lines[element];
    return transform_normal(
        instance.frame, normalize(interpolate_line(shape.normals, line, uv)));
  } else if (is_points(shape)) {
    return transform_normal(
        instance.frame, normalize(shape.normals[shape.points[element]]));
  } else {
    return {0, 0, 0};
  }
}

// Eval texcoord
vec2f eval_texcoord(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.texcoords.empty()) return uv;
  if (is_triangles(shape)) {
    auto triangle = shape.triangles[element];
    return interpolate_triangle(shape.texcoords, triangle, uv);
  } else if (is_quads(shape)) {
    auto quad = shape.quads[element];
    return interpolate_quad(shape.texcoords, quad, uv);
  } else if (is_lines(shape)) {
    auto line = shape.lines[element];
    return interpolate_line(shape.texcoords, line, uv);
  } else if (is_points(shape)) {
    return shape.texcoords[shape.points[element]];
  } else {
    return vec2f{0, 0};
  }
}

// Shape element normal.
pair<vec3f, vec3f> eval_element_tangents(
    const scene_data& scene, const instance_data& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (is_triangles(shape) && !shape.texcoords.empty()) {
    auto triangle = shape.triangles[element];
    auto [tu, tv] = triangle_tangents_fromuv(
        shape.positions, shape.texcoords, triangle);
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else if (is_quads(shape) && !shape.texcoords.empty()) {
    auto quad     = shape.quads[element];
    auto [tu, tv] = quad_tangents_fromuv(
        shape.positions, shape.texcoords, quad, {0, 0});
    return {transform_direction(instance.frame, tu),
        transform_direction(instance.frame, tv)};
  } else {
    return {};
  }
}

vec3f eval_normalmap(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalidid &&
      (is_triangles(shape) || is_quads(shape))) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, true));
    auto [tu, tv]    = eval_element_tangents(scene, instance, element);
    auto frame       = orthonormalize(frame3f{tu, tv, normal, {0, 0, 0}});
    auto flip_v      = dot(frame[1], tv) < 0;
    normalmap[1] *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading position
vec3f eval_shading_position(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape = scene.shapes[instance.shape];
  if (is_triangles(shape) || is_quads(shape)) {
    return eval_position(scene, instance, element, uv);
  } else if (is_lines(shape)) {
    return eval_position(scene, instance, element, uv);
  } else if (is_points(shape)) {
    return eval_position(shape, element, uv);
  } else {
    return {0, 0, 0};
  }
}

// Eval shading normal
vec3f eval_shading_normal(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (is_triangles(shape) || is_quads(shape)) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalidid) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::refractive) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else if (is_lines(shape)) {
    auto normal = eval_normal(scene, instance, element, uv);
    return orthonormalize(outgoing, normal);
  } else if (is_points(shape)) {
    return outgoing;
  } else {
    return {0, 0, 0};
  }
}

// Eval color
vec4f eval_color(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (is_triangles(shape)) {
    auto triangle = shape.triangles[element];
    return interpolate_triangle(shape.colors, triangle, uv);
  } else if (is_quads(shape)) {
    auto quad = shape.quads[element];
    return interpolate_quad(shape.colors, quad, uv);
  } else if (is_lines(shape)) {
    auto line = shape.lines[element];
    return interpolate_line(shape.colors, line, uv);
  } else if (is_points(shape)) {
    return shape.colors[shape.points[element]];
  } else {
    return {0, 0, 0, 0};
  }
}

// Evaluate material
material_point eval_material(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex   = eval_texture(scene, material.emission_tex, texcoord);
  auto color_shp      = eval_color(scene, instance, element, uv);
  auto color_tex      = eval_texture(scene, material.color_tex, texcoord);
  auto roughness_tex  = eval_texture(scene, material.roughness_tex, texcoord);
  auto scattering_tex = eval_texture(scene, material.scattering_tex, texcoord);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex) * xyz(color_shp);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * alpha(color_tex) * alpha(color_shp);
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
    const environment_data& environment, const vec3f& direction) {
  auto texcoord = cartesiany_to_sphericaluv(
      transform_direction(inverse(environment.frame), direction));
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_data& scene, const vec3f& direction) {
  auto emission = vec3f{0, 0, 0};
  for (auto& environment : scene.environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Make a camera
camera_data make_camera(const frame3f& frame, float lens, float aspect,
    float aperture, float focus) {
  return {.frame = frame,
      .lens      = lens,
      .aspect    = aspect,
      .aperture  = aperture,
      .focus     = focus};
}
camera_data make_camera(const vec3f& from, const vec3f& to, const vec3f& up,
    float lens, float aspect, float aperture) {
  return {.frame = lookat_frame(from, to, up),
      .lens      = lens,
      .aspect    = aspect,
      .aperture  = aperture,
      .focus     = distance(from, to)};
}

// Make a texture
texture_data make_texture(
    const array2d<vec4f>& image, bool as_float, bool srgb) {
  return {
      .pixelsf = as_float ? (srgb ? srgb_to_rgb(image) : image)
                          : array2d<vec4f>{},
      .pixelsb = as_float ? array2d<vec4b>{}
                          : (srgb ? float_to_byte(image) : rgb_to_srgbb(image)),
  };
}

// Make a material of different types
material_data make_emissive_material(
    const vec3f& emission, int emission_tex, int normal_tex) {
  return {.type     = material_type::matte,
      .emission     = emission,
      .color        = {0, 0, 0},
      .emission_tex = emission_tex,
      .normal_tex   = normal_tex};
}
material_data make_matte_material(
    const vec3f& color, int color_tex, int normal_tex) {
  return {.type   = material_type::matte,
      .color      = color,
      .color_tex  = color_tex,
      .normal_tex = normal_tex};
}
material_data make_glossy_material(const vec3f& color, float roughness,
    int color_tex, int roughness_tex, int normal_tex) {
  return {.type      = material_type::glossy,
      .color         = color,
      .roughness     = roughness,
      .color_tex     = color_tex,
      .roughness_tex = roughness_tex,
      .normal_tex    = normal_tex};
}
material_data make_reflective_material(const vec3f& color, float roughness,
    int color_tex, int roughness_tex, int normal_tex) {
  return {.type      = material_type::reflective,
      .color         = color,
      .roughness     = roughness,
      .color_tex     = color_tex,
      .roughness_tex = roughness_tex,
      .normal_tex    = normal_tex};
}
material_data make_transparent_material(const vec3f& color, float roughness,
    int color_tex, int roughness_tex, int normal_tex) {
  return {.type      = material_type::transparent,
      .color         = color,
      .roughness     = roughness,
      .color_tex     = color_tex,
      .roughness_tex = roughness_tex,
      .normal_tex    = normal_tex};
}
material_data make_refractive_material(const vec3f& color, float roughness,
    int color_tex, int roughness_tex, int normal_tex, float ior) {
  return {.type      = material_type::refractive,
      .color         = color,
      .roughness     = roughness,
      .ior           = ior,
      .color_tex     = color_tex,
      .roughness_tex = roughness_tex,
      .normal_tex    = normal_tex};
}
material_data make_scattering_material(const vec3f& color,
    const vec3f& scattering, float roughness, int color_tex, int scattering_tex,
    int roughness_tex, int normal_tex, float ior, float scanisotropy,
    float trdepth) {
  return {.type       = material_type::refractive,
      .color          = color,
      .roughness      = roughness,
      .ior            = ior,
      .scattering     = scattering,
      .scanisotropy   = scanisotropy,
      .trdepth        = trdepth,
      .color_tex      = color_tex,
      .roughness_tex  = roughness_tex,
      .scattering_tex = scattering_tex,
      .normal_tex     = normal_tex};
}
material_data make_volumetric_material(const vec3f& color,
    const vec3f& scattering, float scanisotropy, float trdepth) {
  return {.type     = material_type::volumetric,
      .color        = color,
      .scattering   = scattering,
      .scanisotropy = scanisotropy,
      .trdepth      = trdepth};
}

// Make an instance
instance_data make_instance(const frame3f& frame, int shape, int material) {
  return {.frame = frame, .shape = shape, .material = material};
}

// Make an environment
environment_data make_environment(
    const frame3f& frame, const vec3f& emission, int emission_tex) {
  return {.frame = frame, .emission = emission, .emission_tex = emission_tex};
}

// Make a subdiv
subdiv_data make_subdiv(const shape_data& subdiv, int shape, int subdivisions,
    float displacement, int displacement_tex) {
  return {.quadspos     = subdiv.quads,
      .positions        = subdiv.positions,
      .shape            = shape,
      .subdivisions     = subdivisions,
      .displacement     = displacement,
      .displacement_tex = displacement_tex};
}
subdiv_data make_subdiv(const fvshape_data& subdiv, int shape, int subdivisions,
    float displacement, int displacement_tex) {
  return {.quadspos     = subdiv.quadspos,
      .quadsnorm        = subdiv.quadsnorm,
      .quadstexcoord    = subdiv.quadstexcoord,
      .positions        = subdiv.positions,
      .normals          = subdiv.normals,
      .texcoords        = subdiv.texcoords,
      .shape            = shape,
      .subdivisions     = subdivisions,
      .displacement     = displacement,
      .displacement_tex = displacement_tex};
}

// Add a camera
int add_camera(
    scene_data& scene, const string& name, const camera_data& camera) {
  scene.camera_names.push_back(name);
  scene.cameras.push_back(camera);
  return (int)scene.cameras.size() - 1;
}
// Add a shape
int add_shape(scene_data& scene, const string& name, const shape_data& shape) {
  scene.shape_names.push_back(name);
  scene.shapes.push_back(shape);
  return (int)scene.shapes.size() - 1;
}
// Add a texture
int add_texture(
    scene_data& scene, const string& name, const texture_data& texture) {
  scene.texture_names.push_back(name);
  scene.textures.push_back(texture);
  return (int)scene.textures.size() - 1;
}
// Add a material
int add_material(
    scene_data& scene, const string& name, const material_data& material) {
  scene.material_names.push_back(name);
  scene.materials.push_back(material);
  return (int)scene.materials.size() - 1;
}
// Add an instance
int add_instance(
    scene_data& scene, const string& name, const instance_data& instance) {
  scene.instance_names.push_back(name);
  scene.instances.push_back(instance);
  return (int)scene.instances.size() - 1;
}
// Add an environment
int add_environment(scene_data& scene, const string& name,
    const environment_data& environment) {
  scene.environment_names.push_back(name);
  scene.environments.push_back(environment);
  return (int)scene.environments.size() - 1;
}
// Add a subdiv
int add_subdiv(
    scene_data& scene, const string& name, const subdiv_data& subdiv) {
  scene.subdiv_names.push_back(name);
  scene.subdivs.push_back(subdiv);
  return (int)scene.subdivs.size() - 1;
}

// Add a camera
int add_camera(scene_data& scene, const string& name, const frame3f& frame,
    float lens, float aspect, float aperture, float focus) {
  return add_camera(
      scene, name, make_camera(frame, lens, aspect, aperture, focus));
}
int add_camera(scene_data& scene, const string& name, const vec3f& from,
    const vec3f& to, const vec3f& up, float lens, float aspect,
    float aperture) {
  return add_camera(
      scene, name, make_camera(from, to, up, lens, aspect, aperture));
}

// Add a texture
int add_texture(scene_data& scene, const string& name,
    const array2d<vec4f>& image, bool as_float, bool srgb) {
  return add_texture(scene, name, make_texture(image, as_float, srgb));
}

// Add a material
int add_emission_material(scene_data& scene, const string& name,
    const vec3f& emission, int emission_tex, int normal_tex) {
  return add_material(
      scene, name, make_emissive_material(emission, emission_tex, normal_tex));
}
int add_matte_material(scene_data& scene, const string& name,
    const vec3f& color, int color_tex, int normal_tex) {
  return add_material(
      scene, name, make_matte_material(color, color_tex, normal_tex));
}
int add_glossy_material(scene_data& scene, const string& name,
    const vec3f& color, float roughness, int color_tex, int roughness_tex,
    int normal_tex) {
  return add_material(scene, name,
      make_glossy_material(
          color, roughness, color_tex, roughness_tex, normal_tex));
}
int add_reflective_material(scene_data& scene, const string& name,
    const vec3f& color, float roughness, int color_tex, int roughness_tex,
    int normal_tex) {
  return add_material(scene, name,
      make_reflective_material(
          color, roughness, color_tex, roughness_tex, normal_tex));
}
int add_transparent_material(scene_data& scene, const string& name,
    const vec3f& color, float roughness, int color_tex, int roughness_tex,
    int normal_tex) {
  return add_material(scene, name,
      make_transparent_material(
          color, roughness, color_tex, roughness_tex, normal_tex));
}
int add_refractive_material(scene_data& scene, const string& name,
    const vec3f& color, float roughness, int color_tex, int roughness_tex,
    int normal_tex, float ior) {
  return add_material(scene, name,
      make_refractive_material(
          color, roughness, color_tex, roughness_tex, normal_tex, ior));
}
int add_scattering_material(scene_data& scene, const string& name,
    const vec3f& color, const vec3f& scattering, float roughness, int color_tex,
    int scattering_tex, int roughness_tex, int normal_tex, float ior,
    float scanisotropy, float trdepth) {
  return add_material(scene, name,
      make_scattering_material(color, scattering, roughness, color_tex,
          scattering_tex, roughness_tex, normal_tex, ior, scanisotropy,
          trdepth));
}
int add_volumetric_material(scene_data& scene, const string& name,
    const vec3f& color, const vec3f& scattering, float scanisotropy,
    float trdepth) {
  return add_material(scene, name,
      make_volumetric_material(color, scattering, scanisotropy, trdepth));
}

// Add an instance
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    int shape, int material) {
  return add_instance(scene, name, make_instance(frame, shape, material));
}
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    const shape_data& shape, const material_data& material) {
  return add_instance(scene, name,
      make_instance(frame, add_shape(scene, name, shape),
          add_material(scene, name, material)));
}
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    const shape_data& shape, const material_data& material,
    const array2d<vec4f>& color_tex, const array2d<vec4f>& roughness_tex,
    const array2d<vec4f>& normal_tex) {
  auto material_ = material;
  if (!color_tex.empty()) {
    material_.color     = {1, 1, 1};
    material_.color_tex = add_texture(scene, name + "_diffuse", color_tex);
  }
  if (!roughness_tex.empty()) {
    material_.roughness     = 1;
    material_.metallic      = 1;
    material_.roughness_tex = add_texture(
        scene, name + "_roughness", roughness_tex);
  }
  if (!normal_tex.empty()) {
    material_.normal_tex = add_texture(scene, name + "_normal", normal_tex);
  }
  return add_instance(scene, name,
      make_instance(frame, add_shape(scene, name, shape),
          add_material(scene, name, material_)));
}
int add_instance(scene_data& scene, const string& name, const frame3f& frame,
    const shape_data& shape, const material_data& material,
    const texture_data& color_tex, const texture_data& roughness_tex,
    const texture_data& normal_tex) {
  auto material_ = material;
  if (!color_tex.pixelsf.empty() || !color_tex.pixelsb.empty())
    material_.color_tex = add_texture(scene, name + "_diffuse", color_tex);
  if (!roughness_tex.pixelsf.empty() || !roughness_tex.pixelsb.empty())
    material_.roughness_tex = add_texture(
        scene, name + "_roughness", roughness_tex);
  if (!normal_tex.pixelsf.empty() || !normal_tex.pixelsb.empty())
    material_.roughness_tex = add_texture(scene, name + "_normal", normal_tex);
  return add_instance(scene, name,
      make_instance(frame, add_shape(scene, name, shape),
          add_material(scene, name, material_)));
}

// Add an environment
int add_environment(scene_data& scene, const string& name, const frame3f& frame,
    const vec3f& emission, int emission_tex) {
  return add_environment(
      scene, name, make_environment(frame, emission, emission_tex));
}
int add_environment(scene_data& scene, const string& name, const frame3f& frame,
    const vec3f& emission, const array2d<vec4f>& emission_tex) {
  return add_environment(scene, name,
      make_environment(frame, emission,
          add_texture(scene, name, emission_tex, true, false)));
}

// Add a subdiv
int add_subdiv(scene_data& scene, const string& name, const shape_data& subdiv,
    int shape, int subdivisions, float displacement, int displacement_tex) {
  return add_subdiv(scene, name,
      make_subdiv(subdiv, shape, subdivisions, displacement, displacement_tex));
}
int add_subdiv(scene_data& scene, const string& name,
    const fvshape_data& subdiv, int shape, int subdivisions, float displacement,
    int displacement_tex) {
  return add_subdiv(scene, name,
      make_subdiv(subdiv, shape, subdivisions, displacement, displacement_tex));
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
  environment.emission     = {0.25, 0.25, 0.25};
  environment.emission_tex = (int)scene.textures.size() - 1;
}

// get named camera or default if camera is empty
int find_camera(const scene_data& scene, const string& name) {
  if (scene.cameras.empty()) return invalidid;
  if (scene.camera_names.empty()) return 0;
  for (auto&& [idx, cname] : enumerate(scene.camera_names)) {
    if (cname == name) return (int)idx;
  }
  for (auto&& [idx, cname] : enumerate(scene.camera_names)) {
    if (cname == "default") return (int)idx;
  }
  for (auto&& [idx, cname] : enumerate(scene.camera_names)) {
    if (cname == "camera") return (int)idx;
  }
  for (auto&& [idx, cname] : enumerate(scene.camera_names)) {
    if (cname == "camera0") return (int)idx;
  }
  for (auto&& [idx, cname] : enumerate(scene.camera_names)) {
    if (cname == "camera1") return (int)idx;
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
scene_data make_shape_scene(const shape_data& shape, bool sky) {
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
  if (sky) add_sky(scene);
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
            eval_texture(displacement_tex, subdiv.texcoords[qtxt[i]], true));
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

  shape = {};
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, subdiv.quadspos, subdiv.quadsnorm, subdiv.quadstexcoord,
      subdiv.positions, subdiv.normals, subdiv.texcoords);
}

void tesselate_subdivs(scene_data& scene) {
  // tesselate shapes
  for (auto& subdiv : scene.subdivs) {
    tesselate_subdiv(scene.shapes[subdiv.shape], subdiv, scene);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

size_t compute_memory(const scene_data& scene) {
  auto vector_memory = [](auto& values) -> size_t {
    if (values.empty()) return 0;
    return values.size() * sizeof(*values.data());
  };

  auto memory = (size_t)0;
  memory += vector_memory(scene.cameras);
  memory += vector_memory(scene.instances);
  memory += vector_memory(scene.materials);
  memory += vector_memory(scene.shapes);
  memory += vector_memory(scene.textures);
  memory += vector_memory(scene.environments);
  memory += vector_memory(scene.camera_names);
  memory += vector_memory(scene.instance_names);
  memory += vector_memory(scene.material_names);
  memory += vector_memory(scene.shape_names);
  memory += vector_memory(scene.texture_names);
  memory += vector_memory(scene.environment_names);
  for (auto& shape : scene.shapes) {
    memory += vector_memory(shape.points);
    memory += vector_memory(shape.lines);
    memory += vector_memory(shape.triangles);
    memory += vector_memory(shape.quads);
    memory += vector_memory(shape.positions);
    memory += vector_memory(shape.normals);
    memory += vector_memory(shape.texcoords);
    memory += vector_memory(shape.colors);
    memory += vector_memory(shape.triangles);
  }
  for (auto& subdiv : scene.subdivs) {
    memory += vector_memory(subdiv.quadspos);
    memory += vector_memory(subdiv.quadsnorm);
    memory += vector_memory(subdiv.quadstexcoord);
    memory += vector_memory(subdiv.positions);
    memory += vector_memory(subdiv.normals);
    memory += vector_memory(subdiv.texcoords);
  }
  for (auto& texture : scene.textures) {
    memory += vector_memory(texture.pixelsb);
    memory += vector_memory(texture.pixelsf);
  }
  return memory;
}

vector<string> scene_stats(const scene_data& scene, bool verbose) {
  auto accumulate = [](const auto& values, const auto& func) -> size_t {
    auto sum = (size_t)0;
    for (auto& value : values) sum += func(value);
    return sum;
  };
  auto format = [](size_t num) {
    auto str = string{};
    while (num > 0) {
      str = std::to_string(num % 1000) + (str.empty() ? "" : ",") + str;
      num /= 1000;
    }
    if (str.empty()) str = "0";
    while (str.size() < 20) str = " " + str;
    return str;
  };
  auto format3 = [](auto nums) {
    auto str = string{};
    for (auto num : nums) {
      str += (str.empty() ? "" : " ") + std::to_string(num);
    }
    while (str.size() < 48) str = " " + str;
    return str;
  };

  auto bbox = compute_bounds(scene);

  auto stats = vector<string>{};
  stats.push_back("cameras:      " + format(scene.cameras.size()));
  stats.push_back("instances:    " + format(scene.instances.size()));
  stats.push_back("materials:    " + format(scene.materials.size()));
  stats.push_back("shapes:       " + format(scene.shapes.size()));
  stats.push_back("subdivs:      " + format(scene.subdivs.size()));
  stats.push_back("environments: " + format(scene.environments.size()));
  stats.push_back("textures:     " + format(scene.textures.size()));
  stats.push_back("memory:       " + format(compute_memory(scene)));
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
                  format(accumulate(scene.subdivs,
                      [](auto& subdiv) { return subdiv.quadspos.size(); })));
  stats.push_back("texels4b:     " +
                  format(accumulate(scene.textures,
                      [](auto& texture) { return texture.pixelsb.size(); })));
  stats.push_back("texels4f:     " +
                  format(accumulate(scene.textures,
                      [](auto& texture) { return texture.pixelsf.size(); })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("diagonal:     " + format3(diagonal(bbox)));

  return stats;
}

// Checks for validity of the scene.
vector<string> scene_validation(const scene_data& scene, bool notextures) {
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
  auto check_empty_textures = [&errs](const scene_data& scene) {
    for (auto&& [idx, texture] : enumerate(scene.textures)) {
      if (texture.pixelsf.empty() && texture.pixelsb.empty()) {
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
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

scene_data make_cornellbox() {
  auto scene = scene_data{};

  add_camera(scene, "camera", {0, 1, 3.75f}, {0, 1, 0}, {0, 1, 0}, 0.050f, 1);

  add_instance(scene, "floor", identity3x4f,
      {.positions    = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.color = {0.725f, 0.71f, 0.68f}});

  add_instance(scene, "ceiling", identity3x4f,
      {.positions    = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.color = {0.725f, 0.71f, 0.68f}});

  add_instance(scene, "backwall", identity3x4f,
      {.positions    = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.color = {0.725f, 0.71f, 0.68f}});

  add_instance(scene, "rightwall", identity3x4f,
      {.positions    = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.color = {0.14f, 0.45f, 0.091f}});

  add_instance(scene, "leftwall", identity3x4f,
      {.positions    = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.color = {0.63f, 0.065f, 0.05f}});

  add_instance(scene, "shortbox", identity3x4f,
      {.positions    = {{0.53f, 0.6f, 0.75f}, {0.7f, 0.6f, 0.17f},
              {0.13f, 0.6f, 0.0f}, {-0.05f, 0.6f, 0.57f}, {-0.05f, 0.0f, 0.57f},
              {-0.05f, 0.6f, 0.57f}, {0.13f, 0.6f, 0.0f}, {0.13f, 0.0f, 0.0f},
              {0.53f, 0.0f, 0.75f}, {0.53f, 0.6f, 0.75f}, {-0.05f, 0.6f, 0.57f},
              {-0.05f, 0.0f, 0.57f}, {0.7f, 0.0f, 0.17f}, {0.7f, 0.6f, 0.17f},
              {0.53f, 0.6f, 0.75f}, {0.53f, 0.0f, 0.75f}, {0.13f, 0.0f, 0.0f},
              {0.13f, 0.6f, 0.0f}, {0.7f, 0.6f, 0.17f}, {0.7f, 0.0f, 0.17f},
              {0.53f, 0.0f, 0.75f}, {0.7f, 0.0f, 0.17f}, {0.13f, 0.0f, 0.0f},
              {-0.05f, 0.0f, 0.57f}},
          .triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4}, {8, 9, 10},
              {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
              {18, 19, 16}, {20, 21, 22}, {22, 23, 20}}},
      {.color = {0.725f, 0.71f, 0.68f}});

  add_instance(scene, "tallbox", identity3x4f,
      {.positions    = {{-0.53f, 1.2f, 0.09f}, {0.04f, 1.2f, -0.09f},
              {-0.14f, 1.2f, -0.67f}, {-0.71f, 1.2f, -0.49f},
              {-0.53f, 0.0f, 0.09f}, {-0.53f, 1.2f, 0.09f}, {-0.71f, 1.2f, -0.49f},
              {-0.71f, 0.0f, -0.49f}, {-0.71f, 0.0f, -0.49f},
              {-0.71f, 1.2f, -0.49f}, {-0.14f, 1.2f, -0.67f},
              {-0.14f, 0.0f, -0.67f}, {-0.14f, 0.0f, -0.67f},
              {-0.14f, 1.2f, -0.67f}, {0.04f, 1.2f, -0.09f}, {0.04f, 0.0f, -0.09f},
              {0.04f, 0.0f, -0.09f}, {0.04f, 1.2f, -0.09f}, {-0.53f, 1.2f, 0.09f},
              {-0.53f, 0.0f, 0.09f}, {-0.53f, 0.0f, 0.09f}, {0.04f, 0.0f, -0.09f},
              {-0.14f, 0.0f, -0.67f}, {-0.71f, 0.0f, -0.49f}},
          .triangles = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4}, {8, 9, 10},
              {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
              {18, 19, 16}, {20, 21, 22}, {22, 23, 20}}},
      {.color = {0.725f, 0.71f, 0.68f}});

  add_instance(scene, "backwall", identity3x4f,
      {.positions    = {{-0.25f, 1.99f, 0.25f}, {-0.25f, 1.99f, -0.25f},
              {0.25f, 1.99f, -0.25f}, {0.25f, 1.99f, 0.25f}},
          .triangles = {{0, 1, 2}, {2, 3, 0}}},
      {.emission = {17, 12, 4}});

  return scene;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Find the first keyframe value that is greater than the argument.
inline int keyframe_index(const vector<float>& times, const float& time) {
  for (auto i : range(times.size()))
    if (times[i] > time) return (int)i;
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
