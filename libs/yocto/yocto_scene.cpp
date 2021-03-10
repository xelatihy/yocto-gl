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

#include "ext/stb_image_resize.h"
#include "yocto_cli.h"
#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_noise.h"
#include "yocto_parallel.h"
#include "yocto_sceneio.h"
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEXTURE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// pixel access
vec4f lookup_texture(
    const scene_texture& texture, int i, int j, bool as_linear) {
  auto color = vec4f{0, 0, 0, 0};
  if (!texture.pixelsf.empty()) {
    color = texture.pixelsf[j * texture.width + i];
  } else {
    color = byte_to_float(texture.pixelsb[j * texture.width + i]);
  }
  if (as_linear && !texture.linear) {
    return srgb_to_rgb(color);
  } else {
    return color;
  }
}

// Evaluates an image at a point `uv`.
vec4f eval_texture(const scene_texture& texture, const vec2f& uv,
    bool as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (texture.width == 0 || texture.height == 0) return {0, 0, 0, 0};

  // get texture width/height
  auto size = vec2i{texture.width, texture.height};

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

  // handle interpolation
  if (no_interpolation) {
    return lookup_texture(texture, i, j, as_linear);
  } else {
    return lookup_texture(texture, i, j, as_linear) * (1 - u) * (1 - v) +
           lookup_texture(texture, i, jj, as_linear) * (1 - u) * v +
           lookup_texture(texture, ii, j, as_linear) * u * (1 - v) +
           lookup_texture(texture, ii, jj, as_linear) * u * v;
  }
}

// Helpers
vec4f eval_texture(const scene_model& scene, int texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation, bool clamp_to_edge) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(
      scene.textures[texture], uv, ldr_as_linear, no_interpolation);
}

// conversion from image
scene_texture image_to_texture(const color_image& image) {
  auto texture = scene_texture{image.width, image.height, image.linear, {}, {}};
  if (image.linear) {
    texture.pixelsf = image.pixels;
  } else {
    texture.pixelsb.resize(image.pixels.size());
    float_to_byte(texture.pixelsb, image.pixels);
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
material_point eval_material(const scene_model& scene,
    const scene_material& material, const vec2f& texcoord,
    const vec4f& color_shp) {
  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

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
  if (material.type == scene_material_type::glass ||
      material.type == scene_material_type::volume ||
      material.type == scene_material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == scene_material_type::matte ||
      point.type == scene_material_type::metallic ||
      point.type == scene_material_type::plastic) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// check if a material is a delta or volumetric
bool is_delta(const scene_material& material) {
  return (material.type == scene_material_type::metal &&
             material.roughness == 0) ||
         (material.type == scene_material_type::glass &&
             material.roughness == 0) ||
         (material.type == scene_material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == scene_material_type::volume);
}
bool is_volumetric(const scene_material& material) {
  return material.type == scene_material_type::glass ||
         material.type == scene_material_type::volume ||
         material.type == scene_material_type::subsurface;
}

// check if a brdf is a delta
bool is_delta(const material_point& material) {
  return (material.type == scene_material_type::metal &&
             material.roughness == 0) ||
         (material.type == scene_material_type::glass &&
             material.roughness == 0) ||
         (material.type == scene_material_type::thinglass &&
             material.roughness == 0) ||
         (material.type == scene_material_type::volume);
}
bool has_volume(const material_point& material) {
  return material.type == scene_material_type::glass ||
         material.type == scene_material_type::volume ||
         material.type == scene_material_type::subsurface;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Interpolate vertex data
vec3f eval_position(const scene_shape& shape, int element, const vec2f& uv) {
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.positions[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.positions[line.x], shape.positions[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const scene_shape& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return normalize(shape.normals[point]);
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return normalize(
        interpolate_line(shape.normals[line.x], shape.normals[line.y], uv.x));
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return normalize(interpolate_triangle(shape.normals[triangle.x],
        shape.normals[triangle.y], shape.normals[triangle.z], uv));
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec3f eval_tangent(const scene_shape& shape, int element, const vec2f& uv) {
  return eval_normal(shape, element, uv);
}

vec2f eval_texcoord(const scene_shape& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.texcoords[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(
        shape.texcoords[line.x], shape.texcoords[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.texcoords[triangle.x],
        shape.texcoords[triangle.y], shape.texcoords[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

vec4f eval_color(const scene_shape& shape, int element, const vec2f& uv) {
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.colors[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.colors[line.x], shape.colors[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.colors[triangle.x],
        shape.colors[triangle.y], shape.colors[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.colors[quad.x], shape.colors[quad.y],
        shape.colors[quad.z], shape.colors[quad.w], uv);
  } else {
    return {0, 0};
  }
}

float eval_radius(const scene_shape& shape, int element, const vec2f& uv) {
  if (shape.radius.empty()) return 0;
  if (!shape.points.empty()) {
    auto& point = shape.points[element];
    return shape.radius[point];
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return interpolate_line(shape.radius[line.x], shape.radius[line.y], uv.x);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return interpolate_triangle(shape.radius[triangle.x],
        shape.radius[triangle.y], shape.radius[triangle.z], uv);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return interpolate_quad(shape.radius[quad.x], shape.radius[quad.y],
        shape.radius[quad.z], shape.radius[quad.w], uv);
  } else {
    return 0;
  }
}

// Evaluate element normals
vec3f eval_element_normal(const scene_shape& shape, int element) {
  if (!shape.points.empty()) {
    return {0, 0, 1};
  } else if (!shape.lines.empty()) {
    auto& line = shape.lines[element];
    return line_tangent(shape.positions[line.x], shape.positions[line.y]);
  } else if (!shape.triangles.empty()) {
    auto& triangle = shape.triangles[element];
    return triangle_normal(shape.positions[triangle.x],
        shape.positions[triangle.y], shape.positions[triangle.z]);
  } else if (!shape.quads.empty()) {
    auto& quad = shape.quads[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const scene_shape& shape) {
  if (!shape.points.empty()) {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    return lines_tangents(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return triangles_normals(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return quads_normals(shape.quads, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const scene_shape& shape) {
  if (!shape.points.empty()) {
    normals.assign(shape.positions.size(), {0, 0, 1});
  } else if (!shape.lines.empty()) {
    lines_tangents(normals, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    triangles_normals(normals, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    quads_normals(normals, shape.quads, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Shape sampling
vector<float> sample_shape_cdf(const scene_shape& shape) {
  if (!shape.points.empty()) {
    return sample_points_cdf((int)shape.points.size());
  } else if (!shape.lines.empty()) {
    return sample_lines_cdf(shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    return sample_triangles_cdf(shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    return sample_quads_cdf(shape.quads, shape.positions);
  } else {
    return sample_points_cdf((int)shape.positions.size());
  }
}

void sample_shape_cdf(vector<float>& cdf, const scene_shape& shape) {
  if (!shape.points.empty()) {
    sample_points_cdf(cdf, (int)shape.points.size());
  } else if (!shape.lines.empty()) {
    sample_lines_cdf(cdf, shape.lines, shape.positions);
  } else if (!shape.triangles.empty()) {
    sample_triangles_cdf(cdf, shape.triangles, shape.positions);
  } else if (!shape.quads.empty()) {
    sample_quads_cdf(cdf, shape.quads, shape.positions);
  } else {
    sample_points_cdf(cdf, (int)shape.positions.size());
  }
}

shape_point sample_shape(const scene_shape& shape, const vector<float>& cdf,
    float rn, const vec2f& ruv) {
  if (!shape.points.empty()) {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  } else if (!shape.lines.empty()) {
    auto [element, u] = sample_lines(cdf, rn, ruv.x);
    return {element, {u, 0}};
  } else if (!shape.triangles.empty()) {
    auto [element, uv] = sample_triangles(cdf, rn, ruv);
    return {element, uv};
  } else if (!shape.quads.empty()) {
    auto [element, uv] = sample_quads(cdf, rn, ruv);
    return {element, uv};
  } else {
    auto element = sample_points(cdf, rn);
    return {element, {0, 0}};
  }
}

vector<shape_point> sample_shape(
    const scene_shape& shape, int num_samples, uint64_t seed) {
  auto cdf    = sample_shape_cdf(shape);
  auto points = vector<shape_point>(num_samples);
  auto rng    = make_rng(seed);
  for (auto& point : points) {
    point = sample_shape(shape, cdf, rand1f(rng), rand2f(rng));
  }
  return points;
}

// Conversions
scene_shape quads_to_triangles(const scene_shape& shape) {
  auto result = shape;
  quads_to_triangles(result, result);
  return result;
}
void quads_to_triangles(scene_shape& result, const scene_shape& shape) {
  result.triangles = quads_to_triangles(shape.quads);
  result.quads     = {};
}

// Subdivision
scene_shape subdivide_shape(
    const scene_shape& shape, int subdivisions, bool catmullclark) {
  // This should probably be reimplemented in a faster fashion,
  // but how it is not obvious
  auto subdivided = scene_shape{};
  if (!shape.points.empty()) {
    // nothing to do
  } else if (!shape.lines.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_lines(
        shape.lines, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_lines(
        shape.lines, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_lines(
        shape.lines, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_lines(
        shape.lines, shape.radius, subdivisions);
    std::tie(subdivided.lines, subdivided.positions) = subdivide_lines(
        shape.lines, shape.positions, subdivisions);
  } else if (!shape.triangles.empty()) {
    std::tie(std::ignore, subdivided.normals) = subdivide_triangles(
        shape.triangles, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_triangles(
        shape.triangles, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_triangles(
        shape.triangles, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_triangles(
        shape.triangles, shape.radius, subdivisions);
    std::tie(subdivided.triangles, subdivided.positions) = subdivide_triangles(
        shape.triangles, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && !catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_quads(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_quads(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_quads(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_quads(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_quads(
        shape.quads, shape.positions, subdivisions);
  } else if (!shape.quads.empty() && catmullclark) {
    std::tie(std::ignore, subdivided.normals) = subdivide_catmullclark(
        shape.quads, shape.normals, subdivisions);
    std::tie(std::ignore, subdivided.texcoords) = subdivide_catmullclark(
        shape.quads, shape.texcoords, subdivisions);
    std::tie(std::ignore, subdivided.colors) = subdivide_catmullclark(
        shape.quads, shape.colors, subdivisions);
    std::tie(std::ignore, subdivided.radius) = subdivide_catmullclark(
        shape.quads, shape.radius, subdivisions);
    std::tie(subdivided.quads, subdivided.positions) = subdivide_catmullclark(
        shape.quads, shape.positions, subdivisions);
  } else {
    // empty shape
  }
  return subdivided;
}

// Interpolate vertex data
vec3f eval_position(const scene_fvshape& shape, int element, const vec2f& uv) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return interpolate_quad(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w], uv);
  } else {
    return {0, 0, 0};
  }
}

vec3f eval_normal(const scene_fvshape& shape, int element, const vec2f& uv) {
  if (shape.normals.empty()) return eval_element_normal(shape, element);
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadsnorm[element];
    return normalize(
        interpolate_quad(shape.normals[quad.x], shape.normals[quad.y],
            shape.normals[quad.z], shape.normals[quad.w], uv));
  } else {
    return {0, 0, 1};
  }
}

vec2f eval_texcoord(const scene_fvshape& shape, int element, const vec2f& uv) {
  if (shape.texcoords.empty()) return {0, 0};
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadstexcoord[element];
    return interpolate_quad(shape.texcoords[quad.x], shape.texcoords[quad.y],
        shape.texcoords[quad.z], shape.texcoords[quad.w], uv);
  } else {
    return {0, 0};
  }
}

// Evaluate element normals
vec3f eval_element_normal(const scene_fvshape& shape, int element) {
  if (!shape.quadspos.empty()) {
    auto& quad = shape.quadspos[element];
    return quad_normal(shape.positions[quad.x], shape.positions[quad.y],
        shape.positions[quad.z], shape.positions[quad.w]);
  } else {
    return {0, 0, 0};
  }
}

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const scene_fvshape& shape) {
  if (!shape.quadspos.empty()) {
    return quads_normals(shape.quadspos, shape.positions);
  } else {
    return vector<vec3f>(shape.positions.size(), {0, 0, 1});
  }
}
void compute_normals(vector<vec3f>& normals, const scene_fvshape& shape) {
  if (!shape.quadspos.empty()) {
    quads_normals(normals, shape.quadspos, shape.positions);
  } else {
    normals.assign(shape.positions.size(), {0, 0, 1});
  }
}

// Conversions
scene_shape fvshape_to_shape(const scene_fvshape& fvshape, bool as_triangles) {
  auto shape = scene_shape{};
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, fvshape.quadspos, fvshape.quadsnorm,
      fvshape.quadstexcoord, fvshape.positions, fvshape.normals,
      fvshape.texcoords);
  return shape;
}
scene_fvshape shape_to_fvshape(const scene_shape& shape) {
  if (!shape.points.empty() || !shape.lines.empty())
    throw std::invalid_argument{"cannor convert shape"};
  auto fvshape          = scene_fvshape{};
  fvshape.positions     = shape.positions;
  fvshape.normals       = shape.normals;
  fvshape.texcoords     = shape.texcoords;
  fvshape.quadspos      = !shape.quads.empty() ? shape.quads
                                               : triangles_to_quads(shape.triangles);
  fvshape.quadsnorm     = !shape.normals.empty() ? fvshape.quadspos
                                                 : vector<vec4i>{};
  fvshape.quadstexcoord = !shape.texcoords.empty() ? fvshape.quadspos
                                                   : vector<vec4i>{};
  return fvshape;
}

// Subdivision
scene_fvshape subdivide_fvshape(
    const scene_fvshape& shape, int subdivisions, bool catmullclark) {
  auto subdivided = scene_fvshape{};
  if (!catmullclark) {
    std::tie(subdivided.quadspos, subdivided.positions) = subdivide_quads(
        shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_quads(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) = subdivide_quads(
        shape.quadstexcoord, shape.texcoords, subdivisions);
  } else {
    std::tie(subdivided.quadspos, subdivided.positions) =
        subdivide_catmullclark(shape.quadspos, shape.positions, subdivisions);
    std::tie(subdivided.quadsnorm, subdivided.normals) = subdivide_catmullclark(
        shape.quadsnorm, shape.normals, subdivisions);
    std::tie(subdivided.quadstexcoord, subdivided.texcoords) =
        subdivide_catmullclark(
            shape.quadstexcoord, shape.texcoords, subdivisions, true);
  }
  return subdivided;
}

vector<string> shape_stats(const scene_shape& shape, bool verbose) {
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

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("points:       " + format(shape.points.size()));
  stats.push_back("lines:        " + format(shape.lines.size()));
  stats.push_back("triangles:    " + format(shape.triangles.size()));
  stats.push_back("quads:        " + format(shape.quads.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("colors:       " + format(shape.colors.size()));
  stats.push_back("radius:       " + format(shape.radius.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

vector<string> fvshape_stats(const scene_fvshape& shape, bool verbose) {
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

  auto bbox = invalidb3f;
  for (auto& pos : shape.positions) bbox = merge(bbox, pos);

  auto stats = vector<string>{};
  stats.push_back("fvquads:      " + format(shape.quadspos.size()));
  stats.push_back("positions:    " + format(shape.positions.size()));
  stats.push_back("normals:      " + format(shape.normals.size()));
  stats.push_back("texcoords:    " + format(shape.texcoords.size()));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));
  stats.push_back("min:          " + format3(bbox.min));
  stats.push_back("max:          " + format3(bbox.max));

  return stats;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INSTANCE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Eval position
vec3f eval_position(const scene_model& scene, const scene_instance& instance,
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
    const scene_model& scene, const scene_instance& instance, int element) {
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
vec3f eval_normal(const scene_model& scene, const scene_instance& instance,
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
vec2f eval_texcoord(const scene_model& scene, const scene_instance& instance,
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
    const scene_model& scene, const scene_instance& instance, int element) {
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

vec3f eval_normalmap(const scene_model& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalidid &&
      (!shape.triangles.empty() || !shape.quads.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, false));
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
vec3f eval_shading_normal(const scene_model& scene,
    const scene_instance& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty() || !shape.quads.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalidid) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == scene_material_type::glass) return normal;
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
vec4f eval_color(const scene_model& scene, const scene_instance& instance,
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

// Evaluate material
material_point eval_material(const scene_model& scene,
    const scene_instance& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

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
  if (material.type == scene_material_type::glass ||
      material.type == scene_material_type::volume ||
      material.type == scene_material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == scene_material_type::matte ||
      point.type == scene_material_type::metallic ||
      point.type == scene_material_type::plastic) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// check if an instance is volumetric
bool is_volumetric(const scene_model& scene, const scene_instance& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENVIRONMENT PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate environment color.
vec3f eval_environment(const scene_model& scene,
    const scene_environment& environment, const vec3f& direction) {
  auto wl       = transform_direction(inverse(environment.frame), direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
vec3f eval_environment(const scene_model& scene, const vec3f& direction) {
  auto emission = zero3f;
  for (auto environment : scene.environments) {
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
void add_camera(scene_model& scene) {
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

// Add a sky environment
void add_sky(scene_model& scene, float sun_angle) {
  scene.texture_names.emplace_back("sky");
  auto& texture = scene.textures.emplace_back();
  texture       = image_to_texture(make_sunsky(1024, 512, sun_angle));
  scene.environment_names.emplace_back("sky");
  auto& environment        = scene.environments.emplace_back();
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
}

// get named camera or default if camera is empty
int find_camera(const scene_model& scene, const string& name) {
  if (scene.cameras.empty()) return invalidid;
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
    if (scene.camera_names[idx] == "camera0") return idx;
  }
  for (auto idx = 0; idx < (int)scene.camera_names.size(); idx++) {
    if (scene.camera_names[idx] == "camera1") return idx;
  }
  return 0;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene_model& scene) {
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
    scene_shape& shape, scene_subdiv& subdiv_, const scene_model& scene) {
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

  if (subdiv.displacement != 0 && subdiv.displacement_tex != invalidid) {
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
            eval_texture(displacement_tex, subdiv.texcoords[qtxt[i]], false));
        if (!displacement_tex.pixelsb.empty()) disp -= 0.5f;
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
  split_facevarying(shape.quads, shape.positions, shape.normals,
      shape.texcoords, subdiv.quadspos, subdiv.quadsnorm, subdiv.quadstexcoord,
      subdiv.positions, subdiv.normals, subdiv.texcoords);
}

void tesselate_subdivs(scene_model& scene) {
  // handle progress
  auto progress = vec2i{0, (int)scene.subdivs.size() + 1};
  log_progress("tesselate subdivs", progress.x++, progress.y);

  // tesselate shapes
  for (auto& subdiv : scene.subdivs) {
    log_progress("tesselate subdiv", progress.x++, progress.y);
    tesselate_subdiv(scene.shapes[subdiv.shape], subdiv, scene);
  }

  // done
  log_progress("tesselate subdivs", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

size_t compute_memory(const scene_model& scene) {
  auto vector_memory = [](auto& values) -> size_t {
    if (values.empty()) return 0;
    return values.size() * sizeof(values[0]);
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

vector<string> scene_stats(const scene_model& scene, bool verbose) {
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
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
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
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene.
vector<string> scene_validation(const scene_model& scene, bool notextures) {
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
  auto check_empty_textures = [&errs](const scene_model& scene) {
    for (auto idx = 0; idx < (int)scene.textures.size(); idx++) {
      auto& texture = scene.textures[idx];
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
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a plane.
scene_shape make_rect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_rect(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
scene_shape make_bulged_rect(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = scene_shape{};
  make_bulged_rect(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a plane in the xz plane.
scene_shape make_recty(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_recty(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
scene_shape make_bulged_recty(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, float radius) {
  auto shape = scene_shape{};
  make_bulged_recty(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Make a box.
scene_shape make_box(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = scene_shape{};
  make_box(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}
scene_shape make_rounded_box(const vec3i& steps, const vec3f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = scene_shape{};
  make_rounded_box(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, radius);
  return shape;
}

// Make a quad stack
scene_shape make_rect_stack(
    const vec3i& steps, const vec3f& scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_rect_stack(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a floor.
scene_shape make_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}
scene_shape make_bent_floor(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale, float bent) {
  auto shape = scene_shape{};
  make_bent_floor(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, bent);
  return shape;
}

// Make a sphere.
scene_shape make_sphere(int steps, float scale, float uvscale) {
  auto shape = scene_shape{};
  make_sphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere.
scene_shape make_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_uvsphere(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a sphere with slipped caps.
scene_shape make_capped_uvsphere(
    const vec2i& steps, float scale, const vec2f& uvscale, float height) {
  auto shape = scene_shape{};
  make_capped_uvsphere(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, height);
  return shape;
}
// Make a disk
scene_shape make_disk(int steps, float scale, float uvscale) {
  auto shape = scene_shape{};
  make_disk(shape.quads, shape.positions, shape.normals, shape.texcoords, steps,
      scale, uvscale);
  return shape;
}

// Make a bulged disk
scene_shape make_bulged_disk(
    int steps, float scale, float uvscale, float height) {
  auto shape = scene_shape{};
  make_bulged_disk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale, height);
  return shape;
}

// Make a uv disk
scene_shape make_uvdisk(const vec2i& steps, float scale, const vec2f& uvscale) {
  auto shape = scene_shape{};
  make_uvdisk(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a uv cylinder
scene_shape make_uvcylinder(
    const vec3i& steps, const vec2f& scale, const vec3f& uvscale) {
  auto shape = scene_shape{};
  make_uvcylinder(shape.quads, shape.positions, shape.normals, shape.texcoords,
      steps, scale, uvscale);
  return shape;
}

// Make a rounded uv cylinder
scene_shape make_rounded_uvcylinder(const vec3i& steps, const vec2f& scale,
    const vec3f& uvscale, float radius) {
  auto shape = scene_shape{};
  make_rounded_uvcylinder(shape.quads, shape.positions, shape.normals,
      shape.texcoords, steps, scale, uvscale, radius);
  return shape;
}

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
scene_shape make_lines(const vec2i& steps, const vec2f& scale,
    const vec2f& uvscale, const vec2f& rad) {
  auto shape = scene_shape{};
  make_lines(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, steps, scale, uvscale, rad);
  return shape;
}

// Make point primitives. Returns points, pos, norm, texcoord, radius.
scene_shape make_point(float radius) {
  auto shape = scene_shape{};
  make_point(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, radius);
  return shape;
}

scene_shape make_points(int num, float uvscale, float radius) {
  auto shape = scene_shape{};
  make_points(shape.points, shape.positions, shape.normals, shape.texcoords,
      shape.radius, num, uvscale, radius);
  return shape;
}

scene_shape make_random_points(
    int num, const vec3f& size, float uvscale, float radius, uint64_t seed) {
  auto shape = scene_shape{};
  make_random_points(shape.points, shape.positions, shape.normals,
      shape.texcoords, shape.radius, num, size, uvscale, radius, seed);
  return shape;
}

// Make a facevarying rect
scene_fvshape make_fvrect(
    const vec2i& steps, const vec2f& scale, const vec2f& uvscale) {
  auto shape = scene_fvshape{};
  make_fvrect(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying box
scene_fvshape make_fvbox(
    const vec3i& steps, const vec3f& scale, const vec3f& uvscale) {
  auto shape = scene_fvshape{};
  make_fvbox(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Make a facevarying sphere
scene_fvshape make_fvsphere(int steps, float scale, float uvscale) {
  auto shape = scene_fvshape{};
  make_fvsphere(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, steps, scale, uvscale);
  return shape;
}

// Predefined meshes
scene_shape make_monkey(float scale) {
  auto shape = scene_shape{};
  make_monkey(shape.quads, shape.positions, scale);
  return shape;
}
scene_shape make_quad(float scale) {
  auto shape = scene_shape{};
  make_quad(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
scene_shape make_quady(float scale) {
  auto shape = scene_shape{};
  make_quady(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
scene_shape make_cube(float scale) {
  auto shape = scene_shape{};
  make_cube(
      shape.quads, shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
scene_fvshape make_fvcube(float scale) {
  auto shape = scene_fvshape{};
  make_fvcube(shape.quadspos, shape.quadsnorm, shape.quadstexcoord,
      shape.positions, shape.normals, shape.texcoords, scale);
  return shape;
}
scene_shape make_geosphere(float scale) {
  auto shape = scene_shape{};
  make_geosphere(shape.triangles, shape.positions, scale);
  return shape;
}

// Make a hair ball around a shape
scene_shape make_hair(const scene_shape& base, const vec2i& steps,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
  auto shape = scene_shape{};
  make_hair(shape.lines, shape.positions, shape.normals, shape.texcoords,
      shape.radius, base.triangles, base.quads, base.positions, base.normals,
      base.texcoords, steps, len, rad, noise, clump, rotation, seed);
  return shape;
}

// Make a heightfield mesh.
scene_shape make_heightfield(const vec2i& size, const vector<float>& height) {
  auto shape = scene_shape{};
  make_heightfield(shape.quads, shape.positions, shape.normals, shape.texcoords,
      size, height);
  return shape;
}
scene_shape make_heightfield(const vec2i& size, const vector<vec4f>& color) {
  auto shape = scene_shape{};
  make_heightfield(shape.quads, shape.positions, shape.normals, shape.texcoords,
      size, color);
  return shape;
}

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res and without texcoords and normals.
scene_shape points_to_spheres(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = scene_shape{};
  points_to_spheres(shape.quads, shape.positions, shape.normals,
      shape.texcoords, vertices, steps, scale);
  return shape;
}
scene_shape lines_to_cylinders(
    const vector<vec3f>& vertices, int steps, float scale) {
  auto shape = scene_shape{};
  lines_to_cylinders(shape.quads, shape.positions, shape.normals,
      shape.texcoords, vertices, steps, scale);
  return shape;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox(scene_model& scene) {
  auto& camera    = scene.cameras.emplace_back();
  camera.frame    = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 3.9}};
  camera.lens     = 0.035;
  camera.aperture = 0.0;
  camera.focus    = 3.9;
  camera.film     = 0.024;
  camera.aspect   = 1;

  auto& floor_shape       = scene.shapes.emplace_back();
  floor_shape.positions   = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& floor_material    = scene.materials.emplace_back();
  floor_material.color    = {0.725, 0.71, 0.68};
  auto& floor_instance    = scene.instances.emplace_back();
  floor_instance.shape    = (int)scene.shapes.size() - 1;
  floor_instance.material = (int)scene.materials.size() - 1;

  auto& ceiling_shape       = scene.shapes.emplace_back();
  ceiling_shape.positions   = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& ceiling_material    = scene.materials.emplace_back();
  ceiling_material.color    = {0.725, 0.71, 0.68};
  auto& ceiling_instance    = scene.instances.emplace_back();
  ceiling_instance.shape    = (int)scene.shapes.size() - 1;
  ceiling_instance.material = (int)scene.materials.size() - 1;

  auto& backwall_shape     = scene.shapes.emplace_back();
  backwall_shape.positions = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& backwall_material  = scene.materials.emplace_back();
  backwall_material.color  = {0.725, 0.71, 0.68};
  auto& backwall_instance  = scene.instances.emplace_back();
  backwall_instance.shape  = (int)scene.shapes.size() - 1;
  backwall_instance.material = (int)scene.materials.size() - 1;

  auto& rightwall_shape       = scene.shapes.emplace_back();
  rightwall_shape.positions   = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& rightwall_material    = scene.materials.emplace_back();
  rightwall_material.color    = {0.14, 0.45, 0.091};
  auto& rightwall_instance    = scene.instances.emplace_back();
  rightwall_instance.shape    = (int)scene.shapes.size() - 1;
  rightwall_instance.material = (int)scene.materials.size() - 1;

  auto& leftwall_shape     = scene.shapes.emplace_back();
  leftwall_shape.positions = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall_shape.triangles = {{0, 1, 2}, {2, 3, 0}};
  auto& leftwall_material  = scene.materials.emplace_back();
  leftwall_material.color  = {0.63, 0.065, 0.05};
  auto& leftwall_instance  = scene.instances.emplace_back();
  leftwall_instance.shape  = (int)scene.shapes.size() - 1;
  leftwall_instance.material = (int)scene.materials.size() - 1;

  auto& shortbox_shape       = scene.shapes.emplace_back();
  shortbox_shape.positions   = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  shortbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& shortbox_material    = scene.materials.emplace_back();
  shortbox_material.color    = {0.725, 0.71, 0.68};
  auto& shortbox_instance    = scene.instances.emplace_back();
  shortbox_instance.shape    = (int)scene.shapes.size() - 1;
  shortbox_instance.material = (int)scene.materials.size() - 1;

  auto& tallbox_shape       = scene.shapes.emplace_back();
  tallbox_shape.positions   = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
      {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
      {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
      {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
      {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
      {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
      {-0.71, 0.0, -0.49}};
  tallbox_shape.triangles   = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  auto& tallbox_material    = scene.materials.emplace_back();
  tallbox_material.color    = {0.725, 0.71, 0.68};
  auto& tallbox_instance    = scene.instances.emplace_back();
  tallbox_instance.shape    = (int)scene.shapes.size() - 1;
  tallbox_instance.material = (int)scene.materials.size() - 1;

  auto& light_shape       = scene.shapes.emplace_back();
  light_shape.positions   = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
      {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  light_shape.triangles   = {{0, 1, 2}, {2, 3, 0}};
  auto& light_material    = scene.materials.emplace_back();
  light_material.emission = {17, 12, 4};
  auto& light_instance    = scene.instances.emplace_back();
  light_instance.shape    = (int)scene.shapes.size() - 1;
  light_instance.material = (int)scene.materials.size() - 1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
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
