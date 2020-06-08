//
// Implementation for Yocto/Trace.
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

#include "yocto_trace.h"

#include <atomic>
#include <cstring>
#include <deque>
#include <future>
#include <memory>
#include <mutex>

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

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
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

// constant values
static const auto coat_ior       = 1.5f;
static const auto coat_roughness = 0.03f * 0.03f;

// Shape element normal.
static vec3f eval_normal(const scene_shape* shape, int element) {
  auto norm = zero3f;
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    norm   = triangle_normal(
        shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
  } else if (!shape->quads.empty()) {
    auto q = shape->quads[element];
    norm   = quad_normal(shape->positions[q.x], shape->positions[q.y],
        shape->positions[q.z], shape->positions[q.w]);
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    norm   = line_tangent(shape->positions[l.x], shape->positions[l.y]);
  } else if (!shape->points.empty()) {
    norm = {0, 0, 1};
  } else {
    // empty shape
    norm = {0, 0, 1};
  }
  return norm;
}

// Shape element normal.
static std::pair<vec3f, vec3f> eval_tangents(
    const scene_shape* shape, int element, const vec2f& uv) {
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    if (shape->texcoords.empty()) {
      return triangle_tangents_fromuv(shape->positions[t.x],
          shape->positions[t.y], shape->positions[t.z], {0, 0}, {1, 0}, {0, 1});
    } else {
      return triangle_tangents_fromuv(shape->positions[t.x],
          shape->positions[t.y], shape->positions[t.z], shape->texcoords[t.x],
          shape->texcoords[t.y], shape->texcoords[t.z]);
    }
  } else if (!shape->quads.empty()) {
    auto q = shape->quads[element];
    if (shape->texcoords.empty()) {
      return quad_tangents_fromuv(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w], {0, 0}, {1, 0}, {0, 1},
          {1, 1}, uv);
    } else {
      return quad_tangents_fromuv(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w], shape->texcoords[q.x],
          shape->texcoords[q.y], shape->texcoords[q.z], shape->texcoords[q.w],
          uv);
    }
  } else {
    return {zero3f, zero3f};
  }
}

// Shape value interpolated using barycentric coordinates
template <typename T>
static T eval_shape(const scene_shape* shape, const vector<T>& vals,
    int element, const vec2f& uv, const T& def) {
  if (vals.empty()) return def;
  if (!shape->triangles.empty()) {
    auto t = shape->triangles[element];
    return interpolate_triangle(vals[t.x], vals[t.y], vals[t.z], uv);
  } else if (!shape->quads.empty()) {
    auto q = shape->quads[element];
    if (q.w == q.z)
      return interpolate_triangle(vals[q.x], vals[q.y], vals[q.z], uv);
    return interpolate_quad(vals[q.x], vals[q.y], vals[q.z], vals[q.w], uv);
  } else if (!shape->lines.empty()) {
    auto l = shape->lines[element];
    return interpolate_line(vals[l.x], vals[l.y], uv.x);
  } else if (!shape->points.empty()) {
    return vals[shape->points[element]];
  } else {
    return def;
  }
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

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_perspective_camera(
    const scene_camera* camera, const vec2f& image_uv, const vec2f& lens_uv) {
  // focal plane correction --- we skip it since we consider the lens to be
  // the effective focal lens
  // if (camera->focus < flt_max) {
  //   distance = camera->lens * camera->focus / (camera->focus - camera->lens);
  // }
  // point on the image plane
  auto film = camera->aspect >= 1
                  ? vec2f{camera->film, camera->film / camera->aspect}
                  : vec2f{camera->film * camera->aspect, camera->film};
  auto q = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), camera->lens};
  // ray direction through the lens center
  auto dc = -normalize(q);
  // point on the lens
  auto e = vec3f{
      lens_uv.x * camera->aperture / 2, lens_uv.y * camera->aperture / 2, 0};
  // point on the focus plane
  auto p = dc * camera->focus / abs(dc.z);
  // correct ray direction to account for camera focusing
  auto d = normalize(p - e);
  // done
  return ray3f{
      transform_point(camera->frame, e), transform_direction(camera->frame, d)};
  // old implementation that was derived differently --- kept here in case
  // bugs start to show up
  // auto e = vec3f{(lens_uv.x - 0.5f) * camera->aperture,
  //     (lens_uv.y - 0.5f) * camera->aperture, 0};
  // auto q = vec3f{camera->film.x * (0.5f - image_uv.x),
  //     camera->film.y * (image_uv.y - 0.5f), distance};
  // // distance of the image of the point
  // auto distance1 = camera->lens * distance / (distance - camera->lens);
  // auto q1        = -q * distance1 / distance;
  // auto d         = normalize(q1 - e);
  // // auto q1 = - normalize(q) * camera->focus / normalize(q).z;
  // auto ray = ray3f{transform_point(camera->frame, e),
  //     transform_direction(camera->frame, d)};
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_orthographic_camera(
    const scene_camera* camera, const vec2f& image_uv, const vec2f& lens_uv) {
  // point on the image plane
  auto film = camera->aspect >= 1
                  ? vec2f{camera->film, camera->film / camera->aspect}
                  : vec2f{camera->film * camera->aspect, camera->film};
  auto scale = 1 / camera->lens;
  auto q     = vec3f{film.x * (0.5f - image_uv.x) * scale,
      film.y * (image_uv.y - 0.5f) * scale, camera->lens};
  // point on the lens
  auto e = vec3f{-q.x, -q.y, 0} + vec3f{lens_uv.x * camera->aperture / 2,
                                      lens_uv.y * camera->aperture / 2, 0};
  // point on the focus plane
  auto p = vec3f{-q.x, -q.y, -camera->focus};
  // correct ray direction to account for camera focusing
  auto d = normalize(p - e);
  // done
  return ray3f{
      transform_point(camera->frame, e), transform_direction(camera->frame, d)};
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_camera(
    const scene_camera* camera, const vec2f& uv, const vec2f& luv) {
  if (camera->orthographic)
    return eval_orthographic_camera(camera, uv, luv);
  else
    return eval_perspective_camera(camera, uv, luv);
}

// Sample camera
static ray3f sample_camera(const scene_camera* camera, const vec2i& ij,
    const vec2i& image_size, const vec2f& puv, const vec2f& luv, bool tent) {
  if (!tent) {
    auto uv = vec2f{
        (ij.x + puv.x) / image_size.x, (ij.y + puv.y) / image_size.y};
    return eval_camera(camera, uv, sample_disk(luv));
  } else {
    const auto width  = 2.0f;
    const auto offset = 0.5f;
    auto       fuv =
        width *
            vec2f{
                puv.x < 0.5f ? sqrt(2 * puv.x) - 1 : 1 - sqrt(2 - 2 * puv.x),
                puv.y - 0.5f ? sqrt(2 * puv.y) - 1 : 1 - sqrt(2 - 2 * puv.y),
            } +
        offset;
    auto uv = vec2f{
        (ij.x + fuv.x) / image_size.x, (ij.y + fuv.y) / image_size.y};
    return eval_camera(camera, uv, sample_disk(luv));
  }
}

// Point used for tracing
struct trace_point {
  // shape
  int   object   = -1;
  vec3f position = zero3f;
  vec3f normal   = zero3f;
  vec3f gnormal  = zero3f;
  vec2f texcoord = zero2f;
  vec3f color    = zero3f;
  // directions
  vec3f outgoing = zero3f;
  vec3f incoming = zero3f;
  // material
  vec3f emission     = {0, 0, 0};
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f metal        = {0, 0, 0};
  vec3f coat         = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  vec3f translucency = {0, 0, 0};
  vec3f refraction   = {0, 0, 0};
  float roughness    = 0;
  float opacity      = 1;
  float ior          = 1;
  vec3f meta         = {0, 0, 0};
  vec3f metak        = {0, 0, 0};
  // weights
  float diffuse_pdf      = 0;
  float specular_pdf     = 0;
  float metal_pdf        = 0;
  float coat_pdf         = 0;
  float transmission_pdf = 0;
  float translucency_pdf = 0;
  float refraction_pdf   = 0;
};

// Evaluate point
static trace_point eval_point(const scene_model* scene,
    const scene_intersection& intersection, const ray3f& ray) {
  // get data
  auto object   = scene->objects[intersection.object];
  auto shape    = object->shape;
  auto material = object->material;
  auto frame = object->instance->frames[intersection.instance] * object->frame;
  auto element          = intersection.element;
  auto uv               = intersection.uv;
  auto non_rigid_frames = true;

  // initialize point
  auto point     = trace_point{};
  point.object   = intersection.object;
  point.outgoing = -ray.d;
  point.incoming = -ray.d;

  // geometric properties
  point.position = eval_shape(shape, shape->positions, element, uv, zero3f);
  point.gnormal  = eval_normal(shape, element);
  point.normal   = normalize(
      eval_shape(shape, shape->normals, element, uv, point.gnormal));
  point.texcoord = eval_shape(shape, shape->texcoords, element, uv, uv);
  point.color    = eval_shape(shape, shape->colors, element, uv, vec3f{1});

  // apply normal mapping
  if (material->normal_tex &&
      (!shape->triangles.empty() || !shape->quads.empty())) {
    auto normalmap = -1 + 2 * eval_texture(
                                  material->normal_tex, point.texcoord, true);
    auto z      = point.normal;
    auto basis  = identity3x3f;
    auto flip_v = false;
    if (shape->tangents.empty()) {
      auto tangents = eval_tangents(shape, element, uv);
      auto x        = orthonormalize(tangents.first, z);
      auto y        = normalize(cross(z, x));
      basis         = {x, y, z};
      flip_v        = dot(y, tangents.second) < 0;
    } else {
      auto tangsp = eval_shape(
          shape, shape->tangents, element, uv, {0, 0, 0, 1});
      auto x = orthonormalize(xyz(tangsp), z);
      auto y = normalize(cross(z, x));
      basis  = {x, y, z};
      flip_v = tangsp.w < 0;
    }
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    point.normal = normalize(basis * normalmap);
  }

  // transforms
  point.position = transform_point(frame, point.position);
  point.normal   = transform_normal(frame, point.normal, non_rigid_frames);
  point.gnormal  = transform_normal(frame, point.gnormal, non_rigid_frames);

  // correct normals
  if (!shape->points.empty()) {
    point.normal = point.outgoing;
  } else if (!shape->lines.empty()) {
    point.normal = orthonormalize(point.outgoing, point.normal);
  } else if (!shape->triangles.empty()) {
    if (material->thin && dot(point.outgoing, point.normal) < 0)
      point.normal = -point.normal;
  } else if (!shape->quads.empty()) {
    if (material->thin && dot(point.outgoing, point.normal) < 0)
      point.normal = -point.normal;
  }

  // material -------
  // initialize factors
  auto texcoord = point.texcoord;
  auto emission = material->emission *
                  eval_texture(material->emission_tex, texcoord, false);
  auto base = material->color * point.color *
              eval_texture(material->color_tex, texcoord, false);
  auto specular = material->specular *
                  eval_texture(material->specular_tex, texcoord, true).x;
  auto metallic = material->metallic *
                  eval_texture(material->metallic_tex, texcoord, true).x;
  auto roughness = material->roughness *
                   eval_texture(material->roughness_tex, texcoord, true).x;

  auto ior  = material->ior;
  auto coat = material->coat *
              eval_texture(material->coat_tex, texcoord, true).x;
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord, true).x;
  auto translucency =
      material->translucency *
      eval_texture(material->translucency_tex, texcoord, true).x;
  auto opacity = material->opacity *
                 mean(eval_texture(material->opacity_tex, texcoord, true));
  auto thin = material->thin || !material->transmission;

  // factors
  auto weight    = vec3f{1, 1, 1};
  point.emission = weight * emission;
  point.coat     = weight * coat;
  weight *= 1 - point.coat *
                    fresnel_dielectric(coat_ior, point.outgoing, point.normal);
  point.metal = weight * metallic;
  weight *= 1 - metallic;
  point.refraction = thin ? zero3f : (weight * transmission);
  weight *= 1 - (thin ? 0 : transmission);
  point.specular = weight * specular;
  weight *= 1 -
            specular * fresnel_dielectric(ior, point.outgoing, point.normal);
  point.transmission = thin ? (weight * transmission * base) : zero3f;
  weight *= 1 - (thin ? transmission : 0);
  point.translucency = thin ? (weight * translucency * base)
                            : (weight * translucency);
  weight *= 1 - translucency;
  point.diffuse   = weight * base;
  point.meta      = reflectivity_to_eta(base);
  point.metak     = zero3f;
  point.roughness = roughness * roughness;
  point.ior       = ior;
  point.opacity   = opacity;

  // textures
  if (point.diffuse != zero3f || point.translucency != zero3f ||
      point.roughness) {
    point.roughness = clamp(point.roughness, 0.03f * 0.03f, 1.0f);
  }
  if (point.specular == zero3f && point.metal == zero3f &&
      point.transmission == zero3f && point.refraction == zero3f) {
    point.roughness = 1;
  }
  if (point.opacity > 0.999f) point.opacity = 1;

  // weights
  point.diffuse_pdf  = max(point.diffuse);
  point.specular_pdf = max(point.specular * fresnel_dielectric(point.ior,
                                                point.normal, point.outgoing));
  point.metal_pdf = max(point.metal * fresnel_conductor(point.meta, point.metak,
                                          point.normal, point.outgoing));
  point.coat_pdf  = max(
      point.coat * fresnel_dielectric(coat_ior, point.normal, point.outgoing));
  point.transmission_pdf = max(point.transmission);
  point.translucency_pdf = max(point.translucency);
  point.refraction_pdf   = max(point.refraction);
  auto pdf_sum = point.diffuse_pdf + point.specular_pdf + point.metal_pdf +
                 point.coat_pdf + point.transmission_pdf +
                 point.translucency_pdf + point.refraction_pdf;
  if (pdf_sum) {
    point.diffuse_pdf /= pdf_sum;
    point.specular_pdf /= pdf_sum;
    point.metal_pdf /= pdf_sum;
    point.coat_pdf /= pdf_sum;
    point.transmission_pdf /= pdf_sum;
    point.translucency_pdf /= pdf_sum;
    point.refraction_pdf /= pdf_sum;
  }
  return point;
}

// Point used for tracing
struct volume_point {
  // shape
  int   object   = -1;
  vec3f position = zero3f;
  // directions
  vec3f outgoing = zero3f;
  vec3f incoming = zero3f;
  // volume
  vec3f voldensity    = {0, 0, 0};
  vec3f volemission   = {0, 0, 0};
  vec3f volscatter    = {0, 0, 0};
  float volanisotropy = 0;
};

// Evaluate point
static volume_point eval_volume(const scene_model* scene,
    const scene_intersection& intersection, const ray3f& ray) {
  // get data
  auto& object   = scene->objects[intersection.object];
  auto& shape    = object->shape;
  auto& material = object->material;
  auto  frame = object->instance->frames[intersection.instance] * object->frame;
  auto  element = intersection.element;
  auto  uv      = intersection.uv;

  // initialize point
  auto point     = volume_point{};
  point.object   = intersection.object;
  point.outgoing = -ray.d;
  point.incoming = -ray.d;

  // geometric properties
  point.position = eval_shape(shape, shape->positions, element, uv, zero3f);
  point.position = transform_point(frame, point.position);

  // material -------
  // initialize factors
  auto texcoord = eval_shape(shape, shape->texcoords, element, uv, uv);
  auto color    = eval_shape(shape, shape->colors, element, uv, vec3f{1});
  auto base     = material->color * color *
              eval_texture(material->color_tex, texcoord, false);
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord, true).x;
  auto translucency =
      material->translucency *
      eval_texture(material->translucency_tex, texcoord, true).x;
  auto thin = material->thin ||
              (!material->transmission && !material->translucency);
  auto scattering = material->scattering *
                    eval_texture(material->scattering_tex, texcoord, false);
  auto scanisotropy = material->scanisotropy;
  auto trdepth      = material->trdepth;

  // factors
  point.volemission = zero3f;
  point.voldensity  = ((transmission || translucency) && !thin)
                         ? -log(clamp(base, 0.0001f, 1.0f)) / trdepth
                         : zero3f;
  point.volscatter    = scattering;
  point.volanisotropy = scanisotropy;

  return point;
}

// Check if an instance as volume scattering
static bool has_volume(
    const scene_model* scene, const scene_intersection& intersection) {
  auto object = scene->objects[intersection.object];
  return !object->material->thin &&
         (object->material->transmission || object->material->translucency);
}

// Evaluate all environment color.
static vec3f eval_environment(const scene_model* scene, const ray3f& ray) {
  auto emission = zero3f;
  for (auto environment : scene->environments) {
    auto wl       = transform_direction(inverse(environment->frame), ray.d);
    auto texcoord = vec2f{
        atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
    if (texcoord.x < 0) texcoord.x += 1;
    emission += environment->emission *
                eval_texture(environment->emission_tex, texcoord);
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Set non-rigid frames as default
static const bool non_rigid_frames = true;

static vec3f eval_emission(const trace_point& point) { return point.emission; }

static vec3f eval_volemission(const volume_point& point) {
  return point.volemission;
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_brdfcos(const trace_point& point) {
  if (!point.roughness) return zero3f;

  // accumulate the lobes
  auto brdfcos = zero3f;
  if (point.diffuse) {
    brdfcos += point.diffuse * eval_diffuse_reflection(point.normal,
                                   point.outgoing, point.incoming);
  }
  if (point.specular) {
    brdfcos += point.specular * eval_microfacet_reflection(point.ior,
                                    point.roughness, point.normal,
                                    point.outgoing, point.incoming);
  }
  if (point.metal) {
    brdfcos += point.metal * eval_microfacet_reflection(point.meta, point.metak,
                                 point.roughness, point.normal, point.outgoing,
                                 point.incoming);
  }
  if (point.coat) {
    brdfcos += point.coat * eval_microfacet_reflection(coat_ior, coat_roughness,
                                point.normal, point.outgoing, point.incoming);
  }
  if (point.transmission) {
    brdfcos += point.transmission * eval_microfacet_transmission(point.ior,
                                        point.roughness, point.normal,
                                        point.outgoing, point.incoming);
  }
  if (point.translucency) {
    brdfcos += point.translucency * eval_diffuse_transmission(point.normal,
                                        point.outgoing, point.incoming);
  }
  if (point.refraction) {
    brdfcos += point.refraction * eval_microfacet_refraction(point.ior,
                                      point.roughness, point.normal,
                                      point.outgoing, point.incoming);
  }
  return brdfcos;
}

static vec3f eval_delta(const trace_point& point) {
  if (point.roughness) return zero3f;

  auto brdfcos = zero3f;

  if (point.specular && !point.refraction) {
    brdfcos += point.specular * eval_delta_reflection(point.ior, point.normal,
                                    point.outgoing, point.incoming);
  }
  if (point.metal) {
    brdfcos += point.metal * eval_delta_reflection(point.meta, point.metak,
                                 point.normal, point.outgoing, point.incoming);
  }
  if (point.coat) {
    brdfcos += point.coat * eval_delta_reflection(coat_ior, point.normal,
                                point.outgoing, point.incoming);
  }
  if (point.transmission) {
    brdfcos += point.transmission * eval_delta_transmission(point.ior,
                                        point.normal, point.outgoing,
                                        point.incoming);
  }
  if (point.refraction) {
    brdfcos += point.refraction * eval_delta_refraction(point.ior, point.normal,
                                      point.outgoing, point.incoming);
  }

  return brdfcos;
}

// Picks a direction based on the BRDF
static vec3f sample_brdf(const trace_point& point, float rnl, const vec2f& rn) {
  if (!point.roughness) return zero3f;

  auto cdf = 0.0f;

  if (point.diffuse_pdf) {
    cdf += point.diffuse_pdf;
    if (rnl < cdf)
      return sample_diffuse_reflection(point.normal, point.outgoing, rn);
  }

  if (point.specular_pdf && !point.refraction_pdf) {
    cdf += point.specular_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          point.ior, point.roughness, point.normal, point.outgoing, rn);
  }

  if (point.metal_pdf) {
    cdf += point.metal_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(point.meta, point.metak,
          point.roughness, point.normal, point.outgoing, rn);
  }

  if (point.coat_pdf) {
    cdf += point.coat_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          coat_ior, coat_roughness, point.normal, point.outgoing, rn);
  }

  if (point.transmission_pdf) {
    cdf += point.transmission_pdf;
    if (rnl < cdf)
      return sample_microfacet_transmission(
          point.ior, point.roughness, point.normal, point.outgoing, rn);
  }

  if (point.translucency_pdf) {
    cdf += point.translucency_pdf;
    if (rnl < cdf)
      return sample_diffuse_transmission(point.normal, point.outgoing, rn);
  }

  if (point.refraction_pdf) {
    cdf += point.refraction_pdf;
    if (rnl < cdf)
      return sample_microfacet_refraction(
          point.ior, point.roughness, point.normal, point.outgoing, rnl, rn);
  }

  return zero3f;
}

static vec3f sample_delta(const trace_point& point, float rnl) {
  if (point.roughness) return zero3f;

  // keep a weight sum to pick a lobe
  auto cdf = 0.0f;
  cdf += point.diffuse_pdf;

  if (point.specular_pdf && !point.refraction_pdf) {
    cdf += point.specular_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(point.ior, point.normal, point.outgoing);
    }
  }

  if (point.metal_pdf) {
    cdf += point.metal_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(
          point.meta, point.metak, point.normal, point.outgoing);
    }
  }

  if (point.coat_pdf) {
    cdf += point.coat_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(coat_ior, point.normal, point.outgoing);
    }
  }

  if (point.transmission_pdf) {
    cdf += point.transmission_pdf;
    if (rnl < cdf) {
      return sample_delta_transmission(point.ior, point.normal, point.outgoing);
    }
  }

  if (point.refraction_pdf) {
    cdf += point.refraction_pdf;
    if (rnl < cdf) {
      return sample_delta_refraction(
          point.ior, point.normal, point.outgoing, rnl);
    }
  }

  return zero3f;
}

// Compute the weight for sampling the BRDF
static float sample_brdf_pdf(const trace_point& point) {
  if (!point.roughness) return 0;

  auto pdf = 0.0f;

  if (point.diffuse_pdf) {
    pdf += point.diffuse_pdf * sample_diffuse_reflection_pdf(point.normal,
                                   point.outgoing, point.incoming);
  }

  if (point.specular_pdf && !point.refraction_pdf) {
    pdf += point.specular_pdf * sample_microfacet_reflection_pdf(point.ior,
                                    point.roughness, point.normal,
                                    point.outgoing, point.incoming);
  }

  if (point.metal_pdf) {
    pdf += point.metal_pdf * sample_microfacet_reflection_pdf(point.meta,
                                 point.metak, point.roughness, point.normal,
                                 point.outgoing, point.incoming);
  }

  if (point.coat_pdf) {
    pdf += point.coat_pdf * sample_microfacet_reflection_pdf(coat_ior,
                                coat_roughness, point.normal, point.outgoing,
                                point.incoming);
  }

  if (point.transmission_pdf) {
    pdf += point.transmission_pdf *
           sample_microfacet_transmission_pdf(point.ior, point.roughness,
               point.normal, point.outgoing, point.incoming);
  }

  if (point.translucency_pdf) {
    pdf += point.translucency_pdf *
           sample_diffuse_transmission_pdf(
               point.normal, point.outgoing, point.incoming);
  }

  if (point.refraction_pdf) {
    pdf += point.refraction_pdf * sample_microfacet_refraction_pdf(point.ior,
                                      point.roughness, point.normal,
                                      point.outgoing, point.incoming);
  }

  return pdf;
}

static float sample_delta_pdf(const trace_point& point) {
  if (point.roughness) return 0;

  auto pdf = 0.0f;
  if (point.specular_pdf && !point.refraction_pdf) {
    pdf += point.specular_pdf * sample_delta_reflection_pdf(point.ior,
                                    point.normal, point.outgoing,
                                    point.incoming);
  }
  if (point.metal_pdf) {
    pdf += point.metal_pdf * sample_delta_reflection_pdf(point.meta,
                                 point.metak, point.normal, point.outgoing,
                                 point.incoming);
  }
  if (point.coat_pdf) {
    pdf += point.coat_pdf * sample_delta_reflection_pdf(coat_ior, point.normal,
                                point.outgoing, point.incoming);
  }
  if (point.transmission_pdf) {
    pdf += point.transmission_pdf * sample_delta_transmission_pdf(point.ior,
                                        point.normal, point.outgoing,
                                        point.incoming);
  }
  if (point.refraction_pdf) {
    pdf += point.refraction_pdf * sample_delta_refraction_pdf(point.ior,
                                      point.normal, point.outgoing,
                                      point.incoming);
  }
  return pdf;
}

static vec3f eval_scattering(const volume_point& point) {
  if (point.voldensity == zero3f) return zero3f;
  return point.volscatter * point.voldensity *
         eval_phasefunction(
             point.volanisotropy, point.outgoing, point.incoming);
}

static vec3f sample_scattering(
    const volume_point& point, float rnl, const vec2f& rn) {
  if (point.voldensity == zero3f) return zero3f;
  return sample_phasefunction(point.volanisotropy, point.outgoing, rn);
}

static float sample_scattering_pdf(const volume_point& point) {
  if (point.voldensity == zero3f) return 0;
  return sample_phasefunction_pdf(
      point.volanisotropy, point.outgoing, point.incoming);
}

// Sample lights wrt solid angle
static vec3f sample_lights(const scene_model* scene, const vec3f& position,
    float rl, float rel, const vec2f& ruv) {
  auto  light_id = sample_uniform(scene->lights.size(), rl);
  auto& light    = scene->lights[light_id];
  if (light->object) {
    auto& object    = light->object;
    auto  shape     = object->shape;
    auto  frame     = object->instance->frames[light->instance] * object->frame;
    auto  element   = sample_discrete_cdf(shape->elements_cdf, rel);
    auto  uv        = (!shape->triangles.empty()) ? sample_triangle(ruv) : ruv;
    auto  lposition = transform_point(
        frame, eval_shape(shape, shape->positions, element, uv, zero3f));
    return normalize(lposition - position);
  } else if (light->environment) {
    auto& environment = light->environment;
    if (environment->emission_tex) {
      auto emission_tex = environment->emission_tex;
      auto idx          = sample_discrete_cdf(environment->texels_cdf, rel);
      auto size         = texture_size(emission_tex);
      auto uv           = vec2f{
          (idx % size.x + 0.5f) / size.x, (idx / size.x + 0.5f) / size.y};
      return transform_direction(environment->frame,
          {cos(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif),
              sin(uv.x * 2 * pif) * sin(uv.y * pif)});
    } else {
      return sample_sphere(ruv);
    }
  } else {
    return zero3f;
  }
}

// Sample lights pdf
static float sample_lights_pdf(
    const scene_model* scene, const vec3f& position, const vec3f& direction) {
  auto pdf = 0.0f;
  for (auto& light : scene->lights) {
    if (light->object) {
      // check all intersection
      auto  lpdf          = 0.0f;
      auto  next_position = position;
      auto& object        = light->object;
      auto  frame = object->instance->frames[light->instance] * object->frame;
      for (auto bounce = 0; bounce < 100; bounce++) {
        auto intersection = intersect_instance_bvh(
            light->object, light->instance, {next_position, direction});
        if (!intersection.hit) break;
        // accumulate pdf
        auto lposition = transform_point(
            frame, eval_shape(object->shape, object->shape->positions,
                       intersection.element, intersection.uv, zero3f));
        auto lnormal = transform_normal(frame,
            eval_normal(object->shape, intersection.element), non_rigid_frames);
        // prob triangle * area triangle = area triangle mesh
        auto area = object->shape->elements_cdf.back();
        lpdf += distance_squared(lposition, position) /
                (abs(dot(lnormal, direction)) * area);
        // continue
        next_position = lposition + direction * 1e-3f;
      }
      pdf += lpdf;
    } else if (light->environment) {
      auto& environment = light->environment;
      if (environment->emission_tex) {
        auto& cdf          = environment->texels_cdf;
        auto  emission_tex = environment->emission_tex;
        auto  size         = texture_size(emission_tex);
        auto  wl = transform_direction(inverse(environment->frame), direction);
        auto  texcoord = vec2f{atan2(wl.z, wl.x) / (2 * pif),
            acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
        if (texcoord.x < 0) texcoord.x += 1;
        auto i     = clamp((int)(texcoord.x * size.x), 0, size.x - 1);
        auto j     = clamp((int)(texcoord.y * size.y), 0, size.y - 1);
        auto prob  = sample_discrete_cdf_pdf(cdf, j * size.x + i) / cdf.back();
        auto angle = (2 * pif / size.x) * (pif / size.y) *
                     sin(pif * (j + 0.5f) / size.y);
        pdf += prob / angle;
      } else {
        pdf += 1 / (4 * pif);
      }
    }
  }
  pdf *= sample_uniform_pdf(scene->lights.size());
  return pdf;
}

// Recursive path tracing.
static std::pair<vec3f, bool> trace_path(const scene_model* scene,
    const ray3f& ray_, rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance      = zero3f;
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = vector<volume_point>{};
  auto max_roughness = 0.0f;
  auto hit           = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // handle transmission if inside a volume
    auto in_volume = false;
    if (!volume_stack.empty()) {
      auto& point    = volume_stack.back();
      auto  distance = sample_transmittance(
          point.voldensity, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(point.voldensity, distance) /
                sample_transmittance_pdf(
                    point.voldensity, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto point = eval_point(scene, intersection, ray);

      // correct roughness
      if (params.nocaustics) {
        max_roughness   = max(point.roughness, max_roughness);
        point.roughness = max_roughness;
      }

      // handle opacity
      if (point.opacity < 1 && rand1f(rng) >= point.opacity) {
        ray = {point.position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }
      hit = true;

      // accumulate emission
      radiance += weight * eval_emission(point);

      // next direction
      if (point.roughness) {
        if (rand1f(rng) < 0.5f) {
          point.incoming = sample_brdf(point, rand1f(rng), rand2f(rng));
        } else {
          point.incoming = sample_lights(
              scene, point.position, rand1f(rng), rand1f(rng), rand2f(rng));
        }
        weight *= eval_brdfcos(point) /
                  (0.5f * sample_brdf_pdf(point) +
                      0.5f * sample_lights_pdf(
                                 scene, point.position, point.incoming));
      } else {
        point.incoming = sample_delta(point, rand1f(rng));
        weight *= eval_delta(point) / sample_delta_pdf(point);
      }

      // update volume stack
      if (has_volume(scene, intersection) &&
          dot(point.normal, point.outgoing) *
                  dot(point.normal, point.incoming) <
              0) {
        if (volume_stack.empty()) {
          auto volpoint = eval_volume(scene, intersection, ray);
          volume_stack.push_back(volpoint);
        } else {
          volume_stack.pop_back();
        }
      }

      // setup next iteration
      ray = {point.position, point.incoming};
    } else {
      // prepare shading point
      auto point     = volume_stack.back();
      point.outgoing = -ray.d;
      point.position = ray.o + ray.d * intersection.distance;

      // handle opacity
      hit = true;

      // accumulate emission
      radiance += weight * eval_volemission(point);

      // next direction
      if (rand1f(rng) < 0.5f) {
        point.incoming = sample_scattering(point, rand1f(rng), rand2f(rng));
      } else {
        point.incoming = sample_lights(
            scene, point.position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      weight *=
          eval_scattering(point) /
          (0.5f * sample_scattering_pdf(point) +
              0.5f * sample_lights_pdf(scene, point.position, point.incoming));

      // setup next iteration
      ray = {point.position, point.incoming};
    }

    // check weight
    if (weight == zero3f || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit};
}

// Recursive path tracing.
static std::pair<vec3f, bool> trace_naive(const scene_model* scene,
    const ray3f& ray_, rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // prepare shading point
    auto point = eval_point(scene, intersection, ray);

    // handle opacity
    if (point.opacity < 1 && rand1f(rng) >= point.opacity) {
      ray = {point.position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    radiance += weight * eval_emission(point);

    // next direction
    if (point.roughness) {
      point.incoming = sample_brdf(point, rand1f(rng), rand2f(rng));
      weight *= eval_brdfcos(point) / sample_brdf_pdf(point);
    } else {
      point.incoming = sample_delta(point, rand1f(rng));
      weight *= eval_delta(point) / sample_delta_pdf(point);
    }

    // check weight
    if (weight == zero3f || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // setup next iteration
    ray = {point.position, point.incoming};
  }

  return {radiance, hit};
}

// Eyelight for quick previewing.
static std::pair<vec3f, bool> trace_eyelight(const scene_model* scene,
    const ray3f& ray_, rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = false;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_scene_bvh(scene, ray);
    if (!intersection.hit) {
      radiance += weight * eval_environment(scene, ray);
      break;
    }

    // prepare shading point
    auto point = eval_point(scene, intersection, ray);

    // handle opacity
    if (point.opacity < 1 && rand1f(rng) >= point.opacity) {
      ray = {point.position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    point.incoming = point.outgoing;
    radiance += weight * eval_emission(point);

    // brdf * light
    radiance += weight * pif * eval_brdfcos(point);

    // continue path
    if (point.roughness) break;
    point.incoming = sample_delta(point, rand1f(rng));
    weight *= eval_delta(point) / sample_delta_pdf(point);
    if (weight == zero3f || !isfinite(weight)) break;

    // setup next iteration
    ray = {point.position, point.incoming};
  }

  return {radiance, hit};
}

// False color rendering
static std::pair<vec3f, bool> trace_falsecolor(const scene_model* scene,
    const ray3f& ray, rng_state& rng, const trace_params& params) {
  // intersect next point
  auto intersection = intersect_scene_bvh(scene, ray);
  if (!intersection.hit) {
    return {zero3f, false};
  }

  // prepare shading point
  auto point = eval_point(scene, intersection, ray);

  // hash color
  auto hashed_color = [](int id) {
    auto hashed = std::hash<int>()(id);
    auto rng    = make_rng(trace_default_seed, hashed);
    return pow(0.5f + 0.5f * rand3f(rng), 2.2f);
  };

  switch (params.falsecolor) {
    case trace_falsecolor_type::normal: return {point.normal * 0.5f + 0.5f, 1};
    case trace_falsecolor_type::frontfacing:
      return {
          dot(point.normal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0}, 1};
    case trace_falsecolor_type::gnormal:
      return {point.gnormal * 0.5f + 0.5f, 1};
    case trace_falsecolor_type::gfrontfacing:
      return {
          dot(point.gnormal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0}, 1};
    case trace_falsecolor_type::texcoord:
      return {
          {fmod(point.texcoord.x, 1.0f), fmod(point.texcoord.y, 1.0f), 0}, 1};
    case trace_falsecolor_type::color: return {point.color, 1};
    case trace_falsecolor_type::emission: return {point.emission, 1};
    case trace_falsecolor_type::diffuse: return {point.diffuse, 1};
    case trace_falsecolor_type::specular: return {point.specular, 1};
    case trace_falsecolor_type::coat: return {point.coat, 1};
    case trace_falsecolor_type::metal: return {point.metal, 1};
    case trace_falsecolor_type::transmission: return {point.transmission, 1};
    case trace_falsecolor_type::translucency: return {point.translucency, 1};
    case trace_falsecolor_type::refraction: return {point.refraction, 1};
    case trace_falsecolor_type::roughness: return {vec3f{point.roughness}, 1};
    case trace_falsecolor_type::opacity: return {vec3f{point.opacity}, 1};
    case trace_falsecolor_type::ior: return {vec3f{point.ior}, 1};
    case trace_falsecolor_type::element:
      return {hashed_color(intersection.element), 1};
    case trace_falsecolor_type::object:
      return {hashed_color(intersection.object), 1};
    case trace_falsecolor_type::highlight: {
      auto emission = point.emission;
      if (emission == zero3f) emission = {0.2f, 0.2f, 0.2f};
      return {emission * abs(dot(-ray.d, point.normal)), 1};
    } break;
    default: return {zero3f, false};
  }
}

// Trace a single ray from the camera using the given algorithm.
using sampler_func = std::pair<vec3f, bool> (*)(const scene_model* scene,
    const ray3f& ray, rng_state& rng, const trace_params& params);
static sampler_func get_trace_sampler_func(const trace_params& params) {
  switch (params.sampler) {
    case trace_sampler_type::path: return trace_path;
    case trace_sampler_type::naive: return trace_naive;
    case trace_sampler_type::eyelight: return trace_eyelight;
    case trace_sampler_type::falsecolor: return trace_falsecolor;
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params) {
  switch (params.sampler) {
    case trace_sampler_type::path: return true;
    case trace_sampler_type::naive: return true;
    case trace_sampler_type::eyelight: return false;
    case trace_sampler_type::falsecolor: return false;
    default: {
      throw std::runtime_error("sampler unknown");
      return false;
    }
  }
}

// Trace a block of samples
vec4f trace_sample(trace_state* state, const scene_model* scene,
    const scene_camera* camera, const vec2i& ij, const trace_params& params) {
  auto  sampler = get_trace_sampler_func(params);
  auto& pixel   = state->pixels[ij];
  auto  ray = sample_camera(camera, ij, state->pixels.size(), rand2f(pixel.rng),
      rand2f(pixel.rng), params.tentfilter);
  auto [radiance, hit] = sampler(scene, ray, pixel.rng, params);
  if (!hit) {
    if (params.envhidden || scene->environments.empty()) {
      radiance = zero3f;
      hit      = false;
    } else {
      hit = true;
    }
  }
  if (!isfinite(radiance)) radiance = zero3f;
  if (max(radiance) > params.clamp)
    radiance = radiance * (params.clamp / max(radiance));
  pixel.radiance += radiance;
  pixel.hits += hit ? 1 : 0;
  pixel.samples += 1;
  return {pixel.hits ? pixel.radiance / pixel.hits : zero3f,
      (float)pixel.hits / (float)pixel.samples};
}

// Init a sequence of random number generators.
void init_state(trace_state* state, const scene_model* scene,
    const scene_camera* camera, const trace_params& params) {
  auto image_size = (camera->aspect >= 1)
                        ? vec2i{params.resolution,
                              (int)round(params.resolution / camera->aspect)}
                        : vec2i{(int)round(params.resolution * camera->aspect),
                              params.resolution};
  state->pixels.assign(image_size, trace_pixel{});
  state->render.assign(image_size, zero4f);
  auto rng = make_rng(1301081);
  for (auto& pixel : state->pixels) {
    pixel.rng = make_rng(params.seed, rand1i(rng, 1 << 31) / 2 + 1);
  }
}

// Build the bvh acceleration structure.
void init_bvh(scene_model* scene, const trace_params& params,
    progress_callback progress_cb) {
  auto params_       = scene_bvh_params{};
  params_.bvh        = (scene_bvh_type)params.bvh;
  params_.noparallel = params.noparallel;
  init_bvh(scene, params_, progress_cb);
}

// Refit bvh data
void update_bvh(scene_model*       scene,
    const vector<scene_object*>&   updated_objects,
    const vector<scene_shape*>&    updated_shapes,
    const vector<scene_instance*>& updated_instances,
    const trace_params&            params) {
  auto params_       = scene_bvh_params{};
  params_.bvh        = (scene_bvh_type)params.bvh;
  params_.noparallel = params.noparallel;
  update_bvh(
      scene, updated_objects, updated_shapes, updated_instances, params_);
}

// Forward declaration
static scene_light* add_light(scene_model* scene) {
  return scene->lights.emplace_back(new scene_light{});
}

// Init trace lights
void init_lights(scene_model* scene, progress_callback progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);

  for (auto light : scene->lights) delete light;
  scene->lights.clear();

  for (auto object : scene->objects) {
    if (object->material->emission == zero3f) continue;
    auto shape = object->shape;
    if (shape->triangles.empty() && shape->quads.empty()) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    if (!shape->triangles.empty()) {
      shape->elements_cdf = vector<float>(shape->triangles.size());
      for (auto idx = 0; idx < shape->elements_cdf.size(); idx++) {
        auto& t                  = shape->triangles[idx];
        shape->elements_cdf[idx] = triangle_area(shape->positions[t.x],
            shape->positions[t.y], shape->positions[t.z]);
        if (idx) shape->elements_cdf[idx] += shape->elements_cdf[idx - 1];
      }
    }
    if (!shape->quads.empty()) {
      shape->elements_cdf = vector<float>(shape->quads.size());
      for (auto idx = 0; idx < shape->elements_cdf.size(); idx++) {
        auto& t                  = shape->quads[idx];
        shape->elements_cdf[idx] = quad_area(shape->positions[t.x],
            shape->positions[t.y], shape->positions[t.z],
            shape->positions[t.w]);
        if (idx) shape->elements_cdf[idx] += shape->elements_cdf[idx - 1];
      }
    }
    for (auto iidx = 0; iidx < object->instance->frames.size(); iidx++) {
      auto light         = add_light(scene);
      light->object      = object;
      light->instance    = iidx;
      light->environment = nullptr;
    }
  }
  for (auto environment : scene->environments) {
    if (environment->emission == zero3f) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    if (environment->emission_tex) {
      auto texture            = environment->emission_tex;
      auto size               = texture_size(texture);
      environment->texels_cdf = vector<float>(size.x * size.y);
      if (size != zero2i) {
        for (auto i = 0; i < environment->texels_cdf.size(); i++) {
          auto ij                    = vec2i{i % size.x, i / size.x};
          auto th                    = (ij.y + 0.5f) * pif / size.y;
          auto value                 = lookup_texture(texture, ij);
          environment->texels_cdf[i] = max(value) * sin(th);
          if (i) environment->texels_cdf[i] += environment->texels_cdf[i - 1];
        }
      }
    }
    auto light         = add_light(scene);
    light->object      = nullptr;
    light->instance    = -1;
    light->environment = environment;
  }

  // handle progress
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto        futures  = vector<future<void>>{};
  auto        nthreads = std::thread::hardware_concurrency();
  atomic<int> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, size]() {
          while (true) {
            auto j = next_idx.fetch_add(1);
            if (j >= size.y) break;
            for (auto i = 0; i < size.x; i++) func({i, j});
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image(const scene_model* scene, const scene_camera* camera,
    const trace_params& params, progress_callback progress_cb,
    image_callback image_cb) {
  auto state_guard = std::make_unique<trace_state>();
  auto state       = state_guard.get();
  init_state(state, scene, camera, params);

  for (auto sample = 0; sample < params.samples; sample++) {
    if (progress_cb) progress_cb("trace image", sample, params.samples);
    if (params.noparallel) {
      for (auto j = 0; j < state->render.size().y; j++) {
        for (auto i = 0; i < state->render.size().x; i++) {
          state->render[{i, j}] = trace_sample(
              state, scene, camera, {i, j}, params);
        }
      }
    } else {
      parallel_for(state->render.size(),
          [state, scene, camera, &params](const vec2i& ij) {
            state->render[ij] = trace_sample(state, scene, camera, ij, params);
          });
    }
    if (image_cb) image_cb(state->render, sample + 1, params.samples);
  }

  if (progress_cb) progress_cb("trace image", params.samples, params.samples);
  return state->render;
}

// [experimental] Asynchronous interface
void trace_start(trace_state* state, const scene_model* scene,
    const scene_camera* camera, const trace_params& params,
    progress_callback progress_cb, image_callback image_cb,
    async_callback async_cb) {
  init_state(state, scene, camera, params);
  state->worker = {};
  state->stop   = false;

  // render preview
  if (progress_cb) progress_cb("trace preview", 0, params.samples);
  auto pprms = params;
  pprms.resolution /= params.pratio;
  pprms.samples = 1;
  auto preview  = trace_image(scene, camera, pprms);
  for (auto j = 0; j < state->render.size().y; j++) {
    for (auto i = 0; i < state->render.size().x; i++) {
      auto pi               = clamp(i / params.pratio, 0, preview.size().x - 1),
           pj               = clamp(j / params.pratio, 0, preview.size().y - 1);
      state->render[{i, j}] = preview[{pi, pj}];
    }
  }
  if (image_cb) image_cb(state->render, 0, params.samples);

  // start renderer
  state->worker = std::async(std::launch::async, [=]() {
    for (auto sample = 0; sample < params.samples; sample++) {
      if (state->stop) return;
      if (progress_cb) progress_cb("trace image", sample, params.samples);
      parallel_for(state->render.size(), [&](const vec2i& ij) {
        if (state->stop) return;
        state->render[ij] = trace_sample(state, scene, camera, ij, params);
        if (async_cb) async_cb(state->render, sample, params.samples, ij);
      });
      if (image_cb) image_cb(state->render, sample + 1, params.samples);
    }
    if (progress_cb) progress_cb("trace image", params.samples, params.samples);
    if (image_cb) image_cb(state->render, params.samples, params.samples);
  });
}
void trace_stop(trace_state* state) {
  if (!state) return;
  state->stop = true;
  if (state->worker.valid()) state->worker.get();
}

}  // namespace yocto
