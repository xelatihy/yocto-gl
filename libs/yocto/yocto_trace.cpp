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
#include <deque>
#include <future>
#include <memory>
#include <mutex>
using namespace std::string_literals;

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::trace {

// import math symbols for use
using math::abs;
using math::acos;
using math::atan2;
using math::clamp;
using math::cos;
using math::exp;
using math::flt_max;
using math::fmod;
using math::fresnel_conductor;
using math::fresnel_dielectric;
using math::identity3x3f;
using math::invalidb3f;
using math::log;
using math::make_rng;
using math::max;
using math::min;
using math::pif;
using math::pow;
using math::rng_state;
using math::sample_discrete;
using math::sample_discrete_pdf;
using math::sample_uniform;
using math::sample_uniform_pdf;
using math::sin;
using math::sqrt;
using math::zero2f;
using math::zero2i;
using math::zero3f;
using math::zero3i;
using math::zero4f;
using math::zero4i;

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::trace {

static std::pair<float, int> sample_distance(
    const vec3f& density, float rl, float rd) {
  auto channel         = clamp((int)(rl * 3), 0, 2);
  auto density_channel = density[channel];
  if (density_channel == 0 || rd == 0)
    return {flt_max, channel};
  else
    return {-log(rd) / density_channel, channel};
}

static float sample_distance_pdf(
    const vec3f& density, float distance, int channel) {
  auto density_channel = density[channel];
  return exp(-density_channel * distance);
}

static vec3f eval_transmission(const vec3f& density, float distance) {
  return exp(-density * distance);
}

static vec3f sample_phasefunction(float g, const vec2f& u) {
  auto cos_theta = 0.0f;
  if (abs(g) < 1e-3f) {
    cos_theta = 1 - 2 * u.x;
  } else {
    float square = (1 - g * g) / (1 - g + 2 * g * u.x);
    cos_theta    = (1 + g * g - square * square) / (2 * g);
  }

  auto sin_theta = sqrt(max(0.0f, 1 - cos_theta * cos_theta));
  auto phi       = 2 * pif * u.y;
  return {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
}

static float eval_phasefunction(float cos_theta, float g) {
  auto denom = 1 + g * g + 2 * g * cos_theta;
  return (1 - g * g) / (4 * pif * denom * sqrt(denom));
}

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto::trace {

// constant values
static const auto coat_ior       = 1.5;
static const auto coat_roughness = 0.03f * 0.03f;

// Shape element normal.
static vec3f eval_normal(const trc::shape* shape, int element) {
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
  } else {
    throw std::runtime_error("empty shape");
    norm = {0, 0, 1};
  }
  return norm;
}

// Shape element normal.
static std::pair<vec3f, vec3f> eval_tangents(
    const trc::shape* shape, int element, const vec2f& uv) {
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
static T eval_shape(const trc::shape* shape, const std::vector<T>& vals,
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
static vec2i texture_size(const trc::texture* texture) {
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
    const trc::texture* texture, const vec2i& ij, bool ldr_as_linear = false) {
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
static vec3f eval_texture(const trc::texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  // get texture
  if (!texture) return {1, 1, 1};

  // get img::image width/height
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

  // get img::image coordinates and residuals
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

// Generates a ray from a camera for img::image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_perspective_camera(
    const trc::camera* camera, const vec2f& image_uv, const vec2f& lens_uv) {
  // focal plane correction --- we skip it since we consider the lens to be
  // the effective focal lens
  // if (camera->focus < flt_max) {
  //   distance = camera->lens * camera->focus / (camera->focus - camera->lens);
  // }
  // point on the image plane
  auto q = vec3f{camera->film.x * (0.5f - image_uv.x),
      camera->film.y * (image_uv.y - 0.5f), camera->lens};
  // ray direction through the lens center
  auto dc = -normalize(q);
  // point on the lens
  auto e = vec3f{(lens_uv.x - 0.5f) * camera->aperture,
      (lens_uv.y - 0.5f) * camera->aperture, 0};
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
  // // distance of the img::image of the point
  // auto distance1 = camera->lens * distance / (distance - camera->lens);
  // auto q1        = -q * distance1 / distance;
  // auto d         = normalize(q1 - e);
  // // auto q1 = - normalize(q) * camera->focus / normalize(q).z;
  // auto ray = ray3f{transform_point(camera->frame, e),
  //     transform_direction(camera->frame, d)};
}

// Generates a ray from a camera for img::image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_orthographic_camera(
    const trc::camera* camera, const vec2f& image_uv, const vec2f& lens_uv) {
  // point on the image plane
  auto scale = 1 / camera->lens;
  auto q     = vec3f{camera->film.x * (0.5f - image_uv.x) * scale,
      camera->film.y * (image_uv.y - 0.5f) * scale, camera->lens};
  // point on the lens
  auto e = vec3f{-q.x, -q.y, 0} + vec3f{(lens_uv.x - 0.5f) * camera->aperture,
                                      (lens_uv.y - 0.5f) * camera->aperture, 0};
  // point on the focus plane
  auto p = vec3f{-q.x, -q.y, -camera->focus};
  // correct ray direction to account for camera focusing
  auto d = normalize(p - e);
  // done
  return ray3f{
      transform_point(camera->frame, e), transform_direction(camera->frame, d)};
}

// Generates a ray from a camera for img::image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_camera(
    const trc::camera* camera, const vec2f& uv, const vec2f& luv) {
  if (camera->orthographic)
    return eval_orthographic_camera(camera, uv, luv);
  else
    return eval_perspective_camera(camera, uv, luv);
}

// Sample camera
static ray3f sample_camera(const trc::camera* camera, const vec2i& ij,
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
  float refraction_pdf   = 0;
};

// Evaluate point
static trace_point eval_point(const trc::scene* scene,
    const intersection3f& intersection, const ray3f& ray) {
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
  point.refraction = thin ? zero3f : weight * transmission;
  weight *= 1 - (thin ? 0 : transmission);
  point.specular = weight * specular;
  weight *= 1 -
            specular * fresnel_dielectric(ior, point.outgoing, point.normal);
  point.transmission = weight * transmission * base;
  weight *= 1 - transmission;
  point.diffuse   = weight * base;
  point.meta      = reflectivity_to_eta(base);
  point.metak     = zero3f;
  point.roughness = roughness * roughness;
  point.ior       = ior;
  point.opacity   = opacity;

  // textures
  if (point.diffuse != zero3f || point.roughness) {
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
  point.refraction_pdf   = max(point.refraction);
  auto pdf_sum = point.diffuse_pdf + point.specular_pdf + point.metal_pdf +
                 point.coat_pdf + point.transmission_pdf + point.refraction_pdf;
  if (pdf_sum) {
    point.diffuse_pdf /= pdf_sum;
    point.specular_pdf /= pdf_sum;
    point.metal_pdf /= pdf_sum;
    point.coat_pdf /= pdf_sum;
    point.transmission_pdf /= pdf_sum;
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
static volume_point eval_volume(const trc::scene* scene,
    const intersection3f& intersection, const ray3f& ray) {
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
  auto thin       = material->thin || !material->transmission;
  auto scattering = material->scattering *
                    eval_texture(material->scattering_tex, texcoord, false).x;
  auto scanisotropy = material->scanisotropy;
  auto trdepth      = material->trdepth;

  // factors
  point.volemission = zero3f;
  point.voldensity  = (transmission && !thin)
                         ? -log(clamp(base, 0.0001f, 1.0f)) / trdepth
                         : zero3f;
  point.volscatter    = scattering;
  point.volanisotropy = scanisotropy;

  return point;
}

// Check if an instance as volume scattering
static bool has_volume(
    const trc::scene* scene, const intersection3f& intersection) {
  auto object = scene->objects[intersection.object];
  return !object->material->thin && object->material->transmission;
}

// Evaluate all environment color.
static vec3f eval_environment(const trc::scene* scene, const ray3f& ray) {
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

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto::trace {

#ifdef YOCTO_EMBREE
// Get Embree device
std::atomic<ssize_t> embree_memory = 0;
static RTCDevice     embree_device() {
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
static void init_embree_bvh(trc::shape* shape, const trace_params& params) {
  auto edevice = embree_device();
  if (shape->embree_bvh) rtcReleaseScene(shape->embree_bvh);
  shape->embree_bvh = rtcNewScene(edevice);
  auto escene       = shape->embree_bvh;
  if (params.bvh == bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == bvh_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  if (!shape->points.empty()) {
    throw std::runtime_error("embree does not support points");
  } else if (!shape->lines.empty()) {
    auto elines     = std::vector<int>{};
    auto epositions = std::vector<vec4f>{};
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
    if (params.bvh == bvh_type::embree_compact) {
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
    if (params.bvh == bvh_type::embree_compact) {
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

static void init_embree_bvh(trc::scene* scene, const trace_params& params) {
  // scene bvh
  auto edevice = embree_device();
  if (scene->embree_bvh) rtcReleaseScene(scene->embree_bvh);
  scene->embree_bvh = rtcNewScene(edevice);
  auto escene       = scene->embree_bvh;
  if (params.bvh == bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == bvh_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  auto object_id = 0;
  for (auto object : scene->objects) {
    auto instance_id = 0;
    for (auto frame : object->instance->frames) {
      frame          = frame * object->frame;
      auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_INSTANCE);
      rtcSetGeometryInstancedScene(egeometry, object->shape->embree_bvh);
      rtcSetGeometryTransform(
          egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &frame);
      rtcCommitGeometry(egeometry);
      rtcAttachGeometryByID(
          escene, egeometry, (int)scene->embree_instances.size());
      scene->embree_instances.push_back({object_id, instance_id});
      instance_id += 1;
    }
    object_id += 1;
  }
  rtcCommitScene(escene);
}

static void update_embree_bvh(trc::scene* scene,
    const std::vector<trc::object*>&      updated_objects,
    const std::vector<trc::shape*>&       updated_shapes,
    const std::vector<trc::instance*>&    updated_instances,
    const trace_params&                   params) {
  // scene bvh
  auto escene = scene->embree_bvh;
  for (auto& [object_id, instance_id] : scene->embree_instances) {
    auto object    = scene->objects[object_id];
    auto frame     = scene->objects[instance_id]->frame * object->frame;
    auto egeometry = rtcGetGeometry(escene, instance_id);
    rtcSetGeometryInstancedScene(egeometry, object->shape->embree_bvh);
    rtcSetGeometryTransform(
        egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &frame);
    rtcCommitGeometry(egeometry);
  }
  rtcCommitScene(escene);
}

static bool intersect_shape_embree_bvh(trc::shape* shape, const ray3f& ray,
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

static bool intersect_scene_embree_bvh(const trc::scene* scene,
    const ray3f& ray, int& shape, int& instance, int& element, vec2f& uv,
    float& distance, bool find_any) {
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
  shape    = scene->embree_instances[(int)embree_ray.hit.instID[0]].x;
  instance = scene->embree_instances[(int)embree_ray.hit.instID[0]].y;
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}
#endif

// primitive used to sort bvh entries
struct bvh_primitive {
  bbox3f bbox      = invalidb3f;
  vec3f  center    = zero3f;
  vec2i  primitive = {0, 0};
};

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
static std::pair<int, int> split_sah(
    std::vector<bvh_primitive>& primitives, int start, int end) {
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
    std::vector<bvh_primitive>& primitives, int start, int end) {
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

// Splits a BVH node using the middle heutirtic. Returns split position and
// axis.
static std::pair<int, int> split_middle(
    std::vector<bvh_primitive>& primitives, int start, int end) {
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
static std::pair<int, int> split_nodes(
    std::vector<bvh_primitive>& primitives, int start, int end, bvh_type type) {
  switch (type) {
    case bvh_type::default_: return split_middle(primitives, start, end);
    case bvh_type::highquality: return split_sah(primitives, start, end);
    case bvh_type::middle: return split_middle(primitives, start, end);
    case bvh_type::balanced: return split_balanced(primitives, start, end);
    default: throw std::runtime_error("should not have gotten here");
  }
}

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh_serial(std::vector<bvh_node>& nodes,
    std::vector<bvh_primitive>& primitives, bvh_type type) {
  // prepare to build nodes
  nodes.clear();
  nodes.reserve(primitives.size() * 2);

  // queue up first node
  auto queue = std::deque<vec3i>{{0, 0, (int)primitives.size()}};
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
    if (end - start > bvh_max_prims) {
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
    const shared_ptr<bvh_tree>& bvh, std::vector<bbox3f>& bboxes, bvh_type type) {
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
  auto centers = std::vector<vec3f>(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++)
    centers[idx] = center(bboxes[idx]);

  // queue up first node
  auto queue = std::deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // synchronization
  std::atomic<int>          num_processed_prims(0);
  std::mutex                queue_mutex;
  std::vector<std::future<void>> futures;
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
static void update_bvh(bvh_tree* bvh, const std::vector<bbox3f>& bboxes) {
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

static void init_bvh(trc::shape* shape, const trace_params& params) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (params.bvh == bvh_type::embree_default ||
      params.bvh == bvh_type::embree_highquality ||
      params.bvh == bvh_type::embree_compact) {
    return init_embree_bvh(shape, params);
  }
#endif

  // build primitives
  auto primitives = std::vector<bvh_primitive>{};
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < shape->points.size(); idx++) {
      auto& p             = shape->points[idx];
      auto& primitive     = primitives.emplace_back();
      primitive.bbox      = point_bounds(shape->positions[p], shape->radius[p]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {idx, 0};
    }
  } else if (!shape->lines.empty()) {
    for (auto idx = 0; idx < shape->lines.size(); idx++) {
      auto& l         = shape->lines[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {idx, 1};
    }
  } else if (!shape->triangles.empty()) {
    for (auto idx = 0; idx < shape->triangles.size(); idx++) {
      auto& primitive = primitives.emplace_back();
      auto& t         = shape->triangles[idx];
      primitive.bbox  = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {idx, 2};
    }
  } else if (!shape->quads.empty()) {
    for (auto idx = 0; idx < shape->quads.size(); idx++) {
      auto& q         = shape->quads[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {idx, 3};
    }
  }

  // build nodes
  if (shape->bvh) delete shape->bvh;
  shape->bvh = new bvh_tree{};
  build_bvh_serial(shape->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  shape->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    shape->bvh->primitives.push_back(primitive.primitive);
  }
}

void init_bvh(trc::scene* scene, const trace_params& params,
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
  if (params.bvh == bvh_type::embree_default ||
      params.bvh == bvh_type::embree_highquality ||
      params.bvh == bvh_type::embree_compact) {
    return init_embree_bvh(scene, params);
  }
#endif

  // handle progress
  if (progress_cb) progress_cb("build scene bvh", progress.x++, progress.y);

  // instance bboxes
  auto primitives            = std::vector<bvh_primitive>{};
  auto object_id             = 0;
  auto empty_instance_frames = std::vector<frame3f>{identity3x4f};
  for (auto object : scene->objects) {
    auto instance_id = 0;
    for (auto& frame : object->instance->frames) {
      auto& primitive = primitives.emplace_back();
      primitive.bbox  = object->shape->bvh->nodes.empty()
                           ? invalidb3f
                           : transform_bbox(frame * object->frame,
                                 object->shape->bvh->nodes[0].bbox);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {object_id, instance_id};
      instance_id += 1;
    }
    object_id += 1;
  }

  // build nodes
  if (scene->bvh) delete scene->bvh;
  scene->bvh = new bvh_tree{};
  build_bvh_serial(scene->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  scene->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    scene->bvh->primitives.push_back(primitive.primitive);
  }

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

static void update_bvh(trc::shape* shape, const trace_params& params) {
#ifdef YOCTO_EMBREE
  if (shape->embree_bvh) {
    throw std::runtime_error("embree shape update not implemented");
  }
#endif

  // build primitives
  auto bboxes = std::vector<bbox3f>(shape->bvh->primitives.size());
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& p     = shape->points[shape->bvh->primitives[idx].x];
      bboxes[idx] = point_bounds(shape->positions[p], shape->radius[p]);
    }
  } else if (!shape->lines.empty()) {
    bboxes = std::vector<bbox3f>(shape->lines.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& l     = shape->lines[shape->bvh->primitives[idx].x];
      bboxes[idx] = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
    }
  } else if (!shape->triangles.empty()) {
    bboxes = std::vector<bbox3f>(shape->triangles.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& t     = shape->triangles[shape->bvh->primitives[idx].x];
      bboxes[idx] = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
    }
  } else if (!shape->quads.empty()) {
    bboxes = std::vector<bbox3f>(shape->quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape->quads[shape->bvh->primitives[idx].x];
      bboxes[idx] = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
    }
  }

  // update nodes
  update_bvh(shape->bvh, bboxes);
}

void update_bvh(trc::scene*            scene,
    const std::vector<trc::object*>&   updated_objects,
    const std::vector<trc::shape*>&    updated_shapes,
    const std::vector<trc::instance*>& updated_instances,
    const trace_params&                params) {
  for (auto shape : updated_shapes) update_bvh(shape, params);

#ifdef YOCTO_EMBREE
  if (scene->embree_bvh) {
    update_embree_bvh(
        scene, updated_objects, updated_shapes, updated_instances, params);
  }
#endif

  // build primitives
  auto bboxes = std::vector<bbox3f>(scene->bvh->primitives.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto instance = scene->bvh->primitives[idx];
    auto object   = scene->objects[instance.x];
    auto sbvh     = object->shape->bvh;
    bboxes[idx]   = transform_bbox(
        object->instance->frames[instance.y] * object->frame,
        sbvh->nodes[0].bbox);
  }

  // update nodes
  update_bvh(scene->bvh, bboxes);
}

// Intersect ray with a bvh->
static bool intersect_shape_bvh(trc::shape* shape, const ray3f& ray_,
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
        auto& p = shape->points[shape->bvh->primitives[idx].x];
        if (intersect_point(
                ray, shape->positions[p], shape->radius[p], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx].x;
          ray.tmax = distance;
        }
      }
    } else if (!shape->lines.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& l = shape->lines[shape->bvh->primitives[idx].x];
        if (intersect_line(ray, shape->positions[l.x], shape->positions[l.y],
                shape->radius[l.x], shape->radius[l.y], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx].x;
          ray.tmax = distance;
        }
      }
    } else if (!shape->triangles.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& t = shape->triangles[shape->bvh->primitives[idx].x];
        if (intersect_triangle(ray, shape->positions[t.x],
                shape->positions[t.y], shape->positions[t.z], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx].x;
          ray.tmax = distance;
        }
      }
    } else if (!shape->quads.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape->quads[shape->bvh->primitives[idx].x];
        if (intersect_quad(ray, shape->positions[q.x], shape->positions[q.y],
                shape->positions[q.z], shape->positions[q.w], uv, distance)) {
          hit      = true;
          element  = shape->bvh->primitives[idx].x;
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
static bool intersect_scene_bvh(const trc::scene* scene, const ray3f& ray_,
    int& objecct, int& instance, int& element, vec2f& uv, float& distance,
    bool find_any, bool non_rigid_frames) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (scene->embree_bvh) {
    return intersect_scene_embree_bvh(
        scene, ray_, objecct, instance, element, uv, distance, find_any);
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
        auto [object_id, instance_id] = scene->bvh->primitives[idx];
        auto object                   = scene->objects[object_id];
        auto frame   = object->instance->frames[instance_id] * object->frame;
        auto inv_ray = transform_ray(inverse(frame, non_rigid_frames), ray);
        if (intersect_shape_bvh(
                object->shape, inv_ray, element, uv, distance, find_any)) {
          hit      = true;
          objecct  = object_id;
          instance = instance_id;
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
static bool intersect_instance_bvh(const trc::object* object, int instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto frame   = object->instance->frames[instance] * object->frame;
  auto inv_ray = transform_ray(inverse(frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      object->shape, inv_ray, element, uv, distance, find_any);
}

intersection3f intersect_scene_bvh(const trc::scene* scene, const ray3f& ray,
    bool find_any, bool non_rigid_frames) {
  auto intersection = intersection3f{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.object,
      intersection.instance, intersection.element, intersection.uv,
      intersection.distance, find_any, non_rigid_frames);
  return intersection;
}
intersection3f intersect_instance_bvh(const trc::object* object, int instance,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = intersection3f{};
  intersection.hit  = intersect_instance_bvh(object, instance, ray,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto::trace {

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
  return point.volscatter *
         eval_phasefunction(
             dot(point.outgoing, point.incoming), point.volanisotropy);
}

static vec3f sample_scattering(
    const volume_point& point, float rnl, const vec2f& rn) {
  if (point.voldensity == zero3f) return zero3f;
  auto direction = sample_phasefunction(point.volanisotropy, rn);
  return basis_fromz(-point.outgoing) * direction;
}

static float sample_scattering_pdf(const volume_point& point) {
  if (point.voldensity == zero3f) return 0;
  return eval_phasefunction(
      dot(point.outgoing, point.incoming), point.volanisotropy);
}

// Sample lights wrt solid angle
static vec3f sample_lights(const trc::scene* scene, const vec3f& position,
    float rl, float rel, const vec2f& ruv) {
  auto  light_id = sample_uniform(scene->lights.size(), rl);
  auto& light    = scene->lights[light_id];
  if (light->object) {
    auto& object    = light->object;
    auto  shape     = object->shape;
    auto  frame     = object->instance->frames[light->instance] * object->frame;
    auto  element   = sample_discrete(shape->elements_cdf, rel);
    auto  uv        = (!shape->triangles.empty()) ? sample_triangle(ruv) : ruv;
    auto  lposition = transform_point(
        frame, eval_shape(shape, shape->positions, element, uv, zero3f));
    return normalize(lposition - position);
  } else if (light->environment) {
    auto& environment = light->environment;
    if (environment->emission_tex) {
      auto emission_tex = environment->emission_tex;
      auto idx          = sample_discrete(environment->texels_cdf, rel);
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
    const trc::scene* scene, const vec3f& position, const vec3f& direction) {
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
        auto prob  = sample_discrete_pdf(cdf, j * size.x + i) / cdf.back();
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
static std::pair<vec3f, bool> trace_path(const trc::scene* scene,
    const ray3f& ray_, rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance      = zero3f;
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = std::vector<volume_point>{};
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
      auto& point              = volume_stack.back();
      auto [distance, channel] = sample_distance(
          point.voldensity, rand1f(rng), rand1f(rng));
      distance = min(distance, intersection.distance);
      weight *= eval_transmission(point.voldensity, distance) /
                sample_distance_pdf(point.voldensity, distance, channel);
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
    if (max(weight) < 1 && bounce > 6) {
      auto rr_prob = max((float)0.05, 1 - max(weight));
      if (rand1f(rng) > rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit};
}

// Recursive path tracing.
static std::pair<vec3f, bool> trace_naive(const trc::scene* scene,
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
    if (max(weight) < 1 && bounce > 3) {
      auto rr_prob = max((float)0.05, 1 - max(weight));
      if (rand1f(rng) > rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // setup next iteration
    ray = {point.position, point.incoming};
  }

  return {radiance, hit};
}

// Eyelight for quick previewing.
static std::pair<vec3f, bool> trace_eyelight(const trc::scene* scene,
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
static std::pair<vec3f, bool> trace_falsecolor(const trc::scene* scene,
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
    auto rng    = make_rng(default_seed, hashed);
    return pow(0.5f + 0.5f * rand3f(rng), 2.2f);
  };

  switch (params.falsecolor) {
    case falsecolor_type::normal: return {point.normal * 0.5f + 0.5f, 1};
    case falsecolor_type::frontfacing:
      return {
          dot(point.normal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0}, 1};
    case falsecolor_type::gnormal: return {point.gnormal * 0.5f + 0.5f, 1};
    case falsecolor_type::gfrontfacing:
      return {
          dot(point.gnormal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0}, 1};
    case falsecolor_type::texcoord:
      return {
          {fmod(point.texcoord.x, 1.0f), fmod(point.texcoord.y, 1.0f), 0}, 1};
    case falsecolor_type::color: return {point.color, 1};
    case falsecolor_type::emission: return {point.emission, 1};
    case falsecolor_type::diffuse: return {point.diffuse, 1};
    case falsecolor_type::specular: return {point.specular, 1};
    case falsecolor_type::coat: return {point.coat, 1};
    case falsecolor_type::metal: return {point.metal, 1};
    case falsecolor_type::transmission: return {point.transmission, 1};
    case falsecolor_type::refraction: return {point.refraction, 1};
    case falsecolor_type::roughness: return {vec3f{point.roughness}, 1};
    case falsecolor_type::opacity: return {vec3f{point.opacity}, 1};
    case falsecolor_type::element:
      return {hashed_color(intersection.element), 1};
    case falsecolor_type::object: return {hashed_color(intersection.object), 1};
    case falsecolor_type::highlight: {
      auto emission = point.emission;
      if (emission == zero3f) emission = {0.2f, 0.2f, 0.2f};
      return {emission * abs(dot(-ray.d, point.normal)), 1};
    } break;
    default: return {zero3f, false};
  }
}

// Trace a single ray from the camera using the given algorithm.
using sampler_func = std::pair<vec3f, bool> (*)(const trc::scene* scene,
    const ray3f& ray, rng_state& rng, const trace_params& params);
static sampler_func get_trace_sampler_func(const trace_params& params) {
  switch (params.sampler) {
    case sampler_type::path: return trace_path;
    case sampler_type::naive: return trace_naive;
    case sampler_type::eyelight: return trace_eyelight;
    case sampler_type::falsecolor: return trace_falsecolor;
    default: {
      throw std::runtime_error("sampler unknown");
      return nullptr;
    }
  }
}

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params) {
  switch (params.sampler) {
    case sampler_type::path: return true;
    case sampler_type::naive: return true;
    case sampler_type::eyelight: return false;
    case sampler_type::falsecolor: return false;
    default: {
      throw std::runtime_error("sampler unknown");
      return false;
    }
  }
}

// Trace a block of samples
vec4f trace_sample(trc::state* state, const trc::scene* scene,
    const trc::camera* camera, const vec2i& ij, const trace_params& params) {
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
void init_state(trc::state* state, const trc::scene* scene,
    const trc::camera* camera, const trace_params& params) {
  auto image_size =
      (camera->film.x > camera->film.y)
          ? vec2i{params.resolution,
                (int)round(params.resolution * camera->film.y / camera->film.x)}
          : vec2i{
                (int)round(params.resolution * camera->film.x / camera->film.y),
                params.resolution};
  state->pixels.assign(image_size, pixel{});
  state->render.assign(image_size, zero4f);
  auto rng = make_rng(1301081);
  for (auto& pixel : state->pixels) {
    pixel.rng = make_rng(params.seed, rand1i(rng, 1 << 31) / 2 + 1);
  }
}

// Forward declaration
trc::light* add_light(trc::scene* scene);

// Init trace lights
void init_lights(trc::scene* scene, progress_callback progress_cb) {
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
      shape->elements_cdf = std::vector<float>(shape->triangles.size());
      for (auto idx = 0; idx < shape->elements_cdf.size(); idx++) {
        auto& t                  = shape->triangles[idx];
        shape->elements_cdf[idx] = triangle_area(shape->positions[t.x],
            shape->positions[t.y], shape->positions[t.z]);
        if (idx) shape->elements_cdf[idx] += shape->elements_cdf[idx - 1];
      }
    }
    if (!shape->quads.empty()) {
      shape->elements_cdf = std::vector<float>(shape->quads.size());
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
      environment->texels_cdf = std::vector<float>(size.x * size.y);
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

using std::atomic;
using std::deque;
using std::future;

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(const vec2i& size, Func&& func) {
  auto             futures  = std::vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(0);
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
img::image<vec4f> trace_image(const trc::scene* scene,
    const trc::camera* camera, const trace_params& params,
    progress_callback progress_cb, image_callback image_cb) {
  auto state_guard = std::make_unique<state>();
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
void trace_start(trc::state* state, const trc::scene* scene,
    const trc::camera* camera, const trace_params& params,
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
      if (progress_cb) progress_cb("trace img::image", sample, params.samples);
      parallel_for(state->render.size(), [&](const vec2i& ij) {
        if (state->stop) return;
        state->render[ij] = trace_sample(state, scene, camera, ij, params);
        if (async_cb) async_cb(state->render, sample, params.samples, ij);
      });
      if (image_cb) image_cb(state->render, sample + 1, params.samples);
    }
    if (progress_cb)
      progress_cb("trace img::image", params.samples, params.samples);
    if (image_cb) image_cb(state->render, params.samples, params.samples);
  });
}
void trace_stop(trc::state* state) {
  if (!state) return;
  state->stop = true;
  if (state->worker.valid()) state->worker.get();
}

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto::trace {

// cleanup
shape::~shape() {
  if (bvh) delete bvh;
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene(embree_bvh);
#endif
}

// cleanup
scene::~scene() {
  if (bvh) delete bvh;
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene(embree_bvh);
#endif
  for (auto camera : cameras) delete camera;
  for (auto object : objects) delete object;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto instance : instances) delete instance;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
}

// Default instance
static auto default_instance = instance{{identity3x4f}};

// Add element
trc::camera* add_camera(trc::scene* scene) {
  return scene->cameras.emplace_back(new camera{});
}
trc::texture* add_texture(trc::scene* scene) {
  return scene->textures.emplace_back(new texture{});
}
trc::shape* add_shape(trc::scene* scene) {
  return scene->shapes.emplace_back(new shape{});
}
trc::material* add_material(trc::scene* scene) {
  return scene->materials.emplace_back(new material{});
}
trc::instance* add_instance(trc::scene* scene) {
  return scene->instances.emplace_back(new instance{});
}
trc::object* add_object(trc::scene* scene) {
  auto object_      = scene->objects.emplace_back(new object{});
  object_->instance = &default_instance;
  return object_;
}
trc::environment* add_environment(trc::scene* scene) {
  return scene->environments.emplace_back(new environment{});
}
trc::light* add_light(trc::scene* scene) {
  return scene->lights.emplace_back(new light{});
}

// Set cameras
void set_frame(trc::camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(
    trc::camera* camera, float lens, float aspect, float film, bool ortho) {
  camera->lens = lens;
  camera->film = aspect >= 1 ? vec2f{film, film / aspect}
                             : vec2f{film * aspect, film};
  camera->orthographic = ortho;
}
void set_focus(trc::camera* camera, float aperture, float focus) {
  camera->aperture = aperture;
  camera->focus    = focus;
}

// Add texture
void set_texture(trc::texture* texture, const img::image<vec3b>& img) {
  texture->colorb  = img;
  texture->colorf  = {};
  texture->scalarb = {};
  texture->scalarf = {};
}
void set_texture(trc::texture* texture, const img::image<vec3f>& img) {
  texture->colorb  = {};
  texture->colorf  = img;
  texture->scalarb = {};
  texture->scalarf = {};
}
void set_texture(trc::texture* texture, const img::image<byte>& img) {
  texture->colorb  = {};
  texture->colorf  = {};
  texture->scalarb = img;
  texture->scalarf = {};
}
void set_texture(trc::texture* texture, const img::image<float>& img) {
  texture->colorb  = {};
  texture->colorf  = {};
  texture->scalarb = {};
  texture->scalarf = img;
}

// Add shape
void set_points(trc::shape* shape, const std::vector<int>& points) {
  shape->points = points;
}
void set_lines(trc::shape* shape, const std::vector<vec2i>& lines) {
  shape->lines = lines;
}
void set_triangles(trc::shape* shape, const std::vector<vec3i>& triangles) {
  shape->triangles = triangles;
}
void set_quads(trc::shape* shape, const std::vector<vec4i>& quads) {
  shape->quads = quads;
}
void set_positions(trc::shape* shape, const std::vector<vec3f>& positions) {
  shape->positions = positions;
}
void set_normals(trc::shape* shape, const std::vector<vec3f>& normals) {
  shape->normals = normals;
}
void set_texcoords(trc::shape* shape, const std::vector<vec2f>& texcoords) {
  shape->texcoords = texcoords;
}
void set_colors(trc::shape* shape, const std::vector<vec3f>& colors) {
  shape->colors = colors;
}
void set_radius(trc::shape* shape, const std::vector<float>& radius) {
  shape->radius = radius;
}
void set_tangents(trc::shape* shape, const std::vector<vec4f>& tangents) {
  shape->tangents = tangents;
}

// Add object
void set_frame(trc::object* object, const frame3f& frame) {
  object->frame = frame;
}
void set_shape(trc::object* object, trc::shape* shape) {
  object->shape = shape;
}
void set_material(trc::object* object, trc::material* material) {
  object->material = material;
}
void set_instance(trc::object* object, trc::instance* instance) {
  object->instance = instance;
  if (!object->instance) object->instance = &default_instance;
}

// Add instance
void set_frames(trc::instance* instance, const std::vector<frame3f>& frames) {
  instance->frames = frames;
}

// Add material
void set_emission(trc::material* material, const vec3f& emission,
    trc::texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    trc::material* material, const vec3f& color, trc::texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    trc::material* material, float specular, trc::texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_metallic(
    trc::material* material, float metallic, trc::texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_ior(trc::material* material, float ior) { material->ior = ior; }
void set_transmission(trc::material* material, float transmission, bool thin,
    float trdepth, trc::texture* transmission_tex) {
  material->transmission     = transmission;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->transmission_tex = transmission_tex;
}
void set_thin(trc::material* material, bool thin) { material->thin = thin; }
void set_roughness(
    trc::material* material, float roughness, trc::texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    trc::material* material, float opacity, trc::texture* opacity_tex) {
  material->opacity     = opacity;
  material->opacity_tex = opacity_tex;
}
void set_scattering(trc::material* material, const vec3f& scattering,
    float scanisotropy, trc::texture* scattering_tex) {
  material->scattering     = scattering;
  material->scanisotropy   = scanisotropy;
  material->scattering_tex = scattering_tex;
}
void set_normalmap(trc::material* material, trc::texture* normal_tex) {
  material->normal_tex = normal_tex;
}

// Add environment
void set_frame(trc::environment* environment, const frame3f& frame) {
  environment->frame = frame;
}
void set_emission(trc::environment* environment, const vec3f& emission,
    trc::texture* emission_tex) {
  environment->emission     = emission;
  environment->emission_tex = emission_tex;
}

}  // namespace yocto::trace
