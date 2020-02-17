//
// Implementation for Yocto/Trace.
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

#include "yocto_trace.h"

#include <atomic>
#include <deque>
#include <future>
#include <memory>
#include <mutex>

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(const vec2f& ruv) {
  auto z   = ruv.y;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_pdf(const vec3f& direction) {
  return (direction.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(const vec3f& normal, const vec2f& ruv) {
  auto z               = ruv.y;
  auto r               = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec3f{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
inline float sample_hemisphere_pdf(
    const vec3f& normal, const vec3f& direction) {
  return (dot(normal, direction) <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(const vec2f& ruv) {
  auto z   = 2 * ruv.y - 1;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_sphere_pdf(const vec3f& w) { return 1 / (4 * pif); }

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(const vec2f& ruv) {
  auto z   = sqrt(ruv.y);
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cos_pdf(const vec3f& direction) {
  return (direction.z <= 0) ? 0 : direction.z / pif;
}

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(const vec3f& normal, const vec2f& ruv) {
  auto z               = sqrt(ruv.y);
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec3f{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
inline float sample_hemisphere_cos_pdf(
    const vec3f& normal, const vec3f& direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : cosw / pif;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float exponent, const vec2f& ruv) {
  auto z   = pow(ruv.y, 1 / (exponent + 1));
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cospower_pdf(
    float exponent, const vec3f& direction) {
  return (direction.z <= 0)
             ? 0
             : pow(direction.z, exponent) * (exponent + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
inline vec2f sample_disk(const vec2f& ruv) {
  auto r   = sqrt(ruv.y);
  auto phi = 2 * pif * ruv.x;
  return {cos(phi) * r, sin(phi) * r};
}
inline float sample_disk_pdf() { return 1 / pif; }

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(const vec2f& ruv) {
  auto phi = 2 * pif * ruv.x;
  return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
inline float sample_cylinder_pdf() { return 1 / pif; }

// Sample a point uniformly on a triangle returning the baricentric coordinates.
inline vec2f sample_triangle(const vec2f& ruv) {
  return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

// Sample a point uniformly on a triangle.
inline vec3f sample_triangle(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec2f& ruv) {
  auto uv = sample_triangle(ruv);
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
  return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
inline int sample_uniform(int size, float r) {
  return clamp((int)(r * size), 0, size - 1);
}
inline float sample_uniform_pdf(int size) { return (float)1 / (float)size; }

// Sample an index with uniform distribution.
inline float sample_uniform(const vector<float>& elements, float r) {
  if (elements.empty()) return {};
  auto size = (int)elements.size();
  return elements[clamp((int)(r * size), 0, size - 1)];
}
inline float sample_uniform_pdf(const vector<float>& elements) {
  if (elements.empty()) return 0;
  return 1.0f / (int)elements.size();
}

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const vector<float>& cdf, float r) {
  r        = clamp(r * cdf.back(), (float)0, cdf.back() - (float)0.00001);
  auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                   cdf.data());
  return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const vector<float>& cdf, int idx) {
  if (idx == 0) return cdf.at(0);
  return cdf.at(idx) - cdf.at(idx - 1);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& specular, float direction_cosine) {
  if (specular == zero3f) return zero3f;
  return specular + (1 - specular) *
                        pow(clamp(1 - abs(direction_cosine), 0.0f, 1.0f), 5.0f);
}

// Evaluates the GGX distribution and geometric term
float eval_microfacetD(float roughness, const vec3f& normal,
    const vec3f& half_vector, bool ggx = true) {
  auto cosine = dot(normal, half_vector);
  if (cosine <= 0) return 0;
  auto roughness_square = roughness * roughness;
  auto cosine_square    = cosine * cosine;
  auto tangent_square   = clamp(1 - cosine_square, 0.0f, 1.0f) / cosine_square;
  if (ggx) {
    return roughness_square / (pif * cosine_square * cosine_square *
                                  (roughness_square + tangent_square) *
                                  (roughness_square + tangent_square));
  } else {
    return exp(-tangent_square / roughness_square) /
           (pif * roughness_square * cosine_square * cosine_square);
  }
}
float evaluate_microfacetG1(float roughness, const vec3f& normal,
    const vec3f& half_vector, const vec3f& direction, bool ggx = true) {
  auto cosine = dot(normal, direction);
  if (dot(half_vector, direction) * cosine <= 0) return 0;
  auto roughness_square = roughness * roughness;
  auto cosine_square    = cosine * cosine;
  auto tangent_square   = clamp(1 - cosine_square, 0.0f, 1.0f) / cosine_square;
  if (ggx) {
    return 2 / (1 + sqrt(1.0f + roughness_square * tangent_square));
  } else {
    auto tangent       = sqrt(tangent_square);
    auto inv_rt        = 1 / (roughness * tangent);
    auto inv_rt_square = 1 / (roughness_square * tangent_square);
    if (inv_rt < 1.6f) {
      return (3.535f * inv_rt + 2.181f * inv_rt_square) /
             (1.0f + 2.276f * inv_rt + 2.577f * inv_rt_square);
    } else {
      return 1.0f;
    }
  }
}
float eval_microfacetG(float roughness, const vec3f& normal,
    const vec3f& half_vector, const vec3f& outgoing, const vec3f& incoming,
    bool ggx = true) {
  return evaluate_microfacetG1(roughness, normal, half_vector, outgoing, ggx) *
         evaluate_microfacetG1(roughness, normal, half_vector, incoming, ggx);
}
vec3f sample_microfacet(
    float roughness, const vec3f& normal, const vec2f& rn, bool ggx = true) {
  auto phi              = 2 * pif * rn.x;
  auto roughness_square = roughness * roughness;
  auto tangent_square   = 0.0f;
  if (ggx) {
    tangent_square = -roughness_square * log(1 - rn.y);
  } else {
    tangent_square = roughness_square * rn.y / (1 - rn.y);
  }
  auto cosine_square     = 1 / (1 + tangent_square);
  auto cosine            = 1 / sqrt(1 + tangent_square);
  auto radius            = sqrt(clamp(1 - cosine_square, 0.0f, 1.0f));
  auto local_half_vector = vec3f{cos(phi) * radius, sin(phi) * radius, cosine};
  return transform_direction(basis_fromz(normal), local_half_vector);
}
float sample_microfacet_pdf(float roughness, const vec3f& normal,
    const vec3f& half_vector, bool ggx = true) {
  auto cosine = dot(normal, half_vector);
  if (cosine < 0) return 0;
  return eval_microfacetD(roughness, normal, half_vector, ggx) * cosine;
}

// Specular to  eta.
vec3f reflectivity_to_eta(const vec3f& reflectivity_) {
  auto reflectivity = clamp(reflectivity_, 0.0f, 0.99f);
  return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
}

// Specular to fresnel eta.
pair<vec3f, vec3f> reflectivity_to_eta(
    const vec3f& reflectivity, const vec3f& edge_tint) {
  auto r = clamp(reflectivity, 0.0f, 0.99f);
  auto g = edge_tint;

  auto r_sqrt = sqrt(r);
  auto n_min  = (1 - r) / (1 + r);
  auto n_max  = (1 + r_sqrt) / (1 - r_sqrt);

  auto n  = lerp(n_max, n_min, g);
  auto k2 = ((n + 1) * (n + 1) * r - (n - 1) * (n - 1)) / (1 - r);
  k2      = max(k2, 0.0f);
  auto k  = sqrt(k2);
  return {n, k};
}

vec3f eta_to_reflectivity(const vec3f& eta) {
  return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
}
vec3f eta_to_reflectivity(const vec3f& eta, const vec3f& etak) {
  return ((eta - 1) * (eta - 1) + etak * etak) /
         ((eta + 1) * (eta + 1) + etak * etak);
}
vec3f eta_to_edge_tint(const vec3f& eta, const vec3f& etak) {
  auto r     = eta_to_reflectivity(eta, etak);
  auto numer = (1 + sqrt(r)) / (1 - sqrt(r)) - eta;
  auto denom = (1 + sqrt(r)) / (1 - sqrt(r)) - (1 - r) / (1 + r);
  return numer / denom;
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
float fresnel_dielectric(float eta, float cosw) {
  if (cosw < 0) {
    eta  = 1 / eta;
    cosw = -cosw;
  }

  auto sin2 = 1 - cosw * cosw;
  auto eta2 = eta * eta;

  auto cos2t = 1 - sin2 / eta2;
  if (cos2t < 0) return 1;  // tir

  auto t0 = sqrt(cos2t);
  auto t1 = eta * t0;
  auto t2 = eta * cosw;

  auto rs = (cosw - t1) / (cosw + t1);
  auto rp = (t0 - t2) / (t0 + t2);

  return (rs * rs + rp * rp) / 2;
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_dielectric(const vec3f& eta_, float cosw) {
  auto eta = eta_;
  if (cosw < 0) {
    eta  = 1 / eta;
    cosw = -cosw;
  }

  auto sin2 = 1 - cosw * cosw;
  auto eta2 = eta * eta;

  auto cos2t = 1 - sin2 / eta2;
  if (cos2t.x < 0 || cos2t.y < 0 || cos2t.z < 0) return {1, 1, 1};  // tir

  auto t0 = sqrt(cos2t);
  auto t1 = eta * t0;
  auto t2 = eta * cosw;

  auto rs = (cosw - t1) / (cosw + t1);
  auto rp = (t0 - t2) / (t0 + t2);

  return (rs * rs + rp * rp) / 2;
}

// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_conductor(const vec3f& eta, const vec3f& etak, float cosw) {
  if (cosw <= 0) return zero3f;
  if (etak == zero3f) return fresnel_dielectric(eta, cosw);

  cosw       = clamp(cosw, (float)-1, (float)1);
  auto cos2  = cosw * cosw;
  auto sin2  = clamp(1 - cos2, (float)0, (float)1);
  auto eta2  = eta * eta;
  auto etak2 = etak * etak;

  auto t0       = eta2 - etak2 - sin2;
  auto a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
  auto t1       = a2plusb2 + cos2;
  auto a        = sqrt((a2plusb2 + t0) / 2);
  auto t2       = 2 * a * cosw;
  auto rs       = (t1 - t2) / (t1 + t2);

  auto t3 = cos2 * a2plusb2 + sin2 * sin2;
  auto t4 = t2 * sin2;
  auto rp = rs * (t3 - t4) / (t3 + t4);

  return (rp + rs) / 2;
}

pair<float, int> sample_distance(const vec3f& density, float rl, float rd) {
  auto channel         = clamp((int)(rl * 3), 0, 2);
  auto density_channel = density[channel];
  if (density_channel == 0 || rd == 0)
    return {flt_max, channel};
  else
    return {-log(rd) / density_channel, channel};
}

float sample_distance_pdf(const vec3f& density, float distance, int channel) {
  auto density_channel = density[channel];
  return exp(-density_channel * distance);
}

vec3f eval_transmission(const vec3f& density, float distance) {
  return exp(-density * distance);
}

vec3f sample_phasefunction(float g, const vec2f& u) {
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

float eval_phasefunction(float cos_theta, float g) {
  auto denom = 1 + g * g + 2 * g * cos_theta;
  return (1 - g * g) / (4 * pif * denom * sqrt(denom));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SCENE EVALUATION
// -----------------------------------------------------------------------------
namespace yocto {

// constant values
static const auto coat_ior       = 1.5;
static const auto coat_roughness = 0.03f * 0.03f;

// Shape element normal.
static vec3f eval_element_normal(const trace_shape* shape, int element) {
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
  } else if (!shape->quadspos.empty()) {
    auto q = shape->quadspos[element];
    norm   = quad_normal(shape->positions[q.x], shape->positions[q.y],
        shape->positions[q.z], shape->positions[q.w]);
  } else {
    throw std::runtime_error("empty shape");
    norm = {0, 0, 1};
  }
  return norm;
}

// Shape element normal.
static pair<vec3f, vec3f> eval_element_tangents(
    const trace_shape* shape, int element, const vec2f& uv) {
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
  } else if (!shape->quadspos.empty()) {
    auto q = shape->quadspos[element];
    if (shape->texcoords.empty()) {
      return quad_tangents_fromuv(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w], {0, 0}, {1, 0}, {0, 1},
          {1, 1}, uv);
    } else {
      auto qt = shape->quadstexcoord[element];
      return quad_tangents_fromuv(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w], shape->texcoords[qt.x],
          shape->texcoords[qt.y], shape->texcoords[qt.z],
          shape->texcoords[qt.w], uv);
    }
  } else {
    return {zero3f, zero3f};
  }
}

// Shape value interpolated using barycentric coordinates
template <typename T>
static T eval_shape_elem(const trace_shape* shape,
    const vector<vec4i>& facevarying_quads, const vector<T>& vals, int element,
    const vec2f& uv) {
  if (vals.empty()) return {};
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
  } else if (!shape->quadspos.empty()) {
    auto q = facevarying_quads[element];
    if (q.w == q.z)
      return interpolate_triangle(vals[q.x], vals[q.y], vals[q.z], uv);
    return interpolate_quad(vals[q.x], vals[q.y], vals[q.z], vals[q.w], uv);
  } else {
    return {};
  }
}

// Check texture size
static vec2i texture_size(const trace_texture* texture) {
  if (!texture->hdr.empty()) {
    return texture->hdr.size();
  } else if (!texture->ldr.empty()) {
    return texture->ldr.size();
  } else {
    return zero2i;
  }
}

// Evaluate a texture
static vec3f lookup_texture(
    const trace_texture* texture, const vec2i& ij, bool ldr_as_linear = false) {
  if (texture->hdr.empty() && texture->ldr.empty()) return {1, 1, 1};
  if (!texture->hdr.empty()) {
    return xyz(texture->hdr[ij]);
  } else if (!texture->ldr.empty() && ldr_as_linear) {
    return xyz(byte_to_float(texture->ldr[ij]));
  } else if (!texture->ldr.empty() && !ldr_as_linear) {
    return xyz(srgb_to_rgb(byte_to_float(texture->ldr[ij])));
  } else {
    return {1, 1, 1};
  }
}

// Evaluate a texture
static vec3f eval_texture(const trace_texture* texture, const vec2f& uv,
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
static ray3f eval_perspective_camera(const trace_scene* scene, int camera_,
    const vec2f& image_uv, const vec2f& lens_uv) {
  auto& camera   = scene->cameras[camera_];
  auto  distance = camera->lens;
  if (camera->focus < flt_max) {
    distance = camera->lens * camera->focus / (camera->focus - camera->lens);
  }
  if (camera->aperture) {
    auto e = vec3f{(lens_uv.x - 0.5f) * camera->aperture,
        (lens_uv.y - 0.5f) * camera->aperture, 0};
    auto q = vec3f{camera->film.x * (0.5f - image_uv.x),
        camera->film.y * (image_uv.y - 0.5f), distance};
    // distance of the image of the point
    auto distance1 = camera->lens * distance / (distance - camera->lens);
    auto q1        = -q * distance1 / distance;
    auto d         = normalize(q1 - e);
    // auto q1 = - normalize(q) * camera->focus / normalize(q).z;
    auto ray = ray3f{transform_point(camera->frame, e),
        transform_direction(camera->frame, d)};
    return ray;
  } else {
    auto e   = zero3f;
    auto q   = vec3f{camera->film.x * (0.5f - image_uv.x),
        camera->film.y * (image_uv.y - 0.5f), distance};
    auto q1  = -q;
    auto d   = normalize(q1 - e);
    auto ray = ray3f{transform_point(camera->frame, e),
        transform_direction(camera->frame, d)};
    return ray;
  }
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_orthographic_camera(const trace_scene* scene, int camera_,
    const vec2f& image_uv, const vec2f& lens_uv) {
  auto& camera = scene->cameras[camera_];
  if (camera->aperture) {
    auto scale = 1 / camera->lens;
    auto q     = vec3f{camera->film.x * (0.5f - image_uv.x) * scale,
        camera->film.y * (image_uv.y - 0.5f) * scale, scale};
    auto q1    = vec3f{-q.x, -q.y, -camera->focus};
    auto e = vec3f{-q.x, -q.y, 0} + vec3f{(lens_uv.x - 0.5f) * camera->aperture,
                                        (lens_uv.y - 0.5f) * camera->aperture,
                                        0};
    auto d = normalize(q1 - e);
    auto ray = ray3f{transform_point(camera->frame, e),
        transform_direction(camera->frame, d)};
    return ray;
  } else {
    auto scale = 1 / camera->lens;
    auto q     = vec3f{camera->film.x * (0.5f - image_uv.x) * scale,
        camera->film.y * (image_uv.y - 0.5f) * scale, scale};
    auto q1    = -q;
    auto e     = vec3f{-q.x, -q.y, 0};
    auto d     = normalize(q1 - e);
    auto ray   = ray3f{transform_point(camera->frame, e),
        transform_direction(camera->frame, d)};
    return ray;
  }
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
static ray3f eval_camera(
    const trace_scene* scene, int camera_, const vec2f& uv, const vec2f& luv) {
  auto& camera = scene->cameras[camera_];
  if (camera->orthographic)
    return eval_orthographic_camera(scene, camera_, uv, luv);
  else
    return eval_perspective_camera(scene, camera_, uv, luv);
}

// Sample camera
static ray3f sample_camera(const trace_scene* scene, int camera,
    const vec2i& ij, const vec2i& image_size, const vec2f& puv,
    const vec2f& luv, bool tent) {
  if (!tent) {
    auto uv = vec2f{
        (ij.x + puv.x) / image_size.x, (ij.y + puv.y) / image_size.y};
    return eval_camera(scene, camera, uv, sample_disk(luv));
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
    return eval_camera(scene, camera, uv, sample_disk(luv));
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
  vec4f color    = zero4f;
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
static trace_point eval_point(const trace_scene* scene,
    const trace_intersection& intersection, const ray3f& ray) {
  // get data
  auto object   = scene->objects[intersection.object];
  auto shape    = object->shape;
  auto material = object->material;
  auto frame = object->instance->frames[intersection.instance] * object->frame;
  auto element                = intersection.element;
  auto uv                     = intersection.uv;
  auto trace_non_rigid_frames = true;

  // initialize point
  auto point     = trace_point{};
  point.object   = intersection.object;
  point.outgoing = -ray.d;
  point.incoming = -ray.d;

  // geometric properties
  point.position = eval_shape_elem(
      shape, shape->quadspos, shape->positions, element, uv);
  point.position = transform_point(frame, point.position);
  point.normal   = shape->normals.empty()
                     ? eval_element_normal(shape, element)
                     : normalize(eval_shape_elem(shape, shape->quadsnorm,
                           shape->normals, element, uv));
  point.normal  = transform_normal(frame, point.normal, trace_non_rigid_frames);
  point.gnormal = eval_element_normal(shape, element);
  point.gnormal = transform_normal(
      frame, point.gnormal, trace_non_rigid_frames);
  point.texcoord = shape->texcoords.empty()
                       ? uv
                       : eval_shape_elem(shape, shape->quadstexcoord,
                             shape->texcoords, element, uv);
  point.color = shape->colors.empty()
                    ? vec4f{1, 1, 1, 1}
                    : eval_shape_elem(shape, {}, shape->colors, element, uv);

  // correct normal
  if (!shape->points.empty()) {
    point.normal = point.outgoing;
  } else if (!shape->lines.empty()) {
    point.normal = orthonormalize(point.outgoing, point.normal);
  } else if (!material->normal_tex) {
    if (material->thin && dot(point.outgoing, point.normal) < 0)
      point.normal = -point.normal;
  } else {
    auto texcoord = shape->texcoords.empty()
                        ? uv
                        : eval_shape_elem(shape, shape->quadstexcoord,
                              shape->texcoords, element, uv);
    auto normalmap = -1 +
                     2 * eval_texture(material->normal_tex, texcoord, true);
    auto z = shape->normals.empty()
                 ? eval_element_normal(shape, element)
                 : normalize(eval_shape_elem(
                       shape, shape->quadsnorm, shape->normals, element, uv));
    auto basis  = identity3x3f;
    auto flip_v = false;
    if (shape->tangents.empty()) {
      auto tangents = eval_element_tangents(shape, element, uv);
      auto x        = orthonormalize(tangents.first, z);
      auto y        = normalize(cross(z, x));
      basis         = {x, y, z};
      flip_v        = dot(y, tangents.second) < 0;
    } else {
      auto tangsp = eval_shape_elem(shape, {}, shape->tangents, element, uv);
      auto x      = orthonormalize(xyz(tangsp), z);
      auto y      = normalize(cross(z, x));
      basis       = {x, y, z};
      flip_v      = tangsp.w < 0;
    }
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    auto normal  = normalize(basis * normalmap);
    point.normal = transform_normal(frame, normal, trace_non_rigid_frames);
    if (material->thin && dot(point.outgoing, point.normal) < 0)
      point.normal = -point.normal;
  }

  // material -------
  // initialize factors
  auto texcoord = point.texcoord;
  auto emission = material->emission *
                  eval_texture(material->emission_tex, texcoord);
  auto base = material->color * xyz(point.color) *
              eval_texture(material->color_tex, texcoord);
  auto specular = material->specular *
                  eval_texture(material->specular_tex, texcoord).x;
  auto metallic = material->metallic *
                  eval_texture(material->metallic_tex, texcoord).z;
  auto roughness = material->roughness *
                   eval_texture(material->roughness_tex, texcoord).x *
                   (material->gltf_textures
                           ? eval_texture(material->metallic_tex, texcoord).x
                           : 1);
  auto ior  = material->ior;
  auto coat = material->coat * eval_texture(material->coat_tex, texcoord).x;
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord).x;
  auto opacity = material->opacity * point.color.w *
                 mean(eval_texture(material->opacity_tex, texcoord, true));
  auto thin = material->thin || !material->transmission;

  // factors
  auto weight    = vec3f{1, 1, 1};
  point.emission = weight * emission;
  point.coat     = weight * coat;
  weight *= 1 - point.coat * fresnel_dielectric(
                                 coat_ior, dot(point.outgoing, point.normal));
  point.metal = weight * metallic;
  weight *= 1 - metallic;
  point.refraction = thin ? zero3f : weight * transmission;
  weight *= 1 - (thin ? 0 : transmission);
  point.specular = weight * specular;
  weight *= 1 - specular *
                    fresnel_dielectric(ior, dot(point.outgoing, point.normal));
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
  point.specular_pdf = max(
      point.specular *
      fresnel_dielectric(point.ior, dot(point.outgoing, point.normal)));
  point.metal_pdf = max(point.metal * fresnel_conductor(point.meta, point.metak,
                                          dot(point.outgoing, point.normal)));
  point.coat_pdf  = max(point.coat * fresnel_dielectric(coat_ior,
                                        dot(point.outgoing, point.normal)));
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
static volume_point eval_volume(const trace_scene* scene,
    const trace_intersection& intersection, const ray3f& ray) {
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
  point.position = eval_shape_elem(
      shape, shape->quadspos, shape->positions, element, uv);
  point.position = transform_point(frame, point.position);

  // material -------
  // initialize factors
  auto texcoord = shape->texcoords.empty()
                      ? uv
                      : eval_shape_elem(shape, shape->quadstexcoord,
                            shape->texcoords, element, uv);
  auto color = shape->colors.empty()
                   ? vec4f{1, 1, 1, 1}
                   : eval_shape_elem(shape, {}, shape->colors, element, uv);
  auto base = material->color * xyz(color) *
              eval_texture(material->color_tex, texcoord);
  auto transmission = material->transmission *
                      eval_texture(material->emission_tex, texcoord).x;
  auto thin       = material->thin || !material->transmission;
  auto scattering = material->scattering *
                    eval_texture(material->scattering_tex, texcoord).x;
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
    const trace_scene* scene, const trace_intersection& intersection) {
  auto object = scene->objects[intersection.object];
  return !object->material->thin && object->material->transmission;
}

// Evaluate all environment color.
static vec3f eval_environment(const trace_scene* scene, const ray3f& ray) {
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
// IMPLEMENRTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate)
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, vec2f& uv, float& dist) {
  // find parameter for line-point minimum distance
  auto w = p - ray.o;
  auto t = dot(w, ray.d) / dot(ray.d, ray.d);

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return false;

  // test for line-point distance vs point radius
  auto rp  = ray.o + ray.d * t;
  auto prp = p - rp;
  if (dot(prp, prp) > r * r) return false;

  // intersection occurred: set params and exit
  uv   = {0, 0};
  dist = t;
  return true;
}

// Intersect a ray with a line
inline bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, vec2f& uv, float& dist) {
  // setup intersection params
  auto u = ray.d;
  auto v = p1 - p0;
  auto w = ray.o - p0;

  // compute values to solve a linear system
  auto a   = dot(u, u);
  auto b   = dot(u, v);
  auto c   = dot(v, v);
  auto d   = dot(u, w);
  auto e   = dot(v, w);
  auto det = a * c - b * b;

  // check determinant and exit if lines are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return false;

  // compute Parameters on both ray and segment
  auto t = (b * e - c * d) / det;
  auto s = (a * e - b * d) / det;

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return false;

  // clamp segment param to segment corners
  s = clamp(s, (float)0, (float)1);

  // compute segment-segment distance on the closest points
  auto pr  = ray.o + ray.d * t;
  auto pl  = p0 + (p1 - p0) * s;
  auto prl = pr - pl;

  // check with the line radius at the same point
  auto d2 = dot(prl, prl);
  auto r  = r0 * (1 - s) + r1 * s;
  if (d2 > r * r) return {};

  // intersection occurred: set params and exit
  uv   = {s, sqrt(d2) / r};
  dist = t;
  return true;
}

// Intersect a ray with a triangle
inline bool intersect_triangle(const ray3f& ray, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, vec2f& uv, float& dist) {
  // compute triangle edges
  auto edge1 = p1 - p0;
  auto edge2 = p2 - p0;

  // compute determinant to solve a linear system
  auto pvec = cross(ray.d, edge2);
  auto det  = dot(edge1, pvec);

  // check determinant and exit if triangle and ray are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return false;
  auto inv_det = 1.0f / det;

  // compute and check first bricentric coordinated
  auto tvec = ray.o - p0;
  auto u    = dot(tvec, pvec) * inv_det;
  if (u < 0 || u > 1) return false;

  // compute and check second bricentric coordinated
  auto qvec = cross(tvec, edge1);
  auto v    = dot(ray.d, qvec) * inv_det;
  if (v < 0 || u + v > 1) return false;

  // compute and check ray parameter
  auto t = dot(edge2, qvec) * inv_det;
  if (t < ray.tmin || t > ray.tmax) return false;

  // intersection occurred: set params and exit
  uv   = {u, v};
  dist = t;
  return true;
}

// Intersect a ray with a quad.
inline bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, vec2f& uv, float& dist) {
  if (p2 == p3) {
    return intersect_triangle(ray, p0, p1, p3, uv, dist);
  }
  auto hit  = false;
  auto tray = ray;
  if (intersect_triangle(tray, p0, p1, p3, uv, dist)) {
    hit       = true;
    tray.tmax = dist;
  }
  if (intersect_triangle(tray, p2, p3, p1, uv, dist)) {
    hit       = true;
    uv        = 1 - uv;
    tray.tmax = dist;
  }
  return hit;
}

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(const ray3f& ray, const bbox3f& bbox) {
  // determine intersection ranges
  auto invd = 1.0f / ray.d;
  auto t0   = (bbox.min - ray.o) * invd;
  auto t1   = (bbox.max - ray.o) * invd;
  // flip based on range directions
  if (invd.x < 0.0f) swap(t0.x, t1.x);
  if (invd.y < 0.0f) swap(t0.y, t1.y);
  if (invd.z < 0.0f) swap(t0.z, t1.z);
  auto tmin = max(t0.z, max(t0.y, max(t0.x, ray.tmin)));
  auto tmax = min(t1.z, min(t1.y, min(t1.x, ray.tmax)));
  tmax *= 1.00000024f;  // for double: 1.0000000000000004
  return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
inline bool intersect_bbox(
    const ray3f& ray, const vec3f& ray_dinv, const bbox3f& bbox) {
  auto it_min = (bbox.min - ray.o) * ray_dinv;
  auto it_max = (bbox.max - ray.o) * ray_dinv;
  auto tmin   = min(it_min, it_max);
  auto tmax   = max(it_min, it_max);
  auto t0     = max(max(tmin), ray.tmin);
  auto t1     = min(min(tmax), ray.tmax);
  t1 *= 1.00000024f;  // for double: 1.0000000000000004
  return t0 <= t1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

#ifdef YOCTO_EMBREE
// Get Embree device
std::atomic<ssize_t> trace_embree_memory = 0;
static RTCDevice     trace_embree_device() {
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
          trace_embree_memory += bytes;
          return true;
        },
        nullptr);
  }
  return device;
}

// Initialize Embree BVH
static void init_embree_bvh(trace_shape* shape, const trace_params& params) {
  auto edevice = trace_embree_device();
  if (shape->embree_bvh) rtcReleaseScene(shape->embree_bvh);
  shape->embree_bvh = rtcNewScene(edevice);
  auto escene       = shape->embree_bvh;
  if (params.bvh == trace_bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == trace_bvh_type::embree_highquality)
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
    if (params.bvh == trace_bvh_type::embree_compact) {
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
    if (params.bvh == trace_bvh_type::embree_compact) {
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
  } else if (!shape->quadspos.empty()) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    if (params.bvh == trace_bvh_type::embree_compact) {
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape->positions.data(), 0, 3 * 4,
          shape->positions.size());
      rtcSetSharedGeometryBuffer(egeometry, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT4, shape->quadspos.data(), 0, 4 * 4,
          shape->quadspos.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape->positions.size());
      auto embree_quads     = rtcSetNewGeometryBuffer(egeometry,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4,
          shape->quadspos.size());
      memcpy(embree_positions, shape->positions.data(),
          shape->positions.size() * 12);
      memcpy(embree_quads, shape->quadspos.data(), shape->quadspos.size() * 16);
    }
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else {
    throw std::runtime_error("empty shapes not supported");
  }
  rtcCommitScene(escene);
}

static void init_embree_bvh(trace_scene* scene, const trace_params& params) {
  // scene bvh
  auto edevice = trace_embree_device();
  if (scene->embree_bvh) rtcReleaseScene(scene->embree_bvh);
  scene->embree_bvh = rtcNewScene(edevice);
  auto escene       = scene->embree_bvh;
  if (params.bvh == trace_bvh_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == trace_bvh_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  auto shape_id = 0;
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
      scene->embree_instances.push_back({shape_id, instance_id});
      instance_id += 1;
    }
    shape_id += 1;
  }
  rtcCommitScene(escene);
}

static void update_embree_bvh(
    trace_scene* scene, const vector<int>& updated_instances) {
  // scene bvh
  auto escene       = scene->embree_bvh;
  auto update_flags = vector<bool>(scene->shapes.size(), false);
  for (auto shape_id : updated_instances) update_flags[shape_id] = true;
  for (auto& [shape_id, instance_id] : scene->embree_instances) {
    if (!update_flags[shape_id]) continue;
    // FIXME: this is incorrect!!!
    auto shape     = scene->shapes[shape_id];
    auto frame     = scene->objects[instance_id]->frame;
    auto egeometry = rtcGetGeometry(escene, instance_id);
    rtcSetGeometryInstancedScene(egeometry, shape->embree_bvh);
    rtcSetGeometryTransform(
        egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &frame);
    rtcCommitGeometry(egeometry);
  }
  rtcCommitScene(escene);
}

static bool intersect_shape_embree_bvh(trace_shape* shape, const ray3f& ray,
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

static bool intersect_scene_embree_bvh(const trace_scene* scene,
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
struct trace_bvh_primitive {
  bbox3f bbox      = invalidb3f;
  vec3f  center    = zero3f;
  vec2i  primitive = {0, 0};
};

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
static pair<int, int> split_sah(
    vector<trace_bvh_primitive>& primitives, int start, int end) {
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
static pair<int, int> split_balanced(
    vector<trace_bvh_primitive>& primitives, int start, int end) {
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
    throw std::runtime_error("bad bvh split");
    axis = 0;
    mid  = (start + end) / 2;
  }

  return {mid, axis};
}

// Splits a BVH node using the middle heutirtic. Returns split position and
// axis.
static pair<int, int> split_middle(
    vector<trace_bvh_primitive>& primitives, int start, int end) {
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
    throw std::runtime_error("bad bvh split");
    axis = 0;
    mid  = (start + end) / 2;
  }

  return {mid, axis};
}

// Split bvh nodes according to a type
static pair<int, int> split_nodes(vector<trace_bvh_primitive>& primitives,
    int start, int end, trace_bvh_type type) {
  switch (type) {
    case trace_bvh_type::default_: return split_middle(primitives, start, end);
    case trace_bvh_type::highquality: return split_sah(primitives, start, end);
    case trace_bvh_type::middle: return split_middle(primitives, start, end);
    case trace_bvh_type::balanced:
      return split_balanced(primitives, start, end);
    default: throw std::runtime_error("should not have gotten here");
  }
}

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh_serial(vector<trace_bvh_node>& nodes,
    vector<trace_bvh_primitive>& primitives, trace_bvh_type type) {
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
    const shared_ptr<trace_bvh>& bvh, vector<bbox3f>& bboxes, trace_bvh_type type) {
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
  auto queue = std::deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // synchronization
  std::atomic<int>          num_processed_prims(0);
  std::mutex                queue_mutex;
  vector<std::future<void>> futures;
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
static void update_bvh(trace_bvh* bvh, const vector<bbox3f>& bboxes) {
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

static void init_bvh(trace_shape* shape, const trace_params& params) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (params.bvh == trace_bvh_type::embree_default ||
      params.bvh == trace_bvh_type::embree_highquality ||
      params.bvh == trace_bvh_type::embree_compact) {
    return init_embree_bvh(shape, params);
  }
#endif

  // build primitives
  auto primitives = vector<trace_bvh_primitive>{};
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
  } else if (!shape->quadspos.empty()) {
    for (auto idx = 0; idx < shape->quadspos.size(); idx++) {
      auto& q         = shape->quadspos[idx];
      auto& primitive = primitives.emplace_back();
      primitive.bbox = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
      primitive.center    = center(primitive.bbox);
      primitive.primitive = {idx, 4};
    }
  }

  // build nodes
  if (shape->bvh) delete shape->bvh;
  shape->bvh = new trace_bvh{};
  build_bvh_serial(shape->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  shape->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    shape->bvh->primitives.push_back(primitive.primitive);
  }
}

void init_bvh(trace_scene* scene, const trace_params& params,
    trace_progress progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1 + (int)scene->shapes.size()};

  // shapes
  for (auto idx = 0; idx < scene->shapes.size(); idx++) {
    if (progress_cb) progress_cb("build shape bvh", progress.x++, progress.y);
    init_bvh(scene->shapes[idx], params);
  }

  // embree
#ifdef YOCTO_EMBREE
  if (params.bvh == trace_bvh_type::embree_default ||
      params.bvh == trace_bvh_type::embree_highquality ||
      params.bvh == trace_bvh_type::embree_compact) {
    return init_embree_bvh(scene, params);
  }
#endif

  // handle progress
  if (progress_cb) progress_cb("build scene bvh", progress.x++, progress.y);

  // instance bboxes
  auto primitives            = vector<trace_bvh_primitive>{};
  auto object_id             = 0;
  auto empty_instance_frames = vector<frame3f>{identity3x4f};
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
  scene->bvh = new trace_bvh{};
  build_bvh_serial(scene->bvh->nodes, primitives, params.bvh);

  // set bvh primitives
  scene->bvh->primitives.reserve(primitives.size());
  for (auto& primitive : primitives) {
    scene->bvh->primitives.push_back(primitive.primitive);
  }

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

static void update_bvh(trace_shape* shape, const trace_params& params) {
#ifdef YOCTO_EMBREE
  if (shape->embree_bvh) {
    throw std::runtime_error("embree shape update not implemented");
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>(shape->bvh->primitives.size());
  if (!shape->points.empty()) {
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& p     = shape->points[shape->bvh->primitives[idx].x];
      bboxes[idx] = point_bounds(shape->positions[p], shape->radius[p]);
    }
  } else if (!shape->lines.empty()) {
    bboxes = vector<bbox3f>(shape->lines.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& l     = shape->lines[shape->bvh->primitives[idx].x];
      bboxes[idx] = line_bounds(shape->positions[l.x], shape->positions[l.y],
          shape->radius[l.x], shape->radius[l.y]);
    }
  } else if (!shape->triangles.empty()) {
    bboxes = vector<bbox3f>(shape->triangles.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& t     = shape->triangles[shape->bvh->primitives[idx].x];
      bboxes[idx] = triangle_bounds(
          shape->positions[t.x], shape->positions[t.y], shape->positions[t.z]);
    }
  } else if (!shape->quads.empty()) {
    bboxes = vector<bbox3f>(shape->quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape->quads[shape->bvh->primitives[idx].x];
      bboxes[idx] = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
    }
  } else if (!shape->quadspos.empty()) {
    bboxes = vector<bbox3f>(shape->quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape->quads[shape->bvh->primitives[idx].x];
      bboxes[idx] = quad_bounds(shape->positions[q.x], shape->positions[q.y],
          shape->positions[q.z], shape->positions[q.w]);
    }
  }

  // update nodes
  update_bvh(shape->bvh, bboxes);
}

void update_bvh(trace_scene* scene, const vector<int>& updated_instances,
    const vector<int>& updated_shapes, const trace_params& params) {
  // update shapes
  for (auto shape : updated_shapes) update_bvh(scene->shapes[shape], params);

#ifdef YOCTO_EMBREE
  if (scene->embree_bvh) {
    update_embree_bvh(scene, updated_instances);
  }
#endif

  // build primitives
  auto empty_instance_frames = vector<frame3f>{identity3x4f};
  auto bboxes                = vector<bbox3f>(scene->bvh->primitives.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& instance = scene->bvh->primitives[idx];
    auto  object   = scene->objects[instance.x];
    auto& sbvh     = scene->shapes[instance.x]->bvh;
    bboxes[idx]    = transform_bbox(
        object->instance->frames[instance.y] * object->frame,
        sbvh->nodes[0].bbox);
  }

  // update nodes
  update_bvh(scene->bvh, bboxes);
}

// Intersect ray with a bvh->
static bool intersect_shape_bvh(trace_shape* shape, const ray3f& ray_,
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
    } else if (!shape->quadspos.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape->quadspos[shape->bvh->primitives[idx].x];
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
static bool intersect_scene_bvh(const trace_scene* scene, const ray3f& ray_,
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
static bool intersect_instance_bvh(const trace_object* object, int instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto frame   = object->instance->frames[instance] * object->frame;
  auto inv_ray = transform_ray(inverse(frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      object->shape, inv_ray, element, uv, distance, find_any);
}

trace_intersection intersect_scene_bvh(const trace_scene* scene,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = trace_intersection{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.object,
      intersection.instance, intersection.element, intersection.uv,
      intersection.distance, find_any, non_rigid_frames);
  return intersection;
}
trace_intersection intersect_instance_bvh(const trace_object* object,
    int instance, const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = trace_intersection{};
  intersection.hit  = intersect_instance_bvh(object, instance, ray,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Set non-rigid frames as default
static const bool trace_non_rigid_frames = true;

static vec3f eval_emission(const trace_point& point) { return point.emission; }

static vec3f eval_volemission(const volume_point& point) {
  return point.volemission;
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_brdfcos(const trace_point& point) {
  if (!point.roughness) return zero3f;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;
  auto& incoming = point.incoming;

  auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
  auto same_hemi = dot(normal, outgoing) * dot(normal, incoming) > 0;

  auto brdfcos = zero3f;

  if (point.diffuse != zero3f && same_hemi) {
    brdfcos += point.diffuse / pif * abs(dot(normal, incoming));
  }

  if (point.specular != zero3f && point.refraction == zero3f && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    auto F       = fresnel_dielectric(point.ior, dot(halfway, outgoing));
    auto D       = eval_microfacetD(point.roughness, up_normal, halfway);
    auto G       = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, incoming);
    brdfcos += point.specular * F * D * G /
               abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming));
  }

  if (point.metal != zero3f && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    auto F = fresnel_conductor(point.meta, point.metak, dot(halfway, outgoing));
    auto D = eval_microfacetD(point.roughness, up_normal, halfway);
    auto G = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, incoming);
    brdfcos += point.metal * F * D * G /
               abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming));
  }

  if (point.coat != zero3f && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    auto F       = fresnel_dielectric(coat_ior, dot(halfway, outgoing));
    auto D       = eval_microfacetD(coat_roughness, up_normal, halfway);
    auto G       = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, incoming);
    brdfcos += point.coat * F * D * G /
               abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming));
  }

  if (point.transmission != zero3f && !same_hemi) {
    auto ir      = reflect(-incoming, up_normal);
    auto halfway = normalize(ir + outgoing);
    // auto F       = fresnel_schlick(
    //     point.reflectance, abs(dot(halfway, outgoing)), entering);
    auto D = eval_microfacetD(point.roughness, up_normal, halfway);
    auto G = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, ir);
    brdfcos += point.transmission * D * G /
               abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming));
  }

  if (point.refraction != zero3f && !same_hemi) {
    auto halfway_vector = dot(outgoing, normal) > 0
                              ? -(outgoing + point.ior * incoming)
                              : (point.ior * outgoing + incoming);
    auto halfway = normalize(halfway_vector);
    // auto F       = fresnel_dielectric(point.ior, dot(halfway, outgoing));
    auto F = fresnel_dielectric(point.ior, dot(normal, outgoing));
    auto D = eval_microfacetD(point.roughness, up_normal, halfway);
    auto G = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, incoming);

    auto dot_terms = (dot(outgoing, halfway) * dot(incoming, halfway)) /
                     (dot(outgoing, normal) * dot(incoming, normal));

    // [Walter 2007] equation 21
    brdfcos += point.refraction * abs(dot_terms) * (1 - F) * D * G /
               dot(halfway_vector, halfway_vector) * abs(dot(normal, incoming));
  }

  if (point.refraction != zero3f && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    auto F       = fresnel_dielectric(point.ior, dot(halfway, outgoing));
    auto D       = eval_microfacetD(point.roughness, up_normal, halfway);
    auto G       = eval_microfacetG(
        point.roughness, up_normal, halfway, outgoing, incoming);
    brdfcos += point.refraction * F * D * G /
               abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
               abs(dot(normal, incoming));
  }

  return brdfcos;
}

static vec3f eval_delta(const trace_point& point) {
  if (point.roughness) return zero3f;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;
  auto& incoming = point.incoming;

  auto same_hemi = dot(normal, outgoing) * dot(normal, incoming) > 0;

  auto brdfcos = zero3f;

  if (point.specular != zero3f && point.refraction == zero3f && same_hemi) {
    brdfcos += point.specular *
               fresnel_dielectric(point.ior, dot(normal, outgoing));
  }
  if (point.metal != zero3f && same_hemi) {
    brdfcos += point.metal * fresnel_conductor(point.meta, point.metak,
                                 dot(normal, outgoing));
  }
  if (point.coat != zero3f && same_hemi) {
    brdfcos += point.coat * fresnel_dielectric(coat_ior, dot(outgoing, normal));
  }
  if (point.transmission != zero3f && !same_hemi) {
    brdfcos += point.transmission;
  }
  if (point.refraction != zero3f && !same_hemi) {
    brdfcos += point.refraction *
               (1 - fresnel_dielectric(point.ior, dot(normal, outgoing)));
  }
  if (point.refraction != zero3f && same_hemi) {
    brdfcos += point.refraction *
               fresnel_dielectric(point.ior, dot(normal, outgoing));
  }

  return brdfcos;
}

// Picks a direction based on the BRDF
static vec3f sample_brdf(const trace_point& point, float rnl, const vec2f& rn) {
  if (!point.roughness) return zero3f;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;

  auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;

  auto cdf = 0.0f;

  if (point.diffuse_pdf) {
    cdf += point.diffuse_pdf;
    if (rnl < cdf) {
      return sample_hemisphere_cos(up_normal, rn);
    }
  }

  if (point.specular_pdf && !point.refraction_pdf) {
    cdf += point.specular_pdf;
    if (rnl < cdf) {
      auto halfway = sample_microfacet(point.roughness, up_normal, rn);
      return reflect(outgoing, halfway);
    }
  }

  if (point.metal_pdf) {
    cdf += point.metal_pdf;
    if (rnl < cdf) {
      auto halfway = sample_microfacet(point.roughness, up_normal, rn);
      return reflect(outgoing, halfway);
    }
  }

  if (point.coat_pdf) {
    cdf += point.coat_pdf;
    if (rnl < cdf) {
      auto halfway = sample_microfacet(coat_roughness, up_normal, rn);
      return reflect(outgoing, halfway);
    }
  }

  if (point.transmission_pdf) {
    cdf += point.transmission_pdf;
    if (rnl < cdf) {
      auto halfway = sample_microfacet(point.roughness, up_normal, rn);
      auto ir      = reflect(outgoing, halfway);
      return -reflect(ir, up_normal);
    }
  }

  if (point.refraction_pdf) {
    cdf += point.refraction_pdf;
    if (rnl < cdf) {
      // auto nrnl = (rnl - 1 + point.refraction_pdf) /
      // point.refraction_pdf;
      auto nrnl = rnl;
      if (nrnl < fresnel_dielectric(point.ior, dot(normal, outgoing))) {
        auto halfway = sample_microfacet(point.roughness, up_normal, rn);
        return reflect(outgoing, halfway);
      } else {
        auto halfway = sample_microfacet(point.roughness, up_normal, rn);
        return refract_notir(outgoing, halfway,
            dot(normal, outgoing) > 0 ? 1 / point.ior : point.ior);
      }
    }
  }

  return zero3f;
}

static vec3f sample_delta(const trace_point& point, float rnl) {
  if (point.roughness) return zero3f;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;

  auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;

  // keep a weight sum to pick a lobe
  auto cdf = 0.0f;
  cdf += point.diffuse_pdf;

  if (point.specular_pdf && !point.refraction_pdf) {
    cdf += point.specular_pdf;
    if (rnl < cdf) {
      return reflect(outgoing, up_normal);
    }
  }

  if (point.metal_pdf) {
    cdf += point.metal_pdf;
    if (rnl < cdf) {
      return reflect(outgoing, up_normal);
    }
  }

  if (point.coat_pdf) {
    cdf += point.coat_pdf;
    if (rnl < cdf) {
      return reflect(outgoing, up_normal);
    }
  }

  if (point.transmission_pdf) {
    cdf += point.transmission_pdf;
    if (rnl < cdf) {
      return -outgoing;
    }
  }

  if (point.refraction_pdf) {
    cdf += point.refraction_pdf;
    if (rnl < cdf) {
      // auto nrnl = (rnl - 1 + point.refraction_pdf) /
      // point.refraction_pdf;
      auto nrnl = rnl;
      if (nrnl < fresnel_dielectric(point.ior, dot(normal, outgoing))) {
        return reflect(outgoing, up_normal);
      } else {
        return refract_notir(outgoing, up_normal,
            dot(normal, outgoing) > 0 ? 1 / point.ior : point.ior);
      }
    }
  }

  return zero3f;
}

// Compute the weight for sampling the BRDF
static float sample_brdf_pdf(const trace_point& point) {
  if (!point.roughness) return 0;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;
  auto& incoming = point.incoming;

  auto up_normal = dot(normal, outgoing) >= 0 ? normal : -normal;
  auto same_hemi = dot(normal, outgoing) * dot(normal, incoming) > 0;

  auto pdf = 0.0f;

  if (point.diffuse_pdf && same_hemi) {
    pdf += point.diffuse_pdf * sample_hemisphere_cos_pdf(up_normal, incoming);
  }

  if (point.specular_pdf && !point.refraction_pdf && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    pdf += point.specular_pdf *
           sample_microfacet_pdf(point.roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  }

  if (point.metal_pdf && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    pdf += point.metal_pdf *
           sample_microfacet_pdf(point.roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  }

  if (point.coat_pdf && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    pdf += point.coat_pdf *
           sample_microfacet_pdf(coat_roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  }

  if (point.transmission_pdf && !same_hemi) {
    auto up_normal = dot(outgoing, normal) > 0 ? normal : -normal;
    auto ir        = reflect(-incoming, up_normal);
    auto halfway   = normalize(ir + outgoing);
    auto d         = sample_microfacet_pdf(point.roughness, up_normal, halfway);
    pdf += point.transmission_pdf * d / (4 * abs(dot(outgoing, halfway)));
  }

  if (point.refraction_pdf && !same_hemi) {
    auto halfway_vector = dot(outgoing, normal) > 0
                              ? -(outgoing + point.ior * incoming)
                              : (point.ior * outgoing + incoming);
    auto halfway = normalize(halfway_vector);
    // [Walter 2007] equation 17
    pdf += point.refraction_pdf *
           (1 - fresnel_dielectric(point.ior, dot(normal, outgoing))) *
           sample_microfacet_pdf(point.roughness, up_normal, halfway) *
           abs(dot(halfway, incoming)) / dot(halfway_vector, halfway_vector);
  }

  if (point.refraction_pdf && same_hemi) {
    auto halfway = normalize(incoming + outgoing);
    pdf += point.refraction_pdf *
           fresnel_dielectric(point.ior, dot(normal, outgoing)) *
           sample_microfacet_pdf(point.roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  }

  return pdf;
}

static float sample_delta_pdf(const trace_point& point) {
  if (point.roughness) return 0;

  auto& normal   = point.normal;
  auto& outgoing = point.outgoing;
  auto& incoming = point.incoming;

  auto same_hemi = dot(normal, outgoing) * dot(normal, incoming) > 0;

  auto pdf = 0.0f;
  if (point.specular_pdf && !point.refraction_pdf && same_hemi)
    pdf += point.specular_pdf;
  if (point.metal_pdf && same_hemi) pdf += point.metal_pdf;
  if (point.coat_pdf && same_hemi) pdf += point.coat_pdf;
  if (point.transmission_pdf && !same_hemi) pdf += point.transmission_pdf;
  if (point.refraction_pdf && !same_hemi)
    pdf += point.refraction_pdf *
           (1 - fresnel_dielectric(point.ior, dot(normal, outgoing)));
  if (point.refraction_pdf && same_hemi)
    pdf += point.refraction_pdf *
           fresnel_dielectric(point.ior, dot(normal, outgoing));
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
static vec3f sample_lights(const trace_scene* scene, const vec3f& position,
    float rl, float rel, const vec2f& ruv) {
  auto  light_id = sample_uniform(scene->lights.size(), rl);
  auto& light    = scene->lights[light_id];
  if (light->object) {
    auto& object    = light->object;
    auto  shape     = object->shape;
    auto  frame     = object->instance->frames[light->instance] * object->frame;
    auto  element   = sample_discrete(shape->elements_cdf, rel);
    auto  uv        = (!shape->triangles.empty()) ? sample_triangle(ruv) : ruv;
    auto  lposition = transform_point(frame,
        eval_shape_elem(shape, shape->quadspos, shape->positions, element, uv));
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
    const trace_scene* scene, const vec3f& position, const vec3f& direction) {
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
            frame, eval_shape_elem(object->shape, object->shape->quadspos,
                       object->shape->positions, intersection.element,
                       intersection.uv));
        auto lnormal = transform_normal(frame,
            eval_element_normal(object->shape, intersection.element),
            trace_non_rigid_frames);
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
static pair<vec3f, bool> trace_path(const trace_scene* scene, const ray3f& ray_,
    rng_state& rng, const trace_params& params) {
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
static pair<vec3f, bool> trace_naive(const trace_scene* scene,
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
static pair<vec3f, bool> trace_eyelight(const trace_scene* scene,
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
static pair<vec3f, bool> trace_falsecolor(const trace_scene* scene,
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
    case trace_falsecolor_type::color: return {xyz(point.color), 1};
    case trace_falsecolor_type::emission: return {point.emission, 1};
    case trace_falsecolor_type::diffuse: return {point.diffuse, 1};
    case trace_falsecolor_type::specular: return {point.specular, 1};
    case trace_falsecolor_type::coat: return {point.coat, 1};
    case trace_falsecolor_type::metal: return {point.metal, 1};
    case trace_falsecolor_type::transmission: return {point.transmission, 1};
    case trace_falsecolor_type::refraction: return {point.refraction, 1};
    case trace_falsecolor_type::roughness: return {vec3f{point.roughness}, 1};
    case trace_falsecolor_type::opacity: return {vec3f{point.opacity}, 1};
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
using trace_sampler_func = pair<vec3f, bool> (*)(const trace_scene* scene,
    const ray3f& ray, rng_state& rng, const trace_params& params);
static trace_sampler_func get_trace_sampler_func(const trace_params& params) {
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
vec4f trace_sample(trace_state* state, const trace_scene* scene,
    const vec2i& ij, const trace_params& params) {
  auto  sampler = get_trace_sampler_func(params);
  auto& pixel   = state->pixels[ij];
  auto  ray     = sample_camera(scene, params.camera, ij, state->pixels.size(),
      rand2f(pixel.rng), rand2f(pixel.rng), params.tentfilter);
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
void init_state(
    trace_state* state, const trace_scene* scene, const trace_params& params) {
  auto camera = scene->cameras[params.camera];
  auto image_size =
      (camera->film.x > camera->film.y)
          ? vec2i{params.resolution,
                (int)round(params.resolution * camera->film.y / camera->film.x)}
          : vec2i{
                (int)round(params.resolution * camera->film.x / camera->film.y),
                params.resolution};
  state->pixels.assign(image_size, trace_pixel{});
  state->render.assign(image_size, zero4f);
  auto rng = make_rng(1301081);
  for (auto& pixel : state->pixels) {
    pixel.rng = make_rng(params.seed, rand1i(rng, 1 << 31) / 2 + 1);
  }
}

// Init trace lights
void init_lights(trace_scene* scene, trace_progress progress_cb) {
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
      auto light         = scene->lights.emplace_back(new trace_light{});
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
    auto light         = scene->lights.emplace_back(new trace_light{});
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
  auto             futures  = vector<std::future<void>>{};
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
image<vec4f> trace_image(const trace_scene* scene, const trace_params& params,
    trace_progress progress_cb, trace_progress_image progress_image_cb) {
  auto state_guard = make_unique<trace_state>();
  auto state       = state_guard.get();
  init_state(state, scene, params);

  for (auto sample = 0; sample < params.samples; sample++) {
    if (progress_cb) progress_cb("trace image", sample, params.samples);
    if (params.noparallel) {
      for (auto j = 0; j < state->render.size().y; j++) {
        for (auto i = 0; i < state->render.size().x; i++) {
          state->render[{i, j}] = trace_sample(state, scene, {i, j}, params);
        }
      }
    } else {
      parallel_for(
          state->render.size(), [state, scene, &params](const vec2i& ij) {
            state->render[ij] = trace_sample(state, scene, ij, params);
          });
    }
    if (progress_image_cb)
      progress_image_cb(state->render, sample + 1, params.samples);
  }

  if (progress_cb) progress_cb("trace image", params.samples, params.samples);
  return state->render;
}

// [experimental] Asynchronous interface
unique_ptr<trace_state> trace_async_start(const trace_scene* scene,
    const trace_params& params, trace_progress progress_cb,
    trace_progress_image progress_image_cb,
    trace_process_async  progress_async_cb) {
  auto state = make_unique<trace_state>();
  trace_async_start(state.get(), scene, params, progress_cb, progress_image_cb,
      progress_async_cb);
  return state;
}
void trace_async_start(trace_state* state, const trace_scene* scene,
    const trace_params& params, trace_progress progress_cb,
    trace_progress_image progress_image_cb,
    trace_process_async  progress_async_cb) {
  init_state(state, scene, params);
  state->worker = {};
  state->stop   = false;

  // render preview
  if (progress_cb) progress_cb("trace preview", 0, params.samples);
  auto pprms = params;
  pprms.resolution /= params.pratio;
  pprms.samples = 1;
  auto preview  = trace_image(scene, pprms);
  for (auto j = 0; j < state->render.size().y; j++) {
    for (auto i = 0; i < state->render.size().x; i++) {
      auto pi               = clamp(i / params.pratio, 0, preview.size().x - 1),
           pj               = clamp(j / params.pratio, 0, preview.size().y - 1);
      state->render[{i, j}] = preview[{pi, pj}];
    }
  }
  if (progress_image_cb) progress_image_cb(state->render, 0, params.samples);

  // start renderer
  state->worker = std::async(std::launch::async, [=]() {
    for (auto sample = 0; sample < params.samples; sample++) {
      if (state->stop) return;
      if (progress_cb) progress_cb("trace image", sample, params.samples);
      parallel_for(state->render.size(), [&](const vec2i& ij) {
        if (state->stop) return;
        state->render[ij] = trace_sample(state, scene, ij, params);
        if (progress_async_cb)
          progress_async_cb(state->render, sample, params.samples, ij);
      });
      if (progress_image_cb)
        progress_image_cb(state->render, sample + 1, params.samples);
    }
    if (progress_cb) progress_cb("trace image", params.samples, params.samples);
    if (progress_image_cb)
      progress_image_cb(state->render, params.samples, params.samples);
  });
}
void trace_async_stop(trace_state* state) {
  if (!state) return;
  state->stop = true;
  if (state->worker.valid()) state->worker.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// cleanup
trace_shape::~trace_shape() {
  if (bvh) delete bvh;
  if (embree_bvh) rtcReleaseScene(embree_bvh);
}

// cleanup
trace_scene::~trace_scene() {
  for (auto camera : cameras) delete camera;
  for (auto object : objects) delete object;
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto instance : instances) delete instance;
  for (auto texture : textures) delete texture;
  for (auto environment : environments) delete environment;
  if (bvh) delete bvh;
  if (embree_bvh) rtcReleaseScene(embree_bvh);
}

// create scene
unique_ptr<trace_scene> make_trace_scene() {
  return make_unique<trace_scene>();
}

// Add cameras
trace_camera* add_camera(trace_scene* scene) {
  return scene->cameras.emplace_back(new trace_camera{});
}
void set_frame(trace_camera* camera, const frame3f& frame) {
  camera->frame = frame;
}
void set_lens(trace_camera* camera, float lens, float aspect, float film) {
  camera->lens = lens;
  camera->film = aspect >= 1 ? vec2f{film, film / aspect}
                             : vec2f{film * aspect, film};
}
void set_focus(trace_camera* camera, float aperture, float focus) {
  camera->aperture = aperture;
  camera->focus    = focus;
}

// Add texture
trace_texture* add_texture(trace_scene* scene) {
  return scene->textures.emplace_back(new trace_texture{});
}
void set_texture(trace_texture* texture, const image<vec4b>& img) {
  texture->ldr = img;
  texture->hdr = {};
}
void set_texture(trace_texture* texture, const image<vec4f>& img) {
  texture->ldr = {};
  texture->hdr = img;
}

// Add shape
trace_shape* add_shape(trace_scene* scene) {
  return scene->shapes.emplace_back(new trace_shape{});
}
void set_points(trace_shape* shape, const vector<int>& points) {
  shape->points = points;
}
void set_lines(trace_shape* shape, const vector<vec2i>& lines) {
  shape->lines = lines;
}
void set_triangles(trace_shape* shape, const vector<vec3i>& triangles) {
  shape->triangles = triangles;
}
void set_quads(trace_shape* shape, const vector<vec4i>& quads) {
  shape->quads = quads;
}
void set_fvquads(trace_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord) {
  shape->quadspos      = quadspos;
  shape->quadsnorm     = quadsnorm;
  shape->quadstexcoord = quadstexcoord;
}
void set_positions(trace_shape* shape, const vector<vec3f>& positions) {
  shape->positions = positions;
}
void set_normals(trace_shape* shape, const vector<vec3f>& normals) {
  shape->normals = normals;
}
void set_texcoords(trace_shape* shape, const vector<vec2f>& texcoords) {
  shape->texcoords = texcoords;
}
void set_colors(trace_shape* shape, const vector<vec4f>& colors) {
  shape->colors = colors;
}
void set_radius(trace_shape* shape, const vector<float>& radius) {
  shape->radius = radius;
}
void set_tangents(trace_shape* shape, const vector<vec4f>& tangents) {
  shape->tangents = tangents;
}

// Default instance
static auto default_instance = trace_instance{{identity3x4f}};

// Add object
trace_object* add_object(trace_scene* scene) {
  auto object      = scene->objects.emplace_back(new trace_object{});
  object->instance = &default_instance;
  return object;
}
void set_frame(trace_object* object, const frame3f& frame) {
  object->frame = frame;
}
void set_shape(trace_object* object, trace_shape* shape) {
  object->shape = shape;
}
void set_material(trace_object* object, trace_material* material) {
  object->material = material;
}
void set_instance(trace_object* object, trace_instance* instance) {
  object->instance = instance;
  if (!object->instance) object->instance = &default_instance;
}

// Add instance
trace_instance* add_instance(trace_scene* scene) {
  return scene->instances.emplace_back(new trace_instance{});
}
void set_frames(trace_instance* instance, const vector<frame3f>& frames) {
  instance->frames = frames;
}

// Add material
trace_material* add_material(trace_scene* scene) {
  return scene->materials.emplace_back(new trace_material{});
}
void set_emission(trace_material* material, const vec3f& emission,
    trace_texture* emission_tex) {
  material->emission     = emission;
  material->emission_tex = emission_tex;
}
void set_color(
    trace_material* material, const vec3f& color, trace_texture* color_tex) {
  material->color     = color;
  material->color_tex = color_tex;
}
void set_specular(
    trace_material* material, float specular, trace_texture* specular_tex) {
  material->specular     = specular;
  material->specular_tex = specular_tex;
}
void set_metallic(
    trace_material* material, float metallic, trace_texture* metallic_tex) {
  material->metallic     = metallic;
  material->metallic_tex = metallic_tex;
}
void set_ior(trace_material* material, float ior) { material->ior = ior; }
void set_transmission(trace_material* material, float transmission, bool thin,
    float trdepth, trace_texture* transmission_tex) {
  material->transmission     = transmission;
  material->thin             = thin;
  material->trdepth          = trdepth;
  material->transmission_tex = transmission_tex;
}
void set_thin(trace_material* material, bool thin) { material->thin = thin; }
void set_roughness(
    trace_material* material, float roughness, trace_texture* roughness_tex) {
  material->roughness     = roughness;
  material->roughness_tex = roughness_tex;
}
void set_opacity(
    trace_material* material, float opacity, trace_texture* opacity_tex) {
  material->opacity     = opacity;
  material->opacity_tex = opacity_tex;
}
void set_scattering(trace_material* material, const vec3f& scattering,
    float scanisotropy, trace_texture* scattering_tex) {
  material->scattering     = scattering;
  material->scanisotropy   = scanisotropy;
  material->scattering_tex = scattering_tex;
}
void set_normalmap(trace_material* material, trace_texture* normal_tex) {
  material->normal_tex = normal_tex;
}
void set_gltftextures(trace_material* material, bool gltf_textures) {
  material->gltf_textures = gltf_textures;
}

// Add environment
trace_environment* add_environment(trace_scene* scene) {
  return scene->environments.emplace_back(new trace_environment{});
}
void set_frame(trace_environment* environment, const frame3f& frame) {
  environment->frame = frame;
}
void set_emission(trace_environment* environment, const vec3f& emission,
    trace_texture* emission_tex) {
  environment->emission     = emission;
  environment->emission_tex = emission_tex;
}

}  // namespace yocto
