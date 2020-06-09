//
// # Yocto/Shading: Tiny library of shading functions for path tracing.
//
// Yocto/Shading provides shading and sampling functions useful to write path
// tracing shaders.
//
// 1. use `fresnel_dielectric()` or `fresnel_conductor()` to evaluate the
//    fresnel term for dielectrics or conductors; use `fresnel_schlick()` for
//    the Schlick fresnel approximation
// 2. use `eta_to_reflectivity()` and `reflective_to_eta()` to convert eta to
//    reflectivity and vice-versa; use `eta_to_edgetint()` and
//    `edgetint_to_eta()`
// 3. use `microfacet_distribution()` and `microfacet_shadowing()` to evaluate
//    a microfacet distribution and its associated shadowing term
// 4. evaluate BRDF lobes with
//    - `eval_diffuse_reflection()`: diffuse brdf
//    - `eval_microfacet_reflection()`: specular brdf for dielectrics and metals
//    - `eval_microfacet_transmission()`: transmission brdf for thin dielectrics
//    - `eval_microfacet_refraction()`: refraction brdf for dielectrics
// 5. sample BRDF lobes with `sample_XXX()` using the above lobe names
// 6. compute the PDF for BRDF lobe sampling with `sample_XXX_pdf()` using the
//    above lobe names
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

#ifndef _YOCTO_SHADING_H_
#define _YOCTO_SHADING_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_sampling.h"

// -----------------------------------------------------------------------------
// SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Schlick approximation of the Fresnel term.
inline vec3f fresnel_schlick(
    const vec3f& specular, const vec3f& normal, const vec3f& outgoing);
// Compute the fresnel term for dielectrics.
inline float fresnel_dielectric(
    float eta, const vec3f& normal, const vec3f& outgoing);
// Compute the fresnel term for metals.
inline vec3f fresnel_conductor(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing);

// Convert eta to reflectivity
inline vec3f eta_to_reflectivity(const vec3f& eta);
// Convert reflectivity to  eta.
inline vec3f reflectivity_to_eta(const vec3f& reflectivity);
// Convert conductor eta to reflectivity.
inline vec3f eta_to_reflectivity(const vec3f& eta, const vec3f& etak);
// Convert eta to edge tint parametrization.
inline std::pair<vec3f, vec3f> eta_to_edgetint(
    const vec3f& eta, const vec3f& etak);
// Convert reflectivity and edge tint to eta.
inline std::pair<vec3f, vec3f> edgetint_to_eta(
    const vec3f& reflectivity, const vec3f& edgetint);

// Evaluates the microfacet distribution.
inline float microfacet_distribution(float roughness, const vec3f& normal,
    const vec3f& halfway, bool ggx = true);
// Evaluates the microfacet shadowing.
inline float microfacet_shadowing(float roughness, const vec3f& normal,
    const vec3f& halfway, const vec3f& outgoing, const vec3f& incoming,
    bool ggx = true);

// Samples a microfacet distribution.
inline vec3f sample_microfacet(
    float roughness, const vec3f& normal, const vec2f& rn, bool ggx = true);
// Pdf for microfacet distribution sampling.
inline float sample_microfacet_pdf(float roughness, const vec3f& normal,
    const vec3f& halfway, bool ggx = true);

// Samples a microfacet distribution with the distribution of visible normals.
inline vec3f sample_microfacet(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn, bool ggx = true);
// Pdf for microfacet distribution sampling with the distribution of visible
// normals.
inline float sample_microfacet_pdf(float roughness, const vec3f& normal,
    const vec3f& halfway, const vec3f& outgoing, bool ggx = true);

// Evaluates a diffuse BRDF lobe.
inline vec3f eval_diffuse_reflection(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluate a translucent BRDF lobe.
inline vec3f eval_diffuse_transmission(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a specular BRDF lobe.
inline vec3f eval_microfacet_reflection(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a metal BRDF lobe.
inline vec3f eval_microfacet_reflection(const vec3f& eta, const vec3f& etak,
    float roughness, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming);
// Evaluates a transmission BRDF lobe.
inline vec3f eval_microfacet_transmission(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluates a refraction BRDF lobe.
inline vec3f eval_microfacet_refraction(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);

// Sample a diffuse BRDF lobe.
inline vec3f sample_diffuse_reflection(
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a translucency BRDF lobe.
inline vec3f sample_diffuse_transmission(
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a specular BRDF lobe.
inline vec3f sample_microfacet_reflection(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a metal BRDF lobe.
inline vec3f sample_microfacet_reflection(const vec3f& eta, const vec3f& etak,
    float roughness, const vec3f& normal, const vec3f& outgoing,
    const vec2f& rn);
// Sample a transmission BRDF lobe.
inline vec3f sample_microfacet_transmission(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn);
// Sample a refraction BRDF lobe.
inline vec3f sample_microfacet_refraction(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn);

// Pdf for diffuse BRDF lobe sampling.
inline float sample_diffuse_reflection_pdf(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for translucency BRDF lobe sampling.
inline float sample_diffuse_transmission_pdf(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for specular BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for metal BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(const vec3f& eta,
    const vec3f& etak, float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);
// Pdf for transmission BRDF lobe sampling.
inline float sample_microfacet_transmission_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for refraction BRDF lobe sampling.
inline float sample_microfacet_refraction_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);

// Evaluate a delta specular BRDF lobe.
inline vec3f eval_delta_reflection(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta metal BRDF lobe.
inline vec3f eval_delta_reflection(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta transmission BRDF lobe.
inline vec3f eval_delta_transmission(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);
// Evaluate a delta refraction BRDF lobe.
inline vec3f eval_delta_refraction(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);

// Sample a delta specular BRDF lobe.
inline vec3f sample_delta_reflection(
    float ior, const vec3f& normal, const vec3f& outgoing);
// Sample a delta metal BRDF lobe.
inline vec3f sample_delta_reflection(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing);
// Sample a delta transmission BRDF lobe.
inline vec3f sample_delta_transmission(
    float ior, const vec3f& normal, const vec3f& outgoing);
// Sample a delta refraction BRDF lobe.
inline vec3f sample_delta_refraction(
    float ior, const vec3f& normal, const vec3f& outgoing, float rnl);

// Pdf for delta specular BRDF lobe sampling.
inline float sample_delta_reflection_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta metal BRDF lobe sampling.
inline float sample_delta_reflection_pdf(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta transmission BRDF lobe sampling.
inline float sample_delta_transmission_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);
// Pdf for delta refraction BRDF lobe sampling.
inline float sample_delta_refraction_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming);

// Convert mean-free-path to transmission
inline vec3f mfp_to_transmission(const vec3f& mfp, float depth);

// Evaluate transmittance
inline vec3f eval_transmittance(const vec3f& density, float distance);
// Sample a distance proportionally to transmittance
inline float sample_transmittance(
    const vec3f& density, float max_distance, float rl, float rd);
// Pdf for distance sampling
inline float sample_transmittance_pdf(
    const vec3f& density, float distance, float max_distance);

// Evaluate phase function
inline float eval_phasefunction(
    float anisotropy, const vec3f& outgoing, const vec3f& incoming);
// Sample phase function
inline vec3f sample_phasefunction(
    float anisotropy, const vec3f& outgoing, const vec2f& rn);
// Pdf for phase function sampling
inline float sample_phasefunction_pdf(
    float anisotropy, const vec3f& outgoing, const vec3f& incoming);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Schlick approximation of the Fresnel term
inline vec3f fresnel_schlick(
    const vec3f& specular, const vec3f& normal, const vec3f& outgoing) {
  if (specular == zero3f) return zero3f;
  auto cosine = dot(normal, outgoing);
  return specular +
         (1 - specular) * pow(clamp(1 - abs(cosine), 0.0f, 1.0f), 5.0f);
}

// Compute the fresnel term for dielectrics.
inline float fresnel_dielectric(
    float eta, const vec3f& normal, const vec3f& outgoing) {
  // Implementation from
  // https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
  auto cosw = abs(dot(normal, outgoing));

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

// Compute the fresnel term for metals.
inline vec3f fresnel_conductor(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing) {
  // Implementation from
  // https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
  auto cosw = dot(normal, outgoing);
  if (cosw <= 0) return zero3f;

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

// Convert eta to reflectivity
inline vec3f eta_to_reflectivity(const vec3f& eta) {
  return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
}
// Convert reflectivity to  eta.
inline vec3f reflectivity_to_eta(const vec3f& reflectivity_) {
  auto reflectivity = clamp(reflectivity_, 0.0f, 0.99f);
  return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
}
// Convert conductor eta to reflectivity
inline vec3f eta_to_reflectivity(const vec3f& eta, const vec3f& etak) {
  return ((eta - 1) * (eta - 1) + etak * etak) /
         ((eta + 1) * (eta + 1) + etak * etak);
}
// Convert eta to edge tint parametrization
inline std::pair<vec3f, vec3f> eta_to_edgetint(
    const vec3f& eta, const vec3f& etak) {
  auto reflectivity = eta_to_reflectivity(eta, etak);
  auto numer        = (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity)) - eta;
  auto denom        = (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity)) -
               (1 - reflectivity) / (1 + reflectivity);
  auto edgetint = numer / denom;
  return {reflectivity, edgetint};
}
// Convert reflectivity and edge tint to eta.
inline std::pair<vec3f, vec3f> edgetint_to_eta(
    const vec3f& reflectivity, const vec3f& edgetint) {
  auto r = clamp(reflectivity, 0.0f, 0.99f);
  auto g = edgetint;

  auto r_sqrt = sqrt(r);
  auto n_min  = (1 - r) / (1 + r);
  auto n_max  = (1 + r_sqrt) / (1 - r_sqrt);

  auto n  = lerp(n_max, n_min, g);
  auto k2 = ((n + 1) * (n + 1) * r - (n - 1) * (n - 1)) / (1 - r);
  k2      = max(k2, 0.0f);
  auto k  = sqrt(k2);
  return {n, k};
}

// Evaluate microfacet distribution
inline float microfacet_distribution(
    float roughness, const vec3f& normal, const vec3f& halfway, bool ggx) {
  // https://google.github.io/filament/Filament.html#materialsystem/specularbrdf
  // http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html
  auto cosine = dot(normal, halfway);
  if (cosine <= 0) return 0;
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  if (ggx) {
    return roughness2 / (pif * (cosine2 * roughness2 + 1 - cosine2) *
                            (cosine2 * roughness2 + 1 - cosine2));
  } else {
    return exp((cosine2 - 1) / (roughness2 * cosine2)) /
           (pif * roughness2 * cosine2 * cosine2);
  }
}

// Evaluate the microfacet shadowing1
inline float microfacet_shadowing1(float roughness, const vec3f& normal,
    const vec3f& halfway, const vec3f& direction, bool ggx) {
  // https://google.github.io/filament/Filament.html#materialsystem/specularbrdf
  // http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html
  auto cosine  = dot(normal, direction);
  auto cosineh = dot(halfway, direction);
  if (cosine * cosineh <= 0) return 0;
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  if (ggx) {
    return 2 * abs(cosine) /
           (abs(cosine) + sqrt(cosine2 - roughness2 * cosine2 + roughness2));
  } else {
    auto ci = abs(cosine) / (roughness * sqrt(1 - cosine2));
    return ci < 1.6f ? (3.535f * ci + 2.181f * ci * ci) /
                           (1.0f + 2.276f * ci + 2.577f * ci * ci)
                     : 1.0f;
  }
}

// Evaluate microfacet shadowing
inline float microfacet_shadowing(float roughness, const vec3f& normal,
    const vec3f& halfway, const vec3f& outgoing, const vec3f& incoming,
    bool ggx) {
  return microfacet_shadowing1(roughness, normal, halfway, outgoing, ggx) *
         microfacet_shadowing1(roughness, normal, halfway, incoming, ggx);
}

// Sample a microfacet ditribution.
inline vec3f sample_microfacet(
    float roughness, const vec3f& normal, const vec2f& rn, bool ggx) {
  auto phi   = 2 * pif * rn.x;
  auto theta = 0.0f;
  if (ggx) {
    theta = atan(roughness * sqrt(rn.y / (1 - rn.y)));
  } else {
    auto roughness2 = roughness * roughness;
    theta           = atan(sqrt(-roughness2 * log(1 - rn.y)));
  }
  auto local_half_vector = vec3f{
      cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
  return transform_direction(basis_fromz(normal), local_half_vector);
}

// Pdf for microfacet distribution sampling.
inline float sample_microfacet_pdf(
    float roughness, const vec3f& normal, const vec3f& halfway, bool ggx) {
  auto cosine = dot(normal, halfway);
  if (cosine < 0) return 0;
  return microfacet_distribution(roughness, normal, halfway, ggx) * cosine;
}

// Sample a microfacet ditribution with the distribution of visible normals.
inline vec3f sample_microfacet(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn, bool ggx) {
  // http://jcgt.org/published/0007/04/01/
  if (ggx) {
    // move to local coordinate system
    auto basis   = basis_fromz(normal);
    auto Ve      = transform_direction(transpose(basis), outgoing);
    auto alpha_x = roughness, alpha_y = roughness;
    // Section 3.2: transforming the view direction to the hemisphere
    // configuration
    auto Vh = normalize(vec3f{alpha_x * Ve.x, alpha_y * Ve.y, Ve.z});
    // Section 4.1: orthonormal basis (with special case if cross product is
    // zero)
    auto lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    auto T1    = lensq > 0 ? vec3f{-Vh.y, Vh.x, 0} * (1 / sqrt(lensq))
                        : vec3f{1, 0, 0};
    auto T2 = cross(Vh, T1);
    // Section 4.2: parameterization of the projected area
    auto r   = sqrt(rn.y);
    auto phi = 2 * pif * rn.x;
    auto t1  = r * cos(phi);
    auto t2  = r * sin(phi);
    auto s   = 0.5f * (1 + Vh.z);
    t2       = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
    // Section 4.3: reprojection onto hemisphere
    auto Nh = t1 * T1 + t2 * T2 + sqrt(max(0.0f, 1 - t1 * t1 - t2 * t2)) * Vh;
    // Section 3.4: transforming the normal back to the ellipsoid configuration
    auto Ne = normalize(vec3f{alpha_x * Nh.x, alpha_y * Nh.y, max(0.0f, Nh.z)});
    // move to world coordinate
    auto local_halfway = Ne;
    return transform_direction(basis, local_halfway);
  } else {
    throw std::invalid_argument{"not implemented yet"};
  }
}

// Pdf for microfacet distribution sampling with the distribution of visible
// normals.
inline float sample_microfacet_pdf(float roughness, const vec3f& normal,
    const vec3f& halfway, const vec3f& outgoing, bool ggx) {
  // http://jcgt.org/published/0007/04/01/
  if (dot(normal, halfway) < 0) return 0;
  if (dot(halfway, outgoing) < 0) return 0;
  return microfacet_distribution(roughness, normal, halfway, ggx) *
         microfacet_shadowing1(roughness, normal, halfway, outgoing, ggx) *
         max(0.0f, dot(halfway, outgoing)) / abs(dot(normal, outgoing));
}

// Evaluate a diffuse BRDF lobe.
inline vec3f eval_diffuse_reflection(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return zero3f;
  return vec3f{1} / pif * dot(normal, incoming);
}

// Evaluate a translucent BRDF lobe.
inline vec3f eval_diffuse_transmission(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) return zero3f;
  return vec3f{1} / pif * abs(dot(normal, incoming));
}

// Evaluate a specular BRDF lobe.
inline vec3f eval_microfacet_reflection(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return zero3f;
  auto halfway = normalize(incoming + outgoing);
  auto F       = fresnel_dielectric(ior, halfway, incoming);
  auto D       = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(roughness, normal, halfway, outgoing, incoming);
  return vec3f{1} * F * D * G /
         (4 * dot(normal, outgoing) * dot(normal, incoming)) *
         dot(normal, incoming);
}

// Evaluate a metal BRDF lobe.
inline vec3f eval_microfacet_reflection(const vec3f& eta, const vec3f& etak,
    float roughness, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return zero3f;
  auto halfway = normalize(incoming + outgoing);
  auto F       = fresnel_conductor(eta, etak, halfway, incoming);
  auto D       = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(roughness, normal, halfway, outgoing, incoming);
  return F * D * G / (4 * dot(normal, outgoing) * dot(normal, incoming)) *
         dot(normal, incoming);
}

// Evaluate a transmission BRDF lobe.
inline vec3f eval_microfacet_transmission(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return zero3f;
  auto reflected = reflect(-incoming, normal);
  auto halfway   = normalize(reflected + outgoing);
  // auto F       = fresnel_schlick(
  //     point.reflectance, abs(dot(halfway, outgoing)), entering);
  auto D = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(
      roughness, normal, halfway, outgoing, reflected);
  return vec3f{1} * D * G /
         (4 * dot(normal, outgoing) * dot(normal, reflected)) *
         (dot(normal, reflected));
}

// Evaluate a refraction BRDF lobe.
inline vec3f eval_microfacet_refraction(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    auto halfway = normalize(incoming + outgoing);
    auto F       = fresnel_dielectric(rel_ior, halfway, outgoing);
    auto D       = microfacet_distribution(roughness, up_normal, halfway);
    auto G       = microfacet_shadowing(
        roughness, up_normal, halfway, outgoing, incoming);
    return vec3f{1} * F * D * G /
           abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
           abs(dot(normal, incoming));
  } else {
    auto halfway = -normalize(rel_ior * incoming + outgoing) *
                   (entering ? 1 : -1);
    auto F = fresnel_dielectric(rel_ior, halfway, outgoing);
    auto D = microfacet_distribution(roughness, up_normal, halfway);
    auto G = microfacet_shadowing(
        roughness, up_normal, halfway, outgoing, incoming);
    // [Walter 2007] equation 21
    return vec3f{1} *
           abs((dot(outgoing, halfway) * dot(incoming, halfway)) /
               (dot(outgoing, normal) * dot(incoming, normal))) *
           (1 - F) * D * G /
           pow(rel_ior * dot(halfway, incoming) + dot(halfway, outgoing), 2) *
           abs(dot(normal, incoming));
  }
}

// Sample a diffuse BRDF lobe.
inline vec3f sample_diffuse_reflection(
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  return sample_hemisphere_cos(normal, rn);
}

// Sample a translucency BRDF lobe.
inline vec3f sample_diffuse_transmission(
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
  auto up_normal = dot(normal, outgoing) >= 0 ? normal : -normal;
  return sample_hemisphere_cos(-up_normal, rn);
}

// Sample a specular BRDF lobe.
inline vec3f sample_microfacet_reflection(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  // auto halfway = sample_microfacet(roughness, normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, normal, rn);
  return reflect(outgoing, halfway);
}

// Sample a metal BRDF lobe.
inline vec3f sample_microfacet_reflection(const vec3f& eta, const vec3f& etak,
    float roughness, const vec3f& normal, const vec3f& outgoing,
    const vec2f& rn) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  // auto halfway = sample_microfacet(roughness, normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, normal, rn);
  return reflect(outgoing, halfway);
}

// Sample a transmission BRDF lobe.
inline vec3f sample_microfacet_transmission(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  auto halfway = sample_microfacet(roughness, normal, rn);
  // auto halfway   = sample_microfacet(roughness, normal, outgoing, rn);
  auto reflected = reflect(outgoing, halfway);
  return -reflect(reflected, normal);
}

// Sample a refraction BRDF lobe.
inline vec3f sample_microfacet_refraction(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  // auto halfway   = sample_microfacet(roughness, up_normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, up_normal, rn);
  if (rnl < fresnel_dielectric(entering ? ior : (1 / ior), halfway, outgoing)) {
    return reflect(outgoing, halfway);
  } else {
    return refract(outgoing, halfway, entering ? (1 / ior) : ior);
  }
}

// Pdf for diffuse BRDF lobe sampling.
inline float sample_diffuse_reflection_pdf(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return sample_hemisphere_cos_pdf(normal, incoming);
}

// Pdf for translucency BRDF lobe sampling.
inline float sample_diffuse_transmission_pdf(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) return 0;
  auto up_normal = dot(normal, outgoing) >= 0 ? normal : -normal;
  return sample_hemisphere_cos_pdf(-up_normal, incoming);
}

// Pdf for specular BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  auto halfway = normalize(outgoing + incoming);
  // return sample_microfacet_pdf(roughness, normal, halfway, outgoing) /
  return sample_microfacet_pdf(roughness, normal, halfway) /
         (4 * abs(dot(outgoing, halfway)));
}

// Pdf for metal BRDF lobe sampling.
inline float sample_microfacet_reflection_pdf(const vec3f& eta,
    const vec3f& etak, float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  auto halfway = normalize(outgoing + incoming);
  // return sample_microfacet_pdf(roughness, normal, halfway, outgoing) /
  return sample_microfacet_pdf(roughness, normal, halfway) /
         (4 * abs(dot(outgoing, halfway)));
}

// Pdf for transmission BRDF lobe sampling.
inline float sample_microfacet_transmission_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return 0;
  auto reflected = reflect(-incoming, normal);
  auto halfway   = normalize(reflected + outgoing);
  // auto d         = sample_microfacet_pdf(roughness, normal, halfway,
  // outgoing);
  auto d = sample_microfacet_pdf(roughness, normal, halfway);
  return d / (4 * abs(dot(outgoing, halfway)));
}

// Pdf for refraction BRDF lobe sampling.
inline float sample_microfacet_refraction_pdf(float ior, float roughness,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    auto halfway = normalize(incoming + outgoing);
    return fresnel_dielectric(rel_ior, halfway, outgoing) *
           //  sample_microfacet_pdf(roughness, up_normal, halfway, outgoing) /
           sample_microfacet_pdf(roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  } else {
    auto halfway = -normalize(rel_ior * incoming + outgoing) *
                   (entering ? 1 : -1);
    // [Walter 2007] equation 17
    return (1 - fresnel_dielectric(rel_ior, halfway, outgoing)) *
           //  sample_microfacet_pdf(roughness, up_normal, halfway, outgoing) *
           sample_microfacet_pdf(roughness, up_normal, halfway) *
           abs(dot(halfway, outgoing)) /
           pow(rel_ior * dot(halfway, incoming) + dot(halfway, outgoing), 2);
  }
}

// Evaluate a delta specular BRDF lobe.
inline vec3f eval_delta_reflection(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return zero3f;
  return vec3f{1} * fresnel_dielectric(ior, normal, outgoing);
}

// Evaluate a delta metal BRDF lobe.
inline vec3f eval_delta_reflection(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return zero3f;
  return fresnel_conductor(eta, etak, normal, outgoing);
}

// Evaluate a delta transmission BRDF lobe.
inline vec3f eval_delta_transmission(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return zero3f;
  return vec3f{1};
}

// Evaluate a delta refraction BRDF lobe.
inline vec3f eval_delta_refraction(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (abs(ior - 1) < 1e-3)
    return dot(normal, incoming) * dot(normal, outgoing) <= 0 ? vec3f{1}
                                                              : vec3f{0};
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    return vec3f{1} * fresnel_dielectric(rel_ior, up_normal, outgoing);
  } else {
    return vec3f{1} * (1 / (rel_ior * rel_ior)) *
           (1 - fresnel_dielectric(rel_ior, up_normal, outgoing));
  }
}

// Sample a delta specular BRDF lobe.
inline vec3f sample_delta_reflection(
    float ior, const vec3f& normal, const vec3f& outgoing) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  return reflect(outgoing, normal);
}

// Sample a delta metal BRDF lobe.
inline vec3f sample_delta_reflection(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  return reflect(outgoing, normal);
}

// Sample a delta transmission BRDF lobe.
inline vec3f sample_delta_transmission(
    float ior, const vec3f& normal, const vec3f& outgoing) {
  if (dot(normal, outgoing) <= 0) return zero3f;
  return -outgoing;
}

// Sample a delta refraction BRDF lobe.
inline vec3f sample_delta_refraction(
    float ior, const vec3f& normal, const vec3f& outgoing, float rnl) {
  if (abs(ior - 1) < 1e-3) return -outgoing;
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (rnl < fresnel_dielectric(rel_ior, up_normal, outgoing)) {
    return reflect(outgoing, up_normal);
  } else {
    return refract(outgoing, up_normal, 1 / rel_ior);
  }
}

// Pdf for delta specular BRDF lobe sampling.
inline float sample_delta_reflection_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta metal BRDF lobe sampling.
inline float sample_delta_reflection_pdf(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta transmission BRDF lobe sampling.
inline float sample_delta_transmission_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta refraction BRDF lobe sampling.
inline float sample_delta_refraction_pdf(float ior, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (abs(ior - 1) < 1e-3)
    return dot(normal, incoming) * dot(normal, outgoing) < 0 ? 1 : 0;
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    return fresnel_dielectric(rel_ior, up_normal, outgoing);
  } else {
    return (1 - fresnel_dielectric(rel_ior, up_normal, outgoing));
  }
}

// Convert mean-free-path to transmission
inline vec3f mfp_to_transmission(const vec3f& mfp, float depth) {
  return exp(-depth / mfp);
}

// Evaluate transmittance
inline vec3f eval_transmittance(const vec3f& density, float distance) {
  return exp(-density * distance);
}

// Sample a distance proportionally to transmittance
inline float sample_transmittance(
    const vec3f& density, float max_distance, float rl, float rd) {
  auto channel  = clamp((int)(rl * 3), 0, 2);
  auto distance = (density[channel] == 0) ? flt_max
                                          : -log(1 - rd) / density[channel];
  return min(distance, max_distance);
}

// Pdf for distance sampling
inline float sample_transmittance_pdf(
    const vec3f& density, float distance, float max_distance) {
  if (distance < max_distance) {
    return sum(density * exp(-density * distance)) / 3;
  } else {
    return sum(exp(-density * max_distance)) / 3;
  }
}

// Evaluate phase function
inline float eval_phasefunction(
    float anisotropy, const vec3f& outgoing, const vec3f& incoming) {
  auto cosine = -dot(outgoing, incoming);
  auto denom  = 1 + anisotropy * anisotropy - 2 * anisotropy * cosine;
  return (1 - anisotropy * anisotropy) / (4 * pif * denom * sqrt(denom));
}

// Sample phase function
inline vec3f sample_phasefunction(
    float anisotropy, const vec3f& outgoing, const vec2f& rn) {
  auto cos_theta = 0.0f;
  if (abs(anisotropy) < 1e-3f) {
    cos_theta = 1 - 2 * rn.y;
  } else {
    auto square = (1 - anisotropy * anisotropy) /
                  (1 + anisotropy - 2 * anisotropy * rn.y);
    cos_theta = (1 + anisotropy * anisotropy - square * square) /
                (2 * anisotropy);
  }

  auto sin_theta      = sqrt(max(0.0f, 1 - cos_theta * cos_theta));
  auto phi            = 2 * pif * rn.x;
  auto local_incoming = vec3f{
      sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
  return basis_fromz(-outgoing) * local_incoming;
}

// Pdf for phase function sampling
inline float sample_phasefunction_pdf(
    float anisotropy, const vec3f& outgoing, const vec3f& incoming) {
  return eval_phasefunction(anisotropy, outgoing, incoming);
}

}  // namespace yocto

#endif
