//
// Implementation for Yocto/Trace.
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

#include "yocto_trace.h"

#include <algorithm>
#include <cstring>
#include <deque>
#include <memory>
#include <stdexcept>
#include <utility>

#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_json.h"
#include "yocto_parallel.h"
#include "yocto_sampling.h"
#include "yocto_shading.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::deque;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF RAY-SCENE INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// Build the bvh acceleration structure.
void init_bvh(trace_bvh& bvh, const scene_scene& scene,
    const trace_params& params, const progress_callback& progress_cb) {
  init_bvh(bvh, scene,
      bvh_params{(bvh_build_type)params.bvh, params.noparallel}, progress_cb);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// constant values
static const auto coat_ior       = 1.5f;
static const auto coat_roughness = 0.03f * 0.03f;

// Evaluate bsdf
trace_bsdf eval_bsdf(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv, const vec3f& normal, const vec3f& outgoing) {
  auto& material = get_material(scene, instance.material);
  auto  texcoord = eval_texcoord(scene, instance, element, uv);
  auto  color = material.color * xyz(eval_color(scene, instance, element, uv)) *
               xyz(eval_texture(scene, material.color_tex, texcoord, false));
  auto specular = material.specular *
                  eval_texture(scene, material.specular_tex, texcoord, true).x;
  auto metallic = material.metallic *
                  eval_texture(scene, material.metallic_tex, texcoord, true).x;
  auto roughness =
      material.roughness *
      eval_texture(scene, material.roughness_tex, texcoord, true).x;
  auto ior  = material.ior;
  auto coat = material.coat *
              eval_texture(scene, material.coat_tex, texcoord, true).x;
  auto transmission =
      material.transmission *
      eval_texture(scene, material.emission_tex, texcoord, true).x;
  auto translucency =
      material.translucency *
      eval_texture(scene, material.translucency_tex, texcoord, true).x;
  auto thin = material.thin || material.transmission == 0;

  // factors
  auto bsdf   = trace_bsdf{};
  auto weight = vec3f{1, 1, 1};
  bsdf.coat   = weight * coat;
  weight *= 1 - bsdf.coat * fresnel_dielectric(coat_ior, outgoing, normal);
  bsdf.metal = weight * metallic;
  weight *= 1 - metallic;
  bsdf.refraction = thin ? zero3f : (weight * transmission);
  weight *= 1 - (thin ? 0 : transmission);
  bsdf.specular = weight * specular;
  weight *= 1 - specular * fresnel_dielectric(ior, outgoing, normal);
  bsdf.transmission = thin ? (weight * transmission * color) : zero3f;
  weight *= 1 - (thin ? transmission : 0);
  bsdf.translucency = thin ? (weight * translucency * color)
                           : (weight * translucency);
  weight *= 1 - translucency;
  bsdf.diffuse   = weight * color;
  bsdf.meta      = reflectivity_to_eta(color);
  bsdf.metak     = zero3f;
  bsdf.roughness = roughness * roughness;
  bsdf.ior       = ior;

  // textures
  if (bsdf.diffuse != zero3f || bsdf.translucency != zero3f ||
      bsdf.roughness != 0) {
    bsdf.roughness = clamp(bsdf.roughness, coat_roughness, 1.0f);
  }
  if (bsdf.specular == zero3f && bsdf.metal == zero3f &&
      bsdf.transmission == zero3f && bsdf.refraction == zero3f) {
    bsdf.roughness = 1;
  }

  // weights
  bsdf.diffuse_pdf  = max(bsdf.diffuse);
  bsdf.specular_pdf = max(
      bsdf.specular * fresnel_dielectric(bsdf.ior, normal, outgoing));
  bsdf.metal_pdf = max(
      bsdf.metal * fresnel_conductor(bsdf.meta, bsdf.metak, normal, outgoing));
  bsdf.coat_pdf = max(
      bsdf.coat * fresnel_dielectric(coat_ior, normal, outgoing));
  bsdf.transmission_pdf = max(bsdf.transmission);
  bsdf.translucency_pdf = max(bsdf.translucency);
  bsdf.refraction_pdf   = max(bsdf.refraction);
  auto pdf_sum = bsdf.diffuse_pdf + bsdf.specular_pdf + bsdf.metal_pdf +
                 bsdf.coat_pdf + bsdf.transmission_pdf + bsdf.translucency_pdf +
                 bsdf.refraction_pdf;
  if (pdf_sum != 0) {
    bsdf.diffuse_pdf /= pdf_sum;
    bsdf.specular_pdf /= pdf_sum;
    bsdf.metal_pdf /= pdf_sum;
    bsdf.coat_pdf /= pdf_sum;
    bsdf.transmission_pdf /= pdf_sum;
    bsdf.translucency_pdf /= pdf_sum;
    bsdf.refraction_pdf /= pdf_sum;
  }
  return bsdf;
}

// check if a brdf is a delta
bool is_delta(const trace_bsdf& bsdf) { return bsdf.roughness == 0; }

// evaluate volume
trace_vsdf eval_vsdf(const scene_scene& scene, const scene_instance& instance,
    int element, const vec2f& uv) {
  auto& material = get_material(scene, instance.material);
  // initialize factors
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  auto color = material.color * xyz(eval_color(scene, instance, element, uv)) *
               xyz(eval_texture(scene, material.color_tex, texcoord, false));
  auto transmission =
      material.transmission *
      eval_texture(scene, material.emission_tex, texcoord, true).x;
  auto translucency =
      material.translucency *
      eval_texture(scene, material.translucency_tex, texcoord, true).x;
  auto thin = material.thin ||
              (material.transmission == 0 && material.translucency == 0);
  auto scattering =
      material.scattering *
      xyz(eval_texture(scene, material.scattering_tex, texcoord, false));
  auto scanisotropy = material.scanisotropy;
  auto trdepth      = material.trdepth;

  // factors
  auto vsdf       = trace_vsdf{};
  vsdf.density    = ((transmission != 0 || translucency != 0) && !thin)
                        ? -log(clamp(color, 0.0001f, 1.0f)) / trdepth
                        : zero3f;
  vsdf.scatter    = scattering;
  vsdf.anisotropy = scanisotropy;

  return vsdf;
}

// check if we have a volume
bool has_volume(const scene_scene& scene, const scene_instance& instance) {
  auto& material = get_material(scene, instance.material);
  return !material.thin && material.transmission != 0;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluate emission
static vec3f eval_emission(
    const vec3f& emission, const vec3f& normal, const vec3f& outgoing) {
  return emission;
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_bsdfcos(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (bsdf.roughness == 0) return zero3f;

  // accumulate the lobes
  auto brdfcos = zero3f;
  if (bsdf.diffuse != zero3f) {
    brdfcos += bsdf.diffuse *
               eval_diffuse_reflection(normal, outgoing, incoming);
  }
  if (bsdf.specular != zero3f) {
    brdfcos += bsdf.specular * eval_microfacet_reflection(bsdf.ior,
                                   bsdf.roughness, normal, outgoing, incoming);
  }
  if (bsdf.metal != zero3f) {
    brdfcos += bsdf.metal * eval_microfacet_reflection(bsdf.meta, bsdf.metak,
                                bsdf.roughness, normal, outgoing, incoming);
  }
  if (bsdf.coat != zero3f) {
    brdfcos += bsdf.coat * eval_microfacet_reflection(coat_ior, coat_roughness,
                               normal, outgoing, incoming);
  }
  if (bsdf.transmission != zero3f) {
    brdfcos += bsdf.transmission * eval_microfacet_transmission(bsdf.ior,
                                       bsdf.roughness, normal, outgoing,
                                       incoming);
  }
  if (bsdf.translucency != zero3f) {
    brdfcos += bsdf.translucency *
               eval_diffuse_transmission(normal, outgoing, incoming);
  }
  if (bsdf.refraction != zero3f) {
    brdfcos += bsdf.refraction * eval_microfacet_refraction(bsdf.ior,
                                     bsdf.roughness, normal, outgoing,
                                     incoming);
  }
  return brdfcos;
}

static vec3f eval_delta(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (bsdf.roughness != 0) return zero3f;

  auto brdfcos = zero3f;

  if (bsdf.specular != zero3f && bsdf.refraction == zero3f) {
    brdfcos += bsdf.specular *
               eval_delta_reflection(bsdf.ior, normal, outgoing, incoming);
  }
  if (bsdf.metal != zero3f) {
    brdfcos += bsdf.metal * eval_delta_reflection(bsdf.meta, bsdf.metak, normal,
                                outgoing, incoming);
  }
  if (bsdf.coat != zero3f) {
    brdfcos += bsdf.coat *
               eval_delta_reflection(coat_ior, normal, outgoing, incoming);
  }
  if (bsdf.transmission != zero3f) {
    brdfcos += bsdf.transmission *
               eval_delta_transmission(bsdf.ior, normal, outgoing, incoming);
  }
  if (bsdf.refraction != zero3f) {
    brdfcos += bsdf.refraction *
               eval_delta_refraction(bsdf.ior, normal, outgoing, incoming);
  }

  return brdfcos;
}

// Picks a direction based on the BRDF
static vec3f sample_bsdfcos(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  if (bsdf.roughness == 0) return zero3f;

  auto cdf = 0.0f;

  if (bsdf.diffuse_pdf != 0) {
    cdf += bsdf.diffuse_pdf;
    if (rnl < cdf) return sample_diffuse_reflection(normal, outgoing, rn);
  }

  if (bsdf.specular_pdf != 0 && bsdf.refraction_pdf == 0) {
    cdf += bsdf.specular_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          bsdf.ior, bsdf.roughness, normal, outgoing, rn);
  }

  if (bsdf.metal_pdf != 0) {
    cdf += bsdf.metal_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          bsdf.meta, bsdf.metak, bsdf.roughness, normal, outgoing, rn);
  }

  if (bsdf.coat_pdf != 0) {
    cdf += bsdf.coat_pdf;
    if (rnl < cdf)
      return sample_microfacet_reflection(
          coat_ior, coat_roughness, normal, outgoing, rn);
  }

  if (bsdf.transmission_pdf != 0) {
    cdf += bsdf.transmission_pdf;
    if (rnl < cdf)
      return sample_microfacet_transmission(
          bsdf.ior, bsdf.roughness, normal, outgoing, rn);
  }

  if (bsdf.translucency_pdf != 0) {
    cdf += bsdf.translucency_pdf;
    if (rnl < cdf) return sample_diffuse_transmission(normal, outgoing, rn);
  }

  if (bsdf.refraction_pdf != 0) {
    cdf += bsdf.refraction_pdf;
    if (rnl < cdf)
      return sample_microfacet_refraction(
          bsdf.ior, bsdf.roughness, normal, outgoing, rnl, rn);
  }

  return zero3f;
}

static vec3f sample_delta(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, float rnl) {
  if (bsdf.roughness != 0) return zero3f;

  // keep a weight sum to pick a lobe
  auto cdf = 0.0f;
  cdf += bsdf.diffuse_pdf;

  if (bsdf.specular_pdf != 0 && bsdf.refraction_pdf == 0) {
    cdf += bsdf.specular_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(bsdf.ior, normal, outgoing);
    }
  }

  if (bsdf.metal_pdf != 0) {
    cdf += bsdf.metal_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(bsdf.meta, bsdf.metak, normal, outgoing);
    }
  }

  if (bsdf.coat_pdf != 0) {
    cdf += bsdf.coat_pdf;
    if (rnl < cdf) {
      return sample_delta_reflection(coat_ior, normal, outgoing);
    }
  }

  if (bsdf.transmission_pdf != 0) {
    cdf += bsdf.transmission_pdf;
    if (rnl < cdf) {
      return sample_delta_transmission(bsdf.ior, normal, outgoing);
    }
  }

  if (bsdf.refraction_pdf != 0) {
    cdf += bsdf.refraction_pdf;
    if (rnl < cdf) {
      return sample_delta_refraction(bsdf.ior, normal, outgoing, rnl);
    }
  }

  return zero3f;
}

// Compute the weight for sampling the BRDF
static float sample_bsdfcos_pdf(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (bsdf.roughness == 0) return 0;

  auto pdf = 0.0f;

  if (bsdf.diffuse_pdf != 0) {
    pdf += bsdf.diffuse_pdf *
           sample_diffuse_reflection_pdf(normal, outgoing, incoming);
  }

  if (bsdf.specular_pdf != 0 && bsdf.refraction_pdf == 0) {
    pdf += bsdf.specular_pdf * sample_microfacet_reflection_pdf(bsdf.ior,
                                   bsdf.roughness, normal, outgoing, incoming);
  }

  if (bsdf.metal_pdf != 0) {
    pdf += bsdf.metal_pdf * sample_microfacet_reflection_pdf(bsdf.meta,
                                bsdf.metak, bsdf.roughness, normal, outgoing,
                                incoming);
  }

  if (bsdf.coat_pdf != 0) {
    pdf += bsdf.coat_pdf * sample_microfacet_reflection_pdf(coat_ior,
                               coat_roughness, normal, outgoing, incoming);
  }

  if (bsdf.transmission_pdf != 0) {
    pdf += bsdf.transmission_pdf * sample_microfacet_transmission_pdf(bsdf.ior,
                                       bsdf.roughness, normal, outgoing,
                                       incoming);
  }

  if (bsdf.translucency_pdf != 0) {
    pdf += bsdf.translucency_pdf *
           sample_diffuse_transmission_pdf(normal, outgoing, incoming);
  }

  if (bsdf.refraction_pdf != 0) {
    pdf += bsdf.refraction_pdf * sample_microfacet_refraction_pdf(bsdf.ior,
                                     bsdf.roughness, normal, outgoing,
                                     incoming);
  }

  return pdf;
}

static float sample_delta_pdf(const trace_bsdf& bsdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (bsdf.roughness != 0) return 0;

  auto pdf = 0.0f;
  if (bsdf.specular_pdf != 0 && bsdf.refraction_pdf == 0) {
    pdf += bsdf.specular_pdf *
           sample_delta_reflection_pdf(bsdf.ior, normal, outgoing, incoming);
  }
  if (bsdf.metal_pdf != 0) {
    pdf += bsdf.metal_pdf * sample_delta_reflection_pdf(bsdf.meta, bsdf.metak,
                                normal, outgoing, incoming);
  }
  if (bsdf.coat_pdf != 0) {
    pdf += bsdf.coat_pdf *
           sample_delta_reflection_pdf(coat_ior, normal, outgoing, incoming);
  }
  if (bsdf.transmission_pdf != 0) {
    pdf += bsdf.transmission_pdf *
           sample_delta_transmission_pdf(bsdf.ior, normal, outgoing, incoming);
  }
  if (bsdf.refraction_pdf != 0) {
    pdf += bsdf.refraction_pdf *
           sample_delta_refraction_pdf(bsdf.ior, normal, outgoing, incoming);
  }
  return pdf;
}

static vec3f eval_scattering(
    const trace_vsdf& vsdf, const vec3f& outgoing, const vec3f& incoming) {
  if (vsdf.density == zero3f) return zero3f;
  return vsdf.scatter * vsdf.density *
         eval_phasefunction(vsdf.anisotropy, outgoing, incoming);
}

static vec3f sample_scattering(
    const trace_vsdf& vsdf, const vec3f& outgoing, float rnl, const vec2f& rn) {
  if (vsdf.density == zero3f) return zero3f;
  return sample_phasefunction(vsdf.anisotropy, outgoing, rn);
}

static float sample_scattering_pdf(
    const trace_vsdf& vsdf, const vec3f& outgoing, const vec3f& incoming) {
  if (vsdf.density == zero3f) return 0;
  return sample_phasefunction_pdf(vsdf.anisotropy, outgoing, incoming);
}

// Sample camera
static ray3f sample_camera(const scene_camera& camera, const vec2i& ij,
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
                puv.y < 0.5f ? sqrt(2 * puv.y) - 1 : 1 - sqrt(2 - 2 * puv.y),
            } +
        offset;
    auto uv = vec2f{
        (ij.x + fuv.x) / image_size.x, (ij.y + fuv.y) / image_size.y};
    return eval_camera(camera, uv, sample_disk(luv));
  }
}

// Sample lights wrt solid angle
static vec3f sample_lights(const scene_scene& scene, const trace_lights& lights,
    const vec3f& position, float rl, float rel, const vec2f& ruv) {
  auto  light_id = sample_uniform((int)lights.lights.size(), rl);
  auto& light    = lights.lights[light_id];
  if (light.instance != invalid_handle) {
    auto& instance  = get_instance(scene, light.instance);
    auto& shape     = get_shape(scene, instance.shape);
    auto  element   = sample_discrete_cdf(light.elements_cdf, rel);
    auto  uv        = (!shape.triangles.empty()) ? sample_triangle(ruv) : ruv;
    auto  lposition = eval_position(scene, instance, element, uv);
    return normalize(lposition - position);
  } else if (light.environment != invalid_handle) {
    auto& environment = get_environment(scene, light.environment);
    if (environment.emission_tex != invalid_handle) {
      auto& emission_tex = get_texture(scene, environment.emission_tex);
      auto  idx          = sample_discrete_cdf(light.elements_cdf, rel);
      auto  size         = texture_size(emission_tex);
      auto  uv           = vec2f{
          ((idx % size.x) + 0.5f) / size.x, ((idx / size.x) + 0.5f) / size.y};
      return transform_direction(environment.frame,
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
static float sample_lights_pdf(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction) {
  auto pdf = 0.0f;
  for (auto& light : lights.lights) {
    if (light.instance != invalid_handle) {
      auto& instance = get_instance(scene, light.instance);
      // check all intersection
      auto lpdf          = 0.0f;
      auto next_position = position;
      for (auto bounce = 0; bounce < 100; bounce++) {
        auto intersection = intersect_bvh(
            bvh, scene, light.instance, {next_position, direction});
        if (!intersection.hit) break;
        // accumulate pdf
        auto lposition = eval_position(
            scene, instance, intersection.element, intersection.uv);
        auto lnormal = eval_element_normal(
            scene, instance, intersection.element);
        // prob triangle * area triangle = area triangle mesh
        auto area = light.elements_cdf.back();
        lpdf += distance_squared(lposition, position) /
                (abs(dot(lnormal, direction)) * area);
        // continue
        next_position = lposition + direction * 1e-3f;
      }
      pdf += lpdf;
    } else if (light.environment != invalid_handle) {
      auto& environment = get_environment(scene, light.environment);
      if (environment.emission_tex != invalid_handle) {
        auto& emission_tex = get_texture(scene, environment.emission_tex);
        auto  size         = texture_size(emission_tex);
        auto  wl = transform_direction(inverse(environment.frame), direction);
        auto  texcoord = vec2f{atan2(wl.z, wl.x) / (2 * pif),
            acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
        if (texcoord.x < 0) texcoord.x += 1;
        auto i    = clamp((int)(texcoord.x * size.x), 0, size.x - 1);
        auto j    = clamp((int)(texcoord.y * size.y), 0, size.y - 1);
        auto prob = sample_discrete_cdf_pdf(
                        light.elements_cdf, j * size.x + i) /
                    light.elements_cdf.back();
        auto angle = (2 * pif / size.x) * (pif / size.y) *
                     sin(pif * (j + 0.5f) / size.y);
        pdf += prob / angle;
      } else {
        pdf += 1 / (4 * pif);
      }
    }
  }
  pdf *= sample_uniform_pdf((int)lights.lights.size());
  return pdf;
}

// Recursive path tracing.
static vec4f trace_path(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance      = zero3f;
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = vector<trace_vsdf>{};
  auto max_roughness = 0.0f;
  auto hit           = !params.envhidden && !scene.environments.empty();

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // handle transmission if inside a volume
    auto in_volume = false;
    if (!volume_stack.empty()) {
      auto& vsdf     = volume_stack.back();
      auto  distance = sample_transmittance(
          vsdf.density, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(vsdf.density, distance) /
                sample_transmittance_pdf(
                    vsdf.density, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto& instance = scene.instances[intersection.instance];
      auto  element  = intersection.element;
      auto  uv       = intersection.uv;
      auto  position = eval_position(scene, instance, element, uv);
      auto normal = eval_shading_normal(scene, instance, element, uv, outgoing);
      auto emission = eval_emission(
          scene, instance, element, uv, normal, outgoing);
      auto opacity = eval_opacity(
          scene, instance, element, uv, normal, outgoing);
      auto bsdf = eval_bsdf(scene, instance, element, uv, normal, outgoing);

      // correct roughness
      if (params.nocaustics) {
        max_roughness  = max(bsdf.roughness, max_roughness);
        bsdf.roughness = max_roughness;
      }

      // handle opacity
      if (opacity < 1 && rand1f(rng) >= opacity) {
        ray = {position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }
      hit = true;

      // accumulate emission
      radiance += weight * eval_emission(emission, normal, outgoing);

      // next direction
      auto incoming = zero3f;
      if (!is_delta(bsdf)) {
        if (rand1f(rng) < 0.5f) {
          incoming = sample_bsdfcos(
              bsdf, normal, outgoing, rand1f(rng), rand2f(rng));
        } else {
          incoming = sample_lights(
              scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
        }
        weight *= eval_bsdfcos(bsdf, normal, outgoing, incoming) /
                  (0.5f * sample_bsdfcos_pdf(bsdf, normal, outgoing, incoming) +
                      0.5f * sample_lights_pdf(
                                 scene, bvh, lights, position, incoming));
      } else {
        incoming = sample_delta(bsdf, normal, outgoing, rand1f(rng));
        weight *= eval_delta(bsdf, normal, outgoing, incoming) /
                  sample_delta_pdf(bsdf, normal, outgoing, incoming);
      }

      // update volume stack
      if (has_volume(scene, instance) &&
          dot(normal, outgoing) * dot(normal, incoming) < 0) {
        if (volume_stack.empty()) {
          auto vsdf = eval_vsdf(scene, instance, element, uv);
          volume_stack.push_back(vsdf);
        } else {
          volume_stack.pop_back();
        }
      }

      // setup next iteration
      ray = {position, incoming};
    } else {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto  position = ray.o + ray.d * intersection.distance;
      auto& vsdf     = volume_stack.back();

      // handle opacity
      hit = true;

      // accumulate emission
      // radiance += weight * eval_volemission(emission, outgoing);

      // next direction
      auto incoming = zero3f;
      if (rand1f(rng) < 0.5f) {
        incoming = sample_scattering(vsdf, outgoing, rand1f(rng), rand2f(rng));
      } else {
        incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      weight *=
          eval_scattering(vsdf, outgoing, incoming) /
          (0.5f * sample_scattering_pdf(vsdf, outgoing, incoming) +
              0.5f * sample_lights_pdf(scene, bvh, lights, position, incoming));

      // setup next iteration
      ray = {position, incoming};
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

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Recursive path tracing.
static vec4f trace_naive(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = !params.envhidden && !scene.environments.empty();

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto instance = scene.instances[intersection.instance];
    auto element  = intersection.element;
    auto uv       = intersection.uv;
    auto position = eval_position(scene, instance, element, uv);
    auto normal   = eval_shading_normal(scene, instance, element, uv, outgoing);
    auto emission = eval_emission(
        scene, instance, element, uv, normal, outgoing);
    auto opacity = eval_opacity(scene, instance, element, uv, normal, outgoing);
    auto bsdf    = eval_bsdf(scene, instance, element, uv, normal, outgoing);

    // handle opacity
    if (opacity < 1 && rand1f(rng) >= opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    radiance += weight * eval_emission(emission, normal, outgoing);

    // next direction
    auto incoming = zero3f;
    if (bsdf.roughness != 0) {
      incoming = sample_bsdfcos(
          bsdf, normal, outgoing, rand1f(rng), rand2f(rng));
      weight *= eval_bsdfcos(bsdf, normal, outgoing, incoming) /
                sample_bsdfcos_pdf(bsdf, normal, outgoing, incoming);
    } else {
      incoming = sample_delta(bsdf, normal, outgoing, rand1f(rng));
      weight *= eval_delta(bsdf, normal, outgoing, incoming) /
                sample_delta_pdf(bsdf, normal, outgoing, incoming);
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
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// Eyelight for quick previewing.
static vec4f trace_eyelight(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance = zero3f;
  auto weight   = vec3f{1, 1, 1};
  auto ray      = ray_;
  auto hit      = !params.envhidden && !scene.environments.empty();

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_bvh(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto instance = scene.instances[intersection.instance];
    auto element  = intersection.element;
    auto uv       = intersection.uv;
    auto position = eval_position(scene, instance, element, uv);
    auto normal   = eval_shading_normal(scene, instance, element, uv, outgoing);
    auto emission = eval_emission(
        scene, instance, element, uv, normal, outgoing);
    auto opacity = eval_opacity(scene, instance, element, uv, normal, outgoing);
    auto bsdf    = eval_bsdf(scene, instance, element, uv, normal, outgoing);

    // handle opacity
    if (opacity < 1 && rand1f(rng) >= opacity) {
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }
    hit = true;

    // accumulate emission
    auto incoming = outgoing;
    radiance += weight * eval_emission(emission, normal, outgoing);

    // brdf * light
    radiance += weight * pif * eval_bsdfcos(bsdf, normal, outgoing, incoming);

    // continue path
    if (bsdf.roughness != 0) break;
    incoming = sample_delta(bsdf, normal, outgoing, rand1f(rng));
    weight *= eval_delta(bsdf, normal, outgoing, incoming) /
              sample_delta_pdf(bsdf, normal, outgoing, incoming);
    if (weight == zero3f || !isfinite(weight)) break;

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance.x, radiance.y, radiance.z, hit ? 1.0f : 0.0f};
}

// False color rendering
static vec4f trace_falsecolor(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params) {
  // intersect next point
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) {
    return {0, 0, 0, 0};
  }

  // prepare shading point
  auto outgoing = -ray.d;
  auto instance = scene.instances[intersection.instance];
  auto element  = intersection.element;
  auto uv       = intersection.uv;
  auto position = eval_position(scene, instance, element, uv);
  auto normal   = eval_shading_normal(scene, instance, element, uv, outgoing);
  auto gnormal  = eval_element_normal(scene, instance, element);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  auto color    = eval_color(scene, instance, element, uv);
  auto emission = eval_emission(scene, instance, element, uv, normal, outgoing);
  auto opacity  = eval_opacity(scene, instance, element, uv, normal, outgoing);
  auto bsdf     = eval_bsdf(scene, instance, element, uv, normal, outgoing);

  // hash color
  auto hashed_color = [](int id) {
    auto hashed = std::hash<int>()(id);
    auto rng    = make_rng(trace_default_seed, hashed);
    return pow(0.5f + 0.5f * rand3f(rng), 2.2f);
  };

  // make vec4f
  auto make_vec = [](const vec3f& xyz, float w) {
    return vec4f{xyz.x, xyz.y, xyz.z, w};
  };

  switch (params.falsecolor) {
    case trace_falsecolor_type::position:
      return make_vec(position * 0.5f + 0.5f, 1);
    case trace_falsecolor_type::normal:
      return make_vec(normal * 0.5f + 0.5f, 1);
    case trace_falsecolor_type::frontfacing:
      return dot(normal, -ray.d) > 0 ? vec4f{0, 1, 0, 1} : vec4f{1, 0, 0, 1};
    case trace_falsecolor_type::gnormal:
      return make_vec(gnormal * 0.5f + 0.5f, 1);
    case trace_falsecolor_type::gfrontfacing:
      return make_vec(
          dot(gnormal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0}, 1);
    case trace_falsecolor_type::texcoord:
      return {fmod(texcoord.x, 1.0f), fmod(texcoord.y, 1.0f), 0, 1};
    case trace_falsecolor_type::color: return make_vec(xyz(color), 1);
    case trace_falsecolor_type::emission: return make_vec(emission, 1);
    case trace_falsecolor_type::diffuse: return make_vec(bsdf.diffuse, 1);
    case trace_falsecolor_type::specular: return make_vec(bsdf.specular, 1);
    case trace_falsecolor_type::coat: return make_vec(bsdf.coat, 1);
    case trace_falsecolor_type::metal: return make_vec(bsdf.metal, 1);
    case trace_falsecolor_type::transmission:
      return make_vec(bsdf.transmission, 1);
    case trace_falsecolor_type::translucency:
      return make_vec(bsdf.translucency, 1);
    case trace_falsecolor_type::refraction: return make_vec(bsdf.refraction, 1);
    case trace_falsecolor_type::roughness:
      return {bsdf.roughness, bsdf.roughness, bsdf.roughness, 1};
    case trace_falsecolor_type::opacity: return {opacity, opacity, opacity, 1};
    case trace_falsecolor_type::ior: return {bsdf.ior, bsdf.ior, bsdf.ior, 1};
    case trace_falsecolor_type::element:
      return make_vec(hashed_color(intersection.element), 1);
    case trace_falsecolor_type::instance:
      return make_vec(hashed_color(intersection.instance), 1);
    case trace_falsecolor_type::highlight: {
      if (emission == zero3f) emission = {0.2f, 0.2f, 0.2f};
      return make_vec(emission * abs(dot(-ray.d, normal)), 1);
    } break;
    default: return {0, 0, 0, 0};
  }
}

static vec4f trace_albedo(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params, int bounce) {
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) {
    auto radiance = eval_environment(scene, ray.d);
    return {radiance.x, radiance.y, radiance.z, 1};
  }

  // prepare shading point
  auto  outgoing = -ray.d;
  auto& instance = scene.instances[intersection.instance];
  auto  element  = intersection.element;
  auto  uv       = intersection.uv;
  auto& material = get_material(scene, instance.material);
  auto  position = eval_position(scene, instance, element, uv);
  auto  normal   = eval_shading_normal(scene, instance, element, uv, outgoing);
  auto  texcoord = eval_texcoord(scene, instance, element, uv);
  auto  color    = eval_color(scene, instance, element, uv);
  auto emission = eval_emission(scene, instance, element, uv, normal, outgoing);
  auto opacity  = eval_opacity(scene, instance, element, uv, normal, outgoing);
  auto bsdf     = eval_bsdf(scene, instance, element, uv, normal, outgoing);

  if (emission != zero3f) {
    return {emission.x, emission.y, emission.z, 1};
  }

  auto albedo = material.color * xyz(color) *
                xyz(eval_texture(scene, material.color_tex, texcoord, false));

  // handle opacity
  if (opacity < 1.0f) {
    auto blend_albedo = trace_albedo(scene, bvh, lights,
        ray3f{position + ray.d * 1e-2f, ray.d}, rng, params, bounce);
    return lerp(blend_albedo, vec4f{albedo.x, albedo.y, albedo.z, 1}, opacity);
  }

  if (bsdf.roughness < 0.05 && bounce < 5) {
    if (bsdf.transmission != zero3f && material.thin) {
      auto incoming     = -outgoing;
      auto trans_albedo = trace_albedo(scene, bvh, lights,
          ray3f{position, incoming}, rng, params, bounce + 1);

      incoming         = reflect(outgoing, normal);
      auto spec_albedo = trace_albedo(scene, bvh, lights,
          ray3f{position, incoming}, rng, params, bounce + 1);

      auto fresnel = fresnel_dielectric(material.ior, outgoing, normal);
      auto dielectric_albedo = lerp(trans_albedo, spec_albedo, fresnel);
      return dielectric_albedo * vec4f{albedo.x, albedo.y, albedo.z, 1};
    } else if (bsdf.metal != zero3f) {
      auto incoming    = reflect(outgoing, normal);
      auto refl_albedo = trace_albedo(scene, bvh, lights,
          ray3f{position, incoming}, rng, params, bounce + 1);
      return refl_albedo * vec4f{albedo.x, albedo.y, albedo.z, 1};
    }
  }

  return {albedo.x, albedo.y, albedo.z, 1};
}

static vec4f trace_albedo(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params) {
  auto albedo = trace_albedo(scene, bvh, lights, ray, rng, params, 0);
  return clamp(albedo, 0.0, 1.0);
}

static vec4f trace_normal(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params, int bounce) {
  auto intersection = intersect_bvh(bvh, scene, ray);
  if (!intersection.hit) {
    return {0, 0, 0, 1};
  }

  // prepare shading point
  auto  outgoing = -ray.d;
  auto& instance = scene.instances[intersection.instance];
  auto  element  = intersection.element;
  auto  uv       = intersection.uv;
  auto& material = get_material(scene, instance.material);
  auto  position = eval_position(scene, instance, element, uv);
  auto  normal   = eval_shading_normal(scene, instance, element, uv, outgoing);
  auto  opacity  = eval_opacity(scene, instance, element, uv, normal, outgoing);
  auto  bsdf     = eval_bsdf(scene, instance, element, uv, normal, outgoing);

  // handle opacity
  if (opacity < 1.0f) {
    auto normal = trace_normal(scene, bvh, lights,
        ray3f{position + ray.d * 1e-2f, ray.d}, rng, params, bounce);
    return lerp(normal, normal, opacity);
  }

  if (bsdf.roughness < 0.05f && bounce < 5) {
    if (bsdf.transmission != zero3f && material.thin) {
      auto incoming   = -outgoing;
      auto trans_norm = trace_normal(scene, bvh, lights,
          ray3f{position, incoming}, rng, params, bounce + 1);

      incoming       = reflect(outgoing, normal);
      auto spec_norm = trace_normal(scene, bvh, lights,
          ray3f{position, incoming}, rng, params, bounce + 1);

      auto fresnel = fresnel_dielectric(material.ior, outgoing, normal);
      return lerp(trans_norm, spec_norm, fresnel);
    } else if (bsdf.metal != zero3f) {
      auto incoming = reflect(outgoing, normal);
      return trace_normal(scene, bvh, lights, ray3f{position, incoming}, rng,
          params, bounce + 1);
    }
  }

  return {normal.x, normal.y, normal.z, 1};
}

static vec4f trace_normal(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params) {
  return trace_normal(scene, bvh, lights, ray, rng, params, 0);
}

// Trace a single ray from the camera using the given algorithm.
using sampler_func = vec4f (*)(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params);
static sampler_func get_trace_sampler_func(const trace_params& params) {
  switch (params.sampler) {
    case trace_sampler_type::path: return trace_path;
    case trace_sampler_type::naive: return trace_naive;
    case trace_sampler_type::eyelight: return trace_eyelight;
    case trace_sampler_type::falsecolor: return trace_falsecolor;
    case trace_sampler_type::albedo: return trace_albedo;
    case trace_sampler_type::normal: return trace_normal;
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
    case trace_sampler_type::albedo: return false;
    case trace_sampler_type::normal: return false;
    default: {
      throw std::runtime_error("sampler unknown");
      return false;
    }
  }
}

// Trace a block of samples
void trace_sample(trace_state& state, const scene_scene& scene,
    const scene_camera& camera, const trace_bvh& bvh,
    const trace_lights& lights, const vec2i& ij, const trace_params& params) {
  auto sampler = get_trace_sampler_func(params);
  auto ray     = sample_camera(camera, ij, state.render.imsize(),
      rand2f(state.rngs[ij]), rand2f(state.rngs[ij]), params.tentfilter);
  auto sample  = sampler(scene, bvh, lights, ray, state.rngs[ij], params);
  if (!isfinite(xyz(sample))) sample = {0, 0, 0, sample.w};
  if (max(sample) > params.clamp)
    sample = sample * (params.clamp / max(sample));
  state.accumulation[ij] += sample;
  state.samples[ij] += 1;
  auto radiance    = state.accumulation[ij].w != 0
                         ? xyz(state.accumulation[ij]) / state.accumulation[ij].w
                         : zero3f;
  auto coverage    = state.accumulation[ij].w / state.samples[ij];
  state.render[ij] = {radiance.x, radiance.y, radiance.z, coverage};
}

// Init a sequence of random number generators.
void init_state(trace_state& state, const scene_scene& scene,
    const scene_camera& camera, const trace_params& params) {
  auto image_size = (camera.aspect >= 1)
                        ? vec2i{params.resolution,
                              (int)round(params.resolution / camera.aspect)}
                        : vec2i{(int)round(params.resolution * camera.aspect),
                              params.resolution};
  state.render.assign(image_size, zero4f);
  state.accumulation.assign(image_size, zero4f);
  state.samples.assign(image_size, 0);
  state.rngs.assign(image_size, {});
  auto rng_ = make_rng(1301081);
  for (auto& rng : state.rngs) {
    rng = make_rng(params.seed, rand1i(rng_, 1 << 31) / 2 + 1);
  }
}

// Forward declaration
static trace_light& add_light(trace_lights& lights) {
  return lights.lights.emplace_back();
}

// Init trace lights
void init_lights(trace_lights& lights, const scene_scene& scene,
    const trace_params& params, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);

  lights.lights.clear();

  for (auto handle = 0; handle < scene.instances.size(); handle++) {
    auto& instance = get_instance(scene, handle);
    auto& material = get_material(scene, instance.material);
    if (material.emission == zero3f) continue;
    auto& shape = get_shape(scene, instance.shape);
    if (shape.triangles.empty() && shape.quads.empty()) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    auto& light       = add_light(lights);
    light.instance    = handle;
    light.environment = invalid_handle;
    if (!shape.triangles.empty()) {
      light.elements_cdf = vector<float>(shape.triangles.size());
      for (auto idx = 0; idx < light.elements_cdf.size(); idx++) {
        auto& t                 = shape.triangles[idx];
        light.elements_cdf[idx] = triangle_area(
            shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
        if (idx != 0) light.elements_cdf[idx] += light.elements_cdf[idx - 1];
      }
    }
    if (!shape.quads.empty()) {
      light.elements_cdf = vector<float>(shape.quads.size());
      for (auto idx = 0; idx < light.elements_cdf.size(); idx++) {
        auto& t                 = shape.quads[idx];
        light.elements_cdf[idx] = quad_area(shape.positions[t.x],
            shape.positions[t.y], shape.positions[t.z], shape.positions[t.w]);
        if (idx != 0) light.elements_cdf[idx] += light.elements_cdf[idx - 1];
      }
    }
  }
  for (auto handle = 0; handle < scene.environments.size(); handle++) {
    auto& environment = get_environment(scene, handle);
    if (environment.emission == zero3f) continue;
    if (progress_cb) progress_cb("build light", progress.x++, ++progress.y);
    auto& light       = add_light(lights);
    light.instance    = invalid_handle;
    light.environment = handle;
    if (environment.emission_tex != invalid_handle) {
      auto texture       = get_texture(scene, environment.emission_tex);
      auto size          = texture_size(texture);
      light.elements_cdf = vector<float>(size.x * size.y);
      if (size != zero2i) {
        for (auto i = 0; i < light.elements_cdf.size(); i++) {
          auto ij               = vec2i{i % size.x, i / size.x};
          auto th               = (ij.y + 0.5f) * pif / size.y;
          auto value            = lookup_texture(texture, ij);
          light.elements_cdf[i] = max(value) * sin(th);
          if (i != 0) light.elements_cdf[i] += light.elements_cdf[i - 1];
        }
      }
    }
  }

  // handle progress
  if (progress_cb) progress_cb("build light", progress.x++, progress.y);
}

// Progressively computes an image.
image<vec4f> trace_image(const scene_scene& scene, const scene_camera& camera,
    const trace_params& params, const progress_callback& progress_cb,
    const image_callback& image_cb) {
  auto bvh = trace_bvh{};
  init_bvh(bvh, scene, params, progress_cb);

  auto lights = trace_lights{};
  init_lights(lights, scene, params, progress_cb);

  return trace_image(scene, camera, bvh, lights, params, progress_cb, image_cb);
}

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image(const scene_scene& scene, const scene_camera& camera,
    const trace_bvh& bvh, const trace_lights& lights,
    const trace_params& params, const progress_callback& progress_cb,
    const image_callback& image_cb) {
  auto state = trace_state{};
  init_state(state, scene, camera, params);

  for (auto sample = 0; sample < params.samples; sample++) {
    if (progress_cb) progress_cb("trace image", sample, params.samples);
    if (params.noparallel) {
      for (auto j = 0; j < state.render.height(); j++) {
        for (auto i = 0; i < state.render.width(); i++) {
          trace_sample(state, scene, camera, bvh, lights, {i, j}, params);
        }
      }
    } else {
      parallel_for(
          state.render.width(), state.render.height(), [&](int i, int j) {
            trace_sample(state, scene, camera, bvh, lights, {i, j}, params);
          });
    }
    if (image_cb) image_cb(state.render, sample + 1, params.samples);
  }

  if (progress_cb) progress_cb("trace image", params.samples, params.samples);
  return state.render;
}

// [experimental] Asynchronous interface
void trace_start(trace_state& state, const scene_scene& scene,
    const scene_camera& camera, const trace_bvh& bvh,
    const trace_lights& lights, const trace_params& params,
    const progress_callback& progress_cb, const image_callback& image_cb,
    const async_callback& async_cb) {
  init_state(state, scene, camera, params);
  state.worker = {};
  state.stop   = false;

  // render preview
  if (progress_cb) progress_cb("trace preview", 0, params.samples);
  auto pprms = params;
  pprms.resolution /= params.pratio;
  pprms.samples = 1;
  auto preview  = trace_image(scene, camera, bvh, lights, pprms);
  for (auto j = 0; j < state.render.height(); j++) {
    for (auto i = 0; i < state.render.width(); i++) {
      auto pi              = clamp(i / params.pratio, 0, preview.width() - 1),
           pj              = clamp(j / params.pratio, 0, preview.height() - 1);
      state.render[{i, j}] = preview[{pi, pj}];
    }
  }
  if (image_cb) image_cb(state.render, 0, params.samples);

  // start renderer
  state.worker = std::async(
      std::launch::async, [=, &state, &scene, &lights, &camera, &bvh]() {
        for (auto sample = 0; sample < params.samples; sample++) {
          if (state.stop) return;
          if (progress_cb) progress_cb("trace image", sample, params.samples);
          parallel_for(
              state.render.width(), state.render.height(), [&](int i, int j) {
                if (state.stop) return;
                trace_sample(state, scene, camera, bvh, lights, {i, j}, params);
                if (async_cb)
                  async_cb(state.render, sample, params.samples, {i, j});
              });
          if (image_cb) image_cb(state.render, sample + 1, params.samples);
        }
        if (progress_cb)
          progress_cb("trace image", params.samples, params.samples);
        if (image_cb) image_cb(state.render, params.samples, params.samples);
      });
}
void trace_stop(trace_state& state) {
  state.stop = true;
  if (state.worker.valid()) state.worker.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE IO
// -----------------------------------------------------------------------------
namespace yocto {

// clang-format off
    
 void serialize_value(json_mode mode,
    json_value& json, trace_params& value, const string& description) {
  serialize_object(mode, json, value, description);
  serialize_property(mode, json, value.resolution, "resolution", "Image resolution.");
  serialize_property(mode, json, value.sampler, "sampler", "Sampler type.");
  serialize_property(mode, json, value.falsecolor, "falsecolor", "False color type.");
  serialize_property(mode, json, value.samples, "samples", "Number of samples.");
  serialize_property(mode, json, value.bounces, "bounces", "Number of bounces.");
  serialize_property(mode, json, value.clamp, "clamp", "Clamp value.");
  serialize_property(mode, json, value.nocaustics, "nocaustics", "Disable caustics.");
  serialize_property(mode, json, value.envhidden, "envhidden", "Hide environment.");
  serialize_property(mode, json, value.tentfilter, "tentfilter", "Filter image.");
  serialize_property(mode, json, value.seed, "seed", "Random seed.");
  serialize_property(mode, json, value.bvh, "bvh", "Bvh type.");
  serialize_property(mode, json, value.noparallel, "noparallel", "Disable threading.");
  serialize_property(mode, json, value.pratio, "pratio", "Preview ratio.");
  serialize_property(mode, json, value.exposure, "exposure", "Image exposure.");
}

// Json enum conventions
 const vector<pair<trace_bvh_type, string>>& json_enum_labels(
    trace_bvh_type) {
  static const auto trace_bvh_labels = vector<pair<trace_bvh_type, string>>{
      {trace_bvh_type::default_, "default"},
      {trace_bvh_type::highquality, "highquality"},
      {trace_bvh_type::middle, "middle"},
      {trace_bvh_type::balanced, "balanced"},
#ifdef YOCTO_EMBREE
      {trace_bvh_type::embree_default, "embree-default"},
      {trace_bvh_type::embree_highquality, "embree-highquality"},
      {trace_bvh_type::embree_compact, "embree-compact"},
#endif
  };
  return trace_bvh_labels;
}

 const vector<pair<trace_falsecolor_type, string>>& json_enum_labels(
    trace_falsecolor_type) {
  static const auto trace_falsecolor_labels =
      vector<pair<trace_falsecolor_type, string>>{
          {trace_falsecolor_type::position, "position"},
          {trace_falsecolor_type::normal, "normal"},
          {trace_falsecolor_type::frontfacing, "frontfacing"},
          {trace_falsecolor_type::gnormal, "gnormal"},
          {trace_falsecolor_type::gfrontfacing, "gfrontfacing"},
          {trace_falsecolor_type::texcoord, "texcoord"},
          {trace_falsecolor_type::color, "color"},
          {trace_falsecolor_type::emission, "emission"},
          {trace_falsecolor_type::diffuse, "diffuse"},
          {trace_falsecolor_type::specular, "specular"},
          {trace_falsecolor_type::coat, "coat"},
          {trace_falsecolor_type::metal, "metal"},
          {trace_falsecolor_type::transmission, "transmission"},
          {trace_falsecolor_type::translucency, "translucency"},
          {trace_falsecolor_type::refraction, "refraction"},
          {trace_falsecolor_type::roughness, "roughness"},
          {trace_falsecolor_type::opacity, "opacity"},
          {trace_falsecolor_type::ior, "ior"},
          {trace_falsecolor_type::instance, "instance"},
          {trace_falsecolor_type::element, "element"},
          {trace_falsecolor_type::highlight, "highlight"}};
  return trace_falsecolor_labels;
}

const vector<pair<trace_sampler_type, string>>& json_enum_labels(
    trace_sampler_type) {
  static const auto trace_sampler_labels =
      vector<pair<trace_sampler_type, string>>{
          {trace_sampler_type::path, "path"},
          {trace_sampler_type::naive, "naive"},
          {trace_sampler_type::eyelight, "eyelight"},
          {trace_sampler_type::falsecolor, "falsecolor"},
          {trace_sampler_type::albedo, "albedo"},
          {trace_sampler_type::normal, "normal"}};
  return trace_sampler_labels;
}

// clang-format on

}  // namespace yocto
