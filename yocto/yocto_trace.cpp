//
// Implementation for Yocto/Trace.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Trace stats.
atomic<uint64_t> _trace_npaths{0};
atomic<uint64_t> _trace_nrays{0};

// Trace point
struct trace_point {
    int             instance_id  = -1;
    int             element_id   = -1;
    vec2f           element_uv   = zero_vec2f;
    vec3f           position     = zero_vec3f;
    vec3f           normal       = zero_vec3f;
    vec2f           texturecoord = zero_vec2f;
    vec3f           emission     = zero_vec3f;
    microfacet_brdf brdf         = {};
    float           opacity      = 1;
};

// Make a trace point
trace_point make_trace_point(const yocto_scene& scene, int instance_id,
    int element_id, const vec2f& element_uv,
    const vec3f& shading_direction = zero_vec3f) {
    auto& instance    = scene.instances[instance_id];
    auto  point       = trace_point();
    point.instance_id = instance_id;
    point.element_id  = element_id;
    point.element_uv  = element_uv;
    point.position    = evaluate_instance_position(
        scene, instance, element_id, element_uv);
    point.normal = evaluate_instance_shading_normal(
        scene, instance, element_id, element_uv, -shading_direction);
    point.texturecoord = evaluate_instance_texturecoord(
        scene, instance, element_id, element_uv);
    point.emission = evaluate_instance_emission(
        scene, instance, element_id, element_uv);
    point.brdf = evaluate_instance_brdf(scene, instance, element_id, element_uv);
    point.opacity = evaluate_instance_opacity(
        scene, instance, element_id, element_uv);
    return point;
}

// Intersects a ray and returns a point
trace_point trace_ray(const yocto_scene& scene, const bvh_scene& bvh,
    const vec3f& position, const vec3f& direction) {
    auto isec = intersect_scene(scene, bvh, make_ray(position, direction));
    _trace_nrays += 1;
    if (isec.instance_id >= 0) {
        return make_trace_point(scene, isec.instance_id, isec.element_id,
            isec.element_uv, direction);
    } else {
        auto point     = trace_point();
        point.emission = evaluate_environment_emission(scene, direction);
        return point;
    }
}

// Intersects a ray and returns a point accounting for opacity treated as
// coveregae
trace_point trace_ray_with_opacity(const yocto_scene& scene,
    const bvh_scene& bvh, const vec3f& position_, const vec3f& direction,
    rng_state& rng, int max_bounces) {
    auto position = position_;
    for (auto b = 0; b < max_bounces * 10; b++) {
        auto point = trace_ray(scene, bvh, position, direction);
        if (point.instance_id < 0) return point;
        if (point.opacity > 0.999f) return point;
        if (get_random_float(rng) < point.opacity) return point;
        position = point.position + direction * (float)default_ray_eps;
    }
    return {};
}

// Intersect a scene handling opacity.
scene_intersection intersect_scene_with_opacity(const yocto_scene& scene,
    const bvh_scene& bvh, const ray3f& ray_, rng_state& rng, int max_bounces) {
    auto ray = ray_;
    for (auto b = 0; b < max_bounces; b++) {
        _trace_nrays += 1;
        auto isec = intersect_scene(scene, bvh, ray);
        if (isec.instance_id < 0) return isec;
        auto& instance = scene.instances[isec.instance_id];
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[shape.material];
        auto  op       = evaluate_material_opacity(scene, material,
            evaluate_shape_texturecoord(shape, isec.element_id, isec.element_uv),
            evaluate_shape_color(shape, isec.element_id, isec.element_uv));
        if (op > 0.999f) return isec;
        if (get_random_float(rng) < op) return isec;
        ray = make_ray(evaluate_instance_position(
                           scene, instance, isec.element_id, isec.element_uv),
            ray.direction);
    }
    return {};
}

// Sample camera
ray3f sample_camera_ray(const yocto_camera& camera, const vec2i& ij,
    const vec2i& image_size, rng_state& rng) {
    auto puv = get_random_vec2f(rng);  // force order of evaluation with
                                       // assignments
    auto luv = get_random_vec2f(rng);  // force order of evaluation with
                                       // assignments
    return evaluate_camera_ray(camera, ij, image_size, puv, luv);
}

// Check if we are near the mirror direction.
inline bool check_near_mirror(vec3f normal, vec3f outgoing, vec3f incoming) {
    return fabs(dot(incoming,
                    normalize(normal * 2.0f * dot(outgoing, normal) - outgoing)) -
                1) < 0.001f;
}

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& incoming) {
    if (ks == zero_vec3f) return zero_vec3f;
    return ks + (vec3f{1, 1, 1} - ks) *
                    pow(clamp(1.0f - fabs(dot(h, incoming)), 0.0f, 1.0f), 5.0f);
}
vec3f fresnel_schlick(
    const vec3f& ks, const vec3f& h, const vec3f& incoming, float rs) {
    if (ks == zero_vec3f) return zero_vec3f;
    auto fks = fresnel_schlick(ks, fabs(dot(h, incoming)));
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX distribution and geometric term
float evaluate_ggx_dist(float rs, const vec3f& normal, const vec3f& h) {
    auto di = (dot(normal, h) * dot(normal, h)) * (rs * rs - 1) + 1;
    return rs * rs / (pif * di * di);
}
float evaluate_ggx_sm(float rs, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
#if 0
    // evaluate G from Heitz
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    auto Go = (2 * fabs(dot(normal, outgoing))) /
              (fabs(dot(normal, outgoing)) +
                  sqrt(rs * rs + (1 - rs * rs) * dot(normal, outgoing) *
                                     dot(normal, outgoing)));
    auto Gi = (2 * fabs(dot(normal, incoming))) /
              (fabs(dot(normal, incoming)) +
                  sqrt(rs * rs + (1 - rs * rs) * dot(normal, incoming) *
                                     dot(normal, incoming)));
    return Go * Gi;
#endif
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f evaluate_smooth_brdf_cosine(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) return zero_vec3f;
    auto brdf_cosine = zero_vec3f;

    // diffuse
    if (brdf.diffuse != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto h = normalize(incoming + outgoing);
        auto F = fresnel_schlick(brdf.specular, h, outgoing);
        brdf_cosine += brdf.diffuse * (vec3f{1, 1, 1} - F) / pif;
    }

    // specular
    if (brdf.specular != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto h = normalize(incoming + outgoing);
        auto F = fresnel_schlick(brdf.specular, h, outgoing);
        auto D = evaluate_ggx_dist(brdf.roughness, normal, h);
        auto G = evaluate_ggx_sm(brdf.roughness, normal, outgoing, incoming);
        brdf_cosine += F * D * G /
                       (4 * fabs(dot(normal, outgoing)) *
                           fabs(dot(normal, incoming)));
    }

    // transmission (thin sheet)
    if (brdf.transmission != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto ir = (dot(normal, outgoing) >= 0) ? reflect(-incoming, normal) :
                                                 reflect(-incoming, -normal);
        auto h = normalize(ir + outgoing);
        auto F = fresnel_schlick(brdf.specular, h, outgoing);
        auto D = evaluate_ggx_dist(brdf.roughness, normal, h);
        auto G = evaluate_ggx_sm(brdf.roughness, normal, outgoing, ir);
        brdf_cosine += brdf.transmission * (vec3f{1, 1, 1} - F) * D * G /
                       (4 * fabs(dot(normal, outgoing)) * fabs(dot(normal, ir)));
    }

    return brdf_cosine * abs(dot(normal, incoming));
}

// Evaluates the BRDF assuming that it is called only from the directions
// generated by sample_smooth_brdf_direction.
vec3f evaluate_delta_brdf_cosine(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!is_brdf_delta(brdf)) return zero_vec3f;
    auto microfacet_brdf = zero_vec3f;

    // specular
    if (brdf.specular != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto F = fresnel_schlick(brdf.specular, normal, outgoing);
        microfacet_brdf += F;
    }

    // transmission (thin sheet)
    if (brdf.transmission != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto F = fresnel_schlick(brdf.specular, normal, outgoing);
        microfacet_brdf += brdf.transmission * (1 - F);
    }

    return microfacet_brdf;
}

// Evalates the BRDF value. For delta BRDFs this assumes that it is called only
// from the directions generated by sample_delta_brdf_direction.
vec3f evaluate_brdf_cosine(const microfacet_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) {
        return evaluate_delta_brdf_cosine(brdf, normal, outgoing, incoming);
    } else {
        return evaluate_smooth_brdf_cosine(brdf, normal, outgoing, incoming);
    }
}

// Picks a direction based on the BRDF
vec3f sample_smooth_brdf_direction(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn) {
    if (is_brdf_delta(brdf)) return zero_vec3f;
    auto F    = fresnel_schlick(brdf.specular, normal, outgoing);
    auto prob = vec3f{max(brdf.diffuse * (vec3f{1, 1, 1} - F)), max(F),
        max(brdf.transmission * (vec3f{1, 1, 1} - F))};
    if (prob == zero_vec3f) return zero_vec3f;
    prob /= prob[0] + prob[1] + prob[2];

    // sample according to diffuse
    if (brdf.diffuse != zero_vec3f && rnl < prob[0]) {
        auto rz = sqrtf(rn[1]), rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn[0];
        auto il = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        auto fp = dot(normal, outgoing) >= 0 ?
                      make_frame_fromz(zero_vec3f, normal) :
                      make_frame_fromz(zero_vec3f, -normal);
        return transform_direction(fp, il);
    }
    // sample according to specular GGX
    else if (brdf.specular != zero_vec3f && rnl < prob[0] + prob[1]) {
        auto hl = sample_ggx(brdf.roughness, rn);
        auto fp = dot(normal, outgoing) >= 0 ?
                      make_frame_fromz(zero_vec3f, normal) :
                      make_frame_fromz(zero_vec3f, -normal);
        auto h = transform_direction(fp, hl);
        return reflect(outgoing, h);
    }
    // transmission hack
    else if (brdf.transmission != zero_vec3f &&
             rnl < prob[0] + prob[1] + prob[2]) {
        auto hl = sample_ggx(brdf.roughness, rn);
        auto fp = dot(normal, outgoing) >= 0 ?
                      make_frame_fromz(zero_vec3f, normal) :
                      make_frame_fromz(zero_vec3f, -normal);
        auto h  = transform_direction(fp, hl);
        auto ir = reflect(outgoing, h);
        return dot(normal, outgoing) >= 0 ? reflect(-ir, -normal) :
                                            reflect(-ir, normal);
    } else {
        return zero_vec3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_smooth_brdf_direction(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, rng_state& rng) {
    auto rnl = get_random_float(rng);  // force order of evaluation with
                                       // assignments
    auto rni = get_random_vec2f(rng);  // force order of evaluation with
                                       // assignments
    return sample_smooth_brdf_direction(brdf, normal, outgoing, rnl, rni);
}

// Picks a direction based on the BRDF
vec3f sample_delta_brdf_direction(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn) {
    if (!is_brdf_delta(brdf)) return zero_vec3f;
    auto F    = fresnel_schlick(brdf.specular, normal, outgoing);
    auto prob = vec3f{0, max(F), max(brdf.transmission * (vec3f{1, 1, 1} - F))};
    if (prob == zero_vec3f) return zero_vec3f;
    prob /= prob[0] + prob[1] + prob[2];

    // sample according to specular mirror
    if (brdf.specular != zero_vec3f && rnl < prob[0] + prob[1]) {
        return reflect(outgoing, dot(normal, outgoing) >= 0 ? normal : -normal);
    }
    // sample according to transmission
    else if (brdf.transmission != zero_vec3f && !brdf.refract &&
             rnl < prob[0] + prob[1] + prob[2]) {
        return -outgoing;
    }
    // sample according to transmission
    else if (brdf.transmission != zero_vec3f && brdf.refract &&
             rnl < prob[0] + prob[1] + prob[2]) {
        if (dot(normal, outgoing) >= 0) {
            return refract(outgoing, normal, 1 / specular_to_eta(brdf.specular));
        } else {
            return refract(outgoing, -normal, specular_to_eta(brdf.specular));
        }
    }
    // no sampling
    else {
        return zero_vec3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_delta_brdf_direction(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, rng_state& rng) {
    auto rnl = get_random_float(rng);  // force order of evaluation with
                                       // assignments
    auto rni = get_random_vec2f(rng);  // force order of evaluation with
                                       // assignments
    return sample_delta_brdf_direction(brdf, normal, outgoing, rnl, rni);
}

// Picks a direction based on the BRDF
vec3f sample_brdf_direction(const microfacet_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, rng_state& rng) {
    if (is_brdf_delta(brdf)) {
        return sample_delta_brdf_direction(brdf, normal, outgoing, rng);
    } else {
        return sample_smooth_brdf_direction(brdf, normal, outgoing, rng);
    }
}

// Compute the weight for sampling the BRDF
float sample_smooth_brdf_direction_pdf(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) return 0;
    auto F    = fresnel_schlick(brdf.specular, normal, outgoing);
    auto prob = vec3f{max(brdf.diffuse * (vec3f{1, 1, 1} - F)), max(F),
        max(brdf.transmission * (vec3f{1, 1, 1} - F))};
    if (prob == zero_vec3f) return 0;
    prob /= prob[0] + prob[1] + prob[2];

    auto pdf = 0.0f;

    if (brdf.diffuse != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        pdf += prob[0] * fabs(dot(normal, incoming)) / pif;
    }
    if (brdf.specular != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto h = normalize(incoming + outgoing);
        auto d = sample_ggx_pdf(brdf.roughness, fabs(dot(normal, h)));
        pdf += prob[1] * d / (4 * fabs(dot(outgoing, h)));
    }
    if (brdf.transmission != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto ir = (dot(normal, outgoing) >= 0) ? reflect(-incoming, normal) :
                                                 reflect(-incoming, -normal);
        auto h = normalize(ir + outgoing);
        auto d = sample_ggx_pdf(brdf.roughness, fabs(dot(normal, h)));
        pdf += prob[2] * d / (4 * fabs(dot(outgoing, h)));
    }

    return pdf;
}

// Compute the weight for sampling the BRDF
float sample_delta_brdf_direction_pdf(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!is_brdf_delta(brdf)) return 0;
    auto F    = fresnel_schlick(brdf.specular, normal, outgoing);
    auto prob = vec3f{0, max(F), max(brdf.transmission * (vec3f{1, 1, 1} - F))};
    if (prob == zero_vec3f) return 0;
    prob /= prob[0] + prob[1] + prob[2];

    auto pdf = 0.0f;

    if (brdf.specular != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        return prob[1];
    }
    if (brdf.transmission != zero_vec3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        return prob[2];
    }

    return pdf;
}

// Compute the weight for BRDF smapling
float sample_brdf_direction_pdf(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) {
        return sample_delta_brdf_direction_pdf(brdf, normal, outgoing, incoming);
    } else {
        return sample_smooth_brdf_direction_pdf(brdf, normal, outgoing, incoming);
    }
}

// Picks a point on a light.
trace_point sample_instance_point(const yocto_scene& scene,
    const trace_lights& lights, int instance_id, float rel, const vec2f& ruv) {
    auto& instance = scene.instances[instance_id];
    auto  sample   = pair<int, vec2f>();
    if (instance.shape >= 0) {
        auto& shape        = scene.shapes[instance.shape];
        auto& elements_cdf = lights.shape_elements_cdf[instance.shape];
        sample = sample_shape_element(shape, elements_cdf, rel, ruv);
    } else if (instance.surface >= 0) {
        auto& surface      = scene.surfaces[instance.surface];
        auto& elements_cdf = lights.surface_elements_cdf[instance.shape];
        sample = sample_surface_element(surface, elements_cdf, rel, ruv);
    } else {
        log_error("empty instance");
    }
    return make_trace_point(scene, instance_id, sample.first, sample.second);
}

// Sample pdf for a light point.
float sample_instance_point_pdf(const yocto_scene& scene,
    const trace_lights& lights, int instance_id, int element_id,
    const vec2f& element_uv) {
    auto& instance = scene.instances[instance_id];
    if (!is_instance_emissive(scene, instance)) return 0;
    if (instance.shape >= 0) {
        auto& shape        = scene.shapes[instance.shape];
        auto& elements_cdf = lights.shape_elements_cdf[instance.shape];
        return sample_shape_element_pdf(
            shape, elements_cdf, element_id, element_uv);
    } else if (instance.surface >= 0) {
        auto& surface      = scene.shapes[instance.surface];
        auto& elements_cdf = lights.surface_elements_cdf[instance.surface];
        return sample_shape_element_pdf(
            surface, elements_cdf, element_id, element_uv);
    } else {
        return 0;
    }
}

// Sample a point from all shape lights.
trace_point sample_lights_point(const yocto_scene& scene,
    const trace_lights& lights, const vec3f& position, rng_state& rng) {
    if (lights.instances.empty()) return {};
    auto light_id = sample_uniform_index(
        lights.instances.size(), get_random_float(rng));
    auto instance_id = lights.instances[light_id];
    auto rel         = get_random_float(rng);  // force order of evaluation
    auto ruv         = get_random_vec2f(rng);  // force order of evaluation
    auto point = sample_instance_point(scene, lights, instance_id, rel, ruv);
    auto direction = normalize(position - point.position);
    point.normal   = evaluate_instance_shading_normal(scene,
        scene.instances[instance_id], point.element_id, point.element_uv,
        direction);
    return point;
}

// Sample pdf for a light point.
float sample_lights_point_pdf(const yocto_scene& scene,
    const trace_lights& lights, const vec3f& position,
    const trace_point& light_point) {
    auto& instance = scene.instances[light_point.instance_id];
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[shape.material];
    if (lights.instances.empty()) return 0;
    if (material.emission == zero_vec3f) return 0;
    return sample_instance_point_pdf(scene, lights, light_point.instance_id,
               light_point.element_id, light_point.element_uv) *
           sample_uniform_index_pdf(lights.instances.size());
}

// Sample pdf for an environment.
float sample_environment_direction_pdf(const yocto_scene& scene,
    const trace_lights& lights, int environment_id, const vec3f& incoming) {
    auto& environment = scene.environments[environment_id];
    if (environment.emission_texture >= 0) {
        auto& elements_cdf = lights.environment_texture_cdf[environment.emission_texture];
        auto& emission_texture = scene.textures[environment.emission_texture];
        auto  size             = evaluate_texture_size(emission_texture);
        auto texcoord = evaluate_environment_texturecoord(environment, incoming);
        auto i        = (int)(texcoord[0] * size[0]);
        auto j        = (int)(texcoord[1] * size[1]);
        auto idx      = j * size[0] + i;
        auto prob     = sample_discrete_distribution_pdf(elements_cdf, idx) /
                    elements_cdf.back();
        auto angle = (2 * pif / size[0]) * (pif / size[1]) *
                     sin(pif * (j + 0.5f) / size[1]);
        return prob / angle;
    } else {
        return 1 / (4 * pif);
    }
}

// Picks a point on an environment.
vec3f sample_environment_direction(const yocto_scene& scene,
    const trace_lights& lights, int environment_id, float rel, const vec2f& ruv) {
    auto& environment = scene.environments[environment_id];
    if (environment.emission_texture >= 0) {
        auto& elements_cdf = lights.environment_texture_cdf[environment.emission_texture];
        auto& emission_texture = scene.textures[environment.emission_texture];
        auto  idx  = sample_discrete_distribution(elements_cdf, rel);
        auto  size = evaluate_texture_size(emission_texture);
        auto  u    = (idx % size[0] + 0.5f) / size[0];
        auto  v    = (idx / size[0] + 0.5f) / size[1];
        return evaluate_environment_direction(environment, {u, v});
    } else {
        return sample_sphere_direction(ruv);
    }
}

vec3f sample_environment_direction(const yocto_scene& scene,
    const trace_lights& lights, int environment_id, rng_state& rng) {
    auto rel = get_random_float(rng);  // force order of evaluation with
                                       // assignments
    auto ruv = get_random_vec2f(rng);  // force order of evaluation with
                                       // assignments
    return sample_environment_direction(scene, lights, environment_id, rel, ruv);
}

// Picks a point on a light.
vec3f sample_instance_direction(const yocto_scene& scene,
    const trace_lights& lights, int instance_id, const vec3f& p, float rel,
    const vec2f& ruv) {
    auto& instance = scene.instances[instance_id];
    auto  sample   = pair<int, vec2f>();
    if (instance.shape >= 0) {
        auto& shape        = scene.shapes[instance.shape];
        auto& elements_cdf = lights.shape_elements_cdf[instance.shape];
        sample = sample_shape_element(shape, elements_cdf, rel, ruv);
    } else if (instance.surface >= 0) {
        auto& surface      = scene.surfaces[instance.surface];
        auto& elements_cdf = lights.surface_elements_cdf[instance.surface];
        sample = sample_surface_element(surface, elements_cdf, rel, ruv);
    } else {
    }
    return normalize(
        evaluate_instance_position(scene, instance, sample.first, sample.second) -
        p);
}

// Sample pdf for a light point.
float sample_instance_direction_pdf(const yocto_scene& scene,
    const trace_lights& lights, int instance_id, const bvh_scene& bvh,
    const vec3f& position_, const vec3f& direction) {
    auto& instance = scene.instances[instance_id];
    if (!is_instance_emissive(scene, instance)) return 0;
    auto& elements_cdf = instance.shape >= 0 ?
                             lights.shape_elements_cdf[instance.shape] :
                             lights.surface_elements_cdf[instance.surface];
    // check all intersection
    auto pdf      = 0.0f;
    auto position = position_;
    for (auto bounce = 0; bounce < 10; bounce++) {
        auto isec = intersect_scene(
            scene, instance_id, bvh, make_ray(position, direction));
        if (isec.instance_id < 0) break;
        // accumulate pdf
        auto& instance       = scene.instances[isec.instance_id];
        auto  light_position = evaluate_instance_position(
            scene, instance, isec.element_id, isec.element_uv);
        auto light_normal = evaluate_instance_shading_normal(
            scene, instance, isec.element_id, isec.element_uv, direction);
        // prob triangle * area triangle = area triangle mesh
        auto area = elements_cdf.back();
        pdf += distance_squared(light_position, position) /
               (abs(dot(light_normal, direction)) * area);
        // continue
        position = light_position + direction * 1e-3f;
    }
    return pdf;
}

// Sample lights wrt solid angle
vec3f sample_lights_direction(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh, const vec3f& position,
    rng_state& rng) {
    auto light_id = sample_uniform_index(
        lights.instances.size() + lights.environments.size(),
        get_random_float(rng));
    if (light_id < lights.instances.size()) {
        auto instance = lights.instances[light_id];
        auto rel      = get_random_float(rng);  // force order of evaluation
        auto ruv      = get_random_vec2f(rng);  // force order of evaluation
        return sample_instance_direction(
            scene, lights, instance, position, rel, ruv);
    } else {
        auto environment = lights.environments[light_id -
                                               (int)lights.instances.size()];
        auto rel         = get_random_float(rng);  // force order of evaluation
        auto ruv         = get_random_vec2f(rng);  // force order of evaluation
        return sample_environment_direction(scene, lights, environment, rel, ruv);
    }
}

// Sample lights pdf
float sample_lights_direction_pdf(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh, const vec3f& position,
    const vec3f& direction) {
    auto pdf = 0.0f;
    for (auto instance : lights.instances) {
        pdf += sample_instance_direction_pdf(
            scene, lights, instance, bvh, position, direction);
    }
    for (auto environment : lights.environments) {
        pdf += sample_environment_direction_pdf(
            scene, lights, environment, direction);
    }
    pdf *= sample_uniform_index_pdf(
        lights.instances.size() + lights.environments.size());
    return pdf;
}

// Sample a direction accoding to either ligts or brdf
vec3f sample_lights_or_brdf_direction(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh,
    const microfacet_brdf& brdf, const vec3f& position, const vec3f& normal,
    const vec3f& outgoing, rng_state& rng) {
    auto rmode = get_random_float(rng);
    if (rmode < 0.5f) {
        return sample_lights_direction(scene, lights, bvh, position, rng);
    } else {
        return sample_smooth_brdf_direction(brdf, normal, outgoing, rng);
    }
}

// Pdf for direction sampling
float sample_lights_or_brdf_direction_pdf(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh,
    const microfacet_brdf& brdf, const vec3f& position, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    return 0.5f * sample_lights_direction_pdf(
                      scene, lights, bvh, position, incoming) +
           0.5f * sample_smooth_brdf_direction_pdf(
                      brdf, normal, outgoing, incoming);
}

// Russian roulette
bool sample_russian_roulette(
    const vec3f& weight, int bounce, rng_state& rng, int min_bounce = 2) {
    if (bounce <= min_bounce) return false;
    auto rrprob = 1.0f - min(max(weight), 0.95f);
    return get_random_float(rng) < rrprob;
}
float sample_russian_roulette_pdf(
    const vec3f& weight, int bounce, int min_bounce = 2) {
    if (bounce <= min_bounce) return 1;
    auto rrprob = 1.0f - min(max(weight), 0.95f);
    return 1 - rrprob;
}

#if 0
// Test occlusion.
vec3f evaluate_transmission(const yocto_scene& scene, const bvh_scene& bvh,
    const vec3f& from, const vec3f& to, int max_bounces) {
    auto weight = vec3f{1, 1, 1};
    auto p      = from;
    for (auto bounce = 0; bounce < max_bounces; bounce++) {
        auto ray  = make_segment(p, to);
        auto isec = intersect_scene(scene, bvh, ray);
        if (!isec.instance) break;
        auto f = evaluate_brdf(isec.instance.shape, isec.element_id, isec.element_uv);
        auto op = evaluate_opacity(
            isec.instance, isec.element_id, isec.element_uv);
        weight *= brdf.transmission + vec3f{1 - op, 1 - op, 1 - op};
        if (weight == zero_vec3f) break;
        p = evaluate_position(isec.instance, isec.element_id, isec.element_uv);
    }
    return weight;
}

#endif

vec3f evaluate_transmission(const yocto_scene& scene,
    const yocto_material& material, const vec3f& from, const vec3f& dir,
    float distance, int channel, rng_state& rng) {
    auto& vd = material.volume_density;
    if (is_material_volume_homogeneus(material))
        return vec3f{exp(-distance * vd[0]), exp(-distance * vd[1]),
            exp(-distance * vd[2])};

    // ratio tracking
    auto tr = 1.0f, t = 0.0f;
    auto pos = from;
    while (true) {
        auto step = -log(1 - get_random_float(rng)) / vd[channel];
        t += step;
        if (t >= distance) break;
        pos += dir * step;
        auto density = material.volume_density;
        if (material.volume_density_texture >= 0) {
            auto& volume_density_texture = scene.voltextures[material.volume_density_texture];
            density *= evaluate_voltexture(volume_density_texture, pos);
        }
        tr *= 1.0f - max(0.0f, density[channel] / vd[channel]);
    }
    return {tr, tr, tr};
}

float sample_distance(const yocto_scene& scene, const yocto_material& material,
    const vec3f& from, const vec3f& dir, int channel, rng_state& rng) {
    auto pos      = from;
    auto majorant = material.volume_density[channel];
    if (majorant == 0) return float_max;

    // delta tracking
    auto distance = 0.0f;
    while (true) {
        auto r = get_random_float(rng);
        if (r == 0) return float_max;
        auto step = -log(r) / majorant;
        if (is_material_volume_homogeneus(material)) return step;

        pos += dir * step;
        distance += step;
        auto density = material.volume_density;
        if (material.volume_density_texture >= 0) {
            auto& volume_density_texture = scene.voltextures[material.volume_density_texture];
            density *= evaluate_voltexture(volume_density_texture, pos);
        }
        if (density[channel] / majorant >= get_random_float(rng))
            return distance;

        // Escape from volume.
        if (pos[0] > 1 || pos[1] > 1 || pos[2] > 1) return float_max;
        if (pos[0] < -1 || pos[1] < -1 || pos[2] < -1) return float_max;
    }
}

float sample_distance(const yocto_scene& scene, const yocto_instance& instance,
    const bbox3f& bbox, const vec3f& from, const vec3f& dir, int channel,
    rng_state& rng) {
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[shape.material];
    if (material.volume_density == zero_vec3f) return float_max;

    // Transform coordinates so that every position in the bounding box of the
    // instance is mapped to the cube [-1,1]^3 (the same space of volume texture
    // sampling).
    auto scale    = bbox.max - bbox.min;
    auto frame    = instance.frame;
    auto froml    = transform_point_inverse(frame, from) / scale;
    auto dirl     = transform_direction_inverse(frame, dir) / scale;
    auto ll       = length(dirl);
    auto distance = sample_distance(
        scene, material, froml, dirl / ll, channel, rng);
    return distance * ll;
}

vec3f sample_phase_function(float g, const vec2f& u) {
    auto cos_theta = 0.0f;
    if (abs(g) < 1e-3) {
        cos_theta = 1 - 2 * u[0];
    } else {
        float square = (1 - g * g) / (1 - g + 2 * g * u[0]);
        cos_theta    = (1 + g * g - square * square) / (2 * g);
    }

    auto sin_theta = sqrt(max(0.0f, 1 - cos_theta * cos_theta));
    auto phi       = 2 * pif * u[1];
    return {sin_theta * cos(phi), sin_theta * sin(phi), cos_theta};
}

float evaluate_phase_function(float cos_theta, float g) {
    auto denom = 1 + g * g + 2 * g * cos_theta;
    return (1 - g * g) / (4 * pif * denom * sqrt(denom));
}

// Probability of computing direct illumination.
float prob_direct(const microfacet_brdf& brdf) {
    // This is just heuristic. Any other choice is equally correct.
    if (brdf.diffuse + brdf.specular == zero_vec3f) return 0;
    auto kd = max(brdf.diffuse);
    auto ks = max(brdf.specular);
    return (kd + brdf.roughness * ks) / (kd + ks);
}

// Sample a direction of direct illumination from the point p, which is inside
// mediums.back(). pdf and incoming radiance le are returned in reference. It
// works for both surface rendering and volume rendering.
vec3f direct_illumination(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& p, int channel,
    const vector<int>& mediums_, rng_state& rng, float& pdf, vec3f& le) {
    auto  incoming = zero_vec3f;
    vec3f weight   = vec3f{1, 1, 1};
    auto  mediums  = mediums_;

    auto idx = sample_uniform_index(
        lights.instances.size() + lights.environments.size(),
        get_random_float(rng));
    pdf = 1.0f / (lights.instances.size() + lights.environments.size());
    if (idx < lights.instances.size()) {
        auto instance_id = lights.instances[idx];
        incoming = sample_instance_direction(scene, lights, instance_id, p,
            get_random_float(rng), get_random_vec2f(rng));
        auto& instance = scene.instances[instance_id];
        pdf *= 1.0 / lights.shape_elements_cdf[instance.shape].back();
    } else {
        auto environment_id = lights.environments[idx - lights.instances.size()];
        incoming = sample_environment_direction(scene, lights, environment_id,
            get_random_float(rng), get_random_vec2f(rng));
        pdf *= sample_environment_direction_pdf(
            scene, lights, environment_id, incoming);
        auto isec = intersect_scene_with_opacity(
            scene, bvh, make_ray(p, incoming), rng, 10);
        if (isec.instance_id < 0) {
            auto& environment = scene.environments[environment_id];
            le = evaluate_environment_emission(scene, environment, incoming);
            return incoming;
        }
    }

    auto isec = intersect_scene(scene, bvh, make_ray(p, incoming));

    while (isec.instance_id >= 0) {
        auto& isec_instance = scene.instances[isec.instance_id];
        auto& isec_shape    = scene.shapes[isec_instance.shape];
        auto  lp            = evaluate_instance_position(
            scene, isec_instance, isec.element_id, isec.element_uv);
        auto ln = evaluate_instance_shading_normal(
            scene, isec_instance, isec.element_id, isec.element_uv, -incoming);
        auto& isec_material = scene.materials[isec_shape.material];
        auto  emission      = evaluate_material_emission(scene, isec_material,
            evaluate_shape_texturecoord(
                isec_shape, isec.element_id, isec.element_uv),
            evaluate_shape_color(isec_shape, isec.element_id, isec.element_uv));

        auto& medium_instance = scene.instances[mediums.back()];
        auto& medium_shape    = scene.shapes[medium_instance.shape];
        auto& medium_material = scene.materials[medium_shape.material];
        if (medium_material.volume_density != zero_vec3f)
            weight *= evaluate_transmission(scene, medium_material, lp,
                incoming, isec.distance, channel, rng);

        // Hack: Uncomment this or the result will be biased
        // If mediums refracts, the transmission ray won't reach the sampled
        // light point if(isec.instance.mat->refract) break;

        if (emission != zero_vec3f) {
            // Geometric term.
            weight *= fabs(dot(ln, incoming)) / dot(lp - p, lp - p);
            le += weight * emission;
            break;
        }

        auto brdf = evaluate_material_brdf(scene, isec_material,
            evaluate_shape_texturecoord(
                isec_shape, isec.element_id, isec.element_uv),
            evaluate_shape_color(isec_shape, isec.element_id, isec.element_uv));
        if (brdf.transmission == zero_vec3f) {
            le = zero_vec3f;
            break;
        }

        auto  ndi       = dot(incoming, ln);
        float threshold = 0.05;

        if (ndi > threshold) {
            // Exiting from medium.
            if (isec.instance_id != mediums.back()) {  // exiting a different
                                                       // medium??
                pdf = 0;
                return zero_vec3f;
            }
            if (mediums.size() <= 1) {
                pdf = 0;
                return zero_vec3f;
            }
            mediums.pop_back();
        } else if (ndi < -threshold) {
            // Entering new medium.
            if (isec.instance_id == mediums.back()) {  // entering the same
                                                       // medium??
                pdf = 0;
                return zero_vec3f;
            }
            mediums.push_back(isec.instance_id);
        } else {
            pdf = 0;
            return zero_vec3f;
        }

        isec = intersect_scene(scene, bvh, make_ray(lp, incoming));  //@Hack:
                                                                     // 10? Don't
                                                                     // know...
    }

    return incoming;
}

// Recursive path tracing.
pair<vec3f, bool> trace_path(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto weight   = vec3f{1, 1, 1};
    auto outgoing = -direction;

    // trace  path
    for (auto bounce = 0; bounce < max_bounces; bounce++) {
        // direct
        if (!is_brdf_delta(point.brdf) && !empty(lights)) {
            auto light_direction = sample_lights_or_brdf_direction(scene, lights,
                bvh, point.brdf, point.position, point.normal, outgoing, rng);
            auto light_pdf = sample_lights_or_brdf_direction_pdf(scene, lights,
                bvh, point.brdf, point.position, point.normal, outgoing,
                light_direction);
            auto light_point = trace_ray_with_opacity(
                scene, bvh, point.position, light_direction, rng, max_bounces);
            auto brdf_cosine = evaluate_smooth_brdf_cosine(
                point.brdf, point.normal, outgoing, light_direction);
            if (light_pdf)
                radiance += weight * light_point.emission * brdf_cosine /
                            light_pdf;
        }

        // continue path
        auto next_direction = sample_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);

        // accumulate weight
        if (next_pdf == 0) break;
        weight *= brdf_cosine / next_pdf;
        if (weight == zero_vec3f) break;

        // russian roulette
        if (sample_russian_roulette(weight, bounce, rng)) break;
        weight /= sample_russian_roulette_pdf(weight, bounce);

        // intersect next point
        auto next_point = trace_ray_with_opacity(
            scene, bvh, point.position, next_direction, rng, max_bounces);
        if (is_brdf_delta(point.brdf)) radiance += weight * next_point.emission;
        if (next_point.instance_id < 0 || is_brdf_zero(next_point.brdf)) break;

        // setup next iteration
        point    = next_point;
        outgoing = -next_direction;
    }

    return {radiance, true};
}

// Evaluates the weight after sampling distance in a medium.
vec3f evaluate_transmission_div_pdf(const vec3f& vd, float distance, int ch) {
    auto weight = zero_vec3f;

    // For the sampled channel, transmission / pdf == 1.0
    weight[ch] = 1.0;

    // Compute weight for the remaining channels i.
    // In order to avoid numerical nasties (NaNs) transmission / pdf is
    // evaluated. transmission[i] = exp(-distance * vd[i]) pdf             =
    // exp(-distance * vd[channel])
    int i = (ch + 1) % 3, j = (ch + 2) % 3;
    weight[i] = exp(-distance * (vd[i] - vd[ch]));
    weight[j] = exp(-distance * (vd[j] - vd[ch]));
    return weight;
}

// Iterative volume path tracing.
pair<vec3f, bool> trace_volpath(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    if (empty(lights)) return {zero_vec3f, false};

    // initialize
    auto radiance = zero_vec3f;
    auto weight   = vec3f{1, 1, 1};
    auto emission = true;
    auto ray      = make_ray(position, direction);

#if 0
    // @Hack: air volume properties should be set in the scene struct.
    if (air == nullptr) {
        air = new instance();
        air->name = "air";
        air->mat = new material();
        air->mat->vd = vec3f{0.0, 0.0, 0.0};
        air->mat->va = vec3f{0.0, 0.0, 0.0};
        air->mat->vg = 0.0;
    }
#endif

    // List of mediums that contains the path. The path starts in air.
    auto mediums = vector<int>{-1};

    // Sample color channel. This won't matter if there are no heterogeneus
    // materials.
    auto ch             = sample_uniform_index(3, get_random_float(rng));
    auto single_channel = false;

    int bounce = 0;
    while (bounce < max_bounces) {
        auto  medium   = mediums.back();
        auto& instance = scene.instances[medium];
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[shape.material];
        auto  ve       = material.volume_emission;
        auto  va       = material.volume_albedo;
        auto  vd       = material.volume_density;
        auto  vg       = material.volume_phaseg;

        // If medium has color but must use delta tracking, integrate only the
        // sampled spectrum.
        if (!single_channel && is_material_volume_colored(material) &&
            !is_material_volume_homogeneus(material)) {
            weight[ch] *= 3;
            weight[(ch + 1) % 3] = 0;
            weight[(ch + 2) % 3] = 0;
            single_channel       = true;
        }

        // TODO: FIXME REMOVING BBOX
        // Sample distance of next absorption/scattering event in the medium.
        // dist_pdf is unknown due to delta tracking.
        auto bbox = transform_bbox(
            instance.frame, bbox3f{{-1, -1, -1}, {1, 1, 1}});
        auto distance = sample_distance(
            scene, instance, bbox, ray.origin, ray.direction, ch, rng);

        // Create ray and clamp it to make the intersection faster.
        ray       = make_ray(ray.origin, ray.direction);
        ray.tmax  = distance;
        auto isec = intersect_scene_with_opacity(
            scene, bvh, ray, rng, max_bounces);

        // @Hack: When isec.instance == nullptr, we must discern if the ray hit
        // nothing (the environment)
        //        or a medium interaction was sampled. Doing isec.distance ==
        //        maxf doesn't work, why??
        auto scene_size = max(bvh.nodes[0].bbox.max - bvh.nodes[0].bbox.min);

        // environment
        if (isec.instance_id < 0 && distance > scene_size) {
            if (emission) {
                for (auto& environment : scene.environments)
                    radiance += weight * evaluate_environment_emission(
                                             scene, environment, ray.direction);
            }
            return {radiance, false};
        }

        // surface intersection
        if (isec.instance_id >= 0) {
            auto& isec_instance = scene.instances[isec.instance_id];
            auto& isec_shape    = scene.shapes[isec_instance.shape];
            auto  outgoing      = -ray.direction;
            auto& isec_material = scene.materials[isec_shape.material];
            auto  p             = evaluate_instance_position(
                scene, isec_instance, isec.element_id, isec.element_uv);
            auto normal = evaluate_instance_shading_normal(scene, isec_instance,
                isec.element_id, isec.element_uv, outgoing);
            auto brdf   = evaluate_material_brdf(scene, isec_material,
                evaluate_shape_texturecoord(
                    isec_shape, isec.element_id, isec.element_uv),
                evaluate_shape_color(
                    isec_shape, isec.element_id, isec.element_uv));

            // distance sampling pdf is unknown due to delta tracking, but we do
            // know the value of transmission / pdf_dist.
            weight *= evaluate_transmission_div_pdf(vd, isec.distance, ch);

            // emission
            if (emission)
                radiance += weight * evaluate_material_emission(scene,
                                         isec_material,
                                         evaluate_shape_texturecoord(isec_shape,
                                             isec.element_id, isec.element_uv),
                                         evaluate_shape_color(isec_shape,
                                             isec.element_id, isec.element_uv));

            // early exit
            if (brdf.diffuse + brdf.specular + brdf.transmission == zero_vec3f ||
                bounce >= max_bounces - 1)
                break;

            // direct lighting
            if (get_random_float(rng) < prob_direct(brdf)) {
                // With some probabilty, this is a naive path tracer (works
                // great with delta-like brdfs)
                vec3f direct;
                float pdf;
                vec3f incoming = direct_illumination(
                    scene, bvh, lights, p, ch, mediums, rng, pdf, direct);
                if (pdf != 0) {
                    auto brdf_cosine = evaluate_smooth_brdf_cosine(
                        brdf, normal, outgoing, incoming);
                    radiance += weight * direct * brdf_cosine / pdf;
                    emission = false;
                }
            } else
                emission = true;

            // continue path
            vec3f incoming, brdf_cosine;
            float pdf = 0;
            if (!is_brdf_delta(brdf)) {
                incoming = sample_smooth_brdf_direction(
                    brdf, normal, outgoing, rng);
                brdf_cosine = evaluate_smooth_brdf_cosine(
                    brdf, normal, outgoing, incoming);
                pdf = sample_smooth_brdf_direction_pdf(
                    brdf, normal, outgoing, incoming);
            } else {
                incoming = sample_delta_brdf_direction(
                    brdf, normal, outgoing, rng);
                brdf_cosine = evaluate_delta_brdf_cosine(
                    brdf, normal, outgoing, incoming);
                pdf = sample_delta_brdf_direction_pdf(
                    brdf, normal, outgoing, incoming);
            }
            auto ndi = dot(normal, incoming);
            auto ndo = dot(normal, outgoing);

            // accumulate weight
            if (pdf == 0) break;
            weight *= brdf_cosine / pdf;
            if (weight == zero_vec3f) break;
            ray.origin       = p;
            ray.direction    = incoming;
            bool transmitted = (ndi > 0) != (ndo > 0);

            // transmission in medium
            if (transmitted) {
                float tr = 0.05;  // avoid numerical errors
                if (ndo < -tr) {
                    // Exiting from medium.
                    if (isec.instance_id != medium) break;
                    if (mediums.size() <= 1) break;
                    mediums.pop_back();
                } else if (ndo > tr) {
                    // Entering new medium.
                    if (isec.instance_id == medium) break;
                    mediums.push_back(isec.instance_id);
                } else
                    break;
            }
            bounce += 1;
        }
        // medium interaction
        else {
            ray.origin += ray.direction * distance;
            float scattering_prob = va[ch];

            // absorption and emission
            if (get_random_float(rng) >= scattering_prob) {
                weight /= 1 - scattering_prob;
                radiance += weight * ve;
                break;
            }

            // scattering event
            weight /= scattering_prob;
            weight *= evaluate_transmission_div_pdf(vd, distance, ch);

            // direct lighting
            vec3f direct;
            float pdf_direct;
            vec3f l = direct_illumination(scene, bvh, lights, ray.origin, ch,
                mediums, rng, pdf_direct, direct);
            if (pdf_direct != 0) {
                auto f = va * evaluate_phase_function(dot(l, -ray.direction), vg);
                radiance += weight * direct * f / pdf_direct;
                emission = false;
            }

            // indirect
            vec3f incoming = sample_phase_function(vg, get_random_vec2f(rng));
            weight *= va;
            ray.direction = transform_direction(
                make_frame_fromz(zero_vec3f, ray.direction), incoming);
        }

        // russian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (get_random_float(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }
    }

    return {radiance, true};
}

// Recursive path tracing.
pair<vec3f, bool> trace_path_naive(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto weight   = vec3f{1, 1, 1};
    auto outgoing = -direction;

    // trace  path
    for (auto bounce = 0; bounce < max_bounces; bounce++) {
        // continue path
        auto next_direction = sample_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);

        // accumulate weight
        if (next_pdf == 0) break;
        weight *= brdf_cosine / next_pdf;
        if (weight == zero_vec3f) break;

        // russian roulette
        if (sample_russian_roulette(weight, bounce, rng)) break;
        weight /= sample_russian_roulette_pdf(weight, bounce);

        // intersect next point
        auto next_point = trace_ray_with_opacity(
            scene, bvh, point.position, next_direction, rng, max_bounces);
        radiance += weight * next_point.emission;
        if (next_point.instance_id < 0 || is_brdf_zero(next_point.brdf)) break;

        // setup next iteration
        point    = next_point;
        outgoing = -next_direction;
    }

    return {radiance, true};
}

// Recursive path tracing.
pair<vec3f, bool> trace_path_nomis(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto weight   = vec3f{1, 1, 1};
    auto outgoing = -direction;

    // trace  path
    for (auto bounce = 0; bounce < max_bounces; bounce++) {
        // direct
        if (!is_brdf_delta(point.brdf) && !empty(lights)) {
            auto light_point = sample_lights_point(
                scene, lights, point.position, rng);
            auto light_pdf = sample_lights_point_pdf(
                scene, lights, point.position, light_point);
            auto light_direction = normalize(
                light_point.position - point.position);
            auto intersection_point = trace_ray_with_opacity(
                scene, bvh, point.position, light_direction, rng, max_bounces);
            if (light_pdf &&
                light_point.instance_id == intersection_point.instance_id) {
                auto brdf_cosine = evaluate_smooth_brdf_cosine(
                    point.brdf, point.normal, outgoing, light_direction);
                auto geometric_term = abs(dot(light_point.normal,
                                          light_direction)) /
                                      distance_squared(
                                          light_point.position, point.position);
                radiance += weight * light_point.emission * brdf_cosine *
                            geometric_term / light_pdf;
            }
        }

        // continue path
        auto next_direction = sample_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);

        // accumulate weight
        if (next_pdf == 0) break;
        weight *= brdf_cosine / next_pdf;
        if (weight == zero_vec3f) break;

        // russian roulette
        if (sample_russian_roulette(weight, bounce, rng)) break;
        weight /= sample_russian_roulette_pdf(weight, bounce);

        // intersect next point
        auto next_point = trace_ray_with_opacity(
            scene, bvh, point.position, next_direction, rng, max_bounces);
        if (next_point.instance_id < 0 || is_brdf_delta(point.brdf))
            radiance += weight * next_point.emission;
        if (next_point.instance_id < 0 || is_brdf_zero(next_point.brdf)) break;

        // setup next iteration
        point    = next_point;
        outgoing = -next_direction;
    }

    return {radiance, true};
}

// Direct illumination.
pair<vec3f, bool> trace_direct(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto outgoing = -direction;

    // direct
    if (!is_brdf_delta(point.brdf) && !empty(lights)) {
        auto light_direction = sample_lights_or_brdf_direction(scene, lights,
            bvh, point.brdf, point.position, point.normal, outgoing, rng);
        auto light_pdf = sample_lights_or_brdf_direction_pdf(scene, lights, bvh,
            point.brdf, point.position, point.normal, outgoing, light_direction);
        auto light_point = trace_ray_with_opacity(
            scene, bvh, point.position, light_direction, rng, max_bounces);
        auto brdf_cosine = evaluate_smooth_brdf_cosine(
            point.brdf, point.normal, outgoing, light_direction);
        if (light_pdf)
            radiance += light_point.emission * brdf_cosine / light_pdf;
    }

    // deltas
    if (is_brdf_delta(point.brdf) && max_bounces) {
        auto next_direction = sample_delta_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_delta_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_delta_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);
        auto incoming_radiance = trace_direct(scene, bvh, lights,
            point.position, next_direction, rng, max_bounces - 1, true)
                                     .first;
        radiance += brdf_cosine * incoming_radiance / next_pdf;
    }

    // done
    return {radiance, true};
}

// Direct illumination.
pair<vec3f, bool> trace_direct_nomis(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto outgoing = -direction;

    // direct
    if (!is_brdf_delta(point.brdf) && !lights.instances.empty()) {
        auto light_point = sample_lights_point(
            scene, lights, point.position, rng);
        auto light_pdf = sample_lights_point_pdf(
            scene, lights, point.position, light_point);
        auto light_direction = normalize(light_point.position - point.position);
        auto intersection_point = trace_ray_with_opacity(
            scene, bvh, point.position, light_direction, rng, max_bounces);
        if (light_pdf &&
            light_point.instance_id == intersection_point.instance_id) {
            auto brdf_cosine = evaluate_smooth_brdf_cosine(
                point.brdf, point.normal, outgoing, light_direction);
            auto geometric_term = abs(dot(light_point.normal, light_direction)) /
                                  distance_squared(
                                      light_point.position, point.position);
            radiance += light_point.emission * brdf_cosine * geometric_term /
                        light_pdf;
        }
    }

    // environments
    if (!is_brdf_delta(point.brdf) && !lights.environments.empty()) {
        auto next_direction = sample_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);
        auto emission = evaluate_environment_emission(scene, next_direction);
        if (next_pdf &&
            intersect_scene_with_opacity(scene, bvh,
                make_ray(point.position, next_direction), rng, max_bounces)
                    .instance_id < 0)
            radiance += emission * brdf_cosine / next_pdf;
    }

    // deltas
    if (is_brdf_delta(point.brdf)) {
        auto next_direction = sample_delta_brdf_direction(
            point.brdf, point.normal, outgoing, rng);
        auto brdf_cosine = evaluate_delta_brdf_cosine(
            point.brdf, point.normal, outgoing, next_direction);
        auto next_pdf = sample_delta_brdf_direction_pdf(
            point.brdf, point.normal, outgoing, next_direction);
        auto incoming_radiance = trace_direct_nomis(scene, bvh, lights,
            point.position, next_direction, rng, max_bounces - 1, true)
                                     .first;
        radiance += brdf_cosine * incoming_radiance * next_pdf;
    }

    // done
    return {radiance, true};
}

// Environment illumination only with no shadows.
pair<vec3f, bool> trace_environment(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto outgoing = -direction;

    // continue path
    auto next_direction = sample_brdf_direction(
        point.brdf, point.normal, outgoing, rng);
    auto brdf_cosine = evaluate_brdf_cosine(
        point.brdf, point.normal, outgoing, next_direction);
    auto next_pdf = sample_brdf_direction_pdf(
        point.brdf, point.normal, outgoing, next_direction);

    // accumulate environment illumination
    if (next_pdf)
        radiance += brdf_cosine *
                    evaluate_environment_emission(scene, next_direction) /
                    next_pdf;

    // done
    return {radiance, true};
}

// Eyelight for quick previewing.
pair<vec3f, bool> trace_eyelight(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction,
    rng_state& rng, int max_bounces, bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) {
        if (environments_hidden || scene.environments.empty())
            return {zero_vec3f, false};
        return {point.emission, true};
    }

    // initialize
    auto radiance = point.emission;
    auto outgoing = -direction;

    // microfacet_brdf * light
    radiance += evaluate_smooth_brdf_cosine(
                    point.brdf, point.normal, outgoing, outgoing) *
                pif;

    // done
    return {radiance, true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_normal(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    return {point.normal * 0.5f + 0.5f, true};
}

// Debug frontfacing.
pair<vec3f, bool> trace_debug_frontfacing(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    auto outgoing    = -direction;
    auto frontfacing = dot(point.normal, outgoing) > 0 ? vec3f{0, 1, 0} :
                                                         vec3f{1, 0, 0};
    return {frontfacing, true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_albedo(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    auto albedo = point.brdf.diffuse + point.brdf.specular +
                  point.brdf.transmission;
    return {albedo, true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_diffuse(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    return {point.brdf.diffuse, true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_specular(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    return {point.brdf.specular, true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_roughness(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    return {{point.brdf.roughness, point.brdf.roughness, point.brdf.roughness},
        true};
}

// Debug previewing.
pair<vec3f, bool> trace_debug_texcoord(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces,
    bool environments_hidden) {
    // intersect ray
    auto point = trace_ray_with_opacity(
        scene, bvh, position, direction, rng, max_bounces);
    if (point.instance_id < 0) return {zero_vec3f, false};

    // shade
    return {{point.texturecoord[0], point.texturecoord[1], 0}, true};
}

// Trace a single ray from the camera using the given algorithm.
trace_sampler_func get_trace_sampler_func(trace_sampler_type type) {
    switch (type) {
        case trace_sampler_type::path: return trace_path;
        // case trace_type::volpath:
        //     return trace_volpath(
        //         scene, bvh, lights, position, direction, rng, max_bounces);
        case trace_sampler_type::direct: return trace_direct;
        case trace_sampler_type::environment: return trace_environment;
        case trace_sampler_type::eyelight: return trace_eyelight;
        case trace_sampler_type::path_nomis: return trace_path_nomis;
        case trace_sampler_type::path_naive: return trace_path_naive;
        case trace_sampler_type::direct_nomis: return trace_direct_nomis;
        case trace_sampler_type::debug_normal: return trace_debug_normal;
        case trace_sampler_type::debug_albedo: return trace_debug_albedo;
        case trace_sampler_type::debug_texcoord: return trace_debug_texcoord;
        case trace_sampler_type::debug_frontfacing:
            return trace_debug_frontfacing;
        case trace_sampler_type::debug_diffuse: return trace_debug_diffuse;
        case trace_sampler_type::debug_specular: return trace_debug_specular;
        case trace_sampler_type::debug_roughness: return trace_debug_roughness;
    }
    return {};
}

// Trace a block of samples
void trace_image_region(image<vec4f>& rendered_image, image<trace_pixel>& pixels,
    const yocto_scene& scene, const bvh_scene& bvh, const trace_lights& lights,
    const bbox2i& region, int num_samples, const trace_image_options& options) {
    auto& camera  = scene.cameras.at(options.camera_id);
    auto  sampler = get_trace_sampler_func(options.sampler_type);
    for (auto j = region.min[1]; j < region.max[1]; j++) {
        for (auto i = region.min[0]; i < region.max[0]; i++) {
            auto& pixel = pixels[{i, j}];
            for (auto s = 0; s < num_samples; s++) {
                if (options.cancel_flag && *options.cancel_flag) return;
                _trace_npaths += 1;
                auto ray = sample_camera_ray(
                    camera, {i, j}, rendered_image.size(), pixel.rng);
                auto radiance_hit = sampler(scene, bvh, lights, ray.origin,
                    ray.direction, pixel.rng, options.max_bounces,
                    options.environments_hidden);
                auto radiance     = radiance_hit.first;
                auto hit          = radiance_hit.second;
                if (!isfinite(radiance[0]) || !isfinite(radiance[1]) ||
                    !isfinite(radiance[2])) {
                    log_error("NaN detected");
                    radiance = zero_vec3f;
                }
                if (max(radiance) > options.pixel_clamp)
                    radiance = radiance * (options.pixel_clamp / max(radiance));
                pixel.radiance += radiance;
                pixel.hits += hit ? 1 : 0;
                pixel.samples += 1;
            }
            auto radiance = pixel.hits ? pixel.radiance / pixel.hits : zero_vec3f;
            auto coverage          = (float)pixel.hits / (float)pixel.samples;
            rendered_image[{i, j}] = {
                radiance[0], radiance[1], radiance[2], coverage};
        }
    }
}

// Init a sequence of random number generators.
image<trace_pixel> make_trace_pixels(const vec2i& image_size, uint64_t seed) {
    auto pixels = image<trace_pixel>{image_size};
    auto rng    = make_rng(1301081);
    for (auto j = 0; j < pixels.height(); j++) {
        for (auto i = 0; i < pixels.width(); i++) {
            auto& pixel = pixels[{i, j}];
            pixel.rng   = make_rng(seed, get_random_int(rng, 1 << 31) / 2 + 1);
        }
    }
    return pixels;
}

// Init trace lights
trace_lights make_trace_lights(const yocto_scene& scene) {
    auto scope  = log_trace_scoped("making trace lights");
    auto lights = trace_lights{};

    lights.shape_elements_cdf.resize(scene.shapes.size());
    lights.surface_elements_cdf.resize(scene.surfaces.size());
    lights.environment_texture_cdf.resize(scene.textures.size());

    for (auto instance_id = 0; instance_id < scene.instances.size();
         instance_id++) {
        auto& instance = scene.instances[instance_id];
        if (!is_instance_emissive(scene, instance)) continue;
        if (instance.shape >= 0) {
            auto& shape = scene.shapes[instance.shape];
            if (shape.triangles.empty() && shape.quads.empty()) continue;
            lights.instances.push_back(instance_id);
            lights.shape_elements_cdf[instance.shape] = compute_shape_elements_cdf(
                shape);
        } else if (instance.surface >= 0) {
            auto& surface = scene.surfaces[instance.surface];
            if (surface.quads_positions.empty()) continue;
            lights.instances.push_back(instance_id);
            lights.surface_elements_cdf[instance.surface] = compute_surface_elements_cdf(
                surface);
        } else {
            continue;
        }
    }

    for (auto environment_id = 0; environment_id < scene.environments.size();
         environment_id++) {
        auto& environment = scene.environments[environment_id];
        if (environment.emission == zero_vec3f) continue;
        lights.environments.push_back(environment_id);
        if (environment.emission_texture >= 0) {
            lights.environment_texture_cdf[environment.emission_texture] = compute_environment_texels_cdf(
                scene, environment);
        }
    }

    if (lights.instances.empty() && lights.environments.empty()) return {};
    return lights;
}

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const trace_image_options& options) {
    auto scope      = log_trace_scoped("tracing image");
    auto image_size = get_camera_image_size(
        scene.cameras.at(options.camera_id), options.image_size);
    auto rendered_image = image<vec4f>{image_size};
    auto pixels         = make_trace_pixels(image_size, options.random_seed);
    auto regions        = make_image_regions(rendered_image.size());

    if (options.run_serially) {
        for (auto& region : regions) {
            if (options.cancel_flag && *options.cancel_flag) break;
            trace_image_region(rendered_image, pixels, scene, bvh, lights,
                region, options.num_samples, options);
        }
    } else {
        auto nthreads = thread::hardware_concurrency();
        auto threads  = vector<thread>();
        for (auto tid = 0; tid < nthreads; tid++) {
            threads.push_back(thread([&, tid]() {
                for (auto region_id = tid; region_id < regions.size();
                     region_id += nthreads) {
                    if (options.cancel_flag && *options.cancel_flag) break;
                    auto& region = regions[region_id];
                    trace_image_region(rendered_image, pixels, scene, bvh,
                        lights, region, options.num_samples, options);
                }
            }));
        }
        for (auto& t : threads) t.join();
    }
    return rendered_image;
}

// Progressively compute an image by calling trace_samples multiple times.
int trace_image_samples(image<vec4f>& rendered_image, image<trace_pixel>& pixels,
    const yocto_scene& scene, const bvh_scene& bvh, const trace_lights& lights,
    int current_sample, const trace_image_options& options) {
    auto regions = make_image_regions(rendered_image.size());
    auto scope   = log_trace_scoped(
        "tracing samples {}-{}", current_sample, options.num_samples);
    auto num_samples = min(
        options.samples_per_batch, options.num_samples - current_sample);
    if (options.run_serially) {
        for (auto& region : regions) {
            if (options.cancel_flag && *options.cancel_flag) break;
            trace_image_region(rendered_image, pixels, scene, bvh, lights,
                region, num_samples, options);
        }
    } else {
        auto nthreads = thread::hardware_concurrency();
        auto threads  = vector<thread>();
        for (auto tid = 0; tid < nthreads; tid++) {
            threads.push_back(thread([&, tid]() {
                for (auto region_id = tid; region_id < regions.size();
                     region_id += nthreads) {
                    if (options.cancel_flag && *options.cancel_flag) break;
                    auto& region = regions[region_id];
                    trace_image_region(rendered_image, pixels, scene, bvh,
                        lights, region, num_samples, options);
                }
            }));
        }
        for (auto& t : threads) t.join();
    }
    return current_sample + num_samples;
}

// Starts an anyncrhounous renderer.
void trace_image_async_start(image<vec4f>& rendered_image,
    image<trace_pixel>& pixels, const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, vector<thread>& threads,
    atomic<int>& current_sample, concurrent_queue<bbox2i>& queue,
    const trace_image_options& options) {
    log_trace("start tracing async");
    auto& camera     = scene.cameras.at(options.camera_id);
    auto  image_size = get_camera_image_size(camera, options.image_size);
    rendered_image   = {image_size, zero_vec4f};
    pixels           = make_trace_pixels(image_size, options.random_seed);
    auto regions     = make_image_regions(rendered_image.size());
    if (options.cancel_flag) *options.cancel_flag = false;

#if 0
    auto nthreads = thread::hardware_concurrency();
    threads.clear();
    for (auto tid = 0; tid < nthreads; tid++) {
        threads.push_back(thread([options, tid, nthreads, regions, &current_sample, &queue, &lights, &scene, &rendered_image, &pixels, &bvh]() {
            for (auto s = 0; s < options.num_samples;
                 s += options.samples_per_batch) {
                if (!tid) current_sample = s;
                auto num_samples = min(options.samples_per_batch,
                    options.num_samples - current_sample);
                for (auto region_id = tid; region_id < regions.size();
                     region_id += nthreads) {
                    if (options.cancel_flag && *options.cancel_flag) break;
                    auto& region = regions[region_id];
                    trace_image_region(rendered_image, pixels, scene, bvh,
                        lights, region, num_samples, options);
                    queue.push(region);
                }
            }
            if (!tid) current_sample = options.num_samples;
        }));
    }
#else
    threads.clear();
    threads.emplace_back([options, regions, &current_sample, &rendered_image,
                             &scene, &lights, &bvh, &pixels, &queue]() {
        for (auto sample = 0; sample < options.num_samples;
             sample += options.samples_per_batch) {
            if (options.cancel_flag && *options.cancel_flag) return;
            current_sample   = sample;
            auto num_samples = min(options.samples_per_batch,
                options.num_samples - current_sample);
            parallel_foreach(regions,
                [num_samples, &options, &rendered_image, &scene, &lights, &bvh,
                    &pixels, &queue](const bbox2i& region) {
                    trace_image_region(rendered_image, pixels, scene, bvh,
                        lights, region, num_samples, options);
                    queue.push(region);
                },
                options.cancel_flag);
        }
        current_sample = options.num_samples;
    });
#endif
}

// Stop the asynchronous renderer.
void trace_image_async_stop(vector<thread>& threads,
    concurrent_queue<bbox2i>& queue, const trace_image_options& options) {
    if (options.cancel_flag) *options.cancel_flag = true;
    for (auto& t : threads) t.join();
    threads.clear();
    queue.clear();
}

// Trace statistics for last run used for fine tuning implementation.
// For now returns number of paths and number of rays.
pair<uint64_t, uint64_t> get_trace_stats() {
    return {_trace_nrays, _trace_npaths};
}
void reset_trace_stats() {
    _trace_nrays  = 0;
    _trace_npaths = 0;
}

}  // namespace yocto
// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace yocto {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float exponent) {
    return sqrtf(2 / (exponent + 2));
}

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es  = {(1 + sqrt(ks[0])) / (1 - sqrt(ks[0])),
        (1 + sqrt(ks[1])) / (1 - sqrt(ks[1])),
        (1 + sqrt(ks[2])) / (1 - sqrt(ks[2]))};
    esk = {0, 0, 0};
}

// Specular to  eta.
float specular_to_eta(const vec3f& ks) {
    auto f0 = (ks[0] + ks[1] + ks[2]) / 3;
    return (1 + sqrt(f0)) / (1 - sqrt(f0));
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_dielectric(float cosw, const vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta  = vec3f{1, 1, 1} / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = vec3f{1, 1, 1} - vec3f{sin2, sin2, sin2} / eta2;
    if (cos2t[0] < 0 || cos2t[1] < 0 || cos2t[2] < 0)
        return vec3f{1, 1, 1};  // tir

    auto t0 = vec3f{sqrt(cos2t[0]), sqrt(cos2t[1]), sqrt(cos2t[2])};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (vec3f{cosw, cosw, cosw} - t1) / (vec3f{cosw, cosw, cosw} + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2.0f;
}

// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak) {
    if (etak == zero_vec3f) return fresnel_dielectric(cosw, eta);

    cosw       = clamp(cosw, (float)-1, (float)1);
    auto cos2  = cosw * cosw;
    auto sin2  = clamp(1 - cos2, (float)0, (float)1);
    auto eta2  = eta * eta;
    auto etak2 = etak * etak;

    auto t0         = eta2 - etak2 - vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2   = vec3f{
        sqrt(a2plusb2_2[0]), sqrt(a2plusb2_2[1]), sqrt(a2plusb2_2[2])};
    auto t1  = a2plusb2 + vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a   = vec3f{sqrt(a_2[0]), sqrt(a_2[1]), sqrt(a_2[2])};
    auto t2  = 2.0f * a * cosw;
    auto rs  = (t1 - t2) / (t1 + t2);

    auto t3 = vec3f{cos2, cos2, cos2} * a2plusb2 +
              vec3f{sin2, sin2, sin2} * vec3f{sin2, sin2, sin2};
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2.0f;
}

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& ks, float cosw) {
    if (ks == zero_vec3f) return zero_vec3f;
    return ks + (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    if (ks == zero_vec3f) return zero_vec3f;
    auto fks = fresnel_schlick(ks, cosw);
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX pdf
float sample_ggx_pdf(float rs, float ndh) {
    auto alpha2 = rs * rs;
    auto di     = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d      = alpha2 / (pif * di * di);
    return d * ndh;
}

// Sample the GGX distribution
vec3f sample_ggx(float rs, const vec2f& rn) {
    auto tan2 = rs * rs * rn[1] / (1 - rn[1]);
    auto rz   = sqrt(1 / (tan2 + 1));
    auto rr   = sqrt(1 - rz * rz);
    auto rphi = 2 * pif * rn[0];
    // set to wh
    auto wh_local = vec3f{rr * cos(rphi), rr * sin(rphi), rz};
    return wh_local;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// NUMERICAL TESTS FOR MONTE CARLO INTEGRATION
// -----------------------------------------------------------------------------
namespace yocto {

float integrate_func_base(
    function<float(float)> f, float a, float b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = get_random_float(rng);
        auto x = a + r * (b - a);
        integral += f(x) * (b - a);
    }
    integral /= nsamples;
    return integral;
}

float integrate_func_stratified(
    function<float(float)> f, float a, float b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = (i + get_random_float(rng)) / nsamples;
        auto x = a + r * (b - a);
        integral += f(x) * (b - a);
    }
    integral /= nsamples;
    return integral;
}

float integrate_func_importance(function<float(float)> f,
    function<float(float)> pdf, function<float(float)> warp, int nsamples,
    rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = get_random_float(rng);
        auto x = warp(r);
        integral += f(x) / pdf(x);
    }
    integral /= nsamples;
    return integral;
}

// compute the integral and error using different Monte Carlo scehems
// example 1: ---------
// auto f = [](double x) { return 1.0 - (3.0 / 4.0) * x * x; };
// auto a = 0.0, b = 1.0;
// auto expected = 3.0 / 4.0;
// auto nsamples = 10000
// example 2: ---------
// auto f = [](double x) { return sin(x); }
// auto a = 0.0, b = (double)M_PI;
// auto expected = (double)M_PI;
void print_integrate_func_test(function<float(float)> f, float a, float b,
    float expected, int nsamples, function<float(float)> pdf,
    function<float(float)> warp) {
    auto rng = rng_state();
    printf("nsamples base base-err stratified-err importance-err\n");
    for (auto ns = 10; ns < nsamples; ns += 10) {
        auto integral_base       = integrate_func_base(f, a, b, ns, rng);
        auto integral_stratified = integrate_func_stratified(f, a, b, ns, rng);
        auto integral_importance = integrate_func_importance(
            f, pdf, warp, ns, rng);
        auto error_base       = fabs(integral_base - expected) / expected;
        auto error_stratified = fabs(integral_stratified - expected) / expected;
        auto error_importance = fabs(integral_importance - expected) / expected;
        printf("%d %g %g %g %g\n", ns, integral_base, error_base,
            error_stratified, error_importance);
    }
}

float integrate_func2_base(
    function<float(vec2f)> f, vec2f a, vec2f b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = get_random_vec2f(rng);
        auto x = a + r * (b - a);
        integral += f(x) * (b[0] - a[0]) * (b[1] - a[1]);
    }
    integral /= nsamples;
    return integral;
}

float integrate_func2_stratified(
    function<float(vec2f)> f, vec2f a, vec2f b, int nsamples, rng_state& rng) {
    auto integral  = 0.0f;
    auto nsamples2 = (int)sqrt(nsamples);
    for (auto i = 0; i < nsamples2; i++) {
        for (auto j = 0; j < nsamples2; j++) {
            auto r = vec2f{(i + get_random_float(rng)) / nsamples2,
                (j + get_random_float(rng)) / nsamples2};
            auto x = a + r * (b - a);
            integral += f(x) * (b[0] - a[0]) * (b[1] - a[1]);
        }
    }
    integral /= nsamples2 * nsamples2;
    return integral;
}

float integrate_func2_importance(function<float(vec2f)> f,
    function<float(vec2f)> pdf, function<vec2f(vec2f)> warp, int nsamples,
    rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = get_random_vec2f(rng);
        auto x = warp(r);
        integral += f(x) / pdf(x);
    }
    integral /= nsamples;
    return integral;
}

// compute the integral and error using different Monte Carlo scehems
// example 1: ---------
// auto f = [](double x) { return 1.0 - (3.0 / 4.0) * x * x; };
// auto a = 0.0, b = 1.0;
// auto expected = 3.0 / 4.0;
// auto nsamples = 10000
void print_integrate_func2_test(function<float(vec2f)> f, vec2f a, vec2f b,
    float expected, int nsamples, function<float(vec2f)> pdf,
    function<vec2f(vec2f)> warp) {
    auto rng = rng_state();
    printf("nsamples base base-err stratified-err importance-err\n");
    for (auto ns = 10; ns < nsamples; ns += 10) {
        auto integral_base       = integrate_func2_base(f, a, b, ns, rng);
        auto integral_stratified = integrate_func2_stratified(f, a, b, ns, rng);
        auto integral_importance = integrate_func2_importance(
            f, pdf, warp, ns, rng);
        auto error_base       = fabs(integral_base - expected) / expected;
        auto error_stratified = fabs(integral_stratified - expected) / expected;
        auto error_importance = fabs(integral_importance - expected) / expected;
        printf("%d %g %g %g %g\n", ns, integral_base, error_base,
            error_stratified, error_importance);
    }
}

}  // namespace yocto
