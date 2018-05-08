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
#include "yocto_bvh.h"
#include "yocto_image.h"
#include "yocto_utils.h"

#include <cassert>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a scene handling opacity.
scene_intersection intersect_ray_cutout(
    const scene* scn, const ray3f& ray_, rng_state& rng, int nbounces) {
    auto ray = ray_;
    for (auto b = 0; b < nbounces; b++) {
        auto isec = intersect_ray(scn, ray);
        if (!isec.ist) return isec;
        auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
        if (op > 0.999f) return isec;
        if (rand1f(rng) < op) return isec;
        ray = make_ray(eval_pos(isec.ist, isec.ei, isec.uv), ray.d);
    }
    return {};
}

// Evaluates emission.
vec3f eval_emission(const vec3f& ke, const vec3f& n, const vec3f& o) {
    if (dot(n, o) >= 0) return ke;
    return zero3f;
}

// Check if we are near the mirror direction.
inline bool check_near_mirror(const vec3f& n, const vec3f& o, const vec3f& i) {
    return fabs(dot(i, normalize(n * 2.0f * dot(o, n) - o)) - 1) < 0.001f;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_brdfcos(
    const brdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    auto brdfcos = zero3f;

    auto ndo = dot(n, o), ndi = dot(n, i);

    if (f.kd != zero3f && ndi > 0 && ndo > 0) { brdfcos += f.kd * ndi / pi; }

    if (f.ks != zero3f && f.rs && ndi > 0 && ndo > 0) {
        auto h = normalize(o + i);
        auto ndh = clamp(dot(h, n), -1.0f, 1.0f);
        auto dg = eval_ggx(f.rs, ndh, ndi, ndo);
        auto odh = clamp(dot(o, h), 0.0f, 1.0f);
        auto ks = fresnel_schlick(f.ks, odh, f.rs);
        brdfcos += ks * ndi * dg / (4 * ndi * ndo);
    }

    if (f.kt != zero3f && f.rs && ndo > 0 && ndi < 0) {
        auto ir = i - 2 * dot(i, n) * n;
        auto h = normalize(o + ir);
        auto ndh = clamp(dot(h, n), -1.0f, 1.0f);
        auto dg = eval_ggx(f.rs, ndh, -ndi, ndo);
        auto odh = clamp(dot(o, h), 0.0f, 1.0f);
        auto kt = f.kt * (vec3f{1, 1, 1} - fresnel_schlick(f.ks, odh, f.rs));
        brdfcos += kt * ndi * dg / (4 * ndi * ndo);
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_delta_brdfcos(
    const brdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    auto brdfcos = zero3f;

    auto ndo = dot(n, o), ndi = dot(n, i);

    if (f.ks != zero3f && !f.rs && ndo > 0 && check_near_mirror(n, o, i)) {
        auto ks = fresnel_schlick(f.ks, ndo);
        brdfcos += ks;
    }

    if (f.kt != zero3f && !f.rs && o == -i) {
        auto kt = f.kt * (vec3f{1, 1, 1} - fresnel_schlick(f.ks, ndo));
        brdfcos += kt;
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Compute the weight for sampling the BRDF
float sample_brdf_pdf(
    const brdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    auto prob_kd = max(f.kd), prob_ks = max(f.ks), prob_kt = max(f.kt);
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto ndo = dot(n, o), ndi = dot(n, i);

    auto pdf = 0.0f;

    if (prob_kd && ndo >= 0 && ndi >= 0) { pdf += prob_kd * ndi / pi; }

    if (prob_ks && f.rs && ndo >= 0 && ndi >= 0) {
        auto h = normalize(i + o);
        auto ndh = dot(n, h);
        auto d = sample_ggx_pdf(f.rs, ndh);
        auto hdo = dot(o, h);
        pdf += prob_ks * d / (4 * hdo);
    }

    if (prob_kt && f.rs && ndo >= 0 && ndi <= 0) {
        auto wir = i - 2 * dot(i, n) * n;
        auto h = normalize(o + wir);
        auto ndh = dot(n, h);
        auto d = sample_ggx_pdf(f.rs, ndh);
        auto hdo = dot(o, h);
        pdf += prob_kt * d / (4 * hdo);
    }

    return pdf;
}

// Compute the weight for sampling the BRDF
float sample_delta_brdf_pdf(
    const brdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    if (f.rs) return 0;
    auto prob_ks = max(f.ks), prob_kt = max(f.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto ndo = dot(n, o), ndi = dot(n, i);

    auto pdf = 0.0f;

    if (prob_ks && !f.rs && ndo > 0 && check_near_mirror(n, o, i)) {
        pdf += prob_ks;
    }

    if (prob_kt && !f.rs && o == -i) { pdf += prob_kt; }

    return pdf;
}

// Picks a direction based on the BRDF
vec3f sample_brdf(
    const brdf& f, const vec3f& n, const vec3f& o, float rnl, const vec2f& rn) {
    auto prob_kd = max(f.kd), prob_ks = (f.rs) ? max(f.ks) : 0,
         prob_kt = (f.rs) ? max(f.kt) : 0;
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto ndo = dot(n, o);
    if (ndo <= 0) return zero3f;

    // sample according to diffuse
    if (rnl < prob_kd) {
        auto fp = make_frame_fromz(zero3f, n);
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
        auto il = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return transform_direction(fp, il);
    }
    // sample according to specular GGX
    else if (rnl < prob_kd + prob_ks && f.rs) {
        auto fp = make_frame_fromz(zero3f, n);
        auto hl = sample_ggx(f.rs, rn);
        auto h = transform_direction(fp, hl);
        return normalize(h * 2.0f * dot(o, h) - o);
    }
    // transmission hack
    else if (rnl < prob_kd + prob_ks + prob_kt && f.rs) {
        auto fp = make_frame_fromz(zero3f, n);
        auto hl = sample_ggx(f.rs, rn);
        auto h = transform_direction(fp, hl);
        auto i = normalize(h * 2.0f * dot(o, h) - o);
        return normalize(i - 2 * dot(i, n) * n);
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_delta_brdf(
    const brdf& f, const vec3f& n, const vec3f& o, float rnl, const vec2f& rn) {
    if (f.rs) return zero3f;
    auto prob_ks = max(f.ks), prob_kt = max(f.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto ndo = dot(n, o);
    if (ndo <= 0) return zero3f;

    // sample according to specular mirror
    if (rnl < prob_ks && !f.rs) {
        return normalize(n * 2.0f * dot(o, n) - o);
    }
    // transmission hack
    else if (rnl < prob_ks + prob_kt && !f.rs) {
        return -o;
    } else {
        return zero3f;
    }
}

float sample_delta_prob(const brdf& f, const vec3f& n, const vec3f& o) {
    if (f.rs) return 0;
    auto dw = max(f.ks) + max(f.kt);
    return dw / (dw + max(f.kd));
}

// Sample pdf for an environment.
float sample_environment_pdf(const environment* env, const vec3f& i) {
    auto txt = env->ke_txt.txt;
    if (!env->elem_cdf.empty() && txt) {
        auto texcoord = eval_texcoord(env, i);
        auto i = (int)(texcoord.x * txt->width);
        auto j = (int)(texcoord.y * txt->height);
        auto idx = j * txt->width + i;
        auto prob =
            sample_discrete_pdf(env->elem_cdf, idx) / env->elem_cdf.back();
        auto angle = (2 * pi / txt->width) * (pi / txt->height) *
                     sin(pi * (j + 0.5f) / txt->height);
        return prob / angle;
    } else {
        return 1 / (4 * pi);
    }
}

// Picks a point on an environment.
vec3f sample_environment(const environment* env, float rel, const vec2f& ruv) {
    auto txt = env->ke_txt.txt;
    if (!env->elem_cdf.empty() && txt) {
        auto idx = sample_discrete(env->elem_cdf, rel);
        auto u = (idx % txt->width + 0.5f) / txt->width;
        auto v = (idx / txt->width + 0.5f) / txt->height;
        return eval_direction(env, {u, v});
    } else {
        return sample_sphere(ruv);
    }
}

// Sample pdf for a light point.
float sample_light_pdf(const instance* ist, int ei, const vec3f& i,
    const vec3f& lnorm, float dist) {
    if (!ist) return 0;
    auto shp = ist->shp;
    if (shp->elem_cdf.empty()) return 0;
    if (!shp->triangles.empty()) {
        // prob triangle * area triangle = area triangle mesh
        auto area = shp->elem_cdf.back();
        return (dist * dist) / (fabs(dot(lnorm, i)) * area);
    } else if (!shp->lines.empty()) {
        throw std::runtime_error("not implemented yet");
        return 0;
    } else if (!shp->pos.empty()) {
        auto num = shp->elem_cdf.back();
        return (dist * dist) / num;
    } else {
        return 0;
    }
}

// Picks a point on a light.
vec3f sample_light(
    const instance* ist, const vec3f& pos, float rel, const vec2f& ruv) {
    auto sample = sample_shape(ist->shp, rel, ruv);
    return normalize(eval_pos(ist, sample.first, sample.second) - pos);
}

// Test occlusion.
vec3f eval_transmission(
    const scene* scn, const vec3f& from, const vec3f& to, int nbounces) {
    auto weight = vec3f{1, 1, 1};
    auto pos = from;
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        auto ray = make_segment(pos, to);
        auto isec = intersect_ray(scn, ray);
        if (!isec.ist) break;
        auto f = eval_brdf(isec.ist, isec.ei, isec.uv);
        auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
        weight *= f.kt + vec3f{1 - op, 1 - op, 1 - op};
        if (weight == zero3f) break;
        pos = eval_pos(isec.ist, isec.ei, isec.uv);
    }
    return weight;
}

// Recursive path tracing.
vec3f trace_path(
    const scene* scn, const ray3f& ray_, rng_state& rng, int nbounces) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto emission = true;
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            if (emission) {
                for (auto env : scn->environments)
                    l += weight * eval_environment(env, ray.d);
            }
            break;
        }

        // point
        auto o = -ray.d;
        auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

        // emission
        if (emission) {
            auto em = eval_emission(isec.ist, isec.ei, isec.uv);
            l += weight * eval_emission(em, n, o);
        }

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // choose delta
        auto delta_prob = sample_delta_prob(f, n, o);
        auto delta = (delta_prob && rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);

        // direct
        if (!delta && (!scn->lights.empty() || !scn->environments.empty())) {
            auto i = zero3f;
            auto nlights = (int)(scn->lights.size() + scn->environments.size());
            if (rand1f(rng) < 0.5f) {
                auto idx = sample_index(nlights, rand1f(rng));
                if (idx < scn->lights.size()) {
                    auto lgt = scn->lights[idx];
                    i = sample_light(lgt, pos, rand1f(rng), rand2f(rng));
                } else {
                    auto env = scn->environments[idx - scn->lights.size()];
                    i = sample_environment(env, rand1f(rng), rand2f(rng));
                }
            } else {
                i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            }
            auto isec =
                intersect_ray_cutout(scn, make_ray(pos, i), rng, nbounces);
            auto pdf = 0.5f * sample_brdf_pdf(f, n, o, i);
            auto le = zero3f;
            if (isec.ist) {
                auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
                auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
                pdf += 0.5f *
                       sample_light_pdf(
                           isec.ist, isec.ei, i, ln, length(lp - pos)) *
                       sample_index_pdf(nlights);
                auto em = eval_emission(isec.ist, isec.ei, isec.uv);
                le += eval_emission(em, ln, -i);
            } else {
                for (auto env : scn->environments) {
                    pdf += 0.5f * sample_environment_pdf(env, i) *
                           sample_index_pdf(nlights);
                    le += eval_environment(env, i);
                }
            }
            auto brdfcos = eval_brdfcos(f, n, o, i);
            if (pdf != 0) l += weight * le * brdfcos / pdf;
        }

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!delta) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_brdfcos(f, n, o, i);
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdfcos(f, n, o, i);
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(pos, i);
        emission = delta;
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_naive(
    const scene* scn, const ray3f& ray_, rng_state& rng, int nbounces) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            for (auto env : scn->environments)
                l += weight * eval_environment(env, ray.d);
            break;
        }

        // point
        auto o = -ray.d;
        auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

        // emission
        auto em = eval_emission(isec.ist, isec.ei, isec.uv);
        l += weight * eval_emission(em, n, o);

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // choose delta
        auto delta_prob = sample_delta_prob(f, n, o);
        auto delta = (delta_prob && rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!delta) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_brdfcos(f, n, o, i);
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdfcos(f, n, o, i);
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(pos, i);
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_nomis(
    const scene* scn, const ray3f& ray_, rng_state& rng, int nbounces) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto emission = true;
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            for (auto env : scn->environments)
                l += weight * eval_environment(env, ray.d);
            break;
        }

        // point
        auto o = -ray.d;
        auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

        // emission
        if (emission) {
            auto em = eval_emission(isec.ist, isec.ei, isec.uv);
            l += weight * eval_emission(em, n, o);
        }

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // choose delta
        auto delta_prob = sample_delta_prob(f, n, o);
        auto delta = (delta_prob && rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);

        // direct
        if (!delta && !scn->lights.empty()) {
            auto lgt =
                scn->lights[sample_index(scn->lights.size(), rand1f(rng))];
            auto i = sample_light(lgt, pos, rand1f(rng), rand2f(rng));
            auto isec =
                intersect_ray_cutout(scn, make_ray(pos, i), rng, nbounces);
            if (isec.ist && isec.ist->mat->ke != zero3f) {
                auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
                auto pdf =
                    sample_light_pdf(isec.ist, isec.ei, i, ln, isec.dist) *
                    sample_index_pdf(scn->lights.size());
                auto em = eval_emission(isec.ist, isec.ei, isec.uv);
                auto le = eval_emission(em, ln, -i);
                auto brdfcos = eval_brdfcos(f, n, o, i);
                if (pdf != 0) l += weight * le * brdfcos / pdf;
            }
        }

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!delta) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_brdfcos(f, n, o, i);
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdfcos(f, n, o, i);
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(pos, i);
        emission = delta;
    }

    return l;
}

// Direct illumination.
vec3f trace_direct(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }

    // point
    auto o = -ray.d;
    auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

    // emission
    auto em = eval_emission(isec.ist, isec.ei, isec.uv);
    l += eval_emission(em, n, o);

    // direct lights
    for (auto lgt : scn->lights) {
        auto i = zero3f;
        if (rand1f(rng) < 0.5f) {
            i = sample_light(lgt, pos, rand1f(rng), rand2f(rng));
        } else {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
        }
        auto isec = intersect_ray(scn, make_ray(pos, i));
        if (lgt != isec.ist) continue;
        auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
        auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
        auto pdf = 0.5f * sample_light_pdf(
                              isec.ist, isec.ei, i, ln, length(lp - pos)) +
                   0.5f * sample_brdf_pdf(f, n, o, i);
        auto em = eval_emission(isec.ist, isec.ei, isec.uv);
        auto le = eval_emission(em, ln, -i);
        auto brdfcos = eval_brdfcos(f, n, o, i);
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // direct environments
    for (auto env : scn->environments) {
        auto i = zero3f;
        if (rand1f(rng) < 0.5f) {
            i = sample_environment(env, rand1f(rng), rand2f(rng));
        } else {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
        }
        auto isec = intersect_ray(scn, make_ray(pos, i));
        if (isec.ist) continue;
        auto pdf = 0.5f * sample_environment_pdf(env, i) +
                   0.5f * sample_brdf_pdf(f, n, o, i);
        auto le = eval_environment(env, i);
        auto brdfcos = eval_brdfcos(f, n, o, i);
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // exit if needed
    if (nbounces <= 0) return l;

    // reflection
    if (f.ks != zero3f && !f.rs) {
        auto i = reflect(o, n);
        l += f.ks * trace_direct(scn, make_ray(pos, i), rng, nbounces - 1);
    }

    // refraction
    if (f.kt != zero3f) {
        l += f.kt * trace_direct(scn, make_ray(pos, -o), rng, nbounces - 1);
    }

    // opacity
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l +
            (1 - op) * trace_direct(scn, make_ray(pos, -o), rng, nbounces - 1);
    }

    // done
    return l;
}

// Direct illumination.
vec3f trace_direct_nomis(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }

    // point
    auto o = -ray.d;
    auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

    // emission
    auto em = eval_emission(isec.ist, isec.ei, isec.uv);
    l += eval_emission(em, n, o);

    // direct lights
    for (auto lgt : scn->lights) {
        auto i = sample_light(lgt, pos, rand1f(rng), rand2f(rng));
        auto isec = intersect_ray(scn, make_ray(pos, i));
        if (lgt != isec.ist) continue;
        auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
        auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
        auto pdf = sample_light_pdf(isec.ist, isec.ei, i, ln, length(lp - pos));
        auto em = eval_emission(isec.ist, isec.ei, isec.uv);
        auto le = eval_emission(em, ln, -i);
        auto brdfcos = eval_brdfcos(f, n, o, i);
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // direct environments
    for (auto env : scn->environments) {
        auto i = sample_environment(env, rand1f(rng), rand2f(rng));
        auto isec = intersect_ray(scn, make_ray(pos, i));
        if (isec.ist) continue;
        auto pdf = sample_environment_pdf(env, i);
        auto le = eval_environment(env, i);
        auto brdfcos = eval_brdfcos(f, n, o, i);
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // exit if needed
    if (nbounces <= 0) return l;

    // reflection
    if (f.ks != zero3f && !f.rs) {
        auto i = reflect(o, n);
        l += f.ks * trace_direct(scn, make_ray(pos, i), rng, nbounces - 1);
    }

    // opacity
    if (f.kt != zero3f) {
        l += f.kt * trace_direct(scn, make_ray(pos, -o), rng, nbounces - 1);
    }

    // opacity
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l +
            (1 - op) * trace_direct(scn, make_ray(pos, -o), rng, nbounces - 1);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
vec3f trace_eyelight(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }

    // point
    auto o = -ray.d;
    auto pos = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

    // emission
    auto em = eval_emission(isec.ist, isec.ei, isec.uv);
    l += eval_emission(em, n, o);

    // brdf*light
    l += eval_brdfcos(f, n, o, o) * pi;

    // opacity
    if (nbounces <= 0) return l;
    if (f.kt != zero3f) {
        l += f.kt * trace_eyelight(scn, make_ray(pos, -o), rng, nbounces - 1);
    }
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l + (1 - op) * trace_eyelight(
                                    scn, make_ray(pos, -o), rng, nbounces - 1);
    }

    // done
    return l;
}

// Debug previewing.
vec3f trace_debug_normal(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;

    // point
    auto o = -ray.d;
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);

    // shade
    return n * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug frontfacing.
vec3f trace_debug_frontfacing(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;

    // point
    auto o = -ray.d;
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);

    // shade
    return dot(n, o) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
}

// Debug previewing.
vec3f trace_debug_albedo(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;

    // point
    auto f = eval_brdf(isec.ist, isec.ei, isec.uv);

    // shade
    return f.kd + f.ks + f.kt;
}

// Debug previewing.
vec3f trace_debug_texcoord(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;

    // point
    auto texcoord = eval_texcoord(isec.ist, isec.ei, isec.uv);

    // shade
    return {texcoord.x, texcoord.y, 0};
}

// Trace a single sample
vec4f trace_sample(const scene* scn, const camera* cam, int i, int j, int width,
    int height, rng_state& rng, trace_func tracer, int nbounces,
    float pixel_clamp = 100) {
    auto crn = rand2f(rng);
    auto lrn = rand2f(rng);
    auto uv = vec2f{(i + crn.x) / width, 1 - (j + crn.y) / height};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto l = tracer(scn, ray, rng, nbounces);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        l = zero3f;
    }
    if (max(l) > pixel_clamp) l = l * (pixel_clamp / max(l));
    // return {l.x, l.y, l.z, (pt.ist || pt.env) ? 1.0f : 0.0f};
    return {l.x, l.y, l.z, 1};
}

// Trace the next nsamples.
void trace_samples(const scene* scn, const camera* cam, int width, int height,
    std::vector<vec4f>& img, std::vector<rng_state>& rngs, int sample,
    int nsamples, trace_func tracer, int nbounces, float pixel_clamp) {
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto pid = i + j * width;
            img[pid] *= sample;
            for (auto s = 0; s < nsamples; s++)
                img[pid] += trace_sample(scn, cam, i, j, width, height,
                    rngs[pid], tracer, nbounces, pixel_clamp);
            img[pid] /= sample + nsamples;
        }
    }
}

// Trace the next nsamples.
void trace_samples_mt(const scene* scn, const camera* cam, int width,
    int height, std::vector<vec4f>& img, std::vector<rng_state>& rngs,
    int sample, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp) {
    auto nthreads = std::thread::hardware_concurrency();
    auto threads = std::vector<std::thread>();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &rngs]() {
            for (auto j = tid; j < height; j += nthreads) {
                for (auto i = 0; i < width; i++) {
                    auto pid = i + j * width;
                    img[pid] *= sample;
                    for (auto s = 0; s < nsamples; s++)
                        img[pid] += trace_sample(scn, cam, i, j, width, height,
                            rngs[pid], tracer, nbounces, pixel_clamp);
                    img[pid] /= sample + nsamples;
                }
            }
        }));
    }
    for (auto& t : threads) t.join();
    threads.clear();
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, int width,
    int height, std::vector<vec4f>& img, std::vector<rng_state>& rngs,
    int nsamples, trace_func tracer, int nbounces,
    std::vector<std::thread>& threads, bool& stop_flag, int& cur_sample,
    float pixel_clamp) {
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < nthreads; tid++) {
        threads.push_back(
            std::thread([=, &img, &rngs, &stop_flag, &cur_sample]() {
                for (auto sample = 0; sample < nsamples; sample++) {
                    if (!tid) cur_sample = sample;
                    for (auto j = tid; j < height; j += nthreads) {
                        for (auto i = 0; i < width; i++) {
                            auto pid = i + j * width;
                            if (stop_flag) return;
                            auto l = trace_sample(scn, cam, i, j, width, height,
                                rngs[pid], tracer, nbounces);
                            img[pid] = (img[pid] * sample + l) / (sample + 1);
                        }
                    }
                    if (!tid) cur_sample = nsamples;
                }
            }));
    }
}

// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag) {
    stop_flag = true;
    for (auto& t : threads) t.join();
    stop_flag = false;
    threads.clear();
}

}  // namespace ygl
// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness. Public API, see above.
float specular_exponent_to_roughness(float n) { return sqrtf(2 / (n + 2)); }

// Specular to fresnel eta. Public API, see above.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es = {(1 + sqrt(ks.x)) / (1 - sqrt(ks.x)),
        (1 + sqrt(ks.y)) / (1 - sqrt(ks.y)),
        (1 + sqrt(ks.z)) / (1 - sqrt(ks.z))};
    esk = {0, 0, 0};
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_dielectric(float cosw, const vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1.0f / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = vec3f{1, 1, 1} - sin2 / eta2;
    if (cos2t.x < 0 || cos2t.y < 0 || cos2t.z < 0)
        return vec3f{1, 1, 1};  // tir

    auto t0 = vec3f{sqrt(cos2t.x), sqrt(cos2t.y), sqrt(cos2t.z)};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (vec3f{cosw, cosw, cosw} - t1) / (vec3f{cosw, cosw, cosw} + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2.0f;
}

// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak) {
    if (etak == zero3f) return fresnel_dielectric(cosw, eta);

    cosw = clamp(cosw, (float)-1, (float)1);
    auto cos2 = cosw * cosw;
    auto sin2 = clamp(1 - cos2, (float)0, (float)1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2 =
        vec3f{sqrt(a2plusb2_2.x), sqrt(a2plusb2_2.y), sqrt(a2plusb2_2.z)};
    auto t1 = a2plusb2 + vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a = vec3f{sqrt(a_2.x), sqrt(a_2.y), sqrt(a_2.z)};
    auto t2 = 2.0f * a * cosw;
    auto rs = (t1 - t2) / (t1 + t2);

    auto t3 = vec3f{cos2, cos2, cos2} * a2plusb2 +
              vec3f{sin2, sin2, sin2} * vec3f{sin2, sin2, sin2};
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2.0f;
}

// Schlick approximation of Fresnel term
vec3f fresnel_schlick(const vec3f& ks, float cosw) {
    if (ks == zero3f) return zero3f;
    return ks +
           (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

// Schlick approximation of Fresnel term weighted by roughness.
// This is a hack, but works better than not doing it.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    if (ks == zero3f) return zero3f;
    auto fks = fresnel_schlick(ks, cosw);
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX distribution and geometric term
float eval_ggx(float rs, float ndh, float ndi, float ndo) {
    // evaluate D
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pi * di * di);
#if 0
    // evaluate G from Heitz
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    // evaluate G from Smith
    auto go = (2 * ndo) / (ndo + sqrt(alpha2 + (1 - alpha2) * ndo * ndo));
    auto gi = (2 * ndi) / (ndi + sqrt(alpha2 + (1 - alpha2) * ndi * ndi));
    auto g = go * gi;
#endif
    return d * g;
}

// Evaluates the GGX pdf
float sample_ggx_pdf(float rs, float ndh) {
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pi * di * di);
    return d * ndh;
}

// Sample the GGX distribution
vec3f sample_ggx(float rs, const vec2f& rn) {
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz = sqrt(1 / (tan2 + 1));
    auto rr = sqrt(1 - rz * rz);
    auto rphi = 2 * pi * rn.x;
    // set to wh
    auto wh_local = vec3f{rr * cos(rphi), rr * sin(rphi), rz};
    return wh_local;
}

}  // namespace ygl
