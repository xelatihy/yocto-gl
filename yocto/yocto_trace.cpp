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
#include "yocto_utils.h"
#include "yocto_image.h"

#include <cassert>

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Generates a 1-dimensional sample.
float sample_next1f(trace_pixel& pxl) {
    return clamp(next_rand1f(pxl.rng), 0.0f, 1 - flt_eps);
}

// Generates a 2-dimensional sample.
vec2f sample_next2f(trace_pixel& pxl) {
    return {next_rand1f(pxl.rng), next_rand1f(pxl.rng)};
}

// Type of point.
enum struct trace_point_type {
    none = 0,         // unitialized point
    point = 1,        // point
    curve = 2,        // curve
    surface = 3,      // surface
    environment = 4,  // environment
};

// Surface point with geometry and material data. Supports point on
// envmap too. This is the key data manipulated in the path tracer.
struct trace_point {
    const instance* ist = nullptr;                   // shape instance
    const environment* env = nullptr;                // environment
    trace_point_type type = trace_point_type::none;  // type
    vec3f pos = zero3f;                              // pos
    vec3f norm = {0, 0, 1};                          // norm
    vec2f texcoord = zero2f;                         // texcoord
    vec3f ke = zero3f;                               // emission
    vec3f kd = {0, 0, 0};                            // diffuse
    vec3f ks = {0, 0, 0};                            // specular
    float rs = 0;                                    // specular roughness
    vec3f kt = {0, 0, 0};                            // thin glass transmission
    float op = 1;                                    // opacity
    bool double_sided = false;                       // double sided
};

// Create a point for an environment map. Resolves material with textures.
trace_point eval_point(const environment* env, const vec2f& uv) {
    auto pt = trace_point();
    pt.env = env;
    pt.type = trace_point_type::environment;

    pt.pos = eval_pos(env, uv);
    pt.norm = eval_norm(env, uv);
    pt.texcoord = eval_texcoord(env, uv);

    pt.pos = transform_point(env->frame, pt.pos);
    pt.norm = transform_direction(env->frame, pt.norm);

    pt.ke = env->ke;
    if (env->ke_txt.txt) {
        auto txt = eval_texture(env->ke_txt, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }
    return pt;
}

// Create a point for a shape. Resolves geometry and material with textures.
trace_point eval_point(
    const instance* ist, int eid, const vec2f& euv, bool force_double_sided) {
    // default material
    static auto def_material = (material*)nullptr;
    if (!def_material) {
        def_material = new material();
        def_material->kd = {0.2f, 0.2f, 0.2f};
        def_material->rs = 1;
    }

    // point
    auto pt = trace_point();

    // shape
    pt.ist = ist;
    switch (get_shape_type(ist->shp)) {
        case shape_elem_type::none: {
            return pt;
        } break;
        case shape_elem_type::points:
        case shape_elem_type::vertices: {
            pt.type = trace_point_type::point;
        } break;
        case shape_elem_type::lines:
        case shape_elem_type::beziers: {
            pt.type = trace_point_type::curve;
        } break;
        case shape_elem_type::triangles:
        case shape_elem_type::quads:
        case shape_elem_type::facevarying: {
            pt.type = trace_point_type::surface;
        } break;
    }

    // shape values
    pt.pos = eval_pos(ist->shp, eid, euv);
    pt.norm = eval_norm(ist->shp, eid, euv);
    pt.texcoord = eval_texcoord(ist->shp, eid, euv);

    // shortcuts
    auto mat = (ist->mat) ? ist->mat : def_material;

    // handle normal map
    if (mat->norm_txt.txt) {
        auto tangsp = eval_tangsp(ist->shp, eid, euv);
        auto txt = eval_texture(mat->norm_txt, pt.texcoord, false) * 2.0f -
                   vec4f{1, 1, 1, 1};
        auto ntxt = normalize(vec3f{txt.x, -txt.y, txt.z});
        auto frame = make_frame_fromzx(
            {0, 0, 0}, pt.norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp.w;
        pt.norm = transform_direction(frame, ntxt);
    }

    // move to world coordinates
    pt.pos = transform_point(ist->frame, pt.pos);
    pt.norm = transform_direction(ist->frame, pt.norm);

    // double-sided
    pt.double_sided = mat->double_sided || force_double_sided;

    // initialized material values
    auto kx = vec3f{1, 1, 1};
    auto op = 1.0f;
    if (!ist->shp->color.empty()) {
        auto col = eval_color(ist->shp, eid, euv);
        kx *= {col.x, col.y, col.z};
        op *= col.w;
    }

    // handle occlusion
    if (mat->occ_txt.txt) {
        auto txt = eval_texture(mat->occ_txt, pt.texcoord);
        kx *= {txt.x, txt.y, txt.z};
    }

    // handle opacity
    pt.op = mat->op * op;
    if (mat->op_txt.txt) {
        auto txt = eval_texture(mat->op_txt, pt.texcoord);
        pt.op *= (txt.x + txt.y + txt.z) / 3;
    }

    // sample emission
    pt.ke = mat->ke * kx;
    if (mat->ke_txt.txt) {
        auto txt = eval_texture(mat->ke_txt, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }

    // sample reflectance
    switch (mat->type) {
        case material_type::specular_roughness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
            }
            pt.kt = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
        case material_type::metallic_roughness: {
            auto kb = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                kb *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            auto km = mat->ks.x;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                km *= txt.y;
                pt.rs *= txt.z;
            }
            pt.kd = kb * (1 - km);
            pt.ks = kb * km + vec3f{0.04f, 0.04f, 0.04f} * (1 - km);
        } break;
        case material_type::specular_glossiness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
                pt.rs *= txt.w;
            }
            pt.rs = 1 - pt.rs;  // glossiness -> roughnes
            pt.kt = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
    }

    // set up final values
    if (pt.ks != zero3f && pt.rs < 0.9999f) {
        pt.rs = pt.rs * pt.rs;
        if (pt.rs < 0.03f * 0.03f)
            pt.rs = 0;
        else
            pt.rs = clamp(pt.rs, 0.03f * 0.03f, 1.0f);
    } else {
        pt.ks = zero3f;
        pt.rs = 1;
    }
    if (pt.kt != zero3f) pt.double_sided = true;
    if (pt.op > 0.999f) pt.op = 1;

    // done
    return pt;
}

// Check if it has reflectance.
bool has_reflectance(const trace_point& pt) {
    return pt.kd + pt.ks + pt.kt != zero3f;
}

// Evaluates emission.
vec3f eval_emission(const trace_point& pt, const vec3f& wo) {
    if (pt.type != trace_point_type::surface || pt.double_sided ||
        dot(pt.norm, wo) >= 0)
        return pt.ke;
    return zero3f;
}

// Check if we are near the mirror direction.
inline bool check_near_mirror(
    const vec3f& wn, const vec3f& wo, const vec3f& wi) {
    return fabs(dot(wi, normalize(wn * 2.0f * dot(wo, wn) - wo)) - 1) < 0.001f;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_surface_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    if (pt.kd != zero3f && ndi > 0 && ndo > 0) { brdfcos += pt.kd * ndi / pi; }

    if (pt.ks != zero3f && pt.rs && ndi > 0 && ndo > 0) {
        auto wh = normalize(wo + wi);
        auto ndh = clamp(dot(wh, wn), -1.0f, 1.0f);
        auto dg = eval_ggx(pt.rs, ndh, ndi, ndo);
        auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
        auto ks = fresnel_schlick(pt.ks, odh, pt.rs);
        brdfcos += ks * ndi * dg / (4 * ndi * ndo);
    }

    if (pt.kt != zero3f && pt.rs && ndo > 0 && ndi < 0) {
        auto wir = wi - 2 * dot(wi, wn) * wn;
        auto wh = normalize(wo + wir);
        auto ndh = clamp(dot(wh, wn), -1.0f, 1.0f);
        auto dg = eval_ggx(pt.rs, ndh, -ndi, ndo);
        auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
        auto kt = pt.kt * (vec3f{1, 1, 1} - fresnel_schlick(pt.ks, odh, pt.rs));
        brdfcos += kt * ndi * dg / (4 * ndi * ndo);
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_surface_delta_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    if (pt.ks != zero3f && !pt.rs && ndo > 0 && check_near_mirror(wn, wo, wi)) {
        auto ks = fresnel_schlick(pt.ks, ndo);
        brdfcos += ks;
    }

    if (pt.kt != zero3f && !pt.rs && wo == -wi) {
        auto kt = pt.kt * (vec3f{1, 1, 1} - fresnel_schlick(pt.ks, ndo));
        brdfcos += kt;
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses Kajiya-Kay for hair
vec3f eval_curve_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);
    auto wn = pt.norm;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi),
         ndh = clamp(dot(wn, wh), 0.0f, 1.0f);
    auto so = sqrt(clamp(1 - ndo * ndo, 0.0001f, 1.0f)),
         si = sqrt(clamp(1 - ndi * ndi, 0.0001f, 1.0f)),
         sh = sqrt(clamp(1 - ndh * ndh, 0.0001f, 1.0f));

    if (pt.kd != zero3f) { brdfcos += pt.kd * si / pi; }

    if (pt.ks != zero3f && pt.rs) {
        auto ns = 2 / (pt.rs * pt.rs) - 2;
        auto d = (ns + 2) * pow(sh, ns) / (2 + pi);
        brdfcos += pt.ks * si * d / (4.0f * si * so);
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses a hack for points
vec3f eval_point_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;

    auto ido = dot(wo, wi);
    brdfcos += pt.kd * (2 * ido + 1) / (2 * pi);

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!has_reflectance(pt)) return zero3f;
    switch (pt.type) {
        case trace_point_type::none: return zero3f;
        case trace_point_type::surface:
            return (!delta) ? eval_surface_brdfcos(pt, wo, wi) :
                              eval_surface_delta_brdfcos(pt, wo, wi);
        case trace_point_type::curve: return eval_curve_brdfcos(pt, wo, wi);
        case trace_point_type::point: return eval_point_brdfcos(pt, wo, wi);
        case trace_point_type::environment: return zero3f;
    }
}

// Compute the weight for sampling the BRDF
float weight_surface_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto prob_kd = max(pt.kd), prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    // normal
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    auto pdf = 0.0f;

    if (prob_kd && ndo > 0 && ndi > 0) { pdf += prob_kd * ndi / pi; }

    if (prob_ks && pt.rs && ndo > 0 && ndi > 0) {
        auto wh = normalize(wi + wo);
        auto ndh = dot(wn, wh);
        auto d = sample_ggx_pdf(pt.rs, ndh);
        auto hdo = dot(wo, wh);
        pdf += prob_ks * d / (4 * hdo);
    }

    if (prob_kt && pt.rs && ndo > 0 && ndi < 0) {
        auto wir = wi - 2 * dot(wi, wn) * wn;
        auto wh = normalize(wo + wir);
        auto ndh = dot(wn, wh);
        auto d = sample_ggx_pdf(pt.rs, ndh);
        auto hdo = dot(wo, wh);
        pdf += prob_kt * d / (4 * hdo);
    }

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_surface_delta_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    if (pt.rs) return 0;
    auto prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    // normal
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    auto pdf = 0.0f;

    if (prob_ks && !pt.rs && ndo > 0 && check_near_mirror(wn, wo, wi)) {
        pdf += prob_ks;
    }

    if (prob_kt && !pt.rs && wo == -wi) { pdf += prob_kt; }

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_curve_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto pdf = 1 / (4 * pi);
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_point_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto pdf = 1 / (4 * pi);
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_brdf(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!has_reflectance(pt)) return 0;
    switch (pt.type) {
        case trace_point_type::none: return 0;
        case trace_point_type::surface:
            return (!delta) ? weight_surface_brdf(pt, wo, wi) :
                              weight_surface_delta_brdf(pt, wo, wi);
        case trace_point_type::curve: return weight_curve_brdf(pt, wo, wi);
        case trace_point_type::point: return weight_point_brdf(pt, wo, wi);
        case trace_point_type::environment: return 0;
    }
}

// Picks a direction based on the BRDF
vec3f sample_surface_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto prob_kd = max(pt.kd), prob_ks = (pt.rs) ? max(pt.ks) : 0,
         prob_kt = (pt.rs) ? max(pt.kt) : 0;
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo);
    if (ndo <= 0) return zero3f;

    // sample according to diffuse
    if (rnl < prob_kd) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return transform_direction(fp, wi_local);
    }
    // sample according to specular GGX
    else if (rnl < prob_kd + prob_ks && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        return normalize(wh * 2.0f * dot(wo, wh) - wo);
    }
    // transmission hack
    else if (rnl < prob_kd + prob_ks + prob_kt && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        auto wi = normalize(wh * 2.0f * dot(wo, wh) - wo);
        return normalize(wi - 2 * dot(wi, wn) * wn);
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_surface_delta_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    if (pt.rs) return zero3f;
    auto prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo);
    if (ndo <= 0) return zero3f;

    // sample according to specular mirror
    if (rnl < prob_ks && !pt.rs) {
        return normalize(wn * 2.0f * dot(wo, wn) - wo);
    }
    // transmission hack
    else if (rnl < prob_ks + prob_kt && !pt.rs) {
        return -wo;
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_curve_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto wn = pt.norm;
    auto fp = make_frame_fromz(pt.pos, wn);
    auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
    auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
    return transform_direction(fp, wi_local);
}

// Picks a direction based on the BRDF
vec3f sample_point_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto wn = pt.norm;
    auto fp = make_frame_fromz(pt.pos, wn);
    auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
    auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
    return transform_direction(fp, wi_local);
}

// Picks a direction based on the BRDF
vec3f sample_brdf(const trace_point& pt, const vec3f& wo, float rnl,
    const vec2f& rn, bool delta = false) {
    if (!has_reflectance(pt)) return zero3f;
    switch (pt.type) {
        case trace_point_type::none: return zero3f;
        case trace_point_type::surface:
            return (!delta) ? sample_surface_brdf(pt, wo, rnl, rn) :
                              sample_surface_delta_brdf(pt, wo, rnl, rn);
        case trace_point_type::curve: return sample_curve_brdf(pt, wo, rnl, rn);
        case trace_point_type::point: return sample_point_brdf(pt, wo, rnl, rn);
        case trace_point_type::environment: return zero3f;
    }
}

float sample_delta_prob(const trace_point& pt, const vec3f& wo) {
    if (pt.type != trace_point_type::surface) return 0;
    if (pt.rs) return 0;
    auto dw = max(pt.ks) + max(pt.kt);
    return dw / (dw + max(pt.kd));
}

// Sample weight for a light point.
float weight_light(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    switch (lpt.type) {
        case trace_point_type::none: {
            throw std::runtime_error("should not have gotten here");
        } break;
        case trace_point_type::point: {
            auto area = lights.shape_distribs.at(lpt.ist->shp).back();
            auto dist = length(lpt.pos - pt.pos);
            return area / (dist * dist);
        } break;
        case trace_point_type::curve: {
            throw std::runtime_error("not implemented yet");
        } break;
        case trace_point_type::surface: {
            auto area = lights.shape_distribs.at(lpt.ist->shp).back();
            auto dist = length(lpt.pos - pt.pos);
            return area * fabs(dot(lpt.norm, normalize(lpt.pos - pt.pos))) /
                   (dist * dist);
        } break;
        case trace_point_type::environment: {
            return 4 * pi;
        } break;
    }
    return 0;
}

// Sample weight for a light point.
float weight_lights(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    return lights.lights.size() * weight_light(lights, lpt, pt);
}

// Picks a point on a light.
trace_point sample_light(const trace_lights& lights, const trace_light& lgt,
    const trace_point& pt, float rel, const vec2f& ruv,
    const trace_params& params) {
    if (lgt.ist) {
        auto& dst = lights.shape_distribs.at(lgt.ist->shp);
        auto sample = sample_shape(lgt.ist->shp, dst, rel, ruv);
        return eval_point(
            lgt.ist, sample.first, sample.second, params.double_sided);
    } else if (lgt.env) {
        // BUG: this is not uniform sampling
        return eval_point(lgt.env, ruv);
    } else {
        throw std::runtime_error("should not have gotten here");
    }
    return {};
}

// Picks a point on a light.
trace_point sample_lights(const trace_lights& lights, const trace_point& pt,
    float rnl, float rne, const vec2f& ruv, const trace_params& params) {
    auto lidx = sample_discrete(lights.light_distrib, rnl);
    auto& lgt = lights.lights.at(lidx);
    return sample_light(lights, lgt, pt, rne, ruv, params);
}

// Intersects a ray with the scn and return the point (or env point).
trace_point intersect_scene(const scene* scn, const bvh_tree* bvh,
    const ray3f& ray, const trace_params& params) {
    auto iid = 0, eid = 0;
    auto euv = zero2f;
    auto ray_t = 0.0f;
    if (intersect_bvh(bvh, ray, false, ray_t, iid, eid, euv)) {
        return eval_point(scn->instances[iid], eid, euv, params.double_sided);
    } else if (!scn->environments.empty()) {
        return eval_point(
            scn->environments[0], eval_uv(scn->environments[0],
                                      transform_direction_inverse(
                                          scn->environments[0]->frame, ray.d)));
    } else {
        return {};
    }
}

// Test occlusion.
vec3f eval_transmission(const scene* scn, const bvh_tree* bvh,
    const trace_point& pt, const trace_point& lpt, const trace_params& params) {
    if (params.notransmission) {
        auto ray = make_segment(pt.pos, lpt.pos);
        auto iid = 0, eid = 0; auto ray_t = 0.0f; auto euv = zero2f;
        return (intersect_bvh(bvh, ray, true, ray_t, iid, eid, euv)) ? zero3f : vec3f{1, 1, 1};
    } else {
        auto cpt = pt;
        auto weight = vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            auto ray = make_segment(cpt.pos, lpt.pos);
            cpt = intersect_scene(scn, bvh, ray, params);
            if (!cpt.ist) break;
            weight *= cpt.kt + vec3f{1 - cpt.op, 1 - cpt.op, 1 - cpt.op};
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight.
float weight_mis(float w0, float w1) {
    if (w0 == 0) return w1;
    if (w1 == 0) return w0;
    return 1 / (1 / w0 + 1 / w1);
}

// Recursive path tracing.
vec3f trace_path(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    rng_state& rng, const trace_params& params) {
    if (lights.lights.empty()) return zero3f;
    auto pt = pt_;
    auto wo = wo_;

    // initialize
    auto l = eval_emission(pt, wo);
    auto weight = vec3f{1, 1, 1};

    // trace  path
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // opacity
        if (pt.op != 1) {
            if (next_rand1f(rng) < 1 - pt.op) {
                pt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
                l += weight * eval_emission(pt, wo);
                continue;
            }
        }

        // early exit
        if (!has_reflectance(pt)) break;

        // direct – light
        auto& lgt = lights.lights[sample_discrete(
            lights.light_distrib, next_rand1f(rng))];
        auto lpt = sample_light(
            lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
        auto lw = weight_light(lights, lpt, pt) * lights.light_distrib.back();
        auto lwi = normalize(lpt.pos - pt.pos);
        auto lke = eval_emission(lpt, -lwi);
        auto lbc = eval_brdfcos(pt, wo, lwi);
        auto lld = lke * lbc;
        if (lld != zero3f) {
            l += weight * lld * eval_transmission(scn, bvh, pt, lpt, params) *
                 weight_mis(lw, weight_brdf(pt, wo, lwi));
        }

        // direct – brdf
        auto bwi = sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng));
        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        auto bw = weight_brdf(pt, wo, bwi);
        auto bke = eval_emission(bpt, -bwi);
        auto bbc = eval_brdfcos(pt, wo, bwi);
        auto bld = bke * bbc;
        if (bld != zero3f) {
            l += weight * bld *
                 weight_mis(bw, weight_light(lights, bpt, pt) *
                                    lights.light_distrib.back());
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // deltas
        auto delta_prob = sample_delta_prob(pt, wo);
        auto delta = (delta_prob && next_rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);
        if (delta) {
            bwi =
                sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng), delta);
            bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
            l += weight * eval_emission(bpt, -bwi) *
                 eval_brdfcos(pt, wo, bwi, delta) *
                 weight_brdf(pt, wo, bwi, delta);
        }

        // continue path
        weight *=
            eval_brdfcos(pt, wo, bwi, delta) * weight_brdf(pt, wo, bwi, delta);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (next_rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        wo = -bwi;
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_nomis(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    rng_state& rng, const trace_params& params) {
    if (lights.lights.empty()) return zero3f;
    auto pt = pt_;
    auto wo = wo_;

    // initialize
    auto l = eval_emission(pt, wo);
    auto weight = vec3f{1, 1, 1};

    // trace  path
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // opacity
        if (pt.op != 1) {
            if (next_rand1f(rng) < 1 - pt.op) {
                pt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
                l += weight * eval_emission(pt, wo);
                continue;
            }
        }

        // early exit
        if (!has_reflectance(pt)) break;

        // direct
        if (!lights.lights.empty()) {
            auto& lgt = lights.lights[sample_discrete(
                lights.light_distrib, next_rand1f(rng))];
            auto lpt = sample_light(
                lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
            auto lwi = normalize(lpt.pos - pt.pos);
            auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                      weight_light(lights, lpt, pt) *
                      lights.light_distrib.back();
            if (ld != zero3f) {
                l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
            }
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // choose delta
        auto delta_prob = sample_delta_prob(pt, wo);
        auto delta = (delta_prob && next_rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);

        // continue path
        auto bwi =
            sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng), delta);
        weight *=
            eval_brdfcos(pt, wo, bwi, delta) * weight_brdf(pt, wo, bwi, delta);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (next_rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        pt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        wo = -bwi;
        if (delta) l += weight * eval_emission(pt, wo);
    }

    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, rng_state& rng, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);

    // ambient
    l += params.ambient * pt.kd;

    // direct
    for (auto& lgt : lights.lights) {
        auto lpt = sample_light(
            lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                  weight_light(lights, lpt, pt);
        if (ld == zero3f) continue;
        l += ld * eval_transmission(scn, bvh, pt, lpt, params);
    }

    // exit if needed
    if (bounce >= params.max_depth) return l;

    // reflection
    if (pt.ks != zero3f && !pt.rs) {
        auto wi = reflect(wo, pt.norm);
        auto rpt = intersect_scene(scn, bvh, make_ray(pt.pos, wi), params);
        l += pt.ks *
             trace_direct(scn, bvh, lights, rpt, -wi, bounce + 1, rng, params);
    }

    // opacity
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += pt.kt *
             trace_direct(scn, bvh, lights, opt, wo, bounce + 1, rng, params);
    }

    // opacity
    if (pt.op != 1) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l = pt.op * l + (1 - pt.op) * trace_direct(scn, bvh, lights, opt, wo,
                                          bounce + 1, rng, params);
    }

    // done
    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return trace_direct(scn, bvh, lights, pt, wo, 0, rng, params);
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, rng_state& rng, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);

    // brdf*light
    l += eval_brdfcos(pt, wo, wo) * pi;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += pt.kt *
             trace_eyelight(scn, bvh, lights, opt, wo, bounce + 1, rng, params);
    }
    if (pt.op != 1) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l = pt.op * l + (1 - pt.op) * trace_eyelight(scn, bvh, lights, opt, wo,
                                          bounce + 1, rng, params);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return trace_eyelight(scn, bvh, lights, pt, wo, 0, rng, params);
}

// Debug previewing.
vec3f trace_debug_normal(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return wn * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug frontfacing.
vec3f trace_debug_frontfacing(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return dot(wn, wo) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
}

// Debug previewing.
vec3f trace_debug_albedo(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return pt.kd + pt.ks + pt.kt;
}

// Debug previewing.
vec3f trace_debug_texcoord(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return {pt.texcoord.x, pt.texcoord.y, 0};
}

// Trace function
vec3f trace_func(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    switch (params.tracer) {
        case trace_type::eyelight:
            return trace_eyelight(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::direct:
            return trace_direct(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::pathtrace:
            return trace_path(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::pathtrace_nomis:
            return trace_path_nomis(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_albedo:
            return trace_debug_albedo(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_normal:
            return trace_debug_normal(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_frontfacing:
            return trace_debug_frontfacing(
                scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_texcoord:
            return trace_debug_texcoord(scn, bvh, lights, pt, wo, rng, params);
    }
}

// Trace a single sample
void trace_sample(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, trace_pixel& pxl, const trace_params& params) {
    auto crn = next_rand2f(pxl.rng);
    auto lrn = next_rand2f(pxl.rng);
    auto uv = vec2f{(pxl.i + crn.x) / (cam->aspect * params.resolution),
        1 - (pxl.j + crn.y) / params.resolution};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto pt = intersect_scene(scn, bvh, ray, params);
    if (!pt.ist && params.envmap_invisible) return;
    auto l = trace_func(scn, bvh, lights, pt, -ray.d, pxl.rng, params);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (max(l) > params.pixel_clamp) l = l * (params.pixel_clamp / max(l));
    pxl.acc += {l.x, l.y, l.z, 1};
    pxl.sample += 1;
}

// Trace the next nsamples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, int width, int height, std::vector<vec4f>& img, std::vector<trace_pixel>& pixels,
    int nsamples, const trace_params& params) {
    if (params.noparallel) {
        for (auto pid = 0; pid < pixels.size(); pid ++) {
            for (auto s = 0; s < params.nsamples; s++)
                trace_sample(scn, cam, bvh, lights, pixels[pid], params);
            img[pid] = pixels[pid].acc / pixels[pid].sample;
        }
    } else {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = std::vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(std::thread([=, &img, &pixels, &params]() {
                for (auto pid = tid; pid < pixels.size(); pid += nthreads) {
                    for (auto s = 0; s < nsamples; s++)
                        trace_sample(scn, cam, bvh, lights, pixels[pid], params);
                    img[pid] = pixels[pid].acc / pixels[pid].sample;
                }
            }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    }
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, int width, int height, std::vector<vec4f>& img, std::vector<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params, const std::function<void(int, int)>& callback) {
    pixels = make_trace_pixels(width, height, params);

    // render preview
    if (params.preview_resolution) {
        auto pparams = params;
        pparams.resolution = params.preview_resolution;
        pparams.nsamples = 1;
        auto pwidth = (int)std::round(cam->aspect * pparams.resolution);
        auto pheight = pparams.resolution;
        auto pimg = std::vector<vec4f>(pwidth*pheight);
        auto ppixels = make_trace_pixels(pwidth, pheight, pparams);
        trace_samples(scn, cam, bvh, lights, pwidth, pheight, pimg, ppixels, 1, pparams);
        img = resize_image(pwidth, pheight, pimg, width, height, ygl::resize_filter::box);
    } else {
        for (auto& p : img) p = zero4f;
    }
    if (callback) callback(0, 0);

    // start rendering
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &pixels, &stop_flag]() {
            for (auto s = 0; s < params.nsamples; s++) {
                for (auto pid = tid; pid < pixels.size(); pid += nthreads) {
                    if (stop_flag) return;
                    trace_sample(scn, cam, bvh, lights, pixels[pid], params);
                    img[pid] = pixels[pid].acc / pixels[pid].sample;
                    if (!pixels[pid].i && callback) callback(s, pixels[pid].j);
                }
            }
            if (!tid && callback) callback(params.nsamples, 0);
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

// Initialize trace lights
trace_lights make_trace_lights(const scene* scn) {
    auto lights = trace_lights();

    for (auto ist : scn->instances) {
        if (!ist->mat || ist->mat->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.ist = ist;
        lights.lights.push_back(lgt);
        if (lights.shape_distribs.find(ist->shp) == lights.shape_distribs.end())
            lights.shape_distribs[ist->shp] = sample_shape_cdf(ist->shp);
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.env = env;
        lights.lights.push_back(lgt);
    }

    lights.light_distrib = std::vector<float>(lights.lights.size());
    for (auto i = 0; i < lights.lights.size(); i++)
        lights.light_distrib[i] = i + 1;

    return lights;
}

// Initialize a rendering state
std::vector<trace_pixel> make_trace_pixels(int width, int height,
    const trace_params& params) {
    auto pixels = std::vector<trace_pixel>(width * height);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto pxl = trace_pixel();
            pxl.i = i;
            pxl.j = j;
            pxl.rng = make_rng(params.seed, (j * width + i) * 2 + 1);
            pixels[j * width + i] = pxl;
        }
    }
    return pixels;
}

}  // namespace ygl
