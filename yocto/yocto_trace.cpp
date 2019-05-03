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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& specular, float direction_cosine) {
    if (specular == zero3f) return zero3f;
    return specular +
           (1 - specular) *
               pow(clamp(1 - abs(direction_cosine), 0.0f, 1.0f), 5.0f);
}

// Evaluates the GGX distribution and geometric term
float eval_microfacetD(
    float roughness, const vec3f& normal, const vec3f& half_vector, bool ggx) {
    auto cosine = dot(normal, half_vector);
    if (cosine <= 0) return 0;
    auto roughness_square = roughness * roughness;
    auto cosine_square    = cosine * cosine;
    auto tangent_square   = clamp01(1 - cosine_square) / cosine_square;
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
    const vec3f& half_vector, const vec3f& direction, bool ggx) {
    auto cosine = dot(normal, direction);
    if (dot(half_vector, direction) * cosine <= 0) return 0;
    auto roughness_square = roughness * roughness;
    auto cosine_square    = cosine * cosine;
    auto tangent_square   = clamp01(1 - cosine_square) / cosine_square;
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
    bool ggx) {
    return evaluate_microfacetG1(
               roughness, normal, half_vector, outgoing, ggx) *
           evaluate_microfacetG1(roughness, normal, half_vector, incoming, ggx);
}
vec3f sample_microfacet(
    float roughness, const vec3f& normal, const vec2f& rn, bool ggx) {
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
    auto radius            = sqrt(clamp01(1 - cosine_square));
    auto local_half_vector = vec3f{
        cos(phi) * radius, sin(phi) * radius, cosine};
    return transform_direction(make_basis_fromz(normal), local_half_vector);
}
float sample_microfacet_pdf(
    float roughness, const vec3f& normal, const vec3f& half_vector, bool ggx) {
    auto cosine = dot(normal, half_vector);
    if (cosine < 0) return 0;
    return eval_microfacetD(roughness, normal, half_vector, ggx) * cosine;
}

// Phong exponent to roughness.
float exponent_to_roughness(float exponent) {
    return sqrtf(2 / (exponent + 2));
}

// Specular to  eta.
vec3f reflectivity_to_eta(const vec3f& reflectivity) {
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
        return {float_max, channel};
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
    if (abs(g) < 1e-3) {
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

// Tabulated ior for metals
// https://github.com/tunabrain/tungsten
const unordered_map<string, pair<vec3f, vec3f>> metal_ior_table = {
    {"a-C", {{2.9440999183f, 2.2271502925f, 1.9681668794f},
                {0.8874329109f, 0.7993216383f, 0.8152862927f}}},
    {"Ag", {{0.1552646489f, 0.1167232965f, 0.1383806959f},
               {4.8283433224f, 3.1222459278f, 2.1469504455f}}},
    {"Al", {{1.6574599595f, 0.8803689579f, 0.5212287346f},
               {9.2238691996f, 6.2695232477f, 4.8370012281f}}},
    {"AlAs", {{3.6051023902f, 3.2329365777f, 2.2175611545f},
                 {0.0006670247f, -0.0004999400f, 0.0074261204f}}},
    {"AlSb", {{-0.0485225705f, 4.1427547893f, 4.6697691348f},
                 {-0.0363741915f, 0.0937665154f, 1.3007390124f}}},
    {"Au", {{0.1431189557f, 0.3749570432f, 1.4424785571f},
               {3.9831604247f, 2.3857207478f, 1.6032152899f}}},
    {"Be", {{4.1850592788f, 3.1850604423f, 2.7840913457f},
               {3.8354398268f, 3.0101260162f, 2.8690088743f}}},
    {"Cr", {{4.3696828663f, 2.9167024892f, 1.6547005413f},
               {5.2064337956f, 4.2313645277f, 3.7549467933f}}},
    {"CsI", {{2.1449030413f, 1.7023164587f, 1.6624194173f},
                {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
    {"Cu", {{0.2004376970f, 0.9240334304f, 1.1022119527f},
               {3.9129485033f, 2.4528477015f, 2.1421879552f}}},
    {"Cu2O", {{3.5492833755f, 2.9520622449f, 2.7369202137f},
                 {0.1132179294f, 0.1946659670f, 0.6001681264f}}},
    {"CuO", {{3.2453822204f, 2.4496293965f, 2.1974114493f},
                {0.5202739621f, 0.5707372756f, 0.7172250613f}}},
    {"d-C", {{2.7112524747f, 2.3185812849f, 2.2288565009f},
                {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
    {"Hg", {{2.3989314904f, 1.4400254917f, 0.9095512090f},
               {6.3276269444f, 4.3719414152f, 3.4217899270f}}},
    {"HgTe", {{4.7795267752f, 3.2309984581f, 2.6600252401f},
                 {1.6319827058f, 1.5808189339f, 1.7295753852f}}},
    {"Ir", {{3.0864098394f, 2.0821938440f, 1.6178866805f},
               {5.5921510077f, 4.0671757150f, 3.2672611269f}}},
    {"K", {{0.0640493070f, 0.0464100621f, 0.0381842017f},
              {2.1042155920f, 1.3489364357f, 0.9132113889f}}},
    {"Li", {{0.2657871942f, 0.1956102432f, 0.2209198538f},
               {3.5401743407f, 2.3111306542f, 1.6685930000f}}},
    {"MgO", {{2.0895885542f, 1.6507224525f, 1.5948759692f},
                {0.0000000000f, -0.0000000000f, 0.0000000000f}}},
    {"Mo", {{4.4837010280f, 3.5254578255f, 2.7760769438f},
               {4.1111307988f, 3.4208716252f, 3.1506031404f}}},
    {"Na", {{0.0602665320f, 0.0561412435f, 0.0619909494f},
               {3.1792906496f, 2.1124800781f, 1.5790940266f}}},
    {"Nb", {{3.4201353595f, 2.7901921379f, 2.3955856658f},
               {3.4413817900f, 2.7376437930f, 2.5799132708f}}},
    {"Ni", {{2.3672753521f, 1.6633583302f, 1.4670554172f},
               {4.4988329911f, 3.0501643957f, 2.3454274399f}}},
    {"Rh", {{2.5857954933f, 1.8601866068f, 1.5544279524f},
               {6.7822927110f, 4.7029501026f, 3.9760892461f}}},
    {"Se-e", {{5.7242724833f, 4.1653992967f, 4.0816099264f},
                 {0.8713747439f, 1.1052845009f, 1.5647788766f}}},
    {"Se", {{4.0592611085f, 2.8426947380f, 2.8207582835f},
               {0.7543791750f, 0.6385150558f, 0.5215872029f}}},
    {"SiC", {{3.1723450205f, 2.5259677964f, 2.4793623897f},
                {0.0000007284f, -0.0000006859f, 0.0000100150f}}},
    {"SnTe", {{4.5251865890f, 1.9811525984f, 1.2816819226f},
                 {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
    {"Ta", {{2.0625846607f, 2.3930915569f, 2.6280684948f},
               {2.4080467973f, 1.7413705864f, 1.9470377016f}}},
    {"Te-e", {{7.5090397678f, 4.2964603080f, 2.3698732430f},
                 {5.5842076830f, 4.9476231084f, 3.9975145063f}}},
    {"Te", {{7.3908396088f, 4.4821028985f, 2.6370708478f},
               {3.2561412892f, 3.5273908133f, 3.2921683116f}}},
    {"ThF4", {{1.8307187117f, 1.4422274283f, 1.3876488528f},
                 {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
    {"TiC", {{3.7004673762f, 2.8374356509f, 2.5823030278f},
                {3.2656905818f, 2.3515586388f, 2.1727857800f}}},
    {"TiN", {{1.6484691607f, 1.1504482522f, 1.3797795097f},
                {3.3684596226f, 1.9434888540f, 1.1020123347f}}},
    {"TiO2-e", {{3.1065574823f, 2.5131551146f, 2.5823844157f},
                   {0.0000289537f, -0.0000251484f, 0.0001775555f}}},
    {"TiO2", {{3.4566203131f, 2.8017076558f, 2.9051485020f},
                 {0.0001026662f, -0.0000897534f, 0.0006356902f}}},
    {"VC", {{3.6575665991f, 2.7527298065f, 2.5326814570f},
               {3.0683516659f, 2.1986687713f, 1.9631816252f}}},
    {"VN", {{2.8656011588f, 2.1191817791f, 1.9400767149f},
               {3.0323264950f, 2.0561075580f, 1.6162930914f}}},
    {"V", {{4.2775126218f, 3.5131538236f, 2.7611257461f},
              {3.4911844504f, 2.8893580874f, 3.1116965117f}}},
    {"W", {{4.3707029924f, 3.3002972445f, 2.9982666528f},
              {3.5006778591f, 2.6048652781f, 2.2731930614f}}},
};

// Get a complex ior table with keys the metal name and values (eta, etak)
pair<vec3f, vec3f> get_conductor_eta(const string& name) {
    if (metal_ior_table.find(name) == metal_ior_table.end())
        return {zero3f, zero3f};
    return metal_ior_table.at(name);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Set non-rigid frames as default
constexpr bool trace_non_rigid_frames = true;

// Emission lobe
struct trace_emission {
    enum struct type_t { diffuse };
    type_t type   = type_t::diffuse;
    vec3f  weight = zero3f;
};
using trace_emissions = short_vector<trace_emission, 8>;
// Scattering lobe
struct trace_bsdf {
    enum struct type_t {
        diffuse,
        reflection,
        transmission,
        transparency,
        translucency
    };
    type_t type      = type_t::diffuse;
    vec3f  weight    = zero3f;
    vec3f  eta       = zero3f;
    vec3f  etak      = zero3f;
    float  roughness = 0;
    float  samplew   = 0;
    float  pdf       = 0;
};
using trace_bsdfs = short_vector<trace_bsdf, 8>;
// Delta lobe
struct trace_delta {
    enum struct type_t { reflection, transmission, transparency, passthrough };
    type_t type    = type_t::reflection;
    vec3f  weight  = zero3f;
    vec3f  eta     = zero3f;
    vec3f  etak    = zero3f;
    float  samplew = 0;
    float  pdf     = 0;
};
using trace_deltas = short_vector<trace_delta, 8>;
// Scattering medium
struct trace_medium {
    enum struct type_t { phaseg };
    type_t type     = type_t::phaseg;
    vec3f  weight   = zero3f;
    vec3f  emission = zero3f;
    vec3f  density  = zero3f;
    vec3f  albedo   = zero3f;
    float  phaseg   = 0;
    float  samplew  = 0;
    float  pdf      = 0;
};
using trace_mediums = short_vector<trace_medium, 8>;

void eval_material(trace_emissions& emissions, trace_bsdfs& bsdfs,
    trace_deltas& deltas, trace_mediums& mediums, const material_point& point_,
    const vec4f& shape_color, const vec3f& normal, const vec3f& outgoing) {
    auto point = point_;
    point.emission_color *= xyz(shape_color);
    point.base_color *= xyz(shape_color);
    point.opacity_factor *= shape_color.w;
    point.specular_roughness = point.specular_roughness *
                               point.specular_roughness;
    if (point.specular_roughness)
        point.specular_roughness = clamp(
            point.specular_roughness, 0.03f * 0.03f, 1.0f);
    point.coat_roughness = point.coat_roughness * point.coat_roughness;
    if (point.coat_roughness)
        point.coat_roughness = clamp(point.coat_roughness, 0.03f * 0.03f, 1.0f);
    if (point.opacity_factor > 0.999f) point.opacity_factor = 1;
    auto weight = vec3f{1};
    if (point.opacity_factor < 0.999f) {
        auto lweight = vec3f{1 - point.opacity_factor};
        deltas.push_back({trace_delta::type_t::passthrough, lweight, zero3f,
            zero3f, max(lweight), 0});
        weight *= point.opacity_factor;
    }
    if (point.coat_factor) {
        auto roughness = point.coat_roughness;
        auto eta       = vec3f{point.coat_ior};
        auto fresnel   = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
        auto lweight   = weight;
        if (lweight != zero3f) {
            if (roughness) {
                bsdfs.push_back({trace_bsdf::type_t::reflection, lweight, eta,
                    zero3f, roughness, max(lweight * fresnel), 0});
            } else {
                deltas.push_back({trace_delta::type_t::reflection, lweight, eta,
                    zero3f, max(lweight * fresnel), 0});
            }
        }
        weight *= point.coat_color * point.coat_factor * (1 - fresnel);
    }
    if (point.emission_factor) {
        auto lweight = weight * point.emission_factor * point.emission_color;
        if (lweight != zero3f) {
            emissions.push_back({trace_emission::type_t::diffuse, lweight});
        }
    }
    if (point.metallic_factor) {
        auto roughness = point.specular_roughness;
        auto eta = reflectivity_to_eta(point.base_color), etak = zero3f;
        // auto [eta1, etak1] = reflectivity_to_eta(point.base_color,
        // point.specular_color); auto eta2 = vec3f{0.1431189557f,
        // 0.3749570432f, 1.4424785571f}; auto etak2 =
        // vec3f{3.9831604247f, 2.3857207478f, 1.6032152899f}; etak = zero3f;
        // print("1: {} {}\n 2: {} {}\n 3: {} {}\n", eta1, etak1, eta2, etak2,
        // eta3, etak3); auto eta = eta2, etak = etak3; auto fresnel =
        // fresnel_conductor(eta, etak, abs(dot(normal, outgoing)));
        auto fresnel = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
        auto lweight = weight * point.metallic_factor;
        if (lweight != zero3f) {
            if (roughness) {
                bsdfs.push_back({trace_bsdf::type_t::reflection, lweight, eta,
                    etak, roughness, max(lweight * fresnel), 0});
            } else {
                deltas.push_back({trace_delta::type_t::reflection, lweight, eta,
                    etak, max(lweight * fresnel), 0});
            }
        }
        weight *= 1 - point.metallic_factor;
    }
    if (point.transmission_factor) {
        auto roughness = point.specular_roughness;
        auto eta       = vec3f{point.specular_ior};
        if (point.thin_walled && !point.specular_factor) eta = {1, 1, 1};
        auto fresnel = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
        auto lweight = weight * point.transmission_factor;
        if (point.thin_walled || !point.transmission_depth) {
            lweight *= point.transmission_color;
        } else {
            auto density = -log(
                               clamp(point.transmission_color, 0.0001f, 1.0f)) /
                           point.transmission_depth;
            // auto density = point.transmission_color *
            // point.transmission_depth;
            mediums.push_back(
                {trace_medium::type_t::phaseg, {1, 1, 1}, zero3f, density,
                    point.transmission_scatter, point.transmission_anisotropy,
                    max(point.transmission_scatter)});
        }
        if (lweight != zero3f) {
            if (roughness) {
                bsdfs.push_back(
                    {point.thin_walled ? trace_bsdf::type_t::transparency
                                       : trace_bsdf::type_t::transmission,
                        lweight, eta, zero3f, roughness,
                        max(lweight * (1 - fresnel)), 0});
            } else {
                deltas.push_back(
                    {point.thin_walled ? trace_delta::type_t::transparency
                                       : trace_delta::type_t::transmission,
                        lweight, eta, zero3f, max(lweight * (1 - fresnel)), 0});
            }
        }
    }
    if (point.specular_factor) {
        auto roughness = point.specular_roughness;
        auto eta       = point.specular_ior_from_color
                       ? reflectivity_to_eta(
                             point.specular_factor * point.specular_color)
                       : vec3f{point.specular_ior};
        auto fresnel = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
        auto lweight = weight *
                       (point.specular_ior_from_color
                               ? vec3f{1, 1, 1}
                               : point.specular_factor * point.specular_color);
        if (lweight != zero3f) {
            if (roughness) {
                bsdfs.push_back({trace_bsdf::type_t::reflection, lweight, eta,
                    zero3f, roughness, max(lweight * fresnel), 0});
            } else {
                deltas.push_back({trace_delta::type_t::reflection, lweight, eta,
                    zero3f, max(lweight * fresnel), 0});
            }
        }
        weight *= 1 - fresnel;
    }
    if (point.subsurface_factor) {
        auto roughness = point.subsurface_factor;
        auto lweight   = weight * point.subsurface_factor;
        if (point.thin_walled) {
            lweight *= point.subsurface_color;
            bsdfs.push_back({trace_bsdf::type_t::translucency, lweight, zero3f,
                zero3f, roughness, max(lweight), 0});
        } else {
            auto density = 1 /
                           (point.subsurface_radius * point.subsurface_scale);
            mediums.push_back({trace_medium::type_t::phaseg, {1, 1, 1},
                point.subsurface_emission, density, point.subsurface_color,
                point.subsurface_anisotropy, max(point.subsurface_color)});
            deltas.push_back({trace_delta::type_t::transparency, lweight,
                zero3f, zero3f, max(lweight), 0});
        }
        weight *= 1 - point.subsurface_factor;
    }
    if (point.diffuse_factor) {
        auto lweight = weight * point.diffuse_factor * point.base_color;
        if (lweight != zero3f) {
            bsdfs.push_back({trace_bsdf::type_t::diffuse, lweight, zero3f,
                zero3f, 0, max(lweight), 0});
        }
    }
    auto normalize_weights = [](auto& lobes) {
        auto weight = 0.0f;
        for (auto& lobe : lobes) weight += lobe.samplew;
        for (auto& lobe : lobes) lobe.pdf = lobe.samplew / weight;
    };
    normalize_weights(bsdfs);
    normalize_weights(deltas);
    normalize_weights(mediums);
}  // namespace yocto

vec3f eval_emission(const trace_emissions& emissions, const vec3f& normal,
    const vec3f& outgoing) {
    auto emission = zero3f;
    for (auto& lobe : emissions) emission += lobe.weight;
    return emission;
}
vec3f eval_emission(
    const trace_mediums& mediums, const vec3f& normal, const vec3f& outgoing) {
    auto emission = zero3f;
    for (auto& lobe : mediums) emission += lobe.weight * lobe.emission;
    return emission;
}

// Check if we are on the same side of the hemisphere
bool same_hemisphere(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    return dot(normal, outgoing) * dot(normal, incoming) > 0;
}
bool other_hemisphere(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    return dot(normal, outgoing) * dot(normal, incoming) < 0;
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
vec3f eval_diffuse_reflection(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return zero3f;
    return vec3f{abs(dot(normal, incoming)) / pif};
}
vec3f eval_diffuse_translucency(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    return vec3f{abs(dot(normal, incoming)) / pif};
}
vec3f eval_microfacet_reflection(float roughness, const vec3f& eta,
    const vec3f& etak, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return zero3f;
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    auto halfway   = normalize(incoming + outgoing);
    auto F         = etak == zero3f
                 ? fresnel_dielectric(eta, abs(dot(halfway, outgoing)))
                 : fresnel_conductor(eta, etak, abs(dot(halfway, outgoing)));
    auto D = eval_microfacetD(roughness, up_normal, halfway);
    auto G = eval_microfacetG(
        roughness, up_normal, halfway, outgoing, incoming);
    return F * D * G / abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
           abs(dot(normal, incoming));
}
vec3f eval_delta_reflection(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return zero3f;
    auto F = etak == zero3f
                 ? fresnel_dielectric(eta, abs(dot(normal, outgoing)))
                 : fresnel_conductor(eta, etak, abs(dot(normal, outgoing)));
    return F;
}
vec3f eval_microfacet_transmission(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    auto up_normal      = dot(outgoing, normal) > 0 ? normal : -normal;
    auto halfway_vector = dot(outgoing, normal) > 0
                              ? -(outgoing + eta * incoming)
                              : (eta * outgoing + incoming);
    auto halfway = normalize(halfway_vector);
    auto F       = fresnel_dielectric(eta, abs(dot(halfway, outgoing)));
    auto D       = eval_microfacetD(roughness, up_normal, halfway);
    auto G       = eval_microfacetG(
        roughness, up_normal, halfway, outgoing, incoming);

    auto dot_terms = (dot(outgoing, halfway) * dot(incoming, halfway)) /
                     (dot(outgoing, normal) * dot(incoming, normal));

    auto numerator   = (1 - F) * D * G;
    auto denominator = dot(halfway_vector, halfway_vector);

    // [Walter 2007] equation 21
    return abs(dot_terms) * numerator / denominator *
           abs(dot(normal, incoming));
}
vec3f eval_delta_transmission(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    auto F = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
    return (1 - F);
}
vec3f eval_microfacet_transparency(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    auto up_normal = dot(outgoing, normal) > 0 ? normal : -normal;
    auto ir        = reflect(-incoming, up_normal);
    auto halfway   = normalize(ir + outgoing);
    auto F         = fresnel_dielectric(eta, abs(dot(halfway, outgoing)));
    auto D         = eval_microfacetD(roughness, up_normal, halfway);
    auto G = eval_microfacetG(roughness, up_normal, halfway, outgoing, ir);
    return (1 - F) * D * G /
           abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
           abs(dot(normal, incoming));
}
vec3f eval_delta_transparency(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    if (eta == zero3f || eta == vec3f{1, 1, 1}) {
        return {1, 1, 1};
    } else {
        auto F = fresnel_dielectric(eta, abs(dot(normal, outgoing)));
        return (1 - F);
    }
}
vec3f eval_delta_passthrough(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return zero3f;
    return {1, 1, 1};
}
vec3f eval_volume_scattering(const vec3f& albedo, float phaseg,
    const vec3f& outgoing, const vec3f& incoming) {
    return albedo * eval_phasefunction(dot(outgoing, incoming), phaseg);
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
vec3f eval_scattering(const trace_bsdfs& bsdfs, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto scattering = zero3f;
    for (auto& lobe : bsdfs) {
        if (lobe.weight == zero3f) continue;
        switch (lobe.type) {
            case trace_bsdf::type_t::diffuse: {
                scattering += lobe.weight *
                              eval_diffuse_reflection(
                                  lobe.roughness, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::translucency: {
                scattering += lobe.weight *
                              eval_diffuse_translucency(
                                  lobe.roughness, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::reflection: {
                scattering += lobe.weight *
                              eval_microfacet_reflection(lobe.roughness,
                                  lobe.eta, zero3f, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::transmission: {
                scattering += lobe.weight *
                              eval_microfacet_transmission(lobe.roughness,
                                  lobe.eta, normal, outgoing, incoming);

            } break;
            case trace_bsdf::type_t::transparency: {
                scattering += lobe.weight *
                              eval_microfacet_transparency(lobe.roughness,
                                  lobe.eta, normal, outgoing, incoming);
            } break;
        }
    }
    return scattering;
}
vec3f eval_scattering(const trace_deltas& deltas, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto scattering = zero3f;
    for (auto& lobe : deltas) {
        if (lobe.weight == zero3f) continue;
        switch (lobe.type) {
            case trace_delta::type_t::reflection: {
                scattering += lobe.weight * eval_delta_reflection(lobe.eta,
                                                zero3f, normal, outgoing,
                                                incoming);
            } break;
            case trace_delta::type_t::transmission: {
                scattering += lobe.weight * eval_delta_transmission(lobe.eta,
                                                normal, outgoing, incoming);

            } break;
            case trace_delta::type_t::transparency: {
                scattering += lobe.weight * eval_delta_transparency(lobe.eta,
                                                normal, outgoing, incoming);
            } break;
            case trace_delta::type_t::passthrough: {
                scattering += lobe.weight * eval_delta_passthrough(
                                                normal, outgoing, incoming);
            } break;
            default: break;
        }
    }
    return scattering;
}
vec3f eval_scattering(const trace_mediums& mediums, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto scattering = zero3f;
    for (auto& lobe : mediums) {
        if (lobe.weight == zero3f) continue;
        switch (lobe.type) {
            case trace_medium::type_t::phaseg: {
                scattering += lobe.weight * eval_volume_scattering(lobe.albedo,
                                                lobe.phaseg, outgoing,
                                                incoming);
            } break;
        }
    }
    return scattering;
}

// Picks a direction based on the BRDF
vec3f sample_diffuse_reflection(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    return sample_hemisphere(up_normal, rn);
}
vec3f sample_diffuse_translucency(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    return sample_hemisphere(-up_normal, rn);
}
vec3f sample_microfacet_reflection(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    auto halfway   = sample_microfacet(roughness, up_normal, rn);
    return reflect(outgoing, halfway);
}
vec3f sample_delta_reflection(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    return reflect(outgoing, up_normal);
}
vec3f sample_microfacet_transmission(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    auto halfway   = sample_microfacet(roughness, up_normal, rn);
    return refract(outgoing, halfway,
        dot(normal, outgoing) > 0 ? 1 / mean(eta) : mean(eta));
}
vec3f sample_delta_transmission(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    return refract(outgoing, up_normal,
        dot(normal, outgoing) > 0 ? 1 / mean(eta) : mean(eta));
}
vec3f sample_microfacet_transparency(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
    auto up_normal = dot(normal, outgoing) > 0 ? normal : -normal;
    auto halfway   = sample_microfacet(roughness, up_normal, rn);
    auto ir        = reflect(outgoing, halfway);
    return -reflect(ir, up_normal);
}
vec3f sample_delta_transparency(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec2f& rn) {
    return -outgoing;
}
vec3f sample_delta_passthrough(
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
    return -outgoing;
}
vec3f sample_volume_scattering(const vec3f& albedo, float phaseg,
    const vec3f& normal, const vec3f& outgoing, const vec2f& rn) {
    auto direction = sample_phasefunction(phaseg, rn);
    return make_basis_fromz(-outgoing) * direction;
}

// Picks a direction based on the BRDF
vec3f sample_scattering(const trace_bsdfs& bsdfs, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
    // keep a weight sum to pick a lobe
    auto weight_sum = 0.0f;
    for (auto& lobe : bsdfs) {
        if (lobe.pdf == 0) continue;

        // check if we pick this lobe
        weight_sum += lobe.pdf;
        if (rnl > weight_sum) continue;

        switch (lobe.type) {
            case trace_bsdf::type_t::diffuse: {
                return sample_diffuse_reflection(
                    lobe.roughness, normal, outgoing, rn);
            } break;
            case trace_bsdf::type_t::translucency: {
                return sample_diffuse_translucency(
                    lobe.roughness, normal, outgoing, rn);
            } break;
            case trace_bsdf::type_t::reflection: {
                return sample_microfacet_reflection(
                    lobe.roughness, lobe.eta, normal, outgoing, rn);
            } break;
            case trace_bsdf::type_t::transmission: {
                return sample_microfacet_transmission(
                    lobe.roughness, lobe.eta, normal, outgoing, rn);
            } break;
            case trace_bsdf::type_t::transparency: {
                return sample_microfacet_transparency(
                    lobe.roughness, lobe.eta, normal, outgoing, rn);
            } break;
        }
    }

    // something went wrong if we got here
    return zero3f;
}
vec3f sample_scattering(const trace_deltas& deltas, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
    // keep a weight sum to pick a lobe
    auto weight_sum = 0.0f;
    for (auto& lobe : deltas) {
        if (lobe.pdf == 0) continue;

        // check if we pick this lobe
        weight_sum += lobe.pdf;
        if (rnl > weight_sum) continue;

        switch (lobe.type) {
            case trace_delta::type_t::reflection: {
                return sample_delta_reflection(lobe.eta, normal, outgoing, rn);
            } break;
            case trace_delta::type_t::transmission: {
                return sample_delta_transmission(
                    lobe.eta, normal, outgoing, rn);
            } break;
            case trace_delta::type_t::transparency: {
                return sample_delta_transparency(
                    lobe.eta, normal, outgoing, rn);
            } break;
            case trace_delta::type_t::passthrough: {
                return sample_delta_passthrough(normal, outgoing, rn);
            } break;
        }
    }

    // something went wrong if we got here
    return zero3f;
}
vec3f sample_scattering(const trace_mediums& mediums, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
    // keep a weight sum to pick a lobe
    auto weight_sum = 0.0f;
    for (auto& lobe : mediums) {
        if (lobe.pdf == 0) continue;

        // check if we pick this lobe
        weight_sum += lobe.pdf;
        if (rnl > weight_sum) continue;

        switch (lobe.type) {
            case trace_medium::type_t::phaseg: {
                return sample_volume_scattering(
                    lobe.albedo, lobe.phaseg, normal, outgoing, rn);
            } break;
        }
    }

    // something went wrong if we got here
    return zero3f;
}

// Compute the weight for sampling the BRDF
float sample_diffuse_reflection_pdf(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return 0;
    return abs(dot(normal, incoming)) / pif;
}
float sample_diffuse_translucency_pdf(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    return abs(dot(normal, incoming)) / pif;
}
float sample_microfacet_reflection_pdf(float roughness, const vec3f& eta,
    const vec3f& etak, const vec3f& normal, const vec3f& outgoing,
    const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return 0;
    auto up_normal = dot(normal, outgoing) >= 0 ? normal : -normal;
    auto halfway   = normalize(incoming + outgoing);
    auto d         = sample_microfacet_pdf(roughness, up_normal, halfway);
    auto jacobian  = 0.25f / abs(dot(outgoing, halfway));
    return d * jacobian;
}
float sample_delta_reflection_pdf(const vec3f& eta, const vec3f& etak,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!same_hemisphere(normal, outgoing, incoming)) return 0;
    return 1;
}
float sample_microfacet_transmission_pdf(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    auto up_normal      = dot(outgoing, normal) > 0 ? normal : -normal;
    auto halfway_vector = dot(outgoing, normal) > 0
                              ? -(outgoing + eta * incoming)
                              : (eta * outgoing + incoming);
    auto halfway = normalize(halfway_vector);
    auto d       = sample_microfacet_pdf(roughness, up_normal, halfway);
    // [Walter 2007] equation 17
    auto jacobian = abs(dot(halfway, incoming)) /
                    dot(halfway_vector, halfway_vector);
    return d * jacobian;
}
float sample_delta_transmission_pdf(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    return 1;
}
float sample_microfacet_transparency_pdf(float roughness, const vec3f& eta,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    auto up_normal = dot(outgoing, normal) > 0 ? normal : -normal;
    auto ir        = reflect(-incoming, up_normal);
    auto halfway   = normalize(ir + outgoing);
    auto d         = sample_microfacet_pdf(roughness, up_normal, halfway);
    auto jacobian  = 0.25f / abs(dot(outgoing, halfway));
    return d * jacobian;
}
float sample_delta_transparency_pdf(const vec3f& eta, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    return 1;
}
float sample_delta_passthrough_pdf(
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!other_hemisphere(normal, outgoing, incoming)) return 0;
    return 1;
}
float sample_volume_scattering_pdf(const vec3f& albedo, float phaseg,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    return eval_phasefunction(dot(outgoing, incoming), phaseg);
}

// Compute the weight for sampling the BRDF
float sample_scattering_pdf(const trace_bsdfs& bsdfs, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto pdf = 0.0f;
    for (auto& lobe : bsdfs) {
        if (lobe.pdf == 0) continue;

        switch (lobe.type) {
            case trace_bsdf::type_t::diffuse: {
                pdf += lobe.pdf * sample_diffuse_reflection_pdf(lobe.roughness,
                                      normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::translucency: {
                pdf += lobe.pdf *
                       sample_diffuse_translucency_pdf(
                           lobe.roughness, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::reflection: {
                pdf += lobe.pdf *
                       sample_microfacet_reflection_pdf(lobe.roughness,
                           lobe.eta, lobe.etak, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::transmission: {
                pdf += lobe.pdf *
                       sample_microfacet_transmission_pdf(lobe.roughness,
                           lobe.eta, normal, outgoing, incoming);
            } break;
            case trace_bsdf::type_t::transparency: {
                pdf += lobe.pdf *
                       sample_microfacet_transparency_pdf(lobe.roughness,
                           lobe.eta, normal, outgoing, incoming);
            } break;
        }
    }

    return pdf;
}
float sample_scattering_pdf(const trace_deltas& deltas, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto pdf = 0.0f;
    for (auto& lobe : deltas) {
        if (lobe.pdf == 0) continue;

        switch (lobe.type) {
            case trace_delta::type_t::reflection: {
                pdf += lobe.pdf * sample_delta_reflection_pdf(lobe.eta,
                                      lobe.etak, normal, outgoing, incoming);
            } break;
            case trace_delta::type_t::transmission: {
                pdf += lobe.pdf * sample_delta_transmission_pdf(
                                      lobe.eta, normal, outgoing, incoming);
            } break;
            case trace_delta::type_t::transparency: {
                pdf += lobe.pdf * sample_delta_transparency_pdf(
                                      lobe.eta, normal, outgoing, incoming);
            } break;
            case trace_delta::type_t::passthrough: {
                pdf += lobe.pdf *
                       sample_delta_passthrough_pdf(normal, outgoing, incoming);
            } break;
        }
    }

    return pdf;
}
float sample_scattering_pdf(const trace_mediums& mediums, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    auto pdf = 0.0f;
    for (auto& lobe : mediums) {
        if (lobe.pdf == 0) continue;

        switch (lobe.type) {
            case trace_medium::type_t::phaseg: {
                pdf += lobe.pdf * sample_volume_scattering_pdf(lobe.albedo,
                                      lobe.phaseg, normal, outgoing, incoming);
            } break;
        }
    }

    return pdf;
}

// Sample pdf for an environment.
float sample_environment_pdf(const yocto_scene& scene,
    const trace_lights& lights, int environment_id, const vec3f& incoming) {
    auto& environment = scene.environments[environment_id];
    if (environment.emission_texture >= 0) {
        auto& elements_cdf =
            lights.environment_cdfs[environment.emission_texture];
        auto& emission_texture = scene.textures[environment.emission_texture];
        auto  size             = texture_size(emission_texture);
        auto  texcoord         = eval_texcoord(environment, incoming);
        auto  i    = clamp((int)(texcoord.x * size.x), 0, size.x - 1);
        auto  j    = clamp((int)(texcoord.y * size.y), 0, size.y - 1);
        auto  prob = sample_discrete_pdf(elements_cdf, j * size.x + i) /
                    elements_cdf.back();
        auto angle = (2 * pif / size.x) * (pif / size.y) *
                     sin(pif * (j + 0.5f) / size.y);
        return prob / angle;
    } else {
        return 1 / (4 * pif);
    }
}

// Picks a point on an environment.
vec3f sample_environment(const yocto_scene& scene, const trace_lights& lights,
    int environment_id, float rel, const vec2f& ruv) {
    auto& environment = scene.environments[environment_id];
    if (environment.emission_texture >= 0) {
        auto& elements_cdf =
            lights.environment_cdfs[environment.emission_texture];
        auto& emission_texture = scene.textures[environment.emission_texture];
        auto  idx              = sample_discrete(elements_cdf, rel);
        auto  size             = texture_size(emission_texture);
        auto  u                = (idx % size.x + 0.5f) / size.x;
        auto  v                = (idx / size.x + 0.5f) / size.y;
        return eval_direction(environment, {u, v});
    } else {
        return sample_sphere(ruv);
    }
}

// Picks a point on a light.
vec3f sample_light(const yocto_scene& scene, const trace_lights& lights,
    int instance_id, const vec3f& p, float rel, const vec2f& ruv) {
    auto& instance                = scene.instances[instance_id];
    auto& shape                   = scene.shapes[instance.shape];
    auto& elements_cdf            = lights.shape_cdfs[instance.shape];
    auto [element_id, element_uv] = sample_shape(shape, elements_cdf, rel, ruv);
    return normalize(
        eval_position(scene, instance, element_id, element_uv) - p);
}

// Sample pdf for a light point.
float sample_light_pdf(const yocto_scene& scene, const trace_lights& lights,
    int instance_id, const bvh_scene& bvh, const vec3f& position,
    const vec3f& direction) {
    auto& instance = scene.instances[instance_id];
    auto& material = scene.materials[instance.material];
    if (!material.emission_factor) return 0;
    auto& elements_cdf = lights.shape_cdfs[instance.shape];
    // check all intersection
    auto pdf           = 0.0f;
    auto next_position = position;
    for (auto bounce = 0; bounce < 100; bounce++) {
        auto isec = bvh_intersection{};
        if (!intersect_bvh(scene, bvh, instance_id,
                make_ray(next_position, direction), isec))
            break;
        // accumulate pdf
        auto& instance       = scene.instances[isec.instance_id];
        auto  light_position = eval_position(
            scene, instance, isec.element_id, isec.element_uv);
        auto light_normal = eval_normal(scene, instance, isec.element_id,
            isec.element_uv, trace_non_rigid_frames);
        // prob triangle * area triangle = area triangle mesh
        auto area = elements_cdf.back();
        pdf += distance_squared(light_position, position) /
               (abs(dot(light_normal, direction)) * area);
        // continue
        next_position = light_position + direction * 1e-3f;
    }
    return pdf;
}

// Sample lights wrt solid angle
vec3f sample_lights(const yocto_scene& scene, const trace_lights& lights,
    const bvh_scene& bvh, const vec3f& position, float rl, float rel,
    const vec2f& ruv) {
    auto light_id = sample_uniform(
        lights.instances.size() + lights.environments.size(), rl);
    if (light_id < lights.instances.size()) {
        auto instance = lights.instances[light_id];
        return sample_light(scene, lights, instance, position, rel, ruv);
    } else {
        auto environment =
            lights.environments[light_id - (int)lights.instances.size()];
        return sample_environment(scene, lights, environment, rel, ruv);
    }
}

// Sample lights pdf
float sample_lights_pdf(const yocto_scene& scene, const trace_lights& lights,
    const bvh_scene& bvh, const vec3f& position, const vec3f& direction) {
    auto pdf = 0.0f;
    for (auto instance : lights.instances) {
        pdf += sample_light_pdf(
            scene, lights, instance, bvh, position, direction);
    }
    for (auto environment : lights.environments) {
        pdf += sample_environment_pdf(scene, lights, environment, direction);
    }
    pdf *= sample_uniform_pdf<float>(
        lights.instances.size() + lights.environments.size());
    return pdf;
}

// Russian roulette mode: 0: weight (pbrt), 1: albedo
#define YOCTO_RUSSIAN_ROULETTE_MODE 0
// Russian roulette on weight causes extremely high noise with prefectly
// transmissive materials because when reaching transparent surfaces (where you
// just want to continue path tracing) a path may still have low weight, hence
// high probability of beign terminated by Russian roulette.
// This causes to throw away almost every path that passes through multiple
// transparent surfaces: the path gets terminated very often, but it doesn't
// carry any radiance yet (since we don't split to trace lights) hence gives no
// contribution.
// The current fix is to sample Russian roulette on albedo, which is a local
// property and does not depend on the path. In this way we ensure that a path
// has always low propabilty of being terminated on transmissive surfaces.

// Russian roulette. Returns a weight and whether to stop.
#if YOCTO_RUSSIAN_ROULETTE_MODE == 0
pair<float, bool> sample_roulette(const vec3f& weight, int bounce, float rr,
    int min_bounce = 4, float rr_threadhold = 1) {
    if (max(weight) < rr_threadhold && bounce > min_bounce) {
        auto rr_prob = max((float)0.05, 1 - max(weight));
        return {1 - rr_prob, rr >= rr_prob};
    } else {
        return {1, false};
    }
}
#elif YOCTO_RUSSIAN_ROULETTE_MODE == 1
pair<float, bool> sample_roulette(const vec3f& albedo, const vec3f& weight,
    int bounce, float rr, int min_bounce = 4, float rr_threadhold = 1) {
    if (bounce <= min_bounce) return 1;
    auto rr_prob = clamp(1.0f - max(albedo), 0.0f, 0.95f);
    return rr < rr_prob ? 1 - rr_prob : 0;
}
#endif

pair<float, vec2i> sample_distance(
    const trace_mediums& mediums, float rl, float rd) {
    if (mediums.empty()) return {0, {-1, -1}};
    auto idx = sample_uniform((int)mediums.size(), rl);
    rl       = clamp(rl * (int)mediums.size() - idx, 0.0f, 1.0f);
    auto [distance, channel] = sample_distance(mediums[idx].density, rl, rd);
    return {distance, {idx, channel}};
}

float sample_distance_pdf(
    const trace_mediums& mediums, float distance, const vec2i& channel) {
    if (mediums.empty() || channel.x < 0) return 0;
    return sample_distance_pdf(mediums[channel.x].density, distance, channel.y);
}

vec3f eval_transmission(const trace_mediums& mediums, float distance) {
    auto transmission = vec3f{1, 1, 1};
    for (auto& lobe : mediums)
        transmission *= eval_transmission(lobe.density, distance);
    return transmission;
}

// Sample next direction. Returns weight and direction.
template <typename Lobe>
tuple<vec3f, vec3f> sample_direction(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh,
    const short_vector<Lobe, 8>& lobes, const vec3f& position,
    const vec3f& normal, const vec3f& outgoing, float rl, const vec2f& ruv) {
    // continue path
    auto incoming     = sample_scattering(lobes, normal, outgoing, rl, ruv);
    auto scattering   = eval_scattering(lobes, normal, outgoing, incoming);
    auto incoming_pdf = sample_scattering_pdf(
        lobes, normal, outgoing, incoming);
    if (incoming == zero3f || incoming_pdf == 0) {
        return {zero3f, zero3f};
    } else {
        return {scattering / incoming_pdf, incoming};
    }
}

template <typename Lobe>
tuple<vec3f, vec3f> sample_direction_mis(const yocto_scene& scene,
    const trace_lights& lights, const bvh_scene& bvh,
    const short_vector<Lobe, 8>& lobes, const vec3f& position,
    const vec3f& normal, const vec3f& outgoing, float rl, const vec2f& ruv,
    float re, float rmis) {
    // continue path
    auto incoming = (rmis < 0.5f)
                        ? sample_scattering(lobes, normal, outgoing, rl, ruv)
                        : sample_lights(
                              scene, lights, bvh, position, re, rl, ruv);
    auto scattering = eval_scattering(lobes, normal, outgoing, incoming);
    auto incoming_pdf =
        0.5f * sample_scattering_pdf(lobes, normal, outgoing, incoming) +
        0.5f * sample_lights_pdf(scene, lights, bvh, position, incoming);
    if (incoming == zero3f || incoming_pdf == 0) {
        return {zero3f, zero3f};
    } else {
        return {scattering / incoming_pdf, incoming};
    }
}

// Sample next mode. Returns weight and whether it is a delta.
pair<float, bool> sample_delta(
    const trace_bsdfs& bsdfs, const trace_deltas& deltas, float rn) {
    auto compute_weight = [](const auto& lobes) {
        auto weight = 0.0f;
        for (auto& lobe : lobes) weight += lobe.samplew;
        return weight;
    };
    auto bsdf_weight  = compute_weight(bsdfs);
    auto delta_weight = compute_weight(deltas);
    if (bsdf_weight == 0 && delta_weight == 0) return {0, false};
    auto delta_pdf = delta_weight / (bsdf_weight + delta_weight);
    if (rn < delta_pdf) {
        return {1 / delta_pdf, true};
    } else {
        return {1 / (1 - delta_pdf), false};
    }
}

// Returns weight and distance
pair<vec3f, float> sample_transmission(
    const trace_mediums& mediums, float max_distance, float rl, float rn) {
    if (mediums.empty()) return {vec3f{1}, max_distance};
    // clamp ray if inside a volume
    auto [distance, channel] = sample_distance(mediums, rl, rn);
    distance                 = min(distance, max_distance);
    auto pdf                 = sample_distance_pdf(mediums, distance, channel);
    auto transmission        = eval_transmission(mediums, distance);
    if (transmission == zero3f || pdf == 0) return {zero3f, 0};
    return {transmission / pdf, distance};
}

// Updates trace volume stack upon entering and exiting a surface
using trace_volume_stack = short_vector<pair<trace_mediums, int>, 32>;
void update_volume_stack(trace_volume_stack& volume_stack,
    trace_mediums& medium, const trace_mediums& material, int instance,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (material.empty()) return;
    if (dot(incoming, normal) > 0 != dot(outgoing, normal) > 0) {
        if (volume_stack.empty()) {
            volume_stack.push_back({material, instance});
            medium = material;
        } else {
            volume_stack.pop_back();
            medium = {};
            if (!volume_stack.empty()) medium = volume_stack.back().first;
        }
    }
}

// Trace point
struct trace_point {
    vec3f           position  = zero3f;
    vec3f           normal    = zero3f;
    trace_emissions emissions = {};
    trace_bsdfs     bsdfs     = {};
    trace_deltas    deltas    = {};
    trace_mediums   mediums   = {};
};

// Make a trace point
trace_point make_surface_point(const yocto_scene& scene,
    const bvh_intersection& intersection, const vec3f& shading_direction) {
    auto& instance = scene.instances[intersection.instance_id];
    auto& shape    = scene.shapes[instance.shape];
    auto& material = scene.materials[instance.material];
    auto  point    = trace_point{};
    auto& [position, normal, emissions, bsdfs, deltas, mediums] = point;
    position                                                    = eval_position(
        scene, instance, intersection.element_id, intersection.element_uv);
    normal         = eval_normal(scene, instance, intersection.element_id,
        intersection.element_uv, trace_non_rigid_frames);
    auto texcoords = eval_texcoord(
        shape, intersection.element_id, intersection.element_uv);
    auto color = eval_color(
        shape, intersection.element_id, intersection.element_uv);
    auto material_point = eval_material(scene, material, texcoords);
    eval_material(emissions, bsdfs, deltas, mediums, material_point, color,
        normal, -shading_direction);
    if (!shape.lines.empty()) {
        normal = orthonormalize(-shading_direction, normal);
    } else if (!shape.points.empty()) {
        normal = -shading_direction;
    } else {
        if (material.normal_texture >= 0) {
            normal = eval_perturbed_normal(scene, instance,
                intersection.element_id, intersection.element_uv,
                material_point.normal_map, trace_non_rigid_frames);
        }
    }
    return point;
}

trace_point make_volume_point(const trace_mediums& last_mediums,
    const vec3f& last_position, const vec3f& last_direction, float distance) {
    auto point                                                  = trace_point{};
    auto& [position, normal, emissions, bsdfs, deltas, mediums] = point;
    position = last_position + last_direction * distance;
    normal   = last_direction;
    mediums  = last_mediums;
    return point;
}

// Trace stats.
atomic<uint64_t> _trace_npaths{0};
atomic<uint64_t> _trace_nrays{0};

// Intersects a ray and returns a point
bool trace_ray(const yocto_scene& scene, const bvh_scene& bvh,
    const vec3f& position, const vec3f& direction,
    bvh_intersection& intersection) {
    _trace_nrays += 1;
    return intersect_bvh(
        scene, bvh, make_ray(position, direction), intersection);
}

// Sample camera
ray3f sample_camera(const yocto_camera& camera, const vec2i& ij,
    const vec2i& image_size, const vec2f& puv, const vec2f& luv) {
    return eval_camera(camera, ij, image_size, puv, sample_disk(luv));
}

// Recursive path tracing.
pair<vec3f, bool> trace_path(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& origin_, const vec3f& direction_,
    rng_state& rng, const trace_params& params) {
    // initialize
    auto radiance      = zero3f;
    auto weight        = vec3f{1, 1, 1};
    auto last_position = origin_;
    auto last_incoming = direction_;
    auto volume_stack  = trace_volume_stack{};
    auto last_mediums  = trace_mediums{};
    auto hit           = false;

    // trace  path
    for (auto bounce = 0; bounce < params.max_bounces; bounce++) {
        // intersect next point
        auto intersection = bvh_intersection{};
        if (!trace_ray(
                scene, bvh, last_position, last_incoming, intersection)) {
            radiance += weight * eval_environments(scene, last_incoming);
            break;
        }
        hit = true;

        // prepare to shade
        auto outgoing = -last_incoming;

        // clamp ray if inside a volume
        auto [transmission, distance] = sample_transmission(
            last_mediums, intersection.distance, rand1f(rng), rand1f(rng));
        weight *= transmission;
        auto on_surface = distance >= intersection.distance;

        // prepare shading point
        auto [position, normal, emissions, bsdfs, deltas, mediums] =
            on_surface ? make_surface_point(scene, intersection, last_incoming)
                       : make_volume_point(last_mediums, last_position,
                             last_incoming, distance);

        // accumulate emission
        radiance += on_surface
                        ? weight * eval_emission(emissions, normal, outgoing)
                        : weight * eval_emission(mediums, normal, outgoing);
        if (bsdfs.empty() && deltas.empty() && mediums.empty()) break;

        // russian roulette
        auto [rr_weight, rr_stop] = sample_roulette(
            weight, bounce, rand1f(rng));
        if (rr_stop) break;
        weight *= rr_weight;
        if (weight == zero3f) break;

        // pick delta
        auto [delta_weight, delta] =
            on_surface ? sample_delta(bsdfs, deltas, rand1f(rng))
                       : pair{1.0f, false};
        weight *= delta_weight;

        // next direction
        auto [scattering, incoming] =
            on_surface
                ? (delta ? sample_direction(scene, lights, bvh, deltas,
                               position, normal, outgoing, rand1f(rng),
                               rand2f(rng))
                         : sample_direction_mis(scene, lights, bvh, bsdfs,
                               position, normal, outgoing, rand1f(rng),
                               rand2f(rng), rand1f(rng), rand1f(rng)))
                : sample_direction_mis(scene, lights, bvh, mediums, position,
                      normal, outgoing, rand1f(rng), rand2f(rng), rand1f(rng),
                      rand1f(rng));
        weight *= scattering;
        if (weight == zero3f) break;

        // setup next iteration
        last_position = position;
        last_incoming = incoming;
        if (on_surface) {
            update_volume_stack(volume_stack, last_mediums, mediums,
                intersection.instance_id, normal, outgoing, incoming);
        }
    }

    return {radiance, hit};
}

// Recursive path tracing.
pair<vec3f, bool> trace_naive(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& origin_, const vec3f& direction_,
    rng_state& rng, const trace_params& params) {
    // initialize
    auto radiance      = zero3f;
    auto weight        = vec3f{1, 1, 1};
    auto last_position = origin_;
    auto last_incoming = direction_;
    auto hit           = false;

    // trace  path
    for (auto bounce = 0; bounce < params.max_bounces; bounce++) {
        // intersect next point
        auto intersection = bvh_intersection{};
        if (!trace_ray(
                scene, bvh, last_position, last_incoming, intersection)) {
            radiance += weight * eval_environments(scene, last_incoming);
            break;
        }
        hit = true;

        // prepare shading point
        auto outgoing = -last_incoming;
        auto [position, normal, emissions, bsdfs, deltas, mediums] =
            make_surface_point(scene, intersection, last_incoming);

        // accumulate emission
        radiance += weight * eval_emission(emissions, normal, outgoing);
        if (bsdfs.empty() && deltas.empty()) break;

        // russian roulette
        auto [rr_weight, rr_stop] = sample_roulette(
            weight, bounce, rand1f(rng));
        if (rr_stop) break;
        weight *= rr_weight;
        if (weight == zero3f) break;

        // pick delta
        auto [delta_weight, delta] = sample_delta(bsdfs, deltas, rand1f(rng));
        weight *= delta_weight;

        // next direction
        auto [scattering, incoming] =
            delta ? sample_direction(scene, lights, bvh, deltas, position,
                        normal, outgoing, rand1f(rng), rand2f(rng))
                  : sample_direction(scene, lights, bvh, bsdfs, position,
                        normal, outgoing, rand1f(rng), rand2f(rng));
        weight *= scattering;
        if (weight == zero3f) break;

        // setup next iteration
        last_position = position;
        last_incoming = incoming;
    }

    return {radiance, hit};
}

// Eyelight for quick previewing.
pair<vec3f, bool> trace_eyelight(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const vec3f& origin_, const vec3f& direction_,
    rng_state& rng, const trace_params& params) {
    // initialize
    auto radiance      = zero3f;
    auto weight        = vec3f{1, 1, 1};
    auto last_position = origin_;
    auto last_incoming = direction_;
    auto hit           = false;

    // trace  path
    for (auto bounce = 0; bounce < max(params.max_bounces, 4); bounce++) {
        // intersect next point
        auto intersection = bvh_intersection{};
        if (!trace_ray(
                scene, bvh, last_position, last_incoming, intersection)) {
            radiance += weight * eval_environments(scene, last_incoming);
            break;
        }
        hit = true;

        // prepare shading point
        auto outgoing = -last_incoming;
        auto [position, normal, emissions, bsdfs, deltas, mediums] =
            make_surface_point(scene, intersection, last_incoming);

        // accumulate emission
        radiance += weight * eval_emission(emissions, normal, outgoing);
        if (bsdfs.empty() && deltas.empty()) break;

        // brdf * light
        radiance += weight * pif *
                    eval_scattering(bsdfs, normal, outgoing, outgoing);

        // exit if needed
        if (weight == zero3f) break;

        // continue path
        auto [scattering, incoming] = sample_direction(scene, lights, bvh,
            deltas, position, normal, outgoing, rand1f(rng), zero2f);

        // exit if no hit
        weight *= scattering;
        if (weight == zero3f) break;

        // setup next iteration
        last_position = position;
        last_incoming = incoming;
    }

    return {radiance, hit};
}

// False color rendering
pair<vec3f, bool> trace_falsecolor(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& origin,
    const vec3f& direction, rng_state& rng, const trace_params& params) {
    // intersect next point
    auto intersection = bvh_intersection{};
    if (!trace_ray(scene, bvh, origin, direction, intersection)) {
        return {zero3f, false};
    }

    // get scene elements
    auto& instance = scene.instances[intersection.instance_id];
    auto& shape    = scene.shapes[instance.shape];
    // auto& material = scene.materials[instance.material];

    // prepare shading point
    auto outgoing = -direction;
    auto [position, normal, emissions, bsdfs, deltas, mediums] =
        make_surface_point(scene, intersection, direction);

    switch (params.falsecolor_type) {
        case trace_falsecolor_type::normal: {
            return {normal * 0.5f + 0.5f, 1};
        }
        case trace_falsecolor_type::frontfacing: {
            auto frontfacing = dot(normal, outgoing) > 0 ? vec3f{0, 1, 0}
                                                         : vec3f{1, 0, 0};
            return {frontfacing, 1};
        }
        case trace_falsecolor_type::gnormal: {
            auto normal = eval_element_normal(
                scene, instance, intersection.element_id, true);
            return {normal * 0.5f + 0.5f, 1};
        }
        case trace_falsecolor_type::gfrontfacing: {
            auto normal = eval_element_normal(
                scene, instance, intersection.element_id, true);
            auto frontfacing = dot(normal, outgoing) > 0 ? vec3f{0, 1, 0}
                                                         : vec3f{1, 0, 0};
            return {frontfacing, 1};
        }
        case trace_falsecolor_type::albedo: {
            auto weight = zero3f;
            for (auto& lobe : bsdfs) {
                weight += lobe.weight;
            }
            for (auto& lobe : deltas) {
                weight += lobe.weight;
            }
            return {weight, 1};
        }
        case trace_falsecolor_type::texcoord: {
            auto texturecoord = eval_texcoord(
                shape, intersection.element_id, intersection.element_uv);
            return {{texturecoord.x, texturecoord.y, 0}, 1};
        }
        case trace_falsecolor_type::color: {
            auto color = eval_color(
                shape, intersection.element_id, intersection.element_uv);
            return {xyz(color), 1};
        }
        case trace_falsecolor_type::emission: {
            auto weight = zero3f;
            for (auto& lobe : emissions) {
                weight += lobe.weight;
            }
            return {weight, 1};
        }
        case trace_falsecolor_type::diffuse: {
            auto weight = zero3f;
            for (auto& lobe : bsdfs) {
                if (lobe.type != trace_bsdf::type_t::diffuse) continue;
                weight += lobe.weight;
            }
            return {weight, 1};
        }
        case trace_falsecolor_type::specular: {
            auto weight = zero3f;
            for (auto& lobe : bsdfs) {
                if (lobe.type != trace_bsdf::type_t::reflection) continue;
                weight += lobe.weight;
            }
            for (auto& lobe : deltas) {
                if (lobe.type != trace_delta::type_t::reflection) continue;
                weight += lobe.weight;
            }
            return {weight, 1};
        }
        case trace_falsecolor_type::transmission: {
            auto weight = zero3f;
            for (auto& lobe : bsdfs) {
                if (lobe.type != trace_bsdf::type_t::transmission) continue;
                weight += lobe.weight;
            }
            for (auto& lobe : deltas) {
                if (lobe.type != trace_delta::type_t::transmission) continue;
                weight += lobe.weight;
            }
            return {weight, 1};
        }
        case trace_falsecolor_type::roughness: {
            auto roughness = zero3f;
            for (auto& lobe : bsdfs) {
                if (lobe.type != trace_bsdf::type_t::reflection) continue;
                roughness += lobe.roughness;
            }
            return {roughness, 1};
        }
        case trace_falsecolor_type::lobes: {
            return {vec3f{float(bsdfs.size() + deltas.size()) / 4}, 1};
        }
        case trace_falsecolor_type::material: {
            auto hashed = std::hash<int>()(instance.material);
            auto rng_   = make_rng(trace_default_seed, hashed);
            return {pow(0.5f + 0.5f * rand3f(rng_), 2.2f), 1};
        }
        case trace_falsecolor_type::shape: {
            auto hashed = std::hash<int>()(instance.shape);
            auto rng_   = make_rng(trace_default_seed, hashed);
            return {pow(0.5f + 0.5f * rand3f(rng_), 2.2f), 1};
        }
        case trace_falsecolor_type::instance: {
            auto hashed = std::hash<int>()(intersection.instance_id);
            auto rng_   = make_rng(trace_default_seed, hashed);
            return {pow(0.5f + 0.5f * rand3f(rng_), 2.2f), 1};
        }
        case trace_falsecolor_type::highlight: {
            auto emission = zero3f;
            for (auto& lobe : emissions) {
                emission += lobe.weight;
            }
            auto outgoing = -direction;
            if (emission == zero3f) emission = {0.2f, 0.2f, 0.2f};
            return {emission * abs(dot(outgoing, normal)), 1};
        }
        default: {
            return {zero3f, false};
        }
    }
}

// Trace a single ray from the camera using the given algorithm.
using trace_sampler_func = pair<vec3f, bool> (*)(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, const trace_params& params);
trace_sampler_func get_trace_sampler_func(const trace_params& params) {
    switch (params.sampler_type) {
        case trace_sampler_type::path: return trace_path;
        case trace_sampler_type::naive: return trace_naive;
        case trace_sampler_type::eyelight: return trace_eyelight;
        case trace_sampler_type::falsecolor: return trace_falsecolor;
        default: {
            throw runtime_error("sampler unknown");
            return nullptr;
        }
    }
}

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params) {
    switch (params.sampler_type) {
        case trace_sampler_type::path: return true;
        case trace_sampler_type::naive: return true;
        case trace_sampler_type::eyelight: return true;
        case trace_sampler_type::falsecolor: return true;
        default: {
            throw runtime_error("sampler unknown");
            return false;
        }
    }
}

// Get trace pixel
trace_pixel& get_trace_pixel(trace_state& state, int i, int j) {
    return state.pixels[j * state.image_size.x + i];
}

// Trace a block of samples
void trace_region(image<vec4f>& image, trace_state& state,
    const yocto_scene& scene, const bvh_scene& bvh, const trace_lights& lights,
    const image_region& region, int num_samples, const trace_params& params) {
    auto& camera  = scene.cameras.at(params.camera_id);
    auto  sampler = get_trace_sampler_func(params);
    for (auto j = region.min.y; j < region.max.y; j++) {
        for (auto i = region.min.x; i < region.max.x; i++) {
            auto& pixel = get_trace_pixel(state, i, j);
            for (auto s = 0; s < num_samples; s++) {
                if (params.cancel_flag && *params.cancel_flag) return;
                _trace_npaths += 1;
                auto ray = sample_camera(camera, {i, j}, image.size(),
                    rand2f(pixel.rng), rand2f(pixel.rng));
                auto [radiance, hit] = sampler(
                    scene, bvh, lights, ray.o, ray.d, pixel.rng, params);
                if (!hit) {
                    if (params.environments_hidden ||
                        scene.environments.empty()) {
                        radiance = zero3f;
                        hit      = false;
                    } else {
                        hit = true;
                    }
                }
                if (!isfinite(radiance)) {
                    // printf("NaN detected\n");
                    radiance = zero3f;
                }
                if (max(radiance) > params.pixel_clamp)
                    radiance = radiance * (params.pixel_clamp / max(radiance));
                pixel.radiance += radiance;
                pixel.hits += hit ? 1 : 0;
                pixel.samples += 1;
            }
            auto radiance = pixel.hits ? pixel.radiance / pixel.hits : zero3f;
            auto coverage = (float)pixel.hits / (float)pixel.samples;
            image[{i, j}] = {radiance.x, radiance.y, radiance.z, coverage};
        }
    }
}

// Init a sequence of random number generators.
void init_trace_state(
    trace_state& state, const vec2i& image_size, uint64_t seed) {
    state    = trace_state{image_size,
        vector<trace_pixel>(image_size.x * image_size.y, trace_pixel{})};
    auto rng = make_rng(1301081);
    for (auto j = 0; j < state.image_size.y; j++) {
        for (auto i = 0; i < state.image_size.x; i++) {
            auto& pixel = get_trace_pixel(state, i, j);
            pixel.rng   = make_rng(seed, rand1i(rng, 1 << 31) / 2 + 1);
        }
    }
}

// Init trace lights
void init_trace_lights(trace_lights& lights, const yocto_scene& scene) {
    lights = {};

    lights.shape_cdfs.resize(scene.shapes.size());
    lights.environment_cdfs.resize(scene.textures.size());

    for (auto instance_id = 0; instance_id < scene.instances.size();
         instance_id++) {
        auto& instance = scene.instances[instance_id];
        auto& shape    = scene.shapes[instance.shape];
        auto& material = scene.materials[instance.material];
        if (!material.emission_factor || material.emission_color == zero3f)
            continue;
        if (shape.triangles.empty() && shape.quads.empty()) continue;
        lights.instances.push_back(instance_id);
        sample_shape_cdf(shape, lights.shape_cdfs[instance.shape]);
    }

    for (auto environment_id = 0; environment_id < scene.environments.size();
         environment_id++) {
        auto& environment = scene.environments[environment_id];
        if (!environment.emission_factor ||
            environment.emission_color == zero3f)
            continue;
        lights.environments.push_back(environment_id);
        if (environment.emission_texture >= 0) {
            sample_environment_cdf(scene, environment,
                lights.environment_cdfs[environment.emission_texture]);
        }
    }
}

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image(const yocto_scene& scene, const bvh_scene& bvh,
    const trace_lights& lights, const trace_params& params) {
    auto image_size = camera_image_size(
        scene.cameras.at(params.camera_id), params.image_size);
    auto image = yocto::image{image_size, zero4f};
    auto state = trace_state{};
    init_trace_state(state, image_size, params.random_seed);
    auto regions = vector<image_region>{};
    make_regions(regions, image.size(), params.region_size, true);

    parallel_foreach(regions, [&image, &state, &scene, &bvh, &lights, &params](
                                  const image_region& region) {
        trace_region(image, state, scene, bvh, lights, region,
            params.num_samples, params);
    });

    return image;
}

// Progressively compute an image by calling trace_samples multiple times.
int trace_samples(image<vec4f>& image, trace_state& state,
    const yocto_scene& scene, const bvh_scene& bvh, const trace_lights& lights,
    int current_sample, const trace_params& params) {
    auto regions = vector<image_region>{};
    make_regions(regions, image.size(), params.region_size, true);
    auto num_samples = min(
        params.samples_per_batch, params.num_samples - current_sample);
    parallel_foreach(
        regions, [&image, &state, &scene, &bvh, &lights, num_samples, &params](
                     const image_region& region) {
            trace_region(
                image, state, scene, bvh, lights, region, num_samples, params);
        });
    return current_sample + num_samples;
}

// Starts an anyncrhounous renderer.
void trace_async_start(image<vec4f>& image, trace_state& state,
    const yocto_scene& scene, const bvh_scene& bvh, const trace_lights& lights,
    vector<future<void>>& futures, atomic<int>& current_sample,
    concurrent_queue<image_region>& queue, const trace_params& params) {
    auto& camera     = scene.cameras.at(params.camera_id);
    auto  image_size = camera_image_size(camera, params.image_size);
    image            = {image_size, zero4f};
    state            = trace_state{};
    init_trace_state(state, image_size, params.random_seed);
    auto regions = vector<image_region>{};
    make_regions(regions, image.size(), params.region_size, true);
    if (params.cancel_flag) *params.cancel_flag = false;

    futures.clear();
    futures.emplace_back(async([params, regions, &current_sample, &image,
                                   &scene, &lights, &bvh, &state, &queue]() {
        for (auto sample = 0; sample < params.num_samples;
             sample += params.samples_per_batch) {
            if (params.cancel_flag && *params.cancel_flag) return;
            current_sample   = sample;
            auto num_samples = min(
                params.samples_per_batch, params.num_samples - current_sample);
            parallel_foreach(
                regions,
                [num_samples, &params, &image, &scene, &lights, &bvh, &state,
                    &queue](const image_region& region) {
                    trace_region(image, state, scene, bvh, lights, region,
                        num_samples, params);
                    queue.push(region);
                },
                params.cancel_flag, params.run_serially);
        }
        current_sample = params.num_samples;
    }));
}

// Stop the asynchronous renderer.
void trace_async_stop(vector<future<void>>& futures,
    concurrent_queue<image_region>& queue, const trace_params& params) {
    if (params.cancel_flag) *params.cancel_flag = true;
    for (auto& f : futures) f.get();
    futures.clear();
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
// UNUSED CODE KEPT FOR REFERENCE
// -----------------------------------------------------------------------------
namespace yocto {

#if YOCTO_TRACE_THINSHEET

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(
    const vec3f& specular, const vec3f& half_vector, const vec3f& incoming) {
    if (specular == zero3f) return zero3f;
    return specular +
           (vec3f{1, 1, 1} - specular) *
               pow(clamp(1.0f - fabs(dot(half_vector, incoming)), 0.0f, 1.0f),
                   5.0f);
}
vec3f fresnel_schlick(const vec3f& specular, const vec3f& half_vector,
    const vec3f& incoming, float roughness) {
    if (specular == zero3f) return zero3f;
    auto fks = fresnel_schlick(specular, half_vector, incoming);
    return specular +
           (fks - specular) * (1 - sqrt(clamp(roughness, 0.0f, 1.0f)));
}

// Evaluates the GGX distribution and geometric term
float evaluate_ggx_distribution(
    float roughness, const vec3f& normal, const vec3f& half_vector) {
    auto di = (dot(normal, half_vector) * dot(normal, half_vector)) *
                  (roughness * roughness - 1) +
              1;
    return roughness * roughness / (pif * di * di);
}
float evaluate_ggx_shadowing(float roughness, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
#if 0
    // evaluate G from Heitz
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    auto Go = (2 * fabs(dot(normal, outgoing))) /
              (fabs(dot(normal, outgoing)) +
                  sqrt(roughness * roughness + (1 - roughness * roughness) *
                                                   dot(normal, outgoing) *
                                                   dot(normal, outgoing)));
    auto Gi = (2 * fabs(dot(normal, incoming))) /
              (fabs(dot(normal, incoming)) +
                  sqrt(roughness * roughness + (1 - roughness * roughness) *
                                                   dot(normal, incoming) *
                                                   dot(normal, incoming)));
    return Go * Gi;
#endif
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
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz   = sqrt(1 / (tan2 + 1));
    auto rr   = sqrt(1 - rz * rz);
    auto rphi = 2 * pif * rn.x;
    // set to wh
    auto wh_local = vec3f{rr * cos(rphi), rr * sin(rphi), rz};
    return wh_local;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f evaluate_brdf_cosine(const microfacet_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) return zero3f;
    auto brdf_cosine = zero3f;

    // diffuse
    if (brdf.diffuse != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto half_vector = normalize(incoming + outgoing);
        auto fresnel     = (brdf.fresnel) ? fresnel_schlick(brdf.specular,
                                            half_vector, outgoing)
                                      : brdf.specular;
        brdf_cosine += brdf.diffuse * (1 - fresnel) / pif;
    }

    // specular
    if (brdf.specular != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto half_vector = normalize(incoming + outgoing);
        auto fresnel     = (brdf.fresnel) ? fresnel_schlick(brdf.specular,
                                            half_vector, outgoing)
                                      : brdf.specular;
        auto D = evaluate_ggx_distribution(brdf.roughness, normal, half_vector);
        auto G = evaluate_ggx_shadowing(
            brdf.roughness, normal, outgoing, incoming);
        brdf_cosine +=
            fresnel * D * G /
            (4 * fabs(dot(normal, outgoing)) * fabs(dot(normal, incoming)));
    }

    // transmission (thin sheet)
    if (brdf.transmission != zero3f && !brdf.refract &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto ir = (dot(normal, outgoing) >= 0) ? reflect(-incoming, normal)
                                               : reflect(-incoming, -normal);
        auto half_vector = normalize(ir + outgoing);
        auto fresnel     = (brdf.fresnel) ? fresnel_schlick(brdf.specular,
                                            half_vector, outgoing)
                                      : brdf.specular;
        auto D = evaluate_ggx_distribution(brdf.roughness, normal, half_vector);
        auto G = evaluate_ggx_shadowing(brdf.roughness, normal, outgoing, ir);
        brdf_cosine +=
            brdf.transmission * (1 - fresnel) * D * G /
            (4 * fabs(dot(normal, outgoing)) * fabs(dot(normal, ir)));
    }

    // refraction through rough surface
    if (brdf.transmission != zero3f && brdf.refract &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto eta            = specular_to_eta(brdf.specular);
        auto halfway_vector = dot(normal, outgoing) > 0
                                  ? -(outgoing + eta * incoming)
                                  : (eta * outgoing + incoming);
        auto halfway = normalize(halfway_vector);

        auto fresnel = fresnel_schlick(brdf.specular, halfway, outgoing);
        auto D = evaluate_ggx_distribution(brdf.roughness, normal, halfway);
        auto G = evaluate_ggx_shadowing(
            brdf.roughness, normal, outgoing, incoming);

        auto dots = dot(outgoing, halfway) * dot(incoming, halfway) /
                    (dot(outgoing, normal) * dot(incoming, normal));

        auto numerator   = (1 - fresnel) * D * G;
        auto denominator = dot(halfway_vector, halfway_vector);

        brdf_cosine += abs(dots) * numerator / denominator;
    }
    return brdf_cosine * abs(dot(normal, incoming));
}

// Evaluates the BRDF assuming that it is called only from the directions
// generated by sample_brdf_direction.
vec3f evaluate_delta_brdf_cosine(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!is_brdf_delta(brdf)) return zero3f;
    auto microfacet_brdf = zero3f;

    // specular
    if (brdf.specular != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto fresnel = (brdf.fresnel)
                           ? fresnel_schlick(brdf.specular, normal, outgoing)
                           : brdf.specular;
        microfacet_brdf += fresnel;
    }

    // transmission (thin sheet)
    if (brdf.transmission != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto fresnel = (brdf.fresnel)
                           ? fresnel_schlick(brdf.specular, normal, outgoing)
                           : brdf.specular;
        microfacet_brdf += brdf.transmission * (1 - fresnel);
    }

    return microfacet_brdf;
}

// Picks a direction based on the BRDF
vec3f sample_brdf_direction(const microfacet_brdf& brdf, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
    if (is_brdf_delta(brdf)) return zero3f;
    auto fresnel = (brdf.fresnel)
                       ? fresnel_schlick(brdf.specular, normal, outgoing)
                       : brdf.specular;
    auto prob = vec3f{max(brdf.diffuse * (1 - fresnel)), max(fresnel),
        max(brdf.transmission * (1 - fresnel))};
    if (prob == zero3f) return zero3f;
    prob /= prob.x + prob.y + prob.z;

    // sample according to diffuse
    if (brdf.diffuse != zero3f && rnl < prob.x) {
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto il = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        auto fp = dot(normal, outgoing) >= 0
                      ? make_frame_fromz(zero3f, normal)
                      : make_frame_fromz(zero3f, -normal);
        return transform_direction(fp, il);
    }
    // sample according to specular GGX
    else if (brdf.specular != zero3f && rnl < prob.x + prob.y) {
        auto hl = sample_ggx(brdf.roughness, rn);
        auto fp = dot(normal, outgoing) >= 0
                      ? make_frame_fromz(zero3f, normal)
                      : make_frame_fromz(zero3f, -normal);
        auto half_vector = transform_direction(fp, hl);
        return reflect(outgoing, half_vector);
    }
    // transmission hack
    else if (brdf.transmission != zero3f && !brdf.refract &&
             rnl < prob.x + prob.y + prob.z) {
        auto hl = sample_ggx(brdf.roughness, rn);
        auto fp = dot(normal, outgoing) >= 0
                      ? make_frame_fromz(zero3f, normal)
                      : make_frame_fromz(zero3f, -normal);
        auto half_vector = transform_direction(fp, hl);
        auto ir          = reflect(outgoing, half_vector);
        return dot(normal, outgoing) >= 0 ? reflect(-ir, -normal)
                                          : reflect(-ir, normal);
    }
    // sample according to rough refraction
    else if (brdf.transmission != zero3f && brdf.refract &&
             rnl < prob.x + prob.y + prob.z) {
        auto hl = sample_ggx(brdf.roughness, rn);
        auto fp = dot(normal, outgoing) >= 0
                      ? make_frame_fromz(zero3f, normal)
                      : make_frame_fromz(zero3f, -normal);
        auto halfway = transform_direction(fp, hl);

        auto eta = specular_to_eta(brdf.specular);
        return refract(
            outgoing, halfway, dot(normal, outgoing) >= 0 ? 1 / eta : eta);
    }

    else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_delta_brdf_direction(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, float rnl, const vec2f& rn) {
    if (!is_brdf_delta(brdf)) return zero3f;
    auto fresnel = (brdf.fresnel)
                       ? fresnel_schlick(brdf.specular, normal, outgoing)
                       : brdf.specular;
    auto prob = vec3f{0, max(fresnel), max(brdf.transmission * (1 - fresnel))};
    if (prob == zero3f) return zero3f;
    prob /= prob.x + prob.y + prob.z;

    // sample according to specular mirror
    if (brdf.specular != zero3f && rnl < prob.x + prob.y) {
        return reflect(outgoing, dot(normal, outgoing) >= 0 ? normal : -normal);
    }
    // sample according to transmission
    else if (brdf.transmission != zero3f && !brdf.refract &&
             rnl < prob.x + prob.y + prob.z) {
        return -outgoing;
    }
    // sample according to transmission
    else if (brdf.transmission != zero3f && brdf.refract &&
             rnl < prob.x + prob.y + prob.z) {
        if (dot(normal, outgoing) >= 0) {
            return refract(
                outgoing, normal, 1 / specular_to_eta(brdf.specular));
        } else {
            return refract(outgoing, -normal, specular_to_eta(brdf.specular));
        }
    }
    // no sampling
    else {
        return zero3f;
    }
}

// Compute the weight for sampling the BRDF
float sample_brdf_direction_pdf(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (is_brdf_delta(brdf)) return 0;
    auto fresnel = (brdf.fresnel)
                       ? fresnel_schlick(brdf.specular, normal, outgoing)
                       : brdf.specular;
    auto prob = vec3f{max(brdf.diffuse * (1 - fresnel)), max(fresnel),
        max(brdf.transmission * (1 - fresnel))};
    if (prob == zero3f) return 0;
    prob /= prob.x + prob.y + prob.z;

    auto pdf = 0.0f;

    if (brdf.diffuse != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        pdf += prob.x * fabs(dot(normal, incoming)) / pif;
    }
    if (brdf.specular != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        auto half_vector = normalize(incoming + outgoing);
        auto d = sample_ggx_pdf(brdf.roughness, fabs(dot(normal, half_vector)));
        pdf += prob.y * d / (4 * fabs(dot(outgoing, half_vector)));
    }
    if (brdf.transmission != zero3f && !brdf.refract &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto ir = (dot(normal, outgoing) >= 0) ? reflect(-incoming, normal)
                                               : reflect(-incoming, -normal);
        auto half_vector = normalize(ir + outgoing);
        auto d = sample_ggx_pdf(brdf.roughness, fabs(dot(normal, half_vector)));
        pdf += prob.z * d / (4 * fabs(dot(outgoing, half_vector)));
    }

    // refraction through rough surface
    if (brdf.transmission != zero3f && brdf.refract &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        auto eta            = specular_to_eta(brdf.specular);
        auto outgoing_up    = dot(normal, outgoing) > 0;
        auto halfway_vector = outgoing_up ? -(outgoing + eta * incoming)
                                          : (eta * outgoing + incoming);
        auto halfway = normalize(halfway_vector);

        auto d = sample_ggx_pdf(brdf.roughness, fabs(dot(normal, halfway)));

        auto jacobian = fabs(dot(halfway, incoming)) /
                        dot(halfway_vector, halfway_vector);
        pdf += prob.z * d * jacobian;
    }

    return pdf;
}

// Compute the weight for sampling the BRDF
float sample_delta_brdf_direction_pdf(const microfacet_brdf& brdf,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
    if (!is_brdf_delta(brdf)) return 0;
    auto fresnel = (brdf.fresnel)
                       ? fresnel_schlick(brdf.specular, normal, outgoing)
                       : brdf.specular;
    auto prob = vec3f{0, max(fresnel), max(brdf.transmission * (1 - fresnel))};
    if (prob == zero3f) return 0;
    prob /= prob.x + prob.y + prob.z;

    auto pdf = 0.0f;

    if (brdf.specular != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) > 0) {
        return prob.y;
    }
    if (brdf.transmission != zero3f &&
        dot(normal, outgoing) * dot(normal, incoming) < 0) {
        return prob.z;
    }

    return pdf;
}

#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// NUMERICAL TESTS FOR MONTE CARLO INTEGRATION
// -----------------------------------------------------------------------------
namespace yocto {

template <typename Func>
float integrate_func_base(
    const Func& f, float a, float b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = rand1f(rng);
        auto x = a + r * (b - a);
        integral += f(x) * (b - a);
    }
    integral /= nsamples;
    return integral;
}

template <typename Func>
float integrate_func_stratified(
    const Func& f, float a, float b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = (i + rand1f(rng)) / nsamples;
        auto x = a + r * (b - a);
        integral += f(x) * (b - a);
    }
    integral /= nsamples;
    return integral;
}

template <typename Func>
float integrate_func_importance(const Func& f, const Func& pdf,
    const Func& warp, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = rand1f(rng);
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
template <typename Func>
void print_integrate_func_test(const Func& f, float a, float b, float expected,
    int nsamples, const Func& pdf, const Func& warp) {
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

template <typename Func>
float integrate_func2_base(
    const Func& f, vec2f a, vec2f b, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = rand2f(rng);
        auto x = a + r * (b - a);
        integral += f(x) * (b.x - a.x) * (b.y - a.y);
    }
    integral /= nsamples;
    return integral;
}

template <typename Func>
float integrate_func2_stratified(
    const Func& f, vec2f a, vec2f b, int nsamples, rng_state& rng) {
    auto integral  = 0.0f;
    auto nsamples2 = (int)sqrt(nsamples);
    for (auto i = 0; i < nsamples2; i++) {
        for (auto j = 0; j < nsamples2; j++) {
            auto r = vec2f{
                (i + rand1f(rng)) / nsamples2, (j + rand1f(rng)) / nsamples2};
            auto x = a + r * (b - a);
            integral += f(x) * (b.x - a.x) * (b.y - a.y);
        }
    }
    integral /= nsamples2 * nsamples2;
    return integral;
}

template <typename Func, typename Func2>
float integrate_func2_importance(const Func& f, const Func& pdf,
    const Func2& warp, int nsamples, rng_state& rng) {
    auto integral = 0.0f;
    for (auto i = 0; i < nsamples; i++) {
        auto r = rand2f(rng);
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
template <typename Func, typename Func2>
void print_integrate_func2_test(const Func& f, vec2f a, vec2f b, float expected,
    int nsamples, const Func& pdf, const Func2& warp) {
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
