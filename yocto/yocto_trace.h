//
// YOCTO_TRACE: path tracer implementation for support for textured mesh
// lights, GGX/Phong materials, environment mapping. The interface supports
// progressive parallel execution with any splitting strategy, by
// generating per-sample random number sequences or fully deterministic
// hash-based sampling.
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. define your scene setting up only the public data
// - init the scene
//     scene = ytrace::scene()
// - define shapes, materials and textures (we store pointers, not copies)
//     foreach shape: scene.shapes.push_back({shape data})
//     foreach material: scene.materials.push_back({mat data})
//     foreach texture: scene.textures.push_back({texture data})
// - define cameras and environments
//     foreach camera: scene.camera.push_back({cam data})
//     foreach environments: scene.environments.push_back({env data})
// - intersection routines are handled outside of this library using callbacks
//     scn.intersect_first = <callback>
//     scn.intersect_any = <callback>
//     - can use yocto_bvh
// 2. prepare for rendering
//    init_lights(scene)
// 3. define rendering params
// 4. render blocks of samples
//    render_block(scn, pixels, image size, block)
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by array of vertex indices for
// points, lines and triangles, and arrays of vertex data. Only one primitive
// type can be non-empty for each shape.
//
// Materials are represented as sums of an emission term, a diffuse term and
// a specular microfacet term (GGX or Phong). Only opaque for now. We pick
// a proper material type for each shape element type (points, lines,
// triangles).
//
// Lights are defined as any shape with a material emission term. Additionally
// one can also add an environment map. But even if you can, you might want to
// add a large triangle mesh with inward normals instead. The latter is more
// general (you can even more an arbitrary shape sun). For now only the first
// env is used.
//
// We generate our own random numbers guarantying that there is one random
// sequence per path. This means you can rul the path tracer in any order
// serially or in parallel.
//
// For now, we support a straightforward path tracer with explicit direct
// illumination using MIS.
//

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use in C++.
// To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.6: minor API change for blocks
// - v 0.5: [major API change] move to modern C++ interface
// - v 0.4: C++ API
// - v 0.3: removal of C interface
// - v 0.2: use of STL containers
// - v 0.1: C++ implementation
// - v 0.0: initial release in C99
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YTRACE_H_
#define _YTRACE_H_

#ifndef YGL_DECLARATION
#define YGL_API inline
#else
#define YGL_API
#endif

#include <functional>
#include "yocto_math.h"

// -----------------------------------------------------------------------------
// C++ INTERFACE
// -----------------------------------------------------------------------------

namespace ytrace {

//
// Using directives
//
using namespace ym;

//
// Type of rendering algorithm (shader)
//
enum struct shader_type {
    def = 0,    // default renderer
    eyelight,   // eye hight for quick previews
    direct,     // direct illumination
    pathtrace,  // path tracing
};

//
// Random number generator type
//
enum struct rng_type {
    def = 0,     // default generator
    uniform,     // uniform random numbers
    stratified,  // stratified random numbers
    cmjs,        // correlated multi-jittered sampling
};

//
// Camera
//
struct camera {
    frame3f xform = identity_frame3f;  // local-to-world transform
    float yfov = pif / 3;              // field of view
    float aspect = 1;                  // aspect ratio
    float aperture = 0;                // lens aperture
    float focus = 1;                   // focus plane distance (cannot be zero)
};

//
// Texture
//
struct texture {
    image_view<vec4f> hdr;
    image_view<vec4b> ldr;
};

//
// Material
//
struct material {
    // material values
    vec3f ke = zero3f;  // emission, term
    vec3f kd = zero3f;  // diffuse term
    vec3f ks = zero3f;  // specular term
    float rs = 0.1;     // specular roughness API

    // fresnel
    vec3f es = zero3f;   // eta
    vec3f eks = zero3f;  // etak (metals only)

    // textures
    int ke_txt = -1;
    int kd_txt = -1;
    int ks_txt = -1;
    int rs_txt = -1;

    // material flags
    bool use_phong = false;  // whether to use phong
};

//
// Shape
//
struct shape {
    frame3f xform = identity_frame3f;  // local-to-world rigid transform
    int matid = -1;                    // material id

    // element data [only one enabled at any given time]
    array_view<vec1i> points;     // elem data
    array_view<vec2i> lines;      // elem data
    array_view<vec3i> triangles;  // elem data

    // vertex data
    array_view<vec3f> pos;       // vertex data
    array_view<vec3f> norm;      // vertex data
    array_view<vec2f> texcoord;  // vertex data
    array_view<vec3f> color;     // vertex data
    array_view<float> radius;    // vertex data
};

//
// Environment
//
struct environment {
    frame3f xform = identity_frame3f;  // local-to-world rigid transform
    vec3f ke = zero3f;                 // emission
    int ke_txt = -1;                   // emission texture
};

//
// Ray-scn Intersection.
//
struct intersect_point {
    float dist = 0;      // ray distance
    int sid = -1;        // shape index
    int eid = -1;        // element index
    vec3f euv = zero3f;  // element baricentric coordinates

    // check whether it was a hit
    operator bool() const { return eid >= 0; }
};

//
// Ray-scn closest intersection callback
//
// Parameters:
// - ray: ray
//
// Return:
// - intersection point
//
using intersect_first_cb = function<intersect_point(const ray3f& ray)>;

//
// Ray-scn intersection callback
//
// Parameters:
// - ray: ray
//
// Return:
// - whether we intersect or not
//
using intersect_any_cb = function<bool(const ray3f& ray)>;

//
// Light (either shape or environment).
// This is only used internally and should not be created.
//
struct light {
    int shape_id = -1;  // shape
    int env_id = -1;    // environment
    vector<float> cdf;  // for shape, cdf of shape elements for sampling
    float area = 0;     // for shape, shape area
};

//
// Scene
//
struct scene {
    // intersection callbaks
    intersect_first_cb intersect_first;  // ray intersection callback
    intersect_any_cb intersect_any;      // ray hit callback

    // scn data
    vector<camera> cameras;            // camera
    vector<environment> environments;  // env
    vector<shape> shapes;              // shapes
    vector<material> materials;        // materials
    vector<texture> textures;          // textures

    // [private] light sources
    vector<light> _lights;  // lights [private]
};

//
// Rendering params
//
struct render_params {
    shader_type stype = shader_type::def;  // smp type
    rng_type rtype = rng_type::def;        // random type
    vec3f amb = zero3f;                    // ambient lighting
    int min_depth = 3;                     // min ray depth
    int max_depth = 8;                     // mas ray depth
    float pixel_clamp = 100;               // final pixel clamping
    float ray_eps = 1e-2f;                 // ray intersection epsilon
};

//
// Convert a Phong exponent to GGX/Phong roughness
//
YGL_API float specular_exponent_to_roughness(float n);

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YGL_API void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);

//
// Initialize rendering.
//
// Parameters:
// - scn: trace scene
//
YGL_API void init_lights(scene& scn);

//
// Renders a block of sample
//
// Parameters:
// - scn: trace scene
// - cid: camera id
// - img: pixel data in RGBA format
// - ns: number of samples
// - block: image block to render [xmin, xmax, ymin, ymax];
// max values are excluded
// - sampples: sample block to render [sample_min, sample_max];
// max values are excluded
// - accumulate: whether to accumulate the results of subsequent calls.
//
// Notes: It is safe to call the function in parallel one different blocks.
// But two threads should not access the same pixels at the same time.
// Also blocks with different samples should be called sequentially if
// accumulate is true.
//
YGL_API void trace_block(const scene& scn, int cid, image_view<vec4f> img,
                         int ns, const vec2i& xy, const vec2i& wh,
                         const vec2i& samples, const render_params& params,
                         bool accumulate = false);

//
// Convenience function to call trace_block with all sample at once.
//
YGL_API void trace_image(const scene& scn, int cid, image_view<vec4f> img,
                         int ns, const render_params& params);

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

#include <cassert>
#include <cfloat>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace ytrace {

//
// Compute shape element cdf for shape sampling.
//
template <typename T, typename Weight_callback>
static inline vector<float> _compute_weight_cdf(
    const array_view<T>& elem, float& total_weight,
    const Weight_callback& weight_cb) {
    // prepare return
    auto cdf = vector<float>();
    cdf.reserve(elem.size());

    // compute weights
    total_weight = 0;
    for (auto&& e : elem) {
        auto w = weight_cb(e);
        total_weight += w;
        cdf.push_back(total_weight);
    }

    // normalize
    for (auto& c : cdf) c /= total_weight;

    // done
    return cdf;
}

//
// Init lights. Public API, see above.
//
YGL_API void init_lights(scene& scn) {
    // clear old lights
    scn._lights.resize(0);

    for (int sid = 0; sid < scn.shapes.size(); sid++) {
        auto& shp = scn.shapes[sid];
        auto& mat = scn.materials[shp.matid];
        if (mat.ke == zero3f) continue;
        scn._lights.push_back(light());
        auto& light = scn._lights.back();
        light.shape_id = sid;
        if (!shp.points.empty()) {
            light.cdf = _compute_weight_cdf(shp.points, light.area,
                                            [&shp](auto e) { return 1; });
        } else if (!shp.lines.empty()) {
            light.cdf =
                _compute_weight_cdf(shp.lines, light.area, [&shp](auto e) {
                    return length(shp.pos[e[1]] - shp.pos[e[0]]);
                });
        } else if (!shp.triangles.empty()) {
            light.cdf =
                _compute_weight_cdf(shp.triangles, light.area, [&shp](auto e) {
                    return triangle_area(shp.pos[e[0]], shp.pos[e[1]],
                                         shp.pos[e[2]]);
                });
        }
    }

    for (int envid = 0; envid < scn.environments.size(); envid++) {
        auto& env = scn.environments[envid];
        if (env.ke == zero3f) continue;
        scn._lights.push_back(light());
        auto& light = scn._lights.back();
        light.env_id = envid;
    }
}

//
// Phong exponent to roughness. Public API, see above.
//
YGL_API float specular_exponent_to_roughness(float n) {
    return sqrtf(2 / (n + 2));
}

//
// Specular to fresnel eta. Public API, see above.
//
YGL_API void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es = {(1 + sqrt(ks[0])) / (1 - sqrt(ks[0])),
          (1 + sqrt(ks[1])) / (1 - sqrt(ks[1])),
          (1 + sqrt(ks[2])) / (1 - sqrt(ks[2]))};
    esk = zero3f;
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// Random number smp. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
//
struct _sampler {
    rng_pcg32 rng;   // rnumber number state
    int i, j;        // pixel coordinates
    int s, d;        // sample and dimension indices
    int ns;          // number of samples
    rng_type rtype;  // random number type
};

//
// Initialize a smp ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
//
static inline _sampler _make_sampler(int i, int j, int s, int ns,
                                     rng_type rtype) {
    // we use various hashes to scramble the pixel values
    _sampler smp = {{0, 0}, i, j, s, 0, ns, rtype};
    uint64_t sample_id = ((uint64_t)(i + 1)) << 0 | ((uint64_t)(j + 1)) << 15 |
                         ((uint64_t)(s + 1)) << 30;
    uint64_t initseq = hash_uint64(sample_id);
    uint64_t initstate = hash_uint64(sample_id * 3202034522624059733ull + 1ull);
    rng_init(smp.rng, initstate, initseq);
    return smp;
}

//
// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
//
static inline float _sample_next1f(_sampler& smp) {
    float rn = 0;
    switch (smp.rtype) {
        case rng_type::def:
        case rng_type::uniform: {
            rn = rng_nextf(smp.rng);
        } break;
        case rng_type::stratified: {
            uint32_t p = hash_uint64_32(((uint64_t)(smp.i + 1)) << 0 |
                                        ((uint64_t)(smp.j + 1)) << 15 |
                                        ((uint64_t)(smp.d + 1)) << 30);
            int s = hash_permute(smp.s, smp.ns, p);
            rn = (s + rng_nextf(smp.rng)) / smp.ns;
        } break;
        case rng_type::cmjs: {
            uint32_t p = hash_uint64_32(((uint64_t)(smp.i + 1)) << 0 |
                                        ((uint64_t)(smp.j + 1)) << 15 |
                                        ((uint64_t)(smp.d + 1)) << 30);
            int s = hash_permute(smp.s, smp.ns, p);
            rn = (s + hash_randfloat(s, p * 0xa399d265)) / smp.ns;
        } break;
        default: assert(false);
    }

    smp.d += 1;

    // make sure all sampled numbers are below 1
    // TODO: use numeric_limits
    if (rn >= 1) rn = 1 - FLT_EPSILON;

    return rn;
}

//
// Generates a 1-dimensional sample.
//
// Implementation notes: see above. Note that using deterministic keyed
// permutaton we can use stratified sampling without preallocating samples.
//
static inline vec2f _sample_next2f(_sampler& smp) {
    vec2f rn = {0, 0};
    switch (smp.rtype) {
        case rng_type::def:
        case rng_type::uniform: {
            rn[0] = rng_nextf(smp.rng);
            rn[1] = rng_nextf(smp.rng);
        } break;
        case rng_type::stratified: {
            uint32_t ns2 = (uint32_t)round(sqrt(smp.ns));
            uint32_t p = hash_uint64_32(((uint64_t)(smp.i + 1)) << 0 |
                                        ((uint64_t)(smp.j + 1)) << 15 |
                                        ((uint64_t)(smp.d + 1)) << 30);
            int s = hash_permute(smp.s, smp.ns, p);
            rn[0] = (s % ns2 + rng_nextf(smp.rng)) / ns2;
            rn[1] = (s / ns2 + rng_nextf(smp.rng)) / ns2;
        } break;
        case rng_type::cmjs: {
            uint32_t ns2 = (uint32_t)round(sqrt(smp.ns));
            uint32_t p = hash_uint64_32(((uint64_t)(smp.i + 1)) << 0 |
                                        ((uint64_t)(smp.j + 1)) << 15 |
                                        ((uint64_t)(smp.d + 1)) << 30);
            int s = hash_permute(smp.s, smp.ns, p);
            int sx = hash_permute(s % ns2, ns2, p * 0xa511e9b3);
            int sy = hash_permute(s / ns2, ns2, p * 0x63d83595);
            float jx = hash_randfloat(s, p * 0xa399d265);
            float jy = hash_randfloat(s, p * 0x711ad6a5);
            rn[0] = (s % ns2 + (sy + jx) / ns2) / ns2;
            rn[1] = (s / ns2 + (sx + jy) / ns2) / ns2;
        } break;
        default: assert(false);
    }

    smp.d += 2;

    // make sure all sampled numbers are below 1
    if (rn[0] >= 1) rn[0] = 1 - FLT_EPSILON;
    if (rn[1] >= 1) rn[1] = 1 - FLT_EPSILON;

    return rn;
}

//
// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
//
struct point {
    // point type -----------------------------------
    enum struct type {
        none = -1,     // invalid
        env = 0,       // environment
        point = 1,     // points
        line = 2,      // lines
        triangle = 3,  // triangle
    };
    type ptype = type::none;  // element type

    // light id -----------------------------
    int light_id = -1;  // light id used for MIS

    // direction ----------------------------
    vec3f wo = zero3f;  // outgoing direction

    // resolved geometry (shape) ------------
    frame3f frame = identity_frame3f;  // local frame

    // shading ------------------------------
    vec3f ke = zero3f;       // material values
    vec3f kd = zero3f;       // material values
    vec3f ks = zero3f;       // material values
    float rs = 0;            // material values
    vec3f es = zero3f;       // material values
    vec3f eks = zero3f;      // material values
    bool use_phong = false;  // material values
};

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static inline ray3f _eval_camera(const camera& cam, const vec2f& uv,
                                 const vec2f& luv) {
    auto h = 2 * tan(cam.yfov / 2);
    auto w = h * cam.aspect;
    vec3f o = vec3f{luv[0] * cam.aperture, luv[1] * cam.aperture, 0};
    vec3f q = {w * cam.focus * (uv[0] - 0.5f), h * cam.focus * (uv[1] - 0.5f),
               -cam.focus};
    return ray3f(transform_point(cam.xform, o),
                 transform_direction(cam.xform, normalize(q - o)));
}

//
// Wrapper for above function
//
static inline vec4f _eval_texture(const texture& txt, const vec2f& texcoord) {
    assert(!txt.hdr.empty() || !txt.ldr.empty());

    // get image width/height
    auto wh = (!txt.ldr.empty()) ? txt.ldr.size() : txt.hdr.size();

    // get coordinates normalized for tiling
    auto st = vec2f{fmod(texcoord[0], 1.0f), fmod(texcoord[1], 1.0f)} *
              vec2f(wh[0], wh[1]);
    if (st[0] < 0) st[0] += wh[0];
    if (st[1] < 0) st[1] += wh[1];

    // get image coordinates and residuals
    auto ij = clamp(vec2i(st[0], st[1]), {0, 0}, wh);
    auto uv = st - vec2f(ij[0], ij[1]);

    // get interpolation weights and indices
    vec2i idx[4] = {ij,
                    {ij[0], (ij[1] + 1) % wh[1]},
                    {(ij[0] + 1) % wh[0], ij[1]},
                    {(ij[0] + 1) % wh[0], (ij[1] + 1) % wh[1]}};
    float w[4] = {(1 - uv[0]) * (1 - uv[1]), (1 - uv[0]) * uv[1],
                  uv[0] * (1 - uv[1]), uv[0] * uv[1]};

    // handle interpolation
    if (!txt.ldr.empty()) {
        return (srgb_to_linear(txt.ldr[idx[0]]) * w[0] +
                srgb_to_linear(txt.ldr[idx[1]]) * w[1] +
                srgb_to_linear(txt.ldr[idx[2]]) * w[2] +
                srgb_to_linear(txt.ldr[idx[3]]) * w[3]);
    } else if (!txt.hdr.empty()) {
        return (txt.hdr[idx[0]] * w[0] + txt.hdr[idx[1]] * w[1] +
                txt.hdr[idx[2]] * w[2] + txt.hdr[idx[3]] * w[3]);
    } else {
        assert(false);
    }

    // should not have gotten here
    return zero4f;
}

//
// Evaluates emission.
//
static inline vec3f _eval_emission(const point& pt) {
    if (pt.ke == zero3f) return zero3f;
    switch (pt.ptype) {
        case point::type::env: return pt.ke;
        case point::type::point: return pt.ke;
        case point::type::line: return pt.ke;
        case point::type::triangle:
            return (dot(pt.frame[2], pt.wo) > 0) ? pt.ke : zero3f;
        default: {
            assert(false);
            return zero3f;
        }
    }
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline vec3f _eval_fresnel_dielectric(float cosw, const vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1 / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = 1 - sin2 / eta2;
    if (cos2t[0] < 0 || cos2t[1] < 0 || cos2t[2] < 0)
        return vec3f{1, 1, 1};  // tir

    auto t0 = vec3f{sqrt(cos2t[0]), sqrt(cos2t[1]), sqrt(cos2t[2])};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (cosw - t1) / (cosw + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2;
}

//
// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline vec3f _eval_fresnel_metal(float cosw, const vec3f& eta,
                                        const vec3f& etak) {
    if (etak == zero3f) return _eval_fresnel_dielectric(cosw, eta);

    cosw = clamp(cosw, (float)-1, (float)1);
    auto cos2 = cosw * cosw;
    auto sin2 = clamp(1 - cos2, (float)0, (float)1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - sin2;
    auto a2plusb2_2 = t0 * t0 + 4 * eta2 * etak2;
    auto a2plusb2 =
        vec3f{sqrt(a2plusb2_2[0]), sqrt(a2plusb2_2[1]), sqrt(a2plusb2_2[2])};
    auto t1 = a2plusb2 + cos2;
    auto a_2 = (a2plusb2 + t0) / 2;
    auto a = vec3f{sqrt(a_2[0]), sqrt(a_2[1]), sqrt(a_2[2])};
    auto t2 = 2 * a * cosw;
    auto rs = (t1 - t2) / (t1 + t2);

    auto t3 = cos2 * a2plusb2 + sin2 * sin2;
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2;
}

//
// Evaluates the BRDF scaled by the cosine of the incoming direction.
//
// Implementation notes:
// - ggx from [Heitz 2014] and [Walter 2007]
// "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
// http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
// - uses Kajiya-Kay for hair
// - uses a hack for points
//
static inline vec3f _eval_brdfcos(const point& pt, const vec3f& wi) {
    // summation over multiple terms
    auto brdfcos = zero3f;

    // exit if not needed
    if (pt.kd == zero3f && pt.ks == zero3f) return brdfcos;

    // save wo
    auto wo = pt.wo;

    // compute wh
    auto wh = normalize(wo + wi);

    // compute dot products
    auto ndo = dot(pt.frame[2], wo), ndi = dot(pt.frame[2], wi),
         ndh = clamp(dot(wh, pt.frame[2]), (float)0, (float)1);

    switch (pt.ptype) {
        case point::type::point: {
            // diffuse term (hack for now)
            if (pt.kd != zero3f) {
                auto ido = dot(wo, wi);
                auto diff = pt.kd * (2 * ido + 1) / (2 * pif);
                brdfcos += diff;
            }
        } break;
        case point::type::line: {
            // take sines
            auto so = sqrt(clamp(1 - ndo * ndo, (float)0, (float)1)),
                 si = sqrt(clamp(1 - ndi * ndi, (float)0, (float)1)),
                 sh = sqrt(clamp(1 - ndh * ndh, (float)0, (float)1));

            // diffuse term (Kajiya-Kay)
            if (si > 0 && so > 0 && pt.kd != zero3f) {
                auto diff = pt.kd * si / pif;
                brdfcos += diff;
            }

            // specular term (Kajiya-Kay)
            if (si > 0 && so > 0 && sh > 0 && pt.ks != zero3f) {
                auto ns = 2 / (pt.rs * pt.rs) - 2;
                auto d = (ns + 2) * pow(sh, ns) / (2 + pif);
                auto spec = pt.ks * si * d / (4 * si * so);
                brdfcos += spec;
            }
        } break;
        case point::type::triangle: {
            // diffuse term
            if (ndi > 0 && ndo && pt.kd != zero3f) {
                auto diff = pt.kd * ndi / pif;
                brdfcos += diff;
            }

            // specular term (GGX)
            if (ndi > 0 && ndo > 0 && ndh > 0 && pt.ks != zero3f) {
                if (!pt.use_phong) {
                    // evaluate GGX
                    auto cos2 = ndh * ndh;
                    auto tan2 = (1 - cos2) / cos2;
                    auto alpha2 = pt.rs * pt.rs;
                    auto d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) *
                                       (alpha2 + tan2));
                    auto lambda_o =
                        (-1 + sqrtf(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    auto lambda_i =
                        (-1 + sqrtf(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    auto g = 1 / (1 + lambda_o + lambda_i);
                    auto spec = pt.ks * ndi * d * g / (4 * ndi * ndo);
                    if (pt.es != zero3f) {
                        spec *= _eval_fresnel_metal(ndh, pt.es, pt.eks);
                    }
                    brdfcos += spec;
                } else {
                    // evaluate Blinn-Phong
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    auto d = (ns + 2) * powf(ndh, ns) / (2 + pif);
                    auto spec = pt.ks * ndi * d / (4 * ndi * ndo);
                    brdfcos += spec;
                }
            }
        } break;
        default: { assert(false); }
    }

    return brdfcos;
}

//
// Compute the weight for sampling the BRDF
//
static inline float _weight_brdfcos(const point& pt, const vec3f& wi) {
    // skip if no component
    if (pt.kd == zero3f && pt.ks == zero3f) return 0;

    // save wo
    auto wo = pt.wo;

    // compute wh
    auto wh = normalize(wi + wo);

    // compute dot products
    auto ndo = dot(pt.frame[2], wo), ndi = dot(pt.frame[2], wi),
         ndh = dot(pt.frame[2], wh);

    // check to make sure we are above the surface
    // updated this for refraction
    if (ndo <= 0 || ndi <= 0) return 0;

    // pick from a sum
    auto wall = mean(pt.kd) + mean(pt.ks);
    auto wd = mean(pt.kd) / wall;
    auto ws = mean(pt.ks) / wall;

    // accumulate probability
    auto pdf = 0.0f;

    switch (pt.ptype) {
        case point::type::point: {
        } break;
        case point::type::line: {
            // diffuse term
            if (wall) {
                // homepherical cosine probability
                pdf += 1 / (4 * pif);
            }
        } break;
        case point::type::triangle: {
            // diffuse term
            if (wd && ndi > 0) {
                // homepherical cosine probability
                pdf += wd * ndi / pif;
            }

            // specular term (GGX or Phong)
            if (ws && ndi > 0 && ndo > 0 && ndh > 0) {
                if (!pt.use_phong) {
                    // probability proportional to d * ndh
                    auto cos2 = ndh * ndh;
                    auto tan2 = (1 - cos2) / cos2;
                    auto alpha2 = pt.rs * pt.rs;
                    auto d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) *
                                       (alpha2 + tan2));
                    auto hdo = dot(wo, wh);
                    pdf += ws * d * ndh / (4 * hdo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // compute wh
                    auto wh = normalize(wi + wo);
                    auto ndh = dot(pt.frame[2], wh);
                    // homerispherical cosine power probability
                    pdf += ws * powf(ndh, ns) * (ns + 1) / (2 * pif);
                }
            }
        } break;
        default: { assert(false); }
    }

    // done
    return (pdf) ? 1 / pdf : 0;
}

//
// Picks a direction based on the BRDF
//
static inline vec3f _sample_brdfcos(const point& pt, float rnl,
                                    const vec2f& rn) {
    // skip if no component
    if (pt.kd == zero3f && pt.ks == zero3f) return zero3f;

    // save wo
    auto wo = pt.wo;

    // compute cosine
    auto ndo = dot(pt.frame[2], wo);

    // check to make sure we are above the surface
    // update this for refraction
    if (ndo <= 0) return zero3f;

    // pick from a sum
    auto wall = mean(pt.kd) + mean(pt.ks);
    auto wd = mean(pt.kd) / wall;
    auto ws = mean(pt.ks) / wall;

    switch (pt.ptype) {
        // TODO: point color
        case point::type::point:
        case point::type::line: {
            if (wall > 0) {
                // sample wi with uniform spherical distribution
                auto rz = rn[1], rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn[0];
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(pt.frame, wi_local);
            }
        } break;
        case point::type::triangle: {
            // sample according to diffuse
            if (rnl < wd) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn[1]), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn[0];
                // set to wi
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return transform_direction(pt.frame, wi_local);
            }

            // sample according to specular (GGX or Phong)
            if (rnl >= wd && rnl < wd + ws) {
                if (!pt.use_phong) {
                    // sample wh with ggx distribution
                    auto tan2 = pt.rs * pt.rs * rn[1] / (1 - rn[1]);
                    auto rz = sqrtf(1 / (tan2 + 1)), rr = sqrtf(1 - rz * rz),
                         rphi = 2 * pif * rn[0];
                    // set to wh
                    auto wh_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                    auto wh = transform_direction(pt.frame, wh_local);
                    // compute wi
                    return normalize(wh * 2 * dot(wo, wh) - wo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // sample wh with hemispherical cosine power distribution
                    auto rz = powf(rn[1], 1 / (ns + 1)),
                         rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn[0];
                    // set to wh
                    auto wh_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                    auto wh = transform_direction(pt.frame, wh_local);
                    // compute wi
                    return normalize(wh * 2 * dot(wo, wh) - wo);
                }
            }
        } break;
        default: { assert(false); }
    }

    // should not have gotten here
    assert(false);
    return zero3f;
}

//
// Create a point for an environment map. Resolves material with textures.
//
static inline point _eval_envpoint(const scene& scn, int env_id,
                                   const vec3f& wo) {
    // set shape data
    auto pt = point();

    // check if null point
    if (env_id < 0) return pt;
    auto& env = scn.environments[env_id];

    // env params
    pt.ptype = point::type::env;

    // direction
    pt.wo = wo;

    // maerial
    pt.ke = env.ke;

    // textures
    if (env.ke_txt >= 0) {
        auto w = transform_direction(inverse(env.xform), -wo);
        auto theta = 1 - (acos(clamp(w[1], (float)-1, (float)1)) / pif);
        auto phi = atan2f(w[2], w[0]) / (2 * pif);
        auto texcoord = vec2f{phi, theta};
        if (env.ke_txt >= 0) {
            auto txt = _eval_texture(scn.textures[env.ke_txt], texcoord);
            pt.ke = lerp({txt[0], txt[1], txt[2]}, pt.ke, txt[3]);
        }
    }

    // done
    return pt;
}

//
// Interpolate a value over an element
//
template <typename T, size_t N>
static inline T _interpolate_value(const array_view<T>& vals,
                                   const array_view<vec<int, N>>& elems,
                                   int eid, const vec3f& euv) {
    auto ret = T();
    if (vals.empty()) return ret;
    auto& elem = elems[eid];
    for (auto i = 0; i < N; i++) ret += vals[elem[i]] * euv[i];
    return ret;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static inline point _eval_shapepoint(const scene& scn, int shape_id, int eid,
                                     const vec3f& euv, const vec3f& wo) {
    // set shape data
    auto pt = point();

    // check if null point
    if (shape_id < 0) return pt;
    auto& shp = scn.shapes[shape_id];

    // direction
    pt.wo = wo;

    // compute points and weights
    auto pos = zero3f, norm = zero3f, color = zero3f;
    auto texcoord = zero2f;
    if (!shp.points.empty()) {
        pt.ptype = point::type::point;
        pos = _interpolate_value(shp.pos, shp.points, eid, euv);
        norm = normalize(_interpolate_value(shp.norm, shp.points, eid, euv));
        texcoord = _interpolate_value(shp.texcoord, shp.points, eid, euv);
        color = _interpolate_value(shp.color, shp.points, eid, euv);
    } else if (!shp.lines.empty()) {
        pt.ptype = point::type::line;
        pos = _interpolate_value(shp.pos, shp.lines, eid, euv);
        norm = normalize(_interpolate_value(shp.norm, shp.lines, eid, euv));
        texcoord = _interpolate_value(shp.texcoord, shp.lines, eid, euv);
        color = _interpolate_value(shp.color, shp.lines, eid, euv);
    } else if (!shp.triangles.empty()) {
        pt.ptype = point::type::triangle;
        pos = _interpolate_value(shp.pos, shp.triangles, eid, euv);
        norm = _interpolate_value(shp.norm, shp.triangles, eid, euv);
        texcoord = _interpolate_value(shp.texcoord, shp.triangles, eid, euv);
        color = _interpolate_value(shp.color, shp.triangles, eid, euv);
    }

    // creating frame
    pt.frame = make_frame3(pos, norm);

    // transform to world space
    pt.frame.o() = transform_point(shp.xform, pt.frame.o());
    pt.frame[2] = transform_direction(shp.xform, pt.frame[2]);
    pt.frame[0] = transform_direction(shp.xform, pt.frame[0]);
    pt.frame[1] = transform_direction(shp.xform, pt.frame[1]);

    // sample material data
    auto& mat = scn.materials[shp.matid];
    pt.ke = mat.ke;
    pt.kd = mat.kd;
    pt.ks = mat.ks;
    pt.rs = mat.rs;
    pt.use_phong = mat.use_phong;

    // handle surface color
    if (!shp.color.empty()) {
        pt.ke *= color;
        pt.kd *= color;
        pt.ks *= color;
    }

    // handle textures
    if (!shp.texcoord.empty()) {
        if (mat.ke_txt >= 0) {
            auto txt = _eval_texture(scn.textures[mat.ke_txt], texcoord);
            pt.ke = lerp(pt.ke, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.kd_txt >= 0) {
            auto txt = _eval_texture(scn.textures[mat.kd_txt], texcoord);
            pt.kd = lerp(pt.kd, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.ks_txt >= 0) {
            auto txt = _eval_texture(scn.textures[mat.ks_txt], texcoord);
            pt.ks = lerp(pt.ks, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.rs_txt >= 0) {
            auto txt = _eval_texture(scn.textures[mat.rs_txt], texcoord);
            pt.rs = lerp(pt.rs, txt[0], txt[3]);
        }
    }

    return pt;
}

//
// Sample weight for a light point.
//
static inline float _weight_light(const scene& scn, int light_id,
                                  const point& lpt, const point& pt) {
    switch (lpt.ptype) {
        case point::type::env: {
            return 4 * pif;
        } break;
        case point::type::point: {
            auto& light = scn._lights[light_id];
            auto d = dist(lpt.frame.o(), pt.frame.o());
            return light.area / (d * d);
        } break;
        case point::type::line: {
            assert(false);
            return 0;
        } break;
        case point::type::triangle: {
            auto& light = scn._lights[light_id];
            auto d = dist(lpt.frame.o(), pt.frame.o());
            return light.area * fabsf(dot(lpt.frame[2], lpt.wo)) / (d * d);
        } break;
        default: {
            assert(false);
            return 0;
        } break;
    }
}

//
// Picks a point on a light.
//
static inline point _sample_light(const scene& scn, int lid, const point& pt,
                                  float rne, const vec2f& rn) {
    auto& light = scn._lights[lid];
    if (light.shape_id >= 0) {
        auto& shp = scn.shapes[light.shape_id];
        auto eid = (int)(lower_bound(light.cdf.begin(), light.cdf.end(), rne) -
                         light.cdf.begin());
        if (eid > light.cdf.size() - 1) eid = (int)light.cdf.size() - 1;

        auto euv = zero3f;
        if (!shp.triangles.empty()) {
            euv = {sqrt(rn[0]) * (1 - rn[1]), 1 - sqrt(rn[0]),
                   rn[1] * sqrt(rn[0])};
        } else if (!shp.lines.empty()) {
            euv = {1 - rn[0], rn[0], 0};
        } else if (!shp.points.empty()) {
            euv = {1, 0, 0};
        } else
            assert(false);

        auto lpt = _eval_shapepoint(scn, light.shape_id, eid, euv, zero3f);
        lpt.wo = normalize(pt.frame.o() - lpt.frame.o());
        lpt.light_id = lid;
        return lpt;
    } else if (light.env_id >= 0) {
        auto z = -1 + 2 * rn[1];
        auto rr = sqrt(clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * pif * rn[0];
        auto wo = vec3f{cosf(phi) * rr, z, sinf(phi) * rr};
        auto lpt = _eval_envpoint(scn, light.env_id, wo);
        lpt.light_id = lid;
        return lpt;
    } else {
        assert(false);
    }
    return point();
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ray3f _offset_ray(const scene& scn, const point& pt,
                                const vec3f& w, const render_params& params) {
    return ray3f(pt.frame.o() + pt.frame[2] * params.ray_eps, w,
                 params.ray_eps);
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ray3f _offset_ray(const scene& scn, const point& pt,
                                const point& pt2, const render_params& params) {
    auto ray_dist = (pt2.ptype != point::type::env)
                        ? dist(pt.frame.o(), pt2.frame.o())
                        : FLT_MAX;
    return ray3f(pt.frame.o() + pt.frame[2] * params.ray_eps, -pt2.wo,
                 params.ray_eps, ray_dist - 2 * params.ray_eps);
}

//
// Intersects a ray with the scn and return the point (or env point).
//
static inline point _intersect_scene(const scene& scn, const ray3f& ray) {
    auto isec = scn.intersect_first(ray);
    if (isec) {
        return _eval_shapepoint(scn, isec.sid, isec.eid, isec.euv, -ray.d);
    } else if (!scn.environments.empty()) {
        return _eval_envpoint(scn, 0, -ray.d);
    } else {
        return {};
    }
}

//
// Evalutes direct illumination using MIS.
//
static inline vec3f _eval_direct(const scene& scn, int lid, const point& pt,
                                 _sampler& smp, const render_params& params) {
    // select whether it goes in all light mode
    auto all_lights = (lid < 0);

    // pick a light if not there
    auto nlweight = 0.0f;
    if (all_lights) {
        lid = _sample_next1f(smp) * scn._lights.size();
        if (lid > scn._lights.size() - 1) lid = (int)scn._lights.size() - 1;
        nlweight = scn._lights.size();
    } else {
        nlweight = 1;
    }

    // sample light according to area
    auto lpt =
        _sample_light(scn, lid, pt, _sample_next1f(smp), _sample_next2f(smp));
    auto lld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo);
    auto lweight = _weight_light(scn, lid, lpt, pt) * nlweight;
    lld *= lweight;
    if (lld != zero3f) {
        auto shadow_ray = _offset_ray(scn, pt, lpt, params);
        if (scn.intersect_any(shadow_ray)) lld = zero3f;
    }

    // check if mis is necessary
    if (pt.ptype == point::type::point || pt.ptype == point::type::line)
        return lld;

    // check if mis is necessary
    auto& light = scn._lights[lid];
    if (light.shape_id < 0) return lld;
    if (lpt.ptype == point::type::point || lpt.ptype == point::type::line) {
        return lld;
    }

    // sample the brdf
    auto bwi = _sample_brdfcos(pt, _sample_next1f(smp), _sample_next2f(smp));
    auto bweight = 0.0f;
    auto bld = zero3f;
    auto bpt = _intersect_scene(scn, _offset_ray(scn, pt, bwi, params));
    if (lid == bpt.light_id || all_lights) {
        bweight = _weight_brdfcos(pt, bwi);
        bld = _eval_emission(bpt) * _eval_brdfcos(pt, bwi) * bweight;
    }

    // accumulate the value with mis
    if (lld != zero3f) {
        auto bweight = _weight_brdfcos(pt, -lpt.wo);
        // float weight =
        //     (1 / lweight) * (1 / lweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        auto weight = (1 / lweight) / ((1 / lweight) + (1 / bweight));
        lld *= weight;
    }
    if (bld != zero3f) {
        auto lweight = _weight_light(scn, lid, bpt, pt) * nlweight;
        // float weight =
        //     (1 / bweight) * (1 / bweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        auto weight = (1 / bweight) / ((1 / lweight) + (1 / bweight));
        bld *= weight;
    }

    // return weighted sum
    return lld + bld;
}

//
// Recursive path tracing.
//
static inline vec4f _shade_pathtrace_recd(const scene& scn, const ray3f& ray,
                                          _sampler& smp, int ray_depth,
                                          const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none) return zero4f;

    // init
    auto la = vec4f{0, 0, 0, 1};
    auto& l = *(vec3f*)la.data();

    // emission
    if (ray_depth == 0) l += _eval_emission(pt);
    if (pt.ptype == point::type::env) return la;

    // check early exit
    if (pt.kd == zero3f && pt.ks == zero3f) return la;

    // direct
    l += _eval_direct(scn, 0, pt, smp, params);

    // roussian roulette
    if (ray_depth >= params.max_depth) return la;
    auto rrweight = 1.0f;
    if (ray_depth >= params.min_depth) {
        auto rrrn = _sample_next1f(smp);
        auto wrr = min(mean(pt.kd) + mean(pt.ks), 0.95f);
        if (rrrn >= wrr) return {l[0], l[1], l[2], 1};
        rrweight /= wrr;
    }

    // continue path
    auto bwi = _sample_brdfcos(pt, _sample_next1f(smp), _sample_next2f(smp));
    if (bwi == zero3f) return la;
    auto bweight = _weight_brdfcos(pt, bwi);
    if (!bweight) return la;
    auto bbrdfcos = _eval_brdfcos(pt, bwi);
    if (bbrdfcos == zero3f) return la;
    auto ble = _shade_pathtrace_recd(scn, _offset_ray(scn, pt, bwi, params),
                                     smp, ray_depth + 1, params);
    l += vec3f{ble[0], ble[1], ble[2]} * bbrdfcos * rrweight;

    return la;
}

//
// Shader interface for the above function.
//
static inline vec4f _shade_pathtrace(const scene& scn, const ray3f& ray,
                                     _sampler& smp,
                                     const render_params& params) {
    return _shade_pathtrace_recd(scn, ray, smp, 0, params);
}

//
// Direct illuination.
//
static inline vec4f _shade_direct(const scene& scn, const ray3f& ray,
                                  _sampler& smp, const render_params& params) {
    // scn intersection
    auto pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none) return zero4f;

    // init
    auto la = vec4f{0, 0, 0, 1};
    auto& l = *(vec3f*)la.data();

    // emission
    l += _eval_emission(pt);
    if (pt.ptype == point::type::env) return la;

    // early exit
    if (pt.kd == zero3f && pt.ks == zero3f) return la;

    // ambient
    l += params.amb * pt.kd;

    // direct
    for (int lid = 0; lid < scn._lights.size(); lid++) {
        l += _eval_direct(scn, lid, pt, smp, params);
    }

    // done
    return la;
}

//
// Eyelight for quick previewing.
//
static inline vec4f _shade_eyelight(const scene& scn, const ray3f& ray,
                                    _sampler& smp,
                                    const render_params& params) {
    // intersection
    point pt = _intersect_scene(scn, ray);
    if (pt.ptype == point::type::none) return zero4f;

    // init
    auto la = vec4f{0, 0, 0, 1};
    auto& l = *(vec3f*)la.data();

    // emission
    l += _eval_emission(pt);
    if (pt.ptype == point::type::env) return la;

    // brdf*light
    l += _eval_brdfcos(pt, pt.wo) * pif;

    return la;
}

//
// Shader function callback.
//
using shade_fn = function<vec4f(const scene& scn, const ray3f& ray,
                                _sampler& smp, const render_params& params)>;

//
// Renders a block of pixels. Public API, see above.
//
YGL_API void trace_block(const scene& scn, int cid, image_view<vec4f> img,
                         int ns, const vec2i& xy, const vec2i& wh,
                         const vec2i& samples, const render_params& params,
                         bool accumulate) {
    auto& cam = scn.cameras[cid];
    shade_fn shade;
    switch (params.stype) {
        case shader_type::eyelight: shade = _shade_eyelight; break;
        case shader_type::def:
        case shader_type::direct: shade = _shade_direct; break;
        case shader_type::pathtrace: shade = _shade_pathtrace; break;
        default: assert(false); return;
    }
    for (auto j = xy[1]; j < xy[1] + wh[1]; j++) {
        for (auto i = xy[0]; i < xy[0] + wh[0]; i++) {
            auto ij = vec2i(i, j);
            auto saved = img[ij];
            img[ij] = zero4f;
            for (auto s = samples[0]; s < samples[1]; s++) {
                auto smp = _make_sampler(i, j, s, ns, params.rtype);
                auto rn = _sample_next2f(smp);
                auto uv = vec2f{(i + rn[0]) / img.size()[0],
                                1 - (j + rn[1]) / img.size()[1]};
                auto ray = _eval_camera(cam, uv, _sample_next2f(smp));
                auto l = shade(scn, ray, smp, params);
                if (!isfinite(l[0]) || !isfinite(l[1]) || !isfinite(l[2]))
                    continue;
                if (params.pixel_clamp > 0)
                    *(vec3f*)&l = clamplen(*(vec3f*)&l, params.pixel_clamp);
                img[ij] += l;
            }
            if (accumulate && samples[0]) {
                img[ij] += saved * samples[0];
                img[ij] /= samples[1];
            } else {
                img[ij] /= samples[1] - samples[0];
            }
        }
    }
}

//
// Renders the whole image. Public API, see above.
//
YGL_API void trace_image(const scene& scn, int cid, image_view<vec4f> img,
                         int ns, const render_params& params) {
    trace_block(scn, cid, img, ns, {0, 0}, img.size(), {0, ns}, params);
}

}  // namespace

#endif

#endif
