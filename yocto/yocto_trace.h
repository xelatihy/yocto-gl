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
//     scene.intersect_first = <callback>
//     scene.intersect_any = <callback>
//     - can use yocto_bvh
// 2. prepare for rendering
//    init_lights(scene)
// 3. define rendering params
// 4. render blocks of samples
//    render_block(scene, pixels, image size, block)
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

#ifndef _YTRACE_H_
#define _YTRACE_H_

#ifndef YGL_DECLARATION
#define YGL_API inline
#else
#define YGL_API
#endif

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// C++ INTERFACE
// -----------------------------------------------------------------------------

namespace ytrace {

//
// Type of rendering algorithm (shader)
//
enum struct stype {
    def = 0,    // default renderer
    eyelight,   // eye hight for quick previews
    direct,     // direct illumination
    pathtrace,  // path tracing
};

//
// Random number generator type
//
enum struct rtype {
    def = 0,     // default generator
    uniform,     // uniform random numbers
    stratified,  // stratified random numbers
    cmjs,        // correlated multi-jittered sampling
};

//
// Camera
//
struct camera {
    ym::frame3f xform = ym::identity_frame3f;  // local-to-world transform
    float yfov = ym::pif / 3;                  // field of view
    float aspect = 1;                          // aspect ratio
    float aperture = 0;                        // lens aperture
    float focus = 1;  // focus plane distance (cannot be zero)
};

//
// Texture
//
struct texture {
    ym::image_view<ym::vec4f> hdr;
    ym::image_view<ym::vec4b> ldr;
};

//
// Material
//
struct material {
    // material values
    ym::vec3f ke = ym::zero3f;  // emission, term
    ym::vec3f kd = ym::zero3f;  // diffuse term
    ym::vec3f ks = ym::zero3f;  // specular term
    float rs = 0.1;             // specular roughness API

    // fresnel
    ym::vec3f es = ym::zero3f;   // eta
    ym::vec3f eks = ym::zero3f;  // etak (metals only)

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
    ym::frame3f xform = ym::identity_frame3f;  // local-to-world rigid transform
    int matid = -1;                            // material id

    // element data [only one enabled at any given time]
    ym::array_view<ym::vec1i> points;     // elem data
    ym::array_view<ym::vec2i> lines;      // elem data
    ym::array_view<ym::vec3i> triangles;  // elem data

    // vertex data
    ym::array_view<ym::vec3f> pos;       // vertex data
    ym::array_view<ym::vec3f> norm;      // vertex data
    ym::array_view<ym::vec2f> texcoord;  // vertex data
    ym::array_view<ym::vec3f> color;     // vertex data
    ym::array_view<float> radius;        // vertex data
};

//
// Environment
//
struct environment {
    ym::frame3f xform = ym::identity_frame3f;  // local-to-world rigid transform
    ym::vec3f ke = ym::zero3f;                 // emission
    int ke_txt = -1;                           // emission texture
};

//
// Light (either shape or environment).
// This is only used internally and should not be created.
//
struct _light {
    int shape_id = -1;       // shape
    int env_id = -1;         // environment
    std::vector<float> cdf;  // for shape, cdf of shape elements for sampling
    float area = 0;          // for shape, shape area
};

//
// Ray-scene Intersection.
//
struct intersect_point {
    float dist = 0;              // ray distance
    int sid = -1;                // shape index
    int eid = -1;                // element index
    ym::vec3f euv = ym::zero3f;  // element baricentric coordinates
    bool hit = false;            // whether we hit

    // check whether it was a hit
    operator bool() const { return hit; }
};

//
// Ray-scene closest intersection callback
//
// Parameters:
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Out Parameters:
// - ray_t: hit distance
// - sid: hit shape index
// - eid: hit element index
// - euv: hit element parameters
//
// Return:
// - whether we intersect or not
//
using intersect_first_cb = std::function<intersect_point(const ym::ray3f& ray)>;

//
// Ray-scene intersection callback
//
// Parameters:
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Return:
// - whether we intersect or not
//
using intersect_any_cb = std::function<bool(const ym::ray3f& ray)>;

//
// Scene
//
struct scene {
    // intersection callbaks
    intersect_first_cb intersect_first;  // ray intersection callback
    intersect_any_cb intersect_any;      // ray hit callback

    // scene data
    std::vector<camera> cameras;            // camera
    std::vector<environment> environments;  // env
    std::vector<shape> shapes;              // shapes
    std::vector<material> materials;        // materials
    std::vector<texture> textures;          // textures

    // [private] light sources
    std::vector<_light> _lights;  // lights [private]
};

//
// Rendering params
//
struct render_params {
    stype stype = stype::def;    // sampler type
    rtype rtype = rtype::def;    // random type
    ym::vec3f amb = ym::zero3f;  // ambient lighting
    int min_depth = 3;           // min ray depth
    int max_depth = 8;           // mas ray depth
    float pixel_clamp = 100;     // final pixel clamping
    float ray_eps = 1e-2f;       // ray intersection epsilon
};

//
// Convert a Phong exponent to GGX/Phong roughness
//
YGL_API float specular_exponent_to_roughness(float n);

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YGL_API void specular_fresnel_from_ks(const ym::vec3f& ks, ym::vec3f& es,
                                      ym::vec3f& esk);

//
// Initialize rendering.
//
// Parameters:
// - scene: trace scene
//
YGL_API void init_lights(scene& scene);

//
// Renders a block of sample
//
// Parameters:
// - scene: trace scene
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
YGL_API void trace_block(const scene& scene, int cid,
                         ym::image_view<ym::vec4f> img, int ns,
                         const ym::range2i& block, const ym::range1i& samples,
                         const render_params& params, bool accumulate = false);

//
// Convenience function to call trace_block with all sample at once.
//
YGL_API void trace_image(const scene& scene, int cid,
                         ym::image_view<ym::vec4f> img, int ns,
                         const render_params& params);

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace ytrace {

//
// Compute shape element cdf for shape sampling.
//
template <typename T, typename Weight_callback>
static inline std::vector<float> _compute_weight_cdf(
    const ym::array_view<T>& elem, float& total_weight,
    const Weight_callback& weight_cb) {
    // prepare return
    auto cdf = std::vector<float>();
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
YGL_API void init_lights(scene& scene) {
    // clear old lights
    scene._lights.resize(0);

    for (int sid = 0; sid < scene.shapes.size(); sid++) {
        auto& shape = scene.shapes[sid];
        auto& mat = scene.materials[shape.matid];
        if (mat.ke == ym::zero3f) continue;
        scene._lights.push_back(_light());
        auto& light = scene._lights.back();
        light.shape_id = sid;
        if (!shape.points.empty()) {
            light.cdf = _compute_weight_cdf(shape.points, light.area,
                                            [&shape](auto e) { return 1; });
        } else if (!shape.lines.empty()) {
            light.cdf =
                _compute_weight_cdf(shape.lines, light.area, [&shape](auto e) {
                    return ym::length(shape.pos[e[1]] - shape.pos[e[0]]);
                });
        } else if (!shape.triangles.empty()) {
            light.cdf = _compute_weight_cdf(
                shape.triangles, light.area, [&shape](auto e) {
                    return ym::triangle_area(shape.pos[e[0]], shape.pos[e[1]],
                                             shape.pos[e[2]]);
                });
        }
    }

    for (int envid = 0; envid < scene.environments.size(); envid++) {
        auto& env = scene.environments[envid];
        if (env.ke == ym::zero3f) continue;
        scene._lights.push_back(_light());
        auto& light = scene._lights.back();
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
YGL_API void specular_fresnel_from_ks(const ym::vec3f& ks, ym::vec3f& es,
                                      ym::vec3f& esk) {
    es = (ym::one3f + ym::sqrt(ks)) / (ym::one3f - ym::sqrt(ks));
    esk = ym::zero3f;
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// Random number sampler. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
//
struct _sampler {
    ym::rng_pcg32 rng;  // rnumber number state
    int i, j;           // pixel coordinates
    int s, d;           // sample and dimension indices
    int ns;             // number of samples
    rtype rtype;        // random number type
};

//
// Initialize a sampler ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
//
static inline _sampler _make_sampler(int i, int j, int s, int ns, rtype rtype) {
    // we use various hashes to scramble the pixel values
    _sampler sampler = {{0, 0}, i, j, s, 0, ns, rtype};
    uint64_t sample_id = ((uint64_t)(i + 1)) << 0 | ((uint64_t)(j + 1)) << 15 |
                         ((uint64_t)(s + 1)) << 30;
    uint64_t initseq = ym::hash_uint64(sample_id);
    uint64_t initstate =
        ym::hash_uint64(sample_id * 3202034522624059733ull + 1ull);
    ym::rng_init(sampler.rng, initstate, initseq);
    return sampler;
}

//
// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
//
static inline float _sample_next1f(_sampler& sampler) {
    float rn = 0;
    switch (sampler.rtype) {
        case rtype::def:
        case rtype::uniform: {
            rn = ym::rng_nextf(sampler.rng);
        } break;
        case rtype::stratified: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(sampler.i + 1)) << 0 |
                                            ((uint64_t)(sampler.j + 1)) << 15 |
                                            ((uint64_t)(sampler.d + 1)) << 30);
            int s = ym::hash_permute(sampler.s, sampler.ns, p);
            rn = (s + ym::rng_nextf(sampler.rng)) / sampler.ns;
        } break;
        case rtype::cmjs: {
            uint32_t p = ym::hash_uint64_32(((uint64_t)(sampler.i + 1)) << 0 |
                                            ((uint64_t)(sampler.j + 1)) << 15 |
                                            ((uint64_t)(sampler.d + 1)) << 30);
            int s = ym::hash_permute(sampler.s, sampler.ns, p);
            rn = (s + ym::hash_randfloat(s, p * 0xa399d265)) / sampler.ns;
        } break;
        default: assert(false);
    }

    sampler.d += 1;

    // make sure all sampled numbers are below 1
    // TODO: use std::numeric_limits
    if (rn >= 1) rn = 1 - FLT_EPSILON;

    return rn;
}

//
// Generates a 1-dimensional sample.
//
// Implementation notes: see above. Note that using deterministic keyed
// permutaton we can use stratified sampling without preallocating samples.
//
static inline ym::vec2f _sample_next2f(_sampler& sampler) {
    ym::vec2f rn = {0, 0};
    switch (sampler.rtype) {
        case rtype::def:
        case rtype::uniform: {
            rn[0] = ym::rng_nextf(sampler.rng);
            rn[1] = ym::rng_nextf(sampler.rng);
        } break;
        case rtype::stratified: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler.ns));
            uint32_t p = ym::hash_uint64_32(((uint64_t)(sampler.i + 1)) << 0 |
                                            ((uint64_t)(sampler.j + 1)) << 15 |
                                            ((uint64_t)(sampler.d + 1)) << 30);
            int s = ym::hash_permute(sampler.s, sampler.ns, p);
            rn[0] = (s % ns2 + ym::rng_nextf(sampler.rng)) / ns2;
            rn[1] = (s / ns2 + ym::rng_nextf(sampler.rng)) / ns2;
        } break;
        case rtype::cmjs: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler.ns));
            uint32_t p = ym::hash_uint64_32(((uint64_t)(sampler.i + 1)) << 0 |
                                            ((uint64_t)(sampler.j + 1)) << 15 |
                                            ((uint64_t)(sampler.d + 1)) << 30);
            int s = ym::hash_permute(sampler.s, sampler.ns, p);
            int sx = ym::hash_permute(s % ns2, ns2, p * 0xa511e9b3);
            int sy = ym::hash_permute(s / ns2, ns2, p * 0x63d83595);
            float jx = ym::hash_randfloat(s, p * 0xa399d265);
            float jy = ym::hash_randfloat(s, p * 0x711ad6a5);
            rn[0] = (s % ns2 + (sy + jx) / ns2) / ns2;
            rn[1] = (s / ns2 + (sx + jy) / ns2) / ns2;
        } break;
        default: assert(false);
    }

    sampler.d += 2;

    // make sure all sampled numbers are below 1
    if (rn[0] >= 1) rn[0] = 1 - FLT_EPSILON;
    if (rn[1] >= 1) rn[1] = 1 - FLT_EPSILON;

    return rn;
}

//
// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
//
struct _point {
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
    ym::vec3f wo = ym::zero3f;  // outgoing direction

    // resolved geometry (shape) ------------
    ym::frame3f frame = ym::identity_frame3f;  // local frame

    // shading ------------------------------
    ym::vec3f ke = ym::zero3f;   // material values
    ym::vec3f kd = ym::zero3f;   // material values
    ym::vec3f ks = ym::zero3f;   // material values
    float rs = 0;                // material values
    ym::vec3f es = ym::zero3f;   // material values
    ym::vec3f eks = ym::zero3f;  // material values
    bool use_phong = false;      // material values
};

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static inline ym::ray3f _eval_camera(const camera& cam, const ym::vec2f& uv,
                                     const ym::vec2f& luv) {
    auto h = 2 * std::tan(cam.yfov / 2);
    auto w = h * cam.aspect;
    ym::vec3f o = ym::vec3f{luv[0] * cam.aperture, luv[1] * cam.aperture, 0};
    ym::vec3f q = {w * cam.focus * (uv[0] - 0.5f),
                   h * cam.focus * (uv[1] - 0.5f), -cam.focus};
    return ym::ray3f(ym::transform_point(cam.xform, o),
                     ym::transform_direction(cam.xform, ym::normalize(q - o)));
}

//
// Wrapper for above function
//
static inline ym::vec4f _eval_texture(const texture& txt,
                                      const ym::vec2f& texcoord) {
    assert(!txt.hdr.empty() || !txt.ldr.empty());

    // get image width/height
    auto wh = (!txt.ldr.empty()) ? txt.ldr.size() : txt.hdr.size();

    // get coordinates normalized for tiling
    auto st =
        ym::vec2f{std::fmod(texcoord[0], 1.0f), std::fmod(texcoord[1], 1.0f)} *
        ym::vec2f(wh);
    if (st[0] < 0) st[0] += wh[0];
    if (st[1] < 0) st[1] += wh[1];

    // get image coordinates and residuals
    auto ij = ym::clamp(ym::vec2i(st), {0, 0}, wh);
    auto uv = st - ym::vec2f(ij);

    // get interpolation weights and indices
    ym::vec2i idx[4] = {ij,
                        {ij[0], (ij[1] + 1) % wh[1]},
                        {(ij[0] + 1) % wh[0], ij[1]},
                        {(ij[0] + 1) % wh[0], (ij[1] + 1) % wh[1]}};
    float w[4] = {(1 - uv[0]) * (1 - uv[1]), (1 - uv[0]) * uv[1],
                  uv[0] * (1 - uv[1]), uv[0] * uv[1]};

    // handle interpolation
    if (!txt.ldr.empty()) {
        return (ym::srgb_to_linear(txt.ldr[idx[0]]) * w[0] +
                ym::srgb_to_linear(txt.ldr[idx[1]]) * w[1] +
                ym::srgb_to_linear(txt.ldr[idx[2]]) * w[2] +
                ym::srgb_to_linear(txt.ldr[idx[3]]) * w[3]);
    } else if (!txt.hdr.empty()) {
        return (txt.hdr[idx[0]] * w[0] + txt.hdr[idx[1]] * w[1] +
                txt.hdr[idx[2]] * w[2] + txt.hdr[idx[3]] * w[3]);
    } else {
        assert(false);
    }

    // should not have gotten here
    return ym::zero4f;
}

//
// Evaluates emission.
//
static inline ym::vec3f _eval_emission(const _point& pt) {
    if (pt.ke == ym::zero3f) return ym::zero3f;
    switch (pt.ptype) {
        case _point::type::env: return pt.ke;
        case _point::type::point: return pt.ke;
        case _point::type::line: return pt.ke;
        case _point::type::triangle:
            return (ym::dot(pt.frame.z(), pt.wo) > 0) ? pt.ke : ym::zero3f;
        default: {
            assert(false);
            return ym::zero3f;
        }
    }
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline ym::vec3f _eval_fresnel_dielectric(float cosw,
                                                 const ym::vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1 / eta;
        cosw = -cosw;
    }

    auto sin2 = ym::vec3f(1 - cosw * cosw);
    auto eta2 = eta * eta;

    auto cos2t = ym::one3f - sin2 / eta2;
    if (cos2t[0] < 0 || cos2t[1] < 0 || cos2t[2] < 0) return ym::one3f;  // tir

    auto t0 = sqrt(cos2t);
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (ym::vec3f(cosw) - t1) / (ym::vec3f(cosw) + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2;
}

//
// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline ym::vec3f _eval_fresnel_metal(float cosw, const ym::vec3f& eta,
                                            const ym::vec3f& etak) {
    if (etak == ym::zero3f) return _eval_fresnel_dielectric(cosw, eta);

    cosw = ym::clamp(cosw, -1, 1);
    auto cos2 = ym::vec3f(cosw * cosw);
    auto sin2 = ym::clamp(ym::one3f - cos2, 0, 1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - sin2;
    auto a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
    auto t1 = a2plusb2 + cos2;
    auto a = sqrt(0.5f * (a2plusb2 + t0));
    auto t2 = 2 * a * cosw;
    auto rs = (t1 - t2) / (t1 + t2);

    auto t3 = ym::vec3f(cos2) * a2plusb2 + sin2 * sin2;
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
static inline ym::vec3f _eval_brdfcos(const _point& pt, const ym::vec3f& wi) {
    // summation over multiple terms
    auto brdfcos = ym::zero3f;

    // exit if not needed
    if (pt.kd == ym::zero3f && pt.ks == ym::zero3f) return brdfcos;

    // save wo
    auto wo = pt.wo;

    // compute wh
    auto wh = ym::normalize(wo + wi);

    // compute dot products
    auto ndo = ym::dot(pt.frame.z(), wo), ndi = ym::dot(pt.frame.z(), wi),
         ndh = ym::clamp(ym::dot(wh, pt.frame.z()), 0, 1);

    switch (pt.ptype) {
        case _point::type::point: {
            // diffuse term (hack for now)
            if (pt.kd != ym::zero3f) {
                auto ido = ym::dot(wo, wi);
                auto diff = pt.kd * (2 * ido + 1) / (2 * ym::pif);
                brdfcos += diff;
            }
        } break;
        case _point::type::line: {
            // take sines
            auto so = sqrtf(ym::clamp(1 - ndo * ndo, 0, 1)),
                 si = sqrtf(ym::clamp(1 - ndi * ndi, 0, 1)),
                 sh = sqrtf(ym::clamp(1 - ndh * ndh, 0, 1));

            // diffuse term (Kajiya-Kay)
            if (si > 0 && so > 0 && pt.kd != ym::zero3f) {
                auto diff = pt.kd * si / ym::pif;
                brdfcos += diff;
            }

            // specular term (Kajiya-Kay)
            if (si > 0 && so > 0 && sh > 0 && pt.ks != ym::zero3f) {
                auto ns = 2 / (pt.rs * pt.rs) - 2;
                auto d = (ns + 2) * powf(sh, ns) / (2 + ym::pif);
                auto spec = pt.ks * si * d / (4 * si * so);
                brdfcos += spec;
            }
        } break;
        case _point::type::triangle: {
            // diffuse term
            if (ndi > 0 && ndo && pt.kd != ym::zero3f) {
                auto diff = pt.kd * ndi / ym::pif;
                brdfcos += diff;
            }

            // specular term (GGX)
            if (ndi > 0 && ndo > 0 && ndh > 0 && pt.ks != ym::zero3f) {
                if (!pt.use_phong) {
                    // evaluate GGX
                    auto cos2 = ndh * ndh;
                    auto tan2 = (1 - cos2) / cos2;
                    auto alpha2 = pt.rs * pt.rs;
                    auto d = alpha2 / (ym::pif * cos2 * cos2 * (alpha2 + tan2) *
                                       (alpha2 + tan2));
                    auto lambda_o =
                        (-1 + sqrtf(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    auto lambda_i =
                        (-1 + sqrtf(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    auto g = 1 / (1 + lambda_o + lambda_i);
                    auto spec = pt.ks * ndi * d * g / (4 * ndi * ndo);
                    if (pt.es != ym::zero3f) {
                        spec *= _eval_fresnel_metal(ndh, pt.es, pt.eks);
                    }
                    brdfcos += spec;
                } else {
                    // evaluate Blinn-Phong
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    auto d = (ns + 2) * powf(ndh, ns) / (2 + ym::pif);
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
static inline float _weight_brdfcos(const _point& pt, const ym::vec3f& wi) {
    // skip if no component
    if (pt.kd == ym::zero3f && pt.ks == ym::zero3f) return 0;

    // save wo
    auto wo = pt.wo;

    // compute wh
    auto wh = ym::normalize(wi + wo);

    // compute dot products
    auto ndo = ym::dot(pt.frame.z(), wo), ndi = ym::dot(pt.frame.z(), wi),
         ndh = ym::dot(pt.frame.z(), wh);

    // check to make sure we are above the surface
    // updated this for refraction
    if (ndo <= 0 || ndi <= 0) return 0;

    // pick from a sum
    auto wall = ym::mean(pt.kd) + ym::mean(pt.ks);
    auto wd = ym::mean(pt.kd) / wall;
    auto ws = ym::mean(pt.ks) / wall;

    // accumulate probability
    auto pdf = 0.0f;

    switch (pt.ptype) {
        case _point::type::point: {
        } break;
        case _point::type::line: {
            // diffuse term
            if (wall) {
                // homepherical cosine probability
                pdf += 1 / (4 * ym::pif);
            }
        } break;
        case _point::type::triangle: {
            // diffuse term
            if (wd && ndi > 0) {
                // homepherical cosine probability
                pdf += wd * ndi / ym::pif;
            }

            // specular term (GGX or Phong)
            if (ws && ndi > 0 && ndo > 0 && ndh > 0) {
                if (!pt.use_phong) {
                    // probability proportional to d * ndh
                    auto cos2 = ndh * ndh;
                    auto tan2 = (1 - cos2) / cos2;
                    auto alpha2 = pt.rs * pt.rs;
                    auto d = alpha2 / (ym::pif * cos2 * cos2 * (alpha2 + tan2) *
                                       (alpha2 + tan2));
                    auto hdo = ym::dot(wo, wh);
                    pdf += ws * d * ndh / (4 * hdo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // compute wh
                    auto wh = ym::normalize(wi + wo);
                    auto ndh = ym::dot(pt.frame.z(), wh);
                    // homerispherical cosine power probability
                    pdf += ws * powf(ndh, ns) * (ns + 1) / (2 * ym::pif);
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
static inline ym::vec3f _sample_brdfcos(const _point& pt, float rnl,
                                        const ym::vec2f& rn) {
    // skip if no component
    if (pt.kd == ym::zero3f && pt.ks == ym::zero3f) return ym::zero3f;

    // save wo
    auto wo = pt.wo;

    // compute cosine
    auto ndo = ym::dot(pt.frame.z(), wo);

    // check to make sure we are above the surface
    // update this for refraction
    if (ndo <= 0) return ym::zero3f;

    // pick from a sum
    auto wall = ym::mean(pt.kd) + ym::mean(pt.ks);
    auto wd = ym::mean(pt.kd) / wall;
    auto ws = ym::mean(pt.ks) / wall;

    switch (pt.ptype) {
        // TODO: point color
        case _point::type::point:
        case _point::type::line: {
            if (wall > 0) {
                // sample wi with uniform spherical distribution
                auto rz = rn[1], rr = sqrtf(1 - rz * rz),
                     rphi = 2 * ym::pif * rn[0];
                auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return ym::transform_direction(pt.frame, wi_local);
            }
        } break;
        case _point::type::triangle: {
            // sample according to diffuse
            if (rnl < wd) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn[1]), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * ym::pif * rn[0];
                // set to wi
                auto wi_local = ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return ym::transform_direction(pt.frame, wi_local);
            }

            // sample according to specular (GGX or Phong)
            if (rnl >= wd && rnl < wd + ws) {
                if (!pt.use_phong) {
                    // sample wh with ggx distribution
                    auto tan2 = pt.rs * pt.rs * rn[1] / (1 - rn[1]);
                    auto rz = sqrtf(1 / (tan2 + 1)), rr = sqrtf(1 - rz * rz),
                         rphi = 2 * ym::pif * rn[0];
                    // set to wh
                    auto wh_local =
                        ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                    auto wh = ym::transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym::normalize(wh * 2 * ym::dot(wo, wh) - wo);
                } else {
                    // get phong exponent
                    auto ns = 2 / (pt.rs * pt.rs) - 2;
                    // sample wh with hemispherical cosine power distribution
                    auto rz = powf(rn[1], 1 / (ns + 1)),
                         rr = sqrtf(1 - rz * rz), rphi = 2 * ym::pif * rn[0];
                    // set to wh
                    auto wh_local =
                        ym::vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                    auto wh = ym::transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym::normalize(wh * 2 * ym::dot(wo, wh) - wo);
                }
            }
        } break;
        default: { assert(false); }
    }

    // should not have gotten here
    assert(false);
    return ym::zero3f;
}

//
// Create a point for an environment map. Resolves material with textures.
//
static inline _point _eval_envpoint(const scene& scene, int env_id,
                                    const ym::vec3f& wo) {
    // set shape data
    auto pt = _point();

    // check if null point
    if (env_id < 0) return pt;
    auto& env = scene.environments[env_id];

    // env params
    pt.ptype = _point::type::env;

    // direction
    pt.wo = wo;

    // maerial
    pt.ke = env.ke;

    // textures
    if (env.ke_txt >= 0) {
        auto w = ym::transform_direction(ym::inverse(env.xform), -wo);
        auto theta = 1 - (acosf(ym::clamp(w[1], -1, 1)) / ym::pif);
        auto phi = atan2f(w.z(), w.x()) / (2 * ym::pif);
        auto texcoord = ym::vec2f{phi, theta};
        if (env.ke_txt >= 0) {
            auto txt = _eval_texture(scene.textures[env.ke_txt], texcoord);
            pt.ke = ym::lerp({txt[0], txt[1], txt[2]}, pt.ke, txt[3]);
        }
    }

    // done
    return pt;
}

//
// Interpolate a value over an element
//
template <typename T, int N>
static inline T _interpolate_value(const ym::array_view<T>& vals,
                                   const ym::array_view<ym::vec<int, N>>& elems,
                                   int eid, const ym::vec3f& euv) {
    auto ret = T();
    if (vals.empty()) return ret;
    auto& elem = elems[eid];
    for (auto i = 0; i < N; i++) ret += vals[elem[i]] * euv[i];
    return ret;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static inline _point _eval_shapepoint(const scene& scene, int shape_id, int eid,
                                      const ym::vec3f& euv,
                                      const ym::vec3f& wo) {
    // set shape data
    auto pt = _point();

    // check if null point
    if (shape_id < 0) return pt;
    auto& shape = scene.shapes[shape_id];

    // direction
    pt.wo = wo;

    // compute points and weights
    auto pos = ym::zero3f, norm = ym::zero3f, color = ym::zero3f;
    auto texcoord = ym::zero2f;
    if (!shape.points.empty()) {
        pt.ptype = _point::type::point;
        pos = _interpolate_value(shape.pos, shape.points, eid, euv);
        norm = ym::normalize(
            _interpolate_value(shape.norm, shape.points, eid, euv));
        texcoord = _interpolate_value(shape.texcoord, shape.points, eid, euv);
        color = _interpolate_value(shape.color, shape.points, eid, euv);
    } else if (!shape.lines.empty()) {
        pt.ptype = _point::type::line;
        pos = _interpolate_value(shape.pos, shape.lines, eid, euv);
        norm = ym::normalize(
            _interpolate_value(shape.norm, shape.lines, eid, euv));
        texcoord = _interpolate_value(shape.texcoord, shape.lines, eid, euv);
        color = _interpolate_value(shape.color, shape.lines, eid, euv);
    } else if (!shape.triangles.empty()) {
        pt.ptype = _point::type::triangle;
        pos = _interpolate_value(shape.pos, shape.triangles, eid, euv);
        norm = _interpolate_value(shape.norm, shape.triangles, eid, euv);
        texcoord =
            _interpolate_value(shape.texcoord, shape.triangles, eid, euv);
        color = _interpolate_value(shape.color, shape.triangles, eid, euv);
    }

    // creating frame
    pt.frame = ym::make_frame3(pos, norm);

    // transform to world space
    pt.frame.o() = ym::transform_point(shape.xform, pt.frame.o());
    pt.frame.z() = ym::transform_direction(shape.xform, pt.frame.z());
    pt.frame.x() = ym::transform_direction(shape.xform, pt.frame.x());
    pt.frame.y() = ym::transform_direction(shape.xform, pt.frame.y());

    // sample material data
    auto& mat = scene.materials[shape.matid];
    pt.ke = mat.ke;
    pt.kd = mat.kd;
    pt.ks = mat.ks;
    pt.rs = mat.rs;
    pt.use_phong = mat.use_phong;

    // handle surface color
    if (!shape.color.empty()) {
        pt.ke *= color;
        pt.kd *= color;
        pt.ks *= color;
    }

    // handle textures
    if (!shape.texcoord.empty()) {
        if (mat.ke_txt >= 0) {
            auto txt = _eval_texture(scene.textures[mat.ke_txt], texcoord);
            pt.ke = ym::lerp(pt.ke, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.kd_txt >= 0) {
            auto txt = _eval_texture(scene.textures[mat.kd_txt], texcoord);
            pt.kd = ym::lerp(pt.kd, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.ks_txt >= 0) {
            auto txt = _eval_texture(scene.textures[mat.ks_txt], texcoord);
            pt.ks = ym::lerp(pt.ks, {txt[0], txt[1], txt[2]}, txt[3]);
        }
        if (mat.rs_txt >= 0) {
            auto txt = _eval_texture(scene.textures[mat.rs_txt], texcoord);
            pt.rs = ym::lerp(pt.rs, txt[0], txt[3]);
        }
    }

    return pt;
}

//
// Sample weight for a light point.
//
static inline float _weight_light(const scene& scene, int light_id,
                                  const _point& lpt, const _point& pt) {
    switch (lpt.ptype) {
        case _point::type::env: {
            return 4 * ym::pif;
        } break;
        case _point::type::point: {
            auto& light = scene._lights[light_id];
            auto d = ym::dist(lpt.frame.o(), pt.frame.o());
            return light.area / (d * d);
        } break;
        case _point::type::line: {
            assert(false);
            return 0;
        } break;
        case _point::type::triangle: {
            auto& light = scene._lights[light_id];
            auto d = ym::dist(lpt.frame.o(), pt.frame.o());
            return light.area * fabsf(ym::dot(lpt.frame.z(), lpt.wo)) / (d * d);
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
static inline _point _sample_light(const scene& scene, int lid,
                                   const _point& pt, float rne,
                                   const ym::vec2f& rn) {
    auto& light = scene._lights[lid];
    if (light.shape_id >= 0) {
        auto& shape = scene.shapes[light.shape_id];
        auto eid =
            (int)(std::lower_bound(light.cdf.begin(), light.cdf.end(), rne) -
                  light.cdf.begin());
        if (eid > light.cdf.size() - 1) eid = (int)light.cdf.size() - 1;

        auto euv = ym::zero3f;
        if (!shape.triangles.empty()) {
            euv = {std::sqrt(rn[0]) * (1 - rn[1]), 1 - std::sqrt(rn[0]),
                   rn[1] * std::sqrt(rn[0])};
        } else if (!shape.lines.empty()) {
            euv = {1 - rn[0], rn[0], 0};
        } else if (!shape.points.empty()) {
            euv = {1, 0, 0};
        } else
            assert(false);

        auto lpt =
            _eval_shapepoint(scene, light.shape_id, eid, euv, ym::zero3f);
        lpt.wo = ym::normalize(pt.frame.o() - lpt.frame.o());
        lpt.light_id = lid;
        return lpt;
    } else if (light.env_id >= 0) {
        auto z = -1 + 2 * rn[1];
        auto rr = sqrtf(ym::clamp(1 - z * z, 0, 1));
        auto phi = 2 * ym::pif * rn[0];
        auto wo = ym::vec3f{cosf(phi) * rr, z, sinf(phi) * rr};
        auto lpt = _eval_envpoint(scene, light.env_id, wo);
        lpt.light_id = lid;
        return lpt;
    } else {
        assert(false);
    }
    return _point();
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f _offset_ray(const scene& scene, const _point& pt,
                                    const ym::vec3f& w,
                                    const render_params& params) {
    return ym::ray3f(pt.frame.o() + pt.frame.z() * params.ray_eps, w,
                     params.ray_eps);
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym::ray3f _offset_ray(const scene& scene, const _point& pt,
                                    const _point& pt2,
                                    const render_params& params) {
    auto ray_dist = (pt2.ptype != _point::type::env)
                        ? ym::dist(pt.frame.o(), pt2.frame.o())
                        : FLT_MAX;
    return ym::ray3f(pt.frame.o() + pt.frame.z() * params.ray_eps, -pt2.wo,
                     params.ray_eps, ray_dist - 2 * params.ray_eps);
}

//
// Intersects a ray with the scene and return the point (or env point).
//
static inline _point _intersect_scene(const scene& scene,
                                      const ym::ray3f& ray) {
    auto isec = scene.intersect_first(ray);
    if (isec) {
        return _eval_shapepoint(scene, isec.sid, isec.eid, isec.euv, -ray.d);
    } else if (!scene.environments.empty()) {
        return _eval_envpoint(scene, 0, -ray.d);
    } else {
        return {};
    }
}

//
// Evalutes direct illumination using MIS.
//
static inline ym::vec3f _eval_direct(const scene& scene, int lid,
                                     const _point& pt, _sampler& sampler,
                                     const render_params& params) {
    // select whether it goes in all light mode
    auto all_lights = (lid < 0);

    // pick a light if not there
    auto nlweight = 0.0f;
    if (all_lights) {
        lid = _sample_next1f(sampler) * scene._lights.size();
        if (lid > scene._lights.size() - 1) lid = (int)scene._lights.size() - 1;
        nlweight = scene._lights.size();
    } else {
        nlweight = 1;
    }

    // sample light according to area
    auto lpt = _sample_light(scene, lid, pt, _sample_next1f(sampler),
                             _sample_next2f(sampler));
    auto lld = _eval_emission(lpt) * _eval_brdfcos(pt, -lpt.wo);
    auto lweight = _weight_light(scene, lid, lpt, pt) * nlweight;
    lld *= lweight;
    if (lld != ym::zero3f) {
        auto shadow_ray = _offset_ray(scene, pt, lpt, params);
        if (scene.intersect_any(shadow_ray)) lld = ym::zero3f;
    }

    // check if mis is necessary
    if (pt.ptype == _point::type::point || pt.ptype == _point::type::line)
        return lld;

    // check if mis is necessary
    auto& light = scene._lights[lid];
    if (light.shape_id < 0) return lld;
    if (lpt.ptype == _point::type::point || lpt.ptype == _point::type::line) {
        return lld;
    }

    // sample the brdf
    auto bwi =
        _sample_brdfcos(pt, _sample_next1f(sampler), _sample_next2f(sampler));
    auto bweight = 0.0f;
    auto bld = ym::zero3f;
    auto bpt = _intersect_scene(scene, _offset_ray(scene, pt, bwi, params));
    if (lid == bpt.light_id || all_lights) {
        bweight = _weight_brdfcos(pt, bwi);
        bld = _eval_emission(bpt) * _eval_brdfcos(pt, bwi) * bweight;
    }

    // accumulate the value with mis
    if (lld != ym::zero3f) {
        auto bweight = _weight_brdfcos(pt, -lpt.wo);
        // float weight =
        //     (1 / lweight) * (1 / lweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        auto weight = (1 / lweight) / ((1 / lweight) + (1 / bweight));
        lld *= weight;
    }
    if (bld != ym::zero3f) {
        auto lweight = _weight_light(scene, lid, bpt, pt) * nlweight;
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
static inline ym::vec4f _shade_pathtrace_recd(const scene& scene,
                                              const ym::ray3f& ray,
                                              _sampler& sampler, int ray_depth,
                                              const render_params& params) {
    // scene intersection
    auto pt = _intersect_scene(scene, ray);
    if (pt.ptype == _point::type::none) return ym::zero4f;

    // init
    auto l = ym::vec4f(0, 0, 0, 1);

    // emission
    if (ray_depth == 0) l.xyz() += _eval_emission(pt);
    if (pt.ptype == _point::type::env) return l;

    // check early exit
    if (pt.kd == ym::zero3f && pt.ks == ym::zero3f) return l;

    // direct
    l.xyz() += _eval_direct(scene, 0, pt, sampler, params);

    // roussian roulette
    if (ray_depth >= params.max_depth) return l;
    auto rrweight = 1.0f;
    if (ray_depth >= params.min_depth) {
        auto rrrn = _sample_next1f(sampler);
        auto wrr = std::min(ym::mean(pt.kd) + ym::mean(pt.ks), 0.95f);
        if (rrrn >= wrr) return l;
        rrweight /= wrr;
    }

    // continue path
    auto bwi =
        _sample_brdfcos(pt, _sample_next1f(sampler), _sample_next2f(sampler));
    if (bwi == ym::zero3f) return l;
    auto bweight = _weight_brdfcos(pt, bwi);
    if (!bweight) return l;
    auto bbrdfcos = _eval_brdfcos(pt, bwi);
    if (bbrdfcos == ym::zero3f) return l;
    auto ble = _shade_pathtrace_recd(scene, _offset_ray(scene, pt, bwi, params),
                                     sampler, ray_depth + 1, params);
    l.xyz() += ble.xyz() * bbrdfcos * rrweight;

    return l;
}

//
// Shader interface for the above function.
//
static inline ym::vec4f _shade_pathtrace(const scene& scene,
                                         const ym::ray3f& ray,
                                         _sampler& sampler,
                                         const render_params& params) {
    return _shade_pathtrace_recd(scene, ray, sampler, 0, params);
}

//
// Direct illuination.
//
static inline ym::vec4f _shade_direct(const scene& scene, const ym::ray3f& ray,
                                      _sampler& sampler,
                                      const render_params& params) {
    // scene intersection
    auto pt = _intersect_scene(scene, ray);
    if (pt.ptype == _point::type::none) return ym::zero4f;

    // init
    auto l = ym::vec4f(0, 0, 0, 1);

    // emission
    l.xyz() += _eval_emission(pt);
    if (pt.ptype == _point::type::env) return l;

    // early exit
    if (pt.kd == ym::zero3f && pt.ks == ym::zero3f) return l;

    // ambient
    l.xyz() += params.amb * pt.kd;

    // direct
    for (int lid = 0; lid < scene._lights.size(); lid++) {
        l.xyz() += _eval_direct(scene, lid, pt, sampler, params);
    }

    // done
    return l;
}

//
// Eyelight for quick previewing.
//
static inline ym::vec4f _shade_eyelight(const scene& scene,
                                        const ym::ray3f& ray, _sampler& sampler,
                                        const render_params& params) {
    // intersection
    _point pt = _intersect_scene(scene, ray);
    if (pt.ptype == _point::type::none) return ym::zero4f;

    // init
    auto l = ym::vec4f(0, 0, 0, 1);

    // emission
    l.xyz() += _eval_emission(pt);
    if (pt.ptype == _point::type::env) return l;

    // brdf*light
    l.xyz() += _eval_brdfcos(pt, pt.wo) * ym::pif;

    return l;
}

//
// Shader function callback.
//
using shade_fn =
    std::function<ym::vec4f(const scene& scene, const ym::ray3f& ray,
                            _sampler& sampler, const render_params& params)>;

//
// Renders a block of pixels. Public API, see above.
//
YGL_API void trace_block(const scene& scene, int cid,
                         ym::image_view<ym::vec4f> img, int ns,
                         const ym::range2i& block, const ym::range1i& samples,
                         const render_params& params, bool accumulate) {
    auto& cam = scene.cameras[cid];
    shade_fn shade;
    switch (params.stype) {
        case stype::eyelight: shade = _shade_eyelight; break;
        case stype::def:
        case stype::direct: shade = _shade_direct; break;
        case stype::pathtrace: shade = _shade_pathtrace; break;
        default: assert(false); return;
    }
    for (auto j = block.min[1]; j < block.max[1]; j++) {
        for (auto i = block.min[0]; i < block.max[0]; i++) {
            auto ij = ym::vec2i(i, j);
            auto saved = img[ij];
            img[ij] = ym::zero4f;
            for (auto s = samples.min; s < samples.max; s++) {
                auto sampler = _make_sampler(i, j, s, ns, params.rtype);
                auto rn = _sample_next2f(sampler);
                auto uv = ym::vec2f{(i + rn[0]) / img.size()[0],
                                    1 - (j + rn[1]) / img.size()[1]};
                auto ray = _eval_camera(cam, uv, _sample_next2f(sampler));
                auto l = shade(scene, ray, sampler, params);
                if (!ym::isfinite(l)) continue;
                if (params.pixel_clamp > 0)
                    l.xyz() = ym::clamplen(l.xyz(), params.pixel_clamp);
                img[ij] += l;
            }
            if (accumulate && samples.min) {
                img[ij] += saved * samples.min;
                img[ij] /= samples.max;
            } else {
                img[ij] /= (float)samples.size();
            }
        }
    }
}

//
// Renders the whole image. Public API, see above.
//
YGL_API void trace_image(const scene& scene, int cid,
                         ym::image_view<ym::vec4f> img, int ns,
                         const render_params& params) {
    trace_block(scene, cid, img, ns, {{0, 0}, {img.size()[0], img.size()[1]}},
                {0, ns}, params);
}

}  // namespace

#endif

#endif
