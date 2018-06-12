//
// # Yocto/Trace: Tiny C++ Library for Physically-based Path Tracing
//
// Yocto/Trace is a tiny, but fully featured, path tracer with support for
// textured mesh area lights, GGX materials, environment mapping. The algorithm
// makes heavy use of MIS for fast convergence.
// The interface supports progressive parallel execution both synchronously,
// for CLI applications, and asynchronously for interactive viewing.
//
// Materials are represented as sums of an emission term, a diffuse term and
// a specular microfacet term (GGX or Phong), and a transmission term for
// this sheet glass.
// Lights are defined as any shape with a material emission term. Additionally
// one can also add environment maps. But even if you can, you might want to
// add a large triangle mesh with inward normals instead. The latter is more
// general (you can even more an arbitrary shape sun). For now only the first
// env is used.
//
// # Usage
//
// 1. prepare the scene for tracing
//    - build the ray-tracing acceleration structure with `update_bvh()`
//     - prepare lights for rendering with `update_lights()`
// 2. create the inmage buffer and random number generators `make_trace_rngs()`
// 3. render blocks of samples with `trace_samples()`
// 4. you can also start an asynchronous renderer with `trace_asynch_start()`
//
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
//

#ifndef _YGL_TRACE_H_
#define _YGL_TRACE_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_scene.h"

#include <functional>
#include <thread>
#include <unordered_map>

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Trace evaluation function.
using trace_func = std::function<vec3f(const std::shared_ptr<scene>& scn,
    const ray3f& ray, rng_state& rng, int nbounces, bool* hit)>;

// Init a sequence of random number generators.
inline image<rng_state> make_trace_rngs(int w, int h, uint64_t seed) {
    auto rngs = image<rng_state>{w, h};
    int rseed = 1301081;  // large prime
    for (auto i = 0; i < w * h; i++) {
        rngs[i] = make_rng(seed, rseed + 1);
        rseed = (rseed * 1103515245 + 12345) & ((1U << 31) - 1);  // bsd rand
    }
    return rngs;
}

// Trace the next nsamples samples. Assumes that the
// image contains cur_samples already. Returns true when done.
void trace_samples(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int cur_samples, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp = 100);
// Like before but with multiplthreading.
void trace_samples_mt(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int cur_samples, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp = 100);

// Starts an anyncrhounous renderer.
void trace_async_start(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int nsamples, trace_func tracer, int nbounces,
    std::vector<std::thread>& threads, bool& stop_flag, int& cur_sample,
    float pixel_clamp = 100);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

// Trace function - path tracer.
vec3f trace_path(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - path tracer without mis.
vec3f trace_path_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - naive path tracer.
vec3f trace_path_naive(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - direct illumination.
vec3f trace_direct(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - direct illumination without mis.
vec3f trace_direct_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - pure environment illumination with no shadows.
vec3f trace_environment(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - eyelight rendering.
vec3f trace_eyelight(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - normal debug visualization.
vec3f trace_debug_normal(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - faceforward debug visualization.
vec3f trace_debug_frontfacing(const std::shared_ptr<scene>& scn,
    const ray3f& ray, rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - albedo debug visualization.
vec3f trace_debug_albedo(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - diffuse debug visualization.
vec3f trace_debug_diffuse(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - specular debug visualization.
vec3f trace_debug_specular(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - roughness debug visualization.
vec3f trace_debug_roughness(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - texcoord debug visualization.
vec3f trace_debug_texcoord(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);
float specular_to_eta(const vec3f& ks);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);
// Compute the fresnel term for metals.
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);
// Schlick approximation of Fresnel term, optionally weighted by rs;
vec3f fresnel_schlick(const vec3f& ks, float cosw);
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& o);
vec3f fresnel_schlick(
    const vec3f& ks, const vec3f& h, const vec3f& o, float rs);

// Evaluates the GGX distribution and geometric term.
float eval_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
float sample_ggx_pdf(float rs, float ndh);

// Evaluates the GGX distribution and geometric term.
float eval_ggx_dist(float rs, const vec3f& n, const vec3f& h);
float eval_ggx_sm(float rs, const vec3f& n, const vec3f& o, const vec3f& i);

}  // namespace ygl

#endif
