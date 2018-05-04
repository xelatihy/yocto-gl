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
// 2. create the inmage buffer and random number generators `make_rng_seq()`
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

#include <thread>
#include <unordered_map>

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Trace evaluation function.
using trace_func = vec3f (*)(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);

// Trace the next nsamples samples. Assumes that the
// image contains cur_samples already. Returns true when done.
void trace_samples(const scene* scn, const camera* cam, int width, int height,
    std::vector<vec4f>& img, std::vector<rng_state>& rngs, int cur_samples,
    int nsamples, trace_func tracer, int nbounces, float pixel_clamp = 100);
// Like before but with multiplthreading.
void trace_samples_mt(const scene* scn, const camera* cam, int width,
    int height, std::vector<vec4f>& img, std::vector<rng_state>& rngs,
    int cur_samples, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp = 100);

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, int width,
    int height, std::vector<vec4f>& img, std::vector<rng_state>& rngs,
    int nsamples, trace_func tracer, int nbounces,
    std::vector<std::thread>& threads, bool& stop_flag, int& cur_sample,
    float pixel_clamp = 100);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

// Trace function - path tracer.
vec3f trace_path(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - path tracer without mis.
vec3f trace_path_nomis(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - naive path tracer.
vec3f trace_path_naive(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - direct illumination.
vec3f trace_direct(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - direct illumination without mis.
vec3f trace_direct_nomis(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - direct illumination without mis.
vec3f trace_direct_nomis(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - eyelight rendering.
vec3f trace_eyelight(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - normal debug visualization.
vec3f trace_debug_normal(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - faceforward debug visualization.
vec3f trace_debug_frontfacing(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - albedo debug visualization.
vec3f trace_debug_albedo(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);
// Trace function - texcoord debug visualization.
vec3f trace_debug_texcoord(
    const scene* scn, const ray3f& ray, rng_state& rng, int nbounces);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);
// Compute the fresnel term for metals.
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);
// Schlick approximation of Fresnel term.
vec3f fresnel_schlick(const vec3f& ks, float cosw);
// Schlick approximation of Fresnel term weighted by roughness.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);

// Evaluates the GGX distribution and geometric term.
float eval_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
// Evaluates the GGX pdf.
float sample_ggx_pdf(float rs, float ndh);

}  // namespace ygl

#endif
