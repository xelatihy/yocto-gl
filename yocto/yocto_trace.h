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
// 1. build the ray-tracing acceleration structure with Yocto/Bvh
// 2. prepare lights for rendering `update_lights()`
// 3. define rendering params with the `trace_params` structure
// 4. initialize buffers of pixel data with `make_trace_pixels()`
// 5. render blocks of samples with `trace_samples()`
//
// The code can also run in fully asynchronous mode to preview images in a
// window.
//
// 1. build the ray-tracing acceleration structure with `make_bvh()`
// 2. prepare lights for rendering `update_lights()`
// 3. define rendering params with the `trace_params` structure
// 4. initialize buffers of pixel data with `make_trace_pixels()`
// 5. start the progressive renderer with `trace_async_start()`
// 6. stop the progressive renderer with `trace_async_stop()`
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

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// #codegen begin refl-trace

// Type of rendering algorithm.
enum struct trace_type {
    pathtrace = 0,
    eyelight,
    direct,
    pathtrace_nomis,
    debug_normal,
    debug_albedo,
    debug_texcoord,
    debug_frontfacing
};

// Rendering params.
struct trace_params {
    int resolution = 512;                       // image vertical resolution
    int nsamples = 256;                         // number of samples
    trace_type tracer = trace_type::pathtrace;  // trace type
    bool notransmission = false;    // whether to test transmission in shadows
    bool double_sided = false;      // force double sided rendering
    vec3f ambient = {0, 0, 0};      // ambient lighting
    bool envmap_invisible = false;  // view environment map
    int min_depth = 3;              // minimum ray depth
    int max_depth = 8;              // maximum ray depth
    float pixel_clamp = 100;        // final pixel clamping
    float ray_eps = 1e-4f;          // ray intersection epsilon
    bool parallel = true;           // parallel execution
    int seed = 0;                   // seed for the random number generators
    int preview_resolution = 64;    // preview resolution for async rendering
    int batch_size = 16;            // sample batch size
};

// #codegen end refl-trace

// Trace pixel state. Handles image accumulation and random number generation
// for uniform and stratified sequences. The members are not part of the
// the public API.
struct trace_pixel {
    // Accumulated radiance and coverage
    vec4f acc = zero4f;
    // Random number state
    rng_state rng = rng_state();
    // Pixel coordinates
    int i = 0, j = 0;
    // Number of samples computed
    int sample = 0;
};

// Trace light as either an instance or an environment.
struct trace_light {
    const instance* ist = nullptr;     // instance for the light
    const environment* env = nullptr;  // environment for the light
};

// Trace lights. Handles sampling of illumination.
struct trace_lights {
    std::vector<trace_light> lights;  // lights
    std::unordered_map<const shape*, std::vector<float>>
        shape_distribs;                // shape dist
    std::vector<float> light_distrib;  // light distribution
};

// Initialize trace pixels.
std::vector<trace_pixel> make_trace_pixels(
    const image4f& img, const trace_params& params);
// Initialize trace lights.
trace_lights make_trace_lights(const scene* scn);

// Trace the next `nsamples` samples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    int nsamples, const trace_params& params);

// Trace the next `nsamples` samples with image filtering.
void trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    std::vector<trace_pixel>& pixels, int nsamples, const trace_params& params);

// Trace the whole image.
inline void trace_image(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    std::vector<trace_pixel>& pixels, const trace_params& params,
    const std::function<void(int)>& callback) {
    for (auto& p : img.pixels) p = zero4f;
    for (auto cur_sample = 0; cur_sample < params.nsamples;
         cur_sample += params.batch_size) {
        if (callback) callback(cur_sample);
        trace_samples(scn, cam, bvh, lights, img, pixels,
            std::min(params.batch_size, params.nsamples - cur_sample), params);
    }
}

// Trace the whole image.
inline image4f trace_image(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_params& params,
    const std::function<void(int)>& callback) {
    auto img = make_image4f(
        (int)std::round(cam->aspect * params.resolution), params.resolution);
    auto pixels = make_trace_pixels(img, params);
    auto lights = make_trace_lights(scn);
    trace_image(scn, cam, bvh, lights, img, pixels, params, callback);
    return img;
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params, const std::function<void(int, int)>& callback);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

// Names of enum values.
inline const std::map<trace_type, std::string>& trace_type_names() {
    static auto names = std::map<trace_type, std::string>{
        {trace_type::pathtrace, "pathtrace"},
        {trace_type::eyelight, "eyelight"},
        {trace_type::direct, "direct"},
        {trace_type::pathtrace_nomis, "pathtrace_nomis"},
        {trace_type::debug_normal, "debug_normal"},
        {trace_type::debug_albedo, "debug_albedo"},
        {trace_type::debug_texcoord, "debug_texcoord"},
        {trace_type::debug_frontfacing, "debug_frontfacing"},
    };
    return names;
}

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
