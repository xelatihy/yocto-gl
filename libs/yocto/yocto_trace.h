//
// # Yocto/Trace: Path tracing
//
// Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
// Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef _YOCTO_TRACE_H_
#define _YOCTO_TRACE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <atomic>
#include <future>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_bvh.h"
#include "yocto_ebvh.h"
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_sampling.h"
#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Type of tracing algorithm
enum struct trace_sampler_type {
  path,           // path tracing
  pathdirect,     // path tracing with direct
  pathmis,        // path tracing with mis
  pathtest,       // path tracing test
  lightsampling,  // path tracing with light sampling
  naive,          // naive path tracing
  eyelight,       // eyelight rendering
  furnace,        // furnace test
  falsecolor,     // false color rendering
};
// Type of false color visualization
enum struct trace_falsecolor_type {
  // clang-format off
  position, normal, frontfacing, gnormal, gfrontfacing, texcoord, mtype, color,
  emission, roughness, opacity, metallic, delta, instance, shape, material, 
  element, highlight
  // clang-format on
};

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Options for trace functions
struct trace_params {
  int                   camera         = 0;
  int                   resolution     = 1280;
  trace_sampler_type    sampler        = trace_sampler_type::path;
  trace_falsecolor_type falsecolor     = trace_falsecolor_type::color;
  int                   samples        = 512;
  int                   bounces        = 8;
  float                 clamp          = 100;
  bool                  nocaustics     = false;
  bool                  envhidden      = false;
  bool                  tentfilter     = false;
  uint64_t              seed           = trace_default_seed;
  bool                  embreebvh      = false;
  bool                  highqualitybvh = false;
  bool                  noparallel     = false;
  int                   pratio         = 8;
  bool                  denoise        = false;
  int                   batch          = 1;
};

// Progressively computes an image.
image<vec4f> trace_image(const scene_data& scene, const trace_params& params);

// Convenience helper
image<vec4f> trace_image(const scene_data& scene,
    trace_sampler_type type = trace_sampler_type::path, int resolution = 1440,
    int samples = 256, int bounces = 8);

// Forward compatibility
inline image<vec4f> pathtrace_image(const scene_data& scene,
    trace_sampler_type type = trace_sampler_type::path, int resolution = 1440,
    int samples = 256, int bounces = 8) {
  return trace_image(scene, type, resolution, samples, bounces);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Scene lights used during rendering. These are created automatically.
struct trace_light {
  int           instance     = invalidid;
  int           environment  = invalidid;
  vector<float> elements_cdf = {};
};

// Scene lights
struct trace_lights {
  vector<trace_light> lights = {};
};

// Trace Bvh, a wrapper of a Yocto/Bvh and an Embree one
struct trace_bvh {
  scene_bvh  bvh  = {};
  scene_ebvh ebvh = {};
};

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// Trace state
struct trace_state {
  image<vec4f>     render   = {};
  image<vec3f>     albedo   = {};
  image<vec3f>     normal   = {};
  image<int>       hits     = {};
  image<rng_state> rngs     = {};
  image<vec4f>     denoised = {};
  int              samples  = 0;

  vec2i size() const { return render.size(); }
};

// Initialize state.
trace_state make_trace_state(
    const scene_data& scene, const trace_params& params);

// Initialize lights.
trace_lights make_trace_lights(
    const scene_data& scene, const trace_params& params);

// Build the bvh acceleration structure.
trace_bvh make_trace_bvh(const scene_data& scene, const trace_params& params);

// Progressively computes an image.
void trace_samples(trace_state& state, const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights,
    const trace_params& params);
void trace_sample(trace_state& state, const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights, int i, int j, int sample,
    const trace_params& params);

// Get resulting render, denoised if requested
image<vec4f> get_image(const trace_state& state);
void         get_image(image<vec4f>& image, const trace_state& state);

// Get internal images from state
image<vec4f> get_rendered_image(const trace_state& state);
void         get_rendered_image(image<vec4f>& image, const trace_state& state);
image<vec4f> get_denoised_image(const trace_state& state);
void         get_denoised_image(image<vec4f>& image, const trace_state& state);
image<vec3f> get_albedo_image(const trace_state& state);
void         get_albedo_image(image<vec3f>& image, const trace_state& state);
image<vec3f> get_normal_image(const trace_state& state);
void         get_normal_image(image<vec3f>& image, const trace_state& state);

// Denoise image
image<vec4f> denoise_image(const image<vec4f>& render,
    const image<vec3f>& albedo, const image<vec3f>& normal);
void         denoise_image(image<vec4f>& denoised, const image<vec4f>& render,
            const image<vec3f>& albedo, const image<vec3f>& normal);
void         denoise_image(vector<vec4f>& denoised, int width, int height,
            const vector<vec4f>& render, const vector<vec3f>& albedo,
            const vector<vec3f>& normal);

// Async implementation
struct trace_context {
  std::future<void> worker = {};
  std::atomic<bool> done   = false;
  std::atomic<bool> stop   = false;
};

// Trace context
trace_context make_trace_context(const trace_params& params);

// Async start
void trace_start(trace_context& context, trace_state& state,
    const scene_data& scene, const trace_bvh& bvh, const trace_lights& lights,
    const trace_params& params);

// Async cancel
void trace_cancel(trace_context& context);

// Async done
void trace_done(trace_context& context);

// Async preview
void trace_preview(image<vec4f>& image, trace_context& context,
    trace_state& state, const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const trace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENUM LABELS
// -----------------------------------------------------------------------------
namespace yocto {

// trace sampler names
inline const auto trace_sampler_names = vector<string>{"path", "pathdirect",
    "pathmis", "pathtest", "lightsampling", "naive", "eyelight", "furnace",
    "falsecolor"};

// false color names
inline const auto trace_falsecolor_names = vector<string>{"position", "normal",
    "frontfacing", "gnormal", "gfrontfacing", "texcoord", "mtype", "color",
    "emission", "roughness", "opacity", "metallic", "delta", "instance",
    "shape", "material", "element", "highlight"};

// trace sampler labels
inline const auto trace_sampler_labels =
    vector<pair<trace_sampler_type, string>>{{trace_sampler_type::path, "path"},
        {trace_sampler_type::pathdirect, "pathdirect"},
        {trace_sampler_type::pathmis, "pathmis"},
        {trace_sampler_type::pathtest, "pathtest"},
        {trace_sampler_type::lightsampling, "lightsampling"},
        {trace_sampler_type::naive, "naive"},
        {trace_sampler_type::eyelight, "eyelight"},
        {trace_sampler_type::furnace, "furnace"},
        {trace_sampler_type::falsecolor, "falsecolor"}};

// false color labels
inline const auto trace_falsecolor_labels =
    vector<pair<trace_falsecolor_type, string>>{
        {trace_falsecolor_type::position, "position"},
        {trace_falsecolor_type::normal, "normal"},
        {trace_falsecolor_type::frontfacing, "frontfacing"},
        {trace_falsecolor_type::gnormal, "gnormal"},
        {trace_falsecolor_type::gfrontfacing, "gfrontfacing"},
        {trace_falsecolor_type::texcoord, "texcoord"},
        {trace_falsecolor_type::mtype, "mtype"},
        {trace_falsecolor_type::color, "color"},
        {trace_falsecolor_type::emission, "emission"},
        {trace_falsecolor_type::roughness, "roughness"},
        {trace_falsecolor_type::opacity, "opacity"},
        {trace_falsecolor_type::metallic, "metallic"},
        {trace_falsecolor_type::delta, "delta"},
        {trace_falsecolor_type::instance, "instance"},
        {trace_falsecolor_type::shape, "shape"},
        {trace_falsecolor_type::material, "material"},
        {trace_falsecolor_type::element, "element"},
        {trace_falsecolor_type::highlight, "highlight"}};

}  // namespace yocto

#endif
