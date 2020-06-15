//
// # Yocto/Trace: Tiny path tracer
//
//
// Yocto/Trace is a simple path tracing library with support for microfacet
// materials, area and environment lights, and advacned sampling.
//
//
// ## Physically-based Path Tracing
//
// Yocto/Trace includes a tiny, but fully featured, path tracer with support for
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
// environment is used.
//
// 1. prepare the ray-tracing acceleration structure with `build_bvh()`
// 2. prepare lights for rendering with `init_trace_lights()`
// 3. create the random number generators with `init_trace_state()`
// 4. render blocks of samples with `trace_samples()`
// 5. you can also start an asynchronous renderer with `trace_asynch_start()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

// TODO: flatten state

#ifndef _YOCTO_TRACE_H_
#define _YOCTO_TRACE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <atomic>
#include <future>
#include <memory>

#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_sampling.h"
#include "yocto_scene.h"

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::function;
using std::future;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Type of tracing algorithm
enum struct trace_sampler_type {
  path,        // path tracing
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  falsecolor,  // false color rendering
  albedo,      // renders the (approximate) albedo of objects for denoising
  normal,      // renders the normals of objects for denoising
};
// Type of false color visualization
enum struct trace_falsecolor_type {
  // clang-format off
  position, normal, frontfacing, gnormal, gfrontfacing, texcoord, color, 
  emission, diffuse, specular, coat, metal, transmission, translucency, 
  refraction, roughness, opacity, ior, instance, element, highlight
  // clang-format on
};
// Strategy used to build the bvh
enum struct trace_bvh_type {
  default_,
  highquality,
  middle,
  balanced,
#ifdef YOCTO_EMBREE
  embree_default,
  embree_highquality,
  embree_compact  // only for copy interface
#endif
};

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Options for trace functions
struct trace_params {
  int                   resolution = 1280;
  trace_sampler_type    sampler    = trace_sampler_type::path;
  trace_falsecolor_type falsecolor = trace_falsecolor_type::diffuse;
  int                   samples    = 512;
  int                   bounces    = 8;
  float                 clamp      = 100;
  bool                  nocaustics = false;
  bool                  envhidden  = false;
  bool                  tentfilter = false;
  uint64_t              seed       = trace_default_seed;
  trace_bvh_type        bvh        = trace_bvh_type::default_;
  bool                  noparallel = false;
  int                   pratio     = 8;
  float                 exposure   = 0;
};

const auto trace_sampler_names = std::vector<std::string>{
    "path", "naive", "eyelight", "falsecolor", "dalbedo", "dnormal"};

const auto trace_falsecolor_names = vector<string>{"position", "normal",
    "frontfacing", "gnormal", "gfrontfacing", "texcoord", "color", "emission",
    "diffuse", "specular", "coat", "metal", "transmission", "translucency",
    "refraction", "roughness", "opacity", "ior", "instance", "element",
    "highlight"};
const auto bvh_names              = vector<string>{
    "default", "highquality", "middle", "balanced",
#ifdef YOCTO_EMBREE
    "embree-default", "embree-highquality", "embree-compact"
#endif
};

// Progress report callback
using progress_callback =
    function<void(const string& message, int current, int total)>;
// Callback used to report partially computed image
using image_callback =
    function<void(const image<vec4f>& render, int current, int total)>;

// Initialize lights.
void init_lights(scene_model* scene, progress_callback progress_cb = {});

// Build the bvh acceleration structure.
void init_bvh(scene_model* scene, const trace_params& params,
    progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(scene_model*       scene,
    const vector<scene_instance*>& updated_objects,
    const vector<scene_shape*>& updated_shapes, const trace_params& params);

// Progressively computes an image.
image<vec4f> trace_image(const scene_model* scene, const scene_camera* camera,
    const trace_params& params, progress_callback progress_cb = {},
    image_callback image_cb = {});

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// [experimental] Asynchronous state
struct trace_state {
  image<vec4f>     render       = {};
  image<vec4f>     accumulation = {};
  image<int>       samples      = {};
  image<rng_state> rngs         = {};
  future<void>     worker       = {};  // async
  atomic<bool>     stop         = {};  // async
};

// [experimental] Callback used to report partially computed image
using async_callback = function<void(
    const image<vec4f>& render, int current, int total, const vec2i& ij)>;

// [experimental] Asynchronous interface
struct trace_state;
void trace_start(trace_state* state, const scene_model* scene,
    const scene_camera* camera, const trace_params& params,
    progress_callback progress_cb = {}, image_callback image_cb = {},
    async_callback async_cb = {});
void trace_stop(trace_state* state);

}  // namespace yocto

#endif
