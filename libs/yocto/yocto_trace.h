//
// # Yocto/Trace: Path tracing
//
// Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
// Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_sampling.h"
#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::function;
using std::future;
using std::pair;
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

// Camera based on a simple lens model. The camera is placed using a frame.
// Camera projection is described in photographic terms. In particular,
// we specify film size (35mm by default), film aspect ration,
// the lens' focal length, the focus distance and the lens aperture.
// All values are in meters. Here are some common aspect ratios used in video
// and still photography.
// 3:2    on 35 mm:  0.036 x 0.024
// 16:9   on 35 mm:  0.036 x 0.02025 or 0.04267 x 0.024
// 2.35:1 on 35 mm:  0.036 x 0.01532 or 0.05640 x 0.024
// 2.39:1 on 35 mm:  0.036 x 0.01506 or 0.05736 x 0.024
// 2.4:1  on 35 mm:  0.036 x 0.015   or 0.05760 x 0.024 (approx. 2.39 : 1)
// To compute good apertures, one can use the F-stop number from photography
// and set the aperture to focal length over f-stop.
using trace_camera = scene_camera;

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
using trace_texture = scene_texture;

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
using trace_material = scene_material;

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
using trace_shape = scene_shape;

// Object.
using trace_instance = scene_instance;

// Environment map.
using trace_environment = scene_environment;

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
using trace_scene = scene_scene;

}  // namespace yocto

// -----------------------------------------------------------------------------
// RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

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
  emission, roughness, opacity, instance, shape, material, element, highlight
  // clang-format on
};

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Options for trace functions
struct trace_params {
  int                   camera     = 0;
  int                   resolution = 1280;
  trace_sampler_type    sampler    = trace_sampler_type::path;
  trace_falsecolor_type falsecolor = trace_falsecolor_type::color;
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

inline const auto trace_sampler_names = std::vector<std::string>{
    "path", "naive", "eyelight", "falsecolor", "dalbedo", "dnormal"};

inline const auto trace_falsecolor_names = vector<string>{"position", "normal",
    "frontfacing", "gnormal", "gfrontfacing", "texcoord", "color", "emission",
    "roughness", "opacity", "instance", "shape", "material", "element",
    "highlight"};
const auto        trace_bvh_names        = vector<string>{
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
    function<void(const image_data& render, int current, int total)>;

// Progressively computes an image.
image_data trace_image(const scene_scene& scene, const trace_params& params,
    const progress_callback& progress_cb = {},
    const image_callback&    image_cb    = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Scene lights used during rendering. These are created automatically.
struct trace_light {
  instance_handle    instance     = invalid_handle;
  environment_handle environment  = invalid_handle;
  vector<float>      elements_cdf = {};
};

// Scene lights
struct trace_lights {
  vector<trace_light> lights = {};
};

// Initialize lights.
trace_lights make_lights(const scene_scene& scene, const trace_params& params,
    const progress_callback& progress_cb = {});

// Define BVH
using trace_bvh = bvh_scene;

// Build the bvh acceleration structure.
trace_bvh make_bvh(const scene_scene& scene, const trace_params& params,
    const progress_callback& progress_cb = {});

// Refit bvh data
void update_bvh(trace_bvh& bvh, const scene_scene& scene,
    const vector<trace_instance*>& updated_instances,
    const vector<trace_shape*>& updated_shapes, const trace_params& params,
    const progress_callback& progress_cb = {});

// Progressively computes an image.
image_data trace_image(const scene_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const trace_params& params,
    const progress_callback& progress_cb = {},
    const image_callback&    image_cb    = {});

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// Trace state
struct trace_state {
  // final rendered image
  image_data image = {};
  // computing buffers
  int               width        = 0;
  int               height       = 0;
  vector<vec4f>     accumulation = {};
  vector<int>       samples      = {};
  vector<rng_state> rngs         = {};
};

// [experimental] Asynchronous state
struct trace_worker {
  future<void> worker = {};  // async
  atomic<bool> stop   = {};  // async
};

// [experimental] Callback used to report partially computed image
using async_callback = function<void(
    const image_data& render, int current, int total, const vec2i& ij)>;

// [experimental] Asynchronous interface
struct trace_state;
void trace_start(trace_worker& worker, trace_state& state,
    const scene_scene& scene, const trace_bvh& bvh, const trace_lights& lights,
    const trace_params& params, const progress_callback& progress_cb = {},
    const image_callback& image_cb = {}, const async_callback& async_cb = {});
void trace_stop(trace_worker& worker);

}  // namespace yocto

#endif
