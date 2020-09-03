//
// # Yocto/Trace: Path tracing
//
// Yocto/Trace is a simple path tracer written on the Yocto/Scene model.
// Yocto/Trace is implemented in `yocto_trace.h` and `yocto_trace.cpp`.
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
struct trace_camera {
  frame3f frame        = identity3x4f;
  bool    orthographic = false;
  float   lens         = 0.050;
  float   film         = 0.036;
  float   aspect       = 1.500;
  float   focus        = 10000;
  float   aperture     = 0;
};

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
struct trace_texture {
  image<vec4f> hdr = {};
  image<vec4b> ldr = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct trace_material {
  // material
  vec3f emission     = {0, 0, 0};
  vec3f color        = {0, 0, 0};
  float specular     = 0;
  float roughness    = 0;
  float metallic     = 0;
  float ior          = 1.5;
  vec3f spectint     = {1, 1, 1};
  float coat         = 0;
  float transmission = 0;
  float translucency = 0;
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  bool  thin         = true;

  // textures
  trace_texture* emission_tex     = nullptr;
  trace_texture* color_tex        = nullptr;
  trace_texture* specular_tex     = nullptr;
  trace_texture* metallic_tex     = nullptr;
  trace_texture* roughness_tex    = nullptr;
  trace_texture* transmission_tex = nullptr;
  trace_texture* translucency_tex = nullptr;
  trace_texture* spectint_tex     = nullptr;
  trace_texture* scattering_tex   = nullptr;
  trace_texture* coat_tex         = nullptr;
  trace_texture* opacity_tex      = nullptr;
  trace_texture* normal_tex       = nullptr;
};

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support face-varying primitives where
// each vertex data has its own topology.
struct trace_shape {
  // primitives
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // face-varying primitives
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};

  // subdivision data [experimental]
  int  subdivisions = 0;
  bool catmullclark = true;
  bool smooth       = true;

  // displacement data [experimental]
  float          displacement     = 0;
  trace_texture* displacement_tex = nullptr;

  // shape is assigned at creation
  int shape_id = -1;
};

// Object.
struct trace_instance {
  frame3f         frame    = identity3x4f;
  trace_shape*    shape    = nullptr;
  trace_material* material = nullptr;

  // instance id assigned at creation
  int instance_id = -1;
};

// Environment map.
struct trace_environment {
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {0, 0, 0};
  trace_texture* emission_tex = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct trace_scene {
  // scene elements
  vector<trace_camera*>      cameras      = {};
  vector<trace_instance*>    instances    = {};
  vector<trace_environment*> environments = {};
  vector<trace_shape*>       shapes       = {};
  vector<trace_texture*>     textures     = {};
  vector<trace_material*>    materials    = {};

  // cleanup
  ~trace_scene();
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// add element to a scene
trace_camera*      add_camera(trace_scene* scene);
trace_environment* add_environment(trace_scene* scene);
trace_instance*    add_instance(trace_scene* scene);
trace_material*    add_material(trace_scene* scene);
trace_shape*       add_shape(trace_scene* scene);
trace_texture*     add_texture(trace_scene* scene);
trace_instance*    add_complete_instance(trace_scene* scene);

// set camera properties
void set_frame(trace_camera* camera, const frame3f& frame);
void set_lens(trace_camera* camera, float lens, float aspect, float film,
    bool ortho = false);
void set_focus(trace_camera* camera, float aperture, float focus);

// set instance properties
void set_frame(trace_instance* instance, const frame3f& frame);
void set_material(trace_instance* instance, trace_material* material);
void set_shape(trace_instance* instance, trace_shape* shape);

// set texture properties
void set_texture(trace_texture* texture, const image<vec4b>& img);
void set_texture(trace_texture* texture, const image<vec4f>& img);

// set material properties
void set_emission(trace_material* material, const vec3f& emission,
    trace_texture* emission_tex = nullptr);
void set_color(trace_material* material, const vec3f& color,
    trace_texture* color_tex = nullptr);
void set_specular(trace_material* material, float specular = 1,
    trace_texture* specular_tex = nullptr);
void set_ior(trace_material* material, float ior);
void set_metallic(trace_material* material, float metallic,
    trace_texture* metallic_tex = nullptr);
void set_transmission(trace_material* material, float transmission, bool thin,
    float trdepth, trace_texture* transmission_tex = nullptr);
void set_translucency(trace_material* material, float translucency, bool thin,
    float trdepth, trace_texture* translucency_tex = nullptr);
void set_roughness(trace_material* material, float roughness,
    trace_texture* roughness_tex = nullptr);
void set_opacity(trace_material* material, float opacity,
    trace_texture* opacity_tex = nullptr);
void set_thin(trace_material* material, bool thin);
void set_scattering(trace_material* material, const vec3f& scattering,
    float scanisotropy, trace_texture* scattering_tex = nullptr);
void set_normalmap(trace_material* material, trace_texture* normal_tex);

// set shape properties
void set_points(trace_shape* shape, const vector<int>& points);
void set_lines(trace_shape* shape, const vector<vec2i>& lines);
void set_triangles(trace_shape* shape, const vector<vec3i>& triangles);
void set_quads(trace_shape* shape, const vector<vec4i>& quads);
void set_fvquads(trace_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord);
void set_positions(trace_shape* shape, const vector<vec3f>& positions);
void set_normals(trace_shape* shape, const vector<vec3f>& normals);
void set_texcoords(trace_shape* shape, const vector<vec2f>& texcoords);
void set_colors(trace_shape* shape, const vector<vec4f>& colors);
void set_radius(trace_shape* shape, const vector<float>& radius);
void set_tangents(trace_shape* shape, const vector<vec4f>& tangents);
void set_subdivision(trace_shape* shape, int subdivisions, bool catmullclark,
    bool smooth = true);
void set_displacement(
    trace_shape* shape, float displacement, trace_texture* displacement_tex);

// set environment properties
void set_frame(trace_environment* environment, const frame3f& frame);
void set_emission(trace_environment* environment, const vec3f& emission,
    trace_texture* emission_tex = nullptr);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Generates a ray from a camera.
ray3f eval_camera(
    const trace_camera* camera, const vec2f& image_uv, const vec2f& lens_uv);

// Evaluates a texture
vec2i texture_size(const trace_texture* texture);
vec4f lookup_texture(
    const trace_texture* texture, const vec2i& ij, bool ldr_as_linear = false);
vec4f eval_texture(const trace_texture* texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false);

// Evaluate instance properties
vec3f eval_position(
    const trace_instance* instance, int element, const vec2f& uv);
vec3f eval_element_normal(const trace_instance* instance, int element);
vec3f eval_normal(const trace_instance* instance, int element, const vec2f& uv);
vec2f eval_texcoord(
    const trace_instance* instance, int element, const vec2f& uv);
pair<vec3f, vec3f> eval_element_tangents(
    const trace_instance* instance, int element);
vec3f eval_normalmap(
    const trace_instance* instance, int element, const vec2f& uv);
vec3f eval_shading_normal(const trace_instance* instance, int element,
    const vec2f& uv, const vec3f& outgoing);
vec4f eval_color(const trace_instance* instance, int element, const vec2f& uv);

// Environment
vec3f eval_environment(
    const trace_environment* environment, const vec3f& direction);
vec3f eval_environment(const trace_scene* scene, const vec3f& direction);

// Material sample
struct trace_material_sample {
  vec3f emission     = {0, 0, 0};
  vec3f color        = {0, 0, 0};
  float specular     = 0;
  float roughness    = 0;
  float metallic     = 0;
  float ior          = 1.5;
  vec3f spectint     = {1, 1, 1};
  float coat         = 0;
  float transmission = 0;
  float translucency = 0;
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  bool  thin         = true;
  vec3f normalmap    = {0, 0, 1};
};

// Evaluates material and textures
trace_material_sample eval_material(
    const trace_material* material, const vec2f& texcoord);

// Material Bsdf parameters
struct trace_bsdf {
  // brdf lobes
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f metal        = {0, 0, 0};
  vec3f coat         = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  vec3f translucency = {0, 0, 0};
  vec3f refraction   = {0, 0, 0};
  float roughness    = 0;
  float ior          = 1;
  vec3f meta         = {0, 0, 0};
  vec3f metak        = {0, 0, 0};
  // weights
  float diffuse_pdf      = 0;
  float specular_pdf     = 0;
  float metal_pdf        = 0;
  float coat_pdf         = 0;
  float transmission_pdf = 0;
  float translucency_pdf = 0;
  float refraction_pdf   = 0;
};

// Eval material to obtain emission, brdf and opacity.
vec3f eval_emission(const trace_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
// Eval material to obatain emission, brdf and opacity.
trace_bsdf eval_bsdf(const trace_instance* instance, int element,
    const vec2f& uv, const vec3f& normal, const vec3f& outgoing);
float eval_opacity(const trace_instance* instance, int element, const vec2f& uv,
    const vec3f& normal, const vec3f& outgoing);
// check if a brdf is a delta
bool is_delta(const trace_bsdf& bsdf);

// Material volume parameters
struct trace_vsdf {
  vec3f density    = {0, 0, 0};
  vec3f scatter    = {0, 0, 0};
  float anisotropy = 0;
};

// check if we have a volume
bool has_volume(const trace_instance* instance);
// evaluate volume
trace_vsdf eval_vsdf(
    const trace_instance* instance, int element, const vec2f& uv);

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
  emission, diffuse, specular, coat, metal, transmission, translucency,
  refraction, roughness, opacity, ior, instance, element, highlight
  // clang-format on
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
const auto trace_bvh_names        = vector<string>{
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

// Apply subdivision and displacement rules.
void tesselate_shapes(
    trace_scene* scene, const progress_callback& progress_cb = {});
void tesselate_shape(trace_scene* shape);

// Progressively computes an image.
image<vec4f> trace_image(const trace_scene* scene, const trace_camera* camera,
    const trace_params& params, const progress_callback& progress_cb = {},
    const image_callback& image_cb = {});

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// Scene lights used during rendering. These are created automatically.
struct trace_light {
  trace_instance*    instance     = nullptr;
  trace_environment* environment  = nullptr;
  vector<float>      elements_cdf = {};
};

// Scene lights
struct trace_lights {
  // light elements
  vector<trace_light*> lights = {};

  // cleanup
  ~trace_lights();
};

// Initialize lights.
void init_lights(trace_lights* lights, const trace_scene* scene,
    const trace_params& params, const progress_callback& progress_cb = {});

// Define BVH
using trace_bvh = bvh_scene;

// Build the bvh acceleration structure.
void init_bvh(trace_bvh* bvh, const trace_scene* scene,
    const trace_params& params, const progress_callback& progress_cb = {});

// Refit bvh data
void update_bvh(trace_bvh* bvh, const trace_scene* scene,
    const vector<trace_instance*>& updated_instances,
    const vector<trace_shape*>& updated_shapes, const trace_params& params);

// Progressively computes an image.
image<vec4f> trace_image(const trace_scene* scene, const trace_camera* camera,
    const trace_bvh* bvh, const trace_lights* lights,
    const trace_params& params, const progress_callback& progress_cb = {},
    const image_callback& image_cb = {});

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
void trace_start(trace_state* state, const trace_scene* scene,
    const trace_camera* camera, const trace_bvh* bvh,
    const trace_lights* lights, const trace_params& params,
    const progress_callback& progress_cb = {},
    const image_callback& image_cb = {}, const async_callback& async_cb = {});
void trace_stop(trace_state* state);

}  // namespace yocto

#endif
