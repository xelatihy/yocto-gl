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

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// HIGH LEVEL API
// -----------------------------------------------------------------------------
namespace yocto {

// Using directives
using std::function;
using std::make_unique;
using std::unique_ptr;

// Trace scene
struct trace_scene;
struct trace_camera;
struct trace_environment;
struct trace_shape;
struct trace_texture;
struct trace_material;
struct trace_instance;
struct trace_object;

// Create scene.
unique_ptr<trace_scene> make_trace_scene();

// Add scene elements
trace_camera*      add_camera(trace_scene* scene);
trace_object*      add_object(trace_scene* scene);
trace_texture*     add_texture(trace_scene* scene);
trace_material*    add_material(trace_scene* scene);
trace_shape*       add_shape(trace_scene* scene);
trace_instance*    add_instance(trace_scene* scene);
trace_environment* add_environment(trace_scene* scene);

// camera properties
void set_frame(trace_camera* camera, const frame3f& frame);
void set_lens(trace_camera* camera, float lens, float aspect, float film);
void set_focus(trace_camera* camera, float aperture, float focus);

// object properties
void set_frame(trace_object* object, const frame3f& frame);
void set_material(trace_object* object, trace_material* material);
void set_shape(trace_object* object, trace_shape* shape);
void set_instance(trace_object* object, trace_instance* instance);

// texture properties
void set_texture(trace_texture* texture, const image<vec4b>& img);
void set_texture(trace_texture* texture, const image<vec4f>& img);

// material properties
void set_emission(trace_material* material, const vec3f& emission,
    trace_texture* emission_txt = nullptr);
void set_color(trace_material* material, const vec3f& color,
    trace_texture* color_txt = nullptr);
void set_specular(trace_material* material, float specular = 1,
    trace_texture* specular_txt = nullptr);
void set_ior(trace_material* material, float ior);
void set_metallic(trace_material* material, float metallic,
    trace_texture* metallic_txt = nullptr);
void set_transmission(trace_material* material, float transmission, bool thin,
    float trdepth, trace_texture* transmission_txt = nullptr);
void set_roughness(trace_material* material, float roughness,
    trace_texture* roughness_txt = nullptr);
void set_opacity(trace_material* material, float opacity,
    trace_texture* opacity_txt = nullptr);
void set_thin(trace_material* material, bool thin);
void set_scattering(trace_material* material, const vec3f& scattering,
    float scanisotropy, trace_texture* scattering_tex = nullptr);
void set_normalmap(trace_material* material, trace_texture* normal_txt);
void set_gltftextures(trace_material* material, bool gltf_textures);

// shape properties
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

// instance properties
void set_frames(trace_instance* instance, const vector<frame3f>& frames);

// environment properties
void set_frame(trace_environment* environment, const frame3f& frame);
void set_emission(trace_environment* environment, const vec3f& emission,
    trace_texture* emission_map = nullptr);

// Type of tracing algorithm
enum struct trace_sampler_type {
  path,        // path tracing
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  falsecolor,  // false color rendering
};
// Type of false color visualization
enum struct trace_falsecolor_type {
  // clang-format off
  normal, frontfacing, gnormal, gfrontfacing, texcoord, color, emission,    
  diffuse, specular, coat, metal, transmission, refraction, roughness, 
  object, element, highlight
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
  int                   camera     = 0;
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

const auto trace_sampler_names = vector<string>{
    "path", "naive", "eyelight", "falsecolor"};

const auto trace_falsecolor_names = vector<string>{"normal", "frontfacing",
    "gnormal", "gfrontfacing", "texcoord", "color", "emission", "diffuse",
    "specular", "coat", "metal", "transmission", "refraction", "roughness",
    "object", "element", "highlight"};
const auto trace_bvh_names        = vector<string>{
    "default", "highquality", "middle", "balanced",
#ifdef YOCTO_EMBREE
    "embree-default", "embree-highquality", "embree-compact"
#endif
};

// Progress report callback
using trace_progress =
    function<void(const string& message, int current, int total)>;
// Callback used to report partially computed image
using trace_progress_image =
    function<void(const image<vec4f>& render, int current, int total)>;

// Initialize lights.
void init_lights(trace_scene* scene, trace_progress progress_cb = {});

// Build the bvh acceleration structure.
void init_bvh(trace_scene* scene, const trace_params& params,
    trace_progress progress_cb = {});

// Refit bvh data
void update_bvh(trace_scene* scene, const vector<int>& updated_instances,
    const vector<int>& updated_shapes, const trace_params& params);

// Progressively computes an image.
image<vec4f> trace_image(const trace_scene* scene, const trace_params& params,
    trace_progress       progress_cb       = {},
    trace_progress_image progress_image_cb = {});

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// [experimental] Callback used to report partially computed image
using trace_process_async = function<void(
    const image<vec4f>& render, int current, int total, const vec2i& ij)>;

// [experimental] Asynchronous interface
struct trace_state;
unique_ptr<trace_state> trace_async_start(const trace_scene* scene,
    const trace_params& params, trace_progress progress_cb = {},
    trace_progress_image progress_image_cb = {},
    trace_process_async  progress_async_cb = {});
void trace_async_start(trace_state* state, const trace_scene* scene,
    const trace_params& params, trace_progress progress_cb = {},
    trace_progress_image progress_image_cb = {},
    trace_process_async  progress_async_cb = {});
void trace_async_stop(trace_state* state);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE AND RENDERING DATA
// -----------------------------------------------------------------------------
namespace yocto {

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct trace_bvh_node {
  bbox3f bbox;
  int    start;
  short  num;
  bool   internal;
  byte   axis;
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// Application data is not stored explicitly.
struct trace_bvh {
  vector<trace_bvh_node> nodes      = {};
  vector<vec2i>          primitives = {};
};

// Camera based on a simple lens model. The camera is placed using a frame.
// Camera projection is described in photorgaphics terms. In particular,
// we specify fil size (35mm by default), the lens' focal length, the focus
// distance and the lens aperture. All values are in meters.
// Here are some common aspect ratios used in video and still photography.
// 3:2    on 35 mm:  0.036 x 0.024
// 16:9   on 35 mm:  0.036 x 0.02025 or 0.04267 x 0.024
// 2.35:1 on 35 mm:  0.036 x 0.01532 or 0.05640 x 0.024
// 2.39:1 on 35 mm:  0.036 x 0.01506 or 0.05736 x 0.024
// 2.4:1  on 35 mm:  0.036 x 0.015   or 0.05760 x 0.024 (approx. 2.39 : 1)
// To compute good apertures, one can use the F-stop number from phostography
// and set the aperture to focal_leangth/f_stop.
struct trace_camera {
  frame3f frame        = identity3x4f;
  bool    orthographic = false;
  float   lens         = 0.050;
  vec2f   film         = {0.036, 0.024};
  float   focus        = flt_max;
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
  vec3f scattering   = {0, 0, 0};
  float scanisotropy = 0;
  float trdepth      = 0.01;
  float opacity      = 1;
  bool  thin         = false;

  // textures
  trace_texture* emission_tex     = nullptr;
  trace_texture* color_tex        = nullptr;
  trace_texture* specular_tex     = nullptr;
  trace_texture* metallic_tex     = nullptr;
  trace_texture* roughness_tex    = nullptr;
  trace_texture* transmission_tex = nullptr;
  trace_texture* spectint_tex     = nullptr;
  trace_texture* scattering_tex   = nullptr;
  trace_texture* coat_tex         = nullptr;
  trace_texture* opacity_tex      = nullptr;
  trace_texture* normal_tex       = nullptr;
  bool           gltf_textures    = false;  // glTF packed textures
};

// Shape data represented as an indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support faceavarying primitives where
// each verftex data has its own topology.
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

  // computed properties
  trace_bvh* bvh = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif

  // element cdf for sampling
  vector<float> elements_cdf = {};

  // cleanup
  ~trace_shape();
};

// Instances.
struct trace_instance {
  vector<frame3f> frames = {};
};

// Object.
struct trace_object {
  frame3f         frame    = identity3x4f;
  trace_shape*    shape    = nullptr;
  trace_material* material = nullptr;
  trace_instance* instance = nullptr;
};

// Environment map.
struct trace_environment {
  frame3f        frame        = identity3x4f;
  vec3f          emission     = {0, 0, 0};
  trace_texture* emission_tex = nullptr;
  vector<float>  texels_cdf   = {};
};

// Trace lights used during rendering. These are created automatically.
struct trace_light {
  trace_object*      object      = nullptr;
  int                instance    = -1;
  trace_environment* environment = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct trace_scene {
  vector<trace_camera*>      cameras      = {};
  vector<trace_object*>      objects      = {};
  vector<trace_shape*>       shapes       = {};
  vector<trace_material*>    materials    = {};
  vector<trace_instance*>    instances    = {};
  vector<trace_texture*>     textures     = {};
  vector<trace_environment*> environments = {};

  // computed properties
  vector<trace_light*> lights = {};
  trace_bvh*           bvh    = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene      embree_bvh       = nullptr;
  vector<vec2i> embree_instances = {};
#endif

  // cleanup
  ~trace_scene();
};

// State of a pixel during tracing
struct trace_pixel {
  vec3f     radiance = zero3f;
  int       hits     = 0;
  int       samples  = 0;
  rng_state rng      = {};
};

// [experimental] Asynchronous state
struct trace_state {
  image<vec4f>       render = {};
  image<trace_pixel> pixels = {};
  std::future<void>  worker = {};  // async
  std::atomic<bool>  stop   = {};  // async
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// Results of intersect functions that include hit flag, the instance id,
// the shape element id, the shape element uv and intersection distance.
// Results values are set only if hit is true.
struct trace_intersection {
  int   object   = -1;
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
trace_intersection intersect_scene_bvh(const trace_scene* scene,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
trace_intersection intersect_instance_bvh(const trace_object* object,
    int instance, const ray3f& ray, bool find_any = false,
    bool non_rigid_frames = true);

}  // namespace yocto

#endif
