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

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::trace {

// Namespace aliases
namespace trc = yocto::trace;
namespace img = yocto::image;

// Math defitions
using math::bbox3f;
using math::byte;
using math::frame3f;
using math::identity3x4f;
using math::ray3f;
using math::rng_state;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// HIGH LEVEL API
// -----------------------------------------------------------------------------
namespace yocto::trace {

// Trace scene
struct scene;
struct camera;
struct environment;
struct shape;
struct texture;
struct material;
struct instance;
struct object;

// Add scene elements
trc::camera*      add_camera(trc::scene* scene);
trc::object*      add_object(trc::scene* scene);
trc::texture*     add_texture(trc::scene* scene);
trc::material*    add_material(trc::scene* scene);
trc::shape*       add_shape(trc::scene* scene);
trc::instance*    add_instance(trc::scene* scene);
trc::environment* add_environment(trc::scene* scene);

// camera properties
void set_frame(trc::camera* camera, const frame3f& frame);
void set_lens(trc::camera* camera, float lens, float aspect, float film,
    bool ortho = false);
void set_focus(trc::camera* camera, float aperture, float focus);

// object properties
void set_frame(trc::object* object, const frame3f& frame);
void set_material(trc::object* object, trc::material* material);
void set_shape(trc::object* object, trc::shape* shape);
void set_instance(trc::object* object, trc::instance* instance);

// texture properties
void set_texture(trc::texture* texture, const img::image<vec3b>& img);
void set_texture(trc::texture* texture, const img::image<vec3f>& img);
void set_texture(trc::texture* texture, const img::image<byte>& img);
void set_texture(trc::texture* texture, const img::image<float>& img);

// material properties
void set_emission(trc::material* material, const vec3f& emission,
    trc::texture* emission_tex = nullptr);
void set_color(trc::material* material, const vec3f& color,
    trc::texture* color_tex = nullptr);
void set_specular(trc::material* material, float specular = 1,
    trc::texture* specular_tex = nullptr);
void set_ior(trc::material* material, float ior);
void set_metallic(trc::material* material, float metallic,
    trc::texture* metallic_tex = nullptr);
void set_transmission(trc::material* material, float transmission, bool thin,
    float trdepth, trc::texture* transmission_tex = nullptr);
void set_translucency(trc::material* material, float translucency, bool thin,
    float trdepth, trc::texture* translucency_tex = nullptr);
void set_roughness(trc::material* material, float roughness,
    trc::texture* roughness_tex = nullptr);
void set_opacity(trc::material* material, float opacity,
    trc::texture* opacity_tex = nullptr);
void set_thin(trc::material* material, bool thin);
void set_scattering(trc::material* material, const vec3f& scattering,
    float scanisotropy, trc::texture* scattering_tex = nullptr);
void set_normalmap(trc::material* material, trc::texture* normal_tex);

// shape properties
void set_points(trc::shape* shape, const std::vector<int>& points);
void set_lines(trc::shape* shape, const std::vector<vec2i>& lines);
void set_triangles(trc::shape* shape, const std::vector<vec3i>& triangles);
void set_quads(trc::shape* shape, const std::vector<vec4i>& quads);
void set_positions(trc::shape* shape, const std::vector<vec3f>& positions);
void set_normals(trc::shape* shape, const std::vector<vec3f>& normals);
void set_texcoords(trc::shape* shape, const std::vector<vec2f>& texcoords);
void set_colors(trc::shape* shape, const std::vector<vec3f>& colors);
void set_radius(trc::shape* shape, const std::vector<float>& radius);
void set_tangents(trc::shape* shape, const std::vector<vec4f>& tangents);

// instance properties
void set_frames(trc::instance* instance, const std::vector<frame3f>& frames);

// environment properties
void set_frame(trc::environment* environment, const frame3f& frame);
void set_emission(trc::environment* environment, const vec3f& emission,
    trc::texture* emission_tex = nullptr);

// Type of tracing algorithm
enum struct sampler_type {
  path,        // path tracing
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  falsecolor,  // false color rendering
};
// Type of false color visualization
enum struct falsecolor_type {
  // clang-format off
  normal, frontfacing, gnormal, gfrontfacing, texcoord, color, emission,    
  diffuse, specular, coat, metal, transmission, translucency, refraction, 
  roughness, opacity, ior, object, element, highlight
  // clang-format on
};
// Strategy used to build the bvh
enum struct bvh_type {
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
const auto default_seed = 961748941ull;

// Options for trace functions
struct trace_params {
  int             resolution = 1280;
  sampler_type    sampler    = sampler_type::path;
  falsecolor_type falsecolor = falsecolor_type::diffuse;
  int             samples    = 512;
  int             bounces    = 8;
  float           clamp      = 100;
  bool            nocaustics = false;
  bool            envhidden  = false;
  bool            tentfilter = false;
  uint64_t        seed       = default_seed;
  bvh_type        bvh        = bvh_type::default_;
  bool            noparallel = false;
  int             pratio     = 8;
  float           exposure   = 0;
};

const auto sampler_names = std::vector<std::string>{
    "path", "naive", "eyelight", "falsecolor"};

const auto falsecolor_names = std::vector<std::string>{"normal", "frontfacing",
    "gnormal", "gfrontfacing", "texcoord", "color", "emission", "diffuse",
    "specular", "coat", "metal", "transmission", "translucency", "refraction",
    "roughness", "opacity", "ior", "object", "element", "highlight"};
const auto bvh_names        = std::vector<std::string>{
    "default", "highquality", "middle", "balanced",
#ifdef YOCTO_EMBREE
    "embree-default", "embree-highquality", "embree-compact"
#endif
};

// Progress report callback
using progress_callback =
    std::function<void(const std::string& message, int current, int total)>;
// Callback used to report partially computed image
using image_callback = std::function<void(
    const img::image<vec4f>& render, int current, int total)>;

// Initialize lights.
void init_lights(trc::scene* scene, progress_callback progress_cb = {});

// Build the bvh acceleration structure.
void init_bvh(trc::scene* scene, const trace_params& params,
    progress_callback progress_cb = {});

// Refit bvh data
void update_bvh(trc::scene*            scene,
    const std::vector<trc::object*>&   updated_objects,
    const std::vector<trc::shape*>&    updated_shapes,
    const std::vector<trc::instance*>& updated_instances,
    const trace_params&                params);

// Progressively computes an image.
img::image<vec4f> trace_image(const trc::scene* scene,
    const trc::camera* camera, const trace_params& params,
    progress_callback progress_cb = {}, image_callback image_cb = {});

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

// [experimental] Callback used to report partially computed image
using async_callback = std::function<void(
    const img::image<vec4f>& render, int current, int total, const vec2i& ij)>;

// [experimental] Asynchronous interface
struct state;
void trace_start(state* state, const trc::scene* scene,
    const trc::camera* camera, const trace_params& params,
    progress_callback progress_cb = {}, image_callback image_cb = {},
    async_callback async_cb = {});
void trace_stop(state* state);

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// SCENE AND RENDERING DATA
// -----------------------------------------------------------------------------
namespace yocto::trace {

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct bvh_node {
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
struct bvh_tree {
  std::vector<bvh_node> nodes      = {};
  std::vector<vec2i>    primitives = {};
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
struct camera {
  frame3f frame        = identity3x4f;
  bool    orthographic = false;
  float   lens         = 0.050;
  vec2f   film         = {0.036, 0.024};
  float   focus        = 10000;
  float   aperture     = 0;
};

// Texture containing either an LDR or HDR image. HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB.
struct texture {
  img::image<vec3f> colorf  = {};
  img::image<vec3b> colorb  = {};
  img::image<float> scalarf = {};
  img::image<byte>  scalarb = {};
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct material {
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
  bool  thin         = false;

  // textures
  trc::texture* emission_tex     = nullptr;
  trc::texture* color_tex        = nullptr;
  trc::texture* specular_tex     = nullptr;
  trc::texture* metallic_tex     = nullptr;
  trc::texture* roughness_tex    = nullptr;
  trc::texture* transmission_tex = nullptr;
  trc::texture* translucency_tex = nullptr;
  trc::texture* spectint_tex     = nullptr;
  trc::texture* scattering_tex   = nullptr;
  trc::texture* coat_tex         = nullptr;
  trc::texture* opacity_tex      = nullptr;
  trc::texture* normal_tex       = nullptr;
};

// Shape data represented as an indexed meshes of elements.
// May contain either points, lines, triangles and quads.
// Additionally, we support faceavarying primitives where
// each verftex data has its own topology.
struct shape {
  // primitives
  std::vector<int>   points    = {};
  std::vector<vec2i> lines     = {};
  std::vector<vec3i> triangles = {};
  std::vector<vec4i> quads     = {};

  // vertex data
  std::vector<vec3f> positions = {};
  std::vector<vec3f> normals   = {};
  std::vector<vec2f> texcoords = {};
  std::vector<vec3f> colors    = {};
  std::vector<float> radius    = {};
  std::vector<vec4f> tangents  = {};

  // computed properties
  bvh_tree* bvh = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif

  // element cdf for sampling
  std::vector<float> elements_cdf = {};

  // cleanup
  ~shape();
};

// Instances.
struct instance {
  std::vector<frame3f> frames = {};
};

// Object.
struct object {
  frame3f        frame    = identity3x4f;
  trc::shape*    shape    = nullptr;
  trc::material* material = nullptr;
  trc::instance* instance = nullptr;
};

// Environment map.
struct environment {
  frame3f            frame        = identity3x4f;
  vec3f              emission     = {0, 0, 0};
  trc::texture*      emission_tex = nullptr;
  std::vector<float> texels_cdf   = {};
};

// Trace lights used during rendering. These are created automatically.
struct light {
  trc::object*      object      = nullptr;
  int               instance    = -1;
  trc::environment* environment = nullptr;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene {
  std::vector<trc::camera*>      cameras      = {};
  std::vector<trc::object*>      objects      = {};
  std::vector<trc::shape*>       shapes       = {};
  std::vector<trc::material*>    materials    = {};
  std::vector<trc::instance*>    instances    = {};
  std::vector<trc::texture*>     textures     = {};
  std::vector<trc::environment*> environments = {};

  // computed properties
  std::vector<trc::light*> lights = {};
  bvh_tree*                bvh    = nullptr;
#ifdef YOCTO_EMBREE
  RTCScene           embree_bvh       = nullptr;
  std::vector<vec2i> embree_instances = {};
#endif

  // cleanup
  ~scene();
};

// State of a pixel during tracing
struct pixel {
  vec3f     radiance = {0, 0, 0};
  int       hits     = 0;
  int       samples  = 0;
  rng_state rng      = {};
};

// [experimental] Asynchronous state
struct state {
  img::image<vec4f> render = {};
  img::image<pixel> pixels = {};
  std::future<void> worker = {};  // async
  std::atomic<bool> stop   = {};  // async
};

}  // namespace yocto::trace

// -----------------------------------------------------------------------------
// INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto::trace {

// Results of intersect functions that include hit flag, the instance id,
// the shape element id, the shape element uv and intersection distance.
// Results values are set only if hit is true.
struct intersection3f {
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
intersection3f intersect_scene_bvh(const trc::scene* scene, const ray3f& ray,
    bool find_any = false, bool non_rigid_frames = true);
intersection3f intersect_instance_bvh(const trc::object* object, int instance,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto::trace

#endif
