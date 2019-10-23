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

#ifndef _YOCTO_TRACE_H_
#define _YOCTO_TRACE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"
#include "yocto_common.h"
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_random.h"
#include "yocto_sceneio.h"

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace yocto {

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
  // lobes
  vec3f emission        = {0, 0, 0};
  vec3f diffuse         = {0, 0, 0};
  vec3f specular        = {0, 0, 0};
  float roughness       = 0;
  float metallic        = 0;
  vec3f coat            = {0, 0, 0};
  vec3f transmission    = {0, 0, 0};
  vec3f voltransmission = {0, 0, 0};
  vec3f volmeanfreepath = {0, 0, 0};
  vec3f volemission     = {0, 0, 0};
  vec3f volscatter      = {0, 0, 0};
  float volanisotropy   = 0;
  float volscale        = 0.01;
  float opacity         = 1;
  bool  refract         = false;

  // textures
  int  emission_tex     = -1;
  int  diffuse_tex      = -1;
  int  specular_tex     = -1;
  int  metallic_tex     = -1;
  int  roughness_tex    = -1;
  int  transmission_tex = -1;
  int  subsurface_tex   = -1;
  int  coat_tex         = -1;
  int  opacity_tex      = -1;
  int  normal_tex       = -1;
  bool gltf_textures    = false;  // glTF packed textures
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
};

// Instance of a visible shape in the scene.
struct trace_instance {
  frame3f frame    = identity3x4f;
  int     shape    = -1;
  int     material = -1;
};

// Environment map.
struct trace_environment {
  frame3f frame        = identity3x4f;
  vec3f   emission     = {0, 0, 0};
  int     emission_tex = -1;
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct trace_scene {
  vector<trace_camera>      cameras      = {};
  vector<trace_shape>       shapes       = {};
  vector<trace_instance>    instances    = {};
  vector<trace_material>    materials    = {};
  vector<trace_texture>     textures     = {};
  vector<trace_environment> environments = {};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE CREATION
// -----------------------------------------------------------------------------
namespace yocto {

// Construct a scene from io
void make_trace_scene(trace_scene& scene, const scene_model& ioscene);

// Update a value from io
void update_camera(trace_camera& camera, const scene_camera& iocamera);
void update_texture(trace_texture& texture, const scene_texture& iotexture);
void update_material(
    trace_material& material, const scene_material& iomaterial);
void update_shape(
    trace_shape& shape, const scene_shape& ioshape, const scene_model& ioscene);
void update_instance(
    trace_instance& instance, const scene_instance& ioinstance);
void update_environment(
    trace_environment& environment, const scene_environment& ioenvironment);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace yocto {

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Trace bvh
using trace_bvh = bvh_shared_scene;

// Build/refit the bvh acceleration structure.
void make_bvh(
    bvh_scene& bvh, const trace_scene& scene, const bvh_params& params);
void update_bvh(bvh_scene& bvh, const trace_scene& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes,
    const bvh_params& params);
void make_bvh(
    bvh_shared_scene& bvh, const trace_scene& scene, const bvh_params& params);
void update_bvh(bvh_shared_scene& bvh, const trace_scene& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes,
    const bvh_params& params);

// Trace lights used during rendering.
struct trace_lights {
  vector<int>           instances        = {};
  vector<int>           environments     = {};
  vector<vector<float>> shape_cdfs       = {};
  vector<vector<float>> environment_cdfs = {};

  bool empty() const { return instances.empty() && environments.empty(); }
};

// Initialize lights.
trace_lights make_trace_lights(const trace_scene& scene);
void         make_trace_lights(trace_lights& lights, const trace_scene& scene);

// State of a pixel during tracing
struct trace_pixel {
  vec3f     radiance = zero3f;
  int       hits     = 0;
  int       samples  = 0;
  rng_state rng      = {};
};
struct trace_state {
  vec2i               image_size = {0, 0};
  vector<trace_pixel> pixels     = {};
};

// Initialize state of the renderer.
trace_state make_trace_state(
    const vec2i& image_size, uint64_t random_seed = trace_default_seed);
void  make_trace_state(trace_state& state, const vec2i& image_size,
     uint64_t random_seed = trace_default_seed);
vec2i camera_resolution(const trace_camera& camera, int resolution);

// Options for trace functions
struct trace_params {
  // clang-format off
  // Type of tracing algorithm to use
  enum struct sampler_type {
    path,        // path tracing
    naive,       // naive path tracing
    eyelight,    // eyelight rendering
    falsecolor,  // false color rendering
  };
  enum struct falsecolor_type {
    normal, frontfacing, gnormal, gfrontfacing, texcoord, color, emission,    
    diffuse, specular, transmission, roughness, material, shape, instance, 
    element, highlight };
  // clang-format on

  int             camera     = 0;
  int             resolution = 1280;
  sampler_type    sampler    = sampler_type::path;
  falsecolor_type falsecolor = falsecolor_type::diffuse;
  int             samples    = 512;
  int             bounces    = 8;
  int             batch      = 16;
  int             region     = 16;
  float           clamp      = 10;
  bool            envhidden  = false;
  bool            tentfilter = false;
  uint64_t        seed       = trace_default_seed;
  bool            noparallel = false;
};

const auto trace_sampler_names = vector<string>{
    "path", "naive", "eyelight", "falsecolor"};

const auto trace_falsecolor_names = vector<string>{"normal", "frontfacing",
    "gnormal", "gfrontfacing", "texcoord", "color", "emission", "diffuse",
    "specular", "transmission", "roughness", "material", "shape", "instance",
    "element", "highlight"};

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image(const trace_scene& scene, const trace_bvh& bvh,
    const trace_lights& lights, const trace_params& params);

// Progressively compute an image by calling trace_samples multiple times.
// Start with an empty state and then successively call this function to
// render the next batch of samples.
int trace_samples(image<vec4f>& image, trace_state& state,
    const trace_scene& scene, const trace_bvh& bvh, const trace_lights& lights,
    int current_sample, const trace_params& params);

// Progressively compute an image by calling trace_region multiple times.
// Compared to `trace_samples` this always runs serially and is helpful
// when building async applications.
void trace_region(image<vec4f>& image, trace_state& state,
    const trace_scene& scene, const trace_bvh& bvh, const trace_lights& lights,
    const image_region& region, int num_samples, const trace_params& params);

// Check is a sampler requires lights
bool is_sampler_lit(const trace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace yocto {

// Phong exponent to roughness.
float exponent_to_roughness(float n);

// Specular to fresnel eta.
vec3f              reflectivity_to_eta(const vec3f& reflectivity);
vec3f              eta_to_reflectivity(const vec3f& eta);
pair<vec3f, vec3f> reflectivity_to_eta(
    const vec3f& reflectivity, const vec3f& edge_tint);
vec3f eta_to_reflectivity(const vec3f& eta, const vec3f& etak);
vec3f eta_to_edge_tint(const vec3f& eta, const vec3f& etak);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(const vec3f& eta, float direction_cosine);
// Compute the fresnel term for metals.
vec3f fresnel_conductor(
    const vec3f& eta, const vec3f& etak, float direction_cosine);
// Schlick approximation of Fresnel term, optionally weighted by roughness;
vec3f fresnel_schlick(const vec3f& specular, float direction_cosine);
vec3f fresnel_schlick(
    const vec3f& specular, float direction_cosine, float roughness);

// Evaluates the microfacet distribution and geometric term (ggx or beckman).
float eval_microfacetD(float roughness, const vec3f& normal,
    const vec3f& half_vector, bool ggx = true);
float eval_microfacetG(float roughness, const vec3f& normal,
    const vec3f& half_vector, const vec3f& outgoing, const vec3f& incoming,
    bool ggx = true);
vec3f sample_microfacet(
    float roughness, const vec3f& normal, const vec2f& rn, bool ggx = true);
float sample_microfacet_pdf(float roughness, const vec3f& normal,
    const vec3f& half_vector, bool ggx = true);

// Evaluate and sample volume phase function.
vec3f sample_phasefunction(float vg, const vec2f& u);
float eval_phasefunction(float cos_theta, float vg);

// Get complex ior from metal names (eta, etak).
// Return zeros if not available.
pair<vec3f, vec3f> get_conductor_eta(const string& element);

// Get subsurface params
pair<vec3f, vec3f> get_subsurface_params(const string& name);

}  // namespace yocto

#endif
