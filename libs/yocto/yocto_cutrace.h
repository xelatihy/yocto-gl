//
// # Yocto/CuTrace: Path tracing on Cuda/Optix
//
// Yocto/CuTrace is a simple path tracer written on the Yocto/Scene model.
// Yocto/CuTrace is implemented in `yocto_cutrace.h`, `yocto_cutrace.cpp`,
// and `yocto_cutrace.cu`.
//
// THIS IS AN EXPERIMENTAL LIBRARY THAT IS NOT READY FOR PRIME TIME
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

#ifndef _YOCTO_CUTRACE_H_
#define _YOCTO_CUTRACE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "yocto_scene.h"
#include "yocto_trace.h"

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
enum struct cutrace_sampler_type {
  path,        // path tracing
  pathdirect,  // path tracing with direct
  pathmis,     // path tracing with mis
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  eyelightao,  // eyelight with ambient occlusion
  furnace,     // furnace test
  falsecolor,  // false color rendering
};
// Type of false color visualization
enum struct cutrace_falsecolor_type {
  // clang-format off
  position, normal, frontfacing, gnormal, gfrontfacing, texcoord, mtype, color,
  emission, roughness, opacity, metallic, delta, instance, shape, material, 
  element, highlight
  // clang-format on
};

// Default trace seed
const auto cutrace_default_seed = 961748941ull;

// Options for trace functions
struct cutrace_params {
  int                     camera         = 0;
  int                     resolution     = 1280;
  cutrace_sampler_type    sampler        = cutrace_sampler_type::path;
  cutrace_falsecolor_type falsecolor     = cutrace_falsecolor_type::color;
  int                     samples        = 512;
  int                     bounces        = 8;
  float                   clamp          = 10;
  bool                    nocaustics     = false;
  bool                    envhidden      = false;
  bool                    tentfilter     = false;
  uint64_t                seed           = cutrace_default_seed;
  bool                    embreebvh      = false;
  bool                    highqualitybvh = false;
  bool                    noparallel     = false;
  int                     pratio         = 8;
  float                   exposure       = 0;
  bool                    filmic         = false;
  bool                    denoise        = false;
  int                     batch          = 1;
};

// Progressively computes an image.
image_data cutrace_image(const scene_data& scene, const cutrace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// forward declarations
struct cutrace_scene;
struct cutrace_sceneext;
struct cubvh_data;
struct cutrace_state;
struct cutrace_lights;
struct cutrace_context;

// Initialize GPU context.
cutrace_context make_cutrace_context(const cutrace_params& params);

// Upload the scene to the GPU.
cutrace_sceneext make_cutrace_scene(
    const scene_data& scene, const cutrace_params& params);

// Build the bvh acceleration structure.
cubvh_data make_cutrace_bvh(cutrace_context& context, cutrace_sceneext& cuscene,
    const scene_data& scene, const cutrace_params& params);

// Initialize state.
cutrace_state make_cutrace_state(
    const scene_data& scene, const cutrace_params& params);

// Initialize lights.
cutrace_lights make_cutrace_lights(
    const scene_data& scene, const cutrace_params& params);

// Start rendering an image.
void trace_start(cutrace_context& context, cutrace_state& state,
    const cutrace_scene& cuscene, const cubvh_data& bvh,
    const cutrace_lights& lights, const scene_data& scene,
    const cutrace_params& params);

// Progressively computes an image.
void trace_samples(cutrace_context& context, cutrace_state& state,
    const cutrace_scene& cuscene, const cubvh_data& bvh,
    const cutrace_lights& lights, const scene_data& scene,
    const cutrace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENUM LABELS
// -----------------------------------------------------------------------------
namespace yocto {

// trace sampler names
inline const auto cutrace_sampler_names = vector<string>{"path", "pathdirect",
    "pathmis", "naive", "eyelight", "eyelightao", "furnace", "falsecolor"};

// false color names
inline const auto cutrace_falsecolor_names = vector<string>{"position",
    "normal", "frontfacing", "gnormal", "gfrontfacing", "texcoord", "mtype",
    "color", "emission", "roughness", "opacity", "metallic", "delta",
    "instance", "shape", "material", "element", "highlight"};

// trace sampler labels
inline const auto cutrace_sampler_labels =
    vector<pair<cutrace_sampler_type, string>>{
        {cutrace_sampler_type::path, "path"},
        {cutrace_sampler_type::pathdirect, "pathdirect"},
        {cutrace_sampler_type::pathmis, "pathmis"},
        {cutrace_sampler_type::naive, "naive"},
        {cutrace_sampler_type::eyelight, "eyelight"},
        {cutrace_sampler_type::eyelightao, "eyelightao"},
        {cutrace_sampler_type::furnace, "furnace"},
        {cutrace_sampler_type::falsecolor, "falsecolor"}};

// false color labels
inline const auto cutrace_falsecolor_labels =
    vector<pair<cutrace_falsecolor_type, string>>{
        {cutrace_falsecolor_type::position, "position"},
        {cutrace_falsecolor_type::normal, "normal"},
        {cutrace_falsecolor_type::frontfacing, "frontfacing"},
        {cutrace_falsecolor_type::gnormal, "gnormal"},
        {cutrace_falsecolor_type::gfrontfacing, "gfrontfacing"},
        {cutrace_falsecolor_type::texcoord, "texcoord"},
        {cutrace_falsecolor_type::mtype, "mtype"},
        {cutrace_falsecolor_type::color, "color"},
        {cutrace_falsecolor_type::emission, "emission"},
        {cutrace_falsecolor_type::roughness, "roughness"},
        {cutrace_falsecolor_type::opacity, "opacity"},
        {cutrace_falsecolor_type::metallic, "metallic"},
        {cutrace_falsecolor_type::delta, "delta"},
        {cutrace_falsecolor_type::instance, "instance"},
        {cutrace_falsecolor_type::shape, "shape"},
        {cutrace_falsecolor_type::material, "material"},
        {cutrace_falsecolor_type::element, "element"},
        {cutrace_falsecolor_type::highlight, "highlight"}};

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#if YOCTO_CUDA

// do not reorder
#include <cuda.h>
// do not reorder
#include <optix.h>

// -----------------------------------------------------------------------------
// DATA DEFINITIONS
// -----------------------------------------------------------------------------
namespace yocto {

// buffer view
template <typename T>
struct cubuffer {
  size_t      size() const { return _size; }
  CUdeviceptr device_ptr() const { return _data; }
  size_t      size_in_bytes() const { return _size * sizeof(T); }

  CUdeviceptr _data = 0;
  size_t      _size = 0;
};

// device params
struct cutrace_camera {
  frame3f frame;
  float   lens;
  float   film;
  float   aspect;
  float   focus;
  float   aperture;
  bool    orthographic;
};

struct cutrace_texture {
  CUarray     array;
  CUtexObject texture;
  int         width  = 0;
  int         height = 0;
  bool        linear = false;
};

struct cutrace_material {
  material_type type         = material_type::matte;
  vec3f         emission     = {0, 0, 0};
  vec3f         color        = {0, 0, 0};
  float         roughness    = 0;
  float         metallic     = 0;
  float         ior          = 1.5f;
  vec3f         scattering   = {0, 0, 0};
  float         scanisotropy = 0;
  float         trdepth      = 0.01f;
  float         opacity      = 1;

  int emission_tex   = invalidid;
  int color_tex      = invalidid;
  int roughness_tex  = invalidid;
  int scattering_tex = invalidid;
  int normal_tex     = invalidid;
};

struct cutrace_instance {
  frame3f frame;
  int     shape;
  int     material;
};

struct cutrace_shape {
  cubuffer<vec3f> positions = {};
  cubuffer<vec3f> normals   = {};
  cubuffer<vec2f> texcoords = {};
  cubuffer<vec4f> colors    = {};
  cubuffer<vec3i> triangles = {};
};

struct cutrace_environment {
  frame3f frame        = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
  vec3f   emission     = {0, 0, 0};
  int     emission_tex = invalidid;
};

struct cutrace_scene {
  cubuffer<cutrace_camera>      cameras      = {};
  cubuffer<cutrace_texture>     textures     = {};
  cubuffer<cutrace_material>    materials    = {};
  cubuffer<cutrace_shape>       shapes       = {};
  cubuffer<cutrace_instance>    instances    = {};
  cubuffer<cutrace_environment> environments = {};
};

struct cutrace_sceneext : cutrace_scene {
  vector<cutrace_texture> cutextures = {};
  vector<cutrace_shape>   cushapes   = {};
};

struct cubvh_tree {
  cubuffer<byte>         buffer = {};
  OptixTraversableHandle handle;
};

struct cubvh_data {
  cubuffer<OptixInstance> instances = {};
  cubvh_tree              instances_bvh;
  vector<cubvh_tree>      shapes_bvhs;
};

// state
struct cutrace_state {
  int                 width   = 0;
  int                 height  = 0;
  int                 samples = 0;
  cubuffer<vec4f>     image   = {};
  cubuffer<vec3f>     albedo  = {};
  cubuffer<vec3f>     normal  = {};
  cubuffer<int>       hits    = {};
  cubuffer<rng_state> rngs    = {};
  cubuffer<vec4f>     display = {};
};

// params
struct cutrace_dparams {
  int                     camera         = 0;
  int                     resolution     = 1280;
  cutrace_sampler_type    sampler        = cutrace_sampler_type::path;
  cutrace_falsecolor_type falsecolor     = cutrace_falsecolor_type::color;
  int                     samples        = 512;
  int                     bounces        = 8;
  float                   clamp          = 10;
  bool                    nocaustics     = false;
  bool                    envhidden      = false;
  bool                    tentfilter     = false;
  uint64_t                seed           = cutrace_default_seed;
  bool                    embreebvh      = false;
  bool                    highqualitybvh = false;
  bool                    noparallel     = false;
  int                     pratio         = 8;
  float                   exposure       = 0;
  bool                    filmic         = false;
  bool                    denoise        = false;
  int                     batch          = 1;
};

// light
struct cutrace_light {
  int             instance     = invalidid;
  int             environment  = invalidid;
  cubuffer<float> elements_cdf = {};
};

// lights
struct cutrace_lights {
  cubuffer<cutrace_light> lights = {};
};

// device params
struct cutrace_globals {
  cutrace_state          state  = {};
  cutrace_scene          scene  = {};
  OptixTraversableHandle bvh    = {};
  cutrace_lights         lights = {};
  cutrace_dparams        params = {};
};

// empty stb record
struct __declspec(align(OPTIX_SBT_RECORD_ALIGNMENT)) cutrace_stbrecord {
  __declspec(align(
      OPTIX_SBT_RECORD_ALIGNMENT)) char header[OPTIX_SBT_RECORD_HEADER_SIZE];
};

struct cutrace_context {
  // context
  CUcontext          cuda_context  = nullptr;
  CUstream           cuda_stream   = nullptr;
  OptixDeviceContext optix_context = nullptr;

  // pipeline
  OptixPipeline optix_pipeline = nullptr;
  OptixModule   optix_module   = nullptr;

  // programs
  OptixProgramGroup raygen_program   = nullptr;
  OptixProgramGroup miss_program     = nullptr;
  OptixProgramGroup hitgroup_program = nullptr;

  // stb
  cubuffer<cutrace_stbrecord> raygen_records   = {};
  cubuffer<cutrace_stbrecord> miss_records     = {};
  cubuffer<cutrace_stbrecord> hitgroup_records = {};
  OptixShaderBindingTable     binding_table    = {};

  // global buffer
  cubuffer<cutrace_globals> globals_buffer = {};
};

}  // namespace yocto

#else

#endif

#endif
