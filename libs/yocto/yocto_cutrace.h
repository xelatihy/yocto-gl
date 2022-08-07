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

#ifndef _YOCTO_CUTRACE_H_
#define _YOCTO_CUTRACE_H_

#ifdef YOCTO_CUDA

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

// Progressively computes an image.
image_data cutrace_image(const scene_data& scene, const trace_params& params);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOWER-LEVEL RENDERING API
// -----------------------------------------------------------------------------
namespace yocto {

// forward declarations
struct cuscene_data;
struct cuscene_bvh;
struct cutrace_state;
struct cutrace_lights;
struct cutrace_context;
using cutrace_bvh = cuscene_bvh;

// Initialize GPU context.
cutrace_context make_cutrace_context(const trace_params& params);

// Upload the scene to the GPU.
cuscene_data make_cutrace_scene(cutrace_context& context,
    const scene_data& scene, const trace_params& params);
void update_cutrace_cameras(cutrace_context& context, cuscene_data& cuscene,
    const scene_data& scene, const trace_params& params);

// Build the bvh acceleration structure.
cutrace_bvh make_cutrace_bvh(cutrace_context& context,
    const cuscene_data& cuscene, const trace_params& params);

// Initialize state.
cutrace_state make_cutrace_state(cutrace_context& context,
    const scene_data& scene, const trace_params& params);
void reset_cutrace_state(cutrace_context& context, cutrace_state& state,
    const scene_data& scene, const trace_params& params);

// Initialize lights.
cutrace_lights make_cutrace_lights(cutrace_context& context,
    const scene_data& scene, const trace_params& params);

// Start rendering an image.
void trace_start(cutrace_context& context, cutrace_state& state,
    const cuscene_data& cuscene, const cutrace_bvh& bvh,
    const cutrace_lights& lights, const scene_data& scene,
    const trace_params& params);

// Progressively computes an image.
void trace_samples(cutrace_context& context, cutrace_state& state,
    const cuscene_data& cuscene, const cutrace_bvh& bvh,
    const cutrace_lights& lights, const scene_data& scene,
    const trace_params& params);

void trace_preview(image_data& image, cutrace_context& context,
    cutrace_state& pstate, const cuscene_data& cuscene, const cutrace_bvh& bvh,
    const cutrace_lights& lights, const scene_data& scene,
    const trace_params& params);

// Get resulting render, denoised if requested
image_data get_image(const cutrace_state& state);
void       get_image(image_data& image, const cutrace_state& state);

// Get internal images from state
image_data get_rendered_image(const cutrace_state& state);
void       get_rendered_image(image_data& image, const cutrace_state& state);
image_data get_denoised_image(const cutrace_state& state);
void       get_denoised_image(image_data& image, const cutrace_state& state);
image_data get_albedo_image(const cutrace_state& state);
void       get_albedo_image(image_data& image, const cutrace_state& state);
image_data get_normal_image(const cutrace_state& state);
void       get_normal_image(image_data& image, const cutrace_state& state);

// denoise image
void denoise_image(cutrace_context& context, cutrace_state& state);

// check if display
bool is_display(const cutrace_context& context);

}  // namespace yocto

// -----------------------------------------------------------------------------
// ENUM LABELS
// -----------------------------------------------------------------------------
namespace yocto {

// trace sampler names
inline const auto cutrace_sampler_names = trace_sampler_names;

// false color names
inline const auto cutrace_falsecolor_names = trace_falsecolor_names;

// trace sampler labels
inline const auto cutrace_sampler_labels = trace_sampler_labels;

// false color labels
inline const auto cutrace_falsecolor_labels = trace_falsecolor_labels;

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// do not reorder
#include <cuda.h>
// do not reorder
#include <optix.h>

// -----------------------------------------------------------------------------
// DATA DEFINITIONS
// -----------------------------------------------------------------------------
namespace yocto {

// cuda buffer
template <typename T>
struct cuspan {
  bool        empty() const { return _size == 0; }
  size_t      size() const { return _size; }
  CUdeviceptr device_ptr() const { return _data; }
  size_t      size_in_bytes() const { return _size * sizeof(T); }
  void        swap(cuspan& other) {
           std::swap(_data, other._data);
           std::swap(_size, other._size);
  }

  CUdeviceptr _data = 0;
  size_t      _size = 0;
};

// cuda array
template <typename T>
struct cuarray {
  CUarray device_array() const { return _array; }
  CUarray _array = nullptr;
};

// device params
struct cucamera_data {
  frame3f frame        = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
  float   lens         = 0.050f;
  float   film         = 0.036f;
  float   aspect       = 1.500f;
  float   focus        = 10000;
  float   aperture     = 0;
  bool    orthographic = false;
};

struct cutexture_data {
  int         width   = 0;
  int         height  = 0;
  bool        linear  = false;
  CUtexObject texture = 0;
  CUarray     array   = nullptr;
};

struct cumaterial_data {
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

struct cuinstance_data {
  frame3f frame    = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
  int     shape    = invalidid;
  int     material = invalidid;
};

struct cushape_data {
  cuspan<vec3f> positions = {};
  cuspan<vec3f> normals   = {};
  cuspan<vec2f> texcoords = {};
  cuspan<vec4f> colors    = {};
  cuspan<vec3i> triangles = {};
};

struct cuenvironment_data {
  frame3f frame        = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
  vec3f   emission     = {0, 0, 0};
  int     emission_tex = invalidid;
};

struct cuscene_data {
  cuspan<cucamera_data>      cameras      = {};
  cuspan<cutexture_data>     textures     = {};
  cuspan<cumaterial_data>    materials    = {};
  cuspan<cushape_data>       shapes       = {};
  cuspan<cuinstance_data>    instances    = {};
  cuspan<cuenvironment_data> environments = {};

  cuscene_data() {}
  cuscene_data(cuscene_data&&);
  cuscene_data& operator=(cuscene_data&&);
  ~cuscene_data();
};

struct cubvh_tree {
  cuspan<byte>           buffer = {};
  OptixTraversableHandle handle = 0;

  cubvh_tree() {}
  cubvh_tree(cubvh_tree&&);
  cubvh_tree& operator=(cubvh_tree&&);
  ~cubvh_tree();
};

struct cushape_bvh {
  cubvh_tree bvh = {};

  cushape_bvh() {}
  cushape_bvh(cushape_bvh&&);
  cushape_bvh& operator=(cushape_bvh&&);
  ~cushape_bvh();
};

struct cuscene_bvh {
  cubvh_tree            bvh       = {};
  vector<cushape_bvh>   shapes    = {};
  cuspan<OptixInstance> instances = {};

  cuscene_bvh() {}
  cuscene_bvh(cuscene_bvh&&);
  cuscene_bvh& operator=(cuscene_bvh&&);
  ~cuscene_bvh();
};

// state
struct cutrace_state {
  int               width            = 0;
  int               height           = 0;
  int               samples          = 0;
  cuspan<vec4f>     image            = {};
  cuspan<vec3f>     albedo           = {};
  cuspan<vec3f>     normal           = {};
  cuspan<int>       hits             = {};
  cuspan<rng_state> rngs             = {};
  cuspan<vec4f>     denoised         = {};
  cuspan<byte>      denoiser_state   = {};
  cuspan<byte>      denoiser_scratch = {};

  cutrace_state() {}
  cutrace_state(cutrace_state&&);
  cutrace_state& operator=(cutrace_state&&);
  ~cutrace_state();
};

// light
struct cutrace_light {
  int           instance     = invalidid;
  int           environment  = invalidid;
  cuspan<float> elements_cdf = {};
};

// lights
struct cutrace_lights {
  cuspan<cutrace_light> lights = {};

  cutrace_lights() {}
  cutrace_lights(cutrace_lights&&);
  cutrace_lights& operator=(cutrace_lights&&);
  ~cutrace_lights();
};

// device params
struct cutrace_globals {
  cutrace_state          state  = {};
  cuscene_data           scene  = {};
  OptixTraversableHandle bvh    = {};
  cutrace_lights         lights = {};
  trace_params           params = {};
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
  cuspan<cutrace_stbrecord> raygen_records   = {};
  cuspan<cutrace_stbrecord> miss_records     = {};
  cuspan<cutrace_stbrecord> hitgroup_records = {};
  OptixShaderBindingTable   binding_table    = {};

  // global buffer
  cuspan<cutrace_globals> globals_buffer = {};

  // denoiser
  OptixDenoiser denoiser = nullptr;

  cutrace_context() {}
  cutrace_context(cutrace_context&&);
  cutrace_context& operator=(cutrace_context&&);
  ~cutrace_context();
};

}  // namespace yocto

#endif

#endif
