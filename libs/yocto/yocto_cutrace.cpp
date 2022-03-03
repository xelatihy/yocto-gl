//
// Implementation for Yocto/Trace.
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
// TODO: upload
// TODO: path tracing
// TODO: environments
// TODO: interactive rendering
//

#include "yocto_cutrace.h"

#include "yocto_sampling.h"

#if YOCTO_CUDA

// do not reorder
#include <cuda.h>
// do not reorder
#include <optix.h>
#include <optix_function_table_definition.h>
#include <optix_stubs.h>

// -----------------------------------------------------------------------------
// CUDA HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

static void check_result(CUresult result) {
  if (result != CUDA_SUCCESS) {
    const char* error_name;
    cuGetErrorName(result, &error_name);
    throw std::runtime_error{"Cuda error: " + string{error_name}};
  }
}

static void check_cusync() {
  check_result(cuStreamSynchronize(nullptr));  // TODO: cuda_stream
}

static void check_result(OptixResult result) {
  if (result != OPTIX_SUCCESS) {
    throw std::runtime_error{"Optix error"};
  }
}

// buffer view
template <typename T>
struct cubuffer {
  size_t      size() const { return _size; }
  CUdeviceptr device_ptr() const { return _data; }
  size_t      size_in_bytes() const { return _size * sizeof(T); }

  CUdeviceptr _data = 0;
  size_t      _size = 0;
};

// make a buffer
template <typename T>
static cubuffer<T> make_buffer(size_t size, const T* data) {
  auto buffer  = cubuffer<T>{};
  buffer._size = size;
  check_result(cuMemAlloc(&buffer._data, buffer.size_in_bytes()));
  if (data) {
    check_result(
        cuMemcpyHtoD(buffer.device_ptr(), data, buffer.size_in_bytes()));
  }
  return buffer;
}
template <typename T>
static cubuffer<T> make_buffer(const vector<T>& data) {
  return make_buffer(data.size(), data.data());
}
template <typename T>
static cubuffer<T> make_buffer(const T& data) {
  return make_buffer(1, &data);
}

// update a buffer
template <typename T>
static void update_buffer(cubuffer<T>& buffer, size_t size, const T* data) {
  if (buffer.size() != size) throw std::runtime_error{"Cuda buffer error"};
  check_result(cuMemcpyHtoD(buffer.device_ptr(), data, buffer.size_in_bytes()));
}
template <typename T>
static void update_buffer(cubuffer<T>& buffer, const vector<T>& data) {
  return update_buffer(buffer, data.size(), data.data());
}
template <typename T>
static void update_buffer(cubuffer<T>& buffer, const T& data) {
  return update_buffer(buffer, 1, &data);
}

// download buffer
template <typename T>
static void download_buffer(
    const cubuffer<T>& buffer, size_t size, void* data) {
  if (buffer.size() != size) throw std::runtime_error{"Cuda download error"};
  check_result(cuMemcpyDtoH(data, buffer.device_ptr(), buffer.size_in_bytes()));
}
template <typename T>
static void download_buffer(const cubuffer<T>& buffer, vector<T>& data) {
  return download_buffer(buffer, data.size(), data.data());
}
template <typename T>
static void download_buffer(const cubuffer<T>& buffer, T& data) {
  return download_buffer(buffer, 1, &data);
}
template <typename T>
static vector<T> download_buffer_vector(const cubuffer<T>& buffer) {
  auto data = vector<T>(buffer.size());
  download_buffer(buffer, data.size(), data.data());
  return data;
}
template <typename T>
static T download_buffer_value(const cubuffer<T>& buffer) {
  if (buffer.size() != 1) throw std::runtime_error{"Cuda download error"};
  auto data = T{};
  download_buffer(buffer, 1, &data);
  return data;
}

// free buffer
template <typename T>
static void clear_buffer(cubuffer<T>& buffer) {
  if (buffer.device_ptr() == 0) return;
  check_result(cuMemFree(buffer.device_ptr()));
  buffer._data = 0;
  buffer._size = 0;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HACKS
// -----------------------------------------------------------------------------
namespace yocto {

extern "C" char yocto_cutrace_ptx[];

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

// device params
struct cutrace_globals {
  cutrace_state          state  = {};
  cutrace_scene          scene  = {};
  OptixTraversableHandle bvh    = {};
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

// init cuda and optix context
static cutrace_context make_cutrace_context(const cutrace_params& params) {
  // context
  auto context = cutrace_context{};

  // init cuda
  check_result(cuInit(0));
  auto device = CUdevice{0};
  check_result(cuCtxCreate(&context.cuda_context, CU_CTX_SCHED_SPIN, device));

  // init optix
  check_result(optixInit());

  // init cuda device
  check_result(cuStreamCreate(&context.cuda_stream, CU_STREAM_DEFAULT));

  // init optix device
  check_result(cuCtxGetCurrent(&context.cuda_context));
  check_result(optixDeviceContextCreate(
      context.cuda_context, 0, &context.optix_context));

  // options
  auto module_options             = OptixModuleCompileOptions{};
  module_options.maxRegisterCount = 50;
  module_options.optLevel         = OPTIX_COMPILE_OPTIMIZATION_DEFAULT;
  module_options.debugLevel       = OPTIX_COMPILE_DEBUG_LEVEL_NONE;
  auto compile_options            = OptixPipelineCompileOptions{};
  compile_options.traversableGraphFlags =
      OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_LEVEL_INSTANCING;
  compile_options.usesMotionBlur                   = false;
  compile_options.numPayloadValues                 = 2;
  compile_options.numAttributeValues               = 2;
  compile_options.exceptionFlags                   = OPTIX_EXCEPTION_FLAG_NONE;
  compile_options.pipelineLaunchParamsVariableName = "globals";
  auto link_options                                = OptixPipelineLinkOptions{};
  link_options.maxTraceDepth                       = 2;

  // optix_module
  auto log_buffer = std::array<char, 2048>();
  auto log_size   = log_buffer.size();
  auto ptx_code   = string{yocto_cutrace_ptx};
  check_result(optixModuleCreateFromPTX(context.optix_context, &module_options,
      &compile_options, ptx_code.c_str(), ptx_code.size(), log_buffer.data(),
      &log_size, &context.optix_module));

  // raygen program
  auto raygen_options                        = OptixProgramGroupOptions{};
  auto raygen_descriptor                     = OptixProgramGroupDesc{};
  raygen_descriptor.kind                     = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
  raygen_descriptor.raygen.module            = context.optix_module;
  raygen_descriptor.raygen.entryFunctionName = "__raygen__trace_pixel";
  check_result(optixProgramGroupCreate(context.optix_context,
      &raygen_descriptor, 1, &raygen_options, log_buffer.data(), &log_size,
      &context.raygen_program));

  // miss program
  auto miss_options                      = OptixProgramGroupOptions{};
  auto miss_descriptor                   = OptixProgramGroupDesc{};
  miss_descriptor.kind                   = OPTIX_PROGRAM_GROUP_KIND_MISS;
  miss_descriptor.miss.module            = context.optix_module;
  miss_descriptor.miss.entryFunctionName = "__miss__intersect_scene";
  check_result(optixProgramGroupCreate(context.optix_context, &miss_descriptor,
      1, &miss_options, log_buffer.data(), &log_size, &context.miss_program));

  // hitgroup program
  auto hitgroup_options                 = OptixProgramGroupOptions{};
  auto hitgroup_descriptor              = OptixProgramGroupDesc{};
  hitgroup_descriptor.kind              = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
  hitgroup_descriptor.hitgroup.moduleCH = context.optix_module;
  hitgroup_descriptor.hitgroup.entryFunctionNameCH =
      "__closesthit__intersect_scene";
  hitgroup_descriptor.hitgroup.moduleAH = context.optix_module;
  hitgroup_descriptor.hitgroup.entryFunctionNameAH =
      "__anyhit__intersect_scene";
  check_result(optixProgramGroupCreate(context.optix_context,
      &hitgroup_descriptor, 1, &hitgroup_options, log_buffer.data(), &log_size,
      &context.hitgroup_program));

  // pipeline
  auto program_groups = array<OptixProgramGroup, 3>{
      context.raygen_program, context.miss_program, context.hitgroup_program};
  check_result(optixPipelineCreate(context.optix_context, &compile_options,
      &link_options, program_groups.data(), (int)program_groups.size(),
      log_buffer.data(), &log_size, &context.optix_pipeline));
  check_result(optixPipelineSetStackSize(
      context.optix_pipeline, 3 * 1024, 3 * 1024, 3 * 1024, 2));

  // stb raygen
  auto raygen_record = cutrace_stbrecord{};
  check_result(
      optixSbtRecordPackHeader(context.raygen_program, &raygen_record));
  context.raygen_records             = make_buffer(raygen_record);
  context.binding_table.raygenRecord = context.raygen_records.device_ptr();

  // stb miss
  auto miss_record = cutrace_stbrecord{};
  check_result(optixSbtRecordPackHeader(context.miss_program, &miss_record));
  context.miss_records                 = make_buffer(miss_record);
  context.binding_table.missRecordBase = context.miss_records.device_ptr();
  context.binding_table.missRecordStrideInBytes = sizeof(cutrace_stbrecord);
  context.binding_table.missRecordCount         = 1;

  // stb hitgroup
  auto hitgroup_record = cutrace_stbrecord{};
  check_result(
      optixSbtRecordPackHeader(context.hitgroup_program, &hitgroup_record));
  context.hitgroup_records = make_buffer(hitgroup_record);
  context.binding_table.hitgroupRecordBase =
      context.hitgroup_records.device_ptr();
  context.binding_table.hitgroupRecordStrideInBytes = sizeof(cutrace_stbrecord);
  context.binding_table.hitgroupRecordCount         = 1;

  // globals
  context.globals_buffer = make_buffer(cutrace_globals{});

  return context;
}

// start a new render
static void render_start(cutrace_context& context, cutrace_state& state,
    const cutrace_scene& cuscene, const cubvh_data& bvh,
    const scene_data& scene, const cutrace_params& params) {
  auto globals   = cutrace_globals{};
  globals.state  = state;
  globals.scene  = cuscene;
  globals.bvh    = bvh.instances_bvh.handle;
  globals.params = (const cutrace_dparams&)params;
  update_buffer(context.globals_buffer, globals);
  // sync so we can get the frame
  check_cusync();
}

// render a batch of samples
static void render_samples(cutrace_context& context, cutrace_state& state,
    const cutrace_scene& cuscene, const cubvh_data& bvh,
    const scene_data& scene, const cutrace_params& params) {
  check_result(optixLaunch(context.optix_pipeline, context.cuda_stream,
      context.globals_buffer.device_ptr(),
      context.globals_buffer.size_in_bytes(), &context.binding_table,
      state.width, state.height, 1));
  // sync so we can get the frame
  check_cusync();
}

static cutrace_sceneext make_cutrace_scene(
    const scene_data& scene, const cutrace_params& params) {
  auto cuscene = cutrace_sceneext{};

  auto cucameras = vector<cutrace_camera>{};
  for (auto& camera : scene.cameras) {
    auto& cucamera        = cucameras.emplace_back();
    cucamera.frame        = camera.frame;
    cucamera.lens         = camera.lens;
    cucamera.aspect       = camera.aspect;
    cucamera.film         = camera.film;
    cucamera.aperture     = camera.aperture;
    cucamera.focus        = camera.focus;
    cucamera.orthographic = (int)camera.orthographic;
  }
  cuscene.cameras = make_buffer(cucameras);

  // shapes
  for (auto& shape : scene.shapes) {
    auto& cushape     = cuscene.cushapes.emplace_back();
    cushape.positions = make_buffer(shape.positions);
    cushape.triangles = make_buffer(shape.triangles);
    if (!shape.normals.empty()) cushape.normals = make_buffer(shape.normals);
    if (!shape.texcoords.empty())
      cushape.texcoords = make_buffer(shape.texcoords);
    if (!shape.colors.empty()) cushape.colors = make_buffer(shape.colors);
  }
  cuscene.shapes = make_buffer(cuscene.cushapes);

  // textures
  for (auto& texture : scene.textures) {
    auto& cutexture  = cuscene.cutextures.emplace_back();
    cutexture.width  = texture.width;
    cutexture.height = texture.height;
    cutexture.linear = texture.linear;

    auto array_descriptor        = CUDA_ARRAY_DESCRIPTOR{};
    array_descriptor.Width       = texture.width;
    array_descriptor.Height      = texture.height;
    array_descriptor.NumChannels = 4;
    array_descriptor.Format      = CU_AD_FORMAT_UNSIGNED_INT8;
    check_result(cuArrayCreate(&cutexture.array, &array_descriptor));

    auto memcpy_descriptor          = CUDA_MEMCPY2D{};
    memcpy_descriptor.dstMemoryType = CU_MEMORYTYPE_ARRAY;
    memcpy_descriptor.dstArray      = cutexture.array;
    memcpy_descriptor.dstXInBytes   = 0;
    memcpy_descriptor.dstY          = 0;
    memcpy_descriptor.srcMemoryType = CU_MEMORYTYPE_HOST;
    memcpy_descriptor.srcHost       = nullptr;
    memcpy_descriptor.srcXInBytes   = 0;
    memcpy_descriptor.srcY          = 0;
    memcpy_descriptor.srcPitch      = texture.width * 4;
    memcpy_descriptor.WidthInBytes  = texture.width * 4;
    memcpy_descriptor.Height        = texture.height;
    if (!texture.pixelsb.empty()) {
      memcpy_descriptor.srcHost = texture.pixelsb.data();
      check_result(cuMemcpy2D(&memcpy_descriptor));
    }
    if (!texture.pixelsf.empty()) {
      auto pixelsb = vector<vec4b>(texture.pixelsf.size());
      rgb_to_srgb(pixelsb, texture.pixelsf);
      memcpy_descriptor.srcHost = pixelsb.data();
      check_result(cuMemcpy2D(&memcpy_descriptor));
    }

    auto resource_descriptor               = CUDA_RESOURCE_DESC{};
    resource_descriptor.resType            = CU_RESOURCE_TYPE_ARRAY;
    resource_descriptor.res.array.hArray   = cutexture.array;
    auto texture_descriptor                = CUDA_TEXTURE_DESC{};
    texture_descriptor.addressMode[0]      = CU_TR_ADDRESS_MODE_WRAP;
    texture_descriptor.addressMode[1]      = CU_TR_ADDRESS_MODE_WRAP;
    texture_descriptor.addressMode[2]      = CU_TR_ADDRESS_MODE_WRAP;
    texture_descriptor.filterMode          = CU_TR_FILTER_MODE_LINEAR;
    texture_descriptor.flags               = CU_TRSF_NORMALIZED_COORDINATES;
    texture_descriptor.maxAnisotropy       = 1;
    texture_descriptor.maxMipmapLevelClamp = 99;
    texture_descriptor.minMipmapLevelClamp = 0;
    texture_descriptor.mipmapFilterMode    = CU_TR_FILTER_MODE_POINT;
    texture_descriptor.borderColor[0]      = 1.0f;
    texture_descriptor.borderColor[1]      = 1.0f;
    texture_descriptor.borderColor[2]      = 1.0f;
    texture_descriptor.borderColor[3]      = 1.0f;
    check_result(cuTexObjectCreate(&cutexture.texture, &resource_descriptor,
        &texture_descriptor, nullptr));
  }
  cuscene.textures = make_buffer(cuscene.cutextures);

  auto materials = vector<cutrace_material>{};
  for (auto& material : scene.materials) {
    auto& cumaterial      = materials.emplace_back();
    cumaterial.type       = material.type;
    cumaterial.emission   = material.emission;
    cumaterial.color      = material.color;
    cumaterial.roughness  = material.roughness;
    cumaterial.metallic   = material.metallic;
    cumaterial.ior        = material.ior;
    cumaterial.scattering = material.scattering;
    cumaterial.trdepth    = material.trdepth;
    cumaterial.opacity    = material.opacity;

    cumaterial.emission_tex   = material.emission_tex;
    cumaterial.color_tex      = material.color_tex;
    cumaterial.roughness_tex  = material.roughness_tex;
    cumaterial.scattering_tex = material.scattering_tex;
    cumaterial.normal_tex     = material.normal_tex;
  }
  cuscene.materials = make_buffer(materials);

  auto instances = vector<cutrace_instance>{};
  for (auto& instance : scene.instances) {
    auto& cuinstance    = instances.emplace_back();
    cuinstance.frame    = instance.frame;
    cuinstance.shape    = instance.shape;
    cuinstance.material = instance.material;
  }
  cuscene.instances = make_buffer(instances);

  auto environments = vector<cutrace_environment>{};
  for (auto& environment : scene.environments) {
    auto& cuenvironment        = environments.emplace_back();
    cuenvironment.frame        = environment.frame;
    cuenvironment.emission     = environment.emission;
    cuenvironment.emission_tex = environment.emission_tex;
  }
  cuscene.environments = make_buffer(environments);

  return cuscene;
}

static cubvh_data make_cutrace_bvh(cutrace_context& context,
    cutrace_sceneext& cuscene, const scene_data& scene,
    const cutrace_params& params) {
  auto bvh = cubvh_data{};

  // shapes
  bvh.shapes_bvhs.resize(scene.shapes.size());
  for (auto shape_id = (size_t)0; shape_id < scene.shapes.size(); shape_id++) {
    auto& shape   = scene.shapes[shape_id];
    auto& cushape = cuscene.cushapes[shape_id];

    // input
    auto built_input                       = OptixBuildInput{};
    built_input.type                       = OPTIX_BUILD_INPUT_TYPE_TRIANGLES;
    built_input.triangleArray.vertexFormat = OPTIX_VERTEX_FORMAT_FLOAT3;
    built_input.triangleArray.vertexStrideInBytes = sizeof(vec3f);
    built_input.triangleArray.numVertices         = (int)shape.positions.size();
    auto vertex_buffer                      = cushape.positions.device_ptr();
    built_input.triangleArray.vertexBuffers = &vertex_buffer;
    built_input.triangleArray.indexFormat = OPTIX_INDICES_FORMAT_UNSIGNED_INT3;
    built_input.triangleArray.indexStrideInBytes = sizeof(vec3i);
    built_input.triangleArray.numIndexTriplets   = (int)shape.triangles.size();
    auto index_buffer                       = cushape.triangles.device_ptr();
    built_input.triangleArray.indexBuffer   = index_buffer;
    auto input_flags                        = (unsigned int)0;
    built_input.triangleArray.flags         = &input_flags;
    built_input.triangleArray.numSbtRecords = 1;
    built_input.triangleArray.sbtIndexOffsetBuffer        = 0;
    built_input.triangleArray.sbtIndexOffsetSizeInBytes   = 0;
    built_input.triangleArray.sbtIndexOffsetStrideInBytes = 0;

    // setup
    auto accelerator_options       = OptixAccelBuildOptions{};
    accelerator_options.buildFlags = OPTIX_BUILD_FLAG_NONE |
                                     OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accelerator_options.motionOptions.numKeys = 1;
    accelerator_options.operation             = OPTIX_BUILD_OPERATION_BUILD;

    auto accelerator_sizes = OptixAccelBufferSizes{};
    check_result(optixAccelComputeMemoryUsage(context.optix_context,
        &accelerator_options, &built_input, (int)1, &accelerator_sizes));

    auto compacted_size_buffer = make_buffer(1, (uint64_t*)nullptr);
    auto readback_descriptor   = OptixAccelEmitDesc{};
    readback_descriptor.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    readback_descriptor.result = compacted_size_buffer.device_ptr();

    // build
    auto temporary_buffer = make_buffer(
        accelerator_sizes.tempSizeInBytes, (byte*)nullptr);
    auto bvh_buffer = make_buffer(
        accelerator_sizes.outputSizeInBytes, (byte*)nullptr);
    auto& sbvh = bvh.shapes_bvhs[shape_id];
    check_result(optixAccelBuild(context.optix_context,
        /* cuda_stream */ 0, &accelerator_options, &built_input, (int)1,
        temporary_buffer.device_ptr(), temporary_buffer.size_in_bytes(),
        bvh_buffer.device_ptr(), bvh_buffer.size_in_bytes(), &sbvh.handle,
        &readback_descriptor, 1));
    check_cusync();

    // compact
    auto compacted_size = download_buffer_value(compacted_size_buffer);
    sbvh.buffer         = make_buffer(compacted_size, (byte*)nullptr);
    check_result(optixAccelCompact(context.optix_context,
        /*cuda_stream:*/ 0, sbvh.handle, sbvh.buffer.device_ptr(),
        sbvh.buffer.size_in_bytes(), &sbvh.handle));
    check_cusync();

    // cleanup
    clear_buffer(bvh_buffer);
    clear_buffer(temporary_buffer);
    clear_buffer(compacted_size_buffer);
  }

  // instances
  {
    // upload data
    auto instances = vector<OptixInstance>(scene.instances.size());
    for (auto instance_id = 0; instance_id < scene.instances.size();
         instance_id++) {
      auto& instance   = scene.instances[instance_id];
      auto& cuinstance = instances[instance_id];
      auto  transform  = transpose(frame_to_mat(instance.frame));
      memcpy(cuinstance.transform, &transform, sizeof(float) * 12);
      cuinstance.sbtOffset         = 0;
      cuinstance.instanceId        = instance_id;
      cuinstance.traversableHandle = bvh.shapes_bvhs[instance.shape].handle;
      cuinstance.flags             = OPTIX_INSTANCE_FLAG_NONE;
      cuinstance.visibilityMask    = 0xff;
    }
    bvh.instances = make_buffer(instances);

    // config
    auto build_input                       = OptixBuildInput{};
    build_input.type                       = OPTIX_BUILD_INPUT_TYPE_INSTANCES;
    build_input.instanceArray.instances    = bvh.instances.device_ptr();
    build_input.instanceArray.numInstances = (int)scene.instances.size();

    auto accelerator_options       = OptixAccelBuildOptions{};
    accelerator_options.buildFlags = OPTIX_BUILD_FLAG_NONE |
                                     OPTIX_BUILD_FLAG_ALLOW_COMPACTION;
    accelerator_options.motionOptions.numKeys = 1;
    accelerator_options.operation             = OPTIX_BUILD_OPERATION_BUILD;

    auto accelerator_sizes = OptixAccelBufferSizes{};
    check_result(optixAccelComputeMemoryUsage(context.optix_context,
        &accelerator_options, &build_input, (int)1, &accelerator_sizes));

    auto compacted_size_buffer = make_buffer(1, (uint64_t*)nullptr);
    auto readback_descriptor   = OptixAccelEmitDesc{};
    readback_descriptor.type   = OPTIX_PROPERTY_TYPE_COMPACTED_SIZE;
    readback_descriptor.result = compacted_size_buffer.device_ptr();

    // build
    auto temporary_buffer = make_buffer(
        accelerator_sizes.tempSizeInBytes, (byte*)nullptr);
    auto bvh_buffer = make_buffer(
        accelerator_sizes.outputSizeInBytes, (byte*)nullptr);

    auto& ibvh = bvh.instances_bvh;
    check_result(optixAccelBuild(context.optix_context,
        /* cuda_stream */ 0, &accelerator_options, &build_input, (int)1,
        temporary_buffer.device_ptr(), temporary_buffer.size_in_bytes(),
        bvh_buffer.device_ptr(), bvh_buffer.size_in_bytes(), &ibvh.handle,
        &readback_descriptor, 1));
    check_cusync();

    // compact
    auto compacted_size = download_buffer_value(compacted_size_buffer);

    ibvh.buffer = make_buffer(compacted_size, (byte*)nullptr);
    check_result(optixAccelCompact(context.optix_context,
        /*cuda_stream:*/ 0, ibvh.handle, ibvh.buffer.device_ptr(),
        ibvh.buffer.size_in_bytes(), &ibvh.handle));
    check_cusync();

    // cleanup
    clear_buffer(bvh_buffer);
    clear_buffer(temporary_buffer);
    clear_buffer(compacted_size_buffer);
  }

  // done
  return bvh;
}

static cutrace_state make_cutrace_state(
    const scene_data& scene, const cutrace_params& params) {
  auto& camera = scene.cameras[params.camera];
  auto  state  = cutrace_state{};
  if (camera.aspect >= 1) {
    state.width  = params.resolution;
    state.height = (int)round(params.resolution / camera.aspect);
  } else {
    state.height = params.resolution;
    state.width  = (int)round(params.resolution * camera.aspect);
  }
  state.samples = 0;
  state.image   = make_buffer(state.width * state.height, (vec4f*)nullptr);
  state.albedo  = make_buffer(state.width * state.height, (vec3f*)nullptr);
  state.normal  = make_buffer(state.width * state.height, (vec3f*)nullptr);
  state.hits    = make_buffer(state.width * state.height, (int*)nullptr);
  state.rngs    = make_buffer(state.width * state.height, (rng_state*)nullptr);
  state.display = make_buffer(state.width * state.height, (vec4f*)nullptr);
  return state;
};

image_data cutrace_image(
    const scene_data& scene, const cutrace_params& params) {
  // initialization
  auto context = make_cutrace_context(params);
  auto cuscene = make_cutrace_scene(scene, params);
  auto bvh     = make_cutrace_bvh(context, cuscene, scene, params);
  auto state   = make_cutrace_state(scene, params);

  // rendering
  render_start(context, state, cuscene, bvh, scene, params);
  // for (auto sample = 0; sample < params.samples; sample++) {
  render_samples(context, state, cuscene, bvh, scene, params);
  // }

  // copy back image
  auto image = make_image(state.width, state.height, false);
  download_buffer(state.image, image.pixels);

  // cleanup

  // done
  return image;
}

}  // namespace yocto

#else

// -----------------------------------------------------------------------------
// NO CUDA
// -----------------------------------------------------------------------------
namespace yocto {

static void exit_nocuda() { throw std::runtime_error{"cuda not linked\n"}; }

image_data cutrace_image(
    const scene_data& scene, const cutrace_params& params) {
  exit_nocuda();
}

}  // namespace yocto

#endif