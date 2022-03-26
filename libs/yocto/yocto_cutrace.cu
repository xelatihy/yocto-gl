//
// # Yocto/CuTrace: Path tracing on Cuda/Optix
//
// Yocto/CuTrace is a simple path tracer written on the Yocto/Scene model.
// Yocto/CuTrace is implemented in `yocto_cutrace.h`, `yocto_cutrace.cpp`,
// and `yocto_cutrace.cu`.
// This library includes a stand-alone implementaton of the PCG32 random number
// generator by M.E. O'Neill.
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
// LICENSE OF INCLUDED SOFTWARE for Pcg random number generator
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//

#include <optix_device.h>

#include "yocto_color.h"
#include "yocto_geometry.h"
#include "yocto_math.h"
#include "yocto_sampling.h"
#include "yocto_shading.h"

// HACK TO ALLOW CUT&PASTING FROM YOCTO'S CODE
#define inline __forceinline__ __device__
#define static static __forceinline__ __device__
#define optix_shader extern "C" __global__
#define optix_constant extern "C" __constant__

// whether to use builtin compound types or yocto's ones
#define CUTRACE_BUILTIN_VECS 0

// -----------------------------------------------------------------------------
// SUBSTITUTES FOR STD TYPES
// -----------------------------------------------------------------------------
namespace yocto {

// pair
template <typename T1, typename T2>
struct pair_ {
  T1 first;
  T2 second;
};

template <typename T>
struct span {
  inline bool   empty() const { return _size == 0; }
  inline size_t size() const { return _size; }

  inline T&       operator[](int idx) { return _data[idx]; }
  inline const T& operator[](int idx) const { return _data[idx]; }
  inline T&       at(int idx) { return _data[idx]; }
  inline const T& at(int idx) const { return _data[idx]; }

  inline T*       begin() { return _data; }
  inline T*       end() { return _data + _size; }
  inline const T* begin() const { return _data; }
  inline const T* end() const { return _data + _size; }

  inline T&       front() { return *_data; }
  inline T&       back() { return *(_data + _size - 1); }
  inline const T& front() const { return *_data; }
  inline const T& back() const { return *(_data + _size - 1); }

  inline T*       data() { return _data; }
  inline const T* data() const { return _data; }

  T*     _data = nullptr;
  size_t _size = 0;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// simplified version of possible implementation from cpprenference.com
template <class T>
static const T* _upper_bound(const T* first, const T* last, const T& value) {
  const T*  it;
  ptrdiff_t count, step;
  count = last - first;

  while (count > 0) {
    it   = first;
    step = count / 2;
    it += step;
    if (!(value < *it)) {
      first = ++it;
      count -= step + 1;
    } else
      count = step;
  }
  return first;
}

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const span<float>& cdf, float r) {
  r = clamp(r * cdf.back(), (float)0, cdf.back() - (float)0.00001);
  auto idx =
      (int)(_upper_bound(cdf.data(), cdf.data() + cdf.size(), r) - cdf.data());
  return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const span<float>& cdf, int idx) {
  if (idx == 0) return cdf.at(0);
  return cdf.at(idx) - cdf.at(idx - 1);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
struct cuspan {
  inline bool     empty() const { return _size == 0; }
  inline size_t   size() const { return _size; }
  inline T&       operator[](int idx) { return _data[idx]; }
  inline const T& operator[](int idx) const { return _data[idx]; }

  inline T*       begin() { return _data; }
  inline T*       end() { return _data + _size; }
  inline const T* begin() const { return _data; }
  inline const T* end() const { return _data + _size; }

  inline T&       front() { return *_data; }
  inline T&       back() { return *(_data + _size - 1); }
  inline const T& front() const { return *_data; }
  inline const T& back() const { return *(_data + _size - 1); }

  inline operator span<T>() const { return {_data, _size}; }

  T*     _data = nullptr;
  size_t _size = 0;
};

template <typename T, size_t Size = 16>
struct svector {
  inline bool     empty() const { return _size == 0; }
  inline size_t   size() const { return _size; }
  inline T&       operator[](int idx) { return _data[idx]; }
  inline const T& operator[](int idx) const { return _data[idx]; }

  inline T*       begin() { return _data; }
  inline T*       end() { return _data + _size; }
  inline const T* begin() const { return _data; }
  inline const T* end() const { return _data + _size; }

  inline T&       front() { return *_data; }
  inline T&       back() { return *(_data + _size - 1); }
  inline const T& front() const { return *_data; }
  inline const T& back() const { return *(_data + _size - 1); }

  inline void push_back(const T& value) { _data[_size++] = value; }
  inline void pop_back() { _size--; }

  T      _data[Size] = {};
  size_t _size       = 0;
};

inline void* unpackPointer(uint32_t i0, uint32_t i1) {
  const uint64_t uptr = static_cast<uint64_t>(i0) << 32 | i1;
  void*          ptr  = reinterpret_cast<void*>(uptr);
  return ptr;
}

inline void packPointer(void* ptr, uint32_t& i0, uint32_t& i1) {
  const uint64_t uptr = reinterpret_cast<uint64_t>(ptr);
  i0                  = uptr >> 32;
  i1                  = uptr & 0x00000000ffffffff;
}

template <typename T>
inline T* getPRD() {
  const uint32_t u0 = optixGetPayload_0();
  const uint32_t u1 = optixGetPayload_1();
  return reinterpret_cast<T*>(unpackPointer(u0, u1));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUTRACE TYPES
// -----------------------------------------------------------------------------
namespace yocto {

constexpr int invalidid = -1;

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
};

struct cucamera_data {
  frame3f frame        = {};
  float   lens         = {};
  float   film         = {};
  float   aspect       = {};
  float   focus        = {};
  float   aperture     = {};
  bool    orthographic = {};
};

struct cutexture_data {
  int                 width   = 0;
  int                 height  = 0;
  bool                linear  = false;
  cudaTextureObject_t texture = 0;
  cudaArray_t         array   = nullptr;
};

enum struct material_type {
  // clang-format off
  matte, glossy, reflective, transparent, refractive, subsurface, volumetric, 
  gltfpbr
  // clang-format on
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
  int     shape    = -1;
  int     material = -1;
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
};

// Type of tracing algorithm
enum struct trace_sampler_type {
  path,        // path tracing
  pathdirect,  // path tracing with direct
  pathmis,     // path tracing with mis
  pathtest,    // path tracing for testing
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  furnace,     // furnace test
  falsecolor,  // false color rendering
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
constexpr auto trace_default_seed = 961748941ull;

// params
struct trace_params {
  int                   camera         = 0;
  int                   resolution     = 1280;
  trace_sampler_type    sampler        = trace_sampler_type::path;
  trace_falsecolor_type falsecolor     = trace_falsecolor_type::color;
  int                   samples        = 512;
  int                   bounces        = 8;
  float                 clamp          = 10;
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

using cutrace_bvh = OptixTraversableHandle;

// light
struct cutrace_light {
  int           instance     = invalidid;
  int           environment  = invalidid;
  cuspan<float> elements_cdf = {};
};

// lights
struct cutrace_lights {
  cuspan<cutrace_light> lights = {};
};

struct cutrace_globals {
  cutrace_state          state  = {};
  cuscene_data           scene  = {};
  OptixTraversableHandle bvh    = 0;
  cutrace_lights         lights = {};
  trace_params           params = {};
};

// global data
optix_constant cutrace_globals globals;

// compatibility aliases
using trace_bvh    = cutrace_bvh;
using trace_lights = cutrace_lights;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// compatibility aliases
using scene_data       = cuscene_data;
using camera_data      = cucamera_data;
using material_data    = cumaterial_data;
using texture_data     = cutexture_data;
using instance_data    = cuinstance_data;
using shape_data       = cushape_data;
using environment_data = cuenvironment_data;

// constant values
constexpr auto min_roughness = 0.03f * 0.03f;

// Evaluates an image at a point `uv`.
static vec4f eval_texture(const texture_data& texture, const vec2f& texcoord,
    bool as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  auto fromTexture = tex2D<float4>(texture.texture, texcoord.x, texcoord.y);
  auto color       = vec4f{
      fromTexture.x, fromTexture.y, fromTexture.z, fromTexture.w};
  if (as_linear && !texture.linear) {
    return srgb_to_rgb(color);
  } else {
    return color;
  }
}

// Helpers
static vec4f eval_texture(const scene_data& scene, int texture, const vec2f& uv,
    bool ldr_as_linear = false, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(
      scene.textures[texture], uv, ldr_as_linear, no_interpolation);
}

// Material parameters evaluated at a point on the surface
struct material_point {
  material_type type         = material_type::gltfpbr;
  vec3f         emission     = {0, 0, 0};
  vec3f         color        = {0, 0, 0};
  float         opacity      = 1;
  float         roughness    = 0;
  float         metallic     = 0;
  float         ior          = 1;
  vec3f         density      = {0, 0, 0};
  vec3f         scattering   = {0, 0, 0};
  float         scanisotropy = 0;
  float         trdepth      = 0.01f;
};

// Eval position
static vec3f eval_position(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_point(
        instance.frame, interpolate_triangle(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z], uv));
  } else {
    return {0, 0, 0};
  }
}

// Shape element normal.
static vec3f eval_element_normal(
    const scene_data& scene, const instance_data& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, triangle_normal(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z]));
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
static vec3f eval_normal(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.normals.empty())
    return eval_element_normal(scene, instance, element);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, normalize(interpolate_triangle(shape.normals[t.x],
                            shape.normals[t.y], shape.normals[t.z], uv)));
  } else {
    return {0, 0, 0};
  }
}

// Eval texcoord
static vec2f eval_texcoord(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.texcoords.empty()) return uv;
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.texcoords[t.x], shape.texcoords[t.y], shape.texcoords[t.z], uv);
  } else {
    return {0, 0};
  }
}

// Shape element normal.
static pair_<vec3f, vec3f> eval_element_tangents(
    const scene_data& scene, const instance_data& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty() && !shape.texcoords.empty()) {
    auto t   = shape.triangles[element];
    auto tuv = triangle_tangents_fromuv(shape.positions[t.x],
        shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
        shape.texcoords[t.y], shape.texcoords[t.z]);
    return {transform_direction(instance.frame, tuv.first),
        transform_direction(instance.frame, tuv.second)};
  } else {
    return {};
  }
}

static vec3f eval_normalmap(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalidid && (!shape.triangles.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, false));
    auto  tuv        = eval_element_tangents(scene, instance, element);
    auto  frame      = frame3f{tuv.first, tuv.second, normal, {0, 0, 0}};
    frame.x          = orthonormalize(frame.x, frame.z);
    frame.y          = normalize(cross(frame.z, frame.x));
    auto flip_v      = dot(frame.y, tuv.second) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading position
static vec3f eval_shading_position(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    return eval_position(scene, instance, element, uv);
  } else {
    return {0, 0, 0};
  }
}

// Eval shading normal
static vec3f eval_shading_normal(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalidid) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::refractive) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else {
    return {0, 0, 0};
  }
}

// Eval color
static vec4f eval_color(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.colors[t.x], shape.colors[t.y], shape.colors[t.z], uv);
  } else {
    return {0, 0, 0, 0};
  }
}

// Evaluate material
static material_point eval_material(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::refractive ||
      material.type == material_type::volumetric ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::gltfpbr ||
      point.type == material_type::glossy) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  } else if (material.type == material_type::volumetric) {
    point.roughness = 0;
  } else {
    if (point.roughness < min_roughness) point.roughness = 0;
  }

  return point;
}

static bool is_volumetric(const material_data& material) {
  return material.type == material_type::refractive ||
         material.type == material_type::volumetric ||
         material.type == material_type::subsurface;
}

// check if an instance is volumetric
static bool is_volumetric(
    const scene_data& scene, const instance_data& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

// check if a brdf is a delta
static bool is_delta(const material_point& material) {
  return (material.type == material_type::reflective &&
             material.roughness == 0) ||
         (material.type == material_type::refractive &&
             material.roughness == 0) ||
         (material.type == material_type::transparent &&
             material.roughness == 0) ||
         (material.type == material_type::volumetric);
}

static ray3f eval_camera(
    const cucamera_data& camera, const vec2f& image_uv, const vec2f& lens_uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  auto q    = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), camera.lens};
  // ray direction through the lens center
  auto dc = -normalize(q);
  // point on the lens
  auto e = vec3f{
      lens_uv.x * camera.aperture / 2, lens_uv.y * camera.aperture / 2, 0};
  // point on the focus plane
  auto p = dc * camera.focus / abs(dc.z);
  // correct ray direction to account for camera focusing
  auto d = normalize(p - e);
  // done
  return ray3f{
      transform_point(camera.frame, e), transform_direction(camera.frame, d)};
}

// Evaluate environment color.
static vec3f eval_environment(const scene_data& scene,
    const environment_data& environment, const vec3f& direction) {
  auto wl       = transform_direction_inverse(environment.frame, direction);
  auto texcoord = vec2f{
      atan2(wl.z, wl.x) / (2 * pif), acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
  if (texcoord.x < 0) texcoord.x += 1;
  return environment.emission *
         xyz(eval_texture(scene, environment.emission_tex, texcoord));
}

// Evaluate all environment color.
static vec3f eval_environment(const scene_data& scene, const vec3f& direction) {
  auto emission = vec3f{0, 0, 0};
  for (auto& environment : scene.environments) {
    emission += eval_environment(scene, environment, direction);
  }
  return emission;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-SCENE INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// intersection result
struct scene_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
  float _pad     = 0;
};

// closest hit
optix_shader void __closesthit__intersect_scene() {
  auto& intersection    = *getPRD<scene_intersection>();
  intersection.instance = optixGetInstanceIndex();
  intersection.element  = optixGetPrimitiveIndex();
  intersection.uv       = {
      optixGetTriangleBarycentrics().x, optixGetTriangleBarycentrics().y};
  intersection.distance = optixGetRayTmax();
  intersection.hit      = true;
}

// anyhit shader
optix_shader void __anyhit__intersect_scene() {}

// miss shader
optix_shader void __miss__intersect_scene() {
  auto& intersection    = *getPRD<scene_intersection>();
  intersection.instance = 0;
  intersection.element  = 0;
  intersection.uv       = {0, 0};
  intersection.distance = optixGetRayTmax();
  intersection.hit      = false;
}

// scene intersection via shaders
static scene_intersection intersect_scene(
    const trace_bvh& bvh, const cuscene_data& scene, const ray3f& ray) {
  auto     intersection = scene_intersection{};
  uint32_t u0, u1;
  packPointer(&intersection, u0, u1);
  optixTrace(bvh, {ray.o.x, ray.o.y, ray.o.z}, {ray.d.x, ray.d.y, ray.d.z},
      ray.tmin, ray.tmax, 0.0f, OptixVisibilityMask(255),
      OPTIX_RAY_FLAG_DISABLE_ANYHIT, 0, 0, 0, u0, u1);
  return intersection;
}

// instance intersection, for now manual
static scene_intersection intersect_instance(const trace_bvh& bvh,
    const cuscene_data& scene, int instance_id, const ray3f& ray) {
  auto& instance     = scene.instances[instance_id];
  auto& shape        = scene.shapes[instance.shape];
  auto  intersection = scene_intersection{};
  auto  tray         = ray3f{transform_point_inverse(instance.frame, ray.o),
      transform_vector_inverse(instance.frame, ray.d)};
  for (auto element = 0; element < shape.triangles.size(); element++) {
    auto& triangle = shape.triangles[element];
    auto  isec     = intersect_triangle(tray, shape.positions[triangle.x],
             shape.positions[triangle.y], shape.positions[triangle.z]);
    if (!isec.hit) continue;
    intersection.hit      = true;
    intersection.instance = instance_id;
    intersection.element  = element;
    intersection.uv       = isec.uv;
    intersection.distance = isec.distance;
    tray.tmax             = isec.distance;
  }
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience functions
[[maybe_unused]] static vec3f eval_position(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_element_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_element_normal(
      scene, scene.instances[intersection.instance], intersection.element);
}
[[maybe_unused]] static vec3f eval_shading_position(const scene_data& scene,
    const scene_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec3f eval_shading_normal(const scene_data& scene,
    const scene_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec2f eval_texcoord(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_texcoord(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static material_point eval_material(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_material(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static bool is_volumetric(
    const scene_data& scene, const scene_intersection& intersection) {
  return is_volumetric(scene, scene.instances[intersection.instance]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_emission(const material_point& material, const vec3f& normal,
    const vec3f& outgoing) {
  return dot(normal, outgoing) >= 0 ? material.emission : vec3f{0, 0, 0};
}

// Evaluates/sample the BRDF scaled by the cosine of the incoming direction.
static vec3f eval_bsdfcos(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (material.roughness == 0) return {0, 0, 0};

  if (material.type == material_type::matte) {
    return eval_matte(material.color, normal, outgoing, incoming);
  } else if (material.type == material_type::glossy) {
    return eval_glossy(material.color, material.ior, material.roughness, normal,
        outgoing, incoming);
  } else if (material.type == material_type::reflective) {
    return eval_reflective(
        material.color, material.roughness, normal, outgoing, incoming);
  } else if (material.type == material_type::transparent) {
    return eval_transparent(material.color, material.ior, material.roughness,
        normal, outgoing, incoming);
  } else if (material.type == material_type::refractive) {
    return eval_refractive(material.color, material.ior, material.roughness,
        normal, outgoing, incoming);
  } else if (material.type == material_type::subsurface) {
    return eval_refractive(material.color, material.ior, material.roughness,
        normal, outgoing, incoming);
  } else if (material.type == material_type::gltfpbr) {
    return eval_gltfpbr(material.color, material.ior, material.roughness,
        material.metallic, normal, outgoing, incoming);
  } else {
    return {0, 0, 0};
  }
}

static vec3f eval_delta(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, const vec3f& incoming) {
  if (material.roughness != 0) return {0, 0, 0};

  if (material.type == material_type::reflective) {
    return eval_reflective(material.color, normal, outgoing, incoming);
  } else if (material.type == material_type::transparent) {
    return eval_transparent(
        material.color, material.ior, normal, outgoing, incoming);
  } else if (material.type == material_type::refractive) {
    return eval_refractive(
        material.color, material.ior, normal, outgoing, incoming);
  } else if (material.type == material_type::volumetric) {
    return eval_passthrough(material.color, normal, outgoing, incoming);
  } else {
    return {0, 0, 0};
  }
}

// Picks a direction based on the BRDF
static vec3f sample_bsdfcos(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  if (material.roughness == 0) return {0, 0, 0};

  if (material.type == material_type::matte) {
    return sample_matte(material.color, normal, outgoing, rn);
  } else if (material.type == material_type::glossy) {
    return sample_glossy(material.color, material.ior, material.roughness,
        normal, outgoing, rnl, rn);
  } else if (material.type == material_type::reflective) {
    return sample_reflective(
        material.color, material.roughness, normal, outgoing, rn);
  } else if (material.type == material_type::transparent) {
    return sample_transparent(material.color, material.ior, material.roughness,
        normal, outgoing, rnl, rn);
  } else if (material.type == material_type::refractive) {
    return sample_refractive(material.color, material.ior, material.roughness,
        normal, outgoing, rnl, rn);
  } else if (material.type == material_type::subsurface) {
    return sample_refractive(material.color, material.ior, material.roughness,
        normal, outgoing, rnl, rn);
  } else if (material.type == material_type::gltfpbr) {
    return sample_gltfpbr(material.color, material.ior, material.roughness,
        material.metallic, normal, outgoing, rnl, rn);
  } else {
    return {0, 0, 0};
  }
}

static vec3f sample_delta(const material_point& material, const vec3f& normal,
    const vec3f& outgoing, float rnl) {
  if (material.roughness != 0) return {0, 0, 0};

  if (material.type == material_type::reflective) {
    return sample_reflective(material.color, normal, outgoing);
  } else if (material.type == material_type::transparent) {
    return sample_transparent(
        material.color, material.ior, normal, outgoing, rnl);
  } else if (material.type == material_type::refractive) {
    return sample_refractive(
        material.color, material.ior, normal, outgoing, rnl);
  } else if (material.type == material_type::volumetric) {
    return sample_passthrough(material.color, normal, outgoing);
  } else {
    return {0, 0, 0};
  }
}

// Compute the weight for sampling the BRDF
static float sample_bsdfcos_pdf(const material_point& material,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (material.roughness == 0) return 0;

  if (material.type == material_type::matte) {
    return sample_matte_pdf(material.color, normal, outgoing, incoming);
  } else if (material.type == material_type::glossy) {
    return sample_glossy_pdf(material.color, material.ior, material.roughness,
        normal, outgoing, incoming);
  } else if (material.type == material_type::reflective) {
    return sample_reflective_pdf(
        material.color, material.roughness, normal, outgoing, incoming);
  } else if (material.type == material_type::transparent) {
    return sample_tranparent_pdf(material.color, material.ior,
        material.roughness, normal, outgoing, incoming);
  } else if (material.type == material_type::refractive) {
    return sample_refractive_pdf(material.color, material.ior,
        material.roughness, normal, outgoing, incoming);
  } else if (material.type == material_type::subsurface) {
    return sample_refractive_pdf(material.color, material.ior,
        material.roughness, normal, outgoing, incoming);
  } else if (material.type == material_type::gltfpbr) {
    return sample_gltfpbr_pdf(material.color, material.ior, material.roughness,
        material.metallic, normal, outgoing, incoming);
  } else {
    return 0;
  }
}

static float sample_delta_pdf(const material_point& material,
    const vec3f& normal, const vec3f& outgoing, const vec3f& incoming) {
  if (material.roughness != 0) return 0;

  if (material.type == material_type::reflective) {
    return sample_reflective_pdf(material.color, normal, outgoing, incoming);
  } else if (material.type == material_type::transparent) {
    return sample_tranparent_pdf(
        material.color, material.ior, normal, outgoing, incoming);
  } else if (material.type == material_type::refractive) {
    return sample_refractive_pdf(
        material.color, material.ior, normal, outgoing, incoming);
  } else if (material.type == material_type::volumetric) {
    return sample_passthrough_pdf(material.color, normal, outgoing, incoming);
  } else {
    return 0;
  }
}

static vec3f eval_scattering(const material_point& material,
    const vec3f& outgoing, const vec3f& incoming) {
  if (material.density == vec3f{0, 0, 0}) return {0, 0, 0};
  return material.scattering * material.density *
         eval_phasefunction(material.scanisotropy, outgoing, incoming);
}

static vec3f sample_scattering(const material_point& material,
    const vec3f& outgoing, float rnl, const vec2f& rn) {
  if (material.density == vec3f{0, 0, 0}) return {0, 0, 0};
  return sample_phasefunction(material.scanisotropy, outgoing, rn);
}

static float sample_scattering_pdf(const material_point& material,
    const vec3f& outgoing, const vec3f& incoming) {
  if (material.density == vec3f{0, 0, 0}) return 0;
  return sample_phasefunction_pdf(material.scanisotropy, outgoing, incoming);
}

// Sample camera
static ray3f sample_camera(const camera_data& camera, const vec2i& ij,
    const vec2i& image_size, const vec2f& puv, const vec2f& luv, bool tent) {
  if (!tent) {
    auto uv = vec2f{
        (ij.x + puv.x) / image_size.x, (ij.y + puv.y) / image_size.y};
    return eval_camera(camera, uv, sample_disk(luv));
  } else {
    const auto width  = 2.0f;
    const auto offset = 0.5f;
    auto       fuv =
        width *
            vec2f{
                puv.x < 0.5f ? sqrt(2 * puv.x) - 1 : 1 - sqrt(2 - 2 * puv.x),
                puv.y < 0.5f ? sqrt(2 * puv.y) - 1 : 1 - sqrt(2 - 2 * puv.y),
            } +
        offset;
    auto uv = vec2f{
        (ij.x + fuv.x) / image_size.x, (ij.y + fuv.y) / image_size.y};
    return eval_camera(camera, uv, sample_disk(luv));
  }
}

// Sample lights wrt solid angle
static vec3f sample_lights(const scene_data& scene, const trace_lights& lights,
    const vec3f& position, float rl, float rel, const vec2f& ruv) {
  auto  light_id = sample_uniform((int)lights.lights.size(), rl);
  auto& light    = lights.lights[light_id];
  if (light.instance != invalidid) {
    auto& instance  = scene.instances[light.instance];
    auto& shape     = scene.shapes[instance.shape];
    auto  element   = sample_discrete(light.elements_cdf, rel);
    auto  uv        = (!shape.triangles.empty()) ? sample_triangle(ruv) : ruv;
    auto  lposition = eval_position(scene, instance, element, uv);
    return normalize(lposition - position);
  } else if (light.environment != invalidid) {
    auto& environment = scene.environments[light.environment];
    if (environment.emission_tex != invalidid) {
      auto& emission_tex = scene.textures[environment.emission_tex];
      auto  idx          = sample_discrete(light.elements_cdf, rel);
      auto  uv = vec2f{((idx % emission_tex.width) + 0.5f) / emission_tex.width,
          ((idx / emission_tex.width) + 0.5f) / emission_tex.height};
      return transform_direction(environment.frame,
          {cos(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif),
              sin(uv.x * 2 * pif) * sin(uv.y * pif)});
    } else {
      return sample_sphere(ruv);
    }
  } else {
    return {0, 0, 0};
  }
}

// Sample lights pdf
static float sample_lights_pdf(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const vec3f& position, const vec3f& direction) {
  auto pdf = 0.0f;
  for (auto& light : lights.lights) {
    if (light.instance != invalidid) {
      auto& instance = scene.instances[light.instance];
      // check all intersection
      auto lpdf          = 0.0f;
      auto next_position = position;
      for (auto bounce = 0; bounce < 100; bounce++) {
        auto intersection = intersect_instance(
            bvh, scene, light.instance, {next_position, direction});
        if (!intersection.hit) break;
        // accumulate pdf
        auto lposition = eval_position(
            scene, instance, intersection.element, intersection.uv);
        auto lnormal = eval_element_normal(
            scene, instance, intersection.element);
        // prob triangle * area triangle = area triangle mesh
        auto area = light.elements_cdf.back();
        lpdf += distance_squared(lposition, position) /
                (abs(dot(lnormal, direction)) * area);
        // continue
        next_position = lposition + direction * 1e-3f;
      }
      pdf += lpdf;
    } else if (light.environment != invalidid) {
      auto& environment = scene.environments[light.environment];
      if (environment.emission_tex != invalidid) {
        auto& emission_tex = scene.textures[environment.emission_tex];
        auto  wl = transform_direction_inverse(environment.frame, direction);
        auto  texcoord = vec2f{atan2(wl.z, wl.x) / (2 * pif),
            acos(clamp(wl.y, -1.0f, 1.0f)) / pif};
        if (texcoord.x < 0) texcoord.x += 1;
        auto i = clamp(
            (int)(texcoord.x * emission_tex.width), 0, emission_tex.width - 1);
        auto j    = clamp((int)(texcoord.y * emission_tex.height), 0,
               emission_tex.height - 1);
        auto prob = sample_discrete_pdf(
                        light.elements_cdf, j * emission_tex.width + i) /
                    light.elements_cdf.back();
        auto angle = (2 * pif / emission_tex.width) *
                     (pif / emission_tex.height) *
                     sin(pif * (j + 0.5f) / emission_tex.height);
        pdf += prob / angle;
      } else {
        pdf += 1 / (4 * pif);
      }
    }
  }
  pdf *= sample_uniform_pdf((int)lights.lights.size());
  return pdf;
}

struct trace_result {
  vec3f radiance = {0, 0, 0};
  bool  hit      = false;
  vec3f albedo   = {0, 0, 0};
  vec3f normal   = {0, 0, 0};
};

// Recursive path tracing.
static trace_result trace_path(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance      = vec3f{0, 0, 0};
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = svector<material_point>{};
  auto max_roughness = 0.0f;
  auto hit           = false;
  auto hit_albedo    = vec3f{0, 0, 0};
  auto hit_normal    = vec3f{0, 0, 0};
  auto opbounce      = 0;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // handle transmission if inside a volume
    auto in_volume = false;
    if (!volume_stack.empty()) {
      auto& vsdf     = volume_stack.back();
      auto  distance = sample_transmittance(
           vsdf.density, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(vsdf.density, distance) /
                sample_transmittance_pdf(
                    vsdf.density, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto outgoing = -ray.d;
      auto position = eval_shading_position(scene, intersection, outgoing);
      auto normal   = eval_shading_normal(scene, intersection, outgoing);
      auto material = eval_material(scene, intersection);

      // correct roughness
      if (params.nocaustics) {
        max_roughness      = max(material.roughness, max_roughness);
        material.roughness = max_roughness;
      }

      // handle opacity
      if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
        if (opbounce++ > 128) break;
        ray = {position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }

      // set hit variables
      if (bounce == 0) {
        hit        = true;
        hit_albedo = material.color;
        hit_normal = normal;
      }

      // accumulate emission
      radiance += weight * eval_emission(material, normal, outgoing);

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (!is_delta(material)) {
        if (rand1f(rng) < 0.5f) {
          incoming = sample_bsdfcos(
              material, normal, outgoing, rand1f(rng), rand2f(rng));
        } else {
          incoming = sample_lights(
              scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
        }
        if (incoming == vec3f{0, 0, 0}) break;
        weight *=
            eval_bsdfcos(material, normal, outgoing, incoming) /
            (0.5f * sample_bsdfcos_pdf(material, normal, outgoing, incoming) +
                0.5f *
                    sample_lights_pdf(scene, bvh, lights, position, incoming));
      } else {
        incoming = sample_delta(material, normal, outgoing, rand1f(rng));
        weight *= eval_delta(material, normal, outgoing, incoming) /
                  sample_delta_pdf(material, normal, outgoing, incoming);
      }

      // update volume stack
      if (is_volumetric(scene, intersection) &&
          dot(normal, outgoing) * dot(normal, incoming) < 0) {
        if (volume_stack.empty()) {
          auto material = eval_material(scene, intersection);
          volume_stack.push_back(material);
        } else {
          volume_stack.pop_back();
        }
      }

      // setup next iteration
      ray = {position, incoming};
    } else {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto  position = ray.o + ray.d * intersection.distance;
      auto& vsdf     = volume_stack.back();

      // accumulate emission
      // radiance += weight * eval_volemission(emission, outgoing);

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (rand1f(rng) < 0.5f) {
        incoming = sample_scattering(vsdf, outgoing, rand1f(rng), rand2f(rng));
      } else {
        incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      if (incoming == vec3f{0, 0, 0}) break;
      weight *=
          eval_scattering(vsdf, outgoing, incoming) /
          (0.5f * sample_scattering_pdf(vsdf, outgoing, incoming) +
              0.5f * sample_lights_pdf(scene, bvh, lights, position, incoming));

      // setup next iteration
      ray = {position, incoming};
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Recursive path tracing.
static trace_result trace_pathdirect(const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights, const ray3f& ray_,
    rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance      = vec3f{0, 0, 0};
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = svector<material_point>{};
  auto max_roughness = 0.0f;
  auto hit           = false;
  auto hit_albedo    = vec3f{0, 0, 0};
  auto hit_normal    = vec3f{0, 0, 0};
  auto next_emission = true;
  auto opbounce      = 0;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if ((bounce > 0 || !params.envhidden) && next_emission)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // handle transmission if inside a volume
    auto in_volume = false;
    if (!volume_stack.empty()) {
      auto& vsdf     = volume_stack.back();
      auto  distance = sample_transmittance(
           vsdf.density, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(vsdf.density, distance) /
                sample_transmittance_pdf(
                    vsdf.density, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto outgoing = -ray.d;
      auto position = eval_shading_position(scene, intersection, outgoing);
      auto normal   = eval_shading_normal(scene, intersection, outgoing);
      auto material = eval_material(scene, intersection);

      // correct roughness
      if (params.nocaustics) {
        max_roughness      = max(material.roughness, max_roughness);
        material.roughness = max_roughness;
      }

      // handle opacity
      if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
        if (opbounce++ > 128) break;
        ray = {position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }

      // set hit variables
      if (bounce == 0) {
        hit        = true;
        hit_albedo = material.color;
        hit_normal = normal;
      }

      // accumulate emission
      if (next_emission)
        radiance += weight * eval_emission(material, normal, outgoing);

      // direct
      if (!is_delta(material)) {
        auto incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
        auto pdf = sample_lights_pdf(scene, bvh, lights, position, incoming);
        auto bsdfcos = eval_bsdfcos(material, normal, outgoing, incoming);
        if (bsdfcos != vec3f{0, 0, 0} && pdf > 0) {
          auto intersection = intersect_scene(bvh, scene, {position, incoming});
          auto emission =
              !intersection.hit
                  ? eval_environment(scene, incoming)
                  : eval_emission(eval_material(scene,
                                      scene.instances[intersection.instance],
                                      intersection.element, intersection.uv),
                        eval_shading_normal(scene,
                            scene.instances[intersection.instance],
                            intersection.element, intersection.uv, -incoming),
                        -incoming);
          radiance += weight * bsdfcos * emission / pdf;
        }
        next_emission = false;
      } else {
        next_emission = true;
      }

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (!is_delta(material)) {
        if (rand1f(rng) < 0.5f) {
          incoming = sample_bsdfcos(
              material, normal, outgoing, rand1f(rng), rand2f(rng));
        } else {
          incoming = sample_lights(
              scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
        }
        if (incoming == vec3f{0, 0, 0}) break;
        weight *=
            eval_bsdfcos(material, normal, outgoing, incoming) /
            (0.5f * sample_bsdfcos_pdf(material, normal, outgoing, incoming) +
                0.5f *
                    sample_lights_pdf(scene, bvh, lights, position, incoming));
      } else {
        incoming = sample_delta(material, normal, outgoing, rand1f(rng));
        if (incoming == vec3f{0, 0, 0}) break;
        weight *= eval_delta(material, normal, outgoing, incoming) /
                  sample_delta_pdf(material, normal, outgoing, incoming);
      }

      // update volume stack
      if (is_volumetric(scene, intersection) &&
          dot(normal, outgoing) * dot(normal, incoming) < 0) {
        if (volume_stack.empty()) {
          auto material = eval_material(scene, intersection);
          volume_stack.push_back(material);
        } else {
          volume_stack.pop_back();
        }
      }

      // setup next iteration
      ray = {position, incoming};
    } else {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto  position = ray.o + ray.d * intersection.distance;
      auto& vsdf     = volume_stack.back();

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (rand1f(rng) < 0.5f) {
        incoming = sample_scattering(vsdf, outgoing, rand1f(rng), rand2f(rng));
      } else {
        incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      if (incoming == vec3f{0, 0, 0}) break;
      weight *=
          eval_scattering(vsdf, outgoing, incoming) /
          (0.5f * sample_scattering_pdf(vsdf, outgoing, incoming) +
              0.5f * sample_lights_pdf(scene, bvh, lights, position, incoming));

      // setup next iteration
      ray = {position, incoming};
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Recursive path tracing with MIS.
static trace_result trace_pathmis(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance      = vec3f{0, 0, 0};
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto volume_stack  = svector<material_point>{};
  auto max_roughness = 0.0f;
  auto hit           = false;
  auto hit_albedo    = vec3f{0, 0, 0};
  auto hit_normal    = vec3f{0, 0, 0};
  auto opbounce      = 0;

  // MIS helpers
  auto mis_heuristic = [](float this_pdf, float other_pdf) {
    return (this_pdf * this_pdf) /
           (this_pdf * this_pdf + other_pdf * other_pdf);
  };
  auto next_emission     = true;
  auto next_intersection = scene_intersection{};

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = next_emission ? intersect_scene(bvh, scene, ray)
                                      : next_intersection;
    if (!intersection.hit) {
      if ((bounce > 0 || !params.envhidden) && next_emission)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // handle transmission if inside a volume
    auto in_volume = false;
    if (!volume_stack.empty()) {
      auto& vsdf     = volume_stack.back();
      auto  distance = sample_transmittance(
           vsdf.density, intersection.distance, rand1f(rng), rand1f(rng));
      weight *= eval_transmittance(vsdf.density, distance) /
                sample_transmittance_pdf(
                    vsdf.density, distance, intersection.distance);
      in_volume             = distance < intersection.distance;
      intersection.distance = distance;
    }

    // switch between surface and volume
    if (!in_volume) {
      // prepare shading point
      auto outgoing = -ray.d;
      auto position = eval_shading_position(scene, intersection, outgoing);
      auto normal   = eval_shading_normal(scene, intersection, outgoing);
      auto material = eval_material(scene, intersection);

      // correct roughness
      if (params.nocaustics) {
        max_roughness      = max(material.roughness, max_roughness);
        material.roughness = max_roughness;
      }

      // handle opacity
      if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
        if (opbounce++ > 128) break;
        ray = {position + ray.d * 1e-2f, ray.d};
        bounce -= 1;
        continue;
      }

      // set hit variables
      if (bounce == 0) {
        hit        = true;
        hit_albedo = material.color;
        hit_normal = normal;
      }

      // accumulate emission
      if (next_emission) {
        radiance += weight * eval_emission(material, normal, outgoing);
      }

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (!is_delta(material)) {
        // direct with MIS --- light
        for (auto sample_light : {true, false}) {
          incoming = sample_light ? sample_lights(scene, lights, position,
                                        rand1f(rng), rand1f(rng), rand2f(rng))
                                  : sample_bsdfcos(material, normal, outgoing,
                                        rand1f(rng), rand2f(rng));
          if (incoming == vec3f{0, 0, 0}) break;
          auto bsdfcos   = eval_bsdfcos(material, normal, outgoing, incoming);
          auto light_pdf = sample_lights_pdf(
              scene, bvh, lights, position, incoming);
          auto bsdf_pdf = sample_bsdfcos_pdf(
              material, normal, outgoing, incoming);
          auto mis_weight = sample_light
                                ? mis_heuristic(light_pdf, bsdf_pdf) / light_pdf
                                : mis_heuristic(bsdf_pdf, light_pdf) / bsdf_pdf;
          if (bsdfcos != vec3f{0, 0, 0} && mis_weight != 0) {
            auto intersection = intersect_scene(
                bvh, scene, {position, incoming});
            if (!sample_light) next_intersection = intersection;
            auto emission = vec3f{0, 0, 0};
            if (!intersection.hit) {
              emission = eval_environment(scene, incoming);
            } else {
              auto material = eval_material(scene,
                  scene.instances[intersection.instance], intersection.element,
                  intersection.uv);
              emission      = eval_emission(material,
                       eval_shading_normal(scene,
                           scene.instances[intersection.instance],
                           intersection.element, intersection.uv, -incoming),
                       -incoming);
            }
            radiance += weight * bsdfcos * emission * mis_weight;
          }
        }

        // indirect
        weight *= eval_bsdfcos(material, normal, outgoing, incoming) /
                  sample_bsdfcos_pdf(material, normal, outgoing, incoming);
        next_emission = false;
      } else {
        incoming = sample_delta(material, normal, outgoing, rand1f(rng));
        weight *= eval_delta(material, normal, outgoing, incoming) /
                  sample_delta_pdf(material, normal, outgoing, incoming);
        next_emission = true;
      }

      // update volume stack
      if (is_volumetric(scene, intersection) &&
          dot(normal, outgoing) * dot(normal, incoming) < 0) {
        if (volume_stack.empty()) {
          auto material = eval_material(scene, intersection);
          volume_stack.push_back(material);
        } else {
          volume_stack.pop_back();
        }
      }

      // setup next iteration
      ray = {position, incoming};
    } else {
      // prepare shading point
      auto  outgoing = -ray.d;
      auto  position = ray.o + ray.d * intersection.distance;
      auto& vsdf     = volume_stack.back();

      // next direction
      auto incoming = vec3f{0, 0, 0};
      if (rand1f(rng) < 0.5f) {
        incoming = sample_scattering(vsdf, outgoing, rand1f(rng), rand2f(rng));
        next_emission = true;
      } else {
        incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
        next_emission = true;
      }
      weight *=
          eval_scattering(vsdf, outgoing, incoming) /
          (0.5f * sample_scattering_pdf(vsdf, outgoing, incoming) +
              0.5f * sample_lights_pdf(scene, bvh, lights, position, incoming));

      // setup next iteration
      ray = {position, incoming};
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Recursive path tracing.
static trace_result trace_pathtest(const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights, const ray3f& ray_,
    rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance      = vec3f{0, 0, 0};
  auto weight        = vec3f{1, 1, 1};
  auto ray           = ray_;
  auto max_roughness = 0.0f;
  auto hit           = false;
  auto hit_albedo    = vec3f{0, 0, 0};
  auto hit_normal    = vec3f{0, 0, 0};
  auto opbounce      = 0;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);
    material.type = material_type::matte;

    // set hit variables
    if (bounce == 0) {
      hit        = true;
      hit_albedo = material.color;
      hit_normal = normal;
    }

    // accumulate emission
    radiance += weight * eval_emission(material, normal, outgoing);

    // next direction
    auto incoming = vec3f{0, 0, 0};
    if (!is_delta(material)) {
      if (rand1f(rng) < 0.5f) {
        incoming = sample_bsdfcos(
            material, normal, outgoing, rand1f(rng), rand2f(rng));
      } else {
        incoming = sample_lights(
            scene, lights, position, rand1f(rng), rand1f(rng), rand2f(rng));
      }
      if (incoming == vec3f{0, 0, 0}) break;
      weight *=
          eval_bsdfcos(material, normal, outgoing, incoming) /
          (0.5f * sample_bsdfcos_pdf(material, normal, outgoing, incoming) +
              0.5f * sample_lights_pdf(scene, bvh, lights, position, incoming));
    } else {
      incoming = sample_delta(material, normal, outgoing, rand1f(rng));
      weight *= eval_delta(material, normal, outgoing, incoming) /
                sample_delta_pdf(material, normal, outgoing, incoming);
    }

    // setup next iteration
    ray = {position, incoming};

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Recursive path tracing.
static trace_result trace_naive(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance   = vec3f{0, 0, 0};
  auto weight     = vec3f{1, 1, 1};
  auto ray        = ray_;
  auto hit        = false;
  auto hit_albedo = vec3f{0, 0, 0};
  auto hit_normal = vec3f{0, 0, 0};
  auto opbounce   = 0;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      if (opbounce++ > 128) break;
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    // set hit variables
    if (bounce == 0) {
      hit        = true;
      hit_albedo = material.color;
      hit_normal = normal;
    }

    // accumulate emission
    radiance += weight * eval_emission(material, normal, outgoing);

    // next direction
    auto incoming = vec3f{0, 0, 0};
    if (material.roughness != 0) {
      incoming = sample_bsdfcos(
          material, normal, outgoing, rand1f(rng), rand2f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_bsdfcos(material, normal, outgoing, incoming) /
                sample_bsdfcos_pdf(material, normal, outgoing, incoming);
    } else {
      incoming = sample_delta(material, normal, outgoing, rand1f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_delta(material, normal, outgoing, incoming) /
                sample_delta_pdf(material, normal, outgoing, incoming);
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Eyelight for quick previewing.
static trace_result trace_eyelight(const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights, const ray3f& ray_,
    rng_state& rng, const trace_params& params) {
  // initialize
  auto radiance   = vec3f{0, 0, 0};
  auto weight     = vec3f{1, 1, 1};
  auto ray        = ray_;
  auto hit        = false;
  auto hit_albedo = vec3f{0, 0, 0};
  auto hit_normal = vec3f{0, 0, 0};
  auto opbounce   = 0;

  // trace  path
  for (auto bounce = 0; bounce < max(params.bounces, 4); bounce++) {
    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto outgoing = -ray.d;
    auto position = eval_shading_position(scene, intersection, outgoing);
    auto normal   = eval_shading_normal(scene, intersection, outgoing);
    auto material = eval_material(scene, intersection);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      if (opbounce++ > 128) break;
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    // set hit variables
    if (bounce == 0) {
      hit        = true;
      hit_albedo = material.color;
      hit_normal = normal;
    }

    // accumulate emission
    auto incoming = outgoing;
    radiance += weight * eval_emission(material, normal, outgoing);

    // brdf * light
    radiance += weight * pif *
                eval_bsdfcos(material, normal, outgoing, incoming);

    // continue path
    if (!is_delta(material)) break;
    incoming = sample_delta(material, normal, outgoing, rand1f(rng));
    if (incoming == vec3f{0, 0, 0}) break;
    weight *= eval_delta(material, normal, outgoing, incoming) /
              sample_delta_pdf(material, normal, outgoing, incoming);
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // setup next iteration
    ray = {position, incoming};
  }

  return {radiance, hit, hit_albedo, hit_normal};
}

// Furnace test.
static trace_result trace_furnace(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray_, rng_state& rng,
    const trace_params& params) {
  // initialize
  auto radiance   = vec3f{0, 0, 0};
  auto weight     = vec3f{1, 1, 1};
  auto ray        = ray_;
  auto hit        = false;
  auto hit_albedo = vec3f{0, 0, 0};
  auto hit_normal = vec3f{0, 0, 0};
  auto opbounce   = 0;
  auto in_volume  = false;

  // trace  path
  for (auto bounce = 0; bounce < params.bounces; bounce++) {
    // exit loop
    if (bounce > 0 && !in_volume) {
      radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // intersect next point
    auto intersection = intersect_scene(bvh, scene, ray);
    if (!intersection.hit) {
      if (bounce > 0 || !params.envhidden)
        radiance += weight * eval_environment(scene, ray.d);
      break;
    }

    // prepare shading point
    auto  outgoing = -ray.d;
    auto& instance = scene.instances[intersection.instance];
    auto  element  = intersection.element;
    auto  uv       = intersection.uv;
    auto  position = eval_position(scene, instance, element, uv);
    auto  normal = eval_shading_normal(scene, instance, element, uv, outgoing);
    auto  material = eval_material(scene, instance, element, uv);

    // handle opacity
    if (material.opacity < 1 && rand1f(rng) >= material.opacity) {
      if (opbounce++ > 128) break;
      ray = {position + ray.d * 1e-2f, ray.d};
      bounce -= 1;
      continue;
    }

    // set hit variables
    if (bounce == 0) {
      hit        = true;
      hit_albedo = material.color;
      hit_normal = normal;
    }

    // accumulate emission
    radiance += weight * eval_emission(material, normal, outgoing);

    // next direction
    auto incoming = vec3f{0, 0, 0};
    if (material.roughness != 0) {
      incoming = sample_bsdfcos(
          material, normal, outgoing, rand1f(rng), rand2f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_bsdfcos(material, normal, outgoing, incoming) /
                sample_bsdfcos_pdf(material, normal, outgoing, incoming);
    } else {
      incoming = sample_delta(material, normal, outgoing, rand1f(rng));
      if (incoming == vec3f{0, 0, 0}) break;
      weight *= eval_delta(material, normal, outgoing, incoming) /
                sample_delta_pdf(material, normal, outgoing, incoming);
    }

    // check weight
    if (weight == vec3f{0, 0, 0} || !isfinite(weight)) break;

    // russian roulette
    if (bounce > 3) {
      auto rr_prob = min((float)0.99, max(weight));
      if (rand1f(rng) >= rr_prob) break;
      weight *= 1 / rr_prob;
    }

    // update volume stack
    if (dot(normal, outgoing) * dot(normal, incoming) < 0)
      in_volume = !in_volume;

    // setup next iteration
    ray = {position, incoming};
  }

  // done
  return {radiance, hit, hit_albedo, hit_normal};
}

// False color rendering
static trace_result trace_falsecolor(const scene_data& scene,
    const trace_bvh& bvh, const trace_lights& lights, const ray3f& ray,
    rng_state& rng, const trace_params& params) {
  // intersect next point
  auto intersection = intersect_scene(bvh, scene, ray);
  if (!intersection.hit) return {};

  // prepare shading point
  auto outgoing = -ray.d;
  auto position = eval_shading_position(scene, intersection, outgoing);
  auto normal   = eval_shading_normal(scene, intersection, outgoing);
  auto gnormal  = eval_element_normal(scene, intersection);
  auto texcoord = eval_texcoord(scene, intersection);
  auto material = eval_material(scene, intersection);
  auto delta    = is_delta(material) ? 1.0f : 0.0f;

  // hash color
  auto hashed_color = [](int id) {
    auto rng = make_rng(trace_default_seed, id * 2 + 1);
    return pow(0.5f + 0.5f * rand3f(rng), 2.2f);
  };

  // compute result
  auto result = vec3f{0, 0, 0};
  switch (params.falsecolor) {
    case trace_falsecolor_type::position:
      result = position * 0.5f + 0.5f;
      break;
    case trace_falsecolor_type::normal: result = normal * 0.5f + 0.5f; break;
    case trace_falsecolor_type::frontfacing:
      result = dot(normal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
      break;
    case trace_falsecolor_type::gnormal: result = gnormal * 0.5f + 0.5f; break;
    case trace_falsecolor_type::gfrontfacing:
      result = dot(gnormal, -ray.d) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
      break;
    case trace_falsecolor_type::mtype:
      result = hashed_color((int)material.type);
      break;
    case trace_falsecolor_type::texcoord:
      result = {fmod(texcoord.x, 1.0f), fmod(texcoord.y, 1.0f), 0};
      break;
    case trace_falsecolor_type::color: result = material.color; break;
    case trace_falsecolor_type::emission: result = material.emission; break;
    case trace_falsecolor_type::roughness:
      result = {material.roughness, material.roughness, material.roughness};
      break;
    case trace_falsecolor_type::opacity:
      result = {material.opacity, material.opacity, material.opacity};
      break;
    case trace_falsecolor_type::metallic:
      result = {material.metallic, material.metallic, material.metallic};
      break;
    case trace_falsecolor_type::delta: result = {delta, delta, delta}; break;
    case trace_falsecolor_type::element:
      result = hashed_color(intersection.element);
      break;
    case trace_falsecolor_type::instance:
      result = hashed_color(intersection.instance);
      break;
    case trace_falsecolor_type::shape:
      result = hashed_color(scene.instances[intersection.instance].shape);
      break;
    case trace_falsecolor_type::material:
      result = hashed_color(scene.instances[intersection.instance].material);
      break;
    case trace_falsecolor_type::highlight: {
      if (material.emission == vec3f{0, 0, 0})
        material.emission = {0.2f, 0.2f, 0.2f};
      result = material.emission * abs(dot(-ray.d, normal));
      break;
    } break;
    default: result = {0, 0, 0};
  }

  // done
  return {srgb_to_rgb(result), true, material.color, normal};
}

static trace_result trace_sampler(const scene_data& scene, const trace_bvh& bvh,
    const trace_lights& lights, const ray3f& ray, rng_state& rng,
    const trace_params& params) {
  switch (params.sampler) {
    case trace_sampler_type::path:
      return trace_path(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::pathdirect:
      return trace_pathdirect(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::pathmis:
      return trace_pathmis(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::pathtest:
      return trace_pathtest(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::naive:
      return trace_naive(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::eyelight:
      return trace_eyelight(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::furnace:
      return trace_furnace(scene, bvh, lights, ray, rng, params);
    case trace_sampler_type::falsecolor:
      return trace_falsecolor(scene, bvh, lights, ray, rng, params);
    default: {
      return {};
    }
  }
}

static void trace_sample(cutrace_state& state, const cuscene_data& scene,
    const cutrace_bvh& bvh, const cutrace_lights& lights, int i, int j,
    int sample, const trace_params& params) {
  auto& camera = scene.cameras[params.camera];
  // auto  sampler = get_trace_sampler_func(params);
  auto idx    = state.width * j + i;
  auto ray    = sample_camera(camera, {i, j}, {state.width, state.height},
         rand2f(state.rngs[idx]), rand2f(state.rngs[idx]), params.tentfilter);
  auto result = trace_sampler(scene, bvh, lights, ray, state.rngs[idx], params);
  // auto [radiance, hit, albedo, normal] = sampler(
  //    scene, bvh, lights, ray, state.rngs[idx], params);
  auto radiance = result.radiance;
  auto hit      = result.hit;
  auto albedo   = result.albedo;
  auto normal   = result.normal;
  if (!isfinite(radiance)) radiance = {0, 0, 0};
  if (max(radiance) > params.clamp)
    radiance = radiance * (params.clamp / max(radiance));
  auto weight = 1.0f / (sample + 1);
  if (hit) {
    state.image[idx] = lerp(
        state.image[idx], {radiance.x, radiance.y, radiance.z, 1}, weight);
    state.albedo[idx] = lerp(state.albedo[idx], albedo, weight);
    state.normal[idx] = lerp(state.normal[idx], normal, weight);
    state.hits[idx] += 1;
  } else if (!params.envhidden && !scene.environments.empty()) {
    state.image[idx] = lerp(
        state.image[idx], {radiance.x, radiance.y, radiance.z, 1}, weight);
    state.albedo[idx] = lerp(state.albedo[idx], {1, 1, 1}, weight);
    state.normal[idx] = lerp(state.normal[idx], -ray.d, weight);
    state.hits[idx] += 1;
  } else {
    state.image[idx]  = lerp(state.image[idx], {0, 0, 0, 0}, weight);
    state.albedo[idx] = lerp(state.albedo[idx], {0, 0, 0}, weight);
    state.normal[idx] = lerp(state.normal[idx], -ray.d, weight);
  }
}

// raygen shader
optix_shader void __raygen__trace_pixel() {
  // pixel index
  auto ij  = optixGetLaunchIndex();
  auto idx = ij.y * globals.state.width + ij.x;

  // initialize state on first sample
  if (globals.state.samples == 0) {
    globals.state.image[idx] = {0, 0, 0, 0};
    globals.state.rngs[idx]  = make_rng(98273987, idx * 2 + 1);
  }

  // run shading
  auto ssample  = globals.state.samples;
  auto nsamples = globals.params.batch;
  for (auto sample = ssample; sample < ssample + nsamples; sample++) {
    trace_sample(globals.state, globals.scene, globals.bvh, globals.lights,
        optixGetLaunchIndex().x, optixGetLaunchIndex().y, sample,
        globals.params);
  }
}

}  // namespace yocto