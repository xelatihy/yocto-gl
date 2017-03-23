///
/// YOCTO_TRACE: path tracer implementation for support for textured mesh
/// lights, GGX/Phong materials, environment mapping. The interface supports
/// progressive parallel execution with any splitting strategy, by
/// generating per-sample random number sequences or fully deterministic
/// hash-based sampling.
///
///
/// USAGE:
///
/// 1. define your scene setting up only the public data
/// - init the scene with make_scene()
/// - define cameras with set_camera()
/// - defile shapes with set set_triagle_shape(), set_line_shape(),
/// set_point_shape()
/// - define materials with set_material()
/// - define textures with set_texture()
/// - define environments with set_environment()
/// - set intersection routines with set_intersection_callbacks()
///     - can use yocto_bvh
/// 2. prepare for rendering with init_lights()
/// 3. define rendering params in render_params
/// 4. render blocks of samples with trace_block() or the whole image with
/// trace_image()
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// Shapes are indexed meshes and are described by array of vertex indices for
/// points, lines and triangles, and arrays of vertex data. Only one primitive
/// type can be non-empty for each shape.
///
/// Materials are represented as sums of an emission term, a diffuse term and
/// a specular microfacet term (GGX or Phong). Only opaque for now. We pick
/// a proper material type for each shape element type (points, lines,
/// triangles).
///
/// Lights are defined as any shape with a material emission term. Additionally
/// one can also add an environment map. But even if you can, you might want to
/// add a large triangle mesh with inward normals instead. The latter is more
/// general (you can even more an arbitrary shape sun). For now only the first
/// env is used.
///
/// We generate our own random numbers guarantying that there is one random
/// sequence per path. This means you can rul the path tracer in any order
/// serially or in parallel.
///
/// For now, we support a straightforward path tracer with explicit direct
/// illumination using MIS.
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
///
/// The .cpp file depends on yocto_math.h.
///
///
/// HISTORY:
/// - v 1.14: normal mapping
/// - v 1.13: simpler Fresnel handling
/// - v 1.12: significantly better path tracing
/// - v 1.11: add progressive sampling to rendering params
/// - v 0.10: switch to .h/.cpp pair
/// - v 0.9: doxygen comments
/// - v 0.8: opaque API (allows for changing internals without altering API)
/// - v 0.7: internally use pointers for performance transparency
/// - v 0.6: minor API change for blocks
/// - v 0.5: [major API change] move to modern C++ interface
/// - v 0.4: C++ API
/// - v 0.3: removal of C interface
/// - v 0.2: use of STL containers
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace ytrace {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YTRACE_H_
#define _YTRACE_H_

#ifdef YTRACE_INLINE
#define YTRACE_API inline
#else
#define YTRACE_API
#endif

#include <array>
#include <cstdarg>
#include <functional>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ytrace {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using byte = unsigned char;
using byte2 = std::array<byte, 2>;
using byte3 = std::array<byte, 3>;
using byte4 = std::array<byte, 4>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;

///
/// Trace scene.
///
struct scene;

///
/// Initialize the scene with the proper number of objects.
///
YTRACE_API scene* make_scene(int ncameras, int nshapes, int nmaterials,
    int ntextures, int nenvironments);

///
/// Free scene.
///
YTRACE_API void free_scene(scene* scn);

///
/// Sets a camera in the scene.
///
/// Parameters:
/// - scn: scene
/// - cid: camera id
/// - frame: local-to-world frame (x, y, z, o in column major order)
/// - yfov: field of view
/// - aspect: aspect ratio
/// - aperture: lens aperture
/// - focus: focus plane distance (cannot be zero)
///
YTRACE_API void set_camera(scene* scn, int cid, const float3x4& frame,
    float yfov, float aspect, float aperture = 0, float focus = 1);

///
/// Sets a texture in the scene.
///
/// Parameters:
/// - scn: scene
/// - tid: texture id
/// - width: width
/// - height: height
/// - ncomp: number of components (1-4)
/// - hdr: hdr pixels
/// - ldr: ldr pixels (sRGB)
///
YTRACE_API void set_texture(
    scene* scn, int tid, int width, int height, int ncomp, const float* hdr);

///
/// Sets a texture in the scene.
///
/// Parameters:
/// - scn: scene
/// - tid: texture id
/// - width: width
/// - height: height
/// - ncomp: number of components (1-4)
/// - hdr: hdr pixels
/// - ldr: ldr pixels (sRGB)
///
YTRACE_API void set_texture(
    scene* scn, int tid, int width, int height, int ncomp, const byte* ldr);

///
/// Sets a material in the scene.
///
/// Parameters:
/// - scn: scene
/// - mid: material id
/// - ke: emission, term
/// - kd: diffuse term
/// - ks: specular term
/// - rs: specular roughness
/// - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for none)
/// - use_phong: whether to use phong
///
YTRACE_API void set_material(scene* scn, int mid, const float3& ke,
    const float3& kd, const float3& ks, const float3& kt, float rs = 0.1,
    int ke_txt = -1, int kd_txt = -1, int ks_txt = -1, int kt_txt = -1,
    int rs_txt = -1, int norm_txt = -1, bool use_phong = false);

///
/// Sets an environment in the scene.
///
/// Parameters:
/// - scn: scene
/// - mid: material id
/// - frame: local-to-world frame (x, y, z, o in column major order)
/// - ke: emission
/// - ke_txt: emission texture (-1 for none)
///
YTRACE_API void set_environment(scene* scn, int eid, const float3x4& frame,
    const float3& ke, int txt_id = -1);

///
/// Sets a shape in the scene.
///
/// Parameters:
/// - scn: scene
/// - sid: shape id
/// - frame: local-to-world frame (x, y, z, o in column major order)
/// - mid: material id
/// - npoints/nlines/ntriangles: number of elements
/// - points/lines/tiangles: elem data
/// - nverts: number of vertices
/// - pos: vertex positions
/// - norm/tang: vertex normals/tangents
/// - texcoord: vertex texcoord
/// - color: vertex color
/// - tangsp: tangent space for normal and bump mapping
///
YTRACE_API void set_triangle_shape(scene* scn, int sid, const float3x4& frame,
    int mid, int ntriangles, const int3* triangles, int nverts,
    const float3* pos, const float3* norm, const float2* texcoord = nullptr,
    const float3* color = nullptr, const float4* tangsp = nullptr);

///
/// Sets a shape in the scene.
///
/// Parameters:
/// - scn: scene
/// - sid: shape id
/// - frame: local-to-world frame (x, y, z, o in column major order)
/// - mid: material id
/// - npoints/nlines/ntriangles: number of elements
/// - points/lines/tiangles: elem data
/// - nverts: number of vertices
/// - pos: vertex positions
/// - norm/tang: vertex normals/tangents
/// - texcoord: vertex texcoord
/// - color: vertex color
/// - radius: vertex radius
///
YTRACE_API void set_point_shape(scene* scn, int sid, const float3x4& frame,
    int mid, int npoints, const int* points, int nverts, const float3* pos,
    const float3* norm, const float2* texcoord = nullptr,
    const float3* color = nullptr, const float* radius = nullptr);

///
/// Sets a shape in the scene.
///
/// Parameters:
/// - scn: scene
/// - sid: shape id
/// - frame: local-to-world frame (x, y, z, o in column major order)
/// - mid: material id
/// - npoints/nlines/ntriangles: number of elements
/// - points/lines/tiangles: elem data
/// - nverts: number of vertices
/// - pos: vertex positions
/// - norm/tang: vertex normals/tangents
/// - texcoord: vertex texcoord
/// - color: vertex color
/// - radius: vertex radius
///
YTRACE_API void set_line_shape(scene* scn, int sid, const float3x4& frame,
    int mid, int nlines, const int2* lines, int nverts, const float3* pos,
    const float3* tang, const float2* texcoord = nullptr,
    const float3* color = nullptr, const float* radius = nullptr);

///
/// Sets per-vertex material properties.
///
/// Parameters:
/// - scn: scene
/// - sid: shape id
/// - ke: per-vertex emission
/// - kd: per-vertex diffuse
/// - ks: per-vertex specular
/// - rs: per-vertex roughness
///
YTRACE_API void set_vert_material(scene* scn, int sid, const float3* ke,
    const float3* kd, const float3* ks, const float* rs);

///
/// Convert a Phong exponent to GGX/Phong roughness
///
YTRACE_API float specular_exponent_to_roughness(float n);

///
/// Estimates the fresnel coefficient es from ks at normal incidence
///
YTRACE_API void specular_fresnel_from_ks(
    const float3& ks, float3& es, float3& esk);

///
/// Ray-scene Intersection.
///
struct intersect_point {
    /// ray distance
    float dist = 0;
    /// shape index
    int sid = -1;
    /// element distance
    int eid = -1;
    /// element baricentric coordinates
    float3 euv = {0, 0, 0};

    /// check whether it was a hit
    operator bool() const { return eid >= 0; }
};

///
/// Ray-scene closest intersection callback.
///
/// Parameters:
/// - ctx: context
/// - o: ray origin
/// - d: ray direction
/// - tmin/tmax: ray min/max distance
///
/// Return:
/// - intersection point
///
using intersect_first_cb = intersect_point (*)(
    void* ctx, const float3& o, const float3& d, float tmin, float tmax);

///
/// Ray-scene intersection callback
///
/// Parameters:
/// - ctx: context
/// - o: ray origin
/// - d: ray direction
/// - tmin/tmax: ray min/max distance
///
/// Return:
/// - whether we intersect or not
///
using intersect_any_cb = bool (*)(
    void* ctx, const float3& o, const float3& d, float tmin, float tmax);

///
/// Sets the intersection callbacks
///
YTRACE_API void set_intersection_callbacks(scene* scn, void* ctx,
    intersect_first_cb intersect_first, intersect_any_cb intersect_any);

///
/// Logger callback
///
using logging_msg_cb = void (*)(
    int level, const char* name, const char* msg, va_list args);

///
/// Logger
///
YTRACE_API void set_logging_callbacks(
    scene* scn, void* ctx, logging_msg_cb log);

///
/// Initialize rendering.
///
/// Parameters:
/// - scn: trace scene
///
YTRACE_API void init_lights(scene* scn);

///
/// Type of rendering algorithm (shader)
///
enum struct shader_type {
    /// default renderer
    def = 0,
    /// eye hight for quick previews
    eyelight,
    /// direct illumination
    direct,
    /// direct illumination with ambient occlusion
    direct_ao,
    /// pathtrace
    pathtrace,
};

///
/// Random number generator type
///
enum struct rng_type {
    /// default generator
    def = 0,
    /// uniform random numbers
    uniform,
    /// stratified random numbers
    stratified,
    /// correlated multi-jittered sampling
    cmjs,
};

///
/// Rendering params
///
struct render_params {
    /// camera id
    int camera_id = 0;
    /// number of samples
    int nsamples = 256;
    /// progressive rendering
    bool progressive = true;
    /// sampler type
    shader_type stype = shader_type::def;
    /// random number generation type
    rng_type rtype = rng_type::def;
    /// ambient lighting
    float3 amb = {0, 0, 0};
    /// view environment map
    bool envmap_invisible = false;
    /// minimum ray depth
    int min_depth = 3;
    /// maximum ray depth
    int max_depth = 8;
    /// final pixel clamping
    float pixel_clamp = 10;
    /// ray intersection epsilon
    float ray_eps = 1e-4f;
};

///
/// Renders a block of sample
///
/// Parameters:
/// - scn: trace scene
/// - cid: camera id
/// - width: image width
/// - height: image height
/// - img: pixel data in RGBA format
/// - nsamples: number of samples
/// - block_x, block_y: block corner
/// - block_width, block_height: block width and height
/// - samples_min, samples_max: sample block to render [sample_min, sample_max];
/// max values are excluded
///
/// Notes: It is safe to call the function in parallel one different blocks.
/// But two threads should not access the same pixels at the same time.
/// Also blocks with different samples should be called sequentially if
/// accumulate is true.
///
YTRACE_API void trace_block(const scene* scn, int width, int height,
    float4* img, int block_x, int block_y, int block_width, int block_height,
    int samples_min, int samples_max, const render_params& params);

///
/// Convenience function to call trace_block with all sample at once.
///
YTRACE_API void trace_image(const scene* scn, const int width, int height,
    float4* img, const render_params& params);

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YTRACE_INLINE
#include "yocto_trace.cpp"
#endif

#endif
