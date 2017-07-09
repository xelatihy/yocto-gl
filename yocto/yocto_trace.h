///
/// # Yocto/Trace
///
/// Path tracer implementation for support for textured mesh
/// lights, GGX/Phong materials, environment mapping. The interface supports
/// progressive parallel execution with any splitting strategy, by
/// generating per-sample random number sequences or fully deterministic
/// hash-based sampling.
///
/// The raytraced scene is a list of instances of basic shapes. Each shape
/// is a collection of points, lines or triangles with associated normals.
/// Shapes are instanced by creating instances with soecific local-to-world
/// trasforms. Instancing shares memory so large scenes can be created easily.
/// Shaope data is in fact shared with the application and not copied
/// internally.
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
/// sequence per path. This means you can run the path tracer in any order
/// serially or in parallel.
///
/// For now, we support a straightforward path tracer with explicit direct
/// illumination using MIS. Also added simpler shaders for a quick preview
/// and a direct-only renderer.
///
/// The library can support raytracing either by building an internal
/// acceleration structure with Yocto/Bvh or with user supplied intersection
/// routines for custom intersection.
///
/// This library depends in yocto_math.h. Optionally depend on yocto_bvh.h/.cpp
/// for internal acceleration. Disable this by setting YTRACE_NO_BVH.
///
///
/// ## Usage for Scene Creation
///
/// 1. create a scene with `make_scene()`
/// 2. add cameras with `add_camera()`, `set_camera()`
/// 3. add add texture with `add_texture()`
/// 4. create material with `add_XXX_material()`
/// 5. add shapes with `add_XXX_shape()`
/// 6. add instances with `add_instance()`
/// 7. add environment maps with `add_environment()`
///
/// ## Usage for Rendering
///
/// 1. either build the ray-tracing acceleration structure with
///   `init_intersection()` or supply your own with
///   `set_intersection_callbacks()`
/// 2. if desired, add logging with `set_logging_callbacks()`
/// 3. prepare lights for rendering `init_lights()`
/// 4. define rendering params with the `render_params` structure
/// 5. render blocks of samples with `trace_block()` or the whole image with
///     `trace_image()`
///
/// ## History
///
/// - v 0.17: move to rgba per-vertex color
/// - v 0.16: use yocto_math in the interface and remove inline compilation
/// - v 0.15: move to add api
/// - v 1.16: internally use yocto_bvh if desired
/// - v 1.15: added gltf/generic material properties (deprecate old interface)
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

#include <array>
#include <cstdarg>
#include <functional>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// Path tracer with MIS support.
///
namespace ytrace {

// convenient typedef for bytes
using byte = unsigned char;

///
/// Trace scene.
///
struct scene;

///
/// Initialize the scene with the proper number of objects.
///
scene* make_scene();

///
/// Free scene.
///
void free_scene(scene* scn);

///
/// Adds a camera in the scene.
///
/// - Parameters:
///     - scn: scene
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - yfov: field of view
///     - aspect: aspect ratio
///     - aperture: lens aperture
///     - focus: focus plane distance (cannot be zero)
/// - Returns:
///     - camera id
///
int add_camera(scene* scn, const ym::frame3f& frame, float yfov, float aspect,
    float aperture = 0, float focus = 1);

///
/// Sets a camera in the scene.
///
/// - Parameters:
///     - scn: scene
///     - cid: camera id
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - yfov: field of view
///     - aspect: aspect ratio
///     - aperture: lens aperture
///     - focus: focus plane distance (cannot be zero)
///
void set_camera(scene* scn, int cid, const ym::frame3f& frame, float yfov,
    float aspect, float aperture = 0, float focus = 1);

///
/// Adds a texture in the scene.
///
/// - Parameters:
///     - scn: scene
///     - tid: texture id
///     - width: width
///     - height: height
///     - ncomp: number of components (1-4)
///     - hdr: hdr pixels
///     - ldr: ldr pixels (sRGB)
/// - Returns:
///     - texture id
///
int add_texture(scene* scn, int width, int height, int ncomp, const float* hdr);

///
/// Sets a texture in the scene.
///
/// - Parameters:
///     - scn: scene
///     - tid: texture id
///     - width: width
///     - height: height
///     - ncomp: number of components (1-4)
///     - hdr: hdr pixels
///     - ldr: ldr pixels (sRGB)
/// - Returns:
///     - texture id
///
int add_texture(scene* scn, int width, int height, int ncomp, const byte* ldr);

///
/// Sets a material in the scene. [DEPRECATED]
///
/// - Parameters:
///     - scn: scene
///     - ke: emission, term
///     - kd: diffuse term
///     - ks: specular term
///     - rs: specular roughness
///     - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
///     none) - use_phong: whether to use phong
/// - Returns:
///     - material id
///
int add_material(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs = 0.1, int ke_txt = -1,
    int kd_txt = -1, int ks_txt = -1, int kt_txt = -1, int rs_txt = -1,
    int norm_txt = -1, bool use_phong = false);

///
/// Sets a material in the scene with the most customization possible.
///
/// - Parameters:
///     - scn: scene
///     - ke: emission, term
///     - kd: diffuse term
///     - ks: specular term
///     - rs: specular roughness
///     - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
///     none) - use_phong: whether to use phong
/// - Returns:
///     - material id
///
int add_material_generic(scene* scn, const ym::vec3f& ke, const ym::vec3f& kd,
    const ym::vec3f& ks, const ym::vec3f& kt, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int kt_txt, int rs_txt, int op_txt, int norm_txt,
    int occ_txt, bool use_phong = false);

///
/// Sets a gltf metallic roughness material.
///
/// - Parameters:
///     - scn: scene
///     - ke: emission, term
///     - kd: diffuse term
///     - ks: specular term
///     - rs: specular roughness
///     - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
///     none)
///     - use_phong: whether to use phong
/// - Returns:
///     - material id
///
int add_material_gltf_metallic_roughness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, float ks, float rs, float op, int ke_txt, int kd_txt,
    int ks_txt, int norm_txt, int occ_txt, bool use_phong = false);

///
/// Sets a gltf metallic specular glossiness.
///
/// - Parameters:
///     - scn: scene
///     - ke: emission, term
///     - kd: diffuse term
///     - ks: specular term
///     - rs: specular roughness
///     - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
///     none)
///    - use_phong: whether to use phong
/// - Returns:
///     - material id
///
int add_material_gltf_specular_glossiness(scene* scn, const ym::vec3f& ke,
    const ym::vec3f& kd, const ym::vec3f& ks, float rs, float op, int ke_txt,
    int kd_txt, int ks_txt, int norm_txt, int occ_txt, bool use_phong = false);

///
/// Sets a gltf emission only material.
///
/// - Parameters:
///     - scn: scene
///     - mid: material id
///     - ke: emission, term
///     - kd: diffuse term
///     - ks: specular term
///     - rs: specular roughness
///     - ke_txt, kd_txt, ks_txt, rs_txt, norm_txt: texture indices (-1 for
///     none)
///     - use_phong: whether to use phong
/// - Returns:
///     - material id
///
int add_material_emission_only(
    scene* scn, const ym::vec3f& ke, int ke_txt, int norm_txt, int occ_txt);

///
/// Sets an environment in the scene.
///
/// - Parameters:
///     - scn: scene
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - ke: emission
///     - ke_txt: emission texture (-1 for none)
/// - Returns:
///     - environment id
///
int add_environment(
    scene* scn, const ym::frame3f& frame, const ym::vec3f& ke, int txt_id = -1);

///
/// Sets a shape in the scene.
///
/// - Parameters:
///     - scn: scene
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - mid: material id
///     - npoints/nlines/ntriangles: number of elements
///     - points/lines/tiangles: elem data
///     - nverts: number of vertices
///     - pos: vertex positions
///     - norm/tang: vertex normals/tangents
///     - texcoord: vertex texcoord
///     - color: vertex color
///     - tangsp: tangent space for normal and bump mapping
/// - Returns:
///     - shape id
///
int add_triangle_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const ym::vec4f* tangsp = nullptr);

///
/// Sets a shape in the scene.
///
/// - Parameters:
///     - scn: scene
///     - sid: shape id
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - mid: material id
///     - npoints/nlines/ntriangles: number of elements
///     - points/lines/tiangles: elem data
///     - nverts: number of vertices
///     - pos: vertex positions
///     - norm/tang: vertex normals/tangents
///     - texcoord: vertex texcoord
///     - color: vertex color
///     - radius: vertex radius
/// - Returns:
///     - shape id
///
int add_point_shape(scene* scn, int npoints, const int* points, int nverts,
    const ym::vec3f* pos, const ym::vec3f* norm,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const float* radius = nullptr);

///
/// Sets a shape in the scene.
///
/// - Parameters:
///     - scn: scene
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - mid: material id
///     - npoints/nlines/ntriangles: number of elements
///     - points/lines/tiangles: elem data
///     - nverts: number of vertices
///     - pos: vertex positions
///     - norm/tang: vertex normals/tangents
///     - texcoord: vertex texcoord
///     - color: vertex color
///     - radius: vertex radius
/// - Returns:
///     - shape id
///
int add_line_shape(scene* scn, int nlines, const ym::vec2i* lines, int nverts,
    const ym::vec3f* pos, const ym::vec3f* tang,
    const ym::vec2f* texcoord = nullptr, const ym::vec4f* color = nullptr,
    const float* radius = nullptr);

///
/// Adds an instance in the scene.
///
/// - Parameters:
///     - scn: scene
///     - iid: instance id
///     - frame: local-to-world frame (x, y, z, o in column major order)
///     - sid: shape id
///     - mid: material id
/// - Returns:
///     - instance id
///
int add_instance(scene* scn, const ym::frame3f& frame, int sid, int mid);

///
/// Sets per-vertex material properties.
///
/// - Parameters:
///     - scn: scene
///     - sid: shape id
///     - ke: per-vertex emission
///     - kd: per-vertex diffuse
///     - ks: per-vertex specular
///     - rs: per-vertex roughness
///
void set_vert_material(scene* scn, int sid, const ym::vec3f* ke,
    const ym::vec3f* kd, const ym::vec3f* ks, const float* rs);

///
/// Convert a Phong exponent to GGX/Phong roughness
///
float specular_exponent_to_roughness(float n);

///
/// Estimates the fresnel coefficient es from ks at normal incidence
///
void specular_fresnel_from_ks(
    const ym::vec3f& ks, ym::vec3f& es, ym::vec3f& esk);

///
/// Ray-scene Intersection.
///
struct intersect_point {
    /// ray distance
    float dist = 0;
    /// instance index
    int iid = -1;
    /// shape index
    int sid = -1;
    /// element distance
    int eid = -1;
    /// element baricentric coordinates
    ym::vec3f euv = {0, 0, 0};

    /// check whether it was a hit
    operator bool() const { return eid >= 0; }
};

///
/// Ray-scene closest intersection callback.
///
/// - Parameters:
///     - ctx: context
///     - ray: ray
/// - Return:
///     - intersection point
///
using intersect_first_cb = intersect_point (*)(void* ctx, const ym::ray3f& ray);

///
/// Ray-scene intersection callback
///
/// - Parameters:
///     - ctx: context
///     - ray: ray
/// - Return:
///     - whether we intersect or not
///
using intersect_any_cb = bool (*)(void* ctx, const ym::ray3f& o);

///
/// Sets the intersection callbacks
///
void set_intersection_callbacks(scene* scn, void* ctx,
    intersect_first_cb intersect_first, intersect_any_cb intersect_any);

///
/// Initialize acceleration structure.
///
/// - Parameters:
///     - scn: trace scene
///
void init_intersection(scene* scn);

///
/// Logger callback
///
using logging_msg_cb = void (*)(
    int level, const char* name, const char* msg, va_list args);

///
/// Logger
///
void set_logging_callbacks(scene* scn, void* ctx, logging_msg_cb log);

///
/// Initialize lighting.
///
/// - Parameters:
///     - scn: trace scene
///
void init_lights(scene* scn);

///
/// Type of rendering algorithm (shader)
///
enum struct shader_type {
    /// pathtrace
    pathtrace = 0,
    /// eye hight for quick previews
    eyelight,
    /// direct illumination
    direct,
    /// direct illumination with ambient occlusion
    direct_ao,
};

///
/// Random number generator type
///
enum struct rng_type {
    /// uniform random numbers
    uniform = 0,
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
    shader_type stype = shader_type::pathtrace;
    /// random number generation type
    rng_type rtype = rng_type::stratified;
    /// ambient lighting
    ym::vec3f amb = {0, 0, 0};
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
/// Notes: It is safe to call the function in parallel one different blocks.
/// But two threads should not access the same pixels at the same time.
/// Also blocks with different samples should be called sequentially if
/// accumulate is true.
///
/// - Parameters:
///     - scn: trace scene
///     - cid: camera id
///     - width: image width
///     - height: image height
///     - img: pixel data in RGBA format
///     - nsamples: number of samples
///     - block_x, block_y: block corner
///     - block_width, block_height: block width and height
///     - samples_min, samples_max: sample block to render
///       [sample_min, sample_max]; max values are excluded
///
void trace_block(const scene* scn, int width, int height, ym::vec4f* img,
    int block_x, int block_y, int block_width, int block_height,
    int samples_min, int samples_max, const render_params& params);

///
/// Convenience function to call trace_block() with all samples at once.
///
void trace_image(const scene* scn, const int width, int height, ym::vec4f* img,
    const render_params& params);

}  // namespace ytrace

#endif
