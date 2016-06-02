//
// YOCTO_TRACE: path tracer implementation for support for textured mesh
// lights, GGX/Phong materials, environment mapping. The interface supports
// progressive parallel execution with any splitting strategy, by
// generating per-sample random number sequences or fully deterministic
// hash-based sampling.
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. define your scene
// - init the scene
//     scene = yt_make_scene(num_shapes, num_materials, num_textures)
// - define shapes, materials and textures (we store pointers, not copies)
//     for(int i = 0; i < nshapes; i ++) yt_set_shape(scene, shape data)
//     for(int i = 0; i < nmaterials; i ++) yt_set_material(scene, mat data)
//     for(int i = 0; i < ntextures; i ++) yt_set_texture(scene, tex data)
// - define camera and environment
//     yt_set_camera(scene, camera data)
//     yt_set_env(scene, env data)
// 2. intersection routines are handled outside of this library, for example
//    in yocto_bvh; to include it in this library, pass two funtion pointers
//    to find the closest intersection and any intersection
//    yt_set_interseciton(scene, intersect_context, intersect_func, hit_func)
// 3. prepare for rendering
//    yt_init_lights(scene)
// 4. render blocks of samples
//    yt_render_block(scene, pixels, image size, block)
// 5. cleanup with yt_free_scene
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles),
// an array of vertex positions, and an array of vertex radius
// (for points and lines).
//
// Materials are represented as sums of an emission term, a diffuse term and
// a specular microfacet term (GGX or Phong). Only opaque for now. We pick
// a proper material type for each shape element type (points, lines,
// triangles).
//
// Lights are defined as any shape with a material emission term. Additionally
// one can also add an environment map. But even if you can, you might want to
// add a large triangle mesh with inward normals instead. The latter is more
// general (you can even more an arbitrary shape sun).
//
// We generate our own random numbers guarantying that there is one random
// sequence per path. This means you can rul the path tracer in any order
// serially or in parallel.
//
// For now, we support a straightforward path tracer with explicit direct
// illumination using MIS.
//

//
// COMPILATION:
//
// The library has two APIs. The default one is usable directly from C++,
// while the other is usable from both C and C++. To use from C, compile the
// library into a static or dynamic lib using a C++ and then include/link from
// C using the C API.
//
// All functions in this library are inlined by default for ease of use in C++.
// To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.1: C++ implementation
// - v 0.0: initial release in C99
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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

#ifndef _YT_H_
#define _YT_H_

#ifdef __cplusplus
#ifndef YGL_DECLARATION
#define YGL_API inline
#define YGLC_API inline
#else
#define YGL_API
#define YGLC_API extern "C"
#endif
#include "yocto_math.h"
#endif

#ifndef __cplusplus
#define YGLC_API extern
#include <stdbool.h>
#endif

// -----------------------------------------------------------------------------
// C++ INTERFACE
// -----------------------------------------------------------------------------

#ifdef __cplusplus

//
// Shape element types
//
enum {
    yt_etype_point = 1,     // points
    yt_etype_line = 2,      // lines
    yt_etype_triangle = 3,  // triangle
};

//
// Type of rendering algorithm (shader)
//
enum {
    yt_stype_default = 0,  // default renderer
    yt_stype_eyelight,     // eye hight for quick previews
    yt_stype_direct,       // direct illumination
    yt_stype_pathtrace,    // path tracing
    yt_stype_max           // max number of renderer
};

//
// Random number generator type
//
enum {
    yt_rtype_default = 0,  // default generator
    yt_rtype_uniform,      // uniform random numbers
    yt_rtype_stratified,   // stratified random numbers
    yt_rtype_cmjs,         // correlated multi-jittered sampling
    yt_rtype_max           // max number of generators
};

//
// Opaque scene data structure
//
typedef struct yt_scene yt_scene;

//
// Rendering params
//
struct yt_render_params {
    int stype = yt_stype_default;  // sampler type
    int rtype = yt_rtype_default;  // random type
    ym_vec3f amb = ym_zero3f;      // ambient lighting
    int min_depth = 3;             // min ray depth
    int max_depth = 8;             // mas ray depth
    float pixel_clamp = 100;       // final pixel clamping
    float ray_eps = 1e-2f;         // ray intersection epsilon
};

//
// Ray-scene closest intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Out Parameters:
// - ray_t: hit distance
// - sid: hit shape index
// - eid: hit element index
// - euv: hit element parameters
//
// Return:
// - whether we intersect or not
//
typedef bool (*yt_intersect_ray)(const void* ctx, const ym_ray3f& ray,
                                 float* ray_t, int* sid, int* eid,
                                 ym_vec2f* euv);

//
// Ray-scene closest intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Return:
// - whether we intersect or not
//
typedef bool (*yt_hit_ray)(const void* ctx, const ym_ray3f& ray);

//
// Initialize the scene data with ncameras, nshapes, nmaterials, ntextures,
// nenvmaps.
//
YGL_API yt_scene* yt_make_scene(int ncameras, int nshapes, int nmaterials,
                                int ntextures, int nenvs);

//
// Cleanup scene data
//
YGL_API void yt_free_scene(yt_scene* scene);

//
// Set intersection callback. See callback description above.
//
YGL_API void yt_set_intersection(yt_scene* scene, void* ray_ctx,
                                 yt_intersect_ray intersect_ray,
                                 yt_hit_ray hit_ray);

//
// Set rendering camera.
//
// Parameters:
// - scene: trace scene
// - cid: camera id
// - xform: camera local to world **rigid** transform
// - width, height: width and height of image plane
// - aperture: aperture
// - focus: focus distance
//
YGL_API void yt_set_camera(yt_scene* scene, int cid, const ym_frame3f& xform,
                           float width, float height, float aperture,
                           float focus);

//
// Set rendering environment map.
//
// Parameters:
// - scene: trace scene
// - eid: envmap id
// - xform: env local to world **rigid** transform
// - ke: emission
// - ke_txt: ke texture index (-1 for no texture)
//
YGL_API void yt_set_env(yt_scene* scene, int eid, const ym_frame3f& xform,
                        const ym_vec3f& ke, int ke_txt);

//
// Set a shape.
//
// Parameters:
// - scene: trace scene
// - sid: shape index
// - mid: material index for this shape
// - xform: shape local to world **rigid** transform
// - nelems: number of elements
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - nverts: number of vertices
// - pos: array of vertex positions
// - norm: array of vertex normals
// - texcoord: optional array of vertex texture coordinates
// - color: optional array of vertex colors (rgb)
// - radius: optional array of vertex radius
//
YGL_API void yt_set_shape(yt_scene* scene, int sid, int mid,
                          const ym_frame3f& xform, int nelems, const int* elem,
                          int etype, int nverts, const ym_vec3f* pos,
                          const ym_vec3f* norm, const ym_vec2f* texcoord,
                          const ym_vec3f* color, const float* radius);

//
// Set a material.
//
// Parameters:
// - scene: trace scene
// - mid: material index for this shape
// - ke: emission
// - kd: diffuse
// - ks: specular
// - rs: specular roughness
// - es: optional index of refraction (eta) that controls the fresnel term;
// NULL to disable fresnel; for materals use eta[0-2] for real part and
// eta[3-5] for complex part; for dielectics use eta[0-2] for ior and eta[3-5]=0
// - ke_txt: ke texture index (-1 for no texture)
// - kd_txt: kd texture index (-1 for no texture)
// - ks_txt: ks texture index (-1 for no texture)
// - rs_txt: rs texture index (-1 for no texture)
// - use_phong: whether to use Phong (normaly use false here)
//
YGL_API void yt_set_material(yt_scene* scene, int mid, const ym_vec3f& ke,
                             const ym_vec3f& kd, const ym_vec3f& ks, float rs,
                             const ym_vec3f& es, const ym_vec3f& esk,
                             int ke_txt, int kd_txt, int ks_txt, int rs_txt,
                             bool use_phong);

//
// Convert a Phong exponent to GGX/Phong roughness
//
YGL_API float yt_specular_exponent_to_roughness(float n);

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YGL_API void yt_specular_fresnel_from_ks(const ym_vec3f& ks, ym_vec3f* es,
                                         ym_vec3f* esk);

//
// Set a texture.
//
// Parameters:
// - scene: trace scene
// - tid: texture id
// - pixel: pixel data
// - w, h: width and height
// - nc: number of components (3-4 supported for now)
//
YGL_API void yt_set_texture(yt_scene* scene, int tid, const float* pixels,
                            int w, int h, int nc);

//
// Set rendering camera.
//
// Parameters:
// - scene: trace scene
// - xform: camera local to world **rigid** transform
// - width, height: width and height of image plane
// - aperture: aperture
// - focus: focus distance
//
YGL_API void yt_init_lights(yt_scene* scene);

//
// Renders a block of sample
//
// Parameters:
// - scene: trace scene
// - cid: camera id
// - pixels: pixel data in RGBA format
// - w, h: image width and height
// - ns: number of samples
// - block: image block to render [xmin, xmax, ymin, ymax];
// max values are excluded
// - sampples: sample block to render [sample_min, sample_max];
// max values are excluded
//
// Notes: It is safe to call the function in parallel one different blocks.
// But two threads should not access the same pixels at the same time.
// Also blocks with different samples should be called sequentially.
//
YGL_API void yt_trace_block(const yt_scene* scene, int cid, ym_vec4f* pixels,
                            int w, int h, int ns, const ym_range2i& block,
                            const ym_range1i& samples,
                            const yt_render_params& params);

//
// Convenience function to call yt_trace_block with all sample at once.
//
YGL_API void yt_trace_image(const yt_scene* scene, int cid, ym_vec4f* pixels,
                            int w, int h, int ns,
                            const yt_render_params& params);

#endif  // __cplusplus

// -----------------------------------------------------------------------------
// C/C++ INTERFACE
// -----------------------------------------------------------------------------

//
// Rendering params
//
typedef struct ytc_render_params {
    int stype;          // sampler type
    int rtype;          // random type
    float amb[3];       // ambient lighting
    int min_depth;      // min ray depth
    int max_depth;      // mas ray depth
    float pixel_clamp;  // final pixel clamping
    float ray_eps;      // ray intersection epsilon
} ytc_render_params;

//
// Ray-scene closest intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Out Parameters:
// - ray_t: hit distance
// - sid: hit shape index
// - eid: hit element index
// - euv: hit element parameters
//
// Return:
// - whether we intersect or not
//
typedef bool (*ytc_intersect_ray)(void* ctx, const float ray_o[3],
                                  const float ray_d[3], float ray_min,
                                  float ray_max, int ray_mask, float* ray_t,
                                  int* sid, int* eid, float* euv);

//
// Ray-scene closest intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Return:
// - whether we intersect or not
//
typedef bool (*ytc_hit_ray)(void* ctx, const float ray_o[3],
                            const float ray_d[3], float ray_min, float ray_max,
                            int ray_mask);

//
// Initialize the scene data with nshapes, nmaterials, ntextures.
//
YGLC_API yt_scene* ytc_make_scene(int ncameras, int nshapes, int nmaterials,
                                  int ntextures, int nenvs);

//
// Cleanup scene data
//
YGLC_API void ytc_free_scene(yt_scene* scene);

//
// Set intersection callback. See callback description above.
//
YGLC_API void ytc_set_intersection(yt_scene* scene, void* ray_ctx,
                                   yt_intersect_ray intersect_ray,
                                   yt_hit_ray hit_ray);

//
// Set rendering camera.
//
// Parameters:
// - scene: trace scene
// - cid: camera id
// - xform: camera local to world **rigid** transform
// - width, height: width and height of image plane
// - aperture: aperture
// - focus: focus distance
//
YGLC_API void ytc_set_camera(yt_scene* scene, int cid, float xform[16],
                             float width, float height, float aperture,
                             float focus);

//
// Set rendering environment map.
//
// Parameters:
// - scene: trace scene
// - eid: envmap id
// - xform: env local to world **rigid** transform
// - ke: emission
// - ke_txt: ke texture index (-1 for no texture)
//
YGL_API void ytc_set_env(yt_scene* scene, int eid, float xform[16], float ke[3],
                         int ke_txt);

//
// Set rendering camera.
//
// Parameters:
// - scene: trace scene
// - xform: camera local to world **rigid** transform
// - width, height: width and height of image plane
// - aperture: aperture
// - focus: focus distance
//
YGL_API void ytc_init_lights(yt_scene* scene);

//
// Set a shape.
//
// Parameters:
// - scene: trace scene
// - sid: shape index
// - mid: material index for this shape
// - xform: shape local to world **rigid** transform
// - nelems: number of elements
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - nverts: number of vertices
// - pos: array of vertex positions
// - norm: array of vertex normals
// - texcoord: optional array of vertex texture coordinates
// - color: optional array of vertex colors (rgb)
// - radius: optional array of vertex radius
//
YGL_API void ytc_set_shape(yt_scene* scene, int sid, int mid, float xform[16],
                           int nelems, int* elem, int etype, int nverts,
                           float* pos, float* norm, float* texcoord,
                           float* color, float* radius);

//
// Set a material.
//
// Parameters:
// - scene: trace scene
// - mid: material index for this shape
// - ke: emission
// - kd: diffuse
// - ks: specular
// - rs: specular roughness
// - es: optional index of refraction (eta) that controls the fresnel term;
// NULL to disable fresnel; for materals use eta[0-2] for real part and
// eta[3-5] for complex part; for dielectics use eta[0-2] for ior and eta[3-5]=0
// - ke_txt: ke texture index (-1 for no texture)
// - kd_txt: kd texture index (-1 for no texture)
// - ks_txt: ks texture index (-1 for no texture)
// - rs_txt: rs texture index (-1 for no texture)
// - use_phong: whether to use Phong (normaly use false here)
//
YGL_API void ytc_set_material(yt_scene* scene, int mid, float ke[3],
                              float kd[3], float ks[3], float rs, float es[3],
                              float esk[3], int ke_txt, int kd_txt, int ks_txt,
                              int rs_txt, bool use_phong);

//
// Convert a Phong exponent to GGX/Phong roughness
//
YGL_API float ytc_specular_exponent_to_roughness(float n);

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YGL_API void ytc_specular_fresnel_from_ks(const float ks[3], float es[3],
                                          float esk[3]);

//
// Set a texture.
//
// Parameters:
// - scene: trace scene
// - tid: texture id
// - pixel: pixel data
// - w, h: width and height
// - nc: number of components (3-4 supported for now)
//
YGL_API void ytc_set_texture(yt_scene* scene, int tid, float* pixels, int w,
                             int h, int nc);

//
// Set rendering parameters
//
// Parameters:
// - scene: trace scene
// - stype: rendering algorithm type
// - rtype: random number generator type
// - amb: ambient lighting for direct rendering algorithm
//
YGL_API ytc_render_params ytc_make_rendering_params();

//
// Renders a block of sample
//
// Parameters:
// - scene: trace scene
// - camera: camera id
// - pixels: pixel data in RGBA format
// - w, h: image width and height
// - ns: number of samples
// - block: image block to render [xmin, xmax, ymin, ymax];
// max values are excluded
// - sampples: sample block to render [sample_min, sample_max];
// max values are excluded
//
// Notes: It is safe to call the function in parallel one different blocks.
// But two threads should not access the same pixels at the same time.
// Also blocks with different samples should be called sequentially.
//
YGL_API void ytc_trace_block(const yt_scene* scene, int cid, float* pixels,
                             int w, int h, int ns, const int block[4],
                             const int samples[2],
                             const ytc_render_params* params);

//
// Convenience function to call yt_trace_block with all sample at once.
//
YGL_API void ytc_trace_image(const yt_scene* scene, int cid, float* pixels,
                             int w, int h, int ns,
                             const ytc_render_params* params);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if defined(__cplusplus) &&                                                    \
    (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

//
// Camera
//
struct yt__camera {
    ym_frame3f xform = ym_identity_frame3f;  // local-to-world transform
    float width = 1, height = 1;             // image width and height
    float aperture = 0, focus = 1;  // lens aperture and focus plane distance
};

//
// Texture
//
struct yt__texture {
    int w = 0, h = 0, nc = 0;       // width, height, number of components
    const float* pixels = nullptr;  // pixel data
};

//
// Material
//
struct yt__material {
    ym_vec3f ke = ym_zero3f;   // emission, term
    ym_vec3f kd = ym_zero3f;   // diffuse term
    ym_vec3f ks = ym_zero3f;   // specular term
    float rs = 0.1;            // specular roughness API
    ym_vec3f es = ym_zero3f;   // eta
    ym_vec3f eks = ym_zero3f;  // etak (metals only)
    yt__texture *ke_txt = nullptr, *kd_txt = nullptr, *ks_txt = nullptr,
                *rs_txt = nullptr;  // textures
    bool use_phong = false;         // whether to use phong
};

//
// Shape
//
struct yt__shape {
    ym_frame3f xform = ym_identity_frame3f;  // local-to-world rigid transform
    yt__material* mat = nullptr;             // material
    int nelems = 0, etype = 0, nverts = 0;   // element data
    const int* elem = nullptr;               // element data
    const ym_vec3f* pos = nullptr;           // vertex data
    const ym_vec3f* norm = nullptr;          // vertex data
    const ym_vec2f* texcoord = nullptr;      // vertex data
    const ym_vec3f* color = nullptr;         // vertex data
    const float* radius = nullptr;           // vertex data
};

//
// Environment
//
struct yt__env {
    ym_frame3f xform = ym_identity_frame3f;  // local-to-world rigid transform
    ym_frame3f xform_inv =
        ym_identity_frame3f;        // local-to-world rigid transform
    bool xformed = false;           // transformed
    ym_vec3f ke = ym_zero3f;        // emission
    yt__texture* ke_txt = nullptr;  // emission texture
};

//
// Light (either shape or environment)
//
struct yt__light {
    yt__shape* shape = nullptr;  // shape
    yt__env* env = nullptr;      // environment
    ym_vector<float> cdf;  // for shape, cdf of shape elements for sampling
    float area = 0;        // for shape, shape area
};

//
// Scene
//
struct yt_scene {
    const void* ray_ctx = nullptr;             // ray context
    yt_intersect_ray intersect_ray = nullptr;  // ray intersection callback
    yt_hit_ray hit_ray = nullptr;              // ray hit callback

    ym_vector<yt__camera> camera;       // camera
    ym_vector<yt__env> env;             // env
    ym_vector<yt__shape> shapes;        // shapes
    ym_vector<yt__material> materials;  // materials
    ym_vector<yt__texture> textures;    // textures

    ym_vector<yt__light> lights;  // lights
};

//
// Init scene. Public API, see above.
//
YGL_API yt_scene* yt_make_scene(int ncameras, int nshapes, int nmaterials,
                                int ntextures, int nenvs) {
    yt_scene* scene = new yt_scene();
    scene->camera.resize(ncameras);
    scene->shapes.resize(nshapes);
    scene->materials.resize(nmaterials);
    scene->textures.resize(ntextures);
    scene->env.resize(nenvs);
    return scene;
}

//
// Free scene. Public API, see above.
//
YGL_API void yt_free_scene(yt_scene* scene) { delete scene; }

//
// Set intersection callbacks. Public API, see above.
//
YGL_API void yt_set_intersection(yt_scene* scene, void* ray_ctx,
                                 yt_intersect_ray intersect_ray,
                                 yt_hit_ray hit_ray) {
    scene->ray_ctx = ray_ctx;
    scene->intersect_ray = intersect_ray;
    scene->hit_ray = hit_ray;
}

//
// Set camera. Public API, see above.
//
YGL_API void yt_set_camera(yt_scene* scene, int cid, const ym_frame3f& xform,
                           float width, float height, float aperture,
                           float focus) {
    yt__camera* cam = &scene->camera[cid];
    cam->xform = xform;
    cam->width = width;
    cam->height = height;
    cam->aperture = aperture;
    cam->focus = (focus) ? focus : 1;
}

//
// Set shape. Public API, see above.
//
YGL_API void yt_set_shape(yt_scene* scene, int sid, int matid,
                          const ym_frame3f& xform, int nelems, const int* elem,
                          int etype, int nverts, const ym_vec3f* pos,
                          const ym_vec3f* norm, const ym_vec2f* texcoord,
                          const ym_vec3f* color, const float* radius) {
    assert(matid >= 0 && matid < scene->materials.size());
    yt__shape* shape = &scene->shapes[sid];
    shape->xform = xform;
    shape->mat = (matid < 0) ? 0 : &scene->materials[matid];
    shape->nelems = nelems;
    shape->elem = elem;
    shape->etype = etype;
    shape->nverts = nverts;
    shape->pos = pos;
    shape->norm = norm;
    shape->texcoord = texcoord;
    shape->color = color;
    shape->radius = radius;
}

//
// Compute shape element cdf for shape sampling.
//
static inline void yt__compute_shape_cdf(int nelems, const int* elem, int etype,
                                         const ym_vec3f* pos, float* cdf,
                                         float* area) {
    switch (etype) {
        case yt_etype_point: {
            for (int i = 0; i < nelems; i++) cdf[i] = i + 1;
        } break;
        case yt_etype_line: {
            for (int i = 0; i < nelems; i++) {
                const int* l = elem + i * 2;
                cdf[i] = ym_dist(pos[l[0]], pos[l[1]]);
                if (i != 0) cdf[i] += cdf[i - 1];
            }
        } break;
        case yt_etype_triangle: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 3;
                const ym_vec3f *v0 = pos + f[0], *v1 = pos + f[1],
                               *v2 = pos + f[2];
                ym_vec3f e1 = *v1 - *v0;
                ym_vec3f e2 = *v2 - *v0;
                ym_vec3f a = ym_cross(e1, e2);
                cdf[i] = ym_length(a) / 2;
                if (i != 0) cdf[i] += cdf[i - 1];
            }
        } break;
        default: assert(false); break;
    }

    // normalize
    *area = cdf[nelems - 1];
    for (int i = 0; i < nelems; i++) cdf[i] /= *area;
}

//
// Init lights. Public API, see above.
//
YGL_API void yt_init_lights(yt_scene* scene) {
    // clear old lights
    scene->lights.resize(0);

    for (int sid = 0; sid < scene->shapes.size(); sid++) {
        yt__shape* shape = &scene->shapes[sid];
        if (ym_iszero(shape->mat->ke)) continue;
        scene->lights.push_back(yt__light());
        yt__light* light = &scene->lights.back();
        light->shape = shape;
        light->cdf.resize(shape->nelems);
        yt__compute_shape_cdf(shape->nelems, shape->elem, shape->etype,
                              (ym_vec3f*)shape->pos, light->cdf.data(),
                              &light->area);
    }

    for (int sid = 0; sid < scene->env.size(); sid++) {
        yt__env* env = &scene->env[sid];
        if (ym_iszero(env->ke)) continue;
        scene->lights.push_back(yt__light());
        yt__light* light = &scene->lights.back();
        light->env = &scene->env[0];
    }
}

//
// Set material. Public API, see above.
//
YGL_API void yt_set_material(yt_scene* scene, int mid, const ym_vec3f& ke,
                             const ym_vec3f& kd, const ym_vec3f& ks, float rs,
                             const ym_vec3f& es, const ym_vec3f& esk,
                             int ke_txt, int kd_txt, int ks_txt, int rs_txt,
                             bool use_phong) {
    yt__material* material = &scene->materials[mid];
    material->ke = ke;
    material->kd = kd;
    material->ks = ks;
    material->es = es;
    material->eks = esk;
    material->rs = rs;
    material->ke_txt = (ke_txt < 0) ? 0 : &scene->textures[ke_txt];
    material->kd_txt = (kd_txt < 0) ? 0 : &scene->textures[kd_txt];
    material->ks_txt = (ks_txt < 0) ? 0 : &scene->textures[ks_txt];
    material->rs_txt = (rs_txt < 0) ? 0 : &scene->textures[rs_txt];
    material->use_phong = use_phong;
}

//
// Phong exponent to roughness. Public API, see above.
//
YGL_API float yt_specular_exponent_to_roughness(float n) {
    return sqrtf(2 / (n + 2));
}

//
// Specular to fresnel eta. Public API, see above.
//
YGL_API void yt_specular_fresnel_from_ks(const ym_vec3f& ks, ym_vec3f* es,
                                         ym_vec3f* esk) {
    *es = ym_vec3f{(1 + sqrtf(ks.x)) / (1 - sqrtf(ks.x)),
                   (1 + sqrtf(ks.y)) / (1 - sqrtf(ks.y)),
                   (1 + sqrtf(ks.z)) / (1 - sqrtf(ks.z))};
    *esk = ym_zero3f;
}

//
// Set env. Public API, see above.
//
YGL_API void yt_set_env(yt_scene* scene, int eid, const ym_frame3f& xform,
                        const ym_vec3f& ke, int ke_txt) {
    yt__env* env = &scene->env[eid];
    env->xform = xform;
    env->xform_inv = ym_inverse(env->xform);
    env->ke = ke;
    env->ke_txt = (ke_txt < 0) ? 0 : &scene->textures[ke_txt];
}

//
// Set texture. Public API, see above.
//
YGL_API void yt_set_texture(yt_scene* scene, int tid, const float* pixels,
                            int w, int h, int nc) {
    yt__texture* txt = &scene->textures[tid];
    txt->w = w;
    txt->h = h;
    txt->nc = nc;
    txt->pixels = pixels;
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// Random number sampler. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
//
struct yt__sampler {
    ym_rng_pcg32 rng;  // rnumber number state
    int i, j;          // pixel coordinates
    int s, d;          // sample and dimension indices
    int ns;            // number of samples
    int rtype;         // random number type
};

//
// Initialize a sampler ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
//
static inline yt__sampler yt__make_sampler(int i, int j, int s, int ns,
                                           int rtype) {
    // we use various hashes to scramble the pixel values
    yt__sampler sampler = {{0, 0}, i, j, s, 0, ns, rtype};
    uint64_t sample_id = ((uint64_t)(i + 1)) << 0 | ((uint64_t)(j + 1)) << 15 |
                         ((uint64_t)(s + 1)) << 30;
    uint64_t initseq = ym_hash_uint64(sample_id);
    uint64_t initstate =
        ym_hash_uint64(sample_id * 3202034522624059733ull + 1ull);
    ym_rng_init(&sampler.rng, initstate, initseq);
    return sampler;
}

//
// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
//
static inline float yt__sample_next1f(yt__sampler* sampler) {
    float rn = 0;
    switch (sampler->rtype) {
        case yt_rtype_default:
        case yt_rtype_uniform: {
            rn = ym_rng_nextf(&sampler->rng);
        } break;
        case yt_rtype_stratified: {
            uint32_t p = ym_hash_uint64_32(((uint64_t)(sampler->i + 1)) << 0 |
                                           ((uint64_t)(sampler->j + 1)) << 15 |
                                           ((uint64_t)(sampler->d + 1)) << 30);
            int s = ym_hash_permute(sampler->s, sampler->ns, p);
            rn = (s + ym_rng_nextf(&sampler->rng)) / sampler->ns;
        } break;
        case yt_rtype_cmjs: {
            uint32_t p = ym_hash_uint64_32(((uint64_t)(sampler->i + 1)) << 0 |
                                           ((uint64_t)(sampler->j + 1)) << 15 |
                                           ((uint64_t)(sampler->d + 1)) << 30);
            int s = ym_hash_permute(sampler->s, sampler->ns, p);
            rn = (s + ym_hash_randfloat(s, p * 0xa399d265)) / sampler->ns;
        } break;
        default: assert(false);
    }

    sampler->d += 1;

    // make sure all sampled numbers are below 1
    if (rn >= 1) rn = 1 - FLT_EPSILON;

    return rn;
}

//
// Generates a 1-dimensional sample.
//
// Implementation notes: see above. Note that using deterministic keyed
// permutaton we can use stratified sampling without preallocating samples.
//
static inline ym_vec2f yt__sample_next2f(yt__sampler* sampler) {
    ym_vec2f rn = {0, 0};
    switch (sampler->rtype) {
        case yt_rtype_default:
        case yt_rtype_uniform: {
            rn.x = ym_rng_nextf(&sampler->rng);
            rn.y = ym_rng_nextf(&sampler->rng);
        } break;
        case yt_rtype_stratified: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler->ns));
            uint32_t p = ym_hash_uint64_32(((uint64_t)(sampler->i + 1)) << 0 |
                                           ((uint64_t)(sampler->j + 1)) << 15 |
                                           ((uint64_t)(sampler->d + 1)) << 30);
            int s = ym_hash_permute(sampler->s, sampler->ns, p);
            rn.x = (s % ns2 + ym_rng_nextf(&sampler->rng)) / ns2;
            rn.y = (s / ns2 + ym_rng_nextf(&sampler->rng)) / ns2;
        } break;
        case yt_rtype_cmjs: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler->ns));
            uint32_t p = ym_hash_uint64_32(((uint64_t)(sampler->i + 1)) << 0 |
                                           ((uint64_t)(sampler->j + 1)) << 15 |
                                           ((uint64_t)(sampler->d + 1)) << 30);
            int s = ym_hash_permute(sampler->s, sampler->ns, p);
            int sx = ym_hash_permute(s % ns2, ns2, p * 0xa511e9b3);
            int sy = ym_hash_permute(s / ns2, ns2, p * 0x63d83595);
            float jx = ym_hash_randfloat(s, p * 0xa399d265);
            float jy = ym_hash_randfloat(s, p * 0x711ad6a5);
            rn.x = (s % ns2 + (sy + jx) / ns2) / ns2;
            rn.y = (s / ns2 + (sx + jy) / ns2) / ns2;
        } break;
        default: assert(false);
    }

    sampler->d += 2;

    // make sure all sampled numbers are below 1
    if (rn.x >= 1) rn.x = 1 - FLT_EPSILON;
    if (rn.y >= 1) rn.y = 1 - FLT_EPSILON;

    return rn;
}

//
// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
//
struct yt__point {
    // shape data -------------------------
    const yt__shape* shape;
    int eid;       // element id
    ym_vec2f euv;  // element baricentric coordinates
    // env data ---------------------------
    const yt__env* env;  // env
    // direction --------------------------
    ym_vec3f wo;  // outgoing direction
    // resolved geometry (shape) ----------
    int etype;         // element type
    ym_frame3f frame;  // local frame
    // shading ----------------------------
    ym_vec3f ke, kd, ks;  // material values
    float rs;             // material values
    ym_vec3f es, eks;     // material values
    bool use_phong;       // material values
    // sampling data ----------------------
    float wkd, wks, wrr;  // weights of kd and ks and total
};

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static inline ym_ray3f yt__eval_camera(const yt__camera* cam,
                                       const ym_vec2f& uv,
                                       const ym_vec2f& luv) {
    ym_vec3f o = ym_vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    ym_vec3f q = {cam->width * cam->focus * (uv.x - 0.5f),
                  cam->height * cam->focus * (uv.y - 0.5f), -cam->focus};
    return ym_ray3f(ym_transform_point(cam->xform, o),
                    ym_transform_direction(cam->xform, ym_normalize(q - o)));
}

//
// Evaluates a texture txt at texture coordinates texcoord. Uses bilinear
// interpolation and tiling. If the texture is semitrasparent, uses c
// as the case color.
//
static inline ym_vec3f yt__eval_texture(const ym_vec3f& c,
                                        const yt__texture* txt,
                                        const ym_vec2f& texcoord) {
    if (!txt) return c;
    assert(txt->nc == 3 || txt->nc == 4);
    ym_vec2f tc = {fmodf(texcoord.x, 1), fmodf(texcoord.y, 1)};
    if (tc.x < 0) tc.x += 1;
    if (tc.y < 0) tc.y += 1;
    int i = tc.x * txt->w, j = tc.y * txt->h;
    if (i > txt->w - 1) i = txt->w - 1;
    if (j > txt->h - 1) i = txt->h - 1;
    int ii = (i + 1) % txt->w, jj = (j + 1) % txt->h;
    float s = tc.x * txt->w - i, t = tc.y * txt->h - j;
    int idx[4] = {j * txt->w + i, jj * txt->w + i, j * txt->w + ii,
                  jj * txt->w + ii};
    float w[4] = {(1 - s) * (1 - t), (1 - s) * t, s * (1 - t), s * t};
    if (txt->nc == 3) {
        ym_vec3f* p3 = (ym_vec3f*)txt->pixels;
        return p3[idx[0]] * w[0] + p3[idx[1]] * w[1] + p3[idx[2]] * w[2] +
               p3[idx[3]] * w[3];
    } else if (txt->nc == 4) {
        ym_vec4f* p4 = (ym_vec4f*)txt->pixels;
        ym_vec4f t = p4[idx[0]] * w[0] + p4[idx[1]] * w[1] + p4[idx[2]] * w[2] +
                     p4[idx[3]] * w[3];
        return ym_vec3f{t.x, t.y, t.z} * t.w + c * (1 - t.w);
    } else {
        assert(false);
        return ym_zero3f;
    }
}

//
// Evaluates emission.
//
static inline ym_vec3f yt__eval_emission(const yt__point& pt) {
    if (ym_iszero(pt.ke)) return ym_zero3f;
    if (pt.shape) {
        switch (pt.shape->etype) {
            case yt_etype_point: return pt.ke;
            case yt_etype_line: return pt.ke;
            case yt_etype_triangle:
                if (ym_dot(pt.frame.norm, pt.wo) > 0)
                    return pt.ke;
                else
                    return ym_zero3f;
            default: assert(false); return ym_zero3f;
        }
    } else if (pt.env) {
        return pt.ke;
    } else {
        return ym_zero3f;
    }
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline float yt__eval_fresnel_dielectric(float cosw, float eta) {
    if (cosw < 0) {
        eta = 1 / eta;
        cosw = -cosw;
    }

    float sin2 = 1 - cosw * cosw;
    float eta2 = eta * eta;

    float cos2t = 1 - sin2 / eta2;
    if (cos2t < 0) return 1;  // tir

    float t0 = sqrt(cos2t);
    float t1 = eta * t0;
    float t2 = eta * cosw;

    float rs = (cosw - t1) / (cosw + t1);
    float rp = (t0 - t2) / (t0 + t2);

    return 0.5f * (rs * rs + rp * rp);
}

//
// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline float yt__eval_fresnel_metal(float cosw, float eta, float etak) {
    // TODO: vec3f
    if (!etak) return yt__eval_fresnel_dielectric(cosw, eta);

    cosw = ym_clamp(cosw, -1, 1);
    float cos2 = cosw * cosw;
    float sin2 = ym_clamp(1 - cos2, 0, 1);
    float eta2 = eta * eta;
    float etak2 = etak * etak;

    float t0 = eta2 - etak2 - sin2;
    float a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
    float t1 = a2plusb2 + cos2;
    float a = sqrt(0.5f * (a2plusb2 + t0));
    float t2 = 2 * a * cosw;
    float rs = (t1 - t2) / (t1 + t2);

    float t3 = cos2 * a2plusb2 + sin2 * sin2;
    float t4 = t2 * sin2;
    float rp = rs * (t3 - t4) / (t3 + t4);

    return 0.5f * (rp + rs);
}

//
// Evaluates the BRDF scaled by the cosine of the incoming direction.
//
// Implementation notes:
// - ggx from [Heitz 2014] and [Walter 2007]
// "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
// http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
// - uses Kajiya-Kay for hair
// - uses a hack for points
//
static inline ym_vec3f yt__eval_brdfcos(const yt__point& pt,
                                        const ym_vec3f& wi) {
    // summation over multiple terms
    ym_vec3f brdfcos = ym_zero3f;

    // exit if not needed
    if (ym_iszero(pt.kd) && ym_iszero(pt.ks)) return brdfcos;

    // save wo
    ym_vec3f wo = pt.wo;

    // compute wh
    ym_vec3f wh = ym_normalize(wo + wi);

    // compute dot products
    float ndo = ym_dot(pt.frame.norm, wo), ndi = ym_dot(pt.frame.norm, wi),
          ndh = ym_clamp(ym_dot(wh, pt.frame.norm), 0, 1);

    switch (pt.shape->etype) {
        case yt_etype_point: {
            // diffuse term (hack for now)
            if (!ym_iszero(pt.kd)) {
                float ido = ym_dot(wo, wi);
                ym_vec3f diff = pt.kd * (2 * ido + 1) / (2 * ym_pif);
                brdfcos += diff;
            }
        } break;
        case yt_etype_line: {
            // take sines
            float so = sqrtf(ym_clamp(1 - ndo * ndo, 0, 1)),
                  si = sqrtf(ym_clamp(1 - ndi * ndi, 0, 1)),
                  sh = sqrtf(ym_clamp(1 - ndh * ndh, 0, 1));

            // diffuse term (Kajiya-Kay)
            if (si > 0 && so > 0 && !ym_iszero(pt.kd)) {
                ym_vec3f diff = pt.kd * si / ym_pif;
                brdfcos += diff;
            }

            // specular term (Kajiya-Kay)
            if (si > 0 && so > 0 && sh > 0 && !ym_iszero(pt.ks)) {
                float ns = 2 / (pt.rs * pt.rs) - 2;
                float d = (ns + 2) * powf(sh, ns) / (2 + ym_pif);
                ym_vec3f spec = pt.ks * si * d / (4 * si * so);
                brdfcos += spec;
            }
        } break;
        case yt_etype_triangle: {
            // diffuse term
            if (ndi > 0 && ndo && !ym_iszero(pt.kd)) {
                ym_vec3f diff = pt.kd * ndi / ym_pif;
                brdfcos += diff;
            }

            // specular term (GGX)
            if (ndi > 0 && ndo > 0 && ndh > 0 && !ym_iszero(pt.ks)) {
                if (!pt.use_phong) {
                    // evaluate GGX
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = pt.rs * pt.rs;
                    float d = alpha2 / (ym_pif * cos2 * cos2 * (alpha2 + tan2) *
                                        (alpha2 + tan2));
                    float lambda_o =
                        (-1 + sqrtf(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i =
                        (-1 + sqrtf(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    ym_vec3f spec = pt.ks * ndi * d * g / (4 * ndi * ndo);
                    if (pt.es.x) {
                        // TODO: fixme
                        spec.x *=
                            yt__eval_fresnel_metal(ndh, pt.es.x, pt.eks.x);
                        spec.y *=
                            yt__eval_fresnel_metal(ndh, pt.es.y, pt.eks.y);
                        spec.z *=
                            yt__eval_fresnel_metal(ndh, pt.es.z, pt.eks.z);
                    }
                    brdfcos += spec;
                } else {
                    // evaluate Blinn-Phong
                    float ns = 2 / (pt.rs * pt.rs) - 2;
                    float d = (ns + 2) * powf(ndh, ns) / (2 + ym_pif);
                    ym_vec3f spec = pt.ks * ndi * d / (4 * ndi * ndo);
                    brdfcos += spec;
                }
            }
        } break;
        default: { assert(false); }
    }

    return brdfcos;
}

//
// Compute the weight for smapling the BRDF
//
static inline float yt__weight_brdfcos(const yt__point& pt,
                                       const ym_vec3f& wi) {
    // skip if no component
    if (ym_iszero(pt.kd) && ym_iszero(pt.ks)) return 0;

    // save wo
    ym_vec3f wo = pt.wo;

    // compute wh
    ym_vec3f wh = ym_normalize(wi + wo);

    // compute dot products
    float ndo = ym_dot(pt.frame.norm, wo), ndi = ym_dot(pt.frame.norm, wi),
          ndh = ym_dot(pt.frame.norm, wh);

    // check to make sure we are above the surface
    // updated this for refraction
    if (ndo <= 0 || ndi <= 0) return 0;

    // pick from a sum
    float wall = ym_mean(pt.kd) + ym_mean(pt.ks);
    float wd = ym_mean(pt.kd) / wall;
    float ws = ym_mean(pt.ks) / wall;

    // accumulate probability
    float pdf = 0;

    switch (pt.shape->etype) {
        case yt_etype_point:
        case yt_etype_line: {
            // diffuse term
            if (wall) {
                // homepherical cosine probability
                pdf += 1 / (4 * ym_pif);
            }
        } break;
        case yt_etype_triangle: {
            // diffuse term
            if (wd && ndi > 0) {
                // homepherical cosine probability
                pdf += wd * ndi / ym_pif;
            }

            // specular term (GGX or Phong)
            if (ws && ndi > 0 && ndo > 0 && ndh > 0) {
                if (!pt.use_phong) {
                    // probability proportional to d * ndh
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = pt.rs * pt.rs;
                    float d = alpha2 / (ym_pif * cos2 * cos2 * (alpha2 + tan2) *
                                        (alpha2 + tan2));
                    float hdo = ym_dot(wo, wh);
                    pdf += ws * d * ndh / (4 * hdo);
                } else {
                    // get phong exponent
                    float ns = 2 / (pt.rs * pt.rs) - 2;
                    // compute wh
                    ym_vec3f wh = ym_normalize(wi + wo);
                    float ndh = ym_dot(pt.frame.norm, wh);
                    // homerispherical cosine power probability
                    pdf += ws * powf(ndh, ns) * (ns + 1) / (2 * ym_pif);
                }
            }
        } break;
        default: { assert(false); }
    }

    // done
    return (pdf) ? 1 / pdf : 0;
}

//
// Picks a direction based on the BRDF
//
static inline ym_vec3f yt__sample_brdfcos(const yt__point& pt, float rnl,
                                          const ym_vec2f& rn) {
    // skip if no component
    if (ym_iszero(pt.kd) && ym_iszero(pt.ks)) return ym_zero3f;

    // save wo
    ym_vec3f wo = pt.wo;

    // compute cosine
    float ndo = ym_dot(pt.frame.norm, wo);

    // check to make sure we are above the surface
    // update this for refraction
    if (ndo <= 0) return ym_zero3f;

    // pick from a sum
    float wall = ym_mean(pt.kd) + ym_mean(pt.ks);
    float wd = ym_mean(pt.kd) / wall;
    float ws = ym_mean(pt.ks) / wall;

    switch (pt.shape->etype) {
        case yt_etype_point:
        case yt_etype_line: {
            if (wall > 0) {
                // sample wi with uniform spherical distribution
                float rz = rn.y, rr = sqrtf(1 - rz * rz),
                      rphi = 2 * ym_pif * rn.x;
                ym_vec3f wi_local = {rr * cosf(rphi), rr * sinf(rphi), rz};
                return ym_transform_direction(pt.frame, wi_local);
            }
        } break;
        case yt_etype_triangle: {
            // sample according to diffuse
            if (rnl < wd) {
                // sample wi with hemispherical cosine distribution
                float rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                      rphi = 2 * ym_pif * rn.x;
                // set to wi
                ym_vec3f wi_local = {rr * cosf(rphi), rr * sinf(rphi), rz};
                return ym_transform_direction(pt.frame, wi_local);
            }

            // sample according to specular (GGX or Phong)
            if (rnl >= wd && rnl < wd + ws) {
                if (!pt.use_phong) {
                    // sample wh with ggx distribution
                    float tan2 = pt.rs * pt.rs * rn.y / (1 - rn.y);
                    float rz = sqrtf(1 / (tan2 + 1)), rr = sqrtf(1 - rz * rz),
                          rphi = 2 * ym_pif * rn.x;
                    // set to wh
                    ym_vec3f wh_local = {rr * cosf(rphi), rr * sinf(rphi), rz};
                    ym_vec3f wh = ym_transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym_normalize(wh * 2 * ym_dot(wo, wh) - wo);
                } else {
                    // get phong exponent
                    float ns = 2 / (pt.rs * pt.rs) - 2;
                    // sample wh with hemispherical cosine power distribution
                    float rz = powf(rn.y, 1 / (ns + 1)),
                          rr = sqrtf(1 - rz * rz), rphi = 2 * ym_pif * rn.x;
                    // set to wh
                    ym_vec3f wh_local = {rr * cosf(rphi), rr * sinf(rphi), rz};
                    ym_vec3f wh = ym_transform_direction(pt.frame, wh_local);
                    // compute wi
                    return ym_normalize(wh * 2 * ym_dot(wo, wh) - wo);
                }
            }
        } break;
        default: { assert(false); }
    }

    assert(false);

    return ym_zero3f;
}

//
// Create a point for an environment map. Resolves material with textures.
//
static inline yt__point yt__eval_envpoint(const yt__env* env,
                                          const ym_vec3f& wo) {
    // set shape data
    yt__point pt;
    memset(&pt, 0, sizeof(pt));

    // check if null point
    if (!env) return pt;

    // shape params
    pt.env = env;

    // direction
    pt.wo = wo;

    // maerial
    pt.ke = env->ke;

    // textures
    if (env->ke_txt) {
        ym_vec3f w = ym_transform_direction(env->xform_inv, -wo);
        float theta = 1 - (acosf(ym_clamp(w.y, -1, 1)) / ym_pif);
        float phi = atan2f(w.z, w.x) / (2 * ym_pif);
        ym_vec2f texcoord = {phi, theta};
        pt.ke =
            pt.ke * yt__eval_texture(ym_vec3f{1, 1, 1}, env->ke_txt, texcoord);
    }

    return pt;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static inline yt__point yt__eval_shapepoint(const yt__shape* shape, int eid,
                                            const ym_vec2f& euv,
                                            const ym_vec3f& wo) {
    // set shape data
    yt__point pt;
    memset(&pt, 0, sizeof(pt));

    // check if null point
    if (!shape) return pt;

    // shape params
    pt.shape = shape;
    pt.eid = eid;
    pt.euv = euv;
    pt.etype = shape->etype;

    // direction
    pt.wo = wo;

    // compute weights
    float w[3];
    int idx[3], nidx = 0;
    switch (pt.etype) {
        case yt_etype_triangle: {
            const int* f = shape->elem + eid * 3;
            w[0] = 1 - euv.x - euv.y;
            w[1] = euv.x;
            w[2] = euv.y;
            idx[0] = f[0];
            idx[1] = f[1];
            idx[2] = f[2];
            nidx = 3;
        }; break;
        case yt_etype_line: {
            const int* f = shape->elem + eid * 2;
            w[0] = 1 - euv.x;
            w[1] = euv.x;
            idx[0] = f[0];
            idx[1] = f[1];
            nidx = 2;
        }; break;
        case yt_etype_point: {
            const int f = shape->elem[eid];
            w[0] = 1;
            idx[0] = f;
            nidx = 1;
        }; break;
    }

    // interpolate geometry positions
    pt.frame.pos = ym_zero3f;
    for (int i = 0; i < nidx; i++) pt.frame.pos += shape->pos[idx[i]] * w[i];

    pt.frame.norm = ym_zero3f;
    if (shape->norm) {
        for (int i = 0; i < nidx && shape->norm; i++)
            pt.frame.norm += shape->norm[idx[i]] * w[i];
        pt.frame.norm = ym_normalize(pt.frame.norm);
    }

    // creating frame
    pt.frame = ym_make_frame(pt.frame.pos, pt.frame.norm);

    // transform to world space
    // TODO: frame transform once we have them
    pt.frame.pos = ym_transform_point(shape->xform, pt.frame.pos);
    pt.frame.norm = ym_transform_direction(shape->xform, pt.frame.norm);
    pt.frame.tangu = ym_transform_direction(shape->xform, pt.frame.tangu);
    pt.frame.tangv = ym_transform_direction(shape->xform, pt.frame.tangv);

    // sample material data
    yt__material* mat = shape->mat;
    pt.ke = mat->ke;
    pt.kd = mat->kd;
    pt.ks = mat->ks;
    pt.rs = mat->rs;
    pt.use_phong = mat->use_phong;

    // handle surface color
    if (shape->color) {
        ym_vec3f color = ym_zero3f;
        for (int i = 0; i < nidx; i++) color += shape->color[idx[i]] * w[i];

        pt.ke *= color;
        pt.kd *= color;
        pt.ks *= color;
    }

    // handle textures
    if (shape->texcoord) {
        ym_vec2f texcoord = {0, 0};
        for (int i = 0; i < nidx; i++)
            texcoord += shape->texcoord[idx[i]] * w[i];

        if (mat->ke_txt) {
            pt.ke *= yt__eval_texture(ym_vec3f{1, 1, 1}, mat->ke_txt, texcoord);
        }
        if (mat->kd_txt) {
            pt.kd *= yt__eval_texture(ym_vec3f{1, 1, 1}, mat->kd_txt, texcoord);
        }
        if (mat->ks_txt) {
            pt.ks *= yt__eval_texture(ym_vec3f{1, 1, 1}, mat->ks_txt, texcoord);
        }
        if (mat->rs_txt) {
            pt.rs *=
                yt__eval_texture(ym_vec3f{1, 1, 1}, mat->rs_txt, texcoord).x;
        }
    }

    return pt;
}

//
// Sample weight for a light point.
//
static inline float yt__weight_light(const yt__light* light,
                                     const yt__point& lpt,
                                     const yt__point& pt) {
    if (light->shape) {
        switch (light->shape->etype) {
            case yt_etype_point: {
                float d = ym_dist(lpt.frame.pos, pt.frame.pos);
                return light->area / (d * d);
            } break;
            case yt_etype_line: {
                assert(false);
                return 0;
            } break;
            case yt_etype_triangle: {
                float d = ym_dist(lpt.frame.pos, pt.frame.pos);
                return light->area * fabsf(ym_dot(lpt.frame.norm, lpt.wo)) /
                       (d * d);
            } break;
            default: { assert(false); } break;
        }
    } else if (light->env) {
        return 4 * ym_pif;
    } else {
        assert(false);
    }
    return 0;
}

//
// Picks a point on a light.
//
static inline yt__point yt__sample_light(const yt__light* light,
                                         const yt__point& pt, float rne,
                                         const ym_vec2f& rn) {
    if (light->shape) {
        int eid;
        for (eid = 0; eid < light->shape->nelems && light->cdf[eid] < rne;
             eid += 1) {
        }
        if (eid > light->shape->nelems - 1) eid = light->shape->nelems - 1;
        ym_vec2f euv = rn;
        if (light->shape->etype == yt_etype_triangle) {
            euv.x = 1 - sqrt(rn.x);
            euv.y = rn.y * sqrt(rn.x);
        }
        yt__point lpt = yt__eval_shapepoint(light->shape, eid, euv, ym_zero3f);
        lpt.wo = ym_normalize(pt.frame.pos - lpt.frame.pos);
        return lpt;
    } else if (light->env) {
        float z = -1 + 2 * rn.y;
        float rr = sqrtf(ym_clamp(1 - z * z, 0, 1));
        float phi = 2 * ym_pif * rn.x;
        ym_vec3f wo = {cosf(phi) * rr, z, sinf(phi) * rr};
        yt__point lpt = yt__eval_envpoint(light->env, wo);
        return lpt;
    } else {
        assert(false);
    }
    return yt__point();
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym_ray3f yt__offset_ray(const yt_scene* scene,
                                      const yt__point& pt, const ym_vec3f& w,
                                      const yt_render_params& params) {
    return ym_ray3f(pt.frame.pos + pt.frame.norm * params.ray_eps, w,
                    params.ray_eps);
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline ym_ray3f yt__offset_ray(const yt_scene* scene,
                                      const yt__point& pt, const yt__point& pt2,
                                      const yt_render_params& params) {
    float ray_dist =
        (pt2.shape) ? ym_dist(pt.frame.pos, pt2.frame.pos) : FLT_MAX;
    return ym_ray3f(pt.frame.pos + pt.frame.norm * params.ray_eps, -pt2.wo,
                    params.ray_eps, ray_dist - 2 * params.ray_eps);
}

//
// Intersects a ray with the scene and return the point (or env point).
//
static inline yt__point yt__intersect_scene(const yt_scene* scene,
                                            const ym_vec3f& ray_o,
                                            const ym_vec3f& ray_d,
                                            const yt_render_params& params) {
    int sid, eid;
    ym_vec2f euv;
    float ray_t;
    ym_ray3f ray = ym_ray3f(ray_o, ray_d, params.ray_eps, HUGE_VALF);
    bool hit =
        scene->intersect_ray(scene->ray_ctx, ray, &ray_t, &sid, &eid, &euv);
    ym_vec3f wo = -ray_d;
    if (hit) {
        return yt__eval_shapepoint(&scene->shapes[sid], eid, euv, wo);
    } else {
        return yt__eval_envpoint(&scene->env[0], wo);
    }
}

//
// Intersects a ray with the scene and return the point (or env point).
//
static inline yt__point yt__intersect_scene(const yt_scene* scene,
                                            const ym_ray3f& ray) {
    int sid, eid;
    ym_vec2f euv;
    float ray_t;
    bool hit =
        scene->intersect_ray(scene->ray_ctx, ray, &ray_t, &sid, &eid, &euv);
    ym_vec3f wo = -ray.d;
    if (hit) {
        return yt__eval_shapepoint(&scene->shapes[sid], eid, euv, wo);
    } else {
        return yt__eval_envpoint(&scene->env[0], wo);
    }
}

//
// Evalutes direct illumination using MIS.
//
static inline ym_vec3f yt__eval_direct(const yt_scene* scene,
                                       const yt__light* light,
                                       const yt__point& pt,
                                       yt__sampler* sampler,
                                       const yt_render_params& params) {
    // select whether it goes in all light mode
    bool all_lights = !light;

    // pick a light if not there
    float nlweight = 0;
    if (all_lights) {
        int lid = yt__sample_next1f(sampler) * scene->lights.size();
        if (lid > scene->lights.size() - 1) lid = (int)scene->lights.size() - 1;
        light = &scene->lights[lid];
        nlweight = scene->lights.size();
    } else {
        nlweight = 1;
    }

    // sample light according to area
    yt__point lpt = yt__sample_light(light, pt, yt__sample_next1f(sampler),
                                     yt__sample_next2f(sampler));
    ym_vec3f lld = yt__eval_emission(lpt) * yt__eval_brdfcos(pt, -lpt.wo);
    float lweight = yt__weight_light(light, lpt, pt) * nlweight;
    lld *= lweight;
    if (!ym_iszero(lld)) {
        ym_ray3f shadow_ray = yt__offset_ray(scene, pt, lpt, params);
        if (scene->hit_ray(scene->ray_ctx, shadow_ray)) lld = ym_zero3f;
    }

    // check if mis is necessary
    if (pt.shape &&
        (pt.shape->etype == yt_etype_point || pt.shape->etype == yt_etype_line))
        return lld;

    // check if mis is necessary
    if (light->shape && (light->shape->etype == yt_etype_point ||
                         light->shape->etype == yt_etype_line))
        return lld;

    // sample the brdf
    ym_vec3f bwi = yt__sample_brdfcos(pt, yt__sample_next1f(sampler),
                                      yt__sample_next2f(sampler));
    float bweight = 0;
    ym_vec3f bld = ym_zero3f;
    yt__point bpt =
        yt__intersect_scene(scene, yt__offset_ray(scene, pt, bwi, params));
    if (light->shape == bpt.shape || all_lights) {
        bweight = yt__weight_brdfcos(pt, bwi);
        bld = yt__eval_emission(bpt) * yt__eval_brdfcos(pt, bwi) * bweight;
    }

    // accumulate the value with mis
    if (!ym_iszero(lld)) {
        float bweight = yt__weight_brdfcos(pt, -lpt.wo);
        // float weight =
        //     (1 / lweight) * (1 / lweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        float weight = (1 / lweight) / ((1 / lweight) + (1 / bweight));
        lld *= weight;
    }
    if (!ym_iszero(bld)) {
        float lweight = yt__weight_light(light, bpt, pt) * nlweight;
        // float weight =
        //     (1 / bweight) * (1 / bweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        float weight = (1 / bweight) / ((1 / lweight) + (1 / bweight));
        bld *= weight;
    }

    // return weighted sum
    return lld + bld;
}

//
// Recursive path tracing.
//
static inline ym_vec3f yt__shade_pathtrace_recd(
    const yt_scene* scene, const ym_ray3f& ray, yt__sampler* sampler, bool* hit,
    int ray_depth, const yt_render_params& params) {
    yt__point pt = yt__intersect_scene(scene, ray);
    if (hit) *hit = pt.shape;

    ym_vec3f l = (!ray_depth) ? yt__eval_emission(pt) : ym_zero3f;
    if (!pt.shape) return l;

    if (ym_iszero(pt.kd) && ym_iszero(pt.ks)) return l;

    ym_vec3f ld = yt__eval_direct(scene, 0, pt, sampler, params);
    l += ld;

    if (ray_depth >= params.max_depth) return l;

    // roussian roulette
    float rrweight = 1;
    if (ray_depth >= params.min_depth) {
        float rrrn = yt__sample_next1f(sampler);
        if (rrrn >= fmin(pt.wrr, 0.95f)) return l;
        rrweight = 1.0f / fmin(pt.wrr, 0.95f);
    }

    // continue path
    ym_vec3f bwi = yt__sample_brdfcos(pt, yt__sample_next1f(sampler),
                                      yt__sample_next2f(sampler));
    if (ym_iszero(bwi)) return l;
    float bweight = yt__weight_brdfcos(pt, bwi);
    if (!bweight) return l;
    ym_vec3f bbrdfcos = yt__eval_brdfcos(pt, bwi);
    if (ym_iszero(bbrdfcos)) return l;
    ym_vec3f ble =
        yt__shade_pathtrace_recd(scene, yt__offset_ray(scene, pt, bwi, params),
                                 sampler, 0, ray_depth + 1, params);
    l += ble * bbrdfcos * rrweight;

    return l;
}

//
// Shader interface for the above function.
//
static inline ym_vec3f yt__shade_pathtrace(const yt_scene* scene,
                                           const ym_ray3f& ray,
                                           yt__sampler* sampler, bool* hit,
                                           const yt_render_params& params) {
    return yt__shade_pathtrace_recd(scene, ray, sampler, hit, 0, params);
}

//
// Direct illuination.
//
static inline ym_vec3f yt__shade_direct(const yt_scene* scene,
                                        const ym_ray3f& ray,
                                        yt__sampler* sampler, bool* hit,
                                        const yt_render_params& params) {
    yt__point pt = yt__intersect_scene(scene, ray);
    if (hit) *hit = pt.shape;

    ym_vec3f l = yt__eval_emission(pt);
    if (!pt.shape) return l;

    if (ym_iszero(pt.kd) && ym_iszero(pt.ks)) return l;

    if (!ym_iszero(params.amb)) {
        l += params.amb * pt.kd;
    }

    for (int lid = 0; lid < scene->lights.size(); lid++) {
        const yt__light* light = &scene->lights[lid];
        l += yt__eval_direct(scene, light, pt, sampler, params);
    }

    return l;
}

//
// Eyelight for quick previewing.
//
static inline ym_vec3f yt__shade_eyelight(const yt_scene* scene,
                                          const ym_ray3f& ray,
                                          yt__sampler* sampler, bool* hit,
                                          const yt_render_params& params) {
    yt__point pt = yt__intersect_scene(scene, ray);
    if (hit) *hit = pt.shape;

    ym_vec3f l = yt__eval_emission(pt);
    if (!pt.shape) return l;

    ym_vec3f brdfcos = yt__eval_brdfcos(pt, pt.wo) * ym_pif;
    l += brdfcos;

    return l;
}

//
// Shader function callback.
//
typedef ym_vec3f (*shade_fn)(const yt_scene* scene, const ym_ray3f& ray,
                             yt__sampler* sampler, bool* hit,
                             const yt_render_params& params);

//
// Renders a block of pixels. Public API, see above.
//
YGL_API void yt_trace_block(const yt_scene* scene, int cid,
                            ym_vec4f* img_pixels, int img_w, int img_h, int ns,
                            const ym_range2i& block, const ym_range1i& samples,
                            const yt_render_params& params) {
    const yt__camera* cam = &scene->camera[cid];
    shade_fn shade;
    switch (params.stype) {
        case yt_stype_eyelight: shade = yt__shade_eyelight; break;
        case yt_stype_default:
        case yt_stype_direct: shade = yt__shade_direct; break;
        case yt_stype_pathtrace: shade = yt__shade_pathtrace; break;
        default: assert(false); return;
    }
    for (int j = block.min.y; j < block.max.y; j++) {
        for (int i = block.min.x; i < block.max.x; i++) {
            ym_vec4f* pixel = img_pixels + (img_w * j + i);
            *pixel = ym_zero4f;
            for (int s = samples.min; s < samples.max; s++) {
                yt__sampler sampler =
                    yt__make_sampler(i, j, s, ns, params.rtype);
                ym_vec2f rn = yt__sample_next2f(&sampler);
                ym_vec2f uv = {(i + rn.x) / img_w, 1 - (j + rn.y) / img_h};
                bool hit;
                ym_ray3f ray =
                    yt__eval_camera(cam, uv, yt__sample_next2f(&sampler));
                ym_vec3f l = shade(scene, ray, &sampler, &hit, params);
                if (!ym_isfinite(l)) continue;
                if (params.pixel_clamp > 0) ym_clamplen(l, params.pixel_clamp);
                *pixel += {l, 1};
            }
            *pixel /= (float)samples.size();
        }
    }
}

//
// Renders the whole image. Public API, see above.
//
YGL_API void yt_trace_image(const yt_scene* scene, int cid,
                            ym_vec4f* img_pixels, int img_w, int img_h, int ns,
                            const yt_render_params& params) {
    yt_trace_block(scene, cid, img_pixels, img_w, img_h, ns,
                   {{0, 0}, {img_w, img_h}}, {0, ns}, params);
}

// -----------------------------------------------------------------------------
// C API IMPLEMENTATION
// -----------------------------------------------------------------------------

//
// Initialize the scene data with nshapes, nmaterials, ntextures.
//
YGLC_API yt_scene* ytc_make_scene(int ncameras, int nshapes, int nmaterials,
                                  int ntextures, int nenvs) {
    return yt_make_scene(ncameras, nshapes, nmaterials, ntextures, nenvs);
}

//
// Cleanup scene data
//
YGLC_API void ytc_free_scene(yt_scene* scene) { return yt_free_scene(scene); }

//
// Set intersection callback. See callback description above.
//
YGLC_API void ytc_set_intersection(yt_scene* scene, void* ray_ctx,
                                   yt_intersect_ray intersect_ray,
                                   yt_hit_ray hit_ray);

//
// Set rendering camera.
//
YGLC_API void ytc_set_camera(yt_scene* scene, int cid, float xform[16],
                             float width, float height, float aperture,
                             float focus) {
    return yt_set_camera(scene, cid, ym_frame3f((ym_mat4f)xform), width, height,
                         aperture, focus);
}

//
// Set rendering environment map.
//
YGL_API void ytc_set_env(yt_scene* scene, int eid, float xform[16], float ke[3],
                         int ke_txt) {
    return yt_set_env(scene, eid, ym_frame3f((ym_mat4f)xform), (ym_vec3f)ke,
                      ke_txt);
}

//
// Set rendering camera.
//
YGL_API void ytc_init_lights(yt_scene* scene) { return yt_init_lights(scene); }

//
// Set a shape.
//
YGL_API void ytc_set_shape(yt_scene* scene, int sid, int mid, float xform[16],
                           int nelems, int* elem, int etype, int nverts,
                           float* pos, float* norm, float* texcoord,
                           float* color, float* radius) {
    return yt_set_shape(scene, sid, mid, ym_frame3f(*(ym_mat4f*)xform), nelems,
                        elem, etype, nverts, (ym_vec3f*)pos, (ym_vec3f*)norm,
                        (ym_vec2f*)texcoord, (ym_vec3f*)color, radius);
}

//
// Set a material.
//
YGL_API void ytc_set_material(yt_scene* scene, int mid, float ke[3],
                              float kd[3], float ks[3], float rs, float es[3],
                              float esk[3], int ke_txt, int kd_txt, int ks_txt,
                              int rs_txt, bool use_phong) {
    return yt_set_material(scene, mid, *(ym_vec3f*)ke, *(ym_vec3f*)kd,
                           *(ym_vec3f*)ks, rs, *(ym_vec3f*)es, *(ym_vec3f*)esk,
                           ke_txt, kd_txt, ks_txt, rs_txt, use_phong);
}

//
// Convert a Phong exponent to GGX/Phong roughness
//
YGL_API float ytc_specular_exponent_to_roughness(float n) {
    return yt_specular_exponent_to_roughness(n);
}

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YGL_API void ytc_specular_fresnel_from_ks(const float ks[3], float es[3],
                                          float esk[3]) {
    return yt_specular_fresnel_from_ks(*(ym_vec3f*)ks, (ym_vec3f*)es,
                                       (ym_vec3f*)esk);
}

//
// Set a texture.
//
YGL_API void ytc_set_texture(yt_scene* scene, int tid, float* pixels, int w,
                             int h, int nc) {
    return yt_set_texture(scene, tid, pixels, w, h, nc);
}

//
// Set rendering parameters
//
YGL_API ytc_render_params ytc_make_rendering_params() {
    yt_render_params cpp_params;
    ytc_render_params params;
    memcpy(&params, &cpp_params, sizeof(params));
    return params;
}

//
// Renders a block of sample
//
YGL_API void ytc_trace_block(const yt_scene* scene, int cid, float* pixels,
                             int w, int h, int ns, const int block[4],
                             const int samples[2],
                             const ytc_render_params* params) {
    return yt_trace_block(scene, cid, (ym_vec4f*)pixels, w, h, ns,
                          {{block[0], block[1]}, {block[2], block[3]}},
                          {samples[0], samples[1]},
                          *(const yt_render_params*)params);
}

//
// Convenience function to call yt_trace_block with all sample at once.
//
YGL_API void ytc_trace_image(const yt_scene* scene, int cid, float* pixels,
                             int w, int h, int ns,
                             const ytc_render_params* params) {
    return yt_trace_image(scene, cid, (ym_vec4f*)pixels, w, h, ns,
                          *(const yt_render_params*)params);
}

#endif

#endif
