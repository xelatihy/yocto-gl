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
// the primitive type (points, lines, triangles, quads),
// an array of vertex positions, and an array of vertex radius
// (for points and lines).
//
// Quad meshes are experimental and might go away in future realeases. If
// you can, please use triangles. Quads are treated as two triangles (v0,v1,v3)
// and (v2,v3,v1). Quads with v2 == v3 are degenerate and represent one
// triangle, thus quads meshes can also represent mixtures of triangle
// and quads. This follows Intel's Embree API. But please note that irregular
// quads will be sampled very very poorly. Really you are better off with
// triangles.
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
// All functions in this library are inlined by default for ease of use.
// To use the library as a .h/.c pair do the following:
// - to use as a .h, just #define YT_NOINLINE before including this file
// - to use as a .c, just #define YT_IMPLEMENTATION before including this file
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

//
//  LICENSE of included software
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//

#ifndef _YT_H_
#define _YT_H_

#ifndef YT_NOINLINE
#define YT_API static inline
#else
#ifdef __cplusplus
#define YT_API extern "C"
#else
#define YT_API
#endif
#endif

#include <stdbool.h>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Shape element types
//
enum {
    yt_etype_point = 1,     // points
    yt_etype_line = 2,      // lines
    yt_etype_triangle = 3,  // triangle
    yt_etype_quad = 4,      // quads
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
typedef bool (*yt_intersect_ray)(void* ctx, const float ray_o[3],
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
typedef bool (*yt_hit_ray)(void* ctx, const float ray_o[3],
                           const float ray_d[3], float ray_min, float ray_max,
                           int ray_mask);

//
// Initialize the scene data with nshapes, nmaterials, ntextures.
//
YT_API yt_scene*
yt_make_scene(int nshapes, int nmaterials, int ntextures);

//
// Cleanup scene data
//
YT_API void
yt_free_scene(yt_scene* scene);

//
// Set intersection callback. See callback description above.
//
YT_API void
yt_set_intersection(yt_scene* scene, void* ray_ctx,
                    yt_intersect_ray intersect_ray, yt_hit_ray hit_ray);

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
YT_API void
yt_set_camera(yt_scene* scene, float xform[16], float width, float height,
              float aperture, float focus);

//
// Set rendering environment map.
//
// Parameters:
// - scene: trace scene
// - xform: env local to world **rigid** transform
// - ke: emission
// - ke_txt: ke texture index (-1 for no texture)
//
YT_API void
yt_set_env(yt_scene* scene, float xform[16], float ke[3], int ke_txt);

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
YT_API void
yt_init_lights(yt_scene* scene);

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
YT_API void
yt_set_shape(yt_scene* scene, int sid, int mid, float xform[16], int nelems,
             int* elem, int etype, int nverts, float* pos, float* norm,
             float* texcoord, float* color, float* radius);

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
YT_API void
yt_set_material(yt_scene* scene, int mid, float ke[3], float kd[3], float ks[3],
                float rs, float es[6], int ke_txt, int kd_txt, int ks_txt,
                int rs_txt, bool use_phong);

//
// Convert a Phong exponent to GGX/Phong roughness
//
YT_API float
yt_specular_exponent_to_roughness(float n);

//
// Estimates the fresnel coefficient es from ks at normal incidence
//
YT_API void
yt_specular_fresnel_from_ks(const float ks[3], float es[6]);

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
YT_API void
yt_set_texture(yt_scene* scene, int tid, float* pixels, int w, int h, int nc);

//
// Set rendering parameters
//
// Parameters:
// - scene: trace scene
// - stype: rendering algorithm type
// - rtype: random number generator type
// - amb: ambient lighting for direct rendering algorithm
//
YT_API void
yt_set_rendering_params(yt_scene* scene, int stype, int rtype, float amb[3]);

//
// Set advanced rendering parameters (normally not used)
//
// Parameters:
// - scene: trace scene
// - max_depth: max ray depth in path tracer
// - pixel_clamp: maximum sample value (used to discard very large values)
// - ray_eps: ray epsilon to avoid self-intersection
//
YT_API void
yt_set_adv_rendering_params(yt_scene* scene, int max_depth, float pixel_clamp,
                            float ray_eps);

//
// Renders a block of sample
//
// Parameters:
// - scene: trace scene
// - pixels: pixel data in RGBA format
// - w, h: image width and height
// - ns: number of samples
// - block: sample block, respecively [xmin, xmax, ymin, ymax, sample_min,
// sample_max]; max values are excluded
//
// Notes: It is safe to call the function in parallel one different blocks.
// But two threads should not access the same pixels at the same time.
// Also blocks with different samples should be called sequentially.
//
YT_API void
yt_trace_block(const yt_scene* scene, float* pixels, int w, int h, int ns,
               int block[6]);

//
// Convenience function to call yt_trace_block with all sample at once.
//
YT_API void
yt_trace_image(const yt_scene* scene, float* pixels, int w, int h, int ns);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YT_NOINLINE) || defined(YT_IMPLEMENTATION)

#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// MATH SUPPORT
// -----------------------------------------------------------------------------

#define yt__pif 3.14159265f

// clamps float values
static inline float
yt__fclampf(float x, float m, float M) {
    return fminf(M, fmaxf(m, x));
}

// 2d vectors
typedef struct { float x, y; } yt__vec2f;

// 3d vectors
typedef struct { float x, y, z; } yt__vec3f;

// 4d vectors
typedef struct { float x, y, z, w; } yt__vec4f;

// 4x4 rigid transform matrix
typedef struct { float m[16]; } yt__mat4f;

// check if zero
static inline bool
yt__iszero3f(const yt__vec3f v) {
    return v.x == 0 && v.y == 0 && v.z == 0;
}

// check if zero
static inline bool
yt__iszero3fv(const float v[3]) {
    return v[0] == 0 && v[1] == 0 && v[2] == 0;
}

// check if finite
static inline bool
yt__isfinite3fv(const float v[3]) {
    return isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
}

// check if finite
static inline bool
yt__isfinite3f(const yt__vec3f v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

// zero vector
static inline yt__vec3f
yt__zero3f() {
    return (yt__vec3f){ 0, 0, 0 };
}

// vector negate
static inline yt__vec3f
yt__neg3f(const yt__vec3f a) {
    return (yt__vec3f){ -a.x, -a.y, -a.z };
}

// vector sum
static inline yt__vec3f
yt__sum3f(const yt__vec3f a, const yt__vec3f b) {
    return (yt__vec3f){ a.x + b.x, a.y + b.y, a.z + b.z };
}

// vector subtraction
static inline yt__vec3f
yt__sub_vec3f(const yt__vec3f a, const yt__vec3f b) {
    return (yt__vec3f){ a.x - b.x, a.y - b.y, a.z - b.z };
}

// vector multiplication
static inline yt__vec3f
yt__mul3f(const yt__vec3f a, const yt__vec3f b) {
    return (yt__vec3f){ a.x * b.x, a.y * b.y, a.z * b.z };
}

// vector scalar multiply
static inline yt__vec3f
yt__smul3f(const yt__vec3f a, float b) {
    return (yt__vec3f){ a.x * b, a.y * b, a.z * b };
}

// vector dot product
static inline float
yt__dot3f(const yt__vec3f a, const yt__vec3f b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

// vector length
static inline float
yt__length3f(const yt__vec3f a) {
    return (float)sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

// vector normalization
static inline yt__vec3f
yt__normalize3f(const yt__vec3f a) {
    return yt__smul3f(a, (1 / yt__length3f(a)));
}

// point distance
static inline float
yt__dist3f(const yt__vec3f a, const yt__vec3f b) {
    return (float)sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                       (a.z - b.z) * (a.z - b.z));
}

// vector cross product
static inline yt__vec3f
yt__cross3f(const yt__vec3f a, const yt__vec3f b) {
    return (yt__vec3f){ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x };
}

// orthogonal vector
static inline yt__vec3f
yt__orthogonal3f(const yt__vec3f v) {
    return fabs(v.x) > fabs(v.z) ? (yt__vec3f){ -v.y, v.x, 0 }
                                 : (yt__vec3f){ 0, -v.z, v.y };
}

// vector mean
static inline float
yt__mean3f(const yt__vec3f v) {
    return (v.x + v.y + v.z) / 3;
}

// vector clamp
static inline yt__vec3f
yt__clamp3f(const yt__vec3f v, float m) {
    float l = yt__length3f(v);
    if (l > m)
        return yt__smul3f(v, m / l);
    else
        return v;
}

// vector sum
static inline yt__vec4f
yt__sum4f(const yt__vec4f a, const yt__vec4f b) {
    return (yt__vec4f){ a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

// vector scale
static inline yt__vec4f
yt__smul4f(const yt__vec4f a, float b) {
    return (yt__vec4f){ a.x * b, a.y * b, a.z * b, a.w * b };
}

// identity matrix
static inline yt__mat4f
yt__identity4f() {
    return (yt__mat4f){ { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 } };
}

// matric inverse
// code from http://stackoverflow.com/questions/2624422/
// efficient-4x4-matrix-inverse-affine-transform
static inline yt__mat4f
yt__inverse4f(const yt__mat4f m) {
    float a[4][4];
    *(yt__mat4f*)a = m;

    float s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    float s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
    float s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3];
    float s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
    float s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3];
    float s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];

    float c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3];
    float c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
    float c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2];
    float c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
    float c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2];
    float c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];

    // TODO: Should check for 0 determinant
    float invdet =
        1.0f / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

    float b[4][4];

    b[0][0] = (a[1][1] * c5 - a[1][2] * c4 + a[1][3] * c3) * invdet;
    b[0][1] = (-a[0][1] * c5 + a[0][2] * c4 - a[0][3] * c3) * invdet;
    b[0][2] = (a[3][1] * s5 - a[3][2] * s4 + a[3][3] * s3) * invdet;
    b[0][3] = (-a[2][1] * s5 + a[2][2] * s4 - a[2][3] * s3) * invdet;

    b[1][0] = (-a[1][0] * c5 + a[1][2] * c2 - a[1][3] * c1) * invdet;
    b[1][1] = (a[0][0] * c5 - a[0][2] * c2 + a[0][3] * c1) * invdet;
    b[1][2] = (-a[3][0] * s5 + a[3][2] * s2 - a[3][3] * s1) * invdet;
    b[1][3] = (a[2][0] * s5 - a[2][2] * s2 + a[2][3] * s1) * invdet;

    b[2][0] = (a[1][0] * c4 - a[1][1] * c2 + a[1][3] * c0) * invdet;
    b[2][1] = (-a[0][0] * c4 + a[0][1] * c2 - a[0][3] * c0) * invdet;
    b[2][2] = (a[3][0] * s4 - a[3][1] * s2 + a[3][3] * s0) * invdet;
    b[2][3] = (-a[2][0] * s4 + a[2][1] * s2 - a[2][3] * s0) * invdet;

    b[3][0] = (-a[1][0] * c3 + a[1][1] * c1 - a[1][2] * c0) * invdet;
    b[3][1] = (a[0][0] * c3 - a[0][1] * c1 + a[0][2] * c0) * invdet;
    b[3][2] = (-a[3][0] * s3 + a[3][1] * s1 - a[3][2] * s0) * invdet;
    b[3][3] = (a[2][0] * s3 - a[2][1] * s1 + a[2][2] * s0) * invdet;

    return *(yt__mat4f*)b;
}

// trasform point (rigid transform)
static inline yt__vec3f
yt__transform_point3f(const yt__mat4f a, const yt__vec3f b) {
    return (yt__vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z + a.m[3],
                        a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z + a.m[7],
                        a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z + a.m[11] };
}

// trasform direction (rigid transform)
static inline yt__vec3f
yt__transform_direction3f(const yt__mat4f a, const yt__vec3f b) {
    return yt__normalize3f(
        (yt__vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z,
                     a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z,
                     a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z });
}

//
// Camera
//
typedef struct yt__camera {
    yt__mat4f xform;        // local-to-world transform
    float width, height;    // image width and height
    float aperture, focus;  // lens aperture and focus plane distance
} yt__camera;

//
// Texture
//
typedef struct yt__texture {
    int w, h, nc;   // width, height, number of components
    float* pixels;  // pixel data
} yt__texture;

//
// Material
//
typedef struct yt__material {
    yt__vec3f ke, kd, ks;  // emission, diffuse, specular term
    float rs;              // specular roughness API
    yt__vec3f es, eks;     // eta
    yt__texture *ke_txt, *kd_txt, *ks_txt, *rs_txt;  // textures
    bool use_phong;                                  // whether to use phong
} yt__material;

//
// Shape
//
typedef struct yt__shape {
    yt__mat4f xform;                   // local-to-world rigid transform
    bool xformed;                      // whether it is transformed
    yt__material* mat;                 // material
    int nelems, *elem, etype, nverts;  // element data
    float *pos, *norm, *texcoord, *color, *radius;  // vertex data
} yt__shape;

//
// Environment
//
typedef struct yt__env {
    yt__mat4f xform;      // local-to-world rigid transform
    yt__mat4f xform_inv;  // world-to-local rigid transform
    bool xformed;         // transformed
    yt__vec3f ke;         // emission
    yt__texture* ke_txt;  // emission texture
} yt__env;

//
// Light (either shape or environment)
//
typedef struct yt__light {
    yt__shape* shape;  // shape
    yt__env* env;      // environment
    float* cdf;        // for shape, cdf of shape elements for sampling
    float area;        // for shape, shape area
} yt__light;

//
// Scene
//
struct yt_scene {
    void* ray_ctx;                   // ray context
    yt_intersect_ray intersect_ray;  // ray intersection callback
    yt_hit_ray hit_ray;              // ray hit callback

    yt__camera* camera;                           // camera
    yt__env* env;                                 // env
    int nshapes, nmaterials, ntextures, nlights;  // number of elements
    yt__shape* shapes;                            // shapes
    yt__material* materials;                      // materials
    yt__texture* textures;                        // textures
    yt__light* lights;                            // lights

    int stype;      // rendering algorithm type
    int rtype;      // random number generator type
    yt__vec3f amb;  // ambient illumination for direct

    int min_ray_depth, max_ray_depth;  // min/max ray depth
    float pixel_clamp;                 // pixel clamping

    float ray_eps;  // ray interseciton epsilon
};

//
// Init scene. Public API, see above.
//
YT_API yt_scene*
yt_make_scene(int nshapes, int nmaterials, int ntextures) {
    yt_scene* scene = (yt_scene*)calloc(1, sizeof(yt_scene));
    scene->camera = (yt__camera*)calloc(1, sizeof(yt__camera));
    scene->camera->xform = yt__identity4f();
    scene->camera->width = 1;
    scene->camera->height = 1;
    scene->camera->aperture = 0;
    scene->nshapes = nshapes;
    scene->shapes = (yt__shape*)calloc(nshapes, sizeof(yt__shape));
    scene->nmaterials = nmaterials;
    scene->materials = (yt__material*)calloc(nmaterials, sizeof(yt__material));
    scene->ntextures = ntextures;
    scene->textures = (yt__texture*)calloc(ntextures, sizeof(yt__texture));
    scene->stype = 0;
    scene->rtype = 0;
    scene->amb = yt__zero3f();
    scene->ray_eps = 1e-2f;
    scene->min_ray_depth = 3;
    scene->max_ray_depth = 8;
    scene->pixel_clamp = 100;
    return scene;
}

//
// Free scene. Public API, see above.
//
YT_API void
yt_free_scene(yt_scene* scene) {
    if (scene->camera) free(scene->camera);
    if (scene->shapes) free(scene->shapes);
    if (scene->materials) free(scene->materials);
    if (scene->textures) free(scene->textures);
    if (scene->lights) free(scene->lights);
    if (scene->env) free(scene->env);
}

//
// Set rendering params. Public API, see above.
//
YT_API void
yt_set_rendering_params(yt_scene* scene, int stype, int rtype, float amb[3]) {
    scene->stype = stype;
    scene->rtype = rtype;
    if (amb) {
        scene->amb = (yt__vec3f){ amb[0], amb[1], amb[2] };
    } else {
        scene->amb = (yt__vec3f){ 0, 0, 0 };
    }
}

//
// Set advanced rendering params. Public API, see above.
//
YT_API void
yt_set_adv_rendering_params(yt_scene* scene, int ray_depth, float pixel_clamp,
                            float ray_eps) {
    scene->max_ray_depth = ray_depth ? ray_depth : 8;
    scene->min_ray_depth =
        (scene->min_ray_depth > 3) ? 3 : scene->min_ray_depth;
    scene->pixel_clamp = pixel_clamp;
    scene->ray_eps = ray_eps;
}

//
// Set intersection callbacks. Public API, see above.
//
YT_API void
yt_set_intersection(yt_scene* scene, void* ray_ctx,
                    yt_intersect_ray intersect_ray, yt_hit_ray hit_ray) {
    scene->ray_ctx = ray_ctx;
    scene->intersect_ray = intersect_ray;
    scene->hit_ray = hit_ray;
}

//
// Set camera. Public API, see above.
//
YT_API void
yt_set_camera(yt_scene* scene, float xform[16], float width, float height,
              float aperture, float focus) {
    yt__camera* cam = scene->camera;
    cam->xform = *(yt__mat4f*)xform;
    cam->width = width;
    cam->height = height;
    cam->aperture = aperture;
    cam->focus = (focus) ? focus : 1;
}

//
// Set shape. Public API, see above.
//
YT_API void
yt_set_shape(yt_scene* scene, int sid, int matid, float xform[16], int nelems,
             int* elem, int etype, int nverts, float* pos, float* norm,
             float* texcoord, float* color, float* radius) {
    assert(matid >= 0 && matid < scene->nmaterials);
    yt__shape* shape = scene->shapes + sid;
    shape->xformed = (bool)xform;
    shape->xform = (xform) ? (*(yt__mat4f*)xform) : yt__identity4f();
    shape->mat = (matid < 0) ? 0 : scene->materials + matid;
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
static inline void
yt__compute_shape_cdf(int nelems, int* elem, int etype, yt__vec3f* pos,
                      float* cdf, float* area) {
    switch (etype) {
        case yt_etype_point: {
            for (int i = 0; i < nelems; i++) cdf[i] = i + 1;
        } break;
        case yt_etype_line: {
            for (int i = 0; i < nelems; i++) {
                int* l = elem + i * 2;
                cdf[i] = yt__dist3f(pos[l[0]], pos[l[1]]);
                if (i != 0) cdf[i] += cdf[i - 1];
            }
        } break;
        case yt_etype_triangle: {
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i * 3;
                yt__vec3f *v0 = pos + f[0], *v1 = pos + f[1], *v2 = pos + f[2];
                yt__vec3f e1 = yt__sub_vec3f(*v1, *v0);
                yt__vec3f e2 = yt__sub_vec3f(*v2, *v0);
                yt__vec3f a = yt__cross3f(e1, e2);
                cdf[i] = yt__length3f(a) / 2;
                if (i != 0) cdf[i] += cdf[i - 1];
            }
        } break;
        case yt_etype_quad: {
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i * 4;
                yt__vec3f *v0 = pos + f[0], *v1 = pos + f[1], *v2 = pos + f[2],
                          *v3 = pos + f[3];
                yt__vec3f e1 = yt__sub_vec3f(*v1, *v0);
                yt__vec3f e2 = yt__sub_vec3f(*v3, *v0);
                yt__vec3f a = yt__cross3f(e1, e2);
                cdf[i] = yt__length3f(a) / 2;
                if (f[2] != f[3]) {
                    e1 = yt__sub_vec3f(*v3, *v2);
                    e2 = yt__sub_vec3f(*v1, *v2);
                    a = yt__cross3f(e1, e2);
                    cdf[i] += yt__length3f(a) / 2;
                }
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
YT_API void
yt_init_lights(yt_scene* scene) {
    // clear old lights
    if (scene->nlights) {
        for (int i = 0; i < scene->nlights; i++) {
            yt__light* light = scene->lights + i;
            if (light->cdf) free(light->cdf);
        }
        free(scene->lights);
    }

    int nenvs = (scene->env && !yt__iszero3f(scene->env->ke)) ? 1 : 0;
    int nshapes = 0;
    for (int sid = 0; sid < scene->nshapes; sid++) {
        yt__shape* shape = scene->shapes + sid;
        if (!yt__iszero3f(shape->mat->ke)) nshapes++;
    }

    scene->nlights = nshapes + nenvs;
    if (!scene->nlights) {
        scene->lights = 0;
        return;
    }

    scene->lights = (yt__light*)calloc(scene->nlights, sizeof(yt__light));
    yt__light* light = scene->lights;

    for (int sid = 0; sid < scene->nshapes; sid++) {
        yt__shape* shape = scene->shapes + sid;
        if (yt__iszero3f(shape->mat->ke)) continue;
        light->shape = shape;
        light->cdf = (float*)calloc(shape->nelems, sizeof(float));
        yt__compute_shape_cdf(shape->nelems, shape->elem, shape->etype,
                              (yt__vec3f*)shape->pos, light->cdf, &light->area);
        light++;
    }

    if (nenvs) {
        light->env = scene->env;
    }
}

//
// Set material. Public API, see above.
//
YT_API void
yt_set_material(yt_scene* scene, int mid, float ke[3], float kd[3], float ks[3],
                float rs, float es[6], int ke_txt, int kd_txt, int ks_txt,
                int rs_txt, bool use_phong) {
    yt__material* material = scene->materials + mid;
    material->ke = *(yt__vec3f*)ke;
    material->kd = *(yt__vec3f*)kd;
    material->ks = *(yt__vec3f*)ks;
    if (es) {
        material->es = *(yt__vec3f*)es;
        material->eks = *(yt__vec3f*)(es + 3);
    } else {
        material->es = yt__zero3f();
        material->eks = yt__zero3f();
    }
    material->rs = rs;
    material->ke_txt = (ke_txt < 0) ? 0 : scene->textures + ke_txt;
    material->kd_txt = (kd_txt < 0) ? 0 : scene->textures + kd_txt;
    material->ks_txt = (ks_txt < 0) ? 0 : scene->textures + ks_txt;
    material->rs_txt = (rs_txt < 0) ? 0 : scene->textures + rs_txt;
    material->use_phong = use_phong;
}

//
// Phong exponent to roughness. Public API, see above.
//
YT_API float
yt_specular_exponent_to_roughness(float n) {
    return sqrtf(2 / (n + 2));
}

//
// Specular to fresnel eta. Public API, see above.
//
YT_API void
yt_specular_fresnel_from_ks(const float ks[3], float es[6]) {
    for (int c = 0; c < 3; c++) {
        es[c] = (1 + sqrtf(ks[c])) / (1 - sqrtf(ks[c]));
        es[c + 3] = 0;
    }
}

//
// Set env. Public API, see above.
//
YT_API void
yt_set_env(yt_scene* scene, float xform[16], float ke[3], int ke_txt) {
    if (!scene->env) scene->env = (yt__env*)calloc(1, sizeof(yt__env));
    yt__env* env = scene->env;
    env->xform = (xform) ? *(yt__mat4f*)xform : yt__identity4f();
    env->xform_inv =
        (xform) ? yt__inverse4f(*(yt__mat4f*)(xform)) : yt__identity4f();
    env->ke = *(yt__vec3f*)ke;
    env->ke_txt = (ke_txt < 0) ? 0 : scene->textures + ke_txt;
}

//
// Set texture. Public API, see above.
//
YT_API void
yt_set_texture(yt_scene* scene, int tid, float* pixels, int w, int h, int nc) {
    yt__texture* txt = scene->textures + tid;
    txt->w = w;
    txt->h = h;
    txt->nc = nc;
    txt->pixels = pixels;
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// Computes the i-th term of a permutation of l values keyed by p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
static inline uint32_t
yt__cmjs_permute(uint32_t i, uint32_t l, uint32_t p) {
    uint32_t w = l - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;
    do {
        i ^= p;
        i *= 0xe170893du;
        i ^= p >> 16;
        i ^= (i & w) >> 4;
        i ^= p >> 8;
        i *= 0x0929eb3f;
        i ^= p >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | p >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;
    } while (i >= l);
    return (i + p) % l;
}

//
// Computes a float value by hashing i with a key p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
static inline float
yt__cmjs_randfloat(uint32_t i, uint32_t p) {
    i ^= p;
    i ^= i >> 17;
    i ^= i >> 10;
    i *= 0xb36534e5;
    i ^= i >> 12;
    i ^= i >> 21;
    i *= 0x93fc4795;
    i ^= 0xdf6e307f;
    i ^= i >> 17;
    i *= 1 | p >> 18;
    return i * (1.0f / 4294967808.0f);
}

//
// 64 bit integer hash. Public domain code.
//
static inline uint64_t
yt__hashint64(uint64_t a) {
    a = (~a) + (a << 21);  // a = (a << 21) - a - 1;
    a ^= (a >> 24);
    a += (a << 3) + (a << 8);  // a * 265
    a ^= (a >> 14);
    a += (a << 2) + (a << 4);  // a * 21
    a ^= (a >> 28);
    a += (a << 31);
    return a;
}

//
// 64-to-32 bit integer hash. Public domain code.
//
static inline uint32_t
yt__hashint64_to_32(uint64_t a) {
    a = (~a) + (a << 18);  // a = (a << 18) - a - 1;
    a ^= (a >> 31);
    a *= 21;  // a = (a + (a << 2)) + (a << 4);
    a ^= (a >> 11);
    a += (a << 6);
    a ^= (a >> 22);
    return (uint32_t)a;
}

//
// PCG random numbers. A family of random number generators that supports
// multiple sequences. In our code, we allocate one sequence for each sample.
// PCG32 from http://www.pcg-random.org/
//

//
// Random number state
//
typedef struct { uint64_t state, inc; } yt__pcg32;

//
// Next random number
//
static inline uint32_t
yt__pcg32_next(yt__pcg32* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ull + (rng->inc | 1u);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

//
// Init a random number generator with a state state from the sequence seq.
//
static inline void
yt__pcg32_init(yt__pcg32* rng, uint64_t state, uint64_t seq) {
    rng->state = 0U;
    rng->inc = (seq << 1u) | 1u;
    yt__pcg32_next(rng);
    rng->state += state;
    yt__pcg32_next(rng);
}

//
// Next random float in [0,1).
//
static inline double
yt__pcg32_nextf(yt__pcg32* rng) {
    return (float)ldexp(yt__pcg32_next(rng), -32);
}

//
// Random number sampler. Handles random number generation for stratified
// sampling and correlated multi-jittered sampling.
//
typedef struct yt__sampler {
    yt__pcg32 rng;  // rnumber number state
    int i, j;       // pixel coordinates
    int s, d;       // sample and dimension indices
    int ns;         // number of samples
    int rtype;      // random number type
} yt__sampler;

//
// Initialize a sampler ot type rtype for pixel i, j with ns total samples.
//
// Implementation Notes: we use hash functions to scramble the pixel ids
// to avoid introducing unwanted correlation between pixels. These should not
// around according to the RNG documentaion, but we still found bad cases.
// Scrambling avoids it.
//
static inline yt__sampler
yt__make_sampler(int i, int j, int s, int ns, int rtype) {
    // we use various hashes to scramble the pixel values
    yt__sampler sampler = { { 0, 0 }, i, j, s, 0, ns, rtype };
    uint64_t sample_id = ((uint64_t)(i + 1)) << 0 | ((uint64_t)(j + 1)) << 15 |
                         ((uint64_t)(s + 1)) << 30;
    uint64_t initseq = yt__hashint64(sample_id);
    uint64_t initstate =
        yt__hashint64(sample_id * 3202034522624059733ull + 1ull);
    yt__pcg32_init(&sampler.rng, initstate, initseq);
    return sampler;
}

//
// Generates a 1-dimensional sample.
//
// Implementation Notes: For deterministic sampling (stratified and cmjs) we
// compute a 64bit sample and use hashing to avoid correlation. Then permutation
// are computed with CMJS procedures.
//
static inline float
yt__sample_next1f(yt__sampler* sampler) {
    float rn = 0;
    switch (sampler->rtype) {
        case yt_rtype_default:
        case yt_rtype_uniform: {
            rn = yt__pcg32_nextf(&sampler->rng);
        } break;
        case yt_rtype_stratified: {
            uint32_t p =
                yt__hashint64_to_32(((uint64_t)(sampler->i + 1)) << 0 |
                                    ((uint64_t)(sampler->j + 1)) << 15 |
                                    ((uint64_t)(sampler->d + 1)) << 30);
            int s = yt__cmjs_permute(sampler->s, sampler->ns, p);
            rn = (s + yt__pcg32_nextf(&sampler->rng)) / sampler->ns;
        } break;
        case yt_rtype_cmjs: {
            uint32_t p =
                yt__hashint64_to_32(((uint64_t)(sampler->i + 1)) << 0 |
                                    ((uint64_t)(sampler->j + 1)) << 15 |
                                    ((uint64_t)(sampler->d + 1)) << 30);
            int s = yt__cmjs_permute(sampler->s, sampler->ns, p);
            rn = (s + yt__cmjs_randfloat(s, p * 0xa399d265)) / sampler->ns;
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
static inline yt__vec2f
yt__sample_next2f(yt__sampler* sampler) {
    yt__vec2f rn = { 0, 0 };
    switch (sampler->rtype) {
        case yt_rtype_default:
        case yt_rtype_uniform: {
            rn.x = yt__pcg32_nextf(&sampler->rng);
            rn.y = yt__pcg32_nextf(&sampler->rng);
        } break;
        case yt_rtype_stratified: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler->ns));
            uint32_t p =
                yt__hashint64_to_32(((uint64_t)(sampler->i + 1)) << 0 |
                                    ((uint64_t)(sampler->j + 1)) << 15 |
                                    ((uint64_t)(sampler->d + 1)) << 30);
            int s = yt__cmjs_permute(sampler->s, sampler->ns, p);
            rn.x = (s % ns2 + yt__pcg32_nextf(&sampler->rng)) / ns2;
            rn.y = (s / ns2 + yt__pcg32_nextf(&sampler->rng)) / ns2;
        } break;
        case yt_rtype_cmjs: {
            uint32_t ns2 = (uint32_t)round(sqrt(sampler->ns));
            uint32_t p =
                yt__hashint64_to_32(((uint64_t)(sampler->i + 1)) << 0 |
                                    ((uint64_t)(sampler->j + 1)) << 15 |
                                    ((uint64_t)(sampler->d + 1)) << 30);
            int s = yt__cmjs_permute(sampler->s, sampler->ns, p);
            int sx = yt__cmjs_permute(s % ns2, ns2, p * 0xa511e9b3);
            int sy = yt__cmjs_permute(s / ns2, ns2, p * 0x63d83595);
            float jx = yt__cmjs_randfloat(s, p * 0xa399d265);
            float jy = yt__cmjs_randfloat(s, p * 0x711ad6a5);
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
typedef struct yt__point {
    // shape data -------------------------
    const yt__shape* shape;
    int eid;        // element id
    yt__vec2f euv;  // element baricentric coordinates
    // env data ---------------------------
    const yt__env* env;  // env
    // direction --------------------------
    yt__vec3f wo;  // outgoing direction
    // resolved geometry (shape) ----------
    int etype;               // element type
    yt__vec3f pos, norm;     // position and normal
    yt__vec3f tangx, tangy;  // tangent frame
    // shading ----------------------------
    yt__vec3f ke, kd, ks;  // material values
    float rs;              // material values
    yt__vec3f es, eks;     // material values
    bool use_phong;        // material values
    // sampling data ----------------------
    float wkd, wks, wrr;  // weights of kd and ks and total
} yt__point;

//
// Generates a ray ray_o, ray_d from a camera cam for image plane coordinate
// uv and the lens coordinates luv.
//
static inline void
yt__eval_camera(const yt__camera* cam, yt__vec2f uv, yt__vec2f luv,
                yt__vec3f* ray_o, yt__vec3f* ray_d) {
    yt__vec3f o =
        (yt__vec3f){ luv.x * cam->aperture, luv.y * cam->aperture, 0 };
    yt__vec3f q = { cam->width * cam->focus * (uv.x - 0.5f),
                    cam->height * cam->focus * (uv.y - 0.5f), -cam->focus };
    *ray_o = yt__transform_point3f(cam->xform, o);
    *ray_d = yt__normalize3f(
        yt__sub_vec3f(yt__transform_point3f(cam->xform, q), *ray_o));
}

//
// Evaluates a texture txt at texture coordinates texcoord. Uses bilinear
// interpolation and tiling. If the texture is semitrasparent, uses c
// as the case color.
//
static inline yt__vec3f
yt__eval_texture(const yt__vec3f c, const yt__texture* txt,
                 const yt__vec2f texcoord) {
    if (!txt) return c;
    assert(txt->nc == 3 || txt->nc == 4);
    yt__vec2f tc = { fmodf(texcoord.x, 1), fmodf(texcoord.y, 1) };
    if (tc.x < 0) tc.x += 1;
    if (tc.y < 0) tc.y += 1;
    int i = tc.x * txt->w, j = tc.y * txt->h;
    if (i > txt->w - 1) i = txt->w - 1;
    if (j > txt->h - 1) i = txt->h - 1;
    int ii = (i + 1) % txt->w, jj = (j + 1) % txt->h;
    float s = tc.x * txt->w - i, t = tc.y * txt->h - j;
    int idx[4] = { j * txt->w + i, jj * txt->w + i, j * txt->w + ii,
                   jj * txt->w + ii };
    float w[4] = { (1 - s) * (1 - t), (1 - s) * t, s * (1 - t), s * t };
    if (txt->nc == 3) {
        yt__vec3f* p3 = (yt__vec3f*)txt->pixels;
        yt__vec3f t = yt__zero3f();
        t = yt__sum3f(t, yt__smul3f(p3[idx[0]], w[0]));
        t = yt__sum3f(t, yt__smul3f(p3[idx[1]], w[1]));
        t = yt__sum3f(t, yt__smul3f(p3[idx[2]], w[2]));
        t = yt__sum3f(t, yt__smul3f(p3[idx[3]], w[3]));
        return t;
    } else if (txt->nc == 4) {
        yt__vec4f* p4 = (yt__vec4f*)txt->pixels;
        yt__vec4f t = { 0, 0, 0, 0 };
        t = yt__sum4f(t, yt__smul4f(p4[idx[0]], w[0]));
        t = yt__sum4f(t, yt__smul4f(p4[idx[1]], w[1]));
        t = yt__sum4f(t, yt__smul4f(p4[idx[2]], w[2]));
        t = yt__sum4f(t, yt__smul4f(p4[idx[3]], w[3]));
        return yt__sum3f(yt__smul3f((yt__vec3f){ t.x, t.y, t.z }, t.w),
                         yt__smul3f(c, 1 - t.w));
    } else {
        assert(false);
    }
}

//
// Evaluates emission.
//
static inline yt__vec3f
yt__eval_emission(const yt__point* pt) {
    if (yt__iszero3f(pt->ke)) return yt__zero3f();
    if (pt->shape) {
        switch (pt->shape->etype) {
            case yt_etype_point: return pt->ke;
            case yt_etype_line: return pt->ke;
            case yt_etype_triangle:
            case yt_etype_quad:
                if (yt__dot3f(pt->norm, pt->wo) > 0)
                    return pt->ke;
                else
                    return yt__zero3f();
            default: assert(false); return yt__zero3f();
        }
    } else if (pt->env) {
        return pt->ke;
    } else {
        return yt__zero3f();
    }
}

//
// Local-to-world direction transform for point frame.
//
static inline yt__vec3f
yt__transform_dirtoworld(const yt__point* pt, yt__vec3f wl) {
    yt__vec3f w = yt__zero3f();
    w = yt__sum3f(w, yt__smul3f(pt->tangx, wl.x));
    w = yt__sum3f(w, yt__smul3f(pt->tangy, wl.y));
    w = yt__sum3f(w, yt__smul3f(pt->norm, wl.z));
    w = yt__normalize3f(w);
    return w;
}

//
// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
//
static inline float
yt__eval_fresnel_dielectric(float cosw, float eta) {
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
static inline float
yt__eval_fresnel_metal(float cosw, float eta, float etak) {
    // TODO: vec3f
    if (!etak) return yt__eval_fresnel_dielectric(cosw, eta);

    cosw = yt__fclampf(cosw, -1, 1);
    float cos2 = cosw * cosw;
    float sin2 = yt__fclampf(1 - cos2, 0, 1);
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
static inline yt__vec3f
yt__eval_brdfcos(const yt__point* pt, yt__vec3f wi) {
    // summation over multiple terms
    yt__vec3f brdfcos = yt__zero3f();

    // exit if not needed
    if (yt__iszero3f(pt->kd) && yt__iszero3f(pt->ks)) return brdfcos;

    // save wo
    yt__vec3f wo = pt->wo;

    // compute wh
    yt__vec3f wh = yt__normalize3f(yt__sum3f(wo, wi));

    // compute dot products
    float ndo = yt__dot3f(pt->norm, wo), ndi = yt__dot3f(pt->norm, wi),
          ndh = yt__fclampf(yt__dot3f(wh, pt->norm), 0, 1);

    switch (pt->shape->etype) {
        case yt_etype_point: {
            // diffuse term (hack for now)
            if (!yt__iszero3f(pt->kd)) {
                float ido = yt__dot3f(wo, wi);
                yt__vec3f diff =
                    yt__smul3f(pt->kd, (2 * ido + 1) / (2 * yt__pif));
                brdfcos = yt__sum3f(brdfcos, diff);
            }
        } break;
        case yt_etype_line: {
            // take sines
            float so = sqrtf(yt__fclampf(1 - ndo * ndo, 0, 1)),
                  si = sqrtf(yt__fclampf(1 - ndi * ndi, 0, 1)),
                  sh = sqrtf(yt__fclampf(1 - ndh * ndh, 0, 1));

            // diffuse term (Kajiya-Kay)
            if (si > 0 && so > 0 && !yt__iszero3f(pt->kd)) {
                yt__vec3f diff = yt__smul3f(pt->kd, si / yt__pif);
                brdfcos = yt__sum3f(brdfcos, diff);
            }

            // specular term (Kajiya-Kay)
            if (si > 0 && so > 0 && sh > 0 && !yt__iszero3f(pt->ks)) {
                float ns = 2 / (pt->rs * pt->rs) - 2;
                float d = (ns + 2) * powf(sh, ns) / (2 + yt__pif);
                yt__vec3f spec = yt__smul3f(pt->ks, si * d / (4 * si * so));
                brdfcos = yt__sum3f(brdfcos, spec);
            }
        } break;
        case yt_etype_triangle:
        case yt_etype_quad: {
            // diffuse term
            if (ndi > 0 && ndo && !yt__iszero3f(pt->kd)) {
                yt__vec3f diff = yt__smul3f(pt->kd, ndi / yt__pif);
                brdfcos = yt__sum3f(brdfcos, diff);
            }

            // specular term (GGX)
            if (ndi > 0 && ndo > 0 && ndh > 0 && !yt__iszero3f(pt->ks)) {
                if (!pt->use_phong) {
                    // evaluate GGX
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = pt->rs * pt->rs;
                    float d = alpha2 / (yt__pif * cos2 * cos2 *
                                        (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o =
                        (-1 + sqrtf(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i =
                        (-1 + sqrtf(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    yt__vec3f spec =
                        yt__smul3f(pt->ks, ndi * d * g / (4 * ndi * ndo));
                    if (pt->es.x) {
                        // TODO: fixme
                        spec.x *=
                            yt__eval_fresnel_metal(ndh, pt->es.x, pt->eks.x);
                        spec.y *=
                            yt__eval_fresnel_metal(ndh, pt->es.y, pt->eks.y);
                        spec.z *=
                            yt__eval_fresnel_metal(ndh, pt->es.z, pt->eks.z);
                    }
                    brdfcos = yt__sum3f(brdfcos, spec);
                } else {
                    // evaluate Blinn-Phong
                    float ns = 2 / (pt->rs * pt->rs) - 2;
                    float d = (ns + 2) * powf(ndh, ns) / (2 + yt__pif);
                    yt__vec3f spec =
                        yt__smul3f(pt->ks, ndi * d / (4 * ndi * ndo));
                    brdfcos = yt__sum3f(brdfcos, spec);
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
static inline float
yt__weight_brdfcos(const yt__point* pt, yt__vec3f wi) {
    // skip if no component
    if (yt__iszero3f(pt->kd) && yt__iszero3f(pt->ks)) return 0;

    // save wo
    yt__vec3f wo = pt->wo;

    // compute wh
    yt__vec3f wh = yt__normalize3f(yt__sum3f(wi, wo));

    // compute dot products
    float ndo = yt__dot3f(pt->norm, wo), ndi = yt__dot3f(pt->norm, wi),
          ndh = yt__dot3f(pt->norm, wh);

    // check to make sure we are above the surface
    // updated this for refraction
    if (ndo <= 0 || ndi <= 0) return 0;

    // pick from a sum
    float wall = yt__mean3f(pt->kd) + yt__mean3f(pt->ks);
    float wd = yt__mean3f(pt->kd) / wall;
    float ws = yt__mean3f(pt->ks) / wall;

    // accumulate probability
    float pdf = 0;

    switch (pt->shape->etype) {
        case yt_etype_point:
        case yt_etype_line: {
            // diffuse term
            if (wall) {
                // homepherical cosine probability
                pdf += 1 / (4 * yt__pif);
            }
        } break;
        case yt_etype_triangle:
        case yt_etype_quad: {
            // diffuse term
            if (wd && ndi > 0) {
                // homepherical cosine probability
                pdf += wd * ndi / yt__pif;
            }

            // specular term (GGX or Phong)
            if (ws && ndi > 0 && ndo > 0 && ndh > 0) {
                if (!pt->use_phong) {
                    // probability proportional to d * ndh
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = pt->rs * pt->rs;
                    float d = alpha2 / (yt__pif * cos2 * cos2 *
                                        (alpha2 + tan2) * (alpha2 + tan2));
                    float hdo = yt__dot3f(wo, wh);
                    pdf += ws * d * ndh / (4 * hdo);
                } else {
                    // get phong exponent
                    float ns = 2 / (pt->rs * pt->rs) - 2;
                    // compute wh
                    yt__vec3f wh = yt__sum3f(wi, wo);
                    wh = yt__normalize3f(wh);
                    float ndh = yt__dot3f(pt->norm, wh);
                    // homerispherical cosine power probability
                    pdf += ws * powf(ndh, ns) * (ns + 1) / (2 * yt__pif);
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
static inline yt__vec3f
yt__sample_brdfcos(const yt__point* pt, float rnl, yt__vec2f rn) {
    // skip if no component
    if (yt__iszero3f(pt->kd) && yt__iszero3f(pt->ks)) return yt__zero3f();

    // save wo
    yt__vec3f wo = pt->wo;

    // compute cosine
    float ndo = yt__dot3f(pt->norm, wo);

    // check to make sure we are above the surface
    // update this for refraction
    if (ndo <= 0) return yt__zero3f();

    // pick from a sum
    float wall = yt__mean3f(pt->kd) + yt__mean3f(pt->ks);
    float wd = yt__mean3f(pt->kd) / wall;
    float ws = yt__mean3f(pt->ks) / wall;

    switch (pt->shape->etype) {
        case yt_etype_point:
        case yt_etype_line: {
            if (wall > 0) {
                // sample wi with uniform spherical distribution
                float rz = rn.y, rr = sqrtf(1 - rz * rz),
                      rphi = 2 * yt__pif * rn.x;
                yt__vec3f wi_local = { rr * cosf(rphi), rr * sinf(rphi), rz };
                return yt__transform_dirtoworld(pt, wi_local);
            }
        } break;
        case yt_etype_triangle:
        case yt_etype_quad: {
            // sample according to diffuse
            if (rnl < wd) {
                // sample wi with hemispherical cosine distribution
                float rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                      rphi = 2 * yt__pif * rn.x;
                // set to wi
                yt__vec3f wi_local = { rr * cosf(rphi), rr * sinf(rphi), rz };
                return yt__transform_dirtoworld(pt, wi_local);
            }

            // sample according to specular (GGX or Phong)
            if (rnl >= wd && rnl < wd + ws) {
                if (!pt->use_phong) {
                    // sample wh with ggx distribution
                    float tan2 = pt->rs * pt->rs * rn.y / (1 - rn.y);
                    float rz = sqrtf(1 / (tan2 + 1)), rr = sqrtf(1 - rz * rz),
                          rphi = 2 * yt__pif * rn.x;
                    // set to wh
                    yt__vec3f wh_local = { rr * cosf(rphi), rr * sinf(rphi),
                                           rz };
                    yt__vec3f wh = yt__transform_dirtoworld(pt, wh_local);
                    // compute wi
                    yt__vec3f wi = yt__smul3f(wh, 2 * yt__dot3f(wo, wh));
                    wi = yt__sub_vec3f(wi, wo);
                    wi = yt__normalize3f(wi);
                    return wi;
                } else {
                    // get phong exponent
                    float ns = 2 / (pt->rs * pt->rs) - 2;
                    // sample wh with hemispherical cosine power distribution
                    float rz = powf(rn.y, 1 / (ns + 1)),
                          rr = sqrtf(1 - rz * rz), rphi = 2 * yt__pif * rn.x;
                    // set to wh
                    yt__vec3f wh_local = { rr * cosf(rphi), rr * sinf(rphi),
                                           rz };
                    yt__vec3f wh = yt__transform_dirtoworld(pt, wh_local);
                    // compute wi
                    yt__vec3f wi = yt__smul3f(wh, 2 * yt__dot3f(wo, wh));
                    wi = yt__sub_vec3f(wi, wo);
                    wi = yt__normalize3f(wi);
                    return wi;
                }
            }
        } break;
        default: { assert(false); }
    }

    assert(false);

    return yt__zero3f();
}

//
// Create a point for an environment map. Resolves material with textures.
//
static inline yt__point
yt__eval_envpoint(const yt__env* env, yt__vec3f wo) {
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
        yt__vec3f w = yt__transform_direction3f(env->xform_inv, yt__neg3f(wo));
        float theta = 1 - (acosf(yt__fclampf(w.y, -1, 1)) / yt__pif);
        float phi = atan2f(w.z, w.x) / (2 * yt__pif);
        yt__vec2f texcoord = { phi, theta };
        pt.ke = yt__mul3f(pt.ke, yt__eval_texture((yt__vec3f){ 1, 1, 1 },
                                                  env->ke_txt, texcoord));
    }

    return pt;
}

//
// Create a point for a shape. Resolves geometry and material with textures.
//
static inline yt__point
yt__eval_shapepoint(const yt__shape* shape, int eid, yt__vec2f euv,
                    yt__vec3f wo) {
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
            int* f = shape->elem + eid * 3;
            w[0] = 1 - euv.x - euv.y;
            w[1] = euv.x;
            w[2] = euv.y;
            idx[0] = f[0];
            idx[1] = f[1];
            idx[2] = f[2];
            nidx = 3;
        }; break;
        case yt_etype_quad: {
            int* f = shape->elem + eid * 4;
            if (f[2] == f[3] || euv.x + euv.y <= 1) {
                idx[0] = f[0];
                idx[1] = f[1];
                idx[2] = f[3];
                w[0] = 1 - euv.x - euv.y;
                w[1] = euv.x;
                w[2] = euv.y;
            } else {
                idx[0] = f[2];
                idx[1] = f[3];
                idx[2] = f[1];
                w[0] = 1 - (1 - euv.x) - (1 - euv.y);
                w[1] = 1 - euv.x;
                w[2] = 1 - euv.y;
            }
            nidx = 3;
        }; break;
        case yt_etype_line: {
            int* f = shape->elem + eid * 2;
            w[0] = 1 - euv.x;
            w[1] = euv.x;
            idx[0] = f[0];
            idx[1] = f[1];
            nidx = 2;
        }; break;
        case yt_etype_point: {
            int f = shape->elem[eid];
            w[0] = 1;
            idx[0] = f;
            nidx = 1;
        }; break;
    }

    // interpolate geometry positions
    pt.pos = yt__zero3f();
    for (int i = 0; i < nidx; i++) {
        pt.pos.x += shape->pos[idx[i] * 3 + 0] * w[i];
        pt.pos.y += shape->pos[idx[i] * 3 + 1] * w[i];
        pt.pos.z += shape->pos[idx[i] * 3 + 2] * w[i];
    }

    pt.norm = yt__zero3f();
    if (shape->norm) {
        for (int i = 0; i < nidx && shape->norm; i++) {
            pt.norm.x += shape->norm[idx[i] * 3 + 0] * w[i];
            pt.norm.y += shape->norm[idx[i] * 3 + 1] * w[i];
            pt.norm.z += shape->norm[idx[i] * 3 + 2] * w[i];
        }
        pt.norm = yt__normalize3f(pt.norm);
    }

    // creating frame
    pt.tangx = yt__orthogonal3f(pt.norm);
    pt.tangx = yt__normalize3f(pt.tangx);
    pt.tangy = yt__cross3f(pt.norm, pt.tangx);
    pt.tangy = yt__normalize3f(pt.tangy);

    // transform to world space
    if (shape->xformed) {
        pt.pos = yt__transform_point3f(shape->xform, pt.pos);
        pt.norm = yt__transform_direction3f(shape->xform, pt.norm);
        pt.tangx = yt__transform_direction3f(shape->xform, pt.tangx);
        pt.tangy = yt__transform_direction3f(shape->xform, pt.tangy);
    }

    // sample material data
    yt__material* mat = shape->mat;
    pt.ke = mat->ke;
    pt.kd = mat->kd;
    pt.ks = mat->ks;
    pt.rs = mat->rs;
    pt.use_phong = mat->use_phong;

    // handle surface color
    if (shape->color) {
        yt__vec3f color = yt__zero3f();
        for (int i = 0; i < nidx; i++) {
            color.x += shape->color[idx[i] * 3 + 0] * w[i];
            color.y += shape->color[idx[i] * 3 + 1] * w[i];
            color.z += shape->color[idx[i] * 3 + 2] * w[i];
        }

        pt.ke = yt__mul3f(pt.ke, color);
        pt.kd = yt__mul3f(pt.kd, color);
        pt.ks = yt__mul3f(pt.ks, color);
    }

    // handle textures
    if (shape->texcoord) {
        yt__vec2f texcoord = { 0, 0 };
        for (int i = 0; i < nidx; i++) {
            texcoord.x += shape->texcoord[idx[i] * 2 + 0] * w[i];
            texcoord.y += shape->texcoord[idx[i] * 2 + 1] * w[i];
        }

        if (mat->ke_txt) {
            pt.ke = yt__mul3f(pt.ke, yt__eval_texture((yt__vec3f){ 1, 1, 1 },
                                                      mat->ke_txt, texcoord));
        }
        if (mat->kd_txt) {
            pt.kd = yt__mul3f(pt.kd, yt__eval_texture((yt__vec3f){ 1, 1, 1 },
                                                      mat->kd_txt, texcoord));
        }
        if (mat->ks_txt) {
            pt.ks = yt__mul3f(pt.ks, yt__eval_texture((yt__vec3f){ 1, 1, 1 },
                                                      mat->ks_txt, texcoord));
        }
        if (mat->rs_txt) {
            pt.rs *=
                yt__eval_texture((yt__vec3f){ 1, 1, 1 }, mat->rs_txt, texcoord)
                    .x;
        }
    }

    return pt;
}

//
// Sample weight for a light point.
//
static inline float
yt__weight_light(const yt__light* light, const yt__point* lpt,
                 const yt__point* pt) {
    if (light->shape) {
        switch (light->shape->etype) {
            case yt_etype_point: {
                float d = yt__dist3f(lpt->pos, pt->pos);
                return light->area / (d * d);
            } break;
            case yt_etype_line: {
                assert(false);
                return 0;
            } break;
            case yt_etype_triangle:
            case yt_etype_quad: {
                float d = yt__dist3f(lpt->pos, pt->pos);
                return light->area * fabsf(yt__dot3f(lpt->norm, lpt->wo)) /
                       (d * d);
            } break;
            default: { assert(false); } break;
        }
    } else if (light->env) {
        return 4 * yt__pif;
    } else {
        assert(false);
    }
    return 0;
}

//
// Picks a point on a light.
//
static inline yt__point
yt__sample_light(const yt__light* light, const yt__point* pt, float rne,
                 yt__vec2f rn) {
    if (light->shape) {
        int eid;
        for (eid = 0; eid < light->shape->nelems && light->cdf[eid] < rne;
             eid += 1) {
        }
        if (eid > light->shape->nelems - 1) eid = light->shape->nelems - 1;
        yt__vec2f euv = rn;
        if (light->shape->etype == yt_etype_triangle) {
            euv.x = 1 - sqrt(rn.x);
            euv.y = rn.y * sqrt(rn.x);
        }
        yt__point lpt =
            yt__eval_shapepoint(light->shape, eid, euv, yt__zero3f());
        lpt.wo = yt__normalize3f(yt__sub_vec3f(pt->pos, lpt.pos));
        return lpt;
    } else if (light->env) {
        float z = -1 + 2 * rn.y;
        float rr = sqrtf(yt__fclampf(1 - z * z, 0, 1));
        float phi = 2 * yt__pif * rn.x;
        yt__vec3f wo = { cosf(phi) * rr, z, sinf(phi) * rr };
        yt__point lpt = yt__eval_envpoint(light->env, wo);
        return lpt;
    } else {
        assert(false);
    }
    return (yt__point){ 0 };
}

//
// Offsets a ray origin to avoid self-intersection.
//
static inline yt__vec3f
yt__offset_ray(const yt_scene* scene, const yt__point* pt) {
    return yt__sum3f(pt->pos, yt__smul3f(pt->norm, scene->ray_eps));
}

//
// Intersects a ray with the scene and return the point (or env point).
//
static inline yt__point
yt__intersect_scene(const yt_scene* scene, const yt__vec3f ray_o,
                    const yt__vec3f ray_d) {
    int sid, eid;
    yt__vec2f euv;
    float ray_t;
    bool hit =
        scene->intersect_ray(scene->ray_ctx, &ray_o.x, &ray_d.x, scene->ray_eps,
                             HUGE_VALF, 0, &ray_t, &sid, &eid, &euv.x);
    yt__vec3f wo = yt__neg3f(ray_d);
    if (hit) {
        return yt__eval_shapepoint(scene->shapes + sid, eid, euv, wo);
    } else {
        return yt__eval_envpoint(scene->env, wo);
    }
}

//
// Evalutes direct illumination using MIS.
//
static inline yt__vec3f
yt__eval_direct(const yt_scene* scene, const yt__light* light,
                const yt__point* pt, yt__sampler* sampler) {
    // select whether it goes in all light mode
    bool all_lights = !light;

    // pick a light if not there
    float nlweight = 0;
    if (all_lights) {
        int lid = yt__sample_next1f(sampler) * scene->nlights;
        if (lid > scene->nlights - 1) lid = scene->nlights - 1;
        light = scene->lights + lid;
        nlweight = scene->nlights;
    } else {
        nlweight = 1;
    }

    // sample light according to area
    yt__point lpt = yt__sample_light(light, pt, yt__sample_next1f(sampler),
                                     yt__sample_next2f(sampler));
    yt__vec3f lld = yt__mul3f(yt__eval_emission(&lpt),
                              yt__eval_brdfcos(pt, yt__neg3f(lpt.wo)));
    float lweight = yt__weight_light(light, &lpt, pt) * nlweight;
    lld = yt__smul3f(lld, lweight);
    if (!yt__iszero3f(lld)) {
        yt__vec3f shadow_ray_o = yt__offset_ray(scene, pt);
        yt__vec3f shadow_ray_d = yt__neg3f(lpt.wo);
        float shadow_ray_dist =
            (lpt.shape) ? yt__dist3f(shadow_ray_o, lpt.pos) : HUGE_VALF;
        if (scene->hit_ray(scene->ray_ctx, &shadow_ray_o.x, &shadow_ray_d.x,
                           scene->ray_eps, shadow_ray_dist - 2 * scene->ray_eps,
                           0))
            lld = yt__zero3f();
    }

    // check if mis is necessary
    if (pt->shape && (pt->shape->etype == yt_etype_point ||
                      pt->shape->etype == yt_etype_line))
        return lld;

    // check if mis is necessary
    if (light->shape && (light->shape->etype == yt_etype_point ||
                         light->shape->etype == yt_etype_line))
        return lld;

    // sample the brdf
    yt__vec3f bwi = yt__sample_brdfcos(pt, yt__sample_next1f(sampler),
                                       yt__sample_next2f(sampler));
    float bweight = 0;
    yt__vec3f bld = yt__zero3f();
    yt__point bpt = yt__intersect_scene(scene, yt__offset_ray(scene, pt), bwi);
    if (light->shape == bpt.shape || all_lights) {
        bweight = yt__weight_brdfcos(pt, bwi);
        bld = yt__mul3f(yt__eval_emission(&bpt), yt__eval_brdfcos(pt, bwi));
        bld = yt__smul3f(bld, bweight);
    }

    // accumulate the value with mis
    if (!yt__iszero3f(lld)) {
        float bweight = yt__weight_brdfcos(pt, yt__neg3f(lpt.wo));
        // float weight =
        //     (1 / lweight) * (1 / lweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        float weight = (1 / lweight) / ((1 / lweight) + (1 / bweight));
        lld = yt__smul3f(lld, weight);
    }
    if (!yt__iszero3f(bld)) {
        float lweight = yt__weight_light(light, &bpt, pt) * nlweight;
        // float weight =
        //     (1 / bweight) * (1 / bweight) /
        //     ((1 / lweight) * (1 / lweight) + (1 / bweight) * (1 / bweight));
        float weight = (1 / bweight) / ((1 / lweight) + (1 / bweight));
        bld = yt__smul3f(bld, weight);
    }

    // return weighted sum
    return yt__sum3f(lld, bld);
}

//
// Recursive path tracing.
//
static inline yt__vec3f
yt__shade_pathtrace_recd(const yt_scene* scene, yt__vec3f ray_o,
                         yt__vec3f ray_d, yt__sampler* sampler, bool* hit,
                         int ray_depth) {
    yt__point pt = yt__intersect_scene(scene, ray_o, ray_d);
    if (hit) *hit = pt.shape;

    yt__vec3f l = (!ray_depth) ? yt__eval_emission(&pt) : yt__zero3f();
    if (!pt.shape) return l;

    if (yt__iszero3f(pt.kd) && yt__iszero3f(pt.ks)) return l;

    yt__vec3f ld = yt__eval_direct(scene, 0, &pt, sampler);
    l = yt__sum3f(l, ld);

    if (ray_depth >= scene->max_ray_depth) return l;

    // roussian roulette
    float rrweight = 1;
    if (ray_depth >= scene->min_ray_depth) {
        float rrrn = yt__sample_next1f(sampler);
        if (rrrn >= fmin(pt.wrr, 0.95f)) return l;
        rrweight = 1.0f / fmin(pt.wrr, 0.95f);
    }

    // continue path
    yt__vec3f bwi = yt__sample_brdfcos(&pt, yt__sample_next1f(sampler),
                                       yt__sample_next2f(sampler));
    if (yt__iszero3f(bwi)) return l;
    float bweight = yt__weight_brdfcos(&pt, bwi);
    if (!bweight) return l;
    yt__vec3f bbrdfcos = yt__eval_brdfcos(&pt, bwi);
    if (yt__iszero3f(bbrdfcos)) return l;
    yt__vec3f ble = yt__shade_pathtrace_recd(scene, yt__offset_ray(scene, &pt),
                                             bwi, sampler, 0, ray_depth + 1);
    ble = yt__mul3f(ble, bbrdfcos);
    ble = yt__smul3f(ble, bweight * rrweight);
    l = yt__sum3f(l, ble);

    return l;
}

//
// Shader interface for the above function.
//
static inline yt__vec3f
yt__shade_pathtrace(const yt_scene* scene, yt__vec3f ray_o, yt__vec3f ray_d,
                    yt__sampler* sampler, bool* hit) {
    return yt__shade_pathtrace_recd(scene, ray_o, ray_d, sampler, hit, 0);
}

//
// Direct illuination.
//
static inline yt__vec3f
yt__shade_direct(const yt_scene* scene, yt__vec3f ray_o, yt__vec3f ray_d,
                 yt__sampler* sampler, bool* hit) {
    yt__point pt = yt__intersect_scene(scene, ray_o, ray_d);
    if (hit) *hit = pt.shape;

    yt__vec3f l = yt__eval_emission(&pt);
    if (!pt.shape) return l;

    if (yt__iszero3f(pt.kd) && yt__iszero3f(pt.ks)) return l;

    if (!yt__iszero3f(scene->amb)) {
        yt__vec3f ale = yt__mul3f(scene->amb, pt.kd);
        l = yt__sum3f(l, ale);
    }

    for (int lid = 0; lid < scene->nlights; lid++) {
        yt__light* light = scene->lights + lid;
        yt__vec3f ld = yt__eval_direct(scene, light, &pt, sampler);
        l = yt__sum3f(l, ld);
    }

    return l;
}

//
// Eyelight for quick previewing.
//
static inline yt__vec3f
yt__shade_eyelight(const yt_scene* scene, yt__vec3f ray_o, yt__vec3f ray_d,
                   yt__sampler* sampler, bool* hit) {
    yt__point pt = yt__intersect_scene(scene, ray_o, ray_d);
    if (hit) *hit = pt.shape;

    yt__vec3f l = yt__eval_emission(&pt);
    if (!pt.shape) return l;

    yt__vec3f brdfcos = yt__smul3f(yt__eval_brdfcos(&pt, pt.wo), yt__pif);
    l = yt__sum3f(l, brdfcos);

    return l;
}

//
// Shader function callback.
//
typedef yt__vec3f (*shade_fn)(const yt_scene* scene, yt__vec3f ray_o,
                              yt__vec3f ray_d, yt__sampler* sampler, bool* hit);

//
// Renders a block of pixels. Public API, see above.
//
YT_API void
yt_trace_block(const yt_scene* scene, float* img_pixels, int img_w, int img_h,
               int ns, int block[6]) {
    const yt__camera* cam = scene->camera;
    shade_fn shade;
    switch (scene->stype) {
        case yt_stype_eyelight: shade = yt__shade_eyelight; break;
        case yt_stype_default:
        case yt_stype_direct: shade = yt__shade_direct; break;
        case yt_stype_pathtrace: shade = yt__shade_pathtrace; break;
        default: assert(false); return;
    }
    for (int j = block[2]; j < block[3]; j++) {
        for (int i = block[0]; i < block[1]; i++) {
            float* pixel = img_pixels + (img_w * j + i) * 4;
            for (int i = 0; i < 4; i++) pixel[i] = 0;
            for (int s = block[4]; s < block[5]; s++) {
                yt__sampler sampler =
                    yt__make_sampler(i, j, s, ns, scene->rtype);
                yt__vec2f rn = yt__sample_next2f(&sampler);
                yt__vec2f uv = { (i + rn.x) / img_w, 1 - (j + rn.y) / img_h };
                yt__vec3f ray_o, ray_d;
                bool hit;
                yt__eval_camera(cam, uv, yt__sample_next2f(&sampler), &ray_o,
                                &ray_d);
                yt__vec3f l = shade(scene, ray_o, ray_d, &sampler, &hit);
                if (!yt__isfinite3f(l)) continue;
                if (scene->pixel_clamp > 0) yt__clamp3f(l, scene->pixel_clamp);
                pixel[0] += l.x;
                pixel[1] += l.y;
                pixel[2] += l.z;
                // pixel[3] += (hit) ? 1.0f : 0.0f;
                pixel[3] += 1;
            }
            for (int i = 0; i < 4; i++) pixel[i] /= block[5] - block[4];
        }
    }
}

//
// Renders the whole image. Public API, see above.
//
YT_API void
yt_trace_image(const yt_scene* scene, float* img_pixels, int img_w, int img_h,
               int ns) {
    int block[6] = { 0, img_w, 0, img_h, 0, ns };
    yt_trace_block(scene, img_pixels, img_w, img_h, ns, block);
}

#endif

#endif
