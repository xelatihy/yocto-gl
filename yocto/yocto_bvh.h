//
// YOCTO_BVH: ray-intersection and closet-point routines supporting points,
// lines, triangles accelerated by a two-level bounding volume
// hierarchy (BVH)
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. create the scene bvh
//      scene = yb_init_scene(mnshapes)
// 2. for each shape, add shape data
//     for(int i = 0; i < nshapes; i ++)
//         yb_set_shape(scene, pass shape data)
// 3. build the bvh with yb_build_scene
// 4.a. perform ray-interseciton tests
//     - use yb_intersect_first if you want to know the hit point
//         hit = yb_intersect_first(scene, ray data, out primitive intersection)
//     - use yb_intersect_any if you only need to know whether there is a hit
//         hit = yb_intersect_any(scene, ray data)
// 4.b. perform closet-point tests
//     - use yb_overlap_first to get the closet element to a point
//       bounded by a radius
//         hit = yb_overlap_first(scene, pos, radius, out primitive)
//     - use yb_overlap_any to check whether there is an element that overlaps
//         a point within a given radius
// 5. use yb_interpolate_vertex to get interpolated vertex values from the
//    intersection data
//    yb_interpolate_vertex(intersection data, shape data, out interpolated val)
// 6. use yb_refit_bvh to recompute the bvh bounds if objects move (you should
//    rebuild the bvh for large changes)
// 7. cleanup with yb_free_bvh
//   - you have to do this for each shape bvh and the scene bvh
//   yb_free_scene(scene)
//
// The interface for each function is described in details in the interface
// section of this file.
//
// You can also just use one untransformed shape by calling yb_init_shape_bvh
// and yb_build_shape_bvh, etc.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles),
// an array of vertex positions, and an array of vertex radius
// (for points and lines).
//
// Implementation notes:
// - using high precision bvh intersection by default
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
// - v 0.3: cleaner interface based on opaque yb_scene
// - v 0.2: better internal intersections and cleaner bvh stats
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
//

#ifndef _YB_H_
#define _YB_H_

// compilation options
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
// Element types for shapes
//
enum {
    yb_etype_point = 1,     // points
    yb_etype_line = 2,      // lines
    yb_etype_triangle = 3,  // triangles
};

//
// Heuristic stratigy for bvh build
//
enum {
    yb_htype_default = 0,  // default strategy (use this for ray-casting)
    yb_htype_equalnum,     // balanced binary tree
    yb_htype_equalsize,    // splitting space binary tree
    yb_htype_sah,          // surface area heuristic (full sweep for accuracy)
    yb_htype_binned_sah,   // surface area heuristic (binned for speed)
    yb_htype_max           // total number of strategies
};

//
// Scene data structure. Implementation details hidden.
//
struct yb_scene;

//
// Create a BVH for a collection of transformed shapes (scene).
//
// Shapes' BVHs can be transformed with transformations matrices. The code
// supports only affine transforms.
//
// Parameters:
// - nshapes: number of shape bvhs
//
YGL_API yb_scene* yb_init_scene(int nshapes);

//
// Sets shape data for scene BVH. Equivalent to calling init_shape_bvh.
//
// Parameters:
// - bvh: scene bvh
// - sid: shape id
// - xform: shape transform
// - nelems: number of elements
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - pos: array of 3D vertex positions
// - radius: array of vertex radius
//
YGL_API void yb_set_shape(yb_scene* scene, int sid, const ym_affine3f& xform,
                          int nelems, const int* elem, int etype, int nverts,
                          const ym_vec3f* pos, const float* radius);

//
// Builds a shape BVH.
//
// Parameters:
// - bvh: bvh to build
// - heuristic: heristic use to build the bvh (0 for deafult)
//
YGL_API void yb_build_bvh(yb_scene* scene, int heuristic);

//
// Clears BVH memory.
//
// This will only clear the memory for the given bvh. For scene bvhs, you are
// still responsible for clearing all shape bvhs.
//
// Parameters:
// - bvh: bvh to clean
//
YGL_API void yb_free_scene(yb_scene* scene);

//
// Refit the bounds of each shape for moving objects. Use this only to avoid
// a rebuild, but note that queries are likely slow if objects move a lot.
//
// Parameters:
// - scene: scene to refit
// - bvhs: array of pointers to shape bvhs
//
YGL_API void yb_refit_scene_bvh(yb_scene* scene, const ym_affine3f* xforms);

//
// Intersect the scene with a ray finding the closest intersection.
//
// Parameters:
// - scene: scene to intersect (scene or shape)
// - ray: ray origin, direction and min/max distance
// - req_shape: [optional] whether to restrict to only a particular shape
//   (-1 for all)
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
YGL_API bool yb_intersect_first(const yb_scene* scene, const ym_ray3f& ray,
                                float* ray_t, int* sid, int* eid, ym_vec2f* euv,
                                int req_shape = -1);

//
// Intersect the scene with a ray finding any intersection (good for shadow
// rays).
//
// Parameters:
// - scene: scene to intersect (scene or shape)
// - ray: ray origin, direction and min/max distance
// - req_shape: [optional] whether to restrict to only a particular shape
//   (-1 for all)
//
// Return:
// - whether we intersect or not
//
YGL_API bool yb_intersect_any(const yb_scene* scene, const ym_ray3f& ray,
                              int req_shape = -1);

//
// Returns a list of shape pairs that can possibly overlap by checking only they
// axis aligned bouds. This is only a conservative check useful for collision
// detection.
//
// Parameters:
// - scene: scene
// - exclude_self: whether to exlude self intersections
// - overlap_cb: function called for each overlap
//
// Return:
// - number of overlaps
//
// Notes:
// - intersections are duplicated, so if (i,j) overlaps than both (i,j) and
// (j,i)
//   will be present; this makes it easier to apply asymmetric checks
// - to remove symmetric checks, just skip all pairs with i > j
//
YGL_API int yb_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                    void* overlap_ctx,
                                    void (*overlap_cb)(const ym_vec2i& v));

//
// Returns a list of shape pairs that can possibly overlap by checking only they
// axis aligned bouds. This is only a conservative check useful for collision
// detection.
//
// Parameters:
// - scene: scene
// - exclude_self: whether to exlude self intersections
// - overlaps: possible shape overlaps
//
// Return:
// - number of overlaps
//
// Notes: See above.
//
YGL_API int yb_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                    ym_vector<ym_vec2i>* overlaps);

//
// Finds the closest element that overlaps a point within a given radius.
//
// Parameters:
// - scene: scene to check
// - pt: ray origin
// - max_dist: max point distance
// - req_shape: [optional] required shape (-1 for all)
//
// Out Parameters:
// - dist: distance
// - sid: shape index
// - eid: element index
// - euv: element parameters
//
// Return:
// - whether we overlap or not
//
YGL_API bool yb_overlap_first(const yb_scene* scene, const ym_vec3f& pt,
                              float max_dist, float* dist, int* sid, int* eid,
                              ym_vec2f* euv, int req_shape = -1);

//
// Finds any element that overlaps a point wihtin a given radius.
//
// Parameters:
// - scene: scene to check
// - pt: ray origin
// - max_dist: max point distance
// - req_shape: [optional] required shape (-1 for all)
//
// Return:
// - whether we overlap or not
//
YGL_API bool yb_overlap_any(const yb_scene* scene, const ym_vec3f& pt,
                            float max_dist, int req_shape = -1);

//
// Interpolates a vertex property from the given intersection data. Uses
// linear interpolation for lines, baricentric for triangles and quads
// (to correctly treat as two triangles) and copies values for points.
//
// Parameters:
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - eid: hit element index
// - euv: hit element parameters
// - vsize: number of floats in the vertex property
// - vert: shape vertex data (contigous array of vsize-sized values)
//
// Out Parameters:
// - v: interpolated vertex (array of vsize size)
//
YGL_API void yb_interpolate_vert(const int* elem, int etype, int eid,
                                 const ym_vec2f& euv, int vsize,
                                 const float* vert, float* v);

// -----------------------------------------------------------------------------
// PERFORMANCE TUNING INTERFACE
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
YGL_API void yb_compute_bvh_stats(const yb_scene* scene, bool include_shapes,
                                  int* nprims, int* ninternals, int* nleaves,
                                  int* min_depth, int* max_depth,
                                  int req_shape = -1);

//
// Logging (not thread safe to avoid affecting speed)
//
#ifdef YGL_BVH_LOG_RAYS
YGL_API void yb_get_ray_log(int* nrays, int* nbbox_inters, int* npoint_inters,
                            int* nline_inters, int* ntriangle_inters);
#endif

#endif  // __cplusplus

// -----------------------------------------------------------------------------
// C/C++ INTERFACE
// -----------------------------------------------------------------------------

//
// Scene data structure. Implementation details hidden.
//
typedef struct yb_scene yb_scene;

//
// Create a BVH for a collection of transformed shapes (scene).
//
// Shapes' BVHs can be transformed with transformations matrices. The code
// supports only affine transforms.
//
// Parameters:
// - nshapes: number of shapes
//
YGLC_API yb_scene* ybc_init_scene(int nshapes);

//
// Sets shape data for scene BVH. Equivalent to calling init_shape_bvh.
//
// Parameters:
// - scene: scene
// - sid: shape id
// - xform: shape transform
// - nelems: number of elements
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - pos: array of 3D vertex positions
// - radius: array of vertex radius
//
YGLC_API void ybc_set_shape(yb_scene* scene, int sid, const float xform[16],
                            int nelems, const int* elem, int etype, int nverts,
                            const float* pos, const float* radius);

//
// Clears BVH memory.
//
// This will only clear the memory for the given bvh. For scene bvhs, you are
// still responsible for clearing all shape bvhs.
//
// Parameters:
// - scene: scene to clean
//
YGLC_API void ybc_free_scene(yb_scene* scene);

//
// Refit the bounds of each shape for moving objects. Use this only to avoid
// a rebuild, but note that queries are likely slow if objects move a lot.
//
// Parameters:
// - scene: scene to refit
// - bvhs: array of pointers to shape bvhs
//
YGLC_API void ybc_refit_bvh(yb_scene* scene, const float* xforms[16]);

//
// Intersect the scene with a ray finding the closest intersection.
//
// Parameters:
// - scene: scene to intersect
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - req_shape: required shape (or -1 for all scene)
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
YGLC_API bool ybc_intersect_first(const yb_scene* scene, const float ray_o[3],
                                  const float ray_d[3], float ray_tmin,
                                  float ray_tmax, float* ray_t, int* sid,
                                  int* eid, float euv[2], int req_shape);

//
// Intersect the scene with a ray finding any intersection.
//
// Parameters:
// - bvh: bvh to intersect (scene or shape)
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - req_shape: required shape (or -1 for all scene)
//
// Return:
// - whether we intersect or not
//
YGLC_API bool ybc_intersect_any(const yb_scene* scene, const float ray_o[3],
                                const float ray_d[3], float ray_tmin,
                                float ray_tmax, int req_shape);

//
// Returns a list of shape pairs that can possibly overlap by checking only they
// axis aligned bouds. This is only a conservative check useful for collision
// detection.
//
// Parameters:
// - bvh: scene bvh
// - exclude_self: whether to exlude self intersections
// - overlaps: pointer to an array of shape pairs (consecutive ints)
// - ocapacity: capacity of the overlap array to avoid reallocations
//
// Return:
// - number of intersections
//
// Notes:
// - intersections are duplicated, so if (i,j) overlaps than both (i,j) and
// (j,i)
//   will be present; this makes it easier to apply asymmetric checks
// - to remove symmetric checks, just skip all pairs with i > j
//
YGLC_API int ybc_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                      void* overlap_ctx,
                                      void (*overlap_cb)(void* ctx,
                                                         const int ij[2]));

//
// Finds the closest element to a point wihtin a given radius.
//
// Parameters:
// - scene: scene to overlap (scene or shape)
// - pt: ray origin
// - max_dist: max point distance
// - req_shape: required shape (-1 for all)
//
// Out Parameters:
// - dist: distance
// - sid: shape index
// - eid: element index
// - euv: element parameters
//
// Return:
// - whether we overlap or not
//
YGLC_API bool ybc_overlap_first(const yb_scene* scene, const float pt[3],
                                float max_dist, float* dist, int* sid, int* eid,
                                float euv[2], int req_shape);

//
// Finds if any element overlaps a point wihtin a given radius.
//
// Parameters:
// - scene: scene to overlap
// - pt: origin
// - max_dist: max point distance
// - req_shape: required shape (-1 for all)
//
// Out Parameters:
// - dist: distance
// - sid: shape index
// - eid: element index
// - euv: element parameters
//
// Return:
// - whether we overlap or not
//
YGLC_API bool ybc_overlap_any(const yb_scene* scene, const float pt[3],
                              float max_dist, int req_shape);

//
// Interpolates a vertex property from the given intersection data. Uses
// linear interpolation for lines, baricentric for triangles and quads
// (to correctly treat as two triangles) and copies values for points.
//
// Parameters:
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - eid: hit element index
// - euv: hit element parameters
// - vsize: number of floats in the vertex property
// - vert: shape vertex data (contigous array of vsize-sized values)
//
// Out Parameters:
// - v: interpolated vertex (array of vsize size)
//
YGLC_API void ybc_interpolate_vert(const int* elem, int etype, int eid,
                                   const float euv[2], int vsize,
                                   const float* vert, float* v);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

// handles compilation options
#if defined(__cplusplus) &&                                                    \
    (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

#include <algorithm>
#include <cstdio>

// -----------------------------------------------------------------------------
// MATH FUNCTIONS SUPPORT
// -----------------------------------------------------------------------------

//
// Transforms a bbox via an affine matrix
//
static inline ym_range3f yb__transform_bbox(const ym_affine3f& a,
                                            const ym_range3f& bbox) {
    ym_vec3f corners[8] = {
        {bbox.min.x, bbox.min.y, bbox.min.z},
        {bbox.min.x, bbox.min.y, bbox.max.z},
        {bbox.min.x, bbox.max.y, bbox.min.z},
        {bbox.min.x, bbox.max.y, bbox.max.z},
        {bbox.max.x, bbox.min.y, bbox.min.z},
        {bbox.max.x, bbox.min.y, bbox.max.z},
        {bbox.max.x, bbox.max.y, bbox.min.z},
        {bbox.max.x, bbox.max.y, bbox.max.z},
    };
    ym_range3f xformed = ym_invalid_range3f;
    for (int j = 0; j < 8; j++) xformed += ym_transform_point(a, corners[j]);
    return xformed;
}

//
// Compute the point on the ray ray_o, ray_d at distance t
//
static inline ym_vec3f yb__eval_ray(const ym_vec3f ray_o, const ym_vec3f ray_d,
                                    float t) {
    return {ray_o.x + ray_d.x * t, ray_o.y + ray_d.y * t,
            ray_o.z + ray_d.z * t};
}

// -----------------------------------------------------------------------------
// ELEMENT-WISE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Intersect a ray with a point (approximate)
//
// Parameters:
// - ray: ray origin and direction, parameter min, max range
// - p: point position
// - r: point radius
//
// Out Parameters:
// - ray_t: ray parameter at the intersection point
// - euv: primitive uv ( {0,0} for points )
//
// Returns:
// - whether the intersection occurred
//
// Iplementation Notes:
// - out parameters and only writtent o if an intersection occurs
// - algorithm finds the closest point on the ray segment to the point and
//    test their distance with the point radius
// - based on http://geomalgorithms.com/a02-lines.html.
//
static inline bool yb__intersect_point(const ym_ray3f& ray, const ym_vec3f& p,
                                       float r, float* ray_t, ym_vec2f* euv) {
    // find parameter for line-point minimum distance
    ym_vec3f w = p - ray.o;
    float t = ym_dot(w, ray.d) / ym_dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    ym_vec3f rp = ray.eval(t);
    ym_vec3f prp = p - rp;
    if (ym_dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = ym_vec2f{0, 0};

    return true;
}

//
// Intersect a ray with a line
//
// Parameters:
// - ray: ray origin and direction, parameter min, max range
// - v0, v1: line segment points
// - r0, r1: line segment radia
//
// Out Parameters:
// - ray_t: ray parameter at the intersection point
// - euv: euv[0] is the line parameter at the intersection ( euv[1] is zero )
//
// Returns:
// - whether the intersection occurred
//
// Notes:
// - out parameters and only writtent o if an intersection occurs
// - algorithm find the closest points on line and ray segment and test
//   their distance with the line radius at that location
// - based on http://geomalgorithms.com/a05-intersect-1.html
// - based on http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
//
static inline bool yb__intersect_line(const ym_ray3f& ray, const ym_vec3f& v0,
                                      const ym_vec3f& v1, float r0, float r1,
                                      float* ray_t, ym_vec2f* euv) {
    // setup intersection params
    ym_vec3f u = ray.d;
    ym_vec3f v = v1 - v0;
    ym_vec3f w = ray.o - v0;

    // compute values to solve a linear system
    float a = ym_dot(u, u);
    float b = ym_dot(u, v);
    float c = ym_dot(v, v);
    float d = ym_dot(u, w);
    float e = ym_dot(v, w);
    float det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute parameters on both ray and segment
    float t = (b * e - c * d) / det;
    float s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // clamp segment param to segment corners
    s = ym_clamp(s, 0, 1);

    // compute segment-segment distance on the closest points
    ym_vec3f p0 = yb__eval_ray(ray.o, ray.d, t);
    ym_vec3f p1 = yb__eval_ray(v0, v1 - v0, s);
    ym_vec3f p01 = p0 - p1;

    // check with the line radius at the same point
    float r = r0 * (1 - s) + r1 * s;
    if (ym_dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = ym_vec2f{s, 0};

    return true;
}

//
// Intersect a ray with a triangle
//
// Parameters:
// - ray: ray origin and direction, parameter min, max range
// - v0, v1, v2: triangle vertices
//
// Out Parameters:
// - ray_t: ray parameter at the intersection point
// - euv: baricentric coordinates of the intersection
//
// Returns:
// - whether the intersection occurred
//
// Notes:
// - out parameters and only writtent o if an intersection occurs
// - algorithm based on Muller-Trombone intersection test
//
static inline bool yb__intersect_triangle(const ym_ray3f& ray,
                                          const ym_vec3f& v0,
                                          const ym_vec3f& v1,
                                          const ym_vec3f& v2, float* ray_t,
                                          ym_vec2f* euv) {
    // compute triangle edges
    ym_vec3f edge1 = v1 - v0;
    ym_vec3f edge2 = v2 - v0;

    // compute determinant to solve a linear system
    ym_vec3f pvec = ym_cross(ray.d, edge2);
    float det = ym_dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    float inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    ym_vec3f tvec = ray.o - v0;
    float u = ym_dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    ym_vec3f qvec = ym_cross(tvec, edge1);
    float v = ym_dot(ray.d, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0) return false;

    // compute and check ray parameter
    float t = ym_dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = ym_vec2f{u, v};

    return true;
}

//
// Intersect a ray with a axis-aligned bounding box
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
// - bbox_min, bbox_max: bounding box min/max bounds
//
// Returns:
// - whether the intersection occurred
//
static inline bool yb__intersect_check_bbox(const ym_ray3f& ray,
                                            const ym_range3f& bbox) {
    // set up convenient pointers for looping over axes
    const float* _ray_o = &ray.o.x;
    const float* _ray_d = &ray.d.x;
    const float* _bbox_min = &bbox.min.x;
    const float* _bbox_max = &bbox.max.x;

    float tmin = ray.tmin, tmax = ray.tmax;

    // for each axis, clip intersection against the bounding planes
    for (int i = 0; i < 3; i++) {
        // determine intersection ranges
        float invd = 1.0f / _ray_d[i];
        float t0 = (_bbox_min[i] - _ray_o[i]) * invd;
        float t1 = (_bbox_max[i] - _ray_o[i]) * invd;
        // flip based on range directions
        if (invd < 0.0f) {
            float a = t0;
            t0 = t1;
            t1 = a;
        }
        // clip intersection
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        // if intersection is empty, exit
        if (tmin > tmax) return false;
    }

    // passed all planes, then intersection occurred
    return true;
}

//
// Min/max used in BVH traversal. Copied here since the traversal code relies
// on the specific behaviour wrt NaNs.
//
template <typename T>
static inline const T& yb__safemin(const T& a, const T& b) {
    return (a < b) ? a : b;
}
template <typename T>
static inline const T& yb__safemax(const T& a, const T& b) {
    return (a > b) ? a : b;
}

//
// Intersect a ray with a axis-aligned bounding box
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
// - ray_dinv: ray inverse direction
// - ray_dsign: ray direction sign
// - bbox_min, bbox_max: bounding box min/max bounds
//
// Returns:
// - whether the intersection occurred
//
// Implementation Notes:
// - based on "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
//
static inline bool yb__intersect_check_bbox(const ym_ray3f& ray,
                                            const ym_vec3f& ray_dinv,
                                            const ym_vec3i& ray_dsign,
                                            const ym_range3f& bbox) {
    float txmin = (bbox[ray_dsign[0]].x - ray.o.x) * ray_dinv.x;
    float txmax = (bbox[1 - ray_dsign[0]].x - ray.o.x) * ray_dinv.x;
    float tymin = (bbox[ray_dsign[1]].y - ray.o.y) * ray_dinv.y;
    float tymax = (bbox[1 - ray_dsign[1]].y - ray.o.y) * ray_dinv.y;
    float tzmin = (bbox[ray_dsign[2]].z - ray.o.z) * ray_dinv.z;
    float tzmax = (bbox[1 - ray_dsign[2]].z - ray.o.z) * ray_dinv.z;
    float tmin =
        yb__safemax(tzmin, yb__safemax(tymin, yb__safemax(txmin, ray.tmin)));
    float tmax =
        yb__safemin(tzmax, yb__safemin(tymax, yb__safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// -----------------------------------------------------------------------------
// ELEMENT-WISE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------

// TODO: documentation
static inline bool yb__distance_point(const ym_vec3f& pos, float dist_max,
                                      const ym_vec3f& p, float r, float* dist,
                                      ym_vec2f* euv) {
    float d2 = ym_distsqr(pos, p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    *dist = sqrtf(d2);
    *euv = ym_vec2f{0, 0};
    return true;
}

// TODO: documentation
static inline float yb__closestuv_line(const ym_vec3f& pos, const ym_vec3f& v0,
                                       const ym_vec3f& v1) {
    ym_vec3f ab = v1 - v0;
    float d = ym_dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    float u = ym_dot(pos - v0, ab) / d;
    u = ym_clamp(u, 0, 1);
    return u;
}

// TODO: documentation
static inline bool yb__distance_line(const ym_vec3f& pos, float dist_max,
                                     const ym_vec3f& v0, const ym_vec3f& v1,
                                     float r0, float r1, float* dist,
                                     ym_vec2f* euv) {
    float u = yb__closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    ym_vec3f p = ym_lerp(v0, v1, u);
    float r = ym_lerp(r0, r1, u);
    float d2 = ym_distsqr(pos, p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    *dist = sqrtf(d2);
    *euv = ym_vec2f{u, 0};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
static inline ym_vec2f yb__closestuv_triangle(const ym_vec3f& pos,
                                              const ym_vec3f& v0,
                                              const ym_vec3f& v1,
                                              const ym_vec3f& v2) {
    ym_vec3f ab = v1 - v0;
    ym_vec3f ac = v2 - v0;
    ym_vec3f ap = pos - v0;

    float d1 = ym_dot(ab, ap);
    float d2 = ym_dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return ym_vec2f{0, 0};

    ym_vec3f bp = pos - v1;
    float d3 = ym_dot(ab, bp);
    float d4 = ym_dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return ym_vec2f{1, 0};

    float vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return ym_vec2f{d1 / (d1 - d3), 0};

    ym_vec3f cp = pos - v2;
    float d5 = ym_dot(ab, cp);
    float d6 = ym_dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return ym_vec2f{0, 1};

    float vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) return ym_vec2f{0, d2 / (d2 - d6)};

    float va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return ym_vec2f{1 - w, w};
    }

    // face case
    float denom = 1 / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return ym_vec2f{v, w};
}

// TODO: documentation
static inline bool yb__distance_triangle(const ym_vec3f& pos, float dist_max,
                                         const ym_vec3f& v0, const ym_vec3f& v1,
                                         const ym_vec3f& v2, float r0, float r1,
                                         float r2, float* dist, ym_vec2f* euv) {
    ym_vec2f uv = yb__closestuv_triangle(pos, v0, v1, v2);
    ym_vec3f p = ym_blerp(v0, v1, v2, uv.x, uv.y);
    float r = ym_blerp(r0, r1, r2, uv.x, uv.y);
    float dd = ym_distsqr(p, pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    *dist = sqrtf(dd);
    *euv = uv;
    return true;
}

// TODO: documentation
static inline bool yb__distance_check_bbox(const ym_vec3f pos, float dist_max,
                                           const ym_vec3f bbox_min,
                                           const ym_vec3f bbox_max) {
    // set up convenient pointers for looping over axes
    const float* _pos = &pos.x;
    const float* _bbox_min = &bbox_min.x;
    const float* _bbox_max = &bbox_max.x;

    // computing distance
    float dd = 0.0f;
    // For each axis count any excess distance outside box extents
    for (int i = 0; i < 3; i++) {
        float v = _pos[i];
        if (v < _bbox_min[i]) dd += (_bbox_min[i] - v) * (_bbox_min[i] - v);
        if (v > _bbox_max[i]) dd += (v - _bbox_max[i]) * (v - _bbox_max[i]);
    }

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
static inline bool yb__overlap_bbox(const ym_vec3f bbox1_min,
                                    const ym_vec3f bbox1_max,
                                    const ym_vec3f bbox2_min,
                                    const ym_vec3f bbox2_max) {
    if (bbox1_max.x < bbox2_min.x || bbox1_min.x > bbox2_max.x) return false;
    if (bbox1_max.y < bbox2_min.y || bbox1_min.y > bbox2_max.y) return false;
    if (bbox1_max.z < bbox2_min.z || bbox1_min.z > bbox2_max.z) return false;
    return true;
}

// -----------------------------------------------------------------------------
// BVH DATA STRUCTURE
// -----------------------------------------------------------------------------

// number of primitives to avoid splitting on
#define YB__BVH_MINPRIMS 4

//
// Types of primitives contained in the entire bvh or its nodes.
//
enum {
    yb__ntype_internal = 0,  // node is internal (not leaf)
    yb__ntype_point = 1,     // points
    yb__ntype_line = 2,      // lines
    yb__ntype_triangle = 3,  // triangles
    yb__ntype_bvh = 9,       // pointers to other bvhs
};

//
// BVH tree node containing its bounds, indices to the BVH arrays of either
// sorted primitives or internal nodes, whether its a leaf or an internal node,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See yb_bvh for more details.
//
// Implemenetation Notes:
// - Padded to 32 bytes for cache fiendly access
//
struct yb__bvhn {
    ym_range3f bbox;  // bounding box
    uint32_t start;   // index to the first sorted primitive/node
    uint16_t count;   // number of primitives/nodes
    uint8_t isleaf;   // whether it is a leaf
    uint8_t axis;     // plit axis
};

//
// BVH tree, stored as a node array. The tree structure is encoded using array
// indices instead of pointers, both for speed but also to simplify code.
// BVH nodes indices refer to either the node array, for internal nodes,
// or a primitive array, for leaf nodes. BVH trees may contain only one type
// of geometric primitive, like points, lines, triangle or shape other BVHs.
// We handle multiple primitive types and transformed primitices by building
// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
// BVHs, shape BVHs, each of which of a uniform primitive type.
//
// Implementation Notes:
// - we handle three types of BVHs, shape, shape and scene, using a union
// - shape BVHs contain sorted geometric primitives for fast intersection
// - shape BVHs contain pointers to external geometric data, to save memory
// - scene BVHs contain other BVHs to build the two level hiearchy
// - all BVH types shares the same type of internal nodes
//
struct yb__bvh {
    // bvh data
    ym_vector<yb__bvhn> nodes;   // sorted array of internal nodes
    ym_vector<int> sorted_prim;  // sorted elements
};

//
// Shape Data
//
struct yb_shape {
    // shape data
    int nelems = 0;  // shape elements
    int etype = 0;   // shape element type
    union {
        const int* elem = nullptr;  // shape elements
        const int* point;
        const ym_vec2i* line;
        const ym_vec3i* triangle;
    };
    const ym_vec3f* pos = nullptr;  // vertex pos
    const float* radius = nullptr;  // vertex radius
    // bvh data
    yb__bvh* bvh = nullptr;  // bvh
};

//
// Scene Data
//
struct yb_scene {
    // scene data
    ym_vector<yb_shape> shapes;  // shapes
    ym_vector<ym_affine3f> xforms;
    ym_vector<ym_affine3f> inv_xforms;
    // bvh data
    yb__bvh* bvh = nullptr;
};

// -----------------------------------------------------------------------------
// BVH BUILD FUNCTIONS
// -----------------------------------------------------------------------------

//
// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
//
struct yb__bound_prim {
    ym_range3f bbox;       // bounding box
    ym_vec3f center;       // bounding box center (for faster sort)
    int pid;               // primitive id
    float sah_cost_left;   // buffer for sah heuristic costs
    float sah_cost_right;  // buffer for sah heuristic costs
};

//
// Comparison function for each axis
//
struct yb__bound_prim_comp {
    int axis;
    float middle;
    yb__bound_prim_comp(int a, float m = 0) : axis(a), middle(m) {}
    bool operator()(const yb__bound_prim& a, const yb__bound_prim& b) const {
        return a.center[axis] < b.center[axis];
    }
    bool operator()(const yb__bound_prim& a) const {
        return a.center[axis] < middle;
    }
};

//
// Given an array sorted_prim of primitives to split between the elements
// start and end, determines the split axis axis, split primitive index mid
// based on the heuristic heuristic. Supports balanced tree (equalnum) and
// Surface-Area Heuristic.
//
static inline bool yb__partition_prims(yb__bound_prim* sorted_prim, int start,
                                       int end, int* axis, int* mid,
                                       int heuristic) {
    const float __box_eps = 1e-12f;
#define __bbox_area(r)                                                         \
    (2 * ((r.x + __box_eps) * (r.y + __box_eps) +                              \
          (r.y + __box_eps) * (r.z + __box_eps) +                              \
          (r.z + __box_eps) * (r.x + __box_eps)))

    // init to default values
    *axis = 0;
    *mid = (start + end) / 2;

    // compute primintive bounds and size
    ym_range3f centroid_bbox = ym_invalid_range3f;
    for (int i = start; i < end; i++) centroid_bbox += sorted_prim[i].center;
    ym_vec3f centroid_size = ym_rsize(centroid_bbox);

    // check if it is not possible to split
    if (centroid_size == ym_zero3f) return false;

    // split along largest
    int largest_axis = ym_max_element(centroid_size);

    // check heuristic
    switch (heuristic) {
        // balanced tree split: find the largest axis of the bounding box
        // and split along this one right in the middle
        case yb_htype_equalnum: {
            *axis = largest_axis;
            *mid = (start + end) / 2;
            std::nth_element(sorted_prim + start, sorted_prim + (*mid),
                             sorted_prim + end,
                             yb__bound_prim_comp(largest_axis));
        } break;
        // split the space in the middle along the largest axis
        case yb_htype_default:
        case yb_htype_equalsize: {
            *axis = largest_axis;
            *mid = (int)(std::partition(
                             sorted_prim + start, sorted_prim + end,
                             yb__bound_prim_comp(
                                 largest_axis,
                                 ym_rcenter(centroid_bbox)[largest_axis])) -
                         sorted_prim);
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case yb_htype_sah: {
            float min_cost = HUGE_VALF;
            int count = end - start;
            for (int a = 0; a < 3; a++) {
                std::sort(sorted_prim + start, sorted_prim + end,
                          yb__bound_prim_comp(a));
                ym_range3f sbbox = ym_invalid_range3f;
                // to avoid an O(n^2) computation, use sweaps to compute the
                // cost,
                // first smallest to largest, then largest to smallest
                for (int i = 0; i < count; i++) {
                    sbbox += sorted_prim[start + i].bbox;
                    ym_vec3f sbbox_size = ym_rsize(sbbox);
                    sorted_prim[start + i].sah_cost_left =
                        __bbox_area(sbbox_size);
                    sorted_prim[start + i].sah_cost_left *= i + 1;
                }
                // the other sweep
                sbbox = ym_invalid_range3f;
                for (int i = 0; i < count; i++) {
                    sbbox += sorted_prim[end - 1 - i].bbox;
                    ym_vec3f sbbox_size = ym_rsize(sbbox);
                    sorted_prim[end - 1 - i].sah_cost_right =
                        __bbox_area(sbbox_size);
                    sorted_prim[end - 1 - i].sah_cost_right *= i + 1;
                }
                // find the minimum cost
                for (int i = start + 2; i <= end - 2; i++) {
                    float cost = sorted_prim[i - 1].sah_cost_left +
                                 sorted_prim[i].sah_cost_right;
                    if (min_cost > cost) {
                        min_cost = cost;
                        *axis = a;
                        *mid = i;
                    }
                }
            }
            std::nth_element(sorted_prim + start, sorted_prim + (*mid),
                             sorted_prim + end, yb__bound_prim_comp(*axis));
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case yb_htype_binned_sah: {
            // allocate bins
            const int nbins = 16;
            ym_range3f bins_bbox[nbins];
            int bins_count[nbins];
            for (int b = 0; b < nbins; b++) {
                bins_bbox[b] = ym_invalid_range3f;
                // bins_bbox[b] = centroid_bbox;
                // bins_bbox[b].min[largest_axis] += b *
                // centroid_size[largest_axis] / nbins;
                // bins_bbox[b].max[largest_axis] -= (nbins - 1 - b) *
                // centroid_size[largest_axis] / nbins;
                bins_count[b] = 0;
            }
            for (int i = start; i < end; i++) {
                int b = (int)(nbins * (sorted_prim[i].center[largest_axis] -
                                       centroid_bbox.min[largest_axis]) /
                              centroid_size[largest_axis]);
                b = ym_clamp(b, 0, nbins - 1);
                bins_count[b] += 1;
                bins_bbox[b] += sorted_prim[i].bbox;
            }
            float min_cost = HUGE_VALF;
            int bin_idx = -1;
            for (int b = 1; b < nbins; b++) {
                if (!bins_count[b]) continue;
                ym_range3f left_bbox = ym_invalid_range3f,
                           right_bbox = ym_invalid_range3f;
                int left_count = 0, right_count = 0;
                for (int j = 0; j < b; j++) {
                    if (!bins_count[j]) continue;
                    left_count += bins_count[j];
                    left_bbox += bins_bbox[j];
                }
                for (int j = b; j < nbins; j++) {
                    if (!bins_count[j]) continue;
                    right_count += bins_count[j];
                    right_bbox += bins_bbox[j];
                }
                ym_vec3f left_bbox_size = ym_rsize(left_bbox);
                ym_vec3f right_bbox_size = ym_rsize(right_bbox);
                float cost = __bbox_area(left_bbox_size) * left_count +
                             __bbox_area(right_bbox_size) * right_count;
                if (min_cost > cost) {
                    min_cost = cost;
                    bin_idx = b;
                }
            }
            assert(bin_idx >= 0);
            *axis = largest_axis;
            *mid = start;
            for (int b = 0; b < bin_idx; b++) {
                *mid += bins_count[b];
            }
            assert(*axis >= 0 && *mid > 0);
            std::nth_element(sorted_prim + start, sorted_prim + (*mid),
                             sorted_prim + end,
                             yb__bound_prim_comp(largest_axis));
        } break;
        default: assert(false); break;
    }
    assert(*axis >= 0 && *mid > 0);
    assert(*mid > start && *mid < end);
    return true;
}

//
// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
//
static inline void yb__make_node(yb__bvhn* node, int* nnodes, yb__bvhn* nodes,
                                 yb__bound_prim* sorted_prims, int start,
                                 int end, int heuristic) {
    // compute node bounds
    node->bbox = ym_invalid_range3f;
    for (int i = start; i < end; i++) {
        node->bbox += sorted_prims[i].bbox;
    }

    // decide whether to create a leaf
    if (end - start <= YB__BVH_MINPRIMS) {
        // makes a leaf node
        node->isleaf = true;
        node->start = start;
        node->count = end - start;
    } else {
        // choose the split axis and position
        int axis, mid;
        bool split = yb__partition_prims(sorted_prims, start, end, &axis, &mid,
                                         heuristic);
        if (!split) {
            // we failed to split for some reasons
            node->isleaf = true;
            node->start = start;
            node->count = end - start;
        } else {
            assert(mid > start && mid < end);
            // makes an internal node
            node->isleaf = false;
            // perform the splits by preallocating the child nodes and recurring
            node->axis = axis;
            node->start = *nnodes;
            node->count = 2;
            *nnodes += 2;
            // build child nodes
            yb__make_node(nodes + node->start, nnodes, nodes, sorted_prims,
                          start, mid, heuristic);
            yb__make_node(nodes + node->start + 1, nnodes, nodes, sorted_prims,
                          mid, end, heuristic);
        }
    }
}

//
// Makes a scene BVH. Public function whose interface is described above.
//
YGL_API yb_scene* yb_init_scene(int nshapes) {
    // allocates the BVH and set its primitive type
    yb_scene* scene = new yb_scene();
    scene->shapes.assign(nshapes, yb_shape());
    scene->xforms.assign(nshapes, ym_identity_affine3f);
    scene->inv_xforms.assign(nshapes, ym_identity_affine3f);
    scene->bvh = nullptr;
    return scene;
}

//
// Sets shape data for scene.
//
YGL_API void yb_set_shape(yb_scene* scene, int sid, const ym_affine3f& xform,
                          int nelems, const int* elem, int etype, int nverts,
                          const ym_vec3f* pos, const float* radius) {
    scene->shapes[sid].nelems = nelems;
    scene->shapes[sid].etype = etype;
    scene->shapes[sid].elem = elem;
    scene->shapes[sid].pos = pos;
    scene->shapes[sid].radius = radius;
    scene->xforms[sid] = xform;
    if (xform == ym_identity_affine3f) {
        scene->inv_xforms[sid] = xform;
    } else {
        scene->inv_xforms[sid] = ym_inverse(xform);
    }
    if (scene->shapes[sid].bvh) delete scene->shapes[sid].bvh;
    scene->shapes[sid].bvh = nullptr;
}

//
// Build a BVH from a set of primitives.
//
YGL_API yb__bvh* yb__build_bvh(ym_vector<yb__bound_prim>& bound_prims,
                               int heuristic) {
    // create bvh
    yb__bvh* bvh = new yb__bvh();

    // allocate nodes (over-allocate now then shrink)
    int nnodes = 0;
    bvh->nodes.resize(bound_prims.size() * 2);

    // start recursive splitting
    nnodes = 1;
    yb__make_node(bvh->nodes.data(), &nnodes, bvh->nodes.data(),
                  bound_prims.data(), 0, (int)bound_prims.size(), heuristic);

    // shrink back
    bvh->nodes.resize(nnodes);
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(bound_prims.size());
    for (int i = 0; i < bound_prims.size(); i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }
    bvh->sorted_prim.shrink_to_fit();

    // done
    return bvh;
}

//
// Build a shape BVH by creating a bound_prims array for all prims and
// passing it to the bvh building code.
//
YGL_API void yb__build_shape_bvh(yb_shape* shape, int heuristic) {
    // create bounded primitivrs used in BVH build
    ym_vector<yb__bound_prim> bound_prims =
        ym_vector<yb__bound_prim>(shape->nelems);
    assert(shape->pos && shape->elem);
    // compute bounds for each primitive type
    switch (shape->etype) {
        case yb_etype_point: {
            // point bounds are computed as small spheres
            for (int i = 0; i < shape->nelems; i++) {
                const int& f = shape->point[i];
                yb__bound_prim* prim = &bound_prims[i];
                prim->pid = i;
                prim->bbox = ym_invalid_range3f;
                float r0 = (shape->radius) ? shape->radius[f] : 0;
                prim->bbox += shape->pos[f] - ym_vec3f(r0);
                prim->bbox += shape->pos[f] + ym_vec3f(r0);
                prim->center = ym_rcenter(prim->bbox);
            }
        } break;
        case yb_etype_line: {
            // line bounds are computed as thick rods
            assert(shape->radius);
            for (int i = 0; i < shape->nelems; i++) {
                const ym_vec2i& f = shape->line[i];
                yb__bound_prim* prim = &bound_prims[i];
                prim->pid = i;
                prim->bbox = ym_invalid_range3f;
                float r0 = (shape->radius) ? shape->radius[f.x] : 0;
                float r1 = (shape->radius) ? shape->radius[f.y] : 0;
                prim->bbox += shape->pos[f.x] - ym_vec3f(r0);
                prim->bbox += shape->pos[f.x] + ym_vec3f(r0);
                prim->bbox += shape->pos[f.y] - ym_vec3f(r1);
                prim->bbox += shape->pos[f.y] + ym_vec3f(r1);
                prim->center = ym_rcenter(prim->bbox);
            }
        } break;
        case yb_etype_triangle: {
            // triangle bounds are computed by including their vertices
            for (int i = 0; i < shape->nelems; i++) {
                const ym_vec3i& f = shape->triangle[i];
                yb__bound_prim* prim = &bound_prims[i];
                prim->pid = i;
                prim->bbox = ym_invalid_range3f;
                prim->bbox += shape->pos[f.x];
                prim->bbox += shape->pos[f.y];
                prim->bbox += shape->pos[f.z];
                prim->center = ym_rcenter(prim->bbox);
            }
        } break;
        default: assert(false);
    }

    // build bvh
    shape->bvh = yb__build_bvh(bound_prims, heuristic);
}

//
// Build a scene BVH by wrapping each shape into a bounding box, transforming it
// and creating a bvh of these.
//
YGL_API yb__bvh* yb__build_scene_bvh(const ym_vector<ym_range3f>& bboxes,
                                     const ym_vector<ym_affine3f>& xforms,
                                     int heuristic) {
    // compute bounds for all transformed shape bvhs
    ym_vector<yb__bound_prim> bound_prims =
        ym_vector<yb__bound_prim>(bboxes.size());
    for (int i = 0; i < bboxes.size(); i++) {
        yb__bound_prim* prim = &bound_prims[i];
        ym_affine3f xform = xforms[i];
        prim->pid = i;
        prim->bbox = bboxes[i];
        // apply transformations if needed (estimate trasform BVH bounds as
        // corners bounds, overestimate, but best effort for now)
        ym_vec3f corners[8] = {
            {prim->bbox.min.x, prim->bbox.min.y, prim->bbox.min.z},
            {prim->bbox.min.x, prim->bbox.min.y, prim->bbox.max.z},
            {prim->bbox.min.x, prim->bbox.max.y, prim->bbox.min.z},
            {prim->bbox.min.x, prim->bbox.max.y, prim->bbox.max.z},
            {prim->bbox.max.x, prim->bbox.min.y, prim->bbox.min.z},
            {prim->bbox.max.x, prim->bbox.min.y, prim->bbox.max.z},
            {prim->bbox.max.x, prim->bbox.max.y, prim->bbox.min.z},
            {prim->bbox.max.x, prim->bbox.max.y, prim->bbox.max.z},
        };
        prim->bbox = ym_invalid_range3f;
        for (int j = 0; j < 8; j++)
            prim->bbox += ym_transform_point(xform, corners[j]);
        prim->center = ym_rcenter(prim->bbox);
    }

    // make bvh
    return yb__build_bvh(bound_prims, heuristic);
}

//
// Build a scene BVH. Public function whose interface is described above.
//
YGL_API void yb_build_bvh(yb_scene* scene, int heuristic) {
    // recursively build shape bvhs
    ym_vector<ym_range3f> shape_bounds(scene->shapes.size());
    for (int sid = 0; sid < scene->shapes.size(); sid++) {
        yb__build_shape_bvh(&scene->shapes[sid], heuristic);
        shape_bounds[sid] = scene->shapes[sid].bvh->nodes[0].bbox;
    }

    // build scene bvh
    scene->bvh = yb__build_scene_bvh(shape_bounds, scene->xforms, heuristic);
}

//
// Recursively recomputes the node bounds for a scene bvh
//
static inline void yb__recompute_scene_bounds(yb_scene* scene, int nodeid) {
    yb__bvh* bvh = scene->bvh;
    yb__bvhn* node = &bvh->nodes[nodeid];
    if (node->isleaf) {
        node->bbox = ym_invalid_range3f;
        for (int i = 0; i < node->count; i++) {
            int idx = bvh->sorted_prim[node->start + i];
            ym_range3f bbox = scene->shapes[idx].bvh->nodes[0].bbox;
            ym_affine3f xform = scene->xforms[idx];
            ym_vec3f corners[8] = {
                {bbox.min.x, bbox.min.y, bbox.min.z},
                {bbox.min.x, bbox.min.y, bbox.max.z},
                {bbox.min.x, bbox.max.y, bbox.min.z},
                {bbox.min.x, bbox.max.y, bbox.max.z},
                {bbox.max.x, bbox.min.y, bbox.min.z},
                {bbox.max.x, bbox.min.y, bbox.max.z},
                {bbox.max.x, bbox.max.y, bbox.min.z},
                {bbox.max.x, bbox.max.y, bbox.max.z},
            };
            for (int j = 0; j < 8; j++)
                node->bbox += ym_transform_point(xform, corners[j]);
        }
    } else {
        node->bbox = ym_invalid_range3f;
        for (int i = 0; i < node->count; i++) {
            yb__recompute_scene_bounds(scene, node->start + i);
            node->bbox += bvh->nodes[node->start + i].bbox;
        }
    }
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
YGL_API void yb_refit_scene_bvh(yb_scene* scene, const ym_affine3f* xforms) {
    // updated xforms
    for (int i = 0; i < scene->shapes.size(); i++) {
        scene->xforms[i] = xforms[i];
        if (xforms[i] == ym_identity_affine3f) {
            scene->inv_xforms[i] = ym_identity_affine3f;
        } else {
            scene->inv_xforms[i] = ym_inverse(xforms[i]);
        }
    }

    // recompute bvh bounds
    yb__recompute_scene_bounds(scene, 0);
}

//
// Free BVH data.
//
YGL_API void yb_free_scene(yb_scene* scene) { delete scene; }

// -----------------------------------------------------------------------------
// BVH INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Logging global variables
//
#ifdef YGL_BVH_LOG_RAYS
int yb__log_nrays = 0;
int yb__log_nbbox_inters = 0;
int yb__log_npoint_inters = 0;
int yb__log_nline_inters = 0;
int yb__log_ntriangle_inters = 0;

YGL_API void yb_get_ray_log(int* nrays, int* nbbox_inters, int* npoint_inters,
                            int* nline_inters, int* ntriangle_inters) {
    *nrays = yb__log_nrays;
    *nbbox_inters = yb__log_nbbox_inters;
    *npoint_inters = yb__log_npoint_inters;
    *nline_inters = yb__log_nline_inters;
    *ntriangle_inters = yb__log_ntriangle_inters;
}
#endif

//
// Intersect ray with a bvh. Similar to the generic public function whose
// interface is described above. See yb_intersect_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of recursive
// calls; this follows general conventions and stragely makes the code shorter
// - The walk is simplified for first hit by obeserving that if we update
// the ray_tmax liminit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequence farthest iterations will be
// rejected in the tmax tests
//
// TODO: update doc
//
static inline bool yb__intersect_bvh(const yb_scene* scene, int shape_id,
                                     bool early_exit, const ym_ray3f& ray_,
                                     float* ray_t, int* sid, int* eid,
                                     ym_vec2f* euv) {
    // node stack
    int node_stack[64];
    int node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    bool hit = false;

    // init shape and bvh
    const yb_shape* shape =
        (shape_id >= 0) ? &scene->shapes[shape_id] : nullptr;
    const yb__bvh* bvh = (shape) ? shape->bvh : scene->bvh;

    // init ray
    ym_ray3f ray = ray_;
    if (shape) {
        ray.o = ym_transform_point(scene->inv_xforms[shape_id], ray.o);
        ray.d = ym_transform_vector(scene->inv_xforms[shape_id], ray.d);
    }

#ifndef YGL_BVH_LEGACY_INTERSECTION
    // prepare ray for fast queries
    ym_vec3f ray_dinv = {1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
    ym_vec3i ray_dsign = {(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
                          (ray_dinv.z < 0) ? 1 : 0};
#endif

#ifdef YGL_BVH_LOG_RAYS
    if (!shape) yb__log_nrays++;
#endif

    // walking stack
    while (node_cur && !(early_exit && hit)) {
        // exit if needed
        if (early_exit && hit) return hit;

        // grab node
        const yb__bvhn* node = &bvh->nodes[node_stack[--node_cur]];

#ifdef YGL_BVH_LOG_RAYS
        yb__log_nbbox_inters++;
#endif

// intersect bbox
#ifndef YGL_BVH_LEGACY_INTERSECTION
        if (!yb__intersect_check_bbox(ray, ray_dinv, ray_dsign, node->bbox))
            continue;
#else
        if (!yb__intersect_check_bbox(ray, node->bbox)) continue;
#endif

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node->isleaf) {
// for internal nodes, attempts to proceed along the
// split axis from smallest to largest nodes
#ifndef YGL_BVH_LEGACY_INTERSECTION
            bool reverse = ray_dsign[node->axis];
#else
            bool reverse = (&ray.d.x)[node->axis] < 0;
#endif
            if (reverse) {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            } else {
                for (int i = node->count - 1; i >= 0; i--) {
                    int idx = node->start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            }
        } else if (shape) {
            if (shape->etype == yb_etype_triangle) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const ym_vec3i& f = shape->triangle[idx];
                    if (yb__intersect_triangle(ray, shape->pos[f.x],
                                               shape->pos[f.y], shape->pos[f.z],
                                               &ray.tmax, euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
                }
#ifdef YGL_BVH_LOG_RAYS
                yb__log_ntriangle_inters += node->count;
#endif

            } else if (shape->etype == yb_etype_line) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const ym_vec2i& f = shape->line[idx];
                    if (yb__intersect_line(ray, shape->pos[f.x],
                                           shape->pos[f.y], shape->radius[f.x],
                                           shape->radius[f.y], &ray.tmax,
                                           euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
                }
#ifdef YGL_BVH_LOG_RAYS
                yb__log_nline_inters += node->count;
#endif
            } else if (shape->etype == yb_etype_point) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const int& f = shape->point[idx];
                    if (yb__intersect_point(ray, shape->pos[f],
                                            shape->radius[f], &ray.tmax, euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
#ifdef YGL_BVH_LOG_RAYS
                    yb__log_npoint_inters += node->count;
#endif
                }
            } else {
                assert(false);
            }
        } else {
            for (int i = 0; i < node->count; i++) {
                int idx = bvh->sorted_prim[node->start + i];
                if (yb__intersect_bvh(scene, idx, early_exit, ray, &ray.tmax,
                                      sid, eid, euv)) {
                    hit = true;
                }
            }
        }
    }

    // update ray distance for outside
    if (hit) *ray_t = ray.tmax;

    return hit;
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool yb_intersect_first(const yb_scene* scene, const ym_ray3f& ray,
                                float* ray_t, int* sid, int* eid, ym_vec2f* euv,
                                int req_shape) {
    return yb__intersect_bvh(scene, req_shape, false, ray, ray_t, sid, eid,
                             euv);
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool yb_intersect_any(const yb_scene* scene, const ym_ray3f& ray,
                              int req_shape) {
    int sid, eid;
    float ray_t;
    ym_vec2f euv;
    return yb__intersect_bvh(scene, req_shape, true, ray, &ray_t, &sid, &eid,
                             &euv);
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP
// -----------------------------------------------------------------------------

//
// Finds the closest element with a bvh. Similar to the generic public function
// whose
// interface is described above. See yb_nightbor_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of recursive
// calls; this follows general conventions and stragely makes the code shorter
// - The walk is simplified for first hit by obeserving that if we update
// the dist_max limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequent farthest iterations will be
// rejected in the tmax tests
//
static inline bool yb__overlap_bvh(const yb_scene* scene, int shape_id,
                                   bool early_exit, const ym_vec3f& pt_,
                                   float dist_max, float* dist, int* sid,
                                   int* eid, ym_vec2f* euv) {
    // node stack
    int node_stack[64];
    int node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    bool hit = false;

    // init shape and bvh
    const yb_shape* shape =
        (shape_id >= 0) ? &scene->shapes[shape_id] : nullptr;
    const yb__bvh* bvh = (shape) ? shape->bvh : scene->bvh;

    // init ray
    ym_vec3f pt = pt_;
    if (shape) {
        pt = ym_transform_point(scene->inv_xforms[shape_id], pt);
    }

    // walking stack
    while (node_cur && !(early_exit && hit)) {
        // exit if needed
        if (early_exit && hit) return hit;

        // grab node
        const yb__bvhn* node = &bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!yb__distance_check_bbox(pt, dist_max, node->bbox.min,
                                     node->bbox.max))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node->isleaf) {
            // internal node
            for (int i = 0; i < node->count; i++) {
                int idx = node->start + i;
                node_stack[node_cur++] = idx;
                assert(node_cur < 64);
            }
        } else if (shape) {
            if (shape->etype == yb_etype_triangle) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const ym_vec3i& f = shape->triangle[idx];
                    if (yb__distance_triangle(
                            pt, dist_max, shape->pos[f.x], shape->pos[f.y],
                            shape->pos[f.z],
                            (shape->radius) ? shape->radius[f.x] : 0,
                            (shape->radius) ? shape->radius[f.y] : 0,
                            (shape->radius) ? shape->radius[f.z] : 0, &dist_max,
                            euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
                }
            } else if (shape->etype == yb_etype_line) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const ym_vec2i& f = shape->line[idx];
                    if (yb__distance_line(
                            pt, dist_max, shape->pos[f.x], shape->pos[f.y],
                            (shape->radius) ? shape->radius[f.x] : 0,
                            (shape->radius) ? shape->radius[f.y] : 0, &dist_max,
                            euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
                }
            } else if (shape->etype == yb_etype_point) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    const int& f = shape->point[idx];
                    if (yb__distance_point(pt, dist_max, shape->pos[f],
                                           (shape->radius) ? shape->radius[f]
                                                           : 0,
                                           &dist_max, euv)) {
                        hit = true;
                        *eid = idx;
                        *sid = shape_id;
                    }
                }

            } else {
                assert(false);
            }
        } else {
            for (int i = 0; i < node->count; i++) {
                int idx = bvh->sorted_prim[node->start + i];
                if (yb__overlap_bvh(scene, idx, early_exit, pt, dist_max,
                                    &dist_max, sid, eid, euv)) {
                    hit = true;
                }
            }
        }
    }

    // update ray distance for outside
    if (hit) *dist = dist_max;

    return hit;
}

//
// Find the closest element to a given point within a maximum distance.
// A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool yb_overlap_first(const yb_scene* scene, const ym_vec3f& pt,
                              float max_dist, float* dist, int* sid, int* eid,
                              ym_vec2f* euv, int req_shape) {
    return yb__overlap_bvh(scene, req_shape, false, pt, max_dist, dist, sid,
                           eid, euv);
}

//
// Find if any element overlaps given point within a maximum distance.
// A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool yb_overlap_any(const yb_scene* scene, const ym_vec3f& pt,
                            float max_dist, int req_shape) {
    float dist;
    int sid, eid;
    ym_vec2f euv;
    return yb__overlap_bvh(scene, req_shape, true, pt, max_dist, &dist, &sid,
                           &eid, &euv);
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
YGL_API int yb_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                    void* overlap_ctx,
                                    void (*overlap_cb)(void* ctx,
                                                       const ym_vec2i& v)) {
    // node stack
    int node_stack[256];
    int node_cur = 0;
    node_stack[node_cur++] = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    int hits = 0;

    // init shape and bvh
    const yb__bvh* bvh = scene->bvh;

    // walking stack
    while (node_cur) {
        // grab node
        int node1_idx = node_stack[--node_cur];
        int node2_idx = node_stack[--node_cur];
        const yb__bvhn* node1 = &bvh->nodes[node1_idx];
        const yb__bvhn* node2 = &bvh->nodes[node2_idx];

        // intersect bbox
        if (!yb__overlap_bbox(node1->bbox.min, node1->bbox.max, node2->bbox.min,
                              node2->bbox.max))
            continue;

        // check for leaves
        if (node1->isleaf && node2->isleaf) {
            // collide primitives
            for (int i1 = 0; i1 < node1->count; i1++) {
                for (int i2 = 0; i2 < node2->count; i2++) {
                    int idx1 = bvh->sorted_prim[node1->start + i1];
                    int idx2 = bvh->sorted_prim[node2->start + i2];
                    if (exclude_self && idx1 == idx2) continue;
                    const ym_range3f bbox1 = yb__transform_bbox(
                        scene->xforms[idx1],
                        scene->shapes[idx1].bvh->nodes[0].bbox);
                    const ym_range3f bbox2 = yb__transform_bbox(
                        scene->xforms[idx2],
                        scene->shapes[idx2].bvh->nodes[0].bbox);
                    if (!yb__overlap_bbox(bbox1.min, bbox1.max, bbox2.min,
                                          bbox2.max))
                        continue;
                    hits += 1;
                    overlap_cb(overlap_ctx, {idx1, idx2});
                }
            }
        } else {
            // descend
            if (node1->isleaf) {
                for (int i = 0; i < node2->count; i++) {
                    int idx = node2->start + i;
                    node_stack[node_cur++] = node1_idx;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 256);
                }
            } else {
                for (int i = 0; i < node1->count; i++) {
                    int idx = node1->start + i;
                    node_stack[node_cur++] = idx;
                    node_stack[node_cur++] = node2_idx;
                    assert(node_cur < 256);
                }
            }
        }
    }

    // done
    return hits;
}

//
// Overlap callback for ym_vector
//
static inline void yb__overlap_cb_vector(void* ctx, const ym_vec2i& v) {
    ym_vector<ym_vec2i>* vec = (ym_vector<ym_vec2i>*)ctx;
    vec->push_back(v);
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
YGL_API int yb_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                    ym_vector<ym_vec2i>* overlaps) {
    overlaps->clear();
    return yb_overlap_shape_bounds(scene, exclude_self, overlaps,
                                   yb__overlap_cb_vector);
}

// -----------------------------------------------------------------------------
// VERTEX PROPERTY INTERPOLATION
// -----------------------------------------------------------------------------

//
// Interpolate vertex properties at the intersection.
// Public function whose interface is described above.
//
YGL_API void yb_interpolate_vert(const int* elem, int etype, int eid,
                                 const ym_vec2f& euv, int vsize,
                                 const float* vert, float* v) {
    for (int c = 0; c < vsize; c++) v[c] = 0;
    switch (etype) {
        case yb_etype_point: {
            const int* f = elem + eid;
            const float* vf = vert + f[0] * vsize;
            for (int c = 0; c < vsize; c++) v[c] += vf[c];
        } break;
        case yb_etype_line: {
            const int* f = elem + eid * 2;
            const float* vf[] = {vert + f[0] * vsize, vert + f[1] * vsize};
            const float w[2] = {1 - euv[0], euv[0]};
            for (int i = 0; i < 2; i++) {
                for (int c = 0; c < vsize; c++) v[c] += w[i] * vf[i][c];
            }
        } break;
        case yb_etype_triangle: {
            const int* f = elem + eid * 3;
            const float* vf[] = {vert + f[0] * vsize, vert + f[1] * vsize,
                                 vert + f[2] * vsize};
            const float w[3] = {1 - euv[0] - euv[1], euv[0], euv[1]};
            for (int i = 0; i < 3; i++) {
                for (int c = 0; c < vsize; c++) v[c] += w[i] * vf[i][c];
            }
        } break;
        default: assert(false);
    }
}

// -----------------------------------------------------------------------------
// STATISTICS FOR DEBUGGING (probably not helpful to all)
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
YGL_API void yb__compute_bvh_stats(const yb_scene* scene, int shape_id,
                                   bool include_shapes, int depth, int* nprims,
                                   int* ninternals, int* nleaves,
                                   int* min_depth, int* max_depth) {
    // get bvh
    const yb__bvh* bvh =
        (shape_id >= 0) ? scene->shapes[shape_id].bvh : scene->bvh;

    // node stack
    ym_vec2i node_stack[128];  // node and depth
    int node_cur = 0;
    node_stack[node_cur++] = ym_vec2i{0, depth};

    // walk the stack
    while (node_cur) {
        // get node and depth
        ym_vec2i node_depth = node_stack[--node_cur];
        const yb__bvhn* node = &bvh->nodes[node_depth.x];

        // update stats
        if (!node->isleaf) {
            *ninternals += 1;
            for (int i = 0; i < node->count; i++) {
                node_stack[node_cur++] =
                    ym_vec2i(node->start + i, node_depth.y + 1);
            }
        } else if (shape_id >= 0) {
            *nleaves += 1;
            *nprims += node->count;
            *min_depth = ym_min(*min_depth, node_depth.y);
            *max_depth = ym_max(*max_depth, node_depth.y);
        } else {
            if (include_shapes) {
                for (int i = 0; i < node->count; i++) {
                    int idx = bvh->sorted_prim[node->start + i];
                    yb__compute_bvh_stats(scene, idx, true, node_depth.y + 1,
                                          nprims, ninternals, nleaves,
                                          min_depth, max_depth);
                }
            } else {
                *nleaves += 1;
                *nprims += node->count;
                *min_depth = ym_min(*min_depth, node_depth.y);
                *max_depth = ym_max(*max_depth, node_depth.y);
            }
        }
    }
}

//
// Compute BVH stats.
//
YGL_API void yb_compute_bvh_stats(const yb_scene* scene, bool include_shapes,
                                  int* nprims, int* ninternals, int* nleaves,
                                  int* min_depth, int* max_depth,
                                  int req_shape) {
    // init out variables
    *nprims = 0;
    *ninternals = 0;
    *nleaves = 0;
    *min_depth = 10000;
    *max_depth = -1;

    yb__compute_bvh_stats(scene, req_shape, include_shapes, 0, nprims,
                          ninternals, nleaves, min_depth, max_depth);
}

// -----------------------------------------------------------------------------
// C/C++ INTERFACE
// -----------------------------------------------------------------------------

//
// Create a BVH for a collection of transformed shapes (scene).
//
YGLC_API yb_scene* ybc_init_scene(int nshapes) {
    return yb_init_scene(nshapes);
}

//
// Sets the data for a given shape.
//
YGLC_API void ybc_set_shape(yb_scene* scene, int sid, const float xform[16],
                            int nelems, const int* elem, int etype, int nverts,
                            const float* pos, const float* radius) {
    return yb_set_shape(scene, sid, ym_affine3f(ym_mat4f(xform)), nelems, elem,
                        etype, nverts, (const ym_vec3f*)pos, radius);
}

//
// Builds a scene BVH.
//
YGLC_API void ybc_build_bvh(yb_scene* scene, int heuristic) {
    return yb_build_bvh(scene, heuristic);
}

//
// Clears BVH memory.
//
YGLC_API void ybc_free_scene(yb_scene* scene) { return yb_free_scene(scene); }

//
// Refit the bounds of each shape for moving objects. Use this only to avoid
// a rebuild, but note that queries are likely slow if objects move a lot.
//
YGLC_API void ybc_refit_scene_bvh(yb_scene* scene, const float* xforms[16]) {
    ym_vector<ym_affine3f> affines;
    if (xforms) {
        affines.resize(scene->shapes.size());
        for (int i = 0; i < affines.size(); i++)
            affines[i] = ym_affine3f(ym_mat4f(xforms[i]));
    }
    yb_refit_scene_bvh(scene, affines.data());
}

//
// Intersect the scene with a ray finding the closest intersection.
//
YGLC_API bool ybc_intersect_first(const yb_scene* scene, const float ray_o[3],
                                  const float ray_d[3], float ray_tmin,
                                  float ray_tmax, float* ray_t, int* sid,
                                  int* eid, float euv[2], int req_shape) {
    return yb_intersect_first(
        scene, ym_ray3f(ym_vec3f(ray_o), ym_vec3f(ray_d), ray_tmin, ray_tmax),
        ray_t, sid, eid, (ym_vec2f*)euv, req_shape);
}

//
// Intersect the scene with a ray finding any intersection.
//
YGLC_API bool ybc_intersect_any(const yb_scene* scene, const float ray_o[3],
                                const float ray_d[3], float ray_tmin,
                                float ray_tmax, int req_shape) {
    return yb_intersect_any(
        scene, ym_ray3f(ym_vec3f(ray_o), ym_vec3f(ray_d), ray_tmin, ray_tmax),
        req_shape);
}

//
// Returns a list of shape pairs that can possibly overlap by checking only they
// axis aligned bouds. This is only a conservative check useful for collision
// detection.
//
YGLC_API int ybc_overlap_shape_bounds(const yb_scene* scene, bool exclude_self,
                                      void* overlap_ctx,
                                      void (*overlap_cb)(void* ctx,
                                                         const int ij[2])) {
    return 0;
// TODO: fixme
#if 0
    return yb_overlap_shape_bounds(bvh, exclude_self, overlap_ctx, (void*)overlap_cb);
#endif
}

//
// Finds the closest element to a point wihtin a given radius.
//
YGLC_API bool ybc_overlap_first(const yb_scene* scene, const float pt[3],
                                float max_dist, float* dist, int* sid, int* eid,
                                float euv[2], int req_shape) {
    return yb_overlap_first(scene, ym_vec3f(pt), max_dist, dist, sid, eid,
                            (ym_vec2f*)euv, req_shape);
}

//
// Finds if any element overlaps a point wihtin a given radius.
//
YGLC_API bool ybc_overlap_any(const yb_scene* scene, const float pt[3],
                              float max_dist, int req_shape) {
    return yb_overlap_any(scene, ym_vec3f(pt), max_dist, req_shape);
}

//
// Interpolates a vertex property from the given intersection data. Uses
// linear interpolation for lines, baricentric for triangles and quads
// (to correctly treat as two triangles) and copies values for points.
//
YGLC_API void ybc_interpolate_vert(const int* elem, int etype, int eid,
                                   const float euv[2], int vsize,
                                   const float* vert, float* v) {
    return yb_interpolate_vert(elem, etype, eid, ym_vec2f(euv), vsize, vert, v);
}

#endif

#endif
