//
// YOCTO_BVH: ray intersection routines supporting points, lines, triangles and
// quads accelerated by a two-level bounding volume hierarchy (BVH)
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. create a bvh for each shape in your scene with yb_make_shape_bvh
//     for(int i = 0; i < nshapes; i ++)
//         shape_bvhs[i] = yb_make_shape_bvh(pass shape data)
// 2. create the scene bvh by combining shape bvhs with yb_make_scene_bvh
//      bvh = yb_make_scene_bvh(mshapes, shape_bvhs, other shape data)
// 3. perform ray-interseciton tests
//     - use yb_intersect_bvh if you want to know the hit point
//         hit = yb_intersect_bvh(bvh, ray data, out primitive intersection)
//     - use yb_hit_bvh if you only need to know whether there is a hit
//         hit = yb_hit_bvh(bvh, ray data)
// 4. use yb_interpolate_vertex to get interpolated vertex values from the
//    intersection data
//    yb_interpolate_vertex(intersection data, shape data, out interpolated val)
// 5. cleanup with yb_free_bvh
//   - you have to do this for each shape bvh and the scene bvh
//   yb_free_bvh(bvh)
//   for(int i = 0; i < nshapes; i ++) yb_free_bvh(shape_bvhs[i])
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
// and quads. This follows Intel's Embree API.
//
//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use.
// To use the library as a .h/.c pair do the following:
// - to use as a .h, just #define YB_NOINLINE before including this file
// - to use as a .c, just #define YB_IMPLEMENTATION before including this file
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
#ifndef YB_NOINLINE
#define YB_API static inline
#else
#ifdef __cplusplus
#define YB_API extern "C"
#else
#define YB_API
#endif
#endif

#include <stdbool.h>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Element types for shapes
//
enum {
    yb_etype_point = 1,     // points
    yb_etype_line = 2,      // lines
    yb_etype_triangle = 3,  // triangles
    yb_etype_quad = 4       // quads
};

//
// Heuristic stratigy for bvh build
//
enum {
    yb_htype_default = 0,  // default strategy (use this for ray-casting)
    yb_htype_equalnum,     // balanced binary tree
    yb_htype_sah,          // surface area heuristic
    yb_htype_max           // total number of strategies
};

//
// BVH data structure. Implementation details hidden.
//
typedef struct yb_bvh yb_bvh;

//
// Create a BVH for a given shape.
//
// Shapes are indexed meshes with 1, 2, 3, 4 indices respectively for points,
// lines, triangles and quads. Vertices have positions and radius, the latter
// required only for points and lines.
//
// Parameters:
// - nelems: number of elements
// - elem: array of vertex indices
// - etype: shape element type (as per previous enum)
// - pos: array of 3D vertex positions
// - radius: array of vertex radius
// - heuristic: heristic use to build the bvh (0 for deafult)
// - share_mem: whether to share or copy shape data (prefer false if you can)
//
// Return:
// - shape bvh
//
YB_API yb_bvh*
yb_make_shape_bvh(int nelems, int* elem, int etype, int nverts, float* pos,
                  float* radius, int heuristic, bool share_mem);

//
// Create a BVH for a collection of transformed shapes (scene).
//
// Shapes' BVHs can be transformed with transformations matrices. The code
// supports only affine transforms.
//
// Shapes can be tagged with ray_masks, useful for filtering intersections.
// During intersection, if a ray is marked with a mask, shapes are
// intersected only if the bitwise and between the ray masks is true.
//
// Parameters:
// - nbvhs: number of shape bvhs
// - bvhs: array of pointers to shape bvhs
// - xforms: arrays of transform matrices (null for no xforms)
// - ray_masks: array of shape ray masks
// - heuristic: heristic use to build the bvh (0 for deafult)
//
YB_API yb_bvh*
yb_make_scene_bvh(int nbvhs, yb_bvh** bvhs, float* xforms[16], int* ray_masks,
                  int heuristic);

//
// Clears BVH memory.
//
// This will only clear the memory for the given bvh. For scene bvhs, you are
// still responsible for clearing all shape bvhs.
//
// Parameters:
// - bvh: bvh to clean
//
YB_API void
yb_free_bvh(yb_bvh* bvh);

//
// Intersect the scene with a ray finding the closest intersection.
//
// Parameters:
// - bvh: bvh to intersect (scene or shape)
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
YB_API bool
yb_intersect_bvh(const yb_bvh* bvh, const float ray_o[3], const float ray_d[3],
                 float ray_tmin, float ray_tmax, int ray_mask, float* ray_t,
                 int* sid, int* eid, float euv[2]);

//
// Intersect the scene with a ray finding any intersection.
//
// Parameters:
// - bvh: bvh to intersect (scene or shape)
// - ray_o: ray origin
// - ray_d: ray direction
// - ray_o: minimal distance along the ray to consider (0 for all)
// - ray_o: maximal distance along the ray to consider (HUGE_VALF for all)
// - ray_mask: ray mask for use (0 for no mask); see yb_make_scene_bvh
//
// Return:
// - whether we intersect or not
//
YB_API bool
yb_hit_bvh(const yb_bvh* bvh, const float ray_o[3], const float ray_d[3],
           float ray_tmin, float ray_tmax, int ray_mask);

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
YB_API void
yb_interpolate_vert(const int* elem, int etype, int eid, const float euv[2],
                    int vsize, const float* vert, float* v);

//
// Print stats for the given BVH. Mostly useful for debugging or performance
// analysis. Ignore for general purpose.
//
// Parameters:
// - bvh: bvh to print stats for
// - whether to print the tree structure
//
YB_API void
yb_print_bvh_stats(const yb_bvh* bvh, bool print_tree);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

// handles compilation options
#if !defined(YB_NOINLINE) || defined(YB_IMPLEMENTATION)

// include standard files
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// MATH FUNCTIONS SUPPORT
// -----------------------------------------------------------------------------

//
// 2d float vectors
//
typedef struct yb__vec2f {
    float x, y;  // vector components
} yb__vec2f;

//
// 3d float vectors
//
typedef struct yb__vec3f {
    float x, y, z;  // vector components
} yb__vec3f;

//
// Packed vertex data (4 float for position and radius)
//
typedef struct yb__vert4f {
    float x, y, z, r;  // vector components
} yb__vert4f;

//
// 4x4 float matrix
//
typedef struct yb__mat4f {
    float m[16];  // matrix components in row-major order
} yb__mat4f;

//
// Axis-aligned bounding box (AABB) storing min, max corners
//
typedef struct yb__range3f {
    yb__vec3f min, max;  // bouding box min and max corners
} yb__range3f;

//
// Vector subtraction
//
static inline yb__vec3f
yb__sub3f(const yb__vec3f a, const yb__vec3f b) {
    return (yb__vec3f){ a.x - b.x, a.y - b.y, a.z - b.z };
}

//
// Vector dot product
//
static inline float
yb__dot3f(const yb__vec3f a, const yb__vec3f b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

//
// Vector cross product
//
static inline yb__vec3f
yb__cross3f(const yb__vec3f a, const yb__vec3f b) {
    return (yb__vec3f){ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x };
}

//
// Init AABB to an invalid box
//
static inline yb__range3f
yb__rinit3f() {
    return (yb__range3f){ { HUGE_VALF, HUGE_VALF, HUGE_VALF },
                          { -HUGE_VALF, -HUGE_VALF, -HUGE_VALF } };
}

//
// Expand AABB to include a point
//
static inline yb__range3f
yb__rexpand3f(const yb__range3f a, const yb__vec3f b) {
    return (yb__range3f){
        { fminf(a.min.x, b.x), fminf(a.min.y, b.y), fminf(a.min.z, b.z) },
        { fmaxf(a.max.x, b.x), fmaxf(a.max.y, b.y), fmaxf(a.max.z, b.z) }
    };
}

//
// Takes the union of two AABBs
//
static inline yb__range3f
yb__runion3f(const yb__range3f a, const yb__range3f b) {
    return (yb__range3f){ { fminf(a.min.x, b.min.x), fminf(a.min.y, b.min.y),
                            fminf(a.min.z, b.min.z) },
                          { fmaxf(a.max.x, b.max.x), fmaxf(a.max.y, b.max.y),
                            fmaxf(a.max.z, b.max.z) } };
}

//
// AABB center
//
static inline yb__vec3f
yb__rcenter3f(const yb__range3f a) {
    return (yb__vec3f){ (a.max.x + a.min.x) * 0.5f, (a.max.y + a.min.y) * 0.5f,
                        (a.max.z + a.min.z) * 0.5f };
}

//
// AABB diagonal size
//
static inline yb__vec3f
yb__rsize3f(const yb__range3f a) {
    return (yb__vec3f){ a.max.x - a.min.x, a.max.y - a.min.y,
                        a.max.z - a.min.z };
}

//
// Inervert a 4x4 matrix
//
// Implementation Notes:
// - based on the code on this link http://stackoverflow.com/questions/2624422/
//     efficient-4x4-matrix-inverse-affine-transform
//
static inline yb__mat4f
yb__inverse4f(const yb__mat4f m) {
    float a[4][4];
    *(yb__mat4f*)a = m;

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

    return *(yb__mat4f*)b;
}

//
// Transforms a point via an affine matrix
//
static inline yb__vec3f
yb__transform_point3f(const yb__mat4f a, const yb__vec3f b) {
    return (yb__vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z + a.m[3],
                        a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z + a.m[7],
                        a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z + a.m[11] };
}

//
// Transforms a vector via an affine matrix
//
static inline yb__vec3f
yb__transform_vector3f(const yb__mat4f a, const yb__vec3f b) {
    return (yb__vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z,
                        a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z,
                        a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z };
}

//
// Clamp a value x between m, M inclusive
//
static inline float
yb__fclampf(float x, float m, float M) {
    return fmin(M, fmax(m, x));
}

//
// Compute the point on the ray ray_o, ray_d at distance t
//
static inline yb__vec3f
yb__eval_ray(const yb__vec3f ray_o, const yb__vec3f ray_d, float t) {
    return (yb__vec3f){ ray_o.x + ray_d.x * t, ray_o.y + ray_d.y * t,
                        ray_o.z + ray_d.z * t };
}

// -----------------------------------------------------------------------------
// ELEMENT-WISE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Intersect a ray with a point (approximate)
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
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
static inline bool
yb__intersect_point(const yb__vec3f ray_o, const yb__vec3f ray_d,
                    float ray_tmin, float ray_tmax, const yb__vec3f p, float r,
                    float* ray_t, yb__vec2f* euv) {
    // find parameter for line-point minimum distance
    yb__vec3f w = yb__sub3f(p, ray_o);
    float t = yb__dot3f(w, ray_d) / yb__dot3f(ray_d, ray_d);

    // exit if not within bounds
    if (t < ray_tmin || t > ray_tmax) return false;

    // test for line-point distance vs point radius
    yb__vec3f rp = yb__eval_ray(ray_o, ray_d, t);
    yb__vec3f prp = yb__sub3f(p, rp);
    if (yb__dot3f(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = (yb__vec2f){ 0, 0 };

    return true;
}

//
// Intersect a ray with a line
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
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
static inline bool
yb__intersect_line(const yb__vec3f ray_o, const yb__vec3f ray_d, float ray_tmin,
                   float ray_tmax, const yb__vec3f v0, const yb__vec3f v1,
                   float r0, float r1, float* ray_t, yb__vec2f* euv) {
    // setup intersection params
    yb__vec3f u = ray_d;
    yb__vec3f v = yb__sub3f(v1, v0);
    yb__vec3f w = yb__sub3f(ray_o, v0);

    // compute values to solve a linear system
    float a = yb__dot3f(u, u);
    float b = yb__dot3f(u, v);
    float c = yb__dot3f(v, v);
    float d = yb__dot3f(u, w);
    float e = yb__dot3f(v, w);
    float det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute parameters on both ray and segment
    float t = (b * e - c * d) / det;
    float s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray_tmin || t > ray_tmax) return false;

    // clamp segment param to segment corners
    s = yb__fclampf(s, 0, 1);

    // compute segment-segment distance on the closest points
    yb__vec3f p0 = yb__eval_ray(ray_o, ray_d, t);
    yb__vec3f p1 = yb__eval_ray(v0, yb__sub3f(v1, v0), s);
    yb__vec3f p01 = yb__sub3f(p0, p1);

    // check with the line radius at the same point
    float r = r0 * (1 - s) + r1 * s;
    if (yb__dot3f(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = (yb__vec2f){ s, 0 };

    return true;
}

//
// Intersect a ray with a triangle
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
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
static inline bool
yb__intersect_triangle(const yb__vec3f ray_o, const yb__vec3f ray_d,
                       float ray_tmin, float ray_tmax, const yb__vec3f v0,
                       const yb__vec3f v1, const yb__vec3f v2, float* ray_t,
                       yb__vec2f* euv) {
    // compute triangle edges
    yb__vec3f edge1 = yb__sub3f(v1, v0);
    yb__vec3f edge2 = yb__sub3f(v2, v0);

    // compute determinant to solve a linear system
    yb__vec3f pvec = yb__cross3f(ray_d, edge2);
    float det = yb__dot3f(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    float inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    yb__vec3f tvec = yb__sub3f(ray_o, v0);
    float u = yb__dot3f(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    yb__vec3f qvec = yb__cross3f(tvec, edge1);
    float v = yb__dot3f(ray_d, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0) return false;

    // compute and check ray parameter
    float t = yb__dot3f(edge2, qvec) * inv_det;
    if (t < ray_tmin || t > ray_tmax) return false;

    // intersection occurred: set params and exit
    *ray_t = t;
    *euv = (yb__vec2f){ u, v };

    return true;
}

//
// Intersect a ray with a quad, represented as two triangles
//
// Parameters:
// - ray_o, ray_d: ray origin and direction
// - ray_tmin, ray_tmax: ray parameter min, max range
// - v0, v1, v2, v3: quad vertices
//
// Out Parameters:
// - ray_t: ray parameter at the intersection point
// - euv: baricentric coordinates (as described above)
//
// Returns:
// - whether the intersection occurred
//
// Notes:
// - out parameters and only writtent o if an intersection occurs
// - might become deprecated in the future
//
static inline bool
yb__intersect_quad(const yb__vec3f ray_o, const yb__vec3f ray_d, float ray_tmin,
                   float ray_tmax, const yb__vec3f v0, const yb__vec3f v1,
                   const yb__vec3f v2, const yb__vec3f v3, float* ray_t,
                   yb__vec2f* euv) {
    // test first triangle
    bool hit = false;
    if (yb__intersect_triangle(ray_o, ray_d, ray_tmin, ray_tmax, v0, v1, v3,
                               &ray_tmax, euv)) {
        hit = true;
        *ray_t = ray_tmax;
    }

    // test second triangle
    if (yb__intersect_triangle(ray_o, ray_d, ray_tmin, ray_tmax, v2, v3, v1,
                               &ray_tmax, euv)) {
        hit = true;
        *ray_t = ray_tmax;
        // flip coordinates to map to [0,0]x[1,1]
        *euv = (yb__vec2f){ 1 - euv->x, 1 - euv->y };
    }

    return hit;
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
static inline bool
yb__intersect_bbox(const yb__vec3f ray_o, const yb__vec3f ray_d, float ray_tmin,
                   float ray_tmax, const yb__vec3f bbox_min,
                   const yb__vec3f bbox_max) {
    // set up convenient pointers for looping over axes
    const float* _ray_o = &ray_o.x;
    const float* _ray_d = &ray_d.x;
    const float* _bbox_min = &bbox_min.x;
    const float* _bbox_max = &bbox_max.x;

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
        ray_tmin = t0 > ray_tmin ? t0 : ray_tmin;
        ray_tmax = t1 < ray_tmax ? t1 : ray_tmax;
        // if intersection is empty, exit
        if (ray_tmin > ray_tmax) return false;
    }

    // passed all planes, then intersection occurred
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
    yb__ntype_internal = 0,         // node is internal (not leaf)
    yb__ntype_point = 1,            // points
    yb__ntype_line = 2,             // lines
    yb__ntype_triangle = 3,         // triangles
    yb__ntype_quad = 4,             // quads
    yb__ntype_shared_point = 5,     // pointer to external points
    yb__ntype_shared_line = 6,      // pointer to external lines
    yb__ntype_shared_triangle = 7,  // pointer to external triangles
    yb__ntype_shared_quad = 8,      // pointer to external quads
    yb__ntype_bvh = 9,              // pointers to other bvhs
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
typedef struct yb__bvhn {
    yb__range3f bbox;  // bounding box
    uint32_t start;    // index to the first sorted primitive/node
    uint16_t count;    // number of primitives/nodes
    uint8_t isleaf;    // whether it is a leaf
    uint8_t axis;      // plit axis
} yb__bvhn;

//
// Structure that contains a reference to an external trnaformed BVH.
// See yb_bvh for more details.
//
typedef struct yb__xformed_bvh {
    const yb_bvh* bvh;    // external BVH pointer
    bool xformed;         // whether it is tranformed
    int pid;              // BVH primitive id to return upon intersection
    int ray_mask;         // BVH ray mask to filter intersections
    yb__mat4f xform;      // transform
    yb__mat4f inv_xform;  // inverse transform
} yb__xformed_bvh;

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
// - we handle three types of BVHs, shape, shape_ref and scene, using a union
// - shape BVHs contain sorted geometric primitives for fast intersection
// - shape_ref BVHs contain pointers to external geometric data, to save memory
// - scene BVHs contain other BVHs to build the two level hiearchy
// - all BVH types shares the same type of internal nodes
//
struct yb_bvh {
    // bvh data
    int ntype;        // type of primitive data contained by bvh
    int nprims;       // number of primitivrs
    int nnodes;       // number of internal nodes
    yb__bvhn* nodes;  // sorted array of internal nodes
    // primitive data
    union {
        // shared elem data
        struct yb__shaperef_bvh {
            int* sorted_prim;  // sorted elements
            int* elems;        // shape elements
            yb__vec3f* pos;    // vertex pos
            float* radius;     // vertex radius
        } shape_ref;
        // elem data
        struct yb__shape_bvh {
            int* pid;          // elem ids
            int* elems;        // shape elements
            yb__vert4f* vert;  // vertex
        } shape;
        // scene bvh data
        struct yb__scene_bvh {
            yb__xformed_bvh* shape_bvhs;  // array of shape BVHs
        } scene;
    };
};

// -----------------------------------------------------------------------------
// BVH BUILD FUNCTIONS
// -----------------------------------------------------------------------------

//
// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
//
typedef struct yb__bound_prim {
    yb__range3f bbox;      // bounding box
    yb__vec3f center;      // bounding box center (for faster sort)
    int pid;               // primitive id
    float sah_cost_left;   // buffer for sah heuristic costs
    float sah_cost_right;  // buffer for sah heuristic costs
} yb__bound_prim;

//
// Shell sort BVH array a of length n along the axis axis. Implementation based
// on http://rosettacode.org/wiki/Sorting_algorithms/Shell_sort
//
static inline void
yb__shellsort_bound_prim(yb__bound_prim* a, int n, int axis) {
    int h, i, j;
    yb__bound_prim t;
    for (h = n; h /= 2;) {
        for (i = h; i < n; i++) {
            t = a[i];
            for (j = i;
                 j >= h && (&t.center.x)[axis] < (&a[j - h].center.x)[axis];
                 j -= h) {
                a[j] = a[j - h];
            }
            a[j] = t;
        }
    }
}

//
// Heap sort BVH array a of length n along the axis axis. Implementation based
// on https://en.wikipedia.org/wiki/Heapsort
//
static inline void
yb__heapsort_bound_prim(yb__bound_prim* arr, int count, int axis) {
    if (!count) return;

    yb__bound_prim t;
    unsigned int n = count, parent = count / 2, index, child;
    while (true) {
        if (parent > 0) {
            t = arr[--parent];
        } else {
            n--;
            if (n == 0) return;
            t = arr[n];
            arr[n] = arr[0];
        }
        index = parent;
        child = index * 2 + 1;
        while (child < n) {
            if (child + 1 < n &&
                (&arr[child + 1].center.x)[axis] > (&arr[child].center.x)[axis])
                child++;
            if ((&arr[child].center.x)[axis] > (&t.center.x)[axis]) {
                arr[index] = arr[child];
                index = child;
                child = index * 2 + 1;
            } else
                break;
        }
        arr[index] = t;
    }
}

// choose a sort algorithm
#define yb__sort_bound_prim yb__heapsort_bound_prim

//
// Given an array sorted_prim of primitives to split between the elements
// start and end, determines the split axis axis, split primitive index mid
// based on the heuristic heuristic. Supports balanced tree (equalnum) and
// Surface-Area Heuristic.
//
static inline void
yb__split_axis(yb__bound_prim* sorted_prim, int start, int end, int* axis,
               int* mid, int heuristic) {
    *axis = -1;
    *mid = -1;
    switch (heuristic) {
        // balanced tree split: find the largest axis of the bounding box
        // and split along this one right in the middle
        case yb_htype_equalnum: {
            // compute primintive bounds
            yb__range3f bbox = yb__rinit3f();
            for (int i = start; i < end; i++) {
                bbox = yb__rexpand3f(bbox, sorted_prim[i].center);
            }
            // split along largest
            yb__vec3f bvh_size = yb__rsize3f(bbox);
            if (bvh_size.x >= bvh_size.y && bvh_size.x >= bvh_size.z)
                *axis = 0;
            else if (bvh_size.y >= bvh_size.x && bvh_size.y >= bvh_size.z)
                *axis = 1;
            else if (bvh_size.z >= bvh_size.x && bvh_size.z >= bvh_size.y)
                *axis = 2;
            else
                assert(false);
            *mid = (start + end) / 2;
        } break;
        case yb_htype_default:
        case yb_htype_sah: {
            // surface area heuristic: estimate the cost of splitting
            // along each of the axis and pick the one with best expected
            // performance
            float min_cost = HUGE_VALF;
            int count = end - start;
            for (int a = 0; a < 3; a++) {
                yb__sort_bound_prim(sorted_prim + start, end - start, a);
                yb__range3f sbbox = yb__rinit3f();
                // to avoid an O(n^2) computation, use sweaps to compute the
                // cost,
                // first smallest to largest, then largest to smallest
                for (int i = 0; i < count; i++) {
                    sbbox = yb__runion3f(sbbox, sorted_prim[start + i].bbox);
                    yb__vec3f sbbox_size = yb__rsize3f(sbbox);
                    sorted_prim[start + i].sah_cost_left =
                        sbbox_size.x * sbbox_size.y +
                        sbbox_size.x * sbbox_size.z +
                        sbbox_size.y * sbbox_size.z;
                    sorted_prim[start + i].sah_cost_left *= i + 1;
                }
                // the other sweep
                sbbox = yb__rinit3f();
                for (int i = 0; i < count; i++) {
                    sbbox = yb__runion3f(sbbox, sorted_prim[end - 1 - i].bbox);
                    yb__vec3f sbbox_size = yb__rsize3f(sbbox);
                    sorted_prim[end - 1 - i].sah_cost_right =
                        sbbox_size.x * sbbox_size.y +
                        sbbox_size.x * sbbox_size.z +
                        sbbox_size.y * sbbox_size.z;
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
        } break;
        default: assert(false); break;
    }
    assert(*axis >= 0 && *mid > 0);
}

//
// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
//
static inline void
yb__make_node(yb__bvhn* node, int* nnodes, yb__bvhn* nodes,
              yb__bound_prim* sorted_prims, int start, int end, int heuristic) {
    // compute node bounds
    node->bbox = yb__rinit3f();
    for (int i = start; i < end; i++) {
        node->bbox = yb__runion3f(node->bbox, sorted_prims[i].bbox);
    }

    // decide whether to create a leaf
    if (end - start <= YB__BVH_MINPRIMS) {
        // makes a leaf node
        node->isleaf = true;
        node->start = start;
        node->count = end - start;
    } else {
        // makes an internal node
        node->isleaf = false;
        // choose the split axis and position
        int axis, mid;
        yb__split_axis(sorted_prims, start, end, &axis, &mid, heuristic);
        // sort primitives along the given axis
        yb__sort_bound_prim(sorted_prims + start, end - start, axis);
        // perform the splits by preallocating the child nodes and recurring
        node->axis = axis;
        node->start = *nnodes;
        node->count = 2;
        *nnodes += 2;
        // build child nodes
        yb__make_node(nodes + node->start, nnodes, nodes, sorted_prims, start,
                      mid, heuristic);
        yb__make_node(nodes + node->start + 1, nnodes, nodes, sorted_prims, mid,
                      end, heuristic);
    }
}

//
// Makes a shape BVH. Public function whose interface is described above.
//
YB_API yb_bvh*
yb_make_shape_bvh(int nelems, int* elem, int etype, int nverts, float* pos,
                  float* radius, int heuristic, bool share_mem) {
    // create bounded primitivrs used in BVH build
    yb__bound_prim* bound_prims =
        (yb__bound_prim*)calloc(nelems, sizeof(yb__bound_prim));
    assert(pos && elem);
    // compute bounds for each primitive type
    switch (etype) {
        case yb_etype_point: {
            // point bounds are computed as small spheres
            assert(radius);
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i;
                yb__bound_prim* prim = bound_prims + i;
                prim->pid = i;
                prim->bbox = yb__rinit3f();
                yb__vec3f p = ((yb__vec3f*)pos)[f[0]];
                float r = radius[f[0]];
                yb__vec3f pm = { p.x - r, p.y - r, p.z - r };
                yb__vec3f pM = { p.x + r, p.y + r, p.z + r };
                prim->bbox = yb__rexpand3f(prim->bbox, pm);
                prim->bbox = yb__rexpand3f(prim->bbox, pM);
                prim->center = yb__rcenter3f(prim->bbox);
            }
        } break;
        case yb_etype_line: {
            // line bounds are computed as thick rods
            assert(radius);
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i * 2;
                yb__bound_prim* prim = bound_prims + i;
                prim->pid = i;
                prim->bbox = yb__rinit3f();
                for (int c = 0; c < 2; c++) {
                    yb__vec3f p = ((yb__vec3f*)pos)[f[c]];
                    float r = radius[f[c]];
                    yb__vec3f pm = { p.x - r, p.y - r, p.z - r };
                    yb__vec3f pM = { p.x + r, p.y + r, p.z + r };
                    prim->bbox = yb__rexpand3f(prim->bbox, pm);
                    prim->bbox = yb__rexpand3f(prim->bbox, pM);
                }
                prim->center = yb__rcenter3f(prim->bbox);
            }
        } break;
        case yb_etype_triangle: {
            // triangle bounds are computed by including their vertices
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i * 3;
                yb__bound_prim* prim = bound_prims + i;
                prim->pid = i;
                prim->bbox = yb__rinit3f();
                for (int c = 0; c < 3; c++) {
                    prim->bbox =
                        yb__rexpand3f(prim->bbox, ((yb__vec3f*)pos)[f[c]]);
                }
                prim->center = yb__rcenter3f(prim->bbox);
            }
        } break;
        case yb_etype_quad: {
            // quad bounds are computed by including their vertices
            for (int i = 0; i < nelems; i++) {
                int* f = elem + i * 4;
                yb__bound_prim* prim = bound_prims + i;
                prim->pid = i;
                prim->bbox = yb__rinit3f();
                for (int c = 0; c < 4; c++) {
                    prim->bbox =
                        yb__rexpand3f(prim->bbox, ((yb__vec3f*)pos)[f[c]]);
                }
                prim->center = yb__rcenter3f(prim->bbox);
            }
        } break;
        default: assert(false);
    }

    // allocates the BVH and set its primitive type
    yb_bvh* bvh = (yb_bvh*)calloc(1, sizeof(yb_bvh));
    bvh->nprims = nelems;
    bvh->ntype = etype - yb_etype_point +
                 ((share_mem) ? yb__ntype_shared_point : yb__ntype_point);

    // allocate nodes (over-allocate now then shrink)
    bvh->nnodes = 0;
    bvh->nodes = (yb__bvhn*)calloc(bvh->nprims * 2, sizeof(yb__bvhn));

    // start recursive splitting
    bvh->nnodes = 1;
    yb__make_node(bvh->nodes, &bvh->nnodes, bvh->nodes, bound_prims, 0,
                  bvh->nprims, heuristic);

    // shrink back
    bvh->nodes = (yb__bvhn*)realloc(bvh->nodes, sizeof(yb__bvhn) * bvh->nnodes);

    // init sorted element arrays
    if (share_mem) {
        // for shared memory, stored pointer to the external data
        bvh->shape_ref.elems = elem;
        bvh->shape_ref.pos = (yb__vec3f*)pos;
        bvh->shape_ref.radius = radius;
        // store the sorted primitive order for BVH walk
        bvh->shape_ref.sorted_prim = (int*)calloc(bvh->nprims, sizeof(int));
        for (int i = 0; i < bvh->nprims; i++) {
            bvh->shape_ref.sorted_prim[i] = bound_prims[i].pid;
        }
    } else {
        // for non-shared memory, copy the sorted primitive data in internal
        // sorted arrays
        bvh->shape.vert = (yb__vert4f*)malloc(nverts * sizeof(yb__vert4f));
        for (int i = 0; i < nverts; i++) {
            bvh->shape.vert[i].x = pos[i * 3 + 0];
            bvh->shape.vert[i].y = pos[i * 3 + 1];
            bvh->shape.vert[i].z = pos[i * 3 + 2];
            bvh->shape.vert[i].r = (radius) ? radius[i] : 0;
        }
        // stored sorted primitive data (elems) and indices to external
        // elements to return upon intersection
        bvh->shape.pid = (int*)calloc(bvh->nprims, sizeof(int));
        bvh->shape.elems = (int*)calloc(bvh->nprims, sizeof(int) * etype);
        for (int i = 0; i < bvh->nprims; i++) {
            int pid = bound_prims[i].pid;
            bvh->shape.pid[i] = pid;
            memcpy(bvh->shape.elems + i * etype, elem + pid * etype,
                   sizeof(int) * etype);
        }
    }

    // cleanup
    free(bound_prims);

    // done
    return bvh;
}

//
// Makes a scene BVH. Public function whose interface is described above.
//
YB_API yb_bvh*
yb_make_scene_bvh(int nbvhs, yb_bvh** bvhs, float* xforms[16], int* ray_mask,
                  int heuristic) {
    // compute bounds for all transformed shape bvhs
    yb__bound_prim* bound_prims =
        (yb__bound_prim*)calloc(nbvhs, sizeof(yb__bound_prim));
    assert(bvhs);
    for (int i = 0; i < nbvhs; i++) {
        yb__bound_prim* prim = bound_prims + i;
        prim->pid = i;
        prim->bbox = bvhs[i]->nodes[0].bbox;
        // apply transformations if needed (estimate trasform BVH bounds as
        // corners bounds, overestimate, but best effort for now)
        if (xforms && xforms[i]) {
            yb__vec3f corners[8] = {
                { prim->bbox.min.x, prim->bbox.min.y, prim->bbox.min.z },
                { prim->bbox.min.x, prim->bbox.min.y, prim->bbox.max.z },
                { prim->bbox.min.x, prim->bbox.max.y, prim->bbox.min.z },
                { prim->bbox.min.x, prim->bbox.max.y, prim->bbox.max.z },
                { prim->bbox.max.x, prim->bbox.min.y, prim->bbox.min.z },
                { prim->bbox.max.x, prim->bbox.min.y, prim->bbox.max.z },
                { prim->bbox.max.x, prim->bbox.max.y, prim->bbox.min.z },
                { prim->bbox.max.x, prim->bbox.max.y, prim->bbox.max.z },
            };
            yb__mat4f xform = *(yb__mat4f*)xforms[i];
            prim->bbox = yb__rinit3f();
            for (int j = 0; j < 8; j++) {
                prim->bbox = yb__rexpand3f(
                    prim->bbox, yb__transform_point3f(xform, corners[j]));
            }
        }
        prim->center = yb__rcenter3f(prim->bbox);
    }

    // allocates the BVH and set its primitive type
    yb_bvh* bvh = (yb_bvh*)calloc(1, sizeof(yb_bvh));
    bvh->nprims = nbvhs;
    bvh->ntype = yb__ntype_bvh;

    // allocate nodes (over-allocate now then shrink)
    bvh->nnodes = 0;
    bvh->nodes = (yb__bvhn*)calloc(bvh->nprims * 2, sizeof(yb__bvhn));

    // start recursive splitting
    bvh->nnodes = 1;
    yb__make_node(bvh->nodes, &bvh->nnodes, bvh->nodes, bound_prims, 0,
                  bvh->nprims, heuristic);

    // shrink back
    bvh->nodes = (yb__bvhn*)realloc(bvh->nodes, sizeof(yb__bvhn) * bvh->nnodes);

    // pack sorted primimitives in the internal bvhs array
    bvh->scene.shape_bvhs =
        (yb__xformed_bvh*)calloc(bvh->nprims, sizeof(yb__xformed_bvh));
    for (int i = 0; i < bvh->nprims; i++) {
        int pid = bound_prims[i].pid;
        bvh->scene.shape_bvhs[i].bvh = bvhs[pid];
        bvh->scene.shape_bvhs[i].pid = pid;
        bvh->scene.shape_bvhs[i].ray_mask = (ray_mask) ? ray_mask[pid] : 0;
        bvh->scene.shape_bvhs[i].xformed = xforms && xforms[pid];
        if (bvh->scene.shape_bvhs[i].xformed) {
            bvh->scene.shape_bvhs[i].xform = *(yb__mat4f*)xforms[pid];
            bvh->scene.shape_bvhs[i].inv_xform =
                yb__inverse4f(bvh->scene.shape_bvhs[i].xform);
        }
    }

    // cleanup
    free(bound_prims);

    // done
    return bvh;
}

//
// FreeBVH data.
//
YB_API void
yb_free_bvh(yb_bvh* bvh) {
    free(bvh->nodes);
    switch (bvh->ntype) {
        case yb__ntype_point:
        case yb__ntype_line:
        case yb__ntype_triangle:
        case yb__ntype_quad: {
            free(bvh->shape.pid);
            free(bvh->shape.elems);
            free(bvh->shape.vert);
        } break;
        case yb__ntype_shared_point:
        case yb__ntype_shared_line:
        case yb__ntype_shared_triangle:
        case yb__ntype_shared_quad: {
            free(bvh->shape_ref.sorted_prim);
        } break;
        case yb__ntype_bvh: {
            free(bvh->scene.shape_bvhs);
        } break;
        default: { assert(false); }
    }
}

// -----------------------------------------------------------------------------
// BVH INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

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
static inline bool
yb__intersect_bvh(const yb_bvh* bvh, bool early_exit, const yb__vec3f ray_o,
                  const yb__vec3f ray_d, float ray_tmin, float ray_tmax,
                  int ray_mask, float* ray_t, int* sid, int* eid,
                  yb__vec2f* euv) {
    // node stack
    int node_stack[64];
    int node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    bool hit = false;

    // walking stack
    while (node_cur && !(early_exit && hit)) {
        // exit if needed
        if (early_exit && hit) return hit;

        // grab node
        yb__bvhn* node = bvh->nodes + node_stack[--node_cur];

        // intersect bbox
        if (!yb__intersect_bbox(ray_o, ray_d, ray_tmin, ray_tmax,
                                node->bbox.min, node->bbox.max))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        int ntype = (node->isleaf) ? bvh->ntype : yb__ntype_internal;
        switch (ntype) {
            // internal node
            case yb__ntype_internal: {
                // for internal nodes, attempts to proceed along the
                // split axis from smallest to largest nodes
                if ((&ray_d.x)[node->axis] >= 0) {
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
            } break;
            case yb__ntype_point: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int* f = bvh->shape.elems + idx * 1;
                    yb__vert4f* vert = bvh->shape.vert;
                    yb__vert4f* v0 = vert + f[0];
                    if (yb__intersect_point(ray_o, ray_d, ray_tmin, ray_tmax,
                                            *(yb__vec3f*)v0, v0->r, &ray_tmax,
                                            euv)) {
                        hit = true;
                        *eid = bvh->shape.pid[idx];
                    }
                }
            } break;
            case yb__ntype_line: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int* f = bvh->shape.elems + idx * 2;
                    yb__vert4f* vert = bvh->shape.vert;
                    yb__vert4f *v0 = vert + f[0], *v1 = vert + f[1];
                    if (yb__intersect_line(ray_o, ray_d, ray_tmin, ray_tmax,
                                           *(yb__vec3f*)v0, *(yb__vec3f*)v1,
                                           v0->r, v1->r, &ray_tmax, euv)) {
                        hit = true;
                        *eid = bvh->shape.pid[idx];
                    }
                }
            } break;
            case yb__ntype_triangle: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int* f = bvh->shape.elems + idx * 3;
                    yb__vert4f* vert = bvh->shape.vert;
                    yb__vert4f *v0 = vert + f[0], *v1 = vert + f[1],
                               *v2 = vert + f[2];
                    if (yb__intersect_triangle(
                            ray_o, ray_d, ray_tmin, ray_tmax, *(yb__vec3f*)v0,
                            *(yb__vec3f*)v1, *(yb__vec3f*)v2, &ray_tmax, euv)) {
                        hit = true;
                        *eid = bvh->shape.pid[idx];
                    }
                }
            } break;
            case yb__ntype_quad: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int* f = bvh->shape.elems + idx * 4;
                    yb__vert4f* vert = bvh->shape.vert;
                    yb__vert4f *v0 = vert + f[0], *v1 = vert + f[1],
                               *v2 = vert + f[2], *v3 = vert + f[3];
                    if (yb__intersect_quad(ray_o, ray_d, ray_tmin, ray_tmax,
                                           *(yb__vec3f*)v0, *(yb__vec3f*)v1,
                                           *(yb__vec3f*)v2, *(yb__vec3f*)v3,
                                           &ray_tmax, euv)) {
                        hit = true;
                        *eid = bvh->shape.pid[idx];
                    }
                }
            } break;
            case yb__ntype_shared_point: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int pid = bvh->shape_ref.sorted_prim[idx];
                    int* f = bvh->shape_ref.elems + pid * 1;
                    yb__vec3f* pos = bvh->shape_ref.pos;
                    float* radius = bvh->shape_ref.radius;
                    if (yb__intersect_point(ray_o, ray_d, ray_tmin, ray_tmax,
                                            pos[f[0]], radius[f[0]], &ray_tmax,
                                            euv)) {
                        hit = true;
                        *eid = pid;
                    }
                }
            } break;
            case yb__ntype_shared_line: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int pid = bvh->shape_ref.sorted_prim[idx];
                    int* f = bvh->shape_ref.elems + pid * 2;
                    yb__vec3f* pos = bvh->shape_ref.pos;
                    float* radius = bvh->shape_ref.radius;
                    if (yb__intersect_line(ray_o, ray_d, ray_tmin, ray_tmax,
                                           pos[f[0]], pos[f[1]], radius[f[0]],
                                           radius[f[1]], &ray_tmax, euv)) {
                        hit = true;
                        *eid = pid;
                    }
                }
            } break;
            case yb__ntype_shared_triangle: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int pid = bvh->shape_ref.sorted_prim[idx];
                    int* f = bvh->shape_ref.elems + pid * 3;
                    yb__vec3f* pos = bvh->shape_ref.pos;
                    if (yb__intersect_triangle(ray_o, ray_d, ray_tmin, ray_tmax,
                                               pos[f[0]], pos[f[1]], pos[f[2]],
                                               &ray_tmax, euv)) {
                        hit = true;
                        *eid = pid;
                    }
                }
            } break;
            case yb__ntype_shared_quad: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    int pid = bvh->shape_ref.sorted_prim[idx];
                    int* f = bvh->shape_ref.elems + pid * 4;
                    yb__vec3f* pos = bvh->shape_ref.pos;
                    if (yb__intersect_quad(ray_o, ray_d, ray_tmin, ray_tmax,
                                           pos[f[0]], pos[f[1]], pos[f[2]],
                                           pos[f[3]], &ray_tmax, euv)) {
                        hit = true;
                        *eid = pid;
                    }
                }
            } break;
            case yb__ntype_bvh: {
                for (int i = 0; i < node->count; i++) {
                    int idx = node->start + i;
                    yb__vec3f tray_o, tray_d;
                    const yb__xformed_bvh* prim = &bvh->scene.shape_bvhs[idx];
                    if (ray_mask && !(ray_mask & prim->ray_mask)) continue;
                    if (prim->xformed) {
                        tray_o = yb__transform_point3f(prim->inv_xform, ray_o);
                        tray_d = yb__transform_vector3f(prim->inv_xform, ray_d);
                    } else {
                        tray_o = ray_o;
                        tray_d = ray_d;
                    }
                    if (yb__intersect_bvh(prim->bvh, early_exit, tray_o, tray_d,
                                          ray_tmin, ray_tmax, ray_mask,
                                          &ray_tmax, sid, eid, euv)) {
                        hit = true;
                        *sid = prim->pid;
                    }
                }
            } break;
            default: { assert(false); } break;
        }
    }

    // update ray distance for outside
    if (hit) *ray_t = ray_tmax;

    return hit;
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YB_API bool
yb_intersect_bvh(const yb_bvh* bvh, const float ray_o[3], const float ray_d[3],
                 float ray_tmin, float ray_tmax, int ray_mask, float* ray_t,
                 int* sid, int* eid, float euv[2]) {
    return yb__intersect_bvh(bvh, false, *(yb__vec3f*)ray_o, *(yb__vec3f*)ray_d,
                             ray_tmin, ray_tmax, ray_mask, ray_t, sid, eid,
                             (yb__vec2f*)euv);
}

//
// Find any ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YB_API bool
yb_hit_bvh(const yb_bvh* bvh, const float ray_o[3], const float ray_d[3],
           float ray_tmin, float ray_tmax, int ray_mask) {
    // check intersection
    int sid, eid;
    float ray_t;
    yb__vec2f euv;
    return yb__intersect_bvh(bvh, true, *(yb__vec3f*)ray_o, *(yb__vec3f*)ray_d,
                             ray_tmin, ray_tmax, ray_mask, &ray_t, &sid, &eid,
                             &euv);
}

// -----------------------------------------------------------------------------
// VERTEX PROPERTY INTERPOLATION
// -----------------------------------------------------------------------------

//
// Interpolate vertex properties at the intersection.
// Public function whose interface is described above.
//
YB_API void
yb_interpolate_vert(const int* elem, int etype, int eid, const float euv[2],
                    int vsize, const float* vert, float* v) {
    for (int c = 0; c < vsize; c++) v[c] = 0;
    switch (etype) {
        case yb_etype_point: {
            const int* f = elem + eid;
            const float* vf = vert + f[0] * vsize;
            for (int c = 0; c < vsize; c++) v[c] += vf[c];
        } break;
        case yb_etype_line: {
            const int* f = elem + eid * 2;
            const float* vf[] = { vert + f[0] * vsize, vert + f[1] * vsize };
            const float w[2] = { 1 - euv[0], euv[0] };
            for (int i = 0; i < 2; i++) {
                for (int c = 0; c < vsize; c++) v[c] += w[i] * vf[i][c];
            }
        } break;
        case yb_etype_triangle: {
            const int* f = elem + eid * 3;
            const float* vf[] = { vert + f[0] * vsize, vert + f[1] * vsize,
                                  vert + f[2] * vsize };
            const float w[3] = { 1 - euv[0] - euv[1], euv[0], euv[1] };
            for (int i = 0; i < 3; i++) {
                for (int c = 0; c < vsize; c++) v[c] += w[i] * vf[i][c];
            }
        } break;
        case yb_etype_quad: {
            const int* f = elem + eid * 4;
            const float* vf[3];
            float w[3];
            if (f[2] == f[3] || euv[0] + euv[1] <= 1) {
                vf[0] = vert + f[0] * vsize;
                vf[1] = vert + f[1] * vsize;
                vf[2] = vert + f[3] * vsize;
                w[0] = 1 - euv[0] - euv[1];
                w[1] = euv[0];
                w[2] = euv[1];
            } else {
                vf[0] = vert + f[2] * vsize;
                vf[1] = vert + f[3] * vsize;
                vf[2] = vert + f[1] * vsize;
                w[0] = 1 - (1 - euv[0]) - (1 - euv[1]);
                w[1] = 1 - euv[0];
                w[2] = 1 - euv[1];
            }
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
// Implementation of the function below
//
static inline void
yb__print_bvh_stats_rec(const yb_bvh* bvh, const yb__bvhn* n, int d) {
    for (int i = 0; i < d; i++) printf(" ");
    printf("%s\n", (!n->isleaf) ? "-" : "*");
    if (!n->isleaf) {
        for (int i = 0; i < n->count; i++) {
            int idx = n->start + i;
            yb__print_bvh_stats_rec(bvh, bvh->nodes + idx, d + 1);
        }
    }
}

//
// Implementation of the function below
//
static inline void
yb__collect_bvh_stats(const yb_bvh* bvh, const yb__bvhn* n, int d, int* nl,
                      int* ni, int* sump, int* mind, int* maxd, int* sumd) {
    if (n->isleaf) {
        *mind = (*mind < d) ? *mind : d;
        *maxd = (*maxd > d) ? *maxd : d;
        *sumd += d;
        (*nl)++;
        (*sump) += n->count;
    } else {
        (*ni)++;
        for (int i = 0; i < n->count; i++) {
            int idx = n->start + i;
            yb__collect_bvh_stats(bvh, bvh->nodes + idx, d + 1, nl, ni, sump,
                                  mind, maxd, sumd);
        }
    }
}

//
// Print bvh for debugging
//
YB_API void
yb_print_bvh_stats(const yb_bvh* bvh, bool print_tree) {
    printf("nodes:  %d\nprims: %d\n", bvh->nnodes, bvh->nprims);
    int nl = 0, ni = 0, sump = 0, mind = 0, maxd = 0, sumd = 0;
    yb__collect_bvh_stats(bvh, bvh->nodes, 0, &nl, &ni, &sump, &mind, &maxd,
                          &sumd);
    printf("leaves:  %d\ninternal: %d\nprims per leaf: %f\n", nl, ni,
           (float)sump / (float)nl);
    printf("min depth:  %d\nmax depth: %d\navg depth: %f\n", mind, maxd,
           (float)sumd / (float)nl);
    yb__print_bvh_stats_rec(bvh, bvh->nodes, 0);
    if (bvh->ntype == yb__ntype_bvh) {
        for (int i = 0; i < bvh->nprims; i++)
            yb_print_bvh_stats(bvh->scene.shape_bvhs[i].bvh, print_tree);
    }
}

#endif

#endif
