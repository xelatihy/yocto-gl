//
// YOCTO_BVH: ray-intersection and closet-point routines supporting points,
// lines, triangles accelerated by a two-level bounding volume
// hierarchy (BVH)
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. create a scene
//      scene = yobj::scene()
// 2. for each shape, add shape data and transforms
//      for(int i = 0; i < nshapes; i ++) {
//          shape = yobj::shape()
//          add shape data as views (only one primitive type though)
//          scene.shapes.push_back(shape)
//          scene.xforms.push_back(shape transform)
//          scene.inv_xforms.push_back(shape inverse transform)
//      - the shape index returned by the interface is the same as
//        the position in the shape
// 3. build the bvh with build_bvh using the specified heuristic (or default)
// 4.a. perform ray-interseciton tests
//     - use intersect_first if you want to know the hit point
//         hit = intersect_first(scene, ray data, out primitive intersection)
//     - use intersect_any if you only need to know whether there is a hit
//         hit = intersect_any(scene, ray data)
// 4.b. perform closet-point tests
//     - use overlap_first to get the closet element to a point
//       bounded by a radius
//         hit = overlap_first(scene, pos, radius, out primitive)
//     - use overlap_any to check whether there is an element that overlaps
//         a point within a given radius
// 4.c perform shape overlap queries with overlap_shape_bounds
// 4.d can restrict computation to only one shape if desired
// 5. use interpolate_shape_vert to get interpolated vertex values from the
//    intersection data
//    interpolate_vertex(intersection data, shape data, out interpolated val)
// 6. use refit_bvh to recompute the bvh bounds if transforms changed
//    (you should rebuild the bvh for large changes)
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles),
// an array of vertex positions, and an array of vertex radius
// (for points and lines). Use ym::array_view to pass shared data
// without copying.
//
// Implementation notes:
// - using high precision bvh intersection by default
//

//
// COMPILATION:
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
// - v 0.7: [major API change] move to modern C++ interface
// - v 0.6: internal refactoring
// - v 0.5: removal of C interface
// - v 0.4: use of STL containers
// - v 0.3: cleaner interface based on opaque scene
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

#ifndef _YBVH_H_
#define _YBVH_H_

// compilation options
#ifndef YGL_DECLARATION
#define YGL_API inline
#else
#define YGL_API
#endif

#include "yocto_math.h"

#include <functional>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ybvh {

//
// Heuristic strategy for bvh build
//
enum struct htype {
    def = 0,     // default strategy
    equalnum,    // balanced binary tree
    equalsize,   // splitting space binary tree
    sah,         // surface area heuristic (full sweep for accuracy)
    binned_sah,  // surface area heuristic (binned for speed)
    htype_max    // total number of strategies
};

//
// BVH tree node containing its bounds, indices to the BVH arrays of either
// sorted primitives or internal nodes, whether its a leaf or an internal node,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh for more details.
//
// This is not part of the public interface.
//
// Implemenetation Notes:
// - Padded to 32 bytes for cache fiendly access
//
struct _bvhn {
    ym::bbox3f bbox;  // bounding box
    uint32_t start;   // index to the first sorted primitive/node
    uint16_t count;   // number of primitives/nodes
    uint8_t isleaf;   // whether it is a leaf
    uint8_t axis;     // slit axis
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
// This is not part of the public interface.
//
struct _bvh {
    // bvh data
    std::vector<_bvhn> nodes;      // sorted array of internal nodes
    std::vector<int> sorted_prim;  // sorted elements
};

//
// Shape Data as index mesh. Only one element type is supported at one time.
//
struct shape {
    // shape transform --------------------
    ym::mat4f xform;      // shape transform
    ym::mat4f inv_xform;  // shape inverse transforms

    // elements data ----------------------
    ym::array_view<int> points;           // point indices
    ym::array_view<ym::vec2i> lines;      // line indices
    ym::array_view<ym::vec3i> triangles;  // triangle indices

    // vertex data ------------------------
    ym::array_view<ym::vec3f> pos;  // vertex pos
    ym::array_view<float> radius;   // vertex radius

    // [private] bvh data -----------------
    _bvh bvh;  // bvh [private]
};

//
// Scene Data
//
struct scene {
    // scene data -------------------------
    std::vector<shape> shapes;  // shapes

    // bvh private data -------------------
    _bvh bvh;  // bvh [private]
};

//
// Builds a scene/shape BVH.
//
// Parameters:
// - shape/scene: object to build the bvh for
// - heuristic: heristic use to build the bvh (htype::none for default)
//
YGL_API void build_bvh(scene& scene, htype heuristic = htype::def,
                       bool do_shapes = true);
YGL_API void build_bvh(shape& shape, htype heuristic = htype::def);

//
// Refit the bounds of each shape for moving objects. Use this only to avoid
// a rebuild, but note that queries are likely slow if objects move a lot.
// Before calling refit, set the scene data.
//
// Parameters:
// - scene/shape: scene to refit
//
YGL_API void refit_bvh(scene& scene, bool do_shapes = false);
YGL_API void refit_bvh(shape& shape);

//
// BVH intersection.
//
struct point {
    float dist = 0;              // distance
    int sid = -1;                // shape index
    int eid = -1;                // element index
    ym::vec3f euv = ym::zero3f;  // element baricentric coordinates
    bool hit = false;            // hit

    // check whether it was a hit
    operator bool() const { return hit; }
};

//
// Intersect the scene with a ray finding the closest intersection.
//
// Parameters:
// - scene/shape: scene/shape to intersect
// - ray: ray
//
// Return:
// - intersection point
//
YGL_API point intersect_first(const scene& scene, const ym::ray3f& ray);
YGL_API point intersect_first(const shape& shape, const ym::ray3f& ray);

//
// Intersect the scene with a ray finding any intersection (good for shadow
// rays).
//
// Parameters:
// - scene/shape: scene/shape to intersect
// - ray: ray
//
// Return:
// - whether we intersect or not
//
YGL_API bool intersect_any(const scene& scene, const ym::ray3f& ray);
YGL_API bool intersect_any(const shape& shape, const ym::ray3f& ray);

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
YGL_API int overlap_shape_bounds(
    const scene& scene, bool exclude_self,
    std::function<void(const ym::vec2i& v)> overlap_cb);

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
YGL_API int overlap_shape_bounds(const scene& scene, bool exclude_self,
                                 std::vector<ym::vec2i>& overlaps);

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
YGL_API point overlap_first(const scene& scene, const ym::vec3f& pt,
                            float max_dist);
YGL_API point overlap_first(const shape& shape, const ym::vec3f& pt,
                            float max_dist);

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
YGL_API bool overlap_any(const scene& scene, const ym::vec3f& pt,
                         float max_dist, int req_shape = -1);
YGL_API bool overlap_any(const shape& shape, const ym::vec3f& pt,
                         float max_dist);

//
// Interpolates a vertex property from the given intersection data. Uses
// linear interpolation for lines, baricentric for triangles and copies values
// for points.
//
// Parameters:
// - points/lines/triangles: array of vertex indices (only one filled)
// - vert: vertex property array
// - eid: element index
// - euv: element parameters
//
// Returns:
// - interpolated vertex data
//
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<int>& points,
                                 const ym::array_view<ym::vec2i>& lines,
                                 const ym::array_view<ym::vec3i>& triangles,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv);
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<ym::vec2i>& lines,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv);
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<ym::vec3i>& triangles,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv);

// -----------------------------------------------------------------------------
// PERFORMANCE TUNING INTERFACE
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
YGL_API void compute_bvh_stats(const scene& scene, bool include_shapes,
                               int& nprims, int& ninternals, int& nleaves,
                               int& min_depth, int& max_depth,
                               int req_shape = -1);

//
// Logging (not thread safe to avoid affecting speed)
//
#ifdef YGL_BVH_LOG_RAYS
YGL_API void get_ray_log(int& nrays, int& nbbox_inters, int& npoint_inters,
                         int& nline_inters, int& ntriangle_inters);
#endif

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

// handles compilation options
#if (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

#include <algorithm>
#include <cstdio>

namespace ybvh {

// -----------------------------------------------------------------------------
// ELEMENT-WISE BOUNDS
// -----------------------------------------------------------------------------

//
// Point bounds
//
static inline ym::bbox3f _bound_point(const ym::vec3f& p, float r = 0) {
    return ym::make_bbox({p - r, p - r});
}

//
// Line bounds
//
static inline ym::bbox3f _bound_line(const ym::vec3f& v0, const ym::vec3f& v1,
                                     float r0 = 0, float r1 = 0) {
    return ym::make_bbox({v0 - r0, v0 + r0, v1 - r1, v1 + r1});
}

//
// Triangle bounds
//
static inline ym::bbox3f _bound_triangle(const ym::vec3f& v0,
                                         const ym::vec3f& v1,
                                         const ym::vec3f& v2) {
    return ym::make_bbox({v0, v1, v2});
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
static inline bool _intersect_point(const ym::ray3f& ray, const ym::vec3f& p,
                                    float r, float& ray_t, ym::vec3f& euv) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = ym::dot(w, ray.d) / ym::dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp = ray.eval(t);
    auto prp = p - rp;
    if (ym::dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1, 0, 0};

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
static inline bool _intersect_line(const ym::ray3f& ray, const ym::vec3f& v0,
                                   const ym::vec3f& v1, float r0, float r1,
                                   float& ray_t, ym::vec3f& euv) {
    // setup intersection params
    auto u = ray.d;
    auto v = v1 - v0;
    auto w = ray.o - v0;

    // compute values to solve a linear system
    auto a = ym::dot(u, u);
    auto b = ym::dot(u, v);
    auto c = ym::dot(v, v);
    auto d = ym::dot(u, w);
    auto e = ym::dot(v, w);
    auto det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute parameters on both ray and segment
    auto t = (b * e - c * d) / det;
    auto s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // clamp segment param to segment corners
    s = ym::clamp(s, 0, 1);

    // compute segment-segment distance on the closest points
    auto p0 = ray.eval(t);
    auto p1 = ym::ray3f(v0, v1 - v0).eval(s);
    auto p01 = p0 - p1;

    // check with the line radius at the same point
    auto r = r0 * (1 - s) + r1 * s;
    if (ym::dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - s, s, 0};

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
static inline bool _intersect_triangle(const ym::ray3f& ray,
                                       const ym::vec3f& v0, const ym::vec3f& v1,
                                       const ym::vec3f& v2, float& ray_t,
                                       ym::vec3f& euv) {
    // compute triangle edges
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;

    // compute determinant to solve a linear system
    auto pvec = ym::cross(ray.d, edge2);
    auto det = ym::dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.o - v0;
    auto u = ym::dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = ym::cross(tvec, edge1);
    auto v = ym::dot(ray.d, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = ym::dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - u - v, u, v};

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
static inline bool _intersect_check_bbox(const ym::ray3f& ray,
                                         const ym::bbox3f& bbox) {
    // set up convenient pointers for looping over axes
    auto tmin = ray.tmin, tmax = ray.tmax;

    // for each axis, clip intersection against the bounding planes
    for (int i = 0; i < 3; i++) {
        // determine intersection ranges
        auto invd = 1.0f / ray.d[i];
        auto t0 = (bbox[0][i] - ray.o[i]) * invd;
        auto t1 = (bbox[1][i] - ray.o[i]) * invd;
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
static inline const T& _safemin(const T& a, const T& b) {
    return (a < b) ? a : b;
}
template <typename T>
static inline const T& _safemax(const T& a, const T& b) {
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
static inline bool _intersect_check_bbox(const ym::ray3f& ray,
                                         const ym::vec3f& ray_dinv,
                                         const ym::vec3i& ray_dsign,
                                         const ym::bbox3f& bbox) {
    auto txmin = (bbox[ray_dsign[0]][0] - ray.o[0]) * ray_dinv[0];
    auto txmax = (bbox[1 - ray_dsign[0]][0] - ray.o[0]) * ray_dinv[0];
    auto tymin = (bbox[ray_dsign[1]][1] - ray.o[1]) * ray_dinv[1];
    auto tymax = (bbox[1 - ray_dsign[1]][1] - ray.o[1]) * ray_dinv[1];
    auto tzmin = (bbox[ray_dsign[2]][2] - ray.o[2]) * ray_dinv[2];
    auto tzmax = (bbox[1 - ray_dsign[2]][2] - ray.o[2]) * ray_dinv[2];
    auto tmin = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// -----------------------------------------------------------------------------
// ELEMENT-WISE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------

// TODO: documentation
static inline bool _overlap_point(const ym::vec3f& pos, float dist_max,
                                  const ym::vec3f& p, float r, float& dist,
                                  ym::vec3f& euv) {
    auto d2 = ym::distsqr(pos, p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = std::sqrt(d2);
    euv = {1, 0, 0};
    return true;
}

// TODO: documentation
static inline float _closestuv_line(const ym::vec3f& pos, const ym::vec3f& v0,
                                    const ym::vec3f& v1) {
    auto ab = v1 - v0;
    auto d = ym::dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    auto u = ym::dot(pos - v0, ab) / d;
    u = ym::clamp(u, 0, 1);
    return u;
}

// TODO: documentation
static inline bool _overlap_line(const ym::vec3f& pos, float dist_max,
                                 const ym::vec3f& v0, const ym::vec3f& v1,
                                 float r0, float r1, float& dist,
                                 ym::vec3f& euv) {
    auto u = _closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = ym::lerp(v0, v1, u);
    auto r = ym::lerp(r0, r1, u);
    auto d2 = ym::distsqr(pos, p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = std::sqrt(d2);
    euv = {1 - u, u, 0};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
static inline ym::vec2f _closestuv_triangle(const ym::vec3f& pos,
                                            const ym::vec3f& v0,
                                            const ym::vec3f& v1,
                                            const ym::vec3f& v2) {
    auto ab = v1 - v0;
    auto ac = v2 - v0;
    auto ap = pos - v0;

    auto d1 = ym::dot(ab, ap);
    auto d2 = ym::dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return ym::vec2f{0, 0};

    auto bp = pos - v1;
    auto d3 = ym::dot(ab, bp);
    auto d4 = ym::dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return ym::vec2f{1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0))
        return ym::vec2f{d1 / (d1 - d3), 0};

    auto cp = pos - v2;
    auto d5 = ym::dot(ab, cp);
    auto d6 = ym::dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return ym::vec2f{0, 1};

    auto vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0))
        return ym::vec2f{0, d2 / (d2 - d6)};

    auto va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return ym::vec2f{1 - w, w};
    }

    // face case
    auto denom = 1 / (va + vb + vc);
    auto v = vb * denom;
    auto w = vc * denom;
    return ym::vec2f{v, w};
}

// TODO: documentation
static inline bool _overlap_triangle(const ym::vec3f& pos, float dist_max,
                                     const ym::vec3f& v0, const ym::vec3f& v1,
                                     const ym::vec3f& v2, float r0, float r1,
                                     float r2, float& dist, ym::vec3f& euv) {
    auto uv = _closestuv_triangle(pos, v0, v1, v2);
    auto p = ym::blerp(v0, v1, v2, ym::vec3f{1 - uv[0] - uv[1], uv[0], uv[1]});
    auto r = ym::blerp(r0, r1, r2, ym::vec3f{1 - uv[0] - uv[1], uv[0], uv[1]});
    auto dd = ym::distsqr(p, pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = std::sqrt(dd);
    euv = {1 - uv[0] - uv[1], uv[0], uv[1]};
    return true;
}

// TODO: documentation
static inline bool _distance_check_bbox(const ym::vec3f pos, float dist_max,
                                        const ym::vec3f bbox_min,
                                        const ym::vec3f bbox_max) {
    // computing distance
    auto dd = 0.0f;
    // For each axis count any excess distance outside box extents
    for (int i = 0; i < 3; i++) {
        auto v = pos[i];
        if (v < bbox_min[i]) dd += (bbox_min[i] - v) * (bbox_min[i] - v);
        if (v > bbox_max[i]) dd += (v - bbox_max[i]) * (v - bbox_max[i]);
    }

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
static inline bool _overlap_bbox(const ym::vec3f bbox1_min,
                                 const ym::vec3f bbox1_max,
                                 const ym::vec3f bbox2_min,
                                 const ym::vec3f bbox2_max) {
    if (bbox1_max[0] < bbox2_min[0] || bbox1_min[0] > bbox2_max[0])
        return false;
    if (bbox1_max[1] < bbox2_min[1] || bbox1_min[1] > bbox2_max[1])
        return false;
    if (bbox1_max[2] < bbox2_min[2] || bbox1_min[2] > bbox2_max[2])
        return false;
    return true;
}

// -----------------------------------------------------------------------------
// BVH DATA STRUCTURE
// -----------------------------------------------------------------------------

// number of primitives to avoid splitting on
#define YB__BVH_MINPRIMS 4

// -----------------------------------------------------------------------------
// BVH BUILD FUNCTIONS
// -----------------------------------------------------------------------------

//
// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
//
struct _bound_prim {
    ym::bbox3f bbox;       // bounding box
    ym::vec3f center;      // bounding box center (for faster sort)
    int pid;               // primitive id
    float sah_cost_left;   // buffer for sah heuristic costs
    float sah_cost_right;  // buffer for sah heuristic costs
};

//
// Comparison function for each axis
//
struct _bound_prim_comp {
    int axis;
    float middle;

    _bound_prim_comp(int a, float m = 0) : axis(a), middle(m) {}

    bool operator()(const _bound_prim& a, const _bound_prim& b) const {
        return a.center[axis] < b.center[axis];
    }

    bool operator()(const _bound_prim& a) const {
        return a.center[axis] < middle;
    }
};

//
// Given an array sorted_prim of primitives to split between the elements
// start and end, determines the split axis axis, split primitive index mid
// based on the heuristic heuristic. Supports balanced tree (equalnum) and
// Surface-Area Heuristic.
//
static inline bool _partition_prims(ym::array_view<_bound_prim> sorted_prim,
                                    int start, int end, int& axis, int& mid,
                                    htype heuristic) {
    const auto __box_eps = 1e-12f;
#define __bbox_area(r)                                                         \
    (2 * ((r[0] + __box_eps) * (r[1] + __box_eps) +                            \
          (r[1] + __box_eps) * (r[2] + __box_eps) +                            \
          (r[2] + __box_eps) * (r[0] + __box_eps)))

    // init to default values
    axis = 0;
    mid = (start + end) / 2;

    // compute primintive bounds and size
    auto centroid_bbox = ym::invalid_bbox3f;
    for (auto i = start; i < end; i++) centroid_bbox += sorted_prim[i].center;
    auto centroid_size = ym::diagonal(centroid_bbox);

    // check if it is not possible to split
    if (centroid_size == ym::zero3f) return false;

    // split along largest
    auto largest_axis = ym::max_element(centroid_size);

    // check heuristic
    switch (heuristic) {
        // balanced tree split: find the largest axis of the bounding box
        // and split along this one right in the middle
        case htype::equalnum: {
            axis = largest_axis;
            mid = (start + end) / 2;
            std::nth_element(sorted_prim.data() + start,
                             sorted_prim.data() + mid, sorted_prim.data() + end,
                             _bound_prim_comp(largest_axis));
        } break;
        // split the space in the middle along the largest axis
        case htype::def:
        case htype::equalsize: {
            axis = largest_axis;
            mid =
                (int)(std::partition(
                          sorted_prim.data() + start, sorted_prim.data() + end,
                          _bound_prim_comp(
                              largest_axis,
                              ym::center(centroid_bbox)[largest_axis])) -
                      sorted_prim.data());
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case htype::sah: {
            auto min_cost = HUGE_VALF;
            auto count = end - start;
            for (auto a = 0; a < 3; a++) {
                std::sort(sorted_prim.data() + start, sorted_prim.data() + end,
                          _bound_prim_comp(a));
                auto sbbox = ym::invalid_bbox3f;
                // to avoid an O(n^2) computation, use sweaps to compute the
                // cost,
                // first smallest to largest, then largest to smallest
                for (auto i = 0; i < count; i++) {
                    sbbox += sorted_prim[start + i].bbox;
                    auto sbbox_size = ym::diagonal(sbbox);
                    sorted_prim[start + i].sah_cost_left =
                        __bbox_area(sbbox_size);
                    sorted_prim[start + i].sah_cost_left *= i + 1;
                }
                // the other sweep
                sbbox = ym::invalid_bbox3f;
                for (auto i = 0; i < count; i++) {
                    sbbox += sorted_prim[end - 1 - i].bbox;
                    auto sbbox_size = ym::diagonal(sbbox);
                    sorted_prim[end - 1 - i].sah_cost_right =
                        __bbox_area(sbbox_size);
                    sorted_prim[end - 1 - i].sah_cost_right *= i + 1;
                }
                // find the minimum cost
                for (auto i = start + 2; i <= end - 2; i++) {
                    auto cost = sorted_prim[i - 1].sah_cost_left +
                                sorted_prim[i].sah_cost_right;
                    if (min_cost > cost) {
                        min_cost = cost;
                        axis = a;
                        mid = i;
                    }
                }
            }
            std::nth_element(sorted_prim.data() + start,
                             sorted_prim.data() + (mid),
                             sorted_prim.data() + end, _bound_prim_comp(axis));
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case htype::binned_sah: {
            // allocate bins
            const auto nbins = 16;
            ym::bbox3f bins_bbox[nbins];
            int bins_count[nbins];
            for (int b = 0; b < nbins; b++) {
                bins_bbox[b] = ym::invalid_bbox3f;
                // bins_bbox[b] = centroid_bbox;
                // bins_bbox[b].min[largest_axis] += b *
                // centroid_size[largest_axis] / nbins;
                // bins_bbox[b].max[largest_axis] -= (nbins - 1 - b) *
                // centroid_size[largest_axis] / nbins;
                bins_count[b] = 0;
            }
            for (int i = start; i < end; i++) {
                auto b = (int)(nbins * (sorted_prim[i].center[largest_axis] -
                                        centroid_bbox[0][largest_axis]) /
                               centroid_size[largest_axis]);
                b = ym::clamp(b, 0, nbins - 1);
                bins_count[b] += 1;
                bins_bbox[b] += sorted_prim[i].bbox;
            }
            float min_cost = HUGE_VALF;
            int bin_idx = -1;
            for (int b = 1; b < nbins; b++) {
                if (!bins_count[b]) continue;
                auto left_bbox = ym::invalid_bbox3f,
                     right_bbox = ym::invalid_bbox3f;
                auto left_count = 0, right_count = 0;
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
                auto left_bbox_size = ym::diagonal(left_bbox);
                auto right_bbox_size = ym::diagonal(right_bbox);
                auto cost = __bbox_area(left_bbox_size) * left_count +
                            __bbox_area(right_bbox_size) * right_count;
                if (min_cost > cost) {
                    min_cost = cost;
                    bin_idx = b;
                }
            }
            assert(bin_idx >= 0);
            axis = largest_axis;
            mid = start;
            for (int b = 0; b < bin_idx; b++) {
                mid += bins_count[b];
            }
            assert(axis >= 0 && mid > 0);
            std::nth_element(
                sorted_prim.data() + start, sorted_prim.data() + (mid),
                sorted_prim.data() + end, _bound_prim_comp(largest_axis));
        } break;
        default: assert(false); break;
    }
    assert(axis >= 0 && mid > 0);
    assert(mid > start && mid < end);
    return true;
}

//
// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
//
static inline void _make_node(_bvhn& node, std::vector<_bvhn>& nodes,
                              ym::array_view<_bound_prim> sorted_prims,
                              int start, int end, htype heuristic) {
    // compute node bounds
    node.bbox = ym::invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += sorted_prims[i].bbox;

    // decide whether to create a leaf
    if (end - start <= YB__BVH_MINPRIMS) {
        // makes a leaf node
        node.isleaf = true;
        node.start = start;
        node.count = end - start;
    } else {
        // choose the split axis and position
        int axis, mid;
        auto split =
            _partition_prims(sorted_prims, start, end, axis, mid, heuristic);
        if (!split) {
            // we failed to split for some reasons
            node.isleaf = true;
            node.start = start;
            node.count = end - start;
        } else {
            assert(mid > start && mid < end);
            // makes an internal node
            node.isleaf = false;
            // perform the splits by preallocating the child nodes and recurring
            node.axis = axis;
            node.start = (int)nodes.size();
            node.count = 2;
            nodes.push_back({});
            nodes.push_back({});
            // build child nodes
            _make_node(nodes[node.start], nodes, sorted_prims, start, mid,
                       heuristic);
            _make_node(nodes[node.start + 1], nodes, sorted_prims, mid, end,
                       heuristic);
        }
    }
}

//
// Build a BVH from a set of primitives.
//
template <typename T, typename bbox_func>
YGL_API void _build_bvh(_bvh& bvh, const ym::array_view<T>& prims,
                        htype heuristic, const bbox_func& bbox) {
    // create bounded primitives used in BVH build
    auto bound_prims = std::vector<_bound_prim>();
    for (auto& f : prims) {
        auto prim = _bound_prim();
        prim.pid = (int)bound_prims.size();
        prim.bbox = bbox(f);
        prim.center = ym::center(prim.bbox);
        bound_prims.push_back(prim);
    }

    // clear bvh
    bvh.nodes.clear();
    bvh.sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh.nodes.reserve(bound_prims.size() * 2);

    // start recursive splitting
    bvh.nodes.push_back({});
    _make_node(bvh.nodes[0], bvh.nodes, bound_prims, 0, (int)bound_prims.size(),
               heuristic);

    // shrink back
    bvh.nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh.sorted_prim.resize(bound_prims.size());
    for (int i = 0; i < bound_prims.size(); i++) {
        bvh.sorted_prim[i] = bound_prims[i].pid;
    }
    bvh.sorted_prim.shrink_to_fit();
}

//
// Build a shape BVH. Public function whose interface is described above.
//
// TODO: move to simpler interface in loops
YGL_API void build_bvh(shape& shape, htype heuristic) {
    // compute bounds for each primitive type
    if (!shape.points.empty()) {
        // point bounds are computed as little spheres
        _build_bvh(shape.bvh, shape.points, heuristic, [&shape](const int f) {
            auto r = (!shape.radius.empty()) ? shape.radius[f] : 0;
            return _bound_point(shape.pos[f], r);
        });
    } else if (!shape.lines.empty()) {
        // line bounds are computed as thick rods
        _build_bvh(
            shape.bvh, shape.lines, heuristic, [&shape](const ym::vec2i& f) {
                auto r0 = (!shape.radius.empty()) ? shape.radius[f[0]] : 0;
                auto r1 = (!shape.radius.empty()) ? shape.radius[f[1]] : 0;
                return _bound_line(shape.pos[f[0]], shape.pos[f[1]], r0, r1);
            });
    } else if (!shape.triangles.empty()) {
        // triangle bounds are computed by including their vertices
        _build_bvh(shape.bvh, shape.triangles, heuristic,
                   [&shape](const ym::vec3i& f) {
                       return _bound_triangle(shape.pos[f[0]], shape.pos[f[1]],
                                              shape.pos[f[2]]);
                   });
    } else {
        assert(false);
    }
}

//
// Build a scene BVH. Public function whose interface is described above.
//
YGL_API void build_bvh(scene& scene, htype heuristic, bool do_shapes) {
    // recursively build shape bvhs
    if (do_shapes) {
        for (auto& shape : scene.shapes) build_bvh(shape, heuristic);
    }

    // build scene bvh
    _build_bvh(scene.bvh, ym::array_view<shape>(scene.shapes), heuristic,
               [](const shape& shape) {
                   return ym::transform_bbox(shape.xform,
                                             shape.bvh.nodes[0].bbox);
               });
}

//
// Recursively recomputes the node bounds for a scene bvh
//
template <typename T, typename bbox_func>
static inline void _refit_bvh(_bvh& bvh, const std::vector<T>& prims,
                              _bvhn& node, const bbox_func& bbox) {
    if (node.isleaf) {
        node.bbox = ym::invalid_bbox3f;
        for (auto i = node.start; i < node.start + node.count; i++) {
            auto idx = bvh.sorted_prim[i];
            node.bbox += bbox(prims[idx]);
        }
    } else {
        node.bbox = ym::invalid_bbox3f;
        for (auto i = node.start; i < node.start + node.count; i++) {
            _refit_bvh(bvh, prims, bvh.nodes[i], bbox);
            node.bbox += bvh.nodes[i].bbox;
        }
    }
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
YGL_API void refit_bvh(shape& shape) {}

//
// Refits a scene BVH. Public function whose interface is described above.
//
YGL_API void refit_bvh(scene& scene, bool do_shapes) {
    // do shapes
    if (do_shapes) {
        for (auto& shape : scene.shapes) refit_bvh(shape);
    }

    // recompute bvh bounds
    _refit_bvh(
        scene.bvh, scene.shapes, scene.bvh.nodes[0], [](const shape& shape) {
            return ym::transform_bbox(shape.xform, shape.bvh.nodes[0].bbox);
        });
}

// -----------------------------------------------------------------------------
// BVH INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Logging global variables
//
#ifdef YGL_BVH_LOG_RAYS
int _log_nrays = 0;
int _log_nbbox_inters = 0;
int _log_npoint_inters = 0;
int _log_nline_inters = 0;
int _log_ntriangle_inters = 0;

//
// Logging information
//
YGL_API void get_ray_log(int& nrays, int& nbbox_inters, int& npoint_inters,
                         int& nline_inters, int& ntriangle_inters) {
    nrays = _log_nrays;
    nbbox_inters = _log_nbbox_inters;
    npoint_inters = _log_npoint_inters;
    nline_inters = _log_nline_inters;
    ntriangle_inters = _log_ntriangle_inters;
}

#endif

//
// Intersect ray with a bvh. Similar to the generic public function whose
// interface is described above. See intersect_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of recursive
// calls; this follows general conventions and stragely makes the code shorter
// - The walk is simplified for first hit by obeserving that if we update
// the ray_tmax limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequence farthest iterations will be
// rejected in the tmax tests
//
template <typename Isec_func>
static inline bool _intersect_bvh(const _bvh& bvh, const ym::ray3f& ray_,
                                  bool early_exit, float& ray_t, int& id,
                                  const Isec_func& intersect_prim) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // copy ray
    auto ray = ray_;

    // shared variables
    auto hit = false;

    // prepare ray for fast queries
    auto ray_dinv = ym::vec3f{1, 1, 1} / ray.d;
    auto ray_dsign =
        ym::vec3i{(ray_dinv[0] < 0) ? 1 : 0, (ray_dinv[1] < 0) ? 1 : 0,
                  (ray_dinv[2] < 0) ? 1 : 0};
    auto ray_reverse = ym::vec<bool, 4>{(bool)ray_dsign[0], (bool)ray_dsign[1],
                                        (bool)ray_dsign[2], false};

    // walking stack
    while (node_cur && !(early_exit && hit)) {
        // exit if needed
        if (early_exit && hit) return hit;

        // grab node
        auto& node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!_intersect_check_bbox(ray, ray_dinv, ray_dsign, node.bbox))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // for internal nodes, attempts to proceed along the
            // split axis from smallest to largest nodes
            if (ray_reverse[node.axis]) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            } else {
                for (auto i = node.count - 1; i >= 0; i--) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh.sorted_prim[node.start + i];
                if (intersect_prim(idx, ray, ray_t)) {
                    hit = true;
                    id = idx;
                    ray.tmax = ray_t;
                }
            }
        }
    }

    return hit;
}

//
// Shape intersection
//
static inline bool _intersect_shape(const shape& shape, const ym::ray3f& ray,
                                    bool early_exit, float& ray_t, int& eid,
                                    ym::vec3f& euv) {
    auto tray = transform_ray(shape.inv_xform, ray);
    if (!shape.triangles.empty()) {
        return _intersect_bvh(shape.bvh, tray, early_exit, ray_t, eid,
                              [&](int idx, const ym::ray3f& ray, float& ray_t) {
                                  auto f = shape.triangles[idx];
                                  return _intersect_triangle(
                                      ray, shape.pos[f[0]], shape.pos[f[1]],
                                      shape.pos[f[2]], ray_t, euv);
                              });
    } else if (!shape.lines.empty()) {
        return _intersect_bvh(
            shape.bvh, tray, early_exit, ray_t, eid,
            [&](int idx, const ym::ray3f& ray, float& ray_t) {
                auto f = shape.lines[idx];
                auto r0 = (!shape.radius.empty()) ? shape.radius[f[0]] : 0;
                auto r1 = (!shape.radius.empty()) ? shape.radius[f[1]] : 0;
                return _intersect_line(ray, shape.pos[f[0]], shape.pos[f[1]],
                                       r0, r1, ray_t, euv);
            });
    } else if (!shape.points.empty()) {
        return _intersect_bvh(
            shape.bvh, tray, early_exit, ray_t, eid,
            [&](int idx, const ym::ray3f& ray, float& ray_t) {
                auto f = shape.points[idx];
                auto r = (!shape.radius.empty()) ? shape.radius[f] : 0;
                return _intersect_point(ray, shape.pos[f], r, ray_t, euv);
            });
    } else {
        assert(false);
        return false;
    }
}

//
// Scene intersection
//
static inline bool _intersect_scene(const scene& scene, const ym::ray3f& ray,
                                    bool early_exit, float& ray_t, int& sid,
                                    int& eid, ym::vec3f& euv) {
    return _intersect_bvh(scene.bvh, ray, early_exit, ray_t, sid,
                          [&](int idx, ym::ray3f& ray, float& ray_t) {
                              return _intersect_shape(scene.shapes[idx], ray,
                                                      early_exit, ray_t, eid,
                                                      euv);
                          });
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API point intersect_first(const scene& scene, const ym::ray3f& ray) {
    auto isec = point();
    isec.hit = _intersect_scene(scene, ray, false, isec.dist, isec.sid,
                                isec.eid, isec.euv);
    return isec;
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API point intersect_first(const shape& shape, const ym::ray3f& ray) {
    auto isec = point();
    isec.hit =
        _intersect_shape(shape, ray, false, isec.dist, isec.eid, isec.euv);
    return isec;
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool intersect_any(const scene& scene, const ym::ray3f& ray) {
    auto isec = point();
    return _intersect_scene(scene, ray, true, isec.dist, isec.sid, isec.eid,
                            isec.euv);
}

//
// Find closest ray-scene intersection. A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool intersect_any(const shape& shape, const ym::ray3f& ray) {
    auto isec = point();
    return _intersect_shape(shape, ray, true, isec.dist, isec.eid, isec.euv);
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP
// -----------------------------------------------------------------------------

//
// Finds the closest element with a bvh. Similar to the generic public function
// whose
// interface is described above. See nightbor_bvh for parameter docs.
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
template <typename Overlap_func>
static inline bool _overlap_bvh(const _bvh& bvh, const ym::vec3f& pt,
                                float max_dist, bool early_exit, float& dist,
                                int& id, const Overlap_func& overlap_prim) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // walking stack
    while (node_cur && !(early_exit && hit)) {
        // exit if needed
        if (early_exit && hit) return hit;

        // grab node
        auto& node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!_distance_check_bbox(pt, max_dist, node.bbox[0], node.bbox[1]))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // internal node
            for (auto idx = node.start; idx < node.start + node.count; idx++) {
                node_stack[node_cur++] = idx;
                assert(node_cur < 64);
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh.sorted_prim[node.start + i];
                if (overlap_prim(idx, pt, max_dist, dist)) {
                    hit = true;
                    id = idx;
                    max_dist = dist;
                }
            }
        }
    }

    return hit;
}

//
// Shape overlap
//
static inline bool _overlap_bvh(const shape& shape, const ym::vec3f& pt,
                                float max_dist, bool early_exit, float& dist,
                                int& eid, ym::vec3f& euv) {
    auto tpt = transform_point(shape.inv_xform, pt);
    if (!shape.triangles.empty()) {
        return _overlap_bvh(
            shape.bvh, tpt, max_dist, early_exit, dist, eid,
            [&](int idx, const ym::vec3f& pt, float max_dist, float& dist) {
                auto f = shape.triangles[idx];
                auto r0 = (!shape.radius.empty()) ? shape.radius[f[0]] : 0;
                auto r1 = (!shape.radius.empty()) ? shape.radius[f[1]] : 0;
                auto r2 = (!shape.radius.empty()) ? shape.radius[f[2]] : 0;
                return _overlap_triangle(pt, max_dist, shape.pos[f[0]],
                                         shape.pos[f[1]], shape.pos[f[2]], r0,
                                         r1, r2, dist, euv);
            });
    } else if (!shape.lines.empty()) {
        return _overlap_bvh(
            shape.bvh, tpt, max_dist, early_exit, dist, eid,
            [&](int idx, const ym::vec3f& pt, float max_dist, float& dist) {
                auto f = shape.lines[idx];
                auto r0 = (!shape.radius.empty()) ? shape.radius[f[0]] : 0;
                auto r1 = (!shape.radius.empty()) ? shape.radius[f[1]] : 0;
                return _overlap_line(pt, max_dist, shape.pos[f[0]],
                                     shape.pos[f[1]], r0, r1, dist, euv);
            });
    } else if (!shape.points.empty()) {
        return _overlap_bvh(
            shape.bvh, tpt, max_dist, early_exit, dist, eid,
            [&](int idx, const ym::vec3f& pt, float max_dist, float& dist) {
                auto f = shape.points[idx];
                auto r = (!shape.radius.empty()) ? shape.radius[f] : 0;
                return _overlap_point(pt, max_dist, shape.pos[f], r, dist, euv);
            });
    } else {
        assert(false);
        return false;
    }
}

//
// Scene overlap
//
static inline bool _overlap_bvh(const scene& scene, const ym::vec3f& pt,
                                float max_dist, bool early_exit, float& dist,
                                int& sid, int& eid, ym::vec3f& euv) {
    return _overlap_bvh(
        scene.bvh, pt, max_dist, early_exit, dist, sid,
        [&](int idx, const ym::vec3f& pt, float max_dist, float& dist) {
            return _overlap_bvh(scene.shapes[idx], pt, max_dist, early_exit,
                                dist, eid, euv);
        });
}

//
// Find the closest element to a given point within a maximum distance.
// A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API point overlap_first(const scene& scene, const ym::vec3f& pt,
                            float max_dist) {
    auto overlap = point();
    overlap.hit = _overlap_bvh(scene, pt, max_dist, false, overlap.dist,
                               overlap.sid, overlap.eid, overlap.euv);
    return overlap;
}
YGL_API point overlap_first(const shape& shape, const ym::vec3f& pt,
                            float max_dist) {
    auto overlap = point();
    overlap.hit = _overlap_bvh(shape, pt, max_dist, false, overlap.dist,
                               overlap.eid, overlap.euv);
    return overlap;
}

//
// Find if any element overlaps given point within a maximum distance.
// A convenient wrapper for the above func.
// Public function whose interface is described above.
//
YGL_API bool overlap_any(const scene& scene, const ym::vec3f& pt,
                         float max_dist, int req_shape) {
    auto overlap = point();
    return _overlap_bvh(scene, pt, max_dist, true, overlap.dist, overlap.sid,
                        overlap.eid, overlap.euv);
}
YGL_API bool overlap_any(const shape& shape, const ym::vec3f& pt,
                         float max_dist, int req_shape) {
    auto overlap = point();
    return _overlap_bvh(shape, pt, max_dist, true, overlap.dist, overlap.eid,
                        overlap.euv);
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
template <typename Overlap_func>
YGL_API int _overlap_bounds(const _bvh& bvh, bool exclude_self,
                            const Overlap_func& overlap_prims) {
    // node stack
    ym::vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // shared variables
    auto hits = 0;

    // walking stack
    while (node_cur) {
        // grab node
        auto node_idx = node_stack[--node_cur];
        auto node1 = &bvh.nodes[node_idx[0]];
        auto node2 = &bvh.nodes[node_idx[1]];

        // intersect bbox
        if (!_overlap_bbox(node1->bbox[0], node1->bbox[1], node2->bbox[0],
                           node2->bbox[1]))
            continue;

        // check for leaves
        if (node1->isleaf && node2->isleaf) {
            // collide primitives
            for (auto i1 = node1->start; i1 < node1->start + node1->count;
                 i1++) {
                for (auto i2 = node2->start; i2 < node2->start + node2->count;
                     i2++) {
                    auto idx =
                        ym::vec2i{bvh.sorted_prim[i1], bvh.sorted_prim[i2]};
                    if (exclude_self && idx[0] == idx[1]) continue;
                    hits += overlap_prims(idx);
                }
            }
        } else {
            // descend
            if (node1->isleaf) {
                for (auto idx = node2->start; idx < node2->start + node2->count;
                     idx++) {
                    node_stack[node_cur++] = {node_idx[0], (int)idx};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx = node1->start; idx < node1->start + node1->count;
                     idx++) {
                    node_stack[node_cur++] = {(int)idx, node_idx[1]};
                    assert(node_cur < 128);
                }
            }
        }
    }

    // done
    return hits;
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
YGL_API int overlap_shape_bounds(
    const scene& scene, bool exclude_self,
    std::function<void(const ym::vec2i& v)> overlap_cb) {
    return _overlap_bounds(scene.bvh, exclude_self, [&](const ym::vec2i& idx) {
        auto& shape1 = scene.shapes[idx[0]];
        auto& shape2 = scene.shapes[idx[1]];
        auto bbox1 = ym::transform_bbox(shape1.xform, shape1.bvh.nodes[0].bbox);
        auto bbox2 = ym::transform_bbox(shape2.xform, shape2.bvh.nodes[0].bbox);
        if (!_overlap_bbox(bbox1[0], bbox1[1], bbox2[0], bbox2[1])) return 0;
        overlap_cb(idx);
        return 1;
    });
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
YGL_API int overlap_shape_bounds(const scene& scene, bool exclude_self,
                                 std::vector<ym::vec2i>& overlaps) {
    overlaps.clear();
    return overlap_shape_bounds(
        scene, exclude_self,
        [&overlaps](const ym::vec2i& v) { overlaps.push_back(v); });
}

// -----------------------------------------------------------------------------
// VERTEX PROPERTY INTERPOLATION
// -----------------------------------------------------------------------------

//
// Interpolate vertex properties at the intersection. Public API.
//
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<int>& points,
                                 const ym::array_view<ym::vec2i>& lines,
                                 const ym::array_view<ym::vec3i>& triangles,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv) {
    if (!points.empty())
        return vert[points[eid]];
    else if (!lines.empty())
        return vert[lines[eid][0]] * (1 - euv[0]) +
               vert[lines[eid][1]] * euv[1];
    else if (!triangles.empty())
        return vert[triangles[eid][0]] * (1 - euv[0] - euv[1]) +
               vert[triangles[eid][1]] * euv[0] +
               vert[triangles[eid][2]] * euv[1];
    assert(false);
    return T();
}

//
// Interpolate vertex properties at the intersection. Public API.
//
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<ym::vec2i>& lines,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv) {
    return vert[lines[eid][0]] * (1 - euv[0]) + vert[lines[eid][1]] * euv[1];
}

//
// Interpolate vertex properties at the intersection. Public API.
//
template <typename T>
YGL_API T interpolate_shape_vert(const ym::array_view<ym::vec3i>& triangles,
                                 const ym::array_view<T>& vert, int eid,
                                 const ym::vec2f& euv) {
    return vert[triangles[eid][0]] * (1 - euv[0] - euv[1]) +
           vert[triangles[eid][1]] * euv[0] + vert[triangles[eid][2]] * euv[1];
}

// -----------------------------------------------------------------------------
// STATISTICS FOR DEBUGGING (probably not helpful to all)
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
YGL_API void _compute_bvh_stats(const scene& scene, int shape_id,
                                bool include_shapes, int depth, int& nprims,
                                int& ninternals, int& nleaves, int& min_depth,
                                int& max_depth) {
    // get bvh
    auto bvh = (shape_id >= 0) ? scene.shapes[shape_id].bvh : scene.bvh;

    // node stack
    ym::vec2i node_stack[128];  // node and depth
    auto node_cur = 0;
    node_stack[node_cur++] = ym::vec2i{0, depth};

    // walk the stack
    while (node_cur) {
        // get node and depth
        auto node_depth = node_stack[--node_cur];
        auto& node = bvh.nodes[node_depth[0]];

        // update stats
        if (!node.isleaf) {
            ninternals += 1;
            for (auto i = 0; i < node.count; i++) {
                node_stack[node_cur++] =
                    ym::vec2i(node.start + i, node_depth[1] + 1);
            }
        } else if (shape_id >= 0) {
            nleaves += 1;
            nprims += node.count;
            min_depth = ym::min(min_depth, node_depth[1]);
            max_depth = ym::max(max_depth, node_depth[1]);
        } else {
            if (include_shapes) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh.sorted_prim[node.start + i];
                    _compute_bvh_stats(scene, idx, true, node_depth[1] + 1,
                                       nprims, ninternals, nleaves, min_depth,
                                       max_depth);
                }
            } else {
                nleaves += 1;
                nprims += node.count;
                min_depth = ym::min(min_depth, node_depth[1]);
                max_depth = ym::max(max_depth, node_depth[1]);
            }
        }
    }
}

//
// Compute BVH stats.
//
YGL_API void compute_bvh_stats(const scene& scene, bool include_shapes,
                               int& nprims, int& ninternals, int& nleaves,
                               int& min_depth, int& max_depth, int req_shape) {
    // init out variables
    nprims = 0;
    ninternals = 0;
    nleaves = 0;
    min_depth = 10000;
    max_depth = -1;

    _compute_bvh_stats(scene, req_shape, include_shapes, 0, nprims, ninternals,
                       nleaves, min_depth, max_depth);
}

}  // namespace

#endif

#endif
