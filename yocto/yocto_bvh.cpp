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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF YOCTO_BVH
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"
#include "yocto_math.h"

#include <algorithm>
#include <cstdio>
#include <unordered_map>

namespace ybvh {

// -----------------------------------------------------------------------------
// ELEMENT-WISE BOUNDS
// -----------------------------------------------------------------------------

//
// Point bounds
//
static inline ym::bbox3f _bound_point(const ym::vec3f& p, float r = 0) {
    return ym::bbox3f{p - ym::vec3f{r, r, r}, p + ym::vec3f{r, r, r}};
}

//
// Line bounds
//
static inline ym::bbox3f _bound_line(
    const ym::vec3f& v0, const ym::vec3f& v1, float r0 = 0, float r1 = 0) {
    return ym::make_bbox(
        {v0 - ym::vec3f{r0, r0, r0}, v0 + ym::vec3f{r0, r0, r0},
            v1 - ym::vec3f{r1, r1, r1}, v1 + ym::vec3f{r1, r1, r1}});
}

//
// Triangle bounds
//
static inline ym::bbox3f _bound_triangle(
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2) {
    return ym::make_bbox({v0, v1, v2});
}

//
// Triangle bounds
//
static inline ym::bbox3f _bound_tetrahedron(const ym::vec3f& v0,
    const ym::vec3f& v1, const ym::vec3f& v2, const ym::vec3f& v3) {
    return ym::make_bbox({v0, v1, v2, v3});
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
static inline bool _intersect_point(
    const ym::ray3f& ray, const ym::vec3f& p, float r, float& ray_t) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = ym::dot(w, ray.d) / ym::dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp = ym::eval(ray, t);
    auto prp = p - rp;
    if (ym::dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;

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
    const ym::vec3f& v1, float r0, float r1, float& ray_t, ym::vec2f& euv) {
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
    s = ym::clamp(s, (float)0, (float)1);

    // compute segment-segment distance on the closest points
    auto p0 = ym::eval(ray, t);
    auto p1 = ym::eval(ym::ray3f{v0, v1 - v0}, s);
    auto p01 = p0 - p1;

    // check with the line radius at the same point
    auto r = r0 * (1 - s) + r1 * s;
    if (ym::dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - s, s};

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
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2, float& ray_t,
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
// Intersect a ray with a tetrahedron. Note that we consider only intersection
// wiht the tetrahedra surface and discount intersction with the interior.
//
// Parameters:
// - ray: ray to intersect with
// - v0, v1, v2: triangle vertices
//
// Out Parameters:
// - ray_t: ray parameter at the intersection point
// - euv: baricentric coordinates of the intersection
//
// Returns:
// - whether the intersection occurred
//
// TODO: check order
// TODO: uv
//
static inline bool _intersect_tetrahedron(const ym::ray3f& ray_,
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2,
    const ym::vec3f& v3, float& ray_t, ym::vec4f& euv) {
    // check intersction for each face
    auto hit = false;
    auto ray = ray_;
    auto tuv = ym::zero3f;
    if (_intersect_triangle(ray, v0, v1, v2, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (_intersect_triangle(ray, v0, v1, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (_intersect_triangle(ray, v0, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (_intersect_triangle(ray, v1, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }

    return hit;
}

//
// Intersect a ray with a axis-aligned bounding box
//
// Parameters:
// - ray: ray to intersect with
// - bbox: bounding box min/max bounds
//
// Returns:
// - whether the intersection occurred
//
static inline bool _intersect_check_bbox(
    const ym::ray3f& ray, const ym::bbox3f& bbox) {
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
    const ym::vec3f& ray_dinv, const ym::vec3i& ray_dsign,
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
    const ym::vec3f& p, float r, float& dist) {
    auto d2 = ym::distsqr(pos, p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(d2);
    return true;
}

// TODO: documentation
static inline float _closestuv_line(
    const ym::vec3f& pos, const ym::vec3f& v0, const ym::vec3f& v1) {
    auto ab = v1 - v0;
    auto d = ym::dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    auto u = ym::dot(pos - v0, ab) / d;
    u = ym::clamp(u, (float)0, (float)1);
    return u;
}

// TODO: documentation
static inline bool _overlap_line(const ym::vec3f& pos, float dist_max,
    const ym::vec3f& v0, const ym::vec3f& v1, float r0, float r1, float& dist,
    ym::vec2f& euv) {
    auto u = _closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = ym::lerp(v0, v1, u);
    auto r = ym::lerp(r0, r1, u);
    auto d2 = ym::distsqr(pos, p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = std::sqrt(d2);
    euv = {1 - u, u};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
static inline ym::vec2f _closestuv_triangle(const ym::vec3f& pos,
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2) {
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
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2, float r0,
    float r1, float r2, float& dist, ym::vec3f& euv) {
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
static inline bool _overlap_tetrahedron(const ym::vec3f& pos,
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2,
    const ym::vec3f& v3, ym::vec4f& euv) {
    auto vol = ym::dot(v3 - v0, ym::cross(v3 - v1, v3 - v0));
    if (vol == 0) return false;
    auto u = ym::dot(v3 - v0, ym::cross(v3 - v1, v3 - v0)) / vol;
    if (u < 0 || u > 1) return false;
    auto v = ym::dot(v3 - v0, ym::cross(v3 - v1, v3 - v0)) / vol;
    if (v < 0 || v > 1 || u + v > 1) return false;
    auto w = ym::dot(v3 - v0, ym::cross(v3 - v1, v3 - v0)) / vol;
    if (w < 0 || w > 1 || u + v + w > 1) return false;
    euv = {u, v, w, 1 - u - v - w};
    return true;
}

// TODO: documentation
static inline bool _overlap_tetrahedron(const ym::vec3f& pos, float dist_max,
    const ym::vec3f& v0, const ym::vec3f& v1, const ym::vec3f& v2,
    const ym::vec3f& v3, float r0, float r1, float r2, float r3, float& dist,
    ym::vec4f& euv) {
    // check interior
    if (_overlap_tetrahedron(pos, v0, v1, v2, v3, euv)) {
        dist = 0;
        return true;
    }

    // check faces
    auto hit = false;
    auto tuv = ym::zero3f;
    if (_overlap_triangle(pos, dist_max, v0, v1, v2, r0, r1, r2, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (_overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (_overlap_triangle(pos, dist_max, v0, v2, v3, r0, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (_overlap_triangle(pos, dist_max, v1, v2, v3, r1, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }

    return hit;
}

// TODO: documentation
static inline bool _distance_check_bbox(
    const ym::vec3f& pos, float dist_max, const ym::bbox3f& bbox) {
    // computing distance
    auto dd = 0.0f;

    // For each axis count any excess distance outside box extents
    for (int i = 0; i < 3; i++) {
        auto v = pos[i];
        if (v < bbox[0][i]) dd += (bbox[0][i] - v) * (bbox[0][i] - v);
        if (v > bbox[1][i]) dd += (v - bbox[1][i]) * (v - bbox[1][i]);
    }

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
static inline bool _overlap_bbox(
    const ym::bbox3f& bbox1, const ym::bbox3f& bbox2) {
    if (bbox1[1][0] < bbox2[0][0] || bbox1[0][0] > bbox2[1][0]) return false;
    if (bbox1[1][1] < bbox2[0][1] || bbox1[0][1] > bbox2[1][1]) return false;
    if (bbox1[1][2] < bbox2[0][2] || bbox1[0][2] > bbox2[1][2]) return false;
    return true;
}

// TODO: doc
// from "Real-Time Collision Detection" by Christer Ericson, Sect. 4.4.1
static inline bool _overlap_bbox(const ym::bbox3f& bbox1,
    const ym::bbox3f& bbox2, const ym::frame3f& frame1,
    const ym::frame3f& frame2) {
#define __YBVH_EPSILON__ 1e-5f
    // compute centered frames and extents
    auto cframe1 = ym::make_frame(
        ym::rot(frame1), transform_point(frame1, ym::center(bbox1)));
    auto cframe2 = ym::make_frame(
        ym::rot(frame2), transform_point(frame2, ym::center(bbox2)));
    auto ext1 = ym::diagonal(bbox1) / 2.0f, ext2 = ym::diagonal(bbox2) / 2.0f;

    // compute frame from 2 to 1
    auto cframe2to1 = ym::inverse(cframe1) * cframe2;

    // split frame components and move to row-major
    auto rot = ym::transpose(ym::rot(cframe2to1));
    auto t = ym::pos(cframe2to1);

    // Compute common subexpressions. Add in an epsilon term to
    // counteract arithmetic errors when two edges are parallel and
    // their cross product is (near) null (see text for details)
    ym::mat3f absrot;
    auto parallel_axis = false;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            absrot[i][j] = std::abs(rot[i][j]) + __YBVH_EPSILON__;
            if (absrot[i][j] > 1) parallel_axis = true;
        }

    // Test axes L = A0, L = A1, L = A2
    for (int i = 0; i < 3; i++) {
        if (std::abs(t[i]) > ext1[i] + ym::dot(ext2, absrot[i])) return false;
    }

    // Test axes L = B0, L = B1, L = B2
    for (int i = 0; i < 3; i++) {
        auto ra = ext1[0] * absrot[0][i] + ext1[1] * absrot[1][i] +
                  ext1[2] * absrot[2][i];
        auto rb = ext2[i];
        if (std::abs(t[0] * rot[0][i] + t[1] * rot[1][i] + t[2] * rot[2][i]) >
            ra + rb)
            return false;
    }

    // if axis were nearly parallel, can exit here
    if (parallel_axis) return true;

    // Test axis L = A0 x B0
    auto ra = ext1[1] * absrot[2][0] + ext1[2] * absrot[1][0];
    auto rb = ext2[1] * absrot[0][2] + ext2[2] * absrot[0][1];
    if (std::abs(t[2] * rot[1][0] - t[1] * rot[2][0]) > ra + rb) return false;

    // Test axis L = A0 x B1
    ra = ext1[1] * absrot[2][1] + ext1[2] * absrot[1][1];
    rb = ext2[0] * absrot[0][2] + ext2[2] * absrot[0][0];
    if (std::abs(t[2] * rot[1][1] - t[1] * rot[2][1]) > ra + rb) return false;

    // Test axis L = A0 x B2
    ra = ext1[1] * absrot[2][2] + ext1[2] * absrot[1][2];
    rb = ext2[0] * absrot[0][1] + ext2[1] * absrot[0][0];
    if (std::abs(t[2] * rot[1][2] - t[1] * rot[2][2]) > ra + rb) return false;

    // Test axis L = A1 x B0
    ra = ext1[0] * absrot[2][0] + ext1[2] * absrot[0][0];
    rb = ext2[1] * absrot[1][2] + ext2[2] * absrot[1][1];
    if (std::abs(t[0] * rot[2][0] - t[2] * rot[0][0]) > ra + rb) return false;

    // Test axis L = A1 x B1
    ra = ext1[0] * absrot[2][1] + ext1[2] * absrot[0][1];
    rb = ext2[0] * absrot[1][2] + ext2[2] * absrot[1][0];
    if (std::abs(t[0] * rot[2][1] - t[2] * rot[0][1]) > ra + rb) return false;

    // Test axis L = A1 x B2
    ra = ext1[0] * absrot[2][2] + ext1[2] * absrot[0][2];
    rb = ext2[0] * absrot[1][1] + ext2[1] * absrot[1][0];
    if (std::abs(t[0] * rot[2][2] - t[2] * rot[0][2]) > ra + rb) return false;

    // Test axis L = A2 x B0
    ra = ext1[0] * absrot[1][0] + ext1[1] * absrot[0][0];
    rb = ext2[1] * absrot[2][2] + ext2[2] * absrot[2][1];
    if (std::abs(t[1] * rot[0][0] - t[0] * rot[1][0]) > ra + rb) return false;

    // Test axis L = A2 x B1
    ra = ext1[0] * absrot[1][1] + ext1[1] * absrot[0][1];
    rb = ext2[0] * absrot[2][2] + ext2[2] * absrot[2][0];
    if (std::abs(t[1] * rot[0][1] - t[0] * rot[1][1]) > ra + rb) return false;

    // Test axis L = A2 x B2
    ra = ext1[0] * absrot[1][2] + ext1[1] * absrot[0][2];
    rb = ext2[0] * absrot[2][1] + ext2[1] * absrot[2][0];
    if (std::abs(t[1] * rot[0][2] - t[0] * rot[1][2]) > ra + rb) return false;

    // Since no separating axis is found, the OBBs must be intersecting
    return true;
#undef __YBVH_EPSILON__
}

// this is only a conservative test!
// TODO: rename to something more clear
// TODO: doc
static inline bool _overlap_bbox_conservative(const ym::bbox3f& bbox1,
    const ym::bbox3f& bbox2, const ym::frame3f& frame1,
    const ym::frame3f& frame2) {
    return _overlap_bbox(
        bbox1, ym::transform_bbox(ym::inverse(frame1) * frame2, bbox2));
}

// -----------------------------------------------------------------------------
// BVH DATA STRUCTURE
// -----------------------------------------------------------------------------

// number of primitives to avoid splitting on
#define YBVH__MINPRIMS 4

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
struct bvhn {
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
struct bvh {
    // bvh data
    std::vector<bvhn> nodes;       // sorted array of internal nodes
    std::vector<int> sorted_prim;  // sorted elements
};

//
// Shape Data as index mesh. Only one element type is supported at one time.
//
struct shape {
    // shape id ---------------------------
    int sid = -1;  // shape id

    // shape transform --------------------
    ym::frame3f frame;  // shape transform

    // elements data ----------------------
    int nelems = 0;                       // number of elements
    const int* point = nullptr;           // point indices
    const ym::vec2i* line = nullptr;      // line indices
    const ym::vec3i* triangle = nullptr;  // triangle indices
    const ym::vec4i* tetra = nullptr;     // tetrahedra indices

    // vertex data ------------------------
    int nverts = 0;                  // number of vertices
    const ym::vec3f* pos = nullptr;  // vertex pos
    const float* radius = nullptr;   // vertex radius

    // [private] bvh data -----------------
    bvh* _bvh = nullptr;  // bvh [private]

    // [private] methods ------------------
    float rad(int i) const { return (radius) ? radius[i] : 0; }
    ym::bbox3f local_bbox() const { return _bvh->nodes[0].bbox; }
    ym::bbox3f world_bbox() const {
        return ym::transform_bbox(frame, _bvh->nodes[0].bbox);
    }

    // destructor
    ~shape();
};

//
// Scene Data
//
struct scene {
    // scene data -------------------------
    std::vector<shape*> shapes;  // shapes

    // bvh private data -------------------
    bvh* _bvh = nullptr;  // bvh [private]

    // destructor
    ~scene();
};

//
// Init scene.
//
YBVH_API scene* make_scene(int nshapes) {
    auto scn = new scene();
    scn->shapes.resize(nshapes);
    for (auto& shp : scn->shapes) shp = new shape();
    return scn;
}

///
/// Free scene.
///
YBVH_API void free_scene(scene* scn) {
    if (scn) delete scn;
}

//
// Destructor
//
YBVH_API scene::~scene() {
    for (auto& shp : shapes) delete shp;
    if (_bvh) delete _bvh;
}

YBVH_API shape::~shape() {
    if (_bvh) delete _bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_point_shape(scene* scn, int sid, const float3x4& frame,
    int npoints, const int* point, int nverts, const float3* pos,
    const float* radius) {
    scn->shapes[sid]->sid = sid;
    scn->shapes[sid]->frame = frame;
    scn->shapes[sid]->nelems = npoints;
    scn->shapes[sid]->point = (const int*)point;
    scn->shapes[sid]->line = nullptr;
    scn->shapes[sid]->triangle = nullptr;
    scn->shapes[sid]->tetra = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = (const ym::vec3f*)pos;
    scn->shapes[sid]->radius = radius;
    if (scn->shapes[sid]->_bvh) delete scn->shapes[sid]->_bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_line_shape(scene* scn, int sid, const float3x4& frame,
    int nlines, const int2* lines, int nverts, const float3* pos,
    const float* radius) {
    scn->shapes[sid]->sid = sid;
    scn->shapes[sid]->frame = frame;
    scn->shapes[sid]->nelems = nlines;
    scn->shapes[sid]->point = nullptr;
    scn->shapes[sid]->line = (const ym::vec2i*)lines;
    scn->shapes[sid]->triangle = nullptr;
    scn->shapes[sid]->tetra = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = (const ym::vec3f*)pos;
    scn->shapes[sid]->radius = radius;
    if (scn->shapes[sid]->_bvh) delete scn->shapes[sid]->_bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_triangle_shape(scene* scn, int sid, const float3x4& frame,
    int ntriangles, const int3* triangles, int nverts, const float3* pos,
    const float* radius) {
    scn->shapes[sid]->sid = sid;
    scn->shapes[sid]->frame = frame;
    scn->shapes[sid]->nelems = ntriangles;
    scn->shapes[sid]->point = nullptr;
    scn->shapes[sid]->line = nullptr;
    scn->shapes[sid]->triangle = (const ym::vec3i*)triangles;
    scn->shapes[sid]->tetra = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = (const ym::vec3f*)pos;
    scn->shapes[sid]->radius = radius;
    if (scn->shapes[sid]->_bvh) delete scn->shapes[sid]->_bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_tetra_shape(scene* scn, int sid, const float3x4& frame,
    int ntetra, const int4* tetra, int nverts, const float3* pos,
    const float* radius) {
    scn->shapes[sid]->sid = sid;
    scn->shapes[sid]->frame = frame;
    scn->shapes[sid]->nelems = ntetra;
    scn->shapes[sid]->point = nullptr;
    scn->shapes[sid]->line = nullptr;
    scn->shapes[sid]->triangle = nullptr;
    scn->shapes[sid]->tetra = (const ym::vec4i*)tetra;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = (const ym::vec3f*)pos;
    scn->shapes[sid]->radius = radius;
    if (scn->shapes[sid]->_bvh) delete scn->shapes[sid]->_bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_point_shape(scene* scn, int sid, const float3x4& frame,
    int nverts, const float3* pos, const float* radius) {
    scn->shapes[sid]->sid = sid;
    scn->shapes[sid]->frame = frame;
    scn->shapes[sid]->nelems = nverts;
    scn->shapes[sid]->point = nullptr;
    scn->shapes[sid]->line = nullptr;
    scn->shapes[sid]->triangle = nullptr;
    scn->shapes[sid]->tetra = nullptr;
    scn->shapes[sid]->nverts = nverts;
    scn->shapes[sid]->pos = (const ym::vec3f*)pos;
    scn->shapes[sid]->radius = radius;
    if (scn->shapes[sid]->_bvh) delete scn->shapes[sid]->_bvh;
}

//
// Set shape. Public API.
//
YBVH_API void set_shape_frame(scene* scn, int sid, const float3x4& frame) {
    scn->shapes[sid]->frame = frame;
}

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
static inline bool _partition_prims(_bound_prim* sorted_prim, int start,
    int end, int& axis, int& mid, heuristic_type htype) {
    // internal function
    auto bbox_area = [](auto r) {
        const auto __box_eps = 1e-12f;
        return (2 * ((r[0] + __box_eps) * (r[1] + __box_eps) +
                        (r[1] + __box_eps) * (r[2] + __box_eps) +
                        (r[2] + __box_eps) * (r[0] + __box_eps)));
    };

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
    auto largest_axis = ym::max_element_idx(centroid_size);

    // check heuristic
    switch (htype) {
        // balanced tree split: find the largest axis of the bounding box
        // and split along this one right in the middle
        case heuristic_type::equalnum: {
            axis = largest_axis;
            mid = (start + end) / 2;
            std::nth_element(sorted_prim + start, sorted_prim + mid,
                sorted_prim + end, _bound_prim_comp(largest_axis));
        } break;
        // split the space in the middle along the largest axis
        case heuristic_type::def:
        case heuristic_type::equalsize: {
            axis = largest_axis;
            mid = (int)(std::partition(sorted_prim + start, sorted_prim + end,
                            _bound_prim_comp(largest_axis,
                                ym::center(centroid_bbox)[largest_axis])) -
                        sorted_prim);
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case heuristic_type::sah: {
            auto min_cost = HUGE_VALF;
            auto count = end - start;
            for (auto a = 0; a < 3; a++) {
                std::sort(sorted_prim + start, sorted_prim + end,
                    _bound_prim_comp(a));
                auto sbbox = ym::invalid_bbox3f;
                // to avoid an O(n^2) computation, use sweaps to compute the
                // cost,
                // first smallest to largest, then largest to smallest
                for (auto i = 0; i < count; i++) {
                    sbbox += sorted_prim[start + i].bbox;
                    auto sbbox_size = ym::diagonal(sbbox);
                    sorted_prim[start + i].sah_cost_left =
                        bbox_area(sbbox_size);
                    sorted_prim[start + i].sah_cost_left *= i + 1;
                }
                // the other sweep
                sbbox = ym::invalid_bbox3f;
                for (auto i = 0; i < count; i++) {
                    sbbox += sorted_prim[end - 1 - i].bbox;
                    auto sbbox_size = ym::diagonal(sbbox);
                    sorted_prim[end - 1 - i].sah_cost_right =
                        bbox_area(sbbox_size);
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
            std::nth_element(sorted_prim + start, sorted_prim + mid,
                sorted_prim + end, _bound_prim_comp(axis));
        } break;
        // surface area heuristic: estimate the cost of splitting
        // along each of the axis and pick the one with best expected
        // performance
        case heuristic_type::binned_sah: {
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
                auto cost = bbox_area(left_bbox_size) * left_count +
                            bbox_area(right_bbox_size) * right_count;
                if (min_cost > cost) {
                    min_cost = cost;
                    bin_idx = b;
                }
            }
            assert(bin_idx >= 0);
            axis = largest_axis;
            mid = start;
            for (int b = 0; b < bin_idx; b++) { mid += bins_count[b]; }
            assert(axis >= 0 && mid > 0);
            std::nth_element(sorted_prim + start, sorted_prim + mid,
                sorted_prim + end, _bound_prim_comp(largest_axis));
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
static inline void _make_node(bvhn& node, std::vector<bvhn>& nodes,
    _bound_prim* sorted_prims, int start, int end, heuristic_type htype) {
    // compute node bounds
    node.bbox = ym::invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += sorted_prims[i].bbox;

    // decide whether to create a leaf
    if (end - start <= YBVH__MINPRIMS) {
        // makes a leaf node
        node.isleaf = true;
        node.start = start;
        node.count = end - start;
    } else {
        // choose the split axis and position
        int axis, mid;
        auto split =
            _partition_prims(sorted_prims, start, end, axis, mid, htype);
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
            nodes.emplace_back();
            nodes.emplace_back();
            // build child nodes
            _make_node(
                nodes[node.start], nodes, sorted_prims, start, mid, htype);
            _make_node(
                nodes[node.start + 1], nodes, sorted_prims, mid, end, htype);
        }
    }
}

//
// Build a BVH from a set of primitives.
//
YBVH_API void _build_bvh(
    bvh* bvh, int nprims, _bound_prim* bound_prims, heuristic_type htype) {
    // clear bvh
    bvh->nodes.clear();
    bvh->sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh->nodes.reserve(nprims * 2);

    // start recursive splitting
    bvh->nodes.emplace_back();
    _make_node(bvh->nodes[0], bvh->nodes, bound_prims, 0, nprims, htype);

    // shrink back
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(nprims);
    for (int i = 0; i < nprims; i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }
}

//
// Gets the bbox of a shape element
//
static inline ym::bbox3f _bound_elem(const shape* shp, int eid) {
    if (shp->point) {
        auto f = shp->point[eid];
        return _bound_point(shp->pos[f], shp->rad(f));
    } else if (shp->line) {
        auto f = shp->line[eid];
        return _bound_line(
            shp->pos[f[0]], shp->pos[f[1]], shp->rad(f[0]), shp->rad(f[1]));
    } else if (shp->triangle) {
        auto f = shp->triangle[eid];
        return _bound_triangle(shp->pos[f[0]], shp->pos[f[1]], shp->pos[f[2]]);
    } else if (shp->tetra) {
        auto f = shp->tetra[eid];
        return _bound_tetrahedron(
            shp->pos[f[0]], shp->pos[f[1]], shp->pos[f[2]], shp->pos[f[3]]);
    } else {
        return _bound_point(shp->pos[eid], shp->rad(eid));
    }
}

//
// Build a shape BVH. Public function whose interface is described above.
//
YBVH_API void build_bvh(shape* shp, heuristic_type htype) {
    // create bounded primitives used in BVH build
    auto bound_prims = std::vector<_bound_prim>(shp->nelems);
    for (auto i = 0; i < bound_prims.size(); i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = _bound_elem(shp, i);
        bound_prims[i].center = ym::center(bound_prims[i].bbox);
    }

    // tree bvh
    if (!shp->_bvh) shp->_bvh = new bvh();
    return _build_bvh(
        shp->_bvh, (int)bound_prims.size(), bound_prims.data(), htype);
}

//
// Build a scene BVH. Public function whose interface is described above.
//
YBVH_API void build_bvh(scene* scn, heuristic_type htype, bool do_shapes) {
    // do shapes
    if (do_shapes) {
        for (auto shp : scn->shapes) build_bvh(shp, htype);
    }

    // create bounded primitives used in BVH build
    auto bound_prims = std::vector<_bound_prim>(scn->shapes.size());
    for (auto i = 0; i < bound_prims.size(); i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = scn->shapes[i]->world_bbox();
        bound_prims[i].center = ym::center(bound_prims[i].bbox);
    }

    // tree bvh
    if (!scn->_bvh) scn->_bvh = new bvh();
    _build_bvh(scn->_bvh, (int)bound_prims.size(), bound_prims.data(), htype);
}

//
// Recursively recomputes the node bounds for a shape bvh
//
static inline void _refit_bvh(scene* scn, int sid, int nodeid, bool do_shapes) {
    // get shape and bvh
    auto shp = (sid < 0) ? nullptr : scn->shapes[sid];
    auto bvh = (!shp) ? scn->_bvh : shp->_bvh;

    // refit
    auto node = &bvh->nodes[nodeid];
    node->bbox = ym::invalid_bbox3f;
    if (node->isleaf) {
        if (!shp) {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                if (do_shapes) _refit_bvh(scn, idx, 0, false);
                node->bbox += scn->shapes[idx]->world_bbox();
            }
        } else {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                node->bbox += _bound_elem(shp, idx);
            }
        }
    } else {
        for (auto i = 0; i < node->count; i++) {
            auto idx = node->start + i;
            _refit_bvh(scn, sid, idx, do_shapes);
            node->bbox += bvh->nodes[i].bbox;
        }
    }
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
YBVH_API void refit_bvh(scene* scn, int sid) {
    return _refit_bvh(scn, sid, 0, false);
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
YBVH_API void refit_bvh(scene* scn, bool do_shapes) {
    // recompute bvh bounds
    _refit_bvh(scn, -1, 0, do_shapes);
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
YBVH_API void get_ray_log(int& nrays, int& nbbox_inters, int& npoint_inters,
    int& nline_inters, int& ntriangle_inters) {
    nrays = _log_nrays;
    nbbox_inters = _log_nbbox_inters;
    npoint_inters = _log_npoint_inters;
    nline_inters = _log_nline_inters;
    ntriangle_inters = _log_ntriangle_inters;
}

#endif

//
// Intersect a shape element. Used in the templated function below.
//
static inline point _intersect_elem(
    const shape* shp, int eid, const ym::ray3f& ray, bool early_exit) {
    // initialize point
    auto pt = point();

    // switch over shape type
    if (shp->triangle) {
        auto f = shp->triangle[eid];
        if (!_intersect_triangle(ray, shp->pos[f[0]], shp->pos[f[1]],
                shp->pos[f[2]], pt.dist, (ym::vec3f&)pt.euv))
            return {};
        pt.euv = {pt.euv[0], pt.euv[1], pt.euv[2], 0};
    } else if (shp->line) {
        assert(shp->radius);
        auto f = shp->line[eid];
        if (!_intersect_line(ray, shp->pos[f[0]], shp->pos[f[1]],
                shp->radius[f[0]], shp->radius[f[1]], pt.dist,
                (ym::vec2f&)pt.euv))
            return point{};
        pt.euv = {pt.euv[0], pt.euv[1], 0, 0};
    } else if (shp->point) {
        assert(shp->radius);
        auto f = shp->point[eid];
        if (!_intersect_point(ray, shp->pos[f], shp->radius[f], pt.dist))
            return {};
        pt.euv = {1, 0, 0, 0};
    } else if (shp->tetra) {
        auto f = shp->tetra[eid];
        if (!_intersect_tetrahedron(ray, shp->pos[f[0]], shp->pos[f[1]],
                shp->pos[f[2]], shp->pos[f[3]], pt.dist, (ym::vec4f&)pt.euv))
            return {};
    } else {
        assert(shp->radius);
        if (!_intersect_point(ray, shp->pos[eid], shp->radius[eid], pt.dist))
            return {};
        pt.euv = {1, 0, 0, 0};
    }

    // finilize point
    pt.eid = eid;
    pt.sid = shp->sid;
    return pt;
}

//
// Intersect ray with a bvh-> Similar to the generic public function whose
// interface is described above. See intersect_ray for parameter docs.
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
static inline point _intersect_ray(
    const scene* scn, int sid, const ym::ray3f& ray_, bool early_exit) {
    // get shape and bvh
    auto shp = (sid < 0) ? nullptr : scn->shapes[sid];
    auto bvh = (!shp) ? scn->_bvh : shp->_bvh;

    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // copy ray and transform it if necessary
    auto ray = (!shp) ? ray_ : ym::transform_ray_inverse(shp->frame, ray_);

    // shared variables
    auto pt = point();

    // prepare ray for fast queries
    auto ray_dinv = ym::vec3f{1, 1, 1} / ray.d;
    auto ray_dsign = ym::vec3i{(ray_dinv[0] < 0) ? 1 : 0,
        (ray_dinv[1] < 0) ? 1 : 0, (ray_dinv[2] < 0) ? 1 : 0};
    auto ray_reverse = ym::vec<bool, 4>{
        (bool)ray_dsign[0], (bool)ray_dsign[1], (bool)ray_dsign[2], false};

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

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
            if (!shp) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto pp = point();
                    if ((pp = _intersect_ray(scn, idx, ray, early_exit))) {
                        if (early_exit) return pp;
                        pt = pp;
                        ray.tmax = pt.dist;
                    }
                }
            } else {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto pp = point();
                    if ((pp = _intersect_elem(shp, idx, ray, early_exit))) {
                        if (early_exit) return pp;
                        pt = pp;
                        ray.tmax = pt.dist;
                    }
                }
            }
        }
    }

    return pt;
}

//
// Shape intersection
//
YBVH_API point intersect_ray(const scene* scn, int sid, const float3& ray_o,
    const float3& ray_d, float ray_tmin, float ray_tmax, bool early_exit) {
    return _intersect_ray(
        scn, sid, {ray_o, ray_d, ray_tmin, ray_tmax}, early_exit);
}

//
// Scene intersection
//
YBVH_API point intersect_ray(const scene* scn, const float3& ray_o,
    const float3& ray_d, float ray_tmin, float ray_tmax, bool early_exit) {
    return _intersect_ray(
        scn, -1, {ray_o, ray_d, ray_tmin, ray_tmax}, early_exit);
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP
// -----------------------------------------------------------------------------

//
// Find point overlap for shape elements.
//
static inline point _overlap_elem(const shape* shp, int eid,
    const ym::vec3f& pos, float max_dist, bool early_exit) {
    // initialize point
    auto pt = point();

    // switch over elemenet type
    if (shp->triangle) {
        auto f = shp->triangle[eid];
        if (!_overlap_triangle(pos, max_dist, shp->pos[f[0]], shp->pos[f[1]],
                shp->pos[f[2]], shp->rad(f[0]), shp->rad(f[1]), shp->rad(f[2]),
                pt.dist, (ym::vec3f&)pt.euv))
            return {};
        pt.euv = {pt.euv[0], pt.euv[1], pt.euv[2], 0};
    } else if (shp->line) {
        auto f = shp->line[eid];
        if (!_overlap_line(pos, max_dist, shp->pos[f[0]], shp->pos[f[1]],
                shp->rad(f[0]), shp->rad(f[1]), pt.dist, (ym::vec2f&)pt.euv))
            return {};
        pt.euv = {pt.euv[0], pt.euv[1], 0, 0};
    } else if (shp->point) {
        auto f = shp->point[eid];
        if (!_overlap_point(pos, max_dist, shp->pos[f], shp->rad(f), pt.dist))
            return {};
        pt.euv = {1, 0, 0, 0};
    } else if (shp->tetra) {
        auto f = shp->tetra[eid];
        if (!_overlap_tetrahedron(pos, max_dist, shp->pos[f[0]], shp->pos[f[1]],
                shp->pos[f[2]], shp->pos[f[3]], shp->rad(f[0]), shp->rad(f[1]),
                shp->rad(f[2]), shp->rad(f[3]), pt.dist, (ym::vec4f&)pt.euv))
            return {};
    } else {
        if (!_overlap_point(
                pos, max_dist, shp->pos[eid], shp->rad(eid), pt.dist))
            return {};
        pt.euv = {1, 0, 0, 0};
    };

    // wrap up point
    pt.eid = eid;
    pt.sid = shp->sid;
    return pt;
}

//
// Finds the closest element with a bvh-> Similar to the generic public function
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
static inline point _overlap_point(const scene* scn, int sid,
    const ym::vec3f& pos_, float max_dist, bool early_exit) {
    // get shape and bvh
    auto shp = (sid < 0) ? nullptr : scn->shapes[sid];
    auto bvh = (!shp) ? scn->_bvh : shp->_bvh;

    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // get point
    auto pos = (!shp) ? pos_ : transform_point_inverse(shp->frame, pos_);

    // shared variables
    auto pt = point();

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!_distance_check_bbox(pos, max_dist, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // internal node
            for (auto idx = node.start; idx < node.start + node.count; idx++) {
                node_stack[node_cur++] = idx;
                assert(node_cur < 64);
            }
        } else {
            if (!shp) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto pp = point();
                    if ((pp = _overlap_point(
                             scn, idx, pos, max_dist, early_exit))) {
                        if (early_exit) return pp;
                        pt = pp;
                        max_dist = pt.dist;
                    }
                }
            } else {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto pp = point();
                    if ((pp = _overlap_elem(
                             shp, idx, pos, max_dist, early_exit))) {
                        if (early_exit) return pp;
                        pt = pp;
                        max_dist = pt.dist;
                    }
                }
            }
        }
    }

    return pt;
}

//
// Shape overlap
//
YBVH_API point overlap_point(const scene* scn, int sid, const float3& pos,
    float max_dist, bool early_exit) {
    return _overlap_point(scn, sid, pos, max_dist, early_exit);
}

//
// Scene overlap
//
YBVH_API point overlap_point(
    const scene* scn, const float3& pos, float max_dist, bool early_exit) {
    return _overlap_point(scn, -1, pos, max_dist, early_exit);
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP FOR INTERNAL ELEMENTS
// -----------------------------------------------------------------------------

//
// Shape elements overlap
//
// TODO: avoid duplicate elements
//
static inline void _overlap_elem(const shape* shp1, const shape* shp2, int idx1,
    int idx2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, int2>>* overlaps,
    std::unordered_map<int, int>* closest) {
    // prepare point
    ym::vec4i verts;
    if (shp2->triangle) {
        verts = {shp2->triangle[idx2][0], shp2->triangle[idx2][1],
            shp2->triangle[idx2][2], -1};
    } else if (shp2->line) {
        verts = {shp2->line[idx2][0], shp2->line[idx2][1], -1, -1};
    } else if (shp2->point) {
        verts = {shp2->point[idx2], -1, -1, -1};
    } else if (shp2->tetra) {
        verts = shp2->tetra[idx2];
    } else {
        verts = {idx2, -1, -1, -1};
    }

    // loop over vertices
    for (auto vid : verts) {
        if (vid < 0) continue;

        // transform point
        auto pos = transform_point_inverse(
            shp1->frame, transform_point(shp2->frame, shp2->pos[vid]));
        auto rad = shp2->rad(vid) + radius;
        auto pt = point();

        // switch over shape type
        if (shp1->triangle) {
            auto f = shp1->triangle[idx1];
            if (!_overlap_triangle(pos, rad, shp1->pos[f[0]], shp1->pos[f[1]],
                    shp1->pos[f[2]], shp1->rad(f[0]), shp1->rad(f[1]),
                    shp1->rad(f[2]), pt.dist, (ym::vec3f&)pt.euv))
                return;
            pt.euv = {pt.euv[0], pt.euv[1], pt.euv[2], 0};
        } else if (shp1->line) {
            auto f = shp1->line[idx1];
            if (!_overlap_line(pos, rad, shp1->pos[f[0]], shp1->pos[f[1]],
                    shp1->rad(f[0]), shp1->rad(f[1]), pt.dist,
                    (ym::vec2f&)pt.euv))
                return;
            pt.euv = {pt.euv[0], pt.euv[1], 0, 0};
        } else if (shp1->point) {
            auto f = shp1->point[idx1];
            if (!_overlap_point(pos, rad, shp1->pos[f], shp1->rad(f), pt.dist))
                return;
            pt.euv = {1, 0, 0, 0};
        } else if (shp1->tetra) {
            auto f = shp1->tetra[idx1];
            if (!_overlap_tetrahedron(pos, rad, shp1->pos[f[0]],
                    shp1->pos[f[1]], shp1->pos[f[2]], shp1->pos[f[3]],
                    shp1->rad(f[0]), shp1->rad(f[1]), shp1->rad(f[2]),
                    shp1->rad(f[3]), pt.dist, (ym::vec4f&)pt.euv))
                return;
        } else {
            if (!_overlap_point(
                    pos, rad, shp1->pos[idx1], shp1->rad(idx1), pt.dist))
                return;
            pt.euv = {1, 0, 0, 0};
        }

        // wrap up
        pt.sid = shp1->sid;
        pt.eid = idx1;
        if (first_only) {
            if (closest->find(vid) == closest->end()) {
                overlaps->push_back({pt, {shp2->sid, vid}});
                (*closest)[vid] = (int)overlaps->size() - 1;
            } else {
                auto& overlap = (*overlaps)[(*closest)[vid]];
                if (overlap.first.dist > pt.dist) {
                    overlap = {pt, {shp2->sid, vid}};
                }
            }
        } else {
            overlaps->push_back({pt, {shp2->sid, vid}});
        }
    }
}

//
// Finds the closest element for all pairs of vertices.
// Similar to the generic public function whose
// interface is described above. See nightbor_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of
// recursive
// calls; this follows general conventions and stragely makes the code
// shorter
// - The walk is simplified for first hit by obeserving that if we update
// the dist_max limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequent farthest iterations will be
// rejected in the tmax tests
//
static inline void _overlap_verts(const scene* scn1, const scene* scn2,
    int sid1, int sid2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, int2>>* overlaps,
    std::unordered_map<int, int>* closest) {
    // get shape and bvh
    auto shp1 = (sid1 < 0) ? nullptr : scn1->shapes[sid1];
    auto bvh1 = (!shp1) ? scn1->_bvh : shp1->_bvh;
    auto shp2 = (sid2 < 0) ? nullptr : scn2->shapes[sid2];
    auto bvh2 = (!shp2) ? scn2->_bvh : shp2->_bvh;

    // node stack
    ym::vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // get frames
    auto frame1 = (!shp1) ? ym::identity_frame3f : shp1->frame;
    auto frame2 = (!shp2) ? ym::identity_frame3f : shp2->frame;

    // check if a trasformed test is needed
    auto xformed = frame1 != frame2;

    // walking stack
    while (node_cur) {
        // grab node
        const auto node_idx = node_stack[--node_cur];
        const auto node1 = bvh1->nodes[node_idx[0]];
        const auto node2 = bvh2->nodes[node_idx[1]];

        // intersect bbox
        if (xformed) {
            if (!_overlap_bbox(node1.bbox,
                    {node2.bbox[0] - ym::vec3f{radius, radius, radius},
                        node2.bbox[1] + ym::vec3f{radius, radius, radius}},
                    frame1, frame2))
                continue;
        } else {
            if (!_overlap_bbox(node1.bbox,
                    {node2.bbox[0] - ym::vec3f{radius, radius, radius},
                        node2.bbox[1] + ym::vec3f{radius, radius, radius}}))
                continue;
        }

        // check for leaves
        if (node1.isleaf && node2.isleaf) {
            if (!shp1) {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count;
                     i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (exclude_self && idx1 == idx2) continue;
                        _overlap_verts(scn1, scn2, idx1, idx2, exclude_self,
                            radius, first_only, overlaps, closest);
                    }
                }
            } else {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count;
                     i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (exclude_self && idx1 == idx2) continue;
                        _overlap_elem(shp1, shp2, idx1, idx2, exclude_self,
                            radius, first_only, overlaps, closest);
                    }
                }
            }
        } else {
            // descend
            if (node1.isleaf) {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    node_stack[node_cur++] = {node_idx[0], (int)idx2};
                    assert(node_cur < 128);
                }
            } else if (node2.isleaf) {
                for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                     idx1++) {
                    node_stack[node_cur++] = {(int)idx1, node_idx[1]};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    for (auto idx1 = node1.start;
                         idx1 < node1.start + node1.count; idx1++) {
                        node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        assert(node_cur < 128);
                    }
                }
            }
        }
    }
}

//
// Find the list of overlaps between shapes.
// Public function whose interface is described above.
//
YBVH_API void overlap_verts(const scene* scn1, const scene* scn2, int sid1,
    int sid2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, int2>>* overlaps) {
    std::unordered_map<int, int> closest;
    _overlap_verts(
        scn1, scn2, sid1, sid2, false, radius, first_only, overlaps, &closest);
}

//
// Find the list of overlaps between scenes.
// Public function whose interface is described above.
//
YBVH_API void overlap_verts(const scene* scn1, const scene* scn2,
    bool exclude_self, float radius, bool first_only,
    std::vector<int2>* soverlaps,
    std::vector<std::pair<point, int2>>* overlaps) {
    overlap_shape_bounds(scn1, scn2, false, false, exclude_self, soverlaps);
    for (auto sh : *soverlaps) {
        overlap_verts(scn1, scn2, sh[0], sh[1], exclude_self, radius,
            first_only, overlaps);
    }
}

//
// Finds the overlap between shape bounds.
// Similat interface as the public function.
//
static inline void _overlap_shape_bounds(const scene* scn1, const scene* scn2,
    bool conservative, bool skip_duplicates, bool skip_self,
    std::vector<int2>* overlaps) {
    // get bvhs
    auto bvh1 = scn1->_bvh;
    auto bvh2 = scn2->_bvh;

    // node stack
    ym::vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto node_idx = node_stack[--node_cur];
        const auto node1 = bvh1->nodes[node_idx[0]];
        const auto node2 = bvh2->nodes[node_idx[1]];

        // intersect bbox
        if (!_overlap_bbox(node1.bbox, node2.bbox)) continue;

        // check for leaves
        if (node1.isleaf && node2.isleaf) {
            // collide primitives
            for (auto i1 = node1.start; i1 < node1.start + node1.count; i1++) {
                for (auto i2 = node2.start; i2 < node2.start + node2.count;
                     i2++) {
                    auto idx1 = bvh1->sorted_prim[i1];
                    auto idx2 = bvh2->sorted_prim[i2];
                    auto shp1 = scn1->shapes[idx1];
                    auto shp2 = scn2->shapes[idx2];
                    if (skip_duplicates && shp1->sid > shp2->sid) continue;
                    if (skip_self && shp1->sid == shp2->sid) continue;
                    if (conservative) {
                        if (_overlap_bbox(
                                transform_bbox(shp1->frame, shp1->local_bbox()),
                                transform_bbox(
                                    shp2->frame, shp2->local_bbox())))
                            overlaps->push_back({shp1->sid, shp2->sid});
                    } else {
                        if (_overlap_bbox(shp1->local_bbox(),
                                shp2->local_bbox(), shp1->frame, shp2->frame))
                            overlaps->push_back({shp1->sid, shp2->sid});
                    }
                }
            }
        } else {
            // descend
            if (node1.isleaf) {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    node_stack[node_cur++] = {node_idx[0], (int)idx2};
                    assert(node_cur < 128);
                }
            } else if (node2.isleaf) {
                for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                     idx1++) {
                    node_stack[node_cur++] = {(int)idx1, node_idx[1]};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    for (auto idx1 = node1.start;
                         idx1 < node1.start + node1.count; idx1++) {
                        node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        assert(node_cur < 128);
                    }
                }
            }
        }
    }
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
YBVH_API void overlap_shape_bounds(const scene* scn1, const scene* scn2,
    bool conservative, bool skip_duplicates, bool skip_self,
    std::vector<int2>* overlaps) {
    overlaps->clear();
    _overlap_shape_bounds(
        scn1, scn2, conservative, skip_duplicates, skip_self, overlaps);
}

// -----------------------------------------------------------------------------
// STATISTICS FOR DEBUGGING (probably not helpful to all)
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
YBVH_API void _compute_bvh_stats(const scene* scn, int shape_id,
    bool include_shapes, int depth, int& nprims, int& ninternals, int& nleaves,
    int& min_depth, int& max_depth) {
    // get bvh
    auto bvh = (shape_id >= 0) ? scn->shapes[shape_id]->_bvh : scn->_bvh;

    // node stack
    ym::vec2i node_stack[128];  // node and depth
    auto node_cur = 0;
    node_stack[node_cur++] = ym::vec2i{0, depth};

    // walk the stack
    while (node_cur) {
        // get node and depth
        auto node_depth = node_stack[--node_cur];
        auto node = &bvh->nodes[node_depth[0]];

        // update stats
        if (!node->isleaf) {
            ninternals += 1;
            for (auto i = 0; i < node->count; i++) {
                node_stack[node_cur++] =
                    ym::vec2i{(int)(node->start + i), node_depth[1] + 1};
            }
        } else if (shape_id >= 0) {
            nleaves += 1;
            nprims += node->count;
            min_depth = ym::min(min_depth, node_depth[1]);
            max_depth = ym::max(max_depth, node_depth[1]);
        } else {
            if (include_shapes) {
                for (auto i = 0; i < node->count; i++) {
                    auto idx = bvh->sorted_prim[node->start + i];
                    _compute_bvh_stats(scn, idx, true, node_depth[1] + 1,
                        nprims, ninternals, nleaves, min_depth, max_depth);
                }
            } else {
                nleaves += 1;
                nprims += node->count;
                min_depth = ym::min(min_depth, node_depth[1]);
                max_depth = ym::max(max_depth, node_depth[1]);
            }
        }
    }
}

//
// Compute BVH stats.
//
YBVH_API void compute_bvh_stats(const scene* scn, bool include_shapes,
    int& nprims, int& ninternals, int& nleaves, int& min_depth, int& max_depth,
    int req_shape) {
    // init out variables
    nprims = 0;
    ninternals = 0;
    nleaves = 0;
    min_depth = 10000;
    max_depth = -1;

    _compute_bvh_stats(scn, req_shape, include_shapes, 0, nprims, ninternals,
        nleaves, min_depth, max_depth);
}

}  // namespace
