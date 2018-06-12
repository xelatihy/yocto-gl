//
// # Yocto/Bvh: Tiny C++ Library for Ray Intersection and Closest-Point Queries
//
// Yocto/Bvh provides ray-scene intersection for points, lines and triangles
// accelerated by a two-level BVH data structure. Quads are supported in an
// experimental fashion.
// Our BVH is written for minimal code and not maximum speed, but still gives
// reasonable results. We suggest the use of Intel's Embree as a more efficient
// alternative.
//
// In Yocto/Bvh, shapes are described as collections of indexed primitives
// (points/lines/triangles) like the standard triangle mesh used in real-time
// graphics. A scene if represented as transformed instances of shapes. The
// internal data structure is a two-level BVH, with a BVH for each shape and
// one top-level BVH for the whole scene. This design support ionstancing
// for large scenes and easy BVH refitting for interactive applictions.
//
//
// ## Usage
//
// 1. fill the shape or instance data
// 2. build the BVH with `build_bvh()`
// 3. perform ray-element intersection with `intersect_bvh()`
// 4. perform point overlap queries with `overlap_bvh()`
// 5. refit the BVH with `refit_bvh()` after updating internal data
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
#ifndef _YGL_BVH_H_
#define _YGL_BVH_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <vector>
#include <memory>

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(
    const ray3f& ray, vec3f p, float r, float& dist, vec2f& uv);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, vec3f v0, vec3f v1, float r0, float r1,
    float& dist, vec2f& uv);

// Intersect a ray with a triangle.
bool intersect_triangle(
    const ray3f& ray, vec3f v0, vec3f v1, vec3f v2, float& dist, vec2f& uv);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, vec3f v0, vec3f v1, vec3f v2, vec3f v3,
    float& dist, vec2f& uv);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(
    const ray3f& ray, vec3f ray_dinv, vec3i ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(
    vec3f pos, float dist_max, vec3f v0, float r0, float& dist, vec2f& uv);

// Find closest line point to a position.
float closestuv_line(vec3f pos, vec3f v0, vec3f v1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(vec3f pos, float dist_max, vec3f v0, vec3f v1, float r0,
    float r1, float& dist, vec2f& uv);

// Find closest triangle point to a position.
vec2f closestuv_triangle(vec3f pos, vec3f v0, vec3f v1, vec3f v2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(vec3f pos, float dist_max, vec3f v0, vec3f v1, vec3f v2,
    float r0, float r1, float r2, float& dist, vec2f& uv);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(vec3f pos, float dist_max, vec3f v0, vec3f v1, vec3f v2,
    vec3f v3, float r0, float r1, float r2, float r3, float& dist, vec2f& uv);

// Check if a bouning box overlaps a position within a max distance.
bool overlap_bbox(vec3f pos, float dist_max, const bbox3f& bbox);

// Check if two bouning boxes overlap.
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

}  // namespace ygl

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace ygl {

// forward declaration
struct bvh_tree;

// Type of BVH node.
enum struct bvh_node_type : uint8_t {
    internal,
    point,
    line,
    triangle,
    vertex,
    instance
};

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh_tree for more details.
struct bvh_node {
    bbox3f bbox;                    // bouds
    uint32_t prims[bvh_max_prims];  // primitives
    uint16_t count;                 // number of prims
    bvh_node_type type;             // node type
    uint8_t split_axis;             // split axis
};

// BVH tree, stored as a node array. The tree structure is encoded using array
// indices instead of pointers, both for speed but also to simplify code.
// BVH nodes indices refer to either the node array, for internal nodes,
// or the primitive arrays, for leaf nodes. BVH trees may contain only one type
// of geometric primitive, like points, lines, triangle or instances of other
// BVHs. To handle multiple primitive types and transformed primitives, build
// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
// BVHs, shape BVHs, each of which of a uniform primitive type.
// To build a BVH, first fill in either the shape or instance data, then
// call `build_bvh()`.
struct bvh_tree {
    // data for shape BVH
    std::vector<vec3f> pos;        // Positions for shape BVHs.
    std::vector<float> radius;     // Radius for shape BVHs.
    std::vector<int> points;       // Points for shape BVHs.
    std::vector<vec2i> lines;      // Lines for shape BVHs.
    std::vector<vec3i> triangles;  // Triangles for shape BVHs.

    // data for instance BVH
    std::vector<frame3f> ist_frames;                  // instance frames
    std::vector<frame3f> ist_inv_frames;              // instance inverse frames
    std::vector<std::shared_ptr<bvh_tree>> ist_bvhs;  // instance shape bvhs

    // bvh nodes
    std::vector<bvh_node> nodes;  // Internal nodes.
};

// Build a BVH from the given set of primitives.
void build_bvh(std::shared_ptr<bvh_tree> bvh, bool equalsize);
// Update the node bounds for a shape bvh.
void refit_bvh(std::shared_ptr<bvh_tree> bvh);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance `dist`, the instance
// id `iid`, the shape id `sid`, the shape element index `eid` and the
// shape barycentric coordinates `uv`.
bool intersect_bvh(const std::shared_ptr<bvh_tree> bvh, const ray3f& ray,
    bool find_any, float& dist, int& iid, int& eid, vec2f& uv);

// Find a shape element that overlaps a point within a given distance
// `max_dist`, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance `dist`, the instance id `iid`, the
// shape id `sid`, the shape element index `eid` and the shape barycentric
// coordinates `uv`.
bool overlap_bvh(const std::shared_ptr<bvh_tree> bvh, vec3f pos, float max_dist,
    bool find_any, float& dist, int& iid, int& eid, vec2f& uv);

}  // namespace ygl

#endif
