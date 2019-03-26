//
// # Yocto/BVH: Tiny library for ray-object intersection using a BVH
//
//
// Yocto/BVH is a simple implementation of ray intersection and
// closest queries using a two-level BVH data structure. We also include
// low-level intersection and closet point primitives.
// Alternatively the library also support wrapping Intel's Embree.
//
//
// ## Ray-Scene and Closest-Point Queries
//
// Yocto/GL provides ray-scene intersection for points, lines, triangles and
// quads accelerated by a two-level BVH data structure. Our BVH is written for
// minimal code and not maximum speed, but still gives reasonable results. We
// suggest the use of Intel's Embree as a more efficient alternative.
//
// In Yocto/Bvh, shapes are described as collections of indexed primitives
// (points/lines/triangles/quads) like the standard triangle mesh used in
// real-time graphics. A scene if represented as transformed instances of
// shapes. The internal data structure is a two-level BVH, with a BVH for each
// shape and one top-level BVH for the whole scene. This design support
// instancing for large scenes and easy BVH refitting for interactive
// applications.
//
// In these functions triangles are parameterized with uv written
// w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
// as pairs of two triangles p0,p1,p3 and p2,p3,p1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
// represent triangles correctly, an this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
// Shape and scene data is copied from the application to the BVH by default.
// This is safer and avoids memory corruption problems. We support
// data sharing with the application to reduce memory usage. Since this is
// less safe, we require this option to be explicitly defined during
// initialization. In the case of using Embree and sharing, positions arrays
// need to be slightly overallocated as per Embree documentation.
//
// We support working either on the whole scene or on a single shape. In the
// description below yoi will see this dual API defined.
//
// 1. create shape/scene bvhs using `init_XXX_bvh()`
// 3. build the shape/scene BVH with `build_bvh()`;
//    for the scene build, we will build appropriately the shape bvhs
// 4. perform ray-shape intersection with `intersect_bvh()`
// 5. perform point overlap queries with `overlap_bvh()`
// 6. a shape bvh can be updated if change to its positions occurs;
//    to do so, update the bvh data with `update_vertices()` and then call
//    `refit_bvh()`
// 7. an instance bvh can be updated if changes to instances or shapes occur;
//    to do so update the relevant shape bvhs as in 6. and update the bvh data
//    with `update_vertices()` and then call `refit_bvh()`; you can retrive a
//    shape bvh by index with `get_shape_bvh()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_BVH_H_
#define _YOCTO_BVH_H_

#ifndef YOCTO_EMBREE
#define YOCTO_EMBREE 1
#endif

#ifndef YOCTO_QUADS_AS_TRIANGLES
#define YOCTO_QUADS_AS_TRIANGLES 1
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_utils.h"

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh_scene for more details.
struct bvh_node {
    bbox3f bbox;
    short  num_primitives;
    bool   is_internal;
    byte   split_axis;
    int    primitive_ids[bvh_max_prims];
};

// internal class for shared memory references
template <typename T>
struct bvh_array_view {
    constexpr bvh_array_view() : ptr{nullptr}, count{0} {}
    constexpr bvh_array_view(const T* ptr, int count)
        : ptr{ptr}, count{count} {}
    constexpr bvh_array_view(const vector<T>& v)
        : ptr{v.data()}, count{(int)v.size()} {}

    constexpr bool     empty() const { return count == 0; }
    constexpr int      size() const { return count; }
    constexpr const T& operator[](int idx) const { return ptr[idx]; }
    constexpr const T* data() const { return ptr; }
    constexpr const T* begin() const { return ptr; }
    constexpr const T* end() const { return ptr + count; }

   private:
    const T* ptr   = nullptr;
    int      count = 0;
};

// internal class for shared memory references
template <typename T>
struct bvh_strided_view {
    constexpr bvh_strided_view() : ptr{nullptr}, count{0} {}
    constexpr bvh_strided_view(const void* ptr, int count, int stride)
        : ptr{ptr}, count{count}, stride{stride} {}
    constexpr bvh_strided_view(const vector<T>& v)
        : ptr{v.data()}, count{(int)v.size()}, stride{sizeof(T)} {}

    constexpr bool     empty() const { return count == 0; }
    constexpr int      size() const { return count; }
    constexpr const T& operator[](int idx) const {
        return *(T*)((char*)ptr + (size_t)stride * (size_t)idx);
    }

    struct iterator {
        const void* ptr    = nullptr;
        int         stride = 0;
        bool        operator!=(const iterator& other) const {
            return ptr != other.ptr;
        }
        iterator& operator++() { ptr = (char*)ptr + stride; }
        const T&  operator*() const { return *(T*)(char*)ptr; }
    };
    iterator begin() const { return {ptr, stride}; }
    iterator end() const { return {(char*)ptr + stride * count, stride}; }

   private:
    const void* ptr    = nullptr;
    int         count  = 0;
    int         stride = 0;
};

// Instance for a scene BVH.
struct bvh_instance {
    frame3f frame    = identity_frame3f;
    int     shape_id = -1;
};

// BVH for shapes made of points, lines, triangles or quads. Only one primitive
// type can be used.
// To build, fill in the shape data, then call `init_XXX_bvh()`.
// The BVH is stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
struct bvh_tree {
    // shape data views
    bvh_array_view<vec3f> positions_view;
    bvh_array_view<float> radius_view;
    bvh_array_view<int>   points_view;
    bvh_array_view<vec2i> lines_view;
    bvh_array_view<vec3i> triangles_view;
    bvh_array_view<vec4i> quads_view;

    // shape data copies
    vector<vec3f> positions_data;
    vector<float> radius_data;
    vector<int>   points_data;
    vector<vec2i> lines_data;
    vector<vec3i> triangles_data;
    vector<vec4i> quads_data;

    // instance data views
    bvh_strided_view<bvh_instance> instances_view;

    // instance data copies
    vector<bvh_instance> instances_data;

    // other data for instance bvh
    vector<bvh_tree> shape_bvhs;
    bool             non_rigid_frames = true;

    // bvh internal nodes
    vector<bvh_node> nodes;

#if YOCTO_EMBREE
    // Embree opaque data
    void* embree_bvh       = nullptr;
    bool  embree_flattened = false;
    // Cleanup for embree data
    ~bvh_tree();
#endif
};

// Options for build bvh
struct build_bvh_options {
    bool          share_memory   = false;
    bool          high_quality   = false;
    bool          use_embree     = false;
    bool          flatten_embree = false;
    bool          run_serially   = false;
    atomic<bool>* cancel_flag    = nullptr;
};

// Build a BVH from the given set of shape primitives. Data is copied or shared
// depending on `copy_data`.
void init_points_bvh(bvh_tree& bvh, bvh_array_view<int> points,
    bvh_array_view<vec3f> positions, bvh_array_view<float> radius,
    const build_bvh_options& options = {});
void init_lines_bvh(bvh_tree& bvh, bvh_array_view<vec2i> lines,
    bvh_array_view<vec3f> positions, bvh_array_view<float> radius,
    const build_bvh_options& options = {});
void init_triangles_bvh(bvh_tree& bvh, bvh_array_view<vec3i> triangles,
    bvh_array_view<vec3f> positions, const build_bvh_options& options = {});
void init_quads_bvh(bvh_tree& bvh, bvh_array_view<vec4i> quads,
    bvh_array_view<vec3f> positions, const build_bvh_options& options = {});
void init_instances_bvh(bvh_tree& bvh, bvh_strided_view<bvh_instance> instances,
    int num_shapes, const build_bvh_options& options = {});
bvh_tree& get_shape_bvh(bvh_tree& bvh, int shape_id);

// Build the bvh acceleration structure.
void build_bvh(bvh_tree& bvh, const build_bvh_options& options = {});

// Update the node data for shape and scene bvhs.
void update_vertices(bvh_tree& bvh, bvh_array_view<vec3f> positions,
    bvh_array_view<float> radius = {});
void update_instances(bvh_tree& bvh, bvh_strided_view<bvh_instance> instances);

// Refit bvh data
void refit_bvh(bvh_tree& bvh);

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_intersection {
    int   instance_id = -1;
    int   element_id  = -1;
    vec2f element_uv  = {0, 0};
    float distance    = 0;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bool intersect_bvh(const bvh_tree& bvh, const ray3f& ray,
    bvh_intersection& intersection, bool find_any = false);
bool intersect_bvh(const bvh_tree& bvh, int instance_id, const ray3f& ray,
    bvh_intersection& intersection, bool find_any = false);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bool overlap_bvh(const bvh_tree& bvh, const vec3f& pos, float max_distance,
    bvh_intersection& intersection, bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, vec2f& uv, float& dist);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, vec2f& uv, float& dist);

// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, vec2f& uv, float& dist);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, vec2f& uv, float& dist);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "A Ray-Box Intersection Algorithm and Efficient Dynamic Voxel Rendering" at
// http://jcgt.org/published/0007/03/04/
// but using the Wald implementation
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p0, float r0,
    vec2f& uv, float& dist);

// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& p0, const vec3f& p1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, float r0, float r1, vec2f& uv, float& dist);

// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& p0, const vec3f& p1, const vec3f& p2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, float r0, float r1, float r2, vec2f& uv,
    float& dist);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
    float r2, float r3, vec2f& uv, float& dist);

// Check if a bounding box overlaps a position within a max distance.
bool overlap_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox);

// Check if two bounding boxes overlap.
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print bvh statistics.
string print_bvh_stats(const bvh_tree& bvh);

}  // namespace yocto

#endif
