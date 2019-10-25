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
// Shape and scene data is not copied from the application to the BVH to
// improve memory footprint at the price of convenience. Shape data is
// explixitly passed on evey call, while instance data uses callbacks,
// since each application has its own conventions for storing those.
// To make usage more convenite, we provide `bvh_XXX_data` to hold application
// data and convenience wrappers for all functions.
//
// We support working either on the whole scene or on a single shape. In the
// description below yoi will see this dual API defined.
//
// 1. build the shape/scene BVH with `make_XXX_bvh()`;
// 2. perform ray-shape intersection with `intersect_XXX_bvh()`
// 3. perform point overlap queries with `overlap_XXX_bvh()`
// 4. refit BVH for dynamic applications with `update_XXX_bvh`
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

#include "yocto_common.h"
#include "yocto_math.h"

#if YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct bvh_node {
  bbox3f bbox;
  int    start;
  short  num;
  bool   internal;
  byte   axis;
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// Application data is not stored explicitly.
struct bvh_tree {
  vector<bvh_node> nodes      = {};
  vector<int>      primitives = {};
};

// Make shape bvh
void make_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel);
void make_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel);
void make_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel);
void make_quads_bvh(bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel);

// Updates shape bvh for changes in positions and radia
void update_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
void update_quads_bvh(
    bvh_tree& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions);

// Find a shape element or scene instances that intersects a ray,
// returning either the closest or any overlap depending on `find_any`.
// Returns the point distance, the instance id, the shape element index and
// the element barycentric coordinates.
bool intersect_points_bvh(const bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, int& element, vec2f& uv, float& distance,
    bool find_any = false);
bool intersect_lines_bvh(const bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, int& element, vec2f& uv, float& distance,
    bool find_any = false);
bool intersect_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, int& element, vec2f& uv, float& distance,
    bool find_any = false);
bool intersect_quads_bvh(const bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const ray3f& ray, int& element, vec2f& uv,
    float& distance, bool find_any = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bool overlap_points_bvh(const bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any = false);
bool overlap_lines_bvh(const bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any = false);
bool overlap_triangles_bvh(const bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any = false);
bool overlap_quads_bvh(const bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION USING INTEL's EMBREE
// -----------------------------------------------------------------------------
namespace yocto {

#if YOCTO_EMBREE
// Wrapper to Interl's Embree
struct bvh_embree {
  bvh_embree() {}
  bvh_embree(const bvh_embree&);
  ~bvh_embree();

  bvh_embree& operator=(const bvh_embree&);

  RTCDevice           device    = nullptr;
  RTCScene            scene     = nullptr;
  RTCGeometry         shape     = nullptr;
  vector<RTCGeometry> instances = {};
};

// Make shape bvh with Intel's Embree
void make_lines_embree_bvh(bvh_embree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool compact);
void make_triangles_embree_bvh(bvh_embree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, bool high_quality, bool compact);
void make_quads_embree_bvh(bvh_embree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, bool high_quality, bool compact);

// Updates shape bvh for changes in positions and radia with Intel's Embree
void update_lines_embree_bvh(bvh_embree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_triangles_embree_bvh(bvh_embree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
void update_quads_embree_bvh(bvh_embree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions);

// Intersect a ray with either a shape or a scene
bool intersect_elements_embree_bvh(const bvh_embree& bvh, const ray3f& ray,
    int& element, vec2f& uv, float& distance, bool find_any);
bool intersect_instances_embree_bvh(const bvh_embree& bvh, const ray3f& ray,
    int& instance, int& element, vec2f& uv, float& distance, bool find_any);
#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// HIGH-LEVEL BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_shape {
  // elements
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};
  vector<vec4i> quadspos  = {};

  // vertices
  vector<vec3f> positions = {};
  vector<float> radius    = {};

  // nodes
  bvh_tree bvh = {};

#if YOCTO_EMBREE
  bvh_embree embree = {};
#endif
};

// instance
struct bvh_instance {
frame3f frame = identity3x4f;
int     shape = -1;
};

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_scene {
  // instances and shapes
  vector<bvh_instance> instances = {};
  vector<bvh_shape>    shapes    = {};

  // nodes
  bvh_tree bvh = {};

#if YOCTO_EMBREE
  bvh_embree embree = {};
#endif
};

// bvh build params
struct bvh_params {
  bool high_quality = false;
#if YOCTO_EMBREE
  bool embree  = false;
  bool compact = false;
#endif
  bool noparallel = false;
};

// Build the bvh acceleration structure.
void init_shape_bvh(bvh_shape& bvh, const bvh_params& params);
void init_scene_bvh(bvh_scene& bvh, const bvh_params& params);

// Refit bvh data
void update_shape_bvh(bvh_shape& bvh, const bvh_params& params);
void update_scene_bvh(bvh_scene& bvh, const vector<int>& updated_instances,
    const vector<int>& updated_shapes, const bvh_params& params);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bool intersect_shape_bvh(const bvh_shape& bvh, const ray3f& ray, int& element,
    vec2f& uv, float& distance, bool find_any = false);
bool intersect_scene_bvh(const bvh_scene& bvh, const ray3f& ray, int& instance,
    int& element, vec2f& uv, float& distance, bool find_any = false,
    bool non_rigid_frames = true);
// Intersects a single instance.
bool intersect_instance_bvh(const bvh_scene& bvh, int instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance,
    bool find_any = false, bool non_rigid_frames = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bool overlap_shape_bvh(const bvh_shape& bvh, const vec3f& pos,
    float max_distance, int& element, vec2f& uv, float& distance,
    bool find_any = false);
bool overlap_scene_bvh(const bvh_scene& bvh, const vec3f& pos,
    float max_distance, int& instance, int& element, vec2f& uv, float& distance,
    bool find_any = false, bool non_rigid_frames = true);

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

bvh_intersection intersect_shape_bvh(
    const bvh_shape& bvh, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_scene_bvh(const bvh_scene& bvh, const ray3f& ray,
    bool find_any = false, bool non_rigid_frames = true);
bvh_intersection intersect_instance_bvh(const bvh_scene& bvh, int instance,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

bvh_intersection overlap_shape_bvh(const bvh_shape& bvh, const vec3f& pos,
    float max_distance, bool find_any = false);
bvh_intersection overlap_scene_bvh(const bvh_scene& bvh, const vec3f& pos,
    float max_distance, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print bvh statistics.
vector<string> format_stats(const bvh_shape& bvh);
vector<string> format_stats(const bvh_scene& bvh);

}  // namespace yocto

#endif
