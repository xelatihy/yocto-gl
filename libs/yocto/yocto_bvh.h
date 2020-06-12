//
// # Yocto/Bvh: Tiny Library for ray-intersection and point-overlap queries
//
//
// Yocto/Bvh implements ray-intersection and point-overlap queries using a
// two-level Bvh data strcture.
//
//
// ## Ray-Scene and Closest-Point Queries
//
// Yocto/BVH is a simple implementation of ray intersection and
// closest queries using a two-level BVH data structure. We also include
// low-level intersection and closet point primitives.
// Alternatively the library also support wrapping Intel's Embree.
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

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <memory>
#include <tuple>
#include <unordered_map>

#include "yocto_math.h"

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH, RAY INTERSECTION AND OVERLAP QUERIES
// -----------------------------------------------------------------------------
namespace yocto {

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct bvh_node_ {
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
struct bvh_tree_ {
  vector<bvh_node_> nodes      = {};
  vector<int>       primitives = {};
};

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_shape {
  // elements
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // vertices
  vector<vec3f> positions = {};
  vector<float> radius    = {};

  // nodes
  bvh_tree_ bvh = {};
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif
  ~bvh_shape();
};

// instance
struct bvh_instance {
  frame3f frame = identity3x4f;
  int     shape = -1;
};

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_scene {
  // instances and shapes
  vector<bvh_instance*> instances = {};
  vector<bvh_shape*>    shapes    = {};

  // nodes
  bvh_tree_ bvh = {};
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif
  ~bvh_scene();
};

// Build the bvh acceleration structure.
void init_bvh(bvh_shape* bvh, bool embree = false);
void init_bvh(bvh_scene* bvh, bool embree = false);

// Refit bvh data
void update_bvh(bvh_shape* bvh);
void update_bvh(bvh_scene* bvh, const vector<int>& updated_instances,
    const vector<int>& updated_shapes);

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_shape_intersection {
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_scene_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bvh_shape_intersection intersect_shape_bvh(
    const bvh_shape* bvh, const ray3f& ray, bool find_any = false);
bvh_scene_intersection intersect_scene_bvh(const bvh_scene* bvh,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);
bvh_shape_intersection intersect_instance_bvh(const bvh_scene* bvh,
    int instance, const ray3f& ray, bool find_any = false,
    bool non_rigid_frames = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bvh_shape_intersection overlap_shape_bvh(const bvh_shape* bvh, const vec3f& pos,
    float max_distance, bool find_any = false);
bvh_scene_intersection overlap_scene_bvh(const bvh_scene* bvh, const vec3f& pos,
    float max_distance, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto

#endif
