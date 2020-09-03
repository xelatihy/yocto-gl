//
// # Yocto/Bvh: Accelerated ray-intersection and point-overlap
//
// Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
// by a two-level BVH using an internal or wrapping Embree.
// Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.
// For now, use the Yocto/Scene support for this.

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

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

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
using std::function;
using std::string;
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
  bbox3f  bbox     = invalidb3f;
  int32_t start    = 0;
  int16_t num      = 0;
  int8_t  axis     = 0;
  bool    internal = false;
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
  frame3f    frame = identity3x4f;
  bvh_shape* shape = nullptr;
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

// Create BVH
bvh_shape*    add_shape(bvh_scene* scene);
bvh_instance* add_instance(bvh_scene* scene);

// set shape properties
void set_points(bvh_shape* shape, const vector<int>& points);
void set_lines(bvh_shape* shape, const vector<vec2i>& lines);
void set_triangles(bvh_shape* shape, const vector<vec3i>& triangles);
void set_quads(bvh_shape* shape, const vector<vec4i>& quads);
void set_positions(bvh_shape* shape, const vector<vec3f>& positions);
void set_radius(bvh_shape* shape, const vector<float>& radius);

// set instance properties
void set_frame(bvh_instance* instance, const frame3f& frame);
void set_shape(bvh_instance* instance, bvh_shape* shape);

// Create BVH shortcuts
bvh_shape*    add_shape(bvh_scene* bvh, const vector<int>& points,
       const vector<vec2i>& lines, const vector<vec3i>& triangles,
       const vector<vec4i>& quads, const vector<vec3f>& positions,
       const vector<float>& radius);
bvh_instance* add_instance(
    bvh_scene* bvh, const frame3f& frame, bvh_shape* shape);

// Strategy used to build the bvh
enum struct bvh_build_type {
  default_,
  highquality,
  middle,
  balanced,
#ifdef YOCTO_EMBREE
  embree_default,
  embree_highquality,
  embree_compact  // only for copy interface
#endif
};

const auto bvh_build_names = vector<string>{
    "default", "highquality", "middle", "balanced",
#ifdef YOCTO_EMBREE
    "embree-default", "embree-highquality", "embree-compact"
#endif
};

// Bvh parameters
struct bvh_params {
  bvh_build_type bvh        = bvh_build_type::default_;
  bool           noparallel = false;  // only serial momentarily
};

// Progress report callback
using progress_callback =
    function<void(const string& message, int current, int total)>;

// Build the bvh acceleration structure.
void init_bvh(bvh_shape* bvh, const bvh_params& params);
void init_bvh(bvh_scene* bvh, const bvh_params& params,
    const progress_callback& progress_cb = {});

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
