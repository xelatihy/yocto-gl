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

#include "yocto_geometry.h"
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

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct bvh_node {
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
struct bvh_tree {
  vector<bvh_node> nodes      = {};
  vector<int>      primitives = {};
};

// BVH span to give a view over an array
template <typename T>
struct bvh_span {
  explicit bvh_span(const vector<T>& data)
      : _data{data.data()}, _size{data.size()} {}
  bvh_span()                = default;
  bvh_span(const bvh_span&) = default;
  bvh_span& operator=(const bvh_span&) = default;
  bool      empty() const { return _size == 0; }
  size_t    size() const { return _size; }
  const T&  operator[](int idx) const { return _data[idx]; }
  const T*  begin() const { return _data; }
  const T*  end() const { return _data + _size; }
  const T*  data() const { return _data; }
  const T*  _data = nullptr;
  size_t    _size = 0;
};

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_shape {
  // elements
  bvh_span<int>   points    = {};
  bvh_span<vec2i> lines     = {};
  bvh_span<vec3i> triangles = {};
  bvh_span<vec4i> quads     = {};

  // vertices
  bvh_span<vec3f> positions = {};
  bvh_span<float> radius    = {};

  // owned elements
  vector<int>   points_data    = {};
  vector<vec2i> lines_data     = {};
  vector<vec3i> triangles_data = {};
  vector<vec4i> quads_data     = {};

  // owned vertices
  vector<vec3f> positions_data = {};
  vector<float> radius_data    = {};

  // nodes
  bvh_tree bvh = {};
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

// Callback to get instance properties
using bvh_instance_callback = function<bvh_instance(int)>;

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_scene {
  // instances and shapes
  int                   num_instances  = 0;
  bvh_instance_callback instance_cb    = {};
  vector<bvh_instance>  instances_data = {};
  vector<bvh_shape*>    shapes         = {};

  // nodes
  bvh_tree bvh = {};
#ifdef YOCTO_EMBREE
  RTCScene embree_bvh = nullptr;
#endif
  ~bvh_scene();
};

// Set shapes
int  add_shape(bvh_scene* bvh, const vector<int>& points,
     const vector<vec2i>& lines, const vector<vec3i>& triangles,
     const vector<vec4i>& quads, const vector<vec3f>& positions,
     const vector<float>& radius, bool as_view = false);
void set_shape(bvh_scene* bvh, int shape_id, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<float>& radius, bool as_view = false);

// Set instances
void set_instances(bvh_scene* bvh, int num_instances,
    bvh_instance_callback instance_cb, bool as_view = false);

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
void init_bvh(bvh_scene* bvh, const bvh_params& params,
    const progress_callback& progress_cb = {});

// Refit bvh data
void update_bvh(bvh_scene* bvh, const vector<int>& updated_instances,
    const vector<int>&       updated_shapes,
    const progress_callback& progress_cb = {});

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

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bvh_intersection intersect_bvh(const bvh_scene* bvh, const ray3f& ray,
    bool find_any = false, bool non_rigid_frames = true);
bvh_intersection intersect_bvh(const bvh_scene* bvh, int instance,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bvh_intersection overlap_shape_bvh(const bvh_shape* bvh, const vec3f& pos,
    float max_distance, bool find_any = false);
bvh_intersection overlap_bvh(const bvh_scene* bvh, const vec3f& pos,
    float max_distance, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto

#endif
