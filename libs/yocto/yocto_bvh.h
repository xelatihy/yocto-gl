//
// # Yocto/Bvh: Accelerated ray-intersection and point-overlap
//
// Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
// by a two-level BVH using an internal or wrapping Embree.
// Yocto/Bvh is implemented in `yocto_bvh.h` and `yocto_bvh.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
#include <memory>
#include <string>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_math.h"
#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::unique_ptr;
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
struct shape_bvh {
  vector<bvh_node> nodes      = {};
  vector<int>      primitives = {};
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// We also store the BVH of the contained shapes.
// Application data is not stored explicitly.
struct scene_bvh {
  vector<bvh_node>  nodes      = {};
  vector<int>       primitives = {};
  vector<shape_bvh> shapes     = {};  // shapes
};

// Build the bvh acceleration structure.
shape_bvh make_shape_bvh(const shape_data& shape, bool highquality = false);
scene_bvh make_scene_bvh(
    const scene_data& scene, bool highquality = false, bool noparallel = false);

// Refit bvh data
void update_shape_bvh(shape_bvh& bvh, const shape_data& shape);
void update_scene_bvh(scene_bvh& bvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes);

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct scene_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
shape_intersection intersect_shape_bvh(const shape_bvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false);
scene_intersection intersect_scene_bvh(const scene_bvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any = false);
scene_intersection intersect_instance_bvh(const scene_bvh& bvh,
    const scene_data& scene, int instance, const ray3f& ray,
    bool find_any = false);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
shape_intersection overlap_shape_bvh(const shape_bvh& bvh,
    const shape_data& shape, const vec3f& pos, float max_distance,
    bool find_any = false);
scene_intersection overlap_scene_bvh(const scene_bvh& bvh,
    const scene_data& scene, const vec3f& pos, float max_distance,
    bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EMBREE WRAPPER
// -----------------------------------------------------------------------------
namespace yocto {

// Wrapper for Intel Embree.
struct shape_embree_bvh {
  unique_ptr<void, void (*)(void*)> embree_bvh = {nullptr, nullptr};  // embree
};

// Wrapper for Intel Embree.
struct scene_embree_bvh {
  vector<shape_embree_bvh>          shapes     = {};                  // shapes
  unique_ptr<void, void (*)(void*)> embree_bvh = {nullptr, nullptr};  // embree
};

// Build the bvh acceleration structure.
shape_embree_bvh make_shape_embree_bvh(
    const shape_data& shape, bool highquality = false);
scene_embree_bvh make_scene_embree_bvh(
    const scene_data& scene, bool highquality = false, bool noparallel = false);

// Refit bvh data
void update_shape_embree_bvh(shape_embree_bvh& bvh, const shape_data& shape);
void update_scene_embree_bvh(scene_embree_bvh& bvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
shape_intersection intersect_shape_embree_bvh(const shape_embree_bvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false);
scene_intersection intersect_scene_embree_bvh(const scene_embree_bvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any = false);
scene_intersection intersect_instance_embree_bvh(const scene_embree_bvh& bvh,
    const scene_data& scene, int instance, const ray3f& ray,
    bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARD COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

// backward compatibility
using bvh_shape [[deprecated]] = shape_bvh;
using bvh_scene [[deprecated]] = scene_bvh;

// backward compatibility
[[deprecated]] inline shape_bvh make_bvh(
    const shape_data& shape, bool embree = false, bool highquality = false) {
  return make_shape_bvh(shape, highquality);
}
[[deprecated]] inline scene_bvh make_bvh(const scene_data& scene,
    bool embree = false, bool highquality = false, bool noparallel = false) {
  return make_scene_bvh(scene, highquality, noparallel);
}

// backward compatibility
[[deprecated]] inline void update_bvh(shape_bvh& bvh, const shape_data& shape) {
  return update_shape_bvh(bvh, shape);
}
[[deprecated]] inline void update_bvh(scene_bvh& bvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes) {
  return update_scene_bvh(bvh, scene, updated_instances, updated_shapes);
}

// backward compatibility
[[deprecated]] inline shape_intersection intersect_bvh(const shape_bvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false) {
  return intersect_shape_bvh(bvh, shape, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_bvh(const scene_bvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any = false) {
  return intersect_scene_bvh(bvh, scene, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_bvh(const scene_bvh& bvh,
    const scene_data& scene, int instance, const ray3f& ray,
    bool find_any = false) {
  return intersect_instance_bvh(bvh, scene, instance, ray, find_any);
}

// backward compatibility
[[deprecated]] inline shape_intersection overlap_bvh(const shape_bvh& bvh,
    const shape_data& shape, const vec3f& pos, float max_distance,
    bool find_any = false) {
  return overlap_shape_bvh(bvh, shape, pos, max_distance, find_any);
}
[[deprecated]] inline scene_intersection overlap_bvh(const scene_bvh& bvh,
    const scene_data& scene, const vec3f& pos, float max_distance,
    bool find_any = false) {
  return overlap_scene_bvh(bvh, scene, pos, max_distance, find_any);
}

// backward compatibility
[[deprecated]] inline shape_intersection intersect_shape(const shape_bvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false) {
  return intersect_shape_bvh(bvh, shape, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_scene(const scene_bvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any) {
  return intersect_scene_bvh(bvh, scene, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_instance(
    const scene_bvh& bvh, const scene_data& scene, int instance,
    const ray3f& ray, bool find_any) {
  return intersect_instance_bvh(bvh, scene, instance, ray, find_any);
}

// backward compatibility
[[deprecated]] inline shape_intersection overlap_shape(const shape_bvh& bvh,
    const shape_data& shape, const vec3f& pos, float max_distance,
    bool find_any = false) {
  return overlap_shape_bvh(bvh, shape, pos, max_distance, find_any);
}
[[deprecated]] inline scene_intersection overlap_scene(const scene_bvh& bvh,
    const scene_data& scene, const vec3f& pos, float max_distance,
    bool find_any = false) {
  return overlap_scene_bvh(bvh, scene, pos, max_distance, find_any);
}

}  // namespace yocto

#endif
