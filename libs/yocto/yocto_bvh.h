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
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef YOCTO_BVH_H_
#define YOCTO_BVH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <cstdint>
#include <memory>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_math.h"
#include "yocto_scene.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
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
  bbox3f  bbox     = {};
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

// Shape BVHs are just the bvh for the shape.
struct shape_bvh {
  bvh_tree bvh = {};
};

// Scene BVHs store the bvh for instances and shapes.
// Application data is not stored explicitly.
struct scene_bvh {
  bvh_tree          bvh    = {};
  vector<shape_bvh> shapes = {};
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
// shape element id, shape element uv and intersection distance.
struct shape_intersection {
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
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
// BVH FOR COLLECTION OF PRIMITIVES
// -----------------------------------------------------------------------------
namespace yocto {

// Build elements bvh
bvh_tree make_points_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool highquality = false);
bvh_tree make_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool highquality = false);
bvh_tree make_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, bool highquality = false);
bvh_tree make_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, bool highquality = false);

// Refit elements bvh
void update_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
void update_quads_bvh(
    bvh_tree& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions);

// Bvh intersection
shape_intersection intersect_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
shape_intersection intersect_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
shape_intersection intersect_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);
shape_intersection intersect_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EMBREE WRAPPER
// -----------------------------------------------------------------------------
namespace yocto {

// Wrapper for Intel's Embree
using ebvh_tree = unique_ptr<void, void (*)(void*)>;

// Wrapper for Intel Embree.
struct shape_ebvh {
  ebvh_tree ebvh = {nullptr, nullptr};
};

// Wrapper for Intel Embree.
struct scene_ebvh {
  ebvh_tree          ebvh   = {nullptr, nullptr};  // instances
  vector<shape_ebvh> shapes = {};                  // shapes
};

// Check if embree is supported
bool embree_supported();

// Build the bvh acceleration structure.
shape_ebvh make_shape_ebvh(const shape_data& shape, bool highquality = false);
scene_ebvh make_scene_ebvh(
    const scene_data& scene, bool highquality = false, bool noparallel = false);

// Refit bvh data
void update_shape_ebvh(shape_ebvh& bvh, const shape_data& shape);
void update_scene_ebvh(scene_ebvh& bvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
shape_intersection intersect_shape_ebvh(const shape_ebvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false);
scene_intersection intersect_scene_ebvh(const scene_ebvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any = false);
scene_intersection intersect_instance_ebvh(const scene_ebvh& bvh,
    const scene_data& scene, int instance, const ray3f& ray,
    bool find_any = false);

}  // namespace yocto

#endif
