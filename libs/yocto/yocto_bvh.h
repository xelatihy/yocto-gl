//
// # Yocto/Bvh: Accelerated ray-intersection and point-overlap
//
// Yocto/Bvh provides ray-intersection and point-overlap queries accelerated
// by a two-level BVH.
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
// BVH FOR SHAPE ELEMENTS
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

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct shape_intersection {
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Make shape bvh
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
    const ray3f& ray, bool find_any = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
shape_intersection overlap_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, vec3f pos, float max_distance,
    bool find_any = false);
shape_intersection overlap_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, vec3f pos, float max_distance,
    bool find_any = false);
shape_intersection overlap_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<float>& radius, vec3f pos, float max_distance,
    bool find_any = false);
shape_intersection overlap_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<float>& radius, vec3f pos, float max_distance,
    bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE AND SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

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
    const shape_data& shape, vec3f pos, float max_distance,
    bool find_any = false);
scene_intersection overlap_scene_bvh(const scene_bvh& bvh,
    const scene_data& scene, vec3f pos, float max_distance,
    bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONVENIENCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience functions
inline vec3f eval_position(
    const scene_data& scene, const scene_intersection& intersection);
inline vec3f eval_normal(
    const scene_data& scene, const scene_intersection& intersection);
inline vec3f eval_element_normal(
    const scene_data& scene, const scene_intersection& intersection);
inline vec3f eval_shading_position(const scene_data& scene,
    const scene_intersection& intersection, vec3f outgoing);
inline vec3f eval_shading_normal(const scene_data& scene,
    const scene_intersection& intersection, vec3f outgoing);
inline vec2f eval_texcoord(
    const scene_data& scene, const scene_intersection& intersection);
inline material_point eval_material(
    const scene_data& scene, const scene_intersection& intersection);
inline bool is_volumetric(
    const scene_data& scene, const scene_intersection& intersection);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// CONVENIENCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience functions
inline vec3f eval_position(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
inline vec3f eval_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
inline vec3f eval_element_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_element_normal(
      scene, scene.instances[intersection.instance], intersection.element);
}
inline vec3f eval_shading_position(const scene_data& scene,
    const scene_intersection& intersection, vec3f outgoing) {
  return eval_shading_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
inline vec3f eval_shading_normal(const scene_data& scene,
    const scene_intersection& intersection, vec3f outgoing) {
  return eval_shading_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
inline vec2f eval_texcoord(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_texcoord(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
inline material_point eval_material(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_material(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
inline bool is_volumetric(
    const scene_data& scene, const scene_intersection& intersection) {
  return is_volumetric(scene, scene.instances[intersection.instance]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARD COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

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
    const shape_data& shape, vec3f pos, float max_distance,
    bool find_any = false) {
  return overlap_shape_bvh(bvh, shape, pos, max_distance, find_any);
}
[[deprecated]] inline scene_intersection overlap_bvh(const scene_bvh& bvh,
    const scene_data& scene, vec3f pos, float max_distance,
    bool find_any = false) {
  return overlap_scene_bvh(bvh, scene, pos, max_distance, find_any);
}

// backward compatibility
[[deprecated]] inline shape_intersection intersect_shape(const shape_bvh& bvh,
    const shape_data& shape, const ray3f& ray, bool find_any = false) {
  return intersect_shape_bvh(bvh, shape, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_scene(const scene_bvh& bvh,
    const scene_data& scene, const ray3f& ray, bool find_any = false) {
  return intersect_scene_bvh(bvh, scene, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_instance(
    const scene_bvh& bvh, const scene_data& scene, int instance,
    const ray3f& ray, bool find_any = false) {
  return intersect_instance_bvh(bvh, scene, instance, ray, find_any);
}

// backward compatibility
[[deprecated]] inline shape_intersection overlap_shape(const shape_bvh& bvh,
    const shape_data& shape, vec3f pos, float max_distance,
    bool find_any = false) {
  return overlap_shape_bvh(bvh, shape, pos, max_distance, find_any);
}
[[deprecated]] inline scene_intersection overlap_scene(const scene_bvh& bvh,
    const scene_data& scene, vec3f pos, float max_distance,
    bool find_any = false) {
  return overlap_scene_bvh(bvh, scene, pos, max_distance, find_any);
}

}  // namespace yocto

#endif
