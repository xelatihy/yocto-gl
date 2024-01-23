//
// # Yocto/EBvh: Wrapper for Embree
//
// Yocto/EBvh provides ray-intersection accelerated Embree
// Yocto/EBvh is implemented in `yocto_ebvh.h` and `yocto_ebvh.cpp`.
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

#ifndef _YOCTO_EBVH_H_
#define _YOCTO_EBVH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"

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

// -----------------------------------------------------------------------------
// EMBREE BACKWARD COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

// backward compatibility
using shape_embree_bvh = shape_ebvh;
using scene_embree_bvh = scene_ebvh;

// Build the bvh acceleration structure.
[[deprecated]] inline shape_embree_bvh make_shape_embree_bvh(
    const shape_data& shape, bool highquality = false) {
  return make_shape_ebvh(shape, highquality);
}
[[deprecated]] inline scene_embree_bvh make_scene_embree_bvh(
    const scene_data& scene, bool highquality = false,
    bool noparallel = false) {
  return make_scene_ebvh(scene, highquality, noparallel);
}

// Refit bvh data
[[deprecated]] inline void update_shape_embree_bvh(
    shape_embree_bvh& bvh, const shape_data& shape) {
  return update_shape_ebvh(bvh, shape);
}
[[deprecated]] inline void update_scene_embree_bvh(scene_ebvh& bvh,
    const scene_data& scene, const vector<int>& updated_instances,
    const vector<int>& updated_shapes) {
  return update_scene_ebvh(bvh, scene, updated_instances, updated_shapes);
}

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
[[deprecated]] inline shape_intersection intersect_shape_embree_bvh(
    const shape_ebvh& bvh, const shape_data& shape, const ray3f& ray,
    bool find_any = false) {
  return intersect_shape_ebvh(bvh, shape, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_scene_embree_bvh(
    const scene_ebvh& bvh, const scene_data& scene, const ray3f& ray,
    bool find_any = false) {
  return intersect_scene_ebvh(bvh, scene, ray, find_any);
}
[[deprecated]] inline scene_intersection intersect_instance_embree_bvh(
    const scene_ebvh& bvh, const scene_data& scene, int instance,
    const ray3f& ray, bool find_any = false) {
  return intersect_instance_ebvh(bvh, scene, instance, ray, find_any);
}

}  // namespace yocto

#endif
