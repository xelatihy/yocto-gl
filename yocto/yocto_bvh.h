///
/// # Yocto/Bvh
///
/// Ray-intersection and closet-point routines supporting points,
/// lines, triangles and tetrahedra accelerated by a two-level bounding volume
/// hierarchy (BVH). Tetrahedra support is still a work in progress.
///
/// The scene model is comprised of a set of shapes, each with its own BVH.
/// Shapes are collections of points, lines, triangles or tetrahedras. each
/// shape may contain only one element type. Shapes are organized into a scene
/// by creating shape instances, each of which with its own transform. For now,
/// we only support rigid transform (translation and rotation). But eventually
/// we might be able to support generic transforms. Shape data is shared with
/// application and not copied to save memory. Instance data is instead
/// kept internally to allow for more flexibility.
///
/// The library supports two main queries. Ray-scene intersection for raytracing
/// and point-scene overlap test for collision detection. Also for collision
/// detection we support shape bounds overlap for broad-phase collisions.
/// Queries and BVH build can be performed over the whole single or single
/// shapes.
///
/// This library depends in yocto_math.h
///
///
/// ## Usage
///
/// 1. create a scene with `make_scene()`
/// 2. for each shape, add shape data and transforms with `add_point_shape()`,
///    `add_line_shape()`, `add_triangle_shape()` and `add_tetra_shape()`; to
///    modify the frame call `set_shape_frame()`
/// 3. add shape instances with `add_instance()`
/// 4. build the bvh with `build_scene_bvh()`
/// 5. perform ray-interseciton tests with `intersect_ray()`
///     - use early_exit=false if you want to know the closest hit point
///     - use early_exit=false if you only need to know whether there is a hit
///     - for points and lines, a radius is required
///     - for triangle and tetrahedra, the radius is ignored
/// 6. perform point overlap tests with `overlap_point()` to if a point overlaps
///       with an element within a maximum distance
///     - use early_exit as above
///     - for all primitives, a radius is used if defined, but should
///       be very small compared to the size of the primitive since the radius
///       overlap is approximate
/// 7. perform shape overlap queries with `overlap_shape_bounds()`
/// 8. use `refit_bvh()` to recompute the bvh bounds if transforms or vertices
///    are (you should rebuild the bvh for large changes); update the instances'
///    transforms with `set_instance_frame()`; shapes use shared memory, so
///    no explicit update is necessary
///
///
/// ## History
///
/// - v 0.18: faster internal intersection
/// - v 0.17: removal of SAH build option (better use embree instead)
/// - v 0.16: use yocto_math in the interface and remove inline compilation
/// - v 0.15: move to add api
/// - v 0.14: added instances
/// - v 0.13: switch to .h/.cpp pair
/// - v 0.12: doxygen comments
/// - v 0.11: opaque API (allows for changing internals without altering API)
/// - v 0.10: internally use pointers for performance transparency
/// - v 0.9: [API change] use rigid frames rather than linear transforms
/// - v 0.8: [API change] added support for tetrahdedra and explicit shape ids
/// - v 0.7: [major API change] move to modern C++ interface
/// - v 0.6: internal refactoring
/// - v 0.5: removal of C interface
/// - v 0.4: use of STL containers
/// - v 0.3: cleaner interface based on opaque scene
/// - v 0.2: better internal intersections and cleaner bvh stats
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
///
/// ## Acknowledgements
///
/// This library includes code from "Real-Time Collision Detection"
/// by Christer Ericson and from the Journal of Graphics Techniques.
///
namespace ybvh {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

#ifndef _YBVH_H_
#define _YBVH_H_

#include <array>
#include <functional>
#include <vector>

#include "yocto_math.h"

///
/// Accellerated ray-intersection and closet-point queries.
///
namespace ybvh {

///
/// Opaque BVH scene.
///
struct scene;

///
/// Init scene.
///
scene* make_scene();

///
/// Free scene.
///
/// - Parameters:
///     - scn: scene to free
///
void free_scene(scene* scn);

///
/// Add shape
///
/// - Parameters:
///     - scn: scene
///     - npoints/point: point indices
///     - nverts: number of vertices
///     - pos/radius: vertex position and radius
/// - Returns:
///     - shape id
///
int add_point_shape(scene* scn, int npoints, const int* point, int nverts,
    const ym::vec3f* pos, const float* radius);

///
/// Add shape.
///
/// - Parameters:
///     - scn: scene
///     - nlines/line: line indices
///     - nverts: number of vertices
///     - pos/radius: vertex position and radius
/// - Returns:
///     - shape id
///
int add_line_shape(scene* scn, int nlines, const ym::vec2i* lines, int nverts,
    const ym::vec3f* pos, const float* radius);

///
/// Add shape
///
/// - Parameters:
///     - scn: scene
///     - ntriange/triangles: triangle indices
///     - nverts: number of vertices
///     - pos/radius: vertex position and radius
/// - Returns:
///     - shape id
///
int add_triangle_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos, const float* radius);

///
/// Add shape
///
/// - Parameters:
///     - scn: scene
///     - npoints/point: point indices
///     - ntetra/tetra: tetrahedra indices
///     - nverts: number of vertices
///     - pos/radius: vertex position and radius
/// - Returns:
///     - shape id
///
int add_tetra_shape(scene* scn, int ntetra, const ym::vec4i* tetra, int nverts,
    const ym::vec3f* pos, const float* radius);

///
/// Add shape
///
/// - Parameters:
///     - scn: scene
///     - nverts: number of vertices
///     - pos/radius: vertex position and radius
/// - Returns:
///     - shape id
///
int add_point_shape(
    scene* scn, int nverts, const ym::vec3f* pos, const float* radius);

///
/// Add an instance.
///
/// - Parameters:
///     - scn: scene
///     - frame: instance transform
///     - sid: shape id
/// - Returns:
///     - instance id
///
int add_instance(scene* scn, const ym::frame3f& frame, int sid);

///
/// Set an instance frame.
///
/// - Parameters:
///     - scn: scene
///     - iid: instance id
///     - frame: shape transform
///
void set_instance_frame(scene* scn, int iid, const ym::frame3f& frame);

///
/// Builds a scene BVH.
///
/// - Parameters:
///     - scn: object to build the bvh for
///     - equalsize: space-splitting tree (otherwise balanced tree)
///     - do_shapes: build shapes
///
void build_scene_bvh(scene* scn, bool equalsize = true, bool do_shapes = true);

///
/// Builds a shape BVH.
///
/// - Parameters:
///     - scn: object to build the bvh for
///     - sid: required shape
///     - equalsize: space-splitting tree (otherwise balanced tree)
///
void build_shape_bvh(scene* scn, int sid, bool equalsize = true);

///
/// Refit the bounds of each shape for moving objects. Use this only to avoid
/// a rebuild, but note that queries are likely slow if objects move a lot.
/// Before calling refit, set the scene data.
///
/// - Parameters:
///     - scn: scene to refit
///     - do_instances: refit instances
///     - do_shapes: refit shapes
///
void refit_scene_bvh(scene* scn, bool do_shapes = false);

///
/// Refit the bounds of each shape for moving objects. Use this only to avoid
/// a rebuild, but note that queries are likely slow if objects move a lot.
/// Before calling refit, set the scene data.
///
/// - Parameters:
///     - scn: scene to refit
///     - sid: shape id
///
void refit_shape_bvh(scene* scn, int sid);

///
/// BVH intersection.
///
struct intersection_point {
    /// distance
    float dist = 0;
    /// instance index
    int iid = -1;
    /// shape index
    int sid = -1;
    /// element index
    int eid = -1;
    /// element baricentric coordinates
    ym::vec4f euv = {0, 0, 0, 0};

    /// Check whether it was a hit.
    operator bool() const { return eid >= 0; }
};

///
/// Intersect the scene with a ray. Find any interstion if early_exit, otherwise
/// find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - ray: ray
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - intersection point
///
intersection_point intersect_scene(
    const scene* scn, const ym::ray3f& ray, bool early_exit);

///
/// Intersect the scene with a ray. Find any interstion if early_exit, otherwise
/// find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - sid: shape id
///     - ray: ray
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - intersection point
///
intersection_point intersect_shape(
    const scene* scn, int sid, const ym::ray3f& ray, bool early_exit);

///
/// Intersect the scene with a ray. Find any interstion if early_exit, otherwise
/// find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - sid: shape id
///     - ray: ray
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - intersection point
///
intersection_point intersect_instance(
    const scene* scn, int iid, const ym::ray3f& ray, bool early_exit);

///
/// Returns a list of instance pairs that can possibly overlap by checking only
/// they axis aligned bouds. This is only a conservative check useful for
/// collision detection.
///
/// - Parameters:
///     - scn1, scn2: scenes to overlap
///     - conservative: use conservative checks
///     - skip_self: exlude self intersections
///     - skip_duplicates: exlude intersections (i1,i2) if (i2,i1)
///       is already present
/// - Out Parameters:
///     - overlaps: vectors of shape overlaps
///
void overlap_instance_bounds(const scene* scn1, const scene* scn2,
    bool conservative, bool exclude_duplicates, bool exclude_self,
    std::vector<ym::vec2i>* overlaps);

///
/// Finds the closest element that overlaps a point within a given radius.
///
/// - Parameters:
///     - scene: scene to check
///     - pt: ray origin
///     - max_dist: max point distance
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - overlap point
///
intersection_point overlap_scene(
    const scene* scn, const ym::vec3f& pt, float max_dist, bool early_exit);

///
/// Finds the closest element that overlaps a point within a given radius.
///
/// - Parameters:
///     - scene: scene to check
///     - sid: shape id
///     - pt: ray origin
///     - max_dist: max point distance
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - overlap point
///
intersection_point overlap_shape(const scene* scn, int sid, const ym::vec3f& pt,
    float max_dist, bool early_exit);

///
/// Finds the closest element that overlaps a point within a given radius.
///
/// - Parameters:
///     - scene: scene to check
///     - iid: instance id
///     - pt: ray origin
///     - max_dist: max point distance
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - overlap point
///
intersection_point overlap_instance(const scene* scn, int iid,
    const ym::vec3f& pt, float max_dist, bool early_exit);

///
/// Compute BVH stats.
///
void compute_bvh_stats(const scene* scn, bool include_shapes, int& nprims,
    int& ninternals, int& nleaves, int& min_depth, int& max_depth,
    int req_shape = -1);

}  // namespace ybvh

#endif
