///
/// YOCTO_BVH: ray-intersection and closet-point routines supporting points,
/// lines, triangles and tetrahedra accelerated by a two-level bounding volume
/// hierarchy (BVH). Tetrahedra support is still a work in progress.
///
///
/// USAGE:
///
/// 1. create a scene with make_scene()
/// 2. for each shape, add shape data and transforms with set_point_shape(),
///    set_line_shape(), set_triangle_shape() and set_tetra_shape(); to
///    modify the frame call set_shape_frame()
/// 3. build the bvh with build_bvh() using the specified heuristic (or default)
/// 4. perform ray-interseciton tests with intersect_ray()
///     - use early_exit=false if you want to know the closest hit point
///     - use early_exit=false if you only need to know whether there is a hit
///     - for points and lines, a radius is required
///     - for triangle and tetrahedra, the radius is ignored
/// 5. perform point overlap tests with overlap_point() to if a point overlaps
///       with an element within a maximum distance
///     - use early_exit as above
///     - for all primitives, a radius is used if defined, but should
///       be very small compared to the size of the primitive since the radius
///       overlap is approximate
/// 6. perform shape overlap queries with overlap_shape_bounds() and
/// overlap_verts()
/// 7. use refit_bvh() to recompute the bvh bounds if transforms or vertices are
///    (you should rebuild the bvh for large changes)
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// Shapes are indexed meshes and are described by their
/// number of elements, an array of vertex indices,
/// the primitive type (points, lines, triangles),
/// an array of vertex positions, and an array of vertex radius
/// (for points and lines).
///
/// Implementation notes:
/// - using high precision bvh intersection by default
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
///
///
/// This .cpp file depends on yocto_math.h.
///
///
/// HISTORY:
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
/// ACKNOLEDGEMENTS
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

// compilation options
#ifdef YBVH_INLINE
#define YBVH_API inline
#else
#define YBVH_API
#endif

#include <array>
#include <functional>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ybvh {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;
using int4 = std::array<int, 4>;

///
/// BVH scene.
///
struct scene;

///
/// Init scene.
///
YBVH_API scene* make_scene(int nshapes);

///
/// Free scene.
///
/// - paramsters:
///   - scn: scene to free
///
YBVH_API void free_scene(scene* scn);

///
/// Set shape
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///   - npoints/point: point indices
///   - nverts: number of vertices
///   - pos/radius: vertex position and radius
///
YBVH_API void set_point_shape(scene* scn, int sid, const float3x4& frame,
    int npoints, const int* point, int nverts, const float3* pos,
    const float* radius);

///
/// Set shape
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///   - nlines/line: line indices
///   - nverts: number of vertices
///   - pos/radius: vertex position and radius
///
YBVH_API void set_line_shape(scene* scn, int sid, const float3x4& frame,
    int nlines, const int2* lines, int nverts, const float3* pos,
    const float* radius);

///
/// Set shape
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///   - ntriange/triangles: triangle indices
///   - nverts: number of vertices
///   - pos/radius: vertex position and radius
///
YBVH_API void set_triangle_shape(scene* scn, int sid, const float3x4& frame,
    int ntriangles, const int3* triangles, int nverts, const float3* pos,
    const float* radius);

///
/// Set shape
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///   - npoints/point: point indices
///   - ntetra/tetra: tetrahedra indices
///   - nverts: number of vertices
///   - pos/radius: vertex position and radius
///
YBVH_API void set_tetra_shape(scene* scn, int sid, const float3x4& frame,
    int ntetra, const int4* tetra, int nverts, const float3* pos,
    const float* radius);

///
/// Set shape
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///   - nverts: number of vertices
///   - pos/radius: vertex position and radius
///
YBVH_API void set_point_shape(scene* scn, int sid, const float3x4& frame,
    int nverts, const float3* pos, const float* radius);

///
/// Set a shape frame.
///
/// - parameters:
///   - scn: scene
///   - sid: shape id
///   - frame: shape transform
///
YBVH_API void set_shape_frame(scene* scn, int sid, const float3x4& frame);

///
/// Heuristic strategy for bvh build
///
enum struct heuristic_type {
    /// default strategy
    def = 0,
    /// balanced binary tree
    equalnum,
    /// splitting space binary tree
    equalsize,
    /// surface area heuristic (full sweep for accuracy)
    sah,
    /// surface area heuristic (binned for speed)
    binned_sah,
    /// total number of strategies
    htype_max
};

///
/// Builds a scene BVH.
///
/// - parameters:
///   - scn: object to build the bvh for
///   - heuristic: heristic use to build the bvh (htype::none for default)
///   - do_shapes: build shapes
///
YBVH_API void build_bvh(scene* scn, heuristic_type htype = heuristic_type::def,
    bool do_shapes = true);

///
/// Builds a shape BVH.
///
/// - parameters:
///   - scn: object to build the bvh for
///   - sid: required shape
///   - heuristic: heristic use to build the bvh (htype::none for default)
///
YBVH_API void build_bvh(
    scene* scn, int sid, heuristic_type htype = heuristic_type::def);

///
/// Refit the bounds of each shape for moving objects. Use this only to avoid
/// a rebuild, but note that queries are likely slow if objects move a lot.
/// Before calling refit, set the scene data.
///
/// - parameters:
///   - scn: scene to refit
///   - do_shapes: refit shapes
///
YBVH_API void refit_bvh(scene* scn, bool do_shapes = false);

///
/// Refit the bounds of each shape for moving objects. Use this only to avoid
/// a rebuild, but note that queries are likely slow if objects move a lot.
/// Before calling refit, set the scene data.
///
/// - parameters:
///   - scn: scene to refit
///   - sid: required shape
///
YBVH_API void refit_bvh(scene* scn, int sid);

///
/// BVH intersection.
///
struct point {
    /// distance
    float dist = 0;
    /// shape index
    int sid = -1;
    /// element index
    int eid = -1;
    /// element baricentric coordinates
    float4 euv = {0, 0, 0, 0};

    /// Check whether it was a hit.
    operator bool() const { return eid >= 0; }
};

///
/// Intersect the scene with a ray. Find any interstion if early_exit, otherwise
/// find first intersection.
///
/// - parameters:
///   - scn: scene to intersect
///   - sid: shape id
///   - ray: ray
///   - early_exit: whether to stop at the first found hit
///
/// - Returns:
///   - intersection point
///
YBVH_API point intersect_ray(const scene* scn, const float3& ray_o,
    const float3& ray_d, float ray_tmin, float ray_tmax, bool early_exit);
YBVH_API point intersect_ray(const scene* scn, int sid, const float3& ray_o,
    const float3& ray_d, float ray_tmin, float ray_tmax, bool early_exit);

///
/// Returns a list of shape pairs that can possibly overlap by checking only
/// they axis aligned bouds. This is only a conservative check useful for
/// collision detection.
///
/// - parameters:
///   - scn1, scn2: scenes to overlap
///   - conservative: use conservative checks
///   - skip_self: exlude self intersections
///   - skip_duplicates: exlude intersections (s1,s2) if (s2,s1) is already
///     present
///
/// - out parameters:
///   - overlaps: vectors of shape overlaps
///
YBVH_API void overlap_shape_bounds(const scene* scn1, const scene* scn2,
    bool conservative, bool exclude_duplicates, bool exclude_self,
    std::vector<int2>* overlaps);

///
/// Returns a list of all elements of scene1/shape1 that overlap with the
/// vectices of scene2/shape2 within the give radius. Has the option of
/// filtering only the closest intersection.
///
/// - parameters:
///   - scn1, scn2: scenes to overlap
///   - sid1, sid2: shape ids
///   - exclude_self: whether to exlude self intersections
///   - radius: radius of the distance query (see below)
///   - first_only: finds only the closest overlap for each vertex
///
/// - out parameters:
///   - overlaps: vectors of element overlaps (the second argument is (sid,eid))
///   - soverlaps: vectors of shape overlaps
///
YBVH_API void overlap_verts(const scene* scn1, const scene* scn2,
    bool exclude_self, float radius, bool first_only,
    std::vector<int2>* soverlaps,
    std::vector<std::pair<point, int2>>* overlaps);

///
/// Returns a list of all elements of scene1/shape1 that overlap with the
/// vectices of scene2/shape2 within the give radius. Has the option of
/// filtering
/// only the closest intersection.
///
/// - parameters:
///   - scn1, scn2: scenes to overlap
///   - sid1, sid2: shape ids
///   - exclude_self: whether to exlude self intersections
///   - radius: radius of the distance query (see below)
///   - first_only: finds only the closest overlap for each vertex
///
/// - out parameters:
///   - overlaps: vectors of element overlaps (the second argument is (sid,eid))
///   - soverlaps: vectors of shape overlaps
///
YBVH_API void overlap_verts(const scene* scn1, const scene* scn2, int sid1,
    int sid2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, int2>>* overlaps);

///
/// Finds the closest element that overlaps a point within a given radius.
///
/// - parameters:
///   - scene: scene to check
///   - pt: ray origin
///   - max_dist: max point distance
///   - early_exit: whether to stop at the first found hit
///
/// - returns:
///   - overlap point
///
YBVH_API point overlap_point(
    const scene* scn, const float3& pt, float max_dist, bool early_exit);

///
/// Finds the closest element that overlaps a point within a given radius.
///
/// - parameters:
///   - scene: scene to check
///   - pt: ray origin
///   - max_dist: max point distance
///   - early_exit: whether to stop at the first found hit
///
/// - returns:
///   - overlap point
///
YBVH_API point overlap_point(const scene* scn, int sid, const float3& pt,
    float max_dist, bool early_exit);

// -----------------------------------------------------------------------------
// PERFORMANCE TUNING INTERFACE
// -----------------------------------------------------------------------------

///
/// Compute BVH stats.
///
YBVH_API void compute_bvh_stats(const scene* scn, bool include_shapes,
    int& nprims, int& ninternals, int& nleaves, int& min_depth, int& max_depth,
    int req_shape = -1);

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YBVH_INLINE
#include "yocto_bvh.cpp"
#endif

#endif
