///
/// YOCTO_SYM: Simple rigib body simulator with collision support
/// for convex and concate triangle meshes.
///
///
/// USAGE:
///
/// 1. define the rigib body scene
/// - init the scene with make_scene()
/// - for each rigid body, set shape data with set_body()
/// - set collision callbacks with set_overlap_callbacks()
/// 2. start tha simulation with init_simulation()
/// 3. for each frame, advacne the simulation with advance_simulation()
/// 4. after each frame, retrive or change the rigid body frame, with
/// get_body_frame() and set_body_frame(), and the rigid body velocities with
/// get_body_velocity() and set_body_velocity()
/// 5. if desired, you can explicitly compute rigid body moments with
/// compute_moments()
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// Shapes are indexed meshes and are described by array of vertex indices for
/// triangles, and arrays of vertex data.
///
/// The rigib body code performs collision detection and response in gravity.
/// For collision detection, we use only mesh vertices, so increase object
/// tesselation to make the simulation more accurate. This allows us to support
/// convex and concave objects and make the simulation very stable compared to
/// convex collision detection such as GJK or MPR.
///
/// The solver is based on the sequential impulse techniques, more correctly
/// known as Projected Guass-Sidel. Friction is grossly approximated now,
/// waiting for a refactoring before getting better.
///
///
/// COMPILATION:
///
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
///
/// The .cpp file depends on yocto_math.h.
///
///
/// HISTORY:
/// - v 0.10: switch to .h/.cpp pair
/// - v 0.9: doxygen comments
/// - v 0.8: opaque API (allows for changing internals without altering API)
/// - v 0.7: internally use pointers for performance transparency
/// - v 0.6: new formulation for moment computation (and bug fixes)
/// - v 0.5: faster collision detection
/// - v 0.4: [major API change] move to modern C++ interface
/// - v 0.3: removal of C interface
/// - v 0.2: use of STL containers
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace ysym {}

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

#ifndef _YSYM_H_
#define _YSYM_H_

// compilation options
#ifdef YSYM_INLINE
#define YSYM_API inline
#else
#define YSYM_API
#endif

#include <array>
#include <functional>
#include <tuple>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ysym {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using float3x3 = std::array<std::array<float, 3>, 3>;  // column major
using float3x2 = std::array<std::array<float, 3>, 2>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;
using int4 = std::array<int, 4>;

///
/// Simulation scene.
///
struct scene;

///
/// Initialze a scene.
///
/// Parameters:
/// - nbodies: number of rigid bodies
///
YSYM_API scene* make_scene(int nbodies);

///
/// Free a scene.
///
/// Parameters:
/// - scn: rigid body scene
///
YSYM_API void free_scene(scene* scn);

///
/// Set a rigid body.
///
/// Parameters:
/// - scn: scene
/// - bid: body index
/// - frame: body frame
/// - density: bidy density (zero for not-simulated)
/// - lin_vel: linear velocity
/// - ang_vel: angular velocity
/// - ntriangles: number of triangles
/// - triangles: triangles indices
/// - nverts: number of vertices
/// - pos: vertex positions
///
YSYM_API void set_body(scene* scn, int bid, const float3x4& frame,
    const float3& lin_vel, const float3& ang_vel, float density, int ntriangles,
    const int3* triangles, int nverts, const float3* pos);

///
/// Get a rigid body frame.
///
/// Parameters:
/// - scn: scene
/// - bid: body index
///
/// Returns:
/// - body frame
///
YSYM_API float3x4 get_body_frame(const scene* scn, int bid);

///
/// Get a rigid body linear and angular velocity.
///
/// Parameters:
/// - scn: scene
/// - bid: body index
///
/// Returns:
/// - linear and angular velocity
///
YSYM_API float3x2 get_body_velocity(const scene* scn, int bid);

///
/// Sets a rigid body frame.
///
/// Parameters:
/// - scn: scene
/// - bid: body index
/// - frame: rigid body frame
///
YSYM_API void set_body_frame(scene* scn, int bid, const float3x4& frame);

///
/// Sets a rigid body linear and angular velocity.
///
/// Parameters:
/// - scn: scene
/// - bid: body index
/// - lin_vel: linear velocity
/// - ang_vel: angular velocity
///
YSYM_API void set_body_velocity(
    scene* scn, int bid, const float3& lin_vel, const float3& ang_vel);

///
/// Point-scene overlap.
///
struct overlap_point {
    /// overlap distance
    float dist = 0;
    /// shape index
    int sid = -1;
    /// elements index
    int eid = -1;
    /// element baricentric coordinates
    float4 euv = {0, 0, 0, 0};

    /// check whether it was a hit
    operator bool() const { return eid >= 0; }
};

///
/// Shape-shape intersection (conservative)
///
/// Out Parameters:
/// - ctx: context
/// - overlaps: overlaps array
///
using overlap_shapes_cb = void (*)(void* ctx, std::vector<int2>*);

///
/// Closest element intersection callback
///
/// Parameters:
/// - ctx: context
/// - sid: shape to check
/// - pt: point
/// - max_dist: maximum distance
///
/// Return:
/// - overlap point
///
using overlap_shape_cb = overlap_point (*)(
    void* ctx, int sid, const float3& pt, float max_dist);

///
/// Closest vertex-to-element overlap
///
/// Parameters:
/// - ctx: context
/// - sid1: element shape to check
/// - sid2: vertex shape to check
/// - max_dist: maximum distance from each vert
///
/// Out Params:
/// - overlaps: overlapping elements
///
using overlap_verts_cb = void (*)(void* ctx, int sid1, int sid2, float max_dist,
    std::vector<std::pair<overlap_point, int2>>* overlaps);

///
/// Refit data structure after transform updates
///
/// Parameters:
/// - ctx: context
/// - scn: scene
///
using overlap_refit_cb = void (*)(void* ctx, const scene* scn, int nshapes);

///
/// Set overlap functions
///
YSYM_API void set_overlap_callbacks(scene* scn, void* ctx,
    overlap_shapes_cb overlap_shapes, overlap_shape_cb overlap_shape,
    overlap_verts_cb overlap_verts, overlap_refit_cb overlap_refit);

///
/// Computes the moments of a shape.
///
/// Parameters:
/// - triangles: triangle indices
/// - pos: vertex positions
///
/// Output parameters:
/// - volume: volume
/// - center: center of mass
/// - inertia: inertia tensore (wrt center of mass)
///
YSYM_API void compute_moments(int ntriangles, const int3* triangles, int nverts,
    const float3* pos, float* volume, float3* center, float3x3* inertia);

///
/// Computes the moments of a shape.
///
/// Parameters:
/// - triangles: triangle indices
/// - pos: vertex positions
///
/// Output parameters:
/// - volume: volume
/// - center: center of mass
/// - inertia: inertia tensore (wrt center of mass)
///
YSYM_API void compute_moments(int ntetra, const int4* tetra, int nverts,
    const float3* pos, float* volume, float3* center, float3x3* inertia);

///
/// Initialize the simulation
///
/// Paramaters:
/// - scene: rigib body scene
///
YSYM_API void init_simulation(scene* scn);

///
/// Advance the simulation one step at a time.
///
/// Paramaters:
/// - scene: rigib body scene
/// - dt: time step
///
YSYM_API void advance_simulation(scene* scn, float dt);

}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YSYM_INLINE
#include "yocto_sym.cpp"
#endif

#endif
