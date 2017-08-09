///
/// # Yocto/Sym
///
/// Simple rigib body simulator with collision support for convex and concave
/// triangle meshes.
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
/// This library depends in yocto_math.h Optionally depend on yocto_bvh.h/.cpp
/// for internal acceleration. Disable this by setting YSYM_NO_BVH.
///
///
/// ## Usage for Scene Creation
///
/// 1. create a scene with `make_scene()`
/// 2. add materials with `add_rigid_material()`
/// 3. add shapes with `add_rigid_shape()`
/// 4. add instances with `add_instance()`
///
/// ## Usage for Simulation
///
/// 1. either build the point-overlap acceleration structure with
///   `init_overlap()` or supply your own with `set_overlap_callbacks()`
/// 2. prepare simuation internal data `init_simulation()`
/// 3. define simulation params with the `simulation_params` structure
/// 4. advance the simiulation with `advance_simulation()`
/// 5. set or get rigid body frames and velocities with `set_XXX()` and
/// `get_XXX()` methods
///
/// ## History
///
/// - v 0.16: simpler logging
/// - v 0.15: removal of group overlap
/// - v 0.14: use yocto_math in the interface and remove inline compilation
/// - v 0.13: move to add api
/// - v 1.12: internally use yocto_bvh if desired
/// - v 0.11: added instances
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

#include <array>
#include <functional>
#include <tuple>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

///
/// Simple rigid body simulator
///
namespace ysym {

///
/// Simulation scene.
///
struct scene;

///
/// Initialze a scene.
///
scene* make_scene();

///
/// Free a scene.
///
/// - Parameters:
///     - scn: rigid body scene
///
void free_scene(scene*& scn);

///
/// Sets a rigid shape.
///
/// - Parameters:
///     - scn: scene
///     - ntriangles: number of triangles
///     - triangles: triangles indices
///     - nverts: number of vertices
///     - pos: vertex positions
///
int add_rigid_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos);

///
/// Sets a material.
///
/// - Parameters:
///     - scn: scene
///     - density: body density (zero for not-simulated)
///
int add_rigid_material(scene* scn, float density);

///
/// Set a rigid body.
///
/// - Parameters:
///     - scn: scene
///     - frame: body frame
///     - sid: shape id
///     - mid: material id
///     - lin_vel: linear velocity
///     - ang_vel: angular velocity
///
int add_rigid_body(scene* scn, const ym::frame3f& frame, int sid, int mid,
    const ym::vec3f& lin_vel = {0, 0, 0}, const ym::vec3f& ang_vel = {0, 0, 0});

///
/// Get a rigid body frame.
///
/// - Parameters:
///     - scn: scene
///     - bid: body index
/// - Returns:
///     - body frame
///
ym::frame3f get_rigid_body_frame(const scene* scn, int bid);

///
/// Get a rigid body linear and angular velocity.
///
/// - Parameters:
///     - scn: scene
///     - bid: body index
/// - Returns:
///     - linear and angular velocity
///
std::pair<ym::vec3f, ym::vec3f> get_rigid_body_velocity(
    const scene* scn, int bid);

///
/// Sets a rigid body frame.
///
/// - Parameters:
///     - scn: scene
///     - bid: body index
///     - frame: rigid body frame
///
void set_rigid_body_frame(scene* scn, int bid, const ym::frame3f& frame);

///
/// Sets a rigid body linear and angular velocity.
///
/// - Parameters:
///     - scn: scene
///     - bid: body index
///     - lin_vel: linear velocity
///     - ang_vel: angular velocity
///
void set_rigid_body_velocity(
    scene* scn, int bid, const ym::vec3f& lin_vel, const ym::vec3f& ang_vel);

///
/// Point-scene overlap.
///
struct overlap_point {
    /// overlap distance
    float dist = 0;
    /// instance index
    int iid = -1;
    /// shape index
    int sid = -1;
    /// elements index
    int eid = -1;
    /// element baricentric coordinates
    ym::vec4f euv = {0, 0, 0, 0};

    /// check whether it was a hit
    operator bool() const { return eid >= 0; }
};

///
/// Shape-shape intersection (conservative)
///
/// - Out Parameters:
///     - overlaps: overlaps array
///
using overlap_shapes_cb = std::function<void(std::vector<ym::vec2i>*)>;

///
/// Closest element intersection callback
///
/// - Parameters:
///     - sid: shape to check
///     - pt: point
///     - max_dist: maximum distance
/// - Return:
///     - overlap point
///
using overlap_shape_cb =
    std::function<overlap_point(int sid, const ym::vec3f& pt, float max_dist)>;

///
/// Refit data structure after transform updates
///
/// - Parameters:
///     - scn: scene
///
using overlap_refit_cb = std::function<void(const scene* scn, int nshapes)>;

///
/// Set overlap functions
///
void set_overlap_callbacks(scene* scn, overlap_shapes_cb overlap_shapes,
    overlap_shape_cb overlap_shape, overlap_refit_cb overlap_refit);

///
/// Initialize overlap functions using internal structures.
///
void init_overlap(scene* scn);

///
/// Computes the moments of a shape.
///
/// - Parameters:
///     - triangles: triangle indices
///     - pos: vertex positions
/// - Output Parameters:
///     - volume: volume
///     - center: center of mass
///     - inertia: inertia tensore (wrt center of mass)
///
void compute_moments(int ntriangles, const ym::vec3i* triangles, int nverts,
    const ym::vec3f* pos, float* volume, ym::vec3f* center, ym::mat3f* inertia);

///
/// Computes the moments of a shape.
///
/// - Parameters:
///     - triangles: triangle indices
///     - pos: vertex positions
/// - Output Parameters:
///     - volume: volume
///     - center: center of mass
///     - inertia: inertia tensore (wrt center of mass)
///
void compute_moments(int ntetra, const ym::vec4i* tetra, int nverts,
    const ym::vec3f* pos, float* volume, ym::vec3f* center, ym::mat3f* inertia);

///
/// Initialize the simulation.
///
/// - Paramaters:
///     - scene: rigib body scene
///
void init_simulation(scene* scn);

///
/// Simulation Parameters.
///
struct simulation_params {
    /// delta time
    float dt = 1 / 60.0f;
    /// gravity
    ym::vec3f gravity = {0.f, -9.82f, 0.f};
    /// solver iterations
    int solver_iterations = 20;
    /// global linear velocity drag
    float lin_drag = 0.01;
    /// global angular velocity drag
    float ang_drag = 0.01;
};

///
/// Advance the simulation one step at a time.
///
/// - Paramaters:
///     - scene: rigib body scene
///     - dt: time step
///
void advance_simulation(scene* scn, const simulation_params& params);

}  // namespace ysym

#endif
