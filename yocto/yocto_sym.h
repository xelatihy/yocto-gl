//
// YOCTO_SYM: Simple rigib body simulator with collision support
// for convex and concate triangle meshes.
//

//
// USAGE:
//
// 0. include this file (more compilation options below)
// 1. define the rigib body scene
// - init an scene
//   scene = ysym::scene()
// - for each shape, set shape data (frame,vels,density,triangles,pos)
//   foreach shape: scene.shapes.push_back({...})
// - set collision callbacks
//    - can use yocto_bvh for this
// 3. advance the time at each time step with advance
// 4. can look up updated shape state directly from shape array
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by array of vertex indices for
// triangles, and arrays of vertex data.
//
// The rigib body code performs collision detection and response in gravity.
// For collision detection, we use only mesh vertices, so increase object
// tesselation to make the simulation more accurate. This allows us to support
// convex and concave objects and make the simulation very stable compared to
// convex collision detection such as GJK or MPR.
//
// The solver is based on the sequential impulse techniques, more correctly
// known as Projected Guass-Sidel. Friction is grossly approximated now,
// waiting for a refactoring before getting better.
//

//
// COMPILATION:
//
// The library has two APIs. The default one is usable directly from C++,
// while the other is usable from both C and C++. To use from C, compile the
// library into a static or dynamic lib using a C++ and then include/link from
// C using the C API.
//
// All functions in this library are inlined by default for ease of use in C++.
// To use the library as a .h/.cpp pair do the following:
// - to use as a .h, just #define YGL_DECLARATION before including this file
// - to build as a .cpp, just #define YGL_IMPLEMENTATION before including this
// file into only one file that you can either link directly or pack as a lib.
//
// This file depends on yocto_math.h.
//

//
// HISTORY:
// - v 0.4: [major API change] move to modern C++ interface
// - v 0.3: removal of C interface
// - v 0.2: use of STL containers
// - v 0.1: C++ implementation
// - v 0.0: initial release in C99
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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
#ifndef YGL_DECLARATION
#define YGL_API inline
#else
#define YGL_API
#endif

#include <functional>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

namespace ysym {

//
// Using directives
//
using namespace ym;

//
// Rigid shape
//
struct shape {
    // shape configuration ------------------------
    frame3f frame;  // rigid transform
    vec3f lin_vel;  // linear velocity
    vec3f ang_vel;  // angular velocity

    // physical properties ------------------------
    float density = 1;      // density
    bool simulated = true;  // simulated

    // shape data ----------------------------------
    array_view<vec3i> triangles;  // triangles
    array_view<vec3f> pos;        // vertex positions

    // [private] computed values ------------------
    float _mass = 1;                        // mass
    mat3f _inertia_local = identity_mat3f;  // shape inertia
    vec3f _centroid_local = zero3f;         // shape center
    vec3f _centroid_world = zero3f;         // shape center
    float _mass_inv = 1;                    // mass inverse
    mat3f _inertia_inv_world =
        identity_mat3f;  // inverse of inertia tensor (world-space)
    mat3f _inertia_inv_local =
        identity_mat3f;  // inverse of inertia tensor (local-space)
};

//
// Collision point and response [private]
//
struct _collision {
    vec2i shapes = zero2i;                           // shapes
    frame3f frame = identity_frame3f;                // collision frame
    vec3f impulse = zero3f, local_impulse = zero3f;  // impulses
    vec3f vel_before = zero3f, vel_after = zero3f;   // velocities (for viz)
    vec3f meff_inv = zero3f;                         // effective mass
    float depth = 0;                                 // penetration depth
};

//
// Point-scene overlap.
//
struct overlap_point {
    float dist = 0;      // ray distance
    int sid = -1;        // shape index
    int eid = -1;        // element index
    vec3f euv = zero3f;  // element baricentric coordinates

    // check whether it was a hit
    operator bool() const { return eid >= 0; }
};

//
// Shape-shape intersection (conservative)
//
// Out Parameters:
// - overlaps: overlaps array
//
// Return:
// - number of intersections
//
using overlap_shapes_cb = function<vector<vec2i>()>;

//
// Closest element intersection callback
//
// Parameters:
// - sid: shape to check
// - pt: point
// - max_dist: maximum distance
//
// Out Parameters:
// - dist: distance
// - eid: element id
// - euv: element uv
//
// Return:
// - whether we intersect or not
//
using overlap_shape_cb =
    function<overlap_point(int sid, const vec3f& pt, float max_dist)>;

//
// Refit data structure after transform updates
//
// Parameters:
// - xform: transform array
//
using overlap_refit_cb = function<void()>;

//
// Rigid body scene
//
struct scene {
    // simulation shapes -----------------------
    vector<shape> shapes;  // shapes

    // global simulation values ----------------
    vec3f gravity = {0.f, -9.82f, 0.f};  // gravity
    float lin_drag = 0.01;               // linear drag
    float ang_drag = 0.01;               // angular drag
    int iterations = 20;                 // solver iterations

    // overlap callbacks -----------------------
    overlap_shapes_cb overlap_shapes;  // overlap callbacks
    overlap_shape_cb overlap_shape;    // overlap callbacks
    overlap_refit_cb overlap_refit;    // overlap callbacks

    // [private] collision data ----------------
    vector<_collision> _collisions;  // [private] collisions
    vector<vec2i> _shapecollisions;  // [private] shape collisions
};

//
// Computes the moments of a shape.
//
// Parameters:
// - nelems: number of elements
// - elem: elements
// - etype: element type
// - nverts: number of vertices
// - pos: vertex positions
//
// Output parameters:
// - volume: volume
// - center: center of mass
// - inertia: inertia tensore (wrt center of mass)
//
YGL_API void compute_moments(const array_view<vec3i>& triangles,
                             const array_view<vec3f>& pos, float& volume,
                             vec3f& center, mat3f& inertia);

//
// Initialize the simulation
//
// Paramaters:
// - scene: rigib body scene
//
YGL_API void init_simulation(scene& scn);

//
// Advance the simulation one step at a time.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YGL_API void advance_simulation(scene& scn, float dt);

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

namespace ysym {

//
// Compute moments. Public API, see above.
//
// Implementation notes:
// - follow eberly
// http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
// - re-implement this using
// https://github.com/melax/sandbox/blob/master/include/geometric.h
//
YGL_API void compute_moments(const array_view<vec3i>& triangles,
                             const array_view<vec3f>& pos, float& volume,
                             vec3f& center, mat3f& inertia) {
    // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
    float integral[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (auto f : triangles) {
        auto v0 = pos[f[0]], v1 = pos[f[1]], v2 = pos[f[2]];

        // Get cross product of edges and normal vector.
        auto V1mV0 = v1 - v0;
        auto V2mV0 = v2 - v0;
        vec3f N = cross(V1mV0, V2mV0);

        // Compute integral terms.
        auto tmp0 = v0[0] + v1[0];
        auto f1x = tmp0 + v2[0];
        auto tmp1 = v0[0] * v0[0];
        auto tmp2 = tmp1 + v1[0] * tmp0;
        auto f2x = tmp2 + v2[0] * f1x;
        auto f3x = v0[0] * tmp1 + v1[0] * tmp2 + v2[0] * f2x;
        auto g0x = f2x + v0[0] * (f1x + v0[0]);
        auto g1x = f2x + v1[0] * (f1x + v1[0]);
        auto g2x = f2x + v2[0] * (f1x + v2[0]);

        tmp0 = v0[1] + v1[1];
        auto f1y = tmp0 + v2[1];
        tmp1 = v0[1] * v0[1];
        tmp2 = tmp1 + v1[1] * tmp0;
        auto f2y = tmp2 + v2[1] * f1y;
        auto f3y = v0[1] * tmp1 + v1[1] * tmp2 + v2[1] * f2y;
        auto g0y = f2y + v0[1] * (f1y + v0[1]);
        auto g1y = f2y + v1[1] * (f1y + v1[1]);
        auto g2y = f2y + v2[1] * (f1y + v2[1]);

        tmp0 = v0[2] + v1[2];
        auto f1z = tmp0 + v2[2];
        tmp1 = v0[2] * v0[2];
        tmp2 = tmp1 + v1[2] * tmp0;
        auto f2z = tmp2 + v2[2] * f1z;
        auto f3z = v0[2] * tmp1 + v1[2] * tmp2 + v2[2] * f2z;
        auto g0z = f2z + v0[2] * (f1z + v0[2]);
        auto g1z = f2z + v1[2] * (f1z + v1[2]);
        auto g2z = f2z + v2[2] * (f1z + v2[2]);

        // Update integrals.
        integral[0] += N[0] * f1x / 6;
        integral[1] += N[0] * f2x / 24;
        integral[2] += N[1] * f2y / 24;
        integral[3] += N[2] * f2z / 24;
        integral[4] += N[0] * f3x / 60;
        integral[5] += N[1] * f3y / 60;
        integral[6] += N[2] * f3z / 60;
        integral[7] += N[0] * (v0[1] * g0x + v1[1] * g1x + v2[1] * g2x) / 120;
        integral[8] += N[1] * (v0[2] * g0y + v1[2] * g1y + v2[2] * g2y) / 120;
        integral[9] += N[2] * (v0[0] * g0z + v1[0] * g1z + v2[0] * g2z) / 120;
    }

    // mass
    volume = integral[0];

    // center of mass
    center = vec3f(integral[1], integral[2], integral[3]) / volume;

    // inertia relative to center
    inertia[0][0] = integral[5] + integral[6] -
                    volume * (center[1] * center[1] + center[2] * center[2]);
    inertia[0][1] = -integral[7] + volume * center[0] * center[1];
    inertia[0][2] = -integral[9] + volume * center[2] * center[0];
    inertia[1][0] = inertia[0][1];
    inertia[1][1] = integral[4] + integral[6] -
                    volume * (center[2] * center[2] + center[0] * center[0]);
    inertia[1][2] = -integral[8] + volume * center[1] * center[2];
    inertia[2][0] = inertia[0][2];
    inertia[2][1] = inertia[1][1];
    inertia[2][2] = integral[4] + integral[5] -
                    volume * (center[0] * center[0] + center[1] * center[1]);
}

//
// Compute collisions.
//
static inline void _compute_collision(scene& scn, const vec2i& shapes) {
    auto& shape1 = scn.shapes[shapes[0]];
    auto& shape2 = scn.shapes[shapes[1]];
    for (auto& p2 : shape2.pos) {
        auto p = transform_point(shape2.frame, p2);
        auto overlap =
            scn.overlap_shape(shapes[0], p, std::numeric_limits<float>::max());
        if (!overlap) continue;
        auto triangle = shape1.triangles[overlap.eid];
        auto v0 = shape1.pos[triangle[0]], v1 = shape1.pos[triangle[1]],
             v2 = shape1.pos[triangle[2]];
        auto tp = transform_point(shape1.frame, blerp(v0, v1, v2, overlap.euv));
        auto n = transform_direction(shape1.frame, cross(v1 - v0, v2 - v0));
        const auto eps = -0.01f;
        auto ptp = normalize(p - tp);
        if (dot(n, ptp) > eps) continue;
        scn._collisions.push_back(_collision());
        auto& col = scn._collisions.back();
        col.shapes = shapes;
        col.depth = overlap.dist;
        col.frame = make_frame3(p, n);
    }
}

//
// Compute collisions.
//
static inline void _compute_collisions(scene& scene) {
    // check which shapes might overlap
    scene._shapecollisions = scene.overlap_shapes();
    // test all pair-wise objects
    scene._collisions.resize(0);
    for (auto& sc : scene._shapecollisions) {
        if (!scene.shapes[sc[0]].simulated && !scene.shapes[sc[1]].simulated)
            continue;
        if (scene.shapes[sc[0]].triangles.empty()) continue;
        if (scene.shapes[sc[1]].triangles.empty()) continue;
        _compute_collision(scene, sc);
    }
}

//
// Apply an impulse where the position is relative to the center of mass.
//
static inline void _apply_rel_impulse(shape& shp, const vec3f impulse,
                                      const vec3f rel_pos) {
    if (!shp.simulated) return;
    shp.lin_vel += impulse * shp._mass_inv;
    shp.ang_vel += shp._inertia_inv_world * cross(rel_pos, impulse);
}

//
// Shortcut math function.
//
static inline float _muldot(const vec3f& v, const mat3f& m) {
    return dot(v, m * v);
}

//
// Solve constraints with PGS.
//
YGL_API void _solve_constraints(scene& scn, float dt) {
    // initialize computation
    for (auto& col : scn._collisions) {
        col.local_impulse = zero3f;
        col.impulse = zero3f;
        auto& shape1 = scn.shapes[col.shapes[0]];
        auto& shape2 = scn.shapes[col.shapes[1]];
        auto r1 = col.frame.o() - shape1._centroid_world,
             r2 = col.frame.o() - shape2._centroid_world;
        col.meff_inv = {
            1 / (shape1._mass_inv + shape2._mass_inv +
                 _muldot(cross(r1, col.frame[0]), shape1._inertia_inv_world) +
                 _muldot(cross(r2, col.frame[0]), shape2._inertia_inv_world)),
            1 / (shape1._mass_inv + shape2._mass_inv +
                 _muldot(cross(r1, col.frame[1]), shape1._inertia_inv_world) +
                 _muldot(cross(r2, col.frame[1]), shape2._inertia_inv_world)),
            1 / (shape1._mass_inv + shape2._mass_inv +
                 _muldot(cross(r1, col.frame[2]), shape1._inertia_inv_world) +
                 _muldot(cross(r2, col.frame[2]), shape2._inertia_inv_world))};
    }

    // compute relative velocity for visualization
    for (auto& col : scn._collisions) {
        auto& shape1 = scn.shapes[col.shapes[0]];
        auto& shape2 = scn.shapes[col.shapes[1]];
        auto r1 = col.frame.o() - shape1._centroid_world,
             r2 = col.frame.o() - shape2._centroid_world;
        auto v1 = shape1.lin_vel + cross(shape1.ang_vel, r1),
             v2 = shape2.lin_vel + cross(shape2.ang_vel, r2);
        col.vel_before = v2 - v1;
    }

    // solve constraints
    for (int i = 0; i < scn.iterations; i++) {
        for (auto& col : scn._collisions) {
            auto& shape1 = scn.shapes[col.shapes[0]];
            auto& shape2 = scn.shapes[col.shapes[1]];
            auto r1 = col.frame.o() - shape1._centroid_world,
                 r2 = col.frame.o() - shape2._centroid_world;
            auto v1 = shape1.lin_vel + cross(shape1.ang_vel, r1),
                 v2 = shape2.lin_vel + cross(shape2.ang_vel, r2);
            auto vr = v2 - v1;
            _apply_rel_impulse(shape1, col.impulse, r1);
            _apply_rel_impulse(shape2, -col.impulse, r2);
            // float offset = col.depth*0.8f/dt;
            auto offset = 0.0f;
            vec3f local_impulse =
                col.meff_inv * vec3f{-dot(col.frame[0], vr),
                                     -dot(col.frame[1], vr),
                                     -dot(col.frame[2], vr) + offset};
            col.local_impulse += local_impulse;
            col.local_impulse[2] = clamp(col.local_impulse[2], 0.0f,
                                         std::numeric_limits<float>::max());
            col.local_impulse[0] =
                clamp(col.local_impulse[0], -col.local_impulse[2] * 0.6f,
                      col.local_impulse[2] * 0.6f);
            col.local_impulse[1] =
                clamp(col.local_impulse[1], -col.local_impulse[2] * 0.6f,
                      col.local_impulse[2] - offset * 0.6f);
            col.impulse = col.local_impulse[2] * col.frame[2] +
                          col.local_impulse[0] * col.frame[0] +
                          col.local_impulse[1] * col.frame[1];
            _apply_rel_impulse(shape1, -col.impulse, r1);
            _apply_rel_impulse(shape2, col.impulse, r2);
        }
    }

    // compute relative velocity for visualization
    for (auto& col : scn._collisions) {
        auto& shape1 = scn.shapes[col.shapes[0]];
        auto& shape2 = scn.shapes[col.shapes[1]];
        auto r1 = col.frame.o() - shape1._centroid_world,
             r2 = col.frame.o() - shape2._centroid_world;
        auto v1 = shape1.lin_vel + cross(shape1.ang_vel, r1),
             v2 = shape2.lin_vel + cross(shape2.ang_vel, r2);
        col.vel_after = v2 - v1;
    }

    // recompute total impulse and velocity for visualization
    for (auto& col : scn._collisions) {
        col.impulse = col.local_impulse[2] * col.frame[2] +
                      col.local_impulse[0] * col.frame[0] +
                      col.local_impulse[1] * col.frame[1];
    }
}

//
// Initialize the simulation
//
YGL_API void init_simulation(scene& scn) {
    for (auto& shp : scn.shapes) {
        if (shp.simulated) {
            float volume = 1;
            ysym::compute_moments(shp.triangles, shp.pos, volume,
                                  shp._centroid_local, shp._inertia_local);
            shp._mass = shp.density * volume;
            shp._centroid_world =
                transform_point(shp.frame, shp._centroid_local);
            shp._mass_inv = 1 / shp._mass;
            shp._inertia_inv_local = inverse(shp._inertia_local);
        } else {
            shp._mass = 0;
            shp._mass_inv = 0;
            shp._centroid_local = zero3f;
            shp._centroid_world = zero3f;
            shp._inertia_local = mat3f(zero3f, zero3f, zero3f);
            shp._inertia_inv_local = mat3f(zero3f, zero3f, zero3f);
        }
    }
}

//
// Check function for numerical problems
//
inline bool _isfinite(const vec3f& v) {
    return isfinite(v[0]) && isfinite(v[1]) && isfinite(v[2]);
}

//
// Advance simulation. Public API, see above.
//
YGL_API void advance_simulation(scene& scn, float dt) {
    // update centroid and inertia
    for (auto& shp : scn.shapes) {
        if (!shp.simulated) continue;
        shp._centroid_world = transform_point(shp.frame, shp._centroid_local);
        shp._inertia_inv_world =
            shp.frame.m() * shp._inertia_inv_local * transpose(shp.frame.m());
    }

    // compute collisions
    _compute_collisions(scn);

    // apply external forces
    vec3f gravity_impulse = scn.gravity * dt;
    for (auto& shp : scn.shapes) {
        if (!shp.simulated) continue;
        shp.lin_vel += gravity_impulse;
    }

    // solve constraints
    _solve_constraints(scn, dt);

    // apply drag
    for (auto& shp : scn.shapes) {
        if (!shp.simulated) continue;
        shp.lin_vel *= 1 - scn.lin_drag;
        shp.ang_vel *= 1 - scn.ang_drag;
    }

    // update position and velocity
    for (auto& shp : scn.shapes) {
        if (!shp.simulated) continue;

        // check for nans
        if (!_isfinite(shp.frame.o())) printf("nan detected\n");
        if (!_isfinite(shp.lin_vel)) printf("nan detected\n");
        if (!_isfinite(shp.ang_vel)) printf("nan detected\n");

        // translate the frame to the centroid
        auto centroid = shp.frame.m() * shp._centroid_local + shp.frame.o();
        // update centroid
        centroid += shp.lin_vel * dt;
        float angle = length(shp.ang_vel) * dt;
        if (angle) {
            vec3f axis = normalize(shp.ang_vel);
            shp.frame.m() = rotation_mat3(axis, angle) * shp.frame.m();
            // TODO: if using matrices, I gotta orthonormalize them
        }
        // translate the frame back
        shp.frame.o() = centroid - shp.frame.m() * shp._centroid_local;
    }

    // update acceleartion for collisions
    scn.overlap_refit();
}

}  // namespace

#endif

#endif
