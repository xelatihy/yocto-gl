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
// - v 0.6: new formulation for moment computation (and bug fixes)
// - v 0.5: faster collision detection
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
// Rigid shape.
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
// Point-scene overlap.
//
struct overlap_point {
    float dist = 0;      // ray distance
    int sid = -1;        // shape index
    int eid = -1;        // element index
    vec4f euv = zero4f;  // element baricentric coordinates

    // check whether it was a hit
    operator bool() const { return eid >= 0; }
};

//
// Shape-shape intersection (conservative)
//
// Out Parameters:
// - overlaps: overlaps array
//
using overlap_shapes_cb = function<void(vector<vec2i>&)>;

//
// Closest element intersection callback
//
// Parameters:
// - sid: shape to check
// - pt: point
// - max_dist: maximum distance
//
// Return:
// - overlap point
//
using overlap_shape_cb =
    function<overlap_point(int sid, const vec3f& pt, float max_dist)>;

//
// Closest vertex-to-element overlap
//
// Parameters:
// - sid1: element shape to check
// - sid2: vertex shape to check
// - max_dist: maximum distance from each vert
//
// Out Params:
// - overlaps: overlapping elements
//
using overlap_verts_cb =
    function<void(int sid1, int sid2, float max_dist,
                  vector<pair<overlap_point, vec2i>>& overlaps)>;

//
// Refit data structure after transform updates
//
// Parameters:
// - xform: transform array
//
using overlap_refit_cb = function<void()>;

//
// Collision point and response [private]
//
struct collision {
    vec2i shapes = zero2i;                           // shapes
    frame3f frame = identity_frame3f;                // collision frame
    vec3f impulse = zero3f, local_impulse = zero3f;  // impulses
    vec3f vel_before = zero3f, vel_after = zero3f;   // velocities (for viz)
    vec3f meff_inv = zero3f;                         // effective mass
    float depth = 0;                                 // penetration depth
};

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
    float overlap_max_radius = 0.25;   // maximum vertex overlap distance
    overlap_shapes_cb overlap_shapes;  // overlap callbacks
    overlap_shape_cb overlap_shape;    // overlap callbacks
    overlap_verts_cb overlap_verts;    // overlap callbacks
    overlap_refit_cb overlap_refit;    // overlap callbacks

    // overlap data used for visualization [private] ----
    vector<collision> __collisions;
};

//
// Computes the moments of a shape.
//
// Parameters:
// - triangles: triangle indices
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
// Computes the moments of a shape.
//
// Parameters:
// - tetra: tetrahedra indices
// - pos: vertex positions
//
// Output parameters:
// - volume: volume
// - center: center of mass
// - inertia: inertia tensore (wrt center of mass)
//
YGL_API void compute_moments(const array_view<vec4i>& tetra,
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

#include <iostream>

namespace ysym {

//
// Computes the tetrahedra moment of inertia
//
// Notes:
// - from "Explicit Exact Formulas for the 3-D Tetrahedron Inertia Tensor
// in Terms of its Vertex Coordinates" by F. Tonon published in Journal of
// Mathematics and Statistics 1 (1), 2004.
// - implemented similarly to
// https://github.com/melax/sandbox/blob/master/include/geometric.h
//
static inline mat3f _compute_tetra_inertia(const vec3f& v0, const vec3f& v1,
                                           const vec3f& v2, const vec3f& v3,
                                           const vec3f& center) {
    // volume
    auto volume = tetrahedron_volume(v0, v1, v2, v3);
    // relative vertices
    auto vr0 = v0 - center, vr1 = v1 - center, vr2 = v2 - center,
         vr3 = v3 - center;
    // diagonal elements
    auto diag = zero3f;  // x^2, y^2, z^2
    for (auto j = 0; j < 3; j++) {
        diag[j] = (vr0[j] * vr0[j] + vr1[j] * vr1[j] + vr2[j] * vr2[j] +
                   vr3[j] * vr3[j] + vr0[j] * vr1[j] + vr0[j] * vr2[j] +
                   vr0[j] * vr3[j] + vr1[j] * vr2[j] + vr1[j] * vr3[j] +
                   vr2[j] * vr3[j]) *
                  6 * volume / 60;
    }
    // off-diagonal elements
    auto offd = zero3f;  // x*y, x*z, y*z
    for (auto j = 0; j < 3; j++) {
        auto j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        offd[j] = (2 * vr0[j1] * vr0[j2] + 2 * vr1[j1] * vr1[j2] +
                   2 * vr2[j1] * vr2[j2] + 2 * vr3[j1] * vr3[j2] +
                   vr1[j1] * vr0[j2] + vr2[j1] * vr0[j2] + vr3[j1] * vr0[j2] +
                   vr0[j1] * vr1[j2] + vr2[j1] * vr1[j2] + vr3[j1] * vr1[j2] +
                   vr0[j1] * vr2[j2] + vr1[j1] * vr2[j2] + vr3[j1] * vr2[j2] +
                   vr0[j1] * vr3[j2] + vr1[j1] * vr3[j2] + vr2[j1] * vr3[j2]) *
                  6 * volume / 120;
    }
    // setup inertia
    return {{diag[1] + diag[2], -offd[2], -offd[1]},
            {-offd[2], diag[0] + diag[2], -offd[0]},
            {-offd[1], -offd[0], diag[0] + diag[1]}};
}

//
// Compute moments. Public API, see above.
//
YGL_API void compute_moments(const array_view<vec3i>& triangles,
                             const array_view<vec3f>& pos, float& volume,
                             vec3f& center, mat3f& inertia) {
    // volume and center
    volume = 0;
    center = zero3f;
    for (auto t : triangles) {
        auto tvolume =
            tetrahedron_volume(zero3f, pos[t[0]], pos[t[1]], pos[t[2]]);
        volume += tvolume;
        center += tvolume * (zero3f + pos[t[0]] + pos[t[1]] + pos[t[2]]) / 4;
    }
    center /= volume;
    // inertia
    inertia = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto t : triangles) {
        inertia += _compute_tetra_inertia(zero3f, pos[t[0]], pos[t[1]],
                                          pos[t[2]], center);
    }
    inertia /= volume;
}

//
// Compute moments. Public API, see above.
//
YGL_API void compute_moments(const array_view<vec4i>& tetra,
                             const array_view<vec3f>& pos, float& volume,
                             vec3f& center, mat3f& inertia) {
    // volume and center
    volume = 0;
    center = zero3f;
    for (auto t : tetra) {
        auto tvolume =
            tetrahedron_volume(pos[t[0]], pos[t[1]], pos[t[2]], pos[t[3]]);
        volume += tvolume;
        center += tvolume * (pos[t[0]] + pos[t[1]] + pos[t[2]] + pos[t[3]]) / 4;
    }
    center /= volume;
    // inertia
    inertia = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto t : tetra) {
        inertia += _compute_tetra_inertia(pos[t[0]], pos[t[1]], pos[t[2]],
                                          pos[t[3]], center);
    }
    inertia /= volume;
}

//
// Compute collisions.
//
#if 1
static inline void _compute_collision(scene& scn, const vec2i& shapes,
                                      vector<collision>& collisions) {
    vector<pair<overlap_point, vec2i>> overlaps;
    scn.overlap_verts(shapes[0], shapes[1], scn.overlap_max_radius, overlaps);
    if (overlaps.empty()) return;
    auto& shape1 = scn.shapes[shapes[0]];
    auto& shape2 = scn.shapes[shapes[1]];
    for (auto& overlap : overlaps) {
        auto p = transform_point(shape2.frame, shape2.pos[overlap.second[1]]);
        auto triangle = shape1.triangles[overlap.first.eid];
        auto v0 = shape1.pos[triangle[0]], v1 = shape1.pos[triangle[1]],
             v2 = shape1.pos[triangle[2]];
        auto tp =
            transform_point(shape1.frame, blerp(v0, v1, v2, overlap.first.euv));
        auto n = transform_direction(shape1.frame, triangle_normal(v0, v1, v2));
        const auto eps = -0.01f;
        auto ptp = normalize(p - tp);
        if (dot(n, ptp) > eps) continue;
        collisions.push_back(collision());
        auto& col = collisions.back();
        col.shapes = shapes;
        col.depth = overlap.first.dist;
        col.frame = make_frame3(p, n);
    }
}
#else
static inline void _compute_collision(scene& scn, const vec2i& shapes,
                                      vector<collision>& collisions) {
    auto& shape1 = scn.shapes[shapes[0]];
    auto& shape2 = scn.shapes[shapes[1]];
    for (auto& p2 : shape2.pos) {
        auto p = transform_point(shape2.frame, p2);
        auto overlap = scn.overlap_shape(shapes[0], p, scn.overlap_max_radius);
        if (!overlap) continue;
        auto triangle = shape1.triangles[overlap.eid];
        auto v0 = shape1.pos[triangle[0]], v1 = shape1.pos[triangle[1]],
             v2 = shape1.pos[triangle[2]];
        auto tp = transform_point(shape1.frame, blerp(v0, v1, v2, overlap.euv));
        auto n = transform_direction(shape1.frame, triangle_normal(v0, v1, v2));
        const auto eps = -0.01f;
        auto ptp = normalize(p - tp);
        if (dot(n, ptp) > eps) continue;
        collisions.push_back(collision());
        auto& col = collisions.back();
        col.shapes = shapes;
        col.depth = overlap.dist;
        col.frame = make_frame3(p, n);
    }
}
#endif

//
// Compute collisions.
//
static inline void _compute_collisions(scene& scene,
                                       vector<collision>& collisions) {
    // check which shapes might overlap
    auto shapecollisions = vector<vec2i>();
    scene.overlap_shapes(shapecollisions);
    // test all pair-wise objects
    collisions.clear();
    for (auto& sc : shapecollisions) {
        if (!scene.shapes[sc[0]].simulated && !scene.shapes[sc[1]].simulated)
            continue;
        if (scene.shapes[sc[0]].triangles.empty()) continue;
        if (scene.shapes[sc[1]].triangles.empty()) continue;
        _compute_collision(scene, sc, collisions);
        _compute_collision(scene, {sc[1], sc[0]}, collisions);
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
YGL_API void _solve_constraints(scene& scn, vector<collision>& collisions,
                                float dt) {
    // initialize computation
    for (auto& col : collisions) {
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
    for (auto& col : collisions) {
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
        for (auto& col : collisions) {
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
    for (auto& col : collisions) {
        auto& shape1 = scn.shapes[col.shapes[0]];
        auto& shape2 = scn.shapes[col.shapes[1]];
        auto r1 = col.frame.o() - shape1._centroid_world,
             r2 = col.frame.o() - shape2._centroid_world;
        auto v1 = shape1.lin_vel + cross(shape1.ang_vel, r1),
             v2 = shape2.lin_vel + cross(shape2.ang_vel, r2);
        col.vel_after = v2 - v1;
    }

    // recompute total impulse and velocity for visualization
    for (auto& col : collisions) {
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
    auto collisions = vector<collision>();
    _compute_collisions(scn, collisions);

    // apply external forces
    vec3f gravity_impulse = scn.gravity * dt;
    for (auto& shp : scn.shapes) {
        if (!shp.simulated) continue;
        shp.lin_vel += gravity_impulse;
    }

    // solve constraints
    _solve_constraints(scn, collisions, dt);

    // copy for visualization
    scn.__collisions = collisions;

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
