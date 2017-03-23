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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR YOCTO_SYM
// -----------------------------------------------------------------------------

#include "yocto_sym.h"
#include "yocto_math.h"

#include <iostream>

//
// TODO: cleanup: frame -> pos/rot
// TODO: cleanup: shape -> body/shape/material
//

namespace ysym {

//
// Rigid shape.
//
struct shape {
    // shape configuration ------------------------
    ym::frame3f frame;  // rigid transform
    ym::vec3f lin_vel;  // linear velocity
    ym::vec3f ang_vel;  // angular velocity

    // physical properties ------------------------
    float density = 1;      // density
    bool simulated = true;  // simulated

    // shape data ----------------------------------
    int nelems = 0;                        // numner of elements
    const ym::vec3i* triangles = nullptr;  // triangles
    int nverts = 0;                        // number of vertices
    const ym::vec3f* pos = nullptr;        // vertex positions

    // [private] computed values ------------------
    float _mass = 1;                                // mass
    ym::mat3f _inertia_local = ym::identity_mat3f;  // shape inertia
    ym::vec3f _centroid_local = ym::zero3f;         // shape center
    ym::vec3f _centroid_world = ym::zero3f;         // shape center
    float _mass_inv = 1;                            // mass inverse
    ym::mat3f _inertia_inv_world =
        ym::identity_mat3f;  // inverse of inertia tensor (world-space)
    ym::mat3f _inertia_inv_local =
        ym::identity_mat3f;  // inverse of inertia tensor (local-space)
};

//
// Collision point and response [private]
//
struct collision {
    ym::vec2i shapes = ym::zero2i;             // shapes
    ym::frame3f frame = ym::identity_frame3f;  // collision frame
    ym::vec3f impulse = ym::zero3f, local_impulse = ym::zero3f;  // impulses
    ym::vec3f vel_before = ym::zero3f,
              vel_after = ym::zero3f;  // velocities (for viz)
    ym::vec3f meff_inv = ym::zero3f;   // effective mass
    float depth = 0;                   // penetration depth
};

//
// Rigid body scene
//
struct scene {
    // simulation shapes -----------------------
    std::vector<shape*> shapes;  // shapes

    // global simulation values ----------------
    ym::vec3f gravity = {0.f, -9.82f, 0.f};  // gravity
    float lin_drag = 0.01;                   // linear drag
    float ang_drag = 0.01;                   // angular drag
    int iterations = 20;                     // solver iterations

    // overlap callbacks -----------------------
    float overlap_max_radius = 0.25;  // maximum vertex overlap distance
    void* overlap_ctx = nullptr;      // overlap context
    overlap_shapes_cb overlap_shapes = nullptr;  // overlap callbacks
    overlap_shape_cb overlap_shape = nullptr;    // overlap callbacks
    overlap_verts_cb overlap_verts = nullptr;    // overlap callbacks
    overlap_refit_cb overlap_refit = nullptr;    // overlap callbacks

    // overlap data used for visualization [private] ----
    std::vector<collision> __collisions;

    // destructor
    ~scene();
};

//
// Public API.
//
YSYM_API scene* make_scene(int nbodies) {
    auto scn = new scene();
    scn->shapes.resize(nbodies);
    for (auto& shp : scn->shapes) shp = new shape();
    return scn;
}

//
// Public API.
//
YSYM_API void free_scene(scene* scn) {
    if (scn) delete scn;
}

//
// Destructor
//
YSYM_API scene::~scene() {
    for (auto shp : shapes)
        if (shp) delete shp;
    shapes.clear();
}

//
// Public API.
//
YSYM_API void set_body(scene* scn, int bid, const float3x4& frame,
    const float3& lin_vel, const float3& ang_vel, float density, int ntriangles,
    const int3* triangles, int nverts, const float3* pos) {
    scn->shapes[bid]->frame = frame;
    scn->shapes[bid]->lin_vel = lin_vel;
    scn->shapes[bid]->ang_vel = ang_vel;
    scn->shapes[bid]->density = density;
    scn->shapes[bid]->simulated = density > 0;
    scn->shapes[bid]->nelems = ntriangles;
    scn->shapes[bid]->triangles = (const ym::vec3i*)triangles;
    scn->shapes[bid]->nverts = nverts;
    scn->shapes[bid]->pos = (const ym::vec3f*)pos;
}

//
// Public API.
//
YSYM_API float3x4 get_body_frame(const scene* scn, int bid) {
    return scn->shapes[bid]->frame;
}

//
// Public API.
//
YSYM_API float3x2 get_body_velocity(const scene* scn, int bid) {
    return {scn->shapes[bid]->lin_vel, scn->shapes[bid]->ang_vel};
}

//
// Public API.
//
YSYM_API void set_body_frame(scene* scn, int bid, const float3x4& frame) {
    scn->shapes[bid]->frame = frame;
}

//
// Public API.
//
YSYM_API void set_body_velocity(
    scene* scn, int bid, const float3& lin_vel, const float3& ang_vel) {
    scn->shapes[bid]->lin_vel = lin_vel;
    scn->shapes[bid]->ang_vel = ang_vel;
}

//
// Public API.
//
YSYM_API void set_overlap_callbacks(scene* scn, void* ctx,
    overlap_shapes_cb overlap_shapes, overlap_shape_cb overlap_shape,
    overlap_verts_cb overlap_verts, overlap_refit_cb overlap_refit) {
    scn->overlap_ctx = ctx;
    scn->overlap_shapes = overlap_shapes;
    scn->overlap_shape = overlap_shape;
    scn->overlap_verts = overlap_verts;
    scn->overlap_refit = overlap_refit;
}

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
static inline ym::mat3f _compute_tetra_inertia(const ym::vec3f& v0,
    const ym::vec3f& v1, const ym::vec3f& v2, const ym::vec3f& v3,
    const ym::vec3f& center) {
    // volume
    auto volume = ym::tetrahedron_volume(v0, v1, v2, v3);
    // relative vertices
    auto vr0 = v0 - center, vr1 = v1 - center, vr2 = v2 - center,
         vr3 = v3 - center;
    // diagonal elements
    auto diag = ym::zero3f;  // x^2, y^2, z^2
    for (auto j = 0; j < 3; j++) {
        diag[j] = (vr0[j] * vr0[j] + vr1[j] * vr1[j] + vr2[j] * vr2[j] +
                      vr3[j] * vr3[j] + vr0[j] * vr1[j] + vr0[j] * vr2[j] +
                      vr0[j] * vr3[j] + vr1[j] * vr2[j] + vr1[j] * vr3[j] +
                      vr2[j] * vr3[j]) *
                  6 * volume / 60;
    }
    // off-diagonal elements
    auto offd = ym::zero3f;  // x*y, x*z, y*z
    for (auto j = 0; j < 3; j++) {
        auto j1 = (j + 1) % 3, j2 = (j + 2) % 3;
        offd[j] =
            (2 * vr0[j1] * vr0[j2] + 2 * vr1[j1] * vr1[j2] +
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
// Public API.
//
YSYM_API std::tuple<float, ym::vec3f, ym::mat3f> _compute_moments(
    int ntriangles, const ym::vec3i* triangles, int nverts,
    const ym::vec3f* pos) {
    // volume and center
    auto volume = 0.0f;
    auto center = ym::zero3f;
    for (auto i = 0; i < ntriangles; i++) {
        auto t = triangles[i];
        auto tvolume =
            ym::tetrahedron_volume(ym::zero3f, pos[t[0]], pos[t[1]], pos[t[2]]);
        volume += tvolume;
        center +=
            tvolume * (ym::zero3f + pos[t[0]] + pos[t[1]] + pos[t[2]]) / 4.0f;
    }
    center /= volume;
    // inertia
    auto inertia = ym::mat3f{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto i = 0; i < ntriangles; i++) {
        auto t = triangles[i];
        inertia += _compute_tetra_inertia(
            ym::zero3f, pos[t[0]], pos[t[1]], pos[t[2]], center);
    }
    inertia /= volume;
    return std::make_tuple(volume, center, inertia);
}

//
// Compute moments. Public API, see above.
//
YSYM_API std::tuple<float, ym::vec3f, ym::mat3f> _compute_moments(
    int ntetra, const ym::vec4i* tetra, int nverts, const ym::vec3f* pos) {
    // volume and center
    auto volume = 0.0f;
    auto center = ym::zero3f;
    for (auto i = 0; i < ntetra; i++) {
        auto t = tetra[i];
        auto tvolume =
            ym::tetrahedron_volume(pos[t[0]], pos[t[1]], pos[t[2]], pos[t[3]]);
        volume += tvolume;
        center +=
            tvolume * (pos[t[0]] + pos[t[1]] + pos[t[2]] + pos[t[3]]) / 4.0f;
    }
    center /= volume;
    // inertia
    auto inertia = ym::mat3f{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (auto i = 0; i < ntetra; i++) {
        auto t = tetra[i];
        inertia += _compute_tetra_inertia(
            pos[t[0]], pos[t[1]], pos[t[2]], pos[t[3]], center);
    }
    inertia /= volume;
    return std::make_tuple(volume, center, inertia);
}

//
// Public API.
//
YSYM_API std::tuple<float, float3, float3x3> compute_moments(
    int ntriangles, const int3* triangles, int nverts, const float3* pos) {
    return _compute_moments(
        ntriangles, (const ym::vec3i*)triangles, nverts, (const ym::vec3f*)pos);
}

//
// Public API.
//
YSYM_API std::tuple<float, float3, float3x3> compute_moments(
    int ntetra, const int4* tetra, int nverts, const float3* pos) {
    return _compute_moments(
        ntetra, (const ym::vec4i*)tetra, nverts, (const ym::vec3f*)pos);
}

//
// Compute collisions.
//
#if 1
static inline void _compute_collision(const scene* scn, const ym::vec2i& sids,
    std::vector<collision>* collisions) {
    std::vector<std::pair<overlap_point, int2>> overlaps;
    scn->overlap_verts(
        scn->overlap_ctx, sids[0], sids[1], scn->overlap_max_radius, &overlaps);
    if (overlaps.empty()) return;
    auto shape1 = scn->shapes[sids[0]];
    auto shape2 = scn->shapes[sids[1]];
    for (auto& overlap : overlaps) {
        auto p = transform_point(shape2->frame, shape2->pos[overlap.second[1]]);
        auto triangle = shape1->triangles[overlap.first.eid];
        auto v0 = shape1->pos[triangle[0]], v1 = shape1->pos[triangle[1]],
             v2 = shape1->pos[triangle[2]];
        auto tp = transform_point(
            shape1->frame, ym::blerp(v0, v1, v2, overlap.first.euv));
        auto n =
            transform_direction(shape1->frame, ym::triangle_normal(v0, v1, v2));
        const auto eps = -0.01f;
        auto ptp = ym::normalize(p - tp);
        if (ym::dot(n, ptp) > eps) continue;
        auto col = collision();
        col.shapes = sids;
        col.depth = overlap.first.dist;
        col.frame = ym::make_frame3_fromz(p, n);
        collisions->push_back(col);
    }
}
#else
static inline void _compute_collision(const scene* scn, const ym::vec2i& sids,
    std::vector<collision>* collisions) {
    auto shape1 = scn->shapes[sids[0]];
    auto shape2 = scn->shapes[sids[1]];
    for (auto& p2 : shape2->pos) {
        auto p = transform_point(shape2->frame, p2);
        auto overlap = scn->overlap_shape(sids[0], p, scn->overlap_max_radius);
        if (!overlap) continue;
        auto triangle = shape1->triangles[overlap.eid];
        auto v0 = shape1->pos[triangle[0]], v1 = shape1->pos[triangle[1]],
             v2 = shape1->pos[triangle[2]];
        auto tp =
            transform_point(shape1->frame, blerp(v0, v1, v2, overlap.euv));
        auto n =
            transform_direction(shape1->frame, triangle_normal(v0, v1, v2));
        const auto eps = -0.01f;
        auto ptp = normalize(p - tp);
        if (dot(n, ptp) > eps) continue;
        auto col = collision();
        col.shapes = sids;
        col.depth = overlap.dist;
        col.frame = make_frame3(p, n);
        collisions.push_back(col);
    }
}
#endif

//
// Compute collisions.
//
static inline void _compute_collisions(
    scene* scene, std::vector<collision>* collisions) {
    // check which shapes might overlap
    auto shapecollisions = std::vector<int2>();
    scene->overlap_shapes(scene->overlap_ctx, &shapecollisions);
    // test all pair-wise objects
    collisions->clear();
    for (auto& sc : shapecollisions) {
        if (!scene->shapes[sc[0]]->simulated &&
            !scene->shapes[sc[1]]->simulated)
            continue;
        if (!scene->shapes[sc[0]]->triangles) continue;
        if (!scene->shapes[sc[1]]->triangles) continue;
        _compute_collision(scene, sc, collisions);
        _compute_collision(scene, {sc[1], sc[0]}, collisions);
    }
}

//
// Apply an impulse where the position is relative to the center of mass.
//
static inline void _apply_rel_impulse(
    shape* shp, const ym::vec3f impulse, const ym::vec3f rel_pos) {
    if (!shp->simulated) return;
    shp->lin_vel += impulse * shp->_mass_inv;
    shp->ang_vel += shp->_inertia_inv_world * ym::cross(rel_pos, impulse);
}

//
// Shortcut math function.
//
static inline float _muldot(const ym::vec3f& v, const ym::mat3f& m) {
    return ym::dot(v, m * v);
}

//
// Solve constraints with PGS.
//
YSYM_API void _solve_constraints(
    scene* scn, std::vector<collision>& collisions, float dt) {
    // initialize computation
    for (auto& col : collisions) {
        col.local_impulse = ym::zero3f;
        col.impulse = ym::zero3f;
        auto shape1 = scn->shapes[col.shapes[0]];
        auto shape2 = scn->shapes[col.shapes[1]];
        auto r1 = ym::pos(col.frame) - shape1->_centroid_world,
             r2 = ym::pos(col.frame) - shape2->_centroid_world;
        col.meff_inv = {1 / (shape1->_mass_inv + shape2->_mass_inv +
                                _muldot(ym::cross(r1, col.frame[0]),
                                    shape1->_inertia_inv_world) +
                                _muldot(ym::cross(r2, col.frame[0]),
                                    shape2->_inertia_inv_world)),
            1 / (shape1->_mass_inv + shape2->_mass_inv +
                    _muldot(ym::cross(r1, col.frame[1]),
                        shape1->_inertia_inv_world) +
                    _muldot(ym::cross(r2, col.frame[1]),
                        shape2->_inertia_inv_world)),
            1 / (shape1->_mass_inv + shape2->_mass_inv +
                    _muldot(ym::cross(r1, col.frame[2]),
                        shape1->_inertia_inv_world) +
                    _muldot(ym::cross(r2, col.frame[2]),
                        shape2->_inertia_inv_world))};
    }

    // compute relative velocity for visualization
    for (auto& col : collisions) {
        auto shape1 = scn->shapes[col.shapes[0]];
        auto shape2 = scn->shapes[col.shapes[1]];
        auto r1 = ym::pos(col.frame) - shape1->_centroid_world,
             r2 = ym::pos(col.frame) - shape2->_centroid_world;
        auto v1 = shape1->lin_vel + ym::cross(shape1->ang_vel, r1),
             v2 = shape2->lin_vel + ym::cross(shape2->ang_vel, r2);
        col.vel_before = v2 - v1;
    }

    // solve constraints
    for (int i = 0; i < scn->iterations; i++) {
        for (auto& col : collisions) {
            auto shape1 = scn->shapes[col.shapes[0]];
            auto shape2 = scn->shapes[col.shapes[1]];
            auto r1 = ym::pos(col.frame) - shape1->_centroid_world,
                 r2 = ym::pos(col.frame) - shape2->_centroid_world;
            auto v1 = shape1->lin_vel + ym::cross(shape1->ang_vel, r1),
                 v2 = shape2->lin_vel + ym::cross(shape2->ang_vel, r2);
            auto vr = v2 - v1;
            _apply_rel_impulse(shape1, col.impulse, r1);
            _apply_rel_impulse(shape2, -col.impulse, r2);
            // float offset = col.depth*0.8f/dt;
            auto offset = 0.0f;
            ym::vec3f local_impulse =
                col.meff_inv * ym::vec3f{-ym::dot(col.frame[0], vr),
                                   -ym::dot(col.frame[1], vr),
                                   -ym::dot(col.frame[2], vr) + offset};
            col.local_impulse += local_impulse;
            col.local_impulse[2] = ym::clamp(
                col.local_impulse[2], 0.0f, std::numeric_limits<float>::max());
            col.local_impulse[0] = ym::clamp(col.local_impulse[0],
                -col.local_impulse[2] * 0.6f, col.local_impulse[2] * 0.6f);
            col.local_impulse[1] =
                ym::clamp(col.local_impulse[1], -col.local_impulse[2] * 0.6f,
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
        auto shape1 = scn->shapes[col.shapes[0]];
        auto shape2 = scn->shapes[col.shapes[1]];
        auto r1 = ym::pos(col.frame) - shape1->_centroid_world,
             r2 = ym::pos(col.frame) - shape2->_centroid_world;
        auto v1 = shape1->lin_vel + ym::cross(shape1->ang_vel, r1),
             v2 = shape2->lin_vel + ym::cross(shape2->ang_vel, r2);
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
YSYM_API void init_simulation(scene* scn) {
    for (auto shp : scn->shapes) {
        if (shp->simulated) {
            float volume = 1;
            std::tie(volume, shp->_centroid_local, shp->_inertia_local) =
                ysym::_compute_moments(
                    shp->nelems, shp->triangles, shp->nverts, shp->pos);
            shp->_mass = shp->density * volume;
            shp->_centroid_world =
                transform_point(shp->frame, shp->_centroid_local);
            shp->_mass_inv = 1 / shp->_mass;
            shp->_inertia_inv_local = ym::inverse(shp->_inertia_local);
        } else {
            shp->_mass = 0;
            shp->_mass_inv = 0;
            shp->_centroid_local = ym::zero3f;
            shp->_centroid_world = ym::zero3f;
            shp->_inertia_local = ym::mat3f{ym::zero3f, ym::zero3f, ym::zero3f};
            shp->_inertia_inv_local =
                ym::mat3f{ym::zero3f, ym::zero3f, ym::zero3f};
        }
    }
}

//
// Check function for numerical problems
//
inline bool _isfinite(const ym::vec3f& v) {
    return std::isfinite(v[0]) && std::isfinite(v[1]) && std::isfinite(v[2]);
}

//
// Advance simulation. Public API, see above.
//
YSYM_API void advance_simulation(scene* scn, float dt) {
    // update centroid and inertia
    for (auto shp : scn->shapes) {
        if (!shp->simulated) continue;
        shp->_centroid_world =
            transform_point(shp->frame, shp->_centroid_local);
        shp->_inertia_inv_world = ym::rot(shp->frame) *
                                  shp->_inertia_inv_local *
                                  ym::transpose(ym::rot(shp->frame));
    }

    // compute collisions
    auto collisions = std::vector<collision>();
    _compute_collisions(scn, &collisions);

    // apply external forces
    ym::vec3f gravity_impulse = scn->gravity * dt;
    for (auto shp : scn->shapes) {
        if (!shp->simulated) continue;
        shp->lin_vel += gravity_impulse;
    }

    // solve constraints
    _solve_constraints(scn, collisions, dt);

    // copy for visualization
    scn->__collisions = collisions;

    // apply drag
    for (auto shp : scn->shapes) {
        if (!shp->simulated) continue;
        shp->lin_vel *= 1 - scn->lin_drag;
        shp->ang_vel *= 1 - scn->ang_drag;
    }

    // update position and velocity
    for (auto shp : scn->shapes) {
        if (!shp->simulated) continue;

        // check for nans
        if (!_isfinite(ym::pos(shp->frame))) printf("nan detected\n");
        if (!_isfinite(shp->lin_vel)) printf("nan detected\n");
        if (!_isfinite(shp->ang_vel)) printf("nan detected\n");

        // translate the frame to the centroid
        auto centroid =
            ym::rot(shp->frame) * shp->_centroid_local + ym::pos(shp->frame);
        // update centroid
        centroid += shp->lin_vel * dt;
        float angle = ym::length(shp->ang_vel) * dt;
        if (angle) {
            ym::vec3f axis = ym::normalize(shp->ang_vel);
            ym::rot(shp->frame) =
                ym::rotation_mat3(axis, angle) * ym::rot(shp->frame);
            // TODO: if using matrices, I gotta orthonormalize them
        }
        // translate the frame back
        ym::pos(shp->frame) =
            centroid - ym::rot(shp->frame) * shp->_centroid_local;
    }

    // update acceleartion for collisions
    scn->overlap_refit(scn->overlap_ctx, scn, (int)scn->shapes.size());
}

}  // namespace
