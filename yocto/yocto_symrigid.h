//
// YOCTO_SYMRIGID: Simple rigib body simulator for with collision support
// for convex and concate triangle meshes.
//

//
// USAGE FOR READING:
//
// 0. include this file (more compilation options below)
// 1. create a rigib body scene
//   - init the scene with ysr_init_scene
//   - for each rigib body:
//     - compute volume, center of mass and inertia with
//     ysr_compute_moments(shape params, moments)
//     - scale volume and inertia by density
//     - substract the center of mass from vertex positions
//     - set the rigib body data ysr_set_body(scene, index, data)
// 2. set collision callbacks with ysr_set_collision
//    - can use yocto_bvh for this
// 3. advance the time at each time step with ysr_advance
// 4. get the updated transforms with ysr_get_transforms
// 5. cleanup with ysr_free_scene
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles, quads),
// and arrays for vertex positions. Only triangle meshes are supported now,
// but others might be added later so we keep here the generic interface.
// Quads are likely to removed soon though.
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

#ifndef _YSR_H_
#define _YSR_H_

// compilation options
#ifdef __cplusplus
#ifndef YGL_DECLARATION
#define YGL_API inline
#define YGLC_API inline
#else
#define YGL_API
#define YGLC_API extern "C"
#endif
#include "yocto_math.h"
#endif

#ifndef __cplusplus
#define YGLC_API extern
#include <stdbool.h>
#endif

// -----------------------------------------------------------------------------
// C++ INTERFACE
// -----------------------------------------------------------------------------

#ifdef __cplusplus

//
// Shape element types
//
enum {
    ysr_etype_point = 1,     // points
    ysr_etype_line = 2,      // lines
    ysr_etype_triangle = 3,  // triangle
    ysr_etype_quad = 4,      // quads
};

//
// Shape-shape intersection (conservative)
//
// Parameters:
// - ctx: pointer to an object passed back to the function
//
// Out Parameters:
// - overlaps: overlaps array
// - overlap capacity
//
// Return:
// - number of intersections
//
typedef int (*ysr_overlap_shapes)(void* ctx, ym_vector<ym_vec2i>* overlaps);

//
// Closest element intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
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
typedef bool (*ysr_overlap_shape)(void* ctx, int sid, const ym_vec3f& pt,
                                  float max_dist, float* dist, int* eid,
                                  ym_vec2f* euv);

//
// Refit data structure after transform updates
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - xform: transform array
//
typedef void (*ysr_overlap_refit)(void* ctx, const ym_frame3f* xform);

//
// Scene forward declaration
//
typedef struct ysr_scene ysr_scene;

//
// Initialize a scene.
//
// Parameters:
// - nbodies: number of rigid bodies
//
// Return:
// - scene
//
YGL_API ysr_scene* ysr_init_scene(int nbodies);

//
// Free a scene.
//
// Parameters:
// - scene: scene to clear
//
YGL_API void ysr_free_scene(ysr_scene* scene);

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
YGL_API void ysr_compute_moments(int nelems, const int* elem, int etype,
                                 int nverts, const ym_vec3f* pos, float* volume,
                                 ym_vec3f* center, ym_mat3f* inertia);

//
// Set rigid body shape.
//
// Parameters:
// - scene: scene
// - bid: body to set
// - xform: transform
// - mass: mass (0 to not simulate object)
// - inertia: interia tensor wrt center of mass
// - nelems: number of elements
// - elem: elements
// - etype: element type
// - nverts: number of vertices
// - pos: vertex position wrt center of mass
//
YGL_API void ysr_set_body(ysr_scene* scene, int bid, const ym_frame3f& xform,
                          float mass, const ym_mat3f& inertia, int nelems,
                          const int* elem, int etype, int nverts,
                          const ym_vec3f* pos);

//
// Get rigib body transform.
//
// Paramaters:
// - scene: rigib body scene
// - bid: body index
//
// Output parameters:
// - xform: transform
//
YGL_API ym_frame3f ysr_get_transform(const ysr_scene* scene, int bid);

//
// Set body transform.
//
// Paramaters:
// - scene: rigib body scene
// - bid: body index
// - xform: transform
//
YGL_API void ysr_set_transform(ysr_scene* scene, int bid,
                               const ym_frame3f& xform);

//
// Set collision callbacks.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YGL_API void ysr_set_collision(ysr_scene* scene, void* ctx,
                               ysr_overlap_shapes overlaps,
                               ysr_overlap_shape overlap,
                               ysr_overlap_refit refit);

//
// Advance the simulation one step at a time.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YGL_API void ysr_advance(ysr_scene* scene, float dt);

#endif

// -----------------------------------------------------------------------------
// C/C++ INTERFACE
// -----------------------------------------------------------------------------

//
// Shape-shape intersection (conservative)
//
// Parameters:
// - ctx: pointer to an object passed back to the function
//
// Out Parameters:
// - overlaps: overlaps array
// - overlap capacity
//
// Return:
// - number of intersections
//
typedef int (*ysrc_overlap_shapes)(void* ctx, int** overlaps, int* ocapacity);

//
// Closest element intersection callback
//
// Parameters:
// - ctx: pointer to an object passed back to the function
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
typedef bool (*ysrc_overlap_shape)(void* ctx, int sid, const float pt[3],
                                   float max_dist, float* dist, int* eid,
                                   float* euv);

//
// Refit data structure after transform updates
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - xform: transform array
//
typedef void (*ysrc_overlap_refit)(void* ctx, const float* xform[16]);

//
// Initialize a scene.
//
// Parameters:
// - nbodies: number of rigid bodies
//
// Return:
// - scene
//
YGLC_API ysr_scene* ysrc_init_scene(int nbodies);

//
// Free a scene.
//
// Parameters:
// - scene: scene to clear
//
YGLC_API void ysrc_free_scene(ysr_scene* scene);

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
YGLC_API void ysrc_compute_moments(int nelems, const int* elem, int etype,
                                   int nverts, const float* pos, float* volume,
                                   float center[3], float inertia[3]);

//
// Set rigid body shape.
//
// Parameters:
// - scene: scene
// - bid: body to set
// - xform: transform
// - mass: mass (0 to not simulate object)
// - inertia: interia tensor wrt center of mass
// - nelems: number of elements
// - elem: elements
// - etype: element type
// - nverts: number of vertices
// - pos: vertex position wrt center of mass
//
YGLC_API void ysrc_set_body(ysr_scene* scene, int bid, const float xform[16],
                            float mass, const float inertia[9], int nelems,
                            const int* elem, int etype, int nverts,
                            const float* pos);

//
// Get rigib body transform.
//
// Paramaters:
// - scene: rigib body scene
// - bid: body index
//
// Output parameters:
// - xform: transform
//
YGLC_API void ysrc_get_transform(const ysr_scene* scene, int bid,
                                 float xform[16]);

//
// Set body transform.
//
// Paramaters:
// - scene: rigib body scene
// - bid: body index
// - xform: transform
//
YGLC_API void ysrc_set_transform(ysr_scene* scene, int bid,
                                 const float xform[16]);

//
// Set collision callbacks.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YGLC_API void ysrc_set_collision(ysr_scene* scene, void* ctx,
                                 ysr_overlap_shapes overlaps,
                                 ysr_overlap_shape overlap,
                                 ysr_overlap_refit refit);

//
// Advance the simulation one step at a time.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YGLC_API void ysrc_advance(ysr_scene* scene, float dt);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if defined(__cplusplus) &&                                                    \
    (!defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION))

//
// Rigid body properties
//
struct ysr__body {
    // world space quantities
    float mass = 1;      // mass
    float mass_inv = 1;  // mass inverse
    ym_mat3f inertia_inv =
        ym_identity_mat3f;  // inverse of inertis tensor (world-space)
    ym_mat3f inertia_inv_local =
        ym_identity_mat3f;  // inverse of inertis tensor (local-space)

    // rigid body shape
    struct {
        int nelems = 0;             // number of elemets
        const int* elem = nullptr;  // elements
        int etype = 0;              // element type
        int nverts = 0;             // number of vertices
        ym_vec3f* pos = nullptr;    // vertex positions
    } shape;
};

//
// Collision point and response
//
struct ysr__collision {
    int bi = 0, bj = 0;                                       // bodies
    ym_frame3f frame = ym_identity_frame3f;                   // collision frame
    ym_vec3f impulse = ym_zero3f, local_impulse = ym_zero3f;  // impulses
    ym_vec3f vel_before = ym_zero3f,
             vel_after = ym_zero3f;  // velocities (for viz)
    ym_vec3f meff_inv = ym_zero3f;   // effective mass
    float depth = 0;                 // penetration depth
};

//
// Rigid body scene
//
struct ysr_scene {
    ym_vector<ym_frame3f> frame;  // rigid body rigid transform
    ym_vector<ym_vec3f> lin_vel;  // linear velocity
    ym_vector<ym_vec3f> ang_vel;  // angular velocity

    ym_vector<ysr__body> bodies;           // bodies
    ym_vector<ysr__collision> collisions;  // collisions
    ym_vector<ym_vec2i> shapecollisions;   // shape collisions
    ym_vec3f gravity;                      // gravity
    float lin_drag, ang_drag;              // linear and angular drag
    int iterations;                        // solver iterations
    ysr_overlap_shapes overlap_shapes;     // overlap callbacks
    ysr_overlap_shape overlap_shape;       // overlap callbacks
    ysr_overlap_refit overlap_refit;       // overlap callbacks
    void* overlap_ctx;                     // overlap callbacks
};

//
// Init scene. Public API, see above.
//
YGL_API ysr_scene* ysr_init_scene(int nbodies) {
    ysr_scene* scene = new ysr_scene();
    scene->bodies.resize(nbodies);
    scene->frame.resize(nbodies);
    scene->lin_vel.resize(nbodies);
    scene->ang_vel.resize(nbodies);
    scene->gravity = ym_vec3f(0.f, -9.82f, 0.f);
    scene->lin_drag = 0.01;
    scene->ang_drag = 0.01;
    scene->iterations = 20;
    return scene;
}

//
// Free scene. Public API, see above.
//
YGL_API void ysr_free_scene(ysr_scene* scene) { delete scene; }

//
// Compute moments. Public API, see above.
//
// Implementation notes:
// - follow eberly
// http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
// - re-implement this using
// https://github.com/melax/sandbox/blob/master/include/geometric.h
//
YGL_API void ysr_compute_moments(int nelems, const int* elem, int etype,
                                 int nverts, const ym_vec3f* pos, float* volume,
                                 ym_vec3f* center, ym_mat3f* inertia) {
    // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
    float integral[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0};

    for (int i = 0; i < nelems; i++) {
        const int* f = elem + i * 3;
        ym_vec3f v0 = pos[f[0]], v1 = pos[f[1]], v2 = pos[f[2]];

        // Get cross product of edges and normal vector.
        ym_vec3f V1mV0 = v1 - v0;
        ym_vec3f V2mV0 = v2 - v0;
        ym_vec3f N = ym_cross(V1mV0, V2mV0);

        // Compute integral terms.
        float tmp0, tmp1, tmp2;
        float f1x, f2x, f3x, g0x, g1x, g2x;
        tmp0 = v0.x + v1.x;
        f1x = tmp0 + v2.x;
        tmp1 = v0.x * v0.x;
        tmp2 = tmp1 + v1.x * tmp0;
        f2x = tmp2 + v2.x * f1x;
        f3x = v0.x * tmp1 + v1.x * tmp2 + v2.x * f2x;
        g0x = f2x + v0.x * (f1x + v0.x);
        g1x = f2x + v1.x * (f1x + v1.x);
        g2x = f2x + v2.x * (f1x + v2.x);

        float f1y, f2y, f3y, g0y, g1y, g2y;
        tmp0 = v0.y + v1.y;
        f1y = tmp0 + v2.y;
        tmp1 = v0.y * v0.y;
        tmp2 = tmp1 + v1.y * tmp0;
        f2y = tmp2 + v2.y * f1y;
        f3y = v0.y * tmp1 + v1.y * tmp2 + v2.y * f2y;
        g0y = f2y + v0.y * (f1y + v0.y);
        g1y = f2y + v1.y * (f1y + v1.y);
        g2y = f2y + v2.y * (f1y + v2.y);

        float f1z, f2z, f3z, g0z, g1z, g2z;
        tmp0 = v0.z + v1.z;
        f1z = tmp0 + v2.z;
        tmp1 = v0.z * v0.z;
        tmp2 = tmp1 + v1.z * tmp0;
        f2z = tmp2 + v2.z * f1z;
        f3z = v0.z * tmp1 + v1.z * tmp2 + v2.z * f2z;
        g0z = f2z + v0.z * (f1z + v0.z);
        g1z = f2z + v1.z * (f1z + v1.z);
        g2z = f2z + v2.z * (f1z + v2.z);

        // Update integrals.
        integral[0] += N.x * f1x / 6;
        integral[1] += N.x * f2x / 24;
        integral[2] += N.y * f2y / 24;
        integral[3] += N.z * f2z / 24;
        integral[4] += N.x * f3x / 60;
        integral[5] += N.y * f3y / 60;
        integral[6] += N.z * f3z / 60;
        integral[7] += N.x * (v0.y * g0x + v1.y * g1x + v2.y * g2x) / 120;
        integral[8] += N.y * (v0.z * g0y + v1.z * g1y + v2.z * g2y) / 120;
        integral[9] += N.z * (v0.x * g0z + v1.x * g1z + v2.x * g2z) / 120;
    }

    // mass
    *volume = integral[0];

    // center of mass
    *center = ym_vec3f(integral[1], integral[2], integral[3]) / *volume;

    // inertia relative to center
    inertia->data()[0] =
        integral[5] + integral[6] -
        (*volume) * (center->y * center->y + center->z * center->z);
    inertia->data()[1] = -integral[7] + (*volume) * center->x * center->y;
    inertia->data()[2] = -integral[9] + (*volume) * center->z * center->x;
    inertia->data()[3] = inertia->data()[1];
    inertia->data()[4] =
        integral[4] + integral[6] -
        (*volume) * (center->z * center->z + center->x * center->x);
    inertia->data()[5] = -integral[8] + (*volume) * center->y * center->z;
    inertia->data()[6] = inertia->data()[2];
    inertia->data()[7] = inertia->data()[4];
    inertia->data()[8] =
        integral[4] + integral[5] -
        (*volume) * (center->x * center->x + center->y * center->y);
}

//
// Set rigid body. Public API, see above.
//
YGL_API void ysr_set_body(ysr_scene* scene, int bid, const ym_frame3f& xform,
                          float mass, const ym_mat3f& inertia, int nelems,
                          const int* elem, int etype, int nverts,
                          const ym_vec3f* pos) {
    ysr__body* body = &scene->bodies[bid];
    body->mass = mass;
    body->mass_inv = (mass) ? 1 / mass : 0;
    body->inertia_inv_local = (mass)
                                  ? ym_inverse(inertia)
                                  : ym_mat3f{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    body->shape.nverts = nverts;
    body->shape.pos = (ym_vec3f*)pos;
    body->shape.nelems = nelems;
    body->shape.etype = etype;
    body->shape.elem = elem;
    ysr_set_transform(scene, bid, xform);
}

//
// Set collision. Public API, see above.
//
YGL_API void ysr_set_collision(ysr_scene* scene, void* ctx,
                               ysr_overlap_shapes overlaps,
                               ysr_overlap_shape overlap,
                               ysr_overlap_refit refit) {
    scene->overlap_shapes = overlaps;
    scene->overlap_shape = overlap;
    scene->overlap_refit = refit;
    scene->overlap_ctx = ctx;
}

//
// Get transform. Public API, see above.
//
YGL_API ym_frame3f ysr_get_transform(const ysr_scene* scene, int bid) {
    return scene->frame[bid];
}

//
// Set transform. Public API, see above.
//
YGL_API void ysr_set_transform(ysr_scene* scene, int bid,
                               const ym_frame3f& xform) {
    scene->frame[bid] = xform;
}

//
// Compute collisions.
//
static inline void ysr__compute_collision(ysr_scene* scene, int bi, int bj) {
    ysr__body* b1 = &scene->bodies[bi];
    ysr__body* b2 = &scene->bodies[bj];
    for (int j = 0; j < b2->shape.nverts; j++) {
        ym_vec3f p = ym_transform_point(scene->frame[bj], b2->shape.pos[j]);
        int eid;
        float dist;
        ym_vec2f euv;
        bool hit = scene->overlap_shape(scene->overlap_ctx, bi, p, ym_max_float,
                                        &dist, &eid, &euv);
        if (!hit) continue;
        const int* f = b1->shape.elem + eid * 3;
        ym_vec3f v0 = b1->shape.pos[f[0]], v1 = b1->shape.pos[f[1]],
                 v2 = b1->shape.pos[f[2]];
        ym_vec3f tp = ym_transform_point(scene->frame[bi],
                                         ym_blerp(v0, v1, v2, euv[0], euv[1]));
        ym_vec3f n = ym_transform_direction(scene->frame[bi],
                                            ym_cross(v1 - v0, v2 - v0));
        const float eps = -0.01f;
        ym_vec3f ptp = ym_normalize(p - tp);
        if (ym_dot(n, ptp) > eps) continue;
        scene->collisions.push_back(ysr__collision());
        ysr__collision* col = &scene->collisions.back();
        col->bi = bi;
        col->bj = bj;
        col->depth = dist;
        col->frame = ym_make_frame(p, n);
    }
}

//
// Compute collisions.
//
static inline void ysr__compute_collisions(ysr_scene* scene) {
    // check which shapes might overlap
    scene->overlap_shapes(scene->overlap_ctx, &scene->shapecollisions);
    // test all pair-wise objects
    scene->collisions.resize(0);
    for (int c = 0; c < scene->shapecollisions.size(); c++) {
        int i = scene->shapecollisions[c].x, j = scene->shapecollisions[c].y;
        if (!scene->bodies[i].mass && !scene->bodies[j].mass) continue;
        if (scene->bodies[i].shape.etype != 3) continue;
        if (scene->bodies[j].shape.etype != 3) continue;
        ysr__compute_collision(scene, i, j);
    }
}

//
// Apply an impulse.
//
static inline void ysr__apply_impulse(ysr_scene* scene, int bid,
                                      const ym_vec3f impulse,
                                      const ym_vec3f local_pos) {
    if (!scene->bodies[bid].mass) return;
    scene->lin_vel[bid] += impulse * scene->bodies[bid].mass_inv;
    scene->ang_vel[bid] +=
        scene->bodies[bid].inertia_inv * ym_cross(local_pos, impulse);
}

//
// Shortcut math function.
//
static inline float ysr__muldot(const ym_vec3f& v, const ym_mat3f& m) {
    return ym_dot(v, m * v);
}

//
// Solve constraints with PGS.
//
YGL_API void ysr__solve_constraints(ysr_scene* scene, float dt) {
    // initialize computation
    for (int j = 0; j < scene->collisions.size(); j++) {
        scene->collisions[j].local_impulse = ym_zero3f;
        scene->collisions[j].impulse = ym_zero3f;
        ysr__collision* col = &scene->collisions[j];
        ym_vec3f r1 = col->frame.pos - scene->frame[col->bi].pos,
                 r2 = col->frame.pos - scene->frame[col->bj].pos;
        col->meff_inv =
            ym_vec3f{1 / (scene->bodies[col->bi].mass_inv +
                          scene->bodies[col->bj].mass_inv +
                          ysr__muldot(ym_cross(r1, col->frame.tangu),
                                      scene->bodies[col->bi].inertia_inv) +
                          ysr__muldot(ym_cross(r2, col->frame.tangu),
                                      scene->bodies[col->bj].inertia_inv)),
                     1 / (scene->bodies[col->bi].mass_inv +
                          scene->bodies[col->bj].mass_inv +
                          ysr__muldot(ym_cross(r1, col->frame.tangv),
                                      scene->bodies[col->bi].inertia_inv) +
                          ysr__muldot(ym_cross(r2, col->frame.tangv),
                                      scene->bodies[col->bj].inertia_inv)),
                     1 / (scene->bodies[col->bi].mass_inv +
                          scene->bodies[col->bj].mass_inv +
                          ysr__muldot(ym_cross(r1, col->frame.norm),
                                      scene->bodies[col->bi].inertia_inv) +
                          ysr__muldot(ym_cross(r2, col->frame.norm),
                                      scene->bodies[col->bj].inertia_inv))};
    }

    // compute relative velocity for visualization
    for (int j = 0; j < scene->collisions.size(); j++) {
        ysr__collision* col = &scene->collisions[j];
        ym_vec3f r1 = col->frame.pos - scene->frame[col->bi].pos,
                 r2 = col->frame.pos - scene->frame[col->bj].pos;
        ym_vec3f v1 = scene->lin_vel[col->bi] +
                      ym_cross(scene->ang_vel[col->bi], r1),
                 v2 = scene->lin_vel[col->bj] +
                      ym_cross(scene->ang_vel[col->bj], r2);
        col->vel_before = v2 - v1;
    }

    // solve constraints
    for (int i = 0; i < scene->iterations; i++) {
        for (int j = 0; j < scene->collisions.size(); j++) {
            ysr__collision* col = &scene->collisions[j];
            ym_vec3f r1 = col->frame.pos - scene->frame[col->bi].pos,
                     r2 = col->frame.pos - scene->frame[col->bj].pos;
            ym_vec3f v1 = scene->lin_vel[col->bi] +
                          ym_cross(scene->ang_vel[col->bi], r1),
                     v2 = scene->lin_vel[col->bj] +
                          ym_cross(scene->ang_vel[col->bj], r2);
            ym_vec3f vr = v2 - v1;
            ysr__apply_impulse(scene, col->bi, col->impulse, r1);
            ysr__apply_impulse(scene, col->bj, -col->impulse, r2);
            // float offset = col->depth*0.8f/dt;
            float offset = 0;
            ym_vec3f local_impulse =
                col->meff_inv * ym_vec3f{-ym_dot(col->frame.tangu, vr),
                                         -ym_dot(col->frame.tangv, vr),
                                         -ym_dot(col->frame.norm, vr) + offset};
            col->local_impulse += local_impulse;
            col->local_impulse.z =
                ym_clamp(col->local_impulse.z, 0.0f, ym_max_float);
            col->local_impulse.x =
                ym_clamp(col->local_impulse.x, -col->local_impulse.z * 0.6f,
                         col->local_impulse.z * 0.6f);
            col->local_impulse.y =
                ym_clamp(col->local_impulse.y, -col->local_impulse.z * 0.6f,
                         col->local_impulse.z - offset * 0.6f);
            col->impulse = col->local_impulse.z * col->frame.norm +
                           col->local_impulse.x * col->frame.tangu +
                           col->local_impulse.y * col->frame.tangv;
            ysr__apply_impulse(scene, col->bi, -col->impulse, r1);
            ysr__apply_impulse(scene, col->bj, col->impulse, r2);
        }
    }

    // compute relative velocity for visualization
    for (int j = 0; j < scene->collisions.size(); j++) {
        ysr__collision* col = &scene->collisions[j];
        ym_vec3f r1 = col->frame.pos - scene->frame[col->bi].pos,
                 r2 = col->frame.pos - scene->frame[col->bj].pos;
        ym_vec3f v1 = scene->lin_vel[col->bi] +
                      ym_cross(scene->ang_vel[col->bi], r1),
                 v2 = scene->lin_vel[col->bj] +
                      ym_cross(scene->ang_vel[col->bj], r2);
        col->vel_after = v2 - v1;
    }

    // recompute total impulse and velocity for visualization
    for (int j = 0; j < scene->collisions.size(); j++) {
        ysr__collision* col = &scene->collisions[j];
        col->impulse = col->local_impulse.z * col->frame.norm +
                       col->local_impulse.x * col->frame.tangu +
                       col->local_impulse.y * col->frame.tangv;
    }
}

//
// Advance simulation. Public API, see above.
//
YGL_API void ysr_advance(ysr_scene* scene, float dt) {
    // update centroid and inertia
    for (int bid = 0; bid < scene->bodies.size(); bid++) {
        ysr__body* body = &scene->bodies[bid];
        if (!body->mass) continue;
        body->inertia_inv = scene->frame[bid].rot * body->inertia_inv_local *
                            ym_transpose(scene->frame[bid].rot);
    }

    // compute collisions
    ysr__compute_collisions(scene);

    // apply external forces
    ym_vec3f gravity_impulse = scene->gravity * dt;
    for (int bid = 0; bid < scene->bodies.size(); bid++) {
        if (!scene->bodies[bid].mass) continue;
        scene->lin_vel[bid] += gravity_impulse;
    }

    // solve constraints
    ysr__solve_constraints(scene, dt);

    // apply drag
    for (int bid = 0; bid < scene->bodies.size(); bid++) {
        if (!scene->bodies[bid].mass) continue;
        scene->lin_vel[bid] *= 1 - scene->lin_drag;
        scene->ang_vel[bid] *= 1 - scene->ang_drag;
    }

    // update position
    for (int bid = 0; bid < scene->bodies.size(); bid++) {
        if (!scene->bodies[bid].mass) continue;

        // check for nans
        if (!ym_isfinite(scene->frame[bid].pos)) printf("nan detected\n");
        if (!ym_isfinite(scene->lin_vel[bid])) printf("nan detected\n");
        if (!ym_isfinite(scene->ang_vel[bid])) printf("nan detected\n");

        // update rotation
        scene->frame[bid].pos += scene->lin_vel[bid] * dt;
        float angle = ym_length(scene->ang_vel[bid]) * dt;
        if (angle) {
            ym_vec3f axis = ym_normalize(scene->ang_vel[bid]);
            scene->frame[bid].rot =
                ym_rotation_mat3(axis, angle) * scene->frame[bid].rot;
            // TODO: if using matrices, I gotta orthonormalize them
        }
    }

    // update acceleartion for collisions
    scene->overlap_refit(scene->overlap_ctx, scene->frame.data());
}

// -----------------------------------------------------------------------------
// C API IMPLEMENTATION
// -----------------------------------------------------------------------------

//
// Initialize a scene.
//
YGLC_API ysr_scene* ysrc_init_scene(int nbodies) {
    return ysr_init_scene(nbodies);
}

//
// Free a scene.
//
// Parameters:
// - scene: scene to clear
//
YGLC_API void ysrc_free_scene(ysr_scene* scene) {
    return ysr_free_scene(scene);
}

//
// Computes the moments of a shape.
//
YGLC_API void ysrc_compute_moments(int nelems, const int* elem, int etype,
                                   int nverts, const float* pos, float* volume,
                                   float center[3], float inertia[3]) {
    return ysr_compute_moments(nelems, elem, etype, nverts,
                               (const ym_vec3f*)pos, volume, (ym_vec3f*)center,
                               (ym_mat3f*)inertia);
}

//
// Set rigid body shape.
//
YGL_API void ysrc_set_body(ysr_scene* scene, int bid, const float xform[16],
                           float mass, const float inertia[9], int nelems,
                           const int* elem, int etype, int nverts,
                           const float* pos) {
    return ysr_set_body(scene, bid, ym_frame3f(ym_mat4f(xform)), mass,
                        ym_mat3f(inertia), nelems, elem, etype, nverts,
                        (const ym_vec3f*)pos);
}

//
// Get rigib body transform.
//
YGLC_API void ysrc_get_transform(const ysr_scene* scene, int bid,
                                 float xform[16]) {
    ym_frame3f frame = ysr_get_transform(scene, bid);
    *(ym_mat4f*)xform = ym_mat4f(frame);
}

//
// Set body transform.
//
YGLC_API void ysrc_set_transform(ysr_scene* scene, int bid,
                                 const float xform[16]) {
    return ysr_set_transform(scene, bid, ym_frame3f(ym_mat4f(xform)));
}

//
// Set collision callbacks.
//
YGLC_API void ysrc_set_collision(ysr_scene* scene, void* ctx,
                                 ysr_overlap_shapes overlaps,
                                 ysr_overlap_shape overlap,
                                 ysr_overlap_refit refit) {
    return ysr_set_collision(scene, ctx, overlaps, overlap, refit);
}

//
// Advance the simulation one step at a time.
//
YGLC_API void ysrc_advance(ysr_scene* scene, float dt) {
    return ysr_advance(scene, dt);
}

#endif

#endif
