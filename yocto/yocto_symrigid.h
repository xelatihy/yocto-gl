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
// All functions in this library are inlined by default for ease of use.
// To use the library as a .h/.c pair do the following:
// - to use as a .h, just #define YSR_NOINLINE before including this file
// - to use as a .c, just #define YSR_IMPLEMENTATION before including this file
//

//
// HISTORY:
// - v 0.0: initial release
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

#ifndef YSR_NOINLINE
#define YSR_API static inline
#else
#ifdef __cplusplus
#define YSR_API extern "C"
#else
#define YSR_API
#endif
#endif

#include <stdbool.h>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Shape element types
//
enum {
    yt_etype_point = 1,     // points
    yt_etype_line = 2,      // lines
    yt_etype_triangle = 3,  // triangle
    yt_etype_quad = 4,      // quads
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
typedef int (*ysr_overlap_shapes)(void* ctx, int** overlaps, int* ocapacity);

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
typedef bool (*ysr_overlap_shape)(void* ctx, int sid, const float pt[3],
                                  float max_dist, float* dist, int* eid,
                                  float* euv);

//
// Refit data structure after transform updates
//
// Parameters:
// - ctx: pointer to an object passed back to the function
// - xform: transform array
//
typedef void (*ysr_overlap_refit)(void* ctx, float* xform[16]);

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
YSR_API ysr_scene*
ysr_init_scene(int nbodies);

//
// Free a scene.
//
// Parameters:
// - scene: scene to clear
//
YSR_API void
ysr_free_scene(ysr_scene* scene);

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
YSR_API void
ysr_compute_moments(int nelems, const int* elem, int etype, int nverts,
                    float* pos, float* volume, float center[3],
                    float inertia[3]);

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
YSR_API void
ysr_set_body(ysr_scene* scene, int bid, const float xform[16], float mass,
             const float inertia[9], int nelems, const int* elem, int etype,
             int nverts, const float* pos);

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
YSR_API void
ysr_get_transform(const ysr_scene* scene, int bid, float xform[16]);

//
// Set body transform.
//
// Paramaters:
// - scene: rigib body scene
// - bid: body index
// - xform: transform
//
YSR_API void
ysr_set_transform(ysr_scene* scene, int bid, const float xform[16]);

//
// Set collision callbacks.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YSR_API void
ysr_set_collision(ysr_scene* scene, void* ctx, ysr_overlap_shapes overlaps,
                  ysr_overlap_shape overlap, ysr_overlap_refit refit);

//
// Advance the simulation one step at a time.
//
// Paramaters:
// - scene: rigib body scene
// - dt: time step
//
YSR_API void
ysr_advance(ysr_scene* scene, float dt);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YSR_NOINLINE) || defined(YSR_IMPLEMENTATION)

#include <assert.h>
#include <float.h>
#include <stdbool.h>

#include "yocto_math.h"

//
// check if finite
//
static inline bool
ysr__isfinite3f(const ym_vec3f v) {
    return isfinite(v.x) && isfinite(v.y) && isfinite(v.z);
}

//
// Rigid body properties
//
typedef struct ysr__body {
    // world space quantities
    float mass;                  // mass
    float mass_inv;              // mass inverse
    ym_mat3f inertia_inv;        // inverse of inertis tensor (world-space)
    ym_mat3f inertia_inv_local;  // inverse of inertis tensor (local-space)

    // rigid body shape
    struct {
        int nelems;       // number of elemets
        const int* elem;  // elements
        int etype;        // element type
        int nverts;       // number of vertices
        ym_vec3f* pos;    // vertex positions
    } shape;
} ysr__body;

//
// Collision point and response
//
typedef struct ysr__collision {
    int bi, bj;                       // bodies
    ym_vec3f pos;                     // position
    ym_vec3f norm, tangu, tangv;      // normal and tangets
    ym_vec3f impulse, local_impulse;  // impulses
    ym_vec3f vel_before, vel_after;   // velocities (for viz)
    ym_vec3f meff_inv;                // effective mass
    float depth;                      // penetration depth
} ysr__collision;

//
// Rigid body scene
//
struct ysr_scene {
    ym_vec3f* pos;      // position
    ym_mat3f* rot;      // rotation
    ym_vec3f* lin_vel;  // linear velocity
    ym_vec3f* ang_vel;  // angular velocity

    ym_mat4f* xform;  // computed transforms

    int nbodies;                        // number of bodies
    ysr__body* bodies;                  // bodies
    int ncollisions;                    // number of collisions
    int ccollisions;                    // capacity of collision array
    ysr__collision* collisions;         // collisions
    int nshapecollisions;               // number of shapes collision pairs
    int cshapecollisions;               // capacity oof shapes c.p.
    int* shapecollisions;               // shape collision pairs
    ym_vec3f gravity;                   // gravity
    float lin_drag, ang_drag;           // linear and angular drag
    int iterations;                     // solver iterations
    ysr_overlap_shapes overlap_shapes;  // overlap callbacks
    ysr_overlap_shape overlap_shape;    // overlap callbacks
    ysr_overlap_refit overlap_refit;    // overlap callbacks
    float** overlap_xforms;             // overlap callbacks
    void* overlap_ctx;                  // overlap callbacks
};

//
// Init scene. Public API, see above.
//
YSR_API ysr_scene*
ysr_init_scene(int nbodies) {
    ysr_scene* scene = (ysr_scene*)calloc(1, sizeof(ysr_scene));
    scene->nbodies = nbodies;
    scene->bodies = (ysr__body*)calloc(scene->nbodies, sizeof(ysr__body));
    scene->pos = (ym_vec3f*)calloc(scene->nbodies, sizeof(ym_vec3f));
    scene->rot = (ym_mat3f*)calloc(scene->nbodies, sizeof(ym_mat3f));
    scene->lin_vel = (ym_vec3f*)calloc(scene->nbodies, sizeof(ym_vec3f));
    scene->ang_vel = (ym_vec3f*)calloc(scene->nbodies, sizeof(ym_vec3f));
    scene->xform = (ym_mat4f*)calloc(scene->nbodies, sizeof(ym_mat4f));
    scene->overlap_xforms = (float**)calloc(scene->nbodies, sizeof(float*));
    for (int i = 0; i < nbodies; i++) {
        scene->overlap_xforms[i] = scene->xform[i].m;
    }
    scene->gravity = (ym_vec3f){ 0, -9.82, 0 };
    scene->lin_drag = 0.01;
    scene->ang_drag = 0.01;
    scene->iterations = 20;
    return scene;
}

//
// Free scene. Public API, see above.
//
YSR_API void
ysr_free_scene(ysr_scene* scene) {
    free(scene->bodies);
    free(scene->pos);
    free(scene->rot);
    free(scene->lin_vel);
    free(scene->ang_vel);
    free(scene->xform);
    free(scene->overlap_xforms);
    free(scene);
}

//
// Compute moments. Public API, see above.
//
// Implementation notes:
// - follow eberly
// http://www.geometrictools.com/Documentation/PolyhedralMassProperties.pdf
// - re-implement this using
// https://github.com/melax/sandbox/blob/master/include/geometric.h
//
YSR_API void
ysr_compute_moments(int nelems, const int* elem, int etype, int nverts,
                    float* pos_, float* volume, float center[3],
                    float inertia[3]) {
    ym_vec3f* pos = (ym_vec3f*)pos_;

    // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx
    float integral[10] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for (int i = 0; i < nelems; i++) {
        const int* f = elem + i * 3;
        ym_vec3f v0 = pos[f[0]], v1 = pos[f[1]], v2 = pos[f[2]];

        // Get cross product of edges and normal vector.
        ym_vec3f V1mV0 = ym_sub3f(v1, v0);
        ym_vec3f V2mV0 = ym_sub3f(v2, v0);
        ym_vec3f N = ym_cross3f(V1mV0, V2mV0);

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
    center[0] = integral[1] / *volume;
    center[1] = integral[2] / *volume;
    center[2] = integral[3] / *volume;

    // inertia relative to center
    inertia[0] = integral[5] + integral[6] -
                 (*volume) * (center[1] * center[1] + center[2] * center[2]);
    inertia[1] = -integral[7] + (*volume) * center[0] * center[1];
    inertia[2] = -integral[9] + (*volume) * center[2] * center[0];
    inertia[3] = inertia[1];
    inertia[4] = integral[4] + integral[6] -
                 (*volume) * (center[2] * center[2] + center[0] * center[0]);
    inertia[5] = -integral[8] + (*volume) * center[1] * center[2];
    inertia[6] = inertia[2];
    inertia[7] = inertia[4];
    inertia[8] = integral[4] + integral[5] -
                 (*volume) * (center[0] * center[0] + center[1] * center[1]);
}

//
// Set rigid body. Public API, see above.
//
YSR_API void
ysr_set_body(ysr_scene* scene, int bid, const float xform[16], float mass,
             const float inertia[9], int nelems, const int* elem, int etype,
             int nverts, const float* pos) {
    ysr__body* body = scene->bodies + bid;
    body->mass = mass;
    body->mass_inv = (mass) ? 1 / mass : 0;
    body->inertia_inv_local = (mass)
                                  ? ym_inverse3f(*(ym_mat3f*)inertia)
                                  : (ym_mat3f){ { 0, 0, 0, 0, 0, 0, 0, 0, 0 } };
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
YSR_API void
ysr_set_collision(ysr_scene* scene, void* ctx, ysr_overlap_shapes overlaps,
                  ysr_overlap_shape overlap, ysr_overlap_refit refit) {
    scene->overlap_shapes = overlaps;
    scene->overlap_shape = overlap;
    scene->overlap_refit = refit;
    scene->overlap_ctx = ctx;
}

//
// Compute xform for a body.
//
static inline void
ysr__compute_xforms(ysr_scene* scene, int bid) {
    ym_mat3f m = scene->rot[bid];
    ym_vec3f t = scene->pos[bid];
    scene->xform[bid] =
        (ym_mat4f){ { m.m[0], m.m[1], m.m[2], t.x, m.m[3], m.m[4], m.m[5], t.y,
                      m.m[6], m.m[7], m.m[8], t.z, 0, 0, 0, 1 } };
}

//
// Get transform. Public API, see above.
//
YSR_API void
ysr_get_transform(const ysr_scene* scene, int bid, float xform[16]) {
    *(ym_mat4f*)xform = scene->xform[bid];
}

//
// Set transform. Public API, see above.
//
YSR_API void
ysr_set_transform(ysr_scene* scene, int bid, const float xform[16]) {
    scene->pos[bid] = (ym_vec3f){ xform[3], xform[7], xform[11] };
    scene->rot[bid] =
        (ym_mat3f){ { xform[0], xform[1], xform[2], xform[4], xform[5],
                      xform[6], xform[8], xform[9], xform[10] } };
    ysr__compute_xforms(scene, bid);
}

//
// Compute collisions.
//
static inline void
ysr__compute_collision(ysr_scene* scene, int bi, int bj) {
    ysr__body* b1 = scene->bodies + bi;
    ysr__body* b2 = scene->bodies + bj;
    for (int j = 0; j < b2->shape.nverts; j++) {
        ym_vec3f p = ym_transform_point3f(scene->xform[bj], b2->shape.pos[j]);
        int eid;
        float dist, euv[2];
        bool hit = scene->overlap_shape(scene->overlap_ctx, bi, &p.x, FLT_MAX,
                                        &dist, &eid, euv);
        if (!hit) continue;
        const int* f = b1->shape.elem + eid * 3;
        ym_vec3f v0 = b1->shape.pos[f[0]], v1 = b1->shape.pos[f[1]],
                 v2 = b1->shape.pos[f[2]];
        ym_vec3f tp = ym_transform_point3f(
            scene->xform[bi], ym_blerp3f(v0, v1, v2, euv[0], euv[1]));
        ym_vec3f n = ym_transform_direction3f(
            scene->xform[bi], ym_cross3f(ym_sub3f(v1, v0), ym_sub3f(v2, v0)));
        const float eps = -0.01f;
        ym_vec3f ptp = ym_normalize3f(ym_sub3f(p, tp));
        if (ym_dot3f(n, ptp) > eps) continue;
        if (scene->ncollisions + 1 > scene->ccollisions) {
            scene->ccollisions =
                (scene->ccollisions) ? scene->ccollisions * 2 : 32;
            scene->collisions = (ysr__collision*)realloc(
                scene->collisions, sizeof(ysr__collision) * scene->ccollisions);
        }
        ysr__collision* col = scene->collisions + scene->ncollisions;
        memset(col, 0, sizeof(*col));
        scene->ncollisions += 1;
        col->bi = bi;
        col->bj = bj;
        col->depth = dist;
        col->pos = p;
        col->norm = n;
        col->tangu = ym_orthogonal3f(n);
        col->tangv = ym_cross3f(n, col->tangv);
    }
}

//
// Compute collisions.
//
static inline void
ysr__compute_collisions(ysr_scene* scene) {
    // check which shapes might overlap
    scene->nshapecollisions = scene->overlap_shapes(
        scene->overlap_ctx, &scene->shapecollisions, &scene->cshapecollisions);
    // test all pair-wise objects
    scene->ncollisions = 0;
    for (int c = 0; c < scene->nshapecollisions; c++) {
        int i = scene->shapecollisions[c * 2 + 0],
            j = scene->shapecollisions[c * 2 + 1];
        if (!scene->bodies[i].mass && !scene->bodies[j].mass) continue;
        if (scene->bodies[i].shape.etype != 3) continue;
        if (scene->bodies[j].shape.etype != 3) continue;
        ysr__compute_collision(scene, i, j);
    }
}

//
// Apply an impulse.
//
static inline void
ysr__apply_impulse(ysr_scene* scene, int bid, const ym_vec3f impulse,
                   const ym_vec3f local_pos) {
    if (!scene->bodies[bid].mass) return;
    scene->lin_vel[bid] = ym_sum3f(
        scene->lin_vel[bid], ym_smul3f(impulse, scene->bodies[bid].mass_inv));
    scene->ang_vel[bid] = ym_sum3f(scene->ang_vel[bid],
                                   ym_vmul3f(scene->bodies[bid].inertia_inv,
                                             ym_cross3f(local_pos, impulse)));
}

//
// Shortcut math function.
//
static inline float
ysr__muldot(const ym_vec3f v, const ym_mat3f m) {
    return ym_dot3f(v, ym_vmul3f(m, v));
}

//
// Solve constraints with PGS.
//
YS_API void
ysr__solve_constraints(ysr_scene* scene, float dt) {
    // initialize computation
    for (int j = 0; j < scene->ncollisions; j++) {
        scene->collisions[j].local_impulse = ym_zero3f();
        scene->collisions[j].impulse = ym_zero3f();
        ysr__collision* col = scene->collisions + j;
        ym_vec3f r1 = ym_sub3f(col->pos, scene->pos[col->bi]),
                 r2 = ym_sub3f(col->pos, scene->pos[col->bj]);
        col->meff_inv =
            (ym_vec3f){ 1 / (scene->bodies[col->bi].mass_inv +
                             scene->bodies[col->bj].mass_inv +
                             ysr__muldot(ym_cross3f(r1, col->tangu),
                                         scene->bodies[col->bi].inertia_inv) +
                             ysr__muldot(ym_cross3f(r2, col->tangu),
                                         scene->bodies[col->bj].inertia_inv)),
                        1 / (scene->bodies[col->bi].mass_inv +
                             scene->bodies[col->bj].mass_inv +
                             ysr__muldot(ym_cross3f(r1, col->tangv),
                                         scene->bodies[col->bi].inertia_inv) +
                             ysr__muldot(ym_cross3f(r2, col->tangv),
                                         scene->bodies[col->bj].inertia_inv)),
                        1 / (scene->bodies[col->bi].mass_inv +
                             scene->bodies[col->bj].mass_inv +
                             ysr__muldot(ym_cross3f(r1, col->norm),
                                         scene->bodies[col->bi].inertia_inv) +
                             ysr__muldot(ym_cross3f(r2, col->norm),
                                         scene->bodies[col->bj].inertia_inv)) };
    }

    // compute relative velocity for visualization
    for (int j = 0; j < scene->ncollisions; j++) {
        ysr__collision* col = scene->collisions + j;
        ym_vec3f r1 = ym_sub3f(col->pos, scene->pos[col->bi]),
                 r2 = ym_sub3f(col->pos, scene->pos[col->bj]);
        ym_vec3f v1 = ym_sum3f(scene->lin_vel[col->bi],
                               ym_cross3f(scene->ang_vel[col->bi], r1)),
                 v2 = ym_sum3f(scene->lin_vel[col->bj],
                               ym_cross3f(scene->ang_vel[col->bj], r2));
        col->vel_before = ym_sub3f(v2, v1);
    }

    // solve constraints
    for (int i = 0; i < scene->iterations; i++) {
        for (int j = 0; j < scene->ncollisions; j++) {
            ysr__collision* col = scene->collisions + j;
            ym_vec3f r1 = ym_sub3f(col->pos, scene->pos[col->bi]),
                     r2 = ym_sub3f(col->pos, scene->pos[col->bj]);
            ym_vec3f v1 = ym_sum3f(scene->lin_vel[col->bi],
                                   ym_cross3f(scene->ang_vel[col->bi], r1)),
                     v2 = ym_sum3f(scene->lin_vel[col->bj],
                                   ym_cross3f(scene->ang_vel[col->bj], r2));
            ym_vec3f vr = ym_sub3f(v2, v1);
            ysr__apply_impulse(scene, col->bi, col->impulse, r1);
            ysr__apply_impulse(scene, col->bj, ym_neg3f(col->impulse), r2);
            // float offset = col->depth*0.8f/dt;
            float offset = 0;
            ym_vec3f local_impulse = ym_mul3f(
                col->meff_inv, (ym_vec3f){ -ym_dot3f(col->tangu, vr),
                                           -ym_dot3f(col->tangv, vr),
                                           -ym_dot3f(col->norm, vr) + offset });
            col->local_impulse = ym_sum3f(col->local_impulse, local_impulse);
            col->local_impulse.z = ym_fclampf(col->local_impulse.z, 0, FLT_MAX);
            col->local_impulse.x =
                ym_fclampf(col->local_impulse.x, -col->local_impulse.z * 0.6,
                           col->local_impulse.z * 0.6);
            col->local_impulse.y =
                ym_fclampf(col->local_impulse.y, -col->local_impulse.z * 0.6,
                           col->local_impulse.z - offset * 0.6);
            col->impulse = ym_zero3f();
            col->impulse =
                ym_ssum3f(col->impulse, col->local_impulse.z, col->norm);
            col->impulse =
                ym_ssum3f(col->impulse, col->local_impulse.x, col->tangu);
            col->impulse =
                ym_ssum3f(col->impulse, col->local_impulse.y, col->tangv);
            ysr__apply_impulse(scene, col->bi, ym_neg3f(col->impulse), r1);
            ysr__apply_impulse(scene, col->bj, col->impulse, r2);
        }
    }

    // compute relative velocity for visualization
    for (int j = 0; j < scene->ncollisions; j++) {
        ysr__collision* col = scene->collisions + j;
        ym_vec3f r1 = ym_sub3f(col->pos, scene->pos[col->bi]),
                 r2 = ym_sub3f(col->pos, scene->pos[col->bj]);
        ym_vec3f v1 = ym_sum3f(scene->lin_vel[col->bi],
                               ym_cross3f(scene->ang_vel[col->bi], r1)),
                 v2 = ym_sum3f(scene->lin_vel[col->bj],
                               ym_cross3f(scene->ang_vel[col->bj], r2));
        col->vel_after = ym_sub3f(v2, v1);
    }

    // recompute total impulse and velocity for visualization
    for (int j = 0; j < scene->ncollisions; j++) {
        ysr__collision* col = scene->collisions + j;
        col->impulse = ym_zero3f();
        col->impulse = ym_ssum3f(col->impulse, col->local_impulse.z, col->norm);
        col->impulse =
            ym_ssum3f(col->impulse, col->local_impulse.x, col->tangu);
        col->impulse =
            ym_ssum3f(col->impulse, col->local_impulse.y, col->tangv);
    }
}

//
// Advance simulation. Public API, see above.
//
YSR_API void
ysr_advance(ysr_scene* scene, float dt) {
    // update centroid and inertia
    for (int bid = 0; bid < scene->nbodies; bid++) {
        ysr__body* body = scene->bodies + bid;
        if (!body->mass) continue;
        body->inertia_inv = ym_mmul3f(
            scene->rot[bid], ym_mmul3f(body->inertia_inv_local,
                                       ym_transpose3f(scene->rot[bid])));
    }

    // compute collisions
    ysr__compute_collisions(scene);

    // apply external forces
    ym_vec3f gravity_impulse = ym_smul3f(scene->gravity, dt);
    for (int bid = 0; bid < scene->nbodies; bid++) {
        if (!scene->bodies[bid].mass) continue;
        scene->lin_vel[bid] = ym_sum3f(scene->lin_vel[bid], gravity_impulse);
    }

    // solve constraints
    ysr__solve_constraints(scene, dt);

    // apply drag
    for (int bid = 0; bid < scene->nbodies; bid++) {
        if (!scene->bodies[bid].mass) continue;
        scene->lin_vel[bid] =
            ym_smul3f(scene->lin_vel[bid], 1 - scene->lin_drag);
        scene->ang_vel[bid] =
            ym_smul3f(scene->ang_vel[bid], 1 - scene->ang_drag);
    }

    // update position
    for (int bid = 0; bid < scene->nbodies; bid++) {
        if (!scene->bodies[bid].mass) continue;

        // check for nans
        if (!ysr__isfinite3f(scene->pos[bid])) printf("nan detected\n");
        if (!ysr__isfinite3f(scene->lin_vel[bid])) printf("nan detected\n");
        if (!ysr__isfinite3f(scene->ang_vel[bid])) printf("nan detected\n");

        // update rotation
        scene->pos[bid] =
            ym_sum3f(scene->pos[bid], ym_smul3f(scene->lin_vel[bid], dt));
        float angle = ym_length3f(scene->ang_vel[bid]) * dt;
        if (angle) {
            ym_vec3f axis = ym_normalize3f(scene->ang_vel[bid]);
            scene->rot[bid] =
                ym_mmul3f(ym_rotation3f(angle, axis), scene->rot[bid]);
            // TODO: if using matrices, I gotta orthonormalize them
        }

        // update transforms
        ysr__compute_xforms(scene, bid);
    }

    // update acceleartion for collisions
    scene->overlap_refit(scene->overlap_ctx, scene->overlap_xforms);
}

#endif

#endif
