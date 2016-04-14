//
// YOCTO_MATH: a collection of vector math functions for 2, 3, 4 float vectors
// and 4x4 matrices. Used by YOCTO authors to create demos. This is not
// documented at the moment.
//
// **We strongly remommend to use a more complete math library.**
//
// We developed our own library since we felt that all existing ones are either
// complete, but unreadable of with tons of dependencies, or just as incomplete
// and untested as ours.
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

#ifndef _YM_H_
#define _YM_H_

#include <math.h>
#include <stdbool.h>

// pi
#define ym_pif 3.14159265f
#define ym_pi 3.1415926535897932384626433832795

// a few utility functions
static inline int
ym_min(int x, int y) {
    return (x < y) ? x : y;
}
static inline int
ym_max(int x, int y) {
    return (x > y) ? x : y;
}
static inline int
ym_clamp(int x, int m, int M) {
    return ym_min(M, ym_max(m, x));
}
static inline double
ym_fclamp(double x, double m, double M) {
    return fmin(M, fmax(m, x));
}
static inline float
ym_fclampf(float x, float m, float M) {
    return fminf(M, fmaxf(m, x));
}

typedef struct { float x, y; } ym_vec2f;
typedef struct { float x, y, z; } ym_vec3f;
typedef struct { float x, y, z, w; } ym_vec4f;
typedef struct { float m[16]; } ym_mat4f;

static inline ym_vec3f
ym_zero3f() {
    return (ym_vec3f){ 0, 0, 0 };
}

static inline ym_vec3f
ym_vmake3f(const float* v) {
    return (ym_vec3f){ v[0], v[1], v[2] };
}

static inline ym_vec3f
ym_neg3f(const ym_vec3f a) {
    return (ym_vec3f){ -a.x, -a.y, -a.z };
}

static inline ym_vec3f
ym_sum3f(const ym_vec3f a, const ym_vec3f b) {
    return (ym_vec3f){ a.x + b.x, a.y + b.y, a.z + b.z };
}

static inline ym_vec3f
ym_sub3f(const ym_vec3f a, const ym_vec3f b) {
    return (ym_vec3f){ a.x - b.x, a.y - b.y, a.z - b.z };
}

static inline ym_vec3f
ym_mul3f(const ym_vec3f a, const ym_vec3f b) {
    return (ym_vec3f){ a.x * b.x, a.y * b.y, a.z * b.z };
}

static inline ym_vec3f
ym_smul3f(const ym_vec3f a, float b) {
    return (ym_vec3f){ a.x * b, a.y * b, a.z * b };
}

static inline float
ym_dot3f(const ym_vec3f a, const ym_vec3f b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

static inline float
ym_length3f(const ym_vec3f a) {
    return (float)sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}

static inline ym_vec3f
ym_normalize3f(const ym_vec3f a) {
    return ym_smul3f(a, (1 / ym_length3f(a)));
}

static inline float
ym_dist3f(const ym_vec3f a, const ym_vec3f b) {
    return (float)sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) +
                       (a.z - b.z) * (a.z - b.z));
}

static inline ym_vec3f
ym_cross3f(const ym_vec3f a, const ym_vec3f b) {
    return (ym_vec3f){ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                       a.x * b.y - a.y * b.x };
}

// orthogonal (from:
// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
static inline ym_vec3f
ym_orthogonal3f(const ym_vec3f v) {
    return fabs(v.x) > fabs(v.z) ? (ym_vec3f){ -v.y, v.x, 0 }
                                 : (ym_vec3f){ 0, -v.z, v.y };
}

static inline ym_vec3f
ym_orthonormalize3f(const ym_vec3f a, const ym_vec3f b) {
    return ym_normalize3f(ym_sub3f(a, ym_smul3f(b, ym_dot3f(a, b))));
}

static inline ym_vec2f
ym_zero2f() {
    return (ym_vec2f){ 0, 0 };
}

static inline ym_vec2f
ym_neg2f(const ym_vec2f a) {
    return (ym_vec2f){ -a.x, -a.y };
}

static inline ym_vec2f
ym_sum2f(const ym_vec2f a, const ym_vec2f b) {
    return (ym_vec2f){ a.x + b.x, a.y + b.y };
}

static inline ym_vec2f
ym_sub2f(const ym_vec2f a, const ym_vec2f b) {
    return (ym_vec2f){ a.x - b.x, a.y - b.y };
}

static inline ym_vec2f
ym_mul2f(const ym_vec2f a, const ym_vec2f b) {
    return (ym_vec2f){ a.x * b.x, a.y * b.y };
}

static inline ym_vec2f
ym_smul2f(const ym_vec2f a, float b) {
    return (ym_vec2f){ a.x * b, a.y * b };
}

static inline float
ym_dot2f(const ym_vec2f a, const ym_vec2f b) {
    return (a.x * b.x + a.y * b.y);
}

static inline ym_vec4f
ym_zero4f() {
    return (ym_vec4f){ 0, 0, 0, 0 };
}

static inline ym_vec4f
ym_neg4f(const ym_vec4f a) {
    return (ym_vec4f){ -a.x, -a.y, -a.z, -a.w };
}

static inline ym_vec4f
ym_sum4f(const ym_vec4f a, const ym_vec4f b) {
    return (ym_vec4f){ a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

static inline ym_vec4f
ym_sub4f(const ym_vec4f a, const ym_vec4f b) {
    return (ym_vec4f){ a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
}

static inline ym_vec4f
ym_mul4f(const ym_vec4f a, const ym_vec4f b) {
    return (ym_vec4f){ a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w };
}

static inline ym_vec4f
ym_smul4f(const ym_vec4f a, float b) {
    return (ym_vec4f){ a.x * b, a.y * b, a.z * b, a.w * b };
}

static inline float
ym_dot4f(const ym_vec4f a, const ym_vec4f b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w);
}

static inline ym_mat4f
ym_identity4f() {
    return (ym_mat4f){ { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 } };
}

static inline ym_vec4f
ym_vmul4f(const ym_mat4f a, const ym_vec4f b) {
    const ym_vec4f* mv = (const ym_vec4f*)a.m;
    return (ym_vec4f){
        ym_dot4f(mv[0], b), ym_dot4f(mv[1], b), ym_dot4f(mv[2], b),
        ym_dot4f(mv[3], b),
    };
}

static inline ym_mat4f
ym_mmul4f(const ym_mat4f a, const ym_mat4f b) {
    ym_mat4f c = { { 0 } };
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            c.m[i * 4 + j] = 0;
            for (int k = 0; k < 4; k++) {
                c.m[i * 4 + j] += a.m[i * 4 + k] * b.m[k * 4 + j];
            }
        }
    }
    return c;
}

// http://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform
static inline ym_mat4f
ym_inverse4f(const ym_mat4f m) {
    float a[4][4];
    *(ym_mat4f*)a = m;

    float s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    float s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
    float s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3];
    float s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
    float s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3];
    float s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];

    float c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3];
    float c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
    float c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2];
    float c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
    float c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2];
    float c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];

    // TODO: Should check for 0 determinant
    float invdet =
        1.0f / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

    float b[4][4];

    b[0][0] = (a[1][1] * c5 - a[1][2] * c4 + a[1][3] * c3) * invdet;
    b[0][1] = (-a[0][1] * c5 + a[0][2] * c4 - a[0][3] * c3) * invdet;
    b[0][2] = (a[3][1] * s5 - a[3][2] * s4 + a[3][3] * s3) * invdet;
    b[0][3] = (-a[2][1] * s5 + a[2][2] * s4 - a[2][3] * s3) * invdet;

    b[1][0] = (-a[1][0] * c5 + a[1][2] * c2 - a[1][3] * c1) * invdet;
    b[1][1] = (a[0][0] * c5 - a[0][2] * c2 + a[0][3] * c1) * invdet;
    b[1][2] = (-a[3][0] * s5 + a[3][2] * s2 - a[3][3] * s1) * invdet;
    b[1][3] = (a[2][0] * s5 - a[2][2] * s2 + a[2][3] * s1) * invdet;

    b[2][0] = (a[1][0] * c4 - a[1][1] * c2 + a[1][3] * c0) * invdet;
    b[2][1] = (-a[0][0] * c4 + a[0][1] * c2 - a[0][3] * c0) * invdet;
    b[2][2] = (a[3][0] * s4 - a[3][1] * s2 + a[3][3] * s0) * invdet;
    b[2][3] = (-a[2][0] * s4 + a[2][1] * s2 - a[2][3] * s0) * invdet;

    b[3][0] = (-a[1][0] * c3 + a[1][1] * c1 - a[1][2] * c0) * invdet;
    b[3][1] = (a[0][0] * c3 - a[0][1] * c1 + a[0][2] * c0) * invdet;
    b[3][2] = (-a[3][0] * s3 + a[3][1] * s1 - a[3][2] * s0) * invdet;
    b[3][3] = (a[2][0] * s3 - a[2][1] * s1 + a[2][2] * s0) * invdet;

    return *(ym_mat4f*)b;
}

static inline ym_vec3f
ym_transform_point3f(const ym_mat4f a, const ym_vec3f b) {
    ym_vec4f vb = { b.x, b.y, b.z, 1 };
    ym_vec4f tvb = ym_vmul4f(a, vb);
    return (ym_vec3f){ tvb.x / tvb.w, tvb.y / tvb.w, tvb.z / tvb.w };
}

static inline ym_vec3f
ym_transform_vector3f(const ym_mat4f a, const ym_vec3f b) {
    ym_vec4f vb = { b.x, b.y, b.z, 0 };
    ym_vec4f tvb = ym_vmul4f(a, vb);
    return (ym_vec3f){ tvb.x, tvb.y, tvb.z };
}

static inline ym_vec3f
ym_transform_direction3f(const ym_mat4f a, const ym_vec3f b) {
    return ym_normalize3f(ym_transform_vector3f(a, b));
}

static inline ym_vec3f
ym_affine_transform_point3f(const ym_mat4f a, const ym_vec3f b) {
    return (ym_vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z + a.m[3],
                       a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z + a.m[7],
                       a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z + a.m[11] };
}

static inline ym_vec3f
ym_affine_transform_vector3f(const ym_mat4f a, const ym_vec3f b) {
    return (ym_vec3f){ a.m[0] * b.x + a.m[1] * b.y + a.m[2] * b.z,
                       a.m[4] * b.x + a.m[5] * b.y + a.m[6] * b.z,
                       a.m[8] * b.x + a.m[9] * b.y + a.m[10] * b.z };
}

static inline ym_vec3f
ym_affine_transform_direction3f(const ym_mat4f a, const ym_vec3f b) {
    return ym_normalize3f(ym_affine_transform_vector3f(a, b));
}

static inline ym_mat4f
ym_translation4f(const ym_vec3f a) {
    return (
        ym_mat4f){ { 1, 0, 0, a.x, 0, 1, 0, a.y, 0, 0, 1, a.z, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_scaling4f(const ym_vec3f a) {
    return (
        ym_mat4f){ { a.x, 0, 0, 0, 0, a.y, 0, 0, 0, 0, a.z, 0, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_rotation4f(float a, const ym_vec3f b) {
    float s = sinf(a), c = cosf(a);
    ym_vec3f vv = ym_normalize3f(b);
    return (ym_mat4f){
        { c + (1 - c) * vv.x * vv.x, (1 - c) * vv.x * vv.y - s * vv.z,
          (1 - c) * vv.x * vv.z + s * vv.y, 0, (1 - c) * vv.x * vv.y + s * vv.z,
          c + (1 - c) * vv.y * vv.y, (1 - c) * vv.y * vv.z - s * vv.x, 0,
          (1 - c) * vv.x * vv.z - s * vv.y, (1 - c) * vv.y * vv.z + s * vv.x,
          c + (1 - c) * vv.z * vv.z, 0, 0, 0, 0, 1 }
    };
}

static inline ym_mat4f
ym_rotationx4f(float a) {
    float s = sinf(a), c = cosf(a);
    return (ym_mat4f){ { 1, 0, 0, 0, 0, c, -s, 0, 0, s, c, 0, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_rotationy4f(float a) {
    float s = sinf(a), c = cosf(a);
    return (ym_mat4f){ { c, 0, s, 0, 0, 1, 0, 0, -s, 0, c, 0, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_rotationz4f(float a) {
    float s = sinf(a), c = cosf(a);
    return (ym_mat4f){ { c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_frustum4f(float l, float r, float b, float t, float n, float f) {
    return (ym_mat4f){ { 2 * n / (r - l), 0, (r + l) / (r - l), 0, 0,
                         2 * n / (t - b), (t + b) / (t - b), 0, 0, 0,
                         -(f + n) / (f - n), -2 * f * n / (f - n), 0, 0, -1,
                         0 } };
}

static inline ym_mat4f
ym_ortho4f(float l, float r, float b, float t, float n, float f) {
    return (ym_mat4f){ { 2 / (r - l), 0, 0, -(r + l) / (r - l), 0, 2 / (t - b),
                         0, -(t + b) / (t - b), 0, 0, -2 / (f - n),
                         -(f + n) / (f - n), 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_lookat_view4f(const ym_vec3f eye, const ym_vec3f center, const ym_vec3f up) {
    ym_vec3f w = ym_normalize3f(ym_sub3f(eye, center));
    ym_vec3f u = ym_normalize3f(ym_cross3f(up, w));
    ym_vec3f v = ym_normalize3f(ym_cross3f(w, u));
    return (ym_mat4f){ { u.x, u.y, u.z, -ym_dot3f(u, eye), v.x, v.y, v.z,
                         -ym_dot3f(v, eye), w.x, w.y, w.z, -ym_dot3f(w, eye), 0,
                         0, 0, 1 } };
}

static inline ym_mat4f
ym_lookat_xform4f(const ym_vec3f eye, const ym_vec3f center,
                  const ym_vec3f up) {
    ym_vec3f w = ym_normalize3f(ym_sub3f(eye, center));
    ym_vec3f u = ym_normalize3f(ym_cross3f(up, w));
    ym_vec3f v = ym_normalize3f(ym_cross3f(w, u));
    return (ym_mat4f){ { u.x, v.x, w.x, eye.x, u.y, v.y, w.y, eye.y, u.z, v.z,
                         w.z, eye.z, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_ortho2d4f(float l, float r, float b, float t) {
    return (ym_mat4f){ { 2 / (r - l), 0, 0, -(r + l) / (r - l), 0, 2 / (t - b),
                         0, -(t + b) / (t - b), 0, 0, -1, 0, 0, 0, 0, 1 } };
}

static inline ym_mat4f
ym_perspective4f(float fovy, float aspect, float near, float far) {
    float f = 1 / tanf(fovy / 2);
    return (ym_mat4f){ { f / aspect, 0, 0, 0, 0, f, 0, 0, 0, 0,
                         (far + near) / (near - far),
                         2 * far * near / (near - far), 0, 0, -1, 0 } };
}

typedef struct { ym_vec3f min, max; } ym_range3f;

static inline ym_range3f
ym_rinit3f() {
    return (ym_range3f){ { HUGE_VALF, HUGE_VALF, HUGE_VALF },
                         { -HUGE_VALF, -HUGE_VALF, -HUGE_VALF } };
}

static inline ym_range3f
ym_rexpand3f(const ym_range3f a, const ym_vec3f b) {
    return (ym_range3f){
        { fminf(a.min.x, b.x), fminf(a.min.y, b.y), fminf(a.min.z, b.z) },
        { fmaxf(a.max.x, b.x), fmaxf(a.max.y, b.y), fmaxf(a.max.z, b.z) }
    };
}

static inline ym_range3f
ym_runion3f(const ym_range3f a, const ym_range3f b) {
    return (ym_range3f){ { fminf(a.min.x, b.min.x), fminf(a.min.y, b.min.y),
                           fminf(a.min.z, b.min.z) },
                         { fmaxf(a.max.x, b.max.x), fmaxf(a.max.y, b.max.y),
                           fmaxf(a.max.z, b.max.z) } };
}

static inline ym_vec3f
ym_rcenter3f(const ym_range3f a) {
    return ym_smul3f(ym_sum3f(a.max, a.min), 0.5f);
}

static inline ym_vec3f
ym_rsize3f(const ym_range3f a) {
    return ym_sub3f(a.max, a.min);
}

#endif
