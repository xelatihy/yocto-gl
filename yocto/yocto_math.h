//
// YOCTO_MATH: a collection of vector math functions and simple containers
// used to implement YOCTO. Features include
// - static length float vectors, with specialization for 2, 3, 4 length
// - static length matrices, with specialization for 2x2, 3x3, 4x4
// - affine and rigid transforms
// - linear algebra operations and transforms for fixed length matrices/vecs
// - axis aligned bounding boxes
// - rays
// - random number generation via PCG32
// - a few hash functions
// - vector/string container replacement
// - timer (depends on C++11 chrono)
//
// The containers are only meant to avoid using the STL in YOCTO. These
// containers only support a subset of the members of the STL ones since they
// are meant to be simple replacement rather than complete "Reinvent-The-Wheel"
// efforts.
//
// While we tested this library in the implementation of our other ones, we
// consider this code incomplete and remommend to use a more complete math
// library. So use it at your own peril.
//
// We developed our own library since we felt that all existing ones are either
// complete, but unreadable of with tons of dependencies, or just as incomplete
// and untested as ours.
//

//
// COMPILATION:
//
// This library can only be used as a header only library in C++ since it uses
// templates for its basic types. To use STL containers in YOCTO, instead of
// the built in ones, #define YGL_USESTL before including this file or any
// other YOCTO file that dependes on this.
//

//
// HISTORY:
// - v 0.1: C++ only implementation
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

//
//  LICENSE of included software
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//

#ifndef _YM_H_
#define _YM_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <type_traits>

// -----------------------------------------------------------------------------
// CONSTANTS
// -----------------------------------------------------------------------------

//
// pi
//
const float ym_pif = 3.14159265f;
const double ym_pi = 3.1415926535897932384626433832795;

//
// shortcut for numeric limits
//
const float ym_max_float = std::numeric_limits<float>::max();

// -----------------------------------------------------------------------------
// BASIC MATH FUNCTIONS
// -----------------------------------------------------------------------------

//
// min/max
//
template <typename T>
inline T ym_min(const T& x, const T& y) {
    return (x < y) ? x : y;
}
template <typename T>
inline T ym_max(const T& x, const T& y) {
    return (x > y) ? x : y;
}

//
// swap utility
//
template <typename T>
inline void ym_swap(T& a, T& b) {
    T c = a;
    a = b;
    b = c;
}

//
// Clamp a value between a min and max values
//
template <typename T>
inline T ym_clamp(const T& x, const T& min, const T& max) {
    return ym_min(ym_max(x, min), max);
}
template <typename T, typename T1>
inline T ym_clamp(const T& x, const T1& min, const T1& max) {
    return ym_min(ym_max(x, (T)min), (T)max);
}

//
// Check if a value if finite (will then be specialized on vectors)
//
template <typename T>
inline bool ym_isfinite(const T& v) {
    return std::isfinite(v);
}

//
// Linear interpolation and triangle baricentric interpolation.
//
template <typename T, typename T1>
inline T ym_lerp(const T& a, const T& b, T1 t) {
    return a * (1 - t) + b * t;
}
template <typename T, typename T1>
inline T ym_blerp(const T& a, const T& b, const T& c, T1 u, T1 v) {
    return a * (1 - u - v) + b * u + c * v;
}

// -----------------------------------------------------------------------------
// LINEAR ALGEBRA TYPES (VECTORS and MATRICES)
// -----------------------------------------------------------------------------

//
// Vector of element of compile time length with default initializer,
// constant initialization and initialization from a C array. Data access
// via operator[].
//
template <typename T, int N>
struct ym_vec {
    T v[N];  // elements

    ym_vec() {
        for (int i = 0; i < N; i++) v[i] = 0;
    }
    explicit ym_vec(T x) {
        for (int i = 0; i < N; i++) v[i] = x;
    }
    explicit ym_vec(const T* vv) {
        for (int i = 0; i < N; i++) v[i] = vv[i];
    }

    const T& operator[](int i) const { return v[i]; }
    T& operator[](int i) { return v[i]; }
};

//
// Specialization of ym_vec template to allow data access also via named
// members (x, y). Also allow constraction to/from vector of length-1.
//
template <typename T>
struct ym_vec<T, 2> {
    union {
        struct {
            T x, y;
        };
        T v[2];
    };

    ym_vec() : x(), y() {}
    ym_vec(T x, T y) : x(x), y(y) {}
    explicit ym_vec(T s) : x(s), y(s) {}
    explicit ym_vec(const T* vv) : x(vv[0]), y(vv[1]) {}
    template <class U>
    explicit ym_vec(const ym_vec<U, 2>& v) : x(v.x), y(v.y) {}

    const T& operator[](int i) const { return v[i]; }
    T& operator[](int i) { return v[i]; }
};

//
// Specialization of ym_vec template to allow data access also via named
// members (x, y, z). Also allow constraction to/from vector of length-1.
//
template <typename T>
struct ym_vec<T, 3> {
    union {
        struct {
            T x, y, z;
        };
        T v[3];
    };

    ym_vec() : x(), y(), z() {}
    ym_vec(T x, T y, T z) : x(x), y(y), z(z) {}
    explicit ym_vec(T s) : x(s), y(s), z(s) {}
    explicit ym_vec(const T* vv) : x(vv[0]), y(vv[1]), z(vv[2]) {}
    template <class U>
    explicit ym_vec(const ym_vec<U, 3>& v) : x(v.x), y(v.y), z(v.z) {}

    explicit ym_vec(const ym_vec<T, 2>& sv, T z) : ym_vec(sv.x, sv.y, z) {}
    explicit operator ym_vec<T, 2>() const { return ym_vec<T, 2>(x, y); }

    const T& operator[](int i) const { return v[i]; }
    T& operator[](int i) { return v[i]; }
};

//
// Specialization of ym_vec template to allow data access also via named
// members (x, y, z, w). Also allow constraction to/from vector of length-1.
//
template <typename T>
struct ym_vec<T, 4> {
    union {
        struct {
            T x, y, z, w;
        };
        T v[4];
    };

    ym_vec() : x(), y(), z(), w() {}
    ym_vec(T x, T y, T z, T w) : x(x), y(y), z(z), w(w) {}
    explicit ym_vec(T s) : x(s), y(s), z(s), w(s) {}
    explicit ym_vec(const T* vv) : x(vv[0]), y(vv[1]), z(vv[2]), w(vv[3]) {}
    template <class U>
    explicit ym_vec(const ym_vec<U, 4>& v) : x(v.x), y(v.y), z(v.z), w(v.w) {}

    ym_vec(const ym_vec<T, 3>& sv, T w) : ym_vec(sv.x, sv.y, sv.z, w) {}
    explicit operator ym_vec<T, 3>() const { return ym_vec<T, 3>(x, y, z); }

    const T& operator[](int i) const { return v[i]; }
    T& operator[](int i) { return v[i]; }
};

//
// Matrix of element of compile time dimensions stored in column major format.
// Allows for default initializer, initalization with columns and to/from
// smaller matrices.
// Access to columns via operator [] and raw data pointer data().
//
template <typename T, int N, int M>
struct ym_mat {
    typedef ym_vec<T, N> V;

    V v[M];

    ym_mat() {
        for (int j = 0; j < M; j++) v[j] = V();
    }

    explicit ym_mat(const T* vv) {
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) v[j][i] = vv[j * N + i];
        }
    }

    ym_mat(const std::initializer_list<V>& cols) {
        assert(cols.size() == M);
        int j = 0;
        for (auto&& c : cols) v[j++] = c;
    }

    explicit ym_mat(const ym_mat<T, N - 1, M - 1>& m,
                    const ym_vec<T, N - 1>& t) {
        for (int j = 0; j < M - 1; j++) {
            v[j] = V{m.v[j], 0};
        }
        v[M - 1] = V{t, 1};
    }

    explicit operator ym_mat<T, N - 1, M - 1>() const {
        ym_mat<T, N - 1, M - 1> m;
        for (int j = 0; j < M - 1; j++) m.v[j] = (ym_vec<T, N - 1>)v[j];
        return m;
    }

    const V& operator[](int i) const { return v[i]; }
    V& operator[](int i) { return v[i]; }

    T* data() { return &v[0].x; }
    const T* data() const { return &v[0].x; }
};

// -----------------------------------------------------------------------------
// AFFINE AND RIGID TRANSFORMS
// -----------------------------------------------------------------------------

//
// Affine matrix stored as a linear tranform NxN-matrix m and a translation
// N-vector t. Provides converion to/from (N+1)x(N+1) matrices. Can also
// specify whether the matrix is orthonormal at compile time (this in turns
// speeds up inversion and transforms). If orthonormal, than the matrix
// columns are the axis of the coordinate system, while the translation is
// its origin.
//
template <typename T, int N, bool orthonormal>
struct ym_affine {
    typedef ym_vec<T, N> V;
    typedef ym_mat<T, N, N> M;

    M m;
    V t;

    ym_affine() : m(), t() {}
    ym_affine(const M& m, const V& t) : m(m), t(t) {}

    explicit ym_affine(const ym_mat<T, N + 1, N + 1>& mat) {
        for (int j = 0; j < N; j++) m.v[j] = V(mat.v[j]);
        m.t = V(mat.v[N]);
    }
    operator ym_mat<T, N + 1, N + 1>() const {
        return ym_mat<T, N + 1, N + 1>{m, t};
    }
};

//
// Specializaiton of affine matrix with esplicit access to thr rigid frame
// with origin o and axis x, y. Also support a rotation/position names for
// elements.
//
template <typename T, bool orthonormal>
struct ym_affine<T, 2, orthonormal> {
    typedef ym_vec<T, 2> V;
    typedef ym_mat<T, 2, 2> M;

    union {
        struct {
            M m;
            V t;
        };
        struct {
            M rot;
            V pos;
        };
        struct {
            V x, y, o;
        };
    };

    ym_affine() : x(1, 0), y(0, 1), o(0, 0) {}
    ym_affine(const M& m, const V& t) : m(m), t(t) {}
    ym_affine(const V& x, const V& y, const V& o) : x(x), y(y), o(o) {}

    explicit ym_affine(const ym_mat<T, 3, 3>& mat)
        : x(mat[0]), y(mat[1]), o(mat[2]) {}
    operator ym_mat<T, 3, 3>() const { return ym_mat<T, 3, 3>{m, t}; }
};

//
// Specializaiton of affine matrix with esplicit access to thr rigid frame
// with origin o and axis x, y, z. Also support a rotation/position names for
// elements.
//
template <typename T, bool orthonormal>
struct ym_affine<T, 3, orthonormal> {
    typedef ym_vec<T, 3> V;
    typedef ym_mat<T, 3, 3> M;

    union {
        struct {
            M m;
            V t;
        };
        struct {
            M rot;
            V pos;
        };
        struct {
            V x, y, z, o;
        };
        struct {
            V tangu, tangv, norm, pt;
        };
    };

    ym_affine() : x(1, 0, 0), y(0, 1, 0), z(0, 0, 1), o(0, 0, 0) {}
    ym_affine(const M& m, const V& t) : m(m), t(t) {}
    ym_affine(const V& x, const V& y, const V& z, const V& o)
        : x(x), y(y), z(z), o(o) {}

    explicit ym_affine(const ym_mat<T, 4, 4>& mat)
        : x(mat[0]), y(mat[1]), z(mat[2]), o(mat[3]) {}
    explicit operator ym_mat<T, 4, 4>() const { return ym_mat<T, 4, 4>{m, t}; }

    explicit ym_affine(const ym_affine<T, 3, !orthonormal>& af)
        : m(af.m), t(af.t) {}
    explicit operator ym_affine<T, 3, !orthonormal>() { return {m, t}; }
};

// -----------------------------------------------------------------------------
// LINEAR ALGEBRA TYPEDEFS AND CONSTANTS
// -----------------------------------------------------------------------------

//
// Typedef for standard vectors.
//

typedef ym_vec<float, 2> ym_vec2f;
typedef ym_vec<float, 3> ym_vec3f;
typedef ym_vec<float, 4> ym_vec4f;

typedef ym_vec<int, 2> ym_vec2i;
typedef ym_vec<int, 3> ym_vec3i;
typedef ym_vec<int, 4> ym_vec4i;

typedef ym_vec<unsigned char, 2> ym_vec2b;
typedef ym_vec<unsigned char, 3> ym_vec3b;
typedef ym_vec<unsigned char, 4> ym_vec4b;

//
// Vector Constants.
//

const ym_vec2f ym_zero2f = ym_vec2f();
const ym_vec3f ym_zero3f = ym_vec3f();
const ym_vec4f ym_zero4f = ym_vec4f();

const ym_vec2f ym_x2f = ym_vec2f(1, 0);
const ym_vec2f ym_y2f = ym_vec2f(0, 1);

const ym_vec3f ym_x3f = ym_vec3f(1, 0, 0);
const ym_vec3f ym_y3f = ym_vec3f(0, 1, 0);
const ym_vec3f ym_z3f = ym_vec3f(0, 0, 1);

//
// Typedef for standard matrices.
//

typedef ym_mat<float, 2, 2> ym_mat2f;
typedef ym_mat<float, 3, 3> ym_mat3f;
typedef ym_mat<float, 4, 4> ym_mat4f;

//
// Matrix Constants.
//

const ym_mat2f ym_identity_mat2f = ym_mat2f{{1, 0}, {0, 1}};
const ym_mat3f ym_identity_mat3f = ym_mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const ym_mat4f ym_identity_mat4f =
    ym_mat4f{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

//
// Typedef for standard affine matrices.
//

typedef ym_affine<float, 2, false> ym_affine2f;
typedef ym_affine<float, 3, false> ym_affine3f;

//
// Affine Transforms Constants.
//

const ym_affine2f ym_identity_affine2f = ym_affine2f{{1, 0}, {0, 1}, {0, 0}};
const ym_affine3f ym_identity_affine3f =
    ym_affine3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

//
// Typedef for standard rigid transforms / frames.
//

typedef ym_affine<float, 2, true> ym_frame2f;
typedef ym_affine<float, 3, true> ym_frame3f;

//
// Rigid Transforms/Frames Constants.
//

const ym_frame2f ym_identity_frame2f = ym_frame2f{{1, 0}, {0, 1}, {0, 0}};
const ym_frame3f ym_identity_frame3f =
    ym_frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// -----------------------------------------------------------------------------
// LINEAR ALGEBRA OPERATIONS
// -----------------------------------------------------------------------------

//
// Components-wise comparison
//

template <typename T, int N>
inline bool ym_iszero(const ym_vec<T, N>& a) {
    for (int i = 0; i < N; i++)
        if (a.v[i] != 0) return false;
    return true;
}

template <typename T, int N>
inline bool ym_isfinite(const ym_vec<T, N>& a) {
    for (int i = 0; i < N; i++)
        if (!ym_isfinite(a.v[i])) return false;
    return true;
}

template <typename T, int N>
inline bool operator==(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    for (int i = 0; i < N; i++)
        if (a.v[i] != b.v[i]) return false;
    return true;
}

template <typename T, int N>
inline bool operator<(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    for (int i = 0; i < N; i++)
        if (a.v[i] >= b.v[i]) return false;
    return true;
}

//
// Component-wise arithmentic.
//

template <typename T, int N>
inline ym_vec<T, N> operator-(const ym_vec<T, N>& a) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = -a.v[i];
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> operator+(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] + b.v[i];
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> operator-(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] - b.v[i];
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> operator*(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] * b.v[i];
    return c;
}

template <typename T, int N, typename T1>
inline ym_vec<T, N> operator*(const ym_vec<T, N>& a, const T1& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] * b;
    return c;
}

template <typename T1, typename T, int N>
inline ym_vec<T, N> operator*(const T1& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a * b.v[i];
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> operator/(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] / b.v[i];
    return c;
}

template <typename T, int N, typename T1>
inline ym_vec<T, N> operator/(const ym_vec<T, N>& a, const T1 b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a.v[i] / b;
    return c;
}

template <typename T1, typename T, int N>
inline ym_vec<T, N> operator/(const T1& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = a / b.v[i];
    return c;
}

//
// Component-wise assignment arithmentic.
//

template <typename T, int N>
inline ym_vec<T, N>& operator+=(ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    a = a + b;
    return a;
}

template <typename T, int N>
inline ym_vec<T, N>& operator-=(ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    a = a - b;
    return a;
}

template <typename T, int N>
inline ym_vec<T, N>& operator*=(ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    a = a * b;
    return a;
}

template <typename T, int N, typename T1>
inline ym_vec<T, N>& operator*=(ym_vec<T, N>& a, const T1& b) {
    a = a * b;
    return a;
}

template <typename T, int N>
inline ym_vec<T, N>& operator/=(ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    a = a / b;
    return a;
}

template <typename T, int N, typename T1>
inline ym_vec<T, N>& operator/=(ym_vec<T, N>& a, const T1 b) {
    a = a / b;
    return a;
}

//
// Component-wise min/max and clamping.
//

template <typename T, int N>
inline ym_vec<T, N> ym_min(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = ym_min(a.v[i], b.v[i]);
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> ym_max(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = ym_max(a.v[i], b.v[i]);
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> ym_min(const ym_vec<T, N>& a, const T& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = ym_min(a.v[i], b);
    return c;
}

template <typename T, int N>
inline ym_vec<T, N> ym_max(const ym_vec<T, N>& a, const T& b) {
    ym_vec<T, N> c;
    for (int i = 0; i < N; i++) c.v[i] = ym_max(a.v[i], b);
    return c;
}

template <typename T, int M>
inline ym_vec<T, M> ym_clamp(const ym_vec<T, M>& x, const T& m, const T& N) {
    ym_vec<T, M> c;
    for (int i = 0; i < M; i++) c.v[i] = ym_clamp(x.v[i], m, N);
    return c;
}

template <typename T, int M>
inline ym_vec<T, M> ym_clamplen(const ym_vec<T, M> x, float N) {
    float l = ym_length(x);
    return (l > N) ? x * N / l : x;
}

//
// Element min/max.
//
template <typename T, int N>
inline int ym_min_element(const ym_vec<T, N>& a) {
    T v = std::numeric_limits<T>::max();
    int pos = -1;
    for (int i = 0; i < N; i++) {
        if (v > a.v[i]) {
            v = a.v[i];
            pos = i;
        }
    }
    return pos;
}
template <typename T, int N>
inline int ym_max_element(const ym_vec<T, N>& a) {
    T v = -std::numeric_limits<T>::max();
    int pos = -1;
    for (int i = 0; i < N; i++) {
        if (v < a.v[i]) {
            v = a.v[i];
            pos = i;
        }
    }
    return pos;
}

//
// Vector operations
//
template <typename T, int N>
inline T ym_dot(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    float c = 0;
    for (int i = 0; i < N; i++) c += a.v[i] * b.v[i];
    return c;
}

template <typename T, int N>
inline T ym_length(const ym_vec<T, N>& a) {
    return sqrt(ym_dot(a, a));
}

template <typename T, int N>
inline T ym_lengthsqr(const ym_vec<T, N>& a) {
    return ym_dot(a, a);
}

template <typename T, int N>
inline ym_vec<T, N> ym_normalize(const ym_vec<T, N>& a) {
    T l = ym_length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

template <typename T, int N>
inline T ym_dist(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    return ym_length(a - b);
}

template <typename T, int N>
inline T ym_distsqr(const ym_vec<T, N>& a, const ym_vec<T, N>& b) {
    return ym_lengthsqr(a - b);
}

template <typename T>
inline T ym_cross(const ym_vec<T, 2>& a, const ym_vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}

template <typename T>
inline ym_vec<T, 3> ym_cross(const ym_vec<T, 3>& a, const ym_vec<T, 3>& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
            a.x * b.y - a.y * b.x};
}

// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
template <typename T>
inline ym_vec<T, 3> ym_orthogonal(const ym_vec<T, 3>& v) {
    return fabs(v.x) > fabs(v.z) ? ym_vec3f{-v.y, v.x, 0}
                                 : ym_vec3f{0, -v.z, v.y};
}

template <typename T>
inline ym_vec<T, 3> ym_orthonormalize(const ym_vec<T, 3>& a,
                                      const ym_vec<T, 3>& b) {
    return ym_normalize(a - b * ym_dot(a, b));
}

//
// Aggregate operations.
//
template <typename T, int M>
inline T ym_sum(const ym_vec<T, M>& a) {
    T s = 0;
    for (int i = 0; i < M; i++) s += a.v[i];
    return s;
}
template <typename T, int M>
inline T ym_mean(const ym_vec<T, M>& a) {
    return ym_sum(a) / M;
}

//
// Components-wise comparison
//

template <typename T, int N, int M>
inline bool operator==(const ym_mat<T, N, M>& a, const ym_mat<T, N, M>& b) {
    for (int i = 0; i < M; i++)
        if (!(a.v[i] == b.v[i])) return false;
    return true;
}

//
// Components-wise arithmetic
//

template <typename T, int N, int M>
inline ym_mat<T, M, N> operator-(const ym_mat<T, N, M>& a) {
    ym_mat<T, N, M> c;
    for (int i = 0; i < M; i++) c.v[i] = -a.v[i];
    return c;
}

template <typename T, int N, int M>
inline ym_mat<T, M, N> operator+(const ym_mat<T, N, M>& a,
                                 const ym_mat<T, N, M>& b) {
    ym_mat<T, N, M> c;
    for (int i = 0; i < M; i++) c.v[i] = a.v[i] + b.v[i];
    return c;
}

//
// Linear Algebra operations.
//

template <typename T, int N, int M>
inline ym_mat<T, M, N> ym_transpose(const ym_mat<T, N, M>& a) {
    ym_mat<T, M, N> c;
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            c.v[i][j] = a.v[j][i];
        }
    }
    return c;
}

// http://stackoverflow.com/questions/983999/simple-3x3-matrix-inverse-code-c
template <typename T>
inline ym_mat<T, 3, 3> ym_inverse(const ym_mat<T, 3, 3>& m_) {
    T m[3][3];
    *(ym_mat<T, 3, 3>*)m = m_;

    // computes the inverse of a matrix m
    T det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
            m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
            m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

    T invdet = 1 / det;

    T minv[3][3];  // inverse of matrix m
    minv[0][0] = (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invdet;
    minv[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invdet;
    minv[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invdet;
    minv[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invdet;
    minv[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invdet;
    minv[1][2] = (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invdet;
    minv[2][0] = (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invdet;
    minv[2][1] = (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invdet;
    minv[2][2] = (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invdet;

    return *(ym_mat<T, 3, 3>*)minv;
}

// http://stackoverflow.com/questions/2624422/efficient-4x4-matrix-inverse-affine-transform
template <typename T>
inline ym_mat<T, 4, 4> ym_inverse(const ym_mat<T, 4, 4>& m) {
    T a[4][4];
    *(ym_mat<T, 4, 4>*)a = m;

    T s0 = a[0][0] * a[1][1] - a[1][0] * a[0][1];
    T s1 = a[0][0] * a[1][2] - a[1][0] * a[0][2];
    T s2 = a[0][0] * a[1][3] - a[1][0] * a[0][3];
    T s3 = a[0][1] * a[1][2] - a[1][1] * a[0][2];
    T s4 = a[0][1] * a[1][3] - a[1][1] * a[0][3];
    T s5 = a[0][2] * a[1][3] - a[1][2] * a[0][3];

    T c5 = a[2][2] * a[3][3] - a[3][2] * a[2][3];
    T c4 = a[2][1] * a[3][3] - a[3][1] * a[2][3];
    T c3 = a[2][1] * a[3][2] - a[3][1] * a[2][2];
    T c2 = a[2][0] * a[3][3] - a[3][0] * a[2][3];
    T c1 = a[2][0] * a[3][2] - a[3][0] * a[2][2];
    T c0 = a[2][0] * a[3][1] - a[3][0] * a[2][1];

    // TODO: Should check for 0 determinant
    T invdet =
        1.0f / (s0 * c5 - s1 * c4 + s2 * c3 + s3 * c2 - s4 * c1 + s5 * c0);

    T b[4][4];

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

    return *(ym_mat<T, 4, 4>*)b;
}

//
// Vector/Matrix multiplies
//

template <typename T, int N, int M>
inline ym_mat<T, M, N> operator*(const ym_mat<T, N, M>& a, T b) {
    ym_mat<T, N, M> c;
    for (int i = 0; i < M; i++) c.v[i] = a.v[i] * b;
    return c;
}

template <typename T, int N, int M>
inline ym_vec<T, N> operator*(const ym_mat<T, N, M>& a, const ym_vec<T, M>& b) {
    ym_vec<T, N> c;
    for (int j = 0; j < M; j++) c += a.v[j] * b.v[j];
    return c;
}

template <typename T, int N, int M, int K>
inline ym_mat<T, N, M> operator*(const ym_mat<T, N, K>& a,
                                 const ym_mat<T, K, M>& b) {
    ym_mat<T, N, M> c;
    for (int j = 0; j < M; j++) c.v[j] = a * b.v[j];
    return c;
}

// -----------------------------------------------------------------------------
// TRANSFORM OPERATIONS
// -----------------------------------------------------------------------------

//
// Component-wise comparisons
//

template <typename T, int N, bool ortho>
inline bool operator==(const ym_affine<T, N, ortho>& a,
                       const ym_affine<T, N, ortho>& b) {
    return a.m == b.m && a.t == b.t;
}

//
// Vector/matrix multiplies for transforms.
//

template <typename T, int N, bool ortho>
inline ym_vec<T, N> operator*(const ym_affine<T, N, ortho>& a,
                              const ym_vec<T, N>& b) {
    return a.m * b + a.t;
}

template <typename T, int N, bool ortho>
inline ym_affine<T, N, ortho> operator*(const ym_affine<T, N, ortho>& a,
                                        const ym_affine<T, N, ortho>& b) {
    return {a.m * b.m, a.m * b.t + a.t};
}

//
// Tranform operations.
//

template <typename T, int N>
inline ym_affine<T, N, false> ym_inverse(const ym_affine<T, N, false>& a) {
    ym_mat<T, N, N> minv = ym_inverse(a.m);
    return {minv, -(minv * a.t)};
}

template <typename T, int N>
inline ym_affine<T, N, true> ym_inverse(const ym_affine<T, N, true>& a) {
    ym_mat<T, N, N> minv = ym_inverse(a.m);
    return {minv, -(minv * a.t)};
}

template <typename T, int N>
inline ym_vec<T, N> ym_transform_point(const ym_mat<T, N + 1, N + 1>& a,
                                       const ym_vec<T, N>& b) {
    // make it generic
    ym_vec<T, N + 1> vb = {b, 1};
    ym_vec<T, N + 1> tvb = a * vb;
    return ym_vec<T, N>(tvb) / tvb.w;
}

template <typename T, int N>
inline ym_vec<T, N> ym_transform_vector(const ym_mat<T, N + 1, N + 1>& a,
                                        const ym_vec<T, N>& b) {
    // make it generic
    ym_vec<T, N + 1> vb = {b, 0};
    ym_vec<T, N + 1> tvb = a * vb;
    return ym_vec<T, N>(tvb);
}

template <typename T, int N>
inline ym_vec<T, N> ym_transform_direction(const ym_mat<T, N + 1, N + 1>& a,
                                           const ym_vec<T, N>& b) {
    return ym_normalize(ym_transform_vector(a, b));
}

template <typename T, int N, bool ortho>
inline ym_vec<T, N> ym_transform_point(const ym_affine<T, N, ortho>& a,
                                       const ym_vec<T, N>& b) {
    return a * b;
}

template <typename T, int N, bool ortho>
inline ym_vec<T, N> ym_transform_vector(const ym_affine<T, N, ortho>& a,
                                        const ym_vec<T, N>& b) {
    return a.m * b;
}

template <typename T, int N>
inline ym_vec<T, N> ym_transform_direction(const ym_affine<T, N, false>& a,
                                           const ym_vec<T, N>& b) {
    return ym_normalize(a.m * b);
}

template <typename T, int N>
inline ym_vec<T, N> ym_transform_direction(const ym_affine<T, N, true>& a,
                                           const ym_vec<T, N>& b) {
    return a.m * b;
}

// -----------------------------------------------------------------------------
// TRANSFORM MATRICES
// -----------------------------------------------------------------------------

//
// Rotation matrix crom axis/angle
//
template <typename T>
inline ym_mat<T, 3, 3> ym_rotation_mat3(const ym_vec<T, 3>& axis, T angle) {
    float s = sin(angle), c = cos(angle);
    ym_vec<T, 3> vv = ym_normalize(axis);
    return ym_transpose(ym_mat<T, 3, 3>{
        {c + (1 - c) * vv.x * vv.x, (1 - c) * vv.x * vv.y - s * vv.z,
         (1 - c) * vv.x * vv.z + s * vv.y},
        {(1 - c) * vv.x * vv.y + s * vv.z, c + (1 - c) * vv.y * vv.y,
         (1 - c) * vv.y * vv.z - s * vv.x},
        {(1 - c) * vv.x * vv.z - s * vv.y, (1 - c) * vv.y * vv.z + s * vv.x,
         c + (1 - c) * vv.z * vv.z}});
}

//
// Translation transform
//
template <typename T>
inline ym_affine<T, 3, true> ym_translation_xform3(const ym_vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}

//
// Scaling transform
//
template <typename T>
inline ym_affine<T, 3, false> ym_scaling_xform3(const ym_vec<T, 3>& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}

//
// Rotation transform
//
template <typename T>
inline ym_affine<T, 3, true> ym_rotation_xform3(const ym_vec<T, 3>& axis,
                                                T angle) {
    return ym_affine<T, 3, true>{ym_rotation_mat3(axis, angle), {0, 0, 0}};
}

//
// Lookat tranform
//
template <typename T>
inline ym_affine<T, 3, true> ym_lookat_xform3(const ym_vec<T, 3>& eye,
                                              const ym_vec<T, 3>& center,
                                              const ym_vec<T, 3>& up) {
    ym_vec<T, 3> w = ym_normalize(eye - center);
    ym_vec<T, 3> u = ym_normalize(ym_cross(up, w));
    ym_vec<T, 3> v = ym_normalize(ym_cross(w, u));
    return {u, v, w, eye};
}

//
// Frame/transform for origin and z.
//
template <typename T>
inline ym_affine<T, 3, true> ym_make_frame(const ym_vec<T, 3>& o,
                                           const ym_vec<T, 3>& z_) {
    ym_vec<T, 3> z = ym_normalize(z_);
    ym_vec<T, 3> x = ym_normalize(ym_orthogonal(z));
    ym_vec<T, 3> y = ym_normalize(ym_cross(z, x));
    return {x, y, z, o};
}

//
// OpenGL perspective frustum matrix
//
template <typename T>
inline ym_mat<T, 4, 4> ym_frustum_mat4(T l, T r, T b, T t, T n, T f) {
    return ym_transpose(
        ym_mat<T, 4, 4>{{2 * n / (r - l), 0, (r + l) / (r - l), 0},
                        {0, 2 * n / (t - b), (t + b) / (t - b), 0},
                        {0, 0, -(f + n) / (f - n), -2 * f * n / (f - n)},
                        {0, 0, -1, 0}});
}

//
// OpenGL orthographic matrix
//
template <typename T>
inline ym_mat<T, 4, 4> ym_ortho_mat4(T l, T r, T b, T t, T n, T f) {
    return ym_transpose(
        ym_mat<T, 4, 4>{{2 / (r - l), 0, 0, -(r + l) / (r - l)},
                        {0, 2 / (t - b), 0, -(t + b) / (t - b)},
                        {0, 0, -2 / (f - n), -(f + n) / (f - n)},
                        {0, 0, 0, 1}});
}

//
// OpenGL orthographic 2D matrix
//
template <typename T>
inline ym_mat<T, 4, 4> ym_ortho2d_mat4(T l, T r, T b, T t) {
    return ym_transpose(ym_mat<T, 4, 4>{{2 / (r - l), 0, 0, -(r + l) / (r - l)},
                                        {0, 2 / (t - b), 0, -(t + b) / (t - b)},
                                        {0, 0, -1, 0},
                                        {0, 0, 0, 1}});
}

//
// OpenGL perspective matrix
//
template <typename T>
inline ym_mat<T, 4, 4> ym_perspective_mat4(T fovy, T aspect, T near, T far) {
    T f = 1 / tanf(fovy / 2);
    return ym_transpose(ym_mat<T, 4, 4>{
        {f / aspect, 0, 0, 0},
        {0, f, 0, 0},
        {0, 0, (far + near) / (near - far), 2 * far * near / (near - far)},
        {0, 0, -1, 0}});
}

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------

//
// Axis aligned bounding box.
//
template <typename T>
struct ym_range {
    union {
        struct {
            T min, max;
        };  // bouding box min and max corners
        T d[2];
    };

    ym_range()
        : min(std::numeric_limits<T>::max()),
          max(std::numeric_limits<T>::lowest()) {}
    ym_range(const T& m, const T& N) : min(m), max(N) {}

    T size() const { return max - min; }
    T center() const { return (max + min) / 2; }

    T& operator[](int i) { return d[i]; }
    const T& operator[](int i) const { return d[i]; }
};

//
// Axis aligned bounding box for vec spaces.
//
template <typename T, int M>
struct ym_range<ym_vec<T, M>> {
    union {
        struct {
            ym_vec<T, M> min, max;
        };  // bouding box min and max corners
        ym_vec<T, M> d[2];
    };

    ym_range()
        : min(std::numeric_limits<T>::max()),
          max(std::numeric_limits<T>::lowest()) {}
    ym_range(const ym_vec<T, M>& m, const ym_vec<T, M>& N) : min(m), max(N) {}

    ym_vec<T, M> size() const { return max - min; }
    ym_vec<T, M> center() const { return (max + min) / 2; }

    ym_vec<T, M>& operator[](int i) { return d[i]; }
    const ym_vec<T, M>& operator[](int i) const { return d[i]; }
};

//
// Axis aligned bounding box typedefs.
//

typedef ym_range<ym_vec2f> ym_range2f;
typedef ym_range<ym_vec3f> ym_range3f;

typedef ym_range<int> ym_range1i;
typedef ym_range<ym_vec2i> ym_range2i;
typedef ym_range<ym_vec3i> ym_range3i;

//
// Axis aligned bounding box constants.
//

const ym_range2f ym_invalid_range2f = ym_range2f();
const ym_range3f ym_invalid_range3f = ym_range3f();

//
// Axis aligned bounding box operations.
//

template <typename T>
inline ym_range<T> ym_rexpand(const ym_range<T>& a, const T& b) {
    return {ym_min(a.min, b), ym_max(a.max, b)};
}

template <typename T>
inline ym_range<T> ym_rexpand(const ym_range<T>& a, const ym_range<T>& b) {
    return {ym_min(a.min, b.min), ym_max(a.max, b.max)};
}

template <typename T>
inline ym_range<T> operator+(const ym_range<T>& a, const T& b) {
    return ym_rexpand(a, b);
}

template <typename T>
inline ym_range<T> operator+(const ym_range<T>& a, const ym_range<T>& b) {
    return ym_rexpand(a, b);
}

template <typename T>
inline ym_range<T>& operator+=(ym_range<T>& a, const T& b) {
    a = ym_rexpand(a, b);
    return a;
}

template <typename T>
inline ym_range<T>& operator+=(ym_range<T>& a, const ym_range<T>& b) {
    a = ym_rexpand(a, b);
    return a;
}

template <typename T>
inline T ym_rcenter(const ym_range<T>& a) {
    return (a.min + a.max) / 2;
}

template <typename T>
inline T ym_rsize(const ym_range<T>& a) {
    return a.max - a.min;
}

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------

//
// Rays with origin, direction and min/max t value.
//
template <typename T, int M>
struct ym_ray {
    ym_vec<T, M> o;
    ym_vec<T, M> d;
    T tmin;
    T tmax;

    ym_ray() : o(), d(0, 0, 1), tmin(0), tmax(std::numeric_limits<T>::max()) {}
    ym_ray(const ym_vec<T, M>& o, const ym_vec<T, M>& d, T tmin = 0,
           T tmax = std::numeric_limits<T>::max())
        : o(o), d(d), tmin(tmin), tmax(tmax) {}

    ym_vec<T, M> eval(T t) const { return o + t * d; }
};

//
// Rays typedefs.
//

typedef ym_ray<float, 2> ym_ray2f;
typedef ym_ray<float, 3> ym_ray3f;

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------

//
// Turntable for UI navigation from a from/to/up parametrization of the camera.
//
static inline void ym_turntable(ym_vec3f* from, ym_vec3f* to, ym_vec3f* up,
                                const ym_vec2f& rotate, float dolly,
                                const ym_vec2f& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        ym_vec3f z = ym_normalize(*to - *from);
        float lz = ym_dist(*to, *from);
        float phi = atan2f(z.z, z.x) + rotate[0];
        float theta = acosf(z.y) + rotate[1];
        theta = fmax(0.001f, fmin(theta, ym_pif - 0.001f));
        ym_vec3f nz = {sinf(theta) * cosf(phi) * lz, cosf(theta) * lz,
                       sinf(theta) * sinf(phi) * lz};
        *from = *to - nz;
    }

    // dolly if necessary
    if (dolly) {
        ym_vec3f z = ym_normalize(*to - *from);
        float lz = fmaxf(0.001f, ym_dist(*to, *from) * (1 + dolly));
        z *= lz;
        *from = *to - z;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        ym_vec3f z = ym_normalize(*to - *from);
        ym_vec3f x = ym_normalize(ym_cross(*up, z));
        ym_vec3f y = ym_normalize(ym_cross(z, x));
        ym_vec3f t = {pan[0] * x.x + pan[1] * y.x, pan[0] * x.y + pan[1] * y.y,
                      pan[0] * x.z + pan[1] * y.z};
        *from += t;
        *to += t;
    }
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// PCG random numbers. A family of random number generators that supports
// multiple sequences. In our code, we allocate one sequence for each sample.
// PCG32 from http://www.pcg-random.org/
//

//
// Random number state (PCG32)
//
struct ym_rng_pcg32 {
    uint64_t state, inc;
};

//
// Next random number
//
inline uint32_t ym_rng_next(ym_rng_pcg32* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ull + (rng->inc | 1u);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-((int32_t)rot)) & 31));
}

//
// Init a random number generator with a state state from the sequence seq.
//
inline void ym_rng_init(ym_rng_pcg32* rng, uint64_t state, uint64_t seq) {
    rng->state = 0U;
    rng->inc = (seq << 1u) | 1u;
    ym_rng_next(rng);
    rng->state += state;
    ym_rng_next(rng);
}

//
// Next random float in [0,1).
//
inline float ym_rng_nextf(ym_rng_pcg32* rng) {
    return (float)ldexp(ym_rng_next(rng), -32);
}

//
// Next random float in [0,1)x[0,1).
//
inline ym_vec2f ym_rng_next2f(ym_rng_pcg32* rng) {
    return {ym_rng_nextf(rng), ym_rng_nextf(rng)};
}

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------

//
// Computes the i-th term of a permutation of l values keyed by p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
inline uint32_t ym_hash_permute(uint32_t i, uint32_t n, uint32_t key) {
    uint32_t w = n - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;
    do {
        i ^= key;
        i *= 0xe170893du;
        i ^= key >> 16;
        i ^= (i & w) >> 4;
        i ^= key >> 8;
        i *= 0x0929eb3f;
        i ^= key >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | key >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;
    } while (i >= n);
    return (i + key) % n;
}

//
// Computes a float value by hashing i with a key p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
inline float ym_hash_randfloat(uint32_t i, uint32_t key) {
    i ^= key;
    i ^= i >> 17;
    i ^= i >> 10;
    i *= 0xb36534e5;
    i ^= i >> 12;
    i ^= i >> 21;
    i *= 0x93fc4795;
    i ^= 0xdf6e307f;
    i ^= i >> 17;
    i *= 1 | key >> 18;
    return i * (1.0f / 4294967808.0f);
}

//
// 64 bit integer hash. Public domain code.
//
inline uint64_t ym_hash_uint64(uint64_t a) {
    a = (~a) + (a << 21);  // a = (a << 21) - a - 1;
    a ^= (a >> 24);
    a += (a << 3) + (a << 8);  // a * 265
    a ^= (a >> 14);
    a += (a << 2) + (a << 4);  // a * 21
    a ^= (a >> 28);
    a += (a << 31);
    return a;
}

//
// 64-to-32 bit integer hash. Public domain code.
//
static inline uint32_t ym_hash_uint64_32(uint64_t a) {
    a = (~a) + (a << 18);  // a = (a << 18) - a - 1;
    a ^= (a >> 31);
    a *= 21;  // a = (a + (a << 2)) + (a << 4);
    a ^= (a >> 11);
    a += (a << 6);
    a ^= (a >> 22);
    return (uint32_t)a;
}

// -----------------------------------------------------------------------------
// SIMPLE CONTAINERS
// -----------------------------------------------------------------------------

//
// Dynamic array with semantic similar to std::vector but only a subset of the
// operations. See std::vector for the individual functions.
//
template <typename T>
struct ym_darray {
    ym_darray() : num(0), cap(0), d(0) {}

    ym_darray(size_t n) : num(0), cap(0), d(0) { resize(n); }

    ym_darray(T* start, T* end) : num(0), cap(0), d(0) { assign(start, end); }

    ym_darray(const std::initializer_list<T>& lst) : num(0), cap(0), d(0) {
        assign(lst);
    }

    ym_darray(const ym_darray<T>& a) : num(0), cap(0), d(0) {
        assign(a.data(), a.data() + a.size());
    }

    ym_darray(ym_darray<T>&& a) : num(a.num), cap(a.cap), d(a.d) {
        a.num = a.cap = 0;
        a.d = nullptr;
    }

    ~ym_darray() {
        if (d) delete[] d;
    }

    ym_darray<T>& operator=(const ym_darray<T>& a) {
        if (&a == this) return *this;
        assign(a.data(), a.data() + a.size());
        return *this;
    }

    ym_darray<T>& operator=(ym_darray<T>&& a) {
        if (&a == this) return *this;
        if (d) delete[] d;
        d = a.d;
        num = a.num;
        cap = a.cap;
        a.num = a.cap = 0;
        a.d = nullptr;
        return *this;
    }

    bool operator==(const ym_darray<T>& b) const {
        if (size() != b.size()) return false;
        for (size_t i = 0; i < size(); i++) {
            if (d[i] != b.d[i]) return false;
        }
        return true;
    }

    bool operator<(const ym_darray<T>& b) const {
        for (size_t i = 0; i < ym_min(size(), b.size()); i++) {
            if (d[i] >= b.d[i]) return false;
        }
        if (size() >= b.size()) return false;
        return true;
    }

    void assign(size_t n, const T& v) {
        resize(n);
        _set(v);
    }

    void assign(const T* start, const T* end) {
        resize(end - start);
        _copy(end - start, start);
    }

    void assign(const std::initializer_list<T>& lst) {
        resize(lst.size());
        size_t c = 0;
        for (auto&& v : lst) d[c++] = v;
    }

    void push_back(const T& v) {
        _grow_upto(num + 1);
        d[num++] = v;
    }

    void resize(size_t n) {
        if (n > cap) _realloc(n);
        num = n;
    }

    void reserve(size_t n) {
        if (n > cap) _realloc(n);
    }

    void pop_back() { num--; }

    void append(size_t n, const T* dta) {
        reserve(size() + n);
        for (int i = 0; i < n; i++) d[num++] = dta[i];
    }

    void clear() { num = 0; }

    void swap(ym_darray<T>& other) { _swap(other); }

    void shrink_to_fit() {
        if (num == cap) return;
        if (num == 0) {
            if (d) delete[] d;
            d = nullptr;
            cap = 0;
        } else {
            _realloc(num);
        }
    }

    bool empty() const { return num == 0; }

    size_t size() const { return num; }

    int findpos(const T& v) const {
        for (int i = 0; i < num; i++)
            if (d[i] == v) return (int)i;
        return -1;
    }

    T& at(int i) {
        assert(i >= 0 && i < num);
        return d[i];
    }

    const T& at(size_t i) const {
        assert(i >= 0 && i < num);
        return d[i];
    }

    T& operator[](size_t i) { return d[i]; }
    const T& operator[](size_t i) const { return d[i]; }

    T& front() { return d[0]; }
    const T& front() const { return d[0]; }

    T& back() { return d[num - 1]; }
    const T& back() const { return d[num - 1]; }

    T* data() {
        if (num)
            return d;
        else
            return nullptr;
    }
    const T* data() const {
        if (num)
            return d;
        else
            return nullptr;
    }

    T* begin() { return d; }
    const T* begin() const { return d; }

    T* end() { return d + num; }
    const T* end() const { return d + num; }

   private:
    size_t num = 0;
    size_t cap = 0;
    T* d = nullptr;

    void _realloc(size_t n) {
        cap = n;
        if (n) {
            T* old_d = d;
            d = new T[cap];
            // _copy(num, old_d);
            _move(num, old_d);
            delete[] old_d;
        } else {
            if (d) delete[] d;
            d = nullptr;
        }
    }

    void _set(const T& v) {
        for (size_t i = 0; i < num; i++) d[i] = v;
    }

    void _copy(size_t n, const T* dta) {
        if (std::is_trivially_copyable<T>::value) {
            memcpy(d, dta, sizeof(T) * n);
        } else {
            for (size_t i = 0; i < n; i++) d[i] = dta[i];
        }
    }

    void _move(size_t n, T* dta) {
        if (std::is_trivially_copyable<T>::value) {
            memcpy(d, dta, sizeof(T) * n);
        } else {
            for (size_t i = 0; i < n; i++) d[i] = std::move(dta[i]);
        }
    }

    void _grow_upto(size_t n) {
        if (n == 0) return;
        if (n <= cap) return;
        int c = 1;
        while (c < n) c *= 2;
        _realloc(c);
    }

    void _swap(ym_darray<T>& a) {
        ym_swap(d, a.d);
        ym_swap(num, a.num);
        ym_swap(cap, a.cap);
    }
};

//
// Dynamic string with std::string like semantic. See std::string for member
// documentation.
//
struct ym_str {
    ym_str() {
        _buf.resize(1);
        _buf[0] = 0;
    }

    ym_str(const char* ns) { assign(ns); }

    ym_str(const ym_str& ns) : _buf(ns._buf) {}
    ym_str(ym_str&& ns) : _buf(std::move(ns._buf)) {}

    ym_str& operator=(const ym_str& ns) {
        _buf = ns._buf;
        return *this;
    }

    ym_str& operator=(ym_str&& ns) {
        _buf = std::move(ns._buf);
        return *this;
    }

    bool operator==(const ym_str& ns) { return !strcmp(c_str(), ns.c_str()); }

    bool operator==(const char* ns) { return !strcmp(c_str(), ns); }

    bool operator<(const ym_str& ns) { return strcmp(c_str(), ns.c_str()) < 0; }

    bool empty() const { return _buf.size() < 2; }
    size_t size() const { return _buf.size() - 1; }
    size_t length() const { return _buf.size() - 1; }

    void assign(size_t n, const char* ns) {
        _buf.assign(ns, ns + n);
        _buf.push_back(0);
    }

    void assign(const char* ns) { assign(strlen(ns), ns); }

    void assign(const ym_str& ns) { assign(ns.size(), ns.data()); }

    void append(size_t n, const char* ns) {
        _buf.pop_back();
        _buf.append(n, ns);
        _buf.push_back(0);
    }

    void append(const ym_str& ns) { append(ns.size(), ns.data()); }

    void append(const char* ns) { append(strlen(ns), ns); }

    const char* c_str() const { return _buf.data(); }

    char* data() { return _buf.data(); }
    const char* data() const { return _buf.data(); }

   private:
    ym_darray<char> _buf;
};

//
// Key-value pair for std::pair substitution.
//
template <typename TK, typename TV>
struct ym_keyvalue {
    union {
        TK key;
        TK first;
    };
    union {
        TV value;
        TV second;
    };
};

// -----------------------------------------------------------------------------
// CONTAINER OPERATIONS
// -----------------------------------------------------------------------------

//
// String append operations
//

inline ym_str& operator+=(ym_str& a, const char* b) {
    a.append(b);
    return a;
}

inline ym_str& operator+=(ym_str& a, const ym_str& b) {
    a.append(b);
    return a;
}

inline ym_str operator+(const ym_str& a, const ym_str& b) {
    ym_str c = a;
    c += b;
    return c;
}

inline ym_str operator+(const ym_str& a, const char* b) {
    ym_str c = a;
    c += b;
    return c;
}

// -----------------------------------------------------------------------------
// CONTAINER TYPEDEFS
// -----------------------------------------------------------------------------

//
// This code selects whether to use the builtin YOCTO containers os the STL
// ones.
//

#ifdef YGL_USESTL

#include <string>
#include <vector>

using ym_string = std::string;
template <typename T>
using ym_vector = std::vector<T>;
template <typename TK, typename TV>
using ym_pair = std::pair<TK, TV>;

#else

using ym_string = ym_str;
template <typename T>
using ym_vector = ym_darray<T>;
template <typename TK, typename TV>
using ym_pair = ym_keyvalue<TK, TV>;

#endif

// -----------------------------------------------------------------------------
// CONTAINER OPERATIONS
// -----------------------------------------------------------------------------

//
// Vector append operations using operators (Python-like).
//

template <typename T>
inline ym_vector<T>& operator+=(ym_vector<T>& a, const T& b) {
    a.push_back(b);
    return a;
}

template <typename T>
inline ym_vector<T>& operator+=(ym_vector<T>& a, const ym_vector<T>& b) {
    for (int i = 0; i < b.size(); i++) a.push_back(b[i]);
    return a;
}

template <typename T>
inline ym_vector<T> operator+(const ym_vector<T>& a, const ym_vector<T>& b) {
    ym_vector<T> c = a;
    c += b;
    return c;
}

template <typename T>
inline ym_vector<T> operator+(const ym_vector<T>& a, const T& b) {
    ym_vector<T> c = a;
    c += b;
    return c;
}

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------

#include <chrono>

struct ym_timer {
    ym_timer(bool autostart = true) {
        if (autostart) start();
    }

    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

#endif
