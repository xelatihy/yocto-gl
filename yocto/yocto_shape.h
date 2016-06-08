//
// YOCTO_SHAPE: A set of utilities to manipulate 3D shapes represented as
// collection of elements.
//

//
// INCLUDED UTILITIES:
//
// 0. include this file (more compilation options below)
// 1. smoothed normals
//   1.a compute smoothed vertex normals or tangets
//      ys_compute_normals(shape elements, vertex pos, vertex normals)
// 2. shape tesselation
//   2.a tesselate the shape once by splitting each element along its edges
//      ys_tesselate_shape(in/out elements data, in/out vertex data)
// 3. create shapes parametrically only useful for testing)
//    3.a. make a few stadard shapes
//      ys_make_stdshape(standard shape data, our elements, out vertex)
//    3.b. make a parametric grid of elements and a uv array
//      - useful to generate parametric shape
//      ys_make_uvgrid(grid data, out elements, out uv grid)
//    3.c. make a parametric shape via function callbacks
//      ys_make_uvshape(grid data, out vertex data, out elements)
// 4. pick points on a shape
//    - sample shapes one point at a time with
//      - first, compute shape cdf
//      ys_sample_shape_cdf(shape definition, out cdf)
//      - generate element id and element uv for each shape
//      ys_sample_shape(cdf, random numbers, out element id and uv)
// 5. interpolate linearly vertex data
//    ys_interpolate_vertex(element data, vertex data, out vertex value)
// 6. [support] make a dictionary of unique undirected edges from elements
//    edge_map = ys_make_edge_map(elements)
//    ys_get_edge_index(edge_map, directed edge data)
//    ys_free_edge_map(edge_map);
//
// The interface for each function is described in details in the interface
// section of this file.
//
// Shapes are indexed meshes and are described by their
// number of elements, an array of vertex indices,
// the primitive type (points, lines, triangles, quads),
// an array of arbitrary vertex data. We differ from other library since
// vertex data is always arbitrary and we require only vertex position
// in most functions.
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

#ifndef _YS_H_
#define _YS_H_

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
// Element types
//
enum {
    ys_etype_point = 1,     // points
    ys_etype_line = 2,      // lines
    ys_etype_triangle = 3,  // triangles
};

//
// Test shapes (this is mostly used to create tests)
//
enum {
    ys_stype_uvsphere = 0,
    ys_stype_uvquad,
    ys_stype_uvcube,
    ys_stype_uvflippedsphere,
    ys_stype_uvspherecube,
    ys_stype_uvspherizedcube,
    ys_stype_uvflipcapsphere,
    ys_stype_points,
    ys_stype_rnpoints,
    ys_stype_lines,
    ys_stype_rnlines,
    ys_stype_max
};

//
// Compute smoothed normals or tangents (for lines).
//
// Parameters:
// - nelems: number of elements
// - elem: element array
// - etype: element type
// - nverts: number of vertices
// - pos: array pf vertex positions
// - weighted: whether to use area weighting (typically true)
//
// Out Parameters:
// - norm: array of computed normals (allocated internally if NULL)
//
// Returns:
// - array of vertex normals (either norm or one allocated internally)
//
YGL_API ym_vec3f* ys_compute_normals(int nelems, const int* elem, int etype,
                                     int nverts, const ym_vec3f* pos,
                                     ym_vec3f* norm, bool weighted);

//
// Tesselates a mesh by subdiving along element edges. Will produce a new
// set of elements referrring to new vertex indices in the range [0,nverts),
// for the original vertices and [nverts,nverts+nedges) for vertices in the
// split edges. The edges are returned so new vertices can be created.
// For a simpler interface, see ys_tesselate_stdshape that handles everything
// internally.
//
// Parameters:
// - nelems: number of elements
// - elem: element array
// - etype: element type
// - nverts: number of vertices
//
// Out Parameters:
// - tess_nelems: number of elements after tesselation
// - tess_elem: tesselated elements (preallocated with max nelems*4*etype
// elements)
// - nedges: number of edges
// -
//
YGL_API void ys_tesselate_shape(int nelems, const int* elem, int etype,
                                int nverts, int* tess_nelems, int* tess_elem,
                                int* nedges, ym_vec2i* edges);

//
// Tesselate a shape inplace.
//
// In/Out Parameters:
// - nelems: number of elements
// - elem: element array
// - etype: element type
// - nverts: number of vertices
// - pos: vertex pos
// - norm: vertex norm
// - texcoord: vertex texcoord
// - color: vertex color
//
YGL_API void ys_tesselate_stdshape(int* nelems, ym_vector<int>* elem, int etype,
                                   int* nverts, ym_vector<ym_vec3f>* pos,
                                   ym_vector<ym_vec3f>* norm,
                                   ym_vector<ym_vec2f>* texcoord,
                                   ym_vector<ym_vec3f>* color,
                                   ym_vector<float>* radius);

//
// Makes a parametric uv grid in [0-1]x[0-1] helpful to generate
// parametric surfaces.
//
// Parameters:
// - usteps: subdivisions in u
// - vsteps: subdivisions in v
// - etype: requested element type
//
// Out Parameters:
// - nelems: number of elements
// - elem: element array
// - nverts: number of vertices
// - uv: array of vertex uv
//
YGL_API void ys_make_uvgrid(int usteps, int vsteps, int etype, int* elem,
                            ym_vec2f* uv);

//
// Gets the size of a parametric uv grid.
//
// Parameters:
// - usteps: subdivisions in u
// - vsteps: subdivisions in v
// - etype: requested element type
//
YGL_API void ys_get_uvgrid_size(int usteps, int vsteps, int etype, int* nelems,
                                int* nverts);

//
// Computes the distribution of area of a shape element for sampling. This is
// needed to sample a shape.
//
// Paramaters:
// - nelems: number of elements
// - elem: element array
// - etype: element type
// - pos: vertex positions
//
// Out Parameters:
// - ecdf: array of element cdfs (or NULL if allocated internally)
// - area: total area of the shape (or length for lines)
//
YGL_API void ys_sample_shape_cdf(int nelems, const int* elem, int etype,
                                 const ym_vec3f* pos, float* ecdf, float* area);

//
// Sampels a shape element.
//
// Paramaters:
// - nelems: number of elements
// - ecdf: element cdf from the above function
// - etype: element type
// - ern, uvrn: random numbers, int [0,1) range, for element and uv choices
//
// Out Parameters:
// - eid: element index
// - euv: element baricentric coordinates
//
YGL_API void ys_sample_shape(int nelems, const float* ecdf, int etype,
                             float ern, const ym_vec2f& uvrn, int* eid,
                             ym_vec2f* euv);

//
// Baricentric interpolation of vertex values.
//
// Parameters:
// - elem: element array
// - etype: element type
// - vert: vertex data
// - eid: element id
// - euv: element baricentric uv
// - vsize: number of float of the per-vertex data
//
// Out Parameters:
// - v: interpolated value
//
template <typename T>
YGL_API T ys_interpolate_vertex(const int* elem, int etype, const T* vert,
                                int eid, const ym_vec2f& euv);

//
// Create standard shapes for testing purposes.
//
// Parameters:
// - stype: shape type
// - level: tesselation level (roughly euivalent to creating 2^level splits)
// - etype: element type
// - params: shape params for sample shapes (see code for documentation)
//
// Out Parameters:
// - nelems: number of elements
// - nverts: number of vertices
// - elem: element array
// - pos: vertex positions
// - norm: vertex normals
// - texcoord: vertex texture coordinates
// - radius: vertex radius
//
YGL_API void ys_make_stdshape(int stype, int level, int etype,
                              const ym_vec4f& params, int* nelems,
                              ym_vector<int>* elem, int* nverts,
                              ym_vector<ym_vec3f>* pos,
                              ym_vector<ym_vec3f>* norm,
                              ym_vector<ym_vec2f>* texcoord,
                              ym_vector<float>* radius);

//
// Dictionary from directed edges to undirected edges implemented as a hashmap
//
struct ys_edge_map;

//
// Build an edge map
//
YGL_API ys_edge_map* ys_make_edge_map(int nelems, const int* elem, int etype);

//
// Get a unique edge index from the indices defining a directed edge
//
YGL_API int ys_get_edge_index(const ys_edge_map* em, int v0, int v1);

//
// Get the number of edges
//
YGL_API int ys_get_nedges(const ys_edge_map* em);

//
// Get the number of edges
//
YGL_API const ym_vec2i* ys_get_edges(const ys_edge_map* em);

//
// Free edge map memory
//
YGL_API void ys_free_edge_map(ys_edge_map* em);

#endif

// -----------------------------------------------------------------------------
// C/C++ INTERFACE
// -----------------------------------------------------------------------------

//
// Compute smoothed normals or tangents (for lines).
//
YGLC_API float* ysc_compute_normals(int nelems, const int* elem, int etype,
                                    int nverts, const float* pos, float* norm,
                                    bool weighted);

//
// Tesselates a mesh by subdiving along element edges.
//
YGLC_API void ysc_tesselate_shape(int nelems, const int* elem, int etype,
                                  int nverts, int* tess_nelems, int* tess_elem,
                                  int* nedges, int* edges);

//
// Makes a parametric uv grid in [0-1]x[0-1] helpful to generate
// parametric surfaces.
//
YGLC_API void ysc_make_uvgrid(int usteps, int vsteps, int etype, int* elem,
                              float* uv);

//
// Gets the size of a parametric uv grid.
//
YGLC_API void ysc_get_uvgrid_size(int usteps, int vsteps, int etype,
                                  int* nelems, int* nverts);

//
// Computes the distribution of area of a shape element for sampling. This is
// needed to sample a shape.
//
YGLC_API void ysc_sample_shape_cdf(int nelems, const int* elem, int etype,
                                   const float* pos, float* ecdf, float* area);

//
// Sampels a shape element.
//
YGLC_API void ysc_sample_shape(int nelems, const float* ecdf, int etype,
                               float ern, const float uvrn[2], int* eid,
                               float euv[2]);

//
// Baricentric interpolation of vertex values.
//
YGLC_API void ysc_interpolate_vertex(const int* elem, int etype,
                                     const float* vert, int eid,
                                     const float euv[2], int vsize, float* v);

//
// Dictionary from directed edges to undirected edges implemented as a hashmap
//
struct ys_edge_map {
    ym_vector<ym_vector<ym_vec3i>> hedges;  // map buckets of directed edges
    ym_vector<ym_vec2i> edges;              // array of undirected edges
};

//
// Build an edge map
//
YGLC_API ys_edge_map* ysc_make_edge_map(int nelems, const int* elem, int etype);

//
// Get the number of edges
//
YGLC_API int ysc_get_nedges(const ys_edge_map* em);

//
// Get the number of edges
//
YGLC_API const int* ysc_get_edges(const ys_edge_map* em);

//
// Get a unique edge index from the indices defining a directed edge
//
YGLC_API int ysc_get_edge_index(const ys_edge_map* em, int v0, int v1);

//
// Free edge map memory
//
YGLC_API void ysc_free_edge_map(ys_edge_map* em);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION)

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// NORMAL COMPUTATION
// -----------------------------------------------------------------------------

// line tangent
static inline ym_vec3f ys__compute_line_tangent(const ym_vec3f& v0,
                                                const ym_vec3f& v1,
                                                bool normalize) {
    ym_vec3f n = v1 - v0;
    if (normalize) n = ym_normalize(n);
    return n;
}

// triangle tangent
static inline ym_vec3f ys__compute_triangle_normal(const ym_vec3f& v0,
                                                   const ym_vec3f& v1,
                                                   const ym_vec3f& v2,
                                                   bool normalize) {
    ym_vec3f e1 = v1 - v0;
    ym_vec3f e2 = v2 - v0;
    ym_vec3f n = ym_cross(e1, e2);
    if (normalize) n = ym_normalize(n);
    return n;
}

//
// Normal computation. Public API described above.
//
// Implementation Notes:
// - computes vertex normals as weighted averages over the face normals.
// - quads are treated as two triangles.
//
YGL_API ym_vec3f* ys_compute_normals(int nelems, const int* elem, int etype,
                                     int nverts, const ym_vec3f* pos,
                                     ym_vec3f* norm, bool weighted) {
    // clear normals
    for (int i = 0; i < nverts; i++) norm[i] = ym_zero3f;

    // handle various primitives
    switch (etype) {
        case ys_etype_point: {
            for (int i = 0; i < nverts; i++) {
                norm[i] = ym_vec3f{0, 0, 1};
            }
        } break;
        case ys_etype_line: {
            ym_vec2i* lines = (ym_vec2i*)elem;
            for (int i = 0; i < nelems; i++) {
                const ym_vec2i f = lines[i];
                ym_vec3f n =
                    ys__compute_line_tangent(pos[f.x], pos[f.y], !weighted);
                norm[f.x] += n;
                norm[f.y] += n;
            }
        } break;
        case ys_etype_triangle: {
            ym_vec3i* triangles = (ym_vec3i*)elem;
            for (int i = 0; i < nelems; i++) {
                const ym_vec3i f = triangles[i];
                ym_vec3f n = ys__compute_triangle_normal(pos[f.x], pos[f.y],
                                                         pos[f.z], !weighted);
                norm[f.x] += n;
                norm[f.y] += n;
                norm[f.z] += n;
            }
        } break;
        default: { assert(false); } break;
    }

    // normalize result
    for (int i = 0; i < nverts; i++) norm[i] = ym_normalize(norm[i]);

    // done
    return norm;
}

//
// Shape tesselation. Public API described above.
//
YGL_API void ys_tesselate_shape(int nelems, const int* elem, int etype,
                                int nverts, int* tess_nelems, int* tess_elem,
                                int* nedges, ym_vec2i* edges) {
    switch (etype) {
        case ys_etype_point: {
            *tess_nelems = nelems;
            *nedges = 0;
            for (int i = 0; i < nelems; i++) tess_elem[i] = elem[i];
        } break;
        case ys_etype_line: {
            const ym_vec2i* lines = (const ym_vec2i*)elem;
            ym_vec2i* tess_lines = (ym_vec2i*)tess_elem;
            *nedges = nelems;
            for (int i = 0; i < nelems; i++) edges[i] = lines[i];
            *tess_nelems = 2 * nelems;
            for (int i = 0; i < nelems; i++) {
                tess_lines[i * 2 + 0] = {lines[i].x, nverts + i};
                tess_lines[i * 2 + 1] = {nverts + i, lines[i].y};
            }
        };
        case ys_etype_triangle: {
            const ym_vec3i* triangles = (const ym_vec3i*)elem;
            ym_vec3i* tess_triangles = (ym_vec3i*)tess_elem;
            ys_edge_map* edge_map = ys_make_edge_map(nelems, elem, etype);
            *nedges = ys_get_nedges(edge_map);
            for (int i = 0; i < *nedges; i++) {
                edges[i] = ys_get_edges(edge_map)[i];
            }
            *tess_nelems = 4 * nelems;
            for (int i = 0; i < nelems; i++) {
                ym_vec3i f = triangles[i];
                ym_vec3i e = {nverts + ys_get_edge_index(edge_map, f.x, f.y),
                              nverts + ys_get_edge_index(edge_map, f.y, f.z),
                              nverts + ys_get_edge_index(edge_map, f.z, f.x)};
                tess_triangles[i * 4 + 0] = {f.x, e.x, e.z};
                tess_triangles[i * 4 + 1] = {f.y, e.y, e.x};
                tess_triangles[i * 4 + 2] = {f.z, e.z, e.y};
                tess_triangles[i * 4 + 3] = {e.x, e.y, e.z};
            }
            ys_free_edge_map(edge_map);
        } break;
        default: { assert(false); } break;
    }
}

//
// Tesselate a shape inplace.
//
YGL_API void ys_tesselate_stdshape(int* nelems, ym_vector<int>* elem, int etype,
                                   int* nverts, ym_vector<ym_vec3f>* pos,
                                   ym_vector<ym_vec3f>* norm,
                                   ym_vector<ym_vec2f>* texcoord,
                                   ym_vector<ym_vec3f>* color,
                                   ym_vector<float>* radius) {
    // prepare edges and elements
    int tess_nelems, tess_nedges, tess_nverts;
    ym_vector<int> tess_elem(*nelems * 4 * etype);
    ym_vector<ym_vec2i> tess_edges(*nelems * 4);
    ys_tesselate_shape(*nelems, elem->data(), etype, *nverts, &tess_nelems,
                       tess_elem.data(), &tess_nedges, tess_edges.data());
    tess_elem.resize(tess_nelems * etype);
    tess_edges.resize(tess_nedges);

    // allocate vertex data
    tess_nverts = *nverts + tess_nedges;
    if (pos && pos->size()) pos->resize(tess_nverts);
    if (norm && norm->size()) norm->resize(tess_nverts);
    if (texcoord && texcoord->size()) texcoord->resize(tess_nverts);
    if (color && color->size()) color->resize(tess_nverts);
    if (radius && radius->size()) radius->resize(tess_nverts);

    // interpolate vertex data
    for (int j = 0; j < tess_nedges; j++) {
        ym_vec2i e = tess_edges[j];
        int vid = j + *nverts;
        if (pos && pos->size())
            pos->at(vid) = (pos->at(e.x) + pos->at(e.y)) / 2;
        if (norm && norm->size())
            norm->at(vid) = (norm->at(e.x) + norm->at(e.y)) / 2;
        if (texcoord && texcoord->size())
            texcoord->at(vid) = (texcoord->at(e.x) + texcoord->at(e.y)) / 2;
        if (color && color->size())
            color->at(vid) = (color->at(e.x) + color->at(e.y)) / 2;
        if (radius && radius->size())
            radius->at(vid) = (radius->at(e.x) + radius->at(e.y)) / 2;
    }

    // copy back
    *nverts = tess_nverts;
    *nelems = tess_nelems;
    *elem = tess_elem;
}

//
// Uv grid. Public API described above.
//
YGL_API void ys_make_uvgrid(int usteps, int vsteps, int etype, int* elem,
                            ym_vec2f* uv) {
#define __vid(i, j) ((j) * (usteps + 1) + (i))
    for (int j = 0; j <= vsteps; j++) {
        for (int i = 0; i <= usteps; i++) {
            (*uv)[__vid(i, j) * 2 + 0] = i / (float)usteps;
            (*uv)[__vid(i, j) * 2 + 1] = j / (float)vsteps;
        }
    }

    switch (etype) {
        case ys_etype_point: {
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    elem[j * usteps + i] = __vid(i, j);
                }
            }
        } break;
        case ys_etype_line: {
            for (int j = 0; j <= vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    ((ym_vec2i*)elem)[j * usteps + i] = {__vid(i, j),
                                                         __vid(i + 1, j)};
                }
            }
        } break;
        case ys_etype_triangle: {
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    ym_vec3i& f1 = ((ym_vec3i*)elem)[(j * usteps + i) * 2];
                    ym_vec3i& f2 = ((ym_vec3i*)elem)[(j * usteps + i) * 2 + 1];
                    if ((i + j) % 2) {
                        f1[0] = __vid(i, j);
                        f1[1] = __vid(i + 1, j);
                        f1[2] = __vid(i + 1, j + 1);
                        f2[0] = __vid(i + 1, j + 1);
                        f2[1] = __vid(i, j + 1);
                        f2[2] = __vid(i, j);
                    } else {
                        f1[0] = __vid(i, j);
                        f1[1] = __vid(i + 1, j);
                        f1[2] = __vid(i, j + 1);
                        f2[0] = __vid(i + 1, j + 1);
                        f2[1] = __vid(i, j + 1);
                        f2[2] = __vid(i + 1, j);
                    }
                }
            }
        } break;
        default: { assert(false); } break;
    }
#undef __vid
}

//
// Uv grid. Public API described above.
//
YGL_API void ys_get_uvgrid_size(int usteps, int vsteps, int etype, int* nelems,
                                int* nverts) {
    *nverts = (usteps + 1) * (vsteps + 1);
    switch (etype) {
        case ys_etype_point: {
            *nelems = usteps * vsteps;
        } break;
        case ys_etype_line: {
            *nelems = usteps * (vsteps + 1);
        } break;
        case ys_etype_triangle: {
            *nelems = usteps * vsteps * 2;
        } break;
        default: { assert(false); } break;
    }
}

//
// Shape tesselation. Public API described above.
//
template <typename T>
YGL_API T ys_interpolate_vertex(const int* elem, int etype, const T* vert,
                                int eid, const ym_vec2f& euv) {
    switch (etype) {
        case ys_etype_point: {
            const int* f = elem + eid;
            return vert[f[0]];
        } break;
        case ys_etype_line: {
            const int* f = elem + eid * 2;
            return vert[f[0]] * (1 - euv[0]) + vert[f[1]] * euv[0];
        } break;
        case ys_etype_triangle: {
            const int* f = elem + eid * 3;
            return vert[f[0]] * (1 - euv[0] - euv[1]) + vert[f[1]] * euv[0] +
                   vert[f[2]] * euv[1];
        } break;
        default: {
            assert(false);
            return T();
            break;
        }
    }

    return T();
}

// finds the first array element smaller than the given one
// http://stackoverflow.com/questions/6553970/
// find-the-first-element-in-an-array-that-is-greater-than-the-target
static inline int ys__bsearch_smaller(float x, const float a[], int n) {
    int low = 0, high = n;
    while (low != high) {
        int mid = (low + high) / 2;
        if (a[mid] < x)
            low = mid + 1;
        else
            high = mid;
    }
    return low;
}

//
// Sample cdf. Public API described above.
//
YGL_API void ys_sample_shape_cdf(int nelems, const int* elem, int etype,
                                 const ym_vec3f* pos, float* cdf, float* area) {
    switch (etype) {
        case ys_etype_point: {
            for (int i = 0; i < nelems; i++) {
                cdf[i] = 1;
            }
        } break;
        case ys_etype_line: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 2;
                cdf[i] = ym_length(pos[f[0]] - pos[f[1]]);
            }
        } break;
        case ys_etype_triangle: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 3;
                cdf[i] = ym_length(ym_cross(pos[f[0]] - pos[f[1]],
                                            pos[f[0]] - pos[f[2]])) /
                         2;
            }
        } break;
        default: { assert(false); }
    }
    for (int i = 1; i < nelems; i++) cdf[i] += cdf[i - 1];
    if (area) *area = cdf[nelems - 1];
    for (int i = 0; i < nelems; i++) cdf[i] /= cdf[nelems - 1];
}

//
// Sample element. Public API described above.
//
YGL_API void ys_sample_shape(int nelems, const float* cdf, int etype, float ern,
                             const ym_vec2f& uvrn, int* es, ym_vec2f* euv) {
    *es = ys__bsearch_smaller(ern, cdf, nelems);
    if (etype == ys_etype_triangle) {
        *euv = {1 - sqrtf(uvrn[0]), uvrn[1] * sqrtf(uvrn[0])};
    } else {
        *euv = {uvrn[0], uvrn[1]};
    }
}

//
// Make standard shape. Public API described above.
//
YGL_API void ys_make_stdshape(int stype, int level, int etype,
                              const ym_vec4f& params, int* nelems_,
                              ym_vector<int>* elem_, int* nverts_,
                              ym_vector<ym_vec3f>* pos_,
                              ym_vector<ym_vec3f>* norm_,
                              ym_vector<ym_vec2f>* texcoord_,
                              ym_vector<float>* radius_) {
    static const ym_vec3f px = {1, 0, 0}, py = {0, 1, 0}, pz = {0, 0, 1},
                          nx = {-1, 0, 0}, ny = {0, -1, 0}, nz = {0, 0, -1};

    int& nverts = *nverts_;
    int& nelems = *nelems_;
    ym_vector<int>& elem = *elem_;
    ym_vector<ym_vec3f>& pos = *pos_;
    ym_vector<ym_vec3f>& norm = *norm_;
    ym_vector<ym_vec2f>& texcoord = *texcoord_;
    ym_vector<float>& radius = *radius_;

#define __pow2(x) (1 << (x))
#define __make_shape()                                                         \
    {                                                                          \
        elem.resize(nelems* etype);                                            \
        pos.resize(nverts);                                                    \
        norm.resize(nverts);                                                   \
        texcoord.resize(nverts);                                               \
    }

    switch (stype) {
        case ys_stype_uvsphere:
        case ys_stype_uvflippedsphere: {
            int usteps = __pow2(level + 2), vsteps = __pow2(level + 1);
            ys_get_uvgrid_size(usteps, vsteps, etype, &nelems, &nverts);
            __make_shape();
            ys_make_uvgrid(usteps, vsteps, etype, elem.data(), texcoord.data());
            for (int i = 0; i < nverts; i++) {
                ym_vec2f uv = texcoord[i];
                if (stype == ys_stype_uvsphere) {
                    ym_vec2f a = {2 * ym_pif * uv.x, ym_pif * (1 - uv.y)};
                    pos[i] = ym_vec3f{cosf(a.x) * sinf(a.y),
                                      sinf(a.x) * sinf(a.y), cosf(a.y)};
                    norm[i] = pos[i];
                } else if (stype == ys_stype_uvflippedsphere) {
                    ym_vec2f a = {2 * ym_pif * uv.x, ym_pif * uv.y};
                    pos[i] = ym_vec3f{cosf(a.x) * sinf(a.y),
                                      sinf(a.x) * sinf(a.y), cosf(a.y)};
                    norm[i] = ym_vec3f{-pos[i].x, -pos[i].y, -pos[i].z};
                    texcoord[i].y = 1 - texcoord[i].y;
                } else {
                    assert(false);
                }
            }
        } break;
        case ys_stype_uvquad: {
            int usteps = __pow2(level), vsteps = __pow2(level);
            ys_get_uvgrid_size(usteps, vsteps, etype, &nelems, &nverts);
            __make_shape();
            ys_make_uvgrid(usteps, vsteps, etype, elem.data(), texcoord.data());
            for (int i = 0; i < nverts; i++) {
                ym_vec2f uv = texcoord[i];
                pos[i] = ym_vec3f{-1 + uv.x * 2, -1 + uv.y * 2, 0};
                norm[i] = ym_vec3f{0, 0, 1};
            }
        } break;
        case ys_stype_uvcube: {
            int usteps = __pow2(level), vsteps = __pow2(level);
            int grid_nelems, grid_nverts;
            ys_get_uvgrid_size(usteps, vsteps, etype, &grid_nelems,
                               &grid_nverts);
            ym_vector<int> grid_elem(grid_nelems * etype);
            ym_vector<ym_vec2f> grid_uv(grid_nverts);
            ys_make_uvgrid(usteps, vsteps, etype, grid_elem.data(),
                           grid_uv.data());
            nelems = grid_nelems * 6;
            nverts = grid_nverts * 6;
            __make_shape();
            ym_vec3f frames[6][3] = {{px, py, pz}, {nx, py, nz}, {nx, pz, py},
                                     {px, pz, ny}, {py, pz, px}, {ny, pz, nx}};
            for (int i = 0; i < 6; i++) {
                int eoffset = i * grid_nelems * etype;
                int voffset = i * grid_nverts;
                ym_vec3f* frame = frames[i];
                for (int j = 0; j < grid_nelems * etype; j++) {
                    elem[eoffset + j] = grid_elem[j] + voffset;
                }
                for (int j = 0; j < grid_nverts; j++) {
                    pos[voffset + j] = frame[2] +
                                       frame[0] * (-1 + 2 * grid_uv[j].x) +
                                       frame[1] * (-1 + 2 * grid_uv[j].y);
                    norm[voffset + j] = frame[2];
                    texcoord[voffset + j] = grid_uv[j];
                }
            }
        } break;
        case ys_stype_uvspherecube: {
            ys_make_stdshape(ys_stype_uvcube, level, etype, ym_zero4f, &nelems,
                             &elem, &nverts, &pos, &norm, &texcoord, &radius);
            for (int i = 0; i < nverts; i++) {
                pos[i] = ym_normalize(pos[i]);
                norm[i] = pos[i];
            }
        } break;
        case ys_stype_uvspherizedcube: {
            ys_make_stdshape(ys_stype_uvcube, level, etype, ym_zero4f, &nelems,
                             &elem, &nverts, &pos, &norm, &texcoord, &radius);
            if (params[0] == 0) return;
            for (int i = 0; i < nverts; i++) {
                norm[i] = ym_normalize(pos[i]);
                pos[i] *= 1 - params[0];
                pos[i] += norm[i] * params[0];
            }
            ys_compute_normals(nelems, elem.data(), etype, nverts, pos.data(),
                               norm.data(), true);
        } break;
        case ys_stype_uvflipcapsphere: {
            ys_make_stdshape(ys_stype_uvsphere, level, etype, ym_zero4f,
                             &nelems, &elem, &nverts, &pos, &norm, &texcoord,
                             &radius);
            if (params[0] == 1) return;
            for (int i = 0; i < nverts; i++) {
                if (pos[i].z > params[0]) {
                    pos[i].z = 2 * params[0] - pos[i].z;
                    norm[i].x = -norm[i].x;
                    norm[i].y = -norm[i].y;
                } else if (pos[i].z < -params[0]) {
                    pos[i].z = -2 * params[0] - pos[i].z;
                    norm[i].x = -norm[i].x;
                    norm[i].y = -norm[i].y;
                }
            }
        } break;
        case ys_stype_points: {
            assert(etype == ys_etype_point);
            int npoints = powf(2, level);
            nverts = npoints;
            nelems = npoints;
            __make_shape();
            if (radius_) radius.resize(nverts);
            for (int i = 0; i < npoints; i++) {
                pos[i] = ym_vec3f{0, 0, 0};
                norm[i] = ym_vec3f{0, 0, 1};
                texcoord[i] = ym_vec2f{0.5f, 0.5f};
                if (radius_) radius[i] = params[0];
            }
            for (int i = 0; i < npoints; i++) elem[i] = i;
        } break;
        case ys_stype_rnpoints: {
            ys_make_stdshape(ys_stype_points, level, etype, ym_zero4f, &nelems,
                             &elem, &nverts, &pos, &norm, &texcoord, &radius);
            ym_rng_pcg32 rn;
            for (int i = 0; i < nverts; i++) {
                pos[i] = ym_vec3f{-1 + 2 * ym_rng_nextf(&rn),
                                  -1 + 2 * ym_rng_nextf(&rn),
                                  -1 + 2 * ym_rng_nextf(&rn)};
                norm[i] = ym_normalize(ym_vec3f{-1 + 2 * ym_rng_nextf(&rn),
                                                -1 + 2 * ym_rng_nextf(&rn),
                                                -1 + 2 * ym_rng_nextf(&rn)});
                texcoord[i] = ym_vec2f{ym_rng_nextf(&rn), ym_rng_nextf(&rn)};
                if (radius_) {
                    (*radius_)[i] =
                        params[0] + (params[1] - params[0]) * ym_rng_nextf(&rn);
                }
            }
        } break;
        default: {
            assert(false);
            break;
        }
    }
}

//
// Make edge map. Public API described above.
//
YGL_API ys_edge_map* ys_make_edge_map(int nelems, const int* elem, int etype) {
    assert(etype >= 2 && etype <= 4);
    ys_edge_map* map = new ys_edge_map();

    // find largest vert index
    int nverts = 0;
    for (int i = 0; i < nelems; i++) {
        const int* e = elem + i * etype;
        for (int j = 0; j < etype; j++) {
            if (nverts < e[j] + 1) nverts = e[j] + 1;
        }
    }

    // count the number of edges per vert
    // insert all half edges and set edge counts
    map->hedges.resize(nverts);
    for (int i = 0; i < nelems; i++) {
        const int* e = elem + i * etype;
        for (int j = 0; j < etype; j++) {
            int e1 = e[j], e2 = e[(j + 1) % etype];
            if (e1 == e2) continue;
            map->hedges[e1].push_back({e1, e2, 0});
            int pos = -1;
            for (int k = 0; k < map->hedges[e2].size() && pos < 0; k++) {
                if (map->hedges[e2][k].y == e1) pos = k;
            }
            if (pos < 0) {
                map->edges.push_back({ym_min(e1, e2), ym_max(e1, e2)});
                map->hedges[e1].back().z = (int)map->edges.size() - 1;
            } else {
                map->hedges[e1].back().z = map->hedges[e2][pos].z;
            }
        }
    }

    // done
    return map;
}

//
// Edge array length
//
YGL_API int ys_get_nedges(const ys_edge_map* em) {
    return (int)em->edges.size();
}

//
// Edge array length
//
YGL_API const ym_vec2i* ys_get_edges(const ys_edge_map* em) {
    return em->edges.data();
}

//
// Get edge index. Public API described above.
//
YGL_API int ys_get_edge_index(const ys_edge_map* em, int v0, int v1) {
    for (int i = 0; i < em->hedges[v0].size(); i++) {
        if (em->hedges[v0][i].y == v1) return em->hedges[v0][i].z;
    }
    return -1;
}

//
// Free edge map. Public API described above.
//
YGL_API void ys_free_edge_map(ys_edge_map* em) { delete em; }

// -----------------------------------------------------------------------------
// C API IMPLEMENTATION
// -----------------------------------------------------------------------------

//
// Compute smoothed normals or tangents (for lines).
//
YGLC_API float* ysc_compute_normals(int nelems, const int* elem, int etype,
                                    int nverts, const float* pos, float* norm,
                                    bool weighted) {
    return (float*)ys_compute_normals(nelems, elem, etype, nverts,
                                      (const ym_vec3f*)pos, (ym_vec3f*)norm,
                                      weighted);
}

//
// Tesselates a mesh by subdiving along element edges.
//
YGLC_API void ysc_tesselate_shape(int nelems, const int* elem, int etype,
                                  int nverts, int* tess_nelems, int* tess_elem,
                                  int* nedges, int* edges) {
    ys_tesselate_shape(nelems, elem, etype, nverts, tess_nelems, tess_elem,
                       nedges, (ym_vec2i*)edges);
}

//
// Makes a parametric uv grid in [0-1]x[0-1] helpful to generate
// parametric surfaces.
//
YGLC_API void ysc_make_uvgrid(int usteps, int vsteps, int etype, int* elem,
                              float* uv) {
    ys_make_uvgrid(usteps, vsteps, etype, elem, (ym_vec2f*)uv);
}

//
// Gets the size of a parametric uv grid.
//
YGLC_API void ysc_get_uvgrid_size(int usteps, int vsteps, int etype,
                                  int* nelems, int* nverts) {
    ys_get_uvgrid_size(usteps, vsteps, etype, nelems, nverts);
}

//
// Computes the distribution of area of a shape element for sampling. This is
// needed to sample a shape.
//
YGLC_API void ysc_sample_shape_cdf(int nelems, const int* elem, int etype,
                                   const float* pos, float* ecdf, float* area) {
    ys_sample_shape_cdf(nelems, elem, etype, (ym_vec3f*)pos, ecdf, area);
}

//
// Sampels a shape element.
//
YGLC_API void ysc_sample_shape(int nelems, const float* ecdf, int etype,
                               float ern, const float uvrn[2], int* eid,
                               float euv[2]) {
    ys_sample_shape(nelems, ecdf, etype, ern, ym_vec2f(uvrn), eid,
                    (ym_vec2f*)euv);
}

//
// Baricentric interpolation of vertex values.
//
YGLC_API void ysc_interpolate_vertex(const int* elem, int etype,
                                     const float* vert, int eid,
                                     const float euv[2], int vsize, float* v) {
    switch (vsize) {
        case 1: {
            v[0] = ys_interpolate_vertex(elem, etype, vert, eid, ym_vec2f(euv));
        } break;
        case 2: {
            *(ym_vec2f*)v = ys_interpolate_vertex(
                elem, etype, (const ym_vec2f*)vert, eid, ym_vec2f(euv));
        } break;
        case 3: {
            *(ym_vec3f*)v = ys_interpolate_vertex(
                elem, etype, (const ym_vec3f*)vert, eid, ym_vec2f(euv));
        } break;
        case 4: {
            *(ym_vec4f*)v = ys_interpolate_vertex(
                elem, etype, (const ym_vec4f*)vert, eid, ym_vec2f(euv));
        } break;
        default: {
            assert(false);
            break;
        }
    }
}

//
// Build an edge map
//
YGLC_API ys_edge_map* ysc_make_edge_map(int nelems, const int* elem,
                                        int etype) {
    return ys_make_edge_map(nelems, elem, etype);
}

//
// Get the number of edges
//
YGLC_API int ysc_get_nedges(const ys_edge_map* em) { return ys_get_nedges(em); }

//
// Get the number of edges
//
YGLC_API const int* ysc_get_edges(const ys_edge_map* em) {
    return (const int*)ys_get_edges(em);
}

//
// Get a unique edge index from the indices defining a directed edge
//
YGLC_API int ysc_get_edge_index(const ys_edge_map* em, int v0, int v1) {
    return ys_get_edge_index(em, v0, v1);
}

//
// Free edge map memory
//
YGLC_API void ysc_free_edge_map(ys_edge_map* em) {
    return ys_free_edge_map(em);
}

#endif

#endif
