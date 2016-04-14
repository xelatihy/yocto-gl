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
//    4.a sample num points uniformly over its surface
//      ys_sample_shape(shape data, num points, out vertex data);
//    4.b. sample shapes one point at a time with
//      - first, compute shape cdf
//      ys_sample_shape_cdf(shape definition, out cdf)
//      - generate element id and element uv for each shape
//      ys_sample_shape_elem(cdf, random numbers, out element id and uv)
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
// Quad meshes are experimental and might go away in future realeases. If
// you can, please use triangles. Quads are treated as two triangles (v0,v1,v3)
// and (v2,v3,v1). Quads with v2 == v3 are degenerate and represent one
// triangle, thus quads meshes can also represent mixtures of triangle and
// quads. This follows Intel's Embree API.
//

//
// COMPILATION:
//
// All functions in this library are inlined by default for ease of use.
// To use the library as a .h/.c pair do the following:
// - to use as a .h, just #define YS_NOINLINE before including this file
// - to use as a .c, just #define YS_IMPLEMENTATION before including this file
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

#ifndef YS_NOINLINE
#define YS_API static inline
#else
#ifdef __cplusplus
#define YS_API extern "C"
#else
#define YS_API
#endif
#endif

#include <stdbool.h>

// -----------------------------------------------------------------------------
// INTERFACE
// -----------------------------------------------------------------------------

//
// Element types
//
enum {
    ys_etype_point = 1,     // points
    ys_etype_line = 2,      // lines
    ys_etype_triangle = 3,  // triangles
    ys_etype_quad = 4       // quads
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
YS_API float*
ys_compute_normals(int nelems, const int* elem, int etype, int nverts,
                   const float* pos, float* norm, bool weighted);

//
// Tesselates a mesh by subdiving along element edges.
//
// Parameters:
// - etype: element type
// - nvertprops: number of vertex properties
// - vsize: array of vertex property sizes
//
// In/Out Parameters:
// - nelems: number of elements
// - nverts: number of vertices
// - elem: element array
// - vert: array of pointers to vertex properties
//
YS_API void
ys_tesselate_shape(int etype, int nvertprops, const int vsize[], int* nelems,
                   int** elem, int* nverts, float** vert[]);

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
YS_API void
ys_make_uvgrid(int usteps, int vsteps, int etype, int* nelems, int** elem,
               int* nverts, float** uv);

//
// Function callback for parametric shapes. Computes the vertex value v
//
//
typedef void (*vfunc)(void* ctx, int vid, const float uv[2], float* v);

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
YS_API void
ys_make_uvshape(int usteps, int vsteps, void* ctx, int etype, int nvertprops,
                int* vsize, vfunc* vfuncs, int* nelems, int** elem, int* nverts,
                float*** vert);

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
// Returns:
// - element cdf
//
YS_API float*
ys_sample_shape_cdf(int nelems, const int* elem, int etype, const float* pos,
                    float* ecdf, float* area);

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
YS_API void
ys_sample_shape_elem(int nelems, const float* ecdf, int etype, float ern,
                     const float uvrn[2], int* eid, float euv[2]);

//
// Sampels a shape generating.
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
YS_API void
ys_sample_shape(int nelems, const int* elem, int etype, const float* pos,
                int num, int nvertprops, const int* vsize, const float** vert,
                float*** sampled);

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
YS_API void
ys_interpolate_vertex(const int* elem, int etype, const float* vert, int eid,
                      const float euv[2], int vsize, float* v);

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
YS_API void
ys_make_stdshape(int stype, int level, int etype, const float params[],
                 int* nelems, int** elem, int* nverts, float** pos,
                 float** norm, float** texcoord, float** radius);

//
// Dictionary from directed edges to undirected edges implemented as a hashmap
//
typedef struct ys_edge_map {
    int nverts, nedges;  // number of vertices and edges
    int* vert_size;      // size of the hash map buckets
    int** vert_hedge;    // map buckets of directed edges
    int* edges;          // array of undirected edges
} ys_edge_map;

//
// Build an edge map
//
YS_API ys_edge_map*
ys_make_edge_map(int nelems, int* elem, int etype);

//
// Get a unique edge index from the indices defining a directed edge
//
YS_API int
ys_get_edge_index(const ys_edge_map* em, int v0, int v1);

//
// Free edge map memory
//
YS_API void
ys_free_edge_map(ys_edge_map* em);

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#if !defined(YS_NOINLINE) || defined(YS_IMPLEMENTATION)

#include <assert.h>
#include <stdlib.h>
#include <string.h>

// -----------------------------------------------------------------------------
// MATH SUPPORT
// -----------------------------------------------------------------------------

#define ys__pif 3.14159265f

//
// 2d vectors
//
typedef struct { float x, y; } ys__vec2f;

//
// 3d vectors
//
typedef struct { float x, y, z; } ys__vec3f;

// vector substration
static inline ys__vec3f
ys__sub3f(const ys__vec3f a, const ys__vec3f b) {
    return (ys__vec3f){ a.x - b.x, a.y - b.y, a.z - b.z };
}

// vector scaling
static inline ys__vec3f
ys__smul3f(const ys__vec3f a, float b) {
    return (ys__vec3f){ a.x * b, a.y * b, a.z * b };
}

// vector dot product
static inline float
ys__dot3f(const ys__vec3f a, const ys__vec3f b) {
    return (a.x * b.x + a.y * b.y + a.z * b.z);
}

// vector length
static inline float
ys__length3f(const ys__vec3f a) {
    return sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
}

// vector normalization
static inline ys__vec3f
ys__normalize3f(const ys__vec3f a) {
    float l = sqrtf(a.x * a.x + a.y * a.y + a.z * a.z);
    if (l > 0)
        return (ys__vec3f){ a.x / l, a.y / l, a.z / l };
    else
        return a;
}

// vector cross product
static inline ys__vec3f
ys__cross3f(const ys__vec3f a, const ys__vec3f b) {
    return (ys__vec3f){ a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
                        a.x * b.y - a.y * b.x };
}

// vector sum
static inline ys__vec3f
ys__sum3f(const ys__vec3f a, const ys__vec3f b) {
    return (ys__vec3f){ a.x + b.x, a.y + b.y, a.z + b.z };
}

// -----------------------------------------------------------------------------
// NORMAL COMPUTATION
// -----------------------------------------------------------------------------

// line tangent
static inline ys__vec3f
ys__compute_line_tangent(const ys__vec3f v0, const ys__vec3f v1,
                         bool normalize) {
    ys__vec3f n = ys__sub3f(v1, v0);
    if (normalize) n = ys__normalize3f(n);
    return n;
}

// triangle tangent
static inline ys__vec3f
ys__compute_triangle_normal(const ys__vec3f v0, const ys__vec3f v1,
                            const ys__vec3f v2, bool normalize) {
    ys__vec3f e1 = ys__sub3f(v1, v0);
    ys__vec3f e2 = ys__sub3f(v2, v0);
    ys__vec3f n = ys__cross3f(e1, e2);
    if (normalize) n = ys__normalize3f(n);
    return n;
}

// quad tangent
static inline ys__vec3f
ys__compute_quad_normal(const ys__vec3f v0, const ys__vec3f v1,
                        const ys__vec3f v2, const ys__vec3f v3,
                        bool normalize) {
    ys__vec3f n1 = ys__compute_triangle_normal(v0, v1, v2, normalize);
    ys__vec3f n2 = ys__compute_triangle_normal(v0, v2, v3, normalize);
    ys__vec3f n = ys__sum3f(n1, n2);
    if (normalize) n = ys__normalize3f(n);
    return n;
}

//
// Normal computation. Public API described above.
//
// Implementation Notes:
// - computes vertex normals as weighted averages over the face normals.
// - quads are treated as two triangles.
//
YS_API float*
ys_compute_normals(int nelems, const int* elem, int etype, int nverts,
                   const float* pos_, float* norm_, bool weighted) {
    // convert to internal types
    ys__vec3f *pos = (ys__vec3f *)pos_, *norm = (ys__vec3f *)norm_;

    // allocate if not passed
    if (!norm) norm = (ys__vec3f*)calloc(nverts, sizeof(ys__vec3f));

    // clear normals
    memset(norm, 0, sizeof(ys__vec3f) * nverts);

    // handle various primitives
    switch (etype) {
        case ys_etype_point: {
            for (int i = 0; i < nverts; i++) {
                norm[i] = (ys__vec3f){ 0, 0, 1 };
            }
        } break;
        case ys_etype_line: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 2;
                ys__vec3f *v0 = pos + f[0], *v1 = pos + f[1];
                ys__vec3f n = ys__compute_line_tangent(*v0, *v1, !weighted);
                for (int j = 0; j < 2; j++) {
                    norm[f[j]] = ys__sum3f(norm[f[j]], n);
                }
            }
        } break;
        case ys_etype_triangle: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 3;
                ys__vec3f *v0 = pos + f[0], *v1 = pos + f[1], *v2 = pos + f[2];
                ys__vec3f n =
                    ys__compute_triangle_normal(*v0, *v1, *v2, !weighted);
                for (int j = 0; j < 3; j++) {
                    norm[f[j]] = ys__sum3f(norm[f[j]], n);
                }
            }
        } break;
        case ys_etype_quad: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 4;
                ys__vec3f *v0 = pos + f[0], *v1 = pos + f[1], *v2 = pos + f[2],
                          *v3 = pos + f[3];
                if (f[2] != f[3]) {
                    ys__vec3f n =
                        ys__compute_quad_normal(*v0, *v1, *v2, *v3, !weighted);
                    for (int j = 0; j < 4; j++) {
                        norm[f[j]] = ys__sum3f(norm[f[j]], n);
                    }
                } else {
                    ys__vec3f n =
                        ys__compute_triangle_normal(*v0, *v1, *v2, !weighted);
                    for (int j = 0; j < 3; j++) {
                        norm[f[j]] = ys__sum3f(norm[f[j]], n);
                    }
                }
            }
        } break;
        default: { assert(false); } break;
    }

    // normalize result
    for (int i = 0; i < nverts; i++) norm[i] = ys__normalize3f(norm[i]);

    // done
    return (float*)norm;
}

//
// Shape tesselation. Public API described above.
//
YS_API void
ys_tesselate_shape(int tess_etype, int nvertprops, const int* vsize,
                   int* tess_nelems, int** tess_elem, int* tess_nverts,
                   float*** tess_vert) {
    // save data
    int nverts = *tess_nverts, nelems = *tess_nelems, *elem = *tess_elem,
        etype = tess_etype;

    // get edges
    int nedges = 0, *edge = 0;
    ys_edge_map* edge_map = 0;
    if (etype == ys_etype_triangle || etype == ys_etype_quad) {
        edge_map = ys_make_edge_map(nelems, elem, etype);
        edge = edge_map->edges;
        nedges = edge_map->nedges;
    } else if (etype == ys_etype_line) {
        edge = elem;
        nedges = nelems;
    } else {
        assert(false);
    }

    // create face lists for quads to mark degenerate quads
    int nfacesplits = 0, *facesplit = 0;
    if (etype == ys_etype_quad) {
        // count face splits
        facesplit = (int*)calloc(nelems, sizeof(int));
        for (int i = 0; i < nelems; i++) {
            int* e = elem + i * 4;
            if (e[2] != e[3]) {
                facesplit[i] = nfacesplits;
                nfacesplits++;
            } else {
                facesplit[i] = -1;
            }
        }
    }

    // create vertices
    *tess_nverts = nverts + nedges + nfacesplits;
    // perform op for each vertex properties
    for (int vp = 0; vp < nvertprops; vp++) {
        // skip if needed
        if (!tess_vert[vp] || !(*tess_vert[vp])) continue;
        // grab vert data
        float* vert = *tess_vert[vp];
        float* tess = (float*)calloc(vsize[vp] * (*tess_nverts), sizeof(float));
        *(tess_vert[vp]) = tess;
        // copy old vertices
        memcpy(tess, vert, sizeof(float) * vsize[vp] * nverts);
        // edge vertices
        int voffset = nverts;
        for (int e = 0; e < nedges; e++) {
            int* f = edge + e * 2;
            for (int c = 0; c < vsize[vp]; c++)
                tess[(voffset + e) * vsize[vp] + c] =
                    (vert[f[0] * vsize[vp] + c] + vert[f[1] * vsize[vp] + c]) /
                    2;
        }
        // face vertices
        if (facesplit) {
            voffset += nedges;
            for (int e = 0; e < nelems; e++) {
                if (facesplit[e] < 0) continue;
                int* face = elem + e * etype;
                for (int c = 0; c < vsize[vp]; c++) {
                    float* tv = tess + (voffset + facesplit[e]) * vsize[vp] + c;
                    *tv = 0;
                    for (int ve = 0; ve < etype; ve++)
                        *tv += vert[face[ve] * vsize[vp] + c];
                    *tv /= etype;
                }
            }
        }
        // cleanup
        free(vert);
    }

    // create elems
    if (etype == ys_etype_triangle) {
        *tess_nelems = 4 * nelems;
        *tess_elem = (int*)calloc(3 * 4 * nelems, sizeof(int));
        for (int e = 0; e < nelems; e++) {
            int* f = elem + e * 3;
            int ev[3] = { nverts + ys_get_edge_index(edge_map, f[0], f[1]),
                          nverts + ys_get_edge_index(edge_map, f[1], f[2]),
                          nverts + ys_get_edge_index(edge_map, f[2], f[0]) };
            int new_face[4][3] = { { f[0], ev[0], ev[2] },
                                   { f[1], ev[1], ev[0] },
                                   { f[2], ev[2], ev[1] },
                                   { ev[0], ev[1], ev[2] } };
            memcpy(*tess_elem + e * 12, new_face, sizeof(new_face));
        }
        free(elem);
    } else if (etype == ys_etype_quad) {
        *tess_nelems = 4 * nelems;
        *tess_elem = (int*)calloc(4 * 4 * nelems, sizeof(int));
        for (int e = 0; e < nelems; e++) {
            int* f = elem + e * 4;
            if (f[2] != f[3]) {
                int ev[4] = { nverts + ys_get_edge_index(edge_map, f[0], f[1]),
                              nverts + ys_get_edge_index(edge_map, f[1], f[2]),
                              nverts + ys_get_edge_index(edge_map, f[2], f[3]),
                              nverts +
                                  ys_get_edge_index(edge_map, f[3], f[0]) };
                assert(facesplit[e] >= 0);
                int fv = nverts + nedges + facesplit[e];
                int new_face[4][4] = { { f[0], ev[0], fv, ev[3] },
                                       { f[1], ev[1], fv, ev[0] },
                                       { f[2], ev[2], fv, ev[1] },
                                       { f[3], ev[3], fv, ev[2] } };
                memcpy(*tess_elem + e * 16, new_face, sizeof(new_face));
            } else {
                int ev[3] = { nverts + ys_get_edge_index(edge_map, f[0], f[1]),
                              nverts + ys_get_edge_index(edge_map, f[1], f[2]),
                              nverts +
                                  ys_get_edge_index(edge_map, f[2], f[0]) };
                int new_face[4][4] = { { f[0], ev[0], ev[2], ev[2] },
                                       { f[1], ev[1], ev[0], ev[0] },
                                       { f[2], ev[2], ev[1], ev[1] },
                                       { ev[0], ev[1], ev[2], ev[2] } };
                memcpy(*tess_elem + e * 16, new_face, sizeof(new_face));
            }
        }
        free(elem);
    } else if (etype == ys_etype_line) {
        *tess_nelems = 2 * nelems;
        *tess_elem = (int*)calloc(2 * 2 * nelems, sizeof(int));
        for (int e = 0; e < nelems; e++) {
            int* f = elem + e * 4;
            int ev = nverts + ys_get_edge_index(edge_map, f[0], f[1]);
            int new_face[2][2] = { { f[0], ev }, { ev, f[1] } };
            memcpy(*tess_elem + e * 4, new_face, sizeof(new_face));
        }
        free(elem);
    } else {
        assert(false);
    }

    // cleanup
    if (facesplit) free(facesplit);
    ys_free_edge_map(edge_map);
}

//
// Uv grid. Public API described above.
//
YS_API void
ys_make_uvgrid(int usteps, int vsteps, int etype, int* nelems, int** elem,
               int* nverts, float** uv) {
#define __vid(i, j) ((j) * (usteps + 1) + (i))
    *nverts = (usteps + 1) * (vsteps + 1);
    *uv = (float*)calloc(*nverts * 2, sizeof(float));
    for (int j = 0; j <= vsteps; j++) {
        for (int i = 0; i <= usteps; i++) {
            (*uv)[__vid(i, j) * 2 + 0] = i / (float)usteps;
            (*uv)[__vid(i, j) * 2 + 1] = j / (float)vsteps;
        }
    }

    switch (etype) {
        case ys_etype_point: {
            *nelems = usteps * vsteps;
            *elem = (int*)calloc(*nelems, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i);
                    f[0] = __vid(i, j);
                }
            }
        } break;
        case ys_etype_line: {
            *nelems = usteps * (vsteps + 1);
            *elem = (int*)calloc(*nelems * 2, sizeof(int));
            for (int j = 0; j <= vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i) * 2;
                    f[0] = __vid(i, j);
                    f[1] = __vid(i + 1, j);
                }
            }
        } break;
        case ys_etype_triangle: {
            *nelems = usteps * vsteps * 2;
            *elem = (int*)calloc(*nelems * 3, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f1 = *elem + (j * usteps + i) * 6;
                    int* f2 = f1 + 3;
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
        case ys_etype_quad: {
            *nelems = usteps * vsteps;
            *elem = (int*)calloc(*nelems * 4, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i) * 4;
                    f[0] = __vid(i, j);
                    f[1] = __vid(i + 1, j);
                    f[2] = __vid(i + 1, j + 1);
                    f[3] = __vid(i, j + 1);
                }
            }
        } break;
        default: { assert(false); } break;
    }
#undef __vid
}

//
// Uv shape. Public API described above.
//
YS_API void
ys_make_uvshape(int usteps, int vsteps, void* ctx, int etype, int nvertprops,
                int* vsize, vfunc* vfuncs, int* nelems, int** elem, int* nverts,
                float*** vert) {
#define __vid(i, j) ((j) * (usteps + 1) + (i))
    *nverts = (usteps + 1) * (vsteps + 1);
    for (int v = 0; v < nvertprops; v++) {
        *vert[v] = (float*)calloc(*nverts * vsize[v], sizeof(float));
    }
    for (int j = 0; j <= vsteps; j++) {
        for (int i = 0; i <= usteps; i++) {
            float uv[2] = { i / (float)usteps, j / (float)vsteps };
            for (int v = 0; v < nvertprops; v++) {
                if (vfuncs[v]) {
                    vfuncs[v](ctx, __vid(i, j), uv,
                              *vert[v] + __vid(i, j) * vsize[v]);
                } else {
                    for (int vv = 0; vv < vsize[v]; vv++) {
                        (*vert[v] + __vid(i, j) * vsize[v])[vv] = uv[vv % 2];
                    }
                }
            }
        }
    }

    switch (etype) {
        case ys_etype_point: {
            *nelems = usteps * vsteps;
            *elem = (int*)calloc(*nelems, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i);
                    f[0] = __vid(i, j);
                }
            }
        } break;
        case ys_etype_line: {
            *nelems = usteps * (vsteps + 1);
            *elem = (int*)calloc(*nelems * 2, sizeof(int));
            for (int j = 0; j <= vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i) * 2;
                    f[0] = __vid(i, j);
                    f[1] = __vid(i + 1, j);
                }
            }
        } break;
        case ys_etype_triangle: {
            *nelems = usteps * vsteps * 2;
            *elem = (int*)calloc(*nelems * 3, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f1 = *elem + (j * usteps + i) * 6;
                    int* f2 = f1 + 3;
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
        case ys_etype_quad: {
            *nelems = usteps * vsteps;
            *elem = (int*)calloc(*nelems * 4, sizeof(int));
            for (int j = 0; j < vsteps; j++) {
                for (int i = 0; i < usteps; i++) {
                    int* f = *elem + (j * usteps + i) * 4;
                    f[0] = __vid(i, j);
                    f[1] = __vid(i + 1, j);
                    f[2] = __vid(i + 1, j + 1);
                    f[3] = __vid(i, j + 1);
                }
            }
        } break;
        default: { assert(false); } break;
    }
#undef __vid
}

//
// Shape tesselation. Public API described above.
//
YS_API void
ys_interpolate_vertex(const int* elem, int etype, const float* vert, int eid,
                      const float euv[2], int vsize, float* v) {
    switch (etype) {
        case ys_etype_point: {
            for (int c = 0; c < vsize; c++) {
                const int* f = elem + eid;
                v[c] = vert[f[0] * vsize + c];
            }
        } break;
        case ys_etype_line: {
            for (int c = 0; c < vsize; c++) {
                const int* f = elem + eid * 2;
                v[c] = vert[f[0] * vsize + c] * (1 - euv[0]) +
                       vert[f[1] * vsize + c] * euv[0];
            }
        } break;
        case ys_etype_triangle: {
            for (int c = 0; c < vsize; c++) {
                const int* f = elem + eid * 3;
                v[c] = vert[f[0] * vsize + c] * (1 - euv[0] - euv[1]) +
                       vert[f[1] * vsize + c] * euv[0] +
                       vert[f[2] * vsize + c] * euv[1];
            }
        } break;
        case ys_etype_quad: {
            for (int c = 0; c < vsize; c++) {
                const int* f = elem + eid * 4;
                if (euv[0] + euv[1] < 1) {
                    v[c] = vert[f[0] * vsize + c] * (1 - euv[0] - euv[1]) +
                           vert[f[1] * vsize + c] * euv[0] +
                           vert[f[3] * vsize + c] * euv[1];
                } else {
                    v[c] = vert[f[2] * vsize + c] *
                               (1 - (1 - euv[0]) - (1 - euv[1])) +
                           vert[f[3] * vsize + c] * (1 - euv[0]) +
                           vert[f[1] * vsize + c] * (1 - euv[1]);
                }
            }
        } break;
        default: break;
    }
}

// lcg random number
static inline float
ys__rng_nextf(unsigned int* state) {
    *state = (1103515245U * *state + 12345U) % 2147483648U;
    return fminf((float)*state / (float)2147483648U, 0.999999f);
}

// finds the first array element smaller than the given one
// http://stackoverflow.com/questions/6553970/
// find-the-first-element-in-an-array-that-is-greater-than-the-target
static inline int
ys__bsearch_smaller(float x, const float a[], int n) {
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
YS_API float*
ys_sample_shape_cdf(int nelems, const int* elem, int etype, const float* pos_,
                    float* cdf, float* area) {
    if (!cdf) cdf = (float*)calloc(nelems, sizeof(float));
    ys__vec3f* pos = (ys__vec3f*)pos_;
    switch (etype) {
        case ys_etype_point: {
            for (int i = 0; i < nelems; i++) {
                cdf[i] = 1;
            }
        } break;
        case ys_etype_line: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 2;
                cdf[i] = ys__length3f(ys__sub3f(pos[f[0]], pos[f[1]]));
            }
        } break;
        case ys_etype_triangle: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 3;
                cdf[i] =
                    ys__length3f(ys__cross3f(ys__sub3f(pos[f[0]], pos[f[1]]),
                                             ys__sub3f(pos[f[0]], pos[f[2]]))) /
                    2;
            }
        } break;
        case ys_etype_quad: {
            for (int i = 0; i < nelems; i++) {
                const int* f = elem + i * 4;
                cdf[i] =
                    ys__length3f(ys__cross3f(ys__sub3f(pos[f[0]], pos[f[1]]),
                                             ys__sub3f(pos[f[0]], pos[f[2]]))) /
                        2 +
                    ys__length3f(ys__cross3f(ys__sub3f(pos[f[0]], pos[f[2]]),
                                             ys__sub3f(pos[f[0]], pos[f[3]]))) /
                        2;
            }
        } break;
        default: { assert(false); }
    }
    for (int i = 1; i < nelems; i++) cdf[i] += cdf[i - 1];
    if (area) *area = cdf[nelems - 1];
    for (int i = 0; i < nelems; i++) cdf[i] /= cdf[nelems - 1];
    return cdf;
}

//
// Sample element. Public API described above.
//
YS_API void
ys_sample_shape_elem(int nelems, const float* cdf, int etype, float ern,
                     const float uvrn[2], int* es, float euv[2]) {
    *es = ys__bsearch_smaller(ern, cdf, nelems);
    if (etype == ys_etype_triangle) {
        euv[0] = 1 - sqrtf(uvrn[0]);
        euv[1] = uvrn[1] * sqrtf(uvrn[0]);
    } else {
        euv[0] = uvrn[0];
        euv[1] = uvrn[1];
    }
}

//
// Sample shape. Public API described above.
//
YS_API void
ys_sample_shape(int nelems, const int* elem, int etype, const float* pos,
                int num, int nvertprops, const int* vsize, const float** vert,
                float*** sampled) {
    float* cdf = ys_sample_shape_cdf(nelems, elem, etype, pos, 0, 0);

    for (int p = 0; p < nvertprops; p++) {
        if (!*sampled[p])
            *sampled[p] = (float*)calloc(num * vsize[p], sizeof(float));
    }

    unsigned int rng_state = 0;
    for (int i = 0; i < num; i++) {
        int eid;
        float euv[2];
        ys_sample_shape_elem(
            nelems, cdf, etype, ys__rng_nextf(&rng_state),
            (float[2]){ ys__rng_nextf(&rng_state), ys__rng_nextf(&rng_state) },
            &eid, euv);
        for (int p = 0; p < nvertprops; p++) {
            if (!*sampled[p]) {
                ys_interpolate_vertex(elem, etype, vert[p] + i * vsize[p], eid,
                                      euv, vsize[p],
                                      *sampled[p] + i * vsize[p]);
            }
        }
    }
}

//
// Make standard shape. Public API described above.
//
YS_API void
ys_make_stdshape(int stype, int level, int etype, const float* params,
                 int* nelems, int** elem, int* nverts, float** pos_,
                 float** norm_, float** texcoord_, float** radius_) {
    static const ys__vec3f px = { 1, 0, 0 }, py = { 0, 1, 0 }, pz = { 0, 0, 1 },
                           nx = { -1, 0, 0 }, ny = { 0, -1, 0 },
                           nz = { 0, 0, -1 };

#define __pow2(x) (1 << (x))

    switch (stype) {
        case ys_stype_uvsphere:
        case ys_stype_uvflippedsphere:
        case ys_stype_uvquad: {
            ys_make_uvgrid(__pow2(level + 2), __pow2(level + 1), etype, nelems,
                           elem, nverts, texcoord_);
            *pos_ = (float*)calloc(*nverts * 3, sizeof(float));
            *norm_ = (float*)calloc(*nverts * 3, sizeof(float));
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            ys__vec2f* texcoord = (ys__vec2f*)*texcoord_;
            for (int i = 0; i < *nverts; i++) {
                ys__vec2f uv = texcoord[i];
                if (stype == ys_stype_uvsphere) {
                    ys__vec2f a = { 2 * ys__pif * uv.x, ys__pif * (1 - uv.y) };
                    pos[i] = (ys__vec3f){ cosf(a.x) * sinf(a.y),
                                          sinf(a.x) * sinf(a.y), cosf(a.y) };
                    norm[i] = pos[i];
                } else if (stype == ys_stype_uvflippedsphere) {
                    ys__vec2f a = { 2 * ys__pif * uv.x, ys__pif * uv.y };
                    pos[i] = (ys__vec3f){ cosf(a.x) * sinf(a.y),
                                          sinf(a.x) * sinf(a.y), cosf(a.y) };
                    norm[i] = (ys__vec3f){ -pos[i].x, -pos[i].y, -pos[i].z };
                    texcoord[i].y = 1 - texcoord[i].y;
                } else if (stype == ys_stype_uvquad) {
                    pos[i] = (ys__vec3f){ -1 + uv.x * 2, -1 + uv.y * 2, 0 };
                    norm[i] = (ys__vec3f){ 0, 0, 1 };
                } else
                    assert(false);
            }
        } break;
        case ys_stype_uvcube: {
            int grid_nelems, *grid_elem, grid_nverts;
            ys__vec2f* grid_uv;
            ys_make_uvgrid(__pow2(level + 2), __pow2(level + 1), etype,
                           &grid_nelems, &grid_elem, &grid_nverts,
                           (float**)&grid_uv);
            *nelems = grid_nelems * 6;
            *elem = (int*)calloc(*nelems * 6 * etype, sizeof(int));
            *nverts = grid_nverts * 6;
            *pos_ = (float*)calloc(*nverts * 3, sizeof(float));
            *norm_ = (float*)calloc(*nverts * 3, sizeof(float));
            *texcoord_ = (float*)calloc(*nverts * 2, sizeof(float));
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            ys__vec2f* texcoord = (ys__vec2f*)*texcoord_;
            ys__vec3f frames[6][3] = { { px, py, pz }, { nx, py, nz },
                                       { nx, pz, py }, { px, pz, ny },
                                       { py, pz, px }, { ny, pz, nx } };
            for (int i = 0; i < 6; i++) {
                int eoffset = i * grid_nelems * etype;
                int voffset = i * grid_nverts;
                ys__vec3f* frame = frames[i];
                for (int j = 0; j < grid_nelems * etype; j++) {
                    (*elem)[eoffset + j] = grid_elem[j] + voffset;
                }
                for (int j = 0; j < grid_nverts; j++) {
                    pos[voffset + j] = frame[2];
                    pos[voffset + j] =
                        ys__sum3f(pos[voffset + j],
                                  ys__smul3f(frame[0], -1 + 2 * grid_uv[j].x));
                    pos[voffset + j] =
                        ys__sum3f(pos[voffset + j],
                                  ys__smul3f(frame[1], -1 + 2 * grid_uv[j].y));
                    norm[voffset + j] = frame[2];
                    texcoord[voffset + j] = grid_uv[j];
                }
            }
            free(grid_uv);
            free(grid_elem);
        } break;
        case ys_stype_uvspherecube: {
            ys_make_stdshape(ys_stype_uvcube, level, etype, 0, nelems, elem,
                             nverts, pos_, norm_, texcoord_, radius_);
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            for (int i = 0; i < *nverts; i++) {
                pos[i] = ys__normalize3f(pos[i]);
                norm[i] = pos[i];
            }
        } break;
        case ys_stype_uvspherizedcube: {
            assert(params);
            ys_make_stdshape(ys_stype_uvcube, level, etype, 0, nelems, elem,
                             nverts, pos_, norm_, texcoord_, radius_);
            if (params[0] == 0) return;
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            for (int i = 0; i < *nverts; i++) {
                norm[i] = ys__normalize3f(pos[i]);
                pos[i] = ys__smul3f(pos[i], 1 - params[0]);
                pos[i] = ys__sum3f(pos[i], ys__smul3f(norm[i], params[0]));
            }
            ys_compute_normals(*nelems, *elem, etype, *nverts, *pos_, *norm_,
                               true);
        } break;
        case ys_stype_uvflipcapsphere: {
            assert(params);
            ys_make_stdshape(ys_stype_uvsphere, level, etype, 0, nelems, elem,
                             nverts, pos_, norm_, texcoord_, radius_);
            if (params[0] == 1) return;
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            for (int i = 0; i < *nverts; i++) {
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
            if (radius_) assert(params);
            int npoints = powf(2, level);
            *nverts = npoints;
            *pos_ = (float*)calloc(*nverts * 3, sizeof(float));
            *norm_ = (float*)calloc(*nverts * 3, sizeof(float));
            *texcoord_ = (float*)calloc(*nverts * 2, sizeof(float));
            if (radius_) *radius_ = (float*)calloc(*nverts * 1, sizeof(float));
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            ys__vec2f* texcoord = (ys__vec2f*)*texcoord_;
            for (int i = 0; i < npoints; i++) {
                pos[i] = (ys__vec3f){ 0, 0, 0 };
                norm[i] = (ys__vec3f){ 0, 0, 1 };
                texcoord[i] = (ys__vec2f){ 0.5f, 0.5f };
                if (radius_) (*radius_)[i] = params[0];
            }
            *nelems = npoints;
            *elem = (int*)calloc(npoints, sizeof(int));
            for (int i = 0; i < npoints; i++) {
                (*elem)[i] = i;
            }
        } break;
        case ys_stype_rnpoints: {
            ys_make_stdshape(ys_stype_points, level, etype, 0, nelems, elem,
                             nverts, pos_, norm_, texcoord_, radius_);
            ys__vec3f* pos = (ys__vec3f*)*pos_;
            ys__vec3f* norm = (ys__vec3f*)*norm_;
            ys__vec2f* texcoord = (ys__vec2f*)*texcoord_;
            unsigned int rn = 0;
            for (int i = 0; i < *nverts; i++) {
                pos[i] = (ys__vec3f){ -1 + 2 * ys__rng_nextf(&rn),
                                      -1 + 2 * ys__rng_nextf(&rn),
                                      -1 + 2 * ys__rng_nextf(&rn) };
                norm[i] = ys__normalize3f((ys__vec3f){
                    -1 + 2 * ys__rng_nextf(&rn), -1 + 2 * ys__rng_nextf(&rn),
                    -1 + 2 * ys__rng_nextf(&rn) });
                texcoord[i] =
                    (ys__vec2f){ ys__rng_nextf(&rn), ys__rng_nextf(&rn) };
                if (radius_) {
                    (*radius_)[i] =
                        params[0] +
                        (params[1] - params[0]) * ys__rng_nextf(&rn);
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
YS_API ys_edge_map*
ys_make_edge_map(int nelems, int* elem, int etype) {
    assert(etype >= 2 && etype <= 4);
    ys_edge_map* map = (ys_edge_map*)calloc(1, sizeof(ys_edge_map));

    // find largest vert index
    map->nverts = 0;
    for (int i = 0; i < nelems; i++) {
        int* e = elem + i * etype;
        for (int j = 0; j < etype; j++) {
            if (map->nverts < e[j] + 1) map->nverts = e[j] + 1;
        }
    }

    // count the number of edges per vert
    map->vert_size = (int*)calloc(map->nverts, sizeof(int));
    for (int i = 0; i < nelems; i++) {
        int* e = elem + i * etype;
        for (int j = 0; j < etype; j++) {
            map->vert_size[e[j]]++;
        }
    }

    // prepare half edge data
    map->vert_hedge = (int**)calloc(map->nverts, sizeof(int*));
    for (int i = 0; i < map->nverts; i++) {
        map->vert_hedge[i] = (int*)calloc(map->vert_size[i], 2 * sizeof(int));
    }

    // insert all half edges and set edge counts
    map->nedges = 0;
    memset(map->vert_size, 0, sizeof(int) * map->nverts);
    for (int i = 0; i < nelems; i++) {
        int* e = elem + i * etype;
        for (int j = 0; j < etype; j++) {
            int e1 = e[j], e2 = e[(j + 1) % etype];
            if (j == 2 && e1 == e2) continue;
            int* map_entry = map->vert_hedge[e1] + 2 * map->vert_size[e1];
            map->vert_size[e1] += 1;
            map_entry[0] = e2;
            int pos = -1;
            for (int k = 0; k < map->vert_size[e2] && pos < 0; k++) {
                if (map->vert_hedge[e2][k * 2] == e1) pos = k * 2;
            }
            if (pos < 0) {
                map_entry[1] = map->nedges;
                map->nedges += 1;
            } else {
                map_entry[1] = map->vert_hedge[e2][pos + 1];
            }
        }
    }

    // create edge list
    map->edges = (int*)calloc(map->nedges * 2, sizeof(int));
    for (int i = 0; i < map->nverts; i++) {
        for (int j = 0; j < map->vert_size[i]; j++) {
            int* map_entry = map->vert_hedge[i] + j * 2;
            map->edges[2 * map_entry[1] + 0] = i;
            map->edges[2 * map_entry[1] + 1] = map_entry[0];
        }
    }

    // ensure edge list  is with edges with e1<e2
    for (int i = 0; i < map->nedges; i++) {
        if (map->edges[i * 2] > map->edges[i * 2 + 1]) {
            int aux = map->edges[i * 2 + 1];
            map->edges[i * 2 + 1] = map->edges[i * 2 + 0];
            map->edges[i * 2 + 0] = aux;
        }
    }

    // done
    return map;
}

//
// Get edge index. Public API described above.
//
YS_API int
ys_get_edge_index(const ys_edge_map* em, int v0, int v1) {
    for (int i = 0; i < em->vert_size[v0]; i++) {
        if (em->vert_hedge[v0][i * 2] == v1)
            return em->vert_hedge[v0][i * 2 + 1];
    }
    return -1;
}

//
// Free edge map. Public API described above.
//
YS_API void
ys_free_edge_map(ys_edge_map* em) {
    if (em->edges) free(em->edges);
    if (em->vert_size) free(em->vert_size);
    for (int i = 0; i < em->nverts; i++) {
        if (em->vert_hedge[i]) free(em->vert_hedge[i]);
    }
}

#endif

#endif
