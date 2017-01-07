///
/// YOCTO_SHAPE: A set of utilities to manipulate 3D shapes represented as
/// collection of elements.
///
///
/// INCLUDED UTILITIES:
///
/// 1. include this file (more compilation options below)
/// 2. smoothed normals
///   - compute smoothed vertex normals or tangets with compute_normals()
/// 3. shape tesselation
///   - tesselate the shape by splitting each element along its edges with
///      split_edges()
///   - for a simpler interface, use tesselate_stdshape() that also prepare
///       computes split vertex data for a shape with per-vertex position,
///       normals, texture coordinates, radia and colors
/// 4. create shapes parametrically using callbacks for vertex position, normal
///    and texture coordinates
///    - make a uv surface with make_uvsurface()
///    - make lines with make_lines()
///    - make points with make_points()
/// 5. create shapes a few standard surfaces for testing with make_stdsurface()
/// 6. pick points on a shape
/// 7. sample shapes uniformly
///    - compute shape element distribution with sample_shape_cdf()
///    - generate elements ids and uvs for each point by repeatedly call
///      sample_shape()
/// 8. interpolate vertex data linearly over primitives
///    interpolate_vert()
///
/// The interface for each function is described in details in the interface
/// section of this file.
///
/// Shapes are indexed meshes and are described by array of vertex indices for
/// points, lines and triangles, and arrays of arbitrary vertex data.
/// We differ from other library since vertex data is always arbitrary and
/// we require only vertex position in most functions.
///
///
/// COMPILATION:
///
/// All functions in this library are inlined by default for ease of use in C++.
/// To use the library as a .h/.cpp pair do the following:
/// - to use as a .h, just define YGL_DECLARATION before including this file
/// - to build as a .cpp, just define YGL_IMPLEMENTATION before including this
/// file into only one file that you can either link directly or pack as a lib.
///
/// This file depends on yocto_math.h.
///
///
/// HISTORY:
/// - v 0.7: doxygen comments
/// - v 0.6: opaque API (allows for changing internals without altering API)
/// - v 0.5: internally use pointers for performance transparency
/// - v 0.4: [major API change] move to modern C++ interface
/// - v 0.4: added some cleaner C++ interface
/// - v 0.3: removal of C interface
/// - v 0.2: use of STL containers
/// - v 0.1: C++ implementation
/// - v 0.0: initial release in C99
///
namespace yshape {}

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

#ifndef _YSHAPE_H_
#define _YSHAPE_H_

// compilation options
#ifndef YGL_DECLARATION
#define YGL_API inline
#else
#define YGL_API
#endif

#include <array>
#include <functional>
#include <vector>

// -----------------------------------------------------------------------------
// C++ INTERFACE
// -----------------------------------------------------------------------------

namespace yshape {

//
// Typedefs for vec/mat types
//
using float2 = std::array<float, 2>;
using float3 = std::array<float, 3>;
using float4 = std::array<float, 4>;
using float3x4 = std::array<std::array<float, 3>, 4>;
using int2 = std::array<int, 2>;
using int3 = std::array<int, 3>;
using int4 = std::array<int, 4>;

///
/// Compute smoothed normals or tangents (for lines). Sets points normals to
/// defaults.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - npoints/points: array of point indices
/// - nlines/lines: array of point indices
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - norm: preallocated array of computed normals
///
YGL_API void compute_normals(int npoints, const int* points, int nlines,
                             const int2* lines, int ntriangles,
                             const int3* triangles, int nverts,
                             const float3* pos, float3* norm,
                             bool weighted = true);

///
/// Tesselates a mesh by subdiving along element edges. Will produce a new
/// set of elements referrring to new vertex indices in the range [0,nverts),
/// for the original vertices and [nverts,nverts+nedges) for vertices in the
/// split edges. The edges are returned so new vertices can be created.
/// For a simpler interface, see tesselate_stdshape that handles everything
/// internally.
///
/// Parameters:
/// - nverts: number of vertices
/// - line/triangle: elements to split
///
/// Out Parameters:
/// - tess_line/tess_triangle: split elements
/// - tess_edges: edges created
///
YGL_API void split_edges(int nverts, int nlines, const int2* lines,
                         int ntriangles, const int3* triangles,
                         std::vector<int2>& tess_lines,
                         std::vector<int3>& tess_triangles,
                         std::vector<int2>& tess_edges);

///
/// Tesselate a shape inplace.
///
/// In/Out Parameters:
/// - points, lines, triangles: elems to split
/// - pos, norm, texcoord, color, radius: vertices to split
///
YGL_API void tesselate_stdshape(std::vector<int2>& lines,
                                std::vector<int3>& triangles,
                                std::vector<float3>& pos,
                                std::vector<float3>& norm,
                                std::vector<float2>& texcoord,
                                std::vector<float3>& color,
                                std::vector<float>& radius);

///
/// Generate a parametric surface with callbacks.
///
/// Parameters:
/// - usteps: subdivisions in u
/// - vsteps: subdivisions in v
/// - pos_fn/norm_fn/texcoord_fn: pos/norm/texcoord callbacks
///
/// Out Parameters:
/// - triangles: element array
/// - pos/norm/texcoord: vertex position/normal/texcoords
///
YGL_API void make_uvsurface(int usteps, int vsteps,
                            std::vector<int3>& triangles,
                            std::vector<float3>& pos, std::vector<float3>& norm,
                            std::vector<float2>& texcoord,
                            std::function<float3(const float2&)> pos_fn,
                            std::function<float3(const float2&)> norm_fn,
                            std::function<float2(const float2&)> texcoord_fn);

///
/// Generate parametric lines with callbacks.
///
/// Parameters:
/// - usteps: subdivisions in u
/// - num: number of lines
/// - pos_fn/tang_fn/texcoord_fn/radius_fn: pos/norm/texcoord/radius callbacks
///
/// Out Parameters:
/// - lines: element array
/// - pos/tang/texcoord/radius: vertex position/tangent/texcoords/radius
///
YGL_API void make_lines(int usteps, int num, std::vector<int2>& lines,
                        std::vector<float3>& tang, std::vector<float3>& norm,
                        std::vector<float2>& texcoord,
                        std::vector<float>& radius,
                        std::function<float3(const float2&)> pos_fn,
                        std::function<float3(const float2&)> tang_fn,
                        std::function<float3(const float2&)> texcoord_fn,
                        std::function<float(const float2&)> radius_fn);

///
/// Generate a parametric surface with callbacks.
///
/// Parameters:
/// - num: number of points
/// - pos_fn/norm_fn/texcoord_fn/radius_fn: pos/norm/texcoord/radius callbacks
///
/// Out Parameters:
/// - points: element array
/// - pos/norm/texcoord/radius: vertex position/normal/texcoords/radius
///
YGL_API void make_points(int num, std::vector<int>& points,
                         std::vector<float3>& pos, std::vector<float3>& norm,
                         std::vector<float2>& texcoord,
                         std::vector<float>& radius,
                         std::function<float3(float)> pos_fn,
                         std::function<float3(float)> norm_fn,
                         std::function<float3(float)> texcoord_fn,
                         std::function<float(float)> radius_fn);

///
/// Test shapes (this is mostly used to create tests)
///
enum struct stdsurface_type {
    uvsphere,         ///< uv sphere
    uvquad,           ///< quad
    uvcube,           ///< cube
    uvflippedsphere,  ///< uv sphere flipped inside/out
    uvspherecube,     ///< sphere obtained by a cube tesselation
    uvspherizedcube,  ///< cube tesselation spherized by a radius
    uvflipcapsphere   ///< us sphere with flipped poles
};

///
/// Create standard shapes for testing purposes.
///
/// Parameters:
/// - stype: shape type
/// - level: tesselation level (roughly euivalent to creating 2^level splits)
/// - params: shape params for sample shapes (see code for documentation)
/// - frame: frame
/// - scale: scale
///
/// Out Parameters:
/// - triangles: element array
/// - pos: vertex positions
/// - norm: vertex normals
/// - texcoord: vertex texture coordinates
///
YGL_API void make_stdsurface(
    stdsurface_type stype, int level, const float4& params,
    std::vector<int3>& triangles, std::vector<float3>& pos,
    std::vector<float3>& norm, std::vector<float2>& texcoord,
    const float3x4& frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}},
    float scale = 1);

///
/// Computes the distribution of area of a shape element for sampling. This is
/// needed to sample a shape. In the case of the sample_shape version, only one
/// element array can be non-empty at any given call.
///
/// Paramaters:
/// - npoints/points: number of elements
/// - elem: element array
/// - etype: element type
/// - pos: vertex positions
///
/// Out Parameters:
/// - ecdf: array of element cdfs (or NULL if allocated internally)
/// - area: total area of the shape (or length for lines)
///
YGL_API void sample_shape_cdf(int npoints, const int* points, int nverts,
                              const float3* pos, float* ecdf, float* num);

///
/// Computes the distribution of area of a shape element for sampling. This is
/// needed to sample a shape. In the case of the sample_shape version, only one
/// element array can be non-empty at any given call.
///
/// Paramaters:
/// - npoints/points: number of elements
/// - elem: element array
/// - etype: element type
/// - pos: vertex positions
///
/// Out Parameters:
/// - ecdf: array of element cdfs (or NULL if allocated internally)
/// - area: total area of the shape (or length for lines)
///
YGL_API void sample_shape_cdf(int nlines, const int2* lines, int nverts,
                              const float3* pos, float* ecdf, float* length);

///
/// Computes the distribution of area of a shape element for sampling. This is
/// needed to sample a shape. In the case of the sample_shape version, only one
/// element array can be non-empty at any given call.
///
/// Paramaters:
/// - npoints/points: number of elements
/// - elem: element array
/// - etype: element type
/// - pos: vertex positions
///
/// Out Parameters:
/// - ecdf: array of element cdfs (or NULL if allocated internally)
/// - area: total area of the shape (or length for lines)
///
YGL_API void sample_shape_cdf(int ntriangles, const int2* triangles, int nverts,
                              const float3* pos, float* ecdf, float* length);

///
/// Sampels a shape elements. In the case of the sample_shape version, only one
/// cdf can be non-empty at any given call.
///
/// Paramaters:
/// - nelems: number of elements
/// - point_cdf/line_cdf/triangle_cdf: element cdfs from the above std::function
/// - ern, uvrn: random numbers, int [0,1) range, for element and uv choices
///
/// Out Parameters:
/// - eid: element index
/// - euv: element baricentric coordinates
///
YGL_API void sample_shape(int npoints, const float* point_cdf, float ern,
                          int* eid);

///
/// Sampels a shape elements. In the case of the sample_shape version, only one
/// cdf can be non-empty at any given call.
///
/// Paramaters:
/// - nelems: number of elements
/// - point_cdf/line_cdf/triangle_cdf: element cdfs from the above std::function
/// - ern, uvrn: random numbers, int [0,1) range, for element and uv choices
///
/// Out Parameters:
/// - eid: element index
/// - euv: element baricentric coordinates
///
YGL_API void sample_shape(int nlines, const float* line_cdf, float ern,
                          float uvrn, int* eid, float* euv);

///
/// Sampels a shape elements. In the case of the sample_shape version, only one
/// cdf can be non-empty at any given call.
///
/// Paramaters:
/// - nelems: number of elements
/// - point_cdf/line_cdf/triangle_cdf: element cdfs from the above std::function
/// - ern, uvrn: random numbers, int [0,1) range, for element and uv choices
///
/// Out Parameters:
/// - eid: element index
/// - euv: element baricentric coordinates
///
YGL_API void sample_shape(int ntriangles, const float* triangle_cdf, float ern,
                          const float2& uvrn, int* eid, float2* euv);

///
/// Interpolates a vertex property using baricentric interpolation. Uses
/// linear interpolation for lines, baricentric for triangles and copies values
/// for points.
///
/// Parameters:
/// - points/lines/triangles: array of vertex indices (only one filled)
/// - vert: vertex property array
/// - eid: element index
/// - euv: element parameters
///
/// Returns:
/// - interpolated vertex data
///
template <typename T>
YGL_API T interpolate_vert(int nlines, const int2* lines, int nverts,
                           const T* vert, int eid, const float2& euv);

///
/// Interpolates a vertex property using baricentric interpolation. Uses
/// linear interpolation for lines, baricentric for triangles and copies values
/// for points.
///
/// Parameters:
/// - points/lines/triangles: array of vertex indices (only one filled)
/// - vert: vertex property array
/// - eid: element index
/// - euv: element parameters
///
/// Returns:
/// - interpolated vertex data
///
template <typename T>
YGL_API T interpolate_vert(int ntriangles, const int3* triangles, int nverts,
                           const T* vert, int eid, const float2& euv);

}  // namespace

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------

#include <unordered_map>

#include "yocto_math.h"

#if !defined(YGL_DECLARATION) || defined(YGL_IMPLEMENTATION)

#include <array>
#include <cassert>
#include <cstdlib>
#include <cstring>

namespace yshape {

//
// Dictionary from directed edges to undirected edges implemented as a hashmap.
// To insert values into the data structure use insert. To access elements use
// the overloaded operators [] and at. For convenience, use make_edge_map.
//
struct edge_map {
    // hash type
    struct _hash {
        size_t operator()(const ym::vec2i& v) const { return ym::hash_vec(v); }
    };

    // underlying map type
    using map_t = std::unordered_map<ym::vec2i, int, _hash>;

    // size and iteration
    size_t size() const { return _map.size(); }

    // edge insertion
    void insert(const ym::vec2i& e) {
        if (!has_edge(e)) _map[_edge(e)] = (int)_map.size();
    }

    // edge check
    bool has_edge(const ym::vec2i& e) const {
        return _map.find(_edge(e)) != _map.end();
    }

    // iteration
    map_t::const_iterator begin() const { return _map.begin(); }
    map_t::const_iterator end() const { return _map.end(); }

    // element access (returns -1 if not present)
    int operator[](const ym::vec2i& e) const {
        assert(has_edge(e));
        if (!has_edge(e))
            return -1;
        else
            return _map.at(_edge(e));
    }

   private:
    // edge id [private]
    ym::vec2i _edge(const ym::vec2i& e) const {
        return {ym::min(e[0], e[1]), ym::max(e[0], e[1])};
    }

    // data [private]
    map_t _map;
};

//
// Build an edge map
//
static inline edge_map make_edge_map(
    ym::array_view<const ym::vec2i> lines,
    ym::array_view<const ym::vec3i> triangles) {
    auto map = edge_map();

    for (auto l : lines) {
        map.insert(l);
    }
    for (auto t : triangles) {
        for (auto i = 0; i < 3; i++) {
            map.insert({t[i], t[(i + 1) % 3]});
        }
    }

    // done
    return map;
}

//
// Normal computation. Public API described above.
//
static inline void _compute_normals(ym::array_view<const int> points,
                                    ym::array_view<const ym::vec2i> lines,
                                    ym::array_view<const ym::vec3i> triangles,
                                    ym::array_view<const ym::vec3f> pos,
                                    ym::array_view<ym::vec3f> norm,
                                    bool weighted) {
    // clear normals
    for (auto& n : norm) n = ym::zero3f;

    // handle various primitives
    for (auto p : points) norm[p] += ym::vec3f{0, 0, 1};
    for (auto l : lines) {
        auto n = pos[l[1]] - pos[l[0]];
        if (!weighted) n = ym::normalize(n);
        norm[l[0]] += n;
        norm[l[1]] += n;
    }
    for (auto t : triangles) {
        auto n = ym::cross(pos[t[1]] - pos[t[0]], pos[t[2]] - pos[t[0]]);
        if (!weighted) n = ym::normalize(n);
        norm[t[0]] += n;
        norm[t[1]] += n;
        norm[t[2]] += n;
    }

    // normalize result
    for (auto& n : norm) n = ym::normalize(n);
}

//
// Compute smoothed normals or tangents (for lines). Sets points normals to
// defaults. Public API.
//
YGL_API void compute_normals(int npoints, const int* points, int nlines,
                             const int2* lines, int ntriangles,
                             const int3* triangles, int nverts,
                             const float3* pos, float3* norm, bool weighted) {
    return _compute_normals({(size_t)npoints, points},
                            {(size_t)nlines, (const ym::vec2i*)lines},
                            {(size_t)ntriangles, (const ym::vec3i*)triangles},
                            {(size_t)nverts, (const ym::vec3f*)pos},
                            {(size_t)nverts, (ym::vec3f*)norm}, weighted);
}

//
// Tesselate a shape by splitting its edges
//
YGL_API void _split_edges(int nverts, ym::array_view<const ym::vec2i> lines,
                          ym::array_view<const ym::vec3i> triangles,
                          std::vector<ym::vec2i>& tess_lines,
                          std::vector<ym::vec3i>& tess_triangles,
                          std::vector<ym::vec2i>& edges) {
    // grab edges
    auto em = make_edge_map(lines, triangles);

    // make new elements
    tess_lines.clear();
    tess_lines.reserve(lines.size() * 2);
    for (auto l : lines) {
        tess_lines.push_back({l[0], nverts + em[l]});
        tess_lines.push_back({nverts + em[l], l[1]});
    }
    tess_triangles.clear();
    tess_triangles.reserve(triangles.size() * 4);
    for (auto t : triangles) {
        for (auto i = 0; i < 3; i++) {
            tess_triangles.push_back({t[i], nverts + em[{t[i], t[(i + 1) % 3]}],
                                      nverts + em[{t[i], t[(i + 2) % 3]}]});
        }
        tess_triangles.push_back({nverts + em[{t[0], t[1]}],
                                  nverts + em[{t[1], t[2]}],
                                  nverts + em[{t[2], t[0]}]});
    }

    // returned edges
    edges.resize(em.size());
    for (auto e : em) {
        edges[e.second] = e.first;
    }
}

//
// Public API.
//
YGL_API void split_edges(int nverts, int nlines, const int2* lines,
                         int ntriangles, const int3* triangles,
                         std::vector<int2>& tess_lines,
                         std::vector<int3>& tess_triangles,
                         std::vector<int2>& tess_edges) {
    return _split_edges(nverts, {(size_t)nlines, (const ym::vec2i*)lines},
                        {(size_t)ntriangles, (const ym::vec3i*)triangles},
                        (std::vector<ym::vec2i>&)tess_lines,
                        (std::vector<ym::vec3i>&)tess_triangles,
                        (std::vector<ym::vec2i>&)tess_edges);
}

//
// Tesselate a shape inplace.
//
YGL_API void _tesselate_stdshape(std::vector<ym::vec2i>& lines,
                                 std::vector<ym::vec3i>& triangles,
                                 std::vector<ym::vec3f>& pos,
                                 std::vector<ym::vec3f>& norm,
                                 std::vector<ym::vec2f>& texcoord,
                                 std::vector<ym::vec3f>& color,
                                 std::vector<float>& radius) {
    // get the number of vertices
    auto nverts = (int)pos.size();

    // prepare edges and elements
    std::vector<ym::vec2i> tess_lines;
    std::vector<ym::vec3i> tess_triangles;
    std::vector<ym::vec2i> tess_edges;
    _split_edges(nverts, lines, triangles, tess_lines, tess_triangles,
                 tess_edges);
    lines = tess_lines;
    triangles = tess_triangles;

    // allocate vertex data
    if (!pos.empty()) pos.resize(nverts + tess_edges.size());
    if (!norm.empty()) norm.resize(nverts + tess_edges.size());
    if (!texcoord.empty()) texcoord.resize(nverts + tess_edges.size());
    if (!color.empty()) color.resize(nverts + tess_edges.size());
    if (!radius.empty()) radius.resize(nverts + tess_edges.size());

    // interpolate vertex data
    for (auto j = 0; j < tess_edges.size(); j++) {
        auto e = tess_edges[j];
        if (!pos.empty()) pos[nverts + j] = (pos[e[0]] + pos[e[1]]) / 2.0f;
        if (!norm.empty()) norm[nverts + j] = (norm[e[0]] + norm[e[1]]) / 2.0f;
        if (!texcoord.empty())
            texcoord[nverts + j] = (texcoord[e[0]] + texcoord[e[1]]) / 2.0f;
        if (!color.empty())
            color[nverts + j] = (color[e[0]] + color[e[1]]) / 2.0f;
        if (!radius.empty())
            radius[nverts + j] = (radius[e[0]] + radius[e[1]]) / 2.0f;
    }

    // fix normals
    for (auto& n : norm) n = ym::normalize(n);
}

//
// Tesselate a shape inplace. Public API.
//
YGL_API void tesselate_stdshape(std::vector<int2>& lines,
                                std::vector<int3>& triangles,
                                std::vector<float3>& pos,
                                std::vector<float3>& norm,
                                std::vector<float2>& texcoord,
                                std::vector<float3>& color,
                                std::vector<float>& radius) {
    return _tesselate_stdshape(
        (std::vector<ym::vec2i>&)lines, (std::vector<ym::vec3i>&)triangles,
        (std::vector<ym::vec3f>&)pos, (std::vector<ym::vec3f>&)norm,
        (std::vector<ym::vec2f>&)texcoord, (std::vector<ym::vec3f>&)color,
        radius);
}

//
// Creates a standard surface. Public interface.
//
static inline void _make_uvsurface(
    int usteps, int vsteps, std::vector<ym::vec3i>& triangles,
    std::vector<ym::vec3f>& pos, std::vector<ym::vec3f>& norm,
    std::vector<ym::vec2f>& texcoord,
    std::function<ym::vec3f(const ym::vec2f&)> pos_fn,
    std::function<ym::vec3f(const ym::vec2f&)> norm_fn,
    std::function<ym::vec2f(const ym::vec2f&)> texcoord_fn) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    pos.resize((usteps + 1) * (vsteps + 1));
    norm.resize((usteps + 1) * (vsteps + 1));
    texcoord.resize((usteps + 1) * (vsteps + 1));
    for (auto j = 0; j <= vsteps; j++) {
        for (auto i = 0; i <= usteps; i++) {
            auto uv = ym::vec2f{i / (float)usteps, j / (float)vsteps};
            pos[vid(i, j)] = pos_fn(uv);
            norm[vid(i, j)] = norm_fn(uv);
            texcoord[vid(i, j)] = texcoord_fn(uv);
        }
    }

    triangles.resize(usteps * vsteps * 2);
    for (auto j = 0; j < vsteps; j++) {
        for (auto i = 0; i < usteps; i++) {
            auto& f1 = triangles[(j * usteps + i) * 2 + 0];
            auto& f2 = triangles[(j * usteps + i) * 2 + 1];
            if ((i + j) % 2) {
                f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
                f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
            } else {
                f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
                f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
            }
        }
    }
}

//
// Creates a standard surface. Public interface.
//
YGL_API void make_uvsurface(int usteps, int vsteps,
                            std::vector<int3>& triangles,
                            std::vector<float3>& pos, std::vector<float3>& norm,
                            std::vector<float2>& texcoord,
                            std::function<float3(const float2&)> pos_fn,
                            std::function<float3(const float2&)> norm_fn,
                            std::function<float2(const float2&)> texcoord_fn) {
    _make_uvsurface(usteps, vsteps, (std::vector<ym::vec3i>&)triangles,
                    (std::vector<ym::vec3f>&)pos, (std::vector<ym::vec3f>&)norm,
                    (std::vector<ym::vec2f>&)texcoord,
                    (std::function<ym::vec3f(const ym::vec2f&)>)pos_fn,
                    (std::function<ym::vec3f(const ym::vec2f&)>)norm_fn,
                    (std::function<ym::vec2f(const ym::vec2f&)>)texcoord_fn);
}

//
// Creates lines. Public interface.
//
YGL_API void make_lines(int usteps, int num, std::vector<int2>& lines,
                        std::vector<float3>& pos, std::vector<float3>& norm,
                        std::vector<float2>& texcoord,
                        std::vector<float>& radius,
                        std::function<float3(const float2&)> pos_fn,
                        std::function<float3(const float2&)> norm_fn,
                        std::function<float2(const float2&)> texcoord_fn,
                        std::function<float(const float2&)> radius_fn) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    pos.resize((usteps + 1) * num);
    norm.resize((usteps + 1) * num);
    texcoord.resize((usteps + 1) * num);
    radius.resize((usteps + 1) * num);
    for (auto j = 0; j < num; j++) {
        for (auto i = 0; i <= usteps; i++) {
            auto uv = ym::vec2f{i / (float)usteps, j / (float)(num - 1)};
            pos[vid(i, j)] = pos_fn(uv);
            norm[vid(i, j)] = norm_fn(uv);
            texcoord[vid(i, j)] = texcoord_fn(uv);
            radius[vid(i, j)] = radius_fn(uv);
        }
    }

    lines.resize(usteps * num);
    for (int j = 0; j < num; j++) {
        for (int i = 0; i < usteps; i++) {
            lines[j * usteps + i] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

//
// Tesselates a surface. Public interface.
//
YGL_API void make_points(int num, std::vector<int>& points,
                         std::vector<float3>& pos, std::vector<float3>& norm,
                         std::vector<float2>& texcoord,
                         std::vector<float>& radius,
                         std::function<float3(float)> pos_fn,
                         std::function<float3(float)> norm_fn,
                         std::function<float2(float)> texcoord_fn,
                         std::function<float(float)> radius_fn) {
    pos.resize(num);
    norm.resize(num);
    texcoord.resize(num);
    radius.resize(num);
    for (auto i = 0; i < num; i++) {
        auto u = i / (float)i;
        pos[i] = pos_fn(u);
        norm[i] = norm_fn(u);
        texcoord[i] = texcoord_fn(u);
        radius[i] = radius_fn(u);
    }

    points.resize(num);
    for (auto i = 0; i < num; i++) points[i] = i;
}

//
// Sample cdf. Public API described above.
//
static inline void _sample_shape_cdf(ym::array_view<const int> elems,
                                     ym::array_view<const ym::vec3f> pos,
                                     ym::array_view<float> cdf, float& weight) {
    for (auto i = 0; i < elems.size(); i++) cdf[i] = 1;
    for (auto i = 1; i < elems.size(); i++) cdf[i] += cdf[i - 1];
    weight = cdf.back();
    for (auto i = 0; i < elems.size(); i++) cdf[i] /= cdf.back();
}

//
// Sample cdf. Public API described above.
//
static inline void _sample_shape_cdf(ym::array_view<const ym::vec2i> elems,
                                     ym::array_view<const ym::vec3f> pos,
                                     ym::array_view<float> cdf, float& weight) {
    for (auto i = 0; i < elems.size(); i++) {
        auto& f = elems[i];
        cdf[i] = ym::length(pos[f[0]] - pos[f[1]]);
    }
    for (auto i = 1; i < elems.size(); i++) cdf[i] += cdf[i - 1];
    weight = cdf.back();
    for (auto i = 0; i < elems.size(); i++) cdf[i] /= cdf.back();
}

//
// Sample cdf. Public API described above.
//
static inline void _sample_shape_cdf(ym::array_view<const ym::vec3i> elems,
                                     ym::array_view<const ym::vec3f> pos,
                                     ym::array_view<float> cdf, float& weight) {
    for (auto i = 0; i < elems.size(); i++) {
        auto& f = elems[i];
        cdf[i] = ym::length(
                     ym::cross(pos[f[0]] - pos[f[1]], pos[f[0]] - pos[f[2]])) /
                 2;
    }
    for (auto i = 1; i < elems.size(); i++) cdf[i] += cdf[i - 1];
    weight = cdf.back();
    for (auto i = 0; i < elems.size(); i++) cdf[i] /= cdf.back();
}

// finds the first array element smaller than the given one
// http://stackoverflow.com/questions/6553970/
// find-the-first-element-in-an-array-that-is-greater-than-the-target
static inline size_t _bsearch_smaller(float x, const float a[], size_t n) {
    size_t low = 0, high = n;
    while (low != high) {
        size_t mid = (low + high) / 2;
        if (a[mid] < x)
            low = mid + 1;
        else
            high = mid;
    }
    return low;
}

//
// Sample shape. Public API described above.
//
static inline void _sample_shape(ym::array_view<const float> cdf, float ern,
                                 int& eid) {
    eid = (int)_bsearch_smaller(ern, cdf.data(), cdf.size());
}

//
// Sample shape. Public API described above.
//
static inline void _sample_shape(ym::array_view<const float> cdf, float ern,
                                 float uvrn, int& eid, float& euv) {
    eid = (int)_bsearch_smaller(ern, cdf.data(), cdf.size());
    euv = uvrn;
}

//
// Sample shape. Public API described above.
//
static inline void _sample_shape(ym::array_view<const float> cdf, float ern,
                                 const ym::vec2f& uvrn, int& eid,
                                 ym::vec2f& euv) {
    eid = (int)_bsearch_smaller(ern, cdf.data(), cdf.size());
    euv = {1 - sqrtf(uvrn[0]), uvrn[1] * sqrtf(uvrn[0])};
}

//
// Interpolate vertex properties. Public API.
//
template <typename T>
static inline T _interpolate_vert(ym::array_view<const ym::vec2i> lines,
                                  ym::array_view<const T> vert, int eid,
                                  const ym::vec2f& euv) {
    return vert[lines[eid][0]] * (1 - euv[0]) + vert[lines[eid][1]] * euv[1];
}

//
// Interpolate vertex properties. Public API.
//
template <typename T>
static inline T _interpolate_vert(ym::array_view<const ym::vec3i> triangles,
                                  ym::array_view<const T> vert, int eid,
                                  const ym::vec2f& euv) {
    return vert[triangles[eid][0]] * (1 - euv[0] - euv[1]) +
           vert[triangles[eid][1]] * euv[0] + vert[triangles[eid][2]] * euv[1];
}

//
// Make standard shape. Public API described above.
//
static inline void _make_stdsurface(stdsurface_type stype, int level,
                                    const ym::vec4f& params,
                                    std::vector<ym::vec3i>& triangles,
                                    std::vector<ym::vec3f>& pos,
                                    std::vector<ym::vec3f>& norm,
                                    std::vector<ym::vec2f>& texcoord,
                                    const ym::frame3f& frame, float scale) {
    switch (stype) {
        case stdsurface_type::uvsphere: {
            auto usteps = ym::pow2(level + 2), vsteps = ym::pow2(level + 1);
            return _make_uvsurface(
                usteps, vsteps, triangles, pos, norm, texcoord,
                [frame, scale](const auto uv) {
                    auto a =
                        ym::vec2f{2 * ym::pif * uv[0], ym::pif * (1 - uv[1])};
                    return transform_point(
                        frame, {scale * std::cos(a[0]) * std::sin(a[1]),
                                scale * std::sin(a[0]) * std::sin(a[1]),
                                scale * std::cos(a[1])});
                },
                [frame](const auto uv) {
                    auto a =
                        ym::vec2f{2 * ym::pif * uv[0], ym::pif * (1 - uv[1])};
                    return transform_direction(frame,
                                               {std::cos(a[0]) * std::sin(a[1]),
                                                std::sin(a[0]) * std::sin(a[1]),
                                                std::cos(a[1])});
                },
                [](const auto uv) { return uv; });
        } break;
        case stdsurface_type::uvflippedsphere: {
            auto usteps = ym::pow2(level + 2), vsteps = ym::pow2(level + 1);
            return _make_uvsurface(
                usteps, vsteps, triangles, pos, norm, texcoord,
                [frame, scale](const auto uv) {
                    auto a = ym::vec2f{2 * ym::pif * uv[0], ym::pif * uv[1]};
                    return transform_point(
                        frame, {scale * std::cos(a[0]) * std::sin(a[1]),
                                scale * std::sin(a[0]) * std::sin(a[1]),
                                scale * std::cos(a[1])});
                },
                [frame](const auto uv) {
                    auto a = ym::vec2f{2 * ym::pif * uv[0], ym::pif * uv[1]};
                    return transform_direction(
                        frame,
                        {-std::cos(a[0]) * std::sin(a[1]),
                         -std::sin(a[0]) * std::sin(a[1]), -std::cos(a[1])});
                },
                [](const auto uv) {
                    return ym::vec2f{uv[0], 1 - uv[1]};
                });
        } break;
        case stdsurface_type::uvquad: {
            auto usteps = ym::pow2(level), vsteps = ym::pow2(level);
            return _make_uvsurface(
                usteps, vsteps, triangles, pos, norm, texcoord,
                [frame, scale](const auto uv) {
                    return transform_point(frame, {-1 + uv[0] * 2 * scale,
                                                   -1 + uv[1] * 2 * scale, 0});
                },
                [frame](const auto uv) {
                    return transform_direction(frame, {0, 0, 1});
                },
                [](const auto uv) {
                    return ym::vec2f{uv[0], uv[1]};
                });
        } break;
        case stdsurface_type::uvcube: {
            auto frames = std::array<ym::frame3f, 6>{
                ym::frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}},
                ym::frame3f{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, -1}},
                ym::frame3f{{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0}},
                ym::frame3f{{1, 0, 0}, {0, 0, 1}, {0, -1, 0}, {0, -1, 0}},
                ym::frame3f{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}},
                ym::frame3f{{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}, {-1, 0, 0}}};
            std::vector<ym::vec3f> quad_pos, quad_norm;
            std::vector<ym::vec2f> quad_texcoord;
            std::vector<ym::vec3i> quad_triangles;
            for (auto frame : frames) {
                auto offset = ym::vec3i{(int)pos.size(), (int)pos.size(),
                                        (int)pos.size()};
                _make_stdsurface(stdsurface_type::uvquad, level, params,
                                 quad_triangles, quad_pos, quad_norm,
                                 quad_texcoord, frame, scale);
                for (auto p : quad_pos) pos.push_back(p);
                for (auto n : quad_norm) norm.push_back(n);
                for (auto t : quad_texcoord) texcoord.push_back(t);
                for (auto t : quad_triangles) triangles.push_back(t + offset);
            }
        } break;
        case stdsurface_type::uvspherecube: {
            _make_stdsurface(stdsurface_type::uvcube, level, ym::zero4f,
                             triangles, pos, norm, texcoord,
                             ym::identity_frame3f, 1);
            for (auto i = 0; i < pos.size(); i++) {
                pos[i] = transform_point(frame, scale * ym::normalize(pos[i]));
                norm[i] = ym::normalize(pos[i]);
            }
        } break;
        case stdsurface_type::uvspherizedcube: {
            _make_stdsurface(stdsurface_type::uvcube, level, ym::zero4f,
                             triangles, pos, norm, texcoord,
                             ym::identity_frame3f, 1);
            if (params[0] != 0) {
                for (auto i = 0; i < pos.size(); i++) {
                    norm[i] = ym::normalize(pos[i]);
                    pos[i] *= 1 - params[0];
                    pos[i] += norm[i] * params[0];
                }
                _compute_normals({}, {}, triangles, pos, norm, true);
            }
        } break;
        case stdsurface_type::uvflipcapsphere: {
            _make_stdsurface(stdsurface_type::uvsphere, level, ym::zero4f,
                             triangles, pos, norm, texcoord,
                             ym::identity_frame3f, 1);
            if (params[0] != 1) {
                for (auto i = 0; i < pos.size(); i++) {
                    if (pos[i][2] > params[0]) {
                        pos[i][2] = 2 * params[0] - pos[i][2];
                        norm[i][0] = -norm[i][0];
                        norm[i][1] = -norm[i][1];
                    } else if (pos[i][2] < -params[0]) {
                        pos[i][2] = -2 * params[0] - pos[i][2];
                        norm[i][0] = -norm[i][0];
                        norm[i][1] = -norm[i][1];
                    }
                }
            }
        } break;
        default: { assert(false); } break;
    }
}

//
// Make standard shape. Public API described above.
//
YGL_API void make_stdsurface(stdsurface_type stype, int level,
                             const float4& params, std::vector<int3>& triangles,
                             std::vector<float3>& pos,
                             std::vector<float3>& norm,
                             std::vector<float2>& texcoord,
                             const float3x4& frame, float scale) {
    return _make_stdsurface(
        stype, level, (ym::vec4f)params, (std::vector<ym::vec3i>&)triangles,
        (std::vector<ym::vec3f>&)pos, (std::vector<ym::vec3f>&)norm,
        (std::vector<ym::vec2f>&)texcoord, (ym::frame3f)frame, scale);
}

}  // namespace

#endif

#endif
