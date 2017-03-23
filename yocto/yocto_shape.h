///
/// YOCTO_SHAPE: A set of utilities to manipulate 3D shapes represented as
/// collection of elements.
///
///
/// INCLUDED UTILITIES:
///
/// 1. smoothed normals
///   - compute smoothed vertex normals or tangets with compute_normals()
/// 2. shape tesselation
///   - tesselate the shape by splitting each element along its edges with
///      split_edges()
///   - for a simpler interface, use tesselate_stdshape() that also prepare
///       computes split vertex data for a shape with per-vertex position,
///       normals, texture coordinates, radia and colors
/// 3. create shapes parametrically using callbacks for vertex position, normal
///    and texture coordinates
///    - make a uv surface with make_uvsurface()
///    - make lines with make_lines()
///    - make points with make_points()
/// 4. create shapes a few standard surfaces for testing with make_stdsurface()
/// 5. pick points on a shape
/// 6. sample shapes uniformly
///    - compute shape element distribution with sample_shape_cdf()
///    - generate elements ids and uvs for each point by repeatedly call
///      sample_shape()
/// 7. interpolate vertex data linearly over primitives
///    interpolate_vert()
/// 8. compute tangent frame with compute_tangent_frame()
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
/// To use the library include the .h and compile the .cpp. To use this library
/// as a header-only library, define YBVH_INLINE before including this file.
///
/// The .cpp file depends on yocto_math.h.
///
///
/// HISTORY:
/// - v 0.9: tangent frame
/// - v 0.8: switch to .h/.cpp pair
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
#ifdef YSHAPE_INLINE
#define YSHAPE_API inline
#else
#define YSHAPE_API
#endif

#include <array>
#include <functional>
#include <vector>

// -----------------------------------------------------------------------------
// INTERFACE
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
YSHAPE_API void compute_normals(int npoints, const int* points, int nlines,
    const int2* lines, int ntriangles, const int3* triangles, int nverts,
    const float3* pos, float3* norm, bool weighted = true);

///
/// Compute tangent frame for triangle mesh. Tangent space is defined by
/// a four component vector. The first three components are the tangent
/// with respect to the U texcoord. The fourth component is the sign of the
/// tangent wrt the V texcoord. Tangent frame is useful in normal mapping.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - tangsp: preallocated array of computed tangent space
///
YSHAPE_API void compute_tangent_frame(int ntriangles, const int3* triangles,
    int nverts, const float3* pos, const float3* norm, const float2* texcoord,
    float4* tangsp, bool weighted = true);

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
YSHAPE_API void split_edges(int nverts, int nlines, const int2* lines,
    int ntriangles, const int3* triangles, std::vector<int2>& tess_lines,
    std::vector<int3>& tess_triangles, std::vector<int2>& tess_edges);

///
/// Tesselate a shape inplace.
///
/// In/Out Parameters:
/// - points, lines, triangles: elems to split
/// - pos, norm, texcoord, color, radius: vertices to split
///
YSHAPE_API void tesselate_stdshape(std::vector<int2>& lines,
    std::vector<int3>& triangles, std::vector<float3>& pos,
    std::vector<float3>& norm, std::vector<float2>& texcoord,
    std::vector<float3>& color, std::vector<float>& radius);

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
YSHAPE_API void make_uvsurface(int usteps, int vsteps,
    std::vector<int3>& triangles, std::vector<float3>& pos,
    std::vector<float3>& norm, std::vector<float2>& texcoord,
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
YSHAPE_API void make_lines(int usteps, int num, std::vector<int2>& lines,
    std::vector<float3>& pos, std::vector<float3>& tang,
    std::vector<float2>& texcoord, std::vector<float>& radius,
    std::function<float3(const float2&)> pos_fn,
    std::function<float3(const float2&)> tang_fn,
    std::function<float2(const float2&)> texcoord_fn,
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
YSHAPE_API void make_points(int num, std::vector<int>& points,
    std::vector<float3>& pos, std::vector<float3>& norm,
    std::vector<float2>& texcoord, std::vector<float>& radius,
    std::function<float3(float)> pos_fn, std::function<float3(float)> norm_fn,
    std::function<float2(float)> texcoord_fn,
    std::function<float(float)> radius_fn);

///
/// Test shapes (this is mostly used to create tests)
///
enum struct stdsurface_type {
    uvsphere,             ///< uv sphere
    uvhemisphere,         ///< uv hemisphere
    uvquad,               ///< quad
    uvcube,               ///< cube
    uvflippedsphere,      ///< uv sphere flipped inside/out
    uvflippedhemisphere,  ///< uv hemisphere flipped inside/out
    uvspherecube,         ///< sphere obtained by a cube tesselation
    uvspherizedcube,      ///< cube tesselation spherized by a radius
    uvflipcapsphere,      ///< uv sphere with flipped poles
    uvhollowcutsphere,    ///< uv hollow sphere and cut at the top
    uvhollowcutsphere1,   ///< uv hollow sphere and cut at the top
    uvcutsphere,          ///< uv sphere cut at the top
    uvflippedcutsphere,   ///< uv sphere flipped and cut at the top
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
/// - urange, vrange: texture coordinates range
///
/// Out Parameters:
/// - triangles: element array
/// - pos: vertex positions
/// - norm: vertex normals
/// - texcoord: vertex texture coordinates
///
YSHAPE_API void make_stdsurface(stdsurface_type stype, int level,
    const float4& params, std::vector<int3>& triangles,
    std::vector<float3>& pos, std::vector<float3>& norm,
    std::vector<float2>& texcoord,
    const float3x4& frame = {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}}},
    float scale = 1, const float2& urange = {0, 1},
    const float2& vrange = {0, 1});

///
/// Merge a standard shapes elements into an existing shape.
///
/// In/Out Parameters:
/// - triangles: element array
/// - pos: vertex positions
/// - norm: vertex normals
/// - texcoord: vertex texture coordinates
///
/// Parameters:
/// - mtriangles: element array
/// - mpos: vertex positions
/// - mnorm: vertex normals
/// - mtexcoord: vertex texture coordinates
///
YSHAPE_API void merge_stdsurface(std::vector<int3>& triangles,
    std::vector<float3>& pos, std::vector<float3>& norm,
    std::vector<float2>& texcoord, const std::vector<int3>& mtriangles,
    const std::vector<float3>& mpos, const std::vector<float3>& mnorm,
    const std::vector<float2>& mtexcoord);

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
YSHAPE_API void sample_shape_cdf(int npoints, const int* points, int nverts,
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
YSHAPE_API void sample_shape_cdf(int nlines, const int2* lines, int nverts,
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
YSHAPE_API void sample_shape_cdf(int ntriangles, const int2* triangles,
    int nverts, const float3* pos, float* ecdf, float* length);

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
YSHAPE_API void sample_shape(
    int npoints, const float* point_cdf, float ern, int* eid);

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
YSHAPE_API void sample_shape(int nlines, const float* line_cdf, float ern,
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
YSHAPE_API void sample_shape(int ntriangles, const float* triangle_cdf,
    float ern, const float2& uvrn, int* eid, float2* euv);
}  // namespace

// -----------------------------------------------------------------------------
// INCLUDE FOR HEADER-ONLY MODE
// -----------------------------------------------------------------------------

#ifdef YSHAPE_INLINE
#include "yocto_shape.cpp"
#endif

#endif
