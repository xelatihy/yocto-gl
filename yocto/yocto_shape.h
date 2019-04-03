//
// # Yocto/Shape: Tiny Library for shape operations for graphics
//
//
// Yocto/Shape is a collection of utilities for manipulating shapes in 3D
// graphics, with a focus on triangle and quad meshes. We support both low-level
// geometry operation and whole shape operations.
//
//
// ## Geometry functions
//
// The library supports basic geomtry functions such as computing
// line/triangle/quad normals and areas, picking points on triangles
// and the like. In these functions triangles are parameterized with us written
// w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
// as pairs of two triangles p0,p1,p3 and p2,p3,p1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
// represent triangles correctly, an this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
//
// ## Shape functions
//
// We provide a small number of utilities for shape manipulation for index
// triangle and quad meshes, indexed line and point sets and indexed beziers.
// The utliities collected here are written to support a global illumination
// rendering and not for generic geometry processing. We support operation for
// shape smoothing, shape subdivision (including Catmull-Clark subdivs), and
// example shape creation.
//
// 1. compute line tangents, and triangle and quad areas and normals with
//    `line_tangent()`, `triamgle_normal()`, `quad_normal()` and
//    `line_length()`, `triangle_area()` and `quad_normal()`
// 2. interpolate values over primitives with `interpolate_line()`,
//    `interpolate_triangle()` and `interpolate_quad()`
// 3. evaluate Bezier curves and derivatives with `interpolate_bezier()` and
//    `interpolate_bezier_derivative()`
// 4. compute smooth normals and tangents with `compute_vertex_normals()`
//   `compute_vertex_tangents()`
// 5. compute tangent frames from texture coordinates with
//    `compute_tangent_spaces()`
// 6. compute skinning with `compute_skinning()` and
//    `compute_matrix_skinning()`
// 6. create shapes with `make_box_shape()`, `make_sphere_shape()`,
//    `make_quad_shape()`,  `make_disk_shape()`, `make_box_fvshape()`,
//    `make_hair_shape()`,  `make_suzanne_shape()`, `make_lines_shape()`,
//    `make_points_shape()`,  `make_uvsphere_shape()`,
//    `make_box_rounded_shape()`,  `make_uvsphere_flipcap_shape()`,
//    `make_uvcylinder_shape()`,  `make_uvcylinder_rounded_shape()`,
//    `make_uvdisk_shape()`,  `make_disk_shape()`
// 7. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
// 8. shape sampling with `sample_points_element()`, `sample_lines_element()`,
//    `sample_triangles_element()`; initialize the sampling CDFs with
//    `sample_points_element_cdf()`, `sample_lines_element_cdf()`,
//    `sample_triangles_element_cdf()`
// 9.  sample a could of point over a surface with
// `sample_triangles_element_points()`
// 10. get edges and boundaries with `get_edges()`
// 11. convert quads to triangles with `convert_quads_to_triangles()`
// 12. convert face varying to vertex shared representations with
//     `convert_face_varying()`
// 13. subdivide elements by edge splits with `subdivide_lines()`,
//     `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
// 14. Catmull-Clark subdivision surface with `subdivide_catmullclark()`
//
//
// ## Shape IO
//
// We support reading and writing shapes in OBJ and PLY.
//
// 1. load/save shapes with `load_shape()`/`save_shape()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_SHAPE_H_
#define _YOCTO_SHAPE_H_

#ifndef YOCTO_QUADS_AS_TRIANGLES
#define YOCTO_QUADS_AS_TRIANGLES 1
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"
#include "yocto_random.h"
#include "yocto_utils.h"

#include <unordered_map>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::unordered_map;

}

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
template <typename T>
inline vec<T, 3> line_tangent(const vec<T, 3>& p0, const vec<T, 3>& p1) {
    return normalize(p1 - p0);
}
template <typename T>
inline T line_length(const vec<T, 3>& p0, const vec<T, 3>& p1) {
    return length(p1 - p0);
}

// Triangle properties.
template <typename T>
inline vec<T, 3> triangle_normal(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
    return normalize(cross(p1 - p0, p2 - p0));
}
template <typename T>
inline T triangle_area(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
    return length(cross(p1 - p0, p2 - p0)) / 2;
}

// Quad propeties.
template <typename T>
inline vec<T, 3> quad_normal(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3) {
    return normalize(triangle_normal(p0, p1, p3) + triangle_normal(p2, p3, p1));
}
template <typename T>
inline T quad_area(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3) {
    return triangle_area(p0, p1, p3) + triangle_area(p2, p3, p1);
}

// Triangle tangent and bitangent from uv
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 2>& uv0,
    const vec<T, 2>& uv1, const vec<T, 2>& uv2);

// Quad tangent and bitangent from uv. Note that we pass a current_uv since
// internally we may want to split the quad in two and we need to known where
// to do it. If not interested in the split, just pass zero2f here.
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2,
    const vec<T, 2>& uv3, const vec<T, 2>& current_uv);

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const T& p0, const T& p1, T1 u) {
    return p0 * (1 - u) + p1 * u;
}
// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const T& p0, const T& p1, const T& p2, const vec<T1, 2>& uv) {
    return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T, typename T1>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, const vec<T1, 2>& uv) {
#if YOCTO_QUADS_AS_TRIANGLES
    if (uv.x + uv.y <= 1) {
        return interpolate_triangle(p0, p1, p3, uv);
    } else {
        return interpolate_triangle(p2, p3, p1, 1 - uv);
    }
#else
    return p0 * (1 - uv.x) * (1 - uv.y) + p1 * uv.x * (1 - uv.y) +
           p2 * uv.x * uv.y + p3 * (1 - uv.x) * uv.y;
#endif
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u) {
    return p0 * (1 - u) * (1 - u) * (1 - u) + p1 * 3 * u * (1 - u) * (1 - u) +
           p2 * 3 * u * u * (1 - u) + p3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u) {
    return (p1 - p0) * 3 * (1 - u) * (1 - u) + (p2 - p1) * 6 * u * (1 - u) +
           (p3 - p2) * 3 * u * u;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Load/Save a mesh
void load_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool preserve_facevarying);
void save_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texturecoords, const vector<vec4f>& colors,
    const vector<float>& radius, bool preserve_facevarying, bool ascii = false);

// shapeio error
struct shapeio_error : runtime_error {
    explicit shapeio_error(const char* msg) : runtime_error{msg} {}
    explicit shapeio_error(const std::string& msg) : runtime_error{msg} {}
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex normals/tangents for lines/triangles/quads.
template <typename T>
inline void compute_vertex_tangents(vector<vec<T, 3>>& tangents,
    const vector<vec2i>& lines, const vector<vec<T, 3>>& positions);
template <typename T>
inline void compute_vertex_normals(vector<vec<T, 3>>& normals,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions);
template <typename T>
inline void compute_vertex_normals(vector<vec<T, 3>>& normals,
    const vector<vec4i>& quads, const vector<vec<T, 3>>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
template <typename T>
inline void compute_tangent_spaces(vector<vec<T, 4>>& tangent_spaces,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 2>>& texturecoords);

// Apply skinning to vertex position and normals.
template <typename T>
inline void compute_skinning(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec4i>& joints, const vector<frame<T, 3>>& xforms);
// Apply skinning as specified in Khronos glTF.
template <typename T>
inline void compute_matrix_skinning(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec4i>& joints, const vector<mat<T, 4, 4>>& xforms);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
template <typename T>
inline void flip_vertex_normals(vector<vec<T, 3>>& normals);
// Flip face orientation
template <typename T>
inline void flip_triangles_orientation(vector<vec<T, 3>>& triangles);
template <typename T>
inline void flip_quads_orientation(vector<vec<T, 4>>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
template <typename T>
inline void align_vertices(
    vector<vec<T, 3>>& positions, const vec3i& alignment);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EDGE AND GRID DATA STRUCTURES
// -----------------------------------------------------------------------------
namespace yocto {

// Dictionary to store edge information. edge_index is the index to the edge
// array, edges the array of edges and adj_faces the adjacent faces. We store
// only bidirectional edges to keep the dictionary small. Use the functions
// below to access this data.
struct edge_map {
    unordered_map<vec2i, int> edge_index = {};
    vector<vec2i>             edges      = {};
    vector<int>               num_faces  = {};
};

// Initialize an edge map with elements.
inline void insert_edges(edge_map& emap, const vector<vec3i>& triangles);
inline void insert_edges(edge_map& emap, const vector<vec4i>& quads);
// Insert an edge and return its index
inline int insert_edge(edge_map& emap, const vec2i& edge);
// Get the edge index / insertion count
inline int get_edge_index(const edge_map& emap, const vec2i& edge);
// Get list of edges / boundary edges
inline int  get_num_edges(const edge_map& emap);
inline void get_edges(const edge_map& emap, vector<vec2i>& edges);
inline void get_boundary(const edge_map& emap, vector<vec2i>& edges);
inline void get_edges(const vector<vec3i>& triangles, vector<vec2i>& edges);
inline void get_edges(const vector<vec4i>& quads, vector<vec2i>& edges);

// A sparse grid of cells, containing list of points. Cells are stored in
// a dictionary to get sparsing. Helpful for nearest neighboor lookups.
template <typename T, int N>
struct hash_grid;

template <typename T>
struct hash_grid<T, 3> {
    T                                 cell_size     = 0;
    T                                 cell_inv_size = 0;
    vector<vec<T, 3>>                 positions     = {};
    unordered_map<vec3i, vector<int>> cells         = {};
};

// Create a hash_grid
template <typename T>
inline void init_hash_grid(hash_grid<T, 3>& grid, T cell_size);
template <typename T>
inline void init_hash_grid(
    hash_grid<T, 3>& grid, const vector<vec<T, 3>>& positions, T cell_size);
// Inserts a point into the grid
template <typename T>
inline int insert_vertex(hash_grid<T, 3>& grid, const vec<T, 3>& position);
// Finds the nearest neighboors within a given radius
template <typename T>
inline void find_nearest_neightbors(const hash_grid<T, 3>& grid,
    vector<int>& neighboors, const vec<T, 3>& position, T max_radius);
template <typename T>
inline void find_nearest_neightbors(const hash_grid<T, 3>& grid,
    vector<int>& neighboors, int vertex_id, T max_radius);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
inline void convert_quads_to_triangles(
    vector<vec3i>& triangles, const vector<vec4i>& quads);
// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
inline void convert_quads_to_triangles(
    vector<vec3i>& triangles, const vector<vec4i>& quads, int row_length);
// Convert triangles to quads by creating degenerate quads
inline void convert_triangles_to_quads(
    vector<vec4i>& quads, const vector<vec3i>& triangles);

// Convert beziers to lines using 3 lines for each bezier.
inline void convert_bezier_to_lines(
    vector<vec2i>& lines, const vector<vec4i>& beziers);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm, texcoord and colors.
template <typename T>
inline void convert_facevarying(vector<vec4i>& split_quads,
    vector<vec<T, 3>>& split_positions, vector<vec<T, 3>>& split_normals,
    vector<vec<T, 2>>&   split_texturecoords,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>&     quads_texturecoords,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords);

// Split primitives per id
inline void ungroup_lines(vector<vector<vec2i>>& split_lines,
    const vector<vec2i>& lines, const vector<int>& ids);
inline void ungroup_triangles(vector<vector<vec3i>>& split_triangles,
    const vector<vec3i>& triangles, const vector<int>& ids);
inline void ungroup_quads(vector<vector<vec4i>>& split_quads,
    const vector<vec4i>& quads, const vector<int>& ids);

// Weld vertices within a threshold.
template <typename T>
inline void weld_vertices(
    vector<vec<T, 3>>& positions, vector<int>& indices, T threshold);
template <typename T>
inline void weld_triangles(
    vector<vec3i>& triangles, vector<vec<T, 3>>& positions, T threshold);
template <typename T>
inline void weld_quads(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, T threshold);

// Merge shape elements
inline void merge_lines(
    vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts);
inline void merge_triangles(vector<vec3i>& triangles,
    const vector<vec2i>& merge_triangles, int num_verts);
inline void merge_quads(
    vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts);
template <typename T>
inline void merge_lines(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& tangents, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vector<vec2i>& merge_lines,
    const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_tangents,
    const vector<vec<T, 2>>& merge_texturecoords,
    const vector<T>&         merge_radius);
template <typename T>
inline void merge_triangles(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vector<vec2i>& merge_triangles,
    const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_normals,
    const vector<vec<T, 2>>& merge_texturecoords);
template <typename T>
inline void merge_quads(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vector<vec4i>& merge_quads, const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_normals,
    const vector<vec<T, 2>>& merge_texturecoords);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines by splitting each line in half.
template <typename T>
inline void subdivide_lines(vector<vec2i>& lines, vector<T>& vert);
template <typename T>
inline void subdivide_lines(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<vec<T, 4>>& colors, vector<T>& radius);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
inline void subdivide_triangles(vector<vec3i>& triangles, vector<T>& vert);
template <typename T>
inline void subdivide_triangles(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
inline void subdivide_quads(vector<vec4i>& quads, vector<T>& vert);
template <typename T>
inline void subdivide_quads(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<vec<T, 4>>& colors, vector<T>& radius);
// Subdivide beziers by splitting each segment in two.
template <typename T>
inline void subdivide_beziers(vector<vec4i>& beziers, vector<T>& vert);
template <typename T>
inline void subdivide_beziers(vector<vec4i>& beziers,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
inline void subdivide_catmullclark(
    vector<vec4i>& quads, vector<T>& vert, bool lock_boundary = false);
template <typename T>
inline void subdivide_catmullclark(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
template <typename T>
inline int sample_points_element(int npoints, T re);
template <typename T>
inline void sample_points_element_cdf(vector<T>& cdf, int npoints);
template <typename T>
inline int sample_points_element(const vector<T>& cdf, T re);

// Pick a point on lines uniformly.
template <typename T>
inline void sample_lines_element_cdf(vector<T>& cdf, const vector<vec2i>& lines,
    const vector<vec<T, 3>>& positions);
template <typename T>
inline pair<int, T> sample_lines_element(const vector<T>& cdf, T re, T ru);

// Pick a point on a triangle mesh uniformly.
template <typename T>
inline void sample_triangles_element_cdf(vector<T>& cdf,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions);
template <typename T>
inline pair<int, vec<T, 2>> sample_triangles_element(
    const vector<T>& cdf, T re, const vec<T, 2>& ruv);

// Pick a point on a quad mesh uniformly.
template <typename T>
inline void sample_quads_element_cdf(vector<T>& cdf, const vector<vec4i>& quads,
    const vector<vec<T, 3>>& positions);
template <typename T>
inline pair<int, vec<T, 2>> sample_quads_element(
    const vector<T>& cdf, T re, const vec<T, 2>& ruv);
template <typename T>
inline pair<int, vec<T, 2>> sample_quads_element(const vector<vec4i>& quads,
    const vector<T>& cdf, T re, const vec<T, 2>& ruv);

// Samples a set of points over a triangle/quad mesh uniformly. Returns pos,
// norm and texcoord of the sampled points.
template <typename T>
inline void sample_triangles_points(vector<vec<T, 3>>& sampled_positions,
    vector<vec<T, 3>>&                                 sampled_normals,
    vector<vec<T, 2>>& sampled_texturecoords, const vector<vec3i>& triangles,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords, int npoints, int seed = 7);
template <typename T>
inline void sample_quads_points(vector<vec<T, 3>>& sampled_positions,
    vector<vec<T, 3>>&                             sampled_normals,
    vector<vec<T, 2>>& sampled_texturecoords, const vector<vec4i>& quads,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords, int npoints, int seed = 7);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto {

// Data structure used for geodesic computation
template <typename T>
struct geodesic_solver {
    struct arc_ {
        int node   = 0;
        T   length = 0;
    };
    struct index_ {
        int node  = 0;
        int index = 0;
    };
    vector<vector<arc_>>   graph      = {};
    vector<vector<index_>> edge_index = {};
    vector<vec<T, 3>>      positions  = {};
    vector<vec2i>          edges      = {};
};

// Construct an edge graph
template <typename T>
inline void init_geodesic_solver(geodesic_solver<T>& solver,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions);
template <typename T>
inline void compute_geodesic_distances(geodesic_solver<T>& solver,
    vector<T>& distances, const vector<int>& sources);
template <typename T>
inline void convert_distance_to_color(
    vector<vec<T, 4>>& colors, const vector<T>& distances);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make examples shapes that are not watertight (besides quads).
// Return (triangles, quads, pos, norm, texcoord)
template <typename T>
inline void make_quad_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec2i& steps, const vec<T, 2>& size, const vec<T, 2>& uvsize);
template <typename T>
inline void make_quad_stack_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 2>& uvsize);
template <typename T>
inline void make_floor_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec2i& steps, const vec<T, 2>& size, const vec<T, 2>& uvsize);
template <typename T>
inline void make_floor_bent_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, const vec<T, 2>& size,
    const vec<T, 2>& uvsize, T radius);
template <typename T>
inline void make_box_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec3i& steps, const vec<T, 3>& size, const vec<T, 3>& uvsize);
template <typename T>
inline void make_box_rounded_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 3>& uvsize, T rounded);
template <typename T>
inline void make_uvsphere_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize);
template <typename T>
inline void make_sphere_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize);
template <typename T>
inline void make_uvsphere_flipcap_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize, const vec<T, 2>& zflip);
template <typename T>
inline void make_uvdisk_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize);
template <typename T>
inline void make_disk_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords, int steps,
    T size, T uvsize);
template <typename T>
inline void make_disk_bulged_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize, T height);
template <typename T>
inline void make_quad_bulged_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize, T height);
template <typename T>
inline void make_uvcylinder_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 2>& size,
    const vec<T, 3>& uvsize);
template <typename T>
inline void make_uvcylinder_rounded_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 2>& size,
    const vec<T, 3>& uvsize, T rounded);
template <typename T>
inline void make_geodesic_sphere_shape(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals, int tesselation,
    T size);

// Make examples shapes with are watertight (good for subdivs).
template <typename T>
inline void make_suzanne_shape(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, T size);
template <typename T>
inline void make_box_shape(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, const vec<T, 3>& size);

// Make facevarying example shapes that are watertight (good for subdivs).
template <typename T>
inline void make_box_fvshape(vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 3>& uvsize);
template <typename T>
inline void make_sphere_fvshape(vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
template <typename T>
inline void make_lines_shape(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vec2i& steps, const vec<T, 2>& size,
    const vec<T, 2>& uvsize, const vec<T, 2>& line_radius = {0.001f, 0.001f});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
template <typename T>
inline void make_point_shape(vector<int>& points, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, T point_radius = 0.001f);
template <typename T>
inline void make_points_shape(vector<int>& points, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, int num, T uvsize, T point_radius = 0.001f);
template <typename T>
inline void make_random_points_shape(vector<int>& points,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<T>& radius, int num,
    const vec<T, 3>& size, T uvsize, T point_radius = 0.001f,
    uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
template <typename T>
inline void make_bezier_circle_shape(
    vector<vec4i>& beziers, vector<vec<T, 3>>& positions, T size);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
template <typename T>
inline void make_hair_shape(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vec2i& steps, const vector<vec3i>& striangles,
    const vector<vec4i>& squads, const vector<vec<T, 3>>& spos,
    const vector<vec<T, 3>>& snorm, const vector<vec<T, 2>>& stexcoord,
    const vec<T, 2>& length = {0.1f, 0.1f},
    const vec<T, 2>& rad = {0.001f, 0.001f}, const vec<T, 2>& noise = zero2f,
    const vec<T, 2>& clump = zero2f, const vec<T, 2>& rotation = zero2f,
    int seed = 7);

// Thickens a shape by copy9ing the shape content, rescaling it and flipping its
// normals. Note that this is very much not robust and only useful for trivial
// cases.
template <typename T>
inline void make_shell_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords, T thickness);

// Shape presets used ofr testing.
template <typename T>
inline void make_shape_presets(vector<vec4i>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const string& type);

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Triangle tangent and bitangent from uv
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 2>& uv0,
    const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
    // Follows the definition in http://www.terathon.com/code/tangent.html and
    // https://gist.github.com/aras-p/2843984
    // normal points up from texture space
    auto p   = p1 - p0;
    auto q   = p2 - p0;
    auto s   = vec<T, 2>{uv1.x - uv0.x, uv2.x - uv0.x};
    auto t   = vec<T, 2>{uv1.y - uv0.y, uv2.y - uv0.y};
    auto div = s.x * t.y - s.y * t.x;

    if (div != 0) {
        auto tu = vec<T, 3>{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                      t.y * p.z - t.x * q.z} /
                  div;
        auto tv = vec<T, 3>{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                      s.x * q.z - s.y * p.z} /
                  div;
        return {tu, tv};
    } else {
        return {{1, 0, 0}, {0, 1, 0}};
    }
}

// Quad tangent and bitangent from uv.
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2,
    const vec<T, 2>& uv3, const vec<T, 2>& current_uv) {
#if YOCTO_QUADS_AS_TRIANGLES
    if (current_uv.x + current_uv.y <= 1) {
        return triangle_tangents_fromuv(p0, p1, p3, uv0, uv1, uv3);
    } else {
        return triangle_tangents_fromuv(p2, p3, p1, uv2, uv3, uv1);
    }
#else
    return triangle_tangents_fromuv(p0, p1, p3, uv0, uv1, uv3);
#endif
}

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex tangents for lines.
template <typename T>
inline void compute_vertex_tangents(vector<vec<T, 3>>& tangents,
    const vector<vec2i>& lines, const vector<vec<T, 3>>& positions) {
    if (tangents.size() != positions.size()) {
        throw std::out_of_range("array should be the same length");
    }
    for (auto& tangent : tangents) tangent = zero3f;
    for (auto& l : lines) {
        auto tangent = line_tangent(positions[l.x], positions[l.y]);
        auto length  = line_length(positions[l.x], positions[l.y]);
        tangents[l.x] += tangent * length;
        tangents[l.y] += tangent * length;
    }
    for (auto& tangent : tangents) tangent = normalize(tangent);
}

// Compute per-vertex normals for triangles.
template <typename T>
inline void compute_vertex_normals(vector<vec<T, 3>>& normals,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions) {
    if (normals.size() != positions.size()) {
        throw std::out_of_range("array should be the same length");
    }
    for (auto& normal : normals) normal = zero3f;
    for (auto& t : triangles) {
        auto normal = triangle_normal(
            positions[t.x], positions[t.y], positions[t.z]);
        auto area = triangle_area(
            positions[t.x], positions[t.y], positions[t.z]);
        normals[t.x] += normal * area;
        normals[t.y] += normal * area;
        normals[t.z] += normal * area;
    }
    for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex normals for quads.
template <typename T>
inline void compute_vertex_normals(vector<vec<T, 3>>& normals,
    const vector<vec4i>& quads, const vector<vec<T, 3>>& positions) {
    if (normals.size() != positions.size()) {
        throw std::out_of_range("array should be the same length");
    }
    for (auto& normal : normals) normal = zero3f;
    for (auto& q : quads) {
        auto normal = quad_normal(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        auto area = quad_area(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        normals[q.x] += normal * area;
        normals[q.y] += normal * area;
        normals[q.z] += normal * area;
        if (q.z != q.w) normals[q.w] += normal * area;
    }
    for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
template <typename T>
inline void compute_tangent_spaces(vector<vec<T, 4>>& tangent_spaces,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 2>>& texturecoords) {
    auto tangu = vector<vec<T, 3>>(positions.size(), zero3f);
    auto tangv = vector<vec<T, 3>>(positions.size(), zero3f);
    for (auto t : triangles) {
        auto [tu, tv] = triangle_tangents_fromuv(positions[t.x], positions[t.y],
            positions[t.z], texturecoords[t.x], texturecoords[t.y],
            texturecoords[t.z]);
        tu            = normalize(tu);
        tv            = normalize(tv);
        for (auto vid : {t.x, t.y, t.z}) tangu[vid] += tu;
        for (auto vid : {t.x, t.y, t.z}) tangv[vid] += tv;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    for (auto& tangent : tangent_spaces) tangent = zero4f;
    for (auto i = 0; i < positions.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], normals[i]);
        auto s   = (dot(cross(normals[i], tangu[i]), tangv[i]) < 0) ? -1.0f
                                                                  : 1.0f;
        tangent_spaces[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
}

// Apply skinning
template <typename T>
inline void compute_skinning(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec4i>& joints, const vector<frame<T, 3>>& xforms) {
    if (skinned_positions.size() != positions.size() ||
        skinned_normals.size() != normals.size()) {
        throw std::out_of_range("arrays should be the same size");
    }
    for (auto i = 0; i < positions.size(); i++) {
        skinned_positions[i] =
            transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
    }
    for (auto i = 0; i < normals.size(); i++) {
        skinned_normals[i] = normalize(
            transform_direction(xforms[joints[i].x], normals[i]) *
                weights[i].x +
            transform_direction(xforms[joints[i].y], normals[i]) *
                weights[i].y +
            transform_direction(xforms[joints[i].z], normals[i]) *
                weights[i].z +
            transform_direction(xforms[joints[i].w], normals[i]) *
                weights[i].w);
    }
}

// Apply skinning as specified in Khronos glTF
template <typename T>
inline void compute_matrix_skinning(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec4i>& joints, const vector<mat<T, 4, 4>>& xforms) {
    if (skinned_positions.size() != positions.size() ||
        skinned_normals.size() != normals.size()) {
        throw std::out_of_range("arrays should be the same size");
    }
    for (auto i = 0; i < positions.size(); i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_positions[i] = transform_point(xform, positions[i]);
        skinned_normals[i] = normalize(transform_direction(xform, normals[i]));
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
template <typename T>
inline void flip_vertex_normals(vector<vec<T, 3>>& normals) {
    for (auto& n : normals) n = -n;
}
// Flip face orientation
template <typename T>
inline void flip_triangles_orientation(vector<vec<T, 3>>& triangles) {
    for (auto& t : triangles) swap(t.y, t.z);
}
template <typename T>
inline void flip_quads_orientation(vector<vec<T, 4>>& quads) {
    for (auto& q : quads) swap(q.y, q.w);
}

// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
template <typename T>
inline void align_vertices(
    vector<vec<T, 3>>& positions, const vec3i& alignment) {
    auto bounds = invalid_bbox3f;
    for (auto& p : positions) bounds += p;
    auto offset = vec<T, 3>{0, 0, 0};
    switch (alignment.x) {
        case 1: offset.x = bounds.min.x; break;
        case 2: offset.x = (bounds.min.x + bounds.max.x) / 2; break;
        case 3: offset.x = bounds.max.x; break;
    }
    switch (alignment.y) {
        case 1: offset.y = bounds.min.y; break;
        case 2: offset.y = (bounds.min.y + bounds.max.y) / 2; break;
        case 3: offset.y = bounds.max.y; break;
    }
    switch (alignment.z) {
        case 1: offset.z = bounds.min.z; break;
        case 2: offset.z = (bounds.min.z + bounds.max.z) / 2; break;
        case 3: offset.z = bounds.max.z; break;
    }
    for (auto& p : positions) p -= offset;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF EDGE AND GRID DATA STRUCTURES
// -----------------------------------------------------------------------------
namespace yocto {

// Create key entry for edge_map
inline vec2i make_edgemap_edge(const vec2i& e) {
    return e.x < e.y ? e : vec2i{e.y, e.x};
}

// Initialize an edge map with elements.
inline void insert_edges(edge_map& emap, const vector<vec3i>& triangles) {
    for (int i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        insert_edge(emap, {t.x, t.y});
        insert_edge(emap, {t.y, t.z});
        insert_edge(emap, {t.z, t.x});
    }
}
inline void insert_edges(edge_map& emap, const vector<vec4i>& quads) {
    for (int i = 0; i < quads.size(); i++) {
        auto& q = quads[i];
        insert_edge(emap, {q.x, q.y});
        insert_edge(emap, {q.y, q.z});
        if (q.z != q.w) insert_edge(emap, {q.z, q.w});
        insert_edge(emap, {q.w, q.x});
    }
}
// Insert an edge and return its index
inline int insert_edge(edge_map& emap, const vec2i& edge) {
    auto es = make_edgemap_edge(edge);
    auto it = emap.edge_index.find(es);
    if (it == emap.edge_index.end()) {
        auto idx = (int)emap.edges.size();
        emap.edge_index.insert(it, {es, idx});
        emap.edges.push_back(es);
        emap.num_faces.push_back(1);
        return idx;
    } else {
        auto idx = it->second;
        emap.num_faces[idx] += 1;
        return idx;
    }
}
// Get number of edges
inline int get_num_edges(const edge_map& emap) { return emap.edges.size(); }
// Get the edge index
inline int get_edge_index(const edge_map& emap, const vec2i& edge) {
    auto iterator = emap.edge_index.find(make_edgemap_edge(edge));
    if (iterator == emap.edge_index.end()) return -1;
    return iterator->second;
}
// Get a list of edges, boundary edges, boundary vertices
inline void get_edges(const edge_map& emap, vector<vec2i>& edges) {
    edges = emap.edges;
}
inline void get_boundary(const edge_map& emap, vector<vec2i>& boundary) {
    boundary.clear();
    for (auto edge_index = 0; edge_index < emap.edges.size(); edge_index++) {
        if (emap.num_faces[edge_index] < 2)
            boundary.push_back(emap.edges[edge_index]);
    }
}
inline void get_edges(const vector<vec3i>& triangles, vector<vec2i>& edges) {
    auto emap = edge_map{};
    insert_edges(emap, triangles);
    get_edges(emap, edges);
}
inline void get_edges(const vector<vec4i>& quads, vector<vec2i>& edges) {
    auto emap = edge_map{};
    insert_edges(emap, quads);
    get_edges(emap, edges);
}

// Gets the cell index
template <typename T>
inline vec3i get_cell_index(
    const hash_grid<T, 3>& grid, const vec<T, 3>& position) {
    auto scaledpos = position * grid.cell_inv_size;
    return vec3i{(int)scaledpos.x, (int)scaledpos.y, (int)scaledpos.z};
}

// Create a hash_grid
template <typename T>
inline void init_hash_grid(hash_grid<T, 3>& grid, T cell_size) {
    grid               = {};
    grid.cell_size     = cell_size;
    grid.cell_inv_size = 1 / cell_size;
}
template <typename T>
inline void init_hash_grid(
    hash_grid<T, 3>& grid, const vector<vec<T, 3>>& positions, T cell_size) {
    grid               = {};
    grid.cell_size     = cell_size;
    grid.cell_inv_size = 1 / cell_size;
    for (auto& position : positions) insert_vertex(grid, position);
}
// Inserts a point into the grid
template <typename T>
inline int insert_vertex(hash_grid<T, 3>& grid, const vec<T, 3>& position) {
    auto vertex_id = (int)grid.positions.size();
    auto cell      = get_cell_index(grid, position);
    grid.cells[cell].push_back(vertex_id);
    grid.positions.push_back(position);
    return vertex_id;
}
// Finds the nearest neighboors within a given radius
template <typename T>
inline void find_nearest_neightbors(const hash_grid<T, 3>& grid,
    vector<int>& neighboors, const vec<T, 3>& position, T max_radius,
    int skip_id) {
    auto cell        = get_cell_index(grid, position);
    auto cell_radius = (int)(max_radius * grid.cell_inv_size) + 1;
    neighboors.clear();
    auto max_radius_squared = max_radius * max_radius;
    for (auto k = -cell_radius; k <= cell_radius; k++) {
        for (auto j = -cell_radius; j <= cell_radius; j++) {
            for (auto i = -cell_radius; i <= cell_radius; i++) {
                auto ncell         = cell + vec3i{i, j, k};
                auto cell_iterator = grid.cells.find(ncell);
                if (cell_iterator == grid.cells.end()) continue;
                auto& ncell_vertices = cell_iterator->second;
                for (auto vertex_id : ncell_vertices) {
                    if (distance_squared(grid.positions[vertex_id], position) >
                        max_radius_squared)
                        continue;
                    if (vertex_id == skip_id) continue;
                    neighboors.push_back(vertex_id);
                }
            }
        }
    }
}
template <typename T>
inline void find_nearest_neightbors(const hash_grid<T, 3>& grid,
    vector<int>& neighboors, const vec<T, 3>& position, T max_radius) {
    find_nearest_neightbors(grid, neighboors, position, max_radius, -1);
}
template <typename T>
inline void find_nearest_neightbors(const hash_grid<T, 3>& grid,
    vector<int>& neighboors, int vertex_id, T max_radius) {
    find_nearest_neightbors(
        grid, neighboors, grid.positions[vertex_id], max_radius, vertex_id);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
inline void convert_quads_to_triangles(
    vector<vec3i>& triangles, const vector<vec4i>& quads) {
    triangles.clear();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
}

// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
inline void convert_quads_to_triangles(
    vector<vec3i>& triangles, const vector<vec4i>& quads, int row_length) {
    triangles.clear();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    // triangles.resize(usteps * vsteps * 2);
    // for (auto j = 0; j < vsteps; j++) {
    //     for (auto i = 0; i < usteps; i++) {
    //         auto f1 = triangles[(j * usteps + i) * 2 + 0];
    //         auto f2 = triangles[(j * usteps + i) * 2 + 1];
    //         if ((i + j) % 2) {
    //             f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
    //             f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
    //         } else {
    //             f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
    //             f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
    //         }
    //     }
    // }
}

// Convert triangles to quads by creating degenerate quads
inline void convert_triangles_to_quads(
    vector<vec4i>& quads, const vector<vec3i>& triangles) {
    quads.clear();
    quads.reserve(triangles.size());
    for (auto& t : triangles) quads.push_back({t.x, t.y, t.z, t.z});
}

// Convert beziers to lines using 3 lines for each bezier.
inline void convert_bezier_to_lines(
    vector<vec2i>& lines, const vector<vec4i>& beziers) {
    lines.clear();
    lines.reserve(beziers.size() * 3);
    for (auto b : beziers) {
        lines.push_back({b.x, b.y});
        lines.push_back({b.y, b.z});
        lines.push_back({b.z, b.w});
    }
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
template <typename T>
inline void convert_facevarying(vector<vec4i>& split_quads,
    vector<vec<T, 3>>& split_positions, vector<vec<T, 3>>& split_normals,
    vector<vec<T, 2>>&   split_texturecoords,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>&     quads_texturecoords,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords) {
    // make faces unique
    unordered_map<vec3i, int> vert_map;
    split_quads.resize(quads_positions.size());
    for (auto fid = 0; fid < quads_positions.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec3i{
                (&quads_positions[fid].x)[c],
                (!quads_normals.empty()) ? (&quads_normals[fid].x)[c] : -1,
                (!quads_texturecoords.empty())
                    ? (&quads_texturecoords[fid].x)[c]
                    : -1,
            };
            auto it = vert_map.find(v);
            if (it == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map.insert(it, {v, s});
                (&split_quads[fid].x)[c] = s;
            } else {
                (&split_quads[fid].x)[c] = it->second;
            }
        }
    }

    // fill vert data
    split_positions.clear();
    if (!positions.empty()) {
        split_positions.resize(vert_map.size());
        for (auto& [vert, index] : vert_map) {
            split_positions[index] = positions[vert.x];
        }
    }
    split_normals.clear();
    if (!normals.empty()) {
        split_normals.resize(vert_map.size());
        for (auto& [vert, index] : vert_map) {
            split_normals[index] = normals[vert.y];
        }
    }
    split_texturecoords.clear();
    if (!texturecoords.empty()) {
        split_texturecoords.resize(vert_map.size());
        for (auto& [vert, index] : vert_map) {
            split_texturecoords[index] = texturecoords[vert.z];
        }
    }
}

// Split primitives per id
template <typename T>
inline void ungroup_elems(vector<vector<T>>& split_elems,
    const vector<T>& elems, const vector<int>& ids) {
    auto max_id = *max_element(ids.begin(), ids.end());
    split_elems.resize(max_id + 1);
    for (auto elem_id = 0; elem_id < elems.size(); elem_id++) {
        split_elems[ids[elem_id]].push_back(elems[elem_id]);
    }
}
inline void ungroup_lines(vector<vector<vec2i>>& split_lines,
    const vector<vec2i>& lines, const vector<int>& ids) {
    ungroup_elems(split_lines, lines, ids);
}
inline void ungroup_triangles(vector<vector<vec3i>>& split_triangles,
    const vector<vec3i>& triangles, const vector<int>& ids) {
    ungroup_elems(split_triangles, triangles, ids);
}
inline void ungroup_quads(vector<vector<vec4i>>& split_quads,
    const vector<vec4i>& quads, const vector<int>& ids) {
    ungroup_elems(split_quads, quads, ids);
}

// Weld vertices within a threshold.
template <typename T>
inline void weld_vertices(
    vector<vec<T, 3>>& positions, vector<int>& indices, T threshold) {
    indices.resize(positions.size());
    auto welded_positions = vector<vec<T, 3>>{};
    auto grid             = hash_grid<T, 3>{};
    init_hash_grid(grid, threshold);
    auto neighboors = vector<int>{};
    for (auto vertex_id = 0; vertex_id < positions.size(); vertex_id++) {
        auto& position = positions[vertex_id];
        find_nearest_neightbors(grid, neighboors, position, threshold);
        if (neighboors.empty()) {
            welded_positions.push_back(position);
            indices[vertex_id] = (int)welded_positions.size() - 1;
            insert_vertex(grid, position);
        } else {
            indices[vertex_id] = neighboors.front();
        }
    }
    swap(welded_positions, positions);
    // for (auto i = 0; i < positions.size(); i++) {
    //     welded_indices[i] = (int)welded_positions.size();
    //     for (auto j = 0; j < welded_positions.size(); j++) {
    //         if (length(positions[i] - welded_positions[j]) < threshold) {
    //             welded_indices[i] = j;
    //             break;
    //         }
    //     }
    //     if (welded_indices[i] == (int)welded_positions.size())
    //         welded_positions.push_back(positions[i]);
    // }
}
template <typename T>
inline void weld_triangles(
    vector<vec3i>& triangles, vector<vec<T, 3>>& positions, T threshold) {
    auto indices = vector<int>{};
    weld_vertices(positions, indices, threshold);
    for (auto& t : triangles) t = {indices[t.x], indices[t.y], indices[t.z]};
}
template <typename T>
inline void weld_quads(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, T threshold) {
    auto indices = vector<int>{};
    weld_vertices(positions, indices, threshold);
    auto welded_quads = vector<vec4i>{};
    for (auto& q : quads)
        q = {
            indices[q.x],
            indices[q.y],
            indices[q.z],
            indices[q.w],
        };
}

// Merge shape elements
inline void merge_lines(
    vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts) {
    for (auto& l : merge_lines)
        lines.push_back({l.x + num_verts, l.y + num_verts});
}
inline void merge_triangles(vector<vec3i>& triangles,
    const vector<vec3i>& merge_triangles, int num_verts) {
    for (auto& t : merge_triangles)
        triangles.push_back(
            {t.x + num_verts, t.y + num_verts, t.z + num_verts});
}
inline void merge_quads(
    vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts) {
    for (auto& q : merge_quads)
        quads.push_back({q.x + num_verts, q.y + num_verts, q.z + num_verts,
            q.w + num_verts});
}
template <typename T>
inline void merge_lines(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& tangents, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vector<vec2i>& merge_lines,
    const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_tangents,
    const vector<vec<T, 2>>& merge_texturecoords,
    const vector<T>&         merge_radius) {
    auto merge_verts = (int)positions.size();
    for (auto& l : merge_lines)
        lines.push_back({l.x + merge_verts, l.y + merge_verts});
    positions.insert(
        positions.end(), merge_positions.begin(), merge_positions.end());
    tangents.insert(
        tangents.end(), merge_tangents.begin(), merge_tangents.end());
    texturecoords.insert(texturecoords.end(), merge_texturecoords.begin(),
        merge_texturecoords.end());
    radius.insert(radius.end(), merge_radius.begin(), merge_radius.end());
}
template <typename T>
inline void merge_triangles(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vector<vec3i>& merge_triangles,
    const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_normals,
    const vector<vec<T, 2>>& merge_texturecoords) {
    auto merge_verts = (int)positions.size();
    for (auto& t : merge_triangles)
        triangles.push_back(
            {t.x + merge_verts, t.y + merge_verts, t.z + merge_verts});
    positions.insert(
        positions.end(), merge_positions.begin(), merge_positions.end());
    normals.insert(normals.end(), merge_normals.begin(), merge_normals.end());
    texturecoords.insert(texturecoords.end(), merge_texturecoords.begin(),
        merge_texturecoords.end());
}
template <typename T>
inline void merge_quads(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vector<vec4i>& merge_quads, const vector<vec<T, 3>>& merge_positions,
    const vector<vec<T, 3>>& merge_normals,
    const vector<vec<T, 2>>& merge_texturecoords) {
    auto merge_verts = (int)positions.size();
    for (auto& q : merge_quads)
        quads.push_back({q.x + merge_verts, q.y + merge_verts,
            q.z + merge_verts, q.w + merge_verts});
    positions.insert(
        positions.end(), merge_positions.begin(), merge_positions.end());
    normals.insert(normals.end(), merge_normals.begin(), merge_normals.end());
    texturecoords.insert(texturecoords.end(), merge_texturecoords.begin(),
        merge_texturecoords.end());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines.
template <typename T>
inline void subdivide_lines(vector<vec2i>& lines, vector<T>& vert) {
    // early exit
    if (lines.empty() || vert.empty()) return;
    // sizes
    auto nverts = (int)vert.size();
    auto nlines = (int)lines.size();
    // create vertices
    auto tvert = vector<T>(nverts + nlines);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nlines; i++) {
        auto l            = lines[i];
        tvert[nverts + i] = (vert[l.x] + vert[l.y]) / 2;
    }
    // create lines
    auto tlines = vector<vec2i>(nlines * 2);
    for (auto i = 0; i < nlines; i++) {
        auto l            = lines[i];
        tlines[i * 2 + 0] = {l.x, nverts + i};
        tlines[i * 2 + 0] = {nverts + i, l.y};
    }
    swap(tlines, lines);
    swap(tvert, vert);
}

// Subdivide triangle.
template <typename T>
inline void subdivide_triangles(vector<vec3i>& triangles, vector<T>& vert) {
    // early exit
    if (triangles.empty() || vert.empty()) return;
    // get edges
    auto emap = edge_map{};
    insert_edges(emap, triangles);
    auto edges = vector<vec2i>{};
    get_edges(emap, edges);
    // number of elements
    auto nverts = (int)vert.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)triangles.size();
    // create vertices
    auto tvert = vector<T>(nverts + nedges);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
        auto e            = edges[i];
        tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    // create triangles
    auto ttriangles = vector<vec3i>(nfaces * 4);
    for (auto i = 0; i < nfaces; i++) {
        auto t                = triangles[i];
        ttriangles[i * 4 + 0] = {t.x, nverts + get_edge_index(emap, {t.x, t.y}),
            nverts + get_edge_index(emap, {t.z, t.x})};
        ttriangles[i * 4 + 1] = {t.y, nverts + get_edge_index(emap, {t.y, t.z}),
            nverts + get_edge_index(emap, {t.x, t.y})};
        ttriangles[i * 4 + 2] = {t.z, nverts + get_edge_index(emap, {t.z, t.x}),
            nverts + get_edge_index(emap, {t.y, t.z})};
        ttriangles[i * 4 + 3] = {nverts + get_edge_index(emap, {t.x, t.y}),
            nverts + get_edge_index(emap, {t.y, t.z}),
            nverts + get_edge_index(emap, {t.z, t.x})};
    }
    swap(ttriangles, triangles);
    swap(tvert, vert);
}

// Subdivide quads.
template <typename T>
inline void subdivide_quads(vector<vec4i>& quads, vector<T>& vert) {
    // early exit
    if (quads.empty() || vert.empty()) return;
    // get edges
    auto emap = edge_map{};
    insert_edges(emap, quads);
    auto edges = vector<vec2i>{};
    get_edges(emap, edges);
    // number of elements
    auto nverts = (int)vert.size();
    auto nedges = (int)edges.size();
    auto nfaces = (int)quads.size();
    // create vertices
    auto tvert = vector<T>(nverts + nedges + nfaces);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
        auto e            = edges[i];
        tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    for (auto i = 0; i < nfaces; i++) {
        auto q = quads[i];
        if (q.z != q.w) {
            tvert[nverts + nedges + i] =
                (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4;
        } else {
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.y]) /
                                         3;
        }
    }
    // create quads
    auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
    auto qi     = 0;
    for (auto i = 0; i < nfaces; i++) {
        auto q = quads[i];
        if (q.z != q.w) {
            tquads[qi++] = {q.x, nverts + get_edge_index(emap, {q.x, q.y}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.w, q.x})};
            tquads[qi++] = {q.y, nverts + get_edge_index(emap, {q.y, q.z}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.x, q.y})};
            tquads[qi++] = {q.z, nverts + get_edge_index(emap, {q.z, q.w}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.y, q.z})};
            tquads[qi++] = {q.w, nverts + get_edge_index(emap, {q.w, q.x}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.z, q.w})};
        } else {
            tquads[qi++] = {q.x, nverts + get_edge_index(emap, {q.x, q.y}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.z, q.x})};
            tquads[qi++] = {q.y, nverts + get_edge_index(emap, {q.y, q.z}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.x, q.y})};
            tquads[qi++] = {q.z, nverts + get_edge_index(emap, {q.z, q.x}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.y, q.z})};
        }
    }
    tquads.resize(qi);
    // done
    swap(tquads, quads);
    swap(tvert, vert);
}

// Subdivide beziers.
template <typename T>
inline void subdivide_beziers(vector<vec4i>& beziers, vector<T>& vert) {
    // early exit
    if (beziers.empty() || vert.empty()) return;
    // get edges
    auto vmap     = unordered_map<int, int>();
    auto tvert    = vector<T>();
    auto tbeziers = vector<vec4i>();
    for (auto b : beziers) {
        if (vmap.find(b.x) == vmap.end()) {
            vmap[b.x] = (int)tvert.size();
            tvert.push_back(vert[b.x]);
        }
        if (vmap.find(b.w) == vmap.end()) {
            vmap[b.w] = (int)tvert.size();
            tvert.push_back(vert[b.w]);
        }
        auto bo = (int)tvert.size();
        tbeziers.push_back({vmap.at(b.x), bo + 0, bo + 1, bo + 2});
        tbeziers.push_back({bo + 2, bo + 3, bo + 4, vmap.at(b.w)});
        tvert.push_back(vert[b.x] / 2 + vert[b.y] / 2);
        tvert.push_back(vert[b.x] / 4 + vert[b.y] / 2 + vert[b.z] / 4);
        tvert.push_back(vert[b.x] / 8 + 3 * vert[b.y] / 8 + 3 * vert[b.z] / 8 +
                        vert[b.w] / 8);
        tvert.push_back(vert[b.y] / 4 + vert[b.z] / 2 + vert[b.w] / 4);
        tvert.push_back(vert[b.z] / 2 + vert[b.w] / 2);
    }

    // done
    swap(tbeziers, beziers);
    swap(tvert, vert);
}

// Subdivide catmullclark.
template <typename T>
inline void subdivide_catmullclark(
    vector<vec4i>& quads, vector<T>& vert, bool lock_boundary) {
    // early exit
    if (quads.empty() || vert.empty()) return;
    // get edges
    auto emap = edge_map{};
    insert_edges(emap, quads);
    auto edges = vector<vec2i>{}, boundary = vector<vec2i>{};
    get_edges(emap, edges);
    get_boundary(emap, boundary);
    // number of elements
    auto nverts    = (int)vert.size();
    auto nedges    = (int)edges.size();
    auto nboundary = (int)boundary.size();
    auto nfaces    = (int)quads.size();

    // split elements ------------------------------------
    // create vertices
    auto tvert = vector<T>(nverts + nedges + nfaces);
    for (auto i = 0; i < nverts; i++) tvert[i] = vert[i];
    for (auto i = 0; i < nedges; i++) {
        auto e            = edges[i];
        tvert[nverts + i] = (vert[e.x] + vert[e.y]) / 2;
    }
    for (auto i = 0; i < nfaces; i++) {
        auto q = quads[i];
        if (q.z != q.w) {
            tvert[nverts + nedges + i] =
                (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4;
        } else {
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.y]) /
                                         3;
        }
    }
    // create quads
    auto tquads = vector<vec4i>(nfaces * 4);  // conservative allocation
    auto qi     = 0;
    for (auto i = 0; i < nfaces; i++) {
        auto q = quads[i];
        if (q.z != q.w) {
            tquads[qi++] = {q.x, nverts + get_edge_index(emap, {q.x, q.y}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.w, q.x})};
            tquads[qi++] = {q.y, nverts + get_edge_index(emap, {q.y, q.z}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.x, q.y})};
            tquads[qi++] = {q.z, nverts + get_edge_index(emap, {q.z, q.w}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.y, q.z})};
            tquads[qi++] = {q.w, nverts + get_edge_index(emap, {q.w, q.x}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.z, q.w})};
        } else {
            tquads[qi++] = {q.x, nverts + get_edge_index(emap, {q.x, q.y}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.z, q.x})};
            tquads[qi++] = {q.y, nverts + get_edge_index(emap, {q.y, q.z}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.x, q.y})};
            tquads[qi++] = {q.z, nverts + get_edge_index(emap, {q.z, q.x}),
                nverts + nedges + i, nverts + get_edge_index(emap, {q.y, q.z})};
        }
    }
    tquads.resize(qi);

    // split boundary
    auto tboundary = vector<vec2i>(nboundary * 2);
    for (auto i = 0; i < nboundary; i++) {
        auto e = boundary[i];
        tboundary.push_back({e.x, nverts + get_edge_index(emap, e)});
        tboundary.push_back({nverts + get_edge_index(emap, e), e.y});
    }

    // setup creases -----------------------------------
    auto tcrease_edges = vector<vec2i>();
    auto tcrease_verts = vector<int>();
    if (lock_boundary) {
        for (auto& b : tboundary) {
            tcrease_verts.push_back(b.x);
            tcrease_verts.push_back(b.y);
        }
    } else {
        for (auto& b : tboundary) tcrease_edges.push_back(b);
    }

    // define vertex valence ---------------------------
    auto tvert_val = vector<int>(tvert.size(), 2);
    for (auto& e : tboundary) {
        tvert_val[e.x] = (lock_boundary) ? 0 : 1;
        tvert_val[e.y] = (lock_boundary) ? 0 : 1;
    }

    // averaging pass ----------------------------------
    auto avert  = vector<T>(tvert.size(), T());
    auto acount = vector<int>(tvert.size(), 0);
    for (auto p : tcrease_verts) {
        if (tvert_val[p] != 0) continue;
        avert[p] += tvert[p];
        acount[p] += 1;
    }
    for (auto& e : tcrease_edges) {
        auto c = (tvert[e.x] + tvert[e.y]) / 2.0f;
        for (auto vid : {e.x, e.y}) {
            if (tvert_val[vid] != 1) continue;
            avert[vid] += c;
            acount[vid] += 1;
        }
    }
    for (auto& q : tquads) {
        auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4.0f;
        for (auto vid : {q.x, q.y, q.z, q.w}) {
            if (tvert_val[vid] != 2) continue;
            avert[vid] += c;
            acount[vid] += 1;
        }
    }
    for (auto i = 0; i < tvert.size(); i++) avert[i] /= (T)acount[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i = 0; i < tvert.size(); i++) {
        if (tvert_val[i] != 2) continue;
        avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4.0f / acount[i]);
    }
    tvert = avert;

    // done
    swap(tquads, quads);
    swap(tvert, vert);
}

template <typename ET, typename T, typename SubdivideFunc>
inline void subdivide_elems_impl(vector<ET>& elems,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius, const SubdivideFunc&& subdivied_func) {
    if (!normals.empty()) {
        auto elems_ = elems;
        subdivied_func(elems_, normals);
    }
    if (!texturecoords.empty()) {
        auto elems_ = elems;
        subdivied_func(elems_, texturecoords);
    }
    if (!colors.empty()) {
        auto elems_ = elems;
        subdivied_func(elems_, colors);
    }
    if (!radius.empty()) {
        auto elems_ = elems;
        subdivied_func(elems_, radius);
    }
    subdivied_func(elems, positions);
}

template <typename T>
inline void subdivide_lines(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<vec<T, 4>>& colors, vector<T>& radius) {
    subdivide_elems_impl(lines, positions, normals, texturecoords, colors,
        radius, [](auto& lines, auto& vert) { subdivide_lines(lines, vert); });
}
template <typename T>
inline void subdivide_triangles(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius) {
    subdivide_elems_impl(triangles, positions, normals, texturecoords, colors,
        radius, [](auto& triangles, auto& vert) {
            subdivide_triangles(triangles, vert);
        });
}
template <typename T>
inline void subdivide_quads(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<vec<T, 4>>& colors, vector<T>& radius) {
    subdivide_elems_impl(quads, positions, normals, texturecoords, colors,
        radius, [](auto& quads, auto& vert) { subdivide_quads(quads, vert); });
}
template <typename T>
inline void subdivide_beziers(vector<vec4i>& beziers,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius) {
    subdivide_elems_impl(beziers, positions, normals, texturecoords, colors,
        radius,
        [](auto& beziers, auto& vert) { subdivide_beziers(beziers, vert); });
}
template <typename T>
inline void subdivide_catmullclark(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<vec<T, 4>>& colors,
    vector<T>& radius) {
    subdivide_elems_impl(quads, positions, normals, texturecoords, colors,
        radius,
        [](auto& quads, auto& vert) { subdivide_catmullclark(quads, vert); });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
template <typename T>
inline int sample_points_element(int npoints, T re) {
    return sample_uniform_index(npoints, re);
}
template <typename T>
inline void sample_points_element_cdf(vector<T>& cdf, int npoints) {
    cdf.resize(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
}
template <typename T>
inline int sample_points_element(const vector<T>& cdf, T re) {
    return sample_discrete_distribution(cdf, re);
}

// Pick a point on lines uniformly.
template <typename T>
inline void sample_lines_element_cdf(vector<T>& cdf, const vector<vec2i>& lines,
    const vector<vec<T, 3>>& positions) {
    cdf.resize(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(positions[l.x], positions[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
}
template <typename T>
inline pair<int, T> sample_lines_element(const vector<T>& cdf, T re, T ru) {
    return {sample_discrete_distribution(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
template <typename T>
inline void sample_triangles_element_cdf(vector<T>& cdf,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions) {
    cdf.resize(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
}
template <typename T>
inline pair<int, vec<T, 2>> sample_triangles_element(
    const vector<T>& cdf, T re, const vec<T, 2>& ruv) {
    return {sample_discrete_distribution(cdf, re),
        sample_triangle_coordinates(ruv)};
}

// Pick a point on a quad mesh uniformly.
template <typename T>
inline void sample_quads_element_cdf(vector<T>& cdf, const vector<vec4i>& quads,
    const vector<vec<T, 3>>& positions) {
    cdf.resize(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
}
template <typename T>
inline pair<int, vec<T, 2>> sample_quads_element(
    const vector<T>& cdf, T re, const vec<T, 2>& ruv) {
    return {sample_discrete_distribution(cdf, re), ruv};
}
template <typename T>
inline pair<int, vec<T, 2>> sample_quads_element(const vector<vec4i>& quads,
    const vector<T>& cdf, T re, const vec<T, 2>& ruv) {
    auto element_id = sample_discrete_distribution(cdf, re);
    if (quads[element_id].z == quads[element_id].w) {
        return {element_id, sample_triangle_coordinates(ruv)};
    } else {
        return {element_id, ruv};
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec<T, 3> numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
template <typename T>
inline void sample_triangles_points(vector<vec<T, 3>>& sampled_positions,
    vector<vec<T, 3>>&                                 sampled_normals,
    vector<vec<T, 2>>& sampled_texturecoords, const vector<vec3i>& triangles,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords, int npoints, int seed) {
    sampled_positions.resize(npoints);
    sampled_normals.resize(npoints);
    sampled_texturecoords.resize(npoints);
    auto cdf = vector<T>{};
    sample_triangles_element_cdf(cdf, triangles, positions);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto [triangle_id, triangle_uv] = sample_triangles_element(cdf,
            (T)get_random_float(rng),
            {(T)get_random_float(rng), (T)get_random_float(rng)});
        auto t                          = triangles[triangle_id];
        sampled_positions[i]            = interpolate_triangle(
            positions[t.x], positions[t.y], positions[t.z], triangle_uv);
        if (!sampled_normals.empty()) {
            sampled_normals[i] = normalize(interpolate_triangle(
                normals[t.x], normals[t.y], normals[t.z], triangle_uv));
        } else {
            sampled_normals[i] = triangle_normal(
                positions[t.x], positions[t.y], positions[t.z]);
        }
        if (!sampled_texturecoords.empty()) {
            sampled_texturecoords[i] = interpolate_triangle(texturecoords[t.x],
                texturecoords[t.y], texturecoords[t.z], triangle_uv);
        } else {
            sampled_texturecoords[i] = zero2f;
        }
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec<T, 3> numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
template <typename T>
inline void sample_quads_points(vector<vec<T, 3>>& sampled_positions,
    vector<vec<T, 3>>&                             sampled_normals,
    vector<vec<T, 2>>& sampled_texturecoords, const vector<vec4i>& quads,
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 2>>& texturecoords, int npoints, int seed) {
    sampled_positions.resize(npoints);
    sampled_normals.resize(npoints);
    sampled_texturecoords.resize(npoints);
    auto cdf = vector<T>{};
    sample_quads_element_cdf(cdf, quads, positions);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto [quad_id, quad_uv] = sample_quads_element(cdf,
            (T)get_random_float(rng),
            {(T)get_random_float(rng), (T)get_random_float(rng)});
        auto q                  = quads[quad_id];
        sampled_positions[i] = interpolate_quad(positions[q.x], positions[q.y],
            positions[q.z], positions[q.w], quad_uv);
        if (!sampled_normals.empty()) {
            sampled_normals[i] = normalize(interpolate_quad(normals[q.x],
                normals[q.y], normals[q.z], normals[q.w], quad_uv));
        } else {
            sampled_normals[i] = quad_normal(
                positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        }
        if (!sampled_texturecoords.empty()) {
            sampled_texturecoords[i] = interpolate_quad(texturecoords[q.x],
                texturecoords[q.y], texturecoords[q.z], texturecoords[q.w],
                quad_uv);
        } else {
            sampled_texturecoords[i] = zero2f;
        }
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline void add_node(geodesic_solver<T>& solver, const vec<T, 3>& position) {
    solver.positions.push_back(position);
    solver.graph.push_back({});
    solver.graph.back().reserve(8);
}

template <typename T>
inline void add_directed_arc(geodesic_solver<T>& solver, int from, int to) {
    // assert(from >= 0 && from < solver.graph.size());
    // assert(to >= 0 && to < solver.graph.size());
    T len = length(solver.positions[from] - solver.positions[to]);
    solver.graph[from].push_back({to, len});
}

template <typename T>
inline void add_undirected_arc(geodesic_solver<T>& solver, int na, int nb) {
    add_directed_arc(solver, na, nb);
    add_directed_arc(solver, nb, na);
}

template <typename T>
inline void make_edge_solver_slow(geodesic_solver<T>& solver,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions,
    bool use_steiner_points) {
    solver.graph.reserve(positions.size());

    for (int i = 0; i < positions.size(); i++) add_node(solver, positions[i]);

    auto emap = edge_map{};
    insert_edges(emap, triangles);
    auto edges = vector<vec2i>{};
    get_edges(emap, edges);
    for (auto& edge : edges) {
        add_undirected_arc(solver, edge.x, edge.y);
    }

    if (!use_steiner_points) return;

    solver.graph.reserve(size(positions) + size(edges));
    auto steiner_per_edge = vector<int>(get_num_edges(emap));

    // On each edge, connect the mid vertex with the vertices on th same edge.
    for (auto edge_index = 0; edge_index < size(edges); edge_index++) {
        auto& edge                   = edges[edge_index];
        auto  steiner_idx            = solver.graph.size();
        steiner_per_edge[edge_index] = steiner_idx;
        add_node(solver, (positions[edge.x] + positions[edge.y]) * 0.5f);
        add_directed_arc(solver, steiner_idx, edge.x);
        add_directed_arc(solver, steiner_idx, edge.y);
    }

    // Make connection for each face
    for (int face = 0; face < triangles.size(); ++face) {
        int steiner_idx[3];
        for (int k : {0, 1, 2}) {
            int a          = triangles[face][k];
            int b          = triangles[face][(k + 1) % 3];
            steiner_idx[k] = steiner_per_edge[get_edge_index(emap, {a, b})];
        }

        // Connect each mid-vertex to the opposite mesh vertex in the triangle
        int opp[3] = {triangles[face].z, triangles[face].x, triangles[face].y};
        add_undirected_arc(solver, steiner_idx[0], opp[0]);
        add_undirected_arc(solver, steiner_idx[1], opp[1]);
        add_undirected_arc(solver, steiner_idx[2], opp[2]);

        // Connect mid-verts of the face between them
        add_undirected_arc(solver, steiner_idx[0], steiner_idx[1]);
        add_undirected_arc(solver, steiner_idx[0], steiner_idx[2]);
        add_undirected_arc(solver, steiner_idx[1], steiner_idx[2]);
    }
}

template <typename T>
inline void add_half_edge(geodesic_solver<T>& solver, const vec2i& edge) {
    // check if edge exists already
    for (auto [vert, _] : solver.graph[edge.x]) {
        if (vert == edge.y) return;
    }
    auto len = length(solver.positions[edge.x] - solver.positions[edge.y]);
    auto edge_index = (int)solver.edges.size();
    solver.graph[edge.x].push_back({edge.y, len});
    solver.edge_index[edge.x].push_back({edge.y, edge_index});
    solver.edges.push_back(edge);
}

template <typename T>
inline void add_edge(geodesic_solver<T>& solver, const vec2i& edge) {
    // check if edge exists already
    for (auto [vert, _] : solver.graph[edge.x]) {
        if (vert == edge.y) return;
    }
    auto len = length(solver.positions[edge.x] - solver.positions[edge.y]);
    solver.graph[edge.x].push_back({edge.y, len});
    solver.graph[edge.y].push_back({edge.x, len});
    auto edge_index = (int)solver.edges.size();
    solver.edge_index[edge.x].push_back({edge.y, edge_index});
    solver.edge_index[edge.y].push_back({edge.x, edge_index});
    solver.edges.push_back(edge);
}

template <typename T>
inline int get_edge_index(const geodesic_solver<T>& solver, const vec2i& edge) {
    for (auto [node, index] : solver.edge_index[edge.x])
        if (edge.y == node) return index;
    return -1;
}

template <typename T>
inline void make_edge_solver_fast(geodesic_solver<T>& solver,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions,
    bool use_steiner_points) {
    solver.positions = positions;
    solver.graph.resize(size(positions));
    solver.edge_index.resize(size(positions));

    // fast construction assuming edges are not repeated
    for (auto t : triangles) {
        add_edge(solver, {t.x, t.y});
        add_edge(solver, {t.y, t.z});
        add_edge(solver, {t.z, t.x});
    }

    if (!use_steiner_points) return;

    auto edges = solver.edges;
    solver.graph.resize(size(positions) + size(edges));
    solver.edge_index.resize(size(positions) + size(edges));
    solver.positions.resize(size(positions) + size(edges));
    auto steiner_per_edge = vector<int>(size(edges));

    // On each edge, connect the mid vertex with the vertices on th same edge.
    auto edge_offset = (int)positions.size();
    for (auto edge_index = 0; edge_index < size(edges); edge_index++) {
        auto& edge                   = edges[edge_index];
        auto  steiner_idx            = edge_offset + edge_index;
        steiner_per_edge[edge_index] = steiner_idx;
        solver.positions[steiner_idx] =
            (positions[edge.x] + positions[edge.y]) * 0.5f;
        add_half_edge(solver, {steiner_idx, edge.x});
        add_half_edge(solver, {steiner_idx, edge.y});
    }

    // Make connection for each face
    for (int face = 0; face < triangles.size(); ++face) {
        int steiner_idx[3];
        for (int k : {0, 1, 2}) {
            int a          = triangles[face][k];
            int b          = triangles[face][(k + 1) % 3];
            steiner_idx[k] = steiner_per_edge[get_edge_index(solver, {a, b})];
        }

        // Connect each mid-vertex to the opposite mesh vertex in the triangle
        int opp[3] = {triangles[face].z, triangles[face].x, triangles[face].y};
        add_edge(solver, {steiner_idx[0], opp[0]});
        add_edge(solver, {steiner_idx[1], opp[1]});
        add_edge(solver, {steiner_idx[2], opp[2]});

        // Connect mid-verts of the face between them
        add_edge(solver, {steiner_idx[0], steiner_idx[1]});
        add_edge(solver, {steiner_idx[1], steiner_idx[2]});
        add_edge(solver, {steiner_idx[2], steiner_idx[0]});
    }
}

template <typename T>
inline void log_geodesic_solver_stats(const geodesic_solver<T>& solver) {
    // stats
    auto num_edges     = 0;
    auto min_adjacents = int_max, max_adjacents = int_min;
    auto min_length = type_max<T>, max_length = type_min<T>;
    auto avg_adjacents = 0.0, avg_length = 0.0;
    for (auto& adj : solver.graph) {
        num_edges += (int)adj.size();
        min_adjacents = min(min_adjacents, (int)adj.size());
        max_adjacents = max(max_adjacents, (int)adj.size());
        avg_adjacents += adj.size() / (double)solver.graph.size();
        for (auto& edge : adj) {
            min_length = min(min_length, edge.length);
            max_length = max(max_length, edge.length);
            avg_length += edge.length;
        }
    }
    avg_length /= num_edges;
}

template <typename T>
inline void update_edge_distances(geodesic_solver<T>& solver) {
    for (auto node = 0; node < solver.graph.size(); node++) {
        for (auto& edge : solver.graph[node])
            edge.length = length(
                solver.positions[node] - solver.positions[edge.node]);
    }
}

template <typename T>
inline void init_geodesic_solver(geodesic_solver<T>& solver,
    const vector<vec3i>& triangles, const vector<vec<T, 3>>& positions) {
    make_edge_solver_fast(solver, triangles, positions, true);
    // auto solver = make_edge_solver_slow(triangles, positions, true);
    log_geodesic_solver_stats(solver);
}

template <typename T>
inline void compute_geodesic_distances(geodesic_solver<T>& graph,
    vector<T>& distances, const vector<int>& sources) {
    // preallocated
    distances.resize(graph.positions.size());
    for (auto& d : distances) d = type_max<T>;

    // Small Label Fisrt + Large Label Last
    // https://en.wikipedia.org/wiki/Shortest_Path_Faster_Algorithm
    auto num_nodes = (int)graph.graph.size();
    // assert(distances.size() == num_nodes);
    auto visited = vector<bool>(num_nodes, false);

    // setup queue
    auto queue = deque<int>{};
    for (auto source : sources) {
        distances[source] = 0.0f;
        visited[source]   = true;
        queue.push_back(source);
    }

    // Cumulative weights of elements in queue.
    auto cumulative_weight = 0.0f;
    while (!queue.empty()) {
        auto node           = queue.front();
        auto average_weight = (double)cumulative_weight / queue.size();

        // Large Label Last: until front node weights more than the average, put
        // it on back. Sometimes average_weight is less than envery value due to
        // Ting point errors  (doesn't happen with double precision).
        for (auto tries = 0; tries < queue.size() + 1; tries++) {
            if (distances[node] <= average_weight) break;
            queue.pop_front();
            queue.push_back(node);
            node = queue.front();
        }
        queue.pop_front();
        visited[node] = false;                 // out of queue
        cumulative_weight -= distances[node];  // update average

        const auto offset_distance = distances[node];
        const auto num_neighbors   = (int)graph.graph[node].size();

        for (int neighbor_idx = 0; neighbor_idx < num_neighbors;
             neighbor_idx++) {
            // distance and id to neightbor through this node
            auto new_distance = offset_distance +
                                graph.graph[node][neighbor_idx].length;
            auto neighbor = graph.graph[node][neighbor_idx].node;

            auto old_distance = distances[neighbor];
            if (new_distance >= old_distance) continue;

            if (visited[neighbor]) {
                // if neighbor already in queue, update cumulative weights.
                cumulative_weight = cumulative_weight - old_distance +
                                    new_distance;
            } else {
                // if neighbor not in queue, Small Label first.
                if (queue.empty() || (new_distance < distances[queue.front()]))
                    queue.push_front(neighbor);
                else
                    queue.push_back(neighbor);

                visited[neighbor] = true;
                cumulative_weight = cumulative_weight + new_distance;
            }

            distances[neighbor] = new_distance;
        }
    }
}

template <typename T>
inline void convert_distance_to_color(
    vector<vec<T, 4>>& colors, const vector<T>& distances) {
    colors.resize(distances.size());
    for (auto idx = 0; idx < distances.size(); idx++) {
        auto distance = fmod(distances[idx] * 10, (T)1);
        colors[idx]   = {distance, distance, distance, (T)1};
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a quad.
template <typename T>
inline void make_quad_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec2i& steps, const vec<T, 2>& size, const vec<T, 2>& uvsize) {
    positions.resize((steps.x + 1) * (steps.y + 1));
    normals.resize((steps.x + 1) * (steps.y + 1));
    texturecoords.resize((steps.x + 1) * (steps.y + 1));
    for (auto j = 0; j <= steps.y; j++) {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec<T, 2>{i / (T)steps.x, j / (T)steps.y};
            positions[j * (steps.x + 1) + i] = {
                (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
            normals[j * (steps.x + 1) + i]       = {0, 0, 1};
            texturecoords[j * (steps.x + 1) + i] = uv * uvsize;
        }
    }

    for (auto& uv : texturecoords) uv.y = 1 - uv.y;

    quads.resize(steps.x * steps.y);
    for (auto j = 0; j < steps.y; j++) {
        for (auto i = 0; i < steps.x; i++) {
            quads[j * steps.x + i] = {j * (steps.x + 1) + i,
                j * (steps.x + 1) + i + 1, (j + 1) * (steps.x + 1) + i + 1,
                (j + 1) * (steps.x + 1) + i};
        }
    }
    // shape.triangles.resize(steps.x * steps.y * 2);
    // for (auto j = 0; j < steps.y; j++) {
    //     for (auto i = 0; i < steps.x; i++) {
    //         shape.triangles[(j * steps.x + i) * 2 + 0] = {
    //             j * (steps.x + 1) + i, j * (steps.x + 1) + i + 1,
    //             (j + 1) * (steps.x + 1) + i + 1};
    //         shape.triangles[(j * steps.x + i) * 2 + 1] = {
    //             j * (steps.x + 1) + i, (j + 1) * (steps.x + 1) + i + 1,
    //             (j + 1) * (steps.x + 1) + i};
    //     }
    // }
}

template <typename T>
inline void make_floor_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec2i& steps, const vec<T, 2>& size, const vec<T, 2>& uvsize) {
    make_quad_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    for (auto& p : positions) p = {p.x, p.z, p.y};
    for (auto& normal : normals) normal = {normal.x, normal.z, normal.y};
    for (auto& q : quads) swap(q.y, q.w);
}

// Make a rounded cube.
template <typename T>
inline void make_floor_bent_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, const vec<T, 2>& size,
    const vec<T, 2>& uvsize, T radius) {
    make_floor_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    auto start = (size.y / 2 - radius) / 2;
    auto end   = start + radius;
    for (auto i = 0; i < positions.size(); i++) {
        if (positions[i].z < -end) {
            positions[i] = {
                positions[i].x, -positions[i].z - end + radius, -end};
            normals[i] = {0, 0, 1};
        } else if (positions[i].z < -start && positions[i].z >= -end) {
            auto phi     = (pif / 2) * (-positions[i].z - start) / radius;
            positions[i] = {positions[i].x, -cos(phi) * radius + radius,
                -sin(phi) * radius - start};
            normals[i]   = {0, cos(phi), sin(phi)};
        } else {
        }
    }
}

// Make a stack of quads
template <typename T>
inline void make_quad_stack_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 2>& uvsize) {
    quads.clear();
    positions.clear();
    normals.clear();
    texturecoords.clear();
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec<T, 3>>{};
    auto qnormals       = vector<vec<T, 3>>{};
    auto qtexturecoords = vector<vec<T, 2>>{};
    for (auto i = 0; i <= steps.z; i++) {
        make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
            {steps.x, steps.y}, {size.x, size.y}, uvsize);
        for (auto& p : qpositions) p.z = (-0.5f + (T)i / steps.z) * size.z;
        merge_quads(quads, positions, normals, texturecoords, qquads,
            qpositions, qnormals, qtexturecoords);
    }
}

// Make a cube.
template <typename T>
inline void make_box_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    const vec3i& steps, const vec<T, 3>& size, const vec<T, 3>& uvsize) {
    quads.clear();
    positions.clear();
    normals.clear();
    texturecoords.clear();
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec<T, 3>>{};
    auto qnormals       = vector<vec<T, 3>>{};
    auto qtexturecoords = vector<vec<T, 2>>{};
    // + z
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto& p : qpositions) p = {p.x, p.y, size.z / 2};
    for (auto& n : qnormals) n = {0, 0, 1};
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - z
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto& p : qpositions) p = {-p.x, p.y, -size.z / 2};
    for (auto& n : qnormals) n = {0, 0, -1};
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // + x
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.z, steps.y}, {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto& p : qpositions) p = {size.x / 2, p.y, -p.x};
    for (auto& n : qnormals) n = {1, 0, 0};
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - x
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.z, steps.y}, {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto& p : qpositions) p = {-size.x / 2, p.y, p.x};
    for (auto& n : qnormals) n = {-1, 0, 0};
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // + y
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.z}, {size.x, size.z}, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {qpositions[i].x, size.y / 2, -qpositions[i].y};
        qnormals[i]   = {0, 1, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - y
    make_quad_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.z}, {size.x, size.z}, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {qpositions[i].x, -size.y / 2, qpositions[i].y};
        qnormals[i]   = {0, -1, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
}

// Make a rounded cube.
template <typename T>
inline void make_box_rounded_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 3>& uvsize, T rounded) {
    make_box_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    auto radius = rounded * min(size) / 2;
    auto c      = size / 2 - vec<T, 3>{radius, radius, radius};
    for (auto i = 0; i < positions.size(); i++) {
        auto pc = vec<T, 3>{
            fabs(positions[i].x), fabs(positions[i].y), fabs(positions[i].z)};
        auto ps = vec<T, 3>{positions[i].x < 0 ? -1.0f : 1.0f,
            positions[i].y < 0 ? -1.0f : 1.0f,
            positions[i].z < 0 ? -1.0f : 1.0f};
        if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
            auto pn      = normalize(pc - c);
            positions[i] = c + radius * pn;
            normals[i]   = pn;
        } else if (pc.x >= c.x && pc.y >= c.y) {
            auto pn      = normalize((pc - c) * vec<T, 3>{1, 1, 0});
            positions[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
            normals[i]   = pn;
        } else if (pc.x >= c.x && pc.z >= c.z) {
            auto pn      = normalize((pc - c) * vec<T, 3>{1, 0, 1});
            positions[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
            normals[i]   = pn;
        } else if (pc.y >= c.y && pc.z >= c.z) {
            auto pn      = normalize((pc - c) * vec<T, 3>{0, 1, 1});
            positions[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
            normals[i]   = pn;
        } else {
            continue;
        }
        positions[i] *= ps;
        normals[i] *= ps;
    }
}

// Make a sphere.
template <typename T>
inline void make_uvsphere_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize) {
    make_quad_shape(
        quads, positions, normals, texturecoords, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv = texturecoords[i];
        auto a  = vec<T, 2>{2 * pif * uv.x, pif * (1 - uv.y)};
        auto p  = vec<T, 3>{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        positions[i]     = p * (size / 2);
        normals[i]       = normalize(p);
        texturecoords[i] = uv * uvsize;
    }
}

// Make a spherecube.
template <typename T>
inline void make_sphere_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize) {
    make_box_shape(quads, positions, normals, texturecoords,
        {steps, steps, steps}, {1, 1, 1}, {uvsize, uvsize, uvsize});
    for (auto i = 0; i < positions.size(); i++) {
        auto p       = positions[i];
        positions[i] = normalize(p) * (size / 2);
        normals[i]   = normalize(p);
    }
}

// Make a flipped sphere. This is not watertight.
template <typename T>
inline void make_uvsphere_flipcap_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize, const vec<T, 2>& zflip_) {
    make_uvsphere_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    auto zflip = zflip_ * size / 2;
    for (auto i = 0; i < positions.size(); i++) {
        if (positions[i].z > zflip.y) {
            positions[i].z = 2 * zflip.y - positions[i].z;
            normals[i].x   = -normals[i].x;
            normals[i].y   = -normals[i].y;
        } else if (positions[i].z < zflip.x) {
            positions[i].z = 2 * zflip.x - positions[i].z;
            normals[i].x   = -normals[i].x;
            normals[i].y   = -normals[i].y;
        }
    }
}

// Make a disk.
template <typename T>
inline void make_uvdisk_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, T size,
    const vec<T, 2>& uvsize) {
    make_quad_shape(
        quads, positions, normals, texturecoords, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv      = texturecoords[i];
        auto phi     = 2 * pif * uv.x;
        positions[i] = {
            cos(phi) * uv.y * size / 2, sin(phi) * uv.y * size / 2, 0};
        normals[i]       = {0, 0, 1};
        texturecoords[i] = uv * uvsize;
    }
}

// Make a disk from a quad.
template <typename T>
inline void make_disk_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords, int steps,
    T size, T uvsize) {
    make_quad_shape(quads, positions, normals, texturecoords, {steps, steps},
        {2, 2}, {uvsize, uvsize});
    for (auto i = 0; i < positions.size(); i++) {
        // Analytical Methods for Squaring the Disc, by C. Fong
        // https://arxiv.org/abs/1509.06344
        auto xy = vec<T, 2>{positions[i].x, positions[i].y};
        auto uv = vec<T, 2>{
            xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
        positions[i] = {uv.x * size / 2, uv.y * size / 2, 0};
    }
}

// Make a bulged disk from a quad.
template <typename T>
inline void make_disk_bulged_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize, T height) {
    make_disk_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    if (height == 0) return;
    auto radius = (size * size / 4 + height * height) / (2 * height);
    auto center = vec<T, 3>{0, 0, -radius + height};
    for (auto i = 0; i < positions.size(); i++) {
        auto pn      = normalize(positions[i] - center);
        positions[i] = center + pn * radius;
        normals[i]   = pn;
    }
}

// Make a bulged quad.
template <typename T>
inline void make_quad_bulged_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize, T height) {
    make_quad_shape(quads, positions, normals, texturecoords, {steps, steps},
        {size, size}, {uvsize, uvsize});
    if (height == 0) return;
    auto radius = (size * size / 4 + height * height) / (2 * height);
    auto center = vec<T, 3>{0, 0, -radius + height};
    for (auto i = 0; i < positions.size(); i++) {
        auto pn      = normalize(positions[i] - center);
        positions[i] = center + pn * radius;
        normals[i]   = pn;
    }
}

// Make a cylinder (side-only).
template <typename T>
inline void make_cylinder_side_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec2i& steps, const vec<T, 2>& size,
    const vec<T, 2>& uvsize) {
    make_quad_shape(
        quads, positions, normals, texturecoords, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv          = texturecoords[i];
        auto phi         = 2 * pif * uv.x;
        positions[i]     = {cos(phi) * size.x / 2, sin(phi) * size.x / 2,
            (uv.y - 0.5f) * size.y};
        normals[i]       = {cos(phi), sin(phi), 0};
        texturecoords[i] = uv * uvsize;
    }
}

// Make a cylinder.
template <typename T>
inline void make_uvcylinder_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 2>& size,
    const vec<T, 3>& uvsize) {
    quads.clear();
    positions.clear();
    normals.clear();
    texturecoords.clear();
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec<T, 3>>{};
    auto qnormals       = vector<vec<T, 3>>{};
    auto qtexturecoords = vector<vec<T, 2>>{};
    // side
    make_cylinder_side_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // top
    make_uvdisk_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.z}, size.x, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i].z = size.y / 2;
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // bottom
    make_uvdisk_shape(qquads, qpositions, qnormals, qtexturecoords,
        {steps.x, steps.z}, size.x, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i].z = -size.y / 2;
        qnormals[i]     = -qnormals[i];
    }
    for (auto i = 0; i < qquads.size(); i++) swap(qquads[i].x, qquads[i].z);
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
}

// Make a rounded cylinder.
template <typename T>
inline void make_uvcylinder_rounded_shape(vector<vec4i>& quads,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 2>& size,
    const vec<T, 3>& uvsize, T rounded) {
    make_uvcylinder_shape(
        quads, positions, normals, texturecoords, steps, size, uvsize);
    auto radius = rounded * max(size) / 2;
    auto c      = size / 2 - vec<T, 2>{radius, radius};
    for (auto i = 0; i < positions.size(); i++) {
        auto phi = atan2(positions[i].y, positions[i].x);
        auto r   = length(vec<T, 2>{positions[i].x, positions[i].y});
        auto z   = positions[i].z;
        auto pc  = vec<T, 2>{r, fabs(z)};
        auto ps  = (z < 0) ? -1.0f : 1.0f;
        if (pc.x >= c.x && pc.y >= c.y) {
            auto pn      = normalize(pc - c);
            positions[i] = {cos(phi) * (c.x + radius * pn.x),
                sin(phi) * (c.x + radius * pn.x), ps * (c.y + radius * pn.y)};
            normals[i]   = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
        } else {
            continue;
        }
    }
}

// Make a geodesic sphere.
template <typename T>
inline void make_geodesic_sphere_shape(vector<vec3i>& triangles,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals, int tesselation,
    T size) {
    // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
    const T     X                = 0.525731112119133606f;
    const T     Z                = 0.850650808352039932f;
    static auto sphere_pos       = vector<vec<T, 3>>{{-X, 0.0, Z}, {X, 0.0, Z},
        {-X, 0.0, -Z}, {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X},
        {0.0, -Z, -X}, {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
    static auto sphere_triangles = vector<vec3i>{{0, 1, 4}, {0, 4, 9},
        {9, 4, 5}, {4, 8, 5}, {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3},
        {5, 3, 2}, {2, 3, 7}, {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0},
        {0, 6, 1}, {6, 10, 1}, {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
    positions                    = sphere_pos;
    triangles                    = sphere_triangles;
    for (auto l = 0; l < max(0, tesselation - 2); l++) {
        subdivide_triangles(triangles, positions);
    }
    for (auto& p : positions) p = normalize(p) * size / 2;
    normals = positions;
}

// Make a facevarying cube.
template <typename T>
inline void make_box_fvshape(vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, const vec3i& steps, const vec<T, 3>& size,
    const vec<T, 3>& uvsize) {
    make_box_shape(quads_positions, positions, normals, texturecoords, steps,
        size, uvsize);
    quads_normals       = quads_positions;
    quads_texturecoords = quads_positions;
    auto positions_     = positions;
    weld_quads(quads_positions, positions,
        min(0.1f * size / vec<T, 3>{(T)steps.x, (T)steps.y, (T)steps.z}));
}

// Make a faceavrying spherecube.
template <typename T>
inline void make_sphere_fvshape(vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, int steps, T size, T uvsize) {
    make_box_fvshape(quads_positions, quads_normals, quads_texturecoords,
        positions, normals, texturecoords, {steps, steps, steps}, {1, 1, 1},
        {uvsize, uvsize, uvsize});
    quads_normals = quads_positions;
    normals       = positions;
    for (auto i = 0; i < positions.size(); i++) {
        auto p       = positions[i];
        positions[i] = normalize(p) * (size / 2);
        normals[i]   = normalize(p);
    }
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
template <typename T>
inline void make_suzanne_shape(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, T size) {
    static auto suzanne_pos = vector<vec<T, 3>>{{0.4375, 0.1640625, 0.765625},
        {-0.4375, 0.1640625, 0.765625}, {0.5, 0.09375, 0.6875},
        {-0.5, 0.09375, 0.6875}, {0.546875, 0.0546875, 0.578125},
        {-0.546875, 0.0546875, 0.578125}, {0.3515625, -0.0234375, 0.6171875},
        {-0.3515625, -0.0234375, 0.6171875}, {0.3515625, 0.03125, 0.71875},
        {-0.3515625, 0.03125, 0.71875}, {0.3515625, 0.1328125, 0.78125},
        {-0.3515625, 0.1328125, 0.78125}, {0.2734375, 0.1640625, 0.796875},
        {-0.2734375, 0.1640625, 0.796875}, {0.203125, 0.09375, 0.7421875},
        {-0.203125, 0.09375, 0.7421875}, {0.15625, 0.0546875, 0.6484375},
        {-0.15625, 0.0546875, 0.6484375}, {0.078125, 0.2421875, 0.65625},
        {-0.078125, 0.2421875, 0.65625}, {0.140625, 0.2421875, 0.7421875},
        {-0.140625, 0.2421875, 0.7421875}, {0.2421875, 0.2421875, 0.796875},
        {-0.2421875, 0.2421875, 0.796875}, {0.2734375, 0.328125, 0.796875},
        {-0.2734375, 0.328125, 0.796875}, {0.203125, 0.390625, 0.7421875},
        {-0.203125, 0.390625, 0.7421875}, {0.15625, 0.4375, 0.6484375},
        {-0.15625, 0.4375, 0.6484375}, {0.3515625, 0.515625, 0.6171875},
        {-0.3515625, 0.515625, 0.6171875}, {0.3515625, 0.453125, 0.71875},
        {-0.3515625, 0.453125, 0.71875}, {0.3515625, 0.359375, 0.78125},
        {-0.3515625, 0.359375, 0.78125}, {0.4375, 0.328125, 0.765625},
        {-0.4375, 0.328125, 0.765625}, {0.5, 0.390625, 0.6875},
        {-0.5, 0.390625, 0.6875}, {0.546875, 0.4375, 0.578125},
        {-0.546875, 0.4375, 0.578125}, {0.625, 0.2421875, 0.5625},
        {-0.625, 0.2421875, 0.5625}, {0.5625, 0.2421875, 0.671875},
        {-0.5625, 0.2421875, 0.671875}, {0.46875, 0.2421875, 0.7578125},
        {-0.46875, 0.2421875, 0.7578125}, {0.4765625, 0.2421875, 0.7734375},
        {-0.4765625, 0.2421875, 0.7734375}, {0.4453125, 0.3359375, 0.78125},
        {-0.4453125, 0.3359375, 0.78125}, {0.3515625, 0.375, 0.8046875},
        {-0.3515625, 0.375, 0.8046875}, {0.265625, 0.3359375, 0.8203125},
        {-0.265625, 0.3359375, 0.8203125}, {0.2265625, 0.2421875, 0.8203125},
        {-0.2265625, 0.2421875, 0.8203125}, {0.265625, 0.15625, 0.8203125},
        {-0.265625, 0.15625, 0.8203125}, {0.3515625, 0.2421875, 0.828125},
        {-0.3515625, 0.2421875, 0.828125}, {0.3515625, 0.1171875, 0.8046875},
        {-0.3515625, 0.1171875, 0.8046875}, {0.4453125, 0.15625, 0.78125},
        {-0.4453125, 0.15625, 0.78125}, {0.0, 0.4296875, 0.7421875},
        {0.0, 0.3515625, 0.8203125}, {0.0, -0.6796875, 0.734375},
        {0.0, -0.3203125, 0.78125}, {0.0, -0.1875, 0.796875},
        {0.0, -0.7734375, 0.71875}, {0.0, 0.40625, 0.6015625},
        {0.0, 0.5703125, 0.5703125}, {0.0, 0.8984375, -0.546875},
        {0.0, 0.5625, -0.8515625}, {0.0, 0.0703125, -0.828125},
        {0.0, -0.3828125, -0.3515625}, {0.203125, -0.1875, 0.5625},
        {-0.203125, -0.1875, 0.5625}, {0.3125, -0.4375, 0.5703125},
        {-0.3125, -0.4375, 0.5703125}, {0.3515625, -0.6953125, 0.5703125},
        {-0.3515625, -0.6953125, 0.5703125}, {0.3671875, -0.890625, 0.53125},
        {-0.3671875, -0.890625, 0.53125}, {0.328125, -0.9453125, 0.5234375},
        {-0.328125, -0.9453125, 0.5234375}, {0.1796875, -0.96875, 0.5546875},
        {-0.1796875, -0.96875, 0.5546875}, {0.0, -0.984375, 0.578125},
        {0.4375, -0.140625, 0.53125}, {-0.4375, -0.140625, 0.53125},
        {0.6328125, -0.0390625, 0.5390625}, {-0.6328125, -0.0390625, 0.5390625},
        {0.828125, 0.1484375, 0.4453125}, {-0.828125, 0.1484375, 0.4453125},
        {0.859375, 0.4296875, 0.59375}, {-0.859375, 0.4296875, 0.59375},
        {0.7109375, 0.484375, 0.625}, {-0.7109375, 0.484375, 0.625},
        {0.4921875, 0.6015625, 0.6875}, {-0.4921875, 0.6015625, 0.6875},
        {0.3203125, 0.7578125, 0.734375}, {-0.3203125, 0.7578125, 0.734375},
        {0.15625, 0.71875, 0.7578125}, {-0.15625, 0.71875, 0.7578125},
        {0.0625, 0.4921875, 0.75}, {-0.0625, 0.4921875, 0.75},
        {0.1640625, 0.4140625, 0.7734375}, {-0.1640625, 0.4140625, 0.7734375},
        {0.125, 0.3046875, 0.765625}, {-0.125, 0.3046875, 0.765625},
        {0.203125, 0.09375, 0.7421875}, {-0.203125, 0.09375, 0.7421875},
        {0.375, 0.015625, 0.703125}, {-0.375, 0.015625, 0.703125},
        {0.4921875, 0.0625, 0.671875}, {-0.4921875, 0.0625, 0.671875},
        {0.625, 0.1875, 0.6484375}, {-0.625, 0.1875, 0.6484375},
        {0.640625, 0.296875, 0.6484375}, {-0.640625, 0.296875, 0.6484375},
        {0.6015625, 0.375, 0.6640625}, {-0.6015625, 0.375, 0.6640625},
        {0.4296875, 0.4375, 0.71875}, {-0.4296875, 0.4375, 0.71875},
        {0.25, 0.46875, 0.7578125}, {-0.25, 0.46875, 0.7578125},
        {0.0, -0.765625, 0.734375}, {0.109375, -0.71875, 0.734375},
        {-0.109375, -0.71875, 0.734375}, {0.1171875, -0.8359375, 0.7109375},
        {-0.1171875, -0.8359375, 0.7109375}, {0.0625, -0.8828125, 0.6953125},
        {-0.0625, -0.8828125, 0.6953125}, {0.0, -0.890625, 0.6875},
        {0.0, -0.1953125, 0.75}, {0.0, -0.140625, 0.7421875},
        {0.1015625, -0.1484375, 0.7421875}, {-0.1015625, -0.1484375, 0.7421875},
        {0.125, -0.2265625, 0.75}, {-0.125, -0.2265625, 0.75},
        {0.0859375, -0.2890625, 0.7421875}, {-0.0859375, -0.2890625, 0.7421875},
        {0.3984375, -0.046875, 0.671875}, {-0.3984375, -0.046875, 0.671875},
        {0.6171875, 0.0546875, 0.625}, {-0.6171875, 0.0546875, 0.625},
        {0.7265625, 0.203125, 0.6015625}, {-0.7265625, 0.203125, 0.6015625},
        {0.7421875, 0.375, 0.65625}, {-0.7421875, 0.375, 0.65625},
        {0.6875, 0.4140625, 0.7265625}, {-0.6875, 0.4140625, 0.7265625},
        {0.4375, 0.546875, 0.796875}, {-0.4375, 0.546875, 0.796875},
        {0.3125, 0.640625, 0.8359375}, {-0.3125, 0.640625, 0.8359375},
        {0.203125, 0.6171875, 0.8515625}, {-0.203125, 0.6171875, 0.8515625},
        {0.1015625, 0.4296875, 0.84375}, {-0.1015625, 0.4296875, 0.84375},
        {0.125, -0.1015625, 0.8125}, {-0.125, -0.1015625, 0.8125},
        {0.2109375, -0.4453125, 0.7109375}, {-0.2109375, -0.4453125, 0.7109375},
        {0.25, -0.703125, 0.6875}, {-0.25, -0.703125, 0.6875},
        {0.265625, -0.8203125, 0.6640625}, {-0.265625, -0.8203125, 0.6640625},
        {0.234375, -0.9140625, 0.6328125}, {-0.234375, -0.9140625, 0.6328125},
        {0.1640625, -0.9296875, 0.6328125}, {-0.1640625, -0.9296875, 0.6328125},
        {0.0, -0.9453125, 0.640625}, {0.0, 0.046875, 0.7265625},
        {0.0, 0.2109375, 0.765625}, {0.328125, 0.4765625, 0.7421875},
        {-0.328125, 0.4765625, 0.7421875}, {0.1640625, 0.140625, 0.75},
        {-0.1640625, 0.140625, 0.75}, {0.1328125, 0.2109375, 0.7578125},
        {-0.1328125, 0.2109375, 0.7578125}, {0.1171875, -0.6875, 0.734375},
        {-0.1171875, -0.6875, 0.734375}, {0.078125, -0.4453125, 0.75},
        {-0.078125, -0.4453125, 0.75}, {0.0, -0.4453125, 0.75},
        {0.0, -0.328125, 0.7421875}, {0.09375, -0.2734375, 0.78125},
        {-0.09375, -0.2734375, 0.78125}, {0.1328125, -0.2265625, 0.796875},
        {-0.1328125, -0.2265625, 0.796875}, {0.109375, -0.1328125, 0.78125},
        {-0.109375, -0.1328125, 0.78125}, {0.0390625, -0.125, 0.78125},
        {-0.0390625, -0.125, 0.78125}, {0.0, -0.203125, 0.828125},
        {0.046875, -0.1484375, 0.8125}, {-0.046875, -0.1484375, 0.8125},
        {0.09375, -0.15625, 0.8125}, {-0.09375, -0.15625, 0.8125},
        {0.109375, -0.2265625, 0.828125}, {-0.109375, -0.2265625, 0.828125},
        {0.078125, -0.25, 0.8046875}, {-0.078125, -0.25, 0.8046875},
        {0.0, -0.2890625, 0.8046875}, {0.2578125, -0.3125, 0.5546875},
        {-0.2578125, -0.3125, 0.5546875}, {0.1640625, -0.2421875, 0.7109375},
        {-0.1640625, -0.2421875, 0.7109375}, {0.1796875, -0.3125, 0.7109375},
        {-0.1796875, -0.3125, 0.7109375}, {0.234375, -0.25, 0.5546875},
        {-0.234375, -0.25, 0.5546875}, {0.0, -0.875, 0.6875},
        {0.046875, -0.8671875, 0.6875}, {-0.046875, -0.8671875, 0.6875},
        {0.09375, -0.8203125, 0.7109375}, {-0.09375, -0.8203125, 0.7109375},
        {0.09375, -0.7421875, 0.7265625}, {-0.09375, -0.7421875, 0.7265625},
        {0.0, -0.78125, 0.65625}, {0.09375, -0.75, 0.6640625},
        {-0.09375, -0.75, 0.6640625}, {0.09375, -0.8125, 0.640625},
        {-0.09375, -0.8125, 0.640625}, {0.046875, -0.8515625, 0.6328125},
        {-0.046875, -0.8515625, 0.6328125}, {0.0, -0.859375, 0.6328125},
        {0.171875, 0.21875, 0.78125}, {-0.171875, 0.21875, 0.78125},
        {0.1875, 0.15625, 0.7734375}, {-0.1875, 0.15625, 0.7734375},
        {0.3359375, 0.4296875, 0.7578125}, {-0.3359375, 0.4296875, 0.7578125},
        {0.2734375, 0.421875, 0.7734375}, {-0.2734375, 0.421875, 0.7734375},
        {0.421875, 0.3984375, 0.7734375}, {-0.421875, 0.3984375, 0.7734375},
        {0.5625, 0.3515625, 0.6953125}, {-0.5625, 0.3515625, 0.6953125},
        {0.5859375, 0.2890625, 0.6875}, {-0.5859375, 0.2890625, 0.6875},
        {0.578125, 0.1953125, 0.6796875}, {-0.578125, 0.1953125, 0.6796875},
        {0.4765625, 0.1015625, 0.71875}, {-0.4765625, 0.1015625, 0.71875},
        {0.375, 0.0625, 0.7421875}, {-0.375, 0.0625, 0.7421875},
        {0.2265625, 0.109375, 0.78125}, {-0.2265625, 0.109375, 0.78125},
        {0.1796875, 0.296875, 0.78125}, {-0.1796875, 0.296875, 0.78125},
        {0.2109375, 0.375, 0.78125}, {-0.2109375, 0.375, 0.78125},
        {0.234375, 0.359375, 0.7578125}, {-0.234375, 0.359375, 0.7578125},
        {0.1953125, 0.296875, 0.7578125}, {-0.1953125, 0.296875, 0.7578125},
        {0.2421875, 0.125, 0.7578125}, {-0.2421875, 0.125, 0.7578125},
        {0.375, 0.0859375, 0.7265625}, {-0.375, 0.0859375, 0.7265625},
        {0.4609375, 0.1171875, 0.703125}, {-0.4609375, 0.1171875, 0.703125},
        {0.546875, 0.2109375, 0.671875}, {-0.546875, 0.2109375, 0.671875},
        {0.5546875, 0.28125, 0.671875}, {-0.5546875, 0.28125, 0.671875},
        {0.53125, 0.3359375, 0.6796875}, {-0.53125, 0.3359375, 0.6796875},
        {0.4140625, 0.390625, 0.75}, {-0.4140625, 0.390625, 0.75},
        {0.28125, 0.3984375, 0.765625}, {-0.28125, 0.3984375, 0.765625},
        {0.3359375, 0.40625, 0.75}, {-0.3359375, 0.40625, 0.75},
        {0.203125, 0.171875, 0.75}, {-0.203125, 0.171875, 0.75},
        {0.1953125, 0.2265625, 0.75}, {-0.1953125, 0.2265625, 0.75},
        {0.109375, 0.4609375, 0.609375}, {-0.109375, 0.4609375, 0.609375},
        {0.1953125, 0.6640625, 0.6171875}, {-0.1953125, 0.6640625, 0.6171875},
        {0.3359375, 0.6875, 0.59375}, {-0.3359375, 0.6875, 0.59375},
        {0.484375, 0.5546875, 0.5546875}, {-0.484375, 0.5546875, 0.5546875},
        {0.6796875, 0.453125, 0.4921875}, {-0.6796875, 0.453125, 0.4921875},
        {0.796875, 0.40625, 0.4609375}, {-0.796875, 0.40625, 0.4609375},
        {0.7734375, 0.1640625, 0.375}, {-0.7734375, 0.1640625, 0.375},
        {0.6015625, 0.0, 0.4140625}, {-0.6015625, 0.0, 0.4140625},
        {0.4375, -0.09375, 0.46875}, {-0.4375, -0.09375, 0.46875},
        {0.0, 0.8984375, 0.2890625}, {0.0, 0.984375, -0.078125},
        {0.0, -0.1953125, -0.671875}, {0.0, -0.4609375, 0.1875},
        {0.0, -0.9765625, 0.4609375}, {0.0, -0.8046875, 0.34375},
        {0.0, -0.5703125, 0.3203125}, {0.0, -0.484375, 0.28125},
        {0.8515625, 0.234375, 0.0546875}, {-0.8515625, 0.234375, 0.0546875},
        {0.859375, 0.3203125, -0.046875}, {-0.859375, 0.3203125, -0.046875},
        {0.7734375, 0.265625, -0.4375}, {-0.7734375, 0.265625, -0.4375},
        {0.4609375, 0.4375, -0.703125}, {-0.4609375, 0.4375, -0.703125},
        {0.734375, -0.046875, 0.0703125}, {-0.734375, -0.046875, 0.0703125},
        {0.59375, -0.125, -0.1640625}, {-0.59375, -0.125, -0.1640625},
        {0.640625, -0.0078125, -0.4296875}, {-0.640625, -0.0078125, -0.4296875},
        {0.3359375, 0.0546875, -0.6640625}, {-0.3359375, 0.0546875, -0.6640625},
        {0.234375, -0.3515625, 0.40625}, {-0.234375, -0.3515625, 0.40625},
        {0.1796875, -0.4140625, 0.2578125}, {-0.1796875, -0.4140625, 0.2578125},
        {0.2890625, -0.7109375, 0.3828125}, {-0.2890625, -0.7109375, 0.3828125},
        {0.25, -0.5, 0.390625}, {-0.25, -0.5, 0.390625},
        {0.328125, -0.9140625, 0.3984375}, {-0.328125, -0.9140625, 0.3984375},
        {0.140625, -0.7578125, 0.3671875}, {-0.140625, -0.7578125, 0.3671875},
        {0.125, -0.5390625, 0.359375}, {-0.125, -0.5390625, 0.359375},
        {0.1640625, -0.9453125, 0.4375}, {-0.1640625, -0.9453125, 0.4375},
        {0.21875, -0.28125, 0.4296875}, {-0.21875, -0.28125, 0.4296875},
        {0.2109375, -0.2265625, 0.46875}, {-0.2109375, -0.2265625, 0.46875},
        {0.203125, -0.171875, 0.5}, {-0.203125, -0.171875, 0.5},
        {0.2109375, -0.390625, 0.1640625}, {-0.2109375, -0.390625, 0.1640625},
        {0.296875, -0.3125, -0.265625}, {-0.296875, -0.3125, -0.265625},
        {0.34375, -0.1484375, -0.5390625}, {-0.34375, -0.1484375, -0.5390625},
        {0.453125, 0.8671875, -0.3828125}, {-0.453125, 0.8671875, -0.3828125},
        {0.453125, 0.9296875, -0.0703125}, {-0.453125, 0.9296875, -0.0703125},
        {0.453125, 0.8515625, 0.234375}, {-0.453125, 0.8515625, 0.234375},
        {0.4609375, 0.5234375, 0.4296875}, {-0.4609375, 0.5234375, 0.4296875},
        {0.7265625, 0.40625, 0.3359375}, {-0.7265625, 0.40625, 0.3359375},
        {0.6328125, 0.453125, 0.28125}, {-0.6328125, 0.453125, 0.28125},
        {0.640625, 0.703125, 0.0546875}, {-0.640625, 0.703125, 0.0546875},
        {0.796875, 0.5625, 0.125}, {-0.796875, 0.5625, 0.125},
        {0.796875, 0.6171875, -0.1171875}, {-0.796875, 0.6171875, -0.1171875},
        {0.640625, 0.75, -0.1953125}, {-0.640625, 0.75, -0.1953125},
        {0.640625, 0.6796875, -0.4453125}, {-0.640625, 0.6796875, -0.4453125},
        {0.796875, 0.5390625, -0.359375}, {-0.796875, 0.5390625, -0.359375},
        {0.6171875, 0.328125, -0.5859375}, {-0.6171875, 0.328125, -0.5859375},
        {0.484375, 0.0234375, -0.546875}, {-0.484375, 0.0234375, -0.546875},
        {0.8203125, 0.328125, -0.203125}, {-0.8203125, 0.328125, -0.203125},
        {0.40625, -0.171875, 0.1484375}, {-0.40625, -0.171875, 0.1484375},
        {0.4296875, -0.1953125, -0.2109375},
        {-0.4296875, -0.1953125, -0.2109375}, {0.890625, 0.40625, -0.234375},
        {-0.890625, 0.40625, -0.234375}, {0.7734375, -0.140625, -0.125},
        {-0.7734375, -0.140625, -0.125}, {1.0390625, -0.1015625, -0.328125},
        {-1.0390625, -0.1015625, -0.328125}, {1.28125, 0.0546875, -0.4296875},
        {-1.28125, 0.0546875, -0.4296875}, {1.3515625, 0.3203125, -0.421875},
        {-1.3515625, 0.3203125, -0.421875}, {1.234375, 0.5078125, -0.421875},
        {-1.234375, 0.5078125, -0.421875}, {1.0234375, 0.4765625, -0.3125},
        {-1.0234375, 0.4765625, -0.3125}, {1.015625, 0.4140625, -0.2890625},
        {-1.015625, 0.4140625, -0.2890625}, {1.1875, 0.4375, -0.390625},
        {-1.1875, 0.4375, -0.390625}, {1.265625, 0.2890625, -0.40625},
        {-1.265625, 0.2890625, -0.40625}, {1.2109375, 0.078125, -0.40625},
        {-1.2109375, 0.078125, -0.40625}, {1.03125, -0.0390625, -0.3046875},
        {-1.03125, -0.0390625, -0.3046875}, {0.828125, -0.0703125, -0.1328125},
        {-0.828125, -0.0703125, -0.1328125}, {0.921875, 0.359375, -0.21875},
        {-0.921875, 0.359375, -0.21875}, {0.9453125, 0.3046875, -0.2890625},
        {-0.9453125, 0.3046875, -0.2890625},
        {0.8828125, -0.0234375, -0.2109375},
        {-0.8828125, -0.0234375, -0.2109375}, {1.0390625, 0.0, -0.3671875},
        {-1.0390625, 0.0, -0.3671875}, {1.1875, 0.09375, -0.4453125},
        {-1.1875, 0.09375, -0.4453125}, {1.234375, 0.25, -0.4453125},
        {-1.234375, 0.25, -0.4453125}, {1.171875, 0.359375, -0.4375},
        {-1.171875, 0.359375, -0.4375}, {1.0234375, 0.34375, -0.359375},
        {-1.0234375, 0.34375, -0.359375}, {0.84375, 0.2890625, -0.2109375},
        {-0.84375, 0.2890625, -0.2109375}, {0.8359375, 0.171875, -0.2734375},
        {-0.8359375, 0.171875, -0.2734375}, {0.7578125, 0.09375, -0.2734375},
        {-0.7578125, 0.09375, -0.2734375}, {0.8203125, 0.0859375, -0.2734375},
        {-0.8203125, 0.0859375, -0.2734375}, {0.84375, 0.015625, -0.2734375},
        {-0.84375, 0.015625, -0.2734375}, {0.8125, -0.015625, -0.2734375},
        {-0.8125, -0.015625, -0.2734375}, {0.7265625, 0.0, -0.0703125},
        {-0.7265625, 0.0, -0.0703125}, {0.71875, -0.0234375, -0.171875},
        {-0.71875, -0.0234375, -0.171875}, {0.71875, 0.0390625, -0.1875},
        {-0.71875, 0.0390625, -0.1875}, {0.796875, 0.203125, -0.2109375},
        {-0.796875, 0.203125, -0.2109375}, {0.890625, 0.2421875, -0.265625},
        {-0.890625, 0.2421875, -0.265625}, {0.890625, 0.234375, -0.3203125},
        {-0.890625, 0.234375, -0.3203125}, {0.8125, -0.015625, -0.3203125},
        {-0.8125, -0.015625, -0.3203125}, {0.8515625, 0.015625, -0.3203125},
        {-0.8515625, 0.015625, -0.3203125}, {0.828125, 0.078125, -0.3203125},
        {-0.828125, 0.078125, -0.3203125}, {0.765625, 0.09375, -0.3203125},
        {-0.765625, 0.09375, -0.3203125}, {0.84375, 0.171875, -0.3203125},
        {-0.84375, 0.171875, -0.3203125}, {1.0390625, 0.328125, -0.4140625},
        {-1.0390625, 0.328125, -0.4140625}, {1.1875, 0.34375, -0.484375},
        {-1.1875, 0.34375, -0.484375}, {1.2578125, 0.2421875, -0.4921875},
        {-1.2578125, 0.2421875, -0.4921875}, {1.2109375, 0.0859375, -0.484375},
        {-1.2109375, 0.0859375, -0.484375}, {1.046875, 0.0, -0.421875},
        {-1.046875, 0.0, -0.421875}, {0.8828125, -0.015625, -0.265625},
        {-0.8828125, -0.015625, -0.265625}, {0.953125, 0.2890625, -0.34375},
        {-0.953125, 0.2890625, -0.34375}, {0.890625, 0.109375, -0.328125},
        {-0.890625, 0.109375, -0.328125}, {0.9375, 0.0625, -0.3359375},
        {-0.9375, 0.0625, -0.3359375}, {1.0, 0.125, -0.3671875},
        {-1.0, 0.125, -0.3671875}, {0.9609375, 0.171875, -0.3515625},
        {-0.9609375, 0.171875, -0.3515625}, {1.015625, 0.234375, -0.375},
        {-1.015625, 0.234375, -0.375}, {1.0546875, 0.1875, -0.3828125},
        {-1.0546875, 0.1875, -0.3828125}, {1.109375, 0.2109375, -0.390625},
        {-1.109375, 0.2109375, -0.390625}, {1.0859375, 0.2734375, -0.390625},
        {-1.0859375, 0.2734375, -0.390625}, {1.0234375, 0.4375, -0.484375},
        {-1.0234375, 0.4375, -0.484375}, {1.25, 0.46875, -0.546875},
        {-1.25, 0.46875, -0.546875}, {1.3671875, 0.296875, -0.5},
        {-1.3671875, 0.296875, -0.5}, {1.3125, 0.0546875, -0.53125},
        {-1.3125, 0.0546875, -0.53125}, {1.0390625, -0.0859375, -0.4921875},
        {-1.0390625, -0.0859375, -0.4921875}, {0.7890625, -0.125, -0.328125},
        {-0.7890625, -0.125, -0.328125}, {0.859375, 0.3828125, -0.3828125},
        {-0.859375, 0.3828125, -0.3828125}};
    static auto suzanne_triangles = vector<vec3i>{{60, 64, 48}, {49, 65, 61},
        {62, 64, 60}, {61, 65, 63}, {60, 58, 62}, {63, 59, 61}, {60, 56, 58},
        {59, 57, 61}, {60, 54, 56}, {57, 55, 61}, {60, 52, 54}, {55, 53, 61},
        {60, 50, 52}, {53, 51, 61}, {60, 48, 50}, {51, 49, 61}, {224, 228, 226},
        {227, 229, 225}, {72, 283, 73}, {73, 284, 72}, {341, 347, 383},
        {384, 348, 342}, {299, 345, 343}, {344, 346, 300}, {323, 379, 351},
        {352, 380, 324}, {441, 443, 445}, {446, 444, 442}, {463, 491, 465},
        {466, 492, 464}, {495, 497, 499}, {500, 498, 496}};
    static auto suzanne_quads = vector<vec4i>{{46, 0, 2, 44}, {3, 1, 47, 45},
        {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4}, {7, 9, 3, 5},
        {0, 10, 8, 2}, {9, 11, 1, 3}, {10, 12, 14, 8}, {15, 13, 11, 9},
        {8, 14, 16, 6}, {17, 15, 9, 7}, {14, 20, 18, 16}, {19, 21, 15, 17},
        {12, 22, 20, 14}, {21, 23, 13, 15}, {22, 24, 26, 20}, {27, 25, 23, 21},
        {20, 26, 28, 18}, {29, 27, 21, 19}, {26, 32, 30, 28}, {31, 33, 27, 29},
        {24, 34, 32, 26}, {33, 35, 25, 27}, {34, 36, 38, 32}, {39, 37, 35, 33},
        {32, 38, 40, 30}, {41, 39, 33, 31}, {38, 44, 42, 40}, {43, 45, 39, 41},
        {36, 46, 44, 38}, {45, 47, 37, 39}, {46, 36, 50, 48}, {51, 37, 47, 49},
        {36, 34, 52, 50}, {53, 35, 37, 51}, {34, 24, 54, 52}, {55, 25, 35, 53},
        {24, 22, 56, 54}, {57, 23, 25, 55}, {22, 12, 58, 56}, {59, 13, 23, 57},
        {12, 10, 62, 58}, {63, 11, 13, 59}, {10, 0, 64, 62}, {65, 1, 11, 63},
        {0, 46, 48, 64}, {49, 47, 1, 65}, {88, 173, 175, 90},
        {175, 174, 89, 90}, {86, 171, 173, 88}, {174, 172, 87, 89},
        {84, 169, 171, 86}, {172, 170, 85, 87}, {82, 167, 169, 84},
        {170, 168, 83, 85}, {80, 165, 167, 82}, {168, 166, 81, 83},
        {78, 91, 145, 163}, {146, 92, 79, 164}, {91, 93, 147, 145},
        {148, 94, 92, 146}, {93, 95, 149, 147}, {150, 96, 94, 148},
        {95, 97, 151, 149}, {152, 98, 96, 150}, {97, 99, 153, 151},
        {154, 100, 98, 152}, {99, 101, 155, 153}, {156, 102, 100, 154},
        {101, 103, 157, 155}, {158, 104, 102, 156}, {103, 105, 159, 157},
        {160, 106, 104, 158}, {105, 107, 161, 159}, {162, 108, 106, 160},
        {107, 66, 67, 161}, {67, 66, 108, 162}, {109, 127, 159, 161},
        {160, 128, 110, 162}, {127, 178, 157, 159}, {158, 179, 128, 160},
        {125, 155, 157, 178}, {158, 156, 126, 179}, {123, 153, 155, 125},
        {156, 154, 124, 126}, {121, 151, 153, 123}, {154, 152, 122, 124},
        {119, 149, 151, 121}, {152, 150, 120, 122}, {117, 147, 149, 119},
        {150, 148, 118, 120}, {115, 145, 147, 117}, {148, 146, 116, 118},
        {113, 163, 145, 115}, {146, 164, 114, 116}, {113, 180, 176, 163},
        {176, 181, 114, 164}, {109, 161, 67, 111}, {67, 162, 110, 112},
        {111, 67, 177, 182}, {177, 67, 112, 183}, {176, 180, 182, 177},
        {183, 181, 176, 177}, {134, 136, 175, 173}, {175, 136, 135, 174},
        {132, 134, 173, 171}, {174, 135, 133, 172}, {130, 132, 171, 169},
        {172, 133, 131, 170}, {165, 186, 184, 167}, {185, 187, 166, 168},
        {130, 169, 167, 184}, {168, 170, 131, 185}, {143, 189, 188, 186},
        {188, 189, 144, 187}, {184, 186, 188, 68}, {188, 187, 185, 68},
        {129, 130, 184, 68}, {185, 131, 129, 68}, {141, 192, 190, 143},
        {191, 193, 142, 144}, {139, 194, 192, 141}, {193, 195, 140, 142},
        {138, 196, 194, 139}, {195, 197, 138, 140}, {137, 70, 196, 138},
        {197, 70, 137, 138}, {189, 143, 190, 69}, {191, 144, 189, 69},
        {69, 190, 205, 207}, {206, 191, 69, 207}, {70, 198, 199, 196},
        {200, 198, 70, 197}, {196, 199, 201, 194}, {202, 200, 197, 195},
        {194, 201, 203, 192}, {204, 202, 195, 193}, {192, 203, 205, 190},
        {206, 204, 193, 191}, {198, 203, 201, 199}, {202, 204, 198, 200},
        {198, 207, 205, 203}, {206, 207, 198, 204}, {138, 139, 163, 176},
        {164, 140, 138, 176}, {139, 141, 210, 163}, {211, 142, 140, 164},
        {141, 143, 212, 210}, {213, 144, 142, 211}, {143, 186, 165, 212},
        {166, 187, 144, 213}, {80, 208, 212, 165}, {213, 209, 81, 166},
        {208, 214, 210, 212}, {211, 215, 209, 213}, {78, 163, 210, 214},
        {211, 164, 79, 215}, {130, 129, 71, 221}, {71, 129, 131, 222},
        {132, 130, 221, 219}, {222, 131, 133, 220}, {134, 132, 219, 217},
        {220, 133, 135, 218}, {136, 134, 217, 216}, {218, 135, 136, 216},
        {216, 217, 228, 230}, {229, 218, 216, 230}, {217, 219, 226, 228},
        {227, 220, 218, 229}, {219, 221, 224, 226}, {225, 222, 220, 227},
        {221, 71, 223, 224}, {223, 71, 222, 225}, {223, 230, 228, 224},
        {229, 230, 223, 225}, {182, 180, 233, 231}, {234, 181, 183, 232},
        {111, 182, 231, 253}, {232, 183, 112, 254}, {109, 111, 253, 255},
        {254, 112, 110, 256}, {180, 113, 251, 233}, {252, 114, 181, 234},
        {113, 115, 249, 251}, {250, 116, 114, 252}, {115, 117, 247, 249},
        {248, 118, 116, 250}, {117, 119, 245, 247}, {246, 120, 118, 248},
        {119, 121, 243, 245}, {244, 122, 120, 246}, {121, 123, 241, 243},
        {242, 124, 122, 244}, {123, 125, 239, 241}, {240, 126, 124, 242},
        {125, 178, 235, 239}, {236, 179, 126, 240}, {178, 127, 237, 235},
        {238, 128, 179, 236}, {127, 109, 255, 237}, {256, 110, 128, 238},
        {237, 255, 257, 275}, {258, 256, 238, 276}, {235, 237, 275, 277},
        {276, 238, 236, 278}, {239, 235, 277, 273}, {278, 236, 240, 274},
        {241, 239, 273, 271}, {274, 240, 242, 272}, {243, 241, 271, 269},
        {272, 242, 244, 270}, {245, 243, 269, 267}, {270, 244, 246, 268},
        {247, 245, 267, 265}, {268, 246, 248, 266}, {249, 247, 265, 263},
        {266, 248, 250, 264}, {251, 249, 263, 261}, {264, 250, 252, 262},
        {233, 251, 261, 279}, {262, 252, 234, 280}, {255, 253, 259, 257},
        {260, 254, 256, 258}, {253, 231, 281, 259}, {282, 232, 254, 260},
        {231, 233, 279, 281}, {280, 234, 232, 282}, {66, 107, 283, 72},
        {284, 108, 66, 72}, {107, 105, 285, 283}, {286, 106, 108, 284},
        {105, 103, 287, 285}, {288, 104, 106, 286}, {103, 101, 289, 287},
        {290, 102, 104, 288}, {101, 99, 291, 289}, {292, 100, 102, 290},
        {99, 97, 293, 291}, {294, 98, 100, 292}, {97, 95, 295, 293},
        {296, 96, 98, 294}, {95, 93, 297, 295}, {298, 94, 96, 296},
        {93, 91, 299, 297}, {300, 92, 94, 298}, {307, 308, 327, 337},
        {328, 308, 307, 338}, {306, 307, 337, 335}, {338, 307, 306, 336},
        {305, 306, 335, 339}, {336, 306, 305, 340}, {88, 90, 305, 339},
        {305, 90, 89, 340}, {86, 88, 339, 333}, {340, 89, 87, 334},
        {84, 86, 333, 329}, {334, 87, 85, 330}, {82, 84, 329, 331},
        {330, 85, 83, 332}, {329, 335, 337, 331}, {338, 336, 330, 332},
        {329, 333, 339, 335}, {340, 334, 330, 336}, {325, 331, 337, 327},
        {338, 332, 326, 328}, {80, 82, 331, 325}, {332, 83, 81, 326},
        {208, 341, 343, 214}, {344, 342, 209, 215}, {80, 325, 341, 208},
        {342, 326, 81, 209}, {78, 214, 343, 345}, {344, 215, 79, 346},
        {78, 345, 299, 91}, {300, 346, 79, 92}, {76, 323, 351, 303},
        {352, 324, 76, 303}, {303, 351, 349, 77}, {350, 352, 303, 77},
        {77, 349, 347, 304}, {348, 350, 77, 304}, {304, 347, 327, 308},
        {328, 348, 304, 308}, {325, 327, 347, 341}, {348, 328, 326, 342},
        {295, 297, 317, 309}, {318, 298, 296, 310}, {75, 315, 323, 76},
        {324, 316, 75, 76}, {301, 357, 355, 302}, {356, 358, 301, 302},
        {302, 355, 353, 74}, {354, 356, 302, 74}, {74, 353, 315, 75},
        {316, 354, 74, 75}, {291, 293, 361, 363}, {362, 294, 292, 364},
        {363, 361, 367, 365}, {368, 362, 364, 366}, {365, 367, 369, 371},
        {370, 368, 366, 372}, {371, 369, 375, 373}, {376, 370, 372, 374},
        {313, 377, 373, 375}, {374, 378, 314, 376}, {315, 353, 373, 377},
        {374, 354, 316, 378}, {353, 355, 371, 373}, {372, 356, 354, 374},
        {355, 357, 365, 371}, {366, 358, 356, 372}, {357, 359, 363, 365},
        {364, 360, 358, 366}, {289, 291, 363, 359}, {364, 292, 290, 360},
        {73, 359, 357, 301}, {358, 360, 73, 301}, {283, 285, 287, 289},
        {288, 286, 284, 290}, {283, 289, 359, 73}, {360, 290, 284, 73},
        {293, 295, 309, 361}, {310, 296, 294, 362}, {309, 311, 367, 361},
        {368, 312, 310, 362}, {311, 381, 369, 367}, {370, 382, 312, 368},
        {313, 375, 369, 381}, {370, 376, 314, 382}, {347, 349, 385, 383},
        {386, 350, 348, 384}, {317, 383, 385, 319}, {386, 384, 318, 320},
        {297, 299, 383, 317}, {384, 300, 298, 318}, {299, 343, 341, 383},
        {342, 344, 300, 384}, {313, 321, 379, 377}, {380, 322, 314, 378},
        {315, 377, 379, 323}, {380, 378, 316, 324}, {319, 385, 379, 321},
        {380, 386, 320, 322}, {349, 351, 379, 385}, {380, 352, 350, 386},
        {399, 387, 413, 401}, {414, 388, 400, 402}, {399, 401, 403, 397},
        {404, 402, 400, 398}, {397, 403, 405, 395}, {406, 404, 398, 396},
        {395, 405, 407, 393}, {408, 406, 396, 394}, {393, 407, 409, 391},
        {410, 408, 394, 392}, {391, 409, 411, 389}, {412, 410, 392, 390},
        {409, 419, 417, 411}, {418, 420, 410, 412}, {407, 421, 419, 409},
        {420, 422, 408, 410}, {405, 423, 421, 407}, {422, 424, 406, 408},
        {403, 425, 423, 405}, {424, 426, 404, 406}, {401, 427, 425, 403},
        {426, 428, 402, 404}, {401, 413, 415, 427}, {416, 414, 402, 428},
        {317, 319, 443, 441}, {444, 320, 318, 442}, {319, 389, 411, 443},
        {412, 390, 320, 444}, {309, 317, 441, 311}, {442, 318, 310, 312},
        {381, 429, 413, 387}, {414, 430, 382, 388}, {411, 417, 439, 443},
        {440, 418, 412, 444}, {437, 445, 443, 439}, {444, 446, 438, 440},
        {433, 445, 437, 435}, {438, 446, 434, 436}, {431, 447, 445, 433},
        {446, 448, 432, 434}, {429, 447, 431, 449}, {432, 448, 430, 450},
        {413, 429, 449, 415}, {450, 430, 414, 416}, {311, 447, 429, 381},
        {430, 448, 312, 382}, {311, 441, 445, 447}, {446, 442, 312, 448},
        {415, 449, 451, 475}, {452, 450, 416, 476}, {449, 431, 461, 451},
        {462, 432, 450, 452}, {431, 433, 459, 461}, {460, 434, 432, 462},
        {433, 435, 457, 459}, {458, 436, 434, 460}, {435, 437, 455, 457},
        {456, 438, 436, 458}, {437, 439, 453, 455}, {454, 440, 438, 456},
        {439, 417, 473, 453}, {474, 418, 440, 454}, {427, 415, 475, 463},
        {476, 416, 428, 464}, {425, 427, 463, 465}, {464, 428, 426, 466},
        {423, 425, 465, 467}, {466, 426, 424, 468}, {421, 423, 467, 469},
        {468, 424, 422, 470}, {419, 421, 469, 471}, {470, 422, 420, 472},
        {417, 419, 471, 473}, {472, 420, 418, 474}, {457, 455, 479, 477},
        {480, 456, 458, 478}, {477, 479, 481, 483}, {482, 480, 478, 484},
        {483, 481, 487, 485}, {488, 482, 484, 486}, {485, 487, 489, 491},
        {490, 488, 486, 492}, {463, 475, 485, 491}, {486, 476, 464, 492},
        {451, 483, 485, 475}, {486, 484, 452, 476}, {451, 461, 477, 483},
        {478, 462, 452, 484}, {457, 477, 461, 459}, {462, 478, 458, 460},
        {453, 473, 479, 455}, {480, 474, 454, 456}, {471, 481, 479, 473},
        {480, 482, 472, 474}, {469, 487, 481, 471}, {482, 488, 470, 472},
        {467, 489, 487, 469}, {488, 490, 468, 470}, {465, 491, 489, 467},
        {490, 492, 466, 468}, {391, 389, 503, 501}, {504, 390, 392, 502},
        {393, 391, 501, 499}, {502, 392, 394, 500}, {395, 393, 499, 497},
        {500, 394, 396, 498}, {397, 395, 497, 495}, {498, 396, 398, 496},
        {399, 397, 495, 493}, {496, 398, 400, 494}, {387, 399, 493, 505},
        {494, 400, 388, 506}, {493, 501, 503, 505}, {504, 502, 494, 506},
        {493, 495, 499, 501}, {500, 496, 494, 502}, {313, 381, 387, 505},
        {388, 382, 314, 506}, {313, 505, 503, 321}, {504, 506, 314, 322},
        {319, 321, 503, 389}, {504, 322, 320, 390}};

    positions = suzanne_pos;
    for (auto& p : positions) p *= size / 2;
    quads = suzanne_quads;
    for (auto& t : suzanne_triangles) {
        quads.push_back({t.x, t.y, t.z, t.z});
    }
}

// Watertight cube
template <typename T>
inline void make_box_shape(
    vector<vec4i>& quads, vector<vec<T, 3>>& positions, const vec<T, 3>& size) {
    static auto cube_pos     = vector<vec<T, 3>>{{-1, -1, -1}, {-1, +1, -1},
        {+1, +1, -1}, {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1},
        {+1, -1, +1}};
    static auto cube_quads   = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_quad_uv = vector<vec<T, 2>>{
        {0, 0}, {1, 0}, {1, 1}, {0, 1}};
    positions = cube_pos;
    for (auto& p : positions) p *= size / 2;
    quads = cube_quads;
}

// Generate lines set along a quad.
template <typename T>
inline void make_lines_shape(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vec2i& steps, const vec<T, 2>& size,
    const vec<T, 2>& uvsize, const vec<T, 2>& line_radius) {
    auto nverts = (steps.x + 1) * steps.y;
    auto nlines = steps.x * steps.y;
    auto vid    = [steps](int i, int j) { return j * (steps.x + 1) + i; };
    auto fid    = [steps](int i, int j) { return j * steps.x + i; };

    positions.resize(nverts);
    normals.resize(nverts);
    texturecoords.resize(nverts);
    radius.resize(nverts);
    if (steps.y > 1) {
        for (auto j = 0; j < steps.y; j++) {
            for (auto i = 0; i <= steps.x; i++) {
                auto uv = vec<T, 2>{
                    i / (T)steps.x, j / (T)(steps.y > 1 ? steps.y - 1 : 1)};
                positions[vid(i, j)] = {
                    (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
                normals[vid(i, j)]       = {1, 0, 0};
                texturecoords[vid(i, j)] = uv * uvsize;
            }
        }
    } else {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv                  = vec<T, 2>{i / (T)steps.x, 0};
            positions[vid(i, 0)]     = {(uv.x - 0.5f) * size.x, 0, 0};
            normals[vid(i, 0)]       = {1, 0, 0};
            texturecoords[vid(i, 0)] = uv * uvsize;
        }
    }

    lines.resize(nlines);
    for (int j = 0; j < steps.y; j++) {
        for (int i = 0; i < steps.x; i++) {
            lines[fid(i, j)] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

// Generate a point set with points placed at the origin with texcoords
// varying along u.
template <typename T>
inline void make_points_shape(vector<int>& points, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, int num, T uvsize, T point_radius) {
    points.resize(num);
    for (auto i = 0; i < num; i++) points[i] = i;
    positions.assign(num, {0, 0, 0});
    normals.assign(num, {0, 0, 1});
    texturecoords.assign(num, {0, 0});
    radius.assign(num, point_radius);
    for (auto i = 0; i < texturecoords.size(); i++)
        texturecoords[i] = {(T)i / (T)num, 0};
}

// Generate a point set.
template <typename T>
inline void make_random_points_shape(vector<int>& points,
    vector<vec<T, 3>>& positions, vector<vec<T, 3>>& normals,
    vector<vec<T, 2>>& texturecoords, vector<T>& radius, int num,
    const vec<T, 3>& size, T uvsize, T point_radius, uint64_t seed) {
    make_points_shape(points, positions, normals, texturecoords, radius, num,
        uvsize, point_radius);
    auto rng = make_rng(seed);
    for (auto i = 0; i < positions.size(); i++) {
        positions[i] = (get_random_vec3f(rng) - vec<T, 3>{0.5f, 0.5f, 0.5f}) *
                       size;
    }
}

// Make a point.
template <typename T>
inline void make_point_shape(vector<int>& points, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, T point_radius) {
    points        = {0};
    positions     = {{0, 0, 0}};
    normals       = {{0, 0, 1}};
    texturecoords = {{0, 0}};
    radius        = {point_radius};
}

// Make a bezier circle. Returns bezier, pos.
template <typename T>
inline void make_bezier_circle_shape(
    vector<vec4i>& beziers, vector<vec<T, 3>>& positions, T size) {
    // constant from http://spencermortensen.com/articles/bezier-circle/
    const auto  c          = 0.551915024494f;
    static auto circle_pos = vector<vec<T, 3>>{{1, 0, 0}, {1, c, 0}, {c, 1, 0},
        {0, 1, 0}, {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
        {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
    static auto circle_beziers = vector<vec4i>{
        {0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
    positions = circle_pos;
    beziers   = circle_beziers;
    for (auto& p : positions) p *= size;
}

// Make a hair ball around a shape
template <typename T>
inline void make_hair_shape(vector<vec2i>& lines, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const vec2i& steps, const vector<vec3i>& striangles,
    const vector<vec4i>& squads, const vector<vec<T, 3>>& spos,
    const vector<vec<T, 3>>& snorm, const vector<vec<T, 2>>& stexcoord,
    const vec<T, 2>& len, const vec<T, 2>& rad, const vec<T, 2>& noise,
    const vec<T, 2>& clump, const vec<T, 2>& rotation, int seed) {
    auto alltriangles    = striangles;
    auto quads_triangles = vector<vec3i>{};
    convert_quads_to_triangles(quads_triangles, squads);
    alltriangles.insert(
        alltriangles.end(), quads_triangles.begin(), quads_triangles.end());
    auto bpos      = vector<vec<T, 3>>{};
    auto bnorm     = vector<vec<T, 3>>{};
    auto btexcoord = vector<vec<T, 2>>{};
    sample_triangles_points(bpos, bnorm, btexcoord, alltriangles, spos, snorm,
        stexcoord, steps.y, seed);

    auto rng  = make_rng(seed, 3);
    auto blen = vector<T>(bpos.size());
    for (auto& l : blen) l = lerp(len.x, len.y, (T)get_random_float(rng));

    auto cidx = vector<int>();
    if (clump.x > 0) {
        for (auto bidx = 0; bidx < bpos.size(); bidx++) {
            cidx.push_back(0);
            auto cdist = type_max<T>;
            for (auto c = 0; c < clump.y; c++) {
                auto d = length(bpos[bidx] - bpos[c]);
                if (d < cdist) {
                    cdist       = d;
                    cidx.back() = c;
                }
            }
        }
    }

    make_lines_shape(lines, positions, normals, texturecoords, radius, steps,
        {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto u       = texturecoords[i].x;
        auto bidx    = i / (steps.x + 1);
        positions[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
        normals[i]   = bnorm[bidx];
        radius[i]    = lerp(rad.x, rad.y, u);
        if (clump.x > 0) {
            positions[i] = positions[i] +
                           (positions[i + (cidx[bidx] - bidx) * (steps.x + 1)] -
                               positions[i]) *
                               u * clump.x;
        }
        if (noise.x > 0) {
            auto nx = perlin_noise(
                          positions[i] * noise.y + vec<T, 3>{0, 0, 0}) *
                      noise.x;
            auto ny = perlin_noise(
                          positions[i] * noise.y + vec<T, 3>{3, 7, 11}) *
                      noise.x;
            auto nz = perlin_noise(
                          positions[i] * noise.y + vec<T, 3>{13, 17, 19}) *
                      noise.x;
            positions[i] += {nx, ny, nz};
        }
    }

    if (clump.x > 0 || noise.x > 0 || rotation.x > 0) {
        compute_vertex_tangents(normals, lines, positions);
    }
}

// Thickens a shape by copy9ing the shape content, rescaling it and flipping its
// normals. Note that this is very much not robust and only useful for trivial
// cases.
template <typename T>
inline void make_shell_shape(vector<vec4i>& quads, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords, T thickness) {
    auto bbox = invalid_bbox3f;
    for (auto p : positions) bbox += p;
    auto center              = bbox_center(bbox);
    auto inner_quads         = quads;
    auto inner_positions     = positions;
    auto inner_normals       = normals;
    auto inner_texturecoords = texturecoords;
    for (auto& p : inner_positions) p = (1 - thickness) * (p - center) + center;
    for (auto& n : inner_normals) n = -n;
    merge_quads(quads, positions, normals, texturecoords, inner_quads,
        inner_positions, inner_normals, inner_texturecoords);
}

// Shape presets used ofr testing.
template <typename T>
inline void make_shape_presets(vector<vec4i>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec<T, 3>>& positions,
    vector<vec<T, 3>>& normals, vector<vec<T, 2>>& texturecoords,
    vector<T>& radius, const string& type) {
    if (type == "default-quad") {
        make_quad_shape(
            quads, positions, normals, texturecoords, {1, 1}, {2, 2}, {1, 1});
    } else if (type == "default-quady") {
        make_quad_shape(
            quads, positions, normals, texturecoords, {1, 1}, {2, 2}, {1, 1});
    } else if (type == "default-quad-stack") {
        make_quad_stack_shape(quads, positions, normals, texturecoords,
            {1, 1, 1}, {2, 2, 2}, {1, 1});
    } else if (type == "default-box") {
        make_box_shape(quads, positions, normals, texturecoords, {1, 1, 1},
            {2, 2, 2}, {1, 1, 1});
    } else if (type == "default-box-rounded") {
        make_box_rounded_shape(quads, positions, normals, texturecoords,
            {32, 32, 32}, {2, 2, 2}, {1, 1, 1}, 0.15f);
    } else if (type == "default-uvsphere") {
        make_uvsphere_shape(
            quads, positions, normals, texturecoords, {64, 32}, 2, {1, 1});
    } else if (type == "default-sphere") {
        make_sphere_shape(quads, positions, normals, texturecoords, 32, 2, 1);
    } else if (type == "default-uvsphere-flipcap") {
        make_uvsphere_flipcap_shape(quads, positions, normals, texturecoords,
            {64, 32}, 2, {1, 1}, {(T)-0.75, (T) + 0.75});
    } else if (type == "default-uvdisk") {
        make_uvdisk_shape(
            quads, positions, normals, texturecoords, {32, 16}, 2, {1, 1});
    } else if (type == "default-disk") {
        make_disk_shape(quads, positions, normals, texturecoords, 32, 2, 1);
    } else if (type == "default-disk-bulged") {
        make_disk_bulged_shape(
            quads, positions, normals, texturecoords, 32, 2, 1, (T)0.25);
    } else if (type == "default-quad-bulged") {
        make_quad_bulged_shape(
            quads, positions, normals, texturecoords, 32, 2, 1, (T)0.25);
    } else if (type == "default-uvcylinder") {
        make_uvcylinder_shape(quads, positions, normals, texturecoords,
            {64, 32, 16}, {2, 2}, {1, 1, 1});
    } else if (type == "default-uvcylinder-rounded") {
        make_uvcylinder_rounded_shape(quads, positions, normals, texturecoords,
            {64, 32, 16}, {2, 2}, {1, 1, 1}, (T)0.075);
    } else if (type == "default-sphere-geodesic") {
        make_geodesic_sphere_shape(triangles, positions, normals, 4, 2);
    } else if (type == "default-floor") {
        make_floor_shape(quads, positions, normals, texturecoords, {1, 1},
            {40, 40}, {20, 20});
    } else if (type == "default-floor-bent") {
        make_floor_bent_shape(quads, positions, normals, texturecoords, {1, 40},
            {40, 40}, {20, 20}, 10);
    } else if (type == "default-matball") {
        make_sphere_shape(quads, positions, normals, texturecoords, 32, 2, 1);
    } else if (type == "default-hairball") {
        auto base_quads         = vector<vec4i>{};
        auto base_positions     = vector<vec3f>{};
        auto base_normals       = vector<vec3f>{};
        auto base_texturecoords = vector<vec2f>{};
        make_sphere_shape(base_quads, base_positions, base_normals,
            base_texturecoords, 32, 2 * (T)0.8, 1);
        make_hair_shape(lines, positions, normals, texturecoords, radius,
            {4, 65536}, {}, base_quads, base_positions, base_normals,
            base_texturecoords, {(T)0.2, (T)0.2}, {(T)0.002, (T)0.001}, {0, 0},
            {0, 0});
    } else if (type == "default-hairball-interior") {
        make_sphere_shape(
            quads, positions, normals, texturecoords, 32, 2 * (T)0.8, 1);
    } else if (type == "default-suzanne") {
        make_suzanne_shape(quads, positions, 2);
    } else if (type == "default-cube-posonly") {
        auto ignore1 = vector<vec4i>{};
        auto ignore2 = vector<vec4i>{};
        auto ignore3 = vector<vec3f>{};
        auto ignore4 = vector<vec2f>{};
        make_box_fvshape(quads, ignore1, ignore2, positions, ignore3, ignore4,
            {1, 1, 1}, {2, 2, 2}, {1, 1, 1});
    } else if (type == "default-cube-facevarying") {
        make_box_fvshape(quads_positions, quads_normals, quads_texturecoords,
            positions, normals, texturecoords, {1, 1, 1}, {2, 2, 2}, {1, 1, 1});
    } else if (type == "default-sphere-facevarying") {
        make_sphere_fvshape(quads_positions, quads_normals, quads_texturecoords,
            positions, normals, texturecoords, 32, 2.0f, 1.0f);
    } else if (type == "test-cube") {
        make_box_rounded_shape(quads, positions, normals, texturecoords,
            {32, 32, 32}, {2, 2, 2}, {1, 1, 1}, 0.15f);
    } else if (type == "test-uvsphere") {
        make_uvsphere_shape(
            quads, positions, normals, texturecoords, {64, 32}, 2.0f, {1, 1});
    } else if (type == "test-sphere") {
        make_sphere_shape(
            quads, positions, normals, texturecoords, 32, 2.0f, 1.0f);
    } else if (type == "test-disk") {
        make_disk_shape(
            quads, positions, normals, texturecoords, 32, 2.0f, 1.0f);
    } else if (type == "test-cylinder") {
        make_uvcylinder_rounded_shape(quads, positions, normals, texturecoords,
            {64, 32, 16}, vec2f{2.0f, 2.0f}, {1, 1, 1}, 0.075f);
    } else if (type == "test-floor") {
        make_floor_shape(quads, positions, normals, texturecoords, {1, 1},
            {40, 40}, {20, 20});
    } else if (type == "test-matball") {
        make_sphere_shape(
            quads, positions, normals, texturecoords, 32, 2.0f, 1.0f);
    } else if (type == "test-hairball1") {
        auto base_quads         = vector<vec4i>{};
        auto base_positions     = vector<vec3f>{};
        auto base_normals       = vector<vec3f>{};
        auto base_texturecoords = vector<vec2f>{};
        make_sphere_shape(base_quads, base_positions, base_normals,
            base_texturecoords, 32, 2.0f * 0.8f, 1.0f);
        make_hair_shape(lines, positions, normals, texturecoords, radius,
            {4, 65536}, {}, base_quads, base_positions, base_normals,
            base_texturecoords, {2 * (T)0.1, 2 * (T)0.1},
            {2 * (T)0.001, 2 * (T)0.0005}, {0, 0}, {0, 0});
    } else if (type == "test-hairball2") {
        auto base_quads         = vector<vec4i>{};
        auto base_positions     = vector<vec3f>{};
        auto base_normals       = vector<vec3f>{};
        auto base_texturecoords = vector<vec2f>{};
        make_sphere_shape(base_quads, base_positions, base_normals,
            base_texturecoords, 32, 2.0f * 0.8f, 1.0f);
        make_hair_shape(lines, positions, normals, texturecoords, radius,
            {4, 65536}, {}, base_quads, base_positions, base_normals,
            base_texturecoords, {2 * (T)0.1, 2 * (T)0.1},
            {2 * (T)0.001, 2 * (T)0.0005}, {0, 0}, {0, 0});
    } else if (type == "test-hairball3") {
        auto base_quads         = vector<vec4i>{};
        auto base_positions     = vector<vec3f>{};
        auto base_normals       = vector<vec3f>{};
        auto base_texturecoords = vector<vec2f>{};
        make_sphere_shape(base_quads, base_positions, base_normals,
            base_texturecoords, 32, 2.0f * 0.8f, 1.0f);
        make_hair_shape(lines, positions, normals, texturecoords, radius,
            {4, 65536}, {}, base_quads, base_positions, base_normals,
            base_texturecoords, vec2f{0.1f, 0.1f} * 2.0f,
            vec2f{0.001f, 0.0005f} * 2.0f, {0, 0}, {0, 0});
    } else if (type == "test-hairball-interior") {
        make_sphere_shape(
            quads, positions, normals, texturecoords, 32, 2.0f * 0.8f, 1.0f);
    } else if (type == "test-suzanne-subdiv") {
        make_suzanne_shape(quads, positions, 2.0f);
    } else if (type == "test-cube-subdiv") {
        make_box_fvshape(quads_positions, quads_normals, quads_texturecoords,
            positions, normals, texturecoords, {1, 1, 1}, {2, 2, 2}, {1, 1, 1});
    } else {
        throw std::invalid_argument("unknown procedural type " + type);
    }
}

}  // namespace yocto

#include "ext/happly.h"
#include "yocto_obj.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace yocto {

// merge quads and triangles
inline void merge_triangles_and_quads(
    vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles) {
    if (quads.empty()) return;
    if (force_triangles) {
        auto qtriangles = vector<vec3i>{};
        convert_quads_to_triangles(qtriangles, quads);
        triangles.insert(triangles.end(), qtriangles.begin(), qtriangles.end());
        quads = {};
    } else {
        auto tquads = vector<vec4i>{};
        convert_triangles_to_quads(tquads, triangles);
        quads.insert(quads.end(), tquads.begin(), tquads.end());
        triangles = {};
    }
}

// hack for CyHair data
inline void load_cyhair_shape(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<vec4f>& color, vector<float>& radius,
    bool flip_texcoord = true);

// Load/Save a ply mesh
inline void load_ply_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& color,
    vector<float>& radius, bool flip_texcoord = true);
inline void save_ply_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius,
    bool ascii = false, bool flip_texcoord = true);

// Load/Save an OBJ mesh
inline void load_obj_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    bool preserve_facevarying, bool flip_texcoord = true);
inline void save_obj_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord = true);

// Load ply mesh
inline void load_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool preserve_facevarrying) {
    points              = {};
    lines               = {};
    triangles           = {};
    quads               = {};
    quads_positions     = {};
    quads_normals       = {};
    quads_texturecoords = {};
    positions           = {};
    normals             = {};
    texturecoords       = {};
    colors              = {};
    radius              = {};

    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        load_ply_shape(filename, points, lines, triangles, quads,
            quads_positions, quads_normals, quads_texturecoords, positions,
            normals, texturecoords, colors, radius);
    } else if (ext == "obj" || ext == "OBJ") {
        load_obj_shape(filename, points, lines, triangles, quads,
            quads_positions, quads_normals, quads_texturecoords, positions,
            normals, texturecoords, preserve_facevarrying);
    } else if (ext == "hair" || ext == "HAIR") {
        load_cyhair_shape(
            filename, lines, positions, normals, texturecoords, colors, radius);
    } else {
        throw shapeio_error("unsupported mesh type " + ext);
    }
}

// Save ply mesh
inline void save_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii) {
    auto ext = get_extension(filename);
    if (ext == "ply" || ext == "PLY") {
        return save_ply_shape(filename, points, lines, triangles, quads,
            quads_positions, quads_normals, quads_texturecoords, positions,
            normals, texturecoords, colors, radius, ascii);
    } else if (ext == "obj" || ext == "OBJ") {
        return save_obj_shape(filename, points, lines, triangles, quads,
            quads_positions, quads_normals, quads_texturecoords, positions,
            normals, texturecoords);
    } else {
        throw shapeio_error("unsupported mesh type " + ext);
    }
}

inline void load_ply_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool flip_texcoord) {
    try {
        // load ply
        happly::PLYData ply(filename);

        // copy vertex data
        if (ply.hasElement("vertex")) {
            auto& vertex = ply.getElement("vertex");
            if (vertex.hasProperty("x") && vertex.hasProperty("y") &&
                vertex.hasProperty("z")) {
                auto x = vertex.getProperty<float>("x");
                auto y = vertex.getProperty<float>("y");
                auto z = vertex.getProperty<float>("z");
                positions.resize(x.size());
                for (auto i = 0; i < positions.size(); i++) {
                    positions[i] = {x[i], y[i], z[i]};
                }
            } else {
                throw shapeio_error("vertex positions not present");
            }
            if (vertex.hasProperty("nx") && vertex.hasProperty("ny") &&
                vertex.hasProperty("nz")) {
                auto x = vertex.getProperty<float>("nx");
                auto y = vertex.getProperty<float>("ny");
                auto z = vertex.getProperty<float>("nz");
                normals.resize(x.size());
                for (auto i = 0; i < normals.size(); i++) {
                    normals[i] = {x[i], y[i], z[i]};
                }
            }
            if (vertex.hasProperty("u") && vertex.hasProperty("v")) {
                auto x = vertex.getProperty<float>("u");
                auto y = vertex.getProperty<float>("v");
                texturecoords.resize(x.size());
                for (auto i = 0; i < texturecoords.size(); i++) {
                    texturecoords[i] = {x[i], y[i]};
                }
            }
            if (vertex.hasProperty("s") && vertex.hasProperty("t")) {
                auto x = vertex.getProperty<float>("s");
                auto y = vertex.getProperty<float>("t");
                texturecoords.resize(x.size());
                for (auto i = 0; i < texturecoords.size(); i++) {
                    texturecoords[i] = {x[i], y[i]};
                }
            }
            if (vertex.hasProperty("red") && vertex.hasProperty("green") &&
                vertex.hasProperty("blue")) {
                auto x = vertex.getProperty<float>("red");
                auto y = vertex.getProperty<float>("green");
                auto z = vertex.getProperty<float>("blue");
                colors.resize(x.size());
                for (auto i = 0; i < colors.size(); i++) {
                    colors[i] = {x[i], y[i], z[i], 1};
                }
                if (vertex.hasProperty("alpha")) {
                    auto w = vertex.getProperty<float>("alpha");
                    for (auto i = 0; i < colors.size(); i++) {
                        colors[i].w = w[i];
                    }
                }
            }
            if (vertex.hasProperty("radius")) {
                radius = vertex.getProperty<float>("radius");
            }
        }

        // fix texture coordinated
        if (flip_texcoord && !texturecoords.empty()) {
            for (auto& uv : texturecoords) uv.y = 1 - uv.y;
        }

        // copy face data
        if (ply.hasElement("face")) {
            auto& elements = ply.getElement("face");
            if (!elements.hasProperty("vertex_indices"))
                throw shapeio_error("bad ply faces");
            auto indices = vector<vector<int>>{};
            try {
                indices = elements.getListProperty<int>("vertex_indices");
            } catch (...) {
                (vector<vector<unsigned int>>&)indices =
                    elements.getListProperty<unsigned int>("vertex_indices");
            }
            for (auto& face : indices) {
                if (face.size() == 4) {
                    quads.push_back({face[0], face[1], face[2], face[3]});
                } else {
                    for (auto i = 2; i < face.size(); i++)
                        triangles.push_back({face[0], face[i - 1], face[i]});
                }
            }
        }

        // copy face data
        if (ply.hasElement("line")) {
            auto& elements = ply.getElement("line");
            if (!elements.hasProperty("vertex_indices"))
                throw shapeio_error("bad ply lines");
            auto indices = vector<vector<int>>{};
            try {
                indices = elements.getListProperty<int>("vertex_indices");
            } catch (...) {
                (vector<vector<unsigned int>>&)indices =
                    elements.getListProperty<unsigned int>("vertex_indices");
            }
            for (auto& line : indices) {
                for (auto i = 1; i < line.size(); i++)
                    lines.push_back({line[i], line[i - 1]});
            }
        }

        merge_triangles_and_quads(triangles, quads, false);

    } catch (const std::exception& e) {
        throw shapeio_error("cannot load mesh " + filename + "\n" + e.what());
    }
}

// Save ply mesh
inline void save_ply_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii,
    bool flip_texcoord) {
    if (!quads_positions.empty()) {
        auto split_quads         = vector<vec4i>{};
        auto split_positions     = vector<vec3f>{};
        auto split_normals       = vector<vec3f>{};
        auto split_texturecoords = vector<vec2f>{};
        convert_facevarying(split_quads, split_positions, split_normals,
            split_texturecoords, quads_positions, quads_normals,
            quads_texturecoords, positions, normals, texturecoords);
        return save_ply_shape(filename, {}, {}, {}, split_quads, {}, {}, {},
            positions, normals, texturecoords, {}, {}, ascii, flip_texcoord);
    }

    // empty data
    happly::PLYData ply;

    // add elements
    ply.addElement("vertex", positions.size());
    if (!positions.empty()) {
        auto& vertex = ply.getElement("vertex");
        auto  x      = vector<float>{};
        auto  y      = vector<float>{};
        auto  z      = vector<float>{};
        for (auto& p : positions) {
            x.push_back(p.x);
            y.push_back(p.y);
            z.push_back(p.z);
        }
        vertex.addProperty("x", x);
        vertex.addProperty("y", y);
        vertex.addProperty("z", z);
    }
    if (!normals.empty()) {
        auto& vertex = ply.getElement("vertex");
        auto  x      = vector<float>{};
        auto  y      = vector<float>{};
        auto  z      = vector<float>{};
        for (auto& n : normals) {
            x.push_back(n.x);
            y.push_back(n.y);
            z.push_back(n.z);
        }
        vertex.addProperty("nx", x);
        vertex.addProperty("ny", y);
        vertex.addProperty("nz", z);
    }
    if (!texturecoords.empty()) {
        auto& vertex = ply.getElement("vertex");
        auto  x      = vector<float>{};
        auto  y      = vector<float>{};
        for (auto& t : texturecoords) {
            x.push_back(t.x);
            y.push_back(flip_texcoord ? 1 - t.y : t.y);
        }
        vertex.addProperty("u", x);
        vertex.addProperty("v", y);
    }
    if (!colors.empty()) {
        auto& vertex = ply.getElement("vertex");
        auto  x      = vector<float>{};
        auto  y      = vector<float>{};
        auto  z      = vector<float>{};
        auto  w      = vector<float>{};
        for (auto& c : colors) {
            x.push_back(c.x);
            y.push_back(c.y);
            z.push_back(c.z);
            w.push_back(c.w);
        }
        vertex.addProperty("red", x);
        vertex.addProperty("green", y);
        vertex.addProperty("blue", z);
        vertex.addProperty("alpha", w);
    }
    if (!radius.empty()) {
        auto& vertex = ply.getElement("vertex");
        vertex.addProperty("radius", radius);
    }

    // face date
    if (!triangles.empty() || !quads.empty()) {
        ply.addElement("face", triangles.size() + quads.size());
        auto elements = vector<vector<int>>{};
        for (auto& t : triangles) {
            elements.push_back({t.x, t.y, t.z});
        }
        for (auto& q : quads) {
            if (q.z == q.w) {
                elements.push_back({q.x, q.y, q.z});
            } else {
                elements.push_back({q.x, q.y, q.z, q.w});
            }
        }
        ply.getElement("face").addListProperty("vertex_indices", elements);
    }
    if (!lines.empty()) {
        ply.addElement("line", lines.size());
        auto elements = vector<vector<int>>{};
        for (auto& l : lines) {
            elements.push_back({l.x, l.y});
        }
        ply.getElement("line").addListProperty("vertex_indices", elements);
    }
    if (!points.empty() || !quads.empty()) {
        ply.addElement("point", points.size());
        auto elements = vector<vector<int>>{};
        for (auto& p : points) {
            elements.push_back({p});
        }
        ply.getElement("point").addListProperty("vertex_indices", elements);
    }

    // Write our data
    try {
        ply.write(filename,
            ascii ? happly::DataFormat::ASCII : happly::DataFormat::Binary);
    } catch (const std::exception& e) {
        throw shapeio_error("cannot save mesh " + filename + "\n" + e.what());
    }
}

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.position == b.position && a.texturecoord == b.texturecoord &&
           a.normal == b.normal;
}

struct obj_vertex_hash {
    size_t operator()(const obj_vertex& v) const {
        auto vh = std::hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.position)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

inline void print_value(const output_file& fs, const obj_vertex& value) {
    print_value(fs, value.position);
    if (value.texturecoord) {
        print_value(fs, '/');
        print_value(fs, value.texturecoord);
        if (value.normal) {
            print_value(fs, '/');
            print_value(fs, value.normal);
        }
    } else {
        if (value.normal) {
            print_value(fs, '/');
            print_value(fs, '/');
            print_value(fs, value.normal);
        }
    }
}

// Load ply mesh
inline void load_obj_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    bool preserve_facevarying, bool flip_texcoord) {
    struct parse_callbacks : obj_callbacks {
        vector<int>&   points;
        vector<vec2i>& lines;
        vector<vec3i>& triangles;
        vector<vec4i>& quads;
        vector<vec4i>& quads_positions;
        vector<vec4i>& quads_normals;
        vector<vec4i>& quads_texturecoords;
        vector<vec3f>& positions;
        vector<vec3f>& normals;
        vector<vec2f>& texturecoords;
        bool           facevarying = false;

        // TODO: implement me

        // obj vertices
        std::deque<vec3f> opos      = std::deque<vec3f>();
        std::deque<vec3f> onorm     = std::deque<vec3f>();
        std::deque<vec2f> otexcoord = std::deque<vec2f>();

        // vertex maps
        unordered_map<obj_vertex, int, obj_vertex_hash> vertex_map =
            unordered_map<obj_vertex, int, obj_vertex_hash>();

        // vertex maps
        unordered_map<int, int> pos_map      = unordered_map<int, int>();
        unordered_map<int, int> texcoord_map = unordered_map<int, int>();
        unordered_map<int, int> norm_map     = unordered_map<int, int>();

        parse_callbacks(vector<int>& points, vector<vec2i>& lines,
            vector<vec3i>& triangles, vector<vec4i>& quads,
            vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
            vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
            vector<vec3f>& normals, vector<vec2f>& texturecoords,
            bool facevarying)
            : points{points}
            , lines{lines}
            , triangles{triangles}
            , quads{quads}
            , quads_positions{quads_positions}
            , quads_normals{quads_normals}
            , quads_texturecoords{quads_texturecoords}
            , positions{positions}
            , normals{normals}
            , texturecoords{texturecoords}
            , facevarying{facevarying} {}

        // Add  vertices to the current shape
        void add_verts(const vector<obj_vertex>& verts) {
            for (auto& vert : verts) {
                auto it = vertex_map.find(vert);
                if (it != vertex_map.end()) continue;
                auto nverts = (int)positions.size();
                vertex_map.insert(it, {vert, nverts});
                if (vert.position)
                    positions.push_back(opos.at(vert.position - 1));
                if (vert.texturecoord)
                    texturecoords.push_back(
                        otexcoord.at(vert.texturecoord - 1));
                if (vert.normal) normals.push_back(onorm.at(vert.normal - 1));
            }
        }

        // add vertex
        void add_fvverts(const vector<obj_vertex>& verts) {
            for (auto& vert : verts) {
                if (!vert.position) continue;
                auto pos_it = pos_map.find(vert.position);
                if (pos_it != pos_map.end()) continue;
                auto nverts = (int)positions.size();
                pos_map.insert(pos_it, {vert.position, nverts});
                positions.push_back(opos.at(vert.position - 1));
            }
            for (auto& vert : verts) {
                if (!vert.texturecoord) continue;
                auto texcoord_it = texcoord_map.find(vert.texturecoord);
                if (texcoord_it != texcoord_map.end()) continue;
                auto nverts = (int)texturecoords.size();
                texcoord_map.insert(texcoord_it, {vert.texturecoord, nverts});
                texturecoords.push_back(otexcoord.at(vert.texturecoord - 1));
            }
            for (auto& vert : verts) {
                if (!vert.normal) continue;
                auto norm_it = norm_map.find(vert.normal);
                if (norm_it != norm_map.end()) continue;
                auto nverts = (int)normals.size();
                norm_map.insert(norm_it, {vert.normal, nverts});
                normals.push_back(onorm.at(vert.normal - 1));
            }
        }

        void vert(const vec3f& v) { opos.push_back(v); }
        void norm(const vec3f& v) { onorm.push_back(v); }
        void texcoord(const vec2f& v) { otexcoord.push_back(v); }
        void face(const vector<obj_vertex>& verts) {
            if (!facevarying) {
                add_verts(verts);
                if (verts.size() == 4) {
                    quads.push_back(
                        {vertex_map.at(verts[0]), vertex_map.at(verts[1]),
                            vertex_map.at(verts[2]), vertex_map.at(verts[3])});
                } else {
                    for (auto i = 2; i < verts.size(); i++)
                        triangles.push_back({vertex_map.at(verts[0]),
                            vertex_map.at(verts[i - 1]),
                            vertex_map.at(verts[i])});
                }
            } else {
                add_fvverts(verts);
                if (verts.size() == 4) {
                    if (verts[0].position) {
                        quads_positions.push_back(
                            {pos_map.at(verts[0].position),
                                pos_map.at(verts[1].position),
                                pos_map.at(verts[2].position),
                                pos_map.at(verts[3].position)});
                    }
                    if (verts[0].texturecoord) {
                        quads_texturecoords.push_back(
                            {texcoord_map.at(verts[0].texturecoord),
                                texcoord_map.at(verts[1].texturecoord),
                                texcoord_map.at(verts[2].texturecoord),
                                texcoord_map.at(verts[3].texturecoord)});
                    }
                    if (verts[0].normal) {
                        quads_normals.push_back({norm_map.at(verts[0].normal),
                            norm_map.at(verts[1].normal),
                            norm_map.at(verts[2].normal),
                            norm_map.at(verts[3].normal)});
                    }
                    // quads_materials.push_back(current_material_id);
                } else {
                    if (verts[0].position) {
                        for (auto i = 2; i < verts.size(); i++)
                            quads_positions.push_back(
                                {pos_map.at(verts[0].position),
                                    pos_map.at(verts[1].position),
                                    pos_map.at(verts[i].position),
                                    pos_map.at(verts[i].position)});
                    }
                    if (verts[0].texturecoord) {
                        for (auto i = 2; i < verts.size(); i++)
                            quads_texturecoords.push_back(
                                {texcoord_map.at(verts[0].texturecoord),
                                    texcoord_map.at(verts[1].texturecoord),
                                    texcoord_map.at(verts[i].texturecoord),
                                    texcoord_map.at(verts[i].texturecoord)});
                    }
                    if (verts[0].normal) {
                        for (auto i = 2; i < verts.size(); i++)
                            quads_normals.push_back(
                                {norm_map.at(verts[0].normal),
                                    norm_map.at(verts[1].normal),
                                    norm_map.at(verts[i].normal),
                                    norm_map.at(verts[i].normal)});
                    }
                    //                    for (auto i = 2; i < verts.size();
                    //                    i++)
                    //                        quads_materials.push_back(current_material_id);
                }
            }
        }
        void line(const vector<obj_vertex>& verts) {
            add_verts(verts);
            for (auto i = 1; i < verts.size(); i++)
                lines.push_back(
                    {vertex_map.at(verts[i - 1]), vertex_map.at(verts[i])});
        }
        void point(const vector<obj_vertex>& verts) {
            add_verts(verts);
            for (auto i = 0; i < verts.size(); i++)
                points.push_back(vertex_map.at(verts[i]));
        }
        //        void usemtl(const string& name) {
        //            auto pos = std::find(
        //                                 material_group.begin(),
        //                                 material_group.end(), name);
        //            if (pos == material_group.end()) {
        //                material_group.push_back(name);
        //                current_material_id = (int)material_group.size() - 1;
        //            } else {
        //                current_material_id = (int)(pos -
        //                material_group.begin());
        //            }
        //        }
    };

    try {
        // load obj
        auto obj_options          = load_obj_options();
        obj_options.exit_on_error = false;
        obj_options.geometry_only = true;
        obj_options.flip_texcoord = flip_texcoord;
        auto cb = parse_callbacks{points, lines, triangles, quads,
            quads_positions, quads_normals, quads_texturecoords, positions,
            normals, texturecoords, preserve_facevarying};
        load_obj(filename, cb, obj_options);

        // merging quads and triangles
        if (!preserve_facevarying) {
            merge_triangles_and_quads(triangles, quads, false);
        }

    } catch (const std::exception& e) {
        throw shapeio_error("cannot load mesh " + filename + "\n" + e.what());
    }
}

// Load ply mesh
inline void save_obj_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quads_positions,
    const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord) {
    auto fs = output_file(filename);

    // TODO: implement me

    println_values(
        fs, "# Saved by Yocto/GL - https://github.com/xelatihy/yocto-gl\n");

    for (auto& p : positions) println_values(fs, "v", p);
    for (auto& n : normals) println_values(fs, "vn", n);
    for (auto& t : texturecoords)
        println_values(fs, "vt", vec2f{t.x, (flip_texcoord) ? 1 - t.y : t.y});

    auto mask = obj_vertex{
        1, texturecoords.empty() ? 0 : 1, normals.empty() ? 0 : 1};
    auto vert = [mask](int i) {
        return obj_vertex{(i + 1) * mask.position, (i + 1) * mask.texturecoord,
            (i + 1) * mask.normal};
    };

    for (auto& p : points) {
        println_values(fs, "p", vert(p));
    }
    for (auto& l : lines) {
        println_values(fs, "l", vert(l.x), vert(l.y));
    }
    for (auto& t : triangles) {
        println_values(fs, "f", vert(t.x), vert(t.y), vert(t.z));
    }
    for (auto& q : quads) {
        if (q.z == q.w) {
            println_values(fs, "f", vert(q.x), vert(q.y), vert(q.z));
        } else {
            println_values(fs, "f", vert(q.x), vert(q.y), vert(q.z), vert(q.w));
        }
    }

    auto fvmask = obj_vertex{
        1, texturecoords.empty() ? 0 : 1, normals.empty() ? 0 : 1};
    auto fvvert = [fvmask](int pi, int ti, int ni) {
        return obj_vertex{(pi + 1) * fvmask.position,
            (ti + 1) * fvmask.texturecoord, (ni + 1) * fvmask.normal};
    };

    // auto last_material_id = -1;
    for (auto i = 0; i < quads_positions.size(); i++) {
        //        if (!quads_materials.empty() &&
        //            quads_materials[i] != last_material_id) {
        //            last_material_id = quads_materials[i];
        //            println_values(fs, "usemtl material_{}\n",
        //            last_material_id);
        //        }
        auto qp = quads_positions.at(i);
        auto qt = !quads_texturecoords.empty() ? quads_texturecoords.at(i)
                                               : vec4i{-1, -1, -1, -1};
        auto qn = !quads_normals.empty() ? quads_normals.at(i)
                                         : vec4i{-1, -1, -1, -1};
        if (qp.z != qp.w) {
            println_values(fs, "f", fvvert(qp.x, qt.x, qn.x),
                fvvert(qp.y, qt.y, qn.y), fvvert(qp.z, qt.z, qn.z),
                fvvert(qp.w, qt.w, qn.w));
        } else {
            println_values(fs, "f", fvvert(qp.x, qt.x, qn.x),
                fvvert(qp.y, qt.y, qn.y), fvvert(qp.z, qt.z, qn.z));
        }
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CYHAIR
// -----------------------------------------------------------------------------
namespace yocto {

struct cyhair_strand {
    vector<vec3f> positions;
    vector<float> radius;
    vector<float> transparency;
    vector<vec3f> color;
};

struct cyhair_data {
    vector<cyhair_strand> strands              = {};
    float                 default_thickness    = 0;
    float                 default_transparency = 0;
    vec3f                 default_color        = zero3f;
};

inline void load_cyhair(const string& filename, cyhair_data& hair) {
    // open file
    hair    = {};
    auto fs = input_file(filename, std::ios::binary);

    // Bytes 0-3    Must be "HAIR" in ascii code (48 41 49 52)
    // Bytes 4-7    Number of hair strands as unsigned int
    // Bytes 8-11    Total number of points of all strands as unsigned int
    // Bytes 12-15    Bit array of data in the file
    // Bit-0 is 1 if the file has segments array.
    // Bit-1 is 1 if the file has points array (this bit must be 1).
    // Bit-2 is 1 if the file has radius array.
    // Bit-3 is 1 if the file has transparency array.
    // Bit-4 is 1 if the file has color array.
    // Bit-5 to Bit-31 are reserved for future extension (must be 0).
    // Bytes 16-19    Default number of segments of hair strands as unsigned int
    // If the file does not have a segments array, this default value is used.
    // Bytes 20-23    Default radius hair strands as float
    // If the file does not have a radius array, this default value is used.
    // Bytes 24-27    Default transparency hair strands as float
    // If the file does not have a transparency array, this default value is
    // used. Bytes 28-39    Default color hair strands as float array of size 3
    // If the file does not have a radius array, this default value is used.
    // Bytes 40-127    File information as char array of size 88 in ascii

    // parse header
    hair = cyhair_data{};
    struct cyhair_header {
        char         magic[4]             = {0};
        unsigned int num_strands          = 0;
        unsigned int num_points           = 0;
        unsigned int flags                = 0;
        unsigned int default_segments     = 0;
        float        default_thickness    = 0;
        float        default_transparency = 0;
        vec3f        default_color        = zero3f;
        char         info[88]             = {0};
    };
    static_assert(sizeof(cyhair_header) == 128);
    auto header = cyhair_header{};
    read_value(fs, header);
    if (header.magic[0] != 'H' || header.magic[1] != 'A' ||
        header.magic[2] != 'I' || header.magic[3] != 'R')
        throw shapeio_error("bad cyhair header");

    // set up data
    hair.default_thickness    = header.default_thickness;
    hair.default_transparency = header.default_transparency;
    hair.default_color        = header.default_color;
    hair.strands.resize(header.num_strands);

    // get segments length
    auto segments = vector<unsigned short>();
    if (header.flags & 1) {
        segments.resize(header.num_strands);
        read_values(fs, segments);
    } else {
        segments.assign(header.num_strands, header.default_segments);
    }

    // check segment length
    auto total_length = 0;
    for (auto segment : segments) total_length += segment + 1;
    if (total_length != header.num_points) {
        throw shapeio_error("bad cyhair file");
    }

    // read positions data
    if (header.flags & 2) {
        for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
            auto strand_size = (int)segments[strand_id] + 1;
            hair.strands[strand_id].positions.resize(strand_size);
            read_values(fs, hair.strands[strand_id].positions);
        }
    }
    // read radius data
    if (header.flags & 4) {
        for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
            auto strand_size = (int)segments[strand_id] + 1;
            hair.strands[strand_id].radius.resize(strand_size);
            read_values(fs, hair.strands[strand_id].radius);
        }
    }
    // read transparency data
    if (header.flags & 8) {
        for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
            auto strand_size = (int)segments[strand_id] + 1;
            hair.strands[strand_id].transparency.resize(strand_size);
            read_values(fs, hair.strands[strand_id].transparency);
        }
    }
    // read color data
    if (header.flags & 16) {
        for (auto strand_id = 0; strand_id < header.num_strands; strand_id++) {
            auto strand_size = (int)segments[strand_id] + 1;
            hair.strands[strand_id].color.resize(strand_size);
            read_values(fs, hair.strands[strand_id].color);
        }
    }
}

inline void load_cyhair_shape(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<vec4f>& color, vector<float>& radius,
    bool flip_texcoord) {
    // load hair file
    auto hair = cyhair_data();
    load_cyhair(filename, hair);

    // generate curve data
    for (auto& strand : hair.strands) {
        auto offset = (int)positions.size();
        for (auto segment = 0; segment < (int)strand.positions.size() - 1;
             segment++) {
            lines.push_back({offset + segment, offset + segment + 1});
        }
        positions.insert(
            positions.end(), strand.positions.begin(), strand.positions.end());
        if (strand.radius.empty()) {
            radius.insert(
                radius.end(), strand.positions.size(), hair.default_thickness);
        } else {
            radius.insert(
                radius.end(), strand.radius.begin(), strand.radius.end());
        }
        if (strand.color.empty()) {
            color.insert(color.end(), strand.positions.size(),
                {hair.default_color.x, hair.default_color.y,
                    hair.default_color.z, 1});
        } else {
            for (auto i = 0; i < strand.color.size(); i++) {
                auto scolor = strand.color[i];
                color.push_back({scolor.x, scolor.y, scolor.z, 1});
            }
        }
    }

    // flip yz
    for (auto& p : positions) std::swap(p.y, p.z);

    // compute tangents
    normals.resize(positions.size());
    compute_vertex_tangents(normals, lines, positions);

    // fix colors
    for (auto& c : color) c = srgb_to_linear(c);
}

}  // namespace yocto

#endif
