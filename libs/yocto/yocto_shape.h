//
// # Yocto/Shape: Tiny Library for shape operations for graphics
//
//
// Yocto/Shape is a collection of utilities for manipulating shapes in 3D
// graphics, with a focus on triangle and quad meshes.
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
// 1. compute smooth normals and tangents with `compute_normals()`
//   `compute_tangents()`
// 2. compute tangent frames from texture coordinates with
//    `compute_tangent_spaces()`
// 3. compute skinning with `compute_skinning()` and
//    `compute_matrix_skinning()`
// 4. create shapes with `make_proc_image()`, `make_hair()`,
//    `make_points()`, `make_point()`
// 5. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
// 6. shape sampling with `sample_points()`, `sample_lines()`,
//    `sample_triangles()`; initialize the sampling CDFs with
//    `sample_points_cdf()`, `sample_lines_cdf()`,
//    `sample_triangles_cdf()`
// 7.  sample a could of point over a surface with `sample_triangles()`
// 8. get edges and boundaries with `get_edges()`
// 9. convert quads to triangles with `quads_to_triangles()`
// 10. convert face varying to vertex shared representations with
//     `convert_face_varying()`
// 11. subdivide elements by edge splits with `subdivide_lines()`,
//     `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
// 12. Catmull-Clark subdivision surface with `subdivide_catmullclark()`
//
//
// ## Shape IO
//
// We support reading and writing shapes in OBJ and PLY.
//
// 1. load/save shapes with `load_shape()`/`save_shape()`
//
//
// ## Ray-Scene and Closest-Point Queries
//
// Yocto/BVH is a simple implementation of ray intersection and
// closest queries using a two-level BVH data structure. We also include
// low-level intersection and closet point primitives.
// Alternatively the library also support wrapping Intel's Embree.
//
// Yocto/GL provides ray-scene intersection for points, lines, triangles and
// quads accelerated by a two-level BVH data structure. Our BVH is written for
// minimal code and not maximum speed, but still gives reasonable results. We
// suggest the use of Intel's Embree as a more efficient alternative.
//
// In Yocto/Bvh, shapes are described as collections of indexed primitives
// (points/lines/triangles/quads) like the standard triangle mesh used in
// real-time graphics. A scene if represented as transformed instances of
// shapes. The internal data structure is a two-level BVH, with a BVH for each
// shape and one top-level BVH for the whole scene. This design support
// instancing for large scenes and easy BVH refitting for interactive
// applications.
//
// In these functions triangles are parameterized with uv written
// w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
// as pairs of two triangles p0,p1,p3 and p2,p3,p1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
// represent triangles correctly, an this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
// Shape and scene data is not copied from the application to the BVH to
// improve memory footprint at the price of convenience. Shape data is
// explixitly passed on evey call, while instance data uses callbacks,
// since each application has its own conventions for storing those.
// To make usage more convenite, we provide `bvh_XXX_data` to hold application
// data and convenience wrappers for all functions.
//
// We support working either on the whole scene or on a single shape. In the
// description below yoi will see this dual API defined.
//
// 1. build the shape/scene BVH with `make_XXX_bvh()`;
// 2. perform ray-shape intersection with `intersect_XXX_bvh()`
// 3. perform point overlap queries with `overlap_XXX_bvh()`
// 4. refit BVH for dynamic applications with `update_XXX_bvh`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <memory>
#include <tuple>
#include <unordered_map>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions);
vector<vec3f> compute_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
vector<vec3f> compute_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions);
// Update normals and tangents
void update_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions);
void update_normals(vector<vec3f>& normals, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
void update_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords);

// Apply skinning to vertex position and normals.
std::pair<vector<vec3f>, vector<vec3f>> compute_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
std::pair<vector<vec3f>, vector<vec3f>> compute_matrix_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms);
// Update skinning
void udpate_skinning(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms);
void update_matrix_skinning(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
vector<vec3f> flip_normals(const vector<vec3f>& normals);
// Flip face orientation
vector<vec3i> flip_triangles(const vector<vec3i>& triangles);
vector<vec4i> flip_quads(const vector<vec4i>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
vector<vec3f> align_vertices(
    const vector<vec3f>& positions, const vec3i& alignment);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EDGES AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Dictionary to store edge information. `index` is the index to the edge
// array, `edges` the array of edges and `nfaces` the number of adjacent faces.
// We store only bidirectional edges to keep the dictionary small. Use the
// functions below to access this data.
struct edge_map {
  unordered_map<vec2i, int> index  = {};
  vector<vec2i>             edges  = {};
  vector<int>               nfaces = {};
};

// Initialize an edge map with elements.
edge_map make_edge_map(const vector<vec3i>& triangles);
edge_map make_edge_map(const vector<vec4i>& quads);
void     insert_edges(edge_map& emap, const vector<vec3i>& triangles);
void     insert_edges(edge_map& emap, const vector<vec4i>& quads);
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge);
// Get the edge index
int edge_index(const edge_map& emap, const vec2i& edge);
// Get edges and boundaries
int           num_edges(const edge_map& emap);
vector<vec2i> get_edges(const edge_map& emap);
vector<vec2i> get_boundary(const edge_map& emap);
vector<vec2i> get_edges(const vector<vec3i>& triangles);
vector<vec2i> get_edges(const vector<vec4i>& quads);

// Build adjacencies between faces (sorted counter-clockwise)
vector<vec3i> face_adjacencies(const vector<vec3i>& triangles);

// Build adjacencies between vertices (sorted counter-clockwise)
vector<vector<int>> vertex_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies);

// Compute boundaries as a list of loops (sorted counter-clockwise)
vector<vector<int>> ordered_boundaries(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacency, int num_vertices);

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
vector<vector<int>> vertex_to_faces_adjacencies(
    const vector<vec3i>& triangles, const vector<vec3i>& adjacencies);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BVH, RAY INTERSECTION AND OVERLAP QUERIES
// -----------------------------------------------------------------------------
namespace yocto {

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal nodes.
struct bvh_node {
  bbox3f bbox;
  int    start;
  short  num;
  bool   internal;
  byte   axis;
};

// BVH tree stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
// Application data is not stored explicitly.
struct bvh_tree {
  vector<bvh_node> nodes      = {};
  vector<int>      primitives = {};
};

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_intersection {
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Make shape bvh
void make_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
void make_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
void make_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius);
void make_quads_bvh(bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius);

// Updates shape bvh for changes in positions and radia
void update_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius);
void update_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
void update_quads_bvh(
    bvh_tree& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions);

// Find a shape element or scene instances that intersects a ray,
// returning either the closest or any overlap depending on `find_any`.
// Returns the point distance, the instance id, the shape element index and
// the element barycentric coordinates.
bvh_intersection intersect_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);
bvh_intersection intersect_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const ray3f& ray, bool find_any = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bvh_intersection overlap_points_bvh(const bvh_tree& bvh,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_lines_bvh(const bvh_tree& bvh,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_quads_bvh(const bvh_tree& bvh,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------
namespace yocto {

// A sparse grid of cells, containing list of points. Cells are stored in
// a dictionary to get sparsity. Helpful for nearest neighboor lookups.
struct hash_grid {
  float                             cell_size     = 0;
  float                             cell_inv_size = 0;
  vector<vec3f>                     positions     = {};
  unordered_map<vec3i, vector<int>> cells         = {};
};

// Create a hash_grid
hash_grid make_hash_grid(float cell_size);
hash_grid make_hash_grid(const vector<vec3f>& positions, float cell_size);
// Inserts a point into the grid
int insert_vertex(hash_grid& grid, const vec3f& position);
// Finds the nearest neighbors within a given radius
void find_neighbors(const hash_grid& grid, vector<int>& neighbors,
    const vec3f& position, float max_radius);
void find_neighbors(const hash_grid& grid, vector<int>& neighbors, int vertex,
    float max_radius);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yocto {

// Convert quads to triangles
vector<vec3i> quads_to_triangles(const vector<vec4i>& quads);
// Convert triangles to quads by creating degenerate quads
vector<vec4i> triangles_to_quads(const vector<vec3i>& triangles);
// Convert beziers to lines using 3 lines for each bezier.
vector<vec4i> bezier_to_lines(vector<vec2i>& lines);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm, texcoord and colors.
std::tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
split_facevarying(const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords);

// Split primitives per id
vector<vector<vec2i>> ungroup_lines(
    const vector<vec2i>& lines, const vector<int>& ids);
vector<vector<vec3i>> ungroup_triangles(
    const vector<vec3i>& triangles, const vector<int>& ids);
vector<vector<vec4i>> ungroup_quads(
    const vector<vec4i>& quads, const vector<int>& ids);

// Weld vertices within a threshold.
std::pair<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold);
std::pair<vector<vec3i>, vector<vec3f>> weld_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    float threshold);
std::pair<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold);

// Merge shape elements
void merge_lines(
    vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts);
void merge_triangles(vector<vec3i>& triangles,
    const vector<vec2i>& merge_triangles, int num_verts);
void merge_quads(
    vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts);
void merge_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& tangents, vector<vec2f>& texcoords, vector<float>& radius,
    const vector<vec2i>& merge_lines, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_tangents,
    const vector<vec2f>& merge_texturecoords,
    const vector<float>& merge_radius);
void merge_triangles(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vector<vec2i>& merge_triangles, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals,
    const vector<vec2f>& merge_texturecoords);
void merge_quads(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vector<vec4i>& merge_quads, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals,
    const vector<vec2f>& merge_texturecoords);

// Merge quads and triangles
void merge_triangles_and_quads(
    vector<vec3i>& triangles, vector<vec4i>& quads, bool force_triangles);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines by splitting each line in half.
std::pair<vector<vec2i>, vector<float>> subdivide_lines(
    const vector<vec2i>& lines, const vector<float>& vert, int level);
std::pair<vector<vec2i>, vector<vec2f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec2f>& vert, int level);
std::pair<vector<vec2i>, vector<vec3f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec3f>& vert, int level);
std::pair<vector<vec2i>, vector<vec4f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec4f>& vert, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
std::pair<vector<vec3i>, vector<float>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<float>& vert, int level);
std::pair<vector<vec3i>, vector<vec2f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec2f>& vert, int level);
std::pair<vector<vec3i>, vector<vec3f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& vert, int level);
std::pair<vector<vec3i>, vector<vec4f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec4f>& vert, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
std::pair<vector<vec4i>, vector<float>> subdivide_quads(
    const vector<vec4i>& quads, const vector<float>& vert, int level);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec2f>& vert, int level);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec3f>& vert, int level);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec4f>& vert, int level);
// Subdivide beziers by splitting each segment in two.
std::pair<vector<vec4i>, vector<float>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<float>& vert, int level);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec2f>& vert, int level);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec3f>& vert, int level);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec4f>& vert, int level);
// Subdivide quads using Carmull-Clark subdivision rules.
std::pair<vector<vec4i>, vector<float>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<float>& vert, int level,
    bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec2f>& vert, int level,
    bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec3f>& vert, int level,
    bool lock_boundary = false);
std::pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec4f>& vert, int level,
    bool lock_boundary = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
int           sample_points(int npoints, float re);
int           sample_points(const vector<float>& cdf, float re);
vector<float> sample_points_cdf(int npoints);

// Pick a point on lines uniformly.
std::pair<int, float> sample_lines(
    const vector<float>& cdf, float re, float ru);
vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
std::pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
std::pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv);
std::pair<int, vec2f> sample_quads(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float>         sample_quads_cdf(
            const vector<vec4i>& quads, const vector<vec3f>& positions);

// Samples a set of points over a triangle/quad mesh uniformly. Returns pos,
// norm and texcoord of the sampled points.
void sample_triangles(vector<vec3f>& sampled_positions,
    vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
    int seed = 7);
void sample_quads(vector<vec3f>& sampled_positions,
    vector<vec3f>& sampled_normals, vector<vec2f>& sampled_texcoords,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords, int npoints,
    int seed = 7);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape as indexed meshes
[[nodiscard]] bool load_shape(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, vector<float>& radius, string& error,
    bool flip_texcoords = true);
[[nodiscard]] bool save_shape(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, const vector<float>& radius, string& error,
    bool ascii = false, bool flip_texcoords = true);

// Load/save a facevarying shape
[[nodiscard]] bool load_fvshape(const string& filename, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_fvshape(const string& filename,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords, string& error,
    bool ascii = false, bool flip_texcoords = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Get mesh statistics for printing
vector<string> shape_stats(const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<vec3f>& colors,
    const vector<float>& radius, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a plane.
void make_rect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
void make_bulged_rect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1}, float radius = 0.3);
void make_recty(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
// Make a box.
void make_box(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
void make_rounded_box(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1}, float radius = 0.3);
// Make a quad stack
void make_rect_stack(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
// Make a floor.
void make_floor(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
    const vec2f& uvscale = {10, 10});
void make_bent_floor(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
    const vec2f& uvscale = {10, 10}, float bent = 0.5);
// Make a sphere.
void make_sphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);
// Make a sphere.
void make_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a sphere with slipped caps.
void make_capped_uvsphere(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1}, float height = 0.3);
// Make a disk
void make_disk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);
// Make a bulged disk
void make_bulged_disk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1, float height = 0.3);
// Make a uv disk
void make_uvdisk(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv cylinder
void make_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1},
    const vec3f& uvscale = {1, 1, 1});
// Make a rounded uv cylinder
void make_rounded_uvcylinder(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1},
    const vec3f& uvscale = {1, 1, 1}, float radius = 0.3);
// Make a plane in the xz plane.
void make_yrect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
void make_bulged_yrect(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1}, float radius = 0.3);

// Make a facevarying rect
void make_fvrect(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
// Make a facevarying box
void make_fvbox(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
void make_fvsphere(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
void make_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    int num = 65536, const vec2i& steps = {4, 65536},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1},
    const vec2f& rad = {0.001, 0.001});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
void make_point(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    float point_radius);
void make_points(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    int num = 65536, float uvscale = 1, float point_radius = 0.001);
void make_random_points(vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    int num = 65536, const vec3f& size = {1, 1, 1}, float uvscale = 1,
    float point_radius = 0.001, uint64_t seed = 17);

// Predefined meshes
void make_monkey(
    vector<vec4i>& quads, vector<vec3f>& positions, float scale = 1);
void make_quad(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_quady(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_cube(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_fvcube(vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float scale = 1);
void make_geosphere(
    vector<vec3i>& triangles, vector<vec3f>& positions, float scale = 1);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (strength/number)
// rotation: rotation added to hair (angle/strength)
void make_hair(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<float>& radius,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2i& steps = {8, 65536},
    const vec2f& length = {0.1, 0.1}, const vec2f& rad = {0.001, 0.001},
    const vec2f& noise = {0, 10}, const vec2f& clump = {0, 128},
    const vec2f& rotation = {0, 0}, int seed = 7);

// Thickens a shape by copying the shape content, rescaling it and flipping its
// normals. Note that this is very much not robust and only useful for trivial
// cases.
void make_shell(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, float thickness);

// Make a heightfield mesh.
void make_heightfield(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& size,
    const vector<float>& height);

}  // namespace yocto

#endif
