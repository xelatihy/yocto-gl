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
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Math defitions
using math::bbox3f;
using math::byte;
using math::flt_max;
using math::frame3f;
using math::identity3x4f;
using math::mat3f;
using math::mat4f;
using math::ray3f;
using math::vec2f;
using math::vec2i;
using math::vec3b;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Compute per-vertex normals/tangents for lines/triangles/quads.
std::vector<vec3f> compute_tangents(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& positions);
std::vector<vec3f> compute_normals(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions);
std::vector<vec3f> compute_normals(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions);
// Update normals and tangents
void update_tangents(std::vector<vec3f>& tangents,
    const std::vector<vec2i>& lines, const std::vector<vec3f>& positions);
void update_normals(std::vector<vec3f>& normals,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions);
void update_normals(std::vector<vec3f>& normals,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component std::vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
std::vector<vec4f> compute_tangent_spaces(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords);

// Apply skinning to vertex position and normals.
std::pair<std::vector<vec3f>, std::vector<vec3f>> compute_skinning(
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec4f>& weights, const std::vector<vec4i>& joints,
    const std::vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
std::pair<std::vector<vec3f>, std::vector<vec3f>> compute_matrix_skinning(
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec4f>& weights, const std::vector<vec4i>& joints,
    const std::vector<mat4f>& xforms);
// Update skinning
void udpate_skinning(std::vector<vec3f>& skinned_positions,
    std::vector<vec3f>& skinned_normals, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms);
void update_matrix_skinning(std::vector<vec3f>& skinned_positions,
    std::vector<vec3f>& skinned_normals, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Flip vertex normals
std::vector<vec3f> flip_normals(const std::vector<vec3f>& normals);
// Flip face orientation
std::vector<vec3i> flip_triangles(const std::vector<vec3i>& triangles);
std::vector<vec4i> flip_quads(const std::vector<vec4i>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
std::vector<vec3f> align_vertices(
    const std::vector<vec3f>& positions, const vec3i& alignment);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// EDGEA AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Dictionary to store edge information. `index` is the index to the edge
// array, `edges` the array of edges and `nfaces` the number of adjacent faces.
// We store only bidirectional edges to keep the dictionary small. Use the
// functions below to access this data.
struct edge_map {
  std::unordered_map<vec2i, int> index  = {};
  std::vector<vec2i>             edges  = {};
  std::vector<int>               nfaces = {};
};

// Initialize an edge map with elements.
edge_map make_edge_map(const std::vector<vec3i>& triangles);
edge_map make_edge_map(const std::vector<vec4i>& quads);
void     insert_edges(edge_map& emap, const std::vector<vec3i>& triangles);
void     insert_edges(edge_map& emap, const std::vector<vec4i>& quads);
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge);
// Get the edge index
int edge_index(const edge_map& emap, const vec2i& edge);
// Get edges and boundaries
int                num_edges(const edge_map& emap);
std::vector<vec2i> get_edges(const edge_map& emap);
std::vector<vec2i> get_boundary(const edge_map& emap);
std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles);
std::vector<vec2i> get_edges(const std::vector<vec4i>& quads);

// Build adjacencies between faces (sorted counter-clockwise)
std::vector<vec3i> face_adjacencies(const std::vector<vec3i>& triangles);

// Build adjacencies between vertices (sorted counter-clockwise)
std::vector<std::vector<int>> vertex_adjacencies(
    const std::vector<vec3i>& triangles, const std::vector<vec3i>& adjacencies);

// Compute boundaries as a list of loops (sorted counter-clockwise)
std::vector<std::vector<int>> ordered_boundaries(
    const std::vector<vec3i>& triangles, const std::vector<vec3i>& adjacency,
    int num_vertices);

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
std::vector<std::vector<int>> vertex_to_faces_adjacencies(
    const std::vector<vec3i>& triangles, const std::vector<vec3i>& adjacencies);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// BVH, RAY INTERSECTION AND OVERLAP QUERIES
// -----------------------------------------------------------------------------
namespace yocto::shape {

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
  std::vector<bvh_node> nodes      = {};
  std::vector<int>      primitives = {};
};

// Results of intersect_xxx and overlap_xxx functions that include hit flag,
// instance id, shape element id, shape element uv and intersection distance.
// The values are all set for scene intersection. Shape intersection does not
// set the instance id and element intersections do not set shape element id
// and the instance id. Results values are set only if hit is true.
struct bvh_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
};

// Make shape bvh
void make_points_bvh(bvh_tree& bvh, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);
void make_lines_bvh(bvh_tree& bvh, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);
void make_triangles_bvh(bvh_tree& bvh, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);
void make_quads_bvh(bvh_tree& bvh, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);

// Updates shape bvh for changes in positions and radia
void update_points_bvh(bvh_tree& bvh, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);
void update_lines_bvh(bvh_tree& bvh, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<float>& radius);
void update_triangles_bvh(bvh_tree& bvh, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions);
void update_quads_bvh(bvh_tree& bvh, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions);

// Find a shape element or scene instances that intersects a ray,
// returning either the closest or any overlap depending on `find_any`.
// Returns the point distance, the instance id, the shape element index and
// the element barycentric coordinates.
bvh_intersection intersect_points_bvh(const bvh_tree& bvh,
    const std::vector<int>& points, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_lines_bvh(const bvh_tree& bvh,
    const std::vector<vec2i>& lines, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_triangles_bvh(const bvh_tree& bvh,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    const ray3f& ray, bool find_any = false);
bvh_intersection intersect_quads_bvh(const bvh_tree& bvh,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions,
    const ray3f& ray, bool find_any = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bvh_intersection overlap_points_bvh(const bvh_tree& bvh,
    const std::vector<int>& points, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_lines_bvh(const bvh_tree& bvh,
    const std::vector<vec2i>& lines, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_triangles_bvh(const bvh_tree& bvh,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);
bvh_intersection overlap_quads_bvh(const bvh_tree& bvh,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions,
    const std::vector<float>& radius, const vec3f& pos, float max_distance,
    bool find_any = false);

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_shape {
  // elements
  std::vector<int>   points    = {};
  std::vector<vec2i> lines     = {};
  std::vector<vec3i> triangles = {};
  std::vector<vec4i> quads     = {};

  // vertices
  std::vector<vec3f> positions = {};
  std::vector<float> radius    = {};

  // nodes
  bvh_tree bvh = {};
#ifdef YOCTO_EMBREE
  std::shared_ptr<void> embree_bvh = {};
#endif
};

// instance
struct bvh_instance {
  frame3f frame = identity3x4f;
  int     shape = -1;
};

// BVH data for whole shapes. This interface makes copies of all the data.
struct bvh_scene {
  // instances and shapes
  std::vector<bvh_instance> instances = {};
  std::vector<bvh_shape>    shapes    = {};

  // nodes
  bvh_tree bvh = {};
#ifdef YOCTO_EMBREE
  std::shared_ptr<void> embree_bvh = {};
#endif
};

// Build the bvh acceleration structure.
void init_shape_bvh(bvh_shape& bvh, bool embree = false);
void init_scene_bvh(bvh_scene& bvh, bool embree = false);

// Refit bvh data
void update_shape_bvh(bvh_shape& bvh);
void update_scene_bvh(bvh_scene& bvh, const std::vector<int>& updated_instances,
    const std::vector<int>& updated_shapes);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bvh_intersection intersect_shape_bvh(
    const bvh_shape& bvh, const ray3f& ray, bool find_any = false);
bvh_intersection intersect_scene_bvh(const bvh_scene& bvh, const ray3f& ray,
    bool find_any = false, bool non_rigid_frames = true);
bvh_intersection intersect_instance_bvh(const bvh_scene& bvh, int instance,
    const ray3f& ray, bool find_any = false, bool non_rigid_frames = true);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bvh_intersection overlap_shape_bvh(const bvh_shape& bvh, const vec3f& pos,
    float max_distance, bool find_any = false);
bvh_intersection overlap_scene_bvh(const bvh_scene& bvh, const vec3f& pos,
    float max_distance, bool find_any = false, bool non_rigid_frames = true);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------
namespace yocto::shape {

// A sparse grid of cells, containing list of points. Cells are stored in
// a dictionary to get sparsity. Helpful for nearest neighboor lookups.
struct hash_grid {
  float                                       cell_size     = 0;
  float                                       cell_inv_size = 0;
  std::vector<vec3f>                          positions     = {};
  std::unordered_map<vec3i, std::vector<int>> cells         = {};
};

// Create a hash_grid
hash_grid make_hash_grid(float cell_size);
hash_grid make_hash_grid(const std::vector<vec3f>& positions, float cell_size);
// Inserts a point into the grid
int insert_vertex(hash_grid& grid, const vec3f& position);
// Finds the nearest neighbors within a given radius
void find_neighbors(const hash_grid& grid, std::vector<int>& neighbors,
    const vec3f& position, float max_radius);
void find_neighbors(const hash_grid& grid, std::vector<int>& neighbors,
    int vertex, float max_radius);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Convert quads to triangles
std::vector<vec3i> quads_to_triangles(const std::vector<vec4i>& quads);
// Convert triangles to quads by creating degenerate quads
std::vector<vec4i> triangles_to_quads(const std::vector<vec3i>& triangles);
// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec4i> bezier_to_lines(std::vector<vec2i>& lines);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm, texcoord and colors.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
split_facevarying(const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>&               quadsnorm,
    const std::vector<vec4i>&               quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords);

// Split primitives per id
std::vector<std::vector<vec2i>> ungroup_lines(
    const std::vector<vec2i>& lines, const std::vector<int>& ids);
std::vector<std::vector<vec3i>> ungroup_triangles(
    const std::vector<vec3i>& triangles, const std::vector<int>& ids);
std::vector<std::vector<vec4i>> ungroup_quads(
    const std::vector<vec4i>& quads, const std::vector<int>& ids);

// Weld vertices within a threshold.
std::pair<std::vector<vec3f>, std::vector<int>> weld_vertices(
    const std::vector<vec3f>& positions, float threshold);
std::pair<std::vector<vec3i>, std::vector<vec3f>> weld_triangles(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    float threshold);
std::pair<std::vector<vec4i>, std::vector<vec3f>> weld_quads(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions,
    float threshold);

// Merge shape elements
void merge_lines(std::vector<vec2i>& lines,
    const std::vector<vec2i>& merge_lines, int num_verts);
void merge_triangles(std::vector<vec3i>& triangles,
    const std::vector<vec2i>& merge_triangles, int num_verts);
void merge_quads(std::vector<vec4i>& quads,
    const std::vector<vec4i>& merge_quads, int num_verts);
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& positions,
    std::vector<vec3f>& tangents, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, const std::vector<vec2i>& merge_lines,
    const std::vector<vec3f>& merge_positions,
    const std::vector<vec3f>& merge_tangents,
    const std::vector<vec2f>& merge_texturecoords,
    const std::vector<float>& merge_radius);
void merge_triangles(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const std::vector<vec2i>& merge_triangles,
    const std::vector<vec3f>& merge_positions,
    const std::vector<vec3f>& merge_normals,
    const std::vector<vec2f>& merge_texturecoords);
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const std::vector<vec4i>& merge_quads,
    const std::vector<vec3f>& merge_positions,
    const std::vector<vec3f>& merge_normals,
    const std::vector<vec2f>& merge_texturecoords);

// Merge quads and triangles
void merge_triangles_and_quads(std::vector<vec3i>& triangles,
    std::vector<vec4i>& quads, bool force_triangles);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Subdivide lines by splitting each line in half.
std::pair<std::vector<vec2i>, std::vector<float>> subdivide_lines(
    const std::vector<vec2i>& lines, const std::vector<float>& vert, int level);
std::pair<std::vector<vec2i>, std::vector<vec2f>> subdivide_lines(
    const std::vector<vec2i>& lines, const std::vector<vec2f>& vert, int level);
std::pair<std::vector<vec2i>, std::vector<vec3f>> subdivide_lines(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& vert, int level);
std::pair<std::vector<vec2i>, std::vector<vec4f>> subdivide_lines(
    const std::vector<vec2i>& lines, const std::vector<vec4f>& vert, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
std::pair<std::vector<vec3i>, std::vector<float>> subdivide_triangles(
    const std::vector<vec3i>& triangles, const std::vector<float>& vert,
    int level);
std::pair<std::vector<vec3i>, std::vector<vec2f>> subdivide_triangles(
    const std::vector<vec3i>& triangles, const std::vector<vec2f>& vert,
    int level);
std::pair<std::vector<vec3i>, std::vector<vec3f>> subdivide_triangles(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& vert,
    int level);
std::pair<std::vector<vec3i>, std::vector<vec4f>> subdivide_triangles(
    const std::vector<vec3i>& triangles, const std::vector<vec4f>& vert,
    int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
std::pair<std::vector<vec4i>, std::vector<float>> subdivide_quads(
    const std::vector<vec4i>& quads, const std::vector<float>& vert, int level);
std::pair<std::vector<vec4i>, std::vector<vec2f>> subdivide_quads(
    const std::vector<vec4i>& quads, const std::vector<vec2f>& vert, int level);
std::pair<std::vector<vec4i>, std::vector<vec3f>> subdivide_quads(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& vert, int level);
std::pair<std::vector<vec4i>, std::vector<vec4f>> subdivide_quads(
    const std::vector<vec4i>& quads, const std::vector<vec4f>& vert, int level);
// Subdivide beziers by splitting each segment in two.
std::pair<std::vector<vec4i>, std::vector<float>> subdivide_beziers(
    const std::vector<vec4i>& beziers, const std::vector<float>& vert,
    int level);
std::pair<std::vector<vec4i>, std::vector<vec2f>> subdivide_beziers(
    const std::vector<vec4i>& beziers, const std::vector<vec2f>& vert,
    int level);
std::pair<std::vector<vec4i>, std::vector<vec3f>> subdivide_beziers(
    const std::vector<vec4i>& beziers, const std::vector<vec3f>& vert,
    int level);
std::pair<std::vector<vec4i>, std::vector<vec4f>> subdivide_beziers(
    const std::vector<vec4i>& beziers, const std::vector<vec4f>& vert,
    int level);
// Subdivide quads using Carmull-Clark subdivision rules.
std::pair<std::vector<vec4i>, std::vector<float>> subdivide_catmullclark(
    const std::vector<vec4i>& quads, const std::vector<float>& vert, int level,
    bool lock_boundary = false);
std::pair<std::vector<vec4i>, std::vector<vec2f>> subdivide_catmullclark(
    const std::vector<vec4i>& quads, const std::vector<vec2f>& vert, int level,
    bool lock_boundary = false);
std::pair<std::vector<vec4i>, std::vector<vec3f>> subdivide_catmullclark(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& vert, int level,
    bool lock_boundary = false);
std::pair<std::vector<vec4i>, std::vector<vec4f>> subdivide_catmullclark(
    const std::vector<vec4i>& quads, const std::vector<vec4f>& vert, int level,
    bool lock_boundary = false);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Pick a point in a point set uniformly.
int                sample_points(int npoints, float re);
int                sample_points(const std::vector<float>& cdf, float re);
std::vector<float> sample_points_cdf(int npoints);

// Pick a point on lines uniformly.
std::pair<int, float> sample_lines(
    const std::vector<float>& cdf, float re, float ru);
std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
std::pair<int, vec2f> sample_triangles(
    const std::vector<float>& cdf, float re, const vec2f& ruv);
std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
std::pair<int, vec2f> sample_quads(
    const std::vector<float>& cdf, float re, const vec2f& ruv);
std::pair<int, vec2f> sample_quads(const std::vector<vec4i>& quads,
    const std::vector<float>& cdf, float re, const vec2f& ruv);
std::vector<float>    sample_quads_cdf(
       const std::vector<vec4i>& quads, const std::vector<vec3f>& positions);

// Samples a set of points over a triangle/quad mesh uniformly. Returns pos,
// norm and texcoord of the sampled points.
void sample_triangles(std::vector<vec3f>& sampled_positions,
    std::vector<vec3f>& sampled_normals, std::vector<vec2f>& sampled_texcoords,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec2f>& texcoords,
    int npoints, int seed = 7);
void sample_quads(std::vector<vec3f>& sampled_positions,
    std::vector<vec3f>& sampled_normals, std::vector<vec2f>& sampled_texcoords,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec2f>& texcoords,
    int npoints, int seed = 7);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Data structure used for geodesic computation
struct geodesic_solver {
  static const int min_arcs = 12;
  struct graph_edge {
    int   node   = -1;
    float length = flt_max;
  };
#ifdef YOCTO_ABSEIL
  std::vector<short_vector<adjancency_list, min_arcs>> graph = {};
#else
  std::vector<std::vector<graph_edge>> graph = {};
#endif
};

// Construct a a graph to compute geodesic distances
geodesic_solver make_geodesic_solver(const std::vector<vec3i>& triangles,
    const std::vector<vec3i>& adjacencies, const std::vector<vec3f>& positions);

// Compute geodesic distances
void update_geodesic_distances(std::vector<float>& distances,
    const geodesic_solver& solver, const std::vector<int>& sources,
    float max_distance = flt_max);

std::vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const std::vector<int>& sources, float max_distance = flt_max);

// Compute all shortest paths from source vertices to any other vertex.
// Paths are implicitly represented: each node is assignes its previous node in
// the path. Graph search early exits when reching end_vertex.
std::vector<int> compute_geodesic_paths(const geodesic_solver& solver,
    const std::vector<int>& sources, int end_vertex = -1);

// Sample vertices with a Poisson distribution using geodesic distances.
// Sampling strategy is farthest point sampling (FPS): at every step
// take the farthers point from current sampled set until done.
std::vector<int> sample_vertices_poisson(
    const geodesic_solver& solver, int num_samples);

// Compute the distance field needed to compute a voronoi diagram
std::vector<std::vector<float>> compute_voronoi_fields(
    const geodesic_solver& solver, const std::vector<int>& generators);

// Convert distances to colors
std::vector<vec3f> colors_from_field(const std::vector<float>& field,
    float scale = 1, const vec3f& c0 = {1, 1, 1},
    const vec3f& c1 = {1, 0.1, 0.1});

// Description of a discrete path along the surface of a triangle mesh.
struct surface_path {
  struct vertex {
    vec2i edge  = {0, 0};
    int   face  = 0;
    float alpha = 0;
  };
  int                 start, end;
  std::vector<vertex> vertices;
};

// Trace integral path following the gradient of a scalar field
surface_path integrate_field(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<vec3i>& adjacency,
    const std::vector<int>& tags, int tag, const std::vector<float>& field,
    int from);
surface_path integrate_field(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<vec3i>& adjacency,
    const std::vector<int>& tags, int tag, const std::vector<float>& field,
    int from, int to);

std::vector<vec3f> make_positions_from_path(
    const surface_path& path, const std::vector<vec3f>& mesh_positions);

vec3f compute_gradient(const vec3i& triangle,
    const std::vector<vec3f>& positions, const std::vector<float>& field);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Load/save a shape as indexed meshes
[[nodiscard]] bool load_shape(const std::string& filename,
    std::vector<int>& points, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<vec3f>& colors,
    std::vector<float>& radius, std::string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_shape(const std::string& filename,
    const std::vector<int>& points, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<vec3f>& colors,
    const std::vector<float>& radius, std::string& error, bool ascii = false,
    bool flip_texcoords = true);

// Load/save a facevarying shape
[[nodiscard]] bool load_fvshape(const std::string& filename,
    std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_fvshape(const std::string& filename,
    const std::vector<vec4i>& quadspos, const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, std::string& error, bool ascii = false,
    bool flip_texcoords = true);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Get mesh statistics for printing
std::vector<std::string> shape_stats(const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<vec3f>& colors,
    const std::vector<float>& radius, bool verbose = false);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Make a plane.
void make_rect(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
void make_bulged_rect(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1}, float radius = 0.3);
// Make a box.
void make_box(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
void make_rounded_box(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1}, float radius = 0.3);
// Make a quad stack
void make_rect_stack(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
// Make a floor.
void make_floor(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
    const vec2f& uvscale = {10, 10});
void make_bent_floor(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {10, 10},
    const vec2f& uvscale = {10, 10}, float bent = 0.5);
// Make a sphere.
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);
// Make a sphere.
void make_uvsphere(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a sphere with slipped caps.
void make_capped_uvsphere(std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec2i& steps = {32, 32},
    float scale = 1, const vec2f& uvscale = {1, 1}, float height = 0.3);
// Make a disk
void make_disk(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);
// Make a bulged disk
void make_bulged_disk(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1, float height = 0.3);
// Make a uv disk
void make_uvdisk(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {32, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv cylinder
void make_uvcylinder(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1},
    const vec3f& uvscale = {1, 1, 1});
// Make a rounded uv cylinder
void make_rounded_uvcylinder(std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1},
    float radius = 0.3);
// Make a plane in the xz plane.
void make_yrect(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
void make_bulged_yrect(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1}, float radius = 0.3);

// Make a facevarying rect
void make_fvrect(std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
// Make a facevarying box
void make_fvbox(std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1},
    const vec3f& uvscale = {1, 1, 1});
void make_fvsphere(std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, int steps = 32,
    float scale = 1, float uvscale = 1);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, int num = 65536,
    const vec2i& steps = {4, 65536}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1}, const vec2f& rad = {0.001, 0.001});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
void make_point(std::vector<int>& points, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, float point_radius);
void make_points(std::vector<int>& points, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, int num, float uvscale, float point_radius);
void make_random_points(std::vector<int>& points, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, int num, const vec3f& size, float uvscale,
    float point_radius, uint64_t seed);

// Predefined meshes
void make_monkey(
    std::vector<vec4i>& quads, std::vector<vec3f>& positions, float scale = 1);
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    float scale = 1);
void make_quady(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    float scale = 1);
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    float scale = 1);
void make_fvcube(std::vector<vec4i>& quadspos, std::vector<vec4i>& quadsnorm,
    std::vector<vec4i>& quadstexcoord, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    float scale = 1);
void make_geosphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, float scale = 1);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (strength/number)
// rotation: rotation added to hair (angle/strength)
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<float>& radius, const std::vector<vec3i>& striangles,
    const std::vector<vec4i>& squads, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2i& steps = {8, 65536}, const vec2f& length = {0.1, 0.1},
    const vec2f& rad = {0.001, 0.001}, const vec2f& noise = {0, 10},
    const vec2f& clump = {0, 128}, const vec2f& rotation = {0, 0},
    int seed = 7);

// Thickens a shape by copying the shape content, rescaling it and flipping its
// normals. Note that this is very much not robust and only useful for trivial
// cases.
void make_shell(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    float thickness);

// Make a heightfield mesh.
void make_heightfield(std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    const vec2i& size, const std::vector<float>& height);

}  // namespace yocto::shape

// -----------------------------------------------------------------------------
// PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto::shape {

// Extract isoline from surface scalar field.
void meandering_triangles(const std::vector<float>& field, float isoline,
    int selected_tag, int t0, int t1, std::vector<vec3i>& triangles,
    std::vector<int>& tags, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals);

}  // namespace yocto::shape

#endif
