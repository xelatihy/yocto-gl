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
// `make_points()`
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <tuple>
#include <unordered_map>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yshp {

// Math defitions
using ym::bbox3f;
using ym::byte;
using ym::frame3f;
using ym::mat4f;
using ym::vec2f;
using ym::vec2i;
using ym::vec3b;
using ym::vec3f;
using ym::vec3i;
using ym::vec4f;
using ym::vec4i;

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yshp {

// Flip vertex normals
std::vector<vec3f> flip_normals(const std::vector<vec3f>& normals);
// Flip face orientation
std::vector<vec3i> flip_triangles(const std::vector<vec3i>& triangles);
std::vector<vec4i> flip_quads(const std::vector<vec4i>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
std::vector<vec3f> align_vertices(
    const std::vector<vec3f>& positions, const vec3i& alignment);

}  // namespace yshp

// -----------------------------------------------------------------------------
// EDGEA AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// HASH GRID AND NEAREST NEIGHBORS
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE ELEMENT CONVERSION AND GROUPING
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yshp {

// Data structure used for geodesic computation
struct geodesic_solver {
  static const int min_arcs = 12;
  struct graph_edge {
    int   node   = -1;
    float length = ym::flt_max;
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
    float max_distance = ym::flt_max);

std::vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const std::vector<int>& sources, float max_distance = ym::flt_max);

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yshp {

// Get mesh statistics for printing
std::vector<std::string> shape_stats(const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<vec3f>& colors,
    const std::vector<float>& radius, bool verbose = false);

}  // namespace yshp

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yshp {

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

}  // namespace yshp

// -----------------------------------------------------------------------------
// PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yshp {

// Extract isoline from surface scalar field.
void meandering_triangles(const std::vector<float>& field, float isoline,
    int selected_tag, int t0, int t1, std::vector<vec3i>& triangles,
    std::vector<int>& tags, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals);

}  // namespace yshp

#endif
