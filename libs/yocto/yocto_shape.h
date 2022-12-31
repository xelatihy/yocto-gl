//
// # Yocto/Shape: Shape utilities
//
// Yocto/Shape is a collection of utilities for manipulating shapes in 3D
// graphics, with a focus on triangle and quad meshes.
// Yocto/Shape is implemented in `yocto_shape.h` and `yocto_shape.cpp`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "yocto_geometry.h"
#include "yocto_math.h"
#include "yocto_ndarray.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Shape data represented as indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct shape_data {
  // element data
  vector<int>   points    = {};
  vector<vec2i> lines     = {};
  vector<vec3i> triangles = {};
  vector<vec4i> quads     = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};
};

// Shape type
inline bool is_points(const shape_data& shape) { return !shape.points.empty(); }
inline bool is_lines(const shape_data& shape) { return !shape.lines.empty(); }
inline bool is_triangles(const shape_data& shape) {
  return !shape.triangles.empty();
}
inline bool is_quads(const shape_data& shape) { return !shape.quads.empty(); }

// Interpolate vertex data
vec3f eval_position(const shape_data& shape, int element, const vec2f& uv);
vec3f eval_normal(const shape_data& shape, int element, const vec2f& uv);
vec3f eval_tangent(const shape_data& shape, int element, const vec2f& uv);
vec2f eval_texcoord(const shape_data& shape, int element, const vec2f& uv);
vec4f eval_color(const shape_data& shape, int element, const vec2f& uv);
float eval_radius(const shape_data& shape, int element, const vec2f& uv);

// Evaluate element normals
vec3f eval_element_normal(const shape_data& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const shape_data& shape);

// An unevaluated location on a shape
struct shape_point {
  int   element = 0;
  vec2f uv      = {0, 0};
};

// Shape sampling
vector<float> sample_shape_cdf(const shape_data& shape);
shape_point   sample_shape(const shape_data& shape, const vector<float>& cdf,
      float rn, const vec2f& ruv);
vector<shape_point> sample_shape(
    const shape_data& shape, int num_samples, uint64_t seed = 98729387);

// Conversions
shape_data quads_to_triangles(const shape_data& shape);
void       quads_to_triangles_inplace(shape_data& shape);

// Subdivision
shape_data subdivide_shape(
    const shape_data& shape, int subdivisions, bool catmullclark);

// Transform shape
shape_data transform_shape(
    const shape_data& shape, const frame3f& frame, bool non_rigid = false);
shape_data transform_shape(const shape_data& shape, const frame3f& frame,
    float radius_scale, bool non_rigid = false);
shape_data remove_normals(const shape_data& shape);
shape_data add_normals(const shape_data& shape);

// Merge a shape into another
void merge_shape_inplace(shape_data& shape, const shape_data& merge);

// Shape statistics
vector<string> shape_stats(const shape_data& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FACE-VARYING SHAPE DATA AND UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Shape data stored as a face-varying mesh
struct fvshape_data {
  // element data
  vector<vec4i> quadspos      = {};
  vector<vec4i> quadsnorm     = {};
  vector<vec4i> quadstexcoord = {};

  // vertex data
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
};

// Interpolate vertex data
vec3f eval_position(const fvshape_data& shape, int element, const vec2f& uv);
vec3f eval_normal(const fvshape_data& shape, int element, const vec2f& uv);
vec2f eval_texcoord(const shape_data& shape, int element, const vec2f& uv);

// Evaluate element normals
vec3f eval_element_normal(const fvshape_data& shape, int element);

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_normals(const fvshape_data& shape);
void compute_normals(vector<vec3f>& normals, const fvshape_data& shape);

// Conversions
shape_data fvshape_to_shape(
    const fvshape_data& shape, bool as_triangles = false);
fvshape_data shape_to_fvshape(const shape_data& shape);

// Subdivision
fvshape_data subdivide_fvshape(
    const fvshape_data& shape, int subdivisions, bool catmullclark);

// Transform shape
fvshape_data transform_fvshape(
    const frame3f& frame, const fvshape_data& shape, bool non_rigid = false);

// Shape statistics
vector<string> fvshape_stats(const fvshape_data& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a quad
shape_data make_quad(int subdivisions = 0);
// Make a quad along y
shape_data make_quady(int subdivisions = 0);
// Make a cube
shape_data make_cube(int subdivisions = 0);
// Make a sphere
shape_data make_sphere(int subdivisions = 5);
// Make a geosphere
shape_data make_geosphere(int subdivisions = 3);
// Make a cube that is watertight (only positions)
shape_data make_wtcube(int subdivisions = 0);
// Make a sphere that is watertight (only positions and normals)
shape_data make_wtsphere(int subdivisions = 5);
// Make a disk
shape_data make_disk(int subdivisions = 5);
// Make a floor
shape_data make_floor(int subdivisions = 0, float size = 10);
// Make a monkey
shape_data make_monkey(int subdivisions = 0);
// Make a rounded cube
shape_data make_rounded_cube(int subdivisions = 5, float radius = 0.3f);
// Make a bulged quad
shape_data make_bulged_quad(int subdivisions = 5, float radius = 0.3f);
// Make a bulged disk
shape_data make_bulged_disk(int subdivisions = 5, float height = 0.3f);
// Make a bent disk
shape_data make_bent_floor(
    int subdivisions = 5, float size = 10, float bent = 0.5f);

// Make a rectangle
shape_data make_rect(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1});
// Make a rectangle along y
shape_data make_recty(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1});
// Make a box
shape_data make_box(
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1});
// Make a rectangle stack
shape_data make_rect_stack(int vsteps = 1, const vec2i& steps = {1, 1},
    const vec3f& scale = {1, 1, 1});
// Make a tesselated sphere
shape_data make_tsphere(const vec3i& steps = {32, 32, 32}, float scale = 1);
// Make a uv sphere
shape_data make_uvsphere(const vec2i& steps = {64, 32}, float scale = 1);
// Make a uv sphere aligned to y
shape_data make_uvspherey(const vec2i& steps = {64, 32}, float scale = 1);
// Make a uv disk
shape_data make_uvdisk(const vec2i& steps = {32, 1}, float scale = 1);
// Make a uv cylinder
shape_data make_uvcylinder(
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1});
// Make a uv capsule
shape_data make_uvcapsule(
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1});
// Make a uv cone
shape_data make_uvcone(
    const vec3i& steps = {32, 32, 32}, const vec2f& scale = {1, 1});
// Make a bulged rectangle
shape_data make_bulged_rect(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, float radius = 0.3f);
// Make a bulged rectangle along y
shape_data make_bulged_recty(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, float radius = 0.3f);
// Make a rounded box
shape_data make_rounded_box(const vec3i& steps = {32, 32, 32},
    const vec3f& scale = {1, 1, 1}, float radius = 0.3f);
// Make a capped uv sphere
shape_data make_capped_uvsphere(
    const vec2i& steps = {32, 32}, float scale = 1, float height = 0.3f);
// Make a capped uv spherey
shape_data make_capped_uvspherey(
    const vec2i& steps = {32, 32}, float scale = 1, float height = 0.3f);
// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, float radius = 0.3f);

// Make a face-varying cube
fvshape_data make_fvcube(int subdivisions = 0);

// Make a face-varying rectangle
fvshape_data make_fvrect(
    const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1});
// Make a face-varying box
fvshape_data make_fvbox(
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1});
// Make a face-varying sphere
fvshape_data make_fvsphere(int steps = 32, float scale = 1);

// Make a point.
shape_data make_point(float radius = 0.001f, bool generate_uv = true);
// Make many points at the origin (useful for particle systems).
shape_data make_points(
    int num = 65536, float radius = 0.001f, bool generate_uv = true);
// Make lines along a quad.
shape_data make_lines(int num = 65536, int steps = 4,
    const vec2f& scale = {1, 1}, const vec2f& radius = {0.001f, 0.001f});

// Make a grid of quads
shape_data make_quad_grid(const vec2i& steps = {256, 256}, float scale = 1);
// Make points on a grid.
shape_data make_point_grid(
    const vec2i& steps = {256, 256}, float scale = 1, float radius = 0.001f);
// Make lines on a grid.
shape_data make_line_grid(
    const vec2i& steps = {256, 256}, float scale = 1, float radius = 0.001f);

// Make random points in a cube.
// TODO: switch to quad?
shape_data make_random_points(int num = 65536, float radius = 0.001f,
    bool generate_uv = true, uint64_t seed = 17);
// Make points on a shape
shape_data make_random_points(const shape_data& shape, int num = 65536,
    float radius = 0.001f, bool generate_uv = true, int seed = 7);
// Make lines on a shape
shape_data make_random_lines(const shape_data& shape, int num = 65536,
    int steps = 8, const vec2f& length = {0.1f, 0.1f},
    const vec2f& radius = {0.001f, 0.001f}, bool generate_uv = true,
    int seed = 7);
// Make hairs on a shape
shape_data make_random_hairs(const shape_data& shape, int num = 65536,
    int steps = 8, const vec2f& length = {0.1f, 0.1f},
    const vec2f& radius = {0.001f, 0.001f}, float noise = 0,
    float gravity = 0.001f, bool generate_uv = true, int seed = 7);

// Convert points to small spheres and lines to small cylinders. This is
// intended for making very small primitives for display in interactive
// applications, so the spheres are low res.
shape_data points_to_spheres(
    const vector<vec3f>& vertices, int steps = 2, float scale = 0.01f);
shape_data polyline_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);
shape_data lines_to_cylinders(
    const vector<vec3f>& vertices, int steps = 4, float scale = 0.01f);
shape_data lines_to_cylinders(const vector<vec2i>& lines,
    const vector<vec3f>& positions, int steps = 4, float scale = 0.01f);

// Make a heightfield mesh.
shape_data make_heightfield(const array2d<float>& height);
shape_data make_heightfield(const array2d<vec4f>& height);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex normals/tangents for lines/triangles/quads.
template <typename T, typename I>
inline vector<vec<T, 3>> lines_tangents(
    const vector<vec<I, 2>>& lines, const vector<vec<T, 3>>& positions);
template <typename T, typename I>
inline vector<vec<T, 3>> triangles_normals(
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions);
template <typename T, typename I>
inline vector<vec<T, 3>> quads_normals(
    const vector<vec<I, 4>>& quads, const vector<vec<T, 3>>& positions);
// Update normals and tangents
template <typename T, typename I>
inline void lines_tangents(vector<vec<T, 3>>& tangents,
    const vector<vec<I, 2>>& lines, const vector<vec<T, 3>>& positions);
template <typename T, typename I>
inline void triangles_normals(vector<vec<T, 3>>& normals,
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions);
template <typename T, typename I>
inline void quads_normals(vector<vec<T, 3>>& normals,
    const vector<vec<I, 4>>& quads, const vector<vec<T, 3>>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
template <typename T, typename I>
inline vector<vec<T, 4>> triangle_tangent_spaces(
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 2>>& texcoords);

// Apply skinning to vertex position and normals.
template <typename T, typename I>
inline pair<vector<vec<T, 3>>, vector<vec<T, 3>>> skin_vertices(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 4>>& weights, const vector<vec<I, 4>>& joints,
    const vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
template <typename T, typename I>
inline pair<vector<vec<T, 3>>, vector<vec<T, 3>>> skin_matrices(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 4>>& weights, const vector<vec<I, 4>>& joints,
    const vector<mat4f>& xforms);
// Update skinning
template <typename T, typename I>
inline void skin_vertices(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec<I, 4>>& joints, const vector<frame3f>& xforms);
template <typename T, typename I>
inline void skin_matrices(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec<I, 4>>& joints, const vector<mat4f>& xforms);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
template <typename T>
inline vector<vec<T, 3>> flip_normals(const vector<vec<T, 3>>& normals);
// Flip face orientation
template <typename I>
inline vector<vec<I, 3>> flip_triangles(const vector<vec<I, 3>>& triangles);
template <typename I>
inline vector<vec<I, 4>> flip_quads(const vector<vec<I, 4>>& quads);
// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
template <typename T, typename I>
inline vector<vec<T, 3>> align_vertices(
    const vector<vec<T, 3>>& positions, const vec<I, 3>& alignment);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTOR HASHING
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <typename T, size_t N>
struct hash<yocto::vec<T, N>> {
  size_t operator()(const yocto::vec<T, N>& v) const {
    static const auto hasher = std::hash<T>();
    auto              h      = (size_t)0;
    for (auto i = 0; i < N; i++) {
      h ^= hasher(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    }
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// EDGES AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Data related to edges for an edge_map
template <typename I>
struct edge_map_edge {
  I index, nfaces;
};

// Dictionary to store edge information. `index` is the index to the edge
// array, `edges` the array of edges and `nfaces` the number of adjacent faces.
// We store only bidirectional edges to keep the dictionary small. Use the
// functions below to access this data.
template <typename I>
struct edge_map {
  unordered_map<vec<I, 2>, edge_map_edge<I>> edges = {};
};

// Initialize an edge map with elements.
template <typename I>
inline edge_map<I> make_edge_map(const vector<vec<I, 3>>& triangles);
template <typename I>
inline edge_map<I> make_edge_map(const vector<vec<I, 4>>& quads);
// Insert an edge and return its index
template <typename I>
inline I insert_edge(edge_map<I>& emap, const vec<I, 2>& edge);
// Get the edge index
template <typename I>
inline I edge_index(const edge_map<I>& emap, const vec<I, 2>& edge);
// Get edges and boundaries
template <typename I>
inline I num_edges(const edge_map<I>& emap);
template <typename I>
inline vector<vec<I, 2>> get_edges(const edge_map<I>& emap);
template <typename I>
inline vector<vec<I, 2>> get_boundary(const edge_map<I>& emap);

// Edges helpers
template <typename I>
constexpr inline vec<I, 2> _make_edge(const vec<I, 2>& hedge) {
  return {min(hedge), max(hedge)};
}
template <typename I>
constexpr inline vec<I, 2> _make_edge(I a, I b) {
  return {min(a, b), max(a, b)};
}

// Edge lists
template <typename I>
inline ptrdiff_t search_edge(
    const vector<vec<I, 2>>& edges, const vec<I, 2>& edge) {
  auto pos = std::binary_search(edges.begin(), edges.end(), _make_edge(edge));
  if (pos == edges.end()) return -1;
  return pos - edges.begin();
}
template <typename I>
inline vector<vec<I, 2>> get_edges(
    const vector<vec<I, 3>>& triangles, bool no_map = false) {
  if (!no_map) return get_edges(make_edge_map(triangles));
  auto edges = vector<vec<I, 2>>();
  for (auto& triangle : triangles) {
    edges.push_back(_make_edge(triangle.x, triangle.y));
    edges.push_back(_make_edge(triangle.y, triangle.z));
    edges.push_back(_make_edge(triangle.z, triangle.x));
  }
  return remove_duplicates(edges);
}
template <typename I>
inline vector<vec<I, 2>> get_edges(
    const vector<vec<I, 4>>& quads, bool no_map = false) {
  if (!no_map) return get_edges(make_edge_map(quads));
  auto edges = vector<vec<I, 2>>();
  for (auto& quad : quads) {
    edges.push_back(_make_edge(quad.x, quad.y));
    edges.push_back(_make_edge(quad.y, quad.z));
    if (quad.z != quad.w) edges.push_back(_make_edge(quad.z, quad.w));
    edges.push_back(_make_edge(quad.w, quad.x));
  }
  return remove_duplicates(edges);
}
template <typename I>
inline vector<vec<I, 2>> get_boundary(
    const vector<vec<I, 3>>& triangles, bool no_map = false) {
  if (!no_map) return get_boundary(make_edge_map(triangles));
}
template <typename I>
inline vector<vec<I, 2>> get_boundary(
    const vector<vec<I, 4>>& quads, bool no_map = false) {
  if (!no_map) return get_boundary(make_edge_map(quads));
}

// Build adjacencies between faces (sorted counter-clockwise)
template <typename I>
inline vector<vec<I, 3>> face_adjacencies(const vector<vec<I, 3>>& triangles);

// Build adjacencies between vertices (sorted counter-clockwise)
template <typename I>
inline vector<vector<int>> vertex_adjacencies(
    const vector<vec<I, 3>>& triangles, const vector<vec<I, 3>>& adjacencies);

// Compute boundaries as a list of loops (sorted counter-clockwise)
template <typename I>
inline vector<vector<int>> ordered_boundaries(
    const vector<vec<I, 3>>& triangles, const vector<vec<I, 3>>& adjacency,
    int num_vertices);

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
template <typename I>
inline vector<vector<int>> vertex_to_faces_adjacencies(
    const vector<vec<I, 3>>& triangles, const vector<vec<I, 3>>& adjacencies);

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
// and filled vectors for pos, norm and texcoord.
void split_facevarying(vector<vec4i>& split_quads,
    vector<vec3f>& split_positions, vector<vec3f>& split_normals,
    vector<vec2f>& split_texcoords, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords);

// Weld vertices within a threshold.
pair<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold);
pair<vector<vec3i>, vector<vec3f>> weld_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    float threshold);
pair<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold);

// Merge shape elements
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines by splitting each line in half.
pair<vector<vec2i>, vector<float>> subdivide_lines(
    const vector<vec2i>& lines, const vector<float>& vertex);
pair<vector<vec2i>, vector<vec2f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec2f>& vertex);
pair<vector<vec2i>, vector<vec3f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec3f>& vertex);
pair<vector<vec2i>, vector<vec4f>> subdivide_lines(
    const vector<vec2i>& lines, const vector<vec4f>& vertex);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
pair<vector<vec3i>, vector<float>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<float>& vertex);
pair<vector<vec3i>, vector<vec2f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec2f>& vertex);
pair<vector<vec3i>, vector<vec3f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& vertex);
pair<vector<vec3i>, vector<vec4f>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<vec4f>& vertex);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
pair<vector<vec4i>, vector<float>> subdivide_quads(
    const vector<vec4i>& quads, const vector<float>& vertex);
pair<vector<vec4i>, vector<vec2f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec2f>& vertex);
pair<vector<vec4i>, vector<vec3f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec3f>& vertex);
pair<vector<vec4i>, vector<vec4f>> subdivide_quads(
    const vector<vec4i>& quads, const vector<vec4f>& vertex);
// Subdivide beziers by splitting each segment in two.
pair<vector<vec4i>, vector<float>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<float>& vertex);
pair<vector<vec4i>, vector<vec2f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec2f>& vertex);
pair<vector<vec4i>, vector<vec3f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec3f>& vertex);
pair<vector<vec4i>, vector<vec4f>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<vec4f>& vertex);
// Subdivide quads using Carmull-Clark subdivision rules.
pair<vector<vec4i>, vector<float>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<float>& vertex,
    bool lock_boundary = false);
pair<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec2f>& vertex,
    bool lock_boundary = false);
pair<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec3f>& vertex,
    bool lock_boundary = false);
pair<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<vec4f>& vertex,
    bool lock_boundary = false);

// Subdivide lines by splitting each line in half.
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vertices, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
inline pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vertices, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vertices, int level);
// Subdivide beziers by splitting each segment in two.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vertices, int level);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vertices, int level,
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
void          sample_points_cdf(vector<float>& cdf, int npoints);

// Pick a point on lines uniformly.
pair<int, float> sample_lines(const vector<float>& cdf, float re, float ru);
vector<float>    sample_lines_cdf(
       const vector<vec2i>& lines, const vector<vec3f>& positions);
void sample_lines_cdf(vector<float>& cdf, const vector<vec2i>& lines,
    const vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
void sample_triangles_cdf(vector<float>& cdf, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions);
void sample_quads_cdf(vector<float>& cdf, const vector<vec4i>& quads,
    const vector<vec3f>& positions);

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
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMPUTATION OF PER-VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex tangents for lines.
template <typename T, typename I>
inline vector<vec<T, 3>> lines_tangents(
    const vector<vec<I, 2>>& lines, const vector<vec<T, 3>>& positions) {
  auto tangents = vector<vec<T, 3>>{positions.size()};
  for (auto& tangent : tangents) tangent = {0, 0, 0};
  for (auto& line : lines) {
    auto tangent = line_tangent(positions, line);
    auto length  = line_length(positions, line);
    for (auto vid : line) tangents[vid] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
  return tangents;
}

// Compute per-vertex normals for triangles.
template <typename T, typename I>
inline vector<vec<T, 3>> triangles_normals(
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions) {
  auto normals = vector<vec<T, 3>>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& triangle : triangles) {
    auto normal = triangle_normal(positions, triangle);
    auto area   = triangle_area(positions, triangle);
    for (auto vid : triangle) normals[vid] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Compute per-vertex normals for quads.
template <typename T, typename I>
inline vector<vec<T, 3>> quads_normals(
    const vector<vec<I, 4>>& quads, const vector<vec<T, 3>>& positions) {
  auto normals = vector<vec<T, 3>>{positions.size()};
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& quad : quads) {
    auto normal = quad_normal(positions, quad);
    auto area   = quad_area(positions, quad);
    if (!is_triangle(quad)) {
      for (auto vid : quad) normals[vid] += normal * area;
    } else {
      for (auto vid : as_triangle(quad)) normals[vid] += normal * area;
    }
  }
  for (auto& normal : normals) normal = normalize(normal);
  return normals;
}

// Compute per-vertex tangents for lines.
template <typename T, typename I>
inline void lines_tangents(vector<vec<T, 3>>& tangents,
    const vector<vec<I, 2>>& lines, const vector<vec<T, 3>>& positions) {
  check_same_size(tangents, positions);
  for (auto& tangent : tangents) tangent = {0, 0, 0};
  for (auto& line : lines) {
    auto tangent = line_tangent(positions, line);
    auto length  = line_length(positions, line);
    for (auto vid : line) tangents[vid] += tangent * length;
  }
  for (auto& tangent : tangents) tangent = normalize(tangent);
}

// Compute per-vertex normals for triangles.
template <typename T, typename I>
inline void triangles_normals(vector<vec<T, 3>>& normals,
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions) {
  check_same_size(normals, positions);
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& triangle : triangles) {
    auto normal = triangle_normal(positions, triangle);
    auto area   = triangle_area(positions, triangle);
    for (auto vid : triangle) normals[vid] += normal * area;
  }
  for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex normals for quads.
template <typename T, typename I>
inline void quads_normals(vector<vec<T, 3>>& normals,
    const vector<vec<I, 4>>& quads, const vector<vec<T, 3>>& positions) {
  check_same_size(normals, positions);
  for (auto& normal : normals) normal = {0, 0, 0};
  for (auto& quad : quads) {
    auto normal = quad_normal(positions, quad);
    auto area   = quad_area(positions, quad);
    if (!is_triangle(quad)) {
      for (auto vid : quad) normals[vid] += normal * area;
    } else {
      for (auto vid : as_triangle(quad)) normals[vid] += normal * area;
    }
  }
  for (auto& normal : normals) normal = normalize(normal);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
template <typename T, typename I>
inline vector<vec<T, 4>> triangles_tangent_spaces(
    const vector<vec<I, 3>>& triangles, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 2>>& texcoords) {
  auto tangu = vector<vec<T, 3>>(positions.size(), vec<T, 3>{0, 0, 0});
  auto tangv = vector<vec<T, 3>>(positions.size(), vec<T, 3>{0, 0, 0});
  for (auto& triangle : triangles) {
    auto tutv = triangle_tangents_fromuv(positions[triangle.x],
        positions[triangle.y], positions[triangle.z], texcoords[triangle.x],
        texcoords[triangle.y], texcoords[triangle.z]);
    for (auto vid : triangle) tangu[vid] += normalize(tutv.first);
    for (auto vid : triangle) tangv[vid] += normalize(tutv.second);
  }
  for (auto& tangent : tangu) tangent = normalize(tangent);
  for (auto& tangent : tangv) tangent = normalize(tangent);

  auto tangent_spaces = vector<vec<T, 4>>(positions.size());
  for (auto& tangent : tangent_spaces) tangent = vec<T, 4>{0, 0, 0, 0};
  for (auto vid : range(positions.size())) {
    tangu[vid] = orthonormalize(tangu[vid], normals[vid]);
    auto s     = dot(cross(normals[vid], tangu[vid]), tangv[vid]) < 0 ? -1.0f
                                                                      : 1.0f;
    tangent_spaces[vid] = {tangu[vid], s};
  }
  return tangent_spaces;
}

// Apply skinning
template <typename T, typename I>
inline pair<vector<vec<T, 3>>, vector<vec<T, 3>>> skin_vertices(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 4>>& weights, const vector<vec<I, 4>>& joints,
    const vector<frame3f>& xforms) {
  auto skinned_positions = vector<vec<T, 3>>{positions.size()};
  auto skinned_normals   = vector<vec<T, 3>>{positions.size()};
  for (auto i : range(positions.size())) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning as specified in Khronos glTF
template <typename T, typename I>
inline pair<vector<vec<T, 3>>, vector<vec<T, 3>>> skin_matrices(
    const vector<vec<T, 3>>& positions, const vector<vec<T, 3>>& normals,
    const vector<vec<T, 4>>& weights, const vector<vec<I, 4>>& joints,
    const vector<mat4f>& xforms) {
  auto skinned_positions = vector<vec<T, 3>>{positions.size()};
  auto skinned_normals   = vector<vec<T, 3>>{positions.size()};
  for (auto i : range(positions.size())) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
  return {skinned_positions, skinned_normals};
}

// Apply skinning
template <typename T, typename I>
inline void skin_vertices(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec<I, 4>>& joints, const vector<frame3f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i : range(positions.size())) {
    skinned_positions[i] =
        transform_point(xforms[joints[i].x], positions[i]) * weights[i].x +
        transform_point(xforms[joints[i].y], positions[i]) * weights[i].y +
        transform_point(xforms[joints[i].z], positions[i]) * weights[i].z +
        transform_point(xforms[joints[i].w], positions[i]) * weights[i].w;
    skinned_normals[i] = normalize(
        transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
        transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
        transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
        transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
  }
}

// Apply skinning as specified in Khronos glTF
template <typename T, typename I>
inline void skin_matrices(vector<vec<T, 3>>& skinned_positions,
    vector<vec<T, 3>>& skinned_normals, const vector<vec<T, 3>>& positions,
    const vector<vec<T, 3>>& normals, const vector<vec<T, 4>>& weights,
    const vector<vec<I, 4>>& joints, const vector<mat4f>& xforms) {
  if (skinned_positions.size() != positions.size() ||
      skinned_normals.size() != normals.size()) {
    throw std::out_of_range("arrays should be the same size");
  }
  for (auto i : range(positions.size())) {
    auto xform = xforms[joints[i].x] * weights[i].x +
                 xforms[joints[i].y] * weights[i].y +
                 xforms[joints[i].z] * weights[i].z +
                 xforms[joints[i].w] * weights[i].w;
    skinned_positions[i] = transform_point(xform, positions[i]);
    skinned_normals[i]   = normalize(transform_direction(xform, normals[i]));
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER VERTEX PROPETIES
// -----------------------------------------------------------------------------
namespace yocto {

// Flip vertex normals
template <typename T>
inline vector<vec<T, 3>> flip_normals(const vector<vec<T, 3>>& normals) {
  auto flipped = normals;
  for (auto& n : flipped) n = -n;
  return flipped;
}
// Flip face orientation
template <typename I>
inline vector<vec<I, 3>> flip_triangles(const vector<vec<I, 3>>& triangles) {
  auto flipped = triangles;
  for (auto& [v1, v2, v3] : flipped) swap(v2, v3);
  return flipped;
}
template <typename I>
inline vector<vec<I, 4>> flip_quads(const vector<vec<I, 4>>& quads) {
  auto flipped = quads;
  for (auto& [v1, v2, v3, v4] : flipped) {
    if (v3 != v4) {
      swap(v2, v4);
    } else {
      swap(v2, v3);
      v4 = v3;
    }
  }
  return flipped;
}

// Align vertex positions. Alignment is 0: none, 1: min, 2: max, 3: center.
template <typename T, typename I>
inline vector<vec<T, 3>> align_vertices(
    const vector<vec<T, 3>>& positions, const vec<I, 3>& alignment) {
  auto bounds = invalidb3f;
  for (auto& p : positions) bounds = merge(bounds, p);
  auto offset = vec<T, 3>{0, 0, 0};
  for (auto a : range(3)) {
    switch (alignment[a]) {
      case 0: break;
      case 1: offset[a] = bounds.min[a]; break;
      case 2: offset[a] = (bounds.min[a] + bounds.max[a]) / 2; break;
      case 3: offset[a] = bounds.max[a]; break;
      default: throw std::invalid_argument{"invalid alignment"};
    }
  }
  auto aligned = positions;
  for (auto& p : aligned) p -= offset;
  return aligned;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines by splitting each line in half.
template <typename T>
inline pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vertices, int level) {
  if (level < 1) return {lines, vertices};
  auto tess = pair{lines, vertices};
  for (auto idx : range(level)) tess = subdivide_lines(tess.first, tess.second);
  return tess;
}
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
inline pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vertices, int level) {
  if (level < 1) return {triangles, vertices};
  auto tess = pair{triangles, vertices};
  for (auto idx : range(level))
    tess = subdivide_triangles(tess.first, tess.second);
  return tess;
}
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vertices, int level) {
  if (level < 1) return {quads, vertices};
  auto tess = pair{quads, vertices};
  for (auto idx : range(level)) tess = subdivide_quads(tess.first, tess.second);
  return tess;
}
// Subdivide beziers by splitting each segment in two.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vertices, int level) {
  if (level < 1) return {beziers, vertices};
  auto tess = pair{beziers, vertices};
  for (auto idx : range(level))
    tess = subdivide_beziers(tess.first, tess.second);
  return tess;
}
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
inline pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vertices, int level,
    bool lock_boundary) {
  if (level < 1) return {quads, vertices};
  auto tess = pair{quads, vertices};
  for (auto idx : range(level))
    tess = subdivide_catmullclark(tess.first, tess.second, lock_boundary);
  return tess;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EDGEA AND ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Initialize an edge map with elements.
template <typename I>
inline edge_map<I> make_edge_map(const vector<vec<I, 3>>& triangles) {
  auto emap = edge_map<I>{};
  for (auto& [v1, v2, v3] : triangles) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    insert_edge(emap, {v3, v1});
  }
  return emap;
}
template <typename I>
inline edge_map<I> make_edge_map(const vector<vec<I, 4>>& quads) {
  auto emap = edge_map<I>{};
  for (auto& [v1, v2, v3, v4] : quads) {
    insert_edge(emap, {v1, v2});
    insert_edge(emap, {v2, v3});
    if (v3 != v4) insert_edge(emap, {v3, v4});
    insert_edge(emap, {v4, v1});
  }
  return emap;
}
// Insert an edge and return its index
template <typename I>
inline I insert_edge(edge_map<I>& emap, const vec<I, 2>& edge) {
  auto es = vec<I, 2>{min(edge), max(edge)};
  auto it = emap.edges.find(es);
  if (it == emap.edges.end()) {
    auto index = (I)emap.edges.size();
    emap.edges.insert(it, {es, {index, 1}});
    return index;
  } else {
    auto& data = it->second;
    data.nfaces += 1;
    return data.index;
  }
}
// Get number of edges
template <typename I>
inline I num_edges(const edge_map<I>& emap) {
  return (I)emap.edges.size();
}
// Get the edge index
template <typename I>
inline I edge_index(const edge_map<I>& emap, const vec<I, 2>& edge) {
  auto es       = vec<I, 2>{min(edge), max(edge)};
  auto iterator = emap.edges.find(es);
  if (iterator == emap.edges.end()) return -1;
  return iterator->second.index;
}
// Get a list of edges, boundary edges, boundary vertices
template <typename I>
inline vector<vec<I, 2>> get_edges(const edge_map<I>& emap) {
  auto edges = vector<vec<I, 2>>(emap.edges.size());
  for (auto& [edge, data] : emap.edges) edges[data.index] = edge;
  return edges;
}
template <typename I>
inline vector<vec<I, 2>> get_boundary(const edge_map<I>& emap) {
  auto boundary = vector<vec<I, 2>>{};
  for (auto& [edge, data] : emap.edges) {
    if (data.nfaces < 2) boundary.push_back(edge);
  }
  return boundary;
}

// Build adjacencies between faces (sorted counter-clockwise)
template <typename I>
inline vector<vec<I, 3>> face_adjacencies(const vector<vec<I, 3>>& triangles) {
  auto get_edge = [](const vec<I, 3>& triangle, I i) -> vec<I, 2> {
    auto x = triangle[i], y = triangle[i < 2 ? i + 1 : 0];
    return x < y ? vec<I, 2>{x, y} : vec<I, 2>{y, x};
  };
  auto adjacencies = vector<vec<I, 3>>{triangles.size(), vec<I, 3>{-1, -1, -1}};
  auto edge_map    = unordered_map<vec<I, 2>, I>();
  edge_map.reserve((size_t)(triangles.size() * 1.5));
  for (auto i = 0; i < (I)triangles.size(); ++i) {
    for (auto k = 0; k < 3; ++k) {
      auto edge = get_edge(triangles[i], k);
      auto it   = edge_map.find(edge);
      if (it == edge_map.end()) {
        edge_map.insert(it, {edge, i});
      } else {
        auto neighbor     = it->second;
        adjacencies[i][k] = neighbor;
        for (auto kk = 0; kk < 3; ++kk) {
          auto edge2 = get_edge(triangles[neighbor], kk);
          if (edge2 == edge) {
            adjacencies[neighbor][kk] = i;
            break;
          }
        }
      }
    }
  }
  return adjacencies;
}

// Build adjacencies between vertices (sorted counter-clockwise)
template <typename I>
inline vector<vector<I>> vertex_adjacencies(
    const vector<vec<I, 3>>& triangles, const vector<vec<I, 3>>& adjacencies) {
  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<I>(triangles.size() * 3, -1);

  for (auto i = 0; i < (I)triangles.size(); ++i) {
    for (auto k : range(3)) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<I>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (auto i = 0; i < num_vertices; ++i) {
    result[i].reserve(6);
    auto first_face = face_from_vertex[i];
    if (first_face == -1) continue;

    auto face = first_face;
    while (true) {
      auto k = find_index(triangles[face], i);
      k      = k != 0 ? k - 1 : 2;
      result[i].push_back(triangles[face][k]);
      face = adjacencies[face][k];
      if (face == -1) break;
      if (face == first_face) break;
    }
  }

  return result;
}

// Build adjacencies between each vertex and its adjacent faces.
// Adjacencies are sorted counter-clockwise and have same starting points as
// vertex_adjacencies()
template <typename I>
inline vector<vector<I>> vertex_to_faces_adjacencies(
    const vector<vec<I, 3>>& triangles, const vector<vec<I, 3>>& adjacencies) {
  // For each vertex, find any adjacent face.
  auto num_vertices     = 0;
  auto face_from_vertex = vector<I>(triangles.size() * 3, -1);

  for (auto i = 0; i < (I)triangles.size(); ++i) {
    for (auto k : range(3)) {
      face_from_vertex[triangles[i][k]] = i;
      num_vertices                      = max(num_vertices, triangles[i][k]);
    }
  }

  // Init result.
  auto result = vector<vector<I>>(num_vertices);

  // For each vertex, loop around it and build its adjacency.
  for (auto i = 0; i < num_vertices; ++i) {
    result[i].reserve(6);
    auto first_face = face_from_vertex[i];
    if (first_face == -1) continue;

    auto face = first_face;
    while (true) {
      auto k = find_index(triangles[face], i);
      k      = k != 0 ? k - 1 : 2;
      face   = adjacencies[face][k];
      result[i].push_back(face);
      if (face == -1) break;
      if (face == first_face) break;
    }
  }

  return result;
}

// Compute boundaries as a list of loops (sorted counter-clockwise)
template <typename I>
inline vector<vector<I>> ordered_boundaries(const vector<vec<I, 3>>& triangles,
    const vector<vec<I, 3>>& adjacency, I num_vertices) {
  // map every boundary vertex to its next one
  auto next_vert = vector<I>(num_vertices, -1);
  for (auto i = 0; i < (I)triangles.size(); ++i) {
    for (auto k = 0; k < 3; ++k) {
      if (adjacency[i][k] == -1)
        next_vert[triangles[i][k]] = triangles[i][(k + 1) % 3];
    }
  }

  // result
  auto boundaries = vector<vector<I>>();

  // arrange boundary vertices in loops
  for (auto i : range(next_vert.size())) {
    if (next_vert[i] == -1) continue;

    // add new empty boundary
    boundaries.emplace_back();
    auto current = (I)i;

    while (true) {
      auto next = next_vert[current];
      if (next == -1) {
        return {};
      }
      next_vert[current] = -1;
      boundaries.back().push_back(current);

      // close loop if necessary
      if (next == i)
        break;
      else
        current = next;
    }
  }

  return boundaries;
}

}  // namespace yocto

#endif
