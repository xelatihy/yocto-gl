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

#ifndef YOCTO_SHAPE_H_
#define YOCTO_SHAPE_H_

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
#include "yocto_image.h"
#include "yocto_math.h"
#include "yocto_ndarray.h"
#include "yocto_sampling.h"

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

// Displacement
shape_data displace_shape(const shape_data& shape,
    const array2d<float>& displacement, float height = 1, float offset = 0.5f);
shape_data displace_shape(const shape_data& shape,
    const array2d<vec4f>& displacement, float height = 1, float offset = 0.5f);

// Transform shape
shape_data transform_shape(
    const shape_data& shape, const frame3f& frame, bool non_rigid = false);
shape_data transform_shape(const shape_data& shape, const frame3f& frame,
    float radius_scale, bool non_rigid = false);
shape_data scale_shape(const shape_data& shape, float scale, float uvscale = 1);
shape_data scale_shape(shape_data&& shape, float scale, float uvscale = 1);
shape_data flipyz_shape(const shape_data& shape);

// Manipulate vertex data
shape_data remove_normals(const shape_data& shape);
shape_data add_normals(const shape_data& shape);
shape_data weld_vertices(const shape_data& shape, float threshold,
    bool normals = false, bool others = false);

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
vec2f eval_texcoord(const fvshape_data& shape, int element, const vec2f& uv);

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
    const fvshape_data& shape, const frame3f& frame, bool non_rigid = false);
fvshape_data scale_fvshape(
    const fvshape_data& shape, float scale, float uvscale = 1);
fvshape_data scale_fvshape(
    fvshape_data&& shape, float scale, float uvscale = 1);

// Vertex properties
fvshape_data remove_normals(const fvshape_data& shape);
fvshape_data add_normals(const fvshape_data& shape);

// Shape statistics
vector<string> fvshape_stats(const fvshape_data& shape, bool verbose = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a quad
shape_data make_quad(int steps = 1, float scale = 1, float uvscale = 1);
// Make a quad along y
shape_data make_quady(int steps = 1, float scale = 1, float uvscale = 1);
// Make a rectangle
shape_data make_rect(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
// Make a rectangle along y
shape_data make_recty(const vec2i& steps = {1, 1}, const vec2f& scale = {1, 1},
    const vec2f& uvscale = {1, 1});
// Make a floor
shape_data make_floor(const vec2i& steps = {1, 1},
    const vec2f& scale = {20, 20}, const vec2f& uvscale = {20, 20});
// Make a disk
shape_data make_disk(int steps = 32, float scale = 1, float uvscale = 1);
// Make a cube
shape_data make_cube(int steps = 1, float scale = 1, float uvscale = 1);
// Make a cube that is watertight (only positions)
shape_data make_wtcube(int steps = 1, float scale = 1);
// Make a cube that is watertight and open (only positions)
shape_data make_opcube(int steps = 1, float scale = 1);
// Make a box
shape_data make_box(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a watertight box (positions-only)
shape_data make_wtbox(
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1});
// Make a watertight and open box (positions-only)
shape_data make_opbox(
    const vec3i& steps = {1, 1, 1}, const vec3f& scale = {1, 1, 1});
// Make a tesselated sphere
shape_data make_sphere(int steps = 32, float scale = 1, float uvscale = 1);
// Make a uv sphere
shape_data make_uvsphere(const vec2i& steps = {64, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv sphere aligned to y
shape_data make_uvspherey(const vec2i& steps = {64, 32}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv disk
shape_data make_uvdisk(const vec2i& steps = {32, 1}, float scale = 1,
    const vec2f& uvscale = {1, 1});
// Make a uv cylinder
shape_data make_uvcylinder(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a uv capsule
shape_data make_uvcapsule(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a uv cone
shape_data make_uvcone(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a bulged rectangle
shape_data make_bulged_rect(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, float height = 0.3f,
    const vec2f& uvscale = {1, 1});
// Make a bulged rectangle along y
shape_data make_bulged_recty(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, float height = 0.3f,
    const vec2f& uvscale = {1, 1});
// Make a bent disk
shape_data make_bent_floor(const vec2i& steps = {1, 1},
    const vec2f& scale = {20, 20}, float radius = 0.5f,
    const vec2f& uvscale = {20, 20});
// Make a bulged disk
shape_data make_bulged_disk(
    int steps = 32, float scale = 1, float height = 0.3f, float uvscale = 1);
// Make a rounded box
shape_data make_rounded_box(const vec3i& steps = {32, 32, 32},
    const vec3f& scale = {1, 1, 1}, float radius = 0.3f,
    const vec3f& uvscale = {1, 1, 1});
// Make a capped uv sphere
shape_data make_capped_uvsphere(const vec2i& steps = {32, 32}, float scale = 1,
    float height = 0.3f, const vec2f& uvscale = {1, 1});
// Make a capped uv spherey
shape_data make_capped_uvspherey(const vec2i& steps = {32, 32}, float scale = 1,
    float height = 0.3f, const vec2f& uvscale = {1, 1});
// Make a rounded uv cylinder
shape_data make_rounded_uvcylinder(const vec3i& steps = {32, 32, 32},
    const vec2f& scale = {1, 1}, float radius = 0.3f,
    const vec3f& uvscale = {1, 1, 1});

// Make a geosphere
shape_data make_geosphere(int subdivisions = 3, float scale = 1);
// Make a monkey
shape_data make_monkey(int subdivisions = 2, float scale = 1);
// Make a subdiv cube
shape_data make_sdcube(int subdivisions = 4, float scale = 1);
// Make a subdiv cube with texture coordinates
shape_data make_sdtcube(int subdivisions = 4, float scale = 1);

// Make a face-varying quad
fvshape_data make_fvquad(int steps = 1, float scale = 1, float uvscale = 1);
// Make a face-varying rectangle
fvshape_data make_fvrect(const vec2i& steps = {1, 1},
    const vec2f& scale = {1, 1}, const vec2f& uvscale = {1, 1});
// Make a face-varying cube
fvshape_data make_fvcube(int steps = 1, float scale = 1, float uvscale = 1);
// Make a face-varying box
fvshape_data make_fvbox(const vec3i& steps = {1, 1, 1},
    const vec3f& scale = {1, 1, 1}, const vec3f& uvscale = {1, 1, 1});
// Make a face-varying sphere
fvshape_data make_fvsphere(int steps = 32, float scale = 1, float uvscale = 1);

// Make a point.
shape_data make_point(float radius = 0.001f);
// Make many points at the origin (useful for particle systems).
shape_data make_points(
    int num = 65536, float radius = 0.001f, float uvscale = 1);

// Make a line.
shape_data make_line(int steps = 4, float scale = 1,
    const vec2f& radius = {0.001f, 0.001f}, float uvscale = 1);
// Make lines along a quad.
shape_data make_lines(int num = 65536, int steps = 4,
    const vec2f& scale = {1, 1}, const vec2f& radius = {0.001f, 0.001f},
    const vec2f& uvscale = {1, 1});

// Make a grid of quads
shape_data make_quad_grid(const vec2i& steps = {256, 256}, float scale = 1);
// Make points on a grid.
shape_data make_point_grid(
    const vec2i& steps = {256, 256}, float scale = 1, float radius = 0.001f);
// Make lines on a grid.
shape_data make_line_grid(
    const vec2i& steps = {256, 256}, float scale = 1, float radius = 0.001f);

// Make random points in a qiad or cube.
shape_data make_random_points(int num = 65536, float radius = 0.001f,
    bool on_quad = false, uint64_t seed = 17);
// Make points on a shape
shape_data make_random_points(const shape_data& shape, int num = 65536,
    float radius = 0.001f, uint64_t seed = 17);
// Make lines on a shape
shape_data make_random_lines(const shape_data& shape, int num = 65536,
    int steps = 8, const vec2f& length = {0.1f, 0.1f},
    const vec2f& radius = {0.001f, 0.001f}, uint64_t seed = 17);
// Make hairs on a shape
shape_data make_random_hairs(const shape_data& shape, int num = 65536,
    int steps = 8, const vec2f& length = {0.35f, 0.35f},
    const vec2f& radius = {0.001f, 0.001f}, float noise = 0.1f,
    float gravity = 0.001f, uint64_t seed = 17);

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
shape_data make_heightfield(
    const array2d<float>& height, const vec2f& scale = {1, 1});
shape_data make_heightfield(
    const array2d<vec4f>& height, const vec2f& scale = {1, 1});

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF PER_VERTEX PROPERTIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> lines_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions);
vector<vec3f> triangles_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
vector<vec3f> quads_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions);
// Update normals and tangents
void lines_tangents(vector<vec3f>& tangents, const vector<vec2i>& lines,
    const vector<vec3f>& positions);
void triangles_normals(vector<vec3f>& normals, const vector<vec3i>& triangles,
    const vector<vec3f>& positions);
void quads_normals(vector<vec3f>& normals, const vector<vec4i>& quads,
    const vector<vec3f>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> triangle_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords);

// Apply skinning to vertex position and normals.
pair<vector<vec3f>, vector<vec3f>> skin_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
pair<vector<vec3f>, vector<vec3f>> skin_matrices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms);
// Update skinning
void skin_vertices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms);
void skin_matrices(vector<vec3f>& skinned_positions,
    vector<vec3f>& skinned_normals, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms);

}  // namespace yocto

// -----------------------------------------------------------------------------
// COMPUTATION OF VERTEX PROPERTIES
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
// SHAPE ELEMENT CONVERSION
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
pair<int, float> sample_lines(const vector<float>& cdf, float re, float ru);
vector<float>    sample_lines_cdf(
       const vector<vec2i>& lines, const vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv);
pair<int, vec2f> sample_quads(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv);
vector<float>    sample_quads_cdf(
       const vector<vec4i>& quads, const vector<vec3f>& positions);

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
struct edge_map_edge {
  int index, nfaces;
};

// Dictionary to store edge information. `index` is the index to the edge
// array, `edges` the array of edges and `nfaces` the number of adjacent faces.
// We store only bidirectional edges to keep the dictionary small. Use the
// functions below to access this data.
struct edge_map {
  unordered_map<vec2i, edge_map_edge> edges = {};
};

// Initialize an edge map with elements.
edge_map make_edge_map(const vector<vec3i>& triangles);
edge_map make_edge_map(const vector<vec4i>& quads);

// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge);

// Get the edge index
int edge_index(const edge_map& emap, const vec2i& edge);
// Get edges and boundaries
int           num_edges(const edge_map& emap);
vector<vec2i> get_edges(const edge_map& emap);
vector<vec2i> get_boundary(const edge_map& emap);

// Edges helpers
vector<vec2i> get_edges(const vector<vec3i>& triangles, bool no_map = false);
vector<vec2i> get_edges(const vector<vec4i>& quads, bool no_map = false);
vector<vec2i> get_boundary(const vector<vec3i>& triangles);
vector<vec2i> get_boundary(const vector<vec4i>& quads);

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
// SHAPE SUBDIVISION
// -----------------------------------------------------------------------------
namespace yocto {

// Subdivide lines by splitting each line in half.
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vert);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vert);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vert);
// Subdivide beziers by splitting each segment in two.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec3i>& beziers, const vector<T>& vert);
// Subdivide lines using B-splines subdivision rules.
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vert);
// Subdivide quads using Catmull-Clark subdivision rules.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vert,
    bool lock_boundary = false);

// Subdivide lines by splitting each line in half.
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vertices, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
pair<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vertices, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vertices, int level);
// Subdivide beziers by splitting each segment in two.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vertices, int level);
// Subdivide lines using B-splines subdivision rules.
template <typename T>
pair<vector<vec2i>, vector<T>> subdivide_bspline(
    const vector<vec2i>& lines, const vector<T>& vert, int level);
// Subdivide quads using Catmull-Clark subdivision rules.
template <typename T>
pair<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vertices, int level,
    bool lock_boundary = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE DISPLACEMENT
// -----------------------------------------------------------------------------
namespace yocto {

// Displace vertices
vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const array2d<float>& displacement, float scale = 1, float offset = 0.5);
vector<vec3f> displace_vertices(const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const array2d<vec4f>& displacement, float scale = 1, float offset = 0.5);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE ELEMENT GROUPING
// -----------------------------------------------------------------------------
namespace yocto {

// Weld vertices within a threshold.
pair<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold);
pair<vector<int>, vector<int>> weld_indices(
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

#endif
