//
// # Yocto/Mesh: Tiny Library for mesh operations for computational geometry
//
// Yocto/Mesh is a collection of computational geometry routines on triangle
// meshes.
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

#ifndef _YOCTO_MESH_H_
#define _YOCTO_MESH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// ADJACENCIES
// -----------------------------------------------------------------------------
namespace yocto {

// Returns the list of triangles incident at each vertex in CCW order.
// Note: this works only if the mesh does not have a boundary.
vector<vector<int>> vertex_to_triangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies);

// Face adjacent to t and opposite to vertex vid
int opposite_face(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, int t, int vid);

// Finds the opposite vertex of an edge
int opposite_vertex(const vec3i& triangle, const vec2i& edge);
int opposite_vertex(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, int face, int k);
int common_vertex(const vector<vec3i>& triangles, int pid0, int pid1);

// Finds common edge between triangles
vec2i common_edge(const vec3i& triangle0, const vec3i& triangle1);
vec2i opposite_edge(const vec3i& t, int vid);
vec2i common_edge(const vector<vec3i>& triangles, int pid0, int pid1);

// Triangle fan starting from a face and going towards the k-th neighbor face.
vector<int> triangle_fan(
    const vector<vec3i>& adjacencies, int face, int k, bool clockwise = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE GEODESICS
// -----------------------------------------------------------------------------
namespace yocto {

// Data structure used for geodesic computation
struct geodesic_solver {
  static const int min_arcs = 12;
  struct graph_edge {
    int   node   = -1;
    float length = flt_max;
  };
#ifdef YOCTO_ABSEIL
  vector<short_vector<adjancency_list, min_arcs>> graph = {};
#else
  vector<vector<graph_edge>> graph = {};
#endif
};

// Construct a graph to compute geodesic distances
geodesic_solver make_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3i>& adjacencies, const vector<vec3f>& positions);

// Compute angles in tangent space and total angles of every vertex
vector<vector<float>> compute_angles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& v2t, vector<float>& total_angles,
    bool with_opposite);
// Construct a graph to compute geodesic distances with arcs arranged in
// counterclockwise order by using the vertex-to-face adjacencies
geodesic_solver make_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<vector<int>>& vertex_to_faces);

// Compute geodesic distances
void update_geodesic_distances(vector<float>& distances,
    const geodesic_solver& solver, const vector<int>& sources,
    float max_distance = flt_max);

vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<int>& sources, float max_distance = flt_max);

// Compute all shortest paths from source vertices to any other vertex.
// Paths are implicitly represented: each node is assignes its previous node
// in the path. Graph search early exits when reching end_vertex.
vector<int> compute_geodesic_paths(const geodesic_solver& solver,
    const vector<int>& sources, int end_vertex = -1);

// Sample vertices with a Poisson distribution using geodesic distances.
// Sampling strategy is farthest point sampling (FPS): at every step
// take the farthers point from current sampled set until done.
vector<int> sample_vertices_poisson(
    const geodesic_solver& solver, int num_samples);

// Compute the distance field needed to compute a voronoi diagram
vector<vector<float>> compute_voronoi_fields(
    const geodesic_solver& solver, const vector<int>& generators);

// Convert distances to colors
vector<vec3f> colors_from_field(const vector<float>& field, float scale = 1,
    const vec3f& c0 = {1, 1, 1}, const vec3f& c1 = {1, 0.1f, 0.1f});

// Description of a discrete path along the surface of a triangle mesh.
struct surface_path {
  struct vertex {
    vec2i edge  = {0, 0};
    int   face  = 0;
    float alpha = 0;
  };
  int            start, end;
  vector<vertex> vertices;
};

// Trace integral path following the gradient of a scalar field
surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from);
surface_path integrate_field(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacency,
    const vector<int>& tags, int tag, const vector<float>& field, int from,
    int to);

vector<vec3f> make_positions_from_path(
    const surface_path& path, const vector<vec3f>& mesh_positions);

vec3f compute_gradient(const vec3i& triangle, const vector<vec3f>& positions,
    const vector<float>& field);

struct mesh_point {
  int   face = -1;
  vec2f uv   = {0, 0};
};

// compute geodesic distance from  a source Point to all the vertices of the
// mesh(triangles,positions,adjacencies)
vector<float> compute_geodesic_distances(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<mesh_point>& sources);

// given a mesh(triangles,positions,adjacencies), computes the list of parents
// from Point target to Point source (discrete shortest path in the graph)
vector<int> point_to_point_geodesic_path(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const mesh_point& source,
    const mesh_point& target);

// Dual geodesic graph
struct dual_geodesic_solver {
  struct edge {
    int   node   = -1;
    float length = flt_max;
  };
  vector<array<edge, 3>> graph = {};
};

// Construct a graph to compute geodesic distances
dual_geodesic_solver make_dual_geodesic_solver(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies);

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEODESIC PATH
// -----------------------------------------------------------------------------
namespace yocto {

struct geodesic_path {
  // surface data
  mesh_point    start = {};
  mesh_point    end   = {};
  vector<int>   strip = {};
  vector<float> lerps = {};
};

// compute the shortest path connecting two surface points
// initial guess of the connecting strip must be given
geodesic_path shortest_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& start, const mesh_point& end, const vector<int>& strip);

// compute the straightest path given a surface point and tangent direction
geodesic_path straightest_path(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& start, const vec2f& direction, float path_length);

vector<vec3f> path_positions(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies);
vector<float> path_parameters(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies);
float path_length(const geodesic_path& path, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies);
vector<float> path_parameters(const vector<vec3f>& positions);
float         path_length(const vector<vec3f>& positions);

mesh_point eval_path_point(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, float t);

inline mesh_point eval_path_midpoint(const geodesic_path& path,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies) {
  return eval_path_point(path, triangles, positions, adjacencies, 0.5);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STRIPS
// -----------------------------------------------------------------------------
namespace yocto {

vector<int> strip_on_dual_graph(const dual_geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions, int start,
    int end);

vector<int> strip_ascending_distance_field(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int start, int end);

// returns a strip of triangles such target belongs to the first one and
// source to the last one
// TODO(fabio): the name should change in order to get the call
// more consistent with the output)
vector<int> get_strip(const geodesic_solver& solver,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, const vector<vector<int>>& v2t,
    const vector<vector<float>>& angles, const mesh_point& source,
    const mesh_point& target);

}  // namespace yocto

// -----------------------------------------------------------------------------
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

using unfold_triangle = std::array<vec2f, 3>;

vec3f eval_position(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& sample);
vec3f eval_normal(const vector<vec3i>& triangles, const vector<vec3f>& normals,
    const mesh_point& point);

pair<bool, int>   bary_is_edge(const vec3f& bary, float tol = 1e-2f);
pair<bool, int>   bary_is_vert(const vec3f& bary, float tol = 1e-2f);
pair<bool, int>   point_is_vert(const mesh_point& p, float tol = 1e-2f);
pair<bool, int>   point_is_edge(const mesh_point& p, float tol = 5e-3f);
pair<bool, vec2f> point_in_triangle(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, int tid, const vec3f& point,
    float tol = 5e-3f);

std::array<vec2f, 3> init_flat_triangle(
    const vector<vec3f>& positions, const vec3i& tr);

float length_by_flattening(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const mesh_point& p, const vector<int>& strip);

// intersection of two circles with centers c1, c2 and radii squared R1, R2
vec2f intersect_circles(const vec2f& c2, float R2, const vec2f& c1, float R1);

// given the 2D coordinates in tanget space of a triangle, find the coordinates
// of the k-th neighbor triangle
unfold_triangle unfold_face(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const unfold_triangle& tr, int face, int k);

// assign 2D coordinates to a strip of triangles. point start is at (0, 0)
vector<unfold_triangle> unfold_strip(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3i>& adjacencies,
    const vector<int>& strip, const mesh_point& start);

// assign 2D coordinates to vertices of the triangle containing the mesh point,
// putting the point at (0, 0)
unfold_triangle triangle_coordinates(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const mesh_point& point);

// generic utilities for paths
vec2i get_edge(const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3i>& adjacencies, int f0, int f1);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PROCEDURAL MODELING
// -----------------------------------------------------------------------------
namespace yocto {

// Extract isoline from surface scalar field.
void meandering_triangles(const vector<float>& field, float isoline,
    int selected_tag, int t0, int t1, vector<vec3i>& triangles,
    vector<int>& tags, vector<vec3f>& positions, vector<vec3f>& normals);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a shape as indexed meshes
[[nodiscard]] bool load_mesh(const string& filename, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_mesh(const string& filename,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, string& error, bool ascii = false,
    bool flip_texcoords = true);

// Load/save a set of lines
[[nodiscard]] bool load_lines(const string& filename, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<vec3f>& colors, string& error, bool flip_texcoords = true);
[[nodiscard]] bool save_lines(const string& filename,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<vec3f>& colors, string& error, bool ascii = false,
    bool flip_texcoords = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

// Get mesh statistics for printing
vector<string> mesh_stats(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<vec3f>& colors,
    bool verbose = false);

}  // namespace yocto

#endif
