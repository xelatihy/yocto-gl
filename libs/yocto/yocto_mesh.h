//
// # Yocto/Mesh: Tiny Library for mesh operations for computational geometry
//
// Yocto/Mesh is a collection of computational geometry routines on triangle
// meshes.
//
//
// ## Mesh IO
//
// We support reading and writing shapes in OBJ and PLY.
//
// 1. load/save meshes with `load_mesh()`/`save_mesh()`
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

// Finds common edge between triangles
vec2i common_edge(const vec3i& triangle0, const vec3i& triangle1);

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
    const vec3f& c0 = {1, 1, 1}, const vec3f& c1 = {1, 0.1, 0.1});

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
