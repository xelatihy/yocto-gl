//
// # Yocto/Shape: Tiny C++ Library for 3D Shape utilities
//
// Yocto/Shape is a C++ library for 3D shape manipulations with a focus on
// shapes represented as index triangle and quad meshes, indexed line and
// point sets and indexed bezier sets. The utlities collected here are written
// to support a global illujmination rendering and not for generic geometry
// processing.. We support operation for shape smoothing, shape subdivision
// (including Catmull-Clark subdivs), and example shape creation.
//
// ## Usage
//
// 1. compute line tangents, and triangle and quad areas and normals
//    (with Yocto/Math)
// 2. interpolate values over primitives with `eval_line()`,
//    `eval_triangle()` and `eval_quad()`
// 3. evaluate Bezier curves and derivatives with `eval_bezier()` and
//    `eval_bezier_derivative()`
// 4. compute smooth normals and tangents with `compute_normals()`
//   `compute_tangents()`
// 5. compute tangent frames from texture coordinates with
//    `compute_tangent_space()`
// 6. compute skinning with `compute_skinning()` and
//    `compute_matrix_skinning()`
// 6. create shapes with `make_cube()`, `make_sphere()`, `make_quad()`,
//    `make_fvcube()`, `make_hair()`, `make_suzanne()`, `make_lines()`,
//    `make_points()`, `make_sphere_cube()`, `make_cube_rounded()`,
//    `make_sphere_flipcap()`, `make_cylinder()`, `make_cylinder_rounded()`,
//    `make_disk()`, `make_cylinder_side()`, `make_disk_quad()`
// 7. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
// 8. shape sampling with `sample_points()`, `sample_lines()`,
//    `sample_triangles()`; initialize the sampling CDFs with
//    `sample_points_cdf()`, `sample_lines_cdf()`, `sample_triangles_cdf()`
// 9.  samnple a could of point over a surface with `sample_triangles_points()`
// 10. get edges and boundaries with `get_edges()`
// 11. convert quads to triangles with `convert_quads_to_triangles()`
// 12. convert face varying to vertex shared representations with
//     `convert_face_varying()`
// 13. subdivide elements by edge splits with `subdivide_lines()`,
//     `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
// 14. Catmull-Clark subdivision surface with `subdivide_catmullclark()`
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
#ifndef _YGL_SHAPE_H_
#define _YGL_SHAPE_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS AND INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex normals/tangents for lines/triangles/quads.
void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    bool weighted = true);
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true);
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_space(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted = true);

// Apply skinning to vertex position and normals.
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
// Apply skinning as specified in Khronos glTF.
void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles);
std::vector<vec2i> get_edges(const std::vector<vec4i>& quads);

// Convert quads to triangles
std::vector<vec3i> convert_quads_to_triangles(const std::vector<vec4i>& quads);
// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
std::vector<vec3i> convert_quads_to_triangles(
    const std::vector<vec4i>& quads, int row_length);

// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec2i> convert_bezier_to_lines(const std::vector<vec4i>& beziers);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm and texcoord.
void convert_face_varying(std::vector<vec4i>& qquads, std::vector<vec3f>& qpos,
    std::vector<vec3f>& qnorm, std::vector<vec2f>& qtexcoord,
    std::vector<vec4f>& qcolor, const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec4i>& quads_color, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color);

// Subdivide lines by splitting each line in half.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
void subdivide_triangles(
    std::vector<vec3i>& triangles, std::vector<T>& vert, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, int level);
// Subdivide beziers by splitting each segment in two.
template <typename T>
void subdivide_beziers(
    std::vector<vec4i>& beziers, std::vector<T>& vert, int level);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<T>& vert,
    int level, bool lock_boundary = false);

// Merge lines between shapes.
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);
// Merge triangles between shapes.
void merge_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec3i>& triangles1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);
// Merge quads between shapes.
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec4i>& quads1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);

// Weld vertices within a threshold. For noe the implementation is O(n^2).
std::vector<int> weld_vertices(std::vector<vec3f>& pos, float threshold);
void weld_triangles(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, float threshold);
void weld_quads(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, float threshold);

// Pick a point in a point set uniformly.
inline int sample_points(int npoints, float re) {
    return sample_index(npoints, re);
}
inline std::vector<float> sample_points_cdf(int npoints) {
    auto cdf = std::vector<float>(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
    return cdf;
}
inline int sample_points(const std::vector<float>& cdf, float re) {
    return sample_discrete(cdf, re);
}

// Pick a point on lines uniformly.
inline std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(pos[l.x], pos[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, float> sample_lines(
    const std::vector<float>& cdf, float re, float ru) {
    return {sample_discrete(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
inline std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(pos[t.x], pos[t.y], pos[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, vec2f> sample_triangles(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), sample_triangle(ruv)};
}

// Pick a point on a quad mesh uniformly.
inline std::vector<float> sample_quads_cdf(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(pos[q.x], pos[q.y], pos[q.z], pos[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, vec2f> sample_quads(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), ruv};
}

// Samples a set of points over a triangle mesh uniformly. Returns pos, norm
// and tecoord of the sampled points.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, int seed = 7);

// Make examples triangle shapes with shared vertices (not watertight).
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize);
void make_quad_stack(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec2f& uvsize);
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize);
void make_cube_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize, float radius);
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize);
void make_sphere_cube(std::vector<vec4i>& quadss, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize);
void make_sphere_flipcap(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize, const vec2f& zflip);
void make_disk(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize);
void make_disk_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize);
void make_disk_bulged(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize, float height);
void make_cylinder_side(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize, bool capped);
void make_cylinder(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize);
void make_cylinder_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius);
void make_sphere(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation);
void make_geodesic_sphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, int tesselation, float size);

// Make example watertight quad meshes for subdivision surfaces.
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);
void make_suzanne(std::vector<vec4i>& quads, std::vector<vec3f>& pos);

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius = {0.001f, 0.001f});

// Make point primitives
void make_point(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, float point_radius = 0.001f);
void make_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, float uvsize,
    float point_radius = 0.001f);
void make_random_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, const vec3f& size, float uvsize,
    float point_radius = 0.001f, uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
void make_bezier_circle(std::vector<vec4i>& beziers, std::vector<vec3f>& pos);

// Make a hair ball around a shape.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& tang, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2f& length = {0.1f, 0.1f}, const vec2f& rad = {0.001f, 0.001f},
    const vec2f& noise = zero2f, const vec2f& clump = zero2f,
    const vec2f& rotation = zero2f, int seed = 7);

}  // namespace ygl

#endif
