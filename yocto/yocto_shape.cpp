//
// Implementation for Yocto/Shape.
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

#include "yocto_shape.h"
#include "yocto_random.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Compute per-vertex tangents for lines.
vector<vec3f> compute_vertex_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
    auto tangents = vector<vec3f>{positions.size()};
    for (auto& tangent : tangents) tangent = zero3f;
    for (auto& l : lines) {
        auto tangent = line_tangent(positions[l.x], positions[l.y]);
        auto length  = line_length(positions[l.x], positions[l.y]);
        tangents[l.x] += tangent * length;
        tangents[l.y] += tangent * length;
    }
    for (auto& tangent : tangents) tangent = normalize(tangent);
    return tangents;
}

// Compute per-vertex normals for triangles.
vector<vec3f> compute_vertex_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
    auto normals = vector<vec3f>{positions.size()};
    for (auto& normal : normals) normal = zero3f;
    for (auto& t : triangles) {
        auto normal = triangle_normal(
            positions[t.x], positions[t.y], positions[t.z]);
        auto area = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
        normals[t.x] += normal * area;
        normals[t.y] += normal * area;
        normals[t.z] += normal * area;
    }
    for (auto& normal : normals) normal = normalize(normal);
    return normals;
}

// Compute per-vertex normals for quads.
vector<vec3f> compute_vertex_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
    auto normals = vector<vec3f>{positions.size()};
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
    return normals;
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texturecoords) {
    auto tangu = vector<vec3f>(positions.size(), zero3f);
    auto tangv = vector<vec3f>(positions.size(), zero3f);
    for (auto t : triangles) {
        auto tutv = triangle_tangents_fromuv(positions[t.x], positions[t.y],
            positions[t.z], texturecoords[t.x], texturecoords[t.y],
            texturecoords[t.z]);
        tutv      = {normalize(tutv.first), normalize(tutv.second)};
        for (auto vid : {t.x, t.y, t.z}) tangu[vid] += tutv.first;
        for (auto vid : {t.x, t.y, t.z}) tangv[vid] += tutv.second;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    auto tangentspaces = vector<vec4f>(positions.size(), zero4f);
    for (auto i = 0; i < positions.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], normals[i]);
        auto s = (dot(cross(normals[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangentspaces[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
    return tangentspaces;
}

// Apply skinning
tuple<vector<vec3f>, vector<vec3f>> compute_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms) {
    auto skinned_positions = vector<vec3f>{positions.size()};
    auto skinned_normals   = vector<vec3f>{normals};
    for (auto i = 0; i < positions.size(); i++) {
        skinned_positions[i] = transform_point(xforms[joints[i].x], positions[i]) *
                                   weights[i].x +
                               transform_point(xforms[joints[i].y], positions[i]) *
                                   weights[i].y +
                               transform_point(xforms[joints[i].z], positions[i]) *
                                   weights[i].z +
                               transform_point(xforms[joints[i].w], positions[i]) *
                                   weights[i].w;
    }
    for (auto i = 0; i < positions.size(); i++) {
        skinned_normals[i] = normalize(
            transform_direction(xforms[joints[i].x], normals[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], normals[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], normals[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], normals[i]) * weights[i].w);
    }
    return {skinned_positions, skinned_normals};
}

// Apply skinning as specified in Khronos glTF
tuple<vector<vec3f>, vector<vec3f>> compute_matrix_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms) {
    auto skinned_positions = vector<vec3f>{positions.size()};
    auto skinned_normals   = vector<vec3f>{normals};
    for (auto i = 0; i < positions.size(); i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_positions[i] = transform_point(xform, positions[i]);
        skinned_normals[i] = normalize(transform_direction(xform, normals[i]));
    }
    return {skinned_positions, skinned_normals};
}

// Creates an edge map
edge_map make_edge_map(const vector<vec3i>& triangles) {
    auto emap = edge_map{};
    for (int i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        insert_edge(emap, {t.x, t.y}, i);
        insert_edge(emap, {t.y, t.z}, i);
        insert_edge(emap, {t.z, t.x}, i);
    }
    return emap;
}
edge_map make_edge_map(const vector<vec4i>& quads) {
    auto emap = edge_map{};
    for (int i = 0; i < quads.size(); i++) {
        auto& q = quads[i];
        insert_edge(emap, {q.x, q.y}, i);
        insert_edge(emap, {q.y, q.z}, i);
        if (q.z != q.w) insert_edge(emap, {q.z, q.w}, i);
        insert_edge(emap, {q.w, q.x}, i);
    }
    return emap;
}

// Create key entry for edge_map
vec2i make_edge(const vec2i& e) { return e.x < e.y ? e : vec2i{e.y, e.x}; }

// Initialize an edge map with elements.
void insert_edges(edge_map& emap, const vector<vec3i>& triangles) {
    for (int i = 0; i < triangles.size(); i++) {
        auto& t = triangles[i];
        insert_edge(emap, {t.x, t.y}, i);
        insert_edge(emap, {t.y, t.z}, i);
        insert_edge(emap, {t.z, t.x}, i);
    }
}
void insert_edges(edge_map& emap, const vector<vec4i>& quads) {
    for (int i = 0; i < quads.size(); i++) {
        auto& q = quads[i];
        insert_edge(emap, {q.x, q.y}, i);
        insert_edge(emap, {q.y, q.z}, i);
        if (q.z != q.w) insert_edge(emap, {q.z, q.w}, i);
        insert_edge(emap, {q.w, q.x}, i);
    }
}
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& e, int face) {
    auto es = make_edge(e);
    auto it = emap.find(es);
    if (it == emap.end()) {
        auto idx = (int)emap.size();
        emap.insert(it, {es, {idx, face, -1}});
        return idx;
    } else {
        it->second.z = face;
        return it->second.x;
    }
}
// Get the edge index
int get_edge_index(const edge_map& emap, const vec2i& e) {
    auto es = make_edge(e);
    return emap.at(es).x;
}
// Get the edge count
int get_edge_count(const edge_map& emap, const vec2i& e) {
    auto es = make_edge(e);
    return emap.at(es).z == -1 ? 1 : 2;
}
// Get a list of edges, boundary edges, boundary vertices
vector<vec2i> get_edges(const edge_map& emap) {
    auto edges = vector<vec2i>(emap.size());
    for (auto& kv : emap) edges[kv.second.x] = kv.first;
    return edges;
}
vector<vec2i> get_boundary(const edge_map& emap) {
    auto boundary = vector<vec2i>();
    for (auto& kv : emap)
        if (kv.second.z == -1) boundary.push_back(kv.first);
    return boundary;
}

// Convert quads to triangles
vector<vec3i> convert_quads_to_triangles(const vector<vec4i>& quads) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
}

// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
vector<vec3i> convert_quads_to_triangles(
    const vector<vec4i>& quads, int row_length) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
#if 0
        triangles.resize(usteps * vsteps * 2);
        for (auto j = 0; j < vsteps; j++) {
            for (auto i = 0; i < usteps; i++) {
                auto f1 = triangles[(j * usteps + i) * 2 + 0];
                auto f2 = triangles[(j * usteps + i) * 2 + 1];
                if ((i + j) % 2) {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
                } else {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
                }
            }
        }
#endif
    return triangles;
}

// Convert triangles to quads by creating degenerate quads
vector<vec4i> convert_triangles_to_quads(const vector<vec3i>& triangles) {
    auto quads = vector<vec4i>();
    quads.reserve(triangles.size());
    for (auto& t : triangles) quads.push_back({t.x, t.y, t.z, t.z});
    return quads;
}

// Convert beziers to lines using 3 lines for each bezier.
vector<vec2i> convert_bezier_to_lines(const vector<vec4i>& beziers) {
    auto lines = vector<vec2i>();
    lines.reserve(beziers.size() * 3);
    for (auto b : beziers) {
        lines.push_back({b.x, b.y});
        lines.push_back({b.y, b.z});
        lines.push_back({b.z, b.w});
    }
    return lines;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> convert_face_varying(
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords) {
    // make faces unique
    unordered_map<vec3i, int> vert_map;
    auto split_quads = vector<vec4i>(quads_positions.size());
    for (auto fid = 0; fid < quads_positions.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec3i{
                (&quads_positions[fid].x)[c],
                (!quads_normals.empty()) ? (&quads_normals[fid].x)[c] : -1,
                (!quads_texturecoords.empty()) ? (&quads_texturecoords[fid].x)[c] :
                                                 -1,
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
    auto split_positions = vector<vec3f>{};
    if (!positions.empty()) {
        split_positions.resize(vert_map.size());
        for (auto kv : vert_map) {
            split_positions[kv.second] = positions[kv.first.x];
        }
    }
    auto split_normals = vector<vec3f>{};
    if (!normals.empty()) {
        split_normals.resize(vert_map.size());
        for (auto kv : vert_map) {
            split_normals[kv.second] = normals[kv.first.y];
        }
    }
    auto split_texcturecoords = vector<vec2f>{};
    if (!texturecoords.empty()) {
        split_texcturecoords.resize(vert_map.size());
        for (auto kv : vert_map) {
            split_texcturecoords[kv.second] = texturecoords[kv.first.z];
        }
    }

    return {split_quads, split_positions, split_normals, split_texcturecoords};
}

// Split primitives per id
template <typename T>
vector<vector<T>> ungroup_elems(const vector<T>& elems, const vector<int>& ids) {
    auto max_id      = *std::max_element(ids.begin(), ids.end());
    auto split_elems = vector<vector<T>>(max_id + 1);
    for (auto elem_id = 0; elem_id < elems.size(); elem_id++) {
        split_elems[ids[elem_id]].push_back(elems[elem_id]);
    }
    return split_elems;
}
vector<vector<vec2i>> ungroup_lines(
    const vector<vec2i>& lines, const vector<int>& ids) {
    return ungroup_elems(lines, ids);
}
vector<vector<vec3i>> ungroup_triangles(
    const vector<vec3i>& triangles, const vector<int>& ids) {
    return ungroup_elems(triangles, ids);
}
vector<vector<vec4i>> ungroup_quads(
    const vector<vec4i>& quads, const vector<int>& ids) {
    return ungroup_elems(quads, ids);
}

// Subdivide lines.
template <typename T>
tuple<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vert) {
    // early exit
    if (lines.empty() || vert.empty()) return {lines, vert};
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
    return {tlines, tvert};
}

template tuple<vector<vec2i>, vector<float>> subdivide_lines(
    const vector<vec2i>&, const vector<float>&);
template tuple<vector<vec2i>, vector<vec2f>> subdivide_lines(
    const vector<vec2i>&, const vector<vec2f>&);
template tuple<vector<vec2i>, vector<vec3f>> subdivide_lines(
    const vector<vec2i>&, const vector<vec3f>&);
template tuple<vector<vec2i>, vector<vec4f>> subdivide_lines(
    const vector<vec2i>&, const vector<vec4f>&);

// Subdivide triangle.
template <typename T>
tuple<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vert) {
    // early exit
    if (triangles.empty() || vert.empty()) return {triangles, vert};
    // get edges
    auto emap  = make_edge_map(triangles);
    auto edges = get_edges(emap);
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
    return {ttriangles, tvert};
}

template tuple<vector<vec3i>, vector<float>> subdivide_triangles(
    const vector<vec3i>&, const vector<float>&);
template tuple<vector<vec3i>, vector<vec2f>> subdivide_triangles(
    const vector<vec3i>&, const vector<vec2f>&);
template tuple<vector<vec3i>, vector<vec3f>> subdivide_triangles(
    const vector<vec3i>&, const vector<vec3f>&);
template tuple<vector<vec3i>, vector<vec4f>> subdivide_triangles(
    const vector<vec3i>&, const vector<vec4f>&);

// Subdivide quads.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vert) {
    // early exit
    if (quads.empty() || vert.empty()) return {quads, vert};
    // get edges
    auto emap  = make_edge_map(quads);
    auto edges = get_edges(emap);
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
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.z] +
                                             vert[q.w]) /
                                         4;
        } else {
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.y]) / 3;
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
    return {tquads, tvert};
}

template tuple<vector<vec4i>, vector<float>> subdivide_quads(
    const vector<vec4i>&, const vector<float>&);
template tuple<vector<vec4i>, vector<vec2f>> subdivide_quads(
    const vector<vec4i>&, const vector<vec2f>&);
template tuple<vector<vec4i>, vector<vec3f>> subdivide_quads(
    const vector<vec4i>&, const vector<vec3f>&);
template tuple<vector<vec4i>, vector<vec4f>> subdivide_quads(
    const vector<vec4i>&, const vector<vec4f>&);

// Subdivide beziers.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vert) {
    // early exit
    if (beziers.empty() || vert.empty()) return {beziers, vert};
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
    return {tbeziers, tvert};
}

template tuple<vector<vec4i>, vector<float>> subdivide_beziers(
    const vector<vec4i>&, const vector<float>&);
template tuple<vector<vec4i>, vector<vec2f>> subdivide_beziers(
    const vector<vec4i>&, const vector<vec2f>&);
template tuple<vector<vec4i>, vector<vec3f>> subdivide_beziers(
    const vector<vec4i>&, const vector<vec3f>&);
template tuple<vector<vec4i>, vector<vec4f>> subdivide_beziers(
    const vector<vec4i>&, const vector<vec4f>&);

// Subdivide catmullclark.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_catmullclark(
    const vector<vec4i>& quads, const vector<T>& vert, bool lock_boundary) {
    // early exit
    if (quads.empty() || vert.empty()) return {quads, vert};
    // get edges
    auto emap     = make_edge_map(quads);
    auto edges    = get_edges(emap);
    auto boundary = get_boundary(emap);
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
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.z] +
                                             vert[q.w]) /
                                         4;
        } else {
            tvert[nverts + nedges + i] = (vert[q.x] + vert[q.y] + vert[q.y]) / 3;
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
    for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i = 0; i < tvert.size(); i++) {
        if (tvert_val[i] != 2) continue;
        avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4.0f / acount[i]);
    }
    tvert = avert;

    // done
    return {tquads, tvert};
}

template tuple<vector<vec4i>, vector<float>> subdivide_catmullclark(
    const vector<vec4i>&, const vector<float>&, bool);
template tuple<vector<vec4i>, vector<vec2f>> subdivide_catmullclark(
    const vector<vec4i>&, const vector<vec2f>&, bool);
template tuple<vector<vec4i>, vector<vec3f>> subdivide_catmullclark(
    const vector<vec4i>&, const vector<vec3f>&, bool);
template tuple<vector<vec4i>, vector<vec4f>> subdivide_catmullclark(
    const vector<vec4i>&, const vector<vec4f>&, bool);

// Weld vertices within a threshold. For noe the implementation is O(n^2).
tuple<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold) {
    auto welded_indices   = vector<int>(positions.size());
    auto welded_positions = vector<vec3f>();
    for (auto i = 0; i < positions.size(); i++) {
        welded_indices[i] = (int)welded_positions.size();
        for (auto j = 0; j < welded_positions.size(); j++) {
            if (length(positions[i] - welded_positions[j]) < threshold) {
                welded_indices[i] = j;
                break;
            }
        }
        if (welded_indices[i] == (int)welded_positions.size())
            welded_positions.push_back(positions[i]);
    }
    return {welded_positions, welded_indices};
}
void weld_vertices(
    vector<vec3f>& positions, float threshold, vector<int>& welded_indices) {
    tie(positions, welded_indices) = weld_vertices(positions, threshold);
}
tuple<vector<vec3i>, vector<vec3f>> weld_triangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, float threshold) {
    auto welded_positions                 = vector<vec3f>{};
    auto welded_indices                   = vector<int>{};
    tie(welded_positions, welded_indices) = weld_vertices(positions, threshold);
    auto welded_triangles                 = vector<vec3i>{};
    for (auto& t : triangles) {
        welded_triangles.push_back(
            {welded_indices[t.x], welded_indices[t.y], welded_indices[t.z]});
    }
    return {welded_triangles, welded_positions};
}
tuple<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold) {
    auto welded_positions                 = vector<vec3f>{};
    auto welded_indices                   = vector<int>{};
    tie(welded_positions, welded_indices) = weld_vertices(positions, threshold);
    auto welded_quads                     = vector<vec4i>{};
    for (auto& q : quads) {
        welded_quads.push_back({
            welded_indices[q.x],
            welded_indices[q.y],
            welded_indices[q.z],
            welded_indices[q.w],
        });
    }
    return {welded_quads, welded_positions};
}

// Pick a point in a point set uniformly.
int sample_points_element(int npoints, float re) {
    return sample_uniform_index(npoints, re);
}
vector<float> sample_points_element_cdf(int npoints) {
    auto cdf = vector<float>(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
    return cdf;
}
int sample_points_element(const vector<float>& cdf, float re) {
    return sample_discrete_distribution(cdf, re);
}

// Pick a point on lines uniformly.
vector<float> sample_lines_element_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
    auto cdf = vector<float>(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(positions[l.x], positions[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
pair<int, float> sample_lines_element(
    const vector<float>& cdf, float re, float ru) {
    return {sample_discrete_distribution(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
vector<float> sample_triangles_element_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
    auto cdf = vector<float>(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
pair<int, vec2f> sample_triangles_element(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete_distribution(cdf, re),
        sample_triangle_coordinates(ruv)};
}

// Pick a point on a quad mesh uniformly.
vector<float> sample_quads_element_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
    auto cdf = vector<float>(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
pair<int, vec2f> sample_quads_element(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete_distribution(cdf, re), ruv};
}
pair<int, vec2f> sample_quads_element(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv) {
    auto element_id = sample_discrete_distribution(cdf, re);
    if (quads[element_id].z == quads[element_id].w) {
        return {element_id, sample_triangle_coordinates(ruv)};
    } else {
        return {element_id, ruv};
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_triangles_points(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    int npoints, int seed) {
    auto sampled_positions     = vector<vec3f>(npoints);
    auto sampled_normals       = vector<vec3f>(npoints);
    auto sampled_texturecoords = vector<vec2f>(npoints);
    auto cdf = sample_triangles_element_cdf(triangles, positions);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto sample = sample_triangles_element(cdf, get_random_float(rng),
            {get_random_float(rng), get_random_float(rng)});
        auto t      = triangles[sample.first];
        sampled_positions[i] = interpolate_triangle(
            positions[t.x], positions[t.y], positions[t.z], sample.second);
        if (!sampled_normals.empty()) {
            sampled_normals[i] = normalize(interpolate_triangle(
                normals[t.x], normals[t.y], normals[t.z], sample.second));
        } else {
            sampled_normals[i] = triangle_normal(
                positions[t.x], positions[t.y], positions[t.z]);
        }
        if (!sampled_texturecoords.empty()) {
            sampled_texturecoords[i] = interpolate_triangle(texturecoords[t.x],
                texturecoords[t.y], texturecoords[t.z], sample.second);
        } else {
            sampled_texturecoords[i] = zero2f;
        }
    }
    return {sampled_positions, sampled_normals, sampled_texturecoords};
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_quads_points(
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    int npoints, int seed) {
    auto sampled_positions     = vector<vec3f>(npoints);
    auto sampled_normals       = vector<vec3f>(npoints);
    auto sampled_texturecoords = vector<vec2f>(npoints);
    auto cdf                   = sample_quads_element_cdf(quads, positions);
    auto rng                   = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto sample          = sample_quads_element(cdf, get_random_float(rng),
            {get_random_float(rng), get_random_float(rng)});
        auto q               = quads[sample.first];
        sampled_positions[i] = interpolate_quad(positions[q.x], positions[q.y],
            positions[q.z], positions[q.w], sample.second);
        if (!sampled_normals.empty()) {
            sampled_normals[i] = normalize(interpolate_quad(normals[q.x],
                normals[q.y], normals[q.z], normals[q.w], sample.second));
        } else {
            sampled_normals[i] = quad_normal(
                positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        }
        if (!sampled_texturecoords.empty()) {
            sampled_texturecoords[i] = interpolate_quad(texturecoords[q.x],
                texturecoords[q.y], texturecoords[q.z], texturecoords[q.w],
                sample.second);
        } else {
            sampled_texturecoords[i] = zero2f;
        }
    }
    return {sampled_positions, sampled_normals, sampled_texturecoords};
}

// Merge shape elements
void merge_lines(
    vector<vec2i>& lines, const vector<vec2i>& merge_lines, int num_verts) {
    for (auto& l : merge_lines)
        lines.push_back({l.x + num_verts, l.y + num_verts});
}
void merge_triangles(vector<vec3i>& triangles,
    const vector<vec3i>& merge_triangles, int num_verts) {
    for (auto& t : merge_triangles)
        triangles.push_back({t.x + num_verts, t.y + num_verts, t.z + num_verts});
}
void merge_quads(
    vector<vec4i>& quads, const vector<vec4i>& merge_quads, int num_verts) {
    for (auto& q : merge_quads)
        quads.push_back({q.x + num_verts, q.y + num_verts, q.z + num_verts,
            q.w + num_verts});
}
void merge_lines(vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& tangents, vector<vec2f>& texturecoords,
    vector<float>& radius, const vector<vec2i>& merge_lines,
    const vector<vec3f>& merge_positions, const vector<vec3f>& merge_tangents,
    const vector<vec2f>& merge_texturecoords, const vector<float>& merge_radius) {
    auto merge_verts = (int)positions.size();
    for (auto& l : merge_lines)
        lines.push_back({l.x + merge_verts, l.y + merge_verts});
    positions.insert(
        positions.end(), merge_positions.begin(), merge_positions.end());
    tangents.insert(tangents.end(), merge_tangents.begin(), merge_tangents.end());
    texturecoords.insert(texturecoords.end(), merge_texturecoords.begin(),
        merge_texturecoords.end());
    radius.insert(radius.end(), merge_radius.begin(), merge_radius.end());
}
void merge_triangles(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    const vector<vec3i>& merge_triangles, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals, const vector<vec2f>& merge_texturecoords) {
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
void merge_quads(vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    const vector<vec4i>& merge_quads, const vector<vec3f>& merge_positions,
    const vector<vec3f>& merge_normals, const vector<vec2f>& merge_texturecoords) {
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
// IMPLEMENTATION OF SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace yocto {

// Make a quad.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_quad_shape(
    const vec2i& steps, const vec2f& size, const vec2f& uvsize) {
    auto positions     = vector<vec3f>((steps.x + 1) * (steps.y + 1));
    auto normals       = vector<vec3f>((steps.x + 1) * (steps.y + 1));
    auto texturecoords = vector<vec2f>((steps.x + 1) * (steps.y + 1));
    for (auto j = 0; j <= steps.y; j++) {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
            positions[j * (steps.x + 1) + i] = {
                (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
            normals[j * (steps.x + 1) + i]       = {0, 0, 1};
            texturecoords[j * (steps.x + 1) + i] = uv * uvsize;
        }
    }

    auto quads = vector<vec4i>(steps.x * steps.y);
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
    return {quads, positions, normals, texturecoords};
}

tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_floor_shape(
    const vec2i& steps, const vec2f& size, const vec2f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_quad_shape(
        steps, size, uvsize);
    for (auto& p : positions) p = {p.x, p.z, p.y};
    for (auto& normal : normals) normal = {normal.x, normal.z, normal.y};
    return {quads, positions, normals, texturecoords};
}

// Make a stack of quads
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_quad_stack_shape(
    const vec3i& steps, const vec3f& size, const vec2f& uvsize) {
    auto quads          = vector<vec4i>{};
    auto positions      = vector<vec3f>{};
    auto normals        = vector<vec3f>{};
    auto texturecoords  = vector<vec2f>{};
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec3f>{};
    auto qnormals       = vector<vec3f>{};
    auto qtexturecoords = vector<vec2f>{};
    for (auto i = 0; i <= steps.z; i++) {
        tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
            {steps.x, steps.y}, {size.x, size.y}, uvsize);
        for (auto& p : qpositions) p.z = (-0.5f + (float)i / steps.z) * size.z;
        merge_quads(quads, positions, normals, texturecoords, qquads,
            qpositions, qnormals, qtexturecoords);
    }
    return {quads, positions, normals, texturecoords};
}

// Make a cube.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_cube_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize) {
    auto quads          = vector<vec4i>{};
    auto positions      = vector<vec3f>{};
    auto normals        = vector<vec3f>{};
    auto texturecoords  = vector<vec2f>{};
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec3f>{};
    auto qnormals       = vector<vec3f>{};
    auto qtexturecoords = vector<vec2f>{};
    // + z
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {qpositions[i].x, qpositions[i].y, size.z / 2};
        qnormals[i]   = {0, 0, 1};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - z
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {-qpositions[i].x, qpositions[i].y, -size.z / 2};
        qnormals[i]   = {0, 0, -1};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // + x
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.y, steps.z}, {size.y, size.z}, {uvsize.y, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {size.x / 2, qpositions[i].y, -qpositions[i].x};
        qnormals[i]   = {1, 0, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - x
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.y, steps.z}, {size.y, size.z}, {uvsize.y, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {-size.x / 2, qpositions[i].y, qpositions[i].x};
        qnormals[i]   = {-1, 0, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // + y
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {qpositions[i].x, size.y / 2, -qpositions[i].y};
        qnormals[i]   = {0, 1, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // - y
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_quad_shape(
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i] = {qpositions[i].x, -size.y / 2, qpositions[i].y};
        qnormals[i]   = {0, -1, 0};
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    return {quads, positions, normals, texturecoords};
}

// Make a rounded cube.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_cube_rounded_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize, float radius) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_cube_shape(
        steps, size, uvsize);
    auto c = size / 2 - vec3f{radius, radius, radius};
    for (auto i = 0; i < positions.size(); i++) {
        auto pc = vec3f{
            fabs(positions[i].x), fabs(positions[i].y), fabs(positions[i].z)};
        auto ps = vec3f{positions[i].x < 0 ? -1.0f : 1.0f,
            positions[i].y < 0 ? -1.0f : 1.0f, positions[i].z < 0 ? -1.0f : 1.0f};
        if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
            auto pn      = normalize(pc - c);
            positions[i] = c + radius * pn;
            normals[i]   = pn;
        } else if (pc.x >= c.x && pc.y >= c.y) {
            auto pn      = normalize((pc - c) * vec3f{1, 1, 0});
            positions[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
            normals[i]   = pn;
        } else if (pc.x >= c.x && pc.z >= c.z) {
            auto pn      = normalize((pc - c) * vec3f{1, 0, 1});
            positions[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
            normals[i]   = pn;
        } else if (pc.y >= c.y && pc.z >= c.z) {
            auto pn      = normalize((pc - c) * vec3f{0, 1, 1});
            positions[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
            normals[i]   = pn;
        } else {
            continue;
        }
        positions[i] *= ps;
        normals[i] *= ps;
    }
    return {quads, positions, normals, texturecoords};
}

// Make a sphere.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_sphere_shape(
    const vec2i& steps, float size, const vec2f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_quad_shape(
        steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv = texturecoords[i];
        auto a  = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        auto p  = vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        positions[i]     = p * (size / 2);
        normals[i]       = normalize(p);
        texturecoords[i] = uv * uvsize;
    }
    return {quads, positions, normals, texturecoords};
}

// Make a spherecube.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_sphere_cube_shape(
    int steps, float size, float uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_cube_shape(
        {steps, steps, steps}, {1, 1, 1}, {uvsize, uvsize, uvsize});
    for (auto i = 0; i < positions.size(); i++) {
        auto p       = positions[i];
        positions[i] = normalize(p) * (size / 2);
        normals[i]   = normalize(p);
    }
    return {quads, positions, normals, texturecoords};
}

// Make a flipped sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_sphere_flipcap_shape(
    const vec2i& steps, float size, const vec2f& uvsize, const vec2f& zflip) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_sphere_shape(
        steps, size, uvsize);
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
    return {quads, positions, normals, texturecoords};
}

// Make a disk.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_disk_shape(
    const vec2i& steps, float size, const vec2f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_quad_shape(
        steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv      = texturecoords[i];
        auto phi     = 2 * pif * uv.x;
        positions[i] = {
            cos(phi) * uv.y * size / 2, sin(phi) * uv.y * size / 2, 0};
        normals[i]       = {0, 0, 1};
        texturecoords[i] = uv * uvsize;
    }
    return {quads, positions, normals, texturecoords};
}

// Make a disk from a quad.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_disk_quad_shape(
    int steps, float size, float uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_quad_shape(
        {steps, steps}, {2, 2}, {uvsize, uvsize});
    for (auto i = 0; i < positions.size(); i++) {
        // Analytical Methods for Squaring the Disc, by C. Fong
        // https://arxiv.org/abs/1509.06344
        auto xy = vec2f{positions[i].x, positions[i].y};
        auto uv = vec2f{
            xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
        positions[i] = {uv.x * size / 2, uv.y * size / 2, 0};
    }
    return {quads, positions, normals, texturecoords};
}

// Make a bulged disk from a quad.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_disk_bulged_shape(
    int steps, float size, float uvsize, float height) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_disk_quad_shape(
        steps, size, uvsize);
    if (height == 0) return {quads, positions, normals, texturecoords};
    auto radius = (size * size / 4 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i = 0; i < positions.size(); i++) {
        auto pn      = normalize(positions[i] - center);
        positions[i] = center + pn * radius;
        normals[i]   = pn;
    }
    return {quads, positions, normals, texturecoords};
}

// Make a cylinder (side-only).
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_cylinder_side_shape(
    const vec2i& steps, const vec2f& size, const vec2f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_quad_shape(
        steps, {1, 1}, {1, 1});
    for (auto i = 0; i < positions.size(); i++) {
        auto uv          = texturecoords[i];
        auto phi         = 2 * pif * uv.x;
        positions[i]     = {cos(phi) * size.x / 2, sin(phi) * size.x / 2,
            (uv.y - 0.5f) * size.y};
        normals[i]       = {cos(phi), sin(phi), 0};
        texturecoords[i] = uv * uvsize;
    }
    return {quads, positions, normals, texturecoords};
}

// Make a cylinder.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_cylinder_shape(
    const vec3i& steps, const vec2f& size, const vec3f& uvsize) {
    auto quads          = vector<vec4i>{};
    auto positions      = vector<vec3f>{};
    auto normals        = vector<vec3f>{};
    auto texturecoords  = vector<vec2f>{};
    auto qquads         = vector<vec4i>{};
    auto qpositions     = vector<vec3f>{};
    auto qnormals       = vector<vec3f>{};
    auto qtexturecoords = vector<vec2f>{};
    // side
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_cylinder_side_shape(
        {steps.x, steps.y}, {size.x, size.y}, {uvsize.x, uvsize.y});
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // top
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_disk_shape(
        {steps.x, steps.z}, size.x, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i].z = size.y / 2;
    }
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    // bottom
    tie(qquads, qpositions, qnormals, qtexturecoords) = make_disk_shape(
        {steps.x, steps.z}, size.x, {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpositions.size(); i++) {
        qpositions[i].z = -size.y / 2;
        qnormals[i]     = -qnormals[i];
    }
    for (auto i = 0; i < qquads.size(); i++) swap(qquads[i].x, qquads[i].z);
    merge_quads(quads, positions, normals, texturecoords, qquads, qpositions,
        qnormals, qtexturecoords);
    return {quads, positions, normals, texturecoords};
}

// Make a rounded cylinder.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_cylinder_rounded_shape(
    const vec3i& steps, const vec2f& size, const vec3f& uvsize, float radius) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_cylinder_shape(
        steps, size, uvsize);
    auto c = size / 2 - vec2f{radius, radius};
    for (auto i = 0; i < positions.size(); i++) {
        auto phi = atan2(positions[i].y, positions[i].x);
        auto r   = length(vec2f{positions[i].x, positions[i].y});
        auto z   = positions[i].z;
        auto pc  = vec2f{r, fabs(z)};
        auto ps  = (z < 0) ? -1.0f : 1.0f;
        if (pc.x >= c.x && pc.y >= c.y) {
            auto pn      = normalize(pc - c);
            positions[i] = {cos(phi) * c.x + radius * pn.x,
                sin(phi) * c.x + radius * pn.x, ps * (c.y + radius * pn.y)};
            normals[i]   = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
        } else {
            continue;
        }
    }
    return {quads, positions, normals, texturecoords};
}

// Make a geodesic sphere.
tuple<vector<vec3i>, vector<vec3f>, vector<vec3f>> make_geodesic_sphere_shape(
    int tesselation, float size) {
    // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
    const float X                = 0.525731112119133606f;
    const float Z                = 0.850650808352039932f;
    static auto sphere_pos       = vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z},
        {-X, 0.0, -Z}, {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X},
        {0.0, -Z, -X}, {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
    static auto sphere_triangles = vector<vec3i>{{0, 1, 4}, {0, 4, 9},
        {9, 4, 5}, {4, 8, 5}, {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3},
        {5, 3, 2}, {2, 3, 7}, {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0},
        {0, 6, 1}, {6, 10, 1}, {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
    auto        positions        = sphere_pos;
    auto        triangles        = sphere_triangles;
    for (auto l = 0; l < max(0, tesselation - 2); l++) {
        tie(triangles, positions) = subdivide_triangles(triangles, positions);
    }
    for (auto& p : positions) p = normalize(p) * size / 2;
    auto normals = positions;
    return {triangles, positions, normals};
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
tuple<vector<vec4i>, vector<vec4i>, vector<vec4i>, vector<vec3f>, vector<vec3f>,
    vector<vec2f>>
make_cube_facevarying_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_cube_shape(
        steps, size, uvsize);
    auto quads_positions            = quads;
    auto quads_normals              = quads;
    auto quads_texturecoords        = quads;
    tie(quads_positions, positions) = weld_quads(quads_positions, positions,
        min(0.1f * size / vec3f{(float)steps.x, (float)steps.y, (float)steps.z}));
    return {quads_positions, quads_normals, quads_texturecoords, positions,
        normals, texturecoords};
}
tuple<vector<vec4i>, vector<vec4i>, vector<vec4i>, vector<int>, vector<vec3f>,
    vector<vec3f>, vector<vec2f>>
make_cube_multiplematerials_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize) {
    auto quads_positions     = vector<vec4i>{};
    auto quads_normals       = vector<vec4i>{};
    auto quads_texturecoords = vector<vec4i>{};
    auto positions           = vector<vec3f>{};
    auto normals             = vector<vec3f>{};
    auto texturecoords       = vector<vec2f>{};
    tie(quads_positions, quads_normals, quads_texturecoords, positions, normals,
        texturecoords)       = make_cube_facevarying_shape(steps, size, uvsize);
    auto quads_materials     = vector<int>(quads_positions.size());
    auto quads_per_face      = (int)quads_positions.size() / 6;
    for (auto i = 0; i < quads_positions.size(); i++) {
        quads_materials[i] = i / quads_per_face;
    }
    return {quads_positions, quads_normals, quads_texturecoords,
        quads_materials, positions, normals, texturecoords};
}
tuple<vector<vec4i>, vector<vec3f>> make_cube_posonly_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize) {
    auto quads                                    = vector<vec4i>{};
    auto positions                                = vector<vec3f>{};
    auto normals                                  = vector<vec3f>{};
    auto texturecoords                            = vector<vec2f>{};
    tie(quads, positions, normals, texturecoords) = make_cube_shape(
        steps, size, uvsize);
    tie(quads, positions) = weld_quads(quads, positions,
        min(0.1f * size / vec3f{(float)steps.x, (float)steps.y, (float)steps.z}));
    return {quads, positions};
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
tuple<vector<vec4i>, vector<vec3f>> make_suzanne_shape(float size) {
    static auto suzanne_pos       = vector<vec3f>{{0.4375, 0.1640625, 0.765625},
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
        {-0.9453125, 0.3046875, -0.2890625}, {0.8828125, -0.0234375, -0.2109375},
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

    auto positions = suzanne_pos;
    for (auto& p : positions) p *= size / 2;
    auto quads = suzanne_quads;
    for (auto& t : suzanne_triangles) {
        quads.push_back({t.x, t.y, t.z, t.z});
    }
    return {quads, positions};
}

// Watertight cube
tuple<vector<vec4i>, vector<vec3f>> make_cube_shape(const vec3f& size) {
    static auto cube_pos = vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1},
        {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    static auto cube_quads   = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_quad_uv = vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    auto        positions    = cube_pos;
    for (auto& p : positions) p *= size / 2;
    auto quads = cube_quads;
    return {quads, positions};
}

// Generate lines set along a quad.
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>> make_lines_shape(
    const vec2i& steps, const vec2f& size, const vec2f& uvsize,
    const vec2f& line_radius) {
    auto nverts = (steps.x + 1) * steps.y;
    auto nlines = steps.x * steps.y;
    auto vid    = [steps](int i, int j) { return j * (steps.x + 1) + i; };
    auto fid    = [steps](int i, int j) { return j * steps.x + i; };

    auto positions     = vector<vec3f>(nverts);
    auto normals       = vector<vec3f>(nverts);
    auto texturecoords = vector<vec2f>(nverts);
    auto radius        = vector<float>(nverts);
    if (steps.y > 1) {
        for (auto j = 0; j < steps.y; j++) {
            for (auto i = 0; i <= steps.x; i++) {
                auto uv              = vec2f{i / (float)steps.x,
                    j / (float)(steps.y > 1 ? steps.y - 1 : 1)};
                positions[vid(i, j)] = {
                    (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
                normals[vid(i, j)]       = {1, 0, 0};
                texturecoords[vid(i, j)] = uv * uvsize;
            }
        }
    } else {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv                  = vec2f{i / (float)steps.x, 0};
            positions[vid(i, 0)]     = {(uv.x - 0.5f) * size.x, 0, 0};
            normals[vid(i, 0)]       = {1, 0, 0};
            texturecoords[vid(i, 0)] = uv * uvsize;
        }
    }

    auto lines = vector<vec2i>(nlines);
    for (int j = 0; j < steps.y; j++) {
        for (int i = 0; i < steps.x; i++) {
            lines[fid(i, j)] = {vid(i, j), vid(i + 1, j)};
        }
    }

    return {lines, positions, normals, texturecoords, radius};
}

// Generate a point set with points placed at the origin with texcoords
// varying along u.
tuple<vector<int>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>> make_points_shape(
    int num, float uvsize, float point_radius) {
    auto points = vector<int>(num);
    for (auto i = 0; i < num; i++) points[i] = i;
    auto positions     = vector<vec3f>(num, {0, 0, 0});
    auto normals       = vector<vec3f>(num, {0, 0, 1});
    auto texturecoords = vector<vec2f>(num, {0, 0});
    auto radius        = vector<float>(num, point_radius);
    for (auto i = 0; i < texturecoords.size(); i++)
        texturecoords[i] = {(float)i / (float)num, 0};
    return {points, positions, normals, texturecoords, radius};
}

// Generate a point set.
tuple<vector<int>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>> make_random_points_shape(
    int num, const vec3f& size, float uvsize, float point_radius, uint64_t seed) {
    auto points                                            = vector<int>();
    auto positions                                         = vector<vec3f>();
    auto normals                                           = vector<vec3f>();
    auto texturecoords                                     = vector<vec2f>();
    auto radius                                            = vector<float>();
    tie(points, positions, normals, texturecoords, radius) = make_points_shape(
        num, uvsize, point_radius);
    auto rng = make_rng(seed);
    for (auto i = 0; i < positions.size(); i++) {
        positions[i] = (get_random_vec3f(rng) - vec3f{0.5f, 0.5f, 0.5f}) * size;
    }
    return {points, positions, normals, texturecoords, radius};
}

// Make a point.
tuple<vector<int>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>> make_point_shape(
    float point_radius) {
    auto points        = vector<int>{0};
    auto positions     = vector<vec3f>{{0, 0, 0}};
    auto normals       = vector<vec3f>{{0, 0, 1}};
    auto texturecoords = vector<vec2f>{{0, 0}};
    auto radius        = vector<float>{point_radius};
    return {points, positions, normals, texturecoords, radius};
}

// Make a bezier circle. Returns bezier, pos.
tuple<vector<vec4i>, vector<vec3f>> make_bezier_circle_shape(float size) {
    // constant from http://spencermortensen.com/articles/bezier-circle/
    const auto  c              = 0.551915024494f;
    static auto circle_pos     = vector<vec3f>{{1, 0, 0}, {1, c, 0}, {c, 1, 0},
        {0, 1, 0}, {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
        {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
    static auto circle_beziers = vector<vec4i>{
        {0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
    auto positions = circle_pos;
    auto beziers   = circle_beziers;
    for (auto& p : positions) p *= size;
    return {beziers, positions};
}

// Make a hair ball around a shape
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>> make_hair_shape(
    const vec2i& steps, const vector<vec3i>& striangles,
    const vector<vec4i>& squads, const vector<vec3f>& spos,
    const vector<vec3f>& snorm, const vector<vec2f>& stexcoord,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
    vector<vec3f> bpos;
    vector<vec3f> bnorm;
    vector<vec2f> btexcoord;
    auto          alltriangles    = striangles;
    auto          quads_triangles = convert_quads_to_triangles(squads);
    alltriangles.insert(
        alltriangles.end(), quads_triangles.begin(), quads_triangles.end());
    tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        alltriangles, spos, snorm, stexcoord, steps.y, seed);

    auto rng  = make_rng(seed, 3);
    auto blen = vector<float>(bpos.size());
    for (auto& l : blen) l = lerp(len.x, len.y, get_random_float(rng));

    auto cidx = vector<int>();
    if (clump.x > 0) {
        for (auto bidx = 0; bidx < bpos.size(); bidx++) {
            cidx.push_back(0);
            auto cdist = maxf;
            for (auto c = 0; c < clump.y; c++) {
                auto d = length(bpos[bidx] - bpos[c]);
                if (d < cdist) {
                    cdist       = d;
                    cidx.back() = c;
                }
            }
        }
    }

    auto lines                                            = vector<vec2i>{};
    auto positions                                        = vector<vec3f>{};
    auto normals                                          = vector<vec3f>{};
    auto texturecoords                                    = vector<vec2f>{};
    auto radius                                           = vector<float>{};
    tie(lines, positions, normals, texturecoords, radius) = make_lines_shape(
        steps, {1, 1}, {1, 1});
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
            auto nx = perlin_noise(positions[i] * noise.y + vec3f{0, 0, 0}) *
                      noise.x;
            auto ny = perlin_noise(positions[i] * noise.y + vec3f{3, 7, 11}) *
                      noise.x;
            auto nz = perlin_noise(positions[i] * noise.y + vec3f{13, 17, 19}) *
                      noise.x;
            positions[i] += {nx, ny, nz};
        }
    }

    if (clump.x > 0 || noise.x > 0 || rotation.x > 0)
        normals = compute_vertex_tangents(lines, positions);

    return {lines, positions, normals, texturecoords, radius};
}

}  // namespace yocto
