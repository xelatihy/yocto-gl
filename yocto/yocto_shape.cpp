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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex tangents for lines.
void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto l : lines) {
        auto n = pos[l.y] - pos[l.x];
        if (!weighted) n = normalize(n);
        norm[l.x] += n;
        norm[l.y] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for triangles.
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto t : triangles) {
        auto n = cross(pos[t.y] - pos[t.x], pos[t.z] - pos[t.x]);
        if (!weighted) n = normalize(n);
        for (auto vid : {t.x, t.y, t.z}) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for quads.
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto q : quads) {
        auto n = cross(pos[q.y] - pos[q.x], pos[q.w] - pos[q.x]) +
                 cross(pos[q.w] - pos[q.z], pos[q.x] - pos[q.z]);
        if (!weighted) n = normalize(n);
        for (auto vid : {q.x, q.y, q.z, q.w}) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_frames(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted) {
    auto tangu = std::vector<vec3f>(pos.size(), zero3f);
    auto tangv = std::vector<vec3f>(pos.size(), zero3f);
    for (auto t : triangles) {
        auto tutv = triangle_tangents_fromuv(pos[t.x], pos[t.y], pos[t.z],
            texcoord[t.x], texcoord[t.y], texcoord[t.z]);
        if (!weighted) tutv = {normalize(tutv.first), normalize(tutv.second)};
        for (auto vid : {t.x, t.y, t.z}) tangu[vid] += tutv.first;
        for (auto vid : {t.x, t.y, t.z}) tangv[vid] += tutv.second;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    tangsp.resize(pos.size());
    for (auto& t : tangsp) t = zero4f;
    for (auto i = 0; i < pos.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], norm[i]);
        auto s = (dot(cross(norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangsp[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
}

// Apply skinning
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

// Apply skinning
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

// Apply skinning as specified in Khronos glTF
void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_pos[i] = transform_point(xform, pos[i]);
        skinned_norm[i] = normalize(transform_direction(xform, norm[i]));
    }
}

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles) {
    auto edges = std::vector<vec2i>();
    auto eset = std::unordered_set<vec2i>();
    for (auto t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges.push_back(e);
        }
    }
    return edges;
}

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec4i>& quads) {
    auto edges = std::vector<vec2i>();
    auto eset = std::unordered_set<vec2i>();
    for (auto q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges.push_back(e);
        }
    }
    return edges;
}

// Convert quads to triangles
std::vector<vec3i> convert_quads_to_triangles(const std::vector<vec4i>& quads) {
    auto triangles = std::vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
}

// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
std::vector<vec3i> convert_quads_to_triangles(
    const std::vector<vec4i>& quads, int row_length) {
    auto triangles = std::vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
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

// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec2i> convert_bezier_to_lines(const std::vector<vec4i>& beziers) {
    auto lines = std::vector<vec2i>();
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
void convert_face_varying(std::vector<vec4i>& qquads, std::vector<vec3f>& qpos,
    std::vector<vec3f>& qnorm, std::vector<vec2f>& qtexcoord,
    std::vector<vec4f>& qcolor, const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec4i>& quads_color, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color) {
    // make faces unique
    std::unordered_map<vec4i, int> vert_map;
    qquads = std::vector<vec4i>(quads_pos.size());
    for (auto fid = 0; fid < quads_pos.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec4i{
                (&quads_pos[fid].x)[c],
                (!quads_norm.empty()) ? (&quads_norm[fid].x)[c] : -1,
                (!quads_texcoord.empty()) ? (&quads_texcoord[fid].x)[c] : -1,
                (!quads_color.empty()) ? (&quads_color[fid].x)[c] : -1,
            };
            auto it = vert_map.find(v);
            if (it == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map.insert(it, {v, s});
                (&qquads[fid].x)[c] = s;
            } else {
                (&qquads[fid].x)[c] = it->second;
            }
        }
    }

    // fill vert data
    qpos.clear();
    if (!pos.empty()) {
        qpos.resize(vert_map.size());
        for (auto kv : vert_map) { qpos[kv.second] = pos[kv.first.x]; }
    }
    qnorm.clear();
    if (!norm.empty()) {
        qnorm.resize(vert_map.size());
        for (auto kv : vert_map) { qnorm[kv.second] = norm[kv.first.y]; }
    }
    qtexcoord.clear();
    if (!texcoord.empty()) {
        qtexcoord.resize(vert_map.size());
        for (auto kv : vert_map) {
            qtexcoord[kv.second] = texcoord[kv.first.z];
        }
    }
}

// Subdivide lines.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, int level) {
    if (lines.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto tlines = std::vector<vec2i>();
        for (auto l : lines) {
            tlines.push_back({l.x, (int)tvert.size()});
            tlines.push_back({(int)tvert.size(), l.y});
            tvert.push_back((vert[l.x] + vert[l.y]) / 2);
        }

        std::swap(vert, tvert);
        std::swap(lines, tlines);
    }
}

template void subdivide_lines(std::vector<vec2i>&, std::vector<float>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec2f>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec3f>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec4f>&, int);

// Subdivide triangle.
template <typename T>
void subdivide_triangles(
    std::vector<vec3i>& triangles, std::vector<T>& vert, int level) {
    if (triangles.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto ttriangles = std::vector<vec3i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto t : triangles) {
            for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            ttriangles.push_back(
                {t.x, emap.at({t.x, t.y}), emap.at({t.z, t.x})});
            ttriangles.push_back(
                {t.y, emap.at({t.y, t.z}), emap.at({t.x, t.y})});
            ttriangles.push_back(
                {t.z, emap.at({t.z, t.x}), emap.at({t.y, t.z})});
            ttriangles.push_back({emap.at({t.x, t.y}), emap.at({t.y, t.z}),
                emap.at({t.z, t.x})});
        }

        std::swap(vert, tvert);
        std::swap(triangles, ttriangles);
    }
}

template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<float>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec2f>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec3f>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec4f>&, int);

// Subdivide quads.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, int level) {
    if (quads.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto tquads = std::vector<vec4i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            if (q.z != q.w) {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.w, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.w}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tquads.push_back({q.w, emap.at({q.w, q.x}), (int)tvert.size(),
                    emap.at({q.z, q.w})});
                tvert.push_back(
                    (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4);
            } else {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.z, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.x}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tvert.push_back((vert[q.x] + vert[q.y] + vert[q.y]) / 3);
            }
        }

        std::swap(vert, tvert);
        std::swap(quads, tquads);
    }
}

template void subdivide_quads(std::vector<vec4i>&, std::vector<float>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec2f>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec3f>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec4f>&, int);

// Subdivide beziers.
template <typename T>
void subdivide_beziers(
    std::vector<vec4i>& beziers, std::vector<T>& vert, int level) {
    if (beziers.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto vmap = std::unordered_map<int, int>();
        auto tvert = std::vector<T>();
        auto tbeziers = std::vector<vec4i>();
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
            tvert.push_back(vert[b.x] / 8 + 3 * vert[b.y] / 8 +
                            3 * vert[b.z] / 8 + vert[b.w] / 8);
            tvert.push_back(vert[b.y] / 4 + vert[b.z] / 2 + vert[b.w] / 4);
            tvert.push_back(vert[b.z] / 2 + vert[b.w] / 2);
        }

        std::swap(vert, tvert);
        std::swap(beziers, tbeziers);
    }
}

template void subdivide_beziers(std::vector<vec4i>&, std::vector<float>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec2f>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec3f>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec4f>&, int);

// Subdivide catmullclark.
template <typename T>
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<T>& vert,
    int level, bool lock_boundary) {
    if (quads.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        // split elements ------------------------------------
        auto tvert = vert;
        auto tquads = std::vector<vec4i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            if (q.z != q.w) {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.w, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.w}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tquads.push_back({q.w, emap.at({q.w, q.x}), (int)tvert.size(),
                    emap.at({q.z, q.w})});
                tvert.push_back(
                    (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4);
            } else {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.z, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.x}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tvert.push_back((vert[q.x] + vert[q.y] + vert[q.y]) / 3);
            }
        }

        // find boundary -------------------------------------
        auto eset = std::unordered_set<vec2i>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                eset.insert({e.x, e.y});
            }
        }
        auto boundary = std::vector<vec2i>();
        for (auto e : eset) {
            if (eset.find({e.y, e.x}) != eset.end()) continue;
            boundary.push_back(e);
        }

        // split boundary
        auto tboundary = std::vector<vec2i>();
        for (auto e : boundary) {
            tboundary.push_back({e.x, emap.at(e)});
            tboundary.push_back({emap.at(e), e.y});
        }

        // setup creases -----------------------------------
        auto tcrease_edges = std::vector<vec2i>();
        auto tcrease_verts = std::vector<int>();
        if (lock_boundary) {
            for (auto b : tboundary) {
                tcrease_verts.push_back(b.x);
                tcrease_verts.push_back(b.y);
            }
        } else {
            for (auto b : tboundary) tcrease_edges.push_back(b);
        }

        // define vertex valence ---------------------------
        auto tvert_val = std::vector<int>(tvert.size(), 2);
        for (auto e : tboundary) {
            tvert_val[e.x] = (lock_boundary) ? 0 : 1;
            tvert_val[e.y] = (lock_boundary) ? 0 : 1;
        }

        // averaging pass ----------------------------------
        auto avert = std::vector<T>(tvert.size(), T());
        auto acount = std::vector<int>(tvert.size(), 0);
        for (auto p : tcrease_verts) {
            if (tvert_val[p] != 0) continue;
            avert[p] += tvert[p];
            acount[p] += 1;
        }
        for (auto e : tcrease_edges) {
            auto c = (tvert[e.x] + tvert[e.y]) / 2.0f;
            for (auto vid : {e.x, e.y}) {
                if (tvert_val[vid] != 1) continue;
                avert[vid] += c;
                acount[vid] += 1;
            }
        }
        for (auto q : tquads) {
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

        std::swap(vert, tvert);
        std::swap(quads, tquads);
    }
}

template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<float>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec2f>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec3f>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec4f>&, int, bool);

// Merge lines between shapes.
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    lines.reserve(lines.size() + lines1.size());
    auto nverts = (int)pos.size();
    for (auto l : lines1) lines.push_back({l.x + nverts, l.y + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Merge triangles between shapes.
void merge_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec3i>& triangles1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    triangles.reserve(triangles.size() + triangles1.size());
    auto nverts = (int)pos.size();
    for (auto t : triangles1)
        triangles.push_back({t.x + nverts, t.y + nverts, t.z + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Merge quads between shapes.
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec4i>& quads1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    quads.reserve(quads.size() + quads1.size());
    auto nverts = (int)pos.size();
    for (auto q : quads1)
        quads.push_back(
            {q.x + nverts, q.y + nverts, q.z + nverts, q.w + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Weld vertices within a threshold. For noe the implementation is O(n^2).
std::vector<int> weld_vertices(std::vector<vec3f>& pos, float threshold) {
    auto vid = std::vector<int>(pos.size());
    auto wpos = std::vector<vec3f>();
    for (auto i = 0; i < pos.size(); i++) {
        vid[i] = (int)wpos.size();
        for (auto j = 0; j < wpos.size(); j++) {
            if (length(pos[i] - wpos[j]) < threshold) {
                vid[i] = j;
                break;
            }
        }
        if (vid[i] == (int)wpos.size()) wpos.push_back(pos[i]);
    }
    std::swap(pos, wpos);
    return vid;
}
void weld_triangles(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, float threshold) {
    auto vid = weld_vertices(pos, threshold);
    for (auto t : triangles) {
        t.x = vid[t.x];
        t.y = vid[t.y];
        t.z = vid[t.z];
    }
}
void weld_quads(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, float threshold) {
    auto vid = weld_vertices(pos, threshold);
    for (auto q : quads) {
        q.x = vid[q.x];
        q.y = vid[q.y];
        q.z = vid[q.z];
        q.w = vid[q.w];
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, int seed) {
    auto sampled_pos = std::vector<vec3f>(npoints);
    auto sampled_norm = std::vector<vec3f>(npoints);
    auto sampled_texcoord = std::vector<vec2f>(npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto ei = 0;
        auto uv = zero2f;
        std::tie(ei, uv) =
            sample_triangles(cdf, rand1f(rng), {rand1f(rng), rand1f(rng)});
        auto t = triangles[ei];
        sampled_pos[i] = interpolate_triangle(pos[t.x], pos[t.y], pos[t.z], uv);
        if (!sampled_norm.empty()) {
            sampled_norm[i] = normalize(
                interpolate_triangle(norm[t.x], norm[t.y], norm[t.z], uv));
        } else {
            sampled_norm[i] = triangle_normal(pos[t.x], pos[t.y], pos[t.z]);
        }
        if (!sampled_texcoord.empty()) {
            sampled_texcoord[i] = interpolate_triangle(
                texcoord[t.x], texcoord[t.y], texcoord[t.z], uv);
        } else {
            sampled_texcoord[i] = zero2f;
        }
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

// Make a quad.
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize) {
    auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };

    pos.resize((steps.x + 1) * (steps.y + 1));
    norm.resize((steps.x + 1) * (steps.y + 1));
    texcoord.resize((steps.x + 1) * (steps.y + 1));
    for (auto j = 0; j <= steps.y; j++) {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
            pos[vid(i, j)] = {
                (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
            norm[vid(i, j)] = {0, 0, 1};
            texcoord[vid(i, j)] = uv * uvsize;
        }
    }

    quads.resize(steps.x * steps.y);
    for (auto j = 0; j < steps.y; j++) {
        for (auto i = 0; i < steps.x; i++) {
            quads[j * steps.x + i] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
    }
}

// Make a stack of quads
void make_quad_stack(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec2f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    for (auto i = 0; i <= steps.z; i++) {
        make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
            {size.x, size.y}, uvsize);
        for (auto& p : qpos) p.z = (-0.5f + (float)i / steps.z) * size.z;
        merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    }
}

// Make a cube.
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    // +z
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, qpos[i].y, size.z / 2};
        qnorm[i] = {0, 0, 1};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // -z
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {-qpos[i].x, qpos[i].y, -size.z / 2};
        qnorm[i] = {0, 0, -1};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +x
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.z, steps.y},
        {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {size.x / 2, qpos[i].y, -qpos[i].x};
        qnorm[i] = {1, 0, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // -x
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.z, steps.y},
        {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {-size.x / 2, qpos[i].y, qpos[i].x};
        qnorm[i] = {-1, 0, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +y
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, size.y / 2, -qpos[i].y};
        qnorm[i] = {0, 1, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +y
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, -size.y / 2, qpos[i].y};
        qnorm[i] = {0, -1, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
}

// Make a rounded cube.
void make_cube_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize, float radius) {
    make_cube(quads, pos, norm, texcoord, steps, size, uvsize);
    auto c = size / 2 - vec3f{radius, radius, radius};
    for (auto i = 0; i < pos.size(); i++) {
        auto pc = vec3f{fabs(pos[i].x), fabs(pos[i].y), fabs(pos[i].z)};
        auto ps = vec3f{pos[i].x < 0 ? -1.0f : 1.0f,
            pos[i].y < 0 ? -1.0f : 1.0f, pos[i].z < 0 ? -1.0f : 1.0f};
        if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
            auto pn = normalize(pc - c);
            pos[i] = c + radius * pn;
            norm[i] = pn;
        } else if (pc.x >= c.x && pc.y >= c.y) {
            auto pn = normalize((pc - c) * vec3f{1, 1, 0});
            pos[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
            norm[i] = pn;
        } else if (pc.x >= c.x && pc.z >= c.z) {
            auto pn = normalize((pc - c) * vec3f{1, 0, 1});
            pos[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
            norm[i] = pn;
        } else if (pc.y >= c.y && pc.z >= c.z) {
            auto pn = normalize((pc - c) * vec3f{0, 1, 1});
            pos[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
            norm[i] = pn;
        } else {
            continue;
        }
        pos[i] *= ps;
        norm[i] *= ps;
    }
}

// Make a sphere.
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pi * uv.x, pi * (1 - uv.y)};
        auto p = vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        pos[i] = p * (size / 2);
        norm[i] = normalize(p);
        texcoord[i] = uv * uvsize;
    }
}

// Make a spherecube.
void make_sphere_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize) {
    make_cube(quads, pos, norm, texcoord, {steps, steps, steps}, {1, 1, 1},
        {uvsize, uvsize, uvsize});
    for (auto i = 0; i < pos.size(); i++) {
        auto p = pos[i];
        pos[i] = normalize(p) * (size / 2);
        norm[i] = normalize(p);
    }
}

// Make a flipped sphere. This is not watertight.
void make_sphere_flipcap(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize, const vec2f& zflip) {
    make_sphere(quads, pos, norm, texcoord, steps, size, uvsize);
    for (auto i = 0; i < pos.size(); i++) {
        if (pos[i].z > zflip.y) {
            pos[i].z = 2 * zflip.y - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        } else if (pos[i].z < zflip.x) {
            pos[i].z = 2 * zflip.x - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        }
    }
}

// Make a disk.
void make_disk(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto phi = 2 * pi * uv.x;
        pos[i] = {cos(phi) * uv.y * size / 2, sin(phi) * uv.y * size / 2, 0};
        norm[i] = {0, 0, 1};
        texcoord[i] = uv * uvsize;
    }
}

// Make a disk from a quad.
void make_disk_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize) {
    make_quad(
        quads, pos, norm, texcoord, {steps, steps}, {2, 2}, {uvsize, uvsize});
    for (auto i = 0; i < pos.size(); i++) {
        // Analytical Methods for Squaring the Disc, by C. Fong
        // https://arxiv.org/abs/1509.06344
        auto xy = vec2f{pos[i].x, pos[i].y};
        auto uv = vec2f{
            xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
        pos[i] = {uv.x * size / 2, uv.y * size / 2, 0};
    }
}

// Make a bulged disk from a quad.
void make_disk_bulged(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize, float height) {
    make_disk_quad(quads, pos, norm, texcoord, steps, size, uvsize);
    if (height == 0) return;
    auto radius = (size * size / 4 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i = 0; i < pos.size(); i++) {
        auto pn = normalize(pos[i] - center);
        pos[i] = center + pn * radius;
        norm[i] = pn;
    }
}

// Make a cylinder (side-only).
void make_cylinder_side(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto phi = 2 * pi * uv.x;
        pos[i] = {cos(phi) * size.x / 2, sin(phi) * size.x / 2,
            (uv.y - 0.5f) * size.y};
        norm[i] = {cos(phi), sin(phi), 0};
        texcoord[i] = uv * uvsize;
    }
}

// Make a cylinder.
void make_cylinder(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    // side
    make_cylinder_side(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // top
    make_disk(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.z}, size.x,
        {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpos.size(); i++) { qpos[i].z = size.y / 2; }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // bottom
    make_disk(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.z}, size.x,
        {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i].z = -size.y / 2;
        qnorm[i] = -qnorm[i];
    }
    for (auto i = 0; i < qquads.size(); i++)
        std::swap(qquads[i].x, qquads[i].z);
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
}

// Make a rounded cylinder.
void make_cylinder_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius) {
    make_cylinder(quads, pos, norm, texcoord, steps, size, uvsize);
    auto c = size / 2 - vec2f{radius, radius};
    for (auto i = 0; i < pos.size(); i++) {
        auto phi = atan2(pos[i].y, pos[i].x);
        auto r = length(vec2f{pos[i].x, pos[i].y});
        auto z = pos[i].z;
        auto pc = vec2f{r, fabs(z)};
        auto ps = (z < 0) ? -1.0f : 1.0f;
        if (pc.x >= c.x && pc.y >= c.y) {
            auto pn = normalize(pc - c);
            pos[i] = {cos(phi) * c.x + radius * pn.x,
                sin(phi) * c.x + radius * pn.x, ps * (c.y + radius * pn.y)};
            norm[i] = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
        } else {
            continue;
        }
    }
}

// Make a geodesic sphere.
void make_geodesic_sphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, int tesselation, float size) {
    // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
    const float X = 0.525731112119133606f;
    const float Z = 0.850650808352039932f;
    pos = std::vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z},
        {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
        {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
    triangles = std::vector<vec3i>{{0, 1, 4}, {0, 4, 9}, {9, 4, 5}, {4, 8, 5},
        {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3}, {5, 3, 2}, {2, 3, 7},
        {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0}, {0, 6, 1}, {6, 10, 1},
        {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
    subdivide_triangles(triangles, pos, max(0, tesselation - 2));
    for (auto& p : pos) p = normalize(p) * size / 2;
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    const vec3i& steps, const vec3f& size, const vec3f& uvsize) {
    make_cube(quads_pos, pos, norm, texcoord, steps, size, uvsize);
    quads_norm = quads_pos;
    quads_texcoord = quads_pos;
    weld_quads(quads_pos, pos,
        min(0.1f * size /
            vec3f{(float)steps.x, (float)steps.y, (float)steps.z}));
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
void make_suzanne(std::vector<vec4i>& quads, std::vector<vec3f>& pos) {
    static auto suzanne_pos = std::vector<vec3f>{{0.4375, 0.1640625, 0.765625},
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
    static auto suzanne_triangles = std::vector<vec3i>{{60, 64, 48},
        {49, 65, 61}, {62, 64, 60}, {61, 65, 63}, {60, 58, 62}, {63, 59, 61},
        {60, 56, 58}, {59, 57, 61}, {60, 54, 56}, {57, 55, 61}, {60, 52, 54},
        {55, 53, 61}, {60, 50, 52}, {53, 51, 61}, {60, 48, 50}, {51, 49, 61},
        {224, 228, 226}, {227, 229, 225}, {72, 283, 73}, {73, 284, 72},
        {341, 347, 383}, {384, 348, 342}, {299, 345, 343}, {344, 346, 300},
        {323, 379, 351}, {352, 380, 324}, {441, 443, 445}, {446, 444, 442},
        {463, 491, 465}, {466, 492, 464}, {495, 497, 499}, {500, 498, 496}};
    static auto suzanne_quads = std::vector<vec4i>{{46, 0, 2, 44},
        {3, 1, 47, 45}, {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4},
        {7, 9, 3, 5}, {0, 10, 8, 2}, {9, 11, 1, 3}, {10, 12, 14, 8},
        {15, 13, 11, 9}, {8, 14, 16, 6}, {17, 15, 9, 7}, {14, 20, 18, 16},
        {19, 21, 15, 17}, {12, 22, 20, 14}, {21, 23, 13, 15}, {22, 24, 26, 20},
        {27, 25, 23, 21}, {20, 26, 28, 18}, {29, 27, 21, 19}, {26, 32, 30, 28},
        {31, 33, 27, 29}, {24, 34, 32, 26}, {33, 35, 25, 27}, {34, 36, 38, 32},
        {39, 37, 35, 33}, {32, 38, 40, 30}, {41, 39, 33, 31}, {38, 44, 42, 40},
        {43, 45, 39, 41}, {36, 46, 44, 38}, {45, 47, 37, 39}, {46, 36, 50, 48},
        {51, 37, 47, 49}, {36, 34, 52, 50}, {53, 35, 37, 51}, {34, 24, 54, 52},
        {55, 25, 35, 53}, {24, 22, 56, 54}, {57, 23, 25, 55}, {22, 12, 58, 56},
        {59, 13, 23, 57}, {12, 10, 62, 58}, {63, 11, 13, 59}, {10, 0, 64, 62},
        {65, 1, 11, 63}, {0, 46, 48, 64}, {49, 47, 1, 65}, {88, 173, 175, 90},
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

    pos = suzanne_pos;
    quads = suzanne_quads;
    for (auto t : suzanne_triangles) { quads.push_back({t.x, t.y, t.z, t.z}); }
}

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius) {
    auto nverts = (steps.x + 1) * steps.y;
    auto nlines = steps.x * steps.y;
    auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };
    auto fid = [steps](int i, int j) { return j * steps.x + i; };

    pos.resize(nverts);
    norm.resize(nverts);
    texcoord.resize(nverts);
    radius.resize(nverts);
    if (steps.y > 1) {
        for (auto j = 0; j < steps.y; j++) {
            for (auto i = 0; i <= steps.x; i++) {
                auto uv = vec2f{i / (float)steps.x,
                    j / (float)(steps.y > 1 ? steps.y - 1 : 1)};
                pos[vid(i, j)] = {
                    (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
                norm[vid(i, j)] = {1, 0, 0};
                texcoord[vid(i, j)] = uv * uvsize;
            }
        }
    } else {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, 0};
            pos[vid(i, 0)] = {(uv.x - 0.5f) * size.x, 0, 0};
            norm[vid(i, 0)] = {1, 0, 0};
            texcoord[vid(i, 0)] = uv * uvsize;
        }
    }

    lines.resize(nlines);
    for (int j = 0; j < steps.y; j++) {
        for (int i = 0; i < steps.x; i++) {
            lines[fid(i, j)] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

// Generate a point set with points placed at the origin with texcoords varying
// along u.
void make_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, float uvsize, float point_radius) {
    points.resize(num);
    for (auto i = 0; i < num; i++) points[i] = i;
    pos.assign(num, {0, 0, 0});
    norm.assign(num, {0, 0, 1});
    texcoord.assign(num, {0, 0});
    radius.assign(num, point_radius);
    for (auto i = 0; i < texcoord.size(); i++)
        texcoord[i] = {(float)i / (float)num, 0};
}

// Generate a point set.
void make_random_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, const vec3f& size, float uvsize,
    float point_radius, uint64_t seed) {
    make_points(points, pos, norm, texcoord, radius, num, uvsize, point_radius);
    auto rng = make_rng(seed);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = (rand3f(rng) - vec3f{0.5f, 0.5f, 0.5f}) * size;
    }
}

// Make a point.
void make_point(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, float point_radius) {
    points = {0};
    pos = {{0, 0, 0}};
    norm = {{0, 0, 1}};
    texcoord = {{0, 0}};
    radius = {point_radius};
}

// Make a bezier circle. Returns bezier, pos.
void make_bezier_circle(std::vector<vec4i>& beziers, std::vector<vec3f>& pos) {
    // constant from http://spencermortensen.com/articles/bezier-circle/
    auto c = 0.551915024494f;
    pos = std::vector<vec3f>{{1, 0, 0}, {1, c, 0}, {c, 1, 0}, {0, 1, 0},
        {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
        {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
    beziers = std::vector<vec4i>{
        {0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
}

// Make a hair ball around a shape
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
    std::vector<vec3f> bpos;
    std::vector<vec3f> bnorm;
    std::vector<vec2f> btexcoord;
    std::tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        striangles, spos, snorm, stexcoord, steps.y, seed);

    auto rng = make_rng(seed, 3);
    auto blen = std::vector<float>(bpos.size());
    for (auto& l : blen) l = lerp(len.x, len.y, rand1f(rng));

    auto cidx = std::vector<int>();
    if (clump.x > 0) {
        for (auto bidx = 0; bidx < bpos.size(); bidx++) {
            cidx.push_back(0);
            auto cdist = flt_max;
            for (auto c = 0; c < clump.y; c++) {
                auto d = length(bpos[bidx] - bpos[c]);
                if (d < cdist) {
                    cdist = d;
                    cidx.back() = c;
                }
            }
        }
    }

    make_lines(lines, pos, norm, texcoord, radius, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto u = texcoord[i].x;
        auto bidx = i / (steps.x + 1);
        pos[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
        norm[i] = bnorm[bidx];
        radius[i] = lerp(rad.x, rad.y, u);
        if (clump.x > 0) {
            pos[i] = pos[i] +
                     (pos[i + (cidx[bidx] - bidx) * (steps.x + 1)] - pos[i]) *
                         u * clump.x;
        }
        if (noise.x > 0) {
            auto nx = perlin_noise(pos[i] * noise.y + vec3f{0, 0, 0}) * noise.x;
            auto ny =
                perlin_noise(pos[i] * noise.y + vec3f{3, 7, 11}) * noise.x;
            auto nz =
                perlin_noise(pos[i] * noise.y + vec3f{13, 17, 19}) * noise.x;
            pos[i] += {nx, ny, nz};
        }
    }

    if (clump.x > 0 || noise.x > 0 || rotation.x > 0)
        compute_tangents(lines, pos, norm);
}

}  // namespace ygl
