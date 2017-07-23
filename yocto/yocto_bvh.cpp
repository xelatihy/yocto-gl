//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF YOCTO_BVH
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"

#include <algorithm>
#include <cstdio>
#include <unordered_map>

namespace ybvh {

// -----------------------------------------------------------------------------
// BVH DATA STRUCTURE
// -----------------------------------------------------------------------------

// number of primitives to avoid splitting on
#define YBVH__MINPRIMS 4

//
// BVH tree node containing its bounds, indices to the BVH arrays of either
// sorted primitives or internal nodes, whether its a leaf or an internal node,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh for more details.
//
// This is not part of the public interface.
//
// Implemenetation Notes:
// - Padded to 32 bytes for cache fiendly access
//
struct bvh_node {
    ym::bbox3f bbox;  // bounding box
    uint32_t start;   // index to the first sorted primitive/node
    uint16_t count;   // number of primitives/nodes
    uint8_t isleaf;   // whether it is a leaf
    uint8_t axis;     // slit axis
};

//
// BVH tree, stored as a node array. The tree structure is encoded using array
// indices instead of pointers, both for speed but also to simplify code.
// BVH nodes indices refer to either the node array, for internal nodes,
// or a primitive array, for leaf nodes. BVH trees may contain only one type
// of geometric primitive, like points, lines, triangle or shape other BVHs.
// We handle multiple primitive types and transformed primitices by building
// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
// BVHs, shape BVHs, each of which of a uniform primitive type.
//
// This is not part of the public interface.
//
struct bvh_tree {
    // bvh data
    std::vector<bvh_node> nodes;   // sorted array of internal nodes
    std::vector<int> sorted_prim;  // sorted elements
};

//
// Shape Data as index mesh. Only one element type is supported at one time.
//
struct shape {
    // shape id ---------------------------
    int sid = -1;  // shape id

    // elements data ----------------------
    int nelems = 0;                       // number of elements
    const int* point = nullptr;           // point indices
    const ym::vec2i* line = nullptr;      // line indices
    const ym::vec3i* triangle = nullptr;  // triangle indices
    const ym::vec4i* tetra = nullptr;     // tetrahedra indices

    // vertex data ------------------------
    int nverts = 0;                  // number of vertices
    const ym::vec3f* pos = nullptr;  // vertex pos
    const float* radius = nullptr;   // vertex radius

    // [private] bvh data -----------------
    bvh_tree* bvh = nullptr;               // bvh [private]
    ym::bbox3f bbox = ym::invalid_bbox3f;  // bbox [private]

    // [private] methods ------------------
    float rad(int i) const { return (radius) ? radius[i] : 0; }

    // destructor
    ~shape() {
        if (bvh) delete bvh;
    }
};

//
// Instance
//
struct instance {
    // instance id -----------------------
    int iid = -1;  // instance id

    // instance transform ----------------
    ym::mat4f xform;      // instance transform
    ym::mat4f xform_inv;  // instance transform inverse

    // instance data ---------------------
    shape* shp = nullptr;  // shape

    // [private] bvh data -----------------
    ym::bbox3f bbox = ym::invalid_bbox3f;  // bbox [private]
};

//
// Scene Data
//
struct scene {
    // scene data -------------------------
    std::vector<instance*> instances;  // instances
    std::vector<shape*> shapes;        // shapes

    // bvh private data -------------------
    bvh_tree* bvh = nullptr;  // bvh [private]

    // [private] methods -----------------
    ym::bbox3f bbox() const { return bvh->nodes[0].bbox; }

    // destructor
    ~scene() {
        for (auto& shp : shapes) delete shp;
        for (auto& ist : instances) delete ist;
        if (bvh) delete bvh;
    }
};

// -----------------------------------------------------------------------------
// SCENE SPECIFICION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Init scene.
//
scene* make_scene() { return new scene(); }

///
/// Free scene.
///
void free_scene(scene*& scn) {
    if (scn) delete scn;
    scn = nullptr;
}

//
// Set shape. Public API.
//
int add_point_shape(scene* scn, int npoints, const int* point, int nverts,
    const ym::vec3f* pos, const float* radius) {
    auto shp = new shape();
    shp->sid = (int)scn->shapes.size();
    shp->nelems = npoints;
    shp->point = (const int*)point;
    shp->line = nullptr;
    shp->triangle = nullptr;
    shp->tetra = nullptr;
    shp->nverts = nverts;
    shp->pos = (const ym::vec3f*)pos;
    shp->radius = radius;
    scn->shapes.push_back(shp);
    return shp->sid;
}

//
// Set shape. Public API.
//
int add_line_shape(scene* scn, int nlines, const ym::vec2i* lines, int nverts,
    const ym::vec3f* pos, const float* radius) {
    auto shp = new shape();
    shp->sid = (int)scn->shapes.size();
    shp->nelems = nlines;
    shp->point = nullptr;
    shp->line = (const ym::vec2i*)lines;
    shp->triangle = nullptr;
    shp->tetra = nullptr;
    shp->nverts = nverts;
    shp->pos = (const ym::vec3f*)pos;
    shp->radius = radius;
    scn->shapes.push_back(shp);
    return shp->sid;
}

//
// Set shape. Public API.
//
int add_triangle_shape(scene* scn, int ntriangles, const ym::vec3i* triangles,
    int nverts, const ym::vec3f* pos, const float* radius) {
    auto shp = new shape();
    shp->sid = (int)scn->shapes.size();
    shp->nelems = ntriangles;
    shp->point = nullptr;
    shp->line = nullptr;
    shp->triangle = (const ym::vec3i*)triangles;
    shp->tetra = nullptr;
    shp->nverts = nverts;
    shp->pos = (const ym::vec3f*)pos;
    shp->radius = radius;
    scn->shapes.push_back(shp);
    return shp->sid;
}

//
// Set shape. Public API.
//
int add_tetra_shape(scene* scn, int ntetra, const ym::vec4i* tetra, int nverts,
    const ym::vec3f* pos, const float* radius) {
    auto shp = new shape();
    shp->sid = (int)scn->shapes.size();
    shp->nelems = ntetra;
    shp->point = nullptr;
    shp->line = nullptr;
    shp->triangle = nullptr;
    shp->tetra = (const ym::vec4i*)tetra;
    shp->nverts = nverts;
    shp->pos = (const ym::vec3f*)pos;
    shp->radius = radius;
    scn->shapes.push_back(shp);
    return shp->sid;
}

//
// Set shape. Public API.
//
int add_point_shape(
    scene* scn, int nverts, const ym::vec3f* pos, const float* radius) {
    auto shp = new shape();
    shp->sid = (int)scn->shapes.size();
    shp->nelems = nverts;
    shp->point = nullptr;
    shp->line = nullptr;
    shp->triangle = nullptr;
    shp->tetra = nullptr;
    shp->nverts = nverts;
    shp->pos = (const ym::vec3f*)pos;
    shp->radius = radius;
    scn->shapes.push_back(shp);
    return shp->sid;
}

//
// Set instance frame. Public API.
//
int add_instance(
    scene* scn, const ym::mat4f& xf, const ym::mat4f& xf_inv, int sid) {
    auto ist = new instance();
    ist->iid = (int)scn->instances.size();
    ist->xform = xf;
    ist->xform_inv = xf_inv;
    ist->shp = scn->shapes[sid];
    scn->instances.push_back(ist);
    return ist->iid;
}

//
// Set instance frame. Public API.
//
void set_instance_transform(
    scene* scn, int iid, const ym::mat4f& xform, const ym::mat4f& xform_inv) {
    scn->instances[iid]->xform = xform;
    scn->instances[iid]->xform_inv = xform_inv;
}

// -----------------------------------------------------------------------------
// BVH BUILD FUNCTIONS
// -----------------------------------------------------------------------------

//
// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
//
struct bound_prim {
    ym::bbox3f bbox;       // bounding box
    ym::vec3f center;      // bounding box center (for faster sort)
    int pid;               // primitive id
    float sah_cost_left;   // buffer for sah heuristic costs
    float sah_cost_right;  // buffer for sah heuristic costs
};

//
// Comparison function for each axis
//
struct bound_prim_comp {
    int axis;
    float middle;

    bound_prim_comp(int a, float m = 0) : axis(a), middle(m) {}

    bool operator()(const bound_prim& a, const bound_prim& b) const {
        return a.center[axis] < b.center[axis];
    }

    bool operator()(const bound_prim& a) const {
        return a.center[axis] < middle;
    }
};

//
// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
//
void make_node(bvh_node* node, std::vector<bvh_node>& nodes,
    bound_prim* sorted_prims, int start, int end, bool equalsize) {
    // compute node bounds
    node->bbox = ym::invalid_bbox3f;
    for (auto i = start; i < end; i++) node->bbox += sorted_prims[i].bbox;

    // decide whether to create a leaf
    if (end - start <= YBVH__MINPRIMS) {
        // makes a leaf node
        node->isleaf = true;
        node->start = start;
        node->count = end - start;
    } else {
        // choose the split axis and position
        // init to default values
        auto axis = 0;
        auto mid = (start + end) / 2;

        // compute primintive bounds and size
        auto centroid_bbox = ym::invalid_bbox3f;
        for (auto i = start; i < end; i++)
            centroid_bbox += sorted_prims[i].center;
        auto centroid_size = ym::diagonal(centroid_bbox);

        // check if it is not possible to split
        if (centroid_size == ym::zero3f) {
            // we failed to split for some reasons
            node->isleaf = true;
            node->start = start;
            node->count = end - start;
        } else {
            // split along largest
            auto largest_axis = ym::max_element_idx(centroid_size);

            // check heuristic
            if (equalsize) {
                // split the space in the middle along the largest axis
                axis = largest_axis;
                mid = (int)(std::partition(sorted_prims + start,
                                sorted_prims + end,
                                bound_prim_comp(largest_axis,
                                    ym::center(centroid_bbox)[largest_axis])) -
                            sorted_prims);
            } else {
                // balanced tree split: find the largest axis of the bounding
                // box and split along this one right in the middle
                axis = largest_axis;
                mid = (start + end) / 2;
                std::nth_element(sorted_prims + start, sorted_prims + mid,
                    sorted_prims + end, bound_prim_comp(largest_axis));
            }

            // check correctness
            assert(axis >= 0 && mid > 0);
            assert(mid > start && mid < end);

            // makes an internal node
            node->isleaf = false;
            // perform the splits by preallocating the child nodes and recurring
            node->axis = axis;
            node->start = (int)nodes.size();
            node->count = 2;
            nodes.emplace_back();
            nodes.emplace_back();
            // build child nodes
            make_node(&nodes[node->start], nodes, sorted_prims, start, mid,
                equalsize);
            make_node(&nodes[node->start + 1], nodes, sorted_prims, mid, end,
                equalsize);
        }
    }
}

//
// Build a BVH from a set of primitives.
//
template <typename ElemBbox>
void build_bvh(
    bvh_tree*& bvh, int nprims, bool equalsize, const ElemBbox& elem_bbox) {
    // allocate if needed
    if (bvh) delete bvh;
    bvh = new bvh_tree();

    // prepare prims
    auto bound_prims = std::vector<bound_prim>(nprims);
    for (auto i = 0; i < nprims; i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = elem_bbox(i);
        bound_prims[i].center = ym::center(bound_prims[i].bbox);
    }

    // clear bvh
    bvh->nodes.clear();
    bvh->sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh->nodes.reserve(nprims * 2);

    // start recursive splitting
    bvh->nodes.emplace_back();
    make_node(
        &bvh->nodes[0], bvh->nodes, bound_prims.data(), 0, nprims, equalsize);

    // shrink back
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(nprims);
    for (int i = 0; i < nprims; i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }
}

//
// Build a shape BVH. Public function whose interface is described above.
//
void build_shape_bvh(shape* shp, bool equalsize) {
    if (shp->point) {
        build_bvh(shp->bvh, shp->nelems, equalsize, [shp](int eid) {
            auto f = shp->point[eid];
            return point_bbox(shp->pos[f], shp->rad(f));
        });
    } else if (shp->line) {
        build_bvh(shp->bvh, shp->nelems, equalsize, [shp](int eid) {
            auto f = shp->line[eid];
            return line_bbox(
                shp->pos[f.x], shp->pos[f.y], shp->rad(f.x), shp->rad(f.y));
        });
    } else if (shp->triangle) {
        build_bvh(shp->bvh, shp->nelems, equalsize, [shp](int eid) {
            auto f = shp->triangle[eid];
            return triangle_bbox(shp->pos[f.x], shp->pos[f.y], shp->pos[f.z]);
        });
    } else if (shp->tetra) {
        build_bvh(shp->bvh, shp->nelems, equalsize, [shp](int eid) {
            auto f = shp->tetra[eid];
            return tetrahedron_bbox(
                shp->pos[f.x], shp->pos[f.y], shp->pos[f.z], shp->pos[f.w]);
        });
    } else {
        build_bvh(shp->bvh, shp->nelems, equalsize, [shp](int eid) {
            return point_bbox(shp->pos[eid], shp->rad(eid));
        });
    }
    shp->bbox = shp->bvh->nodes[0].bbox;
}

//
// Build a scene BVH. Public function whose interface is described above.
//
void build_scene_bvh(scene* scn, bool equalsize, bool do_shapes) {
    // do shapes
    if (do_shapes) {
        for (auto shp : scn->shapes) build_shape_bvh(shp, equalsize);
    }

    // update instance bbox
    for (auto ist : scn->instances)
        ist->bbox = ym::transform_bbox(ist->xform, ist->shp->bbox);

    // tree bvh
    build_bvh(scn->bvh, (int)scn->instances.size(), equalsize,
        [scn](int eid) { return scn->instances[eid]->bbox; });
}

//
// Recursively recomputes the node bounds for a shape bvh
//
template <typename ElemBbox>
void refit_bvh(bvh_tree* bvh, int nodeid, const ElemBbox& elem_bbox) {
    // refit
    auto node = &bvh->nodes[nodeid];
    node->bbox = ym::invalid_bbox3f;
    if (node->isleaf) {
        for (auto i = 0; i < node->count; i++) {
            auto idx = bvh->sorted_prim[node->start + i];
            node->bbox += elem_bbox(idx);
        }
    } else {
        for (auto i = 0; i < node->count; i++) {
            auto idx = node->start + i;
            refit_bvh(bvh, idx, elem_bbox);
            node->bbox += bvh->nodes[idx].bbox;
        }
    }
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
void refit_shape_bvh(shape* shp) {
    if (shp->point) {
        refit_bvh(shp->bvh, 0, [shp](int eid) {
            auto f = shp->point[eid];
            return point_bbox(shp->pos[f], shp->rad(f));
        });
    } else if (shp->line) {
        refit_bvh(shp->bvh, 0, [shp](int eid) {
            auto f = shp->line[eid];
            return line_bbox(
                shp->pos[f.x], shp->pos[f.y], shp->rad(f.x), shp->rad(f.y));
        });
    } else if (shp->triangle) {
        refit_bvh(shp->bvh, 0, [shp](int eid) {
            auto f = shp->triangle[eid];
            return triangle_bbox(shp->pos[f.x], shp->pos[f.y], shp->pos[f.z]);
        });
    } else if (shp->tetra) {
        refit_bvh(shp->bvh, 0, [shp](int eid) {
            auto f = shp->tetra[eid];
            return tetrahedron_bbox(
                shp->pos[f.x], shp->pos[f.y], shp->pos[f.z], shp->pos[f.w]);
        });
    } else {
        refit_bvh(shp->bvh, 0, [shp](int eid) {
            return point_bbox(shp->pos[eid], shp->rad(eid));
        });
    }
    shp->bbox = shp->bvh->nodes[0].bbox;
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
void refit_shape_bvh(scene* scn, int sid) {
    return refit_shape_bvh(scn->shapes[sid]);
}

//
// Refits a scene BVH. Public function whose interface is described above.
//
void refit_scene_bvh(scene* scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) refit_shape_bvh(scn, shp->sid);
    }

    // update instance bbox
    for (auto ist : scn->instances)
        ist->bbox = ym::transform_bbox(ist->xform, ist->shp->bbox);

    // recompute bvh bounds
    refit_bvh(
        scn->bvh, 0, [scn](int eid) { return scn->instances[eid]->bbox; });
}

// -----------------------------------------------------------------------------
// BVH INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

//
// Logging global variables
//
#ifdef YGL_BVH_LOG_RAYS
int _log_nrays = 0;
int _log_nbbox_inters = 0;
int _log_npoint_inters = 0;
int _log_nline_inters = 0;
int _log_ntriangle_inters = 0;

//
// Logging information
//
void get_ray_log(int& nrays, int& nbbox_inters, int& npoint_inters,
    int& nline_inters, int& ntriangle_inters) {
    nrays = _log_nrays;
    nbbox_inters = _log_nbbox_inters;
    npoint_inters = _log_npoint_inters;
    nline_inters = _log_nline_inters;
    ntriangle_inters = _log_ntriangle_inters;
}

#endif

//
// Intersect ray with a bvh-> Similar to the generic public function whose
// interface is described above. See intersect_ray for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of recursive
// calls; this follows general conventions and stragely makes the code shorter
// - The walk is simplified for first hit by obeserving that if we update
// the ray_tmax limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequence farthest iterations will be
// rejected in the tmax tests
//
template <typename Isec>
intersection_point intersect_bvh(const bvh_tree* bvh, const ym::ray3f& ray_,
    bool early_exit, const Isec& intersect_elem) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto pt = intersection_point();

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = ym::vec3f{1, 1, 1} / ray.d;
    auto ray_dsign = ym::vec3i{(ray_dinv.x < 0) ? 1 : 0,
        (ray_dinv.y < 0) ? 1 : 0, (ray_dinv.z < 0) ? 1 : 0};
    auto ray_reverse = ym::vec<bool, 4>{
        (bool)ray_dsign.x, (bool)ray_dsign.y, (bool)ray_dsign.z, false};

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!ym::intersect_check_bbox(ray, ray_dinv, ray_dsign, node.bbox))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // for internal nodes, attempts to proceed along the
            // split axis from smallest to largest nodes
            if (ray_reverse[node.axis]) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            } else {
                for (auto i = node.count - 1; i >= 0; i--) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh->sorted_prim[node.start + i];
                auto pp = intersection_point();
                if ((pp = intersect_elem(idx, ray, early_exit))) {
                    if (early_exit) return pp;
                    pt = pp;
                    ray.tmax = pt.dist;
                }
            }
        }
    }

    return pt;
}

//
// Shape intersection
//
intersection_point intersect_shape(
    const shape* shp, const ym::ray3f& ray, bool early_exit) {
    // initialize point
    auto pt = intersection_point();

    // switch over shape type
    if (shp->triangle) {
        pt = intersect_bvh(shp->bvh, ray, early_exit,
            [shp](int eid, const ym::ray3f& ray, bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->triangle[eid];
                if (!ym::intersect_triangle(ray, shp->pos[f.x], shp->pos[f.y],
                        shp->pos[f.z], pt.dist, (ym::vec3f&)pt.euv))
                    return intersection_point{};
                pt.euv = {pt.euv.x, pt.euv.y, pt.euv.z, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->line) {
        assert(shp->radius);
        pt = intersect_bvh(shp->bvh, ray, early_exit,
            [shp](int eid, const ym::ray3f& ray, bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->line[eid];
                if (!ym::intersect_line(ray, shp->pos[f.x], shp->pos[f.y],
                        shp->radius[f.x], shp->radius[f.y], pt.dist,
                        (ym::vec2f&)pt.euv))
                    return intersection_point{};
                pt.euv = {pt.euv.x, pt.euv.y, 0, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->point) {
        assert(shp->radius);
        pt = intersect_bvh(shp->bvh, ray, early_exit,
            [shp](int eid, const ym::ray3f& ray, bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->point[eid];
                if (!ym::intersect_point(
                        ray, shp->pos[f], shp->radius[f], pt.dist))
                    return intersection_point{};
                pt.euv = {1, 0, 0, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->tetra) {
        pt = intersect_bvh(shp->bvh, ray, early_exit,
            [shp](int eid, const ym::ray3f& ray, bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->tetra[eid];
                if (!ym::intersect_tetrahedron(ray, shp->pos[f.x],
                        shp->pos[f.y], shp->pos[f.z], shp->pos[f.w], pt.dist,
                        (ym::vec4f&)pt.euv))
                    return intersection_point{};
                pt.eid = eid;
                return pt;
            });
    } else {
        assert(shp->radius);
        pt = intersect_bvh(shp->bvh, ray, early_exit,
            [shp](int eid, const ym::ray3f& ray, bool early_exit) {
                auto pt = intersection_point();
                if (!ym::intersect_point(
                        ray, shp->pos[eid], shp->radius[eid], pt.dist))
                    return intersection_point{};
                pt.euv = {1, 0, 0, 0};
                pt.eid = eid;
                return pt;
            });
    }
    if (pt.eid >= 0) pt.sid = shp->sid;
    return pt;
}

//
// Instance intersection
//
intersection_point intersect_instance(
    const instance* ist, const ym::ray3f& ray, bool early_exit) {
    auto pt = intersect_shape(
        ist->shp, ym::transform_ray(ist->xform_inv, ray), early_exit);
    if (pt.eid >= 0) pt.iid = ist->iid;
    return pt;
}

//
// Scene intersection
//
intersection_point intersect_scene(
    const scene* scn, const ym::ray3f& ray, bool early_exit) {
    return intersect_bvh(scn->bvh, ray, early_exit,
        [scn](int eid, const ym::ray3f& ray, bool early_exit) {
            return intersect_instance(scn->instances[eid], ray, early_exit);
        });
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP
// -----------------------------------------------------------------------------

//
// Finds the closest element with a bvh.
// Similar to the generic public function whose
// interface is described above. See nightbor_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of recursive
// calls; this follows general conventions and stragely makes the code shorter
// - The walk is simplified for first hit by obeserving that if we update
// the dist_max limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequent farthest iterations will be
// rejected in the tmax tests
//
template <typename OverlapElem>
intersection_point overlap_bvh(const bvh_tree* bvh, const ym::vec3f& pos,
    float max_dist, bool early_exit, const OverlapElem& overlap_elem) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto pt = intersection_point();

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!ym::distance_check_bbox(pos, max_dist, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // internal node
            for (auto idx = node.start; idx < node.start + node.count; idx++) {
                node_stack[node_cur++] = idx;
                assert(node_cur < 64);
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh->sorted_prim[node.start + i];
                auto pp = intersection_point();
                if ((pp = overlap_elem(idx, pos, max_dist, early_exit))) {
                    if (early_exit) return pp;
                    pt = pp;
                    max_dist = pt.dist;
                }
            }
        }
    }

    return pt;
}

//
// Shape overlap
//
intersection_point overlap_shape(
    const shape* shp, const ym::vec3f& pos, float max_dist, bool early_exit) {
    // initialize point
    auto pt = intersection_point();

    // switch over elemenet type
    if (shp->triangle) {
        pt = overlap_bvh(shp->bvh, pos, max_dist, early_exit,
            [shp](int eid, const ym::vec3f& pos, float max_dist,
                bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->triangle[eid];
                if (!ym::overlap_triangle(pos, max_dist, shp->pos[f.x],
                        shp->pos[f.y], shp->pos[f.z], shp->rad(f.x),
                        shp->rad(f.y), shp->rad(f.z), pt.dist,
                        (ym::vec3f&)pt.euv))
                    return intersection_point{};
                pt.euv = {pt.euv.x, pt.euv.y, pt.euv.z, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->line) {
        pt = overlap_bvh(shp->bvh, pos, max_dist, early_exit,
            [shp](int eid, const ym::vec3f& pos, float max_dist,
                bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->line[eid];
                if (!ym::overlap_line(pos, max_dist, shp->pos[f.x],
                        shp->pos[f.y], shp->rad(f.x), shp->rad(f.y), pt.dist,
                        (ym::vec2f&)pt.euv))
                    return intersection_point{};
                pt.euv = {pt.euv.x, pt.euv.y, 0, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->point) {
        pt = overlap_bvh(shp->bvh, pos, max_dist, early_exit,
            [shp](int eid, const ym::vec3f& pos, float max_dist,
                bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->point[eid];
                if (!ym::overlap_point(
                        pos, max_dist, shp->pos[f], shp->rad(f), pt.dist))
                    return intersection_point{};
                pt.euv = {1, 0, 0, 0};
                pt.eid = eid;
                return pt;
            });
    } else if (shp->tetra) {
        pt = overlap_bvh(shp->bvh, pos, max_dist, early_exit,
            [shp](int eid, const ym::vec3f& pos, float max_dist,
                bool early_exit) {
                auto pt = intersection_point();
                auto f = shp->tetra[eid];
                if (!ym::overlap_tetrahedron(pos, max_dist, shp->pos[f.x],
                        shp->pos[f.y], shp->pos[f.z], shp->pos[f.w],
                        shp->rad(f.x), shp->rad(f.y), shp->rad(f.z),
                        shp->rad(f.w), pt.dist, (ym::vec4f&)pt.euv))
                    return intersection_point{};
                pt.eid = eid;
                return pt;
            });
    } else {
        pt = overlap_bvh(shp->bvh, pos, max_dist, early_exit,
            [shp](int eid, const ym::vec3f& pos, float max_dist,
                bool early_exit) {
                auto pt = intersection_point();
                if (!ym::overlap_point(
                        pos, max_dist, shp->pos[eid], shp->rad(eid), pt.dist))
                    return intersection_point{};
                pt.euv = {1, 0, 0, 0};
                pt.eid = eid;
                return pt;
            });
    };

    if (pt.eid >= 0) pt.sid = shp->sid;
    return pt;
}

//
// Shape overlap
//
intersection_point overlap_shape(const scene* scn, int sid,
    const ym::vec3f& pos, float max_dist, bool early_exit) {
    return overlap_shape(scn->shapes[sid], pos, max_dist, early_exit);
}

//
// Instance overlap
//
intersection_point overlap_instance(const instance* ist, const ym::vec3f& pos,
    float max_dist, bool early_exit) {
    auto pt = overlap_shape(ist->shp,
        ym::transform_point(ist->xform_inv, (ym::vec3f)pos), max_dist,
        early_exit);
    if (pt.eid >= 0) pt.iid = ist->iid;
    return pt;
}

//
// Instance overlap
//
intersection_point overlap_instance(const scene* scn, int iid,
    const ym::vec3f& pos, float max_dist, bool early_exit) {
    return overlap_instance(scn->instances[iid], pos, max_dist, early_exit);
}

//
// Scene overlap
//
intersection_point overlap_scene(
    const scene* scn, const ym::vec3f& pos, float max_dist, bool early_exit) {
    return overlap_bvh(scn->bvh, pos, max_dist, early_exit,
        [scn](int eid, const ym::vec3f& pos, float max_dist, bool early_exit) {
            return overlap_instance(
                scn->instances[eid], pos, max_dist, early_exit);
        });
}

// -----------------------------------------------------------------------------
// BVH CLOSEST ELEMENT LOOKUP FOR INTERNAL ELEMENTS (DISABLE FOR REFACTORING)
// -----------------------------------------------------------------------------

#if 0

//
// Shape elements overlap
//
// TODO: avoid duplicate elements
//
void overlap_elem(const shape* shp1, const shape* shp2, int idx1,
    int idx2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, ym::vec2i>>* overlaps,
    std::unordered_map<int, int>* closest) {
    // prepare point
    ym::vec4i verts;
    if (shp2->triangle) {
        verts = {shp2->triangle[idx2].x, shp2->triangle[idx2].y,
            shp2->triangle[idx2].z, -1};
    } else if (shp2->line) {
        verts = {shp2->line[idx2].x, shp2->line[idx2].y, -1, -1};
    } else if (shp2->point) {
        verts = {shp2->point[idx2], -1, -1, -1};
    } else if (shp2->tetra) {
        verts = shp2->tetra[idx2];
    } else {
        verts = {idx2, -1, -1, -1};
    }

    // loop over vertices
    for (auto vid : verts) {
        if (vid < 0) continue;

        // transform point
        auto pos = transform_point_inverse(
            shp1->frame, transform_point(shp2->frame, shp2->pos[vid]));
        auto rad = shp2->rad(vid) + radius;
        auto pt = point();

        // switch over shape type
        if (shp1->triangle) {
            auto f = shp1->triangle[idx1];
            if (!_overlap_triangle(pos, rad, shp1->pos[f.x], shp1->pos[f.y],
                    shp1->pos[f.z], shp1->rad(f.x), shp1->rad(f.y),
                    shp1->rad(f.z), pt.dist, (ym::vec3f&)pt.euv))
                return;
            pt.euv = {pt.euv.x, pt.euv.y, pt.euv.z, 0};
        } else if (shp1->line) {
            auto f = shp1->line[idx1];
            if (!_overlap_line(pos, rad, shp1->pos[f.x], shp1->pos[f.y],
                    shp1->rad(f.x), shp1->rad(f.y), pt.dist,
                    (ym::vec2f&)pt.euv))
                return;
            pt.euv = {pt.euv.x, pt.euv.y, 0, 0};
        } else if (shp1->point) {
            auto f = shp1->point[idx1];
            if (!_overlap_point(pos, rad, shp1->pos[f], shp1->rad(f), pt.dist))
                return;
            pt.euv = {1, 0, 0, 0};
        } else if (shp1->tetra) {
            auto f = shp1->tetra[idx1];
            if (!_overlap_tetrahedron(pos, rad, shp1->pos[f.x],
                    shp1->pos[f.y], shp1->pos[f.z], shp1->pos[f.w],
                    shp1->rad(f.x), shp1->rad(f.y), shp1->rad(f.z),
                    shp1->rad(f.w), pt.dist, (ym::vec4f&)pt.euv))
                return;
        } else {
            if (!_overlap_point(
                    pos, rad, shp1->pos[idx1], shp1->rad(idx1), pt.dist))
                return;
            pt.euv = {1, 0, 0, 0};
        }

        // wrap up
        pt.sid = shp1->sid;
        pt.eid = idx1;
        if (first_only) {
            if (closest->find(vid) == closest->end()) {
                overlaps->push_back({pt, {shp2->sid, vid}});
                (*closest)[vid] = (int)overlaps->size() - 1;
            } else {
                auto& overlap = (*overlaps)[(*closest)[vid]];
                if (overlap.first.dist > pt.dist) {
                    overlap = {pt, {shp2->sid, vid}};
                }
            }
        } else {
            overlaps->push_back({pt, {shp2->sid, vid}});
        }
    }
}

//
// Finds the closest element for all pairs of vertices.
// Similar to the generic public function whose
// interface is described above. See nightbor_bvh for parameter docs.
// With respect to that, only adds early_exit to decide whether we exit
// at the first primitive hit or we find the closest hit.
//
// Implementation Notes:
// - Walks the BVH using an internal stack to avoid the slowness of
// recursive
// calls; this follows general conventions and stragely makes the code
// shorter
// - The walk is simplified for first hit by obeserving that if we update
// the dist_max limit with the closest intersection distance during
// traversal, we will speed up computation significantly while simplifying
// the code; note in fact that all subsequent farthest iterations will be
// rejected in the tmax tests
//
void overlap_verts(const scene* scn1, const scene* scn2,
    int sid1, int sid2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, ym::vec2i>>* overlaps,
    std::unordered_map<int, int>* closest) {
    // get shape and bvh
    auto shp1 = (sid1 < 0) ? nullptr : scn1->shapes[sid1];
    auto bvh1 = (!shp1) ? scn1->bvh : shp1->bvh;
    auto shp2 = (sid2 < 0) ? nullptr : scn2->shapes[sid2];
    auto bvh2 = (!shp2) ? scn2->bvh : shp2->bvh;

    // node stack
    ym::vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // get frames
    auto frame1 = (!shp1) ? ym::identity_frame3f : shp1->frame;
    auto frame2 = (!shp2) ? ym::identity_frame3f : shp2->frame;

    // check if a trasformed test is needed
    auto xformed = frame1 != frame2;

    // walking stack
    while (node_cur) {
        // grab node
        const auto node_idx = node_stack[--node_cur];
        const auto node1 = bvh1->nodes[node_idx.x];
        const auto node2 = bvh2->nodes[node_idx.y];

        // intersect bbox
        if (xformed) {
            if (!_overlap_bbox(node1.bbox,
                    {node2.bbox.min - ym::vec3f{radius, radius, radius},
                        node2.bbox.max + ym::vec3f{radius, radius, radius}},
                    frame1, frame2))
                continue;
        } else {
            if (!_overlap_bbox(node1.bbox,
                    {node2.bbox.min - ym::vec3f{radius, radius, radius},
                        node2.bbox.max + ym::vec3f{radius, radius, radius}}))
                continue;
        }

        // check for leaves
        if (node1.isleaf && node2.isleaf) {
            if (!shp1) {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count;
                     i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (exclude_self && idx1 == idx2) continue;
                        _overlap_verts(scn1, scn2, idx1, idx2, exclude_self,
                            radius, first_only, overlaps, closest);
                    }
                }
            } else {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count;
                     i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (exclude_self && idx1 == idx2) continue;
                        _overlap_elem(shp1, shp2, idx1, idx2, exclude_self,
                            radius, first_only, overlaps, closest);
                    }
                }
            }
        } else {
            // descend
            if (node1.isleaf) {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    node_stack[node_cur++] = {node_idx.x, (int)idx2};
                    assert(node_cur < 128);
                }
            } else if (node2.isleaf) {
                for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                     idx1++) {
                    node_stack[node_cur++] = {(int)idx1, node_idx.y};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    for (auto idx1 = node1.start;
                         idx1 < node1.start + node1.count; idx1++) {
                        node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        assert(node_cur < 128);
                    }
                }
            }
        }
    }
}

//
// Find the list of overlaps between shapes.
// Public function whose interface is described above.
//
 void overlap_verts(const scene* scn1, const scene* scn2, int sid1,
    int sid2, bool exclude_self, float radius, bool first_only,
    std::vector<std::pair<point, ym::vec2i>>* overlaps) {
    std::unordered_map<int, int> closest;
    overlap_verts(
        scn1, scn2, sid1, sid2, false, radius, first_only, overlaps, &closest);
}

//
// Find the list of overlaps between scenes.
// Public function whose interface is described above.
//
 void overlap_verts(const scene* scn1, const scene* scn2,
    bool exclude_self, float radius, bool first_only,
    std::vector<ym::vec2i>* soverlaps,
    std::vector<std::pair<point, ym::vec2i>>* overlaps) {
    overlap_shape_bounds(scn1, scn2, false, false, exclude_self, soverlaps);
    for (auto sh : *soverlaps) {
        overlap_verts(scn1, scn2, sh.x, sh.y, exclude_self, radius,
            first_only, overlaps);
    }
}
#endif

//
// Finds the overlap between shape bounds.
// Similat interface as the public function.
//
void overlap_instance_bounds_(const scene* scn1, const scene* scn2,
    bool skip_duplicates, bool skip_self, std::vector<ym::vec2i>* overlaps) {
    // get bvhs
    auto bvh1 = scn1->bvh;
    auto bvh2 = scn2->bvh;

    // node stack
    ym::vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto node_idx = node_stack[--node_cur];
        const auto node1 = bvh1->nodes[node_idx.x];
        const auto node2 = bvh2->nodes[node_idx.y];

        // intersect bbox
        if (!ym::overlap_bbox(node1.bbox, node2.bbox)) continue;

        // check for leaves
        if (node1.isleaf && node2.isleaf) {
            // collide primitives
            for (auto i1 = node1.start; i1 < node1.start + node1.count; i1++) {
                for (auto i2 = node2.start; i2 < node2.start + node2.count;
                     i2++) {
                    auto idx1 = bvh1->sorted_prim[i1];
                    auto idx2 = bvh2->sorted_prim[i2];
                    auto ist1 = scn1->instances[idx1];
                    auto ist2 = scn2->instances[idx2];
                    if (skip_duplicates && ist1->iid > ist2->iid) continue;
                    if (skip_self && ist1->iid == ist2->iid) continue;
                    if (ym::overlap_bbox(ist1->bbox, ist2->bbox))
                        overlaps->push_back({ist1->iid, ist2->iid});
                }
            }
        } else {
            // descend
            if (node1.isleaf) {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    node_stack[node_cur++] = {node_idx.x, (int)idx2};
                    assert(node_cur < 128);
                }
            } else if (node2.isleaf) {
                for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                     idx1++) {
                    node_stack[node_cur++] = {(int)idx1, node_idx.y};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    for (auto idx1 = node1.start;
                         idx1 < node1.start + node1.count; idx1++) {
                        node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        assert(node_cur < 128);
                    }
                }
            }
        }
    }
}

//
// Find the list of overlaps between shape bounds.
// Public function whose interface is described above.
//
void overlap_instance_bounds(const scene* scn1, const scene* scn2,
    bool skip_duplicates, bool skip_self, std::vector<ym::vec2i>* overlaps) {
    overlaps->clear();
    overlap_instance_bounds_(scn1, scn2, skip_duplicates, skip_self, overlaps);
}

// -----------------------------------------------------------------------------
// STATISTICS FOR DEBUGGING (probably not helpful to all)
// -----------------------------------------------------------------------------

//
// Compute BVH stats.
//
void compute_bvh_stats(const scene* scn, int shape_id, bool include_shapes,
    int depth, int& nprims, int& ninternals, int& nleaves, int& min_depth,
    int& max_depth) {
    // get bvh
    auto bvh = (shape_id >= 0) ? scn->shapes[shape_id]->bvh : scn->bvh;

    // node stack
    ym::vec2i node_stack[128];  // node and depth
    auto node_cur = 0;
    node_stack[node_cur++] = ym::vec2i{0, depth};

    // walk the stack
    while (node_cur) {
        // get node and depth
        auto node_depth = node_stack[--node_cur];
        auto node = &bvh->nodes[node_depth[0]];

        // update stats
        if (!node->isleaf) {
            ninternals += 1;
            for (auto i = 0; i < node->count; i++) {
                node_stack[node_cur++] =
                    ym::vec2i{(int)(node->start + i), node_depth.y + 1};
            }
        } else if (shape_id >= 0) {
            nleaves += 1;
            nprims += node->count;
            min_depth = ym::min(min_depth, node_depth.y);
            max_depth = ym::max(max_depth, node_depth.y);
        } else {
            if (include_shapes) {
                for (auto i = 0; i < node->count; i++) {
                    auto idx = bvh->sorted_prim[node->start + i];
                    compute_bvh_stats(scn, idx, true, node_depth.y + 1, nprims,
                        ninternals, nleaves, min_depth, max_depth);
                }
            } else {
                nleaves += 1;
                nprims += node->count;
                min_depth = ym::min(min_depth, node_depth.y);
                max_depth = ym::max(max_depth, node_depth.y);
            }
        }
    }
}

//
// Compute BVH stats.
//
void compute_bvh_stats(const scene* scn, bool include_shapes, int& nprims,
    int& ninternals, int& nleaves, int& min_depth, int& max_depth,
    int req_shape) {
    // init out variables
    nprims = 0;
    ninternals = 0;
    nleaves = 0;
    min_depth = 10000;
    max_depth = -1;

    compute_bvh_stats(scn, req_shape, include_shapes, 0, nprims, ninternals,
        nleaves, min_depth, max_depth);
}

}  // namespace ybvh
