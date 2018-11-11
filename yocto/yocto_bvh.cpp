//
// Implementation for Yocto/BVH.
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

#include "yocto_bvh.h"
#include "yocto_utils.h"

#if YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect a ray with a point (approximate)
bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& distance, vec2f& uv) {
    // find parameter for line-point minimum distance
    auto w = p - ray.origin;
    auto t = dot(w, ray.direction) / dot(ray.direction, ray.direction);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp  = ray.origin + ray.direction * t;
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    // intersection occurred: set options and exit
    distance = t;
    uv       = {0, 0};

    return true;
}

// Intersect a ray with a line
bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, float& distance, vec2f& uv) {
    // setup intersection options
    auto u = ray.direction;
    auto v = p1 - p0;
    auto w = ray.origin - p0;

    // compute values to solve a linear system
    auto a   = dot(u, u);
    auto b   = dot(u, v);
    auto c   = dot(v, v);
    auto d   = dot(u, w);
    auto e   = dot(v, w);
    auto det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute Parameters on both ray and segment
    auto t = (b * e - c * d) / det;
    auto s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // clamp segment param to segment corners
    s = clamp(s, (float)0, (float)1);

    // compute segment-segment distance on the closest points
    auto pr  = ray.origin + ray.direction * t;
    auto pl  = p0 + (p1 - p0) * s;
    auto prl = pr - pl;

    // check with the line radius at the same point
    auto d2 = dot(prl, prl);
    auto r  = r0 * (1 - s) + r1 * s;
    if (d2 > r * r) return false;

    // intersection occurred: set options and exit
    distance = t;
    uv       = {s, sqrt(d2) / r};

    return true;
}

// Intersect a ray with a triangle
bool intersect_triangle(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, float& distance, vec2f& uv) {
    // compute triangle edges
    auto edge1 = p1 - p0;
    auto edge2 = p2 - p0;

    // compute determinant to solve a linear system
    auto pvec = cross(ray.direction, edge2);
    auto det  = dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.origin - p0;
    auto u    = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = cross(tvec, edge1);
    auto v    = dot(ray.direction, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set options and exit
    distance = t;
    uv       = {u, v};

    return true;
}

// Intersect a ray with a quad.
bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, float& distance, vec2f& uv) {
    auto hit  = false;
    auto tray = ray;
    if (intersect_triangle(tray, p0, p1, p3, distance, uv)) {
        tray.tmax = distance;
        hit       = true;
    }
    if (intersect_triangle(tray, p2, p3, p1, distance, uv)) {
        uv        = {1 - uv[0], 1 - uv[1]};
        tray.tmax = distance;
        hit       = true;
    }
    return hit;
}

// Min/max used in BVH traversal. Copied here since the traversal code
// relies on the specific behaviour wrt NaNs.
static inline const float& _safemin(const float& a, const float& b) {
    return (a < b) ? a : b;
}
// Min/max used in BVH traversal. Copied here since the traversal code
// relies on the specific behaviour wrt NaNs.
static inline const float& _safemax(const float& a, const float& b) {
    return (a > b) ? a : b;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox) {
    // determine intersection ranges
    auto invd = 1.0f / ray.direction;
    auto t0   = (bbox.min - ray.origin) * invd;
    auto t1   = (bbox.max - ray.origin) * invd;
    // flip based on range directions
    if (invd[0] < 0.0f) swap(t0[0], t1[0]);
    if (invd[1] < 0.0f) swap(t0[1], t1[1]);
    if (invd[2] < 0.0f) swap(t0[2], t1[2]);
    auto tmin = _safemax(t0[2], _safemax(t0[1], _safemax(t0[0], ray.tmin)));
    auto tmax = _safemin(t1[2], _safemin(t1[1], _safemin(t1[0], ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox) {
    auto txmin = (bbox[ray_dsign[0]][0] - ray.origin[0]) * ray_dinv[0];
    auto txmax = (bbox[1 - ray_dsign[0]][0] - ray.origin[0]) * ray_dinv[0];
    auto tymin = (bbox[ray_dsign[1]][1] - ray.origin[1]) * ray_dinv[1];
    auto tymax = (bbox[1 - ray_dsign[1]][1] - ray.origin[1]) * ray_dinv[1];
    auto tzmin = (bbox[ray_dsign[2]][2] - ray.origin[2]) * ray_dinv[2];
    auto tzmax = (bbox[1 - ray_dsign[2]][2] - ray.origin[2]) * ray_dinv[2];
    auto tmin  = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax  = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// TODO: documentation
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p, float r,
    float& distance, vec2f& uv) {
    auto d2 = dot(pos - p, pos - p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    distance = sqrt(d2);
    uv       = {0, 0};
    return true;
}

// TODO: documentation
float closestuv_line(const vec3f& pos, const vec3f& p0, const vec3f& p1) {
    auto ab = p1 - p0;
    auto d  = dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b –
    // a)
    auto u = dot(pos - p0, ab) / d;
    u      = clamp(u, (float)0, (float)1);
    return u;
}

// TODO: documentation
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, float r0, float r1, float& distance, vec2f& uv) {
    auto u = closestuv_line(pos, p0, p1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p  = p0 + (p1 - p0) * u;
    auto r  = r0 + (r1 - r0) * u;
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    distance = sqrt(d2);
    uv       = {u, 0};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably "--"+prefix to use a sequence of
// test (triangle body, and 3 edges)
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    auto ab = p1 - p0;
    auto ac = p2 - p0;
    auto ap = pos - p0;

    auto d1 = dot(ab, ap);
    auto d2 = dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return {0, 0};

    auto bp = pos - p1;
    auto d3 = dot(ab, bp);
    auto d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return {1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return {d1 / (d1 - d3), 0};

    auto cp = pos - p2;
    auto d5 = dot(ab, cp);
    auto d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return {0, 1};

    auto vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) return {0, d2 / (d2 - d6)};

    auto va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {1 - w, w};
    }

    // face case
    auto denom = 1 / (va + vb + vc);
    auto u     = vb * denom;
    auto v     = vc * denom;
    return {u, v};
}

// TODO: documentation
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, float r0, float r1, float r2,
    float& distance, vec2f& uv) {
    uv      = closestuv_triangle(pos, p0, p1, p2);
    auto p  = p0 * (1 - uv[0] - uv[1]) + p1 * uv[0] + p2 * uv[1];
    auto r  = r0 * (1 - uv[0] - uv[1]) + r1 * uv[0] + r2 * uv[1];
    auto dd = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    distance = sqrt(dd);
    return true;
}

// TODO: documentation
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
    float r2, float r3, float& distance, vec2f& uv) {
    auto hit = false;
    if (overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r3, distance, uv)) {
        dist_max = distance;
        hit      = true;
    }
    if (overlap_triangle(pos, dist_max, p2, p3, p1, r2, r3, r1, distance, uv)) {
        // dist_max = distance;
        uv  = {1 - uv[0], 1 - uv[1]};
        hit = true;
    }
    return hit;
}

// TODO: documentation
bool distance_check_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox) {
    // computing distance
    auto dd = 0.0f;

    // For each axis count any excess distance outside box extents
    if (pos[0] < bbox.min[0])
        dd += (bbox.min[0] - pos[0]) * (bbox.min[0] - pos[0]);
    if (pos[0] > bbox.max[0])
        dd += (pos[0] - bbox.max[0]) * (pos[0] - bbox.max[0]);
    if (pos[1] < bbox.min[1])
        dd += (bbox.min[1] - pos[1]) * (bbox.min[1] - pos[1]);
    if (pos[1] > bbox.max[1])
        dd += (pos[1] - bbox.max[1]) * (pos[1] - bbox.max[1]);
    if (pos[2] < bbox.min[2])
        dd += (bbox.min[2] - pos[2]) * (bbox.min[2] - pos[2]);
    if (pos[2] > bbox.max[2])
        dd += (pos[2] - bbox.max[2]) * (pos[2] - bbox.max[2]);

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
    if (bbox1.max[0] < bbox2.min[0] || bbox1.min[0] > bbox2.max[0])
        return false;
    if (bbox1.max[1] < bbox2.min[1] || bbox1.min[1] > bbox2.max[1])
        return false;
    if (bbox1.max[2] < bbox2.min[2] || bbox1.min[2] > bbox2.max[2])
        return false;
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH
// -----------------------------------------------------------------------------
namespace yocto {

// Cleanup
void clear_shape_bvh_embree(bvh_shape& bvh) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) {
        rtcReleaseScene((RTCScene)bvh.embree_bvh);
    }
#endif
}
void clear_scene_bvh_embree(bvh_scene& bvh) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) {
        for (auto i = 0; i < max(1, (int)bvh.instances.size()); i++) {
            auto geom = rtcGetGeometry((RTCScene)bvh.embree_bvh, i);
            rtcDetachGeometry((RTCScene)bvh.embree_bvh, i);
            rtcReleaseGeometry(geom);
        }
        rtcReleaseScene((RTCScene)bvh.embree_bvh);
    }
#endif
}

#if YOCTO_EMBREE
void embree_error(void* ctx, RTCError code, const char* str) {
    switch (code) {
        case RTC_ERROR_UNKNOWN: printf("RTC_ERROR_UNKNOWN"); break;
        case RTC_ERROR_INVALID_ARGUMENT:
            printf("RTC_ERROR_INVALID_ARGUMENT");
            break;
        case RTC_ERROR_INVALID_OPERATION:
            printf("RTC_ERROR_INVALID_OPERATION");
            break;
        case RTC_ERROR_OUT_OF_MEMORY: printf("RTC_ERROR_OUT_OF_MEMORY"); break;
        case RTC_ERROR_UNSUPPORTED_CPU:
            printf("RTC_ERROR_UNSUPPORTED_CPU");
            break;
        case RTC_ERROR_CANCELLED: printf("RTC_ERROR_CANCELLED"); break;
        default: printf("invalid error code"); break;
    }
    printf(": %s\n", str);
}

// Get Embree device
RTCDevice get_embree_device() {
    static RTCDevice device = nullptr;
    if (!device) {
        device = rtcNewDevice("");
        rtcSetDeviceErrorFunction(device, embree_error, nullptr);
    }
    return device;
}

// Build a BVH using Embree. Calls `build_scene_bvh()` if Embree is not available.
void build_embree_bvh(bvh_shape& bvh) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if (!bvh.points.empty()) {
        log_error("embree does not support points");
    } else if (!bvh.lines.empty()) {
        log_error("not yet implemented");
    } else if (!bvh.triangles.empty()) {
        auto embree_geom = rtcNewGeometry(
            embree_device, RTC_GEOMETRY_TYPE_TRIANGLE);
        rtcSetGeometryVertexAttributeCount(embree_geom, 1);
        auto vert = rtcSetNewGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX,
            0, RTC_FORMAT_FLOAT3, 3 * 4, bvh.positions.size());
        auto triangles = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
            bvh.triangles.size());
        memcpy(vert, bvh.positions.data(), bvh.positions.size() * 12);
        memcpy(triangles, bvh.triangles.data(), bvh.triangles.size() * 12);
        rtcCommitGeometry(embree_geom);
        rtcAttachGeometryByID(embree_scene, embree_geom, 0);
    } else if (!bvh.quads.empty()) {
        auto embree_geom = rtcNewGeometry(embree_device, RTC_GEOMETRY_TYPE_QUAD);
        rtcSetGeometryVertexAttributeCount(embree_geom, 1);
        auto vert = rtcSetNewGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX,
            0, RTC_FORMAT_FLOAT3, 3 * 4, bvh.positions.size());
        auto quads = rtcSetNewGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_INDEX,
            0, RTC_FORMAT_UINT4, 4 * 4, bvh.quads.size());
        memcpy(vert, bvh.positions.data(), bvh.positions.size() * 12);
        memcpy(quads, bvh.quads.data(), bvh.quads.size() * 16);
        rtcCommitGeometry(embree_geom);
        rtcAttachGeometryByID(embree_scene, embree_geom, 0);
    } else {
        log_error("empty bvh");
    }
    rtcCommitScene(embree_scene);
    bvh.embree_bvh = embree_scene;
}
void build_embree_bvh(bvh_scene& bvh) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if (!bvh.instances.empty()) {
        for (auto instance_id = 0; instance_id < bvh.instances.size();
             instance_id++) {
            auto& instance = bvh.instances[instance_id];
            if (instance.shape_id < 0 && instance.surface_id < 0) continue;
            if (instance.shape_id >= 0 &&
                bvh.shape_bvhs[instance.shape_id].positions.empty())
                continue;
            if (instance.surface_id >= 0 &&
                bvh.surface_bvhs[instance.surface_id].positions.empty())
                continue;
            auto embree_geom = rtcNewGeometry(
                embree_device, RTC_GEOMETRY_TYPE_INSTANCE);
            if (instance.shape_id >= 0) {
                rtcSetGeometryInstancedScene(embree_geom,
                    (RTCScene)bvh.shape_bvhs[instance.shape_id].embree_bvh);
            } else if (instance.surface_id >= 0) {
                rtcSetGeometryInstancedScene(embree_geom,
                    (RTCScene)bvh.surface_bvhs[instance.surface_id].embree_bvh);
            } else {
                log_error("empty instance");
            }
            rtcSetGeometryTransform(embree_geom, 0,
                RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
            rtcCommitGeometry(embree_geom);
            rtcAttachGeometryByID(embree_scene, embree_geom, instance_id);
        }
    }
    rtcCommitScene(embree_scene);
    bvh.embree_bvh = embree_scene;
}
// Refit a BVH using Embree. Calls `refit_scene_bvh()` if Embree is not available.
void refit_embree_bvh(bvh_shape& bvh) { log_error("not yet implemented"); }
void refit_embree_bvh(bvh_scene& bvh) { log_error("not yet implemented"); }
bool intersect_embree_bvh(const bvh_shape& bvh, const ray3f& ray, bool find_any,
    float& distance, int& element_id, vec2f& uv) {
    RTCRayHit embree_ray;
    embree_ray.ray.org_x     = ray.origin[0];
    embree_ray.ray.org_y     = ray.origin[1];
    embree_ray.ray.org_z     = ray.origin[2];
    embree_ray.ray.dir_x     = ray.direction[0];
    embree_ray.ray.dir_y     = ray.direction[1];
    embree_ray.ray.dir_z     = ray.direction[2];
    embree_ray.ray.tnear     = ray.tmin;
    embree_ray.ray.tfar      = ray.tmax;
    embree_ray.ray.flags     = 0;
    embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
    embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    RTCIntersectContext embree_ctx;
    rtcInitIntersectContext(&embree_ctx);
    rtcIntersect1((RTCScene)bvh.embree_bvh, &embree_ctx, &embree_ray);
    if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
    distance   = embree_ray.ray.tfar;
    uv         = {embree_ray.hit.u, embree_ray.hit.v};
    element_id = embree_ray.hit.primID;
    return true;
}
bool intersect_embree_bvh(const bvh_scene& bvh, const ray3f& ray, bool find_any,
    float& distance, int& instance_id, int& element_id, vec2f& uv) {
    RTCRayHit embree_ray;
    embree_ray.ray.org_x     = ray.origin[0];
    embree_ray.ray.org_y     = ray.origin[1];
    embree_ray.ray.org_z     = ray.origin[2];
    embree_ray.ray.dir_x     = ray.direction[0];
    embree_ray.ray.dir_y     = ray.direction[1];
    embree_ray.ray.dir_z     = ray.direction[2];
    embree_ray.ray.tnear     = ray.tmin;
    embree_ray.ray.tfar      = ray.tmax;
    embree_ray.ray.flags     = 0;
    embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
    embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    RTCIntersectContext embree_ctx;
    rtcInitIntersectContext(&embree_ctx);
    rtcIntersect1((RTCScene)bvh.embree_bvh, &embree_ctx, &embree_ray);
    if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
    distance    = embree_ray.ray.tfar;
    uv          = {embree_ray.hit.u, embree_ray.hit.v};
    element_id  = embree_ray.hit.primID;
    instance_id = embree_ray.hit.instID[0];
    return true;
}
#endif

// BVH primitive with its bbox, its center and the index to the primitive
struct bvh_prim {
    bbox3f bbox   = invalid_bbox3f;
    vec3f  center = zero_vec3f;
    int    primid = 0;
};

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
pair<int, int> split_bvh_node_sah(vector<bvh_prim>& prims, int start, int end) {
    // initialize split axis and position
    auto split_axis = 0;
    auto mid        = (start + end) / 2;

    // compute primintive bounds and size
    auto cbbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) cbbox += prims[i].center;
    auto csize = cbbox.max - cbbox.min;
    if (csize == zero_vec3f) return {mid, split_axis};

    // consider N bins, compute their cost and keep the minimum
    const int nbins     = 16;
    auto      middle    = 0.0f;
    auto      min_cost  = float_max;
    auto      bbox_area = [](auto& b) {
        auto size = b.max - b.min;
        return 1e-12f + 2 * size[0] * size[1] + 2 * size[0] * size[2] +
               2 * size[1] * size[2];
    };
    auto min_left = 0, min_right = 0;
    for (auto axis = 0; axis < 3; axis++) {
        for (auto b = 1; b < nbins; b++) {
            auto split     = (&cbbox.min[0])[axis] + b * csize[axis] / nbins;
            auto left_bbox = invalid_bbox3f, right_bbox = invalid_bbox3f;
            auto left_nprims = 0, right_nprims = 0;
            for (auto i = start; i < end; i++) {
                if ((&prims[i].center[0])[axis] < split) {
                    left_bbox += prims[i].bbox;
                    left_nprims += 1;
                } else {
                    right_bbox += prims[i].bbox;
                    right_nprims += 1;
                }
            }
            auto cost = 1 +
                        left_nprims * bbox_area(left_bbox) / bbox_area(cbbox) +
                        right_nprims * bbox_area(right_bbox) / bbox_area(cbbox);
            if (cost < min_cost) {
                min_cost   = cost;
                middle     = split;
                split_axis = axis;
                min_left   = left_nprims;
                min_right  = right_nprims;
            }
        }
    }
    // split
    mid = (int)(std::partition(prims.data() + start, prims.data() + end,
                    [split_axis, middle](
                        auto& a) { return a.center[split_axis] < middle; }) -
                prims.data());

    // if we were not able to split, just break the primitives in half
    if (mid == start || mid == end) {
        log_error("bad bvh split");
        split_axis = 0;
        mid        = (start + end) / 2;
    }

    return {mid, split_axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and axis.
pair<int, int> split_bvh_node_balanced(
    vector<bvh_prim>& prims, int start, int end) {
    // initialize split axis and position
    auto split_axis = 0;
    auto mid        = (start + end) / 2;

    // compute primintive bounds and size
    auto cbbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) cbbox += prims[i].center;
    auto csize = cbbox.max - cbbox.min;
    if (csize == zero_vec3f) return {mid, split_axis};

    // split along largest
    auto largest_axis = 0;
    if (csize[0] >= csize[1] && csize[0] >= csize[2]) largest_axis = 0;
    if (csize[1] >= csize[0] && csize[1] >= csize[2]) largest_axis = 1;
    if (csize[2] >= csize[0] && csize[2] >= csize[1]) largest_axis = 2;

    // balanced tree split: find the largest axis of the
    // bounding box and split along this one right in the middle
    split_axis = largest_axis;
    mid        = (start + end) / 2;
    std::nth_element(prims.data() + start, prims.data() + mid,
        prims.data() + end, [split_axis](auto& a, auto& b) {
            return a.center[split_axis] < b.center[split_axis];
        });

    // if we were not able to split, just break the primitives in half
    if (mid == start || mid == end) {
        log_error("bad bvh split");
        split_axis = 0;
        mid        = (start + end) / 2;
    }

    return {mid, split_axis};
}

// Splits a BVH node using the middle heutirtic. Returns split position and axis.
pair<int, int> split_bvh_node_middle(vector<bvh_prim>& prims, int start, int end) {
    // initialize split axis and position
    auto split_axis = 0;
    auto mid        = (start + end) / 2;

    // compute primintive bounds and size
    auto cbbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) cbbox += prims[i].center;
    auto csize = cbbox.max - cbbox.min;
    if (csize == zero_vec3f) return {mid, split_axis};

    // split along largest
    auto largest_axis = 0;
    if (csize[0] >= csize[1] && csize[0] >= csize[2]) largest_axis = 0;
    if (csize[1] >= csize[0] && csize[1] >= csize[2]) largest_axis = 1;
    if (csize[2] >= csize[0] && csize[2] >= csize[1]) largest_axis = 2;

    // split the space in the middle along the largest axis
    split_axis   = largest_axis;
    auto cmiddle = (cbbox.max + cbbox.min) / 2;
    auto middle  = cmiddle[largest_axis];
    mid = (int)(std::partition(prims.data() + start, prims.data() + end,
                    [split_axis, middle](
                        auto& a) { return a.center[split_axis] < middle; }) -
                prims.data());

    // if we were not able to split, just break the primitives in half
    if (mid == start || mid == end) {
        log_error("bad bvh split");
        split_axis = 0;
        mid        = (start + end) / 2;
    }

    return {mid, split_axis};
}

// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
void make_bvh_node(vector<bvh_node>& nodes, vector<bvh_prim>& prims,
    deque<vec3i>& queue, bool high_quality) {
    // grab node to work on
    auto next = queue.front();
    queue.pop_front();
    auto nodeid = next[0], start = next[1], end = next[2];

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

    // split into two children
    if (end - start > bvh_max_prims) {
        // get split
        auto split = (high_quality) ? split_bvh_node_sah(prims, start, end) :
                                      split_bvh_node_balanced(prims, start, end);
        auto mid = split.first, split_axis = split.second;

        // make an internal node
        node.is_internal      = true;
        node.split_axis       = split_axis;
        node.num_primitives   = 2;
        node.primitive_ids[0] = (int)nodes.size() + 0;
        node.primitive_ids[1] = (int)nodes.size() + 1;
        nodes.emplace_back();
        nodes.emplace_back();
        queue.push_back({node.primitive_ids[0], start, mid});
        queue.push_back({node.primitive_ids[1], mid, end});
    } else {
        // Make a leaf node
        node.is_internal    = false;
        node.num_primitives = end - start;
        for (auto i = 0; i < node.num_primitives; i++)
            node.primitive_ids[i] = prims[start + i].primid;
    }
}

// Build BVH nodes
void build_bvh_nodes_serial(vector<bvh_node>& nodes, vector<bvh_prim>& prims,
    const build_bvh_options& options) {
    // prepare to build nodes
    nodes.clear();
    nodes.reserve(prims.size() * 2);

    // queue up first node
    auto queue = deque<vec3i>{{0, 0, (int)prims.size()}};
    nodes.emplace_back();

    // create nodes until the queue is empty
    while (!queue.empty()) {
        // exit if needed
        if (options.cancel_flag && *options.cancel_flag) return;

        // grab node to work on
        auto next = queue.front();
        queue.pop_front();
        auto nodeid = next[0], start = next[1], end = next[2];

        // grab node
        auto& node = nodes[nodeid];

        // compute bounds
        node.bbox = invalid_bbox3f;
        for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

        // split into two children
        if (end - start > bvh_max_prims) {
            // get split
            auto split = (options.high_quality) ?
                             split_bvh_node_sah(prims, start, end) :
                             split_bvh_node_balanced(prims, start, end);
            auto mid = split.first, split_axis = split.second;

            // make an internal node
            node.is_internal      = true;
            node.split_axis       = split_axis;
            node.num_primitives   = 2;
            node.primitive_ids[0] = (int)nodes.size() + 0;
            node.primitive_ids[1] = (int)nodes.size() + 1;
            nodes.emplace_back();
            nodes.emplace_back();
            queue.push_back({node.primitive_ids[0], start, mid});
            queue.push_back({node.primitive_ids[1], mid, end});
        } else {
            // Make a leaf node
            node.is_internal    = false;
            node.num_primitives = end - start;
            for (auto i = 0; i < node.num_primitives; i++)
                node.primitive_ids[i] = prims[start + i].primid;
        }
    }

    // cleanup
    nodes.shrink_to_fit();
}

// Build BVH nodes
void build_bvh_nodes_parallel(vector<bvh_node>& nodes, vector<bvh_prim>& prims,
    const build_bvh_options& options) {
    // prepare to build nodes
    nodes.clear();
    nodes.reserve(prims.size() * 2);

    // queue up first node
    auto queue = deque<vec3i>{{0, 0, (int)prims.size()}};
    nodes.emplace_back();

    // synchronization
    atomic<int>    num_processed_prims(0);
    mutex          queue_mutex;
    vector<thread> threads;
    auto           nthreads = thread::hardware_concurrency();

    // create nodes until the queue is empty
    for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
        threads.emplace_back([&nodes, &prims, &options, &num_processed_prims,
                                 &queue_mutex, &queue] {
            while (true) {
                // exit if needed
                if (num_processed_prims >= prims.size()) return;
                if (options.cancel_flag && *options.cancel_flag) return;

                // grab node to work on
                auto next = zero_vec3i;
                {
                    lock_guard<mutex> lock{queue_mutex};
                    if (!queue.empty()) {
                        next = queue.front();
                        queue.pop_front();
                    }
                }

                // wait a bit if needed
                if (next == zero_vec3i) {
                    std::this_thread::sleep_for(10us);
                    continue;
                }

                // grab node
                auto  nodeid = next[0], start = next[1], end = next[2];
                auto& node = nodes[nodeid];

                // compute bounds
                node.bbox = invalid_bbox3f;
                for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

                // split into two children
                if (end - start > bvh_max_prims) {
                    // get split
                    auto split = (options.high_quality) ?
                                     split_bvh_node_sah(prims, start, end) :
                                     split_bvh_node_balanced(prims, start, end);
                    auto mid = split.first, split_axis = split.second;

                    // make an internal node
                    {
                        lock_guard<mutex> lock{queue_mutex};
                        node.is_internal      = true;
                        node.split_axis       = split_axis;
                        node.num_primitives   = 2;
                        node.primitive_ids[0] = (int)nodes.size() + 0;
                        node.primitive_ids[1] = (int)nodes.size() + 1;
                        nodes.emplace_back();
                        nodes.emplace_back();
                        queue.push_back({node.primitive_ids[0], start, mid});
                        queue.push_back({node.primitive_ids[1], mid, end});
                    }
                } else {
                    // Make a leaf node
                    node.is_internal    = false;
                    node.num_primitives = end - start;
                    for (auto i = 0; i < node.num_primitives; i++)
                        node.primitive_ids[i] = prims[start + i].primid;
                    num_processed_prims += node.num_primitives;
                }
            }
        });
    }
    for (auto& t : threads) t.join();

    // cleanup
    nodes.shrink_to_fit();
}

// Build a BVH from a set of primitives.
void build_shape_bvh(bvh_shape& bvh, const build_bvh_options& options) {
#if YOCTO_EMBREE
    if (options.use_embree) return build_embree_bvh(bvh);
#endif

    // get the number of primitives and the primitive type
    auto prims = vector<bvh_prim>();
    if (!bvh.points.empty()) {
        for (auto& p : bvh.points) {
            prims.push_back({point_bounds(bvh.positions[p], bvh.radius[p])});
        }
    } else if (!bvh.lines.empty()) {
        for (auto& l : bvh.lines) {
            prims.push_back({line_bounds(bvh.positions[l[0]],
                bvh.positions[l[1]], bvh.radius[l[0]], bvh.radius[l[1]])});
        }
    } else if (!bvh.triangles.empty()) {
        for (auto& t : bvh.triangles) {
            prims.push_back({triangle_bounds(bvh.positions[t[0]],
                bvh.positions[t[1]], bvh.positions[t[2]])});
        }
    } else if (!bvh.quads.empty()) {
        for (auto& q : bvh.quads) {
            prims.push_back({quad_bounds(bvh.positions[q[0]],
                bvh.positions[q[1]], bvh.positions[q[2]], bvh.positions[q[3]])});
        }
    }

    // create an array of primitives to sort
    for (auto i = 0; i < prims.size(); i++) {
        prims[i].center = (prims[i].bbox.min + prims[i].bbox.max) / 2;
        prims[i].primid = i;
    }

    // build nodes
    if (options.run_serially) {
        build_bvh_nodes_serial(bvh.nodes, prims, options);
    } else {
        build_bvh_nodes_parallel(bvh.nodes, prims, options);
    }
}

// Build a BVH from the given set of shape primitives.
bvh_shape make_shape_bvh(const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const build_bvh_options& options) {
    auto bvh      = bvh_shape{};
    bvh.points    = points;
    bvh.positions = positions;
    bvh.radius    = radius;
    build_shape_bvh(bvh, options);
    return bvh;
}
bvh_shape make_shape_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const build_bvh_options& options) {
    auto bvh      = bvh_shape{};
    bvh.lines     = lines;
    bvh.positions = positions;
    bvh.radius    = radius;
    build_shape_bvh(bvh, options);
    return bvh;
}
bvh_shape make_shape_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const build_bvh_options& options) {
    auto bvh      = bvh_shape{};
    bvh.triangles = triangles;
    bvh.positions = positions;
    build_shape_bvh(bvh, options);
    return bvh;
}
bvh_shape make_shape_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& positions, const build_bvh_options& options) {
    auto bvh      = bvh_shape{};
    bvh.quads     = quads;
    bvh.positions = positions;
    build_shape_bvh(bvh, options);
    return bvh;
}

// Build a BVH from a set of primitives.
void build_scene_bvh(bvh_scene& bvh, const build_bvh_options& options) {
#if YOCTO_EMBREE
    if (options.use_embree) return build_embree_bvh(bvh);
#endif

    // get the number of primitives and the primitive type
    auto prims = vector<bvh_prim>();
    if (!bvh.instances.empty()) {
        for (auto& instance : bvh.instances) {
            if (instance.shape_id >= 0) {
                auto& sbvh = bvh.shape_bvhs[instance.shape_id];
                if (sbvh.nodes.empty()) {
                    prims.push_back({invalid_bbox3f});
                } else {
                    prims.push_back(
                        {transform_bbox(instance.frame, sbvh.nodes[0].bbox)});
                }
            } else if (instance.surface_id >= 0) {
                auto& sbvh = bvh.surface_bvhs[instance.surface_id];
                if (sbvh.nodes.empty()) {
                    prims.push_back({invalid_bbox3f});
                } else {
                    prims.push_back(
                        {transform_bbox(instance.frame, sbvh.nodes[0].bbox)});
                }
            } else {
                log_warning("empty instance");
                prims.push_back({invalid_bbox3f});
            }
        }
    }

    // create an array of primitives to sort
    for (auto i = 0; i < prims.size(); i++) {
        prims[i].center = (prims[i].bbox.min + prims[i].bbox.max) / 2;
        prims[i].primid = i;
    }

    // build nodes
    if (options.run_serially) {
        build_bvh_nodes_serial(bvh.nodes, prims, options);
    } else {
        build_bvh_nodes_parallel(bvh.nodes, prims, options);
    }
}

// Build a BVH from the given set of instances.
bvh_scene make_scene_bvh(const vector<bvh_instance>& instances,
    const vector<bvh_shape>& shape_bvhs, const vector<bvh_shape>& surface_bvhs,
    const build_bvh_options& options) {
    auto bvh         = bvh_scene{};
    bvh.instances    = instances;
    bvh.shape_bvhs   = shape_bvhs;
    bvh.surface_bvhs = surface_bvhs;
    build_scene_bvh(bvh, options);
    return bvh;
}

// Recursively recomputes the node bounds for a shape bvh
void refit_shape_bvh(bvh_shape& bvh, int nodeid) {
    // refit
    auto& node = bvh.nodes[nodeid];
    node.bbox  = invalid_bbox3f;
    if (node.is_internal) {
        for (auto i = 0; i < 2; i++) {
            refit_shape_bvh(bvh, node.primitive_ids[i]);
            node.bbox += bvh.nodes[node.primitive_ids[i]].bbox;
        }
    } else if (!bvh.triangles.empty()) {
        for (auto i = 0; i < node.num_primitives; i++) {
            auto& t = bvh.triangles[node.primitive_ids[i]];
            node.bbox += triangle_bounds(
                bvh.positions[t[0]], bvh.positions[t[1]], bvh.positions[t[2]]);
        }
    } else if (!bvh.quads.empty()) {
        for (auto i = 0; i < node.num_primitives; i++) {
            auto& t = bvh.quads[node.primitive_ids[i]];
            node.bbox += quad_bounds(bvh.positions[t[0]], bvh.positions[t[1]],
                bvh.positions[t[2]], bvh.positions[t[3]]);
        }
    } else if (!bvh.lines.empty()) {
        for (auto i = 0; i < node.num_primitives; i++) {
            auto& l = bvh.lines[node.primitive_ids[i]];
            node.bbox += line_bounds(bvh.positions[l[0]], bvh.positions[l[1]],
                bvh.radius[l[0]], bvh.radius[l[1]]);
        }
    } else if (!bvh.points.empty()) {
        for (auto i = 0; i < node.num_primitives; i++) {
            auto& p = bvh.points[node.primitive_ids[i]];
            node.bbox += point_bounds(bvh.positions[p], bvh.radius[p]);
        }
    } else {
        log_warning("empty bvh");
    }
}

// Recursively recomputes the node bounds for a shape bvh
void refit_scene_bvh(bvh_scene& bvh, int nodeid) {
    // refit
    auto& node = bvh.nodes[nodeid];
    node.bbox  = invalid_bbox3f;
    if (node.is_internal) {
        for (auto i = 0; i < 2; i++) {
            refit_scene_bvh(bvh, node.primitive_ids[i]);
            node.bbox += bvh.nodes[node.primitive_ids[i]].bbox;
        }
    } else if (!bvh.instances.empty()) {
        for (auto i = 0; i < node.num_primitives; i++) {
            auto& instance = bvh.instances[node.primitive_ids[i]];
            auto  sbvh     = bvh.shape_bvhs[instance.shape_id];
            node.bbox += transform_bbox(instance.frame, sbvh.nodes[0].bbox);
        }
    } else {
        log_warning("empty bvh");
    }
}

// Recursively recomputes the node bounds for a shape bvh
void refit_shape_bvh(bvh_shape& bvh, const vector<vec3f>& positions) {
    bvh.positions = positions;
    refit_shape_bvh(bvh, 0);
}
void refit_shape_bvh(bvh_shape& bvh, const vector<vec3f>& positions,
    const vector<float>& radius) {
    bvh.positions = positions;
    bvh.radius    = radius;
    refit_shape_bvh(bvh, 0);
}
void refit_scene_bvh(bvh_scene& bvh, const vector<bvh_instance>& instances) {
    bvh.instances = instances;
    refit_scene_bvh(bvh, 0);
}

// Intersect ray with a bvh.
bool intersect_shape_bvh(const bvh_shape& bvh, const ray3f& ray_, bool find_any,
    float& distance, int& element_id, vec2f& element_uv) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh)
        return intersect_embree_bvh(
            bvh, ray_, find_any, distance, element_id, element_uv);
#endif

    // check empty
    if (bvh.nodes.empty()) return false;

    // node stack
    int  node_stack[128];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = vec3f{
        1 / ray.direction[0], 1 / ray.direction[1], 1 / ray.direction[2]};
    auto ray_dsign = vec3i{(ray_dinv[0] < 0) ? 1 : 0, (ray_dinv[1] < 0) ? 1 : 0,
        (ray_dinv[2] < 0) ? 1 : 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (node.is_internal) {
            // for internal nodes, attempts to proceed along the
            // split axis from smallest to largest nodes
            if (ray_dsign[node.split_axis]) {
                node_stack[node_cur++] = node.primitive_ids[0];
                node_stack[node_cur++] = node.primitive_ids[1];
            } else {
                node_stack[node_cur++] = node.primitive_ids[1];
                node_stack[node_cur++] = node.primitive_ids[0];
            }
        } else if (!bvh.triangles.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& t = bvh.triangles[node.primitive_ids[i]];
                if (intersect_triangle(ray, bvh.positions[t[0]],
                        bvh.positions[t[1]], bvh.positions[t[2]], distance,
                        element_uv)) {
                    hit        = true;
                    ray.tmax   = distance;
                    element_id = node.primitive_ids[i];
                }
            }
        } else if (!bvh.quads.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& t = bvh.quads[node.primitive_ids[i]];
                if (intersect_quad(ray, bvh.positions[t[0]],
                        bvh.positions[t[1]], bvh.positions[t[2]],
                        bvh.positions[t[3]], distance, element_uv)) {
                    hit        = true;
                    ray.tmax   = distance;
                    element_id = node.primitive_ids[i];
                }
            }
        } else if (!bvh.lines.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& l = bvh.lines[node.primitive_ids[i]];
                if (intersect_line(ray, bvh.positions[l[0]],
                        bvh.positions[l[1]], bvh.radius[l[0]], bvh.radius[l[1]],
                        distance, element_uv)) {
                    hit        = true;
                    ray.tmax   = distance;
                    element_id = node.primitive_ids[i];
                }
            }
        } else if (!bvh.points.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& p = bvh.points[node.primitive_ids[i]];
                if (intersect_point(ray, bvh.positions[p], bvh.radius[p],
                        distance, element_uv)) {
                    hit        = true;
                    ray.tmax   = distance;
                    element_id = node.primitive_ids[i];
                }
            }
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

// Intersect ray with a bvh.
bool intersect_scene_bvh(const bvh_scene& bvh, const ray3f& ray_, bool find_any,
    float& distance, int& instance_id, int& element_id, vec2f& element_uv) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh)
        return intersect_embree_bvh(
            bvh, ray_, find_any, distance, instance_id, element_id, element_uv);
#endif

    // check empty
    if (bvh.nodes.empty()) return false;

    // node stack
    int  node_stack[128];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = vec3f{
        1 / ray.direction[0], 1 / ray.direction[1], 1 / ray.direction[2]};
    auto ray_dsign = vec3i{(ray_dinv[0] < 0) ? 1 : 0, (ray_dinv[1] < 0) ? 1 : 0,
        (ray_dinv[2] < 0) ? 1 : 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (node.is_internal) {
            // for internal nodes, attempts to proceed along the
            // split axis from smallest to largest nodes
            if (ray_dsign[node.split_axis]) {
                node_stack[node_cur++] = node.primitive_ids[0];
                node_stack[node_cur++] = node.primitive_ids[1];
            } else {
                node_stack[node_cur++] = node.primitive_ids[1];
                node_stack[node_cur++] = node.primitive_ids[0];
            }
        } else if (!bvh.instances.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& instance = bvh.instances[node.primitive_ids[i]];
                if (instance.shape_id >= 0) {
                    if (intersect_shape_bvh(bvh.shape_bvhs[instance.shape_id],
                            transform_ray(instance.frame_inverse, ray),
                            find_any, distance, element_id, element_uv)) {
                        hit         = true;
                        ray.tmax    = distance;
                        instance_id = node.primitive_ids[i];
                    }
                } else if (instance.surface_id >= 0) {
                    if (intersect_shape_bvh(bvh.surface_bvhs[instance.surface_id],
                            transform_ray(instance.frame_inverse, ray),
                            find_any, distance, element_id, element_uv)) {
                        hit         = true;
                        ray.tmax    = distance;
                        instance_id = node.primitive_ids[i];
                    }
                }
            }
        } else {
            log_error("empty bvh");
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

// Finds the closest element with a bvh.
bool overlap_shape_bvh(const bvh_shape& bvh, const vec3f& pos,
    float max_distance, bool find_any, float& distance, int& instance_id,
    int& element_id, vec2f& element_uv) {
    // node stack
    int  node_stack[64];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // hit
    auto hit = false;

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (node.is_internal) {
            // internal node
            node_stack[node_cur++] = node.primitive_ids[0];
            node_stack[node_cur++] = node.primitive_ids[1];
        } else if (!bvh.triangles.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& t = bvh.triangles[node.primitive_ids[i]];
                if (overlap_triangle(pos, max_distance, bvh.positions[t[0]],
                        bvh.positions[t[1]], bvh.positions[t[2]],
                        bvh.radius[t[0]], bvh.radius[t[1]], bvh.radius[t[2]],
                        distance, element_uv)) {
                    hit          = true;
                    max_distance = distance;
                    element_id   = node.primitive_ids[i];
                }
            }
        } else if (!bvh.quads.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& q = bvh.quads[node.primitive_ids[i]];
                if (overlap_quad(pos, max_distance, bvh.positions[q[0]],
                        bvh.positions[q[1]], bvh.positions[q[2]],
                        bvh.positions[q[3]], bvh.radius[q[0]], bvh.radius[q[1]],
                        bvh.radius[q[2]], bvh.radius[q[3]], distance,
                        element_uv)) {
                    hit          = true;
                    max_distance = distance;
                    element_id   = node.primitive_ids[i];
                }
            }
        } else if (!bvh.lines.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& l = bvh.lines[node.primitive_ids[i]];
                if (overlap_line(pos, max_distance, bvh.positions[l[0]],
                        bvh.positions[l[1]], bvh.radius[l[0]], bvh.radius[l[1]],
                        distance, element_uv)) {
                    hit          = true;
                    max_distance = distance;
                    element_id   = node.primitive_ids[i];
                }
            }
        } else if (!bvh.points.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& p = bvh.points[node.primitive_ids[i]];
                if (overlap_point(pos, max_distance, bvh.positions[p],
                        bvh.radius[p], distance, element_uv)) {
                    hit          = true;
                    max_distance = distance;
                    element_id   = node.primitive_ids[i];
                }
            }
        } else {
            log_error("empty bvh");
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

// Finds the closest element with a bvh.
bool overlap_scene_bvh(const bvh_scene& bvh, const vec3f& pos,
    float max_distance, bool find_any, float& distance, int& instance_id,
    int& element_id, vec2f& element_uv) {
    // node stack
    int  node_stack[64];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // hit
    auto hit = false;

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh.nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (node.is_internal) {
            // internal node
            node_stack[node_cur++] = node.primitive_ids[0];
            node_stack[node_cur++] = node.primitive_ids[1];
        } else if (!bvh.instances.empty()) {
            for (auto i = 0; i < node.num_primitives; i++) {
                auto& instance = bvh.instances[node.primitive_ids[i]];
                if (overlap_shape_bvh(bvh.shape_bvhs[instance.shape_id],
                        transform_point(instance.frame_inverse, pos),
                        max_distance, find_any, distance, instance_id,
                        element_id, element_uv)) {
                    hit          = true;
                    max_distance = distance;
                    instance_id  = node.primitive_ids[i];
                }
            }
        } else {
            log_error("empty bvh");
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

#if 0
    // Finds the overlap between BVH leaf nodes.
    template <typename OverlapElem>
    void overlap_bvh_elems(const bvh_scene& bvh1, const bvh_scene& bvh2,
                           bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps,
                           const OverlapElem& overlap_elems) {
        // node stack
        vec2i node_stack[128];
        auto node_cur = 0;
        node_stack[node_cur++] = {0, 0};

        // walking stack
        while (node_cur) {
            // grab node
            auto node_idx = node_stack[--node_cur];
            const auto node1 = bvh1->nodes[node_idx[0]];
            const auto node2 = bvh2->nodes[node_idx[1]];

            // intersect bbox
            if (!overlap_bbox(node1.bbox, node2.bbox)) continue;

            // check for leaves
            if (node1.isleaf && node2.isleaf) {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count; i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (skip_duplicates && idx1 > idx2) continue;
                        if (skip_self && idx1 == idx2) continue;
                        if (overlap_elems(idx1, idx2))
                            overlaps.push_back({idx1, idx2});
                    }
                }
            } else {
                // descend
                if (node1.isleaf) {
                    for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                         idx2++) {
                        node_stack[node_cur++] = {node_idx[0], (int)idx2};
                    }
                } else if (node2.isleaf) {
                    for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                         idx1++) {
                        node_stack[node_cur++] = {(int)idx1, node_idx[1]};
                    }
                } else {
                    for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                         idx2++) {
                        for (auto idx1 = node1.start;
                             idx1 < node1.start + node1.count; idx1++) {
                            node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        }
                    }
                }
            }
        }
    }
#endif

}  // namespace yocto
