//
// Implementation for Yocto/BVH.
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
    const ray3f& ray, const vec3f& p, float r, vec2f& uv, float& dist) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = dot(w, ray.d) / dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp  = ray.o + ray.d * t;
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    // intersection occurred: set options and exit
    uv   = {0, 0};
    dist = t;
    return true;
}

// Intersect a ray with a line
bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, vec2f& uv, float& dist) {
    // setup intersection options
    auto u = ray.d;
    auto v = p1 - p0;
    auto w = ray.o - p0;

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
    auto pr  = ray.o + ray.d * t;
    auto pl  = p0 + (p1 - p0) * s;
    auto prl = pr - pl;

    // check with the line radius at the same point
    auto d2 = dot(prl, prl);
    auto r  = r0 * (1 - s) + r1 * s;
    if (d2 > r * r) return {};

    // intersection occurred: set options and exit
    uv   = {s, sqrt(d2) / r};
    dist = t;
    return true;
}

// Intersect a ray with a triangle
bool intersect_triangle(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, vec2f& uv, float& dist) {
    // compute triangle edges
    auto edge1 = p1 - p0;
    auto edge2 = p2 - p0;

    // compute determinant to solve a linear system
    auto pvec = cross(ray.d, edge2);
    auto det  = dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.o - p0;
    auto u    = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = cross(tvec, edge1);
    auto v    = dot(ray.d, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set options and exit
    uv   = {u, v};
    dist = t;
    return true;
}

// Intersect a ray with a quad.
bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, vec2f& uv, float& dist) {
    auto hit  = false;
    auto tray = ray;
    if (intersect_triangle(tray, p0, p1, p3, uv, dist)) {
        hit       = true;
        tray.tmax = dist;
    }
    if (intersect_triangle(tray, p2, p3, p1, uv, dist)) {
        hit       = true;
        uv        = 1 - uv;
        tray.tmax = dist;
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
    auto invd = 1.0f / ray.d;
    auto t0   = (bbox.min - ray.o) * invd;
    auto t1   = (bbox.max - ray.o) * invd;
    // flip based on range directions
    if (invd.x < 0.0f) swap(t0.x, t1.x);
    if (invd.y < 0.0f) swap(t0.y, t1.y);
    if (invd.z < 0.0f) swap(t0.z, t1.z);
    auto tmin = _safemax(t0.z, _safemax(t0.y, _safemax(t0.x, ray.tmin)));
    auto tmax = _safemin(t1.z, _safemin(t1.y, _safemin(t1.x, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox_) {
    auto bbox  = &bbox_.min;
    auto txmin = (bbox[ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto txmax = (bbox[1 - ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto tymin = (bbox[ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tymax = (bbox[1 - ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tzmin = (bbox[ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tzmax = (bbox[1 - ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tmin  = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax  = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(
    const ray3f& ray, const vec3f& ray_dinv, const bbox3f& bbox) {
    auto it_min = (bbox.min - ray.o) * ray_dinv;
    auto it_max = (bbox.max - ray.o) * ray_dinv;
    auto tmin   = vec3f{_safemin(it_min.x, it_max.x),
        _safemin(it_min.y, it_max.y), _safemin(it_min.z, it_max.z)};
    auto tmax   = vec3f{_safemax(it_min.x, it_max.x),
        _safemax(it_min.y, it_max.y), _safemax(it_min.z, it_max.z)};
    auto t0 = _safemax(tmin.x, _safemax(tmin.y, _safemax(tmin.z, ray.tmin)));
    auto t1 = _safemin(tmax.x, _safemin(tmax.y, _safemin(tmax.z, ray.tmax)));
    t1 *= 1.00000024f;  // for double: 1.0000000000000004
    return t0 <= t1;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// TODO: documentation
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p, float r,
    vec2f& uv, float& dist) {
    auto d2 = dot(pos - p, pos - p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    uv   = {0, 0};
    dist = sqrt(d2);
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
    const vec3f& p1, float r0, float r1, vec2f& uv, float& dist) {
    auto u = closestuv_line(pos, p0, p1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p  = p0 + (p1 - p0) * u;
    auto r  = r0 + (r1 - r0) * u;
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    uv   = {u, 0};
    dist = sqrt(d2);
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
    const vec3f& p1, const vec3f& p2, float r0, float r1, float r2, vec2f& uv,
    float& dist) {
    auto cuv = closestuv_triangle(pos, p0, p1, p2);
    auto p   = p0 * (1 - cuv.x - cuv.y) + p1 * cuv.x + p2 * cuv.y;
    auto r   = r0 * (1 - cuv.x - cuv.y) + r1 * cuv.x + r2 * cuv.y;
    auto dd  = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    uv   = cuv;
    dist = sqrt(dd);
    return true;
}

// TODO: documentation
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
    float r2, float r3, vec2f& uv, float& dist) {
    auto hit = false;
    if (overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r2, uv, dist)) {
        hit      = true;
        dist_max = dist;
    }
    if (!overlap_triangle(pos, dist_max, p2, p3, p1, r2, r3, r1, uv, dist)) {
        hit      = true;
        uv       = 1 - uv;
        dist_max = dist;
    }
    return hit;
}

// TODO: documentation
bool distance_check_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox) {
    // computing distance
    auto dd = 0.0f;

    // For each axis count any excess distance outside box extents
    if (pos.x < bbox.min.x) dd += (bbox.min.x - pos.x) * (bbox.min.x - pos.x);
    if (pos.x > bbox.max.x) dd += (pos.x - bbox.max.x) * (pos.x - bbox.max.x);
    if (pos.y < bbox.min.y) dd += (bbox.min.y - pos.y) * (bbox.min.y - pos.y);
    if (pos.y > bbox.max.y) dd += (pos.y - bbox.max.y) * (pos.y - bbox.max.y);
    if (pos.z < bbox.min.z) dd += (bbox.min.z - pos.z) * (bbox.min.z - pos.z);
    if (pos.z > bbox.max.z) dd += (pos.z - bbox.max.z) * (pos.z - bbox.max.z);

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
    if (bbox1.max.x < bbox2.min.x || bbox1.min.x > bbox2.max.x) return false;
    if (bbox1.max.y < bbox2.min.y || bbox1.min.y > bbox2.max.y) return false;
    if (bbox1.max.z < bbox2.min.z || bbox1.min.z > bbox2.max.z) return false;
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH
// -----------------------------------------------------------------------------
namespace yocto {

#if YOCTO_EMBREE
// Cleanup
bvh_shape::~bvh_shape() {
    if (embree_bvh) {
        // TODO: fix me
        // for (auto i = 0; i < max(1, (int)instances.size()); i++) {
        //     auto geom = rtcGetGeometry((RTCScene)embree_bvh, i);
        //     rtcDetachGeometry((RTCScene)embree_bvh, i);
        //     rtcReleaseGeometry(geom);
        // }
        // rtcReleaseScene((RTCScene)embree_bvh);
    }
}

// Cleanup
bvh_scene::~bvh_scene() {
    if (embree_bvh) {
        // TODO: fix me
        // for (auto i = 0; i < max(1, (int)instances.size()); i++) {
        //     auto geom = rtcGetGeometry((RTCScene)embree_bvh, i);
        //     rtcDetachGeometry((RTCScene)embree_bvh, i);
        //     rtcReleaseGeometry(geom);
        // }
        // rtcReleaseScene((RTCScene)embree_bvh);
    }
}

void embree_error(void* ctx, RTCError code, const char* str) {
    switch (code) {
        case RTC_ERROR_UNKNOWN:
            throw runtime_error("RTC_ERROR_UNKNOWN: "s + str);
            break;
        case RTC_ERROR_INVALID_ARGUMENT:
            throw runtime_error("RTC_ERROR_INVALID_ARGUMENT: "s + str);
            break;
        case RTC_ERROR_INVALID_OPERATION:
            throw runtime_error("RTC_ERROR_INVALID_OPERATION: "s + str);
            break;
        case RTC_ERROR_OUT_OF_MEMORY:
            throw runtime_error("RTC_ERROR_OUT_OF_MEMORY: "s + str);
            break;
        case RTC_ERROR_UNSUPPORTED_CPU:
            throw runtime_error("RTC_ERROR_UNSUPPORTED_CPU: "s + str);
            break;
        case RTC_ERROR_CANCELLED:
            throw runtime_error("RTC_ERROR_CANCELLED: "s + str);
            break;
        default: throw runtime_error("invalid error code"); break;
    }
}

// Embree memory
atomic<ssize_t> embree_memory = 0;
bool            embree_memory_monitor(void* userPtr, ssize_t bytes, bool post) {
    embree_memory += bytes;
    return true;
}

// Get Embree device
RTCDevice get_embree_device() {
    static RTCDevice device = nullptr;
    if (!device) {
        device = rtcNewDevice("");
        rtcSetDeviceErrorFunction(device, embree_error, nullptr);
        rtcSetDeviceMemoryMonitorFunction(
            device, embree_memory_monitor, nullptr);
    }
    return device;
}

// Initialize Embree BVH
void build_embree_points_bvh(bvh_shape& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if(options.embree_shared) {
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_COMPACT);
    }
    // rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.embree_bvh = embree_scene;
    throw runtime_error("embree does not support points");
    rtcCommitScene(embree_scene);
    bvh.embree_flattened = false;
}
void build_embree_lines_bvh(bvh_shape& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if(options.embree_shared) {
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_COMPACT);
    }
    // rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.embree_bvh  = embree_scene;
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto l : lines) {
        if (last_index == l.x) {
            elines.push_back((int)epositions.size() - 1);
            epositions.push_back({positions[l.y], radius[l.y]});
        } else {
            elines.push_back((int)epositions.size());
            epositions.push_back({positions[l.x], radius[l.x]});
            epositions.push_back({positions[l.y], radius[l.y]});
        }
        last_index = l.y;
    }
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_FLAT_LINEAR_CURVE);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, 4 * 4, epositions.size());
    auto embree_lines     = rtcSetNewGeometryBuffer(embree_geom,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, 4, elines.size());
    memcpy(embree_positions, epositions.data(), epositions.size() * 16);
    memcpy(embree_lines, elines.data(), elines.size() * 4);
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, 0);
    rtcCommitScene(embree_scene);
    bvh.embree_flattened = false;
}
void build_embree_triangles_bvh(bvh_shape& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const bvh_build_options& options) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if(options.embree_shared) {
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_COMPACT);
    }
    // rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.embree_bvh   = embree_scene;
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    if (options.embree_shared) {
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX, 0,
            RTC_FORMAT_FLOAT3, positions.data(), 0, 3 * 4, positions.size());
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_INDEX, 0,
            RTC_FORMAT_UINT3, triangles.data(), 0, 3 * 4, triangles.size());
    } else {
        auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
            positions.size());
        auto embree_triangles = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
            triangles.size());
        memcpy(embree_positions, positions.data(), positions.size() * 12);
        memcpy(embree_triangles, triangles.data(), triangles.size() * 12);
    }
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, 0);
    rtcCommitScene(embree_scene);
    bvh.embree_flattened = false;
}
void build_embree_quads_bvh(bvh_shape& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const bvh_build_options& options) {
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if(options.embree_shared) {
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_COMPACT);
    }
    // rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.embree_bvh   = embree_scene;
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    if (options.embree_shared) {
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX, 0,
            RTC_FORMAT_FLOAT3, positions.data(), 0, 3 * 4, positions.size());
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_INDEX, 0,
            RTC_FORMAT_UINT4, quads.data(), 0, 4 * 4, quads.size());
    } else {
        auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
            positions.size());
        auto embree_quads     = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4, quads.size());
        memcpy(embree_positions, positions.data(), positions.size() * 12);
        memcpy(embree_quads, quads.data(), quads.size() * 16);
    }
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, 0);
    rtcCommitScene(embree_scene);
    bvh.embree_flattened = false;
}

// Initialize Embree BVH
void build_embree_points_bvh(RTCScene& embree_scene, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool try_sharing, int id = 0) {
    throw runtime_error("embree does not support points");
}
void build_embree_lines_bvh(RTCScene& embree_scene, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool try_sharing, int id = 0) {
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto l : lines) {
        if (last_index == l.x) {
            elines.push_back((int)epositions.size() - 1);
            epositions.push_back({positions[l.y], radius[l.y]});
        } else {
            elines.push_back((int)epositions.size());
            epositions.push_back({positions[l.x], radius[l.x]});
            epositions.push_back({positions[l.y], radius[l.y]});
        }
        last_index = l.y;
    }
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_FLAT_LINEAR_CURVE);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, 4 * 4, epositions.size());
    auto embree_lines     = rtcSetNewGeometryBuffer(embree_geom,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, 4, elines.size());
    memcpy(embree_positions, epositions.data(), epositions.size() * 16);
    memcpy(embree_lines, elines.data(), elines.size() * 4);
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, id);
}
void build_embree_triangles_bvh(RTCScene& embree_scene,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    bool try_sharing, int id = 0) {
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    if (try_sharing) {
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX, 0,
            RTC_FORMAT_FLOAT3, positions.data(), 0, 3 * 4, positions.size());
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_INDEX, 0,
            RTC_FORMAT_UINT3, triangles.data(), 0, 3 * 4, triangles.size());
    } else {
        auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
            positions.size());
        auto embree_triangles = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
            triangles.size());
        memcpy(embree_positions, positions.data(), positions.size() * 12);
        memcpy(embree_triangles, triangles.data(), triangles.size() * 12);
    }
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, 0);
}
void build_embree_quads_bvh(RTCScene& embree_scene, const vector<vec4i>& quads,
    const vector<vec3f>& positions, bool try_sharing, int id = 0) {
    auto embree_geom = rtcNewGeometry(
        get_embree_device(), RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(embree_geom, 1);
    if (try_sharing) {
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_VERTEX, 0,
            RTC_FORMAT_FLOAT3, positions.data(), 0, 3 * 4, positions.size());
        rtcSetSharedGeometryBuffer(embree_geom, RTC_BUFFER_TYPE_INDEX, 0,
            RTC_FORMAT_UINT4, quads.data(), 0, 4 * 4, quads.size());
    } else {
        auto embree_positions = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
            positions.size());
        auto embree_quads     = rtcSetNewGeometryBuffer(embree_geom,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4, quads.size());
        memcpy(embree_positions, positions.data(), positions.size() * 12);
        memcpy(embree_quads, quads.data(), quads.size() * 16);
    }
    rtcCommitGeometry(embree_geom);
    rtcAttachGeometryByID(embree_scene, embree_geom, id);
}

// Build a BVH using Embree.
void build_embree_instances_bvh(bvh_scene& bvh, int num_instances,
    const void* context, bvh_instance (*get_instance)(const void*, int),
    const bvh_build_options& options) {
    // scene bvh
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    if(options.embree_shared) {
        rtcSetSceneFlags(embree_scene, RTC_SCENE_FLAG_COMPACT);
    }
    // rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.embree_bvh = embree_scene;
    if (!num_instances) {
        rtcCommitScene(embree_scene);
        return;
    }
    for (auto instance_id = 0; instance_id < num_instances; instance_id++) {
        auto instance = get_instance(context, instance_id);
        if (instance.shape_id < 0) throw runtime_error("empty instance");
        auto& shape_bvh = bvh.shapes[instance.shape_id];
        if (!shape_bvh.embree_bvh) throw runtime_error("bvh not built");
        auto embree_geom = rtcNewGeometry(
            embree_device, RTC_GEOMETRY_TYPE_INSTANCE);
        rtcSetGeometryInstancedScene(
            embree_geom, (RTCScene)shape_bvh.embree_bvh);
        rtcSetGeometryTransform(
            embree_geom, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
        rtcCommitGeometry(embree_geom);
        rtcAttachGeometryByID(embree_scene, embree_geom, instance_id);
    }
    rtcCommitScene(embree_scene);
    bvh.embree_flattened = false;
}
void build_embree_flattened_bvh(
    bvh_scene_data& bvh, const bvh_build_options& options) {
    // scene bvh
    auto embree_device = get_embree_device();
    auto embree_scene  = rtcNewScene(embree_device);
    rtcSetSceneBuildQuality(embree_scene, RTC_BUILD_QUALITY_HIGH);
    bvh.bvh_.embree_bvh = embree_scene;
    if (bvh.instances.empty()) {
        rtcCommitScene(embree_scene);
        return;
    }
    for (auto instance_id = 0; instance_id < bvh.instances.size();
         instance_id++) {
        auto& instance  = bvh.instances[instance_id];
        auto& shape_bvh = bvh.shapes[instance.shape_id];
        if (shape_bvh.positions.empty()) continue;
        auto transformed_positions = shape_bvh.positions;
        if (instance.frame != identity_frame3f) {
            for (auto& p : transformed_positions)
                p = transform_point(instance.frame, p);
        }
        if (!shape_bvh.points.empty()) {
            build_embree_points_bvh(embree_scene, shape_bvh.points,
                transformed_positions, shape_bvh.radius, false, instance_id);
        } else if (!shape_bvh.lines.empty()) {
            build_embree_lines_bvh(embree_scene, shape_bvh.lines,
                transformed_positions, shape_bvh.radius, false, instance_id);
        } else if (!shape_bvh.triangles.empty()) {
            build_embree_triangles_bvh(embree_scene, shape_bvh.triangles,
                transformed_positions, false, instance_id);
        } else if (!shape_bvh.quads.empty()) {
            build_embree_quads_bvh(embree_scene, shape_bvh.quads,
                transformed_positions, false, instance_id);
        } else {
            throw runtime_error("empty bvh");
        }
    }
    rtcCommitScene(embree_scene);
    bvh.bvh_.embree_flattened = true;
}
// Refit a BVH using Embree. Calls `refit_scene_bvh()` if Embree is not
// available.
void refit_embree_bvh(bvh_shape& bvh) {
    throw runtime_error("not yet implemented");
}
bool intersect_embree_bvh(const bvh_shape& bvh, const ray3f& ray,
    bvh_intersection& intersection, bool find_any) {
    RTCRayHit embree_ray;
    embree_ray.ray.org_x     = ray.o.x;
    embree_ray.ray.org_y     = ray.o.y;
    embree_ray.ray.org_z     = ray.o.z;
    embree_ray.ray.dir_x     = ray.d.x;
    embree_ray.ray.dir_y     = ray.d.y;
    embree_ray.ray.dir_z     = ray.d.z;
    embree_ray.ray.tnear     = ray.tmin;
    embree_ray.ray.tfar      = ray.tmax;
    embree_ray.ray.flags     = 0;
    embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
    embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    RTCIntersectContext embree_ctx;
    rtcInitIntersectContext(&embree_ctx);
    rtcIntersect1((RTCScene)bvh.embree_bvh, &embree_ctx, &embree_ray);
    if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
    intersection.instance_id = bvh.embree_flattened
                                   ? (int)embree_ray.hit.geomID
                                   : (int)embree_ray.hit.instID[0];
    intersection.element_id = (int)embree_ray.hit.primID;
    intersection.element_uv = {embree_ray.hit.u, embree_ray.hit.v};
    intersection.distance   = embree_ray.ray.tfar;
    return true;
}
bool intersect_embree_bvh(const bvh_scene& bvh, const ray3f& ray,
    bvh_intersection& intersection, bool find_any) {
    RTCRayHit embree_ray;
    embree_ray.ray.org_x     = ray.o.x;
    embree_ray.ray.org_y     = ray.o.y;
    embree_ray.ray.org_z     = ray.o.z;
    embree_ray.ray.dir_x     = ray.d.x;
    embree_ray.ray.dir_y     = ray.d.y;
    embree_ray.ray.dir_z     = ray.d.z;
    embree_ray.ray.tnear     = ray.tmin;
    embree_ray.ray.tfar      = ray.tmax;
    embree_ray.ray.flags     = 0;
    embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
    embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
    RTCIntersectContext embree_ctx;
    rtcInitIntersectContext(&embree_ctx);
    rtcIntersect1((RTCScene)bvh.embree_bvh, &embree_ctx, &embree_ray);
    if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
    intersection.instance_id = bvh.embree_flattened
                                   ? (int)embree_ray.hit.geomID
                                   : (int)embree_ray.hit.instID[0];
    intersection.element_id = (int)embree_ray.hit.primID;
    intersection.element_uv = {embree_ray.hit.u, embree_ray.hit.v};
    intersection.distance   = embree_ray.ray.tfar;
    return true;
}
#endif

// BVH primitive with its bbox, its center and the index to the primitive
struct bvh_prim {
    bbox3f bbox   = invalid_bbox3f;
    vec3f  center = zero3f;
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
    if (csize == zero3f) return {mid, split_axis};

    // consider N bins, compute their cost and keep the minimum
    const int nbins     = 16;
    auto      middle    = 0.0f;
    auto      min_cost  = float_max;
    auto      bbox_area = [](auto& b) {
        auto size = b.max - b.min;
        return 1e-12f + 2 * size.x * size.y + 2 * size.x * size.z +
               2 * size.y * size.z;
    };
    for (auto axis = 0; axis < 3; axis++) {
        for (auto b = 1; b < nbins; b++) {
            auto split     = (&cbbox.min.x)[axis] + b * csize[axis] / nbins;
            auto left_bbox = invalid_bbox3f, right_bbox = invalid_bbox3f;
            auto left_nprims = 0, right_nprims = 0;
            for (auto i = start; i < end; i++) {
                if ((&prims[i].center.x)[axis] < split) {
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
        throw runtime_error("bad bvh split");
        split_axis = 0;
        mid        = (start + end) / 2;
    }

    return {mid, split_axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and
// axis.
pair<int, int> split_bvh_node_balanced(
    vector<bvh_prim>& prims, int start, int end) {
    // initialize split axis and position
    auto split_axis = 0;
    auto mid        = (start + end) / 2;

    // compute primintive bounds and size
    auto cbbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) cbbox += prims[i].center;
    auto csize = cbbox.max - cbbox.min;
    if (csize == zero3f) return {mid, split_axis};

    // split along largest
    auto largest_axis = 0;
    if (csize.x >= csize.y && csize.x >= csize.z) largest_axis = 0;
    if (csize.y >= csize.x && csize.y >= csize.z) largest_axis = 1;
    if (csize.z >= csize.x && csize.z >= csize.y) largest_axis = 2;

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
        throw runtime_error("bad bvh split");
        split_axis = 0;
        mid        = (start + end) / 2;
    }

    return {mid, split_axis};
}

// Splits a BVH node using the middle heutirtic. Returns split position and
// axis.
pair<int, int> split_bvh_node_middle(
    vector<bvh_prim>& prims, int start, int end) {
    // initialize split axis and position
    auto split_axis = 0;
    auto mid        = (start + end) / 2;

    // compute primintive bounds and size
    auto cbbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) cbbox += prims[i].center;
    auto csize = cbbox.max - cbbox.min;
    if (csize == zero3f) return {mid, split_axis};

    // split along largest
    auto largest_axis = 0;
    if (csize.x >= csize.y && csize.x >= csize.z) largest_axis = 0;
    if (csize.y >= csize.x && csize.y >= csize.z) largest_axis = 1;
    if (csize.z >= csize.x && csize.z >= csize.y) largest_axis = 2;

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
        throw runtime_error("bad bvh split");
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
    auto nodeid = next.x, start = next.y, end = next.z;

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

    // split into two children
    if (end - start > bvh_max_prims) {
        // get split
        auto [mid, split_axis] =
            (high_quality) ? split_bvh_node_sah(prims, start, end)
                           : split_bvh_node_balanced(prims, start, end);

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
    const bvh_build_options& options) {
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
        auto nodeid = next.x, start = next.y, end = next.z;

        // grab node
        auto& node = nodes[nodeid];

        // compute bounds
        node.bbox = invalid_bbox3f;
        for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

        // split into two children
        if (end - start > bvh_max_prims) {
            // get split
            auto [mid, split_axis] = (options.high_quality)
                                         ? split_bvh_node_sah(prims, start, end)
                                         : split_bvh_node_balanced(
                                               prims, start, end);

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
    const bvh_build_options& options) {
    // prepare to build nodes
    nodes.clear();
    nodes.reserve(prims.size() * 2);

    // queue up first node
    auto queue = deque<vec3i>{{0, 0, (int)prims.size()}};
    nodes.emplace_back();

    // synchronization
    atomic<int>          num_processed_prims(0);
    mutex                queue_mutex;
    vector<future<void>> futures;
    auto                 nthreads = thread::hardware_concurrency();

    // create nodes until the queue is empty
    for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
        futures.emplace_back(async([&nodes, &prims, &options,
                                       &num_processed_prims, &queue_mutex,
                                       &queue] {
            while (true) {
                // exit if needed
                if (num_processed_prims >= prims.size()) return;
                if (options.cancel_flag && *options.cancel_flag) return;

                // grab node to work on
                auto next = zero3i;
                {
                    lock_guard<mutex> lock{queue_mutex};
                    if (!queue.empty()) {
                        next = queue.front();
                        queue.pop_front();
                    }
                }

                // wait a bit if needed
                if (next == zero3i) {
                    std::this_thread::sleep_for(10us);
                    continue;
                }

                // grab node
                auto  nodeid = next.x, start = next.y, end = next.z;
                auto& node = nodes[nodeid];

                // compute bounds
                node.bbox = invalid_bbox3f;
                for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

                // split into two children
                if (end - start > bvh_max_prims) {
                    // get split
                    auto [mid, split_axis] =
                        (options.high_quality)
                            ? split_bvh_node_sah(prims, start, end)
                            : split_bvh_node_balanced(prims, start, end);

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
        }));
    }
    for (auto& f : futures) f.get();

    // cleanup
    nodes.shrink_to_fit();
}

// Build a BVH from a set of primitives.
template <typename ElemBounds>
void build_bvh_nodes(vector<bvh_node>& nodes, size_t num_elements,
    const ElemBounds& element_bounds, const bvh_build_options& options) {
    // get the number of primitives and the primitive type
    auto prims = vector<bvh_prim>(num_elements);
    for (auto element_id = 0; element_id < num_elements; element_id++) {
        prims[element_id].bbox   = element_bounds(element_id);
        prims[element_id].center = bbox_center(prims[element_id].bbox);
        prims[element_id].primid = element_id;
    }

    // build nodes
    if (options.run_serially) {
        build_bvh_nodes_serial(nodes, prims, options);
    } else {
        build_bvh_nodes_parallel(nodes, prims, options);
    }
}

// Build the bvh acceleration structure.
void build_points_bvh(bvh_shape& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (options.use_embree) {
        return build_embree_points_bvh(bvh, points, positions, radius, options);
    }
#endif

    // get the number of primitives and the primitive type
    return build_bvh_nodes(
        bvh.nodes, points.size(),
        [&points, &positions, &radius](int idx) {
            auto& p = points[idx];
            return point_bounds(positions[p], radius[p]);
        },
        options);
}
void build_lines_bvh(bvh_shape& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (options.use_embree) {
        return build_embree_lines_bvh(bvh, lines, positions, radius, options);
    }
#endif

    return build_bvh_nodes(
        bvh.nodes, lines.size(),
        [&lines, &positions, &radius](int idx) {
            auto& l = lines[idx];
            return line_bounds(
                positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
        },
        options);
}
void build_triangles_bvh(bvh_shape& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const bvh_build_options& options) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (options.use_embree) {
        return build_embree_triangles_bvh(bvh, triangles, positions, options);
    }
#endif

    return build_bvh_nodes(
        bvh.nodes, triangles.size(),
        [&triangles, &positions](int idx) {
            auto& t = triangles[idx];
            return triangle_bounds(
                positions[t.x], positions[t.y], positions[t.z]);
        },
        options);
}
void build_quads_bvh(bvh_shape& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const bvh_build_options& options) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (options.use_embree) {
        return build_embree_quads_bvh(bvh, quads, positions, options);
    }
#endif

    return build_bvh_nodes(
        bvh.nodes, quads.size(),
        [&quads, &positions](int idx) {
            auto& q = quads[idx];
            return quad_bounds(
                positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        },
        options);
}

void build_instances_bvh(bvh_scene& bvh, int num_instances,
    const void* context, bvh_instance (*get_instance)(const void*, int),
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (options.use_embree) {
        return build_embree_instances_bvh(
            bvh, num_instances, context, get_instance, options);
    }
#endif

    if (num_instances) {
        // get the number of primitives and the primitive type
        return build_bvh_nodes(
            bvh.nodes, num_instances,
            [&bvh, context, get_instance](int idx) {
                auto  instance = get_instance(context, idx);
                auto& sbvh     = bvh.shapes[instance.shape_id];
                return sbvh.nodes.empty()
                           ? invalid_bbox3f
                           : transform_bbox(instance.frame, sbvh.nodes[0].bbox);
            },
            options);
    }
}

// Recursively recomputes the node bounds for a shape bvh
template <typename ElemBound>
void refit_bvh_nodes(
    vector<bvh_node>& nodes, int nodeid, const ElemBound& element_bounds) {
    // refit
    auto& node = nodes[nodeid];
    node.bbox  = invalid_bbox3f;
    if (node.is_internal) {
        for (auto i = 0; i < 2; i++) {
            refit_bvh_nodes(nodes, node.primitive_ids[i], element_bounds);
            node.bbox += nodes[node.primitive_ids[i]].bbox;
        }
    } else {
        for (auto elememt_id = 0; elememt_id < node.num_primitives;
             elememt_id++) {
            node.bbox += element_bounds(elememt_id);
        }
    }
}

void refit_points_bvh(bvh_shape& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) throw runtime_error("Embree reftting disabled");
#endif

    return refit_bvh_nodes(
        bvh.nodes, 0, [&points, &positions, &radius](int idx) {
            auto& p = points[idx];
            return point_bounds(positions[p], radius[p]);
        });
}
void refit_lines_bvh(bvh_shape& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) throw runtime_error("Embree reftting disabled");
#endif

    return refit_bvh_nodes(
        bvh.nodes, 0, [&lines, &positions, &radius](int idx) {
            auto& l = lines[idx];
            return line_bounds(
                positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
        });
}
void refit_triangles_bvh(bvh_shape& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) throw runtime_error("Embree reftting disabled");
#endif

    return refit_bvh_nodes(bvh.nodes, 0, [&triangles, &positions](int idx) {
        auto& t = triangles[idx];
        return triangle_bounds(positions[t.x], positions[t.y], positions[t.z]);
    });
}
void refit_quads_bvh(bvh_shape& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) throw runtime_error("Embree reftting disabled");
#endif

    return refit_bvh_nodes(bvh.nodes, 0, [&quads, &positions](int idx) {
        auto& q = quads[idx];
        return quad_bounds(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    });
}
void refit_instances_bvh(bvh_scene& bvh, int num_instances,
    const void* context, bvh_instance (*get_instance)(const void*, int),
    const bvh_build_options& options) {
#if YOCTO_EMBREE
    if (bvh.embree_bvh) throw runtime_error("Embree reftting disabled");
#endif

    if (num_instances) {
        // get the number of primitives and the primitive type
        return refit_bvh_nodes(
            bvh.nodes, 0, [get_instance, context, &bvh](int idx) {
                auto  instance = get_instance(context, idx);
                auto& sbvh     = bvh.shapes[instance.shape_id];
                return sbvh.nodes.empty()
                           ? invalid_bbox3f
                           : transform_bbox(instance.frame, sbvh.nodes[0].bbox);
            });
    } else {
        throw runtime_error("empty instances");
    }
}

// Intersect ray with a bvh.
template <bool IsInstanced, typename PrimitiveIntersect>
bool intersect_bvh_nodes(const vector<bvh_node>& nodes, const ray3f& ray_,
    bvh_intersection& intersection, bool find_any,
    const PrimitiveIntersect& intersect_primitive) {
    // check empty
    if (nodes.empty()) return false;

    // node stack
    int  node_stack[128];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
    auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
        (ray_dinv.z < 0) ? 1 : 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = nodes[node_stack[--node_cur]];

        // intersect bbox
        // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
        if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

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
        } else {
            for (auto i = 0; i < node.num_primitives; i++) {
                if constexpr (!IsInstanced) {
                    if (intersect_primitive(ray, node.primitive_ids[i],
                            intersection.element_uv, intersection.distance)) {
                        hit                     = true;
                        intersection.element_id = node.primitive_ids[i];
                        ray.tmax                = intersection.distance;
                    }
                } else {
                    if (intersect_primitive(ray, node.primitive_ids[i],
                            intersection, find_any)) {
                        hit                      = true;
                        intersection.instance_id = node.primitive_ids[i];
                        ray.tmax                 = intersection.distance;
                    }
                }
            }
        }

        // check for early exit
        if (find_any && hit) return hit;
    }

    return hit;
}

bool intersect_points_bvh(const bvh_shape& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, bvh_intersection& intersection, bool find_any) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh) {
        return intersect_embree_bvh(bvh, ray, intersection, find_any);
    }
#endif

    if (!points.empty()) {
        return intersect_bvh_nodes<false>(bvh.nodes, ray, intersection,
            find_any,
            [&points, &positions, &radius](
                const ray3f& ray, int idx, vec2f& uv, float& dist) {
                auto& p = points[idx];
                return intersect_point(ray, positions[p], radius[p], uv, dist);
            });
    } else {
        return false;
    }
}
bool intersect_lines_bvh(const bvh_shape& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, bvh_intersection& intersection, bool find_any) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh) {
        return intersect_embree_bvh(bvh, ray, intersection, find_any);
    }
#endif

    if (!lines.empty()) {
        return intersect_bvh_nodes<false>(bvh.nodes, ray, intersection,
            find_any,
            [&lines, &positions, &radius](
                const ray3f& ray, int idx, vec2f& uv, float& dist) {
                auto& l = lines[idx];
                return intersect_line(ray, positions[l.x], positions[l.y],
                    radius[l.x], radius[l.y], uv, dist);
            });
    } else {
        return false;
    }
}
bool intersect_triangles_bvh(const bvh_shape& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, bvh_intersection& intersection, bool find_any) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh) {
        return intersect_embree_bvh(bvh, ray, intersection, find_any);
    }
#endif

    if (!triangles.empty()) {
        return intersect_bvh_nodes<false>(bvh.nodes, ray, intersection,
            find_any,
            [&triangles, &positions](
                const ray3f& ray, int idx, vec2f& uv, float& dist) {
                auto& t = triangles[idx];
                return intersect_triangle(ray, positions[t.x], positions[t.y],
                    positions[t.z], uv, dist);
            });
    } else {
        return false;
    }
}
bool intersect_quads_bvh(const bvh_shape& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const ray3f& ray,
    bvh_intersection& intersection, bool find_any) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh) {
        return intersect_embree_bvh(bvh, ray, intersection, find_any);
    }
#endif

    if (!quads.empty()) {
        return intersect_bvh_nodes<false>(bvh.nodes, ray, intersection,
            find_any,
            [&quads, &positions](
                const ray3f& ray, int idx, vec2f& uv, float& dist) {
                auto& q = quads[idx];
                return intersect_quad(ray, positions[q.x], positions[q.y],
                    positions[q.z], positions[q.w], uv, dist);
            });
    } else {
        return false;
    }
}

bool intersect_instances_bvh(const bvh_scene& bvh, int num_instances,
    const void* context, bvh_instance (*get_instance)(const void*, int),
    bool (*intersect_shape)(const bvh_shape&, const void*, int, const ray3f&,
        bvh_intersection&, bool),
    const ray3f& ray, bvh_intersection& intersection, bool find_any,
    bool non_rigid_frames) {
#if YOCTO_EMBREE
    // call Embree if needed
    if (bvh.embree_bvh) {
        return intersect_embree_bvh(bvh, ray, intersection, find_any);
    }
#endif

    if (num_instances) {
        return intersect_bvh_nodes<true>(bvh.nodes, ray, intersection, find_any,
            [&bvh, get_instance, context, intersect_shape, non_rigid_frames](
                const ray3f& ray, int idx, bvh_intersection& intersection,
                bool find_any) {
                auto instance = get_instance(context, idx);
                auto inv_ray  = non_rigid_frames
                                   ? transform_ray(
                                         inverse((affine3f)instance.frame), ray)
                                   : transform_ray_inverse(instance.frame, ray);
                return intersect_shape(bvh.shapes[instance.shape_id], context,
                    instance.shape_id, inv_ray, intersection, find_any);
            });
    } else {
        return false;
    }
}

// Intersect ray with a bvh.
template <bool IsInstanced, typename PrimitiveOverlap>
bool overlap_bvh_nodes(const vector<bvh_node>& nodes, const vec3f& pos,
    float max_distance, bvh_intersection& intersection, bool find_any,
    const PrimitiveOverlap& overlap_primitive) {
    // check if empty
    if (nodes.empty()) return false;

    // node stack
    int  node_stack[64];
    auto node_cur          = 0;
    node_stack[node_cur++] = 0;

    // hit
    auto hit = false;

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (node.is_internal) {
            // internal node
            node_stack[node_cur++] = node.primitive_ids[0];
            node_stack[node_cur++] = node.primitive_ids[1];
        } else {
            for (auto i = 0; i < node.num_primitives; i++) {
                if constexpr (!IsInstanced) {
                    if (overlap_primitive(pos, max_distance,
                            node.primitive_ids[i], intersection.element_uv,
                            intersection.distance)) {
                        hit                     = true;
                        intersection.element_id = node.primitive_ids[i];
                        max_distance            = intersection.distance;
                    }
                } else {
                    if (overlap_primitive(pos, max_distance,
                            node.primitive_ids[i], intersection, find_any)) {
                        hit                      = true;
                        intersection.instance_id = node.primitive_ids[i];
                        max_distance             = intersection.distance;
                    }
                }
            }
        }

        // check for early exit
        if (find_any && hit) return hit;
    }

    return hit;
}

// Finds the closest element with a bvh.
bool overlap_points_bvh(const bvh_shape& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, bvh_intersection& intersection,
    bool find_any) {
    return overlap_bvh_nodes<false>(bvh.nodes, pos, max_distance, intersection,
        find_any,
        [&points, &positions, &radius](const vec3f& pos, float max_distance,
            int idx, vec2f& element_uv, float& distance) {
            auto& p = points[idx];
            return overlap_point(pos, max_distance, positions[p], radius[p],
                element_uv, distance);
        });
}
bool overlap_lines_bvh(const bvh_shape& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, bvh_intersection& intersection,
    bool find_any) {
    return overlap_bvh_nodes<false>(bvh.nodes, pos, max_distance, intersection,
        find_any,
        [&lines, &positions, &radius](const vec3f& pos, float max_distance,
            int idx, vec2f& element_uv, float& distance) {
            auto& l = lines[idx];
            return overlap_line(pos, max_distance, positions[l.x],
                positions[l.y], radius[l.x], radius[l.y], element_uv, distance);
        });
}
bool overlap_triangles_bvh(const bvh_shape& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vec3f& pos, float max_distance,
    bvh_intersection& intersection, bool find_any) {
    return overlap_bvh_nodes<false>(bvh.nodes, pos, max_distance, intersection,
        find_any,
        [&triangles, &positions](const vec3f& pos, float max_distance, int idx,
            vec2f& element_uv, float& distance) {
            auto& t = triangles[idx];
            return overlap_triangle(pos, max_distance, positions[t.x],
                positions[t.y], positions[t.z], 0, 0, 0, element_uv, distance);
            // return overlap_triangle(pos, max_distance, positions[t.x],
            //     positions[t.y], positions[t.z], radius[t.x],
            //     radius[t.y], radius[t.z], element_uv, distance);
        });
}
bool overlap_quads_bvh(const bvh_shape& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vec3f& pos, float max_distance,
    bvh_intersection& intersection, bool find_any) {
    return overlap_bvh_nodes<false>(bvh.nodes, pos, max_distance, intersection,
        find_any,
        [&quads, &positions](const vec3f& pos, float max_distance, int idx,
            vec2f& element_uv, float& distance) {
            auto& q = quads[idx];
            return overlap_quad(pos, max_distance, positions[q.x],
                positions[q.y], positions[q.z], positions[q.w], 0, 0, 0, 0,
                element_uv, distance);
            // return overlap_quad(pos, max_distance, positions[q.x],
            //     positions[q.y], positions[q.z], positions[q.w],
            //     radius[q.x], radius[q.y], radius[q.z],
            //     radius[q.w], element_uv, distance);
        });
}

// Finds the closest element with a bvh.
bool overlap_instances_bvh(const bvh_scene& bvh, int num_instances,
    const void* context, bvh_instance (*get_instance)(const void*, int),
    bool (*overlap_shape)(const bvh_shape&, const void*, int, const vec3f&,
        float, bvh_intersection&, bool),
    const vec3f& pos, float max_distance, bvh_intersection& intersection,
    bool find_any, bool non_rigid_frames) {
    if (num_instances) {
        return overlap_bvh_nodes<true>(bvh.nodes, pos, max_distance,
            intersection, find_any,
            [&bvh, get_instance, context, &overlap_shape, non_rigid_frames](
                const vec3f& pos, float max_distance, int idx,
                bvh_intersection& intersection, bool find_any) {
                auto instance = get_instance(context, idx);
                auto inv_pos  = non_rigid_frames
                                   ? transform_point(
                                         inverse((affine3f)instance.frame), pos)
                                   : transform_point_inverse(
                                         instance.frame, pos);
                return overlap_shape(bvh.shapes[instance.shape_id], context,
                    instance.shape_id, inv_pos, max_distance, intersection,
                    find_any);
            });

    } else {
        return false;
    }
}

#if 0
    // Finds the overlap between BVH leaf nodes.
    template <typename OverlapElem>
    void overlap_bvh_elems(const bvh_scene_data& bvh1, const bvh_scene_data& bvh2,
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
            const auto node1 = bvh1->nodes[node_idx.x];
            const auto node2 = bvh2->nodes[node_idx.y];

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
                        node_stack[node_cur++] = {node_idx.x, (int)idx2};
                    }
                } else if (node2.isleaf) {
                    for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                         idx1++) {
                        node_stack[node_cur++] = {(int)idx1, node_idx.y};
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF BVH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Print bvh statistics.
string print_shape_bvh_stats(const bvh_shape& bvh) {
    // TODO
    auto str = ""s;
    return str;
}
string print_scene_bvh_stats(const bvh_scene& bvh) {
#if 0
    auto num_shapes    = (size_t)0;
    auto num_instances = (size_t)0;

    auto elem_points    = (size_t)0;
    auto elem_lines     = (size_t)0;
    auto elem_triangles = (size_t)0;
    auto elem_quads     = (size_t)0;

    auto vert_pos    = (size_t)0;
    auto vert_radius = (size_t)0;

    auto shape_nodes = (size_t)0;
    auto scene_nodes = (size_t)0;

    auto stored_elem_points    = (size_t)0;
    auto stored_elem_lines     = (size_t)0;
    auto stored_elem_triangles = (size_t)0;
    auto stored_elem_quads     = (size_t)0;

    auto stored_vert_pos    = (size_t)0;
    auto stored_vert_radius = (size_t)0;

    auto memory_elems = (size_t)0;
    auto memory_verts = (size_t)0;

    auto memory_ists = (size_t)0;

    auto memory_shape_nodes = (size_t)0;
    auto memory_scene_nodes = (size_t)0;

    for (auto& sbvh : bvh.shapes) {
        shape_nodes += sbvh.nodes.size();
    }

    num_shapes    = bvh.shapes.size();
    num_instances = bvh.instances.size();
    scene_nodes   = bvh.bvh_.nodes.size();

    memory_elems = stored_elem_points * sizeof(int) +
                   stored_elem_lines * sizeof(vec2i) +
                   stored_elem_triangles * sizeof(vec3i) +
                   stored_elem_quads * sizeof(vec4i);
    memory_verts = stored_vert_pos * sizeof(vec3f) +
                   stored_vert_radius * sizeof(float);

    memory_ists = num_instances * sizeof(bvh_instance);

    memory_shape_nodes = shape_nodes * sizeof(bvh_node);
    memory_scene_nodes = scene_nodes * sizeof(bvh_node);

    auto str = ""s;

    str += "num_shapes: " + std::to_string(num_shapes) + "\n";
    str += "num_instances: " + std::to_string(num_instances) + "\n";

    str += "elem_points: " + std::to_string(elem_points) + "\n";
    str += "elem_lines: " + std::to_string(elem_lines) + "\n";
    str += "elem_triangles: " + std::to_string(elem_triangles) + "\n";
    str += "elem_quads: " + std::to_string(elem_quads) + "\n";
    str += "vert_pos: " + std::to_string(vert_pos) + "\n";
    str += "vert_radius: " + std::to_string(vert_radius) + "\n";

    str += "shape_nodes: " + std::to_string(shape_nodes) + "\n";
    str += "scene_nodes: " + std::to_string(scene_nodes) + "\n";

    str += "memory_elems: " + std::to_string(memory_elems) + "\n";
    str += "memory_verts: " + std::to_string(memory_verts) + "\n";
    str += "memory_ists: " + std::to_string(memory_ists) + "\n";
    str += "memory_shape_nodes: " + std::to_string(memory_shape_nodes) + "\n";
    str += "memory_scene_nodes: " + std::to_string(memory_scene_nodes) + "\n";

#if YOCTO_EMBREE
    str += "memory_embree: " + std::to_string(embree_memory) + "\n";
#endif
#endif
    // TODO
    auto str = ""s;
    return str;
}

}  // namespace yocto
