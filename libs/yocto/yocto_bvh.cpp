//
// Implementation for Yocto/Bvh
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstring>
#include <deque>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "yocto_geometry.h"
#include "yocto_parallel.h"

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::atomic;
using std::deque;
using std::pair;
using std::string;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EMBREE BVH
// -----------------------------------------------------------------------------
namespace yocto {

#ifdef YOCTO_EMBREE

// Get Embree device
atomic<ssize_t>  bvh_embree_memory = 0;
static RTCDevice bvh_embree_device() {
  static RTCDevice device = nullptr;
  if (!device) {
    device = rtcNewDevice("");
    rtcSetDeviceErrorFunction(
        device,
        [](void* ctx, RTCError code, const char* message) {
          auto str = string{message};
          switch (code) {
            case RTC_ERROR_UNKNOWN:
              throw std::runtime_error("RTC_ERROR_UNKNOWN: " + str);
            case RTC_ERROR_INVALID_ARGUMENT:
              throw std::runtime_error("RTC_ERROR_INVALID_ARGUMENT: " + str);
            case RTC_ERROR_INVALID_OPERATION:
              throw std::runtime_error("RTC_ERROR_INVALID_OPERATION: " + str);
            case RTC_ERROR_OUT_OF_MEMORY:
              throw std::runtime_error("RTC_ERROR_OUT_OF_MEMORY: " + str);
            case RTC_ERROR_UNSUPPORTED_CPU:
              throw std::runtime_error("RTC_ERROR_UNSUPPORTED_CPU: " + str);
            case RTC_ERROR_CANCELLED:
              throw std::runtime_error("RTC_ERROR_CANCELLED: " + str);
            default: throw std::runtime_error("invalid error code");
          }
        },
        nullptr);
    rtcSetDeviceMemoryMonitorFunction(
        device,
        [](void* userPtr, ssize_t bytes, bool post) {
          bvh_embree_memory += bytes;
          return true;
        },
        nullptr);
  }
  return device;
}

// Initialize Embree BVH
static void init_embree_bvh(
    bvh_shape& bvh, const scene_shape& shape, const bvh_params& params) {
  auto edevice   = bvh_embree_device();
  bvh.embree_bvh = rtcNewScene(edevice);
  auto escene    = (RTCScene)bvh.embree_bvh;
  if (params.bvh == bvh_build_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == bvh_build_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  if (!shape.points.empty()) {
    throw std::runtime_error("embree does not support points");
  } else if (!shape.lines.empty()) {
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto& l : shape.lines) {
      if (last_index == l.x) {
        elines.push_back((int)epositions.size() - 1);
        auto& posy = shape.positions[l.y];
        auto& rady = shape.radius[l.y];
        epositions.push_back({posy.x, posy.y, posy.z, rady});
      } else {
        elines.push_back((int)epositions.size());
        auto& posx = shape.positions[l.x];
        auto& radx = shape.radius[l.x];
        epositions.push_back({posx.x, posx.y, posx.z, radx});
        auto& posy = shape.positions[l.y];
        auto& rady = shape.radius[l.y];
        epositions.push_back({posy.x, posy.y, posy.z, rady});
      }
      last_index = l.y;
    }
    auto egeometry = rtcNewGeometry(
        edevice, RTC_GEOMETRY_TYPE_FLAT_LINEAR_CURVE);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, 4 * 4, epositions.size());
    auto embree_lines     = rtcSetNewGeometryBuffer(
        egeometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, 4, elines.size());
    memcpy(embree_positions, epositions.data(), epositions.size() * 16);
    memcpy(embree_lines, elines.data(), elines.size() * 4);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else if (!shape.triangles.empty()) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
        shape.positions.size());
    auto embree_triangles = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
        shape.triangles.size());
    memcpy(
        embree_positions, shape.positions.data(), shape.positions.size() * 12);
    memcpy(
        embree_triangles, shape.triangles.data(), shape.triangles.size() * 12);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else if (!shape.quads.empty()) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(egeometry, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
        shape.positions.size());
    auto embree_quads     = rtcSetNewGeometryBuffer(egeometry,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4, shape.quads.size());
    memcpy(
        embree_positions, shape.positions.data(), shape.positions.size() * 12);
    memcpy(embree_quads, shape.quads.data(), shape.quads.size() * 16);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, 0);
  } else {
    throw std::runtime_error("empty shapes not supported");
  }
  rtcCommitScene(escene);
}

static void init_embree_bvh(
    bvh_scene& bvh, const scene_scene& scene, const bvh_params& params) {
  // scene bvh
  auto edevice   = bvh_embree_device();
  bvh.embree_bvh = rtcNewScene(edevice);
  auto escene    = (RTCScene)bvh.embree_bvh;
  if (params.bvh == bvh_build_type::embree_compact)
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  if (params.bvh == bvh_build_type::embree_highquality)
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  for (auto instance_id = 0; instance_id < (int)scene.instances.size();
       instance_id++) {
    auto& instance  = scene.instances[instance_id];
    auto& sbvh      = bvh.shapes[instance.shape];
    auto  egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_INSTANCE);
    rtcSetGeometryInstancedScene(egeometry, (RTCScene)sbvh.embree_bvh);
    rtcSetGeometryTransform(
        egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, instance_id);
  }
  rtcCommitScene(escene);
}

static void update_embree_bvh(bvh_scene& bvh, const scene_scene& scene,
    const vector<int>& updated_instances) {
  // scene bvh
  auto escene = (RTCScene)bvh.embree_bvh;
  for (auto instance_id : updated_instances) {
    auto& instance    = scene.instances[instance_id];
    auto& sbvh        = bvh.shapes[instance.shape];
    auto  embree_geom = rtcGetGeometry(escene, instance_id);
    rtcSetGeometryInstancedScene(embree_geom, (RTCScene)sbvh.embree_bvh);
    rtcSetGeometryTransform(
        embree_geom, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(embree_geom);
  }
  rtcCommitScene(escene);
}

static bool intersect_embree_bvh(const bvh_shape& bvh, const scene_shape& shape,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any) {
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
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}

static bool intersect_embree_bvh(const bvh_scene& bvh, const scene_scene& scene,
    const ray3f& ray, int& instance, int& element, vec2f& uv, float& distance,
    bool find_any) {
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
  instance = (int)embree_ray.hit.instID[0];
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}
#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH CREATION
// -----------------------------------------------------------------------------
namespace yocto {

bvh_shape::bvh_shape(bvh_shape&& other) {
  bvh.nodes.swap(other.bvh.nodes);
  bvh.primitives.swap(other.bvh.primitives);
  std::swap(embree_bvh, other.embree_bvh);
}

bvh_scene::bvh_scene(bvh_scene&& other) {
  bvh.nodes.swap(other.bvh.nodes);
  bvh.primitives.swap(other.bvh.primitives);
  shapes.swap(other.shapes);
  std::swap(embree_bvh, other.embree_bvh);
}

bvh_shape::~bvh_shape() {
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene((RTCScene)embree_bvh);
#endif
}

bvh_scene::~bvh_scene() {
#ifdef YOCTO_EMBREE
  if (embree_bvh) rtcReleaseScene((RTCScene)embree_bvh);
#endif
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH BUILD
// -----------------------------------------------------------------------------
namespace yocto {

// Splits a BVH node using the SAH heuristic. Returns split position and axis.
static pair<int, int> split_sah(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto split_axis = 0;
  auto mid        = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, split_axis};

  // consider N bins, compute their cost and keep the minimum
  const int nbins    = 16;
  auto      middle   = 0.0f;
  auto      min_cost = flt_max;
  auto      area     = [](auto& b) {
    auto size = b.max - b.min;
    return 1e-12f + 2 * size.x * size.y + 2 * size.x * size.z +
           2 * size.y * size.z;
  };
  for (auto saxis = 0; saxis < 3; saxis++) {
    for (auto b = 1; b < nbins; b++) {
      auto split     = cbbox.min[saxis] + b * csize[saxis] / nbins;
      auto left_bbox = invalidb3f, right_bbox = invalidb3f;
      auto left_nprims = 0, right_nprims = 0;
      for (auto i = start; i < end; i++) {
        if (centers[primitives[i]][saxis] < split) {
          left_bbox = merge(left_bbox, bboxes[primitives[i]]);
          left_nprims += 1;
        } else {
          right_bbox = merge(right_bbox, bboxes[primitives[i]]);
          right_nprims += 1;
        }
      }
      auto cost = 1 + left_nprims * area(left_bbox) / area(cbbox) +
                  right_nprims * area(right_bbox) / area(cbbox);
      if (cost < min_cost) {
        min_cost   = cost;
        middle     = split;
        split_axis = saxis;
      }
    }
  }
  // split
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [split_axis, middle, &centers](
                      auto a) { return centers[a][split_axis] < middle; }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    split_axis = 0;
    mid        = (start + end) / 2;
    // throw std::runtime_error{"bad bvh split"};
  }

  return {mid, split_axis};
}

// Splits a BVH node using the balance heuristic. Returns split position and
// axis.
static pair<int, int> split_balanced(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // balanced tree split: find the largest axis of the
  // bounding box and split along this one right in the middle
  mid = (start + end) / 2;
  std::nth_element(primitives.data() + start, primitives.data() + mid,
      primitives.data() + end, [axis, &centers](auto a, auto b) {
        return centers[a][axis] < centers[b][axis];
      });

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    axis = 0;
    mid  = (start + end) / 2;
    // throw std::runtime_error("bad bvh split");
  }

  return {mid, axis};
}

// Splits a BVH node using the middle heutirtic. Returns split position and
// axis.
static pair<int, int> split_middle(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end) {
  // initialize split axis and position
  auto axis = 0;
  auto mid  = (start + end) / 2;

  // compute primintive bounds and size
  auto cbbox = invalidb3f;
  for (auto i = start; i < end; i++)
    cbbox = merge(cbbox, centers[primitives[i]]);
  auto csize = cbbox.max - cbbox.min;
  if (csize == zero3f) return {mid, axis};

  // split along largest
  if (csize.x >= csize.y && csize.x >= csize.z) axis = 0;
  if (csize.y >= csize.x && csize.y >= csize.z) axis = 1;
  if (csize.z >= csize.x && csize.z >= csize.y) axis = 2;

  // split the space in the middle along the largest axis
  auto cmiddle = (cbbox.max + cbbox.min) / 2;
  auto middle  = cmiddle[axis];
  mid = (int)(std::partition(primitives.data() + start, primitives.data() + end,
                  [axis, middle, &centers](
                      auto a) { return centers[a][axis] < middle; }) -
              primitives.data());

  // if we were not able to split, just break the primitives in half
  if (mid == start || mid == end) {
    axis = 0;
    mid  = (start + end) / 2;
    // throw std::runtime_error("bad bvh split");
  }

  return {mid, axis};
}

// Split bvh nodes according to a type
static pair<int, int> split_nodes(vector<int>& primitives,
    const vector<bbox3f>& bboxes, const vector<vec3f>& centers, int start,
    int end, bvh_build_type type) {
  switch (type) {
    case bvh_build_type::default_:
      return split_middle(primitives, bboxes, centers, start, end);
    case bvh_build_type::highquality:
      return split_sah(primitives, bboxes, centers, start, end);
    case bvh_build_type::middle:
      return split_middle(primitives, bboxes, centers, start, end);
    case bvh_build_type::balanced:
      return split_balanced(primitives, bboxes, centers, start, end);
    default: throw std::runtime_error("should not have gotten here");
  }
}

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// Build BVH nodes
static void build_bvh_serial(
    bvh_tree& bvh, const vector<bbox3f>& bboxes, const bvh_params& params) {
  // get values
  auto& nodes      = bvh.nodes;
  auto& primitives = bvh.primitives;

  // prepare to build nodes
  nodes.clear();
  nodes.reserve(bboxes.size() * 2);

  // prepare primitives
  bvh.primitives.resize(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) bvh.primitives[idx] = idx;

  // prepare centers
  auto centers = vector<vec3f>(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++)
    centers[idx] = center(bboxes[idx]);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)bboxes.size()}};
  nodes.emplace_back();

  // create nodes until the queue is empty
  while (!queue.empty()) {
    // grab node to work on
    auto next = queue.front();
    queue.pop_front();
    auto nodeid = next.x, start = next.y, end = next.z;

    // grab node
    auto& node = nodes[nodeid];

    // compute bounds
    node.bbox = invalidb3f;
    for (auto i = start; i < end; i++)
      node.bbox = merge(node.bbox, bboxes[primitives[i]]);

    // split into two children
    if (end - start > bvh_max_prims) {
      // get split
      auto [mid, axis] = split_nodes(
          primitives, bboxes, centers, start, end, params.bvh);

      // make an internal node
      node.internal = true;
      node.axis     = (uint8_t)axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = (int16_t)(end - start);
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

#if 0

// Build BVH nodes
static void build_bvh_parallel(
    bvh_tree_& bvh, const vector<bbox3f>& bboxes, bvh_build_type type) {
  // get values
  auto& nodes      = bvh.nodes;
  auto& primitives = bvh.primitives;

  // prepare to build nodes
  nodes.clear();
  nodes.reserve(bboxes.size() * 2);

  // prepare primitives
  bvh.primitives.resize(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) bvh.primitives[idx] = idx;

  // prepare centers
  auto centers = vector<vec3f>(bboxes.size());
  for (auto idx = 0; idx < bboxes.size(); idx++)
    centers[idx] = center(bboxes[idx]);

  // queue up first node
  auto queue = deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // synchronization
  atomic<int>          num_processed_prims(0);
  std::mutex           queue_mutex;
  vector<future<void>> futures;
  auto                 nthreads = std::thread::hardware_concurrency();

  // create nodes until the queue is empty
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(std::async(
        std::launch::async, [&nodes, &primitives, &bboxes, &centers, &type,
                                &num_processed_prims, &queue_mutex, &queue] {
          while (true) {
            // exit if needed
            if (num_processed_prims >= primitives.size()) return;

            // grab node to work on
            auto next = zero3i;
            {
              std::lock_guard<std::mutex> lock{queue_mutex};
              if (!queue.empty()) {
                next = queue.front();
                queue.pop_front();
              }
            }

            // wait a bit if needed
            if (next == zero3i) {
              std::this_thread::sleep_for(std::chrono::microseconds(10));
              continue;
            }

            // grab node
            auto  nodeid = next.x, start = next.y, end = next.z;
            auto& node = nodes[nodeid];

            // compute bounds
            node.bbox = invalidb3f;
            for (auto i = start; i < end; i++)
              node.bbox = merge(node.bbox, bboxes[primitives[i]]);

            // split into two children
            if (end - start > bvh_max_prims) {
              // get split
              auto [mid, axis] = split_nodes(
                  primitives, bboxes, centers, start, end, type);

              // make an internal node
              {
                std::lock_guard<std::mutex> lock{queue_mutex};
                node.internal = true;
                node.axis     = (int8_t)axis;
                node.num      = 2;
                node.start    = (int)nodes.size();
                nodes.emplace_back();
                nodes.emplace_back();
                queue.push_back({node.start + 0, start, mid});
                queue.push_back({node.start + 1, mid, end});
              }
            } else {
              // Make a leaf node
              node.internal = false;
              node.num      = (int16_t)(end - start);
              node.start    = start;
              num_processed_prims += node.num;
            }
          }
        }));
  }
  for (auto& f : futures) f.get();

  // cleanup
  nodes.shrink_to_fit();
}

#endif

// Update bvh
static void update_bvh(bvh_tree& bvh, const vector<bbox3f>& bboxes) {
  for (auto nodeid = (int)bvh.nodes.size() - 1; nodeid >= 0; nodeid--) {
    auto& node = bvh.nodes[nodeid];
    node.bbox  = invalidb3f;
    if (node.internal) {
      for (auto idx = 0; idx < 2; idx++) {
        node.bbox = merge(node.bbox, bvh.nodes[node.start + idx].bbox);
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        node.bbox = merge(node.bbox, bboxes[bvh.primitives[node.start + idx]]);
      }
    }
  }
}

static void build_bvh(
    bvh_shape& bvh, const scene_shape& shape, const bvh_params& params) {
#ifdef YOCTO_EMBREE
  if (params.bvh == bvh_build_type::embree_default ||
      params.bvh == bvh_build_type::embree_highquality ||
      params.bvh == bvh_build_type::embree_compact) {
    return init_embree_bvh(bvh, shape, params);
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>{};
  if (!shape.points.empty()) {
    bboxes = vector<bbox3f>(shape.points.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& p     = shape.points[idx];
      bboxes[idx] = point_bounds(shape.positions[p], shape.radius[p]);
    }
  } else if (!shape.lines.empty()) {
    bboxes = vector<bbox3f>(shape.lines.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& l     = shape.lines[idx];
      bboxes[idx] = line_bounds(shape.positions[l.x], shape.positions[l.y],
          shape.radius[l.x], shape.radius[l.y]);
    }
  } else if (!shape.triangles.empty()) {
    bboxes = vector<bbox3f>(shape.triangles.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& t     = shape.triangles[idx];
      bboxes[idx] = triangle_bounds(
          shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
    }
  } else if (!shape.quads.empty()) {
    bboxes = vector<bbox3f>(shape.quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape.quads[idx];
      bboxes[idx] = quad_bounds(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w]);
    }
  }

  // build nodes
  build_bvh_serial(bvh.bvh, bboxes, params);
}

static void build_bvh(
    bvh_scene& bvh, const scene_scene& scene, const bvh_params& params) {
  // embree
#ifdef YOCTO_EMBREE
  if (params.bvh == bvh_build_type::embree_default ||
      params.bvh == bvh_build_type::embree_highquality ||
      params.bvh == bvh_build_type::embree_compact) {
    return init_embree_bvh(bvh, scene, params);
  }
#endif

  // instance bboxes
  auto bboxes = vector<bbox3f>(scene.instances.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& instance = scene.instances[idx];
    auto& sbvh     = bvh.shapes[instance.shape];
    bboxes[idx]    = sbvh.bvh.nodes.empty()
                         ? invalidb3f
                         : transform_bbox(instance.frame, sbvh.bvh.nodes[0].bbox);
  }

  // build nodes
  build_bvh_serial(bvh.bvh, bboxes, params);
}

void init_bvh(bvh_shape& bvh, const scene_shape& shape,
    const bvh_params& params, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1};

  // build scene bvh
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
  build_bvh(bvh, shape, params);

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

void init_bvh(bvh_scene& bvh, const scene_scene& scene,
    const bvh_params& params, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1 + (int)scene.shapes.size()};

  // build shape bvh
  bvh.shapes.clear();
  for (auto idx = 0; idx < scene.shapes.size(); idx++) {
    if (progress_cb) progress_cb("build shape bvh", progress.x++, progress.y);
    bvh.shapes.emplace_back();
    build_bvh(bvh.shapes[idx], scene.shapes[idx], params);
  }

  // build scene bvh
  if (progress_cb) progress_cb("build scene bvh", progress.x++, progress.y);
  build_bvh(bvh, scene, params);

  // handle progress
  if (progress_cb) progress_cb("build bvh", progress.x++, progress.y);
}

static void update_bvh(bvh_shape& bvh, const scene_shape& shape) {
#ifdef YOCTO_EMBREE
  if (bvh.embree_bvh) {
    throw std::runtime_error("embree shape refit not supported");
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>{};
  if (!shape.points.empty()) {
    bboxes = vector<bbox3f>(shape.points.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& p     = shape.points[idx];
      bboxes[idx] = point_bounds(shape.positions[p], shape.radius[p]);
    }
  } else if (!shape.lines.empty()) {
    bboxes = vector<bbox3f>(shape.lines.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& l     = shape.lines[idx];
      bboxes[idx] = line_bounds(shape.positions[l.x], shape.positions[l.y],
          shape.radius[l.x], shape.radius[l.y]);
    }
  } else if (!shape.triangles.empty()) {
    bboxes = vector<bbox3f>(shape.triangles.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& t     = shape.triangles[idx];
      bboxes[idx] = triangle_bounds(
          shape.positions[t.x], shape.positions[t.y], shape.positions[t.z]);
    }
  } else if (!shape.quads.empty()) {
    bboxes = vector<bbox3f>(shape.quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape.quads[idx];
      bboxes[idx] = quad_bounds(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w]);
    }
  }

  // update nodes
  update_bvh(bvh.bvh, bboxes);
}

void update_bvh(bvh_scene& bvh, const scene_scene& scene,
    const vector<int>& updated_instances) {
#ifdef YOCTO_EMBREE
  if (bvh.embree_bvh) {
    return update_embree_bvh(bvh, scene, updated_instances);
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>(scene.instances.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& instance = scene.instances[idx];
    auto& sbvh     = bvh.shapes[instance.shape].bvh;
    bboxes[idx]    = transform_bbox(instance.frame, sbvh.nodes[0].bbox);
  }

  // update nodes
  update_bvh(bvh.bvh, bboxes);
}

void update_bvh(bvh_shape& bvh, const scene_shape& shape,
    const bvh_params& params, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1};

  // handle instances
  if (progress_cb) progress_cb("update bvh", progress.x++, progress.y);
  update_bvh(bvh, shape);

  // handle progress
  if (progress_cb) progress_cb("update bvh", progress.x++, progress.y);
}

void update_bvh(bvh_scene& bvh, const scene_scene& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes,
    const bvh_params& params, const progress_callback& progress_cb) {
  // handle progress
  auto progress = vec2i{0, 1 + (int)updated_shapes.size()};

  // update shapes
  for (auto shape : updated_shapes) {
    if (progress_cb) progress_cb("update shape bvh", progress.x++, progress.y);
    update_bvh(bvh.shapes[shape], scene.shapes[shape]);
  }

  // handle instances
  if (progress_cb) progress_cb("update scene bvh", progress.x++, progress.y);
  update_bvh(bvh, scene, updated_instances);

  // handle progress
  if (progress_cb) progress_cb("update bvh", progress.x++, progress.y);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect ray with a bvh.
static bool intersect_bvh(const bvh_shape& bvh, const scene_shape& shape,
    const ray3f& ray_, int& element, vec2f& uv, float& distance,
    bool find_any) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (bvh.embree_bvh) {
    return intersect_embree_bvh(
        bvh, shape, ray_, element, uv, distance, find_any);
  }
#endif

  // check empty
  if (bvh.bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 128>{};
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
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis] != 0) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else if (!shape.points.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& p = shape.points[bvh.bvh.primitives[idx]];
        if (intersect_point(
                ray, shape.positions[p], shape.radius[p], uv, distance)) {
          hit      = true;
          element  = bvh.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.lines.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& l = shape.lines[bvh.bvh.primitives[idx]];
        if (intersect_line(ray, shape.positions[l.x], shape.positions[l.y],
                shape.radius[l.x], shape.radius[l.y], uv, distance)) {
          hit      = true;
          element  = bvh.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.triangles.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& t = shape.triangles[bvh.bvh.primitives[idx]];
        if (intersect_triangle(ray, shape.positions[t.x], shape.positions[t.y],
                shape.positions[t.z], uv, distance)) {
          hit      = true;
          element  = bvh.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.quads.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape.quads[bvh.bvh.primitives[idx]];
        if (intersect_quad(ray, shape.positions[q.x], shape.positions[q.y],
                shape.positions[q.z], shape.positions[q.w], uv, distance)) {
          hit      = true;
          element  = bvh.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh.
static bool intersect_bvh(const bvh_scene& bvh, const scene_scene& scene,
    const ray3f& ray_, int& instance, int& element, vec2f& uv, float& distance,
    bool find_any, bool non_rigid_frames) {
#ifdef YOCTO_EMBREE
  // call Embree if needed
  if (bvh.embree_bvh) {
    return intersect_embree_bvh(
        bvh, scene, ray_, instance, element, uv, distance, find_any);
  }
#endif

  // check empty
  if (bvh.bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 128>{};
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
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis] != 0) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& instance_ = scene.instances[bvh.bvh.primitives[idx]];
        auto& shape     = get_shape(scene, instance_.shape);
        auto& sbvh      = bvh.shapes[instance_.shape];
        auto  inv_ray   = transform_ray(
            inverse(instance_.frame, non_rigid_frames), ray);
        if (intersect_bvh(
                sbvh, shape, inv_ray, element, uv, distance, find_any)) {
          hit      = true;
          instance = bvh.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh.
static bool intersect_bvh(const bvh_scene& bvh, const scene_scene& scene,
    int instance_, const ray3f& ray, int& element, vec2f& uv, float& distance,
    bool find_any, bool non_rigid_frames) {
  auto& instance = scene.instances[instance_];
  auto  shape    = get_shape(scene, instance.shape);
  auto& sbvh     = bvh.shapes[instance.shape];
  auto  inv_ray = transform_ray(inverse(instance.frame, non_rigid_frames), ray);
  return intersect_bvh(sbvh, shape, inv_ray, element, uv, distance, find_any);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH OVERLAP
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect ray with a bvh.
static bool overlap_bvh(const bvh_shape& bvh, const scene_shape& shape,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  // check if empty
  if (bvh.bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 64>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else if (!shape.points.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = bvh.bvh.primitives[node.start + idx];
        auto& p         = shape.points[primitive];
        if (overlap_point(pos, max_distance, shape.positions[p],
                shape.radius[p], uv, distance)) {
          hit          = true;
          element      = primitive;
          max_distance = distance;
        }
      }
    } else if (!shape.lines.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = bvh.bvh.primitives[node.start + idx];
        auto& l         = shape.lines[primitive];
        if (overlap_line(pos, max_distance, shape.positions[l.x],
                shape.positions[l.y], shape.radius[l.x], shape.radius[l.y], uv,
                distance)) {
          hit          = true;
          element      = primitive;
          max_distance = distance;
        }
      }
    } else if (!shape.triangles.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = bvh.bvh.primitives[node.start + idx];
        auto& t         = shape.triangles[primitive];
        if (overlap_triangle(pos, max_distance, shape.positions[t.x],
                shape.positions[t.y], shape.positions[t.z], shape.radius[t.x],
                shape.radius[t.y], shape.radius[t.z], uv, distance)) {
          hit          = true;
          element      = primitive;
          max_distance = distance;
        }
      }
    } else if (!shape.quads.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = bvh.bvh.primitives[node.start + idx];
        auto& q         = shape.quads[primitive];
        if (overlap_quad(pos, max_distance, shape.positions[q.x],
                shape.positions[q.y], shape.positions[q.z],
                shape.positions[q.w], shape.radius[q.x], shape.radius[q.y],
                shape.radius[q.z], shape.radius[q.w], uv, distance)) {
          hit          = true;
          element      = primitive;
          max_distance = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

// Intersect ray with a bvh.
static bool overlap_bvh(const bvh_scene& bvh, const scene_scene& scene,
    const vec3f& pos, float max_distance, int& instance, int& element,
    vec2f& uv, float& distance, bool find_any, bool non_rigid_frames) {
  // check if empty
  if (bvh.bvh.nodes.empty()) return false;

  // node stack
  auto node_stack        = array<int, 64>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = bvh.bvh.primitives[node.start + idx];
        auto& instance_ = scene.instances[primitive];
        auto& shape     = scene.shapes[instance_.shape];
        auto& sbvh      = bvh.shapes[instance_.shape];
        auto  inv_pos   = transform_point(
            inverse(instance_.frame, non_rigid_frames), pos);
        if (overlap_bvh(sbvh, shape, inv_pos, max_distance, element, uv,
                distance, find_any)) {
          hit          = true;
          instance     = primitive;
          max_distance = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

#if 0
// Finds the overlap between BVH leaf nodes.
template <typename OverlapElem>
void overlap_bvh_elems(const bvh_data& bvh1, const bvh_data& bvh2,
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

bvh_intersection intersect_bvh(const bvh_scene& bvh, const scene_scene& scene,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection = bvh_intersection{};
  intersection.hit  = intersect_bvh(bvh, scene, ray, intersection.instance,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  return intersection;
}
bvh_intersection intersect_bvh(const bvh_scene& bvh, const scene_scene& scene,
    int instance, const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection     = bvh_intersection{};
  intersection.hit      = intersect_bvh(bvh, scene, instance, ray,
      intersection.element, intersection.uv, intersection.distance, find_any,
      non_rigid_frames);
  intersection.instance = instance;
  return intersection;
}

bvh_intersection overlap_bvh(const bvh_scene& bvh, const scene_scene& scene,
    const vec3f& pos, float max_distance, bool find_any,
    bool non_rigid_frames) {
  auto intersection = bvh_intersection{};
  intersection.hit  = overlap_bvh(bvh, scene, pos, max_distance,
      intersection.instance, intersection.element, intersection.uv,
      intersection.distance, find_any, non_rigid_frames);
  return intersection;
}

}  // namespace yocto
