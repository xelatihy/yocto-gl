//
// Implementation for Yocto/BVH
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
//

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_bvh.h"

#include <algorithm>
#include <atomic>
#include <deque>
#include <future>
#include <thread>

#if YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EMBREE BVH
// -----------------------------------------------------------------------------
namespace yocto {

#if YOCTO_EMBREE
// copies
bvh_embree::bvh_embree(const bvh_embree& other) { *this = other; }
bvh_embree& bvh_embree::operator=(const bvh_embree& other) {
  shape  = other.shape;
  scene  = other.scene;
  device = other.device;
  if (shape) rtcRetainGeometry(shape);
  if (scene) rtcRetainScene(scene);
  if (device) rtcRetainDevice(device);
  return *this;
}
// cleanup
bvh_embree::~bvh_embree() {
  if (shape) rtcReleaseGeometry(shape);
  if (scene) rtcReleaseScene(scene);
  if (device) rtcReleaseDevice(device);
}

// Get Embree device
std::atomic<ssize_t> bvh_embree_memory = 0;
static RTCDevice     bvh_embree_device() {
  static RTCDevice device = nullptr;
  if (!device) {
    device = rtcNewDevice("");
    rtcSetDeviceErrorFunction(
        device,
        [](void* ctx, RTCError code, const char* str) {
          switch (code) {
            case RTC_ERROR_UNKNOWN:
              throw std::runtime_error("RTC_ERROR_UNKNOWN: "s + str);
            case RTC_ERROR_INVALID_ARGUMENT:
              throw std::runtime_error("RTC_ERROR_INVALID_ARGUMENT: "s + str);
            case RTC_ERROR_INVALID_OPERATION:
              throw std::runtime_error("RTC_ERROR_INVALID_OPERATION: "s + str);
            case RTC_ERROR_OUT_OF_MEMORY:
              throw std::runtime_error("RTC_ERROR_OUT_OF_MEMORY: "s + str);
            case RTC_ERROR_UNSUPPORTED_CPU:
              throw std::runtime_error("RTC_ERROR_UNSUPPORTED_CPU: "s + str);
            case RTC_ERROR_CANCELLED:
              throw std::runtime_error("RTC_ERROR_CANCELLED: "s + str);
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
  rtcRetainDevice(device);
  return device;
}

// Initialize Embree BVH
void init_shape_embree_bvh(bvh_shape& shape, const bvh_params& params) {
  shape.embree.device = bvh_embree_device();
  shape.embree.scene  = rtcNewScene(shape.embree.device);
  if (params.compact)
    rtcSetSceneFlags(shape.embree.scene, RTC_SCENE_FLAG_COMPACT);
  if (params.high_quality)
    rtcSetSceneBuildQuality(shape.embree.scene, RTC_BUILD_QUALITY_HIGH);
  if (!shape.points.empty()) {
    throw std::runtime_error("embree does not support points");
  } else if (!shape.lines.empty()) {
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto& l : shape.lines) {
      if (last_index == l.x) {
        elines.push_back((int)epositions.size() - 1);
        epositions.push_back({shape.positions[l.y], shape.radius[l.y]});
      } else {
        elines.push_back((int)epositions.size());
        epositions.push_back({shape.positions[l.x], shape.radius[l.x]});
        epositions.push_back({shape.positions[l.y], shape.radius[l.y]});
      }
      last_index = l.y;
    }
    shape.embree.shape = rtcNewGeometry(
        shape.embree.device, RTC_GEOMETRY_TYPE_FLAT_LINEAR_CURVE);
    rtcSetGeometryVertexAttributeCount(shape.embree.shape, 1);
    auto embree_positions = rtcSetNewGeometryBuffer(shape.embree.shape,
        RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT4, 4 * 4, epositions.size());
    auto embree_lines     = rtcSetNewGeometryBuffer(shape.embree.shape,
        RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT, 4, elines.size());
    memcpy(embree_positions, epositions.data(), epositions.size() * 16);
    memcpy(embree_lines, elines.data(), elines.size() * 4);
  } else if (!shape.triangles.empty()) {
    shape.embree.shape = rtcNewGeometry(
        shape.embree.device, RTC_GEOMETRY_TYPE_TRIANGLE);
    rtcSetGeometryVertexAttributeCount(shape.embree.shape, 1);
    if (params.compact) {
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape.positions.data(), 0, 3 * 4,
          shape.positions.size());
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT3, shape.triangles.data(), 0, 3 * 4,
          shape.triangles.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape.positions.size());
      auto embree_triangles = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * 4,
          shape.triangles.size());
      memcpy(embree_positions, shape.positions.data(),
          shape.positions.size() * 12);
      memcpy(embree_triangles, shape.triangles.data(),
          shape.triangles.size() * 12);
    }
  } else if (!shape.quads.empty()) {
    shape.embree.shape = rtcNewGeometry(
        shape.embree.device, RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(shape.embree.shape, 1);
    if (params.compact) {
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape.positions.data(), 0, 3 * 4,
          shape.positions.size());
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT4, shape.quads.data(), 0, 4 * 4, shape.quads.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape.positions.size());
      auto embree_quads     = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4,
          shape.quads.size());
      memcpy(embree_positions, shape.positions.data(),
          shape.positions.size() * 12);
      memcpy(embree_quads, shape.quads.data(), shape.quads.size() * 16);
    }
  } else if (!shape.quadspos.empty()) {
    shape.embree.shape = rtcNewGeometry(
        shape.embree.device, RTC_GEOMETRY_TYPE_QUAD);
    rtcSetGeometryVertexAttributeCount(shape.embree.shape, 1);
    if (params.compact) {
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_VERTEX, 0,
          RTC_FORMAT_FLOAT3, shape.positions.data(), 0, 3 * 4,
          shape.positions.size());
      rtcSetSharedGeometryBuffer(shape.embree.shape, RTC_BUFFER_TYPE_INDEX, 0,
          RTC_FORMAT_UINT4, shape.quadspos.data(), 0, 4 * 4,
          shape.quadspos.size());
    } else {
      auto embree_positions = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * 4,
          shape.positions.size());
      auto embree_quads     = rtcSetNewGeometryBuffer(shape.embree.shape,
          RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT4, 4 * 4,
          shape.quadspos.size());
      memcpy(embree_positions, shape.positions.data(),
          shape.positions.size() * 12);
      memcpy(embree_quads, shape.quadspos.data(), shape.quadspos.size() * 16);
    }
  } else {
    throw std::runtime_error("empty shapes not supported");
  }
  rtcCommitGeometry(shape.embree.shape);
  rtcAttachGeometryByID(shape.embree.scene, shape.embree.shape, 0);
  rtcCommitScene(shape.embree.scene);
}
void init_scene_embree_bvh(bvh_scene& scene, const bvh_params& params) {
  // scene bvh
  scene.embree.device = bvh_embree_device();
  scene.embree.scene  = rtcNewScene(scene.embree.device);
  if (params.compact)
    rtcSetSceneFlags(scene.embree.scene, RTC_SCENE_FLAG_COMPACT);
  if (params.high_quality)
    rtcSetSceneBuildQuality(scene.embree.scene, RTC_BUILD_QUALITY_HIGH);
  for (auto instance_id = 0; instance_id < scene.instances.size();
       instance_id++) {
    auto& instance = scene.instances[instance_id];
    auto& sbvh     = scene.shapes[instance.shape].embree;
    if (!sbvh.scene) throw std::runtime_error("bvh not built");
    scene.embree.instances.push_back(
        rtcNewGeometry(scene.embree.device, RTC_GEOMETRY_TYPE_INSTANCE));
    rtcSetGeometryInstancedScene(scene.embree.instances.back(), sbvh.scene);
    rtcSetGeometryTransform(scene.embree.instances.back(), 0,
        RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(scene.embree.instances.back());
    rtcAttachGeometryByID(
        scene.embree.scene, scene.embree.instances.back(), instance_id);
  }
  rtcCommitScene(scene.embree.scene);
}

void update_scene_embree_bvh(
    bvh_scene& scene, const vector<int>& updated_instances) {
  // scene bvh
  for (auto instance_id : updated_instances) {
    auto& instance = scene.instances[instance_id];
    auto& sbvh     = scene.shapes[instance.shape].embree;
    if (!sbvh.scene) throw std::runtime_error("bvh not built");
    auto embree_geom = rtcGetGeometry(scene.embree.scene, instance_id);
    rtcSetGeometryInstancedScene(embree_geom, sbvh.scene);
    rtcSetGeometryTransform(
        embree_geom, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(embree_geom);
  }
  rtcCommitScene(scene.embree.scene);
}

bool intersect_shape_embree_bvh(const bvh_shape& shape, const ray3f& ray,
    int& element, vec2f& uv, float& distance, bool find_any) {
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
  rtcIntersect1(shape.embree.scene, &embree_ctx, &embree_ray);
  if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return false;
  element  = (int)embree_ray.hit.primID;
  uv       = {embree_ray.hit.u, embree_ray.hit.v};
  distance = embree_ray.ray.tfar;
  return true;
}
bool intersect_scene_embree_bvh(const bvh_scene& scene, const ray3f& ray,
    int& instance, int& element, vec2f& uv, float& distance, bool find_any) {
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
  rtcIntersect1(scene.embree.scene, &embree_ctx, &embree_ray);
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
// IMPLEMENTATION FOR BVH
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
    throw std::runtime_error("bad bvh split");
    split_axis = 0;
    mid        = (start + end) / 2;
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
    throw std::runtime_error("bad bvh split");
    axis = 0;
    mid  = (start + end) / 2;
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
    throw std::runtime_error("bad bvh split");
    axis = 0;
    mid  = (start + end) / 2;
  }

  return {mid, axis};
}

// Build BVH nodes
static void build_bvh_serial(
    bvh_tree& bvh, vector<bbox3f>& bboxes, bool high_quality) {
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
  auto queue = std::deque<vec3i>{{0, 0, (int)bboxes.size()}};
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
      auto [mid, axis] =
          high_quality
              ? split_sah(primitives, bboxes, centers, start, end)
              : split_balanced(primitives, bboxes, centers, start, end);

      // make an internal node
      node.internal = true;
      node.axis     = axis;
      node.num      = 2;
      node.start    = (int)nodes.size();
      nodes.emplace_back();
      nodes.emplace_back();
      queue.push_back({node.start + 0, start, mid});
      queue.push_back({node.start + 1, mid, end});
    } else {
      // Make a leaf node
      node.internal = false;
      node.num      = end - start;
      node.start    = start;
    }
  }

  // cleanup
  nodes.shrink_to_fit();
}

// Build BVH nodes
static void build_bvh_parallel(
    bvh_tree& bvh, vector<bbox3f>& bboxes, bool high_quality) {
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
  auto queue = std::deque<vec3i>{{0, 0, (int)primitives.size()}};
  nodes.emplace_back();

  // synchronization
  std::atomic<int>          num_processed_prims(0);
  std::mutex                queue_mutex;
  vector<std::future<void>> futures;
  auto                      nthreads = std::thread::hardware_concurrency();

  // create nodes until the queue is empty
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(std::async(std::launch::async,
        [&nodes, &primitives, &bboxes, &centers, &high_quality,
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
              auto [mid, axis] =
                  high_quality
                      ? split_sah(primitives, bboxes, centers, start, end)
                      : split_balanced(primitives, bboxes, centers, start, end);

              // make an internal node
              {
                std::lock_guard<std::mutex> lock{queue_mutex};
                node.internal = true;
                node.axis     = axis;
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
              node.num      = end - start;
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

// Build shape bvh
void make_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // build nodes
  if (!parallel) {
    build_bvh_serial(bvh, bboxes, high_quality);
  } else {
    build_bvh_parallel(bvh, bboxes, high_quality);
  }
}
void make_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& l     = lines[idx];
    bboxes[idx] = line_bounds(
        positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
  }

  // build nodes
  if (!parallel) {
    build_bvh_serial(bvh, bboxes, high_quality);
  } else {
    build_bvh_parallel(bvh, bboxes, high_quality);
  }
}
void make_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& t     = triangles[idx];
    bboxes[idx] = triangle_bounds(
        positions[t.x], positions[t.y], positions[t.z]);
  }

  // build nodes
  if (!parallel) {
    build_bvh_serial(bvh, bboxes, high_quality);
  } else {
    build_bvh_parallel(bvh, bboxes, high_quality);
  }
}
void make_quads_bvh(bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius,
    bool high_quality, bool parallel) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& q     = quads[idx];
    bboxes[idx] = quad_bounds(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
  }

  // build nodes
  if (!parallel) {
    build_bvh_serial(bvh, bboxes, high_quality);
  } else {
    build_bvh_parallel(bvh, bboxes, high_quality);
  }
}

void update_points_bvh(bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(points.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& p     = points[idx];
    bboxes[idx] = point_bounds(positions[p], radius[p]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_lines_bvh(bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius) {
  // build primitives
  auto bboxes = vector<bbox3f>(lines.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& l     = lines[idx];
    bboxes[idx] = line_bounds(
        positions[l.x], positions[l.y], radius[l.x], radius[l.y]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_triangles_bvh(bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(triangles.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& t     = triangles[idx];
    bboxes[idx] = triangle_bounds(
        positions[t.x], positions[t.y], positions[t.z]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}
void update_quads_bvh(
    bvh_tree& bvh, const vector<vec4i>& quads, const vector<vec3f>& positions) {
  // build primitives
  auto bboxes = vector<bbox3f>(quads.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& q     = quads[idx];
    bboxes[idx] = quad_bounds(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
  }

  // update nodes
  update_bvh(bvh, bboxes);
}

// Intersect ray with a bvh.
template <typename Intersect>
static bool intersect_elements_bvh(const bvh_tree& bvh,
    Intersect&& intersect_element, const ray3f& ray_, int& element, vec2f& uv,
    float& distance, bool find_any) {
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
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        if (intersect_element(primitive, ray, uv, distance)) {
          hit      = true;
          element  = primitive;
          ray.tmax = distance;
        }
      }
    }

    // check for early exit
    if (find_any && hit) return hit;
  }

  return hit;
}

template <typename Frame, typename Intersect>
static bool intersect_elements_bvh(const bvh_tree& bvh, Frame&& instance_frame,
    Intersect&& intersect_shape, const ray3f& ray_, int& instance, int& element,
    vec2f& uv, float& distance, bool find_any, bool non_rigid_frames) {
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
  auto ray_dinv  = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
  auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
      (ray_dinv.z < 0) ? 1 : 0};

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        auto inv_ray   = transform_ray(
            inverse(instance_frame(primitive), non_rigid_frames), ray);
        if (intersect_shape(
                primitive, inv_ray, element, uv, distance, find_any)) {
          hit      = true;
          instance = primitive;
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
bool intersect_points_bvh(const bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&points, &positions, &radius](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& p = points[idx];
        return intersect_point(ray, positions[p], radius[p], uv, distance);
      },
      ray, element, uv, distance, find_any);
}
bool intersect_lines_bvh(const bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&lines, &positions, &radius](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& l = lines[idx];
        return intersect_line(ray, positions[l.x], positions[l.y], radius[l.x],
            radius[l.y], uv, distance);
      },
      ray, element, uv, distance, find_any);
}
bool intersect_triangles_bvh(const bvh_tree& bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&triangles, &positions](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& t = triangles[idx];
        return intersect_triangle(
            ray, positions[t.x], positions[t.y], positions[t.z], uv, distance);
      },
      ray, element, uv, distance, find_any);
}
bool intersect_quads_bvh(const bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const ray3f& ray, int& element, vec2f& uv,
    float& distance, bool find_any) {
  return intersect_elements_bvh(
      bvh,
      [&quads, &positions](
          int idx, const ray3f& ray, vec2f& uv, float& distance) {
        auto& t = quads[idx];
        return intersect_quad(ray, positions[t.x], positions[t.y],
            positions[t.z], positions[t.w], uv, distance);
      },
      ray, element, uv, distance, find_any);
}

// Intersect ray with a bvh.
template <typename Overlap>
static bool overlap_elements_bvh(const bvh_tree& bvh, Overlap&& overlap_element,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  // check if empty
  if (bvh.nodes.empty()) return false;

  // node stack
  int  node_stack[64];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        if (overlap_element(primitive, pos, max_distance, uv, distance)) {
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
template <typename Frame, typename Overlap>
bool overlap_elements_bvh(const bvh_tree& bvh, Frame&& instance_frame,
    Overlap&& overlap_shape, const vec3f& pos, float max_distance,
    int& instance, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  // check if empty
  if (bvh.nodes.empty()) return false;

  // node stack
  int  node_stack[64];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = bvh.primitives[node.start + idx];
        auto inv_pos   = transform_point(
            inverse(instance_frame(primitive), non_rigid_frames), pos);
        if (overlap_shape(primitive, inv_pos, max_distance, element, uv,
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

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bool overlap_points_bvh(const bvh_tree& bvh, const vector<int>& points,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&points, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& p = points[idx];
        return overlap_point(
            pos, max_distance, positions[p], radius[p], uv, distance);
      },
      pos, max_distance, element, uv, distance, find_any);
}
bool overlap_lines_bvh(const bvh_tree& bvh, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&lines, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& l = lines[idx];
        return overlap_line(pos, max_distance, positions[l.x], positions[l.y],
            radius[l.x], radius[l.y], uv, distance);
      },
      pos, max_distance, element, uv, distance, find_any);
}
bool overlap_triangles_bvh(const bvh_tree& bvh, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&triangles, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& t = triangles[idx];
        return overlap_triangle(pos, max_distance, positions[t.x],
            positions[t.y], positions[t.z], radius[t.x], radius[t.y],
            radius[t.z], uv, distance);
      },
      pos, max_distance, element, uv, distance, find_any);
}
bool overlap_quads_bvh(const bvh_tree& bvh, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<float>& radius,
    const vec3f& pos, float max_distance, int& element, vec2f& uv,
    float& distance, bool find_any) {
  return overlap_elements_bvh(
      bvh,
      [&quads, &positions, &radius](int idx, const vec3f& pos,
          float max_distance, vec2f& uv, float& distance) {
        auto& q = quads[idx];
        return overlap_quad(pos, max_distance, positions[q.x], positions[q.y],
            positions[q.z], positions[q.w], radius[q.x], radius[q.y],
            radius[q.z], radius[q.w], uv, distance);
      },
      pos, max_distance, element, uv, distance, find_any);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE/SCENE BVH
// -----------------------------------------------------------------------------
namespace yocto {

void init_shape_bvh(bvh_shape& shape, const bvh_params& params) {
#if YOCTO_EMBREE
  // call Embree if needed
  if (params.embree) {
    return init_shape_embree_bvh(shape, params);
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
  } else if (!shape.quadspos.empty()) {
    bboxes = vector<bbox3f>(shape.quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape.quads[idx];
      bboxes[idx] = quad_bounds(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w]);
    }
  }

  // build nodes
  if (params.noparallel) {
    build_bvh_serial(shape.bvh, bboxes, params.high_quality);
  } else {
    build_bvh_parallel(shape.bvh, bboxes, params.high_quality);
  }
}

void init_scene_bvh(bvh_scene& scene, const bvh_params& params) {
  // Make shape bvh
  for (auto idx = 0; idx < scene.shapes.size(); idx++) {
    init_shape_bvh(scene.shapes[idx], params);
  }

  // embree
#if YOCTO_EMBREE
  if (params.embree) {
    return init_scene_embree_bvh(scene, params);
  }
#endif

  // instance bboxes
  auto bboxes = vector<bbox3f>(scene.instances.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& instance = scene.instances[idx];
    auto& shape    = scene.shapes[instance.shape];
    bboxes[idx]    = shape.bvh.nodes.empty()
                      ? invalidb3f
                      : transform_bbox(instance.frame, shape.bvh.nodes[0].bbox);
  }

  // build nodes
  if (params.noparallel) {
    build_bvh_serial(scene.bvh, bboxes, params.high_quality);
  } else {
    build_bvh_parallel(scene.bvh, bboxes, params.high_quality);
  }
}

void update_shape_bvh(bvh_shape& shape, const bvh_params& params) {
#if YOCTO_EMBREE
  if (params.embree) {
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
  } else if (!shape.quadspos.empty()) {
    bboxes = vector<bbox3f>(shape.quads.size());
    for (auto idx = 0; idx < bboxes.size(); idx++) {
      auto& q     = shape.quads[idx];
      bboxes[idx] = quad_bounds(shape.positions[q.x], shape.positions[q.y],
          shape.positions[q.z], shape.positions[q.w]);
    }
  }

  // update nodes
  update_bvh(shape.bvh, bboxes);
}

void update_scene_bvh(bvh_scene& scene, const vector<int>& updated_instances,
    const vector<int>& updated_shapes, const bvh_params& params) {
  // update shapes
  for (auto shape : updated_shapes)
    update_shape_bvh(scene.shapes[shape], params);

#if YOCTO_EMBREE
  if (params.embree) {
    update_scene_embree_bvh(scene, updated_instances);
  }
#endif

  // build primitives
  auto bboxes = vector<bbox3f>(scene.instances.size());
  for (auto idx = 0; idx < bboxes.size(); idx++) {
    auto& instance = scene.instances[idx];
    auto& sbvh     = scene.shapes[instance.shape].bvh;
    bboxes[idx]    = transform_bbox(instance.frame, sbvh.nodes[0].bbox);
  }

  // update nodes
  update_bvh(scene.bvh, bboxes);
}

// Intersect ray with a bvh.
bool intersect_shape_bvh(const bvh_shape& shape, const ray3f& ray_,
    int& element, vec2f& uv, float& distance, bool find_any) {
#if YOCTO_EMBREE
  // call Embree if needed
  if (shape.embree.scene) {
    return intersect_shape_embree_bvh(
        shape, ray_, element, uv, distance, find_any);
  }
#endif

  // check empty
  if (shape.bvh.nodes.empty()) return false;

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
    auto& node = shape.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else if (!shape.points.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& p = shape.points[shape.bvh.primitives[idx]];
        if (intersect_point(
                ray, shape.positions[p], shape.radius[p], uv, distance)) {
          hit      = true;
          element  = shape.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.lines.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& l = shape.lines[shape.bvh.primitives[idx]];
        if (intersect_line(ray, shape.positions[l.x], shape.positions[l.y],
                shape.radius[l.x], shape.radius[l.y], uv, distance)) {
          hit      = true;
          element  = shape.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.triangles.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& t = shape.triangles[shape.bvh.primitives[idx]];
        if (intersect_triangle(ray, shape.positions[t.x], shape.positions[t.y],
                shape.positions[t.z], uv, distance)) {
          hit      = true;
          element  = shape.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.quads.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape.quads[shape.bvh.primitives[idx]];
        if (intersect_quad(ray, shape.positions[q.x], shape.positions[q.y],
                shape.positions[q.z], shape.positions[q.w], uv, distance)) {
          hit      = true;
          element  = shape.bvh.primitives[idx];
          ray.tmax = distance;
        }
      }
    } else if (!shape.quadspos.empty()) {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& q = shape.quadspos[shape.bvh.primitives[idx]];
        if (intersect_quad(ray, shape.positions[q.x], shape.positions[q.y],
                shape.positions[q.z], shape.positions[q.w], uv, distance)) {
          hit      = true;
          element  = shape.bvh.primitives[idx];
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
bool intersect_scene_bvh(const bvh_scene& scene, const ray3f& ray_,
    int& instance, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
#if YOCTO_EMBREE
  // call Embree if needed
  if (scene.embree.scene) {
    return intersect_scene_embree_bvh(
        scene, ray_, instance, element, uv, distance, find_any);
  }
#endif

  // check empty
  if (scene.bvh.nodes.empty()) return false;

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
    auto& node = scene.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    // if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;
    if (!intersect_bbox(ray, ray_dinv, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // for internal nodes, attempts to proceed along the
      // split axis from smallest to largest nodes
      if (ray_dsign[node.axis]) {
        node_stack[node_cur++] = node.start + 0;
        node_stack[node_cur++] = node.start + 1;
      } else {
        node_stack[node_cur++] = node.start + 1;
        node_stack[node_cur++] = node.start + 0;
      }
    } else {
      for (auto idx = node.start; idx < node.start + node.num; idx++) {
        auto& instance_ = scene.instances[scene.bvh.primitives[idx]];
        auto  inv_ray   = transform_ray(
            inverse(instance_.frame, non_rigid_frames), ray);
        if (intersect_shape_bvh(scene.shapes[instance_.shape], inv_ray, element,
                uv, distance, find_any)) {
          hit      = true;
          instance = scene.bvh.primitives[idx];
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
bool intersect_instance_bvh(const bvh_scene& scene, int instance,
    const ray3f& ray, int& element, vec2f& uv, float& distance, bool find_any,
    bool non_rigid_frames) {
  auto& instance_ = scene.instances[instance];
  auto inv_ray = transform_ray(inverse(instance_.frame, non_rigid_frames), ray);
  return intersect_shape_bvh(
      scene.shapes[instance_.shape], inv_ray, element, uv, distance, find_any);
}

// Intersect ray with a bvh.
bool overlap_shape_bvh(const bvh_shape& shape, const vec3f& pos,
    float max_distance, int& element, vec2f& uv, float& distance,
    bool find_any) {
  // check if empty
  if (shape.bvh.nodes.empty()) return false;

  // node stack
  int  node_stack[64];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = shape.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else if (!shape.points.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = shape.bvh.primitives[node.start + idx];
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
        auto  primitive = shape.bvh.primitives[node.start + idx];
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
        auto  primitive = shape.bvh.primitives[node.start + idx];
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
        auto  primitive = shape.bvh.primitives[node.start + idx];
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
    } else if (!shape.quadspos.empty()) {
      for (auto idx = 0; idx < node.num; idx++) {
        auto  primitive = shape.bvh.primitives[node.start + idx];
        auto& q         = shape.quadspos[primitive];
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
bool overlap_scene_bvh(const bvh_scene& scene, const vec3f& pos,
    float max_distance, int& instance, int& element, vec2f& uv, float& distance,
    bool find_any, bool non_rigid_frames) {
  // check if empty
  if (scene.bvh.nodes.empty()) return false;

  // node stack
  int  node_stack[64];
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // hit
  auto hit = false;

  // walking stack
  while (node_cur) {
    // grab node
    auto& node = scene.bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!distance_check_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx = 0; idx < node.num; idx++) {
        auto primitive = scene.bvh.primitives[node.start + idx];
        auto instance_ = scene.instances[primitive];
        auto inv_pos   = transform_point(
            inverse(instance_.frame, non_rigid_frames), pos);
        if (overlap_shape_bvh(scene.shapes[instance_.shape], inv_pos,
                max_distance, element, uv, distance, find_any)) {
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

bvh_intersection intersect_shape_bvh(
    const bvh_shape& shape, const ray3f& ray, bool find_any) {
  auto intersection = bvh_intersection{};
  intersection.hit  = intersect_shape_bvh(shape, ray, intersection.element,
      intersection.uv, intersection.distance, find_any);
  return intersection;
}
bvh_intersection intersect_scene_bvh(const bvh_scene& scene, const ray3f& ray,
    bool find_any, bool non_rigid_frames) {
  auto intersection = bvh_intersection{};
  intersection.hit  = intersect_scene_bvh(scene, ray, intersection.instance,
      intersection.element, intersection.uv, intersection.distance, find_any);
  return intersection;
}
bvh_intersection intersect_instance_bvh(const bvh_scene& scene, int instance,
    const ray3f& ray, bool find_any, bool non_rigid_frames) {
  auto intersection     = bvh_intersection{};
  intersection.hit      = intersect_instance_bvh(scene, instance, ray,
      intersection.element, intersection.uv, intersection.distance, find_any);
  intersection.instance = instance;
  return intersection;
}

bvh_intersection overlap_shape_bvh(const bvh_shape& shape, const vec3f& pos,
    float max_distance, bool find_any) {
  auto intersection = bvh_intersection{};
  intersection.hit  = overlap_shape_bvh(shape, pos, max_distance,
      intersection.element, intersection.uv, intersection.distance, find_any);
  return intersection;
}
bvh_intersection overlap_scene_bvh(const bvh_scene& scene, const vec3f& pos,
    float max_distance, bool find_any, bool non_rigid_frames) {
  auto intersection = bvh_intersection{};
  intersection.hit  = overlap_scene_bvh(scene, pos, max_distance,
      intersection.instance, intersection.element, intersection.uv,
      intersection.distance, find_any);
  return intersection;
}

}  // namespace yocto
