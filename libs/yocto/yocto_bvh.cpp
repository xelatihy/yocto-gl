//
// Implementation for Yocto/Bvh
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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
#include <climits>
#include <cstring>
#include <future>
#include <stdexcept>

#include "yocto_geometry.h"

#ifdef YOCTO_EMBREE
#include <embree3/rtcore.h>
#endif

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;

// consts
inline const auto uint_max = std::numeric_limits<unsigned int>::max();

}  // namespace yocto

// -----------------------------------------------------------------------------
// PARALLEL HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Sequence1, typename Sequence2, typename Func>
inline void parallel_zip(
    Sequence1&& sequence1, Sequence2&& sequence2, Func&& func) {
  if (std::size(sequence1) != std::size(sequence2))
    throw std::out_of_range{"invalid sequence lengths"};
  auto                num      = std::size(sequence1);
  auto                futures  = vector<std::future<void>>{};
  auto                nthreads = std::thread::hardware_concurrency();
  std::atomic<size_t> next_idx(0);
  std::atomic<bool>   has_error(false);
  for (auto thread_id = 0; thread_id < (int)nthreads; thread_id++) {
    futures.emplace_back(std::async(std::launch::async,
        [&func, &next_idx, &has_error, num, &sequence1, &sequence2]() {
          try {
            while (true) {
              auto idx = next_idx.fetch_add(1);
              if (idx >= num) break;
              if (has_error) break;
              func(std::forward<Sequence1>(sequence1)[idx],
                  std::forward<Sequence2>(sequence2)[idx]);
            }
          } catch (...) {
            has_error = true;
            throw;
          }
        }));
  }
  for (auto& f : futures) f.get();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH BUILD
// -----------------------------------------------------------------------------
namespace yocto {

shape_bvh make_shape_bvh(const shape_data& shape, bool highquality) {
  if (is_points(shape)) {
    return {make_bvh(shape.points, highquality, [&](const int& point) {
      return point_bounds(shape.positions, shape.radius, point);
    })};
  } else if (is_lines(shape)) {
    return {make_bvh(shape.lines, highquality, [&](const vec2i& line) {
      return line_bounds(shape.positions, shape.radius, line);
    })};
  } else if (is_triangles(shape)) {
    return {make_bvh(shape.triangles, highquality, [&](const vec3i& triangle) {
      return triangle_bounds(shape.positions, triangle);
    })};
  } else if (is_quads(shape)) {
    return {make_bvh(shape.quads, highquality,
        [&](const vec4i& quad) { return quad_bounds(shape.positions, quad); })};
  } else {
    return {};
  }
}

scene_bvh make_scene_bvh(
    const scene_data& scene, bool highquality, bool noparallel) {
  // bvh
  auto bvh = scene_bvh{};

  // build shape bvh
  bvh.shapes = vector<shape_bvh>(scene.shapes.size());
  if (noparallel) {
    for (auto&& [sbvh, shape] : zip(bvh.shapes, scene.shapes)) {
      sbvh = make_shape_bvh(shape, highquality);
    }
  } else {
    parallel_zip(bvh.shapes, scene.shapes, [&](auto&& sbvh, auto&& shape) {
      sbvh = make_shape_bvh(shape, highquality);
    });
  }

  // instance bvh
  bvh.bvh = make_bvh(
      scene.instances, highquality, [&](const instance_data& instance) {
        return bvh.shapes[instance.shape].bvh.nodes.empty()
                   ? invalidb3f
                   : transform_bbox(instance.frame,
                         bvh.shapes[instance.shape].bvh.nodes[0].bbox);
      });

  // done
  return bvh;
}

void update_shape_bvh(shape_bvh& sbvh, const shape_data& shape) {
  // build primitives
  if (is_points(shape)) {
    refit_bvh(sbvh.bvh, shape.points, [&](const int& point) {
      return point_bounds(shape.positions, shape.radius, point);
    });
  } else if (is_lines(shape)) {
    refit_bvh(sbvh.bvh, shape.lines, [&](const vec2i& line) {
      return line_bounds(shape.positions, shape.radius, line);
    });
  } else if (is_triangles(shape)) {
    refit_bvh(sbvh.bvh, shape.triangles, [&](const vec3i& triangle) {
      return triangle_bounds(shape.positions, triangle);
    });
  } else if (is_quads(shape)) {
    refit_bvh(sbvh.bvh, shape.quads,
        [&](const vec4i& quad) { return quad_bounds(shape.positions, quad); });
  }
}

void update_scene_bvh(scene_bvh& sbvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes) {
  // update shapes
  for (auto shape : updated_shapes) {
    update_shape_bvh(sbvh.shapes[shape], scene.shapes[shape]);
  }

  // update nodes
  refit_bvh(sbvh.bvh, scene.instances, [&](const instance_data& instance) {
    return transform_bbox(
        instance.frame, sbvh.shapes[instance.shape].bvh.nodes[0].bbox);
  });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

intersection3f intersect_shape_bvh(const shape_bvh& sbvh,
    const shape_data& shape, const ray3f& ray, bool find_any) {
  if (is_points(shape)) {
    return intersect_elements_bvh(sbvh.bvh, shape.points, ray, find_any,
        [&shape](const ray3f& ray, const int& point) {
          return intersect_point(ray, shape.positions, shape.radius, point);
        });
  } else if (is_lines(shape)) {
    return intersect_elements_bvh(sbvh.bvh, shape.lines, ray, find_any,
        [&shape](const ray3f& ray, const vec2i& line) {
          return intersect_line(ray, shape.positions, shape.radius, line);
        });
  } else if (is_triangles(shape)) {
    return intersect_elements_bvh(sbvh.bvh, shape.triangles, ray, find_any,
        [&shape](const ray3f& ray, const vec3i& triangle) {
          return intersect_triangle(ray, shape.positions, triangle);
        });
  } else if (is_quads(shape)) {
    return intersect_elements_bvh(sbvh.bvh, shape.quads, ray, find_any,
        [&shape](const ray3f& ray, const vec4i& quad) {
          return intersect_quad(ray, shape.positions, quad);
        });
  } else {
    return {};
  }
}

intersection3f intersect_scene_bvh(const scene_bvh& sbvh,
    const scene_data& scene, const ray3f& ray, bool find_any) {
  return intersect_instances_bvh(sbvh.bvh, scene.instances, ray, find_any,
      [&](const ray3f& ray, const instance_data& instance, bool find_any) {
        auto inv_ray = transform_ray(inverse(instance.frame, true), ray);
        return intersect_shape_bvh(sbvh.shapes[instance.shape],
            scene.shapes[instance.shape], inv_ray, find_any);
      });
}

intersection3f intersect_instance_bvh(const scene_bvh& sbvh,
    const scene_data& scene, int instance_, const ray3f& ray, bool find_any) {
  auto& instance      = scene.instances[instance_];
  auto  inv_ray       = transform_ray(inverse(instance.frame, true), ray);
  auto  sintersection = intersect_shape_bvh(sbvh.shapes[instance.shape],
       scene.shapes[instance.shape], inv_ray, find_any);
  if (!sintersection.hit) return {};
  return {instance_, sintersection};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH OVERLAP
// -----------------------------------------------------------------------------
namespace yocto {

// Intersect ray with a bvh.
template <typename T, typename Func>
intersection3f overlap_elements_bvh(const bvh_data& bvh,
    const vector<T>& elements, const vec3f& pos, float max_distance,
    bool find_any, Func&& overlap_element) {
  // check if empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 64>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // intersection
  auto intersection = intersection3f{};

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx : range(node.start, node.start + node.num)) {
        auto eintersection = overlap_element(
            pos, max_distance, elements[bvh.primitives[idx]]);
        if (!eintersection.hit) continue;
        intersection = {bvh.primitives[idx], eintersection};
        max_distance = eintersection.distance;
      }
    }

    // check for early exit
    if (find_any && intersection.hit) return intersection;
  }

  return intersection;
}

// Intersect ray with a bvh.
intersection3f overlap_shape_bvh(const shape_bvh& sbvh, const shape_data& shape,
    const vec3f& pos, float max_distance, bool find_any) {
  if (is_points(shape)) {
    return overlap_elements_bvh(sbvh.bvh, shape.points, pos, max_distance,
        find_any, [&](const vec3f& pos, float max_distance, const int& point) {
          return overlap_point(
              pos, max_distance, shape.positions, shape.radius, point);
        });
  } else if (is_lines(shape)) {
    return overlap_elements_bvh(sbvh.bvh, shape.lines, pos, max_distance,
        find_any, [&](const vec3f& pos, float max_distance, const vec2i& line) {
          return overlap_line(
              pos, max_distance, shape.positions, shape.radius, line);
        });
  } else if (is_triangles(shape)) {
    return overlap_elements_bvh(sbvh.bvh, shape.triangles, pos, max_distance,
        find_any,
        [&](const vec3f& pos, float max_distance, const vec3i& triangle) {
          return overlap_triangle(
              pos, max_distance, shape.positions, shape.radius, triangle);
        });
  } else if (is_quads(shape)) {
    return overlap_elements_bvh(sbvh.bvh, shape.quads, pos, max_distance,
        find_any, [&](const vec3f& pos, float max_distance, const vec4i& quad) {
          return overlap_quad(
              pos, max_distance, shape.positions, shape.radius, quad);
        });
  } else {
    return {};
  }
}

// Intersect ray with a bvh.
template <typename T, typename Func>
intersection3f overlap_instances_bvh(const bvh_data& bvh,
    const vector<T>& elements, const vec3f& pos, float max_distance,
    bool find_any, Func&& overlap_element) {
  // check if empty
  if (bvh.nodes.empty()) return {};

  // node stack
  auto node_stack        = array<int, 64>{};
  auto node_cur          = 0;
  node_stack[node_cur++] = 0;

  // intersection
  auto intersection = intersection3f{};

  // walking stack
  while (node_cur != 0) {
    // grab node
    auto& node = bvh.nodes[node_stack[--node_cur]];

    // intersect bbox
    if (!overlap_bbox(pos, max_distance, node.bbox)) continue;

    // intersect node, switching based on node type
    // for each type, iterate over the the primitive list
    if (node.internal) {
      // internal node
      node_stack[node_cur++] = node.start + 0;
      node_stack[node_cur++] = node.start + 1;
    } else {
      for (auto idx : range(node.start, node.start + node.num)) {
        auto sintersection = overlap_element(
            pos, max_distance, elements[bvh.primitives[idx]], find_any);
        if (!sintersection.hit) continue;
        intersection = {bvh.primitives[idx], sintersection};
        max_distance = sintersection.distance;
      }
    }

    // check for early exit
    if (find_any && intersection.hit) return intersection;
  }

  return intersection;
}

// Intersect ray with a bvh.
intersection3f overlap_scene_bvh(const scene_bvh& sbvh, const scene_data& scene,
    const vec3f& pos, float max_distance, bool find_any) {
  return overlap_instances_bvh(sbvh.bvh, scene.instances, pos, max_distance,
      find_any,
      [&](const vec3f& pos, float max_distance, const instance_data& instance,
          bool find_any) {
        auto inv_pos = transform_point(inverse(instance.frame, true), pos);
        return overlap_shape_bvh(sbvh.shapes[instance.shape],
            scene.shapes[instance.shape], inv_pos, max_distance, find_any);
      });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EMBREE WRAPPER
// -----------------------------------------------------------------------------
namespace yocto {

// Embree error
struct embree_error : std::runtime_error {
  embree_error(const string& msg) : std::runtime_error{msg} {}
};

#ifdef YOCTO_EMBREE

// Check if embree is supported
bool embree_supported() { return true; }

// Get Embree device
std::atomic<ssize_t> embree_memory = 0;
static RTCDevice     embree_device() {
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
          embree_memory += bytes;
          return true;
        },
        nullptr);
  }
  return device;
}

// Clear Embree bvh
void clear_ebvh(void* ebvh) {
  if (ebvh) rtcReleaseScene((RTCScene)ebvh);
}

// Build the bvh acceleration structure.
shape_ebvh make_shape_ebvh(const shape_data& shape, bool highquality) {
  auto sbvh    = shape_ebvh{};
  auto edevice = embree_device();
  sbvh.ebvh    = unique_ptr<void, void (*)(void*)>{
      rtcNewScene(edevice), &clear_ebvh};
  auto escene = (RTCScene)sbvh.ebvh.get();
  if (highquality) {
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  } else {
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  }
  if (is_points(shape)) {
    throw std::runtime_error("embree does not support points");
  } else if (is_lines(shape)) {
    auto elines     = vector<int>{};
    auto epositions = vector<vec4f>{};
    auto last_index = -1;
    for (auto& [v1, v2] : shape.lines) {
      if (last_index == v1) {
        elines.push_back((int)epositions.size() - 1);
        epositions.push_back({shape.positions[v2], shape.radius[v2]});
      } else {
        elines.push_back((int)epositions.size());
        epositions.push_back({shape.positions[v1], shape.radius[v1]});
        epositions.push_back({shape.positions[v2], shape.radius[v2]});
      }
      last_index = v2;
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
    rtcReleaseGeometry(egeometry);
  } else if (is_triangles(shape)) {
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
    rtcReleaseGeometry(egeometry);
  } else if (is_quads(shape)) {
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
    rtcReleaseGeometry(egeometry);
  } else {
    throw std::runtime_error("empty shapes not supported");
  }
  rtcCommitScene(escene);
  return sbvh;
}

scene_ebvh make_scene_ebvh(
    const scene_data& scene, bool highquality, bool noparallel) {
  // scene bvh
  auto bvh = scene_ebvh{};

  // shape bvhs
  bvh.shapes = vector<shape_ebvh>(scene.shapes.size());
  if (noparallel) {
    for (auto&& [sshape, shape] : zip(bvh.shapes, scene.shapes)) {
      sshape = make_shape_ebvh(shape, highquality);
    }
  } else {
    parallel_zip(bvh.shapes, scene.shapes, [&](auto&& sbvh, auto&& shape) {
      sbvh = make_shape_ebvh(shape, highquality);
    });
  }

  // scene bvh
  auto edevice = embree_device();
  bvh.ebvh     = unique_ptr<void, void (*)(void*)>{
      rtcNewScene(edevice), &clear_ebvh};
  auto escene = (RTCScene)bvh.ebvh.get();
  if (highquality) {
    rtcSetSceneBuildQuality(escene, RTC_BUILD_QUALITY_HIGH);
  } else {
    rtcSetSceneFlags(escene, RTC_SCENE_FLAG_COMPACT);
  }
  for (auto&& [instance_id, instance] : enumerate(scene.instances)) {
    auto egeometry = rtcNewGeometry(edevice, RTC_GEOMETRY_TYPE_INSTANCE);
    rtcSetGeometryInstancedScene(
        egeometry, (RTCScene)bvh.shapes[instance.shape].ebvh.get());
    rtcSetGeometryTransform(
        egeometry, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(egeometry);
    rtcAttachGeometryByID(escene, egeometry, (uint)instance_id);
    rtcReleaseGeometry(egeometry);
  }
  rtcCommitScene(escene);
  return bvh;
}

// Refit bvh data
void update_shape_ebvh(shape_ebvh& sbvh, const shape_data& shape) {
  throw embree_error{"feature not supported"};
}
void update_scene_ebvh(scene_ebvh& sbvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes) {
  // scene bvh
  auto escene = (RTCScene)sbvh.ebvh.get();
  for (auto instance_id : updated_instances) {
    auto& instance    = scene.instances[instance_id];
    auto  embree_geom = rtcGetGeometry(escene, instance_id);
    rtcSetGeometryInstancedScene(
        embree_geom, (RTCScene)sbvh.shapes[instance.shape].ebvh.get());
    rtcSetGeometryTransform(
        embree_geom, 0, RTC_FORMAT_FLOAT3X4_COLUMN_MAJOR, &instance.frame);
    rtcCommitGeometry(embree_geom);
  }
  rtcCommitScene(escene);
}

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
intersection3f intersect_shape_ebvh(const shape_ebvh& sbvh,
    const shape_data& shape, const ray3f& ray, bool find_any) {
  RTCRayHit embree_ray;
  embree_ray.ray.org_x     = ray.o[0];
  embree_ray.ray.org_y     = ray.o[1];
  embree_ray.ray.org_z     = ray.o[2];
  embree_ray.ray.dir_x     = ray.d[0];
  embree_ray.ray.dir_y     = ray.d[1];
  embree_ray.ray.dir_z     = ray.d[2];
  embree_ray.ray.tnear     = ray.tmin;
  embree_ray.ray.tfar      = ray.tmax;
  embree_ray.ray.flags     = 0;
  embree_ray.ray.mask      = uint_max;
  embree_ray.ray.id        = 0;
  embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
  embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  RTCIntersectContext embree_ctx;
  rtcInitIntersectContext(&embree_ctx);
  rtcIntersect1((RTCScene)sbvh.ebvh.get(), &embree_ctx, &embree_ray);
  if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return {};
  auto element  = (int)embree_ray.hit.primID;
  auto uv       = vec2f{embree_ray.hit.u, embree_ray.hit.v};
  auto distance = embree_ray.ray.tfar;
  return {element, uv, distance};
}

intersection3f intersect_scene_ebvh(const scene_ebvh& sbvh,
    const scene_data& scene, const ray3f& ray, bool find_any) {
  RTCRayHit embree_ray;
  embree_ray.ray.org_x     = ray.o[0];
  embree_ray.ray.org_y     = ray.o[1];
  embree_ray.ray.org_z     = ray.o[2];
  embree_ray.ray.dir_x     = ray.d[0];
  embree_ray.ray.dir_y     = ray.d[1];
  embree_ray.ray.dir_z     = ray.d[2];
  embree_ray.ray.tnear     = ray.tmin;
  embree_ray.ray.tfar      = ray.tmax;
  embree_ray.ray.flags     = 0;
  embree_ray.ray.mask      = uint_max;
  embree_ray.ray.id        = 0;
  embree_ray.hit.geomID    = RTC_INVALID_GEOMETRY_ID;
  embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
  RTCIntersectContext embree_ctx;
  rtcInitIntersectContext(&embree_ctx);
  rtcIntersect1((RTCScene)sbvh.ebvh.get(), &embree_ctx, &embree_ray);
  if (embree_ray.hit.geomID == RTC_INVALID_GEOMETRY_ID) return {};
  auto instance = (int)embree_ray.hit.instID[0];
  auto element  = (int)embree_ray.hit.primID;
  auto uv       = vec2f{embree_ray.hit.u, embree_ray.hit.v};
  auto distance = embree_ray.ray.tfar;
  return {instance, element, uv, distance};
}

intersection3f intersect_instance_ebvh(const scene_ebvh& sbvh,
    const scene_data& scene, int instance_, const ray3f& ray, bool find_any) {
  auto& instance     = scene.instances[instance_];
  auto  inv_ray      = transform_ray(inverse(instance.frame, true), ray);
  auto  intersection = intersect_shape_ebvh(sbvh.shapes[instance.shape],
       scene.shapes[instance.shape], inv_ray, find_any);
  if (!intersection.hit) return {};
  return {instance_, intersection};
}

#else

// Check if embree is supported
bool embree_supported() { return false; }

// Not implemented
shape_ebvh make_shape_ebvh(const shape_data& shape, bool highquality) {
  throw embree_error{"Embree not available"};
}
scene_ebvh make_scene_ebvh(
    const scene_data& scene, bool highquality, bool noparallel) {
  throw embree_error{"Embree not available"};
}

// Not implemented
void update_shape_ebvh(shape_ebvh& sbvh, const shape_data& shape) {
  throw embree_error{"Embree not available"};
}
void update_scene_ebvh(scene_ebvh& sbvh, const scene_data& scene,
    const vector<int>& updated_instances, const vector<int>& updated_shapes) {
  throw embree_error{"Embree not available"};
}

// Not implemented
intersection3f intersect_shape_ebvh(const shape_ebvh& sbvh,
    const shape_data& shape, const ray3f& ray, bool find_any) {
  throw embree_error{"Embree not available"};
}
intersection3f intersect_scene_ebvh(const scene_ebvh& sbvh,
    const scene_data& scene, const ray3f& ray, bool find_any) {
  throw embree_error{"Embree not available"};
}
intersection3f intersect_instance_ebvh(const scene_ebvh& sbvh,
    const scene_data& scene, int instance, const ray3f& ray, bool find_any) {
  throw embree_error{"Embree not available"};
}

#endif

}  // namespace yocto