//
// Implementation for Yocto/GL Input and Output functions.
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
// TODO: update transforms -> should be compute transforms?
// TODO: update tesselation -> should be compute tesselation
// TODO: move out animation utilities
//

#include "yocto_sceneio.h"

#include <atomic>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <deque>
#include <future>
#include <memory>

#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
static inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)

    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get directory name (including '/').
static inline string get_dirname(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
static inline string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
static inline string get_filename(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
static inline string get_noextension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
static inline string get_basename(const string& filename) {
  return get_noextension(get_filename(filename));
}

// Replaces extensions
static inline string replace_extension(
    const string& filename, const string& ext) {
  return get_noextension(filename) + ext;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, Func&& func) {
  auto             futures  = vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(begin);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, end]() {
          while (true) {
            auto idx = next_idx.fetch_add(1);
            if (idx >= end) break;
            func(idx);
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Find the first keyframe value that is greater than the argument.
inline int keyframe_index(const vector<float>& times, const float& time) {
  for (auto i = 0; i < times.size(); i++)
    if (times[i] > time) return i;
  return (int)times.size();
}

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T keyframe_step(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f keyframe_slerp(
    const vector<float>& times, const vector<vec4f>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline T keyframe_linear(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evaluates a keyframed value using Bezier interpolation.
template <typename T>
inline T keyframe_bezier(
    const vector<float>& times, const vector<T>& vals, float time) {
  if (time <= times.front()) return vals.front();
  if (time >= times.back()) return vals.back();
  time     = clamp(time, times.front(), times.back() - 0.001f);
  auto idx = keyframe_index(times, time);
  auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
  return interpolate_bezier(
      vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE STATS AND VALIDATION
// -----------------------------------------------------------------------------
namespace yocto {

vector<string> scene_stats(const sceneio_model& scene, bool verbose) {
  auto accumulate = [](const auto& values, const auto& func) -> size_t {
    auto sum = (size_t)0;
    for (auto& value : values) sum += func(value);
    return sum;
  };
  auto format = [](auto num) {
    auto str = std::to_string(num);
    while (str.size() < 13) str = " " + str;
    return str;
  };
  auto format3 = [](auto num) {
    auto str = std::to_string(num.x) + " " + std::to_string(num.y) + " " +
               std::to_string(num.z);
    while (str.size() < 13) str = " " + str;
    return str;
  };

  auto bbox = compute_bounds(scene);

  auto stats = vector<string>{};
  stats.push_back("cameras:      " + format(scene.cameras.size()));
  stats.push_back("shapes:       " + format(scene.shapes.size()));
  stats.push_back("subdivs:      " + format(scene.subdivs.size()));
  stats.push_back("environments: " + format(scene.environments.size()));
  stats.push_back("textures:     " + format(scene.textures.size()));
  stats.push_back("nodes:        " + format(scene.nodes.size()));
  stats.push_back("animations:   " + format(scene.animations.size()));
  stats.push_back(
      "points:       " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.points.size(); })));
  stats.push_back(
      "lines:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.lines.size(); })));
  stats.push_back("triangles:    " +
                  format(accumulate(scene.shapes,
                      [](auto& shape) { return shape.triangles.size(); })));
  stats.push_back(
      "quads:        " + format(accumulate(scene.shapes,
                             [](auto& shape) { return shape.quads.size(); })));
  stats.push_back(
      "spoints:      " + format(accumulate(scene.subdivs,
                             [](auto& shape) { return shape.points.size(); })));
  stats.push_back(
      "slines:       " + format(accumulate(scene.subdivs,
                             [](auto& shape) { return shape.lines.size(); })));
  stats.push_back("striangles:   " +
                  format(accumulate(scene.subdivs,
                      [](auto& shape) { return shape.triangles.size(); })));
  stats.push_back(
      "squads:       " + format(accumulate(scene.subdivs,
                             [](auto& shape) { return shape.quads.size(); })));
  stats.push_back("sfvquads:     " +
                  format(accumulate(scene.subdivs,
                      [](auto& shape) { return shape.quadspos.size(); })));
  stats.push_back(
      "texels4b:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.ldr.size().x * (size_t)texture.ldr.size().x;
      })));
  stats.push_back(
      "texels4f:     " + format(accumulate(scene.textures, [](auto& texture) {
        return (size_t)texture.hdr.size().x * (size_t)texture.hdr.size().y;
      })));
  stats.push_back("center:       " + format3(center(bbox)));
  stats.push_back("size:         " + format3(size(bbox)));

  return stats;
}

// Checks for validity of the scene.
vector<string> scene_validation(const sceneio_model& scene, bool notextures) {
  auto errs        = vector<string>();
  auto check_names = [&errs](const auto& vals, const string& base) {
    auto used = unordered_map<string, int>();
    used.reserve(vals.size());
    for (auto& value : vals) used[value.name] += 1;
    for (auto& [name, used] : used) {
      if (name == "") {
        errs.push_back("empty " + base + " name");
      } else if (used > 1) {
        errs.push_back("duplicated " + base + " name " + name);
      }
    }
  };
  auto check_empty_textures = [&errs](const vector<sceneio_texture>& vals) {
    for (auto& value : vals) {
      if (value.hdr.empty() && value.ldr.empty()) {
        errs.push_back("empty texture " + value.name);
      }
    }
  };

  check_names(scene.cameras, "camera");
  check_names(scene.shapes, "shape");
  check_names(scene.textures, "texture");
  check_names(scene.environments, "environment");
  check_names(scene.nodes, "node");
  check_names(scene.animations, "animation");
  if (!notextures) check_empty_textures(scene.textures);

  return errs;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const sceneio_model& scene) {
  auto shape_bbox = vector<bbox3f>{};
  auto bbox       = invalidb3f;
  for (auto& shape : scene.shapes) {
    auto sbvh = invalidb3f;
    for (auto p : shape.positions)
      sbvh = merge(sbvh, transform_point(shape.frame, p));
    if (shape.instances.empty()) {
      bbox = merge(bbox, sbvh);
    } else {
      for (auto& frame : shape.instances) {
        bbox = merge(bbox, transform_bbox(frame * shape.frame, sbvh));
      }
    }
  }
  return bbox;
}

// Add missing cameras.
void add_cameras(sceneio_model& scene) {
  if (!scene.cameras.empty()) return;
  auto& camera = scene.cameras.emplace_back();
  camera.name  = "default";
  // TODO: error in camera.lens and camera.film
  camera.orthographic = false;
  camera.film         = 0.036;
  camera.aperture     = 0;
  camera.lens         = 0.050;
  auto bbox           = compute_bounds(scene);
  auto center         = (bbox.max + bbox.min) / 2;
  auto bbox_radius    = length(bbox.max - bbox.min) / 2;
  auto camera_dir     = camera.frame.o - center;
  if (camera_dir == zero3f) camera_dir = {0, 0, 1};
  auto camera_dist = bbox_radius / camera.film;
  auto from        = camera_dir * (camera_dist * 1) + center;
  auto to          = center;
  auto up          = vec3f{0, 1, 0};
  camera.frame     = lookat_frame(from, to, up);
  camera.focus     = length(from - to);
}

// Add missing radius.
void add_radius(sceneio_model& scene, float radius = 0.001f) {
  for (auto& shape : scene.shapes) {
    if (shape.points.empty() && shape.lines.empty()) continue;
    if (!shape.radius.empty()) continue;
    shape.radius.assign(shape.positions.size(), radius);
  }
}

// Add a sky environment
void add_sky(sceneio_model& scene, float sun_angle) {
  auto texture     = sceneio_texture{};
  texture.name     = "sky";
  texture.filename = "textures/sky.hdr";
  texture.hdr      = make_sunsky({1024, 512}, sun_angle);
  scene.textures.push_back(texture);
  auto environment         = sceneio_environment{};
  environment.name         = "sky";
  environment.emission     = {1, 1, 1};
  environment.emission_tex = (int)scene.textures.size() - 1;
  scene.environments.push_back(environment);
}

// Reduce memory usage
void trim_memory(sceneio_model& scene) {
  for (auto& shape : scene.shapes) {
    shape.points.shrink_to_fit();
    shape.lines.shrink_to_fit();
    shape.triangles.shrink_to_fit();
    shape.quads.shrink_to_fit();
    shape.positions.shrink_to_fit();
    shape.normals.shrink_to_fit();
    shape.texcoords.shrink_to_fit();
    shape.colors.shrink_to_fit();
    shape.radius.shrink_to_fit();
    shape.tangents.shrink_to_fit();
  }
  for (auto& subdiv : scene.subdivs) {
    subdiv.points.shrink_to_fit();
    subdiv.lines.shrink_to_fit();
    subdiv.triangles.shrink_to_fit();
    subdiv.quads.shrink_to_fit();
    subdiv.quadspos.shrink_to_fit();
    subdiv.quadsnorm.shrink_to_fit();
    subdiv.quadstexcoord.shrink_to_fit();
    subdiv.positions.shrink_to_fit();
    subdiv.normals.shrink_to_fit();
    subdiv.texcoords.shrink_to_fit();
    subdiv.colors.shrink_to_fit();
    subdiv.radius.shrink_to_fit();
    subdiv.tangents.shrink_to_fit();
  }
  for (auto& texture : scene.textures) {
    texture.ldr.shrink_to_fit();
    texture.hdr.shrink_to_fit();
  }
  scene.cameras.shrink_to_fit();
  scene.shapes.shrink_to_fit();
  scene.textures.shrink_to_fit();
  scene.environments.shrink_to_fit();
  scene.nodes.shrink_to_fit();
  scene.animations.shrink_to_fit();
}

// Apply subdivision and displacement rules.
sceneio_subdiv subdivide_subdiv(const sceneio_subdiv& shape) {
  using std::ignore;
  if (!shape.subdivisions) return shape;
  auto tesselated         = shape;
  tesselated.subdivisions = 0;
  if (!shape.points.empty()) {
    throw std::runtime_error("point subdivision not supported");
  } else if (!shape.lines.empty()) {
    tie(ignore, tesselated.normals) = subdivide_lines(
        tesselated.lines, tesselated.normals, shape.subdivisions);
    tie(ignore, tesselated.texcoords) = subdivide_lines(
        tesselated.lines, tesselated.texcoords, shape.subdivisions);
    tie(ignore, tesselated.colors) = subdivide_lines(
        tesselated.lines, tesselated.colors, shape.subdivisions);
    tie(ignore, tesselated.radius) = subdivide_lines(
        tesselated.lines, tesselated.radius, shape.subdivisions);
    tie(tesselated.lines, tesselated.positions) = subdivide_lines(
        tesselated.lines, tesselated.positions, shape.subdivisions);
    if (shape.smooth)
      tesselated.normals = compute_tangents(
          tesselated.lines, tesselated.positions);
  } else if (!shape.triangles.empty()) {
    tie(ignore, tesselated.normals) = subdivide_triangles(
        tesselated.triangles, tesselated.normals, shape.subdivisions);
    tie(ignore, tesselated.texcoords) = subdivide_triangles(
        tesselated.triangles, tesselated.texcoords, shape.subdivisions);
    tie(ignore, tesselated.colors) = subdivide_triangles(
        tesselated.triangles, tesselated.colors, shape.subdivisions);
    tie(ignore, tesselated.radius) = subdivide_triangles(
        tesselated.triangles, tesselated.radius, shape.subdivisions);
    tie(tesselated.triangles, tesselated.positions) = subdivide_triangles(
        tesselated.triangles, tesselated.positions, shape.subdivisions);
    if (shape.smooth)
      tesselated.normals = compute_normals(
          tesselated.triangles, tesselated.positions);
  } else if (!shape.quads.empty() && !shape.catmullclark) {
    tie(ignore, tesselated.normals) = subdivide_quads(
        tesselated.quads, tesselated.normals, shape.subdivisions);
    tie(ignore, tesselated.texcoords) = subdivide_quads(
        tesselated.quads, tesselated.texcoords, shape.subdivisions);
    tie(ignore, tesselated.colors) = subdivide_quads(
        tesselated.quads, tesselated.colors, shape.subdivisions);
    tie(ignore, tesselated.radius) = subdivide_quads(
        tesselated.quads, tesselated.radius, shape.subdivisions);
    tie(tesselated.quads, tesselated.positions) = subdivide_quads(
        tesselated.quads, tesselated.positions, shape.subdivisions);
    if (tesselated.smooth)
      tesselated.normals = compute_normals(
          tesselated.quads, tesselated.positions);
  } else if (!shape.quads.empty() && shape.catmullclark) {
    tie(ignore, tesselated.normals) = subdivide_catmullclark(
        tesselated.quads, tesselated.normals, shape.subdivisions);
    tie(ignore, tesselated.texcoords) = subdivide_catmullclark(
        tesselated.quads, tesselated.texcoords, shape.subdivisions);
    tie(ignore, tesselated.colors) = subdivide_catmullclark(
        tesselated.quads, tesselated.colors, shape.subdivisions);
    tie(ignore, tesselated.radius) = subdivide_catmullclark(
        tesselated.quads, tesselated.radius, shape.subdivisions);
    tie(tesselated.quads, tesselated.positions) = subdivide_catmullclark(
        tesselated.quads, tesselated.positions, shape.subdivisions);
    if (tesselated.smooth)
      tesselated.normals = compute_normals(
          tesselated.quads, tesselated.positions);
  } else if (!shape.quadspos.empty() && !shape.catmullclark) {
    std::tie(tesselated.quadsnorm, tesselated.normals) = subdivide_quads(
        tesselated.quadsnorm, tesselated.normals, shape.subdivisions);
    std::tie(tesselated.quadstexcoord, tesselated.texcoords) = subdivide_quads(
        tesselated.quadstexcoord, tesselated.texcoords, shape.subdivisions);
    std::tie(tesselated.quadspos, tesselated.positions) = subdivide_quads(
        tesselated.quadspos, tesselated.positions, shape.subdivisions);
    if (tesselated.smooth) {
      tesselated.normals = compute_normals(
          tesselated.quadspos, tesselated.positions);
      tesselated.quadsnorm = tesselated.quadspos;
    }
  } else if (!shape.quadspos.empty() && shape.catmullclark) {
    std::tie(tesselated.quadstexcoord, tesselated.texcoords) =
        subdivide_catmullclark(tesselated.quadstexcoord, tesselated.texcoords,
            shape.subdivisions, true);
    std::tie(tesselated.quadsnorm, tesselated.normals) = subdivide_catmullclark(
        tesselated.quadsnorm, tesselated.normals, shape.subdivisions, true);
    std::tie(tesselated.quadspos, tesselated.positions) =
        subdivide_catmullclark(
            tesselated.quadspos, tesselated.positions, shape.subdivisions);
    if (shape.smooth) {
      tesselated.normals = compute_normals(
          tesselated.quadspos, tesselated.positions);
      tesselated.quadsnorm = tesselated.quadspos;
    } else {
      tesselated.normals   = {};
      tesselated.quadsnorm = {};
    }
  } else {
    throw std::runtime_error("empty shape");
  }
  return tesselated;
}
// Apply displacement to a shape
sceneio_subdiv displace_subdiv(
    const sceneio_model& scene, const sceneio_subdiv& subdiv) {
  // Evaluate a texture
  auto eval_texture = [](const sceneio_texture& texture,
                          const vec2f&          texcoord) -> vec4f {
    if (!texture.hdr.empty()) {
      return eval_image(texture.hdr, texcoord, false, false);
    } else if (!texture.ldr.empty()) {
      return eval_image(texture.ldr, texcoord, true, false, false);
    } else {
      return {1, 1, 1, 1};
    }
  };

  if (!subdiv.displacement || subdiv.displacement_tex < 0) return subdiv;
  auto& displacement = scene.textures[subdiv.displacement_tex];
  if (subdiv.texcoords.empty()) {
    throw std::runtime_error("missing texture coordinates");
    return {};
  }

  auto displaced             = subdiv;
  displaced.displacement     = 0;
  displaced.displacement_tex = -1;

  // simple case
  if (!subdiv.triangles.empty()) {
    auto normals = subdiv.normals;
    if (normals.empty())
      normals = compute_normals(subdiv.triangles, subdiv.positions);
    for (auto vid = 0; vid < subdiv.positions.size(); vid++) {
      auto disp = mean(xyz(eval_texture(displacement, subdiv.texcoords[vid])));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      displaced.positions[vid] += normals[vid] * subdiv.displacement * disp;
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      displaced.normals = compute_normals(
          displaced.triangles, displaced.positions);
    }
  } else if (!subdiv.quads.empty()) {
    auto normals = subdiv.normals;
    if (normals.empty())
      normals = compute_normals(subdiv.triangles, subdiv.positions);
    for (auto vid = 0; vid < subdiv.positions.size(); vid++) {
      auto disp = mean(xyz(eval_texture(displacement, subdiv.texcoords[vid])));
      if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
      displaced.positions[vid] += normals[vid] * subdiv.displacement * disp;
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      displaced.normals = compute_normals(displaced.quads, displaced.positions);
    }
  } else if (!subdiv.quadspos.empty()) {
    // facevarying case
    auto offset = vector<float>(subdiv.positions.size(), 0);
    auto count  = vector<int>(subdiv.positions.size(), 0);
    for (auto fid = 0; fid < subdiv.quadspos.size(); fid++) {
      auto qpos = subdiv.quadspos[fid];
      auto qtxt = subdiv.quadstexcoord[fid];
      for (auto i = 0; i < 4; i++) {
        auto disp = mean(
            xyz(eval_texture(displacement, subdiv.texcoords[qtxt[i]])));
        if (!is_hdr_filename(displacement.filename)) disp -= 0.5f;
        offset[qpos[i]] += subdiv.displacement * disp;
        count[qpos[i]] += 1;
      }
    }
    auto normals = compute_normals(subdiv.quadspos, subdiv.positions);
    for (auto vid = 0; vid < subdiv.positions.size(); vid++) {
      displaced.positions[vid] += normals[vid] * offset[vid] / count[vid];
    }
    if (subdiv.smooth || !subdiv.normals.empty()) {
      displaced.quadsnorm = subdiv.quadspos;
      displaced.normals   = compute_normals(
          displaced.quadspos, displaced.positions);
    }
  }
  return displaced;
}

void tesselate_subdiv(
    sceneio_model& scene, const sceneio_subdiv& subdiv, bool no_quads) {
  auto tesselated = subdiv;
  if (tesselated.subdivisions) tesselated = subdivide_subdiv(tesselated);
  if (tesselated.displacement) tesselated = displace_subdiv(scene, tesselated);
  if (!subdiv.quadspos.empty()) {
    std::tie(tesselated.quads, tesselated.positions, tesselated.normals,
        tesselated.texcoords) = split_facevarying(tesselated.quadspos,
        tesselated.quadsnorm, tesselated.quadstexcoord, tesselated.positions,
        tesselated.normals, tesselated.texcoords);
  }
  if (!tesselated.quads.empty() && no_quads) {
    tesselated.triangles = quads_to_triangles(tesselated.quads);
    tesselated.quads     = {};
  }
  auto& shape     = scene.shapes[tesselated.shape];
  shape.points    = tesselated.points;
  shape.lines     = tesselated.lines;
  shape.triangles = tesselated.triangles;
  shape.quads     = tesselated.quads;
  shape.positions = tesselated.positions;
  shape.normals   = tesselated.normals;
  shape.texcoords = tesselated.texcoords;
  shape.colors    = tesselated.colors;
  shape.radius    = tesselated.radius;
}

// Update animation transforms
void update_transforms(sceneio_model& scene, sceneio_animation& animation,
    float time, const string& anim_group) {
  if (anim_group != "" && anim_group != animation.group) return;

  if (!animation.translations.empty()) {
    auto value = vec3f{0, 0, 0};
    switch (animation.interpolation) {
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.translations, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.translations, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.translations, time);
        break;
      default: throw std::runtime_error("should not have been here");
    }
    for (auto target : animation.targets)
      scene.nodes[target].translation = value;
  }
  if (!animation.rotations.empty()) {
    auto value = vec4f{0, 0, 0, 1};
    switch (animation.interpolation) {
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.rotations, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.rotations, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.rotations, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].rotation = value;
  }
  if (!animation.scales.empty()) {
    auto value = vec3f{1, 1, 1};
    switch (animation.interpolation) {
      case sceneio_animation::interpolation_type::step:
        value = keyframe_step(animation.times, animation.scales, time);
        break;
      case sceneio_animation::interpolation_type::linear:
        value = keyframe_linear(animation.times, animation.scales, time);
        break;
      case sceneio_animation::interpolation_type::bezier:
        value = keyframe_bezier(animation.times, animation.scales, time);
        break;
    }
    for (auto target : animation.targets) scene.nodes[target].scale = value;
  }
}

// Update node transforms
void update_transforms(sceneio_model& scene, sceneio_node& node,
    const frame3f& parent = identity3x4f) {
  auto frame = parent * node.local * translation_frame(node.translation) *
               rotation_frame(node.rotation) * scaling_frame(node.scale);
  if (node.shape >= 0) scene.shapes[node.shape].frame = frame;
  if (node.instance >= 0)
    scene.shapes[node.shape].instances[node.instance] = frame;
  if (node.camera >= 0) scene.cameras[node.camera].frame = frame;
  if (node.environment >= 0) scene.environments[node.environment].frame = frame;
  for (auto child : node.children)
    update_transforms(scene, scene.nodes[child], frame);
}

// Update node transforms
void update_transforms(
    sceneio_model& scene, float time, const string& anim_group) {
  for (auto& agr : scene.animations)
    update_transforms(scene, agr, time, anim_group);
  for (auto& node : scene.nodes) node.children.clear();
  for (auto node_id = 0; node_id < scene.nodes.size(); node_id++) {
    auto& node = scene.nodes[node_id];
    if (node.parent >= 0) scene.nodes[node.parent].children.push_back(node_id);
  }
  for (auto& node : scene.nodes)
    if (node.parent < 0) update_transforms(scene, node);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Helpers for throwing
static void throw_format_error(const string& filename) {
  throw std::runtime_error{filename + ": unknown format"};
}
static void throw_dependent_error(const string& filename, const string& err) {
  throw std::runtime_error{filename + ": error in resource (" + err + ")"};
}
static void throw_emptyshape_error(const string& filename, const string& name) {
  throw std::runtime_error{filename + ": empty shape " + name};
}
static void throw_missing_reference_error(
    const string& filename, const string& type, const string& name) {
  throw std::runtime_error{filename + ": missing " + type + " " + name};
}

// Load/save a scene in the builtin YAML format.
static void load_yaml_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static void save_yaml_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load/save a scene from/to OBJ.
static void load_obj_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static void save_obj_scene(const string& filename, const sceneio_model& scene,
    bool instances, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static void load_ply_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static void save_ply_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load/save a scene from/to glTF.
static void load_gltf_scene(
    const string& filename, sceneio_model& scene, bool noparallel);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static void load_pbrt_scene(
    const string& filename, sceneio_model& scene, bool noparallel);
static void save_pbrt_scene(
    const string& filename, const sceneio_model& scene, bool noparallel);

// Load a scene
void load_scene(const string& filename, sceneio_model& scene, bool noparallel) {
  auto ext = get_extension(filename);
  if (ext == ".yaml" || ext == ".YAML") {
    return load_yaml_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, noparallel);
  } else {
    scene = {};
    throw_format_error(filename);
  }
}

// Save a scene
void save_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  auto ext = get_extension(filename);
  if (ext == ".yaml" || ext == ".YAML") {
    return save_yaml_scene(filename, scene, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, false, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, noparallel);
  } else {
    throw_format_error(filename);
  }
}

void load_texture(const string& filename, sceneio_texture& texture) {
  try {
    if (is_hdr_filename(texture.filename)) {
      load_image(get_dirname(filename) + texture.filename, texture.hdr);
    } else {
      load_imageb(get_dirname(filename) + texture.filename, texture.ldr);
    }
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

void save_texture(const string& filename, const sceneio_texture& texture) {
  try {
    if (!texture.hdr.empty()) {
      save_image(get_dirname(filename) + texture.filename, texture.hdr);
    } else {
      save_imageb(get_dirname(filename) + texture.filename, texture.ldr);
    }
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

void load_shape(const string& filename, sceneio_shape& shape) {
  try {
    load_shape(get_dirname(filename) + shape.filename, shape.points,
        shape.lines, shape.triangles, shape.quads, shape.positions,
        shape.normals, shape.texcoords, shape.colors, shape.radius);
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

void save_shape(const string& filename, const sceneio_shape& shape) {
  try {
    save_shape(get_dirname(filename) + shape.filename, shape.points,
        shape.lines, shape.triangles, shape.quads, shape.positions,
        shape.normals, shape.texcoords, shape.colors, shape.radius);
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

void load_subdiv(const string& filename, sceneio_subdiv& subdiv) {
  try {
    if (!subdiv.facevarying) {
      load_shape(get_dirname(filename) + subdiv.filename, subdiv.points,
          subdiv.lines, subdiv.triangles, subdiv.quads, subdiv.positions,
          subdiv.normals, subdiv.texcoords, subdiv.colors, subdiv.radius);
    } else {
      load_fvshape(get_dirname(filename) + subdiv.filename, subdiv.quadspos,
          subdiv.quadsnorm, subdiv.quadstexcoord, subdiv.positions,
          subdiv.normals, subdiv.texcoords);
    }
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

void save_subdiv(const string& filename, const sceneio_subdiv& subdiv) {
  try {
    if (subdiv.quadspos.empty()) {
      save_shape(get_dirname(filename) + subdiv.filename, subdiv.points,
          subdiv.lines, subdiv.triangles, subdiv.quads, subdiv.positions,
          subdiv.normals, subdiv.texcoords, subdiv.colors, subdiv.radius);
    } else {
      save_fvshape(get_dirname(filename) + subdiv.filename, subdiv.quadspos,
          subdiv.quadsnorm, subdiv.quadstexcoord, subdiv.positions,
          subdiv.normals, subdiv.texcoords);
    }
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }
}

// create and cleanup names and filenames
static string make_safe_name(
    const string& name_, const string& base, int count) {
  auto name = name_;
  if (name.empty()) name = base + std::to_string(count);
  if (name.front() == '-') name = "_" + name;
  if (name.front() >= '0' && name.front() <= '9') name = "_" + name;
  for (auto& c : name) {
    if (c == '-' || c == '_') continue;
    if (c >= '0' && c <= '9') continue;
    if (c >= 'a' && c <= 'z') continue;
    if (c >= 'A' && c <= 'Z') continue;
    c = '_';
  }
  std::transform(name.begin(), name.end(), name.begin(),
      [](unsigned char c) { return std::tolower(c); });
  return name;
}
static inline string make_safe_filename(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename) {
    if (c == ' ') c = '_';
  }
  return filename;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// YAML SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

// Save a scene in the builtin YAML format.
static void load_yaml_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  scene = {};

  // open file
  auto yaml = yaml_model{};
  load_yaml(filename, yaml);

  auto tmap = unordered_map<string, int>{{"", -1}};
  auto vmap = unordered_map<string, int>{{"", -1}};
  auto smap = unordered_map<string, int>{{"", -1}};

  // parse yaml reference
  auto get_yaml_ref = [](const yaml_element& yelment, const string& name,
                          int& value, const unordered_map<string, int>& refs) {
    auto ref = ""s;
    get_yaml_value(yelment, name, ref);
    if (ref == "") {
      value = -1;
    } else {
      if (refs.find(ref) == refs.end())
        throw std::invalid_argument{"missing reference to " + ref};
      value = refs.at(ref);
    }
  };

  // parse yaml reference
  auto get_yaml_texture = [&scene, &tmap](const yaml_element& yelment,
                              const string& name, int& value) {
    auto filename = ""s;
    get_yaml_value(yelment, name, filename);
    if (filename == "") {
      value = -1;
    } else {
      if (tmap.find(filename) == tmap.end()) {
        auto& texture    = scene.textures.emplace_back();
        texture.name     = get_basename(filename);
        texture.filename = filename;
        tmap[filename]   = (int)scene.textures.size() - 1;
      }
      value = tmap.at(filename);
    }
  };

  // parse yaml reference copy
  auto get_yaml_refcopy = [](const yaml_element& yelment, const string& name,
                              auto& value, const auto& values,
                              const unordered_map<string, int>& refs) {
    auto ref = ""s;
    get_yaml_value(yelment, name, ref);
    if (ref == "") return;
    if (refs.find(ref) == refs.end())
      throw std::invalid_argument{"missing reference to " + ref};
    value = values[refs.at(ref)];
  };

  // hacked groups for large models
  struct sceneio_group {
    string          filename = "";
    vector<frame3f> frames   = {};
  };
  auto groups  = vector<sceneio_group>{};
  auto igroups = vector<int>{};

  // hacked lists for backward compatibility
  auto materials = vector<sceneio_material>{};
  auto mmap      = unordered_map<string, int>{{"", -1}};

  // check for conversion errors
  try {
    // cameras
    for (auto& yelement : yaml.elements) {
      if (yelement.name == "cameras") {
        auto& camera = scene.cameras.emplace_back();
        get_yaml_value(yelement, "name", camera.name);
        get_yaml_value(yelement, "frame", camera.frame);
        get_yaml_value(yelement, "orthographic", camera.orthographic);
        get_yaml_value(yelement, "lens", camera.lens);
        get_yaml_value(yelement, "aspect", camera.aspect);
        get_yaml_value(yelement, "film", camera.film);
        get_yaml_value(yelement, "focus", camera.focus);
        get_yaml_value(yelement, "aperture", camera.aperture);
        if (has_yaml_value(yelement, "lookat")) {
          auto lookat = identity3x3f;
          get_yaml_value(yelement, "lookat", lookat);
          camera.frame = lookat_frame(lookat.x, lookat.y, lookat.z);
          camera.focus = length(lookat.x - lookat.y);
        }
      } else if (yelement.name == "environments") {
        auto& environment = scene.environments.emplace_back();
        get_yaml_value(yelement, "name", environment.name);
        get_yaml_value(yelement, "frame", environment.frame);
        get_yaml_value(yelement, "emission", environment.emission);
        get_yaml_texture(yelement, "emission_tex", environment.emission_tex);
        if (has_yaml_value(yelement, "uri")) {
          auto uri = ""s;
          get_yaml_value(yelement, "uri", uri);
          environment.name = get_basename(uri);
        }
        if (has_yaml_value(yelement, "lookat")) {
          auto lookat = identity3x3f;
          get_yaml_value(yelement, "lookat", lookat);
          environment.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
        }
      } else if (yelement.name == "materials") {
        auto& material = materials.emplace_back();
        auto  name     = ""s;
        get_yaml_value(yelement, "name", name);
        get_yaml_value(yelement, "emission", material.emission);
        get_yaml_value(yelement, "base", material.base);
        get_yaml_value(yelement, "metallic", material.metallic);
        get_yaml_value(yelement, "specular", material.specular);
        get_yaml_value(yelement, "roughness", material.roughness);
        get_yaml_value(yelement, "coat", material.coat);
        get_yaml_value(yelement, "transmission", material.transmission);
        get_yaml_value(yelement, "thin", material.thin);
        get_yaml_value(yelement, "ior", material.ior);
        get_yaml_value(yelement, "scattering", material.scattering);
        get_yaml_value(yelement, "phaseg", material.phaseg);
        get_yaml_value(yelement, "radius", material.radius);
        get_yaml_value(yelement, "opacity", material.opacity);
        get_yaml_value(yelement, "coat", material.coat);
        get_yaml_texture(yelement, "emission_tex", material.emission_tex);
        get_yaml_texture(yelement, "base_tex", material.base_tex);
        get_yaml_texture(yelement, "metallic_tex", material.metallic_tex);
        get_yaml_texture(yelement, "specular_tex", material.specular_tex);
        get_yaml_texture(
            yelement, "transmission_tex", material.transmission_tex);
        get_yaml_texture(yelement, "roughness_tex", material.roughness_tex);
        get_yaml_texture(yelement, "scattering_tex", material.scattering_tex);
        get_yaml_texture(yelement, "normal_tex", material.normal_tex);
        get_yaml_texture(yelement, "normal_tex", material.normal_tex);
        get_yaml_value(yelement, "gltf_textures", material.gltf_textures);
        if (has_yaml_value(yelement, "uri")) {
          get_yaml_value(yelement, "uri", name);
          mmap[name] = (int)materials.size() - 1;
          name       = get_basename(name);
        }
        mmap[name] = (int)materials.size() - 1;
      } else if (yelement.name == "shapes") {
        auto& shape = scene.shapes.emplace_back();
        get_yaml_value(yelement, "name", shape.name);
        get_yaml_value(yelement, "shape", shape.filename);
        if (has_yaml_value(yelement, "filename")) {
          get_yaml_value(yelement, "filename", shape.filename);
        }
        get_yaml_value(yelement, "frame", shape.frame);
        if (has_yaml_value(yelement, "lookat")) {
          auto lookat = identity3x3f;
          get_yaml_value(yelement, "lookat", lookat);
          shape.frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
        }
        if (has_yaml_value(yelement, "instances")) {
          auto& group = groups.emplace_back();
          get_yaml_value(yelement, "instances", group.filename);
          while (igroups.size() < scene.shapes.size())
            igroups.emplace_back() = -1;
          igroups.back() = (int)groups.size() - 1;
        }
        get_yaml_refcopy(yelement, "material", shape.material, materials, mmap);
        get_yaml_value(yelement, "emission", shape.material.emission);
        get_yaml_value(yelement, "base", shape.material.base);
        get_yaml_value(yelement, "metallic", shape.material.metallic);
        get_yaml_value(yelement, "specular", shape.material.specular);
        get_yaml_value(yelement, "roughness", shape.material.roughness);
        get_yaml_value(yelement, "coat", shape.material.coat);
        get_yaml_value(yelement, "transmission", shape.material.transmission);
        get_yaml_value(yelement, "thin", shape.material.thin);
        get_yaml_value(yelement, "ior", shape.material.ior);
        get_yaml_value(yelement, "scattering", shape.material.scattering);
        get_yaml_value(yelement, "phaseg", shape.material.phaseg);
        get_yaml_value(yelement, "radius", shape.material.radius);
        get_yaml_value(yelement, "opacity", shape.material.opacity);
        get_yaml_value(yelement, "coat", shape.material.coat);
        get_yaml_texture(yelement, "emission_tex", shape.material.emission_tex);
        get_yaml_texture(yelement, "base_tex", shape.material.base_tex);
        get_yaml_texture(yelement, "metallic_tex", shape.material.metallic_tex);
        get_yaml_texture(yelement, "specular_tex", shape.material.specular_tex);
        get_yaml_texture(
            yelement, "transmission_tex", shape.material.transmission_tex);
        get_yaml_texture(
            yelement, "roughness_tex", shape.material.roughness_tex);
        get_yaml_texture(
            yelement, "scattering_tex", shape.material.scattering_tex);
        get_yaml_texture(yelement, "normal_tex", shape.material.normal_tex);
        get_yaml_texture(yelement, "normal_tex", shape.material.normal_tex);
        get_yaml_value(yelement, "gltf_textures", shape.material.gltf_textures);
        if (has_yaml_value(yelement, "uri")) {
          get_yaml_value(yelement, "uri", shape.filename);
          shape.name           = get_basename(shape.filename);
          smap[shape.filename] = (int)scene.shapes.size() - 1;
        }
        smap[shape.name] = (int)scene.shapes.size() - 1;
      } else if (yelement.name == "subdivs") {
        auto& subdiv = scene.subdivs.emplace_back();
        get_yaml_value(yelement, "name", subdiv.name);
        get_yaml_value(yelement, "subdiv", subdiv.filename);
        if (has_yaml_value(yelement, "filename")) {
          get_yaml_value(yelement, "filename", subdiv.filename);
        }
        get_yaml_ref(yelement, "shape", subdiv.shape, smap);
        get_yaml_value(yelement, "subdivisions", subdiv.subdivisions);
        get_yaml_value(yelement, "catmullclark", subdiv.catmullclark);
        get_yaml_value(yelement, "smooth", subdiv.smooth);
        get_yaml_value(yelement, "facevarying", subdiv.facevarying);
        get_yaml_ref(
            yelement, "displacement_tex", subdiv.displacement_tex, tmap);
        get_yaml_value(yelement, "displacement", subdiv.displacement);
        if (has_yaml_value(yelement, "uri")) {
          get_yaml_value(yelement, "uri", subdiv.filename);
          subdiv.name = get_basename(subdiv.filename);
        }
      }
    }

  } catch (std::invalid_argument& e) {
    throw std::runtime_error{filename + ": parse error [" + e.what() + "]"};
  }

  // instance groups
  if (!groups.empty()) {
    // load groups
    if (noparallel) {
      for (auto& group : groups) {
        auto ply = ply_model{};
        try {
          load_ply(get_dirname(filename) + group.filename, ply);
        } catch (std::exception& e) {
          throw_dependent_error(filename, e.what());
        }
        group.frames = get_ply_values(ply, "frame",
            array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy",
                "zz", "ox", "oy", "oz"});
      }
    } else {
      parallel_foreach(groups, [filename](sceneio_group& group) {
        auto ply = ply_model{};
        try {
          load_ply(get_dirname(filename) + group.filename, ply);
        } catch (std::exception& e) {
          throw_dependent_error(filename, e.what());
        }
        group.frames = get_ply_values(ply, "frame",
            array<string, 12>{"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy",
                "zz", "ox", "oy", "oz"});
      });
    }
    for (auto idx = 0; idx < scene.shapes.size(); idx++) {
      scene.shapes[idx].instances = groups[igroups[idx]].frames;
    }
  }

  // load shapes
  for (auto& shape : scene.shapes) {
    if (!shape.positions.empty()) continue;
    load_shape(filename, shape);
  }

  // load subdivs
  for (auto& subdiv : scene.subdivs) {
    if (!subdiv.positions.empty()) continue;
    load_subdiv(filename, subdiv);
  }

  // load textures
  for (auto& texture : scene.textures) {
    if (!texture.hdr.empty() || !texture.ldr.empty()) continue;
    load_texture(filename, texture);
  }

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  trim_memory(scene);
}

// Save a scene in the builtin YAML format.
static void save_yaml_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // helper
  auto add_yaml_texture = [&scene](yaml_element& yelement, const string& name,
                              int ref) {
    if (ref < 0) return;
    add_yaml_value(yelement, name, scene.textures[ref].filename);
  };

  // save yaml file
  auto yaml = yaml_model{};

  for (auto stat : scene_stats(scene)) yaml.comments.push_back(stat);

  for (auto& camera : scene.cameras) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "cameras";
    add_yaml_value(yelement, "name", camera.name);
    add_yaml_value(yelement, "frame", camera.frame);
    if (camera.orthographic)
      add_yaml_value(yelement, "orthographic", camera.orthographic);
    add_yaml_value(yelement, "lens", camera.lens);
    add_yaml_value(yelement, "aspect", camera.aspect);
    add_yaml_value(yelement, "film", camera.film);
    add_yaml_value(yelement, "focus", camera.focus);
    add_yaml_value(yelement, "aperture", camera.aperture);
  }

  for (auto& environment : scene.environments) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "environments";
    add_yaml_value(yelement, "name", environment.name);
    add_yaml_value(yelement, "frame", environment.frame);
    add_yaml_value(yelement, "emission", environment.emission);
    add_yaml_texture(yelement, "emission_tex", environment.emission_tex);
  }

  for (auto& shape : scene.shapes) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "shapes";
    add_yaml_value(yelement, "name", shape.name);
    add_yaml_value(yelement, "filename", shape.filename);
    add_yaml_value(yelement, "frame", shape.frame);
    add_yaml_value(yelement, "emission", shape.material.emission);
    add_yaml_value(yelement, "base", shape.material.base);
    add_yaml_value(yelement, "specular", shape.material.specular);
    add_yaml_value(yelement, "metallic", shape.material.metallic);
    add_yaml_value(yelement, "coat", shape.material.coat);
    add_yaml_value(yelement, "transmission", shape.material.transmission);
    add_yaml_value(yelement, "roughness", shape.material.roughness);
    add_yaml_value(yelement, "thin", shape.material.thin);
    add_yaml_value(yelement, "ior", shape.material.ior);
    add_yaml_value(yelement, "scattering", shape.material.scattering);
    add_yaml_value(yelement, "phaseg", shape.material.phaseg);
    add_yaml_value(yelement, "radius", shape.material.radius);
    add_yaml_value(yelement, "opacity", shape.material.opacity);
    add_yaml_texture(yelement, "emission_tex", shape.material.emission_tex);
    add_yaml_texture(yelement, "base_tex", shape.material.base_tex);
    add_yaml_texture(yelement, "metallic_tex", shape.material.metallic_tex);
    add_yaml_texture(yelement, "specular_tex", shape.material.specular_tex);
    add_yaml_texture(yelement, "roughness_tex", shape.material.roughness_tex);
    add_yaml_texture(
        yelement, "transmission_tex", shape.material.transmission_tex);
    add_yaml_texture(yelement, "scattering_tex", shape.material.scattering_tex);
    add_yaml_texture(yelement, "coat_tex", shape.material.coat_tex);
    add_yaml_texture(yelement, "opacity_tex", shape.material.opacity_tex);
    add_yaml_texture(yelement, "normal_tex", shape.material.normal_tex);
    if (shape.material.gltf_textures)
      add_yaml_value(yelement, "gltf_textures", shape.material.gltf_textures);
    if (!shape.instances.empty()) {
      // TODO: save instances
    }
  }

  for (auto& subdiv : scene.subdivs) {
    auto& yelement = yaml.elements.emplace_back();
    yelement.name  = "subdivs";
    add_yaml_value(yelement, "name", subdiv.name);
    add_yaml_value(yelement, "filename", subdiv.filename);
    add_yaml_value(yelement, "shape", scene.shapes[subdiv.shape].name);
    add_yaml_value(yelement, "subdivisions", subdiv.subdivisions);
    add_yaml_value(yelement, "catmullclark", subdiv.catmullclark);
    add_yaml_value(yelement, "smooth", subdiv.smooth);
    add_yaml_value(yelement, "facevarying", subdiv.facevarying);
    add_yaml_texture(yelement, "displacement_tex", subdiv.displacement_tex);
    add_yaml_value(yelement, "displacement", subdiv.displacement);
  }

  // save yaml
  save_yaml(filename, yaml);

  // save shapes
  for (auto& shape : scene.shapes) {
    save_shape(filename, shape);
  }

  // save subdivs
  for (auto& subdiv : scene.subdivs) {
    save_subdiv(filename, subdiv);
  }

  // save textures
  for (auto& texture : scene.textures) {
    save_texture(filename, texture);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_obj(const string& filename, sceneio_model& scene) {
  // load obj
  auto obj = obj_model{};
  load_obj(filename, obj, false, true, true);

  // convert cameras
  for (auto& ocam : obj.cameras) {
    auto& camera = scene.cameras.emplace_back();
    camera.name  = make_safe_name(ocam.name, "cam", (int)scene.cameras.size());
    camera.frame = ocam.frame;
    camera.orthographic = ocam.ortho;
    camera.film         = max(ocam.width, ocam.height);
    camera.aspect       = ocam.width / ocam.height;
    camera.focus        = ocam.focus;
    camera.lens         = ocam.lens;
    camera.aperture     = ocam.aperture;
  }

  // helper to create texture maps
  auto texture_map = unordered_map<string, int>{{"", -1}};
  auto get_texture = [&texture_map, &scene](const obj_texture_info& info) {
    if (info.path == "") return -1;
    auto it = texture_map.find(info.path);
    if (it != texture_map.end()) return it->second;
    auto& texture = scene.textures.emplace_back();
    texture.name  = make_safe_name(
        get_basename(info.path), "texture", (int)scene.textures.size());
    texture.filename       = info.path;
    texture_map[info.path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert materials and textures
  auto materials    = vector<sceneio_material>{};
  auto material_map = unordered_map<string, int>{{"", -1}};
  for (auto& omat : obj.materials) {
    auto& material            = materials.emplace_back();
    material.emission         = omat.emission;
    material.base             = omat.diffuse;
    material.specular         = omat.specular != zero3f ? 1 : 0;
    material.roughness        = obj_exponent_to_roughness(omat.exponent);
    material.ior              = omat.ior;
    material.metallic         = omat.pbr_metallic;
    material.coat             = omat.pbr_clearcoat;
    material.transmission     = mean(omat.transmission);
    material.scattering       = omat.vol_scattering;
    material.phaseg           = omat.vol_anisotropy;
    material.radius           = omat.vol_scale;
    material.opacity          = omat.opacity;
    material.thin             = true;
    material.emission_tex     = get_texture(omat.emission_map);
    material.base_tex         = get_texture(omat.diffuse_map);
    material.specular_tex     = get_texture(omat.specular_map);
    material.metallic_tex     = get_texture(omat.pbr_metallic_map);
    material.roughness_tex    = get_texture(omat.pbr_roughness_map);
    material.transmission_tex = get_texture(omat.transmission_map);
    material.coat_tex         = get_texture(omat.pbr_clearcoat_map);
    material.opacity_tex      = get_texture(omat.opacity_map);
    material.normal_tex       = get_texture(omat.normal_map);
    material_map[omat.name]   = (int)materials.size() - 1;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto& oshape : obj.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shape.name  = oshape.name;
    if (shape.name == "") shape.name = "shape";
    shape_name_counts[shape.name] += 1;
    if (shape_name_counts[shape.name] > 1)
      shape.name += std::to_string(shape_name_counts[shape.name]);
    shape.name = make_safe_name(shape.name, "shape", (int)scene.shapes.size());
    shape.filename  = make_safe_filename("shapes/" + shape.name + ".ply");
    auto nmaterials = vector<string>{};
    auto ematerials = vector<int>{};
    auto has_quads  = has_obj_quads(oshape);
    if (!oshape.faces.empty() && !has_quads) {
      get_obj_triangles(obj, oshape, shape.triangles, shape.positions,
          shape.normals, shape.texcoords, nmaterials, ematerials, true);
    } else if (!oshape.faces.empty() && has_quads) {
      get_obj_quads(obj, oshape, shape.quads, shape.positions, shape.normals,
          shape.texcoords, nmaterials, ematerials, true);
    } else if (!oshape.lines.empty()) {
      get_obj_lines(obj, oshape, shape.lines, shape.positions, shape.normals,
          shape.texcoords, nmaterials, ematerials, true);
    } else if (!oshape.points.empty()) {
      get_obj_points(obj, oshape, shape.points, shape.positions, shape.normals,
          shape.texcoords, nmaterials, ematerials, true);
    } else {
      throw_emptyshape_error(filename, oshape.name);
    }
    // get material
    if (oshape.materials.size() != 1) {
      throw_missing_reference_error(filename, "material for", oshape.name);
    }
    if (material_map.find(oshape.materials.at(0)) == material_map.end()) {
      throw_missing_reference_error(
          filename, "material", oshape.materials.at(0));
    }
    shape.material  = materials[material_map.at(oshape.materials.at(0))];
    shape.instances = oshape.instances;
  }

  // convert environments
  for (auto& oenvironment : obj.environments) {
    auto& environment = scene.environments.emplace_back();
    environment.name  = make_safe_name(
        oenvironment.name, "environment", scene.environments.size());
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = get_texture(oenvironment.emission_map);
  }
}

// Loads an OBJ
static void load_obj_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // Parse obj
  load_obj(filename, scene);

  // load textures
  for (auto& texture : scene.textures) {
    if (!texture.hdr.empty() || !texture.ldr.empty()) continue;
    load_texture(filename, texture);
  }

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_radius(scene);
}

static void save_obj(
    const string& filename, const sceneio_model& scene, bool instances) {
  auto obj = obj_model{};

  for (auto stat : scene_stats(scene)) obj.comments.push_back(stat);

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = camera.name;
    ocamera.frame    = camera.frame;
    ocamera.ortho    = camera.orthographic;
    ocamera.width    = camera.film;
    ocamera.height   = camera.film / camera.aspect;
    ocamera.focus    = camera.focus;
    ocamera.lens     = camera.lens;
    ocamera.aperture = camera.aperture;
  }

  // textures
  auto get_texture = [&scene](int tex) {
    if (tex < 0) return obj_texture_info{};
    auto info = obj_texture_info{};
    info.path = scene.textures[tex].filename;
    return info;
  };

  // convert materials and textures
  for (auto& shape : scene.shapes) {
    auto& material     = shape.material;
    auto& omaterial    = obj.materials.emplace_back();
    omaterial.name     = shape.name;
    omaterial.illum    = 2;
    omaterial.emission = material.emission;
    omaterial.diffuse  = material.base;
    omaterial.specular = material.specular ? vec3f{1, 1, 1} : vec3f{0, 0, 0};
    omaterial.exponent = obj_roughness_to_exponent(material.roughness);
    omaterial.pbr_metallic      = material.metallic;
    omaterial.pbr_clearcoat     = material.coat;
    omaterial.transmission      = vec3f{material.transmission};
    omaterial.opacity           = material.opacity;
    omaterial.emission_map      = get_texture(material.emission_tex);
    omaterial.diffuse_map       = get_texture(material.base_tex);
    omaterial.specular_map      = get_texture(material.specular_tex);
    omaterial.pbr_metallic_map  = get_texture(material.metallic_tex);
    omaterial.pbr_roughness_map = get_texture(material.roughness_tex);
    omaterial.transmission_map  = get_texture(material.transmission_tex);
    omaterial.pbr_clearcoat_map = get_texture(material.coat_tex);
    omaterial.opacity_map       = get_texture(material.opacity_tex);
    omaterial.normal_map        = get_texture(material.normal_tex);
  }

  // convert shapes
  for (auto& shape : scene.shapes) {
    if (instances && !shape.instances.empty()) {
      auto positions = shape.positions, normals = shape.normals;
      for (auto& p : positions) p = transform_point(shape.frame, p);
      for (auto& n : normals) n = transform_normal(shape.frame, n);
      if (!shape.triangles.empty()) {
        add_obj_triangles(obj, shape.name, shape.triangles, positions, normals,
            shape.texcoords, {}, {}, true);
      } else if (!shape.quads.empty()) {
        add_obj_quads(obj, shape.name, shape.quads, positions, normals,
            shape.texcoords, {}, {}, true);
      } else if (!shape.lines.empty()) {
        add_obj_lines(obj, shape.name, shape.lines, positions, normals,
            shape.texcoords, {}, {}, true);
      } else if (!shape.points.empty()) {
        add_obj_points(obj, shape.name, shape.points, positions, normals,
            shape.texcoords, {}, {}, true);
      } else {
        throw_emptyshape_error(filename, shape.name);
      }
      // TODO: instances
    } else {
      for (auto& frame : shape.instances) {
        auto materials = vector<string>{shape.name};
        auto positions = shape.positions, normals = shape.normals;
        for (auto& p : positions) p = transform_point(frame * shape.frame, p);
        for (auto& n : normals) n = transform_normal(frame * shape.frame, n);
        if (!shape.triangles.empty()) {
          add_obj_triangles(obj, shape.name, shape.triangles, positions,
              normals, shape.texcoords, materials, {}, true);
        } else if (!shape.quads.empty()) {
          add_obj_quads(obj, shape.name, shape.quads, positions, normals,
              shape.texcoords, materials, {}, true);
        } else if (!shape.lines.empty()) {
          add_obj_lines(obj, shape.name, shape.lines, positions, normals,
              shape.texcoords, materials, {}, true);
        } else if (!shape.points.empty()) {
          add_obj_points(obj, shape.name, shape.points, positions, normals,
              shape.texcoords, materials, {}, true);
        } else {
          throw_emptyshape_error(filename, shape.name);
        }
      }
    }
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = environment.name;
    oenvironment.frame        = environment.frame;
    oenvironment.emission     = environment.emission;
    oenvironment.emission_map = get_texture(environment.emission_tex);
  }

  // save obj
  save_obj(filename, obj);
}

static void save_obj_scene(const string& filename, const sceneio_model& scene,
    bool instances, bool noparallel) {
  save_obj(filename, scene, instances);

  // save textures
  for (auto& texture : scene.textures) {
    save_texture(filename, texture);
  }
}

void print_obj_camera(const sceneio_camera& camera) {
  printf("c %s %d %g %g %g %g %g %g %g %g %g %g%g %g %g %g %g %g %g\n",
      camera.name.c_str(), (int)camera.orthographic, camera.film,
      camera.film / camera.aspect, camera.lens, camera.focus, camera.aperture,
      camera.frame.x.x, camera.frame.x.y, camera.frame.x.z, camera.frame.y.x,
      camera.frame.y.y, camera.frame.y.z, camera.frame.z.x, camera.frame.z.y,
      camera.frame.z.z, camera.frame.o.x, camera.frame.o.y, camera.frame.o.z);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void load_ply_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  scene = {};

  // load ply mesh
  scene.shapes.push_back({});
  auto& shape    = scene.shapes.back();
  shape.name     = "shape";
  shape.filename = get_filename(filename);
  try {
    load_shape(filename, shape.points, shape.lines, shape.triangles,
        shape.quads, shape.positions, shape.normals, shape.texcoords,
        shape.colors, shape.radius);
  } catch (std::exception& e) {
    throw_dependent_error(filename, e.what());
  }

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_radius(scene);
}

static void save_ply_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  if (scene.shapes.empty()) throw_emptyshape_error(filename, "");
  auto& shape = scene.shapes.front();
  save_shape(filename, shape.points, shape.lines, shape.triangles, shape.quads,
      shape.positions, shape.normals, shape.texcoords, shape.colors,
      shape.radius);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static void load_gltf_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // load gltf
  auto gltf = gltf_model{};
  load_gltf(filename, gltf);

  // convert textures
  for (auto& gtexture : gltf.textures) {
    auto& texture = scene.textures.emplace_back();
    if (!gtexture.name.empty()) {
      texture.name = make_safe_name(
          gtexture.name, "texture", (int)scene.textures.size());
    } else {
      texture.name = make_safe_name(get_basename(gtexture.filename), "texture",
          (int)scene.textures.size());
    }
    texture.filename = gtexture.filename;
  }

  // convert materials
  auto materials = vector<sceneio_material>{};
  for (auto& gmaterial : gltf.materials) {
    auto& material        = materials.emplace_back();
    material.emission     = gmaterial.emission;
    material.emission_tex = gmaterial.emission_tex;
    if (gmaterial.has_specgloss) {
      material.base    = xyz(gmaterial.sg_diffuse);
      material.opacity = gmaterial.sg_diffuse.w;
      // TODO: better specular
      // material.specular     = gmaterial.sg_specular;
      material.base_tex = gmaterial.sg_diffuse_tex;
      // material.specular_tex = gmaterial.sg_specular_tex;
    } else if (gmaterial.has_metalrough) {
      material.base         = xyz(gmaterial.mr_base);
      material.opacity      = gmaterial.mr_base.w;
      material.specular     = 1;
      material.base_tex     = gmaterial.mr_base_tex;
      material.metallic_tex = gmaterial.mr_metallic_tex;
    }
    material.normal_tex = gmaterial.normal_tex;
  }

  // convert shapes
  auto shape_indices = vector<vector<int>>{};
  for (auto& gmesh : gltf.meshes) {
    shape_indices.push_back({});
    for (auto& gprim : gmesh.primitives) {
      auto& shape    = scene.shapes.emplace_back();
      shape.material = materials[gprim.material];
      shape_indices.back().push_back((int)scene.shapes.size() - 1);
      shape.name =
          gmesh.name.empty()
              ? ""s
              : (gmesh.name + std::to_string(shape_indices.back().size()));
      make_safe_name(shape.name, "shape", (int)scene.shapes.size());
      shape.filename = make_safe_filename(
          "shapes/shape" + std::to_string(scene.shapes.size()));
      shape.positions = gprim.positions;
      shape.normals   = gprim.normals;
      shape.texcoords = gprim.texcoords;
      shape.colors    = gprim.colors;
      shape.radius    = gprim.radius;
      shape.tangents  = gprim.tangents;
      shape.triangles = gprim.triangles;
      shape.lines     = gprim.lines;
      shape.points    = gprim.points;
    }
  }

  // convert cameras
  auto cameras = vector<sceneio_camera>{};
  for (auto& gcamera : gltf.cameras) {
    auto& camera  = cameras.emplace_back();
    camera.name   = gcamera.name;
    camera.aspect = gcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = gcamera.aspect >= 1
                      ? (2 * camera.aspect * tan(gcamera.yfov / 2))
                      : (2 * tan(gcamera.yfov / 2));
    camera.focus = 10;
  }

  // convert scene nodes
  for (auto& gnode : gltf.nodes) {
    if (gnode.camera >= 0) {
      auto& camera = scene.cameras.emplace_back(cameras[gnode.camera]);
      camera.name  = make_safe_name(
          camera.name, "caemra", (int)scene.cameras.size());
      camera.frame = gnode.frame;
    }
    if (gnode.mesh >= 0) {
      for (auto shape_id : shape_indices[gnode.mesh]) {
        scene.shapes[shape_id].instances.push_back(gnode.frame);
      }
    }
  }

  // load textures
  for (auto& texture : scene.textures) {
    if (!texture.hdr.empty() || !texture.ldr.empty()) continue;
    load_texture(filename, texture);
  }

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_radius(scene);

  // fix cameras
  auto bbox = compute_bounds(scene);
  for (auto& camera : scene.cameras) {
    auto center   = (bbox.min + bbox.max) / 2;
    auto distance = dot(-camera.frame.z, center - camera.frame.o);
    if (distance > 0) camera.focus = distance;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static void load_pbrt_scene(
    const string& filename, sceneio_model& scene, bool noparallel) {
  // load pbrt
  auto pbrt = pbrt_model{};
  load_pbrt(filename, pbrt);

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.name   = make_safe_name("", "camera", (int)scene.cameras.size());
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = pcamera.lens;
    camera.focus  = pcamera.focus;
  }

  // convert textures
  auto texture_map = unordered_map<string, int>{{"", -1}};
  for (auto& ptexture : pbrt.textures) {
    if (ptexture.filename.empty()) continue;
    auto& texture = scene.textures.emplace_back();
    texture.name  = make_safe_name(
        ptexture.name, "texture", (int)scene.textures.size());
    texture.filename           = ptexture.filename;
    texture_map[ptexture.name] = (int)scene.textures.size() - 1;
  }

  // convert materials
  auto get_texture = [&texture_map](const string& name) {
    if (name == "") return -1;
    if (texture_map.find(name) == texture_map.end())
      throw std::runtime_error("cannot find texture " + name);
    return texture_map.at(name);
  };
  auto materials    = vector<sceneio_material>{};
  auto material_map = unordered_map<string, int>{{"", -1}};
  for (auto& pmaterial : pbrt.materials) {
    auto& material = materials.emplace_back();
    material.base  = pmaterial.diffuse;
    // TODO: better conversion
    // material.specular     = pmaterial.specular;
    // material.transmission = pmaterial.transmission;
    material.roughness = mean(pmaterial.roughness);
    material.opacity   = pmaterial.opacity == vec3f{1} ? 1
                                                     : mean(pmaterial.opacity);
    material.base_tex            = get_texture(pmaterial.diffuse_map);
    material_map[pmaterial.name] = (int)materials.size() - 1;
  }

  // convert arealights
  auto arealight_map = unordered_map<string, int>{{"", -1}};
  for (auto& parealight : pbrt.arealights) {
    auto& material                 = materials.emplace_back();
    material.emission              = parealight.emission;
    arealight_map[parealight.name] = (int)materials.size() - 1;
  }

  // convert shapes
  for (auto& pshape : pbrt.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shape.name  = make_safe_name(
        get_basename(shape.filename), "shape", (int)scene.shapes.size());
    if (pshape.filename.empty()) {
      shape.name     = make_safe_name("", "shape", (int)scene.shapes.size());
      shape.filename = make_safe_filename(
          "shapes/shape" + std::to_string(scene.shapes.size()) + ".ply");
    } else {
      shape.filename = pshape.filename;
      shape.name     = make_safe_name(
          get_basename(pshape.filename), "shape", (int)scene.shapes.size());
    }
    shape.frame     = pshape.frame;
    shape.positions = pshape.positions;
    shape.normals   = pshape.normals;
    shape.texcoords = pshape.texcoords;
    shape.triangles = pshape.triangles;
    for (auto& uv : shape.texcoords) uv.y = 1 - uv.y;
    auto material_id  = material_map.at(pshape.material);
    auto arealight_id = arealight_map.at(pshape.arealight);
    shape.material  = materials[arealight_id >= 0 ? arealight_id : material_id];
    shape.instances = pshape.instance_frames;
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment = scene.environments.emplace_back();
    environment.name  = make_safe_name(
        "", "environment", (int)scene.environments.size());
    environment.frame    = penvironment.frame;
    environment.emission = penvironment.emission;
    if (!penvironment.filename.empty()) {
      auto& texture    = scene.textures.emplace_back();
      texture.name     = make_safe_name(get_basename(penvironment.filename),
          "environment", (int)scene.environments.size());
      texture.filename = penvironment.filename;
      environment.emission_tex = (int)scene.textures.size() - 1;
    } else {
      environment.emission_tex = -1;
    }
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape     = scene.shapes.emplace_back();
    shape.name      = make_safe_name("", "light", (int)scene.shapes.size());
    shape.filename  = make_safe_filename("shapes/" + shape.name + ".ply");
    shape.frame     = plight.area_frame;
    shape.triangles = plight.area_triangles;
    shape.positions = plight.area_positions;
    shape.normals   = plight.area_normals;
    shape.material.emission = plight.area_emission;
  }

  // load shapes
  for (auto& shape : scene.shapes) {
    if (!shape.positions.empty()) continue;
    load_shape(filename, shape);
  }

  // load textures
  for (auto& texture : scene.textures) {
    if (!texture.hdr.empty() || !texture.ldr.empty()) continue;
    load_texture(filename, texture);
  }

  // fix scene
  scene.name = get_basename(filename);
  add_cameras(scene);
  add_radius(scene);
}

// Save a pbrt scene
void save_pbrt_scene(
    const string& filename, const sceneio_model& scene, bool noparallel) {
  // save pbrt
  auto pbrt = pbrt_model{};

  // embed data
  for (auto stat : scene_stats(scene)) pbrt.comments.push_back(stat);

  // convert camera
  auto& camera     = scene.cameras.front();
  auto& pcamera    = pbrt.cameras.emplace_back();
  pcamera.frame    = camera.frame;
  pcamera.lens     = camera.lens;
  pcamera.aspect   = camera.aspect;
  auto& pfilm      = pbrt.films.emplace_back();
  pfilm.filename   = "out.png";
  pfilm.resolution = {1280, (int)(1280 / pcamera.aspect)};

  // convert textures
  for (auto& texture : scene.textures) {
    auto& ptexture    = pbrt.textures.emplace_back();
    ptexture.name     = texture.name;
    ptexture.filename = texture.filename;
  }

  // convert materials
  for (auto& shape : scene.shapes) {
    auto& material    = shape.material;
    auto& pmaterial   = pbrt.materials.emplace_back();
    pmaterial.name    = shape.name;
    pmaterial.diffuse = material.base;
    // TODO: better conversion
    // pmaterial.specular     = material.specular;
    // pmaterial.transmission = material.transmission;
    pmaterial.roughness = {material.roughness, material.roughness};
    pmaterial.diffuse_map =
        material.base_tex >= 0 ? scene.textures[material.base_tex].name : ""s;
    auto& parealight    = pbrt.arealights.emplace_back();
    parealight.name     = shape.name;
    parealight.emission = material.emission;
  }

  // convert instances
  for (auto& shape : scene.shapes) {
    auto& material         = shape.material;
    auto& pshape           = pbrt.shapes.emplace_back();
    pshape.filename        = replace_extension(shape.filename, ".ply");
    pshape.frame           = shape.frame;
    pshape.material        = shape.name;
    pshape.arealight       = material.emission == zero3f ? ""s : shape.name;
    pshape.instance_frames = shape.instances;
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment    = pbrt.environments.emplace_back();
    penvironment.emission = environment.emission;
    if (environment.emission_tex >= 0) {
      penvironment.filename = scene.textures[environment.emission_tex].filename;
    }
  }

  // save pbrt
  save_pbrt(filename, pbrt);

  // save meshes
  auto dirname = get_dirname(filename);
  for (auto& shape : scene.shapes) {
    save_shape(replace_extension(dirname + shape.filename, ".ply"),
        shape.points, shape.lines, shape.triangles, shape.quads,
        shape.positions, shape.normals, shape.texcoords, shape.colors,
        shape.radius);
  }

  // save textures
  for (auto& texture : scene.textures) {
    save_texture(filename, texture);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace yocto {

void make_cornellbox_scene(sceneio_model& scene) {
  scene.name                = "cornellbox";
  auto& camera              = scene.cameras.emplace_back();
  camera.name               = "camera";
  camera.frame              = frame3f{{0, 1, 3.9}};
  camera.lens               = 0.035;
  camera.aperture           = 0.0;
  camera.film               = 0.024;
  camera.aspect             = 1;
  auto& floor_shp           = scene.shapes.emplace_back();
  floor_shp.name            = "floor";
  floor_shp.filename        = "shapes/floor.obj";
  floor_shp.positions       = {{-1, 0, 1}, {1, 0, 1}, {1, 0, -1}, {-1, 0, -1}};
  floor_shp.triangles       = {{0, 1, 2}, {2, 3, 0}};
  floor_shp.material.base   = {0.725, 0.71, 0.68};
  auto& ceiling_shp         = scene.shapes.emplace_back();
  ceiling_shp.name          = "ceiling";
  ceiling_shp.name          = "shapes/ceiling.obj";
  ceiling_shp.positions     = {{-1, 2, 1}, {-1, 2, -1}, {1, 2, -1}, {1, 2, 1}};
  ceiling_shp.triangles     = {{0, 1, 2}, {2, 3, 0}};
  ceiling_shp.material.base = {0.725, 0.71, 0.68};
  auto& backwall_shp        = scene.shapes.emplace_back();
  backwall_shp.name         = "backwall";
  backwall_shp.filename     = "shapes/backwall.obj";
  backwall_shp.positions = {{-1, 0, -1}, {1, 0, -1}, {1, 2, -1}, {-1, 2, -1}};
  backwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
  backwall_shp.material.base  = {0.725, 0.71, 0.68};
  auto& rightwall_shp         = scene.shapes.emplace_back();
  rightwall_shp.name          = "rightwall";
  rightwall_shp.filename      = "shapes/rightwall.obj";
  rightwall_shp.positions     = {{1, 0, -1}, {1, 0, 1}, {1, 2, 1}, {1, 2, -1}};
  rightwall_shp.triangles     = {{0, 1, 2}, {2, 3, 0}};
  rightwall_shp.material.base = {0.14, 0.45, 0.091};
  auto& leftwall_shp          = scene.shapes.emplace_back();
  leftwall_shp.name           = "leftwall";
  leftwall_shp.filename       = "shapes/leftwall.obj";
  leftwall_shp.positions = {{-1, 0, 1}, {-1, 0, -1}, {-1, 2, -1}, {-1, 2, 1}};
  leftwall_shp.triangles = {{0, 1, 2}, {2, 3, 0}};
  leftwall_shp.material.base  = {0.63, 0.065, 0.05};
  auto& shortbox_shp          = scene.shapes.emplace_back();
  shortbox_shp.name           = "shortbox";
  shortbox_shp.filename       = "shapes/shortbox.obj";
  shortbox_shp.positions      = {{0.53, 0.6, 0.75}, {0.7, 0.6, 0.17},
      {0.13, 0.6, 0.0}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {-0.05, 0.6, 0.57}, {0.13, 0.6, 0.0}, {0.13, 0.0, 0.0}, {0.53, 0.0, 0.75},
      {0.53, 0.6, 0.75}, {-0.05, 0.6, 0.57}, {-0.05, 0.0, 0.57},
      {0.7, 0.0, 0.17}, {0.7, 0.6, 0.17}, {0.53, 0.6, 0.75}, {0.53, 0.0, 0.75},
      {0.13, 0.0, 0.0}, {0.13, 0.6, 0.0}, {0.7, 0.6, 0.17}, {0.7, 0.0, 0.17},
      {0.53, 0.0, 0.75}, {0.7, 0.0, 0.17}, {0.13, 0.0, 0.0},
      {-0.05, 0.0, 0.57}};
  shortbox_shp.triangles      = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  shortbox_shp.material.base  = {0.725, 0.71, 0.68};
  auto& tallbox_shp           = scene.shapes.emplace_back();
  tallbox_shp.name            = "tallbox";
  tallbox_shp.filename        = "shapes/tallbox.obj";
  tallbox_shp.positions       = {{-0.53, 1.2, 0.09}, {0.04, 1.2, -0.09},
      {-0.14, 1.2, -0.67}, {-0.71, 1.2, -0.49}, {-0.53, 0.0, 0.09},
      {-0.53, 1.2, 0.09}, {-0.71, 1.2, -0.49}, {-0.71, 0.0, -0.49},
      {-0.71, 0.0, -0.49}, {-0.71, 1.2, -0.49}, {-0.14, 1.2, -0.67},
      {-0.14, 0.0, -0.67}, {-0.14, 0.0, -0.67}, {-0.14, 1.2, -0.67},
      {0.04, 1.2, -0.09}, {0.04, 0.0, -0.09}, {0.04, 0.0, -0.09},
      {0.04, 1.2, -0.09}, {-0.53, 1.2, 0.09}, {-0.53, 0.0, 0.09},
      {-0.53, 0.0, 0.09}, {0.04, 0.0, -0.09}, {-0.14, 0.0, -0.67},
      {-0.71, 0.0, -0.49}};
  tallbox_shp.triangles       = {{0, 1, 2}, {2, 3, 0}, {4, 5, 6}, {6, 7, 4},
      {8, 9, 10}, {10, 11, 8}, {12, 13, 14}, {14, 15, 12}, {16, 17, 18},
      {18, 19, 16}, {20, 21, 22}, {22, 23, 20}};
  tallbox_shp.material.base   = {0.725, 0.71, 0.68};
  auto& light_shp             = scene.shapes.emplace_back();
  light_shp.name              = "light";
  light_shp.filename          = "shapes/light.obj";
  light_shp.positions         = {{-0.25, 1.99, 0.25}, {-0.25, 1.99, -0.25},
      {0.25, 1.99, -0.25}, {0.25, 1.99, 0.25}};
  light_shp.triangles         = {{0, 1, 2}, {2, 3, 0}};
  light_shp.material.emission = {17, 12, 4};
}

}  // namespace yocto
