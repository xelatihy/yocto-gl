//
// Implementation for Yocto/Scene Input and Output functions.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include "yocto_sceneio.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <deque>
#include <filesystem>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ext/cgltf.h"
#include "ext/json.hpp"
#include "yocto_commonio.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_parallel.h"
#include "yocto_shape.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::deque;
using std::unique_ptr;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);
static bool save_json_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);
static bool save_obj_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);
static bool save_ply_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);
static bool save_pbrt_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel);

// Load a scene
bool load_scene(const string& filename, scene_model* scene, string& error,
    progress_callback progress_cb, bool noparallel) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return load_json_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return load_obj_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return load_gltf_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return load_pbrt_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return load_ply_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    return format_error();
  }
}

// Save a scene
bool save_scene(const string& filename, const scene_model* scene, string& error,
    progress_callback progress_cb, bool noparallel) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };

  auto ext = path_extension(filename);
  if (ext == ".json" || ext == ".JSON") {
    return save_json_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".obj" || ext == ".OBJ") {
    return save_obj_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".pbrt" || ext == ".PBRT") {
    return save_pbrt_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INDIVIDUAL ELEMENTS
// -----------------------------------------------------------------------------
namespace yocto {

// load instances
static bool load_instance(
    const string& filename, vector<frame3f>& frames, string& error) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    if (!load_ply(filename, &ply, error)) return false;
    get_values(&ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    return true;
  } else {
    return format_error();
  }
}

// save instances
// static
bool save_instance(const string& filename, const vector<frame3f>& frames,
    string& error, bool ascii = false) {
  auto format_error = [filename, &error]() {
    error = filename + ": unknown format";
    return false;
  };
  auto ext = path_extension(filename);
  if (ext == ".ply" || ext == ".PLY") {
    auto ply = ply_model{};
    add_values(&ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    if (!save_ply(filename, &ply, error)) return false;
    return true;
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

using json = nlohmann::json;
using std::array;

// support for json conversions
inline void to_json(json& j, const vec3f& value) {
  nlohmann::to_json(j, (const array<float, 3>&)value);
}
inline void to_json(json& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}

inline void from_json(const json& j, vec3f& value) {
  nlohmann::from_json(j, (array<float, 3>&)value);
}
inline void from_json(const json& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const json& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

using json = nlohmann::json;

// load/save json
static bool load_json(const string& filename, json& js, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    js = json::parse(text);
    return true;
  } catch (std::exception&) {
    return parse_error();
  }
}

static bool save_json(const string& filename, const json& js, string& error) {
  return save_text(filename, js.dump(2), error);
}

// Save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto material_error = [filename, &error](const string& name) {
    error = filename + ": missing material " + name;
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // open file
  auto js = json{};
  if (!load_json(filename, js, error)) return false;

  // gets a json value
  auto get_value = [](const json& ejs, const string& name,
                       auto& value) -> bool {
    if (!ejs.contains(name)) return true;
    try {
      ejs.at(name).get_to(value);
      return true;
    } catch (...) {
      return false;
    }
  };

  // parse json reference
  auto get_ref = [&material_error, &get_value](const json& ejs,
                     const string& name, auto& value,
                     const auto& refs) -> bool {
    if (!ejs.contains(name)) return true;
    auto ref = ""s;
    if (!get_value(ejs, name, ref)) return false;
    if (ref == "") {
      value = nullptr;
    } else {
      if (refs.find(ref) == refs.end()) return material_error(ref);
      value = refs.at(ref);
    }
    return true;
  };

  // parse json reference
  auto ctexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_ctexture = [scene, &ctexture_map, &get_value](const json& ejs,
                          const string& name, scene_texture*& value,
                          const string& dirname = "textures/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) {
      value = it->second;
      return true;
    }
    auto texture       = add_texture(scene, path);
    ctexture_map[path] = texture;
    value              = texture;
    return true;
  };

  // parse json reference
  auto stexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_stexture = [scene, &stexture_map, &get_value](const json& ejs,
                          const string& name, scene_texture*& value,
                          const string& dirname = "textures/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) {
      value = it->second;
      return true;
    }
    auto texture       = add_texture(scene, path);
    stexture_map[path] = texture;
    value              = texture;
    return true;
  };

  // parse json reference
  auto shape_map = unordered_map<string, scene_shape*>{{"", nullptr}};
  auto get_shape = [scene, &shape_map, &get_value](const json& ejs,
                       const string& name, scene_shape*& value,
                       const string& dirname = "shapes/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = shape_map.find(path);
    if (it != shape_map.end()) {
      value = it->second;
      return it->second;
    }
    auto shape      = add_shape(scene, path);
    shape_map[path] = shape;
    value           = shape;
    return true;
  };

  struct ply_instance {
    vector<frame3f> frames = {};
  };

  // load json instance
  auto ply_instances     = vector<unique_ptr<ply_instance>>{};
  auto ply_instance_map  = unordered_map<string, ply_instance*>{{"", nullptr}};
  auto instance_ply      = unordered_map<scene_instance*, ply_instance*>{};
  auto get_ply_instances = [&ply_instances, &ply_instance_map, &instance_ply,
                               &get_value](const json& ejs, const string& name,
                               scene_instance* instance,
                               const string&   dirname = "instances/") -> bool {
    if (!ejs.contains(name)) return true;
    auto path = ""s;
    if (!get_value(ejs, name, path)) return false;
    if (path == "") return true;
    auto it = ply_instance_map.find(path);
    if (it != ply_instance_map.end()) {
      instance_ply[instance] = it->second;
      return true;
    }
    auto ply_instance_ = ply_instances.emplace_back(new ply_instance()).get();
    ply_instance_map[path] = ply_instance_;
    instance_ply[instance] = ply_instance_;
    return true;
  };

  // material map
  auto material_map = unordered_map<string, scene_material*>{{"", nullptr}};

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // asset
  if (js.contains("asset")) {
    auto& ejs = js.at("asset");
    if (!get_value(ejs, "copyright", scene->copyright)) return false;
  }

  // cameras
  if (js.contains("cameras")) {
    for (auto& [name, ejs] : js.at("cameras").items()) {
      auto camera  = add_camera(scene);
      camera->name = name;
      if (!get_value(ejs, "frame", camera->frame)) return false;
      if (!get_value(ejs, "orthographic", camera->orthographic)) return false;
      if (!get_value(ejs, "lens", camera->lens)) return false;
      if (!get_value(ejs, "aspect", camera->aspect)) return false;
      if (!get_value(ejs, "film", camera->film)) return false;
      if (!get_value(ejs, "focus", camera->focus)) return false;
      if (!get_value(ejs, "aperture", camera->aperture)) return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        camera->frame = lookat_frame(lookat.x, lookat.y, lookat.z);
        camera->focus = length(lookat.x - lookat.y);
      }
    }
  }
  if (js.contains("environments")) {
    for (auto& [name, ejs] : js.at("environments").items()) {
      auto environment  = add_environment(scene);
      environment->name = name;
      if (!get_value(ejs, "frame", environment->frame)) return false;
      if (!get_value(ejs, "emission", environment->emission)) return false;
      if (!get_ctexture(
              ejs, "emission_tex", environment->emission_tex, "environments/"))
        return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        environment->frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
    }
  }
  if (js.contains("materials")) {
    for (auto& [name, ejs] : js.at("materials").items()) {
      auto material  = add_material(scene);
      material->name = name;
      if (!get_value(ejs, "emission", material->emission)) return false;
      if (!get_value(ejs, "color", material->color)) return false;
      if (!get_value(ejs, "metallic", material->metallic)) return false;
      if (!get_value(ejs, "specular", material->specular)) return false;
      if (!get_value(ejs, "roughness", material->roughness)) return false;
      if (!get_value(ejs, "coat", material->coat)) return false;
      if (!get_value(ejs, "transmission", material->transmission)) return false;
      if (!get_value(ejs, "translucency", material->translucency)) return false;
      if (!get_value(ejs, "thin", material->thin)) return false;
      if (!get_value(ejs, "ior", material->ior)) return false;
      if (!get_value(ejs, "trdepth", material->trdepth)) return false;
      if (!get_value(ejs, "scattering", material->scattering)) return false;
      if (!get_value(ejs, "scanisotropy", material->scanisotropy)) return false;
      if (!get_value(ejs, "opacity", material->opacity)) return false;
      if (!get_value(ejs, "coat", material->coat)) return false;
      if (!get_ctexture(ejs, "emission_tex", material->emission_tex))
        return false;
      if (!get_ctexture(ejs, "color_tex", material->color_tex)) return false;
      if (!get_stexture(ejs, "metallic_tex", material->metallic_tex))
        return false;
      if (!get_stexture(ejs, "specular_tex", material->specular_tex))
        return false;
      if (!get_stexture(ejs, "transmission_tex", material->transmission_tex))
        return false;
      if (!get_stexture(ejs, "translucency_tex", material->translucency_tex))
        return false;
      if (!get_stexture(ejs, "roughness_tex", material->roughness_tex))
        return false;
      if (!get_ctexture(ejs, "scattering_tex", material->scattering_tex))
        return false;
      if (!get_stexture(ejs, "opacity_tex", material->opacity_tex))
        return false;
      if (!get_ctexture(ejs, "normal_tex", material->normal_tex)) return false;
      material_map[material->name] = material;
    }
  }
  if (js.contains("instances")) {
    for (auto& [name, ejs] : js.at("instances").items()) {
      auto instance  = add_instance(scene);
      instance->name = name;
      if (!get_value(ejs, "frame", instance->frame)) return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        instance->frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
      if (!get_ref(ejs, "material", instance->material, material_map))
        return false;
      if (!get_shape(ejs, "shape", instance->shape)) return false;
      if (!get_ply_instances(ejs, "instance", instance)) return false;
      if (instance->shape) {
        if (!get_value(ejs, "subdivisions", instance->shape->subdivisions))
          return false;
        if (!get_value(ejs, "catmullcark", instance->shape->catmullclark))
          return false;
        if (!get_value(ejs, "smooth", instance->shape->smooth)) return false;
        if (!get_value(ejs, "displacement", instance->shape->displacement))
          return false;
        if (!get_stexture(
                ejs, "displacement_tex", instance->shape->displacement_tex))
          return false;
      }
    }
  }
  if (js.contains("objects")) {
    for (auto& [name, ejs] : js.at("objects").items()) {
      auto instance  = add_instance(scene);
      instance->name = name;
      if (!get_value(ejs, "frame", instance->frame)) return false;
      if (ejs.contains("lookat")) {
        auto lookat = identity3x3f;
        if (!get_value(ejs, "lookat", lookat)) return false;
        instance->frame = lookat_frame(lookat.x, lookat.y, lookat.z, true);
      }
      if (!get_ref(ejs, "material", instance->material, material_map))
        return false;
      if (!get_shape(ejs, "shape", instance->shape)) return false;
      if (!get_ply_instances(ejs, "instance", instance)) return false;
      if (instance->shape) {
        if (!get_value(ejs, "subdivisions", instance->shape->subdivisions))
          return false;
        if (!get_value(ejs, "catmullcark", instance->shape->catmullclark))
          return false;
        if (!get_value(ejs, "smooth", instance->shape->smooth)) return false;
        if (!get_value(ejs, "displacement", instance->shape->displacement))
          return false;
        if (!get_stexture(
                ejs, "displacement_tex", instance->shape->displacement_tex))
          return false;
      }
    }
  }

  // handle progress
  progress.y += scene->shapes.size();
  progress.y += scene->textures.size();
  progress.y += ply_instances.size();

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const vector<string>& extensions) {
    for (auto& extension : extensions) {
      auto filepath = path_join(
          path_dirname(filename), group, name + extension);
      if (path_exists(filepath)) return filepath;
    }
    return path_join(path_dirname(filename), group, name + extensions.front());
  };

  // load shapes
  shape_map.erase("");
  for (auto [name, shape] : shape_map) {
    if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
    auto path = make_filename(name, "shapes", {".ply", ".obj"});
    if (!load_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->quadspos, shape->quadsnorm,
            shape->quadstexcoord, shape->positions, shape->normals,
            shape->texcoords, shape->colors, shape->radius, error,
            shape->catmullclark && shape->subdivisions))
      return dependent_error();
  }
  // load textures
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto path = make_filename(
        name, "textures", {".hdr", ".exr", ".png", ".jpg"});
    if (!load_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }
  // load textures
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto path = make_filename(
        name, "textures", {".hdr", ".exr", ".png", ".jpg"});
    if (!load_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // load instances
  ply_instance_map.erase("");
  for (auto [name, instance] : ply_instance_map) {
    if (progress_cb) progress_cb("load instance", progress.x++, progress.y);
    auto path = make_filename(name, "instances", {".ply"});
    if (!load_instance(path, instance->frames, error)) return dependent_error();
  }

  // apply instances
  if (!ply_instances.empty()) {
    if (progress_cb)
      progress_cb("flatten instances", progress.x++, progress.y++);
    auto instances = scene->instances;
    scene->instances.clear();
    for (auto instance : instances) {
      auto it = instance_ply.find(instance);
      if (it == instance_ply.end()) {
        auto ninstance      = add_instance(scene, instance->name);
        ninstance->frame    = instance->frame;
        ninstance->shape    = instance->shape;
        ninstance->material = instance->material;
      } else {
        auto ply_instance = it->second;
        for (auto& frame : ply_instance->frames) {
          auto ninstance      = add_instance(scene, instance->name);
          ninstance->frame    = frame * instance->frame;
          ninstance->shape    = instance->shape;
          ninstance->material = instance->material;
        }
      }
    }
    for (auto instance : instances) delete instance;
  }

  // fix scene
  if (scene->name == "") scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);
  trim_memory(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // helper
  auto add_opt = [](json& ejs, const string& name, const auto& value,
                     const auto& def) {
    if (value == def) return;
    ejs[name] = value;
  };
  auto add_tex = [](json& ejs, const string& name, scene_texture* texture) {
    if (!texture) return;
    ejs[name] = texture->name;
  };
  auto add_ref = [](json& ejs, const string& name, auto ref) {
    if (!ref) return;
    ejs[name] = ref->name;
  };

  // handle progress
  auto progress = vec2i{
      0, 2 + (int)scene->shapes.size() + (int)scene->textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json file
  auto js = json::object();

  // asset
  {
    auto& ejs        = js["asset"];
    ejs["generator"] = "Yocto/GL - https://github.com/xelatihy/yocto-gl";
    add_opt(ejs, "copyright", scene->copyright, ""s);
  }

  auto def_cam = scene_camera{};
  if (!scene->cameras.empty()) js["cameras"] = json::object();
  for (auto& camera : scene->cameras) {
    auto& ejs = js["cameras"][camera->name];
    add_opt(ejs, "frame", camera->frame, def_cam.frame);
    add_opt(ejs, "ortho", camera->orthographic, def_cam.orthographic);
    add_opt(ejs, "lens", camera->lens, def_cam.lens);
    add_opt(ejs, "aspect", camera->aspect, def_cam.aspect);
    add_opt(ejs, "film", camera->film, def_cam.film);
    add_opt(ejs, "focus", camera->focus, def_cam.focus);
    add_opt(ejs, "aperture", camera->aperture, def_cam.aperture);
  }

  auto def_env = scene_environment{};
  if (!scene->environments.empty()) js["environments"] = json::object();
  for (auto environment : scene->environments) {
    auto& ejs = js["environments"][environment->name];
    add_opt(ejs, "frame", environment->frame, def_env.frame);
    add_opt(ejs, "emission", environment->emission, def_env.emission);
    add_tex(ejs, "emission_tex", environment->emission_tex);
  }

  auto def_material = scene_material{};
  if (!scene->materials.empty()) js["materials"] = json::object();
  for (auto material : scene->materials) {
    auto& ejs = js["materials"][material->name];
    add_opt(ejs, "emission", material->emission, def_material.emission);
    add_opt(ejs, "color", material->color, def_material.color);
    add_opt(ejs, "specular", material->specular, def_material.specular);
    add_opt(ejs, "metallic", material->metallic, def_material.metallic);
    add_opt(ejs, "coat", material->coat, def_material.coat);
    add_opt(ejs, "roughness", material->roughness, def_material.roughness);
    add_opt(ejs, "ior", material->ior, def_material.ior);
    add_opt(
        ejs, "transmission", material->transmission, def_material.transmission);
    add_opt(
        ejs, "translucency", material->translucency, def_material.translucency);
    add_opt(ejs, "trdepth", material->trdepth, def_material.trdepth);
    add_opt(ejs, "scattering", material->scattering, def_material.scattering);
    add_opt(
        ejs, "scanisotropy", material->scanisotropy, def_material.scanisotropy);
    add_opt(ejs, "opacity", material->opacity, def_material.opacity);
    add_opt(ejs, "thin", material->thin, def_material.thin);
    add_tex(ejs, "emission_tex", material->emission_tex);
    add_tex(ejs, "color_tex", material->color_tex);
    add_tex(ejs, "metallic_tex", material->metallic_tex);
    add_tex(ejs, "specular_tex", material->specular_tex);
    add_tex(ejs, "roughness_tex", material->roughness_tex);
    add_tex(ejs, "transmission_tex", material->transmission_tex);
    add_tex(ejs, "translucency_tex", material->translucency_tex);
    add_tex(ejs, "scattering_tex", material->scattering_tex);
    add_tex(ejs, "coat_tex", material->coat_tex);
    add_tex(ejs, "opacity_tex", material->opacity_tex);
    add_tex(ejs, "normal_tex", material->normal_tex);
  }

  auto def_object = scene_instance{};
  auto def_shape  = scene_shape{};
  if (!scene->instances.empty()) js["instances"] = json::object();
  for (auto instance : scene->instances) {
    auto& ejs = js["instances"][instance->name];
    add_opt(ejs, "frame", instance->frame, def_object.frame);
    add_ref(ejs, "shape", instance->shape);
    add_ref(ejs, "material", instance->material);
    if (instance->shape) {
      add_opt(ejs, "subdivisions", instance->shape->subdivisions,
          def_shape.subdivisions);
      add_opt(ejs, "catmullclark", instance->shape->catmullclark,
          def_shape.catmullclark);
      add_opt(ejs, "smooth", instance->shape->smooth, def_shape.smooth);
      add_opt(ejs, "displacement", instance->shape->displacement,
          def_shape.displacement);
      add_tex(ejs, "displacement_tex", instance->shape->displacement_tex);
    }
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json
  if (!save_json(filename, js, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // save shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    auto path = make_filename(shape->name, "shapes",
        (shape->catmullclark && shape->subdivisions) ? ".obj" : ".ply");
    if (!save_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->quadspos, shape->quadsnorm,
            shape->quadstexcoord, shape->positions, shape->normals,
            shape->texcoords, shape->colors, shape->radius, error,
            shape->catmullclark && shape->subdivisions))
      return dependent_error();
  }

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr" : ".png");
    if (!save_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Loads an OBJ
static bool load_obj_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load obj
  auto obj_guard = std::make_unique<obj_model>();
  auto obj       = obj_guard.get();
  if (!load_obj(filename, obj, error, false, true, false)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto ocam : obj->cameras) {
    auto camera = add_camera(scene);
    // camera->name         = make_safe_name("camera", ocam->name);
    camera->frame        = ocam->frame;
    camera->orthographic = ocam->ortho;
    camera->film         = max(ocam->width, ocam->height);
    camera->aspect       = ocam->width / ocam->height;
    camera->focus        = ocam->focus;
    camera->lens         = ocam->lens;
    camera->aperture     = ocam->aperture;
  }

  // helper to create texture maps
  auto ctexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_ctexture = [&ctexture_map, scene](
                          const obj_texture& tinfo) -> scene_texture* {
    auto path = tinfo.path;
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };

  // helper to create texture maps
  auto stexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_stexture = [&stexture_map, scene](
                          const obj_texture& tinfo) -> scene_texture* {
    auto path = tinfo.path;
    if (path == "") return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
    return texture;
  };

  // handler for materials
  auto material_map = unordered_map<obj_material*, scene_material*>{};
  for (auto omat : obj->materials) {
    auto material = add_material(scene);
    // material->name             = make_safe_name("material", omat->name);
    material->emission         = omat->pbr_emission;
    material->color            = omat->pbr_base;
    material->specular         = omat->pbr_specular;
    material->roughness        = omat->pbr_roughness;
    material->ior              = omat->pbr_ior;
    material->metallic         = omat->pbr_metallic;
    material->coat             = omat->pbr_coat;
    material->transmission     = omat->pbr_transmission;
    material->translucency     = omat->pbr_translucency;
    material->scattering       = omat->pbr_volscattering;
    material->scanisotropy     = omat->pbr_volanisotropy;
    material->trdepth          = omat->pbr_volscale;
    material->opacity          = omat->pbr_opacity;
    material->thin             = omat->pbr_thin;
    material->emission_tex     = get_ctexture(omat->pbr_emission_tex);
    material->color_tex        = get_ctexture(omat->pbr_base_tex);
    material->specular_tex     = get_stexture(omat->pbr_specular_tex);
    material->metallic_tex     = get_stexture(omat->pbr_metallic_tex);
    material->roughness_tex    = get_stexture(omat->pbr_roughness_tex);
    material->transmission_tex = get_stexture(omat->pbr_transmission_tex);
    material->translucency_tex = get_stexture(omat->pbr_translucency_tex);
    material->coat_tex         = get_stexture(omat->pbr_coat_tex);
    material->opacity_tex      = get_stexture(omat->pbr_opacity_tex);
    material->normal_tex       = get_ctexture(omat->normal_tex);
    material->scattering_tex   = get_ctexture(omat->pbr_volscattering_tex);
    material_map[omat]         = material;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto oshape : obj->shapes) {
    auto& materials = oshape->materials;
    if (materials.empty()) materials.push_back(nullptr);
    for (auto material_idx = 0; material_idx < materials.size();
         material_idx++) {
      auto shape      = add_shape(scene);
      auto material   = material_map.at(materials[material_idx]);
      auto has_quads_ = has_quads(oshape);
      if (!oshape->faces.empty() && !has_quads_) {
        get_triangles(oshape, material_idx, shape->triangles, shape->positions,
            shape->normals, shape->texcoords, true);
      } else if (!oshape->faces.empty() && has_quads_) {
        get_quads(oshape, material_idx, shape->quads, shape->positions,
            shape->normals, shape->texcoords, true);
      } else if (!oshape->lines.empty()) {
        get_lines(oshape, material_idx, shape->lines, shape->positions,
            shape->normals, shape->texcoords, true);
      } else if (!oshape->points.empty()) {
        get_points(oshape, material_idx, shape->points, shape->positions,
            shape->normals, shape->texcoords, true);
      } else {
        return shape_error();
      }
      if (oshape->instances.empty()) {
        auto instance      = add_instance(scene);
        instance->shape    = shape;
        instance->material = material;
      } else {
        for (auto& frame : oshape->instances) {
          auto instance      = add_instance(scene);
          instance->frame    = frame;
          instance->shape    = shape;
          instance->material = material;
        }
      }
    }
  }

  // convert environments
  for (auto oenvironment : obj->environments) {
    auto environment          = add_environment(scene);
    environment->frame        = oenvironment->frame;
    environment->emission     = oenvironment->emission;
    environment->emission_tex = get_ctexture(oenvironment->emission_tex);
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name) {
    return path_join(path_dirname(filename), name);
  };

  // load textures
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // load textures
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // fix scene
  if (scene->name == "") scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_obj_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2 + (int)scene->textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  auto obj_guard = std::make_unique<obj_model>();
  auto obj       = obj_guard.get();

  // convert cameras
  for (auto camera : scene->cameras) {
    auto ocamera      = add_camera(obj);
    ocamera->name     = path_basename(camera->name);
    ocamera->frame    = camera->frame;
    ocamera->ortho    = camera->orthographic;
    ocamera->width    = camera->film;
    ocamera->height   = camera->film / camera->aspect;
    ocamera->focus    = camera->focus;
    ocamera->lens     = camera->lens;
    ocamera->aperture = camera->aperture;
  }

  // textures
  auto get_texture = [](scene_texture* texture) {
    if (!texture) return obj_texture{};
    auto tinfo = obj_texture{};
    tinfo.path = texture->name;
    return tinfo;
  };

  // convert materials and textures
  auto material_map = unordered_map<scene_material*, obj_material*>{
      {nullptr, nullptr}};
  for (auto material : scene->materials) {
    auto omaterial                  = add_material(obj);
    omaterial->name                 = path_basename(material->name);
    omaterial->illum                = 2;
    omaterial->as_pbr               = true;
    omaterial->pbr_emission         = material->emission;
    omaterial->pbr_base             = material->color;
    omaterial->pbr_specular         = material->specular;
    omaterial->pbr_roughness        = material->roughness;
    omaterial->pbr_metallic         = material->metallic;
    omaterial->pbr_coat             = material->coat;
    omaterial->pbr_transmission     = material->transmission;
    omaterial->pbr_translucency     = material->translucency;
    omaterial->pbr_opacity          = material->opacity;
    omaterial->pbr_emission_tex     = get_texture(material->emission_tex);
    omaterial->pbr_base_tex         = get_texture(material->color_tex);
    omaterial->pbr_specular_tex     = get_texture(material->specular_tex);
    omaterial->pbr_metallic_tex     = get_texture(material->metallic_tex);
    omaterial->pbr_roughness_tex    = get_texture(material->roughness_tex);
    omaterial->pbr_transmission_tex = get_texture(material->transmission_tex);
    omaterial->pbr_translucency_tex = get_texture(material->translucency_tex);
    omaterial->pbr_coat_tex         = get_texture(material->coat_tex);
    omaterial->pbr_opacity_tex      = get_texture(material->opacity_tex);
    omaterial->pbr_normal_tex       = get_texture(material->normal_tex);
    material_map[material]          = omaterial;
  }

  // convert objects
  for (auto instance : scene->instances) {
    auto shape     = instance->shape;
    auto positions = shape->positions, normals = shape->normals;
    for (auto& p : positions) p = transform_point(instance->frame, p);
    for (auto& n : normals) n = transform_normal(instance->frame, n);
    auto oshape       = add_shape(obj);
    oshape->name      = shape->name;
    oshape->materials = {material_map.at(instance->material)};
    if (!shape->triangles.empty()) {
      set_triangles(oshape, shape->triangles, positions, normals,
          shape->texcoords, {}, true);
    } else if (!shape->quads.empty()) {
      set_quads(
          oshape, shape->quads, positions, normals, shape->texcoords, {}, true);
    } else if (!shape->lines.empty()) {
      set_lines(
          oshape, shape->lines, positions, normals, shape->texcoords, {}, true);
    } else if (!shape->points.empty()) {
      set_points(oshape, shape->points, positions, normals, shape->texcoords,
          {}, true);
    } else {
      return shape_error();
    }
  }

  // convert environments
  for (auto environment : scene->environments) {
    auto oenvironment          = add_environment(obj);
    oenvironment->name         = path_basename(environment->name);
    oenvironment->frame        = environment->frame;
    oenvironment->emission     = environment->emission;
    oenvironment->emission_tex = get_texture(environment->emission_tex);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save obj
  if (!save_obj(filename, obj, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr" : ".png");
    if (!save_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

void print_obj_camera(scene_camera* camera) {
  printf("c %s %d %g %g %g %g %g %g %g %g %g %g%g %g %g %g %g %g %g\n",
      camera->name.c_str(), (int)camera->orthographic, camera->film,
      camera->film / camera->aspect, camera->lens, camera->focus,
      camera->aperture, camera->frame.x.x, camera->frame.x.y, camera->frame.x.z,
      camera->frame.y.x, camera->frame.y.y, camera->frame.y.z,
      camera->frame.z.x, camera->frame.z.y, camera->frame.z.z,
      camera->frame.o.x, camera->frame.o.y, camera->frame.o.z);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_ply_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load ply mesh
  auto shape = add_shape(scene);
  if (!load_shape(filename, shape->points, shape->lines, shape->triangles,
          shape->quads, shape->quadspos, shape->quadsnorm, shape->quadstexcoord,
          shape->positions, shape->normals, shape->texcoords, shape->colors,
          shape->radius, error))
    return false;

  // create instance
  auto instance   = add_instance(scene);
  instance->shape = shape;

  // fix scene
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_ply_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  if (scene->shapes.empty()) return shape_error();

  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save shape
  auto shape = scene->shapes.front();
  if (!save_shape(filename, shape->points, shape->lines, shape->triangles,
          shape->quads, shape->quadspos, shape->quadsnorm, shape->quadstexcoord,
          shape->positions, shape->normals, shape->texcoords, shape->colors,
          shape->radius, error))
    return false;

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene
static bool load_gltf_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto primitive_error = [filename, &error]() {
    error = filename + ": primitive error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 3};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load gltf
  auto params = cgltf_options{};
  memset(&params, 0, sizeof(params));
  auto data   = (cgltf_data*)nullptr;
  auto result = cgltf_parse_file(&params, filename.c_str(), &data);
  if (result != cgltf_result_success) return read_error();
  auto gltf = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      data, cgltf_free};

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load buffers
  auto dirname = path_dirname(filename);
  if (dirname != "") dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success)
    return read_error();

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert asset
  {
    auto gast = &gltf->asset;
    if (gast->copyright) scene->copyright = gast->copyright;
  }

  // prepare list of effective nodes
  auto visible_nodes = vector<bool>(gltf->nodes_count, false);
  auto gscene        = gltf->scene ? gltf->scene : gltf->scenes;
  if (gscene) {
    auto node_index = unordered_map<cgltf_node*, int>{};
    node_index.reserve(gltf->nodes_count);
    for (auto nid = 0; nid < gltf->nodes_count; nid++)
      node_index[&gltf->nodes[nid]] = nid;
    auto stack = vector<cgltf_node*>{};
    for (auto nid = 0; nid < gscene->nodes_count; nid++)
      stack.push_back(gscene->nodes[nid]);
    while (!stack.empty()) {
      auto gnde = stack.back();
      stack.pop_back();
      visible_nodes[node_index[gnde]] = true;
      for (auto nid = 0; nid < gnde->children_count; nid++)
        stack.push_back(gnde->children[nid]);
    }
  } else {
    for (auto nid = 0; nid < gltf->nodes_count; nid++)
      visible_nodes[nid] = true;
  }

  // convert cameras
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    if (!visible_nodes[nid]) continue;
    auto gnde = &gltf->nodes[nid];
    if (!gnde->camera) continue;
    auto mat = mat4f{};
    cgltf_node_transform_world(gnde, &mat.x.x);
    auto gcam            = gnde->camera;
    auto camera          = add_camera(scene);
    camera->frame        = mat_to_frame(mat);
    camera->orthographic = gcam->type == cgltf_camera_type_orthographic;
    if (camera->orthographic) {
      auto ortho     = &gcam->data.orthographic;
      camera->aspect = ortho->xmag / ortho->ymag;
      camera->lens   = ortho->ymag;  // this is probably bogus
      camera->film   = 0.036;
    } else {
      auto persp     = &gcam->data.perspective;
      camera->aspect = persp->aspect_ratio;
      if (camera->aspect == 0) camera->aspect = 16.0f / 9.0f;
      camera->film = 0.036;
      if (camera->aspect >= 1) {
        camera->lens = (camera->film / camera->aspect) /
                       (2 * tan(persp->yfov / 2));
      } else {
        camera->lens = camera->film / (2 * tan(persp->yfov / 2));
      }
    }
  }

  // convert color textures
  auto ctexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const cgltf_texture_view& ginfo) -> scene_texture* {
    if (!ginfo.texture || !ginfo.texture->image) return nullptr;
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  // convert color opacity textures
  auto cotexture_map =
      unordered_map<string, pair<scene_texture*, scene_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_cotexture = [&scene, &cotexture_map](const cgltf_texture_view& ginfo)
      -> pair<scene_texture*, scene_texture*> {
    if (!ginfo.texture || !ginfo.texture->image) return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return {nullptr, nullptr};
    auto it = cotexture_map.find(path);
    if (it != cotexture_map.end()) return it->second;
    auto color_texture   = add_texture(scene);
    auto opacity_texture = add_texture(scene);
    cotexture_map[path]  = {color_texture, opacity_texture};
    return {color_texture, opacity_texture};
  };
  // convert textures
  auto mrtexture_map =
      unordered_map<string, pair<scene_texture*, scene_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_mrtexture = [&scene, &mrtexture_map](const cgltf_texture_view& ginfo)
      -> pair<scene_texture*, scene_texture*> {
    if (!ginfo.texture || !ginfo.texture->image) return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path == "") return {nullptr, nullptr};
    auto it = mrtexture_map.find(path);
    if (it != mrtexture_map.end()) return it->second;
    auto metallic_texture  = add_texture(scene);
    auto roughness_texture = add_texture(scene);
    mrtexture_map[path]    = {metallic_texture, roughness_texture};
    return {metallic_texture, roughness_texture};
  };

  // convert materials
  auto material_map = unordered_map<cgltf_material*, scene_material*>{
      {nullptr, nullptr}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto gmaterial         = &gltf->materials[mid];
    auto material          = add_material(scene);
    material->emission     = {gmaterial->emissive_factor[0],
        gmaterial->emissive_factor[1], gmaterial->emissive_factor[2]};
    material->emission_tex = get_ctexture(gmaterial->emissive_texture);
    if (gmaterial->has_pbr_metallic_roughness) {
      auto gmr          = &gmaterial->pbr_metallic_roughness;
      material->color   = {gmr->base_color_factor[0], gmr->base_color_factor[1],
          gmr->base_color_factor[2]};
      material->opacity = gmr->base_color_factor[3];
      material->metallic  = gmr->metallic_factor;
      material->roughness = gmr->roughness_factor;
      material->specular  = 1;
      std::tie(material->color_tex, material->opacity_tex) = get_cotexture(
          gmr->base_color_texture);
      std::tie(material->metallic_tex, material->roughness_tex) = get_mrtexture(
          gmr->metallic_roughness_texture);
    }
    material->normal_tex    = get_ctexture(gmaterial->normal_texture);
    material_map[gmaterial] = material;
  }

  // convert meshes
  auto mesh_map = unordered_map<cgltf_mesh*, vector<scene_instance*>>{
      {nullptr, {}}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh = &gltf->meshes[mid];
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (!gprim->attributes_count) continue;
      auto instance = add_instance(scene);
      mesh_map[gmesh].push_back(instance);
      auto shape         = add_shape(scene);
      instance->shape    = shape;
      instance->material = material_map.at(gprim->material);
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name ? gattr->name : "");
        auto gacc     = gattr->data;
        if (semantic == "POSITION") {
          shape->positions.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->positions[i].x, 3);
        } else if (semantic == "NORMAL") {
          shape->normals.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->normals[i].x, 3);
        } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
          shape->texcoords.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->texcoords[i].x, 2);
        } else if (semantic == "COLOR" || semantic == "COLOR_0") {
          shape->colors.resize(gacc->count);
          if (cgltf_num_components(gacc->type) == 3) {
            for (auto i = 0; i < gacc->count; i++)
              cgltf_accessor_read_float(gacc, i, &shape->colors[i].x, 3);
          } else {
            for (auto i = 0; i < gacc->count; i++) {
              auto color4 = vec4f{0, 0, 0, 0};
              cgltf_accessor_read_float(gacc, i, &color4.x, 4);
              shape->colors[i] = color4;
            }
          }
        } else if (semantic == "TANGENT") {
          shape->tangents.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->tangents[i].x, 4);
          for (auto& t : shape->tangents) t.w = -t.w;
        } else if (semantic == "RADIUS") {
          shape->radius.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->radius[i], 1);
        } else {
          // ignore
        }
      }
      // indices
      if (!gprim->indices) {
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape->triangles.resize(shape->positions.size() / 3);
          for (auto i = 0; i < shape->positions.size() / 3; i++)
            shape->triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape->triangles.resize(shape->positions.size() - 2);
          for (auto i = 2; i < shape->positions.size(); i++)
            shape->triangles[i - 2] = {0, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape->triangles.resize(shape->positions.size() - 2);
          for (auto i = 2; i < shape->positions.size(); i++)
            shape->triangles[i - 2] = {i - 2, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape->lines.resize(shape->positions.size() / 2);
          for (auto i = 0; i < shape->positions.size() / 2; i++)
            shape->lines[i] = {i * 2 + 0, i * 2 + 1};
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape->lines.resize(shape->positions.size());
          for (auto i = 1; i < shape->positions.size(); i++)
            shape->lines[i - 1] = {i - 1, i};
          shape->lines.back() = {(int)shape->positions.size() - 1, 0};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape->lines.resize(shape->positions.size() - 1);
          for (auto i = 1; i < shape->positions.size(); i++)
            shape->lines[i - 1] = {i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          return primitive_error();
        } else {
          return primitive_error();
        }
      } else {
        auto giacc = gprim->indices;
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape->triangles.resize(giacc->count / 3);
          for (auto i = 0; i < giacc->count / 3; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 3 + 0, (uint*)&shape->triangles[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 1, (uint*)&shape->triangles[i].y, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 2, (uint*)&shape->triangles[i].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape->triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, 0 + 0, (uint*)&shape->triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape->triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape->triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape->triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, i - 2, (uint*)&shape->triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape->triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape->triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape->lines.resize(giacc->count / 2);
          for (auto i = 0; i < giacc->count / 2; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 2 + 0, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 2 + 1, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape->lines.resize(giacc->count);
          for (auto i = 0; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape->lines.resize(giacc->count - 1);
          for (auto i = 0; i < giacc->count - 1; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape->lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape->lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          return primitive_error();
        } else {
          return primitive_error();
        }
      }
    }
  }

  // convert nodes
  auto instance_map = unordered_map<cgltf_mesh*, vector<frame3f>>{};
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    if (!visible_nodes[nid]) continue;
    auto gnde = &gltf->nodes[nid];
    if (!gnde->mesh) continue;
    auto mat = mat4f{};
    cgltf_node_transform_world(gnde, &mat.x.x);
    auto frame = mat_to_frame(mat);
    instance_map[gnde->mesh].push_back(frame);
  }
  for (auto& [gmsh, frames] : instance_map) {
    if (frames.size() == 1) {
      for (auto object : mesh_map.at(gmsh)) object->frame = frames.front();
    } else {
      for (auto object : mesh_map.at(gmsh)) {
        scene->instances.erase(std::remove(scene->instances.begin(),
                                   scene->instances.end(), object),
            scene->instances.end());
        for (auto fid = 0; fid < frames.size(); fid++) {
          auto nobject = add_instance(
              scene, object->name + "@" + std::to_string(fid));
          nobject->frame    = frames[fid];
          nobject->shape    = object->shape;
          nobject->material = object->material;
        }
        delete object;
      }
    }
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // load texture
  ctexture_map.erase("");
  for (auto [tpath, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(path_join(path_dirname(filename), tpath), texture->hdr,
            texture->ldr, error))
      return dependent_error();
  }

  // load texture
  cotexture_map.erase("");
  for (auto [tpath, textures] : cotexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto color_opacityf = image<vec4f>{};
    auto color_opacityb = image<vec4b>{};
    if (!load_image(path_join(path_dirname(filename), tpath), color_opacityf,
            color_opacityb, error))
      return dependent_error();
    if (!color_opacityf.empty()) {
      auto [ctexture, otexture] = textures;
      ctexture->hdr.resize(color_opacityf.imsize());
      otexture->hdr.resize(color_opacityf.imsize());
      auto oempty = true;
      for (auto j = 0; j < color_opacityf.height(); j++) {
        for (auto i = 0; i < color_opacityf.width(); i++) {
          auto color            = xyz(color_opacityf[{i, j}]);
          auto opacity          = color_opacityf[{i, j}].w;
          ctexture->hdr[{i, j}] = {color.x, color.y, color.z, opacity};
          otexture->hdr[{i, j}] = {opacity, opacity, opacity, opacity};
          if (opacity != 1) oempty = false;
        }
      }
      if (oempty) otexture->hdr.clear();
    }
    if (!color_opacityb.empty()) {
      auto [ctexture, otexture] = textures;
      ctexture->ldr.resize(color_opacityb.imsize());
      otexture->ldr.resize(color_opacityb.imsize());
      auto oempty = true;
      for (auto j = 0; j < color_opacityb.height(); j++) {
        for (auto i = 0; i < color_opacityb.width(); i++) {
          auto color            = xyz(color_opacityb[{i, j}]);
          auto opacity          = color_opacityb[{i, j}].w;
          ctexture->ldr[{i, j}] = {color.x, color.y, color.z, opacity};
          otexture->ldr[{i, j}] = {opacity, opacity, opacity, opacity};
          if (opacity != 1) oempty = false;
        }
      }
      if (oempty) otexture->ldr.clear();
    }
  }

  // load texture
  mrtexture_map.erase("");
  for (auto [tpath, textures] : mrtexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto metallic_roughnessf = image<vec4f>{};
    auto metallic_roughnessb = image<vec4b>{};
    if (!load_image(path_join(path_dirname(filename), tpath),
            metallic_roughnessf, metallic_roughnessb, error))
      return dependent_error();
    if (!metallic_roughnessf.empty()) {
      auto [mtexture, rtexture] = textures;
      mtexture->hdr.resize(metallic_roughnessf.imsize());
      rtexture->hdr.resize(metallic_roughnessf.imsize());
      for (auto j = 0; j < metallic_roughnessf.height(); j++) {
        for (auto i = 0; i < metallic_roughnessf.width(); i++) {
          auto metallic         = metallic_roughnessf[{i, j}].z;
          auto roughness        = metallic_roughnessf[{i, j}].y;
          mtexture->hdr[{i, j}] = {metallic, metallic, metallic, 1};
          rtexture->hdr[{i, j}] = {roughness, roughness, roughness, 1};
        }
      }
    }
    if (!metallic_roughnessb.empty()) {
      auto [mtexture, rtexture] = textures;
      mtexture->ldr.resize(metallic_roughnessb.imsize());
      rtexture->ldr.resize(metallic_roughnessb.imsize());
      for (auto j = 0; j < metallic_roughnessb.height(); j++) {
        for (auto i = 0; i < metallic_roughnessb.width(); i++) {
          auto metallic         = metallic_roughnessb[{i, j}].z;
          auto roughness        = metallic_roughnessb[{i, j}].y;
          mtexture->ldr[{i, j}] = {metallic, metallic, metallic, 1};
          rtexture->ldr[{i, j}] = {roughness, roughness, roughness, 1};
        }
      }
    }
  }

  // remove empty textures
  for (auto material : scene->materials) {
    if (material->opacity_tex) {
      if (material->opacity_tex->hdr.empty() &&
          material->opacity_tex->ldr.empty())
        material->opacity_tex = nullptr;
    }
  }
  for (auto& texture : scene->textures) {
    if (texture->hdr.empty() && texture->ldr.empty()) {
      delete texture;
      texture = nullptr;
    }
  }
  scene->textures.erase(
      std::remove(scene->textures.begin(), scene->textures.end(), nullptr),
      scene->textures.end());

  // fix scene
  if (scene->name == "") scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // fix cameras
  auto bbox = compute_bounds(scene);
  for (auto camera : scene->cameras) {
    auto center   = (bbox.min + bbox.max) / 2;
    auto distance = dot(-camera->frame.z, center - camera->frame.o);
    if (distance > 0) camera->focus = distance;
  }

  // load done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}  // namespace yocto

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static bool load_pbrt_scene(const string& filename, scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load pbrt
  auto pbrt_guard = std::make_unique<pbrt_model>();
  auto pbrt       = pbrt_guard.get();
  if (!load_pbrt(filename, pbrt, error)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto pcamera : pbrt->cameras) {
    auto camera    = add_camera(scene);
    camera->frame  = pcamera->frame;
    camera->aspect = pcamera->aspect;
    camera->film   = 0.036;
    camera->lens   = pcamera->lens;
    camera->focus  = pcamera->focus;
  }

  // convert materials
  auto ctexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const string& path) -> scene_texture* {
    if (path == "") return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  auto stexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_stexture = [&scene, &stexture_map](
                          const string& path) -> scene_texture* {
    if (path == "") return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
    return texture;
  };
  auto atexture_map = unordered_map<string, scene_texture*>{{"", nullptr}};
  auto get_atexture = [&scene, &atexture_map](
                          const string& path) -> scene_texture* {
    if (path == "") return nullptr;
    auto it = atexture_map.find(path);
    if (it != atexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    atexture_map[path] = texture;
    return texture;
  };

  // convert material
  auto material_map = unordered_map<pbrt_material*, scene_material*>{};
  for (auto pmaterial : pbrt->materials) {
    auto material          = add_material(scene);
    material->emission     = pmaterial->emission;
    material->color        = pmaterial->color;
    material->metallic     = pmaterial->metallic;
    material->specular     = pmaterial->specular;
    material->transmission = pmaterial->transmission;
    material->ior          = pmaterial->ior;
    material->roughness    = pmaterial->roughness;
    material->opacity      = pmaterial->opacity;
    material->thin         = pmaterial->thin;
    material->color_tex    = get_ctexture(pmaterial->color_tex);
    material->opacity_tex  = get_stexture(pmaterial->opacity_tex);
    if (!material->opacity_tex)
      material->opacity_tex = get_atexture(pmaterial->alpha_tex);
    material_map[pmaterial] = material;
  }

  // hack for pbrt empty material
  material_map[nullptr] = add_material(scene);

  // convert shapes
  for (auto pshape : pbrt->shapes) {
    auto shape       = add_shape(scene);
    shape->positions = pshape->positions;
    shape->normals   = pshape->normals;
    shape->texcoords = pshape->texcoords;
    shape->triangles = pshape->triangles;
    for (auto& uv : shape->texcoords) uv.y = 1 - uv.y;
    auto material = material_map.at(pshape->material);
    if (pshape->instances.empty()) {
      auto instance      = add_instance(scene);
      instance->frame    = pshape->frame;
      instance->shape    = shape;
      instance->material = material;
    } else {
      for (auto frame : pshape->instances) {
        auto instance      = add_instance(scene);
        instance->frame    = frame * pshape->frame;
        instance->shape    = shape;
        instance->material = material;
      }
    }
  }

  // convert environments
  for (auto penvironment : pbrt->environments) {
    auto environment          = add_environment(scene);
    environment->frame        = penvironment->frame;
    environment->emission     = penvironment->emission;
    environment->emission_tex = get_ctexture(penvironment->emission_tex);
  }

  // lights
  for (auto plight : pbrt->lights) {
    auto instance                = add_instance(scene);
    instance->shape              = add_shape(scene);
    instance->frame              = plight->area_frame;
    instance->shape->triangles   = plight->area_triangles;
    instance->shape->positions   = plight->area_positions;
    instance->shape->normals     = plight->area_normals;
    instance->material           = add_material(scene);
    instance->material->emission = plight->area_emission;
  }

  // handle progress
  progress.y += (int)scene->textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name) {
    return path_join(path_dirname(filename), name);
  };

  // load texture
  ctexture_map.erase("");
  for (auto [name, texture] : ctexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // load texture
  stexture_map.erase("");
  for (auto [name, texture] : stexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // load alpha
  atexture_map.erase("");
  for (auto [name, texture] : atexture_map) {
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture->hdr, texture->ldr, error))
      return dependent_error();
    for (auto& c : texture->hdr) {
      c = (max(vec3f{c.x, c.y, c.z}) < 0.01) ? vec4f{0, 0, 0, c.w}
                                             : vec4f{1, 1, 1, c.w};
    }
    for (auto& c : texture->ldr) {
      c = (max(vec3i{c.x, c.y, c.z}) < 2) ? vec4b{0, 0, 0, c.w}
                                          : vec4b{255, 255, 255, c.w};
    }
  }

  // fix scene
  if (scene->name == "") scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const scene_model* scene,
    string& error, progress_callback progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt_guard = std::make_unique<pbrt_model>();
  auto pbrt       = pbrt_guard.get();

  // convert camera
  auto camera         = scene->cameras.front();
  auto pcamera        = add_camera(pbrt);
  pcamera->frame      = camera->frame;
  pcamera->lens       = camera->lens;
  pcamera->aspect     = camera->aspect;
  pcamera->resolution = {1280, (int)(1280 / pcamera->aspect)};

  // get texture name
  auto get_texture = [](const scene_texture* texture) {
    return texture ? texture->name : "";
  };

  // convert materials
  auto material_map = unordered_map<scene_material*, pbrt_material*>{};
  for (auto material : scene->materials) {
    auto pmaterial          = add_material(pbrt);
    pmaterial->name         = path_basename(material->name);
    pmaterial->emission     = material->emission;
    pmaterial->color        = material->color;
    pmaterial->metallic     = material->metallic;
    pmaterial->specular     = material->specular;
    pmaterial->transmission = material->transmission;
    pmaterial->roughness    = material->roughness;
    pmaterial->ior          = material->ior;
    pmaterial->opacity      = material->opacity;
    pmaterial->color_tex    = get_texture(material->color_tex);
    pmaterial->opacity_tex  = get_texture(material->opacity_tex);
    material_map[material]  = pmaterial;
  }

  // convert instances
  for (auto instance : scene->instances) {
    auto pshape       = add_shape(pbrt);
    pshape->filename_ = instance->shape->name + ".ply";
    pshape->frame     = instance->frame;
    pshape->frend     = instance->frame;
    pshape->material  = material_map.at(instance->material);
  }

  // convert environments
  for (auto environment : scene->environments) {
    auto penvironment          = add_environment(pbrt);
    penvironment->emission     = environment->emission;
    penvironment->emission_tex = get_texture(environment->emission_tex);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  if (!save_pbrt(filename, pbrt, error)) return false;

  // handle progress
  progress.y += (int)scene->shapes.size() + (int)scene->textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr" : ".png");
    if (!save_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto
