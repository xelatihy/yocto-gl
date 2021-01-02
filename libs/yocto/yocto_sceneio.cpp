//
// Implementation for Yocto/Scene Input and Output functions.
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

#include "yocto_sceneio.h"

#include <algorithm>
#include <cassert>
#include <cctype>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <unordered_map>

#include "ext/cgltf.h"
#include "yocto_color.h"
#include "yocto_commonio.h"
#include "yocto_geometry.h"
#include "yocto_image.h"
#include "yocto_modelio.h"
#include "yocto_parallel.h"
#include "yocto_shading.h"
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
static bool load_json_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_json_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_obj_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_ply_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_stl_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_gltf_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_pbrt_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load a scene
bool load_scene(const string& filename, sceneio_scene* scene, string& error,
    const progress_callback& progress_cb, bool noparallel) {
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
  } else if (ext == ".stl" || ext == ".STL") {
    return load_stl_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    return format_error();
  }
}

// Save a scene
bool save_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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
  } else if (ext == ".gltf" || ext == ".GLTF") {
    return save_gltf_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".ply" || ext == ".PLY") {
    return save_ply_scene(filename, scene, error, progress_cb, noparallel);
  } else if (ext == ".stl" || ext == ".STL") {
    return save_stl_scene(filename, scene, error, progress_cb, noparallel);
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
    return save_ply(filename, &ply, error);
  } else {
    return format_error();
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

#if 0
using njson = nlohmann::json;
using std::array;

// support for json conversions
inline void to_json(njson& j, const vec3f& value) {
  nlohmann::to_json(j, (const array<float, 3>&)value);
}
inline void to_json(njson& j, const vec4f& value) {
  nlohmann::to_json(j, (const array<float, 4>&)value);
}
inline void to_json(njson& j, const frame3f& value) {
  nlohmann::to_json(j, (const array<float, 12>&)value);
}
inline void to_json(njson& j, const mat4f& value) {
  nlohmann::to_json(j, (const array<float, 16>&)value);
}

inline void from_json(const njson& j, vec3f& value) {
  nlohmann::from_json(j, (array<float, 3>&)value);
}
inline void from_json(const njson& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const njson& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}
#endif

// support for json conversions
inline bool set_value(json_tview js, const vec3f& value) {
  return set_value(js, (const array<float, 3>&)value);
}
inline bool set_value(json_tview js, const vec4f& value) {
  return set_value(js, (const array<float, 4>&)value);
}
inline bool set_value(json_tview js, const frame3f& value) {
  return set_value(js, (const array<float, 12>&)value);
}
inline bool set_value(json_tview js, const mat4f& value) {
  return set_value(js, (const array<float, 16>&)value);
}

inline bool get_value(json_ctview js, vec3f& value) {
  return get_value(js, (array<float, 3>&)value);
}
inline bool get_value(json_ctview js, mat3f& value) {
  return get_value(js, (array<float, 9>&)value);
}
inline bool get_value(json_ctview js, frame3f& value) {
  return get_value(js, (array<float, 12>&)value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

// Save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto json_error = [filename]() {
    // error does not need setting
    return false;
  };
  auto parse_error = [filename, &error](string_view message) {
    error = filename + ": parse error (" + string{message} + ")";
    return false;
  };
  auto material_error = [filename, &error](string_view name) {
    error = filename + ": missing material " + string{name};
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
  auto js_tree = json_tree{};
  if (!load_json(filename, js_tree, error)) return json_error();
  auto js = get_croot(js_tree);

  // reference disctionaries
  auto texture_map = unordered_map<string, pair<sceneio_texture*, bool>>{
      {"", {nullptr, true}}};
  auto shape_map = unordered_map<string, pair<sceneio_shape*, bool>>{
      {"", {nullptr, true}}};
  auto material_map = unordered_map<string, pair<sceneio_material*, bool>>{
      {"", {nullptr, true}}};

  // parse json reference
  auto get_shape = [scene, &shape_map](
                       json_ctview js, sceneio_shape*& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = shape_map.find(name);
    if (it != shape_map.end()) {
      value = it->second.first;
      return it->second.first != nullptr;
    }
    auto shape      = add_shape(scene, name);
    shape_map[name] = {shape, false};
    value           = shape;
    return true;
  };

  // parse json reference
  auto get_material = [scene, &material_map](
                          json_ctview js, sceneio_material*& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = material_map.find(name);
    if (it != material_map.end()) {
      value = it->second.first;
      return it->second.first != nullptr;
    }
    auto material      = add_material(scene, name);
    material_map[name] = {material, false};
    value              = material;
    return true;
  };

  // parse json reference
  auto get_texture = [scene, &texture_map](
                         json_ctview js, sceneio_texture*& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = texture_map.find(name);
    if (it != texture_map.end()) {
      value = it->second.first;
      return it->second.first != nullptr;
    }
    auto texture      = add_texture(scene, name);
    texture_map[name] = {texture, false};
    value             = texture;
    return true;
  };

  struct ply_instance {
    vector<frame3f> frames = {};
  };

  // load json instance
  auto ply_instances     = vector<unique_ptr<ply_instance>>{};
  auto ply_instance_map  = unordered_map<string, ply_instance*>{{"", nullptr}};
  auto instance_ply      = unordered_map<sceneio_instance*, ply_instance*>{};
  auto get_ply_instances = [&ply_instances, &ply_instance_map, &instance_ply](
                               json_ctview       js,
                               sceneio_instance* instance) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    if (name.empty()) return true;
    auto it = ply_instance_map.find(name);
    if (it != ply_instance_map.end()) {
      instance_ply[instance] = it->second;
      return true;
    }
    auto ply_instance_ = ply_instances.emplace_back(new ply_instance()).get();
    ply_instance_map[name] = ply_instance_;
    instance_ply[instance] = ply_instance_;
    return true;
  };

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // loop over external dictioaries
  for (auto [gname, group] : iterate_object(js)) {
    if (gname == "asset") {
      for (auto [key, value] : iterate_object(group)) {
        if (key == "copyright") {
          get_value(value, scene->copyright);
        } else if (key == "generator") {
          auto generator = string{};
          get_value(value, generator);
        } else {
          set_error(group, "unknown key " + string{key});
        }
      }
    } else if (gname == "cameras") {
      for (auto [name, element] : iterate_object(group)) {
        auto camera = add_camera(scene, string{name});
        for (auto [key, value] : iterate_object(element)) {
          if (key == "frame") {
            get_value(value, camera->frame);
          } else if (key == "orthographic") {
            get_value(value, camera->orthographic);
          } else if (key == "ortho") {
            // backward compatibility
            get_value(value, camera->orthographic);
          } else if (key == "lens") {
            get_value(value, camera->lens);
          } else if (key == "aspect") {
            get_value(value, camera->aspect);
          } else if (key == "film") {
            get_value(value, camera->film);
          } else if (key == "focus") {
            get_value(value, camera->focus);
          } else if (key == "aperture") {
            get_value(value, camera->aperture);
          } else if (key == "lookat") {
            get_value(value, (mat3f&)camera->frame);
            camera->focus = length(camera->frame.x - camera->frame.y);
            camera->frame = lookat_frame(
                camera->frame.x, camera->frame.y, camera->frame.z);
          } else {
            set_error(js, "unknown key " + string{key});
          }
        }
      }
    } else if (gname == "environments") {
      for (auto [name, element] : iterate_object(group)) {
        auto environment = add_environment(scene, string{name});
        for (auto [key, value] : iterate_object(element)) {
          if (key == "frame") {
            get_value(value, environment->frame);
          } else if (key == "emission") {
            get_value(value, environment->emission);
          } else if (key == "emission_tex") {
            get_texture(value, environment->emission_tex);
          } else if (key == "lookat") {
            get_value(value, (mat3f&)environment->frame);
            environment->frame = lookat_frame(environment->frame.x,
                environment->frame.y, environment->frame.z, true);
          } else {
            set_error(js, "unknown key " + string{key});
          }
        }
      }
    } else if (gname == "materials") {
      for (auto [name, element] : iterate_object(group)) {
        auto material_it = material_map.find(string{name});
        auto material    = (material_it == material_map.end())
                               ? add_material(scene, string{name})
                               : material_it->second.first;
        for (auto [key, value] : iterate_object(element)) {
          if (key == "emission") {
            get_value(value, material->emission);
          } else if (key == "color") {
            get_value(value, material->color);
          } else if (key == "metallic") {
            get_value(value, material->metallic);
          } else if (key == "specular") {
            get_value(value, material->specular);
          } else if (key == "roughness") {
            get_value(value, material->roughness);
          } else if (key == "coat") {
            get_value(value, material->coat);
          } else if (key == "transmission") {
            get_value(value, material->transmission);
          } else if (key == "translucency") {
            get_value(value, material->translucency);
          } else if (key == "thin") {
            get_value(value, material->thin);
          } else if (key == "ior") {
            get_value(value, material->ior);
          } else if (key == "trdepth") {
            get_value(value, material->trdepth);
          } else if (key == "scattering") {
            get_value(value, material->scattering);
          } else if (key == "scanisotropy") {
            get_value(value, material->scanisotropy);
          } else if (key == "opacity") {
            get_value(value, material->opacity);
          } else if (key == "coat") {
            get_value(value, material->coat);
          } else if (key == "emission_tex") {
            get_texture(value, material->emission_tex);
          } else if (key == "color_tex") {
            get_texture(value, material->color_tex);
          } else if (key == "metallic_tex") {
            get_texture(value, material->metallic_tex);
          } else if (key == "specular_tex") {
            get_texture(value, material->specular_tex);
          } else if (key == "transmission_tex") {
            get_texture(value, material->transmission_tex);
          } else if (key == "translucency_tex") {
            get_texture(value, material->translucency_tex);
          } else if (key == "roughness_tex") {
            get_texture(value, material->roughness_tex);
          } else if (key == "scattering_tex") {
            get_texture(value, material->scattering_tex);
          } else if (key == "opacity_tex") {
            get_texture(value, material->opacity_tex);
          } else if (key == "normal_tex") {
            get_texture(value, material->normal_tex);
          } else {
            set_error(element, "unknown key " + string{key});
          }
        }
        material_map[material->name] = {material, true};
      }
    } else if (gname == "instances" || gname == "objects") {
      for (auto [name, element] : iterate_object(group)) {
        auto instance = add_instance(scene, string{name});
        for (auto [key, value] : iterate_object(element)) {
          if (key == "frame") {
            get_value(value, instance->frame);
          } else if (key == "lookat") {
            get_value(value, (mat3f&)instance->frame);
            instance->frame = lookat_frame(
                instance->frame.x, instance->frame.y, instance->frame.z, true);
          } else if (key == "material") {
            get_material(value, instance->material);
          } else if (key == "shape") {
            get_shape(value, instance->shape);
          } else if (key == "instance") {
            get_ply_instances(value, instance);
          } else if (key == "subdivisions") {
            if (instance->shape) {
              get_value(value, instance->shape->subdivisions);
            }
          } else if (key == "catmullcark") {
            if (instance->shape) {
              get_value(value, instance->shape->catmullclark);
            }
          } else if (key == "smooth") {
            if (instance->shape) {
              get_value(value, instance->shape->smooth);
            }
          } else if (key == "displacement") {
            if (instance->shape) {
              get_value(value, instance->shape->displacement);
            }
          } else if (key == "displacement_tex") {
            if (instance->shape) {
              get_texture(value, instance->shape->displacement_tex);
            }
          } else if (key == "instance") {
            get_ply_instances(value, instance);
          } else {
            // set_error(element, "unknown key " + string{key});
          }
        }
      }
    } else {
      set_error(js, "unknown key " + string{gname});
    }
  }

  // check for parsing errors
  if (!is_valid(js)) return parse_error(get_error(js));

  // check materials
  for (auto& [key, value] : material_map) {
    if (!value.second) return material_error(key);
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
  for (auto [name, value] : shape_map) {
    auto shape = value.first;
    if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
    auto path = make_filename(name, "shapes", {".ply", ".obj"});
    if (!load_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->quadspos, shape->quadsnorm,
            shape->quadstexcoord, shape->positions, shape->normals,
            shape->texcoords, shape->colors, shape->radius, error,
            shape->catmullclark && shape->subdivisions > 0))
      return dependent_error();
  }
  // load textures
  texture_map.erase("");
  for (auto [name, value] : texture_map) {
    auto texture = value.first;
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
  if (scene->name.empty()) scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);
  trim_memory(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto conversion_error = [filename, &error](const string& message) {
    // should never happen
    throw std::runtime_error{"programmer error"};
    error = filename + ": conversion error (" + message + ")";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{
      0, 2 + (int)scene->shapes.size() + (int)scene->textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json file
  auto js_tree = json_tree{};
  auto js      = get_root(js_tree);
  set_object(js);

  // asset
  {
    auto element = insert_object(js, "asset");
    insert_value(element, "generator",
        "Yocto/GL - https://github.com/xelatihy/yocto-gl");
    if (!scene->copyright.empty()) {
      insert_value(element, "copyright", scene->copyright);
    }
  }

  auto def_cam = sceneio_camera{};
  if (!scene->cameras.empty()) {
    // auto _ = append_object(js, "cameras");
    auto group = insert_object(js, "cameras");
    for (auto& camera : scene->cameras) {
      auto elemnt = insert_object(group, camera->name);
      if (camera->frame != def_cam.frame) {
        insert_value(elemnt, "frame", camera->frame);
      }
      if (camera->orthographic != def_cam.orthographic) {
        insert_value(elemnt, "orthographic", camera->orthographic);
      }
      if (camera->lens != def_cam.lens) {
        insert_value(elemnt, "lens", camera->lens);
      }
      if (camera->aspect != def_cam.aspect) {
        insert_value(elemnt, "aspect", camera->aspect);
      }
      if (camera->film != def_cam.film) {
        insert_value(elemnt, "film", camera->film);
      }
      if (camera->focus != def_cam.focus) {
        insert_value(elemnt, "focus", camera->focus);
      }
      if (camera->aperture != def_cam.aperture) {
        insert_value(elemnt, "aperture", camera->aperture);
      }
    }
  }

  auto def_env = sceneio_environment{};
  if (!scene->environments.empty()) {
    auto group = insert_object(js, "environments");
    for (auto environment : scene->environments) {
      auto elemnt = insert_object(group, environment->name);
      if (environment->frame != def_env.frame) {
        insert_value(elemnt, "frame", environment->frame);
      }
      if (environment->emission != def_env.emission) {
        insert_value(elemnt, "emission", environment->emission);
      }
      if (environment->emission_tex != nullptr) {
        insert_value(elemnt, "emission_tex", environment->emission_tex->name);
      }
    }
  }

  auto def_material = sceneio_material{};
  if (!scene->materials.empty()) {
    auto group = insert_object(js, "materials");
    for (auto material : scene->materials) {
      auto elment = insert_object(group, material->name);
      if (material->emission != def_material.emission) {
        insert_value(elment, "emission", material->emission);
      }
      if (material->color != def_material.color) {
        insert_value(elment, "color", material->color);
      }
      if (material->specular != def_material.specular) {
        insert_value(elment, "specular", material->specular);
      }
      if (material->metallic != def_material.metallic) {
        insert_value(elment, "metallic", material->metallic);
      }
      if (material->coat != def_material.coat) {
        insert_value(elment, "coat", material->coat);
      }
      if (material->roughness != def_material.roughness) {
        insert_value(elment, "roughness", material->roughness);
      }
      if (material->ior != def_material.ior) {
        insert_value(elment, "ior", material->ior);
      }
      if (material->transmission != def_material.transmission) {
        insert_value(elment, "transmission", material->transmission);
      }
      if (material->translucency != def_material.translucency) {
        insert_value(elment, "translucency", material->translucency);
      }
      if (material->trdepth != def_material.trdepth) {
        insert_value(elment, "trdepth", material->trdepth);
      }
      if (material->scattering != def_material.scattering) {
        insert_value(elment, "scattering", material->scattering);
      }
      if (material->scanisotropy != def_material.scanisotropy) {
        insert_value(elment, "scanisotropy", material->scanisotropy);
      }
      if (material->opacity != def_material.opacity) {
        insert_value(elment, "opacity", material->opacity);
      }
      if (material->thin != def_material.thin) {
        insert_value(elment, "thin", material->thin);
      }
      if (material->emission_tex != nullptr) {
        insert_value(elment, "emission_tex", material->emission_tex->name);
      }
      if (material->color_tex != nullptr) {
        insert_value(elment, "color_tex", material->color_tex->name);
      }
      if (material->metallic_tex != nullptr) {
        insert_value(elment, "metallic_tex", material->metallic_tex->name);
      }
      if (material->specular_tex != nullptr) {
        insert_value(elment, "specular_tex", material->specular_tex->name);
      }
      if (material->roughness_tex != nullptr) {
        insert_value(elment, "roughness_tex", material->roughness_tex->name);
      }
      if (material->transmission_tex != nullptr) {
        insert_value(
            elment, "transmission_tex", material->transmission_tex->name);
      }
      if (material->translucency_tex != nullptr) {
        insert_value(
            elment, "translucency_tex", material->translucency_tex->name);
      }
      if (material->scattering_tex != nullptr) {
        insert_value(elment, "scattering_tex", material->scattering_tex->name);
      }
      if (material->coat_tex != nullptr) {
        insert_value(elment, "coat_tex", material->coat_tex->name);
      }
      if (material->opacity_tex != nullptr) {
        insert_value(elment, "opacity_tex", material->opacity_tex->name);
      }
      if (material->normal_tex != nullptr) {
        insert_value(elment, "normal_tex", material->normal_tex->name);
      }
    }
  }

  auto def_instance = sceneio_instance{};
  auto def_shape    = sceneio_shape{};
  if (!scene->instances.empty()) {
    auto group = insert_object(js, "instances");
    for (auto instance : scene->instances) {
      auto elment = insert_object(group, instance->name);
      if (instance->frame != def_instance.frame) {
        insert_value(elment, "frame", instance->frame);
      }
      if (instance->shape != nullptr) {
        insert_value(elment, "shape", instance->shape->name);
      }
      if (instance->material != nullptr) {
        insert_value(elment, "material", instance->material->name);
      }
      if (instance->shape != nullptr) {
        if (instance->shape->subdivisions != def_shape.subdivisions) {
          insert_value(elment, "subdivisions", instance->shape->subdivisions);
        }
        if (instance->shape->catmullclark != def_shape.catmullclark) {
          insert_value(elment, "catmullclark", instance->shape->catmullclark);
        }
        if (instance->shape->smooth != def_shape.smooth) {
          insert_value(elment, "smooth", instance->shape->smooth);
        }
        if (instance->shape->displacement != def_shape.displacement) {
          insert_value(elment, "displacement", instance->shape->displacement);
        }
        if (instance->shape->displacement_tex != nullptr) {
          insert_value(elment, "displacement_tex",
              instance->shape->displacement_tex->name);
        }
      }
    }
  }

  // chck for errors
  if (!is_valid(js)) return conversion_error(get_error(js));

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json
  if (!save_json(filename, js_tree, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // save shapes
  for (auto shape : scene->shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    auto path = make_filename(shape->name, "shapes",
        (shape->catmullclark && shape->subdivisions > 0) ? ".obj"s : ".ply"s);
    if (!save_shape(path, shape->points, shape->lines, shape->triangles,
            shape->quads, shape->quadspos, shape->quadsnorm,
            shape->quadstexcoord, shape->positions, shape->normals,
            shape->texcoords, shape->colors, shape->radius, error,
            shape->catmullclark && shape->subdivisions > 0))
      return dependent_error();
  }

  // save textures
  for (auto texture : scene->textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr"s : ".png"s);
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
static bool load_obj_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
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

  // load obj
  auto obj_guard = std::make_unique<obj_scene>();
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&ctexture_map, scene](
                          const obj_texture& tinfo) -> sceneio_texture* {
    auto path = tinfo.path;
    if (path.empty()) return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };

  // helper to create texture maps
  auto stexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_stexture = [&stexture_map, scene](
                          const obj_texture& tinfo) -> sceneio_texture* {
    auto path = tinfo.path;
    if (path.empty()) return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
    return texture;
  };

  // handler for materials
  auto material_map = unordered_map<string, sceneio_material*>{};
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
    material_map[omat->name]   = material;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto oshape : obj->shapes) {
    auto& materials = oshape->materials;
    if (materials.empty()) materials.push_back(nullptr);
    for (auto material_idx = 0; material_idx < materials.size();
         material_idx++) {
      auto shape = add_shape(scene);
      if (material_map.find(materials[material_idx]) == material_map.end())
        return material_error(materials[material_idx]);
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
  if (scene->name.empty()) scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_obj_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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

  auto obj_guard = std::make_unique<obj_scene>();
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
  auto get_texture = [](sceneio_texture* texture) {
    if (texture == nullptr) return obj_texture{};
    auto tinfo = obj_texture{};
    tinfo.path = "textures/" + texture->name +
                 (!texture->hdr.empty() ? ".hdr"s : ".png"s);
    return tinfo;
  };

  // convert materials and textures
  auto material_map = unordered_map<sceneio_material*, string>{{nullptr, ""s}};
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
    material_map[material]          = omaterial->name;
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
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr"s : ".png"s);
    if (!save_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

void print_obj_camera(sceneio_camera* camera) {
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

static bool load_ply_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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

static bool save_ply_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_stl_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load stl mesh
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

static bool save_stl_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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
static bool load_gltf_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
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
  if (!dirname.empty()) dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success)
    return read_error();

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert asset
  {
    auto gast = &gltf->asset;
    if (gast->copyright != nullptr) scene->copyright = gast->copyright;
  }

  // prepare list of effective nodes
  auto visible_nodes = vector<bool>(gltf->nodes_count, false);
  auto gscene        = gltf->scene != nullptr ? gltf->scene : gltf->scenes;
  if (gscene != nullptr) {
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
    if (gnde->camera == nullptr) continue;
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const cgltf_texture_view& ginfo) -> sceneio_texture* {
    if (ginfo.texture == nullptr || ginfo.texture->image == nullptr)
      return nullptr;
    auto path = string{ginfo.texture->image->uri};
    if (path.empty()) return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  // convert color opacity textures
  auto cotexture_map =
      unordered_map<string, pair<sceneio_texture*, sceneio_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_cotexture = [&scene, &cotexture_map](const cgltf_texture_view& ginfo)
      -> pair<sceneio_texture*, sceneio_texture*> {
    if (ginfo.texture == nullptr || ginfo.texture->image == nullptr)
      return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path.empty()) return {nullptr, nullptr};
    auto it = cotexture_map.find(path);
    if (it != cotexture_map.end()) return it->second;
    auto color_texture   = add_texture(scene);
    auto opacity_texture = add_texture(scene);
    cotexture_map[path]  = {color_texture, opacity_texture};
    return {color_texture, opacity_texture};
  };
  // convert textures
  auto mrtexture_map =
      unordered_map<string, pair<sceneio_texture*, sceneio_texture*>>{
          {"", {nullptr, nullptr}}};
  auto get_mrtexture = [&scene, &mrtexture_map](const cgltf_texture_view& ginfo)
      -> pair<sceneio_texture*, sceneio_texture*> {
    if (ginfo.texture == nullptr || ginfo.texture->image == nullptr)
      return {nullptr, nullptr};
    auto path = string{ginfo.texture->image->uri};
    if (path.empty()) return {nullptr, nullptr};
    auto it = mrtexture_map.find(path);
    if (it != mrtexture_map.end()) return it->second;
    auto metallic_texture  = add_texture(scene);
    auto roughness_texture = add_texture(scene);
    mrtexture_map[path]    = {metallic_texture, roughness_texture};
    return {metallic_texture, roughness_texture};
  };

  // convert materials
  auto material_map = unordered_map<cgltf_material*, sceneio_material*>{
      {nullptr, nullptr}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto gmaterial         = &gltf->materials[mid];
    auto material          = add_material(scene);
    material->emission     = {gmaterial->emissive_factor[0],
        gmaterial->emissive_factor[1], gmaterial->emissive_factor[2]};
    material->emission_tex = get_ctexture(gmaterial->emissive_texture);
    if (gmaterial->has_pbr_metallic_roughness != 0) {
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
  auto mesh_map = unordered_map<cgltf_mesh*, vector<sceneio_instance*>>{
      {nullptr, {}}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh = &gltf->meshes[mid];
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (gprim->attributes_count == 0) continue;
      auto instance = add_instance(scene);
      mesh_map[gmesh].push_back(instance);
      auto shape         = add_shape(scene);
      instance->shape    = shape;
      instance->material = material_map.at(gprim->material);
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name != nullptr ? gattr->name : "");
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
        } else if (semantic == "_RADIUS") {
          shape->radius.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape->radius[i], 1);
        } else {
          // ignore
        }
      }
      // indices
      if (gprim->indices == nullptr) {
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
    if (gnde->mesh == nullptr) continue;
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
    if (material->opacity_tex != nullptr) {
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
  if (scene->name.empty()) scene->name = path_basename(filename);
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
}

// Load a scene
static bool save_gltf_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };
  auto fvshape_error = [filename, &error]() {
    error = filename + ": face-varying not supported";
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 3};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // setup json
  using json = json_value;

  // convert scene to json
  auto js = json{};
  js      = json::object();

  // asset
  {
    auto& ajs      = js["asset"];
    ajs            = json::object();
    ajs["version"] = "2.0";
    ajs["generator"] =
        "Saved with Yocto/GL --- https://github.com/xelatihy/yocto-gl";
    ajs["copyright"] = scene->copyright;
  }

  // cameras
  if (!scene->cameras.empty()) {
    auto& ajs = js["cameras"];
    ajs       = json::array();
    for (auto camera : scene->cameras) {
      auto& cjs          = ajs.emplace_back();
      cjs                = json::object();
      cjs["name"]        = camera->name;
      cjs["type"]        = "perspective";
      auto& pjs          = cjs["perspective"];
      pjs                = json::object();
      pjs["aspectRatio"] = camera->aspect;
      pjs["yfov"]        = 0.660593;  // TODO(fabio): yfov
      pjs["znear"]       = 0.001;     // TODO(fabio): configurable?
    }
  }

  // materials
  auto textures     = vector<pair<string, image<vec4b>>>{};
  auto texture_map  = unordered_map<string, int>{};
  auto material_map = unordered_map<const sceneio_material*, int>{};
  if (!scene->materials.empty()) {
    auto& ajs = js["materials"];
    ajs       = json::array();
    for (auto material : scene->materials) {
      auto& mjs              = ajs.emplace_back();
      mjs                    = json::object();
      mjs["name"]            = material->name;
      mjs["emissiveFactor"]  = to_json(array<float, 3>{
          material->emission.x, material->emission.y, material->emission.z});
      auto& pjs              = mjs["pbrMetallicRoughness"];
      pjs                    = json::object();
      pjs["baseColorFactor"] = to_json(array<float, 4>{material->color.x,
          material->color.y, material->color.z, material->opacity});
      pjs["metallicFactor"]  = material->metallic;
      pjs["roughnessFactor"] = material->roughness;
      if (material->emission_tex) {
        auto tname = material->emission_tex->name;  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          textures.emplace_back(tname, material->emission_tex->ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        mjs["emissiveTexture"]          = json::object();
        mjs["emissiveTexture"]["index"] = texture_map.at(tname);
      }
      if (material->normal_tex) {
        auto tname = material->normal_tex->name;  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          textures.emplace_back(tname, material->normal_tex->ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        mjs["normalTexture"]          = json::object();
        mjs["normalTexture"]["index"] = texture_map.at(tname);
      }
      if (material->color_tex) {                 // TODO(fabio): opacity
        auto tname = material->color_tex->name;  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          textures.emplace_back(tname, material->color_tex->ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        pjs["baseColorTexture"]          = json::object();
        pjs["baseColorTexture"]["index"] = texture_map.at(tname);
      }
      if (material->roughness_tex) {                 // TODO(fabio): roughness
        auto tname = material->roughness_tex->name;  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          textures.emplace_back(tname, material->roughness_tex->ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        pjs["metallicRoughnessTexture"]          = json::object();
        pjs["metallicRoughnessTexture"]["index"] = texture_map.at(tname);
      }
      material_map[material] = (int)js["materials"].size() - 1;
    }
  }

  // textures
  if (!textures.empty()) {
    js["textures"] = json::array();
    js["samplers"] = json::array();
    js["images"]   = json::array();
    auto& sjs      = js["samplers"].emplace_back();
    sjs            = json::object();
    sjs["name"]    = "sampler";
    for (auto& [name, img] : textures) {
      auto& ijs      = js["images"].emplace_back();
      ijs            = json::object();
      ijs["name"]    = name;
      ijs["uri"]     = "textures/" + name + ".png";
      auto& tjs      = js["textures"].emplace_back();
      tjs            = json::object();
      tjs["name"]    = name;
      tjs["sampler"] = 0;
      tjs["source"]  = (int)js["images"].size() - 1;
    }
  }

  // add an accessor
  auto add_accessor = [](json& js, vector<pair<string, vector<byte>>>& buffers,
                          const void* data, size_t count, size_t size,
                          bool is_index = false) -> int {
    static auto types = unordered_map<size_t, string>{
        {1, "SCALAR"}, {2, "VEC2"}, {3, "VEC3"}, {4, "VEC4"}};
    auto  length         = count * size * 4;
    auto& vjs            = js["bufferViews"].emplace_back();
    vjs                  = json::object();
    vjs["buffer"]        = (int)buffers.size() - 1;
    vjs["byteLength"]    = (uint64_t)length;
    vjs["byteOffset"]    = (uint64_t)buffers.back().second.size();
    vjs["target"]        = is_index ? 34963 : 34962;
    auto& ajs            = js["accessors"].emplace_back();
    ajs                  = json::object();
    ajs["bufferView"]    = (int)js["bufferViews"].size() - 1;
    ajs["byteOffset"]    = 0;
    ajs["componentType"] = is_index ? 5125 : 5126;
    ajs["count"]         = (uint64_t)count;
    ajs["type"]          = types.at(size);
    if (!is_index) {
      auto min_ = vector<float>(size, flt_max);
      auto max_ = vector<float>(size, flt_min);
      for (auto idx = (size_t)0; idx < count; idx++) {
        for (auto channel = (size_t)0; channel < size; channel++) {
          auto value    = (float*)data + idx * size + channel;
          min_[channel] = min(min_[channel], *value);
          max_[channel] = max(max_[channel], *value);
        }
      }
      ajs["min"] = to_json(min_);
      ajs["max"] = to_json(max_);
    }
    buffers.back().second.insert(
        buffers.back().second.end(), (byte*)data, (byte*)data + length);
    return (int)js["accessors"].size() - 1;
  };

  // meshes
  auto buffers        = vector<pair<string, vector<byte>>>{};
  auto primitives_map = unordered_map<const sceneio_shape*, json>{};
  if (!scene->shapes.empty()) {
    js["accessors"]   = json::array();
    js["bufferViews"] = json::array();
    js["buffers"]     = json::array();
    for (auto shape : scene->shapes) {
      auto& buffer = buffers.emplace_back(shape->name, vector<byte>{}).second;
      auto& pjs    = primitives_map[shape];
      pjs          = json::object();
      auto& ajs    = pjs["attributes"];
      ajs          = json::object();
      if (!shape->positions.empty()) {
        ajs["POSITION"] = add_accessor(
            js, buffers, shape->positions.data(), shape->positions.size(), 3);
      }
      if (!shape->normals.empty()) {
        ajs["NORMAL"] = add_accessor(
            js, buffers, shape->normals.data(), shape->normals.size(), 3);
      }
      if (!shape->texcoords.empty()) {
        ajs["TEXCOORD_0"] = add_accessor(
            js, buffers, shape->texcoords.data(), shape->texcoords.size(), 2);
      }
      if (!shape->colors.empty()) {
        ajs["COLOR_0"] = add_accessor(
            js, buffers, shape->colors.data(), shape->colors.size(), 4);
      }
      if (!shape->radius.empty()) {
        ajs["_RADIUS"] = add_accessor(
            js, buffers, shape->radius.data(), shape->radius.size(), 1);
      }
      if (!shape->points.empty()) {
        pjs["indices"] = add_accessor(
            js, buffers, shape->points.data(), shape->points.size(), 1, true);
        pjs["mode"] = 0;
      } else if (!shape->lines.empty()) {
        pjs["indices"] = add_accessor(
            js, buffers, shape->lines.data(), shape->lines.size() * 2, 1, true);
        pjs["mode"] = 1;
      } else if (!shape->triangles.empty()) {
        pjs["indices"] = add_accessor(js, buffers, shape->triangles.data(),
            shape->triangles.size() * 3, 1, true);
        pjs["mode"]    = 4;
      } else if (!shape->quads.empty()) {
        auto triangles = quads_to_triangles(shape->quads);
        pjs["indices"] = add_accessor(
            js, buffers, triangles.data(), triangles.size() * 3, 1, true);
        pjs["mode"] = 4;
      } else if (!shape->quadspos.empty()) {
        return fvshape_error();
      }
      auto& bjs         = js["buffers"].emplace_back();
      bjs               = json::object();
      bjs["byteLength"] = (uint64_t)buffer.size();
      bjs["uri"]        = "shapes/" + shape->name + ".bin";
    }
  }

  // nodes
  js["nodes"] = json::array();
  if (!scene->cameras.empty()) {
    auto camera_id = 0;
    for (auto camera : scene->cameras) {
      auto& njs   = js["nodes"].emplace_back();
      njs         = json::object();
      njs["name"] = camera->name;
      auto matrix = frame_to_mat(camera->frame);  // TODO(fabio): do this better
      njs["matrix"] = to_json((array<float, 16>&)matrix);
      njs["camera"] = camera_id++;
    }
  }
  if (!scene->instances.empty()) {
    js["meshes"]   = json::array();
    using mesh_key = pair<const sceneio_shape*, const sceneio_material*>;
    struct mesh_key_hash {
      size_t operator()(const mesh_key& v) const {
        const std::hash<const void*> hasher = std::hash<const void*>();
        auto                         h      = (size_t)0;
        h ^= hasher(v.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hasher(v.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
      }
    };
    auto mesh_map = unordered_map<mesh_key, int, mesh_key_hash>{};
    for (auto instance : scene->instances) {
      auto& njs   = js["nodes"].emplace_back();
      njs         = json::object();
      njs["name"] = instance->name;
      auto matrix = frame_to_mat(
          instance->frame);  // TODO(fabio): do this better
      njs["matrix"] = to_json((array<float, 16>&)matrix);
      if (mesh_map.find(mesh_key{instance->shape, instance->material}) ==
          mesh_map.end()) {
        auto& mjs   = js["meshes"].emplace_back();
        mjs         = json::object();
        mjs["name"] = instance->shape->name + "_" + instance->material->name;
        mjs["primitives"] = json::array();
        mjs["primitives"].push_back(primitives_map.at(instance->shape));
        mjs["primitives"].back()["material"] = material_map.at(
            instance->material);
        mesh_map[mesh_key{instance->shape, instance->material}] =
            (int)js["meshes"].size() - 1;
      }
      njs["mesh"] = mesh_map.at({instance->shape, instance->material});
    }
  } else {
    js["meshes"] = json::array();
    for (auto& [shape, pjs] : primitives_map) {
      auto& mjs         = js["meshes"].emplace_back();
      mjs               = json::object();
      mjs["name"]       = shape->name;
      mjs["primitives"] = json::array();
      mjs["primitives"].push_back(pjs);
      auto& njs   = js["nodes"].emplace_back();
      njs         = json::object();
      njs["name"] = shape->name;
      njs["mesh"] = (int)js["meshes"].size() - 1;
    }
  }

  // root children
  {
    auto& rjs       = js["nodes"].emplace_back();
    rjs             = json::object();
    rjs["name"]     = "root";
    rjs["children"] = json::array();
    for (auto idx = 0; idx < (int)js["nodes"].size() - 1; idx++)
      rjs["children"].push_back(json{idx});
  }

  // scene
  {
    js["scenes"] = json::array();
    auto& sjs    = js["scenes"].emplace_back();
    sjs          = json::object();
    sjs["nodes"] = json::array();
    sjs["nodes"].push_back(json{(int)js["nodes"].size() - 1});
    js["scene"] = 0;
  }

  // save json
  if (!save_json(filename, js, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // dirname
  auto dirname = path_dirname(filename);

  // save shapes
  for (auto& [name, buffer] : buffers) {
    if (progress_cb) progress_cb("save buffer", progress.x++, progress.y);
    if (!save_binary(
            path_join(dirname, "shapes/" + name + ".bin"), buffer, error))
      return dependent_error();
  }

  // save textures
  for (auto& [name, texture] : textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    if (!save_image(
            path_join(dirname, "textures/" + name + ".png"), texture, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static bool load_pbrt_scene(const string& filename, sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load pbrt
  auto pbrt_guard = std::make_unique<pbrt_scene>();
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
  auto ctexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_ctexture = [&scene, &ctexture_map](
                          const string& path) -> sceneio_texture* {
    if (path.empty()) return nullptr;
    auto it = ctexture_map.find(path);
    if (it != ctexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    ctexture_map[path] = texture;
    return texture;
  };
  auto stexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_stexture = [&scene, &stexture_map](
                          const string& path) -> sceneio_texture* {
    if (path.empty()) return nullptr;
    auto it = stexture_map.find(path);
    if (it != stexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    stexture_map[path] = texture;
    return texture;
  };
  auto atexture_map = unordered_map<string, sceneio_texture*>{{"", nullptr}};
  auto get_atexture = [&scene, &atexture_map](
                          const string& path) -> sceneio_texture* {
    if (path.empty()) return nullptr;
    auto it = atexture_map.find(path);
    if (it != atexture_map.end()) return it->second;
    auto texture       = add_texture(scene);
    atexture_map[path] = texture;
    return texture;
  };

  // convert material
  auto material_map = unordered_map<string, sceneio_material*>{};
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
    if (material->opacity_tex == nullptr)
      material->opacity_tex = get_atexture(pmaterial->alpha_tex);
    material_map[pmaterial->name] = material;
  }

  // hack for pbrt empty material
  material_map[""] = add_material(scene);

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
  if (scene->name.empty()) scene->name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const sceneio_scene* scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt_guard = std::make_unique<pbrt_scene>();
  auto pbrt       = pbrt_guard.get();

  // convert camera
  auto camera         = scene->cameras.front();
  auto pcamera        = add_camera(pbrt);
  pcamera->frame      = camera->frame;
  pcamera->lens       = camera->lens;
  pcamera->aspect     = camera->aspect;
  pcamera->resolution = {1280, (int)(1280 / pcamera->aspect)};

  // get texture name
  auto get_texture = [](const sceneio_texture* texture) {
    return texture != nullptr ? texture->name : "";
  };

  // convert materials
  auto material_map = unordered_map<sceneio_material*, string>{};
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
    material_map[material]  = pmaterial->name;
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
        texture->name, "textures", (!texture->hdr.empty()) ? ".hdr"s : ".png"s);
    if (!save_image(path, texture->hdr, texture->ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto
