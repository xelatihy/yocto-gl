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
#include "ext/cgltf_write.h"
#include "ext/fast_obj.h"
#include "ext/json.hpp"
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
// UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Enumerate
template <typename T>
auto enumerate(const vector<T>& elements) {
  struct item {
    const size_t& idx;
    const T&      element;
  };
  struct iterator {
    size_t    idx;
    const T*  element;
    bool      operator!=(const iterator& other) { return idx != other.idx; }
    iterator& operator++() {
      ++element;
      ++idx;
      return *this;
    }
    item operator*() { return {idx, *element}; }
  };
  struct wrapper {
    const vector<T>& elements;
    iterator         begin() { return {0, elements.data()}; }
    iterator         end() {
      return {elements.size(), elements.data() + elements.size()};
    }
  };
  return wrapper{elements};
}

template <typename T1, typename T2>
auto zip(const vector<T1>& keys, const vector<T2>& values) {
  struct item {
    const T1& key;
    const T2& value;
  };
  struct iterator {
    const T1* key;
    const T2* value;
    bool      operator!=(const iterator& other) {
      return key != other.key || value != other.value;
    }
    void operator++() {
      ++key;
      ++value;
    }
    item operator*() { return item{*key, *value}; }
  };
  struct wrapper {
    const vector<T1>& keys;
    const vector<T2>& values;
    iterator          begin() { return {keys.data(), values.data()}; }
    iterator          end() {
      return {keys.data() + keys.size(), values.data() + values.size()};
    }
  };
  return wrapper{keys, values};
}

// get name
[[maybe_unused]] static string get_camera_name(
    const scene_scene& scene, int idx) {
  if (scene.camera_names.empty()) return "camera" + std::to_string(idx);
  return scene.camera_names[idx];
}
[[maybe_unused]] static string get_environment_name(
    const scene_scene& scene, int idx) {
  if (scene.environment_names.empty())
    return "environment" + std::to_string(idx);
  return scene.environment_names[idx];
}
[[maybe_unused]] static string get_shape_name(
    const scene_scene& scene, int idx) {
  if (scene.shape_names.empty()) return "shape" + std::to_string(idx);
  return scene.shape_names[idx];
}
[[maybe_unused]] static string get_texture_name(
    const scene_scene& scene, int idx) {
  if (scene.texture_names.empty()) return "texture" + std::to_string(idx);
  return scene.texture_names[idx];
}
[[maybe_unused]] static string get_instance_name(
    const scene_scene& scene, int idx) {
  if (scene.instance_names.empty()) return "instance" + std::to_string(idx);
  return scene.instance_names[idx];
}
[[maybe_unused]] static string get_material_name(
    const scene_scene& scene, int idx) {
  if (scene.material_names.empty()) return "material" + std::to_string(idx);
  return scene.material_names[idx];
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_scene& scene, int idx) {
  if (scene.subdiv_names.empty()) return "subdiv" + std::to_string(idx);
  return scene.subdiv_names[idx];
}

[[maybe_unused]] static string get_camera_name(
    const scene_scene& scene, const scene_camera& camera) {
  return get_camera_name(scene, (int)(&camera - scene.cameras.data()));
}
[[maybe_unused]] static string get_environment_name(
    const scene_scene& scene, const scene_environment& environment) {
  return get_environment_name(
      scene, (int)(&environment - scene.environments.data()));
}
[[maybe_unused]] static string get_shape_name(
    const scene_scene& scene, const scene_shape& shape) {
  return get_shape_name(scene, (int)(&shape - scene.shapes.data()));
}
[[maybe_unused]] static string get_texture_name(
    const scene_scene& scene, const scene_texture& texture) {
  return get_texture_name(scene, (int)(&texture - scene.textures.data()));
}
[[maybe_unused]] static string get_instance_name(
    const scene_scene& scene, const scene_instance& instance) {
  return get_instance_name(scene, (int)(&instance - scene.instances.data()));
}
[[maybe_unused]] static string get_material_name(
    const scene_scene& scene, const scene_material& material) {
  return get_material_name(scene, (int)(&material - scene.materials.data()));
}
[[maybe_unused]] static string get_subdiv_name(
    const scene_scene& scene, const scene_subdiv& subdiv) {
  return get_subdiv_name(scene, (int)(&subdiv - scene.subdivs.data()));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TEST SCENES
// -----------------------------------------------------------------------------
namespace yocto {

enum struct test_cameras_type { standard, wide };
enum struct test_environments_type { none, sky, sunsky };
enum struct test_arealights_type { none, standard, large };
enum struct test_floor_type { none, standard };
enum struct test_instance_name_type { material, shape };
enum struct test_shapes_type {
  // clang-format off
  features1, features2, rows, bunny_sphere,
  shapes1, shapes2, shapes3
  // clang-format off
};
enum struct test_materials_type {
  // clang-format off
  features1, features2, uvgrid, hair, plastic_metal,
  materials1, materials2, materials3, materials4, materials5,
  // clang-format on
};

struct test_params {
  test_cameras_type       cameras       = test_cameras_type::standard;
  test_environments_type  environments  = test_environments_type::sky;
  test_arealights_type    arealights    = test_arealights_type::standard;
  test_floor_type         floor         = test_floor_type::standard;
  test_shapes_type        shapes        = test_shapes_type::features1;
  test_materials_type     materials     = test_materials_type::features1;
  test_instance_name_type instance_name = test_instance_name_type::material;
};

// Scene test
void make_test(scene_scene& scene, const test_params& params) {
  // cameras
  switch (params.cameras) {
    case test_cameras_type::standard: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
    // TODO(fabio): fix wide camera
    case test_cameras_type::wide: {
      add_camera(scene, "default", {-0.75, 0.4, 0.9}, {-0.075, 0.05, -0.05},
          {0, 1, 0}, 0.05, 2.4, 0);
    } break;
  }
  // TODO(fabio): port other cameras
  switch (params.environments) {
    case test_environments_type::none: break;
    case test_environments_type::sky: {
      add_environment(scene, "sky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, false, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
    } break;
    case test_environments_type::sunsky: {
      add_environment(scene, "sunsky", identity3x4f, {0.5, 0.5, 0.5},
          add_texture(scene, "sky",
              make_sunsky(
                  {2048, 1024}, pif / 4, 3.0, true, 1.0, 1.0, {0.7, 0.7, 0.7}),
              true));
    } break;
  }
  switch (params.arealights) {
    case test_arealights_type::none: break;
    case test_arealights_type::standard: {
      add_instance(scene, "arealight1",
          lookat_frame({-0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight1", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(
              scene, "arealight1", {20, 20, 20}, invalid_handle));
      add_instance(scene, "arealight2",
          lookat_frame({+0.4, 0.8, 0.8}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "arealight2", make_rect({1, 1}, {0.2, 0.2})),
          add_emission_material(
              scene, "arealight2", {20, 20, 20}, invalid_handle));
    } break;
    case test_arealights_type::large: {
      add_instance(scene, "largearealight1",
          lookat_frame({-0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "largearealight1", make_rect({1, 1}, {0.4, 0.4})),
          add_emission_material(
              scene, "largearealight1", {10, 10, 10}, invalid_handle));
      add_instance(scene, "largearealight2",
          lookat_frame({+0.8, 1.6, 1.6}, {0, 0.1, 0}, {0, 1, 0}, true),
          add_shape(scene, "largearealight2", make_rect({1, 1}, {0.4, 0.4})),
          add_emission_material(
              scene, "largearealight2", {10, 10, 10}, invalid_handle));
    } break;
  }
  switch (params.floor) {
    case test_floor_type::none: break;
    case test_floor_type::standard: {
      add_instance(scene, "floor", identity3x4f,
          add_shape(scene, "floor", make_floor({1, 1}, {2, 2}, {20, 20})),
          add_matte_material(scene, "floor", {1, 1, 1},
              add_texture(scene, "floor", make_grid({1024, 1024}))));
    } break;
  }
  auto shapes = vector<shape_handle>{}, shapesi = vector<shape_handle>{};
  auto subdivs   = vector<subdiv_handle>{};
  auto materials = vector<material_handle>{};
  switch (params.shapes) {
    case test_shapes_type::features1: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::features2: {
      shapes  = {add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
          add_shape(scene, "suzanne", make_monkey(0.075f * 0.8f)),
          add_shape(scene, "hair",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
          add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
          add_shape(scene, "cube",
              make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1, 1},
                  0.3 * 0.075f))};
      shapesi = {invalid_handle, invalid_handle,
          add_shape(scene, "hairi", make_sphere(32, 0.075f * 0.8f, 1)),
          invalid_handle, invalid_handle};
      subdivs = {add_subdiv(scene, "suzanne", make_monkey(0.075f * 0.8f),
                     shapes[1], 2),
          add_subdiv(scene, "displaced", make_sphere(128, 0.075f, 1), shapes[3],
              0, 0.025,
              add_texture(scene, "bumps-displacement", make_bumps({1024, 1024}),
                  false, true))};
    } break;
    case test_shapes_type::rows: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, bunny, bunny, bunny, bunny, sphere, sphere, sphere,
          sphere, sphere};
    } break;
    case test_shapes_type::bunny_sphere: {
      auto bunny  = add_shape(scene, "bunny", make_bunny(0.075));
      auto sphere = add_shape(scene, "sphere", make_sphere(32, 0.075, 1));
      shapes      = {bunny, sphere, bunny, sphere, bunny};
    } break;
    case test_shapes_type::shapes1: {
      shapes = {
          add_shape(scene, "sphere", make_sphere(32, 0.075, 1)),
          add_shape(scene, "uvsphere-flipcap",
              make_capped_uvsphere({32, 32}, 0.075, {1, 1}, 0.3 * 0.075)),
          add_shape(scene, "disk", make_disk(32, 0.075f, 1)),
          add_shape(scene, "uvcylinder",
              make_rounded_uvcylinder(
                  {32, 32, 32}, {0.075, 0.075}, {1, 1, 1}, 0.3 * 0.075)),
          add_shape(scene, "cube",
              make_rounded_box({32, 32, 32}, {0.075, 0.075, 0.075}, {1, 1, 1},
                  0.3 * 0.075f)),
      };
    } break;
    case test_shapes_type::shapes2: {
      shapes = {
          add_shape(scene, "cube-subdiv", make_fvcube(0.075)),
          add_shape(scene, "suzanne-subdiv", make_monkey(0.075)),
          add_shape(scene, "displaced", make_sphere(128, 0.075f, 1)),
          add_shape(scene, "bunny", make_bunny(0.075)),
          add_shape(scene, "teapot", make_sphere(32, 0.075, 1)),
      };
      subdivs = {
          add_subdiv(scene, "cube-subdiv", make_fvcube(0.075), shapes[0], 4),
          add_subdiv(scene, "suzanne-subdiv", make_monkey(0.075), shapes[1], 2),
          add_subdiv(scene, "displaced", make_sphere(128, 0.075f, 1), shapes[2],
              0, 0.025,
              add_texture(scene, "bumps-displacement", make_bumps({1024, 1024}),
                  false, true))};
    } break;
    case test_shapes_type::shapes3: {
      shapes = {
          invalid_handle,
          add_shape(scene, "hair1",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0.03, 100})),
          add_shape(scene, "hair2",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f})),
          add_shape(scene, "hair3",
              make_hair(make_sphere(32, 0.075f * 0.8f, 1), {4, 65536},
                  {0.1f * 0.15f, 0.1f * 0.15f},
                  {0.001f * 0.15f, 0.0005f * 0.15f}, {0, 0}, {0.5, 128})),
          invalid_handle,
      };
    } break;
  }
  switch (params.materials) {
    case test_materials_type::features1: {
      materials = {
          add_plastic_material(scene, "coated", {1, 1, 1}, 0.2,
              add_texture(scene, "uvgrid", make_uvgrid({1024, 1024}))),
          add_glass_material(scene, "glass", {1, 0.5, 0.5}, 0),
          add_glass_material(
              scene, "jade", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
          add_plastic_material(scene, "bumped", {0.5, 0.7, 0.5}, 0.2,
              invalid_handle, invalid_handle,
              add_texture(scene, "bumps-normal",
                  bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true)),
          add_metal_material(scene, "metal", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::features2: {
      auto uvgrid  = add_plastic_material(scene, "uvgrid", {1, 1, 1}, 0.2,
          add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
      auto plastic = add_plastic_material(
          scene, "plastic", {0.5, 0.7, 0.5}, 0.2);
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
      materials = {uvgrid, plastic, hair, plastic, uvgrid};
    } break;
    case test_materials_type::uvgrid: {
      auto uvgrid = add_plastic_material(scene, "uvgrid", {1, 1, 1}, 0.2,
          add_texture(scene, "uvgrid", make_uvgrid({1024, 1024})));
      materials   = {uvgrid, uvgrid, uvgrid, uvgrid, uvgrid};
    } break;
    case test_materials_type::hair: {
      auto hair = add_matte_material(scene, "hair", {0.7, 0.7, 0.7});
      materials = {hair, hair, hair, hair, hair};
    } break;
    case test_materials_type::plastic_metal: {
      materials = {
          add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
          add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
          add_metal_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
          add_metal_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials1: {
      materials = {
          add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01),
          add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_matte_material(scene, "matte", {0.7, 0.7, 0.7}),
          add_plastic_material(scene, "metal1", {0.7, 0.7, 0.7}, 0),
          add_plastic_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials2: {
      materials = {
          add_glass_material(scene, "glass1", {1, 1, 1}, 0),
          add_glass_material(scene, "glass2", {1, 0.7, 0.7}, 0.1),
          add_transparent_material(scene, "transparent", {0.7, 0.5, 0.5}, 0.2),
          add_thinglass_material(scene, "tglass1", {1, 1, 1}, 0),
          add_thinglass_material(scene, "tglass2", {1, 0.7, 0.7}, 0.1),
      };
    } break;
    case test_materials_type::materials3: {
      auto bumps_normal = add_texture(scene, "bumps-normal",
          bump_to_normal(make_bumps({1024, 1024}), 0.05), false, true);
      materials         = {
          add_plastic_material(scene, "plastic1", {0.5, 0.5, 0.7}, 0.01,
              invalid_handle, invalid_handle, bumps_normal),
          add_plastic_material(scene, "plastic2", {0.5, 0.7, 0.5}, 0.2),
          add_metal_material(scene, "metal1", {0.7, 0.7, 0.7}, 0,
              invalid_handle, invalid_handle, bumps_normal),
          add_metal_material(scene, "metal2", {0.66, 0.45, 0.34}, 0.2),
          add_metal_material(scene, "metal3", {0.66, 0.45, 0.34}, 0.2),
      };
    } break;
    case test_materials_type::materials4: {
      materials = {
          add_volume_material(
              scene, "cloud", {0.65, 0.65, 0.65}, {0.9, 0.9, 0.9}, 1),
          add_glass_material(scene, "glass", {1, 0.5, 0.5}, 0),
          add_glass_material(
              scene, "jade", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
          add_glass_material(
              scene, "jade2", {0.5, 0.5, 0.5}, 0, {0.3, 0.6, 0.3}),
          add_volume_material(scene, "smoke", {0.5, 0.5, 0.5}, {0.2, 0.2, 0.2}),
      };
    } break;
    case test_materials_type::materials5: {
      materials = {
          add_glass_material(scene, "skin1a", {0.76, 0.48, 0.23}, 0.25,
              {0.436, 0.227, 0.131}, invalid_handle, invalid_handle,
              invalid_handle, 1.5, -0.8, 0.001),
          add_glass_material(scene, "skin2a", {0.82, 0.55, 0.4}, 0.25,
              {0.623, 0.433, 0.343}, invalid_handle, invalid_handle,
              invalid_handle, 1.5, -0.8, 0.001),
          add_glass_material(scene, "skins", {0.76, 0.48, 0.23}, 0,
              {0.436, 0.227, 0.131}, invalid_handle, invalid_handle,
              invalid_handle, 1.5, -0.8, 0.001),
          add_glass_material(scene, "skin1b", {0.76, 0.48, 0.23}, 0.25,
              {0.436, 0.227, 0.131}, invalid_handle, invalid_handle,
              invalid_handle, 1.5, -0.8, 0.001),
          add_glass_material(scene, "skin2b", {0.82, 0.55, 0.4}, 0.25,
              {0.623, 0.433, 0.343}, invalid_handle, invalid_handle,
              invalid_handle, 1.5, -0.8, 0.001),
      };
    } break;
  }
  for (auto idx = 0; idx < shapes.size(); idx++) {
    if (!shapes[idx]) continue;
    if (shapes.size() > 5) {
      add_instance(scene,
          scene.shape_names[idx] + "-" + scene.shape_names[idx % 5],
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1},
              {0.2f * (idx % 5 - 2), 0.075, -0.4f * (idx / 5)}},
          shapes[idx], materials[idx % 5]);
    } else {
      auto name = params.instance_name == test_instance_name_type::material
                      ? scene.material_names[idx]
                      : scene.shape_names[idx];
      add_instance(scene, name,
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx % 5 - 2), 0.075, 0}},
          shapes[idx], materials[idx]);
    }
    if (!shapesi.empty() && shapesi[idx]) {
      // TODO(fabio): fix name
      add_instance(scene, "",
          {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0.2f * (idx - 2), 0.075, 0}},
          shapesi[idx], materials[idx]);
    }
  }
}

// Scene presets used for testing.
bool make_scene_preset(scene_scene& scene, const string& type, string& error) {
  if (type == "cornellbox") {
    make_cornellbox(scene);
    return true;
  } else if (type == "features1") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features1, test_materials_type::features1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "features2") {
    make_test(
        scene, {test_cameras_type::standard, test_environments_type::sky,
                   test_arealights_type::standard, test_floor_type::standard,
                   test_shapes_type::features2, test_materials_type::features2,
                   test_instance_name_type::shape});
    return true;
  } else if (type == "materials1") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials1,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials2") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials2,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials3") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials3,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials4") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials4,
                   test_instance_name_type::material});
    return true;
  } else if (type == "materials5") {
    make_test(
        scene, {test_cameras_type::wide, test_environments_type::sky,
                   test_arealights_type::large, test_floor_type::standard,
                   test_shapes_type::rows, test_materials_type::materials5,
                   test_instance_name_type::material});
    return true;
  } else if (type == "shapes1") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes1, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes2") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes2, test_materials_type::uvgrid,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "shapes3") {
    make_test(scene, {test_cameras_type::standard, test_environments_type::sky,
                         test_arealights_type::large, test_floor_type::standard,
                         test_shapes_type::shapes3, test_materials_type::hair,
                         test_instance_name_type::shape});
    return true;
  } else if (type == "environments1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "environments2") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::sunsky,
            test_arealights_type::none, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else if (type == "arealights1") {
    make_test(scene,
        {test_cameras_type::standard, test_environments_type::none,
            test_arealights_type::standard, test_floor_type::standard,
            test_shapes_type::bunny_sphere, test_materials_type::plastic_metal,
            test_instance_name_type::material});
    return true;
  } else {
    error = "unknown preset";
    return false;
  }
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GENERIC SCENE LOADING
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_json_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to OBJ.
static bool load_obj_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_obj_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to PLY. Loads/saves only one mesh with no other data.
static bool load_ply_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_ply_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to STL. Loads/saves only one mesh with no other data.
static bool load_stl_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_stl_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to glTF.
static bool load_gltf_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_gltf_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene from/to pbrt-> This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
static bool load_pbrt_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);
static bool save_pbrt_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load/save a scene preset.
static bool load_preset_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel);

// Load a scene
bool load_scene(const string& filename, scene_scene& scene, string& error,
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
  } else if (ext == ".ypreset" || ext == ".YPRESET") {
    return load_preset_scene(filename, scene, error, progress_cb, noparallel);
  } else {
    return format_error();
  }
}

// Save a scene
bool save_scene(const string& filename, const scene_scene& scene, string& error,
    const progress_callback& progress_cb, bool noparallel) {
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

// Load/save a scene preset.
static bool load_preset_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto preset_error = [filename, &error]() {
    error = filename + ": " + error;
    return false;
  };

  // handle progress
  if (progress_cb) progress_cb("make preset", 0, 1);

  // make preset
  if (!make_scene_preset(scene, path_basename(filename), error))
    return preset_error();

  // handle progress
  if (progress_cb) progress_cb("make preset", 0, 1);

  // done
  return true;
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
    if (!load_ply(filename, ply, error)) return false;
    get_values(ply, "instance",
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
    add_values(ply, "instance",
        {"xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz", "ox", "oy",
            "oz"},
        frames);
    return save_ply(filename, ply, error);
  } else {
    return format_error();
  }
}

// Loads/saves a  channel float/byte image in linear/srgb color space.
static bool load_image(const string& filename, image<vec4f>& imgf,
    image<vec4b>& imgb, string& error) {
  if (is_hdr_filename(filename)) {
    return load_image(filename, imgf, error);
  } else {
    return load_image(filename, imgb, error);
  }
}

// Loads/saves a  channel float/byte image in linear/srgb color space.
static bool save_image(const string& filename, const image<vec4f>& imgf,
    const image<vec4b>& imgb, string& error) {
  if (!imgf.empty()) {
    return save_image(filename, imgf, error);
  } else {
    return save_image(filename, imgb, error);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

using njson = nlohmann::json;
using std::array;

// load/save json
[[maybe_unused]] static bool load_json(
    const string& filename, njson& json, string& error) {
  // error helpers
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error in json";
    return false;
  };
  auto text = ""s;
  if (!load_text(filename, text, error)) return false;
  try {
    json = njson::parse(text);
    return true;
  } catch (std::exception&) {
    return parse_error();
  }
}

[[maybe_unused]] static bool save_json(
    const string& filename, const njson& json, string& error) {
  return save_text(filename, json.dump(2), error);
}

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
inline void from_json(const njson& j, vec4f& value) {
  nlohmann::from_json(j, (array<float, 4>&)value);
}
inline void from_json(const njson& j, mat3f& value) {
  nlohmann::from_json(j, (array<float, 9>&)value);
}
inline void from_json(const njson& j, frame3f& value) {
  nlohmann::from_json(j, (array<float, 12>&)value);
}
inline void from_json(const njson& j, mat4f& value) {
  nlohmann::from_json(j, (array<float, 16>&)value);
}

inline void to_json(njson& j, material_type value) {
  j = material_type_names.at((int)value);
}
inline void from_json(const njson& j, material_type& value) {
  auto values = j.get<string>();
  auto pos    = std::find(
      material_type_names.begin(), material_type_names.end(), values);
  if (pos == material_type_names.end())
    throw std::invalid_argument{"unknown value"};
  value = (material_type)(pos - material_type_names.begin());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// JSON IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load a scene in the builtin JSON format.
static bool load_json_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto json_error = [filename]() {
    // error does not need setting
    return false;
  };
  auto parse_error = [filename, &error](const string& patha,
                         const string& pathb = "", const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    error = filename + ": parse error at " + path;
    return false;
  };
  auto key_error = [filename, &error](const string& patha,
                       const string& pathb = "", const string& pathc = "") {
    auto path = patha;
    if (!pathb.empty()) path += "/" + pathb;
    if (!pathc.empty()) path += "/" + pathc;
    error = filename + ": unknow key at " + path;
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
  auto js = njson{};
  if (!load_json(filename, js, error)) return json_error();

  // reference disctionaries
  auto texture_map = unordered_map<string, pair<texture_handle, bool>>{
      {"", {invalid_handle, true}}};
  auto shape_map = unordered_map<string, pair<shape_handle, bool>>{
      {"", {invalid_handle, true}}};
  auto material_map = unordered_map<string, pair<material_handle, bool>>{
      {"", {invalid_handle, true}}};
  auto subdiv_map = unordered_map<string, pair<subdiv_handle, bool>>{
      {"", {invalid_handle, true}}};

  // parse json value
  auto get_value = [](const njson& js, auto& value) -> bool {
    try {
      value = js;
      return true;
    } catch (...) {
      return false;
    }
  };

  // parse json reference
  auto get_shape = [&scene, &shape_map, &get_value](
                       const njson& js, shape_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = shape_map.find(name);
    if (it != shape_map.end()) {
      value = it->second.first;
      return it->second.first != invalid_handle;
    }
    scene.shape_names.emplace_back(name);
    scene.shapes.emplace_back();
    auto shape_id   = (int)scene.shapes.size() - 1;
    shape_map[name] = {shape_id, false};
    value           = shape_id;
    return true;
  };

  // parse json reference
  auto get_material = [&scene, &material_map, &get_value](
                          const njson& js, material_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = material_map.find(name);
    if (it != material_map.end()) {
      value = it->second.first;
      return it->second.first != invalid_handle;
    }
    scene.material_names.emplace_back(name);
    scene.materials.emplace_back();
    auto material_id   = (int)scene.materials.size() - 1;
    material_map[name] = {material_id, false};
    value              = material_id;
    return true;
  };

  // parse json reference
  auto get_texture = [&scene, &texture_map, &get_value](
                         const njson& js, texture_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = texture_map.find(name);
    if (it != texture_map.end()) {
      value = it->second.first;
      return it->second.first != invalid_handle;
    }
    scene.texture_names.emplace_back(name);
    scene.textures.emplace_back();
    auto texture_id   = (int)scene.textures.size() - 1;
    texture_map[name] = {texture_id, false};
    value             = texture_id;
    return true;
  };

  struct ply_instance {
    vector<frame3f> frames = {};
  };

  // load json instance
  using ply_instance_handle = int;
  auto ply_instances        = vector<ply_instance>{};
  auto ply_instance_map     = unordered_map<string, ply_instance_handle>{
      {"", invalid_handle}};
  auto instance_ply = unordered_map<instance_handle, ply_instance_handle>{};
  auto get_ply_instances = [&ply_instances, &ply_instance_map, &instance_ply,
                               &get_value](const njson& js,
                               instance_handle          instance) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    if (name.empty()) return true;
    auto it = ply_instance_map.find(name);
    if (it != ply_instance_map.end()) {
      instance_ply[instance] = it->second;
      return true;
    }
    ply_instances.emplace_back(ply_instance());
    ply_instance_map[name] = (int)ply_instances.size() - 1;
    instance_ply[instance] = (int)ply_instances.size() - 1;
    return true;
  };
  auto get_instance_handle =
      [](const scene_scene&     scene,
          const scene_instance& instance) -> instance_handle {
    return (instance_handle)(&instance - scene.instances.data());
  };

  // helper for iteration
  auto iterate_object = [](const njson& js) { return js.items(); };
  auto check_object   = [](const njson& js) { return js.is_object(); };

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // loop over external dictionaries
  for (auto& [gname, group] : iterate_object(js)) {
    if (gname == "asset") {
      auto& asset = scene.asset;
      if (!check_object(group)) return parse_error(gname);
      for (auto& [key, value] : iterate_object(group)) {
        if (key == "copyright") {
          if (!get_value(value, asset.copyright))
            return parse_error(gname, key);
        } else if (key == "generator") {
          if (!get_value(value, asset.generator))
            return parse_error(gname, key);
        } else {
          return key_error(gname, key);
        }
      }
    } else if (gname == "cameras") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.camera_names.emplace_back(name);
        auto& camera = scene.cameras.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, camera.frame))
              return parse_error(gname, name, key);
          } else if (key == "orthographic") {
            if (!get_value(value, camera.orthographic))
              return parse_error(gname, name, key);
          } else if (key == "ortho") {
            // backward compatibility
            if (!get_value(value, camera.orthographic))
              return parse_error(gname, name, key);
          } else if (key == "lens") {
            if (!get_value(value, camera.lens))
              return parse_error(gname, name, key);
          } else if (key == "aspect") {
            if (!get_value(value, camera.aspect))
              return parse_error(gname, name, key);
          } else if (key == "film") {
            if (!get_value(value, camera.film))
              return parse_error(gname, name, key);
          } else if (key == "focus") {
            if (!get_value(value, camera.focus))
              return parse_error(gname, name, key);
          } else if (key == "aperture") {
            if (!get_value(value, camera.aperture))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)camera.frame))
              return parse_error(gname, name, key);
            camera.focus = length(camera.frame.x - camera.frame.y);
            camera.frame = lookat_frame(
                camera.frame.x, camera.frame.y, camera.frame.z);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "environments") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.environment_names.emplace_back(name);
        auto& environment = scene.environments.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, environment.frame))
              return parse_error(gname, name, key);
          } else if (key == "emission") {
            if (!get_value(value, environment.emission))
              return parse_error(gname, name, key);
          } else if (key == "emission_tex") {
            if (!get_texture(value, environment.emission_tex))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)environment.frame))
              return parse_error(gname, name, key);
            environment.frame = lookat_frame(environment.frame.x,
                environment.frame.y, environment.frame.z, true);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "materials") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        auto material_it = material_map.find(string{name});
        auto handle      = invalid_handle;
        if (material_it == material_map.end()) {
          scene.material_names.emplace_back(name);
          scene.materials.emplace_back();
          handle = (int)scene.materials.size() - 1;
        } else {
          handle = material_it->second.first;
        }
        auto& material = scene.materials[handle];
        material.type  = material_type::metallic;
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "type") {
            if (!get_value(value, material.type))
              return parse_error(gname, name, key);
          } else if (key == "emission") {
            if (!get_value(value, material.emission))
              return parse_error(gname, name, key);
          } else if (key == "color") {
            if (!get_value(value, material.color))
              return parse_error(gname, name, key);
          } else if (key == "metallic") {
            if (!get_value(value, material.metallic))
              return parse_error(gname, name, key);
          } else if (key == "roughness") {
            if (!get_value(value, material.roughness))
              return parse_error(gname, name, key);
          } else if (key == "ior") {
            if (!get_value(value, material.ior))
              return parse_error(gname, name, key);
          } else if (key == "trdepth") {
            if (!get_value(value, material.trdepth))
              return parse_error(gname, name, key);
          } else if (key == "scattering") {
            if (!get_value(value, material.scattering))
              return parse_error(gname, name, key);
          } else if (key == "scanisotropy") {
            if (!get_value(value, material.scanisotropy))
              return parse_error(gname, name, key);
          } else if (key == "opacity") {
            if (!get_value(value, material.opacity))
              return parse_error(gname, name, key);
          } else if (key == "emission_tex") {
            if (!get_texture(value, material.emission_tex))
              return parse_error(gname, name, key);
          } else if (key == "color_tex") {
            if (!get_texture(value, material.color_tex))
              return parse_error(gname, name, key);
          } else if (key == "roughness_tex") {
            if (!get_texture(value, material.roughness_tex))
              return parse_error(gname, name, key);
          } else if (key == "scattering_tex") {
            if (!get_texture(value, material.scattering_tex))
              return parse_error(gname, name, key);
          } else if (key == "normal_tex") {
            if (!get_texture(value, material.normal_tex))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
        material_map[string{name}] = {handle, true};
      }
    } else if (gname == "instances" || gname == "objects") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.instance_names.emplace_back(name);
        auto& instance = scene.instances.emplace_back();
        for (auto [key, value] : iterate_object(element)) {
          if (key == "frame") {
            if (!get_value(value, instance.frame))
              return parse_error(gname, name, key);
          } else if (key == "lookat") {
            if (!get_value(value, (mat3f&)instance.frame))
              return parse_error(gname, name, key);
            instance.frame = lookat_frame(
                instance.frame.x, instance.frame.y, instance.frame.z, true);
          } else if (key == "material") {
            if (!get_material(value, instance.material))
              return parse_error(gname, name, key);
          } else if (key == "shape") {
            if (!get_shape(value, instance.shape))
              return parse_error(gname, name, key);
          } else if (key == "instance") {
            throw std::invalid_argument{"instances not suppotyed yet"};
            // get_ply_instances(value, instance);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "subdivs") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.subdiv_names.emplace_back();
        auto& subdiv = scene.subdivs.emplace_back();
        for (auto& [key, value] : iterate_object(element)) {
          if (key == "shape") {
            if (!get_shape(value, subdiv.shape))
              return parse_error(gname, name, key);
          } else if (key == "subdivisions") {
            if (!get_value(value, subdiv.subdivisions))
              return parse_error(gname, name, key);
          } else if (key == "catmullcark") {
            if (!get_value(value, subdiv.catmullclark))
              return parse_error(gname, name, key);
          } else if (key == "smooth") {
            if (!get_value(value, subdiv.smooth))
              return parse_error(gname, name, key);
          } else if (key == "displacement") {
            if (!get_value(value, subdiv.displacement))
              return parse_error(gname, name, key);
          } else if (key == "displacement_tex") {
            if (!get_texture(value, subdiv.displacement_tex))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
        subdiv_map[string{name}] = {(int)scene.subdivs.size() - 1, false};
      }
    } else {
      return key_error(gname);
    }
  }

  // check materials
  for (auto& [key, value] : material_map) {
    if (!value.second) return material_error(key);
  }

  // handle progress
  progress.y += scene.shapes.size();
  progress.y += scene.textures.size();
  progress.y += scene.subdivs.size();
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
    auto& shape = scene.shapes[value.first];
    if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
    auto path   = make_filename(name, "shapes", {".ply", ".obj"});
    auto lshape = shape_data{};
    if (!load_shape(path, lshape, error, false)) return dependent_error();
    shape.points    = lshape.points;
    shape.lines     = lshape.lines;
    shape.triangles = lshape.triangles;
    shape.quads     = lshape.quads;
    shape.positions = lshape.positions;
    shape.normals   = lshape.normals;
    shape.texcoords = lshape.texcoords;
    shape.colors    = lshape.colors;
    shape.radius    = lshape.radius;
  }

  // load subdivs
  subdiv_map.erase("");
  for (auto [name, value] : subdiv_map) {
    auto& subdiv = scene.subdivs[value.first];
    if (progress_cb) progress_cb("load subdiv", progress.x++, progress.y);
    auto path    = make_filename(name, "subdivs", {".ply", ".obj"});
    auto lsubdiv = fvshape_data{};
    if (!load_fvshape(path, lsubdiv, error, true)) return dependent_error();
    subdiv.quadspos      = lsubdiv.quadspos;
    subdiv.quadsnorm     = lsubdiv.quadsnorm;
    subdiv.quadstexcoord = lsubdiv.quadstexcoord;
    subdiv.positions     = lsubdiv.positions;
    subdiv.normals       = lsubdiv.normals;
    subdiv.texcoords     = lsubdiv.texcoords;
  }

  // load textures
  texture_map.erase("");
  for (auto [name, value] : texture_map) {
    auto& texture = scene.textures[value.first];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    auto path = make_filename(
        name, "textures", {".hdr", ".exr", ".png", ".jpg"});
    if (!load_image(path, texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // load instances
  ply_instance_map.erase("");
  for (auto& [name, ihandle] : ply_instance_map) {
    throw std::invalid_argument{"instances not supported"};
    // auto& instance = ply_instances.at(ihandle);
    // if (progress_cb) progress_cb("load instance", progress.x++, progress.y);
    // auto path = make_filename(name, "instances", {".ply"});
    // if (!load_instance(path, instance->frames, error)) return
    // dependent_error();
  }

  // apply instances
  if (!ply_instances.empty()) {
    throw std::invalid_argument{"not supported for now"};
    // if (progress_cb)
    //   progress_cb("flatten instances", progress.x++, progress.y++);
    // auto instances = scene.instances;
    // scene.instances.clear();
    // for (auto instance : instances) {
    //   auto it = instance_ply.find(instance);
    //   if (it == instance_ply.end()) {
    //     auto ninstance      = add_instance(scene, instance->name);
    //     ninstance->frame    = instance->frame;
    //     ninstance->shape    = instance->shape;
    //     ninstance->material = instance->material;
    //   } else {
    //     auto ply_instance = it->second;
    //     for (auto& frame : ply_instance->frames) {
    //       auto ninstance      = add_instance(scene, instance->name);
    //       ninstance->frame    = frame * instance->frame;
    //       ninstance->shape    = instance->shape;
    //       ninstance->material = instance->material;
    //     }
    //   }
    // }
    // for (auto instance : instances) delete instance;
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);
  trim_memory(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a scene in the builtin JSON format.
static bool save_json_scene(const string& filename, const scene_scene& scene,
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

  // helpers to handel old code paths
  auto insert_object = [](njson& js, const string& name) -> njson& {
    js[name] = njson::object();
    return js[name];
  };
  auto insert_value = [](njson& js, const string& name, const auto& value) {
    js[name] = value;
  };

  // handle progress
  auto progress = vec2i{
      0, 2 + (int)scene.shapes.size() + (int)scene.textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json file
  auto js = njson::object();

  // asset
  {
    auto& asset   = scene.asset;
    auto& element = insert_object(js, "asset");
    if (!asset.copyright.empty()) {
      insert_value(element, "copyright", asset.copyright);
    }
    if (!asset.generator.empty()) {
      insert_value(element, "generator", asset.generator);
    }
  }

  auto def_cam = sceneio_camera{};
  if (!scene.cameras.empty()) {
    auto& group = insert_object(js, "cameras");
    for (auto& camera : scene.cameras) {
      auto& element = insert_object(group, get_camera_name(scene, camera));
      if (camera.frame != def_cam.frame) {
        insert_value(element, "frame", camera.frame);
      }
      if (camera.orthographic != def_cam.orthographic) {
        insert_value(element, "orthographic", camera.orthographic);
      }
      if (camera.lens != def_cam.lens) {
        insert_value(element, "lens", camera.lens);
      }
      if (camera.aspect != def_cam.aspect) {
        insert_value(element, "aspect", camera.aspect);
      }
      if (camera.film != def_cam.film) {
        insert_value(element, "film", camera.film);
      }
      if (camera.focus != def_cam.focus) {
        insert_value(element, "focus", camera.focus);
      }
      if (camera.aperture != def_cam.aperture) {
        insert_value(element, "aperture", camera.aperture);
      }
    }
  }

  auto def_env = sceneio_environment{};
  if (!scene.environments.empty()) {
    auto& group = insert_object(js, "environments");
    for (auto& environment : scene.environments) {
      auto& element = insert_object(
          group, get_environment_name(scene, environment));
      if (environment.frame != def_env.frame) {
        insert_value(element, "frame", environment.frame);
      }
      if (environment.emission != def_env.emission) {
        insert_value(element, "emission", environment.emission);
      }
      if (environment.emission_tex != invalid_handle) {
        insert_value(element, "emission_tex",
            get_texture_name(scene, environment.emission_tex));
      }
    }
  }

  auto def_material = sceneio_material{};
  if (!scene.materials.empty()) {
    auto& group = insert_object(js, "materials");
    for (auto& material : scene.materials) {
      auto& element = insert_object(group, get_material_name(scene, material));
      if (material.type != def_material.type) {
        insert_value(element, "type", material.type);
      }
      if (material.emission != def_material.emission) {
        insert_value(element, "emission", material.emission);
      }
      if (material.color != def_material.color) {
        insert_value(element, "color", material.color);
      }
      if (material.metallic != def_material.metallic) {
        insert_value(element, "metallic", material.metallic);
      }
      if (material.roughness != def_material.roughness) {
        insert_value(element, "roughness", material.roughness);
      }
      if (material.ior != def_material.ior) {
        insert_value(element, "ior", material.ior);
      }
      if (material.trdepth != def_material.trdepth) {
        insert_value(element, "trdepth", material.trdepth);
      }
      if (material.scattering != def_material.scattering) {
        insert_value(element, "scattering", material.scattering);
      }
      if (material.scanisotropy != def_material.scanisotropy) {
        insert_value(element, "scanisotropy", material.scanisotropy);
      }
      if (material.opacity != def_material.opacity) {
        insert_value(element, "opacity", material.opacity);
      }
      if (material.emission_tex != invalid_handle) {
        insert_value(element, "emission_tex",
            get_texture_name(scene, material.emission_tex));
      }
      if (material.color_tex != invalid_handle) {
        insert_value(
            element, "color_tex", get_texture_name(scene, material.color_tex));
      }
      if (material.roughness_tex != invalid_handle) {
        insert_value(element, "roughness_tex",
            get_texture_name(scene, material.roughness_tex));
      }
      if (material.scattering_tex != invalid_handle) {
        insert_value(element, "scattering_tex",
            get_texture_name(scene, material.scattering_tex));
      }
      if (material.normal_tex != invalid_handle) {
        insert_value(element, "normal_tex",
            get_texture_name(scene, material.normal_tex));
      }
    }
  }

  auto def_instance = sceneio_instance{};
  if (!scene.instances.empty()) {
    auto& group = insert_object(js, "instances");
    for (auto& instance : scene.instances) {
      auto& element = insert_object(group, get_instance_name(scene, instance));
      if (instance.frame != def_instance.frame) {
        insert_value(element, "frame", instance.frame);
      }
      if (instance.shape != invalid_handle) {
        insert_value(element, "shape", get_shape_name(scene, instance.shape));
      }
      if (instance.material != invalid_handle) {
        insert_value(
            element, "material", get_material_name(scene, instance.material));
      }
    }
  }

  auto def_subdiv = scene_subdiv{};
  if (!scene.subdivs.empty()) {
    auto& group = insert_object(js, "subdivs");
    for (auto& subdiv : scene.subdivs) {
      auto& element = insert_object(group, get_subdiv_name(scene, subdiv));
      if (subdiv.shape != invalid_handle) {
        insert_value(element, "shape", get_shape_name(scene, subdiv.shape));
      }
      if (subdiv.subdivisions != def_subdiv.subdivisions) {
        insert_value(element, "subdivisions", subdiv.subdivisions);
      }
      if (subdiv.catmullclark != def_subdiv.catmullclark) {
        insert_value(element, "catmullclark", subdiv.catmullclark);
      }
      if (subdiv.smooth != def_subdiv.smooth) {
        insert_value(element, "smooth", subdiv.smooth);
      }
      if (subdiv.displacement != def_subdiv.displacement) {
        insert_value(element, "displacement", subdiv.displacement);
      }
      if (subdiv.displacement_tex != invalid_handle) {
        insert_value(element, "displacement_tex",
            get_texture_name(scene, subdiv.displacement_tex));
      }
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
  for (auto& shape : scene.shapes) {
    if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
    auto path = make_filename(get_shape_name(scene, shape), "shapes", ".ply"s);
    auto sshape      = shape_data{};
    sshape.points    = shape.points;
    sshape.lines     = shape.lines;
    sshape.triangles = shape.triangles;
    sshape.quads     = shape.quads;
    sshape.positions = shape.positions;
    sshape.normals   = shape.normals;
    sshape.texcoords = shape.texcoords;
    sshape.colors    = shape.colors;
    sshape.radius    = shape.radius;
    if (!save_shape(path, sshape, error, false)) return dependent_error();
  }

  // save subdiv
  for (auto& subdiv : scene.subdivs) {
    if (progress_cb) progress_cb("save subdiv", progress.x++, progress.y);
    auto path = make_filename(
        get_subdiv_name(scene, subdiv), "subdivs", ".obj");
    auto ssubdiv          = fvshape_data{};
    ssubdiv.quadspos      = subdiv.quadspos;
    ssubdiv.quadsnorm     = subdiv.quadsnorm;
    ssubdiv.quadstexcoord = subdiv.quadstexcoord;
    ssubdiv.positions     = subdiv.positions;
    ssubdiv.normals       = subdiv.normals;
    ssubdiv.texcoords     = subdiv.texcoords;
    if (!save_fvshape(path, ssubdiv, error, true)) return dependent_error();
  }

  // save textures
  for (auto& texture : scene.textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(get_texture_name(scene, texture), "textures",
        (!texture.hdr.empty()) ? ".hdr"s : ".png"s);
    if (!save_image(path, texture.hdr, texture.ldr, error))
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

#define YOCTO_FASTOBJ

#ifdef YOCTO_FASTOBJ

// Loads an OBJ
static bool load_obj_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto read_error = [filename, &error]() {
    error = filename + ": error reading obj";
    return false;
  };
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
  auto obj = fast_obj_read(filename.c_str());
  if (!obj) return read_error();
  auto obj_guard = std::unique_ptr<fastObjMesh, void (*)(fastObjMesh*)>(
      obj, &fast_obj_destroy);

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // helper to create texture maps
  auto texture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_texture = [&texture_map, &scene](
                         const fastObjTexture& tinfo) -> texture_handle {
    if (tinfo.name == nullptr) return invalid_handle;
    auto path = string{tinfo.name};
    for (auto& c : path)
      if (c == '\\') c = '/';
    if (path.empty()) return invalid_handle;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    scene.textures.emplace_back();
    texture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert between roughness and exponent
  auto exponent_to_roughness = [](float exponent) {
    auto roughness = exponent;
    roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
    if (roughness < 0.01f) roughness = 0;
    if (roughness > 0.99f) roughness = 1;
    return roughness;
  };

  // handler for cameras
  for (auto idx = 0; idx < obj->material_count; idx++) {
    auto omaterial = obj->materials + idx;
    if (omaterial->name == nullptr) continue;
    if (string{omaterial->name}.find("ycamera_") != 0) continue;
    auto& camera = scene.cameras.emplace_back();
    camera.frame = lookat_frame(
        vec3f{omaterial->Kd[0], omaterial->Kd[1], omaterial->Kd[2]},
        vec3f{omaterial->Ks[0], omaterial->Ks[1], omaterial->Ks[2]}, {0, 1, 0});
    camera.focus = distance(
        vec3f{omaterial->Kd[0], omaterial->Kd[1], omaterial->Kd[2]},
        vec3f{omaterial->Ks[0], omaterial->Ks[1], omaterial->Ks[2]});
    camera.aspect = omaterial->Ni;
    camera.lens   = omaterial->Ns / 1000.0f;
  }

  // handler for environments
  for (auto idx = 0; idx < obj->material_count; idx++) {
    auto omaterial = obj->materials + idx;
    if (omaterial->name == nullptr) continue;
    if (string{omaterial->name}.find("yenvironment_") != 0) continue;
    auto& environment = scene.environments.emplace_back();
    environment.frame = lookat_frame(
        vec3f{omaterial->Kd[0], omaterial->Kd[1], omaterial->Kd[2]},
        vec3f{omaterial->Ks[0], omaterial->Ks[1], omaterial->Ks[2]}, {0, 1, 0},
        true);
    environment.emission = vec3f{
        omaterial->Ke[0], omaterial->Ke[1], omaterial->Ke[2]};
    environment.emission_tex = get_texture(omaterial->map_Ke);
  }

  // handler for materials
  auto material_map = unordered_map<string, material_handle>{};
  for (auto idx = 0; idx < obj->material_count; idx++) {
    auto omaterial = obj->materials + idx;
    if (omaterial->name == nullptr) continue;
    if (string{omaterial->name}.find("ycamera_") == 0) continue;
    if (string{omaterial->name}.find("yenvironment_") == 0) continue;
    auto& material    = scene.materials.emplace_back();
    material.type     = material_type::metallic;
    material.emission = vec3f{
        omaterial->Ke[0], omaterial->Ke[1], omaterial->Ke[2]};
    material.emission_tex = get_texture(omaterial->map_Ke);
    if (max(max(omaterial->Kt[0], omaterial->Kt[1]), omaterial->Kt[2]) > 0.1) {
      material.type  = material_type::thinglass;
      material.color = vec3f{
          omaterial->Kt[0], omaterial->Kt[1], omaterial->Kt[2]};
      material.color_tex = get_texture(omaterial->map_Kt);
      material.roughness = exponent_to_roughness(omaterial->Ns);
      material.ior       = 1.5;
      material.metallic  = 0;
      material.opacity   = 1;
    } else if (min(min(omaterial->Tf[0], omaterial->Tf[1]), omaterial->Tf[2]) <
               0.99) {
      material.type      = material_type::thinglass;
      material.color     = vec3f{1, 1, 1};
      material.color_tex = get_texture(omaterial->map_Kt);
      material.roughness = 0;
      material.ior       = 1.5;
      material.metallic  = 0;
      material.opacity   = 1;
    } else if (max(max(omaterial->Ks[0], omaterial->Ks[1]), omaterial->Ks[2]) >
               0.2) {
      material.type  = material_type::metal;
      material.color = vec3f{
          omaterial->Ks[0], omaterial->Ks[1], omaterial->Ks[2]};
      material.color_tex = get_texture(omaterial->map_Ks);
      material.roughness = exponent_to_roughness(omaterial->Ns);
      material.ior       = 1.5;
      material.metallic  = 0;
      material.opacity   = omaterial->d;
    } else if (max(max(omaterial->Ks[0], omaterial->Ks[1]), omaterial->Ks[2]) >
               0) {
      material.type  = material_type::plastic;
      material.color = vec3f{
          omaterial->Kd[0], omaterial->Kd[1], omaterial->Kd[2]};
      material.color_tex = get_texture(omaterial->map_Kd);
      material.roughness = exponent_to_roughness(omaterial->Ns);
      material.ior       = 1.5;
      material.metallic  = 0;
      material.opacity   = omaterial->d;
    } else {
      material.type  = material_type::matte;
      material.color = vec3f{
          omaterial->Kd[0], omaterial->Kd[1], omaterial->Kd[2]};
      material.color_tex = get_texture(omaterial->map_Kd);
      material.roughness = exponent_to_roughness(omaterial->Ns);
      material.ior       = 1.5;
      material.metallic  = 0;
      material.opacity   = omaterial->d;
    }
    material.normal_tex           = get_texture(omaterial->map_bump);
    material_map[omaterial->name] = (int)scene.materials.size() - 1;
  }

  struct index_hash {
    size_t operator()(const vec3i& v) const {
      const std::hash<int> hasher = std::hash<int>();
      auto                 h      = (size_t)0;
      h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  // convert shapes
  auto shape_map  = unordered_map<string, int>{};
  auto shape_vmap = vector<unordered_map<vec3i, int, index_hash>>{};
  for (auto idx = 0; idx < obj->group_count; idx++) {
    auto group = obj->groups + idx;
    if (group->face_count == 0) continue;
    auto cur_mat          = -1;
    auto cur_shape        = -1;
    auto cur_index_offset = group->index_offset;
    for (auto cur_face = group->face_offset;
         cur_face < group->face_offset + group->face_count; cur_face++) {
      if (cur_mat != obj->face_materials[cur_face]) {
        cur_mat       = obj->face_materials[cur_face];
        auto cur_name = string{group->name != nullptr ? group->name : ""} +
                        "@@@" + std::to_string(cur_mat);
        if (shape_map.find(cur_name) != shape_map.end()) {
          cur_shape = shape_map.at(cur_name);
        } else {
          scene.shapes.emplace_back();
          shape_vmap.emplace_back();
          cur_shape           = (int)scene.shapes.size() - 1;
          auto& instance      = scene.instances.emplace_back();
          instance.material   = cur_mat;
          instance.shape      = cur_shape;
          shape_map[cur_name] = cur_shape;
        }
      }
      auto& shape = scene.shapes[cur_shape];
      auto& vmap  = shape_vmap[cur_shape];
      auto  vids  = array<int, 128>{};
      for (auto vidx = 0; vidx < obj->face_vertices[cur_face]; vidx++) {
        auto indices  = obj->indices[cur_index_offset + vidx];
        auto vindices = vec3i{(int)indices.p, (int)indices.n, (int)indices.t};
        auto vert_it  = vmap.find(vindices);
        if (vert_it == vmap.end()) {
          shape.positions.push_back({obj->positions[indices.p * 3 + 0],
              obj->positions[indices.p * 3 + 1],
              obj->positions[indices.p * 3 + 2]});
          if (!shape.normals.empty() || indices.n != 0) {
            shape.normals.push_back({obj->normals[indices.n * 3 + 0],
                obj->normals[indices.n * 3 + 1],
                obj->normals[indices.n * 3 + 2]});
          }
          if (!shape.texcoords.empty() || indices.t != 0) {
            shape.texcoords.push_back({obj->texcoords[indices.t * 2 + 0],
                1 - obj->texcoords[indices.t * 2 + 1]});
          }
          vids[vidx] = (int)shape.positions.size() - 1;
          vmap.insert(vert_it, {vindices, vids[vidx]});
        } else {
          vids[vidx] = vert_it->second;
        }
      }
      if (obj->face_vertices[cur_face] == 3) {
        shape.triangles.push_back({vids[0], vids[1], vids[2]});
      } else if (obj->face_vertices[cur_face] == 4) {
        shape.quads.push_back({vids[0], vids[1], vids[2], vids[3]});
      } else if (obj->face_vertices[cur_face] > 4) {
        for (auto vidx = 2; vidx < obj->face_vertices[cur_face]; vidx++) {
          shape.triangles.push_back({vids[0], vids[vidx - 1], vids[vidx]});
        }
      } else {
        // not supported
      }
      cur_index_offset += obj->face_vertices[cur_face];
    }
  }

  // handle mixed shapes
  for (auto& shape : scene.shapes) {
    if (!shape.quads.empty() && !shape.triangles.empty()) {
      auto tquads = quads_to_triangles(shape.quads);
      shape.triangles.insert(
          shape.triangles.end(), tquads.begin(), tquads.end());
      shape.quads.clear();
    }
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name) {
    return path_join(path_dirname(filename), name);
  };

  // load textures
  texture_map.erase("");
  for (auto [name, thandle] : texture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

#else

// Loads an OBJ
static bool load_obj_scene(const string& filename, scene_scene& scene,
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
  auto obj = obj_scene{};
  if (!load_obj(filename, obj, error, false, true, false)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto& ocamera : obj.cameras) {
    auto& camera        = scene.cameras.emplace_back();
    camera.frame        = ocamera.frame;
    camera.orthographic = ocamera.ortho;
    camera.film         = max(ocamera.width, ocamera.height);
    camera.aspect       = ocamera.width / ocamera.height;
    camera.focus        = ocamera.focus;
    camera.lens         = ocamera.lens;
    camera.aperture     = ocamera.aperture;
  }

  // helper to create texture maps
  auto texture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_texture = [&texture_map, &scene](
                         const obj_texture& tinfo) -> texture_handle {
    auto path = tinfo.path;
    if (path.empty()) return invalid_handle;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    scene.textures.emplace_back();
    texture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert between roughness and exponent
  auto exponent_to_roughness = [](float exponent) {
    auto roughness = exponent;
    roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
    if (roughness < 0.01f) roughness = 0;
    if (roughness > 0.99f) roughness = 1;
    return roughness;
  };

  // handler for materials
  auto material_map = unordered_map<string, material_handle>{};
  for (auto& omaterial : obj.materials) {
    auto& material        = scene.materials.emplace_back();
    material.type         = material_type::metallic;
    material.emission     = omaterial.emission;
    material.emission_tex = get_texture(omaterial.emission_tex);
    if (max(omaterial.transmission) > 0.1) {
      material.type      = material_type::thinglass;
      material.color     = omaterial.transmission;
      material.color_tex = get_texture(omaterial.transmission_tex);
    } else if (max(omaterial.specular) > 0.2) {
      material.type      = material_type::metal;
      material.color     = omaterial.specular;
      material.color_tex = get_texture(omaterial.specular_tex);
    } else if (max(omaterial.specular) > 0) {
      material.type      = material_type::plastic;
      material.color     = omaterial.diffuse;
      material.color_tex = get_texture(omaterial.diffuse_tex);
    } else {
      material.type      = material_type::matte;
      material.color     = omaterial.diffuse;
      material.color_tex = get_texture(omaterial.diffuse_tex);
    }
    material.roughness           = exponent_to_roughness(omaterial.exponent);
    material.ior                 = omaterial.ior;
    material.metallic            = 0;
    material.opacity             = omaterial.opacity;
    material.normal_tex          = get_texture(omaterial.normal_tex);
    material_map[omaterial.name] = (int)scene.materials.size() - 1;
  }

  // convert shapes
  auto shape_name_counts = unordered_map<string, int>{};
  for (auto& oshape : obj.shapes) {
    auto& materials = oshape.materials;
    if (materials.empty()) materials.push_back(nullptr);
    for (auto material_idx = 0; material_idx < materials.size();
         material_idx++) {
      auto& shape = scene.shapes.emplace_back();
      if (material_map.find(materials[material_idx]) == material_map.end())
        return material_error(materials[material_idx]);
      auto material   = material_map.at(materials[material_idx]);
      auto has_quads_ = has_quads(oshape);
      if (!oshape.faces.empty() && !has_quads_) {
        get_triangles(oshape, material_idx, shape.triangles, shape.positions,
            shape.normals, shape.texcoords, true);
      } else if (!oshape.faces.empty() && has_quads_) {
        get_quads(oshape, material_idx, shape.quads, shape.positions,
            shape.normals, shape.texcoords, true);
      } else if (!oshape.lines.empty()) {
        get_lines(oshape, material_idx, shape.lines, shape.positions,
            shape.normals, shape.texcoords, true);
      } else if (!oshape.points.empty()) {
        get_points(oshape, material_idx, shape.points, shape.positions,
            shape.normals, shape.texcoords, true);
      } else {
        return shape_error();
      }
      auto shape_id = (int)scene.shapes.size() - 1;
      if (oshape.instances.empty()) {
        auto& instance    = scene.instances.emplace_back();
        instance.shape    = shape_id;
        instance.material = material;
      } else {
        for (auto& frame : oshape.instances) {
          auto instance     = scene.instances.emplace_back();
          instance.frame    = frame;
          instance.shape    = shape_id;
          instance.material = material;
        }
      }
    }
  }

  // convert environments
  for (auto& oenvironment : obj.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = get_texture(oenvironment.emission_tex);
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name) {
    return path_join(path_dirname(filename), name);
  };

  // load textures
  texture_map.erase("");
  for (auto [name, thandle] : texture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

#endif

static bool save_obj_scene(const string& filename, const scene_scene& scene,
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
  auto progress = vec2i{0, 2 + (int)scene.textures.size()};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // build obj
  auto obj = obj_scene{};

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = get_camera_name(scene, camera);
    ocamera.frame    = camera.frame;
    ocamera.ortho    = camera.orthographic;
    ocamera.width    = camera.film;
    ocamera.height   = camera.film / camera.aspect;
    ocamera.focus    = camera.focus;
    ocamera.lens     = camera.lens;
    ocamera.aperture = camera.aperture;
  }

  // textures
  auto get_texture = [&](texture_handle texture) {
    if (texture == invalid_handle) return obj_texture{};
    auto tinfo  = obj_texture{};
    auto is_hdr = !scene.textures[texture].hdr.empty();
    tinfo.path  = "textures/" + get_texture_name(scene, texture) +
                 (is_hdr ? ".hdr"s : ".png"s);
    return tinfo;
  };

  auto roughness_to_exponent = [](float roughness) -> float {
    if (roughness < 0.01f) return 10000;
    if (roughness > 0.99f) return 10;
    return 2 / pow(roughness, 4.0f) - 2;
  };

  // convert materials and textures
  for (auto& material : scene.materials) {
    auto& omaterial        = obj.materials.emplace_back();
    omaterial.name         = get_material_name(scene, material);
    omaterial.illum        = 2;
    omaterial.emission     = material.emission;
    omaterial.diffuse      = material.color;
    omaterial.specular     = {0, 0, 0};
    omaterial.exponent     = roughness_to_exponent(material.roughness);
    omaterial.opacity      = material.opacity;
    omaterial.emission_tex = get_texture(material.emission_tex);
    omaterial.diffuse_tex  = get_texture(material.color_tex);
    omaterial.normal_tex   = get_texture(material.normal_tex);
  }

  // convert objects
  for (auto& instance : scene.instances) {
    auto& shape     = scene.shapes[instance.shape];
    auto  positions = shape.positions, normals = shape.normals;
    for (auto& p : positions) p = transform_point(instance.frame, p);
    for (auto& n : normals) n = transform_normal(instance.frame, n);
    auto& oshape     = obj.shapes.emplace_back();
    oshape.name      = get_shape_name(scene, shape);
    oshape.materials = {get_material_name(scene, instance.material)};
    if (!shape.triangles.empty()) {
      set_triangles(oshape, shape.triangles, positions, normals,
          shape.texcoords, {}, true);
    } else if (!shape.quads.empty()) {
      set_quads(
          oshape, shape.quads, positions, normals, shape.texcoords, {}, true);
    } else if (!shape.lines.empty()) {
      set_lines(
          oshape, shape.lines, positions, normals, shape.texcoords, {}, true);
    } else if (!shape.points.empty()) {
      set_points(
          oshape, shape.points, positions, normals, shape.texcoords, {}, true);
    } else {
      return shape_error();
    }
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = get_environment_name(scene, environment);
    oenvironment.frame        = environment.frame;
    oenvironment.emission     = environment.emission;
    oenvironment.emission_tex = get_texture(environment.emission_tex);
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
  auto texture_id = 0;
  for (auto& texture : scene.textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(get_texture_name(scene, texture_id++), "textures",
        (!texture.hdr.empty()) ? ".hdr"s : ".png"s);
    if (!save_image(path, texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_ply_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load ply mesh
  auto& shape  = scene.shapes.emplace_back();
  auto  lshape = shape_data{};
  if (!load_shape(filename, lshape, error, false)) return false;
  shape.points    = lshape.points;
  shape.lines     = lshape.lines;
  shape.triangles = lshape.triangles;
  shape.quads     = lshape.quads;
  shape.positions = lshape.positions;
  shape.normals   = lshape.normals;
  shape.texcoords = lshape.texcoords;
  shape.colors    = lshape.colors;
  shape.radius    = lshape.radius;

  // create instance
  auto& instance = scene.instances.emplace_back();
  instance.shape = (int)scene.shapes.size() - 1;

  // fix scene
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_ply_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  if (scene.shapes.empty()) return shape_error();

  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save shape
  auto& shape      = scene.shapes.front();
  auto  sshape     = shape_data{};
  sshape.points    = shape.points;
  sshape.lines     = shape.lines;
  sshape.triangles = shape.triangles;
  sshape.quads     = shape.quads;
  sshape.positions = shape.positions;
  sshape.normals   = shape.normals;
  sshape.texcoords = shape.texcoords;
  sshape.colors    = shape.colors;
  sshape.radius    = shape.radius;
  if (!save_shape(filename, sshape, error)) return false;

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// STL CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool load_stl_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load stl mesh
  auto& shape  = scene.shapes.emplace_back();
  auto  lshape = shape_data{};
  if (!load_shape(filename, lshape, error, false)) return false;
  shape.points    = lshape.points;
  shape.lines     = lshape.lines;
  shape.triangles = lshape.triangles;
  shape.quads     = lshape.quads;
  shape.positions = lshape.positions;
  shape.normals   = lshape.normals;
  shape.texcoords = lshape.texcoords;
  shape.colors    = lshape.colors;
  shape.radius    = lshape.radius;

  // create instance
  auto& instance = scene.instances.emplace_back();
  instance.shape = (int)scene.shapes.size() - 1;

  // fix scene
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

static bool save_stl_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto shape_error = [filename, &error]() {
    error = filename + ": empty shape";
    return false;
  };

  if (scene.shapes.empty()) return shape_error();

  // handle progress
  auto progress = vec2i{0, 1};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save shape
  auto& shape      = scene.shapes.front();
  auto  sshape     = shape_data{};
  sshape.points    = shape.points;
  sshape.lines     = shape.lines;
  sshape.triangles = shape.triangles;
  sshape.quads     = shape.quads;
  sshape.positions = shape.positions;
  sshape.normals   = shape.normals;
  sshape.texcoords = shape.texcoords;
  sshape.colors    = shape.colors;
  sshape.radius    = shape.radius;
  if (!save_shape(filename, sshape, error, false)) return false;

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

// #define YOCTO_CGLTF_IN

#ifdef YOCTO_CGLTF_IN

// Load a scene
static bool load_gltf_scene(const string& filename, scene_scene& scene,
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
    if (gast->copyright != nullptr) scene.asset.copyright = gast->copyright;
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
    auto  gcam          = gnde->camera;
    auto& camera        = scene.cameras.emplace_back();
    camera.frame        = mat_to_frame(mat);
    camera.orthographic = gcam->type == cgltf_camera_type_orthographic;
    if (camera.orthographic) {
      auto ortho    = &gcam->data.orthographic;
      camera.aspect = ortho->xmag / ortho->ymag;
      camera.lens   = ortho->ymag;  // this is probably bogus
      camera.film   = 0.036;
    } else {
      auto persp    = &gcam->data.perspective;
      camera.aspect = persp->aspect_ratio;
      if (camera.aspect == 0) camera.aspect = 16.0f / 9.0f;
      camera.film = 0.036;
      if (camera.aspect >= 1) {
        camera.lens = (camera.film / camera.aspect) /
                      (2 * tan(persp->yfov / 2));
      } else {
        camera.lens = camera.film / (2 * tan(persp->yfov / 2));
      }
      camera.focus = 1;
    }
  }

  // convert color textures
  auto texture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_texture = [&scene, &texture_map](
                         const cgltf_texture_view& ginfo) -> texture_handle {
    if (ginfo.texture == nullptr || ginfo.texture->image == nullptr)
      return invalid_handle;
    auto path = string{ginfo.texture->image->uri};
    if (path.empty()) return invalid_handle;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    scene.textures.emplace_back();
    texture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert materials
  auto material_map = unordered_map<cgltf_material*, material_handle>{
      {nullptr, invalid_handle}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto  gmaterial       = &gltf->materials[mid];
    auto& material        = scene.materials.emplace_back();
    material.type         = material_type::metallic;
    material.emission     = {gmaterial->emissive_factor[0],
        gmaterial->emissive_factor[1], gmaterial->emissive_factor[2]};
    material.emission_tex = get_texture(gmaterial->emissive_texture);
    if (gmaterial->has_pbr_metallic_roughness != 0) {
      auto gmr          = &gmaterial->pbr_metallic_roughness;
      material.color    = {gmr->base_color_factor[0], gmr->base_color_factor[1],
          gmr->base_color_factor[2]};
      material.opacity  = gmr->base_color_factor[3];
      material.metallic = gmr->metallic_factor;
      material.roughness      = gmr->roughness_factor;
      material.color_tex      = get_texture(gmr->base_color_texture);
      material.roughness_tex  = get_texture(gmr->metallic_roughness_texture);
      material.normal_tex     = get_texture(gmaterial->normal_texture);
      material_map[gmaterial] = (int)scene.materials.size() - 1;
    }
  }

  // convert meshes
  auto mesh_map = unordered_map<cgltf_mesh*, vector<sceneio_instance>>{
      {nullptr, {}}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh = &gltf->meshes[mid];
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (gprim->attributes_count == 0) continue;
      auto& shape       = scene.shapes.emplace_back();
      auto  instance    = scene_instance{};
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = material_map.at(gprim->material);
      mesh_map[gmesh].push_back(instance);
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name != nullptr ? gattr->name : "");
        auto gacc     = gattr->data;
        if (semantic == "POSITION") {
          shape.positions.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape.positions[i].x, 3);
        } else if (semantic == "NORMAL") {
          shape.normals.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape.normals[i].x, 3);
        } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
          shape.texcoords.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape.texcoords[i].x, 2);
        } else if (semantic == "COLOR" || semantic == "COLOR_0") {
          shape.colors.resize(gacc->count);
          if (cgltf_num_components(gacc->type) == 3) {
            for (auto i = 0; i < gacc->count; i++)
              cgltf_accessor_read_float(gacc, i, &shape.colors[i].x, 3);
          } else {
            for (auto i = 0; i < gacc->count; i++) {
              auto color4 = vec4f{0, 0, 0, 0};
              cgltf_accessor_read_float(gacc, i, &color4.x, 4);
              shape.colors[i] = color4;
            }
          }
        } else if (semantic == "TANGENT") {
          shape.tangents.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape.tangents[i].x, 4);
          for (auto& t : shape.tangents) t.w = -t.w;
        } else if (semantic == "_RADIUS") {
          shape.radius.resize(gacc->count);
          for (auto i = 0; i < gacc->count; i++)
            cgltf_accessor_read_float(gacc, i, &shape.radius[i], 1);
        } else {
          // ignore
        }
      }
      // indices
      if (gprim->indices == nullptr) {
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.resize(shape.positions.size() / 3);
          for (auto i = 0; i < shape.positions.size() / 3; i++)
            shape.triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.resize(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles[i - 2] = {0, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.resize(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles[i - 2] = {i - 2, i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.resize(shape.positions.size() / 2);
          for (auto i = 0; i < shape.positions.size() / 2; i++)
            shape.lines[i] = {i * 2 + 0, i * 2 + 1};
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.resize(shape.positions.size());
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines[i - 1] = {i - 1, i};
          shape.lines.back() = {(int)shape.positions.size() - 1, 0};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.resize(shape.positions.size() - 1);
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines[i - 1] = {i - 1, i};
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          return primitive_error();
        } else {
          return primitive_error();
        }
      } else {
        auto giacc = gprim->indices;
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.resize(giacc->count / 3);
          for (auto i = 0; i < giacc->count / 3; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 3 + 0, (uint*)&shape.triangles[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 1, (uint*)&shape.triangles[i].y, 1);
            cgltf_accessor_read_uint(
                giacc, i * 3 + 2, (uint*)&shape.triangles[i].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, 0 + 0, (uint*)&shape.triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape.triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape.triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.resize(giacc->count - 2);
          for (auto i = 2; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, i - 2, (uint*)&shape.triangles[i - 2].x, 1);
            cgltf_accessor_read_uint(
                giacc, i - 1, (uint*)&shape.triangles[i - 2].y, 1);
            cgltf_accessor_read_uint(
                giacc, i + 0, (uint*)&shape.triangles[i - 2].z, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.resize(giacc->count / 2);
          for (auto i = 0; i < giacc->count / 2; i++) {
            cgltf_accessor_read_uint(
                giacc, i * 2 + 0, (uint*)&shape.lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, i * 2 + 1, (uint*)&shape.lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.resize(giacc->count);
          for (auto i = 0; i < giacc->count; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape.lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape.lines[i].y, 1);
          }
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.resize(giacc->count - 1);
          for (auto i = 0; i < giacc->count - 1; i++) {
            cgltf_accessor_read_uint(
                giacc, (i + 0) % giacc->count, (uint*)&shape.lines[i].x, 1);
            cgltf_accessor_read_uint(
                giacc, (i + 1) % giacc->count, (uint*)&shape.lines[i].y, 1);
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
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    if (!visible_nodes[nid]) continue;
    auto gnde = &gltf->nodes[nid];
    if (gnde->mesh == nullptr) continue;
    auto mat = mat4f{};
    cgltf_node_transform_world(gnde, &mat.x.x);
    for (auto& prims : mesh_map.at(gnde->mesh)) {
      auto& instance = scene.instances.emplace_back();
      instance       = prims;
      auto mat       = mat4f{};
      cgltf_node_transform_world(gnde, &mat.x.x);
      instance.frame = mat_to_frame(mat);
    }
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // load texture
  texture_map.erase("");
  for (auto [tpath, thandle] : texture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(path_join(path_dirname(filename), tpath), texture.hdr,
            texture.ldr, error))
      return dependent_error();
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // fix cameras
  // auto bbox = compute_bounds(scene);
  // for (auto& camera : scene.cameras) {
  //   auto center   = (bbox.min + bbox.max) / 2;
  //   auto distance = dot(-camera.frame.z, center - camera.frame.o);
  //   if (distance > 0) camera.focus = distance;
  // }

  // load done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

#else

// Load a scene
static bool load_gltf_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto primitive_error = [filename, &error]() {
    error = filename + ": primitive error";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
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
  auto gltf = njson{};
  if (!load_json(filename, gltf, error)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load buffers
  auto bnames = vector<string>{};
  try {
    if (gltf.contains("buffers")) {
      for (auto& gbuffer : gltf.at("buffers")) {
        if (!gbuffer.contains("uri")) return parse_error();
        bnames.push_back(gbuffer.value("uri", ""));
      }
    }
  } catch (...) {
    return parse_error();
  }
  auto buffers = vector<vector<byte>>();
  buffers.reserve(bnames.size());
  auto dirname = path_dirname(filename);
  for (auto& name : bnames) {
    progress_cb("load buffer", progress.x++, progress.y);
    if (!load_binary(path_join(dirname, name), buffers.emplace_back(), error))
      return dependent_error();
  }

  // handle progress
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // convert asset
  if (gltf.contains("asset")) {
    try {
      scene.asset.copyright = gltf.value("copyright", ""s);
    } catch (...) {
      return parse_error();
    }
  }

  // convert cameras
  auto cameras = vector<scene_camera>{};
  if (gltf.contains("cameras")) {
    try {
      for (auto& gcamera : gltf.at("cameras")) {
        auto& camera = cameras.emplace_back();
        auto  type   = gcamera.value("type", "perspective");
        if (type == "orthographic") {
          auto& gortho  = gcamera.at("orthographic");
          auto  xmag    = gortho.value("xmag", 1.0f);
          auto  ymag    = gortho.value("ymag", 1.0f);
          camera.aspect = xmag / ymag;
          camera.lens   = ymag;  // this is probably bogus
          camera.film   = 0.036;
        } else if (type == "perspective") {
          auto& gpersp  = gcamera.at("perspective");
          camera.aspect = gpersp.value("aspectRatio", 0.0f);
          auto yfov     = gpersp.value("yfov", radians(45));
          if (camera.aspect == 0) camera.aspect = 16.0f / 9.0f;
          camera.film = 0.036;
          if (camera.aspect >= 1) {
            camera.lens = (camera.film / camera.aspect) / (2 * tan(yfov / 2));
          } else {
            camera.lens = camera.film / (2 * tan(yfov / 2));
          }
          camera.focus = 1;
        } else {
          return parse_error();
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  // prepare list of effective nodes
  auto visible_nodes = vector<bool>{};
  if (gltf.contains("nodes")) {
    try {
      visible_nodes.assign(gltf.at("nodes").size(), true);
    } catch (...) {
      return parse_error();
    }
  }

  // convert color textures
  auto texture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_texture = [&gltf, &scene, &texture_map](const njson& js,
                         const string& name) -> texture_handle {
    if (!js.contains(name)) return invalid_handle;
    auto& ginfo    = js.at(name);
    auto& gtexture = gltf.at("textures").at(ginfo.value("index", -1));
    auto& gimage   = gltf.at("images").at(gtexture.value("source", -1));
    auto  path     = gimage.value("uri", "");
    if (path.empty()) return invalid_handle;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    scene.textures.emplace_back();
    texture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // convert materials
  if (gltf.contains("materials")) {
    try {
      for (auto& gmaterial : gltf.at("materials")) {
        auto& material    = scene.materials.emplace_back();
        material.type     = material_type::metallic;
        material.emission = gmaterial.value("emissiveFactor", vec3f{0, 0, 0});
        material.emission_tex = get_texture(gmaterial, "emissiveTexture");
        material.normal_tex   = get_texture(gmaterial, "normalTexture");
        if (gmaterial.contains("pbrMetallicRoughness")) {
          auto& gpbr         = gmaterial.at("pbrMetallicRoughness");
          auto  base         = gpbr.value("baseColorFactor", vec4f{1, 1, 1, 1});
          material.color     = xyz(base);
          material.opacity   = base.w;
          material.metallic  = gpbr.value("metallicFactor", 1.0f);
          material.roughness = gpbr.value("roughnessFactor", 1.0f);
          material.color_tex = get_texture(gpbr, "baseColorTexture");
          material.roughness_tex = get_texture(
              gpbr, "metallicRoughnessTexture");
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  auto get_attribute = [](const njson&                 gltf,
                           const vector<vector<byte>>& buffers, auto& values,
                           int index, bool color = false) {
    constexpr auto ncomp     = sizeof(values.front()) / sizeof(float);
    auto&          gaccessor = gltf.at("accessors").at(index);
    if (gaccessor.contains("sparse"))
      throw std::invalid_argument{"not implemented"};
    auto  type      = gaccessor.value("type", "VEC3");
    auto  component = gaccessor.value("type", 5126);
    auto  count     = gaccessor.value("count", (size_t)0);
    auto  offset1   = gaccessor.value("byteOffset", (size_t)0);
    auto& gview = gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
    auto  offset2 = gview.value("byteOffset", (size_t)0);
    auto  stride  = gview.value("byteStride", (size_t)0);
    auto& buffer  = buffers.at(gview.value("buffer", -1));
    values.resize(count);
    auto data    = (float*)values.data();
    auto current = buffer.data() + offset1 + offset2;
    auto nncomp  = 0;
    if (type == "SCALAR") nncomp = 1;
    if (type == "VEC2") nncomp = 2;
    if (type == "VEC3") nncomp = 3;
    if (type == "VEC4") nncomp = 4;
    if (ncomp != nncomp) throw std::invalid_argument{"bad value"};
    if (component == 5121) {
      if (stride == 0) stride += nncomp * 1;
      for (auto idx = (size_t)0; idx < count; idx++, current += stride) {
        for (auto c = 0; c < nncomp; c++) {
          data[idx * ncomp + c] = *(byte*)current / 255.0f;
        }
      }
    } else if (component == 5123) {
      if (stride == 0) stride += nncomp * 2;
      for (auto idx = (size_t)0; idx < count; idx++, current += stride) {
        for (auto c = 0; c < nncomp; c++) {
          data[idx * ncomp + c] = *(ushort*)current / 65535.0f;
        }
      }
    } else if (component == 5126) {
      if (stride == 0) stride += nncomp * 4;
      for (auto idx = (size_t)0; idx < count; idx++, current += stride) {
        for (auto c = 0; c < nncomp; c++) {
          data[idx * ncomp + c] = *(float*)current;
        }
      }
    } else {
      throw std::invalid_argument{"bad value"};
    }
  };

  // convert meshes
  auto mesh_primitives = vector<vector<sceneio_instance>>{};
  if (gltf.contains("meshes")) {
    try {
      auto type_components = unordered_map<string, int>{
          {"SCALAR", 1}, {"VEC2", 2}, {"VEC3", 3}, {"VEC4", 4}};
      for (auto& gmesh : gltf.at("meshes")) {
        auto& primitives = mesh_primitives.emplace_back();
        if (!gmesh.contains("primitives")) continue;
        for (auto& gprimitive : gmesh.at("primitives")) {
          if (!gprimitive.contains("attributes")) continue;
          auto& shape       = scene.shapes.emplace_back();
          auto& instance    = primitives.emplace_back();
          instance.shape    = (int)scene.shapes.size() - 1;
          instance.material = gprimitive.value("material", -1);
          for (auto& [gname, gattribute] :
              gprimitive.at("attributes").items()) {
            auto& gaccessor = gltf.at("accessors").at(gattribute.get<int>());
            if (gaccessor.contains("sparse"))
              throw std::invalid_argument{"sparse accessor"};
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer      = buffers.at(gview.value("buffer", 0));
            auto  components  = type_components.at(gaccessor.value("type", ""));
            auto  dcomponents = components;
            auto  count       = gaccessor.value("count", (size_t)0);
            auto  data        = (float*)nullptr;
            if (gname == "POSITION") {
              if (components != 3)
                throw std::invalid_argument{"invalid accessor"};
              shape.positions.resize(count);
              data = (float*)shape.positions.data();
            } else if (gname == "NORMAL") {
              if (components != 3)
                throw std::invalid_argument{"invalid accessor"};
              shape.normals.resize(count);
              data = (float*)shape.normals.data();
            } else if (gname == "TEXCOORD" || gname == "TEXCOORD_0") {
              if (components != 2)
                throw std::invalid_argument{"invalid accessor"};
              shape.texcoords.resize(count);
              data = (float*)shape.texcoords.data();
            } else if (gname == "COLOR" || gname == "COLOR_0") {
              if (components != 3 && components != 4)
                throw std::invalid_argument{"invalid accessor"};
              shape.colors.resize(count);
              data = (float*)shape.colors.data();
              if (components == 3) {
                dcomponents = 4;
                for (auto& c : shape.colors) c.w = 1;
              }
            } else if (gname == "TANGENT") {
              if (components != 4)
                throw std::invalid_argument{"invalid accessor"};
              shape.tangents.resize(count);
              data = (float*)shape.tangents.data();
            } else if (gname == "RADIUS") {
              if (components != 1)
                throw std::invalid_argument{"invalid accessor"};
              shape.radius.resize(count);
              data = (float*)shape.radius.data();
            } else {
              // ignore
            }
            // convert values
            auto current = buffer.data() +
                           gaccessor.value("byteOffset", (size_t)0) +
                           gview.value("byteOffset", (size_t)0);
            auto stride = gaccessor.value("byteStride", (size_t)0);
            auto ctype  = gaccessor.value("componentType", 0);
            if (ctype == 5121) {
              if (stride == 0) stride = components * 1;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] =
                      *(byte*)(current + comp * 1) / 255.0f;
                }
              }
            } else if (ctype == 5123) {
              if (stride == 0) stride = components * 2;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] =
                      *(ushort*)(current + comp * 2) / 65535.0f;
                }
              }
            } else if (ctype == 5126) {
              if (stride == 0) stride = components * 4;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                for (auto comp = 0; comp < components; comp++) {
                  data[idx * dcomponents + comp] = *(
                      float*)(current + comp * 4);
                }
              }
            } else {
              throw std::invalid_argument{"invalid accessor"};
            }
            // fixes
            if (gname == "TANGENT") {
              for (auto& t : shape.tangents) t.w = -t.w;
            }
          }
          // mode
          auto mode = gprimitive.value("mode", 4);
          // indices
          if (!gprimitive.contains("indices")) {
            if (mode == 4) {  // triangles
              shape.triangles.resize(shape.positions.size() / 3);
              for (auto i = 0; i < shape.positions.size() / 3; i++)
                shape.triangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
            } else if (mode == 6) {  // fans
              shape.triangles.resize(shape.positions.size() - 2);
              for (auto i = 2; i < shape.positions.size(); i++)
                shape.triangles[i - 2] = {0, i - 1, i};
            } else if (mode == 5) {  // strips
              shape.triangles.resize(shape.positions.size() - 2);
              for (auto i = 2; i < shape.positions.size(); i++)
                shape.triangles[i - 2] = {i - 2, i - 1, i};
            } else if (mode == 1) {  // lines
              shape.lines.resize(shape.positions.size() / 2);
              for (auto i = 0; i < shape.positions.size() / 2; i++)
                shape.lines[i] = {i * 2 + 0, i * 2 + 1};
            } else if (mode == 2) {  // lines loops
              shape.lines.resize(shape.positions.size());
              for (auto i = 1; i < shape.positions.size(); i++)
                shape.lines[i - 1] = {i - 1, i};
              shape.lines.back() = {(int)shape.positions.size() - 1, 0};
            } else if (mode == 3) {  // lines strips
              shape.lines.resize(shape.positions.size() - 1);
              for (auto i = 1; i < shape.positions.size(); i++)
                shape.lines[i - 1] = {i - 1, i};
            } else if (mode == 0) {  // points strips
              // points
              return primitive_error();
            } else {
              return primitive_error();
            }
          } else {
            auto& gaccessor =
                gltf.at("accessors").at(gprimitive.value("indices", -1));
            auto& gview =
                gltf.at("bufferViews").at(gaccessor.value("bufferView", -1));
            auto& buffer = buffers.at(gview.value("buffer", 0));
            if (gaccessor.value("type", "") != "SCALAR")
              throw std::invalid_argument{"invalid accessor"};
            auto count   = gaccessor.value("count", (size_t)0);
            auto indices = vector<int>(count);
            // convert values
            auto current = buffer.data() +
                           gaccessor.value("byteOffset", (size_t)0) +
                           gview.value("byteOffset", (size_t)0);
            auto stride = gaccessor.value("byteStride", (size_t)0);
            auto ctype  = gaccessor.value("componentType", 0);
            if (ctype == 5121) {
              if (stride == 0) stride = 1;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(byte*)current;
              }
            } else if (ctype == 5123) {
              if (stride == 0) stride = 2;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(ushort*)current;
              }
            } else if (ctype == 5125) {
              if (stride == 0) stride = 4;
              for (auto idx = 0; idx < count; idx++, current += stride) {
                indices[idx] = (int)*(uint*)current;
              }
            } else {
              throw std::invalid_argument{"invalid accessor"};
            }
            if (mode == 4) {  // triangles
              shape.triangles.resize(indices.size() / 3);
              for (auto i = 0; i < (int)indices.size() / 3; i++) {
                shape.triangles[i] = {
                    indices[i * 3 + 0], indices[i * 3 + 1], indices[i * 3 + 2]};
              }
            } else if (mode == 6) {  // fans
              shape.triangles.resize(indices.size() - 2);
              for (auto i = 2; i < (int)indices.size(); i++) {
                shape.triangles[i - 2] = {
                    indices[0], indices[i - 1], indices[i + 0]};
              }
            } else if (mode == 5) {  // strips
              shape.triangles.resize(indices.size() - 2);
              for (auto i = 2; i < (int)indices.size(); i++) {
                shape.triangles[i - 2] = {
                    indices[i - 2], indices[i - 1], indices[i + 0]};
              }
            } else if (mode == 1) {  // lines
              shape.lines.resize(indices.size() / 2);
              for (auto i = 0; i < (int)indices.size() / 2; i++) {
                shape.lines[i] = {indices[i * 2 + 0], indices[i * 2 + 1]};
              }
            } else if (mode == 2) {  // lines loops
              shape.lines.resize(indices.size());
              for (auto i = 0; i < (int)indices.size(); i++) {
                shape.lines[i] = {
                    indices[i + 0], indices[i + 1] % (int)indices.size()};
              }
            } else if (mode == 3) {  // lines strips
              shape.lines.resize(indices.size() - 1);
              for (auto i = 0; i < (int)indices.size() - 1; i++) {
                shape.lines[i] = {indices[i + 0], indices[i + 1]};
              }
            } else if (mode == 0) {  // points strips
              // points
              return primitive_error();
            } else {
              return primitive_error();
            }
          }
        }
      }
    } catch (...) {
      return parse_error();
    }
  }

  // convert nodes
  if (gltf.contains("nodes")) {
    try {
      auto parents = vector<int>(gltf.at("nodes").size(), -1);
      auto lxforms = vector<frame3f>(gltf.at("nodes").size(), identity3x4f);
      auto node_id = 0;
      for (auto& gnode : gltf.at("nodes")) {
        auto& xform = lxforms.at(node_id);
        if (gnode.contains("matrix")) {
          xform = mat_to_frame(gnode.value("matrix", identity4x4f));
        }
        if (gnode.contains("scale")) {
          xform = scaling_frame(gnode.value("scale", vec3f{1, 1, 1})) * xform;
        }
        if (gnode.contains("rotation")) {
          xform = rotation_frame(gnode.value("rotation", vec4f{0, 0, 0, 1})) *
                  xform;
        }
        if (gnode.contains("translation")) {
          xform = translation_frame(
                      gnode.value("translation", vec3f{0, 0, 0})) *
                  xform;
        }
        if (gnode.contains("children")) {
          for (auto& gchild : gnode.at("children")) {
            parents.at(gchild.get<int>()) = node_id;
          }
        }
        node_id++;
      }
      auto xforms = vector<frame3f>(gltf.at("nodes").size(), identity3x4f);
      node_id     = 0;
      for (auto& gnode : gltf.at("nodes")) {
        if (!gnode.contains("camera") && !gnode.contains("mesh")) {
          node_id++;
          continue;
        }
        auto& xform = xforms.at(node_id);
        xform       = lxforms.at(node_id);
        auto parent = parents.at(node_id);
        while (parent >= 0) {
          xform  = lxforms.at(parent) * xform;
          parent = parents.at(parent);
        }
        if (gnode.contains("camera")) {
          auto& camera = scene.cameras.emplace_back();
          camera       = cameras.at(gnode.value("camera", -1));
          camera.frame = xform;
        }
        if (gnode.contains("mesh")) {
          for (auto& primitive : mesh_primitives.at(gnode.value("mesh", -1))) {
            auto& instance = scene.instances.emplace_back();
            instance       = primitive;
            instance.frame = xform;
          }
        }
        node_id++;
      }
    } catch (...) {
      return parse_error();
    }
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // load texture
  texture_map.erase("");
  for (auto [tpath, thandle] : texture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(path_join(path_dirname(filename), tpath), texture.hdr,
            texture.ldr, error))
      return dependent_error();
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // load done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

#endif

#define YOCTO_CGLTF_OUT

#ifdef YOCTO_CGLTF_OUT

// Load a scene
static bool save_gltf_scene(const string& filename, const scene_scene& scene,
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

  // C helpers
  auto alloc_arrays = [](auto& values, size_t& num, size_t count) {
    num    = count;
    values = (std::remove_reference_t<decltype(values)>)malloc(
        sizeof(*values) * count);
    memset(values, 0, sizeof(*values) * count);
  };
  auto clear_value = [](auto& value) { memset(&value, 0, sizeof(value)); };
  auto copy_string = [](const string& str) {
    if (str.empty()) return (char*)nullptr;
    return strdup(str.c_str());
  };

  // handle progress
  auto progress = vec2i{0, 3};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // convert scene to json
  auto gltf_guard = unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      new cgltf_data{}, &cgltf_free};
  auto& gltf = *gltf_guard.get();
  clear_value(gltf);
  gltf.memory.free = [](void* user, void* ptr) {
    if (ptr) free(ptr);
  };

  // asset
  {
    gltf.asset.generator = copy_string(scene.asset.generator);
    gltf.asset.copyright = copy_string(scene.asset.copyright);
    gltf.asset.version   = copy_string("2.0");
  }

  // cameras
  if (!scene.cameras.empty()) {
    alloc_arrays(gltf.cameras, gltf.cameras_count, scene.cameras.size());
    for (auto& camera : scene.cameras) {
      auto& gcamera       = gltf.cameras[&camera - &scene.cameras.front()];
      gcamera.name        = copy_string(get_camera_name(scene, camera));
      gcamera.type        = cgltf_camera_type_perspective;
      auto& gpersp        = gcamera.data.perspective;
      gpersp.aspect_ratio = camera.aspect;
      gpersp.yfov         = 0.660593;  // TODO(fabio): yfov
      gpersp.znear        = 0.001;     // TODO(fabio): configurable?
      gpersp.zfar         = -1.0f;
    }
  }

  // textures
  if (!scene.textures.empty()) {
    alloc_arrays(gltf.textures, gltf.textures_count, scene.textures.size());
    alloc_arrays(gltf.samplers, gltf.samplers_count, scene.textures.size());
    alloc_arrays(gltf.images, gltf.images_count, scene.textures.size());
    for (auto& texture : scene.textures) {
      auto& gimage = gltf.images[&texture - &scene.textures.front()];
      gimage.name  = copy_string(get_texture_name(scene, texture));
      gimage.uri   = copy_string(
          "textures/" + get_texture_name(scene, texture) + ".png");
      auto& gsampler   = gltf.samplers[&texture - &scene.textures.front()];
      gsampler.wrap_s  = 10497;
      gsampler.wrap_t  = 10497;
      auto& gtexture   = gltf.textures[&texture - &scene.textures.front()];
      gtexture.name    = copy_string(get_texture_name(scene, texture));
      gtexture.image   = &gimage;
      gtexture.sampler = &gsampler;
    }
  }

  // materials
  if (!scene.materials.empty()) {
    alloc_arrays(gltf.materials, gltf.materials_count, scene.materials.size());
    for (auto& material : scene.materials) {
      auto& gmaterial = gltf.materials[&material - &scene.materials.front()];
      gmaterial.name  = copy_string(get_material_name(scene, material));
      gmaterial.alpha_cutoff               = 0.5f;
      gmaterial.emissive_factor[0]         = material.emission.x;
      gmaterial.emissive_factor[1]         = material.emission.y;
      gmaterial.emissive_factor[2]         = material.emission.z;
      gmaterial.has_pbr_metallic_roughness = true;
      auto& gpbr                           = gmaterial.pbr_metallic_roughness;
      gpbr.base_color_factor[0]            = material.color.x;
      gpbr.base_color_factor[1]            = material.color.y;
      gpbr.base_color_factor[2]            = material.color.z;
      gpbr.base_color_factor[3]            = material.opacity;
      gpbr.metallic_factor                 = material.metallic;
      gpbr.roughness_factor                = material.roughness;
      if (material.emission_tex != invalid_handle) {
        gmaterial.emissive_texture.texture =
            &gltf.textures[material.emission_tex];
      }
      if (material.normal_tex != invalid_handle) {
        gmaterial.normal_texture.texture = &gltf.textures[material.normal_tex];
      }
      if (material.color_tex != invalid_handle) {
        gpbr.base_color_texture.texture = &gltf.textures[material.color_tex];
      }
      if (material.roughness_tex != invalid_handle) {
        gpbr.metallic_roughness_texture.texture =
            &gltf.textures[material.roughness_tex];
      }
    }
  }

  // mesh dictionary
  using mesh_key = pair<shape_handle, material_handle>;
  struct mesh_key_hash {
    size_t operator()(const mesh_key& v) const {
      const std::hash<element_handle> hasher = std::hash<element_handle>();
      auto                            h      = (size_t)0;
      h ^= hasher(v.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };
  auto mesh_map = unordered_map<mesh_key, size_t, mesh_key_hash>{};

  // shapes
  if (!scene.shapes.empty()) {
    // meshes
    auto shape_materials = vector<vector<int>>(scene.shapes.size());
    for (auto& instance : scene.instances) {
      if (mesh_map.find({instance.shape, instance.material}) != mesh_map.end())
        continue;
      auto index                                    = mesh_map.size();
      mesh_map[{instance.shape, instance.material}] = index;
      shape_materials[instance.shape].push_back(instance.material);
    }
    // count views
    auto num_views = (size_t)0;
    for (auto& shape : scene.shapes) {
      if (!shape.positions.empty()) num_views += 1;
      if (!shape.normals.empty()) num_views += 1;
      if (!shape.texcoords.empty()) num_views += 1;
      if (!shape.colors.empty()) num_views += 1;
      if (!shape.radius.empty()) num_views += 1;
      if (!shape.points.empty() || !shape.lines.empty() ||
          !shape.triangles.empty() || !shape.quads.empty())
        num_views += 1;
    }
    alloc_arrays(gltf.buffers, gltf.buffers_count, scene.shapes.size());
    alloc_arrays(gltf.buffer_views, gltf.buffer_views_count, num_views);
    alloc_arrays(gltf.accessors, gltf.accessors_count, num_views);
    alloc_arrays(gltf.meshes, gltf.meshes_count, mesh_map.size());
    auto cur_view = (size_t)0, cur_mesh = (size_t)0, cur_accessor = (size_t)0;
    for (auto& shape : scene.shapes) {
      auto& materials = shape_materials[&shape - &scene.shapes.front()];
      auto& gbuffer   = gltf.buffers[&shape - &scene.shapes.front()];
      gbuffer.uri     = copy_string(
          "shapes/" + get_shape_name(scene, shape) + ".bin");
      for (auto mesh_id = cur_mesh; mesh_id < cur_mesh + materials.size();
           mesh_id++) {
        auto& gmesh = gltf.meshes[mesh_id];
        gmesh.name  = copy_string(get_shape_name(scene, shape));
        alloc_arrays(gmesh.primitives, gmesh.primitives_count, 1);
        auto& gprim    = gmesh.primitives[0];
        gprim.material = &gltf.materials[materials[mesh_id - cur_mesh]];
        if (!shape.positions.empty()) gprim.attributes_count += 1;
        if (!shape.normals.empty()) gprim.attributes_count += 1;
        if (!shape.texcoords.empty()) gprim.attributes_count += 1;
        if (!shape.colors.empty()) gprim.attributes_count += 1;
        if (!shape.radius.empty()) gprim.attributes_count += 1;
        alloc_arrays(
            gprim.attributes, gprim.attributes_count, gprim.attributes_count);
      }
      auto cur_attribute = (size_t)0;
      auto add_attribute = [&](const auto& values, const string& name) {
        auto& gview  = gltf.buffer_views[cur_view++];
        gview.buffer = &gbuffer;
        gview.size   = sizeof(values.front()) * values.size();
        gview.offset = gbuffer.size;
        gview.stride = 0;
        gview.type   = cgltf_buffer_view_type_vertices;
        gbuffer.size += gview.size;
        auto& gaccessor          = gltf.accessors[cur_accessor++];
        gaccessor.buffer_view    = &gview;
        gaccessor.component_type = cgltf_component_type_r_32f;
        if (sizeof(values.front()) == sizeof(vec3f)) {
          gaccessor.type = cgltf_type_vec3;
        } else if (sizeof(values.front()) == sizeof(vec2f)) {
          gaccessor.type = cgltf_type_vec2;
        } else if (sizeof(values.front()) == sizeof(vec4f)) {
          gaccessor.type = cgltf_type_vec4;
        } else if (sizeof(values.front()) == sizeof(float)) {
          gaccessor.type = cgltf_type_scalar;
        } else {
          // shoud not get here
          gaccessor.type = cgltf_type_scalar;
        }
        gaccessor.count = values.size();
        if constexpr (sizeof(values.front()) == sizeof(vec3f)) {
          if (name == "POSITION") {
            auto bbox = invalidb3f;
            for (auto& value : values) bbox = merge(bbox, value);
            gaccessor.has_min = true;
            gaccessor.has_max = true;
            gaccessor.min[0]  = bbox.min.x;
            gaccessor.min[1]  = bbox.min.y;
            gaccessor.min[2]  = bbox.min.z;
            gaccessor.max[0]  = bbox.max.x;
            gaccessor.max[1]  = bbox.max.y;
            gaccessor.max[2]  = bbox.max.z;
          }
        }
        for (auto mesh_id = cur_mesh; mesh_id < cur_mesh + materials.size();
             mesh_id++) {
          auto& gattribute =
              gltf.meshes[mesh_id].primitives[0].attributes[cur_attribute];
          gattribute.name  = copy_string(name);
          gattribute.data  = &gaccessor;
          gattribute.index = (int)cur_attribute;
        }
        cur_attribute += 1;
      };
      if (!shape.positions.empty()) add_attribute(shape.positions, "POSITION");
      if (!shape.normals.empty()) add_attribute(shape.normals, "NORMAL");
      if (!shape.texcoords.empty())
        add_attribute(shape.texcoords, "TEXCOORD_0");
      if (!shape.colors.empty()) add_attribute(shape.colors, "COLOR_0");
      if (!shape.radius.empty()) add_attribute(shape.radius, "RADIUS");
      auto add_indices = [&](const auto& values) {
        if (gltf.meshes[cur_mesh].primitives[0].indices != nullptr) return;
        auto& gview  = gltf.buffer_views[cur_view++];
        gview.buffer = &gbuffer;
        gview.size   = sizeof(values.front()) * values.size();
        gview.offset = gbuffer.size;
        gview.type   = cgltf_buffer_view_type_indices;
        gbuffer.size += gview.size;
        auto& gaccessor          = gltf.accessors[cur_accessor++];
        gaccessor.buffer_view    = &gview;
        gaccessor.component_type = cgltf_component_type_r_32u;
        gaccessor.type           = cgltf_type_scalar;
        gaccessor.count = values.size() * sizeof(values.front()) / sizeof(int);
        for (auto mesh_id = cur_mesh; mesh_id < cur_mesh + materials.size();
             mesh_id++) {
          gltf.meshes[mesh_id].primitives[0].indices = &gaccessor;
          if (sizeof(values.front()) == sizeof(int))
            gltf.meshes[mesh_id].primitives[0].type =
                cgltf_primitive_type_points;
          if (sizeof(values.front()) == sizeof(vec2i))
            gltf.meshes[mesh_id].primitives[0].type =
                cgltf_primitive_type_lines;
          if (sizeof(values.front()) == sizeof(vec3i))
            gltf.meshes[mesh_id].primitives[0].type =
                cgltf_primitive_type_triangles;
        }
      };
      if (!shape.points.empty()) add_indices(shape.points);
      if (!shape.lines.empty()) add_indices(shape.lines);
      if (!shape.triangles.empty()) add_indices(shape.triangles);
      if (!shape.quads.empty()) add_indices(quads_to_triangles(shape.quads));
      cur_mesh += materials.size();
    }
  }

  // nodes
  if (!scene.cameras.empty() || !scene.instances.empty()) {
    auto set_matrix = [](cgltf_node& gnode, const frame3f& frame) {
      auto mat         = frame_to_mat(frame);
      gnode.has_matrix = true;
      memcpy(gnode.matrix, &mat, sizeof(mat4f));
    };
    alloc_arrays(gltf.nodes, gltf.nodes_count,
        scene.cameras.size() + scene.instances.size() + 1);
    alloc_arrays(gltf.scenes, gltf.scenes_count, 1);
    auto cur_node = (size_t)0;
    for (auto& camera : scene.cameras) {
      auto& gnode  = gltf.nodes[cur_node++];
      gnode.name   = copy_string(get_camera_name(scene, camera));
      gnode.camera = &gltf.cameras[&camera - &scene.cameras.front()];
      set_matrix(gnode, camera.frame);
    }
    for (auto& instance : scene.instances) {
      auto& gnode = gltf.nodes[cur_node++];
      gnode.name  = copy_string(get_instance_name(scene, instance));
      gnode.mesh  = &gltf.meshes[mesh_map[{instance.shape, instance.material}]];
      set_matrix(gnode, instance.frame);
    }
    auto& groot = gltf.nodes[cur_node++];
    groot.name  = copy_string("root");
    alloc_arrays(groot.children, groot.children_count, gltf.nodes_count - 1);
    for (auto idx = 0; idx < gltf.nodes_count - 1; idx++) {
      groot.children[idx] = &gltf.nodes[idx];
    }
    auto& gscene = gltf.scenes[0];
    gscene.name  = copy_string("scene");
    alloc_arrays(gscene.nodes, gscene.nodes_count, 1);
    gscene.nodes[0] = &groot;
  }

  // save json
  auto goptions = cgltf_options{};
  memset(&goptions, 0, sizeof(goptions));
  goptions.type = cgltf_file_type_gltf;
  if (cgltf_write_file(&goptions, filename.c_str(), &gltf) !=
      cgltf_result_success)
    return write_error();

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // dirname
  auto dirname = path_dirname(filename);

  // save binary shape
  auto save_binshape = [](const string& filename, const scene_shape& shape,
                           string& error) {
    auto open_error = [filename, &error]() {
      error = filename + ": file not found";
      return false;
    };
    auto write_error = [filename, &error]() {
      error = filename + ": write error";
      return false;
    };

    auto fs = open_file(filename, "wb");
    if (!fs) return open_error();

    if (!write_values(fs, shape.positions)) return write_error();
    if (!write_values(fs, shape.normals)) return write_error();
    if (!write_values(fs, shape.texcoords)) return write_error();
    if (!write_values(fs, shape.colors)) return write_error();
    if (!write_values(fs, shape.radius)) return write_error();
    if (!write_values(fs, shape.points)) return write_error();
    if (!write_values(fs, shape.lines)) return write_error();
    if (!write_values(fs, shape.triangles)) return write_error();
    if (!write_values(fs, quads_to_triangles(shape.quads)))
      return write_error();

    return true;
  };

  // save shapes
  for (auto& shape : scene.shapes) {
    if (progress_cb) progress_cb("save buffer", progress.x++, progress.y);
    if (!save_binshape(path_join(dirname,
                           "shapes/" + get_shape_name(scene, shape) + ".bin"),
            shape, error))
      return dependent_error();
  }

  // save textures
  for (auto& texture : scene.textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    if (!save_image(
            path_join(dirname,
                "textures/" + get_texture_name(scene, texture) + ".png"),
            texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save done", progress.x++, progress.y);
  return true;
}

#else

// Load a scene
static bool save_gltf_scene(const string& filename, const scene_scene& scene,
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

  // convert scene to json
  auto js = njson::object();

  // asset
  {
    auto& ajs        = js["asset"];
    ajs              = njson::object();
    ajs["version"]   = "2.0";
    ajs["generator"] = scene.asset.generator;
    ajs["copyright"] = scene.asset.copyright;
  }

  // cameras
  if (!scene.cameras.empty()) {
    auto& ajs = js["cameras"];
    ajs       = njson::array();
    for (auto& camera : scene.cameras) {
      auto& cjs          = ajs.emplace_back();
      cjs                = njson::object();
      cjs["name"]        = get_camera_name(scene, camera);
      cjs["type"]        = "perspective";
      auto& pjs          = cjs["perspective"];
      pjs                = njson::object();
      pjs["aspectRatio"] = camera.aspect;
      pjs["yfov"]        = 0.660593;  // TODO(fabio): yfov
      pjs["znear"]       = 0.001;     // TODO(fabio): configurable?
    }
  }

  // materials
  auto textures    = vector<pair<string, image<vec4b>>>{};
  auto texture_map = unordered_map<string, int>{};
  if (!scene.materials.empty()) {
    auto& ajs = js["materials"];
    ajs       = njson::array();
    for (auto& material : scene.materials) {
      auto& mjs              = ajs.emplace_back();
      mjs                    = njson::object();
      mjs["name"]            = get_material_name(scene, material);
      mjs["emissiveFactor"]  = material.emission;
      auto& pjs              = mjs["pbrMetallicRoughness"];
      pjs                    = njson::object();
      pjs["baseColorFactor"] = vec4f{material.color.x, material.color.y,
          material.color.z, material.opacity};
      pjs["metallicFactor"]  = material.metallic;
      pjs["roughnessFactor"] = material.roughness;
      if (material.emission_tex != invalid_handle) {
        auto tname = get_texture_name(
            scene, material.emission_tex);  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          auto& texture = scene.textures[material.emission_tex];
          textures.emplace_back(tname, texture.ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        mjs["emissiveTexture"]          = njson::object();
        mjs["emissiveTexture"]["index"] = texture_map.at(tname);
      }
      if (material.normal_tex != invalid_handle) {
        auto tname = get_texture_name(
            scene, material.normal_tex);  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          auto& texture = scene.textures[material.normal_tex];
          textures.emplace_back(tname, texture.ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        mjs["normalTexture"]          = njson::object();
        mjs["normalTexture"]["index"] = texture_map.at(tname);
      }
      if (material.color_tex != invalid_handle) {  // TODO(fabio): opacity
        auto tname = get_texture_name(
            scene, material.color_tex);  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          auto& texture = scene.textures[material.color_tex];
          textures.emplace_back(tname, texture.ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        pjs["baseColorTexture"]          = njson::object();
        pjs["baseColorTexture"]["index"] = texture_map.at(tname);
      }
      if (material.roughness_tex != invalid_handle) {  // TODO(fabio): roughness
        auto tname = get_texture_name(
            scene, material.roughness_tex);  // TODO(fabio): ldr
        if (texture_map.find(tname) == texture_map.end()) {
          auto& texture = scene.textures[material.roughness_tex];
          textures.emplace_back(tname, texture.ldr);
          texture_map[tname] = (int)textures.size() - 1;
        }
        pjs["metallicRoughnessTexture"]          = njson::object();
        pjs["metallicRoughnessTexture"]["index"] = texture_map.at(tname);
      }
    }
  }

  // textures
  if (!textures.empty()) {
    js["textures"] = njson::array();
    js["samplers"] = njson::array();
    js["images"]   = njson::array();
    auto& sjs      = js["samplers"].emplace_back();
    sjs            = njson::object();
    sjs["name"]    = "sampler";
    for (auto& [name, img] : textures) {
      auto& ijs      = js["images"].emplace_back();
      ijs            = njson::object();
      ijs["name"]    = name;
      ijs["uri"]     = "textures/" + name + ".png";
      auto& tjs      = js["textures"].emplace_back();
      tjs            = njson::object();
      tjs["name"]    = name;
      tjs["sampler"] = 0;
      tjs["source"]  = (int)js["images"].size() - 1;
    }
  }

  // add an accessor
  auto add_accessor = [](njson& js, vector<pair<string, vector<byte>>>& buffers,
                          const void* data, size_t count, size_t size,
                          bool is_index = false) -> int {
    static auto types = unordered_map<size_t, string>{
        {1, "SCALAR"}, {2, "VEC2"}, {3, "VEC3"}, {4, "VEC4"}};
    auto  length         = count * size * 4;
    auto& vjs            = js["bufferViews"].emplace_back();
    vjs                  = njson::object();
    vjs["buffer"]        = (int)buffers.size() - 1;
    vjs["byteLength"]    = (uint64_t)length;
    vjs["byteOffset"]    = (uint64_t)buffers.back().second.size();
    vjs["target"]        = is_index ? 34963 : 34962;
    auto& ajs            = js["accessors"].emplace_back();
    ajs                  = njson::object();
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
      ajs["min"] = min_;
      ajs["max"] = max_;
    }
    buffers.back().second.insert(
        buffers.back().second.end(), (byte*)data, (byte*)data + length);
    return (int)js["accessors"].size() - 1;
  };

  // meshes
  auto buffers        = vector<pair<string, vector<byte>>>{};
  auto primitives_map = unordered_map<shape_handle, njson>{};
  if (!scene.shapes.empty()) {
    js["accessors"]   = njson::array();
    js["bufferViews"] = njson::array();
    js["buffers"]     = njson::array();
    auto shape_id     = 0;
    for (auto& shape : scene.shapes) {
      auto& buffer =
          buffers.emplace_back(get_shape_name(scene, shape), vector<byte>{})
              .second;
      auto& pjs = primitives_map[shape_id++];
      pjs       = njson::object();
      auto& ajs = pjs["attributes"];
      ajs       = njson::object();
      if (!shape.positions.empty()) {
        ajs["POSITION"] = add_accessor(
            js, buffers, shape.positions.data(), shape.positions.size(), 3);
      }
      if (!shape.normals.empty()) {
        ajs["NORMAL"] = add_accessor(
            js, buffers, shape.normals.data(), shape.normals.size(), 3);
      }
      if (!shape.texcoords.empty()) {
        ajs["TEXCOORD_0"] = add_accessor(
            js, buffers, shape.texcoords.data(), shape.texcoords.size(), 2);
      }
      if (!shape.colors.empty()) {
        ajs["COLOR_0"] = add_accessor(
            js, buffers, shape.colors.data(), shape.colors.size(), 4);
      }
      if (!shape.radius.empty()) {
        ajs["_RADIUS"] = add_accessor(
            js, buffers, shape.radius.data(), shape.radius.size(), 1);
      }
      if (!shape.points.empty()) {
        pjs["indices"] = add_accessor(
            js, buffers, shape.points.data(), shape.points.size(), 1, true);
        pjs["mode"] = 0;
      } else if (!shape.lines.empty()) {
        pjs["indices"] = add_accessor(
            js, buffers, shape.lines.data(), shape.lines.size() * 2, 1, true);
        pjs["mode"] = 1;
      } else if (!shape.triangles.empty()) {
        pjs["indices"] = add_accessor(js, buffers, shape.triangles.data(),
            shape.triangles.size() * 3, 1, true);
        pjs["mode"]    = 4;
      } else if (!shape.quads.empty()) {
        auto triangles = quads_to_triangles(shape.quads);
        pjs["indices"] = add_accessor(
            js, buffers, triangles.data(), triangles.size() * 3, 1, true);
        pjs["mode"] = 4;
      }
      auto& bjs         = js["buffers"].emplace_back();
      bjs               = njson::object();
      bjs["byteLength"] = (uint64_t)buffer.size();
      bjs["uri"]        = "shapes/" + get_shape_name(scene, shape) + ".bin";
    }
  }

  // nodes
  js["nodes"] = njson::array();
  if (!scene.cameras.empty()) {
    for (auto idx = 0; idx < (int)scene.cameras.size(); idx++) {
      auto& camera = scene.cameras[idx];
      auto& njs    = js["nodes"].emplace_back();
      njs          = njson::object();
      njs["name"]  = scene.camera_names[idx];
      auto matrix  = frame_to_mat(camera.frame);  // TODO(fabio): do this better
      njs["matrix"] = matrix;
      njs["camera"] = idx;
    }
  }
  if (!scene.instances.empty()) {
    js["meshes"]   = njson::array();
    using mesh_key = pair<shape_handle, material_handle>;
    struct mesh_key_hash {
      size_t operator()(const mesh_key& v) const {
        const std::hash<element_handle> hasher = std::hash<element_handle>();
        auto                            h      = (size_t)0;
        h ^= hasher(v.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
        h ^= hasher(v.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
      }
    };
    auto mesh_map = unordered_map<mesh_key, int, mesh_key_hash>{};
    for (auto& instance : scene.instances) {
      auto& njs     = js["nodes"].emplace_back();
      njs           = njson::object();
      njs["name"]   = get_instance_name(scene, instance);
      njs["matrix"] = frame_to_mat(instance.frame);
      if (mesh_map.find(mesh_key{instance.shape, instance.material}) ==
          mesh_map.end()) {
        auto& mjs   = js["meshes"].emplace_back();
        mjs         = njson::object();
        mjs["name"] = get_shape_name(scene, instance.shape) + "_" +
                      get_material_name(scene, instance.material);
        mjs["primitives"] = njson::array();
        mjs["primitives"].push_back(primitives_map.at(instance.shape));
        mjs["primitives"].back()["material"] = instance.material;
        mesh_map[mesh_key{instance.shape, instance.material}] =
            (int)js["meshes"].size() - 1;
      }
      njs["mesh"] = mesh_map.at({instance.shape, instance.material});
    }
  } else {
    js["meshes"] = njson::array();
    for (auto& [shape, pjs] : primitives_map) {
      auto& mjs         = js["meshes"].emplace_back();
      mjs               = njson::object();
      mjs["name"]       = get_shape_name(scene, shape);
      mjs["primitives"] = njson::array();
      mjs["primitives"].push_back(pjs);
      auto& njs   = js["nodes"].emplace_back();
      njs         = njson::object();
      njs["name"] = get_shape_name(scene, shape);
      njs["mesh"] = (int)js["meshes"].size() - 1;
    }
  }

  // root children
  {
    auto& rjs       = js["nodes"].emplace_back();
    rjs             = njson::object();
    rjs["name"]     = "root";
    rjs["children"] = njson::array();
    for (auto idx = 0; idx < (int)js["nodes"].size() - 1; idx++)
      rjs["children"].push_back(idx);
  }

  // scene
  {
    js["scenes"] = njson::array();
    auto& sjs    = js["scenes"].emplace_back();
    sjs          = njson::object();
    sjs["nodes"] = njson::array();
    sjs["nodes"].push_back((int)js["nodes"].size() - 1);
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

#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PBRT
// -----------------------------------------------------------------------------
namespace yocto {

// load pbrt scenes
static bool load_pbrt_scene(const string& filename, scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load pbrt
  auto pbrt = pbrt_scene{};
  if (!load_pbrt(filename, pbrt, error)) return false;

  // handle progress
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = pcamera.lens;
    camera.focus  = pcamera.focus;
  }

  // convert materials
  auto texture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_texture = [&scene, &texture_map](
                         const string& path) -> texture_handle {
    if (path.empty()) return invalid_handle;
    auto it = texture_map.find(path);
    if (it != texture_map.end()) return it->second;
    scene.textures.emplace_back();
    texture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };
  auto atexture_map = unordered_map<string, texture_handle>{
      {"", invalid_handle}};
  auto get_atexture = [&scene, &atexture_map](
                          const string& path) -> texture_handle {
    if (path.empty()) return invalid_handle;
    auto it = atexture_map.find(path);
    if (it != atexture_map.end()) return it->second;
    scene.textures.emplace_back();
    atexture_map[path] = (int)scene.textures.size() - 1;
    return (int)scene.textures.size() - 1;
  };

  // material type map
  auto material_type_map = unordered_map<pbrt_material_type, material_type>{
      {pbrt_material_type::matte, material_type::matte},
      {pbrt_material_type::plastic, material_type::plastic},
      {pbrt_material_type::metal, material_type::metal},
      {pbrt_material_type::glass, material_type::glass},
      {pbrt_material_type::thinglass, material_type::thinglass},
      {pbrt_material_type::subsurface, material_type::matte},
  };

  // convert material
  auto material_map = unordered_map<string, material_handle>{};
  for (auto& pmaterial : pbrt.materials) {
    auto& material = scene.materials.emplace_back();
    material.type  = material_type_map.at(pmaterial.type);
    if (pmaterial.emission != zero3f) {
      material.type = material_type::matte;
    }
    material.emission  = pmaterial.emission;
    material.color     = pmaterial.color;
    material.ior       = pmaterial.ior;
    material.roughness = pmaterial.roughness;
    material.opacity   = pmaterial.opacity;
    material.color_tex = get_texture(pmaterial.color_tex);
    // material.opacity_tex  = get_texture(pmaterial.opacity_tex);
    // if (material.opacity_tex == invalid_handle)
    //   material.opacity_tex = get_atexture(pmaterial.alpha_tex);
    material_map[pmaterial.name] = (int)scene.materials.size() - 1;
  }

  // hack for pbrt empty material
  scene.materials.emplace_back();
  material_map[""] = (int)scene.materials.size() - 1;

  // convert shapes
  for (auto& pshape : pbrt.shapes) {
    auto& shape     = scene.shapes.emplace_back();
    shape.positions = pshape.positions;
    shape.normals   = pshape.normals;
    shape.texcoords = pshape.texcoords;
    shape.triangles = pshape.triangles;
    for (auto& uv : shape.texcoords) uv.y = 1 - uv.y;
    auto material = material_map.at(pshape.material);
    if (pshape.instances.empty()) {
      auto& instance    = scene.instances.emplace_back();
      instance.frame    = pshape.frame;
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = material;
    } else {
      for (auto frame : pshape.instances) {
        auto& instance    = scene.instances.emplace_back();
        instance.frame    = frame * pshape.frame;
        instance.shape    = (int)scene.shapes.size() - 1;
        instance.material = material;
      }
    }
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = penvironment.frame;
    environment.emission     = penvironment.emission;
    environment.emission_tex = get_texture(penvironment.emission_tex);
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape       = scene.shapes.emplace_back();
    shape.triangles   = plight.area_triangles;
    shape.positions   = plight.area_positions;
    shape.normals     = plight.area_normals;
    auto& material    = scene.materials.emplace_back();
    material.emission = plight.area_emission;
    auto& instance    = scene.instances.emplace_back();
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = (int)scene.materials.size() - 1;
    instance.frame    = plight.area_frame;
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name) {
    return path_join(path_dirname(filename), name);
  };

  // load texture
  texture_map.erase("");
  for (auto [name, thandle] : texture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // load alpha
  atexture_map.erase("");
  for (auto [name, thandle] : atexture_map) {
    auto& texture = scene.textures[thandle];
    if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
    if (!load_image(make_filename(name), texture.hdr, texture.ldr, error))
      return dependent_error();
    for (auto& c : texture.hdr) {
      c = (max(vec3f{c.x, c.y, c.z}) < 0.01) ? vec4f{0, 0, 0, c.w}
                                             : vec4f{1, 1, 1, c.w};
    }
    for (auto& c : texture.ldr) {
      c = (max(vec3i{c.x, c.y, c.z}) < 2) ? vec4b{0, 0, 0, c.w}
                                          : vec4b{255, 255, 255, c.w};
    }
  }

  // fix scene
  if (scene.asset.name.empty()) scene.asset.name = path_basename(filename);
  add_cameras(scene);
  add_radius(scene);
  add_materials(scene);

  // done
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);
  return true;
}

// Save a pbrt scene
static bool save_pbrt_scene(const string& filename, const scene_scene& scene,
    string& error, const progress_callback& progress_cb, bool noparallel) {
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // handle progress
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt = pbrt_scene{};

  // convert camera
  auto& camera       = scene.cameras.front();
  auto& pcamera      = add_camera(pbrt);
  pcamera.frame      = camera.frame;
  pcamera.lens       = camera.lens;
  pcamera.aspect     = camera.aspect;
  pcamera.resolution = {1280, (int)(1280 / pcamera.aspect)};

  // get texture name
  auto get_texture = [&](texture_handle texture) -> string {
    return texture != invalid_handle ? get_texture_name(scene, texture) : "";
  };

  // material type map
  auto material_type_map = unordered_map<material_type, pbrt_material_type>{
      {material_type::matte, pbrt_material_type::matte},
      {material_type::plastic, pbrt_material_type::plastic},
      {material_type::metal, pbrt_material_type::metal},
      {material_type::glass, pbrt_material_type::glass},
      {material_type::thinglass, pbrt_material_type::thinglass},
      {material_type::subsurface, pbrt_material_type::matte},
      {material_type::volume, pbrt_material_type::matte},
  };

  // convert materials
  auto material_map = unordered_map<material_handle, string>{};
  auto material_id  = 0;
  for (auto& material : scene.materials) {
    auto& pmaterial             = add_material(pbrt);
    pmaterial.name              = get_material_name(scene, material);
    pmaterial.type              = material_type_map.at(material.type);
    pmaterial.emission          = material.emission;
    pmaterial.color             = material.color;
    pmaterial.roughness         = material.roughness;
    pmaterial.ior               = material.ior;
    pmaterial.opacity           = material.opacity;
    pmaterial.color_tex         = get_texture(material.color_tex);
    material_map[material_id++] = pmaterial.name;
  }

  // convert instances
  for (auto& instance : scene.instances) {
    auto& pshape     = add_shape(pbrt);
    pshape.filename_ = get_shape_name(scene, instance.shape) + ".ply";
    pshape.frame     = instance.frame;
    pshape.frend     = instance.frame;
    pshape.material  = material_map.at(instance.material);
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment        = add_environment(pbrt);
    penvironment.emission     = environment.emission;
    penvironment.emission_tex = get_texture(environment.emission_tex);
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  if (!save_pbrt(filename, pbrt, error)) return false;

  // handle progress
  progress.y += (int)scene.shapes.size() + (int)scene.textures.size();

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // save textures
  auto texture_id = 0;
  for (auto& texture : scene.textures) {
    if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
    auto path = make_filename(get_texture_name(scene, texture_id++), "textures",
        (!texture.hdr.empty()) ? ".hdr"s : ".png"s);
    if (!save_image(path, texture.hdr, texture.ldr, error))
      return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto
