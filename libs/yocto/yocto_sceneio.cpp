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

// load texture
static bool load_texture(
    const string& filename, scene_texture& texture, string& error) {
  if (is_hdr_filename(filename)) {
    return load_image(filename, texture.hdr, error);
  } else {
    return load_image(filename, texture.ldr, error);
  }
}

// save texture
static bool save_texture(
    const string& filename, const scene_texture& texture, string& error) {
  if (!texture.hdr.empty()) {
    return save_image(filename, texture.hdr, error);
  } else {
    return save_image(filename, texture.ldr, error);
  }
}

// load shape
static bool load_shape(
    const string& filename, scene_shape& shape, string& error) {
  auto lshape = shape_data{};
  if (!load_shape(filename, lshape, error, true)) return false;
  shape.points    = lshape.points;
  shape.lines     = lshape.lines;
  shape.triangles = lshape.triangles;
  shape.quads     = lshape.quads;
  shape.positions = lshape.positions;
  shape.normals   = lshape.normals;
  shape.texcoords = lshape.texcoords;
  shape.colors    = lshape.colors;
  shape.radius    = lshape.radius;
  return true;
}

// save shape
static bool save_shape(
    const string& filename, const scene_shape& shape, string& error) {
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
  return save_shape(filename, sshape, error, true);
}

// load subdiv
static bool load_subdiv(
    const string& filename, scene_subdiv& subdiv, string& error) {
  auto lsubdiv = fvshape_data{};
  if (!load_fvshape(filename, lsubdiv, error, true)) return false;
  subdiv.quadspos      = lsubdiv.quadspos;
  subdiv.quadsnorm     = lsubdiv.quadsnorm;
  subdiv.quadstexcoord = lsubdiv.quadstexcoord;
  subdiv.positions     = lsubdiv.positions;
  subdiv.normals       = lsubdiv.normals;
  subdiv.texcoords     = lsubdiv.texcoords;
  return true;
}

// save subdiv
static bool save_subdiv(
    const string& filename, const scene_subdiv& subdiv, string& error) {
  auto ssubdiv          = fvshape_data{};
  ssubdiv.quadspos      = subdiv.quadspos;
  ssubdiv.quadsnorm     = subdiv.quadsnorm;
  ssubdiv.quadstexcoord = subdiv.quadstexcoord;
  ssubdiv.positions     = subdiv.positions;
  ssubdiv.normals       = subdiv.normals;
  ssubdiv.texcoords     = subdiv.texcoords;
  return save_fvshape(filename, ssubdiv, error, true);
}

// save binary shape
static bool save_binshape(
    const string& filename, const scene_shape& shape, string& error) {
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
  if (!write_values(fs, quads_to_triangles(shape.quads))) return write_error();

  return true;
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
  auto material_error = [filename, &error](const string& name) {
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
  auto shape_map = unordered_map<string, shape_handle>{};
  auto get_shape = [&scene, &shape_map, &get_value](
                       const njson& js, shape_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = shape_map.find(name);
    if (it != shape_map.end()) {
      value = it->second;
      return true;
    }
    scene.shape_names.emplace_back(name);
    scene.shapes.emplace_back();
    auto shape_id   = (int)scene.shapes.size() - 1;
    shape_map[name] = shape_id;
    value           = shape_id;
    return true;
  };

  // parse json reference
  auto material_map = unordered_map<string, material_handle>{};
  auto material_set = vector<bool>{};
  auto get_material = [&scene, &material_map, &material_set, &get_value](
                          const njson& js, material_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = material_map.find(name);
    if (it != material_map.end()) {
      value = it->second;
      return true;
    }
    scene.material_names.emplace_back(name);
    scene.materials.emplace_back();
    auto material_id   = (int)scene.materials.size() - 1;
    material_map[name] = material_id;
    value              = material_id;
    material_set.push_back(false);
    return true;
  };

  // parse json reference
  auto texture_map = unordered_map<string, texture_handle>{};
  auto get_texture = [&scene, &texture_map, &get_value](
                         const njson& js, texture_handle& value) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    auto it = texture_map.find(name);
    if (it != texture_map.end()) {
      value = it->second;
      return true;
    }
    scene.texture_names.emplace_back(name);
    scene.textures.emplace_back();
    auto texture_id   = (int)scene.textures.size() - 1;
    texture_map[name] = texture_id;
    value             = texture_id;
    return true;
  };

  // load json instance
  struct ply_instance {
    vector<frame3f> frames = {};
  };
  using ply_instance_handle = int;
  auto ply_instances        = vector<ply_instance>{};
  auto ply_instances_names  = vector<string>{};
  auto ply_instance_map     = unordered_map<string, ply_instance_handle>{
      {"", invalid_handle}};
  auto instance_ply = unordered_map<instance_handle, ply_instance_handle>{};
  auto get_ply_instances = [&scene, &ply_instances, &ply_instances_names,
                               &ply_instance_map, &instance_ply,
                               &get_value](const njson& js,
                               const scene_instance&    instance) -> bool {
    auto name = ""s;
    if (!get_value(js, name)) return false;
    if (name.empty()) return true;
    auto instance_id = (int)(&instance - scene.instances.data());
    auto it          = ply_instance_map.find(name);
    if (it != ply_instance_map.end()) {
      instance_ply[instance_id] = it->second;
      return true;
    }
    ply_instances_names.emplace_back(name);
    ply_instances.emplace_back(ply_instance());
    auto ply_instance_id      = (int)ply_instances.size() - 1;
    ply_instance_map[name]    = ply_instance_id;
    instance_ply[instance_id] = ply_instance_id;
    return true;
  };
  auto get_ply_instance_name = [&ply_instances, &ply_instances_names](
                                   const scene_scene&  scene,
                                   const ply_instance& instance) -> string {
    return ply_instances_names[&instance - ply_instances.data()];
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
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

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
        if (material_map.find(name) == material_map.end()) {
          scene.material_names.emplace_back(name);
          scene.materials.emplace_back();
          material_map[name] = (int)scene.materials.size() - 1;
          material_set.push_back(false);
        }
        auto& material = scene.materials.at(material_map.at(name));
        material_set[&material - &scene.materials.front()] = true;
        material.type = material_type::metallic;
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
            if (!get_ply_instances(value, instance))
              return parse_error(gname, name, key);
          } else {
            return key_error(gname, name, key);
          }
        }
      }
    } else if (gname == "subdivs") {
      if (!check_object(group)) return parse_error(gname);
      for (auto& [name, element] : iterate_object(group)) {
        if (!check_object(element)) return parse_error(gname, name);
        scene.subdiv_names.emplace_back(name);
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
      }
    } else {
      return key_error(gname);
    }
  }

  // check materials
  for (auto& material : scene.materials) {
    if (!material_set[&material - &scene.materials.front()])
      return material_error(get_material_name(scene, material));
  }

  // handle progress
  progress.y += scene.shapes.size();
  progress.y += scene.textures.size();
  progress.y += scene.subdivs.size();
  progress.y += ply_instances.size();

  // dirname
  auto dirname = path_dirname(filename);

  // get filename from name
  auto find_path = [dirname](const string& name, const string& group,
                       const vector<string>& extensions) {
    for (auto& extension : extensions) {
      auto path = path_join(dirname, group, name + extension);
      if (path_exists(path)) return path_join(group, name + extension);
    }
    return path_join(group, name + extensions.front());
  };

  // load resources
  if (noparallel) {
    // load shapes
    for (auto& shape : scene.shapes) {
      if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
      auto path = find_path(
          get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
      if (!load_shape(path_join(dirname, path), shape, error))
        return dependent_error();
    }
    // load subdivs
    for (auto& subdiv : scene.subdivs) {
      if (progress_cb) progress_cb("load subdiv", progress.x++, progress.y);
      auto path = find_path(
          get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
      if (!load_subdiv(path_join(dirname, path), subdiv, error))
        return dependent_error();
    }
    // load textures
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      auto path = find_path(get_texture_name(scene, texture), "textures",
          {".hdr", ".exr", ".png", ".jpg"});
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
    // load instances
    for (auto& ply_instance : ply_instances) {
      if (progress_cb) progress_cb("load instances", progress.x++, progress.y);
      auto path = find_path(
          get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
      if (!load_instance(path_join(dirname, path), ply_instance.frames, error))
        return dependent_error();
    }
  } else {
    // helpers
    auto mutex = std::mutex{};
    // load shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
      }
      auto path = find_path(
          get_shape_name(scene, shape), "shapes", {".ply", ".obj"});
      auto err = string{};
      if (!load_shape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load subdivs
    parallel_foreach(scene.subdivs, [&](auto& subdiv) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load subdiv", progress.x++, progress.y);
      }
      auto path = find_path(
          get_subdiv_name(scene, subdiv), "subdivs", {".ply", ".obj"});
      auto err = string{};
      if (!load_subdiv(path_join(dirname, path), subdiv, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      }
      auto path = find_path(get_texture_name(scene, texture), "textures",
          {".hdr", ".exr", ".png", ".jpg"});
      auto err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load instances
    parallel_foreach(ply_instances, [&](auto& ply_instance) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb)
          progress_cb("load instances", progress.x++, progress.y);
      }
      auto path = find_path(
          get_ply_instance_name(scene, ply_instance), "instances", {".ply"});
      auto err = string{};
      if (!load_instance(path_join(dirname, path), ply_instance.frames, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // apply instances
  if (!ply_instances.empty()) {
    if (progress_cb)
      progress_cb("flatten instances", progress.x++, progress.y++);
    auto instances      = scene.instances;
    auto instance_names = scene.instance_names;
    scene.instances.clear();
    scene.instance_names.clear();
    for (auto& instance : instances) {
      auto it = instance_ply.find((int)(&instance - instances.data()));
      if (it == instance_ply.end()) {
        auto& ninstance = scene.instances.emplace_back();
        scene.instance_names.emplace_back(
            instance_names[&instance - instances.data()]);
        ninstance.frame    = instance.frame;
        ninstance.shape    = instance.shape;
        ninstance.material = instance.material;
      } else {
        auto& ply_instance = ply_instances[it->second];
        auto  instance_id  = 0;
        for (auto& frame : ply_instance.frames) {
          auto& ninstance = scene.instances.emplace_back();
          scene.instance_names.emplace_back(
              instance_names[&instance - instances.data()] + "_" +
              std::to_string(instance_id++));
          ninstance.frame    = frame * instance.frame;
          ninstance.shape    = instance.shape;
          ninstance.material = instance.material;
        }
      }
    }
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
  auto progress = vec2i{0, 2 + (int)scene.shapes.size() +
                               (int)scene.subdivs.size() +
                               (int)scene.textures.size()};
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

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

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save shapes
    for (auto& shape : scene.shapes) {
      if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      if (!save_shape(path_join(dirname, path), shape, error))
        return dependent_error();
    }

    // save subdiv
    for (auto& subdiv : scene.subdivs) {
      if (progress_cb) progress_cb("save subdiv", progress.x++, progress.y);
      auto path = "subdivs/" + get_subdiv_name(scene, subdiv) + ".obj";
      if (!save_subdiv(path_join(dirname, path), subdiv, error))
        return dependent_error();
    }

    // save textures
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      auto err  = string{};
      if (!save_shape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save subdivs
    parallel_foreach(scene.subdivs, [&](auto& subdiv) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save subdiv", progress.x++, progress.y);
      }
      auto path = "subdivs/" + get_subdiv_name(scene, subdiv) + ".obj";
      auto err  = string{};
      if (!save_subdiv(path_join(dirname, path), subdiv, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

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
  if (!load_obj(filename, obj, error, false, true)) return false;

  // handle progress
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // convert cameras
  for (auto& ocamera : obj.cameras) {
    auto& camera        = scene.cameras.emplace_back();
    camera.frame        = ocamera.frame;
    camera.orthographic = ocamera.ortho;
    camera.film         = ocamera.film;
    camera.aspect       = ocamera.aspect;
    camera.focus        = ocamera.focus;
    camera.lens         = ocamera.lens;
    camera.aperture     = ocamera.aperture;
  }

  // convert between roughness and exponent
  auto exponent_to_roughness = [](float exponent) {
    auto roughness = exponent;
    roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
    if (roughness < 0.01f) roughness = 0;
    if (roughness > 0.99f) roughness = 1;
    return roughness;
  };

  // handler for textures
  auto textures_paths = vector<string>{};
  for (auto& otexture : obj.textures) {
    scene.textures.emplace_back();
    textures_paths.emplace_back(otexture.path);
  }

  // handler for materials
  for (auto& omaterial : obj.materials) {
    auto& material        = scene.materials.emplace_back();
    material.type         = material_type::metallic;
    material.emission     = omaterial.emission;
    material.emission_tex = omaterial.emission_tex;
    if (max(omaterial.transmission) > 0.1) {
      material.type      = material_type::thinglass;
      material.color     = omaterial.transmission;
      material.color_tex = omaterial.transmission_tex;
    } else if (max(omaterial.specular) > 0.2) {
      material.type      = material_type::metal;
      material.color     = omaterial.specular;
      material.color_tex = omaterial.specular_tex;
    } else if (max(omaterial.specular) > 0) {
      material.type      = material_type::plastic;
      material.color     = omaterial.diffuse;
      material.color_tex = omaterial.diffuse_tex;
    } else {
      material.type      = material_type::matte;
      material.color     = omaterial.diffuse;
      material.color_tex = omaterial.diffuse_tex;
    }
    material.roughness  = exponent_to_roughness(omaterial.exponent);
    material.ior        = omaterial.ior;
    material.metallic   = 0;
    material.opacity    = omaterial.opacity;
    material.normal_tex = omaterial.normal_tex;
  }

  // convert shapes
  for (auto& oshape : obj.shapes) {
    if (oshape.elements.empty()) continue;
    auto& shape       = scene.shapes.emplace_back();
    auto& instance    = scene.instances.emplace_back();
    instance.shape    = (int)scene.shapes.size() - 1;
    instance.material = oshape.elements.front().material;
    get_positions(oshape, shape.positions);
    get_normals(oshape, shape.normals);
    get_texcoords(oshape, shape.texcoords, true);
    get_faces(oshape, instance.material, shape.triangles, shape.quads);
    get_lines(oshape, instance.material, shape.lines);
    get_points(oshape, instance.material, shape.points);
  }

  // convert environments
  for (auto& oenvironment : obj.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = oenvironment.frame;
    environment.emission     = oenvironment.emission;
    environment.emission_tex = oenvironment.emission_tex;
  }

  // handle progress
  progress.y += (int)scene.textures.size();

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // load textures
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      auto& path = textures_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      }
      auto& path = textures_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // build obj
  auto obj = obj_scene{};

  // convert cameras
  for (auto& camera : scene.cameras) {
    auto& ocamera    = obj.cameras.emplace_back();
    ocamera.name     = get_camera_name(scene, camera);
    ocamera.frame    = camera.frame;
    ocamera.ortho    = camera.orthographic;
    ocamera.film     = camera.film;
    ocamera.aspect   = camera.aspect;
    ocamera.focus    = camera.focus;
    ocamera.lens     = camera.lens;
    ocamera.aperture = camera.aperture;
  }

  // helper
  auto roughness_to_exponent = [](float roughness) -> float {
    if (roughness < 0.01f) return 10000;
    if (roughness > 0.99f) return 10;
    return 2 / pow(roughness, 4.0f) - 2;
  };

  // convert textures
  for (auto& texture : scene.textures) {
    auto& otexture = obj.textures.emplace_back();
    otexture.path  = "textures/" + get_texture_name(scene, texture) +
                    (!texture.hdr.empty() ? ".hdr"s : ".png"s);
  }

  // convert materials
  for (auto& material : scene.materials) {
    auto& omaterial        = obj.materials.emplace_back();
    omaterial.name         = get_material_name(scene, material);
    omaterial.illum        = 2;
    omaterial.emission     = material.emission;
    omaterial.diffuse      = material.color;
    omaterial.specular     = {0, 0, 0};
    omaterial.exponent     = roughness_to_exponent(material.roughness);
    omaterial.opacity      = material.opacity;
    omaterial.emission_tex = material.emission_tex;
    omaterial.diffuse_tex  = material.color_tex;
    omaterial.normal_tex   = material.normal_tex;
  }

  // convert objects
  for (auto& instance : scene.instances) {
    auto& shape     = scene.shapes[instance.shape];
    auto  positions = shape.positions, normals = shape.normals;
    for (auto& p : positions) p = transform_point(instance.frame, p);
    for (auto& n : normals) n = transform_normal(instance.frame, n);
    auto& oshape = obj.shapes.emplace_back();
    oshape.name  = get_shape_name(scene, shape);
    add_positions(oshape, positions);
    add_normals(oshape, normals);
    add_texcoords(oshape, shape.texcoords, true);
    add_triangles(oshape, shape.triangles, instance.material,
        !shape.normals.empty(), !shape.texcoords.empty());
    add_quads(oshape, shape.quads, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_lines(oshape, shape.lines, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
    add_points(oshape, shape.points, instance.material, !shape.normals.empty(),
        !shape.texcoords.empty());
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& oenvironment        = obj.environments.emplace_back();
    oenvironment.name         = get_environment_name(scene, environment);
    oenvironment.frame        = environment.frame;
    oenvironment.emission     = environment.emission;
    oenvironment.emission_tex = environment.emission_tex;
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save obj
  if (!save_obj(filename, obj, error)) return false;

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save textures
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
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
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace yocto {

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
  auto progress = vec2i{0, 2};
  if (progress_cb) progress_cb("load scene", progress.x++, progress.y);

  // load gltf
  auto gltf = njson{};
  if (!load_json(filename, gltf, error)) return false;

  // parse buffers
  auto buffers_paths = vector<string>{};
  auto buffers       = vector<vector<byte>>();
  try {
    if (gltf.contains("buffers")) {
      for (auto& gbuffer : gltf.at("buffers")) {
        if (!gbuffer.contains("uri")) return parse_error();
        buffers_paths.push_back(gbuffer.value("uri", ""));
        buffers.emplace_back();
      }
    }
  } catch (...) {
    return parse_error();
  }

  // handle progress
  progress.y += (int)buffers_paths.size();

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // load buffers
    for (auto& buffer : buffers) {
      if (progress_cb) progress_cb("load buffer", progress.x++, progress.y);
      auto& path = buffers_paths[&buffer - &buffers.front()];
      if (!load_binary(path_join(dirname, path), buffer, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load buffers
    parallel_foreach(buffers, [&](auto& buffer) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load buffer", progress.x++, progress.y);
      }
      auto& path = buffers_paths[&buffer - &buffers.front()];
      auto  err  = string{};
      if (!load_binary(path_join(dirname, path), buffer, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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

  // convert color textures
  auto get_texture = [&gltf](const njson& js,
                         const string&    name) -> texture_handle {
    if (!js.contains(name)) return invalid_handle;
    auto& ginfo    = js.at(name);
    auto& gtexture = gltf.at("textures").at(ginfo.value("index", -1));
    return gtexture.value("source", -1);
  };

  // convert textures
  auto textures_paths = vector<string>{};
  if (gltf.contains("images")) {
    try {
      for (auto& gimage : gltf.at("images")) {
        scene.textures.emplace_back();
        textures_paths.push_back(gimage.value("uri", ""));
      }
    } catch (...) {
      return parse_error();
    }
  }

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

  if (noparallel) {
    // load texture
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      auto& path = textures_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      }
      auto& path = textures_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
  auto progress = vec2i{
      0, 2 + (int)scene.shapes.size() + (int)scene.textures.size()};
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // convert scene to json
  auto gltf = njson::object();

  // asset
  {
    auto& gasset        = gltf["asset"];
    gasset              = njson::object();
    gasset["version"]   = "2.0";
    gasset["generator"] = scene.asset.generator;
    gasset["copyright"] = scene.asset.copyright;
  }

  // cameras
  if (!scene.cameras.empty()) {
    auto& gcameras = gltf["cameras"];
    gcameras       = njson::array();
    for (auto& camera : scene.cameras) {
      auto& gcamera               = gcameras.emplace_back();
      gcamera                     = njson::object();
      gcamera["name"]             = get_camera_name(scene, camera);
      gcamera["type"]             = "perspective";
      auto& gperspective          = gcamera["perspective"];
      gperspective                = njson::object();
      gperspective["aspectRatio"] = camera.aspect;
      gperspective["yfov"]        = 0.660593;  // TODO(fabio): yfov
      gperspective["znear"]       = 0.001;     // TODO(fabio): configurable?
    }
  }

  // textures
  if (!scene.textures.empty()) {
    auto& gtextures  = gltf["textures"];
    gtextures        = njson::array();
    auto& gsamplers  = gltf["samplers"];
    gsamplers        = njson::array();
    auto& gimages    = gltf["images"];
    gimages          = njson::array();
    auto& gsampler   = gsamplers.emplace_back();
    gsampler         = njson::object();
    gsampler["name"] = "sampler";
    for (auto& texture : scene.textures) {
      auto  name          = get_texture_name(scene, texture);
      auto& gimage        = gimages.emplace_back();
      gimage              = njson::object();
      gimage["name"]      = name;
      gimage["uri"]       = "textures/" + name + ".png";
      auto& gtexture      = gtextures.emplace_back();
      gtexture            = njson::object();
      gtexture["name"]    = name;
      gtexture["sampler"] = 0;
      gtexture["source"]  = (int)gimages.size() - 1;
    }
  }

  // materials
  if (!scene.materials.empty()) {
    auto& gmaterials = gltf["materials"];
    gmaterials       = njson::array();
    for (auto& material : scene.materials) {
      auto& gmaterial             = gmaterials.emplace_back();
      gmaterial                   = njson::object();
      gmaterial["name"]           = get_material_name(scene, material);
      gmaterial["emissiveFactor"] = material.emission;
      auto& gpbr                  = gmaterial["pbrMetallicRoughness"];
      gpbr                        = njson::object();
      gpbr["baseColorFactor"]     = vec4f{material.color.x, material.color.y,
          material.color.z, material.opacity};
      gpbr["metallicFactor"]      = material.metallic;
      gpbr["roughnessFactor"]     = material.roughness;
      if (material.emission_tex != invalid_handle) {
        gmaterial["emissiveTexture"]          = njson::object();
        gmaterial["emissiveTexture"]["index"] = material.emission_tex;
      }
      if (material.normal_tex != invalid_handle) {
        gmaterial["normalTexture"]          = njson::object();
        gmaterial["normalTexture"]["index"] = material.normal_tex;
      }
      if (material.color_tex != invalid_handle) {
        gpbr["baseColorTexture"]          = njson::object();
        gpbr["baseColorTexture"]["index"] = material.color_tex;
      }
      if (material.roughness_tex != invalid_handle) {
        gpbr["metallicRoughnessTexture"]          = njson::object();
        gpbr["metallicRoughnessTexture"]["index"] = material.roughness_tex;
      }
    }
  }

  // add an accessor
  auto set_view = [](njson& gview, njson& gbuffer, const auto& data,
                      size_t buffer_id) {
    gview                 = njson::object();
    gview["buffer"]       = buffer_id;
    gview["byteLength"]   = data.size() * sizeof(data.front());
    gview["byteOffset"]   = gbuffer["byteLength"];
    gbuffer["byteLength"] = gbuffer.value("byteLength", (size_t)0) +
                            data.size() * sizeof(data.front());
  };
  auto set_vaccessor = [](njson& gaccessor, const auto& data, size_t view_id,
                           bool minmax = false) {
    static auto types = unordered_map<size_t, string>{
        {1, "SCALAR"}, {2, "VEC2"}, {3, "VEC3"}, {4, "VEC4"}};
    gaccessor                  = njson::object();
    gaccessor["bufferView"]    = view_id;
    gaccessor["componentType"] = 5126;
    gaccessor["count"]         = data.size();
    gaccessor["type"]          = types.at(sizeof(data.front()) / sizeof(float));
    if constexpr (sizeof(data.front()) == sizeof(vec3f)) {
      if (minmax) {
        auto bbox = invalidb3f;
        for (auto& value : data) bbox = merge(bbox, value);
        gaccessor["min"] = bbox.min;
        gaccessor["max"] = bbox.max;
      }
    }
  };
  auto set_iaccessor = [](njson& gaccessor, const auto& data, size_t view_id,
                           bool minmax = false) {
    gaccessor                  = njson::object();
    gaccessor["bufferView"]    = view_id;
    gaccessor["componentType"] = 5125;
    gaccessor["count"] = data.size() * sizeof(data.front()) / sizeof(int);
    gaccessor["type"]  = "SCALAR";
  };

  // meshes
  auto shape_primitives = vector<njson>();
  shape_primitives.reserve(scene.shapes.size());
  if (!scene.shapes.empty()) {
    auto& gaccessors = gltf["accessors"];
    gaccessors       = njson::array();
    auto& gviews     = gltf["bufferViews"];
    gviews           = njson::array();
    auto& gbuffers   = gltf["buffers"];
    gbuffers         = njson::array();
    for (auto& shape : scene.shapes) {
      auto& gbuffer         = gbuffers.emplace_back();
      gbuffer["uri"]        = "shapes/" + get_shape_name(scene, shape) + ".bin";
      gbuffer["byteLength"] = (size_t)0;
      auto& gprimitive      = shape_primitives.emplace_back();
      gprimitive            = njson::object();
      auto& gattributes     = gprimitive["attributes"];
      gattributes           = njson::object();
      if (!shape.positions.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.positions,
            gbuffers.size() - 1);
        set_vaccessor(gaccessors.emplace_back(), shape.positions,
            gviews.size() - 1, true);
        gattributes["POSITION"] = (int)gaccessors.size() - 1;
      }
      if (!shape.normals.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.normals, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.normals, gviews.size() - 1);
        gattributes["NORMAL"] = (int)gaccessors.size() - 1;
      }
      if (!shape.texcoords.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.texcoords,
            gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.texcoords, gviews.size() - 1);
        gattributes["TEXCOORD_0"] = (int)gaccessors.size() - 1;
      }
      if (!shape.colors.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.colors, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.colors, gviews.size() - 1);
        gattributes["COLOR_0"] = (int)gaccessors.size() - 1;
      }
      if (!shape.radius.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.radius, gbuffers.size() - 1);
        set_vaccessor(
            gaccessors.emplace_back(), shape.radius, gviews.size() - 1);
        gattributes["RADIUS"] = (int)gaccessors.size() - 1;
      }
      if (!shape.points.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.points, gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.points, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 0;
      } else if (!shape.lines.empty()) {
        set_view(
            gviews.emplace_back(), gbuffer, shape.lines, gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.lines, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 1;
      } else if (!shape.triangles.empty()) {
        set_view(gviews.emplace_back(), gbuffer, shape.triangles,
            gbuffers.size() - 1);
        set_iaccessor(
            gaccessors.emplace_back(), shape.triangles, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 4;
      } else if (!shape.quads.empty()) {
        auto triangles = quads_to_triangles(shape.quads);
        set_view(
            gviews.emplace_back(), gbuffer, triangles, gbuffers.size() - 1);
        set_iaccessor(gaccessors.emplace_back(), triangles, gviews.size() - 1);
        gprimitive["indices"] = (int)gaccessors.size() - 1;
        gprimitive["mode"]    = 4;
      }
    }
  }

  // meshes
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
  if (!scene.instances.empty()) {
    auto& gmeshes = gltf["meshes"];
    gmeshes       = njson::array();
    for (auto& instance : scene.instances) {
      auto key = mesh_key{instance.shape, instance.material};
      if (mesh_map.find(key) != mesh_map.end()) continue;
      auto& gmesh   = gmeshes.emplace_back();
      gmesh         = njson::object();
      gmesh["name"] = get_shape_name(scene, instance.shape) + "_" +
                      get_material_name(scene, instance.material);
      gmesh["primitives"] = njson::array();
      gmesh["primitives"].push_back(shape_primitives.at(instance.shape));
      gmesh["primitives"].back()["material"] = instance.material;
      mesh_map[key]                          = gmeshes.size() - 1;
    }
  } else if (!scene.shapes.empty()) {
    auto& gmeshes = gltf["meshes"];
    gmeshes       = njson::array();
    auto shape_id = 0;
    for (auto& primitives : shape_primitives) {
      auto& gmesh         = gmeshes.emplace_back();
      gmesh               = njson::object();
      gmesh["name"]       = get_shape_name(scene, shape_id++);
      gmesh["primitives"] = njson::array();
      gmesh["primitives"].push_back(primitives);
    }
  }

  // nodes
  if (!scene.cameras.empty() || !scene.instances.empty()) {
    auto& gnodes   = gltf["nodes"];
    gnodes         = njson::array();
    auto camera_id = 0;
    for (auto& camera : scene.cameras) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = njson::object();
      gnode["name"]   = get_camera_name(scene, camera);
      gnode["matrix"] = frame_to_mat(camera.frame);
      gnode["camera"] = camera_id++;
    }
    for (auto& instance : scene.instances) {
      auto& gnode     = gnodes.emplace_back();
      gnode           = njson::object();
      gnode["name"]   = get_instance_name(scene, instance);
      gnode["matrix"] = frame_to_mat(instance.frame);
      gnode["mesh"]   = mesh_map.at({instance.shape, instance.material});
    }
    // root children
    auto& groot     = gnodes.emplace_back();
    groot           = njson::object();
    groot["name"]   = "root";
    auto& gchildren = groot["children"];
    gchildren       = njson::array();
    for (auto idx = (size_t)0; idx < gnodes.size() - 1; idx++)
      gchildren.push_back(idx);
    // scene
    auto& gscenes     = gltf["scenes"];
    gscenes           = njson::array();
    auto& gscene      = gscenes.emplace_back();
    gscene            = njson::object();
    auto& gscenenodes = gscene["nodes"];
    gscenenodes       = njson::array();
    gscenenodes.push_back(gnodes.size() - 1);
    gltf["scene"] = 0;
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save json
  if (!save_json(filename, gltf, error)) return false;

  // get filename from name
  auto make_filename = [filename](const string& name, const string& group,
                           const string& extension) {
    return path_join(path_dirname(filename), group, name + extension);
  };

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save shapes
    for (auto& shape : scene.shapes) {
      if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
      if (!save_binshape(path_join(dirname, path), shape, error))
        return dependent_error();
    }
    // save textures
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr" : ".png");
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".bin";
      auto err  = string{};
      if (!save_binshape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

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
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // convert cameras
  for (auto& pcamera : pbrt.cameras) {
    auto& camera  = scene.cameras.emplace_back();
    camera.frame  = pcamera.frame;
    camera.aspect = pcamera.aspect;
    camera.film   = 0.036;
    camera.lens   = pcamera.lens;
    camera.focus  = pcamera.focus;
  }

  // convert material
  auto textures_paths = vector<string>{};
  for (auto& ptexture : pbrt.textures) {
    scene.textures.emplace_back();
    textures_paths.push_back(ptexture.filename);
  }

  // material type map
  auto material_type_map = unordered_map<pbrt_mtype, material_type>{
      {pbrt_mtype::matte, material_type::matte},
      {pbrt_mtype::plastic, material_type::plastic},
      {pbrt_mtype::metal, material_type::metal},
      {pbrt_mtype::glass, material_type::glass},
      {pbrt_mtype::thinglass, material_type::thinglass},
      {pbrt_mtype::subsurface, material_type::matte},
  };

  // convert material
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
    material.color_tex = pmaterial.color_tex;
  }

  // convert shapes
  auto shapes_paths = vector<string>{};
  for (auto& pshape : pbrt.shapes) {
    auto& shape = scene.shapes.emplace_back();
    shapes_paths.emplace_back(pshape.filename_);
    shape.positions = pshape.positions;
    shape.normals   = pshape.normals;
    shape.texcoords = pshape.texcoords;
    shape.triangles = pshape.triangles;
    for (auto& uv : shape.texcoords) uv.y = 1 - uv.y;
    if (pshape.instances.empty()) {
      auto& instance    = scene.instances.emplace_back();
      instance.frame    = pshape.frame;
      instance.shape    = (int)scene.shapes.size() - 1;
      instance.material = pshape.material;
    } else {
      for (auto frame : pshape.instances) {
        auto& instance    = scene.instances.emplace_back();
        instance.frame    = frame * pshape.frame;
        instance.shape    = (int)scene.shapes.size() - 1;
        instance.material = pshape.material;
      }
    }
  }

  // convert environments
  for (auto& penvironment : pbrt.environments) {
    auto& environment        = scene.environments.emplace_back();
    environment.frame        = penvironment.frame;
    environment.emission     = penvironment.emission;
    environment.emission_tex = penvironment.emission_tex;
  }

  // lights
  for (auto& plight : pbrt.lights) {
    auto& shape = scene.shapes.emplace_back();
    shapes_paths.emplace_back();
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
  progress.y += (int)scene.shapes.size();
  progress.y += (int)scene.textures.size();

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // load shape
    for (auto& shape : scene.shapes) {
      if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
      auto& path = shapes_paths[&shape - &scene.shapes.front()];
      if (path.empty()) continue;
      if (!load_shape(path_join(dirname, path), shape, error))
        return dependent_error();
    }
    // load texture
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      auto& path = textures_paths[&texture - &scene.textures.front()];
      if (!load_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // load shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load shape", progress.x++, progress.y);
      }
      auto& path = shapes_paths[&shape - &scene.shapes.front()];
      if (path.empty()) return;
      auto err = string{};
      if (!load_shape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // load textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("load texture", progress.x++, progress.y);
      }
      auto& path = textures_paths[&texture - &scene.textures.front()];
      auto  err  = string{};
      if (!load_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
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
  auto progress = vec2i{
      0, 2 + (int)scene.shapes.size() + (int)scene.textures.size()};
  if (progress_cb) progress_cb("convert scene", progress.x++, progress.y);

  // save pbrt
  auto pbrt = pbrt_scene{};

  // convert camera
  auto& camera       = scene.cameras.front();
  auto& pcamera      = add_camera(pbrt);
  pcamera.frame      = camera.frame;
  pcamera.lens       = camera.lens;
  pcamera.aspect     = camera.aspect;
  pcamera.resolution = {1280, (int)(1280 / pcamera.aspect)};

  // convert textures
  for (auto& texture : scene.textures) {
    auto& ptexture    = pbrt.textures.emplace_back();
    ptexture.filename = "textures/" + get_texture_name(scene, texture) +
                        (!texture.hdr.empty() ? ".hdr" : ".png");
  }

  // material type map
  auto material_type_map = unordered_map<material_type, pbrt_mtype>{
      {material_type::matte, pbrt_mtype::matte},
      {material_type::plastic, pbrt_mtype::plastic},
      {material_type::metal, pbrt_mtype::metal},
      {material_type::glass, pbrt_mtype::glass},
      {material_type::thinglass, pbrt_mtype::thinglass},
      {material_type::subsurface, pbrt_mtype::matte},
      {material_type::volume, pbrt_mtype::matte},
  };

  // convert materials
  for (auto& material : scene.materials) {
    auto& pmaterial     = add_material(pbrt);
    pmaterial.name      = get_material_name(scene, material);
    pmaterial.type      = material_type_map.at(material.type);
    pmaterial.emission  = material.emission;
    pmaterial.color     = material.color;
    pmaterial.roughness = material.roughness;
    pmaterial.ior       = material.ior;
    pmaterial.opacity   = material.opacity;
    pmaterial.color_tex = material.color_tex;
  }

  // convert instances
  for (auto& instance : scene.instances) {
    auto& pshape     = add_shape(pbrt);
    pshape.filename_ = get_shape_name(scene, instance.shape) + ".ply";
    pshape.frame     = instance.frame;
    pshape.frend     = instance.frame;
    pshape.material  = instance.material;
  }

  // convert environments
  for (auto& environment : scene.environments) {
    auto& penvironment        = add_environment(pbrt);
    penvironment.emission     = environment.emission;
    penvironment.emission_tex = environment.emission_tex;
  }

  // handle progress
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);

  // save pbrt
  if (!save_pbrt(filename, pbrt, error)) return false;

  // dirname
  auto dirname = path_dirname(filename);

  if (noparallel) {
    // save textures
    for (auto& shape : scene.shapes) {
      if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      if (!save_shape(path_join(dirname, path), shape, error))
        return dependent_error();
    }
    // save shapes
    for (auto& texture : scene.textures) {
      if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr" : ".png");
      if (!save_texture(path_join(dirname, path), texture, error))
        return dependent_error();
    }
  } else {
    // mutex
    auto mutex = std::mutex{};
    // save shapes
    parallel_foreach(scene.shapes, [&](auto& shape) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save shape", progress.x++, progress.y);
      }
      auto path = "shapes/" + get_shape_name(scene, shape) + ".ply";
      auto err  = string{};
      if (!save_shape(path_join(dirname, path), shape, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
    // save textures
    parallel_foreach(scene.textures, [&](auto& texture) {
      {
        auto lock = std::lock_guard{mutex};
        if (!error.empty()) return;
        if (progress_cb) progress_cb("save texture", progress.x++, progress.y);
      }
      auto path = "textures/" + get_texture_name(scene, texture) +
                  (!texture.hdr.empty() ? ".hdr"s : ".png"s);
      auto err = string{};
      if (!save_texture(path_join(dirname, path), texture, err)) {
        auto lock = std::lock_guard{mutex};
        error     = err;
        return;
      }
    });
    if (!error.empty()) return dependent_error();
  }

  // done
  if (progress_cb) progress_cb("save scene", progress.x++, progress.y);
  return true;
}

}  // namespace yocto
