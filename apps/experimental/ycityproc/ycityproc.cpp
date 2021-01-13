//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

//
// This is  a quick demonsration on how to generate city geometry from geojson.
// It is not efficient since it creates a shape for each building.
// It shiould be rewritten to make it faster and probably reqrite roofs.
//

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

#include <iostream>
#include <string>
#include <tuple>
#include <unordered_map>

#include "ext/earcut.hpp"
#include "ext/json.hpp"

using namespace yocto;
using json = nlohmann::json;
using namespace std::string_literals;
using std::array;

using double2 = array<double, 2>;

enum struct geojson_element_type {
  // clang-format off
  building, highway, pedestrian, water, waterway, grass, forest, sand, tree,
  other
  // clang-format on
};
enum struct geojson_tree_type { standard, palm, oak, pine, cypress };
enum struct geojson_roof_type { missing, flat, gabled };
enum struct geojson_building_type { standard, historic };

struct geojson_element {
  string                  name        = "";
  geojson_element_type    type        = geojson_element_type::other;
  geojson_roof_type       roof        = geojson_roof_type::missing;
  geojson_tree_type       tree        = geojson_tree_type::standard;
  geojson_building_type   building    = geojson_building_type::standard;
  string                  colour      = "";
  int                     level       = 0;
  float                   height      = 0;
  float                   roof_height = 0;
  float                   thickness   = 0;
  vector<double2>         coords      = {};
  vector<double2>         new_coords  = {};
  vector<vector<double2>> holes       = {};
  vector<vector<double2>> new_holes   = {};
};

struct geojson_shape {
  string        name  = "";
  generic_shape shape = {};
};

struct geojson_texture {
  string       name = "";
  image<vec4b> img  = {};
};

struct geojson_scene {
  string                  name      = "";
  string                  copyright = "";
  vector<geojson_element> elements  = {};
  vector<geojson_shape>   shapes    = {};
  vector<geojson_texture> textures  = {};
};

int get_building_level(geojson_element_type type, const json& properties) {
  auto is_integer = [](const string& lev) -> bool {
    for (int i = 0; i < lev.size(); i++) {
      if (lev[i] == '.') return false;
    }
    return false;
  };

  auto is_number = [](const string& lev) -> bool {
    for (int i = 0; i < lev.size(); i++) {
      if ((lev[i] >= 'a' && lev[i] <= 'z') ||
          (lev[i] >= 'A' && lev[i] <= 'B') || lev[i] == ';' || lev[i] == ',')
        return false;
    }
    return true;
  };

  auto is_tall = [](const json& properties) -> bool {
    if (!properties.contains("building")) return false;
    auto building = properties.at("building").get<string>();
    return building == "apartments" || building == "residential" ||
           building == "tower" || building == "hotel";
  };

  auto level  = 1;
  auto height = -1.0f;

  if (properties.contains("building:levels")) {
    auto lev = properties.at("building:levels").get<string>();
    if (is_number(lev)) {
      if (is_integer(lev)) {
        level = (int)round(std::stoi(lev)) + 1;
      } else {
        level = (int)round(std::stof(lev)) + 1;
      }
    } else {
      level = 1;
    }
  }

  // Check if the building:height is given in the GeoJson file
  if (type == geojson_element_type::building && properties.contains("height")) {
    auto h = properties.at("height").get<string>();
    if (is_number(h)) height = std::stof(h);
  }

  if (type == geojson_element_type::building &&
      properties.contains("building:height")) {
    auto h = properties.at("building:height").get<string>();
    if (is_number(h)) height = std::stof(h);
  }

  if (height > -1.0) level = int(float(height) / 3.2);

  if (type == geojson_element_type::building &&
      properties.contains("building") && is_tall(properties))
    level = 3;

  return level;
}

float get_height(const geojson_element& element, float scale) {
  if (element.type == geojson_element_type::building && element.level > 0) {
    return (element.level + scale / 20) / 20;
  } else if (element.type == geojson_element_type::water) {
    return 0.0001f;
  } else if (element.type == geojson_element_type::waterway) {
    return 0.0001f;
  } else if (element.type == geojson_element_type::highway) {
    return 0.0005f;
  } else if (element.type == geojson_element_type::pedestrian) {
    return 0.0004f;
  } else if (element.type == geojson_element_type::grass) {
    return 0.0001f;
  } else if (element.type == geojson_element_type::forest) {
    return 0.0001f;
  } else {
    return 0.0001f;
  }
}

float get_roof_height(const string& roof_h, float scale) {
  if (roof_h != "null") {
    return std::stof(roof_h) / scale;
  } else {
    return 0.109f;
  }
}

vec3f get_color(geojson_element_type type) {
  if (type == geojson_element_type::building) {
    return vec3f{0.79, 0.74, 0.62};
  } else if (type == geojson_element_type::highway) {
    return vec3f{0.26, 0.26, 0.28};
  } else if (type == geojson_element_type::pedestrian) {
    return vec3f{0.45, 0.4, 0.27};  // color = vec3f{0.82, 0.82, 0.82};
  } else if (type == geojson_element_type::water) {
    return vec3f{0.72, 0.95, 1.0};
  } else if (type == geojson_element_type::waterway) {
    return vec3f{0.72, 0.95, 1.0};
  } else if (type == geojson_element_type::sand) {
    return vec3f{0.69, 0.58, 0.43};
  } else if (type == geojson_element_type::forest) {
    return vec3f{0.004, 0.25, 0.16};
  } else if (type == geojson_element_type::grass) {
    return vec3f{0.337, 0.49, 0.274};
  } else {
    return vec3f{0.725, 0.71, 0.68};  // floor
  }
}

vec3f get_building_color(const string& building_color) {
  if (building_color == "yellow") {
    return vec3f{0.882, 0.741, 0.294};
  } else if (building_color == " light yellow") {
    return vec3f{0.922, 0.925, 0.498};
  } else if (building_color == "brown") {
    return vec3f{0.808, 0.431, 0.271};
  } else if (building_color == "light brown") {
    return vec3f{0.8, 0.749, 0.596};
  } else if (building_color == "light orange") {
    return vec3f{0.933, 0.753, 0.416};
  } else {
    return vec3f{1.0, 1.0, 1.0};  // white
  }
}

bool geojson_to_scene(sceneio_scene& scene, const geojson_scene& geojson,
    const string& dirname, string& ioerror) {
  scene.asset.name = geojson.name;
  scene.camera_names.emplace_back("camera");
  auto& camera    = scene.cameras.emplace_back();
  camera.frame    = frame3f{{-0.028f, 0.0f, 1.0f}, {0.764f, 0.645f, 0.022f},
      {-0.645f, 0.764f, -0.018f}, {-13.032f, 16.750f, -1.409f}};
  camera.lens     = 0.035;
  camera.aperture = 0.0;
  camera.focus    = 3.9;
  camera.film     = 0.024;
  camera.aspect   = 1;
  scene.instance_names.emplace_back("floor");
  scene.material_names.emplace_back("floor");
  scene.shape_names.emplace_back("floor");
  auto& floor      = scene.instances.emplace_back();
  auto& floors     = scene.shapes.emplace_back();
  auto& floorm     = scene.materials.emplace_back();
  floor.shape      = (int)scene.shapes.size() - 1;
  floor.material   = (int)scene.materials.size() - 1;
  auto floor_size  = 60.0f;
  floors.positions = {{-floor_size, 0, floor_size}, {floor_size, 0, floor_size},
      {floor_size, 0, -floor_size}, {-floor_size, 0, -floor_size}};
  floors.triangles = {{0, 1, 2}, {2, 3, 0}};
  floorm.color     = {0.725, 0.71, 0.68};

  add_sky(scene);

  // add a tree shape
  auto add_tree_shape = [](sceneio_scene& scene, const geojson_scene& geojson,
                            const string& name) -> shape_handle {
    auto shape_it = std::find_if(geojson.shapes.begin(), geojson.shapes.end(),
        [name](auto& shape) { return shape.name == name; });
    if (shape_it == geojson.shapes.end()) return invalid_handle;
    auto& gshape = *shape_it;
    scene.shape_names.emplace_back(name);
    auto& shape     = scene.shapes.emplace_back();
    shape.triangles = gshape.shape.triangles;
    shape.quads     = gshape.shape.quads;
    shape.positions = gshape.shape.positions;
    shape.normals   = gshape.shape.normals;
    shape.texcoords = gshape.shape.texcoords;
    shape.colors    = gshape.shape.colors;
    return (int)scene.shapes.size() - 1;
  };

  // add a tree shape
  auto add_tree_material = [](sceneio_scene& scene, const vec3f& color,
                               const string& name) -> shape_handle {
    scene.material_names.emplace_back(name);
    auto& material = scene.materials.emplace_back();
    material.color = color;
    return (int)scene.materials.size() - 1;
  };

  // add building texture
  auto add_building_texture = [](sceneio_scene&        scene,
                                  const geojson_scene& geojson,
                                  const string&        name) -> shape_handle {
    auto texture_it = std::find_if(geojson.textures.begin(),
        geojson.textures.end(),
        [name](auto& texture) { return texture.name == name; });
    if (texture_it == geojson.textures.end()) return invalid_handle;
    auto& gtexture = *texture_it;
    scene.texture_names.emplace_back(name);
    auto& texture = scene.textures.emplace_back();
    texture.ldr   = gtexture.img;
    return (int)scene.textures.size() - 1;
  };

  // tree shapes
  auto shape_standard = add_tree_shape(scene, geojson, "standard");
  auto shape_palm     = add_tree_shape(scene, geojson, "palm");
  auto shape_pine     = add_tree_shape(scene, geojson, "pine");
  auto shape_cypress  = add_tree_shape(scene, geojson, "cypress");
  auto shape_oak      = add_tree_shape(scene, geojson, "oak");

  // tree materials
  auto material_standard = add_tree_material(
      scene, {0.002, 0.187, 0.008}, "standard");
  auto material_palm = add_tree_material(scene, {0.224, 0.5, 0.06}, "palm");
  auto material_pine = add_tree_material(scene, {0.145, 0.182, 0.036}, "pine");
  auto material_cypress = add_tree_material(
      scene, {0.019, 0.175, 0.039}, "cypress");
  auto material_oak = add_tree_material(scene, {0.084, 0.193, 0.005}, "oak");

  // buidling textures
  auto texture_1        = add_building_texture(scene, geojson, "1");
  auto texture_2        = add_building_texture(scene, geojson, "2");
  auto texture_3        = add_building_texture(scene, geojson, "3");
  auto texture_4        = add_building_texture(scene, geojson, "4");
  auto texture_5        = add_building_texture(scene, geojson, "5");
  auto texture_6        = add_building_texture(scene, geojson, "6");
  auto texture_7        = add_building_texture(scene, geojson, "7");
  auto texture_8        = add_building_texture(scene, geojson, "8");
  auto texture_8_11     = add_building_texture(scene, geojson, "8_11");
  auto texture_10_41    = add_building_texture(scene, geojson, "10_41");
  auto texture_40_71    = add_building_texture(scene, geojson, "40_71");
  auto texture_70_101   = add_building_texture(scene, geojson, "70_101");
  auto texture_more_101 = add_building_texture(scene, geojson, "more_101");
  auto texture_colosseo = add_building_texture(scene, geojson, "colosseo");

  // Check if exists the element of interest
  auto exist_element = false;
  for (auto& element : geojson.elements) {
    if (element.type == geojson_element_type::building ||
        element.type == geojson_element_type::water ||
        element.type == geojson_element_type::waterway ||
        element.type == geojson_element_type::highway ||
        element.type == geojson_element_type::pedestrian ||
        element.type == geojson_element_type::forest ||
        element.type == geojson_element_type::grass ||
        element.type == geojson_element_type::tree) {
      exist_element = true;
    }
  }

  if (exist_element) {
    for (auto& element : geojson.elements) {
      auto name      = element.name;
      auto type_s    = element.type;
      auto type_roof = element.roof;

      if (type_s == geojson_element_type::tree &&
          element.tree == geojson_tree_type::standard) {
        for (auto& elem : element.new_coords) {
          scene.instance_names.emplace_back(name);
          auto& tree    = scene.instances.emplace_back();
          auto  coord   = vec3f{(float)elem[0], 0, (float)elem[1]};
          auto  x       = coord.x + 0.09f;
          auto  z       = coord.z + 0.09f;
          tree.shape    = shape_standard;
          tree.material = material_standard;
          tree.frame = frame3f{vec3f{1.0f, 0.0f, 0.0f}, vec3f{0.0f, 1.0f, 0.0f},
              vec3f{0.0f, 0.0f, 1.0f}, vec3f{x, coord.y, z}};
        }
      } else if (type_s == geojson_element_type::tree &&
                 element.tree == geojson_tree_type::palm) {
        for (auto& elem : element.new_coords) {
          scene.instance_names.emplace_back(name);
          auto& tree    = scene.instances.emplace_back();
          auto  coord   = vec3f{(float)elem[0], 0, (float)elem[1]};
          tree.shape    = shape_palm;
          tree.material = material_palm;
          tree.frame = frame3f{vec3f{1.0f, 0.0f, 0.0f}, vec3f{0.0f, 1.0f, 0.0f},
              vec3f{0.0f, 0.0f, 1.0f}, vec3f{coord.x, coord.y, coord.z}};
        }
      } else if (type_s == geojson_element_type::tree &&
                 element.tree == geojson_tree_type::cypress) {
        for (auto& elem : element.new_coords) {
          scene.instance_names.emplace_back(name);
          auto& tree    = scene.instances.emplace_back();
          auto  coord   = vec3f{(float)elem[0], 0, (float)elem[1]};
          tree.shape    = shape_cypress;
          tree.material = material_cypress;
          tree.frame = frame3f{vec3f{1.0f, 0.0f, 0.0f}, vec3f{0.0f, 1.0f, 0.0f},
              vec3f{0.0f, 0.0f, 1.0f}, vec3f{coord.x, coord.y, coord.z}};
        }
      } else if (type_s == geojson_element_type::tree &&
                 element.tree == geojson_tree_type::oak) {
        for (auto& elem : element.new_coords) {
          scene.instance_names.emplace_back(name);
          auto& tree    = scene.instances.emplace_back();
          auto  coord   = vec3f{(float)elem[0], 0, (float)elem[1]};
          tree.shape    = shape_oak;
          tree.material = material_oak;
          tree.frame = frame3f{vec3f{1.0f, 0.0f, 0.0f}, vec3f{0.0f, 1.0f, 0.0f},
              vec3f{0.0f, 0.0f, 1.0f}, vec3f{coord.x, coord.y, coord.z}};
        }
      } else if (type_s == geojson_element_type::tree &&
                 element.tree == geojson_tree_type::pine) {
        scene.instance_names.emplace_back(name);
        auto& tree = scene.instances.emplace_back();
        for (auto& elem : element.new_coords) {
          auto coord    = vec3f{(float)elem[0], 0, (float)elem[1]};
          tree.shape    = shape_pine;
          tree.material = material_pine;
          tree.frame = frame3f{vec3f{1.0f, 0.0f, 0.0f}, vec3f{0.0f, 1.0f, 0.0f},
              vec3f{0.0f, 0.0f, 1.0f}, vec3f{coord.x, coord.y, coord.z}};
        }
      } else {
        auto polygon = vector<vector<double2>>{};
        scene.instance_names.push_back(name);
        scene.shape_names.push_back(name);
        scene.material_names.push_back(name);
        auto& instance     = scene.instances.emplace_back();
        auto& shape        = scene.shapes.emplace_back();
        auto& material     = scene.materials.emplace_back();
        instance.shape     = (int)scene.shapes.size() - 1;
        instance.material  = (int)scene.materials.size() - 1;
        auto triangles     = vector<vec3i>{};
        auto positions     = vector<vec3f>{};
        auto vect_building = vector<double2>{};
        auto height        = element.height;
        auto level         = element.level > 0 ? element.level : 0;
        auto type          = element.type;

        for (auto& elem : element.new_coords) {
          auto coord = vec3f{(float)elem[0], height, (float)elem[1]};
          positions.push_back(coord);
          vect_building.push_back({coord.x, coord.z});
        }
        polygon.push_back(vect_building);

        auto vect_hole = vector<double2>{};
        for (auto& list : element.new_holes) {
          for (auto& h : list) {
            auto coord = vec3f{(float)h[0], height, (float)h[1]};
            positions.push_back(coord);
            vect_hole.push_back({coord.x, coord.z});
          }
          polygon.push_back(vect_hole);
          vect_hole = {};
        }

        auto num_holes   = (int)element.new_holes.size();
        auto color_given = false;
        if (element.colour != "null") color_given = true;
        auto color = get_color(type);
        if (type_roof == geojson_roof_type::flat && num_holes == 0) {
          type_roof = geojson_roof_type::gabled;
        } else if (name == "building_relation_1834818") {  // colosseo
          material.color = vec3f{0.725, 0.463, 0.361};
        } else if (type == geojson_element_type::building && level < 3 &&
                   element.building != geojson_building_type::historic) {
          material.color = vec3f{0.538, 0.426, 0.347};
        } else if (element.building == geojson_building_type::historic &&
                   color_given) {
          string building_color = element.colour;
          vec3f  build_color    = get_building_color(building_color);
          material.color        = build_color;
        } else {
          material.color = color;
        }

        auto _polygon = positions;
        auto indices  = mapbox::earcut<int>(polygon);
        for (auto k = 0; k < (int)indices.size() - 2; k += 3) {
          triangles.push_back({indices[k], indices[k + 1], indices[k + 2]});
        }
        if (triangles.empty()) print_info("empty shape " + element.name);

        // Water characteristics
        if (type == geojson_element_type::water ||
            type == geojson_element_type::waterway) {
          material.specular     = 1.0f;
          material.transmission = 0.99f;
          material.metallic     = 0.8f;
          material.roughness    = 0.1f;
        }

        // Road characteristics
        if (type == geojson_element_type::highway) {
          material.roughness = 0.9f;
          material.specular  = 0.7f;
        }

        // Filling buildings
        if (type == geojson_element_type::building) {
          scene.instance_names.push_back(name);
          scene.shape_names.push_back(name);
          scene.material_names.push_back(name);
          auto& instance2    = scene.instances.emplace_back();
          auto& shape2       = scene.shapes.emplace_back();
          auto& material2    = scene.materials.emplace_back();
          instance2.shape    = (int)scene.shapes.size() - 1;
          instance2.material = (int)scene.materials.size() - 1;

          material2.color = color;
          auto _polygon2  = positions;

          // Quads on the building sides
          auto quads = vector<vec4i>{};
          for (auto i = 0; i < (int)positions.size(); i++) {
            auto prev_index = i - 1;
            if (prev_index == -1) {
              prev_index = (int)positions.size() - 1;
            }
            auto index = (int)_polygon2.size();
            _polygon2.push_back({positions[i].x, 0, positions[i].z});
            auto index_2 = (int)_polygon2.size();
            _polygon2.push_back(
                {positions[prev_index].x, 0, positions[prev_index].z});
            quads.push_back({prev_index, i, index, index_2});
          }

          material2.color = color;

          if (element.building == geojson_building_type::historic) {
            if (name == "building_relation_1834818") {  // colosseo
              material2.color_tex = texture_colosseo;
            } else if (element.colour != "null") {
              auto building_color = element.colour;
              auto build_color    = get_building_color(building_color);
              material2.color     = build_color;
            } else {
              material2.color = color;
            }
          } else {
            if (level == 1) {
              material2.color_tex = texture_1;
            } else if (level == 2) {
              material2.color_tex = texture_2;
            } else if (level == 3) {
              material2.color_tex = texture_3;
            } else if (level == 4) {
              material2.color_tex = texture_4;
            } else if (level == 5) {
              material2.color_tex = texture_5;
            } else if (level == 6) {
              material2.color_tex = texture_6;
            } else if (level == 7) {
              material2.color_tex = texture_7;
            } else if (level == 8) {
              material2.color_tex = texture_8;
            } else if (level > 8 && level < 11) {
              material2.color_tex = texture_8_11;
            } else if (level > 10 && level < 41) {
              material2.color_tex = texture_10_41;
            } else if (level > 40 && level < 71) {
              material2.color_tex = texture_40_71;
            } else if (level > 70 && level < 101) {
              material2.color_tex = texture_70_101;
            } else if (level > 101) {
              material2.color_tex = texture_more_101;
            }
          }

          shape2.positions = _polygon2;
          shape2.quads     = quads;
        }

        shape.positions = _polygon;
        shape.triangles = triangles;

        // Gabled roof
        if (type_roof == geojson_roof_type::gabled) {
          auto polygon_roof = vector<vector<double2>>{};
          scene.instance_names.push_back(name + "_roof1");
          scene.shape_names.push_back(name + "_roof1");
          scene.material_names.push_back(name + "_roof1");
          auto& instancer     = scene.instances.emplace_back();
          auto& shaper        = scene.shapes.emplace_back();
          auto& materialr     = scene.materials.emplace_back();
          instancer.shape     = (int)scene.shapes.size() - 1;
          instancer.material  = (int)scene.materials.size() - 1;
          auto triangles_roof = vector<vec3i>{};
          auto positions_roof = vector<vec3f>{};
          auto vect_roof      = vector<double2>{};
          auto height         = element.height;
          auto roof_height    = element.roof_height;
          auto centroid_x = 0.0f, centroid_y = 0.0f;
          auto num_vert  = (int)element.new_coords.size();
          auto num_holes = (int)element.new_holes.size();

          if (num_holes == 0) {
            for (auto& elem : element.new_coords) {
              auto coord = vec3f{(float)elem[0], height, (float)elem[1]};
              positions_roof.push_back(coord);
              vect_roof.push_back({coord.x, coord.z});
              centroid_x += coord.x;
              centroid_y += coord.z;
            }

            centroid_x = centroid_x / num_vert;
            centroid_y = centroid_y / num_vert;

            polygon_roof.push_back(vect_roof);

            auto roof_color = vec3f{0.351, 0.096, 0.091};  // brown/red
            materialr.color = roof_color;

            auto _polygon_roof = positions_roof;
            auto indices_roof  = mapbox::earcut<int>(polygon_roof);
            for (int k = 0; k < indices_roof.size() - 2; k += 3) {
              triangles_roof.push_back(
                  {indices_roof[k], indices_roof[k + 1], indices_roof[k + 2]});
            }

            // Filling roofs
            scene.instance_names.push_back(name + "_roof2");
            scene.shape_names.push_back(name + "_roof2");
            scene.material_names.push_back(name + "_roof2");
            auto& instancer2     = scene.instances.emplace_back();
            auto& shaper2        = scene.shapes.emplace_back();
            auto& materialr2     = scene.materials.emplace_back();
            instancer2.shape     = (int)scene.shapes.size() - 1;
            instancer2.material  = (int)scene.materials.size() - 1;
            materialr2.color     = roof_color;
            auto _polygon2_roof  = positions_roof;
            auto triangles2_roof = vector<vec3i>{};
            for (auto i = 0; i < (int)positions_roof.size(); i++) {
              auto prev_index = i - 1;
              if (prev_index == -1) prev_index = (int)positions_roof.size() - 1;
              auto total_height = height + roof_height;
              auto index        = (int)_polygon2_roof.size();
              _polygon2_roof.push_back({centroid_x, total_height, centroid_y});
              auto index_2 = (int)_polygon2_roof.size();
              _polygon2_roof.push_back({centroid_x, total_height, centroid_y});
              triangles2_roof.push_back({prev_index, i, index});
              triangles2_roof.push_back({index, index_2, prev_index});
            }

            shaper2.positions = _polygon2_roof;
            shaper2.triangles = triangles2_roof;
            shaper.positions  = _polygon_roof;
            shaper.triangles  = triangles_roof;
          }
        }
      }
    }
  }

  return true;
}

vector<double2> compute_area(
    double x, double next_x, double y, double next_y, double road_thickness) {
  auto line_1 = vector<double2>{
      {next_x + road_thickness, next_y + road_thickness},
      {next_x - road_thickness, next_y - road_thickness},
      {x - road_thickness, y - road_thickness},
      {x + road_thickness, y + road_thickness}};

  auto vec_x = vector<double>{}, vec_y = vector<double>{};
  for (auto& couple : line_1) {
    vec_x.push_back(couple[0]);
    vec_y.push_back(couple[1]);
  }

  auto shifted_vec_x = vec_x, shifted_vec_y = vec_y;
  std::rotate(
      shifted_vec_x.begin(), shifted_vec_x.begin() + 3, shifted_vec_x.end());
  std::rotate(
      shifted_vec_y.begin(), shifted_vec_y.begin() + 3, shifted_vec_y.end());

  auto sum_first = 0.0;
  for (auto i = 0; i < vec_x.size(); i++) {
    auto first_prod = vec_x[i] * shifted_vec_y[i];
    sum_first += first_prod;
  }

  auto sum_second = 0.0;
  for (auto i = 0; i < vec_y.size(); i++) {
    auto second_prod = vec_y[i] * shifted_vec_x[i];
    sum_second += second_prod;
  }

  auto area_1 = fabs((float)sum_first - (float)sum_second) / 2;
  auto line_2 = vector<double2>{{next_x + road_thickness, next_y},
      {next_x - road_thickness, next_y}, {x - road_thickness, y},
      {x + road_thickness, y}};

  auto vec_x_2 = vector<double>{}, vec_y_2 = vector<double>{};
  for (auto& couple : line_2) {
    vec_x_2.push_back(couple[0]);
    vec_y_2.push_back(couple[1]);
  }

  auto shifted_vec_x_2 = vec_x_2, shifted_vec_y_2 = vec_y_2;
  std::rotate(shifted_vec_x_2.begin(), shifted_vec_x_2.begin() + 3,
      shifted_vec_x_2.end());
  std::rotate(shifted_vec_y_2.begin(), shifted_vec_y_2.begin() + 3,
      shifted_vec_y_2.end());

  auto sum_first_2 = 0.0;
  for (auto i = 0; i < vec_x_2.size(); i++) {
    auto first_prod_2 = vec_x_2[i] * shifted_vec_y_2[i];
    sum_first_2 += first_prod_2;
  }

  auto sum_second_2 = 0.0;
  for (auto i = 0; i < vec_y_2.size(); i++) {
    auto second_prod_2 = vec_y_2[i] * shifted_vec_x_2[i];
    sum_second_2 += second_prod_2;
  }

  auto area_2 = fabs((float)sum_first_2 - (float)sum_second_2) / 2;
  /*auto line_3 = vector<double2>{{next_x, next_y + road_thickness},
      {next_x, next_y - road_thickness}, {x, y - road_thickness},
      {x, y + road_thickness}};

  auto vec_x_3 = vector<double>{}, vec_y_3 = vector<double>{};
  for (auto& couple : line_3) {
    vec_x_3.push_back(couple[0]);
    vec_y_3.push_back(couple[1]);
  }

  auto shifted_vec_x_3 = vec_x_3, shifted_vec_y_3 = vec_y_3;
  std::rotate(shifted_vec_x_3.begin(), shifted_vec_x_3.begin() + 3,
      shifted_vec_x_3.end());
  std::rotate(shifted_vec_y_3.begin(), shifted_vec_y_3.begin() + 3,
      shifted_vec_y_3.end());

  auto sum_first_3 = 0.0;
  for (auto i = 0; i < vec_x_3.size(); i++) {
    auto first_prod_3 = vec_x_3[i] * shifted_vec_y_3[i];
    sum_first_3 += first_prod_3;
  }

  auto sum_second_3 = 0.0;
  for (auto i = 0; i < vec_y_3.size(); i++) {
    auto second_prod_3 = vec_y_3[i] * shifted_vec_x_3[i];
    sum_second_3 += second_prod_3;
  }

  float area_3 = (float)(0.5f * fabs(sum_first_3 - sum_second_3));
  if (area_2 > area_1) {
    if (area_3 > area_2) {
      return line_3;
    } else {
      return line_2;
    }
  } else {
    if (area_3 > area_1) {
      return line_3;
    } else {
      return line_1;
    }
  }*/

  if (area_2 > area_1) return line_2;
  return line_1;
}

float get_thickness(geojson_element_type type) {
  if (type == geojson_element_type::pedestrian) {
    return 0.00005f;
  } else if (type == geojson_element_type::water ||
             type == geojson_element_type::waterway) {
    return 1.0f;
  } else {
    return 0.0001f;
  }
}

void assign_polygon_type(
    geojson_element& element, const json& properties, float scale) {
  auto is_grass = [](const string& type) -> bool {
    return type == "park" || type == "pitch" || type == "garden" ||
           type == "playground" || type == "greenfield" || type == "scrub" ||
           type == "heath" || type == "farmyard" || type == "grass" ||
           type == "farmland" || type == "village_green" || type == "meadow" ||
           type == "orchard" || type == "vineyard" ||
           type == "recreation_ground" || type == "grassland";
  };

  auto is_pedestrian = [](const json& properties) {
    if (!properties.contains("highway")) return false;
    auto highway_category = properties.at("highway").get<string>();
    return highway_category == "footway" || highway_category == "pedestrian" ||
           highway_category == "track" || highway_category == "steps" ||
           highway_category == "path" || highway_category == "living_street" ||
           highway_category == "pedestrian_area" ||
           highway_category == "pedestrian_line";
  };

  if (properties.contains("building")) {
    element.type = geojson_element_type::building;
    if (properties.contains("roof:shape")) {
      auto roof_shape = properties.at("roof:shape").get<string>();
      if (roof_shape == "gabled" || roof_shape == "onion" ||
          roof_shape == "pyramid") {
        element.roof = geojson_roof_type::gabled;
      } else if (roof_shape == "flat") {
        element.roof = geojson_roof_type::flat;
      } else {
        element.roof = geojson_roof_type::missing;
      }
    }
    if (properties.contains("roof:height")) {
      auto roof_height    = properties.at("roof:height").get<string>();
      element.roof_height = get_roof_height(roof_height, scale);
    }
    if (properties.contains("historic")) {
      element.building = geojson_building_type::historic;
      if (properties.contains("building:colour")) {
        auto color     = properties.at("building:colour").get<string>();
        element.colour = color;
      }
    }
    if (properties.contains("tourism")) {
      auto tourism = properties.at("tourism").get<string>();
      if (tourism == "attraction") {
        element.building = geojson_building_type::historic;
        if (properties.contains("building:colour")) {
          auto color     = properties.at("building:colour").get<string>();
          element.colour = color;
        }
      }
    }
  } else if (properties.contains("water")) {
    element.type = geojson_element_type::water;
  } else if (properties.contains("waterway")) {
    element.type = geojson_element_type::waterway;
  } else if (properties.contains("landuse")) {
    auto landuse = properties.at("landuse").get<string>();
    if (is_grass(landuse)) {
      element.type = geojson_element_type::grass;
    } else if (landuse == "forest") {
      element.type = geojson_element_type::forest;
    } else {
      element.type = geojson_element_type::other;
    }
  } else if (properties.contains("natural")) {
    auto natural = properties.at("natural").get<string>();
    if (natural == "wood") {
      element.type = geojson_element_type::forest;
    } else if (is_grass(natural)) {
      element.type = geojson_element_type::grass;
    } else if (natural == "water") {
      element.type = geojson_element_type::water;
    } else {
      element.type = geojson_element_type::other;
    }
  } else if (properties.contains("leisure")) {
    auto leisure = properties.at("leisure").get<string>();
    if (is_grass(leisure)) {
      element.type = geojson_element_type::grass;
    } else {
      element.type = geojson_element_type::other;
    }
  } else if (properties.contains("highway")) {
    if (is_pedestrian(properties)) {
      element.type = geojson_element_type::pedestrian;
    } else {
      element.type = geojson_element_type::highway;
    }
  } else {
    element.type = geojson_element_type::other;
  }
}

void assign_line_type(
    geojson_element& line, const json& properties, float scale) {
  auto is_pedestrian = [](const json& properties) {
    if (!properties.contains("highway")) return false;
    auto highway = properties.at("highway").get<string>();
    return highway == "footway" || highway == "pedestrian" ||
           highway == "track" || highway == "steps" || highway == "path" ||
           highway == "living_street" || highway == "pedestrian_area" ||
           highway == "pedestrian_line";
  };

  if (properties.contains("highway")) {
    bool pedestrian = is_pedestrian(properties);
    if (pedestrian) {
      line.type = geojson_element_type::pedestrian;
    } else {
      line.type = geojson_element_type::highway;
    }
  } else if (properties.contains("natural")) {
    auto natural = properties.at("natural").get<string>();
    line.type    = geojson_element_type::grass;
  } else {
    line.type = geojson_element_type::other;
  }
}

void assign_multiline_type(
    geojson_element& line, const json& properties, float scale) {
  if (properties.contains("waterway")) {
    line.type = geojson_element_type::waterway;
  } else {
    line.type = geojson_element_type::other;
  }
}

void assign_tree_type(geojson_element& point, const json& properties) {
  if (properties.contains("natural")) {
    auto point_type_nat = properties.at("natural").get<string>();
    if (point_type_nat == "tree") {
      point.type = geojson_element_type::tree;
      if (properties.contains("type")) {
        auto type_tree = properties.at("type").get<string>();
        if (type_tree == "palm") {
          point.tree = geojson_tree_type::palm;
        } else if (type_tree == "pine") {
          point.tree = geojson_tree_type::pine;
        } else if (type_tree == "cypress") {
          point.tree = geojson_tree_type::cypress;
        } else {
          point.tree = geojson_tree_type::standard;
        }
      } else if (properties.contains("tree")) {
        point.type = geojson_element_type::tree;
        point.tree = geojson_tree_type::standard;
      } else if (properties.contains("genus")) {
        point.type      = geojson_element_type::tree;
        auto genus_tree = properties.at("genus").get<string>();
        if (genus_tree == "Quercus") {
          point.tree = geojson_tree_type::oak;
        } else if (genus_tree == "Cupressus") {
          point.tree = geojson_tree_type::cypress;
        } else if (genus_tree == "Pinus") {
          point.tree = geojson_tree_type::pine;
        } else {
          point.tree = geojson_tree_type::standard;
        }
      } else {
        point.type = geojson_element_type::tree;
        point.tree = geojson_tree_type::standard;
      }
    }
  } else {
    point.type = geojson_element_type::other;
  }
}

bool check_valid_type(const geojson_element& element) {
  return element.type == geojson_element_type::building ||
         element.type == geojson_element_type::water ||
         element.type == geojson_element_type::waterway ||
         element.type == geojson_element_type::sand ||
         element.type == geojson_element_type::grass ||
         element.type == geojson_element_type::highway ||
         element.type == geojson_element_type::pedestrian ||
         element.type == geojson_element_type::forest;
}

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

bool load_geojson(const string& filename, geojson_scene& geojson,
    string& ioerror, float scale = 50) {
  // load json
  auto js = json{};
  if (!load_json(filename, js, ioerror)) return false;

  // parse features
  for (auto& feature : js.at("features")) {
    auto geometry   = feature.at("geometry");
    auto properties = feature.at("properties");
    auto id         = properties.at("@id").get<string>();
    std::replace(id.begin(), id.end(), '/', '_');  // replace all '/' to '_'
    auto type = geometry.at("type").get<string>();

    if (type == "Polygon") {
      auto element             = geojson_element{};
      auto count               = 0;
      int  multi_polygon_count = 0;
      assign_polygon_type(element, properties, scale);
      if (element.type == geojson_element_type::other) continue;
      // element.name  = "building_" + id;
      element.level = get_building_level(element.type, properties);
      // auto first    = true;
      for (auto& list_coords : geometry.at("coordinates")) {
        if (count == 0) {  // outer polygon
          element.name = "building_" + id +
                         std::to_string(multi_polygon_count++);
          element.coords = list_coords.get<vector<double2>>();
          // first          = false;
          count++;
        } else {  // analysis of building holes
          element.holes.push_back(list_coords.get<vector<double2>>());
          count++;
        }
        geojson.elements.push_back(element);
        count = 0;
      }
    } else if (type == "MultiPolygon") {
      auto element             = geojson_element{};
      auto multi_polygon_count = 0;
      auto count               = 0;
      assign_polygon_type(element, properties, scale);
      if (element.type == geojson_element_type::other) continue;
      element.level = get_building_level(element.type, properties);
      // auto first    = true;
      for (auto& multi_pol : geometry.at("coordinates")) {
        auto num_lists = multi_pol.size();
        for (auto& list_coords : multi_pol) {
          if (count == 0) {  // outer polygon
            element.name = "building_" + id +
                           std::to_string(multi_polygon_count++);
            element.coords = list_coords.get<vector<double2>>();
            // first          = false;
            count++;
          } else {  // analysis of building holes
            element.holes.push_back(list_coords.get<vector<double2>>());
            count++;
          }
          if (count == num_lists) {
            geojson.elements.push_back(element);
          }
        }
        count = 0;
      }
    } else if (geometry.at("type") == "LineString") {
      auto count = 0;
      for (auto i = 0; i < (int)geometry.at("coordinates").size() - 1; i++) {
        auto [x0, y0]  = geometry.at("coordinates")[i + 0].get<double2>();
        auto [x1, y1]  = geometry.at("coordinates")[i + 1].get<double2>();
        auto thickness = 0.00005f;
        auto area      = compute_area(x0, x1, y0, y1, thickness);
        auto line      = geojson_element{};
        assign_line_type(line, properties, scale);
        if (line.type == geojson_element_type::other) continue;
        line.name      = "line_" + id + std::to_string(count++);
        line.thickness = get_thickness(line.type);
        line.coords    = area;
        geojson.elements.push_back(line);
      }
    } else if (geometry.at("type") == "MultiLineString") {
      auto count = 0;
      for (auto& list_line : geometry.at("coordinates")) {
        for (auto i = 0; i < (int)list_line.size() - 1; i++) {
          auto [x0, y0]  = list_line[i + 0].get<double2>();
          auto [x1, y1]  = list_line[i + 1].get<double2>();
          auto thickness = 0.0004f;
          auto area      = compute_area(x0, x1, y0, y1, thickness);
          auto line      = geojson_element{};
          assign_multiline_type(line, properties, scale);
          if (line.type == geojson_element_type::other) continue;
          line.name      = "multiline_" + id + std::to_string(count++);
          line.thickness = thickness;
          line.coords    = area;
          geojson.elements.push_back(line);
        }
      }
    } else if (geometry.at("type") == "Point") {
      auto point = geojson_element{};
      assign_tree_type(point, properties);
      if (point.type == geojson_element_type::other) continue;
      point.name   = "point_" + id;
      point.coords = {geometry.at("coordinates").get<double2>()};
      geojson.elements.push_back(point);
    }
  }

  // compute bounds
  auto bounds_min = double2{
      std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
  auto bounds_max = double2{std::numeric_limits<double>::lowest(),
      std::numeric_limits<double>::lowest()};
  for (auto& element : geojson.elements) {
    for (auto coord : element.coords) {
      bounds_min = {
          std::min(coord[0], bounds_min[0]), std::min(coord[1], bounds_min[1])};
      bounds_max = {
          std::max(coord[0], bounds_max[0]), std::max(coord[1], bounds_max[1])};
    }
  }

  // scale elements
  for (auto& element : geojson.elements) {
    element.height     = get_height(element, scale);
    element.new_coords = element.coords;
    for (auto& [x, y] : element.new_coords) {
      x = (x - bounds_min[0]) / (bounds_max[0] - bounds_min[0]) * scale -
          (scale / 2);
      y = (y - bounds_min[1]) / (bounds_max[1] - bounds_min[1]) * scale -
          (scale / 2);
    }
    element.new_holes = element.holes;
    for (auto& hole : element.new_holes) {
      for (auto& [x, y] : hole) {
        x = (x - bounds_min[0]) / (bounds_max[0] - bounds_min[0]) * scale -
            (scale / 2);
        y = (y - bounds_min[1]) / (bounds_max[1] - bounds_min[1]) * scale -
            (scale / 2);
      }
    }
  }

  return true;
}

int main(int argc, const char* argv[]) {
  // command line parameters
  auto validate   = false;
  auto info       = false;
  auto copyright  = ""s;
  auto add_skyenv = false;
  auto output     = "out.json"s;
  auto path       = ""s;

  // parse command line
  auto cli = make_cli("ycityproc", "Process scene");
  add_optional(cli, "info", info, "print scene info", {}, "i");
  add_optional(cli, "copyright", copyright, "copyright string", {}, "c");
  add_optional(cli, "validate", validate, "Validate scene");
  add_optional(cli, "output", output, "output scene", {}, "o");
  add_positional(cli, "dirname", path, "input directory");
  parse_cli(cli, argc, argv);

  // load data
  auto geojson = geojson_scene{};
  auto ioerror = ""s;
  print_progress("load geojsons", 0, 1);
  for (auto& filename : list_directory(path)) {
    if (path_extension(filename) != ".geojson") continue;
    if (!load_geojson(filename, geojson, ioerror)) print_fatal(ioerror);
  }
  print_progress("load geojsons", 1, 1);
  print_progress("load shapes", 0, 1);
  for (auto& filename : list_directory(path_join(path, "tree_models"))) {
    if (path_extension(filename) != ".ply") continue;
    auto& shape = geojson.shapes.emplace_back();
    shape.name  = path_basename(filename);
    if (!load_shape(filename, shape.shape, ioerror)) print_fatal(ioerror);
  }
  print_progress("load shapes", 1, 1);
  print_progress("load textures", 0, 1);
  for (auto& filename : list_directory(path_join(path, "buildings_texture"))) {
    if (path_extension(filename) != ".jpg") continue;
    auto& texture = geojson.textures.emplace_back();
    texture.name  = path_basename(filename);
    if (!load_image(filename, texture.img, ioerror)) print_fatal(ioerror);
  }
  print_progress("load textures", 1, 1);

  // Create city
  auto scene = sceneio_scene{};
  print_progress("convert scene", 0, 1);
  if (!geojson_to_scene(scene, geojson, path, ioerror)) print_fatal(ioerror);
  print_progress("convert scene", 1, 1);

  // sky
  if (add_skyenv) add_sky(scene);

  // print info
  if (info) {
    print_info("scene stats ------------");
    for (auto stat : scene_stats(scene)) print_info(stat);
  }

  // make a directory if needed
  if (!make_directory(path_dirname(output), ioerror)) print_fatal(ioerror);
  if (!scene.shapes.empty()) {
    if (!make_directory(path_join(path_dirname(output), "shapes"), ioerror))
      print_fatal(ioerror);
  }
  if (!scene.textures.empty()) {
    if (!make_directory(path_join(path_dirname(output), "textures"), ioerror))
      print_fatal(ioerror);
  }

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // Done
  return 0;
}
