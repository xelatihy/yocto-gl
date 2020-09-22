//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#include <yocto/yocto_commonio.h>
#include <yocto/yocto_sceneio.h>
#include <yocto/yocto_shape.h>
#include <yocto/yocto_trace.h>

#include <atomic>
#include <deque>
#include <fstream>
#include <future>
#include <iostream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <unordered_map>

#include "ext/earcut.hpp"
#include "ext/json.hpp"

using namespace yocto;
using json = nlohmann::json;
using namespace std::string_literals;
using std::array;

int scale = 50;  // 10

struct city_object {
  string                           name        = "";
  string                           type        = "";
  string                           roof_shape  = "";
  string                           colour      = "";
  int                              level       = 0;
  float                            height      = 0;
  float                            roof_height = 0;
  string                           historic    = "";
  float                            thickness   = 0;
  vector<array<double, 2>>         coords      = {};
  vector<array<double, 2>>         new_coords  = {};
  vector<vector<array<double, 2>>> holes       = {};
  vector<vector<array<double, 2>>> new_holes   = {};
};

class Coordinate {
 public:
  double x_minimum = __DBL_MAX__;
  double y_minimum = __DBL_MAX__;
  double x_maximum = __DBL_MIN__;
  double y_maximum = __DBL_MIN__;

  void set_x_min(double x_min) {
    if (x_minimum > x_min) x_minimum = x_min;
  }

  void set_y_min(double y_min) {
    if (y_minimum > y_min) y_minimum = y_min;
  }

  void set_y_max(double y_max) {
    if (y_maximum < y_max) y_maximum = y_max;
  }

  void set_x_max(double x_max) {
    if (x_maximum < x_max) x_maximum = x_max;
  }

  void update(double x, double y) {
    set_x_max(x);
    set_x_min(x);
    set_y_max(y);
    set_y_min(y);
  }
};

//  --------------- FUNCTIONS --------------

bool check_high(json properties) {
  json building_category = properties["building"];
  bool high_building     = false;

  if ((building_category == "apartments") ||
      (building_category == "residential") || (building_category == "tower") ||
      (building_category == "hotel")) {
    high_building = true;
  }
  return high_building;
}

bool check_digit(const string& lev) {
  bool digit = true;
  for (int i = 0; i < lev.size(); i++) {
    // std::cout << typeid(lev[i]).name() << std::endl;
    if ((lev[i] >= 'a' && lev[i] <= 'z') || (lev[i] >= 'A' && lev[i] <= 'B') ||
        lev[i] == ';' || lev[i] == ',')
      digit = false;
  }
  return digit;
}

bool check_int(const string& lev) {
  bool integer = true;
  for (int i = 0; i < lev.size(); i++) {
    if (lev[i] == '.') {
      integer = false;
    }
  }
  return integer;
}

int generate_building_level(const string& footprint_type, json properties) {
  int               level         = 1;
  float             height        = -1.0f;
  bool              high_building = false;
  string::size_type sz;

  if (!properties["building:levels"].empty()) {
    string lev = properties["building:levels"];
    /*std::cout << "level" << std::endl;
    std::cout << lev << std::endl;*/
    bool  digit      = check_digit(lev);
    int   n_levels_i = -1;
    float n_levels_f = -1.0f;

    if (digit) {
      bool integer = check_int(lev);
      if (integer) {
        n_levels_i = std::stoi(lev, &sz);
        level      = (int)round(n_levels_i) + 1;
      } else {
        n_levels_f = std::stof(lev, &sz);
        level      = (int)round(n_levels_f) + 1;
      }
    } else {
      level = 1;  // floor level
    }

    // std::cout << "level" << std::endl;
    // std::cout << level << std::endl;
  }

  // Check if the building:height is given in the GeoJson file
  if (footprint_type == "building" && properties.contains("height")) {
    string h     = properties["height"];
    bool   digit = check_digit(h);
    if (digit) {
      height = std::stof(h, &sz);
    }
  }

  if (footprint_type == "building" && !properties["building:height"].empty()) {
    string h     = properties["building:height"];
    bool   digit = check_digit(h);
    if (digit) {
      height = std::stof(h, &sz);
    }
  }

  if (height > -1.0) {
    level = int(float(height) / 3.2);
  }

  high_building = check_high(properties);
  if (footprint_type == "building" && properties.contains("building") &&
      high_building) {
    level = 3;
  }
  // std::cout << level << std::endl;
  return level;
}

float generate_height(city_object building, int scale) {
  float             height = 0.0001f;
  string::size_type sz;

  if (building.type == "building" && building.level > 0) {
    height = (float)(building.level + (scale / 20.0)) /
             20.0;  //(float) (level + (scale / 20.0)) / 20.0;

  } else if (building.type == "water") {
    height = (float)0.0001f;
  } else if (building.type == "highway") {
    height = (float)0.0005f;
  } else if (building.type == "pedestrian") {
    height = (float)0.0004f;
  }
  // std::cout << height << std::endl;

  return height;
}

float generate_roof_height(const string& roof_h, int scale) {
  float             roof_height = 0.109f;
  string::size_type sz;

  if (roof_h != "null") {
    float roof_hei = (float)std::stof(roof_h, &sz);
    roof_height    = (float)roof_hei / scale;
  }

  return roof_height;
}

bool check_grass_type(const string& building_type) {
  bool grass_area = false;
  if (building_type == "park" || building_type == "pitch" ||
      building_type == "garden" || building_type == "playground" ||
      building_type == "greenfield" || building_type == "scrub" ||
      building_type == "heath" || building_type == "farmyard" ||
      building_type == "grass" || building_type == "farmland" ||
      building_type == "village_green" || building_type == "meadow" ||
      building_type == "orchard" || building_type == "vineyard" ||
      building_type == "recreation_ground" || building_type == "grassland") {
    grass_area = true;
  }
  return grass_area;
}

bool check_pedestrian(json properties) {
  json highway_category = properties["highway"];
  bool is_pedestrian    = false;

  if ((highway_category == "footway") || (highway_category == "pedestrian") ||
      (highway_category == "track") || (highway_category == "steps") ||
      (highway_category == "path") || (highway_category == "living_street") ||
      (highway_category == "pedestrian_area") ||
      (highway_category == "pedestrian_line")) {
    is_pedestrian = true;
  }
  return is_pedestrian;
}

vec3f get_color(const string& type, bool grass_type) {
  vec3f color = {0.725, 0.71, 0.68};  // floor color
  if (type == "building") {
    color = vec3f{0.79, 0.74, 0.62};
  } else if (type == "highway") {
    color = vec3f{0.26, 0.26, 0.28};
  } else if (type == "pedestrian") {
    color = vec3f{0.45, 0.4, 0.27};  // color = vec3f{0.82, 0.82, 0.82};
  } else if (type == "water") {
    color = vec3f{0.72, 0.95, 1.0};
  } else if (type == "sand") {
    color = vec3f{0.69, 0.58, 0.43};
  } else if (type == "forest") {
    color = vec3f{0.004, 0.25, 0.16};
  } else if (grass_type)
    color = vec3f{0.337, 0.49, 0.274};
  return color;
}

vec3f get_building_color(const string& building_color) {
  vec3f color;
  if (building_color == "yellow") {
    color = vec3f{0.882, 0.741, 0.294};
  } else if (building_color == " light yellow") {
    color = vec3f{0.922, 0.925, 0.498};
  } else if (building_color == "brown") {
    color = vec3f{0.808, 0.431, 0.271};
  } else if (building_color == "light brown") {
    color = vec3f{0.8, 0.749, 0.596};
  } else if (building_color == "light orange") {
    color = vec3f{0.933, 0.753, 0.416};
  } else {                         // white
    color = vec3f{1.0, 1.0, 1.0};  // white
  }
  return color;
}

bool create_city_from_json(sceneio_scene* scene,
    vector<city_object> all_geometries, const string& dirname,
    string& ioerror) {
  scene->name      = "cornellbox";
  auto camera      = add_camera(scene);
  camera->frame    = frame3f{{-0.028f, 0.0f, 1.0f}, {0.764f, 0.645f, 0.022f},
      {-0.645f, 0.764f, -0.018f}, {-13.032f, 16.750f, -1.409f}};
  camera->lens     = 0.035;
  camera->aperture = 0.0;
  camera->focus    = 3.9;
  camera->film     = 0.024;
  camera->aspect   = 1;
  auto floor       = add_complete_instance(scene, "floor");
  auto floor_size  = 60.0f;
  floor->shape->positions = {{-floor_size, 0, floor_size},
      {floor_size, 0, floor_size}, {floor_size, 0, -floor_size},
      {-floor_size, 0, -floor_size}};
  floor->shape->triangles = {{0, 1, 2}, {2, 3, 0}};
  floor->material->color  = {0.725, 0.71, 0.68};

  add_sky(scene);

  // standard tree
  string path_standard  = path_join(dirname, "tree_models/standard.ply");
  auto   shape_standard = add_shape(scene, "standard");
  if (!load_shape(path_standard, shape_standard->points, shape_standard->lines,
          shape_standard->triangles, shape_standard->quads,
          shape_standard->quadspos, shape_standard->quadsnorm,
          shape_standard->quadstexcoord, shape_standard->positions,
          shape_standard->normals, shape_standard->texcoords,
          shape_standard->colors, shape_standard->radius, ioerror))
    return false;

  // palm tree
  string path_palm  = path_join(dirname, "tree_models/palm.ply");
  auto   shape_palm = add_shape(scene, "palm");
  if (!load_shape(path_palm, shape_palm->points, shape_palm->lines,
          shape_palm->triangles, shape_palm->quads, shape_palm->quadspos,
          shape_palm->quadsnorm, shape_palm->quadstexcoord,
          shape_palm->positions, shape_palm->normals, shape_palm->texcoords,
          shape_palm->colors, shape_palm->radius, ioerror))
    return false;

  // pine tree
  string path_pine  = path_join(dirname, "tree_models/pine.ply");
  auto   shape_pine = add_shape(scene, "pine");
  if (!load_shape(path_pine, shape_pine->points, shape_pine->lines,
          shape_pine->triangles, shape_pine->quads, shape_pine->quadspos,
          shape_pine->quadsnorm, shape_pine->quadstexcoord,
          shape_pine->positions, shape_pine->normals, shape_pine->texcoords,
          shape_pine->colors, shape_pine->radius, ioerror))
    return false;

  // cypress tree
  string path_cypress  = path_join(dirname, "tree_models/cypress.ply");
  auto   shape_cypress = add_shape(scene, "cypress");
  if (!load_shape(path_cypress, shape_cypress->points, shape_cypress->lines,
          shape_cypress->triangles, shape_cypress->quads,
          shape_cypress->quadspos, shape_cypress->quadsnorm,
          shape_cypress->quadstexcoord, shape_cypress->positions,
          shape_cypress->normals, shape_cypress->texcoords,
          shape_cypress->colors, shape_cypress->radius, ioerror))
    return false;

  // oak tree
  string path_oak  = path_join(dirname, "tree_models/oak.ply");
  auto   shape_oak = add_shape(scene, "oak");
  if (!load_shape(path_oak, shape_oak->points, shape_oak->lines,
          shape_oak->triangles, shape_oak->quads, shape_oak->quadspos,
          shape_oak->quadsnorm, shape_oak->quadstexcoord, shape_oak->positions,
          shape_oak->normals, shape_oak->texcoords, shape_oak->colors,
          shape_oak->radius, ioerror))
    return false;

  // Load textures
  // string texture_path = "./city/buildings_texture/";

  // buidling texture1
  auto   texture_1       = add_texture(scene, "texture1");
  string path_text_1     = path_join(dirname, "buildings_texture/1.jpg");
  auto   build_texture_1 = load_image(path_text_1, texture_1->hdr, ioerror);

  // buidling texture2
  auto   texture_2       = add_texture(scene, "texture2");
  string path_text_2     = path_join(dirname, "buildings_texture/2.jpg");
  auto   build_texture_2 = load_image(path_text_2, texture_2->hdr, ioerror);

  // buidling texture3
  auto   texture_3       = add_texture(scene, "texture3");
  string path_text_3     = path_join(dirname, "buildings_texture/3.jpg");
  auto   build_texture_3 = load_image(path_text_3, texture_3->hdr, ioerror);

  // buidling texture4
  auto   texture_4       = add_texture(scene, "texture4");
  string path_text_4     = path_join(dirname, "buildings_texture/4.jpg");
  auto   build_texture_4 = load_image(path_text_4, texture_4->hdr, ioerror);

  // buidling texture5
  auto   texture_5       = add_texture(scene, "texture5");
  string path_text_5     = path_join(dirname, "buildings_texture/5.jpg");
  auto   build_texture_5 = load_image(path_text_5, texture_5->hdr, ioerror);

  // buidling texture6
  auto   texture_6       = add_texture(scene, "texture6");
  string path_text_6     = path_join(dirname, "buildings_texture/6.jpg");
  auto   build_texture_6 = load_image(path_text_6, texture_6->hdr, ioerror);

  // buidling texture7
  auto   texture_7       = add_texture(scene, "texture7");
  string path_text_7     = path_join(dirname, "buildings_texture/7.jpg");
  auto   build_texture_7 = load_image(path_text_7, texture_7->hdr, ioerror);

  // buidling texture8
  auto   texture_8       = add_texture(scene, "texture8");
  string path_text_8     = path_join(dirname, "buildings_texture/8.jpg");
  auto   build_texture_8 = load_image(path_text_8, texture_8->hdr, ioerror);

  // buidling texture8_11
  auto   texture_8_11       = add_texture(scene, "texture8_11");
  string path_text_8_11     = path_join(dirname, "buildings_texture/8_11.jpg");
  auto   build_texture_8_11 = load_image(
      path_text_8_11, texture_8_11->hdr, ioerror);

  // buidling texture10_41
  auto   texture_10_41   = add_texture(scene, "texture10_41");
  string path_text_10_41 = path_join(dirname, "buildings_texture/10_41.jpg");
  auto   build_texture_10_41 = load_image(
      path_text_10_41, texture_10_41->hdr, ioerror);

  // buidling texture40_71
  auto   texture_40_71   = add_texture(scene, "texture40_71");
  string path_text_40_71 = path_join(dirname, "buildings_texture/40_71.jpg");
  auto   build_texture_40_71 = load_image(
      path_text_40_71, texture_40_71->hdr, ioerror);

  // buidling texture70_101
  auto   texture_70_101   = add_texture(scene, "texture70_101");
  string path_text_70_101 = path_join(dirname, "buildings_texture/70_101.jpg");
  auto   build_texture_70_101 = load_image(
      path_text_70_101, texture_70_101->hdr, ioerror);

  // buidling texturemore_101
  auto   texture_more_101   = add_texture(scene, "texturemore_101");
  string path_text_more_101 = path_join(
      dirname, "buildings_texture/more_101.jpg");
  auto build_texture_more_101 = load_image(
      path_text_more_101, texture_more_101->hdr, ioerror);

  // buidling texture_colosseo
  auto   texture_colosseo   = add_texture(scene, "texture_colosseo");
  string path_text_colosseo = path_join(
      dirname, "buildings_texture/colosseo.jpg");

  // Check if exists the element of interest
  bool exist_element = false;
  for (auto& building_geometry : all_geometries) {
    auto building_type = building_geometry.type;
    bool grass_area    = check_grass_type(building_type);
    if (building_geometry.type == "building" ||
        building_geometry.type == "water" ||
        building_geometry.type == "highway" ||
        building_geometry.type == "pedestrian" ||
        building_geometry.type == "forest" || grass_area ||
        building_geometry.type == "standard" ||
        building_geometry.type == "palm" || building_geometry.type == "pine" ||
        building_geometry.type == "oak" ||
        building_geometry.type == "cypress") {
      exist_element = true;
    }
  }

  if (exist_element) {
    using Coord = double;
    using N     = int32_t;
    using Point = array<Coord, 2>;
    for (auto& element : all_geometries) {
      auto name = element.name;

      string type_s    = element.type;
      string type_roof = "null";
      string historic  = "no";

      if (element.roof_shape != "null") type_roof = element.roof_shape;

      if (element.historic != "no") historic = element.historic;

      if (type_s == "standard") {
        auto  tree  = add_complete_instance(scene, name);
        vec3f coord = {};

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = 0.0;
          coord.z = elem[1];

          auto x = coord.x + 0.09f;
          auto z = coord.z + 0.09f;

          // create TREE object
          tree->shape           = shape_standard;
          tree->material->color = {0.002, 0.187, 0.008};
          tree->frame           = frame3f{vec3f{1.0f, 0.0f, 0.0f},
              vec3f{0.0f, 1.0f, 0.0f}, vec3f{0.0f, 0.0f, 1.0f},
              vec3f{x, coord.y, z}};
          coord                 = {};
          continue;
        }

      } else if (type_s == "palm") {
        auto  tree  = add_complete_instance(scene, name);
        vec3f coord = {};

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = 0.0;
          coord.z = elem[1];

          // create TREE object
          tree->shape           = shape_palm;
          tree->material->color = {0.224, 0.5, 0.06};
          tree->frame           = frame3f{vec3f{1.0f, 0.0f, 0.0f},
              vec3f{0.0f, 1.0f, 0.0f}, vec3f{0.0f, 0.0f, 1.0f},
              vec3f{coord.x, coord.y, coord.z}};

          coord = {};
          continue;
        }
      } else if (type_s == "cypress") {
        auto  tree  = add_complete_instance(scene, name);
        vec3f coord = {};

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = 0.0;
          coord.z = elem[1];

          // create TREE object
          tree->shape           = shape_cypress;
          tree->material->color = {0.019, 0.175, 0.039};
          tree->frame           = frame3f{vec3f{1.0f, 0.0f, 0.0f},
              vec3f{0.0f, 1.0f, 0.0f}, vec3f{0.0f, 0.0f, 1.0f},
              vec3f{coord.x, coord.y, coord.z}};

          coord = {};
          continue;
        }
      } else if (type_s == "oak") {
        auto  tree  = add_complete_instance(scene, name);
        vec3f coord = {};

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = 0.0;
          coord.z = elem[1];

          // create TREE object OAK
          tree->shape           = shape_oak;
          tree->material->color = {0.084, 0.193, 0.005};
          tree->frame           = frame3f{vec3f{1.0f, 0.0f, 0.0f},
              vec3f{0.0f, 1.0f, 0.0f}, vec3f{0.0f, 0.0f, 1.0f},
              vec3f{coord.x, coord.y, coord.z}};

          coord = {};
          continue;
        }
      } else if (type_s == "pine") {
        auto  tree  = add_complete_instance(scene, name);
        vec3f coord = {};

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = 0.0;
          coord.z = elem[1];

          // create TREE object
          tree->shape           = shape_pine;
          tree->material->color = {0.145, 0.182, 0.036};
          tree->frame           = frame3f{vec3f{1.0f, 0.0f, 0.0f},
              vec3f{0.0f, 1.0f, 0.0f}, vec3f{0.0f, 0.0f, 1.0f},
              vec3f{coord.x, coord.y, coord.z}};

          coord = {};
          continue;
        }
      } else {
        vector<vector<Point>> polygon;

        auto          build = add_complete_instance(scene, name);
        vector<vec3i> triangles;
        vector<vec3f> positions;

        vector<Point> vect_building;
        vec3f         coord       = {};
        float         height      = -1.0f;
        float         roof_height = -1.0f;
        int           level       = 0;
        string        type        = "";

        string::size_type sz;

        type = element.type;

        if (element.level > 0) {
          level = element.level;
        }

        height = element.height;

        for (auto& elem : element.new_coords) {
          coord.x = elem[0];
          coord.y = height;
          coord.z = elem[1];

          positions.push_back(coord);
          vect_building.push_back({coord.x, coord.z});

          coord = {};
          continue;
        }
        polygon.push_back(vect_building);

        vector<Point> vect_hole;
        coord = {};

        for (auto& list : element.new_holes) {
          for (auto& h : list) {
            coord.x = h[0];
            coord.z = h[1];
            coord.y = height;

            positions.push_back(coord);
            vect_hole.push_back({coord.x, coord.z});
            coord = {};
          }

          polygon.push_back(vect_hole);
          vect_hole = {};
        }

        int num_holes = element.new_holes.size();

        bool color_given = false;
        if (element.colour != "null") color_given = true;

        bool grass_area = check_grass_type(element.type);
        auto color = get_color(type, grass_area);  // vec3f{0.79, 0.74, 0.62};

        if (type_roof == "flat" && num_holes == 0) {
          type_roof = "gabled";
        } else if (name == "building_relation_1834818") {  // colosseo
          build->material->color = vec3f{0.725, 0.463, 0.361};
        } else if (type == "building" && level < 3 && historic != "no") {
          build->material->color = vec3f{
              0.538, 0.426, 0.347};  // vec3f{0.402,0.319,0.261}; // light brown
        } else if (historic == "yes" && color_given) {
          string building_color  = element.colour;
          vec3f  build_color     = get_building_color(building_color);
          build->material->color = build_color;
        } else {
          build->material->color = color;
        }

        vector<vec3f> _polygon;
        for (auto point : positions) {
          _polygon.push_back(point);
        }

        vector<N> indices = mapbox::earcut<N>(polygon);

        for (int k = 0; k < indices.size() - 2; k += 3) {
          triangles.push_back({indices[k], indices[k + 1], indices[k + 2]});
        }

        // Water characteristics
        if (type == "water") {
          build->material->specular     = 1.0f;
          build->material->transmission = 0.99f;
          build->material->metallic     = 0.8f;
          build->material->roughness    = 0.1f;
        }

        // Road characteristics
        if (type == "highway") {
          build->material->roughness = 0.9f;
          build->material->specular  = 0.7f;
        }

        // Filling buildings
        if (type == "building") {
          auto build2             = add_complete_instance(scene, name + "_1");
          build2->material->color = color;
          vector<vec3f> _polygon2;
          for (auto point : positions) {
            _polygon2.push_back(point);
          }

          // Quads on the building sides
          vector<vec4i> quads;
          for (int i = 0; i < positions.size(); i++) {
            auto prev_index = i - 1;
            if (prev_index == -1) {
              prev_index = positions.size() - 1;
            }
            auto index = (int)_polygon2.size();
            _polygon2.push_back({positions[i].x, 0, positions[i].z});
            auto index_2 = (int)_polygon2.size();
            _polygon2.push_back(
                {positions[prev_index].x, 0, positions[prev_index].z});

            quads.push_back({prev_index, i, index, index_2});
          }

          build2->material->color = color;

          if (historic == "yes") {
            if (name == "building_relation_1834818") {  // colosseo
              auto build_texture_colosseo = load_image(
                  path_text_colosseo, texture_colosseo->hdr, ioerror);
              build2->material->color_tex = texture_colosseo;
            } else if (element.colour != "null") {
              string building_color   = element.colour;
              vec3f  build_color      = get_building_color(building_color);
              build2->material->color = build_color;
            } else {
              build2->material->color = color;
            }
          } else {
            if (level == 1)
              build2->material->color_tex = texture_1;
            else if (level == 2)
              build2->material->color_tex = texture_2;
            else if (level == 3)
              build2->material->color_tex = texture_3;
            else if (level == 4)
              build2->material->color_tex = texture_4;
            else if (level == 5)
              build2->material->color_tex = texture_5;
            else if (level == 6)
              build2->material->color_tex = texture_6;
            else if (level == 7)
              build2->material->color_tex = texture_7;
            else if (level == 8)
              build2->material->color_tex = texture_8;
            else if (level > 8 && level < 11)
              build2->material->color_tex = texture_8_11;
            else if (level > 10 && level < 41)
              build2->material->color_tex = texture_10_41;
            else if (level > 40 && level < 71)
              build2->material->color_tex = texture_40_71;
            else if (level > 70 && level < 101)
              build2->material->color_tex = texture_70_101;
            else if (level > 101)
              build2->material->color_tex = texture_more_101;
          }

          build2->shape->positions = _polygon2;
          build2->shape->quads     = quads;
        }

        build->shape->positions = _polygon;
        build->shape->triangles = triangles;

        // Gabled roof
        if (type_roof == "gabled") {
          vector<vector<Point>> polygon_roof;

          auto          roof = add_complete_instance(scene, name);
          vector<vec3i> triangles_roof;
          vector<vec3f> positions_roof;

          vector<Point> vect_roof;
          float         roof_height = -1.0f;
          vec3f         coord       = {};
          float         height      = -1.0f;

          height      = element.height;
          roof_height = element.roof_height;

          float centroid_x = 0.0f;
          float centroid_y = 0.0f;
          int   num_vert   = (int)element.new_coords.size();
          int   num_holes  = element.new_holes.size();

          if (num_holes == 0) {
            for (auto& elem : element.new_coords) {
              coord.x = (double)elem[0];
              coord.y = height;
              coord.z = (double)elem[1];

              positions_roof.push_back(coord);
              vect_roof.push_back({coord.x, coord.z});

              centroid_x += coord.x;
              centroid_y += coord.z;

              coord = {};
              continue;
            }

            centroid_x = centroid_x / num_vert;
            centroid_y = centroid_y / num_vert;

            polygon_roof.push_back(vect_roof);

            auto roof_color       = vec3f{0.351, 0.096, 0.091};  // brown/red
            roof->material->color = roof_color;

            vector<vec3f> _polygon_roof;
            for (auto point : positions_roof) {
              _polygon_roof.push_back(point);
            }

            vector<N> indices_roof = mapbox::earcut<N>(polygon_roof);

            for (int k = 0; k < indices_roof.size() - 2; k += 3) {
              triangles_roof.push_back(
                  {indices_roof[k], indices_roof[k + 1], indices_roof[k + 2]});
            }

            // Filling roofs
            auto roof2 = add_complete_instance(scene, name + "_roof");
            roof2->material->color = roof_color;
            vector<vec3f> _polygon2_roof;
            for (auto point : positions_roof) {
              _polygon2_roof.push_back(point);
            }
            vector<vec3i> triangles2_roof;
            for (int i = 0; i < positions_roof.size(); i++) {
              auto prev_index = i - 1;
              if (prev_index == -1) {
                prev_index = positions_roof.size() - 1;
              }
              auto total_height = height + roof_height;
              auto index        = (int)_polygon2_roof.size();
              _polygon2_roof.push_back({centroid_x, total_height, centroid_y});
              auto index_2 = (int)_polygon2_roof.size();
              _polygon2_roof.push_back({centroid_x, total_height, centroid_y});
              triangles2_roof.push_back({prev_index, i, index});
              triangles2_roof.push_back({index, index_2, prev_index});
            }

            roof2->shape->positions = _polygon2_roof;
            roof2->shape->triangles = triangles2_roof;

            roof->shape->positions = _polygon_roof;
            roof->shape->triangles = triangles_roof;
          }
        }
      }
    }
  }
  return true;
}

vector<array<double, 2>> compute_area(
    double x, double next_x, double y, double next_y, double road_thickness) {
  vector<array<double, 2>> line_1 = {
      {next_x + road_thickness, next_y + road_thickness},
      {next_x - road_thickness, next_y - road_thickness},
      {x - road_thickness, y - road_thickness},
      {x + road_thickness, y + road_thickness}};

  vector<double> vec_x = {};
  vector<double> vec_y = {};

  for (auto& couple : line_1) {
    vec_x.push_back(couple[0]);
    vec_y.push_back(couple[1]);
  }

  vector<double> shifted_vec_x = vec_x;
  vector<double> shifted_vec_y = vec_y;

  std::rotate(
      shifted_vec_x.begin(), shifted_vec_x.begin() + 3, shifted_vec_x.end());
  std::rotate(
      shifted_vec_y.begin(), shifted_vec_y.begin() + 3, shifted_vec_y.end());

  double sum_first = 0.0f;
  for (int i = 0; i < vec_x.size(); i++) {
    double first_prod = vec_x[i] * shifted_vec_y[i];
    sum_first += first_prod;
  }

  double sum_second = 0.0f;
  for (int i = 0; i < vec_y.size(); i++) {
    double second_prod = vec_y[i] * shifted_vec_x[i];
    sum_second += second_prod;
  }

  float area_1 = (float)(0.5f * fabs(sum_first - sum_second));

  // -----------
  vector<array<double, 2>> line_2 = {{next_x + road_thickness, next_y},
      {next_x - road_thickness, next_y}, {x - road_thickness, y},
      {x + road_thickness, y}};

  vector<double> vec_x_2 = {};
  vector<double> vec_y_2 = {};

  for (auto& couple : line_2) {
    vec_x_2.push_back(couple[0]);
    vec_y_2.push_back(couple[1]);
  }

  vector<double> shifted_vec_x_2 = vec_x_2;
  vector<double> shifted_vec_y_2 = vec_y_2;

  std::rotate(shifted_vec_x_2.begin(), shifted_vec_x_2.begin() + 3,
      shifted_vec_x_2.end());
  std::rotate(shifted_vec_y_2.begin(), shifted_vec_y_2.begin() + 3,
      shifted_vec_y_2.end());

  double sum_first_2 = 0.0f;
  for (int i = 0; i < vec_x_2.size(); i++) {
    double first_prod_2 = vec_x_2[i] * shifted_vec_y_2[i];
    sum_first_2 += first_prod_2;
  }

  double sum_second_2 = 0.0f;
  for (int i = 0; i < vec_y_2.size(); i++) {
    double second_prod_2 = vec_y_2[i] * shifted_vec_x_2[i];
    sum_second_2 += second_prod_2;
  }

  float area_2 = (float)(0.5f * fabs(sum_first_2 - sum_second_2));

  // -----------
  vector<array<double, 2>> line_3 = {{next_x, next_y + road_thickness},
      {next_x, next_y - road_thickness}, {x, y - road_thickness},
      {x, y + road_thickness}};

  vector<double> vec_x_3 = {};
  vector<double> vec_y_3 = {};

  for (auto& couple : line_3) {
    vec_x_3.push_back(couple[0]);
    vec_y_3.push_back(couple[1]);
  }

  vector<double> shifted_vec_x_3 = vec_x_3;
  vector<double> shifted_vec_y_3 = vec_y_3;

  std::rotate(shifted_vec_x_3.begin(), shifted_vec_x_3.begin() + 3,
      shifted_vec_x_3.end());
  std::rotate(shifted_vec_y_3.begin(), shifted_vec_y_3.begin() + 3,
      shifted_vec_y_3.end());

  double sum_first_3 = 0.0f;
  for (int i = 0; i < vec_x_3.size(); i++) {
    double first_prod_3 = vec_x_3[i] * shifted_vec_y_3[i];
    sum_first_3 += first_prod_3;
  }

  double sum_second_3 = 0.0f;
  for (int i = 0; i < vec_y_3.size(); i++) {
    double second_prod_3 = vec_y_3[i] * shifted_vec_x_3[i];
    sum_second_3 += second_prod_3;
  }

  float area_3 = (float)(0.5f * fabs(sum_first_3 - sum_second_3));

  /*if (area_2 > area_1) {
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

float get_thickness(string type) {
  float thickness = 0.0001;
  if (type == "pedestrian") {
    thickness = 0.00005;
  } else if (type == "water") {  // MultiLineString
    thickness = 1.0;
  }
  return thickness;
}

city_object assign_type(city_object building, json properties) {
  if (properties.contains("building")) {
    building.type = "building";
    if (!properties["roof:shape"].empty()) {
      string roof_shape = properties["roof:shape"];
      if (roof_shape == "gabled" || roof_shape == "onion" ||
          roof_shape == "pyramid")
        building.roof_shape = "gabled";
      else if (roof_shape == "flat")
        building.roof_shape = "flat";
    }

    if (!properties["roof:height"].empty()) {
      string roof_h        = properties["roof:height"];
      double roof_height   = generate_roof_height(roof_h, scale);
      building.roof_height = roof_height;
    }

    if (properties.contains("historic")) {
      building.historic = "yes";
      if (!properties["building:colour"].empty()) {
        string build_colour = properties["building:colour"];
        building.colour     = build_colour;
      }
    }

    if (properties.contains("tourism")) {
      string tourism = properties["tourism"];
      if (tourism == "attraction") {
        building.historic = "yes";
        if (!properties["building:colour"].empty()) {
          string build_colour = properties["building:colour"];
          building.colour     = build_colour;
        }
      }
    }
  }

  else if (properties.contains("water")) {
    building.type = "water";
  }

  else if (properties.contains("landuse")) {
    string landuse = properties["landuse"];
    building.type  = landuse;
  }

  else if (properties.contains("natural")) {
    string natural = properties["natural"];
    if (natural == "wood") {
      building.type = "forest";
    } else {
      building.type = natural;
    }
  }

  else if (properties.contains("leisure")) {
    string leisure = properties["leisure"];
    building.type  = leisure;
  }

  else if (properties.contains("highway")) {
    bool pedestrian = check_pedestrian(properties);
    if (pedestrian) {
      building.type = "pedestrian";
    } else {
      building.type = "highway";
    }
  }

  else {
    building.type = "null";
  }

  return building;
}

vector<city_object> assign_tree_type(
    city_object point, json properties, vector<city_object> all_buildings) {
  if (properties.contains("natural")) {
    string point_type_nat = properties["natural"];
    if (point_type_nat == "tree") {
      if (properties.contains("type")) {
        string type_tree = properties["type"];
        if (type_tree == "palm") {
          point.type = "palm";
          all_buildings.push_back(point);
        } else if (type_tree == "pine") {
          point.type = "pine";
          all_buildings.push_back(point);
        } else if (type_tree == "cypress") {
          point.type = "cypress";
          all_buildings.push_back(point);
        } else {
          point.type = "standard";
          all_buildings.push_back(point);
        }
      } else if (properties.contains("tree")) {
        point.type = "standard";
        all_buildings.push_back(point);
      } else if (properties.contains("genus")) {
        string genus_tree = properties["genus"];
        if (genus_tree == "Quercus") {
          point.type = "oak";
          all_buildings.push_back(point);
        } else if (genus_tree == "Cupressus") {
          point.type = "cypress";
          all_buildings.push_back(point);
        } else if (genus_tree == "Pinus") {
          point.type = "pine";
          all_buildings.push_back(point);
        } else {
          point.type = "standard";
          all_buildings.push_back(point);
        }
      } else {
        point.type = "standard";
        all_buildings.push_back(point);
      }
    }
  }

  else {
    point.type = "null";
  }

  return all_buildings;
}

bool check_valid_type(city_object building, json properties) {
  bool valid = false;

  bool grass_area = check_grass_type(building.type);
  if (building.type == "building" || building.type == "water" ||
      building.type == "sand" || grass_area || building.type == "highway" ||
      building.type == "pedestrian" || building.type == "forest")
    valid = true;

  return valid;
}

std::pair<vector<city_object>, Coordinate> data_analysis(json geojson_file,
    vector<city_object> all_buildings, Coordinate class_coord) {
  for (auto& feature : geojson_file["features"]) {
    auto   geometry   = feature["geometry"];
    auto   properties = feature["properties"];
    string id         = properties["@id"];
    std::replace(id.begin(), id.end(), '/', '_');  // replace all '/' to '_'
    int count_list = 0;

    if (geometry["type"] == "Polygon") {
      auto building = city_object();

      building              = assign_type(building, properties);
      string footprint_type = building.type;

      if (footprint_type == "null") {
        continue;
      }

      int level = generate_building_level(footprint_type, properties);

      vector<array<double, 2>>         couple;
      vector<vector<array<double, 2>>> list_holes;
      vector<array<double, 2>>         list_coordinates = {};

      int num_lists = geometry["coordinates"].size();
      for (auto& list_coords : geometry["coordinates"]) {
        if (count_list == 0) {  // outer polygon
          building.level = level;

          for (auto& coord : list_coords) {
            double x = (double)coord[0];
            double y = (double)coord[1];
            class_coord.update(x, y);
            array<double, 2> arr = {x, y};
            list_coordinates.push_back(arr);
          }

          building.coords = list_coordinates;

          string name       = id;
          string build_name = "building_" + name;
          building.name     = build_name;

          count_list++;
        } else {  // analysis of building holes

          for (auto& coord : list_coords) {
            double x = (double)coord[0];
            double y = (double)coord[1];
            couple.push_back({x, y});
          }
          list_holes.push_back(couple);
          couple = {};

          count_list++;
        }
        // std::cout << count << std::endl;

        if (count_list == num_lists) {
          building.holes = list_holes;

          bool valid_type = check_valid_type(building, properties);
          if (valid_type) all_buildings.push_back(building);
        }
      }
      count_list = 0;

    } else if (geometry["type"] == "MultiPolygon") {
      auto building = city_object();

      building              = assign_type(building, properties);
      string footprint_type = building.type;

      if (footprint_type == "null") {
        continue;
      }

      int level = generate_building_level(footprint_type, properties);

      vector<array<double, 2>>         couple;
      vector<vector<array<double, 2>>> list_holes;
      vector<array<double, 2>>         list_coordinates = {};

      for (auto& multi_pol : geometry["coordinates"]) {
        int num_lists = multi_pol.size();
        for (auto& list_coords : multi_pol) {
          // std::cout << geometry["coordinates"].size() << std::endl;

          if (count_list == 0) {  // outer polygon
            building.level = level;

            for (auto& coord : list_coords) {
              double x = (double)coord[0];
              double y = (double)coord[1];
              class_coord.update(x, y);
              class_coord.update(x, y);
              array<double, 2> arr = {x, y};
              list_coordinates.push_back(arr);
            }

            building.coords = list_coordinates;

            string name   = id;
            building.name = "building_" + name;

            count_list++;
          } else {  // analysis of building holes

            for (auto& coord : list_coords) {
              double x = (double)coord[0];
              double y = (double)coord[1];
              couple.push_back({x, y});
            }
            list_holes.push_back(couple);
            couple = {};

            count_list++;
          }

          if (count_list == num_lists) {
            building.holes  = list_holes;
            bool valid_type = check_valid_type(building, properties);
            if (valid_type) all_buildings.push_back(building);
          }
        }
        count_list = 0;
      }

    } else if (geometry["type"] == "LineString") {
      int cont = 0;
      for (int i = 0; i < geometry["coordinates"].size() - 1; i++) {
        auto coord_i      = geometry["coordinates"][i];
        auto coord_i_next = geometry["coordinates"][i + 1];

        double x              = (double)coord_i[0];
        double y              = (double)coord_i[1];
        double next_x         = (double)coord_i_next[0];
        double next_y         = (double)coord_i_next[1];
        double line_thickness = 0.00005f;

        auto area = compute_area(x, next_x, y, next_y, line_thickness);

        auto line = city_object();

        string name = id;
        line.name   = "line_" + name + std::to_string(cont);
        cont++;

        if (properties.contains("highway")) {
          bool pedestrian = check_pedestrian(properties);
          if (pedestrian) {
            line.type = "pedestrian";
          } else {
            line.type = "highway";
          }
        }

        else if (properties.contains("natural")) {
          string natural = properties["natural"];
          line.type      = natural;
        } else {
          continue;
        }

        line.thickness = get_thickness(line.type);
        line.coords    = area;
        all_buildings.push_back(line);

        for (auto& coord : area) {
          double x = (double)coord[0];
          double y = (double)coord[1];
          class_coord.update(x, y);
        }
      }

    } else if (geometry["type"] == "MultiLineString") {
      int cont = 0;
      for (auto& list_line : geometry["coordinates"]) {
        for (int i = 0; i < list_line.size() - 1; i++) {
          auto coord_i      = list_line[i];
          auto coord_i_next = list_line[i + 1];

          double x              = (double)coord_i[0];
          double y              = (double)coord_i[1];
          double next_x         = (double)coord_i_next[0];
          double next_y         = (double)coord_i_next[1];
          double line_thickness = 0.0004f;

          auto area = compute_area(x, next_x, y, next_y, line_thickness);

          auto line = city_object();

          string name = id;
          line.name   = "multiline_" + name + std::to_string(cont);
          cont++;

          if (properties.contains("waterway")) {
            line.type = "water";
          } else {
            continue;
          }

          line.thickness = line_thickness;
          line.coords    = area;
          all_buildings.push_back(line);

          for (auto& coord : area) {
            double x = (double)coord[0];
            double y = (double)coord[1];
            class_coord.update(x, y);
          }
        }
      }

    } else if (geometry["type"] == "Point") {
      vector<array<double, 2>> points;
      points.push_back(geometry["coordinates"]);

      auto point = city_object();

      string name  = id;
      point.name   = "point_" + name;
      point.coords = points;

      all_buildings = assign_tree_type(point, properties, all_buildings);

      if (point.type == "null") {
        continue;
      }

      for (auto& coord : points) {
        double x = (double)coord[0];
        double y = (double)coord[1];
        class_coord.update(x, y);
      }
    }
  }
  // std::cout << all_buildings << std::endl;
  return {all_buildings, class_coord};
}

vector<city_object> generate_new_coordinates(json geojson_file,
    vector<city_object> all_buildings, Coordinate class_coord) {
  std::pair<vector<city_object>, Coordinate> gen_out;

  gen_out       = data_analysis(geojson_file, all_buildings, class_coord);
  all_buildings = gen_out.first;
  class_coord   = gen_out.second;

  vector<city_object> all_objects = {};

  // Scale the CityObject outer polygon in the scene
  for (auto& building_geometry : all_buildings) {
    float height             = generate_height(building_geometry, scale);
    building_geometry.height = height;

    vector<array<double, 2>>         new_coords = {};
    vector<vector<array<double, 2>>> new_holes  = {};

    for (auto& couple : building_geometry.coords) {
      double x = (double)couple[0];
      double y = (double)couple[1];

      double new_x = (x - class_coord.x_minimum) /
                         (class_coord.x_maximum - class_coord.x_minimum) *
                         scale -
                     (scale / 2);
      double new_y = (y - class_coord.y_minimum) /
                         (class_coord.y_maximum - class_coord.y_minimum) *
                         scale -
                     (scale / 2);

      array<double, 2> arr = {new_x, new_y};
      new_coords.push_back(arr);
    }
    building_geometry.new_coords = new_coords;

    // Scale the CityObject holes in the scene
    for (auto& list_hole : building_geometry.holes) {
      vector<array<double, 2>> new_hole_l = {};
      for (auto& hole : list_hole) {
        double x = (double)hole[0];
        double y = (double)hole[1];

        double new_hole_x =
            (x - class_coord.x_minimum) /
                (class_coord.x_maximum - class_coord.x_minimum) * scale -
            (scale / 2);
        double new_hole_y =
            (y - class_coord.y_minimum) /
                (class_coord.y_maximum - class_coord.y_minimum) * scale -
            (scale / 2);
        array<double, 2> new_hole_coords = {new_hole_x, new_hole_y};
        new_hole_l.push_back(new_hole_coords);
      }
      new_holes.push_back(new_hole_l);
    }
    building_geometry.new_holes = new_holes;

    // std::cout << "------" << std::endl;
    // std::cout << building_geometry["new_holes"] << std::endl;
    all_objects.push_back(building_geometry);
  }
  // std::cout << all_objects << std::endl;
  return all_objects;
}

//  ---------------- MAIN FUNCTION --------------------------

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
  add_option(cli, "--info,-i", info, "print scene info");
  add_option(cli, "--copyright,-c", copyright, "copyright string");
  add_option(cli, "--validate/--no-validate", validate, "Validate scene");
  add_option(cli, "--output,-o", output, "output scene");
  add_option(cli, "dirname", path, "input directory", true);
  parse_cli(cli, argc, argv);

  // load data

  // read GeoJson files
  auto buildings   = vector<city_object>{};
  auto class_coord = Coordinate();
  auto ioerror     = ""s;
  print_progress("load geojsons", 0, 1);
  for (auto& filename : list_directory(path)) {
    if (path_extension(filename) != ".geojson") continue;
    auto geojson = json{};
    if (!load_json(filename, geojson, ioerror)) print_fatal(ioerror);
    buildings = generate_new_coordinates(geojson, buildings, class_coord);
  }
  print_progress("load geojsons", 1, 1);

  // Create city
  auto scene_guard = std::make_unique<sceneio_scene>();
  auto scene       = scene_guard.get();
  print_progress("convert scene", 0, 1);
  if (!create_city_from_json(scene, buildings, path, ioerror))
    print_fatal(ioerror);
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
  if (!scene->shapes.empty()) {
    if (!make_directory(path_join(path_dirname(output), "shapes"), ioerror))
      print_fatal(ioerror);
  }
  if (!scene->textures.empty()) {
    if (!make_directory(path_join(path_dirname(output), "textures"), ioerror))
      print_fatal(ioerror);
  }

  // save scene
  if (!save_scene(output, scene, ioerror, print_progress)) print_fatal(ioerror);

  // Done
  return 0;
}
