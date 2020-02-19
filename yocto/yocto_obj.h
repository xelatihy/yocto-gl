//
// # Yocto/Obj: Tiny library for Obj parsing and writing
//
// Yocto/Obj is a tiny library for loading and saving Obj.
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::obj {

// Using directives
using namespace yocto::math;

// OBJ vertex
struct obj_vertex {
  int position = 0;
  int texcoord = 0;
  int normal   = 0;
};

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
  return a.position == b.position && a.texcoord == b.texcoord &&
         a.normal == b.normal;
}

// Obj texture information.
struct obj_texture_info {
  string path  = "";     // file path
  bool   clamp = false;  // clamp to edge
  float  scale = 1;      // scale for bump/displacement

  // Properties not explicitly handled.
  unordered_map<string, vector<float>> props;

  obj_texture_info() {}
  obj_texture_info(const char* path) : path{path} {}
  obj_texture_info(const string& path) : path{path} {}
};

// Obj element
struct obj_element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj material
struct obj_material {
  // material name and type
  string name  = "";
  int    illum = 0;

  // material colors and values
  vec3f emission     = zero3f;
  vec3f ambient      = zero3f;
  vec3f diffuse      = zero3f;
  vec3f specular     = zero3f;
  vec3f reflection   = zero3f;
  vec3f transmission = zero3f;
  float exponent     = 10;
  float ior          = 1.5;
  float opacity      = 1;

  // material textures
  obj_texture_info emission_tex     = {};
  obj_texture_info ambient_tex      = {};
  obj_texture_info diffuse_tex      = {};
  obj_texture_info specular_tex     = {};
  obj_texture_info reflection_tex   = {};
  obj_texture_info transmission_tex = {};
  obj_texture_info exponent_tex     = {};
  obj_texture_info opacity_tex      = {};
  obj_texture_info bump_tex         = {};
  obj_texture_info normal_tex       = {};
  obj_texture_info displacement_tex = {};

  // pbrt extension values
  bool  as_pbr            = false;
  vec3f pbr_emission      = {0, 0, 0};
  vec3f pbr_base          = {0, 0, 0};
  float pbr_specular      = 0;
  float pbr_roughness     = 0;
  float pbr_metallic      = 0;
  float pbr_sheen         = 0;
  float pbr_coat          = 0;
  float pbr_coatroughness = 0;
  float pbr_transmission  = 0;
  float pbr_ior           = 1.5;
  float pbr_opacity       = 1;
  vec3f pbr_volscattering = zero3f;
  float pbr_volanisotropy = 0;
  float pbr_volscale      = 0.01;

  // pbr extension textures
  obj_texture_info pbr_emission_tex      = {};
  obj_texture_info pbr_base_tex          = {};
  obj_texture_info pbr_specular_tex      = {};
  obj_texture_info pbr_roughness_tex     = {};
  obj_texture_info pbr_metallic_tex      = {};
  obj_texture_info pbr_sheen_tex         = {};
  obj_texture_info pbr_coat_tex          = {};
  obj_texture_info pbr_coatroughness_tex = {};
  obj_texture_info pbr_transmission_tex  = {};
  obj_texture_info pbr_opacity_tex       = {};
  obj_texture_info pbr_volscattering_tex = {};
};

// Obj shape
struct obj_shape {
  string                name      = "";
  vector<vec3f>         positions = {};
  vector<vec3f>         normals   = {};
  vector<vec2f>         texcoords = {};
  vector<obj_material*> materials = {};
  vector<obj_vertex>    vertices  = {};
  vector<obj_element>   faces     = {};
  vector<obj_element>   lines     = {};
  vector<obj_element>   points    = {};
  vector<frame3f>       instances = {};
};

// Obj camera
struct obj_camera {
  string  name     = "";
  frame3f frame    = identity3x4f;
  bool    ortho    = false;
  float   width    = 0.036;
  float   height   = 0.028;
  float   lens     = 0.50;
  float   focus    = 0;
  float   aperture = 0;
};

// Obj environment
struct obj_environment {
  string           name         = "";
  frame3f          frame        = identity3x4f;
  vec3f            emission     = zero3f;
  obj_texture_info emission_tex = {};
};

// Obj model
struct obj_model {
  vector<string>           comments     = {};
  vector<obj_shape*>       shapes       = {};
  vector<obj_material*>    materials    = {};
  vector<obj_camera*>      cameras      = {};
  vector<obj_environment*> environments = {};
  ~obj_model();
};

// Load and save obj
inline bool load_obj(const string& filename, obj_model* obj, string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
inline bool save_obj(const string& filename, obj_model* obj, string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
inline void get_triangles(obj_model* obj, obj_shape* shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<obj_material*>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline void get_quads(obj_model* obj, obj_shape* shape, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
inline void get_lines(obj_model* obj, obj_shape* shape, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
inline void get_points(obj_model* obj, obj_shape* shape, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
inline void get_fvquads(obj_model* obj, obj_shape* shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
inline bool has_quads(obj_shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
inline void get_triangles(obj_model* obj, obj_shape* shape, int material,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
inline void get_quads(obj_model* obj, obj_shape* shape, int material,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
inline void get_lines(obj_model* obj, obj_shape* shape, int material,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
inline void get_points(obj_model* obj, obj_shape* shape, int material,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
inline vector<obj_material*> get_materials(obj_model* obj, obj_shape* shape);

// Create OBJ
inline obj_camera*      add_camera(obj_model* obj);
inline obj_material*    add_material(obj_model* obj);
inline obj_environment* add_environment(obj_model* obj);
inline obj_shape*       add_shape(obj_model* obj);

// Add obj shape
inline void add_triangles(obj_model* obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
inline void add_quads(obj_model* obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
inline void add_lines(obj_model* obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
inline void add_points(obj_model* obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
inline void add_fvquads(obj_model* obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);

}  // namespace yocto::obj

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::obj::obj_vertex> {
  size_t operator()(const yocto::obj::obj_vertex& v) const {
    const std::hash<int> hasher = std::hash<int>();
    auto                 h      = (size_t)0;
    h ^= hasher(v.position) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.normal) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.texcoord) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#include <cstdio>
#include <memory>
#include <string_view>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::obj {

using std::string_view;

// utilities
inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

// Parse values from a string
[[nodiscard]] inline bool parse_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) return false;
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
    if (str.front() != '"') return false;
    str.remove_prefix(1);
    if (str.empty()) return false;
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
    if (cpy.empty()) return false;
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    str.remove_prefix(1);
  }
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  if (!parse_value(str, valuei)) return false;
  value = (bool)valuei;
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}

[[nodiscard]] inline bool parse_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

// Formats values to string
inline void format_value(string& str, const string& value) { str += value; }
inline void format_value(string& str, int value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}

// Foramt to file
inline void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
[[nodiscard]] inline bool format_values(
    FILE* fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  if (fputs(str.c_str(), fs) < 0) return false;
  return true;
}
template <typename T>
[[nodiscard]] inline bool format_value(FILE* fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  if (fputs(str.c_str(), fs) < 0) return false;
  return true;
}

inline void remove_obj_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

[[nodiscard]] inline bool parse_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  if (!parse_value(str, value.position)) return false;
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      if (!parse_value(str, value.normal)) return false;
    } else {
      if (!parse_value(str, value.texcoord)) return false;
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        if (!parse_value(str, value.normal)) return false;
      }
    }
  }
  return true;
}

// Input for OBJ textures
[[nodiscard]] inline bool parse_value(
    string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    if (!parse_value(str, token)) return false;
    tokens.push_back(token);
    skip_whitespace(str);
  }
  if (tokens.empty()) return false;

  // texture name
  info.path = tokens.back();
  for (auto& c : info.path)
    if (c == '\\') c = '/';

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }

  return true;
}

// Read obj
[[nodiscard]] inline bool load_mtl(
    const string& filename, obj_model* obj, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // init parsing
  auto material = obj->materials.emplace_back(new obj_material{});

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // grab material
    material = obj->materials.back();

    // possible token values
    if (cmd == "newmtl") {
      auto material = obj->materials.emplace_back(new obj_material{});
      if (!parse_value(str, material->name)) return parse_error();
    } else if (cmd == "illum") {
      if (!parse_value(str, material->illum)) return parse_error();
    } else if (cmd == "Ke") {
      if (!parse_value(str, material->emission)) return parse_error();
    } else if (cmd == "Ka") {
      if (!parse_value(str, material->ambient)) return parse_error();
    } else if (cmd == "Kd") {
      if (!parse_value(str, material->diffuse)) return parse_error();
    } else if (cmd == "Ks") {
      if (!parse_value(str, material->specular)) return parse_error();
    } else if (cmd == "Kt") {
      if (!parse_value(str, material->transmission)) return parse_error();
    } else if (cmd == "Tf") {
      if (!parse_value(str, material->transmission)) return parse_error();
      material->transmission = max(1 - material->transmission, 0.0f);
      if (max(material->transmission) < 0.001) material->transmission = zero3f;
    } else if (cmd == "Tr") {
      if (!parse_value(str, material->opacity)) return parse_error();
      material->opacity = 1 - material->opacity;
    } else if (cmd == "Ns") {
      if (!parse_value(str, material->exponent)) return parse_error();
    } else if (cmd == "d") {
      if (!parse_value(str, material->opacity)) return parse_error();
    } else if (cmd == "map_Ke") {
      if (!parse_value(str, material->emission_tex)) return parse_error();
    } else if (cmd == "map_Ka") {
      if (!parse_value(str, material->ambient_tex)) return parse_error();
    } else if (cmd == "map_Kd") {
      if (!parse_value(str, material->diffuse_tex)) return parse_error();
    } else if (cmd == "map_Ks") {
      if (!parse_value(str, material->specular_tex)) return parse_error();
    } else if (cmd == "map_Tr") {
      if (!parse_value(str, material->transmission_tex)) return parse_error();
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      if (!parse_value(str, material->opacity_tex)) return parse_error();
    } else if (cmd == "map_bump" || cmd == "bump") {
      if (!parse_value(str, material->bump_tex)) return parse_error();
    } else if (cmd == "map_disp" || cmd == "disp") {
      if (!parse_value(str, material->displacement_tex)) return parse_error();
    } else if (cmd == "map_norm" || cmd == "norm") {
      if (!parse_value(str, material->normal_tex)) return parse_error();
    } else if (cmd == "Pe") {
      if (!parse_value(str, material->pbr_emission)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pb") {
      if (!parse_value(str, material->pbr_base)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Ps") {
      if (!parse_value(str, material->pbr_specular)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pm") {
      if (!parse_value(str, material->pbr_metallic)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pr") {
      if (!parse_value(str, material->pbr_roughness)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Ps") {
      if (!parse_value(str, material->pbr_sheen)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pc") {
      if (!parse_value(str, material->pbr_coat)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pcr") {
      if (!parse_value(str, material->pbr_coatroughness)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pt") {
      if (!parse_value(str, material->pbr_transmission)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pn") {
      if (!parse_value(str, material->pbr_ior)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Po") {
      if (!parse_value(str, material->pbr_opacity)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvs") {
      if (!parse_value(str, material->pbr_volscattering)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvg") {
      if (!parse_value(str, material->pbr_volanisotropy)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvr") {
      if (!parse_value(str, material->pbr_volscale)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pe") {
      if (!parse_value(str, material->pbr_emission_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pb") {
      if (!parse_value(str, material->pbr_base_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Ps") {
      if (!parse_value(str, material->pbr_specular_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pm") {
      if (!parse_value(str, material->pbr_metallic_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pr") {
      if (!parse_value(str, material->pbr_roughness_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Ps") {
      if (!parse_value(str, material->pbr_sheen_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pc") {
      if (!parse_value(str, material->pbr_coat_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pcr") {
      if (!parse_value(str, material->pbr_coatroughness_tex))
        return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Po") {
      if (!parse_value(str, material->pbr_opacity_tex)) return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pt") {
      if (!parse_value(str, material->pbr_transmission_tex))
        return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Vs") {
      if (!parse_value(str, material->pbr_volscattering_tex))
        return parse_error();
      material->as_pbr = true;
    } else {
      continue;
    }
  }

  // remove placeholder material
  obj->materials.erase(obj->materials.begin());

  // convert between roughness and exponent
  auto exponent_to_roughness = [](float exponent) {
    auto roughness = exponent;
    roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
    if (roughness < 0.01f) roughness = 0;
    if (roughness > 0.99f) roughness = 1;
    return roughness;
  };

  // convert values when possible
  for (auto material : obj->materials) {
    if (material->as_pbr) continue;
    material->pbr_emission     = material->emission;
    material->pbr_emission_tex = material->emission_tex;
    material->pbr_roughness    = exponent_to_roughness(material->exponent);
    material->pbr_ior          = material->ior;
    material->pbr_opacity      = material->opacity;
    material->pbr_opacity_tex  = material->opacity_tex;
    if (max(material->transmission) > 0.1) {
      material->pbr_base         = material->transmission;
      material->pbr_transmission = 1;
      material->pbr_specular     = 1;
    } else if (max(material->specular) > 0.2) {
      material->pbr_base     = material->specular;
      material->pbr_base_tex = material->specular_tex;
      material->pbr_metallic = 1;
    } else {
      material->pbr_base     = material->diffuse;
      material->pbr_base_tex = material->diffuse_tex;
      material->pbr_specular = max(material->specular) ? 1 : 0;
    }
  }

  return true;
}

// Read obj
[[nodiscard]] inline bool load_objx(
    const string& filename, obj_model* obj, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // shape map for instances
  auto shape_map = unordered_map<string, vector<obj_shape*>>{};
  for (auto shape : obj->shapes) {
    shape_map[shape->name].push_back(shape);
  }

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto camera = obj->cameras.emplace_back(new obj_camera{});
      if (!parse_value(str, camera->name)) return parse_error();
      if (!parse_value(str, camera->ortho)) return parse_error();
      if (!parse_value(str, camera->width)) return parse_error();
      if (!parse_value(str, camera->height)) return parse_error();
      if (!parse_value(str, camera->lens)) return parse_error();
      if (!parse_value(str, camera->focus)) return parse_error();
      if (!parse_value(str, camera->aperture)) return parse_error();
      if (!parse_value(str, camera->frame)) return parse_error();
    } else if (cmd == "e") {
      auto environment = obj->environments.emplace_back(new obj_environment{});
      if (!parse_value(str, environment->name)) return parse_error();
      if (!parse_value(str, environment->emission)) return parse_error();
      auto emission_path = ""s;
      if (!parse_value(str, emission_path)) return parse_error();
      if (emission_path == "\"\"") emission_path = "";
      environment->emission_tex.path = emission_path;
      if (!parse_value(str, environment->frame)) return parse_error();
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      if (!parse_value(str, object)) return parse_error();
      if (!parse_value(str, frame)) return parse_error();
      if (shape_map.find(object) == shape_map.end()) {
        return parse_error();
      }
      for (auto shape : shape_map.at(object)) {
        shape->instances.push_back(frame);
      }
    } else {
      // unused
    }
  }

  return true;
}

inline obj_model::~obj_model() {
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto camera : cameras) delete camera;
  for (auto environment : environments) delete environment;
}

// Make obj
inline obj_camera* add_camera(obj_model* obj) {
  return obj->cameras.emplace_back(new obj_camera{});
}
inline obj_material* add_material(obj_model* obj) {
  return obj->materials.emplace_back(new obj_material{});
}
inline obj_environment* add_environment(obj_model* obj) {
  return obj->environments.emplace_back(new obj_environment{});
}
inline obj_shape* add_shape(obj_model* obj) {
  return obj->shapes.emplace_back(new obj_shape{});
}

// Read obj
inline bool load_obj(const string& filename, obj_model* obj, string& error,
    bool geom_only, bool split_elements, bool split_materials) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto read_error = [filename, &error]() {
    error = filename + ": read error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // parsing state
  auto opositions   = vector<vec3f>{};
  auto onormals     = vector<vec3f>{};
  auto otexcoords   = vector<vec2f>{};
  auto vert_size    = obj_vertex{};
  auto oname        = ""s;
  auto gname        = ""s;
  auto mname        = ""s;
  auto mtllibs      = vector<string>{};
  auto material_map = unordered_map<string, obj_material*>{};

  // initialize obj
  obj->~obj_model();
  obj->cameras.clear();
  obj->environments.clear();
  obj->shapes.clear();
  obj->materials.clear();

  // initialize load
  obj->shapes.emplace_back(new obj_shape{});
  auto empty_material = (obj_material*)nullptr;

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      if (!parse_value(str, opositions.emplace_back(zero3f)))
        return parse_error();
      vert_size.position += 1;
    } else if (cmd == "vn") {
      if (!parse_value(str, onormals.emplace_back(zero3f)))
        return parse_error();
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      if (!parse_value(str, otexcoords.emplace_back(zero2f)))
        return parse_error();
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      // split if split_elements and different primitives
      if (auto shape = obj->shapes.back();
          split_elements && !shape->vertices.empty()) {
        if ((cmd == "f" && (!shape->lines.empty() || !shape->points.empty())) ||
            (cmd == "l" && (!shape->faces.empty() || !shape->points.empty())) ||
            (cmd == "p" && (!shape->faces.empty() || !shape->lines.empty()))) {
          obj->shapes.emplace_back(new obj_shape{});
          obj->shapes.back()->name = oname + gname;
        }
      }
      // split if splt_material and different materials
      if (auto shape = obj->shapes.back();
          !geom_only && split_materials && !shape->materials.empty()) {
        if (shape->materials.size() > 1)
          throw std::runtime_error("should not have happened");
        if (shape->materials.back()->name != mname) {
          obj->shapes.emplace_back(new obj_shape{});
          obj->shapes.back()->name = oname + gname;
        }
      }
      // grab shape and add element
      auto  shape   = obj->shapes.back();
      auto& element = (cmd == "f")
                          ? shape->faces.emplace_back()
                          : (cmd == "l") ? shape->lines.emplace_back()
                                         : shape->points.emplace_back();
      // get element material or add if needed
      if (!geom_only) {
        if (mname.empty() && !empty_material) {
          empty_material   = obj->materials.emplace_back(new obj_material{});
          material_map[""] = empty_material;
        }
        auto mat_idx = -1;
        for (auto midx = 0; midx < shape->materials.size(); midx++)
          if (shape->materials[midx]->name == mname) mat_idx = midx;
        if (mat_idx < 0) {
          shape->materials.push_back(material_map.at(mname));
          mat_idx = shape->materials.size() - 1;
        }
        element.material = (uint8_t)mat_idx;
      }
      // parse vertices
      skip_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        if (!parse_value(str, vert)) return parse_error();
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        shape->vertices.push_back(vert);
        element.size += 1;
        skip_whitespace(str);
      }
    } else if (cmd == "o" || cmd == "g") {
      if (geom_only) continue;
      skip_whitespace(str);
      if (cmd == "o") {
        if (str.empty()) {
          oname = "";
        } else {
          if (!parse_value(str, oname)) return parse_error();
        }
      } else {
        if (str.empty()) {
          gname = "";
        } else {
          if (!parse_value(str, gname)) return parse_error();
        }
      }
      if (!obj->shapes.back()->vertices.empty()) {
        obj->shapes.emplace_back(new obj_shape{});
        obj->shapes.back()->name = oname + gname;
      } else {
        obj->shapes.back()->name = oname + gname;
      }
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      if (!parse_value(str, mname)) return parse_error();
    } else if (cmd == "s") {
      if (geom_only) continue;
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      if (!parse_value(str, mtllib)) return parse_error();
      if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) == mtllibs.end()) {
        mtllibs.push_back(mtllib);
        if (!load_mtl(fs::path(filename).parent_path() / mtllib, obj, error))
          return dependent_error();
        for (auto material : obj->materials)
          material_map[material->name] = material;
      }
    } else {
      // unused
    }
  }

  // fix empty material
  if (empty_material) {
    empty_material->name     = "empty_material";
    empty_material->diffuse  = {0.8, 0.8, 0.8};
    empty_material->pbr_base = {0.8, 0.8, 0.8};
  }

  // convert vertex data
  auto ipositions = vector<int>{};
  auto inormals   = vector<int>{};
  auto itexcoords = vector<int>{};
  for (auto shape : obj->shapes) {
    ipositions.assign(opositions.size() + 1, 0);
    inormals.assign(onormals.size() + 1, 0);
    itexcoords.assign(otexcoords.size() + 1, 0);
    for (auto& vertex : shape->vertices) {
      if (vertex.position && !ipositions[vertex.position]) {
        shape->positions.push_back(opositions[vertex.position - 1]);
        ipositions[vertex.position] = (int)shape->positions.size();
      }
      if (vertex.normal && !inormals[vertex.normal]) {
        shape->normals.push_back(onormals[vertex.normal - 1]);
        inormals[vertex.normal] = (int)shape->normals.size();
      }
      if (vertex.texcoord && !itexcoords[vertex.texcoord]) {
        shape->texcoords.push_back(otexcoords[vertex.texcoord - 1]);
        itexcoords[vertex.texcoord] = (int)shape->texcoords.size();
      }
      vertex.position = ipositions[vertex.position];
      vertex.normal   = inormals[vertex.normal];
      vertex.texcoord = itexcoords[vertex.texcoord];
    }
  }

  // exit if done
  if (geom_only) return true;

  // load extensions
  auto extfilename = fs::path(filename).replace_extension(".objx");
  if (fs::exists(fs::path(extfilename))) {
    if (!load_objx(extfilename, obj, error)) return dependent_error();
  }

  return true;
}

// Format values
inline void format_value(string& str, const obj_texture_info& value) {
  str += value.path.empty() ? "" : value.path;
}
inline void format_value(string& str, const obj_vertex& value) {
  format_value(str, value.position);
  if (value.texcoord) {
    str += "/";
    format_value(str, value.texcoord);
    if (value.normal) {
      str += "/";
      format_value(str, value.normal);
    }
  } else if (value.normal) {
    str += "//";
    format_value(str, value.normal);
  }
}

// Save obj
[[nodiscard]] inline bool save_mtl(
    const string& filename, obj_model* obj, string& error) {
  // throw helpers
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // save comments
  if (!format_values(fs, "#\n")) return write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) return write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return write_error();
  if (!format_values(fs, "#\n\n")) return write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) return write_error();
  }
  if (!format_values(fs, "\n")) return write_error();

  // write material
  for (auto material : obj->materials) {
    if (!format_values(fs, "newmtl {}\n", material->name)) return write_error();
    if (!material->as_pbr) {
      if (!format_values(fs, "illum {}\n", material->illum))
        return write_error();
      if (material->emission != zero3f)
        if (!format_values(fs, "Ke {}\n", material->emission))
          return write_error();
      if (material->ambient != zero3f)
        if (!format_values(fs, "Ka {}\n", material->ambient))
          return write_error();
      if (!format_values(fs, "Kd {}\n", material->diffuse))
        return write_error();
      if (!format_values(fs, "Ks {}\n", material->specular))
        return write_error();
      if (material->reflection != zero3f)
        if (!format_values(fs, "Kr {}\n", material->reflection))
          return write_error();
      if (material->transmission != zero3f)
        if (!format_values(fs, "Kt {}\n", material->transmission))
          return write_error();
      if (!format_values(fs, "Ns {}\n", (int)material->exponent))
        return write_error();
      if (material->opacity != 1)
        if (!format_values(fs, "d {}\n", material->opacity))
          return write_error();
      if (!material->emission_tex.path.empty())
        if (!format_values(fs, "map_Ke {}\n", material->emission_tex))
          return write_error();
      if (!material->diffuse_tex.path.empty())
        if (!format_values(fs, "map_Kd {}\n", material->diffuse_tex))
          return write_error();
      if (!material->specular_tex.path.empty())
        if (!format_values(fs, "map_Ks {}\n", material->specular_tex))
          return write_error();
      if (!material->transmission_tex.path.empty())
        if (!format_values(fs, "map_Kt {}\n", material->transmission_tex))
          return write_error();
      if (!material->reflection_tex.path.empty())
        if (!format_values(fs, "map_Kr {}\n", material->reflection_tex))
          return write_error();
      if (!material->exponent_tex.path.empty())
        if (!format_values(fs, "map_Ns {}\n", material->exponent_tex))
          return write_error();
      if (!material->opacity_tex.path.empty())
        if (!format_values(fs, "map_d {}\n", material->opacity_tex))
          return write_error();
      if (!material->bump_tex.path.empty())
        if (!format_values(fs, "map_bump {}\n", material->bump_tex))
          return write_error();
      if (!material->displacement_tex.path.empty())
        if (!format_values(fs, "map_disp {}\n", material->displacement_tex))
          return write_error();
      if (!material->normal_tex.path.empty())
        if (!format_values(fs, "map_norm {}\n", material->normal_tex))
          return write_error();
    } else {
      if (!format_values(fs, "illum 2\n")) return write_error();
      if (material->pbr_emission != zero3f)
        if (!format_values(fs, "Pe {}\n", material->pbr_emission))
          return write_error();
      if (material->pbr_base != zero3f)
        if (!format_values(fs, "Pb {}\n", material->pbr_base))
          return write_error();
      if (material->pbr_specular)
        if (!format_values(fs, "Psp {}\n", material->pbr_specular))
          return write_error();
      if (material->pbr_roughness)
        if (!format_values(fs, "Pr {}\n", material->pbr_roughness))
          return write_error();
      if (material->pbr_metallic)
        if (!format_values(fs, "Pm {}\n", material->pbr_metallic))
          return write_error();
      if (material->pbr_sheen)
        if (!format_values(fs, "Ps {}\n", material->pbr_sheen))
          return write_error();
      if (material->pbr_coat)
        if (!format_values(fs, "Pc {}\n", material->pbr_coat))
          return write_error();
      if (material->pbr_coatroughness)
        if (!format_values(fs, "Pcr {}\n", material->pbr_coatroughness))
          return write_error();
      if (material->pbr_volscattering != zero3f)
        if (!format_values(fs, "Pvs {}\n", material->pbr_volscattering))
          return write_error();
      if (material->pbr_volanisotropy)
        if (!format_values(fs, "Pvg {}\n", material->pbr_volanisotropy))
          return write_error();
      if (material->pbr_volscale)
        if (!format_values(fs, "Pvr {}\n", material->pbr_volscale))
          return write_error();
      if (!material->pbr_emission_tex.path.empty())
        if (!format_values(fs, "map_Pe {}\n", material->pbr_emission_tex))
          return write_error();
      if (!material->pbr_base_tex.path.empty())
        if (!format_values(fs, "map_Pb {}\n", material->pbr_base_tex))
          return write_error();
      if (!material->pbr_specular_tex.path.empty())
        if (!format_values(fs, "map_Psp {}\n", material->pbr_specular_tex))
          return write_error();
      if (!material->pbr_roughness_tex.path.empty())
        if (!format_values(fs, "map_Pr {}\n", material->pbr_roughness_tex))
          return write_error();
      if (!material->pbr_metallic_tex.path.empty())
        if (!format_values(fs, "map_Pm {}\n", material->pbr_metallic_tex))
          return write_error();
      if (!material->pbr_sheen_tex.path.empty())
        if (!format_values(fs, "map_Ps {}\n", material->pbr_sheen_tex))
          return write_error();
      if (!material->pbr_coat_tex.path.empty())
        if (!format_values(fs, "map_Pc {}\n", material->pbr_coat_tex))
          return write_error();
      if (!material->pbr_coatroughness_tex.path.empty())
        if (!format_values(fs, "map_Pcr {}\n", material->pbr_coatroughness_tex))
          return write_error();
      if (!material->pbr_volscattering_tex.path.empty())
        if (!format_values(fs, "map_Pvs {}\n", material->pbr_volscattering_tex))
          return write_error();
      if (!material->bump_tex.path.empty())
        if (!format_values(fs, "map_bump {}\n", material->bump_tex))
          return write_error();
      if (!material->displacement_tex.path.empty())
        if (!format_values(fs, "map_disp {}\n", material->displacement_tex))
          return write_error();
      if (!material->normal_tex.path.empty())
        if (!format_values(fs, "map_norm {}\n", material->normal_tex))
          return write_error();
    }
    if (!format_values(fs, "\n")) return write_error();
  }
  return true;
}

// Save obj
[[nodiscard]] inline bool save_objx(
    const string& filename, obj_model* obj, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // save comments
  if (!format_values(fs, "#\n")) return write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) return write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return write_error();
  if (!format_values(fs, "#\n\n")) return write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) return write_error();
  }
  if (!format_values(fs, "\n")) return write_error();

  // cameras
  for (auto camera : obj->cameras) {
    if (!format_values(fs, "c {} {} {} {} {} {} {} {}\n", camera->name,
            camera->ortho, camera->width, camera->height, camera->lens,
            camera->focus, camera->aperture, camera->frame))
      return write_error();
  }

  // environments
  for (auto environment : obj->environments) {
    if (!format_values(fs, "e {} {} {} {}\n", environment->name,
            environment->emission,
            environment->emission_tex.path.empty()
                ? "\"\""s
                : environment->emission_tex.path,
            environment->frame))
      return write_error();
  }

  // instances
  for (auto shape : obj->shapes) {
    for (auto& frame : shape->instances) {
      if (!format_values(fs, "i {} {}\n", shape->name, frame))
        return write_error();
    }
  }

  // done
  return true;
}

// Save obj
[[nodiscard]] inline bool save_obj(
    const string& filename, obj_model* obj, string& error) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // save comments
  if (!format_values(fs, "#\n")) return write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) return write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return write_error();
  if (!format_values(fs, "#\n\n")) return write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) return write_error();
  }
  if (!format_values(fs, "\n")) return write_error();

  // save material library
  if (!obj->materials.empty()) {
    if (!format_values(fs, "mtllib {}\n\n",
            fs::path(filename).filename().replace_extension(".mtl")))
      return write_error();
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto shape : obj->shapes) {
    if (!format_values(fs, "o {}\n", shape->name)) return write_error();
    for (auto& p : shape->positions)
      if (!format_values(fs, "v {}\n", p)) return write_error();
    for (auto& n : shape->normals)
      if (!format_values(fs, "vn {}\n", n)) return write_error();
    for (auto& t : shape->texcoords)
      if (!format_values(fs, "vt {}\n", t)) return write_error();
    auto element_labels = vector<string>{"f", "l", "p"};
    auto element_groups = vector<const vector<obj_element>*>{
        &shape->faces, &shape->lines, &shape->points};
    for (auto element_idx = 0; element_idx < 3; element_idx++) {
      auto& label        = element_labels[element_idx];
      auto& elements     = *element_groups[element_idx];
      auto  cur_material = -1, cur_vertex = 0;
      for (auto& element : elements) {
        if (!shape->materials.empty() && cur_material != element.material) {
          if (!format_values(
                  fs, "usemtl {}\n", shape->materials[element.material]->name))
            return write_error();
          cur_material = element.material;
        }
        if (!format_values(fs, "{}", label)) return write_error();
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape->vertices[cur_vertex++];
          if (vert.position) vert.position += vert_size.position;
          if (vert.normal) vert.normal += vert_size.normal;
          if (vert.texcoord) vert.texcoord += vert_size.texcoord;
          if (!format_values(fs, " {}", vert)) return write_error();
        }
        if (!format_values(fs, "\n")) return write_error();
      }
    }
    if (!format_values(fs, "\n")) return write_error();
    vert_size.position += (int)shape->positions.size();
    vert_size.normal += (int)shape->normals.size();
    vert_size.texcoord += (int)shape->texcoords.size();
  }

  // save mtl
  if (!obj->materials.empty()) {
    if (!save_mtl(fs::path(filename).replace_extension(".mtl"), obj, error))
      return dependent_error();
  }

  // save objx
  if (!obj->cameras.empty() || !obj->environments.empty() ||
      std::any_of(obj->shapes.begin(), obj->shapes.end(),
          [](auto shape) { return !shape->instances.empty(); })) {
    if (!save_objx(fs::path(filename).replace_extension(".objx"), obj, error))
      return dependent_error();
  }

  // done
  return true;
}

// Get obj vertices
inline void get_vertices(obj_shape* shape, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<int>& vindex,
    bool flipv) {
  auto vmap = unordered_map<obj_vertex, int>{};
  vmap.reserve(shape->vertices.size());
  vindex.reserve(shape->vertices.size());
  for (auto& vert : shape->vertices) {
    auto it = vmap.find(vert);
    if (it != vmap.end()) {
      vindex.push_back(it->second);
      continue;
    }
    auto nverts = (int)positions.size();
    vindex.push_back(nverts);
    vmap.insert(it, {vert, nverts});
    if (!shape->positions.empty() && vert.position)
      positions.push_back(shape->positions[vert.position - 1]);
    if (!shape->normals.empty() && vert.normal)
      normals.push_back(shape->normals[vert.normal - 1]);
    if (!shape->texcoords.empty() && vert.texcoord)
      texcoords.push_back(shape->texcoords[vert.texcoord - 1]);
  }
  if (flipv) {
    for (auto& texcoord : texcoords) texcoord.y = 1 - texcoord.y;
  }
}
inline vector<vec2f> flip_obj_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get obj shape
inline void get_triangles(obj_model* obj, obj_shape* shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<obj_material*>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape->materials;
  triangles.reserve(shape->faces.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& face : shape->faces) {
    for (auto c = 2; c < face.size; c++) {
      triangles.push_back(
          {vindex[cur + 0], vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(face.material);
    }
    cur += face.size;
  }
}
inline void get_quads(obj_model* obj, obj_shape* shape, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape->materials;
  quads.reserve(shape->faces.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& face : shape->faces) {
    if (face.size == 4) {
      quads.push_back(
          {vindex[cur + 0], vindex[cur + 1], vindex[cur + 2], vindex[cur + 3]});
      if (!materials.empty()) ematerials.push_back(face.material);
    } else {
      for (auto c = 2; c < face.size; c++) {
        quads.push_back({vindex[cur + 0], vindex[cur + c - 1], vindex[cur + c],
            vindex[cur + c]});
        if (!materials.empty()) ematerials.push_back(face.material);
      }
    }
    cur += face.size;
  }
}
inline void get_lines(obj_model* obj, obj_shape* shape, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials, bool flipv) {
  if (shape->lines.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape->materials;
  lines.reserve(shape->lines.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& str : shape->lines) {
    for (auto c = 1; c < str.size; c++) {
      lines.push_back({vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(str.material);
    }
    cur += str.size;
  }
}
inline void get_points(obj_model* obj, obj_shape* shape, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials, bool flipv) {
  if (shape->points.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape->materials;
  points.reserve(shape->points.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& point : shape->points) {
    for (auto c = 0; c < point.size; c++) {
      points.push_back({vindex[cur + 0]});
      if (!materials.empty()) ematerials.push_back(point.material);
    }
    cur += point.size;
  }
}
inline void get_fvquads(obj_model* obj, obj_shape* shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  positions = shape->positions;
  normals   = shape->normals;
  texcoords = flipv ? flip_obj_texcoord(shape->texcoords) : shape->texcoords;
  materials = shape->materials;
  if (shape->vertices[0].position) quadspos.reserve(shape->faces.size());
  if (shape->vertices[0].normal) quadsnorm.reserve(shape->faces.size());
  if (shape->vertices[0].texcoord) quadstexcoord.reserve(shape->faces.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& face : shape->faces) {
    if (face.size == 4) {
      if (shape->vertices[0].position)
        quadspos.push_back({shape->vertices[cur + 0].position - 1,
            shape->vertices[cur + 1].position - 1,
            shape->vertices[cur + 2].position - 1,
            shape->vertices[cur + 3].position - 1});
      if (shape->vertices[0].normal)
        quadsnorm.push_back({shape->vertices[cur + 0].normal - 1,
            shape->vertices[cur + 1].normal - 1,
            shape->vertices[cur + 2].normal - 1,
            shape->vertices[cur + 3].normal - 1});
      if (shape->vertices[0].texcoord)
        quadstexcoord.push_back({shape->vertices[cur + 0].texcoord - 1,
            shape->vertices[cur + 1].texcoord - 1,
            shape->vertices[cur + 2].texcoord - 1,
            shape->vertices[cur + 3].texcoord - 1});
      if (!materials.empty()) ematerials.push_back(face.material);
    } else {
      for (auto c = 2; c < face.size; c++) {
        if (shape->vertices[0].position)
          quadspos.push_back({shape->vertices[cur + 0].position - 1,
              shape->vertices[cur + c - 1].position - 1,
              shape->vertices[cur + c].position - 1,
              shape->vertices[cur + c].position - 1});
        if (shape->vertices[0].normal)
          quadsnorm.push_back({shape->vertices[cur + 0].normal - 1,
              shape->vertices[cur + c - 1].normal - 1,
              shape->vertices[cur + c].normal - 1,
              shape->vertices[cur + c].normal - 1});
        if (shape->vertices[0].texcoord)
          quadstexcoord.push_back({shape->vertices[cur + 0].texcoord - 1,
              shape->vertices[cur + c - 1].texcoord - 1,
              shape->vertices[cur + c].texcoord - 1,
              shape->vertices[cur + c].texcoord - 1});
        if (!materials.empty()) ematerials.push_back(face.material);
      }
    }
    cur += face.size;
  }
}

inline bool has_quads(obj_shape* shape) {
  for (auto& face : shape->faces)
    if (face.size == 4) return true;
  return false;
}

// Get obj vertices
inline void get_vertices(obj_shape* shape, int material,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<int>& vindex, bool flipv) {
  auto used_vertices = vector<bool>(shape->vertices.size(), false);
  auto count         = 0;
  for (auto& elem : shape->faces) {
    if (elem.material == material) {
      for (auto vid = count; vid < count + elem.size; vid++)
        used_vertices[vid] = true;
    }
    count += elem.size;
  }
  for (auto& elem : shape->lines) {
    if (elem.material == material) {
      for (auto vid = count; vid < count + elem.size; vid++)
        used_vertices[vid] = true;
    }
    count += elem.size;
  }
  auto vmap = unordered_map<obj_vertex, int>{};
  vmap.reserve(shape->vertices.size());
  vindex.resize(shape->vertices.size());
  for (auto vid = 0; vid < shape->vertices.size(); vid++) {
    if (!used_vertices[vid]) {
      vindex[vid] = -1;
      continue;
    }
    auto& vert = shape->vertices[vid];
    auto  it   = vmap.find(vert);
    if (it != vmap.end()) {
      vindex[vid] = it->second;
      continue;
    }
    auto nverts = (int)positions.size();
    vindex[vid] = nverts;
    vmap.insert(it, {vert, nverts});
    if (!shape->positions.empty() && vert.position)
      positions.push_back(shape->positions[vert.position - 1]);
    if (!shape->normals.empty() && vert.normal)
      normals.push_back(shape->normals[vert.normal - 1]);
    if (!shape->texcoords.empty() && vert.texcoord)
      texcoords.push_back(shape->texcoords[vert.texcoord - 1]);
  }
  if (flipv) {
    for (auto& texcoord : texcoords) texcoord.y = 1 - texcoord.y;
  }
}

// Get obj shape
inline void get_triangles(obj_model* obj, obj_shape* shape, int material,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, material, positions, normals, texcoords, vindex, flipv);
  triangles.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& elem : shape->faces) {
    if (elem.material == material) {
      for (auto c = 2; c < elem.size; c++) {
        triangles.push_back(
            {vindex[cur + 0], vindex[cur + c - 1], vindex[cur + c]});
      }
    }
    cur += elem.size;
  }
  triangles.shrink_to_fit();
}
inline void get_quads(obj_model* obj, obj_shape* shape, int material,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, material, positions, normals, texcoords, vindex, flipv);
  quads.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& elem : shape->faces) {
    if (elem.material == material) {
      if (elem.size == 4) {
        quads.push_back({vindex[cur + 0], vindex[cur + 1], vindex[cur + 2],
            vindex[cur + 3]});
      } else {
        for (auto c = 2; c < elem.size; c++) {
          quads.push_back({vindex[cur + 0], vindex[cur + c - 1],
              vindex[cur + c], vindex[cur + c]});
        }
      }
    }
    cur += elem.size;
  }
  quads.shrink_to_fit();
}
inline void get_lines(obj_model* obj, obj_shape* shape, int material,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flipv) {
  if (shape->lines.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, material, positions, normals, texcoords, vindex, flipv);
  lines.reserve(shape->lines.size());
  auto cur = 0;
  for (auto& elem : shape->lines) {
    if (elem.material == material) {
      for (auto c = 1; c < elem.size; c++) {
        lines.push_back({vindex[cur + c - 1], vindex[cur + c]});
      }
    }
    cur += elem.size;
  }
  lines.shrink_to_fit();
}
inline void get_points(obj_model* obj, obj_shape* shape, int material,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flipv) {
  if (shape->points.empty()) return;
  auto vindex = vector<int>{};
  get_vertices(shape, material, positions, normals, texcoords, vindex, flipv);
  points.reserve(shape->points.size());
  auto cur = 0;
  for (auto& elem : shape->points) {
    if (elem.material == material) {
      for (auto c = 0; c < elem.size; c++) {
        points.push_back({vindex[cur + 0]});
      }
    }
    cur += elem.size;
  }
}
inline vector<obj_material*> get_materials(obj_model* obj, obj_shape* shape) {
  return shape->materials;
}

// Add obj shape
inline void add_triangles(obj_model* obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials, const vector<int>& ematerials,
    const vector<frame3f>& instances, bool flipv) {
  auto shape       = obj->shapes.emplace_back(new obj_shape{});
  shape->name      = name;
  shape->materials = materials;
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape->instances = instances;
  shape->vertices.reserve(triangles.size() * 3);
  for (auto idx = 0; idx < triangles.size(); idx++) {
    auto& triangle = triangles[idx];
    for (auto c = 0; c < 3; c++) {
      shape->vertices.push_back({
          positions.empty() ? 0 : triangle[c] + 1,
          texcoords.empty() ? 0 : triangle[c] + 1,
          normals.empty() ? 0 : triangle[c] + 1,
      });
    }
    shape->faces.push_back(
        {3, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_quads(obj_model* obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials, const vector<int>& ematerials,
    const vector<frame3f>& instances, bool flipv) {
  auto shape       = obj->shapes.emplace_back(new obj_shape{});
  shape->name      = name;
  shape->materials = materials;
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape->instances = instances;
  shape->vertices.reserve(quads.size() * 4);
  for (auto idx = 0; idx < quads.size(); idx++) {
    auto& quad = quads[idx];
    auto  nv   = quad.z == quad.w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape->vertices.push_back({
          positions.empty() ? 0 : quad[c] + 1,
          texcoords.empty() ? 0 : quad[c] + 1,
          normals.empty() ? 0 : quad[c] + 1,
      });
    }
    shape->faces.push_back({(uint8_t)nv,
        ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_lines(obj_model* obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials, const vector<int>& ematerials,
    const vector<frame3f>& instances, bool flipv) {
  auto shape       = obj->shapes.emplace_back(new obj_shape{});
  shape->name      = name;
  shape->materials = materials;
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape->instances = instances;
  shape->vertices.reserve(lines.size() * 2);
  for (auto idx = 0; idx < lines.size(); idx++) {
    auto& str = lines[idx];
    for (auto c = 0; c < 2; c++) {
      shape->vertices.push_back({
          positions.empty() ? 0 : str[c] + 1,
          texcoords.empty() ? 0 : str[c] + 1,
          normals.empty() ? 0 : str[c] + 1,
      });
    }
    shape->lines.push_back(
        {2, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_points(obj_model* obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials, const vector<int>& ematerials,
    const vector<frame3f>& instances, bool flipv) {
  auto shape       = obj->shapes.emplace_back(new obj_shape{});
  shape->name      = name;
  shape->materials = materials;
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape->instances = instances;
  shape->vertices.reserve(points.size());
  for (auto idx = 0; idx < points.size(); idx++) {
    auto& point = points[idx];
    shape->vertices.push_back({
        positions.empty() ? 0 : point + 1,
        texcoords.empty() ? 0 : point + 1,
        normals.empty() ? 0 : point + 1,
    });
    shape->faces.push_back(
        {1, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_fvquads(obj_model* obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials, const vector<int>& ematerials,
    const vector<frame3f>& instances, bool flipv) {
  auto shape       = obj->shapes.emplace_back(new obj_shape{});
  shape->name      = name;
  shape->materials = materials;
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape->instances = instances;
  shape->vertices.reserve(quadspos.size() * 4);
  for (auto idx = 0; idx < quadspos.size(); idx++) {
    auto nv = quadspos[idx].z == quadspos[idx].w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape->vertices.push_back({
          quadspos.empty() ? 0 : quadspos[idx][c] + 1,
          quadstexcoord.empty() ? 0 : quadstexcoord[idx][c] + 1,
          quadsnorm.empty() ? 0 : quadsnorm[idx][c] + 1,
      });
    }
    shape->faces.push_back({(uint8_t)nv,
        ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}

}  // namespace yocto::obj

#endif
