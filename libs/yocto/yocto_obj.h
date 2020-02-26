//
// # Yocto/Obj: Tiny library for Obj parsing and writing
//
// Yocto/Obj is a tiny library for loading and saving Obj.
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::obj {

// Math defitions
using math::frame3f;
using math::identity3x4f;
using math::vec2f;
using math::vec2i;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;
using math::zero2f;
using math::zero3f;

}  // namespace yocto::obj

// -----------------------------------------------------------------------------
// OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::obj {

// OBJ vertex
struct vertex {
  int position = 0;
  int texcoord = 0;
  int normal   = 0;
};

inline bool operator==(const vertex& a, const vertex& b) {
  return a.position == b.position && a.texcoord == b.texcoord &&
         a.normal == b.normal;
}

// Obj texture information.
struct texture {
  std::string path  = "";     // file path
  bool        clamp = false;  // clamp to edge
  float       scale = 1;      // scale for bump/displacement

  // Properties not explicitly handled.
  std::unordered_map<std::string, std::vector<float>> props;

  texture() {}
  texture(const char* path) : path{path} {}
  texture(const std::string& path) : path{path} {}
};

// Obj element
struct element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj material
struct material {
  // material name and type
  std::string name  = "";
  int         illum = 0;

  // material colors and values
  vec3f emission     = {0, 0, 0};
  vec3f ambient      = {0, 0, 0};
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f reflection   = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  float exponent     = 10;
  float ior          = 1.5;
  float opacity      = 1;

  // material textures
  texture emission_tex     = {};
  texture ambient_tex      = {};
  texture diffuse_tex      = {};
  texture specular_tex     = {};
  texture reflection_tex   = {};
  texture transmission_tex = {};
  texture exponent_tex     = {};
  texture opacity_tex      = {};
  texture bump_tex         = {};
  texture normal_tex       = {};
  texture displacement_tex = {};

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
  vec3f pbr_volscattering = {0, 0, 0};
  float pbr_volanisotropy = 0;
  float pbr_volscale      = 0.01;

  // pbr extension textures
  texture pbr_emission_tex      = {};
  texture pbr_base_tex          = {};
  texture pbr_specular_tex      = {};
  texture pbr_roughness_tex     = {};
  texture pbr_metallic_tex      = {};
  texture pbr_sheen_tex         = {};
  texture pbr_coat_tex          = {};
  texture pbr_coatroughness_tex = {};
  texture pbr_transmission_tex  = {};
  texture pbr_opacity_tex       = {};
  texture pbr_volscattering_tex = {};
};

// Obj shape
struct shape {
  std::string            name      = "";
  std::vector<vec3f>     positions = {};
  std::vector<vec3f>     normals   = {};
  std::vector<vec2f>     texcoords = {};
  std::vector<material*> materials = {};
  std::vector<vertex>    vertices  = {};
  std::vector<element>   faces     = {};
  std::vector<element>   lines     = {};
  std::vector<element>   points    = {};
  std::vector<frame3f>   instances = {};
};

// Obj camera
struct camera {
  std::string name     = "";
  frame3f     frame    = identity3x4f;
  bool        ortho    = false;
  float       width    = 0.036;
  float       height   = 0.028;
  float       lens     = 0.50;
  float       focus    = 0;
  float       aperture = 0;
};

// Obj environment
struct environment {
  std::string name         = "";
  frame3f     frame        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  texture     emission_tex = {};
};

// Obj model
struct model {
  std::vector<std::string>       comments     = {};
  std::vector<obj::shape*>       shapes       = {};
  std::vector<obj::material*>    materials    = {};
  std::vector<obj::camera*>      cameras      = {};
  std::vector<obj::environment*> environments = {};
  ~model();
};

// Load and save obj
inline bool load_obj(const std::string& filename, obj::model* obj,
    std::string& error, bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
inline bool save_obj(
    const std::string& filename, obj::model* obj, std::string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
inline void get_triangles(const obj::shape* shape,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<obj::material*>& materials, std::vector<int>& ematerials,
    bool flip_texcoord = false);
inline void get_quads(const obj::shape* shape, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
inline void get_lines(const obj::shape* shape, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
inline void get_points(const obj::shape* shape, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
inline void get_fvquads(const obj::shape* shape, std::vector<vec4i>& quadspos,
    std::vector<vec4i>& quadsnorm, std::vector<vec4i>& quadstexcoord,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
inline bool has_quads(obj::shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
inline void get_triangles(const obj::shape* shape, int material,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);
inline void get_quads(const obj::shape* shape, int material,
    std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);
inline void get_lines(const obj::shape* shape, int material,
    std::vector<vec2i>& lines, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);
inline void get_points(const obj::shape* shape, int material,
    std::vector<int>& points, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);

// Create OBJ
inline obj::camera*      add_camera(obj::model* obj);
inline obj::material*    add_material(obj::model* obj);
inline obj::environment* add_environment(obj::model* obj);
inline obj::shape*       add_shape(obj::model* obj);

// Add obj shape
inline void set_triangles(obj::shape* shape,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
inline void add_quads(obj::shape* shape, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
inline void set_lines(obj::shape* shape, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>&          texcoords,
    const std::vector<obj::material*>& materials = {},
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
inline void set_points(obj::shape* shape, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
inline void set_fvquads(obj::shape* shape, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
inline void set_materials(
    obj::shape* shape, const std::vector<obj::material*>& materials);
inline void set_instances(
    obj::shape* shape, const std::vector<frame3f>& instances);

}  // namespace yocto::obj

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for std::vector for use with hash_map
template <>
struct hash<yocto::obj::vertex> {
  size_t operator()(const yocto::obj::vertex& v) const {
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
namespace sfs = ghc::filesystem;

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::obj {

// string literals
using namespace std::string_literals;

// utilities
inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_whitespace(std::string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

// Parse values from a std::string
[[nodiscard]] inline bool parse_value(
    std::string_view& str, std::string_view& value) {
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
[[nodiscard]] inline bool parse_value(
    std::string_view& str, std::string& value) {
  auto valuev = std::string_view{};
  if (!parse_value(str, valuev)) return false;
  value = std::string{valuev};
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, bool& value) {
  auto valuei = 0;
  if (!parse_value(str, valuei)) return false;
  value = (bool)valuei;
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}

[[nodiscard]] inline bool parse_value(std::string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

// Formats values to std::string
inline void format_value(std::string& str, const std::string& value) {
  str += value;
}
inline void format_value(std::string& str, int value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(std::string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_value(std::string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(std::string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(std::string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}

// Foramt to file
inline void format_values(std::string& str, const std::string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != std::string::npos)
    throw std::runtime_error("bad format std::string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_values(std::string& str, const std::string& fmt,
    const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == std::string::npos)
    throw std::invalid_argument("bad format std::string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
[[nodiscard]] inline bool format_values(
    FILE* fs, const std::string& fmt, const Args&... args) {
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

inline void remove_comment(std::string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

[[nodiscard]] inline bool parse_value(std::string_view& str, vertex& value) {
  value = vertex{0, 0, 0};
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
[[nodiscard]] inline bool parse_value(std::string_view& str, texture& info) {
  // initialize
  info = texture();

  // get tokens
  auto tokens = std::vector<std::string>();
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
  auto last = std::string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }

  return true;
}

// Read obj
[[nodiscard]] inline bool load_mtl(
    const std::string& filename, obj::model* obj, std::string& error) {
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
  add_material(obj);

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = std::string_view{buffer};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // grab material
    auto material = obj->materials.back();

    // possible token values
    if (cmd == "newmtl") {
      auto material = add_material(obj);
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
      if (max(material->transmission) < 0.001)
        material->transmission = {0, 0, 0};
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
    const std::string& filename, obj::model* obj, std::string& error) {
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
  auto shape_map = std::unordered_map<std::string, std::vector<obj::shape*>>{};
  for (auto shape : obj->shapes) {
    shape_map[shape->name].push_back(shape);
  }

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = std::string_view{buffer};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto camera = add_camera(obj);
      if (!parse_value(str, camera->name)) return parse_error();
      if (!parse_value(str, camera->ortho)) return parse_error();
      if (!parse_value(str, camera->width)) return parse_error();
      if (!parse_value(str, camera->height)) return parse_error();
      if (!parse_value(str, camera->lens)) return parse_error();
      if (!parse_value(str, camera->focus)) return parse_error();
      if (!parse_value(str, camera->aperture)) return parse_error();
      if (!parse_value(str, camera->frame)) return parse_error();
    } else if (cmd == "e") {
      auto environment = add_environment(obj);
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

inline model::~model() {
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto camera : cameras) delete camera;
  for (auto environment : environments) delete environment;
}

// Make obj
inline obj::camera* add_camera(obj::model* obj) {
  return obj->cameras.emplace_back(new camera{});
}
inline obj::material* add_material(obj::model* obj) {
  return obj->materials.emplace_back(new material{});
}
inline obj::environment* add_environment(obj::model* obj) {
  return obj->environments.emplace_back(new environment{});
}
inline obj::shape* add_shape(obj::model* obj) {
  return obj->shapes.emplace_back(new shape{});
}

// Read obj
inline bool load_obj(const std::string& filename, obj::model* obj,
    std::string& error, bool geom_only, bool split_elements,
    bool split_materials) {
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
  auto opositions   = std::vector<vec3f>{};
  auto onormals     = std::vector<vec3f>{};
  auto otexcoords   = std::vector<vec2f>{};
  auto vert_size    = vertex{};
  auto oname        = ""s;
  auto gname        = ""s;
  auto mname        = ""s;
  auto mtllibs      = std::vector<std::string>{};
  auto material_map = std::unordered_map<std::string, obj::material*>{};

  // initialize obj
  obj->~model();
  obj->cameras.clear();
  obj->environments.clear();
  obj->shapes.clear();
  obj->materials.clear();

  // initialize load
  obj->shapes.emplace_back(new shape{});
  auto empty_material = (obj::material*)nullptr;

  // read the file str by str
  char buffer[4096];
  while (fgets(buffer, sizeof(buffer), fs)) {
    // str
    auto str = std::string_view{buffer};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      if (!parse_value(str, opositions.emplace_back())) return parse_error();
      vert_size.position += 1;
    } else if (cmd == "vn") {
      if (!parse_value(str, onormals.emplace_back())) return parse_error();
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      if (!parse_value(str, otexcoords.emplace_back())) return parse_error();
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      // split if split_elements and different primitives
      if (auto shape = obj->shapes.back();
          split_elements && !shape->vertices.empty()) {
        if ((cmd == "f" && (!shape->lines.empty() || !shape->points.empty())) ||
            (cmd == "l" && (!shape->faces.empty() || !shape->points.empty())) ||
            (cmd == "p" && (!shape->faces.empty() || !shape->lines.empty()))) {
          add_shape(obj);
          obj->shapes.back()->name = oname + gname;
        }
      }
      // split if splt_material and different materials
      if (auto shape = obj->shapes.back();
          !geom_only && split_materials && !shape->materials.empty()) {
        if (shape->materials.size() > 1)
          throw std::runtime_error("should not have happened");
        if (shape->materials.back()->name != mname) {
          add_shape(obj);
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
          empty_material   = obj->materials.emplace_back(new material{});
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
        auto vert = vertex{};
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
        obj->shapes.emplace_back(new shape{});
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
        if (!load_mtl(sfs::path(filename).parent_path() / mtllib, obj, error))
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
  auto ipositions = std::vector<int>{};
  auto inormals   = std::vector<int>{};
  auto itexcoords = std::vector<int>{};
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
  auto extfilename = sfs::path(filename).replace_extension(".objx");
  if (sfs::exists(sfs::path(extfilename))) {
    if (!load_objx(extfilename, obj, error)) return dependent_error();
  }

  return true;
}

// Format values
inline void format_value(std::string& str, const texture& value) {
  str += value.path.empty() ? "" : value.path;
}
inline void format_value(std::string& str, const vertex& value) {
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
    const std::string& filename, obj::model* obj, std::string& error) {
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
    const std::string& filename, obj::model* obj, std::string& error) {
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
    const std::string& filename, obj::model* obj, std::string& error) {
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
            sfs::path(filename).filename().replace_extension(".mtl")))
      return write_error();
  }

  // save objects
  auto vert_size = vertex{0, 0, 0};
  for (auto shape : obj->shapes) {
    if (!format_values(fs, "o {}\n", shape->name)) return write_error();
    for (auto& p : shape->positions)
      if (!format_values(fs, "v {}\n", p)) return write_error();
    for (auto& n : shape->normals)
      if (!format_values(fs, "vn {}\n", n)) return write_error();
    for (auto& t : shape->texcoords)
      if (!format_values(fs, "vt {}\n", t)) return write_error();
    auto element_labels = std::vector<std::string>{"f", "l", "p"};
    auto element_groups = std::vector<const std::vector<element>*>{
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
    if (!save_mtl(sfs::path(filename).replace_extension(".mtl"), obj, error))
      return dependent_error();
  }

  // save objx
  if (!obj->cameras.empty() || !obj->environments.empty() ||
      std::any_of(obj->shapes.begin(), obj->shapes.end(),
          [](auto shape) { return !shape->instances.empty(); })) {
    if (!save_objx(sfs::path(filename).replace_extension(".objx"), obj, error))
      return dependent_error();
  }

  // done
  return true;
}

// Get obj vertices
inline void get_vertices(const obj::shape* shape, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<int>& vindex, bool flipv) {
  auto vmap = std::unordered_map<vertex, int>{};
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
inline std::vector<vec2f> flip_texcoord(const std::vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get obj shape
inline void get_triangles(const obj::shape* shape,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    std::vector<obj::material*>& materials, std::vector<int>& ematerials,
    bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_quads(const obj::shape* shape, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_lines(const obj::shape* shape, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flipv) {
  if (shape->lines.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_points(const obj::shape* shape, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flipv) {
  if (shape->points.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_fvquads(const obj::shape* shape, std::vector<vec4i>& quadspos,
    std::vector<vec4i>& quadsnorm, std::vector<vec4i>& quadstexcoord,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  positions = shape->positions;
  normals   = shape->normals;
  texcoords = flipv ? flip_texcoord(shape->texcoords) : shape->texcoords;
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

inline bool has_quads(obj::shape* shape) {
  for (auto& face : shape->faces)
    if (face.size == 4) return true;
  return false;
}

// Get obj vertices
inline void get_vertices(const obj::shape* shape, int material,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<int>& vindex, bool flipv) {
  auto used_vertices = std::vector<bool>(shape->vertices.size(), false);
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
  auto vmap = std::unordered_map<vertex, int>{};
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
inline void get_triangles(const obj::shape* shape, int material,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_quads(const obj::shape* shape, int material,
    std::vector<vec4i>& quads, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, bool flipv) {
  if (shape->faces.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_lines(const obj::shape* shape, int material,
    std::vector<vec2i>& lines, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, bool flipv) {
  if (shape->lines.empty()) return;
  auto vindex = std::vector<int>{};
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
inline void get_points(const obj::shape* shape, int material,
    std::vector<int>& points, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords, bool flipv) {
  if (shape->points.empty()) return;
  auto vindex = std::vector<int>{};
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

// Add obj shape
inline void set_triangles(obj::shape* shape,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& positions,
    const std::vector<vec3f>& normals, const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
inline void set_quads(obj::shape* shape, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<int>& ematerials,
    bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
inline void set_lines(obj::shape* shape, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<int>& ematerials,
    bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
inline void set_points(obj::shape* shape, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<int>& ematerials,
    bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
inline void set_fvquads(obj::shape* shape, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<int>& ematerials,
    bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
inline void set_materials(
    obj::shape* shape, const std::vector<obj::material*>& materials) {
  shape->materials = materials;
}
inline void set_instances(
    obj::shape* shape, const std::vector<frame3f>& instances) {
  shape->instances = instances;
}

}  // namespace yocto::obj

#endif
