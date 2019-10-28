//
// # Yocto/Obj: Tiny library for Obj parsing and writing
//
// Yocto/Obj is a tiny library for loading and saving Wavefront Obj. Yocto/Obj
// supports two interfaces: a simple interface where all Obj data is loaded
// and saved at once and a low-level interface where Obj values are read
// and written one at a time.
// Error reporting is done by throwing `std::runtime_error` exceptions.
//
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

#include "yocto_math.h"

#include <algorithm>

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

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

// Obj shape
struct obj_shape {
  string              name      = "";
  vector<vec3f>       positions = {};
  vector<vec3f>       normals   = {};
  vector<vec2f>       texcoords = {};
  vector<string>      materials = {};
  vector<obj_vertex>  vertices  = {};
  vector<obj_element> faces     = {};
  vector<obj_element> lines     = {};
  vector<obj_element> points    = {};
  vector<frame3f>     instances = {};
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
  obj_texture_info emission_map     = {};
  obj_texture_info ambient_map      = {};
  obj_texture_info diffuse_map      = {};
  obj_texture_info specular_map     = {};
  obj_texture_info reflection_map   = {};
  obj_texture_info transmission_map = {};
  obj_texture_info exponent_map     = {};
  obj_texture_info opacity_map      = {};
  obj_texture_info bump_map         = {};
  obj_texture_info normal_map       = {};
  obj_texture_info displacement_map = {};

  // pbrt extension values
  float pbr_roughness     = 0;
  float pbr_metallic      = 0;
  float pbr_sheen         = 0;
  float pbr_clearcoat     = 0;
  float pbr_coatroughness = 0;

  // pbr extension textures
  obj_texture_info pbr_roughness_map     = {};
  obj_texture_info pbr_metallic_map      = {};
  obj_texture_info pbr_sheen_map         = {};
  obj_texture_info pbr_clearcoat_map     = {};
  obj_texture_info pbr_coatroughness_map = {};

  // volume extension colors and values
  vec3f vol_emission     = zero3f;
  vec3f vol_transmission = zero3f;
  vec3f vol_meanfreepath = zero3f;
  vec3f vol_scattering   = zero3f;
  float vol_anisotropy   = 0;
  float vol_scale        = 0.01;

  // volument textures
  obj_texture_info vol_scattering_map = {};
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
  obj_texture_info emission_map = {};
};

// Obj model
struct obj_model {
  vector<string>          comments     = {};
  vector<obj_shape>       shapes       = {};
  vector<obj_material>    materials    = {};
  vector<obj_camera>      cameras      = {};
  vector<obj_environment> environments = {};
};

// Load and save obj
inline void load_obj(const string& filename, obj_model& obj,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
inline void save_obj(const string& filename, const obj_model& obj);

// convert between roughness and exponent
inline float obj_exponent_to_roughness(float exponent);
inline float obj_roughness_to_exponent(float roughness);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
inline void get_obj_triangles(const obj_model& obj, const obj_shape& shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline void get_obj_quads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline void get_obj_lines(const obj_model& obj, const obj_shape& shape,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline void get_obj_points(const obj_model& obj, const obj_shape& shape,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline void get_obj_fvquads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
inline bool has_obj_quads(const obj_shape& shape);

// Add obj shape
inline void add_obj_triangles(obj_model& obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
inline void add_obj_quads(obj_model& obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
inline void add_obj_lines(obj_model& obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
inline void add_obj_points(obj_model& obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
inline void add_obj_fvquads(obj_model& obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------
namespace yocto {

// A class that wraps a C file ti handle safe opening/closgin with RIIA.
struct obj_file {
  obj_file() {}
  obj_file(obj_file&& other);
  obj_file(const obj_file&) = delete;
  obj_file& operator=(const obj_file&) = delete;
  ~obj_file();

  FILE*  fs       = nullptr;
  string filename = "";
  string mode     = "rt";
  int    linenum  = 0;
};

// open a file
inline obj_file open_obj(const string& filename, const string& mode = "rt");
inline void     open_obj(
        obj_file& fs, const string& filename, const string& mode = "rt");
inline void close_obj(obj_file& fs);

// Obj/Mtl/Objx command
enum struct obj_command {
  // clang-format off
  vertex, normal, texcoord,         // data in value
  face, line, point,                // data in vertices
  object, group, usemtl, smoothing, // data in name
  mtllib, objxlib,                  // data in name
  // clang-format on
};
enum struct mtl_command { material };
enum struct objx_command { camera, environment, instance };

// Obj instance
struct obj_instance {
  string  object = "";
  frame3f frame  = identity3x4f;
};

// Read obj/mtl/objx elements
inline bool read_obj_command(obj_file& fs, obj_command& command, string& name,
    vec3f& value, vector<obj_vertex>& vertices, obj_vertex& vert_size);
inline bool read_mtl_command(obj_file& fs, mtl_command& command,
    obj_material& material, bool fliptr = true);
inline bool read_objx_command(obj_file& fs, objx_command& command,
    obj_camera& camera, obj_environment& environment, obj_instance& instance);

// Write obj/mtl/objx elements
inline void write_obj_comment(obj_file& fs, const string& comment);
inline void write_obj_command(obj_file& fs, obj_command command,
    const string& name, const vec3f& value,
    const vector<obj_vertex>& vertices = {});
inline void write_mtl_command(obj_file& fs, mtl_command command,
    obj_material& material, const obj_texture_info& texture = {});
inline void write_objx_command(obj_file& fs, objx_command command,
    const obj_camera& camera, const obj_environment& environment,
    const obj_instance& instance);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::obj_vertex> {
  size_t operator()(const yocto::obj_vertex& v) const {
    static const std::hash<int> hasher = std::hash<int>();
    auto                        h      = (size_t)0;
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

#include <string_view>

// -----------------------------------------------------------------------------
// LOW-LEVEL FILE HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
inline obj_file::obj_file(obj_file&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
inline obj_file::~obj_file() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
inline void open_obj(obj_file& fs, const string& filename, const string& mode) {
  close_obj(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
inline obj_file open_obj(const string& filename, const string& mode) {
  auto fs = obj_file{};
  open_obj(fs, filename, mode);
  return fs;
}
inline void close_obj(obj_file& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

inline bool read_obj_line(obj_file& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

inline void write_ply_text(obj_file& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
inline void write_ply_text(obj_file& fs, const char* value) {
  if (fputs(value, fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOAD-LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;

// utilities
inline bool is_obj_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_obj_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_obj_whitespace(string_view& str) {
  while (!str.empty() && is_obj_space(str.front())) str.remove_prefix(1);
}

inline void remove_obj_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_obj_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
inline void parse_obj_value(string_view& str, string_view& value) {
  skip_obj_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  auto cpy = str;
  while (!cpy.empty() && !is_obj_space(cpy.front())) cpy.remove_prefix(1);
  value = str;
  value.remove_suffix(cpy.size());
  str.remove_prefix(str.size() - cpy.size());
}
inline void parse_obj_value(string_view& str, string& value) {
  auto valuev = string_view{};
  parse_obj_value(str, valuev);
  value = string{valuev};
}
inline void parse_obj_value(string_view& str, int& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_obj_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_obj_value(str, valuei);
  value = (bool)valuei;
}
inline void parse_obj_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
inline void parse_obj_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_obj_value(str, values[i]);
}

inline void parse_obj_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++) parse_obj_value(str, value[i]);
}
inline void parse_obj_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++) parse_obj_value(str, value[i]);
}
inline void parse_obj_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++) parse_obj_value(str, value[i]);
}

// Parse values from a string
inline void parse_obj_value_or_empty(string_view& str, string& value) {
  skip_obj_whitespace(str);
  if (str.empty()) {
    value = "";
  } else {
    parse_obj_value(str, value);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL PRINTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_obj_value(string& str, const string& value) { str += value; }
inline void format_obj_value(string& str, const char* value) { str += value; }
inline void format_obj_value(string& str, int value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_obj_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_obj_value(string& str, const vec2f& value) {
  for(auto i = 0; i < 2; i ++) format_obj_value(str, value[i]);
}
inline void format_obj_value(string& str, const vec3f& value) {
  for(auto i = 0; i < 3; i ++) format_obj_value(str, value[i]);
}
inline void format_obj_value(string& str, const frame3f& value) {
  for(auto i = 0; i < 4; i ++) format_obj_value(str, value[i]);
}

// Foramt to file
inline void format_obj_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_obj_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_obj_value(str, arg);
  format_obj_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
inline void format_obj_values(
    obj_file& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_obj_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
inline void format_obj_value(obj_file& fs, const T& value) {
  auto str = ""s;
  format_obj_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

inline void parse_obj_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  parse_obj_value(str, value.position);
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      parse_obj_value(str, value.normal);
    } else {
      parse_obj_value(str, value.texcoord);
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        parse_obj_value(str, value.normal);
      }
    }
  }
}

// Input for OBJ textures
inline void parse_obj_value(string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_obj_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    parse_obj_value(str, token);
    tokens.push_back(token);
    skip_obj_whitespace(str);
  }
  if (tokens.empty()) throw std::runtime_error("cannot parse value");

  // texture name
  info.path = normalize_path(tokens.back());

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }
}

// Read obj
inline void load_mtl(
    const string& filename, obj_model& obj, bool fliptr = true) {
  // open file
  auto fs = open_obj(filename, "rt");

  // init parsing
  obj.materials.emplace_back();

  // read the file line by line
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      obj.materials.emplace_back();
      parse_obj_value(line, obj.materials.back().name);
    } else if (cmd == "illum") {
      parse_obj_value(line, obj.materials.back().illum);
    } else if (cmd == "Ke") {
      parse_obj_value(line, obj.materials.back().emission);
    } else if (cmd == "Ka") {
      parse_obj_value(line, obj.materials.back().ambient);
    } else if (cmd == "Kd") {
      parse_obj_value(line, obj.materials.back().diffuse);
    } else if (cmd == "Ks") {
      parse_obj_value(line, obj.materials.back().specular);
    } else if (cmd == "Kt") {
      parse_obj_value(line, obj.materials.back().transmission);
    } else if (cmd == "Tf") {
      obj.materials.back().transmission = vec3f{-1};
      parse_obj_value(line, obj.materials.back().transmission);
      if (obj.materials.back().transmission.y < 0)
        obj.materials.back().transmission = vec3f{
            obj.materials.back().transmission.x};
      if (fliptr)
        obj.materials.back().transmission = 1 -
                                            obj.materials.back().transmission;
    } else if (cmd == "Tr") {
      parse_obj_value(line, obj.materials.back().opacity);
      if (fliptr)
        obj.materials.back().opacity = 1 - obj.materials.back().opacity;
    } else if (cmd == "Ns") {
      parse_obj_value(line, obj.materials.back().exponent);
    } else if (cmd == "d") {
      parse_obj_value(line, obj.materials.back().opacity);
    } else if (cmd == "map_Ke") {
      parse_obj_value(line, obj.materials.back().emission_map);
    } else if (cmd == "map_Ka") {
      parse_obj_value(line, obj.materials.back().ambient_map);
    } else if (cmd == "map_Kd") {
      parse_obj_value(line, obj.materials.back().diffuse_map);
    } else if (cmd == "map_Ks") {
      parse_obj_value(line, obj.materials.back().specular_map);
    } else if (cmd == "map_Tr") {
      parse_obj_value(line, obj.materials.back().transmission_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_obj_value(line, obj.materials.back().opacity_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_obj_value(line, obj.materials.back().bump_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_obj_value(line, obj.materials.back().displacement_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_obj_value(line, obj.materials.back().normal_map);
    } else if (cmd == "Pm") {
      parse_obj_value(line, obj.materials.back().pbr_metallic);
    } else if (cmd == "Pr") {
      parse_obj_value(line, obj.materials.back().pbr_roughness);
    } else if (cmd == "Ps") {
      parse_obj_value(line, obj.materials.back().pbr_sheen);
    } else if (cmd == "Pc") {
      parse_obj_value(line, obj.materials.back().pbr_clearcoat);
    } else if (cmd == "Pcr") {
      parse_obj_value(line, obj.materials.back().pbr_coatroughness);
    } else if (cmd == "map_Pm") {
      parse_obj_value(line, obj.materials.back().pbr_metallic_map);
    } else if (cmd == "map_Pr") {
      parse_obj_value(line, obj.materials.back().pbr_roughness_map);
    } else if (cmd == "map_Ps") {
      parse_obj_value(line, obj.materials.back().pbr_sheen_map);
    } else if (cmd == "map_Pc") {
      parse_obj_value(line, obj.materials.back().pbr_clearcoat_map);
    } else if (cmd == "map_Pcr") {
      parse_obj_value(line, obj.materials.back().pbr_coatroughness_map);
    } else if (cmd == "Vt") {
      parse_obj_value(line, obj.materials.back().vol_transmission);
    } else if (cmd == "Vp") {
      parse_obj_value(line, obj.materials.back().vol_meanfreepath);
    } else if (cmd == "Ve") {
      parse_obj_value(line, obj.materials.back().vol_emission);
    } else if (cmd == "Vs") {
      parse_obj_value(line, obj.materials.back().vol_scattering);
    } else if (cmd == "Vg") {
      parse_obj_value(line, obj.materials.back().vol_anisotropy);
    } else if (cmd == "Vr") {
      parse_obj_value(line, obj.materials.back().vol_scale);
    } else if (cmd == "map_Vs") {
      parse_obj_value(line, obj.materials.back().vol_scattering_map);
    } else {
      continue;
    }
  }

  // remove placeholder material
  obj.materials.erase(obj.materials.begin());
}

// Read obj
inline void load_objx(const string& filename, obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "rt");

  // shape map for instances
  auto shape_map = unordered_map<string, vector<int>>{};
  for (auto idx = 0; idx < obj.shapes.size(); idx++) {
    shape_map[obj.shapes[idx].name].push_back(idx);
  }

  // read the file line by line
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto& camera = obj.cameras.emplace_back();
      parse_obj_value(line, camera.name);
      parse_obj_value(line, camera.ortho);
      parse_obj_value(line, camera.width);
      parse_obj_value(line, camera.height);
      parse_obj_value(line, camera.lens);
      parse_obj_value(line, camera.focus);
      parse_obj_value(line, camera.aperture);
      parse_obj_value(line, camera.frame);
    } else if (cmd == "e") {
      auto& environment = obj.environments.emplace_back();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.emission);
      auto emission_path = ""s;
      parse_obj_value(line, emission_path);
      if (emission_path == "\"\"") emission_path = "";
      environment.emission_map.path = emission_path;
      parse_obj_value(line, environment.frame);
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      parse_obj_value(line, object);
      parse_obj_value(line, frame);
      if (shape_map.find(object) == shape_map.end()) {
        throw std::runtime_error("cannot find object " + object);
      }
      for (auto idx : shape_map.at(object)) {
        obj.shapes[idx].instances.push_back(frame);
      }
    } else {
      // unused
    }
  }
}

// Read obj
inline void load_obj(const string& filename, obj_model& obj, bool geom_only,
    bool split_elements, bool split_materials) {
  // open file
  auto fs = open_obj(filename, "rt");

  // parsing state
  auto opositions = vector<vec3f>{};
  auto onormals   = vector<vec3f>{};
  auto otexcoords = vector<vec2f>{};
  auto vert_size  = obj_vertex{};
  auto oname      = ""s;
  auto gname      = ""s;
  auto mname      = ""s;
  auto mtllibs    = vector<string>{};

  // initialize obj
  obj = {};
  obj.shapes.emplace_back();

  // read the file line by line
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      parse_obj_value(line, opositions.emplace_back());
      vert_size.position += 1;
    } else if (cmd == "vn") {
      parse_obj_value(line, onormals.emplace_back());
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      parse_obj_value(line, otexcoords.emplace_back());
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      // split if split_elements and different primitives
      if (auto& shape = obj.shapes.back();
          split_elements && !shape.vertices.empty()) {
        if ((cmd == "f" && (!shape.lines.empty() || !shape.points.empty())) ||
            (cmd == "l" && (!shape.faces.empty() || !shape.points.empty())) ||
            (cmd == "p" && (!shape.faces.empty() || !shape.lines.empty()))) {
          obj.shapes.emplace_back();
          obj.shapes.back().name = oname + gname;
        }
      }
      // split if splt_material and different materials
      if (auto& shape = obj.shapes.back();
          !geom_only && split_materials && !shape.materials.empty()) {
        if (shape.materials.size() > 1)
          throw std::runtime_error("should not have happened");
        if (shape.materials.back() != mname) {
          obj.shapes.emplace_back();
          obj.shapes.back().name = oname + gname;
        }
      }
      // grab shape and add element
      auto& shape   = obj.shapes.back();
      auto& element = (cmd == "f") ? shape.faces.emplace_back()
                                   : (cmd == "l") ? shape.lines.emplace_back()
                                                  : shape.points.emplace_back();
      // get element material or add if needed
      if (!geom_only) {
        auto mat_idx = -1;
        for (auto midx = 0; midx < shape.materials.size(); midx++)
          if (shape.materials[midx] == mname) mat_idx = midx;
        if (mat_idx < 0) {
          shape.materials.push_back(mname);
          mat_idx = shape.materials.size() - 1;
        }
        element.material = (uint8_t)mat_idx;
      }
      // parse vertices
      skip_obj_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_obj_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        shape.vertices.push_back(vert);
        element.size += 1;
        skip_obj_whitespace(line);
      }
    } else if (cmd == "o" || cmd == "g") {
      if (geom_only) continue;
      parse_obj_value_or_empty(line, cmd == "o" ? oname : gname);
      if (!obj.shapes.back().vertices.empty()) {
        obj.shapes.emplace_back();
        obj.shapes.back().name = oname + gname;
      } else {
        obj.shapes.back().name = oname + gname;
      }
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      parse_obj_value_or_empty(line, mname);
    } else if (cmd == "s") {
      if (geom_only) continue;
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      parse_obj_value(line, mtllib);
      if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) == mtllibs.end()) {
        mtllibs.push_back(mtllib);
      }
    } else {
      // unused
    }
  }

  // convert vertex data
  auto ipositions = vector<int>{};
  auto inormals   = vector<int>{};
  auto itexcoords = vector<int>{};
  for (auto& shape : obj.shapes) {
    ipositions.assign(opositions.size() + 1, 0);
    inormals.assign(onormals.size() + 1, 0);
    itexcoords.assign(otexcoords.size() + 1, 0);
    for (auto& vertex : shape.vertices) {
      if (vertex.position && !ipositions[vertex.position]) {
        shape.positions.push_back(opositions[vertex.position - 1]);
        ipositions[vertex.position] = (int)shape.positions.size();
      }
      if (vertex.normal && !inormals[vertex.normal]) {
        shape.normals.push_back(onormals[vertex.normal - 1]);
        inormals[vertex.normal] = (int)shape.normals.size();
      }
      if (vertex.texcoord && !itexcoords[vertex.texcoord]) {
        shape.texcoords.push_back(otexcoords[vertex.texcoord - 1]);
        itexcoords[vertex.texcoord] = (int)shape.texcoords.size();
      }
      vertex.position = ipositions[vertex.position];
      vertex.normal   = inormals[vertex.normal];
      vertex.texcoord = itexcoords[vertex.texcoord];
    }
  }

  // exit if done
  if (geom_only) return;

  // load materials
  auto dirname = get_dirname(filename);
  for (auto& mtllib : mtllibs) {
    load_mtl(dirname + mtllib, obj);
  }

  // load extensions
  auto extfilename = replace_extension(filename, ".objx");
  if (exists_file(extfilename)) {
    load_objx(extfilename, obj);
  }
}

// Format values
inline void format_obj_value(string& str, const obj_texture_info& value) {
  str += value.path.empty() ? "" : value.path;
}
inline void format_obj_value(string& str, const obj_vertex& value) {
  format_obj_value(str, value.position);
  if (value.texcoord) {
    str += "/";
    format_obj_value(str, value.texcoord);
    if (value.normal) {
      str += "/";
      format_obj_value(str, value.normal);
    }
  } else if (value.normal) {
    str += "//";
    format_obj_value(str, value.normal);
  }
}

// Save obj
inline void save_mtl(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");

  // save comments
  format_obj_values(fs, "#\n");
  format_obj_values(fs, "# Written by Yocto/GL\n");
  format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_obj_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_obj_values(fs, "# {}\n", comment);
  }
  format_obj_values(fs, "\n");

  // write material
  for (auto& material : obj.materials) {
    format_obj_values(fs, "newmtl {}\n", material.name);
    format_obj_values(fs, "illum {}\n", material.illum);
    if (material.emission != zero3f)
      format_obj_values(fs, "Ke {}\n", material.emission);
    if (material.ambient != zero3f)
      format_obj_values(fs, "Ka {}\n", material.ambient);
    format_obj_values(fs, "Kd {}\n", material.diffuse);
    format_obj_values(fs, "Ks {}\n", material.specular);
    if (material.reflection != zero3f)
      format_obj_values(fs, "Kr {}\n", material.reflection);
    if (material.transmission != zero3f)
      format_obj_values(fs, "Kt {}\n", material.transmission);
    format_obj_values(fs, "Ns {}\n", (int)material.exponent);
    if (material.opacity != 1)
      format_obj_values(fs, "d {}\n", material.opacity);
    if (!material.emission_map.path.empty())
      format_obj_values(fs, "map_Ke {}\n", material.emission_map);
    if (!material.diffuse_map.path.empty())
      format_obj_values(fs, "map_Kd {}\n", material.diffuse_map);
    if (!material.specular_map.path.empty())
      format_obj_values(fs, "map_Ks {}\n", material.specular_map);
    if (!material.transmission_map.path.empty())
      format_obj_values(fs, "map_Kt {}\n", material.transmission_map);
    if (!material.reflection_map.path.empty())
      format_obj_values(fs, "map_Kr {}\n", material.reflection_map);
    if (!material.exponent_map.path.empty())
      format_obj_values(fs, "map_Ns {}\n", material.exponent_map);
    if (!material.opacity_map.path.empty())
      format_obj_values(fs, "map_d {}\n", material.opacity_map);
    if (!material.bump_map.path.empty())
      format_obj_values(fs, "map_bump {}\n", material.bump_map);
    if (!material.displacement_map.path.empty())
      format_obj_values(fs, "map_disp {}\n", material.displacement_map);
    if (!material.normal_map.path.empty())
      format_obj_values(fs, "map_norm {}\n", material.normal_map);
    if (material.pbr_roughness)
      format_obj_values(fs, "Pr {}\n", material.pbr_roughness);
    if (material.pbr_metallic)
      format_obj_values(fs, "Pm {}\n", material.pbr_metallic);
    if (material.pbr_sheen)
      format_obj_values(fs, "Ps {}\n", material.pbr_sheen);
    if (material.pbr_clearcoat)
      format_obj_values(fs, "Pc {}\n", material.pbr_clearcoat);
    if (material.pbr_coatroughness)
      format_obj_values(fs, "Pcr {}\n", material.pbr_coatroughness);
    if (!material.pbr_roughness_map.path.empty())
      format_obj_values(fs, "map_Pr {}\n", material.pbr_roughness_map);
    if (!material.pbr_metallic_map.path.empty())
      format_obj_values(fs, "map_Pm {}\n", material.pbr_metallic_map);
    if (!material.pbr_sheen_map.path.empty())
      format_obj_values(fs, "map_Ps {}\n", material.pbr_sheen_map);
    if (!material.pbr_clearcoat_map.path.empty())
      format_obj_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map);
    if (!material.pbr_coatroughness_map.path.empty())
      format_obj_values(fs, "map_Pcr {}\n", material.pbr_coatroughness_map);
    if (material.vol_transmission != zero3f)
      format_obj_values(fs, "Vt {}\n", material.vol_transmission);
    if (material.vol_meanfreepath != zero3f)
      format_obj_values(fs, "Vp {}\n", material.vol_meanfreepath);
    if (material.vol_emission != zero3f)
      format_obj_values(fs, "Ve {}\n", material.vol_emission);
    if (material.vol_scattering != zero3f)
      format_obj_values(fs, "Vs {}\n", material.vol_scattering);
    if (material.vol_anisotropy)
      format_obj_values(fs, "Vg {}\n", material.vol_anisotropy);
    if (material.vol_scale)
      format_obj_values(fs, "Vr {}\n", material.vol_scale);
    if (!material.vol_scattering_map.path.empty())
      format_obj_values(fs, "map_Vs {}\n", material.vol_scattering_map);
    format_obj_values(fs, "\n");
  }
}

// Save obj
inline void save_objx(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");

  // save comments
  format_obj_values(fs, "#\n");
  format_obj_values(fs, "# Written by Yocto/GL\n");
  format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_obj_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_obj_values(fs, "# {}\n", comment);
  }
  format_obj_values(fs, "\n");

  // cameras
  for (auto& camera : obj.cameras) {
    format_obj_values(fs, "c {} {} {} {} {} {} {} {}\n", camera.name,
        camera.ortho, camera.width, camera.height, camera.lens, camera.focus,
        camera.aperture, camera.frame);
  }

  // environments
  for (auto& environment : obj.environments) {
    format_obj_values(fs, "e {} {} {} {}\n", environment.name,
        environment.emission,
        environment.emission_map.path.empty() ? "\"\""s
                                              : environment.emission_map.path,
        environment.frame);
  }

  // instances
  for (auto& shape : obj.shapes) {
    for (auto& frame : shape.instances) {
      format_obj_values(fs, "i {} {}\n", shape.name, frame);
    }
  }
}

// Save obj
inline void save_obj(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");

  // save comments
  format_obj_values(fs, "#\n");
  format_obj_values(fs, "# Written by Yocto/GL\n");
  format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_obj_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_obj_values(fs, "# {}\n", comment);
  }
  format_obj_values(fs, "\n");

  // save material library
  if (!obj.materials.empty()) {
    format_obj_values(
        fs, "mtllib {}\n\n", replace_extension(get_filename(filename), ".mtl"));
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto& shape : obj.shapes) {
    format_obj_values(fs, "o {}\n", shape.name);
    for (auto& p : shape.positions) format_obj_values(fs, "v {}\n", p);
    for (auto& n : shape.normals) format_obj_values(fs, "vn {}\n", n);
    for (auto& t : shape.texcoords) format_obj_values(fs, "vt {}\n", t);
    auto element_labels = vector<string>{"f", "l", "p"};
    auto element_groups = vector<const vector<obj_element>*>{
        &shape.faces, &shape.lines, &shape.points};
    for (auto element_idx = 0; element_idx < 3; element_idx++) {
      auto& label        = element_labels[element_idx];
      auto& elements     = *element_groups[element_idx];
      auto  cur_material = -1, cur_vertex = 0;
      for (auto& element : elements) {
        if (!shape.materials.empty() && cur_material != element.material) {
          format_obj_values(
              fs, "usemtl {}\n", shape.materials[element.material]);
          cur_material = element.material;
        }
        format_obj_values(fs, "{}", label);
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape.vertices[cur_vertex++];
          if (vert.position) vert.position += vert_size.position;
          if (vert.normal) vert.normal += vert_size.normal;
          if (vert.texcoord) vert.texcoord += vert_size.texcoord;
          format_obj_values(fs, " {}", vert);
        }
        format_obj_values(fs, "\n");
      }
    }
    format_obj_values(fs, "\n");
    vert_size.position += (int)shape.positions.size();
    vert_size.normal += (int)shape.normals.size();
    vert_size.texcoord += (int)shape.texcoords.size();
  }

  // save mtl
  if (!obj.materials.empty())
    save_mtl(replace_extension(filename, ".mtl"), obj);

  // save objx
  if (!obj.cameras.empty() || !obj.environments.empty() ||
      std::any_of(obj.shapes.begin(), obj.shapes.end(),
          [](auto& shape) { return !shape.instances.empty(); }))
    save_objx(replace_extension(filename, ".objx"), obj);
}

// convert between roughness and exponent
inline float obj_exponent_to_roughness(float exponent) {
  auto roughness = exponent;
  roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
  if (roughness < 0.01f) roughness = 0;
  if (roughness > 0.99f) roughness = 1;
  return roughness;
}
inline float obj_roughness_to_exponent(float roughness) {
  return (int)clamp(
      2 / pow(clamp(roughness, 0.0f, 0.99f) + 1e-10f, 4.0f) - 2, 0.0f, 1.0e9f);
}

// Get obj vertices
inline void get_obj_vertices(const obj_shape& shape, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<int>& vindex,
    bool flipv) {
  auto vmap = unordered_map<obj_vertex, int>{};
  vmap.reserve(shape.vertices.size());
  vindex.reserve(shape.vertices.size());
  for (auto& vert : shape.vertices) {
    auto it = vmap.find(vert);
    if (it != vmap.end()) {
      vindex.push_back(it->second);
      continue;
    }
    auto nverts = (int)positions.size();
    vindex.push_back(nverts);
    vmap.insert(it, {vert, nverts});
    if (!shape.positions.empty() && vert.position)
      positions.push_back(shape.positions[vert.position - 1]);
    if (!shape.normals.empty() && vert.normal)
      normals.push_back(shape.normals[vert.normal - 1]);
    if (!shape.texcoords.empty() && vert.texcoord)
      texcoords.push_back(shape.texcoords[vert.texcoord - 1]);
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
inline void get_obj_triangles(const obj_model& obj, const obj_shape& shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.faces.empty()) return;
  auto vindex = vector<int>{};
  get_obj_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape.materials;
  triangles.reserve(shape.faces.size());
  if (!materials.empty()) ematerials.reserve(shape.faces.size());
  auto cur = 0;
  for (auto& face : shape.faces) {
    for (auto c = 2; c < face.size; c++) {
      triangles.push_back(
          {vindex[cur + 0], vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(face.material);
    }
    cur += face.size;
  }
}
inline void get_obj_quads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.faces.empty()) return;
  auto vindex = vector<int>{};
  get_obj_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape.materials;
  quads.reserve(shape.faces.size());
  if (!materials.empty()) ematerials.reserve(shape.faces.size());
  auto cur = 0;
  for (auto& face : shape.faces) {
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
inline void get_obj_lines(const obj_model& obj, const obj_shape& shape,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.lines.empty()) return;
  auto vindex = vector<int>{};
  get_obj_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape.materials;
  lines.reserve(shape.lines.size());
  if (!materials.empty()) ematerials.reserve(shape.faces.size());
  auto cur = 0;
  for (auto& line : shape.lines) {
    for (auto c = 1; c < line.size; c++) {
      lines.push_back({vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(line.material);
    }
    cur += line.size;
  }
}
inline void get_obj_points(const obj_model& obj, const obj_shape& shape,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.points.empty()) return;
  auto vindex = vector<int>{};
  get_obj_vertices(shape, positions, normals, texcoords, vindex, flipv);
  materials = shape.materials;
  points.reserve(shape.points.size());
  if (!materials.empty()) ematerials.reserve(shape.faces.size());
  auto cur = 0;
  for (auto& point : shape.points) {
    for (auto c = 0; c < point.size; c++) {
      points.push_back({vindex[cur + 0]});
      if (!materials.empty()) ematerials.push_back(point.material);
    }
    cur += point.size;
  }
}
inline void get_obj_fvquads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.faces.empty()) return;
  positions = shape.positions;
  normals   = shape.normals;
  texcoords = flipv ? flip_obj_texcoord(shape.texcoords) : shape.texcoords;
  materials = shape.materials;
  if (shape.vertices[0].position) quadspos.reserve(shape.faces.size());
  if (shape.vertices[0].normal) quadsnorm.reserve(shape.faces.size());
  if (shape.vertices[0].texcoord) quadstexcoord.reserve(shape.faces.size());
  if (!materials.empty()) ematerials.reserve(shape.faces.size());
  auto cur = 0;
  for (auto& face : shape.faces) {
    if (face.size == 4) {
      if (shape.vertices[0].position)
        quadspos.push_back({shape.vertices[cur + 0].position - 1,
            shape.vertices[cur + 1].position - 1,
            shape.vertices[cur + 2].position - 1,
            shape.vertices[cur + 3].position - 1});
      if (shape.vertices[0].normal)
        quadsnorm.push_back({shape.vertices[cur + 0].normal - 1,
            shape.vertices[cur + 1].normal - 1,
            shape.vertices[cur + 2].normal - 1,
            shape.vertices[cur + 3].normal - 1});
      if (shape.vertices[0].texcoord)
        quadstexcoord.push_back({shape.vertices[cur + 0].texcoord - 1,
            shape.vertices[cur + 1].texcoord - 1,
            shape.vertices[cur + 2].texcoord - 1,
            shape.vertices[cur + 3].texcoord - 1});
      if (!materials.empty()) ematerials.push_back(face.material);
    } else {
      for (auto c = 2; c < face.size; c++) {
        if (shape.vertices[0].position)
          quadspos.push_back({shape.vertices[cur + 0].position - 1,
              shape.vertices[cur + c - 1].position - 1,
              shape.vertices[cur + c].position - 1,
              shape.vertices[cur + c].position - 1});
        if (shape.vertices[0].normal)
          quadsnorm.push_back({shape.vertices[cur + 0].normal - 1,
              shape.vertices[cur + c - 1].normal - 1,
              shape.vertices[cur + c].normal - 1,
              shape.vertices[cur + c].normal - 1});
        if (shape.vertices[0].texcoord)
          quadstexcoord.push_back({shape.vertices[cur + 0].texcoord - 1,
              shape.vertices[cur + c - 1].texcoord - 1,
              shape.vertices[cur + c].texcoord - 1,
              shape.vertices[cur + c].texcoord - 1});
        if (!materials.empty()) ematerials.push_back(face.material);
      }
    }
    cur += face.size;
  }
}

inline bool has_obj_quads(const obj_shape& shape) {
  for (auto& face : shape.faces)
    if (face.size == 4) return true;
  return false;
}

// Add obj shape
inline void add_obj_triangles(obj_model& obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials, const vector<int>& ematerials,
    bool flipv) {
  auto& shape     = obj.shapes.emplace_back();
  shape.name      = name;
  shape.materials = materials;
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape.vertices.reserve(triangles.size() * 3);
  for (auto idx = 0; idx < triangles.size(); idx++) {
    auto& triangle = triangles[idx];
    for (auto c = 0; c < 3; c++) {
      shape.vertices.push_back({
          positions.empty() ? 0 : triangle[c] + 1,
          texcoords.empty() ? 0 : triangle[c] + 1,
          normals.empty() ? 0 : triangle[c] + 1,
      });
    }
    shape.faces.push_back(
        {3, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_obj_quads(obj_model& obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials, const vector<int>& ematerials,
    bool flipv) {
  auto& shape     = obj.shapes.emplace_back();
  shape.name      = name;
  shape.materials = materials;
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape.vertices.reserve(quads.size() * 4);
  for (auto idx = 0; idx < quads.size(); idx++) {
    auto& quad = quads[idx];
    auto  nv   = quad.z == quad.w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          positions.empty() ? 0 : quad[c] + 1,
          texcoords.empty() ? 0 : quad[c] + 1,
          normals.empty() ? 0 : quad[c] + 1,
      });
    }
    shape.faces.push_back({(uint8_t)nv,
        ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_obj_lines(obj_model& obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials, const vector<int>& ematerials,
    bool flipv) {
  auto& shape     = obj.shapes.emplace_back();
  shape.name      = name;
  shape.materials = materials;
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape.vertices.reserve(lines.size() * 2);
  for (auto idx = 0; idx < lines.size(); idx++) {
    auto& line = lines[idx];
    for (auto c = 0; c < 2; c++) {
      shape.vertices.push_back({
          positions.empty() ? 0 : line[c] + 1,
          texcoords.empty() ? 0 : line[c] + 1,
          normals.empty() ? 0 : line[c] + 1,
      });
    }
    shape.lines.push_back(
        {2, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_obj_points(obj_model& obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials, const vector<int>& ematerials,
    bool flipv) {
  auto& shape     = obj.shapes.emplace_back();
  shape.name      = name;
  shape.materials = materials;
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape.vertices.reserve(points.size());
  for (auto idx = 0; idx < points.size(); idx++) {
    auto& point = points[idx];
    shape.vertices.push_back({
        positions.empty() ? 0 : point + 1,
        texcoords.empty() ? 0 : point + 1,
        normals.empty() ? 0 : point + 1,
    });
    shape.faces.push_back(
        {1, ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}
inline void add_obj_fvquads(obj_model& obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials, const vector<int>& ematerials,
    bool flipv) {
  auto& shape     = obj.shapes.emplace_back();
  shape.name      = name;
  shape.materials = materials;
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
  shape.vertices.reserve(quadspos.size() * 4);
  for (auto idx = 0; idx < quadspos.size(); idx++) {
    auto nv = quadspos[idx].z == quadspos[idx].w ? 3 : 4;
    for (auto c = 0; c < nv; c++) {
      shape.vertices.push_back({
          quadspos.empty() ? 0 : quadspos[idx][c] + 1,
          quadstexcoord.empty() ? 0 : quadstexcoord[idx][c] + 1,
          quadsnorm.empty() ? 0 : quadsnorm[idx][c] + 1,
      });
    }
    shape.faces.push_back({(uint8_t)nv,
        ematerials.empty() ? (uint8_t)0 : (uint8_t)ematerials[idx]});
  }
}  // namespace yocto

// Read obj
inline bool read_obj_command(obj_file& fs, obj_command& command, string& name,
    vec3f& value, vector<obj_vertex>& vertices, obj_vertex& vert_size) {
  // read the file line by line
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      parse_obj_value(line, value);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      parse_obj_value(line, value);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
      parse_obj_value(line, (vec2f&)value);
      value.z = 0;
      vert_size.texcoord += 1;
      return true;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      skip_obj_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_obj_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        vertices.push_back(vert);
        skip_obj_whitespace(line);
      }
      if (cmd == "f") command = obj_command::face;
      if (cmd == "l") command = obj_command::line;
      if (cmd == "p") command = obj_command::point;
      return true;
    } else if (cmd == "o") {
      command = obj_command::object;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      parse_obj_value_or_empty(line, name);
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      parse_obj_value(line, name);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
inline bool read_mtl_command(
    obj_file& fs, mtl_command& command, obj_material& material, bool fliptr) {
  material = {};

  // read the file line by line
  auto pos   = ftell(fs.fs);
  auto found = false;
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      if (found) {
        command = mtl_command::material;
        fseek(fs.fs, pos, SEEK_SET);
        return true;
      } else {
        found = true;
      }
      parse_obj_value(line, material.name);
    } else if (cmd == "illum") {
      parse_obj_value(line, material.illum);
    } else if (cmd == "Ke") {
      parse_obj_value(line, material.emission);
    } else if (cmd == "Ka") {
      parse_obj_value(line, material.ambient);
    } else if (cmd == "Kd") {
      parse_obj_value(line, material.diffuse);
    } else if (cmd == "Ks") {
      parse_obj_value(line, material.specular);
    } else if (cmd == "Kt") {
      parse_obj_value(line, material.transmission);
    } else if (cmd == "Tf") {
      material.transmission = vec3f{-1};
      parse_obj_value(line, material.transmission);
      if (material.transmission.y < 0)
        material.transmission = vec3f{material.transmission.x};
      if (fliptr) material.transmission = 1 - material.transmission;
    } else if (cmd == "Tr") {
      parse_obj_value(line, material.opacity);
      if (fliptr) material.opacity = 1 - material.opacity;
    } else if (cmd == "Ns") {
      parse_obj_value(line, material.exponent);
    } else if (cmd == "d") {
      parse_obj_value(line, material.opacity);
    } else if (cmd == "map_Ke") {
      parse_obj_value(line, material.emission_map);
    } else if (cmd == "map_Ka") {
      parse_obj_value(line, material.ambient_map);
    } else if (cmd == "map_Kd") {
      parse_obj_value(line, material.diffuse_map);
    } else if (cmd == "map_Ks") {
      parse_obj_value(line, material.specular_map);
    } else if (cmd == "map_Tr") {
      parse_obj_value(line, material.transmission_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_obj_value(line, material.opacity_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_obj_value(line, material.bump_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_obj_value(line, material.displacement_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_obj_value(line, material.normal_map);
    } else if (cmd == "Pm") {
      parse_obj_value(line, material.pbr_metallic);
    } else if (cmd == "Pr") {
      parse_obj_value(line, material.pbr_roughness);
    } else if (cmd == "Ps") {
      parse_obj_value(line, material.pbr_sheen);
    } else if (cmd == "Pc") {
      parse_obj_value(line, material.pbr_clearcoat);
    } else if (cmd == "Pcr") {
      parse_obj_value(line, material.pbr_coatroughness);
    } else if (cmd == "map_Pm") {
      parse_obj_value(line, material.pbr_metallic_map);
    } else if (cmd == "map_Pr") {
      parse_obj_value(line, material.pbr_roughness_map);
    } else if (cmd == "map_Ps") {
      parse_obj_value(line, material.pbr_sheen_map);
    } else if (cmd == "map_Pc") {
      parse_obj_value(line, material.pbr_clearcoat_map);
    } else if (cmd == "map_Pcr") {
      parse_obj_value(line, material.pbr_coatroughness_map);
    } else if (cmd == "Vt") {
      parse_obj_value(line, material.vol_transmission);
    } else if (cmd == "Vp") {
      parse_obj_value(line, material.vol_meanfreepath);
    } else if (cmd == "Ve") {
      parse_obj_value(line, material.vol_emission);
    } else if (cmd == "Vs") {
      parse_obj_value(line, material.vol_scattering);
    } else if (cmd == "Vg") {
      parse_obj_value(line, material.vol_anisotropy);
    } else if (cmd == "Vr") {
      parse_obj_value(line, material.vol_scale);
    } else if (cmd == "map_Vs") {
      parse_obj_value(line, material.vol_scattering_map);
    } else {
      continue;
    }
    pos = ftell(fs.fs);
  }

  if (found) {
    command = mtl_command::material;
    return true;
  }

  return false;
}

// Read objx
inline bool read_objx_command(obj_file& fs, objx_command& command,
    obj_camera& camera, obj_environment& environment, obj_instance& instance) {
  // read the file line by line
  char buffer[4096];
  auto found = false;
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_obj_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      command = objx_command::camera;
      parse_obj_value(line, camera.name);
      parse_obj_value(line, camera.ortho);
      parse_obj_value(line, camera.width);
      parse_obj_value(line, camera.height);
      parse_obj_value(line, camera.lens);
      parse_obj_value(line, camera.focus);
      parse_obj_value(line, camera.aperture);
      parse_obj_value(line, camera.frame);
      return true;
    } else if (cmd == "e") {
      command = objx_command::environment;
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.emission);
      parse_obj_value(line, environment.emission_map);
      parse_obj_value(line, environment.frame);
      return true;
    } else if (cmd == "i") {
      command = objx_command::instance;
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.frame);
      return true;
    }
  }

  if (found) return true;

  return false;
}

inline vector<string> split_obj_string(const string& str, const string& delim) {
  auto tokens = vector<string>{};
  auto last = (size_t)0, next = (size_t)0;
  while ((next = str.find(delim, last)) != string::npos) {
    tokens.push_back(str.substr(last, next - last));
    last = next + delim.size();
  }
  if (last < str.size()) tokens.push_back(str.substr(last));
  return tokens;
}

// Write obj elements
inline void write_obj_comment(obj_file& fs, const string& comment) {
  auto lines = split_obj_string(comment, "\n");
  for (auto& line : lines) {
    format_obj_values(fs, "# {}\n", line);
  }
  format_obj_values(fs, "\n");
}

inline void write_obj_command(obj_file& fs, obj_command command,
    const string& name, const vec3f& value,
    const vector<obj_vertex>& vertices) {
  switch (command) {
    case obj_command::vertex: format_obj_values(fs, "v {}\n", value); break;
    case obj_command::normal: format_obj_values(fs, "vn {}\n", value); break;
    case obj_command::texcoord: format_obj_values(fs, "vt {}\n", value); break;
    case obj_command::face:
    case obj_command::line:
    case obj_command::point:
      if (command == obj_command::face) format_obj_values(fs, "f ");
      if (command == obj_command::line) format_obj_values(fs, "l ");
      if (command == obj_command::point) format_obj_values(fs, "p ");
      for (auto& vert : vertices) format_obj_values(fs, " {}", vert);
      format_obj_values(fs, "\n");
      break;
    case obj_command::object: format_obj_values(fs, "o {}\n", name); break;
    case obj_command::group: format_obj_values(fs, "g {}\n", name); break;
    case obj_command::usemtl: format_obj_values(fs, "usemtl {}\n", name); break;
    case obj_command::smoothing: format_obj_values(fs, "s {}\n", name); break;
    case obj_command::mtllib: format_obj_values(fs, "mtllib {}\n", name); break;
    case obj_command::objxlib: break;
  }
}

inline void write_mtl_command(
    obj_file& fs, mtl_command command, const obj_material& material) {
  // write material
  format_obj_values(fs, "newmtl {}\n", material.name);
  format_obj_values(fs, "illum {}\n", material.illum);
  if (material.emission != zero3f)
    format_obj_values(fs, "Ke {}\n", material.emission);
  if (material.ambient != zero3f)
    format_obj_values(fs, "Ka {}\n", material.ambient);
  format_obj_values(fs, "Kd {}\n", material.diffuse);
  format_obj_values(fs, "Ks {}\n", material.specular);
  if (material.reflection != zero3f)
    format_obj_values(fs, "Kr {}\n", material.reflection);
  if (material.transmission != zero3f)
    format_obj_values(fs, "Kt {}\n", material.transmission);
  format_obj_values(fs, "Ns {}\n", (int)material.exponent);
  if (material.opacity != 1) format_obj_values(fs, "d {}\n", material.opacity);
  if (!material.emission_map.path.empty())
    format_obj_values(fs, "map_Ke {}\n", material.emission_map);
  if (!material.diffuse_map.path.empty())
    format_obj_values(fs, "map_Kd {}\n", material.diffuse_map);
  if (!material.specular_map.path.empty())
    format_obj_values(fs, "map_Ks {}\n", material.specular_map);
  if (!material.transmission_map.path.empty())
    format_obj_values(fs, "map_Kt {}\n", material.transmission_map);
  if (!material.reflection_map.path.empty())
    format_obj_values(fs, "map_Kr {}\n", material.reflection_map);
  if (!material.exponent_map.path.empty())
    format_obj_values(fs, "map_Ns {}\n", material.exponent_map);
  if (!material.opacity_map.path.empty())
    format_obj_values(fs, "map_d {}\n", material.opacity_map);
  if (!material.bump_map.path.empty())
    format_obj_values(fs, "map_bump {}\n", material.bump_map);
  if (!material.displacement_map.path.empty())
    format_obj_values(fs, "map_disp {}\n", material.displacement_map);
  if (!material.normal_map.path.empty())
    format_obj_values(fs, "map_norm {}\n", material.normal_map);
  if (material.pbr_roughness)
    format_obj_values(fs, "Pr {}\n", material.pbr_roughness);
  if (material.pbr_metallic)
    format_obj_values(fs, "Pm {}\n", material.pbr_metallic);
  if (material.pbr_sheen) format_obj_values(fs, "Ps {}\n", material.pbr_sheen);
  if (material.pbr_clearcoat)
    format_obj_values(fs, "Pc {}\n", material.pbr_clearcoat);
  if (material.pbr_coatroughness)
    format_obj_values(fs, "Pcr {}\n", material.pbr_coatroughness);
  if (!material.pbr_roughness_map.path.empty())
    format_obj_values(fs, "map_Pr {}\n", material.pbr_roughness_map);
  if (!material.pbr_metallic_map.path.empty())
    format_obj_values(fs, "map_Pm {}\n", material.pbr_metallic_map);
  if (!material.pbr_sheen_map.path.empty())
    format_obj_values(fs, "map_Ps {}\n", material.pbr_sheen_map);
  if (!material.pbr_clearcoat_map.path.empty())
    format_obj_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map);
  if (!material.pbr_coatroughness_map.path.empty())
    format_obj_values(fs, "map_Pcr {}\n", material.pbr_coatroughness_map);
  if (material.vol_transmission != zero3f)
    format_obj_values(fs, "Vt {}\n", material.vol_transmission);
  if (material.vol_meanfreepath != zero3f)
    format_obj_values(fs, "Vp {}\n", material.vol_meanfreepath);
  if (material.vol_emission != zero3f)
    format_obj_values(fs, "Ve {}\n", material.vol_emission);
  if (material.vol_scattering != zero3f)
    format_obj_values(fs, "Vs {}\n", material.vol_scattering);
  if (material.vol_anisotropy)
    format_obj_values(fs, "Vg {}\n", material.vol_anisotropy);
  if (material.vol_scale) format_obj_values(fs, "Vr {}\n", material.vol_scale);
  if (!material.vol_scattering_map.path.empty())
    format_obj_values(fs, "map_Vs {}\n", material.vol_scattering_map);
  format_obj_values(fs, "\n");
}

inline void write_objx_command(obj_file& fs, objx_command command,
    const obj_camera& camera, const obj_environment& environment,
    const obj_instance& instance) {
  switch (command) {
    case objx_command::camera: {
      format_obj_values(fs, "c {} {} {} {} {} {} {} {}\n", camera.name,
          camera.ortho, camera.width, camera.height, camera.lens, camera.focus,
          camera.aperture, camera.frame);
    } break;
    case objx_command::environment: {
      format_obj_values(fs, "e {} {} {} {}\n", environment.name,
          environment.emission,
          environment.emission_map.path.empty() ? "\"\""s
                                                : environment.emission_map.path,
          environment.frame);
    } break;
    case objx_command::instance: {
      format_obj_values(fs, "i {} {}\n", instance.object, instance.frame);
    } break;
  }
}

}  // namespace yocto

#endif
