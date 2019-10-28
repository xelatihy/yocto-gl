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

  operator bool() const { return (bool)fs; }

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
  face, str, point,                // data in vertices
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

inline void write_obj_text(obj_file& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
inline void write_obj_text(obj_file& fs, const char* value) {
  if (fputs(value, fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
inline string normalize_obj_path(const string& filename_) {
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
inline string get_obj_dirname(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
inline string get_obj_extension(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
inline string get_obj_filename(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
inline string get_obj_noextension(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
inline string get_obj_basename(const string& filename) {
  return get_obj_noextension(get_obj_filename(filename));
}

// Replaces extensions
inline string replace_obj_extension(const string& filename, const string& ext) {
  return get_obj_noextension(filename) + ext;
}

// Check if a file can be opened for reading.
inline bool exists_obj_file(const string& filename) {
  auto fs = fopen(filename.c_str(), "r");
  if (fs) {
    fclose(fs);
    return true;
  } else {
    return false;
  }
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
[[nodiscard]] inline string_view skip_obj_whitespace(string_view str) {
  while (!str.empty() && is_obj_space(str.front())) str.remove_prefix(1);
  return str;
}

[[nodiscard]] inline string_view remove_obj_comment(
    string_view str, char comment_char = '#') {
  while (!str.empty() && is_obj_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
  return str;
}

// Parse values from a string
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, string_view& value) {
  str = skip_obj_whitespace(str);
  if (str.empty()) return {};
  auto cpy = str;
  while (!cpy.empty() && !is_obj_space(cpy.front())) cpy.remove_prefix(1);
  value = str;
  value.remove_suffix(cpy.size());
  str.remove_prefix(str.size() - cpy.size());
  return str;
}
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, string& value) {
  auto valuev = string_view{};
  str         = parse_obj_value(str, valuev);
  if (!str.data()) return {};
  value = string{valuev};
  return str;
}
[[nodiscard]] inline string_view parse_obj_value(string_view str, int& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return {};
  str.remove_prefix(end - str.data());
  return str;
}
[[nodiscard]] inline string_view parse_obj_value(string_view str, bool& value) {
  auto valuei = 0;
  str         = parse_obj_value(str, valuei);
  if (!str.data()) return {};
  value = (bool)valuei;
  return str;
}
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return {};
  str.remove_prefix(end - str.data());
  return str;
}
template <typename T, size_t N>
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, array<T, N>& value) {
  for (auto i = 0; i < value.size(); i++) {
    str = parse_obj_value(str, value[i]);
    if (!str.data()) return {};
  }
  return str;
}

[[nodiscard]] inline string_view parse_obj_value(
    string_view str, vec2f& value) {
  return parse_obj_value(str, reinterpret_cast<array<float, 2>&>(value));
}
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, vec3f& value) {
  return parse_obj_value(str, reinterpret_cast<array<float, 3>&>(value));
}
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, frame3f& value) {
  return parse_obj_value(str, reinterpret_cast<array<float, 12>&>(value));
}

// Parse values from a string
[[nodiscard]] inline string_view parse_obj_value_or_empty(
    string_view& str, string& value) {
  str = skip_obj_whitespace(str);
  if (str.empty()) {
    value = "";
    return str;
  } else {
    return parse_obj_value(str, value);
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
  for (auto i = 0; i < 2; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
}
inline void format_obj_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
}
inline void format_obj_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
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

[[nodiscard]] inline string_view parse_obj_value(
    string_view str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  str   = parse_obj_value(str, value.position);
  if (!str.data()) return {};
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      str = parse_obj_value(str, value.normal);
      if (!str.data()) return {};
    } else {
      str = parse_obj_value(str, value.texcoord);
      if (!str.data()) return {};
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        str = parse_obj_value(str, value.normal);
        if (!str.data()) return {};
      }
    }
  }
  return str;
}

// Input for OBJ textures
[[nodiscard]] inline string_view parse_obj_value(
    string_view str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  str         = skip_obj_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    str        = parse_obj_value(str, token);
    tokens.push_back(token);
    str = skip_obj_whitespace(str);
  }
  if (tokens.empty()) throw std::runtime_error("cannot parse value");

  // texture name
  info.path = normalize_obj_path(tokens.back());

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }

  return str;
}

// Read obj
inline void load_mtl(
    const string& filename, obj_model& obj, bool fliptr = true) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) throw std::runtime_error("cannot open " + filename);

  // init parsing
  obj.materials.emplace_back();

  // initialize parsing
  auto parse_error = [&filename](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
    if (parse_error(str)) return;
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      obj.materials.emplace_back();
      str = parse_obj_value(str, obj.materials.back().name);
      if (parse_error(str)) return;
    } else if (cmd == "illum") {
      str = parse_obj_value(str, obj.materials.back().illum);
      if (parse_error(str)) return;
    } else if (cmd == "Ke") {
      str = parse_obj_value(str, obj.materials.back().emission);
      if (parse_error(str)) return;
    } else if (cmd == "Ka") {
      str = parse_obj_value(str, obj.materials.back().ambient);
      if (parse_error(str)) return;
    } else if (cmd == "Kd") {
      str = parse_obj_value(str, obj.materials.back().diffuse);
      if (parse_error(str)) return;
    } else if (cmd == "Ks") {
      str = parse_obj_value(str, obj.materials.back().specular);
      if (parse_error(str)) return;
    } else if (cmd == "Kt") {
      str = parse_obj_value(str, obj.materials.back().transmission);
      if (parse_error(str)) return;
    } else if (cmd == "Tf") {
      obj.materials.back().transmission = vec3f{-1};
      str = parse_obj_value(str, obj.materials.back().transmission);
      if (obj.materials.back().transmission.y < 0)
        obj.materials.back().transmission = vec3f{
            obj.materials.back().transmission.x};
      if (fliptr)
        obj.materials.back().transmission = 1 -
                                            obj.materials.back().transmission;
      if (parse_error(str)) return;
    } else if (cmd == "Tr") {
      str = parse_obj_value(str, obj.materials.back().opacity);
      if (fliptr)
        obj.materials.back().opacity = 1 - obj.materials.back().opacity;
      if (parse_error(str)) return;
    } else if (cmd == "Ns") {
      str = parse_obj_value(str, obj.materials.back().exponent);
      if (parse_error(str)) return;
    } else if (cmd == "d") {
      str = parse_obj_value(str, obj.materials.back().opacity);
      if (parse_error(str)) return;
    } else if (cmd == "map_Ke") {
      str = parse_obj_value(str, obj.materials.back().emission_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Ka") {
      str = parse_obj_value(str, obj.materials.back().ambient_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Kd") {
      str = parse_obj_value(str, obj.materials.back().diffuse_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Ks") {
      str = parse_obj_value(str, obj.materials.back().specular_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Tr") {
      str = parse_obj_value(str, obj.materials.back().transmission_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      str = parse_obj_value(str, obj.materials.back().opacity_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_bump" || cmd == "bump") {
      str = parse_obj_value(str, obj.materials.back().bump_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_disp" || cmd == "disp") {
      str = parse_obj_value(str, obj.materials.back().displacement_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_norm" || cmd == "norm") {
      str = parse_obj_value(str, obj.materials.back().normal_map);
      if (parse_error(str)) return;
    } else if (cmd == "Pm") {
      str = parse_obj_value(str, obj.materials.back().pbr_metallic);
      if (parse_error(str)) return;
    } else if (cmd == "Pr") {
      str = parse_obj_value(str, obj.materials.back().pbr_roughness);
      if (parse_error(str)) return;
    } else if (cmd == "Ps") {
      str = parse_obj_value(str, obj.materials.back().pbr_sheen);
      if (parse_error(str)) return;
    } else if (cmd == "Pc") {
      str = parse_obj_value(str, obj.materials.back().pbr_clearcoat);
      if (parse_error(str)) return;
    } else if (cmd == "Pcr") {
      str = parse_obj_value(str, obj.materials.back().pbr_coatroughness);
      if (parse_error(str)) return;
    } else if (cmd == "map_Pm") {
      str = parse_obj_value(str, obj.materials.back().pbr_metallic_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Pr") {
      str = parse_obj_value(str, obj.materials.back().pbr_roughness_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Ps") {
      str = parse_obj_value(str, obj.materials.back().pbr_sheen_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Pc") {
      str = parse_obj_value(str, obj.materials.back().pbr_clearcoat_map);
      if (parse_error(str)) return;
    } else if (cmd == "map_Pcr") {
      str = parse_obj_value(str, obj.materials.back().pbr_coatroughness_map);
      if (parse_error(str)) return;
    } else if (cmd == "Vt") {
      str = parse_obj_value(str, obj.materials.back().vol_transmission);
      if (parse_error(str)) return;
    } else if (cmd == "Vp") {
      str = parse_obj_value(str, obj.materials.back().vol_meanfreepath);
      if (parse_error(str)) return;
    } else if (cmd == "Ve") {
      str = parse_obj_value(str, obj.materials.back().vol_emission);
      if (parse_error(str)) return;
    } else if (cmd == "Vs") {
      str = parse_obj_value(str, obj.materials.back().vol_scattering);
      if (parse_error(str)) return;
    } else if (cmd == "Vg") {
      str = parse_obj_value(str, obj.materials.back().vol_anisotropy);
      if (parse_error(str)) return;
    } else if (cmd == "Vr") {
      str = parse_obj_value(str, obj.materials.back().vol_scale);
      if (parse_error(str)) return;
    } else if (cmd == "map_Vs") {
      str = parse_obj_value(str, obj.materials.back().vol_scattering_map);
      if (parse_error(str)) return;
    } else {
      continue;
    }
  }

  // check error
  if (read_error(fs)) return;

  // remove placeholder material
  obj.materials.erase(obj.materials.begin());
}

// Read obj
inline void load_objx(const string& filename, obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) throw std::runtime_error("cannot open " + filename);

  // shape map for instances
  auto shape_map = unordered_map<string, vector<int>>{};
  for (auto idx = 0; idx < obj.shapes.size(); idx++) {
    shape_map[obj.shapes[idx].name].push_back(idx);
  }

  // initialize parsing
  auto parse_error = [&filename](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
    if (parse_error(str)) return;
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto& camera = obj.cameras.emplace_back();
      str          = parse_obj_value(str, camera.name);
      str          = parse_obj_value(str, camera.ortho);
      str          = parse_obj_value(str, camera.width);
      str          = parse_obj_value(str, camera.height);
      str          = parse_obj_value(str, camera.lens);
      str          = parse_obj_value(str, camera.focus);
      str          = parse_obj_value(str, camera.aperture);
      str          = parse_obj_value(str, camera.frame);
      if (parse_error(str)) return;
    } else if (cmd == "e") {
      auto& environment  = obj.environments.emplace_back();
      str                = parse_obj_value(str, environment.name);
      str                = parse_obj_value(str, environment.emission);
      auto emission_path = ""s;
      str                = parse_obj_value(str, emission_path);
      if (emission_path == "\"\"") emission_path = "";
      environment.emission_map.path = emission_path;
      str                           = parse_obj_value(str, environment.frame);
      if (parse_error(str)) return;
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      str         = parse_obj_value(str, object);
      str         = parse_obj_value(str, frame);
      if (shape_map.find(object) == shape_map.end()) {
        throw std::runtime_error("cannot find object " + object);
      }
      for (auto idx : shape_map.at(object)) {
        obj.shapes[idx].instances.push_back(frame);
      }
      if (parse_error(str)) return;
    } else {
      // unused
    }
  }

  // check error
  if (read_error(fs)) return;
}

// Read obj
inline void load_obj(const string& filename, obj_model& obj, bool geom_only,
    bool split_elements, bool split_materials) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) throw std::runtime_error("cannot open " + filename);

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

  // initialize parsing
  auto parse_error = [&filename](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
    if (parse_error(str)) return;
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      str = parse_obj_value(str, opositions.emplace_back());
      vert_size.position += 1;
      if (parse_error(str)) return;
    } else if (cmd == "vn") {
      str = parse_obj_value(str, onormals.emplace_back());
      vert_size.normal += 1;
      if (parse_error(str)) return;
    } else if (cmd == "vt") {
      str = parse_obj_value(str, otexcoords.emplace_back());
      vert_size.texcoord += 1;
      if (parse_error(str)) return;
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
      str = skip_obj_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        str       = parse_obj_value(str, vert);
        if (parse_error(str)) return;
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        shape.vertices.push_back(vert);
        element.size += 1;
        str = skip_obj_whitespace(str);
      }
      if (parse_error(str)) return;
    } else if (cmd == "o" || cmd == "g") {
      if (geom_only) continue;
      str = parse_obj_value_or_empty(str, cmd == "o" ? oname : gname);
      if (!obj.shapes.back().vertices.empty()) {
        obj.shapes.emplace_back();
        obj.shapes.back().name = oname + gname;
      } else {
        obj.shapes.back().name = oname + gname;
      }
      if (parse_error(str)) return;
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      str = parse_obj_value_or_empty(str, mname);
      if (parse_error(str)) return;
    } else if (cmd == "s") {
      if (geom_only) continue;
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      str         = parse_obj_value(str, mtllib);
      if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) == mtllibs.end()) {
        mtllibs.push_back(mtllib);
      }
      if (parse_error(str)) return;
    } else {
      // unused
    }
  }

  // check error
  if (read_error(fs)) return;

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
  auto dirname = get_obj_dirname(filename);
  for (auto& mtllib : mtllibs) {
    load_mtl(dirname + mtllib, obj);
  }

  // load extensions
  auto extfilename = replace_obj_extension(filename, ".objx");
  if (exists_obj_file(extfilename)) {
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
  if (!fs) throw std::runtime_error("cannot open " + filename);

  auto write_error = [&filename](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };

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

  // check error
  if (write_error(fs)) return;
}

// Save obj
inline void save_objx(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");
  if (!fs) throw std::runtime_error("cannot open " + filename);

  auto write_error = [&filename](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };

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

  // check error
  if (write_error(fs)) return;
}

// Save obj
inline void save_obj(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");
  if (!fs) throw std::runtime_error("cannot open " + filename);

  auto write_error = [&filename](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + filename);
    return true;
  };

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
    format_obj_values(fs, "mtllib {}\n\n",
        replace_obj_extension(get_obj_filename(filename), ".mtl"));
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
    save_mtl(replace_obj_extension(filename, ".mtl"), obj);

  // save objx
  if (!obj.cameras.empty() || !obj.environments.empty() ||
      std::any_of(obj.shapes.begin(), obj.shapes.end(),
          [](auto& shape) { return !shape.instances.empty(); }))
    save_objx(replace_obj_extension(filename, ".objx"), obj);

  // check error
  if (write_error(fs)) return;
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
  for (auto& str : shape.lines) {
    for (auto c = 1; c < str.size; c++) {
      lines.push_back({vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(str.material);
    }
    cur += str.size;
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
    auto& str = lines[idx];
    for (auto c = 0; c < 2; c++) {
      shape.vertices.push_back({
          positions.empty() ? 0 : str[c] + 1,
          texcoords.empty() ? 0 : str[c] + 1,
          normals.empty() ? 0 : str[c] + 1,
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
  // initialize parsing
  auto parse_error = [&fs](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
    if (parse_error(str)) return false;
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      str     = parse_obj_value(str, value);
      if (parse_error(str)) return false;
      vert_size.position += 1;
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      str     = parse_obj_value(str, value);
      if (parse_error(str)) return false;
      vert_size.normal += 1;
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
      str     = parse_obj_value(str, (vec2f&)value);
      if (parse_error(str)) return false;
      value.z = 0;
      vert_size.texcoord += 1;
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      str = skip_obj_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        str       = parse_obj_value(str, vert);
        if (parse_error(str)) return false;
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        vertices.push_back(vert);
        str = skip_obj_whitespace(str);
      }
      if (cmd == "f") command = obj_command::face;
      if (cmd == "l") command = obj_command::str;
      if (cmd == "p") command = obj_command::point;
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "o") {
      command = obj_command::object;
      str     = parse_obj_value_or_empty(str, name);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      str     = parse_obj_value_or_empty(str, name);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      str     = parse_obj_value_or_empty(str, name);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      str     = parse_obj_value_or_empty(str, name);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      str     = parse_obj_value(str, name);
      if (parse_error(str)) return false;
      if (parse_error(str)) return false;
      return true;
    } else {
      // unused
    }
  }

  // check error
  if (read_error(fs)) return false;

  return false;
}

// Read mtl
inline bool read_mtl_command(
    obj_file& fs, mtl_command& command, obj_material& material, bool fliptr) {
  material = {};

  // initialize parsing
  auto parse_error = [&fs](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  auto pos   = ftell(fs.fs);
  auto found = false;
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
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
      str = parse_obj_value(str, material.name);
      if (parse_error(str)) return false;
    } else if (cmd == "illum") {
      str = parse_obj_value(str, material.illum);
      if (parse_error(str)) return false;
    } else if (cmd == "Ke") {
      str = parse_obj_value(str, material.emission);
      if (parse_error(str)) return false;
    } else if (cmd == "Ka") {
      str = parse_obj_value(str, material.ambient);
      if (parse_error(str)) return false;
    } else if (cmd == "Kd") {
      str = parse_obj_value(str, material.diffuse);
      if (parse_error(str)) return false;
    } else if (cmd == "Ks") {
      str = parse_obj_value(str, material.specular);
      if (parse_error(str)) return false;
    } else if (cmd == "Kt") {
      str = parse_obj_value(str, material.transmission);
      if (parse_error(str)) return false;
    } else if (cmd == "Tf") {
      material.transmission = vec3f{-1};
      str                   = parse_obj_value(str, material.transmission);
      if (parse_error(str)) return false;
      if (material.transmission.y < 0)
        material.transmission = vec3f{material.transmission.x};
      if (fliptr) material.transmission = 1 - material.transmission;
    } else if (cmd == "Tr") {
      str = parse_obj_value(str, material.opacity);
      if (parse_error(str)) return false;
      if (fliptr) material.opacity = 1 - material.opacity;
    } else if (cmd == "Ns") {
      str = parse_obj_value(str, material.exponent);
      if (parse_error(str)) return false;
    } else if (cmd == "d") {
      str = parse_obj_value(str, material.opacity);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Ke") {
      str = parse_obj_value(str, material.emission_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Ka") {
      str = parse_obj_value(str, material.ambient_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Kd") {
      str = parse_obj_value(str, material.diffuse_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Ks") {
      str = parse_obj_value(str, material.specular_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Tr") {
      str = parse_obj_value(str, material.transmission_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      str = parse_obj_value(str, material.opacity_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_bump" || cmd == "bump") {
      str = parse_obj_value(str, material.bump_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_disp" || cmd == "disp") {
      str = parse_obj_value(str, material.displacement_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_norm" || cmd == "norm") {
      str = parse_obj_value(str, material.normal_map);
      if (parse_error(str)) return false;
    } else if (cmd == "Pm") {
      str = parse_obj_value(str, material.pbr_metallic);
      if (parse_error(str)) return false;
    } else if (cmd == "Pr") {
      str = parse_obj_value(str, material.pbr_roughness);
      if (parse_error(str)) return false;
    } else if (cmd == "Ps") {
      str = parse_obj_value(str, material.pbr_sheen);
      if (parse_error(str)) return false;
    } else if (cmd == "Pc") {
      str = parse_obj_value(str, material.pbr_clearcoat);
      if (parse_error(str)) return false;
    } else if (cmd == "Pcr") {
      str = parse_obj_value(str, material.pbr_coatroughness);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Pm") {
      str = parse_obj_value(str, material.pbr_metallic_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Pr") {
      str = parse_obj_value(str, material.pbr_roughness_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Ps") {
      str = parse_obj_value(str, material.pbr_sheen_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Pc") {
      str = parse_obj_value(str, material.pbr_clearcoat_map);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Pcr") {
      str = parse_obj_value(str, material.pbr_coatroughness_map);
      if (parse_error(str)) return false;
    } else if (cmd == "Vt") {
      str = parse_obj_value(str, material.vol_transmission);
      if (parse_error(str)) return false;
    } else if (cmd == "Vp") {
      str = parse_obj_value(str, material.vol_meanfreepath);
      if (parse_error(str)) return false;
    } else if (cmd == "Ve") {
      str = parse_obj_value(str, material.vol_emission);
      if (parse_error(str)) return false;
    } else if (cmd == "Vs") {
      str = parse_obj_value(str, material.vol_scattering);
      if (parse_error(str)) return false;
    } else if (cmd == "Vg") {
      str = parse_obj_value(str, material.vol_anisotropy);
      if (parse_error(str)) return false;
    } else if (cmd == "Vr") {
      str = parse_obj_value(str, material.vol_scale);
      if (parse_error(str)) return false;
    } else if (cmd == "map_Vs") {
      str = parse_obj_value(str, material.vol_scattering_map);
      if (parse_error(str)) return false;
    } else {
      continue;
    }
    pos = ftell(fs.fs);
  }

  if (found) {
    command = mtl_command::material;
    return true;
  }

  // check error
  if (read_error(fs)) return false;

  return false;
}

// Read objx
inline bool read_objx_command(obj_file& fs, objx_command& command,
    obj_camera& camera, obj_environment& environment, obj_instance& instance) {
  // initialize parsing
  auto parse_error = [&fs](string_view str) {
    if (str.data()) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };
  auto read_error = [](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

  // read the file str by str
  char buffer[4096];
  auto found = false;
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    str      = remove_obj_comment(str);
    str      = skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    str      = parse_obj_value(str, cmd);
    if (parse_error(str)) return false;
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      command = objx_command::camera;
      str     = parse_obj_value(str, camera.name);
      str     = parse_obj_value(str, camera.ortho);
      str     = parse_obj_value(str, camera.width);
      str     = parse_obj_value(str, camera.height);
      str     = parse_obj_value(str, camera.lens);
      str     = parse_obj_value(str, camera.focus);
      str     = parse_obj_value(str, camera.aperture);
      str     = parse_obj_value(str, camera.frame);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "e") {
      command = objx_command::environment;
      str     = parse_obj_value(str, environment.name);
      str     = parse_obj_value(str, environment.emission);
      str     = parse_obj_value(str, environment.emission_map);
      str     = parse_obj_value(str, environment.frame);
      if (parse_error(str)) return false;
      return true;
    } else if (cmd == "i") {
      command = objx_command::instance;
      str     = parse_obj_value(str, instance.object);
      str     = parse_obj_value(str, instance.frame);
      if (parse_error(str)) return false;
      return true;
    }
  }

  if (found) return true;

  // check error
  if (read_error(fs)) return false;

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
  for (auto& str : lines) {
    format_obj_values(fs, "# {}\n", str);
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
    case obj_command::str:
    case obj_command::point:
      if (command == obj_command::face) format_obj_values(fs, "f ");
      if (command == obj_command::str) format_obj_values(fs, "l ");
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
  auto write_error = [&](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

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

  // check error
  if (write_error(fs)) return;
}

inline void write_objx_command(obj_file& fs, objx_command command,
    const obj_camera& camera, const obj_environment& environment,
    const obj_instance& instance) {
  auto write_error = [&](obj_file& fs) {
    if (!ferror(fs.fs)) return false;
    throw std::runtime_error("cannot parse " + fs.filename);
    return true;
  };

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

  // check error
  if (write_error(fs)) return;
}

}  // namespace yocto

#endif
