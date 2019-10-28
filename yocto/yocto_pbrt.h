//
// # Yocto/Pbrt: Tiny library for Obj parsing and writing
//
// Yocto/Pbrt is a tiny library for loading and saving Pbrt file. Yocto/Pbrt
// supports two interfaces: a simple interface where all Pbrt data is loaded
// and saved at once and a low-level interface where Pbrt commands are read
// and written one at a time. In the high-level interface we provice a
// simplified representation for Pbrt types such as cameras, mnaterials and
// shapes.
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

#ifndef _YOCTO_PBRT_H_
#define _YOCTO_PBRT_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt value type
enum struct pbrt_value_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color, 
  point2, vector2, spectrum
  // clang-format on
};

// Pbrt value
struct pbrt_value {
  string          name     = "";
  pbrt_value_type type     = pbrt_value_type::real;
  int             value1i  = 0;
  float           value1f  = 0;
  vec2f           value2f  = {0, 0};
  vec3f           value3f  = {0, 0, 0};
  bool            value1b  = false;
  string          value1s  = "";
  vector<float>   vector1f = {};
  vector<vec2f>   vector2f = {};
  vector<vec3f>   vector3f = {};
  vector<int>     vector1i = {};
};

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // camera approximation
  float width    = 0;
  float height   = 0;
  float lens     = 0;
  float aspect   = 0;
  float focus    = 0;
  float aperture = 0;
};

// Pbrt texture
struct pbrt_texture {
  // texture parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  // texture approximation
  vec3f  constant = vec3f{1, 1, 1};
  string filename = "";
};

// Pbrt material
struct pbrt_material {
  // material parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  // material approximation
  vec3f  diffuse          = zero3f;
  vec3f  specular         = zero3f;
  vec3f  transmission     = zero3f;
  vec2f  roughness        = zero2f;
  vec3f  opacity          = vec3f{1};
  vec3f  eta              = zero3f;
  vec3f  etak             = zero3f;
  vec3f  sspecular        = zero3f;  // specular scaled by fresnel
  string diffuse_map      = "";
  string specular_map     = "";
  string transmission_map = "";
  string roughness_map    = "";
  string opacity_map      = "";
  string eta_map          = "";
  string etak_map         = "";
  vec3f  volmeanfreepath  = vec3f{0};
  vec3f  volscatter       = vec3f{0};
  float  volscale         = 0.01;
  bool   refract          = false;
};

// Pbrt medium
struct pbrt_medium {
  // medium parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
};

// Pbrt shape
struct pbrt_shape {
  // shape parameters
  string             type            = "";
  vector<pbrt_value> values          = {};
  frame3f            frame           = identity3x4f;
  frame3f            frend           = identity3x4f;
  string             material        = "";
  string             arealight       = "";
  string             interior        = "";
  string             exterior        = "";
  bool               is_instanced    = false;
  vector<frame3f>    instance_frames = {};
  vector<frame3f>    instance_frends = {};
  // shape approximation
  string        filename  = "";
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec3i> triangles = {};
  float         radius    = 0;  // radius for sphere, cylinder, disk
};

// Pbrt lights
struct pbrt_light {
  // light parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // light approximation
  vec3f emission = zero3f;
  vec3f from     = zero3f;
  vec3f to       = zero3f;
  bool  distant  = false;
  // arealight approximation
  vec3f         area_emission  = zero3f;
  frame3f       area_frame     = identity3x4f;
  frame3f       area_frend     = identity3x4f;
  vector<vec3i> area_triangles = {};
  vector<vec3f> area_positions = {};
  vector<vec3f> area_normals   = {};
};
struct pbrt_arealight {
  // arealight parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // arealight approximation
  vec3f emission = zero3f;
};
struct pbrt_environment {
  // shape parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // environment approximation
  vec3f  emission = zero3f;
  string filename = "";
};

// Other pbrt elements
struct pbrt_integrator {
  // integrator parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_film {
  // film parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  // film approximation
  string filename   = "";
  vec2i  resolution = zero2i;
};
struct pbrt_filter {
  // filter parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_accelerator {
  // accelerator parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_sampler {
  // sampler parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};

// Pbrt model
struct pbrt_model {
  vector<string>           comments     = {};
  vector<pbrt_camera>      cameras      = {};
  vector<pbrt_shape>       shapes       = {};
  vector<pbrt_texture>     textures     = {};
  vector<pbrt_material>    materials    = {};
  vector<pbrt_medium>      mediums      = {};
  vector<pbrt_environment> environments = {};
  vector<pbrt_arealight>   arealights   = {};
  vector<pbrt_light>       lights       = {};
  vector<pbrt_integrator>  integrators  = {};
  vector<pbrt_film>        films        = {};
  vector<pbrt_filter>      filters      = {};
  vector<pbrt_sampler>     samplers     = {};
  vector<pbrt_accelerator> accelerators = {};
};

// Load/save pbrt
inline void load_pbrt(const string& filename, pbrt_model& pbrt);
inline void save_pbrt(const string& filename, const pbrt_model& pbrt);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------
namespace yocto {

// A class that wraps a C file ti handle safe opening/closgin with RIIA.
struct pbrt_file {
  pbrt_file() {}
  pbrt_file(pbrt_file&& other);
  pbrt_file(const pbrt_file&) = delete;
  pbrt_file& operator=(const pbrt_file&) = delete;
  ~pbrt_file();

  operator bool() const { return (bool)fs; }

  FILE*  fs       = nullptr;
  string filename = "";
  string mode     = "rt";
  int    linenum  = 0;
};

// open a file
inline pbrt_file open_pbrt(const string& filename, const string& mode = "rt");
inline void      open_pbrt(
         pbrt_file& fs, const string& filename, const string& mode = "rt");
inline void close_pbrt(pbrt_file& fs);

// Pbrt command
enum struct pbrt_command {
  // clang-format off
  world_begin, world_end, attribute_begin, attribute_end,
  transform_begin, transform_end, reverse_orientation,
  set_transform, concat_transform, lookat_transform,
  object_instance, object_begin, object_end, include,      
  sampler, integrator, accelerator, film, filter, camera, shape, light,
  material, arealight, named_texture, named_medium, named_material,             
  use_material, medium_interface, active_transform,
  coordinate_system_set, coordinate_system_transform
  // clang-format on
};

// Read pbrt commands
inline bool read_pbrt_command(pbrt_file& fs, pbrt_command& command,
    string& name, string& type, frame3f& xform, vector<pbrt_value>& values);
inline bool read_pbrt_command(pbrt_file& fs, pbrt_command& command,
    string& name, string& type, frame3f& xform, vector<pbrt_value>& values,
    string& buffer);

// Write pbrt commands
inline void write_pbrt_comment(pbrt_file& fs, const string& comment);
inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_as_float = false);
inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name = "", const frame3f& xform = identity3x4f);
inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name, const string& type, const vector<pbrt_value>& values,
    bool texture_as_float = false);

// type-cheked pbrt value access
inline void get_pbrt_value(const pbrt_value& pbrt, string& value);
inline void get_pbrt_value(const pbrt_value& pbrt, bool& value);
inline void get_pbrt_value(const pbrt_value& pbrt, int& value);
inline void get_pbrt_value(const pbrt_value& pbrt, float& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vec2f& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vec3f& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vector<float>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vector<int>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value);
inline void get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value);
template <typename T>
inline void get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& value, T def) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      get_pbrt_value(p, value);
      return;
    }
  }
  value = def;
}
template <typename T>
inline T get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T def) {
  auto value = T{};
  get_pbrt_value(pbrt, name, value, def);
  return value;
}

// pbrt value construction
inline pbrt_value make_pbrt_value(const string& name, const string& value,
    pbrt_value_type type = pbrt_value_type::string);
inline pbrt_value make_pbrt_value(const string& name, bool value,
    pbrt_value_type type = pbrt_value_type::boolean);
inline pbrt_value make_pbrt_value(const string& name, int value,
    pbrt_value_type type = pbrt_value_type::integer);
inline pbrt_value make_pbrt_value(const string& name, float value,
    pbrt_value_type type = pbrt_value_type::real);
inline pbrt_value make_pbrt_value(const string& name, const vec2f& value,
    pbrt_value_type type = pbrt_value_type::point2);
inline pbrt_value make_pbrt_value(const string& name, const vec3f& value,
    pbrt_value_type type = pbrt_value_type::color);
inline pbrt_value make_pbrt_value(const string& name,
    const vector<vec2f>& value, pbrt_value_type type = pbrt_value_type::point2);
inline pbrt_value make_pbrt_value(const string& name,
    const vector<vec3f>& value, pbrt_value_type type = pbrt_value_type::point);
inline pbrt_value make_pbrt_value(const string& name,
    const vector<vec3i>&                        value,
    pbrt_value_type type = pbrt_value_type::integer);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#include <string_view>

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
inline string normalize_pbrt_path(const string& filename_) {
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
inline string get_pbrt_dirname(const string& filename_) {
  auto filename = normalize_pbrt_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
inline string get_pbrt_extension(const string& filename_) {
  auto filename = normalize_pbrt_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
inline string get_pbrt_filename(const string& filename_) {
  auto filename = normalize_pbrt_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
inline string get_pbrt_noextension(const string& filename_) {
  auto filename = normalize_pbrt_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
inline string get_pbrt_basename(const string& filename) {
  return get_pbrt_noextension(get_pbrt_filename(filename));
}

// Replaces extensions
inline string replace_pbrt_extension(
    const string& filename, const string& ext) {
  return get_pbrt_noextension(filename) + ext;
}

// Check if a file can be opened for reading.
inline bool exists_pbrt_file(const string& filename) {
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
// LOW-LEVEL FILE HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
inline pbrt_file ::pbrt_file(pbrt_file&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
inline pbrt_file ::~pbrt_file() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
inline void open_pbrt(
    pbrt_file& fs, const string& filename, const string& mode) {
  close_pbrt(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
inline pbrt_file open_pbrt(const string& filename, const string& mode) {
  auto fs = pbrt_file{};
  open_pbrt(fs, filename, mode);
  return fs;
}
inline void close_pbrt(pbrt_file& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

inline bool read_pbrt_line(pbrt_file& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

inline void write_ply_text(pbrt_file& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
inline void write_ply_text(pbrt_file& fs, const char* value) {
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
inline bool is_pbrt_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_pbrt_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_pbrt_whitespace(string_view& str) {
  while (!str.empty() && is_pbrt_space(str.front())) str.remove_prefix(1);
}

inline void remove_pbrt_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_pbrt_newline(str.back())) str.remove_suffix(1);
  auto cpy       = str;
  auto in_string = false;
  while (!cpy.empty()) {
    if (cpy.front() == '"') in_string = !in_string;
    if (cpy.front() == comment_char && !in_string) break;
    cpy.remove_prefix(1);
  }
  str.remove_suffix(cpy.size());
}

inline bool read_pbrt_line(obj_file& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

// Read a pbrt command from file
inline bool read_pbrt_cmdline(pbrt_file& fs, string& cmd, int& line_num) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs.fs);
  while (read_pbrt_line(fs, buffer, sizeof(buffer))) {
    // line
    line_num += 1;
    auto line = string_view{buffer};
    remove_pbrt_comment(line);
    skip_pbrt_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs.fs, pos, SEEK_SET);
        line_num -= 1;
        return true;
      } else {
        found = true;
      }
    } else if (!found) {
      throw std::runtime_error("bad pbrt command");
    }
    cmd += line;
    cmd += " ";
    pos = ftell(fs.fs);
  }
  return found;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL PRINTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_pbrt_value(string& str, const string& value) {
  str += value;
}
inline void format_pbrt_value(string& str, const char* value) { str += value; }
inline void format_pbrt_value(string& str, int value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_pbrt_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_pbrt_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) { if(i) str += " "; format_pbrt_value(str, value[i]); }
}
inline void format_pbrt_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) { if(i) str += " "; format_pbrt_value(str, value[i]); }
}
inline void format_pbrt_value(string& str, const vec4f& value) {
  for (auto i = 0; i < 4; i++) { if(i) str += " "; format_pbrt_value(str, value[i]); }
}
inline void format_pbrt_value(string& str, const mat4f& value) {
  for (auto i = 0; i < 4; i++) { if(i) str += " "; format_pbrt_value(str, value[i]); }
}

// Foramt to file
inline void format_pbrt_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_pbrt_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_pbrt_value(str, arg);
  format_pbrt_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
inline void format_pbrt_values(
    pbrt_file& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_pbrt_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
inline void format_pbrt_value(pbrt_file& fs, const T& value) {
  auto str = ""s;
  format_pbrt_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// parse a quoted string
inline void parse_pbrt_value(string_view& str, string_view& value) {
  skip_pbrt_whitespace(str);
  if (str.front() != '"') throw std::runtime_error("cannot parse value");
  str.remove_prefix(1);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != '"') cpy.remove_prefix(1);
  if (cpy.empty()) throw std::runtime_error("cannot parse value");
  value = str;
  value.remove_suffix(cpy.size());
  str.remove_prefix(str.size() - cpy.size());
  str.remove_prefix(1);
}

inline void parse_pbrt_value(string_view& str, string& value) {
  auto view = string_view{};
  parse_pbrt_value(str, view);
  value = string{view};
}

// parse a quoted string
inline void parse_pbrt_command(string_view& str, string& value) {
  skip_pbrt_whitespace(str);
  if (!isalpha((int)str.front())) {
    throw std::runtime_error("bad command");
  }
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == string_view::npos) {
    value.assign(str);
    str.remove_prefix(str.size());
  } else {
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
  }
}

// parse a number
inline void parse_pbrt_value(string_view& str, float& value) {
  skip_pbrt_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtof(str.data(), &next);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}

// parse a number
inline void parse_pbrt_value(string_view& str, int& value) {
  skip_pbrt_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtol(str.data(), &next, 10);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}
template <typename T>
inline void parse_pbrt_value(
    string_view& str, T& value, unordered_map<string, T>& value_names) {
  auto value_name = ""s;
  parse_pbrt_value(str, value_name);
  try {
    value = value_names.at(value_name);
  } catch (std::out_of_range&) {
    throw std::runtime_error("expected enum value");
  }
}

// parse a vec type
inline void parse_pbrt_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++) parse_pbrt_value(str, value[i]);
}
inline void parse_pbrt_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}
inline void parse_pbrt_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
inline void parse_pbrt_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}

// parse pbrt value with optional parens
template <typename T>
inline void parse_pbrt_param(string_view& str, T& value) {
  skip_pbrt_whitespace(str);
  auto parens = !str.empty() && str.front() == '[';
  if (parens) str.remove_prefix(1);
  parse_pbrt_value(str, value);
  if (parens) {
    skip_pbrt_whitespace(str);
    if (!str.empty() && str.front() == '[')
      throw std::runtime_error("bad pbrt param");
    str.remove_prefix(1);
  }
}

// parse a quoted string
inline void parse_pbrt_nametype(string_view& str_, string& name, string& type) {
  auto value = ""s;
  parse_pbrt_value(str_, value);
  auto str  = string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == string_view::npos) {
    throw std::runtime_error("bad type " + value);
  }
  type = string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == string_view::npos) {
    throw std::runtime_error("bad type " + value);
  }
  str.remove_prefix(pos2);
  name = string(str);
}

inline pair<vec3f, vec3f> get_pbrt_etak(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> metal_ior_table = {
      {"a-C", {{2.9440999183f, 2.2271502925f, 1.9681668794f},
                  {0.8874329109f, 0.7993216383f, 0.8152862927f}}},
      {"Ag", {{0.1552646489f, 0.1167232965f, 0.1383806959f},
                 {4.8283433224f, 3.1222459278f, 2.1469504455f}}},
      {"Al", {{1.6574599595f, 0.8803689579f, 0.5212287346f},
                 {9.2238691996f, 6.2695232477f, 4.8370012281f}}},
      {"AlAs", {{3.6051023902f, 3.2329365777f, 2.2175611545f},
                   {0.0006670247f, -0.0004999400f, 0.0074261204f}}},
      {"AlSb", {{-0.0485225705f, 4.1427547893f, 4.6697691348f},
                   {-0.0363741915f, 0.0937665154f, 1.3007390124f}}},
      {"Au", {{0.1431189557f, 0.3749570432f, 1.4424785571f},
                 {3.9831604247f, 2.3857207478f, 1.6032152899f}}},
      {"Be", {{4.1850592788f, 3.1850604423f, 2.7840913457f},
                 {3.8354398268f, 3.0101260162f, 2.8690088743f}}},
      {"Cr", {{4.3696828663f, 2.9167024892f, 1.6547005413f},
                 {5.2064337956f, 4.2313645277f, 3.7549467933f}}},
      {"CsI", {{2.1449030413f, 1.7023164587f, 1.6624194173f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Cu", {{0.2004376970f, 0.9240334304f, 1.1022119527f},
                 {3.9129485033f, 2.4528477015f, 2.1421879552f}}},
      {"Cu2O", {{3.5492833755f, 2.9520622449f, 2.7369202137f},
                   {0.1132179294f, 0.1946659670f, 0.6001681264f}}},
      {"CuO", {{3.2453822204f, 2.4496293965f, 2.1974114493f},
                  {0.5202739621f, 0.5707372756f, 0.7172250613f}}},
      {"d-C", {{2.7112524747f, 2.3185812849f, 2.2288565009f},
                  {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Hg", {{2.3989314904f, 1.4400254917f, 0.9095512090f},
                 {6.3276269444f, 4.3719414152f, 3.4217899270f}}},
      {"HgTe", {{4.7795267752f, 3.2309984581f, 2.6600252401f},
                   {1.6319827058f, 1.5808189339f, 1.7295753852f}}},
      {"Ir", {{3.0864098394f, 2.0821938440f, 1.6178866805f},
                 {5.5921510077f, 4.0671757150f, 3.2672611269f}}},
      {"K", {{0.0640493070f, 0.0464100621f, 0.0381842017f},
                {2.1042155920f, 1.3489364357f, 0.9132113889f}}},
      {"Li", {{0.2657871942f, 0.1956102432f, 0.2209198538f},
                 {3.5401743407f, 2.3111306542f, 1.6685930000f}}},
      {"MgO", {{2.0895885542f, 1.6507224525f, 1.5948759692f},
                  {0.0000000000f, -0.0000000000f, 0.0000000000f}}},
      {"Mo", {{4.4837010280f, 3.5254578255f, 2.7760769438f},
                 {4.1111307988f, 3.4208716252f, 3.1506031404f}}},
      {"Na", {{0.0602665320f, 0.0561412435f, 0.0619909494f},
                 {3.1792906496f, 2.1124800781f, 1.5790940266f}}},
      {"Nb", {{3.4201353595f, 2.7901921379f, 2.3955856658f},
                 {3.4413817900f, 2.7376437930f, 2.5799132708f}}},
      {"Ni", {{2.3672753521f, 1.6633583302f, 1.4670554172f},
                 {4.4988329911f, 3.0501643957f, 2.3454274399f}}},
      {"Rh", {{2.5857954933f, 1.8601866068f, 1.5544279524f},
                 {6.7822927110f, 4.7029501026f, 3.9760892461f}}},
      {"Se-e", {{5.7242724833f, 4.1653992967f, 4.0816099264f},
                   {0.8713747439f, 1.1052845009f, 1.5647788766f}}},
      {"Se", {{4.0592611085f, 2.8426947380f, 2.8207582835f},
                 {0.7543791750f, 0.6385150558f, 0.5215872029f}}},
      {"SiC", {{3.1723450205f, 2.5259677964f, 2.4793623897f},
                  {0.0000007284f, -0.0000006859f, 0.0000100150f}}},
      {"SnTe", {{4.5251865890f, 1.9811525984f, 1.2816819226f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"Ta", {{2.0625846607f, 2.3930915569f, 2.6280684948f},
                 {2.4080467973f, 1.7413705864f, 1.9470377016f}}},
      {"Te-e", {{7.5090397678f, 4.2964603080f, 2.3698732430f},
                   {5.5842076830f, 4.9476231084f, 3.9975145063f}}},
      {"Te", {{7.3908396088f, 4.4821028985f, 2.6370708478f},
                 {3.2561412892f, 3.5273908133f, 3.2921683116f}}},
      {"ThF4", {{1.8307187117f, 1.4422274283f, 1.3876488528f},
                   {0.0000000000f, 0.0000000000f, 0.0000000000f}}},
      {"TiC", {{3.7004673762f, 2.8374356509f, 2.5823030278f},
                  {3.2656905818f, 2.3515586388f, 2.1727857800f}}},
      {"TiN", {{1.6484691607f, 1.1504482522f, 1.3797795097f},
                  {3.3684596226f, 1.9434888540f, 1.1020123347f}}},
      {"TiO2-e", {{3.1065574823f, 2.5131551146f, 2.5823844157f},
                     {0.0000289537f, -0.0000251484f, 0.0001775555f}}},
      {"TiO2", {{3.4566203131f, 2.8017076558f, 2.9051485020f},
                   {0.0001026662f, -0.0000897534f, 0.0006356902f}}},
      {"VC", {{3.6575665991f, 2.7527298065f, 2.5326814570f},
                 {3.0683516659f, 2.1986687713f, 1.9631816252f}}},
      {"VN", {{2.8656011588f, 2.1191817791f, 1.9400767149f},
                 {3.0323264950f, 2.0561075580f, 1.6162930914f}}},
      {"V", {{4.2775126218f, 3.5131538236f, 2.7611257461f},
                {3.4911844504f, 2.8893580874f, 3.1116965117f}}},
      {"W", {{4.3707029924f, 3.3002972445f, 2.9982666528f},
                {3.5006778591f, 2.6048652781f, 2.2731930614f}}},
  };
  return metal_ior_table.at(name);
}

inline void parse_pbrt_params(string_view& str, vector<pbrt_value>& values) {
  auto parse_pbrt_pvalues = [](string_view& str, auto& value, auto& values) {
    values.clear();
    skip_pbrt_whitespace(str);
    if (str.empty()) throw std::runtime_error("bad pbrt value");
    if (str.front() == '[') {
      str.remove_prefix(1);
      skip_pbrt_whitespace(str);
      if (str.empty()) throw std::runtime_error("bad pbrt value");
      while (!str.empty()) {
        auto& val = values.empty() ? value : values.emplace_back();
        parse_pbrt_value(str, val);
        skip_pbrt_whitespace(str);
        if (str.empty()) break;
        if (str.front() == ']') break;
        if (values.empty()) values.push_back(value);
      }
      if (str.empty()) throw std::runtime_error("bad pbrt value");
      if (str.front() != ']') throw std::runtime_error("bad pbrt value");
      str.remove_prefix(1);
    } else {
      parse_pbrt_value(str, value);
    }
  };

  values.clear();
  skip_pbrt_whitespace(str);
  while (!str.empty()) {
    auto& value = values.emplace_back();
    auto  type  = ""s;
    parse_pbrt_nametype(str, value.name, type);
    skip_pbrt_whitespace(str);
    if (str.empty()) throw std::runtime_error("expected value");
    if (type == "float") {
      value.type = pbrt_value_type::real;
      parse_pbrt_pvalues(str, value.value1f, value.vector1f);
    } else if (type == "integer") {
      value.type = pbrt_value_type::integer;
      parse_pbrt_pvalues(str, value.value1i, value.vector1i);
    } else if (type == "string") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::string;
      parse_pbrt_pvalues(str, value.value1s, vector1s);
      if (!vector1s.empty())
        throw std::runtime_error("do not support pbrt string array");
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::boolean;
      parse_pbrt_pvalues(str, value1s, vector1s);
      if (!vector1s.empty())
        throw std::runtime_error("do not support pbrt string array");
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::texture;
      parse_pbrt_pvalues(str, value.value1s, vector1s);
      if (!vector1s.empty())
        throw std::runtime_error("do not support pbrt string array");
    } else if (type == "point" || type == "point3") {
      value.type = pbrt_value_type::point;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "normal" || type == "normal3") {
      value.type = pbrt_value_type::normal;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "vector" || type == "vector3") {
      value.type = pbrt_value_type::vector;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "point2") {
      value.type = pbrt_value_type::point2;
      parse_pbrt_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "vector2") {
      value.type = pbrt_value_type::vector2;
      parse_pbrt_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "blackbody") {
      value.type     = pbrt_value_type::color;
      auto blackbody = zero2f;
      auto vector2f  = vector<vec2f>{};
      parse_pbrt_pvalues(str, blackbody, vector2f);
      if (!vector2f.empty())
        throw std::runtime_error("bad pbrt " + type + " property");
      value.value3f = blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "color" || type == "rgb") {
      value.type = pbrt_value_type::color;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "xyz") {
      value.type = pbrt_value_type::color;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
      throw std::runtime_error("xyz conversion");
    } else if (type == "spectrum") {
      auto is_string = false;
      auto str1      = str;
      skip_pbrt_whitespace(str1);
      if (!str1.empty() && str1.front() == '"')
        is_string = true;
      else if (!str1.empty() && str1.front() == '[') {
        str1.remove_prefix(1);
        skip_pbrt_whitespace(str1);
        if (!str1.empty() && str1.front() == '"') is_string = true;
      }
      if (is_string) {
        value.type     = pbrt_value_type::color;
        auto filename  = ""s;
        auto filenames = vector<string>{};
        parse_pbrt_value(str, filename);
        auto filenamep = get_pbrt_filename(filename);
        if (get_pbrt_extension(filenamep) == ".spd") {
          filenamep = replace_pbrt_extension(filenamep, "");
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (get_pbrt_extension(filenamep) == ".eta") {
            auto eta =
                get_pbrt_etak(replace_pbrt_extension(filenamep, "")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (get_pbrt_extension(filenamep) == ".k") {
            auto k =
                get_pbrt_etak(replace_pbrt_extension(filenamep, "")).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            throw std::runtime_error("unknown spectrum file " + filename);
          }
        } else {
          throw std::runtime_error("unsupported spectrum format");
        }
      } else {
        value.type = pbrt_value_type::spectrum;
        parse_pbrt_pvalues(str, value.value1f, value.vector1f);
      }
    } else {
      throw std::runtime_error("unknown pbrt type");
    }
    skip_pbrt_whitespace(str);
  }
}

// convert pbrt films
inline void convert_pbrt_films(vector<pbrt_film>& films, bool verbose = false) {
  for (auto& film : films) {
    auto& values = film.values;
    if (film.type == "image") {
      film.resolution = {
          get_pbrt_value(values, "xresolution", 512),
          get_pbrt_value(values, "yresolution", 512),
      };
      film.filename = get_pbrt_value(values, "filename", "out.png"s);
    } else {
      throw std::runtime_error("unsupported Film type " + film.type);
    }
  }
}

// convert pbrt elements
inline void convert_pbrt_cameras(vector<pbrt_camera>& cameras,
    const vector<pbrt_film>& films, bool verbose = false) {
  auto film_aspect = 1.0f;
  for (auto& film : films) {
    film_aspect = (float)film.resolution.x / (float)film.resolution.y;
  }
  for (auto& camera : cameras) {
    auto& values   = camera.values;
    camera.frame   = inverse((frame3f)camera.frame);
    camera.frame.z = -camera.frame.z;
    if (camera.type == "perspective") {
      auto fov    = get_pbrt_value(values, "fov", 90.0f);
      camera.lens = 2 * tan(fov / 2) * 0.036;
      // auto lensradius = get_pbrt_value(values, "lensradius", 0.0f);
      camera.aspect = get_pbrt_value(values, "frameaspectratio", film_aspect);
      camera.focus  = get_pbrt_value(values, "focaldistance", 10.0f);
      if (!camera.aspect) camera.aspect = 1;
      if (!camera.focus) camera.focus = 10;
    } else if (camera.type == "realistic") {
      auto lensfile   = get_pbrt_value(values, "lensfile", ""s);
      lensfile        = lensfile.substr(0, lensfile.size() - 4);
      lensfile        = lensfile.substr(lensfile.find('.') + 1);
      lensfile        = lensfile.substr(0, lensfile.size() - 2);
      auto lens       = max(std::atof(lensfile.c_str()), 35.0f) * 0.001f;
      camera.lens     = 2 * atan(0.036f / (2 * lens));
      camera.aperture = get_pbrt_value(values, "aperturediameter", 0.0f);
      camera.focus    = get_pbrt_value(values, "focusdistance", 10.0f);
      camera.aspect   = film_aspect;
    } else {
      throw std::runtime_error("unsupported Camera type " + camera.type);
    }
  }
}

// convert pbrt textures
inline void convert_pbrt_textures(
    vector<pbrt_texture>& textures, bool verbose = false) {
  auto texture_map = unordered_map<string, int>{};
  for (auto& texture : textures) {
    auto index                = (int)texture_map.size();
    texture_map[texture.name] = index;
  }
  auto get_filename = [&textures, &texture_map](const string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return textures[pos->second].filename;
  };
  auto make_placeholder = [verbose](pbrt_texture& texture,
                              const vec3f&        color = {1, 0, 0}) {
    texture.constant = color;
    if (verbose)
      printf("texture %s not supported well\n", texture.type.c_str());
  };

  for (auto& texture : textures) {
    auto& values = texture.values;
    if (texture.type == "imagemap") {
      texture.filename = get_pbrt_value(values, "filename", ""s);
    } else if (texture.type == "constant") {
      texture.constant = get_pbrt_value(values, "value", vec3f{1});
    } else if (texture.type == "bilerp") {
      make_placeholder(texture, {1, 0, 0});
    } else if (texture.type == "checkerboard") {
      // auto tex1     = get_pbrt_value(values, "tex1", pair{vec3f{1}, ""s});
      // auto tex2     = get_pbrt_value(values, "tex2", pair{vec3f{0}, ""s});
      // auto rgb1     = tex1.second == "" ? tex1.first : vec3f{0.4f, 0.4f,
      // 0.4f}; auto rgb2     = tex1.second == "" ? tex2.first : vec3f{0.6f,
      // 0.6f, 0.6f}; auto params   = proc_image_params{}; params.type   =
      // proc_image_params::type_t::checker; params.color0 = {rgb1.x, rgb1.y,
      // rgb1.z, 1}; params.color1 = {rgb2.x, rgb2.y, rgb2.z, 1}; params.scale
      // = 2; make_proc_image(texture.hdr, params); float_to_byte(texture.ldr,
      // texture.hdr); texture.hdr = {};
      make_placeholder(texture, vec3f{0.5});
    } else if (texture.type == "dots") {
      make_placeholder(texture, vec3f{0.5});
    } else if (texture.type == "fbm") {
      make_placeholder(texture, vec3f{0.5});
    } else if (texture.type == "marble") {
      make_placeholder(texture, vec3f{0.5});
    } else if (texture.type == "mix") {
      auto tex1 = get_pbrt_value(values, "tex1", pair{vec3f{0}, ""s});
      auto tex2 = get_pbrt_value(values, "tex2", pair{vec3f{1}, ""s});
      if (!get_filename(tex1.second).empty()) {
        texture.filename = get_filename(tex1.second);
      } else if (!get_filename(tex2.second).empty()) {
        texture.filename = get_filename(tex2.second);
      } else {
        make_placeholder(texture);
      }
    } else if (texture.type == "scale") {
      auto tex1 = get_pbrt_value(values, "tex1", pair{vec3f{1}, ""s});
      auto tex2 = get_pbrt_value(values, "tex2", pair{vec3f{1}, ""s});
      if (!get_filename(tex1.second).empty()) {
        texture.filename = get_filename(tex1.second);
      } else if (!get_filename(tex2.second).empty()) {
        texture.filename = get_filename(tex2.second);
      } else {
        make_placeholder(texture);
      }
    } else if (texture.type == "uv") {
      make_placeholder(texture);
    } else if (texture.type == "windy") {
      make_placeholder(texture);
    } else if (texture.type == "wrinkled") {
      make_placeholder(texture);
    } else {
      throw std::runtime_error("unsupported texture type " + texture.type);
    }
  }
}

// convert pbrt materials
inline void convert_pbrt_materials(vector<pbrt_material>& materials,
    const vector<pbrt_texture>& textures, bool verbose = false) {
  // add constant textures
  auto constants = unordered_map<string, vec3f>{};
  for (auto& texture : textures) {
    if (!texture.filename.empty()) continue;
    constants[texture.name] = texture.constant;
  }

  // helpers
  auto get_scaled_texture = [&](const vector<pbrt_value>& values,
                                const string& name, vec3f& color,
                                string& texture, const vec3f& def) {
    auto textured = get_pbrt_value(values, name, pair{def, ""s});
    if (textured.second == "") {
      color   = textured.first;
      texture = "";
    } else if (constants.find(textured.second) != constants.end()) {
      color   = constants.at(textured.second);
      texture = "";
    } else {
      color   = {1, 1, 1};
      texture = textured.second;
    }
  };

  auto get_pbrt_roughness = [&](const vector<pbrt_value>& values,
                                vec2f& roughness, float def = 0.1) {
    auto roughness_ = get_pbrt_value(
        values, "roughness", pair{vec3f{def}, ""s});
    auto uroughness     = get_pbrt_value(values, "uroughness", roughness_);
    auto vroughness     = get_pbrt_value(values, "vroughness", roughness_);
    auto remaproughness = get_pbrt_value(values, "remaproughness", true);

    roughness = zero2f;
    if (uroughness.first == zero3f || vroughness.first == zero3f) return;
    roughness = vec2f{mean(uroughness.first), mean(vroughness.first)};
    // from pbrt code
    if (remaproughness) {
      roughness = max(roughness, 1e-3f);
      auto x    = log(roughness);
      roughness = 1.62142f + 0.819955f * x + 0.1734f * x * x +
                  0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }
    roughness = sqrt(roughness);
  };

  auto eta_to_reflectivity = [](const vec3f&  eta,
                                 const vec3f& etak = zero3f) -> vec3f {
    return ((eta - 1) * (eta - 1) + etak * etak) /
           ((eta + 1) * (eta + 1) + etak * etak);
  };

  // convert materials
  for (auto& material : materials) {
    auto& values = material.values;
    if (material.type == "uber") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.25});
      get_scaled_texture(
          values, "Ks", material.specular, material.specular_map, vec3f{0.25});
      get_scaled_texture(values, "Kt", material.transmission,
          material.transmission_map, vec3f{0});
      get_scaled_texture(
          values, "opacity", material.opacity, material.opacity_map, vec3f{1});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0.1f);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
    } else if (material.type == "plastic") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.25});
      get_scaled_texture(
          values, "Ks", material.specular, material.specular_map, vec3f{0.25});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0.1);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
    } else if (material.type == "translucent") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.25});
      get_scaled_texture(
          values, "Ks", material.specular, material.specular_map, vec3f{0.25});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0.1);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
    } else if (material.type == "matte") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.5});
      material.roughness = vec2f{1};
    } else if (material.type == "mirror") {
      get_scaled_texture(
          values, "Kr", material.specular, material.specular_map, vec3f{0.9});
      material.eta       = zero3f;
      material.etak      = zero3f;
      material.roughness = zero2f;
      material.sspecular = material.specular;
    } else if (material.type == "metal") {
      get_scaled_texture(
          values, "Kr", material.specular, material.specular_map, vec3f{1});
      get_scaled_texture(values, "eta", material.eta, material.eta_map,
          vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f});
      get_scaled_texture(values, "k", material.etak, material.etak_map,
          vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f});
      get_pbrt_roughness(values, material.roughness, 0.01);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta, material.etak);
    } else if (material.type == "substrate") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.5});
      get_scaled_texture(
          values, "Ks", material.specular, material.specular_map, vec3f{0.5});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0.1);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
    } else if (material.type == "glass") {
      get_scaled_texture(
          values, "Kr", material.specular, material.specular_map, vec3f{1});
      get_scaled_texture(values, "Kt", material.transmission,
          material.transmission_map, vec3f{1});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
      material.refract = true;
    } else if (material.type == "hair") {
      get_scaled_texture(
          values, "color", material.diffuse, material.diffuse_map, vec3f{0});
      material.roughness = {1, 1};
      if (verbose) printf("hair material not properly supported\n");
    } else if (material.type == "disney") {
      get_scaled_texture(
          values, "color", material.diffuse, material.diffuse_map, vec3f{0.5});
      material.roughness = {1, 1};
      if (verbose) printf("disney material not properly supported\n");
    } else if (material.type == "kdsubsurface") {
      get_scaled_texture(
          values, "Kd", material.diffuse, material.diffuse_map, vec3f{0.5});
      get_scaled_texture(
          values, "Kr", material.specular, material.specular_map, vec3f{1});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
      if (verbose) printf("kdsubsurface material not properly supported\n");
    } else if (material.type == "subsurface") {
      get_scaled_texture(
          values, "Kr", material.specular, material.specular_map, vec3f{1});
      get_scaled_texture(values, "Kt", material.transmission,
          material.transmission_map, vec3f{1});
      get_scaled_texture(
          values, "eta", material.eta, material.eta_map, vec3f{1.5});
      get_pbrt_roughness(values, material.roughness, 0);
      material.sspecular = material.specular *
                           eta_to_reflectivity(material.eta);
      auto scale        = get_pbrt_value(values, "scale", 1.0f);
      material.volscale = 1 / scale;
      auto sigma_a = zero3f, sigma_s = zero3f;
      auto sigma_a_tex = ""s, sigma_s_tex = ""s;
      get_scaled_texture(
          values, "sigma_a", sigma_a, sigma_a_tex, vec3f{0011, .0024, .014});
      get_scaled_texture(values, "sigma_prime_s", sigma_s, sigma_s_tex,
          vec3f{2.55, 3.12, 3.77});
      material.volmeanfreepath = 1 / (sigma_a + sigma_s);
      material.volscatter      = sigma_s / (sigma_a + sigma_s);
      if (verbose) printf("subsurface material not properly supported\n");
    } else if (material.type == "mix") {
      auto namedmaterial1 = get_pbrt_value(values, "namedmaterial1", ""s);
      auto namedmaterial2 = get_pbrt_value(values, "namedmaterial2", ""s);
      auto matname        = (!namedmaterial1.empty()) ? namedmaterial1
                                               : namedmaterial2;
      auto matit = std::find_if(materials.begin(), materials.end(),
          [&matname](auto& material) { return material.name == matname; });
      if (matit == materials.end())
        throw std::runtime_error("cannot find material " + matname);
      auto saved_name = material.name;
      material        = *matit;
      material.name   = saved_name;
      if (verbose) printf("mix material not properly supported\n");
    } else if (material.type == "fourier") {
      auto bsdffile = get_pbrt_value(values, "bsdffile", ""s);
      if (bsdffile.rfind("/") != string::npos)
        bsdffile = bsdffile.substr(bsdffile.rfind("/") + 1);
      if (bsdffile == "paint.bsdf") {
        material.diffuse   = {0.6f, 0.6f, 0.6f};
        material.specular  = {1, 1, 1};
        material.eta       = vec3f{1.5};
        material.roughness = vec2f{0.2};
        // material.roughness = get_pbrt_roughnessf(0.2f, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta);
      } else if (bsdffile == "ceramic.bsdf") {
        material.diffuse   = {0.6f, 0.6f, 0.6f};
        material.specular  = {1, 1, 1};
        material.eta       = vec3f{1.5};
        material.roughness = vec2f{0.25};
        // material.roughness = get_pbrt_roughnessf(0.25, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta);
      } else if (bsdffile == "leather.bsdf") {
        material.diffuse   = {0.6f, 0.57f, 0.48f};
        material.specular  = {1, 1, 1};
        material.eta       = vec3f{1.5};
        material.roughness = vec2f{0.3};
        // material.roughness = get_pbrt_roughnessf(0.3, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta);
      } else if (bsdffile == "coated_copper.bsdf") {
        material.specular  = vec3f{1};
        material.eta       = vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f};
        material.etak      = vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f};
        material.roughness = vec2f{0.01};
        // material.roughness = get_pbrt_roughnessf(0.01, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta, material.etak);
      } else if (bsdffile == "roughglass_alpha_0.2.bsdf") {
        material.specular     = {1, 1, 1};
        material.eta          = vec3f{1.5};
        material.transmission = {1, 1, 1};
        material.roughness    = vec2f{0.2};
        // material.roughness = get_pbrt_roughness(0.2, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta);
      } else if (bsdffile == "roughgold_alpha_0.2.bsdf") {
        material.specular  = vec3f{1, 1, 1};
        material.eta       = vec3f{0.1431189557f, 0.3749570432f, 1.4424785571f};
        material.etak      = vec3f{3.9831604247f, 2.3857207478f, 1.6032152899f};
        material.roughness = vec2f{0.2};
        // material.roughness = get_pbrt_roughness(0.2, true);
        material.sspecular = material.specular *
                             eta_to_reflectivity(material.eta, material.etak);
      } else {
        throw std::runtime_error("unsupported bsdffile " + bsdffile);
      }
    } else {
      throw std::runtime_error("unsupported material type" + material.type);
    }
  }
}

// Make a triangle shape from a quad grid
template <typename PositionFunc, typename NormalFunc>
inline void make_pbrt_shape(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    const PositionFunc& position_func, const NormalFunc& normal_func) {
  auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };
  auto tid = [steps](int i, int j, int c) { return (j * steps.x + i) * 2 + c; };
  positions.resize((steps.x + 1) * (steps.y + 1));
  normals.resize((steps.x + 1) * (steps.y + 1));
  texcoords.resize((steps.x + 1) * (steps.y + 1));
  for (auto j = 0; j < steps.y + 1; j++) {
    for (auto i = 0; i < steps.x + 1; i++) {
      auto uv              = vec2f{i / (float)steps.x, j / (float)steps.y};
      positions[vid(i, j)] = position_func(uv);
      normals[vid(i, j)]   = normal_func(uv);
      texcoords[vid(i, j)] = uv;
    }
  }
  triangles.resize(steps.x * steps.y * 2);
  for (auto j = 0; j < steps.y; j++) {
    for (auto i = 0; i < steps.x; i++) {
      triangles[tid(i, j, 0)] = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
      triangles[tid(i, j, 1)] = {vid(i, j), vid(i + 1, j + 1), vid(i, j + 1)};
    }
  }
}

// pbrt sphere
inline void make_pbrt_sphere(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_pbrt_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return radius *
               vec3f{cos(pt.x) * sin(pt.y), sin(pt.x) * sin(pt.y), cos(pt.y)};
      },
      [](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return vec3f{cos(pt.x) * cos(pt.y), sin(pt.x) * cos(pt.y), sin(pt.y)};
      });
}
inline void make_pbrt_disk(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_pbrt_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto a = 2 * pif * uv.x;
        return radius * (1 - uv.y) * vec3f{cos(a), sin(a), 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}
inline void make_pbrt_quad(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_pbrt_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        return vec3f{(uv.x - 0.5f) * radius, (uv.y - 0.5f) * radius, 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}

// Convert pbrt shapes
inline void convert_pbrt_shapes(
    vector<pbrt_shape>& shapes, bool verbose = false) {
  for (auto& shape : shapes) {
    auto& values = shape.values;
    if (shape.type == "trianglemesh") {
      get_pbrt_value(values, "P", shape.positions, {});
      get_pbrt_value(values, "N", shape.normals, {});
      get_pbrt_value(values, "uv", shape.texcoords, {});
      for (auto& uv : shape.texcoords) uv.y = (1 - uv.y);
      get_pbrt_value(values, "indices", shape.triangles, {});
    } else if (shape.type == "loopsubdiv") {
      get_pbrt_value(values, "P", shape.positions, {});
      get_pbrt_value(values, "indices", shape.triangles, {});
      shape.normals.resize(shape.positions.size());
      // compute_normals(shape.normals, shape.triangles, shape.positions);
    } else if (shape.type == "plymesh") {
      shape.filename = get_pbrt_value(values, "filename", ""s);
    } else if (shape.type == "sphere") {
      auto radius = get_pbrt_value(values, "radius", 1.0f);
      make_pbrt_sphere(shape.triangles, shape.positions, shape.normals,
          shape.texcoords, {32, 16}, radius);
    } else if (shape.type == "disk") {
      auto radius = get_pbrt_value(values, "radius", 1.0f);
      make_pbrt_disk(shape.triangles, shape.positions, shape.normals,
          shape.texcoords, {32, 1}, radius);
    } else {
      throw std::runtime_error("unsupported shape type " + shape.type);
    }
  }
}

// Convert pbrt arealights
inline void convert_pbrt_arealights(
    vector<pbrt_arealight>& lights, bool verbose = false) {
  for (auto& light : lights) {
    auto& values = light.values;
    if (light.type == "diffuse") {
      light.emission = get_pbrt_value(values, "L", vec3f{1, 1, 1}) *
                       get_pbrt_value(values, "scale", vec3f{1, 1, 1});
    } else {
      throw std::runtime_error("unsupported arealight type " + light.type);
    }
  }
}

// Convert pbrt lights
inline void convert_pbrt_lights(
    vector<pbrt_light>& lights, bool verbose = false) {
  for (auto& light : lights) {
    auto& values = light.values;
    if (light.type == "distant") {
      light.emission = get_pbrt_value(values, "scale", vec3f{1, 1, 1}) *
                       get_pbrt_value(values, "L", vec3f{1, 1, 1});
      light.from          = get_pbrt_value(values, "from", vec3f{0, 0, 0});
      light.to            = get_pbrt_value(values, "to", vec3f{0, 0, 1});
      light.distant       = true;
      auto distant_dist   = 100;
      auto size           = distant_dist * sin(5 * pif / 180);
      light.area_emission = light.emission * (distant_dist * distant_dist) /
                            (size * size);
      light.area_frame =
          light.frame *
          lookat_frame(normalize(light.from - light.to) * distant_dist, zero3f,
              {0, 1, 0}, true);
      light.area_frend =
          light.frend *
          lookat_frame(normalize(light.from - light.to) * distant_dist, zero3f,
              {0, 1, 0}, true);
      auto texcoords = vector<vec2f>{};
      make_pbrt_quad(light.area_triangles, light.area_positions,
          light.area_normals, texcoords, {4, 2}, size);
    } else if (light.type == "point" || light.type == "goniometric" ||
               light.type == "spot") {
      light.emission = get_pbrt_value(values, "scale", vec3f{1, 1, 1}) *
                       get_pbrt_value(values, "I", vec3f{1, 1, 1});
      light.from          = get_pbrt_value(values, "from", vec3f{0, 0, 0});
      light.area_emission = light.emission;
      light.area_frame    = light.frame * translation_frame(light.from);
      light.area_frend    = light.frend * translation_frame(light.from);
      auto texcoords      = vector<vec2f>{};
      make_pbrt_sphere(light.area_triangles, light.area_positions,
          light.area_normals, texcoords, {4, 2}, 0.0025f);
    } else {
      throw std::runtime_error("unsupported light type " + light.type);
    }
  }
}

inline void convert_pbrt_environments(vector<pbrt_environment>& environments,
    vector<pbrt_texture>& textures, bool verbose = false) {
  for (auto& light : environments) {
    auto& values = light.values;
    if (light.type == "infinite") {
      light.emission = get_pbrt_value(values, "scale", vec3f{1, 1, 1}) *
                       get_pbrt_value(values, "L", vec3f{1, 1, 1});
      light.filename = get_pbrt_value(values, "mapname", ""s);
      // environment.frame =
      // frame3f{{1,0,0},{0,0,-1},{0,-1,0},{0,0,0}}
      // * stack.back().frame;
      light.frame = light.frame *
                    frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
      light.frend = light.frend *
                    frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
    } else {
      throw std::runtime_error("unsupported environment type " + light.type);
    }
  }
}

// pbrt stack ctm
struct pbrt_context {
  frame3f transform_start        = identity3x4f;
  frame3f transform_end          = identity3x4f;
  string  material               = "";
  string  arealight              = "";
  string  medium_interior        = "";
  string  medium_exterior        = "";
  bool    reverse                = false;
  bool    active_transform_start = true;
  bool    active_transform_end   = true;
};

// load pbrt
inline void load_pbrt(const string& filename, pbrt_model& pbrt) {
  auto files = vector<pbrt_file>{};
  open_pbrt(files.emplace_back(), filename);
  if(!files.back()) throw std::runtime_error("cannot open " + filename);

  // parser state
  auto   stack      = vector<pbrt_context>{};
  string cur_object = "";

  // objects and coords
  unordered_map<string, pbrt_context> coordsys = {};
  unordered_map<string, vector<int>>  objects  = {};

  // helpers
  auto set_transform = [](pbrt_context& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](pbrt_context& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= xform;
    if (ctx.active_transform_end) ctx.transform_end *= xform;
  };

  // init stack
  if (stack.empty()) stack.emplace_back();

  // parse command by command
  while (!files.empty()) {
    auto line     = ""s;
    auto line_num = 0;
    while (read_pbrt_cmdline(files.back(), line, line_num)) {
      auto str = string_view{line};
      // get command
      auto cmd = ""s;
      parse_pbrt_command(str, cmd);
      if (cmd == "WorldBegin") {
        stack.push_back({});
      } else if (cmd == "WorldEnd") {
        if (stack.empty()) throw std::runtime_error("bad pbrt stack");
        stack.pop_back();
        if (stack.size() != 1) throw std::runtime_error("bad stack");
      } else if (cmd == "AttributeBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "AttributeEnd") {
        if (stack.empty()) throw std::runtime_error("bad pbrt stack");
        stack.pop_back();
      } else if (cmd == "TransformBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "TransformEnd") {
        if (stack.empty()) throw std::runtime_error("bad pbrt stack");
        stack.pop_back();
      } else if (cmd == "ObjectBegin") {
        stack.push_back(stack.back());
        parse_pbrt_param(str, cur_object);
        objects[cur_object] = {};
      } else if (cmd == "ObjectEnd") {
        stack.pop_back();
        cur_object = "";
      } else if (cmd == "ObjectInstance") {
        auto object = ""s;
        parse_pbrt_param(str, object);
        if (objects.find(object) == objects.end())
          throw std::runtime_error("cannot find object " + object);
        for (auto shape_id : objects.at(object)) {
          auto& shape = pbrt.shapes[shape_id];
          shape.instance_frames.push_back(stack.back().transform_start);
          shape.instance_frends.push_back(stack.back().transform_end);
        }
      } else if (cmd == "ActiveTransform") {
        auto name = ""s;
        parse_pbrt_command(str, name);
        if (name == "StartTime") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = false;
        } else if (name == "EndTime") {
          stack.back().active_transform_start = false;
          stack.back().active_transform_end   = true;
        } else if (name == "All") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = true;
        } else {
          throw std::runtime_error("bad active transform");
        }
      } else if (cmd == "Transform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        set_transform(stack.back(), frame3f{xf});
      } else if (cmd == "ConcatTransform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        concat_transform(stack.back(), frame3f{xf});
      } else if (cmd == "Scale") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), scaling_frame(v));
      } else if (cmd == "Translate") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), translation_frame(v));
      } else if (cmd == "Rotate") {
        auto v = zero4f;
        parse_pbrt_param(str, v);
        concat_transform(
            stack.back(), rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
      } else if (cmd == "LookAt") {
        auto from = zero3f, to = zero3f, up = zero3f;
        parse_pbrt_param(str, from);
        parse_pbrt_param(str, to);
        parse_pbrt_param(str, up);
        auto frame = lookat_frame(from, to, up, true);
        concat_transform(stack.back(), inverse(frame));
      } else if (cmd == "ReverseOrientation") {
        stack.back().reverse = !stack.back().reverse;
      } else if (cmd == "CoordinateSystem") {
        auto name = ""s;
        parse_pbrt_param(str, name);
        coordsys[name].transform_start = stack.back().transform_start;
        coordsys[name].transform_end   = stack.back().transform_end;
      } else if (cmd == "CoordSysTransform") {
        auto name = ""s;
        parse_pbrt_param(str, name);
        if (coordsys.find(name) != coordsys.end()) {
          stack.back().transform_start = coordsys.at(name).transform_start;
          stack.back().transform_end   = coordsys.at(name).transform_end;
        }
      } else if (cmd == "Integrator") {
        auto& integrator = pbrt.integrators.emplace_back();
        parse_pbrt_param(str, integrator.type);
        parse_pbrt_params(str, integrator.values);
      } else if (cmd == "Sampler") {
        auto& sampler = pbrt.samplers.emplace_back();
        parse_pbrt_param(str, sampler.type);
        parse_pbrt_params(str, sampler.values);
      } else if (cmd == "PixelFilter") {
        auto& filter = pbrt.filters.emplace_back();
        parse_pbrt_param(str, filter.type);
        parse_pbrt_params(str, filter.values);
      } else if (cmd == "Film") {
        auto& film = pbrt.films.emplace_back();
        parse_pbrt_param(str, film.type);
        parse_pbrt_params(str, film.values);
      } else if (cmd == "Accelerator") {
        auto& accelerator = pbrt.accelerators.emplace_back();
        parse_pbrt_param(str, accelerator.type);
        parse_pbrt_params(str, accelerator.values);
      } else if (cmd == "Camera") {
        auto& camera = pbrt.cameras.emplace_back();
        parse_pbrt_param(str, camera.type);
        parse_pbrt_params(str, camera.values);
        camera.frame = stack.back().transform_start;
        camera.frend = stack.back().transform_end;
      } else if (cmd == "Texture") {
        auto& texture  = pbrt.textures.emplace_back();
        auto  comptype = ""s;
        parse_pbrt_param(str, texture.name);
        parse_pbrt_param(str, comptype);
        parse_pbrt_param(str, texture.type);
        parse_pbrt_params(str, texture.values);
      } else if (cmd == "Material") {
        static auto material_id = 0;
        auto&       material    = pbrt.materials.emplace_back();
        material.name           = "material_" + std::to_string(material_id++);
        parse_pbrt_param(str, material.type);
        parse_pbrt_params(str, material.values);
        if (material.type == "") {
          stack.back().material = "";
          pbrt.materials.pop_back();
        } else {
          stack.back().material = material.name;
        }
      } else if (cmd == "MakeNamedMaterial") {
        auto& material = pbrt.materials.emplace_back();
        parse_pbrt_param(str, material.name);
        parse_pbrt_params(str, material.values);
        material.type = "";
        for (auto& value : material.values)
          if (value.name == "type") material.type = value.value1s;
      } else if (cmd == "NamedMaterial") {
        parse_pbrt_param(str, stack.back().material);
      } else if (cmd == "Shape") {
        auto& shape = pbrt.shapes.emplace_back();
        parse_pbrt_param(str, shape.type);
        parse_pbrt_params(str, shape.values);
        shape.frame     = stack.back().transform_start;
        shape.frend     = stack.back().transform_end;
        shape.material  = stack.back().material;
        shape.arealight = stack.back().arealight;
        shape.interior  = stack.back().medium_interior;
        shape.exterior  = stack.back().medium_exterior;
        if (cur_object != "") {
          shape.is_instanced = true;
          objects[cur_object].push_back((int)pbrt.shapes.size() - 1);
        } else {
          shape.instance_frames.push_back(identity3x4f);
          shape.instance_frends.push_back(identity3x4f);
        }
      } else if (cmd == "AreaLightSource") {
        static auto arealight_id = 0;
        auto&       arealight    = pbrt.arealights.emplace_back();
        arealight.name = "arealight_" + std::to_string(arealight_id++);
        parse_pbrt_param(str, arealight.type);
        parse_pbrt_params(str, arealight.values);
        arealight.frame        = stack.back().transform_start;
        arealight.frend        = stack.back().transform_end;
        stack.back().arealight = arealight.name;
      } else if (cmd == "LightSource") {
        auto& light = pbrt.lights.emplace_back();
        parse_pbrt_param(str, light.type);
        parse_pbrt_params(str, light.values);
        light.frame = stack.back().transform_start;
        light.frend = stack.back().transform_end;
        if (light.type == "infinite") {
          auto& environment  = pbrt.environments.emplace_back();
          environment.type   = light.type;
          environment.values = light.values;
          environment.frame  = light.frame;
          environment.frend  = light.frend;
          pbrt.lights.pop_back();
        }
      } else if (cmd == "MakeNamedMedium") {
        auto& medium = pbrt.mediums.emplace_back();
        parse_pbrt_param(str, medium.name);
        parse_pbrt_params(str, medium.values);
        medium.type = "";
        for (auto& value : medium.values)
          if (value.name == "type") medium.type = value.value1s;
      } else if (cmd == "MediumInterface") {
        parse_pbrt_param(str, stack.back().medium_interior);
        parse_pbrt_param(str, stack.back().medium_exterior);
      } else if (cmd == "Include") {
        auto includename = ""s;
        parse_pbrt_param(str, includename);
        open_pbrt(
            files.emplace_back(), get_pbrt_dirname(filename) + includename);
        if(!files.back()) throw std::runtime_error("cannot open " + filename);
      } else {
        throw std::runtime_error("unknown command " + cmd);
      }
    }
    files.pop_back();
  }

  // convert objects
  convert_pbrt_films(pbrt.films);
  convert_pbrt_cameras(pbrt.cameras, pbrt.films);
  convert_pbrt_textures(pbrt.textures);
  convert_pbrt_materials(pbrt.materials, pbrt.textures);
  convert_pbrt_shapes(pbrt.shapes);
  convert_pbrt_lights(pbrt.lights);
  convert_pbrt_arealights(pbrt.arealights);
  convert_pbrt_environments(pbrt.environments, pbrt.textures);

  // remove_pbrt_materials(pbrt.materials, pbrt.shapes);
  // remove_pbrt_textures(pbrt.textures, pbrt.materials);
}

inline static void format_pbrt_value(string& str, const pbrt_value& value) {
  static auto type_labels = unordered_map<pbrt_value_type, string>{
      {pbrt_value_type::real, "float"},
      {pbrt_value_type::integer, "integer"},
      {pbrt_value_type::boolean, "bool"},
      {pbrt_value_type::string, "string"},
      {pbrt_value_type::point, "point"},
      {pbrt_value_type::normal, "normal"},
      {pbrt_value_type::vector, "vector"},
      {pbrt_value_type::texture, "texture"},
      {pbrt_value_type::color, "rgb"},
      {pbrt_value_type::point2, "point2"},
      {pbrt_value_type::vector2, "vector2"},
      {pbrt_value_type::spectrum, "spectrum"},
  };

  auto format_vector = [](string& str, auto& values) {
    str += "[ ";
    for (auto& value : values) {
      str += " ";
      format_pbrt_value(str, value);
    }
    str += " ]";
  };

  format_pbrt_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case pbrt_value_type::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_pbrt_value(str, value.value1f);
      }
      break;
    case pbrt_value_type::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_pbrt_value(str, value.value1i);
      }
      break;
    case pbrt_value_type::boolean:
      format_pbrt_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case pbrt_value_type::string:
    case pbrt_value_type::texture:
      format_pbrt_values(str, "\"{}\"", value.value1s);
      break;
    case pbrt_value_type::point:
    case pbrt_value_type::vector:
    case pbrt_value_type::normal:
    case pbrt_value_type::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_pbrt_values(str, "[ {} ]", value.value3f);
      }
      break;
    case pbrt_value_type::spectrum: format_vector(str, value.vector1f); break;
    case pbrt_value_type::point2:
    case pbrt_value_type::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_pbrt_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

inline static void format_pbrt_value(
    string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_pbrt_value(str, value);
  }
}

inline void save_pbrt(const string& filename, const pbrt_model& pbrt) {
  auto fs = open_pbrt(filename, "wt");

  // save comments
  format_pbrt_values(fs, "#\n");
  format_pbrt_values(fs, "# Written by Yocto/GL\n");
  format_pbrt_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_pbrt_values(fs, "#\n\n");
  for (auto& comment : pbrt.comments) {
    format_pbrt_values(fs, "# {}\n", comment);
  }
  format_pbrt_values(fs, "\n");

  for (auto& camera_ : pbrt.cameras) {
    auto camera = camera_;
    if (camera.type == "") {
      camera.type = "perspective";
      camera.values.push_back(make_pbrt_value(
          "fov", 2 * tan(0.036f / (2 * camera.lens)) * 180 / pif));
    }
    format_pbrt_values(fs, "LookAt {} {} {}\n", camera.frame.o,
        camera.frame.o - camera.frame.z, camera.frame.y);
    format_pbrt_values(fs, "Camera \"{}\" {}\n", camera.type, camera.values);
  }
  for (auto& film_ : pbrt.films) {
    auto film = film_;
    if (film.type == "") {
      film.type = "image";
      film.values.push_back(make_pbrt_value("xresolution", film.resolution.x));
      film.values.push_back(make_pbrt_value("yresolution", film.resolution.y));
      film.values.push_back(make_pbrt_value("filename", film.filename));
    }
    format_pbrt_values(fs, "Film \"{}\" {}\n", film.type, film.values);
  }
  for (auto& integrator_ : pbrt.integrators) {
    auto integrator = integrator_;
    format_pbrt_values(
        fs, "Integrator \"{}\" {}\n", integrator.type, integrator.values);
  }
  for (auto& sampler_ : pbrt.samplers) {
    auto sampler = sampler_;
    format_pbrt_values(fs, "Sampler \"{}\" {}\n", sampler.type, sampler.values);
  }
  for (auto& filter_ : pbrt.filters) {
    auto filter = filter_;
    format_pbrt_values(
        fs, "PixelFilter \"{}\" {}\n", filter.type, filter.values);
  }
  for (auto& accelerator_ : pbrt.accelerators) {
    auto accelerator = accelerator_;
    format_pbrt_values(
        fs, "Accelerator \"{}\" {}\n", accelerator.type, accelerator.values);
  }

  format_pbrt_values(fs, "\nWorldBegin\n\n");

  for (auto& texture_ : pbrt.textures) {
    auto texture = texture_;
    if (texture.type == "") {
      if (texture.filename.empty()) {
        texture.type = "constant";
        texture.values.push_back(make_pbrt_value("value", texture.constant));
      } else {
        texture.type = "imagemap";
        texture.values.push_back(make_pbrt_value("filename", texture.filename));
      }
    }
    format_pbrt_values(fs, "Texture \"{}\" \"color\" \"{}\" {}\n", texture.name,
        texture.type, texture.values);
  }

  auto reflectivity_to_eta = [](const vec3f& reflectivity) {
    return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
  };
  for (auto& material_ : pbrt.materials) {
    auto material = material_;
    if (material.type == "") {
      if (material.specular != zero3f && material.transmission != zero3f &&
          material.refract) {
        material.type = "glass";
        material.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
        material.values.push_back(
            make_pbrt_value("roughness", pow(mean(material.roughness), 2)));
        material.values.push_back(make_pbrt_value(
            "eta", mean(reflectivity_to_eta(material.specular))));
        material.values.push_back(make_pbrt_value("remaproughness", false));
      } else if (mean(material.specular) > 0.1f &&
                 material.transmission == zero3f) {
        material.type = "metal";
        material.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
        material.values.push_back(
            make_pbrt_value("roughness", pow(mean(material.roughness), 2)));
        material.values.push_back(
            make_pbrt_value("eta", reflectivity_to_eta(material.specular)));
        material.values.push_back(make_pbrt_value("remaproughness", false));
      } else {
        material.type = "uber";
        if (material.diffuse_map.empty()) {
          material.values.push_back(make_pbrt_value("Kd", material.diffuse));
        } else if (material.diffuse != zero3f) {
          material.values.push_back(make_pbrt_value(
              "Kd", material.diffuse_map, pbrt_value_type::texture));
        }
        if (material.specular != zero3f) {
          material.values.push_back(make_pbrt_value("Ks", vec3f{1, 1, 1}));
          material.values.push_back(
              make_pbrt_value("roughness", pow(mean(material.roughness), 2)));
          material.values.push_back(make_pbrt_value(
              "eta", mean(reflectivity_to_eta(material.specular))));
          material.values.push_back(make_pbrt_value("remaproughness", false));
        }
        if (material.transmission != zero3f) {
          material.values.push_back(
              make_pbrt_value("Kt", material.transmission));
        }
        if (!material.opacity_map.empty()) {
          material.values.push_back(make_pbrt_value(
              "opacity", material.opacity_map, pbrt_value_type::texture));
        } else if (material.opacity != vec3f{1}) {
          material.values.push_back(
              make_pbrt_value("opacity", material.opacity));
        }
      }
    }
    format_pbrt_values(fs,
        "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n", material.name,
        material.type, material.values);
  }

  for (auto& medium_ : pbrt.mediums) {
    auto medium = medium_;
    format_pbrt_values(fs, "MakeNamedMedium \"{}\" \"string type\" \"{}\" {}\n",
        medium.name, medium.type, medium.values);
  }

  for (auto& light_ : pbrt.lights) {
    auto light = light_;
    if (light.type == "") {
      if (light.distant) {
        light.type = "distance";
        light.values.push_back(make_pbrt_value("L", light.emission));
      } else {
        light.type = "point";
        light.values.push_back(make_pbrt_value("I", light.emission));
      }
    }
    format_pbrt_values(fs, "AttributeBegin\n");
    format_pbrt_values(fs, "Transform {}\n", (mat4f)light.frame);
    format_pbrt_values(fs, "LightSource \"{}\" {}\n", light.type, light.values);
    format_pbrt_values(fs, "AttributeEnd\n");
  }

  for (auto& environment_ : pbrt.environments) {
    auto environment = environment_;
    if (environment.type == "") {
      environment.type = "infinite";
      environment.values.push_back(make_pbrt_value("L", environment.emission));
      environment.values.push_back(
          make_pbrt_value("mapname", environment.filename));
    }
    format_pbrt_values(fs, "AttributeBegin\n");
    format_pbrt_values(fs, "Transform {}\n", (mat4f)environment.frame);
    format_pbrt_values(
        fs, "LightSource \"{}\" {}\n", environment.type, environment.values);
    format_pbrt_values(fs, "AttributeEnd\n");
  }

  auto arealights_map = unordered_map<string, string>{};
  for (auto& arealight_ : pbrt.arealights) {
    auto arealight = arealight_;
    if (arealight.type == "") {
      arealight.type = "diffuse";
      arealight.values.push_back(make_pbrt_value("L", arealight.emission));
    }
    format_pbrt_values(arealights_map[arealight.name],
        "AreaLightSource \"{}\" {}\n", arealight.type, arealight.values);
  }

  auto object_id = 0;
  for (auto& shape_ : pbrt.shapes) {
    auto shape = shape_;
    if (shape.type == "") {
      if (!shape.filename.empty()) {
        shape.type = "plymesh";
        shape.values.push_back(make_pbrt_value("filename", shape.filename));
      } else {
        shape.type = "trianglemesh";
        shape.values.push_back(make_pbrt_value("indices", shape.triangles));
        shape.values.push_back(
            make_pbrt_value("P", shape.positions, pbrt_value_type::point));
        if (!shape.normals.empty())
          shape.values.push_back(
              make_pbrt_value("N", shape.triangles, pbrt_value_type::normal));
        if (!shape.texcoords.empty())
          shape.values.push_back(make_pbrt_value("uv", shape.texcoords));
      }
    }
    auto object = "object" + std::to_string(object_id++);
    if (shape.is_instanced)
      format_pbrt_values(fs, "ObjectBegin \"{}\"\n", object);
    format_pbrt_values(fs, "AttributeBegin\n");
    format_pbrt_values(fs, "Transform {}\n", (mat4f)shape.frame);
    format_pbrt_values(fs, "NamedMaterial \"{}\"\n", shape.material);
    if (shape.arealight != "")
      format_pbrt_values(fs, arealights_map.at(shape.arealight));
    format_pbrt_values(fs, "Shape \"{}\" {}\n", shape.type, shape.values);
    format_pbrt_values(fs, "AttributeEnd\n");
    if (shape.is_instanced) format_pbrt_values(fs, "ObjectEnd\n");
    for (auto& iframe : shape.instance_frames) {
      format_pbrt_values(fs, "AttributeBegin\n");
      format_pbrt_values(fs, "Transform {}\n", (mat4f)iframe);
      format_pbrt_values(fs, "ObjectInstance \"{}\"\n", object);
      format_pbrt_values(fs, "AttributeEnd\n");
    }
  }

  format_pbrt_values(fs, "\nWorldEnd\n\n");
}

// Read pbrt commands
inline bool read_pbrt_command(pbrt_file& fs, pbrt_command& command,
    string& name, string& type, frame3f& xform, vector<pbrt_value>& values,
    string& line) {
  // parse command by command
  auto line_num = 0;
  while (read_pbrt_cmdline(fs, line, line_num)) {
    auto str = string_view{line};
    // get command
    auto cmd = ""s;
    parse_pbrt_command(str, cmd);
    if (cmd == "WorldBegin") {
      command = pbrt_command::world_begin;
      return true;
    } else if (cmd == "WorldEnd") {
      command = pbrt_command::world_end;
      return true;
    } else if (cmd == "AttributeBegin") {
      command = pbrt_command::attribute_begin;
      return true;
    } else if (cmd == "AttributeEnd") {
      command = pbrt_command::attribute_end;
      return true;
    } else if (cmd == "TransformBegin") {
      command = pbrt_command::transform_begin;
      return true;
    } else if (cmd == "TransformEnd") {
      command = pbrt_command::transform_end;
      return true;
    } else if (cmd == "ObjectBegin") {
      parse_pbrt_param(str, name);
      command = pbrt_command::object_begin;
      return true;
    } else if (cmd == "ObjectEnd") {
      command = pbrt_command::object_end;
      return true;
    } else if (cmd == "ObjectInstance") {
      parse_pbrt_param(str, name);
      command = pbrt_command::object_instance;
      return true;
    } else if (cmd == "ActiveTransform") {
      parse_pbrt_command(str, name);
      command = pbrt_command::active_transform;
      return true;
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      xform   = frame3f{xf};
      command = pbrt_command::set_transform;
      return true;
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      xform   = frame3f{xf};
      command = pbrt_command::concat_transform;
      return true;
    } else if (cmd == "Scale") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      xform   = scaling_frame(v);
      command = pbrt_command::concat_transform;
      return true;
    } else if (cmd == "Translate") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      xform   = translation_frame(v);
      command = pbrt_command::concat_transform;
      return true;
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      parse_pbrt_param(str, v);
      xform   = rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x));
      command = pbrt_command::concat_transform;
      return true;
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      parse_pbrt_param(str, from);
      parse_pbrt_param(str, to);
      parse_pbrt_param(str, up);
      xform   = {from, to, up, zero3f};
      command = pbrt_command::lookat_transform;
      return true;
    } else if (cmd == "ReverseOrientation") {
      command = pbrt_command::reverse_orientation;
      return true;
    } else if (cmd == "CoordinateSystem") {
      parse_pbrt_param(str, name);
      command = pbrt_command::coordinate_system_set;
      return true;
    } else if (cmd == "CoordSysTransform") {
      parse_pbrt_param(str, name);
      command = pbrt_command::coordinate_system_transform;
      return true;
    } else if (cmd == "Integrator") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::integrator;
      return true;
    } else if (cmd == "Sampler") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::sampler;
      return true;
    } else if (cmd == "PixelFilter") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::filter;
      return true;
    } else if (cmd == "Film") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::film;
      return true;
    } else if (cmd == "Accelerator") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::accelerator;
      return true;
    } else if (cmd == "Camera") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::camera;
      return true;
    } else if (cmd == "Texture") {
      auto comptype = ""s;
      parse_pbrt_param(str, name);
      parse_pbrt_param(str, comptype);
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::named_texture;
      return true;
    } else if (cmd == "Material") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::material;
      return true;
    } else if (cmd == "MakeNamedMaterial") {
      parse_pbrt_param(str, name);
      parse_pbrt_params(str, values);
      type = "";
      for (auto& value : values)
        if (value.name == "type") type = value.value1s;
      command = pbrt_command::named_material;
      return true;
    } else if (cmd == "NamedMaterial") {
      parse_pbrt_param(str, name);
      command = pbrt_command::use_material;
      return true;
    } else if (cmd == "Shape") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::shape;
      return true;
    } else if (cmd == "AreaLightSource") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::arealight;
      return true;
    } else if (cmd == "LightSource") {
      parse_pbrt_param(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command::light;
      return true;
    } else if (cmd == "MakeNamedMedium") {
      parse_pbrt_param(str, name);
      parse_pbrt_params(str, values);
      type = "";
      for (auto& value : values)
        if (value.name == "type") type = value.value1s;
      command = pbrt_command::named_medium;
      return true;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      parse_pbrt_param(str, interior);
      parse_pbrt_param(str, exterior);
      name    = interior + "####" + exterior;
      command = pbrt_command::medium_interface;
      return true;
    } else if (cmd == "Include") {
      parse_pbrt_param(str, name);
      command = pbrt_command::include;
      return true;
    } else {
      throw std::runtime_error("unknown command " + cmd);
    }
  }
  return false;
}
inline bool read_pbrt_command(pbrt_file& fs, pbrt_command& command,
    string& name, string& type, frame3f& xform, vector<pbrt_value>& values) {
  auto command_buffer = ""s;
  return read_pbrt_command(
      fs, command, name, type, xform, values, command_buffer);
}

inline vector<string> split_pbrt_string(
    const string& str, const string& delim) {
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
inline void write_pbrt_comment(pbrt_file& fs, const string& comment) {
  auto lines = split_pbrt_string(comment, "\n");
  for (auto& line : lines) {
    format_pbrt_values(fs, "# {}\n", line);
  }
  format_pbrt_values(fs, "\n");
}

inline void write_pbrt_values(pbrt_file& fs, const vector<pbrt_value>& values) {
  static auto type_labels = unordered_map<pbrt_value_type, string>{
      {pbrt_value_type::real, "float"},
      {pbrt_value_type::integer, "integer"},
      {pbrt_value_type::boolean, "bool"},
      {pbrt_value_type::string, "string"},
      {pbrt_value_type::point, "point"},
      {pbrt_value_type::normal, "normal"},
      {pbrt_value_type::vector, "vector"},
      {pbrt_value_type::texture, "texture"},
      {pbrt_value_type::color, "rgb"},
      {pbrt_value_type::point2, "point2"},
      {pbrt_value_type::vector2, "vector2"},
      {pbrt_value_type::spectrum, "spectrum"},
  };
  for (auto& value : values) {
    format_pbrt_values(
        fs, " \"{} {}\" ", type_labels.at(value.type), value.name);
    switch (value.type) {
      case pbrt_value_type::real:
        if (!value.vector1f.empty()) {
          format_pbrt_values(fs, "[ ");
          for (auto& v : value.vector1f) format_pbrt_values(fs, " {}", v);
          format_pbrt_values(fs, " ]");
        } else {
          format_pbrt_values(fs, "{}", value.value1f);
        }
        break;
      case pbrt_value_type::integer:
        if (!value.vector1f.empty()) {
          format_pbrt_values(fs, "[ ");
          for (auto& v : value.vector1i) format_pbrt_values(fs, " {}", v);
          format_pbrt_values(fs, " ]");
        } else {
          format_pbrt_values(fs, "{}", value.value1i);
        }
        break;
      case pbrt_value_type::boolean:
        format_pbrt_values(fs, "\"{}\"", value.value1b ? "true" : "false");
        break;
      case pbrt_value_type::string:
        format_pbrt_values(fs, "\"{}\"", value.value1s);
        break;
      case pbrt_value_type::point:
      case pbrt_value_type::vector:
      case pbrt_value_type::normal:
      case pbrt_value_type::color:
        if (!value.vector3f.empty()) {
          format_pbrt_values(fs, "[ ");
          for (auto& v : value.vector3f) format_pbrt_values(fs, " {}", v);
          format_pbrt_values(fs, " ]");
        } else {
          format_pbrt_values(fs, "[ {} ]", value.value3f);
        }
        break;
      case pbrt_value_type::spectrum:
        format_pbrt_values(fs, "[ ");
        for (auto& v : value.vector1f) format_pbrt_values(fs, " {}", v);
        format_pbrt_values(fs, " ]");
        break;
      case pbrt_value_type::texture:
        format_pbrt_values(fs, "\"{}\"", value.value1s);
        break;
      case pbrt_value_type::point2:
      case pbrt_value_type::vector2:
        if (!value.vector2f.empty()) {
          format_pbrt_values(fs, "[ ");
          for (auto& v : value.vector2f) format_pbrt_values(fs, " {}", v);
          format_pbrt_values(fs, " ]");
        } else {
          format_pbrt_values(fs, "[ {} ]", value.value2f);
        }
        break;
    }
  }
  format_pbrt_values(fs, "\n");
}

inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_float) {
  switch (command) {
    case pbrt_command::world_begin:
      format_pbrt_values(fs, "WorldBegin\n");
      break;
    case pbrt_command::world_end: format_pbrt_values(fs, "WorldEnd\n"); break;
    case pbrt_command::attribute_begin:
      format_pbrt_values(fs, "AttributeBegin\n");
      break;
    case pbrt_command::attribute_end:
      format_pbrt_values(fs, "AttributeEnd\n");
      break;
    case pbrt_command::transform_begin:
      format_pbrt_values(fs, "TransformBegin\n");
      break;
    case pbrt_command::transform_end:
      format_pbrt_values(fs, "TransformEnd\n");
      break;
    case pbrt_command::object_begin:
      format_pbrt_values(fs, "ObjectBegin \"{}\"\n", name);
      break;
    case pbrt_command::object_end: format_pbrt_values(fs, "ObjectEnd\n"); break;
    case pbrt_command::object_instance:
      format_pbrt_values(fs, "ObjectInstance \"{}\"\n", name);
      break;
    case pbrt_command::sampler:
      format_pbrt_values(fs, "Sampler \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::integrator:
      format_pbrt_values(fs, "Integrator \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::accelerator:
      format_pbrt_values(fs, "Accelerator \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::film:
      format_pbrt_values(fs, "Film \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::filter:
      format_pbrt_values(fs, "Filter \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::camera:
      format_pbrt_values(fs, "Camera \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::shape:
      format_pbrt_values(fs, "Shape \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::light:
      format_pbrt_values(fs, "LightSource \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::material:
      format_pbrt_values(fs, "Material \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::arealight:
      format_pbrt_values(fs, "AreaLightSource \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_texture:
      format_pbrt_values(fs, "Texture \"{}\" \"{}\" \"{}\"", name,
          texture_float ? "float" : "rgb", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_medium:
      format_pbrt_values(
          fs, "MakeNamedMedium \"{}\" \"string type\" \"{}\"", name, type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_material:
      format_pbrt_values(
          fs, "MakeNamedMaterial \"{}\" \"string type\" \"{}\"", name, type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::include:
      format_pbrt_values(fs, "Include \"{}\"\n", name);
      break;
    case pbrt_command::reverse_orientation:
      format_pbrt_values(fs, "ReverseOrientation\n");
      break;
    case pbrt_command::set_transform:
      format_pbrt_values(fs, "Transform {}\n", (mat4f)xform);
      break;
    case pbrt_command::concat_transform:
      format_pbrt_values(fs, "ConcatTransform {}\n", (mat4f)xform);
      break;
    case pbrt_command::lookat_transform:
      format_pbrt_values(fs, "LookAt {} {} {}\n", xform.x, xform.y, xform.z);
      break;
    case pbrt_command::use_material:
      format_pbrt_values(fs, "NamedMaterial \"{}\"\n", name);
      break;
    case pbrt_command::medium_interface: {
      auto interior = ""s, exterior = ""s;
      auto found = false;
      for (auto c : name) {
        if (c == '#') {
          found = true;
          continue;
        }
        if (found)
          exterior.push_back(c);
        else
          interior.push_back(c);
      }
      format_pbrt_values(
          fs, "MediumInterface \"{}\" \"{}\"\n", interior, exterior);
    } break;
    case pbrt_command::active_transform:
      format_pbrt_values(fs, "ActiveTransform \"{}\"\n", name);
      break;
    case pbrt_command::coordinate_system_set:
      format_pbrt_values(fs, "CoordinateSystem \"{}\"\n", name);
      break;
    case pbrt_command::coordinate_system_transform:
      format_pbrt_values(fs, "CoordinateSysTransform \"{}\"\n", name);
      break;
  }
}

inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name, const frame3f& xform) {
  return write_pbrt_command(fs, command, name, "", xform, {});
}
inline void write_pbrt_command(pbrt_file& fs, pbrt_command command,
    const string& name, const string& type, const vector<pbrt_value>& values,
    bool texture_as_float) {
  return write_pbrt_command(
      fs, command, name, type, identity3x4f, values, texture_as_float);
}

// get pbrt value
inline void get_pbrt_value(const pbrt_value& pbrt, string& value) {
  if (pbrt.type == pbrt_value_type::string ||
      pbrt.type == pbrt_value_type::texture) {
    value = pbrt.value1s;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, bool& value) {
  if (pbrt.type == pbrt_value_type::boolean) {
    value = pbrt.value1b;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, int& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    value = pbrt.value1i;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, float& value) {
  if (pbrt.type == pbrt_value_type::real) {
    value = pbrt.value1f;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vec2f& value) {
  if (pbrt.type == pbrt_value_type::point2 ||
      pbrt.type == pbrt_value_type::vector2) {
    value = pbrt.value2f;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vec3f& value) {
  if (pbrt.type == pbrt_value_type::point ||
      pbrt.type == pbrt_value_type::vector ||
      pbrt.type == pbrt_value_type::normal ||
      pbrt.type == pbrt_value_type::color) {
    value = pbrt.value3f;
  } else if (pbrt.type == pbrt_value_type::real) {
    value = vec3f{pbrt.value1f};
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vector<float>& value) {
  if (pbrt.type == pbrt_value_type::real) {
    if (!pbrt.vector1f.empty()) {
      value = pbrt.vector1f;
    } else {
      value = {pbrt.value1f};
    }
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& value) {
  if (pbrt.type == pbrt_value_type::point2 ||
      pbrt.type == pbrt_value_type::vector2) {
    if (!pbrt.vector2f.empty()) {
      value = pbrt.vector2f;
    } else {
      value = {pbrt.value2f};
    }
  } else if (pbrt.type == pbrt_value_type::real) {
    if (pbrt.vector1f.empty() || pbrt.vector1f.size() % 2)
      throw std::runtime_error("bad pbrt type");
    value.resize(pbrt.vector1f.size() / 2);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1f[i * 2 + 0], pbrt.vector1f[i * 2 + 1]};
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& value) {
  if (pbrt.type == pbrt_value_type::point ||
      pbrt.type == pbrt_value_type::vector ||
      pbrt.type == pbrt_value_type::normal ||
      pbrt.type == pbrt_value_type::color) {
    if (!pbrt.vector3f.empty()) {
      value = pbrt.vector3f;
    } else {
      value = {pbrt.value3f};
    }
  } else if (pbrt.type == pbrt_value_type::real) {
    if (pbrt.vector1f.empty() || pbrt.vector1f.size() % 3)
      throw std::runtime_error("bad pbrt type");
    value.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}

inline void get_pbrt_value(const pbrt_value& pbrt, vector<int>& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    if (!pbrt.vector1i.empty()) {
      value = pbrt.vector1i;
    } else {
      value = {pbrt.vector1i};
    }
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    if (pbrt.vector1i.empty() || pbrt.vector1i.size() % 3)
      throw std::runtime_error("bad pbrt type");
    value.resize(pbrt.vector1i.size() / 3);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1i[i * 3 + 0], pbrt.vector1i[i * 3 + 1],
          pbrt.vector1i[i * 3 + 2]};
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value) {
  if (pbrt.type == pbrt_value_type::string) {
    value.first = 0;
    get_pbrt_value(pbrt, value.second);
  } else {
    get_pbrt_value(pbrt, value.first);
    value.second = "";
  }
}
inline void get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value) {
  if (pbrt.type == pbrt_value_type::string ||
      pbrt.type == pbrt_value_type::texture) {
    value.first = zero3f;
    get_pbrt_value(pbrt, value.second);
  } else {
    get_pbrt_value(pbrt, value.first);
    value.second = "";
  }
}

// pbrt value construction
inline pbrt_value make_pbrt_value(
    const string& name, const string& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, bool value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, int value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, float value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec2f& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec3f& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vector<vec2f>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vector<vec3f>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = value;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vector<vec3i>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {(int*)value.data(), (int*)value.data() + value.size() * 3};
  return pbrt;
}

// old code --- maintained here in case we want to integrate back
#if 0
inline void approximate_fourier_material(pbrt_material::fourier_t& fourier) {
  auto filename = get_filename(fourier.bsdffile);
  if (filename == "paint.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.6f, 0.6f};
    // plastic.Ks = {0.4f, 0.4f, 0.4f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.2f;
    plastic.vroughness = 0.2f;
  } else if (filename == "ceramic.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.6f, 0.6f};
    // plastic.Ks = {0.1f, 0.1f, 0.1f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.025f;
    plastic.vroughness = 0.025f;
  } else if (filename == "leather.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::plastic;
    auto& plastic       = fourier.approx_plastic;
    plastic.Kd          = {0.6f, 0.57f, 0.48f};
    // plastic.Ks = {0.1f, 0.1f, 0.1f};
    plastic.Ks         = {1.0f, 1.0f, 1.0f};
    plastic.uroughness = 0.3f;
    plastic.vroughness = 0.3f;
  } else if (filename == "coated_copper.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::metal;
    auto& metal         = fourier.approx_metal;
    auto  etak          = get_pbrt_etak("Cu");
    metal.eta           = {etak.first.x, etak.first.y, etak.first.z};
    metal.k             = {etak.second.x, etak.second.y, etak.second.z};
    metal.uroughness    = 0.01f;
    metal.vroughness    = 0.01f;
  } else if (filename == "roughglass_alpha_0.2.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::glass;
    auto& glass         = fourier.approx_glass;
    glass.uroughness    = 0.2f;
    glass.vroughness    = 0.2f;
    glass.Kr            = {1, 1, 1};
    glass.Kt            = {1, 1, 1};
  } else if (filename == "roughgold_alpha_0.2.bsdf") {
    fourier.approx_type = pbrt_material::fourier_t::approx_type_t::metal;
    auto& metal         = fourier.approx_metal;
    auto  etak          = get_pbrt_etak("Au");
    metal.eta           = {etak.first.x, etak.first.y, etak.first.z};
    metal.k             = {etak.second.x, etak.second.y, etak.second.z};
    metal.uroughness    = 0.2f;
    metal.vroughness    = 0.2f;
  } else {
    throw std::runtime_error("unknown pbrt bsdf filename " + fourier.bsdffile);
  }
}

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
inline pair<vec3f, vec3f> parse_pbrt_subsurface(const string& name) {
  static const unordered_map<string, pair<vec3f, vec3f>> params = {
      // From "A Practical Model for Subsurface Light Transport"
      // Jensen, Marschner, Levoy, Hanrahan
      // Proc SIGGRAPH 2001
      {"Apple", {{2.29, 2.39, 1.97}, {0.0030, 0.0034, 0.046}}},
      {"Chicken1", {{0.15, 0.21, 0.38}, {0.015, 0.077, 0.19}}},
      {"Chicken2", {{0.19, 0.25, 0.32}, {0.018, 0.088, 0.20}}},
      {"Cream", {{7.38, 5.47, 3.15}, {0.0002, 0.0028, 0.0163}}},
      {"Ketchup", {{0.18, 0.07, 0.03}, {0.061, 0.97, 1.45}}},
      {"Marble", {{2.19, 2.62, 3.00}, {0.0021, 0.0041, 0.0071}}},
      {"Potato", {{0.68, 0.70, 0.55}, {0.0024, 0.0090, 0.12}}},
      {"Skimmilk", {{0.70, 1.22, 1.90}, {0.0014, 0.0025, 0.0142}}},
      {"Skin1", {{0.74, 0.88, 1.01}, {0.032, 0.17, 0.48}}},
      {"Skin2", {{1.09, 1.59, 1.79}, {0.013, 0.070, 0.145}}},
      {"Spectralon", {{11.6, 20.4, 14.9}, {0.00, 0.00, 0.00}}},
      {"Wholemilk", {{2.55, 3.21, 3.77}, {0.0011, 0.0024, 0.014}}},
      // From "Acquiring Scattering Properties of Participating Media by
      // Dilution",
      // Narasimhan, Gupta, Donner, Ramamoorthi, Nayar, Jensen
      // Proc SIGGRAPH 2006
      {"Lowfat Milk", {{0.89187, 1.5136, 2.532}, {0.002875, 0.00575, 0.0115}}},
      {"Reduced Milk",
          {{2.4858, 3.1669, 4.5214}, {0.0025556, 0.0051111, 0.012778}}},
      {"Regular Milk",
          {{4.5513, 5.8294, 7.136}, {0.0015333, 0.0046, 0.019933}}},
      {"Espresso", {{0.72378, 0.84557, 1.0247}, {4.7984, 6.5751, 8.8493}}},
      {"Mint Mocha Coffee",
          {{0.31602, 0.38538, 0.48131}, {3.772, 5.8228, 7.82}}},
      {"Lowfat Soy Milk",
          {{0.30576, 0.34233, 0.61664}, {0.0014375, 0.0071875, 0.035937}}},
      {"Regular Soy Milk",
          {{0.59223, 0.73866, 1.4693}, {0.0019167, 0.0095833, 0.065167}}},
      {"Lowfat Chocolate Milk",
          {{0.64925, 0.83916, 1.1057}, {0.0115, 0.0368, 0.1564}}},
      {"Regular Chocolate Milk",
          {{1.4585, 2.1289, 2.9527}, {0.010063, 0.043125, 0.14375}}},
      {"Coke", {{8.9053e-05, 8.372e-05, 0}, {0.10014, 0.16503, 0.2468}}},
      {"Pepsi", {{6.1697e-05, 4.2564e-05, 0}, {0.091641, 0.14158, 0.20729}}},
      {"Sprite", {{6.0306e-06, 6.4139e-06, 6.5504e-06},
                     {0.001886, 0.0018308, 0.0020025}}},
      {"Gatorade",
          {{0.0024574, 0.003007, 0.0037325}, {0.024794, 0.019289, 0.008878}}},
      {"Chardonnay", {{1.7982e-05, 1.3758e-05, 1.2023e-05},
                         {0.010782, 0.011855, 0.023997}}},
      {"White Zinfandel", {{1.7501e-05, 1.9069e-05, 1.288e-05},
                              {0.012072, 0.016184, 0.019843}}},
      {"Merlot", {{2.1129e-05, 0, 0}, {0.11632, 0.25191, 0.29434}}},
      {"Budweiser Beer", {{2.4356e-05, 2.4079e-05, 1.0564e-05},
                             {0.011492, 0.024911, 0.057786}}},
      {"Coors Light Beer",
          {{5.0922e-05, 4.301e-05, 0}, {0.006164, 0.013984, 0.034983}}},
      {"Clorox",
          {{0.0024035, 0.0031373, 0.003991}, {0.0033542, 0.014892, 0.026297}}},
      {"Apple Juice",
          {{0.00013612, 0.00015836, 0.000227}, {0.012957, 0.023741, 0.052184}}},
      {"Cranberry Juice", {{0.00010402, 0.00011646, 7.8139e-05},
                              {0.039437, 0.094223, 0.12426}}},
      {"Grape Juice", {{5.382e-05, 0, 0}, {0.10404, 0.23958, 0.29325}}},
      {"Ruby Grapefruit Juice",
          {{0.011002, 0.010927, 0.011036}, {0.085867, 0.18314, 0.25262}}},
      {"White Grapefruit Juice",
          {{0.22826, 0.23998, 0.32748}, {0.0138, 0.018831, 0.056781}}},
      {"Shampoo",
          {{0.0007176, 0.0008303, 0.0009016}, {0.014107, 0.045693, 0.061717}}},
      {"Strawberry Shampoo",
          {{0.00015671, 0.00015947, 1.518e-05}, {0.01449, 0.05796, 0.075823}}},
      {"Head & Shoulders Shampoo",
          {{0.023805, 0.028804, 0.034306}, {0.084621, 0.15688, 0.20365}}},
      {"Lemon Tea Powder",
          {{0.040224, 0.045264, 0.051081}, {2.4288, 4.5757, 7.2127}}},
      {"Orange Powder", {{0.00015617, 0.00017482, 0.0001762},
                            {0.001449, 0.003441, 0.007863}}},
      {"Pink Lemonade Powder", {{0.00012103, 0.00013073, 0.00012528},
                                   {0.001165, 0.002366, 0.003195}}},
      {"Cappuccino Powder",
          {{1.8436, 2.5851, 2.1662}, {35.844, 49.547, 61.084}}},
      {"Salt Powder",
          {{0.027333, 0.032451, 0.031979}, {0.28415, 0.3257, 0.34148}}},
      {"Sugar Powder",
          {{0.00022272, 0.00025513, 0.000271}, {0.012638, 0.031051, 0.050124}}},
      {"Suisse Mocha Powder",
          {{2.7979, 3.5452, 4.3365}, {17.502, 27.004, 35.433}}},
      {"Pacific Ocean Surface Water", {{0.0001764, 0.00032095, 0.00019617},
                                          {0.031845, 0.031324, 0.030147}}},
  };
  return params.at(name);
}
#endif

}  // namespace yocto

#endif
