//
// # Yocto/Pbrt: Tiny library for Pbrt parsing and writing
//
// Yocto/Pbrt is a tiny library for loading and saving Pbrt.
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

#ifndef _YOCTO_PBRT_H_
#define _YOCTO_PBRT_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::pbrt {

// Math defitions
using math::frame3f;
using math::identity3x4f;
using math::identity4x4f;
using math::mat4f;
using math::pif;
using math::vec2f;
using math::vec2i;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;
using math::zero2f;
using math::zero2i;
using math::zero3f;
using math::zero3i;
using math::zero4f;

}  // namespace yocto::pbrt

// -----------------------------------------------------------------------------
// PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::pbrt {

// Pbrt camera
struct camera {
  // camera parameters
  frame3f frame      = identity3x4f;
  frame3f frend      = identity3x4f;
  vec2i   resolution = {0, 0};
  float   lens       = 0;
  float   aspect     = 0;
  float   focus      = 0;
  float   aperture   = 0;
};

// Pbrt material
struct material {
  // material parameters
  std::string name            = "";
  vec3f       emission        = {0, 0, 0};
  vec3f       color           = {0, 0, 0};
  float       specular        = 0;
  float       metallic        = 0;
  float       transmission    = 0;
  float       roughness       = 0;
  float       ior             = 1.5;
  float       opacity         = 1;
  std::string color_tex       = "";
  std::string opacity_tex     = "";
  std::string alpha_tex       = "";
  bool        thin            = true;
  vec3f       volmeanfreepath = {0, 0, 0};
  vec3f       volscatter      = {0, 0, 0};
  float       volscale        = 0.01;
};

// Pbrt shape
struct shape {
  // frames
  frame3f              frame     = identity3x4f;
  frame3f              frend     = identity3x4f;
  std::vector<frame3f> instances = {};
  std::vector<frame3f> instaends = {};
  // shape
  std::string        filename_ = "";
  std::vector<vec3f> positions = {};
  std::vector<vec3f> normals   = {};
  std::vector<vec2f> texcoords = {};
  std::vector<vec3i> triangles = {};
  // material
  pbrt::material* material = nullptr;
};

// Pbrt lights
struct light {
  // light parameters
  frame3f frame    = identity3x4f;
  frame3f frend    = identity3x4f;
  vec3f   emission = {0, 0, 0};
  vec3f   from     = {0, 0, 0};
  vec3f   to       = {0, 0, 0};
  bool    distant  = false;
  // arealight approximation
  vec3f              area_emission  = {0, 0, 0};
  frame3f            area_frame     = identity3x4f;
  frame3f            area_frend     = identity3x4f;
  std::vector<vec3i> area_triangles = {};
  std::vector<vec3f> area_positions = {};
  std::vector<vec3f> area_normals   = {};
};
struct environment {
  // environment approximation
  frame3f     frame        = identity3x4f;
  frame3f     frend        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  std::string emission_tex = "";
};

// Pbrt model
struct model {
  // pbrt data
  std::vector<std::string>        comments     = {};
  std::vector<pbrt::camera*>      cameras      = {};
  std::vector<pbrt::shape*>       shapes       = {};
  std::vector<pbrt::environment*> environments = {};
  std::vector<pbrt::light*>       lights       = {};
  std::vector<pbrt::material*>    materials    = {};

  // cleanup
  ~model();
};

// Load/save pbrt
inline bool load_pbrt(
    const std::string& filename, pbrt::model* pbrt, std::string& error);
inline bool save_pbrt(const std::string& filename, pbrt::model* pbrt,
    std::string& error, bool ply_meshes = false);

// Create pbrt
inline pbrt::camera*      add_camera(pbrt::model* pbrt);
inline pbrt::shape*       add_shape(pbrt::model* pbrt);
inline pbrt::material*    add_material(pbrt::model* pbrt);
inline pbrt::environment* add_environment(pbrt::model* pbrt);
inline pbrt::light*       add_light(pbrt::model* pbrt);

}  // namespace yocto::pbrt

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
#include <unordered_set>

#include "ext/filesystem.hpp"
#include "yocto_ply.h"
namespace sfs = ghc::filesystem;

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::pbrt {

// string literals
using namespace std::string_literals;

}  // namespace yocto::pbrt

// -----------------------------------------------------------------------------
// PBRT PARSING
// -----------------------------------------------------------------------------
namespace yocto::pbrt {

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
[[nodiscard]] inline bool parse_value(std::string_view& str, int& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
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
[[nodiscard]] inline bool parse_value(std::string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

// Formats values to std::string
inline void format_value(std::string& str, const std::string& value) {
  str += value;
}
inline void format_value(std::string& str, const char* value) { str += value; }
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
inline void format_value(std::string& str, const vec4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(std::string& str, const mat4f& value) {
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

// Pbrt value
struct value {
  // Pbrt value type
  enum struct type_t {
    // clang-format off
    real, integer, boolean, string, point, normal, vector, texture, color, 
    point2, vector2, spectrum
    // clang-format on
  };

  std::string        name     = "";
  type_t             type     = type_t::real;
  int                value1i  = 0;
  float              value1f  = 0;
  vec2f              value2f  = {0, 0};
  vec3f              value3f  = {0, 0, 0};
  bool               value1b  = false;
  std::string        value1s  = "";
  std::vector<float> vector1f = {};
  std::vector<vec2f> vector2f = {};
  std::vector<vec3f> vector3f = {};
  std::vector<int>   vector1i = {};
};

// Pbrt command
struct command {
  std::string        name   = "";
  std::string        type   = "";
  std::vector<value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
};

// get pbrt value
[[nodiscard]] inline bool get_value(const value& pbrt, std::string& val) {
  if (pbrt.type == value::type_t::string ||
      pbrt.type == value::type_t::texture) {
    val = pbrt.value1s;
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(const value& pbrt, bool& val) {
  if (pbrt.type == value::type_t::boolean) {
    val = pbrt.value1b;
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(const value& pbrt, int& val) {
  if (pbrt.type == value::type_t::integer) {
    val = pbrt.value1i;
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(const value& pbrt, float& val) {
  if (pbrt.type == value::type_t::real) {
    val = pbrt.value1f;
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(const value& pbrt, vec2f& val) {
  if (pbrt.type == value::type_t::point2 ||
      pbrt.type == value::type_t::vector2) {
    val = pbrt.value2f;
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(const value& pbrt, vec3f& val) {
  if (pbrt.type == value::type_t::point || pbrt.type == value::type_t::vector ||
      pbrt.type == value::type_t::normal || pbrt.type == value::type_t::color) {
    val = pbrt.value3f;
    return true;
  } else if (pbrt.type == value::type_t::real) {
    val = vec3f{pbrt.value1f};
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::vector<float>& val) {
  if (pbrt.type == value::type_t::real) {
    if (!pbrt.vector1f.empty()) {
      val = pbrt.vector1f;
    } else {
      val = {pbrt.value1f};
    }
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::vector<vec2f>& val) {
  if (pbrt.type == value::type_t::point2 ||
      pbrt.type == value::type_t::vector2) {
    if (!pbrt.vector2f.empty()) {
      val = pbrt.vector2f;
    } else {
      val = {pbrt.value2f};
    }
    return true;
  } else if (pbrt.type == value::type_t::real) {
    if (pbrt.vector1f.empty() || pbrt.vector1f.size() % 2)
      throw std::runtime_error("bad pbrt type");
    val.resize(pbrt.vector1f.size() / 2);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1f[i * 2 + 0], pbrt.vector1f[i * 2 + 1]};
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::vector<vec3f>& val) {
  if (pbrt.type == value::type_t::point || pbrt.type == value::type_t::vector ||
      pbrt.type == value::type_t::normal || pbrt.type == value::type_t::color) {
    if (!pbrt.vector3f.empty()) {
      val = pbrt.vector3f;
    } else {
      val = {pbrt.value3f};
    }
    return true;
  } else if (pbrt.type == value::type_t::real) {
    if (pbrt.vector1f.empty() || pbrt.vector1f.size() % 3)
      throw std::invalid_argument{"expected float3 array"};
    val.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}

[[nodiscard]] inline bool get_value(const value& pbrt, std::vector<int>& val) {
  if (pbrt.type == value::type_t::integer) {
    if (!pbrt.vector1i.empty()) {
      val = pbrt.vector1i;
    } else {
      val = {pbrt.vector1i};
    }
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::vector<vec3i>& val) {
  if (pbrt.type == value::type_t::integer) {
    if (pbrt.vector1i.empty() || pbrt.vector1i.size() % 3)
      throw std::invalid_argument{"expected int3 array"};
    val.resize(pbrt.vector1i.size() / 3);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1i[i * 3 + 0], pbrt.vector1i[i * 3 + 1],
          pbrt.vector1i[i * 3 + 2]};
    return true;
  } else {
    return false;
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::pair<float, std::string>& val) {
  if (pbrt.type == value::type_t::string ||
      pbrt.type == value::type_t::texture) {
    val.first = 0;
    return get_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_value(pbrt, val.first);
  }
}
[[nodiscard]] inline bool get_value(
    const value& pbrt, std::pair<vec3f, std::string>& val) {
  if (pbrt.type == value::type_t::string ||
      pbrt.type == value::type_t::texture) {
    val.first = zero3f;
    return get_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_value(pbrt, val.first);
  }
}
template <typename T>
[[nodiscard]] inline bool get_value(
    const std::vector<value>& pbrt, const std::string& name, T& val) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_value(p, val);
    }
  }
  return true;
}

// pbrt value construction
inline value make_value(const std::string& name, const std::string& val,
    value::type_t type = value::type_t::string) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = val;
  return pbrt;
}
inline value make_value(const std::string& name, bool val,
    value::type_t type = value::type_t::boolean) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = val;
  return pbrt;
}
inline value make_value(const std::string& name, int val,
    value::type_t type = value::type_t::integer) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = val;
  return pbrt;
}
inline value make_value(const std::string& name, float val,
    value::type_t type = value::type_t::real) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = val;
  return pbrt;
}
inline value make_value(const std::string& name, const vec2f& val,
    value::type_t type = value::type_t::point2) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = val;
  return pbrt;
}
inline value make_value(const std::string& name, const vec3f& val,
    value::type_t type = value::type_t::color) {
  auto pbrt    = value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = val;
  return pbrt;
}
inline value make_value(const std::string& name, const std::vector<vec2f>& val,
    value::type_t type = value::type_t::point2) {
  auto pbrt     = value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = val;
  return pbrt;
}
inline value make_value(const std::string& name, const std::vector<vec3f>& val,
    value::type_t type = value::type_t::point) {
  auto pbrt     = value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = val;
  return pbrt;
}
inline value make_value(const std::string& name, const std::vector<vec3i>& val,
    value::type_t type = value::type_t::integer) {
  auto pbrt     = value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {(int*)val.data(), (int*)val.data() + val.size() * 3};
  return pbrt;
}

inline void remove_comment(std::string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy       = str;
  auto in_string = false;
  while (!cpy.empty()) {
    if (cpy.front() == '"') in_string = !in_string;
    if (cpy.front() == comment_char && !in_string) break;
    cpy.remove_prefix(1);
  }
  str.remove_suffix(cpy.size());
}

// Read a pbrt command from file
[[nodiscard]] inline bool read_cmdline(FILE* fs, std::string& cmd) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs);
  while (fgets(buffer, sizeof(buffer), fs)) {
    // line
    auto line = std::string_view{buffer};
    remove_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs, pos, SEEK_SET);
        // line_num -= 1;
        return true;
      } else {
        found = true;
      }
    } else if (!found) {
      return false;
    }
    cmd += line;
    cmd += " ";
    pos = ftell(fs);
  }
  return found;
}

// parse a quoted std::string
[[nodiscard]] inline bool parse_command(
    std::string_view& str, std::string& value) {
  skip_whitespace(str);
  if (!isalpha((int)str.front())) return false;
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == std::string_view::npos) {
    value.assign(str);
    str.remove_prefix(str.size());
  } else {
    value.assign(str.substr(0, pos));
    str.remove_prefix(pos + 1);
  }
  return true;
}

// parse pbrt value with optional parens
template <typename T>
[[nodiscard]] inline bool parse_param(std::string_view& str, T& value) {
  skip_whitespace(str);
  auto parens = !str.empty() && str.front() == '[';
  if (parens) str.remove_prefix(1);
  if (!parse_value(str, value)) return false;
  if (!str.data()) return false;
  if (parens) {
    skip_whitespace(str);
    if (!str.empty() && str.front() == '[') return false;
    str.remove_prefix(1);
  }
  return true;
}

// parse a quoted std::string
[[nodiscard]] inline bool parse_nametype(
    std::string_view& str_, std::string& name, std::string& type) {
  auto value = ""s;
  if (!parse_value(str_, value)) return false;
  if (!str_.data()) return false;
  auto str  = std::string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == std::string_view::npos) return false;
  type = std::string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == std::string_view::npos) return false;
  str.remove_prefix(pos2);
  name = std::string(str);
  return true;
}

inline std::pair<vec3f, vec3f> get_etak(const std::string& name) {
  static const std::unordered_map<std::string, std::pair<vec3f, vec3f>>
      metal_ior_table = {
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

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
inline std::pair<vec3f, vec3f> get_subsurface(const std::string& name) {
  static const std::unordered_map<std::string, std::pair<vec3f, vec3f>> params =
      {
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
          {"Lowfat Milk",
              {{0.89187, 1.5136, 2.532}, {0.002875, 0.00575, 0.0115}}},
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
          {"Pepsi",
              {{6.1697e-05, 4.2564e-05, 0}, {0.091641, 0.14158, 0.20729}}},
          {"Sprite", {{6.0306e-06, 6.4139e-06, 6.5504e-06},
                         {0.001886, 0.0018308, 0.0020025}}},
          {"Gatorade", {{0.0024574, 0.003007, 0.0037325},
                           {0.024794, 0.019289, 0.008878}}},
          {"Chardonnay", {{1.7982e-05, 1.3758e-05, 1.2023e-05},
                             {0.010782, 0.011855, 0.023997}}},
          {"White Zinfandel", {{1.7501e-05, 1.9069e-05, 1.288e-05},
                                  {0.012072, 0.016184, 0.019843}}},
          {"Merlot", {{2.1129e-05, 0, 0}, {0.11632, 0.25191, 0.29434}}},
          {"Budweiser Beer", {{2.4356e-05, 2.4079e-05, 1.0564e-05},
                                 {0.011492, 0.024911, 0.057786}}},
          {"Coors Light Beer",
              {{5.0922e-05, 4.301e-05, 0}, {0.006164, 0.013984, 0.034983}}},
          {"Clorox", {{0.0024035, 0.0031373, 0.003991},
                         {0.0033542, 0.014892, 0.026297}}},
          {"Apple Juice", {{0.00013612, 0.00015836, 0.000227},
                              {0.012957, 0.023741, 0.052184}}},
          {"Cranberry Juice", {{0.00010402, 0.00011646, 7.8139e-05},
                                  {0.039437, 0.094223, 0.12426}}},
          {"Grape Juice", {{5.382e-05, 0, 0}, {0.10404, 0.23958, 0.29325}}},
          {"Ruby Grapefruit Juice",
              {{0.011002, 0.010927, 0.011036}, {0.085867, 0.18314, 0.25262}}},
          {"White Grapefruit Juice",
              {{0.22826, 0.23998, 0.32748}, {0.0138, 0.018831, 0.056781}}},
          {"Shampoo", {{0.0007176, 0.0008303, 0.0009016},
                          {0.014107, 0.045693, 0.061717}}},
          {"Strawberry Shampoo", {{0.00015671, 0.00015947, 1.518e-05},
                                     {0.01449, 0.05796, 0.075823}}},
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
          {"Sugar Powder", {{0.00022272, 0.00025513, 0.000271},
                               {0.012638, 0.031051, 0.050124}}},
          {"Suisse Mocha Powder",
              {{2.7979, 3.5452, 4.3365}, {17.502, 27.004, 35.433}}},
          {"Pacific Ocean Surface Water", {{0.0001764, 0.00032095, 0.00019617},
                                              {0.031845, 0.031324, 0.030147}}},
      };
  return params.at(name);
}

[[nodiscard]] inline bool parse_params(
    std::string_view& str, std::vector<value>& values) {
  auto parse_pvalues = [](std::string_view& str, auto& value,
                           auto& values) -> bool {
    values.clear();
    skip_whitespace(str);
    if (str.empty()) return false;
    if (str.front() == '[') {
      str.remove_prefix(1);
      skip_whitespace(str);
      if (str.empty()) return false;
      while (!str.empty()) {
        auto& val = values.empty() ? value : values.emplace_back();
        if (!parse_value(str, val)) return false;
        if (!str.data()) return false;
        skip_whitespace(str);
        if (str.empty()) break;
        if (str.front() == ']') break;
        if (values.empty()) values.push_back(value);
      }
      if (str.empty()) return false;
      if (str.front() != ']') return false;
      str.remove_prefix(1);
      return true;
    } else {
      return parse_value(str, value);
    }
  };

  values.clear();
  skip_whitespace(str);
  while (!str.empty()) {
    auto& value = values.emplace_back();
    auto  type  = ""s;
    if (!parse_nametype(str, value.name, type)) return false;
    skip_whitespace(str);
    if (str.empty()) return false;
    if (type == "float") {
      value.type = value::type_t::real;
      if (!parse_pvalues(str, value.value1f, value.vector1f)) return false;
    } else if (type == "integer") {
      value.type = value::type_t::integer;
      if (!parse_pvalues(str, value.value1i, value.vector1i)) return false;
    } else if (type == "string") {
      auto vector1s = std::vector<std::string>{};
      value.type    = value::type_t::string;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = std::vector<std::string>{};
      value.type    = value::type_t::boolean;
      if (!parse_pvalues(str, value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = std::vector<std::string>{};
      value.type    = value::type_t::texture;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "point" || type == "point3") {
      value.type = value::type_t::point;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "normal" || type == "normal3") {
      value.type = value::type_t::normal;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "std::vector" || type == "vector3") {
      value.type = value::type_t::vector;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "point2") {
      value.type = value::type_t::point2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "vector2") {
      value.type = value::type_t::vector2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "blackbody") {
      value.type     = value::type_t::color;
      auto blackbody = zero2f;
      auto vector2f  = std::vector<vec2f>{};
      parse_pvalues(str, blackbody, vector2f);
      if (!vector2f.empty()) return false;
      value.value3f = math::blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "color" || type == "rgb") {
      value.type = value::type_t::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "xyz") {
      value.type = value::type_t::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
      // xyz conversion
      return false;
    } else if (type == "spectrum") {
      auto is_string = false;
      auto str1      = str;
      skip_whitespace(str1);
      if (!str1.empty() && str1.front() == '"')
        is_string = true;
      else if (!str1.empty() && str1.front() == '[') {
        str1.remove_prefix(1);
        skip_whitespace(str1);
        if (!str1.empty() && str1.front() == '"') is_string = true;
      }
      if (is_string) {
        value.type     = value::type_t::color;
        auto filename  = ""s;
        auto filenames = std::vector<std::string>{};
        if (!parse_value(str, filename)) return false;
        if (!str.data()) return false;
        auto filenamep = sfs::path(filename).filename();
        if (sfs::path(filenamep).extension() == ".spd") {
          filenamep = sfs::path(filenamep).replace_extension("").string();
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (sfs::path(filenamep).extension() == ".eta") {
            auto eta =
                get_etak(sfs::path(filenamep).replace_extension("")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (sfs::path(filenamep).extension() == ".k") {
            auto k =
                get_etak(sfs::path(filenamep).replace_extension("")).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            return false;
          }
        } else {
          return false;
        }
      } else {
        value.type = value::type_t::spectrum;
        if (!parse_pvalues(str, value.value1f, value.vector1f)) return false;
      }
    } else {
      return false;
    }
    skip_whitespace(str);
  }
  return true;
}

// Other pbrt elements
struct film {
  // film approximation
  std::string filename   = "";
  vec2i       resolution = {0, 0};
};

// Pbrt texture
struct texture {
  // texture parameters
  std::string name     = "";
  vec3f       constant = {1, 1, 1};
  std::string filename = "";
};

// Pbrt area light
struct arealight {
  // arealight parameters
  std::string name     = "";
  vec3f       emission = {0, 0, 0};
};

// Pbrt medium. Not parsed at the moment.
struct medium {
  // medium parameters
  std::string name = "";
};

// convert pbrt films
inline bool convert_film(pbrt::film* film, const command& command,
    const std::string& filename, std::string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  if (command.type == "image") {
    film->resolution = {512, 512};
    if (!get_value(command.values, "xresolution", film->resolution.x))
      return parse_error();
    if (!get_value(command.values, "yresolution", film->resolution.y))
      return parse_error();
    film->filename = "out.png"s;
    if (!get_value(command.values, "filename", film->filename))
      return parse_error();
    return true;
  } else {
    return type_error();
  }
}

// convert pbrt elements
inline bool convert_camera(pbrt::camera* pcamera, const command& command,
    const vec2i& resolution, const std::string& filename, std::string& error,
    bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  pcamera->frame      = command.frame;
  pcamera->frend      = command.frend;
  pcamera->frame      = inverse((frame3f)pcamera->frame);
  pcamera->frame.z    = -pcamera->frame.z;
  pcamera->resolution = resolution;
  auto film_aspect =
      (resolution == zero2i) ? 1 : (float)resolution.x / (float)resolution.y;
  if (command.type == "perspective") {
    auto fov = 90.0f;
    if (!get_value(command.values, "fov", fov)) return parse_error();
    // auto lensradius = if(!get_value(values, "lensradius", 0.0f);
    pcamera->aspect = film_aspect;
    if (pcamera->aspect >= 1) {
      pcamera->lens = (0.036 / pcamera->aspect) /
                      (2 * math::tan(math::radians(fov) / 2));
    } else {
      pcamera->lens = (0.036 * pcamera->aspect) /
                      (2 * math::tan(math::radians(fov) / 2));
    }
    if (!get_value(command.values, "frameaspectratio", pcamera->aspect))
      return parse_error();
    pcamera->focus = 10.0f;
    if (!get_value(command.values, "focaldistance", pcamera->focus))
      return parse_error();
    return true;
  } else if (command.type == "realistic") {
    auto lensfile = ""s;
    if (!get_value(command.values, "lensfile", lensfile)) return parse_error();
    lensfile          = lensfile.substr(0, lensfile.size() - 4);
    lensfile          = lensfile.substr(lensfile.find('.') + 1);
    lensfile          = lensfile.substr(0, lensfile.size() - 2);
    auto lens         = math::max(std::atof(lensfile.c_str()), 35.0f) * 0.001f;
    pcamera->lens     = 2 * atan(0.036f / (2 * lens));
    pcamera->aperture = 0.0f;
    if (!get_value(command.values, "aperturediameter", pcamera->aperture))
      return parse_error();
    pcamera->focus = 10.0f;
    if (!get_value(command.values, "focusdistance", pcamera->focus))
      return parse_error();
    pcamera->aspect = film_aspect;
    return true;
  } else {
    return type_error();
  }
}

// convert pbrt textures
inline bool convert_texture(pbrt::texture* ptexture, const command& command,
    std::unordered_map<std::string, texture>& texture_map,
    const std::string& filename, std::string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  auto get_filename = [&texture_map](const std::string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return pos->second.filename;
  };

  ptexture->name = command.name;
  if (command.type == "imagemap") {
    ptexture->filename = "";
    if (!get_value(command.values, "filename", ptexture->filename))
      return parse_error();
    return true;
  } else if (command.type == "constant") {
    ptexture->constant = vec3f{1};
    if (!get_value(command.values, "value", ptexture->constant))
      return parse_error();
    return true;
  } else if (command.type == "bilerp") {
    ptexture->constant = {1, 0, 0};
    return true;
  } else if (command.type == "checkerboard") {
    // auto tex1     = if(!get_value(command.values, "tex1", std::pair{vec3f{1},
    // ""s}); auto tex2     = if(!get_value(command.values, "tex2",
    //  std::pair{vec3f{0}, ""s}); auto rgb1     = tex1.second == "" ?
    //  tex1.first :
    // vec3f{0.4f, 0.4f, 0.4f}; auto rgb2     = tex1.second == "" ? tex2.first :
    // vec3f{0.6f, 0.6f, 0.6f}; auto params   = proc_image_params{}; params.type
    // = proc_image_params::type_t::checker; params.color0 = {rgb1.x, rgb1.y,
    // rgb1.z, 1}; params.color1 = {rgb2.x, rgb2.y, rgb2.z, 1}; params.scale
    // = 2; make_proc_image(texture.hdr, params); float_to_byte(texture.ldr,
    // texture.hdr); texture.hdr = {};
    ptexture->constant = {0.5, 0.5, 0.5};
    return true;
  } else if (command.type == "dots") {
    ptexture->constant = {0.5, 0.5, 0.5};
    return true;
  } else if (command.type == "fbm") {
    ptexture->constant = {0.5, 0.5, 0.5};
    return true;
  } else if (command.type == "marble") {
    ptexture->constant = {0.5, 0.5, 0.5};
    return true;
  } else if (command.type == "mix") {
    auto tex1 = std::pair{vec3f{0}, ""s}, tex2 = std::pair{vec3f{1}, ""s};
    if (!get_value(command.values, "tex1", tex1)) return parse_error();
    if (!get_value(command.values, "tex2", tex2)) return parse_error();
    if (!get_filename(tex1.second).empty()) {
      ptexture->filename = get_filename(tex1.second);
    } else if (!get_filename(tex2.second).empty()) {
      ptexture->filename = get_filename(tex2.second);
    } else {
      ptexture->constant = {1, 0, 0};
    }
    return true;
  } else if (command.type == "scale") {
    auto tex1 = std::pair{vec3f{1}, ""s}, tex2 = std::pair{vec3f{1}, ""s};
    if (!get_value(command.values, "tex1", tex2)) return parse_error();
    if (!get_value(command.values, "tex2", tex1)) return parse_error();
    if (!get_filename(tex1.second).empty()) {
      ptexture->filename = get_filename(tex1.second);
    } else if (!get_filename(tex2.second).empty()) {
      ptexture->filename = get_filename(tex2.second);
    } else {
      ptexture->constant = {1, 0, 0};
    }
    return true;
  } else if (command.type == "uv") {
    ptexture->constant = {1, 0, 0};
    return true;
  } else if (command.type == "windy") {
    ptexture->constant = {1, 0, 0};
    return true;
  } else if (command.type == "wrinkled") {
    ptexture->constant = {1, 0, 0};
    return true;
  } else {
    return type_error();
  }
}

// convert pbrt materials
inline bool convert_material(pbrt::material* pmaterial, const command& command,
    const std::unordered_map<std::string, material>& named_materials,
    const std::unordered_map<std::string, texture>&  named_textures,
    const std::string& filename, std::string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };
  auto material_error = [filename, &error](const std::string& name) {
    error = filename + ": missing material " + name;
    return false;
  };
  auto bsdf_error = [filename, &error](const std::string& name) {
    error = filename + ": missing bsdf " + name;
    return false;
  };

  // helpers
  auto get_texture = [&](const std::vector<value>& values,
                         const std::string& name, vec3f& color,
                         std::string& filename, const vec3f& def) -> bool {
    auto textured = std::pair{def, ""s};
    if (!get_value(values, name, textured)) return parse_error();
    if (textured.second == "") {
      color    = textured.first;
      filename = "";
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color    = texture.constant;
        filename = "";
      } else {
        color    = {1, 1, 1};
        filename = texture.filename;
      }
    }
    return true;
  };
  auto get_scalar = [&](const std::vector<value>& values,
                        const std::string& name, float& scalar,
                        float def) -> bool {
    auto textured = std::pair{vec3f{def}, ""s};
    if (!get_value(values, name, textured)) return parse_error();
    if (textured.second == "") {
      scalar = mean(textured.first);
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        scalar = mean(texture.constant);
      } else {
        scalar = def;
      }
    }
    return true;
  };
  auto get_color = [&](const std::vector<value>& values,
                       const std::string& name, vec3f& color,
                       const vec3f& def) -> bool {
    auto textured = std::pair{def, ""s};
    if (!get_value(values, name, textured)) return parse_error();
    if (textured.second == "") {
      color = textured.first;
    } else {
      auto& texture = named_textures.at(textured.second);
      if (texture.filename.empty()) {
        color = texture.constant;
      } else {
        color = def;
      }
    }
    return true;
  };

  auto get_roughness = [&](const std::vector<value>& values, float& roughness,
                           float def = 0.1) -> bool {
    auto roughness_ = std::pair{vec3f{def}, ""s};
    if (!get_value(values, "roughness", roughness_)) return parse_error();
    auto uroughness = roughness_, vroughness = roughness_;
    auto remaproughness = true;
    if (!get_value(values, "uroughness", uroughness)) return parse_error();
    if (!get_value(values, "vroughness", vroughness)) return parse_error();
    if (!get_value(values, "remaproughness", remaproughness))
      return parse_error();

    roughness = 0;
    if (uroughness.first == zero3f || vroughness.first == zero3f) return true;
    roughness = mean(vec2f{mean(uroughness.first), mean(vroughness.first)});
    // from pbrt code
    if (remaproughness) {
      roughness = math::max(roughness, 1e-3f);
      auto x    = math::log(roughness);
      roughness = 1.62142f + 0.819955f * x + 0.1734f * x * x +
                  0.0171201f * x * x * x + 0.000640711f * x * x * x * x;
    }
    roughness = sqrt(roughness);
    return true;
  };

  auto eta_to_reflectivity = [](const vec3f&  eta,
                                 const vec3f& etak = zero3f) -> vec3f {
    return ((eta - 1) * (eta - 1) + etak * etak) /
           ((eta + 1) * (eta + 1) + etak * etak);
  };

  pmaterial->name = command.name;
  if (command.type == "uber") {
    auto diffuse = zero3f, specular = zero3f, transmission = zero3f;
    auto diffuse_map = ""s, specular_map = ""s, transmission_map = ""s;
    if (!get_texture(command.values, "Kd", diffuse, diffuse_map, vec3f{0.25}))
      return parse_error();
    if (!get_texture(command.values, "Ks", specular, specular_map, vec3f{0.25}))
      return parse_error();
    if (!get_texture(
            command.values, "Kt", transmission, transmission_map, vec3f{0}))
      return parse_error();
    if (max(transmission) > 0.1) {
      pmaterial->color        = transmission;
      pmaterial->color_tex    = transmission_map;
      pmaterial->specular     = 1;
      pmaterial->transmission = 1;
    } else {
      pmaterial->color     = diffuse;
      pmaterial->color_tex = diffuse_map;
      pmaterial->specular  = 1;
    }
    if (!get_scalar(command.values, "opacity", pmaterial->opacity, 1))
      return parse_error();
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    if (!get_roughness(command.values, pmaterial->roughness, 0.1f))
      return parse_error();
    return true;
  } else if (command.type == "plastic") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.25}))
      return parse_error();
    if (!get_scalar(command.values, "Ks", pmaterial->specular, 0.25))
      return parse_error();
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    pmaterial->roughness = 0.1f;
    if (!get_roughness(command.values, pmaterial->roughness, 0.1))
      return parse_error();
    return true;
  } else if (command.type == "translucent") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.25}))
      return parse_error();
    if (!get_scalar(command.values, "Ks", pmaterial->specular, 0.25))
      return parse_error();
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    if (!get_roughness(command.values, pmaterial->roughness, 0.1))
      return parse_error();
    return true;
  } else if (command.type == "matte") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5}))
      return parse_error();
    return true;
  } else if (command.type == "mirror") {
    if (!get_texture(command.values, "Kr", pmaterial->color,
            pmaterial->color_tex, vec3f{0.9}))
      return parse_error();
    pmaterial->metallic  = 1;
    pmaterial->roughness = 0;
    return true;
  } else if (command.type == "metal") {
    // get_texture(
    //     values, "Kr", material->specular, material->specular_tex,
    //     vec3f{1});
    auto eta = zero3f, etak = zero3f;
    if (!get_color(command.values, "eta", eta,
            vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f}))
      return parse_error();
    if (!get_color(command.values, "k", etak,
            vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f}))
      return parse_error();
    pmaterial->color     = eta_to_reflectivity(eta, etak);
    pmaterial->roughness = 0.01f;
    if (!get_roughness(command.values, pmaterial->roughness, 0.01))
      return parse_error();
    return true;
  } else if (command.type == "substrate") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5}))
      return parse_error();
    if (!get_scalar(command.values, "Ks", pmaterial->specular, 0.5))
      return parse_error();
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    pmaterial->roughness = 0.1f;
    if (!get_roughness(command.values, pmaterial->roughness, 0.1))
      return parse_error();
    return true;
  } else if (command.type == "glass") {
    // get_texture(
    //     values, "Kr", material->specular, material->specular_tex,
    //     vec3f{1});
    // get_texture(command.values, "Kt", material->transmission,
    //     material->transmission_tex, vec3f{1});
    pmaterial->color        = {1, 1, 1};
    pmaterial->specular     = 1;
    pmaterial->transmission = 1;
    pmaterial->thin         = false;
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    pmaterial->roughness = 0;
    if (!get_roughness(command.values, pmaterial->roughness, 0))
      return parse_error();
    return true;
  } else if (command.type == "hair") {
    if (!get_texture(command.values, "color", pmaterial->color,
            pmaterial->color_tex, vec3f{0}))
      return parse_error();
    pmaterial->roughness = 1;
    if (verbose) printf("hair material not properly supported\n");
    return true;
  } else if (command.type == "disney") {
    if (!get_texture(command.values, "color", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5}))
      return parse_error();
    pmaterial->roughness = 1;
    if (verbose) printf("disney material not properly supported\n");
    return true;
  } else if (command.type == "kdsubsurface") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5}))
      return parse_error();
    if (!get_scalar(command.values, "Kr", pmaterial->specular, 1))
      return parse_error();
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    pmaterial->roughness = 0;
    if (!get_roughness(command.values, pmaterial->roughness, 0))
      return parse_error();
    if (verbose) printf("kdsubsurface material not properly supported\n");
    return true;
  } else if (command.type == "subsurface") {
    if (!get_scalar(command.values, "Kr", pmaterial->specular, 1))
      return parse_error();
    if (!get_scalar(command.values, "Kt", pmaterial->transmission, 1))
      return parse_error();
    pmaterial->color = {1, 1, 1};
    if (!get_scalar(command.values, "eta", pmaterial->ior, 1.5))
      return parse_error();
    pmaterial->roughness = 0;
    if (!get_roughness(command.values, pmaterial->roughness, 0))
      return parse_error();
    auto scale = 1.0f;
    if (!get_value(command.values, "scale", scale)) return parse_error();
    pmaterial->volscale = 1 / scale;
    auto sigma_a = zero3f, sigma_s = zero3f;
    auto sigma_a_tex = ""s, sigma_s_tex = ""s;
    if (!get_texture(command.values, "sigma_a", sigma_a, sigma_a_tex,
            vec3f{0011, .0024, .014}))
      return parse_error();
    if (!get_texture(command.values, "sigma_prime_s", sigma_s, sigma_s_tex,
            vec3f{2.55, 3.12, 3.77}))
      return parse_error();
    pmaterial->volmeanfreepath = 1 / (sigma_a + sigma_s);
    pmaterial->volscatter      = sigma_s / (sigma_a + sigma_s);
    if (verbose) printf("subsurface material not properly supported\n");
    return true;
  } else if (command.type == "mix") {
    auto namedmaterial1 = ""s, namedmaterial2 = ""s;
    if (!get_value(command.values, "namedmaterial1", namedmaterial1))
      return parse_error();
    if (!get_value(command.values, "namedmaterial2", namedmaterial2))
      return parse_error();
    auto matname = (!namedmaterial1.empty()) ? namedmaterial1 : namedmaterial2;
    auto matit   = named_materials.find(matname);
    if (matit == named_materials.end()) return material_error(matname);
    auto saved_name = pmaterial->name;
    *pmaterial      = matit->second;
    pmaterial->name = saved_name;
    if (verbose) printf("mix material not properly supported\n");
    return true;
  } else if (command.type == "fourier") {
    auto bsdffile = ""s;
    if (!get_value(command.values, "bsdffile", bsdffile)) return parse_error();
    if (bsdffile.rfind("/") != std::string::npos)
      bsdffile = bsdffile.substr(bsdffile.rfind("/") + 1);
    if (bsdffile == "paint.bsdf") {
      pmaterial->color     = {0.6f, 0.6f, 0.6f};
      pmaterial->specular  = 1;
      pmaterial->ior       = 1.5;
      pmaterial->roughness = 0.2;
    } else if (bsdffile == "ceramic.bsdf") {
      pmaterial->color     = {0.6f, 0.6f, 0.6f};
      pmaterial->specular  = 1;
      pmaterial->ior       = 1.5;
      pmaterial->roughness = 0.25;
    } else if (bsdffile == "leather.bsdf") {
      pmaterial->color     = {0.6f, 0.57f, 0.48f};
      pmaterial->specular  = 1;
      pmaterial->ior       = 1.5;
      pmaterial->roughness = 0.3;
    } else if (bsdffile == "coated_copper.bsdf") {
      auto eta             = vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f};
      auto etak            = vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f};
      pmaterial->color     = eta_to_reflectivity(eta, etak);
      pmaterial->metallic  = 1;
      pmaterial->roughness = 0.01;
    } else if (bsdffile == "roughglass_alpha_0.2.bsdf") {
      pmaterial->color        = {1, 1, 1};
      pmaterial->specular     = 1;
      pmaterial->ior          = 1.5;
      pmaterial->transmission = 1;
      pmaterial->roughness    = 0.2;
    } else if (bsdffile == "roughgold_alpha_0.2.bsdf") {
      auto eta             = vec3f{0.1431189557f, 0.3749570432f, 1.4424785571f};
      auto etak            = vec3f{3.9831604247f, 2.3857207478f, 1.6032152899f};
      pmaterial->color     = eta_to_reflectivity(eta, etak);
      pmaterial->metallic  = 1;
      pmaterial->roughness = 0.2;
    } else {
      return bsdf_error(bsdffile);
    }
    return true;
  } else {
    return type_error();
  }
}

// Make a triangle shape from a quad grid
template <typename PositionFunc, typename NormalFunc>
inline void make_shape(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec2i& steps,
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
inline void make_sphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec2i& steps, float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return radius * vec3f{math::cos(pt.x) * math::sin(pt.y),
                            math::sin(pt.x) * math::sin(pt.y), math::cos(pt.y)};
      },
      [](const vec2f& uv) {
        auto pt = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        return vec3f{math::cos(pt.x) * math::cos(pt.y),
            math::sin(pt.x) * math::cos(pt.y), math::sin(pt.y)};
      });
}
inline void make_disk(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec2i& steps, float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto a = 2 * pif * uv.x;
        return radius * (1 - uv.y) * vec3f{math::cos(a), math::sin(a), 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}
inline void make_quad(std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, const vec2i& steps, float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        return vec3f{(uv.x - 0.5f) * radius, (uv.y - 0.5f) * radius, 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}

// Convert pbrt shapes
inline bool convert_shape(pbrt::shape* shape, const command& command,
    std::string&                                    alphamap,
    const std::unordered_map<std::string, texture>& named_textures,
    const std::string& ply_dirname, const std::string& filename,
    std::string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };
  auto dependent_error = [filename, &error]() {
    error = filename + ": error in " + error;
    return false;
  };

  // helpers
  auto get_alpha = [&](const std::vector<value>& values,
                       const std::string& name, std::string& filename) -> bool {
    auto def      = 1.0f;
    auto textured = std::pair{def, ""s};
    if (!get_value(values, name, textured)) return parse_error();
    if (textured.second == "") {
      filename = "";
    } else {
      filename = named_textures.at(textured.second).filename;
    }
    return true;
  };

  shape->frame = command.frame;
  shape->frend = command.frend;
  if (command.type == "trianglemesh") {
    shape->positions = {};
    shape->normals   = {};
    shape->texcoords = {};
    shape->triangles = {};
    if (!get_value(command.values, "P", shape->positions)) return parse_error();
    if (!get_value(command.values, "N", shape->normals)) return parse_error();
    if (!get_value(command.values, "uv", shape->texcoords))
      return parse_error();
    for (auto& uv : shape->texcoords) uv.y = (1 - uv.y);
    if (!get_value(command.values, "indices", shape->triangles))
      return parse_error();
    return true;
  } else if (command.type == "loopsubdiv") {
    shape->positions = {};
    shape->triangles = {};
    if (!get_value(command.values, "P", shape->positions)) return parse_error();
    if (!get_value(command.values, "indices", shape->triangles))
      return parse_error();
    shape->normals.resize(shape->positions.size());
    // compute_normals(shape->normals, shape->triangles, shape->positions);
    return true;
  } else if (command.type == "plymesh") {
    shape->filename_ = ""s;
    if (!get_value(command.values, "filename", shape->filename_))
      return parse_error();
    if (!get_alpha(command.values, "alpha", alphamap)) return parse_error();
    auto ply = std::make_unique<ply::model>();
    if (!load_ply(ply_dirname + shape->filename_, ply.get(), error))
      return dependent_error();
    get_positions(ply.get(), shape->positions);
    get_normals(ply.get(), shape->normals);
    get_texcoords(ply.get(), shape->texcoords);
    get_triangles(ply.get(), shape->triangles);
    return true;
  } else if (command.type == "sphere") {
    auto radius = 1.0f;
    if (!get_value(command.values, "radius", radius)) return parse_error();
    make_sphere(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 16}, radius);
    return true;
  } else if (command.type == "disk") {
    auto radius = 1.0f;
    if (!get_value(command.values, "radius", radius)) return parse_error();
    make_disk(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 1}, radius);
    return true;
  } else {
    return type_error();
  }
}

// Convert pbrt arealights
inline bool convert_arealight(pbrt::arealight* parealight,
    const command& command, const std::string& filename, std::string& error,
    bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  parealight->name = command.name;
  if (command.type == "diffuse") {
    auto l = vec3f{1}, scale = vec3f{1};
    if (!get_value(command.values, "L", l)) return parse_error();
    if (!get_value(command.values, "scale", scale)) return parse_error();
    parealight->emission = l * scale;
    return true;
  } else {
    return type_error();
  }
}

// Convert pbrt lights
inline bool convert_light(pbrt::light* plight, const command& command,
    const std::string& filename, std::string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  plight->frame = command.frame;
  plight->frend = command.frend;
  if (command.type == "distant") {
    auto l = vec3f{1}, scale = vec3f{1};
    if (!get_value(command.values, "L", l)) return parse_error();
    if (!get_value(command.values, "scale", scale)) return parse_error();
    plight->emission = l * scale;
    plight->from     = vec3f{0, 0, 0};
    plight->to       = vec3f{0, 0, 1};
    if (!get_value(command.values, "from", plight->from)) return parse_error();
    if (!get_value(command.values, "to", plight->to)) return parse_error();
    plight->distant       = true;
    auto distant_dist     = 100;
    auto size             = distant_dist * math::sin(5 * pif / 180);
    plight->area_emission = plight->emission * (distant_dist * distant_dist) /
                            (size * size);
    plight->area_frame =
        plight->frame *
        lookat_frame(normalize(plight->from - plight->to) * distant_dist,
            {0, 0, 0}, {0, 1, 0}, true);
    plight->area_frend =
        plight->frend *
        lookat_frame(normalize(plight->from - plight->to) * distant_dist,
            {0, 0, 0}, {0, 1, 0}, true);
    auto texcoords = std::vector<vec2f>{};
    make_quad(plight->area_triangles, plight->area_positions,
        plight->area_normals, texcoords, {4, 2}, size);
    return true;
  } else if (command.type == "point" || command.type == "goniometric" ||
             command.type == "spot") {
    auto i = vec3f{1}, scale = vec3f{1};
    if (!get_value(command.values, "I", i)) return parse_error();
    if (!get_value(command.values, "scale", scale)) return parse_error();
    plight->emission = i * scale;
    plight->from     = zero3f;
    if (!get_value(command.values, "from", plight->from)) return parse_error();
    plight->area_emission = plight->emission;
    plight->area_frame    = plight->frame * translation_frame(plight->from);
    plight->area_frend    = plight->frend * translation_frame(plight->from);
    auto texcoords        = std::vector<vec2f>{};
    make_sphere(plight->area_triangles, plight->area_positions,
        plight->area_normals, texcoords, {4, 2}, 0.0025f);
    return true;
  } else {
    return type_error();
  }
}

inline bool convert_environment(pbrt::environment* penvironment,
    const command& command, const std::string& filename, std::string& error,
    bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  penvironment->frame = command.frame;
  penvironment->frend = command.frend;
  penvironment->frame = penvironment->frame *
                        frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  penvironment->frend = penvironment->frend *
                        frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  if (command.type == "infinite") {
    auto l = vec3f{1}, scale = vec3f{1};
    if (!get_value(command.values, "L", l)) return parse_error();
    if (!get_value(command.values, "scale", scale)) return parse_error();
    penvironment->emission     = scale * l;
    penvironment->emission_tex = ""s;
    if (!get_value(command.values, "mapname", penvironment->emission_tex))
      return parse_error();
    return true;
  } else {
    return type_error();
  }
}

// pbrt stack ctm
struct stack_element {
  frame3f         transform_start        = identity3x4f;
  frame3f         transform_end          = identity3x4f;
  pbrt::material  material               = {};
  pbrt::arealight arealight              = {};
  pbrt::medium    interior               = {};
  pbrt::medium    exterior               = {};
  bool            reverse                = false;
  bool            active_transform_start = true;
  bool            active_transform_end   = true;
};

// pbrt parsing context
struct context {
  std::vector<stack_element>                                 stack      = {};
  std::unordered_map<std::string, stack_element>             coordsys   = {};
  std::unordered_map<std::string, std::vector<pbrt::shape*>> objects    = {};
  std::string                                                cur_object = "";
  vec2i film_resolution = {512, 512};
};

// load pbrt
[[nodiscard]] inline bool load_pbrt(const std::string& filename,
    pbrt::model* pbrt, std::string& error, context& ctx,
    std::unordered_map<std::string, pbrt::material*>& material_map,
    std::unordered_map<std::string, material>&        named_materials,
    std::unordered_map<std::string, texture>&         named_textures,
    std::unordered_map<std::string, medium>&          named_mediums,
    const std::string&                                ply_dirname) {
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
  auto command_error = [filename, &error](const std::string& cmd) {
    error = filename + ": unknown command " + cmd;
    return false;
  };

  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // helpers
  auto set_transform = [](stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= xform;
    if (ctx.active_transform_end) ctx.transform_end *= xform;
  };

  // init stack
  if (ctx.stack.empty()) ctx.stack.emplace_back();

  // parse command by command
  auto line = ""s;
  while (read_cmdline(fs, line)) {
    auto str = std::string_view{line};
    // get command
    auto cmd = ""s;
    if (!parse_command(str, cmd)) return parse_error();
    if (cmd == "WorldBegin") {
      ctx.stack.push_back({});
    } else if (cmd == "WorldEnd") {
      if (ctx.stack.empty())
        throw std::runtime_error{filename + ": parse error [bad stack]"};
      ctx.stack.pop_back();
      if (ctx.stack.size() != 1)
        throw std::runtime_error{filename + ": parse error [bad stack]"};
    } else if (cmd == "AttributeBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "AttributeEnd") {
      if (ctx.stack.empty())
        throw std::runtime_error{filename + ": parse error [bad stack]"};
      ctx.stack.pop_back();
    } else if (cmd == "TransformBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "TransformEnd") {
      if (ctx.stack.empty())
        throw std::runtime_error{filename + ": parse error [bad stack]"};
      ctx.stack.pop_back();
    } else if (cmd == "ObjectBegin") {
      ctx.stack.push_back(ctx.stack.back());
      if (!parse_param(str, ctx.cur_object)) return parse_error();
      ctx.objects[ctx.cur_object] = {};
    } else if (cmd == "ObjectEnd") {
      ctx.stack.pop_back();
      ctx.cur_object = "";
    } else if (cmd == "ObjectInstance") {
      auto object = ""s;
      if (!parse_param(str, object)) return parse_error();
      if (ctx.objects.find(object) == ctx.objects.end())
        throw std::runtime_error{filename + ": parse error [unknown object]"};
      for (auto shape : ctx.objects.at(object)) {
        shape->instances.push_back(ctx.stack.back().transform_start);
        shape->instaends.push_back(ctx.stack.back().transform_end);
      }
    } else if (cmd == "ActiveTransform") {
      auto name = ""s;
      if (!parse_command(str, name)) return parse_error();
      if (name == "StartTime") {
        ctx.stack.back().active_transform_start = true;
        ctx.stack.back().active_transform_end   = false;
      } else if (name == "EndTime") {
        ctx.stack.back().active_transform_start = false;
        ctx.stack.back().active_transform_end   = true;
      } else if (name == "All") {
        ctx.stack.back().active_transform_start = true;
        ctx.stack.back().active_transform_end   = true;
      } else {
        throw std::runtime_error{filename + ": parse error [bad coordsys]"};
      }
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      if (!parse_param(str, xf)) return parse_error();
      set_transform(ctx.stack.back(), frame3f{xf});
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      if (!parse_param(str, xf)) return parse_error();
      concat_transform(ctx.stack.back(), frame3f{xf});
    } else if (cmd == "Scale") {
      auto v = zero3f;
      if (!parse_param(str, v)) return parse_error();
      concat_transform(ctx.stack.back(), scaling_frame(v));
    } else if (cmd == "Translate") {
      auto v = zero3f;
      if (!parse_param(str, v)) return parse_error();
      concat_transform(ctx.stack.back(), translation_frame(v));
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      if (!parse_param(str, v)) return parse_error();
      concat_transform(ctx.stack.back(),
          rotation_frame(vec3f{v.y, v.z, v.w}, math::radians(v.x)));
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      if (!parse_param(str, from)) return parse_error();
      if (!parse_param(str, to)) return parse_error();
      if (!parse_param(str, up)) return parse_error();
      auto frame = lookat_frame(from, to, up, true);
      concat_transform(ctx.stack.back(), inverse(frame));
    } else if (cmd == "ReverseOrientation") {
      ctx.stack.back().reverse = !ctx.stack.back().reverse;
    } else if (cmd == "CoordinateSystem") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      ctx.coordsys[name].transform_start = ctx.stack.back().transform_start;
      ctx.coordsys[name].transform_end   = ctx.stack.back().transform_end;
    } else if (cmd == "CoordSysTransform") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      if (ctx.coordsys.find(name) != ctx.coordsys.end()) {
        ctx.stack.back().transform_start =
            ctx.coordsys.at(name).transform_start;
        ctx.stack.back().transform_end = ctx.coordsys.at(name).transform_end;
      }
    } else if (cmd == "Integrator") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Sampler") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "PixelFilter") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Film") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      auto film = pbrt::film{};
      if (!convert_film(&film, command, filename, error)) return false;
      ctx.film_resolution = film.resolution;
    } else if (cmd == "Accelerator") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Camera") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto camera   = add_camera(pbrt);
      if (!convert_camera(
              camera, command, ctx.film_resolution, filename, error))
        return false;
    } else if (cmd == "Texture") {
      auto command  = pbrt::command{};
      auto comptype = ""s;
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_param(str, comptype)) return parse_error();
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (!convert_texture(&named_textures[command.name], command,
              named_textures, filename, error))
        return false;
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        command     = pbrt::command{};
      command.name = "__unnamed__material__" + std::to_string(material_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (command.type == "") {
        ctx.stack.back().material = {};
      } else {
        ctx.stack.back().material = {};
        if (!convert_material(&ctx.stack.back().material, command,
                named_materials, named_textures, filename, error))
          return false;
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (value.name == "type") command.type = value.value1s;
      if (!convert_material(&named_materials[command.name], command,
              named_materials, named_textures, filename, error))
        return false;
    } else if (cmd == "NamedMaterial") {
      auto name = ""s;
      if (!parse_param(str, name)) return parse_error();
      ctx.stack.back().material = named_materials.at(name);
    } else if (cmd == "Shape") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto shape    = add_shape(pbrt);
      auto alphamap = ""s;
      if (!convert_shape(shape, command, alphamap, named_textures, ply_dirname,
              filename, error))
        return false;
      auto matkey = "?!!!?" + ctx.stack.back().material.name + "?!!!?" +
                    ctx.stack.back().arealight.name + "?!!!?" + alphamap;
      if (material_map.find(matkey) == material_map.end()) {
        auto material  = add_material(pbrt);
        (*material)    = ctx.stack.back().material;
        material->name = "material" + std::to_string(pbrt->materials.size());
        material->emission   = ctx.stack.back().arealight.emission;
        material->alpha_tex  = alphamap;
        material_map[matkey] = material;
      }
      shape->material = material_map.at(matkey);
      if (ctx.cur_object != "") {
        ctx.objects[ctx.cur_object].push_back(shape);
      }
    } else if (cmd == "AreaLightSource") {
      static auto arealight_id = 0;
      auto        command      = pbrt::command{};
      command.name = "__unnamed__arealight__" + std::to_string(arealight_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (!convert_arealight(
              &ctx.stack.back().arealight, command, filename, error))
        return false;
    } else if (cmd == "LightSource") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (command.type == "infinite") {
        auto environment = add_environment(pbrt);
        if (!convert_environment(environment, command, filename, error))
          return false;
      } else {
        auto light = add_light(pbrt);
        if (!convert_light(light, command, filename, error)) return false;
      }
    } else if (cmd == "MakeNamedMedium") {
      auto command = pbrt::command{};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (command.name == "type") command.type = value.value1s;
      auto medium                 = pbrt::medium{};
      named_mediums[command.name] = medium;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      if (!parse_param(str, interior)) return parse_error();
      if (!parse_param(str, exterior)) return parse_error();
      ctx.stack.back().interior = named_mediums.at(interior);
      ctx.stack.back().exterior = named_mediums.at(exterior);
    } else if (cmd == "Include") {
      auto includename = ""s;
      if (!parse_param(str, includename)) return parse_error();
      if (!load_pbrt(sfs::path(filename).parent_path() / includename, pbrt,
              error, ctx, material_map, named_materials, named_textures,
              named_mediums, ply_dirname))
        return dependent_error();
    } else {
      return command_error(cmd);
    }
  }
  return true;
}

inline model::~model() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto environment : environments) delete environment;
  for (auto light : lights) delete light;
  for (auto material : materials) delete material;
}

// Make pbrt
inline pbrt::camera* add_camera(pbrt::model* pbrt) {
  return pbrt->cameras.emplace_back(new camera{});
}
inline pbrt::shape* add_shape(pbrt::model* pbrt) {
  return pbrt->shapes.emplace_back(new shape{});
}
inline pbrt::material* add_material(pbrt::model* pbrt) {
  return pbrt->materials.emplace_back(new material{});
}
inline pbrt::environment* add_environment(pbrt::model* pbrt) {
  return pbrt->environments.emplace_back(new environment{});
}
inline pbrt::light* add_light(pbrt::model* pbrt) {
  return pbrt->lights.emplace_back(new light{});
}

// load pbrt
[[nodiscard]] inline bool load_pbrt(
    const std::string& filename, pbrt::model* pbrt, std::string& error) {
  auto ctx             = context{};
  auto material_map    = std::unordered_map<std::string, pbrt::material*>{};
  auto named_materials = std::unordered_map<std::string, material>{{"", {}}};
  auto named_mediums   = std::unordered_map<std::string, medium>{{"", {}}};
  auto named_textures  = std::unordered_map<std::string, texture>{{"", {}}};
  auto dirname         = sfs::path(filename).parent_path().string();
  if (dirname != "") dirname += "/";
  if (!load_pbrt(filename, pbrt, error, ctx, material_map, named_materials,
          named_textures, named_mediums, dirname))
    return false;

  // remove unused materials
  auto used_materials = std::unordered_set<pbrt::material*>{};
  for (auto shape : pbrt->shapes) used_materials.insert(shape->material);
  pbrt->materials.erase(
      std::remove_if(pbrt->materials.begin(), pbrt->materials.end(),
          [&used_materials](auto material) {
            auto found = used_materials.find(material) != used_materials.end();
            if (!found) delete material;
            return !found;
          }),
      pbrt->materials.end());

  return true;
}

inline void format_value(std::string& str, const value& value) {
  static auto type_labels = std::unordered_map<value::type_t, std::string>{
      {value::type_t::real, "float"},
      {value::type_t::integer, "integer"},
      {value::type_t::boolean, "bool"},
      {value::type_t::string, "std::string"},
      {value::type_t::point, "point"},
      {value::type_t::normal, "normal"},
      {value::type_t::vector, "std::vector"},
      {value::type_t::texture, "texture"},
      {value::type_t::color, "rgb"},
      {value::type_t::point2, "point2"},
      {value::type_t::vector2, "vector2"},
      {value::type_t::spectrum, "spectrum"},
  };

  auto format_vector = [](std::string& str, auto& values) {
    str += "[ ";
    for (auto& value : values) {
      str += " ";
      format_value(str, value);
    }
    str += " ]";
  };

  format_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case value::type_t::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_value(str, value.value1f);
      }
      break;
    case value::type_t::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_value(str, value.value1i);
      }
      break;
    case value::type_t::boolean:
      format_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case value::type_t::string:
    case value::type_t::texture:
      format_values(str, "\"{}\"", value.value1s);
      break;
    case value::type_t::point:
    case value::type_t::vector:
    case value::type_t::normal:
    case value::type_t::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_values(str, "[ {} ]", value.value3f);
      }
      break;
    case value::type_t::spectrum: format_vector(str, value.vector1f); break;
    case value::type_t::point2:
    case value::type_t::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

inline void format_value(std::string& str, const std::vector<value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

[[nodiscard]] inline bool save_pbrt(const std::string& filename,
    pbrt::model* pbrt, std::string& error, bool ply_meshes) {
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
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // save comments
  if (!format_values(fs, "#\n")) return write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) return write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return write_error();
  if (!format_values(fs, "#\n\n")) return write_error();
  for (auto& comment : pbrt->comments) {
    if (!format_values(fs, "# {}\n", comment)) return write_error();
  }
  if (!format_values(fs, "\n")) return write_error();

  for (auto camera : pbrt->cameras) {
    auto command = pbrt::command{};
    command.type = "image";
    command.values.push_back(make_value("xresolution", camera->resolution.x));
    command.values.push_back(make_value("yresolution", camera->resolution.y));
    command.values.push_back(make_value("filename", "image.exr"s));
    if (!format_values(fs, "Film \"{}\" {}\n", command.type, command.values))
      return write_error();
  }

  for (auto camera : pbrt->cameras) {
    auto command  = pbrt::command{};
    command.type  = "perspective";
    command.frame = camera->frame;
    command.values.push_back(make_value(
        "fov", 2 * math::tan(0.036f / (2 * camera->lens)) * 180 / pif));
    if (!format_values(fs, "LookAt {} {} {}\n", command.frame.o,
            command.frame.o - command.frame.z, command.frame.y))
      return write_error();
    if (!format_values(fs, "Camera \"{}\" {}\n", command.type, command.values))
      return write_error();
  }

  if (!format_values(fs, "\nWorldBegin\n\n")) return write_error();

  for (auto light : pbrt->lights) {
    auto command  = pbrt::command{};
    command.frame = light->frame;
    if (light->distant) {
      command.type = "distance";
      command.values.push_back(make_value("L", light->emission));
    } else {
      command.type = "point";
      command.values.push_back(make_value("I", light->emission));
    }
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)command.frame))
      return write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
  }

  for (auto environment : pbrt->environments) {
    auto command  = pbrt::command{};
    command.frame = environment->frame;
    command.type  = "infinite";
    command.values.push_back(make_value("L", environment->emission));
    command.values.push_back(make_value("mapname", environment->emission_tex));
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)command.frame))
      return write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
  }

  auto reflectivity_to_eta = [](const vec3f& reflectivity) {
    return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
  };

  for (auto material : pbrt->materials) {
    auto command = pbrt::command{};
    if (material->specular != 0 && material->transmission != 0 &&
        !material->thin) {
      command.type = "glass";
      command.values.push_back(make_value("Kr", vec3f{1, 1, 1}));
      command.values.push_back(make_value("Kt", vec3f{1, 1, 1}));
      command.values.push_back(
          make_value("roughness", math::pow(material->roughness, 2)));
      command.values.push_back(make_value("eta", material->ior));
      command.values.push_back(make_value("remaproughness", false));
    } else if (material->metallic > 0.1f) {
      command.type = "metal";
      command.values.push_back(make_value("Kr", vec3f{1, 1, 1}));
      command.values.push_back(
          make_value("roughness", math::pow(material->roughness, 2)));
      command.values.push_back(
          make_value("eta", reflectivity_to_eta(material->color)));
      command.values.push_back(make_value("remaproughness", false));
    } else {
      command.type = "uber";
      if (material->color_tex.empty()) {
        command.values.push_back(make_value("Kd", material->color));
      } else if (material->color != zero3f) {
        command.values.push_back(
            make_value("Kd", material->color_tex, value::type_t::texture));
      }
      if (material->specular != 0) {
        command.values.push_back(make_value("Ks", vec3f{material->specular}));
        command.values.push_back(
            make_value("roughness", math::pow(material->roughness, 2)));
        command.values.push_back(make_value("eta", material->ior));
        command.values.push_back(make_value("remaproughness", false));
      }
      if (material->transmission != 0) {
        command.values.push_back(
            make_value("Kt", vec3f{material->transmission}));
      }
      if (!material->opacity_tex.empty()) {
        command.values.push_back(make_value(
            "opacity", material->opacity_tex, value::type_t::texture));
      } else if (material->opacity != 1) {
        command.values.push_back(make_value("opacity", material->opacity));
      }
    }
    if (!format_values(fs,
            "MakeNamedMaterial \"{}\" \"std::string type\" \"{}\" {}\n",
            material->name, command.type, command.values))
      return write_error();
  }

  auto object_id = 0;
  for (auto shape : pbrt->shapes) {
    auto command  = pbrt::command{};
    command.frame = shape->frame;
    if (ply_meshes) {
      command.type = "plymesh";
      command.values.push_back(make_value("filename", shape->filename_));
    } else {
      command.type = "trianglemesh";
      command.values.push_back(make_value("indices", shape->triangles));
      command.values.push_back(
          make_value("P", shape->positions, value::type_t::point));
      if (!shape->normals.empty())
        command.values.push_back(
            make_value("N", shape->triangles, value::type_t::normal));
      if (!shape->texcoords.empty())
        command.values.push_back(make_value("uv", shape->texcoords));
    }
    if (ply_meshes) {
      auto ply_guard = std::make_unique<ply::model>();
      auto ply       = ply_guard.get();
      add_positions(ply, shape->positions);
      add_normals(ply, shape->normals);
      add_texcoords(ply, shape->texcoords);
      add_triangles(ply, shape->triangles);
      if (!save_ply(
              sfs::path(filename).parent_path() / shape->filename_, ply, error))
        return dependent_error();
    }
    auto object = "object" + std::to_string(object_id++);
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectBegin \"{}\"\n", object))
        return write_error();
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)shape->frame))
      return write_error();
    if (shape->material->emission != zero3f) {
      auto acommand = pbrt::command{};
      acommand.type = "diffuse";
      acommand.values.push_back(make_value("L", shape->material->emission));
      if (!format_values(fs, "AreaLightSource \"{}\" {}\n", acommand.type,
              acommand.values))
        return write_error();
    }
    if (!format_values(fs, "NamedMaterial \"{}\"\n", shape->material->name))
      return write_error();
    if (!format_values(fs, "Shape \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectEnd\n")) return write_error();
    for (auto& iframe : shape->instances) {
      if (!format_values(fs, "AttributeBegin\n")) return write_error();
      if (!format_values(fs, "Transform {}\n", (mat4f)iframe))
        return write_error();
      if (!format_values(fs, "ObjectInstance \"{}\"\n", object))
        return write_error();
      if (!format_values(fs, "AttributeEnd\n")) return write_error();
    }
  }

  if (!format_values(fs, "\nWorldEnd\n\n")) return write_error();

  // done
  return true;
}

}  // namespace yocto::pbrt

#endif
