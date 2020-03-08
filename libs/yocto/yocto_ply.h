//
// # Yocto/Ply: Tiny library for Ply parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Ply.
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

#ifndef _YOCTO_PLY_H_
#define _YOCTO_PLY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// ALIASES
// -----------------------------------------------------------------------------
namespace yocto::ply {

// Math defitions
using math::byte;
using math::frame3f;
using math::vec2f;
using math::vec2i;
using math::vec3f;
using math::vec3i;
using math::vec4f;
using math::vec4i;

}  // namespace yocto::ply

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::ply {

// Ply property
struct property {
  // property types
  enum struct type_t { i8, i16, i32, i64, u8, u16, u32, u64, f32, f64 };

  // description
  std::string name    = "";
  bool        is_list = false;
  type_t      type    = type_t::f32;

  // data
  std::vector<int8_t>   data_i8  = {};
  std::vector<int16_t>  data_i16 = {};
  std::vector<int32_t>  data_i32 = {};
  std::vector<int64_t>  data_i64 = {};
  std::vector<uint8_t>  data_u8  = {};
  std::vector<uint16_t> data_u16 = {};
  std::vector<uint32_t> data_u32 = {};
  std::vector<uint64_t> data_u64 = {};
  std::vector<float>    data_f32 = {};
  std::vector<double>   data_f64 = {};

  // list length
  std::vector<uint8_t> ldata_u8 = {};
};

// Ply elements
struct element {
  // element content
  std::string                 name       = "";
  size_t                      count      = 0;
  std::vector<ply::property*> properties = {};

  // cleanup
  ~element();
};

// Ply model
struct model {
  // Format type
  enum struct format_t { ascii, binary_little_endian, binary_big_endian };

  // ply content
  format_t                   format   = format_t::binary_little_endian;
  std::vector<std::string>   comments = {};
  std::vector<ply::element*> elements = {};

  // cleanup
  ~model();
};

// Load and save ply
inline bool load_ply(
    const std::string& filename, ply::model* ply, std::string& error);
inline bool save_ply(
    const std::string& filename, ply::model* ply, std::string& error);

// Get ply properties
inline bool has_property(
    ply::model* ply, const std::string& element, const std::string& property);
inline ply::property* get_property(
    ply::model* ply, const std::string& element, const std::string& property);

inline bool get_value(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<float>& values);
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties, std::vector<vec4f>& values);
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties, std::vector<vec3f>& values);
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties, std::vector<vec4f>& values);
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    std::vector<frame3f>&              values);

inline bool get_lists(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<std::vector<int>>& lists);
inline bool get_list_sizes(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<byte>& sizes);
inline bool get_list_values(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<int>& values);

// Get ply properties for meshes
inline bool get_positions(ply::model* ply, std::vector<vec3f>& values);
inline bool get_normals(ply::model* ply, std::vector<vec3f>& values);
inline bool get_texcoords(
    ply::model* ply, std::vector<vec2f>& values, bool flipv = false);
inline bool get_colors(ply::model* ply, std::vector<vec3f>& values);
inline bool get_radius(ply::model* ply, std::vector<float>& values);
inline bool get_faces(ply::model* ply, std::vector<std::vector<int>>*& values);
inline bool get_lines(ply::model* ply, std::vector<vec2i>& values);
inline bool get_points(ply::model* ply, std::vector<int>& values);
inline bool get_triangles(ply::model* ply, std::vector<vec3i>& values);
inline bool get_quads(ply::model* ply, std::vector<vec4i>& values);
inline bool has_quads(ply::model* ply);

// Add ply properties
inline bool add_value(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<float>& values);
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties,
    const std::vector<vec2f>&         values);
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties,
    const std::vector<vec3f>&         values);
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties,
    const std::vector<vec4f>&         values);
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    const std::vector<frame3f>&        values);

inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<std::vector<int>>& values);
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<byte>& sizes,
    const std::vector<int>& values);
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<int>& values);
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec2i>& values);
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec3i>& values);
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec4i>& values);

// Add ply properties for meshes
inline bool add_positions(ply::model* ply, const std::vector<vec3f>& values);
inline bool add_normals(ply::model* ply, const std::vector<vec3f>& values);
inline bool add_texcoords(
    ply::model* ply, const std::vector<vec2f>& values, bool flipv = false);
inline bool add_colors(ply::model* ply, const std::vector<vec3f>& values);
inline bool add_radius(ply::model* ply, const std::vector<float>& values);
inline bool add_faces(
    ply::model* ply, const std::vector<std::vector<int>>& values);
inline bool add_faces(ply::model* ply, const std::vector<vec3i>& tvalues,
    const std::vector<vec4i>& qvalues);
inline bool add_triangles(ply::model* ply, const std::vector<vec3i>& values);
inline bool add_quads(ply::model* ply, const std::vector<vec4i>& values);
inline bool add_lines(ply::model* ply, const std::vector<vec2i>& values);
inline bool add_points(ply::model* ply, const std::vector<int>& values);

}  // namespace yocto::ply

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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::ply {

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
[[nodiscard]] inline bool parse_value(std::string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(std::string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
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
[[nodiscard]] inline bool parse_value(std::string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#ifdef __APPLE__
[[nodiscard]] inline bool parse_value(std::string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#endif

// Formats values to std::string
inline void format_value(std::string& str, const std::string& value) {
  str += value;
}
inline void format_value(std::string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(std::string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(std::string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(std::string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
inline void format_value(std::string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(std::string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(std::string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(std::string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
inline void format_value(std::string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_value(std::string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
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

template <typename T>
inline T swap_endian(T value) {
  // https://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
  static_assert(CHAR_BIT == 8, "CHAR_BIT != 8");
  union {
    T             value;
    unsigned char bytes[sizeof(T)];
  } source, dest;
  source.value = value;
  for (auto k = (size_t)0; k < sizeof(T); k++)
    dest.bytes[k] = source.bytes[sizeof(T) - k - 1];
  return dest.value;
}

template <typename T>
[[nodiscard]] inline bool read_value(FILE* fs, T& value, bool big_endian) {
  if (fread(&value, sizeof(value), 1, fs) != 1) return false;
  if (big_endian) value = swap_endian(value);
  return true;
}

template <typename T>
[[nodiscard]] inline bool write_value(
    FILE* fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  if (fwrite(&value, sizeof(value), 1, fs) != 1) return false;
  return true;
}

inline void remove_comment(std::string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

inline element::~element() {
  for (auto property : properties) delete property;
}
inline model::~model() {
  for (auto element : elements) delete element;
}

// Make ply
inline ply::element* add_property(ply::model* ply) {
  return ply->elements.emplace_back(new element{});
}
inline ply::property* add_property(ply::element* element) {
  return element->properties.emplace_back(new property{});
}

// Load ply
inline bool load_ply(
    const std::string& filename, ply::model* ply, std::string& error) {
  // ply type names
  static auto type_map = std::unordered_map<std::string, property::type_t>{
      {"char", property::type_t::i8}, {"short", property::type_t::i16},
      {"int", property::type_t::i32}, {"long", property::type_t::i64},
      {"uchar", property::type_t::u8}, {"ushort", property::type_t::u16},
      {"uint", property::type_t::u32}, {"ulong", property::type_t::u64},
      {"float", property::type_t::f32}, {"double", property::type_t::f64},
      {"int8", property::type_t::i8}, {"int16", property::type_t::i16},
      {"int32", property::type_t::i32}, {"int64", property::type_t::i64},
      {"uint8", property::type_t::u8}, {"uint16", property::type_t::u16},
      {"uint32", property::type_t::u32}, {"uint64", property::type_t::u64},
      {"float32", property::type_t::f32}, {"float64", property::type_t::f64}};

  // initialize data
  ply->comments.clear();
  ply->elements.clear();

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
  auto fs = fopen(filename.c_str(), "rb");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // read header ---------------------------------------------
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

    // check magic number
    if (first_line) {
      if (cmd != "ply") return parse_error();
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) return parse_error();
    } else if (cmd == "format") {
      auto fmt = ""s;
      if (!parse_value(str, fmt)) return parse_error();
      if (fmt == "ascii") {
        ply->format = model::format_t::ascii;
      } else if (fmt == "binary_little_endian") {
        ply->format = model::format_t::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        ply->format = model::format_t::binary_big_endian;
      } else {
        return parse_error();
      }
    } else if (cmd == "comment") {
      skip_whitespace(str);
      ply->comments.push_back(std::string{str});
    } else if (cmd == "obj_info") {
      skip_whitespace(str);
      // comment is the rest of the str
    } else if (cmd == "element") {
      auto elem = ply->elements.emplace_back(new element{});
      if (!parse_value(str, elem->name)) return parse_error();
      if (!parse_value(str, elem->count)) return parse_error();
    } else if (cmd == "property") {
      if (ply->elements.empty()) return parse_error();
      auto prop = ply->elements.back()->properties.emplace_back(new property{});
      auto tname = ""s;
      if (!parse_value(str, tname)) return parse_error();
      if (tname == "list") {
        prop->is_list = true;
        if (!parse_value(str, tname)) return parse_error();
        auto itype = type_map.at(tname);
        if (itype != property::type_t::u8) return parse_error();
        if (!parse_value(str, tname)) return parse_error();
        if (type_map.find(tname) == type_map.end()) return parse_error();
        prop->type = type_map.at(tname);
      } else {
        prop->is_list = false;
        if (type_map.find(tname) == type_map.end()) return parse_error();
        prop->type = type_map.at(tname);
      }
      if (!parse_value(str, prop->name)) return parse_error();
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      return parse_error();
    }
  }

  // check exit
  if (!end_header) return parse_error();

  // allocate data ---------------------------------
  for (auto element : ply->elements) {
    for (auto property : element->properties) {
      auto count = property->is_list ? element->count * 3 : element->count;
      switch (property->type) {
        case property::type_t::i8: property->data_i8.reserve(count); break;
        case property::type_t::i16: property->data_i16.reserve(count); break;
        case property::type_t::i32: property->data_i32.reserve(count); break;
        case property::type_t::i64: property->data_i64.reserve(count); break;
        case property::type_t::u8: property->data_u8.reserve(count); break;
        case property::type_t::u16: property->data_u16.reserve(count); break;
        case property::type_t::u32: property->data_u32.reserve(count); break;
        case property::type_t::u64: property->data_u64.reserve(count); break;
        case property::type_t::f32: property->data_f32.reserve(count); break;
        case property::type_t::f64: property->data_f64.reserve(count); break;
      }
      if (property->is_list) property->ldata_u8.reserve(element->count);
    }
  }

  // read data -------------------------------------
  if (ply->format == model::format_t::ascii) {
    char buffer[4096];
    for (auto elem : ply->elements) {
      for (auto idx = 0; idx < elem->count; idx++) {
        if (!fgets(buffer, sizeof(buffer), fs)) return read_error();
        auto str = std::string_view{buffer};
        for (auto prop : elem->properties) {
          if (prop->is_list) {
            if (!parse_value(str, prop->ldata_u8.emplace_back()))
              return parse_error();
          }
          auto vcount = prop->is_list ? prop->ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case property::type_t::i8:
                if (!parse_value(str, prop->data_i8.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::i16:
                if (!parse_value(str, prop->data_i16.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::i32:
                if (!parse_value(str, prop->data_i32.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::i64:
                if (!parse_value(str, prop->data_i64.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::u8:
                if (!parse_value(str, prop->data_u8.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::u16:
                if (!parse_value(str, prop->data_u16.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::u32:
                if (!parse_value(str, prop->data_u32.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::u64:
                if (!parse_value(str, prop->data_u64.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::f32:
                if (!parse_value(str, prop->data_f32.emplace_back()))
                  return parse_error();
                break;
              case property::type_t::f64:
                if (!parse_value(str, prop->data_f64.emplace_back()))
                  return parse_error();
                break;
            }
          }
        }
      }
    }
  } else {
    auto big_endian = ply->format == model::format_t::binary_big_endian;
    for (auto elem : ply->elements) {
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto prop : elem->properties) {
          if (prop->is_list) {
            if (!read_value(fs, prop->ldata_u8.emplace_back(), big_endian))
              return read_error();
          }
          auto vcount = prop->is_list ? prop->ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case property::type_t::i8:
                if (!read_value(fs, prop->data_i8.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::i16:
                if (!read_value(fs, prop->data_i16.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::i32:
                if (!read_value(fs, prop->data_i32.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::i64:
                if (!read_value(fs, prop->data_i64.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::u8:
                if (!read_value(fs, prop->data_u8.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::u16:
                if (!read_value(fs, prop->data_u16.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::u32:
                if (!read_value(fs, prop->data_u32.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::u64:
                if (!read_value(fs, prop->data_u64.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::f32:
                if (!read_value(fs, prop->data_f32.emplace_back(), big_endian))
                  return read_error();
                break;
              case property::type_t::f64:
                if (!read_value(fs, prop->data_f64.emplace_back(), big_endian))
                  return read_error();
                break;
            }
          }
        }
      }
    }
  }
  return true;
}

// Save ply
inline bool save_ply(
    const std::string& filename, ply::model* ply, std::string& error) {
  // ply type names
  static auto type_map = std::unordered_map<property::type_t, std::string>{
      {property::type_t::i8, "char"}, {property::type_t::i16, "short"},
      {property::type_t::i32, "int"}, {property::type_t::i64, "uint"},
      {property::type_t::u8, "uchar"}, {property::type_t::u16, "ushort"},
      {property::type_t::u32, "uint"}, {property::type_t::u64, "ulong"},
      {property::type_t::f32, "float"}, {property::type_t::f64, "double"}};
  static auto format_map = std::unordered_map<model::format_t, std::string>{
      {model::format_t::ascii, "ascii"},
      {model::format_t::binary_little_endian, "binary_little_endian"},
      {model::format_t::binary_big_endian, "binary_big_endian"}};

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
  auto fs = fopen(filename.c_str(), "wb");
  if (!fs) return open_error();
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // header
  if (!format_values(fs, "ply\n")) return write_error();
  if (!format_values(fs, "format {} 1.0\n", format_map.at(ply->format)))
    return write_error();
  if (!format_values(fs, "comment Written by Yocto/GL\n")) return write_error();
  if (!format_values(fs, "comment https://github.com/xelatihy/yocto-gl\n"))
    return write_error();
  for (auto& comment : ply->comments)
    if (!format_values(fs, "comment {}\n", comment)) return write_error();
  for (auto elem : ply->elements) {
    if (!format_values(
            fs, "element {} {}\n", elem->name, (uint64_t)elem->count))
      return write_error();
    for (auto prop : elem->properties) {
      if (prop->is_list) {
        if (!format_values(fs, "property list uchar {} {}\n",
                type_map[prop->type], prop->name))
          return write_error();
      } else {
        if (!format_values(
                fs, "property {} {}\n", type_map[prop->type], prop->name))
          return write_error();
      }
    }
  }

  if (!format_values(fs, "end_header\n")) return write_error();

  // properties
  if (ply->format == model::format_t::ascii) {
    for (auto elem : ply->elements) {
      auto cur = std::vector<size_t>(elem->properties.size(), 0);
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto pidx = 0; pidx < elem->properties.size(); pidx++) {
          auto prop = elem->properties[pidx];
          if (prop->is_list)
            if (!format_values(fs, "{} ", (int)prop->ldata_u8[idx]))
              return write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case property::type_t::i8:
                if (!format_values(fs, "{} ", prop->data_i8[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::i16:
                if (!format_values(fs, "{} ", prop->data_i16[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::i32:
                if (!format_values(fs, "{} ", prop->data_i32[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::i64:
                if (!format_values(fs, "{} ", prop->data_i64[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::u8:
                if (!format_values(fs, "{} ", prop->data_u8[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::u16:
                if (!format_values(fs, "{} ", prop->data_u16[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::u32:
                if (!format_values(fs, "{} ", prop->data_u32[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::u64:
                if (!format_values(fs, "{} ", prop->data_u64[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::f32:
                if (!format_values(fs, "{} ", prop->data_f32[cur[idx]++]))
                  return write_error();
                break;
              case property::type_t::f64:
                if (!format_values(fs, "{} ", prop->data_f64[cur[idx]++]))
                  return write_error();
                break;
            }
          }
          if (!format_values(fs, "\n")) return write_error();
        }
      }
    }
  } else {
    auto big_endian = ply->format == model::format_t::binary_big_endian;
    for (auto elem : ply->elements) {
      auto cur = std::vector<size_t>(elem->properties.size(), 0);
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto pidx = 0; pidx < elem->properties.size(); pidx++) {
          auto prop = elem->properties[pidx];
          if (prop->is_list)
            if (!write_value(fs, prop->ldata_u8[idx], big_endian))
              return write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case property::type_t::i8:
                if (!write_value(fs, prop->data_i8[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::i16:
                if (!write_value(fs, prop->data_i16[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::i32:
                if (!write_value(fs, prop->data_i32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::i64:
                if (!write_value(fs, prop->data_i64[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::u8:
                if (!write_value(fs, prop->data_u8[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::u16:
                if (!write_value(fs, prop->data_u16[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::u32:
                if (!write_value(fs, prop->data_u32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::u64:
                if (!write_value(fs, prop->data_u64[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::f32:
                if (!write_value(fs, prop->data_f32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case property::type_t::f64:
                if (!write_value(fs, prop->data_f64[cur[pidx]++], big_endian))
                  return write_error();
                break;
            }
          }
        }
      }
    }
  }

  return true;
}

// Get ply properties
inline bool has_property(
    ply::model* ply, const std::string& element, const std::string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return true;
    }
  }
  return false;
}
inline ply::property* get_property(
    ply::model* ply, const std::string& element, const std::string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
template <typename T, typename T1>
inline bool convert_property(
    const std::vector<T1>& prop, std::vector<T>& values) {
  values = std::vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return true;
}
template <typename T>
inline bool convert_property(ply::property* prop, std::vector<T>& values) {
  switch (prop->type) {
    case property::type_t::i8: return convert_property(prop->data_i8, values);
    case property::type_t::i16: return convert_property(prop->data_i16, values);
    case property::type_t::i32: return convert_property(prop->data_i32, values);
    case property::type_t::i64: return convert_property(prop->data_i64, values);
    case property::type_t::u8: return convert_property(prop->data_u8, values);
    case property::type_t::u16: return convert_property(prop->data_u16, values);
    case property::type_t::u32: return convert_property(prop->data_u32, values);
    case property::type_t::u64: return convert_property(prop->data_u64, values);
    case property::type_t::f32: return convert_property(prop->data_f32, values);
    case property::type_t::f64: return convert_property(prop->data_f64, values);
  }
  // return here to silence warnings
  std::runtime_error("should not have gotten here");
  return false;
}
inline bool get_value(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<float>& values) {
  values.clear();
  if (!has_property(ply, element, property)) return false;
  auto prop = get_property(ply, element, property);
  if (prop->is_list) return false;
  if (!convert_property(prop, values)) return false;
  return true;
}
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties, std::vector<vec2f>& values) {
  values.clear();
  auto x = std::vector<float>{}, y = std::vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  values = std::vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return true;
}
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties, std::vector<vec3f>& values) {
  values.clear();
  auto x = std::vector<float>{}, y = std::vector<float>{},
       z = std::vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  values = std::vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return true;
}
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties, std::vector<vec4f>& values) {
  values.clear();
  auto x = std::vector<float>{}, y = std::vector<float>{},
       z = std::vector<float>{}, w = std::vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  if (!get_value(ply, element, properties[3], w)) return false;
  values = std::vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w[i]};
  return true;
}
inline bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    std::vector<frame3f>&              values) {
  values.clear();
  auto coords = std::array<std::vector<float>, 12>{};
  for (auto idx = 0; idx < 12; idx++)
    if (!get_value(ply, element, properties[idx], coords[idx])) return false;
  values = std::vector<frame3f>(coords[0].size());
  for (auto i = (size_t)0; i < values.size(); i++) {
    for (auto c = 0; c < 12; c++) (&values[i].x.x)[c] = coords[c][i];
  }
  return true;
}
inline bool get_lists(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<std::vector<int>>& lists) {
  lists.clear();
  if (!has_property(ply, element, property)) return false;
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return false;
  auto& sizes  = prop->ldata_u8;
  auto  values = std::vector<int>{};
  if (!convert_property(prop, values)) return false;
  lists    = std::vector<std::vector<int>>(sizes.size());
  auto cur = (size_t)0;
  for (auto i = (size_t)0; i < lists.size(); i++) {
    lists[i].resize(sizes[i]);
    for (auto c = 0; c < sizes[i]; c++) {
      lists[i][c] = values[cur++];
    }
  }
  return true;
}
inline bool get_list_sizes(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<byte>& sizes) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  sizes = prop->ldata_u8;
  return true;
}
inline bool get_list_values(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<int>& values) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return convert_property<int>(prop, values);
}

inline std::vector<vec2f> flip_texcoord(const std::vector<vec2f>& texcoords) {
  auto flipped = texcoords;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
inline bool get_positions(ply::model* ply, std::vector<vec3f>& positions) {
  return get_values(ply, "vertex", {"x", "y", "z"}, positions);
}
inline bool get_normals(ply::model* ply, std::vector<vec3f>& normals) {
  return get_values(ply, "vertex", {"nx", "ny", "nz"}, normals);
}
inline bool get_texcoords(
    ply::model* ply, std::vector<vec2f>& texcoords, bool flipv) {
  if (has_property(ply, "vertex", "u")) {
    if (!get_values(ply, "vertex", {"u", "v"}, texcoords)) return false;
  } else {
    if (!get_values(ply, "vertex", {"s", "t"}, texcoords)) return false;
  }
  if (flipv) {
    for (auto& uv : texcoords) uv.y = 1 - uv.y;
  }
  return true;
}
inline bool get_colors(ply::model* ply, std::vector<vec3f>& colors) {
  return get_values(ply, "vertex", {"red", "green", "blue"}, colors);
}
inline bool get_radius(ply::model* ply, std::vector<float>& radius) {
  return get_value(ply, "vertex", "radius", radius);
}
inline bool get_faces(ply::model* ply, std::vector<std::vector<int>>& faces) {
  return get_lists(ply, "face", "vertex_indices", faces);
}
inline bool get_triangles(ply::model* ply, std::vector<vec3i>& triangles) {
  triangles.clear();
  auto indices = std::vector<int>{};
  auto sizes   = std::vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  triangles = std::vector<vec3i>{};
  triangles.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 2; c < size; c++) {
      triangles.push_back(
          {indices[cur + 0], indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return true;
}
inline bool get_quads(ply::model* ply, std::vector<vec4i>& quads) {
  quads.clear();
  auto indices = std::vector<int>{};
  auto sizes   = std::vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  quads = std::vector<vec4i>{};
  quads.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    if (size == 4) {
      quads.push_back({indices[cur + 0], indices[cur + 1], indices[cur + 2],
          indices[cur + 3]});
    } else {
      for (auto c = 2; c < size; c++) {
        quads.push_back({indices[cur + 0], indices[cur + c - 1],
            indices[cur + c], indices[cur + c]});
      }
    }
    cur += size;
  }
  return true;
}
inline bool get_lines(ply::model* ply, std::vector<vec2i>& lines) {
  auto indices = std::vector<int>{};
  auto sizes   = std::vector<uint8_t>{};
  if (!get_list_values(ply, "line", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "line", "vertex_indices", sizes)) return false;
  lines = std::vector<vec2i>{};
  lines.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 1; c < size; c++) {
      lines.push_back({indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return true;
}
inline bool get_points(ply::model* ply, std::vector<int>& values) {
  return get_list_values(ply, "point", "vertex_indices", values);
}
inline bool has_quads(ply::model* ply) {
  auto sizes = std::vector<uint8_t>{};
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
inline ply::element* add_element(
    ply::model* ply, const std::string& element_name, size_t count) {
  for (auto elem : ply->elements) {
    if (elem->name == element_name) return elem;
  }
  auto elem   = ply->elements.emplace_back(new element{});
  elem->name  = element_name;
  elem->count = count;
  return elem;
}
inline ply::property* add_property(ply::model* ply,
    const std::string& element_name, const std::string& property_name,
    size_t count, property::type_t type, bool is_list) {
  if (!add_element(ply, element_name, count)) return nullptr;
  for (auto elem : ply->elements) {
    if (elem->name != element_name) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property_name) return prop;
    }
    auto prop     = elem->properties.emplace_back(new property{});
    prop->name    = property_name;
    prop->type    = type;
    prop->is_list = is_list;
    return prop;
  }
  return nullptr;
}
template <typename T>
inline std::vector<T> make_vector(const T* value, size_t count, int stride) {
  auto ret = std::vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

inline bool add_values(ply::model* ply, const float* values, size_t count,
    const std::string& element, const std::string* properties, int nprops) {
  if (!values) return false;
  for (auto p = 0; p < nprops; p++) {
    if (!add_property(
            ply, element, properties[p], count, property::type_t::f32, false))
      return false;
    auto prop = get_property(ply, element, properties[p]);
    prop->data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_f32[i] = values[p + i * nprops];
  }
  return true;
}

inline bool add_value(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<float>& values) {
  if (values.empty()) return false;
  auto properties = std::vector{property};
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties,
    const std::vector<vec2f>&         values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties,
    const std::vector<vec3f>&         values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties,
    const std::vector<vec4f>&         values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}
inline bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    const std::vector<frame3f>&        values) {
  if (values.empty()) return false;
  return add_values(ply, (float*)values.data(), values.size(), element,
      properties.data(), properties.size());
}

inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<std::vector<int>>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, values.size(), property::type_t::i32, true))
    return false;
  auto prop = get_property(ply, element, property);
  prop->data_i32.reserve(values.size() * 4);
  prop->ldata_u8.reserve(values.size());
  for (auto& value : values) {
    prop->data_i32.insert(prop->data_i32.end(), value.begin(), value.end());
    prop->ldata_u8.push_back((uint8_t)value.size());
  }
  return true;
}
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<byte>& sizes,
    const std::vector<int>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, sizes.size(), property::type_t::i32, true))
    return false;
  auto prop      = get_property(ply, element, property);
  prop->data_i32 = values;
  prop->ldata_u8 = sizes;
  return true;
}
inline bool add_lists(ply::model* ply, const int* values, size_t count,
    int size, const std::string& element, const std::string& property) {
  if (!values) return false;
  if (!add_property(ply, element, property, count, property::type_t::i32, true))
    return false;
  auto prop = get_property(ply, element, property);
  prop->data_i32.assign(values, values + count * size);
  prop->ldata_u8.assign(count, size);
  return true;
}
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<int>& values) {
  if (values.empty()) return false;
  return add_lists(ply, values.data(), values.size(), 1, element, property);
}
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec2i>& values) {
  if (values.empty()) return false;
  return add_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
inline bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
}

// Add ply properties for meshes
inline bool add_positions(ply::model* ply, const std::vector<vec3f>& values) {
  return add_values(ply, "vertex", {"x", "y", "z"}, values);
}
inline bool add_normals(ply::model* ply, const std::vector<vec3f>& values) {
  return add_values(ply, "vertex", {"nx", "ny", "nz"}, values);
}
inline bool add_texcoords(
    ply::model* ply, const std::vector<vec2f>& values, bool flipv) {
  return add_values(
      ply, "vertex", {"u", "v"}, flipv ? flip_texcoord(values) : values);
}
inline bool add_colors(ply::model* ply, const std::vector<vec3f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue"}, values);
}
inline bool add_radius(ply::model* ply, const std::vector<float>& values) {
  return add_value(ply, "vertex", "radius", values);
}
inline bool add_faces(
    ply::model* ply, const std::vector<std::vector<int>>& values) {
  return add_lists(ply, "face", "vertex_indices", values);
}
inline bool add_faces(ply::model* ply, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads) {
  if (triangles.empty() && quads.empty()) return false;
  if (quads.empty()) {
    return add_lists(ply, "face", "vertex_indices", triangles);
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const vec4i& q) { return q.z != q.w; })) {
    return add_lists(ply, "face", "vertex_indices", quads);
  } else {
    auto sizes   = std::vector<uint8_t>();
    auto indices = std::vector<int>{};
    sizes.reserve(triangles.size() + quads.size());
    indices.reserve(triangles.size() * 3 + quads.size() * 4);
    for (auto& t : triangles) {
      sizes.push_back(3);
      indices.push_back(t.x);
      indices.push_back(t.y);
      indices.push_back(t.z);
    }
    for (auto& q : quads) {
      sizes.push_back(q.z == q.w ? 3 : 4);
      indices.push_back(q.x);
      indices.push_back(q.y);
      indices.push_back(q.z);
      if (q.z != q.w) indices.push_back(q.w);
    }
    return add_lists(ply, "face", "vertex_indices", sizes, indices);
  }
}
inline bool add_triangles(ply::model* ply, const std::vector<vec3i>& values) {
  return add_faces(ply, values, {});
}
inline bool add_quads(ply::model* ply, const std::vector<vec4i>& values) {
  return add_faces(ply, {}, values);
}
inline bool add_lines(ply::model* ply, const std::vector<vec2i>& values) {
  return add_lists(ply, "line", "vertex_indices", values);
}
inline bool add_points(ply::model* ply, const std::vector<int>& values) {
  return add_lists(ply, "point", "vertex_indices", values);
}

}  // namespace yocto::ply

#endif
