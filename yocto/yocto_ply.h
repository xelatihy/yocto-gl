//
// # Yocto/Ply: Tiny library for Ply parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Ply.
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

#ifndef _YOCTO_PLY_H_
#define _YOCTO_PLY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::ply {

// Using directives
using namespace yocto::math;

// Type of ply file. For best performance, choose binary_little_endian when
// writing ply files.
enum struct ply_format { ascii, binary_little_endian, binary_big_endian };

// Type of Ply data
enum struct ply_type { i8, i16, i32, i64, u8, u16, u32, u64, f32, f64 };

// Ply property
struct ply_property {
  // description
  string   name    = "";
  bool     is_list = false;
  ply_type type    = ply_type::f32;

  // data if property is loaded
  vector<int8_t>   data_i8  = {};
  vector<int16_t>  data_i16 = {};
  vector<int32_t>  data_i32 = {};
  vector<int64_t>  data_i64 = {};
  vector<uint8_t>  data_u8  = {};
  vector<uint16_t> data_u16 = {};
  vector<uint32_t> data_u32 = {};
  vector<uint64_t> data_u64 = {};
  vector<float>    data_f32 = {};
  vector<double>   data_f64 = {};

  // list length
  vector<uint8_t> ldata_u8 = {};
};

// Ply elements
struct ply_element {
  string                name       = "";
  size_t                count      = 0;
  vector<ply_property*> properties = {};
  ~ply_element();
};

// Ply model
struct ply_model {
  ply_format           format   = ply_format::binary_little_endian;
  vector<string>       comments = {};
  vector<ply_element*> elements = {};
  ~ply_model();
};

// Load and save ply
inline bool load_ply(const string& filename, ply_model* ply, string& error);
inline bool save_ply(const string& filename, ply_model* ply, string& error);

// Get ply properties
inline bool has_property(
    ply_model* ply, const string& element, const string& property);
inline ply_property* get_property(
    ply_model* ply, const string& element, const string& property);

inline vector<float> get_values(
    ply_model* ply, const string& element, const string& property);
inline vector<vec2f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2);
inline vector<vec3f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3);
inline vector<vec4f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      const string& property4);
inline vector<vec4f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      float property4);
inline vector<frame3f> get_values(
    ply_model* ply, const string& element, const array<string, 12>& properties);

inline vector<vector<int>> get_lists(
    ply_model* ply, const string& element, const string& property);
inline vector<byte> get_list_sizes(
    ply_model* ply, const string& element, const string& property);
inline vector<int> get_list_values(
    ply_model* ply, const string& element, const string& property);
inline vec2i get_list_minxmax(
    ply_model* ply, const string& element, const string& property);

// Get ply properties for meshes
inline vector<vec3f>       get_positions(ply_model* ply);
inline vector<vec3f>       get_normals(ply_model* ply);
inline vector<vec2f>       get_texcoords(ply_model* ply, bool flipv = false);
inline vector<vec3f>       get_colors(ply_model* ply);
inline vector<float>       get_radius(ply_model* ply);
inline vector<vector<int>> get_faces(ply_model* ply);
inline vector<vec2i>       get_lines(ply_model* ply);
inline vector<int>         get_points(ply_model* ply);
inline vector<vec3i>       get_triangles(ply_model* ply);
inline vector<vec4i>       get_quads(ply_model* ply);
inline bool                has_quads(ply_model* ply);

// Add ply properties
inline void add_values(ply_model* ply, const vector<float>& values,
    const string& element, const string& property);
inline void add_values(ply_model* ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2);
inline void add_values(ply_model* ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3);
inline void add_values(ply_model* ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4);
inline void add_values(ply_model* ply, const vector<frame3f>& values,
    const string& element, const array<string, 12>& properties);

inline void add_lists(ply_model* ply, const vector<vector<int>>& values,
    const string& element, const string& property);
inline void add_lists(ply_model* ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property);
inline void add_lists(ply_model* ply, const vector<int>& values,
    const string& element, const string& property);
inline void add_lists(ply_model* ply, const vector<vec2i>& values,
    const string& element, const string& property);
inline void add_lists(ply_model* ply, const vector<vec3i>& values,
    const string& element, const string& property);
inline void add_lists(ply_model* ply, const vector<vec4i>& values,
    const string& element, const string& property);

// Add ply properties for meshes
inline void add_positions(ply_model* ply, const vector<vec3f>& values);
inline void add_normals(ply_model* ply, const vector<vec3f>& values);
inline void add_texcoords(
    ply_model* ply, const vector<vec2f>& values, bool flipv = false);
inline void add_colors(ply_model* ply, const vector<vec3f>& values);
inline void add_radius(ply_model* ply, const vector<float>& values);
inline void add_faces(ply_model* ply, const vector<vector<int>>& values);
inline void add_faces(
    ply_model* ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
inline void add_triangles(ply_model* ply, const vector<vec3i>& values);
inline void add_quads(ply_model* ply, const vector<vec4i>& values);
inline void add_lines(ply_model* ply, const vector<vec2i>& values);
inline void add_points(ply_model* ply, const vector<int>& values);

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
[[nodiscard]] inline bool parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] inline bool parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#ifdef __APPLE__
[[nodiscard]] inline bool parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#endif

// Formats values to string
inline void format_value(string& str, const string& value) { str += value; }
inline void format_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
inline void format_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
inline void format_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
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

inline void remove_ply_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

inline ply_element::~ply_element() {
  for (auto property : properties) delete property;
}
inline ply_model::~ply_model() {
  for (auto element : elements) delete element;
}

// Make ply
inline ply_element* add_property(ply_model* ply) {
  return ply->elements.emplace_back(new ply_element{});
}
inline ply_property* add_property(ply_element* element) {
  return element->properties.emplace_back(new ply_property{});
}

// Load ply
inline bool load_ply(const string& filename, ply_model* ply, string& error) {
  // ply type names
  static auto type_map = unordered_map<string, ply_type>{{"char", ply_type::i8},
      {"short", ply_type::i16}, {"int", ply_type::i32}, {"long", ply_type::i64},
      {"uchar", ply_type::u8}, {"ushort", ply_type::u16},
      {"uint", ply_type::u32}, {"ulong", ply_type::u64},
      {"float", ply_type::f32}, {"double", ply_type::f64},
      {"int8", ply_type::i8}, {"int16", ply_type::i16},
      {"int32", ply_type::i32}, {"int64", ply_type::i64},
      {"uint8", ply_type::u8}, {"uint16", ply_type::u16},
      {"uint32", ply_type::u32}, {"uint64", ply_type::u64},
      {"float32", ply_type::f32}, {"float64", ply_type::f64}};

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
    auto str = string_view{buffer};
    remove_ply_comment(str);
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
        ply->format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        ply->format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        ply->format = ply_format::binary_big_endian;
      } else {
        return parse_error();
      }
    } else if (cmd == "comment") {
      skip_whitespace(str);
      ply->comments.push_back(string{str});
    } else if (cmd == "obj_info") {
      skip_whitespace(str);
      // comment is the rest of the str
    } else if (cmd == "element") {
      auto elem = ply->elements.emplace_back(new ply_element{});
      if (!parse_value(str, elem->name)) return parse_error();
      if (!parse_value(str, elem->count)) return parse_error();
    } else if (cmd == "property") {
      if (ply->elements.empty()) return parse_error();
      auto prop = ply->elements.back()->properties.emplace_back(
          new ply_property{});
      auto tname = ""s;
      if (!parse_value(str, tname)) return parse_error();
      if (tname == "list") {
        prop->is_list = true;
        if (!parse_value(str, tname)) return parse_error();
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8) return parse_error();
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
        case ply_type::i8: property->data_i8.reserve(count); break;
        case ply_type::i16: property->data_i16.reserve(count); break;
        case ply_type::i32: property->data_i32.reserve(count); break;
        case ply_type::i64: property->data_i64.reserve(count); break;
        case ply_type::u8: property->data_u8.reserve(count); break;
        case ply_type::u16: property->data_u16.reserve(count); break;
        case ply_type::u32: property->data_u32.reserve(count); break;
        case ply_type::u64: property->data_u64.reserve(count); break;
        case ply_type::f32: property->data_f32.reserve(count); break;
        case ply_type::f64: property->data_f64.reserve(count); break;
      }
      if (property->is_list) property->ldata_u8.reserve(element->count);
    }
  }

  // read data -------------------------------------
  if (ply->format == ply_format::ascii) {
    char buffer[4096];
    for (auto elem : ply->elements) {
      for (auto idx = 0; idx < elem->count; idx++) {
        if (!fgets(buffer, sizeof(buffer), fs)) return read_error();
        auto str = string_view{buffer};
        for (auto prop : elem->properties) {
          if (prop->is_list) {
            if (!parse_value(str, prop->ldata_u8.emplace_back()))
              return parse_error();
          }
          auto vcount = prop->is_list ? prop->ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!parse_value(str, prop->data_i8.emplace_back()))
                  return parse_error();
                break;
              case ply_type::i16:
                if (!parse_value(str, prop->data_i16.emplace_back()))
                  return parse_error();
                break;
              case ply_type::i32:
                if (!parse_value(str, prop->data_i32.emplace_back()))
                  return parse_error();
                break;
              case ply_type::i64:
                if (!parse_value(str, prop->data_i64.emplace_back()))
                  return parse_error();
                break;
              case ply_type::u8:
                if (!parse_value(str, prop->data_u8.emplace_back()))
                  return parse_error();
                break;
              case ply_type::u16:
                if (!parse_value(str, prop->data_u16.emplace_back()))
                  return parse_error();
                break;
              case ply_type::u32:
                if (!parse_value(str, prop->data_u32.emplace_back()))
                  return parse_error();
                break;
              case ply_type::u64:
                if (!parse_value(str, prop->data_u64.emplace_back()))
                  return parse_error();
                break;
              case ply_type::f32:
                if (!parse_value(str, prop->data_f32.emplace_back()))
                  return parse_error();
                break;
              case ply_type::f64:
                if (!parse_value(str, prop->data_f64.emplace_back()))
                  return parse_error();
                break;
            }
          }
        }
      }
    }
  } else {
    auto big_endian = ply->format == ply_format::binary_big_endian;
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
              case ply_type::i8:
                if (!read_value(fs, prop->data_i8.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::i16:
                if (!read_value(fs, prop->data_i16.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::i32:
                if (!read_value(fs, prop->data_i32.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::i64:
                if (!read_value(fs, prop->data_i64.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::u8:
                if (!read_value(fs, prop->data_u8.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::u16:
                if (!read_value(fs, prop->data_u16.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::u32:
                if (!read_value(fs, prop->data_u32.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::u64:
                if (!read_value(fs, prop->data_u64.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::f32:
                if (!read_value(fs, prop->data_f32.emplace_back(), big_endian))
                  return read_error();
                break;
              case ply_type::f64:
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
inline bool save_ply(const string& filename, ply_model* ply, string& error) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{{ply_type::i8, "char"},
      {ply_type::i16, "short"}, {ply_type::i32, "int"}, {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"}, {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"}, {ply_type::u64, "ulong"},
      {ply_type::f32, "float"}, {ply_type::f64, "double"}};
  static auto format_map = unordered_map<ply_format, string>{
      {ply_format::ascii, "ascii"},
      {ply_format::binary_little_endian, "binary_little_endian"},
      {ply_format::binary_big_endian, "binary_big_endian"}};

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
  if (ply->format == ply_format::ascii) {
    for (auto elem : ply->elements) {
      auto cur = vector<size_t>(elem->properties.size(), 0);
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto pidx = 0; pidx < elem->properties.size(); pidx++) {
          auto prop = elem->properties[pidx];
          if (prop->is_list)
            if (!format_values(fs, "{} ", (int)prop->ldata_u8[idx]))
              return write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!format_values(fs, "{} ", prop->data_i8[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::i16:
                if (!format_values(fs, "{} ", prop->data_i16[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::i32:
                if (!format_values(fs, "{} ", prop->data_i32[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::i64:
                if (!format_values(fs, "{} ", prop->data_i64[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::u8:
                if (!format_values(fs, "{} ", prop->data_u8[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::u16:
                if (!format_values(fs, "{} ", prop->data_u16[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::u32:
                if (!format_values(fs, "{} ", prop->data_u32[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::u64:
                if (!format_values(fs, "{} ", prop->data_u64[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::f32:
                if (!format_values(fs, "{} ", prop->data_f32[cur[idx]++]))
                  return write_error();
                break;
              case ply_type::f64:
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
    auto big_endian = ply->format == ply_format::binary_big_endian;
    for (auto elem : ply->elements) {
      auto cur = vector<size_t>(elem->properties.size(), 0);
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto pidx = 0; pidx < elem->properties.size(); pidx++) {
          auto prop = elem->properties[pidx];
          if (prop->is_list)
            if (!write_value(fs, prop->ldata_u8[idx], big_endian))
              return write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!write_value(fs, prop->data_i8[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::i16:
                if (!write_value(fs, prop->data_i16[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::i32:
                if (!write_value(fs, prop->data_i32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::i64:
                if (!write_value(fs, prop->data_i64[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::u8:
                if (!write_value(fs, prop->data_u8[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::u16:
                if (!write_value(fs, prop->data_u16[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::u32:
                if (!write_value(fs, prop->data_u32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::u64:
                if (!write_value(fs, prop->data_u64[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::f32:
                if (!write_value(fs, prop->data_f32[cur[pidx]++], big_endian))
                  return write_error();
                break;
              case ply_type::f64:
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
    ply_model* ply, const string& element, const string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return true;
    }
  }
  return false;
}
inline ply_property* get_property(
    ply_model* ply, const string& element, const string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
template <typename T, typename T1>
inline vector<T> convert_ply_property(const vector<T1>& prop) {
  auto values = vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return values;
}
template <typename T>
inline vector<T> convert_ply_property(ply_property* prop) {
  switch (prop->type) {
    case ply_type::i8: return convert_ply_property<T>(prop->data_i8);
    case ply_type::i16: return convert_ply_property<T>(prop->data_i16);
    case ply_type::i32: return convert_ply_property<T>(prop->data_i32);
    case ply_type::i64: return convert_ply_property<T>(prop->data_i64);
    case ply_type::u8: return convert_ply_property<T>(prop->data_u8);
    case ply_type::u16: return convert_ply_property<T>(prop->data_u16);
    case ply_type::u32: return convert_ply_property<T>(prop->data_u32);
    case ply_type::u64: return convert_ply_property<T>(prop->data_u64);
    case ply_type::f32: return convert_ply_property<T>(prop->data_f32);
    case ply_type::f64: return convert_ply_property<T>(prop->data_f64);
  }
  // return here to silence warnings
  std::runtime_error("should not have gotten here");
  return {};
}
inline vector<float> get_values(
    ply_model* ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (prop->is_list) return {};
  return convert_ply_property<float>(prop);
}
inline vector<vec2f> get_values(ply_model* ply, const string& element,
    const string& property1, const string& property2) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return values;
}
inline vector<vec3f> get_values(ply_model* ply, const string& element,
    const string& property1, const string& property2, const string& property3) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto z      = get_values(ply, element, property3);
  auto values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return values;
}
inline vector<vec4f> get_values(ply_model* ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    const string& property4) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto z      = get_values(ply, element, property3);
  auto w      = get_values(ply, element, property4);
  auto values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w[i]};
  return values;
}
inline vector<vec4f> get_values(ply_model* ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    float property4) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto z      = get_values(ply, element, property3);
  auto w      = property4;
  auto values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w};
  return values;
}
inline vector<frame3f> get_values(ply_model* ply, const string& element,
    const array<string, 12>& properties) {
  auto coords = array<vector<float>, 12>{};
  for (auto idx = 0; idx < 12; idx++)
    coords[idx] = get_values(ply, element, properties[idx]);
  auto values = vector<frame3f>(coords[0].size());
  for (auto i = (size_t)0; i < values.size(); i++) {
    for (auto c = 0; c < 12; c++) (&values[i].x.x)[c] = coords[c][i];
  }
  return values;
}
inline vector<vector<int>> get_lists(
    ply_model* ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  auto& sizes  = prop->ldata_u8;
  auto  values = convert_ply_property<int>(prop);
  auto  lists  = vector<vector<int>>(sizes.size());
  auto  cur    = (size_t)0;
  for (auto i = (size_t)0; i < lists.size(); i++) {
    lists[i].resize(sizes[i]);
    for (auto c = 0; c < sizes[i]; c++) {
      lists[i][c] = values[cur++];
    }
  }
  return lists;
}
inline vector<byte> get_list_sizes(
    ply_model* ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return prop->ldata_u8;
}
inline vector<int> get_list_values(
    ply_model* ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return convert_ply_property<int>(prop);
}

inline vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
inline vector<vec3f> get_positions(ply_model* ply) {
  return get_values(ply, "vertex", "x", "y", "z");
}
inline vector<vec3f> get_normals(ply_model* ply) {
  return get_values(ply, "vertex", "nx", "ny", "nz");
}
inline vector<vec2f> get_texcoords(ply_model* ply, bool flipv) {
  auto texcoord = has_property(ply, "vertex", "u")
                      ? get_values(ply, "vertex", "u", "v")
                      : get_values(ply, "vertex", "s", "t");
  return flipv ? flip_ply_texcoord(texcoord) : texcoord;
}
inline vector<vec3f> get_colors(ply_model* ply) {
  return get_values(ply, "vertex", "red", "green", "blue");
}
inline vector<float> get_radius(ply_model* ply) {
  return get_values(ply, "vertex", "radius");
}
inline vector<vector<int>> get_faces(ply_model* ply) {
  return get_lists(ply, "face", "vertex_indices");
}
inline vector<vec3i> get_triangles(ply_model* ply) {
  auto indices   = get_list_values(ply, "face", "vertex_indices");
  auto sizes     = get_list_sizes(ply, "face", "vertex_indices");
  auto triangles = vector<vec3i>{};
  triangles.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 2; c < size; c++) {
      triangles.push_back(
          {indices[cur + 0], indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return triangles;
}
inline vector<vec4i> get_quads(ply_model* ply) {
  auto indices = get_list_values(ply, "face", "vertex_indices");
  auto sizes   = get_list_sizes(ply, "face", "vertex_indices");
  auto quads   = vector<vec4i>{};
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
  return quads;
}
inline vector<vec2i> get_lines(ply_model* ply) {
  auto indices = get_list_values(ply, "str", "vertex_indices");
  auto sizes   = get_list_sizes(ply, "str", "vertex_indices");
  auto lines   = vector<vec2i>{};
  lines.reserve(sizes.size());
  auto cur = 0;
  for (auto size : sizes) {
    for (auto c = 1; c < size; c++) {
      lines.push_back({indices[cur + c - 1], indices[cur + c]});
    }
    cur += size;
  }
  return lines;
}
inline vector<int> get_points(ply_model* ply) {
  return get_list_values(ply, "point", "vertex_indices");
}
inline bool has_quads(ply_model* ply) {
  auto sizes = get_list_sizes(ply, "face", "vertex_indices");
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
inline void add_element(ply_model* ply, const string& element, size_t count) {
  for (auto elem : ply->elements) {
    if (elem->name == element) return;
  }
  auto elem   = ply->elements.emplace_back(new ply_element{});
  elem->name  = element;
  elem->count = count;
}
inline void add_property(ply_model* ply, const string& element,
    const string& property, size_t count, ply_type type, bool is_list) {
  add_element(ply, element, count);
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property)
        throw std::runtime_error("property already added");
    }
    auto prop     = elem->properties.emplace_back(new ply_property{});
    prop->name    = property;
    prop->type    = type;
    prop->is_list = is_list;
    return;
  }
}
template <typename T>
inline vector<T> make_ply_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

inline void add_values(ply_model* ply, const float* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (!values) return;
  for (auto p = 0; p < nprops; p++) {
    add_property(ply, element, properties[p], count, ply_type::f32, false);
    auto prop = get_property(ply, element, properties[p]);
    prop->data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_f32[i] = values[p + i * nprops];
  }
}

inline void add_values(ply_model* ply, const vector<float>& values,
    const string& element, const string& property) {
  auto properties = vector{property};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
inline void add_values(ply_model* ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2) {
  auto properties = vector{property1, property2};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
inline void add_values(ply_model* ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3) {
  auto properties = vector{property1, property2, property3};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
inline void add_values(ply_model* ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4) {
  auto properties = vector{property1, property2, property3, property4};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}
inline void add_values(ply_model* ply, const vector<frame3f>& values,
    const string& element, const array<string, 12>& properties) {
  add_values(ply, (float*)values.data(), values.size(), element,
      properties.data(), properties.size());
}

inline void add_lists(ply_model* ply, const vector<vector<int>>& values,
    const string& element, const string& property) {
  if (values.empty()) return;
  add_property(ply, element, property, values.size(), ply_type::i32, true);
  auto prop = get_property(ply, element, property);
  prop->data_i32.reserve(values.size() * 4);
  prop->ldata_u8.reserve(values.size());
  for (auto& value : values) {
    prop->data_i32.insert(prop->data_i32.end(), value.begin(), value.end());
    prop->ldata_u8.push_back((uint8_t)value.size());
  }
}
inline void add_lists(ply_model* ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property) {
  if (values.empty()) return;
  add_property(ply, element, property, sizes.size(), ply_type::i32, true);
  auto prop      = get_property(ply, element, property);
  prop->data_i32 = values;
  prop->ldata_u8 = sizes;
}
inline void add_lists(ply_model* ply, const int* values, size_t count, int size,
    const string& element, const string& property) {
  if (!values) return;
  add_property(ply, element, property, count, ply_type::i32, true);
  auto prop = get_property(ply, element, property);
  prop->data_i32.assign(values, values + count * size);
  prop->ldata_u8.assign(count, size);
}
inline void add_lists(ply_model* ply, const vector<int>& values,
    const string& element, const string& property) {
  return add_lists(ply, values.data(), values.size(), 1, element, property);
}
inline void add_lists(ply_model* ply, const vector<vec2i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
inline void add_lists(ply_model* ply, const vector<vec3i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
inline void add_lists(ply_model* ply, const vector<vec4i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
}

// Add ply properties for meshes
inline void add_positions(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, values, "vertex", "x", "y", "z");
}
inline void add_normals(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, values, "vertex", "nx", "ny", "nz");
}
inline void add_texcoords(
    ply_model* ply, const vector<vec2f>& values, bool flipv) {
  return add_values(
      ply, flipv ? flip_ply_texcoord(values) : values, "vertex", "u", "v");
}
inline void add_colors(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, values, "vertex", "red", "green", "blue");
}
inline void add_radius(ply_model* ply, const vector<float>& values) {
  return add_values(ply, values, "vertex", "radius");
}
inline void add_faces(ply_model* ply, const vector<vector<int>>& values) {
  return add_lists(ply, values, "face", "vertex_indices");
}
inline void add_faces(ply_model* ply, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
  if (triangles.empty() && quads.empty()) return;
  if (quads.empty()) {
    return add_lists(ply, triangles, "face", "vertex_indices");
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const vec4i& q) { return q.z != q.w; })) {
    return add_lists(ply, quads, "face", "vertex_indices");
  } else {
    auto sizes   = vector<uint8_t>();
    auto indices = vector<int>{};
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
    return add_lists(ply, sizes, indices, "face", "vertex_indices");
  }
}
inline void add_triangles(ply_model* ply, const vector<vec3i>& values) {
  return add_faces(ply, values, {});
}
inline void add_quads(ply_model* ply, const vector<vec4i>& values) {
  return add_faces(ply, {}, values);
}
inline void add_lines(ply_model* ply, const vector<vec2i>& values) {
  return add_lists(ply, values, "str", "vertex_indices");
}
inline void add_points(ply_model* ply, const vector<int>& values) {
  return add_lists(ply, values, "point", "vertex_indices");
}

}  // namespace yocto::ply

#endif
