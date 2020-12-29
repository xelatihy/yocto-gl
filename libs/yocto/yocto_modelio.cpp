//
// Implementation for Yocto/Ply.
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

#include "yocto_modelio.h"

#include <cstdio>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "yocto_color.h"
#include "yocto_commonio.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string_view;
using std::unordered_map;
using std::unordered_set;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const vec4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}
inline void format_value(string& str, const mat4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i != 0) str += " ";
    format_value(str, value[i]);
  }
}

inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

inline void remove_comment(
    string_view& str, char comment_char = '#', bool handle_quotes = false) {
  if (!handle_quotes) {
    while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
    auto cpy = str;
    while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
    str.remove_suffix(cpy.size());
  } else {
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
}

// Parse values from a string
inline bool parse_value(string_view& str, string_view& value) {
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
inline bool parse_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
inline bool parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
inline bool parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#ifdef __APPLE__
inline bool parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#endif
inline bool parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  if (!parse_value(str, valuei)) return false;
  value = (bool)valuei;
  return true;
}

inline bool parse_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
inline bool parse_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
inline bool parse_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
inline bool parse_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
inline bool parse_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

ply_element::~ply_element() {
  for (auto property : properties) delete property;
}
ply_model::~ply_model() {
  for (auto element : elements) delete element;
}

// Load ply
bool load_ply(const string& filename, ply_model* ply, string& error) {
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
  auto fs = open_file(filename, "rb");
  if (!fs) return open_error();

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // read header ---------------------------------------------
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    // str
    auto str = string_view{buffer.data()};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd.empty()) continue;

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
      ply->comments.emplace_back(str);
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
    auto buffer = array<char, 4096>{};
    for (auto elem : ply->elements) {
      for (auto idx = 0; idx < elem->count; idx++) {
        if (!read_line(fs, buffer)) return read_error();
        auto str = string_view{buffer.data()};
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
bool save_ply(const string& filename, const ply_model* ply, string& error) {
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
  auto fs = open_file(filename, "wb");
  if (!fs) return open_error();

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
        for (auto prop : elem->properties) {
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
bool has_property(
    ply_model* ply, const string& element, const string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return true;
    }
  }
  return false;
}
ply_property* get_property(
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
inline bool convert_property(const vector<T1>& prop, vector<T>& values) {
  values = vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return true;
}
template <typename T>
inline bool convert_property(ply_property* prop, vector<T>& values) {
  switch (prop->type) {
    case ply_type::i8: return convert_property(prop->data_i8, values);
    case ply_type::i16: return convert_property(prop->data_i16, values);
    case ply_type::i32: return convert_property(prop->data_i32, values);
    case ply_type::i64: return convert_property(prop->data_i64, values);
    case ply_type::u8: return convert_property(prop->data_u8, values);
    case ply_type::u16: return convert_property(prop->data_u16, values);
    case ply_type::u32: return convert_property(prop->data_u32, values);
    case ply_type::u64: return convert_property(prop->data_u64, values);
    case ply_type::f32: return convert_property(prop->data_f32, values);
    case ply_type::f64: return convert_property(prop->data_f64, values);
  }
  // return here to silence warnings
  throw std::runtime_error{"should not have gotten here"};
  return false;
}
bool get_value(ply_model* ply, const string& element, const string& property,
    vector<float>& values) {
  values.clear();
  if (!has_property(ply, element, property)) return false;
  auto prop = get_property(ply, element, property);
  if (prop->is_list) return false;
  if (!convert_property(prop, values)) return false;
  return true;
}
bool get_values(ply_model* ply, const string& element,
    const array<string, 2>& properties, vector<vec2f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return true;
}
bool get_values(ply_model* ply, const string& element,
    const array<string, 3>& properties, vector<vec3f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{}, z = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return true;
}
bool get_values(ply_model* ply, const string& element,
    const array<string, 4>& properties, vector<vec4f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{}, z = vector<float>{},
       w = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  if (!get_value(ply, element, properties[2], z)) return false;
  if (!get_value(ply, element, properties[3], w)) return false;
  values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w[i]};
  return true;
}
bool get_values(ply_model* ply, const string& element,
    const array<string, 12>& properties, vector<frame3f>& values) {
  values.clear();
  auto coords = array<vector<float>, 12>{};
  for (auto idx = 0; idx < 12; idx++)
    if (!get_value(ply, element, properties[idx], coords[idx])) return false;
  values = vector<frame3f>(coords[0].size());
  for (auto i = (size_t)0; i < values.size(); i++) {
    for (auto c = 0; c < 12; c++) (&values[i].x.x)[c] = coords[c][i];
  }
  return true;
}
bool get_lists(ply_model* ply, const string& element, const string& property,
    vector<vector<int>>& lists) {
  lists.clear();
  if (!has_property(ply, element, property)) return false;
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return false;
  auto& sizes  = prop->ldata_u8;
  auto  values = vector<int>{};
  if (!convert_property(prop, values)) return false;
  lists    = vector<vector<int>>(sizes.size());
  auto cur = (size_t)0;
  for (auto i = (size_t)0; i < lists.size(); i++) {
    lists[i].resize(sizes[i]);
    for (auto c = 0; c < sizes[i]; c++) {
      lists[i][c] = values[cur++];
    }
  }
  return true;
}
bool get_list_sizes(ply_model* ply, const string& element,
    const string& property, vector<byte>& sizes) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  sizes = prop->ldata_u8;
  return true;
}
bool get_list_values(ply_model* ply, const string& element,
    const string& property, vector<int>& values) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return convert_property<int>(prop, values);
}

inline vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoords) {
  auto flipped = texcoords;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
bool get_positions(ply_model* ply, vector<vec3f>& positions) {
  return get_values(ply, "vertex", {"x", "y", "z"}, positions);
}
bool get_normals(ply_model* ply, vector<vec3f>& normals) {
  return get_values(ply, "vertex", {"nx", "ny", "nz"}, normals);
}
bool get_texcoords(ply_model* ply, vector<vec2f>& texcoords, bool flipv) {
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
bool get_colors(ply_model* ply, vector<vec3f>& colors) {
  return get_values(ply, "vertex", {"red", "green", "blue"}, colors);
}
bool get_colors(ply_model* ply, vector<vec4f>& colors) {
  if (has_property(ply, "vertex", "alpha")) {
    return get_values(ply, "vertex", {"red", "green", "blue", "alpha"}, colors);
  } else {
    auto colors3 = vector<vec3f>{};
    if (!get_values(ply, "vertex", {"red", "green", "blue"}, colors3))
      return false;
    colors.resize(colors3.size());
    for (auto i = 0; i < colors.size(); i++)
      colors[i] = {colors3[i].x, colors3[i].y, colors3[i].z, 1};
    return true;
  }
}
bool get_radius(ply_model* ply, vector<float>& radius) {
  return get_value(ply, "vertex", "radius", radius);
}
bool get_faces(ply_model* ply, vector<vector<int>>& faces) {
  return get_lists(ply, "face", "vertex_indices", faces);
}
bool get_triangles(ply_model* ply, vector<vec3i>& triangles) {
  triangles.clear();
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  triangles = vector<vec3i>{};
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
bool get_quads(ply_model* ply, vector<vec4i>& quads) {
  quads.clear();
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "face", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  quads = vector<vec4i>{};
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
bool get_lines(ply_model* ply, vector<vec2i>& lines) {
  auto indices = vector<int>{};
  auto sizes   = vector<uint8_t>{};
  if (!get_list_values(ply, "line", "vertex_indices", indices)) return false;
  if (!get_list_sizes(ply, "line", "vertex_indices", sizes)) return false;
  lines = vector<vec2i>{};
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
bool get_points(ply_model* ply, vector<int>& values) {
  return get_list_values(ply, "point", "vertex_indices", values);
}
bool has_quads(ply_model* ply) {
  auto sizes = vector<uint8_t>{};
  if (!get_list_sizes(ply, "face", "vertex_indices", sizes)) return false;
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
inline ply_element* add_element(
    ply_model* ply, const string& element_name, size_t count) {
  for (auto elem : ply->elements) {
    if (elem->name == element_name) return elem;
  }
  auto elem   = ply->elements.emplace_back(new ply_element{});
  elem->name  = element_name;
  elem->count = count;
  return elem;
}
inline ply_property* add_property(ply_model* ply, const string& element_name,
    const string& property_name, size_t count, ply_type type, bool is_list) {
  if (add_element(ply, element_name, count) == nullptr) return nullptr;
  for (auto elem : ply->elements) {
    if (elem->name != element_name) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property_name) return prop;
    }
    auto prop     = elem->properties.emplace_back(new ply_property{});
    prop->name    = property_name;
    prop->type    = type;
    prop->is_list = is_list;
    return prop;
  }
  return nullptr;
}
template <typename T>
inline vector<T> make_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

inline bool add_values(ply_model* ply, const float* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (values == nullptr) return false;
  for (auto p = 0; p < nprops; p++) {
    if (add_property(ply, element, properties[p], count, ply_type::f32,
            false) == nullptr)
      return false;
    auto prop = get_property(ply, element, properties[p]);
    prop->data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_f32[i] = values[p + i * nprops];
  }
  return true;
}

inline bool add_values(ply_model* ply, const int* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (values == nullptr) return false;
  for (auto p = 0; p < nprops; p++) {
    if (add_property(ply, element, properties[p], count, ply_type::i32,
            false) == nullptr)
      return false;
    auto prop = get_property(ply, element, properties[p]);
    prop->data_i32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_i32[i] = values[p + i * nprops];
  }
  return true;
}

bool add_value(ply_model* ply, const string& element, const string& property,
    const vector<float>& values) {
  if (values.empty()) return false;
  auto properties = vector{property};
  return add_values(
      ply, values.data(), values.size(), element, properties.data(), 1);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 2>& properties, const vector<vec2f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 2);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 3>& properties, const vector<vec3f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 3);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 4>& properties, const vector<vec4f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 4);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 12>& properties, const vector<frame3f>& values) {
  if (values.empty()) return false;
  return add_values(ply, &values.front().x.x, values.size(), element,
      properties.data(), (int)properties.size());
}

bool add_value(ply_model* ply, const string& element, const string& property,
    const vector<int>& values) {
  if (values.empty()) return false;
  auto properties = vector{property};
  return add_values(
      ply, values.data(), values.size(), element, properties.data(), 1);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 2>& properties, const vector<vec2i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 2);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 3>& properties, const vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 3);
}
bool add_values(ply_model* ply, const string& element,
    const array<string, 4>& properties, const vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, &values.front().x, values.size(), element, properties.data(), 4);
}

bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vector<int>>& values) {
  if (values.empty()) return false;
  if (add_property(ply, element, property, values.size(), ply_type::i32,
          true) == nullptr)
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
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<byte>& sizes, const vector<int>& values) {
  if (values.empty()) return false;
  if (add_property(ply, element, property, sizes.size(), ply_type::i32, true) ==
      nullptr)
    return false;
  auto prop      = get_property(ply, element, property);
  prop->data_i32 = values;
  prop->ldata_u8 = sizes;
  return true;
}
bool add_lists(ply_model* ply, const int* values, size_t count, int size,
    const string& element, const string& property) {
  if (values == nullptr) return false;
  if (add_property(ply, element, property, count, ply_type::i32, true) ==
      nullptr)
    return false;
  auto prop = get_property(ply, element, property);
  prop->data_i32.assign(values, values + count * size);
  prop->ldata_u8.assign(count, (byte)size);
  return true;
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<int>& values) {
  if (values.empty()) return false;
  return add_lists(ply, values.data(), values.size(), 1, element, property);
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec2i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 2, element, property);
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 3, element, property);
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_lists(ply, &values.front().x, values.size(), 4, element, property);
}

// Add ply properties for meshes
bool add_positions(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"x", "y", "z"}, values);
}
bool add_normals(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"nx", "ny", "nz"}, values);
}
bool add_texcoords(ply_model* ply, const vector<vec2f>& values, bool flipv) {
  return add_values(
      ply, "vertex", {"u", "v"}, flipv ? flip_ply_texcoord(values) : values);
}
bool add_colors(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue"}, values);
}
bool add_colors(ply_model* ply, const vector<vec4f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue", "alpha"}, values);
}
bool add_radius(ply_model* ply, const vector<float>& values) {
  return add_value(ply, "vertex", "radius", values);
}
bool add_faces(ply_model* ply, const vector<vector<int>>& values) {
  return add_lists(ply, "face", "vertex_indices", values);
}
bool add_faces(ply_model* ply, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
  if (triangles.empty() && quads.empty()) return false;
  if (quads.empty()) {
    return add_lists(ply, "face", "vertex_indices", triangles);
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const vec4i& q) { return q.z != q.w; })) {
    return add_lists(ply, "face", "vertex_indices", quads);
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
    return add_lists(ply, "face", "vertex_indices", sizes, indices);
  }
}
bool add_triangles(ply_model* ply, const vector<vec3i>& values) {
  return add_faces(ply, values, {});
}
bool add_quads(ply_model* ply, const vector<vec4i>& values) {
  return add_faces(ply, {}, values);
}
bool add_lines(ply_model* ply, const vector<vec2i>& values) {
  return add_lists(ply, "line", "vertex_indices", values);
}
bool add_points(ply_model* ply, const vector<int>& values) {
  return add_lists(ply, "point", "vertex_indices", values);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

inline bool parse_value(string_view& str, obj_vertex& value) {
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
inline bool parse_value(string_view& str, obj_texture& info) {
  // initialize
  info = obj_texture();

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
inline bool load_mtl(const string& filename, obj_scene* obj, string& error) {
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
  auto fs = open_file(filename, "rt");
  if (!fs) return open_error();

  // init parsing
  add_material(obj);

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    // str
    auto str = string_view{buffer.data()};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd.empty()) continue;

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
    } else if (cmd == "Psh") {
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
    } else if (cmd == "Pss") {
      if (!parse_value(str, material->pbr_translucency)) return parse_error();
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
    } else if (cmd == "Pthin") {
      if (!parse_value(str, material->pbr_thin)) return parse_error();
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
    } else if (cmd == "map_Psh") {
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
    } else if (cmd == "map_Pss") {
      if (!parse_value(str, material->pbr_translucency_tex))
        return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pvs") {
      if (!parse_value(str, material->pbr_volscattering_tex))
        return parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pdisp") {
      if (!parse_value(str, material->pbr_displacement_tex))
        return parse_error();
    } else if (cmd == "map_Pnorm") {
      if (!parse_value(str, material->pbr_normal_tex)) return parse_error();
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
      material->pbr_specular = max(material->specular) != 0 ? 1 : 0;
    }
    material->pbr_bump_tex         = material->bump_tex;
    material->pbr_normal_tex       = material->normal_tex;
    material->pbr_displacement_tex = material->displacement_tex;
  }

  return true;
}

// Read obj
inline bool load_objx(const string& filename, obj_scene* obj, string& error) {
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
  auto fs = open_file(filename, "rt");
  if (!fs) return open_error();

  // shape map for instances
  auto shape_map = unordered_map<string, vector<obj_shape*>>{};
  for (auto shape : obj->shapes) {
    shape_map[shape->name].push_back(shape);
  }

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    // str
    auto str = string_view{buffer.data()};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd.empty()) continue;

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

obj_scene::~obj_scene() {
  for (auto shape : shapes) delete shape;
  for (auto material : materials) delete material;
  for (auto camera : cameras) delete camera;
  for (auto environment : environments) delete environment;
}

// Make obj
obj_camera* add_camera(obj_scene* obj) {
  return obj->cameras.emplace_back(new obj_camera{});
}
obj_material* add_material(obj_scene* obj) {
  return obj->materials.emplace_back(new obj_material{});
}
obj_environment* add_environment(obj_scene* obj) {
  return obj->environments.emplace_back(new obj_environment{});
}
obj_shape* add_shape(obj_scene* obj) {
  return obj->shapes.emplace_back(new obj_shape{});
}

// Read obj
bool load_obj(const string& filename, obj_scene* obj, string& error,
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
  auto fs = open_file(filename, "rt");
  if (!fs) return open_error();

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
  obj->~obj_scene();
  obj->cameras.clear();
  obj->environments.clear();
  obj->shapes.clear();
  obj->materials.clear();

  // initialize load
  obj->shapes.emplace_back(new obj_shape{});
  auto empty_material = (obj_material*)nullptr;

  // read the file str by str
  auto buffer = array<char, 4096>{};
  while (read_line(fs, buffer)) {
    // str
    auto str = string_view{buffer.data()};
    remove_comment(str);
    skip_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_value(str, cmd)) return parse_error();
    if (cmd.empty()) continue;

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
        if (shape->materials.back() != mname) {
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
        if (mname.empty() && empty_material == nullptr) {
          empty_material   = obj->materials.emplace_back(new obj_material{});
          material_map[""] = empty_material;
        }
        auto mat_idx = -1;
        for (auto midx = 0; midx < shape->materials.size(); midx++)
          if (shape->materials[midx] == mname) mat_idx = midx;
        if (mat_idx < 0) {
          shape->materials.push_back(mname);
          mat_idx = (int)shape->materials.size() - 1;
        }
        element.material = (uint8_t)mat_idx;
      }
      // parse vertices
      skip_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        if (!parse_value(str, vert)) return parse_error();
        if (vert.position == 0) break;
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
        if (!load_mtl(path_join(path_dirname(filename), mtllib), obj, error))
          return dependent_error();
        for (auto material : obj->materials)
          material_map[material->name] = material;
      }
    } else {
      // unused
    }
  }

  // fix empty material
  if (empty_material != nullptr) {
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
      if (vertex.position != 0 && ipositions[vertex.position] == 0) {
        shape->positions.push_back(opositions[vertex.position - 1]);
        ipositions[vertex.position] = (int)shape->positions.size();
      }
      if (vertex.normal != 0 && inormals[vertex.normal] == 0) {
        shape->normals.push_back(onormals[vertex.normal - 1]);
        inormals[vertex.normal] = (int)shape->normals.size();
      }
      if (vertex.texcoord != 0 && itexcoords[vertex.texcoord] == 0) {
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
  auto extfilename = replace_extension(filename, ".objx");
  if (path_exists(extfilename)) {
    if (!load_objx(extfilename, obj, error)) return dependent_error();
  }

  return true;
}

// Format values
inline void format_value(string& str, const obj_texture& value) {
  str += value.path.empty() ? "" : value.path;
}
inline void format_value(string& str, const obj_vertex& value) {
  format_value(str, value.position);
  if (value.texcoord != 0) {
    str += "/";
    format_value(str, value.texcoord);
    if (value.normal != 0) {
      str += "/";
      format_value(str, value.normal);
    }
  } else if (value.normal != 0) {
    str += "//";
    format_value(str, value.normal);
  }
}

// Save obj
inline bool save_mtl(
    const string& filename, const obj_scene* obj, string& error) {
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
  auto fs = open_file(filename, "wt");
  if (!fs) return open_error();

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
      if (material->pbr_specular != 0)
        if (!format_values(fs, "Ps {}\n", material->pbr_specular))
          return write_error();
      if (material->pbr_roughness != 0)
        if (!format_values(fs, "Pr {}\n", material->pbr_roughness))
          return write_error();
      if (material->pbr_metallic != 0)
        if (!format_values(fs, "Pm {}\n", material->pbr_metallic))
          return write_error();
      if (material->pbr_sheen != 0)
        if (!format_values(fs, "Psh {}\n", material->pbr_sheen))
          return write_error();
      if (material->pbr_transmission != 0)
        if (!format_values(fs, "Pt {}\n", material->pbr_transmission))
          return write_error();
      if (material->pbr_translucency != 0)
        if (!format_values(fs, "Pss {}\n", material->pbr_translucency))
          return write_error();
      if (material->pbr_coat != 0)
        if (!format_values(fs, "Pc {}\n", material->pbr_coat))
          return write_error();
      if (material->pbr_coatroughness != 0)
        if (!format_values(fs, "Pcr {}\n", material->pbr_coatroughness))
          return write_error();
      if (material->pbr_volscattering != zero3f)
        if (!format_values(fs, "Pvs {}\n", material->pbr_volscattering))
          return write_error();
      if (material->pbr_volanisotropy != 0)
        if (!format_values(fs, "Pvg {}\n", material->pbr_volanisotropy))
          return write_error();
      if (material->pbr_volscale != 0)
        if (!format_values(fs, "Pvr {}\n", material->pbr_volscale))
          return write_error();
      if (!material->pbr_emission_tex.path.empty())
        if (!format_values(fs, "map_Pe {}\n", material->pbr_emission_tex))
          return write_error();
      if (!material->pbr_base_tex.path.empty())
        if (!format_values(fs, "map_Pb {}\n", material->pbr_base_tex))
          return write_error();
      if (!material->pbr_specular_tex.path.empty())
        if (!format_values(fs, "map_Ps {}\n", material->pbr_specular_tex))
          return write_error();
      if (!material->pbr_roughness_tex.path.empty())
        if (!format_values(fs, "map_Pr {}\n", material->pbr_roughness_tex))
          return write_error();
      if (!material->pbr_metallic_tex.path.empty())
        if (!format_values(fs, "map_Pm {}\n", material->pbr_metallic_tex))
          return write_error();
      if (!material->pbr_sheen_tex.path.empty())
        if (!format_values(fs, "map_Psh {}\n", material->pbr_sheen_tex))
          return write_error();
      if (!material->pbr_transmission_tex.path.empty())
        if (!format_values(fs, "map_Pt {}\n", material->pbr_transmission_tex))
          return write_error();
      if (!material->pbr_translucency_tex.path.empty())
        if (!format_values(fs, "map_Pss {}\n", material->pbr_translucency_tex))
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
        if (!format_values(fs, "map_Pbump {}\n", material->pbr_bump_tex))
          return write_error();
      if (!material->displacement_tex.path.empty())
        if (!format_values(
                fs, "map_Pdisp {}\n", material->pbr_displacement_tex))
          return write_error();
      if (!material->normal_tex.path.empty())
        if (!format_values(fs, "map_Pnorm {}\n", material->pbr_normal_tex))
          return write_error();
    }
    if (!format_values(fs, "\n")) return write_error();
  }
  return true;
}

// Save obj
inline bool save_objx(
    const string& filename, const obj_scene* obj, string& error) {
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
  auto fs = open_file(filename, "wt");
  if (!fs) return open_error();

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
bool save_obj(const string& filename, const obj_scene* obj, string& error) {
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
  auto fs = open_file(filename, "wt");
  if (!fs) return open_error();

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
            replace_extension(path_filename(filename), ".mtl")))
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
                  fs, "usemtl {}\n", shape->materials[element.material]))
            return write_error();
          cur_material = element.material;
        }
        if (!format_values(fs, "{}", label)) return write_error();
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape->vertices[cur_vertex++];
          if (vert.position != 0) vert.position += vert_size.position;
          if (vert.normal != 0) vert.normal += vert_size.normal;
          if (vert.texcoord != 0) vert.texcoord += vert_size.texcoord;
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
    if (!save_mtl(replace_extension(filename, ".mtl"), obj, error))
      return dependent_error();
  }

  // save objx
  if (!obj->cameras.empty() || !obj->environments.empty() ||
      std::any_of(obj->shapes.begin(), obj->shapes.end(),
          [](auto shape) { return !shape->instances.empty(); })) {
    if (!save_objx(replace_extension(filename, ".objx"), obj, error))
      return dependent_error();
  }

  // done
  return true;
}

// Get obj vertices
void get_vertices(const obj_shape* shape, vector<vec3f>& positions,
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
    if (!shape->positions.empty() && vert.position != 0)
      positions.push_back(shape->positions[vert.position - 1]);
    if (!shape->normals.empty() && vert.normal != 0)
      normals.push_back(shape->normals[vert.normal - 1]);
    if (!shape->texcoords.empty() && vert.texcoord != 0)
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
void get_triangles(const obj_shape* shape, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<string>& materials, vector<int>& ematerials, bool flipv) {
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
void get_quads(const obj_shape* shape, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<string>& materials, vector<int>& ematerials, bool flipv) {
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
void get_lines(const obj_shape* shape, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<string>& materials, vector<int>& ematerials, bool flipv) {
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
void get_points(const obj_shape* shape, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<string>& materials, vector<int>& ematerials, bool flipv) {
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
void get_fvquads(const obj_shape* shape, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<string>& materials, vector<int>& ematerials, bool flipv) {
  if (shape->faces.empty()) return;
  positions = shape->positions;
  normals   = shape->normals;
  texcoords = flipv ? flip_obj_texcoord(shape->texcoords) : shape->texcoords;
  materials = shape->materials;
  if (shape->vertices[0].position != 0) quadspos.reserve(shape->faces.size());
  if (shape->vertices[0].normal != 0) quadsnorm.reserve(shape->faces.size());
  if (shape->vertices[0].texcoord != 0)
    quadstexcoord.reserve(shape->faces.size());
  if (!materials.empty()) ematerials.reserve(shape->faces.size());
  auto cur = 0;
  for (auto& face : shape->faces) {
    if (face.size == 4) {
      if (shape->vertices[0].position != 0)
        quadspos.push_back({shape->vertices[cur + 0].position - 1,
            shape->vertices[cur + 1].position - 1,
            shape->vertices[cur + 2].position - 1,
            shape->vertices[cur + 3].position - 1});
      if (shape->vertices[0].normal != 0)
        quadsnorm.push_back({shape->vertices[cur + 0].normal - 1,
            shape->vertices[cur + 1].normal - 1,
            shape->vertices[cur + 2].normal - 1,
            shape->vertices[cur + 3].normal - 1});
      if (shape->vertices[0].texcoord != 0)
        quadstexcoord.push_back({shape->vertices[cur + 0].texcoord - 1,
            shape->vertices[cur + 1].texcoord - 1,
            shape->vertices[cur + 2].texcoord - 1,
            shape->vertices[cur + 3].texcoord - 1});
      if (!materials.empty()) ematerials.push_back(face.material);
    } else {
      for (auto c = 2; c < face.size; c++) {
        if (shape->vertices[0].position != 0)
          quadspos.push_back({shape->vertices[cur + 0].position - 1,
              shape->vertices[cur + c - 1].position - 1,
              shape->vertices[cur + c].position - 1,
              shape->vertices[cur + c].position - 1});
        if (shape->vertices[0].normal != 0)
          quadsnorm.push_back({shape->vertices[cur + 0].normal - 1,
              shape->vertices[cur + c - 1].normal - 1,
              shape->vertices[cur + c].normal - 1,
              shape->vertices[cur + c].normal - 1});
        if (shape->vertices[0].texcoord != 0)
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

bool has_quads(obj_shape* shape) {
  for (auto& face : shape->faces)
    if (face.size == 4) return true;
  return false;
}

// Get obj vertices
void get_vertices(const obj_shape* shape, int material,
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
    if (!shape->positions.empty() && vert.position != 0)
      positions.push_back(shape->positions[vert.position - 1]);
    if (!shape->normals.empty() && vert.normal != 0)
      normals.push_back(shape->normals[vert.normal - 1]);
    if (!shape->texcoords.empty() && vert.texcoord != 0)
      texcoords.push_back(shape->texcoords[vert.texcoord - 1]);
  }
  if (flipv) {
    for (auto& texcoord : texcoords) texcoord.y = 1 - texcoord.y;
  }
}

// Get obj shape
void get_triangles(const obj_shape* shape, int material,
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
void get_quads(const obj_shape* shape, int material, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flipv) {
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
void get_lines(const obj_shape* shape, int material, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flipv) {
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
void get_points(const obj_shape* shape, int material, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flipv) {
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

// Add obj shape
void set_triangles(obj_shape* shape, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
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
void set_quads(obj_shape* shape, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
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
void set_lines(obj_shape* shape, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
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
void set_points(obj_shape* shape, const vector<int>& points,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
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
void set_fvquads(obj_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape->positions = positions;
  shape->normals   = normals;
  shape->texcoords = flipv ? flip_obj_texcoord(texcoords) : texcoords;
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
void set_materials(obj_shape* shape, const vector<string>& materials) {
  shape->materials = materials;
}
void set_instances(obj_shape* shape, const vector<frame3f>& instances) {
  shape->instances = instances;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR StL
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::vec3f> {
  size_t operator()(const yocto::vec3f& v) const {
    const std::hash<float> hasher = std::hash<float>();
    auto                   h      = (size_t)0;
    h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// STL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// cleanup
stl_model::~stl_model() {
  for (auto shape : shapes) delete shape;
}

// Load/save stl
bool load_stl(const string& filename, stl_model* stl, string& error,
    bool unique_vertices) {
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

  for (auto shape : stl->shapes) delete shape;
  stl->shapes.clear();

  // open file
  auto fs = open_file(filename, "rb");
  if (!fs) return open_error();

  // assume it is binary and read hader
  auto header = array<char, 80>{};
  if (!read_value(fs, header)) return read_error();

  // check if binary
  auto binary = header[0] != 's' || header[1] != 'o' || header[2] != 'l' ||
                header[3] != 'i' || header[4] != 'd';

  // check size in case the binary had a bad header
  if (!binary) {
    auto ntriangles = (uint32_t)0;
    if (!read_value(fs, ntriangles)) return read_error();
    fseek(fs.fs, 0, SEEK_SET);
    fseek(fs.fs, 0, SEEK_END);
    auto length = ftell(fs.fs);
    fseek(fs.fs, 0, SEEK_SET);
    auto size = 80 + 4 + (4 * 12 + 2) * (size_t)ntriangles;
    binary    = length == size;
  }

  // close file
  close_file(fs);

  // switch on type
  if (binary) {
      // open file
      auto fs = open_file(filename, "rb");
      if (!fs) return open_error();
      
    // skip header
    auto header = array<char, 80>{};
    if (!read_value(fs, header)) return read_error();

    // read shapes until the end
    auto ntriangles = (uint32_t)0;
    while (read_value(fs, ntriangles)) {
      // append shape
      auto shape = stl->shapes.emplace_back(new stl_shape{});

      // resize buffers
      shape->fnormals.resize(ntriangles);
      shape->triangles.resize(ntriangles);
      shape->positions.resize(ntriangles * 3);

      // read all data
      for (auto triangle_id = 0; triangle_id < ntriangles; triangle_id++) {
        // read triangle data
        if (!read_value(fs, shape->fnormals[triangle_id])) return read_error();
        if (!read_value(fs, shape->positions[triangle_id * 3 + 0]))
          return read_error();
        if (!read_value(fs, shape->positions[triangle_id * 3 + 1]))
          return read_error();
        if (!read_value(fs, shape->positions[triangle_id * 3 + 2]))
          return read_error();
        shape->triangles[triangle_id] = {
            triangle_id * 3 + 0, triangle_id * 3 + 1, triangle_id * 3 + 2};
        // read unused attrobute count
        auto attribute_count = (uint16_t)0;
        if (!read_value(fs, attribute_count)) return read_error();
        // if (attribute_count != 0) return parse_error();
      }
    }

    // check if read at least one
    if (stl->shapes.empty()) return read_error();
  } else {
    // if ascii, re-open the file as text
    auto fs = open_file(filename, "rt");
    if (!fs) return open_error();

    // parse state
    auto in_solid = false, in_facet = false, in_loop = false;
    // raed all lines
    auto buffer = array<char, 4096>{};
    while (read_line(fs, buffer)) {
      // str
      auto str = string_view{buffer.data()};
      remove_comment(str);
      skip_whitespace(str);
      if (str.empty()) continue;

      // get command
      auto cmd = ""s;
      if (!parse_value(str, cmd)) return parse_error();
      if (cmd.empty()) continue;

      // switch over command
      if (cmd == "solid") {
        if (in_solid) return parse_error();
        in_solid = true;
        stl->shapes.emplace_back(new stl_shape{});
      } else if (cmd == "endsolid") {
        if (!in_solid) return parse_error();
        in_solid = false;
      } else if (cmd == "facet") {
        if (!in_solid || in_facet) return parse_error();
        in_facet = true;
        // next command
        if (!parse_value(str, cmd)) return parse_error();
        if (cmd != "normal") return parse_error();
        // vertex normal
        if (!parse_value(str, stl->shapes.back()->fnormals.emplace_back()))
          return parse_error();
      } else if (cmd == "endfacet") {
        if (!in_solid || !in_facet || in_loop) return parse_error();
        in_facet = false;
        // check that it was a triangle
        auto last_pos = (int)stl->shapes.back()->positions.size() - 3;
        if (stl->shapes.back()->triangles.empty() && last_pos != 0)
            return parse_error();
        if (!stl->shapes.back()->triangles.empty() && last_pos != stl->shapes.back()->triangles.back().z + 1)
          return parse_error();
        // add triangle
        stl->shapes.back()->triangles.push_back(
            {last_pos + 0, last_pos + 1, last_pos + 2});
      } else if (cmd == "outer") {
        if (!in_solid || !in_facet || in_loop) return parse_error();
        in_loop = true;
        // next command
        if (!parse_value(str, cmd)) return parse_error();
        if (cmd != "loop") return parse_error();
      } else if (cmd == "endloop") {
        if (!in_solid || !in_facet || !in_loop) return parse_error();
        in_loop = false;
      } else if (cmd == "vertex") {
        // vertex position
        if (!parse_value(str, stl->shapes.back()->positions.emplace_back()))
          return parse_error();
      } else {
        return parse_error();
      }
    }
  }

  // make unique vertices
  if (unique_vertices) {
    for (auto& shape : stl->shapes) {
      auto vertex_map       = unordered_map<vec3f, int>{};
      auto unique_positions = vector<vec3f>{};
      for (auto& triangle : shape->triangles) {
        for (auto& vertex_id : triangle) {
          auto vertex_it = vertex_map.find(shape->positions[vertex_id]);
          if (vertex_it == vertex_map.end()) {
            auto new_vertex_id = (int)unique_positions.size();
            unique_positions.push_back(shape->positions[vertex_id]);
            vertex_map.insert(
                vertex_it, {unique_positions.back(), new_vertex_id});
            vertex_id = new_vertex_id;
          } else {
            vertex_id = vertex_it->second;
          }
        }
      }
      std::swap(unique_positions, shape->positions);
    }
  }

  // done
  return true;
}

bool save_stl(
    const string& filename, const stl_model* stl, string& error, bool ascii) {
  // error helpers
  auto open_error = [filename, &error]() {
    error = filename + ": file not found";
    return false;
  };
  auto write_error = [filename, &error]() {
    error = filename + ": write error";
    return false;
  };

  // helper
  auto triangle_normal = [](const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return normalize(cross(p1 - p0, p2 - p0));
  };

  // open file
  auto fs = open_file(filename, ascii ? "wt" : "wb");
  if (!fs) return open_error();

  // switch on format
  if (!ascii) {
    // header
    auto header = array<char, 80>{0};
    snprintf(header.data(), header.size(), "Binary STL - Written by Yocto/GL");
    if (!write_value(fs, header)) return write_error();

    // write shapes
    for (auto& shape : stl->shapes) {
      auto ntriangles = (uint32_t)shape->triangles.size();
      if (!write_value(fs, ntriangles)) return write_error();
      for (auto triangle_idx = 0; triangle_idx < shape->triangles.size();
           triangle_idx++) {
        auto& triangle = shape->triangles[triangle_idx];
        auto  fnormal  = !shape->fnormals.empty()
                           ? shape->fnormals[triangle_idx]
                           : triangle_normal(shape->positions[triangle.x],
                                 shape->positions[triangle.y],
                                 shape->positions[triangle.z]);
        if (!write_value(fs, fnormal)) return write_error();
        if (!write_value(fs, shape->positions[triangle.x]))
          return write_error();
        if (!write_value(fs, shape->positions[triangle.y]))
          return write_error();
        if (!write_value(fs, shape->positions[triangle.z]))
          return write_error();
        auto attribute_count = (uint16_t)0;
        if (!write_value(fs, attribute_count)) return write_error();
      }
    }
  } else {
    for (auto& shape : stl->shapes) {
      if (!format_values(fs, "solid \n")) return write_error();
      for (auto triangle_idx = 0; triangle_idx < shape->triangles.size();
           triangle_idx++) {
        auto& triangle = shape->triangles[triangle_idx];
        auto  fnormal  = !shape->fnormals.empty()
                           ? shape->fnormals[triangle_idx]
                           : triangle_normal(shape->positions[triangle.x],
                                 shape->positions[triangle.y],
                                 shape->positions[triangle.z]);
        if (!format_values(fs, "facet normal {}\n", fnormal))
          return write_error();
        if (!format_values(fs, "outer loop\n")) return write_error();
        if (!format_values(fs, "vertex {}\n", shape->positions[triangle.x]))
          return write_error();
        if (!format_values(fs, "vertex {}\n", shape->positions[triangle.y]))
          return write_error();
        if (!format_values(fs, "vertex {}\n", shape->positions[triangle.z]))
          return write_error();
        if (!format_values(fs, "endloop\n")) return write_error();
        if (!format_values(fs, "endfacet\n")) return write_error();
      }
      if (!format_values(fs, "endsolid \n")) return write_error();
    }
  }

  // done
  return true;
}

// Get/set data
bool get_triangles(const stl_model* stl, int shape_id, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& fnormals) {
  if (shape_id < 0 || shape_id >= stl->shapes.size()) return false;
  auto shape = stl->shapes.at(shape_id);
  triangles  = shape->triangles;
  positions  = shape->positions;
  fnormals   = shape->fnormals;
  return true;
}
void add_triangles(stl_model* stl, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& fnormals) {
  auto shape       = stl->shapes.emplace_back(new stl_shape{});
  shape->triangles = triangles;
  shape->positions = positions;
  shape->fnormals  = fnormals;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT PARSING
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt type
enum struct pbrt_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color,
  point2, vector2, spectrum
  // clang-format on
};

// Pbrt value
struct pbrt_value {
  string        name     = "";
  pbrt_type     type     = pbrt_type::real;
  int           value1i  = 0;
  float         value1f  = 0;
  vec2f         value2f  = {0, 0};
  vec3f         value3f  = {0, 0, 0};
  bool          value1b  = false;
  string        value1s  = "";
  vector<float> vector1f = {};
  vector<vec2f> vector2f = {};
  vector<vec3f> vector3f = {};
  vector<int>   vector1i = {};
};

// Pbrt command
struct pbrt_command {
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
};

// get pbrt value
inline bool get_pbrt_value(const pbrt_value& pbrt, string& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val = pbrt.value1s;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, bool& val) {
  if (pbrt.type == pbrt_type::boolean) {
    val = pbrt.value1b;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, int& val) {
  if (pbrt.type == pbrt_type::integer) {
    val = pbrt.value1i;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, float& val) {
  if (pbrt.type == pbrt_type::real) {
    val = pbrt.value1f;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vec2f& val) {
  if (pbrt.type == pbrt_type::point2 || pbrt.type == pbrt_type::vector2) {
    val = pbrt.value2f;
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vec3f& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    val = pbrt.value3f;
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    val = vec3f{pbrt.value1f, pbrt.value1f, pbrt.value1f};
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<float>& val) {
  if (pbrt.type == pbrt_type::real) {
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
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& val) {
  if (pbrt.type == pbrt_type::point2 || pbrt.type == pbrt_type::vector2) {
    if (!pbrt.vector2f.empty()) {
      val = pbrt.vector2f;
    } else {
      val = {pbrt.value2f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 2) != 0)
      throw std::runtime_error("bad pbrt type");
    val.resize(pbrt.vector1f.size() / 2);
    for (auto i = 0; i < val.size(); i++)
      val[i] = {pbrt.vector1f[i * 2 + 0], pbrt.vector1f[i * 2 + 1]};
    return true;
  } else {
    return false;
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& val) {
  if (pbrt.type == pbrt_type::point || pbrt.type == pbrt_type::vector ||
      pbrt.type == pbrt_type::normal || pbrt.type == pbrt_type::color) {
    if (!pbrt.vector3f.empty()) {
      val = pbrt.vector3f;
    } else {
      val = {pbrt.value3f};
    }
    return true;
  } else if (pbrt.type == pbrt_type::real) {
    if (pbrt.vector1f.empty() || (pbrt.vector1f.size() % 3) != 0)
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

inline bool get_pbrt_value(const pbrt_value& pbrt, vector<int>& val) {
  if (pbrt.type == pbrt_type::integer) {
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
inline bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& val) {
  if (pbrt.type == pbrt_type::integer) {
    if (pbrt.vector1i.empty() || (pbrt.vector1i.size() % 3) != 0)
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
inline bool get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = 0;
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
inline bool get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& val) {
  if (pbrt.type == pbrt_type::string || pbrt.type == pbrt_type::texture) {
    val.first = zero3f;
    return get_pbrt_value(pbrt, val.second);
  } else {
    val.second = "";
    return get_pbrt_value(pbrt, val.first);
  }
}
template <typename T>
inline bool get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& val) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_pbrt_value(p, val);
    }
  }
  return true;
}

// pbrt value construction
inline pbrt_value make_pbrt_value(
    const string& name, const string& val, pbrt_type type = pbrt_type::string) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, bool val, pbrt_type type = pbrt_type::boolean) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, int val, pbrt_type type = pbrt_type::integer) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, float val, pbrt_type type = pbrt_type::real) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec2f& val, pbrt_type type = pbrt_type::point2) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(
    const string& name, const vec3f& val, pbrt_type type = pbrt_type::color) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec2f>& val,
    pbrt_type type = pbrt_type::point2) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec3f>& val,
    pbrt_type type = pbrt_type::point) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = val;
  return pbrt;
}
inline pbrt_value make_pbrt_value(const string& name, const vector<vec3i>& val,
    pbrt_type type = pbrt_type::integer) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {&val.front().x, &val.front().x + val.size() * 3};
  return pbrt;
}

inline void remove_pbrt_comment(string_view& str, char comment_char = '#') {
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
inline bool read_pbrt_cmdline(file_stream& fs, string& cmd) {
  auto buffer = array<char, 4096>{};
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs.fs);
  while (read_line(fs, buffer)) {
    // line
    auto line = string_view{buffer.data()};
    remove_comment(line, '#', true);
    skip_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs.fs, pos, SEEK_SET);
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
    pos = ftell(fs.fs);
  }
  return found;
}

// parse a quoted string
inline bool parse_command(string_view& str, string& value) {
  skip_whitespace(str);
  if (!isalpha((int)str.front())) return false;
  auto pos = str.find_first_not_of(
      "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
  if (pos == string_view::npos) {
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
inline bool parse_param(string_view& str, T& value) {
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

// parse a quoted string
inline bool parse_nametype(string_view& str_, string& name, string& type) {
  auto value = ""s;
  if (!parse_value(str_, value)) return false;
  if (str_.empty()) return false;
  auto str  = string_view{value};
  auto pos1 = str.find(' ');
  if (pos1 == string_view::npos) return false;
  type = string(str.substr(0, pos1));
  str.remove_prefix(pos1);
  auto pos2 = str.find_first_not_of(' ');
  if (pos2 == string_view::npos) return false;
  str.remove_prefix(pos2);
  name = string(str);
  return true;
}

inline pair<vec3f, vec3f> get_etak(const string& name) {
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

// Pbrt measure subsurface parameters (sigma_prime_s, sigma_a in mm^-1)
// from pbrt code at pbrt/code/medium.cpp
inline pair<vec3f, vec3f> get_subsurface(const string& name) {
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

inline bool parse_params(string_view& str, vector<pbrt_value>& values) {
  auto parse_pvalues = [](string_view& str, auto& value, auto& values) -> bool {
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
      value.type = pbrt_type::real;
      if (!parse_pvalues(str, value.value1f, value.vector1f)) return false;
    } else if (type == "integer") {
      value.type = pbrt_type::integer;
      if (!parse_pvalues(str, value.value1i, value.vector1i)) return false;
    } else if (type == "string") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::string;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::boolean;
      if (!parse_pvalues(str, value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_type::texture;
      if (!parse_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "point" || type == "point3") {
      value.type = pbrt_type::point;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "normal" || type == "normal3") {
      value.type = pbrt_type::normal;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "vector" || type == "vector3") {
      value.type = pbrt_type::vector;
      parse_pvalues(str, value.value3f, value.vector3f);
    } else if (type == "point2") {
      value.type = pbrt_type::point2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "vector2") {
      value.type = pbrt_type::vector2;
      parse_pvalues(str, value.value2f, value.vector2f);
    } else if (type == "blackbody") {
      value.type     = pbrt_type::color;
      auto blackbody = zero2f;
      auto vector2f  = vector<vec2f>{};
      parse_pvalues(str, blackbody, vector2f);
      if (!vector2f.empty()) return false;
      value.value3f = blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "color" || type == "rgb") {
      value.type = pbrt_type::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "xyz") {
      value.type = pbrt_type::color;
      if (!parse_pvalues(str, value.value3f, value.vector3f)) return false;
      // xyz conversion
      return false;
    } else if (type == "spectrum") {
      auto is_string = false;
      auto str1      = str;
      skip_whitespace(str1);
      if (!str1.empty() && str1.front() == '"') {
        is_string = true;
      } else if (!str1.empty() && str1.front() == '[') {
        str1.remove_prefix(1);
        skip_whitespace(str1);
        if (!str1.empty() && str1.front() == '"') is_string = true;
      }
      if (is_string) {
        value.type     = pbrt_type::color;
        auto filename  = ""s;
        auto filenames = vector<string>{};
        if (!parse_value(str, filename)) return false;
        if (str.empty()) return false;
        auto filenamep = path_filename(filename);
        if (path_extension(filenamep) == ".spd") {
          filenamep = replace_extension(filenamep, "");
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (path_extension(filenamep) == ".eta") {
            auto eta      = get_etak(replace_extension(filenamep, "")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (path_extension(filenamep) == ".k") {
            auto k        = get_etak(replace_extension(filenamep, "")).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            return false;
          }
        } else {
          return false;
        }
      } else {
        value.type = pbrt_type::spectrum;
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
struct pbrt_film {
  // film approximation
  string filename   = "";
  vec2i  resolution = {0, 0};
};

// Pbrt texture
struct pbrt_texture {
  // texture parameters
  string name     = "";
  vec3f  constant = {1, 1, 1};
  string filename = "";
};

// Pbrt area light
struct pbrt_arealight {
  // arealight parameters
  string name     = "";
  vec3f  emission = {0, 0, 0};
};

// Pbrt medium. Not parsed at the moment.
struct pbrt_medium {
  // medium parameters
  string name = "";
};

// convert pbrt films
inline bool convert_film(pbrt_film* film, const pbrt_command& command,
    const string& filename, string& error, bool verbose = false) {
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
    if (!get_pbrt_value(command.values, "xresolution", film->resolution.x))
      return parse_error();
    if (!get_pbrt_value(command.values, "yresolution", film->resolution.y))
      return parse_error();
    film->filename = "out.png"s;
    if (!get_pbrt_value(command.values, "filename", film->filename))
      return parse_error();
    return true;
  } else {
    return type_error();
  }
}

// convert pbrt elements
inline bool convert_camera(pbrt_camera* pcamera, const pbrt_command& command,
    const vec2i& resolution, const string& filename, string& error,
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
    if (!get_pbrt_value(command.values, "fov", fov)) return parse_error();
    // auto lensradius = if(!get_pbrt_value(values, "lensradius", 0.0f);
    pcamera->aspect = film_aspect;
    if (pcamera->aspect >= 1) {
      pcamera->lens = (0.036 / pcamera->aspect) / (2 * tan(radians(fov) / 2));
    } else {
      pcamera->lens = (0.036 * pcamera->aspect) / (2 * tan(radians(fov) / 2));
    }
    if (!get_pbrt_value(command.values, "frameaspectratio", pcamera->aspect))
      return parse_error();
    pcamera->focus = 10.0f;
    if (!get_pbrt_value(command.values, "focaldistance", pcamera->focus))
      return parse_error();
    return true;
  } else if (command.type == "realistic") {
    auto lensfile = ""s;
    if (!get_pbrt_value(command.values, "lensfile", lensfile))
      return parse_error();
    lensfile          = lensfile.substr(0, lensfile.size() - 4);
    lensfile          = lensfile.substr(lensfile.find('.') + 1);
    lensfile          = lensfile.substr(0, lensfile.size() - 2);
    auto lens         = max((float)std::atof(lensfile.c_str()), 35.0f) * 0.001f;
    pcamera->lens     = 2 * atan(0.036f / (2 * lens));
    pcamera->aperture = 0.0f;
    if (!get_pbrt_value(command.values, "aperturediameter", pcamera->aperture))
      return parse_error();
    pcamera->focus = 10.0f;
    if (!get_pbrt_value(command.values, "focusdistance", pcamera->focus))
      return parse_error();
    pcamera->aspect = film_aspect;
    return true;
  } else {
    return type_error();
  }
}

// convert pbrt textures
inline bool convert_texture(pbrt_texture* ptexture, const pbrt_command& command,
    unordered_map<string, pbrt_texture>& texture_map, const string& filename,
    string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };

  auto make_filename = [&texture_map](const string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return pos->second.filename;
  };

  ptexture->name = command.name;
  if (command.type == "imagemap") {
    ptexture->filename = "";
    if (!get_pbrt_value(command.values, "filename", ptexture->filename))
      return parse_error();
    return true;
  } else if (command.type == "constant") {
    ptexture->constant = vec3f{1, 1, 1};
    if (!get_pbrt_value(command.values, "value", ptexture->constant))
      return parse_error();
    return true;
  } else if (command.type == "bilerp") {
    ptexture->constant = {1, 0, 0};
    return true;
  } else if (command.type == "checkerboard") {
    // auto tex1     = if(!get_pbrt_value(command.values, "tex1",
    // pair{vec3f{1},
    // ""s}); auto tex2     = if(!get_pbrt_value(command.values, "tex2",
    //  pair{vec3f{0}, ""s}); auto rgb1     = tex1.second == "" ?
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
    auto tex1 = pair{vec3f{0, 0, 0}, ""s}, tex2 = pair{vec3f{1, 1, 1}, ""s};
    if (!get_pbrt_value(command.values, "tex1", tex1)) return parse_error();
    if (!get_pbrt_value(command.values, "tex2", tex2)) return parse_error();
    if (!make_filename(tex1.second).empty()) {
      ptexture->filename = make_filename(tex1.second);
    } else if (!make_filename(tex2.second).empty()) {
      ptexture->filename = make_filename(tex2.second);
    } else {
      ptexture->constant = {1, 0, 0};
    }
    return true;
  } else if (command.type == "scale") {
    auto tex1 = pair{vec3f{1, 1, 1}, ""s}, tex2 = pair{vec3f{1, 1, 1}, ""s};
    if (!get_pbrt_value(command.values, "tex1", tex2)) return parse_error();
    if (!get_pbrt_value(command.values, "tex2", tex1)) return parse_error();
    if (!make_filename(tex1.second).empty()) {
      ptexture->filename = make_filename(tex1.second);
    } else if (!make_filename(tex2.second).empty()) {
      ptexture->filename = make_filename(tex2.second);
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
inline bool convert_material(pbrt_material*     pmaterial,
    const pbrt_command&                         command,
    const unordered_map<string, pbrt_material>& named_materials,
    const unordered_map<string, pbrt_texture>&  named_textures,
    const string& filename, string& error, bool verbose = false) {
  auto parse_error = [filename, &error]() {
    error = filename + ": parse error";
    return false;
  };
  auto type_error = [filename, &error, &command]() {
    error = filename + ": unknown type " + command.type;
    return false;
  };
  auto material_error = [filename, &error](const string& name) {
    error = filename + ": missing material " + name;
    return false;
  };
  auto bsdf_error = [filename, &error](const string& name) {
    error = filename + ": missing bsdf " + name;
    return false;
  };

  // helpers
  auto get_texture = [&](const vector<pbrt_value>& values, const string& name,
                         vec3f& color, string& filename,
                         const vec3f& def) -> bool {
    auto textured = pair{def, ""s};
    if (!get_pbrt_value(values, name, textured)) return parse_error();
    if (textured.second.empty()) {
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
  auto get_scalar = [&](const vector<pbrt_value>& values, const string& name,
                        float& scalar, float def) -> bool {
    auto textured = pair{vec3f{def, def, def}, ""s};
    if (!get_pbrt_value(values, name, textured)) return parse_error();
    if (textured.second.empty()) {
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
  auto get_color = [&](const vector<pbrt_value>& values, const string& name,
                       vec3f& color, const vec3f& def) -> bool {
    auto textured = pair{def, ""s};
    if (!get_pbrt_value(values, name, textured)) return parse_error();
    if (textured.second.empty()) {
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

  auto get_roughness = [&](const vector<pbrt_value>& values, float& roughness,
                           float def = 0.1) -> bool {
    auto roughness_ = pair{vec3f{def, def, def}, ""s};
    if (!get_pbrt_value(values, "roughness", roughness_)) return parse_error();
    auto uroughness = roughness_, vroughness = roughness_;
    auto remaproughness = true;
    if (!get_pbrt_value(values, "uroughness", uroughness)) return parse_error();
    if (!get_pbrt_value(values, "vroughness", vroughness)) return parse_error();
    if (!get_pbrt_value(values, "remaproughness", remaproughness))
      return parse_error();

    roughness = 0;
    if (uroughness.first == zero3f || vroughness.first == zero3f) return true;
    roughness = mean(vec2f{mean(uroughness.first), mean(vroughness.first)});
    // from pbrt code
    if (remaproughness) {
      roughness = max(roughness, 1e-3f);
      auto x    = log(roughness);
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
    if (!get_texture(command.values, "Kd", diffuse, diffuse_map,
            vec3f{0.25, 0.25, 0.25}))
      return parse_error();
    if (!get_texture(command.values, "Ks", specular, specular_map,
            vec3f{0.25, 0.25, 0.25}))
      return parse_error();
    if (!get_texture(command.values, "Kt", transmission, transmission_map,
            vec3f{0, 0, 0}))
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
    if (!get_roughness(command.values, pmaterial->roughness, 0.1))
      return parse_error();
    return true;
  } else if (command.type == "plastic") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.25, 0.25, 0.25}))
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
            pmaterial->color_tex, vec3f{0.25, 0.25, 0.25}))
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
            pmaterial->color_tex, vec3f{0.5, 0.5, 0.5}))
      return parse_error();
    return true;
  } else if (command.type == "mirror") {
    if (!get_texture(command.values, "Kr", pmaterial->color,
            pmaterial->color_tex, vec3f{0.9, 0.9, 0.9}))
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
            pmaterial->color_tex, vec3f{0.5, 0.5, 0.5}))
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
            pmaterial->color_tex, vec3f{0, 0, 0}))
      return parse_error();
    pmaterial->roughness = 1;
    if (verbose) printf("hair material not properly supported\n");
    return true;
  } else if (command.type == "disney") {
    if (!get_texture(command.values, "color", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5, 0.5, 0.5}))
      return parse_error();
    pmaterial->roughness = 1;
    if (verbose) printf("disney material not properly supported\n");
    return true;
  } else if (command.type == "kdsubsurface") {
    if (!get_texture(command.values, "Kd", pmaterial->color,
            pmaterial->color_tex, vec3f{0.5, 0.5, 0.5}))
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
    if (!get_pbrt_value(command.values, "scale", scale)) return parse_error();
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
    if (!get_pbrt_value(command.values, "namedmaterial1", namedmaterial1))
      return parse_error();
    if (!get_pbrt_value(command.values, "namedmaterial2", namedmaterial2))
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
    if (!get_pbrt_value(command.values, "bsdffile", bsdffile))
      return parse_error();
    if (bsdffile.rfind('/') != string::npos)
      bsdffile = bsdffile.substr(bsdffile.rfind('/') + 1);
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
inline void make_shape(vector<vec3i>& triangles, vector<vec3f>& positions,
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
inline void make_sphere(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_shape(
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
inline void make_disk(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
  make_shape(
      triangles, positions, normals, texcoords, steps,
      [radius](const vec2f& uv) {
        auto a = 2 * pif * uv.x;
        return radius * (1 - uv.y) * vec3f{cos(a), sin(a), 0};
      },
      [](const vec2f& uv) {
        return vec3f{0, 0, 1};
      });
}
inline void make_quad(vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, const vec2i& steps,
    float radius) {
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
inline bool convert_shape(pbrt_shape* shape, const pbrt_command& command,
    string& alphamap, const unordered_map<string, pbrt_texture>& named_textures,
    const string& ply_dirname, const string& filename, string& error,
    bool verbose = false) {
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
  auto get_alpha = [&](const vector<pbrt_value>& values, const string& name,
                       string& filename) -> bool {
    auto def      = 1.0f;
    auto textured = pair{def, ""s};
    if (!get_pbrt_value(values, name, textured)) return parse_error();
    if (textured.second.empty()) {
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
    if (!get_pbrt_value(command.values, "P", shape->positions))
      return parse_error();
    if (!get_pbrt_value(command.values, "N", shape->normals))
      return parse_error();
    if (!get_pbrt_value(command.values, "uv", shape->texcoords))
      return parse_error();
    for (auto& uv : shape->texcoords) uv.y = (1 - uv.y);
    if (!get_pbrt_value(command.values, "indices", shape->triangles))
      return parse_error();
    return true;
  } else if (command.type == "loopsubdiv") {
    shape->positions = {};
    shape->triangles = {};
    if (!get_pbrt_value(command.values, "P", shape->positions))
      return parse_error();
    if (!get_pbrt_value(command.values, "indices", shape->triangles))
      return parse_error();
    shape->normals.resize(shape->positions.size());
    // compute_normals(shape->normals, shape->triangles, shape->positions);
    return true;
  } else if (command.type == "plymesh") {
    shape->filename_ = ""s;
    if (!get_pbrt_value(command.values, "filename", shape->filename_))
      return parse_error();
    if (!get_alpha(command.values, "alpha", alphamap)) return parse_error();
    auto ply = std::make_unique<ply_model>();
    if (!load_ply(path_join(ply_dirname, shape->filename_), ply.get(), error))
      return dependent_error();
    get_positions(ply.get(), shape->positions);
    get_normals(ply.get(), shape->normals);
    get_texcoords(ply.get(), shape->texcoords);
    get_triangles(ply.get(), shape->triangles);
    return true;
  } else if (command.type == "sphere") {
    auto radius = 1.0f;
    if (!get_pbrt_value(command.values, "radius", radius)) return parse_error();
    make_sphere(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 16}, radius);
    return true;
  } else if (command.type == "disk") {
    auto radius = 1.0f;
    if (!get_pbrt_value(command.values, "radius", radius)) return parse_error();
    make_disk(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 1}, radius);
    return true;
  } else {
    return type_error();
  }
}

// Convert pbrt arealights
inline bool convert_arealight(pbrt_arealight* parealight,
    const pbrt_command& command, const string& filename, string& error,
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
    auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
    if (!get_pbrt_value(command.values, "L", l)) return parse_error();
    if (!get_pbrt_value(command.values, "scale", scale)) return parse_error();
    parealight->emission = l * scale;
    return true;
  } else {
    return type_error();
  }
}

// Convert pbrt lights
inline bool convert_light(pbrt_light* plight, const pbrt_command& command,
    const string& filename, string& error, bool verbose = false) {
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
    auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
    if (!get_pbrt_value(command.values, "L", l)) return parse_error();
    if (!get_pbrt_value(command.values, "scale", scale)) return parse_error();
    plight->emission = l * scale;
    plight->from     = vec3f{0, 0, 0};
    plight->to       = vec3f{0, 0, 1};
    if (!get_pbrt_value(command.values, "from", plight->from))
      return parse_error();
    if (!get_pbrt_value(command.values, "to", plight->to)) return parse_error();
    plight->distant       = true;
    auto distant_dist     = 100;
    auto size             = distant_dist * sin(5 * pif / 180);
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
    auto texcoords = vector<vec2f>{};
    make_quad(plight->area_triangles, plight->area_positions,
        plight->area_normals, texcoords, {4, 2}, size);
    return true;
  } else if (command.type == "point" || command.type == "goniometric" ||
             command.type == "spot") {
    auto i = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
    if (!get_pbrt_value(command.values, "I", i)) return parse_error();
    if (!get_pbrt_value(command.values, "scale", scale)) return parse_error();
    plight->emission = i * scale;
    plight->from     = zero3f;
    if (!get_pbrt_value(command.values, "from", plight->from))
      return parse_error();
    plight->area_emission = plight->emission;
    plight->area_frame    = plight->frame * translation_frame(plight->from);
    plight->area_frend    = plight->frend * translation_frame(plight->from);
    auto texcoords        = vector<vec2f>{};
    make_sphere(plight->area_triangles, plight->area_positions,
        plight->area_normals, texcoords, {4, 2}, 0.0025f);
    return true;
  } else {
    return type_error();
  }
}

inline bool convert_environment(pbrt_environment* penvironment,
    const pbrt_command& command, const string& filename, string& error,
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
    auto l = vec3f{1, 1, 1}, scale = vec3f{1, 1, 1};
    if (!get_pbrt_value(command.values, "L", l)) return parse_error();
    if (!get_pbrt_value(command.values, "scale", scale)) return parse_error();
    penvironment->emission     = scale * l;
    penvironment->emission_tex = ""s;
    if (!get_pbrt_value(command.values, "mapname", penvironment->emission_tex))
      return parse_error();
    return true;
  } else {
    return type_error();
  }
}

// pbrt stack ctm
struct pbrt_stack_element {
  frame3f        transform_start        = identity3x4f;
  frame3f        transform_end          = identity3x4f;
  pbrt_material  material               = {};
  pbrt_arealight arealight              = {};
  pbrt_medium    interior               = {};
  pbrt_medium    exterior               = {};
  bool           reverse                = false;
  bool           active_transform_start = true;
  bool           active_transform_end   = true;
};

// pbrt parsing context
struct pbrt_context {
  vector<pbrt_stack_element>                 stack           = {};
  unordered_map<string, pbrt_stack_element>  coordsys        = {};
  unordered_map<string, vector<pbrt_shape*>> objects         = {};
  string                                     cur_object      = "";
  vec2i                                      film_resolution = {512, 512};
};

// load pbrt
inline bool load_pbrt(const string& filename, pbrt_scene* pbrt, string& error,
    pbrt_context& ctx, unordered_map<string, string>& material_map,
    unordered_map<string, pbrt_material>& named_materials,
    unordered_map<string, pbrt_texture>&  named_textures,
    unordered_map<string, pbrt_medium>&   named_mediums,
    const string&                         ply_dirname) {
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
  auto command_error = [filename, &error](const string& cmd) {
    error = filename + ": unknown command " + cmd;
    return false;
  };
  auto stack_error = [filename, &error]() {
    error = filename + ": parse error (bad stack)";
    return false;
  };
  auto object_error = [filename, &error](const string& obj) {
    error = filename + ": unknown object " + obj;
    return false;
  };

  // open file
  auto fs = open_file(filename, "rt");
  if (!fs) return open_error();

  // helpers
  auto set_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= xform;
    if (ctx.active_transform_end) ctx.transform_end *= xform;
  };

  // init stack
  if (ctx.stack.empty()) ctx.stack.emplace_back();

  // parse command by command
  auto line = ""s;
  while (read_pbrt_cmdline(fs, line)) {
    auto str = string_view{line};
    // get command
    auto cmd = ""s;
    if (!parse_command(str, cmd)) return parse_error();
    if (cmd == "WorldBegin") {
      ctx.stack.push_back({});
    } else if (cmd == "WorldEnd") {
      if (ctx.stack.empty()) return stack_error();
      ctx.stack.pop_back();
      if (ctx.stack.size() != 1) return stack_error();
    } else if (cmd == "AttributeBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "AttributeEnd") {
      if (ctx.stack.empty()) return stack_error();
      ctx.stack.pop_back();
    } else if (cmd == "TransformBegin") {
      ctx.stack.push_back(ctx.stack.back());
    } else if (cmd == "TransformEnd") {
      if (ctx.stack.empty()) return stack_error();
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
        return object_error(object);
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
        return parse_error();
      }
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      if (!parse_param(str, xf)) return parse_error();
      set_transform(ctx.stack.back(), mat_to_frame(xf));
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      if (!parse_param(str, xf)) return parse_error();
      concat_transform(ctx.stack.back(), mat_to_frame(xf));
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
      concat_transform(
          ctx.stack.back(), rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
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
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Sampler") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "PixelFilter") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Film") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      auto film = pbrt_film{};
      if (!convert_film(&film, command, filename, error)) return false;
      ctx.film_resolution = film.resolution;
    } else if (cmd == "Accelerator") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
    } else if (cmd == "Camera") {
      auto command = pbrt_command{};
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto camera   = add_camera(pbrt);
      if (!convert_camera(
              camera, command, ctx.film_resolution, filename, error))
        return false;
    } else if (cmd == "Texture") {
      auto command  = pbrt_command{};
      auto comptype = ""s;
      auto str_     = string{str};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_param(str, comptype)) return parse_error();
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (!convert_texture(&named_textures[command.name], command,
              named_textures, filename, error))
        return false;
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        command     = pbrt_command{};
      command.name = "__unnamed__material__" + std::to_string(material_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      if (command.type.empty()) {
        ctx.stack.back().material = {};
      } else {
        ctx.stack.back().material = {};
        if (!convert_material(&ctx.stack.back().material, command,
                named_materials, named_textures, filename, error))
          return false;
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto command = pbrt_command{};
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
      auto command = pbrt_command{};
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
        material_map[matkey] = material->name;
      }
      shape->material = material_map.at(matkey);
      if (!ctx.cur_object.empty()) {
        ctx.objects[ctx.cur_object].push_back(shape);
      }
    } else if (cmd == "AreaLightSource") {
      static auto arealight_id = 0;
      auto        command      = pbrt_command{};
      command.name = "__unnamed__arealight__" + std::to_string(arealight_id++);
      if (!parse_param(str, command.type)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (!convert_arealight(
              &ctx.stack.back().arealight, command, filename, error))
        return false;
    } else if (cmd == "LightSource") {
      auto command = pbrt_command{};
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
      auto command = pbrt_command{};
      if (!parse_param(str, command.name)) return parse_error();
      if (!parse_params(str, command.values)) return parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (command.name == "type") command.type = value.value1s;
      auto medium                 = pbrt_medium{};
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
      if (!load_pbrt(path_join(path_dirname(filename), includename), pbrt,
              error, ctx, material_map, named_materials, named_textures,
              named_mediums, ply_dirname))
        return dependent_error();
    } else {
      return command_error(cmd);
    }
  }
  return true;
}

pbrt_scene::~pbrt_scene() {
  for (auto camera : cameras) delete camera;
  for (auto shape : shapes) delete shape;
  for (auto environment : environments) delete environment;
  for (auto light : lights) delete light;
  for (auto material : materials) delete material;
}

// Make pbrt
pbrt_camera* add_camera(pbrt_scene* pbrt) {
  return pbrt->cameras.emplace_back(new pbrt_camera{});
}
pbrt_shape* add_shape(pbrt_scene* pbrt) {
  return pbrt->shapes.emplace_back(new pbrt_shape{});
}
pbrt_material* add_material(pbrt_scene* pbrt) {
  return pbrt->materials.emplace_back(new pbrt_material{});
}
pbrt_environment* add_environment(pbrt_scene* pbrt) {
  return pbrt->environments.emplace_back(new pbrt_environment{});
}
pbrt_light* add_light(pbrt_scene* pbrt) {
  return pbrt->lights.emplace_back(new pbrt_light{});
}

// load pbrt
bool load_pbrt(const string& filename, pbrt_scene* pbrt, string& error) {
  auto ctx             = pbrt_context{};
  auto material_map    = unordered_map<string, string>{};
  auto named_materials = unordered_map<string, pbrt_material>{{"", {}}};
  auto named_mediums   = unordered_map<string, pbrt_medium>{{"", {}}};
  auto named_textures  = unordered_map<string, pbrt_texture>{{"", {}}};
  if (!load_pbrt(filename, pbrt, error, ctx, material_map, named_materials,
          named_textures, named_mediums, path_dirname(filename)))
    return false;

  // remove unused materials
  auto used_materials = unordered_set<string>{};
  for (auto shape : pbrt->shapes) used_materials.insert(shape->material);
  pbrt->materials.erase(
      std::remove_if(pbrt->materials.begin(), pbrt->materials.end(),
          [&used_materials](auto material) {
            auto found = used_materials.find(material->name) !=
                         used_materials.end();
            if (!found) delete material;
            return !found;
          }),
      pbrt->materials.end());

  return true;
}

inline void format_value(string& str, const pbrt_value& value) {
  static auto type_labels = unordered_map<pbrt_type, string>{
      {pbrt_type::real, "float"},
      {pbrt_type::integer, "integer"},
      {pbrt_type::boolean, "bool"},
      {pbrt_type::string, "string"},
      {pbrt_type::point, "point"},
      {pbrt_type::normal, "normal"},
      {pbrt_type::vector, "vector"},
      {pbrt_type::texture, "texture"},
      {pbrt_type::color, "rgb"},
      {pbrt_type::point2, "point2"},
      {pbrt_type::vector2, "vector2"},
      {pbrt_type::spectrum, "spectrum"},
  };

  auto format_vector = [](string& str, auto& values) {
    str += "[ ";
    for (auto& value : values) {
      str += " ";
      format_value(str, value);
    }
    str += " ]";
  };

  format_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case pbrt_type::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_value(str, value.value1f);
      }
      break;
    case pbrt_type::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_value(str, value.value1i);
      }
      break;
    case pbrt_type::boolean:
      format_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case pbrt_type::string:
    case pbrt_type::texture: format_values(str, "\"{}\"", value.value1s); break;
    case pbrt_type::point:
    case pbrt_type::vector:
    case pbrt_type::normal:
    case pbrt_type::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_values(str, "[ {} ]", value.value3f);
      }
      break;
    case pbrt_type::spectrum: format_vector(str, value.vector1f); break;
    case pbrt_type::point2:
    case pbrt_type::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

inline void format_value(string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

bool save_pbrt(const string& filename, const pbrt_scene* pbrt, string& error,
    bool ply_meshes) {
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
  auto fs = open_file(filename, "wt");
  if (!fs) return open_error();

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
    auto command = pbrt_command{};
    command.type = "image";
    command.values.push_back(
        make_pbrt_value("xresolution", camera->resolution.x));
    command.values.push_back(
        make_pbrt_value("yresolution", camera->resolution.y));
    command.values.push_back(make_pbrt_value("filename", "image.exr"s));
    if (!format_values(fs, "Film \"{}\" {}\n", command.type, command.values))
      return write_error();
  }

  for (auto camera : pbrt->cameras) {
    auto command  = pbrt_command{};
    command.type  = "perspective";
    command.frame = camera->frame;
    command.values.push_back(make_pbrt_value(
        "fov", 2 * tan(0.036f / (2 * camera->lens)) * 180 / pif));
    if (!format_values(fs, "LookAt {} {} {}\n", command.frame.o,
            command.frame.o - command.frame.z, command.frame.y))
      return write_error();
    if (!format_values(fs, "Camera \"{}\" {}\n", command.type, command.values))
      return write_error();
  }

  if (!format_values(fs, "\nWorldBegin\n\n")) return write_error();

  for (auto light : pbrt->lights) {
    auto command  = pbrt_command{};
    command.frame = light->frame;
    if (light->distant) {
      command.type = "distance";
      command.values.push_back(make_pbrt_value("L", light->emission));
    } else {
      command.type = "point";
      command.values.push_back(make_pbrt_value("I", light->emission));
    }
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", frame_to_mat(command.frame)))
      return write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
  }

  for (auto environment : pbrt->environments) {
    auto command  = pbrt_command{};
    command.frame = environment->frame;
    command.type  = "infinite";
    command.values.push_back(make_pbrt_value("L", environment->emission));
    command.values.push_back(
        make_pbrt_value("mapname", environment->emission_tex));
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", frame_to_mat(command.frame)))
      return write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
  }

  auto reflectivity_to_eta = [](const vec3f& reflectivity) {
    return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
  };

  auto material_map = unordered_map<string, pbrt_material*>{};
  for (auto material : pbrt->materials) {
    auto command = pbrt_command{};
    if (material->specular != 0 && material->transmission != 0 &&
        !material->thin) {
      command.type = "glass";
      command.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
      command.values.push_back(make_pbrt_value("Kt", vec3f{1, 1, 1}));
      command.values.push_back(
          make_pbrt_value("roughness", pow(material->roughness, 2)));
      command.values.push_back(make_pbrt_value("eta", material->ior));
      command.values.push_back(make_pbrt_value("remaproughness", false));
    } else if (material->metallic > 0.1f) {
      command.type = "metal";
      command.values.push_back(make_pbrt_value("Kr", vec3f{1, 1, 1}));
      command.values.push_back(
          make_pbrt_value("roughness", pow(material->roughness, 2)));
      command.values.push_back(
          make_pbrt_value("eta", reflectivity_to_eta(material->color)));
      command.values.push_back(make_pbrt_value("remaproughness", false));
    } else {
      command.type = "uber";
      if (material->color_tex.empty()) {
        command.values.push_back(make_pbrt_value("Kd", material->color));
      } else if (material->color != zero3f) {
        command.values.push_back(
            make_pbrt_value("Kd", material->color_tex, pbrt_type::texture));
      }
      if (material->specular != 0) {
        command.values.push_back(make_pbrt_value("Ks",
            vec3f{material->specular, material->specular, material->specular}));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material->roughness, 2)));
        command.values.push_back(make_pbrt_value("eta", material->ior));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      }
      if (material->transmission != 0) {
        command.values.push_back(make_pbrt_value(
            "Kt", vec3f{material->transmission, material->transmission,
                      material->transmission}));
      }
      if (!material->opacity_tex.empty()) {
        command.values.push_back(make_pbrt_value(
            "opacity", material->opacity_tex, pbrt_type::texture));
      } else if (material->opacity != 1) {
        command.values.push_back(make_pbrt_value("opacity", material->opacity));
      }
    }
    if (!format_values(fs,
            "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n",
            material->name, command.type, command.values))
      return write_error();
  }

  auto object_id = 0;
  for (auto shape : pbrt->shapes) {
    auto material = material_map.at(shape->material);
    auto command  = pbrt_command{};
    command.frame = shape->frame;
    if (ply_meshes) {
      command.type = "plymesh";
      command.values.push_back(make_pbrt_value("filename", shape->filename_));
    } else {
      command.type = "trianglemesh";
      command.values.push_back(make_pbrt_value("indices", shape->triangles));
      command.values.push_back(
          make_pbrt_value("P", shape->positions, pbrt_type::point));
      if (!shape->normals.empty())
        command.values.push_back(
            make_pbrt_value("N", shape->triangles, pbrt_type::normal));
      if (!shape->texcoords.empty())
        command.values.push_back(make_pbrt_value("uv", shape->texcoords));
    }
    if (ply_meshes) {
      auto ply_guard = std::make_unique<ply_model>();
      auto ply       = ply_guard.get();
      add_positions(ply, shape->positions);
      add_normals(ply, shape->normals);
      add_texcoords(ply, shape->texcoords);
      add_triangles(ply, shape->triangles);
      if (!save_ply(
              path_dirname(filename) + "/" + shape->filename_, ply, error))
        return dependent_error();
    }
    auto object = "object" + std::to_string(object_id++);
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectBegin \"{}\"\n", object))
        return write_error();
    if (!format_values(fs, "AttributeBegin\n")) return write_error();
    if (!format_values(fs, "Transform {}\n", frame_to_mat(shape->frame)))
      return write_error();
    if (material->emission != zero3f) {
      auto acommand = pbrt_command{};
      acommand.type = "diffuse";
      acommand.values.push_back(make_pbrt_value("L", material->emission));
      if (!format_values(fs, "AreaLightSource \"{}\" {}\n", acommand.type,
              acommand.values))
        return write_error();
    }
    if (!format_values(fs, "NamedMaterial \"{}\"\n", material->name))
      return write_error();
    if (!format_values(fs, "Shape \"{}\" {}\n", command.type, command.values))
      return write_error();
    if (!format_values(fs, "AttributeEnd\n")) return write_error();
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectEnd\n")) return write_error();
    for (auto& iframe : shape->instances) {
      if (!format_values(fs, "AttributeBegin\n")) return write_error();
      if (!format_values(fs, "Transform {}\n", frame_to_mat(iframe)))
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

}  // namespace yocto
