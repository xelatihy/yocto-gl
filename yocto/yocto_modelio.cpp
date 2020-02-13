//
// Implementation for Yocto/ModelIO
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_modelio.h"

#include <algorithm>
#include <memory>
#include <string_view>

#define CGLTF_IMPLEMENTATION
#include "ext/cgltf.h"
#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// -----------------------------------------------------------------------------
// LOAD-LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;

// utilities
static bool is_newline(char c) { return c == '\r' || c == '\n'; }
static bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}

// Parse values from a string
[[nodiscard]] static bool parse_value(string_view& str, string_view& value) {
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
[[nodiscard]] static bool parse_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  if (!parse_value(str, valuei)) return false;
  value = (bool)valuei;
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#ifdef __APPLE__
[[nodiscard]] static bool parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#endif

[[nodiscard]] static bool parse_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}
[[nodiscard]] static bool parse_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_value(str, value[i])) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL PRINTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
static void format_value(string& str, const string& value) { str += value; }
static void format_value(string& str, const char* value) { str += value; }
static void format_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
static void format_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
static void format_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static void format_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}

static void format_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
static void format_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
static void format_value(string& str, const vec4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
static void format_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}
static void format_value(string& str, const mat4f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_value(str, value[i]);
  }
}

// Foramt to file
static void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::invalid_argument("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
[[nodiscard]] static bool format_values(
    FILE* fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  if (fputs(str.c_str(), fs) < 0) return false;
  return true;
}
template <typename T>
[[nodiscard]] static bool format_value(FILE* fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  if (fputs(str.c_str(), fs) < 0) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL BINARY HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
static T swap_endian(T value) {
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
[[nodiscard]] bool read_value(FILE* fs, T& value, bool big_endian) {
  if (fread(&value, sizeof(value), 1, fs) != 1) return false;
  if (big_endian) value = swap_endian(value);
  return true;
}

template <typename T>
[[nodiscard]] bool write_value(FILE* fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  if (fwrite(&value, sizeof(value), 1, fs) != 1) return false;
  return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void remove_ply_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Read ply
shared_ptr<ply_model> make_ply() { return make_shared<ply_model>(); }
shared_ptr<ply_model> load_ply(const string& filename) {
  auto ply = make_ply();
  load_ply(filename, ply);
  return ply;
}

// Load ply
void load_ply(const string& filename, shared_ptr<ply_model> ply) {
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

  // open file
  auto fs = fopen(filename.c_str(), "rb");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_parse_error = [filename]() {
    throw std::runtime_error{filename + ": parse error"};
  };
  auto throw_read_error = [filename]() {
    throw std::runtime_error{filename + ": read error"};
  };

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
    if (!parse_value(str, cmd)) throw_parse_error();
    if (cmd == "") continue;

    // check magic number
    if (first_line) {
      if (cmd != "ply") throw_parse_error();
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) throw_parse_error();
    } else if (cmd == "format") {
      auto fmt = ""s;
      if (!parse_value(str, fmt)) throw_parse_error();
      if (fmt == "ascii") {
        ply->format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        ply->format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        ply->format = ply_format::binary_big_endian;
      } else {
        throw std::runtime_error{filename + ": parse error [bad header]"};
      }
    } else if (cmd == "comment") {
      skip_whitespace(str);
      ply->comments.push_back(string{str});
    } else if (cmd == "obj_info") {
      skip_whitespace(str);
      // comment is the rest of the str
    } else if (cmd == "element") {
      auto elem = ply->elements.emplace_back(make_shared<ply_element>());
      if (!parse_value(str, elem->name)) throw_parse_error();
      if (!parse_value(str, elem->count)) throw_parse_error();
    } else if (cmd == "property") {
      if (ply->elements.empty()) throw_parse_error();
      auto prop = ply->elements.back()->properties.emplace_back(
          make_shared<ply_property>());
      auto tname = ""s;
      if (!parse_value(str, tname)) throw_parse_error();
      if (tname == "list") {
        prop->is_list = true;
        if (!parse_value(str, tname)) throw_parse_error();
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8) throw_parse_error();
        if (!parse_value(str, tname)) throw_parse_error();
        if (type_map.find(tname) == type_map.end()) throw_parse_error();
        prop->type = type_map.at(tname);
      } else {
        prop->is_list = false;
        if (type_map.find(tname) == type_map.end()) throw_parse_error();
        prop->type = type_map.at(tname);
      }
      if (!parse_value(str, prop->name)) throw_parse_error();
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      throw_parse_error();
    }
  }

  // check exit
  if (!end_header) throw_parse_error();

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
        if (!fgets(buffer, sizeof(buffer), fs)) throw_read_error();
        auto str = string_view{buffer};
        for (auto prop : elem->properties) {
          if (prop->is_list) {
            if (!parse_value(str, prop->ldata_u8.emplace_back()))
              throw_parse_error();
          }
          auto vcount = prop->is_list ? prop->ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!parse_value(str, prop->data_i8.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::i16:
                if (!parse_value(str, prop->data_i16.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::i32:
                if (!parse_value(str, prop->data_i32.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::i64:
                if (!parse_value(str, prop->data_i64.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::u8:
                if (!parse_value(str, prop->data_u8.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::u16:
                if (!parse_value(str, prop->data_u16.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::u32:
                if (!parse_value(str, prop->data_u32.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::u64:
                if (!parse_value(str, prop->data_u64.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::f32:
                if (!parse_value(str, prop->data_f32.emplace_back()))
                  throw_parse_error();
                break;
              case ply_type::f64:
                if (!parse_value(str, prop->data_f64.emplace_back()))
                  throw_parse_error();
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
              throw_read_error();
          }
          auto vcount = prop->is_list ? prop->ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!read_value(fs, prop->data_i8.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::i16:
                if (!read_value(fs, prop->data_i16.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::i32:
                if (!read_value(fs, prop->data_i32.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::i64:
                if (!read_value(fs, prop->data_i64.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::u8:
                if (!read_value(fs, prop->data_u8.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::u16:
                if (!read_value(fs, prop->data_u16.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::u32:
                if (!read_value(fs, prop->data_u32.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::u64:
                if (!read_value(fs, prop->data_u64.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::f32:
                if (!read_value(fs, prop->data_f32.emplace_back(), big_endian))
                  throw_read_error();
                break;
              case ply_type::f64:
                if (!read_value(fs, prop->data_f64.emplace_back(), big_endian))
                  throw_read_error();
                break;
            }
          }
        }
      }
    }
  }
}

// Save ply
void save_ply(const string& filename, shared_ptr<ply_model> ply) {
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

  // open file
  auto fs = fopen(filename.c_str(), "wb");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_write_error = [filename]() {
    throw std::runtime_error{filename + ": write error"};
  };

  // header
  if (!format_values(fs, "ply\n")) throw_write_error();
  if (!format_values(fs, "format {} 1.0\n", format_map.at(ply->format)))
    throw_write_error();
  if (!format_values(fs, "comment Written by Yocto/GL\n")) throw_write_error();
  if (!format_values(fs, "comment https://github.com/xelatihy/yocto-gl\n"))
    throw_write_error();
  for (auto& comment : ply->comments)
    if (!format_values(fs, "comment {}\n", comment)) throw_write_error();
  for (auto elem : ply->elements) {
    if (!format_values(
            fs, "element {} {}\n", elem->name, (uint64_t)elem->count))
      throw_write_error();
    for (auto prop : elem->properties) {
      if (prop->is_list) {
        if (!format_values(fs, "property list uchar {} {}\n",
                type_map[prop->type], prop->name))
          throw_write_error();
      } else {
        if (!format_values(
                fs, "property {} {}\n", type_map[prop->type], prop->name))
          throw_write_error();
      }
    }
  }

  if (!format_values(fs, "end_header\n")) throw_write_error();

  // properties
  if (ply->format == ply_format::ascii) {
    for (auto elem : ply->elements) {
      auto cur = vector<size_t>(elem->properties.size(), 0);
      for (auto idx = 0; idx < elem->count; idx++) {
        for (auto pidx = 0; pidx < elem->properties.size(); pidx++) {
          auto prop = elem->properties[pidx];
          if (prop->is_list)
            if (!format_values(fs, "{} ", (int)prop->ldata_u8[idx]))
              throw_write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!format_values(fs, "{} ", prop->data_i8[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::i16:
                if (!format_values(fs, "{} ", prop->data_i16[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::i32:
                if (!format_values(fs, "{} ", prop->data_i32[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::i64:
                if (!format_values(fs, "{} ", prop->data_i64[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::u8:
                if (!format_values(fs, "{} ", prop->data_u8[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::u16:
                if (!format_values(fs, "{} ", prop->data_u16[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::u32:
                if (!format_values(fs, "{} ", prop->data_u32[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::u64:
                if (!format_values(fs, "{} ", prop->data_u64[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::f32:
                if (!format_values(fs, "{} ", prop->data_f32[cur[idx]++]))
                  throw_write_error();
                break;
              case ply_type::f64:
                if (!format_values(fs, "{} ", prop->data_f64[cur[idx]++]))
                  throw_write_error();
                break;
            }
          }
          if (!format_values(fs, "\n")) throw_write_error();
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
              throw_write_error();
          auto vcount = prop->is_list ? prop->ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop->type) {
              case ply_type::i8:
                if (!write_value(fs, prop->data_i8[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::i16:
                if (!write_value(fs, prop->data_i16[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::i32:
                if (!write_value(fs, prop->data_i32[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::i64:
                if (!write_value(fs, prop->data_i64[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::u8:
                if (!write_value(fs, prop->data_u8[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::u16:
                if (!write_value(fs, prop->data_u16[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::u32:
                if (!write_value(fs, prop->data_u32[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::u64:
                if (!write_value(fs, prop->data_u64[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::f32:
                if (!write_value(fs, prop->data_f32[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
              case ply_type::f64:
                if (!write_value(fs, prop->data_f64[cur[pidx]++], big_endian))
                  throw_write_error();
                break;
            }
          }
        }
      }
    }
  }
}

// Get ply properties
bool has_property(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return true;
    }
  }
  return false;
}
shared_ptr<ply_property> get_property(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
template <typename T, typename T1>
static vector<T> convert_ply_property(const vector<T1>& prop) {
  auto values = vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return values;
}
template <typename T>
static vector<T> convert_ply_property(shared_ptr<ply_property> prop) {
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
vector<float> get_values(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (prop->is_list) return {};
  return convert_ply_property<float>(prop);
}
vector<vec2f> get_values(shared_ptr<ply_model> ply, const string& element,
    const string& property1, const string& property2) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return values;
}
vector<vec3f> get_values(shared_ptr<ply_model> ply, const string& element,
    const string& property1, const string& property2, const string& property3) {
  auto x      = get_values(ply, element, property1);
  auto y      = get_values(ply, element, property2);
  auto z      = get_values(ply, element, property3);
  auto values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return values;
}
vector<vec4f> get_values(shared_ptr<ply_model> ply, const string& element,
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
vector<vec4f> get_values(shared_ptr<ply_model> ply, const string& element,
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
vector<frame3f> get_values(shared_ptr<ply_model> ply, const string& element,
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
vector<vector<int>> get_lists(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
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
vector<byte> get_list_sizes(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return prop->ldata_u8;
}
vector<int> get_list_values(
    shared_ptr<ply_model> ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return {};
  auto prop = get_property(ply, element, property);
  if (!prop->is_list) return {};
  return convert_ply_property<int>(prop);
}

static vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
vector<vec3f> get_positions(shared_ptr<ply_model> ply) {
  return get_values(ply, "vertex", "x", "y", "z");
}
vector<vec3f> get_normals(shared_ptr<ply_model> ply) {
  return get_values(ply, "vertex", "nx", "ny", "nz");
}
vector<vec2f> get_texcoords(shared_ptr<ply_model> ply, bool flipv) {
  auto texcoord = has_property(ply, "vertex", "u")
                      ? get_values(ply, "vertex", "u", "v")
                      : get_values(ply, "vertex", "s", "t");
  return flipv ? flip_ply_texcoord(texcoord) : texcoord;
}
vector<vec4f> get_colors(shared_ptr<ply_model> ply) {
  if (has_property(ply, "vertex", "alpha")) {
    return get_values(ply, "vertex", "red", "green", "blue", "alpha");
  } else {
    return get_values(ply, "vertex", "red", "green", "blue", 1);
  }
}
vector<float> get_radius(shared_ptr<ply_model> ply) {
  return get_values(ply, "vertex", "radius");
}
vector<vector<int>> get_faces(shared_ptr<ply_model> ply) {
  return get_lists(ply, "face", "vertex_indices");
}
vector<vec3i> get_triangles(shared_ptr<ply_model> ply) {
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
vector<vec4i> get_quads(shared_ptr<ply_model> ply) {
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
vector<vec2i> get_lines(shared_ptr<ply_model> ply) {
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
vector<int> get_points(shared_ptr<ply_model> ply) {
  return get_list_values(ply, "point", "vertex_indices");
}
bool has_quads(shared_ptr<ply_model> ply) {
  auto sizes = get_list_sizes(ply, "face", "vertex_indices");
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
static void add_element(
    shared_ptr<ply_model> ply, const string& element, size_t count) {
  for (auto elem : ply->elements) {
    if (elem->name == element) return;
  }
  auto elem   = ply->elements.emplace_back(make_shared<ply_element>());
  elem->name  = element;
  elem->count = count;
}
static void add_property(shared_ptr<ply_model> ply, const string& element,
    const string& property, size_t count, ply_type type, bool is_list) {
  add_element(ply, element, count);
  for (auto elem : ply->elements) {
    if (elem->name != element) continue;
    for (auto prop : elem->properties) {
      if (prop->name == property)
        throw std::runtime_error("property already added");
    }
    auto prop     = elem->properties.emplace_back(make_shared<ply_property>());
    prop->name    = property;
    prop->type    = type;
    prop->is_list = is_list;
    return;
  }
}
template <typename T>
static vector<T> make_ply_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

static void add_values(shared_ptr<ply_model> ply, const float* values,
    size_t count, const string& element, const string* properties, int nprops) {
  if (!values) return;
  for (auto p = 0; p < nprops; p++) {
    add_property(ply, element, properties[p], count, ply_type::f32, false);
    auto prop = get_property(ply, element, properties[p]);
    prop->data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_f32[i] = values[p + i * nprops];
  }
}

void add_values(shared_ptr<ply_model> ply, const vector<float>& values,
    const string& element, const string& property) {
  auto properties = vector{property};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
void add_values(shared_ptr<ply_model> ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2) {
  auto properties = vector{property1, property2};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
void add_values(shared_ptr<ply_model> ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3) {
  auto properties = vector{property1, property2, property3};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
void add_values(shared_ptr<ply_model> ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4) {
  auto properties = vector{property1, property2, property3, property4};
  add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}
void add_values(shared_ptr<ply_model> ply, const vector<frame3f>& values,
    const string& element, const array<string, 12>& properties) {
  add_values(ply, (float*)values.data(), values.size(), element,
      properties.data(), properties.size());
}

void add_lists(shared_ptr<ply_model> ply, const vector<vector<int>>& values,
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
void add_lists(shared_ptr<ply_model> ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property) {
  if (values.empty()) return;
  add_property(ply, element, property, sizes.size(), ply_type::i32, true);
  auto prop      = get_property(ply, element, property);
  prop->data_i32 = values;
  prop->ldata_u8 = sizes;
}
void add_lists(shared_ptr<ply_model> ply, const int* values, size_t count,
    int size, const string& element, const string& property) {
  if (!values) return;
  add_property(ply, element, property, count, ply_type::i32, true);
  auto prop = get_property(ply, element, property);
  prop->data_i32.assign(values, values + count * size);
  prop->ldata_u8.assign(count, size);
}
void add_lists(shared_ptr<ply_model> ply, const vector<int>& values,
    const string& element, const string& property) {
  return add_lists(ply, values.data(), values.size(), 1, element, property);
}
void add_lists(shared_ptr<ply_model> ply, const vector<vec2i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
void add_lists(shared_ptr<ply_model> ply, const vector<vec3i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
void add_lists(shared_ptr<ply_model> ply, const vector<vec4i>& values,
    const string& element, const string& property) {
  return add_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
}

// Add ply properties for meshes
void add_positions(shared_ptr<ply_model> ply, const vector<vec3f>& values) {
  return add_values(ply, values, "vertex", "x", "y", "z");
}
void add_normals(shared_ptr<ply_model> ply, const vector<vec3f>& values) {
  return add_values(ply, values, "vertex", "nx", "ny", "nz");
}
void add_texcoords(
    shared_ptr<ply_model> ply, const vector<vec2f>& values, bool flipv) {
  return add_values(
      ply, flipv ? flip_ply_texcoord(values) : values, "vertex", "u", "v");
}
void add_colors(shared_ptr<ply_model> ply, const vector<vec4f>& values) {
  return add_values(ply, values, "vertex", "red", "green", "blue", "alpha");
}
void add_radius(shared_ptr<ply_model> ply, const vector<float>& values) {
  return add_values(ply, values, "vertex", "radius");
}
void add_faces(shared_ptr<ply_model> ply, const vector<vector<int>>& values) {
  return add_lists(ply, values, "face", "vertex_indices");
}
void add_faces(shared_ptr<ply_model> ply, const vector<vec3i>& triangles,
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
void add_triangles(shared_ptr<ply_model> ply, const vector<vec3i>& values) {
  return add_faces(ply, values, {});
}
void add_quads(shared_ptr<ply_model> ply, const vector<vec4i>& values) {
  return add_faces(ply, {}, values);
}
void add_lines(shared_ptr<ply_model> ply, const vector<vec2i>& values) {
  return add_lists(ply, values, "str", "vertex_indices");
}
void add_points(shared_ptr<ply_model> ply, const vector<int>& values) {
  return add_lists(ply, values, "point", "vertex_indices");
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static void remove_obj_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

[[nodiscard]] static bool parse_value(string_view& str, obj_vertex& value) {
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
[[nodiscard]] static bool parse_value(
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
static void load_mtl(
    const string& filename, shared_ptr<obj_model> obj, bool fliptr = true) {
  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_parse_error = [filename]() {
    throw std::runtime_error{filename + ": parse error"};
  };
  auto throw_read_error = [filename]() {
    throw std::runtime_error{filename + ": read error"};
  };

  // init parsing
  auto material = obj->materials.emplace_back(make_shared<obj_material>());

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
    if (!parse_value(str, cmd)) throw_parse_error();
    if (cmd == "") continue;

    // grab material
    material = obj->materials.back();

    // possible token values
    if (cmd == "newmtl") {
      auto material = obj->materials.emplace_back(make_shared<obj_material>());
      if (!parse_value(str, material->name)) throw_parse_error();
    } else if (cmd == "illum") {
      if (!parse_value(str, material->illum)) throw_parse_error();
    } else if (cmd == "Ke") {
      if (!parse_value(str, material->emission)) throw_parse_error();
    } else if (cmd == "Ka") {
      if (!parse_value(str, material->ambient)) throw_parse_error();
    } else if (cmd == "Kd") {
      if (!parse_value(str, material->diffuse)) throw_parse_error();
    } else if (cmd == "Ks") {
      if (!parse_value(str, material->specular)) throw_parse_error();
    } else if (cmd == "Kt") {
      if (!parse_value(str, material->transmission)) throw_parse_error();
    } else if (cmd == "Tf") {
      material->transmission = vec3f{-1};
      if (!parse_value(str, material->transmission)) throw_parse_error();
      if (material->transmission.y < 0)
        material->transmission = vec3f{material->transmission.x};
      if (fliptr) material->transmission = 1 - material->transmission;
    } else if (cmd == "Tr") {
      if (!parse_value(str, material->opacity)) throw_parse_error();
      if (fliptr) material->opacity = 1 - material->opacity;
    } else if (cmd == "Ns") {
      if (!parse_value(str, material->exponent)) throw_parse_error();
    } else if (cmd == "d") {
      if (!parse_value(str, material->opacity)) throw_parse_error();
    } else if (cmd == "map_Ke") {
      if (!parse_value(str, material->emission_map)) throw_parse_error();
    } else if (cmd == "map_Ka") {
      if (!parse_value(str, material->ambient_map)) throw_parse_error();
    } else if (cmd == "map_Kd") {
      if (!parse_value(str, material->diffuse_map)) throw_parse_error();
    } else if (cmd == "map_Ks") {
      if (!parse_value(str, material->specular_map)) throw_parse_error();
    } else if (cmd == "map_Tr") {
      if (!parse_value(str, material->transmission_map)) throw_parse_error();
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      if (!parse_value(str, material->opacity_map)) throw_parse_error();
    } else if (cmd == "map_bump" || cmd == "bump") {
      if (!parse_value(str, material->bump_map)) throw_parse_error();
    } else if (cmd == "map_disp" || cmd == "disp") {
      if (!parse_value(str, material->displacement_map)) throw_parse_error();
    } else if (cmd == "map_norm" || cmd == "norm") {
      if (!parse_value(str, material->normal_map)) throw_parse_error();
    } else if (cmd == "Pe") {
      if (!parse_value(str, material->pbr_emission)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pb") {
      if (!parse_value(str, material->pbr_base)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Ps") {
      if (!parse_value(str, material->pbr_specular)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pm") {
      if (!parse_value(str, material->pbr_metallic)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pr") {
      if (!parse_value(str, material->pbr_roughness)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Ps") {
      if (!parse_value(str, material->pbr_sheen)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pc") {
      if (!parse_value(str, material->pbr_coat)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pcr") {
      if (!parse_value(str, material->pbr_coatroughness)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pt") {
      if (!parse_value(str, material->pbr_transmission)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pn") {
      if (!parse_value(str, material->pbr_ior)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Po") {
      if (!parse_value(str, material->pbr_opacity)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvs") {
      if (!parse_value(str, material->pbr_volscattering)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvg") {
      if (!parse_value(str, material->pbr_volanisotropy)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "Pvr") {
      if (!parse_value(str, material->pbr_volscale)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pe") {
      if (!parse_value(str, material->pbr_emission_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pb") {
      if (!parse_value(str, material->pbr_base_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Ps") {
      if (!parse_value(str, material->pbr_specular_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pm") {
      if (!parse_value(str, material->pbr_metallic_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pr") {
      if (!parse_value(str, material->pbr_roughness_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Ps") {
      if (!parse_value(str, material->pbr_sheen_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pc") {
      if (!parse_value(str, material->pbr_coat_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pcr") {
      if (!parse_value(str, material->pbr_coatroughness_map))
        throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Po") {
      if (!parse_value(str, material->pbr_opacity_map)) throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Pt") {
      if (!parse_value(str, material->pbr_transmission_map))
        throw_parse_error();
      material->as_pbr = true;
    } else if (cmd == "map_Vs") {
      if (!parse_value(str, material->pbr_volscattering_map))
        throw_parse_error();
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
    material->pbr_emission_map = material->emission_map;
    material->pbr_roughness    = exponent_to_roughness(material->exponent);
    material->pbr_ior          = material->ior;
    material->pbr_opacity      = material->opacity;
    material->pbr_opacity_map  = material->opacity_map;
    if (max(material->transmission) > 0.1) {
      material->pbr_base         = material->transmission;
      material->pbr_transmission = 1;
      material->pbr_specular     = 1;
    } else if (max(material->specular) > 0.2) {
      material->pbr_base     = material->specular;
      material->pbr_base_map = material->specular_map;
      material->pbr_metallic = 1;
    } else {
      material->pbr_base     = material->diffuse;
      material->pbr_base_map = material->diffuse_map;
      material->pbr_specular = max(material->specular) ? 1 : 0;
    }
  }
}

// Read obj
static void load_objx(const string& filename, shared_ptr<obj_model> obj) {
  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_parse_error = [filename]() {
    throw std::runtime_error{filename + ": parse error"};
  };
  auto throw_read_error = [filename]() {
    throw std::runtime_error{filename + ": read error"};
  };

  // shape map for instances
  auto shape_map = unordered_map<string, vector<shared_ptr<obj_shape>>>{};
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
    if (!parse_value(str, cmd)) throw_parse_error();
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto camera = obj->cameras.emplace_back(make_shared<obj_camera>());
      if (!parse_value(str, camera->name)) throw_parse_error();
      if (!parse_value(str, camera->ortho)) throw_parse_error();
      if (!parse_value(str, camera->width)) throw_parse_error();
      if (!parse_value(str, camera->height)) throw_parse_error();
      if (!parse_value(str, camera->lens)) throw_parse_error();
      if (!parse_value(str, camera->focus)) throw_parse_error();
      if (!parse_value(str, camera->aperture)) throw_parse_error();
      if (!parse_value(str, camera->frame)) throw_parse_error();
    } else if (cmd == "e") {
      auto environment = obj->environments.emplace_back(
          make_shared<obj_environment>());
      if (!parse_value(str, environment->name)) throw_parse_error();
      if (!parse_value(str, environment->emission)) throw_parse_error();
      auto emission_path = ""s;
      if (!parse_value(str, emission_path)) throw_parse_error();
      if (emission_path == "\"\"") emission_path = "";
      environment->emission_map.path = emission_path;
      if (!parse_value(str, environment->frame)) throw_parse_error();
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      if (!parse_value(str, object)) throw_parse_error();
      if (!parse_value(str, frame)) throw_parse_error();
      if (shape_map.find(object) == shape_map.end()) {
        throw std::runtime_error{filename + ": parse error [unknown object]"};
      }
      for (auto shape : shape_map.at(object)) {
        shape->instances.push_back(frame);
      }
    } else {
      // unused
    }
  }
}

// Read obj
shared_ptr<obj_model> make_obj() { return make_shared<obj_model>(); }
shared_ptr<obj_model> load_obj(const string& filename, bool geom_only,
    bool split_elements, bool split_materials) {
  auto obj = make_obj();
  load_obj(filename, obj, geom_only, split_elements, split_materials);
  return obj;
}

// Read obj
void load_obj(const string& filename, shared_ptr<obj_model> obj, bool geom_only,
    bool split_elements, bool split_materials) {
  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_parse_error = [filename]() {
    throw std::runtime_error{filename + ": parse error"};
  };
  auto throw_read_error = [filename]() {
    throw std::runtime_error{filename + ": read error"};
  };

  // parsing state
  auto opositions   = vector<vec3f>{};
  auto onormals     = vector<vec3f>{};
  auto otexcoords   = vector<vec2f>{};
  auto vert_size    = obj_vertex{};
  auto oname        = ""s;
  auto gname        = ""s;
  auto mname        = ""s;
  auto mtllibs      = vector<string>{};
  auto material_map = unordered_map<string, shared_ptr<obj_material>>{};

  // initialize obj
  obj->~obj_model();
  obj->cameras.clear();
  obj->environments.clear();
  obj->shapes.clear();
  obj->materials.clear();

  // initialize load
  obj->shapes.emplace_back(make_shared<obj_shape>());
  auto empty_material = shared_ptr<obj_material>{};

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
    if (!parse_value(str, cmd)) throw_parse_error();
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      if (!parse_value(str, opositions.emplace_back(zero3f)))
        throw_parse_error();
      vert_size.position += 1;
    } else if (cmd == "vn") {
      if (!parse_value(str, onormals.emplace_back(zero3f))) throw_parse_error();
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      if (!parse_value(str, otexcoords.emplace_back(zero2f)))
        throw_parse_error();
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      // split if split_elements and different primitives
      if (auto shape = obj->shapes.back();
          split_elements && !shape->vertices.empty()) {
        if ((cmd == "f" && (!shape->lines.empty() || !shape->points.empty())) ||
            (cmd == "l" && (!shape->faces.empty() || !shape->points.empty())) ||
            (cmd == "p" && (!shape->faces.empty() || !shape->lines.empty()))) {
          obj->shapes.emplace_back(make_shared<obj_shape>());
          obj->shapes.back()->name = oname + gname;
        }
      }
      // split if splt_material and different materials
      if (auto shape = obj->shapes.back();
          !geom_only && split_materials && !shape->materials.empty()) {
        if (shape->materials.size() > 1)
          throw std::runtime_error("should not have happened");
        if (shape->materials.back()->name != mname) {
          obj->shapes.emplace_back(make_shared<obj_shape>());
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
        if(mname.empty() && !empty_material) {
          empty_material = obj->materials.emplace_back(make_shared<obj_material>());
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
        if (!parse_value(str, vert)) throw_parse_error();
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
          if (!parse_value(str, oname)) throw_parse_error();
        }
      } else {
        if (str.empty()) {
          gname = "";
        } else {
          if (!parse_value(str, gname)) throw_parse_error();
        }
      }
      if (!obj->shapes.back()->vertices.empty()) {
        obj->shapes.emplace_back(make_shared<obj_shape>());
        obj->shapes.back()->name = oname + gname;
      } else {
        obj->shapes.back()->name = oname + gname;
      }
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      if (!parse_value(str, mname)) throw_parse_error();
    } else if (cmd == "s") {
      if (geom_only) continue;
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      if (!parse_value(str, mtllib)) throw_parse_error();
      if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) == mtllibs.end()) {
        mtllibs.push_back(mtllib);
        try {
          load_mtl(fs::path(filename).parent_path() / mtllib, obj);
        } catch (std::exception& e) {
          throw std::runtime_error{
              filename + ": error in resource (" + e.what() + ")"};
        }
        for (auto material : obj->materials)
          material_map[material->name] = material;
      }
    } else {
      // unused
    }
  }

  // fix empty material
  if(empty_material) {
    empty_material->name = "empty_material";
    empty_material->diffuse = {0.8, 0.8, 0.8};
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
  if (geom_only) return;

  // load extensions
  auto extfilename = fs::path(filename).replace_extension(".objx");
  if (fs::exists(fs::path(extfilename))) {
    try {
      load_objx(extfilename, obj);
    } catch (std::exception& e) {
      throw std::runtime_error{
          filename + ": error in resource (" + e.what() + ")"};
    }
  }
}

// Format values
static void format_value(string& str, const obj_texture_info& value) {
  str += value.path.empty() ? "" : value.path;
}
static void format_value(string& str, const obj_vertex& value) {
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
static void save_mtl(const string& filename, shared_ptr<obj_model> obj) {
  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_write_error = [filename]() {
    throw std::runtime_error{filename + ": write error"};
  };

  // save comments
  if (!format_values(fs, "#\n")) throw_write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) throw_write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    throw_write_error();
  if (!format_values(fs, "#\n\n")) throw_write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) throw_write_error();
  }
  if (!format_values(fs, "\n")) throw_write_error();

  // write material
  for (auto material : obj->materials) {
    if (!format_values(fs, "newmtl {}\n", material->name)) throw_write_error();
    if (!material->as_pbr) {
      if (!format_values(fs, "illum {}\n", material->illum))
        throw_write_error();
      if (material->emission != zero3f)
        if (!format_values(fs, "Ke {}\n", material->emission))
          throw_write_error();
      if (material->ambient != zero3f)
        if (!format_values(fs, "Ka {}\n", material->ambient))
          throw_write_error();
      if (!format_values(fs, "Kd {}\n", material->diffuse)) throw_write_error();
      if (!format_values(fs, "Ks {}\n", material->specular))
        throw_write_error();
      if (material->reflection != zero3f)
        if (!format_values(fs, "Kr {}\n", material->reflection))
          throw_write_error();
      if (material->transmission != zero3f)
        if (!format_values(fs, "Kt {}\n", material->transmission))
          throw_write_error();
      if (!format_values(fs, "Ns {}\n", (int)material->exponent))
        throw_write_error();
      if (material->opacity != 1)
        if (!format_values(fs, "d {}\n", material->opacity))
          throw_write_error();
      if (!material->emission_map.path.empty())
        if (!format_values(fs, "map_Ke {}\n", material->emission_map))
          throw_write_error();
      if (!material->diffuse_map.path.empty())
        if (!format_values(fs, "map_Kd {}\n", material->diffuse_map))
          throw_write_error();
      if (!material->specular_map.path.empty())
        if (!format_values(fs, "map_Ks {}\n", material->specular_map))
          throw_write_error();
      if (!material->transmission_map.path.empty())
        if (!format_values(fs, "map_Kt {}\n", material->transmission_map))
          throw_write_error();
      if (!material->reflection_map.path.empty())
        if (!format_values(fs, "map_Kr {}\n", material->reflection_map))
          throw_write_error();
      if (!material->exponent_map.path.empty())
        if (!format_values(fs, "map_Ns {}\n", material->exponent_map))
          throw_write_error();
      if (!material->opacity_map.path.empty())
        if (!format_values(fs, "map_d {}\n", material->opacity_map))
          throw_write_error();
      if (!material->bump_map.path.empty())
        if (!format_values(fs, "map_bump {}\n", material->bump_map))
          throw_write_error();
      if (!material->displacement_map.path.empty())
        if (!format_values(fs, "map_disp {}\n", material->displacement_map))
          throw_write_error();
      if (!material->normal_map.path.empty())
        if (!format_values(fs, "map_norm {}\n", material->normal_map))
          throw_write_error();
    } else {
      if (!format_values(fs, "illum 2\n")) throw_write_error();
      if (material->pbr_emission != zero3f)
        if (!format_values(fs, "Pe {}\n", material->pbr_emission))
          throw_write_error();
      if (material->pbr_base != zero3f)
        if (!format_values(fs, "Pb {}\n", material->pbr_base))
          throw_write_error();
      if (material->pbr_specular)
        if (!format_values(fs, "Psp {}\n", material->pbr_specular))
          throw_write_error();
      if (material->pbr_roughness)
        if (!format_values(fs, "Pr {}\n", material->pbr_roughness))
          throw_write_error();
      if (material->pbr_metallic)
        if (!format_values(fs, "Pm {}\n", material->pbr_metallic))
          throw_write_error();
      if (material->pbr_sheen)
        if (!format_values(fs, "Ps {}\n", material->pbr_sheen))
          throw_write_error();
      if (material->pbr_coat)
        if (!format_values(fs, "Pc {}\n", material->pbr_coat))
          throw_write_error();
      if (material->pbr_coatroughness)
        if (!format_values(fs, "Pcr {}\n", material->pbr_coatroughness))
          throw_write_error();
      if (material->pbr_volscattering != zero3f)
        if (!format_values(fs, "Pvs {}\n", material->pbr_volscattering))
          throw_write_error();
      if (material->pbr_volanisotropy)
        if (!format_values(fs, "Pvg {}\n", material->pbr_volanisotropy))
          throw_write_error();
      if (material->pbr_volscale)
        if (!format_values(fs, "Pvr {}\n", material->pbr_volscale))
          throw_write_error();
      if (!material->pbr_emission_map.path.empty())
        if (!format_values(fs, "map_Pe {}\n", material->pbr_emission_map))
          throw_write_error();
      if (!material->pbr_base_map.path.empty())
        if (!format_values(fs, "map_Pb {}\n", material->pbr_base_map))
          throw_write_error();
      if (!material->pbr_specular_map.path.empty())
        if (!format_values(fs, "map_Psp {}\n", material->pbr_specular_map))
          throw_write_error();
      if (!material->pbr_roughness_map.path.empty())
        if (!format_values(fs, "map_Pr {}\n", material->pbr_roughness_map))
          throw_write_error();
      if (!material->pbr_metallic_map.path.empty())
        if (!format_values(fs, "map_Pm {}\n", material->pbr_metallic_map))
          throw_write_error();
      if (!material->pbr_sheen_map.path.empty())
        if (!format_values(fs, "map_Ps {}\n", material->pbr_sheen_map))
          throw_write_error();
      if (!material->pbr_coat_map.path.empty())
        if (!format_values(fs, "map_Pc {}\n", material->pbr_coat_map))
          throw_write_error();
      if (!material->pbr_coatroughness_map.path.empty())
        if (!format_values(fs, "map_Pcr {}\n", material->pbr_coatroughness_map))
          throw_write_error();
      if (!material->pbr_volscattering_map.path.empty())
        if (!format_values(fs, "map_Pvs {}\n", material->pbr_volscattering_map))
          throw_write_error();
      if (!material->bump_map.path.empty())
        if (!format_values(fs, "map_bump {}\n", material->bump_map))
          throw_write_error();
      if (!material->displacement_map.path.empty())
        if (!format_values(fs, "map_disp {}\n", material->displacement_map))
          throw_write_error();
      if (!material->normal_map.path.empty())
        if (!format_values(fs, "map_norm {}\n", material->normal_map))
          throw_write_error();
    }
    if (!format_values(fs, "\n")) throw_write_error();
  }
}

// Save obj
static void save_objx(const string& filename, shared_ptr<obj_model> obj) {
  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_write_error = [filename]() {
    throw std::runtime_error{filename + ": write error"};
  };

  // save comments
  if (!format_values(fs, "#\n")) throw_write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) throw_write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    throw_write_error();
  if (!format_values(fs, "#\n\n")) throw_write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) throw_write_error();
  }
  if (!format_values(fs, "\n")) throw_write_error();

  // cameras
  for (auto camera : obj->cameras) {
    if (!format_values(fs, "c {} {} {} {} {} {} {} {}\n", camera->name,
            camera->ortho, camera->width, camera->height, camera->lens,
            camera->focus, camera->aperture, camera->frame))
      throw_write_error();
  }

  // environments
  for (auto environment : obj->environments) {
    if (!format_values(fs, "e {} {} {} {}\n", environment->name,
            environment->emission,
            environment->emission_map.path.empty()
                ? "\"\""s
                : environment->emission_map.path,
            environment->frame))
      throw_write_error();
  }

  // instances
  for (auto shape : obj->shapes) {
    for (auto& frame : shape->instances) {
      if (!format_values(fs, "i {} {}\n", shape->name, frame))
        throw_write_error();
    }
  }
}

// Save obj
void save_obj(const string& filename, shared_ptr<obj_model> obj) {
  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_write_error = [filename]() {
    throw std::runtime_error{filename + ": write error"};
  };

  // save comments
  if (!format_values(fs, "#\n")) throw_write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) throw_write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    throw_write_error();
  if (!format_values(fs, "#\n\n")) throw_write_error();
  for (auto& comment : obj->comments) {
    if (!format_values(fs, "# {}\n", comment)) throw_write_error();
  }
  if (!format_values(fs, "\n")) throw_write_error();

  // save material library
  if (!obj->materials.empty()) {
    if (!format_values(fs, "mtllib {}\n\n",
            fs::path(filename).filename().replace_extension(".mtl")))
      throw_write_error();
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto shape : obj->shapes) {
    if (!format_values(fs, "o {}\n", shape->name)) throw_write_error();
    for (auto& p : shape->positions)
      if (!format_values(fs, "v {}\n", p)) throw_write_error();
    for (auto& n : shape->normals)
      if (!format_values(fs, "vn {}\n", n)) throw_write_error();
    for (auto& t : shape->texcoords)
      if (!format_values(fs, "vt {}\n", t)) throw_write_error();
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
            throw_write_error();
          cur_material = element.material;
        }
        if (!format_values(fs, "{}", label)) throw_write_error();
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape->vertices[cur_vertex++];
          if (vert.position) vert.position += vert_size.position;
          if (vert.normal) vert.normal += vert_size.normal;
          if (vert.texcoord) vert.texcoord += vert_size.texcoord;
          if (!format_values(fs, " {}", vert)) throw_write_error();
        }
        if (!format_values(fs, "\n")) throw_write_error();
      }
    }
    if (!format_values(fs, "\n")) throw_write_error();
    vert_size.position += (int)shape->positions.size();
    vert_size.normal += (int)shape->normals.size();
    vert_size.texcoord += (int)shape->texcoords.size();
  }

  // save mtl
  if (!obj->materials.empty()) {
    save_mtl(fs::path(filename).replace_extension(".mtl"), obj);
  }

  // save objx
  if (!obj->cameras.empty() || !obj->environments.empty() ||
      std::any_of(obj->shapes.begin(), obj->shapes.end(),
          [](auto shape) { return !shape->instances.empty(); })) {
    save_objx(fs::path(filename).replace_extension(".objx"), obj);
  }
}

// Get obj vertices
static void get_vertices(shared_ptr<obj_shape> shape, vector<vec3f>& positions,
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
static vector<vec2f> flip_obj_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get obj shape
void get_triangles(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<shared_ptr<obj_material>>& materials,
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
void get_quads(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<shared_ptr<obj_material>>& materials,
    vector<int>& ematerials, bool flipv) {
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
void get_lines(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<shared_ptr<obj_material>>& materials,
    vector<int>& ematerials, bool flipv) {
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
void get_points(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<shared_ptr<obj_material>>& materials,
    vector<int>& ematerials, bool flipv) {
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
void get_fvquads(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<shared_ptr<obj_material>>& materials, vector<int>& ematerials,
    bool flipv) {
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

bool has_quads(shared_ptr<obj_shape> shape) {
  for (auto& face : shape->faces)
    if (face.size == 4) return true;
  return false;
}

// Get obj vertices
static void get_vertices(shared_ptr<obj_shape> shape, int material,
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
void get_triangles(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    int material, vector<vec3i>& triangles, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, bool flipv) {
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
void get_quads(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    int material, vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, bool flipv) {
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
void get_lines(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    int material, vector<vec2i>& lines, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, bool flipv) {
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
void get_points(shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape,
    int material, vector<int>& points, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, bool flipv) {
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
vector<shared_ptr<obj_material>> get_materials(
    shared_ptr<obj_model> obj, shared_ptr<obj_shape> shape) {
  return shape->materials;
}

// Add obj shape
void add_triangles(shared_ptr<obj_model> obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<shared_ptr<obj_material>>& materials,
    const vector<int>& ematerials, const vector<frame3f>& instances,
    bool flipv) {
  auto shape       = obj->shapes.emplace_back(make_shared<obj_shape>());
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
void add_quads(shared_ptr<obj_model> obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<shared_ptr<obj_material>>& materials,
    const vector<int>& ematerials, const vector<frame3f>& instances,
    bool flipv) {
  auto shape       = obj->shapes.emplace_back(make_shared<obj_shape>());
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
void add_lines(shared_ptr<obj_model> obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<shared_ptr<obj_material>>& materials,
    const vector<int>& ematerials, const vector<frame3f>& instances,
    bool flipv) {
  auto shape       = obj->shapes.emplace_back(make_shared<obj_shape>());
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
void add_points(shared_ptr<obj_model> obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<shared_ptr<obj_material>>& materials,
    const vector<int>& ematerials, const vector<frame3f>& instances,
    bool flipv) {
  auto shape       = obj->shapes.emplace_back(make_shared<obj_shape>());
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
void add_fvquads(shared_ptr<obj_model> obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<shared_ptr<obj_material>>& materials,
    const vector<int>& ematerials, const vector<frame3f>& instances,
    bool flipv) {
  auto shape       = obj->shapes.emplace_back(make_shared<obj_shape>());
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOAD-LEVEL PARSING
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

// Pbrt command
struct pbrt_command {
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
};

// get pbrt value
void get_pbrt_value(const pbrt_value& pbrt, string& value) {
  if (pbrt.type == pbrt_value_type::string ||
      pbrt.type == pbrt_value_type::texture) {
    value = pbrt.value1s;
  } else {
    throw std::invalid_argument{"expected string"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, bool& value) {
  if (pbrt.type == pbrt_value_type::boolean) {
    value = pbrt.value1b;
  } else {
    throw std::invalid_argument{"expected bool"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, int& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    value = pbrt.value1i;
  } else {
    throw std::invalid_argument{"expected int"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, float& value) {
  if (pbrt.type == pbrt_value_type::real) {
    value = pbrt.value1f;
  } else {
    throw std::invalid_argument{"expected float"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vec2f& value) {
  if (pbrt.type == pbrt_value_type::point2 ||
      pbrt.type == pbrt_value_type::vector2) {
    value = pbrt.value2f;
  } else {
    throw std::invalid_argument{"expected float2"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vec3f& value) {
  if (pbrt.type == pbrt_value_type::point ||
      pbrt.type == pbrt_value_type::vector ||
      pbrt.type == pbrt_value_type::normal ||
      pbrt.type == pbrt_value_type::color) {
    value = pbrt.value3f;
  } else if (pbrt.type == pbrt_value_type::real) {
    value = vec3f{pbrt.value1f};
  } else {
    throw std::invalid_argument{"expected float3"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vector<float>& value) {
  if (pbrt.type == pbrt_value_type::real) {
    if (!pbrt.vector1f.empty()) {
      value = pbrt.vector1f;
    } else {
      value = {pbrt.value1f};
    }
  } else {
    throw std::invalid_argument{"expected float array"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& value) {
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
    throw std::invalid_argument{"expected float2 array"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& value) {
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
      throw std::invalid_argument{"expected float3 array"};
    value.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
  } else {
    throw std::invalid_argument{"expected float3 array"};
  }
}

void get_pbrt_value(const pbrt_value& pbrt, vector<int>& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    if (!pbrt.vector1i.empty()) {
      value = pbrt.vector1i;
    } else {
      value = {pbrt.vector1i};
    }
  } else {
    throw std::invalid_argument{"expected int array"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    if (pbrt.vector1i.empty() || pbrt.vector1i.size() % 3)
      throw std::invalid_argument{"expected int3 array"};
    value.resize(pbrt.vector1i.size() / 3);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1i[i * 3 + 0], pbrt.vector1i[i * 3 + 1],
          pbrt.vector1i[i * 3 + 2]};
  } else {
    throw std::invalid_argument{"expected int3 array"};
  }
}
void get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value) {
  if (pbrt.type == pbrt_value_type::string) {
    value.first = 0;
    return get_pbrt_value(pbrt, value.second);
  } else {
    value.second = "";
    return get_pbrt_value(pbrt, value.first);
  }
}
void get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value) {
  if (pbrt.type == pbrt_value_type::string ||
      pbrt.type == pbrt_value_type::texture) {
    value.first = zero3f;
    return get_pbrt_value(pbrt, value.second);
  } else {
    value.second = "";
    return get_pbrt_value(pbrt, value.first);
  }
}
template <typename T>
inline void get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& value) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_pbrt_value(p, value);
    }
  }
}

// pbrt value construction
pbrt_value make_pbrt_value(const string& name, const string& value,
    pbrt_value_type type = pbrt_value_type::string) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, bool value,
    pbrt_value_type type = pbrt_value_type::boolean) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, int value,
    pbrt_value_type type = pbrt_value_type::integer) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, float value,
    pbrt_value_type type = pbrt_value_type::real) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, const vec2f& value,
    pbrt_value_type type = pbrt_value_type::point2) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, const vec3f& value,
    pbrt_value_type type = pbrt_value_type::color) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, const vector<vec2f>& value,
    pbrt_value_type type = pbrt_value_type::point2) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, const vector<vec3f>& value,
    pbrt_value_type type = pbrt_value_type::point) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(const string& name, const vector<vec3i>& value,
    pbrt_value_type type = pbrt_value_type::integer) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {(int*)value.data(), (int*)value.data() + value.size() * 3};
  return pbrt;
}

using std::string_view;

static void remove_pbrt_comment(string_view& str, char comment_char = '#') {
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
[[nodiscard]] static bool read_pbrt_cmdline(FILE* fs, string& cmd) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs);
  while (fgets(buffer, sizeof(buffer), fs)) {
    // line
    auto line = string_view{buffer};
    remove_pbrt_comment(line);
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

// parse a quoted string
[[nodiscard]] static bool parse_pbrt_command(string_view& str, string& value) {
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
[[nodiscard]] static bool parse_pbrt_param(string_view& str, T& value) {
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
[[nodiscard]] static bool parse_pbrt_nametype(
    string_view& str_, string& name, string& type) {
  auto value = ""s;
  if (!parse_value(str_, value)) return false;
  if (!str_.data()) return false;
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

static pair<vec3f, vec3f> get_pbrt_etak(const string& name) {
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
// static
pair<vec3f, vec3f> get_pbrt_subsurface(const string& name) {
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

[[nodiscard]] static bool parse_pbrt_params(
    string_view& str, vector<pbrt_value>& values) {
  auto parse_pbrt_pvalues = [](string_view& str, auto& value,
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
    if (!parse_pbrt_nametype(str, value.name, type)) return false;
    skip_whitespace(str);
    if (str.empty()) return false;
    if (type == "float") {
      value.type = pbrt_value_type::real;
      if (!parse_pbrt_pvalues(str, value.value1f, value.vector1f)) return false;
    } else if (type == "integer") {
      value.type = pbrt_value_type::integer;
      if (!parse_pbrt_pvalues(str, value.value1i, value.vector1i)) return false;
    } else if (type == "string") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::string;
      if (!parse_pbrt_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
    } else if (type == "bool") {
      auto value1s  = ""s;
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::boolean;
      if (!parse_pbrt_pvalues(str, value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
      value.value1b = value1s == "true";
    } else if (type == "texture") {
      auto vector1s = vector<string>{};
      value.type    = pbrt_value_type::texture;
      if (!parse_pbrt_pvalues(str, value.value1s, vector1s)) return false;
      if (!vector1s.empty()) return false;
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
      if (!vector2f.empty()) return false;
      value.value3f = blackbody_to_rgb(blackbody.x) * blackbody.y;
    } else if (type == "color" || type == "rgb") {
      value.type = pbrt_value_type::color;
      if (!parse_pbrt_pvalues(str, value.value3f, value.vector3f)) return false;
    } else if (type == "xyz") {
      value.type = pbrt_value_type::color;
      if (!parse_pbrt_pvalues(str, value.value3f, value.vector3f)) return false;
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
        value.type     = pbrt_value_type::color;
        auto filename  = ""s;
        auto filenames = vector<string>{};
        if (!parse_value(str, filename)) return false;
        if (!str.data()) return false;
        auto filenamep = fs::path(filename).filename();
        if (fs::path(filenamep).extension() == ".spd") {
          filenamep = fs::path(filenamep).replace_extension("").string();
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (fs::path(filenamep).extension() == ".eta") {
            auto eta =
                get_pbrt_etak(fs::path(filenamep).replace_extension("")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (fs::path(filenamep).extension() == ".k") {
            auto k =
                get_pbrt_etak(fs::path(filenamep).replace_extension("")).second;
            value.value3f = {k.x, k.y, k.z};
          } else {
            return false;
          }
        } else {
          return false;
        }
      } else {
        value.type = pbrt_value_type::spectrum;
        if (!parse_pbrt_pvalues(str, value.value1f, value.vector1f))
          return false;
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
  vec2i  resolution = zero2i;
};

// Pbrt texture
struct pbrt_texture {
  // texture parameters
  string name     = "";
  vec3f  constant = vec3f{1, 1, 1};
  string filename = "";
};

// convert pbrt films
static shared_ptr<pbrt_film> convert_film(
    const pbrt_command& command, bool verbose = false) {
  auto film = make_shared<pbrt_film>();
  if (command.type == "image") {
    film->resolution = {512, 512};
    get_pbrt_value(command.values, "xresolution", film->resolution.x);
    get_pbrt_value(command.values, "yresolution", film->resolution.y);
    film->filename = "out.png"s;
    get_pbrt_value(command.values, "filename", film->filename);
  } else {
    throw std::invalid_argument{"unknown film " + command.type};
  }
  return film;
}

// convert pbrt elements
static shared_ptr<pbrt_camera> convert_camera(const pbrt_command& command,
    const vec2i& resolution, bool verbose = false) {
  auto camera        = make_shared<pbrt_camera>();
  camera->frame      = command.frame;
  camera->frend      = command.frend;
  camera->frame      = inverse((frame3f)camera->frame);
  camera->frame.z    = -camera->frame.z;
  camera->resolution = resolution;
  auto film_aspect =
      (resolution == zero2i) ? 1 : (float)resolution.x / (float)resolution.y;
  if (command.type == "perspective") {
    auto fov = 90.0f;
    get_pbrt_value(command.values, "fov", fov);
    // auto lensradius = get_pbrt_value(values, "lensradius", 0.0f);
    camera->aspect = film_aspect;
    if (camera->aspect >= 1) {
      camera->lens = (0.036 / camera->aspect) / (2 * tan(radians(fov) / 2));
    } else {
      camera->lens = (0.036 * camera->aspect) / (2 * tan(radians(fov) / 2));
    }
    get_pbrt_value(command.values, "frameaspectratio", camera->aspect);
    camera->focus = 10.0f;
    get_pbrt_value(command.values, "focaldistance", camera->focus);
  } else if (command.type == "realistic") {
    auto lensfile = ""s;
    get_pbrt_value(command.values, "lensfile", lensfile);
    lensfile         = lensfile.substr(0, lensfile.size() - 4);
    lensfile         = lensfile.substr(lensfile.find('.') + 1);
    lensfile         = lensfile.substr(0, lensfile.size() - 2);
    auto lens        = max(std::atof(lensfile.c_str()), 35.0f) * 0.001f;
    camera->lens     = 2 * atan(0.036f / (2 * lens));
    camera->aperture = 0.0f;
    get_pbrt_value(command.values, "aperturediameter", camera->aperture);
    camera->focus = 10.0f;
    get_pbrt_value(command.values, "focusdistance", camera->focus);
    camera->aspect = film_aspect;
  } else {
    throw std::invalid_argument{"unknown camera " + command.type};
  }
  return camera;
}

// convert pbrt textures
static void convert_texture(pbrt_texture& texture, const pbrt_command& command,
    unordered_map<string, pbrt_texture>& texture_map, bool verbose = false) {
  auto get_filename = [&texture_map](const string& name) {
    if (name.empty()) return ""s;
    auto pos = texture_map.find(name);
    if (pos == texture_map.end()) return ""s;
    return pos->second.filename;
  };

  texture.name = command.name;
  if (command.type == "imagemap") {
    texture.filename = "";
    get_pbrt_value(command.values, "filename", texture.filename);
  } else if (command.type == "constant") {
    texture.constant = vec3f{1};
    get_pbrt_value(command.values, "value", texture.constant);
  } else if (command.type == "bilerp") {
    texture.constant = {1, 0, 0};
  } else if (command.type == "checkerboard") {
    // auto tex1     = get_pbrt_value(command.values, "tex1", pair{vec3f{1},
    // ""s}); auto tex2     = get_pbrt_value(command.values, "tex2",
    // pair{vec3f{0}, ""s}); auto rgb1     = tex1.second == "" ? tex1.first :
    // vec3f{0.4f, 0.4f, 0.4f}; auto rgb2     = tex1.second == "" ? tex2.first :
    // vec3f{0.6f, 0.6f, 0.6f}; auto params   = proc_image_params{}; params.type
    // = proc_image_params::type_t::checker; params.color0 = {rgb1.x, rgb1.y,
    // rgb1.z, 1}; params.color1 = {rgb2.x, rgb2.y, rgb2.z, 1}; params.scale
    // = 2; make_proc_image(texture.hdr, params); float_to_byte(texture.ldr,
    // texture.hdr); texture.hdr = {};
    texture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "dots") {
    texture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "fbm") {
    texture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "marble") {
    texture.constant = {0.5, 0.5, 0.5};
  } else if (command.type == "mix") {
    auto tex1 = pair{vec3f{0}, ""s}, tex2 = pair{vec3f{1}, ""s};
    get_pbrt_value(command.values, "tex1", tex1);
    get_pbrt_value(command.values, "tex2", tex2);
    if (!get_filename(tex1.second).empty()) {
      texture.filename = get_filename(tex1.second);
    } else if (!get_filename(tex2.second).empty()) {
      texture.filename = get_filename(tex2.second);
    } else {
      texture.constant = {1, 0, 0};
    }
  } else if (command.type == "scale") {
    auto tex1 = pair{vec3f{1}, ""s}, tex2 = pair{vec3f{1}, ""s};
    get_pbrt_value(command.values, "tex1", tex2);
    get_pbrt_value(command.values, "tex2", tex1);
    if (!get_filename(tex1.second).empty()) {
      texture.filename = get_filename(tex1.second);
    } else if (!get_filename(tex2.second).empty()) {
      texture.filename = get_filename(tex2.second);
    } else {
      texture.constant = {1, 0, 0};
    }
  } else if (command.type == "uv") {
    texture.constant = {1, 0, 0};
  } else if (command.type == "windy") {
    texture.constant = {1, 0, 0};
  } else if (command.type == "wrinkled") {
    texture.constant = {1, 0, 0};
  } else {
    throw std::invalid_argument{"unknown texture " + command.type};
  }
}

// convert pbrt materials
static shared_ptr<pbrt_material> convert_material(const pbrt_command& command,
    const unordered_map<string, shared_ptr<pbrt_material>>& material_map,
    const unordered_map<string, pbrt_texture>&              texture_map,
    bool                                                    verbose = false) {
  // helpers
  auto get_texture = [&](const vector<pbrt_value>& values, const string& name,
                         vec3f& color, string& filename,
                         const vec3f& def) -> void {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second == "") {
      color    = textured.first;
      filename = "";
    } else {
      auto& texture = texture_map.at(textured.second);
      if (texture.filename.empty()) {
        color    = texture.constant;
        filename = "";
      } else {
        color    = {1, 1, 1};
        filename = texture.filename;
      }
    }
  };
  auto get_scalar = [&](const vector<pbrt_value>& values, const string& name,
                        float& scalar, float def) -> void {
    auto textured = pair{vec3f{def}, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second == "") {
      scalar = mean(textured.first);
    } else {
      auto& texture = texture_map.at(textured.second);
      if (texture.filename.empty()) {
        scalar = mean(texture.constant);
      } else {
        scalar = def;
      }
    }
  };
  auto get_color = [&](const vector<pbrt_value>& values, const string& name,
                       vec3f& color, const vec3f& def) -> void {
    auto textured = pair{def, ""s};
    get_pbrt_value(values, name, textured);
    if (textured.second == "") {
      color = textured.first;
    } else {
      auto& texture = texture_map.at(textured.second);
      if (texture.filename.empty()) {
        color = texture.constant;
      } else {
        color = def;
      }
    }
  };

  auto get_roughness = [&](const vector<pbrt_value>& values, float& roughness,
                           float def = 0.1) -> void {
    auto roughness_ = pair{vec3f{def}, ""s};
    get_pbrt_value(values, "roughness", roughness_);
    auto uroughness = roughness_, vroughness = roughness_;
    auto remaproughness = true;
    get_pbrt_value(values, "uroughness", uroughness);
    get_pbrt_value(values, "vroughness", vroughness);
    get_pbrt_value(values, "remaproughness", remaproughness);

    roughness = 0;
    if (uroughness.first == zero3f || vroughness.first == zero3f) return;
    roughness = mean(vec2f{mean(uroughness.first), mean(vroughness.first)});
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

  auto material  = make_shared<pbrt_material>();
  material->name = command.name;
  if (command.type == "uber") {
    auto diffuse = zero3f, specular = zero3f, transmission = zero3f;
    auto diffuse_map = ""s, specular_map = ""s, transmission_map = ""s;
    get_texture(command.values, "Kd", diffuse, diffuse_map, vec3f{0.25});
    get_texture(command.values, "Ks", specular, specular_map, vec3f{0.25});
    get_texture(command.values, "Kt", transmission, transmission_map, vec3f{0});
    if (max(transmission) > 0.1) {
      material->color        = transmission;
      material->color_map    = transmission_map;
      material->specular     = 1;
      material->transmission = 1;
    } else {
      material->color     = diffuse;
      material->color_map = diffuse_map;
      material->specular  = 1;
    }
    get_scalar(command.values, "opacity", material->opacity, 1);
    get_scalar(command.values, "eta", material->ior, 1.5);
    get_roughness(command.values, material->roughness, 0.1f);
  } else if (command.type == "plastic") {
    get_texture(command.values, "Kd", material->color, material->color_map,
        vec3f{0.25});
    get_scalar(command.values, "Ks", material->specular, 0.25);
    get_scalar(command.values, "eta", material->ior, 1.5);
    material->roughness = 0.1f;
    get_roughness(command.values, material->roughness, 0.1);
  } else if (command.type == "translucent") {
    get_texture(command.values, "Kd", material->color, material->color_map,
        vec3f{0.25});
    get_scalar(command.values, "Ks", material->specular, 0.25);
    get_scalar(command.values, "eta", material->ior, 1.5);
    get_roughness(command.values, material->roughness, 0.1);
  } else if (command.type == "matte") {
    get_texture(
        command.values, "Kd", material->color, material->color_map, vec3f{0.5});
  } else if (command.type == "mirror") {
    get_texture(
        command.values, "Kr", material->color, material->color_map, vec3f{0.9});
    material->metallic  = 1;
    material->roughness = 0;
  } else if (command.type == "metal") {
    // get_texture(
    //     values, "Kr", material->specular, material->specular_map,
    //     vec3f{1});
    auto eta = zero3f, etak = zero3f;
    get_color(command.values, "eta", eta,
        vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f});
    get_color(command.values, "k", etak,
        vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f});
    material->color     = eta_to_reflectivity(eta, etak);
    material->roughness = 0.01f;
    get_roughness(command.values, material->roughness, 0.01);
  } else if (command.type == "substrate") {
    get_texture(
        command.values, "Kd", material->color, material->color_map, vec3f{0.5});
    get_scalar(command.values, "Ks", material->specular, 0.5);
    get_scalar(command.values, "eta", material->ior, 1.5);
    material->roughness = 0.1f;
    get_roughness(command.values, material->roughness, 0.1);
  } else if (command.type == "glass") {
    // get_texture(
    //     values, "Kr", material->specular, material->specular_map,
    //     vec3f{1});
    // get_texture(command.values, "Kt", material->transmission,
    //     material->transmission_map, vec3f{1});
    material->color        = {1, 1, 1};
    material->specular     = 1;
    material->transmission = 1;
    material->thin         = false;
    get_scalar(command.values, "eta", material->ior, 1.5);
    material->roughness = 0;
    get_roughness(command.values, material->roughness, 0);
  } else if (command.type == "hair") {
    get_texture(command.values, "color", material->color, material->color_map,
        vec3f{0});
    material->roughness = 1;
    if (verbose) printf("hair material not properly supported\n");
  } else if (command.type == "disney") {
    get_texture(command.values, "color", material->color, material->color_map,
        vec3f{0.5});
    material->roughness = 1;
    if (verbose) printf("disney material not properly supported\n");
  } else if (command.type == "kdsubsurface") {
    get_texture(
        command.values, "Kd", material->color, material->color_map, vec3f{0.5});
    get_scalar(command.values, "Kr", material->specular, 1);
    get_scalar(command.values, "eta", material->ior, 1.5);
    material->roughness = 0;
    get_roughness(command.values, material->roughness, 0);
    if (verbose) printf("kdsubsurface material not properly supported\n");
  } else if (command.type == "subsurface") {
    get_scalar(command.values, "Kr", material->specular, 1);
    get_scalar(command.values, "Kt", material->transmission, 1);
    material->color = {1, 1, 1};
    get_scalar(command.values, "eta", material->ior, 1.5);
    material->roughness = 0;
    get_roughness(command.values, material->roughness, 0);
    auto scale = 1.0f;
    get_pbrt_value(command.values, "scale", scale);
    material->volscale = 1 / scale;
    auto sigma_a = zero3f, sigma_s = zero3f;
    auto sigma_a_tex = ""s, sigma_s_tex = ""s;
    get_texture(command.values, "sigma_a", sigma_a, sigma_a_tex,
        vec3f{0011, .0024, .014});
    get_texture(command.values, "sigma_prime_s", sigma_s, sigma_s_tex,
        vec3f{2.55, 3.12, 3.77});
    material->volmeanfreepath = 1 / (sigma_a + sigma_s);
    material->volscatter      = sigma_s / (sigma_a + sigma_s);
    if (verbose) printf("subsurface material not properly supported\n");
  } else if (command.type == "mix") {
    auto namedmaterial1 = ""s, namedmaterial2 = ""s;
    get_pbrt_value(command.values, "namedmaterial1", namedmaterial1);
    get_pbrt_value(command.values, "namedmaterial2", namedmaterial2);
    auto matname = (!namedmaterial1.empty()) ? namedmaterial1 : namedmaterial2;
    auto matit   = material_map.find(matname);
    if (matit == material_map.end())
      throw std::invalid_argument("cannot find material " + matname);
    auto saved_name = material->name;
    *material       = *matit->second;
    material->name  = saved_name;
    if (verbose) printf("mix material not properly supported\n");
  } else if (command.type == "fourier") {
    auto bsdffile = ""s;
    get_pbrt_value(command.values, "bsdffile", bsdffile);
    if (bsdffile.rfind("/") != string::npos)
      bsdffile = bsdffile.substr(bsdffile.rfind("/") + 1);
    if (bsdffile == "paint.bsdf") {
      material->color     = {0.6f, 0.6f, 0.6f};
      material->specular  = 1;
      material->ior       = 1.5;
      material->roughness = 0.2;
    } else if (bsdffile == "ceramic.bsdf") {
      material->color     = {0.6f, 0.6f, 0.6f};
      material->specular  = 1;
      material->ior       = 1.5;
      material->roughness = 0.25;
    } else if (bsdffile == "leather.bsdf") {
      material->color     = {0.6f, 0.57f, 0.48f};
      material->specular  = 1;
      material->ior       = 1.5;
      material->roughness = 0.3;
    } else if (bsdffile == "coated_copper.bsdf") {
      auto eta            = vec3f{0.2004376970f, 0.9240334304f, 1.1022119527f};
      auto etak           = vec3f{3.9129485033f, 2.4528477015f, 2.1421879552f};
      material->color     = eta_to_reflectivity(eta, etak);
      material->metallic  = 1;
      material->roughness = 0.01;
    } else if (bsdffile == "roughglass_alpha_0.2.bsdf") {
      material->color        = {1, 1, 1};
      material->specular     = 1;
      material->ior          = 1.5;
      material->transmission = 1;
      material->roughness    = 0.2;
    } else if (bsdffile == "roughgold_alpha_0.2.bsdf") {
      auto eta            = vec3f{0.1431189557f, 0.3749570432f, 1.4424785571f};
      auto etak           = vec3f{3.9831604247f, 2.3857207478f, 1.6032152899f};
      material->color     = eta_to_reflectivity(eta, etak);
      material->metallic  = 1;
      material->roughness = 0.2;
    } else {
      throw std::invalid_argument{"unknown bsdffile " + bsdffile};
    }
  } else {
    throw std::invalid_argument{"unknown material type " + command.type};
  }
  return material;
}

// Make a triangle shape from a quad grid
template <typename PositionFunc, typename NormalFunc>
static void make_pbrt_shape(vector<vec3i>& triangles, vector<vec3f>& positions,
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
static void make_pbrt_sphere(vector<vec3i>& triangles, vector<vec3f>& positions,
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
static void make_pbrt_disk(vector<vec3i>& triangles, vector<vec3f>& positions,
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
static void make_pbrt_quad(vector<vec3i>& triangles, vector<vec3f>& positions,
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
static shared_ptr<pbrt_shape> convert_shape(const pbrt_command& command,
    const string& filename, const string& ply_dirname, bool verbose = false) {
  auto shape   = make_shared<pbrt_shape>();
  shape->frame = command.frame;
  shape->frend = command.frend;
  if (command.type == "trianglemesh") {
    shape->positions = {};
    shape->normals   = {};
    shape->texcoords = {};
    shape->triangles = {};
    get_pbrt_value(command.values, "P", shape->positions);
    get_pbrt_value(command.values, "N", shape->normals);
    get_pbrt_value(command.values, "uv", shape->texcoords);
    for (auto& uv : shape->texcoords) uv.y = (1 - uv.y);
    get_pbrt_value(command.values, "indices", shape->triangles);
  } else if (command.type == "loopsubdiv") {
    shape->positions = {};
    shape->triangles = {};
    get_pbrt_value(command.values, "P", shape->positions);
    get_pbrt_value(command.values, "indices", shape->triangles);
    shape->normals.resize(shape->positions.size());
    // compute_normals(shape->normals, shape->triangles, shape->positions);
  } else if (command.type == "plymesh") {
    shape->filename_ = ""s;
    get_pbrt_value(command.values, "filename", shape->filename_);
    try {
      auto ply = make_shared<ply_model>();
      load_ply(ply_dirname + shape->filename_, ply);
      shape->positions = get_positions(ply);
      shape->normals   = get_normals(ply);
      shape->texcoords = get_texcoords(ply);
      shape->triangles = get_triangles(ply);
    } catch (std::exception& e) {
      throw std::runtime_error{
          filename + ": error in resource (" + e.what() + ")"};
    }
  } else if (command.type == "sphere") {
    auto radius = 1.0f;
    get_pbrt_value(command.values, "radius", radius);
    make_pbrt_sphere(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 16}, radius);
  } else if (command.type == "disk") {
    auto radius = 1.0f;
    get_pbrt_value(command.values, "radius", radius);
    make_pbrt_disk(shape->triangles, shape->positions, shape->normals,
        shape->texcoords, {32, 1}, radius);
  } else {
    throw std::invalid_argument{"unknown shape " + command.type};
  }
  return shape;
}

// Convert pbrt arealights
static shared_ptr<pbrt_arealight> convert_arealight(
    const pbrt_command& command, bool verbose = false) {
  auto light  = make_shared<pbrt_arealight>();
  light->name = command.name;
  if (command.type == "diffuse") {
    auto l = vec3f{1}, scale = vec3f{1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    light->emission = l * scale;
  } else {
    throw std::invalid_argument{"unknown arealight " + command.type};
  }
  return light;
}

// Convert pbrt lights
static shared_ptr<pbrt_light> convert_light(
    const pbrt_command& command, bool verbose = false) {
  auto light   = make_shared<pbrt_light>();
  light->frame = command.frame;
  light->frend = command.frend;
  if (command.type == "distant") {
    auto l = vec3f{1}, scale = vec3f{1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    light->emission = l * scale;
    light->from     = zero3f;
    light->to       = vec3f{0, 0, 1};
    get_pbrt_value(command.values, "from", light->from);
    get_pbrt_value(command.values, "to", light->to);
    light->distant       = true;
    auto distant_dist    = 100;
    auto size            = distant_dist * sin(5 * pif / 180);
    light->area_emission = light->emission * (distant_dist * distant_dist) /
                           (size * size);
    light->area_frame =
        light->frame *
        lookat_frame(normalize(light->from - light->to) * distant_dist, zero3f,
            {0, 1, 0}, true);
    light->area_frend =
        light->frend *
        lookat_frame(normalize(light->from - light->to) * distant_dist, zero3f,
            {0, 1, 0}, true);
    auto texcoords = vector<vec2f>{};
    make_pbrt_quad(light->area_triangles, light->area_positions,
        light->area_normals, texcoords, {4, 2}, size);
  } else if (command.type == "point" || command.type == "goniometric" ||
             command.type == "spot") {
    auto i = vec3f{1}, scale = vec3f{1};
    get_pbrt_value(command.values, "I", i);
    get_pbrt_value(command.values, "scale", scale);
    light->emission = i * scale;
    light->from     = zero3f;
    get_pbrt_value(command.values, "from", light->from);
    light->area_emission = light->emission;
    light->area_frame    = light->frame * translation_frame(light->from);
    light->area_frend    = light->frend * translation_frame(light->from);
    auto texcoords       = vector<vec2f>{};
    make_pbrt_sphere(light->area_triangles, light->area_positions,
        light->area_normals, texcoords, {4, 2}, 0.0025f);
  } else {
    throw std::invalid_argument{"unknown light " + command.type};
  }
  return light;
}

static shared_ptr<pbrt_environment> convert_environment(
    const pbrt_command& command, bool verbose = false) {
  auto environment   = make_shared<pbrt_environment>();
  environment->frame = command.frame;
  environment->frend = command.frend;
  environment->frame = environment->frame *
                       frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  environment->frend = environment->frend *
                       frame3f{{1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 0, 0}};
  if (command.type == "infinite") {
    auto l = vec3f{1}, scale = vec3f{1};
    get_pbrt_value(command.values, "L", l);
    get_pbrt_value(command.values, "scale", scale);
    environment->emission     = scale * l;
    environment->emission_map = ""s;
    get_pbrt_value(command.values, "mapname", environment->emission_map);
  } else {
    throw std::invalid_argument{"unknown light " + command.type};
  }
  return environment;
}

// pbrt stack ctm
struct pbrt_stack_element {
  frame3f                    transform_start        = identity3x4f;
  frame3f                    transform_end          = identity3x4f;
  shared_ptr<pbrt_material>  material               = nullptr;
  shared_ptr<pbrt_arealight> arealight              = nullptr;
  shared_ptr<pbrt_medium>    interior               = nullptr;
  shared_ptr<pbrt_medium>    exterior               = nullptr;
  bool                       reverse                = false;
  bool                       active_transform_start = true;
  bool                       active_transform_end   = true;
};

// pbrt parsing context
struct pbrt_context {
  vector<pbrt_stack_element>                            stack      = {};
  unordered_map<string, pbrt_stack_element>             coordsys   = {};
  unordered_map<string, vector<shared_ptr<pbrt_shape>>> objects    = {};
  string                                                cur_object = "";
  vec2i film_resolution                                            = {512, 512};
};

// load pbrt
void load_pbrt(const string& filename, shared_ptr<pbrt_model> pbrt,
    pbrt_context&                                     ctx,
    unordered_map<string, shared_ptr<pbrt_material>>& material_map,
    unordered_map<string, shared_ptr<pbrt_medium>>&   medium_map,
    unordered_map<string, pbrt_texture>&              texture_map,
    const string&                                     ply_dirname) {
  // open file
  auto fs = fopen(filename.c_str(), "rt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // helpers
  auto set_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = xform;
    if (ctx.active_transform_end) ctx.transform_end = xform;
  };
  auto concat_transform = [](pbrt_stack_element& ctx, const frame3f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= xform;
    if (ctx.active_transform_end) ctx.transform_end *= xform;
  };

  // throw helpers
  auto throw_parse_error = [filename]() {
    throw std::runtime_error{filename + ": parse error"};
  };
  auto throw_read_error = [filename]() {
    throw std::runtime_error{filename + ": read error"};
  };

  // init stack
  if (ctx.stack.empty()) ctx.stack.emplace_back();

  // parse command by command
  auto line = ""s;
  while (read_pbrt_cmdline(fs, line)) {
    auto str = string_view{line};
    // get command
    auto cmd = ""s;
    if (!parse_pbrt_command(str, cmd)) throw_parse_error();
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
      if (!parse_pbrt_param(str, ctx.cur_object)) throw_parse_error();
      ctx.objects[ctx.cur_object] = {};
    } else if (cmd == "ObjectEnd") {
      ctx.stack.pop_back();
      ctx.cur_object = "";
    } else if (cmd == "ObjectInstance") {
      auto object = ""s;
      if (!parse_pbrt_param(str, object)) throw_parse_error();
      if (ctx.objects.find(object) == ctx.objects.end())
        throw std::runtime_error{filename + ": parse error [unknown object]"};
      for (auto shape : ctx.objects.at(object)) {
        shape->instances.push_back(ctx.stack.back().transform_start);
        shape->instaends.push_back(ctx.stack.back().transform_end);
      }
    } else if (cmd == "ActiveTransform") {
      auto name = ""s;
      if (!parse_pbrt_command(str, name)) throw_parse_error();
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
      if (!parse_pbrt_param(str, xf)) throw_parse_error();
      set_transform(ctx.stack.back(), frame3f{xf});
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      if (!parse_pbrt_param(str, xf)) throw_parse_error();
      concat_transform(ctx.stack.back(), frame3f{xf});
    } else if (cmd == "Scale") {
      auto v = zero3f;
      if (!parse_pbrt_param(str, v)) throw_parse_error();
      concat_transform(ctx.stack.back(), scaling_frame(v));
    } else if (cmd == "Translate") {
      auto v = zero3f;
      if (!parse_pbrt_param(str, v)) throw_parse_error();
      concat_transform(ctx.stack.back(), translation_frame(v));
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      if (!parse_pbrt_param(str, v)) throw_parse_error();
      concat_transform(
          ctx.stack.back(), rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      if (!parse_pbrt_param(str, from)) throw_parse_error();
      if (!parse_pbrt_param(str, to)) throw_parse_error();
      if (!parse_pbrt_param(str, up)) throw_parse_error();
      auto frame = lookat_frame(from, to, up, true);
      concat_transform(ctx.stack.back(), inverse(frame));
    } else if (cmd == "ReverseOrientation") {
      ctx.stack.back().reverse = !ctx.stack.back().reverse;
    } else if (cmd == "CoordinateSystem") {
      auto name = ""s;
      if (!parse_pbrt_param(str, name)) throw_parse_error();
      ctx.coordsys[name].transform_start = ctx.stack.back().transform_start;
      ctx.coordsys[name].transform_end   = ctx.stack.back().transform_end;
    } else if (cmd == "CoordSysTransform") {
      auto name = ""s;
      if (!parse_pbrt_param(str, name)) throw_parse_error();
      if (ctx.coordsys.find(name) != ctx.coordsys.end()) {
        ctx.stack.back().transform_start =
            ctx.coordsys.at(name).transform_start;
        ctx.stack.back().transform_end = ctx.coordsys.at(name).transform_end;
      }
    } else if (cmd == "Integrator") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
    } else if (cmd == "Sampler") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
    } else if (cmd == "PixelFilter") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
    } else if (cmd == "Film") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      auto cfilm          = convert_film(command);
      ctx.film_resolution = cfilm->resolution;
    } else if (cmd == "Accelerator") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
    } else if (cmd == "Camera") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto camera   = convert_camera(command, ctx.film_resolution);
      pbrt->cameras.push_back(camera);

    } else if (cmd == "Texture") {
      auto command  = pbrt_command{};
      auto comptype = ""s;
      if (!parse_pbrt_param(str, command.name)) throw_parse_error();
      if (!parse_pbrt_param(str, comptype)) throw_parse_error();
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      texture_map[command.name] = {};
      convert_texture(texture_map[command.name], command, texture_map);
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        command     = pbrt_command{};
      command.name            = "material_" + std::to_string(material_id++);
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      if (command.type == "") {
        ctx.stack.back().material = nullptr;
      } else {
        auto material = convert_material(command, material_map, texture_map);
        pbrt->materials.push_back(material);
        ctx.stack.back().material = material;
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.name)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (value.name == "type") command.type = value.value1s;
      auto material = convert_material(command, material_map, texture_map);
      pbrt->materials.push_back(material);
      material_map[command.name] = material;
    } else if (cmd == "NamedMaterial") {
      auto name = ""s;
      if (!parse_pbrt_param(str, name)) throw_parse_error();
      ctx.stack.back().material = material_map.at(name);
    } else if (cmd == "Shape") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      auto shape    = convert_shape(command, filename, ply_dirname);
      pbrt->shapes.push_back(shape);
      shape->material  = ctx.stack.back().material;
      shape->arealight = ctx.stack.back().arealight;
      if (ctx.cur_object != "") {
        ctx.objects[ctx.cur_object].push_back(shape);
      }
    } else if (cmd == "AreaLightSource") {
      static auto arealight_id = 0;
      auto        command      = pbrt_command{};
      command.name             = "arealight_" + std::to_string(arealight_id++);
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.frame  = ctx.stack.back().transform_start;
      command.frend  = ctx.stack.back().transform_end;
      auto arealight = convert_arealight(command);
      pbrt->arealights.push_back(arealight);
      ctx.stack.back().arealight = arealight;
    } else if (cmd == "LightSource") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.type)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.frame = ctx.stack.back().transform_start;
      command.frend = ctx.stack.back().transform_end;
      if (command.type == "infinite") {
        auto environment = convert_environment(command);
        pbrt->environments.push_back(environment);
      } else {
        auto light = convert_light(command);
        pbrt->lights.push_back(light);
      }
    } else if (cmd == "MakeNamedMedium") {
      auto command = pbrt_command{};
      if (!parse_pbrt_param(str, command.name)) throw_parse_error();
      if (!parse_pbrt_params(str, command.values)) throw_parse_error();
      command.type = "";
      for (auto& value : command.values)
        if (command.name == "type") command.type = value.value1s;
      auto medium = pbrt->mediums.emplace_back(make_shared<pbrt_medium>());
      medium_map[command.name] = medium;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      if (!parse_pbrt_param(str, interior)) throw_parse_error();
      if (!parse_pbrt_param(str, exterior)) throw_parse_error();
      ctx.stack.back().interior = medium_map.at(interior);
      ctx.stack.back().exterior = medium_map.at(exterior);
    } else if (cmd == "Include") {
      auto includename = ""s;
      if (!parse_pbrt_param(str, includename)) throw_parse_error();
      try {
        load_pbrt(fs::path(filename).parent_path() / includename, pbrt, ctx,
            material_map, medium_map, texture_map, ply_dirname);
      } catch (std::exception& e) {
        throw std::runtime_error{
            filename + ": error in resource (" + e.what() + ")"};
      }
    } else {
      throw std::runtime_error{filename + ": parse error [unknown command]"};
    }
  }
}

// Read obj
shared_ptr<pbrt_model> make_pbrt() { return make_shared<pbrt_model>(); }
shared_ptr<pbrt_model> load_pbrt(const string& filename) {
  auto pbrt = make_pbrt();
  load_pbrt(filename, pbrt);
  return pbrt;
}

// load pbrt
void load_pbrt(const string& filename, shared_ptr<pbrt_model> pbrt) {
  auto ctx          = pbrt_context{};
  auto material_map = unordered_map<string, shared_ptr<pbrt_material>>{
      {"", {}}};
  auto medium_map  = unordered_map<string, shared_ptr<pbrt_medium>>{{"", {}}};
  auto texture_map = unordered_map<string, pbrt_texture>{{"", {}}};
  auto dirname     = fs::path(filename).parent_path().string();
  if (dirname != "") dirname += "/";
  load_pbrt(
      filename, pbrt, ctx, material_map, medium_map, texture_map, dirname);
}

static void format_value(string& str, const pbrt_value& value) {
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
      format_value(str, value);
    }
    str += " ]";
  };

  format_values(str, "\"{} {}\" ", type_labels.at(value.type), value.name);
  switch (value.type) {
    case pbrt_value_type::real:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1f);
      } else {
        format_value(str, value.value1f);
      }
      break;
    case pbrt_value_type::integer:
      if (!value.vector1f.empty()) {
        format_vector(str, value.vector1i);
      } else {
        format_value(str, value.value1i);
      }
      break;
    case pbrt_value_type::boolean:
      format_values(str, "\"{}\"", value.value1b ? "true" : "false");
      break;
    case pbrt_value_type::string:
    case pbrt_value_type::texture:
      format_values(str, "\"{}\"", value.value1s);
      break;
    case pbrt_value_type::point:
    case pbrt_value_type::vector:
    case pbrt_value_type::normal:
    case pbrt_value_type::color:
      if (!value.vector3f.empty()) {
        format_vector(str, value.vector3f);
      } else {
        format_values(str, "[ {} ]", value.value3f);
      }
      break;
    case pbrt_value_type::spectrum: format_vector(str, value.vector1f); break;
    case pbrt_value_type::point2:
    case pbrt_value_type::vector2:
      if (!value.vector2f.empty()) {
        format_vector(str, value.vector2f);
      } else {
        format_values(str, "[ {} ]", value.value2f);
      }
      break;
  }
}

static void format_value(string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

void save_pbrt(
    const string& filename, shared_ptr<pbrt_model> pbrt, bool ply_meshes) {
  // open file
  auto fs = fopen(filename.c_str(), "wt");
  if (!fs) throw std::runtime_error{filename + ": file not found"};
  auto fs_guard = std::unique_ptr<FILE, decltype(&fclose)>{fs, fclose};

  // throw helpers
  auto throw_write_error = [filename]() {
    throw std::runtime_error{filename + ": write error"};
  };

  // save comments
  if (!format_values(fs, "#\n")) throw_write_error();
  if (!format_values(fs, "# Written by Yocto/GL\n")) throw_write_error();
  if (!format_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    throw_write_error();
  if (!format_values(fs, "#\n\n")) throw_write_error();
  for (auto& comment : pbrt->comments) {
    if (!format_values(fs, "# {}\n", comment)) throw_write_error();
  }
  if (!format_values(fs, "\n")) throw_write_error();

  for (auto camera : pbrt->cameras) {
    auto command = pbrt_command{};
    command.type = "image";
    command.values.push_back(
        make_pbrt_value("xresolution", camera->resolution.x));
    command.values.push_back(
        make_pbrt_value("yresolution", camera->resolution.y));
    command.values.push_back(make_pbrt_value("filename", "image.exr"s));
    if (!format_values(fs, "Film \"{}\" {}\n", command.type, command.values))
      throw_write_error();
  }

  for (auto camera : pbrt->cameras) {
    auto command  = pbrt_command{};
    command.type  = "perspective";
    command.frame = camera->frame;
    command.values.push_back(make_pbrt_value(
        "fov", 2 * tan(0.036f / (2 * camera->lens)) * 180 / pif));
    if (!format_values(fs, "LookAt {} {} {}\n", command.frame.o,
            command.frame.o - command.frame.z, command.frame.y))
      throw_write_error();
    if (!format_values(fs, "Camera \"{}\" {}\n", command.type, command.values))
      throw_write_error();
  }

  if (!format_values(fs, "\nWorldBegin\n\n")) throw_write_error();

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
    if (!format_values(fs, "AttributeBegin\n")) throw_write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)command.frame))
      throw_write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      throw_write_error();
    if (!format_values(fs, "AttributeEnd\n")) throw_write_error();
  }

  for (auto environment : pbrt->environments) {
    auto command  = pbrt_command{};
    command.frame = environment->frame;
    command.type  = "infinite";
    command.values.push_back(make_pbrt_value("L", environment->emission));
    command.values.push_back(
        make_pbrt_value("mapname", environment->emission_map));
    if (!format_values(fs, "AttributeBegin\n")) throw_write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)command.frame))
      throw_write_error();
    if (!format_values(
            fs, "LightSource \"{}\" {}\n", command.type, command.values))
      throw_write_error();
    if (!format_values(fs, "AttributeEnd\n")) throw_write_error();
  }

  auto reflectivity_to_eta = [](const vec3f& reflectivity) {
    return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
  };

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
      if (material->color_map.empty()) {
        command.values.push_back(make_pbrt_value("Kd", material->color));
      } else if (material->color != zero3f) {
        command.values.push_back(make_pbrt_value(
            "Kd", material->color_map, pbrt_value_type::texture));
      }
      if (material->specular != 0) {
        command.values.push_back(
            make_pbrt_value("Ks", vec3f{material->specular}));
        command.values.push_back(
            make_pbrt_value("roughness", pow(material->roughness, 2)));
        command.values.push_back(make_pbrt_value("eta", material->ior));
        command.values.push_back(make_pbrt_value("remaproughness", false));
      }
      if (material->transmission != 0) {
        command.values.push_back(
            make_pbrt_value("Kt", vec3f{material->transmission}));
      }
      if (!material->opacity_map.empty()) {
        command.values.push_back(make_pbrt_value(
            "opacity", material->opacity_map, pbrt_value_type::texture));
      } else if (material->opacity != 1) {
        command.values.push_back(make_pbrt_value("opacity", material->opacity));
      }
    }
    if (!format_values(fs,
            "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n",
            material->name, command.type, command.values))
      throw_write_error();
  }

  auto object_id = 0;
  for (auto shape : pbrt->shapes) {
    auto command  = pbrt_command{};
    command.frame = shape->frame;
    if (ply_meshes) {
      command.type = "plymesh";
      command.values.push_back(make_pbrt_value("filename", shape->filename_));
    } else {
      command.type = "trianglemesh";
      command.values.push_back(make_pbrt_value("indices", shape->triangles));
      command.values.push_back(
          make_pbrt_value("P", shape->positions, pbrt_value_type::point));
      if (!shape->normals.empty())
        command.values.push_back(
            make_pbrt_value("N", shape->triangles, pbrt_value_type::normal));
      if (!shape->texcoords.empty())
        command.values.push_back(make_pbrt_value("uv", shape->texcoords));
    }
    auto acommand = pbrt_command{};
    if (shape->arealight) {
      acommand.type = "diffuse";
      acommand.values.push_back(
          make_pbrt_value("L", shape->arealight->emission));
    }
    if (ply_meshes) {
      try {
        auto ply = make_shared<ply_model>();
        add_positions(ply, shape->positions);
        add_normals(ply, shape->normals);
        add_texcoords(ply, shape->texcoords);
        add_triangles(ply, shape->triangles);
        save_ply(fs::path(filename).parent_path() / shape->filename_, ply);
      } catch (std::exception& e) {
        throw std::runtime_error{
            filename + ": error in resource (" + e.what() + ")"};
      }
    }
    auto object = "object" + std::to_string(object_id++);
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectBegin \"{}\"\n", object))
        throw_write_error();
    if (!format_values(fs, "AttributeBegin\n")) throw_write_error();
    if (!format_values(fs, "Transform {}\n", (mat4f)shape->frame))
      throw_write_error();
    if (shape->arealight) {
      if (!format_values(fs, "AreaLightSource \"{}\" {}\n", acommand.type,
              acommand.values))
        throw_write_error();
    }
    if (!format_values(fs, "NamedMaterial \"{}\"\n", shape->material->name))
      throw_write_error();
    if (!format_values(fs, "Shape \"{}\" {}\n", command.type, command.values))
      throw_write_error();
    if (!format_values(fs, "AttributeEnd\n")) throw_write_error();
    if (!shape->instances.empty())
      if (!format_values(fs, "ObjectEnd\n")) throw_write_error();
    for (auto& iframe : shape->instances) {
      if (!format_values(fs, "AttributeBegin\n")) throw_write_error();
      if (!format_values(fs, "Transform {}\n", (mat4f)iframe))
        throw_write_error();
      if (!format_values(fs, "ObjectInstance \"{}\"\n", object))
        throw_write_error();
      if (!format_values(fs, "AttributeEnd\n")) throw_write_error();
    }
  }

  if (!format_values(fs, "\nWorldEnd\n\n")) throw_write_error();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE GLTF LOADER
// -----------------------------------------------------------------------------
namespace yocto {

void update_transforms(
    gltf_model& scene, gltf_node& node, const frame3f& parent = identity3x4f) {
  auto frame = parent * node.local * translation_frame(node.translation) *
               rotation_frame(node.rotation) * scaling_frame(node.scale);
  for (auto child : node.children)
    update_transforms(scene, scene.nodes[child], frame);
}

// convert gltf to scene
void load_gltf(const string& filename, gltf_model& scene) {
  // load gltf
  auto params = cgltf_options{};
  memset(&params, 0, sizeof(params));
  auto data   = (cgltf_data*)nullptr;
  auto result = cgltf_parse_file(&params, filename.c_str(), &data);
  if (result != cgltf_result_success) {
    throw std::runtime_error{filename + ": read error"};
  }
  auto gltf = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      data, cgltf_free};
  auto dirname = fs::path(filename).parent_path().string();
  if (dirname != "") dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success) {
    throw std::runtime_error(filename + ": error reading buffers");
  }

  // convert textures
  auto _startswith = [](string_view str, string_view substr) {
    if (str.size() < substr.size()) return false;
    return str.substr(0, substr.size()) == substr;
  };
  auto imap = unordered_map<cgltf_image*, int>{};
  for (auto tid = 0; tid < gltf->images_count; tid++) {
    auto gimg        = &gltf->images[tid];
    auto texture     = gltf_texture{};
    texture.name     = gimg->name ? gimg->name : "";
    texture.filename = (_startswith(gimg->uri, "data:"))
                           ? string("[glTF-static inline].png")
                           : gimg->uri;
    scene.textures.push_back(texture);
    imap[gimg] = tid;
  }

  // add a texture
  auto get_texture = [&imap](const cgltf_texture_view& ginfo) {
    if (!ginfo.texture || !ginfo.texture->image) return -1;
    auto gtxt = ginfo.texture;
    return imap.at(gtxt->image);
  };

  // convert materials
  auto mmap = unordered_map<cgltf_material*, int>{{nullptr, -1}};
  for (auto mid = 0; mid < gltf->materials_count; mid++) {
    auto gmat             = &gltf->materials[mid];
    mmap[gmat]            = mid;
    auto& material        = scene.materials.emplace_back();
    material.name         = gmat->name ? gmat->name : "";
    material.emission     = {gmat->emissive_factor[0], gmat->emissive_factor[1],
        gmat->emissive_factor[2]};
    material.emission_tex = get_texture(gmat->emissive_texture);
    if (gmat->has_pbr_specular_glossiness) {
      material.has_specgloss = true;
      auto gsg               = &gmat->pbr_specular_glossiness;
      material.sg_diffuse    = vec4f{gsg->diffuse_factor[0],
          gsg->diffuse_factor[1], gsg->diffuse_factor[2],
          gsg->diffuse_factor[3]};
      material.sg_specular = {gsg->specular_factor[0], gsg->specular_factor[1],
          gsg->specular_factor[2]};
      material.sg_glossiness   = gsg->glossiness_factor;
      material.sg_diffuse_tex  = get_texture(gsg->diffuse_texture);
      material.sg_specular_tex = get_texture(gsg->specular_glossiness_texture);
    } else if (gmat->has_pbr_metallic_roughness) {
      material.has_metalrough  = true;
      auto gmr                 = &gmat->pbr_metallic_roughness;
      material.mr_base         = vec4f{gmr->base_color_factor[0],
          gmr->base_color_factor[1], gmr->base_color_factor[2],
          gmr->base_color_factor[3]};
      material.mr_metallic     = gmr->metallic_factor;
      material.mr_roughness    = gmr->roughness_factor;
      material.mr_base_tex     = get_texture(gmr->base_color_texture);
      material.mr_metallic_tex = get_texture(gmr->metallic_roughness_texture);
    }
    material.normal_tex = get_texture(gmat->normal_texture);
  }

  // get values from accessors
  auto accessor_values =
      [](const cgltf_accessor* gacc,
          bool normalize = false) -> vector<std::array<double, 4>> {
    auto gview       = gacc->buffer_view;
    auto data        = (byte*)gview->buffer->data;
    auto offset      = gacc->offset + gview->offset;
    auto stride      = gview->stride;
    auto compTypeNum = gacc->component_type;
    auto count       = gacc->count;
    auto type        = gacc->type;
    auto ncomp       = 0;
    if (type == cgltf_type_scalar) ncomp = 1;
    if (type == cgltf_type_vec2) ncomp = 2;
    if (type == cgltf_type_vec3) ncomp = 3;
    if (type == cgltf_type_vec4) ncomp = 4;
    auto compSize = 1;
    if (compTypeNum == cgltf_component_type_r_16 ||
        compTypeNum == cgltf_component_type_r_16u) {
      compSize = 2;
    }
    if (compTypeNum == cgltf_component_type_r_32u ||
        compTypeNum == cgltf_component_type_r_32f) {
      compSize = 4;
    }
    if (!stride) stride = compSize * ncomp;
    auto vals = vector<std::array<double, 4>>(count, {{0.0, 0.0, 0.0, 1.0}});
    for (auto i = 0; i < count; i++) {
      auto d = data + offset + i * stride;
      for (auto c = 0; c < ncomp; c++) {
        if (compTypeNum == cgltf_component_type_r_8) {  // char
          vals[i][c] = (double)(*(char*)d);
          if (normalize) vals[i][c] /= SCHAR_MAX;
        } else if (compTypeNum == cgltf_component_type_r_8u) {  // byte
          vals[i][c] = (double)(*(byte*)d);
          if (normalize) vals[i][c] /= UCHAR_MAX;
        } else if (compTypeNum == cgltf_component_type_r_16) {  // short
          vals[i][c] = (double)(*(short*)d);
          if (normalize) vals[i][c] /= SHRT_MAX;
        } else if (compTypeNum ==
                   cgltf_component_type_r_16u) {  // unsigned short
          vals[i][c] = (double)(*(unsigned short*)d);
          if (normalize) vals[i][c] /= USHRT_MAX;
        } else if (compTypeNum == cgltf_component_type_r_32u) {  // unsigned int
          vals[i][c] = (double)(*(unsigned int*)d);
          if (normalize) vals[i][c] /= UINT_MAX;
        } else if (compTypeNum == cgltf_component_type_r_32f) {  // float
          vals[i][c] = (*(float*)d);
        }
        d += compSize;
      }
    }
    return vals;
  };

  // convert meshes
  auto smap = unordered_map<cgltf_mesh*, int>{{nullptr, -1}};
  for (auto mid = 0; mid < gltf->meshes_count; mid++) {
    auto gmesh  = &gltf->meshes[mid];
    smap[gmesh] = mid;
    auto& mesh  = scene.meshes.emplace_back();
    mesh.name   = gmesh->name ? gmesh->name : "";
    for (auto sid = 0; sid < gmesh->primitives_count; sid++) {
      auto gprim = &gmesh->primitives[sid];
      if (!gprim->attributes_count) continue;
      auto& shape = mesh.primitives.emplace_back();
      for (auto aid = 0; aid < gprim->attributes_count; aid++) {
        auto gattr    = &gprim->attributes[aid];
        auto semantic = string(gattr->name ? gattr->name : "");
        auto gacc     = gattr->data;
        auto vals     = accessor_values(gacc);
        if (semantic == "POSITION") {
          shape.positions.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.positions.push_back(
                {(float)vals[i][0], (float)vals[i][1], (float)vals[i][2]});
        } else if (semantic == "NORMAL") {
          shape.normals.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.normals.push_back(
                {(float)vals[i][0], (float)vals[i][1], (float)vals[i][2]});
        } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
          shape.texcoords.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.texcoords.push_back({(float)vals[i][0], (float)vals[i][1]});
        } else if (semantic == "COLOR" || semantic == "COLOR_0") {
          shape.colors.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.colors.push_back({(float)vals[i][0], (float)vals[i][1],
                (float)vals[i][2], (float)vals[i][3]});
        } else if (semantic == "TANGENT") {
          shape.tangents.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.tangents.push_back({(float)vals[i][0], (float)vals[i][1],
                (float)vals[i][2], (float)vals[i][3]});
          for (auto& t : shape.tangents) t.w = -t.w;
        } else if (semantic == "RADIUS") {
          shape.radius.reserve(vals.size());
          for (auto i = 0; i < vals.size(); i++)
            shape.radius.push_back((float)vals[i][0]);
        } else {
          // ignore
        }
      }
      // indices
      if (!gprim->indices) {
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.reserve(shape.positions.size() / 3);
          for (auto i = 0; i < shape.positions.size() / 3; i++)
            shape.triangles.push_back({i * 3 + 0, i * 3 + 1, i * 3 + 2});
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.reserve(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles.push_back({0, i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.reserve(shape.positions.size() - 2);
          for (auto i = 2; i < shape.positions.size(); i++)
            shape.triangles.push_back({i - 2, i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.reserve(shape.positions.size() / 2);
          for (auto i = 0; i < shape.positions.size() / 2; i++)
            shape.lines.push_back({i * 2 + 0, i * 2 + 1});
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.reserve(shape.positions.size());
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines.push_back({i - 1, i});
          shape.lines.back() = {(int)shape.positions.size() - 1, 0};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.reserve(shape.positions.size() - 1);
          for (auto i = 1; i < shape.positions.size(); i++)
            shape.lines.push_back({i - 1, i});
        } else if (gprim->type == cgltf_primitive_type_points) {
          // points
          throw std::runtime_error("points not supported");
        } else {
          throw std::runtime_error("unknown primitive type");
        }
      } else {
        auto indices = accessor_values(gprim->indices);
        if (gprim->type == cgltf_primitive_type_triangles) {
          shape.triangles.reserve(indices.size() / 3);
          for (auto i = 0; i < indices.size() / 3; i++)
            shape.triangles.push_back({(int)indices[i * 3 + 0][0],
                (int)indices[i * 3 + 1][0], (int)indices[i * 3 + 2][0]});
        } else if (gprim->type == cgltf_primitive_type_triangle_fan) {
          shape.triangles.reserve(indices.size() - 2);
          for (auto i = 2; i < indices.size(); i++)
            shape.triangles.push_back({(int)indices[0][0],
                (int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_triangle_strip) {
          shape.triangles.reserve(indices.size() - 2);
          for (auto i = 2; i < indices.size(); i++)
            shape.triangles.push_back({(int)indices[i - 2][0],
                (int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_lines) {
          shape.lines.reserve(indices.size() / 2);
          for (auto i = 0; i < indices.size() / 2; i++)
            shape.lines.push_back(
                {(int)indices[i * 2 + 0][0], (int)indices[i * 2 + 1][0]});
        } else if (gprim->type == cgltf_primitive_type_line_loop) {
          shape.lines.reserve(indices.size());
          for (auto i = 1; i < indices.size(); i++)
            shape.lines.push_back({(int)indices[i - 1][0], (int)indices[i][0]});
          shape.lines.back() = {
              (int)indices[indices.size() - 1][0], (int)indices[0][0]};
        } else if (gprim->type == cgltf_primitive_type_line_strip) {
          shape.lines.reserve(indices.size() - 1);
          for (auto i = 1; i < indices.size(); i++)
            shape.lines.push_back({(int)indices[i - 1][0], (int)indices[i][0]});
        } else if (gprim->type == cgltf_primitive_type_points) {
          throw std::runtime_error("points not supported");
        } else {
          throw std::runtime_error("unknown primitive type");
        }
      }
    }
  }

  // convert cameras
  auto cmap = unordered_map<cgltf_camera*, int>{{nullptr, -1}};
  for (auto cid = 0; cid < gltf->cameras_count; cid++) {
    auto gcam    = &gltf->cameras[cid];
    cmap[gcam]   = cid;
    auto& camera = scene.cameras.emplace_back();
    camera.name  = gcam->name ? gcam->name : "";
    camera.ortho = gcam->type == cgltf_camera_type_orthographic;
    if (camera.ortho) {
      // throw std::runtime_error("orthographic not supported well");
      auto ortho    = &gcam->orthographic;
      camera.yfov   = ortho->ymag;
      camera.aspect = ortho->xmag / ortho->ymag;
    } else {
      auto persp    = &gcam->perspective;
      camera.yfov   = persp->yfov;
      camera.aspect = persp->aspect_ratio;
    }
    scene.cameras.push_back(camera);
    cmap[gcam] = (int)scene.cameras.size() - 1;
  }

  // convert nodes
  auto nmap = unordered_map<cgltf_node*, int>{{nullptr, -1}};
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde  = &gltf->nodes[nid];
    nmap[gnde] = nid;
    auto& node = scene.nodes.emplace_back();
    node.name  = gnde->name ? gnde->name : "";
    if (gnde->camera) node.camera = cmap.at(gnde->camera);
    if (gnde->mesh) node.mesh = smap.at(gnde->mesh);
    if (gnde->has_translation) {
      node.translation = {
          gnde->translation[0], gnde->translation[1], gnde->translation[2]};
    }
    if (gnde->has_rotation) {
      node.rotation = {gnde->rotation[0], gnde->rotation[1], gnde->rotation[2],
          gnde->rotation[3]};
    }
    if (gnde->has_scale) {
      node.scale = {gnde->scale[0], gnde->scale[1], gnde->scale[2]};
    }
    if (gnde->has_matrix) {
      auto m     = gnde->matrix;
      node.local = frame3f(
          mat4f{{m[0], m[1], m[2], m[3]}, {m[4], m[5], m[6], m[7]},
              {m[8], m[9], m[10], m[11]}, {m[12], m[13], m[14], m[15]}});
    }
  }

  // set up parent pointers
  for (auto nid = 0; nid < gltf->nodes_count; nid++) {
    auto gnde = &gltf->nodes[nid];
    if (!gnde->children_count) continue;
    for (auto cid = 0; cid < gnde->children_count; cid++) {
      scene.nodes[nid].children.push_back(nmap.at(gnde->children[cid]));
      scene.nodes[nmap.at(gnde->children[cid])].parent = nid;
    }
  }

  // set up scenes
  for (auto sid = 0; sid < gltf->scenes_count; sid++) {
    auto  gscn = &gltf->scenes[sid];
    auto& scn  = scene.scenes.emplace_back();
    scn.name   = gscn->name ? gscn->name : "";
    for (auto nid = 0; nid < gscn->nodes_count; nid++) {
      scn.nodes.push_back(nmap.at(gscn->nodes[nid]));
    }
  }

  // update transforms
  for (auto& node : scene.nodes)
    if (node.parent < 0) update_transforms(scene, node);

#if 0
  // hasher for later
  struct sampler_map_hash {
    size_t operator()(
        const pair<cgltf_animation_sampler*, cgltf_animation_path_type>& value)
        const {
      auto hasher1 = std::hash<cgltf_animation_sampler*>();
      auto hasher2 = std::hash<int>();
      auto h       = (size_t)0;
      h ^= hasher1(value.first) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher2(value.second) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };

  // convert animations
  for (auto gid = 0; gid < gltf->animations_count; gid++) {
    auto ganm = &gltf->animations[gid];
    auto aid  = 0;
    auto sampler_map =
        unordered_map<pair<cgltf_animation_sampler*, cgltf_animation_path_type>,
            int, sampler_map_hash>();
    for (auto cid = 0; cid < ganm->channels_count; cid++) {
      auto gchannel = &ganm->channels[cid];
      auto path     = gchannel->target_path;
      if (sampler_map.find({gchannel->sampler, path}) == sampler_map.end()) {
        auto gsampler  = gchannel->sampler;
        auto animation = gltf_animation{};
        animation.uri  = (ganm->name ? ganm->name : "anim") +
                        std::to_string(aid++);
        animation.group = ganm->name ? ganm->name : "";
        auto input_view = accessor_values(gsampler->input);
        animation.times.resize(input_view.size());
        for (auto i = 0; i < input_view.size(); i++)
          animation.times[i] = input_view[i][0];
        switch (gsampler->interpolation) {
          case cgltf_interpolation_type_linear:
            animation.interpolation =
                gltf_animation::interpolation_type::linear;
            break;
          case cgltf_interpolation_type_step:
            animation.interpolation = gltf_animation::interpolation_type::step;
            break;
          case cgltf_interpolation_type_cubic_spline:
            animation.interpolation =
                gltf_animation::interpolation_type::bezier;
            break;
        }
        auto output_view = accessor_values(gsampler->output);
        switch (path) {
          case cgltf_animation_path_type_translation: {
            animation.translations.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.translations.push_back({(float)output_view[i][0],
                  (float)output_view[i][1], (float)output_view[i][2]});
          } break;
          case cgltf_animation_path_type_rotation: {
            animation.rotations.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.rotations.push_back(
                  {(float)output_view[i][0], (float)output_view[i][1],
                      (float)output_view[i][2], (float)output_view[i][3]});
          } break;
          case cgltf_animation_path_type_scale: {
            animation.scales.reserve(output_view.size());
            for (auto i = 0; i < output_view.size(); i++)
              animation.scales.push_back({(float)output_view[i][0],
                  (float)output_view[i][1], (float)output_view[i][2]});
          } break;
          case cgltf_animation_path_type_weights: {
            throw std::runtime_error("weights not supported for now");
                    // // get a node that it refers to
                    // auto ncomp = 0;
                    // auto gnode = gltf->get(gchannel->target->node);
                    // auto gmesh = gltf->get(gnode->mesh);
                    // if (gmesh) {
                    //     for (auto gshp : gmesh->primitives) {
                    //         ncomp = max((int)gshp->targets.size(), ncomp);
                    //     }
                    // }
                    // if (ncomp) {
                    //     auto values = vector<float>();
                    //     values.reserve(output_view.size());
                    //     for (auto i = 0; i < output_view.size(); i++)
                    //         values.push_back(output_view.get(i));
                    //     animation.weights.resize(values.size() / ncomp);
                    //     for (auto i = 0; i < animation.weights.size(); i++) {
                    //         animation.weights[i].resize(ncomp);
                    //         for (auto j = 0; j < ncomp; j++)
                    //             animation.weights[i][j] = values[i * ncomp + j];
                    //     }
                    // }
          } break;
          default: {
            throw std::runtime_error("bad gltf animation");
          }
        }
        sampler_map[{gchannel->sampler, path}] = (int)scene.animations.size();
        scene.animations.push_back(animation);
      }
      scene.animations[sampler_map.at({gchannel->sampler, path})]
          .targets.push_back(nmap.at(gchannel->target_node));
    }
  }
#endif
}

}  // namespace yocto
