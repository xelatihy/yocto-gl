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

#include "yocto_ply.h"

#include <cstdio>
#include <memory>
#include <string_view>
#include <unordered_map>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string_view;
using std::unordered_map;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// string literals
using namespace std::string_literals;

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

inline void remove_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

ply_element::~ply_element() {
  for (auto property : properties) delete property;
}
ply_model::~ply_model() {
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
bool save_ply(const string& filename, ply_model* ply, string& error) {
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
  std::runtime_error("should not have gotten here");
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
    const std::array<string, 2>& properties, vector<vec2f>& values) {
  values.clear();
  auto x = vector<float>{}, y = vector<float>{};
  if (!get_value(ply, element, properties[0], x)) return false;
  if (!get_value(ply, element, properties[1], y)) return false;
  values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return true;
}
bool get_values(ply_model* ply, const string& element,
    const std::array<string, 3>& properties, vector<vec3f>& values) {
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
    const std::array<string, 4>& properties, vector<vec4f>& values) {
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
    const std::array<string, 12>& properties, vector<frame3f>& values) {
  values.clear();
  auto coords = std::array<vector<float>, 12>{};
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

inline vector<vec2f> flip_texcoord(const vector<vec2f>& texcoords) {
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
  if (!add_element(ply, element_name, count)) return nullptr;
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
  if (!values) return false;
  for (auto p = 0; p < nprops; p++) {
    if (!add_property(ply, element, properties[p], count, ply_type::f32, false))
      return false;
    auto prop = get_property(ply, element, properties[p]);
    prop->data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop->data_f32[i] = values[p + i * nprops];
  }
  return true;
}

bool add_value(ply_model* ply, const string& element, const string& property,
    const vector<float>& values) {
  if (values.empty()) return false;
  auto properties = vector{property};
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 2>& properties, const vector<vec2f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 3>& properties, const vector<vec3f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 4>& properties, const vector<vec4f>& values) {
  if (values.empty()) return false;
  return add_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 12>& properties, const vector<frame3f>& values) {
  if (values.empty()) return false;
  return add_values(ply, (float*)values.data(), values.size(), element,
      properties.data(), properties.size());
}

bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vector<int>>& values) {
  if (values.empty()) return false;
  if (!add_property(ply, element, property, values.size(), ply_type::i32, true))
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
  if (!add_property(ply, element, property, sizes.size(), ply_type::i32, true))
    return false;
  auto prop      = get_property(ply, element, property);
  prop->data_i32 = values;
  prop->ldata_u8 = sizes;
  return true;
}
bool add_lists(ply_model* ply, const int* values, size_t count, int size,
    const string& element, const string& property) {
  if (!values) return false;
  if (!add_property(ply, element, property, count, ply_type::i32, true))
    return false;
  auto prop = get_property(ply, element, property);
  prop->data_i32.assign(values, values + count * size);
  prop->ldata_u8.assign(count, size);
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
  return add_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec3i>& values) {
  if (values.empty()) return false;
  return add_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec4i>& values) {
  if (values.empty()) return false;
  return add_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
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
      ply, "vertex", {"u", "v"}, flipv ? flip_texcoord(values) : values);
}
bool add_colors(ply_model* ply, const vector<vec3f>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue"}, values);
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
