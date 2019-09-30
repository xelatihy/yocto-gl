//
// Implementation for Yocto/ModelIO.
//

//
// TODO: fov/aspect lens/film everywhere
// TODO: pbrt design, split elements from approximations
// TODO: add uvdisk and uvsphere to pbrt shape loading
// TODO: remove obj instances
// TODO: remove obj procedurals
// TODO: simplify objx parsing using single line commands
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

#include "yocto_image.h"
#include "yocto_utils.h"

#include <algorithm>
#include <cinttypes>
#include <climits>
#include <cstdarg>
#include <limits>
#include <string_view>

#define CGLTF_IMPLEMENTATION
#include "ext/cgltf.h"

// -----------------------------------------------------------------------------
// FILE AND PROPERTY HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
file_wrapper::file_wrapper(file_wrapper&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
file_wrapper::~file_wrapper() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
void open_file(file_wrapper& fs, const string& filename, const string& mode) {
  close_file(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
file_wrapper open_file(const string& filename, const string& mode) {
  auto fs = file_wrapper{};
  open_file(fs, filename, mode);
  return fs;
}
void close_file(file_wrapper& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;
using namespace std::literals::string_view_literals;

template <typename T>
static inline T swap_endian(T value) {
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

// Read a line
static inline bool read_line(file_wrapper& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

static inline bool is_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_newline(char c) { return c == '\r' || c == '\n'; }
static inline bool is_digit(char c) { return c >= '0' && c <= '9'; }
static inline bool is_alpha(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

static inline void skip_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
}
static inline void trim_whitespace(string_view& str) {
  while (!str.empty() && is_space(str.front())) str.remove_prefix(1);
  while (!str.empty() && is_space(str.back())) str.remove_suffix(1);
}

static inline bool is_whitespace(string_view str) {
  while (!str.empty()) {
    if (!is_space(str.front())) return false;
    str.remove_prefix(1);
  }
  return true;
}

static inline vector<string> split_string(
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

static inline vector<vec2f> flip_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Parse values from a string
static inline void parse_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
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
}
static inline void parse_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_value(str, valuev);
  value = string{valuev};
}
static inline void parse_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_value(str, valuei);
  value = (bool)valuei;
}
static inline void parse_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
static inline void parse_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_value(str, values[i]);
}

static inline void parse_value(string_view& str, vec2f& value) {
  parse_value(str, &value.x, 2);
}
static inline void parse_value(string_view& str, vec3f& value) {
  parse_value(str, &value.x, 3);
}
static inline void parse_value(string_view& str, frame3f& value) {
  parse_value(str, &value.x.x, 12);
}
#ifdef __APPLE__
static inline void parse_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
#endif

// Parse values from a string
template <typename T>
static inline void parse_value_or_empty(string_view& str, T& value) {
  skip_whitespace(str);
  if (str.empty()) {
    value = T{};
  } else {
    parse_value(str, value);
  }
}

// Formats values to string
static inline void format_value(string& str, const string& value) {
  str += value;
}
static inline void format_value(string& str, const char* value) {
  str += value;
}
static inline void format_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static inline void format_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
static inline void format_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static inline void format_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
static inline void format_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static inline void format_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static inline void format_value(string& str, const vec2f& value) {
  char buf[256];
  sprintf(buf, "%g %g", value.x, value.y);
  str += buf;
}
static inline void format_value(string& str, const vec3f& value) {
  char buf[256];
  sprintf(buf, "%g %g %g", value.x, value.y, value.z);
  str += buf;
}
#if 0
static inline void format_value(string& str, const vec2i& value) {
  char buf[256];
  sprintf(buf, "%d %d", value.x, value.y);
  str += buf;
}
static inline void format_value(string& str, const vec3i& value) {
  char buf[256];
  sprintf(buf, "%d %d %d", value.x, value.y, value.z);
  str += buf;
}
#endif
static inline void format_value(string& str, const frame3f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g", value.x.x, value.x.y,
      value.x.z, value.y.x, value.y.y, value.y.z, value.z.x, value.z.y,
      value.z.z, value.o.x, value.o.y, value.o.z);
  str += buf;
}
static inline void format_value(string& str, const mat4f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", value.x.x,
      value.x.y, value.x.z, value.x.w, value.y.x, value.y.y, value.y.z,
      value.y.w, value.z.x, value.z.y, value.z.z, value.z.w, value.w.x,
      value.w.y, value.w.z, value.w.w);
  str += buf;
}

// Foramt to file
static inline void format_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static inline void format_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_value(str, arg);
  format_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
static inline void format_values(
    file_wrapper& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
static inline void format_value(file_wrapper& fs, const T& value) {
  auto str = ""s;
  format_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

static inline void write_text(file_wrapper& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
static inline void write_text(file_wrapper& fs, const char* value) {
  if (fputs(value, fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}

template <typename T>
static inline void write_value(file_wrapper& fs, const T& value) {
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}
template <typename T>
static inline void write_value(
    file_wrapper& fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_endian(value_) : value_;
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}

template <typename T>
static inline void read_value(file_wrapper& fs, T& value) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
}
template <typename T>
static inline void read_value(file_wrapper& fs, T& value, bool big_endian) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
  if (big_endian) value = swap_endian(value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static inline void remove_ply_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Load ply
void load_ply(const string& filename, ply_model& ply) {
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
  ply = {};

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // open file
  auto fs = open_file(filename, "rb");

  // read header ---------------------------------------------
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_ply_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // check magic number
    if (first_line) {
      if (cmd != "ply") throw std::runtime_error{"bad ply file"};
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) throw std::runtime_error{"bad ply file"};
    } else if (cmd == "format") {
      auto fmt = ""sv;
      parse_value(line, fmt);
      if (fmt == "ascii") {
        ply.format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        ply.format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        ply.format = ply_format::binary_big_endian;
      } else {
        throw std::runtime_error{"unknown ply format"};
      }
    } else if (cmd == "comment") {
      skip_whitespace(line);
      ply.comments.push_back(string{line});
    } else if (cmd == "obj_info") {
      skip_whitespace(line);
      // comment is the rest of the line
    } else if (cmd == "element") {
      auto& elem = ply.elements.emplace_back();
      parse_value(line, elem.name);
      parse_value(line, elem.count);
    } else if (cmd == "property") {
      if (ply.elements.empty()) throw std::runtime_error{"bad ply header"};
      auto& prop  = ply.elements.back().properties.emplace_back();
      auto  tname = ""s;
      parse_value(line, tname);
      if (tname == "list") {
        prop.is_list = true;
        parse_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          throw std::runtime_error{"unsupported list size type " + tname};
        parse_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      }
      parse_value(line, prop.name);
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      throw std::runtime_error{"unknown ply command"};
    }
  }

  // check exit
  if (!end_header) throw std::runtime_error{"bad ply header"};

  // allocate data ---------------------------------
  for (auto& element : ply.elements) {
    for (auto& property : element.properties) {
      auto count = property.is_list ? element.count * 3 : element.count;
      switch (property.type) {
        case ply_type::i8: property.data_i8.reserve(count); break;
        case ply_type::i16: property.data_i16.reserve(count); break;
        case ply_type::i32: property.data_i32.reserve(count); break;
        case ply_type::i64: property.data_i64.reserve(count); break;
        case ply_type::u8: property.data_u8.reserve(count); break;
        case ply_type::u16: property.data_u16.reserve(count); break;
        case ply_type::u32: property.data_u32.reserve(count); break;
        case ply_type::u64: property.data_u64.reserve(count); break;
        case ply_type::f32: property.data_f32.reserve(count); break;
        case ply_type::f64: property.data_f64.reserve(count); break;
      }
      if (property.is_list) property.ldata_u8.reserve(element.count);
    }
  }

  // read data -------------------------------------
  if (ply.format == ply_format::ascii) {
    for (auto& elem : ply.elements) {
      for (auto idx = 0; idx < elem.count; idx++) {
        if (!read_line(fs, buffer, sizeof(buffer)))
          throw std::runtime_error("cannot read ply");
        auto line = string_view{buffer};
        for (auto& prop : elem.properties) {
          if (prop.is_list) parse_value(line, prop.ldata_u8.emplace_back());
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                parse_value(line, prop.data_i8.emplace_back());
                break;
              case ply_type::i16:
                parse_value(line, prop.data_i16.emplace_back());
                break;
              case ply_type::i32:
                parse_value(line, prop.data_i32.emplace_back());
                break;
              case ply_type::i64:
                parse_value(line, prop.data_i64.emplace_back());
                break;
              case ply_type::u8:
                parse_value(line, prop.data_u8.emplace_back());
                break;
              case ply_type::u16:
                parse_value(line, prop.data_u16.emplace_back());
                break;
              case ply_type::u32:
                parse_value(line, prop.data_u32.emplace_back());
                break;
              case ply_type::u64:
                parse_value(line, prop.data_u64.emplace_back());
                break;
              case ply_type::f32:
                parse_value(line, prop.data_f32.emplace_back());
                break;
              case ply_type::f64:
                parse_value(line, prop.data_f64.emplace_back());
                break;
            }
          }
        }
      }
    }
  } else {
    auto big_endian = ply.format == ply_format::binary_big_endian;
    for (auto& elem : ply.elements) {
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto& prop : elem.properties) {
          if (prop.is_list)
            read_value(fs, prop.ldata_u8.emplace_back(), big_endian);
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                read_value(fs, prop.data_i8.emplace_back(), big_endian);
                break;
              case ply_type::i16:
                read_value(fs, prop.data_i16.emplace_back(), big_endian);
                break;
              case ply_type::i32:
                read_value(fs, prop.data_i32.emplace_back(), big_endian);
                break;
              case ply_type::i64:
                read_value(fs, prop.data_i64.emplace_back(), big_endian);
                break;
              case ply_type::u8:
                read_value(fs, prop.data_u8.emplace_back(), big_endian);
                break;
              case ply_type::u16:
                read_value(fs, prop.data_u16.emplace_back(), big_endian);
                break;
              case ply_type::u32:
                read_value(fs, prop.data_u32.emplace_back(), big_endian);
                break;
              case ply_type::u64:
                read_value(fs, prop.data_u64.emplace_back(), big_endian);
                break;
              case ply_type::f32:
                read_value(fs, prop.data_f32.emplace_back(), big_endian);
                break;
              case ply_type::f64:
                read_value(fs, prop.data_f64.emplace_back(), big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// Save ply
void save_ply(const string& filename, const ply_model& ply) {
  auto fs = open_file(filename, "wb");

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

  // header
  format_values(fs, "ply\n");
  format_values(fs, "format {} 1.0\n", format_map.at(ply.format));
  format_values(fs, "comment Written by Yocto/GL\n");
  format_values(fs, "comment https://github.com/xelatihy/yocto-gl\n");
  for (auto& comment : ply.comments) format_values(fs, "comment {}\n", comment);
  for (auto& elem : ply.elements) {
    format_values(fs, "element {} {}\n", elem.name, (uint64_t)elem.count);
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        format_values(
            fs, "property list uchar {} {}\n", type_map[prop.type], prop.name);
      } else {
        format_values(fs, "property {} {}\n", type_map[prop.type], prop.name);
      }
    }
  }
  format_values(fs, "end_header\n");

  // properties
  if (ply.format == ply_format::ascii) {
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto pidx = 0; pidx < elem.properties.size(); pidx++) {
          auto& prop = elem.properties[pidx];
          if (prop.is_list) format_values(fs, "{} ", (int)prop.ldata_u8[idx]);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                format_values(fs, "{} ", prop.data_i8[cur[idx]++]);
                break;
              case ply_type::i16:
                format_values(fs, "{} ", prop.data_i16[cur[idx]++]);
                break;
              case ply_type::i32:
                format_values(fs, "{} ", prop.data_i32[cur[idx]++]);
                break;
              case ply_type::i64:
                format_values(
                    fs, "{} ", prop.data_i64[cur[idx]++]);
                break;
              case ply_type::u8:
                format_values(fs, "{} ", prop.data_i8[cur[idx]++]);
                break;
              case ply_type::u16:
                format_values(fs, "{} ", prop.data_i16[cur[idx]++]);
                break;
              case ply_type::u32:
                format_values(fs, "{} ", prop.data_u32[cur[idx]++]);
                break;
              case ply_type::u64:
                format_values(
                    fs, "{} ", prop.data_u64[cur[idx]++]);
                break;
              case ply_type::f32:
                format_values(fs, "{}", prop.data_f32[cur[idx]++]);
                break;
              case ply_type::f64:
                format_values(fs, "{}", prop.data_f64[cur[idx]++]);
                break;
            }
          }
          format_values(fs, "\n");
        }
      }
    }
  } else {
    auto big_endian = ply.format == ply_format::binary_big_endian;
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto pidx = 0; pidx < elem.properties.size(); pidx++) {
          auto& prop = elem.properties[pidx];
          if (prop.is_list) write_value(fs, prop.ldata_u8[idx], big_endian);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                write_value(fs, prop.data_i8[cur[pidx]++], big_endian);
                break;
              case ply_type::i16:
                write_value(fs, prop.data_i16[cur[pidx]++], big_endian);
                break;
              case ply_type::i32:
                write_value(fs, prop.data_i32[cur[pidx]++], big_endian);
                break;
              case ply_type::i64:
                write_value(fs, prop.data_i64[cur[pidx]++], big_endian);
                break;
              case ply_type::u8:
                write_value(fs, prop.data_i8[cur[pidx]++], big_endian);
                break;
              case ply_type::u16:
                write_value(fs, prop.data_i16[cur[pidx]++], big_endian);
                break;
              case ply_type::u32:
                write_value(fs, prop.data_u32[cur[pidx]++], big_endian);
                break;
              case ply_type::u64:
                write_value(fs, prop.data_u64[cur[pidx]++], big_endian);
                break;
              case ply_type::f32:
                write_value(fs, prop.data_f32[cur[pidx]++], big_endian);
                break;
              case ply_type::f64:
                write_value(fs, prop.data_f64[cur[pidx]++], big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// Get ply properties
bool has_ply_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return true;
    }
  }
  return false;
}
const ply_property& get_ply_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
ply_property& get_ply_property(
    ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
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
inline vector<T> convert_ply_property(const ply_property& prop) {
  switch (prop.type) {
    case ply_type::i8: return convert_ply_property<T>(prop.data_i8);
    case ply_type::i16: return convert_ply_property<T>(prop.data_i16);
    case ply_type::i32: return convert_ply_property<T>(prop.data_i32);
    case ply_type::i64: return convert_ply_property<T>(prop.data_i64);
    case ply_type::u8: return convert_ply_property<T>(prop.data_u8);
    case ply_type::u16: return convert_ply_property<T>(prop.data_u16);
    case ply_type::u32: return convert_ply_property<T>(prop.data_u32);
    case ply_type::u64: return convert_ply_property<T>(prop.data_u64);
    case ply_type::f32: return convert_ply_property<T>(prop.data_f32);
    case ply_type::f64: return convert_ply_property<T>(prop.data_f64);
  }
}
vector<float> get_ply_values(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (prop.is_list) return {};
  return convert_ply_property<float>(prop);
}
vector<vec2f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return values;
}
vector<vec3f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto z      = get_ply_values(ply, element, property3);
  auto values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return values;
}
vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    const string& property4) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto z      = get_ply_values(ply, element, property3);
  auto w      = get_ply_values(ply, element, property4);
  auto values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w[i]};
  return values;
}
vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    float property4) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto z      = get_ply_values(ply, element, property3);
  auto w      = property4;
  auto values = vector<vec4f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i], w};
  return values;
}
vector<vector<int>> get_ply_lists(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (!prop.is_list) return {};
  auto& sizes  = prop.ldata_u8;
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
vector<byte> get_ply_list_sizes(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (!prop.is_list) return {};
  return prop.ldata_u8;
}
vector<int> get_ply_list_values(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (!prop.is_list) return {};
  return convert_ply_property<int>(prop);
}

// Get ply properties for meshes
vector<vec3f> get_ply_positions(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "x", "y", "z");
}
vector<vec3f> get_ply_normals(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "nx", "ny", "nz");
}
vector<vec2f> get_ply_texcoords(const ply_model& ply, bool flipv) {
  auto texcoord = has_ply_property(ply, "vertex", "u")
                      ? get_ply_values(ply, "vertex", "u", "v")
                      : get_ply_values(ply, "vertex", "s", "t");
  return flipv ? flip_texcoord(texcoord) : texcoord;
}
vector<vec4f> get_ply_colors(const ply_model& ply) {
  if (has_ply_property(ply, "vertex", "alpha")) {
    return get_ply_values(ply, "vertex", "red", "green", "blue", "alpha");
  } else {
    return get_ply_values(ply, "vertex", "red", "green", "blue", 1);
  }
}
vector<float> get_ply_radius(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "radius");
}
vector<vector<int>> get_ply_faces(const ply_model& ply) {
  return get_ply_lists(ply, "face", "vertex_indices");
}
vector<vec3i> get_ply_triangles(const ply_model& ply) {
  auto indices   = get_ply_list_values(ply, "face", "vertex_indices");
  auto sizes     = get_ply_list_sizes(ply, "face", "vertex_indices");
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
vector<vec4i> get_ply_quads(const ply_model& ply) {
  auto indices = get_ply_list_values(ply, "face", "vertex_indices");
  auto sizes   = get_ply_list_sizes(ply, "face", "vertex_indices");
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
vector<vec2i> get_ply_lines(const ply_model& ply) {
  auto indices = get_ply_list_values(ply, "line", "vertex_indices");
  auto sizes   = get_ply_list_sizes(ply, "line", "vertex_indices");
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
vector<int> get_ply_points(const ply_model& ply) {
  return get_ply_list_values(ply, "point", "vertex_indices");
}
bool has_ply_quads(const ply_model& ply) {
  auto sizes = get_ply_list_sizes(ply, "face", "vertex_indices");
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
void add_ply_element(ply_model& ply, const string& element, size_t count) {
  for (auto& elem : ply.elements) {
    if (elem.name == element) return;
  }
  auto& elem = ply.elements.emplace_back();
  elem.name  = element;
  elem.count = count;
}
void add_ply_property(ply_model& ply, const string& element,
    const string& property, size_t count, ply_type type, bool is_list) {
  add_ply_element(ply, element, count);
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property)
        throw std::runtime_error("property already added");
    }
    auto& prop   = elem.properties.emplace_back();
    prop.name    = property;
    prop.type    = type;
    prop.is_list = is_list;
    return;
  }
}
template <typename T>
vector<T> make_ply_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

void add_ply_values(ply_model& ply, const float* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (!values) return;
  for (auto p = 0; p < nprops; p++) {
    add_ply_property(ply, element, properties[p], count, ply_type::f32, false);
    auto& prop = get_ply_property(ply, element, properties[p]);
    prop.data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop.data_f32[i] = values[p + i * nprops];
  }
}

void add_ply_values(ply_model& ply, const vector<float>& values,
    const string& element, const string& property) {
  auto properties = vector{property};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
void add_ply_values(ply_model& ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2) {
  auto properties = vector{property1, property2};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
void add_ply_values(ply_model& ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3) {
  auto properties = vector{property1, property2, property3};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
void add_ply_values(ply_model& ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4) {
  auto properties = vector{property1, property2, property3, property4};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}

void add_ply_lists(ply_model& ply, const vector<vector<int>>& values,
    const string& element, const string& property) {
  if (values.empty()) return;
  add_ply_property(ply, element, property, values.size(), ply_type::i32, true);
  auto& prop = get_ply_property(ply, element, property);
  prop.data_i32.reserve(values.size() * 4);
  prop.ldata_u8.reserve(values.size());
  for (auto& value : values) {
    prop.data_i32.insert(prop.data_i32.end(), value.begin(), value.end());
    prop.ldata_u8.push_back((uint8_t)value.size());
  }
}
void add_ply_lists(ply_model& ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property) {
  if (values.empty()) return;
  add_ply_property(ply, element, property, sizes.size(), ply_type::i32, true);
  auto& prop    = get_ply_property(ply, element, property);
  prop.data_i32 = values;
  prop.ldata_u8 = sizes;
}
void add_ply_lists(ply_model& ply, const int* values, size_t count, int size,
    const string& element, const string& property) {
  if (!values) return;
  add_ply_property(ply, element, property, count, ply_type::i32, true);
  auto& prop = get_ply_property(ply, element, property);
  prop.data_i32.assign(values, values + count * size);
  prop.ldata_u8.assign(count, size);
}
void add_ply_lists(ply_model& ply, const vector<int>& values,
    const string& element, const string& property) {
  return add_ply_lists(ply, values.data(), values.size(), 1, element, property);
}
void add_ply_lists(ply_model& ply, const vector<vec2i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
void add_ply_lists(ply_model& ply, const vector<vec3i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
void add_ply_lists(ply_model& ply, const vector<vec4i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
}

// Add ply properties for meshes
void add_ply_positions(ply_model& ply, const vector<vec3f>& values) {
  return add_ply_values(ply, values, "vertex", "x", "y", "z");
}
void add_ply_normals(ply_model& ply, const vector<vec3f>& values) {
  return add_ply_values(ply, values, "vertex", "nx", "ny", "nz");
}
void add_ply_texcoords(
    ply_model& ply, const vector<vec2f>& values, bool flipv) {
  return add_ply_values(
      ply, flipv ? flip_texcoord(values) : values, "vertex", "u", "v");
}
void add_ply_colors(ply_model& ply, const vector<vec4f>& values) {
  return add_ply_values(ply, values, "vertex", "red", "green", "blue", "alpha");
}
void add_ply_radius(ply_model& ply, const vector<float>& values) {
  return add_ply_values(ply, values, "vertex", "radius");
}
void add_ply_faces(ply_model& ply, const vector<vector<int>>& values) {
  return add_ply_lists(ply, values, "face", "vertex_indices");
}
void add_ply_faces(ply_model& ply, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
  if (triangles.empty() && quads.empty()) return;
  if (quads.empty()) {
    return add_ply_lists(ply, triangles, "face", "vertex_indices");
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const vec4i& q) { return q.z != q.w; })) {
    return add_ply_lists(ply, quads, "face", "vertex_indices");
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
    return add_ply_lists(ply, sizes, indices, "face", "vertex_indices");
  }
}
void add_ply_triangles(ply_model& ply, const vector<vec3i>& values) {
  return add_ply_faces(ply, values, {});
}
void add_ply_quads(ply_model& ply, const vector<vec4i>& values) {
  return add_ply_faces(ply, {}, values);
}
void add_ply_lines(ply_model& ply, const vector<vec2i>& values) {
  return add_ply_lists(ply, values, "line", "vertex_indices");
}
void add_ply_points(ply_model& ply, const vector<int>& values) {
  return add_ply_lists(ply, values, "point", "vertex_indices");
}

// get ply value either ascii or binary
template <typename T, typename VT>
static inline void read_ply_prop(file_wrapper& fs, VT& value, bool big_endian) {
  auto tvalue = T{};
  read_value(fs, tvalue, big_endian);
  value = (VT)tvalue;
}
template <typename VT>
static inline void read_ply_prop(
    file_wrapper& fs, ply_type type, VT& value, bool big_endian) {
  switch (type) {
    case ply_type::i8: read_ply_prop<int8_t>(fs, value, big_endian); break;
    case ply_type::i16: read_ply_prop<int16_t>(fs, value, big_endian); break;
    case ply_type::i32: read_ply_prop<int32_t>(fs, value, big_endian); break;
    case ply_type::i64: read_ply_prop<int64_t>(fs, value, big_endian); break;
    case ply_type::u8: read_ply_prop<uint8_t>(fs, value, big_endian); break;
    case ply_type::u16: read_ply_prop<uint16_t>(fs, value, big_endian); break;
    case ply_type::u32: read_ply_prop<uint32_t>(fs, value, big_endian); break;
    case ply_type::u64: read_ply_prop<uint64_t>(fs, value, big_endian); break;
    case ply_type::f32: read_ply_prop<float>(fs, value, big_endian); break;
    case ply_type::f64: read_ply_prop<double>(fs, value, big_endian); break;
  }
}

template <typename T, typename VT>
static inline void parse_ply_prop(string_view& str, VT& value) {
  auto tvalue = T{};
  parse_value(str, tvalue);
  value = (VT)tvalue;
}
template <typename VT>
static inline void parse_ply_prop(string_view& str, ply_type type, VT& value) {
  switch (type) {
    case ply_type::i8: parse_ply_prop<int8_t>(str, value); break;
    case ply_type::i16: parse_ply_prop<int16_t>(str, value); break;
    case ply_type::i32: parse_ply_prop<int32_t>(str, value); break;
    case ply_type::i64: parse_ply_prop<int64_t>(str, value); break;
    case ply_type::u8: parse_ply_prop<uint8_t>(str, value); break;
    case ply_type::u16: parse_ply_prop<uint16_t>(str, value); break;
    case ply_type::u32: parse_ply_prop<uint32_t>(str, value); break;
    case ply_type::u64: parse_ply_prop<uint64_t>(str, value); break;
    case ply_type::f32: parse_ply_prop<float>(str, value); break;
    case ply_type::f64: parse_ply_prop<double>(str, value); break;
  }
}

// Load ply data
void read_ply_header(file_wrapper& fs, ply_format& format,
    vector<ply_element>& elements, vector<string>& comments) {
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

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // prepare elements
  elements.clear();

  // read the file header line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_ply_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // check magic number
    if (first_line) {
      if (cmd != "ply") throw std::runtime_error{"bad ply file"};
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) throw std::runtime_error{"bad ply file"};
    } else if (cmd == "format") {
      auto fmt = ""sv;
      parse_value(line, fmt);
      if (fmt == "ascii") {
        format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        format = ply_format::binary_big_endian;
      } else {
        throw std::runtime_error{"unknown ply format"};
      }
    } else if (cmd == "comment") {
      skip_whitespace(line);
      comments.push_back(string{line});
    } else if (cmd == "obj_info") {
      skip_whitespace(line);
      // comment is the rest of the line
    } else if (cmd == "element") {
      auto& elem = elements.emplace_back();
      parse_value(line, elem.name);
      parse_value(line, elem.count);
    } else if (cmd == "property") {
      if (elements.empty()) throw std::runtime_error{"bad ply header"};
      auto& prop  = elements.back().properties.emplace_back();
      auto  tname = ""s;
      parse_value(line, tname);
      if (tname == "list") {
        prop.is_list = true;
        parse_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          throw std::runtime_error{"unsupported list size type " + tname};
        parse_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      }
      parse_value(line, prop.name);
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      throw std::runtime_error{"unknown ply command"};
    }
  }

  if (!end_header) throw std::runtime_error{"bad ply header"};
}

template <typename VT, typename LT>
void read_ply_value_generic(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<VT>& values, vector<vector<LT>>& lists) {
  // prepare properties
  if (values.size() != element.properties.size()) {
    values.resize(element.properties.size());
  }
  if (lists.size() != element.properties.size()) {
    lists.resize(element.properties.size());
  }
  for (auto& list : lists) list.clear();

  // read property values
  if (format == ply_format::ascii) {
    char buffer[4096];
    if (!read_line(fs, buffer, sizeof(buffer)))
      throw std::runtime_error("cannot read ply");
    auto line = string_view{buffer};
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop  = element.properties[pidx];
      auto& value = values[pidx];
      auto& list  = lists[pidx];
      if (!prop.is_list) {
        parse_ply_prop(line, prop.type, value);
      } else {
        parse_ply_prop(line, ply_type::u8, value);
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          parse_ply_prop(line, prop.type, list[i]);
      }
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop  = element.properties[pidx];
      auto& value = values[pidx];
      auto& list  = lists[pidx];
      if (!prop.is_list) {
        read_ply_prop(
            fs, prop.type, value, format == ply_format::binary_big_endian);
      } else {
        read_ply_prop(
            fs, ply_type::u8, value, format == ply_format::binary_big_endian);
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          read_ply_prop(
              fs, prop.type, list[i], format == ply_format::binary_big_endian);
      }
    }
  }
}

template <typename VT>
static inline void format_ply_prop(file_wrapper& fs, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: format_value(fs, (int8_t)value); break;
    case ply_type::i16: format_value(fs, (int16_t)value); break;
    case ply_type::i32: format_value(fs, (int32_t)value); break;
    case ply_type::i64: format_value(fs, (int64_t)value); break;
    case ply_type::u8: format_value(fs, (uint8_t)value); break;
    case ply_type::u16: format_value(fs, (uint16_t)value); break;
    case ply_type::u32: format_value(fs, (uint32_t)value); break;
    case ply_type::u64: format_value(fs, (uint64_t)value); break;
    case ply_type::f32: format_value(fs, (float)value); break;
    case ply_type::f64: format_value(fs, (double)value); break;
  }
}

template <typename VT>
static inline void write_ply_prop(
    file_wrapper& fs, ply_type type, VT value, bool big_endian) {
  switch (type) {
    case ply_type::i8: write_value(fs, (int8_t)value, big_endian); break;
    case ply_type::i16: write_value(fs, (int16_t)value, big_endian); break;
    case ply_type::i32: write_value(fs, (int32_t)value, big_endian); break;
    case ply_type::i64: write_value(fs, (int64_t)value, big_endian); break;
    case ply_type::u8: write_value(fs, (uint8_t)value, big_endian); break;
    case ply_type::u16: write_value(fs, (uint16_t)value, big_endian); break;
    case ply_type::u32: write_value(fs, (uint32_t)value, big_endian); break;
    case ply_type::u64: write_value(fs, (uint64_t)value, big_endian); break;
    case ply_type::f32: write_value(fs, (float)value, big_endian); break;
    case ply_type::f64: write_value(fs, (double)value, big_endian); break;
  }
}

// Write Ply functions
void write_ply_header(file_wrapper& fs, ply_format format,
    const vector<ply_element>& elements, const vector<string>& comments) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{{ply_type::i8, "char"},
      {ply_type::i16, "short"}, {ply_type::i32, "int"}, {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"}, {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"}, {ply_type::u64, "ulong"},
      {ply_type::f32, "float"}, {ply_type::f64, "double"}};

  write_text(fs, "ply\n");
  switch (format) {
    case ply_format::ascii: write_text(fs, "format ascii 1.0\n"); break;
    case ply_format::binary_little_endian:
      write_text(fs, "format binary_little_endian 1.0\n");
      break;
    case ply_format::binary_big_endian:
      write_text(fs, "format binary_big_endian 1.0\n");
      break;
  }
  for (auto& comment : comments) write_text(fs, "comment " + comment + "\n");
  for (auto& elem : elements) {
    write_text(
        fs, "element " + elem.name + " " + std::to_string(elem.count) + "\n");
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        write_text(fs, "property list uchar " + type_map[prop.type] + " " +
                           prop.name + "\n");
      } else {
        write_text(
            fs, "property " + type_map[prop.type] + " " + prop.name + "\n");
      }
    }
  }
  write_text(fs, "end_header\n");
}

template <typename VT, typename LT>
void write_ply_value_generic(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<VT>& values, vector<vector<LT>>& lists) {
  if (format == ply_format::ascii) {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (pidx) format_value(fs, " ");
      if (!prop.is_list) {
        format_ply_prop(fs, prop.type, values[pidx]);
      } else {
        format_ply_prop(fs, ply_type::u8, values[pidx]);
        for (auto i = 0; i < (int)lists[pidx].size(); i++) {
          if (i) format_value(fs, " ");
          format_ply_prop(fs, prop.type, lists[pidx][i]);
        }
      }
      format_value(fs, "\n");
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (!prop.is_list) {
        write_ply_prop(fs, prop.type, values[pidx],
            format == ply_format::binary_big_endian);
      } else {
        write_ply_prop(fs, ply_type::u8, values[pidx],
            format == ply_format::binary_big_endian);
        for (auto i = 0; i < (int)lists[pidx].size(); i++)
          write_ply_prop(fs, prop.type, lists[pidx][i],
              format == ply_format::binary_big_endian);
      }
    }
  }
}

void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  write_ply_value_generic(fs, format, element, values, lists);
}
void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  write_ply_value_generic(fs, format, element, values, lists);
}

void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  read_ply_value_generic(fs, format, element, values, lists);
}
void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  read_ply_value_generic(fs, format, element, values, lists);
}

int find_ply_element(const vector<ply_element>& elements, const string& name) {
  for (auto idx = 0; idx < elements.size(); idx++)
    if (elements[idx].name == name) return idx;
  return -1;
}
int find_ply_property(const ply_element& element, const string& name) {
  for (auto idx = 0; idx < element.properties.size(); idx++)
    if (element.properties[idx].name == name) return idx;
  return -1;
}
vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2) {
  auto ids = vec2i{
      find_ply_property(element, name1),
      find_ply_property(element, name2),
  };
  if (ids.x < 0 || ids.y < 0) return vec2i{-1};
  return ids;
}
vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3) {
  auto ids = vec3i{
      find_ply_property(element, name1),
      find_ply_property(element, name2),
      find_ply_property(element, name3),
  };
  if (ids.x < 0 || ids.y < 0 || ids.z < 0) return vec3i{-1};
  return ids;
}
vec4i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3, const string& name4) {
  auto ids = vec4i{
      find_ply_property(element, name1),
      find_ply_property(element, name2),
      find_ply_property(element, name3),
      find_ply_property(element, name4),
  };
  if (ids.x < 0 || ids.y < 0 || ids.z < 0 || ids.w < 0) return vec4i{-1};
  return ids;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static inline void remove_obj_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

static inline void parse_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  parse_value(str, value.position);
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      parse_value(str, value.normal);
    } else {
      parse_value(str, value.texcoord);
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        parse_value(str, value.normal);
      }
    }
  }
}

// Input for OBJ textures
static inline void parse_value(string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    parse_value(str, token);
    tokens.push_back(token);
    skip_whitespace(str);
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
void load_mtl(const string& filename, obj_model& obj, bool fliptr = true) {
  // open file
  auto fs = open_file(filename, "rt");

  // init parsing
  obj.materials.emplace_back();

  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      obj.materials.emplace_back();
      parse_value(line, obj.materials.back().name);
    } else if (cmd == "illum") {
      parse_value(line, obj.materials.back().illum);
    } else if (cmd == "Ke") {
      parse_value(line, obj.materials.back().emission);
    } else if (cmd == "Ka") {
      parse_value(line, obj.materials.back().ambient);
    } else if (cmd == "Kd") {
      parse_value(line, obj.materials.back().diffuse);
    } else if (cmd == "Ks") {
      parse_value(line, obj.materials.back().specular);
    } else if (cmd == "Kt") {
      parse_value(line, obj.materials.back().transmission);
    } else if (cmd == "Tf") {
      obj.materials.back().transmission = vec3f{-1};
      parse_value(line, obj.materials.back().transmission);
      if (obj.materials.back().transmission.y < 0)
        obj.materials.back().transmission = vec3f{
            obj.materials.back().transmission.x};
      if (fliptr)
        obj.materials.back().transmission = 1 -
                                            obj.materials.back().transmission;
    } else if (cmd == "Tr") {
      parse_value(line, obj.materials.back().opacity);
      if (fliptr)
        obj.materials.back().opacity = 1 - obj.materials.back().opacity;
    } else if (cmd == "Ns") {
      parse_value(line, obj.materials.back().exponent);
    } else if (cmd == "d") {
      parse_value(line, obj.materials.back().opacity);
    } else if (cmd == "map_Ke") {
      parse_value(line, obj.materials.back().emission_map);
    } else if (cmd == "map_Ka") {
      parse_value(line, obj.materials.back().ambient_map);
    } else if (cmd == "map_Kd") {
      parse_value(line, obj.materials.back().diffuse_map);
    } else if (cmd == "map_Ks") {
      parse_value(line, obj.materials.back().specular_map);
    } else if (cmd == "map_Tr") {
      parse_value(line, obj.materials.back().transmission_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_value(line, obj.materials.back().opacity_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_value(line, obj.materials.back().bump_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_value(line, obj.materials.back().displacement_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_value(line, obj.materials.back().normal_map);
    } else if (cmd == "Pm") {
      parse_value(line, obj.materials.back().pbr_metallic);
    } else if (cmd == "Pr") {
      parse_value(line, obj.materials.back().pbr_roughness);
    } else if (cmd == "Ps") {
      parse_value(line, obj.materials.back().pbr_sheen);
    } else if (cmd == "Pc") {
      parse_value(line, obj.materials.back().pbr_clearcoat);
    } else if (cmd == "Pcr") {
      parse_value(line, obj.materials.back().pbr_coatroughness);
    } else if (cmd == "map_Pm") {
      parse_value(line, obj.materials.back().pbr_metallic_map);
    } else if (cmd == "map_Pr") {
      parse_value(line, obj.materials.back().pbr_roughness_map);
    } else if (cmd == "map_Ps") {
      parse_value(line, obj.materials.back().pbr_sheen_map);
    } else if (cmd == "map_Pc") {
      parse_value(line, obj.materials.back().pbr_clearcoat_map);
    } else if (cmd == "map_Pcr") {
      parse_value(line, obj.materials.back().pbr_coatroughness_map);
    } else if (cmd == "Vt") {
      parse_value(line, obj.materials.back().vol_transmission);
    } else if (cmd == "Vp") {
      parse_value(line, obj.materials.back().vol_meanfreepath);
    } else if (cmd == "Ve") {
      parse_value(line, obj.materials.back().vol_emission);
    } else if (cmd == "Vs") {
      parse_value(line, obj.materials.back().vol_scattering);
    } else if (cmd == "Vg") {
      parse_value(line, obj.materials.back().vol_anisotropy);
    } else if (cmd == "Vr") {
      parse_value(line, obj.materials.back().vol_scale);
    } else if (cmd == "map_Vs") {
      parse_value(line, obj.materials.back().vol_scattering_map);
    } else {
      continue;
    }
  }

  // remove placeholder material
  obj.materials.erase(obj.materials.begin());
}

// Read obj
void load_objx(const string& filename, obj_model& obj) {
  // open file
  auto fs = open_file(filename, "rt");

  // shape map for instances
  auto shape_map = unordered_map<string, vector<int>>{};
  for (auto idx = 0; idx < obj.shapes.size(); idx++) {
    shape_map[obj.shapes[idx].name].push_back(idx);
  }

  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto& camera = obj.cameras.emplace_back();
      parse_value(line, camera.name);
      parse_value(line, camera.ortho);
      parse_value(line, camera.width);
      parse_value(line, camera.height);
      parse_value(line, camera.lens);
      parse_value(line, camera.focus);
      parse_value(line, camera.aperture);
      parse_value(line, camera.frame);
    } else if (cmd == "e") {
      auto& environment = obj.environments.emplace_back();
      parse_value(line, environment.name);
      parse_value(line, environment.emission);
      auto emission_path = ""s;
      parse_value(line, emission_path);
      if (emission_path == "\"\"") emission_path = "";
      environment.emission_map.path = emission_path;
      parse_value(line, environment.frame);
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      parse_value(line, object);
      parse_value(line, frame);
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
void load_obj(const string& filename, obj_model& obj, bool geom_only,
    bool split_elements, bool split_materials) {
  // open file
  auto fs = open_file(filename, "rt");

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
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      parse_value(line, opositions.emplace_back());
      vert_size.position += 1;
    } else if (cmd == "vn") {
      parse_value(line, onormals.emplace_back());
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      parse_value(line, otexcoords.emplace_back());
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
      skip_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        shape.vertices.push_back(vert);
        element.size += 1;
        skip_whitespace(line);
      }
    } else if (cmd == "o" || cmd == "g") {
      if (geom_only) continue;
      parse_value_or_empty(line, cmd == "o" ? oname : gname);
      if (!obj.shapes.back().vertices.empty()) {
        obj.shapes.emplace_back();
        obj.shapes.back().name = oname + gname;
      } else {
        obj.shapes.back().name = oname + gname;
      }
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      parse_value_or_empty(line, mname);
    } else if (cmd == "s") {
      if (geom_only) continue;
      // TODO: smoothing
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      parse_value(line, mtllib);
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
static inline void format_value(string& str, const obj_texture_info& value) {
  str += value.path.empty() ? "" : value.path;
}
static inline void format_value(string& str, const obj_vertex& value) {
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
void save_mtl(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // write material
  for (auto& material : obj.materials) {
    format_values(fs, "newmtl {}\n", material.name);
    format_values(fs, "illum {}\n", material.illum);
    if (material.emission != zero3f)
      format_values(fs, "Ke {}\n", material.emission);
    if (material.ambient != zero3f)
      format_values(fs, "Ka {}\n", material.ambient);
    format_values(fs, "Kd {}\n", material.diffuse);
    format_values(fs, "Ks {}\n", material.specular);
    if (material.reflection != zero3f)
      format_values(fs, "Kr {}\n", material.reflection);
    if (material.transmission != zero3f)
      format_values(fs, "Kt {}\n", material.transmission);
    format_values(fs, "Ns {}\n", (int)material.exponent);
    if (material.opacity != 1) format_values(fs, "d {}\n", material.opacity);
    if (!material.emission_map.path.empty())
      format_values(fs, "map_Ke {}\n", material.emission_map);
    if (!material.diffuse_map.path.empty())
      format_values(fs, "map_Kd {}\n", material.diffuse_map);
    if (!material.specular_map.path.empty())
      format_values(fs, "map_Ks {}\n", material.specular_map);
    if (!material.transmission_map.path.empty())
      format_values(fs, "map_Kt {}\n", material.transmission_map);
    if (!material.reflection_map.path.empty())
      format_values(fs, "map_Kr {}\n", material.reflection_map);
    if (!material.exponent_map.path.empty())
      format_values(fs, "map_Ns {}\n", material.exponent_map);
    if (!material.opacity_map.path.empty())
      format_values(fs, "map_d {}\n", material.opacity_map);
    if (!material.bump_map.path.empty())
      format_values(fs, "map_bump {}\n", material.bump_map);
    if (!material.displacement_map.path.empty())
      format_values(fs, "map_disp {}\n", material.displacement_map);
    if (!material.normal_map.path.empty())
      format_values(fs, "map_norm {}\n", material.normal_map);
    if (material.pbr_roughness)
      format_values(fs, "Pr {}\n", material.pbr_roughness);
    if (material.pbr_metallic)
      format_values(fs, "Pm {}\n", material.pbr_metallic);
    if (material.pbr_sheen) format_values(fs, "Ps {}\n", material.pbr_sheen);
    if (material.pbr_clearcoat)
      format_values(fs, "Pc {}\n", material.pbr_clearcoat);
    if (material.pbr_coatroughness)
      format_values(fs, "Pcr {}\n", material.pbr_coatroughness);
    if (!material.pbr_roughness_map.path.empty())
      format_values(fs, "map_Pr {}\n", material.pbr_roughness_map);
    if (!material.pbr_metallic_map.path.empty())
      format_values(fs, "map_Pm {}\n", material.pbr_metallic_map);
    if (!material.pbr_sheen_map.path.empty())
      format_values(fs, "map_Ps {}\n", material.pbr_sheen_map);
    if (!material.pbr_clearcoat_map.path.empty())
      format_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map);
    if (!material.pbr_coatroughness_map.path.empty())
      format_values(fs, "map_Pcr {}\n", material.pbr_coatroughness_map);
    if (material.vol_transmission != zero3f)
      format_values(fs, "Vt {}\n", material.vol_transmission);
    if (material.vol_meanfreepath != zero3f)
      format_values(fs, "Vp {}\n", material.vol_meanfreepath);
    if (material.vol_emission != zero3f)
      format_values(fs, "Ve {}\n", material.vol_emission);
    if (material.vol_scattering != zero3f)
      format_values(fs, "Vs {}\n", material.vol_scattering);
    if (material.vol_anisotropy)
      format_values(fs, "Vg {}\n", material.vol_anisotropy);
    if (material.vol_scale) format_values(fs, "Vr {}\n", material.vol_scale);
    if (!material.vol_scattering_map.path.empty())
      format_values(fs, "map_Vs {}\n", material.vol_scattering_map);
    format_values(fs, "\n");
  }
}

// Save obj
void save_objx(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // cameras
  for (auto& camera : obj.cameras) {
    format_values(fs, "c {} {} {} {} {} {} {} {}\n", camera.name, camera.ortho,
        camera.width, camera.height, camera.lens, camera.focus, camera.aperture,
        camera.frame);
  }

  // environments
  for (auto& environment : obj.environments) {
    format_values(fs, "e {} {} {} {}\n", environment.name, environment.emission,
        environment.emission_map.path.empty() ? "\"\""s
                                              : environment.emission_map.path,
        environment.frame);
  }

  // instances
  for (auto& shape : obj.shapes) {
    for (auto& frame : shape.instances) {
      format_values(fs, "i {} {}\n", shape.name, frame);
    }
  }
}

// Save obj
void save_obj(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : obj.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  // save material library
  if (!obj.materials.empty()) {
    format_values(
        fs, "mtllib {}\n\n", replace_extension(get_filename(filename), ".mtl"));
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto& shape : obj.shapes) {
    format_values(fs, "o {}\n", shape.name);
    for (auto& p : shape.positions) format_values(fs, "v {}\n", p);
    for (auto& n : shape.normals) format_values(fs, "vn {}\n", n);
    for (auto& t : shape.texcoords) format_values(fs, "vt {}\n", t);
    auto element_labels = vector<string>{"f", "l", "p"};
    auto element_groups = vector<const vector<obj_element>*>{
        &shape.faces, &shape.lines, &shape.points};
    for (auto element_idx = 0; element_idx < 3; element_idx++) {
      auto& label        = element_labels[element_idx];
      auto& elements     = *element_groups[element_idx];
      auto  cur_material = -1, cur_vertex = 0;
      for (auto& element : elements) {
        if (!shape.materials.empty() && cur_material != element.material) {
          format_values(fs, "usemtl {}\n", shape.materials[element.material]);
          cur_material = element.material;
        }
        format_values(fs, "{}", label);
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape.vertices[cur_vertex++];
          if (vert.position) vert.position += vert_size.position;
          if (vert.normal) vert.normal += vert_size.normal;
          if (vert.texcoord) vert.texcoord += vert_size.texcoord;
          format_values(fs, " {}", vert);
        }
        format_values(fs, "\n");
      }
    }
    format_values(fs, "\n");
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
float obj_exponent_to_roughness(float exponent) {
  auto roughness = exponent;
  roughness      = pow(2 / (roughness + 2), 1 / 4.0f);
  if (roughness < 0.01f) roughness = 0;
  if (roughness > 0.99f) roughness = 1;
  return roughness;
}
float obj_roughness_to_exponent(float roughness) {
  return (int)clamp(
      2 / pow(clamp(roughness, 0.0f, 0.99f) + 1e-10f, 4.0f) - 2, 0.0f, 1.0e9f);
}

// Get obj vertices
void get_obj_vertices(const obj_shape& shape, vector<vec3f>& positions,
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

// Get obj shape
void get_obj_triangles(const obj_model& obj, const obj_shape& shape,
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
void get_obj_quads(const obj_model& obj, const obj_shape& shape,
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
void get_obj_lines(const obj_model& obj, const obj_shape& shape,
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
void get_obj_points(const obj_model& obj, const obj_shape& shape,
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
void get_obj_fvquads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flipv) {
  if (shape.faces.empty()) return;
  positions = shape.positions;
  normals   = shape.normals;
  texcoords = flipv ? flip_texcoord(shape.texcoords) : shape.texcoords;
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

bool has_obj_quads(const obj_shape& shape) {
  for (auto& face : shape.faces)
    if (face.size == 4) return true;
  return false;
}

// Add obj shape
void add_obj_triangles(obj_model& obj, obj_shape& shape,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<int>& ematerials, bool flipv) {
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
void add_obj_quads(obj_model& obj, obj_shape& shape, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
void add_obj_lines(obj_model& obj, obj_shape& shape, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
void add_obj_points(obj_model& obj, obj_shape& shape, const vector<int>& points,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials, bool flipv) {
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
void add_obj_fvquads(obj_model& obj, obj_shape& shape,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<int>& ematerials, bool flipv) {
  shape.positions = positions;
  shape.normals   = normals;
  shape.texcoords = flipv ? flip_texcoord(texcoords) : texcoords;
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
bool read_obj_command(file_wrapper& fs, obj_command& command, string& name,
    vec3f& value, vector<obj_vertex>& vertices, obj_vertex& vert_size) {
  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      parse_value(line, value);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      parse_value(line, value);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
      parse_value(line, (vec2f&)value);
      value.z = 0;
      vert_size.texcoord += 1;
      return true;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      skip_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_value(line, vert);
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        vertices.push_back(vert);
        skip_whitespace(line);
      }
      if (cmd == "f") command = obj_command::face;
      if (cmd == "l") command = obj_command::line;
      if (cmd == "p") command = obj_command::point;
      return true;
    } else if (cmd == "o") {
      command = obj_command::object;
      parse_value_or_empty(line, name);
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      parse_value_or_empty(line, name);
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      parse_value_or_empty(line, name);
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      parse_value_or_empty(line, name);
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      parse_value(line, name);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
bool read_mtl_command(file_wrapper& fs, mtl_command& command,
    obj_material& material, bool fliptr) {
  material = {};

  // read the file line by line
  auto pos   = ftell(fs.fs);
  auto found = false;
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
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
      parse_value(line, material.name);
    } else if (cmd == "illum") {
      parse_value(line, material.illum);
    } else if (cmd == "Ke") {
      parse_value(line, material.emission);
    } else if (cmd == "Ka") {
      parse_value(line, material.ambient);
    } else if (cmd == "Kd") {
      parse_value(line, material.diffuse);
    } else if (cmd == "Ks") {
      parse_value(line, material.specular);
    } else if (cmd == "Kt") {
      parse_value(line, material.transmission);
    } else if (cmd == "Tf") {
      material.transmission = vec3f{-1};
      parse_value(line, material.transmission);
      if (material.transmission.y < 0)
        material.transmission = vec3f{material.transmission.x};
      if (fliptr) material.transmission = 1 - material.transmission;
    } else if (cmd == "Tr") {
      parse_value(line, material.opacity);
      if (fliptr) material.opacity = 1 - material.opacity;
    } else if (cmd == "Ns") {
      parse_value(line, material.exponent);
    } else if (cmd == "d") {
      parse_value(line, material.opacity);
    } else if (cmd == "map_Ke") {
      parse_value(line, material.emission_map);
    } else if (cmd == "map_Ka") {
      parse_value(line, material.ambient_map);
    } else if (cmd == "map_Kd") {
      parse_value(line, material.diffuse_map);
    } else if (cmd == "map_Ks") {
      parse_value(line, material.specular_map);
    } else if (cmd == "map_Tr") {
      parse_value(line, material.transmission_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_value(line, material.opacity_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_value(line, material.bump_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_value(line, material.displacement_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_value(line, material.normal_map);
    } else if (cmd == "Pm") {
      parse_value(line, material.pbr_metallic);
    } else if (cmd == "Pr") {
      parse_value(line, material.pbr_roughness);
    } else if (cmd == "Ps") {
      parse_value(line, material.pbr_sheen);
    } else if (cmd == "Pc") {
      parse_value(line, material.pbr_clearcoat);
    } else if (cmd == "Pcr") {
      parse_value(line, material.pbr_coatroughness);
    } else if (cmd == "map_Pm") {
      parse_value(line, material.pbr_metallic_map);
    } else if (cmd == "map_Pr") {
      parse_value(line, material.pbr_roughness_map);
    } else if (cmd == "map_Ps") {
      parse_value(line, material.pbr_sheen_map);
    } else if (cmd == "map_Pc") {
      parse_value(line, material.pbr_clearcoat_map);
    } else if (cmd == "map_Pcr") {
      parse_value(line, material.pbr_coatroughness_map);
    } else if (cmd == "Vt") {
      parse_value(line, material.vol_transmission);
    } else if (cmd == "Vp") {
      parse_value(line, material.vol_meanfreepath);
    } else if (cmd == "Ve") {
      parse_value(line, material.vol_emission);
    } else if (cmd == "Vs") {
      parse_value(line, material.vol_scattering);
    } else if (cmd == "Vg") {
      parse_value(line, material.vol_anisotropy);
    } else if (cmd == "Vr") {
      parse_value(line, material.vol_scale);
    } else if (cmd == "map_Vs") {
      parse_value(line, material.vol_scattering_map);
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
bool read_objx_command(file_wrapper& fs, objx_command& command,
    obj_camera& camera, obj_environment& environment, obj_instance& instance) {
  // read the file line by line
  char buffer[4096];
  auto found = false;
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_value(line, cmd);
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      command = objx_command::camera;
      parse_value(line, camera.name);
      parse_value(line, camera.ortho);
      parse_value(line, camera.width);
      parse_value(line, camera.height);
      parse_value(line, camera.lens);
      parse_value(line, camera.focus);
      parse_value(line, camera.aperture);
      parse_value(line, camera.frame);
      return true;
    } else if (cmd == "e") {
      command = objx_command::environment;
      parse_value(line, environment.name);
      parse_value(line, environment.emission);
      parse_value(line, environment.emission_map);
      parse_value(line, environment.frame);
      return true;
    } else if (cmd == "i") {
      command = objx_command::instance;
      parse_value(line, instance.object);
      parse_value(line, instance.frame);
      return true;
    }
  }

  if (found) return true;

  return false;
}

// Write obj elements
void write_obj_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    format_values(fs, "# {}\n", line);
  }
  format_values(fs, "\n");
}

void write_obj_command(file_wrapper& fs, obj_command command,
    const string& name, const vec3f& value,
    const vector<obj_vertex>& vertices) {
  switch (command) {
    case obj_command::vertex: format_values(fs, "v {}\n", value); break;
    case obj_command::normal: format_values(fs, "vn {}\n", value); break;
    case obj_command::texcoord: format_values(fs, "vt {}\n", value); break;
    case obj_command::face:
    case obj_command::line:
    case obj_command::point:
      if (command == obj_command::face) format_values(fs, "f ");
      if (command == obj_command::line) format_values(fs, "l ");
      if (command == obj_command::point) format_values(fs, "p ");
      for (auto& vert : vertices) format_values(fs, " {}", vert);
      format_values(fs, "\n");
      break;
    case obj_command::object: format_values(fs, "o {}\n", name); break;
    case obj_command::group: format_values(fs, "g {}\n", name); break;
    case obj_command::usemtl: format_values(fs, "usemtl {}\n", name); break;
    case obj_command::smoothing: format_values(fs, "s {}\n", name); break;
    case obj_command::mtllib: format_values(fs, "mtllib {}\n", name); break;
    case obj_command::objxlib: break;
  }
}

void write_mtl_command(
    file_wrapper& fs, mtl_command command, const obj_material& material) {
  // write material
  format_values(fs, "newmtl {}\n", material.name);
  format_values(fs, "illum {}\n", material.illum);
  if (material.emission != zero3f)
    format_values(fs, "Ke {}\n", material.emission);
  if (material.ambient != zero3f)
    format_values(fs, "Ka {}\n", material.ambient);
  format_values(fs, "Kd {}\n", material.diffuse);
  format_values(fs, "Ks {}\n", material.specular);
  if (material.reflection != zero3f)
    format_values(fs, "Kr {}\n", material.reflection);
  if (material.transmission != zero3f)
    format_values(fs, "Kt {}\n", material.transmission);
  format_values(fs, "Ns {}\n", (int)material.exponent);
  if (material.opacity != 1) format_values(fs, "d {}\n", material.opacity);
  if (!material.emission_map.path.empty())
    format_values(fs, "map_Ke {}\n", material.emission_map);
  if (!material.diffuse_map.path.empty())
    format_values(fs, "map_Kd {}\n", material.diffuse_map);
  if (!material.specular_map.path.empty())
    format_values(fs, "map_Ks {}\n", material.specular_map);
  if (!material.transmission_map.path.empty())
    format_values(fs, "map_Kt {}\n", material.transmission_map);
  if (!material.reflection_map.path.empty())
    format_values(fs, "map_Kr {}\n", material.reflection_map);
  if (!material.exponent_map.path.empty())
    format_values(fs, "map_Ns {}\n", material.exponent_map);
  if (!material.opacity_map.path.empty())
    format_values(fs, "map_d {}\n", material.opacity_map);
  if (!material.bump_map.path.empty())
    format_values(fs, "map_bump {}\n", material.bump_map);
  if (!material.displacement_map.path.empty())
    format_values(fs, "map_disp {}\n", material.displacement_map);
  if (!material.normal_map.path.empty())
    format_values(fs, "map_norm {}\n", material.normal_map);
  if (material.pbr_roughness)
    format_values(fs, "Pr {}\n", material.pbr_roughness);
  if (material.pbr_metallic)
    format_values(fs, "Pm {}\n", material.pbr_metallic);
  if (material.pbr_sheen) format_values(fs, "Ps {}\n", material.pbr_sheen);
  if (material.pbr_clearcoat)
    format_values(fs, "Pc {}\n", material.pbr_clearcoat);
  if (material.pbr_coatroughness)
    format_values(fs, "Pcr {}\n", material.pbr_coatroughness);
  if (!material.pbr_roughness_map.path.empty())
    format_values(fs, "map_Pr {}\n", material.pbr_roughness_map);
  if (!material.pbr_metallic_map.path.empty())
    format_values(fs, "map_Pm {}\n", material.pbr_metallic_map);
  if (!material.pbr_sheen_map.path.empty())
    format_values(fs, "map_Ps {}\n", material.pbr_sheen_map);
  if (!material.pbr_clearcoat_map.path.empty())
    format_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map);
  if (!material.pbr_coatroughness_map.path.empty())
    format_values(fs, "map_Pcr {}\n", material.pbr_coatroughness_map);
  if (material.vol_transmission != zero3f)
    format_values(fs, "Vt {}\n", material.vol_transmission);
  if (material.vol_meanfreepath != zero3f)
    format_values(fs, "Vp {}\n", material.vol_meanfreepath);
  if (material.vol_emission != zero3f)
    format_values(fs, "Ve {}\n", material.vol_emission);
  if (material.vol_scattering != zero3f)
    format_values(fs, "Vs {}\n", material.vol_scattering);
  if (material.vol_anisotropy)
    format_values(fs, "Vg {}\n", material.vol_anisotropy);
  if (material.vol_scale) format_values(fs, "Vr {}\n", material.vol_scale);
  if (!material.vol_scattering_map.path.empty())
    format_values(fs, "map_Vs {}\n", material.vol_scattering_map);
  format_values(fs, "\n");
}

void write_objx_command(file_wrapper& fs, objx_command command,
    const obj_camera& camera, const obj_environment& environment,
    const obj_instance& instance) {
  switch (command) {
    case objx_command::camera: {
      format_values(fs, "c {} {} {} {} {} {} {} {}\n", camera.name,
          camera.ortho, camera.width, camera.height, camera.lens, camera.focus,
          camera.aperture, camera.frame);
    } break;
    case objx_command::environment: {
      format_values(fs, "e {} {} {} {}\n", environment.name,
          environment.emission,
          environment.emission_map.path.empty() ? "\"\""s
                                                : environment.emission_map.path,
          environment.frame);
    } break;
    case objx_command::instance: {
      format_values(fs, "i {} {}\n", instance.object, instance.frame);
    } break;
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// YAML SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

static inline void remove_yaml_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

static inline void parse_yaml_varname(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (!is_alpha(str.front())) throw std::runtime_error("cannot parse value");
  auto pos = 0;
  while (is_alpha(str[pos]) || str[pos] == '_' || is_digit(str[pos])) {
    pos += 1;
    if (pos >= str.size()) break;
  }
  value = str.substr(0, pos);
  str.remove_prefix(pos);
}
static inline void parse_yaml_varname(string_view& str, string& value) {
  auto view = ""sv;
  parse_yaml_varname(str, view);
  value = string{view};
}

inline void parse_yaml_value(string_view& str, string_view& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
  } else {
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
}
inline void parse_yaml_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_yaml_value(str, valuev);
  value = string{valuev};
}
inline void parse_yaml_value(string_view& str, int& value) {
  skip_whitespace(str);
  char* end = nullptr;
  value     = (int)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_yaml_value(string_view& str, float& value) {
  skip_whitespace(str);
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_yaml_value(string_view& str, double& value) {
  skip_whitespace(str);
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}

// parse yaml value
void get_yaml_value(const yaml_value& yaml, string& value) {
  if (yaml.type != yaml_value_type::string)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.string_;
}
void get_yaml_value(const yaml_value& yaml, bool& value) {
  if (yaml.type != yaml_value_type::boolean)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.boolean;
}
void get_yaml_value(const yaml_value& yaml, int& value) {
  if (yaml.type != yaml_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (int)yaml.number;
}
void get_yaml_value(const yaml_value& yaml, float& value) {
  if (yaml.type != yaml_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (float)yaml.number;
}
void get_yaml_value(const yaml_value& yaml, vec2f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 2)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1]};
}
void get_yaml_value(const yaml_value& yaml, vec3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 3)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1], (float)yaml.array_[2]};
}
void get_yaml_value(const yaml_value& yaml, mat3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 9)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 9; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}
void get_yaml_value(const yaml_value& yaml, frame3f& value) {
  if (yaml.type != yaml_value_type::array || yaml.number != 12)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 12; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}

// construction
yaml_value make_yaml_value(const string& value) {
  return {yaml_value_type::string, 0, false, value};
}
yaml_value make_yaml_value(bool value) {
  return {yaml_value_type::boolean, 0, value};
}
yaml_value make_yaml_value(int value) {
  return {yaml_value_type::number, (double)value};
}
yaml_value make_yaml_value(float value) {
  return {yaml_value_type::number, (double)value};
}
yaml_value make_yaml_value(const vec2f& value) {
  return {
      yaml_value_type::array, 2, false, "", {(double)value.x, (double)value.y}};
}
yaml_value make_yaml_value(const vec3f& value) {
  return {yaml_value_type::array, 3, false, "",
      {(double)value.x, (double)value.y, (double)value.z}};
}
yaml_value make_yaml_value(const mat3f& value) {
  auto yaml = yaml_value{yaml_value_type::array, 9};
  for (auto i = 0; i < 9; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
}
yaml_value make_yaml_value(const frame3f& value) {
  auto yaml = yaml_value{yaml_value_type::array, 12};
  for (auto i = 0; i < 12; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
}

void parse_yaml_value(string_view& str, yaml_value& value) {
  trim_whitespace(str);
  if (str.empty()) throw std::runtime_error("bad yaml");
  if (str.front() == '[') {
    str.remove_prefix(1);
    value.type   = yaml_value_type::array;
    value.number = 0;
    while (!str.empty()) {
      skip_whitespace(str);
      if (str.empty()) throw std::runtime_error("bad yaml");
      if (str.front() == ']') {
        str.remove_prefix(1);
        break;
      }
      if (value.number >= 16) throw std::runtime_error("array too large");
      parse_yaml_value(str, value.array_[(int)value.number]);
      value.number += 1;
      skip_whitespace(str);
      if (str.front() == ',') {
        str.remove_prefix(1);
        continue;
      } else if (str.front() == ']') {
        str.remove_prefix(1);
        break;
      } else {
        throw std::runtime_error("bad yaml");
      }
    }
  } else if (is_digit(str.front()) || str.front() == '-' ||
             str.front() == '+') {
    value.type = yaml_value_type::number;
    parse_yaml_value(str, value.number);
  } else {
    value.type = yaml_value_type::string;
    parse_yaml_value(str, value.string_);
    if (value.string_ == "true" || value.string_ == "false") {
      value.type    = yaml_value_type::boolean;
      value.boolean = value.string_ == "true";
    }
  }
  skip_whitespace(str);
  if (!str.empty() && !is_whitespace(str)) throw std::runtime_error("bad yaml");
}

// Load/save yaml
void load_yaml(const string& filename, yaml_model& yaml) {
  auto fs = open_file(filename, "rt");

  // read the file line by line
  auto group = ""s;
  auto key   = ""s;
  auto value = yaml_value{};
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_yaml_comment(line);
    if (line.empty()) continue;
    if (is_whitespace(line)) continue;

    // peek commands
    if (is_space(line.front())) {
      // indented property
      if (group == "") throw std::runtime_error("bad yaml");
      skip_whitespace(line);
      if (line.empty()) throw std::runtime_error("bad yaml");
      if (line.front() == '-') {
        auto& element = yaml.elements.emplace_back();
        element.name  = group;
        line.remove_prefix(1);
        skip_whitespace(line);
      } else if (yaml.elements.empty() || yaml.elements.back().name != group) {
        auto& element = yaml.elements.emplace_back();
        element.name  = group;
      }
      parse_yaml_varname(line, key);
      skip_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      parse_yaml_value(line, value);
      yaml.elements.back().key_values.push_back({key, value});
    } else if (is_alpha(line.front())) {
      // new group
      parse_yaml_varname(line, key);
      skip_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      if (!line.empty() && !is_whitespace(line)) {
        group = "";
        if (yaml.elements.empty() || yaml.elements.back().name != group) {
          auto& element = yaml.elements.emplace_back();
          element.name  = group;
        }
        parse_yaml_value(line, value);
        yaml.elements.back().key_values.push_back({key, value});
      } else {
        group = key;
        key   = "";
      }
    } else {
      throw std::runtime_error("bad yaml");
    }
  }
}

static inline void format_value(string& str, const yaml_value& value) {
  switch (value.type) {
    case yaml_value_type::number: format_value(str, value.number); break;
    case yaml_value_type::boolean:
      format_value(str, value.boolean ? "true" : "false");
      break;
    case yaml_value_type::string:
      if (value.string_.empty() || is_digit(value.string_.front())) {
        format_values(str, "\"{}\"", value.string_);
      } else {
        format_values(str, "{}", value.string_);
      }
      break;
    case yaml_value_type::array:
      format_value(str, "[ ");
      for (auto i = 0; i < value.number; i++) {
        if (i) format_value(str, ", ");
        format_value(str, value.array_[i]);
      }
      format_value(str, " ]");
      break;
  }
}

void save_yaml(const string& filename, const yaml_model& yaml) {
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : yaml.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  auto group = ""s;
  for (auto& element : yaml.elements) {
    if (group != element.name) {
      group = element.name;
      if (group != "") {
        format_values(fs, "\n{}:\n", group);
      } else {
        format_values(fs, "\n");
      }
      auto first = true;
      for (auto& [key, value] : element.key_values) {
        if (group != "") {
          format_values(fs, "  {} {}: {}\n", first ? "-" : " ", key, value);
          first = false;
        } else {
          format_values(fs, "{}: {}\n", key, value);
        }
      }
    }
  }
}

bool read_yaml_property(file_wrapper& fs, string& group, string& key,
    bool& newobj, yaml_value& value) {
  // read the file line by line
  char buffer[4096];
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_yaml_comment(line);
    if (line.empty()) continue;
    if (is_whitespace(line)) continue;

    // peek commands
    if (is_space(line.front())) {
      // indented property
      if (group == "") throw std::runtime_error("bad yaml");
      skip_whitespace(line);
      if (line.empty()) throw std::runtime_error("bad yaml");
      if (line.front() == '-') {
        newobj = true;
        line.remove_prefix(1);
        skip_whitespace(line);
      } else {
        newobj = false;
      }
      parse_yaml_varname(line, key);
      skip_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      parse_yaml_value(line, value);
      return true;
    } else if (is_alpha(line.front())) {
      // new group
      parse_yaml_varname(line, key);
      skip_whitespace(line);
      if (line.empty() || line.front() != ':')
        throw std::runtime_error("bad yaml");
      line.remove_prefix(1);
      if (!line.empty() && !is_whitespace(line)) {
        group = "";
        parse_yaml_value(line, value);
        return true;
      } else {
        group = key;
        key   = "";
        return true;
      }
    } else {
      throw std::runtime_error("bad yaml");
    }
  }
  return false;
}

void write_yaml_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    format_values(fs, "# {}\n", line);
  }
  format_values(fs, "\n");
}

// Save yaml property
void write_yaml_property(file_wrapper& fs, const string& object,
    const string& key, bool newobj, const yaml_value& value) {
  if (key.empty()) {
    format_values(fs, "\n{}:\n", object);
  } else {
    if (!object.empty()) {
      format_values(fs, "  {} {}: {}\n", newobj ? "-" : " ", key, value);
    } else {
      format_values(fs, "{}: {}\n", key, value);
    }
  }
}

void write_yaml_object(file_wrapper& fs, const string& object) {
  format_values(fs, "\n{}:\n", object);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static inline void remove_pbrt_comment(
    string_view& str, char comment_char = '#') {
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
bool read_pbrt_cmdline(file_wrapper& fs, string& cmd, int& line_num) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs.fs);
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    line_num += 1;
    auto line = string_view{buffer};
    remove_pbrt_comment(line);
    skip_whitespace(line);
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

// parse a quoted string
static inline void parse_pbrt_value(string_view& str, string_view& value) {
  skip_whitespace(str);
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

static inline void parse_pbrt_value(string_view& str, string& value) {
  auto view = ""sv;
  parse_pbrt_value(str, view);
  value = string{view};
}

// parse a quoted string
static inline void parse_pbrt_command(string_view& str, string& value) {
  skip_whitespace(str);
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
static inline void parse_pbrt_value(string_view& str, float& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtof(str.data(), &next);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}

// parse a number
static inline void parse_pbrt_value(string_view& str, int& value) {
  skip_whitespace(str);
  if (str.empty()) throw std::runtime_error("number expected");
  auto next = (char*)nullptr;
  value     = strtol(str.data(), &next, 10);
  if (str.data() == next) throw std::runtime_error("number expected");
  str.remove_prefix(next - str.data());
}
template <typename T>
static inline void parse_pbrt_value(
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
static inline void parse_pbrt_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}

// parse pbrt value with optional parens
template <typename T>
static inline void parse_pbrt_param(string_view& str, T& value) {
  skip_whitespace(str);
  auto parens = !str.empty() && str.front() == '[';
  if (parens) str.remove_prefix(1);
  parse_pbrt_value(str, value);
  if (parens) {
    skip_whitespace(str);
    if (!str.empty() && str.front() == '[')
      throw std::runtime_error("bad pbrt param");
    str.remove_prefix(1);
  }
}

// parse a quoted string
static inline void parse_pbrt_nametype(
    string_view& str_, string& name, string& type) {
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

static inline pair<vec3f, vec3f> get_pbrt_etak(const string& name) {
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

static inline void parse_pbrt_params(
    string_view& str, vector<pbrt_value>& values) {
  auto parse_pbrt_pvalues = [](string_view& str, auto& value, auto& values) {
    values.clear();
    skip_whitespace(str);
    if (str.empty()) throw std::runtime_error("bad pbrt value");
    if (str.front() == '[') {
      str.remove_prefix(1);
      skip_whitespace(str);
      if (str.empty()) throw std::runtime_error("bad pbrt value");
      while (!str.empty()) {
        auto& val = values.empty() ? value : values.emplace_back();
        parse_pbrt_value(str, val);
        skip_whitespace(str);
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
  skip_whitespace(str);
  while (!str.empty()) {
    auto& value = values.emplace_back();
    auto  type  = ""s;
    parse_pbrt_nametype(str, value.name, type);
    skip_whitespace(str);
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
      // TODO: xyz conversion
      value.type = pbrt_value_type::color;
      parse_pbrt_pvalues(str, value.value3f, value.vector3f);
      throw std::runtime_error("xyz conversion");
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
        parse_pbrt_value(str, filename);
        auto filenamep = get_filename(filename);
        if (get_extension(filenamep) == ".spd") {
          filenamep = replace_extension(filenamep, "");
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (get_extension(filenamep) == ".eta") {
            auto eta = get_pbrt_etak(replace_extension(filenamep, "")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (get_extension(filenamep) == ".k") {
            auto k = get_pbrt_etak(replace_extension(filenamep, "")).second;
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
    skip_whitespace(str);
  }
}

// convert pbrt films
static void convert_pbrt_films(vector<pbrt_film>& films, bool verbose = false) {
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
static void convert_pbrt_cameras(vector<pbrt_camera>& cameras,
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
      camera.fov = get_pbrt_value(values, "fov", 90.0f);
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
      camera.fov      = 2 * atan(0.036f / (2 * lens));
      camera.aperture = get_pbrt_value(values, "aperturediameter", 0.0f);
      camera.focus    = get_pbrt_value(values, "focusdistance", 10.0f);
      camera.aspect   = film_aspect;
    } else {
      throw std::runtime_error("unsupported Camera type " + camera.type);
    }
  }
}

// convert pbrt textures
static void convert_pbrt_textures(
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
static void convert_pbrt_materials(vector<pbrt_material>& materials,
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
static inline void make_pbrt_shape(vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps, const PositionFunc& position_func,
    const NormalFunc& normal_func) {
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
static inline void make_pbrt_sphere(vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps, float radius) {
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
static inline void make_pbrt_disk(vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps, float radius) {
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
static inline void make_pbrt_quad(vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    const vec2i& steps, float radius) {
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
static void convert_pbrt_shapes(
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
static void convert_pbrt_arealights(
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
static void convert_pbrt_lights(
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

static void convert_pbrt_environments(vector<pbrt_environment>& environments,
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
void load_pbrt(const string& filename, pbrt_model& pbrt) {
  auto files = vector<file_wrapper>{};
  open_file(files.emplace_back(), filename);

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
        open_file(files.emplace_back(), get_dirname(filename) + includename);
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

inline static void format_value(string& str, const pbrt_value& value) {
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

inline static void format_value(string& str, const vector<pbrt_value>& values) {
  for (auto& value : values) {
    str += " ";
    format_value(str, value);
  }
}

void save_pbrt(const string& filename, const pbrt_model& pbrt) {
  auto fs = open_file(filename, "wt");

  // save comments
  format_values(fs, "#\n");
  format_values(fs, "# Written by Yocto/GL\n");
  format_values(fs, "# https://github.com/xelatihy/yocto-gl\n");
  format_values(fs, "#\n\n");
  for (auto& comment : pbrt.comments) {
    format_values(fs, "# {}\n", comment);
  }
  format_values(fs, "\n");

  for (auto& camera_ : pbrt.cameras) {
    auto camera = camera_;
    if (camera.type == "") {
      camera.type = "perspective";
      camera.values.push_back(make_pbrt_value("fov", camera.fov * 180 / pif));
    }
    format_values(fs, "LookAt {} {} {}\n", camera.frame.o,
        camera.frame.o - camera.frame.z, camera.frame.y);
    format_values(fs, "Camera \"{}\" {}\n", camera.type, camera.values);
  }
  for (auto& film_ : pbrt.films) {
    auto film = film_;
    if (film.type == "") {
      film.type = "image";
      film.values.push_back(make_pbrt_value("xresolution", film.resolution.x));
      film.values.push_back(make_pbrt_value("yresolution", film.resolution.y));
      film.values.push_back(make_pbrt_value("filename", film.filename));
    }
    format_values(fs, "Film \"{}\" {}\n", film.type, film.values);
  }
  for (auto& integrator_ : pbrt.integrators) {
    auto integrator = integrator_;
    format_values(
        fs, "Integrator \"{}\" {}\n", integrator.type, integrator.values);
  }
  for (auto& sampler_ : pbrt.samplers) {
    auto sampler = sampler_;
    format_values(fs, "Sampler \"{}\" {}\n", sampler.type, sampler.values);
  }
  for (auto& filter_ : pbrt.filters) {
    auto filter = filter_;
    format_values(fs, "PixelFilter \"{}\" {}\n", filter.type, filter.values);
  }
  for (auto& accelerator_ : pbrt.accelerators) {
    auto accelerator = accelerator_;
    format_values(
        fs, "Accelerator \"{}\" {}\n", accelerator.type, accelerator.values);
  }

  format_values(fs, "\nWorldBegin\n\n");

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
    format_values(fs, "Texture \"{}\" \"color\" \"{}\" {}\n", texture.name,
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
    format_values(fs, "MakeNamedMaterial \"{}\" \"string type\" \"{}\" {}\n",
        material.name, material.type, material.values);
  }

  for (auto& medium_ : pbrt.mediums) {
    auto medium = medium_;
    format_values(fs, "MakeNamedMedium \"{}\" \"string type\" \"{}\" {}\n",
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
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", (mat4f)light.frame);
    format_values(fs, "LightSource \"{}\" {}\n", light.type, light.values);
    format_values(fs, "AttributeEnd\n");
  }

  for (auto& environment_ : pbrt.environments) {
    auto environment = environment_;
    if (environment.type == "") {
      environment.type = "infinite";
      environment.values.push_back(make_pbrt_value("L", environment.emission));
      environment.values.push_back(
          make_pbrt_value("mapname", environment.filename));
    }
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", (mat4f)environment.frame);
    // TODO: should we correct frames?
    format_values(
        fs, "LightSource \"{}\" {}\n", environment.type, environment.values);
    format_values(fs, "AttributeEnd\n");
  }

  auto arealights_map = unordered_map<string, string>{};
  for (auto& arealight_ : pbrt.arealights) {
    auto arealight = arealight_;
    if (arealight.type == "") {
      arealight.type = "diffuse";
      arealight.values.push_back(make_pbrt_value("L", arealight.emission));
    }
    format_values(arealights_map[arealight.name], "AreaLightSource \"{}\" {}\n",
        arealight.type, arealight.values);
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
    if (shape.is_instanced) format_values(fs, "ObjectBegin \"{}\"\n", object);
    format_values(fs, "AttributeBegin\n");
    format_values(fs, "Transform {}\n", (mat4f)shape.frame);
    format_values(fs, "NamedMaterial \"{}\"\n", shape.material);
    if (shape.arealight != "")
      format_values(fs, arealights_map.at(shape.arealight));
    format_values(fs, "Shape \"{}\" {}\n", shape.type, shape.values);
    format_values(fs, "AttributeEnd\n");
    if (shape.is_instanced) format_values(fs, "ObjectEnd\n");
    for (auto& iframe : shape.instance_frames) {
      format_values(fs, "AttributeBegin\n");
      format_values(fs, "Transform {}\n", (mat4f)iframe);
      format_values(fs, "ObjectInstance \"{}\"\n", object);
      format_values(fs, "AttributeEnd\n");
    }
  }

  format_values(fs, "\nWorldEnd\n\n");
}

// Read pbrt commands
bool read_pbrt_command(file_wrapper& fs, pbrt_command& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values, string& line) {
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
bool read_pbrt_command(file_wrapper& fs, pbrt_command& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values) {
  auto command_buffer = ""s;
  return read_pbrt_command(
      fs, command, name, type, xform, values, command_buffer);
}

// Write obj elements
void write_pbrt_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    format_values(fs, "# {}\n", line);
  }
  format_values(fs, "\n");
}

void write_pbrt_values(file_wrapper& fs, const vector<pbrt_value>& values) {
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
    format_values(fs, " \"{} {}\" ", type_labels.at(value.type), value.name);
    switch (value.type) {
      case pbrt_value_type::real:
        if (!value.vector1f.empty()) {
          format_values(fs, "[ ");
          for (auto& v : value.vector1f) format_values(fs, " {}", v);
          format_values(fs, " ]");
        } else {
          format_values(fs, "{}", value.value1f);
        }
        break;
      case pbrt_value_type::integer:
        if (!value.vector1f.empty()) {
          format_values(fs, "[ ");
          for (auto& v : value.vector1i) format_values(fs, " {}", v);
          format_values(fs, " ]");
        } else {
          format_values(fs, "{}", value.value1i);
        }
        break;
      case pbrt_value_type::boolean:
        format_values(fs, "\"{}\"", value.value1b ? "true" : "false");
        break;
      case pbrt_value_type::string:
        format_values(fs, "\"{}\"", value.value1s);
        break;
      case pbrt_value_type::point:
      case pbrt_value_type::vector:
      case pbrt_value_type::normal:
      case pbrt_value_type::color:
        if (!value.vector3f.empty()) {
          format_values(fs, "[ ");
          for (auto& v : value.vector3f) format_values(fs, " {}", v);
          format_values(fs, " ]");
        } else {
          format_values(fs, "[ {} ]", value.value3f);
        }
        break;
      case pbrt_value_type::spectrum:
        format_values(fs, "[ ");
        for (auto& v : value.vector1f) format_values(fs, " {}", v);
        format_values(fs, " ]");
        break;
      case pbrt_value_type::texture:
        format_values(fs, "\"{}\"", value.value1s);
        break;
      case pbrt_value_type::point2:
      case pbrt_value_type::vector2:
        if (!value.vector2f.empty()) {
          format_values(fs, "[ ");
          for (auto& v : value.vector2f) format_values(fs, " {}", v);
          format_values(fs, " ]");
        } else {
          format_values(fs, "[ {} ]", value.value2f);
        }
        break;
    }
  }
  format_values(fs, "\n");
}

void write_pbrt_command(file_wrapper& fs, pbrt_command command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_float) {
  switch (command) {
    case pbrt_command::world_begin: format_values(fs, "WorldBegin\n"); break;
    case pbrt_command::world_end: format_values(fs, "WorldEnd\n"); break;
    case pbrt_command::attribute_begin:
      format_values(fs, "AttributeBegin\n");
      break;
    case pbrt_command::attribute_end:
      format_values(fs, "AttributeEnd\n");
      break;
    case pbrt_command::transform_begin:
      format_values(fs, "TransformBegin\n");
      break;
    case pbrt_command::transform_end:
      format_values(fs, "TransformEnd\n");
      break;
    case pbrt_command::object_begin:
      format_values(fs, "ObjectBegin \"{}\"\n", name);
      break;
    case pbrt_command::object_end: format_values(fs, "ObjectEnd\n"); break;
    case pbrt_command::object_instance:
      format_values(fs, "ObjectInstance \"{}\"\n", name);
      break;
    case pbrt_command::sampler:
      format_values(fs, "Sampler \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::integrator:
      format_values(fs, "Integrator \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::accelerator:
      format_values(fs, "Accelerator \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::film:
      format_values(fs, "Film \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::filter:
      format_values(fs, "Filter \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::camera:
      format_values(fs, "Camera \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::shape:
      format_values(fs, "Shape \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::light:
      format_values(fs, "LightSource \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::material:
      format_values(fs, "Material \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::arealight:
      format_values(fs, "AreaLightSource \"{}\"", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_texture:
      format_values(fs, "Texture \"{}\" \"{}\" \"{}\"", name,
          texture_float ? "float" : "rgb", type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_medium:
      format_values(
          fs, "MakeNamedMedium \"{}\" \"string type\" \"{}\"", name, type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::named_material:
      format_values(
          fs, "MakeNamedMaterial \"{}\" \"string type\" \"{}\"", name, type);
      write_pbrt_values(fs, values);
      break;
    case pbrt_command::include:
      format_values(fs, "Include \"{}\"\n", name);
      break;
    case pbrt_command::reverse_orientation:
      format_values(fs, "ReverseOrientation\n");
      break;
    case pbrt_command::set_transform:
      format_values(fs, "Transform {}\n", (mat4f)xform);
      break;
    case pbrt_command::concat_transform:
      format_values(fs, "ConcatTransform {}\n", (mat4f)xform);
      break;
    case pbrt_command::lookat_transform:
      format_values(fs, "LookAt {} {} {}\n", xform.x, xform.y, xform.z);
      break;
    case pbrt_command::use_material:
      format_values(fs, "NamedMaterial \"{}\"\n", name);
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
      format_values(fs, "MediumInterface \"{}\" \"{}\"\n", interior, exterior);
    } break;
    case pbrt_command::active_transform:
      format_values(fs, "ActiveTransform \"{}\"\n", name);
      break;
    case pbrt_command::coordinate_system_set:
      format_values(fs, "CoordinateSystem \"{}\"\n", name);
      break;
    case pbrt_command::coordinate_system_transform:
      format_values(fs, "CoordinateSysTransform \"{}\"\n", name);
      break;
  }
}

void write_pbrt_command(file_wrapper& fs, pbrt_command command,
    const string& name, const frame3f& xform) {
  return write_pbrt_command(fs, command, name, "", xform, {});
}
void write_pbrt_command(file_wrapper& fs, pbrt_command command,
    const string& name, const string& type, const vector<pbrt_value>& values,
    bool texture_as_float) {
  return write_pbrt_command(
      fs, command, name, type, identity3x4f, values, texture_as_float);
}

// get pbrt value
void get_pbrt_value(const pbrt_value& pbrt, string& value) {
  if (pbrt.type == pbrt_value_type::string ||
      pbrt.type == pbrt_value_type::texture) {
    value = pbrt.value1s;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
void get_pbrt_value(const pbrt_value& pbrt, bool& value) {
  if (pbrt.type == pbrt_value_type::boolean) {
    value = pbrt.value1b;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
void get_pbrt_value(const pbrt_value& pbrt, int& value) {
  if (pbrt.type == pbrt_value_type::integer) {
    value = pbrt.value1i;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
void get_pbrt_value(const pbrt_value& pbrt, float& value) {
  if (pbrt.type == pbrt_value_type::real) {
    value = pbrt.value1f;
  } else {
    throw std::runtime_error("bad pbrt type");
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vec2f& value) {
  if (pbrt.type == pbrt_value_type::point2 ||
      pbrt.type == pbrt_value_type::vector2) {
    value = pbrt.value2f;
  } else {
    throw std::runtime_error("bad pbrt type");
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
    throw std::runtime_error("bad pbrt type");
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
    throw std::runtime_error("bad pbrt type");
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
    throw std::runtime_error("bad pbrt type");
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
      throw std::runtime_error("bad pbrt type");
    value.resize(pbrt.vector1f.size() / 3);
    for (auto i = 0; i < value.size(); i++)
      value[i] = {pbrt.vector1f[i * 3 + 0], pbrt.vector1f[i * 3 + 1],
          pbrt.vector1f[i * 3 + 2]};
  } else {
    throw std::runtime_error("bad pbrt type");
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
    throw std::runtime_error("bad pbrt type");
  }
}
void get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value) {
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
void get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value) {
  if (pbrt.type == pbrt_value_type::string) {
    value.first = 0;
    get_pbrt_value(pbrt, value.second);
  } else {
    get_pbrt_value(pbrt, value.first);
    value.second = "";
  }
}
void get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value) {
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
pbrt_value make_pbrt_value(
    const string& name, const string& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1s = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, bool value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1b = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, int value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1i = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, float value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value1f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, const vec2f& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value2f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, const vec3f& value, pbrt_value_type type) {
  auto pbrt    = pbrt_value{};
  pbrt.name    = name;
  pbrt.type    = type;
  pbrt.value3f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, const vector<vec2f>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector2f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, const vector<vec3f>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector3f = value;
  return pbrt;
}
pbrt_value make_pbrt_value(
    const string& name, const vector<vec3i>& value, pbrt_value_type type) {
  auto pbrt     = pbrt_value{};
  pbrt.name     = name;
  pbrt.type     = type;
  pbrt.vector1i = {(int*)value.data(), (int*)value.data() + value.size() * 3};
  return pbrt;
}

// old code --- maintained here in case we want to integrate back
#if 0
void approximate_fourier_material(pbrt_material::fourier_t& fourier) {
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
static inline pair<vec3f, vec3f> parse_pbrt_subsurface(const string& name) {
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
    throw std::runtime_error("could not load " + filename);
  }
  auto gltf = std::unique_ptr<cgltf_data, void (*)(cgltf_data*)>{
      data, cgltf_free};
  auto dirname = get_dirname(filename);
  if (dirname != "") dirname += "/";
  if (cgltf_load_buffers(&params, data, dirname.c_str()) !=
      cgltf_result_success) {
    throw std::runtime_error("could not load gltf buffers " + filename);
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
