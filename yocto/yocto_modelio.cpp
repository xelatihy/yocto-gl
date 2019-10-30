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
#include <string_view>

// -----------------------------------------------------------------------------
// LOW-LEVEL FILE HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
ply_file::ply_file(ply_file&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
ply_file::~ply_file() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
void open_ply(ply_file& fs, const string& filename, const string& mode) {
  close_ply(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
ply_file open_ply(const string& filename, const string& mode) {
  auto fs = ply_file{};
  open_ply(fs, filename, mode);
  return fs;
}
void close_ply(ply_file& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

// Read a line
bool read_ply_line(ply_file& fs, char* buffer, size_t size) {
  if (fgets(buffer, size, fs.fs)) {
    fs.linenum += 1;
    return true;
  } else {
    return false;
  }
}

template <typename T>
static T swap_ply_endian(T value) {
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
bool write_ply_value(ply_file& fs, const T& value) {
  return fwrite(&value, sizeof(value), 1, fs.fs) == 1;
}
template <typename T>
bool write_ply_value(ply_file& fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_ply_endian(value_) : value_;
  return fwrite(&value, sizeof(value), 1, fs.fs) == 1;
}

template <typename T>
bool read_ply_value(ply_file& fs, T& value) {
  return fread(&value, sizeof(value), 1, fs.fs) == 1;
}
template <typename T>
bool read_ply_value(ply_file& fs, T& value, bool big_endian) {
  auto ok = fread(&value, sizeof(value), 1, fs.fs) == 1;
  if (big_endian) value = swap_ply_endian(value);
  return ok;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOAD-LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;

// utilities
static bool is_ply_newline(char c) { return c == '\r' || c == '\n'; }
static bool is_ply_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static void skip_ply_whitespace(string_view& str) {
  while (!str.empty() && is_ply_space(str.front())) str.remove_prefix(1);
}

static void remove_ply_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_ply_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
static bool parse_ply_value(string_view& str, string_view& value) {
  skip_ply_whitespace(str);
  if (str.empty()) return false;
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_ply_space(cpy.front())) cpy.remove_prefix(1);
    value = str;
    value.remove_suffix(cpy.size());
    str.remove_prefix(str.size() - cpy.size());
    return true;
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
    return true;
  }
}
static bool parse_ply_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_ply_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
static bool parse_ply_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_ply_value(str, valuei);
  value = (bool)valuei;
  return true;
}
static bool parse_ply_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_ply_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#ifdef __APPLE__
static bool parse_ply_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL PRINTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
static void format_ply_value(string& str, const string& value) { str += value; }
static void format_ply_value(string& str, const char* value) { str += value; }
static void format_ply_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_ply_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_ply_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_ply_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
static void format_ply_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_ply_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_ply_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
static void format_ply_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
static void format_ply_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static void format_ply_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}

// Foramt to file
static void format_ply_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static void format_ply_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_ply_value(str, arg);
  format_ply_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
static bool format_ply_values(
    ply_file& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_ply_values(str, fmt, args...);
  return fputs(str.c_str(), fs.fs) >= 0;
}
template <typename T>
static bool format_ply_value(ply_file& fs, const T& value) {
  auto str = ""s;
  format_ply_value(str, value);
  return fputs(str.c_str(), fs.fs) >= 0;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply
plyio_status load_ply(const string& filename, ply_model& ply) {
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
  auto fs = open_ply(filename, "rb");
  if (!fs) return {filename + ": file not found"};

  // read header ---------------------------------------------
  char buffer[4096];
  while (read_ply_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_ply_comment(str);
    skip_ply_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_ply_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // check magic number
    if (first_line) {
      if (cmd != "ply") return {filename + ": not a ply file"};
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) return {filename + ": corrupt header"};
    } else if (cmd == "format") {
      auto fmt = ""s;
      if (!parse_ply_value(str, fmt)) return {filename + ": parse error"};
      if (fmt == "ascii") {
        ply.format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        ply.format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        ply.format = ply_format::binary_big_endian;
      } else {
        return {filename + ": unknown format " + fmt};
      }
    } else if (cmd == "comment") {
      skip_ply_whitespace(str);
      ply.comments.push_back(string{str});
    } else if (cmd == "obj_info") {
      skip_ply_whitespace(str);
      // comment is the rest of the str
    } else if (cmd == "element") {
      auto& elem = ply.elements.emplace_back();
      if (!parse_ply_value(str, elem.name)) return {filename + ": parse error"};
      if (!parse_ply_value(str, elem.count))
        return {filename + ": parse error"};
    } else if (cmd == "property") {
      if (ply.elements.empty()) return {filename + ": corrupt header"};
      auto& prop  = ply.elements.back().properties.emplace_back();
      auto  tname = ""s;
      if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
      if (tname == "list") {
        prop.is_list = true;
        if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          return {filename + ": unsupported list size type " + tname};
        if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
        if (type_map.find(tname) == type_map.end())
          return {filename + ": unknown type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          return {filename + ": unknown type " + tname};
        prop.type = type_map.at(tname);
      }
      if (!parse_ply_value(str, prop.name)) return {filename + ": parse error"};
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      return {filename + ": unknown command " + cmd};
    }
  }

  // check exit
  if (!end_header) return {filename + ": incomplete header"};

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
        if (!read_ply_line(fs, buffer, sizeof(buffer)))
          return {filename + ": read error"};
        auto str = string_view{buffer};
        for (auto& prop : elem.properties) {
          if (prop.is_list) {
            if (!parse_ply_value(str, prop.ldata_u8.emplace_back()))
              return {filename + ": parse error"};
          }
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                if (!parse_ply_value(str, prop.data_i8.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::i16:
                if (!parse_ply_value(str, prop.data_i16.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::i32:
                if (!parse_ply_value(str, prop.data_i32.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::i64:
                if (!parse_ply_value(str, prop.data_i64.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::u8:
                if (!parse_ply_value(str, prop.data_u8.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::u16:
                if (!parse_ply_value(str, prop.data_u16.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::u32:
                if (!parse_ply_value(str, prop.data_u32.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::u64:
                if (!parse_ply_value(str, prop.data_u64.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::f32:
                if (!parse_ply_value(str, prop.data_f32.emplace_back()))
                  return {filename + ": parse error"};
                break;
              case ply_type::f64:
                if (!parse_ply_value(str, prop.data_f64.emplace_back()))
                  return {filename + ": parse error"};
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
          if (prop.is_list) {
            if (!read_ply_value(fs, prop.ldata_u8.emplace_back(), big_endian))
              return {filename + ": read error"};
          }
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                if (!read_ply_value(
                        fs, prop.data_i8.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::i16:
                if (!read_ply_value(
                        fs, prop.data_i16.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::i32:
                if (!read_ply_value(
                        fs, prop.data_i32.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::i64:
                if (!read_ply_value(
                        fs, prop.data_i64.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::u8:
                if (!read_ply_value(
                        fs, prop.data_u8.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::u16:
                if (!read_ply_value(
                        fs, prop.data_u16.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::u32:
                if (!read_ply_value(
                        fs, prop.data_u32.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::u64:
                if (!read_ply_value(
                        fs, prop.data_u64.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::f32:
                if (!read_ply_value(
                        fs, prop.data_f32.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
              case ply_type::f64:
                if (!read_ply_value(
                        fs, prop.data_f64.emplace_back(), big_endian))
                  return {filename + ": read error"};
                break;
            }
          }
        }
      }
    }
  }

  return {};
}

// Save ply
plyio_status save_ply(const string& filename, const ply_model& ply) {
  auto fs = open_ply(filename, "wb");
  if (!fs) return {filename + ": file not found"};

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
  if (!format_ply_values(fs, "ply\n")) return {filename + ": write error"};
  if (!format_ply_values(fs, "format {} 1.0\n", format_map.at(ply.format)))
    return {filename + ": write error"};
  if (!format_ply_values(fs, "comment Written by Yocto/GL\n"))
    return {filename + ": write error"};
  if (!format_ply_values(fs, "comment https://github.com/xelatihy/yocto-gl\n"))
    return {filename + ": write error"};
  for (auto& comment : ply.comments)
    if (!format_ply_values(fs, "comment {}\n", comment))
      return {filename + ": write error"};
  for (auto& elem : ply.elements) {
    if (!format_ply_values(
            fs, "element {} {}\n", elem.name, (uint64_t)elem.count))
      return {filename + ": write error"};
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        if (!format_ply_values(fs, "property list uchar {} {}\n",
                type_map[prop.type], prop.name))
          return {filename + ": write error"};
      } else {
        if (!format_ply_values(
                fs, "property {} {}\n", type_map[prop.type], prop.name))
          return {filename + ": write error"};
      }
    }
  }
  if (!format_ply_values(fs, "end_header\n"))
    return {filename + ": write error"};

  // properties
  if (ply.format == ply_format::ascii) {
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto pidx = 0; pidx < elem.properties.size(); pidx++) {
          auto& prop = elem.properties[pidx];
          if (prop.is_list)
            if (!format_ply_values(fs, "{} ", (int)prop.ldata_u8[idx]))
              return {filename + ": write error"};
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                if (!format_ply_values(fs, "{} ", prop.data_i8[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::i16:
                if (!format_ply_values(fs, "{} ", prop.data_i16[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::i32:
                if (!format_ply_values(fs, "{} ", prop.data_i32[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::i64:
                if (!format_ply_values(fs, "{} ", prop.data_i64[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::u8:
                if (!format_ply_values(fs, "{} ", prop.data_u8[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::u16:
                if (!format_ply_values(fs, "{} ", prop.data_u16[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::u32:
                if (!format_ply_values(fs, "{} ", prop.data_u32[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::u64:
                if (!format_ply_values(fs, "{} ", prop.data_u64[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::f32:
                if (!format_ply_values(fs, "{} ", prop.data_f32[cur[idx]++]))
                  return {filename + ": write error"};
                break;
              case ply_type::f64:
                if (!format_ply_values(fs, "{} ", prop.data_f64[cur[idx]++]))
                  return {filename + ": write error"};
                break;
            }
          }
          if (!format_ply_values(fs, "\n")) return {filename + ": write error"};
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
          if (prop.is_list)
            if (!write_ply_value(fs, prop.ldata_u8[idx], big_endian))
              return {filename + ": write error"};
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                if (!write_ply_value(fs, prop.data_i8[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::i16:
                if (!write_ply_value(
                        fs, prop.data_i16[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::i32:
                if (!write_ply_value(
                        fs, prop.data_i32[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::i64:
                if (!write_ply_value(
                        fs, prop.data_i64[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::u8:
                if (!write_ply_value(fs, prop.data_u8[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::u16:
                if (!write_ply_value(
                        fs, prop.data_u16[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::u32:
                if (!write_ply_value(
                        fs, prop.data_u32[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::u64:
                if (!write_ply_value(
                        fs, prop.data_u64[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::f32:
                if (!write_ply_value(
                        fs, prop.data_f32[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
              case ply_type::f64:
                if (!write_ply_value(
                        fs, prop.data_f64[cur[pidx]++], big_endian))
                  return {filename + ": write error"};
                break;
            }
          }
        }
      }
    }
  }

  return {};
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
static vector<T> convert_ply_property(const vector<T1>& prop) {
  auto values = vector<T>(prop.size());
  for (auto i = (size_t)0; i < prop.size(); i++) values[i] = (T)prop[i];
  return values;
}
template <typename T>
static vector<T> convert_ply_property(const ply_property& prop) {
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
    default: throw std::runtime_error("should not be here");
  }
  // return here to silence warnings
  std::runtime_error("should not have gotten here");
  return {};
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

static vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
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
  return flipv ? flip_ply_texcoord(texcoord) : texcoord;
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
  auto indices = get_ply_list_values(ply, "str", "vertex_indices");
  auto sizes   = get_ply_list_sizes(ply, "str", "vertex_indices");
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
static void add_ply_element(
    ply_model& ply, const string& element, size_t count) {
  for (auto& elem : ply.elements) {
    if (elem.name == element) return;
  }
  auto& elem = ply.elements.emplace_back();
  elem.name  = element;
  elem.count = count;
}
static void add_ply_property(ply_model& ply, const string& element,
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
static vector<T> make_ply_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

static void add_ply_values(ply_model& ply, const float* values, size_t count,
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
void add_ply_lists(ply_model& ply, const int* values, size_t count,
    int size, const string& element, const string& property) {
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
      ply, flipv ? flip_ply_texcoord(values) : values, "vertex", "u", "v");
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
  return add_ply_lists(ply, values, "str", "vertex_indices");
}
void add_ply_points(ply_model& ply, const vector<int>& values) {
  return add_ply_lists(ply, values, "point", "vertex_indices");
}

// get ply value either ascii or binary
template <typename T, typename VT>
[[nodiscard]] static bool read_ply_prop(
    ply_file& fs, VT& value, bool big_endian) {
  auto tvalue = T{};
  auto ok     = read_ply_value(fs, tvalue, big_endian);
  value       = (VT)tvalue;
  return ok;
}
template <typename VT>
[[nodiscard]] static bool read_ply_prop(
    ply_file& fs, ply_type type, VT& value, bool big_endian) {
  switch (type) {
    case ply_type::i8: return read_ply_prop<int8_t>(fs, value, big_endian);
    case ply_type::i16: return read_ply_prop<int16_t>(fs, value, big_endian);
    case ply_type::i32: return read_ply_prop<int32_t>(fs, value, big_endian);
    case ply_type::i64: return read_ply_prop<int64_t>(fs, value, big_endian);
    case ply_type::u8: return read_ply_prop<uint8_t>(fs, value, big_endian);
    case ply_type::u16: return read_ply_prop<uint16_t>(fs, value, big_endian);
    case ply_type::u32: return read_ply_prop<uint32_t>(fs, value, big_endian);
    case ply_type::u64: return read_ply_prop<uint64_t>(fs, value, big_endian);
    case ply_type::f32: return read_ply_prop<float>(fs, value, big_endian);
    case ply_type::f64: return read_ply_prop<double>(fs, value, big_endian);
  }
}

template <typename T, typename VT>
static bool parse_ply_prop(string_view& str, VT& value) {
  auto tvalue = T{};
  if (!parse_ply_value(str, tvalue)) return false;
  value = (VT)tvalue;
  return false;
}
template <typename VT>
static bool parse_ply_prop(string_view& str, ply_type type, VT& value) {
  switch (type) {
    case ply_type::i8: return parse_ply_prop<int8_t>(str, value);
    case ply_type::i16: return parse_ply_prop<int16_t>(str, value);
    case ply_type::i32: return parse_ply_prop<int32_t>(str, value);
    case ply_type::i64: return parse_ply_prop<int64_t>(str, value);
    case ply_type::u8: return parse_ply_prop<uint8_t>(str, value);
    case ply_type::u16: return parse_ply_prop<uint16_t>(str, value);
    case ply_type::u32: return parse_ply_prop<uint32_t>(str, value);
    case ply_type::u64: return parse_ply_prop<uint64_t>(str, value);
    case ply_type::f32: return parse_ply_prop<float>(str, value);
    case ply_type::f64: return parse_ply_prop<double>(str, value);
    default: return false;
  }
}

// Load ply data
plyio_status read_ply_header(const string& filename, ply_file& fs,
    ply_format& format, vector<ply_element>& elements,
    vector<string>& comments) {
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

  // read the file header str by str
  char buffer[4096];
  while (read_ply_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_ply_comment(str);
    skip_ply_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_ply_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // check magic number
    if (first_line) {
      if (cmd != "ply") return {filename + ": format not ply"};
      first_line = false;
      continue;
    }

    // possible token values
    if (cmd == "ply") {
      if (!first_line) return {filename + ": corrupt header"};
    } else if (cmd == "format") {
      auto fmt = string_view{};
      if (!parse_ply_value(str, fmt)) return {filename + ": parse error"};
      if (fmt == "ascii") {
        format = ply_format::ascii;
      } else if (fmt == "binary_little_endian") {
        format = ply_format::binary_little_endian;
      } else if (fmt == "binary_big_endian") {
        format = ply_format::binary_big_endian;
      } else {
        return {filename + ": unknown format"};
      }
    } else if (cmd == "comment") {
      skip_ply_whitespace(str);
      comments.push_back(string{str});
    } else if (cmd == "obj_info") {
      skip_ply_whitespace(str);
      // comment is the rest of the str
    } else if (cmd == "element") {
      auto& elem = elements.emplace_back();
      if (!parse_ply_value(str, elem.name)) return {filename + ": parse error"};
      if (!parse_ply_value(str, elem.count))
        return {filename + ": parse error"};
    } else if (cmd == "property") {
      if (elements.empty()) throw std::runtime_error{"bad ply header"};
      auto& prop  = elements.back().properties.emplace_back();
      auto  tname = ""s;
      if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
      if (tname == "list") {
        prop.is_list = true;
        if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
        if (type_map.find(tname) == type_map.end())
          return {filename + ": unknown type " + tname};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          throw std::runtime_error{"unsupported list size type " + tname};
        if (!parse_ply_value(str, tname)) return {filename + ": parse error"};
        if (type_map.find(tname) == type_map.end())
          return {filename + ": unknown type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          return {filename + ": unknown type " + tname};
        prop.type = type_map.at(tname);
      }
      if (!parse_ply_value(str, prop.name)) return {filename + ": parse error"};
    } else if (cmd == "end_header") {
      end_header = true;
      break;
    } else {
      return {filename + ": unknown command " + cmd};
    }
  }

  if (!end_header) return {filename + ": incomplete header"};

  return {};
}

template <typename VT, typename LT>
static plyio_status read_ply_value_generic(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<VT>& values,
    vector<vector<LT>>& lists) {
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
    if (!read_ply_line(fs, buffer, sizeof(buffer)))
      return {filename + ": read error"};
    auto str = string_view{buffer};
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop  = element.properties[pidx];
      auto& value = values[pidx];
      auto& list  = lists[pidx];
      if (!prop.is_list) {
        if (!parse_ply_prop(str, prop.type, value))
          return {filename + ": parse error"};
      } else {
        if (!parse_ply_prop(str, ply_type::u8, value))
          return {filename + ": parse error"};
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          if (!parse_ply_prop(str, prop.type, list[i]))
            return {filename + ": parse error"};
      }
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop  = element.properties[pidx];
      auto& value = values[pidx];
      auto& list  = lists[pidx];
      if (!prop.is_list) {
        if (!read_ply_prop(
                fs, prop.type, value, format == ply_format::binary_big_endian))
          return {filename + ": read error"};
      } else {
        if (!read_ply_prop(fs, ply_type::u8, value,
                format == ply_format::binary_big_endian))
          return {filename + ": read error"};
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          if (!read_ply_prop(fs, prop.type, list[i],
                  format == ply_format::binary_big_endian))
            return {filename + ": read error"};
      }
    }
  }

  return {};
}

template <typename VT>
static bool format_ply_prop(ply_file& fs, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: return format_ply_value(fs, (int8_t)value);
    case ply_type::i16: return format_ply_value(fs, (int16_t)value);
    case ply_type::i32: return format_ply_value(fs, (int32_t)value);
    case ply_type::i64: return format_ply_value(fs, (int64_t)value);
    case ply_type::u8: return format_ply_value(fs, (uint8_t)value);
    case ply_type::u16: return format_ply_value(fs, (uint16_t)value);
    case ply_type::u32: return format_ply_value(fs, (uint32_t)value);
    case ply_type::u64: return format_ply_value(fs, (uint64_t)value);
    case ply_type::f32: return format_ply_value(fs, (float)value);
    case ply_type::f64: return format_ply_value(fs, (double)value);
  }
}

template <typename VT>
static bool write_ply_prop(
    ply_file& fs, ply_type type, VT value, bool big_endian) {
  switch (type) {
    case ply_type::i8: return write_ply_value(fs, (int8_t)value, big_endian);
    case ply_type::i16: return write_ply_value(fs, (int16_t)value, big_endian);
    case ply_type::i32: return write_ply_value(fs, (int32_t)value, big_endian);
    case ply_type::i64: return write_ply_value(fs, (int64_t)value, big_endian);
    case ply_type::u8: return write_ply_value(fs, (uint8_t)value, big_endian);
    case ply_type::u16: return write_ply_value(fs, (uint16_t)value, big_endian);
    case ply_type::u32: return write_ply_value(fs, (uint32_t)value, big_endian);
    case ply_type::u64: return write_ply_value(fs, (uint64_t)value, big_endian);
    case ply_type::f32: return write_ply_value(fs, (float)value, big_endian);
    case ply_type::f64: return write_ply_value(fs, (double)value, big_endian);
    default: return false;
  }
}

// Write Ply functions
plyio_status write_ply_header(const string& filename, ply_file& fs,
    ply_format format, const vector<ply_element>& elements,
    const vector<string>& comments) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{{ply_type::i8, "char"},
      {ply_type::i16, "short"}, {ply_type::i32, "int"}, {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"}, {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"}, {ply_type::u64, "ulong"},
      {ply_type::f32, "float"}, {ply_type::f64, "double"}};

  if (!format_ply_values(fs, "ply\n")) return {filename + ": write error"};
  switch (format) {
    case ply_format::ascii:
      if (!format_ply_values(fs, "format ascii 1.0\n"))
        return {filename + ": write error"};
      break;
    case ply_format::binary_little_endian:
      if (!format_ply_values(fs, "format binary_little_endian 1.0\n"))
        return {filename + ": write error"};
      break;
    case ply_format::binary_big_endian:
      if (!format_ply_values(fs, "format binary_big_endian 1.0\n"))
        return {filename + ": write error"};
      break;
  }
  for (auto& comment : comments)
    if (!format_ply_values(fs, "comment " + comment + "\n"))
      return {filename + ": write error"};
  for (auto& elem : elements) {
    if (!format_ply_values(fs,
            "element " + elem.name + " " + std::to_string(elem.count) + "\n"))
      return {filename + ": write error"};
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        if (!format_ply_values(fs, "property list uchar " +
                                       type_map[prop.type] + " " + prop.name +
                                       "\n"))
          return {filename + ": write error"};
      } else {
        if (!format_ply_values(
                fs, "property " + type_map[prop.type] + " " + prop.name + "\n"))
          return {filename + ": write error"};
      }
    }
  }
  if (!format_ply_values(fs, "end_header\n"))
    return {filename + ": write error"};

  return {};
}

template <typename VT, typename LT>
static plyio_status write_ply_value_generic(const string& filename,
    ply_file& fs, ply_format format, const ply_element& element,
    vector<VT>& values, vector<vector<LT>>& lists) {
  if (format == ply_format::ascii) {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (pidx)
        if (!format_ply_value(fs, " ")) return {filename + ": write error"};
      if (!prop.is_list) {
        if (!format_ply_prop(fs, prop.type, values[pidx]))
          return {filename + ": write error"};
      } else {
        if (!format_ply_prop(fs, ply_type::u8, values[pidx]))
          return {filename + ": write error"};
        for (auto i = 0; i < (int)lists[pidx].size(); i++) {
          if (i)
            if (!format_ply_value(fs, " ")) return {filename + ": write error"};
          if (!format_ply_prop(fs, prop.type, lists[pidx][i]))
            return {filename + ": write error"};
        }
      }
      if (!format_ply_value(fs, "\n")) return {filename + ": write error"};
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (!prop.is_list) {
        if (!write_ply_prop(fs, prop.type, values[pidx],
                format == ply_format::binary_big_endian))
          return {filename + ": write error"};
      } else {
        if (!write_ply_prop(fs, ply_type::u8, values[pidx],
                format == ply_format::binary_big_endian))
          return {filename + ": write error"};
        for (auto i = 0; i < (int)lists[pidx].size(); i++)
          if (!write_ply_prop(fs, prop.type, lists[pidx][i],
                  format == ply_format::binary_big_endian))
            return {filename + ": write error"};
      }
    }
  }
  return {};
}

plyio_status write_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  return write_ply_value_generic(filename, fs, format, element, values, lists);
}
plyio_status write_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  return write_ply_value_generic(filename, fs, format, element, values, lists);
}

plyio_status read_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  return read_ply_value_generic(filename, fs, format, element, values, lists);
}
plyio_status read_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  return read_ply_value_generic(filename, fs, format, element, values, lists);
}

int find_ply_element(
    const vector<ply_element>& elements, const string& name) {
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
// LOW-LEVEL FILE HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
obj_file::obj_file(obj_file&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
obj_file::~obj_file() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
void open_obj(obj_file& fs, const string& filename, const string& mode) {
  close_obj(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
obj_file open_obj(const string& filename, const string& mode) {
  auto fs = obj_file{};
  open_obj(fs, filename, mode);
  return fs;
}
void close_obj(obj_file& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

// Read a line
static bool read_obj_line(obj_file& fs, char* buffer, size_t size) {
  if (fgets(buffer, size, fs.fs)) {
    fs.linenum += 1;
    return true;
  } else {
    return false;
  }
}

// Check for errors
static bool has_obj_error(obj_file& fs) { return ferror(fs.fs); }

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// Utility to normalize a path
static string normalize_obj_path(const string& filename_) {
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
static string get_obj_dirname(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
static string get_obj_extension(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos);
}

// Get filename without directory.
static string get_obj_filename(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
static string get_obj_noextension(const string& filename_) {
  auto filename = normalize_obj_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
static string get_obj_basename(const string& filename) {
  return get_obj_noextension(get_obj_filename(filename));
}

// Replaces extensions
static string replace_obj_extension(const string& filename, const string& ext) {
  return get_obj_noextension(filename) + ext;
}

// Check if a file can be opened for reading.
static bool exists_obj_file(const string& filename) {
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
static bool is_obj_newline(char c) { return c == '\r' || c == '\n'; }
static bool is_obj_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static void skip_obj_whitespace(string_view& str) {
  while (!str.empty() && is_obj_space(str.front())) str.remove_prefix(1);
}

static void remove_obj_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_obj_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
static bool parse_obj_value(string_view& str, string_view& value) {
  skip_obj_whitespace(str);
  if (str.empty()) return false;
  auto cpy = str;
  while (!cpy.empty() && !is_obj_space(cpy.front())) cpy.remove_prefix(1);
  value = str;
  value.remove_suffix(cpy.size());
  str.remove_prefix(str.size() - cpy.size());
  return true;
}
static bool parse_obj_value(string_view& str, string& value) {
  auto valuev = string_view{};
  if (!parse_obj_value(str, valuev)) return false;
  value = string{valuev};
  return true;
}
static bool parse_obj_value(string_view& str, int& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}
static bool parse_obj_value(string_view& str, bool& value) {
  auto valuei = 0;
  if (!parse_obj_value(str, valuei)) return false;
  value = (bool)valuei;
  return true;
}
static bool parse_obj_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str.data() == end) return false;
  str.remove_prefix(end - str.data());
  return true;
}

static bool parse_obj_value(string_view& str, vec2f& value) {
  for (auto i = 0; i < 2; i++)
    if (!parse_obj_value(str, value[i])) return false;
  return true;
}
static bool parse_obj_value(string_view& str, vec3f& value) {
  for (auto i = 0; i < 3; i++)
    if (!parse_obj_value(str, value[i])) return false;
  return true;
}
static bool parse_obj_value(string_view& str, frame3f& value) {
  for (auto i = 0; i < 4; i++)
    if (!parse_obj_value(str, value[i])) return false;
  return true;
}

// Parse values from a string
static bool parse_obj_value_or_empty(string_view& str, string& value) {
  skip_obj_whitespace(str);
  if (str.empty()) {
    value = "";
    return true;
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
static void format_obj_value(string& str, const string& value) { str += value; }
static void format_obj_value(string& str, const char* value) { str += value; }
static void format_obj_value(string& str, int value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
static void format_obj_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
static void format_obj_value(string& str, const vec2f& value) {
  for (auto i = 0; i < 2; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
}
static void format_obj_value(string& str, const vec3f& value) {
  for (auto i = 0; i < 3; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
}
static void format_obj_value(string& str, const frame3f& value) {
  for (auto i = 0; i < 4; i++) {
    if (i) str += " ";
    format_obj_value(str, value[i]);
  }
}

// Foramt to file
static void format_obj_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
static void format_obj_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_obj_value(str, arg);
  format_obj_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
static bool format_obj_values(
    obj_file& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_obj_values(str, fmt, args...);
  return fputs(str.c_str(), fs.fs) >= 0;
}
template <typename T>
static bool format_obj_value(obj_file& fs, const T& value) {
  auto str = ""s;
  format_obj_value(str, value);
  return fputs(str.c_str(), fs.fs) >= 0;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static bool parse_obj_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  if (!parse_obj_value(str, value.position)) return false;
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      if (!parse_obj_value(str, value.normal)) return false;
    } else {
      if (!parse_obj_value(str, value.texcoord)) return false;
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        if (!parse_obj_value(str, value.normal)) return false;
      }
    }
  }
  return true;
}

// Input for OBJ textures
static bool parse_obj_value(string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_obj_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    if (!parse_obj_value(str, token)) return false;
    tokens.push_back(token);
    skip_obj_whitespace(str);
  }
  if (tokens.empty()) return false;

  // texture name
  info.path = normalize_obj_path(tokens.back());

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }

  return true;
}

// Read obj
static objio_status load_mtl(
    const string& filename, obj_model& obj, bool fliptr = true) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) return {filename + ": file not found"};

  // init parsing
  obj.materials.emplace_back();

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      obj.materials.emplace_back();
      if (!parse_obj_value(str, obj.materials.back().name))
        return {filename + ": parse error"};
    } else if (cmd == "illum") {
      if (!parse_obj_value(str, obj.materials.back().illum))
        return {filename + ": parse error"};
    } else if (cmd == "Ke") {
      if (!parse_obj_value(str, obj.materials.back().emission))
        return {filename + ": parse error"};
    } else if (cmd == "Ka") {
      if (!parse_obj_value(str, obj.materials.back().ambient))
        return {filename + ": parse error"};
    } else if (cmd == "Kd") {
      if (!parse_obj_value(str, obj.materials.back().diffuse))
        return {filename + ": parse error"};
    } else if (cmd == "Ks") {
      if (!parse_obj_value(str, obj.materials.back().specular))
        return {filename + ": parse error"};
    } else if (cmd == "Kt") {
      if (!parse_obj_value(str, obj.materials.back().transmission))
        return {filename + ": parse error"};
    } else if (cmd == "Tf") {
      obj.materials.back().transmission = vec3f{-1};
      if (!parse_obj_value(str, obj.materials.back().transmission))
        return {filename + ": parse error"};
      if (obj.materials.back().transmission.y < 0)
        obj.materials.back().transmission = vec3f{
            obj.materials.back().transmission.x};
      if (fliptr)
        obj.materials.back().transmission = 1 -
                                            obj.materials.back().transmission;
    } else if (cmd == "Tr") {
      if (!parse_obj_value(str, obj.materials.back().opacity))
        return {filename + ": parse error"};
      if (fliptr)
        obj.materials.back().opacity = 1 - obj.materials.back().opacity;
    } else if (cmd == "Ns") {
      if (!parse_obj_value(str, obj.materials.back().exponent))
        return {filename + ": parse error"};
    } else if (cmd == "d") {
      if (!parse_obj_value(str, obj.materials.back().opacity))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ke") {
      if (!parse_obj_value(str, obj.materials.back().emission_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ka") {
      if (!parse_obj_value(str, obj.materials.back().ambient_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Kd") {
      if (!parse_obj_value(str, obj.materials.back().diffuse_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ks") {
      if (!parse_obj_value(str, obj.materials.back().specular_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Tr") {
      if (!parse_obj_value(str, obj.materials.back().transmission_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      if (!parse_obj_value(str, obj.materials.back().opacity_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_bump" || cmd == "bump") {
      if (!parse_obj_value(str, obj.materials.back().bump_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_disp" || cmd == "disp") {
      if (!parse_obj_value(str, obj.materials.back().displacement_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_norm" || cmd == "norm") {
      if (!parse_obj_value(str, obj.materials.back().normal_map))
        return {filename + ": parse error"};
    } else if (cmd == "Pm") {
      if (!parse_obj_value(str, obj.materials.back().pbr_metallic))
        return {filename + ": parse error"};
    } else if (cmd == "Pr") {
      if (!parse_obj_value(str, obj.materials.back().pbr_roughness))
        return {filename + ": parse error"};
    } else if (cmd == "Ps") {
      if (!parse_obj_value(str, obj.materials.back().pbr_sheen))
        return {filename + ": parse error"};
    } else if (cmd == "Pc") {
      if (!parse_obj_value(str, obj.materials.back().pbr_clearcoat))
        return {filename + ": parse error"};
    } else if (cmd == "Pcr") {
      if (!parse_obj_value(str, obj.materials.back().pbr_coatroughness))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pm") {
      if (!parse_obj_value(str, obj.materials.back().pbr_metallic_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pr") {
      if (!parse_obj_value(str, obj.materials.back().pbr_roughness_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ps") {
      if (!parse_obj_value(str, obj.materials.back().pbr_sheen_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pc") {
      if (!parse_obj_value(str, obj.materials.back().pbr_clearcoat_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pcr") {
      if (!parse_obj_value(str, obj.materials.back().pbr_coatroughness_map))
        return {filename + ": parse error"};
    } else if (cmd == "Vt") {
      if (!parse_obj_value(str, obj.materials.back().vol_transmission))
        return {filename + ": parse error"};
    } else if (cmd == "Vp") {
      if (!parse_obj_value(str, obj.materials.back().vol_meanfreepath))
        return {filename + ": parse error"};
    } else if (cmd == "Ve") {
      if (!parse_obj_value(str, obj.materials.back().vol_emission))
        return {filename + ": parse error"};
    } else if (cmd == "Vs") {
      if (!parse_obj_value(str, obj.materials.back().vol_scattering))
        return {filename + ": parse error"};
    } else if (cmd == "Vg") {
      if (!parse_obj_value(str, obj.materials.back().vol_anisotropy))
        return {filename + ": parse error"};
    } else if (cmd == "Vr") {
      if (!parse_obj_value(str, obj.materials.back().vol_scale))
        return {filename + ": parse error"};
    } else if (cmd == "map_Vs") {
      if (!parse_obj_value(str, obj.materials.back().vol_scattering_map))
        return {filename + ": parse error"};
    } else {
      continue;
    }
  }

  // check error
  if (has_obj_error(fs)) return {filename + ": read error"};

  // remove placeholder material
  obj.materials.erase(obj.materials.begin());

  return {};
}

// Read obj
static objio_status load_objx(const string& filename, obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) return {filename + ": file not found"};

  // shape map for instances
  auto shape_map = unordered_map<string, vector<int>>{};
  for (auto idx = 0; idx < obj.shapes.size(); idx++) {
    shape_map[obj.shapes[idx].name].push_back(idx);
  }

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      auto& camera = obj.cameras.emplace_back();
      if (!parse_obj_value(str, camera.name))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.ortho))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.width))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.height))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.lens))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.focus))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.aperture))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.frame))
        return {filename + ": parse error"};
    } else if (cmd == "e") {
      auto& environment = obj.environments.emplace_back();
      if (!parse_obj_value(str, environment.name))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, environment.emission))
        return {filename + ": parse error"};
      auto emission_path = ""s;
      if (!parse_obj_value(str, emission_path))
        return {filename + ": parse error"};
      if (emission_path == "\"\"") emission_path = "";
      environment.emission_map.path = emission_path;
      if (!parse_obj_value(str, environment.frame))
        return {filename + ": parse error"};
    } else if (cmd == "i") {
      auto object = ""s;
      auto frame  = identity3x4f;
      if (!parse_obj_value(str, object)) return {filename + ": parse error"};
      if (!parse_obj_value(str, frame)) return {filename + ": parse error"};
      if (shape_map.find(object) == shape_map.end()) {
        return {filename + ": missing object " + object};
      }
      for (auto idx : shape_map.at(object)) {
        obj.shapes[idx].instances.push_back(frame);
      }
    } else {
      // unused
    }
  }

  // check error
  if (has_obj_error(fs)) return {filename + "read error"};

  return {};
}

// Read obj
objio_status load_obj(const string& filename, obj_model& obj,
    bool geom_only, bool split_elements, bool split_materials) {
  // open file
  auto fs = open_obj(filename, "rt");
  if (!fs) return {filename + ": file not found"};

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

  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      if (!parse_obj_value(str, opositions.emplace_back()))
        return {filename + ": parse error"};
      vert_size.position += 1;
    } else if (cmd == "vn") {
      if (!parse_obj_value(str, onormals.emplace_back()))
        return {filename + ": parse error"};
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      if (!parse_obj_value(str, otexcoords.emplace_back()))
        return {filename + ": parse error"};
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
      skip_obj_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        if (!parse_obj_value(str, vert)) return {filename + ": parse error"};
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        shape.vertices.push_back(vert);
        element.size += 1;
        skip_obj_whitespace(str);
      }
    } else if (cmd == "o" || cmd == "g") {
      if (geom_only) continue;
      if (!parse_obj_value_or_empty(str, cmd == "o" ? oname : gname))
        return {filename + ": parse error"};
      if (!obj.shapes.back().vertices.empty()) {
        obj.shapes.emplace_back();
        obj.shapes.back().name = oname + gname;
      } else {
        obj.shapes.back().name = oname + gname;
      }
    } else if (cmd == "usemtl") {
      if (geom_only) continue;
      if (!parse_obj_value_or_empty(str, mname))
        return {filename + ": parse error"};
    } else if (cmd == "s") {
      if (geom_only) continue;
    } else if (cmd == "mtllib") {
      if (geom_only) continue;
      auto mtllib = ""s;
      if (!parse_obj_value(str, mtllib)) return {filename + ": parse error"};
      if (std::find(mtllibs.begin(), mtllibs.end(), mtllib) == mtllibs.end()) {
        mtllibs.push_back(mtllib);
      }
    } else {
      // unused
    }
  }

  // check error
  if (has_obj_error(fs)) return {filename + ": read error"};

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
  if (geom_only) return {};

  // load materials
  auto dirname = get_obj_dirname(filename);
  for (auto& mtllib : mtllibs) {
    if (auto ret = load_mtl(dirname + mtllib, obj); !ret)
      return {filename + ": mtl error (" + ret.error + ")"};
  }

  // load extensions
  auto extfilename = replace_obj_extension(filename, ".objx");
  if (exists_obj_file(extfilename)) {
    if (auto ret = load_objx(extfilename, obj); !ret)
      return {filename + ": objx error (" + ret.error + ")"};
  }

  return {};
}

// Format values
static void format_obj_value(string& str, const obj_texture_info& value) {
  str += value.path.empty() ? "" : value.path;
}
static void format_obj_value(string& str, const obj_vertex& value) {
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
static objio_status save_mtl(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");
  if (!fs) return {filename + ": file not found"};

  // save comments
  if (!format_obj_values(fs, "#\n")) return {filename + ": write error"};
  if (!format_obj_values(fs, "# Written by Yocto/GL\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "#\n\n")) return {filename + ": write error"};
  for (auto& comment : obj.comments) {
    if (!format_obj_values(fs, "# {}\n", comment))
      return {filename + ": write error"};
  }
  if (!format_obj_values(fs, "\n")) return {filename + ": write error"};

  // write material
  for (auto& material : obj.materials) {
    if (!format_obj_values(fs, "newmtl {}\n", material.name))
      return {filename + ": write error"};
    if (!format_obj_values(fs, "illum {}\n", material.illum))
      return {filename + ": write error"};
    if (material.emission != zero3f)
      if (!format_obj_values(fs, "Ke {}\n", material.emission))
        return {filename + ": write error"};
    if (material.ambient != zero3f)
      if (!format_obj_values(fs, "Ka {}\n", material.ambient))
        return {filename + ": write error"};
    if (!format_obj_values(fs, "Kd {}\n", material.diffuse))
      return {filename + ": write error"};
    if (!format_obj_values(fs, "Ks {}\n", material.specular))
      return {filename + ": write error"};
    if (material.reflection != zero3f)
      if (!format_obj_values(fs, "Kr {}\n", material.reflection))
        return {filename + ": write error"};
    if (material.transmission != zero3f)
      if (!format_obj_values(fs, "Kt {}\n", material.transmission))
        return {filename + ": write error"};
    format_obj_values(fs, "Ns {}\n", (int)material.exponent);
    if (material.opacity != 1)
      if (!format_obj_values(fs, "d {}\n", material.opacity))
        return {filename + ": write error"};
    if (!material.emission_map.path.empty())
      if (!format_obj_values(fs, "map_Ke {}\n", material.emission_map))
        return {filename + ": write error"};
    if (!material.diffuse_map.path.empty())
      if (!format_obj_values(fs, "map_Kd {}\n", material.diffuse_map))
        return {filename + ": write error"};
    if (!material.specular_map.path.empty())
      if (!format_obj_values(fs, "map_Ks {}\n", material.specular_map))
        return {filename + ": write error"};
    if (!material.transmission_map.path.empty())
      if (!format_obj_values(fs, "map_Kt {}\n", material.transmission_map))
        return {filename + ": write error"};
    if (!material.reflection_map.path.empty())
      if (!format_obj_values(fs, "map_Kr {}\n", material.reflection_map))
        return {filename + ": write error"};
    if (!material.exponent_map.path.empty())
      if (!format_obj_values(fs, "map_Ns {}\n", material.exponent_map))
        return {filename + ": write error"};
    if (!material.opacity_map.path.empty())
      if (!format_obj_values(fs, "map_d {}\n", material.opacity_map))
        return {filename + ": write error"};
    if (!material.bump_map.path.empty())
      if (!format_obj_values(fs, "map_bump {}\n", material.bump_map))
        return {filename + ": write error"};
    if (!material.displacement_map.path.empty())
      if (!format_obj_values(fs, "map_disp {}\n", material.displacement_map))
        return {filename + ": write error"};
    if (!material.normal_map.path.empty())
      if (!format_obj_values(fs, "map_norm {}\n", material.normal_map))
        return {filename + ": write error"};
    if (material.pbr_roughness)
      if (!format_obj_values(fs, "Pr {}\n", material.pbr_roughness))
        return {filename + ": write error"};
    if (material.pbr_metallic)
      if (!format_obj_values(fs, "Pm {}\n", material.pbr_metallic))
        return {filename + ": write error"};
    if (material.pbr_sheen)
      if (!format_obj_values(fs, "Ps {}\n", material.pbr_sheen))
        return {filename + ": write error"};
    if (material.pbr_clearcoat)
      if (!format_obj_values(fs, "Pc {}\n", material.pbr_clearcoat))
        return {filename + ": write error"};
    if (material.pbr_coatroughness)
      if (!format_obj_values(fs, "Pcr {}\n", material.pbr_coatroughness))
        return {filename + ": write error"};
    if (!material.pbr_roughness_map.path.empty())
      if (!format_obj_values(fs, "map_Pr {}\n", material.pbr_roughness_map))
        return {filename + ": write error"};
    if (!material.pbr_metallic_map.path.empty())
      if (!format_obj_values(fs, "map_Pm {}\n", material.pbr_metallic_map))
        return {filename + ": write error"};
    if (!material.pbr_sheen_map.path.empty())
      if (!format_obj_values(fs, "map_Ps {}\n", material.pbr_sheen_map))
        return {filename + ": write error"};
    if (!material.pbr_clearcoat_map.path.empty())
      if (!format_obj_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map))
        return {filename + ": write error"};
    if (!material.pbr_coatroughness_map.path.empty())
      if (!format_obj_values(
              fs, "map_Pcr {}\n", material.pbr_coatroughness_map))
        return {filename + ": write error"};
    if (material.vol_transmission != zero3f)
      if (!format_obj_values(fs, "Vt {}\n", material.vol_transmission))
        return {filename + ": write error"};
    if (material.vol_meanfreepath != zero3f)
      if (!format_obj_values(fs, "Vp {}\n", material.vol_meanfreepath))
        return {filename + ": write error"};
    if (material.vol_emission != zero3f)
      if (!format_obj_values(fs, "Ve {}\n", material.vol_emission))
        return {filename + ": write error"};
    if (material.vol_scattering != zero3f)
      if (!format_obj_values(fs, "Vs {}\n", material.vol_scattering))
        return {filename + ": write error"};
    if (material.vol_anisotropy)
      if (!format_obj_values(fs, "Vg {}\n", material.vol_anisotropy))
        return {filename + ": write error"};
    if (material.vol_scale)
      if (!format_obj_values(fs, "Vr {}\n", material.vol_scale))
        return {filename + ": write error"};
    if (!material.vol_scattering_map.path.empty())
      if (!format_obj_values(fs, "map_Vs {}\n", material.vol_scattering_map))
        return {filename + ": write error"};
    if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
  }

  return {};
}

// Save obj
static objio_status save_objx(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");
  if (!fs) return {filename + ": file not found"};

  // save comments
  if (!format_obj_values(fs, "#\n")) return {filename + ": write error"};
  if (!format_obj_values(fs, "# Written by Yocto/GL\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "#\n\n")) return {filename + ": write error"};
  for (auto& comment : obj.comments) {
    if (!format_obj_values(fs, "# {}\n", comment))
      return {filename + ": write error"};
  }
  if (!format_obj_values(fs, "\n")) return {filename + ": write error"};

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

  return {};
}

// Save obj
objio_status save_obj(const string& filename, const obj_model& obj) {
  // open file
  auto fs = open_obj(filename, "wt");
  if (!fs) return {filename + ": file not found"};

  // save comments
  if (!format_obj_values(fs, "#\n")) return {filename + ": write error"};
  if (!format_obj_values(fs, "# Written by Yocto/GL\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "# https://github.com/xelatihy/yocto-gl\n"))
    return {filename + ": write error"};
  if (!format_obj_values(fs, "#\n\n")) return {filename + ": write error"};
  for (auto& comment : obj.comments) {
    if (!format_obj_values(fs, "# {}\n", comment))
      return {filename + ": write error"};
  }
  if (!format_obj_values(fs, "\n")) return {filename + ": write error"};

  // save material library
  if (!obj.materials.empty()) {
    format_obj_values(fs, "mtllib {}\n\n",
        replace_obj_extension(get_obj_filename(filename), ".mtl"));
  }

  // save objects
  auto vert_size = obj_vertex{0, 0, 0};
  for (auto& shape : obj.shapes) {
    if (!format_obj_values(fs, "o {}\n", shape.name))
      return {filename + ": write error"};
    for (auto& p : shape.positions)
      if (!format_obj_values(fs, "v {}\n", p))
        return {filename + ": write error"};
    for (auto& n : shape.normals)
      if (!format_obj_values(fs, "vn {}\n", n))
        return {filename + ": write error"};
    for (auto& t : shape.texcoords)
      if (!format_obj_values(fs, "vt {}\n", t))
        return {filename + ": write error"};
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
        if (!format_obj_values(fs, "{}", label))
          return {filename + ": write error"};
        for (auto c = 0; c < element.size; c++) {
          auto vert = shape.vertices[cur_vertex++];
          if (vert.position) vert.position += vert_size.position;
          if (vert.normal) vert.normal += vert_size.normal;
          if (vert.texcoord) vert.texcoord += vert_size.texcoord;
          if (!format_obj_values(fs, " {}", vert))
            return {filename + ": write error"};
        }
        if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
      }
    }
    if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
    vert_size.position += (int)shape.positions.size();
    vert_size.normal += (int)shape.normals.size();
    vert_size.texcoord += (int)shape.texcoords.size();
  }

  // save mtl
  if (!obj.materials.empty()) {
    if (auto ret = save_mtl(replace_obj_extension(filename, ".mtl"), obj); !ret)
      return ret;
  }

  // save objx
  if (!obj.cameras.empty() || !obj.environments.empty() ||
      std::any_of(obj.shapes.begin(), obj.shapes.end(),
          [](auto& shape) { return !shape.instances.empty(); })) {
    if (auto ret = save_objx(replace_obj_extension(filename, ".objx"), obj);
        !ret)
      return ret;
  }

  return {};
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
static void get_obj_vertices(const obj_shape& shape, vector<vec3f>& positions,
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
static vector<vec2f> flip_obj_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
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
  for (auto& str : shape.lines) {
    for (auto c = 1; c < str.size; c++) {
      lines.push_back({vindex[cur + c - 1], vindex[cur + c]});
      if (!materials.empty()) ematerials.push_back(str.material);
    }
    cur += str.size;
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

bool has_obj_quads(const obj_shape& shape) {
  for (auto& face : shape.faces)
    if (face.size == 4) return true;
  return false;
}

// Add obj shape
void add_obj_triangles(obj_model& obj, const string& name,
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
void add_obj_quads(obj_model& obj, const string& name,
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
void add_obj_lines(obj_model& obj, const string& name,
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
void add_obj_points(obj_model& obj, const string& name,
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
void add_obj_fvquads(obj_model& obj, const string& name,
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
objio_status read_obj_command(const string& filename, obj_file& fs,
    obj_command& command, string& name, vec3f& value,
    vector<obj_vertex>& vertices, obj_vertex& vert_size) {
  // read the file str by str
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      if (!parse_obj_value(str, value)) return {filename + ": parse error"};
      vert_size.position += 1;
      return {};
    } else if (cmd == "vn") {
      command = obj_command::normal;
      if (!parse_obj_value(str, value)) return {filename + ": parse error"};
      vert_size.normal += 1;
      return {};
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
      if (!parse_obj_value(str, (vec2f&)value))
        return {filename + ": parse error"};
      value.z = 0;
      vert_size.texcoord += 1;
      return {};
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      skip_obj_whitespace(str);
      while (!str.empty()) {
        auto vert = obj_vertex{};
        if (!parse_obj_value(str, vert)) return {filename + ": parse error"};
        if (!vert.position) break;
        if (vert.position < 0)
          vert.position = vert_size.position + vert.position + 1;
        if (vert.texcoord < 0)
          vert.texcoord = vert_size.texcoord + vert.texcoord + 1;
        if (vert.normal < 0) vert.normal = vert_size.normal + vert.normal + 1;
        vertices.push_back(vert);
        skip_obj_whitespace(str);
      }
      if (cmd == "f") command = obj_command::face;
      if (cmd == "l") command = obj_command::str;
      if (cmd == "p") command = obj_command::point;
      return {};
    } else if (cmd == "o") {
      command = obj_command::object;
      if (!parse_obj_value_or_empty(str, name))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      if (!parse_obj_value_or_empty(str, name))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "g") {
      command = obj_command::group;
      if (!parse_obj_value_or_empty(str, name))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      if (!parse_obj_value_or_empty(str, name))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      if (!parse_obj_value(str, name)) return {filename + ": parse error"};
      return {};
    } else {
      // unused
    }
  }

  // check error
  if (has_obj_error(fs)) return {filename + ": read error"};

  return {"eof"};
}

// Read mtl
objio_status read_mtl_command(const string& filename, obj_file& fs,
    mtl_command& command, obj_material& material, bool fliptr) {
  material = {};

  // read the file str by str
  auto pos   = ftell(fs.fs);
  auto found = false;
  char buffer[4096];
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      if (found) {
        command = mtl_command::material;
        fseek(fs.fs, pos, SEEK_SET);
        return {};
      } else {
        found = true;
      }
      if (!parse_obj_value(str, material.name))
        return {filename + ": parse error"};
    } else if (cmd == "illum") {
      if (!parse_obj_value(str, material.illum))
        return {filename + ": parse error"};
    } else if (cmd == "Ke") {
      if (!parse_obj_value(str, material.emission))
        return {filename + ": parse error"};
    } else if (cmd == "Ka") {
      if (!parse_obj_value(str, material.ambient))
        return {filename + ": parse error"};
    } else if (cmd == "Kd") {
      if (!parse_obj_value(str, material.diffuse))
        return {filename + ": parse error"};
    } else if (cmd == "Ks") {
      if (!parse_obj_value(str, material.specular))
        return {filename + ": parse error"};
    } else if (cmd == "Kt") {
      if (!parse_obj_value(str, material.transmission))
        return {filename + ": parse error"};
    } else if (cmd == "Tf") {
      material.transmission = vec3f{-1};
      if (!parse_obj_value(str, material.transmission))
        return {filename + ": parse error"};
      if (material.transmission.y < 0)
        material.transmission = vec3f{material.transmission.x};
      if (fliptr) material.transmission = 1 - material.transmission;
    } else if (cmd == "Tr") {
      if (!parse_obj_value(str, material.opacity))
        return {filename + ": parse error"};
      if (fliptr) material.opacity = 1 - material.opacity;
    } else if (cmd == "Ns") {
      if (!parse_obj_value(str, material.exponent))
        return {filename + ": parse error"};
    } else if (cmd == "d") {
      if (!parse_obj_value(str, material.opacity))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ke") {
      if (!parse_obj_value(str, material.emission_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ka") {
      if (!parse_obj_value(str, material.ambient_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Kd") {
      if (!parse_obj_value(str, material.diffuse_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ks") {
      if (!parse_obj_value(str, material.specular_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Tr") {
      if (!parse_obj_value(str, material.transmission_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      if (!parse_obj_value(str, material.opacity_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_bump" || cmd == "bump") {
      if (!parse_obj_value(str, material.bump_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_disp" || cmd == "disp") {
      if (!parse_obj_value(str, material.displacement_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_norm" || cmd == "norm") {
      if (!parse_obj_value(str, material.normal_map))
        return {filename + ": parse error"};
    } else if (cmd == "Pm") {
      if (!parse_obj_value(str, material.pbr_metallic))
        return {filename + ": parse error"};
    } else if (cmd == "Pr") {
      if (!parse_obj_value(str, material.pbr_roughness))
        return {filename + ": parse error"};
    } else if (cmd == "Ps") {
      if (!parse_obj_value(str, material.pbr_sheen))
        return {filename + ": parse error"};
    } else if (cmd == "Pc") {
      if (!parse_obj_value(str, material.pbr_clearcoat))
        return {filename + ": parse error"};
    } else if (cmd == "Pcr") {
      if (!parse_obj_value(str, material.pbr_coatroughness))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pm") {
      if (!parse_obj_value(str, material.pbr_metallic_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pr") {
      if (!parse_obj_value(str, material.pbr_roughness_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Ps") {
      if (!parse_obj_value(str, material.pbr_sheen_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pc") {
      if (!parse_obj_value(str, material.pbr_clearcoat_map))
        return {filename + ": parse error"};
    } else if (cmd == "map_Pcr") {
      if (!parse_obj_value(str, material.pbr_coatroughness_map))
        return {filename + ": parse error"};
    } else if (cmd == "Vt") {
      if (!parse_obj_value(str, material.vol_transmission))
        return {filename + ": parse error"};
    } else if (cmd == "Vp") {
      if (!parse_obj_value(str, material.vol_meanfreepath))
        return {filename + ": parse error"};
    } else if (cmd == "Ve") {
      if (!parse_obj_value(str, material.vol_emission))
        return {filename + ": parse error"};
    } else if (cmd == "Vs") {
      if (!parse_obj_value(str, material.vol_scattering))
        return {filename + ": parse error"};
    } else if (cmd == "Vg") {
      if (!parse_obj_value(str, material.vol_anisotropy))
        return {filename + ": parse error"};
    } else if (cmd == "Vr") {
      if (!parse_obj_value(str, material.vol_scale))
        return {filename + ": parse error"};
    } else if (cmd == "map_Vs") {
      if (!parse_obj_value(str, material.vol_scattering_map))
        return {filename + ": parse error"};
    } else {
      continue;
    }
    pos = ftell(fs.fs);
  }

  if (found) {
    command = mtl_command::material;
    return {};
  }

  // check error
  if (has_obj_error(fs)) return {filename + "read error"};

  return {"eof"};
}

// Read objx
objio_status read_objx_command(const string& filename, obj_file& fs,
    objx_command& command, obj_camera& camera, obj_environment& environment,
    obj_instance& instance) {
  // read the file str by str
  char buffer[4096];
  auto found = false;
  while (read_obj_line(fs, buffer, sizeof(buffer))) {
    // str
    auto str = string_view{buffer};
    remove_obj_comment(str);
    skip_obj_whitespace(str);
    if (str.empty()) continue;

    // get command
    auto cmd = ""s;
    if (!parse_obj_value(str, cmd)) return {filename + ": parse error"};
    if (cmd == "") continue;

    // read values
    if (cmd == "c") {
      command = objx_command::camera;
      if (!parse_obj_value(str, camera.name))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.ortho))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.width))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.height))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.lens))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.focus))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.aperture))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, camera.frame))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "e") {
      command = objx_command::environment;
      if (!parse_obj_value(str, environment.name))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, environment.emission))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, environment.emission_map))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, environment.frame))
        return {filename + ": parse error"};
      return {};
    } else if (cmd == "i") {
      command = objx_command::instance;
      if (!parse_obj_value(str, instance.object))
        return {filename + ": parse error"};
      if (!parse_obj_value(str, instance.frame))
        return {filename + ": parse error"};
      return {};
    }
  }

  if (found) return {};

  // check error
  if (has_obj_error(fs)) return {filename + ": read error"};

  return {"eof"};
}

static vector<string> split_obj_string(const string& str, const string& delim) {
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
objio_status write_obj_comment(
    const string& filename, obj_file& fs, const string& comment) {
  auto lines = split_obj_string(comment, "\n");
  for (auto& str : lines) {
    if (!format_obj_values(fs, "# {}\n", str))
      return {filename + ": write error"};
  }
  if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
  return {};
}

objio_status write_obj_command(const string& filename, obj_file& fs,
    obj_command command, const string& name, const vec3f& value,
    const vector<obj_vertex>& vertices) {
  switch (command) {
    case obj_command::vertex:
      if (!format_obj_values(fs, "v {}\n", value))
        return {filename + ": write error"};
      break;
    case obj_command::normal:
      if (!format_obj_values(fs, "vn {}\n", value))
        return {filename + ": write error"};
      break;
    case obj_command::texcoord:
      if (!format_obj_values(fs, "vt {}\n", value))
        return {filename + ": write error"};
      break;
    case obj_command::face:
    case obj_command::str:
    case obj_command::point:
      if (command == obj_command::face)
        if (!format_obj_values(fs, "f ")) return {filename + ": write error"};
      if (command == obj_command::str)
        if (!format_obj_values(fs, "l ")) return {filename + ": write error"};
      if (command == obj_command::point)
        if (!format_obj_values(fs, "p ")) return {filename + ": write error"};
      for (auto& vert : vertices)
        if (!format_obj_values(fs, " {}", vert))
          return {filename + ": write error"};
      if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
      break;
    case obj_command::object:
      if (!format_obj_values(fs, "o {}\n", name))
        return {filename + ": write error"};
      break;
    case obj_command::group:
      if (!format_obj_values(fs, "g {}\n", name))
        return {filename + ": write error"};
      break;
    case obj_command::usemtl:
      if (!format_obj_values(fs, "usemtl {}\n", name))
        return {filename + ": write error"};
      break;
    case obj_command::smoothing:
      if (!format_obj_values(fs, "s {}\n", name))
        return {filename + ": write error"};
      break;
    case obj_command::mtllib:
      if (!format_obj_values(fs, "mtllib {}\n", name))
        return {filename + ": write error"};
      break;
    case obj_command::objxlib: break;
    case obj_command::error: break;
  }

  return {};
}

static objio_status write_mtl_command(const string& filename, obj_file& fs,
    mtl_command command, const obj_material& material) {
  // write material
  switch (command) {
    case mtl_command::material:
      if (!format_obj_values(fs, "newmtl {}\n", material.name))
        return {filename + ": write error"};
      if (!format_obj_values(fs, "illum {}\n", material.illum))
        return {filename + ": write error"};
      if (material.emission != zero3f)
        if (!format_obj_values(fs, "Ke {}\n", material.emission))
          return {filename + ": write error"};
      if (material.ambient != zero3f)
        if (!format_obj_values(fs, "Ka {}\n", material.ambient))
          return {filename + ": write error"};
      if (!format_obj_values(fs, "Kd {}\n", material.diffuse))
        return {filename + ": write error"};
      if (!format_obj_values(fs, "Ks {}\n", material.specular))
        return {filename + ": write error"};
      if (material.reflection != zero3f)
        if (!format_obj_values(fs, "Kr {}\n", material.reflection))
          return {filename + ": write error"};
      if (material.transmission != zero3f)
        if (!format_obj_values(fs, "Kt {}\n", material.transmission))
          return {filename + ": write error"};
      format_obj_values(fs, "Ns {}\n", (int)material.exponent);
      if (material.opacity != 1)
        if (!format_obj_values(fs, "d {}\n", material.opacity))
          return {filename + ": write error"};
      if (!material.emission_map.path.empty())
        if (!format_obj_values(fs, "map_Ke {}\n", material.emission_map))
          return {filename + ": write error"};
      if (!material.diffuse_map.path.empty())
        if (!format_obj_values(fs, "map_Kd {}\n", material.diffuse_map))
          return {filename + ": write error"};
      if (!material.specular_map.path.empty())
        if (!format_obj_values(fs, "map_Ks {}\n", material.specular_map))
          return {filename + ": write error"};
      if (!material.transmission_map.path.empty())
        if (!format_obj_values(fs, "map_Kt {}\n", material.transmission_map))
          return {filename + ": write error"};
      if (!material.reflection_map.path.empty())
        if (!format_obj_values(fs, "map_Kr {}\n", material.reflection_map))
          return {filename + ": write error"};
      if (!material.exponent_map.path.empty())
        if (!format_obj_values(fs, "map_Ns {}\n", material.exponent_map))
          return {filename + ": write error"};
      if (!material.opacity_map.path.empty())
        if (!format_obj_values(fs, "map_d {}\n", material.opacity_map))
          return {filename + ": write error"};
      if (!material.bump_map.path.empty())
        if (!format_obj_values(fs, "map_bump {}\n", material.bump_map))
          return {filename + ": write error"};
      if (!material.displacement_map.path.empty())
        if (!format_obj_values(fs, "map_disp {}\n", material.displacement_map))
          return {filename + ": write error"};
      if (!material.normal_map.path.empty())
        if (!format_obj_values(fs, "map_norm {}\n", material.normal_map))
          return {filename + ": write error"};
      if (material.pbr_roughness)
        if (!format_obj_values(fs, "Pr {}\n", material.pbr_roughness))
          return {filename + ": write error"};
      if (material.pbr_metallic)
        if (!format_obj_values(fs, "Pm {}\n", material.pbr_metallic))
          return {filename + ": write error"};
      if (material.pbr_sheen)
        if (!format_obj_values(fs, "Ps {}\n", material.pbr_sheen))
          return {filename + ": write error"};
      if (material.pbr_clearcoat)
        if (!format_obj_values(fs, "Pc {}\n", material.pbr_clearcoat))
          return {filename + ": write error"};
      if (material.pbr_coatroughness)
        if (!format_obj_values(fs, "Pcr {}\n", material.pbr_coatroughness))
          return {filename + ": write error"};
      if (!material.pbr_roughness_map.path.empty())
        if (!format_obj_values(fs, "map_Pr {}\n", material.pbr_roughness_map))
          return {filename + ": write error"};
      if (!material.pbr_metallic_map.path.empty())
        if (!format_obj_values(fs, "map_Pm {}\n", material.pbr_metallic_map))
          return {filename + ": write error"};
      if (!material.pbr_sheen_map.path.empty())
        if (!format_obj_values(fs, "map_Ps {}\n", material.pbr_sheen_map))
          return {filename + ": write error"};
      if (!material.pbr_clearcoat_map.path.empty())
        if (!format_obj_values(fs, "map_Pc {}\n", material.pbr_clearcoat_map))
          return {filename + ": write error"};
      if (!material.pbr_coatroughness_map.path.empty())
        if (!format_obj_values(
                fs, "map_Pcr {}\n", material.pbr_coatroughness_map))
          return {filename + ": write error"};
      if (material.vol_transmission != zero3f)
        if (!format_obj_values(fs, "Vt {}\n", material.vol_transmission))
          return {filename + ": write error"};
      if (material.vol_meanfreepath != zero3f)
        if (!format_obj_values(fs, "Vp {}\n", material.vol_meanfreepath))
          return {filename + ": write error"};
      if (material.vol_emission != zero3f)
        if (!format_obj_values(fs, "Ve {}\n", material.vol_emission))
          return {filename + ": write error"};
      if (material.vol_scattering != zero3f)
        if (!format_obj_values(fs, "Vs {}\n", material.vol_scattering))
          return {filename + ": write error"};
      if (material.vol_anisotropy)
        if (!format_obj_values(fs, "Vg {}\n", material.vol_anisotropy))
          return {filename + ": write error"};
      if (material.vol_scale)
        if (!format_obj_values(fs, "Vr {}\n", material.vol_scale))
          return {filename + ": write error"};
      if (!material.vol_scattering_map.path.empty())
        if (!format_obj_values(fs, "map_Vs {}\n", material.vol_scattering_map))
          return {filename + ": write error"};
      if (!format_obj_values(fs, "\n")) return {filename + ": write error"};
      break;
    case mtl_command::error: break;
  }

  return {};
}

objio_status write_objx_command(const string& filename, obj_file& fs,
    objx_command command, const obj_camera& camera,
    const obj_environment& environment, const obj_instance& instance) {
  switch (command) {
    case objx_command::camera: {
      if (!format_obj_values(fs, "c {} {} {} {} {} {} {} {}\n", camera.name,
              camera.ortho, camera.width, camera.height, camera.lens,
              camera.focus, camera.aperture, camera.frame))
        return {filename + ": write error"};
    } break;
    case objx_command::environment: {
      if (!format_obj_values(fs, "e {} {} {} {}\n", environment.name,
              environment.emission,
              environment.emission_map.path.empty()
                  ? "\"\""s
                  : environment.emission_map.path,
              environment.frame))
        return {filename + ": write error"};
    } break;
    case objx_command::instance: {
      if (!format_obj_values(fs, "i {} {}\n", instance.object, instance.frame))
        return {filename + ": write error"};
    } break;
    case objx_command::error: break;
  }

  return {};
}

}  // namespace yocto
