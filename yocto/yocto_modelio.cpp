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

// Check for errors
static bool has_ply_error(ply_file& fs) { return ferror(fs.fs); }

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
