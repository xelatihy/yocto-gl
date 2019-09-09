//
// Implementation for Yocto/ModelIO.
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

#include <algorithm>
#include <cinttypes>
#include <climits>
#include <cstdarg>
#include <limits>
#include <string_view>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

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

static inline void checked_fprintf(file_wrapper& fs, const char* fmt, ...) {
  va_list args1;
  va_start(args1, fmt);
  if (vfprintf(fs.fs, fmt, args1) < 0)
    throw std::runtime_error("cannot write to file");
  va_end(args1);
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

// Parse values from a string
static inline void parse_ply_value(string_view& str, string_view& value) {
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
static inline void parse_ply_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_ply_value(str, valuev);
  value = string{valuev};
}
static inline void parse_ply_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}

// get ply value either ascii or binary
template <typename T>
static inline T read_ply_value(file_wrapper& fs, bool big_endian) {
  auto value = (T)0;
  if (fread(&value, sizeof(T), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read value");
  if (big_endian) value = swap_endian(value);
  return value;
}
template <typename VT>
static inline void read_ply_prop(
    file_wrapper& fs, bool big_endian, ply_type type, VT& value) {
  switch (type) {
    case ply_type::i8:
      value = (VT)read_ply_value<int8_t>(fs, big_endian);
      break;
    case ply_type::i16:
      value = (VT)read_ply_value<int16_t>(fs, big_endian);
      break;
    case ply_type::i32:
      value = (VT)read_ply_value<int32_t>(fs, big_endian);
      break;
    case ply_type::i64:
      value = (VT)read_ply_value<int64_t>(fs, big_endian);
      break;
    case ply_type::u8:
      value = (VT)read_ply_value<uint8_t>(fs, big_endian);
      break;
    case ply_type::u16:
      value = (VT)read_ply_value<uint16_t>(fs, big_endian);
      break;
    case ply_type::u32:
      value = (VT)read_ply_value<uint32_t>(fs, big_endian);
      break;
    case ply_type::u64:
      value = (VT)read_ply_value<uint64_t>(fs, big_endian);
      break;
    case ply_type::f32:
      value = (VT)read_ply_value<float>(fs, big_endian);
      break;
    case ply_type::f64:
      value = (VT)read_ply_value<double>(fs, big_endian);
      break;
  }
}

template <typename VT>
static inline void parse_ply_prop(string_view& str, ply_type type, VT& value) {
  char* end = nullptr;
  switch (type) {
    case ply_type::i8: value = (VT)strtol(str.data(), &end, 10); break;
    case ply_type::i16: value = (VT)strtol(str.data(), &end, 10); break;
    case ply_type::i32: value = (VT)strtol(str.data(), &end, 10); break;
    case ply_type::i64: value = (VT)strtoll(str.data(), &end, 10); break;
    case ply_type::u8: value = (VT)strtoul(str.data(), &end, 10); break;
    case ply_type::u16: value = (VT)strtoul(str.data(), &end, 10); break;
    case ply_type::u32: value = (VT)strtoul(str.data(), &end, 10); break;
    case ply_type::u64: value = (VT)strtoull(str.data(), &end, 10); break;
    case ply_type::f32: value = (VT)strtof(str.data(), &end); break;
    case ply_type::f64: value = (VT)strtod(str.data(), &end); break;
  }
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}

// Load ply data
void read_ply_header(file_wrapper& fs, ply_format& format,
    vector<ply_element>& elements, vector<string>& comments) {
  // ply type names
  static auto type_map = unordered_map<string, ply_type>{
      {"char", ply_type::i8},
      {"short", ply_type::i16},
      {"int", ply_type::i32},
      {"long", ply_type::i64},
      {"uchar", ply_type::u8},
      {"ushort", ply_type::u16},
      {"uint", ply_type::u32},
      {"ulong", ply_type::u64},
      {"float", ply_type::f32},
      {"double", ply_type::f64},
      {"int8", ply_type::i8},
      {"int16", ply_type::i16},
      {"int32", ply_type::i32},
      {"int64", ply_type::i64},
      {"uint8", ply_type::u8},
      {"uint16", ply_type::u16},
      {"uint32", ply_type::u32},
      {"uint64", ply_type::u64},
      {"float32", ply_type::f32},
      {"float64", ply_type::f64},
  };

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
    parse_ply_value(line, cmd);
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
      parse_ply_value(line, fmt);
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
      parse_ply_value(line, elem.name);
      parse_ply_value(line, elem.count);
    } else if (cmd == "property") {
      if (elements.empty()) throw std::runtime_error{"bad ply header"};
      auto& prop  = elements.back().properties.emplace_back();
      auto  tname = ""s;
      parse_ply_value(line, tname);
      if (tname == "list") {
        prop.is_list = true;
        auto ename   = ""s;
        parse_ply_value(line, tname);
        try {
          prop.value_type = type_map.at(tname);
        } catch (...) {
          throw std::runtime_error{"unknown ply type " + tname};
        }
        parse_ply_value(line, tname);
        try {
          prop.list_type = type_map.at(tname);
        } catch (...) {
          throw std::runtime_error{"unknown ply type " + tname};
        }
      } else {
        prop.is_list = false;
        try {
          prop.value_type = type_map.at(tname);
        } catch (...) {
          throw std::runtime_error{"unknown ply type " + tname};
        }
      }
      parse_ply_value(line, prop.name);
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
void read_ply_value_impl(file_wrapper& fs, ply_format format,
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
      parse_ply_prop(line, prop.value_type, value);
      if (prop.is_list) {
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          parse_ply_prop(line, prop.list_type, list[i]);
      }
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop  = element.properties[pidx];
      auto& value = values[pidx];
      auto& list  = lists[pidx];
      read_ply_prop(
          fs, format == ply_format::binary_big_endian, prop.value_type, value);
      if (prop.is_list) {
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          read_ply_prop(fs, format == ply_format::binary_big_endian,
              prop.list_type, list[i]);
      }
    }
  }
}

// Write text to file
static inline void write_ply_text(file_wrapper& fs, const char* value) {
  checked_fprintf(fs, "%s", value);
}
static inline void write_ply_text(file_wrapper& fs, const string& value) {
  checked_fprintf(fs, "%s", value.c_str());
}

template <typename VT>
static inline void write_ply_prop(file_wrapper& fs, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: checked_fprintf(fs, "%d", (int)value); break;
    case ply_type::i16: checked_fprintf(fs, "%d", (int)value); break;
    case ply_type::i32: checked_fprintf(fs, "%d", (int)value); break;
    case ply_type::i64: checked_fprintf(fs, "%lld", (long long)value); break;
    case ply_type::u8: checked_fprintf(fs, "%u", (unsigned)value); break;
    case ply_type::u16: checked_fprintf(fs, "%u", (unsigned)value); break;
    case ply_type::u32: checked_fprintf(fs, "%u", (unsigned)value); break;
    case ply_type::u64:
      checked_fprintf(fs, "%llu", (unsigned long long)value);
      break;
    case ply_type::f32: checked_fprintf(fs, "%g", (float)value); break;
    case ply_type::f64: checked_fprintf(fs, "%g", (double)value); break;
  }
}

template <typename T, typename VT>
static inline void write_ply_binprop(
    file_wrapper& fs, bool big_endian, VT value) {
  auto typed_value = (T)value;
  if (big_endian) typed_value = swap_endian(typed_value);
  if (fwrite(&typed_value, sizeof(T), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to file");
}

template <typename VT>
static inline void write_ply_binprop(
    file_wrapper& fs, bool big_endian, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: write_ply_binprop<int8_t>(fs, big_endian, value); break;
    case ply_type::i16:
      write_ply_binprop<int16_t>(fs, big_endian, value);
      break;
    case ply_type::i32:
      write_ply_binprop<int32_t>(fs, big_endian, value);
      break;
    case ply_type::i64:
      write_ply_binprop<int64_t>(fs, big_endian, value);
      break;
    case ply_type::u8: write_ply_binprop<uint8_t>(fs, big_endian, value); break;
    case ply_type::u16:
      write_ply_binprop<uint16_t>(fs, big_endian, value);
      break;
    case ply_type::u32:
      write_ply_binprop<uint32_t>(fs, big_endian, value);
      break;
    case ply_type::u64:
      write_ply_binprop<uint64_t>(fs, big_endian, value);
      break;
    case ply_type::f32: write_ply_binprop<float>(fs, big_endian, value); break;
    case ply_type::f64: write_ply_binprop<double>(fs, big_endian, value); break;
  }
}

// Write Ply functions
void write_ply_header(file_wrapper& fs, ply_format format,
    const vector<ply_element>& elements, const vector<string>& comments) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{
      {ply_type::i8, "char"},
      {ply_type::i16, "short"},
      {ply_type::i32, "int"},
      {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"},
      {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"},
      {ply_type::u64, "ulong"},
      {ply_type::f32, "float"},
      {ply_type::f64, "double"},
  };

  write_ply_text(fs, "ply\n");
  switch (format) {
    case ply_format::ascii: write_ply_text(fs, "format ascii 1.0\n"); break;
    case ply_format::binary_little_endian:
      write_ply_text(fs, "format binary_little_endian 1.0\n");
      break;
    case ply_format::binary_big_endian:
      write_ply_text(fs, "format binary_big_endian 1.0\n");
      break;
  }
  for (auto& comment : comments)
    write_ply_text(fs, "comment " + comment + "\n");
  for (auto& elem : elements) {
    write_ply_text(
        fs, "element " + elem.name + " " + std::to_string(elem.count) + "\n");
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        write_ply_text(fs, "property list " + type_map[prop.value_type] + " " +
                               type_map[prop.list_type] + " " + prop.name +
                               "\n");
      } else {
        write_ply_text(fs,
            "property " + type_map[prop.value_type] + " " + prop.name + "\n");
      }
    }
  }
  write_ply_text(fs, "end_header\n");
}

template <typename VT, typename LT>
void write_ply_value_impl(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<VT>& values, vector<vector<LT>>& lists) {
  if (format == ply_format::ascii) {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (pidx) write_ply_text(fs, " ");
      write_ply_prop(fs, prop.value_type, values[pidx]);
      if (prop.is_list) {
        for (auto i = 0; i < (int)lists[pidx].size(); i++) {
          if (i) write_ply_text(fs, " ");
          write_ply_prop(fs, prop.list_type, lists[pidx][i]);
        }
      }
      write_ply_text(fs, "\n");
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      write_ply_binprop(fs, format == ply_format::binary_big_endian,
          prop.value_type, values[pidx]);
      if (prop.is_list) {
        for (auto i = 0; i < (int)lists[pidx].size(); i++)
          write_ply_binprop(fs, format == ply_format::binary_big_endian,
              prop.list_type, lists[pidx][i]);
      }
    }
  }
}

void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  write_ply_value_impl(fs, format, element, values, lists);
}
void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  write_ply_value_impl(fs, format, element, values, lists);
}

void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  read_ply_value_impl(fs, format, element, values, lists);
}
void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  read_ply_value_impl(fs, format, element, values, lists);
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

// Parse values from a string
static inline void parse_obj_value(string_view& str, string_view& value) {
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
static inline void parse_obj_value(string_view& str, string& value) {
  auto valuev = ""sv;
  parse_obj_value(str, valuev);
  value = string{valuev};
}
static inline void parse_obj_value(string_view& str, int& value) {
  char* end = nullptr;
  value     = (int)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_obj_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
static inline void parse_obj_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_obj_value(str, values[i]);
}

static inline void parse_obj_value(string_view& str, vec2f& value) {
  parse_obj_value(str, &value.x, 2);
}
static inline void parse_obj_value(string_view& str, vec3f& value) {
  parse_obj_value(str, &value.x, 3);
}
static inline void parse_obj_value(string_view& str, frame3f& value) {
  parse_obj_value(str, &value.x.x, 12);
}

template <typename T>
static inline void parse_obj_value_or_empty(string_view& str, T& value) {
  skip_whitespace(str);
  if (str.empty()) {
    value = T{};
  } else {
    parse_obj_value(str, value);
  }
}

static inline void parse_obj_value(string_view& str, obj_vertex& value) {
  value = obj_vertex{0, 0, 0};
  parse_obj_value(str, value.position);
  if (!str.empty() && str.front() == '/') {
    str.remove_prefix(1);
    if (!str.empty() && str.front() == '/') {
      str.remove_prefix(1);
      parse_obj_value(str, value.normal);
    } else {
      parse_obj_value(str, value.texcoord);
      if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        parse_obj_value(str, value.normal);
      }
    }
  }
}

// Input for OBJ textures
static inline void parse_obj_value(string_view& str, obj_texture_info& info) {
  // initialize
  info = obj_texture_info();

  // get tokens
  auto tokens = vector<string>();
  skip_whitespace(str);
  while (!str.empty()) {
    auto token = ""s;
    parse_obj_value(str, token);
    tokens.push_back(token);
    skip_whitespace(str);
  }
  if (tokens.empty()) throw std::runtime_error("cannot parse value");

  // texture name
  info.path = fs::path(tokens.back()).generic_string();

  // texture params
  auto last = string();
  for (auto i = 0; i < tokens.size() - 1; i++) {
    if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
    if (tokens[i] == "-clamp") info.clamp = true;
  }
}

void parse_obj_value(string_view& str, obj_value& value, obj_value_type type,
    int array_size = 3) {
  switch (type) {
    case obj_value_type::number: {
      auto value_ = 0.0f;
      parse_obj_value(str, value_);
      value = make_obj_value(value_);
    } break;
    case obj_value_type::string: {
      auto value_ = ""s;
      parse_obj_value(str, value_);
      value = make_obj_value(value_);
    } break;
    case obj_value_type::array: {
      if (array_size == 2) {
        auto value_ = zero2f;
        parse_obj_value(str, value_);
        value = make_obj_value(value_);
      } else if (array_size == 3) {
        auto value_ = zero3f;
        parse_obj_value(str, value_);
        value = make_obj_value(value_);
      } else if (array_size == 12) {
        auto value_ = identity3x4f;
        parse_obj_value(str, value_);
        value = make_obj_value(value_);
      } else {
        throw std::runtime_error("should not have gotten here");
      }
    } break;
    case obj_value_type::boolean: {
      auto value_ = 0;
      parse_obj_value(str, value_);
      value = make_obj_value((bool)value_);
    } break;
  }
}

static inline void parse_obj_value_or_empty(
    string_view& str, obj_value& value) {
  skip_whitespace(str);
  if (str.empty()) {
    value = make_obj_value(""s);
  } else {
    parse_obj_value(str, value, obj_value_type::string);
  }
}

// Read obj
bool read_obj_command(file_wrapper& fs, obj_command& command, obj_value& value,
    vector<obj_vertex>& vertices, obj_vertex& vert_size) {
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
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "v") {
      command = obj_command::vertex;
      parse_obj_value(line, value, obj_value_type::array);
      vert_size.position += 1;
      return true;
    } else if (cmd == "vn") {
      command = obj_command::normal;
      parse_obj_value(line, value, obj_value_type::array);
      vert_size.normal += 1;
      return true;
    } else if (cmd == "vt") {
      command = obj_command::texcoord;
      parse_obj_value(line, value, obj_value_type::array, 2);
      vert_size.texcoord += 1;
      return true;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      vertices.clear();
      skip_whitespace(line);
      while (!line.empty()) {
        auto vert = obj_vertex{};
        parse_obj_value(line, vert);
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
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "usemtl") {
      command = obj_command::usemtl;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "g") {
      command = obj_command::group;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "s") {
      command = obj_command::smoothing;
      parse_obj_value_or_empty(line, value);
      return true;
    } else if (cmd == "mtllib") {
      command = obj_command::mtllib;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else {
      // unused
    }
  }
  return false;
}

// Read mtl
bool read_mtl_command(file_wrapper& fs, mtl_command& command, obj_value& value,
    obj_texture_info& texture, bool fliptr) {
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
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // possible token values
    if (cmd == "newmtl") {
      command = mtl_command::material;
      parse_obj_value(line, value, obj_value_type::string);
    } else if (cmd == "illum") {
      command = mtl_command::illum;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Ke") {
      command = mtl_command::emission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Kd") {
      command = mtl_command::diffuse;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Ks") {
      command = mtl_command::specular;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Kt") {
      command = mtl_command::transmission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Tf") {
      command    = mtl_command::transmission;
      auto color = vec3f{-1};
      value      = make_obj_value(color);
      parse_obj_value(line, value, obj_value_type::array);
      get_obj_value(value, color);
      if (color.y < 0) color = vec3f{color.x};
      if (fliptr) color = 1 - color;
      value = make_obj_value(color);
    } else if (cmd == "Tr") {
      command = mtl_command::opacity;
      parse_obj_value(line, value, obj_value_type::number);
      if (fliptr) value.number = 1 - value.number;
    } else if (cmd == "Ns") {
      command = mtl_command::exponent;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "d") {
      command = mtl_command::opacity;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Ke") {
      command = mtl_command::emission_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Kd") {
      command = mtl_command::diffuse_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Ks") {
      command = mtl_command::specular_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Tr") {
      command = mtl_command::transmission_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      command = mtl_command::opacity_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_bump" || cmd == "bump") {
      command = mtl_command::bump_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_disp" || cmd == "disp") {
      command = mtl_command::displacement_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_norm" || cmd == "norm") {
      command = mtl_command::normal_map;
      parse_obj_value(line, texture);
    } else if (cmd == "Pm") {
      command = mtl_command::pbr_metallic;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pr") {
      command = mtl_command::pbr_roughness;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Ps") {
      command = mtl_command::pbr_sheen;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pc") {
      command = mtl_command::pbr_clearcoat;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Pcr") {
      command = mtl_command::pbr_coatroughness;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Pm") {
      command = mtl_command::pbr_metallic_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pr") {
      command = mtl_command::pbr_roughness_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Ps") {
      command = mtl_command::pbr_sheen_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pc") {
      command = mtl_command::pbr_clearcoat_map;
      parse_obj_value(line, texture);
    } else if (cmd == "map_Pcr") {
      command = mtl_command::pbr_coatroughness_map;
      parse_obj_value(line, texture);
    } else if (cmd == "Vt") {
      command = mtl_command::vol_transmission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vp") {
      command = mtl_command::vol_meanfreepath;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Ve") {
      command = mtl_command::vol_emission;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vs") {
      command = mtl_command::vol_scattering;
      parse_obj_value(line, value, obj_value_type::array);
    } else if (cmd == "Vg") {
      command = mtl_command::vol_anisotropy;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "Vr") {
      command = mtl_command::vol_scale;
      parse_obj_value(line, value, obj_value_type::number);
    } else if (cmd == "map_Vs") {
      command = mtl_command::vol_scattering_map;
      parse_obj_value(line, texture);
    } else {
      continue;
    }

    return true;
  }

  return false;
}

// Read objx
bool read_objx_command(file_wrapper& fs, objx_command& command,
    obj_value& value, obj_texture_info& texture) {
  // read the file line by line
  char buffer[4096];
  auto pos = ftell(fs.fs);
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_obj_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // get command
    auto cmd = ""s;
    parse_obj_value(line, cmd);
    if (cmd == "") continue;

    // read values
    if (cmd == "newcam") {
      command = objx_command::camera;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newenv") {
      command = objx_command::environment;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newist") {
      command = objx_command::instance;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "newproc") {
      command = objx_command::procedural;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "frame") {
      command = objx_command::frame;
      parse_obj_value(line, value, obj_value_type::array, 12);
      return true;
    } else if (cmd == "obj") {
      command = objx_command::object;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "mat") {
      command = objx_command::material;
      parse_obj_value(line, value, obj_value_type::string);
      return true;
    } else if (cmd == "ortho") {
      command = objx_command::ortho;
      parse_obj_value(line, value, obj_value_type::boolean);
      return true;
    } else if (cmd == "width") {
      command = objx_command::width;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "height") {
      command = objx_command::height;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "lens") {
      command = objx_command::lens;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "aperture") {
      command = objx_command::aperture;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "focus") {
      command = objx_command::focus;
      parse_obj_value(line, value, obj_value_type::number);
      return true;
    } else if (cmd == "Ke") {
      command = objx_command::emission;
      parse_obj_value(line, value, obj_value_type::array);
      return true;
    } else if (cmd == "map_Ke") {
      command = objx_command::emission_map;
      parse_obj_value(line, texture);
      return true;
    }
    // backward compatibility
    else if (cmd == "c") {
      auto oname = value.string_;
      auto name = obj_value{}, ortho = obj_value{}, width = obj_value{},
           height = obj_value{}, lens = obj_value{}, aperture = obj_value{},
           focus = obj_value{}, frame = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, ortho, obj_value_type::boolean);
      parse_obj_value(line, width, obj_value_type::number);
      parse_obj_value(line, height, obj_value_type::number);
      parse_obj_value(line, lens, obj_value_type::number);
      parse_obj_value(line, focus, obj_value_type::number);
      parse_obj_value(line, aperture, obj_value_type::number);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::camera && oname != "") {
        command = objx_command::ortho;
        value   = ortho;
      } else if (command == objx_command::ortho) {
        command = objx_command::width;
        value   = width;
      } else if (command == objx_command::width) {
        command = objx_command::height;
        value   = height;
      } else if (command == objx_command::height) {
        command = objx_command::lens;
        value   = lens;
      } else if (command == objx_command::lens) {
        command = objx_command::focus;
        value   = focus;
      } else if (command == objx_command::focus) {
        command = objx_command::aperture;
        value   = aperture;
      } else if (command == objx_command::aperture) {
        command = objx_command::frame;
        value   = frame;
      } else {
        command = objx_command::camera;
        value   = name;
      }
      if (command != objx_command::frame) fseek(fs.fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "e") {
      auto name = obj_value{}, frame = obj_value{}, emission = obj_value{},
           emission_map = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, emission, obj_value_type::array);
      parse_obj_value(line, emission_map, obj_value_type::string);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (emission_map.string_ == "\"\"") emission_map.string_ = "";
      if (command == objx_command::environment) {
        command = objx_command::emission;
        value   = emission;
      } else if (command == objx_command::emission) {
        command = objx_command::emission_map;
        get_obj_value(emission_map, texture.path);
      } else if (command == objx_command::emission_map) {
        command = objx_command::frame;
        value   = frame;
      } else {
        command = objx_command::environment;
        value   = name;
      }
      if (command != objx_command::frame) fseek(fs.fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "i") {
      auto name = obj_value{}, frame = obj_value{}, object = obj_value{},
           material = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, object, obj_value_type::string);
      parse_obj_value(line, material, obj_value_type::string);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::instance) {
        command = objx_command::object;
        value   = object;
      } else if (command == objx_command::object) {
        command = objx_command::material;
        value   = material;
      } else if (command == objx_command::material) {
        command = objx_command::frame;
        value   = frame;
      } else {
        command = objx_command::instance;
        value   = name;
      }
      if (command != objx_command::frame) fseek(fs.fs, pos, SEEK_SET);
      return true;
    } else if (cmd == "po") {
      auto name = obj_value{}, frame = obj_value{}, type = obj_value{},
           material = obj_value{}, size = obj_value{}, level = obj_value{};
      parse_obj_value(line, name, obj_value_type::string);
      parse_obj_value(line, type, obj_value_type::string);
      parse_obj_value(line, material, obj_value_type::string);
      parse_obj_value(line, size, obj_value_type::number);
      parse_obj_value(line, level, obj_value_type::number);
      parse_obj_value(line, frame, obj_value_type::array, 12);
      if (command == objx_command::procedural) {
        command = objx_command::object;
        value   = type;
      } else if (command == objx_command::object) {
        command = objx_command::material;
        value   = material;
      } else if (command == objx_command::material) {
        command = objx_command::frame;
        value   = frame;
      } else {
        command = objx_command::procedural;
        value   = name;
      }
      if (command != objx_command::frame) fseek(fs.fs, pos, SEEK_SET);
      return true;
    } else {
      // unused
    }
  }

  return false;
}

// Write obj elements
void write_obj_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    checked_fprintf(fs, "# %s\n", line.c_str());
  }
  checked_fprintf(fs, "\n");
}

void write_obj_command(file_wrapper& fs, obj_command command,
    const obj_value& value_, const vector<obj_vertex>& vertices) {
  auto& name  = value_.string_;
  auto& value = value_.array_;
  switch (command) {
    case obj_command::vertex:
      checked_fprintf(fs, "v %g %g %g\n", value[0], value[1], value[2]);
      break;
    case obj_command::normal:
      checked_fprintf(fs, "vn %g  %g %g\n", value[0], value[1], value[2]);
      break;
    case obj_command::texcoord:
      checked_fprintf(fs, "vt %g %g\n", value[0], value[1]);
      break;
    case obj_command::face:
    case obj_command::line:
    case obj_command::point:
      if (command == obj_command::face) checked_fprintf(fs, "f ");
      if (command == obj_command::line) checked_fprintf(fs, "l ");
      if (command == obj_command::point) checked_fprintf(fs, "p ");
      for (auto& vert : vertices) {
        checked_fprintf(fs, " ");
        checked_fprintf(fs, "%d", vert.position);
        if (vert.texcoord) {
          checked_fprintf(fs, "/%d", vert.texcoord);
          if (vert.normal) {
            checked_fprintf(fs, "/%d", vert.normal);
          }
        } else if (vert.normal) {
          checked_fprintf(fs, "//%d", vert.normal);
        }
      }
      checked_fprintf(fs, "\n");
      break;
    case obj_command::object:
      checked_fprintf(fs, "o %s\n", name.c_str());
      break;
    case obj_command::group: checked_fprintf(fs, "g %s\n", name.c_str()); break;
    case obj_command::usemtl:
      checked_fprintf(fs, "usemtl %s\n", name.c_str());
      break;
    case obj_command::smoothing:
      checked_fprintf(fs, "s %s\n", name.c_str());
      break;
    case obj_command::mtllib:
      checked_fprintf(fs, "mtllib %s\n", name.c_str());
      break;
    case obj_command::objxlib: break;
  }
}

void write_mtl_command(file_wrapper& fs, mtl_command command,
    const obj_value& value_, const obj_texture_info& texture) {
  auto& name  = value_.string_;
  auto  value = value_.number;
  auto& color = value_.array_;
  switch (command) {
    case mtl_command::material:
      checked_fprintf(fs, "\nnewmtl %s\n", name.c_str());
      break;
    case mtl_command::illum:
      checked_fprintf(fs, "  illum %d\n", (int)value);
      break;
    case mtl_command::emission:
      checked_fprintf(fs, "  Ke %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::ambient:
      checked_fprintf(fs, "  Ka %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::diffuse:
      checked_fprintf(fs, "  Kd %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::specular:
      checked_fprintf(fs, "  Ks %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::reflection:
      checked_fprintf(fs, "  Kr %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::transmission:
      checked_fprintf(fs, "  Kt %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::exponent:
      checked_fprintf(fs, "  Ns %d\n", (int)value);
      break;
    case mtl_command::opacity: checked_fprintf(fs, "  d %g\n", value); break;
    case mtl_command::ior: checked_fprintf(fs, "  Ni %g\n", value); break;
    case mtl_command::emission_map:
      checked_fprintf(fs, "  map_Ke %s\n", texture.path.c_str());
      break;
    case mtl_command::ambient_map:
      checked_fprintf(fs, "  map_Ka %s\n", texture.path.c_str());
      break;
    case mtl_command::diffuse_map:
      checked_fprintf(fs, "  map_Kd %s\n", texture.path.c_str());
      break;
    case mtl_command::specular_map:
      checked_fprintf(fs, "  map_Ks %s\n", texture.path.c_str());
      break;
    case mtl_command::reflection_map:
      checked_fprintf(fs, "  map_Kr %s\n", texture.path.c_str());
      break;
    case mtl_command::transmission_map:
      checked_fprintf(fs, "  map_Kt %s\n", texture.path.c_str());
      break;
    case mtl_command::opacity_map:
      checked_fprintf(fs, "  map_d %s\n", texture.path.c_str());
      break;
    case mtl_command::exponent_map:
      checked_fprintf(fs, "  map_Ni %s\n", texture.path.c_str());
      break;
    case mtl_command::bump_map:
      checked_fprintf(fs, "  map_bump %s\n", texture.path.c_str());
      break;
    case mtl_command::normal_map:
      checked_fprintf(fs, "  map_norm %s\n", texture.path.c_str());
      break;
    case mtl_command::displacement_map:
      checked_fprintf(fs, "  map_disp %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_roughness:
      checked_fprintf(fs, "  Pr %g\n", value);
      break;
    case mtl_command::pbr_metallic:
      checked_fprintf(fs, "  Pm %g\n", value);
      break;
    case mtl_command::pbr_sheen: checked_fprintf(fs, "  Ps %g\n", value); break;
    case mtl_command::pbr_clearcoat:
      checked_fprintf(fs, "  Pc %g\n", value);
      break;
    case mtl_command::pbr_coatroughness:
      checked_fprintf(fs, "  Pcr %g\n", value);
      break;
    case mtl_command::pbr_roughness_map:
      checked_fprintf(fs, "  Pr_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_metallic_map:
      checked_fprintf(fs, "  Pm_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_sheen_map:
      checked_fprintf(fs, "  Ps_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_clearcoat_map:
      checked_fprintf(fs, "  Pc_map %s\n", texture.path.c_str());
      break;
    case mtl_command::pbr_coatroughness_map:
      checked_fprintf(fs, "  Pcr_map %s\n", texture.path.c_str());
      break;
    case mtl_command::vol_transmission:
      checked_fprintf(fs, "  Vt %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_meanfreepath:
      checked_fprintf(fs, "  Vp %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_emission:
      checked_fprintf(fs, "  Ve %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_scattering:
      checked_fprintf(fs, "  Vs %g %g %g\n", color[0], color[1], color[2]);
      break;
    case mtl_command::vol_anisotropy:
      checked_fprintf(fs, "  Vg %g\n", value);
      break;
    case mtl_command::vol_scale: checked_fprintf(fs, "  Vr %g\n", value); break;
    case mtl_command::vol_scattering_map:
      checked_fprintf(fs, "  Vs_map %s\n", texture.path.c_str());
  }
}

void write_objx_command(file_wrapper& fs, objx_command command,
    const obj_value& value_, const obj_texture_info& texture) {
  auto& name  = value_.string_;
  auto  value = value_.number;
  auto& color = value_.array_;
  auto& frame = value_.array_;
  switch (command) {
    case objx_command::camera:
      checked_fprintf(fs, "\nnewcam %s\n", name.c_str());
      break;
    case objx_command::environment:
      checked_fprintf(fs, "\nnewenv %s\n", name.c_str());
      break;
    case objx_command::instance:
      checked_fprintf(fs, "\nnewist %s\n", name.c_str());
      break;
    case objx_command::procedural:
      checked_fprintf(fs, "\nnewproc %s\n", name.c_str());
      break;
    case objx_command::frame:
      checked_fprintf(fs, "  frame %g %g %g %g %g %g %g %g %g %g %g %g\n",
          frame[0], frame[1], frame[2], frame[3], frame[4], frame[5], frame[6],
          frame[7], frame[8], frame[9], frame[10], frame[11]);
      break;
    case objx_command::object:
      checked_fprintf(fs, "  obj %s\n", name.c_str());
      break;
    case objx_command::material:
      checked_fprintf(fs, "  mat %s\n", name.c_str());
      break;
    case objx_command::ortho: checked_fprintf(fs, "  ortho %g\n", value); break;
    case objx_command::width: checked_fprintf(fs, "  width %g\n", value); break;
    case objx_command::height:
      checked_fprintf(fs, "  height %g\n", value);
      break;
    case objx_command::lens: checked_fprintf(fs, "  lens %g\n", value); break;
    case objx_command::aperture:
      checked_fprintf(fs, "  aperture %g\n", value);
      break;
    case objx_command::focus: checked_fprintf(fs, "  focus %g\n", value); break;
    case objx_command::emission:
      checked_fprintf(fs, "  Ke %g %g %g\n", color[0], color[1], color[2]);
      break;
    case objx_command::emission_map:
      checked_fprintf(fs, "  map_Ke %s\n", texture.path.c_str());
      break;
  }
}

// typesafe access of obj value
void get_obj_value(const obj_value& yaml, string& value) {
  if (yaml.type != obj_value_type::string)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.string_;
}
void get_obj_value(const obj_value& yaml, bool& value) {
  if (yaml.type != obj_value_type::boolean)
    throw std::runtime_error("error parsing yaml value");
  value = yaml.boolean;
}
void get_obj_value(const obj_value& yaml, int& value) {
  if (yaml.type != obj_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (int)yaml.number;
}
void get_obj_value(const obj_value& yaml, float& value) {
  if (yaml.type != obj_value_type::number)
    throw std::runtime_error("error parsing yaml value");
  value = (float)yaml.number;
}
void get_obj_value(const obj_value& yaml, vec2f& value) {
  if (yaml.type != obj_value_type::array || yaml.number != 2)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1]};
}
void get_obj_value(const obj_value& yaml, vec3f& value) {
  if (yaml.type != obj_value_type::array || yaml.number != 3)
    throw std::runtime_error("error parsing yaml value");
  value = {(float)yaml.array_[0], (float)yaml.array_[1], (float)yaml.array_[2]};
}
void get_obj_value(const obj_value& yaml, mat3f& value) {
  if (yaml.type != obj_value_type::array || yaml.number != 9)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 9; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}
void get_obj_value(const obj_value& yaml, frame3f& value) {
  if (yaml.type != obj_value_type::array || yaml.number != 12)
    throw std::runtime_error("error parsing yaml value");
  for (auto i = 0; i < 12; i++) (&value.x.x)[i] = (float)yaml.array_[i];
}

// typesafe access of obj value
obj_value make_obj_value(const string& value) {
  return {obj_value_type::string, 0, false, value};
}
obj_value make_obj_value(bool value) {
  return {obj_value_type::boolean, 0, value};
}
obj_value make_obj_value(int value) {
  return {obj_value_type::number, (double)value};
}
obj_value make_obj_value(float value) {
  return {obj_value_type::number, (double)value};
}
obj_value make_obj_value(const vec2f& value) {
  return {
      obj_value_type::array, 2, false, "", {(double)value.x, (double)value.y}};
}
obj_value make_obj_value(const vec3f& value) {
  return {obj_value_type::array, 3, false, "",
      {(double)value.x, (double)value.y, (double)value.z}};
}
obj_value make_obj_value(const mat3f& value) {
  auto yaml = obj_value{obj_value_type::array, 9};
  for (auto i = 0; i < 9; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
}
obj_value make_obj_value(const frame3f& value) {
  auto yaml = obj_value{obj_value_type::array, 12};
  for (auto i = 0; i < 12; i++) yaml.array_[i] = (double)(&value.x.x)[i];
  return yaml;
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
    checked_fprintf(fs, "# %s\n", line.c_str());
  }
  checked_fprintf(fs, "\n");
}

// Save yaml property
void write_yaml_property(file_wrapper& fs, const string& object,
    const string& key, bool newobj, const yaml_value& value) {
  if (key.empty()) {
    checked_fprintf(fs, "\n%s:\n", object.c_str());
  } else {
    if (!object.empty()) {
      checked_fprintf(fs, (newobj ? "  - " : "    "));
    }
    checked_fprintf(fs, "%s: ", key.c_str());
    switch (value.type) {
      case yaml_value_type::number:
        checked_fprintf(fs, "%g", value.number);
        break;
      case yaml_value_type::boolean:
        checked_fprintf(fs, "%s", value.boolean ? "true" : "false");
        break;
      case yaml_value_type::string:
        checked_fprintf(fs, "%s", value.string_.c_str());
        break;
      case yaml_value_type::array:
        checked_fprintf(fs, "[ ");
        for (auto i = 0; i < value.number; i++) {
          if (i) checked_fprintf(fs, ", ");
          checked_fprintf(fs, "%g", value.array_[i]);
        }
        checked_fprintf(fs, " ]");
        break;
    }
    checked_fprintf(fs, "\n", key.c_str());
  }
}

void write_yaml_object(file_wrapper& fs, const string& object) {
  checked_fprintf(fs, "\n%s:\n", object.c_str());
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
bool read_pbrt_cmdline(file_wrapper& fs, string& cmd) {
  char buffer[4096];
  cmd.clear();
  auto found = false;
  auto pos   = ftell(fs.fs);
  while (read_line(fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_pbrt_comment(line);
    skip_whitespace(line);
    if (line.empty()) continue;

    // check if command
    auto is_cmd = line[0] >= 'A' && line[0] <= 'Z';
    if (is_cmd) {
      if (found) {
        fseek(fs.fs, pos, SEEK_SET);
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
        auto filenamep = fs::path(filename).filename();
        if (filenamep.extension() == ".spd") {
          filenamep = filenamep.replace_extension("");
          if (filenamep == "SHPS") {
            value.value3f = {1, 1, 1};
          } else if (filenamep.extension() == ".eta") {
            auto eta = get_pbrt_etak(filenamep.replace_extension("")).first;
            value.value3f = {eta.x, eta.y, eta.z};
          } else if (filenamep.extension() == ".k") {
            auto k = get_pbrt_etak(filenamep.replace_extension("")).second;
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

// Read pbrt commands
bool read_pbrt_command(file_wrapper& fs, pbrt_command_& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values, string& line) {
  // parse command by command
  while (read_pbrt_cmdline(fs, line)) {
    auto str = string_view{line};
    // get command
    auto cmd = ""s;
    parse_pbrt_command(str, cmd);
    if (cmd == "WorldBegin") {
      command = pbrt_command_::world_begin;
      return true;
    } else if (cmd == "WorldEnd") {
      command = pbrt_command_::world_end;
      return true;
    } else if (cmd == "AttributeBegin") {
      command = pbrt_command_::attribute_begin;
      return true;
    } else if (cmd == "AttributeEnd") {
      command = pbrt_command_::attribute_end;
      return true;
    } else if (cmd == "TransformBegin") {
      command = pbrt_command_::transform_begin;
      return true;
    } else if (cmd == "TransformEnd") {
      command = pbrt_command_::transform_end;
      return true;
    } else if (cmd == "ObjectBegin") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::object_begin;
      return true;
    } else if (cmd == "ObjectEnd") {
      command = pbrt_command_::object_end;
      return true;
    } else if (cmd == "ObjectInstance") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::object_instance;
      return true;
    } else if (cmd == "ActiveTransform") {
      parse_pbrt_command(str, name);
      command = pbrt_command_::active_transform;
      return true;
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      parse_pbrt_value(str, xf);
      xform   = frame3f{xf};
      command = pbrt_command_::set_transform;
      return true;
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      parse_pbrt_value(str, xf);
      xform   = frame3f{xf};
      command = pbrt_command_::concat_transform;
      return true;
    } else if (cmd == "Scale") {
      auto v = zero3f;
      parse_pbrt_value(str, v);
      xform   = scaling_frame(v);
      command = pbrt_command_::concat_transform;
      return true;
    } else if (cmd == "Translate") {
      auto v = zero3f;
      parse_pbrt_value(str, v);
      xform   = translation_frame(v);
      command = pbrt_command_::concat_transform;
      return true;
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      parse_pbrt_value(str, v);
      xform   = rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x));
      command = pbrt_command_::concat_transform;
      return true;
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      parse_pbrt_value(str, from);
      parse_pbrt_value(str, to);
      parse_pbrt_value(str, up);
      xform   = {from, to, up, zero3f};
      command = pbrt_command_::lookat_transform;
      return true;
    } else if (cmd == "ReverseOrientation") {
      command = pbrt_command_::reverse_orientation;
      return true;
    } else if (cmd == "CoordinateSystem") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::coordinate_system_set;
      return true;
    } else if (cmd == "CoordSysTransform") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::coordinate_system_transform;
      return true;
    } else if (cmd == "Integrator") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::integrator;
      return true;
    } else if (cmd == "Sampler") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::sampler;
      return true;
    } else if (cmd == "PixelFilter") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::filter;
      return true;
    } else if (cmd == "Film") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::film;
      return true;
    } else if (cmd == "Accelerator") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::accelerator;
      return true;
    } else if (cmd == "Camera") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::camera;
      return true;
    } else if (cmd == "Texture") {
      auto comptype = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_value(str, comptype);
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::named_texture;
      return true;
    } else if (cmd == "Material") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::material;
      return true;
    } else if (cmd == "MakeNamedMaterial") {
      parse_pbrt_value(str, name);
      parse_pbrt_params(str, values);
      type = "";
      for (auto& value : values)
        if (value.name == "type") type = value.value1s;
      command = pbrt_command_::named_material;
      return true;
    } else if (cmd == "NamedMaterial") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::use_material;
      return true;
    } else if (cmd == "Shape") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::shape;
      return true;
    } else if (cmd == "AreaLightSource") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::arealight;
      return true;
    } else if (cmd == "LightSource") {
      parse_pbrt_value(str, type);
      parse_pbrt_params(str, values);
      command = pbrt_command_::light;
      return true;
    } else if (cmd == "MakeNamedMedium") {
      parse_pbrt_value(str, name);
      parse_pbrt_params(str, values);
      type = "";
      for (auto& value : values)
        if (value.name == "type") type = value.value1s;
      command = pbrt_command_::named_medium;
      return true;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      parse_pbrt_value(str, interior);
      parse_pbrt_value(str, exterior);
      name    = interior + "####" + exterior;
      command = pbrt_command_::medium_interface;
      return true;
    } else if (cmd == "Include") {
      parse_pbrt_value(str, name);
      command = pbrt_command_::include;
      return true;
    } else {
      throw std::runtime_error("unknown command " + cmd);
    }
  }
  return false;
}
bool read_pbrt_command(file_wrapper& fs, pbrt_command_& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values) {
  auto command_buffer = ""s;
  return read_pbrt_command(
      fs, command, name, type, xform, values, command_buffer);
}

// Write obj elements
void write_pbrt_comment(file_wrapper& fs, const string& comment) {
  auto lines = split_string(comment, "\n");
  for (auto& line : lines) {
    checked_fprintf(fs, "# %s\n", line.c_str());
  }
  checked_fprintf(fs, "\n");
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
    checked_fprintf(fs, " \"%s %s\" ", type_labels.at(value.type).c_str());
    switch (value.type) {
      case pbrt_value_type::real:
        if (value.vector1f.empty()) {
          checked_fprintf(fs, "[ ");
          for (auto& v : value.vector1f) checked_fprintf(fs, " %g", v);
          checked_fprintf(fs, " ]");
        } else {
          checked_fprintf(fs, "%g", value.value1f);
        }
        break;
      case pbrt_value_type::integer:
        if (value.vector1f.empty()) {
          checked_fprintf(fs, "[ ");
          for (auto& v : value.vector1i) checked_fprintf(fs, " %d", v);
          checked_fprintf(fs, " ]");
        } else {
          checked_fprintf(fs, "%d", value.value1i);
        }
        break;
      case pbrt_value_type::boolean:
        checked_fprintf(fs, "\"%s\"", value.value1b ? "true" : "false");
        break;
      case pbrt_value_type::string:
        checked_fprintf(fs, "\"%s\"", value.value1b ? "true" : "false");
        break;
      case pbrt_value_type::point:
      case pbrt_value_type::vector:
      case pbrt_value_type::normal:
      case pbrt_value_type::color:
        if (!value.vector3f.empty()) {
          checked_fprintf(fs, "[ ");
          for (auto& v : value.vector3f)
            checked_fprintf(fs, " %g %g %g", v.x, v.y, v.z);
          checked_fprintf(fs, " ]");
        } else {
          checked_fprintf(fs, "[ %g %g %g ]", value.value3f.x, value.value3f.y,
              value.value3f.z);
        }
        break;
      case pbrt_value_type::spectrum:
        checked_fprintf(fs, "[ ");
        for (auto& v : value.vector1f) checked_fprintf(fs, " %g", v);
        checked_fprintf(fs, " ]");
        break;
      case pbrt_value_type::texture:
        checked_fprintf(fs, "\"%s\"", value.value1s.c_str());
        break;
      case pbrt_value_type::point2:
      case pbrt_value_type::vector2:
        if (!value.vector2f.empty()) {
          checked_fprintf(fs, "[ ");
          for (auto& v : value.vector2f)
            checked_fprintf(fs, " %g %g", v.x, v.y);
          checked_fprintf(fs, " ]");
        } else {
          checked_fprintf(fs, "[ %g %g ]", value.value2f.x, value.value2f.x);
        }
        break;
    }
  }
  checked_fprintf(fs, "\n");
}

void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_float) {
  switch (command) {
    case pbrt_command_::world_begin: checked_fprintf(fs, "WorldBegin\n"); break;
    case pbrt_command_::world_end: checked_fprintf(fs, "WorldEnd\n"); break;
    case pbrt_command_::attribute_begin:
      checked_fprintf(fs, "AttributeBegin\n");
      break;
    case pbrt_command_::attribute_end:
      checked_fprintf(fs, "AttributeEnd\n");
      break;
    case pbrt_command_::transform_begin:
      checked_fprintf(fs, "TransformBegin\n");
      break;
    case pbrt_command_::transform_end:
      checked_fprintf(fs, "TransformEnd\n");
      break;
    case pbrt_command_::object_begin:
      checked_fprintf(fs, "ObjectBegin \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::object_end: checked_fprintf(fs, "ObjectEnd\n"); break;
    case pbrt_command_::object_instance:
      checked_fprintf(fs, "ObjectInstance \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::sampler:
      checked_fprintf(fs, "Sampler \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::integrator:
      checked_fprintf(fs, "Integrator \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::accelerator:
      checked_fprintf(fs, "Accelerator \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::film:
      checked_fprintf(fs, "Film \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::filter:
      checked_fprintf(fs, "Filter \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::camera:
      checked_fprintf(fs, "Camera \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::shape:
      checked_fprintf(fs, "Shape \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::light:
      checked_fprintf(fs, "Light \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::material:
      checked_fprintf(fs, "Material \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::arealight:
      checked_fprintf(fs, "AreaLight \"%s\"", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::named_texture:
      checked_fprintf(fs, "Texture \"%s\" \"%s\" \"%s\"", name.c_str(),
          texture_float ? "float" : "rgb", type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::named_medium:
      checked_fprintf(fs, "MakeNamedMedium \"%s\" \"string type\" \"%s\"",
          name.c_str(), type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::named_material:
      checked_fprintf(fs, "MakeNamedMaterial \"%s\" \"string type\" \"%s\"",
          name.c_str(), type.c_str());
      write_pbrt_values(fs, values);
      break;
    case pbrt_command_::include:
      checked_fprintf(fs, "Include \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::reverse_orientation:
      checked_fprintf(fs, "ReverseOrientation\n");
      break;
    case pbrt_command_::set_transform:
      checked_fprintf(fs,
          "Transform %g %g %g 0 %g %g %g 0 %g %g %g 0 %g %g %g 1\n", xform.x.x,
          xform.x.y, xform.x.z, xform.y.x, xform.y.y, xform.y.z, xform.z.x,
          xform.z.y, xform.z.z, xform.o.x, xform.o.y, xform.o.z);
      break;
    case pbrt_command_::concat_transform:
      checked_fprintf(fs,
          "ConcatTransform %g %g %g 0 %g %g %g 0 %g %g %g 0 %g %g %g 1\n",
          xform.x.x, xform.x.y, xform.x.z, xform.y.x, xform.y.y, xform.y.z,
          xform.z.x, xform.z.y, xform.z.z, xform.o.x, xform.o.y, xform.o.z);
      break;
    case pbrt_command_::lookat_transform:
      checked_fprintf(fs, "LookAt %g %g %g %g %g %g %g %g %g\n", xform.x.x,
          xform.x.y, xform.x.z, xform.y.x, xform.y.y, xform.y.z, xform.z.x,
          xform.z.y, xform.z.z);
      break;
    case pbrt_command_::use_material:
      checked_fprintf(fs, "NamedMaterial \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::medium_interface: {
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
      checked_fprintf(fs, "MediumInterface \"%s\" \"%s\"\n", interior.c_str(),
          exterior.c_str());
    } break;
    case pbrt_command_::active_transform:
      checked_fprintf(fs, "ActiveTransform \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::coordinate_system_set:
      checked_fprintf(fs, "CoordinateSystem \"%s\"\n", name.c_str());
      break;
    case pbrt_command_::coordinate_system_transform:
      checked_fprintf(fs, "CoordinateSysTransform \"%s\"\n", name.c_str());
      break;
  }
}

void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
    const string& name, const frame3f& xform) {
  return write_pbrt_command(fs, command, name, "", xform, {});
}
void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
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
  pbrt.value1b = value;
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

// old code --- maintained here in case we want to integrate back
#if 0
void approximate_fourier_material(pbrt_material::fourier_t& fourier) {
  auto filename = fs::path(fourier.bsdffile).filename().string();
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
