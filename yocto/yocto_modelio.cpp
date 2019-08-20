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
void open_file(
    file_wrapper& fs, const string& filename, const string& mode) {
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

// Check if a file can be opened for reading.
static inline bool exists_file(const string& filename) {
  auto f = fopen(filename.c_str(), "r");
  if (!f) return false;
  fclose(f);
  return true;
}

// Read a line
static inline bool read_line(file_wrapper& fs, char* buffer, size_t size) {
  return fgets(buffer, size, fs.fs) != nullptr;
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
static inline void write_ply_binprop(file_wrapper& fs, bool big_endian, VT value) {
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

void write_ply_value(file_wrapper& fs, ply_format format, const ply_element& element,
    vector<double>& values, vector<vector<double>>& lists) {
  write_ply_value_impl(fs, format, element, values, lists);
}
void write_ply_value(file_wrapper& fs, ply_format format, const ply_element& element,
    vector<float>& values, vector<vector<int>>& lists) {
  write_ply_value_impl(fs, format, element, values, lists);
}

void read_ply_value(file_wrapper& fs, ply_format format, const ply_element& element,
    vector<double>& values, vector<vector<double>>& lists) {
  read_ply_value_impl(fs, format, element, values, lists);
}
void read_ply_value(file_wrapper& fs, ply_format format, const ply_element& element,
    vector<float>& values, vector<vector<int>>& lists) {
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
static inline void parse_obj_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_obj_value(str, valuei);
  value = (bool)valuei;
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

// Load obj materials
void load_mtl(const string& filename, obj_callbacks& cb, bool fliptr) {
  // open file
  auto fs = open_file(filename);

  // currently parsed material
  auto material = obj_material{};
  auto first    = true;

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
      if (!first) cb.material(material);
      first    = false;
      material = obj_material{};
      parse_obj_value(line, material.name);
    } else if (cmd == "illum") {
      parse_obj_value(line, material.illum);
    } else if (cmd == "Ke") {
      parse_obj_value(line, material.ke);
    } else if (cmd == "Kd") {
      parse_obj_value(line, material.kd);
    } else if (cmd == "Ks") {
      parse_obj_value(line, material.ks);
    } else if (cmd == "Kt") {
      parse_obj_value(line, material.kt);
    } else if (cmd == "Tf") {
      material.kt = {-1, -1, -1};
      parse_obj_value(line, material.kt);
      if (material.kt.y < 0)
        material.kt = {material.kt.x, material.kt.x, material.kt.x};
      if (fliptr) material.kt = vec3f{1, 1, 1} - material.kt;
    } else if (cmd == "Tr") {
      parse_obj_value(line, material.op);
      if (fliptr) material.op = 1 - material.op;
    } else if (cmd == "Ns") {
      parse_obj_value(line, material.ns);
      material.pr = pow(2 / (material.ns + 2), 1 / 4.0f);
      if (material.pr < 0.01f) material.pr = 0;
      if (material.pr > 0.99f) material.pr = 1;
    } else if (cmd == "d") {
      parse_obj_value(line, material.op);
    } else if (cmd == "map_Ke") {
      parse_obj_value(line, material.ke_map);
    } else if (cmd == "map_Kd") {
      parse_obj_value(line, material.kd_map);
    } else if (cmd == "map_Ks") {
      parse_obj_value(line, material.ks_map);
    } else if (cmd == "map_Tr") {
      parse_obj_value(line, material.kt_map);
    } else if (cmd == "map_d" || cmd == "map_Tr") {
      parse_obj_value(line, material.op_map);
    } else if (cmd == "map_bump" || cmd == "bump") {
      parse_obj_value(line, material.bump_map);
    } else if (cmd == "map_occ" || cmd == "occ") {
      parse_obj_value(line, material.occ_map);
    } else if (cmd == "map_disp" || cmd == "disp") {
      parse_obj_value(line, material.disp_map);
    } else if (cmd == "map_norm" || cmd == "norm") {
      parse_obj_value(line, material.norm_map);
    } else if (cmd == "Pm") {
      parse_obj_value(line, material.pm);
    } else if (cmd == "Pr") {
      parse_obj_value(line, material.pr);
    } else if (cmd == "Ps") {
      parse_obj_value(line, material.ps);
    } else if (cmd == "Pc") {
      parse_obj_value(line, material.pc);
    } else if (cmd == "Pcr") {
      parse_obj_value(line, material.pcr);
    } else if (cmd == "map_Pm") {
      parse_obj_value(line, material.pm_map);
    } else if (cmd == "map_Pr") {
      parse_obj_value(line, material.pr_map);
    } else if (cmd == "map_Ps") {
      parse_obj_value(line, material.ps_map);
    } else if (cmd == "map_Pc") {
      parse_obj_value(line, material.pc_map);
    } else if (cmd == "map_Pcr") {
      parse_obj_value(line, material.pcr_map);
    } else if (cmd == "Vt") {
      parse_obj_value(line, material.vt);
    } else if (cmd == "Vp") {
      parse_obj_value(line, material.vp);
    } else if (cmd == "Ve") {
      parse_obj_value(line, material.ve);
    } else if (cmd == "Vs") {
      parse_obj_value(line, material.vs);
    } else if (cmd == "Vg") {
      parse_obj_value(line, material.vg);
    } else if (cmd == "Vr") {
      parse_obj_value(line, material.vr);
    } else if (cmd == "map_Vs") {
      parse_obj_value(line, material.vs_map);
    }
  }

  // issue current material
  if (!first) cb.material(material);
}

// Load obj extensions
void load_objx(const string& filename, obj_callbacks& cb) {
  // open file
  auto fs = open_file(filename);

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
    if (cmd == "c") {
      auto camera = obj_camera();
      parse_obj_value(line, camera.name);
      parse_obj_value(line, camera.ortho);
      parse_obj_value(line, camera.width);
      parse_obj_value(line, camera.height);
      parse_obj_value(line, camera.lens);
      parse_obj_value(line, camera.focus);
      parse_obj_value(line, camera.aperture);
      parse_obj_value(line, camera.frame);
      cb.camera(camera);
    } else if (cmd == "e") {
      auto environment = obj_environment();
      parse_obj_value(line, environment.name);
      parse_obj_value(line, environment.ke);
      parse_obj_value(line, environment.ke_txt.path);
      parse_obj_value(line, environment.frame);
      if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
      cb.environmnet(environment);
    } else if (cmd == "i") {
      auto instance = obj_instance();
      parse_obj_value(line, instance.name);
      parse_obj_value(line, instance.object);
      parse_obj_value(line, instance.material);
      parse_obj_value(line, instance.frame);
      cb.instance(instance);
    } else if (cmd == "po") {
      auto procedural = obj_procedural();
      parse_obj_value(line, procedural.name);
      parse_obj_value(line, procedural.type);
      parse_obj_value(line, procedural.material);
      parse_obj_value(line, procedural.size);
      parse_obj_value(line, procedural.level);
      parse_obj_value(line, procedural.frame);
      cb.procedural(procedural);
    } else {
      // unused
    }
  }
}

// Load obj scene
void load_obj(const string& filename, obj_callbacks& cb, bool nomaterials,
    bool flipv, bool fliptr) {
  // open file
  auto fs = open_file(filename);

  // track vertex size
  auto vert_size = obj_vertex();
  auto verts     = vector<obj_vertex>();  // buffer to avoid reallocation

  // material libraries read already
  auto mlibs = vector<string>{};

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
      auto vert = zero3f;
      parse_obj_value(line, vert);
      cb.vert(vert);
      vert_size.position += 1;
    } else if (cmd == "vn") {
      auto vert = zero3f;
      parse_obj_value(line, vert);
      cb.norm(vert);
      vert_size.normal += 1;
    } else if (cmd == "vt") {
      auto vert = zero2f;
      parse_obj_value(line, vert);
      if (flipv) vert.y = 1 - vert.y;
      cb.texcoord(vert);
      vert_size.texcoord += 1;
    } else if (cmd == "f" || cmd == "l" || cmd == "p") {
      verts.clear();
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
        verts.push_back(vert);
        skip_whitespace(line);
      }
      if (cmd == "f") cb.face(verts);
      if (cmd == "l") cb.line(verts);
      if (cmd == "p") cb.point(verts);
    } else if (cmd == "o") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.object(name);
    } else if (cmd == "usemtl") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.usemtl(name);
    } else if (cmd == "g") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.group(name);
    } else if (cmd == "s") {
      auto name = ""s;
      parse_obj_value_or_empty(line, name);
      cb.smoothing(name);
    } else if (cmd == "mtllib") {
      if (nomaterials) continue;
      auto mtlname = ""s;
      parse_obj_value(line, mtlname);
      cb.mtllib(mtlname);
      if (std::find(mlibs.begin(), mlibs.end(), mtlname) != mlibs.end())
        continue;
      mlibs.push_back(mtlname);
      auto mtlpath = fs::path(filename).parent_path() / mtlname;
      load_mtl(mtlpath, cb, fliptr);
    } else {
      // unused
    }
  }

  // parse extensions if presents
  if (!nomaterials) {
    auto extname    = fs::path(filename).replace_extension(".objx");
    auto ext_exists = exists_file(extname);
    if (ext_exists) {
      load_objx(extname, cb);
    }
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
bool read_objx_command(file_wrapper& fs, objx_command& command, obj_value& value,
    obj_texture_info& texture) {
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
      auto oname = value.string;
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
      if (emission_map.string == "\"\"") emission_map.string = "";
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

void write_obj_command(file_wrapper& fs, obj_command command, const obj_value& value_,
    const vector<obj_vertex>& vertices) {
  auto& name  = value_.string;
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

void write_mtl_command(file_wrapper& fs, mtl_command command, const obj_value& value_,
    const obj_texture_info& texture) {
  auto& name  = value_.string;
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

void write_objx_command(file_wrapper& fs, objx_command command, const obj_value& value_,
    const obj_texture_info& texture) {
  auto& name  = value_.string;
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
  value = yaml.string;
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
// IMPLEMENTATION OF LOW LEVEL PARSING
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
static inline void parse_pbrt_value(string_view& str, bool& value) {
  auto value_name = ""s;
  parse_pbrt_value(str, value_name);
  if (value_name == "true") {
    value = true;
  } else if (value_name == "false") {
    value = false;
  } else {
    throw std::runtime_error("expected boolean");
  }
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
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::bilerp_t::mapping_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::bilerp_t::mapping_type>{
          {"uv", pbrt_texture::bilerp_t::mapping_type::uv},
          {"spherical", pbrt_texture::bilerp_t::mapping_type::spherical},
          {"cylindrical", pbrt_texture::bilerp_t::mapping_type::cylindrical},
          {"planar", pbrt_texture::bilerp_t::mapping_type::planar},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::checkerboard_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::dots_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::imagemap_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::uv_t::mapping_type& value) {
  return parse_pbrt_value(str, (pbrt_texture::bilerp_t::mapping_type&)value);
}

static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::checkerboard_t::aamode_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::checkerboard_t::aamode_type>{
          {"closedform", pbrt_texture::checkerboard_t::aamode_type::closedform},
          {"none", pbrt_texture::checkerboard_t::aamode_type::none},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_texture::imagemap_t::wrap_type& value) {
  static auto value_names =
      unordered_map<string, pbrt_texture::imagemap_t::wrap_type>{
          {"repeat", pbrt_texture::imagemap_t::wrap_type::repeat},
          {"clamp", pbrt_texture::imagemap_t::wrap_type::clamp},
          {"black", pbrt_texture::imagemap_t::wrap_type::black},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_shape::curve_t::basis_t& value) {
  static auto value_names = unordered_map<string, pbrt_shape::curve_t::basis_t>{
      {"bezier", pbrt_shape::curve_t::basis_t::bezier},
      {"bspline", pbrt_shape::curve_t::basis_t::bspline},
  };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_shape::curve_t::type_t& value) {
  static auto value_names = unordered_map<string, pbrt_shape::curve_t::type_t>{
      {"flat", pbrt_shape::curve_t::type_t::flat},
      {"cylinder", pbrt_shape::curve_t::type_t::cylinder},
      {"ribbon", pbrt_shape::curve_t::type_t::ribbon},
  };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_accelerator::bvh_t::splitmethod_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_accelerator::bvh_t::splitmethod_t>{
          {"sah", pbrt_accelerator::bvh_t::splitmethod_t::sah},
          {"equal", pbrt_accelerator::bvh_t::splitmethod_t::equal},
          {"middle", pbrt_accelerator::bvh_t::splitmethod_t::middle},
          {"hlbvh", pbrt_accelerator::bvh_t::splitmethod_t::hlbvh},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::path_t::lightsamplestrategy_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_integrator::path_t::lightsamplestrategy_t>{
          {"power", pbrt_integrator::path_t::lightsamplestrategy_t::power},
          {"spatial", pbrt_integrator::path_t::lightsamplestrategy_t::spatial},
          {"uniform", pbrt_integrator::path_t::lightsamplestrategy_t::uniform},
      };
  return parse_pbrt_value(str, value, value_names);
}
static inline void parse_pbrt_value(string_view&       str,
    pbrt_integrator::volpath_t::lightsamplestrategy_t& value) {
  return parse_pbrt_value(
      str, (pbrt_integrator::path_t::lightsamplestrategy_t&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::bdpt_t::lightsamplestrategy_t& value) {
  return parse_pbrt_value(
      str, (pbrt_integrator::path_t::lightsamplestrategy_t&)value);
}
static inline void parse_pbrt_value(
    string_view& str, pbrt_integrator::directlighting_t::strategy_t& value) {
  static auto value_names =
      unordered_map<string, pbrt_integrator::directlighting_t::strategy_t>{
          {"all", pbrt_integrator::directlighting_t::strategy_t::all},
          {"one", pbrt_integrator::directlighting_t::strategy_t::one},
      };
  return parse_pbrt_value(str, value, value_names);
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
static inline void parse_pbrt_value(string_view& str, vec3i& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, vec4i& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, mat4f& value) {
  for (auto i = 0; i < 4; i++) parse_pbrt_value(str, value[i]);
}
static inline void parse_pbrt_value(string_view& str, pbrt_spectrum3f& value) {
  for (auto i = 0; i < 3; i++) parse_pbrt_value(str, value[i]);
}

// Check next
static inline bool is_pbrt_string(string_view& str) {
  skip_whitespace(str);
  return !str.empty() && str.front() == '"';
}
static inline bool is_open_bracket(string_view& str) {
  skip_whitespace(str);
  return !str.empty() && str.front() == '[';
}
static inline bool is_close_bracket(string_view& str) {
  skip_whitespace(str);
  return !str.empty() && str.front() == ']';
}
static inline bool is_pbrt_param(string_view& str) {
  skip_whitespace(str);
  return is_pbrt_string(str);
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

static inline void skip_pbrt_open_bracket(string_view& str) {
  if (!is_open_bracket(str)) throw std::runtime_error("expected bracket");
  str.remove_prefix(1);
  skip_whitespace(str);
}
static inline void skip_pbrt_close_bracket(string_view& str) {
  if (!is_close_bracket(str)) throw std::runtime_error("expected bracket");
  str.remove_prefix(1);
  skip_whitespace(str);
}

template <typename T>
static inline void parse_pbrt_param(string_view& str, T& value) {
  auto has_brackets = is_open_bracket(str);
  if (has_brackets) skip_pbrt_open_bracket(str);
  parse_pbrt_value(str, value);
  if (has_brackets) skip_pbrt_close_bracket(str);
}

template <typename T>
static inline void parse_pbrt_param(string_view& str, vector<T>& values) {
  skip_pbrt_open_bracket(str);
  values.clear();
  while (!is_close_bracket(str)) {
    values.push_back({});
    parse_pbrt_value(str, values.back());
  }
  skip_pbrt_close_bracket(str);
}

template <typename T>
static inline bool is_pbrt_type_compatible(const string& type) {
  if constexpr (std::is_same<T, int>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, float>::value) {
    return type == "float";
  } else if constexpr (std::is_same<T, bool>::value) {
    return type == "bool";
  } else if constexpr (std::is_same<T, string>::value) {
    return type == "string";
  } else if constexpr (std::is_same<T, vec2f>::value) {
    return type == "point2" || type == "vector2" || type == "float";
  } else if constexpr (std::is_same<T, vec3f>::value) {
    return type == "point3" || type == "vector3" || type == "normal3" ||
           type == "point" || type == "vector" || type == "normal" ||
           type == "float";
  } else if constexpr (std::is_same<T, vec4f>::value) {
    return type == "float";
  } else if constexpr (std::is_same<T, pbrt_spectrum3f>::value) {
    return type == "rgb" || type == "pbrt_spectrum" || type == "blackbody";
  } else if constexpr (std::is_same<T, vec3i>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, vec4i>::value) {
    return type == "integer";
  } else if constexpr (std::is_same<T, bbox2f>::value) {
    return type == "float";
  } else if constexpr (std::is_enum<T>::value) {
    return type == "string";
  } else {
    return false;
  }
}

template <typename T>
static inline void parse_pbrt_param(
    string_view& str, const string& type, T& value) {
  if (!is_pbrt_type_compatible<T>(type)) {
    throw std::runtime_error("incompatible type " + type);
  }
  parse_pbrt_param(str, value);
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

static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_spectrum3f& value) {
  bool verbose = false;
  if (type == "rgb") {
    parse_pbrt_param(str, value);
  } else if (type == "color") {
    parse_pbrt_param(str, value);
  } else if (type == "float") {
    auto valuef = 0.0f;
    parse_pbrt_param(str, valuef);
    value = {valuef, valuef, valuef};
  } else if (type == "blackbody") {
    auto blackbody = zero2f;
    parse_pbrt_param(str, blackbody);
    (vec3f&)value = blackbody_to_rgb(blackbody.x) * blackbody.y;
  } else if (type == "spectrum" && is_pbrt_string(str)) {
    if (verbose) printf("spectrum  not well supported\n");
    auto filename = ""s;
    parse_pbrt_param(str, filename);
    auto filenamep = fs::path(filename).filename();
    if (filenamep.extension() == ".spd") {
      filenamep = filenamep.replace_extension("");
      if (filenamep == "SHPS") {
        value = {1, 1, 1};
      } else if (filenamep.extension() == ".eta") {
        auto eta = get_pbrt_etak(filenamep.replace_extension("")).first;
        value    = {eta.x, eta.y, eta.z};
      } else if (filenamep.extension() == ".k") {
        auto k = get_pbrt_etak(filenamep.replace_extension("")).second;
        value  = {k.x, k.y, k.z};
      } else {
        throw std::runtime_error("unknown spectrum file " + filename);
      }
    } else {
      throw std::runtime_error("unsupported spectrum format");
      // value = {1, 0, 0};
    }
  } else if (type == "spectrum" && !is_pbrt_string(str)) {
    if (verbose) printf("spectrum  not well supported\n");
    auto values = vector<float>{};
    parse_pbrt_param(str, values);
    value = {1, 0, 0};
  } else {
    throw std::runtime_error("unsupported spectrum type");
  }
}

template <typename T>
static inline void parse_pbrt_param(
    string_view& str, const string& type, vector<T>& value) {
  if (!is_pbrt_type_compatible<T>(type)) {
    throw std::runtime_error("incompatible type " + type);
  }
  parse_pbrt_param(str, value);
}

static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_textured3f& value) {
  if (type == "texture") {
    parse_pbrt_param(str, value.texture);
  } else {
    parse_pbrt_param(str, type, value.value);
  }
}
static inline void parse_pbrt_param(
    string_view& str, const string& type, pbrt_textured1f& value) {
  if (type == "texture") {
    parse_pbrt_param(str, value.texture);
  } else {
    parse_pbrt_param(str, type, value.value);
  }
}

static inline void skip_pbrt_value(string_view& str) {
  skip_whitespace(str);
  if (str.front() == '"') {
    str.remove_prefix(1);
    str.remove_prefix(str.find('"') + 1);
  } else {
    str.remove_prefix(str.find_first_of(" \n\t\r],\""));
  }
  skip_whitespace(str);
}

static inline void skip_pbrt_param(string_view& str) {
  if (is_open_bracket(str)) {
    skip_pbrt_open_bracket(str);
    while (!is_close_bracket(str)) skip_pbrt_value(str);
    skip_pbrt_close_bracket(str);
  } else {
    skip_pbrt_value(str);
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
  value = yaml.string;
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
    parse_yaml_value(str, value.string);
    if (value.string == "true" || value.string == "false") {
      value.type    = yaml_value_type::boolean;
      value.boolean = value.string == "true";
    }
  }
  skip_whitespace(str);
  if (!str.empty() && !is_whitespace(str)) throw std::runtime_error("bad yaml");
}

bool read_yaml_property(
    file_wrapper& fs, string& group, string& key, bool& newobj, yaml_value& value) {
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
void write_yaml_property(file_wrapper& fs, const string& object, const string& key,
    bool newobj, const yaml_value& value) {
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
        checked_fprintf(fs, "%s", value.string.c_str());
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

// Parse Accelerator
static inline void parse_pbrt_accelerator(
    string_view& str, const string& type, pbrt_accelerator& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "bvh") {
    auto tvalue = pbrt_accelerator::bvh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxnodeprims") {
        parse_pbrt_param(str, ptype, tvalue.maxnodeprims);
      } else if (pname == "splitmethod") {
        parse_pbrt_param(str, ptype, tvalue.splitmethod);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_accelerator::type_t::bvh;
    value.bvh  = tvalue;
  } else if (type == "kdtree") {
    auto tvalue = pbrt_accelerator::kdtree_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "intersectcost") {
        parse_pbrt_param(str, ptype, tvalue.intersectcost);
      } else if (pname == "traversalcost") {
        parse_pbrt_param(str, ptype, tvalue.traversalcost);
      } else if (pname == "emptybonus") {
        parse_pbrt_param(str, ptype, tvalue.emptybonus);
      } else if (pname == "maxprims") {
        parse_pbrt_param(str, ptype, tvalue.maxprims);
      } else if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_accelerator::type_t::kdtree;
    value.kdtree = tvalue;
  } else {
    throw std::runtime_error("unknown Accelerator " + type);
  }
}

// Parse Integrator
static inline void parse_pbrt_integrator(
    string_view& str, const string& type, pbrt_integrator& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "path") {
    auto tvalue = pbrt_integrator::path_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "rrthreshold") {
        parse_pbrt_param(str, ptype, tvalue.rrthreshold);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
      // parse_pbrt_optional_param(str, "lightsamplestrategy",
      // tvalue.lightsamplestrategy); // TODO: enums
    }
    value.type = pbrt_integrator::type_t::path;
    value.path = tvalue;
  } else if (type == "volpath") {
    auto tvalue = pbrt_integrator::volpath_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "rrthreshold") {
        parse_pbrt_param(str, ptype, tvalue.rrthreshold);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_integrator::type_t::volpath;
    value.volpath = tvalue;
  } else if (type == "directlighting") {
    auto tvalue = pbrt_integrator::directlighting_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "strategy") {
        parse_pbrt_param(str, ptype, tvalue.strategy);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type           = pbrt_integrator::type_t::directlighting;
    value.directlighting = tvalue;
  } else if (type == "bdpt") {
    auto tvalue = pbrt_integrator::bdpt_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "lightsamplestrategy") {
        parse_pbrt_param(str, ptype, tvalue.lightsamplestrategy);
      } else if (pname == "visualizestrategies") {
        parse_pbrt_param(str, ptype, tvalue.visualizestrategies);
      } else if (pname == "visualizeweights") {
        parse_pbrt_param(str, ptype, tvalue.visualizeweights);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::bdpt;
    value.bdpt = tvalue;
  } else if (type == "mlt") {
    auto tvalue = pbrt_integrator::mlt_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "bootstrapsamples") {
        parse_pbrt_param(str, ptype, tvalue.bootstrapsamples);
      } else if (pname == "chains") {
        parse_pbrt_param(str, ptype, tvalue.chains);
      } else if (pname == "mutationsperpixel") {
        parse_pbrt_param(str, ptype, tvalue.mutationsperpixel);
      } else if (pname == "largestepprobability") {
        parse_pbrt_param(str, ptype, tvalue.largestepprobability);
      } else if (pname == "sigma") {
        parse_pbrt_param(str, ptype, tvalue.sigma);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::mlt;
    value.mlt  = tvalue;
  } else if (type == "sppm") {
    auto tvalue = pbrt_integrator::sppm_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else if (pname == "iterations") {
        parse_pbrt_param(str, ptype, tvalue.iterations);
      } else if (pname == "numiterations") {
        parse_pbrt_param(str, ptype, tvalue.iterations);
      } else if (pname == "photonsperiteration") {
        parse_pbrt_param(str, ptype, tvalue.photonsperiteration);
      } else if (pname == "imagewritefrequency") {
        parse_pbrt_param(str, ptype, tvalue.imagewritefrequency);
      } else if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_integrator::type_t::sppm;
    value.sppm = tvalue;
  } else if (type == "whitted") {
    auto tvalue = pbrt_integrator::whitted_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "maxdepth") {
        parse_pbrt_param(str, ptype, tvalue.maxdepth);
      } else if (pname == "pixelbounds") {
        parse_pbrt_param(str, ptype, tvalue.pixelbounds);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_integrator::type_t::whitted;
    value.whitted = tvalue;
  } else {
    throw std::runtime_error("unknown Integrator " + type);
  }
}

// Parse Sampler
static inline void parse_pbrt_sampler(
    string_view& str, const string& type, pbrt_sampler& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "random") {
    auto tvalue = pbrt_sampler::random_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_sampler::type_t::random;
    value.random = tvalue;
  } else if (type == "halton") {
    auto tvalue = pbrt_sampler::halton_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_sampler::type_t::halton;
    value.halton = tvalue;
  } else if (type == "sobol") {
    auto tvalue = pbrt_sampler::sobol_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_sampler::type_t::sobol;
    value.sobol = tvalue;
  } else if (type == "02sequence") {
    auto tvalue = pbrt_sampler::zerotwosequence_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type            = pbrt_sampler::type_t::zerotwosequence;
    value.zerotwosequence = tvalue;
  } else if (type == "lowdiscrepancy") {
    auto tvalue = pbrt_sampler::zerotwosequence_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type            = pbrt_sampler::type_t::zerotwosequence;
    value.zerotwosequence = tvalue;
  } else if (type == "maxmindist") {
    auto tvalue = pbrt_sampler::maxmindist_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "pixelsamples") {
        parse_pbrt_param(str, ptype, tvalue.pixelsamples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_sampler::type_t::maxmindist;
    value.maxmindist = tvalue;
  } else if (type == "stratified") {
    auto tvalue = pbrt_sampler::stratified_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xsamples") {
        parse_pbrt_param(str, ptype, tvalue.xsamples);
      } else if (pname == "ysamples") {
        parse_pbrt_param(str, ptype, tvalue.ysamples);
      } else if (pname == "jitter") {
        parse_pbrt_param(str, ptype, tvalue.jitter);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_sampler::type_t::stratified;
    value.stratified = tvalue;
  } else {
    throw std::runtime_error("unknown Sampler " + type);
  }
}

// Parse Filter
static inline void parse_pbrt_filter(
    string_view& str, const string& type, pbrt_filter& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "box") {
    auto tvalue = pbrt_filter::box_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_filter::type_t::box;
    value.box  = tvalue;
  } else if (type == "gaussian") {
    auto tvalue = pbrt_filter::gaussian_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::gaussian;
    value.gaussian = tvalue;
  } else if (type == "mitchell") {
    auto tvalue = pbrt_filter::mitchell_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "B") {
        parse_pbrt_param(str, ptype, tvalue.B);
      } else if (pname == "C") {
        parse_pbrt_param(str, ptype, tvalue.C);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::mitchell;
    value.mitchell = tvalue;
  } else if (type == "sinc") {
    auto tvalue = pbrt_filter::sinc_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else if (pname == "tau") {
        parse_pbrt_param(str, ptype, tvalue.tau);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_filter::type_t::sinc;
    value.sinc = tvalue;
  } else if (type == "triangle") {
    auto tvalue = pbrt_filter::triangle_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xwidth") {
        parse_pbrt_param(str, ptype, tvalue.xwidth);
      } else if (pname == "ywidth") {
        parse_pbrt_param(str, ptype, tvalue.ywidth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_filter::type_t::triangle;
    value.triangle = tvalue;
  } else {
    throw std::runtime_error("unknown PixelFilter " + type);
  }
}

// Parse Filter
static inline void parse_pbrt_film(
    string_view& str, const string& type, pbrt_film& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "image") {
    auto tvalue = pbrt_film::image_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "xresolution") {
        parse_pbrt_param(str, ptype, tvalue.xresolution);
      } else if (pname == "yresolution") {
        parse_pbrt_param(str, ptype, tvalue.yresolution);
      } else if (pname == "yresolution") {
        parse_pbrt_param(str, ptype, tvalue.yresolution);
      } else if (pname == "cropwindow") {
        parse_pbrt_param(str, ptype, tvalue.cropwindow);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "maxsampleluminance") {
        parse_pbrt_param(str, ptype, tvalue.maxsampleluminance);
      } else if (pname == "diagonal") {
        parse_pbrt_param(str, ptype, tvalue.diagonal);
      } else if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_film::type_t::image;
    value.image = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse Camera
static inline void parse_pbrt_camera(
    string_view& str, const string& type, pbrt_camera& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "perspective") {
    auto tvalue = pbrt_camera::perspective_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "fov") {
        parse_pbrt_param(str, ptype, tvalue.fov);
      } else if (pname == "frameaspectratio") {
        parse_pbrt_param(str, ptype, tvalue.frameaspectratio);
      } else if (pname == "lensradius") {
        parse_pbrt_param(str, ptype, tvalue.lensradius);
      } else if (pname == "focaldistance") {
        parse_pbrt_param(str, ptype, tvalue.focaldistance);
      } else if (pname == "screenwindow") {
        parse_pbrt_param(str, ptype, tvalue.screenwindow);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_camera::type_t::perspective;
    value.perspective = tvalue;
  } else if (type == "orthographic") {
    auto tvalue = pbrt_camera::orthographic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "frameaspectratio") {
        parse_pbrt_param(str, ptype, tvalue.frameaspectratio);
      } else if (pname == "lensradius") {
        parse_pbrt_param(str, ptype, tvalue.lensradius);
      } else if (pname == "focaldistance") {
        parse_pbrt_param(str, ptype, tvalue.focaldistance);
      } else if (pname == "screenwindow") {
        parse_pbrt_param(str, ptype, tvalue.screenwindow);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_camera::type_t::orthographic;
    value.orthographic = tvalue;
  } else if (type == "environment") {
    auto tvalue = pbrt_camera::environment_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_camera::type_t::environment;
    value.environment = tvalue;
  } else if (type == "realistic") {
    auto tvalue = pbrt_camera::realistic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "lensfile") {
        parse_pbrt_param(str, ptype, tvalue.lensfile);
        // example: wide.22mm.dat
        auto lensfile = fs::path(tvalue.lensfile).filename().string();
        lensfile      = lensfile.substr(0, lensfile.size() - 4);
        lensfile      = lensfile.substr(lensfile.find('.') + 1);
        lensfile      = lensfile.substr(0, lensfile.size() - 2);
        tvalue.approx_focallength = std::atof(lensfile.c_str());
      } else if (pname == "aperturediameter") {
        parse_pbrt_param(str, ptype, tvalue.aperturediameter);
      } else if (pname == "focusdistance") {
        parse_pbrt_param(str, ptype, tvalue.focusdistance);
      } else if (pname == "simpleweighting") {
        parse_pbrt_param(str, ptype, tvalue.simpleweighting);
      } else if (pname == "shutteropen") {
        parse_pbrt_param(str, ptype, tvalue.shutteropen);
      } else if (pname == "shutterclose") {
        parse_pbrt_param(str, ptype, tvalue.shutterclose);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type      = pbrt_camera::type_t::realistic;
    value.realistic = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse Texture
static inline void parse_pbrt_texture(
    string_view& str, const string& type, pbrt_texture& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "constant") {
    auto tvalue = pbrt_texture::constant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "value") {
        parse_pbrt_param(str, ptype, tvalue.value);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::constant;
    value.constant = tvalue;
  } else if (type == "bilerp") {
    auto tvalue = pbrt_texture::bilerp_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "v00") {
        parse_pbrt_param(str, ptype, tvalue.v00);
      } else if (pname == "v01") {
        parse_pbrt_param(str, ptype, tvalue.v01);
      } else if (pname == "v10") {
        parse_pbrt_param(str, ptype, tvalue.v10);
      } else if (pname == "v11") {
        parse_pbrt_param(str, ptype, tvalue.v11);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_texture::type_t::bilerp;
    value.bilerp = tvalue;
  } else if (type == "checkerboard") {
    auto tvalue = pbrt_texture::checkerboard_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "dimension") {
        parse_pbrt_param(str, ptype, tvalue.dimension);
      } else if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else if (pname == "aamode") {
        parse_pbrt_param(str, ptype, tvalue.aamode);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_texture::type_t::checkerboard;
    value.checkerboard = tvalue;
  } else if (type == "dots") {
    auto tvalue = pbrt_texture::dots_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "inside") {
        parse_pbrt_param(str, ptype, tvalue.inside);
      } else if (pname == "outside") {
        parse_pbrt_param(str, ptype, tvalue.outside);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::dots;
    value.dots = tvalue;
  } else if (type == "imagemap") {
    auto tvalue = pbrt_texture::imagemap_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else if (pname == "wrap") {
        parse_pbrt_param(str, ptype, tvalue.wrap);
      } else if (pname == "maxanisotropy") {
        parse_pbrt_param(str, ptype, tvalue.maxanisotropy);
      } else if (pname == "trilinear") {
        parse_pbrt_param(str, ptype, tvalue.trilinear);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "gamma") {
        parse_pbrt_param(str, ptype, tvalue.gamma);
      } else if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::imagemap;
    value.imagemap = tvalue;
  } else if (type == "mix") {
    auto tvalue = pbrt_texture::mix_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else if (pname == "amount") {
        parse_pbrt_param(str, ptype, tvalue.amount);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::mix;
    value.mix  = tvalue;
  } else if (type == "scale") {
    auto tvalue = pbrt_texture::scale_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "tex1") {
        parse_pbrt_param(str, ptype, tvalue.tex1);
      } else if (pname == "tex2") {
        parse_pbrt_param(str, ptype, tvalue.tex2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_texture::type_t::scale;
    value.scale = tvalue;
  } else if (type == "fbm") {
    auto tvalue = pbrt_texture::fbm_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::fbm;
    value.fbm  = tvalue;
  } else if (type == "wrinkled") {
    auto tvalue = pbrt_texture::wrinkled_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_texture::type_t::wrinkled;
    value.wrinkled = tvalue;
  } else if (type == "windy") {
    auto tvalue = pbrt_texture::windy_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "") {
        // TODO: missing params
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_texture::type_t::windy;
    value.windy = tvalue;
  } else if (type == "marble") {
    auto tvalue = pbrt_texture::marble_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "octaves") {
        parse_pbrt_param(str, ptype, tvalue.octaves);
      } else if (pname == "roughness") {
        parse_pbrt_param(str, ptype, tvalue.roughness);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "variation") {
        parse_pbrt_param(str, ptype, tvalue.variation);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_texture::type_t::marble;
    value.marble = tvalue;
  } else if (type == "uv") {
    auto tvalue = pbrt_texture::uv_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "mapping") {
        parse_pbrt_param(str, ptype, tvalue.mapping);
      } else if (pname == "uscale") {
        parse_pbrt_param(str, ptype, tvalue.uscale);
      } else if (pname == "vscale") {
        parse_pbrt_param(str, ptype, tvalue.vscale);
      } else if (pname == "udelta") {
        parse_pbrt_param(str, ptype, tvalue.udelta);
      } else if (pname == "vdelta") {
        parse_pbrt_param(str, ptype, tvalue.vdelta);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "v2") {
        parse_pbrt_param(str, ptype, tvalue.v2);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_texture::type_t::uv;
    value.uv   = tvalue;
  } else {
    throw std::runtime_error("unknown Texture " + type);
  }
}

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

// Get typename
static inline void parse_pbrt_typeparam(string_view& str, string& value) {
  auto saved = str;
  value      = "";
  auto pname = ""s, ptype = ""s;
  while (is_pbrt_param(str) && value == "") {
    parse_pbrt_nametype(str, pname, ptype);
    if (pname == "type") {
      parse_pbrt_param(str, ptype, value);
    } else {
      skip_pbrt_param(str);
    }
  }
  if (value == "") throw std::runtime_error("type not found");
  str = saved;
}

// Parse param and resolve constant textures
static inline void parse_pbrt_texture(string_view& str, const string& ptype,
    pbrt_textured3f&                              value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  parse_pbrt_param(str, ptype, value);
  if (value.texture == "") return;
  if (constant_values.find(value.texture) == constant_values.end()) return;
  value.value   = constant_values.at(value.texture);
  value.texture = "";
}
static inline void parse_pbrt_texture(string_view& str, const string& ptype,
    pbrt_textured1f&                              value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  parse_pbrt_param(str, ptype, value);
  if (value.texture == "") return;
  if (constant_values.find(value.texture) == constant_values.end()) return;
  auto col      = constant_values.at(value.texture);
  value.value   = (col.x + col.y + col.z) / 3;
  value.texture = "";
}

// Parse Material
static inline void parse_pbrt_material(string_view& str, const string& type,
    pbrt_material&                                value,
    const unordered_map<string, pbrt_spectrum3f>& constant_values) {
  auto pname = ""s, ptype = ""s;
  if (type == "matte") {
    auto tvalue = pbrt_material::matte_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "sigma") {
        parse_pbrt_param(str, ptype, tvalue.sigma);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::matte;
    value.matte = tvalue;
  } else if (type == "mirror") {
    auto tvalue = pbrt_material::mirror_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_material::type_t::mirror;
    value.mirror = tvalue;
  } else if (type == "plastic") {
    auto tvalue = pbrt_material::plastic_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_material::type_t::plastic;
    value.plastic = tvalue;
  } else if (type == "metal") {
    auto tvalue = pbrt_material::metal_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "k") {
        parse_pbrt_texture(str, ptype, tvalue.k, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::metal;
    value.metal = tvalue;
  } else if (type == "glass") {
    auto tvalue = pbrt_material::glass_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_material::type_t::glass;
    value.glass = tvalue;
  } else if (type == "translucent") {
    auto tvalue = pbrt_material::translucent_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "reflect") {
        parse_pbrt_texture(str, ptype, tvalue.reflect, constant_values);
      } else if (pname == "transmit") {
        parse_pbrt_texture(str, ptype, tvalue.transmit, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_material::type_t::translucent;
    value.translucent = tvalue;
  } else if (type == "uber") {
    auto tvalue = pbrt_material::uber_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "opacity") {
        parse_pbrt_texture(str, ptype, tvalue.opacity, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::uber;
    value.uber = tvalue;
  } else if (type == "disney") {
    auto tvalue = pbrt_material::disney_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "color") {
        parse_pbrt_texture(str, ptype, tvalue.color, constant_values);
      } else if (pname == "anisotropic") {
        parse_pbrt_texture(str, ptype, tvalue.anisotropic, constant_values);
      } else if (pname == "clearcoat") {
        parse_pbrt_texture(str, ptype, tvalue.clearcoat, constant_values);
      } else if (pname == "clearcoatgloss") {
        parse_pbrt_texture(str, ptype, tvalue.clearcoatgloss, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "metallic") {
        parse_pbrt_texture(str, ptype, tvalue.metallic, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "scatterdistance") {
        parse_pbrt_texture(str, ptype, tvalue.scatterdistance, constant_values);
      } else if (pname == "sheen") {
        parse_pbrt_texture(str, ptype, tvalue.sheen, constant_values);
      } else if (pname == "sheentint") {
        parse_pbrt_texture(str, ptype, tvalue.sheentint, constant_values);
      } else if (pname == "spectrans") {
        parse_pbrt_texture(str, ptype, tvalue.spectrans, constant_values);
      } else if (pname == "thin") {
        parse_pbrt_param(str, ptype, tvalue.thin);
      } else if (pname == "difftrans") {
        parse_pbrt_texture(str, ptype, tvalue.difftrans, constant_values);
      } else if (pname == "flatness") {
        parse_pbrt_texture(str, ptype, tvalue.flatness, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_material::type_t::disney;
    value.disney = tvalue;
  } else if (type == "hair") {
    auto tvalue = pbrt_material::hair_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "color") {
        parse_pbrt_texture(str, ptype, tvalue.color, constant_values);
      } else if (pname == "sigma_a") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_a, constant_values);
      } else if (pname == "eumelanin") {
        parse_pbrt_texture(str, ptype, tvalue.eumelanin, constant_values);
      } else if (pname == "pheomelanin") {
        parse_pbrt_texture(str, ptype, tvalue.pheomelanin, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "beta_m") {
        parse_pbrt_texture(str, ptype, tvalue.beta_m, constant_values);
      } else if (pname == "beta_n") {
        parse_pbrt_texture(str, ptype, tvalue.beta_n, constant_values);
      } else if (pname == "alpha") {
        parse_pbrt_texture(str, ptype, tvalue.alpha, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::hair;
    value.hair = tvalue;
  } else if (type == "kdsubsurface") {
    auto tvalue = pbrt_material::kdsubsurface_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "mfp") {
        parse_pbrt_texture(str, ptype, tvalue.mfp, constant_values);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_material::type_t::kdsubsurface;
    value.kdsubsurface = tvalue;
  } else if (type == "mix") {
    auto tvalue = pbrt_material::mix_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "amount") {
        parse_pbrt_texture(str, ptype, tvalue.amount, constant_values);
      } else if (pname == "namedmaterial1") {
        parse_pbrt_param(str, ptype, tvalue.namedmaterial1);
      } else if (pname == "namedmaterial2") {
        parse_pbrt_param(str, ptype, tvalue.namedmaterial2);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_material::type_t::mix;
    value.mix  = tvalue;
  } else if (type == "fourier") {
    auto tvalue = pbrt_material::fourier_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "bsdffile") {
        parse_pbrt_param(str, ptype, tvalue.bsdffile);
        approximate_fourier_material(tvalue);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_material::type_t::fourier;
    value.fourier = tvalue;
  } else if (type == "substrate") {
    auto tvalue = pbrt_material::substrate_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "Kd") {
        parse_pbrt_texture(str, ptype, tvalue.Kd, constant_values);
      } else if (pname == "Ks") {
        parse_pbrt_texture(str, ptype, tvalue.Ks, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type      = pbrt_material::type_t::substrate;
    value.substrate = tvalue;
  } else if (type == "subsurface") {
    auto tvalue = pbrt_material::subsurface_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "name") {
        parse_pbrt_param(str, ptype, tvalue.name);
        auto params    = parse_pbrt_subsurface(tvalue.name);
        tvalue.sigma_a = {params.second.x, params.second.y, params.second.z};
        tvalue.sigma_prime_s = {params.first.x, params.first.y, params.first.z};
      } else if (pname == "sigma_a") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_a, constant_values);
      } else if (pname == "sigma_prime_s") {
        parse_pbrt_texture(str, ptype, tvalue.sigma_prime_s, constant_values);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "eta") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "index") {
        parse_pbrt_texture(str, ptype, tvalue.eta, constant_values);
      } else if (pname == "Kr") {
        parse_pbrt_texture(str, ptype, tvalue.Kr, constant_values);
      } else if (pname == "Kt") {
        parse_pbrt_texture(str, ptype, tvalue.Kt, constant_values);
      } else if (pname == "roughness") {
        pbrt_textured1f roughness = 0.01f;
        parse_pbrt_param(str, ptype, roughness);
        tvalue.uroughness = roughness;
        tvalue.vroughness = roughness;
      } else if (pname == "uroughness") {
        parse_pbrt_texture(str, ptype, tvalue.uroughness, constant_values);
      } else if (pname == "vroughness") {
        parse_pbrt_texture(str, ptype, tvalue.vroughness, constant_values);
      } else if (pname == "remaproughness") {
        parse_pbrt_param(str, ptype, tvalue.remaproughness);
      } else if (pname == "bumpmap") {
        parse_pbrt_texture(str, ptype, tvalue.bumpmap, constant_values);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_material::type_t::subsurface;
    value.subsurface = tvalue;
  } else {
    throw std::runtime_error("unknown Material " + type);
  }
}

// Parse Shape
static inline void parse_pbrt_shape(
    string_view& str, const string& type, pbrt_shape& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "trianglemesh") {
    auto tvalue = pbrt_shape::trianglemesh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "indices") {
        parse_pbrt_param(str, ptype, tvalue.indices);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "N") {
        parse_pbrt_param(str, ptype, tvalue.N);
      } else if (pname == "S") {
        parse_pbrt_param(str, ptype, tvalue.S);
      } else if (pname == "uv") {
        parse_pbrt_param(str, ptype, tvalue.uv);
      } else if (pname == "st") {
        parse_pbrt_param(str, ptype, tvalue.uv);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else if (pname == "shadowalpha") {
        parse_pbrt_param(str, ptype, tvalue.shadowalpha);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type         = pbrt_shape::type_t::trianglemesh;
    value.trianglemesh = tvalue;
  } else if (type == "plymesh") {
    auto tvalue = pbrt_shape::plymesh_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "filename") {
        parse_pbrt_param(str, ptype, tvalue.filename);
      } else if (pname == "alpha") {
        parse_pbrt_param(str, ptype, tvalue.alpha);
      } else if (pname == "shadowalpha") {
        parse_pbrt_param(str, ptype, tvalue.shadowalpha);
      } else if (pname == "discarddegenerateUVs") {
        // hack for some files
        auto value = false;
        parse_pbrt_param(str, ptype, value);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_shape::type_t::plymesh;
    value.plymesh = tvalue;
  } else if (type == "curve") {
    auto tvalue = pbrt_shape::curve_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "N") {
        parse_pbrt_param(str, ptype, tvalue.N);
      } else if (pname == "basis") {
        parse_pbrt_param(str, ptype, tvalue.basis);
      } else if (pname == "degree") {
        parse_pbrt_param(str, ptype, tvalue.degree);
      } else if (pname == "type") {
        parse_pbrt_param(str, ptype, tvalue.type);
      } else if (pname == "width") {
        auto width = 1.0f;
        parse_pbrt_param(str, ptype, width);
        tvalue.width0 = width;
        tvalue.width1 = width;
      } else if (pname == "width0") {
        parse_pbrt_param(str, ptype, tvalue.width0);
      } else if (pname == "width1") {
        parse_pbrt_param(str, ptype, tvalue.width1);
      } else if (pname == "splitdepth") {
        parse_pbrt_param(str, ptype, tvalue.splitdepth);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_shape::type_t::curve;
    value.curve = tvalue;
  } else if (type == "loopsubdiv") {
    auto tvalue = pbrt_shape::loopsubdiv_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "indices") {
        parse_pbrt_param(str, ptype, tvalue.indices);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "levels") {
        parse_pbrt_param(str, ptype, tvalue.levels);
      } else if (pname == "nlevels") {
        parse_pbrt_param(str, ptype, tvalue.levels);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_shape::type_t::loopsubdiv;
    value.loopsubdiv = tvalue;
  } else if (type == "nurbs") {
    auto tvalue = pbrt_shape::nurbs_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "nu") {
        parse_pbrt_param(str, ptype, tvalue.nu);
      } else if (pname == "nv") {
        parse_pbrt_param(str, ptype, tvalue.nv);
      } else if (pname == "uknots") {
        parse_pbrt_param(str, ptype, tvalue.uknots);
      } else if (pname == "vknots") {
        parse_pbrt_param(str, ptype, tvalue.vknots);
      } else if (pname == "u0") {
        parse_pbrt_param(str, ptype, tvalue.u0);
      } else if (pname == "v0") {
        parse_pbrt_param(str, ptype, tvalue.v0);
      } else if (pname == "u1") {
        parse_pbrt_param(str, ptype, tvalue.u1);
      } else if (pname == "v1") {
        parse_pbrt_param(str, ptype, tvalue.v1);
      } else if (pname == "P") {
        parse_pbrt_param(str, ptype, tvalue.P);
      } else if (pname == "Pw") {
        parse_pbrt_param(str, ptype, tvalue.Pw);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_shape::type_t::nurbs;
    value.nurbs = tvalue;
  } else if (type == "sphere") {
    auto tvalue = pbrt_shape::sphere_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type   = pbrt_shape::type_t::sphere;
    value.sphere = tvalue;
  } else if (type == "disk") {
    auto tvalue = pbrt_shape::disk_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "height") {
        parse_pbrt_param(str, ptype, tvalue.height);
      } else if (pname == "innerradius") {
        parse_pbrt_param(str, ptype, tvalue.innerradius);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_shape::type_t::disk;
    value.disk = tvalue;
  } else if (type == "cone") {
    auto tvalue = pbrt_shape::cone_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "height") {
        parse_pbrt_param(str, ptype, tvalue.height);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_shape::type_t::cone;
    value.cone = tvalue;
  } else if (type == "cylinder") {
    auto tvalue = pbrt_shape::cylinder_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_shape::type_t::cylinder;
    value.cylinder = tvalue;
  } else if (type == "hyperboloid") {
    auto tvalue = pbrt_shape::hyperboloid_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "p1") {
        parse_pbrt_param(str, ptype, tvalue.p1);
      } else if (pname == "p2") {
        parse_pbrt_param(str, ptype, tvalue.p2);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_shape::type_t::hyperboloid;
    value.hyperboloid = tvalue;
  } else if (type == "paraboloid") {
    auto tvalue = pbrt_shape::paraboloid_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "radius") {
        parse_pbrt_param(str, ptype, tvalue.radius);
      } else if (pname == "zmin") {
        parse_pbrt_param(str, ptype, tvalue.zmin);
      } else if (pname == "zmax") {
        parse_pbrt_param(str, ptype, tvalue.zmax);
      } else if (pname == "phimax") {
        parse_pbrt_param(str, ptype, tvalue.phimax);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_shape::type_t::paraboloid;
    value.paraboloid = tvalue;
  } else if (type == "heightfield") {
    auto tvalue = pbrt_shape::heightfield_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "nu") {
        parse_pbrt_param(str, ptype, tvalue.nu);
      } else if (pname == "nv") {
        parse_pbrt_param(str, ptype, tvalue.nv);
      } else if (pname == "Pz") {
        parse_pbrt_param(str, ptype, tvalue.Pz);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_shape::type_t::heightfield;
    value.heightfield = tvalue;
  } else {
    throw std::runtime_error("unknown Shape " + type);
  }
}

// Parse AreaLightSource
static inline void parse_pbrt_arealight(
    string_view& str, const string& type, pbrt_arealight& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "diffuse") {
    auto tvalue = pbrt_arealight::diffuse_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "twosided") {
        parse_pbrt_param(str, ptype, tvalue.twosided);
      } else if (pname == "samples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "nsamples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_arealight::type_t::diffuse;
    value.diffuse = tvalue;
  } else {
    throw std::runtime_error("unknown Film " + type);
  }
}

// Parse LightSource
static inline void parse_pbrt_light(
    string_view& str, const string& type, pbrt_light& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "distant") {
    auto tvalue = pbrt_light::distant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else if (pname == "to") {
        parse_pbrt_param(str, ptype, tvalue.to);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_light::type_t::distant;
    value.distant = tvalue;
  } else if (type == "goniometric") {
    auto tvalue = pbrt_light::goniometric_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_light::type_t::goniometric;
    value.goniometric = tvalue;
  } else if (type == "infinite") {
    auto tvalue = pbrt_light::infinite_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "samples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "nsamples") {
        parse_pbrt_param(str, ptype, tvalue.samples);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type     = pbrt_light::type_t::infinite;
    value.infinite = tvalue;
  } else if (type == "distant") {
    auto tvalue = pbrt_light::distant_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "L") {
        parse_pbrt_param(str, ptype, tvalue.L);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type    = pbrt_light::type_t::distant;
    value.distant = tvalue;
  } else if (type == "projection") {
    auto tvalue = pbrt_light::projection_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "fov") {
        parse_pbrt_param(str, ptype, tvalue.fov);
      } else if (pname == "mapname") {
        parse_pbrt_param(str, ptype, tvalue.mapname);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type       = pbrt_light::type_t::projection;
    value.projection = tvalue;
  } else if (type == "spot") {
    auto tvalue = pbrt_light::spot_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else if (pname == "to") {
        parse_pbrt_param(str, ptype, tvalue.to);
      } else if (pname == "coneangle") {
        parse_pbrt_param(str, ptype, tvalue.coneangle);
      } else if (pname == "conedeltaangle") {
        parse_pbrt_param(str, ptype, tvalue.conedeltaangle);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type = pbrt_light::type_t::spot;
    value.spot = tvalue;
  } else if (type == "point") {
    auto tvalue = pbrt_light::point_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "I") {
        parse_pbrt_param(str, ptype, tvalue.I);
      } else if (pname == "from") {
        parse_pbrt_param(str, ptype, tvalue.from);
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type  = pbrt_light::type_t::point;
    value.point = tvalue;
  } else {
    throw std::runtime_error("unknown LightSource " + type);
  }
}

// Parse Medium
static inline void parse_pbrt_medium(
    string_view& str, const string& type, pbrt_medium& value) {
  auto pname = ""s, ptype = ""s;
  if (type == "homogeneous") {
    auto tvalue = pbrt_medium::homogeneous_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "sigma_a") {
        parse_pbrt_param(str, ptype, tvalue.sigma_a);
      } else if (pname == "sigma_s") {
        parse_pbrt_param(str, ptype, tvalue.sigma_s);
      } else if (pname == "preset") {
        parse_pbrt_param(str, ptype, tvalue.preset);
      } else if (pname == "g") {
        parse_pbrt_param(str, ptype, tvalue.g);
      } else if (pname == "scale") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type        = pbrt_medium::type_t::homogeneous;
    value.homogeneous = tvalue;
  } else if (type == "heterogeneous") {
    auto tvalue = pbrt_medium::heterogeneous_t{};
    while (is_pbrt_param(str)) {
      parse_pbrt_nametype(str, pname, ptype);
      if (pname == "sigma_a") {
        parse_pbrt_param(str, ptype, tvalue.scale);
      } else if (pname == "sigma_s") {
        parse_pbrt_param(str, ptype, tvalue.sigma_s);
      } else if (pname == "preset") {
        parse_pbrt_param(str, ptype, tvalue.preset);
      } else if (pname == "g") {
        parse_pbrt_param(str, ptype, tvalue.g);
      } else if (pname == "p0") {
        parse_pbrt_param(str, ptype, tvalue.p0);
      } else if (pname == "p1") {
        parse_pbrt_param(str, ptype, tvalue.p1);
      } else if (pname == "nx") {
        parse_pbrt_param(str, ptype, tvalue.nx);
      } else if (pname == "ny") {
        parse_pbrt_param(str, ptype, tvalue.ny);
      } else if (pname == "nz") {
        parse_pbrt_param(str, ptype, tvalue.nz);
      } else if (pname == "density") {
        parse_pbrt_param(str, ptype, tvalue.density);
      } else if (pname == "type") {
        auto ttype = ""s;
        parse_pbrt_param(str, ptype, ttype);
        if (ttype != type) throw std::runtime_error("inconsistent types");
      } else {
        throw std::runtime_error("unknown parameter " + pname);
      }
    }
    value.type          = pbrt_medium::type_t::heterogeneous;
    value.heterogeneous = tvalue;
  } else {
    throw std::runtime_error("unknown Medium " + type);
  }
}

// Load pbrt scene
void load_pbrt(const string& filename, pbrt_callbacks& cb, bool flipv) {
  // start laoding files
  auto files = vector<file_wrapper>{};
  open_file(files.emplace_back(), filename);

  // parsing stack
  auto stack    = vector<pbrt_context>{{}};
  auto object   = pbrt_object{};
  auto coordsys = unordered_map<string, pair<frame3f, frame3f>>{};

  // helpders
  auto set_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end = (frame3f)xform;
  };
  auto concat_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end *= (frame3f)xform;
  };

  // constant values
  unordered_map<string, pbrt_spectrum3f> constant_values = {};

  // parse command by command
  while (!files.empty()) {
    auto& fs   = files.back();
    auto line = ""s;
    auto cmd  = ""s;
    while (read_pbrt_cmdline(fs, line)) {
      auto str = string_view{line};
      // get command
      parse_pbrt_command(str, cmd);
      if (cmd == "WorldBegin") {
        stack.push_back({});
      } else if (cmd == "WorldEnd") {
        stack.pop_back();
        if (stack.size() != 1) throw std::runtime_error("bad stack");
      } else if (cmd == "AttributeBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "AttributeEnd") {
        stack.pop_back();
      } else if (cmd == "TransformBegin") {
        stack.push_back(stack.back());
      } else if (cmd == "TransformEnd") {
        stack.pop_back();
      } else if (cmd == "ObjectBegin") {
        parse_pbrt_value(str, object.name);
        stack.push_back(stack.back());
        cb.begin_object(object, stack.back());
      } else if (cmd == "ObjectEnd") {
        cb.end_object(object, stack.back());
        stack.pop_back();
        object = {};
      } else if (cmd == "ObjectInstance") {
        auto value = pbrt_object{};
        parse_pbrt_value(str, value.name);
        cb.object_instance(value, stack.back());
      } else if (cmd == "ActiveTransform") {
        auto value = ""s;
        parse_pbrt_command(str, value);
        if (value == "StartTime") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = false;
        } else if (value == "EndTime") {
          stack.back().active_transform_start = false;
          stack.back().active_transform_end   = true;
        } else if (value == "All") {
          stack.back().active_transform_start = true;
          stack.back().active_transform_end   = true;
        } else {
          throw std::runtime_error("bad active transform");
        }
      } else if (cmd == "Transform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        set_transform(stack.back(), xf);
      } else if (cmd == "ConcatTransform") {
        auto xf = identity4x4f;
        parse_pbrt_param(str, xf);
        concat_transform(stack.back(), xf);
      } else if (cmd == "Scale") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), (mat4f)scaling_frame(v));
      } else if (cmd == "Translate") {
        auto v = zero3f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(), (mat4f)translation_frame(v));
      } else if (cmd == "Rotate") {
        auto v = zero4f;
        parse_pbrt_param(str, v);
        concat_transform(stack.back(),
            (mat4f)rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
      } else if (cmd == "LookAt") {
        auto from = zero3f, to = zero3f, up = zero3f;
        parse_pbrt_param(str, from);
        parse_pbrt_param(str, to);
        parse_pbrt_param(str, up);
        // from pbrt parser
        auto frame = lookat_frame(from, to, up, true);
        // frame.z = normalize(to-from);
        // frame.x = normalize(cross(frame.z,up));
        // frame.y = cross(frame.x,frame.z);
        // frame.o    = from;
        concat_transform(stack.back(), (mat4f)inverse(frame));
        stack.back().last_lookat_distance = length(from - to);
        // stack.back().focus = length(m.x - m.y);
      } else if (cmd == "ReverseOrientation") {
        stack.back().reverse = !stack.back().reverse;
      } else if (cmd == "CoordinateSystem") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        coordsys[name] = {
            stack.back().transform_start, stack.back().transform_end};
      } else if (cmd == "CoordSysTransform") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        if (coordsys.find(name) != coordsys.end()) {
          stack.back().transform_start = coordsys.at(name).first;
          stack.back().transform_end   = coordsys.at(name).second;
        }
      } else if (cmd == "Integrator") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_integrator{};
        parse_pbrt_integrator(str, type, value);
        cb.integrator(value, stack.back());
      } else if (cmd == "Sampler") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_sampler{};
        parse_pbrt_sampler(str, type, value);
        cb.sampler(value, stack.back());
      } else if (cmd == "PixelFilter") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_filter{};
        parse_pbrt_filter(str, type, value);
        cb.filter(value, stack.back());
      } else if (cmd == "Film") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_film{};
        parse_pbrt_film(str, type, value);
        cb.film(value, stack.back());
      } else if (cmd == "Accelerator") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_accelerator{};
        parse_pbrt_accelerator(str, type, value);
        cb.accelerator(value, stack.back());
      } else if (cmd == "Camera") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_camera{};
        parse_pbrt_camera(str, type, value);
        cb.camera(value, stack.back());
      } else if (cmd == "Texture") {
        auto name = ""s, comptype = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_value(str, comptype);
        parse_pbrt_value(str, type);
        auto value = pbrt_texture{};
        parse_pbrt_texture(str, type, value);
        if (type == "constant") {
          constant_values[name] = value.constant.value.value;
        }
        cb.texture(value, name, stack.back());
      } else if (cmd == "Material") {
        static auto material_id = 0;
        auto        type        = ""s;
        parse_pbrt_value(str, type);
        if (type == "") {
          stack.back().material = "";
        } else {
          auto value = pbrt_material{};
          auto name  = "unnamed_material_" + std::to_string(material_id++);
          parse_pbrt_material(str, type, value, constant_values);
          stack.back().material = name;
          cb.material(value, name, stack.back());
        }
      } else if (cmd == "MakeNamedMaterial") {
        auto name = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_typeparam(str, type);
        auto value = pbrt_material{};
        parse_pbrt_material(str, type, value, constant_values);
        cb.material(value, name, stack.back());
      } else if (cmd == "NamedMaterial") {
        auto name = ""s;
        parse_pbrt_value(str, name);
        stack.back().material = name;
      } else if (cmd == "Shape") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_shape{};
        parse_pbrt_shape(str, type, value);
        cb.shape(value, stack.back());
      } else if (cmd == "AreaLightSource") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        static auto material_id = 0;
        auto        name = "unnamed_arealight_" + std::to_string(material_id++);
        auto        value = pbrt_arealight{};
        parse_pbrt_arealight(str, type, value);
        stack.back().arealight = name;
        cb.arealight(value, name, stack.back());
      } else if (cmd == "LightSource") {
        auto type = ""s;
        parse_pbrt_value(str, type);
        auto value = pbrt_light{};
        parse_pbrt_light(str, type, value);
        cb.light(value, stack.back());
      } else if (cmd == "MakeNamedMedium") {
        auto name = ""s, type = ""s;
        parse_pbrt_value(str, name);
        parse_pbrt_typeparam(str, type);
        auto value = pbrt_medium{};
        parse_pbrt_medium(str, type, value);
        cb.medium(value, name, stack.back());
      } else if (cmd == "MediumInterface") {
        auto interior = ""s, exterior = ""s;
        parse_pbrt_value(str, interior);
        parse_pbrt_value(str, exterior);
        stack.back().medium_interior = interior;
        stack.back().medium_exterior = exterior;
      } else if (cmd == "Include") {
        auto inputname = ""s;
        parse_pbrt_value(str, inputname);
        open_file(
            files.emplace_back(), fs::path(filename).parent_path() / inputname);
      } else {
        throw std::runtime_error("unknown command " + cmd);
      }
    }
    files.pop_back();
  }
}

// Load pbrt scene
bool read_pbrt_element(file_wrapper& fs, pbrt_element& element, string& name,
    pbrt_element_data& data, vector<pbrt_context>& stack,
    pbrt_parser_state& state) {
  // helpders
  auto set_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start = (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end = (frame3f)xform;
  };
  auto concat_transform = [](pbrt_context& ctx, const mat4f& xform) {
    if (ctx.active_transform_start) ctx.transform_start *= (frame3f)xform;
    if (ctx.active_transform_end) ctx.transform_end *= (frame3f)xform;
  };

  // init stack
  if (stack.empty()) stack.emplace_back();

  // parse command by command
  while (read_pbrt_cmdline(fs, state.line)) {
    auto str = string_view{state.line};
    // get command
    auto cmd = ""s;
    parse_pbrt_command(str, cmd);
    if (cmd == "WorldBegin") {
      stack.push_back({});
    } else if (cmd == "WorldEnd") {
      stack.pop_back();
      if (stack.size() != 1) throw std::runtime_error("bad stack");
    } else if (cmd == "AttributeBegin") {
      stack.push_back(stack.back());
    } else if (cmd == "AttributeEnd") {
      stack.pop_back();
    } else if (cmd == "TransformBegin") {
      stack.push_back(stack.back());
    } else if (cmd == "TransformEnd") {
      stack.pop_back();
    } else if (cmd == "ObjectBegin") {
      parse_pbrt_value(str, state.object);
      stack.push_back(stack.back());
      element = pbrt_element::begin_object;
      name    = state.object;
      return true;
    } else if (cmd == "ObjectEnd") {
      element = pbrt_element::end_object;
      name    = state.object;
      stack.pop_back();
      state.object = {};
      return true;
    } else if (cmd == "ObjectInstance") {
      parse_pbrt_value(str, name);
      element = pbrt_element::object_instance;
      return true;
    } else if (cmd == "ActiveTransform") {
      auto value = ""s;
      parse_pbrt_command(str, value);
      if (value == "StartTime") {
        stack.back().active_transform_start = true;
        stack.back().active_transform_end   = false;
      } else if (value == "EndTime") {
        stack.back().active_transform_start = false;
        stack.back().active_transform_end   = true;
      } else if (value == "All") {
        stack.back().active_transform_start = true;
        stack.back().active_transform_end   = true;
      } else {
        throw std::runtime_error("bad active transform");
      }
    } else if (cmd == "Transform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      set_transform(stack.back(), xf);
    } else if (cmd == "ConcatTransform") {
      auto xf = identity4x4f;
      parse_pbrt_param(str, xf);
      concat_transform(stack.back(), xf);
    } else if (cmd == "Scale") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(), (mat4f)scaling_frame(v));
    } else if (cmd == "Translate") {
      auto v = zero3f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(), (mat4f)translation_frame(v));
    } else if (cmd == "Rotate") {
      auto v = zero4f;
      parse_pbrt_param(str, v);
      concat_transform(stack.back(),
          (mat4f)rotation_frame(vec3f{v.y, v.z, v.w}, radians(v.x)));
    } else if (cmd == "LookAt") {
      auto from = zero3f, to = zero3f, up = zero3f;
      parse_pbrt_param(str, from);
      parse_pbrt_param(str, to);
      parse_pbrt_param(str, up);
      // from pbrt parser
      auto frame = lookat_frame(from, to, up, true);
      // frame.z = normalize(to-from);
      // frame.x = normalize(cross(frame.z,up));
      // frame.y = cross(frame.x,frame.z);
      // frame.o    = from;
      concat_transform(stack.back(), (mat4f)inverse(frame));
      stack.back().last_lookat_distance = length(from - to);
      // stack.back().focus = length(m.x - m.y);
    } else if (cmd == "ReverseOrientation") {
      stack.back().reverse = !stack.back().reverse;
    } else if (cmd == "CoordinateSystem") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      state.coordsys[name] = {
          stack.back().transform_start, stack.back().transform_end};
    } else if (cmd == "CoordSysTransform") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      if (state.coordsys.find(name) != state.coordsys.end()) {
        stack.back().transform_start = state.coordsys.at(name).first;
        stack.back().transform_end   = state.coordsys.at(name).second;
      }
    } else if (cmd == "Integrator") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_integrator(str, type, data.intergrator);
      element = pbrt_element::integrator;
      return true;
    } else if (cmd == "Sampler") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_sampler(str, type, data.sampler);
      element = pbrt_element::sampler;
      return true;
    } else if (cmd == "PixelFilter") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_filter(str, type, data.filter);
      element = pbrt_element::filter;
      return true;
    } else if (cmd == "Film") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_film(str, type, data.film);
      element = pbrt_element::film;
      return true;
    } else if (cmd == "Accelerator") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_accelerator(str, type, data.accelerator);
      element = pbrt_element::accelerator;
      return true;
    } else if (cmd == "Camera") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_camera(str, type, data.camera);
      element = pbrt_element::camera;
      return true;
    } else if (cmd == "Texture") {
      auto comptype = ""s, type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_value(str, comptype);
      parse_pbrt_value(str, type);
      parse_pbrt_texture(str, type, data.texture);
      if (type == "constant") {
        state.constant_values[name] = data.texture.constant.value.value;
      }
      element = pbrt_element::texture;
      return true;
    } else if (cmd == "Material") {
      static auto material_id = 0;
      auto        type        = ""s;
      parse_pbrt_value(str, type);
      if (type == "") {
        stack.back().material = "";
      } else {
        name = "unnamed_material_" + std::to_string(material_id++);
        parse_pbrt_material(str, type, data.material, state.constant_values);
        stack.back().material = name;
        element               = pbrt_element::material;
        return true;
      }
    } else if (cmd == "MakeNamedMaterial") {
      auto type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_typeparam(str, type);
      parse_pbrt_material(str, type, data.material, state.constant_values);
      element = pbrt_element::material;
      return true;
    } else if (cmd == "NamedMaterial") {
      auto name = ""s;
      parse_pbrt_value(str, name);
      stack.back().material = name;
    } else if (cmd == "Shape") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_shape(str, type, data.shape);
      element = pbrt_element::shape;
      return true;
    } else if (cmd == "AreaLightSource") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      static auto material_id = 0;
      name = "unnamed_arealight_" + std::to_string(material_id++);
      parse_pbrt_arealight(str, type, data.arealight);
      stack.back().arealight = name;
      element                = pbrt_element::arealight;
      return true;
    } else if (cmd == "LightSource") {
      auto type = ""s;
      parse_pbrt_value(str, type);
      parse_pbrt_light(str, type, data.light);
      element = pbrt_element::light;
      return true;
    } else if (cmd == "MakeNamedMedium") {
      auto type = ""s;
      parse_pbrt_value(str, name);
      parse_pbrt_typeparam(str, type);
      parse_pbrt_medium(str, type, data.medium);
      element = pbrt_element::medium;
      return true;
    } else if (cmd == "MediumInterface") {
      auto interior = ""s, exterior = ""s;
      parse_pbrt_value(str, interior);
      parse_pbrt_value(str, exterior);
      stack.back().medium_interior = interior;
      stack.back().medium_exterior = exterior;
    } else if (cmd == "Include") {
      parse_pbrt_value(str, name);
      element = pbrt_element::include;
      return true;
    } else {
      throw std::runtime_error("unknown command " + cmd);
    }
  }
  return false;
}

}  // namespace yocto
