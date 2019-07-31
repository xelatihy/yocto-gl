//
// Implementation for Yocto/Ply.
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

#include "yocto_ply.h"

#include <algorithm>
#include <cinttypes>
#include <climits>
#include <limits>
#include <string_view>

#include "ext/filesystem.hpp"
namespace fs = ghc::filesystem;

// -----------------------------------------------------------------------------
// PLY CONVERSION
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
static inline bool read_ply_line(FILE* fs, char* buffer, size_t size) {
  return fgets(buffer, size, fs) != nullptr;
}

static inline bool is_ply_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
static inline bool is_ply_newline(char c) { return c == '\r' || c == '\n'; }

static inline void skip_ply_whitespace(string_view& str) {
  while (!str.empty() && is_ply_space(str.front())) str.remove_prefix(1);
}
static inline void remove_ply_comment(
    string_view& str, char comment_char = '#') {
  while (!str.empty() && is_ply_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
static inline void parse_ply_value(string_view& str, string_view& value) {
  skip_ply_whitespace(str);
  if (str.empty()) throw std::runtime_error("cannot parse value");
  if (str.front() != '"') {
    auto cpy = str;
    while (!cpy.empty() && !is_ply_space(cpy.front())) cpy.remove_prefix(1);
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
static inline void parse_ply_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_ply_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_ply_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_ply_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_ply_value(string_view& str, int8_t& value) {
  auto vali = (int32_t)0;
  parse_ply_value(str, vali);
  if (vali < INT8_MIN || vali > INT8_MAX)
    throw std::runtime_error("cannot parse value");
  value = (int8_t)vali;
}
static inline void parse_ply_value(string_view& str, int16_t& value) {
  auto vali = (int32_t)0;
  parse_ply_value(str, vali);
  if (vali < INT16_MIN || vali > INT16_MAX)
    throw std::runtime_error("cannot parse value");
  value = (int16_t)vali;
}
static inline void parse_ply_value(string_view& str, uint8_t& value) {
  auto vali = (uint32_t)0;
  parse_ply_value(str, vali);
  if (vali > UINT8_MAX) throw std::runtime_error("cannot parse value");
  value = (uint8_t)vali;
}
static inline void parse_ply_value(string_view& str, uint16_t& value) {
  auto vali = (uint32_t)0;
  parse_ply_value(str, vali);
  if (vali > UINT16_MAX) throw std::runtime_error("cannot parse value");
  value = (uint16_t)vali;
}
#ifdef __APPLE__
static inline void parse_ply_value(string_view& str, size_t& value) {
  parse_ply_value(str, (uint64_t&)value);
}
#endif
static inline void parse_ply_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
static inline void parse_ply_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}

// get ply value either ascii or binary
template <typename T, typename VT>
static inline void read_ply_prop(FILE* fs, bool big_endian, VT& value) {
  auto typed_value = (T)0;
  if (fread(&typed_value, sizeof(T), 1, fs) != 1)
    throw std::runtime_error("cannot read value");
  if (big_endian) typed_value = swap_endian(typed_value);
  value = (VT)typed_value;
}
template <typename VT>
static inline void read_ply_prop(
    FILE* fs, bool big_endian, ply_type type, VT& value) {
  switch (type) {
    case ply_type::i8: read_ply_prop<int8_t>(fs, big_endian, value); break;
    case ply_type::i16: read_ply_prop<int16_t>(fs, big_endian, value); break;
    case ply_type::i32: read_ply_prop<int32_t>(fs, big_endian, value); break;
    case ply_type::i64: read_ply_prop<int64_t>(fs, big_endian, value); break;
    case ply_type::u8: read_ply_prop<uint8_t>(fs, big_endian, value); break;
    case ply_type::u16: read_ply_prop<uint16_t>(fs, big_endian, value); break;
    case ply_type::u32: read_ply_prop<uint32_t>(fs, big_endian, value); break;
    case ply_type::u64: read_ply_prop<uint64_t>(fs, big_endian, value); break;
    case ply_type::f32: read_ply_prop<float>(fs, big_endian, value); break;
    case ply_type::f64: read_ply_prop<double>(fs, big_endian, value); break;
  }
}

template <typename T, typename VT>
static inline void parse_ply_prop(string_view& str, VT& value) {
  auto typed_value = (T)0;
  parse_ply_value(str, typed_value);
  value = (VT)typed_value;
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
void read_ply_header(
    ply_file& ply, vector<ply_element>& elements, vector<string>& comments) {
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
  };

  // parsing checks
  auto first_line = true;
  auto end_header = false;

  // prepare elements
  elements.clear();

  // read the file header line by line
  char buffer[4096];
  while (read_ply_line(ply.fs, buffer, sizeof(buffer))) {
    // line
    auto line = string_view{buffer};
    remove_ply_comment(line);
    skip_ply_whitespace(line);
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
        ply.ascii = true;
      } else if (fmt == "binary_little_endian") {
        ply.ascii      = false;
        ply.big_endian = false;
      } else if (fmt == "binary_big_endian") {
        ply.ascii      = false;
        ply.big_endian = true;
      } else {
        throw std::runtime_error{"unknown ply format"};
      }
    } else if (cmd == "comment") {
      skip_ply_whitespace(line);
      comments.push_back(string{line});
    } else if (cmd == "obj_info") {
      skip_ply_whitespace(line);
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
void read_ply_value_impl(ply_file& ply, const ply_element& element,
    vector<VT>& values, vector<vector<LT>>& lists) {
  // prepare properties
  if (values.size() != element.properties.size()) {
    values.resize(element.properties.size());
  }
  if (lists.size() != element.properties.size()) {
    lists.resize(element.properties.size());
  }
  for (auto& list : lists) list.clear();

  // read property values
  if (ply.ascii) {
    char buffer[4096];
    if (!read_ply_line(ply.fs, buffer, sizeof(buffer)))
      throw std::runtime_error("cannot read ply " + ply.filename);
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
      read_ply_prop(ply.fs, ply.big_endian, prop.value_type, value);
      if (prop.is_list) {
        list.resize((int)value);
        for (auto i = 0; i < (int)value; i++)
          read_ply_prop(ply.fs, ply.big_endian, prop.list_type, list[i]);
      }
    }
  }
}

// Write text to file
static inline void write_ply_value(FILE* fs, int8_t value) {
  if (fprintf(fs, "%d", (int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, int16_t value) {
  if (fprintf(fs, "%d", (int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, int32_t value) {
  if (fprintf(fs, "%d", (int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, int64_t value) {
  if (fprintf(fs, "%lld", (long long)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, uint8_t value) {
  if (fprintf(fs, "%u", (unsigned int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, uint16_t value) {
  if (fprintf(fs, "%u", (unsigned int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, uint32_t value) {
  if (fprintf(fs, "%u", (unsigned int)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, uint64_t value) {
  if (fprintf(fs, "%llu", (unsigned long long)value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, float value) {
  if (fprintf(fs, "%g", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, double value) {
  if (fprintf(fs, "%g", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, const char* value) {
  if (fprintf(fs, "%s", value) < 0)
    throw std::runtime_error("cannot print value");
}
static inline void write_ply_value(FILE* fs, const string& value) {
  if (fprintf(fs, "%s", value.c_str()) < 0)
    throw std::runtime_error("cannot print value");
}

template <typename T, typename VT>
static inline void write_ply_prop(FILE* fs, VT value) {
  write_ply_value(fs, (T)value);
}

template <typename VT>
static inline void write_ply_prop(FILE* fs, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: write_ply_prop<int8_t>(fs, value); break;
    case ply_type::i16: write_ply_prop<int16_t>(fs, value); break;
    case ply_type::i32: write_ply_prop<int32_t>(fs, value); break;
    case ply_type::i64: write_ply_prop<int64_t>(fs, value); break;
    case ply_type::u8: write_ply_prop<uint8_t>(fs, value); break;
    case ply_type::u16: write_ply_prop<uint16_t>(fs, value); break;
    case ply_type::u32: write_ply_prop<uint32_t>(fs, value); break;
    case ply_type::u64: write_ply_prop<uint64_t>(fs, value); break;
    case ply_type::f32: write_ply_prop<float>(fs, value); break;
    case ply_type::f64: write_ply_prop<double>(fs, value); break;
  }
}

template <typename T, typename VT>
static inline void write_ply_binprop(FILE* fs, bool big_endian, VT value) {
  auto typed_value = (T)value;
  if (big_endian) typed_value = swap_endian(typed_value);
  if (fwrite(&typed_value, sizeof(T), 1, fs) != 1)
    throw std::runtime_error("cannot write to file");
}

template <typename VT>
static inline void write_ply_binprop(
    FILE* fs, bool big_endian, ply_type type, VT value) {
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
void write_ply_header(ply_file& ply, const vector<ply_element>& elements,
    const vector<string>& comments) {
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

  write_ply_value(ply.fs, "ply\n");
  if (ply.ascii) {
    write_ply_value(ply.fs, "format ascii 1.0\n");
  } else if (ply.big_endian) {
    write_ply_value(ply.fs, "format binary_big_endian 1.0\n");
  } else {
    write_ply_value(ply.fs, "format binary_little_endian 1.0\n");
  }
  for (auto& comment : comments)
    write_ply_value(ply.fs, "comment " + comment + "\n");
  for (auto& elem : elements) {
    write_ply_value(ply.fs,
        "element " + elem.name + " " + std::to_string(elem.count) + "\n");
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        write_ply_value(ply.fs, "property list " + type_map[prop.value_type] +
                                    " " + type_map[prop.list_type] + " " +
                                    prop.name + "\n");
      } else {
        write_ply_value(ply.fs,
            "property " + type_map[prop.value_type] + " " + prop.name + "\n");
      }
    }
  }
  write_ply_value(ply.fs, "end_header\n");
}

template <typename VT, typename LT>
void write_ply_value_impl(ply_file& ply, const ply_element& element,
    vector<VT>& values, vector<vector<LT>>& lists) {
  if (ply.ascii) {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (pidx) write_ply_value(ply.fs, " ");
      write_ply_prop(ply.fs, prop.value_type, values[pidx]);
      if (prop.is_list) {
        for (auto i = 0; i < (int)lists[pidx].size(); i++) {
          if (i) write_ply_value(ply.fs, " ");
          write_ply_prop(ply.fs, prop.list_type, lists[pidx][i]);
        }
      }
      write_ply_value(ply.fs, "\n");
    }
  } else {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      write_ply_binprop(ply.fs, ply.big_endian, prop.value_type, values[pidx]);
      if (prop.is_list) {
        for (auto i = 0; i < (int)lists[pidx].size(); i++)
          write_ply_binprop(
              ply.fs, ply.big_endian, prop.list_type, lists[pidx][i]);
      }
    }
  }
}

void write_ply_value(ply_file& ply, const ply_element& element,
    vector<double>& values, vector<vector<double>>& lists) {
  write_ply_value_impl(ply, element, values, lists);
}
void write_ply_value(ply_file& ply, const ply_element& element,
    vector<float>& values, vector<vector<int>>& lists) {
  write_ply_value_impl(ply, element, values, lists);
}

void read_ply_value(ply_file& ply, const ply_element& element,
    vector<double>& values, vector<vector<double>>& lists) {
  read_ply_value_impl(ply, element, values, lists);
}
void read_ply_value(ply_file& ply, const ply_element& element,
    vector<float>& values, vector<vector<int>>& lists) {
  read_ply_value_impl(ply, element, values, lists);
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
