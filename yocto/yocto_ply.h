//
// # Yocto/Ply: Tiny library for Ply parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Stanford Ply. Yocto/Ply
// supports two interfaces: a simple interface where all Ply data is loaded
// and saved at once and a low-level interface where Ply values are read
// and written one at a time.
// Error reporting is done by throwing `std::runtime_error` exceptions.
//
//
// ## Low-Level Ply Loading
//
// Load a PLY by first opening the file and reading its header. Then, for each
// element, read the values of its lists and non-lists properties. Example:
//
//    auto ply = fopen(filename, "rb");                // open for reading
//    auto format = ply_format{};                      // initialize format
//    auto elemnts = vector<ply_element>{};            // initialize elements
//    auto comments = vector<string>{};                // initialize comments
//    read_ply_header(ply, fromat elements, comments); // read ply header
//    for(auto& element : elements) {                  // iterate elements
//      // initialize the element's property values and lists
//      // using either doubles or vector<float> and vector<vector<int>>
//      auto values = vector<double>(element.properties.size());
//      auto lists - vector<vector<double>>(element.properties.size());
//      for(auto i = 0; i < element.count; i ++) {             // iterate values
//        read_ply_value(ply, format, element, values, lists); // read props
//        // values contains values for non-list properties
//        // lists contains the values for list properties
//    }
//
// For convenience during parsing, you can use `find_ply_property()` to
// determine the index of the property you may be interested in.
//
//
// ## Load-Level PLY Saving
//
// Write a PLY by first opening the file for writing and deciding whether to
// use ASCII or binary (we recommend tha letter). Then fill in the elements
// and comments and write its header. Finally, write its values one by one.
// Example:
//
//    auto fs = fopen(filename, "rb");                   // open for writing
//    auto format = ply_format::binary_little_endian;    // initialize format
//    auto elemnts = vector<ply_element>{};              // initialize elements
//    auto comments = vector<string>{};                  // initialize comments
//    // add eleements and comments to the previous lists
//    write_ply_header(ply, format, elements, comments); // read ply header
//    for(auto& element : elements) {                    // iterate elements
//      // initialize the element's property values and lists
//      // using either doubles or vector<float> and vector<vector<int>>
//      auto values = vector<double>(element.properties.size());
//      auto lists - vector<vector<double>>(element.properties.size());
//      for(auto i = 0; i < element.count; i ++) {       // iterate values
//        values = {...}; lists = {...};                 // set values/lists
//        write_ply_value(ply, foramt, element, values, lists); // write props
//    }
//
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

#include "yocto_math.h"

#include <algorithm>

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

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
  string               name       = "";
  size_t               count      = 0;
  vector<ply_property> properties = {};
};

// Ply model
struct ply_model {
  ply_format          format   = ply_format::binary_little_endian;
  vector<string>      comments = {};
  vector<ply_element> elements = {};
};

// Load and save ply
inline void load_ply(const string& filename, ply_model& ply);
inline void save_ply(const string& filename, const ply_model& ply);

// Get ply properties
inline bool has_ply_property(
    const ply_model& ply, const string& element, const string& property);
inline const ply_property& get_ply_property(
    const ply_model& ply, const string& element, const string& property);

inline vector<float> get_ply_values(
    const ply_model& ply, const string& element, const string& property);
inline vector<vec2f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2);
inline vector<vec3f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3);
inline vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    const string& property4);
inline vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    float property4);

inline vector<vector<int>> get_ply_lists(
    const ply_model& ply, const string& element, const string& property);
inline vector<byte> get_ply_list_sizes(
    const ply_model& ply, const string& element, const string& property);
inline vector<int> get_ply_list_values(
    const ply_model& ply, const string& element, const string& property);
inline vec2i get_ply_list_minxmax(
    const ply_model& ply, const string& element, const string& property);

// Get ply properties for meshes
inline vector<vec3f> get_ply_positions(const ply_model& ply);
inline vector<vec3f> get_ply_normals(const ply_model& ply);
inline vector<vec2f> get_ply_texcoords(
    const ply_model& ply, bool flipv = false);
inline vector<vec4f> get_ply_colors(const ply_model& ply);
inline vector<float> get_ply_radius(const ply_model& ply);
inline vector<vector<int>>  get_ply_faces(const ply_model& ply);
inline vector<vec2i> get_ply_lines(const ply_model& ply);
inline vector<int>   get_ply_points(const ply_model& ply);
inline vector<vec3i> get_ply_triangles(const ply_model& ply);
inline vector<vec4i> get_ply_quads(const ply_model& ply);
inline bool          has_ply_quads(const ply_model& ply);

// Add ply properties
inline void add_ply_values(ply_model& ply, const vector<float>& values,
    const string& element, const string& property);
inline void add_ply_values(ply_model& ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2);
inline void add_ply_values(ply_model& ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3);
inline void add_ply_values(ply_model& ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4);

inline void add_ply_lists(ply_model& ply, const vector<vector<int>>& values,
    const string& element, const string& property);
inline void add_ply_lists(ply_model& ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property);
inline void add_ply_lists(ply_model& ply, const vector<int>& values,
    const string& element, const string& property);
inline void add_ply_lists(ply_model& ply, const vector<vec2i>& values,
    const string& element, const string& property);
inline void add_ply_lists(ply_model& ply, const vector<vec3i>& values,
    const string& element, const string& property);
inline void add_ply_lists(ply_model& ply, const vector<vec4i>& values,
    const string& element, const string& property);

// Add ply properties for meshes
inline void add_ply_positions(ply_model& ply, const vector<vec3f>& values);
inline void add_ply_normals(ply_model& ply, const vector<vec3f>& values);
inline void add_ply_texcoords(
    ply_model& ply, const vector<vec2f>& values, bool flipv = false);
inline void add_ply_colors(ply_model& ply, const vector<vec4f>& values);
inline void add_ply_radius(ply_model& ply, const vector<float>& values);
inline void add_ply_faces(ply_model& ply, const vector<vector<int>>& values);
inline void add_ply_faces(
    ply_model& ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
inline void add_ply_triangles(ply_model& ply, const vector<vec3i>& values);
inline void add_ply_quads(ply_model& ply, const vector<vec4i>& values);
inline void add_ply_lines(ply_model& ply, const vector<vec2i>& values);
inline void add_ply_points(ply_model& ply, const vector<int>& values);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW_LEVEL PLY LOADING AND SAVING
// -----------------------------------------------------------------------------
namespace yocto {

// A class that wraps a C file ti handle safe opening/closgin with RIIA.
struct ply_file {
  ply_file() {}
  ply_file(ply_file&& other);
  ply_file(const ply_file&) = delete;
  ply_file& operator=(const ply_file&) = delete;
  ~ply_file();

  FILE*  fs       = nullptr;
  string filename = "";
  string mode     = "rt";
  int    linenum  = 0;
};

// open a file
inline ply_file open_ply(const string& filename, const string& mode = "rt");
inline void     open_ply(
        ply_file& fs, const string& filename, const string& mode = "rt");
inline void close_ply(ply_file& fs);

// Read Ply functions
inline void read_ply_header(ply_file& fs, ply_format& format,
    vector<ply_element>& elements, vector<string>& comments);
inline void read_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
inline void read_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Write Ply functions
inline void write_ply_header(ply_file& fs, ply_format format,
    const vector<ply_element>& elements, const vector<string>& comments);
inline void write_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
inline void write_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Helpers to get element and property indices
inline int find_ply_element(
    const vector<ply_element>& elements, const string& name);
inline int   find_ply_property(const ply_element& element, const string& name);
inline vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2);
inline vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3);
inline vec4i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3, const string& name4);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#include <string_view>

// -----------------------------------------------------------------------------
// LOW-LEVEL FILE HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// copnstrucyor and destructors
inline ply_file::ply_file(ply_file&& other) {
  this->fs       = other.fs;
  this->filename = other.filename;
  other.fs       = nullptr;
}
inline ply_file::~ply_file() {
  if (fs) fclose(fs);
  fs = nullptr;
}

// Opens a file returing a handle with RIIA
inline void open_ply(ply_file& fs, const string& filename, const string& mode) {
  close_ply(fs);
  fs.filename = filename;
  fs.mode     = mode;
  fs.fs       = fopen(filename.c_str(), mode.c_str());
  if (!fs.fs) throw std::runtime_error("could not open file " + filename);
}
inline ply_file open_ply(const string& filename, const string& mode) {
  auto fs = ply_file{};
  open_ply(fs, filename, mode);
  return fs;
}
inline void close_ply(ply_file& fs) {
  if (fs.fs) fclose(fs.fs);
  fs.fs = nullptr;
}

inline bool read_ply_line(ply_file& fs, char* buffer, size_t size) {
  auto ok = fgets(buffer, size, fs.fs) != nullptr;
  if (ok) fs.linenum += 1;
  return ok;
}

template <typename T>
inline T swap_ply_endian(T value) {
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
inline void write_ply_value(ply_file& fs, const T& value) {
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}
template <typename T>
inline void write_ply_value(ply_file& fs, const T& value_, bool big_endian) {
  auto value = big_endian ? swap_ply_endian(value_) : value_;
  if (fwrite(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot write to " + fs.filename);
}

inline void write_ply_text(ply_file& fs, const string& value) {
  if (fputs(value.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}
inline void write_ply_text(ply_file& fs, const char* value) {
  if (fputs(value, fs.fs) < 0)
    throw std::runtime_error("cannot write to " + fs.filename);
}

template <typename T>
inline void read_ply_value(ply_file& fs, T& value) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
}
template <typename T>
inline void read_ply_value(ply_file& fs, T& value, bool big_endian) {
  if (fread(&value, sizeof(value), 1, fs.fs) != 1)
    throw std::runtime_error("cannot read " + fs.filename);
  if (big_endian) value = swap_ply_endian(value);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOAD-LEVEL PARSING
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;

// utilities
inline bool is_ply_newline(char c) { return c == '\r' || c == '\n'; }
inline bool is_ply_space(char c) {
  return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}
inline void skip_ply_whitespace(string_view& str) {
  while (!str.empty() && is_ply_space(str.front())) str.remove_prefix(1);
}

inline void remove_ply_comment(string_view& str, char comment_char = '#') {
  while (!str.empty() && is_ply_newline(str.back())) str.remove_suffix(1);
  auto cpy = str;
  while (!cpy.empty() && cpy.front() != comment_char) cpy.remove_prefix(1);
  str.remove_suffix(cpy.size());
}

// Parse values from a string
inline void parse_ply_value(string_view& str, string_view& value) {
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
inline void parse_ply_value(string_view& str, string& value) {
  auto valuev = string_view{};
  parse_ply_value(str, valuev);
  value = string{valuev};
}
inline void parse_ply_value(string_view& str, int8_t& value) {
  char* end = nullptr;
  value     = (int8_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, int16_t& value) {
  char* end = nullptr;
  value     = (int16_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, int32_t& value) {
  char* end = nullptr;
  value     = (int32_t)strtol(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, int64_t& value) {
  char* end = nullptr;
  value     = (int64_t)strtoll(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, uint8_t& value) {
  char* end = nullptr;
  value     = (uint8_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, uint16_t& value) {
  char* end = nullptr;
  value     = (uint16_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, uint32_t& value) {
  char* end = nullptr;
  value     = (uint32_t)strtoul(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, uint64_t& value) {
  char* end = nullptr;
  value     = (uint64_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, bool& value) {
  auto valuei = 0;
  parse_ply_value(str, valuei);
  value = (bool)valuei;
}
inline void parse_ply_value(string_view& str, float& value) {
  char* end = nullptr;
  value     = strtof(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
inline void parse_ply_value(string_view& str, double& value) {
  char* end = nullptr;
  value     = strtod(str.data(), &end);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
template <typename T>
inline void parse_ply_value(string_view& str, T* values, int num) {
  for (auto i = 0; i < num; i++) parse_ply_value(str, values[i]);
}

inline void parse_ply_value(string_view& str, vec2f& value) {
  parse_ply_value(str, &value.x, 2);
}
inline void parse_ply_value(string_view& str, vec3f& value) {
  parse_ply_value(str, &value.x, 3);
}
inline void parse_ply_value(string_view& str, frame3f& value) {
  parse_ply_value(str, &value.x.x, 12);
}
#ifdef __APPLE__
inline void parse_ply_value(string_view& str, size_t& value) {
  char* end = nullptr;
  value     = (size_t)strtoull(str.data(), &end, 10);
  if (str == end) throw std::runtime_error("cannot parse value");
  str.remove_prefix(end - str.data());
}
#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL PRINTING
// -----------------------------------------------------------------------------
namespace yocto {

// Formats values to string
inline void format_ply_value(string& str, const string& value) { str += value; }
inline void format_ply_value(string& str, const char* value) { str += value; }
inline void format_ply_value(string& str, int8_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_ply_value(string& str, int16_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_ply_value(string& str, int32_t value) {
  char buf[256];
  sprintf(buf, "%d", (int)value);
  str += buf;
}
inline void format_ply_value(string& str, int64_t value) {
  char buf[256];
  sprintf(buf, "%lld", (long long)value);
  str += buf;
}
inline void format_ply_value(string& str, uint8_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_ply_value(string& str, uint16_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_ply_value(string& str, uint32_t value) {
  char buf[256];
  sprintf(buf, "%u", (unsigned)value);
  str += buf;
}
inline void format_ply_value(string& str, uint64_t value) {
  char buf[256];
  sprintf(buf, "%llu", (unsigned long long)value);
  str += buf;
}
inline void format_ply_value(string& str, float value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_ply_value(string& str, double value) {
  char buf[256];
  sprintf(buf, "%g", value);
  str += buf;
}
inline void format_ply_value(string& str, const vec2f& value) {
  char buf[256];
  sprintf(buf, "%g %g", value.x, value.y);
  str += buf;
}
inline void format_ply_value(string& str, const vec3f& value) {
  char buf[256];
  sprintf(buf, "%g %g %g", value.x, value.y, value.z);
  str += buf;
}
#if 0
inline void format_ply_value(string& str, const vec2i& value) {
  char buf[256];
  sprintf(buf, "%d %d", value.x, value.y);
  str += buf;
}
inline void format_ply_value(string& str, const vec3i& value) {
  char buf[256];
  sprintf(buf, "%d %d %d", value.x, value.y, value.z);
  str += buf;
}
#endif
inline void format_ply_value(string& str, const frame3f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g", value.x.x, value.x.y,
      value.x.z, value.y.x, value.y.y, value.y.z, value.z.x, value.z.y,
      value.z.z, value.o.x, value.o.y, value.o.z);
  str += buf;
}
inline void format_ply_value(string& str, const mat4f& value) {
  char buf[512];
  sprintf(buf, "%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g", value.x.x,
      value.x.y, value.x.z, value.x.w, value.y.x, value.y.y, value.y.z,
      value.y.w, value.z.x, value.z.y, value.z.z, value.z.w, value.w.x,
      value.w.y, value.w.z, value.w.w);
  str += buf;
}

// Foramt to file
inline void format_ply_values(string& str, const string& fmt) {
  auto pos = fmt.find("{}");
  if (pos != string::npos) throw std::runtime_error("bad format string");
  str += fmt;
}
template <typename Arg, typename... Args>
inline void format_ply_values(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
  auto pos = fmt.find("{}");
  if (pos == string::npos) throw std::runtime_error("bad format string");
  str += fmt.substr(0, pos);
  format_ply_value(str, arg);
  format_ply_values(str, fmt.substr(pos + 2), args...);
}

template <typename... Args>
inline void format_ply_values(
    ply_file& fs, const string& fmt, const Args&... args) {
  auto str = ""s;
  format_ply_values(str, fmt, args...);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}
template <typename T>
inline void format_ply_value(ply_file& fs, const T& value) {
  auto str = ""s;
  format_ply_value(str, value);
  if (fputs(str.c_str(), fs.fs) < 0)
    throw std::runtime_error("cannor write to " + fs.filename);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

// Load ply
inline void load_ply(const string& filename, ply_model& ply) {
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

  // read header ---------------------------------------------
  char buffer[4096];
  while (read_ply_line(fs, buffer, sizeof(buffer))) {
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
      auto fmt = string_view{};
      parse_ply_value(line, fmt);
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
      skip_ply_whitespace(line);
      ply.comments.push_back(string{line});
    } else if (cmd == "obj_info") {
      skip_ply_whitespace(line);
      // comment is the rest of the line
    } else if (cmd == "element") {
      auto& elem = ply.elements.emplace_back();
      parse_ply_value(line, elem.name);
      parse_ply_value(line, elem.count);
    } else if (cmd == "property") {
      if (ply.elements.empty()) throw std::runtime_error{"bad ply header"};
      auto& prop  = ply.elements.back().properties.emplace_back();
      auto  tname = ""s;
      parse_ply_value(line, tname);
      if (tname == "list") {
        prop.is_list = true;
        parse_ply_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          throw std::runtime_error{"unsupported list size type " + tname};
        parse_ply_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      }
      parse_ply_value(line, prop.name);
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
        if (!read_ply_line(fs, buffer, sizeof(buffer)))
          throw std::runtime_error("cannot read ply");
        auto line = string_view{buffer};
        for (auto& prop : elem.properties) {
          if (prop.is_list) parse_ply_value(line, prop.ldata_u8.emplace_back());
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                parse_ply_value(line, prop.data_i8.emplace_back());
                break;
              case ply_type::i16:
                parse_ply_value(line, prop.data_i16.emplace_back());
                break;
              case ply_type::i32:
                parse_ply_value(line, prop.data_i32.emplace_back());
                break;
              case ply_type::i64:
                parse_ply_value(line, prop.data_i64.emplace_back());
                break;
              case ply_type::u8:
                parse_ply_value(line, prop.data_u8.emplace_back());
                break;
              case ply_type::u16:
                parse_ply_value(line, prop.data_u16.emplace_back());
                break;
              case ply_type::u32:
                parse_ply_value(line, prop.data_u32.emplace_back());
                break;
              case ply_type::u64:
                parse_ply_value(line, prop.data_u64.emplace_back());
                break;
              case ply_type::f32:
                parse_ply_value(line, prop.data_f32.emplace_back());
                break;
              case ply_type::f64:
                parse_ply_value(line, prop.data_f64.emplace_back());
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
            read_ply_value(fs, prop.ldata_u8.emplace_back(), big_endian);
          auto vcount = prop.is_list ? prop.ldata_u8.back() : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                read_ply_value(fs, prop.data_i8.emplace_back(), big_endian);
                break;
              case ply_type::i16:
                read_ply_value(fs, prop.data_i16.emplace_back(), big_endian);
                break;
              case ply_type::i32:
                read_ply_value(fs, prop.data_i32.emplace_back(), big_endian);
                break;
              case ply_type::i64:
                read_ply_value(fs, prop.data_i64.emplace_back(), big_endian);
                break;
              case ply_type::u8:
                read_ply_value(fs, prop.data_u8.emplace_back(), big_endian);
                break;
              case ply_type::u16:
                read_ply_value(fs, prop.data_u16.emplace_back(), big_endian);
                break;
              case ply_type::u32:
                read_ply_value(fs, prop.data_u32.emplace_back(), big_endian);
                break;
              case ply_type::u64:
                read_ply_value(fs, prop.data_u64.emplace_back(), big_endian);
                break;
              case ply_type::f32:
                read_ply_value(fs, prop.data_f32.emplace_back(), big_endian);
                break;
              case ply_type::f64:
                read_ply_value(fs, prop.data_f64.emplace_back(), big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// Save ply
inline void save_ply(const string& filename, const ply_model& ply) {
  auto fs = open_ply(filename, "wb");

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
  format_ply_values(fs, "ply\n");
  format_ply_values(fs, "format {} 1.0\n", format_map.at(ply.format));
  format_ply_values(fs, "comment Written by Yocto/GL\n");
  format_ply_values(fs, "comment https://github.com/xelatihy/yocto-gl\n");
  for (auto& comment : ply.comments)
    format_ply_values(fs, "comment {}\n", comment);
  for (auto& elem : ply.elements) {
    format_ply_values(fs, "element {} {}\n", elem.name, (uint64_t)elem.count);
    for (auto& prop : elem.properties) {
      if (prop.is_list) {
        format_ply_values(
            fs, "property list uchar {} {}\n", type_map[prop.type], prop.name);
      } else {
        format_ply_values(
            fs, "property {} {}\n", type_map[prop.type], prop.name);
      }
    }
  }
  format_ply_values(fs, "end_header\n");

  // properties
  if (ply.format == ply_format::ascii) {
    for (auto& elem : ply.elements) {
      auto cur = vector<size_t>(elem.properties.size(), 0);
      for (auto idx = 0; idx < elem.count; idx++) {
        for (auto pidx = 0; pidx < elem.properties.size(); pidx++) {
          auto& prop = elem.properties[pidx];
          if (prop.is_list)
            format_ply_values(fs, "{} ", (int)prop.ldata_u8[idx]);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                format_ply_values(fs, "{} ", prop.data_i8[cur[idx]++]);
                break;
              case ply_type::i16:
                format_ply_values(fs, "{} ", prop.data_i16[cur[idx]++]);
                break;
              case ply_type::i32:
                format_ply_values(fs, "{} ", prop.data_i32[cur[idx]++]);
                break;
              case ply_type::i64:
                format_ply_values(fs, "{} ", prop.data_i64[cur[idx]++]);
                break;
              case ply_type::u8:
                format_ply_values(fs, "{} ", prop.data_i8[cur[idx]++]);
                break;
              case ply_type::u16:
                format_ply_values(fs, "{} ", prop.data_i16[cur[idx]++]);
                break;
              case ply_type::u32:
                format_ply_values(fs, "{} ", prop.data_u32[cur[idx]++]);
                break;
              case ply_type::u64:
                format_ply_values(fs, "{} ", prop.data_u64[cur[idx]++]);
                break;
              case ply_type::f32:
                format_ply_values(fs, "{}", prop.data_f32[cur[idx]++]);
                break;
              case ply_type::f64:
                format_ply_values(fs, "{}", prop.data_f64[cur[idx]++]);
                break;
            }
          }
          format_ply_values(fs, "\n");
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
          if (prop.is_list) write_ply_value(fs, prop.ldata_u8[idx], big_endian);
          auto vcount = prop.is_list ? prop.ldata_u8[idx] : 1;
          for (auto i = 0; i < vcount; i++) {
            switch (prop.type) {
              case ply_type::i8:
                write_ply_value(fs, prop.data_i8[cur[pidx]++], big_endian);
                break;
              case ply_type::i16:
                write_ply_value(fs, prop.data_i16[cur[pidx]++], big_endian);
                break;
              case ply_type::i32:
                write_ply_value(fs, prop.data_i32[cur[pidx]++], big_endian);
                break;
              case ply_type::i64:
                write_ply_value(fs, prop.data_i64[cur[pidx]++], big_endian);
                break;
              case ply_type::u8:
                write_ply_value(fs, prop.data_i8[cur[pidx]++], big_endian);
                break;
              case ply_type::u16:
                write_ply_value(fs, prop.data_i16[cur[pidx]++], big_endian);
                break;
              case ply_type::u32:
                write_ply_value(fs, prop.data_u32[cur[pidx]++], big_endian);
                break;
              case ply_type::u64:
                write_ply_value(fs, prop.data_u64[cur[pidx]++], big_endian);
                break;
              case ply_type::f32:
                write_ply_value(fs, prop.data_f32[cur[pidx]++], big_endian);
                break;
              case ply_type::f64:
                write_ply_value(fs, prop.data_f64[cur[pidx]++], big_endian);
                break;
            }
          }
        }
      }
    }
  }
}

// Get ply properties
inline bool has_ply_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return true;
    }
  }
  return false;
}
inline const ply_property& get_ply_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
inline ply_property& get_ply_property(
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
    default: throw std::runtime_error("should not be here");
  }
  // return here to silence warnings
  std::runtime_error("should not have gotten here");
  return {};
}
inline vector<float> get_ply_values(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (prop.is_list) return {};
  return convert_ply_property<float>(prop);
}
inline vector<vec2f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto values = vector<vec2f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++) values[i] = {x[i], y[i]};
  return values;
}
inline vector<vec3f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3) {
  auto x      = get_ply_values(ply, element, property1);
  auto y      = get_ply_values(ply, element, property2);
  auto z      = get_ply_values(ply, element, property3);
  auto values = vector<vec3f>(x.size());
  for (auto i = (size_t)0; i < values.size(); i++)
    values[i] = {x[i], y[i], z[i]};
  return values;
}
inline vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
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
inline vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
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
inline vector<vector<int>> get_ply_lists(
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
inline vector<byte> get_ply_list_sizes(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (!prop.is_list) return {};
  return prop.ldata_u8;
}
inline vector<int> get_ply_list_values(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_ply_property(ply, element, property)) return {};
  auto& prop = get_ply_property(ply, element, property);
  if (!prop.is_list) return {};
  return convert_ply_property<int>(prop);
}

inline vector<vec2f> flip_ply_texcoord(const vector<vec2f>& texcoord) {
  auto flipped = texcoord;
  for (auto& uv : flipped) uv.y = 1 - uv.y;
  return flipped;
}

// Get ply properties for meshes
inline vector<vec3f> get_ply_positions(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "x", "y", "z");
}
inline vector<vec3f> get_ply_normals(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "nx", "ny", "nz");
}
inline vector<vec2f> get_ply_texcoords(const ply_model& ply, bool flipv) {
  auto texcoord = has_ply_property(ply, "vertex", "u")
                      ? get_ply_values(ply, "vertex", "u", "v")
                      : get_ply_values(ply, "vertex", "s", "t");
  return flipv ? flip_ply_texcoord(texcoord) : texcoord;
}
inline vector<vec4f> get_ply_colors(const ply_model& ply) {
  if (has_ply_property(ply, "vertex", "alpha")) {
    return get_ply_values(ply, "vertex", "red", "green", "blue", "alpha");
  } else {
    return get_ply_values(ply, "vertex", "red", "green", "blue", 1);
  }
}
inline vector<float> get_ply_radius(const ply_model& ply) {
  return get_ply_values(ply, "vertex", "radius");
}
inline vector<vector<int>> get_ply_faces(const ply_model& ply) {
  return get_ply_lists(ply, "face", "vertex_indices");
}
inline vector<vec3i> get_ply_triangles(const ply_model& ply) {
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
inline vector<vec4i> get_ply_quads(const ply_model& ply) {
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
inline vector<vec2i> get_ply_lines(const ply_model& ply) {
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
inline vector<int> get_ply_points(const ply_model& ply) {
  return get_ply_list_values(ply, "point", "vertex_indices");
}
inline bool has_ply_quads(const ply_model& ply) {
  auto sizes = get_ply_list_sizes(ply, "face", "vertex_indices");
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}

// Add ply properties
inline void add_ply_element(
    ply_model& ply, const string& element, size_t count) {
  for (auto& elem : ply.elements) {
    if (elem.name == element) return;
  }
  auto& elem = ply.elements.emplace_back();
  elem.name  = element;
  elem.count = count;
}
inline void add_ply_property(ply_model& ply, const string& element,
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
inline vector<T> make_ply_vector(const T* value, size_t count, int stride) {
  auto ret = vector<T>(count);
  for (auto idx = (size_t)0; idx < count; idx++) ret[idx] = value[idx * stride];
  return ret;
}

inline void add_ply_values(ply_model& ply, const float* values, size_t count,
    const string& element, const string* properties, int nprops) {
  if (!values) return;
  for (auto p = 0; p < nprops; p++) {
    add_ply_property(ply, element, properties[p], count, ply_type::f32, false);
    auto& prop = get_ply_property(ply, element, properties[p]);
    prop.data_f32.resize(count);
    for (auto i = 0; i < count; i++) prop.data_f32[i] = values[p + i * nprops];
  }
}

inline void add_ply_values(ply_model& ply, const vector<float>& values,
    const string& element, const string& property) {
  auto properties = vector{property};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 1);
}
inline void add_ply_values(ply_model& ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2) {
  auto properties = vector{property1, property2};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 2);
}
inline void add_ply_values(ply_model& ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3) {
  auto properties = vector{property1, property2, property3};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 3);
}
inline void add_ply_values(ply_model& ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4) {
  auto properties = vector{property1, property2, property3, property4};
  add_ply_values(
      ply, (float*)values.data(), values.size(), element, properties.data(), 4);
}

inline void add_ply_lists(ply_model& ply, const vector<vector<int>>& values,
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
inline void add_ply_lists(ply_model& ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property) {
  if (values.empty()) return;
  add_ply_property(ply, element, property, sizes.size(), ply_type::i32, true);
  auto& prop    = get_ply_property(ply, element, property);
  prop.data_i32 = values;
  prop.ldata_u8 = sizes;
}
inline void add_ply_lists(ply_model& ply, const int* values, size_t count,
    int size, const string& element, const string& property) {
  if (!values) return;
  add_ply_property(ply, element, property, count, ply_type::i32, true);
  auto& prop = get_ply_property(ply, element, property);
  prop.data_i32.assign(values, values + count * size);
  prop.ldata_u8.assign(count, size);
}
inline void add_ply_lists(ply_model& ply, const vector<int>& values,
    const string& element, const string& property) {
  return add_ply_lists(ply, values.data(), values.size(), 1, element, property);
}
inline void add_ply_lists(ply_model& ply, const vector<vec2i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 2, element, property);
}
inline void add_ply_lists(ply_model& ply, const vector<vec3i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 3, element, property);
}
inline void add_ply_lists(ply_model& ply, const vector<vec4i>& values,
    const string& element, const string& property) {
  return add_ply_lists(
      ply, (int*)values.data(), values.size(), 4, element, property);
}

// Add ply properties for meshes
inline void add_ply_positions(ply_model& ply, const vector<vec3f>& values) {
  return add_ply_values(ply, values, "vertex", "x", "y", "z");
}
inline void add_ply_normals(ply_model& ply, const vector<vec3f>& values) {
  return add_ply_values(ply, values, "vertex", "nx", "ny", "nz");
}
inline void add_ply_texcoords(
    ply_model& ply, const vector<vec2f>& values, bool flipv) {
  return add_ply_values(
      ply, flipv ? flip_ply_texcoord(values) : values, "vertex", "u", "v");
}
inline void add_ply_colors(ply_model& ply, const vector<vec4f>& values) {
  return add_ply_values(ply, values, "vertex", "red", "green", "blue", "alpha");
}
inline void add_ply_radius(ply_model& ply, const vector<float>& values) {
  return add_ply_values(ply, values, "vertex", "radius");
}
inline void add_ply_faces(ply_model& ply, const vector<vector<int>>& values) {
  return add_ply_lists(ply, values, "face", "vertex_indices");
}
inline void add_ply_faces(ply_model& ply, const vector<vec3i>& triangles,
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
inline void add_ply_triangles(ply_model& ply, const vector<vec3i>& values) {
  return add_ply_faces(ply, values, {});
}
inline void add_ply_quads(ply_model& ply, const vector<vec4i>& values) {
  return add_ply_faces(ply, {}, values);
}
inline void add_ply_lines(ply_model& ply, const vector<vec2i>& values) {
  return add_ply_lists(ply, values, "line", "vertex_indices");
}
inline void add_ply_points(ply_model& ply, const vector<int>& values) {
  return add_ply_lists(ply, values, "point", "vertex_indices");
}

// get ply value either ascii or binary
template <typename T, typename VT>
inline void read_ply_prop(ply_file& fs, VT& value, bool big_endian) {
  auto tvalue = T{};
  read_ply_value(fs, tvalue, big_endian);
  value = (VT)tvalue;
}
template <typename VT>
inline void read_ply_prop(
    ply_file& fs, ply_type type, VT& value, bool big_endian) {
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
inline void parse_ply_prop(string_view& str, VT& value) {
  auto tvalue = T{};
  parse_ply_value(str, tvalue);
  value = (VT)tvalue;
}
template <typename VT>
inline void parse_ply_prop(string_view& str, ply_type type, VT& value) {
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
inline void read_ply_header(ply_file& fs, ply_format& format,
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
  while (read_ply_line(fs, buffer, sizeof(buffer))) {
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
      auto fmt = string_view{};
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
        parse_ply_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        auto itype = type_map.at(tname);
        if (itype != ply_type::u8)
          throw std::runtime_error{"unsupported list size type " + tname};
        parse_ply_value(line, tname);
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
      } else {
        prop.is_list = false;
        if (type_map.find(tname) == type_map.end())
          throw std::runtime_error{"unknown ply type " + tname};
        prop.type = type_map.at(tname);
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
inline void read_ply_value_generic(ply_file& fs, ply_format format,
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
    if (!read_ply_line(fs, buffer, sizeof(buffer)))
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
inline void format_ply_prop(ply_file& fs, ply_type type, VT value) {
  switch (type) {
    case ply_type::i8: format_ply_value(fs, (int8_t)value); break;
    case ply_type::i16: format_ply_value(fs, (int16_t)value); break;
    case ply_type::i32: format_ply_value(fs, (int32_t)value); break;
    case ply_type::i64: format_ply_value(fs, (int64_t)value); break;
    case ply_type::u8: format_ply_value(fs, (uint8_t)value); break;
    case ply_type::u16: format_ply_value(fs, (uint16_t)value); break;
    case ply_type::u32: format_ply_value(fs, (uint32_t)value); break;
    case ply_type::u64: format_ply_value(fs, (uint64_t)value); break;
    case ply_type::f32: format_ply_value(fs, (float)value); break;
    case ply_type::f64: format_ply_value(fs, (double)value); break;
  }
}

template <typename VT>
inline void write_ply_prop(
    ply_file& fs, ply_type type, VT value, bool big_endian) {
  switch (type) {
    case ply_type::i8: write_ply_value(fs, (int8_t)value, big_endian); break;
    case ply_type::i16: write_ply_value(fs, (int16_t)value, big_endian); break;
    case ply_type::i32: write_ply_value(fs, (int32_t)value, big_endian); break;
    case ply_type::i64: write_ply_value(fs, (int64_t)value, big_endian); break;
    case ply_type::u8: write_ply_value(fs, (uint8_t)value, big_endian); break;
    case ply_type::u16: write_ply_value(fs, (uint16_t)value, big_endian); break;
    case ply_type::u32: write_ply_value(fs, (uint32_t)value, big_endian); break;
    case ply_type::u64: write_ply_value(fs, (uint64_t)value, big_endian); break;
    case ply_type::f32: write_ply_value(fs, (float)value, big_endian); break;
    case ply_type::f64: write_ply_value(fs, (double)value, big_endian); break;
  }
}

// Write Ply functions
inline void write_ply_header(ply_file& fs, ply_format format,
    const vector<ply_element>& elements, const vector<string>& comments) {
  // ply type names
  static auto type_map = unordered_map<ply_type, string>{{ply_type::i8, "char"},
      {ply_type::i16, "short"}, {ply_type::i32, "int"}, {ply_type::i64, "uint"},
      {ply_type::u8, "uchar"}, {ply_type::u16, "ushort"},
      {ply_type::u32, "uint"}, {ply_type::u64, "ulong"},
      {ply_type::f32, "float"}, {ply_type::f64, "double"}};

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
        write_ply_text(fs, "property list uchar " + type_map[prop.type] + " " +
                               prop.name + "\n");
      } else {
        write_ply_text(
            fs, "property " + type_map[prop.type] + " " + prop.name + "\n");
      }
    }
  }
  write_ply_text(fs, "end_header\n");
}

template <typename VT, typename LT>
inline void write_ply_value_generic(ply_file& fs, ply_format format,
    const ply_element& element, vector<VT>& values, vector<vector<LT>>& lists) {
  if (format == ply_format::ascii) {
    for (auto pidx = 0; pidx < element.properties.size(); pidx++) {
      auto& prop = element.properties[pidx];
      if (pidx) format_ply_value(fs, " ");
      if (!prop.is_list) {
        format_ply_prop(fs, prop.type, values[pidx]);
      } else {
        format_ply_prop(fs, ply_type::u8, values[pidx]);
        for (auto i = 0; i < (int)lists[pidx].size(); i++) {
          if (i) format_ply_value(fs, " ");
          format_ply_prop(fs, prop.type, lists[pidx][i]);
        }
      }
      format_ply_value(fs, "\n");
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

inline void write_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  write_ply_value_generic(fs, format, element, values, lists);
}
inline void write_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  write_ply_value_generic(fs, format, element, values, lists);
}

inline void read_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists) {
  read_ply_value_generic(fs, format, element, values, lists);
}
inline void read_ply_value(ply_file& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists) {
  read_ply_value_generic(fs, format, element, values, lists);
}

inline int find_ply_element(
    const vector<ply_element>& elements, const string& name) {
  for (auto idx = 0; idx < elements.size(); idx++)
    if (elements[idx].name == name) return idx;
  return -1;
}
inline int find_ply_property(const ply_element& element, const string& name) {
  for (auto idx = 0; idx < element.properties.size(); idx++)
    if (element.properties[idx].name == name) return idx;
  return -1;
}
inline vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2) {
  auto ids = vec2i{
      find_ply_property(element, name1),
      find_ply_property(element, name2),
  };
  if (ids.x < 0 || ids.y < 0) return vec2i{-1};
  return ids;
}
inline vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3) {
  auto ids = vec3i{
      find_ply_property(element, name1),
      find_ply_property(element, name2),
      find_ply_property(element, name3),
  };
  if (ids.x < 0 || ids.y < 0 || ids.z < 0) return vec3i{-1};
  return ids;
}
inline vec4i find_ply_property(const ply_element& element, const string& name1,
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

#endif
