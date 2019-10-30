//
// # Yocto/ModelIO: Tiny library for Ply/Obj/Pbrt/Yaml/glTF parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Ply/Obj/Pbrt/Yaml/glTF. 
// Yocto/ModelIO supports two interfaces: a simple interface where all model 
// data is loaded and saved at once and a low-level interface where single 
// commands values are read and written one at a time.
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

#ifndef _YOCTO_MODELIO_H_
#define _YOCTO_MODELIO_H_

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

// Result of io operations
struct plyio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

// Load and save ply
plyio_status load_ply(const string& filename, ply_model& ply);
plyio_status save_ply(const string& filename, const ply_model& ply);

// Get ply properties
bool has_ply_property(
    const ply_model& ply, const string& element, const string& property);
const ply_property& get_ply_property(
    const ply_model& ply, const string& element, const string& property);

vector<float> get_ply_values(
    const ply_model& ply, const string& element, const string& property);
vector<vec2f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2);
vector<vec3f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3);
vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    const string& property4);
vector<vec4f> get_ply_values(const ply_model& ply, const string& element,
    const string& property1, const string& property2, const string& property3,
    float property4);

vector<vector<int>> get_ply_lists(
    const ply_model& ply, const string& element, const string& property);
vector<byte> get_ply_list_sizes(
    const ply_model& ply, const string& element, const string& property);
vector<int> get_ply_list_values(
    const ply_model& ply, const string& element, const string& property);
vec2i get_ply_list_minxmax(
    const ply_model& ply, const string& element, const string& property);

// Get ply properties for meshes
vector<vec3f> get_ply_positions(const ply_model& ply);
vector<vec3f> get_ply_normals(const ply_model& ply);
vector<vec2f> get_ply_texcoords(
    const ply_model& ply, bool flipv = false);
vector<vec4f>       get_ply_colors(const ply_model& ply);
vector<float>       get_ply_radius(const ply_model& ply);
vector<vector<int>> get_ply_faces(const ply_model& ply);
vector<vec2i>       get_ply_lines(const ply_model& ply);
vector<int>         get_ply_points(const ply_model& ply);
vector<vec3i>       get_ply_triangles(const ply_model& ply);
vector<vec4i>       get_ply_quads(const ply_model& ply);
bool                has_ply_quads(const ply_model& ply);

// Add ply properties
void add_ply_values(ply_model& ply, const vector<float>& values,
    const string& element, const string& property);
void add_ply_values(ply_model& ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2);
void add_ply_values(ply_model& ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3);
void add_ply_values(ply_model& ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4);

void add_ply_lists(ply_model& ply, const vector<vector<int>>& values,
    const string& element, const string& property);
void add_ply_lists(ply_model& ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property);
void add_ply_lists(ply_model& ply, const vector<int>& values,
    const string& element, const string& property);
void add_ply_lists(ply_model& ply, const vector<vec2i>& values,
    const string& element, const string& property);
void add_ply_lists(ply_model& ply, const vector<vec3i>& values,
    const string& element, const string& property);
void add_ply_lists(ply_model& ply, const vector<vec4i>& values,
    const string& element, const string& property);

// Add ply properties for meshes
void add_ply_positions(ply_model& ply, const vector<vec3f>& values);
void add_ply_normals(ply_model& ply, const vector<vec3f>& values);
void add_ply_texcoords(
    ply_model& ply, const vector<vec2f>& values, bool flipv = false);
void add_ply_colors(ply_model& ply, const vector<vec4f>& values);
void add_ply_radius(ply_model& ply, const vector<float>& values);
void add_ply_faces(ply_model& ply, const vector<vector<int>>& values);
void add_ply_faces(
    ply_model& ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
void add_ply_triangles(ply_model& ply, const vector<vec3i>& values);
void add_ply_quads(ply_model& ply, const vector<vec4i>& values);
void add_ply_lines(ply_model& ply, const vector<vec2i>& values);
void add_ply_points(ply_model& ply, const vector<int>& values);

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

  operator bool() const { return (bool)fs; }

  FILE*  fs       = nullptr;
  string filename = "";
  string mode     = "rt";
  int    linenum  = 0;
};

// open a file
ply_file open_ply(const string& filename, const string& mode = "rt");
void     open_ply(
        ply_file& fs, const string& filename, const string& mode = "rt");
void close_ply(ply_file& fs);

// Read Ply functions
plyio_status read_ply_header(const string& filename, ply_file& fs,
    ply_format& format, vector<ply_element>& elements,
    vector<string>& comments);
plyio_status read_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
plyio_status read_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Write Ply functions
plyio_status write_ply_header(const string& filename, ply_file& fs,
    ply_format format, const vector<ply_element>& elements,
    const vector<string>& comments);
plyio_status write_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
plyio_status write_ply_value(const string& filename, ply_file& fs,
    ply_format format, const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Helpers to get element and property indices
int find_ply_element(
    const vector<ply_element>& elements, const string& name);
int   find_ply_property(const ply_element& element, const string& name);
vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2);
vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3);
vec4i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3, const string& name4);

}  // namespace yocto

#endif
