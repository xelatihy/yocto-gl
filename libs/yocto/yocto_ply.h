//
// # Yocto/Ply: Tiny library for Ply parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Ply.
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

#ifndef _YOCTO_PLY_H_
#define _YOCTO_PLY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Ply type
enum struct ply_type { i8, i16, i32, i64, u8, u16, u32, u64, f32, f64 };

// Ply property
struct ply_property {
  // description
  string   name    = "";
  bool     is_list = false;
  ply_type type    = ply_type::f32;

  // data
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
  // element content
  string                name       = "";
  size_t                count      = 0;
  vector<ply_property*> properties = {};

  // cleanup
  ~ply_element();
};

// Ply format
enum struct ply_format { ascii, binary_little_endian, binary_big_endian };

// Ply model
struct ply_model {
  // ply content
  ply_format           format   = ply_format::binary_little_endian;
  vector<string>       comments = {};
  vector<ply_element*> elements = {};

  // cleanup
  ~ply_model();
};

// Load and save ply
bool load_ply(const string& filename, ply_model* ply, string& error);
bool save_ply(const string& filename, ply_model* ply, string& error);

// Get ply properties
bool has_property(
    ply_model* ply, const string& element, const string& property);
ply_property* get_property(
    ply_model* ply, const string& element, const string& property);

bool get_value(ply_model* ply, const string& element, const string& property,
    vector<float>& values);
bool get_values(ply_model* ply, const string& element,
    const std::array<string, 2>& properties, vector<vec4f>& values);
bool get_values(ply_model* ply, const string& element,
    const std::array<string, 3>& properties, vector<vec3f>& values);
bool get_values(ply_model* ply, const string& element,
    const std::array<string, 4>& properties, vector<vec4f>& values);
bool get_values(ply_model* ply, const string& element,
    const std::array<string, 12>& properties, vector<frame3f>& values);

bool get_lists(ply_model* ply, const string& element, const string& property,
    vector<vector<int>>& lists);
bool get_list_sizes(ply_model* ply, const string& element,
    const string& property, vector<byte>& sizes);
bool get_list_values(ply_model* ply, const string& element,
    const string& property, vector<int>& values);

// Get ply properties for meshes
bool get_positions(ply_model* ply, vector<vec3f>& values);
bool get_normals(ply_model* ply, vector<vec3f>& values);
bool get_texcoords(ply_model* ply, vector<vec2f>& values, bool flipv = false);
bool get_colors(ply_model* ply, vector<vec3f>& values);
bool get_radius(ply_model* ply, vector<float>& values);
bool get_faces(ply_model* ply, vector<vector<int>>*& values);
bool get_lines(ply_model* ply, vector<vec2i>& values);
bool get_points(ply_model* ply, vector<int>& values);
bool get_triangles(ply_model* ply, vector<vec3i>& values);
bool get_quads(ply_model* ply, vector<vec4i>& values);
bool has_quads(ply_model* ply);

// Add ply properties
bool add_value(ply_model* ply, const string& element, const string& property,
    const vector<float>& values);
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 2>& properties, const vector<vec2f>& values);
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 3>& properties, const vector<vec3f>& values);
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 4>& properties, const vector<vec4f>& values);
bool add_values(ply_model* ply, const string& element,
    const std::array<string, 12>& properties, const vector<frame3f>& values);

bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vector<int>>& values);
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<byte>& sizes, const vector<int>& values);
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<int>& values);
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec2i>& values);
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec3i>& values);
bool add_lists(ply_model* ply, const string& element, const string& property,
    const vector<vec4i>& values);

// Add ply properties for meshes
bool add_positions(ply_model* ply, const vector<vec3f>& values);
bool add_normals(ply_model* ply, const vector<vec3f>& values);
bool add_texcoords(
    ply_model* ply, const vector<vec2f>& values, bool flipv = false);
bool add_colors(ply_model* ply, const vector<vec3f>& values);
bool add_radius(ply_model* ply, const vector<float>& values);
bool add_faces(ply_model* ply, const vector<vector<int>>& values);
bool add_faces(
    ply_model* ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
bool add_triangles(ply_model* ply, const vector<vec3i>& values);
bool add_quads(ply_model* ply, const vector<vec4i>& values);
bool add_lines(ply_model* ply, const vector<vec2i>& values);
bool add_points(ply_model* ply, const vector<int>& values);

}  // namespace yocto

#endif
