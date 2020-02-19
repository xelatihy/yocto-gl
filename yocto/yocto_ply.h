//
// # Yocto/Ply: Tiny library for Ply parsing and writing
//
// Yocto/Ply is a tiny library for loading and saving Ply. 
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

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::ply {

// Using directives
using namespace yocto::math;

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
  string                name       = "";
  size_t                count      = 0;
  vector<ply_property*> properties = {};
  ~ply_element();
};

// Ply model
struct ply_model {
  ply_format           format   = ply_format::binary_little_endian;
  vector<string>       comments = {};
  vector<ply_element*> elements = {};
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

vector<float> get_values(
    ply_model* ply, const string& element, const string& property);
vector<vec2f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2);
vector<vec3f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3);
vector<vec4f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      const string& property4);
vector<vec4f>   get_values(ply_model* ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      float property4);
vector<frame3f> get_values(
    ply_model* ply, const string& element, const array<string, 12>& properties);

vector<vector<int>> get_lists(
    ply_model* ply, const string& element, const string& property);
vector<byte> get_list_sizes(
    ply_model* ply, const string& element, const string& property);
vector<int> get_list_values(
    ply_model* ply, const string& element, const string& property);
vec2i get_list_minxmax(
    ply_model* ply, const string& element, const string& property);

// Get ply properties for meshes
vector<vec3f>       get_positions(ply_model* ply);
vector<vec3f>       get_normals(ply_model* ply);
vector<vec2f>       get_texcoords(ply_model* ply, bool flipv = false);
vector<vec3f>       get_colors(ply_model* ply);
vector<float>       get_radius(ply_model* ply);
vector<vector<int>> get_faces(ply_model* ply);
vector<vec2i>       get_lines(ply_model* ply);
vector<int>         get_points(ply_model* ply);
vector<vec3i>       get_triangles(ply_model* ply);
vector<vec4i>       get_quads(ply_model* ply);
bool                has_quads(ply_model* ply);

// Add ply properties
void add_values(ply_model* ply, const vector<float>& values,
    const string& element, const string& property);
void add_values(ply_model* ply, const vector<vec2f>& values,
    const string& element, const string& property1, const string& property2);
void add_values(ply_model* ply, const vector<vec3f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3);
void add_values(ply_model* ply, const vector<vec4f>& values,
    const string& element, const string& property1, const string& property2,
    const string& property3, const string& property4);
void add_values(ply_model* ply, const vector<frame3f>& values,
    const string& element, const array<string, 12>& properties);

void add_lists(ply_model* ply, const vector<vector<int>>& values,
    const string& element, const string& property);
void add_lists(ply_model* ply, const vector<byte>& sizes,
    const vector<int>& values, const string& element, const string& property);
void add_lists(ply_model* ply, const vector<int>& values, const string& element,
    const string& property);
void add_lists(ply_model* ply, const vector<vec2i>& values,
    const string& element, const string& property);
void add_lists(ply_model* ply, const vector<vec3i>& values,
    const string& element, const string& property);
void add_lists(ply_model* ply, const vector<vec4i>& values,
    const string& element, const string& property);

// Add ply properties for meshes
void add_positions(ply_model* ply, const vector<vec3f>& values);
void add_normals(ply_model* ply, const vector<vec3f>& values);
void add_texcoords(
    ply_model* ply, const vector<vec2f>& values, bool flipv = false);
void add_colors(ply_model* ply, const vector<vec3f>& values);
void add_radius(ply_model* ply, const vector<float>& values);
void add_faces(ply_model* ply, const vector<vector<int>>& values);
void add_faces(
    ply_model* ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
void add_triangles(ply_model* ply, const vector<vec3i>& values);
void add_quads(ply_model* ply, const vector<vec4i>& values);
void add_lines(ply_model* ply, const vector<vec2i>& values);
void add_points(ply_model* ply, const vector<int>& values);

}  // namespace yocto::ply

#endif
