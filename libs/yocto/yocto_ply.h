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

#include "yocto_math.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::ply {

// Ply property
struct property {
  // property types
  enum struct type_t { i8, i16, i32, i64, u8, u16, u32, u64, f32, f64 };

  // description
  std::string name    = "";
  bool        is_list = false;
  type_t      type    = type_t::f32;

  // data
  std::vector<int8_t>   data_i8  = {};
  std::vector<int16_t>  data_i16 = {};
  std::vector<int32_t>  data_i32 = {};
  std::vector<int64_t>  data_i64 = {};
  std::vector<uint8_t>  data_u8  = {};
  std::vector<uint16_t> data_u16 = {};
  std::vector<uint32_t> data_u32 = {};
  std::vector<uint64_t> data_u64 = {};
  std::vector<float>    data_f32 = {};
  std::vector<double>   data_f64 = {};

  // list length
  std::vector<uint8_t> ldata_u8 = {};
};

// Ply elements
struct element {
  // element content
  std::string                 name       = "";
  size_t                      count      = 0;
  std::vector<ply::property*> properties = {};

  // cleanup
  ~element();
};

// Ply model
struct model {
  // Format type
  enum struct format_t { ascii, binary_little_endian, binary_big_endian };

  // ply content
  format_t                   format   = format_t::binary_little_endian;
  std::vector<std::string>   comments = {};
  std::vector<ply::element*> elements = {};

  // cleanup
  ~model();
};

// Load and save ply
bool load_ply(
    const std::string& filename, ply::model* ply, std::string& error);
bool save_ply(
    const std::string& filename, ply::model* ply, std::string& error);

// Get ply properties
bool has_property(
    ply::model* ply, const std::string& element, const std::string& property);
ply::property* get_property(
    ply::model* ply, const std::string& element, const std::string& property);

bool get_value(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<float>& values);
bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties, std::vector<vec4f>& values);
bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties, std::vector<vec3f>& values);
bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties, std::vector<vec4f>& values);
bool get_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    std::vector<frame3f>&              values);

bool get_lists(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<std::vector<int>>& lists);
bool get_list_sizes(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<byte>& sizes);
bool get_list_values(ply::model* ply, const std::string& element,
    const std::string& property, std::vector<int>& values);

// Get ply properties for meshes
bool get_positions(ply::model* ply, std::vector<vec3f>& values);
bool get_normals(ply::model* ply, std::vector<vec3f>& values);
bool get_texcoords(
    ply::model* ply, std::vector<vec2f>& values, bool flipv = false);
bool get_colors(ply::model* ply, std::vector<vec3f>& values);
bool get_radius(ply::model* ply, std::vector<float>& values);
bool get_faces(ply::model* ply, std::vector<std::vector<int>>*& values);
bool get_lines(ply::model* ply, std::vector<vec2i>& values);
bool get_points(ply::model* ply, std::vector<int>& values);
bool get_triangles(ply::model* ply, std::vector<vec3i>& values);
bool get_quads(ply::model* ply, std::vector<vec4i>& values);
bool has_quads(ply::model* ply);

// Add ply properties
bool add_value(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<float>& values);
bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 2>& properties,
    const std::vector<vec2f>&         values);
bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 3>& properties,
    const std::vector<vec3f>&         values);
bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 4>& properties,
    const std::vector<vec4f>&         values);
bool add_values(ply::model* ply, const std::string& element,
    const std::array<std::string, 12>& properties,
    const std::vector<frame3f>&        values);

bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<std::vector<int>>& values);
bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<byte>& sizes,
    const std::vector<int>& values);
bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<int>& values);
bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec2i>& values);
bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec3i>& values);
bool add_lists(ply::model* ply, const std::string& element,
    const std::string& property, const std::vector<vec4i>& values);

// Add ply properties for meshes
bool add_positions(ply::model* ply, const std::vector<vec3f>& values);
bool add_normals(ply::model* ply, const std::vector<vec3f>& values);
bool add_texcoords(
    ply::model* ply, const std::vector<vec2f>& values, bool flipv = false);
bool add_colors(ply::model* ply, const std::vector<vec3f>& values);
bool add_radius(ply::model* ply, const std::vector<float>& values);
bool add_faces(
    ply::model* ply, const std::vector<std::vector<int>>& values);
bool add_faces(ply::model* ply, const std::vector<vec3i>& tvalues,
    const std::vector<vec4i>& qvalues);
bool add_triangles(ply::model* ply, const std::vector<vec3i>& values);
bool add_quads(ply::model* ply, const std::vector<vec4i>& values);
bool add_lines(ply::model* ply, const std::vector<vec2i>& values);
bool add_points(ply::model* ply, const std::vector<int>& values);

}  // namespace yocto::ply

#endif
