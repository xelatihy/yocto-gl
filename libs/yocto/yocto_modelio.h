//
// # Yocto/ModelIO: Serialization for Obj, Ply, Stl and Pbrt models
//
// Yocto/ModelIO is a collection of utilities for loading and saving scenes
// and meshes in Ply, Obj, Stl and Pbrt formats.
// Yocto/ModelIO is implemented in `yocto_modelio.h` and `yocto_modelio.cpp`,
// and depends on `fast_float.h` for number parsing.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
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
  string               name       = "";
  size_t               count      = 0;
  vector<ply_property> properties = {};
};

// Ply format
enum struct ply_format { ascii, binary_little_endian, binary_big_endian };

// Ply model
struct ply_model {
  // ply content
  ply_format          format   = ply_format::binary_little_endian;
  vector<string>      comments = {};
  vector<ply_element> elements = {};
};

// Load and save ply
[[nodiscard]] bool load_ply(
    const string& filename, ply_model& ply, string& error);
[[nodiscard]] bool save_ply(
    const string& filename, const ply_model& ply, string& error);

// Get ply properties
inline bool has_property(
    const ply_model& ply, const string& element, const string& property);

// Get values
template <typename T>
inline bool get_value(const ply_model& ply, const string& element,
    const string& property, vector<T>& values);
template <typename T, size_t N>
inline bool get_value(const ply_model& ply, const string& element,
    const array<string, N>& properties, vector<array<T, N>>& values);

template <typename T>
inline bool get_lists(const ply_model& ply, const string& element,
    const string& property, vector<vector<T>>& lists);
template <typename T>
inline bool get_list_sizes(const ply_model& ply, const string& element,
    const string& property, vector<T>& sizes);
template <typename T>
inline bool get_list_values(const ply_model& ply, const string& element,
    const string& property, vector<T>& values);

// Get ply properties for meshes
template <typename T>
inline bool get_positions(const ply_model& ply, vector<array<T, 3>>& values);
template <typename T>
inline bool get_normals(const ply_model& ply, vector<array<T, 3>>& values);
template <typename T>
inline bool get_texcoords(
    const ply_model& ply, vector<array<T, 2>>& values, bool flipv = false);
template <typename T>
inline bool get_colors(const ply_model& ply, vector<array<T, 3>>& values);
template <typename T>
inline bool get_colors(const ply_model& ply, vector<array<T, 4>>& values);
template <typename T>
inline bool get_radius(const ply_model& ply, vector<T>& values);
template <typename T>
inline bool get_faces(const ply_model& ply, vector<vector<T>>& faces);
template <typename T>
inline bool get_lines(const ply_model& ply, vector<array<T, 2>>& lines);
template <typename T>
inline bool get_points(const ply_model& ply, vector<T>& points);
template <typename T>
inline bool get_triangles(const ply_model& ply, vector<array<T, 3>>& triangles);
template <typename T>
inline bool get_quads(const ply_model& ply, vector<array<T, 4>>& quads);
template <typename T>
inline bool get_faces(const ply_model& ply, vector<array<T, 3>>& triangles,
    vector<array<T, 4>>& quads);
inline bool has_quads(const ply_model& ply);

// Add ply properties
template <typename T>
inline bool add_value(ply_model& ply, const string& element,
    const string& property, const vector<T>& values);
template <typename T, size_t N>
inline bool add_values(ply_model& ply, const string& element,
    const array<string, N>& properties, const vector<array<T, N>>& values);

template <typename T>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<vector<T>>& values);
template <typename T>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<uint8_t>& sizes,
    const vector<T>& values);
template <typename T, size_t N>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<array<T, N>>& values);

// Add ply properties for meshes
template <typename T>
inline bool add_positions(ply_model& ply, const vector<array<T, 3>>& values);
template <typename T>
inline bool add_normals(ply_model& ply, const vector<array<T, 3>>& values);
template <typename T>
inline bool add_texcoords(
    ply_model& ply, const vector<array<T, 2>>& values, bool flipv = false);
template <typename T>
inline bool add_colors(ply_model& ply, const vector<array<T, 3>>& values);
template <typename T>
inline bool add_colors(ply_model& ply, const vector<array<T, 4>>& values);
template <typename T>
inline bool add_radius(ply_model& ply, const vector<T>& values);
template <typename T>
inline bool add_faces(ply_model& ply, const vector<vector<T>>& values);
template <typename T>
inline bool add_faces(ply_model& ply, const vector<array<T, 3>>& triangles,
    const vector<array<T, 4>>& quads);
template <typename T>
inline bool add_triangles(ply_model& ply, const vector<array<T, 3>>& triangles);
template <typename T>
inline bool add_quads(ply_model& ply, const vector<array<T, 4>>& quads);
template <typename T>
inline bool add_lines(ply_model& ply, const vector<array<T, 2>>& lines);
template <typename T>
inline bool add_points(ply_model& ply, const vector<T>& points);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Obj vertex
struct obj_vertex {
  int position = 0;
  int texcoord = 0;
  int normal   = 0;
};

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
  return a.position == b.position && a.texcoord == b.texcoord &&
         a.normal == b.normal;
}

// Obj element type
enum struct obj_etype : uint16_t { face, line, point };

// Obj element
struct obj_element {
  uint16_t  size     = 0;
  obj_etype etype    = obj_etype::face;
  int       material = 0;
};

// Obj texture information.
struct obj_texture {
  string path  = "";     // file path
  bool   clamp = false;  // clamp to edge
  float  scale = 1;      // scale for bump/displacement

  obj_texture() = default;
  explicit obj_texture(const string& path) : path{path} {}
};

// Obj material
struct obj_material {
  // material name and type
  string name  = "";
  int    illum = 0;

  // material colors and values
  array<float, 3> emission     = {0, 0, 0};
  array<float, 3> ambient      = {0, 0, 0};
  array<float, 3> diffuse      = {0, 0, 0};
  array<float, 3> specular     = {0, 0, 0};
  array<float, 3> reflection   = {0, 0, 0};
  array<float, 3> transmission = {0, 0, 0};
  float           exponent     = 10;
  float           ior          = 1.5;
  float           opacity      = 1;

  // material textures
  int emission_tex     = -1;
  int ambient_tex      = -1;
  int diffuse_tex      = -1;
  int specular_tex     = -1;
  int reflection_tex   = -1;
  int transmission_tex = -1;
  int exponent_tex     = -1;
  int opacity_tex      = -1;
  int bump_tex         = -1;
  int normal_tex       = -1;
  int displacement_tex = -1;
};

// Obj shape
struct obj_shape {
  string                  name      = "";
  vector<array<float, 3>> positions = {};
  vector<array<float, 3>> normals   = {};
  vector<array<float, 2>> texcoords = {};
  vector<obj_vertex>      vertices  = {};
  vector<obj_element>     elements  = {};
};

// Obj camera
struct obj_camera {
  string           name     = "";
  array<float, 12> frame    = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  bool             ortho    = false;
  float            aspect   = 16.0f / 9.0f;
  float            lens     = 0.50f;
  float            film     = 0.036f;
  float            focus    = 0;
  float            aperture = 0;
};

// Obj environment
struct obj_environment {
  string           name         = "";
  array<float, 12> frame        = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 3>  emission     = {0, 0, 0};
  int              emission_tex = -1;
};

// Obj model
struct obj_model {
  vector<string>          comments     = {};
  vector<obj_shape>       shapes       = {};
  vector<obj_material>    materials    = {};
  vector<obj_texture>     textures     = {};
  vector<obj_camera>      cameras      = {};
  vector<obj_environment> environments = {};
};

// Load and save obj shape
obj_shape load_sobj(const string& filename, bool face_varying = false);
void      load_obj(
         const string& filename, obj_shape& obj, bool face_varying = false);
void save_obj(const string& filename, const obj_shape& obj);

// Load and save obj
[[nodiscard]] bool load_obj(const string& filename, obj_model& obj,
    string& error, bool face_varying = false, bool split_materials = false);
[[nodiscard]] bool save_obj(
    const string& filename, const obj_model& obj, string& error);

// Load and save obj shape
[[nodiscard]] bool load_obj(const string& filename, obj_shape& obj,
    string& error, bool face_varying = false);
[[nodiscard]] bool save_obj(
    const string& filename, const obj_shape& obj, string& error);

// Get obj shape.
void get_positions(const obj_shape& obj, vector<array<float, 3>>& positions);
void get_normals(const obj_shape& obj, vector<array<float, 3>>& normals);
void get_texcoords(const obj_shape& obj, vector<array<float, 2>>& texcoords,
    bool flipv = false);
void get_faces(const obj_shape& obj, vector<array<int, 3>>& triangles,
    vector<array<int, 4>>& quads, vector<int>& materials);
void get_triangles(const obj_shape& obj, vector<array<int, 3>>& triangles,
    vector<int>& materials);
void get_quads(
    const obj_shape& obj, vector<array<int, 4>>& quads, vector<int>& materials);
void get_lines(
    const obj_shape& obj, vector<array<int, 2>>& lines, vector<int>& materials);
void get_points(
    const obj_shape& obj, vector<int>& points, vector<int>& materials);
void get_fvquads(const obj_shape& obj, vector<array<int, 4>>& quadspos,
    vector<array<int, 4>>& quadsnorm, vector<array<int, 4>>& quadstexcoord,
    vector<int>& materials);
void get_faces(const obj_shape& obj, int material,
    vector<array<int, 3>>& triangles, vector<array<int, 4>>& quads);
void get_triangles(
    const obj_shape& obj, int material, vector<array<int, 3>>& triangles);
void get_quads(
    const obj_shape& obj, int material, vector<array<int, 4>>& quads);
void get_lines(
    const obj_shape& obj, int material, vector<array<int, 2>>& lines);
void get_points(const obj_shape& obj, int material, vector<int>& points);
bool has_quads(const obj_shape& obj);

// get unique materials from shape
vector<int> get_materials(const obj_shape& obj);

// Add obj shape
void add_positions(obj_shape& obj, const vector<array<float, 3>>& positions);
void add_normals(obj_shape& obj, const vector<array<float, 3>>& normals);
void add_texcoords(obj_shape& obj, const vector<array<float, 2>>& texcoords,
    bool flipv = false);
void add_triangles(obj_shape& obj, const vector<array<int, 3>>& triangles,
    int material, bool has_normals, bool has_texcoord);
void add_quads(obj_shape& obj, const vector<array<int, 4>>& quads, int material,
    bool has_normals, bool has_texcoord);
void add_lines(obj_shape& obj, const vector<array<int, 2>>& lines, int material,
    bool has_normals, bool has_texcoord);
void add_points(obj_shape& obj, const vector<int>& points, int material,
    bool has_normals, bool has_texcoord);
void add_fvquads(obj_shape& obj, const vector<array<int, 4>>& quadspos,
    const vector<array<int, 4>>& quadsnorm,
    const vector<array<int, 4>>& quadstexcoord, int material);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::obj_vertex> {
  size_t operator()(const yocto::obj_vertex& v) const {
    const std::hash<int> hasher = std::hash<int>();
    auto                 h      = (size_t)0;
    h ^= hasher(v.position) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.normal) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.texcoord) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

struct stl_shape {
  vector<array<float, 3>> positions = {};
  vector<array<int, 3>>   triangles = {};
  vector<array<float, 3>> fnormals  = {};
};

struct stl_model {
  vector<stl_shape> shapes = {};
};

// Load/save stl
[[nodiscard]] bool load_stl(const string& filename, stl_model& stl,
    string& error, bool unique_vertices = true);
[[nodiscard]] bool save_stl(const string& filename, const stl_model& stl,
    string& error, bool ascii = false);

// Get/set data
bool get_triangles(const stl_model& stl, int shape_id,
    vector<array<int, 3>>& triangles, vector<array<float, 3>>& positions,
    vector<array<float, 3>>& fnormals);
void add_triangles(stl_model& stl, const vector<array<int, 3>>& triangles,
    const vector<array<float, 3>>& positions,
    const vector<array<float, 3>>& fnormals);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  array<float, 12> frame      = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12> frend      = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<int, 2>    resolution = {0, 0};
  float            lens       = 0;
  float            aspect     = 0;
  float            focus      = 0;
  float            aperture   = 0;
};

// Pbrt material
struct pbrt_texture {
  string          name     = "";
  array<float, 3> constant = {1, 1, 1};
  string          filename = "";
};

// Pbrt material type (simplified and only for the materials that matter here)
enum struct pbrt_mtype {
  // clang-format off
  matte, plastic, metal, glass, thinglass, subsurface
  // clang-format on
};

// Pbrt material
struct pbrt_material {
  string          name            = "";
  pbrt_mtype      type            = pbrt_mtype::matte;
  array<float, 3> emission        = {0, 0, 0};
  array<float, 3> color           = {0, 0, 0};
  float           roughness       = 0;
  float           ior             = 1.5f;
  float           opacity         = 1;
  int             color_tex       = -1;
  array<float, 3> volmeanfreepath = {0, 0, 0};
  array<float, 3> volscatter      = {0, 0, 0};
  float           volscale        = 0.01f;
};

// Pbrt shape
struct pbrt_shape {
  array<float, 12>         frame     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>         frend     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  bool                     instanced = false;
  vector<array<float, 12>> instances = {};
  vector<array<float, 12>> instaends = {};
  int                      material  = -1;
  string                   filename_ = "";
  vector<array<float, 3>>  positions = {};
  vector<array<float, 3>>  normals   = {};
  vector<array<float, 2>>  texcoords = {};
  vector<array<int, 3>>    triangles = {};
};

// Pbrt lights
struct pbrt_light {
  array<float, 12>        frame          = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>        frend          = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 3>         emission       = {0, 0, 0};
  array<float, 3>         from           = {0, 0, 0};
  array<float, 3>         to             = {0, 0, 0};
  bool                    distant        = false;
  array<float, 3>         area_emission  = {0, 0, 0};
  array<float, 12>        area_frame     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12>        area_frend     = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  vector<array<int, 3>>   area_triangles = {};
  vector<array<float, 3>> area_positions = {};
  vector<array<float, 3>> area_normals   = {};
};
struct pbrt_environment {
  array<float, 12> frame        = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 12> frend        = {1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0};
  array<float, 3>  emission     = {0, 0, 0};
  int              emission_tex = -1;
};

// Pbrt model
struct pbrt_model {
  // pbrt data
  vector<string>           comments     = {};
  vector<pbrt_camera>      cameras      = {};
  vector<pbrt_shape>       shapes       = {};
  vector<pbrt_environment> environments = {};
  vector<pbrt_light>       lights       = {};
  vector<pbrt_material>    materials    = {};
  vector<pbrt_texture>     textures     = {};
};

// Load/save pbrt
[[nodiscard]] bool load_pbrt(const string& filename, pbrt_model& pbrt,
    string& error, bool ply_meshes = false);
[[nodiscard]] bool save_pbrt(const string& filename, const pbrt_model& pbrt,
    string& error, bool ply_meshes = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Get ply properties
inline bool has_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return true;
    }
  }
  return false;
}
inline ply_property& get_property(
    ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
inline const ply_property& get_property(
    const ply_model& ply, const string& element, const string& property) {
  for (auto& elem : ply.elements) {
    if (elem.name != element) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property) return prop;
    }
  }
  throw std::runtime_error("property not found");
}
inline size_t get_size(const ply_property& prop) {
  switch (prop.type) {
    case ply_type::i8: return prop.data_i8.size();
    case ply_type::i16: return prop.data_i16.size();
    case ply_type::i32: return prop.data_i32.size();
    case ply_type::i64: return prop.data_i64.size();
    case ply_type::u8: return prop.data_u8.size();
    case ply_type::u16: return prop.data_u16.size();
    case ply_type::u32: return prop.data_u32.size();
    case ply_type::u64: return prop.data_u64.size();
    case ply_type::f32: return prop.data_f32.size();
    case ply_type::f64: return prop.data_f64.size();
  }
  return 0;
}
template <typename T>
inline bool get_value(const ply_property& prop, size_t index, T& value) {
  switch (prop.type) {
    case ply_type::i8: value = (T)prop.data_i8[index]; return true;
    case ply_type::i16: value = (T)prop.data_i16[index]; return true;
    case ply_type::i32: value = (T)prop.data_i32[index]; return true;
    case ply_type::i64: value = (T)prop.data_i64[index]; return true;
    case ply_type::u8: value = (T)prop.data_u8[index]; return true;
    case ply_type::u16: value = (T)prop.data_u16[index]; return true;
    case ply_type::u32: value = (T)prop.data_u32[index]; return true;
    case ply_type::u64: value = (T)prop.data_u64[index]; return true;
    case ply_type::f32: value = (T)prop.data_f32[index]; return true;
    case ply_type::f64: value = (T)prop.data_f64[index]; return true;
  }
  return false;
}
template <typename T>
inline T get_value(const ply_property& prop, size_t index) {
  switch (prop.type) {
    case ply_type::i8: return (T)prop.data_i8[index];
    case ply_type::i16: return (T)prop.data_i16[index];
    case ply_type::i32: return (T)prop.data_i32[index];
    case ply_type::i64: return (T)prop.data_i64[index];
    case ply_type::u8: return (T)prop.data_u8[index];
    case ply_type::u16: return (T)prop.data_u16[index];
    case ply_type::u32: return (T)prop.data_u32[index];
    case ply_type::u64: return (T)prop.data_u64[index];
    case ply_type::f32: return (T)prop.data_f32[index];
    case ply_type::f64: return (T)prop.data_f64[index];
  }
  return 0;
}
template <typename T>
inline bool get_value(const ply_model& ply, const string& element,
    const string& property, vector<T>& values) {
  values.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (prop.is_list) return false;
  values.resize(get_size(prop));
  for (auto index = (size_t)0; index < values.size(); index++) {
    values[index] = get_value<T>(prop, index);
  }
  return true;
}
template <typename T, size_t N>
inline bool get_values(const ply_model& ply, const string& element,
    const array<string, N>& properties, vector<array<T, N>>& values) {
  values.clear();
  for (auto& property : properties) {
    if (!has_property(ply, element, property)) return false;
    auto& prop = get_property(ply, element, property);
    if (prop.is_list) return false;
  }
  values.resize(get_size(get_property(ply, element, properties.front())));
  auto item = (size_t)0;
  for (auto& property : properties) {
    auto& prop = get_property(ply, element, property);
    for (auto index = (size_t)0; index < values.size(); index++) {
      values[index][item] = get_value<T>(prop, index);
    }
    item++;
  }
  return true;
}

template <typename T>
inline bool get_lists(const ply_model& ply, const string& element,
    const string& property, vector<vector<T>>& lists) {
  lists.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes = prop.ldata_u8;
  lists.resize(sizes.size());
  auto list = (size_t)0, current = (size_t)0;
  for (auto size : sizes) {
    lists[list].resize(size);
    for (auto item = (size_t)0; item < size; item++)
      lists[list][item] = get_value<T>(prop, current + item);
    list += 1;
    current += size;
  }
  return true;
}
template <typename T>
inline bool get_list_sizes(const ply_model& ply, const string& element,
    const string& property, vector<T>& sizes) {
  sizes.clear();
  if (!has_property(ply, element, property)) return {};
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  if constexpr (std::is_same_v<T, uint8_t>) {
    sizes = prop.ldata_u8;
  } else {
    sizes.resize(prop.ldata_u8.size());
    for (auto index = (size_t)0; index < sizes.size(); index++) {
      sizes[index] = (T)prop.ldata_u8[index];
    }
  }
  return true;
}
template <typename T>
inline bool get_list_values(const ply_model& ply, const string& element,
    const string& property, vector<T>& values) {
  values.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  values.resize(get_size(prop));
  for (auto index = (size_t)0; index < values.size(); index++) {
    values[index] = get_value<T>(prop, index);
  }
  return true;
}
template <typename T>
inline bool get_triangles(const ply_model& ply, const string& element,
    const string& property, vector<array<T, 3>>& triangles) {
  triangles.clear();
  if (!has_property(ply, element, property)) return {};
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes = prop.ldata_u8;
  triangles.clear();
  triangles.reserve(sizes.size());
  auto current = (size_t)0;
  for (auto size : sizes) {
    if (size == 0) {
      triangles.push_back({(T)-1, (T)-1, (T)-1});
    } else if (size == 1) {
      triangles.push_back({get_value<T>(prop, current + 0), (T)-1, (T)-1});
    } else if (size == 2) {
      triangles.push_back({get_value<T>(prop, current + 0),
          get_value<T>(prop, current + 1), (T)-1});
    } else if (size == 3) {
      triangles.push_back({get_value<T>(prop, current + 0),
          get_value<T>(prop, current + 1), get_value<T>(prop, current + 2)});
    } else {
      for (auto item = (size_t)2; item < size; item++) {
        triangles.push_back({get_value<T>(prop, current + 0),
            get_value<T>(prop, current + item - 1),
            get_value<T>(prop, current + item)});
      }
    }
    current += size;
  }
  return true;
}
template <typename T>
inline bool get_quads(const ply_model& ply, const string& element,
    const string& property, vector<array<T, 4>>& quads) {
  quads.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes = prop.ldata_u8;
  quads.clear();
  quads.reserve(sizes.size());
  auto current = (size_t)0;
  for (auto size : sizes) {
    if (size == 0) {
      quads.push_back({(T)-1, (T)-1, (T)-1, (T)-1});
    } else if (size == 1) {
      quads.push_back({get_value<T>(prop, current + 0), (T)-1, (T)-1, (T)-1});
    } else if (size == 2) {
      quads.push_back({get_value<T>(prop, current + 0),
          get_value<T>(prop, current + 1), (T)-1, (T)-1});
    } else if (size == 3) {
      quads.push_back({get_value<T>(prop, current + 0),
          get_value<T>(prop, current + 1), get_value<T>(prop, current + 2),
          get_value<T>(prop, current + 2)});
    } else if (size == 4) {
      quads.push_back({get_value<T>(prop, current + 0),
          get_value<T>(prop, current + 1), get_value<T>(prop, current + 2),
          get_value<T>(prop, current + 3)});
    } else {
      for (auto item = (size_t)2; item < size; item++) {
        quads.push_back({get_value<T>(prop, current + 0),
            get_value<T>(prop, current + item - 1),
            get_value<T>(prop, current + item),
            get_value<T>(prop, current + item)});
      }
    }
    current += size;
  }
  return true;
}
inline bool has_quads(
    const ply_model& ply, const string& element, const string& property) {
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes = prop.ldata_u8;
  for (auto size : sizes)
    if (size == 4) return true;
  return false;
}
template <typename T>
inline bool get_faces(const ply_model& ply, const string& element,
    const string& property, vector<array<T, 3>>& triangles,
    vector<array<T, 4>>& quads) {
  if (has_quads(ply, element, property)) {
    return get_quads(ply, element, property, quads);
  } else {
    return get_triangles(ply, element, property, triangles);
  }
}
template <typename T>
inline bool get_lines(const ply_model& ply, const string& element,
    const string& property, vector<array<T, 2>>& lines) {
  lines.clear();
  if (!has_property(ply, element, property)) return false;
  auto& prop = get_property(ply, element, property);
  if (!prop.is_list) return false;
  auto& sizes = prop.ldata_u8;
  lines.clear();
  lines.reserve(sizes.size());
  auto current = (size_t)0;
  for (auto size : sizes) {
    if (size == 0) {
      lines.push_back({(T)-1, (T)-1});
    } else if (size == 1) {
      lines.push_back({get_value<T>(prop, current + 0), (T)-1});
    } else if (size == 2) {
      lines.push_back(
          {get_value<T>(prop, current + 0), get_value<T>(prop, current + 1)});
    } else {
      for (auto item = (size_t)1; item < size; item++) {
        lines.push_back({get_value<T>(prop, current + item - 1),
            get_value<T>(prop, current + item)});
      }
    }
    current += size;
  }
  return true;
}
template <typename T>
inline bool get_points(const ply_model& ply, const string& element,
    const string& property, vector<T>& values) {
  return get_list_values(ply, element, property, values);
}

// Get ply properties for meshes
template <typename T>
inline bool get_positions(
    const ply_model& ply, vector<array<T, 3>>& positions) {
  return get_values(ply, "vertex", {"x", "y", "z"}, positions);
}
template <typename T>
inline bool get_normals(const ply_model& ply, vector<array<T, 3>>& normals) {
  return get_values(ply, "vertex", {"nx", "ny", "nz"}, normals);
}
template <typename T>
inline bool get_texcoords(
    const ply_model& ply, vector<array<T, 2>>& texcoords, bool flipv) {
  if (has_property(ply, "vertex", "u")) {
    if (!get_values(ply, "vertex", {"u", "v"}, texcoords)) return false;
  } else {
    if (!get_values(ply, "vertex", {"s", "t"}, texcoords)) return false;
  }
  if (flipv) {
    for (auto& uv : texcoords) uv = {uv[0], 1 - uv[1]};
  }
  return true;
}
template <typename T>
inline bool get_colors(const ply_model& ply, vector<array<T, 3>>& colors) {
  return get_values(ply, "vertex", {"red", "green", "blue"}, colors);
}
template <typename T>
inline bool get_colors(const ply_model& ply, vector<array<T, 4>>& colors) {
  if (has_property(ply, "vertex", "alpha")) {
    return get_values(ply, "vertex", {"red", "green", "blue", "alpha"}, colors);
  } else {
    auto colors3 = vector<array<T, 3>>{};
    if (!get_values(ply, "vertex", {"red", "green", "blue"}, colors3))
      return false;
    colors.resize(colors3.size());
    for (auto i = 0; i < colors.size(); i++)
      colors[i] = {colors3[i][0], colors3[i][1], colors3[i][2], 1};
    return true;
  }
}
template <typename T>
inline bool get_radius(const ply_model& ply, vector<T>& radius) {
  return get_value(ply, "vertex", "radius", radius);
}
template <typename T>
inline bool get_faces(const ply_model& ply, vector<vector<T>>& faces) {
  return get_lists(ply, "face", "vertex_indices", faces);
}
template <typename T>
inline bool get_triangles(
    const ply_model& ply, vector<array<T, 3>>& triangles) {
  return get_triangles(ply, "face", "vertex_indices", triangles);
}
template <typename T>
inline bool get_quads(const ply_model& ply, vector<array<T, 4>>& quads) {
  return get_quads(ply, "face", "vertex_indices", quads);
}
template <typename T>
inline bool get_faces(const ply_model& ply, vector<array<T, 3>>& triangles,
    vector<array<T, 4>>& quads) {
  return get_faces(ply, "face", "vertex_indices", triangles, quads);
}
template <typename T>
inline bool get_lines(const ply_model& ply, vector<array<T, 2>>& lines) {
  return get_lines(ply, "line", "vertex_indices", lines);
}
template <typename T>
inline bool get_points(const ply_model& ply, vector<T>& points) {
  return get_points(ply, "point", "vertex_indices", points);
}
inline bool has_quads(const ply_model& ply) {
  return has_quads(ply, "face", "vertex_indices");
}

// Add ply properties
template <typename T>
inline ply_type get_ply_type() {
  static_assert(std::is_arithmetic_v<T>, "not a supported type");
  if constexpr (std::is_same_v<T, int8_t>) return ply_type::i8;
  if constexpr (std::is_same_v<T, int16_t>) return ply_type::i16;
  if constexpr (std::is_same_v<T, int32_t>) return ply_type::i32;
  if constexpr (std::is_same_v<T, int64_t>) return ply_type::i64;
  if constexpr (std::is_same_v<T, uint8_t>) return ply_type::u8;
  if constexpr (std::is_same_v<T, uint16_t>) return ply_type::u16;
  if constexpr (std::is_same_v<T, uint32_t>) return ply_type::u32;
  if constexpr (std::is_same_v<T, uint64_t>) return ply_type::u64;
  if constexpr (std::is_same_v<T, float>) return ply_type::f32;
  if constexpr (std::is_same_v<T, double>) return ply_type::f64;
}
inline bool add_element(
    ply_model& ply, const string& element_name, size_t count) {
  for (auto& elem : ply.elements) {
    if (elem.name == element_name) return true;
  }
  auto& elem = ply.elements.emplace_back();
  elem.name  = element_name;
  elem.count = count;
  return true;
}
inline bool add_property(ply_model& ply, const string& element_name,
    const string& property_name, size_t count, ply_type type, bool is_list) {
  if (!add_element(ply, element_name, count)) return false;
  for (auto& elem : ply.elements) {
    if (elem.name != element_name) continue;
    for (auto& prop : elem.properties) {
      if (prop.name == property_name) return true;
    }
    auto& prop   = elem.properties.emplace_back();
    prop.name    = property_name;
    prop.type    = type;
    prop.is_list = is_list;
    return true;
  }
  return false;
}
template <typename T>
inline bool set_value(ply_property& prop, size_t index, T value) {
  switch (prop.type) {
    case ply_type::i8: prop.data_i8[index] = (int8_t)value; return true;
    case ply_type::i16: prop.data_i16[index] = (int16_t)value; return true;
    case ply_type::i32: prop.data_i32[index] = (int32_t)value; return true;
    case ply_type::i64: prop.data_i64[index] = (int64_t)value; return true;
    case ply_type::u8: prop.data_u8[index] = (uint8_t)value; return true;
    case ply_type::u16: prop.data_u16[index] = (uint16_t)value; return true;
    case ply_type::u32: prop.data_u32[index] = (uint32_t)value; return true;
    case ply_type::u64: prop.data_u64[index] = (uint64_t)value; return true;
    case ply_type::f32: prop.data_f32[index] = (float)value; return true;
    case ply_type::f64: prop.data_f64[index] = (double)value; return true;
  }
  return false;
}
inline bool resize_values(ply_property& prop, size_t size) {
  switch (prop.type) {
    case ply_type::i8: prop.data_i8.resize(size); return true;
    case ply_type::i16: prop.data_i16.resize(size); return true;
    case ply_type::i32: prop.data_i32.resize(size); return true;
    case ply_type::i64: prop.data_i64.resize(size); return true;
    case ply_type::u8: prop.data_u8.resize(size); return true;
    case ply_type::u16: prop.data_u16.resize(size); return true;
    case ply_type::u32: prop.data_u32.resize(size); return true;
    case ply_type::u64: prop.data_u64.resize(size); return true;
    case ply_type::f32: prop.data_f32.resize(size); return true;
    case ply_type::f64: prop.data_f64.resize(size); return true;
  }
  return false;
}

template <typename T>
inline bool add_value(ply_model& ply, const string& element,
    const string& property, const vector<T>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, values.size(), get_ply_type<T>(), false))
    return false;
  auto& prop = get_property(ply, element, property);
  resize_values(prop, values.size());
  for (auto index = (size_t)0; index < values.size(); index++) {
    if (!set_value(prop, index, values[index])) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool add_values(ply_model& ply, const string& element,
    const array<string, N>& properties, const vector<array<T, N>>& values) {
  if (values.empty()) return false;
  for (auto& property : properties) {
    if (!add_property(
            ply, element, property, values.size(), get_ply_type<T>(), false))
      return false;
  }
  auto item = (size_t)0;
  for (auto& property : properties) {
    auto& prop = get_property(ply, element, property);
    resize_values(prop, values.size());
    for (auto index = (size_t)0; index < values.size(); index++) {
      if (!set_value(prop, index, values[index][item])) return false;
    }
    item++;
  }
  return true;
}

template <typename T>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<vector<T>>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, values.size(), get_ply_type<T>(), true))
    return false;
  auto& prop = get_property(ply, element, property);
  prop.ldata_u8.resize(values.size());
  auto count = (size_t)0;
  for (auto list = 0; list < values.size(); list++) {
    prop.ldata_u8[list] = values[list].size();
    count += values[list].size();
  }
  resize_values(prop, count);
  auto current = (size_t)0;
  for (auto list = (size_t)0; list < values.size(); list++) {
    for (auto item = (size_t)0; item < values[list].size(); item++) {
      if (!set_value(prop, current++, values[list][item])) return false;
    }
  }
  return true;
}
template <typename T>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<byte>& sizes,
    const vector<T>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, values.size(), get_ply_type<T>(), true))
    return false;
  auto& prop    = get_property(ply, element, property);
  prop.ldata_u8 = sizes;
  resize_values(prop, values.size());
  for (auto index = (size_t)0; index < values.size(); index++) {
    if (!set_value(prop, index, values[index])) return false;
  }
  return true;
}
template <typename T, size_t N>
inline bool add_lists(ply_model& ply, const string& element,
    const string& property, const vector<array<T, N>>& values) {
  if (values.empty()) return false;
  if (!add_property(
          ply, element, property, values.size(), get_ply_type<T>(), true))
    return false;
  auto& prop = get_property(ply, element, property);
  prop.ldata_u8.assign(values.size(), (uint8_t)N);
  resize_values(prop, values.size() * N);
  auto current = (size_t)0;
  for (auto list = (size_t)0; list < values.size(); list++) {
    for (auto item = (size_t)0; item < values[list].size(); item++) {
      if (!set_value(prop, current++, values[list][item])) return false;
    }
  }
  return true;
}

// Flip tex coords
template <typename T>
inline vector<array<T, 2>> flip_ply_texcoord(
    const vector<array<T, 2>>& texcoords) {
  auto flipped = texcoords;
  for (auto& uv : flipped) uv = {uv[0], 1 - uv[1]};
  return flipped;
}

// Add ply properties for meshes
template <typename T>
inline bool add_positions(ply_model& ply, const vector<array<T, 3>>& values) {
  return add_values(ply, "vertex", {"x", "y", "z"}, values);
}
template <typename T>
inline bool add_normals(ply_model& ply, const vector<array<T, 3>>& values) {
  return add_values(ply, "vertex", {"nx", "ny", "nz"}, values);
}
template <typename T>
inline bool add_texcoords(
    ply_model& ply, const vector<array<T, 2>>& values, bool flipv) {
  return add_values(
      ply, "vertex", {"u", "v"}, flipv ? flip_ply_texcoord(values) : values);
}
template <typename T>
inline bool add_colors(ply_model& ply, const vector<array<T, 3>>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue"}, values);
}
template <typename T>
inline bool add_colors(ply_model& ply, const vector<array<T, 4>>& values) {
  return add_values(ply, "vertex", {"red", "green", "blue", "alpha"}, values);
}
template <typename T>
inline bool add_radius(ply_model& ply, const vector<T>& values) {
  return add_value(ply, "vertex", "radius", values);
}
template <typename T>
inline bool add_faces(ply_model& ply, const vector<vector<T>>& values) {
  return add_lists(ply, "face", "vertex_indices", values);
}
template <typename T>
inline bool add_faces(ply_model& ply, const vector<array<T, 3>>& triangles,
    const vector<array<T, 4>>& quads) {
  if (triangles.empty() && quads.empty()) return false;
  if (quads.empty()) {
    return add_lists(ply, "face", "vertex_indices", triangles);
  } else if (triangles.empty() &&
             std::all_of(quads.begin(), quads.end(),
                 [](const array<T, 4>& q) { return q[2] != q[3]; })) {
    return add_lists(ply, "face", "vertex_indices", quads);
  } else {
    auto sizes   = vector<uint8_t>();
    auto indices = vector<int>{};
    sizes.reserve(triangles.size() + quads.size());
    indices.reserve(triangles.size() * 3 + quads.size() * 4);
    for (auto& t : triangles) {
      sizes.push_back(3);
      indices.push_back(t[0]);
      indices.push_back(t[1]);
      indices.push_back(t[2]);
    }
    for (auto& q : quads) {
      sizes.push_back(q[2] == q[3] ? 3 : 4);
      indices.push_back(q[0]);
      indices.push_back(q[1]);
      indices.push_back(q[2]);
      if (q[2] != q[3]) indices.push_back(q[3]);
    }
    return add_lists(ply, "face", "vertex_indices", sizes, indices);
  }
}
template <typename T>
inline bool add_triangles(ply_model& ply, const vector<array<T, 3>>& values) {
  return add_faces(ply, values, {});
}
template <typename T>
inline bool add_quads(ply_model& ply, const vector<array<T, 4>>& values) {
  return add_faces(ply, {}, values);
}
template <typename T>
inline bool add_lines(ply_model& ply, const vector<array<T, 2>>& values) {
  return add_lists(ply, "line", "vertex_indices", values);
}
template <typename T>
inline bool add_points(ply_model& ply, const vector<T>& values) {
  return add_lists(
      ply, "point", "vertex_indices", (const vector<array<T, 1>>&)values);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PLY LOADER AND WRITER --- DEPRECATED
// -----------------------------------------------------------------------------
namespace yocto {

[[deprecated]] inline bool get_values(const ply_model& ply,
    const string& element, const array<string, 2>& properties,
    vector<vec2f>& values) {
  return get_values(ply, element, properties, (vector<array<float, 2>>&)values);
}
[[deprecated]] inline bool get_values(const ply_model& ply,
    const string& element, const array<string, 3>& properties,
    vector<vec3f>& values) {
  return get_values(ply, element, properties, (vector<array<float, 3>>&)values);
}
[[deprecated]] inline bool get_values(const ply_model& ply,
    const string& element, const array<string, 4>& properties,
    vector<vec4f>& values) {
  return get_values(ply, element, properties, (vector<array<float, 4>>&)values);
}
[[deprecated]] inline bool get_values(const ply_model& ply,
    const string& element, const array<string, 12>& properties,
    vector<frame3f>& values) {
  return get_values(
      ply, element, properties, (vector<array<float, 12>>&)values);
}

// Get ply properties for meshes
[[deprecated]] inline bool get_positions(
    const ply_model& ply, vector<vec3f>& positions) {
  return get_positions(ply, (vector<array<float, 3>>&)positions);
}
[[deprecated]] inline bool get_normals(
    const ply_model& ply, vector<vec3f>& normals) {
  return get_normals(ply, (vector<array<float, 3>>&)normals);
}
[[deprecated]] inline bool get_texcoords(
    const ply_model& ply, vector<vec2f>& texcoords, bool flipv = false) {
  return get_texcoords(ply, (vector<array<float, 2>>&)texcoords, flipv);
}
[[deprecated]] inline bool get_colors(
    const ply_model& ply, vector<vec3f>& colors) {
  return get_colors(ply, (vector<array<float, 3>>&)colors);
}
[[deprecated]] inline bool get_colors(
    const ply_model& ply, vector<vec4f>& colors) {
  return get_colors(ply, (vector<array<float, 4>>&)colors);
}
[[deprecated]] inline bool get_triangles(
    const ply_model& ply, vector<vec3i>& triangles) {
  return get_triangles(ply, (vector<array<float, 3>>&)triangles);
}
[[deprecated]] inline bool get_quads(
    const ply_model& ply, vector<vec4i>& quads) {
  return get_quads(ply, (vector<array<int, 4>>&)quads);
}
[[deprecated]] inline bool get_faces(
    const ply_model& ply, vector<vec3i>& triangles, vector<vec4i>& quads) {
  return get_faces(
      ply, (vector<array<int, 3>>&)triangles, (vector<array<int, 4>>&)quads);
}
[[deprecated]] inline bool get_lines(
    const ply_model& ply, vector<vec2i>& lines) {
  return get_lines(ply, (vector<array<int, 2>>&)lines);
}

[[deprecated]] inline bool add_values(ply_model& ply, const string& element,
    const array<string, 2>& properties, const vector<vec2f>& values) {
  return add_values(
      ply, element, properties, (const vector<array<float, 2>>&)values);
}
[[deprecated]] inline bool add_values(ply_model& ply, const string& element,
    const array<string, 3>& properties, const vector<vec3f>& values) {
  return add_values(
      ply, element, properties, (const vector<array<float, 3>>&)values);
}
[[deprecated]] inline bool add_values(ply_model& ply, const string& element,
    const array<string, 4>& properties, const vector<vec4f>& values) {
  return add_values(
      ply, element, properties, (const vector<array<float, 4>>&)values);
}
[[deprecated]] inline bool add_values(ply_model& ply, const string& element,
    const array<string, 12>& properties, const vector<frame3f>& values) {
  return add_values(
      ply, element, properties, (const vector<array<float, 12>>&)values);
}

// Add ply properties for meshes
[[deprecated]] inline bool add_positions(
    ply_model& ply, const vector<vec3f>& values) {
  return add_positions(ply, (const vector<array<float, 3>>&)values);
}
[[deprecated]] inline bool add_normals(
    ply_model& ply, const vector<vec3f>& values) {
  return add_normals(ply, (const vector<array<float, 3>>&)values);
}
[[deprecated]] inline bool add_texcoords(
    ply_model& ply, const vector<vec2f>& values, bool flipv = false) {
  return add_texcoords(ply, (const vector<array<float, 2>>&)values, flipv);
}
[[deprecated]] inline bool add_colors(
    ply_model& ply, const vector<vec3f>& values) {
  return add_colors(ply, (const vector<array<float, 3>>&)values);
}
[[deprecated]] inline bool add_colors(
    ply_model& ply, const vector<vec4f>& values) {
  return add_colors(ply, (const vector<array<float, 4>>&)values);
}
[[deprecated]] inline bool add_faces(ply_model& ply,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
  return add_faces(ply, (const vector<array<int, 3>>&)triangles,
      (const vector<array<int, 4>>&)quads);
}
[[deprecated]] inline bool add_triangles(
    ply_model& ply, const vector<vec3i>& values) {
  return add_triangles(ply, (const vector<array<int, 3>>&)values);
}
[[deprecated]] inline bool add_quads(
    ply_model& ply, const vector<vec4i>& values) {
  return add_quads(ply, (const vector<array<int, 4>>&)values);
}
[[deprecated]] inline bool add_lines(
    ply_model& ply, const vector<vec2i>& values) {
  return add_lines(ply, (const vector<array<int, 2>>&)values);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// OBJ LOADER AND WRITER --- DEPRECATED
// -----------------------------------------------------------------------------
namespace yocto {

// Get obj shape.
[[deprecated]] inline void get_positions(
    const obj_shape& obj, vector<vec3f>& positions) {
  return get_positions(obj, (vector<array<float, 3>>&)positions);
}
[[deprecated]] inline void get_normals(
    const obj_shape& obj, vector<vec3f>& normals) {
  return get_normals(obj, (vector<array<float, 3>>&)normals);
}
[[deprecated]] inline void get_texcoords(
    const obj_shape& obj, vector<vec2f>& texcoords, bool flipv = false) {
  return get_texcoords(obj, (vector<array<float, 2>>&)texcoords);
}
[[deprecated]] inline void get_faces(const obj_shape& obj,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<int>& materials) {
  return get_faces(obj, (vector<array<int, 3>>&)triangles,
      (vector<array<int, 4>>&)quads, materials);
}
[[deprecated]] inline void get_triangles(
    const obj_shape& obj, vector<vec3i>& triangles, vector<int>& materials) {
  return get_triangles(obj, (vector<array<int, 3>>&)triangles, materials);
}
[[deprecated]] inline void get_quads(
    const obj_shape& obj, vector<vec4i>& quads, vector<int>& materials) {
  return get_quads(obj, (vector<array<int, 4>>&)quads, materials);
}
[[deprecated]] inline void get_lines(
    const obj_shape& obj, vector<vec2i>& lines, vector<int>& materials) {
  return get_lines(obj, (vector<array<int, 2>>&)lines, materials);
}
[[deprecated]] inline void get_fvquads(const obj_shape& obj,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<int>& materials) {
  return get_fvquads(obj, (vector<array<int, 4>>&)quadspos,
      (vector<array<int, 4>>&)quadsnorm, (vector<array<int, 4>>&)quadstexcoord,
      materials);
}
[[deprecated]] inline void get_faces(const obj_shape& obj, int material,
    vector<vec3i>& triangles, vector<vec4i>& quads) {
  return get_faces(obj, material, (vector<array<int, 3>>&)triangles,
      (vector<array<int, 4>>&)quads);
}
[[deprecated]] inline void get_triangles(
    const obj_shape& obj, int material, vector<vec3i>& triangles) {
  return get_triangles(obj, material, (vector<array<int, 3>>&)triangles);
}
[[deprecated]] inline void get_quads(
    const obj_shape& obj, int material, vector<vec4i>& quads) {
  return get_quads(obj, material, (vector<array<int, 4>>&)quads);
}
[[deprecated]] inline void get_lines(
    const obj_shape& obj, int material, vector<vec2i>& lines) {
  return get_lines(obj, material, (vector<array<int, 2>>&)lines);
}

// Add obj shape
[[deprecated]] inline void add_positions(
    obj_shape& obj, const vector<vec3f>& positions) {
  return add_positions(obj, (vector<array<float, 3>>&)positions);
}
[[deprecated]] inline void add_normals(
    obj_shape& obj, const vector<vec3f>& normals) {
  return add_normals(obj, (vector<array<float, 3>>&)normals);
}
[[deprecated]] inline void add_texcoords(
    obj_shape& obj, const vector<vec2f>& texcoords, bool flipv = false) {
  return add_texcoords(obj, (vector<array<float, 2>>&)texcoords, flipv);
}
[[deprecated]] inline void add_triangles(obj_shape& obj,
    const vector<vec3i>& triangles, int material, bool has_normals,
    bool has_texcoord) {
  return add_triangles(obj, (vector<array<int, 3>>&)triangles, material,
      has_normals, has_texcoord);
}
[[deprecated]] inline void add_quads(obj_shape& obj, const vector<vec4i>& quads,
    int material, bool has_normals, bool has_texcoord) {
  return add_quads(
      obj, (vector<array<int, 4>>&)quads, material, has_normals, has_texcoord);
}
[[deprecated]] inline void add_lines(obj_shape& obj, const vector<vec2i>& lines,
    int material, bool has_normals, bool has_texcoord) {
  return add_lines(
      obj, (vector<array<int, 2>>&)lines, material, has_normals, has_texcoord);
}
[[deprecated]] inline void add_fvquads(obj_shape& obj,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, int material) {
  return add_fvquads(obj, (vector<array<int, 4>>&)quadspos,
      (vector<array<int, 4>>&)quadsnorm, (vector<array<int, 4>>&)quadstexcoord,
      material);
}

}  // namespace yocto

#endif
