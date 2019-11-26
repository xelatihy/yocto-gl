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

#include <algorithm>

#include "yocto_math.h"

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
vector<vec2f>   get_ply_values(const ply_model& ply, const string& element,
      const string& property1, const string& property2);
vector<vec3f>   get_ply_values(const ply_model& ply, const string& element,
      const string& property1, const string& property2, const string& property3);
vector<vec4f>   get_ply_values(const ply_model& ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      const string& property4);
vector<vec4f>   get_ply_values(const ply_model& ply, const string& element,
      const string& property1, const string& property2, const string& property3,
      float property4);
vector<frame3f> get_ply_values(const ply_model& ply, const string& element,
    const array<string, 12>& properties);

vector<vector<int>> get_ply_lists(
    const ply_model& ply, const string& element, const string& property);
vector<byte> get_ply_list_sizes(
    const ply_model& ply, const string& element, const string& property);
vector<int> get_ply_list_values(
    const ply_model& ply, const string& element, const string& property);
vec2i get_ply_list_minxmax(
    const ply_model& ply, const string& element, const string& property);

// Get ply properties for meshes
vector<vec3f>       get_ply_positions(const ply_model& ply);
vector<vec3f>       get_ply_normals(const ply_model& ply);
vector<vec2f>       get_ply_texcoords(const ply_model& ply, bool flipv = false);
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
void add_ply_values(ply_model& ply, const vector<frame3f>& values,
    const string& element, const array<string, 12>& properties);

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

// Read Ply functions
plyio_status read_ply_header(const string& filename, FILE* fs,
    ply_format& format, vector<ply_element>& elements,
    vector<string>& comments);
plyio_status read_ply_value(const string& filename, FILE* fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
plyio_status read_ply_value(const string& filename, FILE* fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Write Ply functions
plyio_status write_ply_header(const string& filename, FILE* fs,
    ply_format format, const vector<ply_element>& elements,
    const vector<string>& comments);
plyio_status write_ply_value(const string& filename, FILE* fs,
    ply_format format, const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
plyio_status write_ply_value(const string& filename, FILE* fs,
    ply_format format, const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Helpers to get element and property indices
int   find_ply_element(const vector<ply_element>& elements, const string& name);
int   find_ply_property(const ply_element& element, const string& name);
vec2i find_ply_property(
    const ply_element& element, const string& name1, const string& name2);
vec3i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3);
vec4i find_ply_property(const ply_element& element, const string& name1,
    const string& name2, const string& name3, const string& name4);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// OBJ vertex
struct obj_vertex {
  int position = 0;
  int texcoord = 0;
  int normal   = 0;
};

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
  return a.position == b.position && a.texcoord == b.texcoord &&
         a.normal == b.normal;
}

// Obj texture information.
struct obj_texture_info {
  string path  = "";     // file path
  bool   clamp = false;  // clamp to edge
  float  scale = 1;      // scale for bump/displacement

  // Properties not explicitly handled.
  unordered_map<string, vector<float>> props;

  obj_texture_info() {}
  obj_texture_info(const char* path) : path{path} {}
  obj_texture_info(const string& path) : path{path} {}
};

// Obj element
struct obj_element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj shape
struct obj_shape {
  string              name      = "";
  vector<vec3f>       positions = {};
  vector<vec3f>       normals   = {};
  vector<vec2f>       texcoords = {};
  vector<string>      materials = {};
  vector<obj_vertex>  vertices  = {};
  vector<obj_element> faces     = {};
  vector<obj_element> lines     = {};
  vector<obj_element> points    = {};
  vector<frame3f>     instances = {};
};

// Obj material
struct obj_material {
  // material name and type
  string name  = "";
  int    illum = 0;

  // material colors and values
  vec3f emission     = zero3f;
  vec3f ambient      = zero3f;
  vec3f diffuse      = zero3f;
  vec3f specular     = zero3f;
  vec3f reflection   = zero3f;
  vec3f transmission = zero3f;
  float exponent     = 10;
  float ior          = 1.5;
  float opacity      = 1;

  // material textures
  obj_texture_info emission_map     = {};
  obj_texture_info ambient_map      = {};
  obj_texture_info diffuse_map      = {};
  obj_texture_info specular_map     = {};
  obj_texture_info reflection_map   = {};
  obj_texture_info transmission_map = {};
  obj_texture_info exponent_map     = {};
  obj_texture_info opacity_map      = {};
  obj_texture_info bump_map         = {};
  obj_texture_info normal_map       = {};
  obj_texture_info displacement_map = {};

  // pbrt extension values
  float pbr_roughness     = 0;
  float pbr_metallic      = 0;
  float pbr_sheen         = 0;
  float pbr_clearcoat     = 0;
  float pbr_coatroughness = 0;

  // pbr extension textures
  obj_texture_info pbr_roughness_map     = {};
  obj_texture_info pbr_metallic_map      = {};
  obj_texture_info pbr_sheen_map         = {};
  obj_texture_info pbr_clearcoat_map     = {};
  obj_texture_info pbr_coatroughness_map = {};

  // volume extension colors and values
  vec3f vol_emission     = zero3f;
  vec3f vol_transmission = zero3f;
  vec3f vol_meanfreepath = zero3f;
  vec3f vol_scattering   = zero3f;
  float vol_anisotropy   = 0;
  float vol_scale        = 0.01;

  // volument textures
  obj_texture_info vol_scattering_map = {};
};

// Obj camera
struct obj_camera {
  string  name     = "";
  frame3f frame    = identity3x4f;
  bool    ortho    = false;
  float   width    = 0.036;
  float   height   = 0.028;
  float   lens     = 0.50;
  float   focus    = 0;
  float   aperture = 0;
};

// Obj environment
struct obj_environment {
  string           name         = "";
  frame3f          frame        = identity3x4f;
  vec3f            emission     = zero3f;
  obj_texture_info emission_map = {};
};

// Obj model
struct obj_model {
  vector<string>          comments     = {};
  vector<obj_shape>       shapes       = {};
  vector<obj_material>    materials    = {};
  vector<obj_camera>      cameras      = {};
  vector<obj_environment> environments = {};
};

// Result of io operations
struct objio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

// Load and save obj
objio_status load_obj(const string& filename, obj_model& obj,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
objio_status save_obj(const string& filename, const obj_model& obj);

// convert between roughness and exponent
float obj_exponent_to_roughness(float exponent);
float obj_roughness_to_exponent(float roughness);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_obj_triangles(const obj_model& obj, const obj_shape& shape,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
void get_obj_quads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
void get_obj_lines(const obj_model& obj, const obj_shape& shape,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
void get_obj_points(const obj_model& obj, const obj_shape& shape,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
void get_obj_fvquads(const obj_model& obj, const obj_shape& shape,
    vector<vec4i>& quadspos, vector<vec4i>& quadsnorm,
    vector<vec4i>& quadstexcoord, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texcoords, vector<string>& materials,
    vector<int>& ematerials, bool flip_texcoord = false);
bool has_obj_quads(const obj_shape& shape);

// Add obj shape
void add_obj_triangles(obj_model& obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void add_obj_quads(obj_model& obj, const string& name,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void add_obj_lines(obj_model& obj, const string& name,
    const vector<vec2i>& lines, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void add_obj_points(obj_model& obj, const string& name,
    const vector<int>& points, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void add_obj_fvquads(obj_model& obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<string>& materials = {}, const vector<int>& ematerials = {},
    bool flip_texcoord = false);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------
namespace yocto {

// Obj/Mtl/Objx command
enum struct obj_command {
  // clang-format off
  vertex, normal, texcoord,         // data in value
  face, str, point,                 // data in vertices
  object, group, usemtl, smoothing, // data in name
  mtllib, objxlib,                  // data in name
  error                             // no data
  // clang-format on
};
enum struct mtl_command { material, error };
enum struct objx_command { camera, environment, instance, error };

// Obj instance
struct obj_instance {
  string  object = "";
  frame3f frame  = identity3x4f;
};

// Read obj/mtl/objx elements
objio_status read_obj_command(const string& filename, FILE* fs,
    obj_command& command, string& name, vec3f& value,
    vector<obj_vertex>& vertices, obj_vertex& vert_size);
objio_status read_mtl_command(const string& filename, FILE* fs,
    mtl_command& command, obj_material& material, bool fliptr = true);
objio_status read_objx_command(const string& filename, FILE* fs,
    objx_command& command, obj_camera& camera, obj_environment& environment,
    obj_instance& instance);

// Write obj/mtl/objx elements
objio_status write_obj_comment(
    const string& filename, FILE* fs, const string& comment);
objio_status write_obj_command(const string& filename, FILE* fs,
    obj_command command, const string& name, const vec3f& value,
    const vector<obj_vertex>& vertices = {});
objio_status write_mtl_command(const string& filename, FILE* fs,
    mtl_command command, obj_material& material,
    const obj_texture_info& texture = {});
objio_status write_objx_command(const string& filename, FILE* fs,
    objx_command command, const obj_camera& camera,
    const obj_environment& environment, const obj_instance& instance);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::obj_vertex> {
  size_t operator()(const yocto::obj_vertex& v) const {
    static const std::hash<int> hasher = std::hash<int>();
    auto                        h      = (size_t)0;
    h ^= hasher(v.position) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.normal) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.texcoord) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// LOW-LEVEL YAML DECLARATIONS
// -----------------------------------------------------------------------------
namespace yocto {

using std::string_view;

// Yaml value type
enum struct yaml_value_type { number, boolean, string, array };

// Yaml value
struct yaml_value {
  yaml_value_type   type    = yaml_value_type::number;
  double            number  = 0;
  bool              boolean = false;
  string            string_ = "";
  array<double, 16> array_  = {};
};

// Yaml element
struct yaml_element {
  string                           name       = "";
  vector<pair<string, yaml_value>> key_values = {};
};

// Yaml model
struct yaml_model {
  vector<string>       comments = {};
  vector<yaml_element> elements = {};
};

// // Result of io operations
struct yamlio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

// Load/save yaml
yamlio_status load_yaml(const string& filename, yaml_model& yaml);
yamlio_status save_yaml(const string& filename, const yaml_model& yaml);

// Load Yaml properties
yamlio_status read_yaml_property(const string& filename, FILE* fs,
    string& group, string& key, bool& newobj, bool& done, yaml_value& value);

// Write Yaml properties
yamlio_status write_yaml_comment(
    const string& filename, FILE* fs, const string& comment);
yamlio_status write_yaml_property(const string& filename, FILE* fs,
    const string& object, const string& key, bool newobj,
    const yaml_value& value);
yamlio_status write_yaml_object(
    const string& filename, FILE* fs, const string& object);

// type-cheked yaml value access
bool get_yaml_value(const yaml_value& yaml, string& value);
bool get_yaml_value(const yaml_value& yaml, bool& value);
bool get_yaml_value(const yaml_value& yaml, int& value);
bool get_yaml_value(const yaml_value& yaml, float& value);
bool get_yaml_value(const yaml_value& yaml, vec2f& value);
bool get_yaml_value(const yaml_value& yaml, vec3f& value);
bool get_yaml_value(const yaml_value& yaml, mat3f& value);
bool get_yaml_value(const yaml_value& yaml, frame3f& value);
template <typename T>
inline bool get_yaml_value(
    const yaml_element& element, const string& name, const T& value);
bool has_yaml_value(const yaml_element& element, const string& name);

// yaml value construction
yaml_value make_yaml_value(const string& value);
yaml_value make_yaml_value(bool value);
yaml_value make_yaml_value(int value);
yaml_value make_yaml_value(float value);
yaml_value make_yaml_value(const vec2f& value);
yaml_value make_yaml_value(const vec3f& value);
yaml_value make_yaml_value(const mat3f& value);
yaml_value make_yaml_value(const frame3f& value);
template <typename T>
inline bool add_yaml_value(
    yaml_element& element, const string& name, const T& value);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt value type
enum struct pbrt_value_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color, 
  point2, vector2, spectrum
  // clang-format on
};

// Pbrt value
struct pbrt_value {
  string          name     = "";
  pbrt_value_type type     = pbrt_value_type::real;
  int             value1i  = 0;
  float           value1f  = 0;
  vec2f           value2f  = {0, 0};
  vec3f           value3f  = {0, 0, 0};
  bool            value1b  = false;
  string          value1s  = "";
  vector<float>   vector1f = {};
  vector<vec2f>   vector2f = {};
  vector<vec3f>   vector3f = {};
  vector<int>     vector1i = {};
};

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // camera approximation
  float width    = 0;
  float height   = 0;
  float lens     = 0;
  float aspect   = 0;
  float focus    = 0;
  float aperture = 0;
};

// Pbrt texture
struct pbrt_texture {
  // texture parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  // texture approximation
  vec3f  constant = vec3f{1, 1, 1};
  string filename = "";
};

// Pbrt material
struct pbrt_material {
  // material parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  // material approximation
  vec3f  diffuse          = zero3f;
  vec3f  specular         = zero3f;
  vec3f  transmission     = zero3f;
  vec2f  roughness        = zero2f;
  vec3f  opacity          = vec3f{1};
  vec3f  eta              = zero3f;
  vec3f  etak             = zero3f;
  vec3f  sspecular        = zero3f;  // specular scaled by fresnel
  string diffuse_map      = "";
  string specular_map     = "";
  string transmission_map = "";
  string roughness_map    = "";
  string opacity_map      = "";
  string eta_map          = "";
  string etak_map         = "";
  vec3f  volmeanfreepath  = vec3f{0};
  vec3f  volscatter       = vec3f{0};
  float  volscale         = 0.01;
  bool   refract          = false;
};

// Pbrt medium
struct pbrt_medium {
  // medium parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
};

// Pbrt shape
struct pbrt_shape {
  // shape parameters
  string             type            = "";
  vector<pbrt_value> values          = {};
  frame3f            frame           = identity3x4f;
  frame3f            frend           = identity3x4f;
  string             material        = "";
  string             arealight       = "";
  string             interior        = "";
  string             exterior        = "";
  bool               is_instanced    = false;
  vector<frame3f>    instance_frames = {};
  vector<frame3f>    instance_frends = {};
  // shape approximation
  string        filename  = "";
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec3i> triangles = {};
  float         radius    = 0;  // radius for sphere, cylinder, disk
};

// Pbrt lights
struct pbrt_light {
  // light parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // light approximation
  vec3f emission = zero3f;
  vec3f from     = zero3f;
  vec3f to       = zero3f;
  bool  distant  = false;
  // arealight approximation
  vec3f         area_emission  = zero3f;
  frame3f       area_frame     = identity3x4f;
  frame3f       area_frend     = identity3x4f;
  vector<vec3i> area_triangles = {};
  vector<vec3f> area_positions = {};
  vector<vec3f> area_normals   = {};
};
struct pbrt_arealight {
  // arealight parameters
  string             name   = "";
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // arealight approximation
  vec3f emission = zero3f;
};
struct pbrt_environment {
  // shape parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  frame3f            frame  = identity3x4f;
  frame3f            frend  = identity3x4f;
  // environment approximation
  vec3f  emission = zero3f;
  string filename = "";
};

// Other pbrt elements
struct pbrt_integrator {
  // integrator parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_film {
  // film parameters
  string             type   = "";
  vector<pbrt_value> values = {};
  // film approximation
  string filename   = "";
  vec2i  resolution = zero2i;
};
struct pbrt_filter {
  // filter parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_accelerator {
  // accelerator parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};
struct pbrt_sampler {
  // sampler parameters
  string             type   = "";
  vector<pbrt_value> values = {};
};

// Pbrt model
struct pbrt_model {
  vector<string>           comments     = {};
  vector<pbrt_camera>      cameras      = {};
  vector<pbrt_shape>       shapes       = {};
  vector<pbrt_texture>     textures     = {};
  vector<pbrt_material>    materials    = {};
  vector<pbrt_medium>      mediums      = {};
  vector<pbrt_environment> environments = {};
  vector<pbrt_arealight>   arealights   = {};
  vector<pbrt_light>       lights       = {};
  vector<pbrt_integrator>  integrators  = {};
  vector<pbrt_film>        films        = {};
  vector<pbrt_filter>      filters      = {};
  vector<pbrt_sampler>     samplers     = {};
  vector<pbrt_accelerator> accelerators = {};
};

// Result of io operations
struct pbrtio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

// Load/save pbrt
pbrtio_status load_pbrt(const string& filename, pbrt_model& pbrt);
pbrtio_status save_pbrt(const string& filename, const pbrt_model& pbrt);

}  // namespace yocto

// -----------------------------------------------------------------------------
// LOW-LEVEL INTERFACE
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt command
enum struct pbrt_command {
  // clang-format off
  world_begin, world_end, attribute_begin, attribute_end,
  transform_begin, transform_end, reverse_orientation,
  set_transform, concat_transform, lookat_transform,
  object_instance, object_begin, object_end, include,      
  sampler, integrator, accelerator, film, filter, camera, shape, light,
  material, arealight, named_texture, named_medium, named_material,             
  use_material, medium_interface, active_transform,
  coordinate_system_set, coordinate_system_transform,
  error
  // clang-format on
};

// Read pbrt commands
pbrtio_status read_pbrt_command(const string& filename, FILE* fs,
    pbrt_command& command, string& name, string& type, frame3f& xform,
    vector<pbrt_value>& values);
pbrtio_status read_pbrt_command(const string& filename, FILE* fs,
    pbrt_command& command, string& name, string& type, frame3f& xform,
    vector<pbrt_value>& values, string& buffer);

// Write pbrt commands
pbrtio_status write_pbrt_comment(
    const string& filename, FILE* fs, const string& comment);
pbrtio_status write_pbrt_command(const string& filename, FILE* fs,
    pbrt_command command, const string& name, const string& type,
    const frame3f& xform, const vector<pbrt_value>& values,
    bool texture_as_float = false);
pbrtio_status write_pbrt_command(const string& filename, FILE* fs,
    pbrt_command command, const string& name = "",
    const frame3f& xform = identity3x4f);
pbrtio_status write_pbrt_command(const string& filename, FILE* fs,
    pbrt_command command, const string& name, const string& type,
    const vector<pbrt_value>& values, bool texture_as_float = false);

// type-cheked pbrt value access
bool get_pbrt_value(const pbrt_value& pbrt, string& value);
bool get_pbrt_value(const pbrt_value& pbrt, bool& value);
bool get_pbrt_value(const pbrt_value& pbrt, int& value);
bool get_pbrt_value(const pbrt_value& pbrt, float& value);
bool get_pbrt_value(const pbrt_value& pbrt, vec2f& value);
bool get_pbrt_value(const pbrt_value& pbrt, vec3f& value);
bool get_pbrt_value(const pbrt_value& pbrt, vector<float>& value);
bool get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& value);
bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& value);
bool get_pbrt_value(const pbrt_value& pbrt, vector<int>& value);
bool get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value);
bool get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value);
bool get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value);
template <typename T>
inline bool get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& value);

// pbrt value construction
pbrt_value make_pbrt_value(const string& name, const string& value,
    pbrt_value_type type = pbrt_value_type::string);
pbrt_value make_pbrt_value(const string& name, bool value,
    pbrt_value_type type = pbrt_value_type::boolean);
pbrt_value make_pbrt_value(const string& name, int value,
    pbrt_value_type type = pbrt_value_type::integer);
pbrt_value make_pbrt_value(const string& name, float value,
    pbrt_value_type type = pbrt_value_type::real);
pbrt_value make_pbrt_value(const string& name, const vec2f& value,
    pbrt_value_type type = pbrt_value_type::point2);
pbrt_value make_pbrt_value(const string& name, const vec3f& value,
    pbrt_value_type type = pbrt_value_type::color);
pbrt_value make_pbrt_value(const string& name, const vector<vec2f>& value,
    pbrt_value_type type = pbrt_value_type::point2);
pbrt_value make_pbrt_value(const string& name, const vector<vec3f>& value,
    pbrt_value_type type = pbrt_value_type::point);
pbrt_value make_pbrt_value(const string& name, const vector<vec3i>& value,
    pbrt_value_type type = pbrt_value_type::integer);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE GLTF LOADER DECLARATIONS
// -----------------------------------------------------------------------------
namespace yocto {

struct gltf_camera {
  string name   = "";
  bool   ortho  = false;
  float  yfov   = 45 * pif / 180;
  float  aspect = 1;
};
struct gltf_texture {
  string name     = "";
  string filename = "";
};
struct gltf_material {
  string name            = "";
  vec3f  emission        = zero3f;
  int    emission_tex    = -1;
  int    normal_tex      = -1;
  bool   has_metalrough  = false;
  vec4f  mr_base         = zero4f;
  float  mr_metallic     = 0;
  float  mr_roughness    = 1;
  int    mr_base_tex     = -1;
  int    mr_metallic_tex = -1;
  bool   has_specgloss   = false;
  vec4f  sg_diffuse      = zero4f;
  vec3f  sg_specular     = zero3f;
  float  sg_glossiness   = 1;
  int    sg_diffuse_tex  = -1;
  int    sg_specular_tex = -1;
};
struct gltf_primitive {
  int           material  = -1;
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec4f> colors    = {};
  vector<float> radius    = {};
  vector<vec4f> tangents  = {};
  vector<vec3i> triangles = {};
  vector<vec2i> lines     = {};
  vector<int>   points    = {};
};
struct gltf_mesh {
  string                 name       = "";
  vector<gltf_primitive> primitives = {};
};
struct gltf_node {
  string      name        = "";
  frame3f     frame       = {};
  vec3f       translation = zero3f;
  vec4f       rotation    = vec4f{0, 0, 0, 1};
  vec3f       scale       = vec3f{1};
  frame3f     local       = identity3x4f;
  int         camera      = -1;
  int         mesh        = -1;
  int         parent      = -1;
  vector<int> children    = {};
};
struct gltf_scene {
  string      name  = "";
  vector<int> nodes = {};
};
struct gltf_model {
  vector<gltf_camera>   cameras   = {};
  vector<gltf_mesh>     meshes    = {};
  vector<gltf_texture>  textures  = {};
  vector<gltf_material> materials = {};
  vector<gltf_node>     nodes     = {};
  vector<gltf_scene>    scenes    = {};
};

// Result of io operations
struct gltfio_status {
  string   error = {};
  explicit operator bool() const { return error.empty(); }
};

gltfio_status load_gltf(const string& filename, gltf_model& gltf);

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline bool get_yaml_value(
    const yaml_element& element, const string& name, T& value) {
  for (auto& [key, value_] : element.key_values) {
    if (key == name) return get_yaml_value(value_, value);
  }
  return true;
}

template <typename T>
inline bool add_yaml_value(
    yaml_element& element, const string& name, const T& value) {
  for (auto& [key, value] : element.key_values)
    if (key == name) return false;
  element.key_values.push_back({name, make_yaml_value(value)});
  return true;
}

template <typename T>
inline bool get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& value) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      return get_pbrt_value(p, value);
    }
  }
  return true;
}

}  // namespace yocto

#endif
