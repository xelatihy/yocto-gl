//
// # Yocto/ModelIO: L<ow-level library for PLY/OBJ/Pbrt/Yaml parsing and writing
//
// Yocto/ModelIO is a collection of simple parsers for Stanford Ply,
// Wavefront Obj, Pbrt, and Yaml formats. The prasers are designed for large
// files and do keep a copy of the model in memory.
// Yocto/ModelIO provides fast/low-level access to model data and requires some
// familiarity with the formats to use effectively. For a higher level
// interface, consider using Yocto/Shape's `load_shape()` and `save_shape()`,
// or Yocto/SceneIO's `load_scene()` and `save_scene()`
//
// Yocto/ModelIO also support writing Ply/Obj/Yaml files again without keeping
// a copy of the model but instead writing elements directly after each call.
// Error reporting is done by throwing `std::runtime_error` exceptions.
//
//
//
// ## Load PLY
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
// ## Write PLY
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
// FILE AND PROPERTY HANDLING
// -----------------------------------------------------------------------------
namespace yocto {

// A class that wraps a C file ti handle safe opening/closgin with RIIA.
struct file_wrapper {
  file_wrapper() {}
  file_wrapper(file_wrapper&& other);
  file_wrapper(const file_wrapper&) = delete;
  file_wrapper& operator=(const file_wrapper&) = delete;
  ~file_wrapper();

  FILE*  fs       = nullptr;
  string filename = "";
  string mode     = "rt";
  int    linenum  = 0;
};

// open a file
file_wrapper open_file(const string& filename, const string& mode = "rt");
void         open_file(
            file_wrapper& fs, const string& filename, const string& mode = "rt");
void close_file(file_wrapper& fs);

}  // namespace yocto

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
  string   name       = "";
  bool     is_list    = false;
  ply_type value_type = ply_type::f32;
  ply_type list_type  = ply_type::f32;
};

// Ply elements
struct ply_element {
  string               name       = "";
  size_t               count      = 0;
  vector<ply_property> properties = {};
};

// Read Ply functions
void read_ply_header(file_wrapper& fs, ply_format& format,
    vector<ply_element>& elements, vector<string>& comments);
void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
void read_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
    vector<vector<int>>& lists);

// Write Ply functions
void write_ply_header(file_wrapper& fs, ply_format format,
    const vector<ply_element>& elements, const vector<string>& comments);
void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<double>& values,
    vector<vector<double>>& lists);
void write_ply_value(file_wrapper& fs, ply_format format,
    const ply_element& element, vector<float>& values,
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

// Obj command
enum struct obj_command {
  // clang-format off
  vertex, normal, texcoord,         // data in value
  face, line, point,                // data in vertices
  object, group, usemtl, smoothing, // data in name
  mtllib, objxlib,                  // data in name
  // clang-format on
};

// Mtl command
enum struct mtl_command {
  // clang-format off
  // material name and type (value)
  material, illum,
  // material colors
  emission, ambient, diffuse, specular, reflection, transmission,
  // material values
  exponent, ior, opacity,
  // material textures
  emission_map, ambient_map, diffuse_map, specular_map, reflection_map,  
  transmission_map, exponent_map, opacity_map, bump_map, normal_map, 
  displacement_map,                  
  // pbrt extension values
  pbr_roughness, pbr_metallic, pbr_sheen, pbr_clearcoat, pbr_coatroughness,
  // pbr extension textures
  pbr_roughness_map, pbr_metallic_map, pbr_sheen_map,
  pbr_clearcoat_map, pbr_coatroughness_map,
  // volume extension colors
  vol_transmission, vol_meanfreepath, vol_scattering, vol_emission,
  // volume extension values
  vol_anisotropy, vol_scale,
  // volument textures
  vol_scattering_map
  // clang-format on
};

// Objx command
enum struct objx_command {
  // clang-format off
  // object names
  camera, environment, instance, procedural,
  // object frames
  frame,
  // camera values
  ortho, width, height, lens, aperture, focus,
  // environment values
  emission, emission_map,
  // instance/procedural values
  object, material
  // clang-format on
};

// Obj value type
enum struct obj_value_type { number, boolean, string, array };

// Obj value
struct obj_value {
  obj_value_type    type    = obj_value_type::number;
  double            number  = 0;
  bool              boolean = false;
  string            string_ = "";
  array<double, 16> array_  = {};
};

// Read obj elements
bool read_obj_command(file_wrapper& fs, obj_command& command, obj_value& value,
    vector<obj_vertex>& vertices, obj_vertex& vert_size);
bool read_mtl_command(file_wrapper& fs, mtl_command& command, obj_value& value,
    obj_texture_info& texture, bool fliptr = true);
bool read_objx_command(file_wrapper& fs, objx_command& command,
    obj_value& value, obj_texture_info& texture);

// Write obj elements
void write_obj_comment(file_wrapper& fs, const string& comment);
void write_obj_command(file_wrapper& fs, obj_command command,
    const obj_value& value, const vector<obj_vertex>& vertices = {});
void write_mtl_command(file_wrapper& fs, mtl_command command,
    const obj_value& value, const obj_texture_info& texture = {});
void write_objx_command(file_wrapper& fs, objx_command command,
    const obj_value& value, const obj_texture_info& texture = {});

// typesafe access of obj value
void get_obj_value(const obj_value& yaml, string& value);
void get_obj_value(const obj_value& yaml, bool& value);
void get_obj_value(const obj_value& yaml, int& value);
void get_obj_value(const obj_value& yaml, float& value);
void get_obj_value(const obj_value& yaml, vec2f& value);
void get_obj_value(const obj_value& yaml, vec3f& value);
void get_obj_value(const obj_value& yaml, mat3f& value);
void get_obj_value(const obj_value& yaml, frame3f& value);

// typesafe access of obj value
obj_value make_obj_value(const string& value);
obj_value make_obj_value(bool value);
obj_value make_obj_value(int value);
obj_value make_obj_value(float value);
obj_value make_obj_value(const vec2f& value);
obj_value make_obj_value(const vec3f& value);
obj_value make_obj_value(const mat3f& value);
obj_value make_obj_value(const frame3f& value);

}  // namespace yocto

// -----------------------------------------------------------------------------
// OLD OBJ INTERFACE
// -----------------------------------------------------------------------------
namespace yocto {

// Obj material.
struct obj_material {
  string name  = "";  // name
  int    illum = 0;   // MTL illum mode

  // base values
  vec3f ke  = {0, 0, 0};  // emission color
  vec3f ka  = {0, 0, 0};  // ambient color
  vec3f kd  = {0, 0, 0};  // diffuse color
  vec3f ks  = {0, 0, 0};  // specular color
  vec3f kr  = {0, 0, 0};  // reflection color
  vec3f kt  = {0, 0, 0};  // transmission color
  float ns  = 0;          // Phong exponent color
  float ior = 1;          // index of refraction
  float op  = 1;          // opacity

  // textures
  obj_texture_info ke_map   = "";  // emission texture
  obj_texture_info ka_map   = "";  // ambient texture
  obj_texture_info kd_map   = "";  // diffuse texture
  obj_texture_info ks_map   = "";  // specular texture
  obj_texture_info kr_map   = "";  // reflection texture
  obj_texture_info kt_map   = "";  // transmission texture
  obj_texture_info ns_map   = "";  // Phong exponent texture
  obj_texture_info op_map   = "";  // opacity texture
  obj_texture_info ior_map  = "";  // ior texture
  obj_texture_info bump_map = "";  // bump map
  obj_texture_info norm_map = "";  // normal map
  obj_texture_info disp_map = "";  // displacement map
  obj_texture_info occ_map  = "";  // occlusion map

  // pbr values
  float pr  = 0;  // roughness
  float pm  = 0;  // metallic
  float ps  = 0;  // sheen
  float pc  = 0;  // coat
  float pcr = 0;  // coat roughness

  // textures
  obj_texture_info pr_map  = "";  // roughness texture
  obj_texture_info pm_map  = "";  // metallic texture
  obj_texture_info ps_map  = "";  // sheen texture
  obj_texture_info pc_map  = "";  // coat texture
  obj_texture_info pcr_map = "";  // coat roughness texture

  // volume values
  vec3f vt = {0, 0, 0};  // volumetric transmission
  vec3f vp = {0, 0, 0};  // volumetric mean-free-path
  vec3f ve = {0, 0, 0};  // volumetric emission
  vec3f vs = {0, 0, 0};  // volumetric scattering
  float vg = 0;          // volumetric anisotropy (phase g)
  float vr = 0.01;       // volumetric scale

  // textures
  obj_texture_info vs_map = "";  // scattering texture

  // Properties not explicitly handled.
  unordered_map<string, vector<string>> props;
};

// Obj camera [extension].
struct obj_camera {
  string  name     = "";            // name
  frame3f frame    = identity3x4f;  // transform
  bool    ortho    = false;         // orthographic
  float   width    = 0.036f;        // film size (default to 35mm)
  float   height   = 0.024f;        // film size (default to 35mm)
  float   lens     = 0.050f;        // focal length
  float   aperture = 0;             // lens aperture
  float   focus    = flt_max;       // focus distance
};

// Obj environment [extension].
struct obj_environment {
  string           name  = "";            // name
  frame3f          frame = identity3x4f;  // transform
  vec3f            ke    = zero3f;        // emission color
  obj_texture_info ke_txt;                // emission texture
};

// Obj procedural object [extension].
struct obj_procedural {
  string  name     = "";            // name
  frame3f frame    = identity3x4f;  // transform
  string  type     = "";            // type
  string  material = "";            // material
  float   size     = 2;             // size
  int     level    = -1;            // level of subdivision (-1 default)
};

// Obj instance [extension]
struct obj_instance {
  string  name     = "";            // name
  frame3f frame    = identity3x4f;  // transform
  string  object   = "";            // object name
  string  material = "";            // material name
};

// Obj callbacks
struct obj_callbacks {
  virtual void vert(const vec3f&) {}
  virtual void norm(const vec3f&) {}
  virtual void texcoord(const vec2f&) {}
  virtual void face(const vector<obj_vertex>&) {}
  virtual void line(const vector<obj_vertex>&) {}
  virtual void point(const vector<obj_vertex>&) {}
  virtual void object(const string&) {}
  virtual void group(const string&) {}
  virtual void usemtl(const string&) {}
  virtual void smoothing(const string&) {}
  virtual void mtllib(const string&) {}
  virtual void material(const obj_material&) {}
  virtual void camera(const obj_camera&) {}
  virtual void environmnet(const obj_environment&) {}
  virtual void instance(const obj_instance&) {}
  virtual void procedural(const obj_procedural&) {}
};

// Load obj scene
void load_obj(const string& filename, obj_callbacks& cb,
    bool nomaterials = false, bool flipv = true, bool fliptr = true);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with unordered_map
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
// SIMPLE YAML LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

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

// Load Yaml properties
bool read_yaml_property(file_wrapper& fs, string& group, string& key,
    bool& newobj, yaml_value& value);

// Write Yaml properties
void write_yaml_comment(file_wrapper& fs, const string& comment);
void write_yaml_property(file_wrapper& fs, const string& object,
    const string& key, bool newobj, const yaml_value& value);
void write_yaml_object(file_wrapper& fs, const string& object);

// type-cheked yaml value access
void get_yaml_value(const yaml_value& yaml, string& value);
void get_yaml_value(const yaml_value& yaml, bool& value);
void get_yaml_value(const yaml_value& yaml, int& value);
void get_yaml_value(const yaml_value& yaml, float& value);
void get_yaml_value(const yaml_value& yaml, vec2f& value);
void get_yaml_value(const yaml_value& yaml, vec3f& value);
void get_yaml_value(const yaml_value& yaml, mat3f& value);
void get_yaml_value(const yaml_value& yaml, frame3f& value);

// yaml value construction
yaml_value make_yaml_value(const string& value);
yaml_value make_yaml_value(bool value);
yaml_value make_yaml_value(int value);
yaml_value make_yaml_value(float value);
yaml_value make_yaml_value(const vec2f& value);
yaml_value make_yaml_value(const vec3f& value);
yaml_value make_yaml_value(const mat3f& value);
yaml_value make_yaml_value(const frame3f& value);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE PBRT LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// Pbrt value type
enum struct pbrt_value_type {
  // clang-format off
  real, integer, boolean, string, point, normal, vector, texture, color, 
  point2, vector2, spectrum
  // clang-format on
};

// Yaml value
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

// Pbrt command
enum struct pbrt_command_ {
  // clang-format off
  world_begin, world_end, attribute_begin, attribute_end,
  transform_begin, transform_end, reverse_orientation,
  set_transform, concat_transform, lookat_transform,
  object_instance, object_begin, object_end, include,      
  sampler, integrator, accelerator, film, filter, camera, shape, light,
  material, arealight, named_texture, named_medium, named_material,             
  use_material, medium_interface, active_transform,
  coordinate_system_set, coordinate_system_transform
  // clang-format on
};

// Read pbrt commands
bool read_pbrt_command(file_wrapper& fs, pbrt_command_& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values);
bool read_pbrt_command(file_wrapper& fs, pbrt_command_& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values, string& buffer);

// Write pbrt commands
void write_pbrt_comment(file_wrapper& fs, const string& comment);
void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_as_float = false);
void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
    const string& name = "", const frame3f& xform = identity3x4f);
void write_pbrt_command(file_wrapper& fs, pbrt_command_ command,
    const string& name, const string& type, const vector<pbrt_value>& values,
    bool texture_as_float = false);

// type-cheked pbrt value access
void get_pbrt_value(const pbrt_value& pbrt, string& value);
void get_pbrt_value(const pbrt_value& pbrt, bool& value);
void get_pbrt_value(const pbrt_value& pbrt, int& value);
void get_pbrt_value(const pbrt_value& pbrt, float& value);
void get_pbrt_value(const pbrt_value& pbrt, vec2f& value);
void get_pbrt_value(const pbrt_value& pbrt, vec3f& value);
void get_pbrt_value(const pbrt_value& pbrt, vector<float>& value);
void get_pbrt_value(const pbrt_value& pbrt, vector<vec2f>& value);
void get_pbrt_value(const pbrt_value& pbrt, vector<vec3f>& value);
void get_pbrt_value(const pbrt_value& pbrt, vector<int>& value);
void get_pbrt_value(const pbrt_value& pbrt, vector<vec3i>& value);
void get_pbrt_value(const pbrt_value& pbrt, pair<float, string>& value);
void get_pbrt_value(const pbrt_value& pbrt, pair<vec3f, string>& value);
template <typename T>
inline void get_pbrt_value(const vector<pbrt_value>& pbrt, const string& name,
    T& value, T def) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      get_pbrt_value(p, value);
      return;
    }
  }
  value = def;
}
template <typename T>
inline T get_pbrt_value(const vector<pbrt_value>& pbrt, const string& name,
    T def) {
  auto value = T{};
  get_pbrt_value(pbrt, name, value, def);
  return value;
}

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

}  // namespace yocto

#endif
