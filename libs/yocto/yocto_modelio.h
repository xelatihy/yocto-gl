//
// # Yocto/ModelIO: Serialization for Obj, Ply and Pbrt models
//
// Yocto/ModelIO is a collection of utilities for loading and saving scenes
// and meshes in Ply, Obj and Pbrt formats.
// Yocto/ModelIO is implemented in `yocto_modelio.h` and `yocto_modelio.cpp`.
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

#ifndef _YOCTO_MODELIO_H_
#define _YOCTO_MODELIO_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <memory>
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
    const array<string, 2>& properties, vector<vec4f>& values);
bool get_values(ply_model* ply, const string& element,
    const array<string, 3>& properties, vector<vec3f>& values);
bool get_values(ply_model* ply, const string& element,
    const array<string, 4>& properties, vector<vec4f>& values);
bool get_values(ply_model* ply, const string& element,
    const array<string, 12>& properties, vector<frame3f>& values);

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
bool get_colors(ply_model* ply, vector<vec4f>& values);
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
    const array<string, 2>& properties, const vector<vec2f>& values);
bool add_values(ply_model* ply, const string& element,
    const array<string, 3>& properties, const vector<vec3f>& values);
bool add_values(ply_model* ply, const string& element,
    const array<string, 4>& properties, const vector<vec4f>& values);
bool add_values(ply_model* ply, const string& element,
    const array<string, 12>& properties, const vector<frame3f>& values);

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
bool add_colors(ply_model* ply, const vector<vec4f>& values);
bool add_radius(ply_model* ply, const vector<float>& values);
bool add_faces(ply_model* ply, const vector<vector<int>>& values);
bool add_faces(
    ply_model* ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
bool add_triangles(ply_model* ply, const vector<vec3i>& values);
bool add_quads(ply_model* ply, const vector<vec4i>& values);
bool add_lines(ply_model* ply, const vector<vec2i>& values);
bool add_points(ply_model* ply, const vector<int>& values);

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

// Obj texture information.
struct obj_texture {
  string path  = "";     // file path
  bool   clamp = false;  // clamp to edge
  float  scale = 1;      // scale for bump/displacement

  obj_texture() {}
  explicit obj_texture(const char* path) : path{path} {}
  explicit obj_texture(const string& path) : path{path} {}
};

// Obj element
struct obj_element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj material
struct obj_material {
  // material name and type
  string name  = "";
  int    illum = 0;

  // material colors and values
  vec3f emission     = {0, 0, 0};
  vec3f ambient      = {0, 0, 0};
  vec3f diffuse      = {0, 0, 0};
  vec3f specular     = {0, 0, 0};
  vec3f reflection   = {0, 0, 0};
  vec3f transmission = {0, 0, 0};
  float exponent     = 10;
  float ior          = 1.5;
  float opacity      = 1;

  // material textures
  obj_texture emission_tex     = {};
  obj_texture ambient_tex      = {};
  obj_texture diffuse_tex      = {};
  obj_texture specular_tex     = {};
  obj_texture reflection_tex   = {};
  obj_texture transmission_tex = {};
  obj_texture exponent_tex     = {};
  obj_texture opacity_tex      = {};
  obj_texture bump_tex         = {};
  obj_texture normal_tex       = {};
  obj_texture displacement_tex = {};

  // pbrt extension values
  bool  as_pbr            = false;
  vec3f pbr_emission      = {0, 0, 0};
  vec3f pbr_base          = {0, 0, 0};
  float pbr_specular      = 0;
  float pbr_roughness     = 0;
  float pbr_metallic      = 0;
  float pbr_sheen         = 0;
  float pbr_coat          = 0;
  float pbr_coatroughness = 0;
  float pbr_transmission  = 0;
  float pbr_translucency  = 0;
  float pbr_ior           = 1.5;
  float pbr_opacity       = 1;
  vec3f pbr_volscattering = {0, 0, 0};
  float pbr_volanisotropy = 0;
  float pbr_volscale      = 0.01;
  bool  pbr_thin          = true;

  // pbr extension textures
  obj_texture pbr_emission_tex      = {};
  obj_texture pbr_base_tex          = {};
  obj_texture pbr_specular_tex      = {};
  obj_texture pbr_roughness_tex     = {};
  obj_texture pbr_metallic_tex      = {};
  obj_texture pbr_sheen_tex         = {};
  obj_texture pbr_coat_tex          = {};
  obj_texture pbr_coatroughness_tex = {};
  obj_texture pbr_transmission_tex  = {};
  obj_texture pbr_translucency_tex  = {};
  obj_texture pbr_opacity_tex       = {};
  obj_texture pbr_volscattering_tex = {};
  obj_texture pbr_bump_tex          = {};
  obj_texture pbr_normal_tex        = {};
  obj_texture pbr_displacement_tex  = {};
};

// Obj shape
struct obj_shape {
  string                name      = "";
  vector<vec3f>         positions = {};
  vector<vec3f>         normals   = {};
  vector<vec2f>         texcoords = {};
  vector<obj_material*> materials = {};
  vector<obj_vertex>    vertices  = {};
  vector<obj_element>   faces     = {};
  vector<obj_element>   lines     = {};
  vector<obj_element>   points    = {};
  vector<frame3f>       instances = {};
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
  string      name         = "";
  frame3f     frame        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  obj_texture emission_tex = {};
};

// Obj model
struct obj_scene {
  vector<string>           comments     = {};
  vector<obj_shape*>       shapes       = {};
  vector<obj_material*>    materials    = {};
  vector<obj_camera*>      cameras      = {};
  vector<obj_environment*> environments = {};
  ~obj_scene();
};

// Load and save obj
bool load_obj(const string& filename, obj_scene* obj, string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
bool save_obj(const string& filename, obj_scene* obj, string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_triangles(const obj_shape* shape, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_quads(const obj_shape* shape, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_lines(const obj_shape* shape, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_points(const obj_shape* shape, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_fvquads(const obj_shape* shape, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
bool has_quads(obj_shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
void get_triangles(const obj_shape* shape, int material,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_quads(const obj_shape* shape, int material, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flip_texcoord = false);
void get_lines(const obj_shape* shape, int material, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flip_texcoord = false);
void get_points(const obj_shape* shape, int material, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    bool flip_texcoord = false);

// Create OBJ
obj_camera*      add_camera(obj_scene* obj);
obj_material*    add_material(obj_scene* obj);
obj_environment* add_environment(obj_scene* obj);
obj_shape*       add_shape(obj_scene* obj);

// Add obj shape
void set_triangles(obj_shape* shape, const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_quads(obj_shape* shape, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_lines(obj_shape* shape, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_points(obj_shape* shape, const vector<int>& points,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_fvquads(obj_shape* shape, const vector<vec4i>& quadspos,
    const vector<vec4i>& quadsnorm, const vector<vec4i>& quadstexcoord,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_materials(obj_shape* shape, const vector<obj_material*>& materials);
void set_instances(obj_shape* shape, const vector<frame3f>& instances);

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

// Pbrt camera
struct pbrt_camera {
  // camera parameters
  frame3f frame      = identity3x4f;
  frame3f frend      = identity3x4f;
  vec2i   resolution = {0, 0};
  float   lens       = 0;
  float   aspect     = 0;
  float   focus      = 0;
  float   aperture   = 0;
};

// Pbrt material
struct pbrt_material {
  // material parameters
  string name            = "";
  vec3f  emission        = {0, 0, 0};
  vec3f  color           = {0, 0, 0};
  float  specular        = 0;
  float  metallic        = 0;
  float  transmission    = 0;
  float  roughness       = 0;
  float  ior             = 1.5;
  float  opacity         = 1;
  string color_tex       = "";
  string opacity_tex     = "";
  string alpha_tex       = "";
  bool   thin            = true;
  vec3f  volmeanfreepath = {0, 0, 0};
  vec3f  volscatter      = {0, 0, 0};
  float  volscale        = 0.01;
};

// Pbrt shape
struct pbrt_shape {
  // frames
  frame3f         frame     = identity3x4f;
  frame3f         frend     = identity3x4f;
  vector<frame3f> instances = {};
  vector<frame3f> instaends = {};
  // shape
  string        filename_ = "";
  vector<vec3f> positions = {};
  vector<vec3f> normals   = {};
  vector<vec2f> texcoords = {};
  vector<vec3i> triangles = {};
  // material
  pbrt_material* material = nullptr;
};

// Pbrt lights
struct pbrt_light {
  // light parameters
  frame3f frame    = identity3x4f;
  frame3f frend    = identity3x4f;
  vec3f   emission = {0, 0, 0};
  vec3f   from     = {0, 0, 0};
  vec3f   to       = {0, 0, 0};
  bool    distant  = false;
  // arealight approximation
  vec3f         area_emission  = {0, 0, 0};
  frame3f       area_frame     = identity3x4f;
  frame3f       area_frend     = identity3x4f;
  vector<vec3i> area_triangles = {};
  vector<vec3f> area_positions = {};
  vector<vec3f> area_normals   = {};
};
struct pbrt_environment {
  // environment approximation
  frame3f frame        = identity3x4f;
  frame3f frend        = identity3x4f;
  vec3f   emission     = {0, 0, 0};
  string  emission_tex = "";
};

// Pbrt model
struct pbrt_scene {
  // pbrt data
  vector<string>            comments     = {};
  vector<pbrt_camera*>      cameras      = {};
  vector<pbrt_shape*>       shapes       = {};
  vector<pbrt_environment*> environments = {};
  vector<pbrt_light*>       lights       = {};
  vector<pbrt_material*>    materials    = {};

  // cleanup
  ~pbrt_scene();
};

// Load/save pbrt
bool load_pbrt(const string& filename, pbrt_scene* pbrt, string& error);
bool save_pbrt(const string& filename, pbrt_scene* pbrt, string& error,
    bool ply_meshes = false);

// Create pbrt
pbrt_camera*      add_camera(pbrt_scene* pbrt);
pbrt_shape*       add_shape(pbrt_scene* pbrt);
pbrt_material*    add_material(pbrt_scene* pbrt);
pbrt_environment* add_environment(pbrt_scene* pbrt);
pbrt_light*       add_light(pbrt_scene* pbrt);

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARDS COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

using obj_model [[deprecated]]  = obj_scene;
using pbrt_model [[deprecated]] = pbrt_scene;

}  // namespace yocto

#endif
