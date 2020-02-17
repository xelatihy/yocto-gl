//
// # Yocto/ModelIO: Tiny library for Ply/Obj/Pbrt/Yaml/glTF parsing and writing
//
// Yocto/ModelIO is a tiny library for loading and saving
// Ply/Obj/Pbrt/glTF. In Yocto/ModelIO, all model data is loaded and saved
// at once. Each format is parsed in a manner that is as close as possible to
// the original. Data can be accessed directly or via converters.
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
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto {

// Using directives
using std::function;
using std::make_unique;
using std::unique_ptr;

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
unique_ptr<ply_model> load_ply(const string& filename, string& error);
[[nodiscard]] bool    load_ply(
       const string& filename, ply_model* ply, string& error);
[[nodiscard]] bool save_ply(
    const string& filename, ply_model* ply, string& error);

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
vector<vec4f>       get_colors(ply_model* ply);
vector<float>       get_radius(ply_model* ply);
vector<vector<int>> get_faces(ply_model* ply);
vector<vec2i>       get_lines(ply_model* ply);
vector<int>         get_points(ply_model* ply);
vector<vec3i>       get_triangles(ply_model* ply);
vector<vec4i>       get_quads(ply_model* ply);
bool                has_quads(ply_model* ply);

// Create PLY
unique_ptr<ply_model> make_ply();

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
void add_colors(ply_model* ply, const vector<vec4f>& values);
void add_radius(ply_model* ply, const vector<float>& values);
void add_faces(ply_model* ply, const vector<vector<int>>& values);
void add_faces(
    ply_model* ply, const vector<vec3i>& tvalues, const vector<vec4i>& qvalues);
void add_triangles(ply_model* ply, const vector<vec3i>& values);
void add_quads(ply_model* ply, const vector<vec4i>& values);
void add_lines(ply_model* ply, const vector<vec2i>& values);
void add_points(ply_model* ply, const vector<int>& values);

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
  obj_texture_info emission_tex     = {};
  obj_texture_info ambient_tex      = {};
  obj_texture_info diffuse_tex      = {};
  obj_texture_info specular_tex     = {};
  obj_texture_info reflection_tex   = {};
  obj_texture_info transmission_tex = {};
  obj_texture_info exponent_tex     = {};
  obj_texture_info opacity_tex      = {};
  obj_texture_info bump_tex         = {};
  obj_texture_info normal_tex       = {};
  obj_texture_info displacement_tex = {};

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
  float pbr_ior           = 1.5;
  float pbr_opacity       = 1;
  vec3f pbr_volscattering = zero3f;
  float pbr_volanisotropy = 0;
  float pbr_volscale      = 0.01;

  // pbr extension textures
  obj_texture_info pbr_emission_tex      = {};
  obj_texture_info pbr_base_tex          = {};
  obj_texture_info pbr_specular_tex      = {};
  obj_texture_info pbr_roughness_tex     = {};
  obj_texture_info pbr_metallic_tex      = {};
  obj_texture_info pbr_sheen_tex         = {};
  obj_texture_info pbr_coat_tex          = {};
  obj_texture_info pbr_coatroughness_tex = {};
  obj_texture_info pbr_transmission_tex  = {};
  obj_texture_info pbr_opacity_tex       = {};
  obj_texture_info pbr_volscattering_tex = {};
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
  string           name         = "";
  frame3f          frame        = identity3x4f;
  vec3f            emission     = zero3f;
  obj_texture_info emission_tex = {};
};

// Obj model
struct obj_model {
  vector<string>           comments     = {};
  vector<obj_shape*>       shapes       = {};
  vector<obj_material*>    materials    = {};
  vector<obj_camera*>      cameras      = {};
  vector<obj_environment*> environments = {};
  ~obj_model();
};

// Load and save obj
unique_ptr<obj_model> load_obj(const string& filename, string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
[[nodiscard]] bool    load_obj(const string& filename, obj_model* obj,
       string& error, bool geom_only = false, bool split_elements = true,
       bool split_materials = false);
[[nodiscard]] bool    save_obj(
       const string& filename, obj_model* obj, string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_triangles(obj_model* obj, obj_shape* shape, vector<vec3i>& triangles,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_quads(obj_model* obj, obj_shape* shape, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_lines(obj_model* obj, obj_shape* shape, vector<vec2i>& lines,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_points(obj_model* obj, obj_shape* shape, vector<int>& points,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
void get_fvquads(obj_model* obj, obj_shape* shape, vector<vec4i>& quadspos,
    vector<vec4i>& quadsnorm, vector<vec4i>& quadstexcoord,
    vector<vec3f>& positions, vector<vec3f>& normals, vector<vec2f>& texcoords,
    vector<obj_material*>& materials, vector<int>& ematerials,
    bool flip_texcoord = false);
bool has_quads(obj_shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
void get_triangles(obj_model* obj, obj_shape* shape, int material,
    vector<vec3i>& triangles, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_quads(obj_model* obj, obj_shape* shape, int material,
    vector<vec4i>& quads, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_lines(obj_model* obj, obj_shape* shape, int material,
    vector<vec2i>& lines, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_points(obj_model* obj, obj_shape* shape, int material,
    vector<int>& points, vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texcoords, bool flip_texcoord = false);
vector<obj_material*> get_materials(obj_model* obj, obj_shape* shape);

// Create OBJ
unique_ptr<obj_model> make_obj();
obj_camera*           add_camera(obj_model* obj);
obj_material*         add_material(obj_model* obj);
obj_environment*      add_environment(obj_model* obj);
obj_shape*            add_shape(obj_model* obj);

// Add obj shape
void add_triangles(obj_model* obj, const string& name,
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
void add_quads(obj_model* obj, const string& name, const vector<vec4i>& quads,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
void add_lines(obj_model* obj, const string& name, const vector<vec2i>& lines,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
void add_points(obj_model* obj, const string& name, const vector<int>& points,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texcoords, const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);
void add_fvquads(obj_model* obj, const string& name,
    const vector<vec4i>& quadspos, const vector<vec4i>& quadsnorm,
    const vector<vec4i>& quadstexcoord, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texcoords,
    const vector<obj_material*>& materials = {},
    const vector<int>& ematerials = {}, const vector<frame3f>& instances = {},
    bool flip_texcoord = false);

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
// SIMPLE PBRT LOADER AND WRITER
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
  vec3f  emission        = zero3f;
  vec3f  color           = zero3f;
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
  vec3f  volmeanfreepath = zero3f;
  vec3f  volscatter      = zero3f;
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
  vec3f   emission = zero3f;
  vec3f   from     = zero3f;
  vec3f   to       = zero3f;
  bool    distant  = false;
  // arealight approximation
  vec3f         area_emission  = zero3f;
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
  vec3f   emission     = zero3f;
  string  emission_tex = "";
};

// Pbrt model
struct pbrt_model {
  vector<string>            comments     = {};
  vector<pbrt_camera*>      cameras      = {};
  vector<pbrt_shape*>       shapes       = {};
  vector<pbrt_environment*> environments = {};
  vector<pbrt_light*>       lights       = {};
  vector<pbrt_material*>    materials    = {};
  ~pbrt_model();
};

// Load/save pbrt
unique_ptr<pbrt_model> load_pbrt(const string& filename, string& error);
[[nodiscard]] bool     load_pbrt(
        const string& filename, pbrt_model* pbrt, string& error);
[[nodiscard]] bool save_pbrt(const string& filename, pbrt_model* pbrt,
    string& error, bool ply_meshes = false);

// Create pbrt
unique_ptr<pbrt_model> make_pbrt();
pbrt_camera*           add_camera(pbrt_model* pbrt);
pbrt_shape*            add_shape(pbrt_model* pbrt);
pbrt_material*         add_material(pbrt_model* pbrt);
pbrt_environment*      add_environment(pbrt_model* pbrt);
pbrt_light*            add_light(pbrt_model* pbrt);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE GLTF LOADER DECLARATIONS
// -----------------------------------------------------------------------------
namespace yocto {

struct gltf_camera {
  string          name   = "";
  vector<frame3f> frames = {};
  bool            ortho  = false;
  float           yfov   = 45 * pif / 180;
  float           aspect = 1;
};
struct gltf_texture {
  string name     = "";
  string filename = "";
};
struct gltf_material {
  string        name         = "";
  vec3f         emission     = {0, 0, 0};
  vec3f         color        = {0, 0, 0};
  float         opacity      = 1;
  float         metallic     = 0;
  float         roughness    = 1;
  gltf_texture* emission_tex = nullptr;
  gltf_texture* color_tex    = nullptr;
  gltf_texture* metallic_tex = nullptr;
  gltf_texture* normal_tex   = nullptr;
};
struct gltf_primitive {
  gltf_material* material  = nullptr;
  vector<vec3f>  positions = {};
  vector<vec3f>  normals   = {};
  vector<vec2f>  texcoords = {};
  vector<vec4f>  colors    = {};
  vector<float>  radius    = {};
  vector<vec4f>  tangents  = {};
  vector<vec3i>  triangles = {};
  vector<vec2i>  lines     = {};
  vector<int>    points    = {};
};
struct gltf_mesh {
  string                  name       = "";
  vector<frame3f>         frames     = {};
  vector<gltf_primitive*> primitives = {};
};
struct gltf_model {
  vector<gltf_camera*>    cameras    = {};
  vector<gltf_mesh*>      meshes     = {};
  vector<gltf_primitive*> primitives = {};
  vector<gltf_texture*>   textures   = {};
  vector<gltf_material*>  materials  = {};
  ~gltf_model();
};

// Load gltf file.
unique_ptr<gltf_model> load_gltf(const string& filename, string& error);
[[nodiscard]] bool     load_gltf(
        const string& filename, gltf_model* gltf, string& error);

}  // namespace yocto

#endif
