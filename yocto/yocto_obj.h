//
// # Yocto/Obj: Tiny library for Obj parsing and writing
//
// Yocto/Obj is a tiny library for loading and saving Obj.
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER AND WRITER
// -----------------------------------------------------------------------------
namespace yocto::obj {

// Using directives
using namespace yocto::math;

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
bool load_obj(const string& filename, obj_model* obj, string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
bool save_obj(const string& filename, obj_model* obj, string& error);

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
obj_camera*      add_camera(obj_model* obj);
obj_material*    add_material(obj_model* obj);
obj_environment* add_environment(obj_model* obj);
obj_shape*       add_shape(obj_model* obj);

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

}  // namespace yocto::obj

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for vector for use with hash_map
template <>
struct hash<yocto::obj::obj_vertex> {
  size_t operator()(const yocto::obj::obj_vertex& v) const {
    static const std::hash<int> hasher = std::hash<int>();
    auto                        h      = (size_t)0;
    h ^= hasher(v.position) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.normal) + 0x9e3779b9 + (h << 6) + (h >> 2);
    h ^= hasher(v.texcoord) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
  }
};

}  // namespace std

#endif
