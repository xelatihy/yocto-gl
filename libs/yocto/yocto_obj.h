//
// # Yocto/Obj: Tiny library for Obj parsing and writing
//
// Yocto/Obj is a tiny library for loading and saving Obj.
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "yocto_math.h"

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
  std::string path  = "";     // file path
  bool        clamp = false;  // clamp to edge
  float       scale = 1;      // scale for bump/displacement

  // Properties not explicitly handled.
  std::unordered_map<std::string, std::vector<float>> props;

  obj_texture() {}
  obj_texture(const char* path) : path{path} {}
  obj_texture(const std::string& path) : path{path} {}
};

// Obj element
struct obj_element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj material
struct obj_material {
  // material name and type
  std::string name  = "";
  int         illum = 0;

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
  std::string                name      = "";
  std::vector<vec3f>         positions = {};
  std::vector<vec3f>         normals   = {};
  std::vector<vec2f>         texcoords = {};
  std::vector<obj_material*> materials = {};
  std::vector<obj_vertex>    vertices  = {};
  std::vector<obj_element>   faces     = {};
  std::vector<obj_element>   lines     = {};
  std::vector<obj_element>   points    = {};
  std::vector<frame3f>       instances = {};
};

// Obj camera
struct obj_camera {
  std::string name     = "";
  frame3f     frame    = identity3x4f;
  bool        ortho    = false;
  float       width    = 0.036;
  float       height   = 0.028;
  float       lens     = 0.50;
  float       focus    = 0;
  float       aperture = 0;
};

// Obj environment
struct obj_environment {
  std::string name         = "";
  frame3f     frame        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  obj_texture emission_tex = {};
};

// Obj model
struct obj_model {
  std::vector<std::string>      comments     = {};
  std::vector<obj_shape*>       shapes       = {};
  std::vector<obj_material*>    materials    = {};
  std::vector<obj_camera*>      cameras      = {};
  std::vector<obj_environment*> environments = {};
  ~obj_model();
};

// Load and save obj
bool load_obj(const std::string& filename, obj_model* obj, std::string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
bool save_obj(const std::string& filename, obj_model* obj, std::string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_triangles(const obj_shape* shape, std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj_material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_quads(const obj_shape* shape, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj_material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_lines(const obj_shape* shape, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj_material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_points(const obj_shape* shape, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj_material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_fvquads(const obj_shape* shape, std::vector<vec4i>& quadspos,
    std::vector<vec4i>& quadsnorm, std::vector<vec4i>& quadstexcoord,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj_material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
bool has_quads(obj_shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
void get_triangles(const obj_shape* shape, int material,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);
void get_quads(const obj_shape* shape, int material, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_lines(const obj_shape* shape, int material, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_points(const obj_shape* shape, int material, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);

// Create OBJ
obj_camera*      add_camera(obj_model* obj);
obj_material*    add_material(obj_model* obj);
obj_environment* add_environment(obj_model* obj);
obj_shape*       add_shape(obj_model* obj);

// Add obj shape
void set_triangles(obj_shape* shape, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_quads(obj_shape* shape, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_lines(obj_shape* shape, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_points(obj_shape* shape, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_fvquads(obj_shape* shape, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_materials(
    obj_shape* shape, const std::vector<obj_material*>& materials);
void set_instances(obj_shape* shape, const std::vector<frame3f>& instances);

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for std::vector for use with hash_map
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

#endif
