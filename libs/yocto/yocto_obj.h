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
namespace yocto::obj {

// OBJ vertex
struct vertex {
  int position = 0;
  int texcoord = 0;
  int normal   = 0;
};

inline bool operator==(const vertex& a, const vertex& b) {
  return a.position == b.position && a.texcoord == b.texcoord &&
         a.normal == b.normal;
}

// Obj texture information.
struct texture {
  std::string path  = "";     // file path
  bool        clamp = false;  // clamp to edge
  float       scale = 1;      // scale for bump/displacement

  // Properties not explicitly handled.
  std::unordered_map<std::string, std::vector<float>> props;

  texture() {}
  texture(const char* path) : path{path} {}
  texture(const std::string& path) : path{path} {}
};

// Obj element
struct element {
  uint8_t size     = 0;
  uint8_t material = 0;
};

// Obj material
struct material {
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
  texture emission_tex     = {};
  texture ambient_tex      = {};
  texture diffuse_tex      = {};
  texture specular_tex     = {};
  texture reflection_tex   = {};
  texture transmission_tex = {};
  texture exponent_tex     = {};
  texture opacity_tex      = {};
  texture bump_tex         = {};
  texture normal_tex       = {};
  texture displacement_tex = {};

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
  texture pbr_emission_tex      = {};
  texture pbr_base_tex          = {};
  texture pbr_specular_tex      = {};
  texture pbr_roughness_tex     = {};
  texture pbr_metallic_tex      = {};
  texture pbr_sheen_tex         = {};
  texture pbr_coat_tex          = {};
  texture pbr_coatroughness_tex = {};
  texture pbr_transmission_tex  = {};
  texture pbr_translucency_tex  = {};
  texture pbr_opacity_tex       = {};
  texture pbr_volscattering_tex = {};
  texture pbr_bump_tex          = {};
  texture pbr_normal_tex        = {};
  texture pbr_displacement_tex  = {};
};

// Obj shape
struct shape {
  std::string            name      = "";
  std::vector<vec3f>     positions = {};
  std::vector<vec3f>     normals   = {};
  std::vector<vec2f>     texcoords = {};
  std::vector<material*> materials = {};
  std::vector<vertex>    vertices  = {};
  std::vector<element>   faces     = {};
  std::vector<element>   lines     = {};
  std::vector<element>   points    = {};
  std::vector<frame3f>   instances = {};
};

// Obj camera
struct camera {
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
struct environment {
  std::string name         = "";
  frame3f     frame        = identity3x4f;
  vec3f       emission     = {0, 0, 0};
  texture     emission_tex = {};
};

// Obj model
struct model {
  std::vector<std::string>       comments     = {};
  std::vector<obj::shape*>       shapes       = {};
  std::vector<obj::material*>    materials    = {};
  std::vector<obj::camera*>      cameras      = {};
  std::vector<obj::environment*> environments = {};
  ~model();
};

// Load and save obj
bool load_obj(const std::string& filename, obj::model* obj, std::string& error,
    bool geom_only = false, bool split_elements = true,
    bool split_materials = false);
bool save_obj(const std::string& filename, obj::model* obj, std::string& error);

// Get obj shape. Obj is a facevarying format, so vertices might be duplicated.
// to ensure that no duplication occurs, either use the facevarying interface,
// or set `no_vertex_duplication`. In the latter case, the code will fallback
// to position only if duplication occurs.
void get_triangles(const obj::shape* shape, std::vector<vec3i>& triangles,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_quads(const obj::shape* shape, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_lines(const obj::shape* shape, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_points(const obj::shape* shape, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
void get_fvquads(const obj::shape* shape, std::vector<vec4i>& quadspos,
    std::vector<vec4i>& quadsnorm, std::vector<vec4i>& quadstexcoord,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, std::vector<obj::material*>& materials,
    std::vector<int>& ematerials, bool flip_texcoord = false);
bool has_quads(obj::shape* shape);

// Get obj shape by extracting the elements beloing to only one material.
void get_triangles(const obj::shape* shape, int material,
    std::vector<vec3i>& triangles, std::vector<vec3f>& positions,
    std::vector<vec3f>& normals, std::vector<vec2f>& texcoords,
    bool flip_texcoord = false);
void get_quads(const obj::shape* shape, int material, std::vector<vec4i>& quads,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_lines(const obj::shape* shape, int material, std::vector<vec2i>& lines,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);
void get_points(const obj::shape* shape, int material, std::vector<int>& points,
    std::vector<vec3f>& positions, std::vector<vec3f>& normals,
    std::vector<vec2f>& texcoords, bool flip_texcoord = false);

// Create OBJ
obj::camera*      add_camera(obj::model* obj);
obj::material*    add_material(obj::model* obj);
obj::environment* add_environment(obj::model* obj);
obj::shape*       add_shape(obj::model* obj);

// Add obj shape
void set_triangles(obj::shape* shape, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_quads(obj::shape* shape, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_lines(obj::shape* shape, const std::vector<vec2i>& lines,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords, const std::vector<int>& ematerials = {},
    bool flip_texcoord = false);
void set_points(obj::shape* shape, const std::vector<int>& points,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_fvquads(obj::shape* shape, const std::vector<vec4i>& quadspos,
    const std::vector<vec4i>& quadsnorm,
    const std::vector<vec4i>& quadstexcoord,
    const std::vector<vec3f>& positions, const std::vector<vec3f>& normals,
    const std::vector<vec2f>& texcoords,
    const std::vector<int>& ematerials = {}, bool flip_texcoord = false);
void set_materials(
    obj::shape* shape, const std::vector<obj::material*>& materials);
void set_instances(obj::shape* shape, const std::vector<frame3f>& instances);

}  // namespace yocto::obj

// -----------------------------------------------------------------------------
// HELPER FOR DICTIONARIES
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for std::vector for use with hash_map
template <>
struct hash<yocto::obj::vertex> {
  size_t operator()(const yocto::obj::vertex& v) const {
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
