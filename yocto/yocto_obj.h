//
// # Yocto/Obj: Tiny library for OBJ parsing/writing
//
// Yocto/Obj is a simple Wavefront OBJ parser that works with callbacks with
// support for a few extensions such as camera, instances and environments.
// The praser is designed for large files and does keep a copy of the model.
//
// Yocto/Obj also support writing OBJ files again without keeping a copy of the
// model but instead writing elements directly after each call.
//
// Error reporting is done by throwing `std::runtime_error` exceptions.
//
// ## Parse an OBJ file
//
// 1. define callbacks in `obj_callbacks` structure
// 2. run the parse with `load_obj()`
//
// ## Write amn OBJ file
//
// 1. use `init_obj_streams()` to initialize the file streams for weriting
// 2. use the `write_obj_XXX()` function to write single Obj elements
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

#ifndef _YOCTO_OBJ_H_
#define _YOCTO_OBJ_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER
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

// Obj command
enum struct obj_command {
  // clang-format off
  vertex, normal, texcoord,         // data in value
  face, line, point,                // data in elements
  object, group, usemtl, smoothing, // data in name
  mtllib, objxlib,                  // data in name
  // clang-format on
};
// Mtl command
enum struct mtl_command {
  // clang-format off
  material,         // data in material
  // clang-format on
};
// Objx command
enum struct objx_command {
  // clang-format off
  camera,       // data in camera
  environment,  // data in environment
  instance,     // data in instance
  procedural,   // data in procedural
  // clang-format on
};

// Write obj elements
void write_obj_comment(FILE* fs, const string& comment);
void write_obj_command(FILE* fs, obj_command command, const vec3f& value,
    const string& name, const vector<obj_vertex>& elements);
void write_mtl_command(FILE* fs, mtl_command command, 
  const obj_material& material);
void write_objx_command(FILE* fs, objx_command command, 
  const obj_camera& camera, const obj_environment& environment, 
  const obj_instance& instance, const obj_procedural& procedural);

// Write obj elements
void write_obj_vertex(FILE* fs, const vec3f& p);
void write_obj_normal(FILE* fs, const vec3f& n);
void write_obj_texcoord(FILE* fs, const vec2f& t);
void write_obj_face(FILE* fs, const vector<obj_vertex>& verts);
void write_obj_face(FILE* fs, const obj_vertex& vert1, const obj_vertex& vert2,
    const obj_vertex& vert3);
void write_obj_face(FILE* fs, const obj_vertex& vert1, const obj_vertex& vert2,
    const obj_vertex& vert3, const obj_vertex& vert4);
void write_obj_line(FILE* fs, const vector<obj_vertex>& verts);
void write_obj_line(FILE* fs, const obj_vertex& vert1, const obj_vertex& vert2);
void write_obj_point(FILE* fs, const vector<obj_vertex>& verts);
void write_obj_point(FILE* fs, const obj_vertex& vert);
void write_obj_object(FILE* fs, const string& name);
void write_obj_group(FILE* fs, const string& name);
void write_obj_usemtl(FILE* fs, const string& name);
void write_obj_smoothing(FILE* fs, const string& name);
void write_obj_mtllib(FILE* fs, const string& name);
void write_mtl_material(FILE* fs, const obj_material& material);
void write_objx_camera(FILE* fs, const obj_camera& camera);
void write_objx_environmnet(FILE* fs, const obj_environment& environment);
void write_objx_instance(FILE* fs, const obj_instance& instamce);
void write_obj_procedural(FILE* fs, const obj_procedural& procedural);

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

#endif
