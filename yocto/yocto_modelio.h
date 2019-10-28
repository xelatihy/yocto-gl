//
// # Yocto/ModelIO: L<ow-level library for OBJ/Pbrt/Yaml parsing and writing
//
// Yocto/ModelIO is a collection of simple parsers for
// Wavefront Obj, Pbrt, and Yaml formats. The prasers are designed for large
// files and do keep a copy of the model in memory.
// Yocto/ModelIO provides fast/low-level access to model data and requires some
// familiarity with the formats to use effectively. For a higher level
// interface, consider using Yocto/Shape's `load_shape()` and `save_shape()`,
// or Yocto/SceneIO's `load_scene()` and `save_scene()`
//
// Yocto/ModelIO also support writing Obj/Yaml files again without keeping
// a copy of the model but instead writing elements directly after each call.
// Error reporting is done by throwing `std::runtime_error` exceptions.
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

// Load/save yaml
void load_yaml(const string& filename, yaml_model& yaml);
void save_yaml(const string& filename, const yaml_model& yaml);

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
// SIMPLE GLTF LOADER
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

void load_gltf(const string& filename, gltf_model& gltf);

}  // namespace yocto

#endif
