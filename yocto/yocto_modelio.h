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

// Load/save pbrt
void load_pbrt(const string& filename, pbrt_model& pbrt);
void save_pbrt(const string& filename, const pbrt_model& pbrt);

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
  coordinate_system_set, coordinate_system_transform
  // clang-format on
};

// Read pbrt commands
bool read_pbrt_command(file_wrapper& fs, pbrt_command& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values);
bool read_pbrt_command(file_wrapper& fs, pbrt_command& command, string& name,
    string& type, frame3f& xform, vector<pbrt_value>& values, string& buffer);

// Write pbrt commands
void write_pbrt_comment(file_wrapper& fs, const string& comment);
void write_pbrt_command(file_wrapper& fs, pbrt_command command,
    const string& name, const string& type, const frame3f& xform,
    const vector<pbrt_value>& values, bool texture_as_float = false);
void write_pbrt_command(file_wrapper& fs, pbrt_command command,
    const string& name = "", const frame3f& xform = identity3x4f);
void write_pbrt_command(file_wrapper& fs, pbrt_command command,
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
inline void get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T& value, T def) {
  for (auto& p : pbrt) {
    if (p.name == name) {
      get_pbrt_value(p, value);
      return;
    }
  }
  value = def;
}
template <typename T>
inline T get_pbrt_value(
    const vector<pbrt_value>& pbrt, const string& name, T def) {
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
pbrt_value make_pbrt_value(const string& name, const vector<vec2f>& value,
    pbrt_value_type type = pbrt_value_type::point2);
pbrt_value make_pbrt_value(const string& name, const vector<vec3f>& value,
    pbrt_value_type type = pbrt_value_type::point);
pbrt_value make_pbrt_value(const string& name, const vector<vec3i>& value,
    pbrt_value_type type = pbrt_value_type::integer);

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
