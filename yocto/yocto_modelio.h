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

// pbrt pbrt_spectrum as rgb color
struct pbrt_spectrum3f {
  float x = 0;
  float y = 0;
  float z = 0;

  pbrt_spectrum3f() : x{0}, y{0}, z{0} {}
  pbrt_spectrum3f(float x, float y, float z) : x{x}, y{y}, z{z} {}
  explicit pbrt_spectrum3f(float v) : x{v}, y{v}, z{v} {}
  explicit operator vec3f() const { return {x, y, z}; };

  float&       operator[](int i) { return (&x)[i]; }
  const float& operator[](int i) const { return (&x)[i]; }
};

// pbrt cameras
struct pbrt_camera {
  struct perspective_t {
    float fov              = 90;
    float frameaspectratio = -1;  // or computed from film
    float lensradius       = 0;
    float focaldistance    = 1e30;
    vec4f screenwindow     = {-1, 1, -1, 1};
    float shutteropen      = 0;
    float shutterclose     = 1;
  };
  struct orthographic_t {
    float frameaspectratio = -1;  // or computed from film
    float lensradius       = 0;
    float focaldistance    = 1e30;
    vec4f screenwindow     = {-1, 1, -1, 1};
    float shutteropen      = 0;
    float shutterclose     = 1;
  };
  struct environment_t {
    float shutteropen  = 0;
    float shutterclose = 1;
  };
  struct realistic_t {
    string lensfile           = "";
    float  aperturediameter   = 1;
    float  focusdistance      = 10;
    bool   simpleweighting    = true;
    float  shutteropen        = 0;
    float  shutterclose       = 1;
    float  approx_focallength = 0;
  };
  enum struct type_t { perspective, orthographic, environment, realistic };
  type_t         type         = type_t::perspective;
  perspective_t  perspective  = {};
  orthographic_t orthographic = {};
  environment_t  environment  = {};
  realistic_t    realistic    = {};
};

// pbrt samplers
struct pbrt_sampler {
  struct random_t {
    int pixelsamples = 16;
  };
  struct halton_t {
    int pixelsamples = 16;
  };
  struct sobol_t {
    int pixelsamples = 16;
  };
  struct zerotwosequence_t {
    int pixelsamples = 16;
  };
  struct maxmindist_t {
    int pixelsamples = 16;
  };
  struct stratified_t {
    bool jitter   = true;
    int  xsamples = 2;
    int  ysamples = 2;
  };
  enum struct type_t {
    random,
    halton,
    sobol,
    zerotwosequence,
    maxmindist,
    stratified
  };
  type_t            type            = type_t::random;
  random_t          random          = {};
  halton_t          halton          = {};
  sobol_t           sobol           = {};
  zerotwosequence_t zerotwosequence = {};
  maxmindist_t      maxmindist      = {};
  stratified_t      stratified      = {};
};

// pbrt film
struct pbrt_film {
  struct image_t {
    int    xresolution        = 640;
    int    yresolution        = 480;
    vec4f  cropwindow         = {0, 1, 0, 1};
    float  scale              = 1;
    float  maxsampleluminance = flt_max;
    float  diagonal           = 35;
    string filename           = "pbrt.exr";
  };
  enum struct type_t { image };
  type_t  type  = type_t::image;
  image_t image = {};
};

// pbrt filters
struct pbrt_filter {
  struct box_t {
    float xwidth = 0.5;
    float ywidth = 0.5;
  };
  struct gaussian_t {
    float xwidth = 2;
    float ywidth = 2;
    float alpha  = 2;
  };
  struct mitchell_t {
    float xwidth = 2;
    float ywidth = 2;
    float B      = 1.0f / 3.0f;
    float C      = 1.0f / 3.0f;
  };
  struct sinc_t {
    float xwidth = 4;
    float ywidth = 4;
    float tau    = 3;
  };
  struct triangle_t {
    float xwidth = 2;
    float ywidth = 2;
  };
  enum struct type_t { box, gaussian, mitchell, sinc, triangle };
  type_t     type     = type_t::box;
  box_t      box      = {};
  gaussian_t gaussian = {};
  mitchell_t mitchell = {};
  sinc_t     sinc     = {};
  triangle_t triangle = {};
};

// pbrt integrators
struct pbrt_integrator {
  struct path_t {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int                   maxdepth            = 5;
    vec4i                 pixelbounds         = {0, 0, int_max, int_max};
    float                 rrthreshold         = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
  };
  struct volpath_t {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int                   maxdepth            = 5;
    vec4i                 pixelbounds         = {0, 0, int_max, int_max};
    float                 rrthreshold         = 1;
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::spatial;
  };
  struct bdpt_t {
    enum struct lightsamplestrategy_t { uniform, power, spatial };
    int                   maxdepth            = 5;
    vec4i                 pixelbounds         = {0, 0, int_max, int_max};
    lightsamplestrategy_t lightsamplestrategy = lightsamplestrategy_t::power;
    bool                  visualizestrategies = false;
    bool                  visualizeweights    = false;
  };
  struct directlighting_t {
    enum struct strategy_t { all, one };
    strategy_t strategy    = strategy_t::all;
    int        maxdepth    = 5;
    vec4i      pixelbounds = {0, 0, int_max, int_max};
  };
  struct mlt_t {
    int   maxdepth             = 5;
    vec4i pixelbounds          = {0, 0, int_max, int_max};
    int   bootstrapsamples     = 100000;
    int   chains               = 1000;
    int   mutationsperpixel    = 100;
    float largestepprobability = 0.3;
    float sigma                = 0.01;
  };
  struct sppm_t {
    int   maxdepth            = 5;
    vec4i pixelbounds         = {0, 0, int_max, int_max};
    int   iterations          = 64;
    int   photonsperiteration = -1;
    int   imagewritefrequency = pow2(31);
    float radius              = 5;
  };
  struct whitted_t {
    int   maxdepth    = 5;
    vec4i pixelbounds = {0, 0, int_max, int_max};
  };
  enum struct type_t {
    path,
    volpath,
    bdpt,
    directlighting,
    mlt,
    sppm,
    whitted
  };
  type_t           type           = type_t::path;
  path_t           path           = {};
  volpath_t        volpath        = {};
  bdpt_t           bdpt           = {};
  directlighting_t directlighting = {};
  mlt_t            mlt            = {};
  sppm_t           sppm           = {};
  whitted_t        whitted        = {};
};

// pbrt accellerators
struct pbrt_accelerator {
  struct bvh_t {
    enum struct splitmethod_t { sah, equal, middle, hlbvh };
    int           maxnodeprims = 4;
    splitmethod_t splitmethod  = splitmethod_t::sah;
  };
  struct kdtree_t {
    int   intersectcost = 80;
    int   traversalcost = 1;
    float emptybonus    = 0.2;
    int   maxprims      = 1;
    int   maxdepth      = -1;
  };
  enum struct type_t { bvh, kdtree };
  type_t   type   = type_t::bvh;
  bvh_t    bvh    = {};
  kdtree_t kdtree = {};
};

// pbrt texture or value
struct pbrt_textured1f {
  float  value   = 0;
  string texture = "";
  pbrt_textured1f() : value{0}, texture{} {}
  pbrt_textured1f(float v) : value{v}, texture{} {}
};
struct pbrt_textured3f {
  pbrt_spectrum3f value   = {0, 0, 0};
  string          texture = "";
  pbrt_textured3f() : value{0, 0, 0}, texture{} {}
  pbrt_textured3f(float x, float y, float z) : value{x, y, z}, texture{} {}
};

// pbrt textures
struct pbrt_texture {
  struct constant_t {
    pbrt_textured3f value = {1, 1, 1};
  };
  struct bilerp_t {
    pbrt_textured3f v00 = {0, 0, 0};
    pbrt_textured3f v01 = {1, 1, 1};
    pbrt_textured3f v10 = {0, 0, 0};
    pbrt_textured3f v11 = {1, 1, 1};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
  };
  struct checkerboard_t {
    enum struct aamode_type { closedform, none };
    int             dimension = 2;
    pbrt_textured3f tex1      = {1, 1, 1};
    pbrt_textured3f tex2      = {0, 0, 0};
    aamode_type     aamode    = aamode_type::closedform;
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
  };
  struct dots_t {
    pbrt_textured3f inside  = {1, 1, 1};
    pbrt_textured3f outside = {0, 0, 0};
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
  };
  struct fbm_t {
    int   octaves   = 8;
    float roughness = 0.5;
  };
  struct imagemap_t {
    enum wrap_type { repeat, black, clamp };
    string    filename      = "";
    wrap_type wrap          = wrap_type::repeat;
    float     maxanisotropy = 8;
    bool      trilinear     = false;
    float     scale         = 1;
    bool      gamma         = true;
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
  };
  struct marble_t {
    int   octaves   = 8;
    float roughness = 0.5;
    float scale     = 1;
    float variation = 0.2;
  };
  struct mix_t {
    pbrt_textured3f tex1   = {1, 1, 1};
    pbrt_textured3f tex2   = {1, 1, 1};
    pbrt_textured1f amount = 0.5;
  };
  struct scale_t {
    pbrt_textured3f tex1 = {1, 1, 1};
    pbrt_textured3f tex2 = {1, 1, 1};
  };
  struct uv_t {
    enum struct mapping_type { uv, spherical, cylindrical, planar };
    mapping_type mapping = mapping_type::uv;
    float        uscale  = 1;
    float        vscale  = 1;
    float        udelta  = 0;
    float        vdelta  = 0;
    vec3f        v1      = {1, 0, 0};
    vec3f        v2      = {0, 1, 0};
  };
  struct windy_t {
    // TODO: missing parameters
  };
  struct wrinkled_t {
    int   octaves   = 8;
    float roughness = 0.5;
  };
  enum struct type_t {
    constant,
    bilerp,
    checkerboard,
    dots,
    fbm,
    imagemap,
    marble,
    mix,
    scale,
    uv,
    windy,
    wrinkled
  };
  type_t         type         = type_t::constant;
  constant_t     constant     = {};
  bilerp_t       bilerp       = {};
  checkerboard_t checkerboard = {};
  dots_t         dots         = {};
  fbm_t          fbm          = {};
  imagemap_t     imagemap     = {};
  marble_t       marble       = {};
  mix_t          mix          = {};
  scale_t        scale        = {};
  uv_t           uv           = {};
  windy_t        windy        = {};
  wrinkled_t     wrinkled     = {};
};

// pbrt materials
struct pbrt_material {
  struct matte_t {
    pbrt_textured3f Kd      = {0.5, 0.5, 0.5};
    pbrt_textured1f sigma   = 0;
    pbrt_textured1f bumpmap = 0;
  };
  struct mirror_t {
    pbrt_textured3f Kr      = {0.9, 0.9, 0.9};
    pbrt_textured1f bumpmap = 0;
  };
  struct plastic_t {
    pbrt_textured3f Kd             = {0.25, 0.25, 0.25};
    pbrt_textured3f Ks             = {0.25, 0.25, 0.25};
    pbrt_textured1f uroughness     = 0.1;
    pbrt_textured1f vroughness     = 0.1;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct metal_t {
    pbrt_textured3f eta        = {0.2004376970f, 0.9240334304f, 1.1022119527f};
    pbrt_textured3f k          = {3.9129485033f, 2.4528477015f, 2.1421879552f};
    pbrt_textured1f uroughness = 0.01;
    pbrt_textured1f vroughness = 0.01;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct glass_t {
    pbrt_textured3f Kr             = {1, 1, 1};
    pbrt_textured3f Kt             = {1, 1, 1};
    pbrt_textured1f eta            = 1.5;
    pbrt_textured1f uroughness     = 0;
    pbrt_textured1f vroughness     = 0;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct translucent_t {
    pbrt_textured3f Kd             = {0.25, 0.25, 0.25};
    pbrt_textured3f Ks             = {0.25, 0.25, 0.25};
    pbrt_textured3f reflect        = {0.5, 0.5, 0.5};
    pbrt_textured3f transmit       = {0.5, 0.5, 0.5};
    pbrt_textured1f uroughness     = 0.1;
    pbrt_textured1f vroughness     = 0.1;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct uber_t {
    pbrt_textured3f Kd             = {0.25, 0.25, 0.25};
    pbrt_textured3f Ks             = {0.25, 0.25, 0.25};
    pbrt_textured3f Kr             = {0, 0, 0};
    pbrt_textured3f Kt             = {0, 0, 0};
    pbrt_textured1f uroughness     = 0.1;
    pbrt_textured1f vroughness     = 0.1;
    pbrt_textured1f eta            = 1.5;
    pbrt_textured3f opacity        = {1, 1, 1};
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct disney_t {
    pbrt_textured3f color           = {0.5, 0.5, 0.5};
    pbrt_textured1f anisotropic     = 0;
    pbrt_textured1f clearcoat       = 0;
    pbrt_textured1f clearcoatgloss  = 1;
    pbrt_textured1f eta             = 1.5;
    pbrt_textured1f metallic        = 0;
    pbrt_textured1f uroughness      = 0.5;
    pbrt_textured1f vroughness      = 0.5;
    pbrt_textured3f scatterdistance = {0, 0, 0};
    pbrt_textured1f sheen           = 0;
    pbrt_textured1f sheentint       = 0.5;
    pbrt_textured1f spectrans       = 0;
    pbrt_textured1f speculartint    = 0;
    bool            thin            = false;
    pbrt_textured3f difftrans       = {1, 1, 1};
    pbrt_textured3f flatness        = {0, 0, 0};
    bool            remaproughness  = true;
    pbrt_textured1f bumpmap         = 0;
  };
  struct fourier_t {
    string          bsdffile = "";
    pbrt_textured1f bumpmap  = 0;
    enum struct approx_type_t { plastic, metal, glass };
    approx_type_t approx_type    = approx_type_t::plastic;
    plastic_t     approx_plastic = {};
    metal_t       approx_metal   = {};
    glass_t       approx_glass   = {};
  };
  struct hair_t {
    pbrt_textured3f color       = {0, 0, 0};  // TODO: missing default
    pbrt_textured3f sigma_a     = {0, 0, 0};  // TODO: missing default
    pbrt_textured1f eumelanin   = 0;          // TODO: missing default
    pbrt_textured1f pheomelanin = 0;          // TODO: missing default
    pbrt_textured1f eta         = 1.55f;
    pbrt_textured1f beta_m      = 0.3f;
    pbrt_textured1f beta_n      = 0.3f;
    pbrt_textured1f alpha       = 2;
    pbrt_textured1f bumpmap     = 0;
  };
  struct kdsubsurface_t {
    pbrt_textured3f Kd             = {0.5, 0.5, 0.5};
    pbrt_textured3f mfp            = {1, 1, 1};
    pbrt_textured1f eta            = 1.3;
    pbrt_textured3f Kr             = {1, 1, 1};
    pbrt_textured3f Kt             = {1, 1, 1};
    pbrt_textured1f uroughness     = 0;
    pbrt_textured1f vroughness     = 0;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct mix_t {
    pbrt_textured3f amount         = {0, 0, 0};
    string          namedmaterial1 = "";
    string          namedmaterial2 = "";
    pbrt_textured1f bumpmap        = 0;
  };
  struct substrate_t {
    pbrt_textured3f Kd             = {0.5, 0.5, 0.5};
    pbrt_textured3f Ks             = {0.5, 0.5, 0.5};
    pbrt_textured1f uroughness     = 0.1;
    pbrt_textured1f vroughness     = 0.1;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  struct subsurface_t {
    string          name           = "";
    pbrt_textured3f sigma_a        = {.0011, .0024, .014};
    pbrt_textured3f sigma_prime_s  = {2.55, 3.12, 3.77};
    float           scale          = 1;
    pbrt_textured1f eta            = 1;
    pbrt_textured3f Kr             = {1, 1, 1};
    pbrt_textured3f Kt             = {1, 1, 1};
    pbrt_textured1f uroughness     = 0;
    pbrt_textured1f vroughness     = 0;
    bool            remaproughness = true;
    pbrt_textured1f bumpmap        = 0;
  };
  enum struct type_t {
    matte,
    mirror,
    plastic,
    metal,
    glass,
    translucent,
    uber,
    disney,
    fourier,
    hair,
    kdsubsurface,
    mix,
    substrate,
    subsurface
  };
  type_t         type         = type_t::matte;
  matte_t        matte        = {};
  mirror_t       mirror       = {};
  plastic_t      plastic      = {};
  metal_t        metal        = {};
  glass_t        glass        = {};
  translucent_t  translucent  = {};
  uber_t         uber         = {};
  disney_t       disney       = {};
  fourier_t      fourier      = {};
  hair_t         hair         = {};
  kdsubsurface_t kdsubsurface = {};
  mix_t          mix          = {};
  substrate_t    substrate    = {};
  subsurface_t   subsurface{};
};

// pbrt shapes
struct pbrt_shape {
  struct trianglemesh_t {
    vector<vec3i>   indices     = {};
    vector<vec3f>   P           = {};
    vector<vec3f>   N           = {};
    vector<vec3f>   S           = {};
    vector<vec2f>   uv          = {};
    pbrt_textured1f alpha       = 1;
    pbrt_textured1f shadowalpha = 1;
  };
  struct plymesh_t {
    string          filename    = {};
    pbrt_textured1f alpha       = 1;
    pbrt_textured1f shadowalpha = 1;
  };
  struct curve_t {
    enum struct type_t { flat, ribbon, cylinder };
    enum struct basis_t { bezier, bspline };
    vector<vec3f> P          = {};
    basis_t       basis      = basis_t::bezier;
    int           degree     = 3;
    type_t        type       = type_t::flat;
    vector<vec3f> N          = {};
    float         width0     = 1;
    float         width1     = 1;
    int           splitdepth = 3;
  };
  struct loopsubdiv_t {
    int           levels  = 3;
    vector<vec3i> indices = {};
    vector<vec3f> P       = {};
  };
  struct nurbs_t {
    int           nu     = -1;
    int           nv     = -1;
    vector<float> uknots = {};
    vector<float> vknots = {};
    float         u0     = -1;
    float         v0     = -1;
    float         u1     = -1;
    float         v1     = -1;
    vector<vec3f> P      = {};
    vector<float> Pw     = {};
  };
  struct sphere_t {
    float radius = 1;
    float zmin   = -radius;
    float zmax   = radius;
    float phimax = 360;
  };
  struct disk_t {
    float height      = 0;
    float radius      = 1;
    float innerradius = 0;
    float phimax      = 360;
  };
  struct cone_t {
    float radius = 1;
    float height = 1;
    float phimax = 360;
  };
  struct cylinder_t {
    float radius = 1;
    float zmin   = -1;
    float zmax   = 1;
    float phimax = 360;
  };
  struct hyperboloid_t {
    vec3f p1     = {0, 0, 0};
    vec3f p2     = {1, 1, 1};
    float phimax = 360;
  };
  struct paraboloid_t {
    float radius = 1;
    float zmin   = 0;
    float zmax   = 1;
    float phimax = 360;
  };
  struct heightfield_t {
    int           nu = 0;
    int           nv = 0;
    vector<float> Pz = {};
  };
  enum struct type_t {
    trianglemesh,
    plymesh,
    curve,
    loopsubdiv,
    nurbs,
    sphere,
    disk,
    cone,
    cylinder,
    hyperboloid,
    paraboloid,
    heightfield
  };
  type_t         type         = type_t::trianglemesh;
  trianglemesh_t trianglemesh = {};
  plymesh_t      plymesh      = {};
  curve_t        curve        = {};
  loopsubdiv_t   loopsubdiv   = {};
  nurbs_t        nurbs        = {};
  sphere_t       sphere       = {};
  disk_t         disk         = {};
  cone_t         cone         = {};
  cylinder_t     cylinder     = {};
  hyperboloid_t  hyperboloid  = {};
  paraboloid_t   paraboloid   = {};
  heightfield_t  heightfield  = {};
};

// pbrt lights
struct pbrt_light {
  struct distant_t {
    pbrt_spectrum3f scale = {1, 1, 1};
    pbrt_spectrum3f L     = {1, 1, 1};
    vec3f           from  = {0, 0, 0};
    vec3f           to    = {0, 0, 1};
  };
  struct goniometric_t {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f I       = {1, 1, 1};
    string          mapname = "";
  };
  struct infinite_t {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f L       = {1, 1, 1};
    int             samples = 1;
    string          mapname = "";
  };
  struct point_t {
    pbrt_spectrum3f scale = {1, 1, 1};
    pbrt_spectrum3f I     = {1, 1, 1};
    vec3f           from  = {0, 0, 0};
  };
  struct projection_t {
    pbrt_spectrum3f scale   = {1, 1, 1};
    pbrt_spectrum3f I       = {1, 1, 1};
    float           fov     = 45;
    string          mapname = "";
  };
  struct spot_t {
    pbrt_spectrum3f scale          = {1, 1, 1};
    pbrt_spectrum3f I              = {1, 1, 1};
    vec3f           from           = {0, 0, 0};
    vec3f           to             = {0, 0, 1};
    float           coneangle      = 30;
    float           conedeltaangle = 5;
  };
  enum struct type_t {
    distant,
    goniometric,
    infinite,
    point,
    projection,
    spot
  };
  type_t        type        = type_t::distant;
  distant_t     distant     = {};
  goniometric_t goniometric = {};
  infinite_t    infinite    = {};
  point_t       point       = {};
  projection_t  projection  = {};
  spot_t        spot        = {};
};

// pbrt area lights
struct pbrt_arealight {
  struct none_t {};
  struct diffuse_t {
    pbrt_spectrum3f scale    = {1, 1, 1};
    pbrt_spectrum3f L        = {1, 1, 1};
    bool            twosided = false;
    int             samples  = 1;
  };
  enum struct type_t { none, diffuse };
  type_t    type    = type_t::none;
  none_t    none    = {};
  diffuse_t diffuse = {};
};

// pbrt mediums
struct pbrt_medium {
  struct homogeneous_t {
    pbrt_spectrum3f sigma_a = {0.0011, 0.0024, 0.014};
    pbrt_spectrum3f sigma_s = {2.55, 3.21, 3.77};
    string          preset  = "";
    float           g       = 0;
    float           scale   = 1;
  };
  struct heterogeneous_t {
    pbrt_spectrum3f sigma_a = {0.0011f, 0.0024f, 0.014f};
    pbrt_spectrum3f sigma_s = {2.55f, 3.21f, 3.77f};
    string          preset  = "";
    float           g       = 0;
    float           scale   = 1;
    vec3f           p0      = {0, 0, 0};
    vec3f           p1      = {1, 1, 1};
    int             nx      = 1;
    int             ny      = 1;
    int             nz      = 1;
    vector<float>   density = {};
  };
  enum struct type_t { homogeneous, heterogeneous };
  type_t          type          = type_t::homogeneous;
  homogeneous_t   homogeneous   = {};
  heterogeneous_t heterogeneous = {};
};

// pbrt medium interface
struct pbrt_mediuminterface {
  string interior = "";
  string exterior = "";
};

// pbrt insstance
struct pbrt_object {
  string name = "";
};

// pbrt include
struct pbrt_include {
  string path = "";
};

// pbrt stack ctm
struct pbrt_context {
  frame3f transform_start        = identity3x4f;
  frame3f transform_end          = identity3x4f;
  string  material               = "";
  string  arealight              = "";
  string  medium_interior        = "";
  string  medium_exterior        = "";
  bool    reverse                = false;
  bool    active_transform_start = true;
  bool    active_transform_end   = true;
  float   last_lookat_distance   = 0;
};

// pbrt callbacks
struct pbrt_callbacks {
  virtual void sampler(const pbrt_sampler& value, const pbrt_context& ctx) {}
  virtual void integrator(
      const pbrt_integrator& value, const pbrt_context& ctx) {}
  virtual void accelerator(
      const pbrt_accelerator& value, const pbrt_context& ctx) {}
  virtual void film(const pbrt_film& value, const pbrt_context& ctx) {}
  virtual void filter(const pbrt_filter& value, const pbrt_context& ctx) {}
  virtual void camera(const pbrt_camera& value, const pbrt_context& ctx) {}
  virtual void texture(
      const pbrt_texture& value, const string& name, const pbrt_context& ctx) {}
  virtual void material(const pbrt_material& value, const string& name,
      const pbrt_context& ctx) {}
  virtual void medium(
      const pbrt_medium& value, const string& name, const pbrt_context& ctx) {}
  virtual void shape(const pbrt_shape& value, const pbrt_context& ctx) {}
  virtual void light(const pbrt_light& value, const pbrt_context& ctx) {}
  virtual void arealight(const pbrt_arealight& value, const string& name,
      const pbrt_context& ctx) {}
  virtual void object_instance(
      const pbrt_object& value, const pbrt_context& ctx) {}
  virtual void begin_object(const pbrt_object& value, const pbrt_context& ctx) {
  }
  virtual void end_object(const pbrt_object& value, const pbrt_context& ctx) {}
};

// Load pbrt scene
void load_pbrt(const string& filename, pbrt_callbacks& cb, bool flipv = true);

// Pbrt element
enum struct pbrt_element {
  // clang-format off
  sampler, integrator, accelerator, film, filter, camera, shape, light, // value
  texture, material, medium, arealight, // name and value
  object_instance, begin_object, end_object, include // name
  // clang-format on
};

// Pbrt element data
struct pbrt_element_data {
  pbrt_sampler     sampler;
  pbrt_integrator  intergrator;
  pbrt_accelerator accelerator;
  pbrt_film        film;
  pbrt_filter      filter;
  pbrt_camera      camera;
  pbrt_texture     texture;
  pbrt_material    material;
  pbrt_medium      medium;
  pbrt_shape       shape;
  pbrt_light       light;
  pbrt_arealight   arealight;
};

// Pbrt parser state. Used only internally.
struct pbrt_parser_state {
  unordered_map<string, pair<frame3f, frame3f>> coordsys        = {};
  unordered_map<string, pbrt_spectrum3f>        constant_values = {};
  string                                        object          = "";
  string                                        line            = "";
};

// Read a pbrt element
bool read_pbrt_element(file_wrapper& fs, pbrt_element& element, string& name,
    pbrt_element_data& data, vector<pbrt_context>& stack,
    pbrt_parser_state& state);

}  // namespace yocto

#endif
