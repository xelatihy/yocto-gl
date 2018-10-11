//
// # Yocto/GLIO: Inpuot and output functions for Yocto/GL
//
// Yocto/GLIO provides loading and saving functionality for images and scenes
// in Yocto/GL. For images we support PNG, JPG, TGA, HDR, EXR formats. For
// scene we support a simple to use JSON format, PLY, OBJ and glTF.
// The JSON serialization is a straight copy of the in-memory scene data.
// To speed up testing, we also support a binary format that is a dump of
// the current scene. This format should not be use for archival though.
//
// We do not use exception as the API for returning errors, although they might
// be used internally in the implementastion of the methods. In load functions,
// as error is signaled by returning an empty object or a null pointer. In
// save functions, errors are returned with the supplied boolean. In the future,
// we will also provide return types with error codes.
//
// ## File Loading and Saving
//
// 1. manipulate paths withe path utilities
// 2. load and save text files with `load_text()` and `save_text()`
// 3. load and save binary files with `load_binary()` and `save_binary()`
// 4. load and save images with `load_image4f()` and `save_image4f()`
// 5. load a scene with `load_json_scene()` and save it with `save_json_scene()`
// 6. load and save OBJs with `load_obj_scene()` and `save_obj_scene()`
// 7. load and save glTFs with `load_gltf_scene()` and `save_gltf_scene()`
// 8. load and save binary dumps with `load_ybin_scene()`and `save_ybin_scene()`
// 9. if desired, the function `load_scene()` and `save_scene()` will either
//    load using the internal format or convert on the fly using on the
//    supported conversions
//
//
// ## Command-Line Parsing
//
// We also provide a simple, immediate-mode, command-line parser. The parser
// works in an immediate-mode manner since it reads each value as you call each
// function, rather than building a data structure and parsing offline. We
// support option and position arguments, automatic help generation, and
// error checking.
//
// 1. initialize the parser with `make_cmdline_parser(argc, argv, help)`
// 2. read a value with `value = parse_arg(parser, name, default, help)`
//    - is name starts with '--' or '-' then it is an option
//    - otherwise it is a positional arguments
//    - options and arguments may be intermixed
//    - the type of each option is determined by the default value `default`
//    - the value is parsed on the stop
// 3. finished parsing with `check_cmdline(parser)`
//    - if an error occurred, the parser will exit and print a usage message
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#ifndef _YGLIO_H_
#define _YGLIO_H_

#include "ygl.h"

#include <array>
#include <chrono>
#include <cstdio>
#include <cstdlib>

#ifndef YGL_FASTPARSE
#define YGL_FASTPARSE 1
#endif

// -----------------------------------------------------------------------------
// STRING/TIME UTILITIES FOR CLI APPLICATIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline std::string format(const std::string& fmt, const Args&... args);

// Converts to string.
template <typename T>
inline std::string to_string(const T& val);

// Prints a formatted string to stdout or file.
template <typename... Args>
inline bool print(FILE* fs, const std::string& fmt, const Args&... args);
template <typename... Args>
inline bool print(const std::string& fmt, const Args&... args) {
    return print(stdout, fmt, args...);
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const std::string& str, Args&... args);
template <typename... Args>
inline bool parse(FILE* fs, Args&... args);

// Exit to the console with an error.
template <typename... Args>
inline void exit_error(const std::string& fmt, const Args&... args) {
    print(fmt + "\n", args...);
    exit(1);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TIMER UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// STRING FORMAT UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Format duration string from nanoseconds
std::string format_duration(int64_t duration);
// Format a large integer number in human readable form
std::string format_num(uint64_t num);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Normalize path delimiters.
std::string normalize_path(const std::string& filename);
// Get directory name (not including '/').
std::string get_dirname(const std::string& filename);
// Get extension (not including '.').
std::string get_extension(const std::string& filename);
// Get filename without directory.
std::string get_filename(const std::string& filename);
// Replace extension.
std::string replace_extension(
    const std::string& filename, const std::string& ext);

// Check if a file can be opened for reading.
bool exists_file(const std::string& filename);

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRIVIAL COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace ygl {

// Command line parser data. All data should be considered private.
struct cmdline_parser {
    std::vector<std::string> args      = {};  // command line arguments
    std::string              usage_cmd = "";  // program name
    std::string              usage_hlp = "";  // program help
    std::string              usage_opt = "";  // options help
    std::string              usage_arg = "";  // arguments help
    std::string              error     = "";  // current parse error
};

// Initialize a command line parser.
cmdline_parser make_cmdline_parser(int argc, char** argv,
    const std::string& usage, const std::string& cmd = "");
// check if any error occurred and exit appropriately
void check_cmdline(cmdline_parser& parser);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line). Booleans are flags.
bool  parse_arg(cmdline_parser& parser, const std::string& name, bool def,
     const std::string& usage);
int   parse_arg(cmdline_parser& parser, const std::string& name, int def,
      const std::string& usage, bool req = false);
float parse_arg(cmdline_parser& parser, const std::string& name, float def,
    const std::string& usage, bool req = false);
vec2f parse_arg(cmdline_parser& parser, const std::string& name,
    const vec2f& def, const std::string& usage, bool req = false);
vec3f parse_arg(cmdline_parser& parser, const std::string& name,
    const vec3f& def, const std::string& usage, bool req = false);
std::string parse_arg(cmdline_parser& parser, const std::string& name,
    const std::string& def, const std::string& usage, bool req = false);
std::string parse_arg(cmdline_parser& parser, const std::string& name,
    const char* def, const std::string& usage, bool req = false);
// Parse all arguments left on the command line.
std::vector<std::string> parse_args(cmdline_parser& parser,
    const std::string& name, const std::vector<std::string>& def,
    const std::string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
int parse_arge(cmdline_parser& parser, const std::string& name, int def,
    const std::string& usage, const std::vector<std::string>& labels,
    bool req = false);

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Load/save a text file
bool load_text(const std::string& filename, std::string& str);
bool save_text(const std::string& filename, const std::string& str);

// Load/save a binary file
bool load_binary(const std::string& filename, std::vector<byte>& data);
bool save_binary(const std::string& filename, const std::vector<byte>& data);

// Load text/binary file with shorthand call. Returns empty values on error.
std::string       load_text(const std::string& filename);
std::vector<byte> load_binary(const std::string& filename);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

// Loads/saves a 4 channel float image in linear color space.
bool load_image4f(const std::string& filename, image<vec4f>& img);
bool save_image4f(const std::string& filename, const image<vec4f>& img);
bool load_image4f_from_memory(
    const byte* data, int data_size, image<vec4f>& img);

// Loads/saves a 4 channel byte image in sRGB color space.
bool load_image4b(const std::string& filename, image<vec4b>& img);
bool save_image4b(const std::string& filename, const image<vec4b>& img);
bool load_image4b_from_memory(
    const byte* data, int data_size, image<vec4b>& img);

// Load 4 channel images with shortened api. Returns empty image on error;
image<vec4f> load_image4f(const std::string& filename);
image<vec4f> load_image4f_from_memory(const byte* data, int data_size);
image<vec4b> load_image4b(const std::string& filename);
image<vec4b> load_image4b_from_memory(const byte* data, int data_size);

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const std::string& filename, const image<vec4f>& hdr,
    float exposure = 0, bool filmic = false, bool srgb = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Loads/saves a 1 channel volume.
bool load_volume1f(const std::string& filename, volume<float>& vol);
bool save_volume1f(const std::string& filename, const volume<float>& vol);

// Loads a 1 channel volume for shortened api.
volume<float> load_volume1f(const std::string& filename);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/save a scene in the supported formats.
bool load_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true);
bool save_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene in the builtin JSON format.
bool load_json_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true);
bool save_json_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene from/to OBJ.
bool load_obj_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true,
    bool split_shapes = true);
bool save_obj_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene from/to glTF.
bool load_gltf_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true);
bool save_gltf_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
bool load_pbrt_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true);
bool save_pbrt_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load/save a binary dump useful for very fast scene IO. This format is not
// an archival format and should only be used as an intermediate format.
bool load_ybin_scene(const std::string& filename, scene*& scn,
    bool save_textures = true, bool skip_missing = true);
bool save_ybin_scene(const std::string& filename, const scene* scn,
    bool save_textures = true, bool skip_missing = true);

// Load a scene in the supported formats using shortened API.
scene* load_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
scene* load_json_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
scene* load_obj_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true, bool split_shapes = true);
scene* load_gltf_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
scene* load_pbrt_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);
scene* load_ybin_scene(const std::string& filename, bool load_textures = true,
    bool skip_missing = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// MESH IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/Save a mesh
bool load_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<vec4f>& color,
    std::vector<float>& radius);
bool save_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, const std::vector<vec4f>& color,
    const std::vector<float>& radius, bool ascii = false);

// Load/Save a ply mesh
bool load_ply_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<vec4f>& color,
    std::vector<float>& radius);
bool save_ply_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, const std::vector<vec4f>& color,
    const std::vector<float>& radius, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_mesh(const std::string& filename, std::vector<int>& points,
    std::vector<vec2i>& lines, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, bool flip_texcoord = true);
bool save_obj_mesh(const std::string& filename, const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, bool flip_texcoord = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// FACE-VARYING IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/Save a mesh
bool load_fvmesh(const std::string& filename, std::vector<vec4i>& quads_pos,
    std::vector<vec3f>& pos, std::vector<vec4i>& quads_norm,
    std::vector<vec3f>& norm, std::vector<vec4i>& quads_texcoord,
    std::vector<vec2f>& texcoord, std::vector<vec4i>& quads_color,
    std::vector<vec4f>& color);
bool save_fvmesh(const std::string& filename,
    const std::vector<vec4i>& quads_pos, const std::vector<vec3f>& pos,
    const std::vector<vec4i>& quads_norm, const std::vector<vec3f>& norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec2f>& texcoord, const std::vector<vec4i>& quads_color,
    const std::vector<vec4f>& color, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_fvmesh(const std::string& filename, std::vector<vec4i>& quads_pos,
    std::vector<vec3f>& pos, std::vector<vec4i>& quads_norm,
    std::vector<vec3f>& norm, std::vector<vec4i>& quads_texcoord,
    std::vector<vec2f>& texcoord, bool flip_texcoord = true);
bool save_obj_fvmesh(const std::string& filename,
    const std::vector<vec4i>& quads_pos, const std::vector<vec3f>& pos,
    const std::vector<vec4i>& quads_norm, const std::vector<vec3f>& norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec2f>& texcoord, bool flip_texcoord = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER
// -----------------------------------------------------------------------------
namespace ygl {

// OBJ vertex
struct obj_vertex {
    int pos      = 0;
    int texcoord = 0;
    int norm     = 0;
};

// Obj texture information.
struct obj_texture_info {
    std::string path  = "";     // file path
    bool        clamp = false;  // clamp to edge
    float       scale = 1;      // scale for bump/displacement
    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<float>> props;
};

// Obj material.
struct obj_material {
    std::string name;       // name
    int         illum = 0;  // MTL illum mode

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
    float rs  = -1;         // roughness (-1 not defined)
    float km  = -1;         // metallic  (-1 not defined)

    // textures
    obj_texture_info ke_txt;    // emission texture
    obj_texture_info ka_txt;    // ambient texture
    obj_texture_info kd_txt;    // diffuse texture
    obj_texture_info ks_txt;    // specular texture
    obj_texture_info kr_txt;    // reflection texture
    obj_texture_info kt_txt;    // transmission texture
    obj_texture_info ns_txt;    // Phong exponent texture
    obj_texture_info op_txt;    // opacity texture
    obj_texture_info rs_txt;    // roughness texture
    obj_texture_info km_txt;    // metallic texture
    obj_texture_info ior_txt;   // ior texture
    obj_texture_info occ_txt;   // occlusion map
    obj_texture_info bump_txt;  // bump map
    obj_texture_info disp_txt;  // displacement map
    obj_texture_info norm_txt;  // normal map

    // volume data [extension]
    vec3f ve = {0, 0, 0};  // volume emission
    vec3f va = {0, 0, 0};  // albedo: scattering / (absorption + scattering)
    vec3f vd = {0, 0, 0};  // density: absorption + scattering
    float vg = 0;          // phase function shape

    // volume textures [extension]
    obj_texture_info vd_txt;  // density

    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Obj camera [extension].
struct obj_camera {
    std::string name;                         // name
    frame3f     frame    = identity_frame3f;  // transform
    bool        ortho    = false;             // orthographic
    vec2f       film     = {0.036f, 0.024f};  // film size (default to 35mm)
    float       focal    = 0.050f;            // focal length
    float       aspect   = 16.0f / 9.0f;      // aspect ratio
    float       aperture = 0;                 // lens aperture
    float       focus    = maxf;              // focus distance
};

// Obj environment [extension].
struct obj_environment {
    std::string      name;                      // name
    frame3f          frame = identity_frame3f;  // transform
    vec3f            ke    = zero3f;            // emission color
    obj_texture_info ke_txt;                    // emission texture
};

// Obj callbacks
struct obj_callbacks {
    std::function<void(const vec3f&)>                   vert        = {};
    std::function<void(const vec3f&)>                   norm        = {};
    std::function<void(const vec2f&)>                   texcoord    = {};
    std::function<void(const std::vector<obj_vertex>&)> face        = {};
    std::function<void(const std::vector<obj_vertex>&)> line        = {};
    std::function<void(const std::vector<obj_vertex>&)> point       = {};
    std::function<void(const std::string& name)>        object      = {};
    std::function<void(const std::string& name)>        group       = {};
    std::function<void(const std::string& name)>        usemtl      = {};
    std::function<void(const std::string& name)>        smoothing   = {};
    std::function<void(const std::string& name)>        mtllib      = {};
    std::function<void(const obj_material&)>            material    = {};
    std::function<void(const obj_camera&)>              camera      = {};
    std::function<void(const obj_environment&)>         environmnet = {};
};

// Load obj scene
bool load_obj(const std::string& filename, const obj_callbacks& cb,
    bool flip_texcoord = true, bool flip_tr = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER
// -----------------------------------------------------------------------------
namespace ygl {

// ply type
enum struct ply_type { ply_uchar, ply_int, ply_float, ply_int_list };

// ply property
struct ply_property {
    std::string                     name    = "";
    ply_type                        type    = ply_type::ply_float;
    std::vector<float>              scalars = {};
    std::vector<std::array<int, 8>> lists   = {};
};

// ply element
struct ply_element {
    std::string               name       = "";
    int                       count      = 0;
    std::vector<ply_property> properties = {};
};

// simple ply api data
struct ply_data {
    std::vector<ply_element> elements = {};
};

// Load ply mesh
bool load_ply(const std::string& filename, ply_data& ply);

// Load ply mesh with shortened API. Returns empty object on error.
ply_data load_ply(const std::string& filename);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING/TIME UTILITIES FOR CLI APPLICATIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Prints basic types
inline bool print_value(std::string& str, const std::string& val) {
    str += val;
    return true;
}
inline bool print_value(std::string& str, const char* val) {
    str += val;
    return true;
}
inline bool print_value(std::string& str, int val) {
    str += std::to_string(val);
    return true;
}
inline bool print_value(std::string& str, float val) {
    str += std::to_string(val);
    return true;
}
inline bool print_value(std::string& str, double val) {
    str += std::to_string(val);
    return true;
}

// Print compound types.
template <typename T>
inline bool print_value(std::string& str, const vec2<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const vec3<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.z)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const vec4<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.z)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.w)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const mat2<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const mat3<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.z)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const mat4<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.z)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.w)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const frame2<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.o)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const frame3<T>& v) {
    if (!print_value(str, v.x)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.y)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.z)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.o)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const bbox1<T>& v) {
    if (!print_value(str, v.min)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.max)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const bbox2<T>& v) {
    if (!print_value(str, v.min)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.max)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const bbox3<T>& v) {
    if (!print_value(str, v.min)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.max)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const bbox4<T>& v) {
    if (!print_value(str, v.min)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.max)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const ray2<T>& v) {
    if (!print_value(str, v.o)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.d)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.tmin)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.tmax)) return false;
    return true;
}
template <typename T>
inline bool print_value(std::string& str, const ray3<T>& v) {
    if (!print_value(str, v.o)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.d)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.tmin)) return false;
    if (!print_value(str, " ")) return false;
    if (!print_value(str, v.tmax)) return false;
    return true;
}

// Prints a string.
inline bool print_next(std::string& str, const std::string& fmt) {
    return print_value(str, fmt);
}
template <typename Arg, typename... Args>
inline bool print_next(std::string& str, const std::string& fmt, const Arg& arg,
    const Args&... args) {
    auto pos = fmt.find("{}");
    if (pos == std::string::npos) return print_value(str, fmt);
    if (!print_value(str, fmt.substr(0, pos))) return false;
    if (!print_value(str, arg)) return false;
    return print_next(str, fmt.substr(pos + 2), args...);
}

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline std::string format(const std::string& fmt, const Args&... args) {
    auto str = std::string();
    print_next(str, fmt, args...);
    return str;
}

// Prints a string.
template <typename... Args>
inline bool print(FILE* fs, const std::string& fmt, const Args&... args) {
    auto str = format(fmt, args...);
    return fprintf(fs, "%s", str.c_str()) >= 0;
}

// Converts to string.
template <typename T>
inline std::string to_string(const T& val) {
    auto str = std::string();
    print_value(str, val);
    return str;
}

// Prints basic types to string
inline bool _parse(const char*& str, std::string& val) {
    auto n = 0;
    char buf[4096];
    if (sscanf(str, "%4095s%n", buf, &n) != 1) return false;
    val = buf;
    str += n;
    return true;
}
inline bool _parse(const char*& str, int& val) {
    auto n = 0;
    if (sscanf(str, "%d%n", &val, &n) != 1) return false;
    str += n;
    return true;
}
inline bool _parse(const char*& str, float& val) {
    auto n = 0;
    if (sscanf(str, "%f%n", &val, &n) != 1) return false;
    str += n;
    return true;
}
inline bool _parse(const char*& str, double& val) {
    auto n = 0;
    if (sscanf(str, "%lf%n", &val, &n) != 1) return false;
    str += n;
    return true;
}

// Prints basic types to stream
inline bool _parse(FILE* fs, std::string& val) {
    char buf[4096];
    if (fscanf(fs, "%4095s", buf) != 1) return false;
    val = buf;
    return true;
}
inline bool _parse(FILE* fs, int& val) { return fscanf(fs, "%d", &val) == 1; }
inline bool _parse(FILE* fs, float& val) { return fscanf(fs, "%f", &val) == 1; }
inline bool _parse(FILE* fs, double& val) {
    return fscanf(fs, "%lf", &val) == 1;
}

// Print compound types
template <typename Archive, typename T, size_t N>
inline bool _parse(Archive& ar, std::array<T, N>& val) {
    for (auto i = 0; i < N; i++) {
        if (!_parse(ar, val[i])) return false;
    }
    return true;
}
// Data acess
template <typename Archive, typename T>
inline bool _parse(Archive& ar, vec2<T>& v) {
    return _parse(ar, (std::array<T, 2>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, vec3<T>& v) {
    return _parse(ar, (std::array<T, 3>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, vec4<T>& v) {
    return _parse(ar, (std::array<T, 4>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, mat2<T>& v) {
    return _parse(ar, (std::array<T, 4>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, mat3<T>& v) {
    return _parse(ar, (std::array<T, 9>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, mat4<T>& v) {
    return _parse(ar, (std::array<T, 16>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, frame2<T>& v) {
    return _parse(ar, (std::array<T, 6>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, frame3<T>& v) {
    return _parse(ar, (std::array<T, 12>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, bbox1<T>& v) {
    return _parse(ar, (std::array<T, 1>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, bbox2<T>& v) {
    return _parse(ar, (std::array<T, 4>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, bbox3<T>& v) {
    return _parse(ar, (std::array<T, 6>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, bbox4<T>& v) {
    return _parse(ar, (std::array<T, 8>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, ray2<T>& v) {
    return _parse(ar, (std::array<T, 6>&)v);
}
template <typename Archive, typename T>
inline bool _parse(Archive& ar, ray3<T>& v) {
    return _parse(ar, (std::array<T, 8>&)v);
}

// Prints a string.
template <typename Archive>
inline bool _parse_next(Archive& ar) {
    return true;
}
template <typename Archive, typename Arg, typename... Args>
inline bool _parse_next(Archive& ar, Arg& arg, Args&... args) {
    if (!_parse(ar, arg)) return false;
    return _parse_next(ar, args...);
}

// Returns trus if this is white space
inline bool _is_whitespace(const char* str) {
    while (*str == ' ' || *str == '\t' || *str == '\r' || *str == '\n') str++;
    return *str == 0;
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const std::string& str, Args&... args) {
    auto str_ = str.c_str();
    if (!_parse_next(str_, args...)) return false;
    return _is_whitespace(str_);
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(FILE* fs, Args&... args) {
    return _parse_next(fs, args...);
}

}  // namespace ygl

#endif
