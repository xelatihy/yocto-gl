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
// 4. load and save images with `load_image()` and `save_image()`
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
// 2. read a value with `value = parse_argument(parser, name, default, help)`
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

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Normalize path delimiters.
string normalize_path(const string& filename);
// Get directory name (not including '/').
string get_dirname(const string& filename);
// Get extension (not including '.').
string get_extension(const string& filename);
// Get filename without directory.
string get_filename(const string& filename);
// Replace extension.
string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
bool exists_file(const string& filename);

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRIVIAL COMMAND LINE PARSING
// -----------------------------------------------------------------------------
namespace ygl {

// Forward declaration
struct cmdline_parser;

// Initialize a command line parser.
void init_cmdline_parser(cmdline_parser& parser, int argc, char** argv,
    const string& usage, const string& cmd = "");
// check if any error occurred and exit appropriately
void check_cmdline(cmdline_parser& parser);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line). Booleans are flags.
template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, bool req = false);
// Parse all arguments left on the command line.
template <typename T>
inline vector<T> parse_arguments(cmdline_parser& parser, const string& name,
    const vector<T>& def, const string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, const vector<string>& labels, bool req = false);

// Parse an int, float, string, vecXX and bool option or positional argument.
// Options's names starts with "--" or "-", otherwise they are arguments.
// vecXX options use space-separated values but all in one argument
// (use " or ' from the common line). Booleans are flags.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& val, const string& usage, bool req = false);
// Parse all arguments left on the command line.
template <typename T>
inline bool parse_arguments_ref(cmdline_parser& parser, const string& name,
    vector<T>& val, const string& usage, bool req = false);
// Parse a labeled enum, with enum values that are successive integers.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& val, const string& usage, const vector<string>& labels, bool req = false);

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Load/save a text file
bool load_text(const string& filename, string& str);
bool save_text(const string& filename, const string& str);

// Load/save a binary file
bool load_binary(const string& filename, vector<byte>& data);
bool save_binary(const string& filename, const vector<byte>& data);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Check if an image is HDR based on filename.
bool is_hdr_filename(const string& filename);

// Loads/saves a 4 channel float image in linear color space.
bool load_image(const string& filename, image<vec4f>& img);
bool save_image(const string& filename, const image<vec4f>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4f>& img);

// Loads/saves a 4 channel byte image in sRGB color space.
bool load_image(const string& filename, image<vec4b>& img);
bool save_image(const string& filename, const image<vec4b>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4b>& img);
bool load_image_from_memory(const byte* data, int data_size, image<vec4b>& img);

// Convenience helper that saves an HDR images as wither a linear HDR file or
// a tonemapped LDR file depending on file name
bool save_tonemapped_image(const string& filename, const image<vec4f>& hdr,
    float exposure = 0, bool filmic = false, bool srgb = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// VOLUME IMAGE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Loads/saves a 1 channel volume.
bool load_volume1f(const string& filename, volume<float>& vol);
bool save_volume1f(const string& filename, const volume<float>& vol);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Scene load options
struct load_scene_options {
    bool          skip_textures             = false;
    bool          exit_on_error             = false;
    bool          obj_split_shapes          = true;
    bool          obj_preserve_face_varying = false;
    bool          assign_texture_opacity    = true;
    atomic<bool>* cancel_flag               = nullptr;
    bool          run_serially              = false;
};
// Scene save options
struct save_scene_options {
    bool          skip_textures = false;
    bool          exit_on_error = false;
    atomic<bool>* cancel_flag   = nullptr;
    bool          run_serially  = false;
};

// Load/save a scene in the supported formats.
bool load_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene in the builtin JSON format.
bool load_json_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_json_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to OBJ.
bool load_obj_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_obj_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to glTF.
bool load_gltf_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_gltf_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a scene from/to pbrt. This is not robust at all and only
// works on scene that have been previously adapted since the two renderers
// are too different to match.
bool load_pbrt_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_pbrt_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

// Load/save a binary dump useful for very fast scene IO. This format is not
// an archival format and should only be used as an intermediate format.
bool load_ybin_scene(const string& filename, yocto_scene& scene,
    const load_scene_options& options = {});
bool save_ybin_scene(const string& filename, const yocto_scene& scene,
    const save_scene_options& options = {});

}  // namespace ygl

// -----------------------------------------------------------------------------
// MESH IO FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Load/Save a mesh
bool load_mesh(const string& filename, vector<int>& points, vector<vec2i>& lines,
    vector<vec3i>& triangles, vector<vec4i>& quads, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords, vector<vec4f>& colors,
    vector<float>& radius, bool force_triangles);
bool save_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii = false);

// Load/Save a ply mesh
bool load_ply_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<vec4f>& color, vector<float>& radius,
    bool force_triangles);
bool save_ply_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<vec4f>& colors, const vector<float>& radius, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_mesh(const string& filename, vector<int>& points,
    vector<vec2i>& lines, vector<vec3i>& triangles, vector<vec4i>& quads,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, bool force_triangles,
    bool flip_texcoord = true);
bool save_obj_mesh(const string& filename, const vector<int>& points,
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    bool flip_texcoord = true);

// Load/Save a face-varying mesh
bool load_facevarying_mesh(const string& filename, vector<vec4i>& quads_positions,
    vector<vec4i>& quads_normals, vector<vec4i>& quads_texturecoords,
    vector<vec3f>& positions, vector<vec3f>& normals,
    vector<vec2f>& texturecoords, vector<int>& quads_materials);
bool save_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool ascii = false);

// Load/Save an OBJ mesh
bool load_obj_facevarying_mesh(const string& filename,
    vector<vec4i>& quads_positions, vector<vec4i>& quads_normals,
    vector<vec4i>& quads_texturecoords, vector<vec3f>& positions,
    vector<vec3f>& normals, vector<vec2f>& texturecoords,
    vector<int>& quads_materials, bool flip_texcoord = true);
bool save_obj_facevarying_mesh(const string& filename,
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    const vector<int>& quads_materials, bool flip_texcoord = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER
// -----------------------------------------------------------------------------
namespace ygl {

// OBJ vertex
struct obj_vertex {
    int position     = 0;
    int texturecoord = 0;
    int normal       = 0;
};

// Obj texture information.
struct obj_texture_info {
    string path  = "";     // file path
    bool   clamp = false;  // clamp to edge
    float  scale = 1;      // scale for bump/displacement
    // Properties not explicitly handled.
    unordered_map<string, vector<float>> props;
};

// Obj material.
struct obj_material {
    string name;       // name
    int    illum = 0;  // MTL illum mode

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
    unordered_map<string, vector<string>> props;
};

// Obj camera [extension].
struct obj_camera {
    string  name;                         // name
    frame3f frame    = identity_frame3f;  // transform
    bool    ortho    = false;             // orthographic
    vec2f   film     = {0.036f, 0.024f};  // film size (default to 35mm)
    float   focal    = 0.050f;            // focal length
    float   aspect   = 16.0f / 9.0f;      // aspect ratio
    float   aperture = 0;                 // lens aperture
    float   focus    = maxf;              // focus distance
};

// Obj environment [extension].
struct obj_environment {
    string           name;                      // name
    frame3f          frame = identity_frame3f;  // transform
    vec3f            ke    = zero3f;            // emission color
    obj_texture_info ke_txt;                    // emission texture
};

// Obj callbacks
struct obj_callbacks {
    function<void(const vec3f&)>              vert        = {};
    function<void(const vec3f&)>              norm        = {};
    function<void(const vec2f&)>              texcoord    = {};
    function<void(const vector<obj_vertex>&)> face        = {};
    function<void(const vector<obj_vertex>&)> line        = {};
    function<void(const vector<obj_vertex>&)> point       = {};
    function<void(const string& name)>        object      = {};
    function<void(const string& name)>        group       = {};
    function<void(const string& name)>        usemtl      = {};
    function<void(const string& name)>        smoothing   = {};
    function<void(const string& name)>        mtllib      = {};
    function<void(const obj_material&)>       material    = {};
    function<void(const obj_camera&)>         camera      = {};
    function<void(const obj_environment&)>    environmnet = {};
};

// Load obj options
struct load_obj_options {
    bool exit_on_error = false;
    bool geometry_only = false;
    bool flip_texcoord = true;
    bool flip_tr       = true;
};

// Load obj scene
bool load_obj(const string& filename, const obj_callbacks& cb,
    const load_obj_options& options = {});

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE PLY LOADER
// -----------------------------------------------------------------------------
namespace ygl {

// ply type
enum struct ply_type { ply_uchar, ply_int, ply_float, ply_int_list };

// ply property
struct ply_property {
    string                name    = "";
    ply_type              type    = ply_type::ply_float;
    vector<float>         scalars = {};
    vector<array<int, 8>> lists   = {};
};

// ply element
struct ply_element {
    string               name       = "";
    int                  count      = 0;
    vector<ply_property> properties = {};
};

// simple ply api data
struct ply_data {
    vector<ply_element> elements = {};
};

// Load ply mesh
bool load_ply(const string& filename, ply_data& ply);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF COMMAND-LINE PARSING
// -----------------------------------------------------------------------------
namespace ygl {

// Command line parser data. All data should be considered private.
struct cmdline_parser {
    vector<string> args      = {};  // command line arguments
    string         usage_cmd = "";  // program name
    string         usage_hlp = "";  // program help
    string         usage_opt = "";  // options help
    string         usage_arg = "";  // arguments help
    string         error     = "";  // current parse error
};

// initialize a command line parser
inline void init_cmdline_parser(cmdline_parser& parser, int argc, char** argv,
    const string& usage, const string& cmd) {
    parser           = cmdline_parser{};
    parser.args      = {argv + 1, argv + argc};
    parser.usage_cmd = (cmd.empty()) ? argv[0] : cmd;
    parser.usage_hlp = usage;
}

// check if option or argument
inline bool is_option(const string& name) {
    return name.size() > 1 && name.front() == '-';
}

// get names from string
inline vector<string> get_option_names(const string& name_) {
    auto names = vector<string>();
    auto name  = name_;
    while (name.find(',') != name.npos) {
        names.push_back(name.substr(0, name.find(',')));
        name = name.substr(name.find(',') + 1);
    }
    names.push_back(name);
    return names;
}

// add help
inline string get_option_usage(const string& name, const string& usage,
    const string& def_, const vector<string>& choices) {
    auto def = def_;
    if (def != "") def = "[" + def + "]";
    auto namevar = name;
    if (name != "") namevar += " " + name;
    char buffer[4096];
    sprintf(
        buffer, "  %-24s %s %s\n", namevar.c_str(), usage.c_str(), def.c_str());
    auto usagelines = string(buffer);
    if (!choices.empty()) {
        usagelines += "        accepted values:";
        for (auto& c : choices) usagelines += " " + c;
        usagelines += "\n";
    }
    return usagelines;
}

// print cmdline help
inline void print_cmdline_usage(const cmdline_parser& parser) {
    printf("%s: %s\n", parser.usage_cmd.c_str(), parser.usage_hlp.c_str());
    printf("usage: %s %s %s\n\n", parser.usage_cmd.c_str(),
        (parser.usage_opt.empty()) ? "" : "[options]",
        (parser.usage_arg.empty()) ? "" : "arguments");
    if (!parser.usage_opt.empty()) {
        printf("options:\n");
        printf("%s\n", parser.usage_opt.c_str());
    }
    if (!parser.usage_arg.empty()) {
        printf("arguments:\n");
        printf("%s\n", parser.usage_arg.c_str());
    }
}

// Parse a flag. Name should start with either "--" or "-".
inline bool parse_flag_argument(cmdline_parser& parser, const string& name,
    bool& value, const string& usage);

// check if any error occurred and exit appropriately
inline void check_cmdline(cmdline_parser& parser) {
    auto help = false;
    if (parse_flag_argument(parser, "--help,-?", help, "print help")) {
        print_cmdline_usage(parser);
        exit(0);
    }
    if (!parser.args.empty()) parser.error += "unmatched arguments remaining\n";
    if (!parser.error.empty()) {
        printf("error: %s", parser.error.c_str());
        print_cmdline_usage(parser);
        exit(1);
    }
}

// Parse an option string. Name should start with "--" or "-".
template <typename T>
inline bool parse_option_argument(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req, const vector<string>& choices) {
    parser.usage_opt += get_option_usage(name, usage, to_string(value), choices);
    if (parser.error != "") return false;
    auto names = get_option_names(name);
    auto pos   = parser.args.end();
    for (auto& name : names) {
        pos = std::min(
            pos, std::find(parser.args.begin(), parser.args.end(), name));
    }
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name;
        return false;
    }
    if (pos == parser.args.end() - 1) {
        parser.error += "missing value for " + name;
        return false;
    }
    auto vals = *(pos + 1);
    parser.args.erase(pos, pos + 2);
    if (!choices.empty() &&
        std::find(choices.begin(), choices.end(), vals) == choices.end()) {
        parser.error += "bad value for " + name;
        return false;
    }
    auto new_value = value;
    if (!parse(vals, new_value)) {
        parser.error += "bad value for " + name;
        return false;
    }
    value = new_value;
    return true;
}

// Parse an argument string. Name should not start with "--" or "-".
template <typename T>
inline bool parse_positional_argument(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req, const vector<string>& choices) {
    parser.usage_arg += get_option_usage(name, usage, to_string(value), choices);
    if (parser.error != "") return false;
    auto pos = std::find_if(parser.args.begin(), parser.args.end(),
        [](auto& v) { return v[0] != '-'; });
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name;
        return false;
    }
    auto vals = *pos;
    parser.args.erase(pos);
    if (!choices.empty() &&
        std::find(choices.begin(), choices.end(), vals) == choices.end()) {
        parser.error += "bad value for " + name;
        return false;
    }
    auto new_value = value;
    if (!parse(vals, new_value)) {
        parser.error += "bad value for " + name;
        return false;
    }
    value = new_value;
    return true;
}

// Parse all left argument strings. Name should not start with "--" or "-".
template <typename T>
inline bool parse_positional_arguments(cmdline_parser& parser,
    const string& name, vector<T>& values, const string& usage, bool req) {
    auto defs = string();
    for (auto& d : values) defs += " " + d;
    parser.usage_arg += get_option_usage(name, usage, defs, {});
    if (parser.error != "") return false;
    auto pos = std::find_if(parser.args.begin(), parser.args.end(),
        [](auto& v) { return v[0] != '-'; });
    if (pos == parser.args.end()) {
        if (req) parser.error += "missing value for " + name;
        return false;
    }
    auto vals = vector<string>{pos, parser.args.end()};
    parser.args.erase(pos, parser.args.end());
    auto new_values = values;
    new_values.resize(vals.size());
    for (auto i = 0; i < vals.size(); i++) {
        if (!parse(vals[i], new_values[i])) {
            parser.error += "bad value for " + name;
            return false;
        }
    }
    values = new_values;
    return true;
}

// Parse a flag. Name should start with either "--" or "-".
inline bool parse_flag_argument(cmdline_parser& parser, const string& name,
    bool& value, const string& usage) {
    parser.usage_opt += get_option_usage(name, usage, "", {});
    if (parser.error != "") return false;
    auto names = get_option_names(name);
    auto pos   = parser.args.end();
    for (auto& name : names)
        pos = std::min(
            pos, std::find(parser.args.begin(), parser.args.end(), name));
    if (pos == parser.args.end()) return false;
    parser.args.erase(pos);
    value = !value;
    return true;
}

// Parse an integer, float, string. If name starts with "--" or "-", then it is
// an option, otherwise it is a position argument.
template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& value, const string& usage, bool req) {
    return is_option(name) ?
               parse_option_argument(parser, name, value, usage, req, {}) :
               parse_positional_argument(parser, name, value, usage, req, {});
}
template <>
inline bool parse_argument_ref<bool>(cmdline_parser& parser, const string& name,
    bool& value, const string& usage, bool req) {
    return parse_flag_argument(parser, name, value, usage);
}

template <typename T>
inline bool parse_argument_ref(cmdline_parser& parser, const string& name,
    T& value, const string& usage, const vector<string>& labels, bool req) {
    auto values = labels.at((int)value);
    auto parsed = is_option(name) ? parse_option_argument(parser, name, values,
                                        usage, req, labels) :
                                    parse_positional_argument(parser, name,
                                        values, usage, req, labels);
    if (!parsed) return false;
    auto pos = std::find(labels.begin(), labels.end(), values);
    if (pos == labels.end()) return false;
    value = (T)(pos - labels.begin());
    return true;
}

// Parser an argument
template <typename T>
inline bool parse_arguments_ref(cmdline_parser& parser, const string& name,
    vector<T>& values, const string& usage, bool req) {
    return parse_positional_arguments(parser, name, values, usage, req);
}

// Parse an integer, float, string. If name starts with "--" or "-", then it is
// an option, otherwise it is a position argument.
template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, bool req) {
    auto value = def;
    if (!parse_argument_ref(parser, name, value, usage, req)) return def;
    return value;
}
template <>
inline bool parse_argument<bool>(cmdline_parser& parser, const string& name,
    bool def, const string& usage, bool req) {
    auto value = def;
    if (!parse_flag_argument(parser, name, value, usage)) return def;
    return value;
}

template <typename T>
inline T parse_argument(cmdline_parser& parser, const string& name, T def,
    const string& usage, const vector<string>& labels, bool req) {
    auto value = def;
    if (!parse_argument_ref(parser, name, value, usage, labels, req))
        return def;
    return value;
}

// Parser an argument
template <typename T>
inline vector<T> parse_arguments(cmdline_parser& parser, const string& name,
    const vector<T>& def, const string& usage, bool req) {
    auto values = vector<T>{};
    if (!parse_arguments_ref(parser, name, values, usage, req)) return def;
    return values;
}

}  // namespace ygl

#endif
