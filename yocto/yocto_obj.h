//
// # Yocto/OBJ: Tiny library for OBJ parsing
//
// Yocto/OBJ is a simple Wavefront OBJ parser that works with callbacks.
// We make no attempt to provide a simple interface for OBJ but just the
// low level parsing code. We support a few extensions such as camera and
// environment map loading.
//
// Error reporting is done through exceptions using the `io_error` exception.
//
// ## Parse an OBJ file
//
// 1. define callbacks in `obj_callback` structure using lambda with capture
//    if desired
// 2. run the parse with `load_obj()`
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
#include "yocto_utils.h"

#include <memory>
#include <unordered_map>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::function;
using std::unique_ptr;
using std::unordered_map;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIMPLE OBJ LOADER
// -----------------------------------------------------------------------------
namespace yocto {

// OBJ vertex
struct obj_vertex {
    int position     = 0;
    int texturecoord = 0;
    int normal       = 0;
};

inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.position == b.position && a.texturecoord == b.texturecoord &&
           a.normal == b.normal;
}

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
    float   width    = 0.036f;            // film size (default to 35mm)
    float   height   = 0.024f;            // film size (default to 35mm)
    float   focal    = 0.050f;            // focal length
    float   aspect   = 16.0f / 9.0f;      // aspect ratio
    float   aperture = 0;                 // lens aperture
    float   focus    = float_max;         // focus distance
};

// Obj environment [extension].
struct obj_environment {
    string           name;                      // name
    frame3f          frame = identity_frame3f;  // transform
    vec3f            ke    = zero3f;            // emission color
    obj_texture_info ke_txt;                    // emission texture
};

// Obj procedural object [extension].
struct obj_procedural {
    string  name;                      // name
    frame3f frame = identity_frame3f;  // transform
    string  type;                      // type
    string  material;                  // material
    float   size  = 2;                 // size
    int     level = -1;                // level of subdivision (-1 default)
};

// Obj callbacks
struct obj_callbacks {
    void vert(const vec3f&) {}
    void norm(const vec3f&) {}
    void texcoord(const vec2f&) {}
    void face(const vector<obj_vertex>&) {}
    void line(const vector<obj_vertex>&) {}
    void point(const vector<obj_vertex>&) {}
    void object(const string&) {}
    void group(const string&) {}
    void usemtl(const string&) {}
    void smoothing(const string&) {}
    void mtllib(const string&) {}
    void material(const obj_material&) {}
    void camera(const obj_camera&) {}
    void environmnet(const obj_environment&) {}
    void procedural(const obj_procedural&) {}
};

// Load obj options
struct load_obj_options {
    bool exit_on_error = false;
    bool geometry_only = false;
    bool flip_texcoord = true;
    bool flip_tr       = true;
};

// Load obj scene
template <typename Callbacks>
inline void load_obj(const string& filename, Callbacks& cb,
    const load_obj_options& options = {});

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

inline void parse_value(string_view& str, obj_vertex& value) {
    value = obj_vertex{0, 0, 0};
    parse_value(str, value.position);
    if (!str.empty() && str.front() == '/') {
        str.remove_prefix(1);
        if (!str.empty() && str.front() == '/') {
            str.remove_prefix(1);
            parse_value(str, value.normal);
        } else {
            parse_value(str, value.texturecoord);
            if (!str.empty() && str.front() == '/') {
                str.remove_prefix(1);
                parse_value(str, value.normal);
            }
        }
    }
}

// Input for OBJ textures
inline void parse_value(string_view& str, obj_texture_info& info) {
    // initialize
    info = obj_texture_info();

    // get tokens
    auto tokens = vector<string>();
    skip_whitespace(str);
    while (!str.empty()) {
        auto token = ""s;
        parse_value(str, token);
        tokens.push_back(token);
        skip_whitespace(str);
    }
    if (tokens.empty()) throw io_error("cannot parse value");

    // texture name
    info.path = normalize_path(tokens.back());

    // texture options
    auto last = string();
    for (auto i = 0; i < tokens.size() - 1; i++) {
        if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
        if (tokens[i] == "-clamp") info.clamp = true;
    }
}

// Load obj materials
template <typename Callbacks>
void load_mtl(
    const string& filename, Callbacks& cb, const load_obj_options& options) {
    // open file
    auto fs = input_file(filename);

    // currently parsed material
    auto material = obj_material();
    auto first    = true;

    // read the file line by line
    char buffer[4096];
    while (read_line(fs, buffer, sizeof(buffer))) {
        // line
        auto line = string_view{buffer};
        remove_comment_and_newline(line);
        skip_whitespace(line);
        if (line.empty()) continue;

        // get command
        auto cmd = ""s;
        parse_value(line, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "newmtl") {
            if (!first) cb.material(material);
            first    = false;
            material = obj_material();
            parse_value(line, material.name);
        } else if (cmd == "illum") {
            parse_value(line, material.illum);
        } else if (cmd == "Ke") {
            parse_value(line, material.ke);
        } else if (cmd == "Kd") {
            parse_value(line, material.kd);
        } else if (cmd == "Ks") {
            parse_value(line, material.ks);
        } else if (cmd == "Kt") {
            parse_value(line, material.kt);
        } else if (cmd == "Tf") {
            material.kt = {-1, -1, -1};
            parse_value(line, material.kt);
            if (material.kt.y < 0)
                material.kt = {material.kt.x, material.kt.x, material.kt.x};
            if (options.flip_tr) material.kt = vec3f{1, 1, 1} - material.kt;
        } else if (cmd == "Tr") {
            parse_value(line, material.op);
            if (options.flip_tr) material.op = 1 - material.op;
        } else if (cmd == "Ns") {
            parse_value(line, material.ns);
            material.rs = pow(2 / (material.ns + 2), 1 / 4.0f);
            if (material.rs < 0.01f) material.rs = 0;
            if (material.rs > 0.99f) material.rs = 1;
        } else if (cmd == "d") {
            parse_value(line, material.op);
        } else if (cmd == "Pr" || cmd == "rs") {
            parse_value(line, material.rs);
        } else if (cmd == "map_Ke") {
            parse_value(line, material.ke_txt);
        } else if (cmd == "map_Kd") {
            parse_value(line, material.kd_txt);
        } else if (cmd == "map_Ks") {
            parse_value(line, material.ks_txt);
        } else if (cmd == "map_Tr") {
            parse_value(line, material.kt_txt);
        } else if (cmd == "map_d" || cmd == "map_Tr") {
            parse_value(line, material.op_txt);
        } else if (cmd == "map_Pr" || cmd == "map_rs") {
            parse_value(line, material.rs_txt);
        } else if (cmd == "map_occ" || cmd == "occ") {
            parse_value(line, material.occ_txt);
        } else if (cmd == "map_bump" || cmd == "bump") {
            parse_value(line, material.bump_txt);
        } else if (cmd == "map_disp" || cmd == "disp") {
            parse_value(line, material.disp_txt);
        } else if (cmd == "map_norm" || cmd == "norm") {
            parse_value(line, material.norm_txt);
        } else if (cmd == "Ve") {
            parse_value(line, material.ve);
        } else if (cmd == "Va") {
            parse_value(line, material.va);
        } else if (cmd == "Vd") {
            parse_value(line, material.vd);
        } else if (cmd == "Vg") {
            parse_value(line, material.vg);
        } else if (cmd == "map_Vd") {
            parse_value(line, material.vd_txt);
        }
    }

    // issue current material
    if (!first) cb.material(material);
}

// Load obj extensions
template <typename Callbacks>
inline void load_objx(
    const string& filename, Callbacks& cb, const load_obj_options& options) {
    // open file
    auto fs = input_file(filename);

    // read the file line by line
    char buffer[4096];
    while (read_line(fs, buffer, sizeof(buffer))) {
        // line
        auto line = string_view{buffer};
        remove_comment_and_newline(line);
        skip_whitespace(line);
        if (line.empty()) continue;

        // get command
        auto cmd = ""s;
        parse_value(line, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "c") {
            auto camera = obj_camera();
            parse_value(line, camera.name);
            parse_value(line, camera.ortho);
            parse_value(line, camera.width);
            parse_value(line, camera.height);
            parse_value(line, camera.focal);
            parse_value(line, camera.focus);
            parse_value(line, camera.aperture);
            parse_value(line, camera.frame);
            cb.camera(camera);
        } else if (cmd == "e") {
            auto environment = obj_environment();
            parse_value(line, environment.name);
            parse_value(line, environment.ke);
            parse_value(line, environment.ke_txt.path);
            parse_value(line, environment.frame);
            if (environment.ke_txt.path == "\"\"") environment.ke_txt.path = "";
            cb.environmnet(environment);
        } else if (cmd == "po") {
            auto procedural = obj_procedural();
            parse_value(line, procedural.name);
            parse_value(line, procedural.type);
            parse_value(line, procedural.material);
            parse_value(line, procedural.size);
            parse_value(line, procedural.level);
            parse_value(line, procedural.frame);
            cb.procedural(procedural);
        } else {
            // unused
        }
    }
}

// Load obj scene
template <typename Callbacks>
inline void load_obj(
    const string& filename, Callbacks& cb, const load_obj_options& options) {
    // open file
    auto fs = input_file(filename);

    // track vertex size
    auto vert_size = obj_vertex();
    auto verts     = vector<obj_vertex>();  // buffer to avoid reallocation

    // read the file line by line
    char buffer[4096];
    while (read_line(fs, buffer, sizeof(buffer))) {
        // line
        auto line = string_view{buffer};
        remove_comment_and_newline(line);
        skip_whitespace(line);
        if (line.empty()) continue;

        // get command
        auto cmd = ""s;
        parse_value(line, cmd);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            auto vert = zero3f;
            parse_value(line, vert);
            cb.vert(vert);
            vert_size.position += 1;
        } else if (cmd == "vn") {
            auto vert = zero3f;
            parse_value(line, vert);
            cb.norm(vert);
            vert_size.normal += 1;
        } else if (cmd == "vt") {
            auto vert = zero2f;
            parse_value(line, vert);
            if (options.flip_texcoord) vert.y = 1 - vert.y;
            cb.texcoord(vert);
            vert_size.texturecoord += 1;
        } else if (cmd == "f" || cmd == "l" || cmd == "p") {
            verts.clear();
            skip_whitespace(line);
            while (!line.empty()) {
                auto vert = obj_vertex{};
                parse_value(line, vert);
                if (!vert.position) break;
                if (vert.position < 0)
                    vert.position = vert_size.position + vert.position + 1;
                if (vert.texturecoord < 0)
                    vert.texturecoord = vert_size.texturecoord +
                                        vert.texturecoord + 1;
                if (vert.normal < 0)
                    vert.normal = vert_size.normal + vert.normal + 1;
                verts.push_back(vert);
                skip_whitespace(line);
            }
            if (cmd == "f") cb.face(verts);
            if (cmd == "l") cb.line(verts);
            if (cmd == "p") cb.point(verts);
        } else if (cmd == "o") {
            auto name = ""s;
            parse_value_or_empty(line, name);
            cb.object(name);
        } else if (cmd == "usemtl") {
            auto name = ""s;
            parse_value_or_empty(line, name);
            cb.usemtl(name);
        } else if (cmd == "g") {
            auto name = ""s;
            parse_value_or_empty(line, name);
            cb.group(name);
        } else if (cmd == "s") {
            auto name = ""s;
            parse_value_or_empty(line, name);
            cb.smoothing(name);
        } else if (cmd == "mtllib") {
            if (options.geometry_only) continue;
            auto mtlname = ""s;
            parse_value(line, mtlname);
            cb.mtllib(mtlname);
            auto mtlpath = get_dirname(filename) + mtlname;
            load_mtl(mtlpath, cb, options);
        } else {
            // unused
        }
    }

    // parse extensions if presents
    if (!options.geometry_only) {
        auto extname    = get_noextension(filename) + ".objx";
        auto ext_exists = exists_file(extname);
        if (ext_exists) {
            load_objx(extname, cb, options);
        }
    }
}

}  // namespace yocto

namespace std {

// Hash functor for vector for use with unordered_map
template <>
struct hash<yocto::obj_vertex> {
    static constexpr std::hash<int> hasher = std::hash<int>();

    size_t operator()(const yocto::obj_vertex& v) const {
        auto h = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= hasher((&v.position)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

#endif
