//
// # Yocto/Obj: Tiny C++ Library for OBJ loading/saving
//
// Yocto/Obj is a library for loading and saving OBJ files for 3D applications,
// support for points, lines, triangles and general polygons and all materials
// properties. Yocto/Obj defined also a few extensions to easily create demos
// such as per-vertex color and radius, cameras, environment maps and instances.
// Can use either a low-level OBJ representation, from this files,
// or a high level flattened representation included in Yocto/Scene.
//
// Both in reading and writing, OBJ has no clear convention on the orientation
// of textures Y axis. So in many cases textures appears flipped. To handle
// that, use the option to flip textures coordinates on either saving or
// loading. By default texture coordinates are flipped since this seems
// the convention found on test cases collected on the web. The value Tr
// has similar problems, since its relation to opacity is software specific.
// Again we let the user choose the conversion and set the default to the
// one found on the web.
//
// Yocto/Obj stores all parsed data in its data structures. This design
// decision does not scale particularly well with very large models. The main
// reasons is the OBJ file format vertices are global for the entire scene and
// not single objects. To avoid memory allocation problems, we internally store
// vertices in queues during reading and make them unique per object right after
// parsing.
//
// In the high level interface, shapes are indexed meshes and are described
// by arrays of vertex indices for points/lines/triangles and arrays for vertex
// positions, normals, texcoords, color and radius. The latter two as
// extensions. Since OBJ is a complex formats that does not match well with
// current GPU rendering / path tracing algorithms, we adopt a simplification
// similar to other single file libraries:
// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
//   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
//   vertex duplication happens thought for same triplets
// 2. we split shapes on changes to groups and materials, instead of keeping
//   per-face group/material data; this makes the data usable right away in
//   a GPU viewer; this is not a major limitation if we accept the previous
//   point that already changes shapes topology.
//
// ## Usage
//
// 1. load a obj data with `load_obj()`; can load also textues
// 2. look at the `obj_XXX` data structures for access to individual elements
// 3. use obj back to disk with `save_obj()`; can also save textures
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
//

#ifndef _YOBJ_H_
#define _YOBJ_H_

// enable fast OBJ parsing
#ifndef YGL_FASTOBJ
#define YGL_FASTOBJ 1
#endif

#include "ygl.h"

#include <algorithm>
#include <array>
#include <deque>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// WAVEFRONT OBJ SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

// Obj element type.
enum struct obj_element_type : uint8_t {
    point = 1,
    line = 2,
    face = 3,
    bezier = 4
};

// Obj element (point/polyline/polygon)
struct obj_element {
    uint32_t start;          // starting vertex index
    uint8_t size;            // number of vertices
    obj_element_type type;   // element type
    uint16_t material = 0;   // material id
    uint16_t group = 0;      // group id
    uint16_t smoothing = 0;  // smoothing group id
};

// Obj group properties.
struct obj_group_props {
    std::string name = "";     // group name
    std::string matname = "";  // material name
    bool faceted = false;      // faceted or smooth
};

// Obj object.
struct obj_object {
    std::string name;                    // name
    std::vector<std::string> materials;  // materials
    std::vector<std::string> groups;     // groups
    std::vector<bool> smoothing;         // smoothing groups
    std::vector<int> verts_pos;          // element pos
    std::vector<int> verts_norm;         // element norm
    std::vector<int> verts_texcoord;     // element texcoord
    std::vector<obj_element> elems;      // elements
    frame3f frame = identity_frame3f;    // frame [extension]
    vec2i subdiv = zero2i;  // type/level of subdivision [extension]
    // Properties not explicitly handled [extension].
    std::unordered_map<std::string, std::vector<float>> props;
};

// Obj texture. Texture data is loaded only if desired.
struct obj_texture {
    std::string path;  // path
    image<vec4f> img;       // image data
};

// Obj texture information.
struct obj_texture_info {
    std::string path = "";  // file path
    bool clamp = false;     // clamp to edge
    float scale = 1;        // scale for bump/displacement
    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<float>> props;
};

// Comparison for texture info.
inline bool operator==(const obj_texture_info& a, const obj_texture_info& b) {
    if (a.path.empty() && b.path.empty()) return true;
    if (a.path != b.path) return false;
    return a.clamp == b.clamp && a.scale == b.scale && a.props == b.props;
}

// Obj material.
struct obj_material {
    std::string name;  // name
    int illum = 0;     // MTL illum mode

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f ka = {0, 0, 0};  // ambient color
    vec3f kd = {0, 0, 0};  // diffuse color
    vec3f ks = {0, 0, 0};  // specular color
    vec3f kr = {0, 0, 0};  // reflection color
    vec3f kt = {0, 0, 0};  // transmission color
    float ns = 0;          // Phong exponent color
    float ior = 1;         // index of refraction
    float op = 1;          // opacity
    float rs = -1;         // roughness (-1 not defined)
    float km = -1;         // metallic  (-1 not defined)

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

    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Obj camera [extension].
struct obj_camera {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    bool ortho = false;                // orthographic
    float width = 0.036f;              // film width (default to 35mm)
    float height = 0.024f;             // film height (default to 35mm)
    float focal = 0.050f;              // focal length
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float aperture = 0;                // lens aperture
    float focus = maxf;                // focus distance
};

// Obj environment [extension].
struct obj_environment {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    vec3f ke = zero3f;                 // emission color
    obj_texture_info ke_txt;           // emission texture
};

// Obj scene.
struct obj_scene {
    std::deque<vec3f> pos;       // vertex positions
    std::deque<vec3f> norm;      // vertex normals
    std::deque<vec2f> texcoord;  // vertex texcoords

    std::vector<obj_object*> objects;            // objects
    std::vector<obj_material*> materials;        // materials
    std::vector<obj_texture*> textures;          // textures
    std::vector<obj_camera*> cameras;            // cameras [extension]
    std::vector<obj_environment*> environments;  // environments [extension]
};

// Load an OBJ from file `filename`. Split shapes at material and group
// boundaries, if `split_shapes` is true.
// Load textures if `load_textures` is true, and report errors only if
// `skip_missing` is false. Texture coordinates and material Tr are flipped
// if `flip_texcoord` and `flip_tp` are respectively true.
inline obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool flip_texcoord = true, bool flip_tr = true);

// Save an OBJ to file `filename`. Save textures if `save_textures` is true,
// and report errors only if `skip_missing` is false.
// Texture coordinates and material Tr are flipped if `flip_texcoord` and
// `flip_tp` are respectively true.
inline void save_obj(const std::string& filename, const obj_scene* obj,
    bool flip_texcoord = true, bool flip_tr = true);

// Load OBJ texture images.
inline void load_obj_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing = true);
// Save OBJ texture images.
inline void save_obj_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing = true);

// Callback object
struct obj_callbacks {
    bool (*object)(void*, const std::string&) = nullptr;  // object callback
    bool (*group)(void*, const std::string&) = nullptr;   // group callback
    bool (*usemat)(void*, const std::string&) = nullptr;  // use material cb
    bool (*smoothing)(void*, bool) = nullptr;             // smoothing callback
    bool (*vertex)(void*, vec3f) = nullptr;               // vertex callback
    bool (*normal)(void*, vec3f) = nullptr;               // normal callback
    bool (*texcoord)(void*, vec2f) = nullptr;             // texcoord callback
    bool (*face)(void*, int, const vec3i*) = nullptr;     // face callback
    bool (*line)(void*, int, const vec3i*) = nullptr;     // line callback
    bool (*point)(void*, int, const vec3i*) = nullptr;    // point callback
    bool (*material)(
        void*, const obj_material&) = nullptr;           // material callback
    bool (*camera)(void*, const obj_camera&) = nullptr;  // camera callback
    bool (*environment)(void*, const obj_environment&) = nullptr;  // envmap cb
};

// Load OBJ with callbacks
inline void load_obj(
    const std::string& filename, const obj_callbacks& callbacks, void* ctx);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

namespace detail {

// skip whitespace
inline void obj_skipws(char*& s) {
    while (*s == ' ') s++;
}

// skip a string if matched
inline bool obj_streq(const char* s, const char* str) {
    while (*s == *str && *s && *str) {
        s++;
        str++;
    }
    return *s == *str;
}

#if YGL_FASTOBJ

// parse base value
inline void obj_parse_base(char*& s, int& val) {
    while (*s == ' ') s++;
    val = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
}

// parse base value
inline void obj_parse_base(char*& s, float& val) {
    while (*s == ' ') s++;
    //    auto ss = s; auto sss = ss;
    auto mantissa = 0, fractional = 0, fractional_length = 0, exponent = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') mantissa = mantissa * 10 + (*s++ - '0');
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            fractional = fractional * 10 + (*s++ - '0');
            fractional_length++;
        }
    }
    mantissa *= sn;
    fractional *= sn;
    if (*s == 'e' || *s == 'E') {
        s++;
        auto en = (*s == '-') ? -1 : 1;
        if (*s == '-' || *s == '+') s++;
        while (*s >= '0' && *s <= '9') exponent = exponent * 10 + (*s++ - '0');
        exponent *= en;
    }
    auto dval = (double)mantissa;
    if (fractional)
        dval += fractional * std::pow(10.0, -(double)fractional_length);
    if (exponent) dval *= std::pow(10.0, (double)exponent);
    val = (float)dval;
#if 0
    auto cval = val;
    sscanf(ss, "%f", &cval);
    if(abs(val - cval) > 0.01f) {
        printf("- %g %g %s\n", val, cval, sss);
    }
    auto len = 0;
    sscanf(s, "%f%n", &val, &len);
    s += len;
#endif
}

// parse base value
inline void obj_parse_base(char*& s, char* val) {
    while (*s == ' ') s++;
    while (*s && *s != ' ') *val++ = *s++;
    *val = 0;
}
// parse base value
inline void obj_parse_base(char*& s, std::string& val) {
    char buf[4096];
    obj_parse_base(s, buf);
    val = buf;
}

#else

// parse base value
inline void obj_parse_base(char*& s, int& val) {
    auto len = 0;
    sscanf(s, "%d%n", &val, &len);
    s += len;
}

// parse base value
inline void obj_parse_base(char*& s, float& val) {
    auto len = 0;
    sscanf(s, "%f%n", &val, &len);
    s += len;
}

// parse base value
inline void obj_parse_base(char*& s, std::string& val) {
    char buf[4096];
    auto len = 0;
    sscanf(s, "%s%n", buf, &len);
    if (len) {
        s += len;
        val = buf;
    } else {
        val = "";
    }
}

// parse base value
inline void obj_parse_base(char*& s, char* val) {
    auto len = 0;
    sscanf(s, "%s%n", val, &len);
    s += len;
}

#endif

// parse value
inline void obj_parse(char*& s, int& val) { obj_parse_base(s, val); }
inline void obj_parse(char*& s, float& val) { obj_parse_base(s, val); }
inline void obj_parse(char*& s, bool& val) {
    auto i = 0;
    obj_parse(s, i);
    val = i;
}
inline void obj_parse(char*& s, std::string& val) { obj_parse_base(s, val); }
inline void obj_parse(char*& s, char* val) { obj_parse_base(s, val); }
inline void obj_parse(char*& s, vec2i& val) {
    for (auto i = 0; i < 2; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec2f& val) {
    for (auto i = 0; i < 2; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec3f& val) {
    for (auto i = 0; i < 3; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, frame3f& val) {
    for (auto i = 0; i < 12; i++) obj_parse(s, (&val.x.x)[i]);
}
inline void obj_parse(char*& s, vec3i& val, vec3i vert_size) {
    char buf[1024];
    obj_parse_base(s, buf);
    val = {-1, -1, -1};
    auto i = 0;
    auto sb = buf;
    while (i < 3 && *sb) {
        obj_parse_base(sb, (&val.x)[i]);
        (&val.x)[i] = ((&val.x)[i] < 0) ? (&vert_size.x)[i] + (&val.x)[i] :
                                          (&val.x)[i] - 1;
        if (*sb != '/') break;
        while (*sb == '/') {
            sb++;
            i++;
        }
    }
}
inline void obj_parse(char*& s, int& num, vec3i* vert_buf, vec3i vert_size) {
    num = 0;
    while (*s) {
        obj_parse(s, vert_buf[num], vert_size);
        num += 1;
    }
}

// clear the whitespace
inline void obj_convertws(char* s) {
    while (*s) {
        if (*s == '\t' || *s == '\r' || *s == '\n') *s = ' ';
        s++;
    }
}

// Parse texture options and name
inline void obj_parse(char*& s, obj_texture_info& info) {
    // initialize
    info = obj_texture_info();

    // get tokens
    auto tokens = std::vector<std::string>();
    obj_skipws(s);
    while (*s) {
        tokens.push_back("");
        obj_parse(s, tokens.back());
        obj_skipws(s);
    }

    // exit if no tokens
    if (tokens.empty()) return;

    // texture name
    info.path = tokens.back();
    for (auto& c : info.path)
        if (c == '\\') c = '/';

    // texture options
    auto last = std::string();
    for (auto& tok : tokens) {
        if (tok == tokens.back()) break;
        if (tok[0] == '-') {
            last = tok;
            info.props[last] = {};
        } else {
            info.props[last].push_back(atof(tok.c_str()));
        }
    }

    // clamp
    if (info.props.find("clamp") != info.props.end()) {
        info.clamp = info.props.at("clamp").front();
        info.props.erase("clamp");
    }

    if (info.props.find("bm") != info.props.end()) {
        info.scale = info.props.at("bm").front();
        info.props.erase("bm");
    }
}

// Load MTL
inline std::vector<obj_material*> load_mtl(
    const std::string& filename, bool flip_tr) {
    // clear materials
    auto materials = std::vector<obj_material*>();

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // add a material preemptively to avoid crashes
    materials.push_back(new obj_material());
    auto mat = materials.back();

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    while (fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "newmtl")) {
            materials.push_back(std::make_shared<obj_material>());
            mat = materials.back();
            obj_parse(ss, mat->name);
        } else if (obj_streq(cmd, "illum")) {
            obj_parse(ss, mat->illum);
        } else if (obj_streq(cmd, "Ke")) {
            obj_parse(ss, mat->ke);
        } else if (obj_streq(cmd, "Ka")) {
            obj_parse(ss, mat->ka);
        } else if (obj_streq(cmd, "Kd")) {
            obj_parse(ss, mat->kd);
        } else if (obj_streq(cmd, "Ks")) {
            obj_parse(ss, mat->ks);
        } else if (obj_streq(cmd, "Kr")) {
            obj_parse(ss, mat->kr);
        } else if (obj_streq(cmd, "Kt")) {
            obj_parse(ss, mat->kt);
        } else if (obj_streq(cmd, "Tf")) {
            auto nchan = 0;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&mat->kt.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) mat->kt = {mat->kt.x, mat->kt.x, mat->kt.x};
            if (flip_tr)
                materials.back()->kt = vec3f{1, 1, 1} - materials.back()->kt;
        } else if (obj_streq(cmd, "Tr")) {
            auto nchan = 0;
            auto tr = zero3f;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&tr.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) tr = {tr.x, tr.x, tr.x};
            materials.back()->op = (tr.x + tr.y + tr.z) / 3;
            if (flip_tr) materials.back()->op = 1 - materials.back()->op;
        } else if (obj_streq(cmd, "Ns")) {
            obj_parse(ss, mat->ns);
        } else if (obj_streq(cmd, "d")) {
            obj_parse(ss, mat->op);
        } else if (obj_streq(cmd, "Ni")) {
            obj_parse(ss, mat->ior);
        } else if (obj_streq(cmd, "Pr") || obj_streq(cmd, "rs")) {
            obj_parse(ss, mat->rs);
        } else if (obj_streq(cmd, "Pm") || obj_streq(cmd, "Km")) {
            obj_parse(ss, mat->km);
        } else if (obj_streq(cmd, "map_Ke")) {
            obj_parse(ss, mat->ke_txt);
        } else if (obj_streq(cmd, "map_Ka")) {
            obj_parse(ss, mat->ka_txt);
        } else if (obj_streq(cmd, "map_Kd")) {
            obj_parse(ss, mat->kd_txt);
        } else if (obj_streq(cmd, "map_Ks")) {
            obj_parse(ss, mat->ks_txt);
        } else if (obj_streq(cmd, "map_Kr")) {
            obj_parse(ss, mat->kr_txt);
        } else if (obj_streq(cmd, "map_Tr")) {
            obj_parse(ss, mat->kt_txt);
        } else if (obj_streq(cmd, "map_Ns")) {
            obj_parse(ss, mat->ns_txt);
        } else if (obj_streq(cmd, "map_d") || obj_streq(cmd, "map_Tr")) {
            obj_parse(ss, mat->op_txt);
        } else if (obj_streq(cmd, "map_Ni")) {
            obj_parse(ss, mat->ior_txt);
        } else if (obj_streq(cmd, "map_Pr") || obj_streq(cmd, "map_rs")) {
            obj_parse(ss, mat->rs_txt);
        } else if (obj_streq(cmd, "map_Pm") || obj_streq(cmd, "map_Km")) {
            obj_parse(ss, mat->km_txt);
        } else if (obj_streq(cmd, "map_occ") || obj_streq(cmd, "occ")) {
            obj_parse(ss, mat->occ_txt);
        } else if (obj_streq(cmd, "map_bump") || obj_streq(cmd, "bump")) {
            obj_parse(ss, mat->bump_txt);
        } else if (obj_streq(cmd, "map_disp") || obj_streq(cmd, "disp")) {
            obj_parse(ss, mat->disp_txt);
        } else if (obj_streq(cmd, "map_norm") || obj_streq(cmd, "norm")) {
            obj_parse(ss, mat->norm_txt);
        } else {
            // copy into strings
            obj_skipws(ss);
            while (*ss) {
                mat->props[cmd].push_back("");
                obj_parse(ss, mat->props[cmd].back());
                obj_skipws(ss);
            }
        }
    }

    // remove first fake material
    materials.erase(materials.begin());

    // clone
    fclose(fs);

    // done
    return materials;
}

// Loads an OBJ
inline obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool flip_texcoord, bool flip_tr) {
    // clear obj
    auto obj = std::make_shared<obj_scene>();

    // splitting policy
    auto split_material = split_shapes;
    auto split_group = split_shapes;
    auto split_smoothing = split_shapes;

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // add object if needed
    auto add_object = [&](obj_scene* obj, std::string name, std::string matname,
                          std::string groupname, bool smoothing) {
        if (obj->objects.empty() || !obj->objects.back()->elems.empty())
            obj->objects.push_back(new obj_object());
        auto oobj = obj->objects.back();
        oobj->name = name;
        if (oobj->materials.empty()) oobj->materials.push_back("");
        if (oobj->groups.empty()) oobj->groups.push_back("");
        if (oobj->smoothing.empty()) oobj->smoothing.push_back(true);
        oobj->materials.back() = matname;
        oobj->groups.back() = groupname;
        oobj->smoothing.back() = smoothing;
        return oobj;
    };

    // current parsing values
    auto matname = std::string();
    auto oname = std::string();
    auto gname = std::string();
    auto smoothing = true;
    auto oobj = add_object(obj, oname, matname, gname, smoothing);

    // properties
    auto oprops = std::unordered_map<std::string,
        std::unordered_map<std::string, std::vector<float>>>();

    // keep track of array lengths
    auto vert_size = zero3i;

    // elem type map
    static auto elem_type_map =
        std::unordered_map<std::string, obj_element_type>{
            {"f", obj_element_type::face}, {"l", obj_element_type::line},
            {"p", obj_element_type::point}, {"b", obj_element_type::bezier}};

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    while (fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "v")) {
            vert_size.x += 1;
            obj->pos.push_back(zero3f);
            obj_parse(ss, obj->pos.back());
        } else if (obj_streq(cmd, "vn")) {
            vert_size.z += 1;
            obj->norm.push_back(zero3f);
            obj_parse(ss, obj->norm.back());
        } else if (obj_streq(cmd, "vt")) {
            vert_size.y += 1;
            obj->texcoord.push_back(zero2f);
            obj_parse(ss, obj->texcoord.back());
            if (flip_texcoord)
                obj->texcoord.back().y = 1 - obj->texcoord.back().y;
        } else if (obj_streq(cmd, "f") || obj_streq(cmd, "l") ||
                   obj_streq(cmd, "p") || obj_streq(cmd, "b")) {
            auto elem = obj_element();
            elem.type = elem_type_map.at(cmd);
            elem.start = (uint32_t)oobj->verts_pos.size();
            elem.size = 0;
            elem.material = (int)oobj->materials.size() - 1;
            elem.group = (int)oobj->groups.size() - 1;
            elem.smoothing = (int)oobj->smoothing.size() - 1;
            oobj->elems.push_back(elem);
            obj_skipws(ss);
            while (*ss) {
                auto vert = vec3i{-1, -1, -1};
                obj_parse(ss, vert, vert_size);
                obj_skipws(ss);
                oobj->verts_pos.push_back(vert.x);
                oobj->verts_norm.push_back(vert.z);
                oobj->verts_texcoord.push_back(vert.y);
                oobj->elems.back().size += 1;
            }
        } else if (obj_streq(cmd, "o")) {
            obj_parse(ss, oname);
            gname = "";
            matname = "";
            smoothing = true;
            oobj = add_object(obj, oname, matname, gname, smoothing);
        } else if (obj_streq(cmd, "usemtl")) {
            obj_parse(ss, matname);
            if (split_material) {
                oobj = add_object(obj, oname, matname, gname, smoothing);
            } else {
                if (oobj->elems.empty()) {
                    oobj->materials.back() = matname;
                } else {
                    oobj->materials.push_back(matname);
                }
            }
        } else if (obj_streq(cmd, "g")) {
            obj_parse(ss, gname);
            if (split_group) {
                oobj = add_object(obj, oname, matname, gname, smoothing);
            } else {
                if (oobj->elems.empty()) {
                    oobj->groups.back() = gname;
                } else {
                    oobj->groups.push_back(gname);
                }
            }
        } else if (obj_streq(cmd, "s")) {
            auto name = std::string();
            obj_parse(ss, name);
            smoothing = (name == "on");
            if (split_smoothing) {
                oobj = add_object(obj, oname, matname, gname, smoothing);
            } else {
                if (oobj->elems.empty()) {
                    oobj->smoothing.back() = smoothing;
                } else {
                    oobj->smoothing.push_back(smoothing);
                }
            }
        } else if (obj_streq(cmd, "of")) {
            obj_parse(ss, oobj->frame);
        } else if (obj_streq(cmd, "os")) {
            obj_parse(ss, oobj->subdiv);
        } else if (obj_streq(cmd, "op")) {
            auto rname = std::string(), pname = std::string();
            obj_parse(ss, rname);
            obj_skipws(ss);
            obj_parse(ss, pname);
            obj_skipws(ss);
            auto& pvalues = oprops[rname][pname];
            while (*ss) {
                auto tok = std::string();
                obj_parse(ss, tok);
                obj_skipws(ss);
                pvalues.push_back((float)std::atof(tok.c_str()));
            }
        } else if (obj_streq(cmd, "mtllib")) {
            auto mtlname = std::string();
            obj_parse(ss, mtlname);
            auto mtlpath = get_dirname(filename) + "/" + mtlname;
            auto mats = load_mtl(mtlpath, flip_tr);
            obj->materials.insert(
                obj->materials.end(), mats.begin(), mats.end());
        } else if (obj_streq(cmd, "c")) {
            auto cam = new obj_camera();
            obj_parse(ss, cam->name);
            obj_parse(ss, cam->ortho);
            obj_parse(ss, cam->width);
            obj_parse(ss, cam->height);
            obj_parse(ss, cam->focal);
            obj_parse(ss, cam->focus);
            obj_parse(ss, cam->aperture);
            obj_parse(ss, cam->frame);
            obj->cameras.push_back(cam);
        } else if (obj_streq(cmd, "e")) {
            auto env = new obj_environment();
            obj_parse(ss, env->name);
            obj_parse(ss, env->ke);
            obj_parse(ss, env->ke_txt.path);
            if (env->ke_txt.path == "\"\"") env->ke_txt.path = "";
            obj_parse(ss, env->frame);
            obj->environments.push_back(env);
        } else {
            // unused
        }
    }

    // cleanup empty
    auto clear_vert_if_unused = [](std::vector<int>& vert) {
        if (vert.empty()) return;
        auto used = false;
        for (auto v : vert)
            if (v >= 0) used = true;
        if (!used) vert.clear();
    };
    for (auto idx = 0; idx < obj->objects.size(); idx++) {
        auto oobj = obj->objects[idx];
        clear_vert_if_unused(oobj->verts_pos);
        clear_vert_if_unused(oobj->verts_norm);
        clear_vert_if_unused(oobj->verts_texcoord);
        if (!oobj->elems.empty() || !oobj->verts_pos.empty()) continue;
        obj->objects.erase(obj->objects.begin() + idx);
        idx--;
    }
    auto end = std::remove_if(obj->objects.begin(), obj->objects.end(),
        [](const obj_object*& x) { return !x; });
    obj->objects.erase(end, obj->objects.end());

    // apply properties
    for (auto oobj : obj->objects) oobj->props = oprops[oobj->name];

    // create textures
    auto txt_set = std::unordered_map<std::string, obj_texture*>();
    auto add_texture = [](obj_scene* obj, auto& txt_set,
                           const obj_texture_info& info) {
        if (info.path == "") return;
        if (txt_set.find(info.path) != txt_set.end()) return;
        auto txt = std::make_shared<obj_texture>();
        txt->path = info.path;
        txt_set[txt->path] = txt;
        obj->textures.push_back(txt);
    };
    for (auto mat : obj->materials) {
        add_texture(obj, txt_set, mat->ke_txt);
        add_texture(obj, txt_set, mat->kd_txt);
        add_texture(obj, txt_set, mat->ks_txt);
        add_texture(obj, txt_set, mat->kr_txt);
        add_texture(obj, txt_set, mat->kt_txt);
        add_texture(obj, txt_set, mat->rs_txt);
        add_texture(obj, txt_set, mat->bump_txt);
        add_texture(obj, txt_set, mat->norm_txt);
        add_texture(obj, txt_set, mat->disp_txt);
        add_texture(obj, txt_set, mat->occ_txt);
    }
    for (auto env : obj->environments) {
        add_texture(obj, txt_set, env->ke_txt);
    }

    // close file
    fclose(fs);

    // done
    return obj;
}

// Load MTL
inline void load_mtl(const std::string& filename,
    const obj_callbacks& callbacks, void* ctx, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // add a material preemptively to avoid crashes
    auto mat = obj_material();
    auto parsed_one = false;

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    auto done = false;
    while (!done && fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "newmtl")) {
            if (parsed_one && callbacks.material)
                done = callbacks.material(ctx, mat);
            parsed_one = true;
            mat = obj_material();
            obj_parse(ss, mat.name);
        } else if (obj_streq(cmd, "illum")) {
            obj_parse(ss, mat.illum);
        } else if (obj_streq(cmd, "Ke")) {
            obj_parse(ss, mat.ke);
        } else if (obj_streq(cmd, "Ka")) {
            obj_parse(ss, mat.ka);
        } else if (obj_streq(cmd, "Kd")) {
            obj_parse(ss, mat.kd);
        } else if (obj_streq(cmd, "Ks")) {
            obj_parse(ss, mat.ks);
        } else if (obj_streq(cmd, "Kr")) {
            obj_parse(ss, mat.kr);
        } else if (obj_streq(cmd, "Kt")) {
            obj_parse(ss, mat.kt);
        } else if (obj_streq(cmd, "Tf")) {
            auto nchan = 0;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&mat.kt.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) mat.kt = {mat.kt.x, mat.kt.x, mat.kt.x};
            if (flip_tr) mat.kt = vec3f{1, 1, 1} - mat.kt;
        } else if (obj_streq(cmd, "Tr")) {
            auto nchan = 0;
            auto tr = zero3f;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&tr.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) tr = {tr.x, tr.x, tr.x};
            mat.op = (tr.x + tr.y + tr.z) / 3;
            if (flip_tr) mat.op = 1 - mat.op;
        } else if (obj_streq(cmd, "Ns")) {
            obj_parse(ss, mat.ns);
        } else if (obj_streq(cmd, "d")) {
            obj_parse(ss, mat.op);
        } else if (obj_streq(cmd, "Ni")) {
            obj_parse(ss, mat.ior);
        } else if (obj_streq(cmd, "Pr") || obj_streq(cmd, "rs")) {
            obj_parse(ss, mat.rs);
        } else if (obj_streq(cmd, "Pm") || obj_streq(cmd, "Km")) {
            obj_parse(ss, mat.km);
        } else if (obj_streq(cmd, "map_Ke")) {
            obj_parse(ss, mat.ke_txt);
        } else if (obj_streq(cmd, "map_Ka")) {
            obj_parse(ss, mat.ka_txt);
        } else if (obj_streq(cmd, "map_Kd")) {
            obj_parse(ss, mat.kd_txt);
        } else if (obj_streq(cmd, "map_Ks")) {
            obj_parse(ss, mat.ks_txt);
        } else if (obj_streq(cmd, "map_Kr")) {
            obj_parse(ss, mat.kr_txt);
        } else if (obj_streq(cmd, "map_Tr")) {
            obj_parse(ss, mat.kt_txt);
        } else if (obj_streq(cmd, "map_Ns")) {
            obj_parse(ss, mat.ns_txt);
        } else if (obj_streq(cmd, "map_d") || obj_streq(cmd, "map_Tr")) {
            obj_parse(ss, mat.op_txt);
        } else if (obj_streq(cmd, "map_Ni")) {
            obj_parse(ss, mat.ior_txt);
        } else if (obj_streq(cmd, "map_Pr") || obj_streq(cmd, "map_rs")) {
            obj_parse(ss, mat.rs_txt);
        } else if (obj_streq(cmd, "map_Pm") || obj_streq(cmd, "map_Km")) {
            obj_parse(ss, mat.km_txt);
        } else if (obj_streq(cmd, "map_occ") || obj_streq(cmd, "occ")) {
            obj_parse(ss, mat.occ_txt);
        } else if (obj_streq(cmd, "map_bump") || obj_streq(cmd, "bump")) {
            obj_parse(ss, mat.bump_txt);
        } else if (obj_streq(cmd, "map_disp") || obj_streq(cmd, "disp")) {
            obj_parse(ss, mat.disp_txt);
        } else if (obj_streq(cmd, "map_norm") || obj_streq(cmd, "norm")) {
            obj_parse(ss, mat.norm_txt);
        } else {
            // copy into strings
            obj_skipws(ss);
            while (*ss) {
                mat.props[cmd].push_back("");
                obj_parse(ss, mat.props[cmd].back());
                obj_skipws(ss);
            }
        }
    }

    // emit last
    if (parsed_one && callbacks.material) done = callbacks.material(ctx, mat);

    // clone
    fclose(fs);
}

// Loads an OBJ
inline void load_obj(const std::string& filename,
    const obj_callbacks& callbacks, void* ctx, bool flip_texcoord,
    bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // keep track of array lengths
    auto vert_size = zero3i;

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    auto done = false;
    vec3i vert_buf[128];
    while (!done && fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "v")) {
            vert_size.x += 1;
            if (callbacks.vertex) {
                auto val = zero3f;
                obj_parse(ss, val);
                done = callbacks.vertex(ctx, val);
            }
        } else if (obj_streq(cmd, "vn")) {
            vert_size.z += 1;
            if (callbacks.normal) {
                auto val = zero3f;
                obj_parse(ss, val);
                done = callbacks.normal(ctx, val);
            }
        } else if (obj_streq(cmd, "vt")) {
            vert_size.y += 1;
            if (callbacks.texcoord) {
                auto val = zero2f;
                if (flip_texcoord) val.y = 1 - val.y;
                done = callbacks.texcoord(ctx, val);
            }
        } else if (obj_streq(cmd, "f")) {
            if (callbacks.face) {
                obj_skipws(ss);
                auto num = 0;
                obj_parse(ss, num, vert_buf, vert_size);
                callbacks.face(ctx, num, vert_buf);
            }
        } else if (obj_streq(cmd, "l")) {
            if (callbacks.line) {
                obj_skipws(ss);
                auto num = 0;
                obj_parse(ss, num, vert_buf, vert_size);
                callbacks.line(ctx, num, vert_buf);
            }
        } else if (obj_streq(cmd, "p")) {
            if (callbacks.point) {
                obj_skipws(ss);
                auto num = 0;
                obj_parse(ss, num, vert_buf, vert_size);
                callbacks.point(ctx, num, vert_buf);
            }
        } else if (obj_streq(cmd, "o")) {
            if (callbacks.object) {
                auto name = std::string();
                obj_parse(ss, name);
                done = callbacks.object(ctx, name);
            }
        } else if (obj_streq(cmd, "usemtl")) {
            if (callbacks.usemat) {
                auto name = std::string();
                obj_parse(ss, name);
                done = callbacks.usemat(ctx, name);
            }
        } else if (obj_streq(cmd, "g")) {
            if (callbacks.group) {
                auto name = std::string();
                obj_parse(ss, name);
                done = callbacks.group(ctx, name);
            }
        } else if (obj_streq(cmd, "s")) {
            if (callbacks.smoothing) {
                auto name = std::string();
                obj_parse(ss, name);
                done = callbacks.smoothing(ctx, name != "off");
            }
        } else if (obj_streq(cmd, "mtllib")) {
            auto mtlname = std::string();
            obj_parse(ss, mtlname);
            auto mtlpath = get_dirname(filename) + "/" + mtlname;
            load_mtl(mtlpath, callbacks, ctx, flip_tr);
        } else if (obj_streq(cmd, "c")) {
            if (callbacks.camera) {
                auto cam = obj_camera();
                obj_parse(ss, cam.name);
                obj_parse(ss, cam.ortho);
                obj_parse(ss, cam.width);
                obj_parse(ss, cam.height);
                obj_parse(ss, cam.focal);
                obj_parse(ss, cam.focus);
                obj_parse(ss, cam.aperture);
                obj_parse(ss, cam.frame);
                done = callbacks.camera(ctx, cam);
            }
        } else if (obj_streq(cmd, "e")) {
            if (callbacks.environment) {
                auto env = obj_environment();
                obj_parse(ss, env.name);
                obj_parse(ss, env.ke);
                obj_parse(ss, env.ke_txt.path);
                if (env.ke_txt.path == "\"\"") env.ke_txt.path = "";
                obj_parse(ss, env.frame);
                done = callbacks.environment(ctx, env);
            }
        } else {
            // unused
        }
    }

    // close file
    fclose(fs);
}

// Dumps a value
inline void obj_dump(char*& s, char* val) {
    while (*val) *s++ = *val++;
}
inline void obj_dump(char*& s, const char* val) {
    while (*val) *s++ = *val++;
}
inline void obj_dump(char*& s, const std::string& val) {
    auto val_ = val.c_str();
    while (*val_) *s++ = *val_++;
}
inline void obj_dump(char*& s, int val) { s += sprintf(s, "%d", val); }
inline void obj_dump(char*& s, float val) { s += sprintf(s, "%g", val); }
inline void obj_dump(char*& s, const int* val, int num) {
    for (auto i = 0; i < num; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, val[i]);
    }
}
inline void obj_dump(char*& s, const float* val, int num) {
    for (auto i = 0; i < num; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, val[i]);
    }
}
inline void obj_dump(char*& s, vec2i val) { obj_dump(s, &val.x, 2); }
inline void obj_dump(char*& s, vec2f val) { obj_dump(s, &val.x, 2); }
inline void obj_dump(char*& s, vec3f val) { obj_dump(s, &val.x, 3); }
inline void obj_dump(char*& s, vec4f val) { obj_dump(s, &val.x, 4); }
inline void obj_dump(char*& s, const frame3f& val) {
    obj_dump(s, &val.x.x, 12);
}
inline void obj_dump(char*& s, const std::array<int, 3>& val) {
    auto nto_write = 0;
    for (auto i = 0; i < 3; i++) {
        if (val[i] >= 0) nto_write = i + 1;
    }
    for (auto i = 0; i < nto_write; i++) {
        if (i) *s++ = '/';
        if (val[i] >= 0) s += sprintf(s, "%d", val[i] + 1);
    }
}
inline void obj_dump(char*& s, const std::vector<std::string>& vals) {
    for (auto i = 0; i < vals.size(); i++) {
        if (i) *s++ = ' ';
        obj_dump(s, vals[i]);
    }
}
inline void obj_dump(char*& s, const obj_texture_info& v) {
    for (auto&& kv : v.props) {
        obj_dump(s, kv.first);
        *s++ = ' ';
        for (auto&& vv : kv.second) {
            obj_dump(s, vv);
            *s++ = ' ';
        }
    }
    if (v.clamp) obj_dump(s, "-clamp on ");
    obj_dump(s, v.path);
}

// Dumps a line
template <typename T>
inline void obj_dump_line(FILE* fs, const char* lbl, const T& val) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, lbl);
    *s++ = ' ';
    obj_dump(s, val);
    *s++ = '\n';
    *s = 0;
    fputs(buf, fs);
}

// Dumps a line
template <typename T>
inline void obj_dump_sp(FILE* fs, const T& val) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, val);
    *s++ = ' ';
    *s = 0;
    fputs(buf, fs);
}

// Dumps a newline
inline void obj_dump_nl(FILE* fs) { fputs("\n", fs); }

// Save an MTL file
inline void save_mtl(const std::string& filename,
    const std::vector<obj_material*>& materials, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // for each material, dump all the values
    for (auto mat : materials) {
        obj_dump_line(fs, "newmtl", mat->name);
        obj_dump_line(fs, "  illum", mat->illum);
        if (mat->ke != zero3f) obj_dump_line(fs, "  Ke", mat->ke);
        if (mat->ka != zero3f) obj_dump_line(fs, "  Ka", mat->ka);
        if (mat->kd != zero3f) obj_dump_line(fs, "  Kd", mat->kd);
        if (mat->ks != zero3f) obj_dump_line(fs, "  Ks", mat->ks);
        if (mat->kr != zero3f) obj_dump_line(fs, "  Kr", mat->kr);
        if (mat->kt != zero3f) obj_dump_line(fs, "  Kt", mat->kt);
        if (mat->ns != 0.0f)
            obj_dump_line(fs, "  Ns", (int)clamp(mat->ns, 0.0f, 1000000000.0f));
        if (mat->op != 1.0f) obj_dump_line(fs, "  d", mat->op);
        if (mat->ior != 1.0f) obj_dump_line(fs, "  Ni", mat->ior);
        if (mat->km != -1.0f) obj_dump_line(fs, "  Pm", mat->km);
        if (mat->rs != -1.0f) obj_dump_line(fs, "  Pr", mat->rs);
        if (mat->ke_txt.path != "") obj_dump_line(fs, "  map_Ke", mat->ke_txt);
        if (mat->ka_txt.path != "") obj_dump_line(fs, "  map_Ka", mat->ka_txt);
        if (mat->kd_txt.path != "") obj_dump_line(fs, "  map_Kd", mat->kd_txt);
        if (mat->ks_txt.path != "") obj_dump_line(fs, "  map_Ks", mat->ks_txt);
        if (mat->kr_txt.path != "") obj_dump_line(fs, "  map_Kr", mat->kr_txt);
        if (mat->kt_txt.path != "") obj_dump_line(fs, "  map_Kt", mat->kt_txt);
        if (mat->ns_txt.path != "") obj_dump_line(fs, "  map_Ns", mat->ns_txt);
        if (mat->op_txt.path != "") obj_dump_line(fs, "  map_d ", mat->op_txt);
        if (mat->ior_txt.path != "")
            obj_dump_line(fs, "  map_Ni", mat->ior_txt);
        if (mat->km_txt.path != "") obj_dump_line(fs, "  map_Pm ", mat->km_txt);
        if (mat->rs_txt.path != "") obj_dump_line(fs, "  map_Pr ", mat->rs_txt);
        if (mat->occ_txt.path != "")
            obj_dump_line(fs, "  map_occ ", mat->occ_txt);
        if (mat->bump_txt.path != "")
            obj_dump_line(fs, "  map_bump", mat->bump_txt);
        if (mat->disp_txt.path != "")
            obj_dump_line(fs, "  map_disp", mat->disp_txt);
        if (mat->norm_txt.path != "")
            obj_dump_line(fs, "  map_norm", mat->norm_txt);
        for (auto&& kv : mat->props) {
            obj_dump_line(fs, ("  " + kv.first + ' ').c_str(), kv.second);
        }
        obj_dump_nl(fs);
    }

    fclose(fs);
}

// Save an OBJ
inline void save_obj(const std::string& filename, const obj_scene*& obj,
    bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // linkup to mtl
    auto mtlname = replace_path_extension(filename, "mtl");
    if (!obj->materials.empty()) { obj_dump_line(fs, "mtllib", mtlname); }

    // save cameras
    for (auto cam : obj->cameras) {
        obj_dump_sp(fs, "c");
        obj_dump_sp(fs, cam->name);
        obj_dump_sp(fs, cam->ortho);
        obj_dump_sp(fs, cam->width);
        obj_dump_sp(fs, cam->height);
        obj_dump_sp(fs, cam->focal);
        obj_dump_sp(fs, cam->focus);
        obj_dump_sp(fs, cam->aperture);
        obj_dump_sp(fs, cam->frame);
        obj_dump_nl(fs);
    }

    // save envs
    for (auto env : obj->environments) {
        obj_dump_sp(fs, "e");
        obj_dump_sp(fs, env->name);
        obj_dump_sp(fs, env->ke);
        obj_dump_sp(fs, (env->ke_txt.path != "") ? env->ke_txt.path : "\"\"");
        obj_dump_sp(fs, env->frame);
        obj_dump_nl(fs);
    }

    // save object properties
    for (auto oobj : obj->objects) {
        for (auto& kv : oobj->props) {
            obj_dump_sp(fs, "op");
            obj_dump_sp(fs, kv.first);
            for (auto v : kv.second) obj_dump_sp(fs, v);
            obj_dump_nl(fs);
        }
    }

    // save all vertex data
    for (auto& v : obj->pos) obj_dump_line(fs, "v", v);
    if (flip_texcoord) {
        for (auto& v : obj->texcoord)
            obj_dump_line(fs, "vt", vec2f{v.x, 1 - v.y});
    } else {
        for (auto& v : obj->texcoord) obj_dump_line(fs, "vt", v);
    }
    for (auto& v : obj->norm) obj_dump_line(fs, "vn", v);

    // save element data
    static auto elem_labels = std::unordered_map<obj_element_type, std::string>{
        {obj_element_type::point, "p"}, {obj_element_type::line, "l"},
        {obj_element_type::face, "f"}, {obj_element_type::bezier, "b"}};
    for (auto oobj : obj->objects) {
        obj_dump_line(fs, "o", oobj->name);
        if (oobj->frame != identity_frame3f)
            obj_dump_line(fs, "of", oobj->frame);
        if (oobj->subdiv != zero2i) obj_dump_line(fs, "os", oobj->subdiv);
        auto last_groupid = -1, last_materialid = -1, last_smoothingid = -1;
        for (auto& elem : oobj->elems) {
            if (last_materialid != elem.material && !oobj->materials.empty()) {
                auto matname = oobj->materials[elem.material];
                if (matname != "") obj_dump_line(fs, "usemtl", matname);
                last_materialid = elem.material;
            }
            if (last_groupid != elem.group && !oobj->groups.empty()) {
                auto groupname = oobj->groups[elem.group];
                if (groupname != "" || elem.group > 0)
                    obj_dump_line(fs, "g", groupname);
                last_groupid = elem.group;
            }
            if (last_smoothingid != elem.smoothing &&
                !oobj->smoothing.empty()) {
                auto smoothing = oobj->smoothing[elem.smoothing];
                if (!smoothing || elem.smoothing > 0)
                    obj_dump_line(fs, "s", smoothing ? "on" : "off");
                last_smoothingid = elem.smoothing;
            }
            obj_dump_sp(fs, elem_labels.at(elem.type).c_str());
            for (auto vid = elem.start; vid < elem.start + elem.size; vid++) {
                auto vert = std::array<int, 3>{{-1, -1, -1}};
                vert[0] = (oobj->verts_pos.empty()) ? -1 : oobj->verts_pos[vid];
                vert[1] = (oobj->verts_texcoord.empty()) ?
                              -1 :
                              oobj->verts_texcoord[vid];
                vert[2] =
                    (oobj->verts_norm.empty()) ? -1 : oobj->verts_norm[vid];
                obj_dump_sp(fs, vert);
            }
            obj_dump_nl(fs);
        }
    }

    fclose(fs);

    // save materials
    if (!obj->materials.empty())
        save_mtl(
            replace_path_extension(filename, "mtl"), obj->materials, flip_tr);
}

// Load OBJ texture images.
inline void load_obj_textures(
    const obj_scene*& obj, const std::string& dirname, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<std::string, float>{{"", 1.0f}};
    for (auto txt : obj->textures) ldr_gamma[txt->path] = 2.2f;
    for (auto mat : obj->materials) {
        ldr_gamma[mat->ke_txt.path] = 2.2f;
        ldr_gamma[mat->kd_txt.path] = 2.2f;
        ldr_gamma[mat->ks_txt.path] = 2.2f;
        ldr_gamma[mat->kr_txt.path] = 2.2f;
        ldr_gamma[mat->kt_txt.path] = 2.2f;
        ldr_gamma[mat->rs_txt.path] = 1;
        ldr_gamma[mat->op_txt.path] = 1;
        ldr_gamma[mat->norm_txt.path] = 1;
        ldr_gamma[mat->disp_txt.path] = 1;
        ldr_gamma[mat->bump_txt.path] = 1;
    }
    for (auto env : obj->environments) { ldr_gamma[env->ke_txt.path] = 2.2f; }

    // load images
    for (auto txt : obj->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            txt->img = load_image4f(filename, ldr_gamma.at(txt->path));
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Save OBJ texture images.
void save_obj_textures(
    const obj_scene*& obj, const std::string& dirname, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<std::string, float>{{"", 1.0f}};
    for (auto txt : obj->textures) ldr_gamma[txt->path] = 2.2f;
    for (auto mat : obj->materials) {
        ldr_gamma[mat->ke_txt.path] = 2.2f;
        ldr_gamma[mat->kd_txt.path] = 2.2f;
        ldr_gamma[mat->ks_txt.path] = 2.2f;
        ldr_gamma[mat->kr_txt.path] = 2.2f;
        ldr_gamma[mat->kt_txt.path] = 2.2f;
        ldr_gamma[mat->rs_txt.path] = 1;
        ldr_gamma[mat->op_txt.path] = 1;
        ldr_gamma[mat->norm_txt.path] = 1;
        ldr_gamma[mat->disp_txt.path] = 1;
        ldr_gamma[mat->bump_txt.path] = 1;
    }
    for (auto env : obj->environments) { ldr_gamma[env->ke_txt.path] = 2.2f; }

    // save images
    for (auto txt : obj->textures) {
        if (txt->empty(img)) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_image4f(filename, txt->img, ldr_gamma.at(txt->path));
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

}  // namespace detail

// Load an OBJ from file `filename`. Split shapes at material and group
// boundaries, if `split_shapes` is true.
// Load textures if `load_textures` is true, and report errors only if
// `skip_missing` is false. Texture coordinates and material Tr are flipped
// if `flip_texcoord` and `flip_tp` are respectively true.
inline obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool flip_texcoord, bool flip_tr) {
    return detail::load_obj(filename, split_shapes, flip_texcoord, flip_tr);
}

// Save an OBJ to file `filename`. Save textures if `save_textures` is true,
// and report errors only if `skip_missing` is false.
// Texture coordinates and material Tr are flipped if `flip_texcoord` and
// `flip_tp` are respectively true.
inline void save_obj(const std::string& filename, const obj_scene*& obj,
    bool flip_texcoord, bool flip_tr) {
    return detail::save_obj(filename, obj, flip_texcoord, flip_tr);
}

// Load OBJ texture images.
inline void load_obj_textures(
    const obj_scene*& obj, const std::string& dirname, bool skip_missing) {
    return detail::load_obj_textures(obj, dirname, skip_missing);
}
// Save OBJ texture images.
inline void save_obj_textures(
    const obj_scene*& obj, const std::string& dirname, bool skip_missing) {
    return detail::save_obj_textures(obj, dirname, skip_missing);
}

// Load OBJ with callbacks
inline void load_obj(
    const std::string& filename, const obj_callbacks& callbacks, void* ctx) {
    return detail::load_obj(filename, callbacks, ctx);
}

}  // namespace ygl

#endif
