//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
// TODO: add tangent space
//

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR YOCTO_OBJ
// -----------------------------------------------------------------------------

#include "yocto_obj.h"

#include <algorithm>
#include <cassert>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <unordered_map>

#ifndef YOBJ_NO_IMAGE
#include "yocto_img.h"
#endif

namespace yobj {

//
// Get extension (including '.').
//
inline std::string get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Get directory name (including '/').
//
inline std::string get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Splits a std::string into an array of strings on whitespace with Python split
// semantic. Modifies original std::string to avoid allocation.
//
inline int splitws(char* str, char** splits, int maxsplits) {
    int n = 0;
    while (*str && n < maxsplits) {
        if (isspace(*str)) {
            *str = 0;
        } else {
            if (n == 0 || !(*(str - 1))) {
                splits[n] = str;
                n++;
            }
        }
        str++;
    }
    splits[n] = nullptr;
    return n;
}

//
// Parses one int.
//
inline int parse_int(char** tok) { return atoi(tok[0]); }

//
// Parses one float.
//
inline float parse_float(char** tok) { return atof(tok[0]); }

//
// Parses two floats.
//
inline ym::vec2f parse_float2(char** tok) {
    return ym::vec2f{(float)atof(tok[0]), (float)atof(tok[1])};
}

//
// Parses three floats.
//
inline ym::vec3f parse_float3(char** tok) {
    return ym::vec3f{
        (float)atof(tok[0]), (float)atof(tok[1]), (float)atof(tok[2])};
}

//
// Parses four floats.
//
inline ym::vec4f parse_float4(char** tok) {
    return ym::vec4f{(float)atof(tok[0]), (float)atof(tok[1]),
        (float)atof(tok[2]), (float)atof(tok[3])};
}

//
// Parses 16 floats.
//
inline ym::mat4f parse_float16(char** tok) {
    ym::mat4f m;
    auto mm = (float*)&m;
    for (auto i = 0; i < 16; i++) mm[i] = (float)atof(tok[i]);
    return m;
}

//
// Parses an OBJ vertex list. Handles negative values.
//
inline void parse_vertlist(char** tok, int ntoks,
    std::vector<obj_vertex>& elems, const obj_vertex& vert_size) {
    elems.clear();
    for (auto i = 0; i < ntoks; i++) {
        // parse triplet
        char* splits[] = {tok[i], 0, 0, 0, 0};
        auto ns = 1;
        while (*tok[i]) {
            if (*tok[i] == '/') {
                *tok[i] = 0;
                if (ns < 5) splits[ns++] = tok[i] + 1;
            }
            tok[i]++;
        }
        auto v = obj_vertex{-1, -1, -1, -1, -1};
        auto v_ptr = &v.pos;
        auto vs_ptr = &vert_size.pos;
        for (auto i = 0; i < 5; i++) {
            if (!splits[i]) {
                v_ptr[i] = -1;
                continue;
            }
            v_ptr[i] = (int)atoi(splits[i]);
            v_ptr[i] = (v_ptr[i] < 0) ? vs_ptr[i] + v_ptr[i] : v_ptr[i] - 1;
        }
        elems.push_back(v);
    }
}

//
// Loads an OBJ
//
obj* load_obj(const std::string& filename, bool flip_texcoord, bool flip_tr,
    std::string* err) {
    // clear obj
    auto asset = std::unique_ptr<obj>(new obj());

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) {
        if (err) *err = "cannot open filename " + filename;
        return nullptr;
    }

    // initializing obj
    asset->objects.push_back({});
    asset->objects.back().groups.push_back({});

    // allocate buffers to avoid re-allocing
    auto cur_elems = std::vector<obj_vertex>();
    auto cur_matname = std::string();
    auto cur_mtllibs = std::vector<std::string>();

    // keep track of array lengths
    auto vert_size = obj_vertex{0, 0, 0, 0, 0};

    // read the file line by line
    char line[4096];
    char* toks[1024];
    auto linenum = 0;
    while (fgets(line, 4096, file)) {
        linenum += 1;
        int ntok = splitws(line, toks, 1024);

        // skip empty and comments
        if (!ntok) continue;
        if (toks[0][0] == '#') continue;

        // set up code
        auto tok_s = std::string(toks[0]);
        auto cur_tok = toks + 1;
        auto cur_ntok = ntok - 1;

        // possible token values
        if (tok_s == "v") {
            vert_size.pos += 1;
            asset->pos.push_back(parse_float3(cur_tok));
        } else if (tok_s == "vn") {
            vert_size.norm += 1;
            asset->norm.push_back(parse_float3(cur_tok));
        } else if (tok_s == "vt") {
            vert_size.texcoord += 1;
            asset->texcoord.push_back(parse_float2(cur_tok));
            if (flip_texcoord)
                asset->texcoord.back()[1] = 1 - asset->texcoord.back()[1];
        } else if (tok_s == "vc") {
            vert_size.color += 1;
            asset->color.push_back(parse_float4(cur_tok));
        } else if (tok_s == "vr") {
            vert_size.radius += 1;
            asset->radius.push_back(parse_float(cur_tok));
        } else if (tok_s == "f") {
            parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::face,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "l") {
            parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::line,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "p") {
            parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::point, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "t") {
            parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::tetra, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "o") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.push_back({name, {}});
            asset->objects.back().groups.push_back({cur_matname, ""});
        } else if (tok_s == "usemtl") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            cur_matname = name;
            asset->objects.back().groups.push_back({cur_matname, ""});
        } else if (tok_s == "g") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.back().groups.push_back({cur_matname, name});
        } else if (tok_s == "s") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            auto smoothing = name == std::string("on");
            if (asset->objects.back().groups.back().smoothing != smoothing) {
                asset->objects.back().groups.push_back(
                    {cur_matname, name, smoothing});
            }
        } else if (tok_s == "mtllib") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            if (name != std::string("")) {
                auto found = false;
                for (auto lib : cur_mtllibs) {
                    if (lib == name) {
                        found = true;
                        break;
                    }
                }
                if (!found) cur_mtllibs.push_back(name);
            }
        } else if (tok_s == "c") {
            asset->cameras.emplace_back();
            auto& cam = asset->cameras.back();
            cam.name = (cur_ntok) ? cur_tok[0] : "";
            cam.ortho = parse_int(cur_tok + 1);
            cam.yfov = parse_float(cur_tok + 2);
            cam.aspect = parse_float(cur_tok + 3);
            cam.aperture = parse_float(cur_tok + 4);
            cam.focus = parse_float(cur_tok + 5);
            cam.translation = parse_float3(cur_tok + 6);
            cam.rotation = (ym::quat4f)parse_float4(cur_tok + 9);
            if (cur_ntok > 13) cam.matrix = parse_float16(cur_tok + 13);
        } else if (tok_s == "e") {
            asset->environments.emplace_back();
            auto& env = asset->environments.back();
            env.name = (cur_ntok) ? cur_tok[0] : "<unnamed>";
            env.matname = (cur_ntok - 1) ? cur_tok[1] : "<unnamed_material>";
            env.rotation = (ym::quat4f)parse_float4(cur_tok + 2);
            if (cur_ntok > 6) env.matrix = parse_float16(cur_tok + 6);
        } else if (tok_s == "i") {
            asset->instances.emplace_back();
            auto& ist = asset->instances.back();
            ist.name = (cur_ntok) ? cur_tok[0] : "<unnamed>";
            ist.meshname = (cur_ntok - 1) ? cur_tok[1] : "<unnamed_mesh>";
            ist.translation = parse_float3(cur_tok + 2);
            ist.rotation = (ym::quat4f)parse_float4(cur_tok + 5);
            ist.scale = parse_float3(cur_tok + 9);
            if (cur_ntok > 12) ist.matrix = parse_float16(cur_tok + 12);
        } else {
            // unused
        }
    }

    // cleanup unused
    for (auto&& o : asset->objects) {
        auto end = std::remove_if(o.groups.begin(), o.groups.end(),
            [](const obj_group& x) { return x.verts.empty(); });
        o.groups.erase(end, o.groups.end());
    }
    auto end = std::remove_if(asset->objects.begin(), asset->objects.end(),
        [](const obj_object& x) { return x.groups.empty(); });
    asset->objects.erase(end, asset->objects.end());

    // parse materials
    for (auto mtllib : cur_mtllibs) {
        auto mtlname = get_dirname(filename) + mtllib;
        std::string errm;
        auto materials = load_mtl(mtlname, flip_tr, &errm);
        if (materials.empty() && !errm.empty()) {
            if (err) *err = errm;
            return nullptr;
        }
        asset->materials.insert(
            asset->materials.end(), materials.begin(), materials.end());
    }

    // done
    return asset.release();
}

//
// Parse texture options and name
//
void parse_texture(char** toks, int ntoks, std::string& path,
    property_map<std::string>& info) {
    // texture name
    if (ntoks > 0) {
        path = toks[ntoks - 1];
        for (auto& c : path)
            if (c == '\\') c = '/';
    }

    // texture options
    if (ntoks > 1) {
        auto cur_ntoks = ntoks - 1;
        auto cur_tok = toks;
        while (cur_ntoks) {
            if (cur_tok[0][0] != '-') break;
            auto name = cur_tok[0];
            info[name] = {};
            cur_ntoks--;
            cur_tok++;
            while (cur_ntoks && cur_tok[0][0] != '-') {
                info[name].push_back(cur_tok[0]);
                cur_ntoks--;
                cur_tok++;
            }
        }
    }
}

//
// Load MTL
//
std::vector<obj_material> load_mtl(
    const std::string& filename, bool flip_tr, std::string* err) {
    // clear materials
    auto materials = std::vector<obj_material>();

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) {
        if (err) *err = "cannot open filename " + filename;
        return {};
    }

    // add a material preemptively to avoid crashes
    materials.emplace_back();

    // read the file line by line
    char line[4096];
    char* toks[1024];
    auto linenum = 0;
    while (fgets(line, 4096, file)) {
        linenum += 1;
        int ntok = splitws(line, toks, 1024);

        // skip empty and comments
        if (!ntok) continue;
        if (toks[0][0] == '#') continue;

        // set up code
        auto tok_s = std::string(toks[0]);
        auto cur_tok = toks + 1;
        auto cur_ntok = ntok - 1;

        // possible token values
        if (tok_s == "newmtl") {
            materials.emplace_back();
            materials.back().name = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "illum") {
            materials.back().illum = parse_int(cur_tok);
        } else if (tok_s == "Ke") {
            materials.back().ke = parse_float3(cur_tok);
        } else if (tok_s == "Ka") {
            materials.back().ka = parse_float3(cur_tok);
        } else if (tok_s == "Kd") {
            materials.back().kd = parse_float3(cur_tok);
        } else if (tok_s == "Ks") {
            materials.back().ks = parse_float3(cur_tok);
        } else if (tok_s == "Kr") {
            materials.back().kr = parse_float3(cur_tok);
        } else if (tok_s == "Kt" || tok_s == "Tf") {
            if (cur_ntok >= 3) {
                materials.back().kt = parse_float3(cur_tok);
            } else {
                auto v = parse_float(cur_tok);
                materials.back().kt = {v, v, v};
            }
        } else if (tok_s == "Tr") {
            if (cur_ntok >= 3) {
                materials.back().kt = parse_float3(cur_tok);
            } else {
                // as tinyobjreader
                if (flip_tr)
                    materials.back().op = 1 - parse_float(cur_tok);
                else
                    materials.back().op = parse_float(cur_tok);
            }
        } else if (tok_s == "Ns") {
            materials.back().ns = parse_float(cur_tok);
        } else if (tok_s == "d") {
            materials.back().op = parse_float(cur_tok);
        } else if (tok_s == "Ni") {
            materials.back().ior = parse_float(cur_tok);
        } else if (tok_s == "map_Ke") {
            parse_texture(cur_tok, cur_ntok, materials.back().ke_txt,
                materials.back().ke_txt_info);
        } else if (tok_s == "map_Ka") {
            parse_texture(cur_tok, cur_ntok, materials.back().ka_txt,
                materials.back().ka_txt_info);
        } else if (tok_s == "map_Kd") {
            parse_texture(cur_tok, cur_ntok, materials.back().kd_txt,
                materials.back().kd_txt_info);
        } else if (tok_s == "map_Ks") {
            parse_texture(cur_tok, cur_ntok, materials.back().ks_txt,
                materials.back().ks_txt_info);
        } else if (tok_s == "map_Kr") {
            parse_texture(cur_tok, cur_ntok, materials.back().kr_txt,
                materials.back().kr_txt_info);
        } else if (tok_s == "map_Tr") {
            parse_texture(cur_tok, cur_ntok, materials.back().kt_txt,
                materials.back().kt_txt_info);
        } else if (tok_s == "map_Ns") {
            parse_texture(cur_tok, cur_ntok, materials.back().ns_txt,
                materials.back().ns_txt_info);
        } else if (tok_s == "map_d") {
            parse_texture(cur_tok, cur_ntok, materials.back().op_txt,
                materials.back().op_txt_info);
        } else if (tok_s == "map_Ni") {
            parse_texture(cur_tok, cur_ntok, materials.back().ior_txt,
                materials.back().ior_txt_info);
        } else if (tok_s == "map_bump" || tok_s == "bump") {
            parse_texture(cur_tok, cur_ntok, materials.back().bump_txt,
                materials.back().bump_txt_info);
        } else if (tok_s == "map_disp" || tok_s == "disp") {
            parse_texture(cur_tok, cur_ntok, materials.back().disp_txt,
                materials.back().disp_txt_info);
        } else if (tok_s == "map_norm" || tok_s == "norm") {
            parse_texture(cur_tok, cur_ntok, materials.back().norm_txt,
                materials.back().norm_txt_info);
        } else {
            // copy into strings
            for (auto i = 0; i < cur_ntok; i++)
                materials.back().unknown_props[tok_s].push_back(cur_tok[i]);
        }
    }

    // remove first fake material
    materials.erase(materials.begin());

    // done
    return materials;
}

//
// write one float prepended by a std::string
//
inline void fwrite_float(
    FILE* file, const char* str, float v, bool newline = true) {
    fprintf(file, "%s %.6g", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write one float prepended by a std::string
//
inline void fwrite_int(
    FILE* file, const char* str, int v, bool newline = true) {
    fprintf(file, "%s %d", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write two floats prepended by a std::string
//
inline void fwrite_float2(
    FILE* file, const char* str, const ym::vec2f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g", str, v[0], v[1]);
    if (newline) fprintf(file, "\n");
}

//
// write three floats prepended by a std::string
//
inline void fwrite_float3(
    FILE* file, const char* str, const ym::vec3f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g %.6g", str, v[0], v[1], v[2]);
    if (newline) fprintf(file, "\n");
}

//
// write four floats prepended by a std::string
//
inline void fwrite_float4(
    FILE* file, const char* str, const ym::vec4f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g %.6g %.6g", str, v[0], v[1], v[2], v[3]);
    if (newline) fprintf(file, "\n");
}

//
// write 16 floats prepended by a std::string
//
inline void fwrite_float16(
    FILE* file, const char* str, const ym::mat4f& v, bool newline = true) {
    const float* vf = (float*)&v;
    fprintf(file, "%s", str);
    for (int i = 0; i < 16; i++) fprintf(file, " %.6g", vf[i]);
    if (newline) fprintf(file, "\n");
}

//
// write a std::string prepended by another if the std::string is not NULL
//
inline void fwrite_str(FILE* file, const char* str, const std::string& s,
    bool force = false, bool newline = true) {
    if (s.empty() && !force) return;
    fprintf(file, "%s %s", str, s.c_str());
    if (newline) fprintf(file, "\n");
}

//
// write a std::string prepended by another if the std::string is not NULL
//
inline void fwrite_str_props(FILE* file, const char* str, const std::string& s,
    const property_map<std::string>& props, bool force = false,
    bool newline = true) {
    if (s.empty() && !force) return;
    if (props.empty()) {
        fprintf(file, "%s %s", str, s.c_str());
    } else {
        auto props_str = std::string();
        for (auto&& prop : props) {
            props_str += prop.first + " ";
            for (auto&& pp : prop.second) props_str += pp + " ";
        }
        fprintf(file, "%s %s %s", str, props_str.c_str(), s.c_str());
    }
    if (newline) fprintf(file, "\n");
}

//
// write one float prepended by a std::string
//
inline void fwrite_float_props(
    FILE* file, const char* str, float v, float def = 0, bool newline = true) {
    if (v == def) return;
    fprintf(file, "%s %.6g", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write three floats prepended by a std::string
//
inline void fwrite_float3_props(FILE* file, const char* str, const ym::vec3f& v,
    const ym::vec3f& def = {0, 0, 0}, bool newline = true) {
    if (v == def) return;
    fprintf(file, "%s %.6g %.6g %.6g", str, v[0], v[1], v[2]);
    if (newline) fprintf(file, "\n");
}

//
// write an OBJ vertex triplet using only the indices that are active
//
inline void fwrite_objverts(FILE* file, const char* str, int nv,
    const obj_vertex* verts, bool newline = true) {
    fprintf(file, "%s", str);
    for (auto v = 0; v < nv; v++) {
        auto vert = verts[v];
        auto vert_ptr = &vert.pos;
        auto nto_write = 0;
        for (auto i = 0; i < 5; i++) {
            if (vert_ptr[i] >= 0) nto_write = i + 1;
        }
        for (auto i = 0; i < nto_write; i++) {
            if (vert_ptr[i] >= 0) {
                fprintf(file, "%c%d", ((i == 0) ? ' ' : '/'), vert_ptr[i] + 1);
            } else {
                fprintf(file, "%c", '/');
            }
        }
    }
    fprintf(file, "\n");
}

//
// Save an OBJ
//
bool save_obj(const std::string& filename, const obj* asset, bool flip_texcoord,
    bool flip_tr, std::string* err) {
    // open file
    auto file = fopen(filename.c_str(), "wt");
    if (!file) {
        if (err) *err = "could not open filename " + filename;
        return false;
    }

    // linkup to mtl
    auto dirname = get_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!asset->materials.empty()) {
        fwrite_str(file, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto& cam : asset->cameras) {
        fwrite_str(file, "c", cam.name, true, false);
        fwrite_int(file, " ", cam.ortho, false);
        fwrite_float(file, " ", cam.yfov, false);
        fwrite_float(file, " ", cam.aspect, false);
        fwrite_float(file, " ", cam.aperture, false);
        fwrite_float(file, " ", cam.focus, false);
        fwrite_float3(file, " ", cam.translation, false);
        fwrite_float4(file, " ", (ym::vec4f)cam.rotation, false);
        if (cam.matrix != ym::identity_mat4f)
            fwrite_float16(file, " ", cam.matrix, true);
    }

    // save envs
    for (auto& env : asset->environments) {
        fwrite_str(file, "e", env.name, true, false);
        fwrite_str(file, " ", env.matname, true, false);
        fwrite_float4(file, " ", (ym::vec4f)env.rotation, false);
        if (env.matrix != ym::identity_mat4f)
            fwrite_float16(file, " ", env.matrix, true);
    }

    // save instances
    for (auto& ist : asset->instances) {
        fwrite_str(file, "i", ist.name, true, false);
        fwrite_str(file, " ", ist.meshname, true, false);
        fwrite_float3(file, " ", ist.translation, false);
        fwrite_float4(file, " ", (ym::vec4f)ist.rotation, false);
        fwrite_float3(file, " ", ist.scale, false);
        if (ist.matrix != ym::identity_mat4f)
            fwrite_float16(file, " ", ist.matrix, true);
        else
            fprintf(file, "\n");
    }

    // save all vertex data
    for (auto& v : asset->pos) fwrite_float3(file, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            fwrite_float2(file, "vt", {v[0], 1 - v[1]});
    } else {
        for (auto& v : asset->texcoord) fwrite_float2(file, "vt", v);
    }
    for (auto& v : asset->norm) fwrite_float3(file, "vn", v);
    for (auto& v : asset->color) fwrite_float4(file, "vc", v);
    for (auto& v : asset->radius) fwrite_float(file, "vr", v);

    // save element data
    const char* elem_labels[] = {"", "p", "l", "f", "t"};
    for (auto& object : asset->objects) {
        fwrite_str(file, "o", object.name, true);
        for (auto& group : object.groups) {
            fwrite_str(file, "usemtl", group.matname);
            fwrite_str(file, "g", group.groupname);
            if (!group.smoothing) fwrite_str(file, "s", "off");
            for (auto elem : group.elems) {
                fwrite_objverts(file, elem_labels[(int)elem.type], elem.size,
                    group.verts.data() + elem.start);
            }
        }
    }

    fclose(file);

    // save materials
    if (!asset->materials.empty()) {
        if (!save_mtl(
                dirname + basename + ".mtl", asset->materials, flip_tr, err))
            return false;
    }

    // done
    return true;
}

//
// Save an MTL file
//
bool save_mtl(const std::string& filename,
    const std::vector<obj_material>& materials, bool flip_tr,
    std::string* err) {
    auto file = fopen(filename.c_str(), "wt");
    if (!file) {
        if (err) *err = "could not open filename " + filename;
        return false;
    }

    // for each material, dump all the values
    for (auto& mat : materials) {
        fwrite_str(file, "newmtl", mat.name, true);
        fwrite_int(file, "  illum", mat.illum);
        fwrite_float3_props(file, "  Ke", mat.ke);
        fwrite_float3_props(file, "  Ka", mat.ka);
        fwrite_float3_props(file, "  Kd", mat.kd);
        fwrite_float3_props(file, "  Ks", mat.ks);
        fwrite_float3_props(file, "  Kr", mat.kr);
        fwrite_float3_props(file, "  Tf", mat.kt);
        fwrite_float_props(file, "  Ns", mat.ns, 0);
        fwrite_float_props(file, "  d", mat.op, 1);
        fwrite_float_props(file, "  Ni", mat.ior, 1);
        fwrite_str_props(file, "  map_Ke", mat.ke_txt, mat.ke_txt_info);
        fwrite_str_props(file, "  map_Ka", mat.ka_txt, mat.ka_txt_info);
        fwrite_str_props(file, "  map_Kd", mat.kd_txt, mat.kd_txt_info);
        fwrite_str_props(file, "  map_Ks", mat.ks_txt, mat.ks_txt_info);
        fwrite_str_props(file, "  map_Kr", mat.kr_txt, mat.kr_txt_info);
        fwrite_str_props(file, "  map_Kt", mat.kt_txt, mat.kt_txt_info);
        fwrite_str_props(file, "  map_Ns", mat.ns_txt, mat.ns_txt_info);
        fwrite_str_props(file, "  map_d", mat.op_txt, mat.op_txt_info);
        fwrite_str_props(file, "  map_Ni", mat.ior_txt, mat.ior_txt_info);
        fwrite_str_props(file, "  map_bump", mat.bump_txt, mat.bump_txt_info);
        fwrite_str_props(file, "  map_disp", mat.disp_txt, mat.disp_txt_info);
        fwrite_str_props(file, "  map_norm", mat.norm_txt, mat.norm_txt_info);
        for (auto&& p : mat.unknown_props) {
            auto s = std::string();
            for (auto&& v : p.second) s += v + " ";
            fwrite_str(file, p.first.c_str(), s.c_str());
        }
        fprintf(file, "\n");
    }

    fclose(file);

    // done
    return true;
}

//
// cleanup
//
mesh::~mesh() {
    for (auto v : shapes)
        if (v) delete v;
}

//
// Cleanup memory
//
scene::~scene() {
    for (auto v : meshes)
        if (v) delete v;
    for (auto v : instances)
        if (v) delete v;
    for (auto v : materials)
        if (v) delete v;
    for (auto v : textures)
        if (v) delete v;
    for (auto v : cameras)
        if (v) delete v;
    for (auto v : environments)
        if (v) delete v;
}

//
// Loads a textures and saves into an array
//
static texture* add_texture(
    const std::string& filename, std::vector<texture*>& txts) {
    if (filename.empty()) return nullptr;
    for (auto txt : txts) {
        if (txt->path == filename) return txt;
    }
    auto txt = new texture();
    txt->path = filename;
    txts.push_back(txt);
    return txt;
}

//
// A hash function for vecs
//
struct vertex_hash {
    std::hash<int> Th;
    size_t operator()(const obj_vertex& vv) const {
        auto v = (const int*)&vv;
        size_t h = 0;
        for (auto i = 0; i < sizeof(obj_vertex) / sizeof(int); i++) {
            // embads hash_combine below
            h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
        }
        return h;
    }
};

//
// Comparison for unordred_map
//
static bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm &&
           a.color == b.color && a.radius == b.radius;
}

//
// Converts texture info
//
texture_info _convert_texture_info(const property_map<std::string>& props) {
    auto info = texture_info();
    if (props.find("-clamp") != props.end() && !props.at("-clamp").empty()) {
        info.clamp =
            props.at("-clamp")[0] == "on" || props.at("-clamp")[0] == "1";
    }
    if (props.find("-bm") != props.end() && !props.at("-bm").empty()) {
        info.clamp = std::atof(props.at("-bm")[0].c_str());
    }
    return info;
}

//
// Flattens an scene
//
scene* obj_to_scene(const obj* asset, bool facet_non_smooth) {
    // clear scene
    auto scn = new scene();

    // convert materials and build textures
    for (auto& omat : asset->materials) {
        auto mat = new material();
        mat->name = omat.name;
        mat->ke = omat.ke;
        mat->kd = omat.kd;
        mat->ks = omat.ks;
        mat->kt = omat.kt;
        mat->rs = sqrt(2 / (omat.ns + 2));
        mat->opacity = omat.op;
        mat->ke_txt = add_texture(omat.ke_txt, scn->textures);
        mat->kd_txt = add_texture(omat.kd_txt, scn->textures);
        mat->ks_txt = add_texture(omat.ks_txt, scn->textures);
        mat->kt_txt = add_texture(omat.kt_txt, scn->textures);
        mat->rs_txt = add_texture(omat.ns_txt, scn->textures);
        mat->norm_txt = add_texture(omat.norm_txt, scn->textures);
        mat->bump_txt = add_texture(omat.bump_txt, scn->textures);
        mat->disp_txt = add_texture(omat.disp_txt, scn->textures);
        mat->ke_txt_info = _convert_texture_info(omat.ke_txt_info);
        mat->kd_txt_info = _convert_texture_info(omat.kd_txt_info);
        mat->ks_txt_info = _convert_texture_info(omat.ks_txt_info);
        mat->kt_txt_info = _convert_texture_info(omat.kt_txt_info);
        mat->rs_txt_info = _convert_texture_info(omat.ns_txt_info);
        mat->norm_txt_info = _convert_texture_info(omat.norm_txt_info);
        mat->bump_txt_info = _convert_texture_info(omat.bump_txt_info);
        mat->disp_txt_info = _convert_texture_info(omat.disp_txt_info);
        mat->unknown_props = omat.unknown_props;
        switch (omat.illum) {
            case 0:  // Color on and Ambient off
            case 1:  // Color on and Ambient on
            case 2:  // Highlight on
            case 3:  // Reflection on and Ray trace on
                mat->opacity = 1;
                mat->kt = {0, 0, 0};
                break;
            case 4:  // Transparency: Glass on
                     // Reflection: Ray trace on
                break;
            case 5:  // Reflection: Fresnel on and Ray trace on
                mat->opacity = 1;
                mat->kt = {0, 0, 0};
                break;
            case 6:  // Transparency: Refraction on
                     // Reflection: Fresnel off and Ray trace on
            case 7:  // Transparency: Refraction on
                     // Reflection: Fresnel on and Ray trace on
                break;
            case 8:  // Reflection on and Ray trace off
                mat->opacity = 1;
                mat->kt = {0, 0, 0};
                break;
            case 9:  // Transparency: Glass on
                     // Reflection: Ray trace off
                break;
        }
        scn->materials.push_back(mat);
    }

    // convert meshes
    for (auto& oshape : asset->objects) {
        auto msh = new mesh();
        msh->name = oshape.name;
        for (auto& group : oshape.groups) {
            if (group.verts.empty()) continue;
            if (group.elems.empty()) continue;
            auto prim = new shape();
            prim->name = group.groupname;
            prim->mat = nullptr;
            for (auto mat : scn->materials) {
                if (mat->name == group.matname) {
                    prim->mat = mat;
                    break;
                }
            }

            // insert all vertices
            std::unordered_map<obj_vertex, int, vertex_hash> vert_map;
            std::vector<int> vert_ids;
            // vert_map.clear();
            // vert_ids.clear();
            for (auto& vert : group.verts) {
                if (vert_map.find(vert) == vert_map.end()) {
                    vert_map[vert] = (int)vert_map.size();
                }
                vert_ids.push_back(vert_map.at(vert));
            }

            // convert elements
            for (auto& elem : group.elems) {
                switch (elem.type) {
                    case obj_element_type::point: {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i++) {
                            prim->points.push_back(vert_ids[i]);
                        }
                    } break;
                    case obj_element_type::line: {
                        for (auto i = elem.start;
                             i < elem.start + elem.size - 1; i++) {
                            prim->lines.push_back(
                                {vert_ids[i], vert_ids[i + 1]});
                        }
                    } break;
                    case obj_element_type::face: {
                        for (auto i = elem.start + 2;
                             i < elem.start + elem.size; i++) {
                            prim->triangles.push_back({vert_ids[elem.start],
                                vert_ids[i - 1], vert_ids[i]});
                        }
                    } break;
                    case obj_element_type::tetra: {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i += 4) {
                            if (i + 3 >= vert_ids.size()) continue;
                            prim->tetras.push_back(
                                {vert_ids[i], vert_ids[i + 1], vert_ids[i + 2],
                                    vert_ids[i + 3]});
                        }
                    } break;
                    default: { assert(false); }
                }
            }

            // check for errors
            // copy vertex data
            auto v = group.verts[0];
            if (v.pos >= 0) prim->pos.resize(vert_map.size());
            if (v.texcoord >= 0) prim->texcoord.resize(vert_map.size());
            if (v.norm >= 0) prim->norm.resize(vert_map.size());
            if (v.color >= 0) prim->color.resize(vert_map.size());
            if (v.radius >= 0) prim->radius.resize(vert_map.size());
            for (auto& kv : vert_map) {
                if (v.pos >= 0 && kv.first.pos >= 0) {
                    prim->pos[kv.second] = asset->pos[kv.first.pos];
                }
                if (v.texcoord >= 0 && kv.first.texcoord >= 0) {
                    prim->texcoord[kv.second] =
                        asset->texcoord[kv.first.texcoord];
                }
                if (v.norm >= 0 && kv.first.norm >= 0) {
                    prim->norm[kv.second] = asset->norm[kv.first.norm];
                }
                if (v.color >= 0 && kv.first.color >= 0) {
                    prim->color[kv.second] = asset->color[kv.first.color];
                }
                if (v.radius >= 0 && kv.first.radius >= 0) {
                    prim->radius[kv.second] = asset->radius[kv.first.radius];
                }
            }

            // fix smoothing
            if (!group.smoothing && facet_non_smooth) {
                auto faceted = new shape();
                faceted->name = prim->name;
                faceted->mat = prim->mat;
                auto pidx = std::vector<int>();
                for (auto point : prim->points) {
                    faceted->points.push_back((int)pidx.size());
                    pidx.push_back(point);
                }
                for (auto line : prim->lines) {
                    faceted->lines.push_back(
                        {(int)pidx.size() + 0, (int)pidx.size() + 1});
                    pidx.push_back(line.x);
                    pidx.push_back(line.y);
                }
                for (auto triangle : prim->triangles) {
                    faceted->triangles.push_back({(int)pidx.size() + 0,
                        (int)pidx.size() + 1, (int)pidx.size() + 2});
                    pidx.push_back(triangle.x);
                    pidx.push_back(triangle.y);
                    pidx.push_back(triangle.z);
                }
                for (auto tetra : prim->tetras) {
                    faceted->tetras.push_back(
                        {(int)pidx.size() + 0, (int)pidx.size() + 1,
                            (int)pidx.size() + 2, (int)pidx.size() + 3});
                    pidx.push_back(tetra.x);
                    pidx.push_back(tetra.y);
                    pidx.push_back(tetra.z);
                    pidx.push_back(tetra.w);
                }
                for (auto idx : pidx) {
                    if (!prim->pos.empty())
                        faceted->pos.push_back(prim->pos[idx]);
                    if (!prim->norm.empty())
                        faceted->norm.push_back(prim->norm[idx]);
                    if (!prim->texcoord.empty())
                        faceted->texcoord.push_back(prim->texcoord[idx]);
                    if (!prim->color.empty())
                        faceted->color.push_back(prim->color[idx]);
                    if (!prim->radius.empty())
                        faceted->radius.push_back(prim->radius[idx]);
                }
                delete prim;
                prim = faceted;
            }

            msh->shapes.push_back(prim);
        }
        scn->meshes.push_back(msh);
    }

    // convert cameras
    for (auto& ocam : asset->cameras) {
        auto cam = new camera();
        cam->name = ocam.name;
        cam->ortho = ocam.ortho;
        cam->yfov = ocam.yfov;
        cam->aspect = ocam.aspect;
        cam->aperture = ocam.aperture;
        cam->focus = ocam.focus;
        cam->translation = ocam.translation;
        cam->rotation = ocam.rotation;
        cam->matrix = ocam.matrix;
        scn->cameras.push_back(cam);
    }

    // convert envs
    for (auto& oenv : asset->environments) {
        auto env = new environment();
        env->name = oenv.name;
        env->mat = nullptr;
        for (auto mat : scn->materials) {
            if (mat->name == oenv.matname) { env->mat = mat; }
        }
        env->rotation = oenv.rotation;
        env->matrix = oenv.matrix;
        scn->environments.push_back(env);
    }

    // convert instances
    for (auto& oist : asset->instances) {
        auto ist = new instance();
        ist->name = oist.name;
        ist->msh = nullptr;
        for (auto mesh : scn->meshes) {
            if (mesh->name == oist.meshname) { ist->msh = mesh; }
        }
        ist->translation = oist.translation;
        ist->rotation = oist.rotation;
        ist->scale = oist.scale;
        ist->matrix = oist.matrix;
        scn->instances.push_back(ist);
    }

    // done
    return scn;
}

//
// Convert texture props
//
property_map<std::string> texture_props(const texture_info& info) {
    auto props = info.unknown_props;
    if (info.clamp) { props["-clamp"] = {"on"}; }
    if (info.bump_scale != 1) {
        props["-bm"] = {std::to_string(info.bump_scale)};
    }
    return props;
}

//
// Save an scene
//
obj* scene_to_obj(const scene* scn) {
    auto asset = new obj();

    // convert materials
    for (auto fl_mat : scn->materials) {
        asset->materials.emplace_back();
        auto mat = &asset->materials.back();
        mat->name = fl_mat->name;
        mat->ke = fl_mat->ke;
        mat->kd = fl_mat->kd;
        mat->ks = fl_mat->ks;
        mat->kt = fl_mat->kt;
        mat->ns = (fl_mat->rs) ? 2 / (fl_mat->rs * fl_mat->rs) - 2 : 1e6;
        mat->op = fl_mat->opacity;
        mat->ke_txt = (fl_mat->ke_txt) ? fl_mat->ke_txt->path : "";
        mat->kd_txt = (fl_mat->kd_txt) ? fl_mat->kd_txt->path : "";
        mat->ks_txt = (fl_mat->ks_txt) ? fl_mat->ks_txt->path : "";
        mat->kt_txt = (fl_mat->kt_txt) ? fl_mat->kt_txt->path : "";
        mat->ns_txt = (fl_mat->rs_txt) ? fl_mat->rs_txt->path : "";
        mat->bump_txt = (fl_mat->bump_txt) ? fl_mat->bump_txt->path : "";
        mat->disp_txt = (fl_mat->disp_txt) ? fl_mat->disp_txt->path : "";
        mat->norm_txt = (fl_mat->norm_txt) ? fl_mat->norm_txt->path : "";
        mat->ke_txt_info = texture_props(fl_mat->ke_txt_info);
        mat->kd_txt_info = texture_props(fl_mat->kd_txt_info);
        mat->ks_txt_info = texture_props(fl_mat->ks_txt_info);
        mat->kt_txt_info = texture_props(fl_mat->kt_txt_info);
        mat->ns_txt_info = texture_props(fl_mat->rs_txt_info);
        mat->bump_txt_info = texture_props(fl_mat->bump_txt_info);
        mat->disp_txt_info = texture_props(fl_mat->disp_txt_info);
        mat->norm_txt_info = texture_props(fl_mat->norm_txt_info);
        mat->unknown_props = fl_mat->unknown_props;
        if (fl_mat->opacity < 1 || fl_mat->kt != ym::zero3f) {
            mat->illum = 4;
        } else {
            mat->illum = 2;
        }
    }

    // convert shapes
    for (auto fl_mesh : scn->meshes) {
        asset->objects.emplace_back();
        auto object = &asset->objects.back();
        object->name = fl_mesh->name;
        for (auto fl_prim : fl_mesh->shapes) {
            auto offset = obj_vertex{(int)asset->pos.size(),
                (int)asset->texcoord.size(), (int)asset->norm.size(),
                (int)asset->color.size(), (int)asset->radius.size()};
            for (auto& v : fl_prim->pos) asset->pos.push_back(v);
            for (auto& v : fl_prim->norm) asset->norm.push_back(v);
            for (auto& v : fl_prim->texcoord) asset->texcoord.push_back(v);
            for (auto& v : fl_prim->color) asset->color.push_back(v);
            for (auto& v : fl_prim->radius) asset->radius.push_back(v);
            object->groups.emplace_back();
            auto group = &object->groups.back();
            group->groupname = fl_prim->name;
            group->matname = (fl_prim->mat) ? fl_prim->mat->name : "";
            for (auto point : fl_prim->points) {
                group->elems.push_back({(uint32_t)group->verts.size(),
                    obj_element_type::point, 1});
                group->verts.push_back({(fl_prim->pos.empty()) ?
                                            -1 :
                                            offset.pos + point,
                    (fl_prim->texcoord.empty()) ? -1 : offset.texcoord + point,
                    (fl_prim->norm.empty()) ? -1 : offset.norm + point,
                    (fl_prim->color.empty()) ? -1 : offset.color + point,
                    (fl_prim->radius.empty()) ? -1 : offset.radius + point});
            }
            for (auto line : fl_prim->lines) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::line, 2});
                for (auto vid : line) {
                    group->verts.push_back(
                        {(fl_prim->pos.empty()) ? -1 : offset.pos + vid,
                            (fl_prim->texcoord.empty()) ? -1 :
                                                          offset.texcoord + vid,
                            (fl_prim->norm.empty()) ? -1 : offset.norm + vid,
                            (fl_prim->color.empty()) ? -1 : offset.color + vid,
                            (fl_prim->radius.empty()) ? -1 :
                                                        offset.radius + vid});
                }
            }
            for (auto triangle : fl_prim->triangles) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::face, 3});
                for (auto vid : triangle) {
                    group->verts.push_back(
                        {(fl_prim->pos.empty()) ? -1 : offset.pos + vid,
                            (fl_prim->texcoord.empty()) ? -1 :
                                                          offset.texcoord + vid,
                            (fl_prim->norm.empty()) ? -1 : offset.norm + vid,
                            (fl_prim->color.empty()) ? -1 : offset.color + vid,
                            (fl_prim->radius.empty()) ? -1 :
                                                        offset.radius + vid});
                }
            }
            for (auto tet : fl_prim->tetras) {
                group->elems.push_back({(uint32_t)group->verts.size(),
                    obj_element_type::tetra, 4});
                for (auto vid : tet) {
                    group->verts.push_back(
                        {(fl_prim->pos.empty()) ? -1 : offset.pos + vid,
                            (fl_prim->texcoord.empty()) ? -1 :
                                                          offset.texcoord + vid,
                            (fl_prim->norm.empty()) ? -1 : offset.norm + vid,
                            (fl_prim->color.empty()) ? -1 : offset.color + vid,
                            (fl_prim->radius.empty()) ? -1 :
                                                        offset.radius + vid});
                }
            }
        }
    }

    // convert cameras
    for (auto fl_cam : scn->cameras) {
        asset->cameras.emplace_back();
        auto cam = &asset->cameras.back();
        cam->name = fl_cam->name;
        cam->ortho = fl_cam->ortho;
        cam->yfov = fl_cam->yfov;
        cam->aspect = fl_cam->aspect;
        cam->focus = fl_cam->focus;
        cam->aperture = fl_cam->aperture;
        cam->translation = fl_cam->translation;
        cam->rotation = fl_cam->rotation;
        cam->matrix = fl_cam->matrix;
    }

    // convert envs
    for (auto fl_env : scn->environments) {
        asset->environments.emplace_back();
        auto env = &asset->environments.back();
        env->name = fl_env->name;
        env->matname = (fl_env->mat) ? fl_env->mat->name : "";
        env->rotation = fl_env->rotation;
        env->matrix = fl_env->matrix;
    }

    // convert instances
    for (auto fl_ist : scn->instances) {
        asset->instances.emplace_back();
        auto ist = &asset->instances.back();
        ist->name = fl_ist->name;
        ist->meshname = (fl_ist->msh) ? fl_ist->msh->name : "<undefined>";
        ist->translation = fl_ist->translation;
        ist->rotation = fl_ist->rotation;
        ist->scale = fl_ist->scale;
        ist->matrix = fl_ist->matrix;
    }

    return asset;
}

//
// Loads textures for an scene.
//
bool load_textures(scene* scn, const std::string& dirname, bool skip_missing,
    std::string* err) {
#ifndef YOBJ_NO_IMAGE
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        if (yimg::is_hdr_filename(filename)) {
            txt->hdr = yimg::load_image4f(filename);
        } else {
            txt->ldr = yimg::load_image4b(filename);
        }
        if (!txt->hdr && !txt->ldr) {
            if (skip_missing) continue;
            if (err) *err = "cannot laod image " + filename;
            return false;
        }
    }
#endif
    return true;
}

//
// Loads textures for an scene.
//
bool save_textures(const scene* scn, const std::string& dirname,
    bool skip_missing, std::string* err) {
#ifndef YOBJ_NO_IMAGE
    for (auto txt : scn->textures) {
        if (!txt->ldr && !txt->hdr) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = true;
        if (txt->ldr) { ok = yimg::save_image4b(filename, txt->ldr); }
        if (txt->hdr) { ok = yimg::save_image4f(filename, txt->hdr); }
        if (!ok) {
            if (skip_missing) continue;
            if (err) *err = "cannot save image " + filename;
            return false;
        }
    }
#endif
    return true;
}

//
// Load scene
//
scene* load_scene(const std::string& filename, bool load_txt, bool skip_missing,
    bool flip_texcoord, bool facet_non_smooth, bool flip_tr, std::string* err) {
    auto oscn =
        std::unique_ptr<obj>(load_obj(filename, flip_texcoord, flip_tr, err));
    if (!oscn) return nullptr;
    auto scn = obj_to_scene(oscn.get(), facet_non_smooth);
    if (!scn) return nullptr;
    if (load_txt)
        if (!load_textures(scn, get_dirname(filename), skip_missing, err))
            return nullptr;
    return scn;
}

//
// Save scene
//
bool save_scene(const std::string& filename, const scene* scn, bool save_txt,
    bool flip_texcoord, bool flip_tr, std::string* err) {
    auto oscn = std::unique_ptr<yobj::obj>(scene_to_obj(scn));
    if (!save_obj(filename, oscn.get(), flip_texcoord, flip_tr, err))
        return false;
    if (save_txt)
        if (!save_textures(scn, get_dirname(filename), true, err)) return false;
    return true;
}

//
// Computes a scene bounding box
//
ym::bbox3f compute_scene_bounds(const scene* scn) {
    auto bbox_meshes = std::map<mesh*, ym::bbox3f>();
    for (auto mesh : scn->meshes) {
        bbox_meshes[mesh] = ym::invalid_bbox3f;
        auto& bbox = bbox_meshes[mesh];
        for (auto shp : mesh->shapes) {
            for (auto p : shp->pos) bbox += ym::vec3f(p);
        }
    }
    auto bbox = ym::invalid_bbox3f;
    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            bbox += ym::transform_bbox(
                ym::mat4f(ist->xform()), bbox_meshes[ist->msh]);
        }
    } else {
        for (auto mesh : scn->meshes) { bbox += bbox_meshes[mesh]; }
    }
    return bbox;
}

//
// Add missing data to the scene.
//
void add_normals(scene* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (!shp->norm.empty()) continue;
            shp->norm.resize(shp->pos.size());
            if (!shp->points.empty()) {
                shp->norm.assign(shp->pos.size(), {0, 0, 1});
            } else if (!shp->lines.empty()) {
                ym::compute_tangents(shp->lines, shp->pos, shp->norm);
            } else if (!shp->triangles.empty()) {
                ym::compute_normals(shp->triangles, shp->pos, shp->norm);
            }
        }
    }
}

//
// Add missing data to the scene.
//
void add_tangent_space(scene* scn) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (!shp->mat) continue;
            if (shp->triangles.empty()) continue;
            if (!shp->tangsp.empty() || shp->texcoord.empty() ||
                !shp->mat->norm_txt)
                continue;
            shp->tangsp.resize(shp->pos.size());
            ym::compute_tangent_frame(shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->tangsp);
        }
    }
}

//
// Add missing data to the scene.
//
void add_radius(scene* scn, float radius) {
    for (auto msh : scn->meshes) {
        for (auto shp : msh->shapes) {
            if (shp->points.empty() && shp->lines.empty()) continue;
            if (!shp->radius.empty()) continue;
            shp->radius.resize(shp->pos.size(), radius);
        }
    }
}

//
// Add missing data to the scene.
//
void add_texture_data(scene* scn) {
    for (auto txt : scn->textures) {
        if (!txt->hdr && !txt->ldr) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->ldr = ym::image4b(1, 1, {255, 255, 255, 255});
        }
    }
}

//
// Add missing data to the scene.
//
void add_instances(scene* scn) {
    if (!scn->instances.empty()) return;
    for (auto mesh : scn->meshes) {
        auto ist = new instance();
        ist->name = mesh->name;
        ist->msh = mesh;
        scn->instances.push_back(ist);
    }
}

//
// Add missing data to the scene.
//
void add_names(scene* scn) {
    auto cid = 0;
    for (auto cam : scn->cameras) {
        if (cam->name.empty())
            cam->name = "<camera " + std::to_string(cid) + ">";
        cid++;
    }

    auto mid = 0;
    for (auto mat : scn->materials) {
        if (mat->name.empty())
            mat->name = "<material " + std::to_string(mid) + ">";
        mid++;
    }

    auto mmid = 0;
    for (auto mesh : scn->meshes) {
        if (mesh->name.empty())
            mesh->name = "<mesh " + std::to_string(mmid) + ">";
        mmid++;
        auto sid = 0;
        for (auto shp : mesh->shapes) {
            if (shp->name.empty())
                shp->name = "<shape " + std::to_string(sid) + ">";
            sid++;
        }
    }

    auto iid = 0;
    for (auto ist : scn->instances) {
        if (ist->name.empty())
            ist->name = "<instance " + std::to_string(iid) + ">";
        iid++;
    }

    auto eid = 0;
    for (auto env : scn->environments) {
        if (env->name.empty())
            env->name = "<environment " + std::to_string(eid) + ">";
        eid++;
    }
}

//
// Add a default camera that views the entire scene.
//
void add_default_camera(scene* scn) {
    auto bbox = ym::bbox3f{compute_scene_bounds(scn)};
    if (scn->cameras.empty()) {
        auto center = ym::center(bbox);
        auto bbox_size = ym::diagonal(bbox);
        auto bbox_msize =
            ym::max(bbox_size[0], ym::max(bbox_size[1], bbox_size[2]));
        // set up camera
        auto cam = new camera();
        auto camera_dir = ym::vec3f{1, 0.4f, 1};
        auto from = camera_dir * bbox_msize + center;
        auto to = center;
        auto up = ym::vec3f{0, 1, 0};
        cam->matrix = ym::to_mat(ym::lookat_frame3(from, to, up));
        cam->ortho = false;
        cam->aspect = 16.0f / 9.0f;
        cam->yfov = 2 * atanf(0.5f);
        cam->aperture = 0;
        cam->focus = ym::length(to - from);
        scn->cameras.push_back(cam);
    }
}

//
// Flatten scene instances into separate meshes.
//
void flatten_instances(scene* scn) {
    if (scn->instances.empty()) return;
    auto meshes = scn->meshes;
    scn->meshes.clear();
    auto instances = scn->instances;
    scn->instances.clear();
    for (auto ist : instances) {
        if (!ist->msh) continue;
        auto xf = ist->xform();
        auto msh = ist->msh;
        auto nmsh = new mesh();
        nmsh->name = ist->name;
        for (auto shp : msh->shapes) {
            auto nshp = new shape(*shp);
            for (auto& p : nshp->pos) p = transform_point(xf, p);
            for (auto& n : nshp->norm) n = transform_direction(xf, n);
            nmsh->shapes.push_back(nshp);
        }
        scn->meshes.push_back(nmsh);
    }
    for (auto e : meshes) delete e;
    for (auto e : instances) delete e;
}

//
// Split meshes into single shapes
//
void split_shapes(scene* scn) {
    auto nmeshes = std::vector<mesh*>();
    for (auto msh : scn->meshes) {
        if (msh->shapes.size() <= 1) {
            nmeshes.push_back(msh);
            continue;
        }
        for (auto shp : msh->shapes) {
            auto nmsh = new mesh();
            nmsh->name = msh->name + shp->name;
            nmsh->shapes.push_back(shp);
            nmeshes.push_back(nmsh);
        }
        msh->shapes.clear();
        delete msh;
    }
    scn->meshes = nmeshes;
}

}  // namespace yobj
