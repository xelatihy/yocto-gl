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
// Get directory name (including '/').
//
static std::string _get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Splits a std::string into an array of strings on whitespace with Python split
// semantic. Modifies original std::string to avoid allocation.
//
static int _splitws(char* str, char** splits, int maxsplits) {
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
static int _parse_int(char** tok) { return atoi(tok[0]); }

//
// Parses one float.
//
static float _parse_float(char** tok) { return atof(tok[0]); }

//
// Parses two floats.
//
static ym::vec2f _parse_float2(char** tok) {
    return ym::vec2f{(float)atof(tok[0]), (float)atof(tok[1])};
}

//
// Parses three floats.
//
static ym::vec3f _parse_float3(char** tok) {
    return ym::vec3f{
        (float)atof(tok[0]), (float)atof(tok[1]), (float)atof(tok[2])};
}

//
// Parses four floats.
//
static ym::vec4f _parse_float4(char** tok) {
    return ym::vec4f{(float)atof(tok[0]), (float)atof(tok[1]),
        (float)atof(tok[2]), (float)atof(tok[3])};
}

//
// Parses 16 floats.
//
static ym::mat4f _parse_float16(char** tok) {
    ym::mat4f m;
    auto mm = (float*)&m;
    for (auto i = 0; i < 16; i++) mm[i] = (float)atof(tok[i]);
    return m;
}

//
// Parses an OBJ vertex list. Handles negative values.
//
static void _parse_vertlist(char** tok, int ntoks,
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
obj* load_obj(const std::string& filename, bool flip_texcoord) {
    // clear obj
    auto asset = std::unique_ptr<obj>(new obj());

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) throw obj_exception("cannot open filename " + filename);

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
        int ntok = _splitws(line, toks, 1024);

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
            asset->pos.push_back(_parse_float3(cur_tok));
        } else if (tok_s == "vn") {
            vert_size.norm += 1;
            asset->norm.push_back(_parse_float3(cur_tok));
        } else if (tok_s == "vt") {
            vert_size.texcoord += 1;
            asset->texcoord.push_back(_parse_float2(cur_tok));
            if (flip_texcoord)
                asset->texcoord.back()[1] = 1 - asset->texcoord.back()[1];
        } else if (tok_s == "vc") {
            vert_size.color += 1;
            asset->color.push_back(_parse_float4(cur_tok));
        } else if (tok_s == "vr") {
            vert_size.radius += 1;
            asset->radius.push_back(_parse_float(cur_tok));
        } else if (tok_s == "f") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::face,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "l") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::line,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "p") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::point, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "t") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::tetra, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "o") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.push_back({name, {}});
            asset->objects.back().groups.push_back({cur_matname, "", {}, {}});
        } else if (tok_s == "usemtl") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            cur_matname = name;
            asset->objects.back().groups.push_back({cur_matname, "", {}, {}});
        } else if (tok_s == "g") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.back().groups.push_back({cur_matname, name, {}, {}});
        } else if (tok_s == "mtllib") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            cur_mtllibs.push_back(name);
        } else if (tok_s == "c") {
            asset->cameras.emplace_back();
            auto& cam = asset->cameras.back();
            cam.name = (cur_ntok) ? cur_tok[0] : "";
            cam.ortho = _parse_int(cur_tok + 1);
            cam.yfov = _parse_float(cur_tok + 2);
            cam.aspect = _parse_float(cur_tok + 3);
            cam.aperture = _parse_float(cur_tok + 4);
            cam.focus = _parse_float(cur_tok + 5);
            cam.xform = _parse_float16(cur_tok + 6);
        } else if (tok_s == "e") {
            asset->environments.emplace_back();
            auto& env = asset->environments.back();
            env.name = (cur_ntok) ? cur_tok[0] : "<unnamed>";
            env.matname = (cur_ntok - 1) ? cur_tok[1] : "<unnamed_material>";
            env.xform = _parse_float16(cur_tok + 2);
        } else if (tok_s == "i") {
            asset->instances.emplace_back();
            auto& ist = asset->instances.back();
            ist.name = (cur_ntok) ? cur_tok[0] : "<unnamed>";
            ist.meshname = (cur_ntok - 1) ? cur_tok[1] : "<unnamed_mesh>";
            ist.xform = _parse_float16(cur_tok + 2);
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
        auto mtlname = _get_dirname(filename) + mtllib;
        auto materials = load_mtl(mtlname);
        asset->materials.insert(
            asset->materials.end(), materials.begin(), materials.end());
    }

    // done
    return asset.release();
}

//
// Parse texture options and name
//
static void parse_texture(char** toks, int ntoks, std::string& path,
    property_map<std::string>& info) {
    // texture name
    if (ntoks > 0) path = toks[ntoks - 1];

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
std::vector<obj_material> load_mtl(const std::string& filename) {
    // clear materials
    auto materials = std::vector<obj_material>();

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) throw(obj_exception("cannot open filename " + filename));

    // add a material preemptively to avoid crashes
    materials.emplace_back();

    // read the file line by line
    char line[4096];
    char* toks[1024];
    auto linenum = 0;
    while (fgets(line, 4096, file)) {
        linenum += 1;
        int ntok = _splitws(line, toks, 1024);

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
            materials.back().illum = _parse_int(cur_tok);
        } else if (tok_s == "Ke") {
            materials.back().ke = _parse_float3(cur_tok);
        } else if (tok_s == "Ka") {
            materials.back().ka = _parse_float3(cur_tok);
        } else if (tok_s == "Kd") {
            materials.back().kd = _parse_float3(cur_tok);
        } else if (tok_s == "Ks") {
            materials.back().ks = _parse_float3(cur_tok);
        } else if (tok_s == "Kr") {
            materials.back().kr = _parse_float3(cur_tok);
        } else if (tok_s == "Kt" || tok_s == "Tf") {
            materials.back().kt = _parse_float3(cur_tok);
        } else if (tok_s == "Tr") {
            if (cur_ntok >= 3) {
                materials.back().kt = _parse_float3(cur_tok);
            } else {
                // as tinyobjreader
                materials.back().op = 1 - _parse_float(cur_tok);
            }
        } else if (tok_s == "Ns") {
            materials.back().ns = _parse_float(cur_tok);
        } else if (tok_s == "d") {
            materials.back().op = _parse_float(cur_tok);
        } else if (tok_s == "Ni") {
            materials.back().ior = _parse_float(cur_tok);
        } else if (tok_s == "phys_stiffness") {
            materials.back().stiffness = _parse_float(cur_tok);
        } else if (tok_s == "phys_density") {
            materials.back().density = _parse_float(cur_tok);
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
static void _fwrite_float(
    FILE* file, const char* str, float v, bool newline = true) {
    fprintf(file, "%s %.6g", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write one float prepended by a std::string
//
static void _fwrite_int(
    FILE* file, const char* str, int v, bool newline = true) {
    fprintf(file, "%s %d", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write two floats prepended by a std::string
//
static void _fwrite_float2(
    FILE* file, const char* str, const ym::vec2f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g", str, v[0], v[1]);
    if (newline) fprintf(file, "\n");
}

//
// write three floats prepended by a std::string
//
static void _fwrite_float3(
    FILE* file, const char* str, const ym::vec3f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g %.6g", str, v[0], v[1], v[2]);
    if (newline) fprintf(file, "\n");
}

//
// write four floats prepended by a std::string
//
static void _fwrite_float4(
    FILE* file, const char* str, const ym::vec4f& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g %.6g %.6g", str, v[0], v[1], v[2], v[3]);
    if (newline) fprintf(file, "\n");
}

//
// write 16 floats prepended by a std::string
//
static void _fwrite_float16(
    FILE* file, const char* str, const ym::mat4f& v, bool newline = true) {
    const float* vf = (float*)&v;
    fprintf(file, "%s", str);
    for (int i = 0; i < 16; i++) fprintf(file, " %.6g", vf[i]);
    if (newline) fprintf(file, "\n");
}

//
// write a std::string prepended by another if the std::string is not NULL
//
static void _fwrite_str(FILE* file, const char* str, const std::string& s,
    bool force = false, bool newline = true) {
    if (s.empty() && !force) return;
    fprintf(file, "%s %s", str, s.c_str());
    if (newline) fprintf(file, "\n");
}

//
// write a std::string prepended by another if the std::string is not NULL
//
static void _fwrite_str_props(FILE* file, const char* str, const std::string& s,
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
// write an OBJ vertex triplet using only the indices that are active
//
static void _fwrite_objverts(FILE* file, const char* str, int nv,
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
void save_obj(
    const std::string& filename, const obj* asset, bool flip_texcoord) {
    // open file
    auto file = fopen(filename.c_str(), "wt");
    if (!file) throw(obj_exception("could not open filename " + filename));

    // linkup to mtl
    auto dirname = _get_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!asset->materials.empty()) {
        _fwrite_str(file, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto& cam : asset->cameras) {
        _fwrite_str(file, "c", cam.name, true, false);
        _fwrite_int(file, " ", cam.ortho, false);
        _fwrite_float(file, " ", cam.yfov, false);
        _fwrite_float(file, " ", cam.aspect, false);
        _fwrite_float(file, " ", cam.aperture, false);
        _fwrite_float(file, " ", cam.focus, false);
        _fwrite_float16(file, " ", cam.xform, true);
    }

    // save envs
    for (auto& env : asset->environments) {
        _fwrite_str(file, "e", env.name, true, false);
        _fwrite_str(file, " ", env.matname, true, false);
        _fwrite_float16(file, " ", env.xform, true);
    }

    // save instances
    for (auto& ist : asset->instances) {
        _fwrite_str(file, "i", ist.name, true, false);
        _fwrite_str(file, " ", ist.meshname, true, false);
        _fwrite_float16(file, " ", ist.xform, true);
    }

    // save all vertex data
    for (auto& v : asset->pos) _fwrite_float3(file, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            _fwrite_float2(file, "vt", {v[0], 1 - v[1]});
    } else {
        for (auto& v : asset->texcoord) _fwrite_float2(file, "vt", v);
    }
    for (auto& v : asset->norm) _fwrite_float3(file, "vn", v);
    for (auto& v : asset->color) _fwrite_float4(file, "vc", v);
    for (auto& v : asset->radius) _fwrite_float(file, "vr", v);

    // save element data
    const char* elem_labels[] = {"", "p", "l", "f", "t"};
    for (auto& object : asset->objects) {
        _fwrite_str(file, "o", object.name, true);
        for (auto& group : object.groups) {
            _fwrite_str(file, "usemtl", group.matname);
            _fwrite_str(file, "g", group.groupname);
            for (auto elem : group.elems) {
                _fwrite_objverts(file, elem_labels[(int)elem.type], elem.size,
                    group.verts.data() + elem.start);
            }
        }
    }

    fclose(file);

    // save materials
    if (!asset->materials.empty()) {
        save_mtl(dirname + basename + ".mtl", asset->materials);
    }
}

//
// Save an MTL file
//
void save_mtl(
    const std::string& filename, const std::vector<obj_material>& materials) {
    auto file = fopen(filename.c_str(), "wt");
    if (!file) throw(obj_exception("could not open filename " + filename));

    // for each material, dump all the values
    for (auto& mat : materials) {
        _fwrite_str(file, "newmtl", mat.name, true);
        _fwrite_int(file, "  illum", mat.illum);
        _fwrite_float3(file, "  Ke", mat.ke);
        _fwrite_float3(file, "  Ka", mat.ka);
        _fwrite_float3(file, "  Kd", mat.kd);
        _fwrite_float3(file, "  Ks", mat.ks);
        _fwrite_float3(file, "  Kr", mat.kr);
        _fwrite_float3(file, "  Kt", mat.kt);
        _fwrite_float(file, "  Ns", mat.ns);
        _fwrite_float(file, "  d", mat.op);
        _fwrite_float(file, "  Ni", mat.ior);
        _fwrite_float(file, "  phys_stiffness", mat.stiffness);
        _fwrite_float(file, "  phys_density", mat.density);
        _fwrite_str_props(file, "  map_Ke", mat.ke_txt, mat.ke_txt_info);
        _fwrite_str_props(file, "  map_Ka", mat.ka_txt, mat.ka_txt_info);
        _fwrite_str_props(file, "  map_Kd", mat.kd_txt, mat.kd_txt_info);
        _fwrite_str_props(file, "  map_Ks", mat.ks_txt, mat.ks_txt_info);
        _fwrite_str_props(file, "  map_Kr", mat.kr_txt, mat.kr_txt_info);
        _fwrite_str_props(file, "  map_Kt", mat.kt_txt, mat.kt_txt_info);
        _fwrite_str_props(file, "  map_Ns", mat.ns_txt, mat.ns_txt_info);
        _fwrite_str_props(file, "  map_d", mat.op_txt, mat.op_txt_info);
        _fwrite_str_props(file, "  map_Ni", mat.ior_txt, mat.ior_txt_info);
        _fwrite_str_props(file, "  map_bump", mat.bump_txt, mat.bump_txt_info);
        _fwrite_str_props(file, "  map_disp", mat.disp_txt, mat.disp_txt_info);
        _fwrite_str_props(file, "  map_norm", mat.norm_txt, mat.norm_txt_info);
        for (auto&& p : mat.unknown_props) {
            auto s = std::string();
            for (auto&& v : p.second) s += v + " ";
            _fwrite_str(file, p.first.c_str(), s.c_str());
        }
        fprintf(file, "\n");
    }

    fclose(file);
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
static texture* _add_texture(
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
scene* obj_to_scene(const obj* asset) {
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
        mat->stiffness = omat.stiffness;
        mat->density = omat.density;
        mat->ke_txt = _add_texture(omat.ke_txt, scn->textures);
        mat->kd_txt = _add_texture(omat.kd_txt, scn->textures);
        mat->ks_txt = _add_texture(omat.ks_txt, scn->textures);
        mat->kt_txt = _add_texture(omat.kt_txt, scn->textures);
        mat->rs_txt = _add_texture(omat.ns_txt, scn->textures);
        mat->norm_txt = _add_texture(omat.norm_txt, scn->textures);
        mat->bump_txt = _add_texture(omat.bump_txt, scn->textures);
        mat->disp_txt = _add_texture(omat.disp_txt, scn->textures);
        mat->ke_txt_info = _convert_texture_info(omat.ke_txt_info);
        mat->kd_txt_info = _convert_texture_info(omat.kd_txt_info);
        mat->ks_txt_info = _convert_texture_info(omat.ks_txt_info);
        mat->kt_txt_info = _convert_texture_info(omat.kt_txt_info);
        mat->rs_txt_info = _convert_texture_info(omat.ns_txt_info);
        mat->norm_txt_info = _convert_texture_info(omat.norm_txt_info);
        mat->bump_txt_info = _convert_texture_info(omat.bump_txt_info);
        mat->disp_txt_info = _convert_texture_info(omat.disp_txt_info);
        mat->unknown_props = omat.unknown_props;
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

            // covert elements
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
        cam->xform = ocam.xform;
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
        env->xform = oenv.xform;
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
        ist->xform = oist.xform;
        scn->instances.push_back(ist);
    }

    // done
    return scn;
}

//
// Convert texture props
//
property_map<std::string> _texture_props(const texture_info& info) {
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
        mat->stiffness = fl_mat->stiffness;
        mat->density = fl_mat->density;
        mat->ke_txt = (fl_mat->ke_txt) ? fl_mat->ke_txt->path : "";
        mat->kd_txt = (fl_mat->kd_txt) ? fl_mat->kd_txt->path : "";
        mat->ks_txt = (fl_mat->ks_txt) ? fl_mat->ks_txt->path : "";
        mat->kt_txt = (fl_mat->kt_txt) ? fl_mat->kt_txt->path : "";
        mat->ns_txt = (fl_mat->rs_txt) ? fl_mat->rs_txt->path : "";
        mat->bump_txt = (fl_mat->bump_txt) ? fl_mat->bump_txt->path : "";
        mat->disp_txt = (fl_mat->disp_txt) ? fl_mat->disp_txt->path : "";
        mat->norm_txt = (fl_mat->norm_txt) ? fl_mat->norm_txt->path : "";
        mat->ke_txt_info = _texture_props(fl_mat->ke_txt_info);
        mat->kd_txt_info = _texture_props(fl_mat->kd_txt_info);
        mat->ks_txt_info = _texture_props(fl_mat->ks_txt_info);
        mat->kt_txt_info = _texture_props(fl_mat->kt_txt_info);
        mat->ns_txt_info = _texture_props(fl_mat->rs_txt_info);
        mat->bump_txt_info = _texture_props(fl_mat->bump_txt_info);
        mat->disp_txt_info = _texture_props(fl_mat->disp_txt_info);
        mat->norm_txt_info = _texture_props(fl_mat->norm_txt_info);
        mat->unknown_props = fl_mat->unknown_props;
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
        cam->xform = fl_cam->xform;
    }

    // convert envs
    for (auto fl_env : scn->environments) {
        asset->environments.emplace_back();
        auto env = &asset->environments.back();
        env->name = fl_env->name;
        env->matname = (fl_env->mat) ? fl_env->mat->name : "";
        env->xform = fl_env->xform;
    }

    // convert instances
    for (auto fl_ist : scn->instances) {
        asset->instances.emplace_back();
        auto ist = &asset->instances.back();
        ist->name = fl_ist->name;
        ist->meshname = (fl_ist->msh) ? fl_ist->msh->name : "<undefined>";
        ist->xform = fl_ist->xform;
    }

    return asset;
}

//
// Loads textures for an scene.
//
void load_textures(scene* scn, const std::string& dirname, bool skip_missing) {
#ifndef YOBJ_NO_IMAGE
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            yimg::load_image(filename, txt->width, txt->height, txt->ncomp,
                txt->dataf, txt->datab);
        } catch (...) {
            if (!skip_missing) throw;
        }
    }
#endif
}

//
// Loads textures for an scene.
//
void save_textures(
    const scene* scn, const std::string& dirname, bool skip_missing) {
#ifndef YOBJ_NO_IMAGE
    for (auto txt : scn->textures) {
        if (txt->datab.empty() && txt->dataf.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        if (!txt->dataf.empty()) {
            yimg::save_image(filename, txt->width, txt->height, txt->ncomp,
                txt->dataf.data());
        }
        if (!txt->datab.empty()) {
            yimg::save_image(filename, txt->width, txt->height, txt->ncomp,
                txt->datab.data());
        }
    }
#endif
}

//
// Load scene
//
scene* load_scene(const std::string& filename, bool load_txt,
    bool flip_texcoord, bool skip_missing) {
    auto oscn = std::unique_ptr<obj>(load_obj(filename, flip_texcoord));
    auto scn = obj_to_scene(oscn.get());
    if (load_txt) load_textures(scn, _get_dirname(filename), skip_missing);
    return scn;
}

//
// Save scene
//
void save_scene(const std::string& filename, const scene* scn, bool save_txt,
    bool flip_texcoord) {
    auto oscn = std::unique_ptr<yobj::obj>(scene_to_obj(scn));
    save_obj(filename, oscn.get(), flip_texcoord);
    if (save_txt) save_textures(scn, _get_dirname(filename), true);
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
                ym::mat4f(ist->xform), bbox_meshes[ist->msh]);
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
        if (txt->dataf.empty() && txt->datab.empty()) {
            printf("unable to load texture %s\n", txt->path.c_str());
            txt->width = 1;
            txt->height = 1;
            txt->ncomp = 4;
            txt->datab = {{255, 255, 255, 255}};
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
        cam->xform = ym::to_mat(ym::lookat_frame3(from, to, up));
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
        auto msh = ist->msh;
        auto nmsh = new mesh();
        nmsh->name = ist->name;
        for (auto shp : msh->shapes) {
            auto nshp = new shape(*shp);
            for (auto& p : nshp->pos) p = transform_point(ist->xform, p);
            for (auto& n : nshp->norm) n = transform_direction(ist->xform, n);
            nmsh->shapes.push_back(nshp);
        }
        scn->meshes.push_back(nmsh);
    }
    for (auto e : meshes) delete e;
    for (auto e : instances) delete e;
}

}  // namespace yobj
