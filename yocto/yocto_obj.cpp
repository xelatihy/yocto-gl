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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR YOCTO_OBJ
// -----------------------------------------------------------------------------

#include "yocto_obj.h"

#include <algorithm>
#include <cassert>
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
static inline std::string _get_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

//
// Get extension (including '.').
//
static inline std::string _get_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

//
// Parse element buffer
//
static inline void _parse_vert_array(char* line_parse,
    std::vector<vert>& element_buffer, const vert& vert_size) {
    element_buffer.clear();
    char vert_buf[2048];
    int advance = 0;
    while (line_parse[0] && sscanf(line_parse, "%s%n", vert_buf, &advance)) {
        line_parse += advance;
        auto v = vert(-1, -1, -1, -1, -1);
        auto vp = &v.pos;
        auto sp = &v.pos;
        auto vert_parse = vert_buf;
        for (auto i = 0; i < 5; i++) {
            if (!vert_parse[0]) break;
            if (vert_parse[0] == '/') {
                vert_parse++;
                continue;
            }
            sscanf(vert_parse, "%d%n", vp + i, &advance);
            vert_parse += advance;
            vp[i] = (vp[i] < 0) ? sp[i] + vp[i] : vp[i] - 1;
            if (vert_parse[0] == '/') vert_parse += 1;
        }
        element_buffer.push_back(v);
    }
}

//
// Splits a std::string into an array of strings on whitespace with Python split
// semantic. Modifies original std::string to avoid allocation.
//
static inline int _splitws(char* str, char** splits, int maxsplits) {
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
static inline int _parse_int(char** tok) { return atoi(tok[0]); }

//
// Parses one float.
//
static inline float _parse_float(char** tok) { return atof(tok[0]); }

//
// Parses two floats.
//
static inline float2 _parse_float2(char** tok) {
    return float2{(float)atof(tok[0]), (float)atof(tok[1])};
}

//
// Parses three floats.
//
static inline float3 _parse_float3(char** tok) {
    return float3{
        (float)atof(tok[0]), (float)atof(tok[1]), (float)atof(tok[2])};
}

//
// Parses 16 floats.
//
static inline float16 _parse_float16(char** tok) {
    float16 m;
    auto mm = (float*)&m;
    for (auto i = 0; i < 16; i++) mm[i] = (float)atof(tok[i]);
    return m;
}

//
// Parses an OBJ vertex list. Handles negative values.
//
static inline void _parse_vertlist(
    char** tok, int ntoks, std::vector<vert>& elems, const vert& vert_size) {
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
        auto v = vert{-1, -1, -1, -1, -1};
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
YOBJ_API obj* load_obj(const std::string& filename, bool flip_texcoord) {
    // clear obj
    auto asset = std::make_unique<obj>();

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) throw obj_exception("cannot open filename " + filename);

    // initializing obj
    asset->objects.emplace_back();
    asset->objects.back().elems.emplace_back();

    // allocate buffers to avoid re-allocing
    auto cur_elems = std::vector<vert>();
    auto cur_matname = std::string();
    auto cur_mtllibs = std::vector<std::string>();

    // keep track of array lengths
    auto vert_size = vert{0, 0, 0, 0, 0};

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
            asset->color.push_back(_parse_float3(cur_tok));
        } else if (tok_s == "vr") {
            vert_size.radius += 1;
            asset->radius.push_back(_parse_float(cur_tok));
        } else if (tok_s == "f") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().elems.back();
            g.elems.push_back({(uint32_t)g.verts.size(), elem::etype::face,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "l") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().elems.back();
            g.elems.push_back({(uint32_t)g.verts.size(), elem::etype::line,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "p") {
            _parse_vertlist(cur_tok, cur_ntok, cur_elems, vert_size);
            auto& g = asset->objects.back().elems.back();
            g.elems.push_back({(uint32_t)g.verts.size(), elem::etype::point,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (tok_s == "o") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.push_back({name, {}});
            asset->objects.back().elems.push_back({cur_matname, "", {}, {}});
        } else if (tok_s == "usemtl") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            cur_matname = name;
            asset->objects.back().elems.push_back({cur_matname, "", {}, {}});
        } else if (tok_s == "g") {
            auto name = (cur_ntok) ? cur_tok[0] : "";
            asset->objects.back().elems.push_back({cur_matname, name, {}, {}});
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
            env.name = (cur_ntok) ? cur_tok[0] : "";
            env.matname = (cur_ntok - 1) ? cur_tok[1] : "";
            env.xform = _parse_float16(cur_tok + 2);
        } else {
            // unused
        }
    }

    // cleanup unused
    for (auto&& o : asset->objects) {
        auto end = std::remove_if(o.elems.begin(), o.elems.end(),
            [](const elem_group& x) { return x.verts.empty(); });
        o.elems.erase(end, o.elems.end());
    }
    auto end = std::remove_if(asset->objects.begin(), asset->objects.end(),
        [](const object& x) { return x.elems.empty(); });
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
// Load MTL
//
YOBJ_API std::vector<material> load_mtl(const std::string& filename) {
    // clear materials
    auto materials = std::vector<material>();

    // open file
    auto file = fopen(filename.c_str(), "rt");
    if (!file) throw(obj_exception("cannot open filename " + filename));

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
        } else if (tok_s == "map_Ke") {
            materials.back().ke_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Ka") {
            materials.back().ka_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Kd") {
            materials.back().kd_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Ks") {
            materials.back().ks_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Kr") {
            materials.back().ke_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Tr") {
            materials.back().kt_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Ns") {
            materials.back().ns_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_d") {
            materials.back().op_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_Ni") {
            materials.back().ior_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_bump" || tok_s == "bump") {
            materials.back().bump_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_disp" || tok_s == "disp") {
            materials.back().disp_txt = (cur_ntok) ? cur_tok[0] : "";
        } else if (tok_s == "map_norm" || tok_s == "norm") {
            materials.back().norm_txt = (cur_ntok) ? cur_tok[0] : "";
        } else {
            // unused
        }
    }

    // done
    return materials;
}

//
// write one float prepended by a std::string
//
static inline void _fwrite_float(
    FILE* file, const char* str, float v, bool newline = true) {
    fprintf(file, "%s %.6g", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write one float prepended by a std::string
//
static inline void _fwrite_int(
    FILE* file, const char* str, int v, bool newline = true) {
    fprintf(file, "%s %d", str, v);
    if (newline) fprintf(file, "\n");
}

//
// write two floats prepended by a std::string
//
static inline void _fwrite_float2(
    FILE* file, const char* str, const float2& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g", str, v[0], v[1]);
    if (newline) fprintf(file, "\n");
}

//
// write three floats prepended by a std::string
//
static inline void _fwrite_float3(
    FILE* file, const char* str, const float3& v, bool newline = true) {
    fprintf(file, "%s %.6g %.6g %.6g", str, v[0], v[1], v[2]);
    if (newline) fprintf(file, "\n");
}

//
// write 16 floats prepended by a std::string
//
static inline void _fwrite_float16(
    FILE* file, const char* str, const float16& v, bool newline = true) {
    const float* vf = (float*)&v;
    fprintf(file, "%s", str);
    for (int i = 0; i < 16; i++) fprintf(file, " %.6g", vf[i]);
    if (newline) fprintf(file, "\n");
}

//
// write a std::string prepended by another if the std::string is not NULL
//
static inline void _fwrite_str(FILE* file, const char* str,
    const std::string& s, bool force = false, bool newline = true) {
    if (s.empty() && !force) return;
    fprintf(file, "%s %s", str, s.c_str());
    if (newline) fprintf(file, "\n");
}

//
// write an OBJ vertex triplet using only the indices that are active
//
static inline void _fwrite_objverts(FILE* file, const char* str, int nv,
    const vert* verts, bool newline = true) {
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
YOBJ_API void save_obj(
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

    // save all vertex data
    for (auto& v : asset->pos) _fwrite_float3(file, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            _fwrite_float2(file, "vt", {v[0], 1 - v[1]});
    } else {
        for (auto& v : asset->texcoord) _fwrite_float2(file, "vt", v);
    }
    for (auto& v : asset->norm) _fwrite_float3(file, "vn", v);
    for (auto& v : asset->color) _fwrite_float3(file, "vc", v);
    for (auto& v : asset->radius) _fwrite_float(file, "vr", v);

    // save element data
    const char* elem_labels[] = {"", "p", "l", "f"};
    for (auto& object : asset->objects) {
        _fwrite_str(file, "o", object.name, true);
        for (auto& elems : object.elems) {
            _fwrite_str(file, "usemtl", elems.matname);
            _fwrite_str(file, "g", elems.groupname);
            for (auto elem : elems.elems) {
                _fwrite_objverts(file, elem_labels[(int)elem.type], elem.size,
                    elems.verts.data() + elem.start);
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
YOBJ_API void save_mtl(
    const std::string& filename, const std::vector<material>& materials) {
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
        _fwrite_str(file, "  map_Ke", mat.ke_txt);
        _fwrite_str(file, "  map_Ka", mat.ke_txt);
        _fwrite_str(file, "  map_Kd", mat.kd_txt);
        _fwrite_str(file, "  map_Ks", mat.ks_txt);
        _fwrite_str(file, "  map_Kr", mat.kr_txt);
        _fwrite_str(file, "  map_Kt", mat.kt_txt);
        _fwrite_str(file, "  map_Ns", mat.ns_txt);
        _fwrite_str(file, "  map_d", mat.op_txt);
        _fwrite_str(file, "  map_Ni", mat.ior_txt);
        _fwrite_str(file, "  map_bump", mat.bump_txt);
        _fwrite_str(file, "  map_disp", mat.disp_txt);
        _fwrite_str(file, "  map_norm", mat.norm_txt);
        fprintf(file, "\n");
    }

    fclose(file);
}

//
// Cleanup memory
//
YOBJ_API fl_obj::~fl_obj() {
    for (auto v : primitives)
        if (v) delete v;
    for (auto v : meshes)
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
static inline int _add_texture(
    const std::string& filename, std::vector<fl_texture*>& txts) {
    if (filename.empty()) return -1;
    for (auto i = 0; i < txts.size(); i++) {
        if (txts[i]->path == filename) return i;
    }
    auto txt = new fl_texture();
    txt->path = filename;
    txts.push_back(txt);
    return (int)txts.size() - 1;
}

//
// A hash function for vecs
//
struct _vert_hash {
    std::hash<int> Th;
    size_t operator()(const vert& vv) const {
        auto v = (const int*)&vv;
        size_t h = 0;
        for (auto i = 0; i < sizeof(vert) / sizeof(int); i++) {
            // embads hash_combine below
            h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
        }
        return h;
    }
};

//
// Comparison for unordred_map
//
static inline bool operator==(const vert& a, const vert& b) {
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm &&
           a.color == b.color && a.radius == b.radius;
}

//
// Flattens an scene
//
YOBJ_API fl_obj* flatten_obj(const obj* asset) {
    // clear scene
    auto scene = new fl_obj();

    // convert materials and build textures
    for (auto& omat : asset->materials) {
        auto mat = new fl_material();
        mat->name = omat.name;
        mat->ke = omat.ke;
        mat->kd = omat.kd;
        mat->ks = omat.ks;
        mat->kt = omat.kt;
        mat->rs = sqrt(2 / (omat.ns + 2));
        mat->ke_txt = _add_texture(omat.ke_txt, scene->textures);
        mat->kd_txt = _add_texture(omat.kd_txt, scene->textures);
        mat->ks_txt = _add_texture(omat.ks_txt, scene->textures);
        mat->kt_txt = _add_texture(omat.kt_txt, scene->textures);
        mat->rs_txt = _add_texture(omat.ns_txt, scene->textures);
        mat->norm_txt = _add_texture(omat.norm_txt, scene->textures);
        scene->materials.push_back(mat);
    }

    // convert shapes
    for (auto& oshape : asset->objects) {
        auto mesh = new fl_mesh();
        mesh->name = oshape.name;
        for (auto& elem_group : oshape.elems) {
            if (elem_group.verts.empty()) continue;
            if (elem_group.elems.empty()) continue;
            auto prim = new fl_primitives();
            prim->name = elem_group.groupname;
            prim->material = -1;
            for (auto i = 0; i < asset->materials.size() && prim->material < 0;
                 i++) {
                if (asset->materials[i].name == elem_group.matname)
                    prim->material = i;
            }

            // insert all vertices
            std::unordered_map<vert, int, _vert_hash> vert_map;
            std::vector<int> vert_ids;
            // vert_map.clear();
            // vert_ids.clear();
            for (auto& vert : elem_group.verts) {
                if (vert_map.find(vert) == vert_map.end()) {
                    vert_map[vert] = (int)vert_map.size();
                }
                vert_ids.push_back(vert_map.at(vert));
            }

            // covert elements
            for (auto& elem : elem_group.elems) {
                switch (elem.type) {
                    case elem::etype::point: {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i++) {
                            prim->points.push_back(vert_ids[i]);
                        }
                    } break;
                    case elem::etype::line: {
                        for (auto i = elem.start;
                             i < elem.start + elem.size - 1; i++) {
                            prim->lines.push_back(
                                {vert_ids[i], vert_ids[i + 1]});
                        }
                    } break;
                    case elem::etype::face: {
                        for (auto i = elem.start + 2;
                             i < elem.start + elem.size; i++) {
                            prim->triangles.push_back({vert_ids[elem.start],
                                vert_ids[i - 1], vert_ids[i]});
                        }
                    } break;
                    default: { assert(false); }
                }
            }

            // check for errors
            // copy vertex data
            auto v = elem_group.verts[0];
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
            scene->primitives.push_back(prim);
            mesh->primitives.push_back((int)scene->primitives.size() - 1);
        }
        scene->meshes.push_back(mesh);
    }

    // convert cameras
    for (auto& ocam : asset->cameras) {
        auto cam = new fl_camera();
        cam->name = ocam.name;
        cam->ortho = ocam.ortho;
        cam->yfov = ocam.yfov;
        cam->aspect = ocam.aspect;
        cam->aperture = ocam.aperture;
        cam->focus = ocam.focus;
        cam->xform = ocam.xform;
        scene->cameras.push_back(cam);
    }

    // convert envs
    for (auto& oenv : asset->environments) {
        auto env = new fl_environment();
        env->name = oenv.name;
        env->matid = -1;
        for (auto i = 0; i < asset->materials.size() && env->matid < 0; i++) {
            if (asset->materials[i].name == oenv.matname) env->matid = i;
        }
        env->xform = oenv.xform;
        scene->environments.push_back(env);
    }

    // done
    return scene;
}

//
// Save an scene
//
YOBJ_API obj* unflatten_obj(const fl_obj* scene) {
    auto asset = new obj();

    // get texture name helper
    auto txt = [](const fl_obj* scene, int id) {
        if (id < 0) return std::string();
        return scene->textures[id]->path;
    };

    // convert materials
    for (auto fl_mat : scene->materials) {
        asset->materials.emplace_back();
        auto mat = &asset->materials.back();
        mat->name = fl_mat->name;
        mat->ke = fl_mat->ke;
        mat->kd = fl_mat->kd;
        mat->ks = fl_mat->ks;
        mat->kt = fl_mat->kt;
        mat->ns = (fl_mat->rs) ? 2 / (fl_mat->rs * fl_mat->rs) - 2 : 1e6;
        mat->ke_txt = txt(scene, fl_mat->ke_txt);
        mat->kd_txt = txt(scene, fl_mat->kd_txt);
        mat->ks_txt = txt(scene, fl_mat->ks_txt);
        mat->kt_txt = txt(scene, fl_mat->kt_txt);
        mat->ns_txt = txt(scene, fl_mat->rs_txt);
        mat->norm_txt = txt(scene, fl_mat->norm_txt);
    }

    // convert shapes
    for (auto fl_mesh : scene->meshes) {
        asset->objects.emplace_back();
        auto object = &asset->objects.back();
        object->name = fl_mesh->name;
        for (auto fl_prim_id : fl_mesh->primitives) {
            auto fl_prim = scene->primitives[fl_prim_id];
            auto offset = vert{(int)asset->pos.size(),
                (int)asset->texcoord.size(), (int)asset->norm.size(),
                (int)asset->color.size(), (int)asset->radius.size()};
            for (auto& v : fl_prim->pos) asset->pos.push_back(v);
            for (auto& v : fl_prim->norm) asset->norm.push_back(v);
            for (auto& v : fl_prim->texcoord) asset->texcoord.push_back(v);
            for (auto& v : fl_prim->color) asset->color.push_back(v);
            for (auto& v : fl_prim->radius) asset->radius.push_back(v);
            object->elems.emplace_back();
            auto elems = &object->elems.back();
            elems->groupname = fl_prim->name;
            elems->matname = (fl_prim->material < 0) ?
                                 "" :
                                 scene->materials[fl_prim->material]->name;
            for (auto point : fl_prim->points) {
                elems->elems.push_back(
                    {(uint32_t)elems->verts.size(), elem::etype::point, 1});
                elems->verts.push_back({(fl_prim->pos.empty()) ?
                                            -1 :
                                            offset.pos + point,
                    (fl_prim->texcoord.empty()) ? -1 : offset.texcoord + point,
                    (fl_prim->norm.empty()) ? -1 : offset.norm + point,
                    (fl_prim->color.empty()) ? -1 : offset.color + point,
                    (fl_prim->radius.empty()) ? -1 : offset.radius + point});
            }
            for (auto line : fl_prim->lines) {
                elems->elems.push_back(
                    {(uint32_t)elems->verts.size(), elem::etype::line, 2});
                elems->verts.push_back(
                    {(fl_prim->pos.empty()) ? -1 : offset.pos + line[0],
                        (fl_prim->texcoord.empty()) ? -1 :
                                                      offset.texcoord + line[0],
                        (fl_prim->norm.empty()) ? -1 : offset.norm + line[0],
                        (fl_prim->color.empty()) ? -1 : offset.color + line[0],
                        (fl_prim->radius.empty()) ? -1 :
                                                    offset.radius + line[0]});
                elems->verts.push_back(
                    {(fl_prim->pos.empty()) ? -1 : offset.pos + line[1],
                        (fl_prim->texcoord.empty()) ? -1 :
                                                      offset.texcoord + line[1],
                        (fl_prim->norm.empty()) ? -1 : offset.norm + line[1],
                        (fl_prim->color.empty()) ? -1 : offset.color + line[1],
                        (fl_prim->radius.empty()) ? -1 :
                                                    offset.radius + line[1]});
            }
            for (auto triangle : fl_prim->triangles) {
                elems->elems.push_back(
                    {(uint32_t)elems->verts.size(), elem::etype::face, 3});
                elems->verts.push_back({(fl_prim->pos.empty()) ?
                                            -1 :
                                            offset.pos + triangle[0],
                    (fl_prim->texcoord.empty()) ? -1 :
                                                  offset.texcoord + triangle[0],
                    (fl_prim->norm.empty()) ? -1 : offset.norm + triangle[0],
                    (fl_prim->color.empty()) ? -1 : offset.color + triangle[0],
                    (fl_prim->radius.empty()) ? -1 :
                                                offset.radius + triangle[0]});
                elems->verts.push_back({(fl_prim->pos.empty()) ?
                                            -1 :
                                            offset.pos + triangle[1],
                    (fl_prim->texcoord.empty()) ? -1 :
                                                  offset.texcoord + triangle[1],
                    (fl_prim->norm.empty()) ? -1 : offset.norm + triangle[1],
                    (fl_prim->color.empty()) ? -1 : offset.color + triangle[1],
                    (fl_prim->radius.empty()) ? -1 :
                                                offset.radius + triangle[1]});
                elems->verts.push_back({(fl_prim->pos.empty()) ?
                                            -1 :
                                            offset.pos + triangle[2],
                    (fl_prim->texcoord.empty()) ? -1 :
                                                  offset.texcoord + triangle[2],
                    (fl_prim->norm.empty()) ? -1 : offset.norm + triangle[2],
                    (fl_prim->color.empty()) ? -1 : offset.color + triangle[2],
                    (fl_prim->radius.empty()) ? -1 :
                                                offset.radius + triangle[2]});
            }
        }
    }

    // convert cameras
    for (auto fl_cam : scene->cameras) {
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
    for (auto fl_env : scene->environments) {
        asset->environments.emplace_back();
        auto env = &asset->environments.back();
        env->name = fl_env->name;
        env->matname =
            (fl_env->matid < 0) ? "" : scene->materials[fl_env->matid]->name;
        env->xform = fl_env->xform;
    }

    return asset;
}

//
// Loads textures for an scene.
//
YOBJ_API void load_textures(
    fl_obj* scene, const std::string& dirname, bool skip_missing) {
#ifndef YOBJ_NO_IMAGE
    for (auto txt : scene->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto img = yimg::load_image(filename);
        txt->width = img->width;
        txt->height = img->height;
        txt->ncomp = img->ncomp;
        if (img->hdr)
            txt->dataf = std::vector<float>(
                img->hdr, img->hdr + img->width * img->height * img->ncomp);
        if (img->ldr)
            txt->datab = std::vector<unsigned char>(
                img->ldr, img->ldr + img->width * img->height * img->ncomp);
        delete img;
    }
#endif
}

}  // namespace
