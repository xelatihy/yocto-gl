//
// Implementation for Yocto/OBJ.
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_obj.h"

// -----------------------------------------------------------------------------
// OBJ CONVERSION
// -----------------------------------------------------------------------------
namespace yocto {

static inline void parse_value(string_view& str, obj_vertex& value) {
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
static inline void parse_value(string_view& str, obj_texture_info& info) {
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

    // texture params
    auto last = string();
    for (auto i = 0; i < tokens.size() - 1; i++) {
        if (tokens[i] == "-bm") info.scale = atof(tokens[i + 1].c_str());
        if (tokens[i] == "-clamp") info.clamp = true;
    }
}

// Load obj materials
void load_mtl(
    const string& filename, obj_callbacks& cb, const obj_params& params) {
    // open file
    auto fs_ = open_input_file(filename);
    auto fs  = fs_.fs;

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
            if (params.flip_tr) material.kt = vec3f{1, 1, 1} - material.kt;
        } else if (cmd == "Tr") {
            parse_value(line, material.op);
            if (params.flip_tr) material.op = 1 - material.op;
        } else if (cmd == "Ns") {
            parse_value(line, material.ns);
            material.pr = pow(2 / (material.ns + 2), 1 / 4.0f);
            if (material.pr < 0.01f) material.pr = 0;
            if (material.pr > 0.99f) material.pr = 1;
        } else if (cmd == "d") {
            parse_value(line, material.op);
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
        } else if (cmd == "map_bump" || cmd == "bump") {
            parse_value(line, material.bump_txt);
        } else if (cmd == "Pm") {
            material.has_pbr = true;
            parse_value(line, material.pm);
        } else if (cmd == "Pr") {
            material.has_pbr = true;
            parse_value(line, material.pr);
        } else if (cmd == "Ps") {
            material.has_pbr = true;
            parse_value(line, material.ps);
        } else if (cmd == "Pc") {
            material.has_pbr = true;
            parse_value(line, material.pc);
        } else if (cmd == "Pcr") {
            material.has_pbr = true;
            parse_value(line, material.pcr);
        } else if (cmd == "map_Pm") {
            material.has_pbr = true;
            parse_value(line, material.pm_txt);
        } else if (cmd == "map_Pr") {
            material.has_pbr = true;
            parse_value(line, material.pr_txt);
        } else if (cmd == "map_Ps") {
            material.has_pbr = true;
            parse_value(line, material.ps_txt);
        } else if (cmd == "map_occ" || cmd == "occ") {
            parse_value(line, material.occ_txt);
        } else if (cmd == "map_disp" || cmd == "disp") {
            parse_value(line, material.disp_txt);
        } else if (cmd == "map_norm" || cmd == "norm") {
            parse_value(line, material.norm_txt);
        }
    }

    // issue current material
    if (!first) cb.material(material);
}

// Load obj extensions
void load_objx(
    const string& filename, obj_callbacks& cb, const obj_params& params) {
    // open file
    auto fs_ = open_input_file(filename);
    auto fs  = fs_.fs;

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
void load_obj(
    const string& filename, obj_callbacks& cb, const obj_params& params) {
    // open file
    auto fs_ = open_input_file(filename);
    auto fs  = fs_.fs;

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
            if (params.flip_texcoord) vert.y = 1 - vert.y;
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
            if (params.geometry_only) continue;
            auto mtlname = ""s;
            parse_value(line, mtlname);
            cb.mtllib(mtlname);
            auto mtlpath = get_dirname(filename) + mtlname;
            load_mtl(mtlpath, cb, params);
        } else {
            // unused
        }
    }

    // parse extensions if presents
    if (!params.geometry_only) {
        auto extname    = get_noextension(filename) + ".objx";
        auto ext_exists = exists_file(extname);
        if (ext_exists) {
            load_objx(extname, cb, params);
        }
    }
}

}  // namespace yocto
