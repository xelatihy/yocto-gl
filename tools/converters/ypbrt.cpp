//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice,
// this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
// this list of conditions and the following disclaimer in the documentation
// and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//

//
// Converts scenes from https://benedikt-bitterli.me/resources/ to Yocto/GL.
//

#include <regex>
#include "../../yocto/yocto_gl.h"
#include "ext/json.hpp"

using namespace ygl;
using namespace nlohmann;
using namespace std::literals;

bool is_cmd(const std::vector<std::string>& tokens, int i) {
    auto& tok = tokens.at(i);
    return !(tok[0] == '[' || tok[0] == ']' || tok[0] == '\"' ||
             tok[0] == '-' || tok[0] == '+' || std::isdigit(tok[0]));
}

bool is_number(const std::vector<std::string>& tokens, int i) {
    auto& tok = tokens.at(i);
    return tok[0] == '-' || tok[0] == '+' || std::isdigit(tok[0]);
}

std::string parse_string(const std::vector<std::string>& tokens, int& i) {
    if (tokens[i][0] != '"') throw std::runtime_error("string expected");
    auto tok = tokens[i++];
    tok = tok.substr(1, tok.size() - 2);
    if (tok.find('|') != tok.npos) tok = tok.substr(tok.find('|') + 1);
    return tok;
}

void parse_param(const std::vector<std::string>& tokens, int& i, json& js) {
    auto open = false, first = true;
    while (i < tokens.size()) {
        if (is_cmd(tokens, i)) {
            break;
        } else if (tokens[i][0] == '[') {
            open = true;
            i++;
        } else if (tokens[i][0] == ']') {
            open = false;
            i++;
            break;
        } else if (tokens[i][0] == '"') {
            if (!first || !open) break;
            first = false;
            js.push_back(tokens[i].substr(1, tokens[i].size() - 2));
            i++;
        } else {
            js.push_back(atof(tokens[i].c_str()));
            i++;
        }
    }
}

void parse_param_list(
    const std::vector<std::string>& tokens, int& i, json& js) {
    while (i < tokens.size()) {
        if (is_cmd(tokens, i)) break;
        auto name = parse_string(tokens, i);
        js[name] = json::array();
        parse_param(tokens, i, js.at(name));
        if (js.at(name).size() == 1) { js.at(name) = js.at(name).at(0); }
    }
}

void parse_param_numbers(
    const std::vector<std::string>& tokens, int& i, json& js) {
    js["values"] = json::array();
    parse_param(tokens, i, js.at("values"));
}

json pbrt_to_json(const std::string& filename) {
    auto pbrt = load_text(filename);
    auto re = std::regex("\"(\\w+)\\s+(\\w+)\"");
    auto tokens = split(std::regex_replace(pbrt, re, "\"$1|$2\""));
    auto js = json::array();
    auto i = 0;
    while (i < tokens.size()) {
        if (!is_cmd(tokens, i)) throw std::runtime_error("command expected");
        auto& tok = tokens[i++];
        auto jcmd = json::object();
        jcmd["cmd"] = tok;
        if (tok == "Transform") {
            parse_param_numbers(tokens, i, jcmd);
        } else if (tok == "Integrator" || tok == "Sampler" ||
                   tok == "PixelFilter" || tok == "Film" || tok == "Camera" ||
                   tok == "Shape" || tok == "AreaLightSource" ||
                   tok == "LightSource") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Texture") {
            jcmd["name"] = parse_string(tokens, i);
            jcmd["value_type"] = parse_string(tokens, i);
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "MakeNamedMaterial") {
            jcmd["name"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "Material") {
            jcmd["type"] = parse_string(tokens, i);
            parse_param_list(tokens, i, jcmd);
        } else if (tok == "NamedMaterial") {
            jcmd["name"] = parse_string(tokens, i);
        } else if (tok == "WorldBegin" || tok == "AttributeBegin" ||
                   tok == "TransformBegin" || tok == "WorldEnd" ||
                   tok == "AttributeEnd" || tok == "TransformEnd") {
        } else {
            throw std::runtime_error("unsupported command " + tok);
        }
        js.push_back(jcmd);
    }
    save_text(filename + ".json", js.dump(2));
    return js;
}

void load_ply(const std::string& filename, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord) {
    auto f = fopen(filename.c_str(), "rb");
    if (!f) throw std::runtime_error("cannot open file " + filename);

    auto nverts = 0, nfaces = 0;
    auto vertex_pos = std::map<std::string, int>{};
    auto vert_size = 0;

    auto in_verts = false, in_faces = false;
    char buf[4096];
    while (fgets(buf, 4096, f)) {
        auto line = std::string(buf);
        auto toks = split(line);
        if (toks[0] == "ply") {
        } else if (toks[0] == "end_header") {
            break;
        } else if (toks[0] == "format") {
            if (toks[1] != "binary_little_endian")
                throw std::runtime_error("bad ply format");
        } else if (toks[0] == "element") {
            if (toks[1] == "vertex") {
                nverts = atoi(toks[2].c_str());
                in_verts = true;
                in_faces = false;
            } else if (toks[1] == "face") {
                nfaces = atoi(toks[2].c_str());
                in_verts = false;
                in_faces = true;
            } else {
                throw std::runtime_error("bad ply element");
            }
        } else if (toks[0] == "property") {
            if (toks[1] == "float") {
                if (!in_verts)
                    throw std::runtime_error("bad ply vertex property");
                vertex_pos[toks[2]] = vert_size++;
            } else if (toks[1] == "list") {
                if (!in_faces || toks[2] != "uint8" || toks[3] != "int" ||
                    toks[4] != "vertex_indices")
                    throw std::runtime_error("bad ply face property");
            } else {
                throw std::runtime_error("bad ply property");
            }
        } else {
            throw std::runtime_error("bad ply header");
        }
    }

    in_verts = true;
    in_faces = false;

    pos.clear();
    norm.clear();
    texcoord.clear();
    for (auto i = 0; i < nverts; i++) {
        float buf[32];
        fread(buf, sizeof(float), vert_size, f);
        if (contains(vertex_pos, "x"s) && contains(vertex_pos, "y"s) &&
            contains(vertex_pos, "z"s)) {
            pos.push_back({buf[vertex_pos["x"]], buf[vertex_pos["y"]],
                buf[vertex_pos["z"]]});
        }
        if (contains(vertex_pos, "nx"s) && contains(vertex_pos, "ny"s) &&
            contains(vertex_pos, "nz"s)) {
            norm.push_back({buf[vertex_pos["nx"]], buf[vertex_pos["ny"]],
                buf[vertex_pos["nz"]]});
        }
        if (contains(vertex_pos, "u"s) && contains(vertex_pos, "v"s)) {
            texcoord.push_back({buf[vertex_pos["u"]], buf[vertex_pos["v"]]});
        }
    }
    triangles.resize(nfaces);
    for (auto i = 0; i < nfaces; i++) {
        auto n = uint8_t(0);
        fread(&n, 1, 1, f);
        if (n != 3) throw std::runtime_error("bad ply face");
        fread(&triangles[i], sizeof(vec3i), 1, f);
    }

    fclose(f);
}

scene* load_pbrt(const std::string& filename) {
    auto js = pbrt_to_json(filename);
    auto dirname = path_dirname(filename);

    struct stack_item {
        frame3f frame = identity_frame3f;
        material* mat = nullptr;
        material* light_mat = nullptr;
    };

    // parse
    auto scn = new scene();
    auto stack = std::vector<stack_item>();
    stack.push_back(stack_item());
    auto txt_map = std::map<std::string, texture*>();
    auto mat_map = std::map<std::string, material*>();
    auto mid = 0;

    auto get_vec3f = [](const json& js) -> vec3f {
        if (js.is_number())
            return {js.get<float>(), js.get<float>(), js.get<float>()};
        if (js.is_array() && js.size() == 3)
            return {js.at(0).get<float>(), js.at(1).get<float>(),
                js.at(2).get<float>()};
        log_error("cannot handle vec3f");
        return zero3f;
    };

    auto get_mat4f = [](const json& js) -> mat4f {
        if (!js.is_array() || js.size() != 16) {
            log_error("cannot handle mat4f");
            return identity_mat4f;
        }
        auto m = identity_mat4f;
        for (auto i = 0; i < 16; i++) (&m.x.x)[i] = js.at(i).get<float>();
        return m;
    };

    auto get_vector_vec3i = [](const json& js) -> std::vector<vec3i> {
        if (!js.is_array() || js.size() % 3) {
            log_error("cannot handle vector<vec3f>");
            return {};
        }
        auto vals = std::vector<vec3i>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = (int)std::round(js.at(i * 3 + 0).get<float>());
            vals[i].y = (int)std::round(js.at(i * 3 + 1).get<float>());
            vals[i].z = (int)std::round(js.at(i * 3 + 2).get<float>());
        }
        return vals;
    };

    auto get_vector_vec3f = [](const json& js) -> std::vector<vec3f> {
        if (!js.is_array() || js.size() % 3) {
            log_error("cannot handle vector<vec3f>");
            return {};
        }
        auto vals = std::vector<vec3f>(js.size() / 3);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 3 + 0).get<float>();
            vals[i].y = js.at(i * 3 + 1).get<float>();
            vals[i].z = js.at(i * 3 + 2).get<float>();
        }
        return vals;
    };

    auto get_vector_vec2f = [](const json& js) -> std::vector<vec2f> {
        if (!js.is_array() || js.size() % 2) {
            log_error("cannot handle vector<vec3f>");
            return {};
        }
        auto vals = std::vector<vec2f>(js.size() / 2);
        for (auto i = 0; i < vals.size(); i++) {
            vals[i].x = js.at(i * 2 + 0).get<float>();
            vals[i].y = js.at(i * 2 + 1).get<float>();
        }
        return vals;
    };

    auto get_scaled_texture =
        [&txt_map, &get_vec3f](const json& js) -> std::pair<vec3f, texture*> {
        if (js.is_string())
            return {{1, 1, 1}, txt_map.at(js.get<std::string>())};
        return {get_vec3f(js), nullptr};
    };

    auto lid = 0;
    for (auto& jcmd : js) {
        auto cmd = jcmd.at("cmd").get<std::string>();
        if (cmd == "Integrator" || cmd == "Sampler" || cmd == "PixelFilter") {
        } else if (cmd == "Transform") {
            auto m = get_mat4f(jcmd.at("values"));
            stack.back().frame = stack.back().frame * mat_to_frame(m);
        } else if (cmd == "Film") {
            if (scn->cameras.empty()) {
                auto cam = new camera();
                cam->name = "cam";
                scn->cameras.push_back(cam);
            }
            scn->cameras.back()->aspect = jcmd.at("xresolution").get<float>() /
                                          jcmd.at("yresolution").get<float>();
        } else if (cmd == "Camera") {
            if (scn->cameras.empty()) {
                auto cam = new camera();
                cam->name = "cam";
                scn->cameras.push_back(cam);
            }
            // scn->cameras.back()->frame = stack.back().frame;
            scn->cameras.back()->frame = inverse(stack.back().frame);
            scn->cameras.back()->frame.z = -scn->cameras.back()->frame.z;
            auto type = jcmd.at("type").get<std::string>();
            if (type == "perspective") {
                scn->cameras.back()->yfov =
                    jcmd.at("fov").get<float>() * pif / 180;
            } else {
                log_error("{} camera not supported", type);
            }
        } else if (cmd == "Texture") {
            auto txt = new texture();
            scn->textures.push_back(txt);
            txt->name = jcmd.at("name").get<std::string>();
            txt_map[txt->name] = txt;
            auto type = jcmd.at("type").get<std::string>();
            if (type == "imagemap") {
                txt->path = jcmd.at("filename").get<std::string>();
            } else {
                log_error("{} texture not supported", type);
            }
        } else if (cmd == "MakeNamedMaterial" || cmd == "Material") {
            auto mat = new material();
            scn->materials.push_back(mat);
            if (cmd == "Material") {
                mat->name = "unnamed_mat" + std::to_string(mid++);
                stack.back().mat = mat;
            } else {
                mat->name = jcmd.at("name").get<std::string>();
                mat_map[mat->name] = mat;
            }
            auto type = jcmd.at("type").get<std::string>();
            if (type == "uber") {
                if (jcmd.count("Kd"))
                    std::tie(mat->kd, mat->kd_txt.txt) =
                        get_scaled_texture(jcmd.at("Kd"));
                if (jcmd.count("Ks"))
                    std::tie(mat->ks, mat->ks_txt.txt) =
                        get_scaled_texture(jcmd.at("Ks"));
                mat->rs = 0;
            } else if (type == "matte") {
                mat->kd = {1,1,1};
                if (jcmd.count("Kd"))
                    std::tie(mat->kd, mat->kd_txt.txt) =
                        get_scaled_texture(jcmd.at("Kd"));
                mat->rs = 1;
            } else if (type == "mirror") {
                mat->ks = {1,1,1};
                mat->rs = 0;
            } else if (type == "metal") {
                auto eta = get_vec3f(jcmd.at("eta"));
                auto k = get_vec3f(jcmd.at("eta"));
                mat->ks = fresnel_metal(1,eta, k);
                mat->rs = 0;
            } else if (type == "substrate") {
                if (jcmd.count("Kd"))
                    std::tie(mat->kd, mat->kd_txt.txt) =
                        get_scaled_texture(jcmd.at("Kd"));
                mat->ks = {0.04,0.04,0.04};
                if (jcmd.count("Ks"))
                    std::tie(mat->ks, mat->ks_txt.txt) =
                        get_scaled_texture(jcmd.at("Ks"));
                mat->rs = 0;
            } else if (type == "glass") {
                mat->ks = {0.04,0.04,0.04};
                if (jcmd.count("Ks"))
                    std::tie(mat->ks, mat->ks_txt.txt) =
                        get_scaled_texture(jcmd.at("Ks"));
                mat->rs = 0;
            } else {
                mat->kd = {1,0,0};
                log_error("{} material not supported", type);
            }
            if(jcmd.count("uroughness")) {
                auto remap = js.count("remaproughness") && js.at("remaproughness").get<bool>();
                if (jcmd.count("uroughness"))
                    mat->rs = jcmd.at("uroughness").get<float>();
                if(!remap) mat->rs = sqrt(mat->rs);
            }
            if(stack.back().light_mat) {
                mat->ke = stack.back().light_mat->ke;
                mat->ke_txt = stack.back().light_mat->ke_txt;
            }
        } else if (cmd == "NamedMaterial") {
            stack.back().mat = mat_map.at(jcmd.at("name").get<std::string>());
            if(stack.back().light_mat) {
                auto mat = new material(*stack.back().mat);
                mat->name += "_" + std::to_string(lid++);
                mat->ke = stack.back().light_mat->ke;
                mat->ke_txt = stack.back().light_mat->ke_txt;
                scn->materials.push_back(mat);
                stack.back().mat = mat;
            }
        } else if (cmd == "Shape") {
            auto shp = new shape();
            scn->shapes.push_back(shp);
            shp->frame = stack.back().frame;
            if (stack.back().mat) {
                shp->groups.push_back({});
                shp->groups.back().mat = stack.back().mat;
            }
            auto type = jcmd.at("type").get<std::string>();
            if (type == "plymesh") {
                auto filename = jcmd.at("filename").get<std::string>();
                shp->name = path_basename(filename);
                load_ply(dirname + filename, shp->triangles, shp->pos,
                    shp->norm, shp->texcoord);
            } else if (type == "trianglemesh") {
                if (jcmd.count("indices"))
                    shp->triangles = get_vector_vec3i(jcmd.at("indices"));
                if (jcmd.count("P")) shp->pos = get_vector_vec3f(jcmd.at("P"));
                if (jcmd.count("N")) shp->norm = get_vector_vec3f(jcmd.at("N"));
                if (jcmd.count("uv"))
                    shp->texcoord = get_vector_vec2f(jcmd.at("uv"));
            } else {
                log_error("{} shape not supported", type);
            }
        } else if (cmd == "AreaLightSource") {
            auto type = jcmd.at("type").get<std::string>();
            if(type == "diffuse") {
            auto lmat = new material();
            lmat->ke = get_vec3f(jcmd.at("L"));
            stack.back().light_mat = lmat;
            } else {
                log_error("{} area light not supported", type);
            }
        } else if (cmd == "LightSource") {
            auto type = jcmd.at("type").get<std::string>();
            if(type == "infinite") {   
                auto env = new environment();
                env->name = "env" + std::to_string(lid++);
                env->frame = stack.back().frame;
                env->ke = {1,1,1};
                if(jcmd.count("mapname")) {
                    auto txt = new texture();
                    txt->path = jcmd.at("mapname").get<std::string>();
                    txt->name = env->name;
                    scn->textures.push_back(txt);
                    env->ke_txt.txt = txt;
                }
                scn->environments.push_back(env);         
            } else {
                log_error("{} light not supported", type);
            }
        } else if (cmd == "WorldBegin") {
            stack.push_back(stack_item());
        } else if (cmd == "AttributeBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "TransformBegin") {
            stack.push_back(stack.back());
        } else if (cmd == "WorldEnd" || cmd == "AttributeEnd" ||
                   cmd == "TransformEnd") {
            stack.pop_back();
        } else {
            log_error("{} command not supported", cmd);
        }
    }
    return scn;
}

int main(int argc, char** argv) {
    // parse command line
    auto parser =
        ygl::make_parser(argc, argv, "ypbrt", "convert pbrt files to yocto");
    auto outfilename = ygl::parse_opt(
        parser, "--output", "-o", "output scene filename", "out.obj"s);
    auto filename =
        ygl::parse_arg(parser, "scene", "input scene filenames", ""s);
    if (ygl::should_exit(parser)) {
        printf("%s\n", get_usage(parser).c_str());
        exit(1);
    }

    // load image
    auto scn = load_pbrt(filename);

    // save scene
    system(("mkdir -p " + path_dirname(outfilename)).c_str());
    auto opts = save_options();
    opts.save_textures = false;
    save_scene(outfilename, scn, opts);

    return 0;
}
