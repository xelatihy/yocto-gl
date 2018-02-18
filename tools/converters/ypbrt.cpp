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

struct stack_node {
    mat4f xform = identity_mat4f;
    std::string matname = "";
};

bool is_cmd(const std::vector<std::string>& tokens, int i) {
    auto& tok = tokens.at(i);
    return !(tok[0] == '[' || tok[0] == ']' || tok[0] == '\"' ||
             tok[0] == '-' || tok[0] == '+' || std::isdigit(tok[0]));
}

std::string parse_string(const std::vector<std::string>& tokens, int& i) {
    if (tokens[i][0] != '"') throw std::runtime_error("string expected");
    return tokens[i++];
}

using param = std::pair<std::vector<std::string>, std::vector<double>>;

param parse_param(const std::vector<std::string>& tokens, int& i) {
    auto open = false, first = false;
    auto ret = std::pair<std::vector<std::string>, std::vector<double>>();
    while (i < tokens.size()) {
        if (is_cmd(tokens, i)) {
            break;
        } else if (tokens[i][0] == '[') {
            open = true;
        } else if (tokens[i][0] == ']') {
            open = false;
            break;
        } else if (tokens[i][0] == '"') {
            if (!first || !open) break;
            first = false;
            ret.first.push_back(tokens[i]);
        } else {
            ret.second.push_back(atof(tokens[i].c_str()));
        }
        i++;
    }
    return ret;
}

std::unordered_map<std::string, param> parse_param_list(
    const std::vector<std::string>& tokens, int& i) {
    auto ret = std::unordered_map<std::string, param>();
    while (i < tokens.size()) {
        if (is_cmd(tokens, i)) break;
        auto name = parse_string(tokens, i);
        auto param = parse_param(tokens, i);
        ret[name] = param;
    }
    return ret;
}

void pbrt_to_json(const std::string& filename) {
    auto pbrt = load_text(filename);
    save_text(filename + ".txt", pbrt);
    auto re = std::regex("\"(\\w+)\\s+(\\w+)\"");
    auto tokens = split(std::regex_replace(pbrt, re, "\"$1|$2\""));
    auto jpbrt = ""s;
    auto stack = std::vector<stack_node>();
    for (auto i = 0; i < tokens.size();) {
        if (!is_cmd(tokens, i)) throw std::runtime_error("command expected");
        auto& tok = tokens[i++];
        if (tok == "Transform") { parse_param(tokens, i); }
        if (tok == "Integrator" || tok == "Sampler" || tok == "PixelFilter" ||
            tok == "Film" || tok == "Camera" || tok == "Shape" ||
            tok == "AreaLightSource") {
            auto type = parse_string(tokens, i);
            auto params = parse_param_list(tokens, i);
            jpbrt += tok + " " + type + "\n";
        }
    }
    save_text(filename + ".j.txt", jpbrt);
    // return js;
}

scene* load_pbrt(const std::string& filename) {
    pbrt_to_json(filename);
    return nullptr;
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
    //    system(("mkdir -p " + path_dirname(outfilename)).c_str());
    //    auto opts = save_options();
    //    opts.save_textures = false;
    //    save_scene(outfilename, scn, opts);

    return 0;
}
