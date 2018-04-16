//
// Implementation for Yocto/Obj. See yocto_gl.h for documentation.
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

#include "yocto_obj.h"
#include "yocto_image.h"
#include "yocto_utils.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

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
    val = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
}

// parse base value
inline void obj_parse_base(char*& s, float& val) {
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
inline void obj_parse(char*& s, int& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, float& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, bool& val) {
    auto i = 0;
    obj_parse(s, i);
    val = i;
}
inline void obj_parse(char*& s, std::string& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, char* val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, vec2f& val) {
    for (auto i = 0; i < 2; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec3f& val) {
    for (auto i = 0; i < 3; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec4f& val) {
    for (auto i = 0; i < 4; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, frame3f& val) {
    for (auto i = 0; i < 12; i++) obj_parse(s, (&val.x.x)[i]);
}
inline void obj_parse(char*& s, obj_vertex& val, const obj_vertex& vert_size) {
    char buf[1024];
    obj_skipws(s);
    obj_parse_base(s, buf);
    val = obj_vertex{-1, -1, -1, -1, -1};
    auto v = &val.pos;
    auto vs = &vert_size.pos;
    auto i = 0;
    auto sb = buf;
    while (i < 5 && *sb) {
        obj_parse_base(sb, v[i]);
        v[i] = (v[i] < 0) ? vs[i] + v[i] : v[i] - 1;
        if (*sb != '/') break;
        while (*sb == '/') {
            sb++;
            i++;
        }
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
            info.props[last].push_back(tok);
        }
    }

    // clamp
    if (info.props.find("-clamp") != info.props.end() &&
        info.props.at("-clamp").size() > 0) {
        auto& clamp_vec = info.props.at("-clamp");
        auto clamp_str = (clamp_vec.empty()) ? "" : clamp_vec.front();
        info.clamp = clamp_str == "on" || clamp_str == "1";
        info.props.erase("-clamp");
    }

    if (info.props.find("-bm") != info.props.end() &&
        info.props.at("-bm").size() > 0) {
        auto& bm_vec = info.props.at("-bm");
        auto bm_str = (bm_vec.empty()) ? "" : bm_vec.front();
        info.scale = std::atof(bm_str.c_str());
        info.props.erase("-bm");
    }
}

// Load MTL
std::vector<obj_material*> load_mtl(const std::string& filename, bool flip_tr,
    std::vector<std::string>& textures) {
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
            materials.push_back(new obj_material());
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
        } else if (obj_streq(cmd, "map_d")) {
            obj_parse(ss, mat->op_txt);
        } else if (obj_streq(cmd, "map_Ni")) {
            obj_parse(ss, mat->ior_txt);
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

    // create texture array
    textures = {};
    auto texture_set = std::unordered_set<std::string>();
    auto add_texture = [&texture_set, &textures](const obj_texture_info& info) {
        if (info.path == "") return;
        if (texture_set.find(info.path) != texture_set.end()) return;
        texture_set.insert(info.path);
        textures.push_back(info.path);
    };
    for (auto mat : materials) {
        add_texture(mat->ke_txt);
        add_texture(mat->ka_txt);
        add_texture(mat->kd_txt);
        add_texture(mat->ks_txt);
        add_texture(mat->kr_txt);
        add_texture(mat->kt_txt);
        add_texture(mat->ns_txt);
        add_texture(mat->op_txt);
        add_texture(mat->ior_txt);
        add_texture(mat->bump_txt);
        add_texture(mat->bump_txt);
        add_texture(mat->disp_txt);
        add_texture(mat->norm_txt);
    }

    // clone
    fclose(fs);

    // done
    return materials;
}

// Loads textures for an scene.
void load_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing) {
    for (auto txt : obj->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
#if YGL_IMAGEIO
        if (is_hdr_filename(filename)) {
            txt->dataf =
                load_imagef(filename, txt->width, txt->height, txt->ncomp);
        } else {
            txt->datab =
                load_imageb(filename, txt->width, txt->height, txt->ncomp);
        }
#endif
        if (txt->datab.empty() && txt->dataf.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}

// Loads an OBJ
obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool load_txt, bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // clear obj
    auto obj = new obj_scene();

    // splitting policy
    auto split_material = split_shapes;
    auto split_group = split_shapes;
    auto split_smoothing = split_shapes;

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // current parsing values
    auto mtllibs = std::vector<std::string>();
    auto oobj = (obj_object*)nullptr;
    auto matname = ""s;
    auto oname = ""s;
    auto faceted = false;
    auto elems = std::vector<obj_vertex>();

    // initializing obj
    obj->objects.push_back(new obj_object());
    oobj = obj->objects.back();
    oobj->groups.push_back({"", "", false});

    // keep track of array lengths
    auto vert_size = obj_vertex{0, 0, 0, 0, 0};

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
            vert_size.pos += 1;
            obj->pos.push_back(zero3f);
            obj_parse(ss, obj->pos.back());
        } else if (obj_streq(cmd, "vn")) {
            vert_size.norm += 1;
            obj->norm.push_back(zero3f);
            obj_parse(ss, obj->norm.back());
        } else if (obj_streq(cmd, "vt")) {
            vert_size.texcoord += 1;
            obj->texcoord.push_back(zero2f);
            obj_parse(ss, obj->texcoord.back());
            if (flip_texcoord)
                obj->texcoord.back().y = 1 - obj->texcoord.back().y;
        } else if (obj_streq(cmd, "vc")) {
            vert_size.color += 1;
            obj->color.push_back(vec4f{0, 0, 0, 1});
            obj_parse(ss, obj->color.back());
        } else if (obj_streq(cmd, "vr")) {
            vert_size.radius += 1;
            obj->radius.push_back(0);
            obj_parse(ss, obj->radius.back());
        } else if (obj_streq(cmd, "f") || obj_streq(cmd, "l") ||
                   obj_streq(cmd, "p") || obj_streq(cmd, "b")) {
            auto elem = obj_element();
            elem.type = elem_type_map.at(cmd);
            elem.start = (uint32_t)oobj->verts.size();
            elem.size = 0;
            elem.groupid = (int)oobj->groups.size() - 1;
            oobj->elems.push_back(elem);
            obj_skipws(ss);
            while (*ss) {
                auto vert = obj_vertex();
                obj_parse(ss, vert, vert_size);
                obj_skipws(ss);
                oobj->verts.push_back(vert);
                oobj->elems.back().size += 1;
            }
        } else if (obj_streq(cmd, "o")) {
            auto name = ""s;
            obj_parse(ss, name);
            obj->objects.push_back(new obj_object());
            oobj = obj->objects.back();
            oobj->name = name;
            oobj->groups.push_back({"", "", false});
            matname = "";
            oname = name;
        } else if (obj_streq(cmd, "usemtl")) {
            obj_parse(ss, matname);
            if (split_material) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname;
                oobj->groups.push_back({"", matname, faceted});
            } else {
                oobj->groups.push_back({oobj->groups.back().name, matname,
                    oobj->groups.back().faceted});
            }
        } else if (obj_streq(cmd, "g")) {
            auto name = ""s;
            obj_parse(ss, name);
            if (split_group) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname + name;
                oobj->groups.push_back({name, matname, faceted});
            } else {
                oobj->groups.push_back({name, oobj->groups.back().matname,
                    oobj->groups.back().faceted});
            }
        } else if (obj_streq(cmd, "s")) {
            auto name = std::string();
            obj_parse(ss, name);
            faceted = (name != "on");
            if (split_smoothing) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname;
                oobj->groups.push_back({"", matname, faceted});
            } else {
                oobj->groups.push_back({oobj->groups.back().name,
                    oobj->groups.back().matname, faceted});
            }
        } else if (obj_streq(cmd, "op")) {
            auto name = std::string();
            obj_parse(ss, name);
            obj_skipws(ss);
            while (*ss) {
                auto tok = std::string();
                obj_parse(ss, tok);
                obj_skipws(ss);
                oobj->props[name].push_back(tok);
            }
        } else if (obj_streq(cmd, "mtllib")) {
            mtllibs.push_back("");
            obj_parse(ss, mtllibs.back());
        } else if (obj_streq(cmd, "c")) {
            auto cam = new obj_camera();
            obj_parse(ss, cam->name);
            obj_parse(ss, cam->ortho);
            obj_parse(ss, cam->yfov);
            obj_parse(ss, cam->aspect);
            obj_parse(ss, cam->aperture);
            obj_parse(ss, cam->focus);
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
        } else if (obj_streq(cmd, "n")) {
            auto nde = new obj_node();
            obj_parse(ss, nde->name);
            obj_parse(ss, nde->parent);
            obj_parse(ss, nde->camname);
            obj_parse(ss, nde->objname);
            obj_parse(ss, nde->envname);
            obj_parse(ss, nde->frame);
            obj_parse(ss, nde->translation);
            obj_parse(ss, nde->rotation);
            obj_parse(ss, nde->scale);
            if (nde->parent == "\"\"") nde->parent = "";
            if (nde->camname == "\"\"") nde->camname = "";
            if (nde->objname == "\"\"") nde->objname = "";
            if (nde->envname == "\"\"") nde->envname = "";
            obj->nodes.push_back(nde);
        } else {
            // unused
        }
    }

    // cleanup empty
    for (auto idx = 0; idx < obj->objects.size(); idx++) {
        if (!obj->objects[idx]->elems.empty()) continue;
        if (!obj->objects[idx]->verts.empty()) continue;
        delete obj->objects[idx];
        obj->objects.erase(obj->objects.begin() + idx);
        idx--;
    }

    // cleanup unused
    for (auto oobj : obj->objects) {
        if (oobj->groups.empty()) continue;
        auto used = std::vector<bool>(oobj->groups.size());
        for (auto& elem : oobj->elems) used[elem.groupid] = true;
        auto emap = std::vector<uint16_t>(oobj->groups.size());
        auto valid = 0;
        for (auto idx = 0; idx < oobj->groups.size(); idx++) {
            emap[idx] = valid;
            if (used[idx]) valid++;
        }
        for (auto& elem : oobj->elems) elem.groupid = emap[elem.groupid];
        for (auto idx = 0; idx < used.size(); idx++) {
            if (used[idx]) continue;
            used.erase(used.begin() + idx);
            oobj->groups.erase(oobj->groups.begin() + idx);
            idx--;
        }
        if (oobj->groups.size() == 1 && oobj->groups[0].name == "" &&
            oobj->groups[0].matname == "" && oobj->groups[0].faceted == false)
            oobj->groups.clear();
    }

    auto end = std::remove_if(obj->objects.begin(), obj->objects.end(),
        [](const obj_object* x) { return !x; });
    obj->objects.erase(end, obj->objects.end());

    // parse materials
    auto mtllibs_set =
        std::unordered_set<std::string>(mtllibs.begin(), mtllibs.end());
    mtllibs = std::vector<std::string>{mtllibs_set.begin(), mtllibs_set.end()};
    auto dirname = path_dirname(filename);
    std::unordered_set<std::string> texture_set;
    for (auto mtllib : mtllibs) {
        auto mtlname = dirname + mtllib;
        std::vector<std::string> textures;
        auto materials = load_mtl(mtlname, flip_tr, textures);
        obj->materials.insert(
            obj->materials.end(), materials.begin(), materials.end());
        for (auto& txt : textures) {
            if (texture_set.find(txt) != texture_set.end()) continue;
            obj->textures.push_back(new obj_texture());
            obj->textures.back()->path = txt;
            texture_set.insert(txt);
        }
    }
    for (auto oenv : obj->environments) {
        auto txt = oenv->ke_txt.path;
        if (txt == "") continue;
        if (texture_set.find(txt) != texture_set.end()) continue;
        obj->textures.push_back(new obj_texture());
        obj->textures.back()->path = txt;
        texture_set.insert(txt);
    }

    // load textures
    if (load_txt) load_textures(obj, dirname, skip_missing);

    // close file
    fclose(fs);

    // done
    return obj;
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
inline void obj_dump(char*& s, const float* val, int num) {
    for (auto i = 0; i < num; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, val[i]);
    }
}
inline void obj_dump(char*& s, const vec2f& val) { obj_dump(s, &val.x, 2); }
inline void obj_dump(char*& s, const vec3f& val) { obj_dump(s, &val.x, 3); }
inline void obj_dump(char*& s, const vec4f& val) { obj_dump(s, &val.x, 4); }
inline void obj_dump(char*& s, const frame3f& val) {
    obj_dump(s, &val.x.x, 12);
}
inline void obj_dump(char*& s, const obj_vertex& val) {
    auto vert_ptr = &val.pos;
    auto nto_write = 0;
    for (auto i = 0; i < 5; i++) {
        if (vert_ptr[i] >= 0) nto_write = i + 1;
    }
    for (auto i = 0; i < nto_write; i++) {
        if (i) *s++ = '/';
        if (vert_ptr[i] >= 0) s += sprintf(s, "%d", vert_ptr[i] + 1);
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
inline void obj_dump_line(
    FILE* fs, const char* lbl, const obj_vertex* vals, int count) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, lbl);
    for (auto i = 0; i < count; i++) {
        *s++ = ' ';
        obj_dump(s, vals[i]);
    }
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

// Save an MTL file
void save_mtl(const std::string& filename,
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
        // if (mat->kt != zero3f) obj_dump_line(fs, "  Tf", mat->kt);
        if (mat->ns != 0.0f)
            obj_dump_line(fs, "  Ns", (int)clamp(mat->ns, 0.0f, 1000000000.0f));
        if (mat->op != 1.0f) obj_dump_line(fs, "  d", mat->op);
        if (mat->ior != 1.0f) obj_dump_line(fs, "  Ni", mat->ior);
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
        if (mat->bump_txt.path != "")
            obj_dump_line(fs, "  map_bump", mat->bump_txt);
        if (mat->disp_txt.path != "")
            obj_dump_line(fs, "  map_disp", mat->disp_txt);
        if (mat->norm_txt.path != "")
            obj_dump_line(fs, "  map_norm", mat->norm_txt);
        for (auto&& kv : mat->props) {
            obj_dump_line(fs, ("  " + kv.first + ' ').c_str(), kv.second);
        }
        fputs("\n", fs);
    }

    fclose(fs);
}

// Loads textures for an scene.
void save_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing) {
    for (auto txt : obj->textures) {
        if (txt->datab.empty() && txt->dataf.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
#if YGL_IMAGEIO
        if (!txt->datab.empty()) {
            ok = save_imageb(filename, txt->width, txt->height, txt->ncomp,
                txt->datab.data());
        }
        if (!txt->dataf.empty()) {
            ok = save_imagef(filename, txt->width, txt->height, txt->ncomp,
                txt->dataf.data());
        }
#endif
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
}

// Save an OBJ
void save_obj(const std::string& filename, const obj_scene* obj, bool save_txt,
    bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // linkup to mtl
    auto dirname = path_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!obj->materials.empty()) {
        obj_dump_line(fs, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto cam : obj->cameras) {
        obj_dump_sp(fs, "c");
        obj_dump_sp(fs, cam->name);
        obj_dump_sp(fs, cam->ortho);
        obj_dump_sp(fs, cam->yfov);
        obj_dump_sp(fs, cam->aspect);
        obj_dump_sp(fs, cam->aperture);
        obj_dump_sp(fs, cam->focus);
        obj_dump_sp(fs, cam->frame);
        obj_dump_sp(fs, "\n");
    }

    // save envs
    for (auto env : obj->environments) {
        obj_dump_sp(fs, "e");
        obj_dump_sp(fs, env->name);
        obj_dump_sp(fs, env->ke);
        obj_dump_sp(fs, (env->ke_txt.path != "") ? env->ke_txt.path : "\"\""s);
        obj_dump_sp(fs, env->frame);
        obj_dump_sp(fs, "\n");
    }

    // save nodes
    for (auto nde : obj->nodes) {
        obj_dump_sp(fs, "n");
        obj_dump_sp(fs, nde->name);
        obj_dump_sp(fs, (nde->parent.empty()) ? "\"\""s : nde->parent);
        obj_dump_sp(fs, (nde->camname.empty()) ? "\"\""s : nde->camname);
        obj_dump_sp(fs, (nde->objname.empty()) ? "\"\""s : nde->objname);
        obj_dump_sp(fs, (nde->envname.empty()) ? "\"\""s : nde->envname);
        obj_dump_sp(fs, nde->frame);
        obj_dump_sp(fs, nde->translation);
        obj_dump_sp(fs, nde->rotation);
        obj_dump_sp(fs, nde->scale);
        obj_dump_sp(fs, "\n");
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
    for (auto& v : obj->color) obj_dump_line(fs, "vc", v);
    for (auto& v : obj->radius) obj_dump_line(fs, "vr", v);

    // save element data
    static auto elem_labels = std::unordered_map<obj_element_type, std::string>{
        {obj_element_type::point, "p"}, {obj_element_type::line, "l"},
        {obj_element_type::face, "f"}, {obj_element_type::bezier, "b"}};
    for (auto oobj : obj->objects) {
        obj_dump_line(fs, "o", oobj->name);
        for (auto& kv : oobj->props) {
            auto nv = kv.second;
            nv.insert(nv.begin(), kv.first);
            obj_dump_line(fs, "op", nv);
        }
        if (!oobj->groups.empty()) {
            auto& ogrp = oobj->groups[0];
            if (ogrp.name != "") obj_dump_line(fs, "g", ogrp.name);
            if (ogrp.matname != "") obj_dump_line(fs, "usemtl", ogrp.matname);
            if (ogrp.faceted) obj_dump_line(fs, "s", "off");
        }
        auto last_groupid = 0;
        for (auto& elem : oobj->elems) {
            if (last_groupid != elem.groupid) {
                auto& lgrp = oobj->groups[last_groupid];
                auto& ogrp = oobj->groups[elem.groupid];
                if (ogrp.name != lgrp.name) obj_dump_line(fs, "g", ogrp.name);
                if (ogrp.matname != "" && ogrp.matname != lgrp.matname)
                    obj_dump_line(fs, "usemtl", ogrp.matname);
                if (ogrp.faceted != lgrp.faceted)
                    obj_dump_line(fs, "s", ogrp.faceted ? "off" : "on");
                last_groupid = elem.groupid;
            }
            obj_dump_line(fs, elem_labels.at(elem.type).c_str(),
                oobj->verts.data() + elem.start, elem.size);
        }
    }

    fclose(fs);

    // save materials
    if (!obj->materials.empty())
        save_mtl(dirname + basename + ".mtl", obj->materials, flip_tr);

    // save textures
    if (save_txt) save_textures(obj, dirname, skip_missing);
}

}  // namespace ygl
