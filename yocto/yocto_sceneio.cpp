//
// Implementation for Yocto/Scene.
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

#include "yocto_sceneio.h"
#include "yocto_image.h"
#include "yocto_shape.h"
#include "yocto_utils.h"

#include <array>
#include <deque>
#include <unordered_map>

#include <fstream>
#include "ext/json.hpp"
#include "yocto_gltf.h"

// -----------------------------------------------------------------------------
// GENERIC METHOD
// -----------------------------------------------------------------------------
namespace ygl {

// Load a scene
scene* load_scene(const std::string& filename, bool skip_missing) {
    auto ext = path_extension(filename);
    auto scn = (scene*)nullptr;
    if (ext == ".json" || ext == ".JSON") {
        scn = load_json_scene(filename);
    } else if (ext == ".obj" || ext == ".OBJ") {
        scn = load_obj_scene(filename, skip_missing);
    } else if (ext == ".gltf" || ext == ".GLTF") {
        scn = load_gltf_scene(filename, skip_missing);
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
    return scn;
}

// Save a scene
void save_scene(
    const std::string& filename, const scene* scn, bool skip_missing) {
    auto ext = path_extension(filename);
    if (ext == ".json" || ext == ".JSON") {
        save_json_scene(filename, scn, skip_missing);
    } else if (ext == ".obj" || ext == ".OBJ") {
        save_obj_scene(filename, scn, skip_missing);
    } else if (ext == ".gltf" || ext == ".GLTF") {
        save_gltf_scene(filename, scn, false, skip_missing);
    } else {
        throw std::runtime_error("unsupported extension " + ext);
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// BUILTIN JSON FORMAT
// -----------------------------------------------------------------------------
namespace ygl {

// Json alias
using json = nlohmann::json;

// Parse int function.
void serialize(int& val, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_number_integer())
            throw std::runtime_error("integer expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse float function.
void serialize(float& val, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_number()) throw std::runtime_error("number expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse bool function.
void serialize(bool& val, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_boolean()) throw std::runtime_error("bool expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse string function.
void serialize(std::string& val, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_string()) throw std::runtime_error("string expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse json function.
void serialize(json& val, json& js, bool reading, void* ctx) {
    if (reading) {
        val = js;
    } else {
        js = val;
    }
}

// Parse support function.
template <typename T>
void serialize(T*& val, json& js, bool reading, void* ctx) {
    if (reading) {
        if (js.is_null()) {
            if (val) delete val;
            val = nullptr;
            return;
        }
        if (!js.is_object()) throw std::runtime_error("object expected");
        if (!val) val = new T();
        serialize(*val, js, reading, ctx);
    } else {
        if (!val) {
            js = nullptr;
            return;
        }
        if (!js.is_object()) js = json::object();
        serialize(*val, js, reading, ctx);
    }
}

// Parse support function.
template <typename T>
void serialize(std::vector<T>& vals, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        vals.resize(js.size());
        for (auto i = 0; i < js.size(); i++) {
            // this is contrived to support for std::vector<bool>
            auto v = T();
            serialize(v, js[i], reading, ctx);
            vals[i] = v;
        }
    } else {
        js = json::array();
        for (auto i = 0; i < vals.size(); i++)
            serialize(vals[i], js[i], reading, ctx);
    }
}

// Parse support function.
template <typename T, size_t N>
void serialize(std::array<T, N>& vals, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        if (N != js.size()) throw std::runtime_error("wrong array size");
        for (auto i = 0; i < N; i++) serialize(vals[i], js.at(i), reading, ctx);
    } else {
        js = json::array();
        for (auto i = 0; i < N; i++) serialize(vals[i], js[i], reading, ctx);
    }
}

// Parse support function.
template <typename T>
void serialize(
    std::map<std::string, T>& vals, json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
        for (auto kv = js.begin(); kv != js.end(); ++kv) {
            serialize(vals[kv.key()], kv.value(), reading, ctx);
        }
    } else {
        js = json::object();
        for (auto& kv : vals) serialize(kv.second, js[kv.first], reading, ctx);
    }
}

// Parse support function.
template <typename T, typename T1>
void serialize(
    T& val, json& js, bool reading, void* ctx, const std::map<T, T1>& table) {
    if (reading) {
        auto v = T1();
        serialize(v, js, reading, ctx);
        auto found = false;
        for (auto& kv : table) {
            if (kv.second == v) {
                val = kv.first;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("bad enum value");
    } else {
        auto found = false;
        auto v = T1();
        for (auto& kv : table) {
            if (kv.first == val) {
                v = kv.second;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("invalid value");
        serialize(v, js, reading, ctx);
    }
}

// Parse support function.
void serialize(vec2f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 2>&)vals, js, reading, ctx);
}
void serialize(vec3f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 3>&)vals, js, reading, ctx);
}
void serialize(vec4f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 4>&)vals, js, reading, ctx);
}
void serialize(vec2i& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<int, 2>&)vals, js, reading, ctx);
}
void serialize(vec3i& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<int, 3>&)vals, js, reading, ctx);
}
void serialize(vec4i& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<int, 4>&)vals, js, reading, ctx);
}
void serialize(mat3f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 9>&)vals, js, reading, ctx);
}
void serialize(mat4f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 16>&)vals, js, reading, ctx);
}
void serialize(frame3f& vals, json& js, bool reading, void* ctx) {
    serialize((std::array<float, 12>&)vals, js, reading, ctx);
}

// Parse support function.
void serialize_obj(json& js, bool reading, void* ctx) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
    } else {
        if (!js.is_object()) js = json::object();
    }
}

// Parse support function.
template <typename T>
void serialize_attr(T& val, json& js, const char* name, bool reading, void* ctx,
    bool required = true, const T& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading, ctx);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading, ctx);
        }
    } else {
        if (required || val != def) serialize(val, js[name], reading, ctx);
    }
}

// Dump support function.
template <typename T>
void serialize_attr(std::vector<T>& val, json& js, const char* name,
    bool reading, void* ctx, bool required = true,
    const std::vector<T>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading, ctx);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading, ctx);
        }
    } else {
        if (required || !val.empty()) serialize(val, js[name], reading, ctx);
    }
}

// Parse support function.
template <typename T>
void serialize_namedref(
    T*& val, json& js, bool reading, void* ctx, const std::vector<T*> values) {
    if (reading) {
        if (js.is_null()) {
            val = nullptr;
        } else if (js.is_string()) {
            auto name = js.get<std::string>();
            if (name == "") {
                val = nullptr;
            } else {
                for (auto v : values)
                    if (v->name == name) val = v;
                if (!val) throw std::runtime_error("undefined reference");
            }
        } else {
            throw std::runtime_error("string or null expected");
        }
    } else {
        if (val) {
            js = val->name;
        } else {
            js = json();
        }
    }
}

// Parse support function.
template <typename T>
void serialize_namedref_attr(T*& val, json& js, const char* name, bool reading,
    void* ctx, const std::vector<T*> values, bool required = true,
    T* def = nullptr) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize_namedref(val, js.at(name), reading, ctx, values);
        } else {
            if (!js.count(name)) {
                val = def;
            } else {
                serialize_namedref(val, js.at(name), reading, ctx, values);
            }
        }
    } else {
        if (required || val != def)
            serialize_namedref(val, js[name], reading, ctx, values);
    }
}

// Parse support function.
template <typename T>
void serialize_namedref(std::vector<T*>& vals, json& js, bool reading,
    void* ctx, const std::vector<T*> values) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        vals.resize(js.size());
        for (auto i = 0; i < js.size(); i++) {
            serialize_namedref(vals[i], js[i], reading, ctx, values);
        }
    } else {
        js = json::array();
        for (auto i = 0; i < vals.size(); i++)
            serialize_namedref(vals[i], js[i], reading, ctx, values);
    }
}

// Dump support function.
template <typename T>
void serialize_namedref_attr(std::vector<T*>& val, json& js, const char* name,
    bool reading, void* ctx, const std::vector<T*> values, bool required = true,
    const std::vector<T*>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize_namedref(val, js.at(name), reading, ctx, values);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize_namedref(val, js.at(name), reading, ctx, values);
        }
    } else {
        if (required || !val.empty())
            serialize_namedref(val, js[name], reading, ctx, values);
    }
}

// Serialize struct
void serialize(camera& val, json& js, bool reading, void* ctx) {
    static auto def = camera();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.frame, js, "frame", reading, ctx, false, def.frame);
    serialize_attr(val.ortho, js, "ortho", reading, ctx, false, def.ortho);
    serialize_attr(val.width, js, "width", reading, ctx, false, def.width);
    serialize_attr(val.height, js, "height", reading, ctx, false, def.height);
    serialize_attr(val.focal, js, "focal", reading, ctx, false, def.focal);
    serialize_attr(val.focus, js, "focus", reading, ctx, false, def.focus);
    serialize_attr(
        val.aperture, js, "aperture", reading, ctx, false, def.aperture);
    serialize_attr(val.near, js, "near", reading, ctx, false, def.near);
    serialize_attr(val.far, js, "far", reading, ctx, false, def.far);
}

// Serialize struct
void serialize(texture& val, json& js, bool reading, void* ctx) {
    static auto def = texture();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.path, js, "path", reading, ctx, false, def.path);
    serialize_attr(val.clamp, js, "clamp", reading, ctx, false, def.clamp);
    serialize_attr(val.scale, js, "scale", reading, ctx, false, def.scale);
    if (val.path == "") {
        serialize_attr(val.width, js, "width", reading, ctx, false, def.width);
        serialize_attr(
            val.height, js, "height", reading, ctx, false, def.height);
        serialize_attr(val.img, js, "img", reading, ctx, false, def.img);
    }
}

// Serialize struct
void serialize(material& val, json& js, bool reading, void* ctx) {
    static auto def = material();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.base_metallic, js, "base_metallic", reading, ctx, false,
        def.base_metallic);
    serialize_attr(val.gltf_textures, js, "gltf_textures", reading, ctx, false,
        def.gltf_textures);
    serialize_attr(val.double_sided, js, "double_sided", reading, ctx, false,
        def.double_sided);
    serialize_attr(val.ke, js, "ke", reading, ctx, false, def.ke);
    serialize_attr(val.kd, js, "kd", reading, ctx, false, def.kd);
    serialize_attr(val.ks, js, "ks", reading, ctx, false, def.ks);
    serialize_attr(val.kt, js, "kt", reading, ctx, false, def.kt);
    serialize_attr(val.rs, js, "rs", reading, ctx, false, def.rs);
    serialize_attr(val.op, js, "op", reading, ctx, false, def.op);
    serialize_attr(
        val.refract, js, "refract", reading, ctx, false, def.refract);
    serialize_namedref_attr(val.ke_txt, js, "ke_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.ke_txt);
    serialize_namedref_attr(val.kd_txt, js, "kd_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.kd_txt);
    serialize_namedref_attr(val.ks_txt, js, "ks_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.ks_txt);
    serialize_namedref_attr(val.kt_txt, js, "kt_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.kt_txt);
    serialize_namedref_attr(val.rs_txt, js, "rs_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.rs_txt);
    serialize_namedref_attr(val.op_txt, js, "op_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.op_txt);
    serialize_namedref_attr(val.occ_txt, js, "occ_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.occ_txt);
    serialize_namedref_attr(val.bump_txt, js, "bump_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.bump_txt);
    serialize_namedref_attr(val.disp_txt, js, "disp_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.disp_txt);
    serialize_namedref_attr(val.norm_txt, js, "norm_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.norm_txt);
}

// Serialize struct
void serialize(shape& val, json& js, bool reading, void* ctx) {
    static auto def = shape();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.path, js, "path", reading, ctx, false, def.path);
    if (val.path == "") {
        serialize_attr(val.pos, js, "pos", reading, ctx, false, def.pos);
        serialize_attr(val.norm, js, "norm", reading, ctx, false, def.norm);
        serialize_attr(
            val.texcoord, js, "texcoord", reading, ctx, false, def.texcoord);
        serialize_attr(val.color, js, "color", reading, ctx, false, def.color);
        serialize_attr(
            val.radius, js, "radius", reading, ctx, false, def.radius);
        serialize_attr(val.lines, js, "lines", reading, ctx, false, def.lines);
        serialize_attr(
            val.triangles, js, "triangles", reading, ctx, false, def.triangles);
    }
}

// Serialize struct
void serialize(subdiv& val, json& js, bool reading, void* ctx) {
    static auto def = subdiv();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.path, js, "path", reading, ctx, false, def.path);
    serialize_attr(val.level, js, "level", reading, ctx, false, def.level);
    serialize_attr(val.catmull_clark, js, "catmull_clark", reading, ctx, false,
        def.catmull_clark);
    serialize_attr(val.compute_normals, js, "compute_normals", reading, ctx,
        false, def.compute_normals);
    if (val.path == "") {
        serialize_attr(val.pos, js, "pos", reading, ctx, false, def.pos);
        serialize_attr(
            val.texcoord, js, "texcoord", reading, ctx, false, def.texcoord);
        serialize_attr(val.color, js, "color", reading, ctx, false, def.color);
        serialize_attr(
            val.quads_pos, js, "quads_pos", reading, ctx, false, def.quads_pos);
        serialize_attr(val.quads_texcoord, js, "quads_texcoord", reading, ctx,
            false, def.quads_texcoord);
        serialize_attr(val.quads_color, js, "quads_color", reading, ctx, false,
            def.quads_color);
    }
}

// Serialize struct
void serialize(instance& val, json& js, bool reading, void* ctx) {
    static auto def = instance();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.frame, js, "frame", reading, ctx, false, def.frame);
    serialize_namedref_attr(val.shp, js, "shp", reading, ctx,
        ((scene*)ctx)->shapes, false, def.shp);
    serialize_namedref_attr(val.mat, js, "mat", reading, ctx,
        ((scene*)ctx)->materials, false, def.mat);
    serialize_namedref_attr(val.sbd, js, "sbd", reading, ctx,
        ((scene*)ctx)->subdivs, false, def.sbd);
}

// Serialize struct
void serialize(environment& val, json& js, bool reading, void* ctx) {
    static auto def = environment();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.frame, js, "frame", reading, ctx, false, def.frame);
    serialize_attr(val.ke, js, "ke", reading, ctx, false, def.ke);
    serialize_namedref_attr(val.ke_txt, js, "ke_txt", reading, ctx,
        ((scene*)ctx)->textures, false, def.ke_txt);
}

// Serialize struct
void serialize(node& val, json& js, bool reading, void* ctx) {
    static auto def = node();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_namedref_attr(val.parent, js, "parent", reading, ctx,
        ((scene*)ctx)->nodes, false, def.parent);
    serialize_attr(val.frame, js, "frame", reading, ctx, false, def.frame);
    serialize_attr(val.translation, js, "translation", reading, ctx, false,
        def.translation);
    serialize_attr(
        val.rotation, js, "rotation", reading, ctx, false, def.rotation);
    serialize_attr(val.scale, js, "scale", reading, ctx, false, def.scale);
    serialize_attr(
        val.weights, js, "weights", reading, ctx, false, def.weights);
    serialize_namedref_attr(val.cam, js, "cam", reading, ctx,
        ((scene*)ctx)->cameras, false, def.cam);
    serialize_namedref_attr(val.ist, js, "ist", reading, ctx,
        ((scene*)ctx)->instances, false, def.ist);
    serialize_namedref_attr(val.env, js, "env", reading, ctx,
        ((scene*)ctx)->environments, false, def.env);
}

// Serialize enum
void serialize(animation_type& val, json& js, bool reading, void* ctx) {
    static auto names = std::map<animation_type, std::string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    serialize(val, js, reading, ctx, names);
}

// Serialize struct
void serialize(animation& val, json& js, bool reading, void* ctx) {
    static auto def = animation();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(val.path, js, "path", reading, ctx, false, def.path);
    serialize_attr(val.group, js, "group", reading, ctx, false, def.group);
    serialize_attr(val.type, js, "type", reading, ctx, false, def.type);
    serialize_namedref_attr(val.targets, js, "targets", reading, ctx,
        ((scene*)ctx)->nodes, false, def.targets);
    if (val.path == "") {
        serialize_attr(val.times, js, "times", reading, ctx, false, def.times);
        serialize_attr(val.translation, js, "translation", reading, ctx, false,
            def.translation);
        serialize_attr(
            val.rotation, js, "rotation", reading, ctx, false, def.rotation);
        serialize_attr(val.scale, js, "scale", reading, ctx, false, def.scale);
    }
}

// Serialize struct
void serialize(scene& val, json& js, bool reading, void* ctx) {
    static auto def = scene();
    serialize_obj(js, reading, ctx);
    serialize_attr(val.name, js, "name", reading, ctx, false, def.name);
    serialize_attr(
        val.cameras, js, "cameras", reading, ctx, false, def.cameras);
    serialize_attr(
        val.textures, js, "textures", reading, ctx, false, def.textures);
    serialize_attr(
        val.materials, js, "materials", reading, ctx, false, def.materials);
    serialize_attr(val.shapes, js, "shapes", reading, ctx, false, def.shapes);
    serialize_attr(
        val.subdivs, js, "subdivs", reading, ctx, false, def.subdivs);
    serialize_attr(
        val.instances, js, "instances", reading, ctx, false, def.instances);
    serialize_attr(val.environments, js, "environments", reading, ctx, false,
        def.environments);
    serialize_attr(val.nodes, js, "nodes", reading, ctx, false, def.nodes);
    serialize_attr(
        val.animations, js, "animations", reading, ctx, false, def.animations);
}

// Load a text file
std::string load_text(const std::string& filename) {
    // https://stackoverflow.com/questions/2912520/read-file-contents-into-a-string-in-c
    auto f = fopen(filename.c_str(), "r");
    if (!f) throw std::runtime_error("could not load " + filename);
    // Determine file size
    fseek(f, 0, SEEK_END);
    auto size = ftell(f);
    rewind(f);
    // read
    auto buf = std::vector<char>(size);
    if (fread(buf.data(), sizeof(char), size, f) != size)
        throw std::runtime_error("could not load " + filename);
    return std::string(buf.data(), buf.size());
}

// load json meshes
void load_json_meshes(
    const std::string& filename, scene* scn, bool skip_missing) {
    // load meshes
    auto dirname = path_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->path == "") continue;
        auto filename = dirname + shp->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            load_mesh(filename, shp->lines, shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Save json meshes
void save_json_meshes(
    const std::string& filename, const scene* scn, bool skip_missing) {
    // load meshes
    auto dirname = path_dirname(filename);
    for (auto shp : scn->shapes) {
        if (shp->path == "") continue;
        auto filename = dirname + shp->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        try {
            save_mesh(filename, shp->lines, shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        } catch (std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// load json textutes
void load_json_textures(
    const std::string& filename, scene* scn, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // load images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        txt->img =
            load_image(filename, txt->width, txt->height, ldr_gamma.at(txt));
        if (txt->img.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}

// Save json texture images.
void save_json_textures(
    const std::string& filename, const scene* scn, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // save images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        if (txt->img.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
        if (!txt->img.empty()) {
            ok = save_image(
                filename, txt->width, txt->height, txt->img, ldr_gamma.at(txt));
        }
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
}

// Save a text file
void save_text(const std::string& filename, const std::string& str) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    auto num = fwrite(str.c_str(), 1, str.size(), f);
    if (num != str.size())
        throw std::runtime_error("cannot write file " + filename);
    fclose(f);
}

// Load a JSON object
json load_json(const std::string& filename) {
    auto txt = load_text(filename);
    return json::parse(txt.begin(), txt.end());
}

// Save a JSON object
void save_json(const std::string& filename, const json& js) {
    save_text(filename, js.dump(2));
}

// Load a scene in the builtin JSON format.
scene* load_json_scene(const std::string& filename, bool skip_missing) {
    // load json
    auto js = load_json(filename);
    auto scn = new scene();
    serialize(scn, js, true, scn);
    // load binaries
    load_json_meshes(filename, scn, skip_missing);
    load_json_textures(filename, scn, skip_missing);
    return scn;
}

// Save a scene in the builtin JSON format.
void save_json_scene(
    const std::string& filename, const scene* scn, bool skip_missing) {
    // save json
    auto js = json();
    if (scn) serialize((scene&)(*scn), js, false, (scene*)scn);
    save_json(filename, js);
    // save meshes and textures
    save_json_meshes(filename, scn, skip_missing);
    save_json_textures(filename, scn, skip_missing);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OBJ CONVESION
// -----------------------------------------------------------------------------
namespace ygl {

// skip whitespace
inline void parse_skipws(char*& s) {
    while (*s == ' ' || *s == '\t' || *s == '\r' || *s == '\n') s++;
}
// check if whitespace
inline bool parse_isws(char*& s) {
    return *s == ' ' || *s == '\t' || *s == '\r' || *s == '\n';
}
// check if whitespace
inline void parse_remove_comment(char*& s, char cmd_char) {
    auto ss = s;
    while (*ss) {
        if (*ss == cmd_char) {
            *ss = 0;
            break;
        }
        ss++;
    }
}

#if YGL_FASTPARSE
// parse base value
inline int parse_int_base(char*& s) {
    parse_skipws(s);
    auto val = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
    return val;
}

// parse base value
inline float parse_float_base(char*& s) {
    parse_skipws(s);
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
    auto val = (float)dval;
    return val;
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
inline std::string parse_string_base(char*& s) {
    parse_skipws(s);
    char buf[4096];
    auto val = buf;
    while (*s && !parse_isws(s)) *val++ = *s++;
    *val = 0;
    return buf;
}

#else

// parse base value
inline int parse_int_base(char*& s) {
    auto len = 0;
    auto val = 0;
    sscanf(s, "%d%n", &val, &len);
    s += len;
    return val;
}

// parse base value
inline float parse_float_base(char*& s) {
    auto len = 0;
    auto val = 0.0f;
    sscanf(s, "%f%n", &val, &len);
    s += len;
    return val;
}

// parse base value
inline std::string parse_string_base(char*& s) {
    char buf[4096];
    auto val = buf;
    auto len = 0;
    sscanf(s, "%s%n", buf, &len);
    if (len) {
        s += len;
        val = buf;
    } else {
        val = "";
    }
    return val;
}

#endif

// parse value
inline int parse_int(char*& s) { return parse_int_base(s); }
inline float parse_float(char*& s) { return parse_float_base(s); }
inline std::string parse_string(char*& s) { return parse_string_base(s); }
inline byte parse_uchar(char*& s) { return (byte)parse_int_base(s); }
inline std::vector<std::string> parse_strings(char*& s) {
    auto vals = std::vector<std::string>();
    parse_skipws(s);
    while (*s) {
        auto val = parse_string(s);
        if (val == "") break;
        vals.push_back(val);
        parse_skipws(s);
    }
    return vals;
}
inline float parse_bool(char*& s) { return parse_int(s); }
inline vec2i parse_vec2i(char*& s) {
    auto x = parse_int(s);
    auto y = parse_int(s);
    return {x, y};
}
inline vec2f parse_vec2f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    return {x, y};
}
inline vec3f parse_vec3f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    auto z = parse_float(s);
    return {x, y, z};
}
inline vec4f parse_vec4f(char*& s) {
    auto x = parse_float(s);
    auto y = parse_float(s);
    auto z = parse_float(s);
    auto w = parse_float(s);
    return {x, y, z, w};
}
inline frame3f parse_frame3f(char*& s) {
    auto val = identity_frame3f;
    for (auto i = 0; i < 12; i++) (&val.x.x)[i] = parse_float(s);
    return val;
}
inline bool parse_getline(FILE* f, char* buf, int max_length) {
    if (fgets(buf, max_length, f)) {
        auto ss = buf;
        while (*ss) {
            if (*ss == '\r' || *ss == '\n') {
                *ss = 0;
                break;
            }
            ss++;
        }
        return true;
    } else {
        return false;
    }
}

inline vec3i parse_objvert(char*& s, const vec3i& vert_size) {
    auto token = parse_string(s);
    auto val = vec3i{-1, -1, -1};
    auto i = 0;
    auto sb = (char*)token.c_str();
    while (i < 3 && *sb) {
        (&val.x)[i] = parse_int(sb);
        (&val.x)[i] = ((&val.x)[i] < 0) ? (&vert_size.x)[i] + (&val.x)[i] :
                                          (&val.x)[i] - 1;
        if (*sb != '/') break;
        while (*sb == '/') {
            sb++;
            i++;
        }
    }
    return val;
}

template <typename T>
inline T read_value(FILE* fs) {
    T val;
    if (fread(&val, sizeof(val), 1, fs) != 1)
        throw std::runtime_error("error reading value");
    return val;
}

inline int read_int(FILE* fs) { return read_value<int>(fs); }
inline float read_float(FILE* fs) { return read_value<float>(fs); }
inline byte read_uchar(FILE* fs) { return read_value<byte>(fs); }

// Load MTL
void load_obj_materials(const std::string& filename, scene* scn,
    std::unordered_map<std::string, material*>& mmap,
    std::unordered_map<std::string, texture*>& tmap, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // Parse texture options and name
    auto add_texture = [scn, &tmap](char*& s) {
        // get tokens
        auto tokens = parse_strings(s);
        if (tokens.empty()) return (texture*)nullptr;

        // texture name
        auto path = tokens.back();
        for (auto& c : path)
            if (c == '\\') c = '/';
        if (tmap.find(path) != tmap.end()) { return tmap.at(path); }

        // create texture
        auto txt = new texture();
        txt->name = path;
        txt->path = path;
        scn->textures.push_back(txt);
        tmap[path] = txt;

        // texture options
        for (auto i = 0; i < tokens.size(); i++) {
            if (tokens[i] == "-clamp") txt->clamp = true;
            // TODO: bump scale
        }

        return txt;
    };

    // add a material preemptively to avoid crashes
    scn->materials.push_back(new material());
    auto mat = scn->materials.back();

    // read the file line by line
    char line[4096];
    while (parse_getline(fs, line, sizeof(line))) {
        // prepare to parse
        auto ss = line;
        parse_skipws(ss);
        parse_remove_comment(ss, '#');

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "newmtl") {
            mat = new material();
            mat->name = parse_string(ss);
            scn->materials.push_back(mat);
            mmap[mat->name] = mat;
        } else if (cmd == "illum") {
            // auto illum = parse_int(ss);
            // TODO: something with illum
        } else if (cmd == "Ke") {
            mat->ke = parse_vec3f(ss);
        } else if (cmd == "Kd") {
            mat->kd = parse_vec3f(ss);
        } else if (cmd == "Ks") {
            mat->ks = parse_vec3f(ss);
        } else if (cmd == "Kt") {
            mat->kt = parse_vec3f(ss);
        } else if (cmd == "Tf") {
            auto nchan = 0;
            while (*ss && nchan < 3) {
                (&mat->kt.x)[nchan++] = parse_float(ss);
                parse_skipws(ss);
            }
            if (nchan < 3) mat->kt = {mat->kt.x, mat->kt.x, mat->kt.x};
            if (flip_tr) mat->kt = vec3f{1, 1, 1} - mat->kt;
        } else if (cmd == "Tr") {
            auto nchan = 0;
            auto tr = zero3f;
            while (*ss && nchan < 3) {
                (&tr.x)[nchan++] = parse_float(ss);
                parse_skipws(ss);
            }
            if (nchan < 3) tr = {tr.x, tr.x, tr.x};
            mat->op = (tr.x + tr.y + tr.z) / 3;
            if (flip_tr) mat->op = 1 - mat->op;
        } else if (cmd == "Ns") {
            auto ns = parse_float(ss);
            mat->rs = pow(2 / (ns + 2), 1 / 4.0f);
            if (mat->rs < 0.01f) mat->rs = 0;
            if (mat->rs > 0.99f) mat->rs = 1;
        } else if (cmd == "d") {
            mat->op = parse_float(ss);
        } else if (cmd == "Pr" || cmd == "rs") {
            mat->rs = parse_float(ss);
        } else if (cmd == "map_Ke") {
            mat->ke_txt = add_texture(ss);
        } else if (cmd == "map_Kd") {
            mat->kd_txt = add_texture(ss);
        } else if (cmd == "map_Ks") {
            mat->ks_txt = add_texture(ss);
        } else if (cmd == "map_Tr") {
            mat->kt_txt = add_texture(ss);
        } else if (cmd == "map_d" || cmd == "map_Tr") {
            mat->op_txt = add_texture(ss);
        } else if (cmd == "map_Pr" || cmd == "map_rs") {
            mat->rs_txt = add_texture(ss);
        } else if (cmd == "map_occ" || cmd == "occ") {
            mat->occ_txt = add_texture(ss);
        } else if (cmd == "map_bump" || cmd == "bump") {
            mat->bump_txt = add_texture(ss);
        } else if (cmd == "map_disp" || cmd == "disp") {
            mat->disp_txt = add_texture(ss);
        } else if (cmd == "map_norm" || cmd == "norm") {
            mat->norm_txt = add_texture(ss);
        }
    }

    // remove first fake material
    scn->materials.erase(scn->materials.begin());

    // clone
    fclose(fs);
}

// Loads an OBJ
void load_obj_geometry(const std::string& filename, scene* scn,
    bool split_shapes, bool flip_texcoord, bool flip_tr) {
    // splitting policy
    auto split_material = split_shapes;
    auto split_group = split_shapes;
    auto split_smoothing = split_shapes;

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // current parsing values
    auto matname = std::string();
    auto oname = std::string();
    auto gname = std::string();
    auto smoothing = true;
    auto ist = (instance*)nullptr;

    // vertices
    auto pos = std::deque<vec3f>();
    auto norm = std::deque<vec3f>();
    auto texcoord = std::deque<vec2f>();

    // object maps
    auto tmap = std::unordered_map<std::string, texture*>();
    auto mmap = std::unordered_map<std::string, material*>();

    // vertex maps
    auto name_map = std::unordered_map<std::string, int>();
    auto vert_map = std::unordered_map<vec3i, int>();
    auto pos_map = std::unordered_map<int, int>();
    auto norm_map = std::unordered_map<int, int>();
    auto texcoord_map = std::unordered_map<int, int>();

    // add object if needed
    auto is_instance_empty = [](instance* ist) {
        if (ist->sbd) {
            return ist->sbd->pos.empty();
        } else if (ist->shp) {
            return ist->shp->pos.empty();
        } else {
            return true;
        }
    };
    auto add_instance = [&](scene* scn, const std::string& objname,
                            const std::string& matname,
                            const std::string& groupname, bool smoothing) {
        if (scn->instances.empty() || objname != scn->instances.back()->name ||
            !is_instance_empty(scn->instances.back())) {
            auto ist = new instance();
            scn->instances.push_back(ist);
            ist->shp = new shape();
            scn->shapes.push_back(ist->shp);
        }
        name_map[objname] += 1;
        auto name = (name_map[objname] == 1) ?
                        objname :
                        objname + "_" + std::to_string(name_map[objname] - 1);
        if (objname == "") name = "object" + name;
        auto ist = scn->instances.back();
        ist->name = name;
        if (ist->shp) ist->shp->name = ist->name;
        if (ist->sbd) ist->sbd->name = ist->name;
        if (matname != "") ist->mat = mmap.at(matname);
        vert_map.clear();
        pos_map.clear();
        norm_map.clear();
        texcoord_map.clear();
        return ist;
    };

    // current objet
    ist = add_instance(scn, "", "", "", true);

    // read the file line by line
    char line[4096];
    while (parse_getline(fs, line, sizeof(line))) {
        // prepare to parse
        auto ss = line;
        parse_skipws(ss);
        parse_remove_comment(ss, '#');

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            pos.push_back(parse_vec3f(ss));
        } else if (cmd == "vn") {
            norm.push_back(parse_vec3f(ss));
        } else if (cmd == "vt") {
            texcoord.push_back(parse_vec2f(ss));
            if (flip_texcoord) texcoord.back().y = 1 - texcoord.back().y;
        } else if (cmd == "f" || cmd == "l") {
            auto num = 0;
            vec3i verts[128];
            int vids[128];
            auto vert_size =
                vec3i{(int)pos.size(), (int)texcoord.size(), (int)norm.size()};
            // elem.material = (int)oobj->materials.size() - 1;
            parse_skipws(ss);
            while (*ss) {
                verts[num++] = parse_objvert(ss, vert_size);
                parse_skipws(ss);
            }
            if (ist->sbd) {
                // TODO: subdivs
            } else if (ist->shp) {
                for (auto i = 0; i < num; i++) {
                    auto it = vert_map.find(verts[i]);
                    if (it == vert_map.end()) {
                        auto nverts = (int)ist->shp->pos.size();
                        vert_map.insert(it, {verts[i], nverts});
                        vids[i] = nverts;
                        if (verts[i].x >= 0)
                            ist->shp->pos.push_back(pos.at(verts[i].x));
                        if (verts[i].y >= 0)
                            ist->shp->texcoord.push_back(
                                texcoord.at(verts[i].y));
                        if (verts[i].z >= 0)
                            ist->shp->norm.push_back(norm.at(verts[i].z));
                    } else {
                        vids[i] = it->second;
                    }
                }
                if (cmd == "f") {
                    for (auto i = 2; i < num; i++)
                        ist->shp->triangles.push_back(
                            {vids[0], vids[i - 1], vids[i]});
                }
                if (cmd == "l") {
                    for (auto i = 1; i < num; i++)
                        ist->shp->lines.push_back({vids[i - 1], vids[i]});
                }
            }
        } else if (cmd == "p") {
            // TODO: support points
        } else if (cmd == "o") {
            oname = parse_string(ss);
            gname = "";
            matname = "";
            smoothing = true;
            ist = add_instance(scn, oname, matname, gname, smoothing);
        } else if (cmd == "usemtl") {
            matname = parse_string(ss);
            if (split_material) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            } else {
                if (matname != "") ist->mat = mmap.at(matname);
            }
        } else if (cmd == "g") {
            gname = parse_string(ss);
            if (split_group) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            }
        } else if (cmd == "s") {
            auto name = parse_string(ss);
            smoothing = (name == "on");
            if (split_smoothing) {
                ist = add_instance(scn, oname, matname, gname, smoothing);
            }
        } else if (cmd == "mtllib") {
            auto mtlname = parse_string(ss);
            auto mtlpath = path_dirname(filename) + mtlname;
            load_obj_materials(mtlpath, scn, mmap, tmap, flip_tr);
        } else if (cmd == "c") {
            auto cam = new camera();
            cam->name = parse_string(ss);
            cam->ortho = parse_bool(ss);
            cam->width = parse_float(ss);
            cam->height = parse_float(ss);
            cam->focal = parse_float(ss);
            cam->focus = parse_float(ss);
            cam->aperture = parse_float(ss);
            cam->frame = parse_frame3f(ss);
            scn->cameras.push_back(cam);
        } else if (cmd == "e") {
            auto env = new environment();
            env->name = parse_string(ss);
            env->ke = parse_vec3f(ss);
            auto ke_txt = parse_string(ss);
            if (ke_txt != "\"\"") {
                if (tmap.find(ke_txt) == tmap.end()) {
                    auto txt = new texture();
                    txt->name = ke_txt;
                    txt->path = ke_txt;
                    tmap[ke_txt] = txt;
                    scn->textures.push_back(txt);
                }
                env->ke_txt = tmap.at(ke_txt);
            }
            env->frame = parse_frame3f(ss);
            scn->environments.push_back(env);
        } else {
            // unused
        }
    }

    // cleanup empty
    for (auto idx = 0; idx < scn->instances.size(); idx++) {
        if (!is_instance_empty(scn->instances[idx])) continue;
        auto ist = scn->instances[idx];
        if (ist->shp) {
            scn->shapes.erase(
                std::find(scn->shapes.begin(), scn->shapes.end(), ist->shp));
            delete ist->shp;
        }
        if (ist->sbd) {
            scn->subdivs.erase(
                std::find(scn->subdivs.begin(), scn->subdivs.end(), ist->sbd));
            delete ist->sbd;
        }
        delete ist;
        scn->instances.erase(scn->instances.begin() + idx);
        idx--;
    }

    // updates
    update_transforms(scn);
    update_bbox(scn);

    // close file
    fclose(fs);
}

// load json textutes
void load_obj_textures(
    const std::string& filename, scene* scn, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // load images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        txt->img =
            load_image(filename, txt->width, txt->height, ldr_gamma.at(txt));
        if (txt->img.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}

// Load a scene
scene* load_obj_scene(const std::string& filename, bool skip_missing) {
    auto ext = path_extension(filename);
    auto scn = new scene();
    load_obj_geometry(filename, scn, true, true, true);
    load_obj_textures(filename, scn, skip_missing);
    scn->name = path_filename(filename);
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
    return scn;
}

void save_obj_geometry(const std::string& filename, const scene* scn) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot save file " + filename);

    // material library
    if (!scn->materials.empty()) {
        auto mtlname = path_basename(filename) + ".mtl";
        fprintf(f, "mtllib %s\n", mtlname.c_str());
    }

    // cameras
    for (auto cam : scn->cameras) {
        fprintf(f, "c");
        fprintf(f, " %s", cam->name.c_str());
        fprintf(f, " %d", (int)cam->ortho);
        fprintf(f, " %g", cam->width);
        fprintf(f, " %g", cam->height);
        fprintf(f, " %g", cam->focal);
        fprintf(f, " %g", cam->focus);
        fprintf(f, " %g", cam->aperture);
        for (auto i = 0; i < 12; i++) fprintf(f, " %g", (&cam->frame.x.x)[i]);
        fprintf(f, "\n");
    }

    // environments
    for (auto env : scn->environments) {
        fprintf(f, "e");
        fprintf(f, " %s", env->name.c_str());
        fprintf(f, " %g %g %g", env->ke.x, env->ke.y, env->ke.z);
        fprintf(f, " %s", (env->ke_txt) ? env->ke_txt->path.c_str() : "\"\"");
        for (auto i = 0; i < 12; i++) fprintf(f, " %g", (&env->frame.x.x)[i]);
        fprintf(f, "\n");
    }

    // shapes
    auto offset = zero3i;
    for (auto ist : scn->instances) {
        fprintf(f, "o %s\n", ist->name.c_str());
        if (ist->mat) fprintf(f, "usemtl %s\n", ist->mat->name.c_str());
        if (!ist->sbd) {
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->shp->pos)
                    fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
                for (auto& n : ist->shp->norm)
                    fprintf(f, "vn %g %g %g\n", n.x, n.y, n.z);
                for (auto& t : ist->shp->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            } else {
                for (auto& p : ist->shp->pos) {
                    auto tp = transform_point(ist->frame, p);
                    fprintf(f, "v %g %g %g\n", tp.x, tp.y, tp.z);
                }
                for (auto& n : ist->shp->norm) {
                    auto tn = transform_direction(ist->frame, n);
                    fprintf(f, "vn %g %g %g\n", tn.x, tn.y, tn.z);
                }
                for (auto& t : ist->shp->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            }
            if (!ist->shp->texcoord.empty() && !ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n",
                        t.x + offset.x + 1, t.x + offset.y + 1,
                        t.x + offset.z + 1, t.y + offset.x + 1,
                        t.y + offset.y + 1, t.y + offset.z + 1,
                        t.z + offset.x + 1, t.z + offset.y + 1,
                        t.z + offset.z + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d/%d/%d %d/%d/%d\n", l.x + offset.x + 1,
                        l.x + offset.y + 1, l.x + offset.z + 1,
                        l.y + offset.x + 1, l.y + offset.y + 1,
                        l.y + offset.z + 1);
            } else if (!ist->shp->texcoord.empty() && ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d/%d %d/%d %d/%d\n", t.x + offset.x + 1,
                        t.x + offset.y + 1, t.y + offset.x + 1,
                        t.y + offset.y + 1, t.z + offset.x + 1,
                        t.z + offset.y + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d/%d %d/%d\n", l.x + offset.x + 1,
                        l.x + offset.y + 1, l.y + offset.x + 1,
                        l.y + offset.y + 1);
            } else if (ist->shp->texcoord.empty() && !ist->shp->norm.empty()) {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d//%d %d//%d %d//%d\n", t.x + offset.x + 1,
                        t.x + offset.z + 1, t.y + offset.x + 1,
                        t.y + offset.z + 1, t.z + offset.x + 1,
                        t.z + offset.z + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(f, "l %d//%d %d//%d\n", l.x + offset.x + 1,
                        l.x + offset.z + 1, l.y + offset.x + 1,
                        l.y + offset.z + 1);
            } else {
                for (auto& t : ist->shp->triangles)
                    fprintf(f, "f %d %d %d\n", t.x + offset.x + 1,
                        t.y + offset.x + 1, t.z + offset.x + 1);
                for (auto& l : ist->shp->lines)
                    fprintf(
                        f, "l %d %d\n", l.x + offset.x + 1, l.y + offset.x + 1);
            }
            offset.x += ist->shp->pos.size();
            offset.y += ist->shp->texcoord.size();
            offset.z += ist->shp->norm.size();
        } else {
            if (ist->frame == identity_frame3f) {
                for (auto& p : ist->sbd->pos)
                    fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
                for (auto& t : ist->sbd->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            } else {
                for (auto& p : ist->sbd->pos) {
                    auto tp = transform_point(ist->frame, p);
                    fprintf(f, "v %g %g %g\n", tp.x, tp.y, tp.z);
                }
                for (auto& t : ist->sbd->texcoord)
                    fprintf(f, "vt %g %g\n", t.x, t.y);
            }
            if (!ist->sbd->texcoord.empty()) {
                for (auto i = 0; i < ist->sbd->quads_pos.size(); i++) {
                    auto qp = ist->sbd->quads_pos[i];
                    auto qt = ist->sbd->quads_texcoord[i];
                    if (qp.z == qp.w) {
                        fprintf(f, "f %d/%d %d/%d %d/%d\n", qp.x + offset.x + 1,
                            qt.x + offset.y + 1, qp.y + offset.x + 1,
                            qt.y + offset.y + 1, qp.z + offset.x + 1,
                            qt.z + offset.y + 1);
                    } else {
                        fprintf(f, "f %d/%d %d/%d %d/%d %d/%d\n",
                            qp.x + offset.x + 1, qt.x + offset.y + 1,
                            qp.y + offset.x + 1, qt.y + offset.y + 1,
                            qp.z + offset.x + 1, qt.z + offset.y + 1,
                            qp.w + offset.x + 1, qt.w + offset.y + 1);
                    }
                }
            } else {
                for (auto& q : ist->sbd->quads_pos) {
                    if (q.z == q.w) {
                        fprintf(f, "f %d %d %d\n", q.x + offset.x + 1,
                            q.y + offset.x + 1, q.z + offset.x + 1);
                    } else {
                        fprintf(f, "f %d %d %d %d\n", q.x + offset.x + 1,
                            q.y + offset.x + 1, q.z + offset.x + 1,
                            q.w + offset.x + 1);
                    }
                }
            }
            offset.x += ist->sbd->pos.size();
            offset.y += ist->sbd->texcoord.size();
        }
    }

    fclose(f);
}

void save_obj_materials(const std::string& filename, const scene* scn) {
    if (scn->materials.empty()) return;

    auto mtlname = replace_path_extension(filename, ".mtl");
    auto f = fopen(mtlname.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot open filename " + mtlname);

    // for each material, dump all the values
    for (auto mat : scn->materials) {
        fprintf(f, "newmtl %s\n", mat->name.c_str());
        fprintf(f, "  illum %d\n", 2);
        if (mat->ke != zero3f)
            fprintf(f, "  Ke %g %g %g\n", mat->ke.x, mat->ke.y, mat->ke.z);
        if (mat->kd != zero3f)
            fprintf(f, "  Kd %g %g %g\n", mat->kd.x, mat->kd.y, mat->kd.z);
        if (mat->ks != zero3f)
            fprintf(f, "  Ks %g %g %g\n", mat->ks.x, mat->ks.y, mat->ks.z);
        if (mat->kt != zero3f)
            fprintf(f, "  Kt %g %g %g\n", mat->kt.x, mat->kt.y, mat->kt.z);
        if (mat->rs != 1.0f)
            fprintf(f, "  Ns %d\n",
                (int)clamp(2 / pow(mat->rs + 1e-10f, 4.0f) - 2, 0.0f, 1.0e12f));
        if (mat->op != 1.0f) fprintf(f, "  d %g\n", mat->op);
        if (mat->rs != -1.0f) fprintf(f, "  Pr %g\n", mat->rs);
        if (mat->ke_txt) fprintf(f, "  map_Ke %s\n", mat->ke_txt->path.c_str());
        if (mat->kd_txt) fprintf(f, "  map_Kd %s\n", mat->kd_txt->path.c_str());
        if (mat->ks_txt) fprintf(f, "  map_Ks %s\n", mat->ks_txt->path.c_str());
        if (mat->kt_txt) fprintf(f, "  map_Kt %s\n", mat->kt_txt->path.c_str());
        if (mat->op_txt) fprintf(f, "  map_d %s\n", mat->op_txt->path.c_str());
        if (mat->rs_txt) fprintf(f, "  map_Pr %s\n", mat->rs_txt->path.c_str());
        if (mat->occ_txt)
            fprintf(f, "  map_occ %s\n", mat->occ_txt->path.c_str());
        if (mat->bump_txt)
            fprintf(f, "  map_bump %s\n", mat->bump_txt->path.c_str());
        if (mat->disp_txt)
            fprintf(f, "  map_disp %s\n", mat->disp_txt->path.c_str());
        if (mat->norm_txt)
            fprintf(f, "  map_norm %s\n", mat->norm_txt->path.c_str());
        fprintf(f, "\n");
    }

    fclose(f);
}

void save_obj_textures(
    const std::string& filename, const scene* scn, bool skip_missing) {
    // set gamma
    auto ldr_gamma = std::unordered_map<texture*, float>{{nullptr, 1.0f}};
    for (auto txt : scn->textures) ldr_gamma[txt] = 2.2f;
    for (auto mat : scn->materials) {
        ldr_gamma[mat->ke_txt] = 2.2f;
        ldr_gamma[mat->kd_txt] = 2.2f;
        ldr_gamma[mat->ks_txt] = 2.2f;
        ldr_gamma[mat->kt_txt] = 2.2f;
        ldr_gamma[mat->rs_txt] = 1;
        ldr_gamma[mat->op_txt] = 1;
        ldr_gamma[mat->norm_txt] = 1;
        ldr_gamma[mat->disp_txt] = 1;
        ldr_gamma[mat->bump_txt] = 1;
    }
    for (auto env : scn->environments) { ldr_gamma[env->ke_txt] = 2.2f; }

    // save images
    auto dirname = path_dirname(filename);
    for (auto txt : scn->textures) {
        if (txt->img.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
        if (!txt->img.empty()) {
            ok = save_image(
                filename, txt->width, txt->height, txt->img, ldr_gamma.at(txt));
        }
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
}

// Save a scene
void save_obj_scene(
    const std::string& filename, const scene* scn, bool skip_missing) {
    save_obj_geometry(filename, scn);
    save_obj_materials(filename, scn);
    save_obj_textures(filename, scn, skip_missing);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GLTF CONVESION
// -----------------------------------------------------------------------------
namespace ygl {

static bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

// Flattens a gltf file into a flattened asset.
scene* gltf_to_scene(const glTF* gltf) {
    // clear asset
    auto scn = new scene();

    // convert textures
    auto tmap = std::map<glTFImage*, texture*>{};
    for (auto gimg : gltf->images) {
        auto txt = new texture();
        txt->name = gimg->name;
        txt->path = (startswith(gimg->uri, "data:")) ?
                        std::string("[glTF-inline].png") :
                        gimg->uri;
        txt->width = gimg->data.width;
        txt->height = gimg->data.height;
        txt->img = gimg->data.img;
        scn->textures.push_back(txt);
        tmap[gimg] = txt;
    }

    // add a texture
    auto add_texture = [&tmap, gltf](const glTFTextureInfo* ginfo,
                           bool normal = false, bool occlusion = false) {
        if (!ginfo) return (texture*)nullptr;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return (texture*)nullptr;
        auto gimg = gltf->get(gtxt->source);
        if (!gimg) return (texture*)nullptr;
        auto txt = tmap.at(gimg);
        auto gsmp = gltf->get(gtxt->sampler);
        if (gsmp) {
            txt->clamp = gsmp->wrapS == glTFSamplerWrapS::ClampToEdge ||
                         gsmp->wrapT == glTFSamplerWrapT::ClampToEdge;
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            txt->scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            txt->scale = ninfo->strength;
        }
        return txt;
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new material();
        mat->name = gmat->name;
        mat->ke = gmat->emissiveFactor;
        mat->ke_txt = add_texture(gmat->emissiveTexture);
        if (gmat->pbrSpecularGlossiness) {
            mat->base_metallic = false;
            mat->gltf_textures = true;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor.x, gsg->diffuseFactor.y,
                gsg->diffuseFactor.z};
            mat->op = gsg->diffuseFactor.w;
            mat->ks = gsg->specularFactor;
            mat->rs = 1 - gsg->glossinessFactor;
            mat->kd_txt = add_texture(gsg->diffuseTexture);
            mat->ks_txt = add_texture(gsg->specularGlossinessTexture);
            mat->rs_txt = mat->ks_txt;
        } else if (gmat->pbrMetallicRoughness) {
            mat->base_metallic = true;
            mat->gltf_textures = true;
            auto gmr = gmat->pbrMetallicRoughness;
            mat->ke = gmat->emissiveFactor;
            mat->kd = {gmr->baseColorFactor.x, gmr->baseColorFactor.y,
                gmr->baseColorFactor.z};
            mat->op = gmr->baseColorFactor.w;
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = add_texture(gmr->baseColorTexture, true);
            mat->ks_txt = add_texture(gmr->metallicRoughnessTexture, true);
            mat->rs_txt = mat->ks_txt;
        }
        mat->occ_txt = add_texture(gmat->occlusionTexture, false, true);
        mat->norm_txt = add_texture(gmat->normalTexture, true, false);
        mat->double_sided = gmat->doubleSided;
        scn->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<std::pair<shape*, material*>>>();
    for (auto gmesh : gltf->meshes) {
        meshes.push_back({});
        auto sid = 0;
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
            shp->name =
                gmesh->name + ((sid) ? std::to_string(sid) : std::string());
            sid++;
            for (auto gattr : gprim->attributes) {
                auto semantic = gattr.first;
                auto vals = accessor_view(gltf, gltf->get(gattr.second));
                if (semantic == "POSITION") {
                    shp->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->pos.push_back(vals.getv3f(i));
                } else if (semantic == "NORMAL") {
                    shp->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->norm.push_back(vals.getv3f(i));
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    shp->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->texcoord.push_back(vals.getv2f(i));
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    shp->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->color.push_back(vals.getv4f(i, {0, 0, 0, 1}));
                } else if (semantic == "TANGENT") {
                    shp->tangsp.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->tangsp.push_back(vals.getv4f(i));
                    for (auto& t : shp->tangsp) t.w = -t.w;
                } else if (semantic == "RADIUS") {
                    shp->radius.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->radius.push_back(vals.get(i, 0));
                } else {
                    // ignore
                }
            }
            // indices
            if (!gprim->indices) {
                if (gprim->mode == glTFMeshPrimitiveMode::Triangles) {
                    shp->triangles.reserve(shp->pos.size() / 3);
                    for (auto i = 0; i < shp->pos.size() / 3; i++)
                        shp->triangles.push_back(
                            {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                } else if (gprim->mode == glTFMeshPrimitiveMode::TriangleFan) {
                    shp->triangles.reserve(shp->pos.size() - 2);
                    for (auto i = 2; i < shp->pos.size(); i++)
                        shp->triangles.push_back({0, i - 1, i});
                } else if (gprim->mode ==
                           glTFMeshPrimitiveMode::TriangleStrip) {
                    shp->triangles.reserve(shp->pos.size() - 2);
                    for (auto i = 2; i < shp->pos.size(); i++)
                        shp->triangles.push_back({i - 2, i - 1, i});
                } else if (gprim->mode == glTFMeshPrimitiveMode::Lines) {
                    shp->lines.reserve(shp->pos.size() / 2);
                    for (auto i = 0; i < shp->pos.size() / 2; i++)
                        shp->lines.push_back({i * 2 + 0, i * 2 + 1});
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineLoop) {
                    shp->lines.reserve(shp->pos.size());
                    for (auto i = 1; i < shp->pos.size(); i++)
                        shp->lines.push_back({i - 1, i});
                    shp->lines.back() = {(int)shp->pos.size() - 1, 0};
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineStrip) {
                    shp->lines.reserve(shp->pos.size() - 1);
                    for (auto i = 1; i < shp->pos.size(); i++)
                        shp->lines.push_back({i - 1, i});
                } else if (gprim->mode == glTFMeshPrimitiveMode::NotSet ||
                           gprim->mode == glTFMeshPrimitiveMode::Points) {
                    log_warning("points not supported");
                } else {
                    throw std::runtime_error("unknown primitive type");
                }
            } else {
                auto indices = accessor_view(gltf, gltf->get(gprim->indices));
                if (gprim->mode == glTFMeshPrimitiveMode::Triangles) {
                    shp->triangles.reserve(indices.size());
                    for (auto i = 0; i < indices.size() / 3; i++)
                        shp->triangles.push_back({indices.geti(i * 3 + 0),
                            indices.geti(i * 3 + 1), indices.geti(i * 3 + 2)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::TriangleFan) {
                    shp->triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shp->triangles.push_back({indices.geti(0),
                            indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode ==
                           glTFMeshPrimitiveMode::TriangleStrip) {
                    shp->triangles.reserve(indices.size() - 2);
                    for (auto i = 2; i < indices.size(); i++)
                        shp->triangles.push_back({indices.geti(i - 2),
                            indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::Lines) {
                    shp->lines.reserve(indices.size() / 2);
                    for (auto i = 0; i < indices.size() / 2; i++)
                        shp->lines.push_back(
                            {indices.geti(i * 2 + 0), indices.geti(i * 2 + 1)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineLoop) {
                    shp->lines.reserve(indices.size());
                    for (auto i = 1; i < indices.size(); i++)
                        shp->lines.push_back(
                            {indices.geti(i - 1), indices.geti(i)});
                    shp->lines.back() = {
                        indices.geti(indices.size() - 1), indices.geti(0)};
                } else if (gprim->mode == glTFMeshPrimitiveMode::LineStrip) {
                    shp->lines.reserve(indices.size() - 1);
                    for (auto i = 1; i < indices.size(); i++)
                        shp->lines.push_back(
                            {indices.geti(i - 1), indices.geti(i)});
                } else if (gprim->mode == glTFMeshPrimitiveMode::NotSet ||
                           gprim->mode == glTFMeshPrimitiveMode::Points) {
                    log_warning("points not supported");
                } else {
                    throw std::runtime_error("unknown primitive type");
                }
            }
            meshes.back().push_back(
                {shp, scn->materials[(int)gprim->material]});
            scn->shapes.push_back(shp);
        }
    }

    // convert cameras
    for (auto gcam : gltf->cameras) {
        auto cam = new camera();
        cam->name = gcam->name;
        cam->ortho = gcam->type == glTFCameraType::Orthographic;
        if (cam->ortho) {
            log_warning("orthographic not supported well");
            auto ortho = gcam->orthographic;
            set_camera_fovy(cam, ortho->ymag, ortho->xmag / ortho->ymag);
            cam->near = ortho->znear;
            cam->far = ortho->zfar;
            cam->focus = flt_max;
            cam->aperture = 0;
        } else {
            auto persp = gcam->perspective;
            set_camera_fovy(cam, persp->yfov, persp->aspectRatio);
            cam->near = persp->znear;
            cam->far = (persp->zfar) ? persp->zfar : flt_max;
            cam->focus = flt_max;
            cam->aperture = 0;
        }
        scn->cameras.push_back(cam);
    }

    // convert nodes
    for (auto gnde : gltf->nodes) {
        auto nde = new node();
        nde->name = gnde->name;
        if (gnde->camera) nde->cam = scn->cameras[(int)gnde->camera];
        nde->translation = gnde->translation;
        nde->rotation = gnde->rotation;
        nde->scale = gnde->scale;
        nde->frame = mat_to_frame(gnde->matrix);
        scn->nodes.push_back(nde);
    }

    // set up parent pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        auto nde = scn->nodes[nid];
        for (auto cid : gnde->children) scn->nodes[(int)cid]->parent = nde;
    }

    // set up instances
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        if (!gnde->mesh) continue;
        auto nde = scn->nodes[nid];
        auto& shps = meshes.at((int)gnde->mesh);
        if (shps.empty()) continue;
        if (shps.size() == 1) {
            nde->ist = new instance();
            nde->ist->name = nde->name;
            nde->ist->shp = shps[0].first;
            nde->ist->mat = shps[0].second;
            scn->instances.push_back(nde->ist);
        } else {
            for (auto shp : shps) {
                auto child = new node();
                child->name = nde->name + "_" + shp.first->name;
                child->parent = nde;
                child->ist = new instance();
                child->ist->name = child->name;
                child->ist->shp = shp.first;
                child->ist->mat = shp.second;
                scn->instances.push_back(child->ist);
            }
        }
    }

    // keyframe type conversion
    static auto keyframe_types =
        std::unordered_map<glTFAnimationSamplerInterpolation, animation_type>{
            {glTFAnimationSamplerInterpolation::NotSet, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Linear, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Step, animation_type::step},
            {glTFAnimationSamplerInterpolation::CubicSpline,
                animation_type::bezier},
        };

    // convert animations
    for (auto ganm : gltf->animations) {
        auto aid = 0;
        auto sampler_map = std::unordered_map<vec2i, int>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto anm = new animation();
                anm->name = ((ganm->name != "") ? ganm->name : "anim") +
                            std::to_string(aid++);
                anm->group = ganm->name;
                auto input_view =
                    accessor_view(gltf, gltf->get(gsampler->input));
                anm->times.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    anm->times[i] = input_view.get(i);
                anm->type = keyframe_types.at(gsampler->interpolation);
                auto output_view =
                    accessor_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        anm->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->translation.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        anm->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->rotation.push_back(output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        anm->scale.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->scale.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Weights: {
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp = max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = std::vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            anm->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < anm->weights.size(); i++) {
                                anm->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    anm->weights[i][j] = values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        throw std::runtime_error("should not have gotten here");
                    }
                }
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = (int)scn->animations.size();
                scn->animations.push_back(anm);
            }
            scn->animations[sampler_map.at({(int)gchannel->sampler,
                                (int)gchannel->target->path})]
                ->targets.push_back(scn->nodes[(int)gchannel->target->node]);
        }
    }

    // compute transforms and bbox
    update_transforms(scn, 0);
    update_bbox(scn);

    // fix camera focus
    for (auto cam : scn->cameras) {
        auto center = (scn->bbox.min + scn->bbox.max) / 2;
        auto dist = dot(-cam->frame.z, center - cam->frame.o);
        if (dist > 0) cam->focus = dist;
    }

    return scn;
}

// Unflattnes gltf
glTF* scene_to_gltf(
    const scene* scn, const std::string& buffer_uri, bool separate_buffers) {
    auto gltf = new glTF();

    // add asset info
    gltf->asset = new glTFAsset();
    gltf->asset->generator = "Yocto/gltf";
    gltf->asset->version = "2.0";

    // convert cameras
    for (auto cam : scn->cameras) {
        auto gcam = new glTFCamera();
        gcam->name = cam->name;
        gcam->type = (cam->ortho) ? glTFCameraType::Orthographic :
                                    glTFCameraType::Perspective;
        if (cam->ortho) {
            auto ortho = new glTFCameraOrthographic();
            ortho->ymag = cam->width;
            ortho->xmag = cam->height;
            ortho->znear = cam->near;
            ortho->zfar = cam->far;
            gcam->orthographic = ortho;
        } else {
            auto persp = new glTFCameraPerspective();
            persp->yfov = eval_camera_fovy(cam);
            persp->aspectRatio = cam->width / cam->height;
            persp->znear = cam->near;
            persp->zfar = (cam->far >= flt_max) ? 0 : cam->far;
            gcam->perspective = persp;
        }
        gltf->cameras.push_back(gcam);
    }

    // convert textures
    auto tmap = std::map<texture*, glTFid<glTFImage>>{};
    for (auto txt : scn->textures) {
        auto gimg = new glTFImage();
        gimg->name = txt->name;
        gimg->uri = txt->path;
        gimg->data.width = txt->width;
        gimg->data.height = txt->height;
        gimg->data.img = txt->img;
        tmap[txt] = glTFid<glTFImage>((int)gltf->images.size());
        gltf->images.push_back(gimg);
    }

    // index of an object
    auto index = [](const auto& vec, auto& val) -> int {
        auto pos = find(vec.begin(), vec.end(), val);
        if (pos == vec.end()) return -1;
        return (int)(pos - vec.begin());
    };

    // add a texture and sampler
    auto add_texture = [gltf, &tmap](texture* txt, bool srgb, bool norm = false,
                           bool occ = false) {
        if (!txt) return (glTFTextureInfo*)nullptr;
        auto gtxt = new glTFTexture();
        gtxt->name = txt->name;
        gtxt->source = tmap.at(txt);

        // check if it is default
        if (txt->clamp) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = glTFSamplerWrapS::ClampToEdge;
            gsmp->wrapT = glTFSamplerWrapT::ClampToEdge;
            gsmp->minFilter = glTFSamplerMinFilter::LinearMipmapLinear;
            gsmp->magFilter = glTFSamplerMagFilter::Linear;
            gtxt->sampler = glTFid<glTFSampler>((int)gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        auto ginfo = (glTFTextureInfo*)nullptr;
        if (norm) {
            ginfo = new glTFMaterialNormalTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ((glTFMaterialNormalTextureInfo*)ginfo)->scale = txt->scale;
        } else if (occ) {
            ginfo = new glTFMaterialOcclusionTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ((glTFMaterialOcclusionTextureInfo*)ginfo)->strength = txt->scale;
        } else {
            ginfo = new glTFTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
        }
        return ginfo;
    };

    // convert materials
    for (auto mat : scn->materials) {
        auto gmat = new glTFMaterial();
        gmat->name = mat->name;
        gmat->emissiveFactor = mat->ke;
        gmat->emissiveTexture = add_texture(mat->ke_txt, true);
        if (!mat->base_metallic) {
            gmat->pbrSpecularGlossiness =
                new glTFMaterialPbrSpecularGlossiness();
            auto gsg = gmat->pbrSpecularGlossiness;
            gsg->diffuseFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
            gsg->specularFactor = mat->ks;
            gsg->glossinessFactor = 1 - mat->rs;
            gsg->diffuseTexture = add_texture(mat->kd_txt, true);
            gsg->specularGlossinessTexture = add_texture(mat->ks_txt, true);
        } else {
            gmat->pbrMetallicRoughness = new glTFMaterialPbrMetallicRoughness();
            auto gmr = gmat->pbrMetallicRoughness;
            gmr->baseColorFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
            gmr->metallicFactor = mat->ks.x;
            gmr->roughnessFactor = mat->rs;
            gmr->baseColorTexture = add_texture(mat->kd_txt, true);
            gmr->metallicRoughnessTexture = add_texture(mat->ks_txt, false);
        }
        gmat->normalTexture = (glTFMaterialNormalTextureInfo*)(add_texture(
            mat->norm_txt, false, true, false));
        gmat->occlusionTexture =
            (glTFMaterialOcclusionTextureInfo*)(add_texture(
                mat->occ_txt, false, false, true));
        gmat->doubleSided = mat->double_sided;
        gltf->materials.push_back(gmat);
    }

    // add buffer
    auto add_buffer = [&gltf](const std::string& buffer_uri) {
        auto gbuffer = new glTFBuffer();
        gltf->buffers.push_back(gbuffer);
        gbuffer->uri = buffer_uri;
        return gbuffer;
    };

    // init buffers
    auto gbuffer_global = add_buffer(buffer_uri);

    // add an optional buffer
    auto add_opt_buffer = [&gbuffer_global, buffer_uri, &add_buffer,
                              separate_buffers](const std::string& uri) {
        if (separate_buffers && uri != "") {
            return add_buffer(uri);
        } else {
            if (!gbuffer_global) gbuffer_global = add_buffer(buffer_uri);
            return gbuffer_global;
        }
    };

    // attribute handling
    auto add_accessor = [&gltf, &index](glTFBuffer* gbuffer,
                            const std::string& name, glTFAccessorType type,
                            glTFAccessorComponentType ctype, int count,
                            int csize, const void* data, bool save_min_max) {
        gltf->bufferViews.push_back(new glTFBufferView());
        auto bufferView = gltf->bufferViews.back();
        bufferView->buffer = glTFid<glTFBuffer>(index(gltf->buffers, gbuffer));
        bufferView->byteOffset = (int)gbuffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = count * csize;
        gbuffer->data.resize(gbuffer->data.size() + bufferView->byteLength);
        gbuffer->byteLength += bufferView->byteLength;
        auto ptr = gbuffer->data.data() + gbuffer->data.size() -
                   bufferView->byteLength;
        bufferView->target = glTFBufferViewTarget::ArrayBuffer;
        memcpy(ptr, data, bufferView->byteLength);
        gltf->accessors.push_back(new glTFAccessor());
        auto accessor = gltf->accessors.back();
        accessor->bufferView =
            glTFid<glTFBufferView>((int)gltf->bufferViews.size() - 1);
        accessor->byteOffset = 0;
        accessor->componentType = ctype;
        accessor->count = count;
        accessor->type = type;
        if (save_min_max && count &&
            ctype == glTFAccessorComponentType::Float) {
            float dmin[4] = {flt_max, flt_max, flt_max, flt_max};
            float dmax[4] = {flt_min, flt_min, flt_min, flt_min};
            auto d = (float*)data;
            auto nc = 0;
            switch (type) {
                case glTFAccessorType::Scalar: nc = 1; break;
                case glTFAccessorType::Vec2: nc = 2; break;
                case glTFAccessorType::Vec3: nc = 3; break;
                case glTFAccessorType::Vec4: nc = 4; break;
                default: break;
            }
            for (auto i = 0; i < count; i++) {
                for (auto c = 0; c < nc; c++) {
                    dmin[c] = min(dmin[c], d[i * nc + c]);
                    dmax[c] = max(dmax[c], d[i * nc + c]);
                }
            }
            for (auto c = 0; c < nc; c++) {
                accessor->min.push_back(dmin[c]);
                accessor->max.push_back(dmax[c]);
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // instances
    auto shape_mats = std::map<shape*, material*>();
    for (auto shp : scn->shapes) shape_mats[shp] = nullptr;
    for (auto ist : scn->instances) {
        if (shape_mats.at(ist->shp) && shape_mats.at(ist->shp) != ist->mat)
            log_error("shapes can only have one material associated");
        else
            shape_mats[ist->shp] = ist->mat;
    }

    // convert meshes
    for (auto shp : scn->shapes) {
        auto gmesh = new glTFMesh();
        gmesh->name = shp->name;
        auto gbuffer = add_opt_buffer(shp->path);
        auto gprim = new glTFMeshPrimitive();
        gprim->material =
            glTFid<glTFMaterial>(index(scn->materials, shape_mats.at(shp)));
        if (!shp->pos.empty())
            gprim->attributes["POSITION"] =
                add_accessor(gbuffer, shp->name + "_pos",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)shp->pos.size(), sizeof(vec3f), shp->pos.data(), true);
        if (!shp->norm.empty())
            gprim->attributes["NORMAL"] = add_accessor(gbuffer,
                shp->name + "_norm", glTFAccessorType::Vec3,
                glTFAccessorComponentType::Float, (int)shp->norm.size(),
                sizeof(vec3f), shp->norm.data(), false);
        if (!shp->texcoord.empty())
            gprim->attributes["TEXCOORD_0"] = add_accessor(gbuffer,
                shp->name + "_texcoord", glTFAccessorType::Vec2,
                glTFAccessorComponentType::Float, (int)shp->texcoord.size(),
                sizeof(vec2f), shp->texcoord.data(), false);
        if (!shp->color.empty())
            gprim->attributes["COLOR_0"] = add_accessor(gbuffer,
                shp->name + "_color", glTFAccessorType::Vec4,
                glTFAccessorComponentType::Float, (int)shp->color.size(),
                sizeof(vec4f), shp->color.data(), false);
        if (!shp->radius.empty())
            gprim->attributes["RADIUS"] = add_accessor(gbuffer,
                shp->name + "_radius", glTFAccessorType::Scalar,
                glTFAccessorComponentType::Float, (int)shp->radius.size(),
                sizeof(float), shp->radius.data(), false);
        if (!shp->lines.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_lines",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->lines.size() * 2, sizeof(int),
                (int*)shp->lines.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Lines;
        } else if (!shp->triangles.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_triangles",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->triangles.size() * 3, sizeof(int),
                (int*)shp->triangles.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Triangles;
        }
        gmesh->primitives.push_back(gprim);
        gltf->meshes.push_back(gmesh);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // shapes
        for (auto ist : scn->instances) {
            auto gnode = new glTFNode();
            gnode->name = ist->name;
            gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, ist->shp));
            gnode->matrix = frame_to_mat(ist->frame);
            gltf->nodes.push_back(gnode);
        }

        // cameras
        for (auto cam : scn->cameras) {
            auto gnode = new glTFNode();
            gnode->name = cam->name;
            gnode->camera = glTFid<glTFCamera>(index(scn->cameras, cam));
            gnode->matrix = frame_to_mat(cam->frame);
            gltf->nodes.push_back(gnode);
        }

        // scenes
        if (!gltf->nodes.empty()) {
            auto gscene = new glTFScene();
            gscene->name = "scene";
            for (auto i = 0; i < gltf->nodes.size(); i++) {
                gscene->nodes.push_back(glTFid<glTFNode>(i));
            }
            gltf->scenes.push_back(gscene);
            gltf->scene = glTFid<glTFScene>(0);
        }
    } else {
        for (auto nde : scn->nodes) {
            auto gnode = new glTFNode();
            gnode->name = nde->name;
            if (nde->cam) {
                gnode->camera =
                    glTFid<glTFCamera>(index(scn->cameras, nde->cam));
            }
            if (nde->ist) {
                gnode->mesh =
                    glTFid<glTFMesh>(index(scn->shapes, nde->ist->shp));
            }
            gnode->matrix = frame_to_mat(nde->frame);
            gnode->translation = nde->translation;
            gnode->rotation = nde->rotation;
            gnode->scale = nde->scale;
            gltf->nodes.push_back(gnode);
        }

        // children
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto gnde = gltf->nodes.at(index(scn->nodes, nde->parent));
            gnde->children.push_back(glTFid<glTFNode>(idx));
        }

        // root nodes
        auto is_root = std::vector<bool>(gltf->nodes.size(), true);
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            auto gnde = gltf->nodes.at(idx);
            for (auto idx1 = 0; idx1 < gnde->children.size(); idx1++) {
                is_root[(int)gnde->children.at(idx1)] = false;
            }
        }

        // scene with root nodes
        auto gscene = new glTFScene();
        gscene->name = "scene";
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            if (is_root[idx]) gscene->nodes.push_back(glTFid<glTFNode>(idx));
        }
        gltf->scenes.push_back(gscene);
        gltf->scene = glTFid<glTFScene>(0);
    }

    // interpolation map
    static const auto interpolation_map =
        std::map<animation_type, glTFAnimationSamplerInterpolation>{
            {animation_type::step, glTFAnimationSamplerInterpolation::Step},
            {animation_type::linear, glTFAnimationSamplerInterpolation::Linear},
            {animation_type::bezier,
                glTFAnimationSamplerInterpolation::CubicSpline},
        };

    // gruop animations
    struct anim_group {
        std::string path;
        std::string name;
        std::vector<animation*> animations;
    };

    std::map<std::string, anim_group> anim_groups;
    for (auto anm : scn->animations) {
        auto agr = &anim_groups[anm->group];
        if (agr->path == "") agr->path = anm->path;
        agr->animations.push_back(anm);
    }

    // animation
    for (auto& agr_kv : anim_groups) {
        auto agr = &agr_kv.second;
        auto ganm = new glTFAnimation();
        ganm->name = agr->name;
        auto gbuffer = add_opt_buffer(agr->path);
        auto count = 0;
        for (auto anm : scn->animations) {
            auto aid = ganm->name + "_" + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input =
                add_accessor(gbuffer, aid + "_time", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)anm->times.size(),
                    sizeof(float), anm->times.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!anm->translation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->translation.size(), sizeof(vec3f),
                    anm->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!anm->rotation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)anm->rotation.size(), sizeof(vec4f),
                    anm->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!anm->scale.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->scale.size(), sizeof(vec3f), anm->scale.data(),
                    false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!anm->weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(anm->weights.size() * anm->weights[0].size());
                for (auto i = 0; i < anm->weights.size(); i++) {
                    values.insert(values.end(), anm->weights[i].begin(),
                        anm->weights[i].end());
                }
                gsmp->output = add_accessor(gbuffer, aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
                throw std::runtime_error("should not have gotten here");
            }
            gsmp->interpolation = interpolation_map.at(anm->type);
            for (auto target : anm->targets) {
                auto gchan = new glTFAnimationChannel();
                gchan->sampler =
                    glTFid<glTFAnimationSampler>{(int)ganm->samplers.size()};
                gchan->target = new glTFAnimationChannelTarget();
                gchan->target->node =
                    glTFid<glTFNode>{index(scn->nodes, target)};
                gchan->target->path = path;
                ganm->channels.push_back(gchan);
            }
            ganm->samplers.push_back(gsmp);
        }

        gltf->animations.push_back(ganm);
    }

    // done
    return gltf;
}

// Load a scene
scene* load_gltf_scene(const std::string& filename, bool skip_missing) {
    auto gscn = load_gltf(filename, true);
    load_gltf_textures(gscn, path_dirname(filename), skip_missing);
    auto scn = gltf_to_scene(gscn);
    delete gscn;
    auto mat = (material*)nullptr;
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_matte_material("<default>", {0.2f, 0.2f, 0.2f});
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
    return scn;
}

// Save a scene
void save_gltf_scene(const std::string& filename, const scene* scn,
    bool separate_buffers, bool skip_missing) {
    auto buffer_uri = path_basename(filename) + ".bin";
    auto gscn = scene_to_gltf(scn, buffer_uri, separate_buffers);
    save_gltf(filename, gscn);
    save_gltf_textures(gscn, path_dirname(filename), skip_missing);
    delete gscn;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE IO
// -----------------------------------------------------------------------------
namespace ygl {

// Load ply mesh
void load_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    auto ext = path_extension(filename);
    if (ext == ".ply" || ext == ".PLY") {
        load_ply_mesh(
            filename, lines, triangles, pos, norm, texcoord, color, radius);
    } else if (ext == ".obj" || ext == ".OBJ") {
        load_obj_mesh(filename, lines, triangles, pos, norm, texcoord);
    } else {
        throw std::runtime_error("unsupported mesh extensions " + ext);
    }
}

// Save ply mesh
void save_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii) {
    auto ext = path_extension(filename);
    if (ext == ".ply" || ext == ".PLY") {
        save_ply_mesh(filename, lines, triangles, pos, norm, texcoord, color,
            radius, ascii);
    } else if (ext == ".obj" || ext == ".OBJ") {
        save_obj_mesh(filename, lines, triangles, pos, norm, texcoord);
    } else {
        throw std::runtime_error("unsupported mesh extensions " + ext);
    }
}

#if 1

// Load ply mesh
void load_ply_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    // clear data
    pos.clear();
    norm.clear();
    texcoord.clear();
    color.clear();
    radius.clear();
    lines.clear();
    triangles.clear();

    // ply header parsed
    enum property_type { uchar_type, int_type, float_type, int_list_type };
    struct property {
        std::string name = "";
        property_type type = float_type;
        std::vector<float> scalars;
        std::vector<std::array<int, 8>> lists;
    };
    struct element {
        std::string name = "";
        int count = 0;
        std::vector<property> props;
    };

    auto fs = fopen(filename.c_str(), "r");
    if (!fs) throw std::runtime_error("could not open file " + filename);

    // parse header
    char line[4096];
    auto ascii = false;
    auto elems = std::vector<element>();
    while (parse_getline(fs, line, sizeof(line))) {
        auto ss = line;
        auto cmd = parse_string(ss);
        if (cmd.empty()) continue;
        if (cmd == "ply") {
        } else if (cmd == "comment") {
        } else if (cmd == "format") {
            auto fmt = parse_string(ss);
            if (fmt != "ascii" && fmt != "binary_little_endian")
                throw std::runtime_error("format not supported");
            ascii = fmt == "ascii";
        } else if (cmd == "element") {
            auto elem = element();
            elem.name = parse_string(ss);
            elem.count = parse_int(ss);
            elems.push_back(elem);
        } else if (cmd == "property") {
            auto prop = property();
            auto type = parse_string(ss);
            if (type == "list") {
                auto count_type = parse_string(ss);
                auto elem_type = parse_string(ss);
                if (count_type != "uchar" && count_type != "uint8")
                    throw std::runtime_error("unsupported ply list type");
                if (elem_type != "int")
                    throw std::runtime_error("unsupported ply list type");
                prop.type = int_list_type;
            } else if (type == "float") {
                prop.type = float_type;
            } else if (type == "uchar" || type == "uint8") {
                prop.type = uchar_type;
            } else if (type == "int") {
                prop.type = int_type;
            } else {
                throw std::runtime_error("unsupported ply type");
            }
            prop.name = parse_string(ss);
            prop.scalars.resize(elems.back().count);
            if (prop.type == int_list_type)
                prop.lists.resize(elems.back().count);
            elems.back().props.push_back(prop);
        } else if (cmd == "end_header") {
            break;
        } else {
            throw std::runtime_error("command not supported " + cmd);
        }
    }

    // parse content
    for (auto& elem : elems) {
        for (auto vid = 0; vid < elem.count; vid++) {
            auto ss = (char*)nullptr;
            if (ascii) {
                if (!parse_getline(fs, line, sizeof(line)))
                    throw std::runtime_error("error reading ply");
                ss = line;
            }
            for (auto pid = 0; pid < elem.props.size(); pid++) {
                auto& prop = elem.props[pid];
                if (prop.type == float_type) {
                    prop.scalars[vid] =
                        (ascii) ? parse_float(ss) : read_float(fs);
                } else if (prop.type == int_type) {
                    prop.scalars[vid] = (ascii) ? parse_int(ss) : read_int(fs);
                } else if (prop.type == uchar_type) {
                    prop.scalars[vid] = (ascii) ? parse_uchar(ss) / 255.0 :
                                                  read_uchar(fs) / 255.0;
                } else if (prop.type == int_list_type) {
                    prop.scalars[vid] =
                        (ascii) ? parse_uchar(ss) : read_uchar(fs);
                    for (auto i = 0; i < (int)prop.scalars[vid]; i++)
                        prop.lists[vid][i] =
                            (ascii) ? parse_int(ss) : read_int(fs);
                } else {
                    throw std::runtime_error("unsupported ply type");
                }
            }
        }
    }

    // copy vertex data
    for (auto& elem : elems) {
        if (elem.name != "vertex") continue;
        auto count = elem.count;
        for (auto& prop : elem.props) {
            auto vals = prop.scalars.data();
            auto copy_floats = [vals, count](auto& vert, const auto& def,
                                   int stride, int offset) {
                if (vert.size() != count) vert.resize(count, def);
                auto dst = (float*)vert.data();
                for (auto i = 0; i < count; i++)
                    dst[i * stride + offset] = vals[i];
            };
            if (prop.name == "x") copy_floats(pos, zero3f, 3, 0);
            if (prop.name == "y") copy_floats(pos, zero3f, 3, 1);
            if (prop.name == "z") copy_floats(pos, zero3f, 3, 2);
            if (prop.name == "nx") copy_floats(norm, zero3f, 3, 0);
            if (prop.name == "ny") copy_floats(norm, zero3f, 3, 1);
            if (prop.name == "nz") copy_floats(norm, zero3f, 3, 2);
            if (prop.name == "u") copy_floats(texcoord, zero2f, 2, 0);
            if (prop.name == "v") copy_floats(texcoord, zero2f, 2, 1);
            if (prop.name == "red") copy_floats(color, vec4f{0, 0, 0, 1}, 4, 0);
            if (prop.name == "green")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 1);
            if (prop.name == "blue")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 2);
            if (prop.name == "alpha")
                copy_floats(color, vec4f{0, 0, 0, 1}, 4, 3);
            if (prop.name == "radius") copy_floats(radius, 0.0f, 1, 0);
        }
    }

    // copy triangle data
    for (auto& elem : elems) {
        if (elem.name != "face") continue;
        auto count = elem.count;
        for (auto& prop : elem.props) {
            if (prop.name == "vertex_indices") {
                for (auto fid = 0; fid < count; fid++) {
                    auto& list = prop.lists[fid];
                    for (auto i = 2; i < (int)prop.scalars[fid]; i++)
                        triangles.push_back({list[0], list[i - 1], list[i]});
                }
            }
        }
    }

    fclose(fs);
}

#else

// Load ply mesh
void load_ply_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    // clear data
    pos.clear();
    norm.clear();
    texcoord.clear();
    color.clear();
    radius.clear();
    lines.clear();
    triangles.clear();

    // open ply and read its header
    auto error_cb = [](p_ply ply, const char* message) {};
    auto ply = ply_open(filename.c_str(), error_cb, 0, nullptr);
    if (!ply) throw std::runtime_error("could not open ply file " + filename);
    if (!ply_read_header(ply))
        throw std::runtime_error("could not read ply file " + filename);

    // read polygons
    static constexpr int poly_size = 8;
    auto polys = std::vector<std::array<int, poly_size>>();
    auto plines = std::vector<std::array<int, poly_size>>();

    // check contained data
    auto elem = (p_ply_element) nullptr;
    while ((elem = ply_get_next_element(ply, elem)) != nullptr) {
        auto elem_name_buf = (const char*)nullptr;
        auto elem_num = (long)0;
        ply_get_element_info(elem, &elem_name_buf, &elem_num);
        auto elem_name = std::string(elem_name_buf);
        if (elem_name == "vertex") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "x") pos.assign(elem_num, zero3f);
                if (prop_name == "nx") norm.assign(elem_num, zero3f);
                if (prop_name == "u" || elem_name == "s" ||
                    elem_name == "texture_u" || elem_name == "texture_s")
                    texcoord.assign(elem_num, zero2f);
                if (prop_name == "red") color.assign(elem_num, {0, 0, 0, 1});
                if (prop_name == "radius") radius.assign(elem_num, 0);
            }
        } else if (elem_name == "face") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "vertex_indices")
                    polys.assign(elem_num, {{0, 0, 0, 0, 0, 0, 0, 0}});
            }
        } else if (elem_name == "line") {
            auto prop = (p_ply_property) nullptr;
            while ((prop = ply_get_next_property(elem, prop)) != nullptr) {
                auto prop_name_buf = (const char*)nullptr;
                ply_get_property_info(
                    prop, &prop_name_buf, nullptr, nullptr, nullptr);
                auto prop_name = std::string(prop_name_buf);
                if (prop_name == "vertex_indices")
                    plines.assign(elem_num, {{0, 0, 0, 0, 0, 0, 0, 0}});
            }
        }
    }

    // vertex data callback
    auto vert_cb = [](p_ply_argument arg) {
        auto data = (float*)nullptr;
        auto idx = (long)0, flags = (long)0;
        ply_get_argument_user_data(arg, (void**)&data, &flags);
        ply_get_argument_element(arg, nullptr, &idx);
        auto stride = (flags & 0x0F0) >> 4;
        auto offset = flags & 0x00F;
        data[idx * stride + offset] = (float)ply_get_argument_value(arg);
        return 1;
    };

    // face data callback
    auto face_cb = [](p_ply_argument arg) {
        auto data = (int*)nullptr;
        auto idx = (long)0, flags = (long)0, len = (long)0, vidx = (long)0;
        ply_get_argument_user_data(arg, (void**)&data, &flags);
        ply_get_argument_element(arg, nullptr, &idx);
        ply_get_argument_property(arg, nullptr, &len, &vidx);
        if (len >= poly_size) throw std::runtime_error("bad face length");
        data[idx * poly_size + (vidx + 1)] = (int)ply_get_argument_value(arg);
        return 1;
    };

    // set up vertex callbacks
    ply_set_read_cb(ply, "vertex", "x", vert_cb, pos.data(), 0x30);
    ply_set_read_cb(ply, "vertex", "y", vert_cb, pos.data(), 0x31);
    ply_set_read_cb(ply, "vertex", "z", vert_cb, pos.data(), 0x32);
    ply_set_read_cb(ply, "vertex", "nx", vert_cb, norm.data(), 0x30);
    ply_set_read_cb(ply, "vertex", "ny", vert_cb, norm.data(), 0x31);
    ply_set_read_cb(ply, "vertex", "nz", vert_cb, norm.data(), 0x32);
    ply_set_read_cb(ply, "vertex", "u", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "v", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "s", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "t", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "texture_u", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "texture_v", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "texture_s", vert_cb, texcoord.data(), 0x20);
    ply_set_read_cb(ply, "vertex", "texture_t", vert_cb, texcoord.data(), 0x21);
    ply_set_read_cb(ply, "vertex", "red", vert_cb, color.data(), 0x40);
    ply_set_read_cb(ply, "vertex", "green", vert_cb, color.data(), 0x41);
    ply_set_read_cb(ply, "vertex", "blue", vert_cb, color.data(), 0x42);
    ply_set_read_cb(ply, "vertex", "alpha", vert_cb, color.data(), 0x43);
    ply_set_read_cb(ply, "vertex", "radius", vert_cb, radius.data(), 0x10);

    // set up triangle and line callbacks
    ply_set_read_cb(ply, "face", "vertex_indices", face_cb, polys.data(), 0x10);
    ply_set_read_cb(
        ply, "line", "vertex_indices", face_cb, plines.data(), 0x10);

    // read file
    if (!ply_read(ply))
        throw std::runtime_error("error reading ply file " + filename);
    ply_close(ply);

    // convert polygons to triangles
    auto ntriangles = 0;
    for (auto& poly : polys) ntriangles += poly[0] - 2;
    triangles.resize(ntriangles);
    auto ti = 0;
    for (auto& poly : polys) {
        for (auto i = 3; i <= poly[0]; i++)
            triangles[ti++] = {poly[1], poly[i - 1], poly[i]};
    }

    // convert polylines to lines
    auto nlines = 0;
    for (auto& poly : plines) nlines += poly[0] - 1;
    lines.resize(nlines);
    auto li = 0;
    for (auto& pline : plines) {
        for (auto i = 1; i <= pline[0]; i++)
            lines[li++] = {pline[i - 1], pline[i]};
    }
}

#endif

// Save ply mesh
void save_ply_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius,
    bool ascii) {
    auto f = fopen(filename.c_str(), "w");
    if (!f) throw std::runtime_error("cannot save file " + filename);

    // header
    fprintf(f, "ply\n");
    fprintf(f, "format %s 1.0\n", (ascii) ? "ascii" : "binary_little_endian");
    fprintf(f, "element vertex %d\n", (int)pos.size());
    if (!pos.empty())
        fprintf(f, "property float x\nproperty float y\nproperty float z\n");
    if (!norm.empty())
        fprintf(f, "property float nx\nproperty float ny\nproperty float nz\n");
    if (!texcoord.empty()) fprintf(f, "property float u\nproperty float v\n");
    if (!color.empty())
        fprintf(f,
            "property float red\nproperty float green\nproperty float "
            "blue\nproperty float alpha\n");
    if (!radius.empty()) fprintf(f, "property float radius\n");
    if (!triangles.empty()) {
        fprintf(f, "element face %d\n", (int)triangles.size());
        fprintf(f, "property list uchar int vertex_indices\n");
    }
    if (!lines.empty()) {
        fprintf(f, "element line %d\n", (int)triangles.size());
        fprintf(f, "property list uchar int vertex_indices\n");
    }
    fprintf(f, "end_header\n");

    // body
    if (ascii) {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty())
                fprintf(f, "%g %g %g ", pos[i].x, pos[i].y, pos[i].z);
            if (!norm.empty())
                fprintf(f, "%g %g %g ", norm[i].x, norm[i].y, norm[i].z);
            if (!texcoord.empty())
                fprintf(f, "%g %g ", texcoord[i].x, texcoord[i].y);
            if (!color.empty())
                fprintf(f, "%g %g %g %g ", color[i].x, color[i].y, color[i].z,
                    color[i].w);
            if (!radius.empty()) fprintf(f, "%g ", radius[i]);
            fprintf(f, "\n");
        }

        // write face data
        for (auto i = 0; i < triangles.size(); i++) {
            fprintf(f, "3 %d %d %d\n", triangles[i].x, triangles[i].y,
                triangles[i].z);
        }
        for (auto i = 0; i < lines.size(); i++) {
            fprintf(f, "2 %d %d\n", lines[i].x, lines[i].y);
        }
    } else {
        // write vertex data
        for (auto i = 0; i < pos.size(); i++) {
            if (!pos.empty()) fwrite(&pos[i], 4, 3, f);
            if (!norm.empty()) fwrite(&norm[i], 4, 3, f);
            if (!texcoord.empty()) fwrite(&texcoord[i], 4, 2, f);
            if (!color.empty()) fwrite(&color[i], 4, 4, f);
            if (!radius.empty()) fwrite(&radius[i], 4, 1, f);
        }

        // write face data
        for (auto i = 0; i < triangles.size(); i++) {
            auto n = (byte)3;
            fwrite(&n, 1, 1, f);
            fwrite(&triangles[i], 4, 3, f);
        }
        for (auto i = 0; i < lines.size(); i++) {
            auto n = (byte)3;
            fwrite(&n, 1, 1, f);
            fwrite(&lines[i], 4, 2, f);
        }
    }

    // done
    fclose(f);
}

// Load ply mesh
void load_obj_mesh(const std::string& filename, std::vector<vec2i>& lines,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    bool flip_texcoord) {
    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    pos.clear();
    norm.clear();
    texcoord.clear();
    lines.clear();
    triangles.clear();

    // vertex maps
    auto vert_map = std::unordered_map<vec3i, int>();

    // read the file line by line
    char line[4096];
    while (parse_getline(fs, line, sizeof(line))) {
        // prepare to parse
        auto ss = line;
        parse_skipws(ss);
        parse_remove_comment(ss, '#');

        // get command
        auto cmd = parse_string(ss);
        if (cmd == "") continue;

        // possible token values
        if (cmd == "v") {
            pos.push_back(parse_vec3f(ss));
        } else if (cmd == "vn") {
            norm.push_back(parse_vec3f(ss));
        } else if (cmd == "vt") {
            texcoord.push_back(parse_vec2f(ss));
            if (flip_texcoord) texcoord.back().y = 1 - texcoord.back().y;
        } else if (cmd == "f" || cmd == "l") {
            auto num = 0;
            vec3i verts[128];
            int vids[128];
            auto vert_size =
                vec3i{(int)pos.size(), (int)texcoord.size(), (int)norm.size()};
            // elem.material = (int)oobj->materials.size() - 1;
            parse_skipws(ss);
            while (*ss) {
                verts[num++] = parse_objvert(ss, vert_size);
                parse_skipws(ss);
            }
            for (auto i = 0; i < num; i++) {
                auto it = vert_map.find(verts[i]);
                if (it == vert_map.end()) {
                    auto nverts = (int)pos.size();
                    vert_map.insert(it, {verts[i], nverts});
                    vids[i] = nverts;
                    if (verts[i].x >= 0) pos.push_back(pos.at(verts[i].x));
                    if (verts[i].y >= 0)
                        texcoord.push_back(texcoord.at(verts[i].y));
                    if (verts[i].z >= 0) norm.push_back(norm.at(verts[i].z));
                } else {
                    vids[i] = it->second;
                }
            }
            if (cmd == "f") {
                for (auto i = 2; i < num; i++)
                    triangles.push_back({vids[0], vids[i - 1], vids[i]});
            }
            if (cmd == "l") {
                for (auto i = 1; i < num; i++)
                    lines.push_back({vids[i - 1], vids[i]});
            }
        }
    }

    // close file
    fclose(fs);
}

// Load ply mesh
void save_obj_mesh(const std::string& filename, const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    bool flip_texcoord) {
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot save file " + filename);
    for (auto& p : pos) fprintf(f, "v %g %g %g\n", p.x, p.y, p.z);
    for (auto& n : norm) fprintf(f, "vn %g %g %g\n", n.x, n.y, n.z);
    if (flip_texcoord)
        for (auto& t : texcoord) fprintf(f, "vt %g %g\n", t.x, 1 - t.y);
    else
        for (auto& t : texcoord) fprintf(f, "vt %g %g\n", t.x, t.y);
    if (!texcoord.empty() && !norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d/%d/%d %d/%d/%d %d/%d/%d\n", t.x + 1, t.x + 1,
                t.x + 1, t.y + 1, t.y + 1, t.y + 1, t.z + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d/%d/%d %d/%d/%d\n", l.x + 1, l.x + 1, l.x + 1,
                l.y + 1, l.y + 1, l.y + 1);
    } else if (!texcoord.empty() && norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d/%d %d/%d %d/%d\n", t.x + 1, t.x + 1, t.y + 1,
                t.y + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d/%d %d/%d\n", l.x + 1, l.x + 1, l.y + 1, l.y + 1);
    } else if (texcoord.empty() && !norm.empty()) {
        for (auto& t : triangles)
            fprintf(f, "f %d//%d %d//%d %d//%d\n", t.x + 1, t.x + 1, t.y + 1,
                t.y + 1, t.z + 1, t.z + 1);
        for (auto& l : lines)
            fprintf(f, "l %d//%d %d//%d\n", l.x + 1, l.x + 1, l.y + 1, l.y + 1);
    } else {
        for (auto& t : triangles)
            fprintf(f, "f %d %d %d\n", t.x + 1, t.y + 1, t.z + 1);
        for (auto& l : lines) fprintf(f, "l %d %d\n", l.x + 1, l.y + 1);
    }
    fclose(f);
}

}  // namespace ygl
